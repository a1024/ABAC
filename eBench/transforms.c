#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include"intercept_malloc.h"
#include"lodepng.h"//
#include<immintrin.h>
static const char file[]=__FILE__;

void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8)
{
	int nlevels=1<<depth;
	int shift=MAXVAR(0, depth-8);
	memset(hist, 0, nlevels*sizeof(int));
	if(hist8)
		memset(hist8, 0, 256*sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int sym=buf[(iw*ky+kx)<<2|kc]+(nlevels>>1);
			sym=CLAMP(0, sym, nlevels-1);
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
	return entropy;//invCR = entropy/(1<<original_depth)
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
void calc_depthfromdata(int *image, int iw, int ih, char *depths)
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
		int nlevels=vmax-vmin;
		depths[kc]=nlevels?floor_log2(nlevels)+1:0;
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

			r*=(int)(4096LL*(nlevels[0]<<(nlevels[0]<0x1000000))/512);//customparam_ct[0]*1000
			g*=(int)(144LL*nlevels[1]/256);
			b*=(int)(1024*(nlevels[2]<<(nlevels[2]<0x1000000))/512);
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
				Y =(double)image->data[k<<2|0]/(144LL*nlevels[1]/256),
				Cb=(double)image->data[k<<2|1]/(1024LL*(nlevels[2]<<(nlevels[2]<0x1000000))/512),
				Cr=(double)image->data[k<<2|2]/(4096LL*(nlevels[0]<<(nlevels[0]<0x1000000))/512);

			double
				L=Y+Cr,
				M=Y-Cr,//g is the sum
				S=Cb+Y;
			L+=cbrt_bias, L*=L*L, L-=bias;
			M+=cbrt_bias, M*=M*M, M-=bias;
			S+=cbrt_bias, S*=S*S, S-=bias;
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

#define RCT_CUSTOM_NITER 100
#define RCT_CUSTOM_DELTAGROUP 2
int rct_custom_params[RCT_CUSTOM_NPARAMS]={0};
void rct_custom(Image *image, int fwd, const int *params)//4 params	fixed 15.16
{
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
}
static void rct_custom_calcloss(Image const *src, Image *dst, int *hist, const int *params, double *loss)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	memcpy(dst->data, src->data, res*sizeof(int[4]));

	rct_custom(dst, 1, params);

	//pred_grad2(dst, 1);
	pred_clampedgrad(dst, 1, 1);
	
	int depths[]=
	{
		src->depth[1],
		(src->depth[2]+1),
		(src->depth[0]+1),
	};
	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(dst->data, dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, depths[kc], hist, 0);
		double entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
		loss[kc]=entropy/src->src_depth[(kc+1)%3];
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
void rct_custom_optimize(Image const *image, int *params)
{
	static int call_idx=0;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();

	Image *im2=0;
	image_copy(&im2, image);
	int inflation[]={1, 0, 1};//original image is expected, RCT permutes depths {r+1, g, b+1} -> {Y:g:+0, Cb:b:+1, Cr:r:+1}
	int maxdepth=calc_maxdepth(image, inflation);
	int maxlevels=1<<maxdepth;
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	if(!im2||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}

	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	int params2[RCT_CUSTOM_NPARAMS];
	memcpy(params2, params, sizeof(int[RCT_CUSTOM_NPARAMS]));

#define CALC_LOSS(L) rct_custom_calcloss(image, im2, hist, params2, L)
	srand((unsigned)__rdtsc());//
	CALC_LOSS(loss_bestsofar);
	memcpy(loss_prev, loss_bestsofar, sizeof(loss_prev));

	int shakethreshold=RCT_CUSTOM_NPARAMS;
	for(int it=0, watchdog=0;it<RCT_CUSTOM_NITER;++it)
	{
		int idx[RCT_CUSTOM_DELTAGROUP]={0}, stuck=0;
		int params_original_selected[RCT_CUSTOM_DELTAGROUP]={0};
		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(params2, params, sizeof(params2));
			for(int k=0;k<RCT_CUSTOM_NPARAMS;++k)
				params2[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<RCT_CUSTOM_DELTAGROUP;++k)//increment params
			{
				int inc=0;
				idx[k]=rand()%RCT_CUSTOM_NPARAMS;
				while(!(inc=rand()-(RAND_MAX>>1)));//reject zero delta
		
				params_original_selected[k]=params2[idx[k]];
				params2[idx[k]]+=inc*(16<<1)/RAND_MAX;
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
					params2[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss_curr[3]<loss_bestsofar[3])//publish if record best
			{
				memcpy(params, params2, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;//bis
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}

		set_window_title(
			"%d %4d/%4d,%d/%d: %lf RGB %lf %lf %lf%s",
			call_idx,
			it+1,
			RCT_CUSTOM_NITER,
			watchdog,
			shakethreshold,
			1/loss_bestsofar[3],
			1/loss_bestsofar[0],
			1/loss_bestsofar[1],
			1/loss_bestsofar[2],
			it+1<RCT_CUSTOM_NITER?"...":" Done."
		);

		//preview
#if 1
		{
			ch_cr[0]=(float)(1/loss_bestsofar[0]);
			ch_cr[1]=(float)(1/loss_bestsofar[1]);
			ch_cr[2]=(float)(1/loss_bestsofar[2]);
			//unsigned char *ptr;
			//addhalf(temp, iw, ih, 3, 4);
			//SWAPVAR(image, temp, ptr);
			io_render();
			//SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	free(hist);
	free(im2);
}



//void addhalf(unsigned char *buf, int iw, int ih, int nch, int bytestride)
//{
//	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
//	{
//		for(int kc=0;kc<nch;++kc)
//			buf[kp+kc]+=128;
//	}
//}
void c_mul_add(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]+=r, dst[1]+=i;
}
void c_mul_sub(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]-=r, dst[1]-=i;
}
double c_abs2(const double *z)
{
	return z[0]*z[0]+z[1]*z[1];
}
void c_div(double *dst, double const *a, double const *b)
{
	double invabsb2=1/c_abs2(b);
	double
		r=(a[0]*b[0]+a[1]*b[1])*invabsb2,
		i=(a[1]*b[0]-a[0]*b[1])*invabsb2;
	dst[0]=r, dst[1]=i;
}
void c_exp(double *dst, double const *x)
{
	double m=exp(x[0]);
	dst[0]=m*cos(x[1]);
	dst[1]=m*sin(x[1]);
}
void c_ln(double *dst, double const *x)
{
	double
		r=log(sqrt(c_abs2(x))),
		i=atan2(x[1], x[0]);
	dst[0]=r, dst[1]=i;
}
void c_cbrt(double *dst, const double *x)//sqrt(x)=exp(0.5lnx)
{
	double temp[2];
	if(x[0]||x[1])
	{
		c_ln(temp, x);
		temp[0]*=1./3, temp[1]*=1./3;
		c_exp(dst, temp);
	}
	else
		dst[0]=x[0], dst[1]=x[1];
}
void impl_solve_cubic(const double *coeffs, double *roots)//finds r[0], r[1], & r[2], the solutions of x^3 + c[2]x^2 + c[1]x + c[0] = 0
{
	//https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
	double p=coeffs[2], q=coeffs[1], r=coeffs[0],
		p2, p3, q2, pq, A[2], B[2];
	double sqrt27=3*sqrt(3), inv3cbrt2=1./(3*cbrt(2)), ninth=1./9;
	double cm[]={-0.5, -sqrt(3)*0.5}, cp[]={-0.5, sqrt(3)*0.5};

	if(!p&&!q)
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
void impl_rref(double *m, short dx, short dy)
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
int impl_nullspace(double *M, int dx, int dy, double *solution, int sstart, char *dep_flags, short *row_idx)
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
	impl_rref(M, dx, dy);//first, get rref(M)
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
			dep_flags[kx]=1, row_idx[ky]=kx;
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
int impl_egvec(double const *M, int n, const double *lambdas, double *S)
{
	int kv, kx, nvec, size=n*n, again;
	double *temp=(double*)malloc(size*sizeof(double));
	char *dep_flags=(char*)malloc(n);
	short *row_idx=(short*)malloc(n*sizeof(short));
	if(!temp||!dep_flags||!row_idx)
	{
		LOG_ERROR("Allocation error");
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
void impl_egval3(double const *M, double *lambdas)//finds the complex eigenvalues of a 3x3 matrix
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
static void rct_adaptive_block(int *src, int iw, int ih, int x, int y, int blocksize, double *axes)
{
	int res=iw*ih;
	double mean[3]={0}, cov[9]={0};
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
	double lambdas[6];
	impl_egval3(cov, lambdas);//get eigenvalues
	lambdas[0]=sqrt(lambdas[0]*lambdas[0]+lambdas[1]*lambdas[1]);
	lambdas[1]=sqrt(lambdas[2]*lambdas[2]+lambdas[3]*lambdas[3]);
	lambdas[2]=sqrt(lambdas[4]*lambdas[4]+lambdas[5]*lambdas[5]);
	double ev[9];
	int nv=impl_egvec(cov, 3, lambdas, ev);//get eigenvectors
	if(nv<2)
		return;

	double d;
#define NORMALIZE(IDX) d=ev[IDX+0]*ev[IDX+0]+ev[IDX+1]*ev[IDX+1]+ev[IDX+2]*ev[IDX+2], d=d?1/sqrt(d):1, ev[IDX+0]*=d, ev[IDX+1]*=d, ev[IDX+2]*=d
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
#define NORMALIZE(X) d=X[0]*X[0]+X[1]*X[1]+X[2]*X[2], d=d?1/sqrt(d):1, X[0]*=d, X[1]*=d, X[2]*=d
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
		int kby1=(ky+(blocksize>>1))/blocksize;
		int kby0=kby1-(kby1>0);
		kby1-=kby1>=by;
		int alphay=(ky+(blocksize>>1))%blocksize;
		for(int kx=0;kx<src->iw;++kx)
		{
			int kbx1=(kx+(blocksize>>1))/blocksize;
			int kbx0=kbx1-(kbx1>0);
			kbx1-=kbx1>=bx;
			int alphax=(kx+(blocksize>>1))%blocksize;
			double
				*axesTL=axes+12*(bx*kby0+kbx0),
				*axesTR=axes+12*(bx*kby0+kbx1),
				*axesBL=axes+12*(bx*kby1+kbx0),
				*axesBR=axes+12*(bx*kby1+kbx1);
			//if(kby0!=kby1&&kbx0!=kbx1)
			//	printf("");
			double localaxes[12];
			for(int k=0;k<12;++k)//alpha blend
			{
				double vTL=axesTL[k], vTR=axesTR[k], vBL=axesBL[k], vBR=axesBR[k];
				double vT=vTL+(vTR-vTL)*alphax/blocksize;
				double vB=vBL+(vBR-vBL)*alphax/blocksize;
				double v=vT+(vB-vT)*alphay/blocksize;
				localaxes[k]=v;
			}
			double d;
#define NORMALIZE(IDX) d=localaxes[IDX+0]*localaxes[IDX+0]+localaxes[IDX+1]*localaxes[IDX+1]+localaxes[IDX+2]*localaxes[IDX+2], d=d?1/sqrt(d):1, localaxes[IDX+0]*=d, localaxes[IDX+1]*=d, localaxes[IDX+2]*=d
			NORMALIZE(0);//normalize axes
			NORMALIZE(3);
			NORMALIZE(6);
#undef  NORMALIZE
			int idx=(src->iw*ky+kx)<<2;
			double r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];
			r-=localaxes[ 9];//subtract mean
			g-=localaxes[10];
			b-=localaxes[11];
			double
				r2=localaxes[6]*r+localaxes[7]*g+localaxes[8]*b,//Cr	project on axes
				g2=localaxes[0]*r+localaxes[1]*g+localaxes[2]*b,//Y
				b2=localaxes[3]*r+localaxes[4]*g+localaxes[5]*b;//Cb
			src->data[idx|0]=(int)g2;
			src->data[idx|1]=(int)b2;
			src->data[idx|2]=(int)r2;
		}
	}
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
#define NORMALIZE(IDX) d=ev[IDX+0]*ev[IDX+0]+ev[IDX+1]*ev[IDX+1]+ev[IDX+2]*ev[IDX+2], d=d?1/sqrt(d):1, ev[IDX+0]*=d, ev[IDX+1]*=d, ev[IDX+2]*=d
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
#define NORMALIZE(X) d=X[0]*X[0]+X[1]*X[1]+X[2]*X[2], d=d?1/sqrt(d):1, X[0]*=d, X[1]*=d, X[2]*=d
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
			for(int k=0;k<_countof(lumas);++k)
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

//clamped gradient / LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS
void pred_clampedgrad(Image *src, int fwd, int enable_ma)
{
	Image *dst=0;
	image_copy(&dst, src);
	int *pixels=fwd?src->data:dst->data, *errors=fwd?dst->data:src->data;
	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int
					NW=kx&&ky?pixels[(idx-src->iw-1)<<2|kc]:0,
					N =ky    ?pixels[(idx-src->iw  )<<2|kc]:0,
					W =kx    ?pixels[(idx          -1)<<2|kc]:0;
				int pred=N+W-NW;
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
	}
	free(dst);
}

//CUSTOM reach-2 predictor
int custom_params[CUSTOM_NPARAMS]={0};
void pred_custom(Image *src, int fwd, int enable_ma, const int *params)
{
	Image *dst=0;
	image_copy(&dst, src);
	int *pixels=fwd?src->data:dst->data, *errors=fwd?dst->data:src->data;
	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[CUSTOM_NNB*2];
				for(int ky2=-CUSTOM_REACH, idx2=0;ky2<=0;++ky2)
				{
					for(int kx2=-CUSTOM_REACH;kx2<=CUSTOM_REACH&&idx2<CUSTOM_NNB*2;++kx2, idx2+=2)
					{
						if((unsigned)(ky+ky2)<(unsigned)src->ih&&(unsigned)(kx+kx2)<(unsigned)src->iw)
						{
							nb[idx2  ]=pixels[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
							nb[idx2+1]=errors[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
						}
						else
						{
							nb[idx2  ]=0;
							nb[idx2+1]=0;
						}
					}
				}

				long long pred=0;
				for(int k=0;k<CUSTOM_NNB*2;++k)
					pred+=(long long)nb[k]*params[k];
				pred+=1LL<<15;
				pred>>=16;
				pred=CLAMP(-(nlevels>>1), pred, (nlevels>>1));

				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=(int)pred;
			}
		}
		params+=24;
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
	}
	free(dst);
}
#define CUSTOM_NITER 100
#define CUSTOM_DELTAGROUP 4
static void pred_custom_calcloss(Image const *src, Image *dst, int *hist, const int *params, double *loss)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	memcpy(dst->data, src->data, res*sizeof(int[4]));

	pred_custom(dst, 1, 1, params);

	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		calc_histogram(dst->data, dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, src->depth[kc], hist, 0);
		double entropy=calc_entropy(hist, nlevels, (int)res);
		loss[kc]=entropy/src->src_depth[(kc+1)%3];
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
void pred_custom_optimize(Image const *image, int *params)
{
	static int call_idx=0;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();

	Image *im2=0;
	image_copy(&im2, image);
	int maxdepth=calc_maxdepth(image, 0);
	int maxlevels=1<<maxdepth;
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	if(!im2||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}

	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	int params2[CUSTOM_NPARAMS];
	memcpy(params2, params, sizeof(params2));

#define CALC_LOSS(L) pred_custom_calcloss(image, im2, hist, params2, L)
	srand((unsigned)__rdtsc());//
	CALC_LOSS(loss_bestsofar);
	memcpy(loss_prev, loss_bestsofar, sizeof(loss_prev));

	int shakethreshold=CUSTOM_NPARAMS;
	for(int it=0, watchdog=0;it<CUSTOM_NITER;++it)
	{
		int idx[CUSTOM_DELTAGROUP]={0}, stuck=0;
		int params_original_selected[CUSTOM_DELTAGROUP]={0};
		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(params2, params, sizeof(params2));
			for(int k=0;k<CUSTOM_NPARAMS;++k)
				params2[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<CUSTOM_DELTAGROUP;++k)//increment params
			{
				int inc=0;
				idx[k]=rand()%CUSTOM_NPARAMS;
				while(!(inc=rand()-(RAND_MAX>>1)));//reject zero delta
		
				params_original_selected[k]=params2[idx[k]];
				params2[idx[k]]+=inc*(16<<1)/RAND_MAX;
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
					params2[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss_curr[3]<loss_bestsofar[3])//publish if record best
			{
				memcpy(params, params2, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;//bis
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}

		set_window_title(
			"%d %4d/%4d,%d/%d: %lf RGB %lf %lf %lf%s",
			call_idx,
			it+1,
			CUSTOM_NITER,
			watchdog,
			shakethreshold,
			1/loss_bestsofar[3],
			1/loss_bestsofar[0],
			1/loss_bestsofar[1],
			1/loss_bestsofar[2],
			it+1<CUSTOM_NITER?"...":" Done."
		);

		//preview
#if 1
		{
			ch_cr[0]=(float)(1/loss_bestsofar[0]);
			ch_cr[1]=(float)(1/loss_bestsofar[1]);
			ch_cr[2]=(float)(1/loss_bestsofar[2]);
			//unsigned char *ptr;
			//addhalf(temp, iw, ih, 3, 4);
			//SWAPVAR(image, temp, ptr);
			io_render();
			//SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	free(hist);
	free(im2);
}




double pw2_errors[PW2_NPRED]={0};//
short pw2_params[PW2_NPARAM*3]=
{
	//0
	
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,-0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C, 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,

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
static int clip(int x, int nlevels)
{
	x=CLAMP(-nlevels, x, nlevels-1);
	return x;
}
void pred_w2_prealloc(const int *src, int iw, int ih, int depth, int kc, short *params, int fwd, int enable_ma, int *dst, int *temp)//temp is (PW2_NPRED+1)*2w
{
	int errorbuflen=iw<<1, rowlen=iw<<2;
	int *error=temp, *pred_errors[PW2_NPRED];
	for(int k=0;k<PW2_NPRED;++k)
		pred_errors[k]=temp+errorbuflen*(k+1);
	int idx=kc;

	int nlevels=1<<depth;
	const int *src2=fwd?src:dst;
	memset(pw2_errors, 0, sizeof(pw2_errors));
	for(int ky=0;ky<ih;++ky)
	{
		//int pred_left=0, pred_left2=0;

		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			//           T3
			//   T2L2    T2    T2R2
			//        TL T  TR TR2
			//L3 L2   L  X
			int
				cT6  =         ky-6>=0?src2[idx-rowlen*6   ]<<8:0,
				cT5  =         ky-5>=0?src2[idx-rowlen*5   ]<<8:0,

				cT4L3=kx-3>=0&&ky-4>=0?src2[idx-rowlen*4-12]<<8:0,
				cT4  =         ky-4>=0?src2[idx-rowlen*4   ]<<8:0,
				cT4R3=kx+3<iw&&ky-4>=0?src2[idx-rowlen*4+12]<<8:0,
				
				cT3L5=kx-5>=0&&ky-3>=0?src2[idx-rowlen*3-20]<<8:0,
				cT3L4=kx-4>=0&&ky-3>=0?src2[idx-rowlen*3-16]<<8:0,
				cT3L2=kx-2>=0&&ky-3>=0?src2[idx-rowlen*3- 8]<<8:0,
				cT3L =kx-1>=0&&ky-3>=0?src2[idx-rowlen*3- 4]<<8:0,
				cT3  =         ky-3>=0?src2[idx-rowlen*3   ]<<8:0,
				cT3R =kx+1<iw&&ky-3>=0?src2[idx-rowlen*3+ 4]<<8:0,
				cT3R2=kx+2<iw&&ky-3>=0?src2[idx-rowlen*3+ 8]<<8:0,
				cT3R3=kx+3<iw&&ky-3>=0?src2[idx-rowlen*3+12]<<8:0,
				cT3R4=kx+4<iw&&ky-3>=0?src2[idx-rowlen*3+16]<<8:0,
				
				cT2L3=kx-3>=0&&ky-2>=0?src2[idx-rowlen*2-12]<<8:0,
				cT2L2=kx-2>=0&&ky-2>=0?src2[idx-rowlen*2- 8]<<8:0,
				cT2L =kx-1>=0&&ky-2>=0?src2[idx-rowlen*2- 4]<<8:0,
				cT2  =         ky-2>=0?src2[idx-rowlen*2   ]<<8:0,
				cT2R =kx+1<iw&&ky-2>=0?src2[idx-rowlen*2+ 4]<<8:0,
				cT2R2=kx+2<iw&&ky-2>=0?src2[idx-rowlen*2+ 8]<<8:0,
				cT2R3=kx+3<iw&&ky-2>=0?src2[idx-rowlen*2+12]<<8:0,
				cT2R4=kx+4<iw&&ky-2>=0?src2[idx-rowlen*2+16]<<8:0,
				
				cTL3 =kx-3>=0&&ky-1>=0?src2[idx-rowlen  -12]<<8:0,
				cTL2 =kx-2>=0&&ky-1>=0?src2[idx-rowlen  - 8]<<8:0,
				cTL  =kx-1>=0&&ky-1>=0?src2[idx-rowlen  - 4]<<8:0,
				cT   =         ky-1>=0?src2[idx-rowlen     ]<<8:0,
				cTR  =kx+1<iw&&ky-1>=0?src2[idx-rowlen  + 4]<<8:0,
				cTR2 =kx+2<iw&&ky-1>=0?src2[idx-rowlen  + 8]<<8:0,
				cTR3 =kx+3<iw&&ky-1>=0?src2[idx-rowlen  +12]<<8:0,
				cTR4 =kx+4<iw&&ky-1>=0?src2[idx-rowlen  +16]<<8:0,
				cTR5 =kx+5<iw&&ky-1>=0?src2[idx-rowlen  +20]<<8:0,
				cTR6 =kx+6<iw&&ky-1>=0?src2[idx-rowlen  +24]<<8:0,
				cTR7 =kx+7<iw&&ky-1>=0?src2[idx-rowlen  +28]<<8:0,

				cL6  =kx-6>=0         ?src2[idx         -24]<<8:0,
				cL5  =kx-5>=0         ?src2[idx         -20]<<8:0,
				cL4  =kx-4>=0         ?src2[idx         -16]<<8:0,
				cL3  =kx-2>=0         ?src2[idx         -12]<<8:0,
				cL2  =kx-2>=0         ?src2[idx         - 8]<<8:0,
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
				int w=
					 (ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0)
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0);
				weights[k]=(params[k]<<8)/(w+1);
			}

			//TL T TR
			//L  X
			int
				eT=ky-1>=0?error[prevrow+kx]:0,
				eL=kx-1>=0?error[currrow+kx-1]:0,
				eTL=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				eTR=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
				eT_L=eT+eL;

			//pred_left=pred_left+((cL-pred_left)*params[PW2_NPRED+11]>>8);
			//pred_left2=pred_left2+((cL-pred_left2)*abs(eL)>>16);

			//int pred_w=(cT*(0x10000-eT)+cL*(0x10000-eL)+cTL*(0x10000-eTL)+cTR*(0x10000-eTR))>>16;

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

			int vmin=cL, vmax=cL, curr;
			if(vmin>cT)vmin=cT;
			if(vmax<cT)vmax=cT;
			//if(vmin>cTR)vmin=cTR;
			//if(vmax<cTR)vmax=cTR;
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

			error[currrow+kx]=curr-(int)pred;
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
	for(int k=0;k<PW2_NPRED;++k)
		pw2_errors[k]/=iw*ih*256;
}
static double pred_w2_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, short *params, int *temp, int *dst, int *hist)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	pred_w2_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	double entropy=calc_entropy(hist, 1<<depth, (int)res);

	double invCR=entropy/src_depth, csize=res*invCR;
	return csize;
}
void pred_w2_opt_v2(Image *src, short *params, int loud)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int maxdepth=calc_maxdepth(src, 0);
	int nlevels=1<<maxdepth;
	int *hist=(int*)malloc(nlevels*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	char title0[256];
	get_window_title(title0, 256);
	int steps[]={128, 64, 32, 16, 8, 4, 2, 1};
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*PW2_NPARAM;
		double csize0=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
		for(int ks=0;ks<_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<PW2_NPARAM;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=step;
				csize=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;

				set_window_title("Ch%d csize %lf [%d/%d %d/%d]...", kc, csize0, kc*_countof(steps)+ks+1, _countof(steps)*3, idx+1, PW2_NPARAM);//
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
			//set_window_title("Ch%d csize %lf [%d/%d]...", kc, csize0, kc*_countof(steps)+ks+1, _countof(steps)*3);//
		}
	}
	free(hist);
	free(temp);
	free(buf3);
	set_window_title("%s", title0);//
}
void pred_w2_apply(Image *src, int fwd, int enable_ma, short *params)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, params             , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[1], 1, params+PW2_NPARAM  , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, params+PW2_NPARAM*2, fwd, enable_ma, buf2, temp);

	for(ptrdiff_t k=0;k<res;++k)
	{
		src->data[k<<2|0]=buf2[k<<2|0];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}


short jxlparams_i16[33]=//signed fixed 7.8 bit
{
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,

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
void pred_jxl_prealloc(const int *src, int iw, int ih, int depth, int kc, short *params, int fwd, int enable_ma, int *dst, int *temp_w10)
{
	int errorbuflen=iw<<1, rowlen=iw<<2;
	int *error=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int nlevels=1<<depth;
	int idx=kc;
	const int *src2=fwd?src:dst;
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred, curr;
			
			char
				ctt      =         ky-2>=0?src2[idx-rowlen*2]:0,
				ctopleft =kx-1>=0&&ky-1>=0?src2[idx-rowlen-4]:0,
				ctop     =kx  <iw&&ky-1>=0?src2[idx-rowlen  ]:0,
				ctopright=kx+1<iw&&ky-1>=0?src2[idx-rowlen+4]:0,
				cleft    =kx-1>=0         ?src2[idx       -4]:0;

			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	kx=iw>>1;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c
			
			int weights[4];//fixed 23.8 bit
			for(int k=0;k<4;++k)
			{
				int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
				weights[k]=(params[k]<<8)/(w+1);
			}

			int
				etop=ky-1>=0?error[prevrow+kx]:0,
				eleft=kx-1>=0?error[currrow+kx-1]:0,
				etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
				etopplusleft=etop+eleft;
			int predictions[]=//fixed 23.8 bit
			{
				(cleft+ctopright-ctop)<<8,
				(ctop<<8)-((etopplusleft+etopright)*params[4]>>8),
				(cleft<<8)-((etopplusleft+etopleft)*params[5]>>8),
				(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
			};

			int sum=weights[0]+weights[1]+weights[2]+weights[3];
			if(sum)
				pred=(predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum;
			else
				pred=predictions[0];

			int vmin=cleft, vmax=cleft;
			if(vmin>ctopright)vmin=ctopright;
			if(vmax<ctopright)vmax=ctopright;

			if(vmin>ctop)vmin=ctop;
			if(vmax<ctop)vmax=ctop;

			vmin<<=8;
			vmax<<=8;

			pred=CLAMP(vmin, pred, vmax);

			if(fwd)
			{
				curr=src[idx]<<8;
				dst[idx]=src[idx]-((pred+127)>>8);
				if(enable_ma)
					dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
			}
			else
			{
				dst[idx]=src[idx]+((pred+127)>>8);
				if(enable_ma)
					dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
				curr=dst[idx]<<8;
			}

			error[currrow+kx]=curr-pred;
			for(int k=0;k<4;++k)
			{
				int e=abs(curr-predictions[k]);
				pred_errors[k][currrow+kx]=e;
				if(kx+1<iw)
					pred_errors[k][prevrow+kx+1]+=e;
			}
		}
	}
}
double pred_jxl_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, short *params, int *temp, int *dst, int *hist)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	pred_jxl_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	//addhalf((unsigned char*)dst+kc, iw, ih, 1, 4);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	double entropy=calc_entropy(hist, 1<<depth, (int)res);
	double invCR=entropy/src_depth, csize=res*invCR;
	return csize;
}
void pred_jxl_opt_v2(Image *src, short *params, int loud)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	int *temp=(int*)malloc((size_t)src->iw*10*sizeof(int));
	int maxdepth=calc_maxdepth(src, 0);
	int nlevels=1<<maxdepth;
	int *hist=(int*)malloc(nlevels*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int steps[]={256, 128, 64, 32, 16, 8, 4, 2, 1};
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*11;
		double csize0=pred_jxl_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
		for(int ks=0;ks<_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<11;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=step;
				csize=pred_jxl_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=step;
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

				param[bestidx]+=beststep;
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
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, params   , fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[1], 1, params+11, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, params+22, fwd, enable_ma, buf2, temp);

	for(int k=0;k<res;++k)
	{
		src->data[k<<2  ]=buf2[k<<2  ];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
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
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, jxlparams_i16        , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc (src->data, src->iw, src->ih, src->depth[1], 1, pw2_params+PW2_NPARAM, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, jxlparams_i16+22     , fwd, enable_ma, buf2, temp);

	for(int k=0;k<res;++k)
	{
		src->data[k<<2  ]=buf2[k<<2  ];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}


//CUSTOM3
#define C3_OPT_NCOMP (C3_NPARAMS>>5)
Custom3Params c3_params={0};
int fast_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-31;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((__m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((__m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
	__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
	s2=_mm_add_epi32(s2, hi);
	s2=_mm_hadd_epi32(s2, s2);
	int s3=_mm_extract_epi32(s2, 0);
	for(;k<count;++k)
		s3+=a[k]*b[k];
	return s3;
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
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
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
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
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
			int count[3], idx2;
			for(int kc=0;kc<3;++kc)
				count[kc]=custom3_loadnb(pixels, errors, iw, ih, kc, kx, ky, nb[kc]);
			
			idx=(iw*ky+kx)<<2;
			idx2=0;
			
			for(int kdst=0;kdst<3;++kdst)
			{
				pred=0;
				for(int kc=0;kc<3;++kc)
					pred+=fast_dot(coeffs[idx2+kc], nb[kc], count[kc]);
				pred+=1<<13;
				pred>>=14;
				pred=CLAMP(-(nlevels[kdst]>>1), pred, (nlevels[kdst]>>1)-1);

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
				//if(fwd)
				//	dst[idx]=src[idx]-pred;
				//else
				//	dst[idx]=src[idx]+pred;

				nb[kdst][C3_NNB  ]=pixels[idx];
				nb[kdst][C3_NNB+1]=errors[idx];
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
		LOG_ERROR("Allocation error");
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

	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(temp, iw, ih, kc, 0, iw, 0, ih, depths[kc], hist, 0);
		double entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
		info->invCR[kc]=entropy/src_depths[kc];
		info->invCR[3]+=info->invCR[kc];
	}
	info->invCR[3]/=3;
}
void custom3_opt(Image const *src, Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();

	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	int *temp=(int*)malloc(res*sizeof(int[4]));
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0xFF, res*sizeof(int[4]));//set alpha for preview
	Custom3OptInfo info;
	memcpy(info.params, srcparams, sizeof(info.params));
#define CALC_LOSS() custom3_calcloss(src->data, src->iw, src->ih, src->depth, src->src_depth, &info, temp, hist)
	if(loud)
		srand((unsigned)__rdtsc());//

	CALC_LOSS();
	double invCR[4], loss0;
	memcpy(invCR, info.invCR, sizeof(invCR));
	loss0=info.invCR[3];
	if(!niter)
		niter=C3_NPARAMS*10;
	int shakethreshold=C3_NPARAMS;
	for(int it=0, watchdog=0;it<niter;++it)
	{
		if(loud)
			set_window_title("%d-%d %4d/%4d,%d/%d: %lf RGB %lf %lf %lf", loud, call_idx, it+1, niter, watchdog, shakethreshold, 1/info.invCR[3], 1/info.invCR[0], 1/info.invCR[1], 1/info.invCR[2]);
		int idx[C3_OPT_NCOMP]={0}, inc[C3_OPT_NCOMP]={0}, stuck=0;
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
				info.params[idx[k]]+=inc[k];
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
					info.params[idx[k]]-=inc[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss0>info.invCR[3])
			{
				memcpy(srcparams, info.params, sizeof(info.params));
				loss0=info.invCR[3];
				--it;//bis
			}
			memcpy(invCR, info.invCR, sizeof(invCR));
			watchdog=0;
		}

		//preview
#if 1
		if(loud)
		{
			ch_cr[0]=(float)(1/info.invCR[0]);
			ch_cr[1]=(float)(1/info.invCR[1]);
			ch_cr[2]=(float)(1/info.invCR[2]);
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
	int *b2=(int*)malloc(res*sizeof(int[4]));
	int *arrN=(int*)malloc(2048*sizeof(int));//count
	int *arrS=(int*)malloc(2048*sizeof(int));//sum
	if(!b2||!arrN||!arrS)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	//const int thresholds[]={5, 15, 25, 42, 60, 85, 140};
	const int *pixels=fwd?src->data:b2, *errors=fwd?b2:src->data;
	int idx, pred;
	memcpy(b2, src->data, (size_t)res*sizeof(int[4]));//copy alpha
	for(int kc=0;kc<3;++kc)//process each channel separately
	{
		int nlevels=1<<src->depth[kc];
		memset(arrN, 0, 2048*sizeof(int));
		memset(arrS, 0, 2048*sizeof(int));
		for(int ky=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
#define LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				int
					NNWW=LOAD(pixels, -2, -2),
					NNW =LOAD(pixels, -1, -2),
					NN  =LOAD(pixels,  0, -2),
					NNE =LOAD(pixels,  1, -2),
					NNEE=LOAD(pixels,  2, -2),
					NWW =LOAD(pixels, -2, -1),
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

				int energy=dx+dy+abs(eN)+abs(eW);
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

				int pattern=(2*W-WW<pred)<<7|(2*N-NN<pred)<<6|(WW<pred)<<5|(NN<pred)<<4|(NE<pred)<<3|(NW<pred)<<2|(W<pred)<<1|(N<pred);

				int context=pattern<<3|energy;

				if(arrN[context]>0)
					pred+=(arrS[context]+arrN[context]/2)/arrN[context];//sum of encountered errors / number of occurrences
				pred=CLAMP(-128, pred, 127);

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
typedef struct SseCellStruct
{
	int count, sum;//, sum7, sum8;
} SseCell;
SseCell g2_SSE[12][256];
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
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *b2=(int*)malloc((size_t)res*sizeof(int[4]));
	int *perrors=(int*)malloc((size_t)src->iw*(G2_NPRED+1)*2*sizeof(int));
	//int *SSE_count=(int*)malloc(256*sizeof(int));
	//int *SSE_sum=(int*)malloc(256*sizeof(int));
	if(!b2||!perrors)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, src->data, res*sizeof(int[4]));//copy alpha
	const int *pixels=fwd?src->data:b2, *errors=fwd?b2:src->data;
	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		int maxerror=0;
		short *params=g2_weights+(_countof(g2_weights)/3)*kc;
		int *hireserror=perrors+src->iw*2*G2_NPRED;
		memset(perrors, 0, 2LL*src->iw*(G2_NPRED+1)*sizeof(int));

		memset(g2_SSE, 0, sizeof(g2_SSE));
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
					NNNNWWWW=LOAD(pixels, -4, -4),
					NNNN    =LOAD(pixels,  0, -4),
					NNNNEEEE=LOAD(pixels,  4, -4),
					NNNWWW  =LOAD(pixels, -3, -3),
					NNN     =LOAD(pixels,  0, -3),
					NNNEEE  =LOAD(pixels,  3, -3),
					NNWW    =LOAD(pixels, -2, -2),
					NNW     =LOAD(pixels, -1, -2),
					NN      =LOAD(pixels,  0, -2),
					NNE     =LOAD(pixels,  1, -2),
					NNEE    =LOAD(pixels,  2, -2),
					NW      =LOAD(pixels, -1, -1),
					N       =LOAD(pixels,  0, -1),
					NE      =LOAD(pixels,  1, -1),
					NEEEE   =LOAD(pixels,  4, -1),
					NEEEEE  =LOAD(pixels,  5, -1),
					NEEEEEE =LOAD(pixels,  6, -1),
					NEEEEEEE=LOAD(pixels,  7, -1),
					NEE     =LOAD(pixels,  2, -1),
					WWWWWW  =LOAD(pixels, -6,  0),
					WWWW    =LOAD(pixels, -4,  0),
					WWW     =LOAD(pixels, -3,  0),
					WW      =LOAD(pixels, -2,  0),
					W       =LOAD(pixels, -1,  0);
#undef  LOAD
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0
				int
					dNNNWWW=LOAD(errors, -3, -3),
					dNNNWW =LOAD(errors, -2, -3),
					dNNNW  =LOAD(errors, -1, -3),
					dNNN   =LOAD(errors,  0, -3),
					dNNNE  =LOAD(errors,  1, -3),
					dNNNEE =LOAD(errors,  2, -3),
					dNNNEEE=LOAD(errors,  3, -3),
					dNNWWW =LOAD(errors, -3, -2),
					dNNWW  =LOAD(errors, -2, -2),
					dNNW   =LOAD(errors, -1, -2),
					dNN    =LOAD(errors,  0, -2),
					dNNE   =LOAD(errors,  1, -2),
					dNNEE  =LOAD(errors,  2, -2),
					dNNEEE =LOAD(errors,  3, -2),
					dNWWW  =LOAD(errors, -3, -1),
					dNWW   =LOAD(errors, -2, -1),
					dNW    =LOAD(errors, -1, -1),
					dN     =LOAD(errors,  0, -1),
					dNE    =LOAD(errors,  1, -1),
					dNEE   =LOAD(errors,  2, -1),
					dNEEE  =LOAD(errors,  3, -1),
					dWWW   =LOAD(errors, -3,  0),
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

				for(int k=0;k<G2_NPRED;++k)
				{
					//eNW + eN + eNE
					prederrs[k]=
						((unsigned)(kx-1)<(unsigned)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx-1]:0)+
						perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx]+
						((unsigned)(kx+1)<(unsigned)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx+1]:0);
				}
				
				//const int den2=2;
				const int correction=0;
				//int correction=((eN+eW)*5+(eNW+eNE)*3)/(16*4);
				int j=-1;
				int vmin, vmax;
#if 1
#define GRAD(pred, N, W, NW, vmin, vmax)\
				do\
				{\
					if(N<W)\
						vmin=N, vmax=W;\
					else\
						vmin=W, vmax=N;\
					pred=N+W-NW;\
					pred=CLAMP(vmin, pred, vmax);\
				}while(0)
#endif
				vmin=N, vmax=N;
				if(vmin>W)vmin=W;
				if(vmax<W)vmax=W;
#if 1
				//the 4 predictors from JPEG XL:
				++j, preds[j]=N-((eNW*params[G2_NPRED]+eN*params[G2_NPRED+1]+eNE*params[G2_NPRED+2]+(NN-N)*params[G2_NPRED+3]+(NW-W)*params[G2_NPRED+4])>>8);
				++j, preds[j]=W-((eN+eW+eNW)*params[G2_NPRED+5]>>8);
				++j, preds[j]=N-((eN+eW+eNE)*params[G2_NPRED+6]>>8);
				++j, preds[j]=W+NE-N;
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
				
				++j, preds[j]=clamp4(N+W -NW +correction, N, W, NW, NE);
				++j, preds[j]=clamp4(W+NE-N  +correction, N, W, NW, NE);
				++j, preds[j]=clamp4(N+NW-NNW+correction, N, W, NW, NE);
				++j, preds[j]=clamp4(N+NE-NNE+correction, N, W, NE, NEE);
				
				++j, preds[j]=(W+NEE)/2+correction;
				++j, preds[j]=NNNNNN+correction;
				++j, preds[j]=(NEEEE+NEEEEEE)/2+correction;
				++j, preds[j]=(WWWW+WWWWWW)/2+correction;
				//++j, preds[j]=W*3-WW*3+WWW;//parabolic
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

				++j, preds[j]=(N+W+NEEEEE+NEEEEEEE)/4+correction;
				++j, preds[j]=clamp4(N*2-NN+correction, N, W, NE, NEE);
				++j, preds[j]=(N+NNN)/2+correction;
				++j, preds[j]=((N+W)*3-NW*2)/4+correction;
#endif

				++j; GRAD(preds[j], N, W, NW, vmin, vmax);
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

				long long num=0;
				int weights[G2_NPRED], den=0;
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
				
				pred=CLAMP(vmin, pred, vmax);
				
				
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

				int idx=(src->iw*ky+kx)<<2|kc;

				//Secondary Symbol Estimation (SSE)
#if 1
				int dnb[]=//error neighbors
				{
					(int)dN   ,//order matters slightly
					(int)dW   ,
					(int)dNW  ,
					(int)dNE  ,
					(int)dNN  ,
					(int)dWW  ,
					(int)dNNWW,
					(int)dNNW ,
					(int)dNNE ,
					(int)dNNEE,
					(int)dNWW ,
					(int)dNEE ,
					
					//(unsigned char)dNNWW,
					//(unsigned char)dNNW ,
					//(unsigned char)dNN  ,
					//(unsigned char)dNNE ,
					//(unsigned char)dNNEE,
					//(unsigned char)dNWW ,
					//(unsigned char)dNW  ,
					//(unsigned char)dN   ,
					//(unsigned char)dNE  ,
					//(unsigned char)dNEE ,
					//(unsigned char)dWW  ,
					//(unsigned char)dW   ,
					
					//(unsigned char)dN     ,
					//(unsigned char)dW     ,
					//(unsigned char)dNW    ,
					//(unsigned char)dNE    ,
					//(unsigned char)dNN    ,
					//(unsigned char)dWW    ,
					//(unsigned char)dNNW   ,
					//(unsigned char)dNNE   ,
					//(unsigned char)dNWW   ,
					//(unsigned char)dNEE   ,
					//(unsigned char)dNNN   ,
					//(unsigned char)dWWW   ,
					//(unsigned char)dNNNW  ,
					//(unsigned char)dNNNE  ,
					//(unsigned char)dNNWW  ,
					//(unsigned char)dNNEE  ,
					//(unsigned char)dNEEE  ,
					//(unsigned char)dNWWW  ,
					//(unsigned char)dNNNWW ,
					//(unsigned char)dNNNEE ,
					//(unsigned char)dNNWWW ,
					//(unsigned char)dNNEEE ,
					//(unsigned char)dNNNWWW,
					//(unsigned char)dNNNEEE,
					
					//(unsigned char)dNNNWWW,
					//(unsigned char)dNNNWW ,
					//(unsigned char)dNNNW  ,
					//(unsigned char)dNNN   ,
					//(unsigned char)dNNNE  ,
					//(unsigned char)dNNNEE ,
					//(unsigned char)dNNNEEE,
					//(unsigned char)dNNWWW ,
					//(unsigned char)dNNWW  ,
					//(unsigned char)dNNW   ,
					//(unsigned char)dNN    ,
					//(unsigned char)dNNE   ,
					//(unsigned char)dNNEE  ,
					//(unsigned char)dNNEEE ,
					//(unsigned char)dNWWW  ,
					//(unsigned char)dNWW   ,
					//(unsigned char)dNW    ,
					//(unsigned char)dN     ,
					//(unsigned char)dNE    ,
					//(unsigned char)dNEE   ,
					//(unsigned char)dNEEE  ,
					//(unsigned char)dWWW   ,
					//(unsigned char)dWW    ,
					//(unsigned char)dW     ,
					
					//eNW,
					//eN ,
					//eNE,
					//eW ,
				};
				SseCell *pc[12];
				unsigned char hashval;

				for(int k2=3;k2>=0;--k2)
				{
					hashval=(pred+128)>>(8+6)&3;
					for(int k=0;k<3;++k)
					{
						hashval<<=3;
						hashval+=(int)(dnb[3*k2+k]*0xFC28)>>(16+6)&3;
					}

					pc[k2]=g2_SSE[k2]+hashval;

					//int c=2, d=pc[k2]->sum8+pc[k2]->sum7;
					//if(pc[k2]->count)
					//{
					//	++c;
					//	d+=(int)((((long long)pc[k2]->sum<<8)+(pc[k2]->count>>1))/pc[k2]->count);
					//}
					//pred+=d/c;
					pred+=pc[k2]->count?(int)((((long long)pc[k2]->sum<<8)+(pc[k2]->count>>1))/pc[k2]->count):0;
					pred=CLAMP(-(nlevels<<7), pred, (nlevels<<7)-1);
				}

				char delta;
				if(fwd)
				{
					delta=src->data[idx]-((pred+128)>>8);
					if(enable_ma)
					{
						delta+=nlevels>>1;
						delta&=nlevels-1;
						delta-=nlevels>>1;
					}
					b2[idx]=delta;
				}
				else
				{
					delta=src->data[idx];
					int pixel=delta+((pred+128)>>8);
					if(enable_ma)
					{
						pixel+=nlevels>>1;
						pixel&=nlevels-1;
						pixel-=nlevels>>1;
					}
					b2[idx]=pixel;
				}
				for(int k=0;k<4;++k)
				{
					if(pc[k]->count+1>640)
					{
						pc[k]->count>>=1;
						pc[k]->sum>>=1;
					}
					++pc[k]->count;
					pc[k]->sum+=delta;
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

				char delta;
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

				char delta;
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
				char delta;
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
				char delta;
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
				
				char delta;//modular arithmetic
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
				b2[idx]=buf[idx]+((spred+128)>>8);
#endif

				int curr=pixels[idx]<<8;
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
	memcpy(src->data, b2, res*sizeof(int[4]));
	free(perrors);
	free(b2);
	//free(SSE_count);
	//free(SSE_sum);
}




//DWTs

ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	ArrayHandle sizes;
	int lw=floor_log2(iw), lh=floor_log2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	int nstages=0;
	DWTSize *p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=iw;
	p->h=ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=w2;
		p->h=h2;
	}
	return sizes;
}

#if 0
void dwt2d_grad_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	//for(int it=sizes_start;it<1;++it)
	{
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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

//lifting-based 8bit Haar DWT
void dwt1d_haar_fwd(int *buffer, int count, int stride, int *b2)
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
void dwt1d_haar_inv(int *buffer, int count, int stride, int *b2)
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit squeeze DWT
char smoothtendency(char B, char a, char n)
{
	char diff=0;
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
void dwt1d_squeeze_fwd(int *buffer, int count, int stride, int *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		char
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
void dwt1d_squeeze_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(int ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
	{
		char
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt1d_cdf53_fwd(int *buffer, int count, int stride, int *b2)
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
void dwt1d_cdf53_inv(int *buffer, int count, int stride, int *b2)
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
static void dwt1d_u8_scale(int *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=buf[k]*coeff>>16;
}
static void dwt1d_u8_unscale(int *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=(buf[k]<<16)/coeff;
}
void dwt1d_cdf97_fwd(int *buffer, int count, int stride, int *b2)
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
void dwt1d_cdf97_inv(int *buffer, int count, int stride, int *b2)
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
		LOG_ERROR("Allocation error");
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
		LOG_ERROR("Allocation error");
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
		LOG_ERROR("Allocation error");
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


//DCTs
static void dct4_fwd_i8(char *x)
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
static void dct4_inv_i8(char *x)
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
void image_dct4_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[ kx>>2           ]=x[0];
				temp[(kx>>2)+(iw>>2)  ]=x[1];
				temp[(kx>>2)+(iw>>2)*2]=x[2];
				temp[(kx>>2)+(iw>>2)*3]=x[3];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[(ky>>2)          ]=x[0];
				temp[(ky>>2)+(ih>>2)  ]=x[1];
				temp[(ky>>2)+(ih>>2)*2]=x[2];
				temp[(ky>>2)+(ih>>2)*3]=x[3];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct4_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>2)          ],
					temp[(ky>>2)+(ih>>2)  ],
					temp[(ky>>2)+(ih>>2)*2],
					temp[(ky>>2)+(ih>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>2)          ],
					temp[(kx>>2)+(iw>>2)  ],
					temp[(kx>>2)+(iw>>2)*2],
					temp[(kx>>2)+(iw>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
}

static void dct8_fwd_i8(char *x)
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
static void dct8_inv_i8(char *x)
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
void image_dct8_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
					image[(idx+4)<<2|kc],
					image[(idx+5)<<2|kc],
					image[(idx+6)<<2|kc],
					image[(idx+7)<<2|kc],
				};

				//char y[8];
				//memcpy(y, x, 8);
				//dct8_fwd_i8(y);
				//dct8_inv_i8(y);
				//if(memcmp(x, y, 8))
				//	x[0]=y[0];

				dct8_fwd_i8(x);

				temp[ kx>>3           ]=x[0];
				temp[(kx>>3)+(iw>>3)  ]=x[1];
				temp[(kx>>3)+(iw>>3)*2]=x[2];
				temp[(kx>>3)+(iw>>3)*3]=x[3];
				temp[(kx>>3)+(iw>>3)*4]=x[4];
				temp[(kx>>3)+(iw>>3)*5]=x[5];
				temp[(kx>>3)+(iw>>3)*6]=x[6];
				temp[(kx>>3)+(iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
					image[(idx+iw*4)<<2|kc],
					image[(idx+iw*5)<<2|kc],
					image[(idx+iw*6)<<2|kc],
					image[(idx+iw*7)<<2|kc],
				};

				dct8_fwd_i8(x);

				temp[(ky>>3)          ]=x[0];
				temp[(ky>>3)+(ih>>3)  ]=x[1];
				temp[(ky>>3)+(ih>>3)*2]=x[2];
				temp[(ky>>3)+(ih>>3)*3]=x[3];
				temp[(ky>>3)+(ih>>3)*4]=x[4];
				temp[(ky>>3)+(ih>>3)*5]=x[5];
				temp[(ky>>3)+(ih>>3)*6]=x[6];
				temp[(ky>>3)+(ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct8_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>3)          ],
					temp[(ky>>3)+(ih>>3)  ],
					temp[(ky>>3)+(ih>>3)*2],
					temp[(ky>>3)+(ih>>3)*3],
					temp[(ky>>3)+(ih>>3)*4],
					temp[(ky>>3)+(ih>>3)*5],
					temp[(ky>>3)+(ih>>3)*6],
					temp[(ky>>3)+(ih>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
				image[(idx+iw*4)<<2|kc]=x[4];
				image[(idx+iw*5)<<2|kc]=x[5];
				image[(idx+iw*6)<<2|kc]=x[6];
				image[(idx+iw*7)<<2|kc]=x[7];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>3)          ],
					temp[(kx>>3)+(iw>>3)  ],
					temp[(kx>>3)+(iw>>3)*2],
					temp[(kx>>3)+(iw>>3)*3],
					temp[(kx>>3)+(iw>>3)*4],
					temp[(kx>>3)+(iw>>3)*5],
					temp[(kx>>3)+(iw>>3)*6],
					temp[(kx>>3)+(iw>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
				image[(idx+4)<<2|kc]=x[4];
				image[(idx+5)<<2|kc]=x[5];
				image[(idx+6)<<2|kc]=x[6];
				image[(idx+7)<<2|kc]=x[7];
			}
		}
#endif
	}
	free(temp);
}


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
int  predict_dct3(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen, const float *dct3, const float *dct4)
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
void pred_dct3_fwd(char *buf, int iw, int ih, int nch, int bytestride)
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
void pred_dct3_inv(char *buf, int iw, int ih, int nch, int bytestride)
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
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, image, (size_t)iw*ih<<2);
	//int maxdim=MAXVAR(iw, ih);
	//char *temp=(char*)malloc(maxdim);
	//if(!temp)
	//{
	//	LOG_ERROR("Allocation error");
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
		LOG_ERROR("Allocation error");
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
		LOG_ERROR("Allocation error");
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
		LOG_ERROR("Allocation error");
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