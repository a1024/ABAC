#include"ebench.h"
#if defined _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
enum
{
	SSIM_SIZE=11
};
void measure_ssim_avg(const Image *ima, const Image *imb, double *ret_ssim)
{
	intptr_t count=0;
	double c1=0.01*255, c2=0.03*255, ssim[4]={0};
	const int *p1, *p2;
	int k, ky, kx, kc, ky2, kx2, idx;

	if(!ima||!imb||ima->iw!=imb->iw||ima->ih!=imb->ih)
	{
		if(loud_transforms&&ima&&imb)
			messagebox(MBOX_OK, "Error", "Dimension mismatch  ima %d*%d  imb %d*%d", ima->iw, ima->ih, imb->iw, imb->ih);
		return;
	}
	c1*=c1;
	c2*=c2;
	for(ky=SSIM_SIZE/2;ky<ima->ih-SSIM_SIZE/2;ky+=5)
	{
		for(kx=SSIM_SIZE/2;kx<ima->iw-SSIM_SIZE/2;kx+=5)
		{
			int16_t patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];

			idx=4*(ima->iw*ky+kx);
			p1=((const int32_t*)HDR_DATC(ima))+idx;
			p2=((const int32_t*)HDR_DATC(imb))+idx;
			for(ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=4*(ima->iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;

					patch1[0][dstidx]=p1[srcidx+0];
					patch2[0][dstidx]=p2[srcidx+0];
					patch1[1][dstidx]=p1[srcidx+1];
					patch2[1][dstidx]=p2[srcidx+1];
					patch1[2][dstidx]=p1[srcidx+2];
					patch2[2][dstidx]=p2[srcidx+2];
				}
			}
			for(kc=0;kc<3;++kc)
			{
				double mean1=0, mean2=0;
				double var1=0, var2=0, cov=0, curr_ssim;

				for(k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[kc][k];
					mean2+=patch2[kc][k];
				}
				mean1/=SSIM_SIZE*SSIM_SIZE;
				mean2/=SSIM_SIZE*SSIM_SIZE;
				for(k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					double d1=patch1[kc][k]-mean1;
					double d2=patch2[kc][k]-mean2;
					var1+=d1*d1;
					var2+=d2*d2;
					cov+=d1*d2;
				}
				var1/=SSIM_SIZE*SSIM_SIZE;
				var2/=SSIM_SIZE*SSIM_SIZE;
				cov/=SSIM_SIZE*SSIM_SIZE;
				curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));
				ssim[kc]+=curr_ssim;
			}
			++count;
		}
	}
	if(count)
	{
		double inv=1./count;

		ssim[0]*=inv;
		ssim[1]*=inv;
		ssim[2]*=inv;
	}
	ssim[3]=(ssim[0]+ssim[1]+ssim[2])/3;
	ret_ssim[0]=ssim[0];
	ret_ssim[1]=ssim[1];
	ret_ssim[2]=ssim[2];
	ret_ssim[3]=ssim[3];
}
void measure_ssim_map(const Image *ima, Image *imb)
{
	double c1=0.01*255, c2=0.03*255;
	Image *im2=0;
	const int *p1, *p2;
	int nlevels[4], k, ky, kx, kc, idx, ky2, kx2;

	nlevels[0]=1<<imb->depth[0];
	nlevels[1]=1<<imb->depth[1];
	nlevels[2]=1<<imb->depth[2];
	nlevels[3]=1<<imb->depth[3];
	if(!ima||!imb||ima->iw!=imb->iw||ima->ih!=imb->ih)
	{
		if(loud_transforms&&ima&&imb)
			messagebox(MBOX_OK, "Error", "Dimension mismatch  ima %d*%d  imb %d*%d", ima->iw, ima->ih, imb->iw, imb->ih);
		return;
	}
	image_copy(&im2, imb);
	if(!im2)
	{
		CRASH("Alloc error");
		return;
	}
	c1*=c1;
	c2*=c2;
	for(ky=SSIM_SIZE/2;ky<ima->ih-SSIM_SIZE/2;++ky)
	{
		for(kx=SSIM_SIZE/2;kx<ima->iw-SSIM_SIZE/2;++kx)
		{
			int16_t patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];

			idx=4*(ima->iw*ky+kx);
			p1=((const int32_t*)HDR_DATC(ima))+idx;
			p2=((const int32_t*)HDR_DATC(im2))+idx;
			for(ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=4*(ima->iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;

					patch1[0][dstidx]=p1[srcidx+0];
					patch2[0][dstidx]=p2[srcidx+0];
					patch1[1][dstidx]=p1[srcidx+1];
					patch2[1][dstidx]=p2[srcidx+1];
					patch1[2][dstidx]=p1[srcidx+2];
					patch2[2][dstidx]=p2[srcidx+2];
				}
			}
			for(kc=0;kc<3;++kc)
			{
				const double brightness=8;
				double mean1=0, mean2=0, var1=0, var2=0, cov=0;
				double curr_ssim, val;

				for(k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[kc][k];
					mean2+=patch2[kc][k];
				}
				mean1/=SSIM_SIZE*SSIM_SIZE;
				mean2/=SSIM_SIZE*SSIM_SIZE;
				for(k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					double d1=patch1[kc][k]-mean1;
					double d2=patch2[kc][k]-mean2;
					var1+=d1*d1;
					var2+=d2*d2;
					cov+=d1*d2;
				}
				var1/=SSIM_SIZE*SSIM_SIZE;
				var2/=SSIM_SIZE*SSIM_SIZE;
				cov/=SSIM_SIZE*SSIM_SIZE;
				curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));

				val=(1-curr_ssim)*brightness*(nlevels[kc]-1)-(nlevels[kc]>>1);
				CLAMP2(val, -(nlevels[kc]>>1), (nlevels[kc]>>1)-1);
				((int32_t*)HDR_DATA(imb))[4*(imb->iw*ky+kx)+kc]=(int)val;
			}
		}
	}
	free(im2);
}