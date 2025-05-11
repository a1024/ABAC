#include"ebench.h"
#include<immintrin.h>
static const char file[]=__FILE__;

#define SSIM_SIZE 11

void measure_ssim_avg(const Image *im0, const Image *im1, double *ret_ssim)
{
	if(!im0||!im1||im0->iw!=im1->iw||im0->ih!=im1->ih)
	{
		if(loud_transforms&&im0&&im1)
			messagebox(MBOX_OK, "Error", "Dimension mismatch  im0 %d*%d  im1 %d*%d", im0->iw, im0->ih, im1->iw, im1->ih);
		return;
	}
	ptrdiff_t count=0;
	double ssim[4]={0};
	double c1=0.01*255, c2=0.03*255;
	c1*=c1;
	c2*=c2;
	for(int ky=SSIM_SIZE/2;ky<im0->ih-SSIM_SIZE/2;ky+=5)
	{
		for(int kx=SSIM_SIZE/2;kx<im0->iw-SSIM_SIZE/2;kx+=5)
		{
			int idx=4*(im0->iw*ky+kx);
			const int *p1=im0->data+idx, *p2=im1->data+idx;
			short patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];
			for(int ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(int kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=4*(im0->iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;
					patch1[0][dstidx]=p1[srcidx+0];
					patch2[0][dstidx]=p2[srcidx+0];
					patch1[1][dstidx]=p1[srcidx+1];
					patch2[1][dstidx]=p2[srcidx+1];
					patch1[2][dstidx]=p1[srcidx+2];
					patch2[2][dstidx]=p2[srcidx+2];
				}
			}
#if 0
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)//RGB->YUV
			{
				int r1=patch1[0][k], g1=patch1[1][k], b1=patch1[2][k];
				int r2=patch2[0][k], g2=patch2[1][k], b2=patch2[2][k];

				r1-=g1;
				b1-=g1;
				g1+=(r1+b1)>>2;
				r2-=g2;
				b2-=g2;
				g2+=(r2+b2)>>2;

				patch1[0][k]=g1;
				patch1[1][k]=b1;
				patch1[2][k]=r1;
				patch2[0][k]=g2;
				patch2[1][k]=b2;
				patch2[2][k]=r2;
			}
#endif
			for(int kc=0;kc<3;++kc)
			{
				double mean1=0, mean2=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[kc][k];
					mean2+=patch2[kc][k];
				}
				mean1/=SSIM_SIZE*SSIM_SIZE;
				mean2/=SSIM_SIZE*SSIM_SIZE;
				double var1=0, var2=0, cov=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
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
				double curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));
				ssim[kc]+=curr_ssim;
			}
			++count;
		}
	}
	if(count)
	{
		ssim[0]/=count;
		ssim[1]/=count;
		ssim[2]/=count;
	}
	ssim[3]=(ssim[0]+ssim[1]+ssim[2])/3;
	ret_ssim[0]=ssim[0];
	ret_ssim[1]=ssim[1];
	ret_ssim[2]=ssim[2];
	ret_ssim[3]=ssim[3];
}
void measure_ssim_map(const Image *im0, Image *im1)
{
	if(!im0||!im1||im0->iw!=im1->iw||im0->ih!=im1->ih)
	{
		if(loud_transforms&&im0&&im1)
			messagebox(MBOX_OK, "Error", "Dimension mismatch  im0 %d*%d  im1 %d*%d", im0->iw, im0->ih, im1->iw, im1->ih);
		return;
	}
	int nlevels[]=
	{
		1<<im1->depth[0],
		1<<im1->depth[1],
		1<<im1->depth[2],
		1<<im1->depth[3],
	};
	Image *im2=0;
	image_copy(&im2, im1);
	if(!im2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	ptrdiff_t count=0;
	double ssim[4]={0};
	double c1=0.01*255, c2=0.03*255;
	c1*=c1;
	c2*=c2;
	//float *dstbuf=(float*)im1->data;
	for(int ky=SSIM_SIZE/2;ky<im0->ih-SSIM_SIZE/2;++ky)
	{
		for(int kx=SSIM_SIZE/2;kx<im0->iw-SSIM_SIZE/2;++kx)
		{
			int idx=4*(im0->iw*ky+kx);
			const int *p1=im0->data+idx, *p2=im2->data+idx;
			short patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];
			for(int ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(int kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=4*(im0->iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;
					patch1[0][dstidx]=p1[srcidx+0];
					patch2[0][dstidx]=p2[srcidx+0];
					patch1[1][dstidx]=p1[srcidx+1];
					patch2[1][dstidx]=p2[srcidx+1];
					patch1[2][dstidx]=p1[srcidx+2];
					patch2[2][dstidx]=p2[srcidx+2];
				}
			}
#if 0
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)//RGB->YUV
			{
				int r1=patch1[0][k], g1=patch1[1][k], b1=patch1[2][k];
				int r2=patch2[0][k], g2=patch2[1][k], b2=patch2[2][k];

				r1-=g1;
				b1-=g1;
				g1+=(r1+b1)>>2;
				r2-=g2;
				b2-=g2;
				g2+=(r2+b2)>>2;

				patch1[0][k]=g1;
				patch1[1][k]=b1;
				patch1[2][k]=r1;
				patch2[0][k]=g2;
				patch2[1][k]=b2;
				patch2[2][k]=r2;
			}
#endif
			for(int kc=0;kc<3;++kc)
			{
				double mean1=0, mean2=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[kc][k];
					mean2+=patch2[kc][k];
				}
				mean1/=SSIM_SIZE*SSIM_SIZE;
				mean2/=SSIM_SIZE*SSIM_SIZE;
				double var1=0, var2=0, cov=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
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
				double curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));

				const double brightness=8;
				double val=(1-curr_ssim)*brightness*(nlevels[kc]-1)-(nlevels[kc]>>1);
				CLAMP2(val, -(nlevels[kc]>>1), (nlevels[kc]>>1)-1);
				im1->data[4*(im1->iw*ky+kx)+kc]=(int)val;
			}
			++count;
		}
	}
	//for(int kc=0;kc<3;++kc)
	//{
	//	float vmin=0, vmax=0;
	//	for(int ky=SSIM_SIZE/2, it=0;ky<im0->ih-SSIM_SIZE/2;++ky)
	//	{
	//		for(int kx=SSIM_SIZE/2;kx<im0->iw-SSIM_SIZE/2;++kx, ++it)
	//		{
	//			float val=dstbuf[4*(im1->iw*ky+kx)+kx];
	//			if(!it||vmin>val)
	//				vmin=val;
	//			if(!it||vmax<val)
	//				vmax=val;
	//		}
	//	}
	//}
	free(im2);
}