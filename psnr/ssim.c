#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
void log_error(const char *file, int line, const char *msg, ...)
{
	printf("\n\n%s(%d): ", file, line);
	if(msg)
	{
		va_list args;
		va_start(args, msg);
		vprintf(msg, args);
		va_end(args);
	}
	{
		int k=0;
		while(!scanf(" %d", &k));
	}
	exit(1);
}
#define LOG_ERROR(MSG, ...) log_error(__FILE__, __LINE__, MSG, ##__VA_ARGS__)
static unsigned char* load_ppm(const char *fn, int *ret_iw, int *ret_ih)
{
	int iw=0, ih=0;
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		LOG_ERROR("Cannot open \"%s\"", fn);
		return 0;
	}
	int c=0;
	ptrdiff_t nread=fread(&c, 1, 2, fsrc);
	if(nread!=2||c!=('P'|'6'<<8))
	{
		LOG_ERROR("Inalid file \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		LOG_ERROR("Invalid PPM file");
		return 0;
	}
	c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	iw=0;
	while((unsigned)(c-'0')<10)
	{
		iw=10*iw+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	ih=0;
	while((unsigned)(c-'0')<10)
	{
		ih=10*ih+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
	{
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *image=(unsigned char*)malloc(size);
	if(!image)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	nread=fread(image, 1, size, fsrc);
	if(nread!=size)
	{
		LOG_ERROR("Truncated file  expected %td  got %td", size, nread);
		free(image);
		return 0;
	}
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	return image;
}
static void measure_ssim_ppm(const char *fn0, const char *fn1, double *ret_ssim)
{
	int iw=0, ih=0;
	unsigned char *im1=0, *im2=0;
	{
		int w2=0, h2=0;
		im1=load_ppm(fn0, &iw, &ih);
		im2=load_ppm(fn1, &w2, &h2);
		if(!im1||!im2||iw!=w2||ih!=h2)
		{
			LOG_ERROR("Dimension mismatch %d*%d vs %d*%d", iw, ih, w2, h2);
			return;
		}
	}
	//printf("WH %d*%d\n", iw, ih);
	ptrdiff_t count=0;
	double ssim[3]={0};
	double c1=0.01*255, c2=0.03*255;
	c1*=c1;
	c2*=c2;
#define SSIM_SIZE 11
	for(int ky=SSIM_SIZE/2;ky<ih-SSIM_SIZE/2;++ky)
	{
		for(int kx=SSIM_SIZE/2;kx<iw-SSIM_SIZE/2;++kx)
		{
			int idx=3*(iw*ky+kx);
			unsigned char *p1=im1+idx, *p2=im2+idx;
			short patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];
			for(int ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(int kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=3*(iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;
					patch1[0][dstidx]=p1[srcidx+0]-128;
					patch2[0][dstidx]=p2[srcidx+0]-128;
					patch1[1][dstidx]=p1[srcidx+1]-128;
					patch2[1][dstidx]=p2[srcidx+1]-128;
					patch1[2][dstidx]=p1[srcidx+2]-128;
					patch2[2][dstidx]=p2[srcidx+2]-128;
				}
			}
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
	free(im1);
	free(im2);
	if(count)
	{
		ssim[0]/=count;
		ssim[1]/=count;
		ssim[2]/=count;
	}
	ret_ssim[0]=ssim[0];
	ret_ssim[1]=ssim[1];
	ret_ssim[2]=ssim[2];
}
int main(int argc, char **argv)
{
	const char *fn1=0, *fn2=0;

	if(argc!=3)
	{
		printf("Usage:  \"%s\"  image1.ppm  image2.ppm\n", argv[0]);
		printf("  Prints SSIM\n");
		return 1;
	}
	fn1=argv[1];
	fn2=argv[2];

	double ssim[4]={0};
	measure_ssim_ppm(fn1, fn2, ssim);
	ssim[3]=(6*ssim[0]+ssim[1]+ssim[2])/8;
	printf("SSIM TYUV %12.9lf %12.9lf %12.9lf %12.9lf\n", ssim[3], ssim[0], ssim[1], ssim[2]);
	//printf("SSIM%11.8lf\n", ssim);
	return 0;
}
