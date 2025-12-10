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
	//{
	//	int k=0;
	//	while(!scanf(" %d", &k));
	//}
	exit(1);
}
#define LOG_ERROR(MSG, ...) log_error(__FILE__, __LINE__, MSG, ##__VA_ARGS__)
static unsigned char* ppm_load(const char *fn, int *ret_iw, int *ret_ih)
{
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		LOG_ERROR("Cannot open \"%s\"", fn);
		return 0;
	}
	int tag=0;
	fread(&tag, 1, 2, fsrc);
	if(tag!=('P'|'6'<<8))
	{
		LOG_ERROR("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	int iw=0, ih=0;
	if(fgetc(fsrc)!='\n')
	{
		LOG_ERROR("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	int nscan=fscanf(fsrc, "%d %d\n", &iw, &ih);
	if(nscan!=2)
	{
		LOG_ERROR("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	int vmax=0;
	nscan=fscanf(fsrc, "%d\n", &vmax);
	if(nscan!=1||vmax!=255)
	{
		LOG_ERROR("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *buf=(unsigned char*)malloc(size+16);
	if(!buf)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	ptrdiff_t nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
	{
		LOG_ERROR("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	fclose(fsrc);
	return buf;
}
int main(int argc, char **argv)
{
	const char *fn1=0, *fn2=0;

	if(argc!=3)
	{
		printf("Usage:  \"%s\"  image1.ppm  image2.ppm\n", argv[0]);
		printf("  Prints RMSE & PSNR\n");
		return 1;
	}
	fn1=argv[1];
	fn2=argv[2];

	int iw=0, ih=0, iw2=0, ih2=0;
	unsigned char *image1=ppm_load(fn1, &iw, &ih);
	unsigned char *image2=ppm_load(fn2, &iw2, &ih2);
	if(iw!=iw2||ih!=ih2)
	{
		LOG_ERROR(
			"Dimension mismatch  im1 %d*%d vs %d*%d"
			, iw, ih
			, iw2, ih2
		);
	}
	long long errors[4]={0};
	ptrdiff_t idx=-1;
	const unsigned char *ptr1=image1, *ptr2=image2;
	for(ptrdiff_t k=0, res=(ptrdiff_t)3*iw*ih;k<res;k+=3)
	{
		int a0=ptr1[0], b0=ptr2[0];
		int a1=ptr1[1], b1=ptr2[1];
		int a2=ptr1[2], b2=ptr2[2];
		if(idx==-1&&((a0^b0)|(a1^b1)|(a2^b2)))
			idx=k;
		int delta;
		delta=a0-b0; errors[0]+=delta*delta;
		delta=a1-b1; errors[1]+=delta*delta;
		delta=a2-b2; errors[2]+=delta*delta;
		ptr1+=3;
		ptr2+=3;
	}
	errors[3]=errors[0]+errors[1]+errors[2];
	double invres=1/((double)iw*ih);
	double rmse[]=
	{
		sqrt((double)errors[0]*invres),
		sqrt((double)errors[1]*invres),
		sqrt((double)errors[2]*invres),
		sqrt((double)errors[3]*invres*(1./3)),
	};
	double psnr[]=
	{
		20*log10(255/rmse[0]),
		20*log10(255/rmse[1]),
		20*log10(255/rmse[2]),
		20*log10(255/rmse[3]),
	};
	printf(
		"RMSE, PSNR\n"
		"T %12.6lf/255  %12.6lf dB\n"
		"R %12.6lf/255  %12.6lf dB\n"
		"G %12.6lf/255  %12.6lf dB\n"
		"B %12.6lf/255  %12.6lf dB\n"
		, rmse[3], psnr[3]
		, rmse[0], psnr[0]
		, rmse[1], psnr[1]
		, rmse[2], psnr[2]
	);
	if(idx!=-1)
	{
		ptrdiff_t idx2=idx;
		int kc, kx, ky;

		kc=(int)(idx2%3);
		idx2/=3;
		kx=(int)(idx2%iw);
		idx2/=iw;
		ky=(int)idx2;
		printf("First diff at XYC %d %d %d\n", kx, ky, kc);
	}
	else
		printf("Perfect match\n");
	free(image1);
	free(image2);
	return 0;
}
