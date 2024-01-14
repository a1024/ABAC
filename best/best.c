#include"best.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

#define ENCODE t45_encode
#define DECODE t45_decode

void compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward)
{
	ptrdiff_t len=(ptrdiff_t)bytestride*iw*ih;
	int inc=bytestride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-bytestride:0;k>=0&&k<len;k+=inc)
	{
		if(memcmp(b1+k, b0+k, symbytes))
		{
			ptrdiff_t idx=k/bytestride, kx=idx%iw, ky=idx/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
			for(int kc=0;kc<symbytes;++kc)
				printf("C%d  0x%02X != 0x%02X    %d != %d\n", kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}
const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
	"pgm",
	"pnm",
};
void batch_test(const char *path)
{
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Start %s\n", g_buf);
	double t_start=time_sec();
	ArrayHandle filenames=get_filenames(path, g_extensions, COUNTOF(g_extensions), 1);
	if(!filenames)
	{
		printf("No images in \'%s\'\n", path);
		return;
	}
	long long
		count_PNG=0, count_JPEG=0,
		sum_cPNGsize=0, sum_cJPEGsize=0,
		sum_uPNGsize=0, sum_uJPEGsize=0,
		sum_testsize=0;

	for(ptrdiff_t k=0;k<(ptrdiff_t)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);

		if(!fn)
		{
			LOG_ERROR("filename read error");
			continue;
		}

		ptrdiff_t formatsize=get_filesize(fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images
			continue;

		int iw=0, ih=0, nch0=3, stride=4;
		long long cycles=__rdtsc();
		unsigned char *buf=stbi_load(fn[0]->data, &iw, &ih, 0, 4);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}

		ptrdiff_t res=(ptrdiff_t)iw*ih, len=res*stride, usize=res*nch0;
		double ratio=(double)usize/formatsize;
		printf("%3lld/%3lld  \"%s\"\tCR %lf (%lf BPP) Dec %lf CPB", k+1, filenames->count, fn[0]->data, ratio, 8/ratio, (double)cycles/usize);
		if(!acme_stricmp(fn[0]->data+fn[0]->count-3, "PNG"))
		{
			sum_cPNGsize+=formatsize;
			sum_uPNGsize+=usize;
			++count_PNG;
		}
		else//assumed
		{
			sum_cJPEGsize+=formatsize;
			sum_uJPEGsize+=usize;
			++count_JPEG;
		}
		unsigned char *b2=(unsigned char*)malloc(len);
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memset(b2, 0, len);
		
		//T34+: ABAC + adaptive Bayesian inference
#if 1
		{
			ArrayHandle cdata=0;
			printf("\n");

			ENCODE(buf, iw, ih, &cdata, 1);
			DECODE(cdata->data, cdata->count, iw, ih, b2, 1);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T45", 0);

			printf("\n");
		}
#endif
		free(buf);
		free(b2);
	}
	printf("Batch elapsed ");
	timedelta2str(0, 0, time_sec()-t_start);
	printf("\n");
	ptrdiff_t totalusize=sum_uPNGsize+sum_uJPEGsize;
	if(totalusize)
	{
		printf("\nOn average:\n");
		printf("BMP     csize %9lld\n", totalusize);
		if(sum_cPNGsize)
			printf("PNG     csize %9lld  CR %lf  (%lld images)\n", sum_cPNGsize, (double)sum_uPNGsize/sum_cPNGsize, count_PNG);
		if(sum_cJPEGsize)
			printf("JPEG    csize %9lld  CR %lf  (%lld images)\n", sum_cJPEGsize, (double)sum_uJPEGsize/sum_cJPEGsize, count_JPEG);
		printf("test    csize %9lld  CR %lf\n", sum_testsize, (double)totalusize/sum_testsize);
	}
	else
		printf("\nNo valid images found\n");

	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);

	printf("\nDone.\n");
	pause();
}
int main(int argc, char **argv)
{
	printf("Best  %s  %s\n", __DATE__, __TIME__);
#if 1
	long long cycles;
	int iw=0, ih=0, nch0=3,
		nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	const char *fn=0;
#ifdef _DEBUG
	//fn="C:/Projects/datasets/CLIC11-crop4-2.PNG";
	//fn="C:/Projects/datasets/CLIC11-small4.PNG";
	//fn="C:/Projects/datasets/dataset-CLIC30/11.png";
	//fn="C:/Projects/datasets/dataset-kodak";
	fn="C:/Projects/datasets/dataset-kodak/kodim13.png";

	//fn="D:/ML/dataset-CLIC30";
	//fn="D:/ML/dataset-kodak";
	//fn="D:/ML/dataset-CLIC30/16.png";//hardest noiseless CLIC30 image
	//fn="D:/ML/dataset-CLIC30/17.png";
	//fn="D:/ML/dataset-kodak/kodim13.png";
	//fn="D:/ML/dataset-kodak-small/13.PNG";
#endif
	if(fn||argc==2)
	{
		if(!fn)
			fn=argv[1];
		ptrdiff_t formatsize=get_filesize(fn);
		if(formatsize==-1)
		{
			LOG_ERROR("Cannot open \"%s\"", fn);
			return 0;
		}
		if(!formatsize)//path
		{
			batch_test(fn);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		cycles=__rdtsc();
		buf=stbi_load(fn, &iw, &ih, 0, 4);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			LOG_ERROR("Couldn't open \"%s\"", fn);
			return 0;
		}
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		printf("Format Dec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", (double)cycles/(resolution*nch0), iw, ih, nch0, formatsize, (double)resolution*nch0/formatsize);
	}
	else if(argc==3)
	{
		const char *fn1=argv[1], *fn2=argv[2];
		int w2, h2;
		buf=stbi_load(fn1, &iw, &ih, 0, 4);
		b2 =stbi_load(fn2, &w2, &h2, 0, 4);
		if(!buf)
		{
			printf("Couldn't open %s\n", fn1);
			return 1;
		}
		if(!b2)
		{
			printf("Couldn't open %s\n", fn2);
			return 1;
		}
		if(iw!=w2||ih!=h2)
		{
			printf("Expected two images of SAME RESOLUTION. %dx%d != %dx%d\n", iw, ih, w2, h2);
			return 1;
		}
		ptrdiff_t formatsize=get_filesize(fn2);
		int res=iw*ih;
		long long sum[3]={0};
		for(int k=0;k<res;++k)
		{
			int dr=buf[k<<2  ]-b2[k<<2  ],
				dg=buf[k<<2|1]-b2[k<<2|1],
				db=buf[k<<2|2]-b2[k<<2|2];
			sum[0]+=dr*dr;
			sum[1]+=dg*dg;
			sum[2]+=db*db;
		}
		double rmse[]=
		{
			sqrt((double)sum[0]/res),
			sqrt((double)sum[1]/res),
			sqrt((double)sum[2]/res),
			sqrt((double)(sum[0]+sum[1]+sum[2])/(res*3)),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		double CR=res*3./formatsize;
		printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)formatsize, CR, 8/CR);
		printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
		printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
		printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
		return 0;
	}
	else
	{
		printf(
			"Usage:\n"
			"[1] best.exe  image\n"
			"    test the algorithm with file\n"
			"\n"
			"[2] best.exe  path\n"
			"    batch-test the algorithm with all images in a folder\n"
			"\n"
			"[3] best.exe  image1  image2\n"
			"    RMSE & PSNR measurement (where image1 is ground truth)\n"
		);
		pause();
		return 0;
	}

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<len;k+=nch)
			buf[k]=0xFF;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	size_t usize=len*nch0>>2;

	printf("\n");
	
	int loud=0;
	
	ArrayHandle cdata=0;

	ENCODE(buf, iw, ih, &cdata, 1);
	DECODE(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T45", 0);
	memset(b2, 0, len);
	printf("\n");


	free(b2);
	free(buf);
#endif
	printf("Done.\n");
	pause();
	return 0;
}
