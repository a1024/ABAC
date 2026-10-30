#include<stdint.h>
#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
void crash(const char *file, int line, const char *msg, ...)
{
	printf("\n\n%s(%d): ", file, line);
	if(msg)
	{
		va_list args;
		va_start(args, msg);
		vprintf(msg, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(MSG, ...) crash(__FILE__, __LINE__, MSG, ##__VA_ARGS__)
static uint8_t* pnm_load(const char *fn, int *ret_iw, int *ret_ih, int *ret_nch)
{
	int iw=0, ih=0, nch=0, c=0;
	ptrdiff_t size, nread;
	FILE *fsrc;

	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", fn);
		return 0;
	}
	c=getc(fsrc);
	c|=getc(fsrc)<<8;
	if(c!=('P'|'6'<<8)&&c!=('P'|'5'<<8))
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	nch=c==('P'|'6'<<8)?3:1;
	c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
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
	iw=0;
	while((uint32_t)(c-'0')<10)
	{
		iw=10*iw+c-'0';
		c=fgetc(fsrc);
	}
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	ih=0;
	while((uint32_t)(c-'0')<10)
	{
		ih=10*ih+c-'0';
		c=fgetc(fsrc);
	}
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	c|=fgetc(fsrc)<<8;
	c|=fgetc(fsrc)<<16;
	if(c!=('2'|'5'<<8|'5'<<16))
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	size=(ptrdiff_t)nch*iw*ih;
	uint8_t *buf=(uint8_t*)malloc(size+32);
	if(!buf)
	{
		CRASH("Alloc error");
		return 0;
	}
	memset(buf, 0, size+32);
	nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
	{
		printf(
			"Truncated \"%s\"\n"
			"CWH %d*%d*%d\n"
			"Requested %12lld bytes\n"
			"Read      %12lld bytes\n"
			"\n"
			, fn
			, nch, iw, ih
			, (int64_t)size
			, (int64_t)nread
		);
	}
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	if(ret_nch)*ret_nch=nch;
	return buf;
}
int main(int argc, char **argv)
{
	const char *fn1=0, *fn2=0;

	if(argc!=3)
	{
		printf("Usage:  \"%s\"  image1.pnm  image2.pnm\n", argv[0]);
		printf("  Prints RMSE & PSNR\n");
		return 1;
	}
	fn1=argv[1];
	fn2=argv[2];

	int iw=0, ih=0, nch=0, iw2=0, ih2=0, nch2=0;
	uint8_t *image1=pnm_load(fn1, &iw, &ih, &nch);
	uint8_t *image2=pnm_load(fn2, &iw2, &ih2, &nch2);
	if(iw!=iw2||ih!=ih2||nch!=nch2)
	{
		CRASH(
			"Dimension mismatch  im1 %d*%d*%d vs %d*%d*%d"
			, nch, iw, ih
			, nch2, iw2, ih2
		);
	}
	{
		long long errors[4]={0};
		ptrdiff_t idx=-1;
		const uint8_t *ptr1=image1, *ptr2=image2;
		for(ptrdiff_t k=0, res=(ptrdiff_t)nch*iw*ih;k<res;k+=3)
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
		{
			double invres=1/((double)iw*ih);
			double rmse[]=
			{
				sqrt((double)errors[0]*invres),
				sqrt((double)errors[1]*invres),
				sqrt((double)errors[2]*invres),
				sqrt((double)errors[3]*invres/nch),
			};
			double psnr[]=
			{
				-20*log10(rmse[0]*(1./255)),
				-20*log10(rmse[1]*(1./255)),
				-20*log10(rmse[2]*(1./255)),
				-20*log10(rmse[3]*(1./255)),
			};
			printf("%12.6f\n\n", psnr[3]);
			if(nch==3)
				printf(
					"RMSE, PSNR\n"
					"T %12.6f/255  %12.6f dB\n"
					"R %12.6f/255  %12.6f dB\n"
					"G %12.6f/255  %12.6f dB\n"
					"B %12.6f/255  %12.6f dB\n"
					, rmse[3], psnr[3]
					, rmse[0], psnr[0]
					, rmse[1], psnr[1]
					, rmse[2], psnr[2]
				);
			else
				printf(
					"RMSE, PSNR\n"
					"G %12.6f/255  %12.6f dB\n"
					, rmse[3], psnr[3]
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
		}
	}
	free(image1);
	free(image2);
	return 0;
}
