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
static ptrdiff_t ppm_save(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return 0;
	}
	ptrdiff_t size=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	ptrdiff_t usize=(ptrdiff_t)3*iw*ih;
	ptrdiff_t nwritten=fwrite(image, 1, usize, fdst);
	if(nwritten!=usize)
	{
		LOG_ERROR("Error saving \"%s\"", fn);
		return 0;
	}
	size+=nwritten;
	return size;
}
int main(int argc, char **argv)
{
	const char *fn1=0, *fn2=0, *dstfn=0;

	if(argc!=4)
	{
		printf("Usage:  \"%s\"  psrc.ppm  nsrc.ppm  dst.ppm\n", argv[0]);
		printf("  dst = psrc - nsrc\n");
		return 1;
	}
	fn1=argv[1];
	fn2=argv[2];
	dstfn=argv[3];

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

	unsigned char *ptr1=image1, *ptr2=image2;
	for(ptrdiff_t k=0, res=(ptrdiff_t)3*iw*ih;k<res;k+=3)
	{
		ptr1[0]=(unsigned char)(ptr1[0]-ptr2[0]+128);
		ptr1[1]=(unsigned char)(ptr1[1]-ptr2[1]+128);
		ptr1[2]=(unsigned char)(ptr1[2]-ptr2[2]+128);
		ptr1+=3;
		ptr2+=3;
	}

	ppm_save(dstfn, image1, iw, ih);
	free(image1);
	free(image2);
}
