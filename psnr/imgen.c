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
	int iw=0, ih=0, nch=0;
	int tag=0;
	int vmax=0;
	int nscan;
	ptrdiff_t size, nread;
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", fn);
		return 0;
	}
	fread(&tag, 1, 2, fsrc);
	if(tag!=('P'|'6'<<8)&&tag!=('P'|'5'<<8))
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	if(fgetc(fsrc)!='\n')
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	nch=tag==('P'|'6'<<8)?3:1;
	nscan=fscanf(fsrc, "%d %d\n", &iw, &ih);
	if(nscan!=2)
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	nscan=fscanf(fsrc, "%d", &vmax);
	if(nscan!=1||vmax!=255)
	{
		CRASH("Invalid PPM file \"%s\"", fn);
		return 0;
	}
	size=(ptrdiff_t)nch*iw*ih;
	uint8_t *buf=(uint8_t*)malloc(size+16);
	if(!buf)
	{
		CRASH("Alloc error");
		return 0;
	}
	memset(buf, 0, size+16);
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
int pnm_save(const char *fn, uint8_t *image, int iw, int ih, int nch)
{
	FILE *fdst=0;

	fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return 1;
	}
	fprintf(fdst, "P%d\n%d %d\n255\n", nch==1?5:6, iw, ih);
	fwrite(image, 1, (size_t)nch*iw*ih, fdst);
	fclose(fdst);
	return 0;
}
int main(int argc, char **argv)
{
	const char *srcfn=0, *dstfn=0;
	int c1=0, c2=0, sign=0, op=0;
	int iw=0, ih=0, nch=0;
	uint8_t *image=0, *imptr=0, *imend=0;
	int64_t size=0;
	int offset=0;

	if(argc!=7)
	{
		printf("Usage:  \"%s\"  SRC.pnm  coeff1  coeff2  s|u  c|m  DST.pnm\n", argv[0]);
		printf("What this does:\n");
		printf("  s:  noise = (uint8_t)rand - 128\n");
		printf("  u:  noise = (uint8_t)rand\n");
		printf("  value = (coeff1*SRC + coeff2*noise)>>8\n");
		printf("  c:  DST  =  Clamp(value,  0,  255)\n");
		printf("  m:  DST  =  value  Mod 256\n");
		return 1;
	}
	srcfn=argv[1];
	c1=atoi(argv[2]);
	c2=atoi(argv[3]);

	sign=argv[4][0];
	if(!argv[4][1]&&(sign&0xDF)=='S')
		sign=1;
	else if(!argv[4][1]&&(sign&0xDF)=='U')
		sign=0;
	else CRASH("Invalid operation \"%s\"", argv[4]);

	op=argv[5][0];
	if(!argv[5][1]&&(op&0xDF)=='C')
		op=1;
	else if(!argv[5][1]&&(op&0xDF)=='M')
		op=0;
	else CRASH("Invalid operation \"%s\"", argv[5]);

	dstfn=argv[6];

	image=pnm_load(srcfn, &iw, &ih, &nch);
	size=(int64_t)nch*iw*ih;
	imptr=image;
	imend=image+size;
	offset=sign?128:0;
	while(imptr<imend)
	{
		int val=(c1 * *imptr + c2*((uint8_t)(rand()>>2)-offset))>>8;
		if(op)
		{
			if(val<0)val=0;
			if(val>255)val=255;
		}
		*imptr=(uint8_t)val;
		++imptr;
	}
	pnm_save(dstfn, image, iw, ih, nch);
	free(image);
	return 0;
}
