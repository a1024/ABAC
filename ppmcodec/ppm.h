#ifndef INC_PPM_H
#define INC_PPM_H
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef enum _CodecID
{
	CODEC_INVALID,
	CODEC_PPM,
	CODEC_PGM,
	CODEC_C01,
} CodecID;

static void ppm_skip(const unsigned char **ptr, const unsigned char *end)
{
	for(;*ptr<end;)
	{
		if(**ptr=='#')
		{
			for(;*ptr<end&&**ptr!='\n';++*ptr);
		}
		else if(**ptr<'!')
			++*ptr;
		else
			return;
	}
}
static int header_read(const unsigned char *src, int len, int *iw, int *ih, CodecID *codec)
{
	const unsigned char *ptr, *end;
	int vmax;

	ptr=src;
	end=src+len;
	if(!memcmp(ptr, "P6", 2))
	{
		*codec=CODEC_PPM;
		ptr+=2;
	}
	else if(!memcmp(ptr, "P5", 2))
	{
		*codec=CODEC_PGM;
		ptr+=2;
	}
	else if(!memcmp(ptr, "C01", 3))
	{
		*codec=CODEC_C01;
		ptr+=3;
	}
	else
		*codec=CODEC_INVALID;
	ppm_skip(&ptr, end);
	*iw=strtol((char*)ptr, (char**)&ptr, 10);
	ppm_skip(&ptr, end);
	*ih=strtol((char*)ptr, (char**)&ptr, 10);
	if(*codec==CODEC_PPM||*codec==CODEC_PGM)
	{
		ppm_skip(&ptr, end);
		vmax=strtol((char*)ptr, (char**)&ptr, 10);//255
		if(vmax!=255)
			*codec=CODEC_INVALID;
	}
	ptr+=*ptr=='\n';//ppm_skip skips binary data < 33 as well
	return (int)(ptr-src);
}
static int compare_bufs_8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud)
{
	ptrdiff_t len=(ptrdiff_t)chstride*iw*ih;
	int inc=chstride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-chstride:0;k>=0&&k<len;k+=inc)
	{
		if(memcmp(b1+k, b0+k, nch*sizeof(char)))
		{
			if(loud)
			{
				ptrdiff_t idx=k/chstride, kx=idx%iw, ky=idx/iw;
				printf("\n%s error IDX %td  XY (%5td, %5td) / %5d x %5d  b1 != b0\n", name, k, kx, ky, iw, ih);
				for(int kc=0;kc<nch;++kc)
				{
					char c=(unsigned char)b1[k+kc]==(unsigned char)b0[k+kc]?'=':'!';
					printf("C%d  0x%04X %c= 0x%04X    %d %c= %d\n",
						kc,
						(unsigned char)b1[k+kc], c, (unsigned char)b0[k+kc],
						b1[k+kc], c, b0[k+kc]
					);
				}
				LOG_ERROR("");
			}
			return 1;
		}
	}
	if(loud)
		printf("%s:\tSUCCESS\n", name);
	return 0;
}

#endif//INC_PPM_H
