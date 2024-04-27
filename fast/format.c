#include"fast.h"
#include<string.h>
#include"lodepng.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

int compare_bufs_16(const short *b1, const short *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud)
{
	ptrdiff_t len=(ptrdiff_t)chstride*iw*ih;
	int inc=chstride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-chstride:0;k>=0&&k<len;k+=inc)
	{
		if(memcmp(b1+k, b0+k, nch*sizeof(short)))
		{
			if(loud)
			{
				ptrdiff_t idx=k/chstride, kx=idx%iw, ky=idx/iw;
				printf("\n%s error IDX %lld  XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, k, kx, ky, iw, ih);
				for(int kc=0;kc<nch;++kc)
				{
					char c=(unsigned short)b1[k+kc]==(unsigned short)b0[k+kc]?'=':'!';
					printf("C%d  0x%04X %c= 0x%04X    %d %c= %d\n",
						kc,
						(unsigned short)b1[k+kc], c, (unsigned short)b0[k+kc],
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
int image_load(const char *fn, Image *image)
{
	int fnlen=(int)strlen(fn);
	if(!_stricmp(fn+fnlen-4, ".DNG"))
	{
		load_dng(fn, image);
		return 1;
	}
	int idx_start=0, idx_end=0;
	get_filetitle(fn, 0, &idx_start, &idx_end);
	image->data=(short*)stbi_load_16(fn, &image->iw, &image->ih, &image->nch, 0);
	if(!image->data)
		return 0;
	ptrdiff_t nvals=(ptrdiff_t)image->iw*image->ih*image->nch;
	image->depth=8;
	for(ptrdiff_t k=0;k<nvals;++k)
	{
		unsigned short val=image->data[k];
		if((val&0xFF)!=(val>>8&0xFF))
		{
			image->depth=16;
			break;
		}
	}
	if(image->depth==8)
	{
		for(ptrdiff_t k=0;k<nvals;++k)
		{
			short val=image->data[k];
			val&=0xFF;
			val-=128;
			image->data[k]=val;
		}
	}
	return 1;
}
unsigned char* image_export8(Image const *src)
{
	ptrdiff_t nvals=(ptrdiff_t)src->iw*src->ih*src->nch;
	unsigned char *dst=(unsigned char*)malloc(nvals*sizeof(char));
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	int nlevels=1<<src->depth, half=nlevels>>1, mask=nlevels-1;
	for(ptrdiff_t k=0;k<nvals;++k)
		dst[k]=(src->data[k]+half)>>(src->depth-8)&mask;
	return dst;
}
int image_save8(const char *fn, Image const *image)
{
	unsigned char *buf=image_export8(image);
	if(!buf)
		return 0;
	LodePNGColorType type=LCT_GREY;
	switch(image->nch)
	{
	case 1:type=LCT_GREY;break;
	case 2:type=LCT_GREY_ALPHA;break;
	case 3:type=LCT_RGB;break;
	case 4:type=LCT_RGBA;break;
	}
	if(!type)
	{
		LOG_ERROR("Invalid number of channels %d", image->nch);
		return 0;
	}
	lodepng_encode_file(fn, buf, image->iw, image->ih, type, 8);
	return 1;
}
int image_snapshot8(Image const *image)
{
	acme_strftime(g_buf, G_BUF_SIZE, "%Y%m%d_%H-%M-%S.PNG");
	return image_save8(g_buf, image);
}
int image_copy(Image *dst, Image const *src)//dst must be initialized to zero, or another image
{
	ptrdiff_t nvals=(ptrdiff_t)src->iw*src->ih*src->nch;
	void *p=realloc(dst->data, nvals*sizeof(short));
	if(!p)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memcpy(dst, src, sizeof(*dst));
	dst->data=(short*)p;
	memcpy(dst->data, src->data, nvals*sizeof(short));
	return 1;
}
int image_copy_nodata(Image *dst, Image const *src)//dst must be initialized to zero, or another image
{
	ptrdiff_t nvals=(ptrdiff_t)src->iw*src->ih*src->nch;
	void *p=realloc(dst->data, nvals*sizeof(short));
	if(!p)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memcpy(dst, src, sizeof(*dst));
	dst->data=(short*)p;
	//memcpy(dst->data, src->data, nvals*sizeof(short));
	return 1;
}
void image_clear(Image *image)
{
	free(image->data);
	memset(image, 0, sizeof(*image));
}
size_t image_getBMPsize(Image const *image)
{
	return ((size_t)image->iw*image->ih*image->nch*image->depth+7)/8;
}