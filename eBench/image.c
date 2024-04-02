#include"ebench.h"
#include"lodepng.h"
#include<stdlib.h>
#include<string.h>
static const char file[]=__FILE__;

Image* image_from_uint8(const unsigned char *src, int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	Image *image=(Image*)malloc(sizeof(Image)+res*sizeof(int[4]));
	if(!image)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	image->iw=iw;
	image->ih=ih;
	image->nch=nch;
	image->src_depth[0]=image->depth[0]=rdepth;
	image->src_depth[1]=image->depth[1]=gdepth;
	image->src_depth[2]=image->depth[2]=bdepth;
	image->src_depth[3]=image->depth[3]=adepth;
	if(src)
	{
		int offset[]=
		{
			image->depth[0]?1<<(image->depth[0]-1):0,
			image->depth[1]?1<<(image->depth[1]-1):0,
			image->depth[2]?1<<(image->depth[2]-1):0,
			image->depth[3]?1<<(image->depth[3]-1):0,
		};
		res<<=2;
		for(ptrdiff_t k=0;k<res;++k)
			image->data[k]=src[k]-offset[k&3];
	}
	else
		memset(image->data, 0, res*sizeof(int[4]));
	return image;
}
Image* image_from_uint16(const unsigned short *src, int iw, int ih, int nch, char *src_depths, char *dst_depths)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	Image *image=(Image*)malloc(sizeof(Image)+res*sizeof(int[4]));
	if(!image)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	image->iw=iw;
	image->ih=ih;
	image->nch=nch;
	image->src_depth[0]=image->depth[0]=dst_depths[0];
	image->src_depth[1]=image->depth[1]=dst_depths[1];
	image->src_depth[2]=image->depth[2]=dst_depths[2];
	image->src_depth[3]=image->depth[3]=dst_depths[3];
	if(src)
	{
		int offset[]=
		{
			image->depth[0]?1<<(image->depth[0]-1):0,
			image->depth[1]?1<<(image->depth[1]-1):0,
			image->depth[2]?1<<(image->depth[2]-1):0,
			image->depth[3]?1<<(image->depth[3]-1):0,
		};
		int shift[]=
		{
			dst_depths[0]-src_depths[0],
			dst_depths[1]-src_depths[1],
			dst_depths[2]-src_depths[2],
			dst_depths[3]-src_depths[3],//this sets alpha to zero when alpha depth is zero
		};
		res<<=2;
		for(ptrdiff_t k=0;k<res;++k)
			image->data[k]=(shift[k&3]<0?src[k]>>-shift[k&3]:src[k]<<shift[k&3])-offset[k&3];
	}
	else
		memset(image->data, 0, res*sizeof(int[4]));
	return image;
}
Image* image_load(const char *fn)
{
	int iw=0, ih=0, nch=0;
	unsigned short *src=stbi_load_16(fn, &iw, &ih, &nch, 4);
	if(!src)
	{
		LOG_ERROR("Cannot open %s", fn);
		return 0;
	}
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	char src_depth[]={16, 16, 16, 16};
	char dst_depth[]={8, 8, 8, 0};
	for(ptrdiff_t k=0;k<res;++k)//detect 8-bit
	{
		int r=src[k<<2|0], g=src[k<<2|1], b=src[k<<2|2], a=src[k<<2|3];
		if(dst_depth[0]<16&&(r>>8)!=(r&0xFF))dst_depth[0]=16;
		if(dst_depth[1]<16&&(g>>8)!=(g&0xFF))dst_depth[1]=16;
		if(dst_depth[2]<16&&(b>>8)!=(b&0xFF))dst_depth[2]=16;
		if(a!=0xFFFF)
		{
			if(dst_depth[3]<16&&(a>>8)!=(a&0xFF))
				dst_depth[3]=16;
			else if(!dst_depth[3])
				dst_depth[3]=8;
		}
	}
	Image *image=image_from_uint16(src, iw, ih, nch, src_depth, dst_depth);
	free(src);
	return image;
}
int image_save_uint8(const char *fn, Image const *image, int override_alpha)
{
	unsigned char *dst=0;
	image_export_uint8(image, &dst, override_alpha, 0);
	if(!dst)
		return 0;
	lodepng_encode_file(fn, dst, image->iw, image->ih, LCT_RGBA, 8);
	free(dst);
	return 1;
}
void image_export_uint8(Image const *image, unsigned char **dst, int override_alpha, int swap_rb)
{
	if(!image)
		return;
	void *p=realloc(*dst, image->iw*image->ih*sizeof(char[4]));
	if(!p)
		return;
	*dst=(unsigned char*)p;
	int shift[]=
	{
		MAXVAR(0, image->depth[0]-8),
		MAXVAR(0, image->depth[1]-8),
		MAXVAR(0, image->depth[2]-8),
		MAXVAR(0, image->depth[3]-8),
	};
	int r, b;
	if(swap_rb)
		r=2, b=0;
	else
		r=0, b=2;
	for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih*4;k<res;k+=4)
	{
		dst[0][k|r]=(image->data[k|0]>>shift[0])+128;
		dst[0][k|1]=(image->data[k|1]>>shift[1])+128;
		dst[0][k|b]=(image->data[k|2]>>shift[2])+128;
		dst[0][k|3]=override_alpha?0xFF:image->data[k|2]>>shift[3];
	}
}
double image_getBMPsize(Image const *image)
{
	int totalbits=0;
	for(int k=0;k<image->nch;++k)
		totalbits+=image->src_depth[k];
	return (double)image->iw*image->ih*totalbits/8;
}
size_t image_getbufsize(Image const *image)
{
	if(!image)
		return 0;
	return sizeof(Image)+(size_t)image->iw*image->ih*sizeof(int[4]);
}
void image_copy_nodata(Image **dst, Image const *src)
{
	size_t srcsize=image_getbufsize(src);
	if(!srcsize)
		return;
	void *p=realloc(*dst, srcsize);
	if(!p)
		return;
	*dst=(Image*)p;
	memcpy(*dst, src, sizeof(Image));
}
void image_copy(Image **dst, Image const *src)
{
	size_t srcsize=image_getbufsize(src);
	if(!srcsize)
		return;
	void *p=realloc(*dst, srcsize);
	if(!p)
		return;
	*dst=(Image*)p;
	memcpy(*dst, src, srcsize);
}