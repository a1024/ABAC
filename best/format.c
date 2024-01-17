#include"best.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"lodepng.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

int calc_maxdepth(Image const *image, int *inflation)
{
	int maxdepth=image->depth[0]+(inflation?inflation[0]:0);
	int next=image->depth[1]+(inflation?inflation[1]:0);
	if(maxdepth<next)
		maxdepth=next;
	next=image->depth[2]+(inflation?inflation[2]:0);
	if(maxdepth<next)
		maxdepth=next;
	return maxdepth;
}
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
void image_export_uint8(Image const *image, unsigned char **dst, int override_alpha)
{
	if(!image)
		return;
	void *p=realloc(*dst, image->iw*image->ih*sizeof(char[4]));
	if(!p)
		return;
	*dst=(unsigned char*)p;
	int shift[]=
	{
		MAXVAR(8, image->depth[0])-8,
		MAXVAR(8, image->depth[1])-8,
		MAXVAR(8, image->depth[2])-8,
		MAXVAR(8, image->depth[3])-8,
	};
	for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih*4;k<res;k+=4)
	{
		dst[0][k|0]=(image->data[k|0]>>shift[0])+128;
		dst[0][k|1]=(image->data[k|1]>>shift[1])+128;
		dst[0][k|2]=(image->data[k|2]>>shift[2])+128;
		dst[0][k|3]=override_alpha?0xFF:image->data[k|2]>>shift[3];
	}
}
void image_export_uint16(Image const *image, unsigned short **dst, int override_alpha, int big_endian)
{
	if(!image)
		return;
	void *p=realloc(*dst, (size_t)image->iw*image->ih*sizeof(short[4]));
	if(!p)
		return;
	*dst=(unsigned short*)p;
	int shift[]=
	{
		16-MAXVAR(16, image->depth[0]),
		16-MAXVAR(16, image->depth[1]),
		16-MAXVAR(16, image->depth[2]),
		16-MAXVAR(16, image->depth[3]),
	};
	for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih*4;k<res;k+=4)
	{
		int val;
		val=SHIFT_LEFT_SIGNED(image->data[k|0], shift[0])+0x8000, dst[0][k|0]=CLAMP(-0x8000, val, 0x7FFF);
		val=SHIFT_LEFT_SIGNED(image->data[k|1], shift[1])+0x8000, dst[0][k|1]=CLAMP(-0x8000, val, 0x7FFF);
		val=SHIFT_LEFT_SIGNED(image->data[k|2], shift[2])+0x8000, dst[0][k|2]=CLAMP(-0x8000, val, 0x7FFF);
		if(override_alpha)
			dst[0][k|3]=0xFFFF;
		else
			val=SHIFT_LEFT_SIGNED(image->data[k|3], shift[3])+0x8000, dst[0][k|3]=CLAMP(-0x8000, val, 0x7FFF);
	}
	if(big_endian)
	{
		for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih*4;k<res;++k)
		{
			unsigned short val=dst[0][k];
			val=val>>8|val<<8;
			dst[0][k]=val;
		}
	}
}
int image_save_uint8(const char *fn, Image const *image, int override_alpha)
{
	unsigned char *dst=0;
	image_export_uint8(image, &dst, override_alpha);
	if(!dst)
		return 0;
	int error=lodepng_encode_file(fn, dst, image->iw, image->ih, LCT_RGBA, 8);
	free(dst);
	return !error;
}
int image_save_native(const char *fn, Image const *image, int override_alpha)
{
	int maxdepth=calc_maxdepth(image, 0);
	if(maxdepth<=8)
		return image_save_uint8(fn, image, override_alpha);
	unsigned short *dst=0;
	image_export_uint16(image, &dst, override_alpha, 1);
	int error=lodepng_encode_file(fn, (unsigned char*)dst, image->iw, image->ih, LCT_RGBA, 16);
	free(dst);
	return !error;
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

#define LSIM_TAG "LSvA"
size_t lsim_writeheader(ArrayHandle *dst, int iw, int ih, int nch, const char *depths, int codec_id)//returns number of bytes written
{
	size_t len=snprintf(g_buf, G_BUF_SIZE, LSIM_TAG " %d %d %d", iw, ih, nch);
	for(int k=0;k<nch;++k)
		len+=snprintf(g_buf+len, G_BUF_SIZE-len, " %d", depths[k]);
	len+=snprintf(g_buf+len, G_BUF_SIZE-len, " c%d\n", codec_id);
	array_append(dst, g_buf, 1, len, 1, 1, 0);
	return len;
}
size_t lsim_readheader(const unsigned char *src, size_t len, LSIMHeader *dst)//returns number of bytes read
{
	const unsigned char *ptr=src, *end=src+len;
#define SKIP_SPACE() for(;ptr<end&&isspace(*ptr);++ptr)
#define GET_INT() strtol((char*)ptr, (char**)&ptr, 10)
	memset(dst, 0, sizeof(*dst));
	if(ptr+4>end||memcmp(ptr, LSIM_TAG, 4))
		return 0;
	ptr+=4;
	SKIP_SPACE();

	dst->iw=GET_INT();
	SKIP_SPACE();

	dst->ih=GET_INT();
	SKIP_SPACE();

	dst->nch=GET_INT();
	SKIP_SPACE();
	
	for(int kc=0;kc<dst->nch;++kc)
	{
		dst->depth[kc]=(char)GET_INT();
		SKIP_SPACE();
	}

	if(ptr<end&&*ptr=='c')
	{
		++ptr;
		dst->codec_id=GET_INT();
	}
	ptr+=ptr<end&&*ptr=='\n';
#undef  SKIP_SPACE
#undef  GET_INT
	return ptr-src;
}
void image_from_lsimheader(Image **dst, LSIMHeader const *src)
{
	size_t size=sizeof(Image)+src->iw*src->ih*sizeof(int[4]);
	void *p=realloc(*dst, size);
	if(!p)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	*dst=(Image*)p;
	memset(*dst, 0, size);
	dst[0]->iw=src->iw;
	dst[0]->ih=src->ih;
	dst[0]->nch=src->nch;
	for(int kc=0;kc<src->nch;++kc)
		dst[0]->src_depth[kc]=dst[0]->depth[kc]=src->depth[kc];
}