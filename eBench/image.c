#include"ebench.h"
#include"lodepng.h"
#include<stdlib.h>
#include<string.h>
#include<tmmintrin.h>
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
Image* image_load(const char *fn, int fnlen)
{
	Image *image;
	int iw=0, ih=0, nch=0;
	ptrdiff_t res;
	if(!_stricmp(fn+fnlen-4, ".PPM"))//shortcut to load PPMs fast
	{
		ArrayHandle src;
		char *ptr, *end;

		src=load_file(fn, 1, 16, 0);
		if(!src)
			return 0;
		ptr=(char*)src->data;
		end=(char*)src->data+src->count;
		if(memcmp(ptr, "P6\n", sizeof(char[3])))
			goto load_other;
		nch=3;
		ptr+=3;
		if(*ptr=='#')//skip comments
		{
			++*ptr;
			for(;ptr!=end&&*ptr!='\n';++ptr);
			++*ptr;//skip newline
		}
		iw=(int)strtol(ptr, &ptr, 10);
		for(;ptr!=end&&isspace(*ptr);++ptr);
		ih=(int)strtol(ptr, &ptr, 10);
		++ptr;//skip newline
		for(;ptr!=end&&*ptr!='\n';++ptr);
		++ptr;//skip newline

		//start of binary data
		res=(ptrdiff_t)iw*ih;
		image=(Image*)malloc(sizeof(Image)+res*sizeof(int[4]));
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image, 0, sizeof(Image));
		image->iw=iw;
		image->ih=ih;
		image->nch=3;
		image->depth[0]=8;
		image->depth[1]=8;
		image->depth[2]=8;
		image->depth[3]=0;
		image->src_depth[0]=8;
		image->src_depth[1]=8;
		image->src_depth[2]=8;
		image->src_depth[3]=0;

		{
			__m128i extract0=_mm_set_epi8(
			//	15,|14, 13, 12,|11, 10,  9,| 8,  7,  6,| 5,  4,  3,| 2,  1,  0
				-1, -1, -1, -1, -1, -1, -1,  2, -1, -1, -1,  1, -1, -1, -1,  0
			);
			__m128i extract1=_mm_set_epi8(
			//	15,|14, 13, 12,|11, 10,  9,| 8,  7,  6,| 5,  4,  3,| 2,  1,  0
				-1, -1, -1, -1, -1, -1, -1,  5, -1, -1, -1,  4, -1, -1, -1,  3
			);
			__m128i extract2=_mm_set_epi8(
			//	15,|14, 13, 12,|11, 10,  9,| 8,  7,  6,| 5,  4,  3,| 2,  1,  0
				-1, -1, -1, -1, -1, -1, -1,  8, -1, -1, -1,  7, -1, -1, -1,  6
			);
			__m128i extract3=_mm_set_epi8(
			//	15,|14, 13, 12,|11, 10,  9,| 8,  7,  6,| 5,  4,  3,| 2,  1,  0
				-1, -1, -1, -1, -1, -1, -1, 11, -1, -1, -1, 10, -1, -1, -1,  9
			);
			__m128i extract4=_mm_set_epi8(
			//	15,|14, 13, 12,|11, 10,  9,| 8,  7,  6,| 5,  4,  3,| 2,  1,  0
				-1, -1, -1, -1, -1, -1, -1, 14, -1, -1, -1, 13, -1, -1, -1, 12
			);
			__m128i half=_mm_set1_epi32(128);
			__m128i *dst=(__m128i*)image->data;
			ptrdiff_t k=0;
			int *dst2;

			for(;k<res-4;k+=5)
			{
				__m128i packed5=_mm_loadu_si128((__m128i*)ptr);
				__m128i e0=_mm_shuffle_epi8(packed5, extract0);
				__m128i e1=_mm_shuffle_epi8(packed5, extract1);
				__m128i e2=_mm_shuffle_epi8(packed5, extract2);
				__m128i e3=_mm_shuffle_epi8(packed5, extract3);
				__m128i e4=_mm_shuffle_epi8(packed5, extract4);
				e0=_mm_sub_epi32(e0, half);
				e1=_mm_sub_epi32(e1, half);
				e2=_mm_sub_epi32(e2, half);
				e3=_mm_sub_epi32(e3, half);
				e4=_mm_sub_epi32(e4, half);
				_mm_storeu_si128(dst+0, e0);
				_mm_storeu_si128(dst+1, e1);
				_mm_storeu_si128(dst+2, e2);
				_mm_storeu_si128(dst+3, e3);
				_mm_storeu_si128(dst+4, e4);
				dst+=5;
				ptr+=3*5;
			}
			dst2=(int*)dst;
			for(;k<res;++k)
			{
				dst2[0]=(unsigned char)ptr[0]-128;
				dst2[1]=(unsigned char)ptr[1]-128;
				dst2[2]=(unsigned char)ptr[2]-128;
				dst2[3]=0;
				dst2+=4;
				ptr+=3;
			}
		}
		array_free(&src);
		return image;
	}
load_other:
	{
		unsigned short *src;
		char src_depth[]={16, 16, 16, 16};
		char dst_depth[]={8, 8, 8, 0};

		src=stbi_load_16(fn, &iw, &ih, &nch, 4);
		if(!src)
		{
			//LOG_ERROR("Cannot open %s", fn);
			return 0;
		}
		res=(ptrdiff_t)iw*ih;
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
		image=image_from_uint16(src, iw, ih, nch, src_depth, dst_depth);
		free(src);
		return image;
	}
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
	{
		void *ptr=realloc(*dst, sizeof(char[4])*image->iw*image->ih);
		if(!ptr)
			return;
		*dst=(unsigned char*)ptr;
	}
	{
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
			dst[0][k|r]=(unsigned char)((image->data[k|0]>>shift[0])+128);
			dst[0][k|1]=(unsigned char)((image->data[k|1]>>shift[1])+128);
			dst[0][k|b]=(unsigned char)((image->data[k|2]>>shift[2])+128);
			dst[0][k|3]=override_alpha?0xFF:(unsigned char)(image->data[k|3]>>shift[3])+128;
		}
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
	{
		void *ptr=realloc(*dst, srcsize);
		if(!ptr)
			return;
		*dst=(Image*)ptr;
	}
	memcpy(*dst, src, sizeof(Image));
}
void image_copy(Image **dst, Image const *src)
{
	size_t srcsize=image_getbufsize(src);
	if(!srcsize)
		return;
	{
		void *ptr=realloc(*dst, srcsize);
		if(!ptr)
			return;
		*dst=(Image*)ptr;
	}
	memcpy(*dst, src, srcsize);
}
