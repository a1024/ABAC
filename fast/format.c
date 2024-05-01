#include"fast.h"
#include<string.h>
#include<ctype.h>
#include"lodepng.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
#ifdef __GNUC__
#define _stricmp strcasecmp
#endif
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
				printf("\n%s error IDX %td  XY (%5td, %5td) / %5d x %5d  b1 != b0\n", name, k, kx, ky, iw, ih);
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
	if(!_stricmp(fn+fnlen-4, ".PPM"))//shortcut to load PPMs fast
	{
		int iw, ih;
		ArrayHandle src=load_file(fn, 1, 16, 0);
		if(!src)
			return 0;
		char *ptr=(char*)src->data, *end=(char*)src->data+src->count;
		if(memcmp(ptr, "P6\n", sizeof(char[3])))
			goto load_other;
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
		ptrdiff_t res=(ptrdiff_t)iw*ih;
		image->data=(short*)malloc(res*sizeof(short[3]));
		if(!image->data)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		image->iw=iw;
		image->ih=ih;
		image->nch=3;
		image->depth=8;

		//serial 16-bit depth code
		res*=3;
		for(ptrdiff_t k=0;k<res;++k)
			image->data[k]=(short)((int)(unsigned char)ptr[k]-128);

		//this is from 32-bit depth code
#if 0
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
		short *dst2=(short*)dst;
		for(;k<res;++k)
		{
			dst2[0]=(unsigned char)(ptr[0]-128);
			dst2[1]=(unsigned char)(ptr[1]-128);
			dst2[2]=(unsigned char)(ptr[2]-128);
			dst2[3]=0;
			dst2+=4;
			ptr+=3;
		}
#endif
		array_free(&src);
		return 1;
	}
load_other:
	if(!_stricmp(fn+fnlen-4, ".DNG"))
		return load_dng(fn, image);
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
		dst[k]=(unsigned char)((src->data[k]+half)>>(src->depth-8)&mask);
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