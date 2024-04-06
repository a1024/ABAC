#include"e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"lodepng.h"
//#define J40_IMPLEMENTATION
//#define J40_CONFIRM_THAT_THIS_IS_EXPERIMENTAL_AND_POTENTIALLY_UNSAFE
//#include"j40.h"
#ifndef _MSC_VER
#include<byteswap.h>
#define _byteswap_ushort bswap_16
#endif
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

#if 0
unsigned char* image_load(const char *filename, int *iw, int *ih)
{
	unsigned char *buf=0;
#ifdef J40_IMPLEMENTATION
	int len=(int)strlen(filename);
	if(len<3)
		return 0;
	if(!acme_stricmp(filename+len-4, ".jxl"))
	{
		j40_image image;
		j40_from_file(&image, filename);
		j40_output_format(&image, J40_RGBA, J40_U8X4);
		
		// JPEG XL supports animation, so `j40_next_frame` calls can be called multiple times
		if(!j40_next_frame(&image))
			LOG_ERROR("Error: %s\n", j40_error_string(&image));
		else
		{
			j40_frame frame = j40_current_frame(&image);
			j40_pixels_u8x4 pixels = j40_frame_pixels_u8x4(&frame, J40_RGBA);
			
			if(j40_error(&image))
			{
				LOG_ERROR("Error: %s\n", j40_error_string(&image));
				return 0;
			}
			*iw=pixels.width;
			*ih=pixels.height;
			buf=(unsigned char*)malloc(pixels.width*pixels.height*4);
			if(!buf)
			{
				LOG_ERROR("Allocation error");
				return 0;
			}
			for(int ky=0;ky<pixels.height;++ky)
			{
				unsigned char *row=(unsigned char*)j40_row_u8x4(pixels, ky);
				memcpy(buf+(pixels.width*ky<<2), row, pixels.width*4);
			}
#if 0
			fprintf(out,
				"P7\n"
				"WIDTH %d\n"
				"HEIGHT %d\n"
				"DEPTH 4\n"
				"MAXVAL 255\n"
				"TUPLTYPE RGB_ALPHA\n"
				"ENDHDR\n",
			pixels.width, pixels.height);
			for(int y=0;y<pixels.height;++y)
				fwrite(j40_row_u8x4(pixels, y), 4, pixels.width, out);
#endif
		}
		j40_free(&image); // also frees all memory associated to j40_frame etc.
	}
	else
#endif
		buf=stbi_load(filename, iw, ih, 0, 4);

	return buf;
}
#endif
int image_save_png_rgba8(const char *filename, const unsigned char *image, int iw, int ih)
{
	return lodepng_encode_file(filename, image, iw, ih, LCT_RGBA, 8);
}
void save_32bit(const char *filename, const int *buf, int iw, int ih, int nch, int saveas8bit)
{
	LodePNGColorType type;
	switch(nch)
	{
	case 1:type=LCT_GREY;break;
	case 3:type=LCT_RGB;break;
	case 4:type=LCT_RGBA;break;
	default:LOG_ERROR("Bitmap type error");return;
	}
	int n=iw*ih*nch;
	if(saveas8bit)
	{
		unsigned char *b3=(unsigned char*)malloc(n);
		if(!b3)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned val=buf[k];
			b3[k]=(unsigned char)(val>>24);
		}
		lodepng_encode_file(filename, b3, iw, ih, type, 8);
		free(b3);
	}
	else
	{
		short *b2=(short*)malloc(n*sizeof(short));
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]>>16;
			b2[k]=_byteswap_ushort(val);
		}
		lodepng_encode_file(filename, (unsigned char*)b2, iw, ih, type, 16);
		free(b2);
	}
}
void save_16bit(const char *filename, const short *buf, const short *sub_b2, int iw, int ih, int nch, int val_offset, int val_shift, int saveas8bit)
{
	LodePNGColorType type;
	switch(nch)
	{
	case 1:type=LCT_GREY;break;
	case 3:type=LCT_RGB;break;
	case 4:type=LCT_RGBA;break;
	default:LOG_ERROR("Bitmap type error");return;
	}
	int n=iw*ih*nch;
	if(saveas8bit)
	{
		unsigned char *b3=(unsigned char*)malloc(n);
		if(!b3)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]+val_offset;
			if(sub_b2)
				val-=sub_b2[k];
			val<<=val_shift;
			b3[k]=(unsigned char)(val>>8);
		}
		lodepng_encode_file(filename, b3, iw, ih, type, 8);
		free(b3);
	}
	else
	{
		short *b2=(short*)malloc(n*sizeof(short));
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]+val_offset;
			if(sub_b2)
				val-=sub_b2[k];
			val<<=val_shift;
			b2[k]=_byteswap_ushort(val);
		}
		lodepng_encode_file(filename, (unsigned char*)b2, iw, ih, type, 16);
		free(b2);
	}
}
void save_mono8(const char *filename, unsigned char *buf, int iw, int ih, int stride)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	unsigned char *b2=(unsigned char*)malloc(res);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(ptrdiff_t k=0;k<res;++k)
		b2[k]=buf[k*stride];
	lodepng_encode_file(filename, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}
void save_channel(unsigned char *buf, int iw, int ih, int stride, int val_offset, const char *format, ...)
{
	va_list args;
	va_start(args, format);
	vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);

	unsigned char *b2=(unsigned char*)malloc((size_t)iw*ih);
	if(!b2)
		return;
	for(int ks=0, kd=0, res=iw*ih;kd<res;ks+=stride, ++kd)
		b2[kd]=buf[ks]+val_offset;
	
	lodepng_encode_file(g_buf, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}

Image* image_alloc(int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth, int initialize, int initval)
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
	image->depth[0]=rdepth;
	image->depth[1]=gdepth;
	image->depth[2]=bdepth;
	image->depth[3]=adepth;
	if(initialize)
		memfill(image->data, &initval, res*sizeof(int[4]), sizeof(int));
	return image;
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
		int val;
		val=(image->data[k|0]>>shift[0])+128, dst[0][k|0]=CLAMP(0, val, 255);
		val=(image->data[k|1]>>shift[1])+128, dst[0][k|1]=CLAMP(0, val, 255);
		val=(image->data[k|2]>>shift[2])+128, dst[0][k|2]=CLAMP(0, val, 255);
		val=(image->data[k|3]>>shift[3])+128, dst[0][k|3]=override_alpha?0xFF:CLAMP(0, val, 255);
		//dst[0][k|0]=(image->data[k|0]>>shift[0])+128;
		//dst[0][k|1]=(image->data[k|1]>>shift[1])+128;
		//dst[0][k|2]=(image->data[k|2]>>shift[2])+128;
		//dst[0][k|3]=override_alpha?0xFF:image->data[k|2]>>shift[3];
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
int image_snapshot(Image const *image)
{
	acme_strftime(g_buf, G_BUF_SIZE, "%Y%m%d_%H-%M-%S.PNG");
	return image_save_uint8(g_buf, image, 1);
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
void calc_depthfromdata(int *image, int iw, int ih, char *depths, const char *src_depths)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	for(int kc=0;kc<3;++kc)
	{
		int vmin=image[0<<2|kc], vmax=image[0<<2|kc];
		for(ptrdiff_t k=1;k<res;++k)
		{
			int sym=image[k<<2|kc];
			if(vmin>sym)
				vmin=sym;
			if(vmax<sym)
				vmax=sym;
		}
		int nlevels=vmax-vmin;
		depths[kc]=nlevels?floor_log2(nlevels)+1:0;
		depths[kc]=MAXVAR(depths[kc], src_depths[kc]);
	}
}
void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8)
{
	int nlevels=1<<depth;
	int shift=MAXVAR(0, depth-8);
	memset(hist, 0, nlevels*sizeof(int));
	if(hist8)
		memset(hist8, 0, 256*sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int sym=buf[(iw*ky+kx)<<2|kc]+(nlevels>>1);
			sym=CLAMP(0, sym, nlevels-1);
			++hist[sym];
			if(hist8)
				++hist8[sym>>shift];
		}
	}
	(void)ih;
}
double calc_entropy(const int *hist, int nlevels, int sum)//nlevels = 1<<current_depth
{
	double entropy=0;
	if(sum<0)
	{
		sum=0;
		for(int k=0;k<nlevels;++k)
			sum+=hist[k];
	}
	if(!sum)
		return 0;//degenerate distribution
	for(int k=0;k<nlevels;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/sum;
			entropy-=p*log2(p);		//Shannon's law
		}
	}
	return entropy;//invCR = entropy/original_depth, csize=res*entropy/8
}
void calc_csize(Image const *image, double *csizes)
{
	int res=image->iw*image->ih;
	int maxdepth=calc_maxdepth(image, 0), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
		int depth=image->depth[kc];
		calc_histogram(image->data, image->iw, image->ih, kc, 0, image->iw, 0, image->ih, depth, hist, 0);
		double entropy=calc_entropy(hist, 1<<depth, res);
		csizes[kc]=res*entropy/8;
		//double invCR=entropy/image->src_depth[kc];
		//csizes[kc]=res*image->src_depth[kc]*entropy/8;
	}
	free(hist);
	//int hist[256]={0};
	//for(int k=0;k<res;++k)
	//{
	//	unsigned char val=src[k<<2|kc]+128;
	//	++hist[val];
	//}
	//double csize=calc_csize_from_hist(hist, 256, 0);
	//return csize;
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