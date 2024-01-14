#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"lodepng.h"
#define J40_IMPLEMENTATION
#define J40_CONFIRM_THAT_THIS_IS_EXPERIMENTAL_AND_POTENTIALLY_UNSAFE
#include"j40.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
#include"battle.h"
static const char file[]=__FILE__;

unsigned char* image_load(const char *filename, int *iw, int *ih)
{
	unsigned char *buf=0;

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
		buf=stbi_load(filename, iw, ih, 0, 4);

	return buf;
}
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