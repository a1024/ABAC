#pragma once
#ifndef INC_IMGCVT_H
#define INC_IMGCVT_H
#include"util.h"


#if 0
typedef struct ImageStruct
{
	int iw, ih, nch, depth;
	short *data;
} Image;
int compare_bufs_16(const short *b1, const short *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);
int image_load(const char *fn, Image *image);
int image_save8(const char *fn, Image const *image);
int image_copy(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_copy_nodata(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_snapshot8(Image const *image);
void image_clear(Image *image);
size_t image_getBMPsize(Image const *image);
#endif


#endif