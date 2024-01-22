#pragma once
#ifndef INC_BEST_H
#define INC_BEST_H
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif


typedef struct ImageStruct
{
	int iw, ih,
		nch;//{greyscale, greyscale+alpha, RGB, RGB+alpha}	alpha can be ignored for now
	char depth[4];
	char src_depth[4];//for entropy calculations
	int data[];//stride always sizeof(int[4])
} Image;
Image* image_load(const char *fn);
int image_save_uint8(const char *fn, Image const *image, int override_alpha);
int image_save_native(const char *fn, Image const *image, int override_alpha);
Image* image_from_uint8(const unsigned char *src, int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth);
Image* image_from_uint16(const unsigned short *src, int iw, int ih, int nch, char *src_depths, char *dst_depths);
void image_export_uint8(Image const *image, unsigned char **dst, int override_alpha);
void image_export_uint16(Image const *image, unsigned short **dst, int override_alpha, int big_endian);
double image_getBMPsize(Image const *image);
size_t image_getbufsize(Image const *image);
void image_copy_nodata(Image **dst, Image const *src);//dst must be 0 or a valid pointed
void image_copy(Image **dst, Image const *src);//dst must be 0 or a valid pointed


typedef struct LSIMHeaderStruct
{
	int iw, ih, nch;
	char depth[4];
	int codec_id;
} LSIMHeader;
size_t lsim_writeheader(ArrayHandle *dst, int iw, int ih, int nch, const char *depths, int codec_id);//returns number of bytes written
size_t lsim_readheader(const unsigned char *src, size_t srclen, LSIMHeader *dst);//returns number of bytes read
void image_from_lsimheader(Image **dst, LSIMHeader const *src);//dst must be 0 or a valid pointed

//auxiliary functions
int calc_maxdepth(Image const *image, int *inflation);
int compare_bufs_32(const int *b1, const int *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);
int get_nch(const char *buf, int res);//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}
int get_nch32(const int *buf, int res);//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}

//transforms
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);
void rct_JPEG2000_32(Image *image, int fwd);
void colortransform_JPEG2000_fwd(char *buf, int iw, int ih);
void colortransform_JPEG2000_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v0_fwd(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v0_inv(char *buf, int iw, int ih);
	
int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

//CALIC
int t45_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t45_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *pixels, int loud);

//SLICv4
int t46_encode(Image const *src, ArrayHandle *data, int loud);
int t46_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//SLICv5
int t47_encode(Image const *src, ArrayHandle *data, int loud);
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);


#ifdef __cplusplus
}
#endif
#endif//INC_BEST_H
