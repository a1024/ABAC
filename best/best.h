#pragma once
#ifndef INC_BEST_H
#define INC_BEST_H
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif

//transforms
int get_nch(const char *buf, int res);//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);
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
	
#ifdef __cplusplus
}
#endif
#endif//INC_BEST_H
