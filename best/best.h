#pragma once
#ifndef INC_BEST_H
#define INC_BEST_H
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif

//transforms
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);
void colortransform_ycocb_fwd(char *buf, int iw, int ih);
void colortransform_ycocb_inv(char *buf, int iw, int ih);
	
int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);
	
#ifdef __cplusplus
}
#endif
#endif/INC_BEST_H
