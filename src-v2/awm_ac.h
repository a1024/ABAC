#pragma once
#ifndef AWM_AC_H
#define AWM_AC_H
#include"awm_util.h"

#define		HAVE_AVX2

void			YCoCg_fwd(int *image, int count);
void			YCoCg_inv(int *image, int count);

//compression testbench

//typedef struct ContextStruct
//{
//	const int *src;//use only src for prediction
//	int *dst;
//	int iw, ih;
//	unsigned short prob[32];
//} Context, *ContextHandle;
//typedef void (*Predictor)(ContextHandle data, int x, int y);//should return 16-bit fixed-point probability, must use causal neighbors

typedef void (*Predictor)(const int *data, int iw, int ih, int bpp, int x, int y, unsigned short *hit_hist, unsigned short *ret_prob);//should return 16-bit fixed-point probability (*bpp), must use causal neighbors

int				test_encode(const int *src, int iw, int ih, int bpp, Predictor predictor, ArrayHandle *output, int loud);
//const void*	test_decode(const void *src_start, const void *src_end, int *dst, int iw, int ih, int loud);

int				abac0a_encode(const unsigned char *src, int count, int bytestride, ArrayHandle *output, int loud);
const void*		abac0a_decode(const void *src_start, const void *src_end, unsigned char *dst, int count, int bytestride, int loud);

int				abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud);
const void*		abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud);

#endif