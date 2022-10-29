#pragma once
#ifndef AWM_AC_H
#define AWM_AC_H
#include"awm_util.h"

//	#define		HAVE_AVX2

//color transforms
void	YCoCg_f32_fwd(const int *image, size_t size, float *bufY, float *bufCo, float *bufCg);
void	YCoCg_f32_inv(int *image, size_t size, const float *bufY, const float *bufCo, const float *bufCg);

void	YCoCg_i32_fwd(int *image, int count);
void	YCoCg_i32_inv(int *image, int count);

//spatial transforms
void	init_fwdDCT(float *matrix, int lgdim);
void	init_invDCT(float *matrix, int lgdim);
void	apply_DCT_1D(const float *matrix, float *data, int count, int stride, float *temp);
void	apply_DCT_2D(const float *hmatrix, const float *vmatrix, float *data, int bw, int bh, float *temp);

void	dwt2_2d_fwd(float *buffer, int bw, int bh, int nstages);
void	dwt2_2d_inv(float *buffer, int bw, int bh, int nstages);


//entropy coders

//compression testbench
typedef void (*Predictor)(const int *data, int iw, int ih, int bpp, int x, int y, unsigned short *hit_hist, unsigned short *ret_prob);//should return 16-bit fixed-point probability (*bpp), must use causal neighbors
int			test_encode(const int *src, int iw, int ih, int bpp, Predictor predictor, ArrayHandle *output, int loud);

int			abac0a_encode(const unsigned char *src, int count, int bytestride, ArrayHandle *output, int loud);
const void*	abac0a_decode(const void *src_start, const void *src_end, unsigned char *dst, int count, int bytestride, int loud);

int			abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud);
const void*	abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud);

#endif