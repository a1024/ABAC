#pragma once
#ifndef INC_LOSSY_H
#define INC_LOSSY_H
#include"util.h"


void measure_distortion(const unsigned char *b1, const unsigned char *b2, int iw, int ih, double *ret_rmse, double *ret_psnr);

void cvt_u8_ps(float *dst, const unsigned char *src, int res, int srcstride, int dststride);
void cvt_ps_u8(unsigned char *dst, const float *src, int res, int srcstride, int dststride);
void dnsample(const float *src, int iw, int ih, float *dst);//src & dst can be the same buffer
void upsample(const float *src, int iw, int ih, float *dst);//src & dst CANNOT be the same buffer


void colortransform_ycocg_fwd(char *buf, int iw, int ih);
void colortransform_ycocg_inv(char *buf, int iw, int ih);
void colortransform_ycocb_fwd(char *buf, int iw, int ih);
void colortransform_ycocb_inv(char *buf, int iw, int ih);

void ycocb_fwd_subsample_separate(const unsigned char *buf, int iw, int ih, short *luma, short *co, short *cb);
void ycocb_inv_upsample_separate(const short *luma, const short *co, const short *cb, int iw, int ih, unsigned char *buf);
void ycocb_inv_upsample_i8(const short *luma, const short *co, const short *cb, int iw, int ih, unsigned char *buf);

void ycocb_fwd_subsample_ps(const unsigned char *buf, int iw, int ih, float *luma, float *co, float *cb);
void ycocb_inv_upsample_ps(const float *luma, const float *co, const float *cb, int iw, int ih, unsigned char *buf);

//void colortransform_ycocb_ps_fwd(float *c0, float *c1, float *c2, int iw, int ih, int stride);
//void colortransform_ycocb_ps_inv(float *c0, float *c1, float *c2, int iw, int ih, int stride);


void DCT2_8x8_ps_buf(float *buf, int iw, int ih);//fwd
void DCT3_8x8_ps_buf(float *buf, int iw, int ih);//inv

//FCT-II/III, O(N*lg(N)), POT size (uses FFT size N)
typedef struct FCT1D_PS_ParamsStruct
{
	int lgsize;
	int *fctp;
	float *re, *im, *rew4N, *imw4N;
} FCT1D_PS_Params;
void FCT1D_ps_gen(int lgsize, FCT1D_PS_Params *p);//lgsize up to 16 (size up to 65536)
void FCT1D_ps_free(FCT1D_PS_Params *p);
void FCT1D_ps_fwd(FCT1D_PS_Params *p, float *data, int stride);
void FCT1D_ps_inv(FCT1D_PS_Params *p, float *data, int stride);

//for quantization: either provide qmatrix, or will mul/div by [de]quantize value before rounding to int
void DCT2_8x8_i16_buf(const short *src, int iw, int ih, short *dst, float quantize, const unsigned char *qmatrix);
void DCT3_8x8_i16_buf(const short *src, int iw, int ih, short *dst, float dequantize, const unsigned char *qmatrix);

void quantize(float *buf, int iw, int ih, const float *qmatrix);


//DWT
typedef struct DWTSizeStruct
{
	unsigned short w, h;
} DWTSize;
ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override);//calculate dimensions of each DWT stage in descending order
void dwt2d_cdf97_fwd(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, short *temp, int qfactor);
void dwt2d_cdf97_inv(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, short *temp, int qfactor);

void cdf97ps2d_fwd(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor);
void cdf97ps2d_inv(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor);

void lg53ps2d_fwd(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor);
void lg53ps2d_inv(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor);


void quantize_dwt_i16_fwd(short *buf, int llw, int llh, int iw, int ih, int Q);
void quantize_dwt_i16_inv(short *buf, int llw, int llh, int iw, int ih, int Q);

void quantize_dwt_ps_i16(const float *src, short *dst, int llw, int llh, int iw, int ih, float Q, float Z);//dead zone Z is at least 1
void dequantize_dwt_i16_ps(const short *src, float *dst, int llw, int llh, int iw, int ih, float Q, float Z);


//tests

//int t44_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T44 pseudo-JPEG float		dims multiple of 16		X  need to use short instead of float

int t45_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T45 pseudo JPEG (slow DWT as complex FFT)
int t45_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud);

//int t46_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T46 pseudo-J2K fixed precision	X  DWT has checkerboard artifacts
//int t46_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud);

int t47_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T47 pseudo-J2K floa
int t47_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud);


#endif//INC_LOSSY_H
