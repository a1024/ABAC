#pragma once
#ifndef INC_E2_H
#define INC_E2_H
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif


//	#define ALLOW_SRAND
//	#define ALLOW_VULKAN
//	#define ALLOW_OPENCL


void calc_histogram(const unsigned char *buf, ptrdiff_t bytesize, ptrdiff_t stride, int *hist);


unsigned char* image_load(const char *filename, int *iw, int *ih);
int image_save_png_rgba8(const char *filename, const unsigned char *image, int iw, int ih);


//	#define ENABLE_GUIDE//debug

//lossless tools
int compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward, int loud);
void compare_bufs_ps(float *b1, float *b0, int iw, int ih, const char *name, int backward);

size_t test16_encode(const unsigned char *src, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, ArrayHandle *data, int loud, int *csizes);
int    test16_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, unsigned char *buf);


void t25_normalize_histogram(const unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF);
int t25_encode(const unsigned char *src, int iw, int ih, int *blockw, int *blockh, int use_ans, ArrayHandle *data, int loud);
int t25_decode(const unsigned char *data, size_t srclen, int iw, int ih, int *blockw, int *blockh, int use_ans, unsigned char *buf, int loud);


typedef struct T26ParamsStruct
{
	unsigned char
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} T26Params;
int t26_encode(const unsigned char *src, int iw, int ih, T26Params const *params, int use_ans, ArrayHandle *data, int loud);
int t26_decode(const unsigned char *data, size_t srclen, int iw, int ih, T26Params const *params, int use_ans, unsigned char *buf, int loud);

//int t27_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T27 X
int t28_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T28: Bayesian inference (histogram to heap)		first good result
int t29_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T29: alpha, slightly better than T28
//int t30_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T30 X  bad so far
int t31_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T31: Adaptive Bayesian inference
//int t32_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T32: Joint adaptive Bayesian inference (needs 128MB RAM)		X  bad
//int t33_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T33: Adaptive Bayesian inference with circular buffer		X  bad
int t34_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);

int t35_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T35: Combines spatial transform with entropy coding
int t35_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

//int t36_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T36: stretch & squeeze	X
//int t36_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t37_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T37: Fixed array as binary tree predictor
int t37_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t38_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T38: Single simple bit predictor
int t38_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T39: Multiple estimators for all maps		NOT GENERALIZED
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

//int t40_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T40 Random generated predictors
//int t40_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

//int t41_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T41 SIMD ABAC
//int t41_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T42: T39 with 'custom2' filter		NOT GENERALIZED
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);
//void t42_explore(void *ctx0);
//void t42_freectx(void **ctx);

//int t43_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T43: Wisdom of the crowd
//int t43_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);


//nch:   must be from {1, 2, 3, 4}.
//depth: must be from [1~16]. If depth<=8, data must be in bytes, otherwise data must be in little-endian uint16's (shorts).
//data:  must be unsigned integers shifted leftmost. For example:
//	A 14-bit subpixel must be stored like this: 0bXXXX_XXXX_XXXX_XX00
//	A 5-bit integer must be stored like this: 0bXXXX_X000
int slic_encode(int iw, int ih, int nch, int depth, const void *pixels, ArrayHandle *data, int loud);
void* slic_decode(const void *data, int len, int *iw, int *ih, int *nch, int *depth);


//nch:   Also the pixel stride. Must be from 1 to 4. The channels are interleaved and packed.
//depth: Must be from [1~16]. If depth<=8, data must be in bytes, otherwise data must be in little-endian uint16's (shorts).
//src:   Must be unsigned integers shifted leftmost. For example:
//	A 5-bit subpixel must be stored like this: 0bXXXX_X000
//	A 14-bit subpixel must be stored like this: 0bXXXX_XXXX_XXXX_XX00
//ret_dummy_alpha:  Tells if image has redundant alpha channel
unsigned char* slic2_encode(int iw, int ih, int nch, int depth, const void *src, int *ret_size);
void*          slic2_decode(const unsigned char *data, int len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha);

//I/O wrappers, return FALSE on error
int   slic2_save(const char *filename, int iw, int ih, int nch, int depth, const void *src);
void* slic2_load(const char *filename, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha);


//unsigned char* slic3_encode(int iw, int ih, int nch, int depth, const void *src, int *ret_size);
//void*          slic3_decode(const unsigned char *src, int len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha);

//T44 paq8pxd
int t44_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t44_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

//void print_ma_test(int testtype);
size_t ma_test(const unsigned char *src, int iw, int ih, int enable_RCT_MA, int RCTtype, int enable_rounding, int loud);

//T45 CALIC
int t45_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);
int t45_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);
double calc_bitsize(unsigned *CDF, int nlevels, int sym);




//transforms
extern short *g_param_ptr;
void apply_transforms_fwd(unsigned char *buf, int bw, int bh);
void apply_transforms_inv(unsigned char *buf, int bw, int bh);

void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);

void pack3_fwd(char *buf, int res);
void pack3_inv(char *buf, int res);
int get_nch(const char *buf, int res);//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}

void colortransform_YCoCg_R_fwd(char *buf, int iw, int ih);
void colortransform_YCoCg_R_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v0_fwd(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v0_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v1_fwd(char *buf, int iw, int ih);//like YCoCg but with green & blue swapped
void colortransform_YCbCr_R_v1_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v2_fwd(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v2_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v3_fwd(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v3_inv(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v4_fwd(char *buf, int iw, int ih);
void colortransform_YCbCr_R_v4_inv(char *buf, int iw, int ih);
void colortransform_JPEG2000_fwd(char *buf, int iw, int ih);
void colortransform_JPEG2000_inv(char *buf, int iw, int ih);
void colortransform_subgreen_fwd(char *buf, int iw, int ih);
void colortransform_subgreen_inv(char *buf, int iw, int ih);


	#define PW2_NPRED 20	//63
	#define PW2_NPARAM (PW2_NPRED+11)
//	#define PW2_NPRED 22
//	#define PW2_NPARAM (PW2_NPRED+8)
extern double pw2_errors[PW2_NPRED];
extern short pw2_params[PW2_NPARAM*3];
void pred_w2_opt_v2(const char *buf2, int iw, int ih, short *params, int loud);
void pred_w2_apply(char *buf, int iw, int ih, short *allparams, int fwd);

void pred_opt_printparam();
void pred_opt_opt_v6(const char *buf2, int iw, int ih, int loud);//multi-threaded grid
void pred_opt_apply(char *buf, int iw, int ih, int fwd);


extern short jxlparams_i16[33];
//extern double jxlpred_params[33];
void pred_jxl_prealloc(const char *src, int iw, int ih, int kc, const short *params, int fwd, char *dst, int *temp_w10);
void pred_jxl_opt_v2(const char *buf2, int iw, int ih, short *params, int loud);
void pred_jxl_apply(char *buf, int iw, int ih, short *allparams, int fwd);

void pred_grad_fwd(char *buf, int iw, int ih, int nch, int bytestride);
void pred_grad_inv(char *buf, int iw, int ih, int nch, int bytestride);
void grad_explore(const unsigned char *buf, int iw, int ih);


//	#define CUSTOM_TRAIN_ON_DOUBLES
	#define CUSTOM_USE_MULHRS

#define CUSTOM_REACH 2
#define CUSTOM_REACH_E 2
#define CUSTOM_NNB (CUSTOM_REACH*(CUSTOM_REACH+1)*2)
#define CUSTOM_NNB_E (CUSTOM_REACH_E*(CUSTOM_REACH_E+1)*2)
#define CUSTOM_NPARAMS (CUSTOM_NNB+CUSTOM_NNB_E)
double opt_custom(const char *buf, int iw, int ih, int kc, int niter, short *params, int loud);
float opt_custom_v2(const char *buf, int iw, int ih, int kc, int niter, short *params, float loss0, int loud);


#ifdef ALLOW_VULKAN
int init_vk();
#endif


#ifdef ALLOW_OPENCL
int init_cl(const char *searchpath);
#endif




void save_32bit(const char *filename, const int *buf, int iw, int ih, int nch, int saveas8bit);
void save_16bit(const char *filename, const short *buf, const short *sub_b2, int iw, int ih, int nch, int val_offset, int val_shift, int saveas8bit);
void save_mono8(const char *filename, unsigned char *buf, int iw, int ih, int stride);
void save_channel(unsigned char *buf, int iw, int ih, int stride, int val_offset, const char *format, ...);


#ifdef __cplusplus
}
#endif
#endif//INC_E2_H
