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
void compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward);
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

int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T39: Multiple estimators for all maps
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t40_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T40 Random generated predictors
int t40_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);

int t41_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud);//T41 SIMD ABAC
int t41_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud);




//transforms
extern short *g_param_ptr;
void apply_transforms_fwd(unsigned char *buf, int bw, int bh);
void apply_transforms_inv(unsigned char *buf, int bw, int bh);

void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);

void colortransform_ycocg_fwd(char *buf, int iw, int ih);
void colortransform_ycocg_inv(char *buf, int iw, int ih);
void colortransform_ycocb_fwd(char *buf, int iw, int ih);//like YCoCg but with green & blue swapped
void colortransform_ycocb_inv(char *buf, int iw, int ih);


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

void pred_grad_fwd     (char *buf, int iw, int ih, int nch, int bytestride);
void pred_grad_inv     (char *buf, int iw, int ih, int nch, int bytestride);


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
