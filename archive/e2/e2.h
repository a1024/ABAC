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
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4200)//no default-constructor for struct with zero-length array
#endif


//	#define ALLOW_SRAND
//	#define ALLOW_VULKAN
//	#define ALLOW_OPENCL


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
int image_snapshot(Image const *image);
int image_save_native(const char *fn, Image const *image, int override_alpha);
Image* image_alloc(int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth, int initialize, int initval);
Image* image_from_uint8(const unsigned char *src, int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth);
Image* image_from_uint16(const unsigned short *src, int iw, int ih, int nch, char *src_depths, char *dst_depths);
void image_export_uint8(Image const *image, unsigned char **dst, int override_alpha);
void image_export_uint16(Image const *image, unsigned short **dst, int override_alpha, int big_endian);
double image_getBMPsize(Image const *image);
size_t image_getbufsize(Image const *image);
void image_copy_nodata(Image **dst, Image const *src);//dst must be 0 or a valid pointed
void image_copy(Image **dst, Image const *src);//dst must be 0 or a valid pointed
void calc_depthfromdata(int *image, int iw, int ih, char *depths, const char *src_depths);
void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8);
double calc_entropy(const int *hist, int nlevels, int sum);
int calc_maxdepth(Image const *image, int *inflation);
void calc_depthfromdata(int *image, int iw, int ih, char *depths, const char *src_depths);
//double calc_csize_from_hist(int *hist, int nlevels, double *ret_usize);//works only with unit-increment histograms initialized with ones
void calc_csize(Image const *image, double *csizes);


typedef struct LSIMHeaderStruct
{
	int iw, ih, nch;
	char depth[4];
	int codec_id;
} LSIMHeader;
size_t lsim_writeheader(ArrayHandle *dst, int iw, int ih, int nch, const char *depths, int codec_id);//returns number of bytes written
size_t lsim_readheader(const unsigned char *src, size_t srclen, LSIMHeader *dst);//returns number of bytes read
void image_from_lsimheader(Image **dst, LSIMHeader const *src);//dst must be 0 or a valid pointed


//void calc_histogram(const unsigned char *buf, ptrdiff_t bytesize, ptrdiff_t stride, int *hist);


//unsigned char* image_load(const char *filename, int *iw, int *ih);
int image_save_png_rgba8(const char *filename, const unsigned char *image, int iw, int ih);


//	#define ENABLE_GUIDE//debug

//lossless tools
int compare_bufs_32(const int *b1, const int *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);
int compare_bufs_uint8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward, int loud);
void compare_bufs_ps(const float *b1, const float *b0, int iw, int ih, const char *name, int backward);

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
int t45_encode(Image const *src, ArrayHandle *data, int loud);
int t45_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

double calc_bitsize(unsigned *CDF, int nlevels, int sym);


//SLIC v4
int t46_encode(Image const *src, ArrayHandle *data, int loud);
int t46_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//SLIC v5
#define RCTLIST\
	RCT(NONE)\
	RCT(SubGreen)\
	RCT(JPEG2000)\
	RCT(YCoCg_R)\
	RCT(YCbCr_R_v1)\
	RCT(A710)\
	RCT(YCbCr_R_v3)\
	RCT(YCbCr_R_v4)\
	RCT(YCbCr_R_v5)\
	RCT(YCbCr_R_v6)\
	RCT(YCbCr_R_v7)\
	RCT(Pei09)
typedef enum RCTTypeEnum
{
#define RCT(X) RCT_##X,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTType;
extern const char *rct_names[RCT_COUNT];

//T47 SLIC5
#define SLIC5_NPREDS 33
#define SLIC5_OPTIMIZE_RCT//extremely slow
typedef struct SLIC5CuriosityStruct
{
#ifdef SLIC5_OPTIMIZE_RCT
#define ORCT_NPARAMS 8//not counting permutation
	char rct_params[ORCT_NPARAMS+1];
#else
	RCTType rct;
	double rct_sizes[RCT_COUNT];
	ptrdiff_t coverage[RCT_COUNT];
#endif
	long long pred_errors[SLIC5_NPREDS];
} SLIC5Curiosity;
extern const char *slic5_prednames[SLIC5_NPREDS];
//extern const char *slic5_orct_permutationnames[6];
void orct_print_compact(const char *params);
int t47_encode(Image const *src, ArrayHandle *data, SLIC5Curiosity *curiosity, int loud);
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);
//RCTType rct_select_best(Image const *src, double *ret_csizes);
//int t47_from_ppm(const char *src, const char *dst);
//int t47_to_ppm(const char *src, const char *dst);
//void t47_analyze_preds(const char *path);

//CABAC
int t48_encode(Image const *src, ArrayHandle *data, int loud);
int t48_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

int t49_encode(Image const *src, ArrayHandle *data, int loud);
int t49_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//T51 ABAC
int t51_encode(Image const *src, ArrayHandle *data, int loud);
int t51_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);
void test_alphaVSbin(Image const *src);

//T52 pseudo-JPEG
void ct_YCbCr_lossy(Image *src, int fwd);
void resample_YUV420(Image *src, int down);
void dct_8x8(Image *dst, Image const *src, int subsampled, int *ymatrix, int *uvmatrix, int fwd);
//void jpeg_quantize(Image *src, int quality, int fwd);//quality [1, 100]
int t52_encode(Image const *src, ArrayHandle *data, int loud);
int t52_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//T53 Adaptive block
int t53_encode(Image const *src, ArrayHandle *data, int loud);
int t53_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//T54 Tries to be fast, slightly sacrificing efficiency
int t54_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define t54_encode(SRC, DATA, LOUD)		t54_codec(SRC, DATA, 0, 0, 0, LOUD)
#define t54_decode(CBUF, CSIZE, DST, LOUD)	t54_codec(0, 0, CBUF, CSIZE, DST, LOUD)
//int t54_encode(Image const *src, ArrayHandle *data, int loud);
//int t54_decode(const unsigned char *data, size_t srclen, Image *dst, int loud);

//T55 Binary
int t55_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define t55_encode(SRC, DATA, LOUD)		t55_codec(SRC, DATA, 0, 0, 0, LOUD)
#define t55_decode(CBUF, CSIZE, DST, LOUD)	t55_codec(0, 0, CBUF, CSIZE, DST, LOUD)




//transforms
extern short *g_param_ptr;
void apply_transforms_fwd(unsigned char *buf, int bw, int bh);
void apply_transforms_inv(unsigned char *buf, int bw, int bh);

void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount);

void pack3_fwd(char *buf, int res);
void pack3_inv(char *buf, int res);
int get_nch32(const int *buf, int res);//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}
int get_nch(const char *buf, int res);

void rct_JPEG2000_32(Image *image, int fwd);
void rct_JPEG2000_32_ma(int *image, int iw, int ih, char *depths, int fwd);
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


void pred_clampgrad(Image *src, int fwd, char *depths);

	#define PW2_NPRED 20	//63
	#define PW2_NPARAM (PW2_NPRED+11)
//	#define PW2_NPRED 22
//	#define PW2_NPARAM (PW2_NPRED+8)
extern double pw2_errors[PW2_NPRED];
extern short pw2_params[PW2_NPARAM*3];
//void pred_w2_opt_v2(const char *buf2, int iw, int ih, short *params, int loud);
void pred_w2_apply(char *buf, int iw, int ih, short *allparams, int fwd);

void pred_opt_printparam();
//void pred_opt_opt_v6(const char *buf2, int iw, int ih, int loud);//multi-threaded grid
void pred_opt_apply(char *buf, int iw, int ih, int fwd);


extern short jxlparams_i16[33];
//extern double jxlpred_params[33];
void pred_jxl_prealloc(const char *src, int iw, int ih, int kc, const short *params, int fwd, char *dst, int *temp_w10);
//void pred_jxl_opt_v2(const char *buf2, int iw, int ih, short *params, int loud);
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
#ifndef __GNUC__
double opt_custom(const char *buf, int iw, int ih, int kc, int niter, short *params, int loud);
float opt_custom_v2(const char *buf, int iw, int ih, int kc, int niter, short *params, float loss0, int loud);
#endif

void image_dct8_fwd(Image *image);
void image_dct8_inv(Image *image);


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



#ifdef _MSC_VER
#pragma warning(pop)
#endif
#ifdef __cplusplus
}
#endif
#endif//INC_E2_H
