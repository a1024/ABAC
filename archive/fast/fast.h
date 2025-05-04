#pragma once
#ifndef INC_FAST_H
#define INC_FAST_H
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif


typedef struct ImageStruct
{
	int iw, ih, nch, depth;
	short *data;
} Image;
int compare_bufs_16(const short *b1, const short *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);
int compare_bufs_8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);
int image_load(const char *fn, Image *image);
int load_dng(const char *fn, Image *image);
unsigned char* image_export8(Image const *src);
int image_save8(const char *fn, Image const *image);
int image_save_native(const char *fn, Image const *image);
int image_save_ppm(const char *fn, Image const *image);
int image_copy(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_copy_nodata(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_snapshot8(Image const *image);
void image_clear(Image *image);
size_t image_getBMPsize(Image const *image);


typedef struct LSIMHeaderStruct
{
	int iw, ih, nch;
	char depth;
	int codec_id;
} LSIMHeader;
size_t lsim_writeheader(ArrayHandle *dst, int iw, int ih, int nch, char depth, int codec_id);//returns number of bytes written
size_t lsim_readheader(const unsigned char *src, size_t srclen, LSIMHeader *dst);//returns number of bytes read
int image_from_lsimheader(Image *dst, LSIMHeader const *src);//dst must be initialized to zero


//RCTs
void rct_JPEG2000_32(Image *image, int fwd);
//void rct_JPEG2000_32_avx2(Image *image, int fwd);


//predictors
void pred_clampgrad(Image *src, int fwd, const char *depths);
void pred_clampgrad_fast(Image *src, int fwd, const char *depths);
void pred_wp_deferred(Image *src, int fwd);
void pred_simd(Image *src, int fwd, const char *depths);

void pred_PU(Image *src, int fwd);


unsigned short* apply_palette_fwd(Image const *src, Image *dst, int *nvals);//nvals[4]
void apply_palette_inv(Image *dst, const unsigned short *palette, int *nvals);


//entropy estimators
void calc_csize(Image const *src, const char *depths, double *ret_csizes, double *ret_invCR);//depths[4];  ret[5] TRGBA
void calc_csize_vlc(Image const *src, const char *depths, double *ret_csizes, double *ret_csizes_vlc);
void calc_csize_bin(Image const *src, const char *depths, double *ret_csizes);
size_t calc_csize_ABAC(Image const *src, const char *depths);


//codecs

//Golomb-Rice coder
int	f01_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f01_encode(SRC, DATA, LOUD)		f01_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f01_decode(CBUF, CSIZE, DST, LOUD)	f01_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//BYPASS coder, for performance evaluation of entropy coders
int	f02_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f02_encode(SRC, DATA, LOUD)		f02_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f02_decode(CBUF, CSIZE, DST, LOUD)	f02_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//planar ABAC
int	f03_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f03_encode(SRC, DATA, LOUD)		f03_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f03_decode(CBUF, CSIZE, DST, LOUD)	f03_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//joint-bit	X
int	f04_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f04_encode(SRC, DATA, LOUD)		f04_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f04_decode(CBUF, CSIZE, DST, LOUD)	f04_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F05 deferred summation (8-bit only)
int	f05_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f05_encode(SRC, DATA, LOUD)		f05_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f05_decode(CBUF, CSIZE, DST, LOUD)	f05_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F06 adaptive 5+4-bit model (8-bit only, AVX2)
int	f06_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f06_encode(SRC, DATA, LOUD)		f06_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f06_decode(CBUF, CSIZE, DST, LOUD)	f06_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F07-20240412-CDFmixSIMD	adaptive 8-bit model, MA, AVX2
int	f07_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f07_encode(SRC, DATA, LOUD)		f07_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f07_decode(CBUF, CSIZE, DST, LOUD)	f07_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F08-adaptive-5+4-deferred-wp		8-bit only, AVX2
int	f08_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f08_encode(SRC, DATA, LOUD)		f08_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f08_decode(CBUF, CSIZE, DST, LOUD)	f08_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F09 ABAC
#define F09_NCTX 37
extern const char *f09_prednames[F09_NCTX];
int	f09_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f09_encode(SRC, DATA, LOUD)		f09_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f09_decode(CBUF, CSIZE, DST, LOUD)	f09_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F10 video test
int f10_mptest(const char *path);

//F11 selects 3 out of 30 channels
int	f11_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f11_encode(SRC, DATA, LOUD)		f11_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f11_decode(CBUF, CSIZE, DST, LOUD)	f11_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//F12 conditional histogram
int f12_statstest(const char *path);

//F13 select 3 from 18 channels, dedicated loops, deferred summation, average of 4 quantized CDFs
int f13_encode(Image const *src, ArrayHandle *data, int loud);
int f13_decode(const unsigned char *data, size_t len, Image *dst, int loud);

//	F14 select 3 from 20 channels, code 4-bit symbols, average of 2 CDFs, AVX2
int	f14_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f14_encode(SRC, DATA, LOUD)		f14_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f14_decode(CBUF, CSIZE, DST, LOUD)	f14_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F15 select 3 from 23 channels, code 4-bit symbols, average of 2 CDFs, AVX2, single-pass decoder
int	f15_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f15_encode(SRC, DATA, LOUD)		f15_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f15_decode(CBUF, CSIZE, DST, LOUD)	f15_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F16 select 3 from 23 channels, code 4-bit symbols, average of 2 CDFs, AVX2, single-pass decoder, uint16 CDFs
int	f16_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f16_encode(SRC, DATA, LOUD)		f16_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f16_decode(CBUF, CSIZE, DST, LOUD)	f16_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F17 faster, with minimum regard for efficiency
int	f17_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f17_encode(SRC, DATA, LOUD)		f17_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f17_decode(CBUF, CSIZE, DST, LOUD)	f17_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F18 array-less adaptive coder - fixed precision
int	f18_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f18_encode(SRC, DATA, LOUD)		f18_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f18_decode(CBUF, CSIZE, DST, LOUD)	f18_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F19 ANS-CDF2sym
int	f19_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f19_encode(SRC, DATA, LOUD)		f19_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f19_decode(CBUF, CSIZE, DST, LOUD)	f19_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F20 Golomb-Rice, SSE2
int	f20_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f20_encode(SRC, DATA, LOUD)		f20_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f20_decode(CBUF, CSIZE, DST, LOUD)	f20_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F21 Golomb-Rice, adaptively selects 3 among N channels
int	f21_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f21_encode(SRC, DATA, LOUD)		f21_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f21_decode(CBUF, CSIZE, DST, LOUD)	f21_codec(0, 0, CBUF, CSIZE, DST, LOUD)

void qoi_test(Image const *src, size_t *csize2, double *enc, double *dec, int *error, int loud);

//	F22 Golomb-Rice, AVX-512 (emulated)
int	f22_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f22_encode(SRC, DATA, LOUD)		f22_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f22_decode(CBUF, CSIZE, DST, LOUD)	f22_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F23 lossless multithreaded block-based  [[NOT Pareto front]]  superseded by C05
int	f23_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f23_encode(SRC, DATA, LOUD)		f23_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f23_decode(CBUF, CSIZE, DST, LOUD)	f23_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F24 lossless multithreaded block-based
int	f24_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f24_encode(SRC, DATA, LOUD)		f24_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f24_decode(CBUF, CSIZE, DST, LOUD)	f24_codec(0, 0, CBUF, CSIZE, DST, LOUD)
void	f24_curiosity();

//	F25 multithreaded, adaptive binary arithmetic coder
int	f25_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f25_encode(SRC, DATA, LOUD)		f25_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f25_decode(CBUF, CSIZE, DST, LOUD)	f25_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F26 multithreaded, 76 RCTs, 3 spatial predictors  [[NOT Pareto front]]  superseded by C01mix4
int	f26_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f26_encode(SRC, DATA, LOUD)		f26_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f26_decode(CBUF, CSIZE, DST, LOUD)	f26_codec(0, 0, CBUF, CSIZE, DST, LOUD)
void	f26_curiosity();

//	F27 multithreaded decorrelation, single-threaded entropy coding
int	f27_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f27_encode(SRC, DATA, LOUD)		f27_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f27_decode(CBUF, CSIZE, DST, LOUD)	f27_codec(0, 0, CBUF, CSIZE, DST, LOUD)
void	f27_curiosity();

//	F28 ORCT from T47
int	f28_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f28_encode(SRC, DATA, LOUD)		f28_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f28_decode(CBUF, CSIZE, DST, LOUD)	f28_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F29
int	f29_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f29_encode(SRC, DATA, LOUD)		f29_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f29_decode(CBUF, CSIZE, DST, LOUD)	f29_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F30 JPEG2000, CG, Q, o2, AC
int	f30_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f30_encode(SRC, DATA, LOUD)		f30_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f30_decode(CBUF, CSIZE, DST, LOUD)	f30_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F31 color NBLIC by WangXuan95
int	f31_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f31_encode(SRC, DATA, LOUD)		f31_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f31_decode(CBUF, CSIZE, DST, LOUD)	f31_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F32 RAW codec
int	f32_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f32_encode(SRC, DATA, LOUD)		f32_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f32_decode(CBUF, CSIZE, DST, LOUD)	f32_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//	F33  MT F26_RCT OLS4 AC_NPOT
int	f33_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f33_encode(SRC, DATA, LOUD)		f33_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f33_decode(CBUF, CSIZE, DST, LOUD)	f33_codec(0, 0, CBUF, CSIZE, DST, LOUD)


#ifdef __cplusplus
}
#endif
#endif//INC_E2_H
