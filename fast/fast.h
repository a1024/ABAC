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
int image_load(const char *fn, Image *image);
int image_save8(const char *fn, Image const *image);
int image_copy(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_copy_nodata(Image *dst, Image const *src);//dst must be initialized to zero, or another image
int image_snapshot8(Image const *image);
void image_clear(Image *image);
size_t image_getBMPsize(Image const *image);


//RCTs
void rct_JPEG2000_32(Image *image, int fwd);


//predictors
void pred_clampgrad(Image *src, int fwd, char *depths);


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

//deferred summation (8-bit only)
int	f05_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f05_encode(SRC, DATA, LOUD)		f05_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f05_decode(CBUF, CSIZE, DST, LOUD)	f05_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//adaptive 5+4-bit model (8-bit only, AVX2)
int	f06_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f06_encode(SRC, DATA, LOUD)		f06_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f06_decode(CBUF, CSIZE, DST, LOUD)	f06_codec(0, 0, CBUF, CSIZE, DST, LOUD)

//adaptive 8-bit model, MA, AVX2
int	f07_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud);
#define f07_encode(SRC, DATA, LOUD)		f07_codec(SRC, DATA, 0, 0, 0, LOUD)
#define f07_decode(CBUF, CSIZE, DST, LOUD)	f07_codec(0, 0, CBUF, CSIZE, DST, LOUD)


#ifdef __cplusplus
}
#endif
#endif//INC_E2_H
