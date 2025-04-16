#ifdef _MSC_VER
#pragma once
#endif
#ifndef INC_CODEC_H
#define INC_CODEC_H
#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifdef __cplusplus
extern "C"
{
#endif
	
	
typedef enum _CodecID
{
	CODEC_INVALID,
	CODEC_PPM,
	CODEC_PGM,
	CODEC_C01,
} CodecID;
int header_read(const unsigned char *src, int len, int *iw, int *ih, CodecID *codec);//returns header size
int compare_bufs_8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud);

int c01_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT AC (mix4)
int c02_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT binary
int c03_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT binary (A2/WG4)
int c04_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT Golomb-Rice WP3
int c05_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT Golomb-Rice 3/7 simple preds
int c06_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT binary
int c07_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT o0 (binary) AC speed test
int c08_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST disk AC5 test
int c09_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST J2K CG o0 disk symbol ANS
int c10_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST J2K CG o0 disk symbol AC
int c11_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST SubG CG o0 disk nibble AC
int c12_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST RCT o1 binary AC
int c13_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT WG4_12
int c14_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT WG4_8 binary AC
int c15_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT WG4_8 binary AC
int c16_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST SubG CG AC o0 separate loops
int c17_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST SubPrev GR (speed test)
int c18_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT SubG AC/GR
int c19_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST CG3D static o0 AC		C10 is better
int c20_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST/3T AVX2 3A3 static o1 rANS
int c21_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST  deferred WP  snapshot-CDF  AVX2 AC/RC 12 prob bits
int c22_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST/3T AVX2 484 static o1 rANS
int c23_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST/3T AVX2 blockwise Shannon, static o1 rANS
int c24_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT AC mix4
int c25_codec(const char *srcfn, const char *dstfn, int nthreads0);//MT AC mix8
int c26_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST ABAC
int c27_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST GR
int c28_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST AC with GR context
int c29_codec(const char *srcfn, const char *dstfn, int nthreads0);//C29 ST interleaved AVX2 cRCT WG4/CG GRctx static-rANS
int c30_codec(const char *srcfn, const char *dstfn, int nthreads0);//ST OLS
int c31_codec(const char *srcfn, const char *dstfn, int nthreads0);//cross-platform C29, slower
int c32_codec(const char *srcfn, const char *dstfn, int nthreads0);//like C29 but 16 coders

	
#ifdef __cplusplus
}
#endif
#endif
