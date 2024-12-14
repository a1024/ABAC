#ifdef _MSC_VER
#pragma once
#endif
#ifndef INC_CODEC_H
#define INC_CODEC_H
#include"util.h"
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

int c01_codec(const char *srcfn, const char *dstfn);//MT AC (mix4)
int c02_codec(const char *srcfn, const char *dstfn);//MT binary
int c03_codec(const char *srcfn, const char *dstfn);//MT binary (A2/WG4)
int c04_codec(const char *srcfn, const char *dstfn);//MT Golomb-Rice WP3
int c05_codec(const char *srcfn, const char *dstfn);//MT Golomb-Rice 3/7 simple preds
int c06_codec(const char *srcfn, const char *dstfn);//binary
int c07_codec(const char *srcfn, const char *dstfn);//o0 (binary) AC speed test
int c08_codec(const char *srcfn, const char *dstfn);//disk AC5 test
int c09_codec(const char *srcfn, const char *dstfn);//J2K CG o0 disk symbol ANS ST
int c10_codec(const char *srcfn, const char *dstfn);//J2K CG o0 disk symbol AC ST
int c11_codec(const char *srcfn, const char *dstfn);//SubG CG o0 disk nibble AC ST
int c12_codec(const char *srcfn, const char *dstfn);//RCT o1 binary AC ST
int c13_codec(const char *srcfn, const char *dstfn);//MT WG4_12
int c14_codec(const char *srcfn, const char *dstfn);//MT WG4_8 binary AC
int c15_codec(const char *srcfn, const char *dstfn);//MT WG4_8 binary AC
int c16_codec(const char *srcfn, const char *dstfn);//ST SubG CG AC o0 separate loops
int c17_codec(const char *srcfn, const char *dstfn);//ST SubPrev GR (speed test)
int c18_codec(const char *srcfn, const char *dstfn);//MT SubG AC/GR
int c19_codec(const char *srcfn, const char *dstfn);//ST CG3D static o0 AC		C10 is better
int c20_codec(const char *srcfn, const char *dstfn);//ST/3T AVX2 3A3 static o1 rANS
int c21_codec(const char *srcfn, const char *dstfn);//ST  deferred WP  snapshot-CDF  AVX2 AC/RC 12 prob bits
int c22_codec(const char *srcfn, const char *dstfn);//ST/3T AVX2 484 static o1 rANS
int c23_codec(const char *srcfn, const char *dstfn);//ST/3T AVX2 blockwise Shannon, static o1 rANS
int c24_codec(const char *srcfn, const char *dstfn);//MT AC mix4
int c25_codec(const char *srcfn, const char *dstfn);//MT AC mix8

	
#ifdef __cplusplus
}
#endif
#endif
