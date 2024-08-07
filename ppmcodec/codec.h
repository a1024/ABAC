﻿#pragma once
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

int c01_codec(const char *srcfn, const char *dstfn);//binary/AC
int c02_codec(const char *srcfn, const char *dstfn);//binary
int c03_codec(const char *srcfn, const char *dstfn);//binary
int c04_codec(const char *srcfn, const char *dstfn);//Golomb-Rice
int c05_codec(const char *srcfn, const char *dstfn);//Golomb-Rice
int c06_codec(const char *srcfn, const char *dstfn);//binary
int c07_codec(const char *srcfn, const char *dstfn);//o0 AC speed test
int c08_codec(const char *srcfn, const char *dstfn);//FILE-only AC5 test
int c09_codec(const char *srcfn, const char *dstfn);//Pei09 CG o0 ST AVX2 ANS

	
#ifdef __cplusplus
}
#endif
#endif
