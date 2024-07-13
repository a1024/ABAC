#pragma once
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

int c01_codec(const char *srcfn, const char *dstfn);
int c02_codec(const char *srcfn, const char *dstfn);
int c03_codec(const char *srcfn, const char *dstfn);

	
#ifdef __cplusplus
}
#endif
#endif
