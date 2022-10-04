#pragma once
#ifndef AWM_AC_H
#define AWM_AC_H
#include"awm_util.h"

int				abac0a_encode(const unsigned char *src, int count, int bytestride, ArrayHandle *output, int loud);
const void*		abac0a_decode(const void *src_start, const void *src_end, unsigned char *dst, int count, int bytestride, int loud);

int				abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud);
const void*		abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud);

#endif