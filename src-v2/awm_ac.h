#pragma once
#ifndef AWM_AC_H
#define AWM_AC_H
#include"awm_util.h"

int				abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud);
const void*		abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud);

#endif