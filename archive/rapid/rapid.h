#pragma once
#ifndef INC_FAST_H
#define INC_FAST_H
#include"util.h"
#ifdef __cplusplus
extern "C"
{
#endif


const char* clerr2str(int error);

int r01_codec(const char *srcfn, const char *dstfn, int fwd, int loud);


#ifdef __cplusplus
}
#endif
#endif//INC_E2_H
