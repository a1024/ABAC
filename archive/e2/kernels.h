#ifdef __OPEN_CL__
#pragma OPENCL EXTENSION cl_khr_fp16 : enable
#else
#include<stdio.h>
#define __kernel
#define __global
#define __constant
#define half float
#endif
//__kernel void 
