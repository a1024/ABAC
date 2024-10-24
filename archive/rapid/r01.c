#include"rapid.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define CL_TARGET_OPENCL_VERSION 300
#include<CL/cl.h>
#ifdef _MSC_VER
#pragma comment(lib, "OpenCL.lib")
#endif
#include<sys/stat.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

//	#define ENABLE_GUIDE
	#define ENABLE_PROFILER		//profiler is not thread-safe

#ifdef __GNUC__
	#define KERNEL_FN "r01_kernels.h"
#else
	#define KERNEL_FN "C:/Projects/rapid/rapid/r01_kernels.h"
#endif

const char* clerr2str(int error)
{
	const char *a=0;
#define 		EC(x)		case x:a=(const char*)#x;break;
#define 		EC2(n, x)	case n:a=(const char*)#x;break;
	switch(error)
	{
	EC(CL_SUCCESS)
	EC(CL_DEVICE_NOT_FOUND)
	EC(CL_DEVICE_NOT_AVAILABLE)
	EC(CL_COMPILER_NOT_AVAILABLE)
	EC(CL_MEM_OBJECT_ALLOCATION_FAILURE)
	EC(CL_OUT_OF_RESOURCES)
	EC(CL_OUT_OF_HOST_MEMORY)
	EC(CL_PROFILING_INFO_NOT_AVAILABLE)
	EC(CL_MEM_COPY_OVERLAP)
	EC(CL_IMAGE_FORMAT_MISMATCH)
	EC(CL_IMAGE_FORMAT_NOT_SUPPORTED)
	EC(CL_BUILD_PROGRAM_FAILURE)
	EC(CL_MAP_FAILURE)
//#ifdef CL_VERSION_1_1
	EC(CL_MISALIGNED_SUB_BUFFER_OFFSET)
	EC(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_COMPILE_PROGRAM_FAILURE)
	EC(CL_LINKER_NOT_AVAILABLE)
	EC(CL_LINK_PROGRAM_FAILURE)
	EC(CL_DEVICE_PARTITION_FAILED)
	EC(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
//#endif
	EC(CL_INVALID_VALUE)
	EC(CL_INVALID_DEVICE_TYPE)
	EC(CL_INVALID_PLATFORM)
	EC(CL_INVALID_DEVICE)
	EC(CL_INVALID_CONTEXT)
	EC(CL_INVALID_QUEUE_PROPERTIES)
	EC(CL_INVALID_COMMAND_QUEUE)
	EC(CL_INVALID_HOST_PTR)
	EC(CL_INVALID_MEM_OBJECT)
	EC(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
	EC(CL_INVALID_IMAGE_SIZE)
	EC(CL_INVALID_SAMPLER)
	EC(CL_INVALID_BINARY)
	EC(CL_INVALID_BUILD_OPTIONS)
	EC(CL_INVALID_PROGRAM)
	EC(CL_INVALID_PROGRAM_EXECUTABLE)
	EC(CL_INVALID_KERNEL_NAME)
	EC(CL_INVALID_KERNEL_DEFINITION)
	EC(CL_INVALID_KERNEL)
	EC(CL_INVALID_ARG_INDEX)
	EC(CL_INVALID_ARG_VALUE)
	EC(CL_INVALID_ARG_SIZE)
	EC(CL_INVALID_KERNEL_ARGS)
	EC(CL_INVALID_WORK_DIMENSION)
	EC(CL_INVALID_WORK_GROUP_SIZE)
	EC(CL_INVALID_WORK_ITEM_SIZE)
	EC(CL_INVALID_GLOBAL_OFFSET)
	EC(CL_INVALID_EVENT_WAIT_LIST)
	EC(CL_INVALID_EVENT)
	EC(CL_INVALID_OPERATION)
	EC(CL_INVALID_GL_OBJECT)
	EC(CL_INVALID_BUFFER_SIZE)
	EC(CL_INVALID_MIP_LEVEL)
	EC(CL_INVALID_GLOBAL_WORK_SIZE)
//#ifdef CL_VERSION_1_1
	EC(CL_INVALID_PROPERTY)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_INVALID_IMAGE_DESCRIPTOR)
	EC(CL_INVALID_COMPILER_OPTIONS)
	EC(CL_INVALID_LINKER_OPTIONS)
	EC(CL_INVALID_DEVICE_PARTITION_COUNT)
//#endif
//#ifdef CL_VERSION_2_0
	EC2(-69, CL_INVALID_PIPE_SIZE)
	EC2(-70, CL_INVALID_DEVICE_QUEUE)
//#endif
//#ifdef CL_VERSION_2_2
	EC2(-71, CL_INVALID_SPEC_ID)
	EC2(-72, CL_MAX_SIZE_RESTRICTION_EXCEEDED)
//#endif
//	EC(CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR)
//	EC(CL_PLATFORM_NOT_FOUND_KHR)
	EC2(-1002, CL_INVALID_D3D10_DEVICE_KHR)
	EC2(-1003, CL_INVALID_D3D10_RESOURCE_KHR)
	EC2(-1004, CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1005, CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1006, CL_INVALID_D3D11_DEVICE_KHR)
	EC2(-1007, CL_INVALID_D3D11_RESOURCE_KHR)
	EC2(-1008, CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1009, CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR)
#ifndef __linux__
	EC2(-1010, CL_INVALID_D3D9_DEVICE_NV_or_CL_INVALID_DX9_DEVICE_INTEL)
	EC2(-1011, CL_INVALID_D3D9_RESOURCE_NV_or_CL_INVALID_DX9_RESOURCE_INTEL)
	EC2(-1012, CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV_or_CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL)
	EC2(-1013, CL_D3D9_RESOURCE_NOT_ACQUIRED_NV_or_CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL)
#endif
	EC2(-1092, CL_EGL_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1093, CL_INVALID_EGL_OBJECT_KHR)
	EC2(-1094, CL_INVALID_ACCELERATOR_INTEL)
	EC2(-1095, CL_INVALID_ACCELERATOR_TYPE_INTEL)
	EC2(-1096, CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL)
	EC2(-1097, CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL)
	EC2(-1098, CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL)
	EC2(-1099, CL_INVALID_VA_API_MEDIA_SURFACE_INTEL)
	EC2(-1101, CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL)
	case 1:a="File failure";break;//
	default:
		a="???";
		break;
	}
#undef			EC
#undef			EC2
	return a;
}
#define CHECKCL(ERR) ASSERT_MSG(!(ERR), "OpenCL error %d: %s\n", ERR, clerr2str(ERR))

#define KERNELLIST\
	KERNEL(prep_planes)\
	KERNEL(pred_planes)\
	KERNEL(calc_entropy)

typedef enum _R01Kernels
{
#define KERNEL(X) OCL_##X,
	KERNELLIST
#undef  KERNEL
	OCL_NKERNELS,
} R01Kernels;
static const char *kernelnames[]=
{
#define KERNEL(X) #X,
	KERNELLIST
#undef  KERNEL
};
typedef struct _R01Context
{
	unsigned nplatforms, ndevices;
	cl_platform_id *platforms;
	cl_device_id *devices;
	size_t maxlocalsize;
	cl_context context;
	cl_command_queue commandqueue;
	cl_program program;
	cl_kernel kernels[OCL_NKERNELS];
} R01Context;
//static R01Context* init_cl(const char *searchpath)
static R01Context* init_cl(int loud)
{
	int error;
	int pl_idx=0;
	size_t temp=0, retlen=0;
	R01Context *ctx=(R01Context*)malloc(sizeof(R01Context));
	if(!ctx)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(ctx, 0, sizeof(R01Context));

	//oclw_loadAPI();

	error=clGetPlatformIDs(0, 0, &ctx->nplatforms);		CHECKCL(error);
	ASSERT_MSG(ctx->nplatforms, "No OpenCL platforms");
	
	ctx->platforms=(cl_platform_id*)malloc(ctx->nplatforms*sizeof(cl_platform_id));
	error=clGetPlatformIDs(ctx->nplatforms, ctx->platforms, 0);	CHECKCL(error);

	ctx->ndevices=0;
	for(;pl_idx<(int)ctx->nplatforms;++pl_idx)
	{
		error=clGetPlatformInfo(ctx->platforms[pl_idx], CL_PLATFORM_VERSION, G_BUF_SIZE, g_buf, &retlen);CHECKCL(error);
		error=clGetDeviceIDs(ctx->platforms[pl_idx], CL_DEVICE_TYPE_GPU, 0, 0, &ctx->ndevices);
		if(loud)
			printf("OpenCL platform: %s, %d device(s)\n", g_buf, ctx->ndevices);
		if(ctx->ndevices)//break before pl_idx increment
			break;
	}
	CHECKCL(error);
	ASSERT_MSG(ctx->ndevices, "No OpenCL devices");

	ctx->devices=(cl_device_id*)malloc(ctx->ndevices*sizeof(cl_device_id));
	error=clGetDeviceIDs(ctx->platforms[pl_idx], CL_DEVICE_TYPE_GPU, ctx->ndevices, ctx->devices, 0);	CHECKCL(error);
	
	//get info
	if(loud)
	{
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_NAME, G_BUF_SIZE, g_buf, &retlen);	CHECKCL(error);
		printf("Device: %s\n", g_buf);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_VENDOR, G_BUF_SIZE, g_buf, &retlen);	CHECKCL(error);
		printf("Vendor: %s\n", g_buf);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &ctx->maxlocalsize, &retlen);	CHECKCL(error);
		printf("CL_DEVICE_MAX_WORK_GROUP_SIZE = %d\n", (int)ctx->maxlocalsize);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF, sizeof(size_t), &temp, &retlen);	CHECKCL(error);
		printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF = %d\n", (int)temp);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(size_t), &temp, &retlen);	CHECKCL(error);
		printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT = %d\n", (int)temp);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_ADDRESS_BITS, sizeof(size_t), &temp, &retlen);	CHECKCL(error);
		printf("CL_DEVICE_ADDRESS_BITS = %d\n", (int)temp);
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_EXTENSIONS, G_BUF_SIZE, g_buf, &retlen);	CHECKCL(error);
		for(int k=0;k<retlen;++k)
		{
			if(g_buf[k]==' ')
				g_buf[k]='\n';
		}
		printf("Extensions:\n%s\n", g_buf);
		//int cl_gl_interop=strstr(g_buf, "cl_khr_gl_sharing")!=0;
		//if(!cl_gl_interop)
		//	printf("\n\tcl_khr_gl_sharing not supported\n\n");
	}
	else
	{
		error=clGetDeviceInfo(ctx->devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &ctx->maxlocalsize, &retlen);	CHECKCL(error);
	}
	
	//create context & command queue
#ifdef CL_GL_INTEROP
	cl_context_properties properties[8]={};
	if(cl_gl_interop)
	{
		auto gl_context=eglGetCurrentContext();//changes when resuming
		auto egl_display=eglGetCurrentDisplay();
		properties[0]=CL_GL_CONTEXT_KHR,	properties[1]=(cl_context_properties)gl_context;//https://stackoverflow.com/questions/26802905/getting-opengl-buffers-using-opencl
		properties[2]=CL_EGL_DISPLAY_KHR,	properties[3]=(cl_context_properties)egl_display;
		properties[4]=CL_CONTEXT_PLATFORM,	properties[5]=(cl_context_properties)platform;
		properties[6]=0, properties[7]=0;
	}
	else
	{
		properties[0]=CL_CONTEXT_PLATFORM, properties[1]=(cl_context_properties)platform;
		properties[2]=0, properties[3]=0;
	}
	context=clCreateContext(properties, 1, &devices[0], nullptr, nullptr, &error);	CHECKCL(error);
#else
	ctx->context=clCreateContext(0, ctx->ndevices, ctx->devices, 0, 0, &error);			CHECKCL(error);
#endif
#if CL_TARGET_OPENCL_VERSION>=200
	//cl_queue_properties properties[]=
	//{
	//	CL_QUEUE_PROPERTIES, 0,
	//	CL_QUEUE_SIZE, 0,
	//	0,
	//};
	ctx->commandqueue=clCreateCommandQueueWithProperties(ctx->context, ctx->devices[0], 0, &error);	CHECKCL(error);
#else
	ctx->commandqueue=clCreateCommandQueue(ctx->context, ctx->devices[0], 0, &error);	CHECKCL(error);
#endif

	//compile kernel
	ArrayHandle srctext=load_file(KERNEL_FN, 0, 0, 0);
	//ArrayHandle filename=searchfor_file(searchpath, KERNEL_FN);
	//if(!filename)
	//{
	//	LOG_ERROR("Cannot find " KERNEL_FN);
	//	return 0;
	//}
	//ArrayHandle srctext=load_file(filename->data, 0, 0, 0);
	//if(!srctext)
	//{
	//	LOG_ERROR("Cannot open " KERNEL_FN);
	//	return 0;
	//}
	//array_free(&filename);
	//ASSERT_MSG(srctext, "Cannot open \'%s\'", srcname);
	const char *k_src=(const char*)srctext->data;
	size_t k_len=srctext->count;
	ctx->program=clCreateProgramWithSource(ctx->context, 1, (const char**)&k_src, &k_len, &error);	CHECKCL(error);
#ifdef PREC_FIXED
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__ -DPREC_FIXED=%d", PREC_FIXED);
#elif defined PREC_HALF
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__ -DPREC_HALF");
#else
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__");
#endif
	error=clBuildProgram(ctx->program, 1, ctx->devices, g_buf, 0, 0);
	if(error)
	{
		error=clGetProgramBuildInfo(ctx->program, ctx->devices[0], CL_PROGRAM_BUILD_LOG, G_BUF_SIZE, g_buf, &retlen);
		if(retlen>G_BUF_SIZE)
		{
			char *buf=(char*)malloc(retlen+10);
			error=clGetProgramBuildInfo(ctx->program, ctx->devices[0], CL_PROGRAM_BUILD_LOG, retlen+10, buf, &retlen);	CHECKCL(error);
			printf("\nOpenCL Compilation failed:\n%s\n", buf);
			free(buf);
			LOG_ERROR("Aborting");
		}
		else
			LOG_ERROR("\nOpenCL Compilation failed:\n%s\n", g_buf);
		return 0;
	}
	array_free(&srctext);

#if 1
	for(int k=0;k<OCL_NKERNELS;++k)
	{
		ctx->kernels[k]=clCreateKernel(ctx->program, kernelnames[k], &error);
		if(error)
		{
			LOG_ERROR("CL error %s: Cannot find kernel %s", clerr2str(error), kernelnames[k]);
			return 0;
		}
	}
#endif
	return ctx;
}


#define BLOCKSIZE 64


#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(GB)\
	OCH(BR)\
	OCH(R3)\
	OCH(G3)\
	OCH(B3)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)\
	OCH(R1)\
	OCH(G1)\
	OCH(B1)\
	OCH(RB)\
	OCH(GR)\
	OCH(BG)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,	0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,	0, 4)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,	4, 0)\
	RCT(R_G_B1,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,	1, 3)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,	2, 2)\
	RCT(R_G_B3,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,	3, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	-1,	4, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	-1,	0, 4)\
	RCT(R_GR_B1,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	-1,	1, 3)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	-1,	2, 2)\
	RCT(R_GR_B3,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	-1,	3, 1)\
	RCT(R_B_G1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,	3, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,	2, 2)\
	RCT(R_B_G3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,	1, 3)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	-1,	0, 4)\
	RCT(R_BR_G1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	-1,	3, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	-1,	2, 2)\
	RCT(R_BR_G3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	-1,	1, 3)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,	4, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,	0, 4)\
	RCT(G_B_R1,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,	1, 3)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,	2, 2)\
	RCT(G_B_R3,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,	3, 1)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	-1,	4, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	-1,	0, 4)\
	RCT(G_BG_R1,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	-1,	1, 3)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	-1,	2, 2)\
	RCT(G_BG_R3,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	-1,	3, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	-1,	0, 4)\
	RCT(G_RG_B1,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	-1,	3, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	-1,	2, 2)\
	RCT(G_RG_B3,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	-1,	1, 3)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,	0, 4)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,	4, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	-1,	4, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	-1,	0, 4)\
	RCT(B_RB_G1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	-1,	1, 3)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	-1,	2, 2)\
	RCT(B_RB_G3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	-1,	3, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	-1,	0, 4)\
	RCT(B_GB_R1,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	-1,	3, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	-1,	2, 2)\
	RCT(B_GB_R3,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	-1,	1, 3)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UMASK1,  VCOEFF0, VCOEFF1) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][9]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UMASK1,  VCOEFF0, VCOEFF1)\
	{YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UMASK1,  VCOEFF0, VCOEFF1},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UMASK1,  VCOEFF0, VCOEFF1) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(N, "N")\
	PRED(W, "W")\
	PRED(AV2, "(N+W)/2")\
	PRED(WG, "(gx*N+gy*W)/(gx+gy)")\
	PRED(CG, "median(N, W, N+W-NW)")\
	PRED(AV3, "[-2 [3];3 ?]/4")\
	PRED(AV4, "[-1 [4] 1;4 ?]/8")\
	PRED(AV5, "[-5 [5] 1;-1 8 ?]/8")\
	PRED(AV6, "[[-1];-5 [6] 1;-1 8 ?]/8")\
	PRED(AV9, "[1 [-2] -1;-1 -9 [10] 4;-2 16 ?]/16")\
	PRED(AVB, "[4 3 [-31] -38;7 -158 [219] 30 19;-42 243 ?]/256")
static const char *pred_desc[]=
{
#define PRED(LABEL, DESC) DESC,
	PREDLIST
#undef  PRED
};
typedef enum _PredIndex
{
#define PRED(LABEL, DESC) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL, DESC) #LABEL,
	PREDLIST
#undef  PRED
};

//https://github.com/samtools/htscodecs
unsigned char *rans_compress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);
unsigned char *rans_compress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);


#ifdef ENABLE_PROFILER//not thread-safe
#define CHECKLIST\
	CHECKPOINT(START)\
	CHECKPOINT(READ)\
	CHECKPOINT(INIT_CL)\
	CHECKPOINT(ANALYSIS)\
	CHECKPOINT(SELECT)\
	CHECKPOINT(PREDICT)\
	CHECKPOINT(ENTROPY_CODING)\
	CHECKPOINT(UNPREDICT)\
	CHECKPOINT(WRITE)
typedef enum _CheckpointType
{
#define CHECKPOINT(X) CHECK_##X,
	CHECKLIST
#undef  CHECKPOINT
	CHECK_COUNT,
} CheckpointType;
static const char *checkpointnames[]=
{
#define CHECKPOINT(X) #X,
	CHECKLIST
#undef  CHECKPOINT
};
static unsigned long long g_c=0;
static double g_t=0;
static unsigned long long g_cycles[CHECK_COUNT]={0};
static double g_seconds[CHECK_COUNT]={0};
FORCE_INLINE void prof_check(int idx)
{
	if(!g_c)
	{
		memset(g_cycles, 0, sizeof(g_cycles));
		g_c=__rdtsc();
	}
	else
	{
		unsigned long long c2=__rdtsc();
		g_cycles[idx]=c2-g_c;
		g_c=c2;
	}

	if(!g_t)
	{
		memset(g_seconds, 0, sizeof(g_seconds));
		g_t=time_sec();
	}
	else
	{
		double t2=time_sec();
		g_seconds[idx]=t2-g_t;
		g_t=t2;
	}
}
static void prof_print(void)
{
	unsigned long long csum=0;
	double ssum=0;
	int len=0;
	for(int k=0;k<CHECK_COUNT;++k)
	{
		int len2=(int)strlen(checkpointnames[k]);
		if(len<len2)
			len=len2;
		csum+=g_cycles[k];
		ssum+=g_seconds[k];
	}
	for(int k=0;k<CHECK_COUNT;++k)
	{
		if(!g_cycles[k])
			continue;
		printf("%-*s  %12.6lf sec %8.4lf%%  %16lld cycles %8.4lf%%\n",
			len, checkpointnames[k],
			g_seconds[k], g_seconds[k]*100./ssum,
			g_cycles[k], (double)g_cycles[k]*100./csum
		);
	}
	memset(g_cycles, 0, sizeof(g_cycles));
	memset(g_seconds, 0, sizeof(g_seconds));
}
#else
#define prof_check(...)
#define prof_print(...)
#endif


#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif


FORCE_INLINE int predict(unsigned char predidx, short **rows)
{
	switch(predidx)
	{
	case PRED_N:	return rows[1][+0*4*2];
	case PRED_W:	return rows[0][-1*4*2];
	case PRED_AV2:	return (rows[1][+0*4*2]+rows[0][-1*4*2]+1)>>1;
	case PRED_WG:
		{
			int
				NN	=rows[2][+0*4*2],
				NNE	=rows[2][+1*4*2],
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				WW	=rows[0][-2*4*2],
				W	=rows[0][-1*4*2];
			int gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1;
			int gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1;
			return (gx*N+gy*W)/(gx+gy);
		}
		break;
	case PRED_CG:
		{
			int
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			pred=N+W-NW;
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
	case PRED_AV3:
		{
			int
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=(3*(N+W)-2*NW+2)>>2;
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	case PRED_AV4:
		{
			int
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=(4*(N+W)+NE-NW+4)>>3;
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	case PRED_AV5:
		{
			int
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				WW	=rows[0][-2*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=W+((5*(N-NW)+NE-WW+4)>>3);
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	case PRED_AV6:
		{
			int
				NN	=rows[2][+0*4*2],
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				WW	=rows[0][-2*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=W+((6*N-5*NW+NE-NN-WW+4)>>3);
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	case PRED_AV9:
		{
			int
				NNW	=rows[2][-1*4*2],
				NN	=rows[2][+0*4*2],
				NNE	=rows[2][+1*4*2],
				NWW	=rows[1][-2*4*2],
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				WW	=rows[0][-2*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=W+((10*N-9*NW+4*NE-2*(NN+WW)-NNE+NNW-NWW+8)>>4);
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	case PRED_AVB:
		{
			int
				NNWW	=rows[2][-2*4*2],
				NNW	=rows[2][-1*4*2],
				NN	=rows[2][+0*4*2],
				NNE	=rows[2][+1*4*2],
				NWW	=rows[1][-2*4*2],
				NW	=rows[1][-1*4*2],
				N	=rows[1][+0*4*2],
				NE	=rows[1][+1*4*2],
				NEE	=rows[1][+2*4*2],
				WW	=rows[0][-2*4*2],
				W	=rows[0][-1*4*2];
			int vmax=N, vmin=W, pred;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			pred=(4*NNWW+3*NNW-31*NN-38*NNE + 7*NWW-158*NW+219*N+30*NE+19*NEE - 42*WW+243*W+128)>>8;
			if(pred<vmin)pred=vmin;
			if(pred<vmax)pred=vmax;
			return pred;
		}
		break;
	}
	return 0;
}
int r01_codec(const char *srcfn, const char *dstfn, int fwd, int loud)
{
	double elapsed=time_sec();
	unsigned long long celapsed=__rdtsc();
	ptrdiff_t srcsize=0;
	unsigned char *srcbuf=0;
	if(loud)
		prof_check(CHECK_START);
	{
		struct stat info={0};
		int ret=stat(srcfn, &info);
		if(ret)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcsize=info.st_size;
	}
	srcbuf=(unsigned char*)malloc(srcsize+16);
	if(!srcbuf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
	}
	unsigned char *srcptr=srcbuf;
	int iw=0, ih=0;
	if(fwd)
	{
		if(memcmp(srcptr, "P6\n", 3))
		{
			LOG_ERROR("Expected a PPM file");
			return 1;
		}
		srcptr+=3;

		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';

		if(*srcptr++!=' ')
		{
			LOG_ERROR("Expected a PPM file");
			return 1;
		}

		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		
		if(memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Expected a PPM file");
			return 1;
		}
		srcptr+=5;
	}
	else
	{
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Invalid image");
		return 1;
	}
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	int xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	int yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	int nblocks=xblocks*yblocks;

	if(loud)
		prof_check(CHECK_START);

	size_t filesize=0;
	if(fwd)//encode
	{
		unsigned char *image=srcptr;
		guide_save(image, iw, ih);

		//Step 0: Init OpenCL
		R01Context *ctx=init_cl(loud);
		if(!ctx)
		{
			LOG_ERROR("Failed to initialize OpenCL");
			return 1;
		}
		if(loud)
			prof_check(CHECK_INIT_CL);

		//Step 1: Analyze
		int error=0;
		size_t
			gpuimagesize=sizeof(char[3])*res,
			gpuindicessize=sizeof(int[4]),
			gpuplanessize=sizeof(short[3])*res,
			gpuhistsize=sizeof(int[3*PRED_COUNT*256])*nblocks,
			gpucsizessize=sizeof(float[3*PRED_COUNT])*nblocks;
		cl_mem gpu_image	=clCreateBuffer(ctx->context, CL_MEM_READ_ONLY, gpuimagesize, 0, &error);	CHECKCL(error);
		cl_mem gpu_indices	=clCreateBuffer(ctx->context, CL_MEM_READ_ONLY, gpuindicessize, 0, &error);	CHECKCL(error);
		cl_mem gpu_planes	=clCreateBuffer(ctx->context, CL_MEM_READ_WRITE, gpuplanessize, 0, &error);	CHECKCL(error);
		cl_mem gpu_hist		=clCreateBuffer(ctx->context, CL_MEM_READ_WRITE, gpuhistsize, 0, &error);	CHECKCL(error);
		cl_mem gpu_csizes	=clCreateBuffer(ctx->context, CL_MEM_READ_WRITE, gpucsizessize, 0, &error);	CHECKCL(error);

		size_t gpumemusage=
			+gpuimagesize
			+gpuindicessize
			+gpuplanessize
			+gpuhistsize
			+gpucsizessize
			;
		if(loud)
			printf("GPU mem usage: %lf MB\n", (double)gpumemusage/(1024*1024));

		const int RCTcoeffs[]=
		{
			0, 0,
			4, 0,
			3, 1,
			2, 2,
			1, 3,
			0, 4,
		};
		clEnqueueWriteBuffer(ctx->commandqueue, gpu_image, CL_FALSE, 0, gpuimagesize, image, 0, 0, 0);
		float *csizes=(float*)malloc(gpucsizessize*6);
		if(!csizes)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		int csizestride=nblocks*(3 * PRED_COUNT);
		error=clSetKernelArg(ctx->kernels[OCL_prep_planes], 0, sizeof(cl_mem), &gpu_indices);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_prep_planes], 1, sizeof(cl_mem), &gpu_image);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_prep_planes], 2, sizeof(cl_mem), &gpu_planes);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_pred_planes], 0, sizeof(cl_mem), &gpu_indices);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_pred_planes], 1, sizeof(cl_mem), &gpu_planes);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_pred_planes], 2, sizeof(cl_mem), &gpu_hist);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_calc_entropy], 0, sizeof(cl_mem), &gpu_indices);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_calc_entropy], 1, sizeof(cl_mem), &gpu_hist);	CHECKCL(error);
		error=clSetKernelArg(ctx->kernels[OCL_calc_entropy], 2, sizeof(cl_mem), &gpu_csizes);	CHECKCL(error);
		for(int k=0;k<6;++k)
		{
			int indices[]=
			{
				iw, ih,
				RCTcoeffs[k<<1|0], RCTcoeffs[k<<1|1],
			};
			size_t worksize=res;
			error=clEnqueueWriteBuffer(ctx->commandqueue, gpu_indices, CL_FALSE, 0, sizeof(indices), indices, 0, 0, 0);	CHECKCL(error);
			error=clFlush(ctx->commandqueue);	CHECKCL(error);
			error=clFinish(ctx->commandqueue);	CHECKCL(error);
			error=clEnqueueNDRangeKernel(ctx->commandqueue, ctx->kernels[OCL_prep_planes], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
			error=clFlush(ctx->commandqueue);	CHECKCL(error);
			error=clFinish(ctx->commandqueue);	CHECKCL(error);

			worksize=3LL*nblocks;
			indices[2]=xblocks;
			indices[3]=yblocks;
			error=clEnqueueWriteBuffer(ctx->commandqueue, gpu_indices, CL_TRUE, 0, sizeof(indices), indices, 0, 0, 0);	CHECKCL(error);
			error=clEnqueueNDRangeKernel(ctx->commandqueue, ctx->kernels[OCL_pred_planes], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
			error=clFlush(ctx->commandqueue);	CHECKCL(error);
			error=clFinish(ctx->commandqueue);	CHECKCL(error);

			worksize=3LL*PRED_COUNT*nblocks;
			indices[0]=(int)(worksize*256*k);
			error=clEnqueueWriteBuffer(ctx->commandqueue, gpu_indices, CL_TRUE, 0, sizeof(indices), indices, 0, 0, 0);	CHECKCL(error);
			error=clEnqueueNDRangeKernel(ctx->commandqueue, ctx->kernels[OCL_calc_entropy], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
			error=clFlush(ctx->commandqueue);	CHECKCL(error);
			error=clFinish(ctx->commandqueue);	CHECKCL(error);

			error=clEnqueueReadBuffer(ctx->commandqueue, gpu_csizes, CL_TRUE, 0, gpucsizessize, csizes+csizestride*k, 0, 0, 0);	CHECKCL(error);
		}
		error=clReleaseMemObject(gpu_image);	CHECKCL(error);
		error=clReleaseMemObject(gpu_indices);	CHECKCL(error);
		error=clReleaseMemObject(gpu_planes);	CHECKCL(error);
		error=clReleaseMemObject(gpu_hist);	CHECKCL(error);
		error=clReleaseMemObject(gpu_csizes);	CHECKCL(error);
		if(loud)
			prof_check(CHECK_ANALYSIS);
		
		//Step 2: Select decorrelation method
		int *RCTinfo=(int*)malloc(nblocks*sizeof(int));
		if(!RCTinfo)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		for(int kb=0;kb<nblocks;++kb)
		{
			int offset=kb*(3*PRED_COUNT);
			float *curr_csizes[]=
			{
				csizes+csizestride*0+offset,//R G B
				csizes+csizestride*1+offset,//RG GB BR
				csizes+csizestride*2+offset,//R3 G3 B3
				csizes+csizestride*3+offset,//R2 G2 B2
				csizes+csizestride*4+offset,//R1 G1 B1
				csizes+csizestride*5+offset,//RB GR BG
			};
			int predsel[3*6]={0};
			for(int kg=0;kg<6;++kg)
			{
				for(int kc=0;kc<3;++kc)
				{
					float *c3_csizes=curr_csizes[kg]+PRED_COUNT*kc;
					int bestpred=0;
					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						if(c3_csizes[bestpred]>c3_csizes[kp])
							bestpred=kp;
					}
					predsel[kg*3+kc]=bestpred;
				}
			}
			int bestrct=0;
			float bestsize=0;
			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *group=rct_combinations[kt];
				float csize=
					+curr_csizes[group[0]%6][group[0]/6*PRED_COUNT+predsel[group[0]]]
					+curr_csizes[group[1]%6][group[1]/6*PRED_COUNT+predsel[group[1]]]
					+curr_csizes[group[2]%6][group[2]/6*PRED_COUNT+predsel[group[2]]]
					;
				if(!kt||bestsize>csize)
				{
					bestsize=csize;
					bestrct=kt;
				}
			}
			{
				const unsigned char *group=rct_combinations[bestrct];
				RCTinfo[kb]=
					predsel[group[0]]<<24|
					predsel[group[1]]<<16|
					predsel[group[0]]<< 8|
					bestrct << 0;//18*18*18*43 = 0x3D398	18/32 bit

#ifdef _DEBUG
				if(loud)
					printf("%-7s %-3s %-3s %-3s%c",
						rct_names[bestrct],
						pred_names[predsel[group[0]]],
						pred_names[predsel[group[1]]],
						pred_names[predsel[group[2]]],
						(kb+1)%xblocks?'\t':'\n'
					);
#endif
			}
		}
		//int predselsize=nblocks*sizeof(int[3*6]);
		//int *predsel=(int*)malloc(predselsize);
		//if(!predsel)
		//{
		//	LOG_ERROR("Alloc error");
		//	return 1;
		//}
		//for(int kc=0;kc<nblocks*3*6;++kc)//select best predictors
		//{
		//	int bestpred=0;
		//	float *curr_csizes=csizes+NPREDS*kc;
		//	for(int kp=1;kp<NPREDS;++kp)
		//	{
		//		if(curr_csizes[bestpred]>curr_csizes[kp])
		//			bestpred=kp;
		//	}
		//	predsel[kc]=bestpred;
		//}
		free(csizes);
		if(loud)
			prof_check(CHECK_SELECT);
		
		//Step 3: Decorrelate
		int psize=(iw+16LL)*sizeof(short[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
		short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
		unsigned char *planes=(unsigned char*)malloc(sizeof(char[3])*res);
		if(!pixels)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(pixels, 0, psize);
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			int by=ky/BLOCKSIZE;
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			int preds[3]={0};
			for(int kx=0;kx<iw;++kx, idx+=3, ++idx2)
			{
				int bx=kx/BLOCKSIZE;
				unsigned char *decorrflags=(unsigned char*)(RCTinfo+xblocks*by+bx);
				const unsigned char *group=rct_combinations[decorrflags[0]];

				int yuv[]=
				{
					image[idx+group[3+0]]-128,
					image[idx+group[3+1]]-128,
					image[idx+group[3+2]]-128,
				};
				for(int kc=0;kc<3;++kc)
				{
					preds[kc]=predict(decorrflags[kc+1], rows);
					rows[0]+=2;
					rows[1]+=2;
					rows[2]+=2;
					rows[3]+=2;
				}
				rows[0]+=2;
				rows[1]+=2;
				rows[2]+=2;
				rows[3]+=2;
				
				rows[0][-1*4*2+0]=yuv[0];
				int offset=yuv[0]&(char)group[6];
				preds[1]+=offset;
				rows[0][-1*4*2+1]=yuv[1]-offset;
				CLAMP2(preds[1], -128, 127);
				offset=(group[7]*yuv[0]+group[8]*yuv[1])>>2;
				preds[2]+=offset;
				rows[0][-1*4*2+2]=yuv[2]-offset;
				CLAMP2(preds[2], -128, 127);

				planes[res*0+idx2]=(unsigned char)(yuv[0]-preds[0]+128);
				planes[res*1+idx2]=(unsigned char)(yuv[1]-preds[1]+128);
				planes[res*2+idx2]=(unsigned char)(yuv[2]-preds[2]+128);
			}
		}
		_mm_free(pixels);
		if(loud)
			prof_check(CHECK_PREDICT);

		//Step 4: Entropy coding
		unsigned char *cdata[3]={0};
		unsigned cdatasize[3]={0};
		cdata[0]=rans_compress_O1_32x16_avx2(planes+res*0, (unsigned)res, cdata[0], cdatasize+0);
		cdata[1]=rans_compress_O1_32x16_avx2(planes+res*1, (unsigned)res, cdata[1], cdatasize+1);
		cdata[2]=rans_compress_O1_32x16_avx2(planes+res*2, (unsigned)res, cdata[2], cdatasize+2);
		free(planes);
		if(loud)
			prof_check(CHECK_ENTROPY_CODING);

		//Step 5: Write file
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			filesize+=fwrite(&iw, 1, 4, fdst);
			filesize+=fwrite(&ih, 1, 4, fdst);
			filesize+=fwrite(cdatasize+0, 1, 4, fdst);
			filesize+=fwrite(cdatasize+1, 1, 4, fdst);
			filesize+=fwrite(RCTinfo, 1, nblocks*sizeof(int), fdst);
			filesize+=fwrite(cdata[0], 1, cdatasize[0], fdst);
			filesize+=fwrite(cdata[1], 1, cdatasize[1], fdst);
			filesize+=fwrite(cdata[2], 1, cdatasize[2], fdst);
			fclose(fdst);
		}
		free(RCTinfo);
		free(cdata[0]);
		free(cdata[1]);
		free(cdata[2]);
		if(loud)
			prof_check(CHECK_WRITE);
	}
	else//decode
	{
		int cdatasizes[3]={0};
		memcpy(cdatasizes+0, srcptr, 4); srcptr+=4;
		memcpy(cdatasizes+1, srcptr, 4); srcptr+=4;
		cdatasizes[2]=(int)(srcbuf+srcsize-(srcptr+cdatasizes[0]+cdatasizes[1]));
		
		int *RCTinfo=(int*)malloc(nblocks*sizeof(int));
		unsigned char *planes=(unsigned char*)malloc(sizeof(char[3])*res);
		unsigned char *image=(unsigned char*)malloc(sizeof(char[3])*res);
		int psize=(iw+16LL)*sizeof(short[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
		short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
		if(!RCTinfo||!planes||!image||!pixels)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memcpy(RCTinfo, srcptr, nblocks*sizeof(int)); srcptr+=nblocks*sizeof(int);
		memset(pixels, 0, psize);
		
		unsigned char *ret[3]={0};
		ret[0]=rans_uncompress_O1_32x16_avx2(srcptr, cdatasizes[0], planes+res*0, (unsigned)res); srcptr+=cdatasizes[0];
		ret[1]=rans_uncompress_O1_32x16_avx2(srcptr, cdatasizes[1], planes+res*1, (unsigned)res); srcptr+=cdatasizes[1];
		ret[2]=rans_uncompress_O1_32x16_avx2(srcptr, cdatasizes[2], planes+res*2, (unsigned)res); //srcptr+=cdatasizes[2];
		if(loud)
			prof_check(CHECK_ENTROPY_CODING);
		if(!ret[0]||!ret[1]||!ret[2])
		{
			LOG_ERROR("HTS Decode error");
			return 1;
		}
		
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			int by=ky/BLOCKSIZE;
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			int preds[3]={0}, yuv[3]={0};
			int kx, bx;
			unsigned char *decorrflags;
			const unsigned char *group;
			for(kx=0;kx<iw;++kx, idx+=3, ++idx2)
			{
				bx=kx/BLOCKSIZE;
				decorrflags=(unsigned char*)(RCTinfo+xblocks*by+bx);
				group=rct_combinations[decorrflags[0]];

				yuv[0]=planes[res*0+idx2]-128;
				yuv[1]=planes[res*1+idx2]-128;
				yuv[2]=planes[res*2+idx2]-128;
				for(int kc=0;kc<3;++kc)
				{
					preds[kc]=predict(decorrflags[kc+1], rows);
					rows[0]+=2;
					rows[1]+=2;
					rows[2]+=2;
					rows[3]+=2;
				}
				rows[0]+=2;
				rows[1]+=2;
				rows[2]+=2;
				rows[3]+=2;

				rows[0][-1*4*2+0]=yuv[0]=(yuv[0]+preds[0])<<24>>24;
				int offset=yuv[0]&(char)group[6];
				preds[1]+=offset;
				CLAMP2(preds[1], -128, 127);
				yuv[1]=(yuv[1]+preds[1])<<24>>24;
				rows[0][-1*4*2+1]=yuv[1]-offset;
				offset=(group[7]*yuv[0]+group[8]*yuv[1])>>2;
				preds[2]+=offset;
				CLAMP2(preds[2], -128, 127);
				yuv[2]=(yuv[2]+preds[2])<<24>>24;
				rows[0][-1*4*2+2]=yuv[2]-offset;
				
				image[idx+group[3+0]]=(unsigned char)(yuv[0]+128);
				image[idx+group[3+1]]=(unsigned char)(yuv[1]+128);
				image[idx+group[3+2]]=(unsigned char)(yuv[2]+128);

				guide_check(image, kx, ky);
			}
		}
		_mm_free(pixels);
		free(planes);
		free(RCTinfo);
		if(loud)
			prof_check(CHECK_UNPREDICT);
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, 3LL*res, fdst);
			fclose(fdst);
		}
		free(image);
		if(loud)
			prof_check(CHECK_WRITE);
	}
	if(loud)
	{
		ptrdiff_t usize=res*3;
		celapsed=__rdtsc()-celapsed;
		elapsed=time_sec()-elapsed;
		if(fwd)
			printf("%12zd/%12td  %lf:1\n", filesize, usize, (double)usize/filesize);
		printf("%c %12.6lf sec  %lf MB/s  %lf CPB\n",
			'D'+fwd,
			elapsed,
			usize/(elapsed*1024*1024),
			(double)celapsed/usize
		);
		prof_print();
		printf("\n");
	}
	return 0;
}