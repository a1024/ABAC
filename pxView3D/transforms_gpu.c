#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
#ifdef ALLOW_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include<CL/cl.h>
#pragma comment(lib, "OpenCL.lib")
#endif
static const char file[]=__FILE__;

#ifdef ALLOW_OPENCL
const char*		clerr2str(int error)
{
	const char *a=0;
#define 		EC(X)		case X:a=(const char*)#X;break;
#define 		EC2(N, X)	case N:a=(const char*)#X;break;
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
#define CHECKCL(E) (!(E)||LOG_ERROR("CL Error %s", clerr2str(E)))

static int clctxinitialized=0;
static cl_platform_id platform=0;
static cl_device_id device=0;
static cl_context context=0;
static cl_command_queue commandqueue=0;
static cl_program program=0;

#define CLKERNELNALELIST CLKERNEL(train)
typedef enum		CLKernelIdxEnum
{
#define				CLKERNEL(LABEL)	OCL_##LABEL,
CLKERNELNALELIST
#undef				CLKERNEL
	OCL_NKERNELS,
} CLKernelIdx;
const char			*kernelnames[]=
{
#define				CLKERNEL(LABEL)	#LABEL,
CLKERNELNALELIST
#undef				CLKERNEL
};
static cl_kernel kernels[OCL_NKERNELS]={0};
int init_ocl()
{
	int error=0;
	unsigned count=0;
	if(clctxinitialized)
		return 1;
	clctxinitialized=1;

	//get platform
	error=clGetPlatformIDs(0, 0, &count);	CHECKCL(error);
	if(!count)
	{
		LOG_ERROR("No OpenCL platforms");
		return 0;
	}
	error=clGetPlatformIDs(1, &platform, 0);	CHECKCL(error);

	//get device
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, 0, &count);	CHECKCL(error);
	if(!count)
	{
		LOG_ERROR("No OpenCL devices");
		return 0;
	}
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, 0);	CHECKCL(error);
	
	//create context
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
	context=clCreateContext(0, 1, &device, 0, 0, &error);	CHECKCL(error);
#endif
	
	//create command queue
#if CL_TARGET_OPENCL_VERSION>=200
	//cl_queue_properties properties[]=
	//{
	//	CL_QUEUE_PROPERTIES, 0,
	//	CL_QUEUE_SIZE, 0,
	//	0,
	//};
	commandqueue=clCreateCommandQueueWithProperties(context, device, 0, &error);	CHECKCL(error);
#else
	commandqueue=clCreateCommandQueue(context, device, 0, &error);	CHECKCL(error);
#endif

	//build kernels
	{
		ArrayHandle srctext=load_file("E:/C/pxView3D/pxView3D/cl_kernels.h", 0, 0);
		const char *k_src=(const char*)srctext->data;
		size_t k_len=srctext->count;
		
		program=clCreateProgramWithSource(context, 1, (const char**)&k_src, &k_len, &error);	CHECKCL(error);
		error=clBuildProgram(program, 1, &device, g_buf, 0, 0);
		if(error)
		{
			size_t retlen=0;
			error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, G_BUF_SIZE, g_buf, &retlen);
			if(retlen>G_BUF_SIZE)
			{
				char *buf=(char*)malloc(retlen+10);
				error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, retlen+10, buf, &retlen);	CHECKCL(error);
				messagebox(MBOX_OK, "OpenCL compilation failed", "%s", buf);
				//printf("\nOpenCL compilation failed:\n%s\n", buf);
				free(buf);
				LOG_ERROR("Aborting");
			}
			else
				LOG_ERROR("OpenCL Compilation failed:\n%s\n", g_buf);
		}
		array_free(&srctext);
	}

	//fetch entry points
	for(int k=0;k<OCL_NKERNELS;++k)
	{
		kernels[k]=clCreateKernel(program, kernelnames[k], &error);
		if(error)
		{
			LOG_ERROR("Couldn't find kernel %s", kernelnames[k]);
			return 0;
		}
	}
	return 1;
}


void pred_learned_gpu(char *buf, int iw, int ih, int fwd)
{
}
#endif