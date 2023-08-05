#include"e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef ALLOW_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include<CL/cl.h>
#pragma comment(lib, "OpenCL.lib")
#endif
static const char file[]=__FILE__;

#ifdef ALLOW_OPENCL
#define KERNEL_FN "E:/C/e2/e2/kernels.h"

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
#define CL_CHECK(ERR) ASSERT_MSG(!(ERR), "OpenCL error %d: %s\n", ERR, clerr2str(ERR))

unsigned nplatforms=0, ndevices=0;
cl_platform_id *platforms=0;
cl_device_id *devices=0;
size_t maxlocalsize=0;
cl_context context=0;
cl_command_queue commandqueue=0;
cl_program program=0;

int init_cl(const char *searchpath)
{
	int error;

	//oclw_loadAPI();

	error=clGetPlatformIDs(0, 0, &nplatforms);		CL_CHECK(error);
	ASSERT_MSG(nplatforms, "No OpenCL platforms");
	
	platforms=(cl_platform_id*)malloc(nplatforms*sizeof(cl_platform_id));
	error=clGetPlatformIDs(nplatforms, platforms, 0);	CL_CHECK(error);

	ndevices=0;
	int pl_idx=0;
	size_t temp=0, retlen=0;
	for(;pl_idx<(int)nplatforms;++pl_idx)
	{
		error=clGetPlatformInfo(platforms[pl_idx], CL_PLATFORM_VERSION, G_BUF_SIZE, g_buf, &retlen);CL_CHECK(error);
		error=clGetDeviceIDs(platforms[pl_idx], CL_DEVICE_TYPE_GPU, 0, 0, &ndevices);
		printf("OpenCL platform: %s, %d device(s)\n", g_buf, ndevices);
		if(ndevices)//break before pl_idx increment
			break;
	}
	CL_CHECK(error);
	ASSERT_MSG(ndevices, "No OpenCL devices");

	devices=(cl_device_id*)malloc(ndevices*sizeof(cl_device_id));
	error=clGetDeviceIDs(platforms[pl_idx], CL_DEVICE_TYPE_GPU, ndevices, devices, 0);	CL_CHECK(error);
	
	//get info
	error=clGetDeviceInfo(devices[0], CL_DEVICE_NAME, G_BUF_SIZE, g_buf, &retlen);	CL_CHECK(error);
	printf("Device: %s\n", g_buf);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_VENDOR, G_BUF_SIZE, g_buf, &retlen);	CL_CHECK(error);
	printf("Vendor: %s\n", g_buf);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &maxlocalsize, &retlen);	CL_CHECK(error);
	printf("CL_DEVICE_MAX_WORK_GROUP_SIZE = %d\n", (int)maxlocalsize);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF, sizeof(size_t), &temp, &retlen);	CL_CHECK(error);
	printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF = %d\n", (int)temp);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(size_t), &temp, &retlen);	CL_CHECK(error);
	printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT = %d\n", (int)temp);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_ADDRESS_BITS, sizeof(size_t), &temp, &retlen);	CL_CHECK(error);
	printf("CL_DEVICE_ADDRESS_BITS = %d\n", (int)temp);
	error=clGetDeviceInfo(devices[0], CL_DEVICE_EXTENSIONS, G_BUF_SIZE, g_buf, &retlen);	CL_CHECK(error);
	for(int k=0;k<retlen;++k)
		if(g_buf[k]==' ')
			g_buf[k]='\n';
	printf("Extensions:\n%s\n", g_buf);
	//int cl_gl_interop=strstr(g_buf, "cl_khr_gl_sharing")!=0;
	//if(!cl_gl_interop)
	//	printf("\n\tcl_khr_gl_sharing not supported\n\n");
	
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
	context=clCreateContext(properties, 1, &devices[0], nullptr, nullptr, &error);	CL_CHECK(error);
#else
	context=clCreateContext(0, ndevices, devices, 0, 0, &error);			CL_CHECK(error);
#endif
#if CL_TARGET_OPENCL_VERSION>=200
	//cl_queue_properties properties[]=
	//{
	//	CL_QUEUE_PROPERTIES, 0,
	//	CL_QUEUE_SIZE, 0,
	//	0,
	//};
	commandqueue=clCreateCommandQueueWithProperties(context, devices[0], 0, &error);	CL_CHECK(error);
#else
	commandqueue=clCreateCommandQueue(context, devices[0], 0, &error);	CL_CHECK(error);
#endif

	//compile kernel
	ArrayHandle filename=searchfor_file(searchpath, "kernels.h");
	if(!filename)
	{
		LOG_ERROR("Cannot open kernels.h");
		return 0;
	}
	ArrayHandle srctext=load_file(filename->data, 0, 0);
	if(!srctext)
	{
		LOG_ERROR("Cannot open kernels.h");
		return 0;
	}
	array_free(&filename);
	//ASSERT_MSG(srctext, "Cannot open \'%s\'", srcname);
	const char *k_src=(const char*)srctext->data;
	size_t k_len=srctext->count;
	program=clCreateProgramWithSource(context, 1, (const char**)&k_src, &k_len, &error);	CL_CHECK(error);
#ifdef PREC_FIXED
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__ -DPREC_FIXED=%d", PREC_FIXED);
#elif defined PREC_HALF
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__ -DPREC_HALF");
#else
	sprintf_s(g_buf, G_BUF_SIZE, "-D__OPEN_CL__");
#endif
	error=clBuildProgram(program, 1, devices, g_buf, 0, 0);
	if(error)
	{
		error=clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, G_BUF_SIZE, g_buf, &retlen);
		if(retlen>G_BUF_SIZE)
		{
			char *buf=(char*)malloc(retlen+10);
			error=clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, retlen+10, buf, &retlen);	CL_CHECK(error);
			printf("\nOpenCL Compilation failed:\n%s\n", buf);
			free(buf);
			LOG_ERROR("Aborting");
		}
		else
			LOG_ERROR("\nOpenCL Compilation failed:\n%s\n", g_buf);
		return 0;
	}
	array_free(&srctext);

#if 0
	for(int k=0;k<OCL_NKERNELS;++k)
	{
		kernels[k]=clCreateKernel(program, kernelnames[k], &error);
		ASSERT_MSG(!error, "Couldn't find kernel %s", kernelnames[k]);
	}
#endif
	return 1;
}
#endif