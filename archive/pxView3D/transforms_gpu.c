#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
#ifdef ALLOW_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include<CL/cl.h>
#pragma comment(lib, "OpenCL.lib")
#endif
static const char file[]=__FILE__;

#ifdef ENABLE_CONSOLE_MAIN_TEST
#include<stdio.h>
#define set_window_title printf
#endif

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


#define CLKERNELNALELIST CLKERNEL(custom3_eval)
//#define CLKERNELNALELIST CLKERNEL(train) CLKERNEL(update) CLKERNEL(predict)

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
static int ctxinitialized=0;
static cl_platform_id platform=0;
static cl_device_id device=0;
static cl_context context=0;
static cl_command_queue commandqueue=0;
static cl_program program=0;
static cl_kernel kernels[OCL_NKERNELS]={0};
int init_ocl()
{
	int error=0;
	unsigned count=0;
	if(ctxinitialized)
		return 1;
	ctxinitialized=1;

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
		ArrayHandle srctext=load_file("cl_kernels2.h", 0, 0, 0);
		if(!srctext)
			srctext=load_file("C:/Projects/pxView3D/pxView3D/cl_kernels2.h", 0, 0, 0);
			//srctext=load_file("E:/C/pxView3D/pxView3D/cl_kernels.h", 0, 0, 0);
		if(!srctext)
		{
			LOG_ERROR("Cannot open cl_kernels2.h");
			return 0;
		}
		const char *k_src=(const char*)srctext->data;
		size_t k_len=srctext->count;
		
		program=clCreateProgramWithSource(context, 1, (const char**)&k_src, &k_len, &error);	CHECKCL(error);
		error=clBuildProgram(program, 1, &device, 0, 0, 0);
		//error=clBuildProgram(program, 1, &device, g_buf, 0, 0);
		if(error)
		{
			size_t retlen=0;
			error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, 0, &retlen);
			//if(retlen>G_BUF_SIZE)
			//{
				char *buf=(char*)malloc(retlen+10);
				error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, retlen+10, buf, &retlen);	CHECKCL(error);
				copy_to_clipboard(buf, (int)retlen);
				messagebox(MBOX_OK, "OpenCL compilation failed", "Errors copied to clipboard");
				//printf("\nOpenCL compilation failed:\n%s\n", buf);
				free(buf);
				return 0;
				//LOG_ERROR("Aborting");
			//}
			//else
			//	LOG_ERROR("OpenCL Compilation failed:\n%s\n", g_buf);
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


#if 0
//SYNC THIS WITH CL_KERNELS.H:

//	#define ZIGZAG_TRAVERSAL

#define V1_XBLOCKS 1//leave XBLOCKS at 1 for optimal performance
#define V1_YBLOCKS 8
#define V1_REACH 7
#define V1_NITER 7
#define V1_XYWRITERS (V1_XBLOCKS*V1_YBLOCKS)
#define V1_TOTALWRITERS (3*V1_XBLOCKS*V1_YBLOCKS)
#define V1_KERNELSIZE (V1_REACH*(V1_REACH+1)*2)
#define V1_NLANES (V1_TOTALWRITERS*V1_KERNELSIZE)

#define NF0_MAIN	55		//55
#define NF0_ERRORS	41		//41
#define NF0_PX		0		//41
#define NF0 (NF0_MAIN+NF0_ERRORS+NF0_PX)
#define NF1 16
#define NF2 16
#define NF3 96	//96
typedef struct TempsStruct
{
	float
		nb[NF0],
		net1[NF1], x1[NF1],
		net2[NF2], x2[NF2],
		net3[NF3], x3[NF3],
		pred, expr, loss;
} Temps;
typedef struct BwdTempsStruct
{
	float dL_dp, dL_dx3[NF3], dL_dx2[NF2], dL_dx1[NF1];
} BwdTemps;
typedef struct ParamsStruct
{
	float
		weight1[NF1*NF0], bias1[NF1],
		weight2[NF2*NF1], bias2[NF2],
		weight3[NF3*NF2], bias3[NF3],
		weight4[1*NF3];
} Params;
static const int NPARAMS=sizeof(Params)/sizeof(float);


static void initialize(float *w, int count, int fan_in)
{
	float gain=1/(0x10000*sqrtf((float)fan_in));
	for(int k=0;k<count;++k)
	{
		int x=(int)(xoroshiro128_next()&0x1FFFF)-0x10000;//[-2^16, 2^16]
		float val=(float)x*gain;//[-2^-4, 2^-4]/sqrt_fan_in
		w[k]=val;
	}
}
static void initialize_all(Params *p, int reset_prng)
{
	if(reset_prng)
		XOROSHIRO128_RESET();
	initialize(p->weight1, _countof(p->weight1), NF0);
	initialize(p->bias1  , _countof(p->bias1  ), NF0);
	initialize(p->weight2, _countof(p->weight2), NF1);
	initialize(p->bias2  , _countof(p->bias2  ), NF1);
	initialize(p->weight3, _countof(p->weight3), NF2);
	initialize(p->bias3  , _countof(p->bias3  ), NF2);
	initialize(p->weight4, _countof(p->weight4), NF3);
}

void pred_learned_gpu(char *buf, int iw, int ih, int fwd)
{
	double t_start=time_ms();
	const int nlanes=V1_TOTALWRITERS*(V1_REACH*((V1_REACH+1)<<1));//
	int res=iw*ih;
	if(!init_ocl())
		return;
	
	//int pretrained=1;//

	int error=0;
	float *buf2=(float*)malloc(res*3LL*sizeof(float));
	cl_mem gpu_pixels=clCreateBuffer(context, CL_MEM_READ_WRITE, res*3LL*sizeof(float), 0, &error);	CHECKCL(error);
	cl_mem gpu_errors=clCreateBuffer(context, CL_MEM_READ_WRITE, res*3LL*sizeof(float), 0, &error);	CHECKCL(error);
	//cl_mem gpu_pixels=clCreateBuffer(context, CL_MEM_READ_WRITE, res*3LL*sizeof(float), fwd?buf2:0, &error);	CHECKCL(error);
	//cl_mem gpu_errors=clCreateBuffer(context, CL_MEM_READ_WRITE, res*3LL*sizeof(float), fwd?0:buf2, &error);	CHECKCL(error);

	Params *params=(Params*)malloc(sizeof(Params)*V1_TOTALWRITERS);
	//Params *params=(Params*)malloc(sizeof(Params));
	cl_mem gpu_params   =clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Params)*V1_TOTALWRITERS, 0, &error);	CHECKCL(error);
	cl_mem gpu_temps    =clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Temps   )*nlanes, 0, &error);		CHECKCL(error);
	cl_mem gpu_bwdtemps =clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(BwdTemps)*nlanes, 0, &error);		CHECKCL(error);
	cl_mem gpu_gradients=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Params  )*nlanes, 0, &error);		CHECKCL(error);
	
	//float *features=(float*)malloc(NF0*nlanes*sizeof(float));
	cl_mem gpu_features=clCreateBuffer(context, CL_MEM_READ_WRITE, NF0*nlanes*sizeof(float), 0, &error);	CHECKCL(error);
	
	int indices[]={0, 0, iw, ih, fwd, 0};
	cl_mem gpu_indices=clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(indices), 0, &error);	CHECKCL(error);
	
	if(!buf2||!gpu_pixels||!gpu_errors||!params||!gpu_params||!gpu_temps||!gpu_bwdtemps||!gpu_gradients||!gpu_features||!gpu_indices)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	initialize_all(params, 1);
	//for(int k=1;k<V1_TOTALWRITERS;++k)
	//	initialize_all(params+k, 0);
	memfill(params+1, params, (V1_TOTALWRITERS-1)*sizeof(*params), sizeof(*params));
	float gain=1.f/128;
	for(int k=0;k<res;++k)
	{
		int idx=k*3;
		buf2[res*0+k]=(buf[k<<2|0]+0.5f)*gain;
		buf2[res*1+k]=(buf[k<<2|1]+0.5f)*gain;
		buf2[res*2+k]=(buf[k<<2|2]+0.5f)*gain;
	}

	error=clEnqueueWriteBuffer(commandqueue, fwd?gpu_pixels:gpu_errors, CL_FALSE, 0, res*3LL*sizeof(float), buf2, 0, 0, 0);		CHECKCL(error);
	//error=clEnqueueFillBuffer(commandqueue, gpu_params, params, sizeof(Params), 0, sizeof(Params)*V1_TOTALWRITERS, 0, 0, 0);	CHECKCL(error);//X  pattern size must be from {1, ...128} bytes
	error=clEnqueueWriteBuffer(commandqueue, gpu_params, CL_FALSE, 0, sizeof(Params)*V1_TOTALWRITERS, params, 0, 0, 0);		CHECKCL(error);

	error=clSetKernelArg(kernels[OCL_train], 0, sizeof(cl_mem), &gpu_pixels);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 1, sizeof(cl_mem), &gpu_errors);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 2, sizeof(cl_mem), &gpu_indices);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 3, sizeof(cl_mem), &gpu_params);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 4, sizeof(cl_mem), &gpu_temps);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 5, sizeof(cl_mem), &gpu_bwdtemps);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_train], 6, sizeof(cl_mem), &gpu_gradients);CHECKCL(error);

	error=clSetKernelArg(kernels[OCL_update], 0, sizeof(cl_mem), &gpu_params);		CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_update], 1, sizeof(cl_mem), &gpu_gradients);	CHECKCL(error);

	error=clSetKernelArg(kernels[OCL_predict], 0, sizeof(cl_mem), &gpu_pixels);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_predict], 1, sizeof(cl_mem), &gpu_errors);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_predict], 2, sizeof(cl_mem), &gpu_indices);CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_predict], 3, sizeof(cl_mem), &gpu_params);	CHECKCL(error);
	error=clSetKernelArg(kernels[OCL_predict], 4, sizeof(cl_mem), &gpu_temps);	CHECKCL(error);
	
	const int blockw=(iw+V1_XBLOCKS-1)/V1_XBLOCKS, blockh=(ih+V1_YBLOCKS-1)/V1_YBLOCKS;

	//for(int kt=0;kt<1+pretrained;++kt)//pretrained
	for(int ky=0;ky<blockh;++ky)
	{
		//error=clEnqueueWriteBuffer(commandqueue, gpu_params, TRUE, 0, sizeof(Params)*V1_TOTALWRITERS, params, 0, 0, 0);			CHECKCL(error);//X
		{
			TimeInfo ti;
			parsetimedelta(time_ms()-t_start, &ti);
			set_window_title("%d/%d, %.2lf%% - %02d-%02d-%06.3f"
#ifdef ENABLE_CONSOLE_MAIN_TEST
				"\n"
#endif
				, ky, blockh, 100.*ky/blockh, ti.hours, ti.mins, ti.secs);
		}
#ifdef ZIGZAG_TRAVERSAL
		int xstart, xend, xstep;
		indices[5]=ky&1;//flip east & west
		if(indices[5])
			xstart=iw-1, xend=0, xstep=-1;
		else
			xstart=0, xend=iw-1, xstep=1;
		for(int kx=xstart;kx!=xend;kx+=xstep)
#else
		for(int kx=0;kx<blockw;++kx)
#endif
		{
			size_t worksize;

			indices[0]=kx;
			indices[1]=ky;
			error=clEnqueueWriteBuffer(commandqueue, gpu_indices, CL_FALSE, 0, sizeof(indices), indices, 0, 0, 0);	CHECKCL(error);
			
			error=clEnqueueBarrierWithWaitList(commandqueue, 0, 0, 0);	CHECKCL(error);
#if 1
			int niter=V1_NITER+(kx?10/kx:5);
			for(int k=0;k<niter;++k)
			{
				worksize=V1_NLANES;
				error=clEnqueueNDRangeKernel(commandqueue, kernels[OCL_train], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
				error=clEnqueueBarrierWithWaitList(commandqueue, 0, 0, 0);	CHECKCL(error);
				//error=clFlush(commandqueue);	CHECKCL(error);
				//error=clFinish(commandqueue);	CHECKCL(error);

				worksize=V1_TOTALWRITERS*NPARAMS;
				error=clEnqueueNDRangeKernel(commandqueue, kernels[OCL_update], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
				error=clEnqueueBarrierWithWaitList(commandqueue, 0, 0, 0);	CHECKCL(error);
				//error=clFlush(commandqueue);	CHECKCL(error);
				//error=clFinish(commandqueue);	CHECKCL(error);
			}
#endif

			worksize=V1_TOTALWRITERS;
			error=clEnqueueNDRangeKernel(commandqueue, kernels[OCL_predict], 1, 0, &worksize, 0, 0, 0, 0);	CHECKCL(error);
			//error=clEnqueueBarrierWithWaitList(commandqueue, 0, 0, 0);	CHECKCL(error);
			error=clFlush(commandqueue);	CHECKCL(error);
			error=clFinish(commandqueue);	CHECKCL(error);
		}
	}

	error=clEnqueueReadBuffer(commandqueue, fwd?gpu_errors:gpu_pixels, CL_TRUE, 0, res*3LL*sizeof(float), buf2, 0, 0, 0);

	for(int k=0;k<res;++k)
	{
		int idx=k*3;
		float val;
		val=floorf(buf2[res*0+k]*128), buf[k<<2|0]=(char)CLAMP(-128, val, 127);
		val=floorf(buf2[res*1+k]*128), buf[k<<2|1]=(char)CLAMP(-128, val, 127);
		val=floorf(buf2[res*2+k]*128), buf[k<<2|2]=(char)CLAMP(-128, val, 127);
	}
	{
		TimeInfo ti;
		parsetimedelta(time_ms()-t_start, &ti);
		set_window_title("Done %02d-%02d-%06.3f"
#ifdef ENABLE_CONSOLE_MAIN_TEST
			"\n"
#endif
			, ti.hours, ti.mins, ti.secs);
	}
	clReleaseMemObject(gpu_indices);
	clReleaseMemObject(gpu_features);
	//free(features);
	clReleaseMemObject(gpu_gradients);
	clReleaseMemObject(gpu_bwdtemps);
	clReleaseMemObject(gpu_temps);
	clReleaseMemObject(gpu_params);
	free(params);
	clReleaseMemObject(gpu_errors);
	clReleaseMemObject(gpu_pixels);
	free(buf2);
}
#ifdef ENABLE_CONSOLE_MAIN_TEST
#include"lodepng.h"
#include"stb_image.h"
int main(int argc, char **argv)
{
	//int arr[]={1, 0, 0, 0, 0};
	//memfill(arr+1, arr, (5-1)*sizeof(int), sizeof(int));
	//
	//pause();
	//return 0;
	int iw, ih, nch0;
	unsigned char *im0=stbi_load(argc==2?argv[1]:

		//"D:/ML/dataset-CLIC/2048x1320_hoach-le-dinh-91879.png"
		//"D:/ML/dataset-kodak/kodim04.png"
		"D:/ML/dataset-kodak/kodim05.png"
		//"D:/ML/dataset-kodak-small/05.PNG"
		//"D:/ML/dataset-kodak/kodim13.png"
		//"D:/ML/dataset-kodak-small/13.PNG"

		, &iw, &ih, &nch0, 4);
	if(!im0)
	{
		printf("Cannot open image\n");
		pause();
		return 0;
	}
	addhalf(im0, iw, ih, 3, 4);
	{
		TimeInfo ti;
		double t_start=time_ms();

		pred_learned_gpu(im0, iw, ih, 1);
		//pred_learned_gpu(im0, iw, ih, 0);

		parsetimedelta(time_ms()-t_start, &ti);
		printf("Done %02d-%02d-%06.3f\n", ti.hours, ti.mins, ti.secs);
	}
	addhalf(im0, iw, ih, 3, 4);

	acme_strftime(g_buf, G_BUF_SIZE, "%Y%m%d-%H%M%S.PNG");
	lodepng_encode_file(g_buf, im0, iw, ih, LCT_RGBA, 8);
	//lodepng_encode_file("result.png", im0, iw, ih, LCT_RGBA, 8);

	free(im0);
	printf("Done.\n");
	pause();
	return 0;
}
#endif//ENABLE_CONSOLE_MAIN_TEST
#endif


//GPU training for CUSTOM3
#if 1
#define C3_NTHREADS 4

static short c3_allparams[C3_NTHREADS*C3_NPARAMS];
static float c3_invCRs[C3_NTHREADS*3];
void custom3_opt_gpu(const char *src, int iw, int ih, Custom3Params *srcparams, int niter, int maskbits, int loud)
{
	static int call_idx=0;
	++call_idx;
	double t_start=time_ms();
	int res=iw*ih;
	if(!init_ocl())
		return;
	
	int error=0;
	cl_mem gpu_indices=clCreateBuffer(context, CL_MEM_READ_ONLY, 3*sizeof(int), 0, &error); CHECKCL(error);
	cl_mem gpu_params=clCreateBuffer(context, CL_MEM_READ_ONLY, C3_NTHREADS*C3_NPARAMS*sizeof(short), 0, &error); CHECKCL(error);
	cl_mem gpu_pixels=clCreateBuffer(context, CL_MEM_READ_ONLY, res*4*sizeof(char), 0, &error); CHECKCL(error);
	cl_mem gpu_allerrors=clCreateBuffer(context, CL_MEM_READ_WRITE, C3_NTHREADS*(C3_REACH+1)*iw*sizeof(char), 0, &error); CHECKCL(error);
	cl_mem gpu_neighbors=clCreateBuffer(context, CL_MEM_READ_WRITE, (C3_NNB+2)*sizeof(short), 0, &error); CHECKCL(error);
	cl_mem gpu_histograms=clCreateBuffer(context, CL_MEM_READ_WRITE, C3_NTHREADS*768*sizeof(int), 0, &error); CHECKCL(error);
	cl_mem gpu_invCRs=clCreateBuffer(context, CL_MEM_WRITE_ONLY, C3_NTHREADS*3*sizeof(float), 0, &error); CHECKCL(error);
	if(!gpu_indices||!gpu_params||!gpu_pixels||!gpu_allerrors||!gpu_neighbors||!gpu_histograms||!gpu_invCRs)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	if(call_idx==1)
	{
		short *p0=(short*)srcparams;
		int sum=0;
		for(int k=0;k<C3_NPARAMS;++k)
			sum+=abs(p0[k]);
		if(!sum)//randomize if starting from all zeros
		{
			for(int k=0;k<C3_NPARAMS;++k)
				p0[k]=(rand()&0x1F)-0x10;
		}
	}
	memfill(c3_allparams, srcparams, C3_NTHREADS*C3_NPARAMS*sizeof(short), C3_NPARAMS*sizeof(short));

	int indices[]={C3_REACH, iw, ih};
	error=clEnqueueWriteBuffer(commandqueue, gpu_indices, CL_FALSE, 0, 3*sizeof(int), indices, 0, 0, 0); CHECKCL(error);
	error=clEnqueueWriteBuffer(commandqueue, gpu_pixels, CL_FALSE, 0, res*4*sizeof(char), src, 0, 0, 0); CHECKCL(error);

	error=clSetKernelArg(kernels[OCL_custom3_eval], 0, sizeof(cl_mem), &gpu_indices);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 1, sizeof(cl_mem), &gpu_params);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 2, sizeof(cl_mem), &gpu_pixels);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 3, sizeof(cl_mem), &gpu_allerrors);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 4, sizeof(cl_mem), &gpu_neighbors);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 5, sizeof(cl_mem), &gpu_histograms);
	error=clSetKernelArg(kernels[OCL_custom3_eval], 6, sizeof(cl_mem), &gpu_invCRs);
	
	if(!niter)
		niter=C3_NPARAMS*10;
	if(loud)
		srand((unsigned)__rdtsc());//
	for(int it=0;it<niter;++it)
	{
		int idx=it%C3_NPARAMS;
		int delta[C3_NTHREADS];
		for(int k=0;k<C3_NTHREADS;++k)
			delta[k]=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1));
		if(!it)
			delta[0]=0;

		for(int k=!it;k<C3_NTHREADS;++k)
			c3_allparams[C3_NPARAMS*k+idx]+=delta[k];
		error=clEnqueueWriteBuffer(commandqueue, gpu_params, CL_FALSE, 0, C3_NTHREADS*C3_NPARAMS*sizeof(short), c3_allparams, 0, 0, 0); CHECKCL(error);
		error=clFlush(commandqueue);	CHECKCL(error);
		error=clFinish(commandqueue);	CHECKCL(error);

		size_t globalsize=C3_NTHREADS;
		error=clEnqueueNDRangeKernel(commandqueue, kernels[OCL_custom3_eval], 1, 0, &globalsize, 0, 0, 0, 0);	CHECKCL(error);
		error=clFlush(commandqueue);	CHECKCL(error);
		error=clFinish(commandqueue);	CHECKCL(error);

		error=clEnqueueReadBuffer(commandqueue, gpu_invCRs, CL_TRUE, 0, C3_NTHREADS*3*sizeof(float), c3_invCRs, 0, 0, 0);
		int bestidx=0;
		float bestloss=c3_invCRs[0]+c3_invCRs[1]+c3_invCRs[2];
		for(int k=1;k<C3_NTHREADS;++k)//select best params
		{
			float loss=c3_invCRs[3*k]+c3_invCRs[3*k+1]+c3_invCRs[3*k+2];
			if(bestloss>loss)
				bestloss=loss, bestidx=k;
		}
		if(bestidx)//populate best params
			memfill(c3_allparams, c3_allparams+C3_NPARAMS*bestidx, bestidx*C3_NPARAMS*sizeof(short), C3_NPARAMS*sizeof(short));
		if(bestidx<C3_NTHREADS-1)
			memfill(c3_allparams+C3_NPARAMS, c3_allparams+C3_NPARAMS*bestidx, (C3_NTHREADS-1-bestidx)*C3_NPARAMS*sizeof(short), C3_NPARAMS*sizeof(short));

		if(loud)
			set_window_title("%4d/%4d: TRGB %lf  %lf %lf %lf  %d", it+1, niter, 3/bestloss, 1/c3_invCRs[3*bestidx], 1/c3_invCRs[3*bestidx+1], 1/c3_invCRs[3*bestidx+2], call_idx);
	}
	memcpy(srcparams, c3_allparams, sizeof(*srcparams));

	error=clReleaseMemObject(gpu_indices); CHECKCL(error);
	error=clReleaseMemObject(gpu_params); CHECKCL(error);
	error=clReleaseMemObject(gpu_pixels); CHECKCL(error);
	error=clReleaseMemObject(gpu_allerrors); CHECKCL(error);
	error=clReleaseMemObject(gpu_neighbors); CHECKCL(error);
	error=clReleaseMemObject(gpu_histograms); CHECKCL(error);
	error=clReleaseMemObject(gpu_invCRs); CHECKCL(error);
}
#endif


#endif//ALLOW_OPENCL