#include"e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#ifdef ALLOW_VULKAN
#include<vulkan/vulkan.h>
#pragma comment(lib, "vulkan-1.lib")
#endif
static const char file[]=__FILE__;

#ifdef ALLOW_VULKAN
static int vulkan_initialized=0;
static VkInstance instance=0;
static VkPhysicalDevice physicaldevice=0;
static VkDevice device=0;

const char* vkerror2str(int error)
{
	switch(error)
	{
#define CASE(E) case E:return #E
	CASE(VK_SUCCESS);
	CASE(VK_NOT_READY);
	CASE(VK_TIMEOUT);
	CASE(VK_EVENT_SET);
	CASE(VK_EVENT_RESET);
	CASE(VK_INCOMPLETE);
	CASE(VK_ERROR_OUT_OF_HOST_MEMORY);
	CASE(VK_ERROR_OUT_OF_DEVICE_MEMORY);
	CASE(VK_ERROR_INITIALIZATION_FAILED);
	CASE(VK_ERROR_DEVICE_LOST);
	CASE(VK_ERROR_MEMORY_MAP_FAILED);
	CASE(VK_ERROR_LAYER_NOT_PRESENT);
	CASE(VK_ERROR_EXTENSION_NOT_PRESENT);
	CASE(VK_ERROR_FEATURE_NOT_PRESENT);
	CASE(VK_ERROR_INCOMPATIBLE_DRIVER);
	CASE(VK_ERROR_TOO_MANY_OBJECTS);
	CASE(VK_ERROR_FORMAT_NOT_SUPPORTED);
	CASE(VK_ERROR_FRAGMENTED_POOL);
	CASE(VK_ERROR_UNKNOWN);
	CASE(VK_ERROR_OUT_OF_POOL_MEMORY);
	CASE(VK_ERROR_INVALID_EXTERNAL_HANDLE);
	CASE(VK_ERROR_FRAGMENTATION);
	CASE(VK_ERROR_INVALID_OPAQUE_CAPTURE_ADDRESS);
	CASE(VK_PIPELINE_COMPILE_REQUIRED);
	CASE(VK_ERROR_SURFACE_LOST_KHR);
	CASE(VK_ERROR_NATIVE_WINDOW_IN_USE_KHR);
	CASE(VK_SUBOPTIMAL_KHR);
	CASE(VK_ERROR_OUT_OF_DATE_KHR);
	CASE(VK_ERROR_INCOMPATIBLE_DISPLAY_KHR);
	CASE(VK_ERROR_VALIDATION_FAILED_EXT);
	CASE(VK_ERROR_INVALID_SHADER_NV);
	CASE(VK_ERROR_IMAGE_USAGE_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_VIDEO_PICTURE_LAYOUT_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_VIDEO_PROFILE_OPERATION_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_VIDEO_PROFILE_FORMAT_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_VIDEO_PROFILE_CODEC_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_VIDEO_STD_VERSION_NOT_SUPPORTED_KHR);
	CASE(VK_ERROR_INVALID_DRM_FORMAT_MODIFIER_PLANE_LAYOUT_EXT);
	CASE(VK_ERROR_NOT_PERMITTED_KHR);
	CASE(VK_ERROR_FULL_SCREEN_EXCLUSIVE_MODE_LOST_EXT);
	CASE(VK_THREAD_IDLE_KHR);
	CASE(VK_THREAD_DONE_KHR);
	CASE(VK_OPERATION_DEFERRED_KHR);
	CASE(VK_OPERATION_NOT_DEFERRED_KHR);
	CASE(VK_ERROR_COMPRESSION_EXHAUSTED_EXT);
	CASE(VK_ERROR_INCOMPATIBLE_SHADER_BINARY_EXT);
#undef CASE
	}
	return "Unknown error";
}
#define CHECKVK(E) (!(E)||LOG_ERROR(vkerror2str(E)))

#if 0
const char hellovksrc[]=
	"#version 460\n"
	"#extension GL_EXT_debug_printf : require\n"
	"void main()\n"
	"{\n"
	"    debugPrintfEXT(\"\\\'Hello World\\\' (said thread %d)\\n\", gl_GlobalInvocationID.x);\n"
	"}\n"
;
#endif


//sheredom/VkComputeSample		https://gist.github.com/sheredom/523f02bbad2ae397d7ed255f3f3b5a7f
#if 0
VkResult getBestTransferQueueNPH(VkPhysicalDevice physicalDevice, uint32_t* queueFamilyIndex)
{
	uint32_t queueFamilyPropertiesCount=0;
	vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyPropertiesCount, 0);
	VkQueueFamilyProperties* const queueFamilyProperties = (VkQueueFamilyProperties*)_malloca(sizeof(VkQueueFamilyProperties) * queueFamilyPropertiesCount);
	if(!queueFamilyProperties)
	{
		LOG_ERROR("");
		return VK_ERROR_OUT_OF_HOST_MEMORY;
	}
	vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyPropertiesCount, queueFamilyProperties);

	for(uint32_t i=0;i<queueFamilyPropertiesCount;++i)//first try and find a queue that has just the transfer bit set
	{
		const VkQueueFlags maskedFlags=queueFamilyProperties[i].queueFlags&~VK_QUEUE_SPARSE_BINDING_BIT;//mask out the sparse binding bit that we aren't caring about (yet!)
		if (!(maskedFlags&(VK_QUEUE_GRAPHICS_BIT|VK_QUEUE_COMPUTE_BIT)) && (maskedFlags&VK_QUEUE_TRANSFER_BIT))
		{
			*queueFamilyIndex=i;
			return VK_SUCCESS;
		}
	}
	for(uint32_t i=0;i<queueFamilyPropertiesCount;++i)//otherwise we'll prefer using a compute-only queue, remember that having compute on the queue implicitly enables transfer!
	{
		const VkQueueFlags maskedFlags=queueFamilyProperties[i].queueFlags&~VK_QUEUE_SPARSE_BINDING_BIT;//mask out the sparse binding bit that we aren't caring about (yet!)
		if(!(maskedFlags&VK_QUEUE_GRAPHICS_BIT) && (maskedFlags&VK_QUEUE_COMPUTE_BIT))
		{
			*queueFamilyIndex=i;
			return VK_SUCCESS;
		}
	}
	for(uint32_t i=0;i<queueFamilyPropertiesCount;++i)//lastly get any queue that'll work for us (graphics, compute or transfer bit set)
	{
		const VkQueueFlags maskedFlags=queueFamilyProperties[i].queueFlags&~VK_QUEUE_SPARSE_BINDING_BIT;// mask out the sparse binding bit that we aren't caring about (yet!)
		if(maskedFlags&(VK_QUEUE_GRAPHICS_BIT|VK_QUEUE_COMPUTE_BIT|VK_QUEUE_TRANSFER_BIT))
		{
			*queueFamilyIndex=i;
			return VK_SUCCESS;
		}
	}
	return VK_ERROR_INITIALIZATION_FAILED;
}
VkResult getBestComputeQueueNPH(VkPhysicalDevice physicalDevice, uint32_t* queueFamilyIndex)
{
	uint32_t queueFamilyPropertiesCount=0;
	vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyPropertiesCount, 0);
	VkQueueFamilyProperties* const queueFamilyProperties = (VkQueueFamilyProperties*)_malloca(sizeof(VkQueueFamilyProperties) * queueFamilyPropertiesCount);
	if(!queueFamilyProperties)
	{
		LOG_ERROR("");
		return VK_ERROR_OUT_OF_HOST_MEMORY;
	}
	vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyPropertiesCount, queueFamilyProperties);

	for(uint32_t i=0;i<queueFamilyPropertiesCount;++i)//first try and find a queue that has just the compute bit set
	{
		const VkQueueFlags maskedFlags=queueFamilyProperties[i].queueFlags&~(VK_QUEUE_TRANSFER_BIT|VK_QUEUE_SPARSE_BINDING_BIT);// mask out the sparse binding bit that we aren't caring about (yet!) and the transfer bit
		if(!(maskedFlags&VK_QUEUE_GRAPHICS_BIT) && (maskedFlags&VK_QUEUE_COMPUTE_BIT))
		{
			*queueFamilyIndex=i;
			return VK_SUCCESS;
		}
	}
	for(uint32_t i=0;i<queueFamilyPropertiesCount;++i)//lastly get any queue that'll work for us
	{
		const VkQueueFlags maskedFlags = queueFamilyProperties[i].queueFlags&~(VK_QUEUE_TRANSFER_BIT|VK_QUEUE_SPARSE_BINDING_BIT);// mask out the sparse binding bit that we aren't caring about (yet!) and the transfer bit
		if(maskedFlags&VK_QUEUE_COMPUTE_BIT)
		{
			*queueFamilyIndex=i;
			return VK_SUCCESS;
		}
	}
	return VK_ERROR_INITIALIZATION_FAILED;
}
#endif

int init_vk()
{
	if(vulkan_initialized)
		return 1;
	
	//Introduction to Vulkan Compute Shaders	https://www.youtube.com/watch?v=KN9nHo9kvZs		https://necrashter.github.io/ceng469/project/compute-shaders-intro
#if 1
	VkResult error;
	VkApplicationInfo appinfo=
	{
		VK_STRUCTURE_TYPE_APPLICATION_INFO, 0,
		"e2", 1,
		"v8", 1,
		VK_API_VERSION_1_3,
	};
	const char *layers[]=
	{
		"VK_LAYER_KHRONOS_validation",
	};
	VkInstanceCreateInfo instancecreateinfo=
	{
		VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, 0, 0,
		&appinfo,
		1, layers,//layers
		0, 0,//extensions
	};
	error=vkCreateInstance(&instancecreateinfo, 0, &instance);	CHECKVK(error);
	if(!instance)
		return 0;


	unsigned count=1;
	error=vkEnumeratePhysicalDevices(instance, &count, &physicaldevice);	CHECKVK(error);
	if(!physicaldevice)
		return 0;
	VkPhysicalDeviceProperties properties;
	vkGetPhysicalDeviceProperties(physicaldevice, &properties);
	printf("Device name: %s\n", properties.deviceName);
	printf("Vulkan version: %d.%d.%d\n", VK_VERSION_MAJOR(properties.apiVersion), VK_VERSION_MINOR(properties.apiVersion), VK_VERSION_PATCH(properties.apiVersion));

	
	int queuefamilyidx=-1;
	count=0;
	vkGetPhysicalDeviceQueueFamilyProperties(physicaldevice, &count, 0);
	ArrayHandle queuefamilies;
	ARRAY_ALLOC(VkQueueFamilyProperties, queuefamilies, 0, count, 0, 0);
	vkGetPhysicalDeviceQueueFamilyProperties(physicaldevice, &count, (VkQueueFamilyProperties*)queuefamilies->data);
	for(int k=0;k<queuefamilies->count;++k)
	{
		VkQueueFamilyProperties *queuefamily=(VkQueueFamilyProperties*)array_at(&queuefamilies, k);
		if(queuefamily->queueFlags&VK_QUEUE_COMPUTE_BIT)
		{
			queuefamilyidx=k;
			break;
		}
	}
	if(queuefamilyidx==-1)
	{
		LOG_ERROR("No compute queue family");
		return 0;
	}

	
	float priority=1;
	VkDeviceQueueCreateInfo queuecreateinfo=
	{
		VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO, 0, 0,
		queuefamilyidx, 1, &priority,
	};
	VkDeviceCreateInfo devicecreateinfo=
	{
		VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO, 0, 0,
		1, &queuecreateinfo,//queue infos
		0, 0,//enabled layers
		0, 0,//enabled extensions
		0,//enabled features
	};
	error=vkCreateDevice(physicaldevice, &devicecreateinfo, 0, &device);	CHECKVK(error);

	
	VkBuffer inbuffer, outbuffer;
	int bufelems=10;
	int bufsize=bufelems*sizeof(int);
	VkBufferCreateInfo buffercreateinfo=
	{
		VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO, 0, 0,
		bufsize,
		VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
		VK_SHARING_MODE_EXCLUSIVE,
		1, &queuefamilyidx,
	};
	error=vkCreateBuffer(device, &buffercreateinfo, 0, &inbuffer);	CHECKVK(error);
	error=vkCreateBuffer(device, &buffercreateinfo, 0, &outbuffer);	CHECKVK(error);

	
	int memorytypeidx=-1;
	VkPhysicalDeviceMemoryProperties memoryproperties;
	vkGetPhysicalDeviceMemoryProperties(physicaldevice, &memoryproperties);
	for(unsigned k=0;k<memoryproperties.memoryTypeCount;++k)
	{
		VkMemoryType *mem=memoryproperties.memoryTypes+k;
		if((mem->propertyFlags&VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT)&&(mem->propertyFlags&VK_MEMORY_PROPERTY_HOST_COHERENT_BIT))
		{
			memorytypeidx=k;
			break;
		}
	}


	//allocate memory
	VkMemoryRequirements2 inmemoryrequirements={VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2}, outmemoryrequirements={VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2};
	VkDeviceBufferMemoryRequirements buffermemoryrequirements=
	{
		VK_STRUCTURE_TYPE_DEVICE_BUFFER_MEMORY_REQUIREMENTS, 0,
		&buffercreateinfo,
	};
	vkGetDeviceBufferMemoryRequirements(device, &buffermemoryrequirements, &inmemoryrequirements);
	vkGetDeviceBufferMemoryRequirements(device, &buffermemoryrequirements, &outmemoryrequirements);

	VkDeviceMemory inmemory=0, outmemory=0;
	VkMemoryAllocateInfo inmemoryallocateinfo=
	{
		VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO, 0,
		inmemoryrequirements.memoryRequirements.size,
		memorytypeidx,
	};
	VkMemoryAllocateInfo outmemoryallocateinfo=
	{
		VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO, 0,
		outmemoryrequirements.memoryRequirements.size,
		memorytypeidx,
	};
	error=vkAllocateMemory(device, &inmemoryallocateinfo, 0, &inmemory);	CHECKVK(error);
	error=vkAllocateMemory(device, &outmemoryallocateinfo, 0, &outmemory);	CHECKVK(error);


	{
		unsigned *ptr=0;
		error=vkMapMemory(device, inmemory, 0, bufsize, 0, (void**)&ptr);	CHECKVK(error);
		for(int k=0;k<bufelems;++k)
			ptr[k]=k;
		vkUnmapMemory(device, inmemory);
	}

	vkBindBufferMemory(device, inbuffer, inmemory, 0);
	vkBindBufferMemory(device, outbuffer, outmemory, 0);

	
	ArrayHandle spirv=0;
	{
		//ptrdiff_t size=get_filesize("comp.spv");
		//if(size==-1)
		//{
			const char *srcpath="shader.comp";
			ptrdiff_t size=get_filesize(srcpath);
			if(size==-1)
			{
				srcpath="E:/C/e2/e2/shader.comp";
				size=get_filesize(srcpath);
				if(size==-1)
				{
					LOG_ERROR("Cannot find \'shader.comp\'");
					return 0;
				}
			}
			snprintf(g_buf, G_BUF_SIZE, "glslc \"%s\" -o comp.spv", srcpath);
			int err2=system(g_buf);
			if(err2==-1)
			{
				LOG_ERROR("System command error %d", errno);
				return 0;
			}
		//}
		spirv=load_file("comp.spv", 1, 0);
		if(!spirv)
		{
			LOG_ERROR("Shader failed to compile");
			return 0;
		}
	}

	
	VkShaderModule shader=0;
	VkShaderModuleCreateInfo shadercreateinfo=
	{
		VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO, 0, 0,
		spirv->count,
		(unsigned*)spirv->data,
	};
	error=vkCreateShaderModule(device, &shadercreateinfo, 0, &shader);	CHECKVK(error);


	VkDescriptorSetLayout descriptorsetlayout=0;
	VkDescriptorSetLayoutBinding descriptorsetlayoutbinding[]=
	{
		{
			0,//binding
			VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			1, VK_SHADER_STAGE_COMPUTE_BIT, 0,
		},
		{
			1,
			VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			1, VK_SHADER_STAGE_COMPUTE_BIT, 0,
		},
	};
	VkDescriptorSetLayoutCreateInfo descriptorsetlayoutcreateinfo=
	{
		VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO, 0, 0,
		2, descriptorsetlayoutbinding,
	};
	vkCreateDescriptorSetLayout(device, &descriptorsetlayoutcreateinfo, 0, &descriptorsetlayout);

	
	VkPipelineLayout pipelinelayout=0;
	VkPipelineLayoutCreateInfo pipelinelayoutcreateinfo=
	{
		VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO, 0, 0,
		1, &descriptorsetlayout,//layouts
		0, 0,//push constant ranges
	};
	error=vkCreatePipelineLayout(device, &pipelinelayoutcreateinfo, 0, &pipelinelayout);	CHECKVK(error);
	
	VkPipelineCache pipelinecache=0;
	VkPipelineCacheCreateInfo pipelinecachecreateinfo=
	{
		VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO, 0, 0,
		0, 0,//initial data
	};
	error=vkCreatePipelineCache(device, &pipelinecachecreateinfo, 0, &pipelinecache);	CHECKVK(error);

	VkPipeline pipeline=0;
	VkComputePipelineCreateInfo computepipelinecreateinfo=
	{
		VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO, 0, 0,
		{
			VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, 0, 0,
			VK_SHADER_STAGE_COMPUTE_BIT,
			shader,
			"main",
			0,//specialization info
		},
		pipelinelayout,
		0, 0,//base pipeline handle & index
	};
	error=vkCreateComputePipelines(device, pipelinecache, 1, &computepipelinecreateinfo, 0, &pipeline);	CHECKVK(error);

	
	VkDescriptorPool descriptorpool=0;
	VkDescriptorPoolSize descriptorpoolsize=
	{
		VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
		2,//descriptor count
	};
	VkDescriptorPoolCreateInfo descriptorpoolcreateinfo=
	{
		VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO, 0, 0,
		1,//maxsets > 0
		1, &descriptorpoolsize,
	};
	error=vkCreateDescriptorPool(device, &descriptorpoolcreateinfo, 0, &descriptorpool);	CHECKVK(error);

	
	VkDescriptorSet descriptorset=0;
	VkDescriptorSetAllocateInfo descriptorsetallocateinfo=
	{
		VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO, 0,
		descriptorpool,
		1, &descriptorsetlayout,
	};
	error=vkAllocateDescriptorSets(device, &descriptorsetallocateinfo, &descriptorset);	CHECKVK(error);


	VkDescriptorBufferInfo inbufferinfo={inbuffer, 0, bufsize};
	VkDescriptorBufferInfo outbufferinfo={outbuffer, 0, bufsize};
	VkWriteDescriptorSet writedescriptorsets[]=
	{
		{
			VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, 0,
			descriptorset,
			0,//binding
			0,//array element
			1,//descriptor count
			VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			0,//imageinfo
			&inbufferinfo,
			0,//texelbufferview
		},
		{
			VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, 0,
			descriptorset,
			1,
			0,
			1,
			VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			0,
			&inbufferinfo,
			0,
		},
	};
	vkUpdateDescriptorSets(device, 2, writedescriptorsets, 0, 0);


	//submit work to GPU (glUseProgram)


	VkCommandPool commandpool=0;
	VkCommandPoolCreateInfo commandpoolcreateinfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO, 0, 0,
		queuefamilyidx,
	};
	error=vkCreateCommandPool(device, &commandpoolcreateinfo, 0, &commandpool);	CHECKVK(error);


	VkCommandBuffer commandbuffer=0;
	VkCommandBufferAllocateInfo commandbufferallocinfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO, 0,
		commandpool,
		VK_COMMAND_BUFFER_LEVEL_PRIMARY,
		1,//command buffer count
	};
	error=vkAllocateCommandBuffers(device, &commandbufferallocinfo, &commandbuffer);	CHECKVK(error);

	//record commands
	VkCommandBufferBeginInfo cmdbufferbegininfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, 0, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
		0,//inheritance info
	};
	error=vkBeginCommandBuffer(commandbuffer, &cmdbufferbegininfo);		CHECKVK(error);
	vkCmdBindPipeline(commandbuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
	vkCmdBindDescriptorSets(commandbuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipelinelayout, 0, 1, &descriptorset, 0, 0);
	vkCmdDispatch(commandbuffer, bufelems, 1, 1);
	error=vkEndCommandBuffer(commandbuffer);	CHECKVK(error);


	//submit work & wait
	VkQueue queue=0;
	vkGetDeviceQueue(device, queuefamilyidx, 0, &queue);

	VkFence fence=0;
	VkFenceCreateInfo fencecreateinfo=
	{
		VK_STRUCTURE_TYPE_FENCE_CREATE_INFO, 0, 0,
	};
	error=vkCreateFence(device, &fencecreateinfo, 0, &fence);	CHECKVK(error);

	VkSubmitInfo submitinfo=
	{
		VK_STRUCTURE_TYPE_SUBMIT_INFO, 0,
		0, 0,//wait semaphores
		0,//wait dst stage mask
		1, &commandbuffer,
		0, 0,//signal semaphores
	};
	error=vkQueueSubmit(queue, 1, &submitinfo, fence);		CHECKVK(error);
	error=vkWaitForFences(device, 1, &fence, 1, (size_t)-1);CHECKVK(error);

	{
		unsigned *ptr=0;
		error=vkMapMemory(device, outmemory, 0, bufsize, 0, &ptr);	CHECKVK(error);
		printf("Result:\n");
		for(int k=0;k<bufelems;++k)
			printf("%d, ", ptr[k]);
		printf("\n");
		vkUnmapMemory(device, outmemory);
	}


	//cleanup
	vkDestroyFence(device, fence, 0);	fence=0;
	vkFreeMemory(device, inmemory, 0);			inmemory=0;
	vkFreeMemory(device, outmemory, 0);			outmemory=0;
	vkDestroyShaderModule(device, shader, 0);	shader=0;
	vkDestroyBuffer(device, inbuffer, 0);		inbuffer=0;
	vkDestroyBuffer(device, outbuffer, 0);		outbuffer=0;
	vkDestroyPipelineCache(device, pipelinecache, 0);	pipelinecache=0;
	vkDestroyPipelineLayout(device, pipelinelayout, 0);	pipelinelayout=0;
	vkDestroyPipeline(device, pipeline, 0);		pipeline=0;
	vkDestroyDescriptorSetLayout(device, descriptorsetlayout, 0);	descriptorsetlayout=0;
	vkDestroyDescriptorPool(device, descriptorpool, 0);	descriptorpool=0;
	vkDestroyCommandPool(device, commandpool, 0);	commandpool=0;
	vkDestroyDevice(device, 0);					device=0;

	printf("Done.\n");
	pause();
	exit(0);
#endif

	//sheredom/VkComputeSample		https://gist.github.com/sheredom/523f02bbad2ae397d7ed255f3f3b5a7f
#if 0
	VkResult error;
	const VkApplicationInfo applicationInfo =
	{
		VK_STRUCTURE_TYPE_APPLICATION_INFO, 0,
		"VKComputeSample", 0,
		"", 0,
		VK_MAKE_VERSION(1, 0, 9),
	};
	const VkInstanceCreateInfo instanceCreateInfo =
	{
		VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, 0, 0,
		&applicationInfo,
		0, 0,//layers
		0, 0,//extensions
	};
	VkInstance instance;
	error=vkCreateInstance(&instanceCreateInfo, 0, &instance);		CHECKVK(error);
	
	uint32_t physicalDeviceCount = 0;
	error=vkEnumeratePhysicalDevices(instance, &physicalDeviceCount, 0);		CHECKVK(error);
	VkPhysicalDevice *physicalDevices=(VkPhysicalDevice*)malloc(physicalDeviceCount*sizeof(VkPhysicalDevice));
	if(!physicalDevices)
	{
		LOG_ERROR("");
		return 0;
	}
	error=vkEnumeratePhysicalDevices(instance, &physicalDeviceCount, physicalDevices);		CHECKVK(error);
	for(uint32_t i=0;i<physicalDeviceCount;++i)
	{
		uint32_t queueFamilyIndex=0;
		error=getBestComputeQueueNPH(physicalDevices[i], &queueFamilyIndex);		CHECKVK(error);
		const float queuePrioritory=1;
		const VkDeviceQueueCreateInfo deviceQueueCreateInfo =
		{
			VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO, 0, 0,
			queueFamilyIndex,
			1,
			&queuePrioritory,
		};
		const VkDeviceCreateInfo deviceCreateInfo =
		{
			VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO, 0, 0,
			1, &deviceQueueCreateInfo,
			0, 0,//layers
			0, 0,//extensions
			0,//features
		};
		VkDevice device;
		error=vkCreateDevice(physicalDevices[i], &deviceCreateInfo, 0, &device);		CHECKVK(error);


		VkPhysicalDeviceMemoryProperties properties;
		vkGetPhysicalDeviceMemoryProperties(physicalDevices[i], &properties);
		const int32_t bufferLength = 16384;
		const uint32_t bufferSize = bufferLength*sizeof(int32_t);
		const VkDeviceSize memorySize = bufferSize * 2;//we are going to need two buffers from this one memory
		uint32_t memoryTypeIndex = VK_MAX_MEMORY_TYPES;//set memoryTypeIndex to an invalid entry in the properties.memoryTypes array
		for (uint32_t k=0;k<properties.memoryTypeCount;k++)
		{
			if ((VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT & properties.memoryTypes[k].propertyFlags) &&
				(VK_MEMORY_PROPERTY_HOST_COHERENT_BIT & properties.memoryTypes[k].propertyFlags) &&
				(memorySize < properties.memoryHeaps[properties.memoryTypes[k].heapIndex].size))
			{
				memoryTypeIndex=k;
				break;
			}
		}
		if(memoryTypeIndex==VK_MAX_MEMORY_TYPES)
			LOG_ERROR("");


		const VkMemoryAllocateInfo memoryAllocateInfo =
		{
			VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO, 0,
			memorySize,
			memoryTypeIndex,
		};
		VkDeviceMemory memory;
		error=vkAllocateMemory(device, &memoryAllocateInfo, 0, &memory);	CHECKVK(error);


		int32_t *payload;
		error=vkMapMemory(device, memory, 0, memorySize, 0, (void*)&payload);	CHECKVK(error);
		for (uint32_t k = 1; k < memorySize / sizeof(int32_t); k++)
			payload[k] = rand();
		vkUnmapMemory(device, memory);


		const VkBufferCreateInfo bufferCreateInfo=
		{
			VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO, 0, 0,
			bufferSize,
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
			VK_SHARING_MODE_EXCLUSIVE,
			1, &queueFamilyIndex,
		};
		VkBuffer in_buffer;
		error=vkCreateBuffer(device, &bufferCreateInfo, 0, &in_buffer);	CHECKVK(error);
		error=vkBindBufferMemory(device, in_buffer, memory, 0);			CHECKVK(error);

		VkBuffer out_buffer;
		error=vkCreateBuffer(device, &bufferCreateInfo, 0, &out_buffer);	CHECKVK(error);
		error=vkBindBufferMemory(device, out_buffer, memory, bufferSize);	CHECKVK(error);


		enum
		{
			RESERVED_ID = 0,
			FUNC_ID,
			IN_ID,
			OUT_ID,
			GLOBAL_INVOCATION_ID,
			VOID_TYPE_ID,
			FUNC_TYPE_ID,
			INT_TYPE_ID,
			INT_ARRAY_TYPE_ID,
			STRUCT_ID,
			POINTER_TYPE_ID,
			ELEMENT_POINTER_TYPE_ID,
			INT_VECTOR_TYPE_ID,
			INT_VECTOR_POINTER_TYPE_ID,
			INT_POINTER_TYPE_ID,
			CONSTANT_ZERO_ID,
			CONSTANT_ARRAY_LENGTH_ID,
			LABEL_ID,
			IN_ELEMENT_ID,
			OUT_ELEMENT_ID,
			GLOBAL_INVOCATION_X_ID,
			GLOBAL_INVOCATION_X_PTR_ID,
			TEMP_LOADED_ID,
			BOUND
		};
		enum
		{
			INPUT = 1,
			UNIFORM = 2,
			BUFFER_BLOCK = 3,
			ARRAY_STRIDE = 6,
			BUILTIN = 11,
			BINDING = 33,
			OFFSET = 35,
			DESCRIPTOR_SET = 34,
			GLOBAL_INVOCATION = 28,
			OP_TYPE_VOID = 19,
			OP_TYPE_FUNCTION = 33,
			OP_TYPE_INT = 21,
			OP_TYPE_VECTOR = 23,
			OP_TYPE_ARRAY = 28,
			OP_TYPE_STRUCT = 30,
			OP_TYPE_POINTER = 32,
			OP_VARIABLE = 59,
			OP_DECORATE = 71,
			OP_MEMBER_DECORATE = 72,
			OP_FUNCTION = 54,
			OP_LABEL = 248,
			OP_ACCESS_CHAIN = 65,
			OP_CONSTANT = 43,
			OP_LOAD = 61,
			OP_STORE = 62,
			OP_RETURN = 253,
			OP_FUNCTION_END = 56,
			OP_CAPABILITY = 17,
			OP_MEMORY_MODEL = 14,
			OP_ENTRY_POINT = 15,
			OP_EXECUTION_MODE = 16,
			OP_COMPOSITE_EXTRACT = 81,
		};
		int32_t shader[]=
		{
			// first is the SPIR-V header
			0x07230203, // magic header ID
			0x00010000, // version 1.0.0
			0,          // generator (optional)
			BOUND,      // bound
			0,          // schema

			(2 << 16) | OP_CAPABILITY, 1,//OpCapability Shader
			(3 << 16) | OP_MEMORY_MODEL, 0, 0,//OpMemoryModel Logical Simple
			(4 << 16) | OP_ENTRY_POINT, 5, FUNC_ID, 'f',// OpEntryPoint GLCompute %FUNC_ID "f" %IN_ID %OUT_ID
			(6 << 16) | OP_EXECUTION_MODE, FUNC_ID, 17, 1, 1, 1,//OpExecutionMode %FUNC_ID LocalSize 1 1 1

			// next declare decorations
			(3 << 16) | OP_DECORATE, STRUCT_ID, BUFFER_BLOCK, 
			(4 << 16) | OP_DECORATE, GLOBAL_INVOCATION_ID, BUILTIN, GLOBAL_INVOCATION,
			(4 << 16) | OP_DECORATE, IN_ID, DESCRIPTOR_SET, 0,
			(4 << 16) | OP_DECORATE, IN_ID, BINDING, 0,
			(4 << 16) | OP_DECORATE, OUT_ID, DESCRIPTOR_SET, 0,
			(4 << 16) | OP_DECORATE, OUT_ID, BINDING, 1,
			(4 << 16) | OP_DECORATE, INT_ARRAY_TYPE_ID, ARRAY_STRIDE, 4,
			(5 << 16) | OP_MEMBER_DECORATE, STRUCT_ID, 0, OFFSET, 0,

			// next declare types
			(2 << 16) | OP_TYPE_VOID, VOID_TYPE_ID,
			(3 << 16) | OP_TYPE_FUNCTION, FUNC_TYPE_ID, VOID_TYPE_ID,
			(4 << 16) | OP_TYPE_INT, INT_TYPE_ID, 32, 1,
			(4 << 16) | OP_CONSTANT, INT_TYPE_ID, CONSTANT_ARRAY_LENGTH_ID, bufferLength,
			(4 << 16) | OP_TYPE_ARRAY, INT_ARRAY_TYPE_ID, INT_TYPE_ID, CONSTANT_ARRAY_LENGTH_ID,
			(3 << 16) | OP_TYPE_STRUCT, STRUCT_ID, INT_ARRAY_TYPE_ID,
			(4 << 16) | OP_TYPE_POINTER, POINTER_TYPE_ID, UNIFORM, STRUCT_ID,
			(4 << 16) | OP_TYPE_POINTER, ELEMENT_POINTER_TYPE_ID, UNIFORM, INT_TYPE_ID,
			(4 << 16) | OP_TYPE_VECTOR, INT_VECTOR_TYPE_ID, INT_TYPE_ID, 3,
			(4 << 16) | OP_TYPE_POINTER, INT_VECTOR_POINTER_TYPE_ID, INPUT, INT_VECTOR_TYPE_ID,
			(4 << 16) | OP_TYPE_POINTER, INT_POINTER_TYPE_ID, INPUT, INT_TYPE_ID,

			// then declare constants
			(4 << 16) | OP_CONSTANT, INT_TYPE_ID, CONSTANT_ZERO_ID, 0,

			// then declare variables
			(4 << 16) | OP_VARIABLE, POINTER_TYPE_ID, IN_ID, UNIFORM,
			(4 << 16) | OP_VARIABLE, POINTER_TYPE_ID, OUT_ID, UNIFORM,
			(4 << 16) | OP_VARIABLE, INT_VECTOR_POINTER_TYPE_ID, GLOBAL_INVOCATION_ID, INPUT,

			// then declare function
			(5 << 16) | OP_FUNCTION, VOID_TYPE_ID, FUNC_ID, 0, FUNC_TYPE_ID,
			(2 << 16) | OP_LABEL, LABEL_ID,
			(5 << 16) | OP_ACCESS_CHAIN, INT_POINTER_TYPE_ID, GLOBAL_INVOCATION_X_PTR_ID, GLOBAL_INVOCATION_ID, CONSTANT_ZERO_ID,
			(4 << 16) | OP_LOAD, INT_TYPE_ID, GLOBAL_INVOCATION_X_ID, GLOBAL_INVOCATION_X_PTR_ID,
			(6 << 16) | OP_ACCESS_CHAIN, ELEMENT_POINTER_TYPE_ID, IN_ELEMENT_ID, IN_ID, CONSTANT_ZERO_ID, GLOBAL_INVOCATION_X_ID,
			(4 << 16) | OP_LOAD, INT_TYPE_ID, TEMP_LOADED_ID, IN_ELEMENT_ID,
			(6 << 16) | OP_ACCESS_CHAIN, ELEMENT_POINTER_TYPE_ID, OUT_ELEMENT_ID, OUT_ID, CONSTANT_ZERO_ID, GLOBAL_INVOCATION_X_ID,
			(3 << 16) | OP_STORE, OUT_ELEMENT_ID, TEMP_LOADED_ID,
			(1 << 16) | OP_RETURN,
			(1 << 16) | OP_FUNCTION_END,
		};
		VkShaderModuleCreateInfo shaderModuleCreateInfo=
		{
			VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO, 0, 0,
			sizeof(shader),
			shader
		};
		VkShaderModule shader_module;
		error=vkCreateShaderModule(device, &shaderModuleCreateInfo, 0, &shader_module);	CHECKVK(error);


		VkDescriptorSetLayoutBinding descriptorSetLayoutBindings[2]=
		{
			{
				0,
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				1,
				VK_SHADER_STAGE_COMPUTE_BIT,
				0
			},
			{
				1,
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				1,
				VK_SHADER_STAGE_COMPUTE_BIT,
				0
			}
		};
		VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfo=
		{
			VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO, 0, 0,
			2, descriptorSetLayoutBindings,
		};
		VkDescriptorSetLayout descriptorSetLayout;
		error=vkCreateDescriptorSetLayout(device, &descriptorSetLayoutCreateInfo, 0, &descriptorSetLayout);	CHECKVK(error);


		VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo=
		{
			VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO, 0, 0,
			1, &descriptorSetLayout,
			0, 0,
		};
		VkPipelineLayout pipelineLayout;
		error=vkCreatePipelineLayout(device, &pipelineLayoutCreateInfo, 0, &pipelineLayout);	CHECKVK(error);


		VkComputePipelineCreateInfo computePipelineCreateInfo=
		{
			VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO, 0, 0,
			{
				VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, 0, 0,
				VK_SHADER_STAGE_COMPUTE_BIT,
				shader_module,
				"f",
				0,
			},
			pipelineLayout,
			0, 0,
		};
		VkPipeline pipeline;
		error=vkCreateComputePipelines(device, 0, 1, &computePipelineCreateInfo, 0, &pipeline);	CHECKVK(error);


		VkDescriptorPool descriptorPool=0;
		VkDescriptorPoolSize descriptorPoolSize=
		{
			VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			2,
		};
		VkDescriptorPoolCreateInfo descriptorPoolCreateInfo=
		{
			VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO, 0, 0,
			1,
			1, &descriptorPoolSize,
		};
		error=vkCreateDescriptorPool(device, &descriptorPoolCreateInfo, 0, &descriptorPool);	CHECKVK(error);


		VkDescriptorSetAllocateInfo descriptorSetAllocateInfo=
		{
			VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO, 0,
			descriptorPool,
			1, &descriptorSetLayout,
		};
		VkDescriptorSet descriptorSet;
		error=vkAllocateDescriptorSets(device, &descriptorSetAllocateInfo, &descriptorSet);	CHECKVK(error);


		VkDescriptorBufferInfo in_descriptorBufferInfo=
		{
			in_buffer,
			0,
			VK_WHOLE_SIZE
		};
		VkDescriptorBufferInfo out_descriptorBufferInfo=
		{
			out_buffer,
			0,
			VK_WHOLE_SIZE
		};
		VkWriteDescriptorSet writeDescriptorSet[2]=
		{
			{
				VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, 0,
				descriptorSet,	//dstSet
				0,		//dstBinding
				0,		//dstArrayElement
				1,		//descriptorCount
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				0,		//imageInfo
				&in_descriptorBufferInfo,
				0,		//texelBufferView
			},
			{
				VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, 0,
				descriptorSet,
				1,
				0,
				1,
				VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
				0,
				&out_descriptorBufferInfo,
				0
			}
		};
		vkUpdateDescriptorSets(device, 2, writeDescriptorSet, 0, 0);

		
		VkCommandPool commandPool=0;
		VkCommandPoolCreateInfo commandPoolCreateInfo=
		{
			VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO, 0, 0,
			queueFamilyIndex,
		};
		error=vkCreateCommandPool(device, &commandPoolCreateInfo, 0, &commandPool);	CHECKVK(error);


		VkCommandBufferAllocateInfo commandBufferAllocateInfo=
		{
			VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO, 0,
			commandPool,
			VK_COMMAND_BUFFER_LEVEL_PRIMARY,
			1,
		};
		VkCommandBuffer commandBuffer;
		error=vkAllocateCommandBuffers(device, &commandBufferAllocateInfo, &commandBuffer);	CHECKVK(error);


		VkCommandBufferBeginInfo commandBufferBeginInfo=
		{
			VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, 0,
			VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
			0,
		};
		error=vkBeginCommandBuffer(commandBuffer, &commandBufferBeginInfo);	CHECKVK(error);


		vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, 0, 1, &descriptorSet, 0, 0);
		vkCmdDispatch(commandBuffer, bufferSize / sizeof(int32_t), 1, 1);
		error=vkEndCommandBuffer(commandBuffer);	CHECKVK(error);


		VkQueue queue;
		vkGetDeviceQueue(device, queueFamilyIndex, 0, &queue);


		VkSubmitInfo submitInfo=
		{
			VK_STRUCTURE_TYPE_SUBMIT_INFO, 0,
			0, 0,//wait semaphores
			0,
			1, &commandBuffer,
			0, 0,//signal semaphores
		};
		error=vkQueueSubmit(queue, 1, &submitInfo, 0);		CHECKVK(error);
		error=vkQueueWaitIdle(queue);		CHECKVK(error);

		error=vkMapMemory(device, memory, 0, memorySize, 0, (void *)&payload);	CHECKVK(error);
		for(uint32_t k=0, e=bufferSize/sizeof(int32_t);k<e;++k)
		{
			if(payload[k+e]!=payload[k])
				LOG_ERROR("VK_ERROR_OUT_OF_HOST_MEMORY");
		}
	}
	printf("Done.\n");
	pause();
	exit(0);
#endif

	//Vulkanised 2023: Transitioning to Vulkan for Compute		https://www.youtube.com/watch?v=0Fv0jDZUyE4
#if 0
	VkResult error;
	VkApplicationInfo appinfo=
	{
		VK_STRUCTURE_TYPE_APPLICATION_INFO, 0,
		"e2", 1,
		"v8", 1,
		VK_API_VERSION_1_3,
	};
	const char *layers[]=
	{
		"VK_LAYER_KHRONOS_validation",
	};
	VkInstanceCreateInfo instancecreateinfo=
	{
		VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, 0, 0,
		&appinfo,
		1, layers,//layers
		0, 0,//extensions
	};
	error=vkCreateInstance(&instancecreateinfo, 0, &instance);	CHECKVK(error);
	if(!instance)
		return 0;


	int count=1;
	error=vkEnumeratePhysicalDevices(instance, &count, &physicaldevice);	CHECKVK(error);
	if(!physicaldevice)
		return 0;


	VkQueueFamilyProperties *family=0;
	int familyidx;
	count=0;
	vkGetPhysicalDeviceQueueFamilyProperties(physicaldevice, &count, 0);
	ArrayHandle queuefamilies;
	ARRAY_ALLOC(VkQueueFamilyProperties, queuefamilies, 0, count, 0, 0);
	vkGetPhysicalDeviceQueueFamilyProperties(physicaldevice, &count, (VkQueueFamilyProperties*)queuefamilies->data);
	for(familyidx=0;familyidx<queuefamilies->count;++familyidx)
	{
		family=(VkQueueFamilyProperties*)array_at(&queuefamilies, familyidx);
		if(family->queueFlags&VK_QUEUE_COMPUTE_BIT)
			break;
	}
	if(familyidx==queuefamilies->count)
	{
		LOG_ERROR("No compute queue family");
		return 0;
	}


	float priority=1;
	VkDeviceQueueCreateInfo queuecreateinfo=
	{
		VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO, 0, 0,
		familyidx,
		1,
		&priority,
	};
	VkDeviceCreateInfo devicecreateinfo=
	{
		VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO, 0, 0,
		1, &queuecreateinfo,//queue infos
		0, 0,//enabled layers
		0, 0,//enabled extensions
		0,//enabled features
	};
	error=vkCreateDevice(physicaldevice, &devicecreateinfo, 0, &device);	CHECKVK(error);


	//vkCompileGlslToSpv
	ArrayHandle spirv=0;
	{
		//ptrdiff_t size=get_filesize("comp.spv");
		//if(size==-1)
		//{
			const char *srcpath="shader.comp";
			ptrdiff_t size=get_filesize(srcpath);
			if(size==-1)
			{
				srcpath="E:/C/e2/e2/shader.comp";
				size=get_filesize(srcpath);
				if(size==-1)
				{
					LOG_ERROR("Cannot find \'shader.comp\'");
					return 0;
				}
			}
			snprintf(g_buf, G_BUF_SIZE, "glslc \"%s\" -o comp.spv", srcpath);
			int err2=system(g_buf);
			if(err2==-1)
			{
				LOG_ERROR("System command error %d", errno);
				return 0;
			}
		//}
		spirv=load_file("comp.spv", 1, 0);
		if(!spirv)
		{
			LOG_ERROR("Shader failed to compile");
			return 0;
		}
	}

	VkShaderModuleCreateInfo shadercreateinfo=
	{
		VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO, 0, 0,
		spirv->count,
		(unsigned*)spirv->data,
	};
	VkShaderModule shader=0;
	error=vkCreateShaderModule(device, &shadercreateinfo, 0, &shader);	CHECKVK(error);


	VkPipelineLayout pipelinelayout=0;
	VkPipelineLayoutCreateInfo pipelinelayoutcreateinfo=
	{
		VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO, 0, 0,
		0, 0,//layouts
		0, 0,//push constant ranges
	};
	error=vkCreatePipelineLayout(device, &pipelinelayoutcreateinfo, 0, &pipelinelayout);	CHECKVK(error);

	VkPipelineCache pipelinecache=0;
	VkPipelineCacheCreateInfo pipelinecachecreateinfo=
	{
		VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO, 0, 0,
		0,//initial data size
		0,//initial data
	};
	error=vkCreatePipelineCache(device, &pipelinecachecreateinfo, 0, &pipelinecache);	CHECKVK(error);
	
	VkPipeline pipeline;
	VkComputePipelineCreateInfo computepipelinecreateinfo=
	{
		VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO, 0, 0,
		{
			VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, 0, 0,
			VK_SHADER_STAGE_COMPUTE_BIT,
			shader,
			"main",
			0,//specialization info
		},
		pipelinelayout,
		0, 0,//base pipeline handle & index
	};
	error=vkCreateComputePipelines(device, pipelinecache, 1, &computepipelinecreateinfo, 0, &pipeline);	CHECKVK(error);


	VkCommandPoolCreateInfo commandpoolcreateinfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO, 0, 0,
		familyidx,
	};
	VkCommandPool commandpool=0;
	error=vkCreateCommandPool(device, &commandpoolcreateinfo, 0, &commandpool);	CHECKVK(error);


	VkCommandBufferAllocateInfo commandbufferallocateinfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO, 0,
		commandpool,
		VK_COMMAND_BUFFER_LEVEL_PRIMARY,
		1,
	};
	VkCommandBuffer commandbuffer=0;
	error=vkAllocateCommandBuffers(device, &commandbufferallocateinfo, &commandbuffer);	CHECKVK(error);


	VkCommandBufferBeginInfo commandbufferbegininfo=
	{
		VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, 0, 0,
		0,//inheritance info
	};
	error=vkBeginCommandBuffer(commandbuffer, &commandbufferbegininfo);	CHECKVK(error);

	vkCmdBindPipeline(commandbuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);

	vkCmdDispatch(commandbuffer, 8, 1, 1);
	error=vkEndCommandBuffer(commandBuffer);	CHECKVK(error);

	VkQueue queue=0;
	vkGetDeviceQueue(device, familyidx, 0, &queue);
	VkSubmitInfo submitinfo=
	{
		VK_STRUCTURE_TYPE_SUBMIT_INFO, 0,
		0, 0,//wait semaphores
		0,//wait stage mask
		1, &commandbuffer,
		0, 0,//signal semaphores
	};
	vkQueueSubmit(queue, 1, &submitinfo, 0);
	vkDeviceWaitIdle(device);


	array_free(&spirv);
	array_free(&queuefamilies);
#endif

	vulkan_initialized=1;
	return 1;
}
#endif