# CMakeList.txt : CMake project for e2, include source and define project specific logic here.


#include_directories(
#	"C:/Projects/a/OpenCL-Headers"
#	"C:/Projects/Vulkan-Headers/include"
#
#	"E:/C/OpenCL-Headers"
#	"E:/C/a/Vulkan-Headers/include"
#)
#link_directories(
#	"C:/Projects/a/OpenCL-Headers"
#	"C:/Projects/Vulkan-Headers"
#
#	"E:/C/OpenCL-Headers"
#	"E:/C/a/Vulkan-Headers"
#)
#find_package(Vulkan REQUIRED)
#find_package(OpenCL REQUIRED)


add_executable(e2
	e2.c
	format.c
	lodepng.c
	#slic.c
	#slic2.c
	#t48.c
	#t49.c		#excluded on 2024-04-06
	t51.c
	#t52.c		#excluded on 2024-04-06
	#t53.c		#excluded on 2024-04-06
	t54.c
	t55.c
	#tests.c
	#tests2.c
	#tests3.c
	tests4.c
	tests5.c
	tests6.c	# == t47
	transforms.c
	#transforms_cpu.c
	#transforms_ocl.c
	#transforms_vk.c
	util.c
)

#if(MSVC)
#	target_compile_options(e2 PRIVATE /W4)
#else()
#	target_compile_options(e2 PRIVATE -Wall -Wextra -Wpedantic)
#endif()

if(!MSVC)
	target_link_libraries(e2 m)
endif()
#target_link_libraries(e2 OpenCL::OpenCL)

if (CMAKE_VERSION VERSION_GREATER 3.12)
	set_property(TARGET e2 PROPERTY C_STANDARD 11)
#	set_property(TARGET e2 PROPERTY CXX_STANDARD 20)
endif()

