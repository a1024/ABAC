# CMakeList.txt : CMake project for e2, include source and define project specific logic here.


include_directories(
	"C:/Projects/a/OpenCL-Headers"
	"C:/Projects/Vulkan-Headers/include"

	"E:/C/OpenCL-Headers"
	"E:/C/a/Vulkan-Headers/include"
)
link_directories(
	"C:/Projects/a/OpenCL-Headers"
	"C:/Projects/Vulkan-Headers"

	"E:/C/OpenCL-Headers"
	"E:/C/a/Vulkan-Headers"
)
#find_package(Vulkan REQUIRED)
#find_package(OpenCL REQUIRED)


add_executable(e2
	e2.c
	format.c
	lodepng.c
	slic.c
	slic2.c
	tests.c
	tests2.c
	tests3.c
	tests4.c
	tests5.c
	transforms.c
	transforms_cpu.c
	transforms_ocl.c
	transforms_vk.c
	util.c
)
#target_link_libraries(e2 OpenCL::OpenCL)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET e2 PROPERTY CXX_STANDARD 20)
endif()

