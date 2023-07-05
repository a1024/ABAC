# CMakeList.txt : CMake project for e2, include source and define
# project specific logic here.
#

#find_package(OpenCL REQUIRED)
#include_directories(
#	"E:/C/OpenCL-Headers"
#)
#link_directories(
#	"E:/C/OpenCL-Headers"
#)
add_executable(e2
	"e2.c"
	"format.c"
	"lodepng.c"
	"tests.c"
	"transforms_cpu.c"
	"transforms_gpu.c"
	"transforms.c"
	"util.c"
)
#target_link_libraries(e2 OpenCL::OpenCL)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET e2 PROPERTY CXX_STANDARD 20)
endif()

