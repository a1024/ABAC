#include_directories(
#	"C:/Projects/a/OpenCL-Headers"
#)
#link_directories(
#	"C:/Projects/a/OpenCL-Headers"
#)
add_executable (eBench WIN32
	bwt.c
	console.c
	custom4.c
	ebench.c
	graphics.c
	image.c
	lodepng.c
	ols4.c
	ols5.c
	ols6.c
	ols7.c
	palette.c
	#transform_c20.c
	transforms.c
	#transforms_avx2.c
	#transforms_gpu.c
	#transforms_learned.c
	#transforms_pt2.c
	util.c
	wg4.c
	wg5.c
	window.c
)

#if(MSVC)
#	target_compile_options(eBench PRIVATE /W4)
#else()
#	target_compile_options(eBench PRIVATE -Wall -Wextra -Wpedantic)
#endif()

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET eBench PROPERTY CXX_STANDARD 20)
#endif()

