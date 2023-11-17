# CMakeList.txt : CMake project for pxView3D, include source and define project specific logic here.

include_directories(
	"C:/Projects/a/OpenCL-Headers"
)
link_directories(
	"C:/Projects/a/OpenCL-Headers"
)
add_executable (pxView3D WIN32
	"console.c"
	"graphics.c"
	"lodepng.c"
	"pxview3d.c"
	"transforms.c"
	"transforms_avx2.c"
	"transforms_gpu.c"
	"transforms_learned.c"
	"transforms_pt2.c"
	"util.c"
	"window.c"
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET pxView3D PROPERTY CXX_STANDARD 20)
endif()

