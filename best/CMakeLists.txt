# CMakeList.txt : CMake project for best, include source and define project specific logic here.

add_executable(best
	best.c
	f02_bypass.c
	f23.c
	format.c
	lodepng.c
	t39.c
	t42.c
	t44_paq8pxd.c
	t45_calic.c
	t46_slic4.c
	t47_slic5.c
	t54.c
	transforms.c
	util.c
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET best PROPERTY CXX_STANDARD 20)
endif()
