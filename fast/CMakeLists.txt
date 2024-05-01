add_executable (fast
	f01.c
	f02.c
	f03.c
	f04.c
	f05.c
	f06.c
	f07.c
	f08.c
	f09.c
	f10.c
	f11.c
	f12.c
	fast.c
	format.c
	huffman.c
	load_dng.cpp
	lodepng.c
	transforms.c
	util.c
)
#if(MSVC)
#	target_compile_options(fast PRIVATE /W4)
#else()
#	target_compile_options(fast PRIVATE -Wall -Wextra -Wpedantic)
#endif()

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET fast PROPERTY CXX_STANDARD 20)
#endif()
