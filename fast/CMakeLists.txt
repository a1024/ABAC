add_executable (fast
	f01.c
	f02.c
	f03.c
	f04.c
	fast.c
	format.c
	huffman.c
	lodepng.c
	transforms.c
	util.c
)

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET fast PROPERTY CXX_STANDARD 20)
#endif()
