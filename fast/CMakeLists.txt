﻿add_executable (fast
#	f01.c
	f02.c
#	f03.c
#	f04.c
#	f05.c
#	f06.c
#	f07.c
	f08.c
#	f09.c
#	f10.c
#	f11.c
	f12.c
#	f13.c
#	f14.c
#	f15.c
	f16.c
#	f17.c
	f18.c
	f19.c
	f20.c
	f21.c
	f22.c
	f23.c
	f24.c
	f25.c
	f26.c
	f27.c
	f28.c
	f29.c
	f30.c
	f31.c
	f32.c
	f33.c
	fast.c
	format.c
#	huffman.c
	load_dng.cpp
	lodepng.c
	qoi.c
	transforms.c
	util.c
)

if(MSVC)
	message("Hello MSVC")
#	target_compile_options(fast PRIVATE /Wall)
else()
	message("Hello GCC")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mavx2 -mlzcnt")
#	message("CMAKE_C_FLAGS = " ${CMAKE_C_FLAGS})
#	target_compile_options(fast PRIVATE -Wall -Wextra -Wpedantic)
endif()

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET fast PROPERTY CXX_STANDARD 20)
#endif()
