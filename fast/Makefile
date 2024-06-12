SRC = f02.c f08.c f20.c f23.c f24.c f25.c f26.c f27.c fast.c format.c lodepng.c qoi.c transforms.c util.c load_dng.cpp
#SRC = f02.c f08.c f12.c f16.c f18.c f19.c f20.c f23.c f24.c fast.c format.c lodepng.c qoi.c transforms.c util.c load_dng.cpp
#SRC = f01.c f02.c f03.c f04.c f05.c f06.c f07.c f08.c f09.c f10.c f11.c f12.c f13.c f14.c f15.c f16.c f17.c f18.c f19.c f20.c fast.c format.c huffman.c lodepng.c transforms.c util.c load_dng.cpp

WARN = -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing -Wpointer-arith -Winit-self -Wshadow\
	-Wstrict-prototypes -Wmissing-prototypes -Wredundant-decls -Wformat=2 -Wfloat-equal -Wundef\
	-Wvla -Wdeclaration-after-statement -Wc++-compat -Winline\
	-Wno-cast-function-type -Wno-redundant-decls -Wno-unused-function -Wno-implicit-fallthrough -Wno-misleading-indentation -Wno-missing-field-initializers

WARN_CLANG = -Wno-cast-align

build:
	g++ -O3 $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o fast
debug:
	g++ -g $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o fast-debug
avx512:
	g++ -O3 $(WARN) -mavx512f -mavx512cd -mbmi -mlzcnt $(SRC) -o fast

use-clang-avx512:
	clang++ -O3 $(WARN) $(WARN_CLANG) -mavx512f -mavx512cd -mbmi -mlzcnt $(SRC) -o fast
