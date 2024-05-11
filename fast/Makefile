SRC = f01.c f02.c f03.c f04.c f05.c f06.c f07.c f08.c f09.c f10.c f11.c f12.c f13.c f14.c f15.c f16.c f17.c f18.c f19.c f20.c fast.c format.c huffman.c lodepng.c transforms.c util.c  load_dng.cpp
build:
	g++ -Wall -mavx2 -mbmi -mlzcnt -O3 $(SRC) -o fast
debug:
	g++ -Wall -mavx2 -mbmi -mlzcnt -g $(SRC) -o fast-debug
