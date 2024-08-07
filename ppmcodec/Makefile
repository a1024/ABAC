SRC = c03.c codec.c util.c
#SRC = c05.c codec.c util.c
#SRC = c01.c c02.c c03.c c04.c c05.c c06.c c07.c codec.c util.c

build:
	gcc -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o codec -static

build-clang:
	clang -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o codec -static
