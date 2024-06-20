SRC = e2.c format.c lodepng.c t51.c t54.c t55.c tests4.c tests5.c tests6.c transforms.c util.c

LINKFLAGS =
ifeq ($(OS),Windows_NT)
LINKFLAGS += -static
else
LINKFLAGS += -lm
endif
build:
	gcc -Wall -O3 -mavx2 -mbmi -mlzcnt $(SRC) -o e2 $(LINKFLAGS)
