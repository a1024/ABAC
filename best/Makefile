SRC = best.c f02_bypass.c f23.c format.c lodepng.c t39.c t42.c t44_paq8pxd.c t45_calic.c t46_slic4.c t47_slic5.c t54.c transforms.c util.c

build:
	gcc -Wall -O3 -mavx2 -mbmi -mlzcnt $(SRC) -o best -static
