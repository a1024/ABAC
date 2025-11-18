SRC =\
	bwt.c\
	console.c\
	custom4.c\
	c18.c\
	ebench.c\
	graphics.c\
	hadamard.c\
	image.c\
	lodepng.c\
	ols4.c\
	ols5.c\
	ols6.c\
	ols7.c\
	ols8.c\
	ols9.c\
	palette.c\
	ssim.c\
	transforms.c\
	util.c\
	wg4.c\
	wg5.c\
	wg6.c\
	wg7.c\
	window.c

WARN = -Wall
#WARN = -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing -Wpointer-arith -Winit-self -Wshadow -Wswitch-enum -Wstrict-prototypes -Wmissing-prototypes -Wredundant-decls -Wformat=2 -Wfloat-equal -Wundef -Wvla -Wdeclaration-after-statement -Wc++-compat -Wno-cast-function-type -Wno-redundant-decls -Wno-switch-enum -Wno-unused-function -Wno-implicit-fallthrough
eBench:
	gcc -O3 $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid -static -s
debug:
	gcc -g $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o ebench-debug -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid -static
