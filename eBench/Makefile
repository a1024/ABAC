SRC = console.c custom4.c ebench.c graphics.c image.c lodepng.c ols4.c ols5.c transforms.c util.c window.c
WARN = -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing -Wpointer-arith -Winit-self -Wshadow -Wswitch-enum -Wstrict-prototypes -Wmissing-prototypes -Wredundant-decls -Wformat=2 -Wfloat-equal -Wundef -Wvla -Wdeclaration-after-statement -Wc++-compat -Wno-cast-function-type -Wno-redundant-decls -Wno-switch-enum -Wno-unused-function
build:
	gcc -O3 $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
debug:
	gcc -g $(WARN) -mavx2 -mbmi -mlzcnt $(SRC) -o ebench-debug -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
