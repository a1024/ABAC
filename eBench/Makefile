SRC = console.c custom4.c ebench.c graphics.c image.c lodepng.c ols4.c ols5.c transforms.c util.c window.c
build:
	gcc -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
debug:
	gcc -g -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o ebench-debug -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
