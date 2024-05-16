SRC = console.c custom4.c ebench.c graphics.c image.c lodepng.c ols4.c ols5.c transforms.c util.c window.c
build:
	gcc -Wall -mavx2 -mbmi -mlzcnt -O3 $(SRC) -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
debug:
	gcc -Wall -mavx2 -mbmi -mlzcnt -g $(SRC) -o ebench-debug -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid
