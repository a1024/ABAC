SRC =\
 bt2.c\
 bwt.c\
 console.c\
 custom4.c\
 c18.c\
 ebench.c\
 graphics.c\
 hadamard.c\
 image.c\
 lodepng.c\
 lpc.c\
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
 window.c\


WARN =\
 -Wall\
 -Wextra\
 -Wpedantic\
 -Wcast-qual\
 -Wcast-align\
 -Wstrict-aliasing\
 -Wpointer-arith\
 -Winit-self\
 -Wshadow\
 -Wswitch-enum\
 -Wstrict-prototypes\
 -Wmissing-prototypes\
 -Wredundant-decls\
 -Wformat=2\
 -Wfloat-equal\
 -Wundef\
 -Wvla\
 -Wdeclaration-after-statement\
 -Wno-cast-function-type\
 -Wno-redundant-decls\
 -Wno-switch-enum\
 -Wno-unused-function\
 -Wno-implicit-fallthrough\


# -std=c89
# -Wc++-compat


OPT =\
 -ffast-math\
 -mavx2\
 -mbmi\
 -mlzcnt\
 -mwindows\
 -static\


LIB =\
 -lcomdlg32\
 -lgdi32\
 -lole32\
 -lopengl32\
 -luuid\


ebench:
	gcc $(WARN) -O3 -s $(OPT) $(SRC) -o ebench-gcc $(LIB)
debug:
	gcc $(WARN) -g $(OPT) $(SRC) -o ebench-debug $(LIB)
