# What this repository contains

These are my attempts and experiments to better understand image compression.
I didn't invent the algorithms used here, but tried to improve the implementations.

Each folder contains a different project.

## Maintained projects

### `best`
The export folder of achievements. Sometimes not up to date.
Reference implementations of the best open-source lossless algorithms that I implemented/encountered.

`gcc -O3 -mavx2  best.c format.c lodepng.c t39.c t42.c t45_calic.c t46_slic4.c t47_slic5.c transforms.c util.c  -o best`

### `e2`
Contains most of my latest experiments on lossless image compression, and other older/obsolete algorithms.

`gcc -O3 -mavx2 -mbmi  e2.c format.c lodepng.c t51.c t54.c t55.c tests4.c tests5.c tests6.c transforms.c util.c  -o e2`

### `eBench`
A GUI tool that evaluates reversible decorrelating transforms.
Supports high bit depth images.
Currently works only on Windows.

`gcc -O3 -mavx2  console.c ebench.c graphics.c image.c lodepng.c ols4.c transforms.c util.c window.c  -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid`

### `fast`
Experiments prioritizing speed.

`gcc -O3 -mavx2 -mbmi -mlzcnt  f02.c f08.c f16.c f18.c f19.c f20.c fast.c format.c lodepng.c transforms.c util.c load_dng.cpp  -o fast`

### `imgcvt`
Batch image converter.

`gcc -O3 imgcvt.c util.c -o imgcvt`


## Archived projects (obsolete)

### `entropy-battle`
Experiments on lossless image compression. Older, messier project, contains even more obsolete techniques.

### `pxView3D`
A GUI visualization and benchmarking tool for reversible image transforms.
Inherently limited to 8-bit images, and uses modular arithmetic, which gives inferior results.

### `lossy`
Experiments on lossy compression (inactive).

### `src` and `src-v2`
Outdated. Only kept in case I need anything from there.
