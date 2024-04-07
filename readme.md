# What this repository contains

These are my attempts and experiments to better understand image compression.
I didn't invent the algorithms used here, but tried to improve the implementations.

Each folder contains a different project.

## Maintained projects

### `best`
The export folder of achievements. Sometimes not up to date.
Reference implementations of the best lossless algorithms that I encountered.

`gcc -Wall -O3 best.c format.c lodepng.c t39.c t42.c t45_calic.c t46_slic4.c t47_slic5.c transforms.c util.c -o best`

### `e2`
Experiments on lossless compression.
Contains my latest experiments, and other obsolete algorithms.

`gcc -Wall -O3 -mavx2 e2.c format.c lodepng.c t51.c t54.c t55.c tests4.c tests5.c tests6.c transforms.c util.c -o e2`

### `eBench`
A GUI tool that evaluates reversible decorrelating transforms.
Supports high bit depth images.
Currently works only on Windows.

`gcc -O3 -mavx2 console.c ebench.c graphics.c image.c lodepng.c transforms.c util.c window.c -o ebench -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid`

### `imgcvt`
Batch image converter.

`gcc -Wall -O3 imgcvt.c util.c -o imgcvt`


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
