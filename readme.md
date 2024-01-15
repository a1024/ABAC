# What this repository contains

These are my attempts and experiments to better understand image compression.
I didn't invent the algorithms used here, but tried to improve the implementations.

Each folder here contains a different project.
The most up-to-date projects are:

### `best`
The export folder of achievements.
Reference implementations of the best lossless algorithms that I comprehended.
This program only prints the efficiency and timing info.
But it doesn't save anything for security reasons.

### `e2`
Experiments on lossless compression.
Like `best` but filled with other obsolete algorithms.

### `eBench`
A GUI tool that evaluates reversible decorrelating transforms.
Supports high bit depth images.


## Archived Projects (Obsolete)

### `entropy-battle`
Experiments on lossless image compression. Older, messier project, contains even more obsolete techniques.

### `pxView3D`
A GUI visualization and benchmarking tool for reversible image transforms.
Inherently limited to 8-bit images, and uses modular arithmetic, which gives inferior results.

### `lossy`
Experiments on lossy compression (inactive).

### `src` and `src-v2`
Outdated. Only kept in case I need anything from there.


# Building
Currently, building these projects in the current state is supported only on Windows.
But the main coding algorithms can be easily ported to any platform.
These projects should compile fine on MSVC 2022.
Not tested yet on GCC.

### `best`
Either create an MSVC 2022 CMake project, or:

`gcc -O -mavx2 best.c t39.c t42.c t45_calic.c transforms.c util.c -o best`

Note that it fails to recover the original image with GCC for unknown reasons [TODO].

### `e2`
Either create an MSVC 2022 CMake project, or:

`gcc -O -mavx2 e2.c format.c lodepng.c tests.c tests2.c tests3.c tests4.c transforms.c util.c -o e2`

Note that it fails to recover the original image with GCC for unknown reasons [TODO].

### `eBench`
Either create an MSVC 2022 CMake project, or:
`gcc -O -mavx2 console.c ebench.c graphics.c image.c lodepng.c transforms.c util.c window.c -o ebench-gcc -lopengl32 -lgdi32 -lole32 -lcomdlg32 -luuid`
