# What this repository contains

These are my attempts and experiments to better understand image compression.
I didn't invent the algorithms used here, but tried to improve the implementations.

Each folder contains a different project.

## Maintained projects

### `best`
The export folder of achievements. Sometimes not up to date.
Reference implementations of the best open-source lossless algorithms that I implemented/encountered.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`make`

### `e2`
Contains most of my latest experiments on lossless image compression, and other older/obsolete algorithms.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`gcc -O3 -mavx2 -mbmi  e2.c format.c lodepng.c t51.c t54.c t55.c tests4.c tests5.c tests6.c transforms.c util.c  -o e2`

### `eBench`
A GUI tool that evaluates reversible decorrelating transforms.
Supports high bit depth images.
Currently works only on Windows.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`make`

### `fast`
Experiments prioritizing speed.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`make`

### `imgcvt`
Batch image converter.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`gcc -O3 imgcvt.c util.c -o imgcvt`

### `ppmcodec`
Plain 8-bit codecs without alpha. PPM only.
To compile, either make an MSVC 2022 CMake project, or to use GCC:

`make`


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
