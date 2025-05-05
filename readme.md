# What this repository contains

These are my attempts and experiments to better understand image compression.
I didn't invent the algorithms used here, but tried to improve the implementations.

Each folder contains a different project.

## Maintained projects

### `cbm`
A codec benchmark that tracks history per-file.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`make`

### `eBench`
A GUI tool that evaluates reversible decorrelating transforms.
Supports high bit depth images.
Currently works only on Windows.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`make`

### `ppmcodec`
Plain 8-bit codecs without alpha. PPM only.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`make c<codecID> [PROF=1]`


## Archived projects (obsolete)

### `best`
[outdated]
The best codecs are placed here. Sometimes not up to date.
Reference implementations of the best open-source lossless algorithms that I implemented/encountered.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`make`

### `e2`
[outdated]
Contains most of my latest experiments on lossless image compression, and other older/obsolete algorithms.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`gcc -O3 -mavx2 -mbmi  e2.c format.c lodepng.c t51.c t54.c t55.c tests4.c tests5.c tests6.c transforms.c util.c  -o e2`

### `entropy-battle`
Experiments on lossless image compression. Older, messier project, contains even more obsolete techniques.

### `fast`
[outdated]
Experiments prioritizing speed.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`make`

### `imgcvt`
[outdated]
Batch image converter.
To build, either create an MSVC 2022 CMake project, or to use GCC:

`gcc -O3 imgcvt.c util.c -o imgcvt`

### `lossy`
Experiments on lossy compression (inactive).

### `pxView3D`
A GUI visualization and benchmarking tool for reversible image transforms.
Inherently limited to 8-bit images, and uses modular arithmetic, which gives inferior results.

### `src` and `src-v2`
Outdated. Only kept in case I need anything from there.
