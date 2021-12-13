ABAC: An Adaptive Binary Arithmetic Coder


Documentation:
ac.cpp: The arithmetic coder. Static and adaptive versions.
	abac_encode/decode: adaptive binary arithmetic coder
	abac_[encode/decode]_[sse2/avx2]: faster SIMD versions,
		incompatible with non-SIMD versions.
	ac_encode/decode: static probability binary arithmetic coder

huffman.cpp: Huffman coder for comparison (uses vector_bool.h)
	huff_encode/decode


How it works:
Each bit-plane is compressed separately.

The probability of next bit is calculated as follows:
	P(next bit=0) = W(0)*confidence + (1/2)*(100%-confidence)
		= 1/2 + (W(0)-1/2)*confidence

where W(0) is the weighted sum:
	W(0) = sum i=1 to LOG_WINDOW_SIZE: (1-bit[-i])*2^-i
which is the complement of the LOG_WINDOW_SIZE (default is 16) previous
bits loaded into an integer such that the MSB is the most recent bit.

while confidence can be a function of the current compression ratio:
	ratio = consumed data size / current output size
	confidence v1 = ratio/(ratio+1)
		= consumed bit count / (produced + consumed bit counts)

v1.5:
Higher compression ratio is achieved with the confidence defined as the
ratio of correctly predicted bits.
	confidence v2 = number of correctly predicted bits
		/ number of data bits consumed

v2:
An even higher compression ratio was achieved.
See abac2_encode/decode in ac.cpp.


The arithmetic coder itself is a modified version of the coder from zpaq 1.10
See:
http://mattmahoney.net/dc/dce.html#Section_32


Evaluation:
Data was compressed without any decorrelating transformations.
SIMD ABAC versions are currently incompatible with normal ABAC,
due to floating point rounding behavior.


Coder		Compression	Encode		Decode
		ratio		cycles		cycles

1) synthetic 1920x1080 image:

ZPAQ0+RVL	5.075552	1011M		1339M
ABAC v2+RVL	6.237550	1094M		 878M
ABAC v2		6.006831	 773M		 876M
ABAC v1.5	5.206832	1176M		1763M
ABAC v1 AVX2	3.608005	 287M		 215M
ABAC v1 SSE4	3.608005	 254M		 352M
ABAC v1		3.607679	 687M		1477M
ZPAQ0 AC	1.379686	1293M		1342M
Huffman		2.41		  92.7M		  49.3M
Simple LZ77	2.006596	 173M		   4.8M
Static prob AC	1.28		 352M		 420M

2) natural 3456x2304 image:

ZPAQ0+RVL	2.235712	4236M		4953M
ABAC v2+RVL	2.161590	2964M		2959M
ABAC v2		1.613031	3555M		2612M
ABAC v1.5	1.508929	4815M		6692M
ABAC v1 AVX2	1.471159	1056M		1070M
ABAC v1 SSE4	1.471159	1347M		1208M
ABAC v1		1.471149	3254M		5400M
ZPAQ0 AC	1.054280	4944M		5468M
Huffman		1.16		 466M		 499M
Simple LZ77	1.003379	1008M		   5.37M
Static prob AC	1.01		2219M		2127M

3) C++ source 1.32 MB:

ABAC v2		1.211128	 710M		 446M
ABAC v1 AVX2	1.223053	 321M		 267M
ABAC v1 SSE4	1.223053	 275M		 285M
ABAC v1		1.223054	 586M		 977M
Huffman		1.43		 132M		  99.7M
Simple LZ77	1.211728	 102M		   3.39M
Static prob AC	1.199		 359M		 437M
LZMA2 (7-zip)	8.5033		~805M		?


Build test:
1) Add the following files to the src directory:
	- stb_image.h
	- lodepng.h
	- lodepng.cpp

lodepng by Lode Vandevenne
https://github.com/lvandeve/lodepng

stb_image by Sean Barrett
https://github.com/nothings/stb/blob/master/stb_image.h

2) On Windows make a Visual Studio project (compiled successfully on MSVC 2013).
   On Linux type 'make release' in terminal.

TODO:
- Port SIMD code to Linux.
