ABAC: An Adaptive Binary Arithmetic Coder


Contents:
ac.cpp: The arithmetic coder. Static and adaptive versions.
huffman.cpp: Huffman coder for comparison (uses vector_bool.h)


How it works:
Each bit-plane is compressed separately.

The probability of next bit is calculated as follows:
	P(next bit=0) = W(0)*confidence + (1/2)*(100%-confidence)

where W(0) is the weighted sum:
	W(0) = sum i=1 to LOG_WINDOW_SIZE: (1-bit[-i])*2^-i
which is the complement of the LOG_WINDOW_SIZE (default is 16) previous
bits loaded into an integer such that the MSB is the most recent bit.

while confidence is a function of the current compression ratio:
	ratio = consumed data size / current output size
	confidence = ratio/(ratio+1)

The arithmetic coder itself is a modified version of the coder from zpaq 1.10
See:
http://mattmahoney.net/dc/dce.html#Section_32


Evaluation:
Compression ratio for data without any transformations with the adaptive coder is as follows:

	~3.6	for a synthetic screenshot
		encode:  594M cycles
		decode: 1124M cycles

	~1.46	for a natural image

Use the functions abac_encode/abac_decode from ac.cpp.


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
