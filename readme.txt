ABAC: An Adaptive Binary Arithmetic Coder


Documentation:
ac.cpp: The arithmetic coder. Static and adaptive versions.
	abac_encode/decode: adaptive binary arithmetic coder
	ac_encode/decode: static probability binary arithmetic coder

huffman.cpp: Huffman coder for comparison (uses vector_bool.h)
	huff_encode/decode


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

Image		Method		compression	encode		decode
				ratio		cycles		cycles


raw synthetic	ABAC		3.61		 609M		1126M
1920x1080	Huffman		2.41		  92.7M		  49.3M
		Static Prob AC	1.28		 352M		 420M

raw natural	ABAC		1.47		2511M		4726M
3456x2304	Huffman		1.16		 466M		 499M
		Static prob AC	1.01		2219M		2127M


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
