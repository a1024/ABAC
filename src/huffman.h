#ifndef HUFFMAN_H
#define HUFFMAN_H
#include"vector_bool.h"

//interface:
void			huff_encode(const short *buffer, int imsize, int depth, int *histogram, vector_bool &bits, bool loud=false);
void			huff_decode(const int *src, long long bitsize, int imsize, int depth, const int *histogram, short *dst, bool loud=false);

//details:
typedef unsigned char byte;
void 			print(const char *format, ...);
void 			print_flush();
void			print_bin(const byte *data, int bytesize);
void			print_histogram(int *histogram, int nlevels, int scanned_count, int *sort_idx);
void			calculate_histogram(const short *image, int size, int *histogram, int nLevels);

#endif