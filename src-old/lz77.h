#ifndef LZ77_H
#define LZ77_H
#include<string>

void			lz77_encode(const void *src, int imsize, std::string &dst, bool loud=false);
void			lz77_decode(const void *src, int csize, int imsize, char *dst, bool loud=false);

#endif