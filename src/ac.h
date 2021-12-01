#ifndef ABAC_H
#define ABAC_H
#include<stdlib.h>
#include<string>

void			ac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, bool loud=false);
void			ac_decode(const char *data, const int *sizes, const int *probabilities, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud);
void			abac_decode(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud);

#ifndef __GNUC__
void			abac_encode_sse2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud);
void			abac_decode_sse2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud);

void			abac_encode_avx2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud);
void			abac_decode_avx2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud);
#endif

void			ac_test_bitplane_differentiation(short *buffer, int imsize, int depth, int &dmask);
void			ac_differentiate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_integrate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_debug(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, short *out, bool loud);

/*#define	SYMBOL_BITS		32

#if SYMBOL_BITS==64
typedef __uint128_t MulType;
typedef __uint64_t Symbol;
#elif SYMBOL_BITS==32
typedef __uint64_t MulType;
typedef __uint32_t Symbol;
#elif SYMBOL_BITS==16
typedef __uint32_t MulType;
typedef __uint16_t Symbol;
#else
typedef __uint16_t MulType;
typedef __uint8_t Symbol;
#endif

void			ac_encode(const short *buffer, int imsize, int depth, std::vector<int> &data, bool loud=false);
void			ac_decode(const int *data, short *buffer, int imsize, int depth, bool loud=false);//*/

#endif