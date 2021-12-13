#ifndef ABAC_H
#define ABAC_H
#include<stdlib.h>
#include<string>

#if defined __GNUC__ ||1
#define			NO_SIMD
#endif

void			ac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, bool loud=false);
void			ac_decode(const char *data, const int *sizes, const int *probabilities, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);

int				abac_estimate(const void *src, int imsize, int bitdepth, int bytestride, bool loud=false, int *sizes=nullptr);

#ifndef NO_SIMD
void			abac_encode_sse2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_sse2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode_avx2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_avx2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);
#endif

void			abac2_encode(const void *src, int imsize, int depth, int bytestride, std::string &out_data, int *out_sizes, int *out_conf, bool loud=false);
void			abac2_decode(const char *data, const int *sizes, const int *conf, void *dst, int imsize, int depth, int bytestride, bool loud=false);

void			abac3_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_conf, bool loud);
void			abac3_decode(const char *data, const int *sizes, const int *conf, short *buffer, int imsize, int depth, bool loud);


void			ac_test_bitplane_differentiation(short *buffer, int imsize, int depth, int &dmask);
void			ac_differentiate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_integrate_bitplanes(short *buffer, int imsize, int depth, int dmask);

#endif