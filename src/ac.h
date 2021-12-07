#ifndef ABAC_H
#define ABAC_H
#include<stdlib.h>
#include<string>

void			ac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, bool loud=false);
void			ac_decode(const char *data, const int *sizes, const int *probabilities, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);

#ifndef __GNUC__
void			abac_encode_sse2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_sse2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode_avx2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_avx2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);
#endif

#if 0
typedef unsigned long long u64;
struct			ABAC_Plane
{
	unsigned start, end;
	int prob, hitcount;
	union
	{
		std::string *str;
		unsigned char *code;
	};
	ABAC_Plane():start(0), end(0xFFFFFFFF), prob(0x8000), hitcount(1){}
};
struct			ABAC_Context
{
	int depth;
	int orig_size, kb;
	ABAC_Plane *planes;
};
void			abac2_encode_init(ABAC_Context &context, int depth, int imsize=0);
void			abac2_encode_finish(ABAC_Context &context, std::string &data, int *out_sizes);
void			abac2_encode(ABAC_Context &context, u64 pixel);

void			abac2_decode_init(ABAC_Context &context, std::string const &data, const int *out_sizes);
u64				abac2_decode(ABAC_Context &context);
#endif

int				abac_estimate(const void *src, int imsize, int bitdepth, int bytestride, bool loud=false, int *sizes=nullptr);


void			ac_test_bitplane_differentiation(short *buffer, int imsize, int depth, int &dmask);
void			ac_differentiate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_integrate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_debug(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, short *out, bool loud);

#endif