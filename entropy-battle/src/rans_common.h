#pragma once
#ifndef INC_RANS_COMMON_H
#define INC_RANS_COMMON_H

//do not change
#define ANS_DEPTH      8
#define ANS_PROB_BITS 16
#define ANS_NLEVELS  256	//(1<<ANS_DEPTH)
#define ANS_L      65536	//(1<<ANS_PROB_BITS)

int rans_calc_histogram(const unsigned char *buffer, int nsymbols, int bytestride, unsigned char *hist, int prob_bits, int integrate);//hist is unsigned char due to alignment issues, but it's 16bit

typedef struct SymbolInfoStruct//32 bytes, size must be a power of two
{
	unsigned short
		freq,
		neg_freq,
	//	shift,
		reserved0;
	unsigned
		CDF,
		bias,
		renorm_limit;
	unsigned long long inv_freq;
	union
	{
		long long invf;
		unsigned invfcomp[2];
	};
} SymbolInfo;
int rans_prep(const void *hist_ptr, int bytespersymbol, SymbolInfo **info, unsigned char **CDF2sym, int loud);


#endif//INC_RANS_COMMON_H
