#pragma once
#ifndef INC_AC_SIMD_H
#define INC_AC_SIMD_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#include<immintrin.h>
#ifdef __cplusplus
extern "C"
{
#endif


#ifdef AC_VALIDATE
static int abac_lane=0;
#endif


#define CMP_LE_U16(X, Y) _mm256_cmpeq_epi16(X, _mm256_min_epu16(X, Y))
#define CMP_GE_U16(X, Y) _mm256_cmpeq_epi16(X, _mm256_max_epu16(X, Y))
inline __m256i abacsimd_helper_condition(__m256i const *lo, __m256i const *hi)
{
	__m256i mone=_mm256_set1_epi32(-1), three=_mm256_set1_epi16(3);
	__m256i mid=_mm256_srli_epi16(_mm256_sub_epi16(*hi, *lo), 8);
	__m256i lo3=_mm256_add_epi16(*lo, three);
	__m256i cmp=CMP_GE_U16(*hi, lo3);
	cmp=_mm256_and_si256(_mm256_xor_si256(CMP_GE_U16(*lo, lo3), mone), cmp);
	cmp=_mm256_and_si256(_mm256_xor_si256(CMP_LE_U16(mid, _mm256_setzero_si256()), mone), cmp);
	return cmp;
}
typedef struct ABACSIMDEncStruct
{
	__m256i lo, hi, cache;
	char nbits[16];
	DList *lists;
} ABACSIMDEnc;
inline void abacsimd_enc_init(ABACSIMDEnc *ctx, DList *lists)
{
	ctx->lo=_mm256_setzero_si256();
	ctx->hi=_mm256_set1_epi32(-1);
	ctx->cache=_mm256_setzero_si256();
	memset(ctx->nbits, 0, 16);
	ctx->lists=lists;
}
inline void abacsimd_enc(ABACSIMDEnc *ctx, __m256i const *p0, __m256i const *bits)
{
	__m256i mid;
	__m256i mone=_mm256_set1_epi32(-1);

	for(;;)//renorm
	{
		__m256i cmp=abacsimd_helper_condition(&ctx->lo, &ctx->hi);
		int result=_mm256_movemask_epi8(cmp);
		if(result==-1)
			break;
		for(int k=0;k<16;++k)
		{
			if(!cmp.m256i_u16[k])//GCC won't like this
			{
				if(ctx->nbits[k]>=(sizeof(short)<<3))
				{
					dlist_push_back(ctx->lists+k, &ctx->cache.m256i_u16[k], 2);
					ctx->cache.m256i_u16[k]=0;
					ctx->nbits[k]=0;
				}
				ctx->cache.m256i_u16[k]|=(ctx->lo.m256i_u16[k]&0x8000)>>ctx->nbits[k];//cache is written MSB -> LSB
				++ctx->nbits[k];

				ctx->lo.m256i_u16[k]<<=1;//shift out MSB
				ctx->hi.m256i_u16[k]<<=1;
				ctx->hi.m256i_u16[k]|=1;
			}
		}
	}

	__m256i prob0=_mm256_slli_epi16(*p0, 8);
	if(_mm256_movemask_epi8(_mm256_cmpeq_epi16(prob0, _mm256_setzero_si256())))//reject degenerate distribution
		LOG_ERROR2("ZPS");

	mid=_mm256_sub_epi16(ctx->hi, ctx->lo);
	mid=_mm256_mulhi_epu16(mid, prob0);
	mid=_mm256_add_epi16(mid, ctx->lo);

#ifdef AC_VALIDATE
	int bit=bits->m256i_u16[abac_lane]&1;
	int prob=prob0.m256i_u16[abac_lane];
	int lo=ctx->lo.m256i_u16[abac_lane];
	int hi=ctx->hi.m256i_u16[abac_lane];
	int midk=mid.m256i_u16[abac_lane];
	acval_enc(bit, bit?prob:0, bit?0x10000-prob:prob, lo, hi, bit?midk:lo, bit?hi:midk-1, ctx->cache.m256i_u16[abac_lane], ctx->nbits[abac_lane]);
#endif

	__m256i bitmask=_mm256_cmpgt_epi16(*bits, _mm256_setzero_si256());
	ctx->lo=_mm256_blendv_epi8(ctx->lo, mid, bitmask);
	__m256i cbits=_mm256_xor_si256(bitmask, mone);
	mid=_mm256_add_epi16(mid, cbits);//OBLIGATORY range leak guard
	ctx->hi=_mm256_blendv_epi8(ctx->hi, mid, cbits);
}
inline void abacsimd_enc_flush(ABACSIMDEnc *ctx)
{
	for(int k=0;k<16;++k)
	{
		int k2=0;
		do//flush
		{
			while(ctx->nbits[k]<(sizeof(short)<<3))
			{
				ctx->cache.m256i_u16[k]|=(ctx->lo.m256i_u16[k]&0x8000)>>ctx->nbits[k];//cache is written MSB -> LSB
				++ctx->nbits[k];
				++k2;

				ctx->lo.m256i_u16[k]<<=1;//shift out MSB
				ctx->hi.m256i_u16[k]<<=1;
				ctx->hi.m256i_u16[k]|=1;
			}
			dlist_push_back(ctx->lists+k, &ctx->cache.m256i_u16[k], 2);
			ctx->cache.m256i_u16[k]=0;
			ctx->nbits[k]=0;
		}while(k2<(sizeof(short)<<3));
	}
}
typedef struct ABACSIMDDecStruct
{
	__m256i lo, hi, cache, code;
	char nbits[16];
	const unsigned char *srcptr[16], *srcend[16];
} ABACSIMDDec;
inline void abacsimd_dec_init(ABACSIMDDec *ctx, const unsigned char **starts, const unsigned char **ends)
{
	ctx->lo=_mm256_setzero_si256();
	ctx->hi=_mm256_set1_epi32(-1);
	memset(ctx->nbits, (sizeof(short)<<3), 16);
	//memcpy((void*)ctx->srcptr, checkpoints, 16*sizeof(void*));
	//memcpy((void*)ctx->srcend, checkpoints+1, 16*sizeof(void*));
	memcpy((void*)ctx->srcptr, starts, 16*sizeof(void*));
	memcpy((void*)ctx->srcend, ends, 16*sizeof(void*));

	for(int kl=0;kl<16;++kl)
	{
		if(ctx->srcend[kl]-ctx->srcptr[kl]<sizeof(short))
			LOG_ERROR2("buffer overflow");
		memcpy(&ctx->code.m256i_u16[kl], ctx->srcptr[kl], sizeof(short));
		ctx->srcptr[kl]+=sizeof(short);

		if(ctx->srcend[kl]-ctx->srcptr[kl]<sizeof(short))
			LOG_ERROR2("buffer overflow");
		memcpy(&ctx->cache.m256i_u16[kl], ctx->srcptr[kl], sizeof(short));
		ctx->srcptr[kl]+=sizeof(short);
	}
}
inline void abacsimd_dec(ABACSIMDDec *ctx, __m256i const *p0, __m256i *bits)
{
	__m256i mone=_mm256_set1_epi32(-1);
	__m256i mid;
	for(;;)//renorm
	{
		__m256i cmp=abacsimd_helper_condition(&ctx->lo, &ctx->hi);
		int result=_mm256_movemask_epi8(cmp);
		if(result==-1)
			break;
		for(int kl=0;kl<16;++kl)
		{
			if(!cmp.m256i_u16[kl])
			{
				if(!ctx->nbits[kl])
				{
					if(ctx->srcend[kl]-ctx->srcptr[kl]<sizeof(short))
					{
	#ifdef AC_VALIDATE
						printf("buffer overflow\n");
						acval_dump();
	#endif
						LOG_ERROR2("buffer overflow");
					}
					memcpy(&ctx->cache.m256i_u16[kl], ctx->srcptr[kl], sizeof(short));
					ctx->srcptr[kl]+=sizeof(short);

					ctx->nbits[kl]=(sizeof(short)<<3);
				}
				--ctx->nbits[kl];
				ctx->code.m256i_u16[kl]<<=1;//shift out MSB		cache is read MSB -> LSB
				ctx->code.m256i_u16[kl]|=(unsigned)(ctx->cache.m256i_u16[kl]>>ctx->nbits[kl]&1);

				ctx->lo.m256i_u16[kl]<<=1;
				ctx->hi.m256i_u16[kl]<<=1;
				ctx->hi.m256i_u16[kl]|=1;
			}
		}
	}

	__m256i prob0=_mm256_slli_epi16(*p0, 8);
	if(_mm256_movemask_epi8(_mm256_cmpeq_epi16(prob0, _mm256_setzero_si256())))//reject degenerate distribution
		LOG_ERROR2("ZPS");
	
	mid=_mm256_sub_epi16(ctx->hi, ctx->lo);
	mid=_mm256_mulhi_epu16(mid, prob0);
	mid=_mm256_add_epi16(mid, ctx->lo);

	*bits=CMP_GE_U16(ctx->code, mid);
	
#ifdef AC_VALIDATE
	int bit=bits->m256i_u16[abac_lane]&1;
	int prob=prob0.m256i_u16[abac_lane];
	int lo=ctx->lo.m256i_u16[abac_lane];
	int hi=ctx->hi.m256i_u16[abac_lane];
	int midk=mid.m256i_u16[abac_lane];
	acval_dec(bit, bit?prob:0, bit?0x10000-prob:prob, lo, hi, bit?midk:lo, bit?hi:midk-1, ctx->cache.m256i_u16[abac_lane], ctx->nbits[abac_lane], ctx->code.m256i_u16[abac_lane]);
#endif
	
	__m256i cbits=_mm256_xor_si256(*bits, mone);
	ctx->lo=_mm256_blendv_epi8(ctx->lo, mid, *bits);
	mid=_mm256_add_epi16(mid, cbits);//OBLIGATORY range leak guard
	ctx->hi=_mm256_blendv_epi8(ctx->hi, mid, cbits);

	*bits=_mm256_and_si256(*bits, _mm256_set1_epi16(1));
}


#ifdef __cplusplus
}
#endif
#endif
