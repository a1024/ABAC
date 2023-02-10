#include"util.h"
#include<stdio.h>//only for debugging
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

//#ifdef USE_SSE4_1
//#include<smmintrin.h>//SSE4.1
//#else
#include<tmmintrin.h>//SSSE3
//#endif

#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
static const char file[]=__FILE__;


//	#define PRINT_SSE2
	#define ENABLE_CONFIDENCE
//	#define ENABLE_GRAYCODE


//#ifdef USE_SSE4_1
//#define EXTRACT_EPU32(REG, IDX)		_mm_extract_epi32(REG, IDX)
//#else
//#define EXTRACT_EPU32(REG, IDX)		(_mm_extract_epi16(REG, (IDX)<<1|1)<<16|_mm_extract_epi16(REG, (IDX)<<1))
//#endif
static const int tag_ac0A='A'|'C'<<8|'0'<<16|'A'<<24;

#ifdef PRINT_SSE2
static void print_reg(__m128i const *r, int n, const char *format, ...)
{
	va_list args;

	if(n==16)
	{
		ALIGN(16) short mem[8];
		_mm_store_si128((__m128i*)mem, r[0]);
		for(int k=7;k>=0;--k)
			printf("%04X ", mem[k]);
		//for(int k=0;k<8;++k)
		//	printf("%04X ", r->m128i_i16[7-k]);
	}
	else if(n==32)
	{
		ALIGN(16) int mem[8];
		_mm_store_si128((__m128i*)mem, r[0]);
		_mm_store_si128((__m128i*)mem+1, r[1]);
		for(int k=7;k>=0;--k)
			printf("%08X ", mem[k]);
		//for(int k2=1;k2>=0;--k2)
		//	for(int k=3;k>=0;--k)
		//		printf("%08X ", r[k2].m128i_i32[k]);
	}
	if(format)
	{
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
}
static void print_prob(__m128i *reg, int idx, int w, const char *format, ...)
{
	va_list args;
	ALIGN(16) short mem[8];
	_mm_store_si128((__m128i*)mem, *reg);

	for(int k=0;k<w;++k)
		printf("%c", '1'-(k*0x7FFF<mem[idx]*w));
	if(format)
	{
		printf(" ");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
}
#else
#define print_reg(...)
#define print_prob(...)
#endif

//abac0a_encode(): Encodes 8-bit symbols. Uses SSSE3 with 15-bit probability.
//If output was initialized, output->esize must be 1.
int abac0a_encode(const unsigned char *src, int count, int bytestride, ArrayHandle *output, int loud)
{
	DList list[8];
	const unsigned char *buffer;
	unsigned char *header;
	size_t offset0;
	long long time1, time2;
	__m128i
		ones=_mm_set1_epi32(-1), ones32=_mm_set1_epi32(1), three=_mm_set1_epi32(3),
		prob_min=_mm_set1_epi16(1), prob_half=_mm_set1_epi16(0x4000), prob_max=_mm_set1_epi16(0x7FFE),
		prob=prob_half,
#ifdef ENABLE_CONFIDENCE
		prob_correct=prob_half,
#endif
		limit=_mm_set1_epi32(0x1000000);
	__m128i shuf_pack=_mm_set_epi8(15, 14, 11, 10, 7, 6, 3, 2, 13, 12, 9, 8, 5, 4, 1, 0);
	__m128i start[2], range[2], r2[2];
	ALIGN(16) unsigned a_start[8], a_range[8];

	if(!src||!count||!bytestride)
		LOG_ERROR("abac4_encode(src=%p, count=%d, stride=%d)", src, count, bytestride);

	time1=__rdtsc();
	if(*output)
	{
		if(output[0]->esize!=1)
			LOG_ERROR("Output array element size must be 1 byte");
		offset0=output[0]->count;
		ARRAY_APPEND(*output, 0, (1+8)*sizeof(int), 1, 0);
	}
	else
	{
		offset0=0;
		ARRAY_ALLOC(char, *output, 0, (1+8)*sizeof(int), 0, 0);//magic + sizes[bitdepth]
	}

	header=(unsigned char*)array_back(output)+1-(1+8)*sizeof(int);
	memcpy(header, &tag_ac0A, sizeof(int));
	header+=sizeof(int);

	for(int k=0;k<8;++k)
		dlist_init(list+k, 1, 1024, 0);

	buffer=(const unsigned char*)src;
	start[0]=_mm_setzero_si128();
	start[1]=_mm_setzero_si128();
	range[0]=_mm_set1_epi32(0x7FFFFFFF);//the MSB is the 2's comp sign, all SSE2 cmp instructions are signed
	range[1]=_mm_set1_epi32(0x7FFFFFFF);
	//unsigned start[8]={0}, r2[8], range[8];
	//for(int k=0;k<8;++k)
	//	range[k]=0xFFFFFFFF;
	for(int ks=0, ks2=0;ks<count;++ks, ks2+=bytestride)
	{
		print_reg(start, 32, "%d start", ks);
		print_reg(range, 32, "%d range", ks);

		//calculate probability of zero bit
#ifdef ENABLE_CONFIDENCE
		__m128i p0=_mm_sub_epi16(prob, prob_half);

		p0=_mm_slli_epi16(p0, 1);
		prob_correct=_mm_slli_epi16(prob_correct, 1);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_mulhi_epu16(p0, p0);
		p0=_mm_srli_epi16(p0, 1);
		prob_correct=_mm_srli_epi16(prob_correct, 1);

		p0=_mm_add_epi16(p0, prob_half);
#else
		__m128i p0=prob;
#endif

		//clamp probability			eg: {0, 1, half, FFFF}
		__m128i cmp=_mm_cmplt_epi16(prob_min, p0);//{0, FFFF, FFFF, FFFF}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		p0=_mm_or_si128(p0, cmp);

		cmp=_mm_cmplt_epi16(p0, prob_max);//{FFFF, FFFF, FFFF, 0}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_max);
		p0=_mm_or_si128(p0, cmp);
		
		//calculate middle (r2)
		__m128i t0=_mm_shuffle_epi8(range[0], shuf_pack);
		__m128i t1=_mm_shuffle_epi8(range[1], shuf_pack);
		__m128i r_lo=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(1, 0, 1, 0)));
		__m128i r_hi=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(3, 2, 3, 2)));
		p0=_mm_slli_epi16(p0, 1);
		__m128i a=_mm_mulhi_epu16(r_lo, p0);
		__m128i b0=_mm_mullo_epi16(r_hi, p0);
		__m128i b1=_mm_mulhi_epu16(r_hi, p0);
		p0=_mm_srli_epi16(p0, 1);

		r2[0]=_mm_unpacklo_epi16(b0, b1);
		r2[1]=_mm_unpackhi_epi16(b0, b1);

		__m128i d0=_mm_unpacklo_epi16(a, _mm_setzero_si128());
		__m128i d1=_mm_unpackhi_epi16(a, _mm_setzero_si128());

		r2[0]=_mm_add_epi32(r2[0], d0);
		r2[1]=_mm_add_epi32(r2[1], d1);

		//avoid degenerate range		r2[k]+=(r2[k]==0)-(r2[k]==range[k]);
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		__m128i cmp2=_mm_cmpeq_epi32(r2[0], range[0]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[0]=_mm_add_epi32(r2[0], cmp);
		
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		cmp2=_mm_cmpeq_epi32(r2[0], range[1]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[1]=_mm_add_epi32(r2[1], cmp);

		print_reg(r2, 32, "%d r2", ks);
		
		print_reg(&prob, 16, "%d prob ", ks);
		print_prob(&prob, 7, 64, "bit 7");
		print_reg(&p0, 16, "%d p0 ", ks);
		print_prob(&p0, 7, 64, "bit 7");


		unsigned char sym=src[ks2];
#ifdef ENABLE_GRAYCODE
		sym^=sym>>1;
#endif
		__m128i bit=_mm_set_epi16(sym>>7&1, sym>>6&1, sym>>5&1, sym>>4&1, sym>>3&1, sym>>2&1, sym>>1&1, sym>>0&1);
		
#ifdef ENABLE_CONFIDENCE
		//update prob_correct			correct=bit^(p0>=PROB_HALF), prob_correct=correct<<15|prob_correct>>1;
		cmp=_mm_cmplt_epi16(p0, prob_half);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		cmp=_mm_xor_si128(cmp, bit);

		cmp=_mm_slli_epi16(cmp, 14);//set pre-MSB, because setting MSB would negate the 2's comp value
		prob_correct=_mm_srli_epi16(prob_correct, 1);
		prob_correct=_mm_or_si128(prob_correct, cmp);
#endif

		//update prob zero				prob=!bit<<15|prob>>1;
		__m128i bit2=_mm_xor_si128(bit, prob_min);
		prob=_mm_srli_epi16(prob, 1);
		prob=_mm_or_si128(prob, _mm_slli_epi16(bit2, 14));
		
		print_reg(&prob, 16, "%d prob2 ", ks);
		print_reg(&bit, 16, "%d bit ", ks);


		//update range
		__m128i bitlo=_mm_unpacklo_epi16(bit, _mm_setzero_si128());
		__m128i bithi=_mm_unpackhi_epi16(bit, _mm_setzero_si128());

		//r2[k]+=bit.m128i_i16[k]-!bit.m128i_i16[k];
		r2[0]=_mm_add_epi32(r2[0], bitlo);
		r2[1]=_mm_add_epi32(r2[1], bithi);
		r2[0]=_mm_sub_epi32(r2[0], _mm_sub_epi32(ones32, bitlo));
		r2[1]=_mm_sub_epi32(r2[1], _mm_sub_epi32(ones32, bithi));

		//start[k]+=r2[k]&-bit.m128i_i16[k];
		bitlo=_mm_cmpeq_epi32(bitlo, ones32);
		bithi=_mm_cmpeq_epi32(bithi, ones32);
		t0=_mm_and_si128(bitlo, r2[0]);
		t1=_mm_and_si128(bithi, r2[1]);
		start[0]=_mm_add_epi32(start[0], t0);
		start[1]=_mm_add_epi32(start[1], t1);
		//__m128i start_lo=_mm_loadu_si128((const __m128i*)start);
		//__m128i start_hi=_mm_loadu_si128((const __m128i*)(start+4));

		//range[k]=((range[k]-r2[k])&-bit.m128i_i16[k])|(r2[k]&-!bit.m128i_i16[k]);
		range[0]=_mm_sub_epi32(range[0], r2[0]);
		range[1]=_mm_sub_epi32(range[1], r2[1]);
		range[0]=_mm_and_si128(range[0], bitlo);
		range[1]=_mm_and_si128(range[1], bithi);
		bitlo=_mm_xor_si128(bitlo, ones);
		bithi=_mm_xor_si128(bithi, ones);
		r2[0]=_mm_and_si128(r2[0], bitlo);
		r2[1]=_mm_and_si128(r2[1], bithi);
		range[0]=_mm_or_si128(range[0], r2[0]);
		range[1]=_mm_or_si128(range[1], r2[1]);

		//renormalize
		unsigned condition;
		for(;;)
		{
			t0=_mm_add_epi32(start[0], range[0]);
			t1=_mm_add_epi32(start[1], range[1]);
			t0=_mm_xor_si128(t0, start[0]);
			t1=_mm_xor_si128(t1, start[1]);
			t0=_mm_cmpeq_epi32(t0, limit);
			t1=_mm_cmpeq_epi32(t1, limit);
			condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
			if(!condition)
				break;
			_mm_store_si128((__m128i*)a_start  , start[0]);
			_mm_store_si128((__m128i*)a_start+1, start[1]);
			_mm_store_si128((__m128i*)a_range  , range[0]);
			_mm_store_si128((__m128i*)a_range+1, range[1]);
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					dlist_push_back1(list+k, (char*)&a_start[k]+3);
#ifdef PRINT_SSE2
					printf("RENORM bit %d <-[%02X]%06X\n", k, ((char*)&a_start[k])[3], a_start[k]&0xFFFFFF);
#endif
					a_start[k]<<=8;
					a_range[k]=a_range[k]<<8|0xFF;
				}
			}
			start[0]=_mm_load_si128((__m128i*)a_start);
			start[1]=_mm_load_si128((__m128i*)a_start+1);
			range[0]=_mm_load_si128((__m128i*)a_range);
			range[1]=_mm_load_si128((__m128i*)a_range+1);
		}
		t0=_mm_cmplt_epi32(range[0], three);
		t1=_mm_cmplt_epi32(range[1], three);
		condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
		if(condition)
		{
			_mm_store_si128((__m128i*)a_start  , start[0]);
			_mm_store_si128((__m128i*)a_start+1, start[1]);
			_mm_store_si128((__m128i*)a_range  , range[0]);
			_mm_store_si128((__m128i*)a_range+1, range[1]);
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					dlist_push_back1(list+k, (char*)&a_start[k]+3);//big endian
					dlist_push_back1(list+k, (char*)&a_start[k]+2);
					dlist_push_back1(list+k, (char*)&a_start[k]+1);
					dlist_push_back1(list+k, &a_start[k]);
#ifdef PRINT_SSE2
					printf("RENORM bit %d <-[%08X]\n", k, a_start[k]);
#endif
					a_start[k]=0;
					a_range[k]=0x7FFFFFFF;//because 1=0.9999...
				}
			}
			start[0]=_mm_load_si128((__m128i*)a_start);
			start[1]=_mm_load_si128((__m128i*)a_start+1);
			range[0]=_mm_load_si128((__m128i*)a_range);
			range[1]=_mm_load_si128((__m128i*)a_range+1);
		}
#ifdef PRINT_SSE2
		printf("\n");
#endif
	}
	_mm_store_si128((__m128i*)a_start  , start[0]);
	_mm_store_si128((__m128i*)a_start+1, start[1]);
	for(int k=0;k<8;++k)
	{
		dlist_push_back1(list+k, (char*)&a_start[k]+3);//big endian
		dlist_push_back1(list+k, (char*)&a_start[k]+2);
		dlist_push_back1(list+k, (char*)&a_start[k]+1);
		dlist_push_back1(list+k, &a_start[k]);
	}
	for(int k=0;k<8;++k)
		memcpy(header+k*sizeof(int), &list[k].nobj, sizeof(int));
	for(int k=0;k<8;++k)//header pointer is now invalid
		dlist_appendtoarray(list+k, output);
	time2=__rdtsc();

	if(loud)
	{
		size_t total_csize=output[0]->count-offset0;
		printf("AC encode:  %lld cycles, %lf c/byte\n", time2-time1, (double)(time2-time1)/count);
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", count, (int)total_csize, (double)count/total_csize, (double)(total_csize<<3)/count);
#ifdef AC_MEASURE_PREDICTION
		//printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", (count+7)>>3);
		for(int k=0;k<8;++k)
			printf("%2d\t%5d\t%lf\n", k, (int)list[k].nobj, (double)count/(list[k].nobj<<3));
		
		printf("Preview:\n");
		int kprint=total_csize<200?(int)total_csize:200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", output[0]->data[offset0+k]&0xFF);
		printf("\n");
	}

	for(int k=0;k<8;++k)
		dlist_clear(list+k);

	return 1;
}
const void* abac0a_decode(const void *src_start, const void *src_end, unsigned char *dst, int count, int bytestride, int loud)
{
	const unsigned char *sizes, *data;
	unsigned char *buffer;
	const unsigned char *ptr[8], *end[8];
	long long time1, time2;
	__m128i
		ones=_mm_set1_epi32(-1), ones32=_mm_set1_epi32(1), three=_mm_set1_epi32(3),
		prob_min=_mm_set1_epi16(1), prob_half=_mm_set1_epi16(0x4000), prob_max=_mm_set1_epi16(0x7FFE),
		prob=prob_half,
#ifdef ENABLE_CONFIDENCE
		prob_correct=prob_half,
#endif
		limit=_mm_set1_epi32(0x1000000);
	__m128i shuf_pack=_mm_set_epi8(15, 14, 11, 10, 7, 6, 3, 2, 13, 12, 9, 8, 5, 4, 1, 0);
	__m128i code[2], start[2], range[2], r2[2];
	ALIGN(16) unsigned a_code[8], a_start[8], a_range[8];

	if(!src_start||!src_end||!dst||!count||!bytestride)
		LOG_ERROR("abac4_decode(src_start=%p, src_end=%p, dst=%p, count=%d, bytestride=%d)", src_start, src_end, dst, count, bytestride);

	time1=__rdtsc();
	data=(const unsigned char*)src_start;
	buffer=(unsigned char*)dst;
	
	int headercount=1+8;
	if((int)(headercount*sizeof(int))>=(unsigned char*)src_end-data)
		LOG_ERROR("File is %d bytes < %d", (int)((unsigned char*)src_end-data), headercount*sizeof(int));
	int tag;
	memcpy(&tag, data, sizeof(int));
	if(tag!=tag_ac0A)
		LOG_ERROR("Invalid tag 0x%08X, expected 0x%08X", tag, tag_ac0A);
	sizes=data+sizeof(int);
	data=sizes+8*sizeof(int);

	size_t offset=0;
	for(int k=0;k<8;++k)
	{
		int size;
		memcpy(&size, sizes+k*sizeof(int), sizeof(int));
		ptr[k]=data+offset;
		offset+=size;
		end[k]=data+offset;
	}
	if((unsigned char*)src_end<end[7])
		LOG_ERROR("File is %d bytes, header requires %d bytes", (int)((size_t)src_end-(size_t)src_start), (int)((size_t)end[7]-(size_t)src_start));
	
	start[0]=_mm_setzero_si128();
	start[1]=_mm_setzero_si128();
	range[0]=_mm_set1_epi32(0x7FFFFFFF);//the MSB is the 2's comp sign, all SSE2 cmp instructions are signed
	range[1]=_mm_set1_epi32(0x7FFFFFFF);
	for(int k=0;k<8;++k)
	{
		ptr[k]+=sizeof(int);
		if(end[k]<ptr[k])
			LOG_ERROR("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
		a_code[k]=ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1];//big endian
	}
	code[0]=_mm_loadu_si128((__m128i*)a_code);
	code[1]=_mm_loadu_si128((__m128i*)a_code+1);
	for(int ks=0;ks<count;++ks, dst+=bytestride)
	{
		print_reg(code, 32, "%d code", ks);
		print_reg(start, 32, "%d start", ks);
		print_reg(range, 32, "%d range", ks);
		
		//calculate probability of zero bit
#ifdef ENABLE_CONFIDENCE
		__m128i p0=_mm_sub_epi16(prob, prob_half);

		p0=_mm_slli_epi16(p0, 1);
		prob_correct=_mm_slli_epi16(prob_correct, 1);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_mulhi_epu16(p0, p0);
		p0=_mm_srli_epi16(p0, 1);
		prob_correct=_mm_srli_epi16(prob_correct, 1);

		p0=_mm_add_epi16(p0, prob_half);
#else
		__m128i p0=prob;
#endif

		//clamp probability			eg: {0, 1, half, FFFF}
		__m128i cmp=_mm_cmplt_epi16(prob_min, p0);//{0, FFFF, FFFF, FFFF}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		p0=_mm_or_si128(p0, cmp);

		cmp=_mm_cmplt_epi16(p0, prob_max);//{FFFF, FFFF, FFFF, 0}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_max);
		p0=_mm_or_si128(p0, cmp);

		//calculate middle (r2)
		__m128i t0=_mm_shuffle_epi8(range[0], shuf_pack);
		__m128i t1=_mm_shuffle_epi8(range[1], shuf_pack);
		__m128i r_lo=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(1, 0, 1, 0)));
		__m128i r_hi=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(3, 2, 3, 2)));
		p0=_mm_slli_epi16(p0, 1);
		__m128i a=_mm_mulhi_epu16(r_lo, p0);
		__m128i b0=_mm_mullo_epi16(r_hi, p0);
		__m128i b1=_mm_mulhi_epu16(r_hi, p0);
		p0=_mm_srli_epi16(p0, 1);

		r2[0]=_mm_unpacklo_epi16(b0, b1);
		r2[1]=_mm_unpackhi_epi16(b0, b1);

		__m128i d0=_mm_unpacklo_epi16(a, _mm_setzero_si128());
		__m128i d1=_mm_unpackhi_epi16(a, _mm_setzero_si128());

		r2[0]=_mm_add_epi32(r2[0], d0);
		r2[1]=_mm_add_epi32(r2[1], d1);

		//avoid degenerate range		r2[k]+=(r2[k]==0)-(r2[k]==range[k]);
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		__m128i cmp2=_mm_cmpeq_epi32(r2[0], range[0]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[0]=_mm_add_epi32(r2[0], cmp);
		
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		cmp2=_mm_cmpeq_epi32(r2[0], range[1]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[1]=_mm_add_epi32(r2[1], cmp);

		__m128i bitlo=_mm_add_epi32(start[0], r2[0]);
		__m128i bithi=_mm_add_epi32(start[1], r2[1]);
		bitlo=_mm_cmplt_epi32(bitlo, code[0]);
		bithi=_mm_cmplt_epi32(bithi, code[1]);

		//pack low 16-bit halves of each 32-bit int
		bitlo=_mm_shuffle_epi8(bitlo, shuf_pack);
		bithi=_mm_shuffle_epi8(bithi, shuf_pack);
		__m128i bit=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(bitlo), _mm_castsi128_ps(bithi), _MM_SHUFFLE(1, 0, 1, 0)));

		//store symbol
		int mask=_mm_movemask_epi8(bit);
		unsigned char sym=mask>>8&128|mask>>7&64|mask>>6&32|mask>>5&16|mask>>4&8|mask>>3&4|mask>>1&2|mask&1;
#ifdef ENABLE_GRAYCODE
		sym^=sym>>4;
		sym^=sym>>2;
		sym^=sym>>1;
#endif
		*dst=sym;
		bit=_mm_and_si128(bit, prob_min);//leave only the LSB

		print_reg(r2, 32, "%d r2", ks);
		
		print_reg(&prob, 16, "%d prob ", ks);
		print_reg(&p0, 16, "%d p0 ", ks);
		
#ifdef ENABLE_CONFIDENCE
		//update prob_correct			prob_correct=correct<<15|prob_correct>>1;
		cmp=_mm_cmplt_epi16(p0, prob_half);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		cmp=_mm_xor_si128(cmp, bit);

		cmp=_mm_slli_epi16(cmp, 14);//set pre-MSB, because setting MSB would negate the 2's comp value
		prob_correct=_mm_srli_epi16(prob_correct, 1);
		prob_correct=_mm_or_si128(prob_correct, cmp);
#endif

		//update prob zero				prob=!bit<<15|prob>>1;
		__m128i bit2=_mm_xor_si128(bit, prob_min);
		prob=_mm_srli_epi16(prob, 1);
		prob=_mm_or_si128(prob, _mm_slli_epi16(bit2, 14));

		print_reg(&prob, 16, "%d prob2 ", ks);
		print_reg(&bit, 16, "%d bit ", ks);

#if 1
		//update range
		bitlo=_mm_unpacklo_epi16(bit, _mm_setzero_si128());
		bithi=_mm_unpackhi_epi16(bit, _mm_setzero_si128());

		//r2[k]+=bit.m128i_i16[k]-!bit.m128i_i16[k];
		r2[0]=_mm_add_epi32(r2[0], bitlo);
		r2[1]=_mm_add_epi32(r2[1], bithi);
		r2[0]=_mm_sub_epi32(r2[0], _mm_sub_epi32(ones32, bitlo));
		r2[1]=_mm_sub_epi32(r2[1], _mm_sub_epi32(ones32, bithi));

		//start[k]+=r2[k]&-bit.m128i_i16[k];
		bitlo=_mm_cmpeq_epi32(bitlo, ones32);
		bithi=_mm_cmpeq_epi32(bithi, ones32);
		t0=_mm_and_si128(bitlo, r2[0]);
		t1=_mm_and_si128(bithi, r2[1]);
		start[0]=_mm_add_epi32(start[0], t0);
		start[1]=_mm_add_epi32(start[1], t1);
		
		//range[k]=((range[k]-r2[k])&-bit.m128i_i16[k])|(r2[k]&-!bit.m128i_i16[k]);
		range[0]=_mm_sub_epi32(range[0], r2[0]);
		range[1]=_mm_sub_epi32(range[1], r2[1]);
		range[0]=_mm_and_si128(range[0], bitlo);
		range[1]=_mm_and_si128(range[1], bithi);
		bitlo=_mm_xor_si128(bitlo, ones);
		bithi=_mm_xor_si128(bithi, ones);
		r2[0]=_mm_and_si128(r2[0], bitlo);
		r2[1]=_mm_and_si128(r2[1], bithi);
		range[0]=_mm_or_si128(range[0], r2[0]);
		range[1]=_mm_or_si128(range[1], r2[1]);

		//renormalize
		unsigned condition;
		for(;;)
		{
			t0=_mm_add_epi32(start[0], range[0]);
			t1=_mm_add_epi32(start[1], range[1]);
			t0=_mm_xor_si128(t0, start[0]);
			t1=_mm_xor_si128(t1, start[1]);
			t0=_mm_cmpeq_epi32(t0, limit);
			t1=_mm_cmpeq_epi32(t1, limit);
			condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
			if(!condition)
				break;
			_mm_store_si128((__m128i*)a_code  , code[0]);
			_mm_store_si128((__m128i*)a_code+1, code[1]);
			_mm_store_si128((__m128i*)a_start  , start[0]);
			_mm_store_si128((__m128i*)a_start+1, start[1]);
			_mm_store_si128((__m128i*)a_range  , range[0]);
			_mm_store_si128((__m128i*)a_range+1, range[1]);
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					++ptr[k];
					if(end[k]<ptr[k])
						LOG_ERROR("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
#ifdef PRINT_SSE2
					printf("RENORM bit %d %08X<-[%02X]\n", k, a_code[k], (int)ptr[k][-1]);
#endif
					a_code[k]=a_code[k]<<8|ptr[k][-1];

					a_start[k]<<=8;
					a_range[k]=a_range[k]<<8|0xFF;
				}
			}
			code[0]=_mm_load_si128((__m128i*)a_code);
			code[1]=_mm_load_si128((__m128i*)a_code+1);
			start[0]=_mm_load_si128((__m128i*)a_start);
			start[1]=_mm_load_si128((__m128i*)a_start+1);
			range[0]=_mm_load_si128((__m128i*)a_range);
			range[1]=_mm_load_si128((__m128i*)a_range+1);
		}
		t0=_mm_cmplt_epi32(range[0], three);
		t1=_mm_cmplt_epi32(range[1], three);
		condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
		if(condition)
		{
			_mm_store_si128((__m128i*)a_code  , code[0]);
			_mm_store_si128((__m128i*)a_code+1, code[1]);
			_mm_store_si128((__m128i*)a_start  , start[0]);
			_mm_store_si128((__m128i*)a_start+1, start[1]);
			_mm_store_si128((__m128i*)a_range  , range[0]);
			_mm_store_si128((__m128i*)a_range+1, range[1]);
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					int reg_idx=k>>2, idx=k&3;

					ptr[k]+=sizeof(int);
					if(end[k]<ptr[k])
						LOG_ERROR("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
#ifdef PRINT_SSE2
					printf("RENORM bit %d %08X<-[%08X]\n", k, a_code[k], ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1]);
#endif
					a_code[k]=ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1];//big endian

					a_start[k]=0;
					a_range[k]=0x7FFFFFFF;//because 1=0.9999...
				}
			}
			code[0]=_mm_load_si128((__m128i*)a_code);
			code[1]=_mm_load_si128((__m128i*)a_code+1);
			start[0]=_mm_load_si128((__m128i*)a_start);
			start[1]=_mm_load_si128((__m128i*)a_start+1);
			range[0]=_mm_load_si128((__m128i*)a_range);
			range[1]=_mm_load_si128((__m128i*)a_range+1);
		}
#ifdef PRINT_SSE2
		printf("\n");
#endif
#endif
	}
	time2=__rdtsc();

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", time2-time1, (double)(time2-time1)/count);
	return ptr[7];
}
