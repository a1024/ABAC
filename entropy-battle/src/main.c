#include"battle.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif

	#define ENABLE_ALL
//	#define UNIFORM

void print_bytes(unsigned char *buf, size_t len)
{
	for(ptrdiff_t k=0;k<(ptrdiff_t)len;++k)
		printf("%d,", buf[k]);
		//printf("%02X-", buf[k]);
	printf("\n");
}
void print_histogram(size_t *hist, int all)
{
	size_t hmax=0;
	for(int sym=0;sym<256;++sym)
	{
		if(hmax<hist[sym])
			hmax=hist[sym];
	}
	for(int sym=0;sym<256;++sym)
	{
		if(all||hist[sym])
		{
			int printed=printf("%3d %7d ", sym, (int)hist[sym]);
			printed=79-printed;
			if(printed>0)
			{
				//printed=(int)((hist[sym]*printed+hmax-1)/hmax);
				printed=(int)(hist[sym]*printed/hmax);
				for(int k2=0;k2<printed;++k2)
					printf("*");
			}
			printf("\n");
		}
	}
	printf("\n");
}
void calc_histogram(unsigned char *buf, ptrdiff_t len, size_t *hist)
{
	memset(hist, 0, 256*sizeof(size_t));
	for(ptrdiff_t k=0;k<len;++k)
		++hist[buf[k]];
}
unsigned hammingweight(unsigned n)
{
#ifdef __GNUC__
	return __builtin_popcount(n);
#elif defined _MSC_VER
	return __popcnt(n);
#else
	//https://helloacm.com/c-coding-exercise-number-of-1-bits-revisited/
	x-=(x>>1)&0x5555555555555555;				//put count of each 2 bits into those 2 bits
	x=(x&0x3333333333333333)+((x>>2)&0x3333333333333333);	//put count of each 4 bits into those 4 bits
	x=(x+(x>>4))&0x0F0F0F0F0F0F0F0F;			//put count of each 8 bits into those 8 bits
	return x*0x0101010101010101>>56;//returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
#endif
}
unsigned sample(unsigned unibits)
{
	unsigned ret=0;
	for(;unibits>7;unibits-=8)
		ret+=hammingweight(rand()&0xFF);
	if(unibits)
		ret+=hammingweight(rand()&((1<<unibits)-1));
	return ret;
}
void fill_uniform(unsigned char *buf, ptrdiff_t len)
{
	for(ptrdiff_t k=0;k<len;++k)
		buf[k]=rand();
}
void fill_halfbinomial(unsigned char *buf, ptrdiff_t len, int unibits)
{
	int mask=(1<<unibits)-1;
	for(ptrdiff_t k=0;k<len;++k)
	{
		int val=sample(unibits)-sample(unibits);
		val^=-(val<0);//no adding one, such that {..., -2->1, -1->0}
		buf[k]=(unsigned char)val;
	}
}
void lowpassfilt(unsigned char *buf, ptrdiff_t len)
{
	for(ptrdiff_t k=0;k<len-1;++k)
		buf[k]=(buf[k]+buf[k+1])>>1;
}
void cvt2graycode(unsigned char *buf, ptrdiff_t len)
{
	for(ptrdiff_t k=0;k<len;++k)
		buf[k]^=buf[k]>>1;
}
void compare_buffers(unsigned char *b1, unsigned char *b2, ptrdiff_t len, const char *name, int backward)//backward is useless: ANS encodes in reverse, decodes forward
{
	int inc=(backward<<1)-1;
	for(ptrdiff_t k=backward?len-1:0;k>=0&&k<len;k+=inc)
	{
		if(b1[k]!=b2[k])
		{
			if(backward)
				printf("%s error at %d - %d: 0x%02X != 0x%02X\n", name, (int)len-1, (int)(len-1-k), b1[k], b2[k]);
			else
				printf("%s error at %d: 0x%02X != 0x%02X\n", name, (int)k, b1[k], b2[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}
size_t hist[256];
int main(int argc, char **argv)
{
	printf("EntropyBattle\n");
#if 1
	int	resolution=1080*1920,//1080*1920	640*480		50
		nch=4;
	size_t len=resolution*nch;
	unsigned char *buf=(unsigned char*)malloc(len),
		*b2=(unsigned char*)malloc(len);
	srand((unsigned)__rdtsc());
	
#ifdef UNIFORM
	printf("Generating test data (uniform)...\n");
	fill_uniform(buf, len);
#else
	int unibits=256;
	printf("Generating test data (%d bit)...\n", unibits);
	fill_halfbinomial(buf, len, unibits);
#endif

	calc_histogram(buf, len, hist);
	//print_histogram(hist, 1);
	printf("\n");
	
	ArrayHandle cdata=0;
	const void *ptr, *end;
	
	long long cycles;
	int loud=0;
	//lowpassfilt(buf, len);//
	//lowpassfilt(buf, len);//
	//cvt2graycode(buf, len);//

	//print_bytes(buf, len);


	double entropy=0;
	for(int k=0;k<256;++k)
	{
		if(hist[k])
		{
			double p=(double)hist[k]/len;
			entropy+=-p*log2(p);
		}
	}
	printf("Entropy = %lf / 8, best possible ratio = %lf\n\n", entropy, 8/entropy);


	//Huffman
#ifdef ENABLE_ALL
	printf("Huffman\n");
	ArrayHandle udata=0;
	cycles=__rdtsc();
	huff_compress(buf, len, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);
	
	cycles=__rdtsc();
	huff_decompress(cdata->data, cdata->count, &udata);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(buf, udata->data, len, "Huffman", 0);
	array_free(&udata);
	printf("\n");
#endif

	//AC
#ifdef ENABLE_ALL
	printf("AC\n");
	cycles=__rdtsc();
	ac0_encode(buf  , len, 0, 8, nch, &cdata, loud);
	ac0_encode(buf+1, len, 0, 8, nch, &cdata, loud);
	ac0_encode(buf+2, len, 0, 8, nch, &cdata, loud);
	ac0_encode(buf+3, len, 0, 8, nch, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);

	cycles=__rdtsc();
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=ac0_decode(ptr, end, b2  , len, 0, 8, nch, loud);
	ptr=ac0_decode(ptr, end, b2+1, len, 0, 8, nch, loud);
	ptr=ac0_decode(ptr, end, b2+2, len, 0, 8, nch, loud);
	ptr=ac0_decode(ptr, end, b2+3, len, 0, 8, nch, loud);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "AC", 0);
	printf("\n");
#endif

	//ABAC
#ifdef ENABLE_ALL
	printf("ABAC\n");
	cycles=__rdtsc();
	abac4_encode(buf  , resolution, 0, 8, nch, &cdata, loud);
	abac4_encode(buf+1, resolution, 0, 8, nch, &cdata, loud);
	abac4_encode(buf+2, resolution, 0, 8, nch, &cdata, loud);
	abac4_encode(buf+3, resolution, 0, 8, nch, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);

	cycles=__rdtsc();
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=abac4_decode(ptr, end, b2  , resolution, 0, 8, nch, loud);
	ptr=abac4_decode(ptr, end, b2+1, resolution, 0, 8, nch, loud);
	ptr=abac4_decode(ptr, end, b2+2, resolution, 0, 8, nch, loud);
	ptr=abac4_decode(ptr, end, b2+3, resolution, 0, 8, nch, loud);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "ABAC", 0);
	printf("\n");
#endif

	//ABAC_SSE2
#ifdef ENABLE_ALL
	printf("ABAC_SSE2\n");
	cycles=__rdtsc();
	abac0a_encode(buf  , resolution, nch, &cdata, loud);
	abac0a_encode(buf+1, resolution, nch, &cdata, loud);
	abac0a_encode(buf+2, resolution, nch, &cdata, loud);
	abac0a_encode(buf+3, resolution, nch, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);
	
	cycles=__rdtsc();
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=abac0a_decode(ptr, end, b2  , resolution, nch, loud);
	ptr=abac0a_decode(ptr, end, b2+1, resolution, nch, loud);
	ptr=abac0a_decode(ptr, end, b2+2, resolution, nch, loud);
	ptr=abac0a_decode(ptr, end, b2+3, resolution, nch, loud);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "ABAC_SSE2", 0);
	printf("\n");
#endif

	//rANS
#ifdef ENABLE_ALL
	printf("rANS\n");
	cycles=__rdtsc();
	rans4_encode(buf, len, nch, 0, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);
	
	cycles=__rdtsc();
	rans4_decode(cdata->data, cdata->count, len, nch, 0, b2, loud);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "rANS", 0);
	//printf("actual ratio %lf\n", (double)len/(cdata->count-(12+(size_t)256*4*sizeof(short))));
	printf("\n");
#endif

	//rANS_SSE2
#ifdef ENABLE_ALL
	printf("rANS_SSE2\n");
	cycles=__rdtsc();
	rans_sse2_encode(buf, len, nch, 0, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);
	
	cycles=__rdtsc();
	rans_sse2_decode(cdata->data, cdata->count, b2, len, nch, 0);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "rANS_SSE2", 0);
	printf("\n");
#endif




#if 0
	printf("ArANS\n");
	cycles=__rdtsc();
	arans_encode(buf, len, 1, &cdata);
	cycles=__rdtsc()-cycles;
	printf("cycles %lld\n%lld/%lld = %lf\n", cycles, len, cdata->count, (double)len/cdata->count);
#endif

	//AC debug
#if 0
	long long pat=0x00FF00FFFF00FF00;//
	memfill(buf, &pat, len, 8);//
	//memset(buf, 0x55, len);//
	loud=1;
	int bitdepth=8;
	size_t enclen=100;
	ac0_encode(buf  , enclen, 0, bitdepth, nch, &cdata, loud);
	ac0_encode(buf+1, enclen, 0, bitdepth, nch, &cdata, loud);
	ac0_encode(buf+2, enclen, 0, bitdepth, nch, &cdata, loud);
	ac0_encode(buf+3, enclen, 0, bitdepth, nch, &cdata, loud);
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=ac0_decode(ptr, end, b2  , enclen, 0, bitdepth, nch, loud, buf  );
	ptr=ac0_decode(ptr, end, b2+1, enclen, 0, bitdepth, nch, loud, buf+1);
	ptr=ac0_decode(ptr, end, b2+2, enclen, 0, bitdepth, nch, loud, buf+2);
	ptr=ac0_decode(ptr, end, b2+3, enclen, 0, bitdepth, nch, loud, buf+3);
#endif

	//uABS debug
#if 0
	size_t enclen=4096;
	uabs_encode_ch(buf, enclen, 0, 8, 1, &cdata);
	printf("uABS ratio %lf\n", (double)enclen/cdata->count);
	uabs_decode_ch(cdata->data, cdata->count, b2, enclen, 0, 8, 1
#ifdef ENABLE_GUIDE
		, buf
#endif
	);
	compare_buffers(b2, buf, resolution, "ABAC_SSE2", 0);
#endif

	//ABAC_SSE2 debug
#if 0
	resolution=64;
	printf("resolution %d\n", resolution);
	abac0a_encode(buf, resolution, 1, &cdata, 0);
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	ptr=abac0a_decode(ptr, end, b2, resolution, 1, 0);
	compare_buffers(b2, buf, resolution, "ABAC_SSE2", 0);
#endif

	//rANS_SSE2 debug
#if 0
	size_t offset=0, len2=len;
	//size_t offset=len/2, len2=len/2;//coincedence permutation
	//size_t offset=len/4, len2=len*3/4;
	//size_t offset=len-16, len2=16;
	rans_sse2_encode(buf+offset, len2, nch, 0, &cdata);
	rans_sse2_decode(cdata->data, cdata->count, b2+offset, len2, nch, 0, buf+offset);

	compare_buffers(b2+offset, buf+offset, len2, "rANS_SSE2", 0);
#endif


	free(buf);
	free(b2);
#endif
	printf("Done.\n");
	pause();
	return 0;
}