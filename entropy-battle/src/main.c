#include"battle.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
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
	int	resolution=1080*1920,//1080*1920	640*480
		nch=4;
	size_t len=resolution*nch;
	unsigned char *buf=(unsigned char*)malloc(len),
		*b2=(unsigned char*)malloc(len);
	//srand((unsigned)__rdtsc());
	printf("Generating test data...\n");
	fill_halfbinomial(buf, len, 16);
	//calc_histogram(buf, len, hist);
	//print_histogram(hist, 1);
	printf("\n");
	
	ArrayHandle cdata=0;
	const void *ptr, *end;

	long long cycles;
	int loud=0;


	//Huffman
#if 1
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

	//ABAC
#if 1
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
#if 1
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
#if 1
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
	printf("\n");
#endif

	//rANS_SSE2
#if 1
	printf("rANS_SSE2\n");
	cycles=__rdtsc();
	rans_sse2_encode(buf, len, nch, 0, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/len, (double)len/cdata->count);
	
	cycles=__rdtsc();
	rans_sse2_decode(cdata->data, cdata->count, b2, len, nch, 0, buf);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/len);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "rANS_SSE2", 0);
	printf("\n");
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