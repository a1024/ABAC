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
#include<sys/stat.h>
#include"lodepng.h"//for testing
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

	#define ENABLE_ALL
	#define UNIFORM

size_t get_filesize(const char *filename)
{
	struct stat info={0};
	int error=stat(filename, &info);
	if(error)
		return 0;
	return info.st_size;
}
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
void calc_histogram(unsigned char *buf, ptrdiff_t len, ptrdiff_t stride, int *hist)
{
	memset(hist, 0, 256*sizeof(int));
	for(ptrdiff_t k=0, end=len-(stride-1);k<end;k+=stride)
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
void set_alpha(unsigned char *buf, ptrdiff_t len)
{
	for(ptrdiff_t k=0;k<len;k+=4)
		buf[k+3]=0xFF;
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
int hist[256*4];
unsigned short pred[]=
{
#if 0
	0xFF, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,//hot zero
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
#endif

#if 0
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,//uniform
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
#endif

#if 0
	0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,//exponential
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
	0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
#endif
	
#if 1
	0xFF, 0xFE, 0xFD, 0xFC, 0xFB, 0xFA, 0xF9, 0xF8, 0xF7, 0xF6, 0xF5, 0xF4, 0xF3, 0xF2, 0xF1, 0xF0,//(bad) ramp prediction
	0xEF, 0xEE, 0xED, 0xEC, 0xEB, 0xEA, 0xE9, 0xE8, 0xE7, 0xE6, 0xE5, 0xE4, 0xE3, 0xE2, 0xE1, 0xE0,
	0xDF, 0xDE, 0xDD, 0xDC, 0xDB, 0xDA, 0xD9, 0xD8, 0xD7, 0xD6, 0xD5, 0xD4, 0xD3, 0xD2, 0xD1, 0xD0,
	0xCF, 0xCE, 0xCD, 0xCC, 0xCB, 0xCA, 0xC9, 0xC8, 0xC7, 0xC6, 0xC5, 0xC4, 0xC3, 0xC2, 0xC1, 0xC0,
	0xBF, 0xBE, 0xBD, 0xBC, 0xBB, 0xBA, 0xB9, 0xB8, 0xB7, 0xB6, 0xB5, 0xB4, 0xB3, 0xB2, 0xB1, 0xB0,
	0xAF, 0xAE, 0xAD, 0xAC, 0xAB, 0xAA, 0xA9, 0xA8, 0xA7, 0xA6, 0xA5, 0xA4, 0xA3, 0xA2, 0xA1, 0xA0,
	0x9F, 0x9E, 0x9D, 0x9C, 0x9B, 0x9A, 0x99, 0x98, 0x97, 0x96, 0x95, 0x94, 0x93, 0x92, 0x91, 0x90,
	0x8F, 0x8E, 0x8D, 0x8C, 0x8B, 0x8A, 0x89, 0x88, 0x87, 0x86, 0x85, 0x84, 0x83, 0x82, 0x81, 0x80,
	0x7F, 0x7E, 0x7D, 0x7C, 0x7B, 0x7A, 0x79, 0x78, 0x77, 0x76, 0x75, 0x74, 0x73, 0x72, 0x71, 0x70,
	0x6F, 0x6E, 0x6D, 0x6C, 0x6B, 0x6A, 0x69, 0x68, 0x67, 0x66, 0x65, 0x64, 0x63, 0x62, 0x61, 0x60,
	0x5F, 0x5E, 0x5D, 0x5C, 0x5B, 0x5A, 0x59, 0x58, 0x57, 0x56, 0x55, 0x54, 0x53, 0x52, 0x51, 0x50,
	0x4F, 0x4E, 0x4D, 0x4C, 0x4B, 0x4A, 0x49, 0x48, 0x47, 0x46, 0x45, 0x44, 0x43, 0x42, 0x41, 0x40,
	0x3F, 0x3E, 0x3D, 0x3C, 0x3B, 0x3A, 0x39, 0x38, 0x37, 0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30,
	0x2F, 0x2E, 0x2D, 0x2C, 0x2B, 0x2A, 0x29, 0x28, 0x27, 0x26, 0x25, 0x24, 0x23, 0x22, 0x21, 0x20,
	0x1F, 0x1E, 0x1D, 0x1C, 0x1B, 0x1A, 0x19, 0x18, 0x17, 0x16, 0x15, 0x14, 0x13, 0x12, 0x11, 0x10,
	0x0F, 0x0E, 0x0D, 0x0C, 0x0B, 0x0A, 0x09, 0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
#endif
};
#if 0
int xhist[65536]={0}, yhist[65536]={0};
#endif
void calc_pairwisehist(unsigned char *buf, int iw, int ih, int bytestride, int *xhist, int *yhist)
{
	memset(xhist, 0, 65536*sizeof(int));
	memset(yhist, 0, 65536*sizeof(int));
	for(int ky=1;ky<ih;++ky)
	{
		for(int kx=1;kx<iw;++kx)
		{
			unsigned char
				s00=buf[bytestride*(iw*(ky-1)+kx-1)],
				s01=buf[bytestride*(iw*(ky-1)+kx)],
				s10=buf[bytestride*(iw* ky   +kx-1)];
			++xhist[s01<<8|s00];
			++xhist[s00<<8|s01];

			++yhist[s10<<8|s00];
			++yhist[s00<<8|s10];
		}
	}
}
double calc_entropy(int *hist, int nlevels)
{
	int sum=0;
	for(int k=0;k<nlevels;++k)
		sum+=hist[k];
	if(!sum)
		return 0;
	double entropy=0;
	for(int k=0;k<nlevels;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/sum;
			entropy-=p*log2(p);
		}
	}
	return entropy;
}
void print_pairwisehist(int *phist)
{
	int vmax=0;
	for(int ky=0;ky<256;++ky)
	{
		for(int kx=0;kx<256;++kx)
		{
			int pfreq=phist[ky<<8|kx];
			if(vmax<pfreq)
				vmax=pfreq;
		}
	}
	if(!vmax)
	{
		printf("Image was empty\n");
		return;
	}
	for(int ky=0;ky<256;++ky)
	{
		for(int kx=0;kx<256;++kx)
		{
			int qfreq=((phist[ky<<8|kx]<<4)+vmax-1)/vmax;
			if(qfreq>15)
				qfreq=15;
			printf("%X", qfreq);
		}
			//printf("%X", ((phist[ky<<8|kx]<<4)+(vmax>>1))/vmax);
		printf("\n");
	}
	printf("\n");
}
void differentiate_image(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				unsigned char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				if(kx||ky)
					sub-=128;
				buf[idx]-=sub;
			}
		}
	}
}
void integrate_image(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				unsigned char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				if(kx||ky)
					sub-=128;
				buf[idx]+=sub;
			}
		}
	}
}
/*int pyramid_getsize(int width){return (width+1)*((width<<1)-1);}
char pyramid_getchar(int width, int idx)
{
	int height=(width<<1)-1;
	++width;
	int x=idx%width, y=idx/width;
	if(x==width-1)
		return '\n';
	if(x>y)
		return ' ';
	if(x>height-1-y)
		return ' ';
	return '*';
}//*/
int main(int argc, char **argv)
{
	//int width=10,
	//	n=pyramid_getsize(width);
	//for(int k=0;k<n;++k)
	//	printf("%c", pyramid_getchar(width, k));

	//DCTtest();

	printf("EntropyBattle\n");
#if 1
	long long cycles;
	int iw=0, ih=0, nch0=0,
		nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	if(argc==2)
	{
		const char *fn=argv[1];
		printf("Opening %s\n", fn);
		cycles=__rdtsc();
		buf=stbi_load(fn, &iw, &ih, &nch0, nch);
		cycles=__rdtsc()-cycles;
		if(!buf)
			LOG_ERROR("Couldn't open %s", fn);
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		size_t formatsize=get_filesize(fn);
		if(formatsize)
			printf("Format Dec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", (double)cycles/(resolution*nch0), iw, ih, nch0, formatsize, (double)resolution*nch0/formatsize);
	}
	else
	{
		iw=1920, ih=1080, nch0=3,//1080*1920*3	640*480		50		4*4*1
			nch=4;
		resolution=(size_t)iw*ih, len=resolution*nch;
		buf=(unsigned char*)malloc(len);
		if(!buf)
			return 0;
		//srand((unsigned)__rdtsc());
	
#ifdef UNIFORM
		printf("Generating test data (uniform)...\n");
		fill_uniform(buf, len);
#else
		int unibits=256;
		printf("Generating test data (%d bit)...\n", unibits);
		fill_halfbinomial(buf, len, unibits);
#endif
	}
	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	size_t usize=len*nch0>>2;

	printf("\n");
	
	ArrayHandle cdata=0;
	const void *ptr, *end;
	
	int loud=0;
	//lowpassfilt(buf, len);//
	//lowpassfilt(buf, len);//
	//cvt2graycode(buf, len);//

	//print_bytes(buf, len);

	//2D differentiation
#if 0
	printf("Differentiating image...\n");

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<len;k+=nch)
			buf[k]=0xFF;
	}

	memcpy(b2, buf, len);//save copy

	differentiate_image(buf, iw, ih, nch0, nch);//add/sub 128 for pleasant image

	//lodepng_encode_file("out.PNG", buf, iw, ih, LCT_RGBA, 8);//save differentiated image

	//integrate_image(buf, iw, ih, nch0, nch);//check
	//compare_buffers(buf, b2, len, "Diff", 0);
#endif

	for(int kc=0;kc<nch0;++kc)
		calc_histogram(buf+kc, len, nch, hist+((size_t)kc<<8));
	//print_histogram(hist, 1);

	double entropy[6]={0};
	for(int kc=0;kc<nch0;++kc)
	{
		for(int k=0;k<256;++k)
		{
			int freq=hist[kc<<8|k];
			double p=(double)freq/(len>>2);
			if(freq)
				entropy[kc]+=-p*log2(p);

			//if(!kc)//
			//	printf("%3d %6d %lf\n", k, freq, p);//
		}
		printf("ch %d E = %lf / 8, optimal ratio = %lf\n", kc, entropy[kc], 8/entropy[kc]);
		entropy[4]+=entropy[kc];
	}
	entropy[4]/=nch0;
	printf("Av. E = %lf / 8, optimal ratio = %lf <- true limit\n", entropy[4], 8/entropy[4]);
#if 0
	for(int k=0;k<256;++k)
	{
		int freq=0;
		for(int kc=0;kc<nch0;++kc)
			freq+=hist[kc<<8|k];
		double p=(double)freq/usize;
		if(freq)
			entropy[5]+=-p*log2(p);
		//printf("%3d %6d %lf\n", k, (int)freq, p);//
	}
	printf("Comb. E = %lf / 8, optimal ratio = %lf\n\n", entropy[5], 8/entropy[5]);

	double xav=0, yav=0;
	for(int kc=0;kc<nch0;++kc)
	{
		calc_pairwisehist(buf+kc, iw, ih, 4, xhist, yhist);
		double
			xpe=calc_entropy(xhist, 65536),
			ype=calc_entropy(yhist, 65536);
		printf("ch %d  XPE %lf / 16  CR %lf  YPE %lf / 16  CR %lf\n", kc, xpe, 16/xpe, ype, 16/ype);
		xav+=xpe;
		yav+=ype;
		//printf("X-hist:\n");
		//print_pairwisehist(xhist);
		//printf("Y-hist:\n");
		//print_pairwisehist(yhist);
	}
	xav/=nch0;
	yav/=nch0;
	printf("Av pairwise entropy:  XE %lf / 16  OCR %lf  YE %lf / 16  OCR %lf\n\n", xav, 16/xav, yav, 16/yav);
#endif

#if 0
	memset(hist, 0, 4*sizeof(int));
	for(int k=0;k<resolution;++k)
		++hist[buf[k<<2]>>6];
	double se2=calc_entropy(hist, 4);
	printf("Entropy of 2 MSB: %lf / 2, CR %lf\n", se2, 2/se2);

	memset(hist, 0, 256*sizeof(size_t));
	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			unsigned char
				s00=buf[(iw* ky   +kx)<<2]>>6, s01=buf[(iw* ky   +kx+1)<<2]>>6,
				s10=buf[(iw*(ky+1)+kx)<<2]>>6, s11=buf[(iw*(ky+1)+kx+1)<<2]>>6;
			++hist[s11<<6|s10<<4|s01<<2|s00];
		}
	}
	double qe2=calc_entropy(hist, 256);
	printf("Quad entropy of 2 MSB: %lf / 8, CR %lf\n\n", qe2, 8/qe2);
	{
		size_t vmax=0;
		for(int k=0;k<256;++k)
		{
			if(vmax<hist[k])
				vmax=hist[k];
		}
		if(vmax)
		{
			for(int ky=0;ky<16;++ky)
			{
				for(int kx=0;kx<16;++kx)
				{
					int qfreq=(int)((((size_t)hist[ky<<4|kx]<<4)+vmax-1)/vmax);
					if(qfreq>15)
						qfreq=15;
					printf("%X ", qfreq);
				}
				printf("\n");
			}
		}
	}
#endif
	printf("\n");
	
#if 0
	printf("Haar wavelet\n");
	int nstages=0;//0
	short *b3=0;
	haar_2d_fwd(buf, iw, ih, nch0, nch, nstages, &b3);

	for(ptrdiff_t k=0;k<len;++k)//9 bit -> 16 bit
	{
		if(!((k+1)&3))
			b3[k]=0xFFFF;
		//else
		//	b3[k]<<=7;
	}
	printf("Saving...\n");
	lodepng_encode_file("out.PNG", b3, iw, ih, LCT_RGBA, 16);
	printf("Done\n");

	//haar_2d_inv(b3, iw, ih, nch0, nch, nstages, &b2);
	//compare_buffers(b2, buf, len, "Haar", 0);
	//memset(b2, 0, len);
#endif

	//test1 & test2
#if 0
	lz2_limit=1;
	//for(lz2_limit=1;lz2_limit<=3;++lz2_limit)
	{
		printf("lz2_limit %d\n", lz2_limit);

		printf("test1\n");
		cycles=test1_encode(buf, iw, ih, nch0, nch, &cdata);
		printf("Enc %lf CPB, ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		array_free(&cdata);
		printf("\n");

		printf("test2\n");
		cycles=test2_encode(buf, iw, ih, nch0, nch, &cdata);
		printf("Enc %lf CPB, ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		array_free(&cdata);
		printf("\n");
	}
#endif
#if 0
	{
		size_t overhead=0;
		printf("LZ\n");
		cycles=lz_encode(buf, iw, ih, nch, &cdata, &overhead);//red channel only
		printf("csize %lld, overhead %lld, diff %lld, ratio %lf\n", cdata->count, overhead, cdata->count-overhead, (double)resolution/(cdata->count-overhead));
		printf("Enc %lf CPB, ratio %lf\n", (double)cycles/resolution, (double)resolution/cdata->count);
		array_free(&cdata);
		printf("\n");
	}
#endif

	//LZ2D
#if 0
	{
		printf("LZ2D3\n");
		ArrayHandle mask=0, rle=0, coeff=0;
		cycles=__rdtsc();
		size_t savedbytes=lz2d3_encode(buf, iw, ih, nch0, nch, &mask, &rle, &coeff, 64);
		//size_t savedbytes=lz2d2_encode(buf, iw, ih, nch0, nch, &mask, &rle, &coeff);
		//size_t savedbytes=lz2d_encode(buf, iw, ih, nch0, nch, &mask, &coeff);
		cycles=__rdtsc()-cycles;
		size_t totalbytes=(size_t)nch0*iw*ih, overhead=rle->count*rle->esize+(coeff?coeff->count*coeff->esize:0), remaining=totalbytes-savedbytes+overhead;
		printf("saved %lld / %lld bytes, csize %lld, ratio %lf,  %lld rle %lld lz, overhead %lld bytes,  %lf CPB\n",
			savedbytes, totalbytes,
			remaining, (double)totalbytes/remaining,
			rle->count, coeff?coeff->count:0, overhead,
			(double)cycles/usize);

		lodepng_encode_file("out.PNG", mask->data, iw, ih, LCT_GREY, 8);//

		if(savedbytes)
		{
			array_free(&mask);
			array_free(&coeff);
		}
	}
#endif

	//test3
#if 1
	//for(int k=8;k<=256;k<<=1)//8=av(6, 10)
	//for(int k=8;k<=32;++k)		//10 works best
	{
		//printf("test3 %3d\t\t", k);
		printf("test3\n");
		cycles=__rdtsc();
		test3_encode(buf, iw, ih, nch0, nch, &cdata, 10);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		array_free(&cdata);
	}
	printf("\n");
#endif

	printf("\n");


	//Huffman
#ifdef ENABLE_ALL
	printf("Huffman (combined channels)\n");
	ArrayHandle udata=0;
	cycles=__rdtsc();
	huff_compress(buf, len, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	
	cycles=__rdtsc();
	huff_decompress(cdata->data, cdata->count, &udata);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

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
	if(nch0==4)
		ac0_encode(buf+3, len, 0, 8, nch, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);

	cycles=__rdtsc();
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=ac0_decode(ptr, end, b2  , len, 0, 8, nch, loud);
	ptr=ac0_decode(ptr, end, b2+1, len, 0, 8, nch, loud);
	ptr=ac0_decode(ptr, end, b2+2, len, 0, 8, nch, loud);
	if(nch0==4)
		ptr=ac0_decode(ptr, end, b2+3, len, 0, 8, nch, loud);
	else
		set_alpha(b2, len);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "AC", 0);
	memset(b2, 0, len);
	printf("\n");
#endif

	//ABAC
#ifdef ENABLE_ALL
	if(!(resolution&7))
	{
		printf("ABAC\n");
		cycles=__rdtsc();
		abac4_encode(buf  , (int)resolution, 0, 8, nch, &cdata, loud);
		abac4_encode(buf+1, (int)resolution, 0, 8, nch, &cdata, loud);
		abac4_encode(buf+2, (int)resolution, 0, 8, nch, &cdata, loud);
		if(nch0==4)
			abac4_encode(buf+3, (int)resolution, 0, 8, nch, &cdata, loud);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);

		cycles=__rdtsc();
		ptr=cdata->data;
		end=cdata->data+cdata->count;
		memset(b2, 0, len);
		ptr=abac4_decode(ptr, end, b2  , (int)resolution, 0, 8, nch, loud);
		ptr=abac4_decode(ptr, end, b2+1, (int)resolution, 0, 8, nch, loud);
		ptr=abac4_decode(ptr, end, b2+2, (int)resolution, 0, 8, nch, loud);
		if(nch0==4)
			ptr=abac4_decode(ptr, end, b2+3, (int)resolution, 0, 8, nch, loud);
		else
			set_alpha(b2, len);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		array_free(&cdata);
		compare_buffers(b2, buf, len, "ABAC", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//ABAC_SSE2
#ifdef ENABLE_ALL
	printf("ABAC_SSE2\n");
	cycles=__rdtsc();
	abac0a_encode(buf  , (int)resolution, nch, &cdata, loud);
	abac0a_encode(buf+1, (int)resolution, nch, &cdata, loud);
	abac0a_encode(buf+2, (int)resolution, nch, &cdata, loud);
	if(nch0==4)
		abac0a_encode(buf+3, (int)resolution, nch, &cdata, loud);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	
	cycles=__rdtsc();
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	memset(b2, 0, len);
	ptr=abac0a_decode(ptr, end, b2  , (int)resolution, nch, loud);
	ptr=abac0a_decode(ptr, end, b2+1, (int)resolution, nch, loud);
	ptr=abac0a_decode(ptr, end, b2+2, (int)resolution, nch, loud);
	if(nch0==4)
		ptr=abac0a_decode(ptr, end, b2+3, (int)resolution, nch, loud);
	else
		set_alpha(b2, len);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "ABAC_SSE2", 0);
	memset(b2, 0, len);
	printf("\n");
#endif

	//rANS
#ifdef ENABLE_ALL
	printf("rANS\n");
	cycles=__rdtsc();
	rans4_encode(buf, len, nch, 0, &cdata, loud, 0);
	//rans4_encode(buf, len, 1, 0, &cdata, loud, pred);
	cycles=__rdtsc()-cycles;
	//double cr=(double)usize/cdata->count, cr0=8/entropy[4];
	//printf("Enc CPB %lf ratio %lf%s\n", (double)cycles/usize, cr, cr>cr0?" IMPOSSIBLE":"");
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	
	cycles=__rdtsc();
	rans4_decode(cdata->data, cdata->count, len, nch, 0, b2, loud, 0);
	//rans4_decode(cdata->data, cdata->count, len, 1, 0, b2, loud, pred);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	array_free(&cdata);
	compare_buffers(b2, buf, len, "rANS", 0);
	memset(b2, 0, len);
	//printf("actual ratio %lf\n", (double)usize/(cdata->count-(12+(size_t)256*4*sizeof(short))));

#if 0
	printf("Adaptive rANS\n");

	Prob *table=rans5_preptable();

	cycles=__rdtsc();
	rans5_encode(buf, iw, ih, nch0, nch, &cdata, table);
	cycles=__rdtsc()-cycles;
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);

	free(table);
#endif

#if 0
	{
		cycles=__rdtsc();
		size_t csize=rans0c_testencode(buf, iw, ih, nch0, nch, &cdata);
		cycles=__rdtsc()-cycles;
		printf("rANS0C: CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/csize);
	}
#endif

#if 0
	entropy=0;
	for(int kp=0;kp<256;++kp)//H(Y|X) = -sum x in X: sum y in Y: p(x, y) lgN(p(x, y)/p(x))		== -sum p lgN(p)
	{
		double prediction=(double)pred[kp]/0x10000;
		for(int ka=0;ka<256;++ka)
		{
			double freq=(double)hist[ka]/len;
			if(freq)
			{
				double term=prediction*freq*log2(freq);
				if(isfinite(term))
					entropy-=term;
			}
		}
	}

	//for(int k=0;k<256;++k)
	//{
	//	if(hist[k])
	//	{
	//		double freq=(double)hist[k]/len, prediction=(double)pred[k]/0x10000;
	//		printf("%3d %lf * %lf\n", k, freq, -log2(prediction));
	//		entropy-=freq*log2(prediction);
	//	}
	//}
	printf("My entropy: %lf\nOptimal compression ratio with ramp predictor: %lf\n", entropy, 8/entropy);
#endif

	printf("\n");

#if 0
	if(!(iw&127)&&!(ih&127))
	{
		printf("Binary block rANS (dimensions divisible by 128)\n");
		unsigned blockdims[]={8, 16, 32, 64, 128};
		unsigned char bitcounts[]={2, 4, 8};
		for(int k=0;k<COUNTOF(blockdims);++k)
		{
			unsigned char blockdim=blockdims[k];
			for(int k2=0;k2<COUNTOF(bitcounts);++k2)
			{
				unsigned char bitcount=bitcounts[k2];
				cycles=__rdtsc();
				rans0b_testencode(buf, iw, ih, 24, 4, blockdim, bitcount, &cdata);
				cycles=__rdtsc()-cycles;
				printf("TestEnc(block %d, bits %d) CPB %lf ratio %lf\n", blockdim, bitcount, (double)cycles/usize, (double)usize/cdata->count);
				array_free(&cdata);
			}
		}
		printf("\n");
	}
#endif
#if 0
	printf("Binary rANS\n");
	int bestprob=0;
	double bestcr=0;
	for(int prob_MPS=0x8000;prob_MPS<0x10000;prob_MPS+=0x0100)
	{
		cycles=__rdtsc();
		rans8_testencode(buf, iw, ih, 24, 4, prob_MPS, &cdata);
		cycles=__rdtsc()-cycles;
		double cr=(double)usize/cdata->count;
		printf("TestEnc(MPS 0x%04X) CPB %lf ratio %lf\n", prob_MPS, (double)cycles/usize, cr);
		if(bestcr<cr)
			bestcr=cr, bestprob=prob_MPS;
		array_free(&cdata);
	}
	printf("\nBest prob_MPS 0x%04X ratio %lf\n\n", bestprob, bestcr);
	//printf("\n");
#endif
#endif

	//rANS_SSE2
#ifdef ENABLE_ALL
	if(!(nch&(nch-1))&&!(len&15))
	{
		printf("rANS_SSE2\n");
		cycles=__rdtsc();
		rans_sse2_encode(buf, len, nch, 0, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	
		cycles=__rdtsc();
		rans_sse2_decode(cdata->data, cdata->count, b2, len, nch, 0);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		array_free(&cdata);
		compare_buffers(b2, buf, len, "rANS_SSE2", 0);
		memset(b2, 0, len);
		printf("\n");
	}
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