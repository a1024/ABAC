#include"battle.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include"lodepng.h"//for testing
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;

//	#define ENABLE_ALL
//	#define UNIFORM

#if 0
size_t get_filesize(const char *filename)
{
	struct stat info={0};
	int error=stat(filename, &info);
	if(error)
		return 0;
	return info.st_size;
}
#endif
void print_shorts2d(short *buf, int bw, int bh)
{
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%d\t", buf[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
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
void image_pred(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				if(!kx)
					kx=0;
				unsigned char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0, sub;

				//sub=left;//diff x
				//sub=(left*3+top*3+topleft*2)>>3;//weighted average
				sub=left+top-topleft;//diff 2d
				//sub=top;//diff y

				//unsigned char p=left+top-topleft;//Paeth
				//sub=left;
				//if(abs(p-top)<abs(p-sub))
				//	sub=top;
				//if(abs(p-topleft)<abs(p-sub))
				//	sub=topleft;

				if(kx||ky)
					sub-=128;
				buf[idx]-=sub;

				//if(!kc)//
				//	buf[idx]=0;//
			}
		}
	}
}
void squeeze_8bit_lossy(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	size_t res=(size_t)iw*ih*bytestride, maxdim=MAXVAR(iw, ih);
	short *b2=(short*)malloc(res*sizeof(short));
	short *temp=(short*)malloc(maxdim*sizeof(short));
	if(!b2||!temp)
		return;
	for(int k=0;k<res;++k)
		b2[k]=buf[k];
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	for(int kc=0;kc<nch;++kc)
		squeeze_2d_fwd(b2+kc, psizes, 0, nsizes, bytestride, 0, temp);//squeeze transform from JPEG XL

	array_free(&sizes);
	for(int k=0;k<res;++k)
		buf[k]=b2[k]+(128&-((k+1)&3));
	free(temp);
	free(b2);
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
//void DCTtest();//
//void DCTtest2();//
//void test4();
void save_32bit(const char *filename, const int *buf, int iw, int ih, int nch, int saveas8bit)
{
	LodePNGColorType type;
	switch(nch)
	{
	case 1:type=LCT_GREY;break;
	case 3:type=LCT_RGB;break;
	case 4:type=LCT_RGBA;break;
	default:LOG_ERROR("Bitmap type error");return;
	}
	int n=iw*ih*nch;
	if(saveas8bit)
	{
		unsigned char *b3=(unsigned char*)malloc(n);
		if(!b3)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned val=buf[k];
			b3[k]=(unsigned char)(val>>24);
		}
		lodepng_encode_file(filename, b3, iw, ih, type, 8);
		free(b3);
	}
	else
	{
		short *b2=(short*)malloc(n*sizeof(short));
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]>>16;
			b2[k]=_byteswap_ushort(val);
		}
		lodepng_encode_file(filename, (unsigned char*)b2, iw, ih, type, 16);
		free(b2);
	}
}
void save_16bit(const char *filename, const short *buf, const short *sub_b2, int iw, int ih, int nch, int val_offset, int val_shift, int saveas8bit)
{
	LodePNGColorType type;
	switch(nch)
	{
	case 1:type=LCT_GREY;break;
	case 3:type=LCT_RGB;break;
	case 4:type=LCT_RGBA;break;
	default:LOG_ERROR("Bitmap type error");return;
	}
	int n=iw*ih*nch;
	if(saveas8bit)
	{
		unsigned char *b3=(unsigned char*)malloc(n);
		if(!b3)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]+val_offset;
			if(sub_b2)
				val-=sub_b2[k];
			val<<=val_shift;
			b3[k]=(unsigned char)(val>>8);
		}
		lodepng_encode_file(filename, b3, iw, ih, type, 8);
		free(b3);
	}
	else
	{
		short *b2=(short*)malloc(n*sizeof(short));
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		for(int k=0;k<n;++k)
		{
			unsigned short val=buf[k]+val_offset;
			if(sub_b2)
				val-=sub_b2[k];
			val<<=val_shift;
			b2[k]=_byteswap_ushort(val);
		}
		lodepng_encode_file(filename, (unsigned char*)b2, iw, ih, type, 16);
		free(b2);
	}
}
void save_mono8(const char *filename, unsigned char *buf, int iw, int ih, int stride)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	unsigned char *b2=(unsigned char*)malloc(res);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(ptrdiff_t k=0;k<res;++k)
		b2[k]=buf[k*stride];
	lodepng_encode_file(filename, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}
void save_squeeze(const char *name, short *buf, DWTSize *sizes, int nsizes, int stride, int maxnbits)
{
	if(nsizes<2)
		return;
	char fn[256]={0};
	int iw=sizes->w, ih=sizes->h;
	int lsbw=(iw>>1)+1, lsbh=(ih>>1)+1, res=lsbw*lsbh;//largest subband size
	unsigned char *b2=(unsigned char*)malloc((size_t)res*4);
	memset(b2, 0xFF, (size_t)res*4);
	//for(int k=0;k<res;++k)
	//	b2[k<<2|3]=0xFF;
	int offset=1<<(maxnbits-1);
	for(int k=1;k<nsizes;++k)
	{
		int px=sizes[k].w, py=sizes[k].h, dx=sizes[k-1].w, dy=sizes[k-1].h, xodd=dx&1, yodd=dy&1;
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
				b2[((px+xodd)*ky+kx)<<2  ]=(buf[stride*(iw* ky    +px+kx)]+offset)>>(maxnbits-8);
		}
		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
				b2[((px+xodd)*ky+kx)<<2|1]=(buf[stride*(iw*(py+ky)+px+kx)]+offset)>>(maxnbits-8);
		}
		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0;kx<px;++kx)
				b2[((px+xodd)*ky+kx)<<2|2]=(buf[stride*(iw*(py+ky)   +kx)]+offset)>>(maxnbits-8);
		}
		snprintf(fn, 128, "%s%02d.PNG", name, k);
		lodepng_encode_file(fn, b2, px+xodd, py+yodd, LCT_RGBA, 8);
	}
	free(b2);
}
void save_channel(unsigned char *buf, int iw, int ih, int stride, const char *format, ...)
{
	va_list args;
	va_start(args, format);
	vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);

	unsigned char *b2=(unsigned char*)malloc((size_t)iw*ih);
	if(!b2)
		return;
	for(int ks=0, kd=0, res=iw*ih;kd<res;ks+=stride, ++kd)
		b2[kd]=buf[ks];

	lodepng_encode_file(g_buf, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}
void save_CDF53(const char *name, unsigned char *buf, DWTSize *sizes, int nsizes, int stride)
{
	if(nsizes<2)
		return;
	char fn[256]={0};
	int iw=sizes->w, ih=sizes->h;
	int lsbw=(iw>>1)+1, lsbh=(ih>>1)+1, res=lsbw*lsbh;//largest subband size
	unsigned char *b2=(unsigned char*)malloc((size_t)res*4);
	memset(b2, 0xFF, (size_t)res*4);
	for(int k=1;k<nsizes;++k)
	{
		int px=sizes[k].w, py=sizes[k].h, dx=sizes[k-1].w, dy=sizes[k-1].h, xodd=dx&1, yodd=dy&1;
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
				b2[((px+xodd)*ky+kx)<<2  ]=buf[stride*(iw* ky    +px+kx)];
		}
		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
				b2[((px+xodd)*ky+kx)<<2|1]=buf[stride*(iw*(py+ky)+px+kx)];
		}
		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0;kx<px;++kx)
				b2[((px+xodd)*ky+kx)<<2|2]=buf[stride*(iw*(py+ky)   +kx)];
		}
		snprintf(fn, 128, "%s%02d.PNG", name, k);
		lodepng_encode_file(fn, b2, px+xodd, py+yodd, LCT_RGBA, 8);
	}
	free(b2);
}

void print_hist32(unsigned *hist, int nlevels, int all)
{
	size_t hmax=0;
	for(int sym=0;sym<nlevels;++sym)
	{
		if(hmax<hist[sym])
			hmax=hist[sym];
	}
	for(int sym=0;sym<nlevels;++sym)
	{
		if(all||hist[sym])
		{
			int printed=printf("%3d %7d ", sym, (int)hist[sym]);
			printed=79-printed;
			if(printed>0)
			{
				//printed=(int)((hist[sym]*printed+hmax-1)/hmax);
				printed=(int)((size_t)hist[sym]*printed/hmax);
				for(int k2=0;k2<printed;++k2)
					printf("*");
			}
			printf("\n");
		}
	}
	printf("\n");
}
void apply_transforms_fwd(unsigned char *buf, int bw, int bh)
{
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char

	//colortransform_ycocg_fwd((char*)buf, bw, bh);
	colortransform_xgz_fwd((char*)buf, bw, bh);
	//colortransform_xyz_fwd((char*)buf, bw, bh);

	for(int kc=0;kc<3;++kc)
		dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

	addbuf(buf, bw, bh, 3, 4, 128);

	free(temp);
	array_free(&sizes);
}
void apply_transforms_inv(unsigned char *buf, int bw, int bh)
{
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	
	addbuf(buf, bw, bh, 3, 4, 128);

	for(int kc=0;kc<3;++kc)
		//dwt2d_haar_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		dwt2d_cdf53_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//dwt2d_cdf97_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

	colortransform_ycocg_inv((char*)buf, bw, bh);
	//colortransform_xgz_inv((char*)buf, bw, bh);
	//colortransform_xyz_inv((char*)buf, bw, bh);

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char

	free(temp);
	array_free(&sizes);
}
long long accumulate_hist(unsigned *hist, int nbits)
{
	int nlevels=1<<nbits;
	long long sum=0;
	for(int sym=0;sym<nlevels;++sym)
		sum+=hist[sym];
	return sum;
}
double estimate_csize_from_hist(unsigned *hist, int nbits, int sum, int bypass)
{
	int nlevels=1<<nbits;
	double entropy=0;
	if(bypass)
		entropy=nbits;
	else
	{
		for(int sym=0;sym<nlevels;++sym)
		{
			unsigned freq=hist[sym];
			if(freq)
			{
				double p=(double)freq/sum;
				//printf("%lf\n", p);
				p*=0x10000-(nlevels-1);
				++p;
				p/=0x10000;
				entropy-=p*log2(p);
			}
		}
	}
	double usize=(double)sum*nbits/8;
	double CR=nbits/entropy;
	double csize=usize/CR;
	//double csize=((double)sum*nbits/8)/(nbits/entropy);
	return csize;
}
double get_mean(const unsigned char *buf, int res, int stride)
{
	long long sum=0;
	for(int k=0;k<res;++k)
		sum+=buf[k];
	double mean=(double)sum/res;
	return mean;
}
double get_var(const unsigned char *buf, int res, int stride, double mean)
{
	double sum=0;
	for(int k=0;k<res;++k)
	{
		double x=buf[k]-mean;
		sum+=x*x;
	}
	sum/=res;
	return sum;
}
double test9_estimate_cr(const unsigned char *buf, int bw, int bh, int loud)
{
	int res=bw*bh;
	double start=time_ms();
	double mean[3]=
	{
		get_mean(buf  , res, 4),
		get_mean(buf+1, res, 4),
		get_mean(buf+2, res, 4),
	};
	double conf[]=
	{
		get_var(buf  , res, 4, mean[0]),
		get_var(buf+1, res, 4, mean[1]),
		get_var(buf+2, res, 4, mean[2]),
	};
	double entropy=0;
	double gain=1/sqrt(8*M_PI*M_PI*M_PI*conf[0]*conf[1]*conf[2]);
	conf[0]=1/sqrt(conf[0]);
	conf[1]=1/sqrt(conf[1]);
	conf[2]=1/sqrt(conf[2]);
	for(int color=0;color<0x1000000;++color)
	{
		double
			r=((color    &255)-mean[0])*conf[0],
			g=((color>> 8&255)-mean[1])*conf[1],
			b=((color>>16&255)-mean[2])*conf[2];
		double prob=exp(-0.5*(r*r+b*b+g*g));
		prob*=gain;
		if(prob)
		{
			double size=-log2(prob);
			if(isfinite(size))
				entropy+=prob*size;
		}
	}
	double CR0=24/entropy, csize=ceil(res*3/CR0)+24, CR=res*3/csize;
	double end=time_ms();
	if(loud)
		printf("Test 9 estimate: usize %d csize %lf CR %lf elapsed %lfms\n\n", res*3, csize, CR, end-start);
	return csize;
}
double test8_estimate_cr(const unsigned char *buf, int bw, int bh, int loud)
{
	const int point=5;

	int headsize=1<<(8-point)*3, remsize=1<<point,
		masklo=remsize-1, maskhi=(1<<(8-point))-1,
		totalhsize=(headsize+(size_t)headsize*remsize*3)*sizeof(unsigned);
	unsigned *hist=(unsigned*)malloc(totalhsize);
	if(!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, totalhsize);
	unsigned *histhead=hist, *histr=histhead+headsize, *histg=histr+headsize*remsize, *histb=histg+headsize*remsize;
	int res=bw*bh;
	for(int k=0;k<res;++k)
	{
		unsigned char r=buf[k<<2]-112, g=buf[k<<2|1]-112, b=buf[k<<2|2]-112;
		//unsigned char r=buf[k<<2]-96, g=buf[k<<2|1]-96, b=buf[k<<2|2]-96;
		//unsigned char r=buf[k<<2]-64, g=buf[k<<2|1]-64, b=buf[k<<2|2]-64;
		//unsigned char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		int head=(b>>point&maskhi)<<((8-point)<<1)|(g>>point&maskhi)<<(8-point)|r>>point&maskhi;
		r&=masklo;
		g&=masklo;
		b&=masklo;
		++histhead[head];
		if(!head)
		{
			++histr[r];
			++histg[g];
			++histb[b];
		}
		else
		{
			++histr[remsize+r];
			++histg[remsize+g];
			++histb[remsize+b];
		}
		//++histr[head<<point|r];
		//++histg[head<<point|g];
		//++histb[head<<point|b];
	}
#if 0
	int inspector;
	print_hist32(histhead,                headsize, 1), inspector=rand()%headsize;
	print_hist32(histr+remsize*inspector, remsize , 1), inspector=rand()%headsize;
	print_hist32(histg+remsize*inspector, remsize , 1), inspector=rand()%headsize;
	print_hist32(histb+remsize*inspector, remsize , 1);
#endif

	int sum[]=
	{
		(int)accumulate_hist(histhead, (8-point)*3),
		(int)accumulate_hist(histr, point),
		(int)accumulate_hist(histg, point),
		(int)accumulate_hist(histb, point),
		(int)accumulate_hist(histr+remsize, point),
		(int)accumulate_hist(histg+remsize, point),
		(int)accumulate_hist(histb+remsize, point),
	};
	double
		s1=estimate_csize_from_hist(histhead, (8-point)*3, sum[0], 0), s2=0, s3=0, s4=0;
	for(int k=0;k<2;++k)
	//for(int k=0;k<headsize;++k)
	{
		s2+=estimate_csize_from_hist(histr+remsize*k, point, sum[1+3*k], k);
		s3+=estimate_csize_from_hist(histg+remsize*k, point, sum[2+3*k], k);
		s4+=estimate_csize_from_hist(histb+remsize*k, point, sum[3+3*k], k);
	}
	int u1=res*(8-point)*3/8, u2=res*point/8;
	//b2/=headsize;
	//b3/=headsize;
	//b4/=headsize;
	int overhead=(headsize+remsize*3)*sizeof(short);
	int total=overhead+(int)ceil(s1+s2+s3+s4);
	double CR=(double)res*3/total;
	if(loud)
	{
		printf("Test 8 estimate:\n");
		printf("\toverhead %d\n", overhead);
		printf("\thead     %14lf / %d = %lf\n", s1, u1, u1/s1);
		printf("\tred      %14lf / %d = %lf\n", s2, u2, u2/s2);
		printf("\tgreen    %14lf / %d = %lf\n", s3, u2, u2/s3);//[sic]
		printf("\tblue     %14lf / %d = %lf\n", s4, u2, u2/s4);
		printf("\n");
		printf("Total:     %d\n", total);
		printf("Unc. Size: %d\n", res*3);
		printf("CR:        %lf\n", CR);
		printf("\n");
	}
	free(hist);
	return total;
}

#if 0
int error_func_p16(int x);
long long error_func_p32(long long x);
void test5()
{
	double maxe16=0, maxe32=0;
	for(int k=0;k<256;++k)
	{
		//if(k==1)//
		//	k=1;//
		double pd=erf((double)k/256);

		int p16=error_func_p16(k<<8);
		long long p32=error_func_p32((long long)k<<24);

		double
			p16_as_pd=(double)p16/0x10000, e16=fabs(pd-p16_as_pd),
			p32_as_pd=(double)p32/0x10000, e32=fabs(pd-p32_as_pd);
		printf("%3d pd %.10lf\tp16 0x%08X = %lf\terr %g\tp32 0x%016llX = %lf\terr %g\n",
			k, pd, p16, p16_as_pd, e16, p32, p32_as_pd, e32);
		if(maxe16<e16)
			maxe16=e16;
		if(maxe32<e32)
			maxe32=e32;
	}
	printf("\n");
	printf("max p16 error 0x%08X = %lf\n", (int)(0x10000*maxe16), maxe16);
	printf("max p32 error 0x%08X = %lf\n", (int)(0x10000*maxe32), maxe32);

	printf("Done.\n");
	pause();
	exit(0);
}
#endif
#if 0
//#define COUNT 256
#define WIDTH 31
#define HEIGHT 17
short g_temp[MAXVAR(WIDTH, HEIGHT)];
short g_b1[]=
{
	195, 196, 196, 192, 195, 195, 192, 196, 195, 195, 193, 195, 192, 193, 193, 191, 192, 195, 193, 193, 193, 191, 192, 192, 192, 189, 192, 193, 192, 192, 193,
	195, 195, 194, 194, 195, 192, 193, 195, 195, 192, 193, 192, 195, 195, 193, 192, 194, 191, 193, 191, 191, 189, 192, 192, 191, 192, 191, 191, 191, 191, 191,
	196, 196, 195, 195, 195, 193, 192, 192, 195, 192, 193, 195, 192, 192, 193, 194, 192, 192, 191, 192, 191, 191, 192, 191, 189, 192, 191, 191, 191, 191, 191,
	195, 192, 194, 195, 196, 192, 192, 192, 192, 192, 192, 194, 191, 194, 192, 189, 192, 192, 189, 192, 191, 191, 192, 189, 189, 191, 192, 189, 191, 191, 191,
	195, 193, 192, 192, 195, 191, 191, 192, 192, 194, 194, 192, 194, 194, 192, 191, 192, 191, 191, 189, 191, 191, 192, 191, 192, 191, 189, 189, 191, 188, 191,
	191, 193, 193, 191, 194, 191, 192, 189, 192, 191, 194, 191, 191, 192, 191, 189, 191, 191, 189, 191, 192, 188, 187, 189, 191, 189, 191, 189, 189, 191, 189,
	193, 192, 193, 193, 195, 195, 192, 192, 192, 193, 191, 194, 192, 194, 191, 191, 192, 191, 191, 189, 189, 189, 188, 189, 191, 189, 188, 189, 188, 191, 189,
	196, 199, 196, 197, 200, 199, 200, 196, 199, 199, 199, 199, 198, 201, 200, 199, 197, 196, 199, 197, 197, 199, 196, 193, 195, 196, 195, 196, 196, 193, 195,
	147, 153, 150, 148, 154, 155, 158, 159, 157, 162, 166, 171, 172, 168, 169, 169, 171, 172, 177, 174, 174, 171, 175, 178, 178, 178, 177, 177, 178, 181, 183,
	143, 140, 140, 139, 149, 143, 134, 140, 153, 141, 129, 148, 158, 137, 130, 132, 132, 134, 136, 133, 133, 133, 136, 136, 133, 134, 132, 132, 130, 136, 136,
	156, 154, 156, 157, 157, 154, 154, 154, 158, 157, 150, 149, 149, 150, 149, 152, 150, 152, 148, 152, 151, 152, 153, 151, 145, 149, 149, 147, 149, 151, 147,
	152, 149, 148, 149, 150, 147, 147, 147, 146, 148, 150, 148, 147, 147, 150, 148, 150, 148, 146, 148, 146, 148, 148, 148, 149, 147, 144, 146, 143, 144, 146,
	143, 139, 142, 140, 142, 147, 138, 143, 143, 140, 142, 142, 142, 142, 142, 142, 142, 141, 141, 141, 142, 141, 140, 143, 139, 139, 140, 141, 137, 141, 142,
	155, 157, 154, 153, 151, 147, 145, 141, 143, 144, 141, 144, 148, 147, 143, 137, 141, 144, 140, 137, 139, 139, 137, 139, 141, 139, 142, 145, 144, 142, 141,
	147, 147, 145, 149, 144, 140, 140, 137, 141, 143, 141, 143, 143, 141, 138, 138, 140, 141, 140, 140, 139, 143, 134, 138, 142, 139, 143, 144, 141, 137, 135,
	135, 129, 136, 141, 134, 134, 130, 136, 138, 140, 137, 136, 137, 134, 137, 143, 140, 137, 136, 136, 137, 138, 135, 137, 139, 133, 140, 141, 140, 139, 139,
	143, 130, 137, 147, 142, 139, 134, 133, 138, 137, 133, 130, 138, 138, 135, 138, 143, 137, 133, 133, 130, 129, 133, 131, 134, 137, 138, 133, 133, 133, 135,

	//0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
	//0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
	//0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
	//0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,
	//0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
	//0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E, 0x2F,
	//0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
	//0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F,
};
void test6()
{
#if 1
	const char fn[]="E:/ML/dataset-kodak/kodim13.png";
	int iw=0, ih=0, nch0=0;
	unsigned char *buf=stbi_load(fn, &iw, &ih, &nch0, 4);
	if(!buf)
		LOG_ERROR("Couldn't open \'%s\'", fn);

	int maxdim=MAXVAR(iw, ih);
	short *temp=(short*)malloc(maxdim*sizeof(short));
	short *b2=(short*)malloc((size_t)iw*ih*sizeof(short));
	if(!temp||!b2)
		LOG_ERROR("Allocation error");

	for(int k=0, res=iw*ih;k<res;++k)
		b2[k]=buf[k<<2];//red channel
		//b2[k]=(buf[k<<2]+buf[k<<2|1]+buf[k<<2|2])/3;//average

	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
	squeeze_2d_fwd(b2, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, 0, temp);

	save_squeeze("squeeze-red", b2, (DWTSize*)sizes->data, (int)sizes->count, 1, 9);
	//save_16bit("squeeze.PNG", b2, 0, iw, ih, 1, 256, 16-9, 1);

	free(b2);
	free(buf);
	free(temp);
#endif

#if 0
	//printf("Fullscreen the console window\n");
	//pause();

	//for(int ky=0;ky<HEIGHT;++ky)
	//{
	//	for(int kx=0;kx<WIDTH;++kx)
	//	{
	//		double x=kx-(WIDTH-1)*0.5, y=ky-(HEIGHT-1)*0.5;
	//		g_b1[WIDTH*ky+kx]=(int)(0+255/(x*x+y*y+1));
	//	}
	//}
	//printf("Original:\n");
	//print_shorts2d(g_b1, WIDTH, HEIGHT);
	save_16bit("blob.PNG", g_b1, 0, WIDTH, HEIGHT, 1, 0, 8, 1);

	ArrayHandle sizes=dwt2d_gensizes(WIDTH, HEIGHT, 3, 3, 0);
	squeeze_2d_fwd(g_b1, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, 0, g_temp);
	//for(int k=0;k<(int)sizes->count-1;++k)
	//{
	//	squeeze_2d_fwd(g_b1, (DWTSize*)sizes->data, k, k+2, 1, 0, g_temp);
	//	
	//	//printf("DWT stage %d:\n", k);
	//	//print_shorts2d(g_b1, WIDTH, HEIGHT);
	//}
	save_16bit("blob_squeeze.PNG", g_b1, 0, WIDTH, HEIGHT, 1, 128, 8, 1);
#endif

	printf("Done.\n");
	pause();
	exit(0);
}
#endif
#if 0
double pdf(int x, double mean, double sdev, int nlevels)
{
	sdev*=sqrt(2.);
	return (erf((x+1-mean)/sdev)-erf((x-mean)/sdev)+1/(nlevels*4096.))/(erf((nlevels-mean)/sdev)-erf((0-mean)/sdev)+2./4096);
}
void test7()//is factorized joint PDF(x, y)=f(x)g(y) any good?
{
	int nlevels=256;

	double
		meanx=64, sdevx=5,
		meany=138, sdevy=10;
	
	double bpp_joint=0, bpp_separate=0;
	for(int ky=0;ky<nlevels;++ky)
	{
		double p=pdf(ky, meany, sdevy, nlevels);
		if(p)
		{
			double lgp=log2(p);
			if(isfinite(lgp))
				bpp_separate-=p*log2(p);
		}
	}
	for(int kx=0;kx<nlevels;++kx)
	{
		double p=pdf(kx, meanx, sdevx, nlevels);
		if(p)
		{
			double lgp=log2(p);
			if(isfinite(lgp))
				bpp_separate-=p*log2(p);
		}
	}

	for(int ky=0;ky<nlevels;++ky)
	{
		double py=pdf(ky, meany, sdevy, nlevels);
		for(int kx=0;kx<nlevels;++kx)
		{
			double px=pdf(kx, meanx, sdevx, nlevels);
			double p=px*py;
			if(p)
			{
				double lgp=log2(p);
				if(isfinite(lgp))
					bpp_joint-=p*log2(p);
			}
		}
	}
	printf("CR_joint\t%lf\nCR_separate\t%lf\n", 16/bpp_joint, 16/bpp_separate);

	printf("Done.\n");
	pause();
	exit(0);
}
#endif
extern unsigned g_conf;

const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
};
void batch_test(const char *path)
{
	ArrayHandle filenames=get_filenames(path, g_extensions, COUNTOF(g_extensions), 1);
	if(!filenames)
	{
		printf("No images in \"%s\"\n", path);
		return;
	}
	long long
		count_PNG=0, count_JPEG=0,
		sum_cPNGsize=0, sum_cJPEGsize=0,
		sum_uPNGsize=0, sum_uJPEGsize=0,
		sum_testsize=0;
		//sum_test3size[2]={0};
	for(ptrdiff_t k=0;k<(ptrdiff_t)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);

		if(!fn)
		{
			LOG_ERROR("filename read error");
			continue;
		}

		ptrdiff_t formatsize=get_filesize(fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images
			continue;

		int iw=0, ih=0, nch0=0, stride=4;
		long long cycles=__rdtsc();
		unsigned char *buf=stbi_load(fn[0]->data, &iw, &ih, &nch0, stride);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}

		ptrdiff_t res=(ptrdiff_t)iw*ih, len=res*stride, usize=res*nch0;
		double ratio=(double)usize/formatsize;
		printf("\"%s\"\tCR %lf (%lf BPP) Dec %lf CPB", fn[0]->data, ratio, 8/ratio, (double)cycles/usize);
		if(!acme_stricmp(fn[0]->data+fn[0]->count-3, "PNG"))
		{
			sum_cPNGsize+=formatsize;
			sum_uPNGsize+=usize;
			++count_PNG;
		}
		else//assumed
		{
			sum_cJPEGsize+=formatsize;
			sum_uJPEGsize+=usize;
			++count_JPEG;
		}

		unsigned char *b2=(unsigned char*)malloc(len);
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memset(b2, 0, len);

#if 1
		apply_transforms_fwd(buf, iw, ih);

		ArrayHandle cdata=0;
		printf("rANS\n");
		cycles=__rdtsc();
		rans4_encode(buf, (ptrdiff_t)iw*ih, nch0, 4, &cdata, 0);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		
		sum_testsize+=cdata->count;

		cycles=__rdtsc();
		rans4_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, nch0, 4, b2, 0);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "rANS", 0);
		memset(b2, 0, len);
#endif

	//test10
#if 0
		ArrayHandle cdata=0;
		//debug_ptr=buf;
		printf("test10\n");
		cycles=__rdtsc();
		test10_encode(buf, iw, ih, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		
		sum_testsize+=cdata->count;

		cycles=__rdtsc();
		test10_decode(cdata->data, cdata->count, iw, ih, b2);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "test10", 0);
		memset(b2, 0, len);

		printf("\n");
#endif

		//test8 / test9
#if 0
		apply_transforms_fwd(buf, iw, ih);
	//	double csize=test9_estimate_cr(buf, iw, ih, 0);
		double csize=test8_estimate_cr(buf, iw, ih, 0);
		double r2=(double)usize/csize;
		printf("\tCR2 %lf", r2);
		if(r2>ratio)
			printf("\t%lf smaller", r2/ratio);
		printf("\n");
		sum_testsize+=(long long)ceil(csize);
#endif

		//transform selection
#if 0
		ArrayHandle cdata=0;
		colortransform_xgz_fwd(buf, iw, ih);

		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		unsigned char *temp=(unsigned char*)malloc(MAXVAR(iw, ih));
		for(int kc=0;kc<3;++kc)
			dwt2d_cdf53_fwd(buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		array_free(&sizes);
		free(temp);

		//image_pred(buf, iw, ih, nch0, 4);

		cycles=__rdtsc();
		rans4_encode(buf, (ptrdiff_t)iw*ih, nch0, 4, &cdata, 0);
		cycles=__rdtsc()-cycles;
		printf("\tCR %lf Enc %lf CPB", (double)usize/cdata->count, (double)cycles/usize);
		double r2=(double)usize/cdata->count;
		sum_testsize+=cdata->count;
	
		cycles=__rdtsc();
		rans4_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, nch0, 4, b2, 0);
		cycles=__rdtsc()-cycles;
		printf("\tDec CPB %lf ", (double)cycles/usize);
		
		if(r2>ratio)
			printf("\t%lf smaller than current", r2/ratio);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "rANS", 0);
		memset(b2, 0, len);
#endif

		//test3
#if 0
		ArrayHandle cdata=0;
		for(int diff=0;diff<3;++diff)//test3s
		{
			const char *testname=0;
			switch(diff)
			{
			case 0:testname="test3s ";break;
			case 1:testname="test3sd";break;
			case 2:testname="test3sD";break;
			}
			printf("%s\t", testname);

			cycles=__rdtsc();
			test3_encode(buf, iw, ih, nch0, stride, &cdata, diff);
			cycles=__rdtsc()-cycles;
			double r2=(double)usize/cdata->count;
			printf("CR %lf Enc %lf", r2, (double)cycles/usize);
		
			cycles=__rdtsc();
			test3_decode(cdata->data, cdata->count, iw, ih, nch0, stride, b2, diff);
			cycles=__rdtsc()-cycles;
			printf(" Dec %lf CPB", (double)cycles/usize);

			//comment
			if(r2>ratio)
				printf(" %lf smaller than current", r2/ratio);
			sum_test3size[diff]+=cdata->count;

			array_free(&cdata);
			printf("\t");
			compare_bufs_uint8(b2, buf, iw, ih, nch0, stride, testname, 1);
			memset(b2, 0, len);
		}
#endif

		//printf("\n");
		free(buf);
		free(b2);
	}
	ptrdiff_t totalusize=sum_uPNGsize+sum_uJPEGsize;
	if(totalusize)
	{
		printf("\nOn average:\n");
		if(sum_cPNGsize)
			printf("PNG     CR %lf  (%lld images)\n", (double)sum_uPNGsize/sum_cPNGsize, count_PNG);
		if(sum_cJPEGsize)
			printf("JPEG    CR %lf  (%lld images)\n", (double)sum_uJPEGsize/sum_cJPEGsize, count_JPEG);
		printf("test    CR %lf\n", (double)totalusize/sum_testsize);
		//printf("test3s  CR %lf\n", (double)totalusize/sum_test3size[0]);
		//printf("test3sd CR %lf\n", (double)totalusize/sum_test3size[1]);
	}
	else
		printf("\nNo valid images found\n");

	array_free(&filenames);

	printf("\nDone.\n");
	pause();
}
int main(int argc, char **argv)
{
	//int width=10,
	//	n=pyramid_getsize(width);
	//for(int k=0;k<n;++k)
	//	printf("%c", pyramid_getchar(width, k));
	//DCTtest();
	//DCTtest2();
	//test4();
	//test5();
	//test6();
	//test7();

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
		ptrdiff_t formatsize=get_filesize(fn);
		if(formatsize==-1)
		{
			LOG_ERROR("Cannot open \"%s\"", fn);
			return 0;
		}
		if(!formatsize)//path
		{
			batch_test(fn);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		cycles=__rdtsc();
		buf=stbi_load(fn, &iw, &ih, &nch0, nch);
		cycles=__rdtsc()-cycles;
		if(!buf)
			LOG_ERROR("Couldn't open \"%s\"", fn);
		resolution=(size_t)iw*ih;
		len=resolution*nch;

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
		printf("Generating test data (%d bit binomial)...\n", unibits);
		fill_halfbinomial(buf, len, unibits);
#endif
	}

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<len;k+=nch)
			buf[k]=0xFF;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	size_t usize=len*nch0>>2;

	printf("\n");
	
	ArrayHandle cdata=0;
	//const void *ptr, *end;
	
	int loud=0;
	//lowpassfilt(buf, len);//
	//lowpassfilt(buf, len);//
	//cvt2graycode(buf, len);//

	//print_bytes(buf, len);

	//test3
#if 0
	//for(int k=8;k<=256;k<<=1)//8=av(6, 10)
	//for(int k=8;k<=32;++k)		//10 works best
	int diff=0;
	//for(int diff=0;diff<3;++diff)
	{
		//buf0=buf;
		//printf("test3 %3d\t\t", k);
		const char *testname=0;
		switch(diff)
		{
		case 0:testname="test3s ";break;
		case 1:testname="test3sd";break;
		case 2:testname="test3sD";break;
		}
		printf("%s\n", testname);

		cycles=__rdtsc();
		size_t savedbytes=test3_encode(buf, iw, ih, nch0, nch, &cdata, diff);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		
		cycles=__rdtsc();
		test3_decode(cdata->data, cdata->count, iw, ih, nch0, nch, b2, diff);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		test3_printsummary(cdata, savedbytes, usize, nch0);//

		//lodepng_encode_file("out.PNG", b2, iw, ih, LCT_RGBA, 8);//DEBUG SAVE

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, testname, 1);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//test4
#if 0
	{
		printf("test4\n");
		cycles=__rdtsc();
		size_t overhead=test4_encode(buf, iw, ih, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  ratio %lf  (%lld + %lld = %lld bytes)\n", (double)cycles/usize, (double)usize/cdata->count, overhead, cdata->count-overhead, cdata->count);


		array_free(&cdata);
		printf("\n");
	}
#endif

	//save squeeze
#if 0
	{
		int maxdim=MAXVAR(iw, ih);
		size_t res=(size_t)iw*ih, nval=res*3;
		short *b3=(short*)malloc(nval*sizeof(short)), *temp=(short*)malloc(maxdim*sizeof(short));
		for(int k=0, ks=0, kd=0;k<res;++k, ks+=4, kd+=3)//copy image to 16 bit buffer
		{
			b3[kd  ]=buf[ks  ];
			b3[kd+1]=buf[ks+1];
			b3[kd+2]=buf[ks+2];
		}
		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		DWTSize *psizes=(DWTSize*)sizes->data;
		int nsizes=(int)sizes->count;
		squeeze_2d_fwd(b3  , psizes, nsizes, 3, 0, temp);
		squeeze_2d_fwd(b3+1, psizes, nsizes, 3, 0, temp);
		squeeze_2d_fwd(b3+2, psizes, nsizes, 3, 0, temp);
		for(int ks=0;ks<(int)sizes->count-1;++ks)//normalize back to 8 bit buffer (lossy)
		{
			int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
			for(int kc=0;kc<3;++kc)
			{
				for(int kq=1;kq<4;++kq)
				{
					int ox, oy;
					short nbits;
					switch(kq)
					{
					case 1:nbits= 9+(kc<2), ox=px, oy= 0;break;
					case 2:nbits= 9+(kc<2), ox= 0, oy=py;break;
					case 3:nbits=10+(kc<2), ox=px, oy=py;break;
					default:LOG_ERROR("");break;
					}
					for(int ky=oy;ky<dy;++ky)
					{
						for(int kx=ox;kx<dx;++kx)
						{
							int idx=iw*ky+kx;
							b2[idx<<2|kc]=(unsigned char)((b3[3*idx+kc]+(1<<(nbits-1)))>>(nbits-8));
						}
					}
				}
			}
		}
		for(int k=0;k<(int)res;++k)//set alpha
			b2[k<<2|3]=0xFF;
		
		lodepng_encode_file("squeeze.PNG", b2, iw, ih, LCT_RGBA, 8);

		array_free(&sizes);
		free(temp);
		free(b3);
	}
#endif

	//test5
#if 0
	printf("test5\n");

	g_conf=0x8000;
	//for(g_conf=0;g_conf<0x10000;g_conf+=256)
	{
		printf("CONF 0x%04X\n", g_conf);
		cycles=__rdtsc();
		test5_encode(buf, iw, ih, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
	}

	array_free(&cdata);
	printf("\n");
#endif

	//test6
#if 0
	printf("test6\n");
	
	extern int g_offset;
	//int best=0x12900;
	//for(g_offset=best-0x800;g_offset<best+0x800;g_offset+=0x80)
	{
		cycles=__rdtsc();
		test6_encode(buf, iw, ih, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		//printf("Enc %lf CPB  ratio %lf  size %lld\t0x%04X\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count, g_offset);

		array_free(&cdata);
	}

	printf("\n");
#endif

	//Exp4.C33
#if 0
	//printf("Exp4.C33\n");
	printf("Exp4.C33 Initializing...\n");
	codec33_init("D:/Share Box/Python/entropybattle3/C33-20230406-005151.txt");
	
	printf("Exp4.C33 Encoding...\n");
	cycles=__rdtsc();
	codec33_encode(buf, iw, ih, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

	array_free(&cdata);
	printf("Exp4.C33\n");
#endif

	//test7
#if 0
	printf("test7: joint PDF\n");
	cycles=__rdtsc();
	test7_encode(buf, iw, ih, 1, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

	array_free(&cdata);
	printf("\n");
#endif

	//test8
#if 0
	printf("test8\n");
	cycles=__rdtsc();
	test8_encode(buf, iw, ih, 1, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
	
	cycles=__rdtsc();
	test8_decode(cdata->data, cdata->count, iw, ih, 1, b2);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	lodepng_encode_file("kodim21-test8.PNG", b2, iw, ih, LCT_RGBA, 8);

	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "test8", 0);
	memset(b2, 0, len);

	printf("\n");
#endif

	//test10
#if 1
	extern const unsigned char *debug_ptr;
	//debug_ptr=buf;
	printf("test10\n");
	cycles=__rdtsc();
	test10_encode(buf, iw, ih, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  ratio %lf  size %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
	
	cycles=__rdtsc();
	test10_decode(cdata->data, cdata->count, iw, ih, b2);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	lodepng_encode_file("kodim21-test10.PNG", b2, iw, ih, LCT_RGBA, 8);

	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "test10", 0);
	memset(b2, 0, len);

	printf("\n");
#endif

	//predict image
#if 1
	printf("Predict image...\n");
	{
		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		unsigned char *temp=(unsigned char*)malloc(MAXVAR(iw, ih));
		
		addbuf(buf, iw, ih, nch0, nch, 128);//unsigned char -> signed char
		
		//colortransform_ycocg_fwd((char*)buf, iw, ih);
		colortransform_xgz_fwd((char*)buf, iw, ih);
		//colortransform_xyz_fwd((char*)buf, iw, ih);

#if 0
		memcpy(b2, buf, len);
		colortransform_xyz_fwd(b2, iw, ih);
		colortransform_xyz_inv(b2, iw, ih);
		for(int kc=0;kc<3;++kc)
		{
			dwt2d_haar_fwd((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, (char*)temp);
			dwt2d_haar_inv((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, (char*)temp);
		}
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "CT", 0);
		printf("\n");
#endif
		
		for(int kc=0;kc<3;++kc)
			//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
			dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
			//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

		addbuf(buf, iw, ih, nch0, nch, 128);
		//save_CDF53("kodim21-haar-stage", buf, (DWTSize*)sizes->data, (int)sizes->count, 4);//
		//lodepng_encode_file("kodim21-haar.PNG", buf, iw, ih, LCT_RGBA, 8);//

		array_free(&sizes);
		free(temp);
	}
	//squeeze_8bit_lossy(buf, iw, ih, nch0, nch);
//	image_pred(buf, iw, ih, nch0, nch);

	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//
#endif

	//2D differentiation
#if 0
	printf("Differentiating image...\n");

	//memcpy(b2, buf, len);//save copy

	differentiate_image(buf, iw, ih, nch0, nch);//add/sub 128 for pleasant image

	//lodepng_encode_file("out.PNG", buf, iw, ih, LCT_RGBA, 8);//save differentiated image

	//integrate_image(buf, iw, ih, nch0, nch);//check
	//compare_bufs_uint8(buf, b2, len, "Diff", 0);
#endif

	test8_estimate_cr(buf, iw, ih, 1);
	//test9_estimate_cr(buf, iw, ih, 1);

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
	//haar_2d_fwd(buf, iw, ih, nch0, nch, nstages, &b3);
	//
	//for(ptrdiff_t k=0;k<len;++k)//9 bit -> 16 bit
	//{
	//	if(!((k+1)&3))
	//	{
	//		b3[k]=0xFFFF;
	//		b2[k]=0xFF;
	//	}
	//	else
	//	{
	//		b3[k]<<=6;
	//		b3[k]+=0x8000;
	//		b2[k]=(unsigned short)b3[k]>>8;
	//		b3[k]=(unsigned short)b3[k]>>8|b3[k]<<8;//PNG is big endian
	//	}
	//}

	memcpy(b2, buf, len);//save diff

	printf("Saving...\n");
	//lodepng_encode_file("out16.PNG", b3, iw, ih, LCT_RGBA, 16);
	lodepng_encode_file("out8.PNG", b2, iw, ih, LCT_RGBA, 8);
	printf("Done\n");

	//haar_2d_inv(b3, iw, ih, nch0, nch, nstages, &b2);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "Haar", 0);
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
	compare_bufs_uint8(buf, udata->data, iw, ih, nch0, nch, "Huffman", 0);
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
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "AC", 0);
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
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "ABAC", 0);
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
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "ABAC_SSE2", 0);
	memset(b2, 0, len);
	printf("\n");
#endif

	//rANS
#if 1
	printf("rANS\n");
	cycles=__rdtsc();
	rans4_encode(buf, (ptrdiff_t)iw*ih, nch0, nch, &cdata, 0);
	//rans4_encode(buf, len, 1, 0, &cdata, loud, pred);
	cycles=__rdtsc()-cycles;
	//double cr=(double)usize/cdata->count, cr0=8/entropy[4];
	//printf("Enc CPB %lf ratio %lf%s\n", (double)cycles/usize, cr, cr>cr0?" IMPOSSIBLE":"");
	printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	
	cycles=__rdtsc();
	rans4_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, nch0, nch, b2, 0);
	//rans4_decode(cdata->data, cdata->count, len, 1, 0, b2, loud, pred);
	cycles=__rdtsc()-cycles;
	printf("Dec CPB %lf\n", (double)cycles/usize);

	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "rANS", 0);
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
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "rANS_SSE2", 0);
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
	compare_bufs_uint8(b2, buf, resolution, "ABAC_SSE2", 0);
#endif

	//ABAC_SSE2 debug
#if 0
	resolution=64;
	printf("resolution %d\n", resolution);
	abac0a_encode(buf, resolution, 1, &cdata, 0);
	ptr=cdata->data;
	end=cdata->data+cdata->count;
	ptr=abac0a_decode(ptr, end, b2, resolution, 1, 0);
	compare_bufs_uint8(b2, buf, resolution, "ABAC_SSE2", 0);
#endif

	//rANS_SSE2 debug
#if 0
	size_t offset=0, len2=len;
	//size_t offset=len/2, len2=len/2;//coincedence permutation
	//size_t offset=len/4, len2=len*3/4;
	//size_t offset=len-16, len2=16;
	rans_sse2_encode(buf+offset, len2, nch, 0, &cdata);
	rans_sse2_decode(cdata->data, cdata->count, b2+offset, len2, nch, 0, buf+offset);

	compare_bufs_uint8(b2+offset, buf+offset, len2, "rANS_SSE2", 0);
#endif


	free(buf);
	free(b2);
#endif
	printf("Done.\n");
	pause();
	return 0;
}