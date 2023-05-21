#include"battle.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
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
		image_save_png_rgba8(fn, b2, px+xodd, py+yodd);
		//lodepng_encode_file(fn, b2, px+xodd, py+yodd, LCT_RGBA, 8);
	}
	free(b2);
}
void save_channel(unsigned char *buf, int iw, int ih, int stride, int val_offset, const char *format, ...)
{
	va_list args;
	va_start(args, format);
	vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);

	unsigned char *b2=(unsigned char*)malloc((size_t)iw*ih);
	if(!b2)
		return;
	for(int ks=0, kd=0, res=iw*ih;kd<res;ks+=stride, ++kd)
		b2[kd]=buf[ks]+val_offset;
	
	image_save_png_rgba8(g_buf, b2, iw, ih);
	//lodepng_encode_file(g_buf, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}
void save_DWT_int8(const char *name, unsigned char *buf, DWTSize *sizes, int nsizes, int stride)
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
		
		image_save_png_rgba8(fn, b2, px+xodd, py+yodd);
		//lodepng_encode_file(fn, b2, px+xodd, py+yodd, LCT_RGBA, 8);
	}
	free(b2);
}
void save_DWT_int8_all(const char *nametag, unsigned char *buf, DWTSize *sizes, int nsizes)
{
	if(nsizes<2)
		return;
	char fn[256]={0};
	int iw=sizes->w, ih=sizes->h;
	int lsbw=(iw>>1)+(iw&1), lsbh=(ih>>1)+(ih&1), res=lsbw*lsbh;//largest subband size
	unsigned char *b2=(unsigned char*)malloc((size_t)res*4);
	memset(b2, 0xFF, (size_t)res*4);
	for(int k=1;k<nsizes;++k)
	{
		int px=sizes[k].w, py=sizes[k].h, dx=sizes[k-1].w, dy=sizes[k-1].h, xodd=dx&1, yodd=dy&1;
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
			{
				int dstidx=(px+xodd)*ky+kx, srcidx=iw*ky+px+kx;
				b2[dstidx<<2  ]=buf[srcidx<<2  ];
				b2[dstidx<<2|1]=buf[srcidx<<2|1];
				b2[dstidx<<2|2]=buf[srcidx<<2|2];
			}
		}
		snprintf(fn, 128, "%s-iter%02d-TR.PNG", nametag, k);
		image_save_png_rgba8(fn, b2, px+xodd, py);
		//lodepng_encode_file(fn, b2, px+xodd, py, LCT_RGBA, 8);

		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0, xend=dx-px;kx<xend;++kx)
			{
				int dstidx=(px+xodd)*ky+kx, srcidx=iw*(py+ky)+px+kx;
				b2[dstidx<<2  ]=buf[srcidx<<2  ];
				b2[dstidx<<2|1]=buf[srcidx<<2|1];
				b2[dstidx<<2|2]=buf[srcidx<<2|2];
			}
		}
		snprintf(fn, 128, "%s-iter%02d-BR.PNG", nametag, k);
		image_save_png_rgba8(fn, b2, px+xodd, py+yodd);
		//lodepng_encode_file(fn, b2, px+xodd, py+yodd, LCT_RGBA, 8);

		for(int ky=0, yend=dy-py;ky<yend;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				int dstidx=px*ky+kx, srcidx=iw*(py+ky)+kx;
				b2[dstidx<<2  ]=buf[srcidx<<2  ];
				b2[dstidx<<2|1]=buf[srcidx<<2|1];
				b2[dstidx<<2|2]=buf[srcidx<<2|2];
			}
		}
		snprintf(fn, 128, "%s-iter%02d-BL.PNG", nametag, k);
		image_save_png_rgba8(fn, b2, px, py+yodd);
		//lodepng_encode_file(fn, b2, px, py+yodd, LCT_RGBA, 8);
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
		sum+=buf[stride*k];
	double mean=(double)sum/res;
	return mean;
}
double get_var(const unsigned char *buf, int res, int stride, double mean)
{
	double sum=0;
	for(int k=0;k<res;++k)
	{
		double x=buf[stride*k]-mean;
		sum+=x*x;
	}
	sum/=res;
	return sum;
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
	//conf[0]*=2;
	//conf[1]*=2;
	//conf[2]*=2;
	double gain=1/sqrt(8*M_PI*M_PI*M_PI*conf[0]*conf[1]*conf[2]);
	conf[0]=1/sqrt(conf[0]);
	conf[1]=1/sqrt(conf[1]);
	conf[2]=1/sqrt(conf[2]);

	double csize=0;
	for(int kp=0;kp<res;++kp)//Zipf's law
	{
		double
			r=(buf[kp<<2  ]-mean[0])*conf[0],
			g=(buf[kp<<2|1]-mean[1])*conf[1],
			b=(buf[kp<<2|2]-mean[2])*conf[2];
		double prob=exp(-0.5*(r*r+g*g+b*b));
		prob*=gain;
		if(prob)
		{
			prob*=0xFFFFFFFF;
			++prob;
			prob/=0x100000000;
			double contribution=-log2(prob);
			if(isfinite(contribution))
				csize+=contribution;
		}
	}
	csize/=8;
	double end=time_ms();
	int usize=res*3;
	double CR=usize/csize;
	if(loud)
		printf("Test 9 estimate: usize %d csize %lf CR %lf elapsed %lfms\n\n", usize, csize, CR, end-start);
	return csize;
#if 0
	double entropy=0;
	for(int color=0;color<0x1000000;++color)//Shannon's law
	{
		//if(color==0x808080)//
		//	color=0x808080;//
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
#endif
}
double test11_estimate_cr(const unsigned char *buf, int bw, int bh, int loud)
{
	const int hsize=0x1000000*sizeof(unsigned);
	unsigned *hist=(unsigned*)malloc(hsize);
	if(!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, hsize);

	int res=bw*bh;
	for(int k=0;k<res;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		++hist[color];
	}
	double csize=0;
	for(int k=0;k<res;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		double p=(double)hist[color]/res;
		p*=0xFFFFFFFF;
		++p;
		p/=0x100000000;
		double bitsize=-log2(p);
		csize+=bitsize;
	}
	free(hist);
	csize/=8;
	
	int usize=res*3;
	double CR=usize/csize;
	if(loud)
		printf("Test 11 estimate: usize %d csize %lf CR %lf\n\n", usize, csize, CR);
	return csize;
}
double test12_estimate_csize(const unsigned char *buf, int bw, int bh, int blocksize, int loud)
{
	int res=bw*bh;
	int blockcount=res/(blocksize*blocksize);
	double csize=0;
	if(blockcount>1000)
	{
		for(int ky=0;ky<bh;ky+=blocksize)
		{
			for(int kx=0;kx<bw;kx+=blocksize)
			{
				int xsize, ysize, count;

				if(ky+blocksize<=bh)
					ysize=blocksize;
				else
					ysize=bh-ky;

				if(kx+blocksize<=bw)
					xsize=blocksize;
				else
					xsize=bw-kx;
				
				count=xsize*ysize;
				for(int ky2=0;ky2<ysize;++ky2)
				{
					for(int kx2=0;kx2<xsize;++kx2)
					{
						unsigned color=((int*)buf)[bw*(ky+ky2)+kx+kx2]&0xFFFFFF;
						int freq=0;
						for(int ky3=0;ky3<ysize;++ky3)
						{
							for(int kx3=0;kx3<xsize;++kx3)
							{
								unsigned c3=((int*)buf)[bw*(ky+ky3)+kx+kx3]&0xFFFFFF;
								freq+=c3==color;
							}
						}
						double p=(double)freq/count, bitsize=-log2(p);
						csize+=bitsize;
					}
				}
			}
		}
	}
	else
	{
		const int hsize=0x1000000*sizeof(unsigned);
		unsigned *hist=(unsigned*)malloc(hsize);
		if(!hist)
		{
			LOG_ERROR("Allocation error");
			return 0;
		}

		for(int ky=0;ky<bh;ky+=blocksize)
		{
			for(int kx=0;kx<bw;kx+=blocksize)
			{
				int xsize, ysize, count;

				if(ky+blocksize<=bh)
					ysize=blocksize;
				else
					ysize=bh-ky;

				if(kx+blocksize<=bw)
					xsize=blocksize;
				else
					xsize=bw-kx;

				count=xsize*ysize;
				memset(hist, 0, hsize);
				for(int ky2=0;ky2<ysize;++ky2)
				{
					for(int kx2=0;kx2<xsize;++kx2)
					{
						unsigned color=((int*)buf)[bw*(ky+ky2)+kx+kx2]&0xFFFFFF;
						++hist[color];
					}
				}
				for(int ky2=0;ky2<ysize;++ky2)
				{
					for(int kx2=0;kx2<xsize;++kx2)
					{
						unsigned color=((int*)buf)[bw*(ky+ky2)+kx+kx2]&0xFFFFFF;
						double p=(double)hist[color]/count, bitsize=-log2(p);
						csize+=bitsize;
					}
				}
			}
		}
		free(hist);
	}
	csize/=8;
	
	int usize=res*3;
	double CR=usize/csize;
	if(loud)
		printf("Test 12 estimate: bsize %4d usize %7d csize %14lf CR %9lf\n", blocksize, usize, csize, CR);
	return csize;
}
double test13_estimate_csize(const unsigned char *buf, int bw, int bh, int blocksize, int loud)
{
	int res=bw*bh;
	int blockcount=res/(blocksize*blocksize);
	double csize=0;
	unsigned *block=(unsigned*)malloc((size_t)blocksize*blocksize*sizeof(unsigned));
	for(int ky=0;ky<bh;ky+=blocksize)
	{
		for(int kx=0;kx<bw;kx+=blocksize)
		{
			int xsize, ysize, count;

			if(ky+blocksize<=bh)
				ysize=blocksize;
			else
				ysize=bh-ky;

			if(kx+blocksize<=bw)
				xsize=blocksize;
			else
				xsize=bw-kx;

			count=xsize*ysize;
			for(int ky2=0;ky2<ysize;++ky2)
				memcpy(block+xsize*ky2, (int*)buf+bw*(ky+ky2)+kx, xsize*sizeof(unsigned));
			double mean[3]=
			{
				get_mean((unsigned char*)block  , count, 4),
				get_mean((unsigned char*)block+1, count, 4),
				get_mean((unsigned char*)block+2, count, 4),
			};
			double conf[]=
			{
				get_var((unsigned char*)block  , count, 4, mean[0]),
				get_var((unsigned char*)block+1, count, 4, mean[1]),
				get_var((unsigned char*)block+2, count, 4, mean[2]),
			};
			double gain=1/sqrt(8*M_PI*M_PI*M_PI*conf[0]*conf[1]*conf[2]);
			conf[0]=1/sqrt(conf[0]);
			conf[1]=1/sqrt(conf[1]);
			conf[2]=1/sqrt(conf[2]);
			for(int k=0;k<count;++k)
			{
				unsigned char *color=(unsigned char*)(block+k);
				double
					r=(color[0]-mean[0])*conf[0],
					g=(color[1]-mean[1])*conf[1],
					b=(color[2]-mean[2])*conf[2];
				double prob=exp(-0.5*(r*r+b*b+g*g));
				prob*=gain;
				if(prob)
				{
					prob*=0xFFFFFFFF;
					++prob;
					prob/=0x100000000;
					double contribution=-log2(prob);
					if(isfinite(contribution))
						csize+=contribution;
				}
			}
		}
	}
	free(block);
	csize/=8;
	
	int usize=res*3;
	double CR=usize/csize;
	//double CR=usize/(csize+24*blockcount);
	if(loud)
		printf("Test 13 estimate: bsize %3d usize %d csize %14lf overhead %7d CR %lf\n", blocksize, usize, csize, 24*blockcount, CR);
	return csize;
}
double test15_estimate_csize(const unsigned char *buf, int bw, int bh, int blocksize, int loud)//per-channel per-block histogram
{
	int res=bw*bh;
	int blockcount=res/(blocksize*blocksize);
	double csize[3]={0};
	if(blockcount>1000)
	{
		for(int kc=0;kc<3;++kc)//for each channel
		{
			for(int ky=0;ky<bh;ky+=blocksize)
			{
				for(int kx=0;kx<bw;kx+=blocksize)//for each block
				{
					int xsize, ysize, count;

					if(ky+blocksize<=bh)
						ysize=blocksize;
					else
						ysize=bh-ky;

					if(kx+blocksize<=bw)
						xsize=blocksize;
					else
						xsize=bw-kx;
				
					count=xsize*ysize;
					for(int ky2=0;ky2<ysize;++ky2)
					{
						for(int kx2=0;kx2<xsize;++kx2)//for each pixel
						{
							unsigned char val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
							int freq=0;
							for(int ky3=0;ky3<ysize;++ky3)
							{
								for(int kx3=0;kx3<xsize;++kx3)
								{
									unsigned v3=buf[(bw*(ky+ky3)+kx+kx3)<<2|kc];
									freq+=v3==val;
								}
							}
							double p=(double)freq/count, bitsize=-log2(p);
							csize[kc]+=bitsize;
						}
					}
				}
			}
		}
	}
	else
	{
		const int hsize=256*sizeof(unsigned);
		unsigned *hist=(unsigned*)malloc(hsize);
		if(!hist)
		{
			LOG_ERROR("Allocation error");
			return 0;
		}

		for(int kc=0;kc<3;++kc)
		{
			for(int ky=0;ky<bh;ky+=blocksize)
			{
				for(int kx=0;kx<bw;kx+=blocksize)
				{
					int xsize, ysize, count;

					if(ky+blocksize<=bh)
						ysize=blocksize;
					else
						ysize=bh-ky;

					if(kx+blocksize<=bw)
						xsize=blocksize;
					else
						xsize=bw-kx;

					count=xsize*ysize;
					memset(hist, 0, hsize);
					for(int ky2=0;ky2<ysize;++ky2)
					{
						for(int kx2=0;kx2<xsize;++kx2)
						{
							unsigned val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
							++hist[val];
						}
					}
					for(int ky2=0;ky2<ysize;++ky2)
					{
						for(int kx2=0;kx2<xsize;++kx2)
						{
							unsigned val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
							double p=(double)hist[val]/count, bitsize=-log2(p);
							csize[kc]+=bitsize;
						}
					}
				}
			}
		}
		free(hist);
	}
	csize[0]/=8;
	csize[1]/=8;
	csize[2]/=8;
	double csize_total=csize[0]+csize[1]+csize[2];
	
	if(loud)
	{
		double CR[3]={res/csize[0], res/csize[1], res/csize[2]}, CR_total=res*3/csize_total;
		printf("bsize %4d csize TRGB %14lf %14lf %14lf %14lf CR TRGB %9lf %9lf %9lf %9lf\n", blocksize, csize_total, csize[0], csize[1], csize[2], CR_total, CR[0], CR[1], CR[2]);
	}
	return csize_total;
}
double test16_estimate_csize(const unsigned char *buf, int bw, int bh, int blocksize, float alpha, int loud)
{
	int res=bw*bh;
	int blockcount=res/(blocksize*blocksize);
	double csize[3]={0};
	const int hsize=256*sizeof(unsigned);
	unsigned *h0=(unsigned*)malloc((size_t)hsize*3);
	if(!h0)
	{
		LOG_ERROR("Allcation error");
		return 0;
	}
	memset(h0, 0, (size_t)hsize*3);
	for(int kc=0;kc<3;++kc)//for each channel
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++h0[kc<<8|sym];
		}
	}
	if(blockcount>1000)
	{
		for(int kc=0;kc<3;++kc)//for each channel
		{
			for(int ky=0;ky<bh;ky+=blocksize)
			{
				for(int kx=0;kx<bw;kx+=blocksize)//for each block
				{
					int xsize, ysize, count;

					if(ky+blocksize<=bh)
						ysize=blocksize;
					else
						ysize=bh-ky;

					if(kx+blocksize<=bw)
						xsize=blocksize;
					else
						xsize=bw-kx;
					count=xsize*ysize;

					int xoffset, yoffset;
					if(kx)
						xoffset=blocksize, yoffset=0;
					else if(ky)
						xoffset=0, yoffset=blocksize;
					else
						xoffset=0, yoffset=0;
					for(int ky2=0;ky2<ysize;++ky2)
					{
						for(int kx2=0;kx2<xsize;++kx2)//for each pixel
						{
							unsigned char val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
							double p=0, p0=0;
							if(kx||ky)
							{
								int freq=0;
								for(int ky3=0;ky3<ysize;++ky3)
								{
									for(int kx3=0;kx3<xsize;++kx3)
									{
										unsigned v3=buf[(bw*(ky-yoffset+ky3)+kx-xoffset+kx3)<<2|kc];
										freq+=v3==val;
									}
								}
								p0=(double)h0[kc<<8|val]/res, p=(double)freq/count;
								p*=0xFFFFFFFF;
								++p;
								p/=0x100000000;
								p=p0+(p-p0)*alpha;
							}
							else//first block: use h0
								p=(double)h0[kc<<8|val]/res;
							double bitsize=-log2(p);
							if(!isfinite(bitsize))
								LOG_ERROR("Infinite size");
							csize[kc]+=bitsize;
						}
					}
				}
			}
		}
	}
	else
	{
		unsigned *hist=(unsigned*)malloc(hsize);
		if(!hist)
		{
			LOG_ERROR("Allocation error");
			return 0;
		}

		for(int kc=0;kc<3;++kc)
		{
			for(int ky=0;ky<bh;ky+=blocksize)
			{
				for(int kx=0;kx<bw;kx+=blocksize)
				{
					int xsize, ysize, count;

					if(ky+blocksize<=bh)
						ysize=blocksize;
					else
						ysize=bh-ky;

					if(kx+blocksize<=bw)
						xsize=blocksize;
					else
						xsize=bw-kx;

					count=xsize*ysize;
					if(kx||ky)//bypass first block
					{
						int xoffset, yoffset;
						if(kx)
							xoffset=blocksize, yoffset=0;
						else
							xoffset=0, yoffset=blocksize;
						memset(hist, 0, hsize);
						for(int ky2=0;ky2<ysize;++ky2)
						{
							for(int kx2=0;kx2<xsize;++kx2)
							{
								unsigned val=buf[(bw*(ky-yoffset+ky2)+kx-xoffset+kx2)<<2|kc];
								++hist[val];
							}
						}
						for(int ky2=0;ky2<ysize;++ky2)
						{
							for(int kx2=0;kx2<xsize;++kx2)
							{
								unsigned val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
								double p0=(double)h0[kc<<8|val]/res, p=(double)hist[val]/count;
								p*=0xFFFFFFFF;
								++p;
								p/=0x100000000;
								p=p0+(p-p0)*alpha;
								double bitsize=-log2(p);
								if(!isfinite(bitsize))
									LOG_ERROR("Infinite size");
								csize[kc]+=bitsize;
							}
						}
					}
					else//first block: use h0
					{
						for(int ky2=0;ky2<ysize;++ky2)
						{
							for(int kx2=0;kx2<xsize;++kx2)//for each pixel
							{
								unsigned char val=buf[(bw*(ky+ky2)+kx+kx2)<<2|kc];
								double p=(double)h0[kc<<8|val]/res;
								double bitsize=-log2(p);
								if(!isfinite(bitsize))
									LOG_ERROR("Infinite size");
								csize[kc]+=bitsize;
							}
						}
					}
				}
			}
		}
		free(hist);
	}
	free(h0);
	csize[0]/=8;
	csize[1]/=8;
	csize[2]/=8;
	double csize_total=csize[0]+csize[1]+csize[2];
	
	if(loud)
	{
		double CR[3]={res/csize[0], res/csize[1], res/csize[2]}, CR_total=res*3/csize_total;
		printf("block %4d %5.2f%% csize TRGB %14lf %14lf %14lf %14lf CR TRGB %9lf %9lf %9lf %9lf\n", blocksize, alpha*100, csize_total, csize[0], csize[1], csize[2], CR_total, CR[0], CR[1], CR[2]);
	}
	return csize_total;
}
#if 0
int t20_accblock(const unsigned char *buf, int bw, int x1, int x2, int y1, int y2, int xstep, int ystep, int attenuation)
{
	unsigned long long sdev=0;
	y2*=ystep;
	x2*=xstep;
	int normal=0;
	for(int ky=y1;ky*ystep<y2;ky+=ystep)
	{
		int weight=0x10000;
		for(int kx=x1;kx*xstep<x2;kx+=xstep)
		{
			char val=buf[bw*ky+kx]-128;
			sdev+=val*val;
		}
	}

}
double test20_estimate_csize(const unsigned char *buf, int bw, int bh, int attenuation, int mixfactor, int blocksize, int margin, int loud)
{
	const int hsize=256*sizeof(unsigned);
	unsigned *hist=(unsigned*)malloc(hsize);
	if(!hist)
	{
		LOG_ERROR("Allcation error");
		return 0;
	}
	int res=bw*bh;
	for(int kc=0;kc<3;++kc)
	{
		memset(hist, 0, hsize);
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++hist[sym];
		}
		for(int ky=0;ky<bh;ky+=blocksize)
		{
			for(int kx=0;kx<bw;kx+=blocksize)//for each block
			{
				int xsize, ysize, count;

				if(ky+blocksize<=bh)
					ysize=blocksize;
				else
					ysize=bh-ky;

				if(kx+blocksize<=bw)
					xsize=blocksize;
				else
					xsize=bw-kx;
				count=xsize*ysize;

				int left=kx-margin, top=ky-margin, right=kx+xsize+margin;
				if(left<0)
					left=0;
				if(top<0)
					top=0;
				if(right>bw)
					right=bw;
				int xend=kx+xsize, yend=ky+ysize;
				if(left<kx)
				{
					for(int ky2=ky;ky2<yend;++ky2)
					{
						int factor=0x10000;
						for(int kx2=kx-1;kx2>=left;--kx2)
						{
						}
					}
				}
				for(int ky2=ky;ky2<yend;++ky2)
				{
					for(int kx2=kx;kx2<xend;++kx2)//for each pixel
					{
					}
				}
			}
		}
	}
}
#endif

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
	unsigned char *buf=image_load(fn, &iw, &ih);
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
void print_ps(float *buf, int count)
{
	const int amplitude=32;
	float vmin=0;
	for(int k=0;k<count;++k)
	{
		if(!k||vmin>buf[k])
			vmin=buf[k];
	}
	
	int c0=(int)round(vmin*amplitude);
	c0=MINVAR(c0, -amplitude);
	for(int k=0;k<count;++k)
	{
		printf("%10.6lf ", buf[k]);
		int c2=(int)roundf(buf[k]*amplitude);
		for(int k2=c0, end=MINVAR(c2, 0);k2<end;++k2)
			printf(" ");
		for(int k2=MINVAR(c2, 0);k2<0;++k2)
			printf("-");
		for(int k2=0;k2<c2;++k2)
			printf("+");
		printf("\n");
	}
	printf("\n");
}
void permute_ps(float *data, int *permutation, int count, float *temp)
{
	memcpy(temp, data, count*sizeof(float));
	for(int k=0;k<count;++k)
		data[permutation[k]]=temp[k];
}
void butterfly_ps(float *a, float *b)
{
	float temp=*a;
	*a+=*b;
	*b=temp-*b;
}
#if 0
void test8()
{
	FCT1D_PS_Params dct8;
	FCT1D_ps_gen(3, &dct8);
	float data[8]={0}, x[8];
	int permutation[]={0, 4, 6, 2, 7, 5, 3, 1};
	float scale[]=
	{
		sinf((float)(M_PI/4))*0.5f,
		sinf((float)(M_PI/4)),
		sinf((float)(M_PI*3/8))*0.5f,
		1/(2*sinf((float)(M_PI*3/8))),
		sinf((float)(M_PI*7/16))*0.5f,
		cosf((float)(M_PI*3/16)),
		1/(2*cosf((float)(M_PI*3/16))),
		1/(2*sinf((float)(M_PI*7/16))),
	};

	float
		p1=13.f/32,//binDCT-C7
		u1=11.f/32,
		p2=22.f/32,
		u2=15.f/32,
		p3= 6.f/32,
		u3= 6.f/32,
		p4=13.f/32,
		u4=22.f/32,
		p5=15.f/32;
	for(int it=0;;it=(it+1)&7)
	{
		memset(data, 0, sizeof(data));
		data[it]=1;
		//srand((unsigned)__rdtsc());
		//for(int k=0;k<8;++k)//initialize
		//	data[k]=(double)rand()/RAND_MAX;
		memcpy(x, data, sizeof(data));

		FCT1D_ps_fwd(&dct8, data, 1);
		
		printf("Data:\n");
		print_ps(x, 8);

		printf("DCT8-II:\n");
		print_ps(data, 8);

		//Chen's factorization
#if 0
		butterfly_ps(x, x+7);
		butterfly_ps(x+1, x+6);
		butterfly_ps(x+2, x+5);
		butterfly_ps(x+3, x+4);

		butterfly_ps(x, x+3);
		butterfly_ps(x+1, x+2);
		
		x[1]=x[0]-x[1], x[0]-=x[1]*0.5f;

		x[2]=p1*x[3]-x[2], x[3]-=u1*x[2];

		x[5]-=p4*x[6], x[6]+=u4*x[5], x[5]=p5*x[6]-x[5];

		x[5]=x[4]-x[5], x[4]-=x[5]*0.5f;
		x[6]=x[7]-x[6], x[7]-=x[6]*0.5f;

		x[4]=p3*x[7]-x[4], x[7]-=u3*x[4];

		x[5]+=p2*x[6], x[6]-=u2*x[5];
#endif

#if 1
		x[7]=x[0]-x[7], x[0]-=x[7]*0.5f;
		x[6]=x[1]-x[6], x[1]-=x[6]*0.5f;
		x[5]=x[2]-x[5], x[2]-=x[5]*0.5f;
		x[4]=x[3]-x[4], x[3]-=x[4]*0.5f;

		x[3]=x[0]-x[3], x[0]-=x[3]*0.5f;
		x[2]=x[1]-x[2], x[1]-=x[2]*0.5f;

		x[1]=x[0]-x[1], x[0]-=x[1]*0.5f;

		x[2]=p1*x[3]-x[2], x[3]-=u1*x[2];

		x[5]-=p4*x[6], x[6]+=u4*x[5], x[5]=p5*x[6]-x[5];

		x[5]=x[4]-x[5], x[4]-=x[5]*0.5f;
		x[6]=x[7]-x[6], x[7]-=x[6]*0.5f;

		x[4]=p3*x[7]-x[4], x[7]-=u3*x[4];

		x[5]+=p2*x[6], x[6]-=u2*x[5];
#endif

		x[0]*=scale[0];
		x[1]*=scale[1];
		x[2]*=scale[2];
		x[3]*=scale[3];
		x[4]*=scale[4];
		x[5]*=scale[5];
		x[6]*=scale[6];
		x[7]*=scale[7];

		//for(int k=0;k<8;++k)
		//	x[k]=(float)permutation[k];

		permute_ps(x, permutation, 8, data);

		printf("Lifting:\n");
		print_ps(x, 8);

		printf("\n\n");

		_getch();
	}
	FCT1D_ps_free(&dct8);
}
#endif
#if 0
void print_pd(double *buf, int count)
{
	const int amplitude=32;
	for(int k=0;k<count;++k)
	{
		printf("%10.6lf ", buf[k]);
		int c2=(int)round(buf[k]*amplitude);
		for(int k2=-amplitude, end=MINVAR(c2, 0);k2<end;++k2)
			printf(" ");
		for(int k2=MINVAR(c2, 0);k2<0;++k2)
			printf("-");
		for(int k2=0;k2<c2;++k2)
			printf("+");
		printf("\n");
	}
	printf("\n");
}
#define DCTSIZE 4
void test8()//useless
{
	float a=cosf((float)(M_PI/4/2)), b=cosf((float)(M_PI/4*3/2)), c=cosf((float)(M_PI/4));
	float DCT[]=
	{
		1, 1, 1, 1,
		a, b, -b, -a,
		c, -c, c, -c,
		b, a, -a, -b,
	};
	float in[]={1.1, 0.9, 0.5, 0.6}, b1[DCTSIZE]={0}, x[DCTSIZE]={0}, temp[DCTSIZE];
	int permutation[]={0, 2, 3, 1};

	float p=13./32, u=11./32;
	for(int it=0;;it=(it+1)&3)
	{
		memset(in, 0, sizeof(in));
		in[it]=1;

#if 0
		printf("QWE: increment P\n");
		printf("ASD: decrement P\n");
		printf("RTY: increment U\n");
		printf("FGH: decrement U\n");
		printf("Z: reset P & U\n");
		printf("Enter: Randomize input\n");
		printf("...: pass\n");
		char c=_getch();
		switch(c)
		{
		case 'Q':case 'q':p+=0.1;  break;
		case 'A':case 'a':p-=0.1;  break;
		case 'W':case 'w':p+=0.01; break;
		case 'S':case 's':p-=0.01; break;
		case 'E':case 'e':p+=0.001;break;
		case 'D':case 'd':p-=0.001;break;
		case 'R':case 'r':u+=0.1;  break;
		case 'F':case 'f':u-=0.1;  break;
		case 'T':case 't':u+=0.01; break;
		case 'G':case 'g':u-=0.01; break;
		case 'Y':case 'y':u+=0.001;break;
		case 'H':case 'h':u-=0.001;break;
		case 'Z':case 'z':p=13./32, u=11./32;break;
		case '\r':case '\n':
			srand((unsigned)__rdtsc());
			for(int k=0;k<DCTSIZE;++k)//initialize
				in[k]=(double)rand()/RAND_MAX;
			break;
		}
		printf("p = %lf\n", p);
		printf("u = %lf\n", u);
#endif
		
		for(int k=0;k<DCTSIZE;++k)//correct DCT
			b1[k]=DCT[k<<2]*in[0]+DCT[k<<2|1]*in[1]+DCT[k<<2|2]*in[2]+DCT[k<<2|3]*in[3];

		memcpy(x, in, sizeof(in));
		x[3]=x[0]-x[3], x[0]-=x[3]*0.5f;
		x[2]=x[1]-x[2], x[1]-=x[2]*0.5f;
		x[1]=x[0]-x[1], x[0]-=x[1]*0.5f;
		x[2]=x[3]*p-x[2], x[3]-=x[2]*u;

		//x[3]-=x[0], x[0]+=x[3]*0.5;
		//x[2]-=x[1], x[1]+=x[2]*0.5;
		//x[1]-=x[0], x[0]+=x[1]*0.5;
		//x[2]-=x[3]*p, x[3]-=x[2]*u;

		permute_ps(x, permutation, 4, temp);

		printf("Input:\n");
		print_ps(in, DCTSIZE);
		printf("DCT:\n");
		print_ps(b1, DCTSIZE);
		printf("Lifting:\n");
		print_ps(x, DCTSIZE);
		printf("\n\n");

		_getch();
	}

	printf("Done.\n");
	pause();
	exit(0);
}
#endif
//extern unsigned g_conf;
#if 0
void test9()
{
	for(int topleft=0;topleft<256;topleft+=16)
	{
		for(int top=0;top<256;top+=16)
		{
			for(int left=0;left<256;left+=16)
			{
				int pred;

				int vmin, vmax;
				if(top<left)
					vmin=top, vmax=left;
				else
					vmin=left, vmax=top;

				if(topleft<vmin)
					pred=vmax;
				else if(topleft>vmax)
					pred=vmin;
				else
					pred=top+left-topleft;

				printf("%3d ", pred);
				//printf("%02X ", pred);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");

	printf("Done.\n");
	pause();
	exit(0);
}
#endif
#if 0
typedef struct E10_CaseStruct
{
	int idx;
	float prob;
} E10_Case;
int e10_cmp_case(const void *p1, const void *p2)
{
	E10_Case const *c1, *c2;

	c1=(E10_Case const*)p1;
	c2=(E10_Case const*)p2;
	return (c1->prob<c2->prob)-(c1->prob>c2->prob);//descending order
}
long long e10_totalpixelcount=0, e10_constcount=0;
int e10_reset(long long *hist)
{
	int histlen=(729LL*9+729LL*3)*sizeof(long long);//9 histograms with 729 cases + 729 sets of means for 3 gradients
	if(hist)
		memset(hist, 0, histlen);
	return histlen;
}
long long* e10_start()
{
	int histlen=e10_reset(0);
	long long *hist=(long long*)malloc(histlen);
	if(!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	e10_reset(hist);
	return hist;
}
void e10_iter(unsigned char *buf, int iw, int ih, int kc, long long *hist)
{
	long long *grad=hist+729LL*9;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char
				topleft =kx-1>=0&&ky-1>=0?buf[(idx-iw-1)<<2|kc]-128:0,
				top     =         ky-1>=0?buf[(idx-iw  )<<2|kc]-128:0,
				topright=kx+1<iw&&ky-1>=0?buf[(idx-iw+1)<<2|kc]-128:0,
				left    =kx-1>=0         ?buf[(idx   -1)<<2|kc]-128:0,
				curr    =                 buf[ idx      <<2|kc]-128  ;
			char nb[]={topleft, top, topright, left};

			int permutation=0;
			char temp;

#define SORT_STEP(A, B)\
			if(nb[A]<nb[B])\
				permutation+=0;\
			else if(nb[A]>nb[B])\
				permutation+=1, temp=nb[A], nb[A]=nb[B], nb[B]=temp;\
			else\
				permutation+=2;

			SORT_STEP(0, 1);
			permutation*=3;
			SORT_STEP(0, 2);
			permutation*=3;
			SORT_STEP(0, 3);

			permutation*=3;
			SORT_STEP(1, 2);
			permutation*=3;
			SORT_STEP(1, 3);

			permutation*=3;
			SORT_STEP(2, 3);
#undef  SORT_STEP

			long long *hk=hist+permutation*9, *gk=grad+permutation*3;
				 if(curr<nb[0])		++hk[0];
			else if(curr==nb[0])	++hk[1];
			else if(curr<nb[1])		++hk[2];
			else if(curr==nb[1])	++hk[3];
			else if(curr<nb[2])		++hk[4];
			else if(curr==nb[2])	++hk[5];
			else if(curr<nb[3])		++hk[6];
			else if(curr==nb[3])	++hk[7];
			else					++hk[8];

			gk[0]+=nb[1]-nb[0];
			gk[1]+=nb[2]-nb[1];
			gk[2]+=nb[3]-nb[2];
			e10_constcount+=nb[0]==nb[1]&&nb[0]==nb[2]&&nb[0]==nb[3];
			++e10_totalpixelcount;
		}
	}
}
double e10_estimate_csize(unsigned char *buf, int iw, int ih, int kc, long long *hist, int loud)//no transforms
{
	int *CDF=(int*)malloc(256*sizeof(int));
	if(!CDF)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	double csize=0;

	double csize_cases=0, csize_rem=0;
	double usize_cases=0, usize_rem=0;
	int casehist[9]={0};

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char
				topleft =kx-1>=0&&ky-1>=0?buf[(idx-iw-1)<<2|kc]-128:0,
				top     =         ky-1>=0?buf[(idx-iw  )<<2|kc]-128:0,
				topright=kx+1<iw&&ky-1>=0?buf[(idx-iw+1)<<2|kc]-128:0,
				left    =kx-1>=0         ?buf[(idx   -1)<<2|kc]-128:0,
				curr    =                 buf[ idx      <<2|kc]-128  ;
			char nb[]={topleft, top, topright, left};

			int permutation=0;
			char temp;

#define SORT_STEP(A, B)\
			if(nb[A]<nb[B])\
				permutation+=0;\
			else if(nb[A]>nb[B])\
				permutation+=1, temp=nb[A], nb[A]=nb[B], nb[B]=temp;\
			else\
				permutation+=2;

			SORT_STEP(0, 1);
			permutation*=3;
			SORT_STEP(0, 2);
			permutation*=3;
			SORT_STEP(0, 3);

			permutation*=3;
			SORT_STEP(1, 2);
			permutation*=3;
			SORT_STEP(1, 3);

			permutation*=3;
			SORT_STEP(2, 3);
#undef  SORT_STEP
			long long *hk=hist+permutation*9;
			long long freq=0, sum=0, sum_hk=0;

			long long intervals[]=
			{
				nb[0] - -128,
				1,
				nb[1] - (nb[0]+1),
				1,
				nb[2] - (nb[1]+1),
				1,
				nb[3] - (nb[2]+1),
				1,
				128 - (nb[3]+1),
			};
			int PDF[9]={0};
			for(int k=0;k<9;++k)
			{
				if(intervals[k])
				{
					sum_hk+=hk[k];
					PDF[k]=hk[k]/intervals[k];
					PDF[k]+=!PDF[k];
					sum+=PDF[k];
				}
			}

			int caseidx;
				 if(curr<nb[0])		caseidx=0;
			else if(curr==nb[0])	caseidx=1;
			else if(curr<nb[1])		caseidx=2;
			else if(curr==nb[1])	caseidx=3;
			else if(curr<nb[2])		caseidx=4;
			else if(curr==nb[2])	caseidx=5;
			else if(curr<nb[3])		caseidx=6;
			else if(curr==nb[3])	caseidx=7;
			else					caseidx=8;
			++casehist[caseidx];
			freq=PDF[caseidx];
			double p=(double)freq/sum, bitsize=-log2(p);
			csize+=bitsize;

			p=(double)hk[caseidx]/sum_hk, bitsize=-log2(p);
			csize_cases+=bitsize;
			p=1./9, bitsize=-log2(p);
			usize_cases+=bitsize;
			if(!(caseidx&1)&&intervals[caseidx]>1)
			{
				p=1./intervals[caseidx], bitsize=-log2(p);
				csize_rem+=bitsize;
			}
		}
	}
	//for(int k=0;k<9;++k)
	//	printf("\tcase %d %7d\n", k, casehist[k]);
	csize/=8;
	csize_cases/=8;
	usize_cases/=8;
	csize_rem/=8;
	if(loud)
	{
		printf("Channel %d estimates:\n", kc);
		printf("naive %14lf  CR %lf\n", csize, iw*ih/csize);
		printf("cases %14lf  CR %lf  U %lf\n", csize_cases, usize_cases/csize_cases, usize_cases);
		printf("rem   %14lf  bypass\n", csize_rem);
		printf("total %14lf  CR %lf\n", csize_cases+csize_rem, iw*ih/(csize_cases+csize_rem));
		printf("\n");
	}
	return csize;
}
long long *e10_sum=0;
int e10_cmp_priority(const void *p1, const void *p2)
{
	int idx1=*(const int*)p1, idx2=*(const int*)p2;
	return (e10_sum[idx1]<e10_sum[idx2])-(e10_sum[idx1]>e10_sum[idx2]);//descending order
}
void e10_print(long long *hist)
{
	#define SORT_BY_PRIORITY
	
	long long *grad=hist+729LL*9;
	int ncases=0;
	long long totalsum=0;
#ifdef SORT_BY_PRIORITY
	long long *sum=(long long*)malloc(729*sizeof(long long));
	int *priority=(int*)malloc(729*sizeof(int));
	if(!sum||!priority)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(sum, 0, 729*sizeof(long long));
	for(int k=0;k<729*9;++k)
		sum[k/9]+=hist[k];
	for(int k=0;k<729;++k)
	{
		priority[k]=k;
		totalsum+=sum[k];
	}
	e10_sum=sum;
	isort(priority, 729, sizeof(int), e10_cmp_priority);
#else
	for(int k=0;k<729*9;++k)
		totalsum+=hist[k];
#endif
	printf("0123 == topleft, top, topright, left\n");
	printf("permutation  order  weight  x<a  x==a  a<x<b  x==b  b<x<c  x==c  c<x<d  x==d  x>d\n\n");
	//printf("               <a           ==a          <b           ==b          <c           ==c          <d           ==d          >d          weight\n\n");
	for(int kp=0;kp<729;++kp)
	{
#ifdef SORT_BY_PRIORITY
		int idx=priority[kp];
#else
		int idx=kp;
#endif
		long long *hk=hist+9*idx;
		char arr[]={'0', '1', '2', '3', 0};
		//char arr[]={'C', 'T', 'R', 'L', 0};

		int permutation=idx, outcome;
		char temp;

#define UNSORT_STEP(A, B)\
			 if(outcome==0)		;\
		else if(outcome==1)		temp=arr[A], arr[A]=arr[B], arr[B]=temp;\
		else					arr[B]=arr[A];
		
		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(2, 3);

		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(1, 3);
		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(1, 2);

		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(0, 3);
		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(0, 2);
		outcome=permutation%3, permutation/=3;
		UNSORT_STEP(0, 1);
#undef  UNSORT_STEP
		
		long long sum=0, kmax=0;
		for(int kr=0;kr<9;++kr)
		{
			sum+=hk[kr];
			if(hk[kmax]<hk[kr])
				kmax=kr;
		}
		
		//if(sum*100./totalsum>0.3)
		if(sum)
		{
			float per[9];
			int kmax=0;
			printf("[%3d] %s %5.2lf  ", idx, arr, sum*100./totalsum);
			for(int kr=0;kr<9;++kr)
			{
				per[kr]=hk[kr]*100.f/sum;
				if(per[kmax]<per[kr])
					kmax=kr;
				if(per[kr])
					printf(" %2.0f", per[kr]);
				else
					printf("   ");
			}
			printf("  -128");
			int start=0, end;
			float CDF=0;
			for(int kr=0;kr<9;++kr)
			{
				CDF+=per[kr];
				end=(int)roundf(CDF);
				char c=kr&1?'a'+(kr>>1):'-';
				//switch(kr)
				//{
				//case 0:c='-';break;
				//case 1:c='a';break;
				//case 2:c='-';break;
				//case 3:c='b';break;
				//case 4:c='-';break;
				//case 5:c='c';break;
				//case 6:c='-';break;
				//case 7:c='d';break;
				//case 8:c='-';break;
				//}
				//char c=kr==kmax?'^':(kr<9-2&&kr&1?'0':'.');

				int e2=start+(end-start+1)/2;
				for(int k2=start;k2<e2;++k2)
					printf("%c", c);
				if(kr>=2&&kr<9-1&&!(kr&1))
				{
					double interval=(double)grad[idx*3+((kr-2)>>1)]/sum;
					if(interval)
						printf("%5.2lf", interval);
					else
						printf("     ");
				}
				for(int k2=e2;k2<end;++k2)
					printf("%c", c);

				start=end;
				if(kr<9-1)
					printf("%c", (per[kr]?'A':'a')+(kr>>1));
				else
					printf("127");
			}
			printf("  [%3d]\n", idx);
			//long long range=grad[idx*3]+grad[idx*3+1]+grad[idx*3+2];
			//printf(" %4.2lf %4.2lf %4.2lf", (double)grad[idx*3]/sum, (double)grad[idx*3+1]/sum, (double)grad[idx*3+2]/sum);
			//printf("\n");

#if 0
			E10_Case h2[9];
			for(int kr=0;kr<9;++kr)
			{
				h2[kr].idx=kr;
				h2[kr].prob=hk[kr]*100.f/sum;
			}
			isort(h2, 9, sizeof(E10_Case), e10_cmp_case);
			printf("[%3d] %s %9lf%%  ", idx, arr, sum*100./totalsum);
			for(int kr=0;kr<9&&h2[kr].prob>0.5f;++kr)
			{
				const char *a=0;
				switch(h2[kr].idx)
				{
				case 0:a="x<a  ";break;
				case 1:a="x==a ";break;
				case 2:a="a<x<b";break;
				case 3:a="x==b ";break;
				case 4:a="b<x<c";break;
				case 5:a="x==c ";break;
				case 6:a="c<x<d";break;
				case 7:a="x==d ";break;
				case 8:a="d<x  ";break;
				}
				printf("  %5.2f%% %s", h2[kr].prob, a);
			}
			printf("\n");
#endif
			++ncases;
		}
	/*	if(sum)
		{
			printf("[%3d] %s", kp, arr);
			for(int kr=0;kr<9;++kr)
			{
				char c1, c2;
				if(kr==kmax)
					c1='[', c2=']';
				else
					c1=c2=' ';
				printf(" %c%9lf%%%c", c1, hk[kr]*100./sum, c2);
			}
			printf("\t%9lf%%\n", sum*100./totalsum);
			++ncases;
		}*/

		//printf("[%3d] %s", kp, arr);
		//if(!sum)
		//	printf(" 0 0 0 0 0 0 0 0 0\n");
		//else
		//{
		//	for(int kr=0;kr<9;++kr)
		//		printf(" %9lf", (double)hk[kr]*100./sum);
		//	printf("\n");
		//}
	}
	printf("\n%d cases with P > 0.5%%\n", ncases);
	printf("const neighbors: %lld / %lld = %lf%%\n", e10_constcount, e10_totalpixelcount, 100.*e10_constcount/e10_totalpixelcount);
	free(sum);
	free(priority);
}
#endif

void estimate_csize_from_transforms(const unsigned char *buf, unsigned char *b2, int iw, int ih, double *csize)
{
	int len=iw*ih<<2;
	memcpy(b2, buf, len);
	apply_transforms_fwd(b2, iw, ih);
	double entropy[3]={0};
	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(b2+kc, len, 4, hist+((size_t)kc<<8));
		for(int k=0;k<256;++k)
		{
			int freq=hist[kc<<8|k];
			double p=(double)freq/(len>>2);
			if(freq)
				entropy[kc]-=p*log2(p);
		}
	}
	csize[0]=entropy[0]*(len>>2)/8;
	csize[1]=entropy[1]*(len>>2)/8;
	csize[2]=entropy[2]*(len>>2)/8;
}


//	#define BATCHTEST_NO_B2

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
#if 0	
	long long *hist=e10_start();//
	if(!hist)
	{
		exit(0);
		return;
	}
#endif

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

		int iw=0, ih=0, nch0=3, stride=4;
		long long cycles=__rdtsc();
		unsigned char *buf=image_load(fn[0]->data, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}

		ptrdiff_t res=(ptrdiff_t)iw*ih, len=res*stride, usize=res*nch0;
		double ratio=(double)usize/formatsize;
#ifndef BATCHTEST_NO_B2
		printf("%3lld/%3lld  \"%s\"\tCR %lf (%lf BPP) Dec %lf CPB", k+1, filenames->count, fn[0]->data, ratio, 8/ratio, (double)cycles/usize);
#else
		printf("%3lld/%3lld  %.2lf%%\r", k+1, filenames->count, (k+1)*100./filenames->count);
#endif
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
#ifndef BATCHTEST_NO_B2
		unsigned char *b2=(unsigned char*)malloc(len);
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memset(b2, 0, len);
#endif
		
	//test 23: test 16 optimizer
#if 0
		{
			ArrayHandle cdata=0;
			printf("T23\n");
			cycles=__rdtsc();
			test23_encode(buf, iw, ih, &cdata, 1, 0);
			cycles=__rdtsc()-cycles;
			printf("Enc %11lf CPB  CR %9lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");

			cycles=__rdtsc();
			test23_decode(cdata->data, cdata->count, iw, ih, b2, 1);
			cycles=__rdtsc()-cycles;
			printf("Dec %11lf CPB\n", (double)cycles/usize);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T23", 0);
			memset(b2, 0, len);
			printf("\n");
		}
#endif

		//E10 codec
#if 0
		{
			ArrayHandle cdata=0;
			cycles=__rdtsc();
			e10_encode_ch(buf, iw, ih, 0, &cdata, 0);
			e10_encode_ch(buf, iw, ih, 1, &cdata, 0);
			e10_encode_ch(buf, iw, ih, 2, &cdata, 0);
			cycles=__rdtsc()-cycles;
			printf(" Enc %lf CPB  CR %lf  csize %lld", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");

			array_free(&cdata);
			printf("\n");
		}
#endif

		//experiment 10: predictor sorts neighbors
#if 0
		addbuf(buf, iw, ih, 3, 4, 128);
		colortransform_ycocb_fwd((char*)buf, iw, ih);
		addbuf(buf, iw, ih, 3, 4, 128);
		//apply_transforms_fwd(buf, iw, ih);

		e10_iter(buf, iw, ih, 1, hist);
#endif

		//test 19
#if 0
		{
			ArrayHandle cdata=0;
			int alpha=0xD3E4, block=30, margin=60;

			printf("T19 a 0x%04X b %3d m %3d ", alpha, block, margin);
			cycles=__rdtsc();
			test19_encode(buf, iw, ih, alpha, block, margin, &cdata, 1);
			cycles=__rdtsc()-cycles;
			printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf("!!! ");

			cycles=__rdtsc();
			test19_decode(cdata->data, cdata->count, iw, ih, alpha, block, margin, b2);
			cycles=__rdtsc()-cycles;
			printf("Dec %11lf CPB ", (double)cycles/usize);

			//printf("\n");
			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T19", 0);
			memset(b2, 0, len);
			//printf("\n");
		}
#endif
		
		//test16
#if 1
		{
			ArrayHandle cdata=0;
			int alpha=0xD3E7,
				blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
				blockh[]={ 1,  1,  1},
				margin[]={26, 37, 26};

			printf("T16\n");
			cycles=__rdtsc();
			test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, 0);
			cycles=__rdtsc()-cycles;
			printf("Enc %lf CPB  CR %lf  csize %lld", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");

			cycles=__rdtsc();
			test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
			cycles=__rdtsc()-cycles;
			printf("Dec %lf CPB\n", (double)cycles/usize);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, 3, 4, "T16", 0);
			memset(b2, 0, len);

			printf("\n");
		}
#endif
		//test16 estimate
#if 0
		apply_transforms_fwd(buf, iw, ih);
		double csize=test16_estimate_csize(buf, iw, ih, 32, 0.6f, 0);
		sum_testsize+=(long long)ceil(csize);
		printf("\tCR2 %f", usize/csize);
		if(csize<formatsize)
			printf(" !!!");
		printf("\n");
#endif

		//transforms + ANS
#if 0
		apply_transforms_fwd(buf, iw, ih);

		ArrayHandle cdata=0;
		printf("rANS\n");
		cycles=__rdtsc();
		rans4_encode(buf, (ptrdiff_t)iw*ih, 3, 4, &cdata, 0);
		cycles=__rdtsc()-cycles;
		printf("Enc CPB %lf ratio %lf\n", (double)cycles/usize, (double)usize/cdata->count);
		
		sum_testsize+=cdata->count;

		cycles=__rdtsc();
		rans4_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, 3, 4, b2, 0);
		cycles=__rdtsc()-cycles;
		printf("Dec CPB %lf\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, 3, 4, "rANS", 0);
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

		//tests 8 / 9 / 11
#if 0
		apply_transforms_fwd(buf, iw, ih);
		double csize=test12_estimate_csize(buf, iw, ih, 16, 0);
	//	double csize=test11_estimate_cr(buf, iw, ih, 0);
	//	double csize=test9_estimate_cr(buf, iw, ih, 0);
	//	double csize=test8_estimate_cr(buf, iw, ih, 0);
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
#ifndef BATCHTEST_NO_B2
		free(b2);
#endif
	}
#if 0
	e10_print(hist);
	free(hist);
#else
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
#endif

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
	//test8();
	//test9();

	printf("EntropyBattle\n");
#if 1
	long long cycles;
	int iw=0, ih=0, nch0=3,
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
		buf=image_load(fn, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			LOG_ERROR("Couldn't open \"%s\"", fn);
			return 0;
		}
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		printf("Format Dec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", (double)cycles/(resolution*nch0), iw, ih, nch0, formatsize, (double)resolution*nch0/formatsize);
	}
	else if(argc==3)
	{
		const char *fn1=argv[1], *fn2=argv[2];
		int w2, h2;
		buf=image_load(fn1, &iw, &ih);
		b2 =image_load(fn2, &w2, &h2);
		if(!buf)
		{
			printf("Couldn't open %s\n", fn1);
			return 1;
		}
		if(!b2)
		{
			printf("Couldn't open %s\n", fn2);
			return 1;
		}
		if(iw!=w2||ih!=h2)
		{
			printf("Expected two images of SAME RESOLUTION. %dx%d != %dx%d\n", iw, ih, w2, h2);
			return 1;
		}
		ptrdiff_t formatsize=get_filesize(fn2);
		int res=iw*ih;
		long long sum[3]={0};
		for(int k=0;k<res;++k)
		{
			int dr=buf[k<<2  ]-b2[k<<2  ],
				dg=buf[k<<2|1]-b2[k<<2|1],
				db=buf[k<<2|2]-b2[k<<2|2];
			sum[0]+=dr*dr;
			sum[1]+=dg*dg;
			sum[2]+=db*db;
		}
		double rmse[]=
		{
			sqrt((double)sum[0]/res),
			sqrt((double)sum[1]/res),
			sqrt((double)sum[2]/res),
			sqrt((double)(sum[0]+sum[1]+sum[2])/(res*3)),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		double CR=res*3./formatsize;
		printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)formatsize, CR, 8/CR);
		printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
		printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
		printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
		return 0;
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
#if 0
	extern const unsigned char *debug_ptr;
	//debug_ptr=buf;
	printf("test10\n");
	cycles=__rdtsc();
	test10_encode(buf, iw, ih, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  CR %lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
	
	cycles=__rdtsc();
	test10_decode(cdata->data, cdata->count, iw, ih, b2);
	cycles=__rdtsc()-cycles;
	printf("Dec %lf CPB\n", (double)cycles/usize);

	//lodepng_encode_file("kodim21-test10.PNG", b2, iw, ih, LCT_RGBA, 8);

	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "test10", 0);
	memset(b2, 0, len);

	printf("\n");
#endif

	//test11_estimate_cr(buf, iw, ih, 1);

	//test16
#if 0
	{
		//debug_ptr=buf;//
		int besta=0, bestb=0, bestm=0, bestc=0;
		int it=0;
		for(int m=30;m<=30;++m)
		{
			for(int b=15;b<=15;++b)
			{
				for(int a=0xD3E7;a<=0xD3E7;++a, ++it)
				{
					int alpha=a,//(60<<16)/100;	0xA51F	0xD3DA	0xD3D6	0xD3EB
						bsize=b,//32 24 26 20
						margin=m;
					printf("T16 a 0x%04X b %3d m %3d ", alpha, bsize, margin);
					cycles=__rdtsc();
					test16_encode(buf, iw, ih, alpha, bsize, margin, &cdata, 1);
					cycles=__rdtsc()-cycles;
					printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

					if(!it||bestc>(int)cdata->count)
						besta=alpha, bestb=bsize, bestm=m, bestc=(int)cdata->count;

					cycles=__rdtsc();
					test16_decode(cdata->data, cdata->count, iw, ih, alpha, bsize, margin, b2);
					cycles=__rdtsc()-cycles;
					printf("Dec %11lf CPB ", (double)cycles/usize);

					array_free(&cdata);
					compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
					memset(b2, 0, len);
				}
				//printf("\n");
			}
		}
		if(it>1)
			printf("T16 best ABMS 0x%04X %2d %2d %d CR %lf\n", besta, bestb, bestm, bestc, (double)usize/bestc);
		else
			printf("\n");
	}
#endif

	//test17_saveconf(buf, iw, ih, 24);

	//test 17
#if 0
	{
		debug_ptr=buf;//
		int bsize=24,//32 24
			alpha=0xA51F;//(60<<16)/100;
		printf("T17 alpha 0x%04X block %d\n", alpha, bsize);
		cycles=__rdtsc();
		test17_encode(buf, iw, ih, bsize, alpha, &cdata, 1);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  CR %lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

		cycles=__rdtsc();
		test17_decode(cdata->data, cdata->count, iw, ih, bsize, alpha, b2);
		cycles=__rdtsc()-cycles;
		printf("Dec %lf CPB\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T17", 0);
		memset(b2, 0, len);

		printf("\n");
	}
#endif

#if 0
	printf("T18 static ANS on DWT subbands & channels\n");
	cycles=__rdtsc();
	test18_encode(buf, iw, ih, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  CR %lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

	array_free(&cdata);
	printf("\n");
#endif

	//test 19
#if 0
	{
		//int a=0xD3E7, b=15, m=30;
		
		int besta=0, bestb=0, bestm=0, bestc=0;
		int it=0;
		for(int m=5;m<=5;++m)
		{
			for(int b=10;b<=10;++b)
			{
				for(int a=0xD3E4;a<=0xD3E4;++a, ++it)
				{
					//debug_ptr=buf;//
					printf("T19 a 0x%04X b %3d m %3d ", a, b, m);
					cycles=__rdtsc();
					test19_encode(buf, iw, ih, a, b, m, &cdata, 0);
					cycles=__rdtsc()-cycles;
					printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

					if(!it||bestc>(int)cdata->count)
						besta=a, bestb=b, bestm=m, bestc=(int)cdata->count;
#if 1
					cycles=__rdtsc();
					test19_decode(cdata->data, cdata->count, iw, ih, a, b, m, b2);
					cycles=__rdtsc()-cycles;
					printf("Dec %11lf CPB ", (double)cycles/usize);

					//printf("\n");
					compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T19", 0);
					memset(b2, 0, len);
#endif
					printf("\n");
					array_free(&cdata);
				}
			}
		}
		if(it>1)
			printf("T19 best ABMS 0x%04X %2d %2d %d CR %lf\n", besta, bestb, bestm, bestc, (double)usize/bestc);
		else
			printf("\n");
	}
#endif
	
	//test20
#if 0
	{
		//debug_ptr=buf;//
		int loud=1;//
		int bestB=0, besta=0, bestb=0, bestm=0, bestc=0;
		int it=0;
		for(int blockcount=8;blockcount<=8;++blockcount)
		{
			for(int margin=30;margin<=30;++margin)
			{
				for(int blocksize=1;blocksize<=1;++blocksize)
				{
					for(int alpha=0;alpha<=0;++alpha, ++it)//0xD3E7
					{
						printf("T20 B %3d a 0x%04X b %3d m %3d ", blockcount, alpha, blocksize, margin);
						cycles=__rdtsc();
						test20_encode(buf, iw, ih, blocksize, margin, alpha, blockcount, &cdata, loud);
						cycles=__rdtsc()-cycles;
						printf("Enc %11lf CPB  CR %9lf  csize %7lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

						if(!it||bestc>(int)cdata->count)
							bestB=blockcount, besta=alpha, bestb=blocksize, bestm=margin, bestc=(int)cdata->count;

						cycles=__rdtsc();
						test20_decode(cdata->data, cdata->count, iw, ih, blocksize, margin, alpha, blockcount, b2);
						cycles=__rdtsc()-cycles;
						printf("Dec %11lf CPB ", (double)cycles/usize);

						array_free(&cdata);
						compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T20", 0);
						memset(b2, 0, len);
					}
					//printf("\n");
				}
			}
		}
		if(it>1)
			printf("T20 best Babms %3d 0x%04X %2d %2d %d CR %lf\n", bestB, besta, bestb, bestm, bestc, (double)usize/bestc);
		else
			printf("\n");
	}
#endif

	//E10 estimate: predictor that sorts neighbors
#if 0
	{
		printf("E10 estimate\n");
		memcpy(b2, buf, len);
		addbuf(b2, iw, ih, 3, 4, 128);
		colortransform_ycocb_fwd((char*)b2, iw, ih);
		//pred_grad_fwd((char*)b2, iw, ih, 3, 4);
		addbuf(b2, iw, ih, 3, 4, 128);

		double csize, totalcsize=0;
		long long *hist=e10_start();
		for(int kc=0;kc<3;++kc)
		{
			e10_iter(b2, iw, ih, kc, hist);
			csize=e10_estimate_csize(b2, iw, ih, kc, hist, 1);
			//printf("ch %d  U %d  C %lf  CR %lf\n", kc, iw*ih, csize, iw*ih/csize);
			totalcsize+=csize;
			//e10_print(hist);
			e10_reset(hist);
			//printf("\n\n");
		}
		free(hist);
		printf("total usize %d  csize %lf  CR %lf\n\n", (int)usize, totalcsize, usize/totalcsize);
		memset(b2, 0, len);
	}
#endif

	//E10 codec
#if 0
	{
		int alpha=0xD3E7, blocksize=15, margin=30;
		printf("E10 sortpred\n");
		cycles=__rdtsc();
		e10dash_encode_ch(buf, iw, ih, 0, alpha, blocksize, margin, &cdata, 1);
		e10dash_encode_ch(buf, iw, ih, 1, alpha, blocksize, margin, &cdata, 1);
		e10dash_encode_ch(buf, iw, ih, 2, alpha, blocksize, margin, &cdata, 1);

		//e10_encode_ch(buf, iw, ih, 0, &cdata, 1);
		//e10_encode_ch(buf, iw, ih, 1, &cdata, 1);
		//e10_encode_ch(buf, iw, ih, 2, &cdata, 1);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  CR %lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		
		cycles=__rdtsc();
		size_t offset=0;
		offset=e10dash_decode_ch(cdata->data, offset, cdata->count, iw, ih, 0, alpha, blocksize, margin, b2);
		offset=e10dash_decode_ch(cdata->data, offset, cdata->count, iw, ih, 1, alpha, blocksize, margin, b2);
		offset=e10dash_decode_ch(cdata->data, offset, cdata->count, iw, ih, 2, alpha, blocksize, margin, b2);
		//offset=e10_decode_ch(cdata->data, offset, cdata->count, iw, ih, 0, b2);
		//offset=e10_decode_ch(cdata->data, offset, cdata->count, iw, ih, 1, b2);
		//offset=e10_decode_ch(cdata->data, offset, cdata->count, iw, ih, 2, b2);
#if 1
		addbuf(b2, iw, ih, 3, 4, 128);
		pred_grad_inv((char*)b2, iw, ih, 3, 4);
		colortransform_ycocb_inv((char*)b2, iw, ih);
		addbuf(b2, iw, ih, 3, 4, 128);
#endif
		cycles=__rdtsc()-cycles;
		printf("Dec %11lf CPB ", (double)cycles/usize);
		
		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "E10", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif


	//test16 codec with jxl predictor optimizer
#if 1
	{
		int alpha=0xD3E7,
			blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
			blockh[]={ 1,  1,  1},
			margin[]={26, 37, 26};

#if 0
		int res=iw*ih;
		double step=0.001, CR0=0, CR, csize[3]={0};
		estimate_csize_from_transforms(buf, b2, iw, ih, csize);
		CR=res*3/(csize[0]+csize[1]+csize[2]);
		printf("%4d TRGB %lf [%lf %lf %lf]\n", 0, CR, res/csize[0], res/csize[1], res/csize[2]);

		for(int k=0;k<256;++k)
		{
			int idx=k%33;
			if(!(k+1)%33)
				step*=0.9;
			do
			{
				CR0=CR;
				jxlpred_params[idx]+=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);

			do
			{
				CR0=CR;
				jxlpred_params[idx]-=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);
		}
#endif

		//for(int k=0;k<3;++k)
		//	printf("%g\t%g\t%g\t%g\n", jxlpred_params[k<<2], jxlpred_params[k<<2|1], jxlpred_params[k<<2|2], jxlpred_params[k<<2|3]);
		//for(int k=12;k<33;k+=7)
		//	printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", jxlpred_params[k], jxlpred_params[k+1], jxlpred_params[k+2], jxlpred_params[k+3], jxlpred_params[k+4], jxlpred_params[k+5], jxlpred_params[k+6]);
		//for(int k=0;k<33;++k)
		//	printf("%3d  %lf\n", k, jxlpred_params[k]);
		//printf("\n");
		
		printf("T16\n");
		int bestcsizes[3]={0}, bestw[3]={0}, besth[3]={0}, bestm[3]={0};
		int it=0;
		//for(int bw=19;bw<25;++bw)//
		{
			//for(int m=31;m<48;++m)
			{
				int csizes[3];
				//blockw[1]=bw, blockh[1]=1, margin[1]=m;
				//blockw[1]=1+k, blockh[1]=1;
				//blockw[1]=4+k%10, blockh[1]=1+k/10;
			
				//blockw[0]=blockw[2]=blockw[1];
				//blockh[0]=blockh[2]=blockh[1];

				cycles=__rdtsc();
				test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, csizes);
				cycles=__rdtsc()-cycles;
				printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

				for(int kc=0;kc<3;++kc)
				{
					if(!it||bestcsizes[kc]>csizes[kc])
						bestcsizes[kc]=csizes[kc], bestw[kc]=blockw[kc], besth[kc]=blockh[kc], bestm[kc]=margin[kc];
				}

				cycles=__rdtsc();
				test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
				cycles=__rdtsc()-cycles;
				printf("Dec %11lf CPB ", (double)cycles/usize);

				array_free(&cdata);
				compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
				memset(b2, 0, len);
				printf("\n");
				++it;
			}
		}
		int res=iw*ih;
		printf("R %7d %lf %dx%d  M %d\n", bestcsizes[0], (double)res/bestcsizes[0], bestw[0], besth[0], bestm[0]);
		printf("G %7d %lf %dx%d  M %d\n", bestcsizes[1], (double)res/bestcsizes[1], bestw[1], besth[1], bestm[1]);
		printf("B %7d %lf %dx%d  M %d\n", bestcsizes[2], (double)res/bestcsizes[2], bestw[2], besth[2], bestm[2]);
		printf("\n");
	}
#endif


	//test 21: reorder blocks
#if 0
	{
		int alpha=0xD3E7,
			blocksize=8;
		//	margin=32;

		//apply_transforms_fwd(buf, iw, ih);

		//memcpy(b2, buf, len);
		//apply_transforms_fwd(b2, iw, ih);
		//apply_transforms_inv(b2, iw, ih);
#if 1
		printf("T21 alpha 0x%04X blocksize %3d\n", alpha, blocksize);
		cycles=__rdtsc();
		test21_encode(buf, iw, ih, alpha, blocksize, &cdata, 1);
		cycles=__rdtsc()-cycles;
		printf("Enc %11lf CPB  CR %9lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

		cycles=__rdtsc();
		test21_decode(cdata->data, cdata->count, iw, ih, alpha, blocksize, b2);
		cycles=__rdtsc()-cycles;
		printf("Dec %11lf CPB\n", (double)cycles/usize);
		array_free(&cdata);

		//lodepng_encode_file("t21.PNG", b2, iw, ih, LCT_RGBA, 8);//
#endif

		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T21", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//test 22: interlaced parallelogram blocks
#if 0
	{
		int alpha=0xD3E7,
			blockw[]={232, 232, 232},
			blockh[]={232, 232, 232};

		//apply_transforms_fwd(buf, iw, ih);

		//memcpy(b2, buf, len);
		//apply_transforms_fwd(b2, iw, ih);
		//apply_transforms_inv(b2, iw, ih);
#if 1
		//float x=0.1f;//1/7
		//printf("%d\n", x==0.1);

		printf("T22\n");
		int bestcsizes[3]={0}, bestw[3]={0}, besth[3]={0};
		for(int k=0;k<1;++k)//
		{
			int csizes[3];
			//blockw[1]=16+k%10*24, blockh[1]=16+k/10*24;
			//
			//blockw[0]=blockw[2]=blockw[1];
			//blockh[0]=blockh[2]=blockh[1];

			//printf("T22 alpha 0x%04X block %dx%d\n", alpha, blockw, blockh);
			cycles=__rdtsc();
			test22_encode(buf, iw, ih, alpha, blockw, blockh, &cdata, 1, csizes);
			cycles=__rdtsc()-cycles;
			printf("Enc %11lf CPB  CR %9lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

			for(int kc=0;kc<3;++kc)
			{
				if(!k||bestcsizes[kc]>csizes[kc])
					bestcsizes[kc]=csizes[kc], bestw[kc]=blockw[kc], besth[kc]=blockh[kc];
			}

			cycles=__rdtsc();
			test22_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, b2);
			cycles=__rdtsc()-cycles;
			printf("Dec %11lf CPB\n", (double)cycles/usize);
			array_free(&cdata);

			compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T22", 0);
			memset(b2, 0, len);
		}
		int res=iw*ih;
		printf("\n");
		printf("R %7d %lf %dx%d\n", bestcsizes[0], (double)res/bestcsizes[0], bestw[0], besth[0]);
		printf("G %7d %lf %dx%d\n", bestcsizes[1], (double)res/bestcsizes[1], bestw[1], besth[1]);
		printf("B %7d %lf %dx%d\n", bestcsizes[2], (double)res/bestcsizes[2], bestw[2], besth[2]);

		//lodepng_encode_file("t21.PNG", b2, iw, ih, LCT_RGBA, 8);//
#endif
		printf("\n");
	}
#endif

	//test 23: test 16 optimizer
#if 0
	{
		printf("T23\n");
		cycles=__rdtsc();
		test23_encode(buf, iw, ih, &cdata, 1, 0);
		cycles=__rdtsc()-cycles;
		printf("Enc %11lf CPB  CR %9lf  csize %lld\n", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

		cycles=__rdtsc();
		test23_decode(cdata->data, cdata->count, iw, ih, b2, 1);
		cycles=__rdtsc()-cycles;
		printf("Dec %11lf CPB\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T23", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//YCoCb
#if 0
	{
		int color=0xFF000000;
		memfill(b2, &color, len, 4);
		for(int ky=0;ky<256;++ky)
		{
			for(int kx=0;kx<256;++kx)
			{
				int idx=iw*ky+kx;
				
				//YCoCg
				b2[idx<<2  ]=kx-128;//Cr
				b2[idx<<2|1]=ky-128;//Cb
				b2[idx<<2|2]=0;

				//YCoCgT
				//b2[idx<<2  ]=kx-128;//Cr
				//b2[idx<<2|1]=0;//Y
				//b2[idx<<2|2]=ky-128;//Cb
			}
		}
		colortransform_ycocg_inv((char*)b2, iw, ih);
		addbuf(b2, iw, ih, 3, 4, 128);
		image_save_png_rgba8("YCoCg.PNG", b2, iw, ih);//
	}
#endif

	//truncation experiment
#if 0
	{
		int *b3=(int*)malloc(len*sizeof(int));
		int *b4=(int*)malloc(len*sizeof(int));
		int *temp=(int*)malloc(iw*10LL*sizeof(int));
		cvt_i8_i32(buf, iw, ih, b3);
		addbuf_i32(b3, iw, ih, 3, 4, -128);
		colortransform_ycocb_fwd_i32(b3, iw, ih);
		for(int kc=0;kc<3;++kc)
			pred_jxl_i32(b3+kc, iw, ih, kc, jxlparams_i32, 1, b4, temp);
		double csize=estimate_csize_i32(b4, iw, ih, 10, 9, 10);
		double csize0=estimate_csize_i32_trunc(b4, iw, ih, 10, 9, 10);
		printf("Overflow test\n");
		printf("est. complete  csize %lf  usize %lld  CR %lf\n", csize, usize, usize/csize);
		printf("est. truncated csize %lf  usize %lld  CR %lf\n", csize0, usize, usize/csize0);
		printf("\n");
		free(b3);
		free(b4);
		free(temp);
	}
#endif

	//experiment 24: adaptive test16 params
#if 1
	{
		printf("E24\n");
		int res=iw*ih;
		double bestsize=0;
		int bestw=0, besti=0, beste=0;
		double csizes[4];
		int it=0;
		for(int kw=15;kw<=15;++kw)
		for(int ki=56;ki<=56;++ki)
		for(int ke=8;ke<=8;++ke)
		{
			unsigned char minwidth[]={kw, 40, 15};//40	kw;//126;
			unsigned char maxinc[]={ki, 84, 56};//84	ki;//121;
			unsigned char encounter_threshold[]={ke, 8, 8};//8	ke;//31;//0xBF		255-100+k
			printf("W %3d I %3d E %3d %9.6lf%%\t", (int)minwidth[2], (int)maxinc[2], encounter_threshold[2], encounter_threshold[2]*100./0xFF);
			csizes[0]=e24_estimate(buf, iw, ih, minwidth, maxinc, encounter_threshold, csizes+1);
			//printf("\n");
			if(!it||bestsize>csizes[0])
				bestsize=csizes[0], bestw=minwidth[2], besti=maxinc[2], beste=encounter_threshold[2];
			++it;
		}
		printf("best W %d  I %d  E %d  CR %lf\n", bestw, besti, beste, usize/bestsize);
		//csizes[0]=e24_estimate(buf, iw, ih, 8, 64, 0xBF, csizes+1);

		//printf("T %14lf  CR %lf\n", csizes[0], usize/csizes[0]);
		//printf("R %14lf  CR %lf\n", csizes[1], res/csizes[1]);
		//printf("G %14lf  CR %lf\n", csizes[2], res/csizes[2]);
		//printf("B %14lf  CR %lf\n", csizes[3], res/csizes[3]);
		printf("\n");
	}
#endif

	//predict image
	apply_transforms_fwd(buf, iw, ih);
	//lodepng_encode_file("kodim21-YCoCgT-unplane.PNG", buf, iw, ih, LCT_RGBA, 8);//
	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//
#if 0
	printf("Predict image...\n");
	{
		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		char *temp=(char*)malloc(MAXVAR(iw, ih));
		
		addbuf(buf, iw, ih, nch0, nch, 128);//unsigned char -> signed char
		
		//colortransform_ycocg_fwd((char*)buf, iw, ih);
		//colortransform_xgz_fwd((char*)buf, iw, ih);
		//colortransform_xyz_fwd((char*)buf, iw, ih);

		//char *b3=(char*)malloc(iw), *b4=(char*)malloc(iw);
		//if(!b3||!b4)
		//	return 0;
		//memcpy(b3, buf, iw);
		//dwt1d_squeeze_fwd(b3, iw, 1, b4);
		//dwt1d_squeeze_inv(b3, iw, 1, b4);
		//compare_bufs_uint8((unsigned char*)b3, buf, iw, 1, 1, 1, "squeeze row", 0);
		//free(b3);
		//free(b4);

#if 1
		memcpy(b2, buf, len);

		colortransform_ycocb_fwd((char*)b2, iw, ih);
		float jxlparams[33]=
		{
			 0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
			 0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
			 0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
		};
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 1);
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 0);

		//colortransform_xyz_fwd(b2, iw, ih);
		//colortransform_xyz_inv(b2, iw, ih);
		//for(int kc=0;kc<3;++kc)
		//{
		//	dwt2d_squeeze_fwd((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//	dwt2d_squeeze_inv((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//}
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "transform", 0);
		printf("\n");
#endif
		
		//for(int kc=0;kc<3;++kc)
		//	//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	dwt2d_squeeze_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, 2, 4, (char*)temp);
		//	//dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

		addbuf(buf, iw, ih, nch0, nch, 128);

		//save_DWT_int8("kodim21-squeeze-stage", buf, (DWTSize*)sizes->data, 2, 4);//
		//lodepng_encode_file("kodim21-cubic.PNG", buf, iw, ih, LCT_RGBA, 8);//

		array_free(&sizes);
		free(temp);
	}
	//squeeze_8bit_lossy(buf, iw, ih, nch0, nch);
//	image_pred(buf, iw, ih, nch0, nch);

	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//
#endif

	//tests 11~13, 15
#if 0
	printf("Test 8 estimate - shifted conditional histogram:\n");
	test8_estimate_cr(buf, iw, ih, 1);

	printf("Test 9 estimate - global mean/conf:\n");
	test9_estimate_cr(buf, iw, ih, 1);

	printf("Test 11 estimate - global joint histogram:\n");
	test11_estimate_cr(buf, iw, ih, 1);

	printf("Test 12 estimate - per block joint histogram:\n");
	for(int bsize=2, bend=MAXVAR(iw, ih);bsize<bend;bsize<<=1)
		test12_estimate_csize(buf, iw, ih, bsize, 1);
	printf("\n");
	
	printf("Test 13 estimate - per block joint mean/var:\n");
	//for(int bsize=1<<floor_log2(MAXVAR(iw, ih));bsize>=2;bsize>>=1)
	for(int bsize=2, bend=MAXVAR(iw, ih);bsize<bend;bsize<<=1)
		test13_estimate_csize(buf, iw, ih, bsize, 1);
	printf("\n");

	printf("Test 15 estimate - per channel per block histogram:\n");
	for(int bsize=2, bend=MAXVAR(iw, ih);bsize<bend;bsize<<=1)
		test15_estimate_csize(buf, iw, ih, bsize, 1);
	printf("\n");
#endif

	//test 14
#if 0
	printf("test14\n");
	cycles=__rdtsc();
	test14_encode(buf, iw, ih, 1, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  CR %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	array_free(&cdata);
	printf("\n");
#endif

	//test 16 estimate
#if 0
	{
		float bestcsize=0;
		int bestalpha=0, bestblock=0;
		printf("Test 16 estimate - per channel previous block histogram:\n");
		for(int bsize=2, bend=MAXVAR(iw, ih);bsize<bend;bsize<<=1)
		{
			printf("block size %d (%d blocks):\n", bsize, iw*ih*3/(bsize*bsize));
			for(int alpha=0;alpha<100;++alpha)
			{
				float csize=(float)test16_estimate_csize(buf, iw, ih, bsize, alpha/100.f, 1);
				if(bsize==2&&!alpha||bestcsize>csize)
					bestcsize=csize, bestblock=bsize, bestalpha=alpha;
			}
			printf("\n");
		}
		printf("Best: block %d alpha %d%% CR %f csize %f\n\n", bestblock, bestalpha, (float)usize/bestcsize, bestcsize);
	}
	//for(int alpha=0;alpha<100;++alpha)
	//{
	//	printf("static + (adaptive-static) * %f%%\n", alpha/100.f);
	//	for(int bsize=2, bend=MAXVAR(iw, ih);bsize<bend;bsize<<=1)
	//		test16_estimate_csize(buf, iw, ih, bsize, alpha/100.f, 1);
	//	printf("\n");
	//}
#endif

	for(int kc=0;kc<nch0;++kc)
		calc_histogram(buf+kc, len, nch, hist+((size_t)kc<<8));
	//print_histogram(hist, 1);

	double entropy[6]={0};
	for(int kc=0;kc<nch0;++kc)
	{
		int freq;
		double p;
		for(int k=0;k<256;++k)
		{
			freq=hist[kc<<8|k];
			if(freq)
			{
				p=(double)freq/(len>>2);
				p*=0x10000-255;
				++p;
				p/=0x10000;
				entropy[kc]+=-p*log2(p);
			}

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
	rans4_encode(buf, (ptrdiff_t)iw*ih, 3, 4, &cdata, 0);
	//rans4_encode(buf, len, 1, 0, &cdata, loud, pred);
	cycles=__rdtsc()-cycles;
	//double cr=(double)usize/cdata->count, cr0=8/entropy[4];
	//printf("Enc CPB %lf ratio %lf%s\n", (double)cycles/usize, cr, cr>cr0?" IMPOSSIBLE":"");
	printf("Enc %lf CPB  csize %d  CR %lf\n", (double)cycles/usize, (int)cdata->count, (double)usize/cdata->count);
	
	cycles=__rdtsc();
	rans4_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, 3, 4, b2, 0);
	//rans4_decode(cdata->data, cdata->count, len, 1, 0, b2, loud, pred);
	cycles=__rdtsc()-cycles;
	printf("Dec %lf CPB\n", (double)cycles/usize);

	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, 3, 4, "rANS", 0);
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

	//ANS128
#if 0
	{
		int depths[]={nch0<<3};
	//	int depths[]={8, 8, 8};
		printf("ANS128\n");
		cycles=__rdtsc();
		size_t puresize=ans128_encode(buf, (ptrdiff_t)iw*ih, depths, COUNTOF(depths), nch, &cdata);
		cycles=__rdtsc()-cycles;
		printf("Enc %lf CPB  CR %lf  CR2 %lf\n", (double)cycles/usize, (double)usize/cdata->count, (double)usize/puresize);
	
		cycles=__rdtsc();
		ans128_decode(cdata->data, cdata->count, (ptrdiff_t)iw*ih, depths, COUNTOF(depths), nch, b2);
		cycles=__rdtsc()-cycles;
		printf("Dec %lf CPB\n", (double)cycles/usize);

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "ANS128", 0);
		memset(b2, 0, len);
	}
#endif

	//ANS16
#if 0
	printf("ANS16\n");
	cycles=__rdtsc();
	ans16_encode(buf, (ptrdiff_t)iw*ih, &cdata);
	cycles=__rdtsc()-cycles;
	printf("Enc %lf CPB  CR %lf\n", (double)cycles/usize, (double)usize/cdata->count);
	array_free(&cdata);
	printf("\n");
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