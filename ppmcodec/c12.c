#include"codec.h"
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
//#include"blist.h"
static const char file[]=__FILE__;


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
//	#define AC_VALIDATE
#endif

//	#define BYPASSALL

#define GRLIMIT 16
#define PROBBITS_STORE 24
#define PROBBITS_USE 9
#define PROBBITS_ADJUST 24
#define PROBBITS_USE2 12
//#define ABAC_PROBBITS 16
//#define ANALYSIS_XSTRIDE 4
//#define ANALYSIS_YSTRIDE 4


#ifdef AC_VALIDATE
#define AC_IMPLEMENTATION
#include"entropy.h"
#endif
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
//	OCH_R2=6,
//	OCH_G2,
//	OCH_B2,

	OCH_COUNT=6,
} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},// 0
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},// 1
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},// 2
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},// 3
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},// 4
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},// 5
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},// 6
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},// 7
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},// 8
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},// 9
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},//10
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},//11
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},//12
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},//13
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},//14
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},//15
};
static unsigned stats1[3][256][14][GRLIMIT+8];
//static unsigned probtable[6][1<<PROBBITS_USE];
#if 0
AWM_INLINE int squash(int x)//sigmoid(x) = 1/(1+exp(-x))		logit sum -> prob
{
#ifdef DISABLE_LOGMIX
	x>>=11;
	x+=1<<ABAC_PROBBITS>>1;
	CLAMP2_32(x, x, 1, (1<<ABAC_PROBBITS)-1);
#else
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(ABAC_PROBBITS-5))-1);
	x=(x>>(ABAC_PROBBITS-5))+16;
	if(x>31)
		return (1<<ABAC_PROBBITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(ABAC_PROBBITS-5))-w)+t[x+1]*w+64)>>(12-5);
#endif
	return x;
}
AWM_INLINE int stretch(int x)//ln(x/(1-x))		probs -> logits
{
#ifndef DISABLE_LOGMIX
	static short t[4096];
	static int initialized=0;
	if(!initialized)
	{
		initialized=1;
		
		int pi=0;
		for(int k=-2047;k<=2047;++k)//invert squash()
		{
			int i=squash(k<<(ABAC_PROBBITS-12))>>(ABAC_PROBBITS-12);
			for(int j=pi;j<=i;++j)
				t[j]=k<<(ABAC_PROBBITS-12);
			pi=i+1;
		}
		t[4095]=(1<<ABAC_PROBBITS>>1)-1;
	}
	x=t[x>>(ABAC_PROBBITS-12)];
#endif
	return x;
}
#endif
typedef struct _C12Info
{
	unsigned long long low, range, code;
	unsigned char *ptr;
	int fwd;
	unsigned p0, bit;
//	unsigned *probtable;
} C12Info;
AWM_INLINE void codebit(C12Info *info, int rcon, int sh)
{
	unsigned long long r2;
	unsigned p0=info->p0;
	unsigned p0e=p0>>(PROBBITS_STORE-PROBBITS_USE);
	//unsigned p0e=((p0+(1<<(PROBBITS_STORE-PROBBITS_USE)>>1))>>(PROBBITS_STORE-PROBBITS_USE));
	//unsigned p0f=info->probtable[p0e];
	//unsigned p0g=(p0f+(1<<(PROBBITS_ADJUST-PROBBITS_USE2)>>1))>>(PROBBITS_ADJUST-PROBBITS_USE2);
	//if(!p02)
	//	LOG_ERROR("");
	//if(p0e<1||p0e>(1<<PROBBITS_USE)-1)//
	//	printf("");//
	CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
//	CLAMP2(p0, 0x4000, 0x10000-0x4000);
	if(info->fwd)
	{		
		if(info->range<(1ULL<<PROBBITS_USE))
		{
			*(unsigned*)info->ptr=(unsigned)(info->low>>32);
			info->ptr+=4;
			info->low=info->low<<32;
			info->range=info->range<<32|0xFFFFFFFF;
			if(info->range>~info->low)
				info->range=~info->low;
		}
#ifdef AC_VALIDATE
		unsigned long long lo0=low, r0=range;
#endif
#ifdef BYPASSALL
		r2=info->range>>1;
#else
		r2=info->range*p0e>>PROBBITS_USE;
#endif
		if(info->bit)
		{
			info->low+=r2;
			info->range-=r2;
		}
		else
			info->range=r2-1;
#ifdef AC_VALIDATE
		acval_enc(bit, 0, p0, lo0, lo0+r0, low, low+range, 0, 0);
#endif
	}
	else
	{
		if(info->range<(1ULL<<PROBBITS_USE))
		{
			info->code=info->code<<32|*(unsigned*)info->ptr;
			info->ptr+=4;
			info->low=info->low<<32;
			info->range=info->range<<32|0xFFFFFFFF;
			if(info->range>~info->low)
				info->range=~info->low;
		}
#ifdef AC_VALIDATE
		unsigned long long lo0=low, r0=range;
#endif
#ifdef BYPASSALL
		r2=info->range>>1;
#else
		r2=info->range*p0e>>PROBBITS_USE;
#endif
		unsigned long long mid=info->low+r2;
		info->range-=r2;
		info->bit=info->code>=mid;
		if(info->bit)
			info->low=mid;
		else
			info->range=r2-1;
#ifdef AC_VALIDATE
		acval_dec(bit, 0, p0, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
	}
#ifndef BYPASSALL
	info->p0=p0+((int)((!info->bit<<PROBBITS_STORE)-p0+rcon)>>sh);
	//info->probtable[p0e]=p0f+((int)((!info->bit<<PROBBITS_ADJUST)-p0f+(1<<8>>1))>>8);
#endif
}
#if 0
AWM_INLINE void bypassbit(C12Info *info)
{
	unsigned long long r2;
	if(info->fwd)
	{		
		if(info->range<2)
		{
			*(unsigned*)info->ptr=(unsigned)(info->low>>32);
			info->ptr+=4;
			info->low=info->low<<32;
			info->range=info->range<<32|0xFFFFFFFF;
			if(info->range>~info->low)
				info->range=~info->low;
		}
#ifdef AC_VALIDATE
		unsigned long long lo0=low, r0=range;
#endif
		r2=info->range>>1;
		if(info->bit)
		{
			info->low+=r2;
			info->range-=r2;
		}
		else
			info->range=r2-1;
#ifdef AC_VALIDATE
		acval_enc(bit, 0, p0, lo0, lo0+r0, low, low+range, 0, 0);
#endif
	}
	else
	{
		if(info->range<2)
		{
			info->code=info->code<<32|*(unsigned*)info->ptr;
			info->ptr+=4;
			info->low=info->low<<32;
			info->range=info->range<<32|0xFFFFFFFF;
			if(info->range>~info->low)
				info->range=~info->low;
		}
#ifdef AC_VALIDATE
		unsigned long long lo0=low, r0=range;
#endif
		r2=info->range>>1;

		unsigned long long mid=info->low+r2;
		info->range-=r2;
		info->bit=info->code>=mid;
		if(info->bit)
			info->low=mid;
		else
			info->range=r2-1;
#ifdef AC_VALIDATE
		acval_dec(bit, 0, p0, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
	}
}
#endif
int c12_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double t=time_sec();
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize=get_filesize(srcfn);
	if(srcsize<0)
	{
		LOG_ERROR("Invalid file \"%s\"", srcfn);
		return 1;
	}
	unsigned char *srcbuf=0;
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+16);
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
	}
	unsigned char *srcptr=srcbuf;
	int iw=0, ih=0;
	int fwd=!memcmp(srcptr, "P6\n", 3);
	if(fwd)
	{
		srcptr+=3;
		if(*srcptr=='#')
		{
			while(*srcptr++!='\n');
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr==' ')
			++srcptr;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Invalid file");
			return 1;
		}
		srcptr+=5;
	}
	else
	{
		if(memcmp(srcptr, "12", 2))
		{
			LOG_ERROR("Invalid file");
			return 1;
		}
		srcptr+=2;
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Invalid file");
		return 1;
	}
	
	ptrdiff_t res=(ptrdiff_t)iw*ih, usize=res*3;
	int bestrct=0, flat=0;
	unsigned char *dstbuf=0, *dstptr=0;
	unsigned long long low=0, range=0xFFFFFFFF, code=0;
	if(fwd)
	{
		guide_save(srcptr, iw, ih);
		dstbuf=(unsigned char*)malloc((ptrdiff_t)4*iw*ih);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		dstptr=dstbuf;
		{//analysis
			long long counters[OCH_COUNT]={0}, flatcount[OCH_COUNT]={0};
			const unsigned char *ptr=srcptr, *end=srcptr+usize-(sizeof(__m256i[3])+6-1);

			__m256i mc[OCH_COUNT], fc[OCH_COUNT];
			__m256i shuf=_mm256_set_epi8(
			//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
				-1, -1, -1, 14, -1, 11, -1,  8, -1, -1, -1, -1, -1,  5, -1,  3,
				-1, -1, -1, 15, -1, 12, -1,  9, -1, -1, -1,  6, -1,  3, -1,  0
			);
			__m256i mmask=_mm256_set_epi32(0, 0, 0, 0xFFFF, 0, 0, 0, 0xFFFF);
			memset(fc, 0, sizeof(fc));
			memset(mc, 0, sizeof(mc));
			while(ptr<end)
			{
				//process 11 pixels (33 SPs) at a time
				__m256i rW	=_mm256_loadu_si256((__m256i*)(ptr+0));
				__m256i gW	=_mm256_loadu_si256((__m256i*)(ptr+1));
				__m256i bW	=_mm256_loadu_si256((__m256i*)(ptr+2));
				__m256i r	=_mm256_loadu_si256((__m256i*)(ptr+3));
				__m256i g	=_mm256_loadu_si256((__m256i*)(ptr+4));
				__m256i b	=_mm256_loadu_si256((__m256i*)(ptr+5));
				ptr+=sizeof(__m256i[3]);
				rW	=_mm256_shuffle_epi8(rW	, shuf);
				gW	=_mm256_shuffle_epi8(gW	, shuf);
				bW	=_mm256_shuffle_epi8(bW	, shuf);
				r	=_mm256_shuffle_epi8(r	, shuf);
				g	=_mm256_shuffle_epi8(g	, shuf);
				b	=_mm256_shuffle_epi8(b	, shuf);

				__m256i d0=_mm256_sad_epu8(r, rW);
				__m256i d1=_mm256_sad_epu8(g, gW);
				__m256i d2=_mm256_sad_epu8(b, bW);
				mc[0]=_mm256_add_epi64(mc[0], d0);
				mc[1]=_mm256_add_epi64(mc[1], d1);
				mc[2]=_mm256_add_epi64(mc[2], d2);

				__m256i rg=_mm256_sub_epi16(r, g);
				__m256i gb=_mm256_sub_epi16(g, b);
				__m256i br=_mm256_sub_epi16(b, r);
				__m256i rgW=_mm256_sub_epi16(rW, gW);
				__m256i gbW=_mm256_sub_epi16(gW, bW);
				__m256i brW=_mm256_sub_epi16(bW, rW);
				__m256i d3=_mm256_sub_epi16(rg, rgW);
				__m256i d4=_mm256_sub_epi16(gb, gbW);
				__m256i d5=_mm256_sub_epi16(br, brW);
				d3=_mm256_abs_epi16(d3);
				d4=_mm256_abs_epi16(d4);
				d5=_mm256_abs_epi16(d5);
				d3=_mm256_hadd_epi16(d3, d3);
				d4=_mm256_hadd_epi16(d4, d4);
				d5=_mm256_hadd_epi16(d5, d5);
				d3=_mm256_hadd_epi16(d3, d3);
				d4=_mm256_hadd_epi16(d4, d4);
				d5=_mm256_hadd_epi16(d5, d5);
				d3=_mm256_hadd_epi16(d3, d3);
				d4=_mm256_hadd_epi16(d4, d4);
				d5=_mm256_hadd_epi16(d5, d5);
				d3=_mm256_and_si256(d3, mmask);
				d4=_mm256_and_si256(d4, mmask);
				d5=_mm256_and_si256(d5, mmask);
				//d3=_mm256_add_epi16(d3, _mm256_srli_epi64(d3, 16));
				//d4=_mm256_add_epi16(d4, _mm256_srli_epi64(d4, 16));
				//d5=_mm256_add_epi16(d5, _mm256_srli_epi64(d5, 16));
				//d3=_mm256_add_epi16(d3, _mm256_srli_epi64(d3, 32));
				//d4=_mm256_add_epi16(d4, _mm256_srli_epi64(d4, 32));
				//d5=_mm256_add_epi16(d5, _mm256_srli_epi64(d5, 32));
				//d3=_mm256_and_si256(d3, _mm256_set1_epi32(0xFFFF));
				//d4=_mm256_and_si256(d4, _mm256_set1_epi32(0xFFFF));
				//d5=_mm256_and_si256(d5, _mm256_set1_epi32(0xFFFF));

				//__m256i r2=_mm256_sub_epi16(r, _mm256_avg_epu16(g, b));
				//__m256i g2=_mm256_sub_epi16(g, _mm256_avg_epu16(b, r));
				//__m256i b2=_mm256_sub_epi16(b, _mm256_avg_epu16(r, g));
				//__m256i rW2=_mm256_sub_epi16(rW, _mm256_avg_epu16(gW, bW));
				//__m256i gW2=_mm256_sub_epi16(gW, _mm256_avg_epu16(bW, rW));
				//__m256i bW2=_mm256_sub_epi16(bW, _mm256_avg_epu16(rW, gW));
				//__m256i d6=_mm256_sub_epi16(r2, rW2);
				//__m256i d7=_mm256_sub_epi16(g2, gW2);
				//__m256i d8=_mm256_sub_epi16(b2, bW2);
				//d6=_mm256_abs_epi16(d6);
				//d7=_mm256_abs_epi16(d7);
				//d8=_mm256_abs_epi16(d8);
				//d6=_mm256_hadd_epi16(d6, d6);
				//d7=_mm256_hadd_epi16(d7, d7);
				//d8=_mm256_hadd_epi16(d8, d8);
				//d6=_mm256_hadd_epi16(d6, d6);
				//d7=_mm256_hadd_epi16(d7, d7);
				//d8=_mm256_hadd_epi16(d8, d8);
				//d6=_mm256_hadd_epi16(d6, d6);
				//d7=_mm256_hadd_epi16(d7, d7);
				//d8=_mm256_hadd_epi16(d8, d8);
				//d6=_mm256_and_si256(d6, mmask);
				//d7=_mm256_and_si256(d7, mmask);
				//d8=_mm256_and_si256(d8, mmask);

				mc[3]=_mm256_add_epi64(mc[3], d3);
				mc[4]=_mm256_add_epi64(mc[4], d4);
				mc[5]=_mm256_add_epi64(mc[5], d5);
				//mc[6]=_mm256_add_epi64(mc[6], d6);
				//mc[7]=_mm256_add_epi64(mc[7], d7);
				//mc[8]=_mm256_add_epi64(mc[8], d8);
				
				d0=_mm256_cmpeq_epi16(r, rW);
				d1=_mm256_cmpeq_epi16(g, gW);
				d2=_mm256_cmpeq_epi16(b, bW);
				d0=_mm256_abs_epi16(d0);
				d1=_mm256_abs_epi16(d1);
				d2=_mm256_abs_epi16(d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_and_si256(d0, mmask);
				d1=_mm256_and_si256(d1, mmask);
				d2=_mm256_and_si256(d2, mmask);
				fc[0]=_mm256_add_epi64(fc[0], d0);
				fc[1]=_mm256_add_epi64(fc[1], d1);
				fc[2]=_mm256_add_epi64(fc[2], d2);

				d0=_mm256_cmpeq_epi16(rg, rgW);
				d1=_mm256_cmpeq_epi16(gb, gbW);
				d2=_mm256_cmpeq_epi16(br, brW);
				d0=_mm256_abs_epi16(d0);
				d1=_mm256_abs_epi16(d1);
				d2=_mm256_abs_epi16(d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_hadd_epi16(d0, d0);
				d1=_mm256_hadd_epi16(d1, d1);
				d2=_mm256_hadd_epi16(d2, d2);
				d0=_mm256_and_si256(d0, mmask);
				d1=_mm256_and_si256(d1, mmask);
				d2=_mm256_and_si256(d2, mmask);
				fc[3]=_mm256_add_epi64(fc[3], d0);
				fc[4]=_mm256_add_epi64(fc[4], d1);
				fc[5]=_mm256_add_epi64(fc[5], d2);
			}
			for(int k=0;k<OCH_COUNT;++k)
			{
				flatcount[k]=
					+_mm256_extract_epi64(fc[k], 0)
					+_mm256_extract_epi64(fc[k], 1)
					+_mm256_extract_epi64(fc[k], 2)
					+_mm256_extract_epi64(fc[k], 3)
				;
				counters[k]=
					+_mm256_extract_epi64(mc[k], 0)
					+_mm256_extract_epi64(mc[k], 1)
					+_mm256_extract_epi64(mc[k], 2)
					+_mm256_extract_epi64(mc[k], 3)
				;
			}
#if 0
			int W[6]={0};
			while(ptr<end)
			{
				int
					r=ptr[0],
					g=ptr[1],
					b=ptr[2],
					rg=r-g,
					gb=g-b,
					br=b-r;
				ptr+=3;
				counters[0]+=abs(r-W[0]);
				counters[1]+=abs(g-W[1]);
				counters[2]+=abs(b-W[2]);
				counters[3]+=abs(rg-W[3]);//abs(r-g-rW+gW)
				counters[4]+=abs(gb-W[4]);//abs(g-b-gW+bW)
				counters[5]+=abs(br-W[5]);//abs(b-r-bW+rW)
				W[0]=r;
				W[1]=g;
				W[2]=b;
				W[3]=rg;
				W[4]=gb;
				W[5]=br;
			}
#endif
			long long minerr=0;
			for(int kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
				long long currerr=
					+counters[rct[0]]
					+counters[rct[1]]
					+counters[rct[2]]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
			const unsigned char *rct=rct_indices[bestrct];
			long long flatarea=
				+flatcount[rct[0]]
				+flatcount[rct[1]]
				+flatcount[rct[2]]
			;
			flat=flatarea>usize/3;
#ifdef LOUD
			long long currerr=
				+counters[rct[0]]
				+counters[rct[1]]
				+counters[rct[2]]
			;
			printf("%s\n", srcfn);
			printf("flat = %d\n", flat);
			printf("F/U = %12lld / %12td = %12.6lf\n", flatarea, usize, (double)flatarea/usize);
			printf("E/U = %12lld / %12td = %12.6lf\n", currerr, usize, (double)currerr/usize);
#endif
#if 0
			int rowstride=iw*3;
			for(int ky=0;ky<ih;ky+=ANALYSIS_YSTRIDE)
			{
				const unsigned char *ptr=srcptr+rowstride*ky, *end=ptr+rowstride;
				int prev2[6]={0}, prev[6]={0};
				while(ptr<end)
				{
					int
						r=ptr[0],
						g=ptr[1],
						b=ptr[2],
						rg=r-g,
						gb=g-b,
						br=b-r;
					ptr+=3*ANALYSIS_XSTRIDE;
					counters[OCH_COUNT*0+0]+=abs(r-prev[0]);
					counters[OCH_COUNT*0+1]+=abs(g-prev[1]);
					counters[OCH_COUNT*0+2]+=abs(b-prev[2]);
					counters[OCH_COUNT*0+3]+=abs(rg-prev[3]);
					counters[OCH_COUNT*0+4]+=abs(gb-prev[4]);
					counters[OCH_COUNT*0+5]+=abs(br-prev[5]);
					counters[OCH_COUNT*1+0]+=abs(r-(2*prev[0]-prev2[0]));
					counters[OCH_COUNT*1+1]+=abs(g-(2*prev[1]-prev2[1]));
					counters[OCH_COUNT*1+2]+=abs(b-(2*prev[2]-prev2[2]));
					counters[OCH_COUNT*1+3]+=abs(rg-(2*prev[3]-prev2[3]));
					counters[OCH_COUNT*1+4]+=abs(gb-(2*prev[4]-prev2[4]));
					counters[OCH_COUNT*1+5]+=abs(br-(2*prev[5]-prev2[5]));
					prev2[0]=prev[0];
					prev2[1]=prev[1];
					prev2[2]=prev[2];
					prev2[3]=prev[3];
					prev2[4]=prev[4];
					prev2[5]=prev[5];
					prev[0]=r;
					prev[1]=g;
					prev[2]=b;
					prev[3]=rg;
					prev[4]=gb;
					prev[5]=br;
				}
			}
			for(int kx=0;kx<iw;kx+=ANALYSIS_XSTRIDE)
			{
				const unsigned char *ptr=srcptr+kx, *end=srcptr+usize;
				int stride=ANALYSIS_YSTRIDE*rowstride;
				int prev2[6]={0}, prev[6]={0};
				while(ptr<end)
				{
					int
						r=ptr[0],
						g=ptr[1],
						b=ptr[2],
						rg=r-g,
						gb=g-b,
						br=b-r;
					ptr+=stride;
					counters[OCH_COUNT*2+0]+=abs(r-prev[0]);
					counters[OCH_COUNT*2+1]+=abs(g-prev[1]);
					counters[OCH_COUNT*2+2]+=abs(b-prev[2]);
					counters[OCH_COUNT*2+3]+=abs(rg-prev[3]);
					counters[OCH_COUNT*2+4]+=abs(gb-prev[4]);
					counters[OCH_COUNT*2+5]+=abs(br-prev[5]);
					counters[OCH_COUNT*3+0]+=abs(r-(2*prev[0]-prev2[0]));
					counters[OCH_COUNT*3+1]+=abs(g-(2*prev[1]-prev2[1]));
					counters[OCH_COUNT*3+2]+=abs(b-(2*prev[2]-prev2[2]));
					counters[OCH_COUNT*3+3]+=abs(rg-(2*prev[3]-prev2[3]));
					counters[OCH_COUNT*3+4]+=abs(gb-(2*prev[4]-prev2[4]));
					counters[OCH_COUNT*3+5]+=abs(br-(2*prev[5]-prev2[5]));
					prev2[0]=prev[0];
					prev2[1]=prev[1];
					prev2[2]=prev[2];
					prev2[3]=prev[3];
					prev2[4]=prev[4];
					prev2[5]=prev[5];
					prev[0]=r;
					prev[1]=g;
					prev[2]=b;
					prev[3]=rg;
					prev[4]=gb;
					prev[5]=br;
				}
			}
			long long minerr=0;
			for(int kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
				long long currerrW0=
					+counters[OCH_COUNT*0+rct[0]]
					+counters[OCH_COUNT*0+rct[1]]
					+counters[OCH_COUNT*0+rct[2]]
				;
				long long currerrW1=
					+counters[OCH_COUNT*1+rct[0]]
					+counters[OCH_COUNT*1+rct[1]]
					+counters[OCH_COUNT*1+rct[2]]
				;
				long long currerrN0=
					+counters[OCH_COUNT*2+rct[0]]
					+counters[OCH_COUNT*2+rct[1]]
					+counters[OCH_COUNT*2+rct[2]]
				;
				long long currerrN1=
					+counters[OCH_COUNT*3+rct[0]]
					+counters[OCH_COUNT*3+rct[1]]
					+counters[OCH_COUNT*3+rct[2]]
				;
				if(!kt||minerrW0>currerrW0)
				{
					minerrW0=currerrW0;
					bestrct=kt;
				}
				if(!kt||minerrW1>currerrW1)
				{
					minerrW1=currerrW1;
					bestrct=kt;
				}
				if(!kt||minerrN0>currerrN0)
				{
					minerrN0=currerrN0;
					bestrct=kt;
				}
				if(!kt||minerrN1>currerrN1)
				{
					minerrN1=currerrN1;
					bestrct=kt;
				}
			}
#endif
		}
		if(nthreads0&1)
			flat=0;
		else if(nthreads0&2)
			flat=1;

		int flag=bestrct<<1|flat;
		low+=range*flag>>5;
		range=(range>>5)-1;
	}
	else
	{
		dstbuf=(unsigned char*)malloc(usize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		dstptr=dstbuf;

		code=*(unsigned*)srcptr;
		srcptr+=4;
		code=code<<4|*(unsigned*)srcptr;
		srcptr+=4;

		int flag=(int)((((code-low+1)<<5)-1)/range);
		low+=range*flag>>5;
		range=(range>>5)-1;

		bestrct=flag>>1;
		flat=flag&1;
	}
#if 0
#define NPREDS 8
#define WP_PREDLIST\
	WP_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WP_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WP_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WP_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WP_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WP_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WP_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WP_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#endif
#if 1
#define NPREDS 8
#define WP_PREDLIST\
	WP_PRED(260, N)\
	WP_PRED(260, W)\
	WP_PRED(150, 3*(N-NN)+NNN)\
	WP_PRED(150, 3*(W-WW)+WWW)\
	WP_PRED(170, W+NE-N)\
	WP_PRED(100, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)\
	WP_PRED(170, N+W-NW)\
	WP_PRED(160, N+NE-NNE)
#endif

	ALIGN(32) int wperrors[3][NPREDS]={0};
	int psize=(iw+16LL)*(int)sizeof(short[4*3*3]);//4 padded rows * 3 channels * {pixels, errors, nbypass}
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	int
		yidx=rct_indices[bestrct][3+0],
		uidx=rct_indices[bestrct][3+1],
		vidx=rct_indices[bestrct][3+2],
		uhelpidx=rct_indices[bestrct][6+0],
		vhelpidx=rct_indices[bestrct][6+1];
	//for(int kc=0;kc<6;++kc)
	//{
	//	for(int k=0;k<1<<PROBBITS_USE;++k)
	//	{
	//	//	//squash(x) = sigmoid(x) = 1/(1+exp(-x))
	//		int input=(k-(1<<PROBBITS_USE>>1))<<(24-PROBBITS_USE)|1<<(24-PROBBITS_USE)>>1;
	//		probtable[kc][k]=(int)((0x1000000ULL<<PROBBITS_ADJUST)/(0x1000000+exp2_fix24(-(input*19>>2))));
	//
	//	//	int input=(k-(1<<PROBBITS_USE>>1))<<(PROBBITS_ADJUST-PROBBITS_USE)|1<<(PROBBITS_ADJUST-PROBBITS_USE)>>1;
	//	//	probtable[kc][k]=(int)((0x1000000ULL<<PROBBITS_ADJUST)/(0x1000000+exp2_fix24(-(input<<(23-PROBBITS_USE)))));
	//	}
	//	//	probtable[kc][k]=k<<(PROBBITS_ADJUST-PROBBITS_USE)|1<<(PROBBITS_ADJUST-PROBBITS_USE)>>1;
	//}
//	memset(stats1, 0, sizeof(stats1));
	FILLMEM((unsigned*)stats1, 1<<PROBBITS_STORE>>1, sizeof(stats1), sizeof(int));
	C12Info info=
	{
		low, range, code,
		fwd?dstptr:srcptr,
		fwd,
		0, 0,
	//	0,
	};
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*3*3,
		};
		char yuv[4]={0};
		int errors[3]={0};
//		ALIGN(32) static const float weights[]=
//		{
//#define WG_PRED(WEIGHT, PRED) WEIGHT,
//			WG_PREDLIST
//#undef  WG_PRED
//		};
//		__m256 mw=_mm256_load_ps(weights);
		//ALIGN(32) float fpreds[8]={0}, fcoeff[8]={0};
		//ptrdiff_t checkpoint=(ptrdiff_t)info.ptr;
		ALIGN(32) int preds[NPREDS]={0};
		int hipred=0, pred=0, ctx=0;
		for(int kx=0;kx<iw;++kx)
		{
			//if(ky==ih/2&&kx==iw/2)//
			//	printf("");//

			short *curr=rows[0]+0*3*3;
			if(fwd)
			{
				yuv[0]=srcptr[yidx]-128;
				yuv[1]=srcptr[uidx]-128;
				yuv[2]=srcptr[vidx]-128;
				srcptr+=3;
				yuv[2]-=yuv[vhelpidx];
				yuv[1]-=yuv[uhelpidx];
				curr[0]=yuv[0];
				curr[3]=yuv[1];
				curr[6]=yuv[2];
			}
			for(int kc=0;kc<3;++kc)
			{
				int//		rows[-Y][E+X*3*3],
					NNN	=rows[3][0+0*3*3],
					NN	=rows[2][0+0*3*3],
					NNE	=rows[2][0+1*3*3],
					NW	=rows[1][0-1*3*3],
					N	=rows[1][0+0*3*3],
					NE	=rows[1][0+1*3*3],
					NEE	=rows[1][0+2*3*3],
					NEEE	=rows[1][0+3*3*3],
					NEEEE	=rows[1][0+4*3*3],
					WWWW	=rows[0][0-4*3*3],
					WWW	=rows[0][0-3*3*3],
					WW	=rows[0][0-2*3*3],
					W	=rows[0][0-1*3*3],
					eNNN	=rows[3][2+0*3*3],
					eNN	=rows[2][2+0*3*3],
					eNNE	=rows[2][2+1*3*3],
					eNW	=rows[1][2-1*3*3],
					eN	=rows[1][2+0*3*3],
					eNE	=rows[1][2+1*3*3],
					eWWW	=rows[0][2-3*3*3],
					eWW	=rows[0][2-2*3*3],
					eW	=rows[0][2-1*3*3],
					bW	=rows[0][1-1*3*3];
				(void)eNNN;
				(void)eNN;
				(void)eNNE;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eWWW;
				(void)eWW;
				(void)eW;
				if(flat)
				{
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					pred=N+W-NW;
					CLAMP2(pred, vmin, vmax);
					hipred=pred<<4;
					ctx=pred&255;
				}
				else
				{
					static const int wg_weights[]=
					{
#define WP_PRED(WEIGHT, EXPR) WEIGHT,
						WP_PREDLIST
						//WG_PREDLIST0
						//WG_PREDLIST1
						//WG_PREDLIST2
#undef  WP_PRED
					};
					int j=0;
#define WP_PRED(WEIGHT, PRED) preds[j++]=PRED;
					WP_PREDLIST
#undef  WP_PRED
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
					static int divlookup[DIVLUTSIZE]={0};
					if(!*divlookup)
					{
						for(int k=0;k<DIVLUTSIZE;++k)
							divlookup[k]=(1<<NUMBITS)/(k+1);
					}
					unsigned coeff[NPREDS], wsum=0;
					int ipred=0;
					for(int kp=0;kp<NPREDS;++kp)
					{
						//			1	1
						//		2	2	1
						//	1	2	?
						coeff[kp]=//maxerror = 18<<depth
							+wperrors[kc][kp]		//+I	*8
						//	+cN	[kp-1*4*NPREDS]*2	//+eNW	*2
						//	+cN	[kp+0*4*NPREDS]*2	//+eN	*2
						//	+ccurr	[kp-1*4*NPREDS]*2	//+eW	*2
						//	+ccurr	[kp+0*4*NPREDS]		//+eNN	*1
						//	+ccurr	[kp+1*4*NPREDS]		//+eNNE	*1
						//	+cN	[kp+1*4*NPREDS]		//+eNE	*1
						//	+ccurr	[kp-2*4*NPREDS]		//+eWW	*1
						;
					}
					int sh[NPREDS];
					sh[0]=31-(DENBITS-1)-_lzcnt_u32(coeff[0]+1);//invert errros
					sh[1]=31-(DENBITS-1)-_lzcnt_u32(coeff[1]+1);
					sh[2]=31-(DENBITS-1)-_lzcnt_u32(coeff[2]+1);
					sh[3]=31-(DENBITS-1)-_lzcnt_u32(coeff[3]+1);
					sh[4]=31-(DENBITS-1)-_lzcnt_u32(coeff[4]+1);
					sh[5]=31-(DENBITS-1)-_lzcnt_u32(coeff[5]+1);
					sh[6]=31-(DENBITS-1)-_lzcnt_u32(coeff[6]+1);
					sh[7]=31-(DENBITS-1)-_lzcnt_u32(coeff[7]+1);
					if(sh[0]<0)sh[0]=0;
					if(sh[1]<0)sh[1]=0;
					if(sh[2]<0)sh[2]=0;
					if(sh[3]<0)sh[3]=0;
					if(sh[4]<0)sh[4]=0;
					if(sh[5]<0)sh[5]=0;
					if(sh[6]<0)sh[6]=0;
					if(sh[7]<0)sh[7]=0;
					coeff[0]=(wg_weights[0]*divlookup[coeff[0]>>sh[0]]>>sh[0])+(1<<DENBITS>>2)/NPREDS;
					coeff[1]=(wg_weights[1]*divlookup[coeff[1]>>sh[1]]>>sh[1])+(1<<DENBITS>>2)/NPREDS;
					coeff[2]=(wg_weights[2]*divlookup[coeff[2]>>sh[2]]>>sh[2])+(1<<DENBITS>>2)/NPREDS;
					coeff[3]=(wg_weights[3]*divlookup[coeff[3]>>sh[3]]>>sh[3])+(1<<DENBITS>>2)/NPREDS;
					coeff[4]=(wg_weights[4]*divlookup[coeff[4]>>sh[4]]>>sh[4])+(1<<DENBITS>>2)/NPREDS;
					coeff[5]=(wg_weights[5]*divlookup[coeff[5]>>sh[5]]>>sh[5])+(1<<DENBITS>>2)/NPREDS;
					coeff[6]=(wg_weights[6]*divlookup[coeff[6]>>sh[6]]>>sh[6])+(1<<DENBITS>>2)/NPREDS;
					coeff[7]=(wg_weights[7]*divlookup[coeff[7]>>sh[7]]>>sh[7])+(1<<DENBITS>>2)/NPREDS;
					wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];//invert coeff sum
					sh[0]=FLOOR_LOG2(wsum)-(DENBITS-2);
					wsum=0;
					coeff[0]>>=sh[0];
					coeff[1]>>=sh[0];
					coeff[2]>>=sh[0];
					coeff[3]>>=sh[0];
					coeff[4]>>=sh[0];
					coeff[5]>>=sh[0];
					coeff[6]>>=sh[0];
					coeff[7]>>=sh[0];
					wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];
					ipred=(wsum>>1)
						+coeff[0]*preds[0]
						+coeff[1]*preds[1]
						+coeff[2]*preds[2]
						+coeff[3]*preds[3]
						+coeff[4]*preds[4]
						+coeff[5]*preds[5]
						+coeff[6]*preds[6]
						+coeff[7]*preds[7]
					;
					hipred=ipred*divlookup[wsum-1];
					pred=hipred>>NUMBITS;
				//	ctx=hipred>>(NUMBITS-2)&1023;
					hipred>>=NUMBITS-4;
#if 0
					__m256i me=_mm256_load_si256((__m256i*)wperrors+kc);
					__m256i mp=_mm256_load_si256((__m256i*)preds);
					me=_mm256_add_epi32(me, _mm256_set1_epi32(1));
					__m256 fe=_mm256_cvtepi32_ps(me);
					__m256 fp=_mm256_cvtepi32_ps(mp);
					fe=_mm256_rcp_ps(fe);
					fe=_mm256_mul_ps(fe, mw);
					fp=_mm256_mul_ps(fp, fe);
					__m128 fe4=_mm_add_ps(_mm256_castps256_ps128(fe), _mm256_extractf128_ps(fe, 1));
					__m128 fp4=_mm_add_ps(_mm256_castps256_ps128(fp), _mm256_extractf128_ps(fp, 1));
					fp4=_mm_hadd_ps(fp4, fp4);
					fe4=_mm_hadd_ps(fe4, fe4);
					fp4=_mm_hadd_ps(fp4, fp4);
					fe4=_mm_hadd_ps(fe4, fe4);
					fp4=_mm_div_ss(fp4, fe4);
					fp4=_mm_castsi128_ps(_mm_add_epi32(_mm_castps_si128(fp4), _mm_set1_epi32(4<<23)));//multiply by 16
					hipred=_mm_cvtt_ss2si(fp4);
					pred=(hipred+(1<<4>>1))>>4;
					ctx=hipred>>2&1023;
#endif
#if 0
					int vmax=N, vmin=W;
					if(N<W)
						vmin=N, vmax=W;
					if(vmin>NW)vmin=NW;
					if(vmax<NW)vmax=NW;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEE)vmin=NEE;
					if(vmax<NEE)vmax=NEE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(cpred, vmin, vmax);
#endif
				}
				int nbypass=FLOOR_LOG2(bW+1);
				ctx=nbypass;
				nbypass-=6;
				if(nbypass<0)
					nbypass=0;
				errors[kc]=0;
				int nzeros=-1, bypass=0, tidx=0;
				if(fwd)
				{
					errors[kc]=(char)(yuv[kc]-pred);
					errors[kc]=errors[kc]<<1^errors[kc]>>31;
					nzeros=errors[kc]>>nbypass;
					bypass=errors[kc]&((1<<nbypass)-1);
				}
				unsigned *statsptr=stats1[kc][(pred+128)&255][ctx];
			//	info.probtable=probtable[kc];
				do
				{
#define GRSHLIST\
	GRSH(6)\
	GRSH(6)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)
					static const int sh[]=
					{
#define GRSH(X) X,
						GRSHLIST
#undef  GRSH
					};
					static const int rcon[]=
					{
#define GRSH(X) 1<<(X)>>1,
						GRSHLIST
#undef  GRSH
					};
					unsigned *p0ptr=statsptr+tidx;
					if(fwd)
						info.bit=nzeros--<=0;
					info.p0=*p0ptr;
					codebit(&info, rcon[tidx], sh[tidx]);
					*p0ptr=info.p0;
					if(!fwd)
						++nzeros;
					++tidx;
				}while(!info.bit&&tidx<GRLIMIT);
				if(tidx==GRLIMIT)
				{
					nbypass=8;
					if(fwd)
						bypass=errors[kc];
				}
			//	info.probtable=probtable[kc+3];
				
				//if(ky==799&&kx==38&&kc==2)//
				//	printf("");//

				for(int kb=nbypass-1, tidx=1;kb>=0;--kb)
				{
					unsigned *p0ptr=statsptr+GRLIMIT+8-nbypass+kb;
					if(fwd)
						info.bit=bypass>>kb&1;
					info.p0=*p0ptr;
					codebit(&info, 1<<9>>1, 9);
					*p0ptr=info.p0;
					if(!fwd)
						bypass|=info.bit<<kb;
					tidx=tidx*2+info.bit;
				}
				if(!fwd)
				{
					errors[kc]=bypass;
					if(tidx<GRLIMIT)
						errors[kc]|=nzeros<<nbypass;
					errors[kc]=errors[kc]>>1^-(errors[kc]&1);
					yuv[kc]=(char)(errors[kc]+pred);
				}
#if 0
				if(fwd)
				{
					errors[kc]=(char)(yuv[kc]-pred);
					errors[kc]=errors[kc]<<1^errors[kc]>>31;
					int nzeros=errors[kc]>>nbypass, bypass=errors[kc]&((1<<nbypass)-1);
					int tidx=0;
					info.bit=0;
					while(nzeros)
					{
						unsigned short *p0ptr=(unsigned short*)&stats1[kc][ctx][tidx];
						info.p0=*p0ptr;
						codebit(&info);
						*p0ptr=info.p0;
						++tidx;
						--nzeros;
					}
					info.bit=1;
					{
						unsigned short *p0ptr=(unsigned short*)&stats1[kc][ctx][tidx];
						info.p0=*p0ptr;
						codebit(&info);
						*p0ptr=info.p0;
						++tidx;
					}
					for(int kb=nbypass-1;kb>=0;--kb)
					{
						unsigned short *p0ptr=(unsigned short*)&stats1[kc][ctx][256-nbypass+kb];
						info.bit=bypass>>kb&1;
						info.p0=*p0ptr;
						codebit(&info);
						*p0ptr=info.p0;
					}
				}
				else
				{
				}
#endif
#if 0
				if(fwd)
					errors[kc]=yuv[kc]-pred;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int kb=7, tidx=1;kb>=0;--kb)
				{
					unsigned short *p0ptr=(unsigned short*)&stats1[kc][ctx][tidx];
					if(fwd)
						info.bit=errors[kc]>>kb&1;
					info.p0=*p0ptr;
					codebit(&info);
					*p0ptr=info.p0;
					if(!fwd)
						errors[kc]|=info.bit<<kb;
					tidx=tidx*2+info.bit;
				}
				if(!fwd)
					yuv[kc]=(char)(errors[kc]+pred);
#endif

				int e=(yuv[kc]<<4)-hipred;
				e=e<<1^e>>31;
				rows[0][1+0*3*3]=(2*bW+(e<<2)+rows[1][1+3*4*3])>>2;
				rows[0][2+0*3*3]=yuv[kc]-pred;
				if(!flat)
				{
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
					for(int k=0;k<8;++k)
					{
						e=yuv[kc]-preds[k];
						e=e<<1^e>>31;
						wperrors[kc][k]+=((e<<8)-wperrors[kc][k]+(1<<3>>1))>>3;
					//	wperrors[kc][k]+=((abs(yuv[kc]-preds[k])<<3)-wperrors[kc][k]+(1<<3>>1))>>3;
					}
				}
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
			}
			if(!fwd)
			{
				curr[0]=(char)yuv[0];
				curr[3]=(char)yuv[1];
				curr[6]=(char)yuv[2];
				yuv[1]+=yuv[uhelpidx];
				yuv[2]+=yuv[vhelpidx];
				dstptr[yidx]=(unsigned char)(yuv[0]+128);
				dstptr[uidx]=(unsigned char)(yuv[1]+128);
				dstptr[vidx]=(unsigned char)(yuv[2]+128);
				dstptr+=3;

				guide_check(dstbuf, kx, ky);//
			}
		}
	//	checkpoint=(ptrdiff_t)info.ptr-checkpoint;
	//	printf("%5td / %5d\n", checkpoint, 3*iw);//
	}
	free(srcbuf);
	_mm_free(pixels);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
			//flush
			*(unsigned*)info.ptr=(unsigned)(info.low>>32);
			info.ptr+=4;
			*(unsigned*)info.ptr=(unsigned)info.low;
			info.ptr+=4;

			fwrite("12", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(dstbuf, 1, info.ptr-dstbuf, fdst);
		}
		else
		{
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}
#if defined _MSC_VER && defined LOUD
	t=time_sec()-t;
	if(fwd)
	{
		ptrdiff_t csize=info.ptr-dstbuf+10;
		printf("%12td/%12td  %12.6lf  %8.4lf%%\n", csize, usize, (double)usize/csize, (double)csize*100/usize);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
	free(dstbuf);
	return 0;
}