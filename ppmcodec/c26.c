#include"codec.h"
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define BYPASS_FLAGS

#ifdef _MSC_VER
	#define LOUD
//	#define ENABLE_GUIDE
#ifndef BYPASS_FLAGS
//	#define AC_VALIDATE
#endif
#endif


#ifdef AC_VALIDATE
#define AC_IMPLEMENTATION
#include"entropy.h"
#endif
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(const unsigned char *image, int iw, int ih)
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
static void guide_check(const unsigned char *image, int kx, int ky)
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

typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_R2,
	OCH_G2,
	OCH_B2,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,

	OCH_COUNT=9,
} OCHIndex;
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_HELP_U,
	II_HELP_V0,
	II_HELP_V1,

	II_COUNT,
} RCTInfoIdx;
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 2)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  2, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  2, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 2)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 2)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  2, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	2,  2, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	2,  0, 2)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	2,  0, 2)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	2,  2, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	2,  0, 2)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	2,  0, 2)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	2,  2, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	2,  0, 2)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	2,  0, 2)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  1, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	2,  1, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  1, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	2,  1, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  1, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	2,  1, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	2,  1, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	2,  1, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	2,  1, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};
unsigned short stats1[3][256][256]={0};
//unsigned short stats0[3][256]={0};
int c26_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double t=time_sec();
	long long bitcount[3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize;
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<0)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
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
		srcbuf[srcsize]=0;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('2'|'6'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	if(fwd)//encode
	{
		if(*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++ - '0';
		while(*srcptr==' ')
			++srcptr;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++ - '0';
		if(memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		srcptr+=5;
	}
	else//decode
	{
		if(srcptr+4*3>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	ptrdiff_t usize=(ptrdiff_t)3*iw*ih;
#ifdef BYPASS_FLAGS
	int nbits=fwd?64:0;
	unsigned long long cache=0;
#else
#ifdef AC_VALIDATE
	unsigned long long lo0=0, r0=0;
#endif
	unsigned long long low=0, range=0xFFFFFFFF, code=0;
#endif
	int bestrct=0;
	ptrdiff_t dstbufsize=0;
	unsigned char *dstbuf=0, *image=0, *streamptr=0;
	if(fwd)
	{
		dstbufsize=(ptrdiff_t)4*iw*ih;
		dstbuf=malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		streamptr=dstbuf;
		image=srcptr;
		guide_save(image, iw, ih);

		long long counters[9]={0};
		int prev[9]={0};
		while(srcptr<srcend)
		{
			int
				r=srcptr[0],
				g=srcptr[1],
				b=srcptr[2],
				rg=r-g,
				gb=g-b,
				br=b-r;
			srcptr+=3;
			counters[0]+=abs(r -prev[0]);
			counters[1]+=abs(g -prev[1]);
			counters[2]+=abs(b -prev[2]);
			counters[3]+=abs(rg-prev[3]);
			counters[4]+=abs(gb-prev[4]);
			counters[5]+=abs(br-prev[5]);
			prev[0]=r;
			prev[1]=g;
			prev[2]=b;
			prev[3]=rg;
			prev[4]=gb;
			prev[5]=br;
		}
		long long minerr=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const unsigned char *rct=rct_combinations[kt];
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
		//no renorm here
#ifdef BYPASS_FLAGS
		nbits-=5;
		cache|=(unsigned long long)bestrct<<nbits;
#else
		low+=range*bestrct/RCT_COUNT;
		range=range/RCT_COUNT-1;
#endif
	}
	else
	{
		dstbufsize=usize;
		dstbuf=malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		streamptr=srcptr;
		image=dstbuf;
		
#ifdef BYPASS_FLAGS
		cache=*(unsigned*)streamptr;
		streamptr+=4;
		nbits-=5;
		bestrct=cache>>nbits&((1ULL<<5)-1);
#else
		code=*(unsigned*)streamptr;
		streamptr+=4;
		code=code<<32|*(unsigned*)streamptr;
		streamptr+=4;
		
		//no renorm here
		bestrct=(int)(((code-low)*RCT_COUNT+RCT_COUNT-1)/range);
		low+=range*bestrct/RCT_COUNT;
		range=range/RCT_COUNT-1;
#endif
	}
	const unsigned char *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	int coeffs[]=
	{
		0, combination[II_HELP_U], combination[II_HELP_V0],
		0, 0, combination[II_HELP_V1],
	};
	int psize=(int)sizeof(short[4*3*2])*(iw+16);//4 padded rows * 3 channels max * {pixels, errors}
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	int paddedwidth=iw+16;
	unsigned char *ptr=image;
	FILLMEM((unsigned short*)stats1, 0x8000, sizeof(stats1), sizeof(short));
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*3*2,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*3*2,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*3*2,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*3*2,
		};
		int yuv[3]={0};
		//int NW[3]={0}, W[3]={0};
		unsigned char error=0, bit=0;
		for(int kx=0;kx<iw;++kx, ptr+=3)
		{
			//short *N=rows[1]+0*3*2;
			if(fwd)
			{
				yuv[0]=ptr[yidx]-128;
				yuv[1]=ptr[uidx]-128;
				yuv[2]=ptr[vidx]-128;
			}
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				int offset=(coeffs[kc+0]*yuv[0]+coeffs[kc+3]*yuv[1])>>1;
				int
					NW	=rows[1][kc-1*3*2+0],
					N	=rows[1][kc+0*3*2+0],
					W	=rows[0][kc-1*3*2+0],
					eN	=rows[1][kc+0*3*2+3],
					eW	=rows[0][kc-1*3*2+3];
				(void)eN;
				(void)eW;
#ifndef BYPASS_FLAGS
				int ctx=(eN+eW)&255;
				unsigned short *curr_stats=stats1[kc][ctx];
#endif
				int pred=N+W-NW;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;

				//int pred=N[kc]+W[kc]-NW[kc];
				//int vmax=N[kc], vmin=W[kc];
				//if(N[kc]<W[kc])vmin=N[kc], vmax=W[kc];

				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, -128, 127);
				if(fwd)
					error=(unsigned char)(yuv[kc]-pred+128);
				else
					error=0;
#if 1
				int sym=128, step=0;
				for(int it=0;;++it)
				{
#ifdef BYPASS_FLAGS
					if(nbits<32)
					{
						if(fwd)
							*(unsigned*)streamptr=(unsigned)(cache>>32);
						cache<<=32;
						if(!fwd)
							cache|=*(unsigned*)streamptr;
						streamptr+=4;
						nbits+=32;
					}
					--nbits;

					bit=cache>>nbits&1;
					if(fwd)
						bit=error>=sym;
					cache|=(unsigned long long)bit<<nbits;
					//if(fwd)
					//{
					//	bit=error>=sym;
					//	cache|=(unsigned long long)bit<<nbits;
					//}
					//else
					//{
					//	bit=cache>>nbits&1;
					//}
#else
					int p0=curr_stats[sym-1];
#ifdef AC_VALIDATE
					lo0=low;
					r0=range;
#endif
					if(range<0x10000)
					{
						if(fwd)
							*(unsigned*)streamptr=(unsigned)(low>>32);
						else
							code=code<<32|*(unsigned*)streamptr;
						streamptr+=4;
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					unsigned long long r2=range*p0>>16;
					unsigned long long mid=low+r2;
					bit=fwd?error>=sym:code>=mid;
					if(bit)
					{
						low=mid;
						range-=r2;
					}
					else
						range=r2-1;
#ifdef AC_VALIDATE
					if(fwd)
						acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, lo0, lo0+r0, low, low+range, 0, 0);
					else
						acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
					int sh=7-(it>1)-(it>3);
					p0+=((!bit<<16)-p0+(1<<sh>>1))>>sh;
					curr_stats[sym-1]=p0;
#endif
#ifdef LOUD
					++bitcount[kc];
#endif

					int step2=2*bit-1;

					if(!(step+step2))//stop bit? (hard to predict branch)
						break;
					step=step2;
					sym+=step;
					if((unsigned)(sym-1)>=255)
						break;
				}
				sym-=sym>=128;
#ifdef _DEBUG
				if(fwd&&error!=sym)//
					LOG_ERROR("Algorithm error");
#endif
				error=(unsigned char)sym;
#endif
#if 0
				unsigned long long r2=0;
				for(int kb=7, tidx=1;kb>=0;--kb)//uniform binary search (inefficient)
				{
					int p0=stats0[kc][tidx];
					if(fwd)
					{
						bit=error>>kb&1;
						
#ifdef AC_VALIDATE
						lo0=low;
						r0=range;
#endif
						if(range<0x10000)
						{
							*(unsigned*)streamptr=(unsigned)(low>>32);
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						r2=range*p0>>16;
						if(bit)
						{
							low+=r2;
							range-=r2;
						}
						else
							range=r2-1;
#ifdef AC_VALIDATE
						acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, lo0, lo0+r0, low, low+range, 0, 0);
#endif
					}
					else
					{
#ifdef AC_VALIDATE
						lo0=low;
						r0=range;
#endif
						if(range<0x10000)
						{
							code=code<<32|*(unsigned*)streamptr;
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						r2=range*p0>>16;
						unsigned long long mid=low+r2;
						bit=code>=mid;
						if(bit)
						{
							low=mid;
							range-=r2;
						}
						else
							range=r2-1;
#ifdef AC_VALIDATE
						acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif

						error|=bit<<kb;
					}
					p0+=((!bit<<16)-p0+(1<<5>>1))>>5;
					stats0[kc][tidx]=p0;
					tidx=2*tidx+bit;
				}
#endif
				if(!fwd)
					yuv[kc]=(char)(error-128+pred);
				rows[0][kc+0]=yuv[kc]-offset;
				rows[0][kc+3]=yuv[kc]-pred;
			}
			if(!fwd)
			{
				ptr[yidx]=yuv[0]+128;
				ptr[uidx]=yuv[1]+128;
				ptr[vidx]=yuv[2]+128;
				guide_check(image, kx, ky);
			}

			//NW[0]=N[0];
			//NW[1]=N[1];
			//NW[2]=N[2];
			//W[0]=rows[0][0];
			//W[1]=rows[0][1];
			//W[2]=rows[0][2];

			rows[0]+=3*2;
			rows[1]+=3*2;
			rows[2]+=3*2;
			rows[3]+=3*2;
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", fdst);
			return 1;
		}
		if(fwd)
		{
#ifdef BYPASS_FLAGS
			*(unsigned*)streamptr=(unsigned)(cache>>32);
			streamptr+=4;
			*(unsigned*)streamptr=(unsigned)cache;
			streamptr+=4;
#else
			*(unsigned*)streamptr=(unsigned)(low>>32);
			streamptr+=4;
			*(unsigned*)streamptr=(unsigned)low;
			streamptr+=4;
#endif

			fwrite("26", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(dstbuf, 1, streamptr-dstbuf, fdst);
		}
		else
		{
			srcsize=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
			srcsize+=usize;
		}
		fclose(fdst);
	}
	free(pixels);
	free(dstbuf);
	free(srcbuf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%s\n", srcfn);
		printf("Flags: (usize %td)\n", srcsize);
		printf("Y %12.2lf\n", bitcount[0]/8.);
		printf("U %12.2lf\n", bitcount[1]/8.);
		printf("V %12.2lf\n", bitcount[2]/8.);
		printf("%8td/%8td bytes\n", 2+4+4+streamptr-dstbuf, srcsize);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, srcsize/(t*1024*1024));
#endif

	(void)nthreads0;
	(void)rct_names;
	return 0;
}