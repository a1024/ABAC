#include"codec.h"
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


//	#define ENABLE_GUIDE
//	#define AC_VALIDATE
//	#define LOUD

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

	OCH_COUNT=6,
} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},
};
static unsigned short stats1[3][256][256];
int c12_codec(const char *srcfn, const char *dstfn)
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
	if(srcsize<1)
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
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr==' ')
			++srcptr;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(iw<1||ih<1||memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Invalid file");
			return 1;
		}
		srcptr+=5;
	}
	else
	{
		if(memcmp(srcptr, "CH", 2))
		{
			LOG_ERROR("Invalid file");
			return 1;
		}
		srcptr+=2;
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	
	ptrdiff_t res=(ptrdiff_t)iw*ih, usize=res*3;
	int bestrct=0;
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
			long long counters[6]={0};
			int W[6]={0};
			for(int k=0;k<usize;k+=3)
			{
				int
					r=srcptr[k+0],
					g=srcptr[k+1],
					b=srcptr[k+2],
					rg=r-g,
					gb=g-b,
					br=b-r;
				counters[0]+=abs(r-W[0]);
				counters[1]+=abs(g-W[1]);
				counters[2]+=abs(b-W[2]);
				counters[3]+=abs(rg-W[3]);
				counters[4]+=abs(gb-W[4]);
				counters[5]+=abs(br-W[5]);
				W[0]=r;
				W[1]=g;
				W[2]=b;
				W[3]=rg;
				W[4]=gb;
				W[5]=br;
			}
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
		}

		low+=range*bestrct>>4;
		range=(range>>4)-1;
	}
	else
	{
		dstbuf=(unsigned char*)malloc(usize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		code=*(unsigned*)srcptr;
		srcptr+=4;
		code=code<<4|*(unsigned*)srcptr;
		srcptr+=4;

		bestrct=(int)((((code-low+1)<<4)-1)/range);
		low+=range*bestrct>>4;
		range=(range>>4)-1;
	}
	int psize=(iw+16LL)*(int)sizeof(short[4*4]);//4 padded rows * 4 channels max
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
	FILLMEM((unsigned short*)stats1, 0x8000, sizeof(stats1), sizeof(short));
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		char yuv[4]={0};
		int bit[3]={0};
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			short
				*NNN	=rows[3]+0*4,
				*NNNE	=rows[3]+1*4,
				*NN	=rows[2]+0*4,
				*NNE	=rows[2]+1*4,
				*NW	=rows[1]-1*4,
				*N	=rows[1]+0*4,
				*NE	=rows[1]+1*4,
				*WWWWW	=rows[0]-5*4,
				*WWWW	=rows[0]-4*4,
				*WWW	=rows[0]-3*4,
				*WW	=rows[0]-2*4,
				*W	=rows[0]-1*4,
				*curr	=rows[0]+0*4;
			unsigned short *stats0[]=
			{
				stats1[0][(4*(N[0]+W[0])+NE[0]-NW[0])>>3&255],
				stats1[1][(4*(N[1]+W[1])+NE[1]-NW[1])>>3&255],
				stats1[2][(4*(N[2]+W[2])+NE[2]-NW[2])>>3&255],
			};
			if(fwd)
			{
				yuv[0]=srcptr[idx+yidx]-128;
				yuv[1]=srcptr[idx+uidx]-128;
				yuv[2]=srcptr[idx+vidx]-128;

				yuv[2]-=yuv[vhelpidx];
				yuv[1]-=yuv[uhelpidx];
				
				curr[0]=yuv[0];
				curr[1]=yuv[1];
				curr[2]=yuv[2];
			}
			else
			{
				yuv[0]=0;
				yuv[1]=0;
				yuv[2]=0;
			}
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
			for(int kb=7, tidx0=1, tidx1=1, tidx2=1;kb>=0;--kb)
			{
				unsigned long long r2;
				unsigned short
					*cell0=stats0[0]+tidx0,
					*cell1=stats0[1]+tidx1,
					*cell2=stats0[2]+tidx2;
				unsigned short p0[]=
				{
					*cell0,
					*cell1,
					*cell2,
				};
				if(fwd)
				{
					bit[0]=yuv[0]>>kb&1;
					bit[1]=yuv[1]>>kb&1;
					bit[2]=yuv[2]>>kb&1;
						
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
					for(int kc=0;kc<3;++kc)
					{
						if(range<0x10000)
						{
							*(unsigned*)dstptr=(unsigned)(low>>32);
							dstptr+=4;
							low=low<<32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
#ifdef AC_VALIDATE
						unsigned long long lo0=low, r0=range;
#endif
						r2=range*p0[kc]>>16;
						if(bit[kc])
						{
							low+=r2;
							range-=r2;
						}
						else
							range=r2-1;
#ifdef AC_VALIDATE
						acval_enc(bit[kc], 0, p0[kc], lo0, lo0+r0, low, low+range, 0, 0);
#endif
					}
				}
				else
				{
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
					for(int kc=0;kc<3;++kc)
					{
						if(range<0x10000)
						{
							code=code<<32|*(unsigned*)srcptr;
							srcptr+=4;
							low=low<<32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
#ifdef AC_VALIDATE
						unsigned long long lo0=low, r0=range;
#endif
						r2=range*p0[kc]>>16;

						unsigned long long mid=low+r2;
						range-=r2;
						bit[kc]=code>=mid;
						if(bit[kc])
							low=mid;
						else
							range=r2-1;
						//bit[kc]=code>=low+r2;
						//if(bit[kc])
						//{
						//	low+=r2;
						//	range-=r2;
						//}
						//else
						//	range=r2-1;
#ifdef AC_VALIDATE
						acval_dec(bit[kc], 0, p0[kc], lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
					}
					yuv[0]|=bit[0]<<kb;
					yuv[1]|=bit[1]<<kb;
					yuv[2]|=bit[2]<<kb;
				}
				*cell0=p0[0]+(((!bit[0]<<16)-p0[0]+(1<<7>>1))>>7);
				*cell1=p0[1]+(((!bit[1]<<16)-p0[1]+(1<<7>>1))>>7);
				*cell2=p0[2]+(((!bit[2]<<16)-p0[2]+(1<<7>>1))>>7);
				tidx0=tidx0*2+bit[0];
				tidx1=tidx1*2+bit[1];
				tidx2=tidx2*2+bit[2];
			}
			if(!fwd)
			{
				curr[0]=(char)yuv[0];
				curr[1]=(char)yuv[1];
				curr[2]=(char)yuv[2];

				yuv[1]+=yuv[uhelpidx];
				yuv[2]+=yuv[vhelpidx];

				dstbuf[idx+yidx]=(unsigned char)(yuv[0]+128);
				dstbuf[idx+uidx]=(unsigned char)(yuv[1]+128);
				dstbuf[idx+vidx]=(unsigned char)(yuv[2]+128);

				guide_check(dstbuf, kx, ky);//
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
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
			*(unsigned*)dstptr=(unsigned)(low>>32);
			dstptr+=4;
			*(unsigned*)dstptr=(unsigned)low;
			dstptr+=4;

			fwrite("CH", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(dstbuf, 1, dstptr-dstbuf, fdst);
		}
		else
		{
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(dstbuf);
#if defined LOUD && !defined __GNUC__
	t=time_sec()-t;
	if(fwd)
	{
		ptrdiff_t csize=dstptr-dstbuf+10;
		printf("%12td/%12td  %12.6lf  %8.4lf%%\n", csize, usize, (double)usize/csize, (double)csize*100/usize);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
	return 0;
}