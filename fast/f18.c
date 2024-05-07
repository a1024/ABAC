#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

#include"ac.h"

#define FRAC_BITS 24
static unsigned exp2_fix24_neg(unsigned x)
{
	//return (unsigned)(exp2(-(x/16777216.))*0x1000000);//53% slower
	/*
	transcendental fractional powers of two
	x					inv(x)
	2^-0x0.000001 = 0x0.FFFFFF4F...		0x1.000000B1... = 2^0x0.000001
	2^-0x0.000002 = 0x0.FFFFFE9D...		0x1.00000163... = 2^0x0.000002
	2^-0x0.000004 = 0x0.FFFFFD3A...		0x1.000002C6... = 2^0x0.000004
	2^-0x0.000008 = 0x0.FFFFFA74...		0x1.0000058C... = 2^0x0.000008
	2^-0x0.000010 = 0x0.FFFFF4E9...		0x1.00000B17... = 2^0x0.000010
	2^-0x0.000020 = 0x0.FFFFE9D2...		0x1.0000162E... = 2^0x0.000020
	2^-0x0.000040 = 0x0.FFFFD3A3...		0x1.00002C5D... = 2^0x0.000040
	2^-0x0.000080 = 0x0.FFFFA747...		0x1.000058B9... = 2^0x0.000080
	2^-0x0.000100 = 0x0.FFFF4E8E...		0x1.0000B172... = 2^0x0.000100
	2^-0x0.000200 = 0x0.FFFE9D1D...		0x1.000162E5... = 2^0x0.000200
	2^-0x0.000400 = 0x0.FFFD3A3B...		0x1.0002C5CD... = 2^0x0.000400
	2^-0x0.000800 = 0x0.FFFA747F...		0x1.00058BA0... = 2^0x0.000800
	2^-0x0.001000 = 0x0.FFF4E91C...		0x1.000B175F... = 2^0x0.001000
	2^-0x0.002000 = 0x0.FFE9D2B3...		0x1.00162F39... = 2^0x0.002000
	2^-0x0.004000 = 0x0.FFD3A752...		0x1.002C605E... = 2^0x0.004000
	2^-0x0.008000 = 0x0.FFA75652...		0x1.0058C86E... = 2^0x0.008000
	2^-0x0.010000 = 0x0.FF4ECB59...		0x1.00B1AFA6... = 2^0x0.010000
	2^-0x0.020000 = 0x0.FE9E115C...		0x1.0163DAA0... = 2^0x0.020000
	2^-0x0.040000 = 0x0.FD3E0C0D...		0x1.02C9A3E7... = 2^0x0.040000
	2^-0x0.080000 = 0x0.FA83B2DB...		0x1.059B0D32... = 2^0x0.080000
	2^-0x0.100000 = 0x0.F5257D15...		0x1.0B5586D0... = 2^0x0.100000
	2^-0x0.200000 = 0x0.EAC0C6E8...		0x1.172B83C8... = 2^0x0.200000
	2^-0x0.400000 = 0x0.D744FCCB...		0x1.306FE0A3... = 2^0x0.400000
	2^-0x0.800000 = 0x0.B504F334...		0x1.6A09E667... = 2^0x0.800000
	*/
	static const unsigned long long frac_pots[]=
	{
		0x100000000,
		0x0FFFFFF4F,//extra 8 bits of precision
		0x0FFFFFE9D,
		0x0FFFFFD3A,
		0x0FFFFFA74,
		0x0FFFFF4E9,
		0x0FFFFE9D2,
		0x0FFFFD3A3,
		0x0FFFFA747,
		0x0FFFF4E8E,
		0x0FFFE9D1D,
		0x0FFFD3A3B,
		0x0FFFA747F,
		0x0FFF4E91C,
		0x0FFE9D2B3,
		0x0FFD3A752,
		0x0FFA75652,
		0x0FF4ECB59,
		0x0FE9E115C,
		0x0FD3E0C0D,
		0x0FA83B2DB,
		0x0F5257D15,
		0x0EAC0C6E8,
		0x0D744FCCB,
		0x0B504F334,
	};
#if 0
	unsigned long long x2=(unsigned long long)x<<1;
	unsigned long long r0=0x1000000, r1=0x1000000, r2=0x1000000, r3=0x1000000;
	for(int k=1;k<=FRAC_BITS;k+=8)
	{
		r0=r0*frac_pots[(k+0)&-(int)(x2>>(k+0)&1)]>>32;
		r1=r1*frac_pots[(k+1)&-(int)(x2>>(k+1)&1)]>>32;
		r2=r2*frac_pots[(k+2)&-(int)(x2>>(k+2)&1)]>>32;
		r3=r3*frac_pots[(k+3)&-(int)(x2>>(k+3)&1)]>>32;
		r0=r0*frac_pots[(k+4)&-(int)(x2>>(k+4)&1)]>>32;
		r1=r1*frac_pots[(k+5)&-(int)(x2>>(k+5)&1)]>>32;
		r2=r2*frac_pots[(k+6)&-(int)(x2>>(k+6)&1)]>>32;
		r3=r3*frac_pots[(k+7)&-(int)(x2>>(k+7)&1)]>>32;
	}
	r2=r2*r3>>24;
	r0=r0*r1>>24;
	r0=r0*r2>>24;
	r0>>=x>>FRAC_BITS;
	return (unsigned)r0;
#endif
#if 1
	unsigned long long x2=(unsigned long long)x<<1;
	unsigned long long r0=0x1000000, r1=0x1000000;
	for(int k=1;k<=FRAC_BITS;k+=2)
	{
		r0=r0*frac_pots[(k+0)&-(int)(x2>>(k+0)&1)]>>32;
		r1=r1*frac_pots[(k+1)&-(int)(x2>>(k+1)&1)]>>32;
		//r0=r0*frac_pots[(k+2)&-(int)(x2>>(k+2)&1)]>>32;
		//r1=r1*frac_pots[(k+3)&-(int)(x2>>(k+3)&1)]>>32;
		//r0=r0*frac_pots[(k+4)&-(int)(x2>>(k+4)&1)]>>32;
		//r1=r1*frac_pots[(k+5)&-(int)(x2>>(k+5)&1)]>>32;
		//r0=r0*frac_pots[(k+6)&-(int)(x2>>(k+6)&1)]>>32;
		//r1=r1*frac_pots[(k+7)&-(int)(x2>>(k+7)&1)]>>32;
	}
	r0*=r1;
	r0>>=(x>>FRAC_BITS)+24;
	return (unsigned)r0;
#endif
#if 0
	unsigned long long result=0x1000000;
	for(int k=0;k<FRAC_BITS;)//up to 24 muls
	{
		int bit=x>>k&1;
		++k;
		result=result*frac_pots[k&-bit]>>32;
	}
	result>>=x>>FRAC_BITS;
	return (unsigned)result;
#endif
}
static unsigned calc_cdf_continuous(int x, int depth, long long mad)
{
	int pos=x>0, posmask=-pos;
	long long xn=((long long)abs(x)<<(48-depth))/mad;
	if(xn>0x18000000)
		xn=0x18000000;
	int e=exp2_fix24_neg((unsigned)xn);
	e^=posmask;//negate then add 1 if x > mean
	e-=posmask;
	e>>=1;
	e+=pos<<24;
	//if(e<0)
	//	LOG_ERROR("");
	return (unsigned)e;
}
static unsigned long long calc_cdf_u(int sym, int depth, int half, int nlevels, long long mad, long long cdf_start, long long cdf_den)
{
	long long num=calc_cdf_continuous(sym-half, depth, mad)-cdf_start;
	num=CLAMP(-0xFFFFFFFFFFFF, num, 0xFFFFFFFFFFFF)*(0x10000LL-nlevels)/cdf_den+sym;
	return num;
}
#define CALC_CDF(X) calc_cdf_u(X, depths[kc], halfs[kc], nlevels[kc], mad, cdf_start, cdf_den)
int f18_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 0x10000, 0);
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	size_t bufsize=sizeof(short[4*4*2])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)_mm_malloc(bufsize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, bufsize);
	char depths[4]={0};
	memfill(depths, &image->depth, image->nch, sizeof(char));
	if(image->depth>=3)
	{
		depths[0]+=image->depth<16;
		depths[2]+=image->depth<16;
	}
	int nlevels[]=
	{
		1<<depths[0],
		1<<depths[1],
		1<<depths[2],
		1<<depths[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int perm[]={1, 2, 0, 3};
	for(int ky=0, idx=0, idx2=0;ky<image->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
		};
		for(int kx=0;kx<image->iw;++kx, idx+=3, ++idx2)
		{
			short
				*NNN	=rows[3]+0*8,
				*NNWW	=rows[2]-2*8,
				*NN	=rows[2]+0*8,
				*NNEE	=rows[2]+2*8,
				*NW	=rows[1]-1*8,
				*N	=rows[1]+0*8,
				*NE	=rows[1]+1*8,
				*WW	=rows[0]-2*8,
				*W	=rows[0]-1*8,
				*curr	=rows[0]+0*8;
			if(fwd)
			{
				memcpy(curr, image->data+idx, sizeof(short[3]));
				if(image->nch>=3)
				{
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
				}
			}
			for(int kc0=0;kc0<image->nch;++kc0)
			{
				int kc=perm[kc0];
				int vmin=MINVAR(N[kc], W[kc]), vmax=MAXVAR(N[kc], W[kc]);
				int pred=N[kc]+W[kc]-NW[kc], delta;
				pred=CLAMP(vmin, pred, vmax);
				long long mad, cdf_start, cdf_den, cdf, freq;
				
				//mad=(((long long)abs(N[kc+4])+(long long)abs(W[kc+4])+(long long)abs(NW[kc+4])+(long long)abs(NE[kc+4]))<<(24-2-depths[kc]))+1;
				//mad=(((long long)abs(N[kc+4])+(long long)abs(W[kc+4]))<<(24-1-depths[kc]))+1;
				//mad=(((long long)abs(NW[kc]-W[kc])+(long long)abs(N[kc]-NW[kc])+(long long)abs(NE[kc]-N[kc])+(long long)abs(W[kc]-NE[kc]))<<(24-3-depths[kc]))+1;
#if 1
				mad=1;
				int count=1;
				const int reach=5;
				for(int k=-reach;k<=reach;++k)
				{
					int weight=0x10000/(k*k+9);
					mad+=(long long)abs(NNN[(k<<3)+kc+4])*weight;
					count+=weight;
				}
				for(int k=-reach;k<=reach;++k)
				{
					int weight=0x10000/(k*k+4);
					mad+=(long long)abs(NN[(k<<3)+kc+4])*weight;
					count+=weight;
				}
				for(int k=-reach;k<=reach+1;++k)
				{
					int weight=0x10000/(k*k+1);
					mad+=(long long)abs(N[(k<<3)+kc+4])*weight;
					count+=weight;
				}
				for(int k=-reach-1;k<0;++k)
				{
					int weight=0x10000/(k*k+0);
					mad+=(long long)abs(curr[(k<<3)+kc+4])*weight;
					count+=weight;
				}
				mad<<=24-depths[kc];
				mad+=count-1LL;
				mad/=count;
#endif

				cdf_start=calc_cdf_continuous(-halfs[kc], depths[kc], mad);
				cdf_den=calc_cdf_continuous(halfs[kc], depths[kc], mad)-cdf_start;
				if(fwd)
				{
					delta=curr[kc]-pred;
					delta+=halfs[kc];
					delta&=nlevels[kc]-1;
					
					if(!cdf_den)
						ac_enc_bypass(&ec, delta, nlevels[kc]);
					else
					{
						cdf=CALC_CDF(delta);
						freq=CALC_CDF(delta+1)-cdf;

						ac_enc_update(&ec, (int)cdf, (int)freq);
					}
					delta-=halfs[kc];
					curr[kc+4]=delta;
				}
				else
				{
					unsigned code=ac_dec_getcdf(&ec);
					int range=nlevels[kc];
					delta=0;
					if(!cdf_den)
						delta=ac_dec_bypass(&ec, nlevels[kc]);
					else
					{
						while(range>1)
						{
							int floorhalf=range>>1;
						
							long long c2=CALC_CDF(delta+floorhalf);
							if(code>=c2)
								delta+=range-floorhalf;
							range=floorhalf;
						}
					
						cdf=CALC_CDF(delta);
						freq=CALC_CDF(delta+1)-cdf;
						ac_dec_update(&ec, (int)cdf, (int)freq);
					}
					
					curr[kc+4]=delta-halfs[kc];
					delta+=pred;
					delta&=nlevels[kc]-1;
					delta-=halfs[kc];
					curr[kc]=delta;
				}
			}
			if(!fwd)
			{
				memcpy(dst->data+idx, curr, sizeof(short[3]));
				if(image->nch>=3)
				{
					dst->data[idx+1]-=(dst->data[idx+0]+dst->data[idx+2])>>2;
					dst->data[idx+2]+=dst->data[idx+1];
					dst->data[idx+0]+=dst->data[idx+1];
				}
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	t0=time_sec()-t0;
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		if(fwd)
		{
			ptrdiff_t csize=list.nobj;
			printf("Memory usage:      %17.2lf KB\n", bufsize/1024.);

			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		}
		printf("%c       \t%16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
		if(fwd)
			printf("\n");
		//prof_print();
	}
	if(fwd)
		dlist_clear(&list);
	_mm_free(pixels);
	return 0;
}