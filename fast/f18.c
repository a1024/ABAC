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

static unsigned calc_cdf_continuous(int x, int depth, long long mad)
{
	int pos=x>0, posmask=-pos;
	long long xn=((long long)abs(x)<<(48-depth))/mad;
	if(xn>0x18000000)
		xn=0x18000000;
	int e=exp2_neg_fix24_avx2((unsigned)xn);
	e^=posmask;//negate then add 1 if x > mean
	e-=posmask;
	e>>=1;//signed shift
	e+=pos<<24;
	//if(e<0)
	//	LOG_ERROR("");
	return (unsigned)e;
}
static unsigned long long calc_cdf_u(unsigned sym, unsigned depth, unsigned half, unsigned nlevels, unsigned long long mad, unsigned long long cdf_start, unsigned long long cdf_den)
{
	unsigned long long num=(unsigned long long)calc_cdf_continuous(sym-half, depth, mad)-cdf_start;
	if(num>0xFFFFFFFFFFFF)
		num=0xFFFFFFFFFFFF;
	num=num*(0x10000LL-nlevels)/cdf_den+sym;
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