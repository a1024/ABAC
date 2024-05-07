#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;

#include"ac.h"

static long long calc_cdf_continuous(int x, long long mad)
{
	if(!x)//this should be hit frequently
		return 0x800000;
	long long pos=x>0, posmask=-pos;
	long long xn=((long long)abs(x)<<48)/mad;
	long long e=exp2_fix24(-(int)CLAMP(-0x7FFFFFFF, xn, 0x7FFFFFFF));
	e^=posmask;//negate then add 1 if x > mean
	e-=posmask;
	e>>=1;
	e+=pos<<24;
	return e;
}
static unsigned long long calc_cdf_u(int sym, int half, int nlevels, long long mad, long long cdf_start, long long cdf_den)
{
	long long num=calc_cdf_continuous(sym-half, mad)-cdf_start;
	num=CLAMP(-0xFFFFFFFFFFFF, num, 0xFFFFFFFFFFFF)*(0x10000LL-nlevels)/cdf_den+sym;
	return num;
}
#define CALC_CDF(X) calc_cdf_u(X, halfs[kc], nlevels[kc], mad, cdf_start, cdf_den)
//#define CALC_CDF(X) (((calc_cdf_continuous((X)-halfs[kc], mad)-cdf_start)*(0x10000-nlevels[kc]))/cdf_den+(X))
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
	//long long mad_clampers[]=
	//{
	//	0x1000000LL<<3>>(depths[0]>>1), 0x1000000LL<<3<<(depths[0]>>1),
	//	0x1000000LL<<3>>(depths[1]>>1), 0x1000000LL<<3<<(depths[1]>>1),
	//	0x1000000LL<<3>>(depths[2]>>1), 0x1000000LL<<3<<(depths[2]>>1),
	//	0x1000000LL<<3>>(depths[3]>>1), 0x1000000LL<<3<<(depths[3]>>1),
	//};
	int perm[]={1, 2, 0, 3};
	long long mads[]={1, 1, 1, 1};
	int counts[]={1, 1, 1, 1};
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

				mad=CLAMP(-0x7FFFFFFFFFFF, mads[kc], 0x7FFFFFFFFFFF);
				mad=(mad<<24)/counts[kc];
#if 0
				mad=1;//this offset decides the ratio between min & max values that the final MAD can take		([0~255]+offset)*59 = [1/256 ~ 255/256]
				int count=1;
				for(int k=-8;k<=8;++k)//(8<<1|1)*3+8 = 59
				{
					int a=abs(NNN[kc+(k<<3)+4]);
					count+=a!=0;
					mad+=(long long)a;
				}
				for(int k=-8;k<=8;++k)
				{
					int a=abs(NN[kc+(k<<3)+4]);
					count+=a!=0;
					mad+=(long long)a;
				}
				for(int k=-8;k<=8;++k)
				{
					int a=abs(N[kc+(k<<3)+4]);
					count+=a!=0;
					mad+=(long long)a;
				}
				for(int k=-8+1;k<=0;++k)
				{
					int a=abs(W[kc+(k<<3)+4]);
					count+=a!=0;
					mad+=(long long)a;
				}
				mad<<=24+depths[kc];
				mad+=count>>1;
				mad/=count;
#endif
				//mad=CLAMP(mad_clampers[kc<<1|0], mad, mad_clampers[kc<<1|1]);

				//if(ky==4&&kx==545)
				//if(ky==0&&kx==2&&kc==2)//
				//	printf("");

				cdf_start=calc_cdf_continuous(-halfs[kc], mad);
				cdf_den=calc_cdf_continuous(halfs[kc], mad)-cdf_start;
				if(fwd)
				{
					delta=curr[kc]-pred;
					delta+=halfs[kc];
					delta&=nlevels[kc]-1;
					
					if(!cdf_den)
						ac_enc_bypass(&ec, delta, nlevels[kc]);
					else
					{
						//if(ky==0&&kx==2&&kc==2)//
						//{
						//	for(int k=0;k<nlevels[kc];++k)
						//	{
						//		printf("%3d  %16lld\n", k, CALC_CDF(k));
						//	}
						//}
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
				mads[kc]+=(long long)abs(curr[kc+4]);
				++counts[kc];
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