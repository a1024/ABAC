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

static unsigned calc_cdf_continuous_v2(long long xn)
{
	int pos=xn>0, posmask=-pos;
	xn=llabs(xn);
	if(xn>0x18000000)
		xn=0x18000000;
	int e=exp2_neg_fix24_avx2((unsigned)xn);
	//int e=exp2_fix24_neg((unsigned)xn);
	e^=posmask;//negate then add 1 if x > mean
	e-=posmask;
	e>>=1;//signed shift
	e+=pos<<24;
	//if(e<0)
	//	LOG_ERROR("");
	return (unsigned)e;
}
static unsigned calc_cdf_v2(int isym, unsigned usym, unsigned long long invmad, unsigned cdf_start, unsigned long long scale)
{
	long long xn=(long long)isym*invmad;
	unsigned cdf=calc_cdf_continuous_v2(xn);
	if(cdf<=cdf_start)
		return usym;
	cdf-=cdf_start;
	cdf=(unsigned)((unsigned long long)cdf*scale>>16);
	cdf+=usym;
	return cdf;
}
#if 0
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
#endif
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
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
			pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
		};
		for(int kx=0;kx<image->iw;++kx, idx+=3)
		{
			short
				*NNN	=rows[3]+0*8,
			//	*NNWW	=rows[2]-2*8,
				*NN	=rows[2]+0*8,
			//	*NNEE	=rows[2]+2*8,
				*NW	=rows[1]-1*8,
				*N	=rows[1]+0*8,
			//	*NE	=rows[1]+1*8,
			//	*WW	=rows[0]-2*8,
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
				unsigned mad, cdf_start, cdf_den, cdf;
				int freq;
				
				//mad=(((long long)abs(N[kc+4])+(long long)abs(W[kc+4])+(long long)abs(NW[kc+4])+(long long)abs(NE[kc+4]))<<(24-2-depths[kc]))+1;
				//mad=(((long long)abs(N[kc+4])+(long long)abs(W[kc+4]))<<(24-1-depths[kc]))+1;
				//mad=(((long long)abs(NW[kc]-W[kc])+(long long)abs(N[kc]-NW[kc])+(long long)abs(NE[kc]-N[kc])+(long long)abs(W[kc]-NE[kc]))<<(24-3-depths[kc]))+1;
#if 1
				mad=1;
				//int count=0;
				const int reach=5;
				for(int k=-reach;k<=reach;++k)
				{
					//count+=(1<<16)/(k*k+3*3)+1;
					mad+=(abs(NNN[(k<<3)+kc+4])<<16)/(k*k+3*3)+1;
				}
				//{
				//	int weight=0x10000/(k*k+9);
				//	mad+=(abs(NNN[(k<<3)+kc+4])+1)*weight;
				//	//count+=weight;
				//}
				for(int k=-reach;k<=reach;++k)
				{
					//count+=(1<<16)/(k*k+2*2)+1;
					mad+=(abs(NN[(k<<3)+kc+4])<<16)/(k*k+2*2)+1;
				}
				//{
				//	int weight=0x10000/(k*k+4);
				//	mad+=(abs(NN[(k<<3)+kc+4])+1)*weight;
				//	//count+=weight;
				//}
				for(int k=-reach;k<=reach+1;++k)
				{
					//count+=(1<<16)/(k*k+1*1)+1;
					mad+=(abs(N[(k<<3)+kc+4])<<16)/(k*k+1*1)+1;
				}
				//{
				//	int weight=0x10000/(k*k+1);
				//	mad+=(abs(N[(k<<3)+kc+4])+1)*weight;
				//	//count+=weight;
				//}
				for(int k=-reach-1;k<0;++k)
				{
					//count+=(1<<16)/(k*k+0*0)+1;
					mad+=(abs(W[(k<<3)+kc+4])<<16)/(k*k+0*0)+1;
				}
				//{
				//	int weight=0x10000/(k*k+0);
				//	mad+=(abs(curr[(k<<3)+kc+4])+1)*weight;
				//	//count+=weight;
				//}
#endif
				//if(idx==0x0001D850)//
				//if(idx==0x000256a7)//
				//if(kx==image->iw/2&&ky==image->ih/2)//
				//if(idx==0x232AD)//
				//	printf("");
				unsigned long long invmad=(0x64046LL<<24)/mad;
				//unsigned long long invmad=((unsigned long long)count<<(48-depths[kc]))/mad;
				//unsigned long long invmad=((unsigned long long)count)/mad;
				//unsigned invmad=(unsigned)((0x6401F0000ULL<<depths[kc])/mad);
				unsigned long long scale;
				long long limit=(long long)halfs[kc]*invmad;
				cdf_start=calc_cdf_continuous_v2(-limit);
				cdf_den=(0x800000-cdf_start)<<1;
				if(cdf_den)
				{
					scale=((0x10000LL-nlevels[kc])<<16)/cdf_den;
					//scale=0xFFFFFFFF/cdf_den;

					//if(idx==0x232AD)//
					//{
					//	for(int ks=-halfs[kc];ks<halfs[kc];++ks)
					//	{
					//		unsigned cdf2=calc_cdf_v2(ks, ks+halfs[kc], invmad, cdf_start, scale);
					//		printf("%4d  %d\n", ks, cdf2);
					//	}
					//	LOG_ERROR("pause");
					//}

					if(fwd)
					{
						delta=curr[kc]-pred;
						delta+=halfs[kc];
						delta&=nlevels[kc]-1;
						int isym=delta-halfs[kc];
						cdf=calc_cdf_v2(isym, delta, invmad, cdf_start, scale);
						freq=calc_cdf_v2(isym+1, delta+1, invmad, cdf_start, scale)-cdf;
						//if(cdf>0x10000||freq<=0)
						//	LOG_ERROR("Invalid CDF");
						ac_enc_update(&ec, cdf, freq);
						curr[kc+4]=isym;
					}
					else
					{
						unsigned code=ac_dec_getcdf(&ec);
						int range=nlevels[kc];
						delta=0;
						while(range>1)
						{
							int floorhalf=range>>1;

							unsigned c2=delta+floorhalf;
							c2=calc_cdf_v2(c2-halfs[kc], c2, invmad, cdf_start, scale);
							if(code>=c2)
								delta+=range-floorhalf;
							range=floorhalf;
						}

						int isym=delta-halfs[kc];
						cdf=calc_cdf_v2(isym, delta, invmad, cdf_start, scale);
						freq=calc_cdf_v2(isym+1, delta+1, invmad, cdf_start, scale)-cdf;
						ac_dec_update(&ec, cdf, freq);
						curr[kc+4]=isym;
						delta+=pred;
						delta&=nlevels[kc]-1;
						delta-=halfs[kc];
						curr[kc]=delta;
					}
				}
				else
				{
					if(fwd)
					{
						delta=curr[kc]-pred;
						delta+=halfs[kc];
						delta&=nlevels[kc]-1;
						ac_enc_bypass(&ec, delta, nlevels[kc]);
						delta-=halfs[kc];
						curr[kc+4]=delta;
					}
					else
					{
						delta=ac_dec_bypass(&ec, nlevels[kc]);
						curr[kc+4]=delta-halfs[kc];
						delta+=pred;
						delta&=nlevels[kc]-1;
						delta-=halfs[kc];
						curr[kc]=delta;
					}
				}
#if 0
				mad<<=24-depths[kc];
				mad+=count-1LL;
				mad/=count;
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
					if(!cdf_den)
						delta=ac_dec_bypass(&ec, nlevels[kc]);
					else
					{
						unsigned code=ac_dec_getcdf(&ec);
						int range=nlevels[kc];
						delta=0;
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
#endif
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