#include"e2.h"
#define AC_IMPLEMENTATION
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;


//	#define FIXEDPREC
//	#define NAN_GUARD


#ifdef NAN_GUARD
float nan_guard(float x)
{
	if(!isfinite(x))
		LOG_ERROR("NaN");
	return x;
}
#else
#define nan_guard(...) __VA_ARGS__
#endif
//Xoroshiro128+ 1.0 by David Blackman and Sebastiano Vigna		https://prng.di.unimi.it/xoroshiro128plus.c
static unsigned long long xoroshiro128_state[2]={0xDF900294D8F554A5, 0x170865DF4B3201FC};
#define XOROSHIRO128_RESET() xoroshiro128_state[0]=0xDF900294D8F554A5, xoroshiro128_state[1]=0x170865DF4B3201FC
static inline unsigned long long rotl(const unsigned long long x, int k)
{
	return (x<<k)|(x>>(64-k));
}
unsigned long long xoroshiro128_next(void)
{
	const unsigned long long s0 = xoroshiro128_state[0];
	unsigned long long s1 = xoroshiro128_state[1];
	const unsigned long long result = s0 + s1;

	s1 ^= s0;
	xoroshiro128_state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	xoroshiro128_state[1] = rotl(s1, 37); // c

	return result;
}

#if 0
static inline void initialize(DataType *w, int count, DataType sqrt_fan_in)
{
	for(int k=0;k<count;++k)
	{
#ifdef FIXEDPREC
		w[k]=(DataType)((xoroshiro128_next()&0x1FFFF)-0x10000)/sqrt_fan_in;
#else
		int x=(int)(xoroshiro128_next()&0x1FFFF)-0x10000;
		w[k]=(DataType)x/(0x10000*sqrt_fan_in);
#endif
	}
	//int x=0;
	//for(int k=0;k<4;++k)
	//{
	//	x+=hamming_weight(xoroshiro128_next());
	//	x-=hamming_weight(xoroshiro128_next());
	//}
	//return (x<<8)*2/(fin+fout);
}
#endif

long long error_func_p32(long long x)
{
	//approximation 3 from Wikipedia
	const unsigned c[]=
	{
		0x120DCCEB,//0.0705230784
		0x0AD2FE74,//0.0422820123
		0x025F8DA3,//0.0092705272
		0x0009F660,//0.0001520143
		0x00122007,//0.0002765672
		0x0002D27E,//0.0000430638
	};
	int neg=x<0;
	unsigned long long x0=(unsigned long long)llabs(x), lo;//32.32 bits
	
	if(x0>0x48E7A0CE1)//erf(4.5565498399) = 0x0.FFFFFFFF800003 ~= 1 in 32.32bit
		lo=0x100000000;
	else
	{
		long long hi;
	#define MUL_P32(DST, A, B)	(lo=_mul128(A, B, &hi), DST=hi<<32|lo>>32)
		MUL_P32(lo, x0, c[5]), lo+=c[4];
		MUL_P32(lo, lo, x0), lo+=c[3];
		MUL_P32(lo, lo, x0), lo+=c[2];
		MUL_P32(lo, lo, x0), lo+=c[1];
		MUL_P32(lo, lo, x0), lo+=c[0];
		MUL_P32(lo, lo, x0), lo+=0x100000000;

		lo=_div128(1, lo>>1, lo, &hi);//round((0x100000000<<32)/lo) = (1<<64|lo>>1)/lo
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);

		lo=0x100000000-lo;
	}
	lo^=-neg;
	lo+=neg;
	lo>>=16;
	return (int)lo;
}
int error_func_p16(int x)
{
	//approximation 1
#if 0
	const int coeff[]=//1 - (1+(c[0]+(c[1]+(c[2]+c[3]*x)*x)*x)*x)^-4, max error=0.0005 = 0x20*2^-16
	{
		0x4745,//0.278393
		0x3AFB,//0.230389
		0x0040,//0.000972
		0x13FF,//0.078108
	};
	int x2=(int)((long long)coeff[3]*x>>16)+coeff[2];
	x2=(int)((long long)x2*x>>16)+coeff[1];
	x2=(int)((long long)x2*x>>16)+coeff[0];
	x2=(int)((long long)x2*x>>16)+0x10000;
	x2=0x10000/x2;
	x2=(int)((long long)x2*x2>>16);
	x2=(int)((long long)x2*x2>>16);
	x2=0x10000-x2;
	return x2;
#endif

	//approximation 3
#if 1
	const unsigned c[]=
	{
		0x120E,//0x120DCCEB,//0.0705230784
		0x0AD3,//0x0AD2FE74,//0.0422820123
		0x0260,//0x025F8DA3,//0.0092705272
		0x000A,//0x0009F660,//0.0001520143
		0x0012,//0x00122007,//0.0002765672
		0x0003,//0x0002D27E,//0.0000430638
	};
	int neg=x<0;
	unsigned long long x0=(unsigned long long)abs(x), res;//16.16 bits

	if(x0>0x32A1F)//erf(3.16453508) = 0x0.FFFF8000003 ~= 1 in 16.16 bit
		res=0x10000;
	else
	{
		res=(x0*c[5]>>16)+c[4];
		res=(res*x0>>16)+c[3];
		res=(res*x0>>16)+c[2];
		res=(res*x0>>16)+c[1];
		res=(res*x0>>16)+c[0];
		res=(res*x0>>16)+0x10000;

		res=(0x100000000|res>>1)/res;
		res=res*res>>16;
		res=res*res>>16;
		res=res*res>>16;
		res=res*res>>16;

		res=0x10000-res;
	}
	res^=-neg;
	res+=neg;
	return (int)res;
#endif
}

//	#define T16_ENABLE_CONF

//test 16: encode blocks using mixed histogram from causal neighbors (margin)
int get_mean_p16(const unsigned char *buf, int bw, int kc, int x1, int x2, int y1, int y2)
{
	long long sum=0;
	int res=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
			sum+=buf[(bw*ky+kx)<<2|kc];
	}
	int mean=(int)(((long long)sum<<16)/res);
	return mean;
}
int get_conf_p16(const unsigned char *buf, int bw, int kc, int x1, int x2, int y1, int y2, int mean, int *trivial)
{
	long long sum=0;
	int res=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			long long x=((long long)buf[(bw*ky+kx)<<2|kc]<<16)-mean;
			sum+=x*x>>16;
		}
	}
	if(!sum)
	{
		*trivial=1;
		return 0;
	}
	sum=(sum<<1)/res;//res is an integer, need variance*2 inside sqrt
	int conf=(int)(65536/sqrt(sum/65536.));//TODO don't use floating point
	return conf;
}
void print_rep(int count, char c)
{
	for(int k2=0;k2<count;++k2)
		printf("%c", c);
}
void print_CDF(const unsigned *CDF, const unsigned char *b2, int bw, int bh, int kc, int x1, int x2, int y1, int y2)
{
	const int graphwidth=32;
	int res=bw*bh;
#if 1
	int *CDF0=(int*)malloc(256*sizeof(int));
	if(!CDF0)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(CDF0, 0, 256*sizeof(int));
	int count=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=b2[(bw*ky+kx)<<2];
			++CDF0[sym];
		}
	}
	int sum=0;
	for(int sym=0;sym<256;++sym)
	{
		int freq=CDF0[sym];
		freq=(int)(((long long)freq*0xFF00)/count)+1;
		CDF0[sym]=sum;
		sum+=freq;
	}
#endif

	int mean=get_mean_p16(b2, bw, kc, x1, x2, y1, y2);
	//int mean=128<<16;
	int bypass=0;
	int conf=get_conf_p16(b2, bw, kc, x1, x2, y1, y2, mean, &bypass);

	//conf*=2.5;//
	//conf<<=1;//
	
	//printf("ACTUAL CDF,  STATIC ALPHA CDF,  CONF CDF\n");
	double entropy[3]={0};

	int start=(int)((long long)(0-mean)*conf>>16),
		end=(int)((long long)((256<<16)-mean)*conf>>16);
	int
		estart=error_func_p16(start),
		eend=error_func_p16(end);
	for(int k=0;k<256;++k)
	{
		int x=(int)((long long)((k<<16)-mean)*conf>>16);
		int e=error_func_p16(x);
		int den=eend-estart;
		int cdf2=(int)(((long long)(e-estart)*0xFF00)/den)+k;//guard

		int x2=(int)((long long)(((k+1)<<16)-mean)*conf>>16);
		int e2=error_func_p16(x2);
		int f2=(int)(((long long)(e-estart)*0xFF00)/den)+1;
		
#if 1
		printf("%3d", k);

		int c0=CDF0[k]*graphwidth>>16;
		printf(" %04X ", CDF0[k]);
		print_rep(c0, '*');
		print_rep(graphwidth-c0, ' ');

		int c1=CDF[k]*graphwidth>>16;
		printf(" %04X ", CDF[k]);
		print_rep(c1, '*');
		print_rep(graphwidth-c1, ' ');

		int c2=cdf2*graphwidth>>16;
		printf(" %04X ", cdf2);
		print_rep(c2, '*');
		print_rep(graphwidth-c2, ' ');

		printf("\n");
#endif

		double p0=((k<255?CDF0[k+1]:0x10000)-CDF0[k])/65536.;
		if(p0)
			entropy[0]-=p0*log2(p0);
		double p1=((k<255?CDF[k+1]:0x10000)-CDF[k])/65536.;
		if(p1)
			entropy[1]-=p1*log2(p1);
		double p2=f2/65536.;
		if(p2)
			entropy[2]-=p2*log2(p2);
	}
	printf("BPP %lf\t%lf\t%lf\n", entropy[0], entropy[1], entropy[2]);
	//printf("CR %lf\t%lf\t%lf\n", 8/entropy[0], 8/entropy[1], 8/entropy[2]);
	//printf("mean 0x%08X conf 0x%08X\n", mean, conf);
	//printf("Actual         %lf\n", 8/entropy[0]);
	//printf("Static Alpha   %lf\n", 8/entropy[1]);
	//printf("Conf           %lf\n", 8/entropy[2]);
	//printf("\n");
	free(CDF0);
}
void t16_prepblock(const unsigned char *b2, const unsigned short *CDF, int bw, int bh, int kc, int bx, int by, int alpha, int blockw, int blockh, int margin, unsigned *CDF2, int *xend, int *yend)
{
	int kx=bx*blockw, ky=by*blockh;

	*yend=ky+blockh<=bh?ky+blockh:bh;
	*xend=kx+blockw<=bw?kx+blockw:bw;
				
	int overflow=0;//CDF overflow can happen only once
	if(!bx&&!by)//first block has zero alpha
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF[sym+1];
			}
		}
		CDF2[256]=0x10000;
	}
	else
	{
		memset(CDF2, 0, 256*sizeof(unsigned));

		int count2=0;

		int left=kx-margin;
		if(left<0)
			left=0;
		int right=kx+blockw+margin;
		if(right>bw)
			right=bw;
		int top=ky-margin;
		if(top<0)
			top=0;

#ifdef T16_ENABLE_CONF
		double sdev=0;
#endif

		if(left<kx)//if left block is available
		{
			for(int ky2=ky;ky2<*yend;++ky2)
			{
				for(int kx2=left;kx2<kx;++kx2)//for each pixel
				{
					int sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx-kx2;
					if(dist<0||dist>margin)
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					//int inc=(kx2>>1)+1;
					//int inc=(kx2*3>>2)+1;
					//int inc=kx2<<1|1;
					//int inc=1;
					//int inc=kx2*kx2+1;
					//int inc=0x10000/(kx2+1);

					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif

					//if(count2<inc)
					//	LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**ysize;
			//count2+=blocksize**ysize;
		}
		if(top<ky)//if top block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=kx;kx2<*xend;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=ky-ky2;
					if(dist<0||dist>margin)
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					//int inc=(ky2>>1)+1;
					//int inc=(ky2*3>>2)+1;
					//int inc=ky2<<1|1;
					//int inc=1;
					//int inc=ky2*ky2+1;
					//int inc=0x10000/(ky2+1);
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif

					if(count2<inc)
						LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**xsize;
			//count2+=blocksize**xsize;
		}
		if(left<kx&&top<ky)//if topleft block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=left;kx2<kx;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx-kx2+ky-ky2;
					//int dist=MAXVAR(kx-kx2, ky-ky2);
					
					if(dist<0||dist>(margin<<1))
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif
				}
			}
		}
		if(right>kx+blockw&&top<ky)//if topright block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=kx+blockw;kx2<right;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx2-(kx+blockw)+ky-ky2;
					//int dist=MAXVAR(kx2-(kx+blockw), ky-ky2);

					if(dist<0||dist>(margin<<1))
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif
				}
			}
		}
#if 0
		int xoffset, yoffset;
		int xsize2, ysize2, count2;
		if(!bx)
		{
			if(!by)
				xoffset=0, yoffset=0;
			else
				xoffset=0, yoffset=blocksize;
		}
		else
			xoffset=blocksize, yoffset=0;
		ysize2=ky-yoffset+blocksize<=bh?blocksize:bh-(ky-yoffset);
		xsize2=kx-xoffset+blocksize<=bw?blocksize:bw-(kx-xoffset);
		count2=xsize2*ysize2;

		for(int ky2=0;ky2<ysize2;++ky2)
		{
			for(int kx2=0;kx2<xsize2;++kx2)//for each pixel
			{
				unsigned char sym=b2[(bw*(ky-yoffset+ky2)+kx-xoffset+kx2)<<2|kc];
				++h2[sym];
			}
		}
#endif
#ifdef T16_ENABLE_CONF
		sdev/=count2;
		sdev=sqrt(sdev);
		double conf=1/(sqrt(2)*sdev);
#endif

		int sum=0;
		for(int sym=0;sym<256;++sym)
		{
			int cdf1=!overflow?CDF[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF[sym+1];
			int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

			int f2=(int)(((long long)CDF2[sym]<<16)/count2);//normalize
			
#ifdef T16_ENABLE_CONF
			double start=erf((0-128)*conf), end=erf((256-128)*conf), x=erf((sym-128)*conf), x2=erf((sym+1-128)*conf);
			int f3=(int)((x2-x)*0x10000/(end-start));
#endif

			//if(f1||f2)//
			//	printf("");

			//int freq=f1+(int)(((long long)f3-f1)*alpha>>16);//blend
			//freq=freq+(int)(((long long)f2-freq)*alpha>>16);//blend

			//int freq=f2;
			int freq=f1+(int)(((long long)f2-f1)*alpha>>16);//blend

#ifdef T16_ENABLE_CONF
			freq=freq+(int)(((long long)f3-freq)*alpha>>16);//blend2
#endif

			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
			if(freq<0||freq>0xFF01)
				LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
				LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
		}
		CDF2[256]=0x10000;
	}
}
size_t test16_encode(const unsigned char *src, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, ArrayHandle *data, int loud, int *csizes)
{
	double t_start=time_sec();
	int res=bw*bh;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);
	apply_transforms_fwd(b2, bw, bh);

	int totalhsize=256*3;
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDF=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!hist||!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	for(int kc=0;kc<3;++kc)
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=b2[k<<2|kc];
			++hist[kc<<8|sym];
		}
	}
	t25_normalize_histogram(hist, 256, res, CDF);//this is just to pack the histogram, CDF is renormalized again with ramp guard
	t25_normalize_histogram(hist+256, 256, res, CDF+256);
	t25_normalize_histogram(hist+512, 256, res, CDF+512);

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, CDF, 768*sizeof(short));

	for(int kc=0;kc<3;++kc)
	{
		int bxcount=(bw+blockw[kc]-1)/blockw[kc],
			bycount=(bh+blockh[kc]-1)/blockh[kc];

		unsigned state=0x10000;
		for(int by=bycount-1;by>=0;--by)
		{
			int ky=by*blockh[kc];
			for(int bx=bxcount-1;bx>=0;--bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;

				int kx=bx*blockw[kc];
				int xend=0, yend=0;
				t16_prepblock(b2, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], margin[kc], h2, &xend, &yend);

				//if(kc==0&&bx==0&&by==0)
				//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

				//encode block
				for(int ky2=yend-1;ky2>=ky;--ky2)
				{
					for(int kx2=xend-1;kx2>=kx;--kx2)//for each pixel
					{
						unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];

						int cdf=h2[sym], freq=h2[sym+1]-cdf;
#if 0
						int cdf=CDF[kc<<8|sym], freq=(sym<255&&cdf<CDF[kc<<8|(sym+1)]?CDF[kc<<8|(sym+1)]:0x10000)-cdf, cdf2, f2;
						if(bx||by)
						{
							cdf2=h2[sym], f2=h2[sym+1]-cdf2;
							cdf2=(cdf2<<16)/h2[256];
							f2=(f2<<16)/h2[256];
							cdf=cdf+(int)(((long long)cdf2-cdf)*alpha>>16);
							freq=freq+(int)(((long long)f2-freq)*alpha>>16);
							cdf=(cdf*0xFF00>>16)+sym;
							freq=(freq*0xFF00>>16)+1;
						}
#endif

						//if(kc==2&&ky2==0&&kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);

						if(!freq)
							LOG_ERROR("ZPS");
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						state=state/freq<<16|(cdf+state%freq);//update
					}
				}
			}
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	int overhead=12+(int)(totalhsize*sizeof(short));
	int ch[]=
	{
		ansbookmarks[0]-overhead,
		ansbookmarks[1]-ansbookmarks[0],
		ansbookmarks[2]-ansbookmarks[1],
	};
	if(csizes)
	{
		csizes[0]=ch[0];
		csizes[1]=ch[1];
		csizes[2]=ch[2];
	}
	if(loud)
	{
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("alpha 0x%04X\n", alpha);
		printf("Total    %7d  %lf\n", ansbookmarks[2], 3.*res/ansbookmarks[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf  %2dx%d  M %d\n", ch[0], (double)res/ch[0], blockw[0], blockh[0], margin[0]);
		printf("Green    %7d  %lf  %2dx%d  M %d\n", ch[1], (double)res/ch[1], blockw[1], blockh[1], margin[1]);
		printf("Blue     %7d  %lf  %2dx%d  M %d\n", ch[2], (double)res/ch[2], blockw[2], blockh[2], margin[2]);
	}

	dlist_clear(&list);
	free(b2);
	free(hist);
	free(CDF);
	free(h2);
	return 1;
}
int threeway_uint32(const void *p1, const void *p2)
{
	unsigned v1=*(const unsigned*)p1, v2=*(const unsigned*)p2;
	return (v1>v2)-(v1<v2);
}
unsigned char *debug_ptr=0;
int test16_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, unsigned char *buf)
{
	const int cdflen=768LL*sizeof(short), overhead=12LL+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	
	unsigned short *CDF=(unsigned short*)malloc(cdflen);
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(CDF, data+12, cdflen);
	
	for(int kc=0;kc<3;++kc)
	{
		int bxcount=(bw+blockw[kc]-1)/blockw[kc],
			bycount=(bh+blockh[kc]-1)/blockh[kc];

		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		for(int by=0;by<bycount;++by)
		{
			int ky=by*blockh[kc];
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;

				int kx=bx*blockw[kc];
				int xend=0, yend=0;
				t16_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], margin[kc], h2, &xend, &yend);
				for(int ky2=ky;ky2<yend;++ky2)
				{
					for(int kx2=kx;kx2<xend;++kx2)//for each pixel
					{
						unsigned c=(unsigned short)state;
						size_t sym=0;
						int found=binary_search(h2, 257, sizeof(int), threeway_uint32, &c, &sym);
						sym-=!found;//binary_search gives insertion index

						//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*ky2+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf=h2[sym], freq=h2[sym+1]-cdf;

						state=freq*(state>>16)+c-cdf;//update
						if(state<0x10000)//renorm
						{
							state<<=16;
							if(srcptr-2>=srcstart)
							{
								srcptr-=2;
								memcpy(&state, srcptr, 2);
							}
						}
					}
				}
			}
		}
	}
	free(CDF);
	free(h2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	apply_transforms_inv(buf, bw, bh);
	return 1;
}

//	#define T25_OPT_PRED

//test 25 (the 4th optimizer for test16): image is divided into large blocks, each lblock has an sblock selected (least compressible part) on which the T16 params are optimized
typedef struct T25ParamsStruct
{
	short
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} T25Params;
typedef struct T25ParamsPackedStruct
{
	unsigned char
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} T25ParamsPacked;
typedef struct RectStruct
{
	int x1, x2, y1, y2;
} Rect;
#define T25_PARAM(P, IDX) ((short*)(P))[IDX]
static T25Params t25_limits={32, 40, 40, 40, 255, 96};
static int t25_ctr=0;
int t25_incparam(T25Params *param, int pidx, int step)
{
	int prevval=0;
	switch(pidx)
	{
	case 0:prevval=param->gwidth, param->gwidth+=step; if(param->gwidth<8)param->gwidth=8; break;
	case 1:prevval=param->mleft , param->mleft +=step; if(param->mleft <0)param->mleft =0; break;
	case 2:prevval=param->mtop  , param->mtop  +=step; if(param->mtop  <0)param->mtop  =0; break;
	case 3:prevval=param->mright, param->mright+=step; if(param->mright<0)param->mright=0; break;
	case 4:prevval=param->alpha , param->alpha +=step; if(param->alpha <0)param->alpha =0; else if(param->alpha>0xFF)param->alpha=0xFF; break;
	case 5:prevval=param->maxinc, param->maxinc+=step; if(param->maxinc<1)param->maxinc=1; break;
	}
	if(T25_PARAM(param, pidx)>T25_PARAM(&t25_limits, pidx))
		T25_PARAM(param, pidx)=T25_PARAM(&t25_limits, pidx);

	if(!BETWEEN_INC(1, param->gwidth, t25_limits.gwidth)||!BETWEEN_INC(0, param->mleft, t25_limits.mleft)||!BETWEEN_INC(0, param->mtop, t25_limits.mtop)||!BETWEEN_INC(0, param->mright, t25_limits.mright)||!BETWEEN_INC(0, param->alpha, t25_limits.alpha)||!BETWEEN_INC(0, param->maxinc, t25_limits.maxinc))
		LOG_ERROR("Invalid params INC  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	
	return prevval;
}
void t25_normalize_histogram(const unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	if(!nsymbols)//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
		return;
	}
	unsigned sum=0, qfreq;
	for(int sym=0;sym<nlevels;++sym)
	{
		qfreq=((long long)srchist[sym]<<16)/nsymbols;
		CDF[sym]=sum;
		sum+=qfreq;
	}
}
void t25_addhist(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0a, int x0b, int y0, int maxinc, unsigned *CDF2)
{
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(iw*ky+kx)<<2|kc];
			int dist=abs(ky-y0);
			if(kx<x0a)
				dist+=abs(kx-x0a);
			else if(kx>x0b)
				dist+=abs(kx-x0b);
			int inc=maxinc-dist;
			if(inc>0)
			{
				CDF2[sym]+=inc;
				CDF2[256]+=inc;
			}
		}
	}
}
int t25_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int iw, int ih, int kc, int x1, int x2, int y, T25Params const *p, unsigned *CDF2, int rec)
{
	int overflow=0;
	int sum, cdf1, f1, f2, freq;
	memset(CDF2, 0, 257*sizeof(unsigned));
	if(p->mtop)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x2+p->mright, y-p->mtop, y, x1, x2, y, p->maxinc, CDF2);
	if(p->mleft)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x1, y, y+1, x1, x2, y, p->maxinc, CDF2);

	if(CDF2[256])
	{
		sum=0;
		if(rec)//
			printf("\nCXXY %d %d %d %d A %d\n", kc, x1, x2, y, p->alpha);
		for(int sym=0;sym<256;++sym)
		{
			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize

			if(f2<f1)
				freq=f2+(int)(((long long)f1-f2)*(0xFF-p->alpha)/0xFF);//blend
			else
				freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

			int f3=freq;//

			freq=(int)((long long)freq*0xFF00>>16)+1;//guard
			//freq=CLAMP(0, freq, 0xFF01);

			if(rec)//
				printf("%3d 0x%04X 0x%04X 0x%04X 0x%04X\n", sym, f1, f2, f3, freq);

			if(freq<0||freq>0xFF01)
			{
				printf("Impossible freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
				return 0;
			}
				//LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000&&sym<255)
			{
				if(!rec)//
				{
					t25_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, y, p, CDF2, 1);//
					//for(int k=0;k<=sym;++k)
					//	printf("%3d 0x%04X 0x%04X\n", k, CDF0[k], CDF2[k]);
					printf("ANS CDF sym 0x%02X sum 0x%04X  freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", sym, sum, freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
					//printf("ANS CDF sym 0x%02X sum 0x%04X freq 0x%04X\n", sym, sum, freq);
				}
				return 0;
			}
		}
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
	}
	CDF2[256]=0x10000;
	return 1;
}
double t25_calcloss(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params const *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double chsize=0;
	if(!BETWEEN_INC(1, param->gwidth, t25_limits.gwidth)||!BETWEEN_INC(0, param->mleft, t25_limits.mleft)||!BETWEEN_INC(0, param->mtop, t25_limits.mtop)||!BETWEEN_INC(0, param->mright, t25_limits.mright)||!BETWEEN_INC(0, param->alpha, t25_limits.alpha)||!BETWEEN_INC(0, param->maxinc, t25_limits.maxinc))
		LOG_ERROR("Invalid params LOSS W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	//if(loud)
	//	printf("W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\r",
	//		param->gwidth,
	//		param->mleft,
	//		param->mtop,
	//		param->mright,
	//		param->alpha,
	//		param->maxinc);
	for(int ky=r->y1;ky<r->y2;++ky)
	{
		for(int kx=r->x1;kx<r->x2;)
		{
			int xend=MINVAR(kx+param->gwidth, r->x2);
			int success=t25_prepblock(buf, CDF0, iw, ih, kc, kx, xend, ky, param, CDF2, 0);
			if(!success)
				return 0;
				
			int kx2=kx;
			for(;kx2<kx+param->gwidth;++kx2)
			{
				unsigned char sym=buf[(iw*ky+kx2)<<2|kc];
				int freq=CDF2[sym+1]-CDF2[sym];
				double prob=(double)freq/0x10000, bitsize=-log2(prob);//Zipf's law
				chsize+=bitsize;
			}
			kx+=param->gwidth;
		}
	}
	++t25_ctr;
	return chsize;
}
#if 0
int t25_opt2(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, const unsigned short *CDF0, unsigned *CDF2, T25Params *param, double *csize, int pidx, int step, int loud)
{
	int went_fwd=0;
	for(int subit=0;subit<20;++subit)
	{
		double csize0=*csize;
		short prevval=t25_incparam(param, pidx, step);
		if(prevval==T25_PARAM(param, pidx))//out of range
			break;
		if(loud)
			printf("W%c%3d  MLTR %c%3d %c%3d %c%3d  A%c0x%02X I%c%3d\r",
				pidx==0?'>':' ', param->gwidth,
				pidx==1?'>':' ', param->mleft,
				pidx==2?'>':' ', param->mtop,
				pidx==3?'>':' ', param->mright,
				pidx==4?'>':' ', param->alpha,
				pidx==5?'>':' ', param->maxinc);
		*csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2);
		if(*csize>csize0)//cancel last change and break
		{
			T25_PARAM(param, pidx)=prevval;
			*csize=csize0;
			break;
		}
		went_fwd=1;
	}
	return went_fwd;
}
double t25_optimize(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double csize00, csize;
	int prevval0;
	int steps[]={16, 8, 4, 2, 1};

	for(int ks=0;ks<_countof(steps);++ks)
	{
		int step=steps[ks];
		for(int it=0, improve=1;it<64&&improve;++it)
		{
			improve=0;
			for(int pidx=0;pidx<sizeof(T25Params)/sizeof(short);++pidx)
			{
				prevval0=T25_PARAM(param, pidx);
				csize00=csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2);

				int went_fwd=t25_opt2(buf, iw, ih, kc, r, CDF0, CDF2, param, &csize, pidx, step, loud);
				if(!went_fwd)
					t25_opt2(buf, iw, ih, kc, r, CDF0, CDF2, param, &csize, pidx, -step, loud);

				if(csize>csize00)//prevent CR from worsening
				{
					T25_PARAM(param, pidx)=prevval0;
					csize=csize00;
				}
			}
		}
	}
	return csize;
}
#endif
double t25_optimize_v2(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double csize0;
	int steps[]={16, 8, 4, 2, 1};
	//	limit[]={ 1,  1, 1, 1, 1, 1};
	//int steps[]={64, 16, 4, 1},
	//	limit[]={4, 4, 4, 4};
	
	csize0=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
	if(!csize0)
		LOG_ERROR("Start W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	for(int ks=0;ks<_countof(steps);++ks)
	{
		int step=steps[ks];
		int bestpidx=0, beststep=0;
		double bestcsize=csize0;
			
		for(int pidx=0;pidx<sizeof(T25Params)/sizeof(short);++pidx)
		{
			double csize;
			short prevval;

			prevval=t25_incparam(param, pidx, step);
			if(T25_PARAM(param, pidx)!=prevval)
			{
				csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
				if(!csize0)
					LOG_ERROR("Plus W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
				T25_PARAM(param, pidx)=prevval;
				if(bestcsize>csize)
					bestcsize=csize, bestpidx=pidx, beststep=step;
			}

			prevval=t25_incparam(param, pidx, -step);
			if(T25_PARAM(param, pidx)!=prevval)
			{
				csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
				if(!csize0)
					LOG_ERROR("Minus W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
				T25_PARAM(param, pidx)=prevval;
				if(bestcsize>csize)
					bestcsize=csize, bestpidx=pidx, beststep=-step;
			}
		}
		if(bestcsize<csize0)
		{
			t25_incparam(param, bestpidx, beststep);
			csize0=bestcsize;
		}
	}
	return csize0;
}
#if 0
int t25_optimizeall(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int loud)
{
	int res=iw*ih;
	unsigned short *CDF0=(unsigned short*)malloc(256LL*3*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++CDF2[sym];
		}
		t25_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}
	
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	
	double csizes[3]={0};
	int steps[]={4, 2, 1};
	int usize=(x2-x1)*(y2-y1);
	for(int kc=0;kc<3;++kc)
	{
		for(int ks=0;ks<_countof(steps);++ks)
		{
			for(int it=0, improve=1;it<64&&improve;++it)
			{
				improve=0;
				for(int pidx=0;pidx<sizeof(t25_params)/sizeof(t25_params->gwidth);++pidx)
				{
					csizes[kc]=t25_optimize(buf, iw, ih, kc, x1, x2, y1, y2, t25_params+kc, pidx, steps[ks], CDF0+((size_t)kc<<8), CDF2);
					if(loud)
						io_render();
				}
			}
		}
		csizes[kc]/=8;
		t25_cr[kc]=usize>0&&csizes[kc]?(x2-x1)*(y2-y1)/csizes[kc]:0;
	}

	free(CDF0);
	free(CDF2);
	return 1;
}
#endif

void t25_calchist(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, unsigned *hist)
{
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=iw*ky+kx;
			unsigned char sym=buf[idx<<2|kc];
			++hist[sym];
		}
	}
}
double t25_calccsize(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, unsigned *hist)
{
	memset(hist, 0, 256*sizeof(int));
	t25_calchist(buf, iw, ih, kc, r->x1, r->x2, r->y1, r->y2, hist);
	//for(int ky=r->y1;ky<r->y2;++ky)
	//{
	//	for(int kx=r->x1;kx<r->x2;++kx)
	//	{
	//		int idx=iw*ky+kx;
	//		unsigned char sym=buf[idx<<2|kc];
	//		++hist[sym];
	//	}
	//}
	int count=(r->x2-r->x1)*(r->y2-r->y1);
	double bitsize=0;
	for(int ky=r->y1;ky<r->y2;++ky)
	{
		for(int kx=r->x1;kx<r->x2;++kx)
		{
			int idx=iw*ky+kx;
			unsigned char sym=buf[idx<<2|kc];
			unsigned freq=hist[sym];
			double prob=(double)freq/count;
			bitsize-=log2(prob);//Zipf's law
		}
	}
	double csize=bitsize/8;
	return csize;
}
void t25_selectsmallblock(const unsigned char *buf, int iw, int ih, int kc, Rect const *lb, Rect *sb, int sbw, int sbh, unsigned *hist)
{
	int bw=lb->x2-lb->x1,
		bh=lb->y2-lb->y1,
		nbx=bw/sbw, nby=bh/sbh;
	if(nbx<=1||nby<=1)
	{
		*sb=*lb;
		return;
	}
	//memset(hist, 0, 256*sizeof(int));
	//for(int ky=lb->y1;ky<lb->y2;++ky)
	//{
	//	for(int kx=lb->x1;kx<lb->x2;++kx)
	//	{
	//		int idx=iw*ky+kx;
	//		unsigned char sym=buf[idx<<2|kc];
	//		++hist[sym];
	//	}
	//}
	double bestcsize=0;
	int bestx=0, besty=0;
	int it=0;
	for(int by=0;by<nby;++by)
	{
		int ky=lb->y1+by*sbh;
		for(int bx=0;bx<nbx;++bx, ++it)
		{
			int kx=lb->x1+bx*sbw;
			Rect r_cand={kx, kx+sbw, ky, ky+sbh};
			double csize=t25_calccsize(buf, iw, ih, kc, &r_cand, hist);
			if(!it||bestcsize>csize)
				bestcsize=csize, bestx=kx, besty=ky;
		}
	}
	sb->x1=bestx;
	sb->x2=bestx+sbw;
	sb->y1=besty;
	sb->y2=besty+sbh;
}

static T25Params t25_params[3]=
{
	{ 8, 26, 26, 26, 0xD4, 32},
	{23, 32, 32, 32, 0xD4, 32},
	{ 8, 26, 26, 26, 0xD4, 32},
};
int t25_encode(const unsigned char *src, int iw, int ih, int *blockw, int *blockh, int use_ans, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifndef T25_OPT_PRED
	apply_transforms_fwd(buf2, iw, ih);
#else
	addbuf(buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd((char*)buf2, iw, ih);
	pred_opt_opt_v6((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v5((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v4((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v3((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v2((char*)buf2, iw, ih);
#endif
	short predparams[22+PW2_NPARAM];
	const short predidx[]={0, 11, 11+PW2_NPARAM}, predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM};
	memcpy(predparams+predidx[0], jxlparams_i16,         predlen[0]*sizeof(short));
	memcpy(predparams+predidx[1], pw2_params+PW2_NPARAM, predlen[1]*sizeof(short));
	memcpy(predparams+predidx[2], jxlparams_i16+22,      predlen[2]*sizeof(short));
#ifdef T25_OPT_PRED
	//pred_jxl_opt_v2((char*)buf2, iw, ih, jxlparams_i16, loud);
#if 1
	if(loud)
		pred_opt_printparam();
	//{
	//	for(int kc=0;kc<3;++kc)//
	//	{
	//		for(int kp=0;kp<predlen[kc];++kp)
	//		{
	//			short val=predparams[predidx[kc]+kp];
	//			printf(" %c0x%04X,", val<0?'-':' ', abs(val));
	//		}
	//		//for(int kp=0;kp<11;++kp)
	//		//{
	//		//	short val=jxlparams_i16[11*kc+kp];
	//		//	printf(" %c0x%04X,", val<0?'-':' ', abs(val));
	//		//}
	//		printf("\n");
	//	}
	//}
#endif
	pred_opt_apply((char*)buf2, iw, ih, 1);
	//pred_jxl_apply((char*)buf2, iw, ih, jxlparams_i16, 1);
	addbuf(buf2, iw, ih, 3, 4, 128);
#endif

	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, CDF2);
		t25_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, predparams, (predlen[3])*sizeof(short));
	//dlist_push_back(&list, jxlparams_i16, 33*sizeof(short));
	dlist_push_back(&list, CDF0, 768*sizeof(short));

	int overhead[4]={0, 0, 0, (int)list.nobj};//

	for(int kc=0;kc<3;++kc)
	{
		int lbw=blockw[kc], lbh=blockh[kc],
		//	sbw=sbsizes[kc<<1], sbh=sbsizes[kc<<1|1],
			lbx=(iw+lbw-1)/lbw, lby=(ih+lbh-1)/lbh;
		Rect lblock;
		//Rect sblock;
		ArrayHandle params;
		ARRAY_ALLOC(T25ParamsPacked, params, 0, 0, (size_t)lbx*lby, 0);
		for(int by=0;by<lby;++by)
		{
			lblock.y1=by*lbh;
			lblock.y2=MINVAR(lblock.y1+lbh, ih);
			for(int bx=0;bx<lbx;++bx)
			{
				lblock.x1=bx*lbw;
				lblock.x2=MINVAR(lblock.x1+lbw, iw);

				t25_ctr=0;

#if 0
				lblock.x1=640;//
				lblock.x2=lblock.x1+lbw;
				lblock.y1=384;
				lblock.y2=lblock.y1+lbh;
				t25_params[kc].gwidth=1;
				t25_params[kc].mleft=30;
				t25_params[kc].mtop=6;
				t25_params[kc].mright=32;
				t25_params[kc].alpha=0xD7;
				t25_params[kc].maxinc=32;
#endif

				t25_optimize_v2(buf2, iw, ih, kc, &lblock, t25_params+kc, CDF0, CDF2, loud);//optimize for whole large block

				//t25_selectsmallblock(buf2, iw, ih, kc, &lblock, &sblock, sbw, sbh, CDF2);
				//t25_optimize(buf2, iw, ih, kc, &sblock, t25_params+kc, CDF0, CDF2);//optimize for small block

				T25ParamsPacked *pp=(T25ParamsPacked*)ARRAY_APPEND(params, 0, 1, 1, 0);
				for(int k=0;k<sizeof(T25ParamsPacked);++k)
					((unsigned char*)pp)[k]=(unsigned char)T25_PARAM(t25_params+kc, k);
				//pp->gwidth=(unsigned char)t25_params[kc].gwidth;
				//pp->mleft =(unsigned char)t25_params[kc].mleft;
				//pp->mtop  =(unsigned char)t25_params[kc].mtop;
				//pp->mright=(unsigned char)t25_params[kc].mright;
				//pp->alpha =(unsigned char)t25_params[kc].alpha;
				//pp->maxinc=(unsigned char)t25_params[kc].maxinc;

				if(loud)
					printf("[%d/3 %2d/%2d %2d/%2d] CXY %d %4d %4d  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d  %d iters\r", kc+1, by+1, lby, bx+1, lbx, kc, lblock.x1, lblock.y1, t25_params[kc].gwidth, t25_params[kc].mleft, t25_params[kc].mtop, t25_params[kc].mright, t25_params[kc].alpha, t25_params[kc].maxinc, t25_ctr);
			}
		}
		if(loud)
			printf("\n");
		dlist_push_back(&list, params->data, params->count*params->esize);

		overhead[kc]=(int)(params->count*params->esize);//

		//if(kc==2)//
		//	kc=2;

		if(use_ans)
		{
			ANSCoder ec;
			ans_enc_init(&ec, &list);
			//unsigned state=0x10000;
			for(int ky=ih-1;ky>=0;--ky)
			{
				int by=ky/lbh;
				lblock.y1=by*lbh;
				lblock.y2=MINVAR(lblock.y1+lbh, ih);
				for(int bx=lbx-1;bx>=0;--bx)
				{
					T25ParamsPacked *pp=(T25ParamsPacked*)array_at(&params, lbx*by+bx);
					T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
					int ng=(lbw+param.gwidth-1)/param.gwidth;
					lblock.x1=bx*lbw;
					lblock.x2=MINVAR(lblock.x1+lbw, iw);
					for(int kg=ng-1;kg>=0;--kg)
					{
						int x1=lblock.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, lblock.x2);
						int success=t25_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2, 0);
						if(!success)
							LOG_ERROR("t25_prepblock error");
						for(int kx=x2-1;kx>=x1;--kx)
						{
							ans_enc(&ec, buf2[(iw*ky+kx)<<2|kc], CDF2, 256);
							//unsigned char sym=buf2[(iw*ky+kx)<<2|kc];
							//ans_enc(&state, CDF2[sym], CDF2[sym+1]-CDF2[sym], &list);
#if 0
							int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

							//if(kc==0&&ky==511&&kx==512)//
							//	printf("CXY %d %d %d  sym 0x%02X cdf 0x%04X freq 0x%04X  state 0x%08X\n", kc, kx, ky, sym, cdf, freq, state);

							if(!freq)
								LOG_ERROR("ZPS");

							//double prob=freq/65536.;
							//csize-=log2(prob);
							//unsigned s0=state;
						
							if(state>=(unsigned)(freq<<16))//renorm
							{
								dlist_push_back(&list, &state, 2);
								state>>=16;
							}
							debug_enc_update(state, cdf, freq, kx, ky, 0, kc, sym);
							state=state/freq<<16|(cdf+state%freq);//update
#endif
						}
					}
				}
			}
			ans_enc_flush(&ec);
			//ans_enc_flush(&list, state);
			//dlist_push_back(&list, &state, 4);
		}
		else
		{
			ArithmeticCoder ec;
			ac_enc_init(&ec, &list);
			for(int ky=0;ky<ih;++ky)
			{
				int by=ky/lbh;
				lblock.y1=by*lbh;
				lblock.y2=MINVAR(lblock.y1+lbh, ih);
				for(int bx=0;bx<lbx;++bx)
				{
					T25ParamsPacked *pp=(T25ParamsPacked*)array_at(&params, lbx*by+bx);
					T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
					int ng=(lbw+param.gwidth-1)/param.gwidth;
					lblock.x1=bx*lbw;
					lblock.x2=MINVAR(lblock.x1+lbw, iw);
					for(int kg=0;kg<ng;++kg)
					{
						int x1=lblock.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, lblock.x2);
						int success=t25_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2, 0);
						if(!success)
							LOG_ERROR("t25_prepblock error");
						for(int kx=x1;kx<x2;++kx)
						{
#if 0
							{//
								static int counter=0;
								unsigned char su=buf2[(iw*ky+kx)<<2|kc];
								int si=(char)(su+128);//signed

								int sz=si;

								int neg=sz<0;//zigzag code
								sz^=-neg;
								sz+=neg;
								sz<<=1;
								sz|=neg;

								int sf=sz;
								int s128=sf==1;
								sf-=sf>1;
								sf+=s128;

								int sg=sf;
								sg^=sg>>1;//to gray code

								if(counter>=10000&&counter<11000)
								{
									if(counter==10000)
										printf("unsigned, signed, zigzag, fixed, grey\n");
									print_bin8((unsigned)su);
									printf(" ");
									print_bin8((unsigned)si);
									printf(" ");
									print_bin8((unsigned)sz);
									printf(" ");
									print_bin8((unsigned)sf);
									printf(" ");
									print_bin8((unsigned)sg);
									printf("\n");
								}
								++counter;
							}
#endif
							ac_enc(&ec, buf2[(iw*ky+kx)<<2|kc], CDF2, 256, 1);
						}
					}
				}
			}
			ac_enc_flush(&ec);
		}
		ansbookmarks[kc]=(int)list.nobj;

		array_free(&params);
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	int chsizes[]=
	{
		ansbookmarks[0]-overhead[3]    -overhead[0],
		ansbookmarks[1]-ansbookmarks[0]-overhead[1],
		ansbookmarks[2]-ansbookmarks[1]-overhead[2],
	};
	if(loud)
	{
		int totaloverhead=overhead[0]+overhead[1]+overhead[2]+overhead[3], totalch=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("Total    %7d  %lf\n", totaloverhead+totalch, 3.*res/list.nobj);
		printf("Overhead %7d\n", totaloverhead);
		printf("Red      %7d  %lf\n", chsizes[0], (double)res/chsizes[0]);
		printf("Green    %7d  %lf\n", chsizes[1], (double)res/chsizes[1]);
		printf("Blue     %7d  %lf\n", chsizes[2], (double)res/chsizes[2]);
	}

	dlist_clear(&list);
	free(buf2);
	free(CDF0);
	free(CDF2);
	return 1;
}
int t25_decode(const unsigned char *data, size_t srclen, int iw, int ih, int *blockw, int *blockh, int use_ans, unsigned char *buf, int loud)
{
	const short predidx[]={0, 11, 11+PW2_NPARAM}, predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM};
	const int cdflen=768LL*sizeof(short), overhead=12LL+predlen[3]*sizeof(short)+cdflen;
	int res=iw*ih;
	
	if(srclen<=overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead||ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short predparams[22+PW2_NPARAM];
	memcpy(predparams, data+12, predlen[3]*sizeof(short));
	memcpy(CDF0, data+12+predlen[3]*sizeof(short), cdflen);
	int lbx[]=
	{
		(iw+blockw[0]-1)/blockw[0],
		(iw+blockw[1]-1)/blockw[1],
		(iw+blockw[2]-1)/blockw[2],
	};
	int lby[]=
	{
		(ih+blockh[0]-1)/blockh[0],
		(ih+blockh[1]-1)/blockh[1],
		(ih+blockh[2]-1)/blockh[2],
	};
	T25ParamsPacked const *params[]=
	{
		(T25ParamsPacked const*)(data+overhead),
		(T25ParamsPacked const*)(data+ansbookmarks[0]),
		(T25ParamsPacked const*)(data+ansbookmarks[1]),
	};
	int pcount[]=
	{
		lbx[0]*lby[0],
		lbx[1]*lby[1],
		lbx[2]*lby[2],
	};
	Rect block;
	for(int kc=0;kc<3;++kc)
	{
		const unsigned char
			*srcstart=(const unsigned char*)(params[kc]+pcount[kc]),
			*srcend=data+ansbookmarks[kc];
		if(use_ans)
		{
			ANSCoder ec;
			ans_dec_init(&ec, srcstart, srcend);
			//ans_dec_init(&ec, data+(kc?ansbookmarks[kc-1]:overhead), data+ansbookmarks[kc]);
			//unsigned state;
			//srcptr=data+ansbookmarks[kc];
			//srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
			//srcptr-=4;
			//if(srcptr<srcstart)
			//	LOG_ERROR("ANS buffer overflow");
			//memcpy(&state, srcptr, 4);

			for(int ky=0;ky<ih;++ky)
			{
				int by=ky/blockh[kc];
				block.y1=by*blockh[kc];
				block.y2=MINVAR(block.y1+blockh[kc], ih);
				for(int bx=0;bx<lbx[kc];++bx)
				{
					T25ParamsPacked const *pp=params[kc]+lbx[kc]*by+bx;
					T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
					int ng=(blockw[kc]+param.gwidth-1)/param.gwidth;
					block.x1=bx*blockw[kc];
					block.x2=MINVAR(block.x1+blockw[kc], iw);
					for(int kg=0;kg<ng;++kg)
					{
						int x1=block.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, block.x2);
						int success=t25_prepblock(buf, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2, 0);
						if(!success)
							LOG_ERROR("t25_prepblock error");
						for(int kx=x1;kx<x2;++kx)
						{
							buf[(iw*ky+kx)<<2|kc]=ans_dec(&ec, CDF2, 256);
#if 0
							unsigned c=(unsigned short)state;
							int sym=0;
					
							//if(kc==0&&ky==2&&kx==512)//
							//	printf("");

							int L=0, R=256, found=0;
							while(L<=R)
							{
								sym=(L+R)>>1;
								if(CDF2[sym]<c)
									L=sym+1;
								else if(CDF2[sym]>c)
									R=sym-1;
								else
								{
									found=1;
									break;
								}
							}
							if(found)
								for(;sym<256-1&&CDF2[sym+1]==c;++sym);
							else
								sym=L+(L<256&&CDF2[L]<c)-(L!=0);
								//sym=L+(L<256&&CDF2[L]<c)-1;

							buf[(iw*ky+kx)<<2|kc]=(unsigned char)sym;

							unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
						
							debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
							state=freq*(state>>16)+c-cdf;//update
							if(state<0x10000)//renorm
							{
								state<<=16;
								if(srcptr-2>=srcstart)
								{
									srcptr-=2;
									memcpy(&state, srcptr, 2);
								}
							}
#endif
						}
					}
				}
			}
		}
		else
		{
			ArithmeticCoder ec;
			ac_dec_init(&ec, srcstart, srcend);

			for(int ky=0;ky<ih;++ky)
			{
				int by=ky/blockh[kc];
				block.y1=by*blockh[kc];
				block.y2=MINVAR(block.y1+blockh[kc], ih);
				for(int bx=0;bx<lbx[kc];++bx)
				{
					T25ParamsPacked const *pp=params[kc]+lbx[kc]*by+bx;
					T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
					int ng=(blockw[kc]+param.gwidth-1)/param.gwidth;
					block.x1=bx*blockw[kc];
					block.x2=MINVAR(block.x1+blockw[kc], iw);
					for(int kg=0;kg<ng;++kg)
					{
						int x1=block.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, block.x2);
						int success=t25_prepblock(buf, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2, 0);
						if(!success)
							LOG_ERROR("t25_prepblock error");
						for(int kx=x1;kx<x2;++kx)
							buf[(iw*ky+kx)<<2|kc]=ac_dec(&ec, CDF2, 256, 1);
					}
				}
			}
		}
	}
	free(CDF0);
	free(CDF2);

#ifndef T25_OPT_PRED
	apply_transforms_inv(buf, iw, ih);
#else
	addbuf(buf, iw, ih, 3, 4, 128);
	memcpy(jxlparams_i16,         predparams+predidx[0], predlen[0]*sizeof(short));
	memcpy(pw2_params+PW2_NPARAM, predparams+predidx[1], predlen[1]*sizeof(short));
	memcpy(jxlparams_i16+22,      predparams+predidx[2], predlen[2]*sizeof(short));
	pred_opt_apply((char*)buf, iw, ih, 0);
	//pred_jxl_apply((char*)buf, iw, ih, jxlparams, 0);
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif

	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	return 1;
}

//	#define T26_OPT_PRED

//test26: T16 with arithmetic range coder
void t25_calchist(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, unsigned *hist);
void t26_normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	if(nsymbols)
	{
		int sum=0;
		for(int sym=0;sym<nlevels;++sym)
		{
			int qfreq=((long long)srchist[sym]*0xFFFF)/nsymbols;
			CDF[sym]=sum;
			sum+=qfreq;
		}
	}
	else//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
	}
}
void t25_addhist(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0a, int x0b, int y0, int maxinc, unsigned *CDF2);
int t26_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int iw, int ih, int kc, int x1, int x2, int y, T26Params const *p, unsigned *CDF2, int rec)
{
	int overflow=0;
	int sum, cdf1, f1, f2, freq;
	memset(CDF2, 0, 257*sizeof(unsigned));
	if(p->mtop)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x2+p->mright, y-p->mtop, y, x1, x2, y, p->maxinc, CDF2);
	if(p->mleft)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x1, y, y+1, x1, x2, y, p->maxinc, CDF2);

	if(CDF2[256])
	{
		sum=0;
		for(int sym=0;sym<256;++sym)
		{
			if(sym==0x80)//
				sym=0x80;

			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize
			
			if(f2<f1)
				freq=f2+(int)(((long long)f1-f2)*(0xFF-p->alpha)/0xFF);//blend
			else
				freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

			int f3=freq;//

			freq=(int)((long long)freq*0xFF00>>16)+1;//guard
			//freq=CLAMP(0, freq, 0xFF01);

			if(rec)//
				printf("%3d 0x%04X 0x%04X 0x%04X 0x%04X\n", sym, f1, f2, f3, freq);

			if(freq<0||freq>0xFF01)
			{
				printf("Impossible freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
				return 0;
			}
				//LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000&&sym<255)
			{
				if(!rec)//
				{
					t26_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, y, p, CDF2, 1);//
					//for(int k=0;k<=sym;++k)
					//	printf("%3d 0x%04X 0x%04X\n", k, CDF0[k], CDF2[k]);
					printf("ANS CDF sym 0x%02X sum 0x%04X  freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", sym, sum, freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
					//printf("ANS CDF sym 0x%02X sum 0x%04X freq 0x%04X\n", sym, sum, freq);
				}
				return 0;
			}
		}
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
	}
	CDF2[256]=0x10000;
	return 1;
}
int t26_encode(const unsigned char *src, int iw, int ih, T26Params const *params, int use_ans, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifndef T26_OPT_PRED
	apply_transforms_fwd(buf2, iw, ih);
#else
	addbuf(buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd((char*)buf2, iw, ih);
	pred_opt_opt_v6((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v5((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v4((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v3((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v2((char*)buf2, iw, ih);
#endif
	short predparams[22+PW2_NPARAM];
	const short
		predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM},
		predidx[]={0, 11, 11+PW2_NPARAM};
	memcpy(predparams+predidx[0], jxlparams_i16,         predlen[0]*sizeof(short));
	memcpy(predparams+predidx[1], pw2_params+PW2_NPARAM, predlen[1]*sizeof(short));
	memcpy(predparams+predidx[2], jxlparams_i16+22,      predlen[2]*sizeof(short));
#ifdef T26_OPT_PRED
	//pred_jxl_opt_v2((char*)buf2, iw, ih, jxlparams_i16, loud);
#if 1
	if(loud)
		pred_opt_printparam();
#endif
	pred_opt_apply((char*)buf2, iw, ih, 1);
	//pred_jxl_apply((char*)buf2, iw, ih, jxlparams_i16, 1);
	addbuf(buf2, iw, ih, 3, 4, 128);
#endif
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, CDF2);
		t26_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int bookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, predparams, (predlen[3])*sizeof(short));
	//dlist_push_back(&list, jxlparams_i16, 33*sizeof(short));
	dlist_push_back(&list, CDF0, 768*sizeof(short));
	
	for(int kc=0;kc<3;++kc)//for each channel
	{
		T26Params const *p=params+kc;
		int gxcount=(iw+p->gwidth-1)/p->gwidth;

		if(use_ans)
		{
			ANSCoder ec;
			ans_enc_init(&ec, &list);

			for(int ky=ih-1;ky>=0;--ky)//for each row
			{
				for(int bx=gxcount-1;bx>=0;--bx)//for each group
				{
					//if(kc==0&&bx==1&&ky==0)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					t26_prepblock(buf2, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);

					//if(kc==0&&bx==0&&by==0)
					//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

					//encode group
					for(int kx=x2-1;kx>=x1;--kx)//for each pixel
					{
						if(kc==0&&ky==0&&bx==0&&kx==0)//
							printf("");

						ans_enc(&ec, buf2[(iw*ky+kx)<<2|kc], CDF2, 256);
					}
				}
			}
			ans_enc_flush(&ec);
		}
		else
		{
			ArithmeticCoder ec;
			ac_enc_init(&ec, &list);

			//unsigned state_lo=0, state_hi=0xFFFFFFFF;
			//unsigned cache=0;
			//int nbits=0;
			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					//if(kc==0&&bx==1&&ky==0)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					t26_prepblock(buf2, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);

					//if(kc==0&&bx==0&&by==0)
					//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

					//encode group
					for(int kx=x1;kx<x2;++kx)//for each pixel
					{
						//if(kc==2&&ky2==0&&kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);
					
						//if(acval&&acval->count==2525)//
						//	printf("");
						//ac_enc_renorm(&state_lo, &state_hi, &cache, &nbits, &list);
	#if 0
						unsigned rlo, rhi, mingap;
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							//if(mingap>0&&state_lo+0x10000>state_hi)
							//	printf("");
							if(state_lo<state_hi&&mingap>0)
								break;
							if(nbits>=32)
							{
								dlist_push_back(&list, &cache, 4);
								cache=0;
								nbits=0;
							}
							cache|=(state_lo&0x80000000)>>nbits;//cache is written MSB -> LSB
							++nbits;
							//cache|=(state_lo>>31)<<32-nbits;

							state_lo<<=1;//shift out MSB
							state_hi<<=1;
							state_hi|=1;

	#if 0
							cache<<=1;
							cache|=state_lo>>31;
							state_lo<<=1;//shift out MSB
							state_hi<<=1;
							state_hi|=1;
							++nbits;
	#endif
						}
	#endif
					
						ac_enc(&ec, buf2[(iw*ky+kx)<<2|kc], CDF2, 256, 1);
					}
				}
			}
			ac_enc_flush(&ec);
#if 0
			int k2=0;
			do//flush
			{
				while(nbits<32)
				{
					cache|=(state_lo&0x80000000)>>nbits;//cache is written MSB -> LSB
					++nbits;
					++k2;

					state_lo<<=1;//shift out MSB
					state_hi<<=1;
					state_hi|=1;
				}
				dlist_push_back(&list, &cache, 4);
				cache=0;
				nbits=0;
			}while(k2<32);
#endif
		}
		bookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, bookmarks, 12);
	
	int overhead=12+(int)(768*sizeof(short));
	int ch[]=
	{
		bookmarks[0]-overhead,
		bookmarks[1]-bookmarks[0],
		bookmarks[2]-bookmarks[1],
	};
	//if(csizes)
	//{
	//	csizes[0]=ch[0];
	//	csizes[1]=ch[1];
	//	csizes[2]=ch[2];
	//}
	if(loud)
	{
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("Total    %7d  %lf\n", bookmarks[2], 3.*res/bookmarks[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf\n", ch[0], (double)res/ch[0]);
		printf("Green    %7d  %lf\n", ch[1], (double)res/ch[1]);
		printf("Blue     %7d  %lf\n", ch[2], (double)res/ch[2]);
	}

	dlist_clear(&list);
	free(CDF2);
	free(CDF0);
	free(buf2);
	return 1;
}
int t26_decode(const unsigned char *data, size_t srclen, int iw, int ih, T26Params const *params, int use_ans, unsigned char *buf, int loud)
{
	const short predidx[]={0, 11, 11+PW2_NPARAM}, predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM};
	const int cdflen=768LL*sizeof(short), overhead=12LL+predlen[3]*sizeof(short)+cdflen;
	int res=iw*ih;
	
	const unsigned char *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}

	unsigned bookmarks[3];
	memcpy(bookmarks, data, 12);
	if(bookmarks[2]<(unsigned)overhead||bookmarks[2]>srclen)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short predparams[22+PW2_NPARAM];
	memcpy(predparams, data+12, predlen[3]*sizeof(short));
	memcpy(CDF0, data+12+predlen[3]*sizeof(short), cdflen);

#ifdef AC_VALIDATE
	acval_idx=0;
#endif
	for(int kc=0;kc<3;++kc)//for each channel
	{
		T26Params const *p=params+kc;
		int gxcount=(iw+p->gwidth-1)/p->gwidth;

		if(use_ans)
		{
			ANSCoder ec;
			ans_dec_init(&ec, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			//unsigned state_lo=0, state_hi=0xFFFFFFFF, code, cache;
			//int nbits=32;
			//srcptr=kc?data+bookmarks[kc-1]:data+overhead;
			//srcend=data+bookmarks[kc];

			//if(ec.srcend-ec.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ec.code, ec.srcptr, 4);
			//ec.srcptr+=4;
			//
			//if(ec.srcend-ec.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ec.cache, ec.srcptr, 4);
			//ec.srcptr+=4;

			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					//if(kc==0&&ky==0&&bx==1)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					int success=t26_prepblock(buf, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);
					if(!success)
						LOG_ERROR("t26_prepblock error");
					for(int kx=x1;kx<x2;++kx)//for each pixel
					{
						buf[(iw*ky+kx)<<2|kc]=ans_dec(&ec, CDF2, 256);
#if 0
						unsigned rlo, rhi, mingap;

						//if(acval_idx==2525)//
						//	printf("");
					
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							if(state_lo<state_hi&&mingap>0)
								break;
							if(!nbits)
							{
								if(srcend-srcptr<4)
								{
#ifdef AC_VALIDATE
									printf("buffer overflow\n");
									acval_dump();
#endif
									LOG_ERROR("buffer overflow");
								}
								memcpy(&cache, srcptr, 4);
								srcptr+=4;

								nbits=32;
							}
							--nbits;
							code<<=1;//shift out MSB		cache is read MSB -> LSB
							code|=(unsigned)(cache>>nbits&1);

							state_lo<<=1;
							state_hi<<=1;
							state_hi|=1;
						}
#if 0
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							if(mingap>0)
								break;
							if(!nbits)
							{
								nbits+=32;
								cache<<=32;

								if(srcptr+4>=srcend)
									LOG_ERROR("buffer overflow");
								memcpy(&cache, srcptr, 4);
								srcptr+=4;
							}
							--nbits;
							code<<=1;//shift out MSB
							code|=(unsigned)(cache>>nbits&1);

							state_lo<<=1;
							state_hi<<=1;
							state_hi|=1;
						}
#endif
						int sym=0;
						int L=0, R=256, found=0;
						unsigned code2;
						while(L<=R)
						{
							sym=(L+R)>>1;
							code2=state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym]>>16);
							if(code2<code)
								L=sym+1;
							else if(code2>code)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(found)
							for(;sym<256-1&&state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym+1]>>16)==code;++sym);
						else
							sym=L+(L<256&&state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym+1]>>16)<code)-(L!=0);
#if 0
						//unsigned c=(unsigned)((((long long)(code-state_lo)<<16)+((state_hi-state_lo)>>1))/(state_hi-state_lo));
						unsigned c=(unsigned)(((long long)(code-state_lo)<<16)/(state_hi-state_lo));
						int sym=0;
					
						//if(kc==0&&ky==2&&kx==512)//
						//	printf("");

						int L=0, R=256, found=0;
						while(L<=R)
						{
							sym=(L+R)>>1;
							if(CDF2[sym]<c)
								L=sym+1;
							else if(CDF2[sym]>c)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(found)
							for(;sym<256-1&&CDF2[sym+1]==c;++sym);
						else
							sym=L+(L<256&&CDF2[L]<c)-(L!=0);
#endif
						buf[(iw*ky+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf_start=CDF2[sym], cdf_end=CDF2[sym+1];
					
						rlo=state_lo+((unsigned long long)(state_hi-state_lo)*cdf_start>>16);
						rhi=state_lo+((unsigned long long)(state_hi-state_lo)*cdf_end  >>16);
						acval_dec(sym, cdf_start, cdf_end, state_lo, state_hi, rlo, rhi, cache, nbits, code);//
						state_lo=rlo;
						state_hi=rhi-1;//OBLIGATORY range leak guard
#endif
						//debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
						//state=freq*(state>>16)+c-cdf;//update
						//if(state<0x10000)//renorm
						//{
						//	state<<=16;
						//	if(srcptr-2>=srcstart)
						//	{
						//		srcptr-=2;
						//		memcpy(&state, srcptr, 2);
						//	}
						//}
					}
				}
			}
		}
		else
		{
			ArithmeticCoder ec;
			ac_dec_init(&ec, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					int success=t26_prepblock(buf, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);
					if(!success)
						LOG_ERROR("t26_prepblock error");
					for(int kx=x1;kx<x2;++kx)//for each pixel
						buf[(iw*ky+kx)<<2|kc]=ac_dec(&ec, CDF2, 256, 1);
				}
			}
		}
	}
	free(CDF0);
	free(CDF2);

#ifndef T25_OPT_PRED
	apply_transforms_inv(buf, iw, ih);
#else
	addbuf(buf, iw, ih, 3, 4, 128);
	memcpy(jxlparams_i16,         predparams+predidx[0], predlen[0]*sizeof(short));
	memcpy(pw2_params+PW2_NPARAM, predparams+predidx[1], predlen[1]*sizeof(short));
	memcpy(jxlparams_i16+22,      predparams+predidx[2], predlen[2]*sizeof(short));
	pred_opt_apply((char*)buf, iw, ih, 0);
	//pred_jxl_apply((char*)buf, iw, ih, jxlparams, 0);
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif

	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
#ifdef AC_VALIDATE
	array_free(&acval);
#endif
	return 1;
}


//T27: bitwise codec	X  bad
#if 0
//static float predctx2[1+24+24];//pred, 24*2 prev bits+hits
//static float weight1[3][11], weight2[24][49];
static DataType t27_weights_pred[24][18], t27_weights_bits[49],//includes bias
	t27_bits[48];//bit & hit history: 3 channels * 8 bit
int t27_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, CDF2);
		t26_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, CDF0, 768*sizeof(short));
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);

	float lr=0.001f;
	
	float csizes[24]={0};
	int hits[24]={0};

	//int ntrains=0;
	//printf("Enter number of training passes: ");
	//while(!scanf("%d", &ntrains));

	//for(int kt=0;kt<ntrains+1;++kt)
	{
		//int test=kt==ntrains;
		unsigned char *ptr=buf2;
		initialize((DataType*)t27_weights_pred, sizeof(t27_weights_pred)/sizeof(DataType), (DataType)sqrtf(14));
		initialize((DataType*)t27_weights_bits, sizeof(t27_weights_bits)/sizeof(DataType), (DataType)sqrtf(24));
		memset(t27_bits, 0, sizeof(t27_bits));
		int w2=iw*2;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ptr+=4)
			{
				if(ky==(ih>>1)&&kx==(iw>>1))//
					printf("");//

				int pnb[]=
				{
					ky-2>=0&&kx-2>=0?((int*)ptr)[-w2-2]:0,
					ky-2>=0&&kx-1>=0?((int*)ptr)[-w2-1]:0,
					ky-2>=0         ?((int*)ptr)[-w2  ]:0,
					ky-2>=0&&kx+1<iw?((int*)ptr)[-w2+1]:0,
					ky-2>=0&&kx+2<iw?((int*)ptr)[-w2+2]:0,

					ky-1>=0&&kx-2>=0?((int*)ptr)[-iw-2]:0,
					ky-1>=0&&kx-1>=0?((int*)ptr)[-iw-1]:0,
					ky-1>=0         ?((int*)ptr)[-iw  ]:0,
					ky-1>=0&&kx+1<iw?((int*)ptr)[-iw+1]:0,
					ky-1>=0&&kx+2<iw?((int*)ptr)[-iw+2]:0,

							 kx-2>=0?((int*)ptr)[   -2]:0,
							 kx-1>=0?((int*)ptr)[   -1]:0,
				};
				unsigned char *cnb=(unsigned char*)(pnb+_countof(pnb));
				const float gain=1.f/128;
#define NORMALIZE(X) (((X)-128)*gain-1)
				DataType predctx1[]=
				{
					0, 0, 0,//current
					(DataType)(ky-(ih>>1))/ih, (DataType)(kx-(iw>>1))/iw,//position
					NORMALIZE(cnb[-1*4]), NORMALIZE(cnb[-1*4+1]), NORMALIZE(cnb[-1*4+2]),//left
					NORMALIZE(cnb[(-5-1)*4]), NORMALIZE(cnb[(-5-1)*4+1]), NORMALIZE(cnb[(-5-1)*4+2]),//topleft
					NORMALIZE(cnb[-5*4]), NORMALIZE(cnb[-5*4+1]), NORMALIZE(cnb[-5*4+2]),//top
					NORMALIZE(cnb[(-5+1)*4]), NORMALIZE(cnb[(-5+1)*4+1]), NORMALIZE(cnb[(-5+1)*4+2]),//topright
				};
				//unsigned char bytes[3]={0};
				for(int kt=0;kt<100;++kt)
				{
				for(int kb=23;kb>=0;--kb)//MSB -> LSB
				{
					int byteidx=kb>>3, bitidx=kb&7;
					DataType fpred=t27_weights_pred[kb][_countof(predctx1)];
					for(int k=0;k<_countof(predctx1);++k)
						fpred+=t27_weights_pred[kb][k]*predctx1[k];
					fpred*=128;
					fpred+=128;
					float net_pred=fpred;
					fpred=CLAMP(0, fpred, 255);
					int pred=(int)fpred;

					DataType fp0=t27_weights_bits[48];
					for(int k=0;k<_countof(t27_bits);++k)
						fp0+=t27_weights_bits[k]*t27_bits[k];
					fp0*=0x8000;
					fp0+=0x8000;
					float net_p0=fp0;
					fp0=CLAMP(1, fp0, 0xFFFF);
					unsigned short p0=(unsigned short)fp0;

					int offset=(1<<bitidx)>>1;
					int bit=ptr[byteidx]-128-pred;
					bit+=offset;
					bit>>=bitidx;
					bit&=1;
					//int bit=(ptr[byteidx]-128-pred+offset)>>bitidx&1;

#if 0
					int temp=p0;
					temp-=0x8000;
					temp=(int)(0x8000*sqrt((double)abs(temp)/0x8000));
					temp=(int)(0x8000*sqrt((double)abs(temp)/0x8000));
					if(p0<0x8000)
						temp=-temp;
					temp+=0x8000;
					p0=(unsigned short)CLAMP(1, temp, 0xFFFF);
#endif
					//if(p0<0x8000)
					//	p0=1;
					//else if(p0>0x8000)
					//	p0=0xFFFF;
					//p0=(unsigned short)(0x10000*erf(10*((double)p0/0x10000-0.5)));//X
				
					if(!kt)
						abac_enc(&ctx, p0, bit);

					//bytes[byteidx]+=(bit<<bitidx);

					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);//loss1
					if(!kt)
						csizes[kb]+=bitsize;//

					int hit=bit?p0<0x8000:p0>0x8000;
					//char hit=(char)(ptr[byteidx]-128-pred+64);//X
					//hit=hit<128;//false if {0, ...127} or true if {128, ...255}//X

					if(!kt)
						hits[kb]+=hit;

					t27_bits[kb]=(DataType)(bit?1:-1);
					t27_bits[kb+24]=(DataType)(hit?1:-1);
					predctx1[byteidx]=(float)(ptr[byteidx]&0xFF<<bitidx)/256;
					//predctx1[byteidx]=(float)bytes[byteidx]/256;


					//backward pass
#if 1
//#define NAN_GUARD(X) (isfinite(X)||LOG_ERROR("%s", #X), X)
					//if(!ntrains||!test)
					{
						float dloss_dnet;
						if(net_p0>=1&&net_p0<0xFFFF)
						{
							float dloss1_dp0=bit?1/(1-fp0):-1/fp0;
							dloss_dnet=dloss1_dp0*0x8000*lr;
							t27_weights_bits[48]-=nan_guard(dloss_dnet);
							for(int k=0;k<_countof(predctx1);++k)
								t27_weights_bits[k]-=nan_guard(dloss_dnet*t27_bits[k]);
						}
						if(net_pred>=0&&net_pred<=255)
						{
							dloss_dnet=fpred-(ptr[byteidx]-128);
							if(dloss_dnet<0)//sgn(pred-px)
								dloss_dnet=-1;
							else if(dloss_dnet>0)
								dloss_dnet=1;
							else
								dloss_dnet=0;

							dloss_dnet*=0x80*lr;
							t27_weights_pred[kb][_countof(predctx1)]-=nan_guard(dloss_dnet);
							for(int k=0;k<_countof(predctx1);++k)
								t27_weights_pred[kb][k]-=nan_guard(dloss_dnet*predctx1[k]);
						}
					}
#endif
				}
				}
			}
		}
	}
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
				printf("C%d\n", k>>3);
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			//if(!((k+1)&7))
			//	printf("\n");
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	free(CDF2);
	free(CDF0);
	free(buf2);
	return 1;
}
#endif


//T28: Bayesian inference		first good result
int t28_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *qtree=(unsigned short*)malloc(765*sizeof(short));
	unsigned *chist=(unsigned*)malloc(256*sizeof(unsigned));
	//unsigned *htree=(unsigned*)malloc(255*sizeof(unsigned));
	if(!buf2||!qtree||!chist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		//isym^=(unsigned char)isym>>1;//to grey code	X

		buf2[k]=isym;
	}

	for(int kc=0;kc<3;++kc)
	{
		memset(chist, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, chist);

		//htree[0]=0;
		//for(int sym=128;sym<255;++sym)
		//	htree[0]+=chist[sym];
		//qtree[0]=htree[0]*0xFFFF/res;

		for(int step=256, idx=0;step>1;step>>=1)
		{
			for(int start=0;start<256;start+=step, ++idx)
			{
				int c0=0, c1=0;
				for(int sym=start;sym<start+(step>>1);++sym)
					c0+=chist[sym];
				for(int sym=start+(step>>1);sym<start+step;++sym)
					c1+=chist[sym];
				if(c0+c1)
					qtree[255*kc+idx]=(unsigned short)((long long)c1*0xFFFE/(c0+c1))+1;
				else
					qtree[255*kc+idx]=0;
			}
		}
	/*	int den=res;
		for(int step=256, idx=0;step>1;step>>=1)
		{
			for(int k=0;k<256;k+=step, ++idx)
			{
				int sum=0;
				for(int sym=k;sym<(step>>1);++sym)
					sum+=chist[sym];
			}
		}*/

		//normalize histogram (to fit in uint16)
		//int nonzerocount=0;
		//for(int sym=0;sym<256;++sym)
		//	nonzerocount+=chist[sym]!=0;
		//int mag=0x10000-nonzerocount;
		//for(int sym=0;sym<256;++sym)
		//	qhist[sym]=((long long)chist[sym]*mag)/res+(chist[sym]!=0);

		//print qtrees
#if 0
		for(int kb=7, idx=0;kb>=0;--kb)
		{
			int step=1<<(kb+1);
			//if(kb==1)
			//	printf("");
			printf("\nbit %d:\n", kb);
			for(int start=0;start<256;start+=step, ++idx)
			{
				int c0=0, c1=0;
				for(int sym=start;sym<start+(step>>1);++sym)
					c0+=chist[sym];
				for(int sym=start+(step>>1);sym<start+step;++sym)
					c1+=chist[sym];
				int p0, p1;
				if(c0+c1)
					p1=(unsigned short)((long long)c1*0xFFFE/(c0+c1))+1, p0=0x10000-p1;
				else
					p0=p1=0;
				//int p1=c0+c1?c1*0xFFFE/(c0+c1)+1:0;
				//int p0=p1?0x10000-p1:0;

				int kb2=0;
				for(;kb2<7-kb;++kb2)
				{
					int bit=start>>(7-kb2)&1;
					printf("%c", '0'+bit);
				}
				printf("?");
				for(;kb2<7;++kb2)
					printf("x");
				printf(" 0x%04X 0x%04X ", p0, p1);
				if(p0+p1)
				{
					int mid=p0*64>>16;
					for(int k=0;k<mid;++k)
						printf("0");
					for(int k=mid;k<64;++k)
						printf("1");
				}
				printf("\n");
			}
		}
#endif

		//print histograms	X
#if 0
		printf("C%d\n", kc);
		int vmax=0;
		for(int sym=0;sym<256;++sym)
		{
			if(vmax<chist[sym])
				vmax=chist[sym];
		}
		for(int sym=0;sym<256;++sym)
		{
			//int gsym=sym^(sym>>1);
			int freq=chist[sym];
			printf("%3d  ", sym);
			print_bin8(sym);
			printf(" 0x%04X ", freq);
			int nstars=freq*64/vmax;
			//int nstars=freq*64>>16;
			for(int k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
#endif
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, qtree, 765*sizeof(short));
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);

	float lr=0.001f;
	
	float csizes[24]={0};
	int hits[24]={0};

	//int test=kt==ntrains;
	unsigned char *ptr=buf2;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))//
			//	printf("");//

			for(int kc=2;kc>=0;--kc)
			{
				const unsigned short *tree=qtree+255*kc;
				int step=1;
				unsigned char sym=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					int p0=0x10000-tree[sym];
					int bit=ptr[kc]>>kb&1;
					ac_enc_bin(&ec, p0, bit);
					
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[kc<<3|kb]+=bitsize;//

					sym<<=1;
					sym|=bit;
					tree+=step;
					step<<=1;
				}
			}
		}
	}
	ac_enc_flush(&ec);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	//free(htree);
	free(chist);
	free(qtree);
	free(buf2);
	return 1;
}


//T29: T28 + alpha
int t29_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *qtree=(unsigned short*)malloc(765*sizeof(short));
	unsigned *chist=(unsigned*)malloc(256*sizeof(unsigned));
	unsigned short *dtree=(unsigned short*)malloc(765*sizeof(short));
	if(!buf2||!qtree||!chist||!dtree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		buf2[k]=isym;
	}

	for(int kc=0;kc<3;++kc)
	{
		memset(chist, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, chist);

		for(int step=256, idx=0;step>1;step>>=1)
		{
			for(int start=0;start<256;start+=step, ++idx)
			{
				int c0=0, c1=0;
				for(int sym=start;sym<start+(step>>1);++sym)
					c0+=chist[sym];
				for(int sym=start+(step>>1);sym<start+step;++sym)
					c1+=chist[sym];
				if(c0+c1)
					qtree[255*kc+idx]=(unsigned short)((long long)c1*0xFFFE/(c0+c1))+1;
				else
					qtree[255*kc+idx]=0;
			}
		}
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, qtree, 765*sizeof(short));
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);

	int alpha[24]=
	{
		//bits 0 ... 7
		0x0000, 0x4000, 0x8000, 0xA000, 0xC000, 0xD000, 0xE000, 0xF000,//C0
		0x0000, 0x0000, 0x0000, 0x0000, 0x1000, 0x2000, 0x3000, 0x4000,//C1
		0x0000, 0x4000, 0x8000, 0xA000, 0xC000, 0xD000, 0xE000, 0xF000,//C2
	};
	//int inc[24]={0};
	
	float csizes[24]={0};
	int hits[24]={0};

	//for(int train=1;train>=0;--train)
	{
		unsigned char *ptr=buf2;
		for(int k=0;k<765;++k)
			dtree[k]=0x8000;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ptr+=4)
			{
				//if(ky==(ih>>1)&&kx==(iw>>1))//
				//	printf("");//

				for(int kc=2;kc>=0;--kc)
				{
					unsigned short *tree=qtree+255*kc, *tree2=dtree+255*kc;
					int step=1;
					unsigned char sym=0;
					for(int kb=7;kb>=0;--kb)//MSB -> LSB
					{
						int bitidx=kc<<3|kb;

						//if(kc==1&&kb==5)//
						//	printf("");

						int p0_1=0x10000-tree[sym], p0_2=0x10000-tree2[sym];
						p0_2+=!p0_2;
						int p0=p0_1+(int)((long long)(p0_2-p0_1)*alpha[bitidx]>>16);

						int bit=ptr[kc]>>kb&1;
						//if(!train)
						{
							ac_enc_bin(&ec, p0, bit);
					
							int hit=bit?p0<0x8000:p0>0x8000;//
							hits[bitidx]+=hit;
							float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
							float bitsize=-log2f(p);
							csizes[bitidx]+=bitsize;//
						}

						if(p0_1!=p0_2)
						{
							if(bit==(p0_1<p0_2))
								alpha[bitidx]+=alpha[bitidx]<0xFFFF;
							else
								alpha[bitidx]-=alpha[bitidx]>1;

						/*	int inc=bit==(p0_1<p0_2);
							int a=alpha[bitidx];
							if(inc)
								a=0x10000-a;
							a=(int)((long long)a*a>>16);
							if(inc)
								a=0x10000-a;
							a=CLAMP(1, a, 0xFFFF);
							alpha[bitidx]=a;*/
						/*	if(bit==(p0_1<p0_2))//decrease alpha
							{
								int a=alpha[bitidx];
								a=(int)((long long)a*a>>16);
								alpha[bitidx]=CLAMP(1, a, 0xFFFF);
							}
							else//increase alpha
							{
								int a=alpha[bitidx];
								a=0x10000-a;
								a=(int)((long long)a*a>>16);
								a=0x10000-a;
								a=CLAMP(1, a, 0xFFFF);
								alpha[bitidx]=a;
							}*/
						/*	if(bit)//p0 should be small
							{
								if(p0_1<p0_2)//p0_1 was better, decrease alpha
									alpha[bitidx]=alpha[bitidx];
								else//p0_2 was better, increase alpha
									alpha[bitidx]=alpha[bitidx];
							}
							else//p0 should be large
							{
								if(p0_1<p0_2)//p0_2 was better, increase alpha
									alpha[bitidx]=alpha[bitidx];
								else//p0_1 was better, decrease alpha
									alpha[bitidx]=alpha[bitidx];
							}*/
						/*	int inc2=bit==(p0_1<p0_2);
							int a=alpha[bitidx];
							if(inc2)
								++inc[bitidx], a+=inc[bitidx];
							else
								--inc[bitidx], a-=inc[bitidx];
							alpha[bitidx]=CLAMP(1, a, 0xFFFF);*/
						/*	int a=alpha[bitidx];
							a+=(p0_1-p0_2)>>6;
							//if(bit)//need small p0
							//{
							//	if(p0_1<p0_2)//p0_1 was better, decrease alpha
							//	{
							//		a-=p0_2-p0_1;
							//	}
							//	else//p0_2 was better, increase alpha
							//	{
							//		a+=p0_1-p0_2;
							//	}
							//}
							//else//need large p0
							//{
							//	if(p0_1<p0_2)
							//	{
							//		a+=p0_1-p0_2;
							//	}
							//	else
							//	{
							//		a-=p0_2-p0_1;
							//	}
							//}
							alpha[bitidx]=CLAMP(1, a, 0xFFFF);//*/
						}

						tree2[sym]=bit<<15|tree2[sym]>>1;
						sym<<=1;
						sym|=bit;
						tree+=step;
						tree2+=step;
						step<<=1;
					}
				}
			}
		}
	}
	ac_enc_flush(&ec);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	free(dtree);
	free(chist);
	free(qtree);
	free(buf2);
	return 1;
}


#if 0
#define T30_NPRED 134
inline unsigned char clip(int px)
{
  if(px>255)return 255;
  if(px<  0)return   0;
  return px;
}
inline unsigned char clamp4(int px, const unsigned char n1, const unsigned char n2, const unsigned char n3, const unsigned char n4)
{
	int maximum=n1, minimum=n1;
	if(maximum<n2)maximum=n2;
	if(maximum<n3)maximum=n3;
	if(maximum<n4)maximum=n4;
	if(minimum>n2)minimum=n2;
	if(minimum>n3)minimum=n3;
	if(minimum>n4)minimum=n4;
	if(px<minimum)return minimum;
	if(px>maximum)return maximum;
	return px;
}
float t30_weight[T30_NPRED], t30_weight2[2];
void t30_init_weight(float *w, int count)
{
	for(int k=0;k<count;++k)
		w[k]=(float)((xoroshiro128_next()&0xFFFF)+0x8000)/0x10000;//[0.5, 1.5]
}
int t30_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *qtree=(unsigned short*)malloc(765*sizeof(short));
	unsigned *chist=(unsigned*)malloc(256*sizeof(unsigned));
	int dtreesize=255*T30_NPRED*2;
	int *dtree=(int*)malloc(dtreesize*sizeof(int));
	if(!buf2||!qtree||!chist||!dtree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		buf2[k]=isym;
	}

	for(int kc=0;kc<3;++kc)
	{
		memset(chist, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, chist);

		for(int step=256, idx=0;step>1;step>>=1)
		{
			for(int start=0;start<256;start+=step, ++idx)
			{
				int c0=0, c1=0;
				for(int sym=start;sym<start+(step>>1);++sym)
					c0+=chist[sym];
				for(int sym=start+(step>>1);sym<start+step;++sym)
					c1+=chist[sym];
				if(c0+c1)
					qtree[255*kc+idx]=(unsigned short)((long long)c1*0xFFFE/(c0+c1))+1;
				else
					qtree[255*kc+idx]=0;
			}
		}
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, qtree, 765*sizeof(short));
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);

	float csizes[24]={0};
	int hits[24]={0};

	float lr=0.001f;

	unsigned char *ptr=buf2;
	XOROSHIRO128_RESET();
	t30_init_weight(t30_weight, _countof(t30_weight));
	t30_init_weight(t30_weight2, _countof(t30_weight2));
	t30_weight2[0]*=10;

	for(int k=0;k<dtreesize;++k)
		dtree[k]=0;

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))//
			//	printf("");//

			for(int kc=2;kc>=0;--kc)
			{
#define ACCESS(SRC, C, X, Y) ((unsigned)(C)<3&&(unsigned)(X)<(unsigned)iw&&(unsigned)(Y)<(unsigned)ih?src[((iw*(Y)+(X))<<2)+C]:0)
				unsigned char
					WWWWWW=ACCESS(buf2, kc, kx-6, ky), WWWWW=ACCESS(buf2, kc, kx-5, ky), WWWW=ACCESS(buf2, kc, kx-4, ky), WWW=ACCESS(buf2, kc, kx-3, ky), WW=ACCESS(buf2, kc, kx-2, ky), W=ACCESS(buf2, kc, kx-1, ky),
					NWWWW=ACCESS(buf2, kc, kx-4, ky-1), NWWW=ACCESS(buf2, kc, kx-3, ky-1), NWW=ACCESS(buf2, kc, kx-2, ky-1), NW=ACCESS(buf2, kc, kx-1, ky-1), N=ACCESS(buf2, kc, kx, ky-1), NE=ACCESS(buf2, kc, kx+1, ky-1), NEE=ACCESS(buf2, kc, kx+2, ky-1), NEEE=ACCESS(buf2, kc, kx+3, ky-1), NEEEE=ACCESS(buf2, kc, kx+4, ky-1),
					NNWWW=ACCESS(buf2, kc, kx-3, ky-2), NNWW=ACCESS(buf2, kc, kx-2, ky-2), NNW=ACCESS(buf2, kc, kx-1, ky-2), NN=ACCESS(buf2, kc, kx, ky-2), NNE=ACCESS(buf2, kc, kx+1, ky-2), NNEE=ACCESS(buf2, kc, kx+2, ky-2), NNEEE=ACCESS(buf2, kc, kx+3, ky-2),
					NNNWW=ACCESS(buf2, kc, kx-2, ky-3), NNNW=ACCESS(buf2, kc, kx-1, ky-3), NNN=ACCESS(buf2, kc, kx, ky-3), NNNE=ACCESS(buf2, kc, kx+1, ky-3), NNNEE=ACCESS(buf2, kc, kx+2, ky-3),
					NNNNW=ACCESS(buf2, kc, kx-1, ky-4), NNNN=ACCESS(buf2, kc, kx, ky-4), NNNNE=ACCESS(buf2, kc, kx+1, ky-4),
					NNNNN=ACCESS(buf2, kc, kx, ky-5), NNNNNN=ACCESS(buf2, kc, kx, ky-6),
					WWp1=ACCESS(buf2, kc-1, kx-2, ky), Wp1=ACCESS(buf2, kc-1, kx-1, ky), p1=ACCESS(buf2, kc-1, kx, ky), NWp1=ACCESS(buf2, kc-1, kx-1, ky+1), Np1=ACCESS(buf2, kc-1, kx, ky+1), NEp1=ACCESS(buf2, kc-1, kx+1, ky+1), NNp1=ACCESS(buf2, kc-1, kx, ky+2),
					WWp2=ACCESS(buf2, kc-2, kx-2, ky), Wp2=ACCESS(buf2, kc-2, kx-1, ky), p2=ACCESS(buf2, kc-2, kx, ky), NWp2=ACCESS(buf2, kc-2, kx-1, ky+1), Np2=ACCESS(buf2, kc-2, kx, ky+1), NEp2=ACCESS(buf2, kc-2, kx+1, ky+1), NNp2=ACCESS(buf2, kc-2, kx, ky+2),
					
					NNNEEE=ACCESS(buf2, kc, kx+3, ky-3),
					NNEp1=ACCESS(buf2, kc-1, kx+1, ky-2),
					NNEp2=ACCESS(buf2, kc-2, kx+1, ky-2),
					NEEp1=ACCESS(buf2, kc-1, kx+2, ky-1),
					NNEEp1=ACCESS(buf2, kc-1, kx+2, ky-2),
					NEEp2=ACCESS(buf2, kc-2, kx+2, ky-1),
					NNEEp2=ACCESS(buf2, kc-2, kx+2, ky-2),
					NNWp1=ACCESS(buf2, kc-1, kx-1, ky-2),
					NNWp2=ACCESS(buf2, kc-2, kx-1, ky-2),
					NNNp1=ACCESS(buf2, kc-1, kx, ky-3),
					NNNp2=ACCESS(buf2, kc-2, kx, ky-3),
					WWWp1=ACCESS(buf2, kc-1, kx-3, ky),
					WWWp2=ACCESS(buf2, kc-2, kx-3, ky),
					NNNWp1=ACCESS(buf2, kc-1, kx-1, ky-3),
					NNNWp2=ACCESS(buf2, kc-2, kx-1, ky-3),
					NNNEp1=ACCESS(buf2, kc-1, kx+1, ky-3),
					NNNEp2=ACCESS(buf2, kc-2, kx+1, ky-3),
					NNNNp1=ACCESS(buf2, kc-1, kx, ky-4),
					NNNNNNp1=ACCESS(buf2, kc-1, kx, ky-6),
					NNNNp2=ACCESS(buf2, kc-2, kx, ky-4),
					NNNNNNp2=ACCESS(buf2, kc-2, kx, ky-6),
					WWWWp1=ACCESS(buf2, kc-1, kx-4, ky),
					WWWWp2=ACCESS(buf2, kc-2, kx-4, ky),
					WWWWWWp1=ACCESS(buf2, kc-1, kx-6, ky),
					WWWWWWp2=ACCESS(buf2, kc-2, kx-6, ky),
					NNNWWW=ACCESS(buf2, kc, kx-3, ky-3),
					NNNWWWW=ACCESS(buf2, kc, kx-4, ky-3),
					NEEEEEE=ACCESS(buf2, kc, kx+6, ky-1),
					NWWp1=ACCESS(buf2, kc+1, kx-2, ky-1),
					NWWp2=ACCESS(buf2, kc+2, kx-2, ky-1),
					NNWWp1=ACCESS(buf2, kc+1, kx-2, ky-2),
					NNWWp2=ACCESS(buf2, kc+2, kx-2, ky-2);

				short cm[T30_NPRED]=
				{
					clamp4(N+p1-Np1, W, NW, N, NE),
					clamp4(N+p2-Np2, W, NW, N, NE),
					(W + clamp4(NE*3 - NNE*3 + NNNE, W, N, NE, NEE)) / 2,
					clamp4((W + clip(NE*2 - NNE))/2, W, NW, N, NE),
					(W + NEE)/2,
					((WWW - 4*WW + 6*W + (NE*4 - NNE*6 + NNNE*4 - NNNNE))/4),
					((-WWWW + 5*WWW - 10*WW + 10*W + clamp4(NE*4 - NNE*6 + NNNE*4 - NNNNE, N, NE, NEE, NEEE))/5),
					((-4*WW + 15*W + 10*(NE*3 - NNE*3 + NNNE) - (NEEE*3 - NNEEE*3 + NNNEEE))/20),
					((-3*WW + 8*W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE))/6),
					((W + (NE*2 - NNE))/2 + p1 - (Wp1 + (NEp1*2 - NNEp1))/2),
					((W + (NE*2 - NNE))/2 + p2 - (Wp2 + (NEp2*2 - NNEp2))/2),
					((-3*WW + 8*W + (NEE*2 - NNEE))/6 + p1 -(-3*WWp1 + 8*Wp1 + (NEEp1*2 - NNEEp1))/6),
					((-3*WW + 8*W + (NEE*2 - NNEE))/6 + p2 -(-3*WWp2 + 8*Wp2 + (NEEp2*2 - NNEEp2))/6),
					((W + NEE)/2 + p1 - (Wp1 + NEEp1)/2),
					((W + NEE)/2 + p2 - (Wp2 + NEEp2)/2),
					((WW + (NEE*2 - NNEE))/2 + p1 - (WWp1 + (NEEp1*2 - NNEEp1))/2),
					((WW + (NEE*2 - NNEE))/2 + p2 - (WWp2 + (NEEp2*2 - NNEEp2))/2),
					(WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1)),
					(WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2)),
					(W + N - NW),
					(W + N - NW + p1 - (Wp1 + Np1 - NWp1)),
					(W + N - NW + p2 - (Wp2 + Np2 - NWp2)),
					(W + NE - N),
					(N + NW - NNW),
					(N + NW - NNW + p1 - (Np1 + NWp1 - NNWp1)),
					(N + NW - NNW + p2 - (Np2 + NWp2 - NNWp2)),
					(N + NE - NNE),
					(N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1)),
					(N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2)),
					(N + NN - NNN),
					(N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1)),
					(N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2)),
					(W + WW - WWW),
					(W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1)),
					(W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2)),
					(W + NEE - NE),
					(W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1)),
					(W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2)),
					(NN + p1 - NNp1),
					(NN + p2 - NNp2),
					(NN + W - NNW),
					(NN + W - NNW + p1 - (NNp1 + Wp1 - NNWp1)),
					(NN + W - NNW + p2 - (NNp2 + Wp2 - NNWp2)),
					(NN + NW - NNNW),
					(NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1)),
					(NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2)),
					(NN + NE - NNNE),
					(NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1)),
					(NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2)),
					(NN + NNNN - NNNNNN),
					(NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1)),
					(NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2)),
					(WW + p1 - WWp1),
					(WW + p2 - WWp2),
					(WW + WWWW - WWWWWW),
					(WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1)),
					(WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2)),
					(N * 2 - NN + p1 - (Np1 * 2 - NNp1)),
					(N * 2 - NN + p2 - (Np2 * 2 - NNp2)),
					(W * 2 - WW + p1 - (Wp1 * 2 - WWp1)),
					(W * 2 - WW + p2 - (Wp2 * 2 - WWp2)),
					(N * 3 - NN * 3 + NNN),
					clamp4(N*3 - NN*3 + NNN, W, NW, N, NE),
					clamp4(W*3 - WW*3 + WWW, W, NW, N, NE),
					clamp4(N*2 - NN, W, NW, N, NE),
					((NNNNN - 6*NNNN + 15*NNN - 20*NN + 15*N + clamp4(W*4 - NWW*6 + NNWWW*4 - NNNWWWW, W, NW, N, NN))/6),
					((NNNEEE - 4*NNEE + 6*NE + (W*4 - NW*6 + NNW*4 - NNNW))/4),
					(((N + 3*NW)/4)*3 - ((NNW + NNWW)/2)*3 + (NNNWW*3 + NNNWWW)/4),
					((W*2 + NW) - (WW + 2*NWW) + NWWW),
					((W*2 - NW) + (W*2 - NWW) + N + NE)/4,
					(N + W + 1) >> 1,
					(NEEEE + NEEEEEE + 1) >> 1,
					(WWWWWW + WWWW + 1) >> 1,
					((W + N) * 3 - NW * 2) >> 2,
					N,
					NN,
					
					N + p1 - Np1,
					N + p2 - Np2,
					W + p1 - Wp1,
					W + p2 - Wp2,
					NW + p1 - NWp1,
					NW + p2 - NWp2,
					NE + p1 - NEp1,
					NE + p2 - NEp2,
					NN + p1 - NNp1,
					NN + p2 - NNp2,
					WW + p1 - WWp1,
					WW + p2 - WWp2,
					W + N - NW,
					W + N - NW + p1 - Wp1 - Np1 + NWp1,
					W + N - NW + p2 - Wp2 - Np2 + NWp2,
					W + NE - N,
					W + NE - N + p1 - Wp1 - NEp1 + Np1,
					W + NE - N + p2 - Wp2 - NEp2 + Np2,
					W + NEE - NE,
					W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1,
					W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2,
					N + NN - NNN,
					N + NN - NNN + p1 - Np1 - NNp1 + NNNp1,
					N + NN - NNN + p2 - Np2 - NNp2 + NNNp2,
					N + NE - NNE,
					N + NE - NNE + p1 - Np1 - NEp1 + NNEp1,
					N + NE - NNE + p2 - Np2 - NEp2 + NNEp2,
					N + NW - NNW,
					N + NW - NNW + p1 - Np1 - NWp1 + NNWp1,
					N + NW - NNW + p2 - Np2 - NWp2 + NNWp2,
					NE + NW - NN,
					NE + NW - NN + p1 - NEp1 - NWp1 + NNp1,
					NE + NW - NN + p2 - NEp2 - NWp2 + NNp2,
					NW + W - NWW,
					NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1,
					NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2,
					W*2 - WW,
					W*2 - WW + p1 - Wp1*2 + WWp1,
					W*2 - WW + p2 - Wp2*2 + WWp2,
					N*2 - NN,
					N*2 - NN + p1 - Np1*2 + NNp1,
					N*2 - NN + p2 - Np2*2 + NNp2,
					NW*2 - NNWW,
					NW*2 - NNWW + p1 - NWp1*2 + NNWWp1,
					NW*2 - NNWW + p2 - NWp2*2 + NNWWp2,
					NE*2 - NNEE,
					NE*2 - NNEE + p1 - NEp1*2 + NNEEp1,
					NE*2 - NNEE + p2 - NEp2*2 + NNEEp2,
					N*3 - NN*3 + NNN + p1 - Np1*3 + NNp1*3 - NNNp1,
					N*3 - NN*3 + NNN + p2 - Np2*3 + NNp2*3 - NNNp2,
					N*3 - NN*3 + NNN,
					(W + NE * 2 - NNE + 1) >> 1,
					(W + NE * 3 - NNE * 3 + NNNE+1) >> 1,
					(W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1*2 - NNEp1)/2,
					(W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2*2 - NNEp2)/2,
					NNE + NE - NNNE,
					NNE + W - NN,
					NNW + W - NNWW,
				};
				for(int k=0;k<T30_NPRED;++k)
					cm[k]=CLAMP(0, cm[k], 255);
				//unsigned char valid[T30_NPRED];
				//memset(valid, 1, sizeof(valid));

				unsigned short *tree=qtree+255*kc, *tree2=dtree;
				int step=1;
				int symlo=0, symhi=256;
				unsigned char sym=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					int bitidx=kc<<3|kb;

					//if(kc==1&&kb==5)//
					//	printf("");
					
					int symmid=(symlo+symhi)>>1;
					float p0_2=0, wsum=0;
					float temp;
					for(int k=0;k<T30_NPRED;++k)
					{
						int s2=cm[k]>>kb;
						int *c=tree2+(k<<(8-kb+1))+s2;
						int bit2=cm[k]>>kb&1;
						++c[bit2];
						float p0_temp=(float)c[0]/(c[0]+c[1]);
						p0_2+=(p0_temp-0.5f)*t30_weight[k];
						wsum+=t30_weight[k];
					}
					if(wsum)
					{
						p0_2*=0x10000/wsum;
						p0_2+=0x8000;
					}
					else
						p0_2=0x8000;

					int p0_1=0x10000-tree[sym];

					float wsum2=t30_weight2[0]+t30_weight2[1];
					float fp0=wsum2?(p0_1*t30_weight2[0]+p0_2*t30_weight2[1])/wsum2:p0_1;
					int p0=CLAMP(1, fp0, 0xFFFF);

					int bit=ptr[kc]>>kb&1;
					abac_enc(&ctx, p0, bit);
					
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[bitidx]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[bitidx]+=bitsize;//

					//backward pass
					if(fp0>=1&&fp0<=0xFFFF)
					{
						float dbs_dp0=bit?1/(1-fp0):-1/fp0;
						float dbs_dp02=dbs_dp0*(p0_2-fp0)/wsum2;

						t30_weight2[0]-=lr*dbs_dp0*(p0_1-fp0)/wsum2;
						t30_weight2[1]-=lr*dbs_dp02;

						for(int k=0;k<T30_NPRED;++k)
						{
							int s2=cm[k]>>kb;
							int *c=tree2+(k<<(8-kb+1))+s2;
							int bit2=cm[k]>>kb&1;
							float p0_temp=(float)c[0]/(c[0]+c[1]);
							t30_weight[k]-=lr*dbs_dp02*(p0_temp-p0_2)/wsum;
						}
					}

					tree2[sym]=bit<<15|tree2[sym]>>1;
					sym<<=1;
					sym|=bit;
					tree+=step;
					tree2+=step;
					step<<=1;
				}
			}
		}
	}
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	free(dtree);
	free(chist);
	free(qtree);
	free(buf2);
	return 1;
}
#endif


//T31: Adaptive Bayesian inference
int t31_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	int *ptree=(int*)malloc(255LL*3*2*sizeof(int));
	//int *ctree=(int*)malloc(255LL*3*2*sizeof(int));
	//unsigned *chist=(unsigned*)malloc(256*sizeof(unsigned));
	if(!buf2||!ptree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		//isym^=(unsigned char)isym>>1;//to grey code	X

		buf2[k]=isym;
	}

#if 0
	for(int kc=0;kc<3;++kc)
	{
		memset(chist, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, chist);
		for(int step=256, idx=0;step>1;step>>=1)
		{
			for(int start=0;start<256;start+=step, ++idx)
			{
				int idx2=(255*kc+idx)<<1;
				int c0=0, c1=0;
				for(int sym=start;sym<start+(step>>1);++sym)
					c0+=chist[sym];
				for(int sym=start+(step>>1);sym<start+step;++sym)
					c1+=chist[sym];
				
				if(c0+c1)
				{
					ctree[idx2  ]=(int)((long long)c0*0xFFFE/(c0+c1))+1;
					ctree[idx2|1]=(int)((long long)c1*0xFFFE/(c0+c1))+1;
				}
				else
				{
					ctree[idx2  ]=0x8000;
					ctree[idx2|1]=0x8000;
				}
				//ctree[idx2]=c0;
				//ctree[idx2|1]=c1;
			}
		}
	}
#endif
	for(int k=0;k<255*3*2;++k)
		ptree[k]=1;
	//double inc=1;

	DList list;
	dlist_init(&list, 1, 1024, 0);

	//dlist_push_back(&list, ctree, 255LL*3*2*sizeof(int));
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float csizes[24]={0};
	int hits[24]={0};

	unsigned char *ptr=buf2;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))//
			//	printf("");//

			for(int kc=2;kc>=0;--kc)
			{
				int *tree=ptree+(255*kc<<1);
				int step=1;
				unsigned char sym=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					int *c=tree+(sym<<1);
					//int p0=(int)(c[0]+c[1]?c[0]*0x10000/(c[0]+c[1]):0x8000);
					//int p0=c[0];
					int p0=c[0]+c[1]?(int)(((long long)c[0]<<16)/(c[0]+c[1])):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=ptr[kc]>>kb&1;
					ac_enc_bin(&ec, p0, bit);

					++c[bit];//CR 2.009992
					//inc*=1.01;//inflation		X  CR 0.074982
					
					//c[bit]+=(ptr-buf2)/(iw*ih);//X  CR 1.5
					
					//++c[bit];

					//if(c[0]+c[1]>0x100)
					//{
					//	c[0]>>=1;
					//	c[1]>>=1;
					//}
					
					//c[bit]-=c[bit]>0;
					
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[kc<<3|kb]+=bitsize;//

					sym<<=1;
					sym|=bit;
					tree+=step<<1;
					step<<=1;
				}
			}
		}
	}
	ac_enc_flush(&ec);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	//free(chist);
	//free(ctree);
	free(ptree);
	free(buf2);
	return 1;
}


//T32: Joint adaptive Bayesian inference	X  bad
#if 0
int t32_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	int treesize=((1<<24)-1)*2;
	int *ptree=(int*)malloc(treesize*sizeof(int));
	//int *ctree=(int*)malloc(255LL*3*2*sizeof(int));
	//unsigned *chist=(unsigned*)malloc(256*sizeof(unsigned));
	if(!buf2||!ptree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		isym^=(unsigned char)isym>>1;//to grey code	X

		buf2[k]=isym;
	}

	memset(ptree, 0, treesize*sizeof(int));
	//for(int k=0;k<treesize;++k)
	//	ptree[k]=1;

	DList list;
	dlist_init(&list, 1, 1024, 0);

	//dlist_push_back(&list, ctree, 255LL*3*2*sizeof(int));
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	int hits[24]={0};

	for(int train=1;train>=0;--train)//CHEAT
	{
		unsigned char *ptr=buf2;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ptr+=4)
			{
				int sym=0;
				int step=1;
				int *tree=ptree;
				//if(ky==(ih>>1)&&kx==(iw>>1))//
				//	printf("");//

				for(int kb=23;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					int *c=tree+(sym<<1);
					int p0=c[0]+c[1]?(int)(((long long)c[0]<<16)/(c[0]+c[1])):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=ptr[kb>>3]>>(kb&7)&1;
					if(!train)
					{
						abac_enc(&ctx, p0, bit);
						--c[bit];
					}
					else
						++c[bit];
				
					if(!train)
					{
						int hit=bit?p0<0x8000:p0>0x8000;//
						hits[kb]+=hit;
						float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
						float bitsize=-log2f(p);
						csizes[kb]+=bitsize;//
					}

					sym<<=1;
					sym|=bit;
					tree+=step<<1;
					step<<=1;
				}
			}
		}
	}
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	//free(chist);
	//free(ctree);
	free(ptree);
	free(buf2);
	return 1;
}
#endif


//T33: Adaptive Bayesian inference with circular buffer		X  bad
#if 0
	#define T33_BITLEN 1

#define T33_BUFLEN ((T33_BITLEN+7)>>3)
typedef struct T33BufferStruct
{
	unsigned short n[2], bitidx;
	unsigned char buf[T33_BUFLEN];
} T33Buffer;
int t33_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	T33Buffer *ptree=(T33Buffer*)malloc(255LL*3*sizeof(T33Buffer));
	if(!buf2||!ptree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	for(int k=0;k<(res<<2);++k)
	{
		char isym=buf2[k]-128;

		int neg=isym<0;//to zigzag code
		isym^=-neg;
		isym+=neg;
		isym<<=1;
		isym|=neg;

		isym-=(isym!=0)+(isym==1);//correction for 0x80 -> 0x01 -> 0xFF

		//isym^=(unsigned char)isym>>1;//to grey code	X

		buf2[k]=isym;
	}

	for(int k=0;k<255*3;++k)
	{
		T33Buffer *b=ptree+k;
		b->n[1]=b->n[0]=T33_BITLEN/2;
		b->bitidx=0;
		for(int k2=0;k2<T33_BUFLEN;++k2)
			b->buf[k2]=0x55;//initialize to alternating ones and zeros
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	//dlist_push_back(&list, ctree, 255LL*3*2*sizeof(int));
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	int hits[24]={0};

	unsigned char *ptr=buf2;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))//
			//	printf("");//

			for(int kc=2;kc>=0;--kc)
			{
				T33Buffer *tree=ptree+255*kc;
				int step=1;
				unsigned char sym=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					T33Buffer *b=tree+sym;
					int p0=b->n[0]+b->n[1]?(int)(((long long)b->n[0]<<16)/(b->n[0]+b->n[1])):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=ptr[kc]>>kb&1;
					abac_enc(&ctx, p0, bit);

					int bit0=b->buf[b->bitidx>>3]>>(b->bitidx&7)&1;
					++b->n[bit0];
					++b->n[bit];
					b->buf[b->bitidx>>3]&=~(1<<(b->bitidx&7));
					b->buf[b->bitidx>>3]|=(bit<<(b->bitidx&7));
					b->bitidx=(b->bitidx+1)%T33_BITLEN;
					
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[kc<<3|kb]+=bitsize;//

					sym<<=1;
					sym|=bit;
					tree+=step;
					step<<=1;
				}
			}
		}
	}
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
	}
	dlist_clear(&list);
	//free(chist);
	//free(ctree);
	free(ptree);
	free(buf2);
	return 1;
}
#endif


#define LR (int)(0.07*0x10000+0.5)


//T34: Adaptive Bayesian inference
#define NESTIMATORS 7
typedef struct NodeStruct
{
	int n[2];				//total counts
	unsigned short rec[6];	//p1 += (revealedbit-p1)/2^i

	//unsigned short rec6;//p1 = p1+(revealedbit-p1)*0.015625
	//unsigned short rec5;//p1 = p1+(revealedbit-p1)*0.03125
	//unsigned short rec4;//p1 = p1+(revealedbit-p1)*0.0625
	//unsigned short rec3;//p1 = p1+(revealedbit-p1)*0.125
	//unsigned short rec2;//p1 = p1+(revealedbit-p1)*0.25
	//unsigned short rec1;//p1 = p1+(revealedbit-p1)*0.5 = 0.5*revealedbit + 0.5*p1_prev	last 16 bits

	//int n65536_0, n65536_1;	//last 65536 symbols
	//int n4096_0, n4096_1;		//last 4096 symbols
	//int n256_0, n256_1;		//last 256 symbols
} Node;
#define HM_W 24
#define HM_H 16
int t34_weights[3*8*NESTIMATORS];
int t34_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	Node *tree=(Node*)malloc(255LL*3*sizeof(Node));
	if(!buf2||!tree)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float ssizes[HM_W*HM_H*3]={0}, ssdevs[HM_W*HM_H*3]={0};
	float csizes[24]={0};
	int hits[24]={0};
	
#if 1
	for(int k=0;k<3*8*NESTIMATORS;++k)//fixed 15.16 bit
		t34_weights[k]=0x8000;
	for(int k=0;k<255*3;++k)
	{
		Node *p=tree+k;
		p->n[0]=1;
		p->n[1]=1;
		for(int k=0;k<NESTIMATORS-1;++k)
			p->rec[k]=0x8000;
		//p->rec6=0x8000;
		//p->rec5=0x8000;
		//p->rec4=0x8000;
		//p->rec3=0x8000;
		//p->rec2=0x8000;
		//p->rec1=0x8000;
	}
	unsigned char *ptr=buf2;
	//for(int ky=ih-1;ky>=0;--ky)
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			for(int kc=2;kc>=0;--kc)
			{
				Node *treek=tree+255*kc;
				int *wk=t34_weights+8*NESTIMATORS*kc;
				int step=1;
				unsigned char sym=0, esym=0;
				for(int kb=7;kb>=0;--kb, wk+=NESTIMATORS)//MSB -> LSB
				{
					Node *node=treek+sym;
					long long sum=node->n[0]+node->n[1];
					int p0a[NESTIMATORS]=
					{
						sum?(int)(((long long)node->n[0]<<16)/sum):0x8000,
						//0x10000-node->rec6,
						//0x10000-node->rec5,
						//0x10000-node->rec4,
						//0x10000-node->rec3,
						//0x10000-node->rec2,
						//0x10000-node->rec1,
					};
					for(int k=0;k<NESTIMATORS-1;++k)
						p0a[k+1]=0x10000-node->rec[k];

					//fwd
					long long wsum=0;
					sum=0;
					for(int k=0;k<NESTIMATORS;++k)
					{
						sum+=(long long)p0a[k]*wk[k];
						wsum+=wk[k];
					}
					int p0=wsum?(int)(sum/wsum):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=ptr[kc]>>kb&1;
					ac_enc_bin(&ec, p0, bit);

					//bwd
					int p_bit=bit?0x10000-p0:p0;
					long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
					dL_dp0^=-bit;
					dL_dp0+=bit;
					for(int k=0;k<NESTIMATORS;++k)
					{
						int diff=p0a[k]-p0;//fixed 15.16 bit
						long long grad = dL_dp0*diff/wsum;
						long long wnew=LR*grad>>16;
						wnew=wk[k]-wnew;
						wnew=CLAMP(1, wnew, 0xFFFF);
						wk[k]=(int)wnew;
					}

					//update
					++node->n[bit];
					for(int k=0;k<NESTIMATORS-1;++k)
					{
						int temp=node->rec[k]+(((bit<<16)-node->rec[k])>>((k+1)<<1));
						node->rec[k]=CLAMP(1, temp, 0xFFFF);
					}
//					int temp;
//#define UPDATE_REC(N) temp=node->rec##N+(((bit<<16)-node->rec##N)>>(N<<1)), node->rec##N=CLAMP(1, temp, 0xFFFF);
//					UPDATE_REC(6);
//					UPDATE_REC(5);
//					UPDATE_REC(4);
//					UPDATE_REC(3);
//					UPDATE_REC(2);
//					UPDATE_REC(1);
					
#if 1
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[kc<<3|kb]+=bitsize;//

					int idx2=HM_W*(HM_H*kc+ky*HM_H/ih)+kx*HM_W/iw;//
					ssizes[idx2]+=bitsize;
					float fp0=(float)(p0-0x8000)/0x10000;
					ssdevs[idx2]+=fp0*fp0;//
#endif

					sym<<=1;
					sym|=bit;
					treek+=step;
					step<<=1;
				}
			}
		}
		if(loud==2)
		{
			static double csize_prev=0;
			double csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("Y%4d  CR%9lf  CR_row%9lf\n", ky, ((ky+1)*iw*3)/(csize+1), (iw*3+1)/(csize-csize_prev+1));
			csize_prev=csize;
		}
	}
#endif
#if 0
	int inc=1;
	int alpha[8];
	for(int k=0;k<8;++k)
		alpha[k]=0x10000;
	for(int k=0;k<255*3*2;++k)
		ptree[k]=1;
	const int blocksize=4;
	for(int y1=0;y1<ih;y1+=blocksize)
	{
		int y2=y1+blocksize;
		if(y2>ih)
			y2=ih;
		for(int x1=0;x1<iw;x1+=blocksize)
		{
			int x2=x1+blocksize;
			if(x2>iw)
				x2=iw;

			for(int k=0;k<255*3*2;++k)//reset dtree
				dtree[k]=1;

			for(int ky=y1;ky<y2;++ky)
			{
				for(int kx=x1;kx<x2;++kx)
				{
					for(int kc=2;kc>=0;--kc)
					{
						int *tree1=ptree+(255*kc<<1);
						int *tree2=dtree+(255*kc<<1);
						int step=1;
						unsigned char sym=buf2[(iw*ky+kx)<<2|kc], seenbits=0;
						int offset=64;
						for(int kb=7;kb>=0;--kb, offset>>=1)//MSB -> LSB
						{
							int *c=tree1+(seenbits<<1);
							int *d=tree2+(seenbits<<1);
							int p0_1=c[0]+c[1]?(int)(((long long)c[0]<<16)/(c[0]+c[1])):0x8000;
							int p0_2=d[0]+d[1]?(int)(((long long)d[0]<<16)/(d[0]+d[1])):0x8000;
							int p0=p0_1+(int)((long long)(p0_2-p0_1)*alpha[kb]>>16);
							p0=CLAMP(1, p0, 0xFFFF);
							int bit=sym>>kb&1;
							abac_enc(&ec, p0, bit);

							++c[bit];
							d[bit]+=inc;

							//int hit1=bit?p0_1<0x8000:p0_1>0x8000,
							//	hit2=bit?p0_2<0x8000:p0_2>0x8000;
							//alpha[kc]+=(hit1<hit2)-(hit1>hit2);

							int hit=bit?p0<0x8000:p0>0x8000;
							hits[kc<<3|kb]+=hit;//
							float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
							float bitsize=-log2f(p);
							csizes[kc<<3|kb]+=bitsize;//

							int idx2=HM_W*(HM_H*kc+ky*HM_H/ih)+kx*HM_W/iw;//
							ssizes[idx2]+=bitsize;
							float fp0=(float)(p0-0x8000)/0x10000;
							ssdevs[idx2]+=fp0*fp0;//

							seenbits<<=1;
							seenbits|=bit;
							tree1+=step<<1;
							tree2+=step<<1;
							step<<=1;
						}
					}
				}
			}
			//inc<<=1;//X
		}
	}
#endif
#if 0
	for(int k=0;k<255*3*2;++k)
		ptree[k]=1;
	unsigned char *ptr=buf2;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=4)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))//
			//	printf("");//

			for(int kc=2;kc>=0;--kc)
			{
				int *tree=ptree+(255*kc<<1);
				int step=1;
				unsigned char sym=0;
				int offset=64;
				for(int kb=7;kb>=0;--kb, offset>>=1)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					int *c=tree+(sym<<1);
					int p0=c[0]+c[1]?(int)(((long long)c[0]<<16)/(c[0]+c[1])):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=(ptr[kc]-offset)>>kb&1;
					abac_enc(&ec, p0, bit);

					++c[bit];
					
					int hit=bit?p0<0x8000:p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					float p=(unsigned short)(bit?0x10000-p0:p0)*(1.f/0x10000);
					float bitsize=-log2f(p);
					csizes[kc<<3|kb]+=bitsize;//

					sym<<=1;
					sym|=bit;
					tree+=step<<1;
					step<<=1;
				}
			}
		}
	}
#endif
	ac_enc_flush(&ec);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");
		
		if(loud==2)
		{
			float ubitsize=8*((float)iw/HM_W)*((float)ih/HM_H);
#if 0
			printf("Sdev maps:\n");
			for(int kc=0;kc<3;++kc)
			{
				for(int ky=0;ky<HM_H;++ky)
				{
					for(int kx=0;kx<HM_W;++kx)
						printf("%8.3f", sqrtf(ssdevs[HM_W*(HM_H*kc+ky)+kx]/ubitsize));
					printf("\n");
				}
				printf("\n");
			}
			printf("\n");
#endif

			printf("CR maps:\n");
			for(int kc=0;kc<3;++kc)
			{
				for(int ky=0;ky<HM_H;++ky)
				{
					for(int kx=0;kx<HM_W;++kx)
						printf("%8.3f", ubitsize/ssizes[HM_W*(HM_H*kc+ky)+kx]);
					printf("\n");
				}
				printf("\n");
			}
			printf("\n");
		}
	}
	dlist_clear(&list);
	free(tree);
	free(buf2);
	return 1;
}


	#define FIXEDPREC

#ifdef FIXEDPREC
typedef long long DataType;
#define FRACBITS 16
#define MAGBITS 12
#define ONE (1<<FRACBITS)
#define ONE_PERCENT ((0x28F5C28+(1<<(32-FRACBITS-1)))>>(32-FRACBITS))//0.01<<32
//#define ONE_PERCENT (int)(0.01*ONE+0.5)
#define ONE_PERMILLE ((0x418937+(1<<(32-FRACBITS-1)))>>(32-FRACBITS))//0.001<<32
#define MAXMAG (1<<(FRACBITS+FRACBITS-MAGBITS))
#define MUL(A, B) (int)((long long)(A)*(B)>>FRACBITS)
#define INITIALIZE(FIN, FOUT) sample(FIN, FOUT)
//#define LOAD(PX) ((PX)<<MAGBITS)
//#define STORE(VAL) (int)(((VAL)+(1<<(MAGBITS-1))-1)>>MAGBITS)
//#define LEARNING_RATE(G) MUL(G, ONE_PERMILLE)
//#define LEARNING_RATE(G) ((G)>>FRACBITS)
#define ABS llabs
#define FROMFLOAT(X) (int)((X)*ONE)
#else
typedef float DataType;
#define ONE 1
#define ONE_PERCENT 0.01f
#define ONE_PERMILLE 0.001f
#define MAXMAG 10
#define MUL(A, B) ((A)*(B))
#define INITIALIZE(FIN, FOUT) (float)sample(FIN, FOUT)/0xFFFF
//#define LOAD(PX) MUL(PX+0.5f, 1.f/256)
//#define STORE(VAL) (char)(MUL(VAL, 256)-0.5f)
//#define LEARNING_RATE(G) ((G)*0.001f)
#define ABS fabsf
#define FROMFLOAT(X) X
#endif
static void add_vec(DataType *dst, const DataType *a, const DataType *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b[k];
}
static void act(DataType *dst, const DataType *src, int count)
{
	for(int k=0;k<count;++k)
	{
		DataType val=src[k];
		DataType negpart=MUL(val, ONE_PERCENT);
		val=val>negpart?val:negpart;
		val=CLAMP(-MAXMAG, val, MAXMAG);
		dst[k]=val;
	}
}
static void act_dash(DataType *dst, const DataType *src, int count)
{
	for(int k=0;k<count;++k)
	{
		DataType val=src[k];
		if(val<-MAXMAG||val>MAXMAG)
			val=0;
		else
			val=val<0?ONE_PERCENT:ONE;
		dst[k]=val;
	}
}
static void mul_vec_scalar(DataType *dst, DataType *vec, DataType scalar, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=MUL(vec[k], scalar);
}
static void mul_vec_ew(DataType *dst, const DataType *v1, const DataType *v2, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=MUL(v1[k], v2[k]);
}
static void mul_vec_outer(DataType *dst, const DataType *left, const DataType *right, int lh, int rw)
{
	for(int ky=0;ky<lh;++ky)
	{
		for(int kx=0;kx<rw;++kx)
			dst[rw*ky+kx]=MUL(left[ky], right[kx]);
	}
}
static void mul_mat(DataType *dst, const DataType *m1, const DataType *m2, int h1, int w1h2, int w2)
{
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			DataType sum=0;
			for(int k=0;k<w1h2;++k)
				sum+=MUL(m1[w1h2*ky+k], m2[w2*k+kx]);
			dst[w2*ky+kx]=sum;
		}
	}
}
static void linear(DataType *dst, const DataType *mat, const DataType *vec, const DataType *bias, int win, int hout)//fixed 15.16 bit
{
	for(int ko=0;ko<hout;++ko)
	{
		DataType temp=bias?bias[ko]:0;
		for(int ki=0;ki<win;++ki)
			temp+=MUL(mat[win*ko+ki], vec[ki]);
		dst[ko]=temp;
	}
}


#if 0
#define T35_PREC 16	//fractional bits
DataType clip_64(DataType x)
{
	x=CLAMP(-(128<<T35_PREC), x, 127<<T35_PREC);
	return x;
}
DataType clamp4_64(DataType x, DataType a, DataType b, DataType c, DataType d)
{
	DataType vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmax<b)vmax=b;
	if(vmin>c)vmin=c;
	if(vmax<c)vmax=c;
	if(vmin>d)vmin=d;
	if(vmax<d)vmax=d;
	x=CLAMP(vmin, x, vmax);
	return x;
}
static void t35_getctx(const char *buf, int iw, int ih, int kx, int ky, int kc, DataType *ctx)
{
	//offsets are in NW direction
#define LOAD(CO, XO, YO) (unsigned)(kx-CO)<3u&&(unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?buf[(iw*(ky-YO)+kx-(XO))<<2|(kx-CO)]<<T35_PREC:0
	int WWWWWW  =LOAD(0,  6, 0),
		WWWWW   =LOAD(0,  5, 0),
		WWWW    =LOAD(0,  4, 0),
		WWW     =LOAD(0,  3, 0),
		WW      =LOAD(0,  2, 0),
		W       =LOAD(0,  1, 0),
		NWWWW   =LOAD(0,  4, 1),
		NWWW    =LOAD(0,  3, 1),
		NWW     =LOAD(0,  2, 1),
		NW      =LOAD(0,  1, 1),
		N       =LOAD(0,  0, 1),
		NE      =LOAD(0, -1, 1),
		NEE     =LOAD(0, -2, 1),
		NEEE    =LOAD(0, -3, 1),
		NEEEE   =LOAD(0, -4, 1),
		NEEEEEE =LOAD(0, -6, 1),
		NNWWW   =LOAD(0,  3, 2),
		NNWW    =LOAD(0,  2, 2),
		NNW     =LOAD(0,  1, 2),
		NN      =LOAD(0,  0, 2),
		NNE     =LOAD(0, -1, 2),
		NNEE    =LOAD(0, -2, 2),
		NNEEE   =LOAD(0, -3, 2),
		NNNWWWW =LOAD(0,  4, 3),
		NNNWWW  =LOAD(0,  3, 3),
		NNNWW   =LOAD(0,  2, 3),
		NNNW    =LOAD(0,  1, 3),
		NNN     =LOAD(0,  0, 3),
		NNNE    =LOAD(0, -1, 3),
		NNNEE   =LOAD(0, -2, 3),
		NNNEEE  =LOAD(0, -3, 3),
		NNNNW   =LOAD(0,  1, 4),
		NNNN    =LOAD(0,  0, 4),
		NNNNE   =LOAD(0, -1, 4),
		NNNNN   =LOAD(0,  0, 5),
		NNNNNN  =LOAD(0,  0, 6),
		WWWWWWp1=LOAD(1,  6, 0),
		WWWWp1  =LOAD(1,  4, 0),
		WWWp1   =LOAD(1,  3, 0),
		WWp1    =LOAD(1,  2, 0),
		Wp1     =LOAD(1,  1, 0),
		p1      =LOAD(1,  0, 0),
		NWWp1   =LOAD(1,  2, 1),
		NWp1    =LOAD(1,  1, 1),
		Np1     =LOAD(1,  0, 1),
		NEp1    =LOAD(1, -1, 1),
		NEEp1   =LOAD(1, -2, 1),
		NNWWp1  =LOAD(1,  2, 2),
		NNp1    =LOAD(1,  0, 2),
		NNEp1   =LOAD(1, -1, 2),
		NNEEp1  =LOAD(1, -2, 2),
		NNNWp1  =LOAD(1,  1, 3),
		NNNp1   =LOAD(1,  0, 3),
		NNNEp1  =LOAD(1, -1, 3),
		NNNNp1  =LOAD(1,  0, 4),
		NNNNNNp1=LOAD(1,  0, 6),
		WWWWWWp2=LOAD(2,  6, 0),
		WWWWp2  =LOAD(2,  4, 0),
		WWWp2   =LOAD(2,  3, 0),
		WWp2    =LOAD(2,  2, 0),
		Wp2     =LOAD(2,  1, 0),
		p2      =LOAD(2,  0, 0),
		NWWp2   =LOAD(2,  2, 1),
		NWp2    =LOAD(2,  1, 1),
		Np2     =LOAD(2,  0, 1),
		NEp2    =LOAD(2, -1, 1),
		NEEp2   =LOAD(2, -2, 1),
		NNWWp2  =LOAD(2,  2, 2),
		NNp2    =LOAD(2,  0, 2),
		NNEp2   =LOAD(2, -1, 2),
		NNEEp2  =LOAD(2, -2, 2),
		NNNWp2  =LOAD(2,  1, 3),
		NNNp2   =LOAD(2,  0, 3),
		NNNEp2  =LOAD(2, -1, 3),
		NNNNp2  =LOAD(2,  0, 4),
		NNNNNNp2=LOAD(2,  0, 6);
#undef LOAD
	int j=-1;
	ctx[++j] = clamp4_64(N + p1 - Np1, W, NW, N, NE);
	ctx[++j] = clamp4_64(N + p2 - Np2, W, NW, N, NE);
	ctx[++j] = (W + clamp4_64(NE * 3 - NNE * 3 + NNNE, W, N, NE, NEE)) / 2;
	ctx[++j] = clamp4_64((W + clip_64(NE * 2 - NNE)) / 2, W, NW, N, NE);
	ctx[++j] = (W + NEE) / 2;
	ctx[++j] = ((WWW - 4 * WW + 6 * W + (NE * 4 - NNE * 6 + NNNE * 4 - NNNNE)) / 4);
	ctx[++j] = ((-WWWW + 5 * WWW - 10 * WW + 10 * W + clamp4_64(NE * 4 - NNE * 6 + NNNE * 4 - NNNNE, N, NE, NEE, NEEE)) / 5);
	ctx[++j] = ((-4 * WW + 15 * W + 10 * (NE * 3 - NNE * 3 + NNNE) - (NEEE * 3 - NNEEE * 3 + NNNEEE)) / 20);
	ctx[++j] = ((-3 * WW + 8 * W + clamp4_64(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);
	ctx[++j] = ((W + (NE * 2 - NNE)) / 2 + p1 - (Wp1 + (NEp1 * 2 - NNEp1)) / 2);
	ctx[++j] = ((W + (NE * 2 - NNE)) / 2 + p2 - (Wp2 + (NEp2 * 2 - NNEp2)) / 2);
	ctx[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p1 -(-3 * WWp1 + 8 * Wp1 + (NEEp1 * 2 - NNEEp1)) / 6);
	ctx[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p2 -(-3 * WWp2 + 8 * Wp2 + (NEEp2 * 2 - NNEEp2)) / 6);
	ctx[++j] = ((W + NEE) / 2 + p1 - (Wp1 + NEEp1) / 2);
	ctx[++j] = ((W + NEE) / 2 + p2 - (Wp2 + NEEp2) / 2);
	ctx[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p1 - (WWp1 + (NEEp1 * 2 - NNEEp1)) / 2);
	ctx[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p2 - (WWp2 + (NEEp2 * 2 - NNEEp2)) / 2);
	ctx[++j] = (WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1));
	ctx[++j] = (WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2));
	ctx[++j] = (W + N - NW);
	ctx[++j] = (W + N - NW + p1 - (Wp1 + Np1 - NWp1));
	ctx[++j] = (W + N - NW + p2 - (Wp2 + Np2 - NWp2));
	ctx[++j] = (W + NE - N);
	ctx[++j] = (N + NW - NNW);
	ctx[++j] = (N + NW - NNW + p1 - (Np1 + NWp1 - NNEp1));
	ctx[++j] = (N + NW - NNW + p2 - (Np2 + NWp2 - NNEp2));
	ctx[++j] = (N + NE - NNE);
	ctx[++j] = (N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1));
	ctx[++j] = (N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2));
	ctx[++j] = (N + NN - NNN);
	ctx[++j] = (N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1));
	ctx[++j] = (N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2));
	ctx[++j] = (W + WW - WWW);
	ctx[++j] = (W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1));
	ctx[++j] = (W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2));
	ctx[++j] = (W + NEE - NE);
	ctx[++j] = (W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1));
	ctx[++j] = (W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2));
	ctx[++j] = (NN + p1 - NNp1);
	ctx[++j] = (NN + p2 - NNp2);
	ctx[++j] = (NN + W - NNW);
	ctx[++j] = (NN + W - NNW + p1 - (NNp1 + Wp1 - NNEp1));
	ctx[++j] = (NN + W - NNW + p2 - (NNp2 + Wp2 - NNEp2));
	ctx[++j] = (NN + NW - NNNW);
	ctx[++j] = (NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1));
	ctx[++j] = (NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2));
	ctx[++j] = (NN + NE - NNNE);
	ctx[++j] = (NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1));
	ctx[++j] = (NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2));
	ctx[++j] = (NN + NNNN - NNNNNN);
	ctx[++j] = (NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1));
	ctx[++j] = (NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2));
	ctx[++j] = (WW + p1 - WWp1);
	ctx[++j] = (WW + p2 - WWp2);
	ctx[++j] = (WW + WWWW - WWWWWW);
	ctx[++j] = (WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1));
	ctx[++j] = (WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2));
	ctx[++j] = (N * 2 - NN + p1 - (Np1 * 2 - NNp1));
	ctx[++j] = (N * 2 - NN + p2 - (Np2 * 2 - NNp2));
	ctx[++j] = (W * 2 - WW + p1 - (Wp1 * 2 - WWp1));
	ctx[++j] = (W * 2 - WW + p2 - (Wp2 * 2 - WWp2));
	ctx[++j] = (N * 3 - NN * 3 + NNN);
	ctx[++j] = clamp4_64(N * 3 - NN * 3 + NNN, W, NW, N, NE);
	ctx[++j] = clamp4_64(W * 3 - WW * 3 + WWW, W, NW, N, NE);
	ctx[++j] = clamp4_64(N * 2 - NN, W, NW, N, NE);
	ctx[++j] = ((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4_64(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW, W, NW, N, NN)) / 6);
	ctx[++j] = ((NNNEEE - 4 * NNEE + 6 * NE + (W * 4 - NW * 6 + NNW * 4 - NNNW)) / 4);
	ctx[++j] = (((N + 3 * NW) / 4) * 3 - ((NNW + NNWW) / 2) * 3 + (NNNWW * 3 + NNNWWW) / 4);
	ctx[++j] = ((W * 2 + NW) - (WW + 2 * NWW) + NWWW);
	ctx[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;
	ctx[++j] = (N + W + 1) >> 1;
	ctx[++j] = (NEEEE + NEEEEEE + 1) >> 1;
	ctx[++j] = (WWWWWW + WWWW + 1) >> 1;
	ctx[++j] = ((W + N) * 3 - NW * 2) >> 2;
	ctx[++j] = N;
	ctx[++j] = NN;
	ctx[++j] = N + p1 - Np1;
	ctx[++j] = N + p2 - Np2;
	ctx[++j] = W + p1 - Wp1;
	ctx[++j] = W + p2 - Wp2;
	ctx[++j] = NW + p1 - NWp1;
	ctx[++j] = NW + p2 - NWp2;
	ctx[++j] = NE + p1 - NEp1;
	ctx[++j] = NE + p2 - NEp2;
	ctx[++j] = NN + p1 - NNp1;
	ctx[++j] = NN + p2 - NNp2;
	ctx[++j] = WW + p1 - WWp1;
	ctx[++j] = WW + p2 - WWp2;
	ctx[++j] = W + N - NW;
	ctx[++j] = W + N - NW + p1 - Wp1 - Np1 + NWp1;
	ctx[++j] = W + N - NW + p2 - Wp2 - Np2 + NWp2;
	ctx[++j] = W + NE - N;
	ctx[++j] = W + NE - N + p1 - Wp1 - NEp1 + Np1;
	ctx[++j] = W + NE - N + p2 - Wp2 - NEp2 + Np2;
	ctx[++j] = W + NEE - NE;
	ctx[++j] = W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1;
	ctx[++j] = W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2;
	ctx[++j] = N + NN - NNN;
	ctx[++j] = N + NN - NNN + p1 - Np1 - NNp1 + NNNp1;
	ctx[++j] = N + NN - NNN + p2 - Np2 - NNp2 + NNNp2;
	ctx[++j] = N + NE - NNE;
	ctx[++j] = N + NE - NNE + p1 - Np1 - NEp1 + NNEp1;
	ctx[++j] = N + NE - NNE + p2 - Np2 - NEp2 + NNEp2;
	ctx[++j] = N + NW - NNW;
	ctx[++j] = N + NW - NNW + p1 - Np1 - NWp1 + NNEp1;
	ctx[++j] = N + NW - NNW + p2 - Np2 - NWp2 + NNEp2;
	ctx[++j] = NE + NW - NN;
	ctx[++j] = NE + NW - NN + p1 - NEp1 - NWp1 + NNp1;
	ctx[++j] = NE + NW - NN + p2 - NEp2 - NWp2 + NNp2;
	ctx[++j] = NW + W - NWW;
	ctx[++j] = NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1;
	ctx[++j] = NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2;
	ctx[++j] = W * 2 - WW;
	ctx[++j] = W * 2 - WW + p1 - Wp1 * 2 + WWp1;
	ctx[++j] = W * 2 - WW + p2 - Wp2 * 2 + WWp2;
	ctx[++j] = N * 2 - NN;
	ctx[++j] = N * 2 - NN + p1 - Np1 * 2 + NNp1;
	ctx[++j] = N * 2 - NN + p2 - Np2 * 2 + NNp2;
	ctx[++j] = NW * 2 - NNWW;
	ctx[++j] = NW * 2 - NNWW + p1 - NWp1 * 2 + NNWWp1;
	ctx[++j] = NW * 2 - NNWW + p2 - NWp2 * 2 + NNWWp2;
	ctx[++j] = NE * 2 - NNEE;
	ctx[++j] = NE * 2 - NNEE + p1 - NEp1 * 2 + NNEEp1;
	ctx[++j] = NE * 2 - NNEE + p2 - NEp2 * 2 + NNEEp2;
	ctx[++j] = N * 3 - NN * 3 + NNN + p1 - Np1 * 3 + NNp1 * 3 - NNNp1;
	ctx[++j] = N * 3 - NN * 3 + NNN + p2 - Np2 * 3 + NNp2 * 3 - NNNp2;
	ctx[++j] = N * 3 - NN * 3 + NNN;
	ctx[++j] = (W + NE * 2 - NNE + 1) >> 1;
	ctx[++j] = (W + NE * 3 - NNE * 3 + NNNE+1) >> 1;
	ctx[++j] = (W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1 * 2 - NNEp1) / 2;
	ctx[++j] = (W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2 * 2 - NNEp2) / 2;
	ctx[++j] = NNE + NE - NNNE;
	ctx[++j] = NNE + W - NN;
	ctx[++j] = NNW + W - NNWW;
}
#endif


//T35: Combines spatial transform with entropy coding		THIS IS THE RECORD HOLDER, DO NOT MODIFY

#define T35_N_REC_ESTIMATORS 6
//#define T35_NMAPS (1+134)
//#define T35_NMAPS 134
#define T35_NMAPS 14
//#define T35_NMAPS 3
#define T35_NESTIMATORS (T35_N_REC_ESTIMATORS+T35_NMAPS)
//#define T35_NESTIMATORS T35_NMAPS
#define T35_PRINT_ESTIMATOR_CR

//#define T35_NF0 134
//#define T35_NF1 7
typedef struct T35PredCtxNodeStruct//starts exactly like T35CtxNode
{
	int key;
	int n[2];
} T35PredCtxNode;
typedef struct T35CtxNodeStruct
{
	int key;
	int n[2];
	unsigned short rec[T35_N_REC_ESTIMATORS];
} T35CtxNode;
static CmpRes t35_cmp_node(const void *key, const void *candidate)
{
	int const *k=(int const*)key;
	T35CtxNode const *c=(T35CtxNode const*)candidate;
	return (*k>c->key)-(*k<c->key);
}
static void t35_debugprinter(RBNodeHandle *node, int depth)
{
	T35CtxNode *p;
	if(node)
	{
		p=(T35CtxNode*)node[0]->data;
		printf("0x%08X %5d %5d\n", p->key, p->n[0], p->n[1]);
	}
}
typedef struct T35CtxStruct
{
	Map maps[24][T35_NMAPS];
	int weights[24][T35_NESTIMATORS];
	int found[T35_NMAPS];
	T35CtxNode *node1;
	//T35PredCtxNode *node2[T35_NMAPS];
	T35PredCtxNode *node2[T35_NMAPS-1];
	int context[T35_NMAPS];
	int p0arr[T35_NESTIMATORS];
	int p0_0, p0;//p0_0 isn't clamped
	long long wsum;
	int nnodes[2];//2 types of nodes
#ifdef T35_PRINT_ESTIMATOR_CR
	float csizes[24*T35_NESTIMATORS];
#endif
} T35Ctx;
void t35_ctx_init(T35Ctx *ctx)
{
	for(int k=0;k<24;++k)//fixed 15.16 bit
		for(int k2=0;k2<T35_NESTIMATORS;++k2)
			ctx->weights[k][k2]=0x8000;
	for(int k=0;k<24;++k)
	{
		//for(int k2=0;k2<T35_NMAPS;++k2)
		//	MAP_INIT(ctx->maps[k]+k2, T35PredCtxNode, t35_cmp_node, 0);
#if 1
		MAP_INIT(ctx->maps[k], T35CtxNode, t35_cmp_node, 0);
		for(int k2=1;k2<T35_NMAPS;++k2)
			MAP_INIT(ctx->maps[k]+k2, T35PredCtxNode, t35_cmp_node, 0);
#endif
	}
}
int clip(int x)
{
	x=CLAMP(-128, x, 127);
	return x;
}
static int clamp4(int x, int a, int b, int c, int d)
{
	int vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmax<b)vmax=b;
	if(vmin>c)vmin=c;
	if(vmax<c)vmax=c;
	if(vmin>d)vmin=d;
	if(vmax<d)vmax=d;
	x=CLAMP(vmin, x, vmax);
	return x;
}
void t35_ctx_get_context(T35Ctx *ctx, char *buf, int iw, int ih, int kc, int kx, int ky)
{
#if 1
	int j=-1;

	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	const int offset=0;
	int W   =kx-1>=0         ?buf[(iw* ky   +kx-1)<<2| kc   ]+offset:0,
		NW  =kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]+offset:0,
		N   =ky-1>=0         ?buf[(iw*(ky-1)+kx  )<<2| kc   ]+offset:0,
		NE  =kx+1<iw&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]+offset:0,
		NN  =ky-2>=0         ?buf[(iw*(ky-2)+kx  )<<2| kc   ]+offset:0,

		m1  =kc-1>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-1)]+offset:0,
		Nm1 =kc-1>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-1)]+offset:0,
		Wm1 =kc-1>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-1)]+offset:0,
		NWm1=kc-1>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-1)]+offset:0,

		m2  =kc-2>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-2)]+offset:0,
		Nm2 =kc-2>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-2)]+offset:0,
		Wm2 =kc-2>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-2)]+offset:0,
		NWm2=kc-2>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-2)]+offset:0;
	//int ctxcount3=0;
	//int NW  =kx-1>=0&&ky-1>=0?(++ctxcount3, (char)buf[(iw*(ky-1)+kx-1)<<2| kc   ]):0,
	//	NE  =kx+1<iw&&ky-1>=0?(++ctxcount3, (char)buf[(iw*(ky-1)+kx-1)<<2| kc   ]):0,
	//	Nm1 =kc-1>=0&&ky-1>=0?(++ctxcount3, (char)buf[(iw*(ky-1)+kx  )<<2|(kc-1)]):0,
	//	Wm1 =kc-1>=0&&kx-1>=0?(++ctxcount3, (char)buf[(iw* ky   +kx-1)<<2|(kc-1)]):0,
	//	NWm1=kc-1>=0&&kx-1>=0&&ky-1>=0?(++ctxcount3, (char)buf[(iw*(ky-1)+kx-1)<<2|(kc-1)]):0;
	//ctxcount3+=ctxcount<<1;
	//int helper3=ctxcount3?(((W+N+m1)<<1)+NW+Nm1+Wm1+NWm1)/ctxcount3:0;
	
	ctx->context[++j]=0;
	ctx->context[++j]=N;
	ctx->context[++j]=W;
	ctx->context[++j]=NW;
	ctx->context[++j]=m1;
	ctx->context[++j]=W+NE-N;
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);
	ctx->context[++j]=NW+NE-NN;
	ctx->context[++j]=(N+W-NW + m1)>>1;
	ctx->context[++j]=m2;
	ctx->context[++j]=(N+W-NW + m2)>>1;

	//ctx->context[++j]=Nm1+Wm1-NWm1;
	//ctx->context[++j]=Nm1;
	//ctx->context[++j]=Wm1;
	//ctx->context[++j]=NWm1;
#endif

#if 0
	//offsets are in NW direction
#define LOAD(CO, XO, YO) (unsigned)(kc-CO)<3u&&(unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?buf[(iw*(ky-YO)+kx-(XO))<<2|(kc-CO)]<<T35_PREC:0
	int WWWWWW  =LOAD(0,  6, 0),
		WWWWW   =LOAD(0,  5, 0),
		WWWW    =LOAD(0,  4, 0),
		WWW     =LOAD(0,  3, 0),
		WW      =LOAD(0,  2, 0),
		W       =LOAD(0,  1, 0),
		NWWWW   =LOAD(0,  4, 1),
		NWWW    =LOAD(0,  3, 1),
		NWW     =LOAD(0,  2, 1),
		NW      =LOAD(0,  1, 1),
		N       =LOAD(0,  0, 1),
		NE      =LOAD(0, -1, 1),
		NEE     =LOAD(0, -2, 1),
		NEEE    =LOAD(0, -3, 1),
		NEEEE   =LOAD(0, -4, 1),
		NEEEEEE =LOAD(0, -6, 1),
		NNWWW   =LOAD(0,  3, 2),
		NNWW    =LOAD(0,  2, 2),
		NNW     =LOAD(0,  1, 2),
		NN      =LOAD(0,  0, 2),
		NNE     =LOAD(0, -1, 2),
		NNEE    =LOAD(0, -2, 2),
		NNEEE   =LOAD(0, -3, 2),
		NNNWWWW =LOAD(0,  4, 3),
		NNNWWW  =LOAD(0,  3, 3),
		NNNWW   =LOAD(0,  2, 3),
		NNNW    =LOAD(0,  1, 3),
		NNN     =LOAD(0,  0, 3),
		NNNE    =LOAD(0, -1, 3),
		NNNEE   =LOAD(0, -2, 3),
		NNNEEE  =LOAD(0, -3, 3),
		NNNNW   =LOAD(0,  1, 4),
		NNNN    =LOAD(0,  0, 4),
		NNNNE   =LOAD(0, -1, 4),
		NNNNN   =LOAD(0,  0, 5),
		NNNNNN  =LOAD(0,  0, 6),
		WWWWWWp1=LOAD(1,  6, 0),
		WWWWp1  =LOAD(1,  4, 0),
		WWWp1   =LOAD(1,  3, 0),
		WWp1    =LOAD(1,  2, 0),
		Wp1     =LOAD(1,  1, 0),
		p1      =LOAD(1,  0, 0),
		NWWp1   =LOAD(1,  2, 1),
		NWp1    =LOAD(1,  1, 1),
		Np1     =LOAD(1,  0, 1),
		NEp1    =LOAD(1, -1, 1),
		NEEp1   =LOAD(1, -2, 1),
		NNWWp1  =LOAD(1,  2, 2),
		NNp1    =LOAD(1,  0, 2),
		NNEp1   =LOAD(1, -1, 2),
		NNEEp1  =LOAD(1, -2, 2),
		NNNWp1  =LOAD(1,  1, 3),
		NNNp1   =LOAD(1,  0, 3),
		NNNEp1  =LOAD(1, -1, 3),
		NNNNp1  =LOAD(1,  0, 4),
		NNNNNNp1=LOAD(1,  0, 6),
		WWWWWWp2=LOAD(2,  6, 0),
		WWWWp2  =LOAD(2,  4, 0),
		WWWp2   =LOAD(2,  3, 0),
		WWp2    =LOAD(2,  2, 0),
		Wp2     =LOAD(2,  1, 0),
		p2      =LOAD(2,  0, 0),
		NWWp2   =LOAD(2,  2, 1),
		NWp2    =LOAD(2,  1, 1),
		Np2     =LOAD(2,  0, 1),
		NEp2    =LOAD(2, -1, 1),
		NEEp2   =LOAD(2, -2, 1),
		NNWWp2  =LOAD(2,  2, 2),
		NNp2    =LOAD(2,  0, 2),
		NNEp2   =LOAD(2, -1, 2),
		NNEEp2  =LOAD(2, -2, 2),
		NNNWp2  =LOAD(2,  1, 3),
		NNNp2   =LOAD(2,  0, 3),
		NNNEp2  =LOAD(2, -1, 3),
		NNNNp2  =LOAD(2,  0, 4),
		NNNNNNp2=LOAD(2,  0, 6);
#undef LOAD
	int j=-1;
	ctx->context[++j]=0;

	ctx->context[++j] = clamp4(N + p1 - Np1, W, NW, N, NE);
	ctx->context[++j] = clamp4(N + p2 - Np2, W, NW, N, NE);
	ctx->context[++j] = (W + clamp4(NE * 3 - NNE * 3 + NNNE, W, N, NE, NEE)) / 2;
	ctx->context[++j] = clamp4((W + clip(NE * 2 - NNE)) / 2, W, NW, N, NE);
	ctx->context[++j] = (W + NEE) / 2;
	ctx->context[++j] = ((WWW - 4 * WW + 6 * W + (NE * 4 - NNE * 6 + NNNE * 4 - NNNNE)) / 4);
	ctx->context[++j] = ((-WWWW + 5 * WWW - 10 * WW + 10 * W + clamp4(NE * 4 - NNE * 6 + NNNE * 4 - NNNNE, N, NE, NEE, NEEE)) / 5);
	ctx->context[++j] = ((-4 * WW + 15 * W + 10 * (NE * 3 - NNE * 3 + NNNE) - (NEEE * 3 - NNEEE * 3 + NNNEEE)) / 20);
	ctx->context[++j] = ((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p1 - (Wp1 + (NEp1 * 2 - NNEp1)) / 2);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p2 - (Wp2 + (NEp2 * 2 - NNEp2)) / 2);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p1 -(-3 * WWp1 + 8 * Wp1 + (NEEp1 * 2 - NNEEp1)) / 6);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p2 -(-3 * WWp2 + 8 * Wp2 + (NEEp2 * 2 - NNEEp2)) / 6);
	ctx->context[++j] = ((W + NEE) / 2 + p1 - (Wp1 + NEEp1) / 2);
	ctx->context[++j] = ((W + NEE) / 2 + p2 - (Wp2 + NEEp2) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p1 - (WWp1 + (NEEp1 * 2 - NNEEp1)) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p2 - (WWp2 + (NEEp2 * 2 - NNEEp2)) / 2);
	ctx->context[++j] = (WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1));
	ctx->context[++j] = (WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2));
	ctx->context[++j] = (W + N - NW);
	ctx->context[++j] = (W + N - NW + p1 - (Wp1 + Np1 - NWp1));
	ctx->context[++j] = (W + N - NW + p2 - (Wp2 + Np2 - NWp2));
	ctx->context[++j] = (W + NE - N);
	ctx->context[++j] = (N + NW - NNW);
	ctx->context[++j] = (N + NW - NNW + p1 - (Np1 + NWp1 - NNEp1));
	ctx->context[++j] = (N + NW - NNW + p2 - (Np2 + NWp2 - NNEp2));
	ctx->context[++j] = (N + NE - NNE);
	ctx->context[++j] = (N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1));
	ctx->context[++j] = (N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2));
	ctx->context[++j] = (N + NN - NNN);
	ctx->context[++j] = (N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1));
	ctx->context[++j] = (N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2));
	ctx->context[++j] = (W + WW - WWW);
	ctx->context[++j] = (W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1));
	ctx->context[++j] = (W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2));
	ctx->context[++j] = (W + NEE - NE);
	ctx->context[++j] = (W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1));
	ctx->context[++j] = (W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2));
	ctx->context[++j] = (NN + p1 - NNp1);
	ctx->context[++j] = (NN + p2 - NNp2);
	ctx->context[++j] = (NN + W - NNW);
	ctx->context[++j] = (NN + W - NNW + p1 - (NNp1 + Wp1 - NNEp1));
	ctx->context[++j] = (NN + W - NNW + p2 - (NNp2 + Wp2 - NNEp2));
	ctx->context[++j] = (NN + NW - NNNW);
	ctx->context[++j] = (NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1));
	ctx->context[++j] = (NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2));
	ctx->context[++j] = (NN + NE - NNNE);
	ctx->context[++j] = (NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1));
	ctx->context[++j] = (NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2));
	ctx->context[++j] = (NN + NNNN - NNNNNN);
	ctx->context[++j] = (NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1));
	ctx->context[++j] = (NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2));
	ctx->context[++j] = (WW + p1 - WWp1);
	ctx->context[++j] = (WW + p2 - WWp2);
	ctx->context[++j] = (WW + WWWW - WWWWWW);
	ctx->context[++j] = (WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1));
	ctx->context[++j] = (WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2));
	ctx->context[++j] = (N * 2 - NN + p1 - (Np1 * 2 - NNp1));
	ctx->context[++j] = (N * 2 - NN + p2 - (Np2 * 2 - NNp2));
	ctx->context[++j] = (W * 2 - WW + p1 - (Wp1 * 2 - WWp1));
	ctx->context[++j] = (W * 2 - WW + p2 - (Wp2 * 2 - WWp2));
	ctx->context[++j] = (N * 3 - NN * 3 + NNN);
	ctx->context[++j] = clamp4(N * 3 - NN * 3 + NNN, W, NW, N, NE);
	ctx->context[++j] = clamp4(W * 3 - WW * 3 + WWW, W, NW, N, NE);
	ctx->context[++j] = clamp4(N * 2 - NN, W, NW, N, NE);
	ctx->context[++j] = ((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW, W, NW, N, NN)) / 6);
	ctx->context[++j] = ((NNNEEE - 4 * NNEE + 6 * NE + (W * 4 - NW * 6 + NNW * 4 - NNNW)) / 4);
	ctx->context[++j] = (((N + 3 * NW) / 4) * 3 - ((NNW + NNWW) / 2) * 3 + (NNNWW * 3 + NNNWWW) / 4);
	ctx->context[++j] = ((W * 2 + NW) - (WW + 2 * NWW) + NWWW);
	ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;
	ctx->context[++j] = (N + W + 1) >> 1;
	ctx->context[++j] = (NEEEE + NEEEEEE + 1) >> 1;
	ctx->context[++j] = (WWWWWW + WWWW + 1) >> 1;
	ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;
	ctx->context[++j] = N;
	ctx->context[++j] = NN;
	ctx->context[++j] = N + p1 - Np1;
	ctx->context[++j] = N + p2 - Np2;
	ctx->context[++j] = W + p1 - Wp1;
	ctx->context[++j] = W + p2 - Wp2;
	ctx->context[++j] = NW + p1 - NWp1;
	ctx->context[++j] = NW + p2 - NWp2;
	ctx->context[++j] = NE + p1 - NEp1;
	ctx->context[++j] = NE + p2 - NEp2;
	ctx->context[++j] = NN + p1 - NNp1;
	ctx->context[++j] = NN + p2 - NNp2;
	ctx->context[++j] = WW + p1 - WWp1;
	ctx->context[++j] = WW + p2 - WWp2;
	ctx->context[++j] = W + N - NW;
	ctx->context[++j] = W + N - NW + p1 - Wp1 - Np1 + NWp1;
	ctx->context[++j] = W + N - NW + p2 - Wp2 - Np2 + NWp2;
	ctx->context[++j] = W + NE - N;
	ctx->context[++j] = W + NE - N + p1 - Wp1 - NEp1 + Np1;
	ctx->context[++j] = W + NE - N + p2 - Wp2 - NEp2 + Np2;
	ctx->context[++j] = W + NEE - NE;
	ctx->context[++j] = W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1;
	ctx->context[++j] = W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2;
	ctx->context[++j] = N + NN - NNN;
	ctx->context[++j] = N + NN - NNN + p1 - Np1 - NNp1 + NNNp1;
	ctx->context[++j] = N + NN - NNN + p2 - Np2 - NNp2 + NNNp2;
	ctx->context[++j] = N + NE - NNE;
	ctx->context[++j] = N + NE - NNE + p1 - Np1 - NEp1 + NNEp1;
	ctx->context[++j] = N + NE - NNE + p2 - Np2 - NEp2 + NNEp2;
	ctx->context[++j] = N + NW - NNW;
	ctx->context[++j] = N + NW - NNW + p1 - Np1 - NWp1 + NNEp1;
	ctx->context[++j] = N + NW - NNW + p2 - Np2 - NWp2 + NNEp2;
	ctx->context[++j] = NE + NW - NN;
	ctx->context[++j] = NE + NW - NN + p1 - NEp1 - NWp1 + NNp1;
	ctx->context[++j] = NE + NW - NN + p2 - NEp2 - NWp2 + NNp2;
	ctx->context[++j] = NW + W - NWW;
	ctx->context[++j] = NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1;
	ctx->context[++j] = NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2;
	ctx->context[++j] = W * 2 - WW;
	ctx->context[++j] = W * 2 - WW + p1 - Wp1 * 2 + WWp1;
	ctx->context[++j] = W * 2 - WW + p2 - Wp2 * 2 + WWp2;
	ctx->context[++j] = N * 2 - NN;
	ctx->context[++j] = N * 2 - NN + p1 - Np1 * 2 + NNp1;
	ctx->context[++j] = N * 2 - NN + p2 - Np2 * 2 + NNp2;
	ctx->context[++j] = NW * 2 - NNWW;
	ctx->context[++j] = NW * 2 - NNWW + p1 - NWp1 * 2 + NNWWp1;
	ctx->context[++j] = NW * 2 - NNWW + p2 - NWp2 * 2 + NNWWp2;
	ctx->context[++j] = NE * 2 - NNEE;
	ctx->context[++j] = NE * 2 - NNEE + p1 - NEp1 * 2 + NNEEp1;
	ctx->context[++j] = NE * 2 - NNEE + p2 - NEp2 * 2 + NNEEp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p1 - Np1 * 3 + NNp1 * 3 - NNNp1;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p2 - Np2 * 3 + NNp2 * 3 - NNNp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN;
	ctx->context[++j] = (W + NE * 2 - NNE + 1) >> 1;
	ctx->context[++j] = (W + NE * 3 - NNE * 3 + NNNE+1) >> 1;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1 * 2 - NNEp1) / 2;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2 * 2 - NNEp2) / 2;
	ctx->context[++j] = NNE + NE - NNNE;
	ctx->context[++j] = NNE + W - NN;
	ctx->context[++j] = NNW + W - NNWW;
#endif
	for(int k=0;k<T35_NMAPS;++k)
	{
		ctx->context[k]-=offset;
		ctx->context[k]<<=8;
	}
}
void t35_ctx_estimate_p0(T35Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	RBNodeHandle *hnode;
#if 1
	hnode=map_insert(ctx->maps[workidx], ctx->context, ctx->found);
	ctx->node1=(T35CtxNode*)hnode[0]->data;
	if(ctx->found[0])
	{
		int k=0;
		sum=ctx->node1->n[0]+ctx->node1->n[1];
		//if(!sum)
		//	LOG_ERROR("");
		ctx->p0arr[p0idx+k]=sum?(int)(((long long)ctx->node1->n[0]<<16)/sum):0x8000;
		++k;
		for(;k<T35_N_REC_ESTIMATORS+1;++k)
			ctx->p0arr[p0idx+k]=ctx->node1->rec[k-1];
		p0idx+=k;
	}
	else
	{
		int k=0;
		for(;k<T35_N_REC_ESTIMATORS+1;++k)
			ctx->p0arr[p0idx+k]=0x8000;
		p0idx+=k;
	}
#endif

	T35PredCtxNode *node;
	for(int k=1;k<T35_NMAPS;++k)
	{
		hnode=map_insert(ctx->maps[workidx]+k, ctx->context+k, ctx->found+k);
		node=ctx->node2[k-1]=(T35PredCtxNode*)hnode[0]->data;
		if(ctx->found[k])
		{
			sum=node->n[0]+node->n[1];
			//if(!sum)
			//	LOG_ERROR("");
			ctx->p0arr[p0idx]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
			++p0idx;
		}
		else
		{
			ctx->p0arr[p0idx]=0x8000;
			++p0idx;
		}
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T35_NESTIMATORS;++k)
	{
		sum+=(long long)ctx->p0arr[k]*wk[k];
		ctx->wsum+=wk[k];
	}
	//ctx->p0=ctx->wsum?(int)((sum+(ctx->wsum>>1))/ctx->wsum):0x8000;//same CR
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t35_ctx_update(T35Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
#ifdef T35_PRINT_ESTIMATOR_CR
	for(int k=0;k<T35_NESTIMATORS;++k)
	{
		int prob=(unsigned short)(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes[T35_NESTIMATORS*workidx+k]+=bitsize;
		}
	}
#endif
	//bwd
	int *wk=ctx->weights[workidx];
	if(ctx->p0_0>=1&&ctx->p0_0<=0xFFFF)
	{
		int p_bit=bit?0x10000-ctx->p0:ctx->p0;
		long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
		dL_dp0^=-bit;
		dL_dp0+=bit;
		for(int k=0;k<T35_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
#if 1
	{
		T35CtxNode *node=ctx->node1;
		if(ctx->found[0])
		{
			++node->n[bit];
			for(int k=0;k<T35_N_REC_ESTIMATORS;++k)
			{
				int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>((k+1)<<1));
				node->rec[k]=CLAMP(1, temp, 0xFFFF);
			}
		}
		else
		{
			node->key=ctx->context[0];
			node->n[0]=1;
			node->n[1]=1;
			for(int k=0;k<T35_N_REC_ESTIMATORS;++k)
				node->rec[k]=0x8000;
			++ctx->nnodes[0];
		}
		ctx->context[0]|=bit<<kb;
	}
#endif

	for(int k=1;k<T35_NMAPS;++k)
	{
		T35PredCtxNode *node=ctx->node2[k-1];
		if(ctx->found[k])
			++node->n[bit];
		else
		{
			node->key=ctx->context[k];
			node->n[0]=1;
			node->n[1]=1;
			++ctx->nnodes[1];
		}
		ctx->context[k]|=bit<<kb;
	}
}
void t35_ctx_clear(T35Ctx *ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T35_NMAPS;++k2)
			MAP_CLEAR(ctx->maps[k]+k2);
	}
	ctx->nnodes[0]=0;
	ctx->nnodes[1]=0;
}
#if 0
static const int t35_bitidx[]=
{//ch, bit
	0, 7,
	1, 7,
	2, 7,
	2, 6,
	1, 6,
	0, 6,
	0, 5,
	1, 5,
	2, 5,
	2, 4,
	1, 4,
	0, 4,
	0, 3,
	1, 3,
	2, 3,
	2, 2,
	1, 2,
	0, 2,
	0, 1,
	1, 1,
	2, 1,
	2, 0,
	1, 0,
	0, 0,
};
static void initialize_fix16(DataType *weights, int count, int fan_in)
{
	DataType gain=(int)(0x10000/fan_in);
	for(int k=0;k<count;++k)
	{
		DataType x=(DataType)(xoroshiro128_next()&0xFFFF);
		weights[k]=x*gain>>16;//[0, 2^16-1]/sqrt_fan_in
	}
}
static void initialize_fp32(float *weights, int count, int fan_in)
{
	float gain=1/(0x10000*sqrtf((float)fan_in));
	for(int k=0;k<count;++k)
	{
		int x=(int)(xoroshiro128_next()&0xFFFF);
		weights[k]=x*gain;//[0, 1]/sqrt_fan_in
	}
}
#endif
T35Ctx t35_context;
int t35_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//colortransform_ycocb_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);

	apply_transforms_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 192);//makes MSB easy

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float csizes[24]={0};
	int hits[24]={0};
	
	t35_ctx_init(&t35_context);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				t35_ctx_get_context(&t35_context, (char*)buf2, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t35_ctx_estimate_p0(&t35_context, kc, kb);

					int bit=buf2[(iw*ky+kx)<<2|kc]>>kb&1;
					ac_enc_bin(&ec, t35_context.p0, bit);
					
					int hit=bit?t35_context.p0<0x8000:t35_context.p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					int prob=bit?0x10000-t35_context.p0:t35_context.p0;
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t35_ctx_update(&t35_context, kc, kb, bit);
				}
			}
		}
		if(loud==2)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("Y %d  CR %f\n", ky, iw*3/(csize-csize_prev));
			csize_prev=csize;
		}
	}
	ac_enc_flush(&ec);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", ((float)t35_context.nnodes[0]*sizeof(T35CtxNode)+(float)t35_context.nnodes[1]*sizeof(T35PredCtxNode))/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");

		for(int kc=0;kc<3;++kc)
		{
			printf("C%d  ", kc);
			for(int kb=7;kb>=0;--kb)
				printf("  B%d %11.3f", kb, iw*ih/csizes[kc<<3|kb]);
			printf("\n");
		}
		printf("\n");

#ifdef T35_PRINT_ESTIMATOR_CR
		printf("Estimator efficiencies:\n");
		int minidx[24]={0}, maxidx[24]={0};
		for(int kb=0;kb<24;++kb)
		{
			float *sizes=t35_context.csizes+T35_NESTIMATORS*kb;
			for(int ke=1;ke<T35_NESTIMATORS;++ke)
			{
				if(sizes[minidx[kb]]>sizes[ke])
					minidx[kb]=ke;
				if(sizes[maxidx[kb]]<sizes[ke])
					maxidx[kb]=ke;
			}
		}
		for(int ke=0;ke<T35_NESTIMATORS;++ke)
		{
			float *sizes=t35_context.csizes+ke;
			printf("E%2d ", ke);
			for(int kb=0;kb<24;++kb)
			{
				char c;
				if(ke==minidx[kb])
					c='*';
				else if(ke==maxidx[kb])
					c='L';
				else
					c=' ';
				printf(" %7.2f%c", sizes[T35_NESTIMATORS*kb]/t35_context.csizes[T35_NESTIMATORS*kb+minidx[kb]], c);
			}
			printf("\n");
		}
#endif
	}
	t35_ctx_clear(&t35_context);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t35_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	//int debug_index=0;

	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t35_ctx_init(&t35_context);

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				t35_ctx_get_context(&t35_context, (char*)buf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t35_ctx_estimate_p0(&t35_context, kc, kb);
					
					int bit=ac_dec_bin(&ec, t35_context.p0);
					buf[(iw*ky+kx)<<2|kc]|=bit<<kb;

					t35_ctx_update(&t35_context, kc, kb, bit);
				}
			}
		}
	}
	t35_ctx_clear(&t35_context);
	
	//addbuf(buf, iw, ih, 3, 4, 128);
	apply_transforms_inv(buf, iw, ih);
	if(loud)
	{
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}


//T36: stretch & squish from paq8px
#if 0
#define T36_NESTIMATORS 25//135
#define T36_LR (long long)(0.001*0x10000+0.5)
#define T36_ACT_ROUGHNESS 9	//do not change
#define T36_SQRT_ONEPLUS_ROUGHNESS 0x3298B
#define T36_PRINT_ESTIMATOR_CR
int sqrt_f16(int v)//https://github.com/chmike/fpsqrt
{
	unsigned t, q, b, r;

	r=v;
	b=0x40000000;
	q=0;
	while(b>0x40)
	{
		t=q+b;
		if(r>=t)
		{
			r-=t;
			q=t+b;//equivalent to q += 2*b
		}
		r<<=1;
		b>>=1;
	}
	q>>=8;
	return q;
}
int t36_stretch(int p0)
{
	//stretch(x) = x/resqrt(1 + sq3*(1-4xx))	[-1/2, 1/2]
	int p2=(int)((long long)p0*p0>>16);
	p2<<=2;
	p2=0x10000-p2;
	p2*=T36_ACT_ROUGHNESS;
	p2+=0x10000;
	p2=sqrt_f16(p2);
	p2=(int)(((long long)p0<<16)/p2);
	return p2;
}
int t36_squish(int x)
{
	//squish(x) = resqrt(1+sq3)*x/resqrt(1+sq(2*3)xx)
	int x2=(int)((long long)x*x>>16);
	x2*=T36_ACT_ROUGHNESS<<2;
	x2+=0x10000;
	x2=sqrt_f16(x2);
	x2=(int)((long long)T36_SQRT_ONEPLUS_ROUGHNESS*x/x2);
	return x2;
}
int t36_squish_dash(int x, int squishx)
{
	//squish'(x) = x ? (1 - sq(2*3)/(1+sq3)*sq(squishx))*(squishx)/x : sqrt(1+sq3)
	int temp;
	if(x)
	{
		temp=(int)((long long)squishx*squishx>>16);
		temp*=T36_ACT_ROUGHNESS<<2;
		temp/=T36_ACT_ROUGHNESS+1;
		temp=0x10000-temp;
		temp=(int)((long long)temp*squishx/x);
	}
	else
		temp=T36_SQRT_ONEPLUS_ROUGHNESS;
	return temp;
}
void t36_stretch_squish_test()
{
	printf("ramp,  stretch, ramp,  squish, ramp\n");
	for(int k=0;k<256;++k)
	{
		int v0=(k<<8)-0x8000;
		int v1=t36_stretch(v0);
		int v2=t36_squish(v1);
		int v3=t36_squish(v0);
		int v4=t36_stretch(v3);
		printf("%12d %12d %12d %12d %12d\n", v0, v1, v2, v3, v4);
		//printf("0x%08X 0x%08X 0x%08X\n", v0, v1, v2);
	}
	pause();
}
typedef struct T36NodeStruct
{
	int key;
	int n[2];
} T36Node;
static CmpRes t36_cmp_node(const void *key, const void *candidate)
{
	int const *k=(int const*)key;
	T36Node const *c=(T36Node const*)candidate;
	return (*k>c->key)-(*k<c->key);
}
//static void t36_debugprinter(RBNodeHandle *node, int depth)
//{
//	T36Node *p;
//	if(node)
//	{
//		p=(T36Node*)node[0]->data;
//		printf("0x%08X %5d %5d\n", p->key, p->n[0], p->n[1]);
//	}
//}
typedef struct T36CtxStruct
{
	Map maps[24][T36_NESTIMATORS];
	int context[T36_NESTIMATORS];
	int found[T36_NESTIMATORS];
	T36Node *node[T36_NESTIMATORS];
	int weights[24][T36_NESTIMATORS];//fixed 15.16 bit
	int t[T36_NESTIMATORS], x, s0, p0;//t[i]=stretch(p0arr[i]),  x = w . t + w0,  s0=squish(x) in [-1/2, 1/2],  p0=CLAMP(s0+1/2)
	long long wsum;
	int nnodes;
	int nestimators;
	int p0arr[T36_NESTIMATORS];
#ifdef T36_PRINT_ESTIMATOR_CR
	float csizes[24*T36_NESTIMATORS];
#endif

	int iw;
	int *ebufs[24*T36_NESTIMATORS];//probability errors
	int ei[T36_NESTIMATORS];

	char *ebuf;//pixel errors
} T36Ctx;
void t36_ctx_init(T36Ctx *ctx, int iw, int ih)
{
	XOROSHIRO128_RESET();
	for(int k=0;k<24;++k)//fixed 15.16 bit
		for(int k2=0;k2<T36_NESTIMATORS;++k2)
			ctx->weights[k][k2]=0x8000+(xoroshiro128_next()&0x1FF)-0x100;
			//ctx->weights[k][k2]=(xoroshiro128_next()&0xFFFF)/T36_NESTIMATORS;
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T36_NESTIMATORS;++k2)
			MAP_INIT(ctx->maps[k]+k2, T36Node, t36_cmp_node, 0);
	}

	ctx->iw=iw;
	int initval=1;
	for(int k=0;k<24*T36_NESTIMATORS;++k)
	{
		ctx->ebufs[k]=(int*)malloc(iw*sizeof(int));
		if(!ctx->ebufs[k])
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memfill(ctx->ebufs[k], &initval, iw*sizeof(int), sizeof(int));
	}

	ctx->ebuf=(char*)malloc((size_t)iw*ih<<2);
	if(!ctx->ebuf)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	initval=0xFF000000;
	memfill(ctx->ebuf, &initval, (size_t)iw*ih<<2, sizeof(int));
}
void t36_ctx_clear(T36Ctx *ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T36_NESTIMATORS;++k2)
			MAP_CLEAR(ctx->maps[k]+k2);
	}
	ctx->nnodes=0;
	
	for(int k=0;k<24*T36_NESTIMATORS;++k)
	{
		free(ctx->ebufs[k]);
		ctx->ebufs[k]=0;
	}

	free(ctx->ebuf);
	ctx->ebuf=0;
}
void t36_ctx_get_context(T36Ctx *ctx, char *buf, int iw, int ih, int kc, int kx, int ky)
{
#if 1
	int j=-1;
#define LOAD(BUF, CO, XO, YO) (unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?BUF[(iw*(ky-YO)+kx-(XO))<<2|(kc+CO)%3]:0
	char
		NNWW=LOAD(buf, 0,  2, 2),
		NNW =LOAD(buf, 0,  1, 2),
		NN  =LOAD(buf, 0,  0, 2),
		NNE =LOAD(buf, 0, -1, 2),
		NNEE=LOAD(buf, 0, -2, 2),
		
		NWW =LOAD(buf, 0,  2, 1),
		NW  =LOAD(buf, 0,  1, 1),
		N   =LOAD(buf, 0,  0, 1),
		NE  =LOAD(buf, 0, -1, 1),
		NEE =LOAD(buf, 0, -2, 1),

		WW  =LOAD(buf, 0,  2, 0),
		W   =LOAD(buf, 0,  1, 0),

		eNNWW=LOAD(ctx->ebuf, 0,  2, 2),
		eNNW =LOAD(ctx->ebuf, 0,  1, 2),
		eNN  =LOAD(ctx->ebuf, 0,  0, 2),
		eNNE =LOAD(ctx->ebuf, 0, -1, 2),
		eNNEE=LOAD(ctx->ebuf, 0, -2, 2),
		
		eNWW =LOAD(ctx->ebuf, 0,  2, 1),
		eNW  =LOAD(ctx->ebuf, 0,  1, 1),
		eN   =LOAD(ctx->ebuf, 0,  0, 1),
		eNE  =LOAD(ctx->ebuf, 0, -1, 1),
		eNEE =LOAD(ctx->ebuf, 0, -2, 1),

		eWW  =LOAD(ctx->ebuf, 0,  2, 0),
		eW   =LOAD(ctx->ebuf, 0,  1, 0),

		Wp1 =LOAD(buf, 1,  1, 0),
		Np1 =LOAD(buf, 1,  0, 1),
		NWp1=LOAD(buf, 1,  1, 1),
		NEp1=LOAD(buf, 1, -1, 1),

		Wp2 =LOAD(buf, 2,  1, 0),
		Np2 =LOAD(buf, 2,  0, 1),
		NWp2=LOAD(buf, 2,  1, 1),
		NEp2=LOAD(buf, 2, -1, 1);
#undef LOAD
	ctx->context[++j]=0;
#if 1
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);
	ctx->context[++j]=(N+W)>>1;
	ctx->context[++j]=N*2-NN;
	ctx->context[++j]=W*2-WW;
	ctx->context[++j]=NW*2-NNWW;
	ctx->context[++j]=NE+NW-NN;
	ctx->context[++j]=NE*2-NNEE;
	ctx->context[++j]=W+NE-N;
	ctx->context[++j]=N*3-NNW-NE;
	ctx->context[++j]=(W*4+NEE*2-(WW+NE*2))/3;
	ctx->context[++j]=N*3-NNE-NW;
	ctx->context[++j]=W+NW-NWW;
	ctx->context[++j]=N+NE-NNE;
	ctx->context[++j]=N+NW-NNW;
	ctx->context[++j]=N-(NE*2+WW)+(W*2+NEE);
	ctx->context[++j]=(N+NW+W+NN)/2-NNW;
	ctx->context[++j]=(N+NE+W+NNEE)/2-NNE;
	ctx->context[++j]=(N+W)-(NNW+NWW)/2;
	ctx->context[++j]=(NNW+W+NW*2)/2-NNWW;
	ctx->context[++j]=(W+NNE+N*2)/2-NN;
	ctx->context[++j]=(N+WW+NW+W)/2-NWW;
	switch(kc)
	{
	case 0://R
		{
			ctx->context[++j]=(
				 0x00AD*NNWW-0x020A*NNW+0x0259*NN-0x0230*NNE+0x0161*NNEE
				-0x034E*NWW +0x01D7*NW +0x0326*N +0x04B7*NE -0x0091*NEE
				+0x0504*WW  +0x0544*W
				+0x002B*eNNWW+0x0189*eNNW+0x00A6*eNN+0x016D*eNNE-0x006B*eNNEE
				+0x0199*eNWW +0x02DD*eNW +0x05F0*eN +0x0148*eNE +0x0117*eNEE
				+0x0024*eWW  +0x068D*eW
			)>>12;
		}
		break;
	case 1://G
		{
			ctx->context[++j]=(
				 0x0033*NNWW-0x0020*NNW+0x00B3*NN+0x0048*NNE+0x00BD*NNEE
				+0x0065*NWW +0x00CD*NW +0x015D*N -0x000D*NE +0x00AD*NEE
				+0x01B7*WW  +0x0928*W
				+0x0058*eNNWW-0x004E*eNNW-0x0132*eNN+0x0083*eNNE-0x0013*eNNEE
				-0x0163*eNWW +0x01A4*eNW +0x0688*eN +0x01E1*eNE +0x004F*eNEE
				-0x03AD*eWW  +0x0347*eW
			)>>12;
		}
		break;
	case 2://B
		{
			ctx->context[++j]=(
				 0x0018*NNWW-0x021A*NNW+0x034A*NN-0x0249*NNE+0x00F7*NNEE
				-0x0261*NWW -0x0008*NW +0x044F*N +0x0347*NE +0x0018*NEE
				+0x0346*WW  +0x0797*W
				+0x0043*eNNWW+0x00EC*eNNW+0x0041*eNN+0x014C*eNNE-0x0067*eNNEE
				+0x009D*eNWW +0x03B9*eNW +0x0627*eN +0x0214*eNE +0x004F*eNEE
				+0x0026*eWW  +0x0543*eW
			)>>12;
		}
		break;
	}
#endif
	ctx->nestimators=j+1;
	int pred=0;
	for(int k=0;k<ctx->nestimators;++k)
	{
		ctx->context[k]=clip(ctx->context[k]);
		pred+=ctx->context[k];
	}
	pred/=ctx->nestimators;
	ctx->ebuf[(iw*ky+kx)<<2|kc]=pred;
#endif

#if 0
	int j=-1;

	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	int W   =kx-1>=0         ?buf[(iw* ky   +kx-1)<<2| kc   ]:0,
		NW  =kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]:0,
		N   =ky-1>=0         ?buf[(iw*(ky-1)+kx  )<<2| kc   ]:0,
		NE  =kx+1<iw&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]:0,
		NN  =ky-2>=0         ?buf[(iw*(ky-2)+kx  )<<2| kc   ]:0,

		m1  =kc-1>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-1)]:0,
		Nm1 =kc-1>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-1)]:0,
		Wm1 =kc-1>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-1)]:0,
		NWm1=kc-1>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-1)]:0,

		m2  =kc-2>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-2)]:0,
		Nm2 =kc-2>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-2)]:0,
		Wm2 =kc-2>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-2)]:0,
		NWm2=kc-2>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-2)]:0;
	
	//ctx->context[++j]=0;
	ctx->context[++j]=N;
	ctx->context[++j]=W;
	ctx->context[++j]=NW;
	ctx->context[++j]=m1;
	//ctx->context[++j]=W+NE-N;
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);
	//ctx->context[++j]=NW+NE-NN;
	//ctx->context[++j]=(N+W-NW + m1)>>1;
	ctx->context[++j]=m2;
	ctx->context[++j]=(N+W-NW + m2)>>1;
	//ctx->context[++j]=Nm1+Wm1-NWm1;
	//ctx->context[++j]=Nm1;
	//ctx->context[++j]=Wm1;
	//ctx->context[++j]=NWm1;

	ctx->nestimators=j+1;
#endif
	
#if 0
	//offsets are in NW direction
#define LOAD(CO, XO, YO) (unsigned)(kc-CO)<3u&&(unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?buf[(iw*(ky-YO)+kx-(XO))<<2|(kc-CO)]:0
	int WWWWWW  =LOAD(0,  6, 0),
		WWWWW   =LOAD(0,  5, 0),
		WWWW    =LOAD(0,  4, 0),
		WWW     =LOAD(0,  3, 0),
		WW      =LOAD(0,  2, 0),
		W       =LOAD(0,  1, 0),
		NWWWW   =LOAD(0,  4, 1),
		NWWW    =LOAD(0,  3, 1),
		NWW     =LOAD(0,  2, 1),
		NW      =LOAD(0,  1, 1),
		N       =LOAD(0,  0, 1),
		NE      =LOAD(0, -1, 1),
		NEE     =LOAD(0, -2, 1),
		NEEE    =LOAD(0, -3, 1),
		NEEEE   =LOAD(0, -4, 1),
		NEEEEEE =LOAD(0, -6, 1),
		NNWWW   =LOAD(0,  3, 2),
		NNWW    =LOAD(0,  2, 2),
		NNW     =LOAD(0,  1, 2),
		NN      =LOAD(0,  0, 2),
		NNE     =LOAD(0, -1, 2),
		NNEE    =LOAD(0, -2, 2),
		NNEEE   =LOAD(0, -3, 2),
		NNNWWWW =LOAD(0,  4, 3),
		NNNWWW  =LOAD(0,  3, 3),
		NNNWW   =LOAD(0,  2, 3),
		NNNW    =LOAD(0,  1, 3),
		NNN     =LOAD(0,  0, 3),
		NNNE    =LOAD(0, -1, 3),
		NNNEE   =LOAD(0, -2, 3),
		NNNEEE  =LOAD(0, -3, 3),
		NNNNW   =LOAD(0,  1, 4),
		NNNN    =LOAD(0,  0, 4),
		NNNNE   =LOAD(0, -1, 4),
		NNNNN   =LOAD(0,  0, 5),
		NNNNNN  =LOAD(0,  0, 6),
		WWWWWWp1=LOAD(1,  6, 0),
		WWWWp1  =LOAD(1,  4, 0),
		WWWp1   =LOAD(1,  3, 0),
		WWp1    =LOAD(1,  2, 0),
		Wp1     =LOAD(1,  1, 0),
		p1      =LOAD(1,  0, 0),
		NWWp1   =LOAD(1,  2, 1),
		NWp1    =LOAD(1,  1, 1),
		Np1     =LOAD(1,  0, 1),
		NEp1    =LOAD(1, -1, 1),
		NEEp1   =LOAD(1, -2, 1),
		NNWWp1  =LOAD(1,  2, 2),
		NNp1    =LOAD(1,  0, 2),
		NNEp1   =LOAD(1, -1, 2),
		NNEEp1  =LOAD(1, -2, 2),
		NNNWp1  =LOAD(1,  1, 3),
		NNNp1   =LOAD(1,  0, 3),
		NNNEp1  =LOAD(1, -1, 3),
		NNNNp1  =LOAD(1,  0, 4),
		NNNNNNp1=LOAD(1,  0, 6),
		WWWWWWp2=LOAD(2,  6, 0),
		WWWWp2  =LOAD(2,  4, 0),
		WWWp2   =LOAD(2,  3, 0),
		WWp2    =LOAD(2,  2, 0),
		Wp2     =LOAD(2,  1, 0),
		p2      =LOAD(2,  0, 0),
		NWWp2   =LOAD(2,  2, 1),
		NWp2    =LOAD(2,  1, 1),
		Np2     =LOAD(2,  0, 1),
		NEp2    =LOAD(2, -1, 1),
		NEEp2   =LOAD(2, -2, 1),
		NNWWp2  =LOAD(2,  2, 2),
		NNp2    =LOAD(2,  0, 2),
		NNEp2   =LOAD(2, -1, 2),
		NNEEp2  =LOAD(2, -2, 2),
		NNNWp2  =LOAD(2,  1, 3),
		NNNp2   =LOAD(2,  0, 3),
		NNNEp2  =LOAD(2, -1, 3),
		NNNNp2  =LOAD(2,  0, 4),
		NNNNNNp2=LOAD(2,  0, 6);
#undef LOAD
	int j=-1;
	ctx->context[++j]=0;
	ctx->context[++j] = clamp4(N + p1 - Np1, W, NW, N, NE);
	ctx->context[++j] = clamp4(N + p2 - Np2, W, NW, N, NE);
#if 1
	ctx->context[++j] = (W + clamp4(NE * 3 - NNE * 3 + NNNE, W, N, NE, NEE)) / 2;
	ctx->context[++j] = clamp4((W + clip(NE * 2 - NNE)) / 2, W, NW, N, NE);
#endif
	ctx->context[++j] = (W + NEE) / 2;
#if 1
	//ctx->context[++j] = ((WWW - 4 * WW + 6 * W + (NE * 4 - NNE * 6 + NNNE * 4 - NNNNE)) / 4);
	//ctx->context[++j] = ((-WWWW + 5 * WWW - 10 * WW + 10 * W + clamp4(NE * 4 - NNE * 6 + NNNE * 4 - NNNNE, N, NE, NEE, NEEE)) / 5);
	//ctx->context[++j] = ((-4 * WW + 15 * W + 10 * (NE * 3 - NNE * 3 + NNNE) - (NEEE * 3 - NNEEE * 3 + NNNEEE)) / 20);
	//ctx->context[++j] = ((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);
	//ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p1 - (Wp1 + (NEp1 * 2 - NNEp1)) / 2);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p2 - (Wp2 + (NEp2 * 2 - NNEp2)) / 2);
	//ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p1 -(-3 * WWp1 + 8 * Wp1 + (NEEp1 * 2 - NNEEp1)) / 6);
	//ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p2 -(-3 * WWp2 + 8 * Wp2 + (NEEp2 * 2 - NNEEp2)) / 6);
#endif
	ctx->context[++j] = ((W + NEE) / 2 + p1 - (Wp1 + NEEp1) / 2);
	ctx->context[++j] = ((W + NEE) / 2 + p2 - (Wp2 + NEEp2) / 2);
#if 1
	//ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p1 - (WWp1 + (NEEp1 * 2 - NNEEp1)) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p2 - (WWp2 + (NEEp2 * 2 - NNEEp2)) / 2);
	ctx->context[++j] = (WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1));
	ctx->context[++j] = (WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2));
	ctx->context[++j] = (W + N - NW);
	ctx->context[++j] = (W + N - NW + p1 - (Wp1 + Np1 - NWp1));
	ctx->context[++j] = (W + N - NW + p2 - (Wp2 + Np2 - NWp2));
	ctx->context[++j] = (W + NE - N);
	ctx->context[++j] = (N + NW - NNW);
	ctx->context[++j] = (N + NW - NNW + p1 - (Np1 + NWp1 - NNEp1));
	ctx->context[++j] = (N + NW - NNW + p2 - (Np2 + NWp2 - NNEp2));
	ctx->context[++j] = (N + NE - NNE);
	ctx->context[++j] = (N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1));
	ctx->context[++j] = (N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2));
	ctx->context[++j] = (N + NN - NNN);
	//ctx->context[++j] = (N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1));
	ctx->context[++j] = (N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2));
	ctx->context[++j] = (W + WW - WWW);
	//ctx->context[++j] = (W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1));
	ctx->context[++j] = (W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2));
	ctx->context[++j] = (W + NEE - NE);
	ctx->context[++j] = (W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1));
	ctx->context[++j] = (W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2));
#endif
	ctx->context[++j] = (NN + p1 - NNp1);
	ctx->context[++j] = (NN + p2 - NNp2);
#if 1
	ctx->context[++j] = (NN + W - NNW);
	ctx->context[++j] = (NN + W - NNW + p1 - (NNp1 + Wp1 - NNEp1));
	ctx->context[++j] = (NN + W - NNW + p2 - (NNp2 + Wp2 - NNEp2));
	ctx->context[++j] = (NN + NW - NNNW);
	//ctx->context[++j] = (NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1));
	ctx->context[++j] = (NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2));
	ctx->context[++j] = (NN + NE - NNNE);
	//ctx->context[++j] = (NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1));
	ctx->context[++j] = (NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2));
	ctx->context[++j] = (NN + NNNN - NNNNNN);
	//ctx->context[++j] = (NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1));
	ctx->context[++j] = (NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2));
#endif
	ctx->context[++j] = (WW + p1 - WWp1);
	ctx->context[++j] = (WW + p2 - WWp2);
#if 1
	ctx->context[++j] = (WW + WWWW - WWWWWW);
	//ctx->context[++j] = (WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1));
	ctx->context[++j] = (WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2));
	//ctx->context[++j] = (N * 2 - NN + p1 - (Np1 * 2 - NNp1));//
	ctx->context[++j] = (N * 2 - NN + p2 - (Np2 * 2 - NNp2));
	//ctx->context[++j] = (W * 2 - WW + p1 - (Wp1 * 2 - WWp1));
	ctx->context[++j] = (W * 2 - WW + p2 - (Wp2 * 2 - WWp2));
	//ctx->context[++j] = (N * 3 - NN * 3 + NNN);
	ctx->context[++j] = clamp4(N * 3 - NN * 3 + NNN, W, NW, N, NE);
	ctx->context[++j] = clamp4(W * 3 - WW * 3 + WWW, W, NW, N, NE);
#endif
	ctx->context[++j] = clamp4(N * 2 - NN, W, NW, N, NE);
	//ctx->context[++j] = ((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW, W, NW, N, NN)) / 6);
	//ctx->context[++j] = ((NNNEEE - 4 * NNEE + 6 * NE + (W * 4 - NW * 6 + NNW * 4 - NNNW)) / 4);
	//ctx->context[++j] = (((N + 3 * NW) / 4) * 3 - ((NNW + NNWW) / 2) * 3 + (NNNWW * 3 + NNNWWW) / 4);
	//ctx->context[++j] = ((W * 2 + NW) - (WW + 2 * NWW) + NWWW);
	//ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;
	ctx->context[++j] = (N + W + 1) >> 1;
	ctx->context[++j] = (NEEEE + NEEEEEE + 1) >> 1;
	ctx->context[++j] = (WWWWWW + WWWW + 1) >> 1;
#if 1
	ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;
#endif
	ctx->context[++j] = N;
	ctx->context[++j] = NN;
	ctx->context[++j] = N + p1 - Np1;
	ctx->context[++j] = N + p2 - Np2;
	ctx->context[++j] = W + p1 - Wp1;
	ctx->context[++j] = W + p2 - Wp2;
	ctx->context[++j] = NW + p1 - NWp1;
	ctx->context[++j] = NW + p2 - NWp2;
	ctx->context[++j] = NE + p1 - NEp1;
	ctx->context[++j] = NE + p2 - NEp2;
	ctx->context[++j] = NN + p1 - NNp1;
	ctx->context[++j] = NN + p2 - NNp2;
	ctx->context[++j] = WW + p1 - WWp1;
	ctx->context[++j] = WW + p2 - WWp2;
#if 1
	ctx->context[++j] = W + N - NW;
	ctx->context[++j] = W + N - NW + p1 - Wp1 - Np1 + NWp1;
	ctx->context[++j] = W + N - NW + p2 - Wp2 - Np2 + NWp2;
	ctx->context[++j] = W + NE - N;
	ctx->context[++j] = W + NE - N + p1 - Wp1 - NEp1 + Np1;
	ctx->context[++j] = W + NE - N + p2 - Wp2 - NEp2 + Np2;
	ctx->context[++j] = W + NEE - NE;
	ctx->context[++j] = W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1;
	ctx->context[++j] = W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2;
	ctx->context[++j] = N + NN - NNN;
	//ctx->context[++j] = N + NN - NNN + p1 - Np1 - NNp1 + NNNp1;
	ctx->context[++j] = N + NN - NNN + p2 - Np2 - NNp2 + NNNp2;
	ctx->context[++j] = N + NE - NNE;
	ctx->context[++j] = N + NE - NNE + p1 - Np1 - NEp1 + NNEp1;
	ctx->context[++j] = N + NE - NNE + p2 - Np2 - NEp2 + NNEp2;
	ctx->context[++j] = N + NW - NNW;
	ctx->context[++j] = N + NW - NNW + p1 - Np1 - NWp1 + NNEp1;
	ctx->context[++j] = N + NW - NNW + p2 - Np2 - NWp2 + NNEp2;
	ctx->context[++j] = NE + NW - NN;
	ctx->context[++j] = NE + NW - NN + p1 - NEp1 - NWp1 + NNp1;
	ctx->context[++j] = NE + NW - NN + p2 - NEp2 - NWp2 + NNp2;
	ctx->context[++j] = NW + W - NWW;
	ctx->context[++j] = NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1;
	ctx->context[++j] = NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2;
	ctx->context[++j] = W * 2 - WW;
	//ctx->context[++j] = W * 2 - WW + p1 - Wp1 * 2 + WWp1;
	ctx->context[++j] = W * 2 - WW + p2 - Wp2 * 2 + WWp2;
	ctx->context[++j] = N * 2 - NN;
	//ctx->context[++j] = N * 2 - NN + p1 - Np1 * 2 + NNp1;
	ctx->context[++j] = N * 2 - NN + p2 - Np2 * 2 + NNp2;
	ctx->context[++j] = NW * 2 - NNWW;
	//ctx->context[++j] = NW * 2 - NNWW + p1 - NWp1 * 2 + NNWWp1;
	ctx->context[++j] = NW * 2 - NNWW + p2 - NWp2 * 2 + NNWWp2;
	ctx->context[++j] = NE * 2 - NNEE;
	//ctx->context[++j] = NE * 2 - NNEE + p1 - NEp1 * 2 + NNEEp1;
	ctx->context[++j] = NE * 2 - NNEE + p2 - NEp2 * 2 + NNEEp2;
	//ctx->context[++j] = N * 3 - NN * 3 + NNN + p1 - Np1 * 3 + NNp1 * 3 - NNNp1;
	//ctx->context[++j] = N * 3 - NN * 3 + NNN + p2 - Np2 * 3 + NNp2 * 3 - NNNp2;
	//ctx->context[++j] = N * 3 - NN * 3 + NNN;
	ctx->context[++j] = (W + NE * 2 - NNE + 1) >> 1;
	//ctx->context[++j] = (W + NE * 3 - NNE * 3 + NNNE+1) >> 1;
	//ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1 * 2 - NNEp1) / 2;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2 * 2 - NNEp2) / 2;
	ctx->context[++j] = NNE + NE - NNNE;
	ctx->context[++j] = NNE + W - NN;
	ctx->context[++j] = NNW + W - NNWW;
#endif

	//if(j+1!=T36_NESTIMATORS)
	//	LOG_ERROR("j %d, estimators %d", j, T36_NESTIMATORS);
	ctx->nestimators=j+1;
	//for(int k=0;k<ctx->nestimators;++k)
	//	ctx->context[k]<<=8;
#endif
}
void t36_ctx_estimate_p0(T36Ctx *ctx, int kc, int kb, int kx)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];
	int p0idx=0;
	long long sum;
	int p0temp;
	RBNodeHandle *hnode;
	T36Node *node;

	for(int k=0;k<ctx->nestimators;++k)
	{
		hnode=map_insert(ctx->maps[workidx]+k, ctx->context+k, ctx->found+k);
		node=ctx->node[k]=(T36Node*)hnode[0]->data;
		if(ctx->found[k])
		{
			sum=node->n[0]+node->n[1];
			p0temp=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
			ctx->p0arr[k]=p0temp;
			p0temp-=0x8000;
			ctx->t[k]=t36_stretch(p0temp);
		}
		else
		{
			ctx->p0arr[k]=0x8000;
			ctx->t[k]=0;
		}
	}

	sum=0;
	ctx->wsum=0;
	int previdx=kx-1;
	MODVAR(previdx, previdx, ctx->iw);
	for(int k=0;k<ctx->nestimators;++k)
	{
		int *ebuf=ctx->ebufs[T36_NESTIMATORS*workidx+k];
		ctx->ei[k]=(ebuf[previdx]+ebuf[kx])>>1;//(left+top)/2

		//int previdx=ctx->eidx+k-24*T36_NESTIMATORS;
		//MODVAR(previdx, previdx, ctx->ebuflen);
		//ctx->ei[k]=(ctx->ebuf[previdx]+ctx->ebuf[(ctx->eidx+k)%ctx->ebuflen])>>1;

		ctx->wsum+=((long long)wk[k]<<16)/ctx->ei[k];
		sum+=ctx->t[k]*((long long)wk[k]<<16)/ctx->ei[k];
	}
	ctx->x=ctx->wsum?(int)(sum/ctx->wsum):0;
	
	//ctx->x=wk[T36_NESTIMATORS];//bias
	//for(int k=0;k<T36_NESTIMATORS;++k)
	//	ctx->x+=(int)((long long)ctx->t[k]*wk[k]>>16);
	
	ctx->s0=t36_squish(ctx->x);
	ctx->p0=ctx->s0+0x8000;
	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t36_ctx_update(T36Ctx *ctx, int kc, int kb, int bit, int kx)
{
	int workidx=kc<<3|kb;
#ifdef T36_PRINT_ESTIMATOR_CR
	for(int k=0;k<ctx->nestimators;++k)
	{
		int prob=(unsigned short)(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			//if(isfinite(bitsize))
				ctx->csizes[T36_NESTIMATORS*workidx+k]+=bitsize;
			//else
			//	LOG_ERROR("p0i = 0x%04X", prob);
		}
	}
#endif

	//bwd
	if((unsigned)(ctx->s0+0x7FFF)<0xFFFF)//1 <= s0 < 0x10000
	{
		int pbit=bit?0x10000-ctx->p0:ctx->p0;
		int dL_dp0=-(int)(0x100000000/pbit);
		dL_dp0^=-bit;
		dL_dp0+=bit;
		int ds0_dx=t36_squish_dash(ctx->x, ctx->s0);
		//int ds0_dx=(int)((long long)ctx->s0*ctx->s0>>16);
		//ds0_dx*=T36_ACT_ROUGHNESS;
		//ds0_dx/=T36_ACT_ROUGHNESS+1;
		//ds0_dx=0x10000-ds0_dx;
		//ds0_dx=(int)(ctx->x ? (long long)ds0_dx*ctx->s0/ctx->x : T36_SQRT_ONEPLUS_ROUGHNESS>>16);
		int grad=(int)((long long)dL_dp0*ds0_dx>>16);
		int delta=(int)(T36_LR*grad>>16);

		int *wk=ctx->weights[workidx];

		int wnew;
		for(int k=0;k<ctx->nestimators;++k)
		{
			wnew=wk[k]-(int)((long long)delta*(ctx->t[k]-ctx->x)/((long long)ctx->ei[k]*ctx->wsum>>16));
			wnew=CLAMP(1, wnew, 0x10000);
			wk[k]=wnew;
		}

		//wk[T36_NESTIMATORS]-=delta;
		//for(int k=0;k<T36_NESTIMATORS;++k)
		//	wk[k]-=(int)((long long)delta*ctx->t[k]>>16);
	}
	int p0_ideal=bit?0:0x10000, *ebuf;
	for(int k=0;k<ctx->nestimators;++k)
	{
		ebuf=ctx->ebufs[T36_NESTIMATORS*workidx+k];
		ebuf[kx]=abs(p0_ideal-ctx->p0arr[k])+1;
	}

	//for(int k=0;k<ctx->nestimators;++k)
	//	ctx->ebuf[(ctx->eidx+k)%ctx->ebuflen]=abs(p0_ideal-ctx->p0arr[k])+1;
	//ctx->eidx+=T36_NESTIMATORS;
	//ctx->eidx%=ctx->ebuflen;

	//update
	T36Node *node;
	for(int k=0;k<ctx->nestimators;++k)
	{
		node=ctx->node[k];
		if(ctx->found[k])
			++node->n[bit];
		else
		{
			node->key=ctx->context[k];
			node->n[0]=1;
			node->n[1]=1;
			++ctx->nnodes;
		}
		ctx->context[k]+=bit<<kb;
		//ctx->context[k]|=bit<<kb;
	}
}
void t36_ctx_print(T36Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];
	for(int k=0;k<ctx->nestimators;++k)
		printf("E%3d ctx 0x%08X found %d p0arr 0x%04X t 0x%04X w 0x%08X\n", k, ctx->context[k], ctx->found[k], ctx->p0arr[k], ctx->t[k], wk[k]);
	printf("x 0x%08X s0 0x%04X p0 0x%04X\n", ctx->x, ctx->s0, ctx->p0);
	//MAP_DEBUGPRINT(ctx->maps[workidx]+1, );
	printf("\n");
}
T36Ctx t36_ctx={0};
int t36_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	//t36_stretch_squish_test();//

	int res=iw*ih;
	double t_start=time_sec();
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);

	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//colortransform_ycocb_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);

	apply_transforms_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 192);//makes MSB easy

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	int hits[24]={0};
	
	t36_ctx_init(&t36_ctx, iw, ih);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				//if(kc==0&&kx==0&&ky==1)//
				//	printf("");
				t36_ctx_get_context(&t36_ctx, (char*)buf2, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==0&&kx==2&&ky==0&&kb==7)//
					//	printf("");

					t36_ctx_estimate_p0(&t36_ctx, kc, kb, kx);
					
					//if(kc==1&&kx==3&&ky==0&&kb==7)//
					//	t36_ctx_print(&t36_ctx, kc, kb);
					//if(t36_ctx.p0!=0x8000)//
					//	printf("");

					int bit=buf2[(iw*ky+kx)<<2|kc]>>kb&1;
					abac_enc(&ctx, t36_ctx.p0, bit);
					
					int hit=bit?t36_ctx.p0<0x8000:t36_ctx.p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					int prob=bit?0x10000-t36_ctx.p0:t36_ctx.p0;
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					//if(t36_ctx.p0!=0x8000)//
					//	printf("");
					
					t36_ctx_update(&t36_ctx, kc, kb, bit, kx);
					//for(int k=0;k<16;++k)//X
					//{
					//	t36_ctx_estimate_p0(&t36_ctx, kc, kb);
					//	t36_ctx_update(&t36_ctx, kc, kb, bit, 0);
					//}
				}
				t36_ctx.ebuf[(iw*ky+kx)<<2|kc]-=buf2[(iw*ky+kx)<<2|kc];
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			TimeInfo ti;
			parsetimedelta(time_sec()-t_start, &ti);
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("Y%5d  CR%10f %8.2f MB  %02d-%02d-%06.3f\n", ky, iw*3/(csize-csize_prev), (float)t36_ctx.nnodes*sizeof(T36Node)/(1024*1024), ti.hours, ti.mins, ti.secs);
			csize_prev=csize;
		}
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t36_ctx.nnodes*sizeof(T36Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");

		for(int kc=0;kc<3;++kc)
		{
			printf("C%d  ", kc);
			for(int kb=7;kb>=0;--kb)
				printf("  B%d %11.3f", kb, iw*ih/csizes[kc<<3|kb]);
			printf("\n");
		}
		printf("\n");

#ifdef T36_PRINT_ESTIMATOR_CR
		printf("Estimator efficiencies:\n");
		int minidx[24]={0}, maxidx[24]={0};
		for(int kb=0;kb<24;++kb)
		{
			float *sizes=t36_ctx.csizes+T36_NESTIMATORS*kb;
			for(int ke=1;ke<t36_ctx.nestimators;++ke)
			{
				if(sizes[minidx[kb]]>sizes[ke])
					minidx[kb]=ke;
				if(sizes[maxidx[kb]]<sizes[ke])
					maxidx[kb]=ke;
			}
		}
		for(int ke=0;ke<t36_ctx.nestimators;++ke)
		{
			float *sizes=t36_ctx.csizes+ke;
			printf("E%3d ", ke);
			for(int kb=0;kb<24;++kb)
			{
				char c;
				if(ke==minidx[kb])
					c='*';
				else if(ke==maxidx[kb])
					c='L';
				else
					c=' ';
				printf(" %7.2f %c", sizes[T36_NESTIMATORS*kb]/t36_ctx.csizes[T36_NESTIMATORS*kb+minidx[kb]], c);
			}
			printf("\n");
		}
#endif
	}
	t36_ctx_clear(&t36_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t36_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	//int debug_index=0;

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t36_ctx_init(&t36_ctx, iw, ih);

	for(int ky=0;ky<ih;++ky)
	{
		//if(loud)
		//	printf("Dec Y%5d  %lf\r", ky, time_sec()-t_start);
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				//if(kc==0&&kx==0&&ky==1)//
				//	printf("");
				t36_ctx_get_context(&t36_ctx, (char*)buf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kx==3&&ky==0&&kb==7)//
					//	printf("");

					t36_ctx_estimate_p0(&t36_ctx, kc, kb, kx);

					//if(kc==1&&kx==3&&ky==0&&kb==7)//
					//	t36_ctx_print(&t36_ctx, kc, kb);
					//if(t36_ctx.p0!=0x8000)//
					//	printf("");
					
					int bit=abac_dec(&ctx, t36_ctx.p0);
					buf[(iw*ky+kx)<<2|kc]|=bit<<kb;

					t36_ctx_update(&t36_ctx, kc, kb, bit, kx);
					//for(int k=0;k<16;++k)
					//{
					//	t36_ctx_estimate_p0(&t36_ctx, kc, kb);
					//	t36_ctx_update(&t36_ctx, kc, kb, bit, 0);
					//}
				}
				t36_ctx.ebuf[(iw*ky+kx)<<2|kc]-=buf[(iw*ky+kx)<<2|kc];
			}
		}
		if(loud)
		{
			TimeInfo ti;
			parsetimedelta(time_sec()-t_start, &ti);
			printf("Y%5d %8.2f MB  %02d-%02d-%06.3f\r", ky, (float)t36_ctx.nnodes*sizeof(T36Node)/(1024*1024), ti.hours, ti.mins, ti.secs);
		}
	}
	if(loud)
		printf("\n");
	t36_ctx_clear(&t36_ctx);
	
	//addbuf((unsigned char*)buf, iw, ih, 3, 4, 128);

	//addbuf((unsigned char*)buf, iw, ih, 3, 4, 128);
	//colortransform_ycocb_fwd(buf, iw, ih);
	//addbuf((unsigned char*)buf, iw, ih, 3, 4, 128);

	apply_transforms_inv(buf, iw, ih);
	if(loud)
	{
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}
#endif


//T37: Fixed array as binary tree predictor

#define T37_NESTIMATORS 14
#define T37_PRINT_ESTIMATOR_CR

void print_weights(int *weights, int count)
{
	int sum=0;
	for(int k=0;k<count;++k)
		sum+=weights[k];
	int sum2=0;
	printf("|");
	for(int k=0, kp=0;k<count;++k)
	{
		sum2+=weights[k];
		int print=sum2*80/sum;
		for(int k2=kp;k2<print;++k2)
			printf("-");
		kp=print;
		printf("|");
	}
}
typedef struct T37CtxStruct
{
	int context[T37_NESTIMATORS], dec;
	unsigned short estimators[3][T37_NESTIMATORS][255];
	//unsigned short estimators[T37_NESTIMATORS][1+4+16+64+256+1024+4096+16384];
	int idx[T37_NESTIMATORS];
	char hintbit[T37_NESTIMATORS];
	int weights[24][T37_NESTIMATORS];
	int p0arr[T37_NESTIMATORS];
	int p0_0, p0;//p0_0 isn't clamped
	long long wsum;
#ifdef T37_PRINT_ESTIMATOR_CR
	float csizes[24*T37_NESTIMATORS];
#endif
} T37Ctx;
void t37_ctx_init(T37Ctx *ctx)
{
	int val=0x8000;
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	//for(int k=0;k<T37_NESTIMATORS;++k)//fixed 15.16 bit
	//	ctx->weights[k]=0x8000;
	memfill(ctx->estimators, &val, sizeof(ctx->estimators), sizeof(short));
}
void t37_ctx_get_context(T37Ctx *ctx, char *buf, int iw, int ih, int kc, int kx, int ky)
{
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	int W   =kx-1>=0         ?buf[(iw* ky   +kx-1)<<2| kc   ]:0,
		NW  =kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]:0,
		N   =ky-1>=0         ?buf[(iw*(ky-1)+kx  )<<2| kc   ]:0,
		NE  =kx+1<iw&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2| kc   ]:0,
		NN  =ky-2>=0         ?buf[(iw*(ky-2)+kx  )<<2| kc   ]:0,

		m1  =kc-1>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-1)]:0,
		Nm1 =kc-1>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-1)]:0,
		Wm1 =kc-1>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-1)]:0,
		NWm1=kc-1>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-1)]:0,

		m2  =kc-2>=0                  ?buf[(iw* ky   +kx  )<<2|(kc-2)]:0,
		Nm2 =kc-2>=0         &&ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|(kc-2)]:0,
		Wm2 =kc-2>=0&&kx-1>=0         ?buf[(iw* ky   +kx-1)<<2|(kc-2)]:0,
		NWm2=kc-2>=0&&kx-1>=0&&ky-1>=0?buf[(iw*(ky-1)+kx-1)<<2|(kc-2)]:0;
	
	int j=-1;

	ctx->context[++j]=0;
	ctx->context[++j]=N;
	ctx->context[++j]=W;
	ctx->context[++j]=NW;
	ctx->context[++j]=m1;
	ctx->context[++j]=W+NE-N;
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);
	ctx->context[++j]=NW+NE-NN;
	ctx->context[++j]=(N+W-NW + m1)>>1;
	ctx->context[++j]=m2;
	ctx->context[++j]=(N+W-NW + m2)>>1;

	//memset(ctx->idx, 0, sizeof(ctx->idx));
	ctx->dec=0;

	for(int k=1;k<T37_NESTIMATORS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
static void t37_access_idx(T37Ctx *ctx, int kb, int ke)
{
	int vmin=ctx->dec, vmax=vmin|((1<<(kb+1))-1);
	int pred=ctx->context[ke];
	if(ke)
		pred=CLAMP(vmin, pred, vmax);
	//if(pred!=ctx->context[ke])
	//	printf("");
	int key=ctx->dec-pred, bit;
	ctx->idx[ke]=0;
	for(int kb2=7;kb2>kb;--kb2)
	{
		bit=key>>kb2&1;
		ctx->idx[ke]=(ctx->idx[ke]<<1)+bit+1;
	}
	ctx->hintbit[ke]=key>>kb&1;
}
void t37_ctx_estimate_p0(T37Ctx *ctx, int kc, int kb)
{
	//static const int factors[]=
	//{
	//	0, 0, 0, 0, 0, 1, 1, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 1, 0, 1, 0,
	//};

	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];
	long long sum=0;
	ctx->wsum=0;
	for(int k=0;k<T37_NESTIMATORS;++k)
	{
		t37_access_idx(ctx, kb, k);
		ctx->p0arr[k]=ctx->estimators[kc][k][ctx->idx[k]];

		if(k)//first estimator has zero-predictor for plain histogram (no hint bit)
		{
			if(ctx->hintbit[k])
				ctx->p0arr[k]-=ctx->p0arr[k]>>2;
			else
				ctx->p0arr[k]+=(0x10000-ctx->p0arr[k])>>2;
		}

		sum+=(long long)ctx->p0arr[k]*wk[k];
		ctx->wsum+=wk[k];
	}
	//ctx->p0_0=ctx->wsum?(int)((sum+(ctx->wsum>>1))/ctx->wsum):0x8000;//same CR
	ctx->p0_0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;

	//if(!factors[workidx])
	//	ctx->p0=ctx->p0arr[0];
	//else
		ctx->p0=ctx->p0_0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t37_ctx_update(T37Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
#ifdef T37_PRINT_ESTIMATOR_CR
	for(int k=0;k<T37_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		prob=CLAMP(1, prob, 0xFFFF);
		float p=(float)prob/0x10000;
		float bitsize=-log2f(p);
		//if(fpclassify(bitsize)!=FP_NORMAL)
		//	LOG_ERROR("");
		ctx->csizes[T37_NESTIMATORS*workidx+k]+=bitsize;
	}
#endif
	//bwd
	int *wk=ctx->weights[kc<<3|kb];
	if(ctx->p0_0>=1&&ctx->p0_0<=0xFFFF)
	{
		int p_bit=bit?0x10000-ctx->p0:ctx->p0;
		long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
		dL_dp0^=-bit;
		dL_dp0+=bit;
		for(int k=0;k<T37_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
	for(int k=0;k<T37_NESTIMATORS;++k)
	{
		unsigned short *p0=ctx->estimators[kc][k]+ctx->idx[k];

		//*p0+=(!bit-bit)<<8;

		//if(!k)
		{
			if(bit)
				*p0-=*p0>>7;
			else
				*p0+=(0x10000-*p0)>>7;
			//*p0=CLAMP(1, *p0, 0xFFFF);
		}

		//int pred=ctx->context[k]>>kb&1;
		//int choice=bit<<1|pred;

		//ctx->idx[k]=(idx<<1)+choice+1;
		//ctx->idx[k]=(idx<<2)+choice+1;
	}
	ctx->dec|=bit<<kb;
}
T37Ctx t37_context;
int t37_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//colortransform_ycocb_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);

	apply_transforms_fwd(buf2, iw, ih);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 192);//makes MSB easy

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float csizes[24]={0};
	int hits[24]={0};
	
	t37_ctx_init(&t37_context);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			if(ky==10&&kx==10)//
				printf("");

			for(int kc=0;kc<3;++kc)
			{
				t37_ctx_get_context(&t37_context, (char*)buf2, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kx==1&&ky==0&&kb==1)
					//	printf("");
					//if((t37_context.dec>>kb)!=(t37_context.context[0]>>kb))//
					//	printf("");

					t37_ctx_estimate_p0(&t37_context, kc, kb);

					//if(t37_context.p0!=0x8000)//
					//	printf("");

					int bit=(buf2[(iw*ky+kx)<<2|kc]+128)>>kb&1;//signed -> unsigned
					ac_enc_bin(&ec, t37_context.p0, bit);
					
					int hit=bit?t37_context.p0<0x8000:t37_context.p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					int prob=bit?0x10000-t37_context.p0:t37_context.p0;
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t37_ctx_update(&t37_context, kc, kb, bit);
				}
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("Y%5d  CR%11f p0 0x%04X ", ky, iw*3/(csize-csize_prev), t37_context.p0);
			print_weights(t37_context.weights[8], T37_NESTIMATORS);
			printf("\n");
			csize_prev=csize;
		}
	}
	ac_enc_flush(&ec);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");

		for(int kc=0;kc<3;++kc)
		{
			printf("C%d  ", kc);
			for(int kb=7;kb>=0;--kb)
				printf("  B%d %11.3f", kb, iw*ih/csizes[kc<<3|kb]);
			printf("\n");
		}
		printf("\n");

#ifdef T37_PRINT_ESTIMATOR_CR
		printf("Estimator efficiencies:\n");
		int minidx[24]={0}, maxidx[24]={0};
		for(int kb=0;kb<24;++kb)
		{
			float *sizes=t37_context.csizes+T37_NESTIMATORS*kb;
			for(int ke=1;ke<T37_NESTIMATORS;++ke)
			{
				if(sizes[minidx[kb]]>sizes[ke])
					minidx[kb]=ke;
				if(sizes[maxidx[kb]]<sizes[ke])
					maxidx[kb]=ke;
			}
		}
		for(int ke=0;ke<T37_NESTIMATORS;++ke)
		{
			float *sizes=t37_context.csizes+ke;
			printf("E%2d ", ke);
			for(int kb=0;kb<24;++kb)
			{
				char c;
				if(ke==minidx[kb])
					c='*';
				else if(ke==maxidx[kb])
					c='L';
				else
					c=' ';
				printf(" %7.2f%c", iw*ih/t37_context.csizes[T37_NESTIMATORS*kb+minidx[kb]], c);
				//printf(" %7.2f%c", sizes[T37_NESTIMATORS*kb]/t37_context.csizes[T37_NESTIMATORS*kb+minidx[kb]], c);
			}
			printf("\n");
		}
#endif
	}
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t37_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	//int debug_index=0;

	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t37_ctx_init(&t37_context);

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				t37_ctx_get_context(&t37_context, (char*)buf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t37_ctx_estimate_p0(&t37_context, kc, kb);
					
					int bit=ac_dec_bin(&ec, t37_context.p0);
					buf[(iw*ky+kx)<<2|kc]|=bit<<kb;//unsigned

					t37_ctx_update(&t37_context, kc, kb, bit);
				}
				buf[(iw*ky+kx)<<2|kc]+=128;//unsigned -> signed
			}
		}
	}
	
	addbuf(buf, iw, ih, 3, 4, 128);
	apply_transforms_inv(buf, iw, ih);
	if(loud)
	{
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}


#define T38_PERSISTENCE 7
int t38_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//colortransform_ycocb_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);

	apply_transforms_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 192);//makes MSB easy

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float csizes[24]={0};
	int hits[24]={0};

	unsigned short prob0[3][255], ctx_idx;
	unsigned short prob_init=0x8000;
	memfill(prob0, &prob_init, sizeof(prob0), 2);
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			if(ky==10&&kx==10)//
				printf("");

			for(int kc=0;kc<3;++kc)
			{
				ctx_idx=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					unsigned short *p0=prob0[kc]+ctx_idx;

					//if(*p0!=0x8000)//
					//	printf("");

					int bit=buf2[(iw*ky+kx)<<2|kc]>>kb&1;
					ac_enc_bin(&ec, *p0, bit);
					
					int hit=bit?*p0<0x8000:*p0>0x8000;//
					hits[kc<<3|kb]+=hit;
					int prob=bit?0x10000-*p0:*p0;
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					if(bit)
						*p0-=*p0>>T38_PERSISTENCE;
					else
						*p0+=(0x10000-*p0)>>T38_PERSISTENCE;
					*p0=CLAMP(1, *p0, 0xFFFF);

					ctx_idx<<=1;
					ctx_idx+=bit+1;
				}
			}
		}
		if(loud==2)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("Y %d  CR %f\n", ky, iw*3/(csize-csize_prev));
			csize_prev=csize;
		}
	}
	ac_enc_flush(&ec);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld  CR %lf  P %7d/8 = %g\n", list.nobj, 3.*iw*ih/list.nobj, iw*ih, iw*ih/8.);
		printf("\n");

		for(int kc=0;kc<3;++kc)
		{
			printf("C%d  ", kc);
			for(int kb=7;kb>=0;--kb)
				printf("  B%d %11.3f", kb, iw*ih/csizes[kc<<3|kb]);
			printf("\n");
		}
		printf("\n");
	}
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t38_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	
	unsigned short prob0[3][255], ctx_idx;
	unsigned short prob_init=0x8000;
	memfill(prob0, &prob_init, sizeof(prob0), 2);

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				ctx_idx=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					unsigned short *p0=prob0[kc]+ctx_idx;
					
					int bit=ac_dec_bin(&ec, *p0);
					buf[(iw*ky+kx)<<2|kc]|=bit<<kb;
					
					if(bit)
						*p0-=*p0>>T38_PERSISTENCE;
					else
						*p0+=(0x10000-*p0)>>T38_PERSISTENCE;
					*p0=CLAMP(1, *p0, 0xFFFF);

					ctx_idx<<=1;
					ctx_idx+=bit+1;
				}
			}
		}
	}
	
	//addbuf(buf, iw, ih, 3, 4, 128);
	apply_transforms_inv(buf, iw, ih);
	if(loud)
	{
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}




//T39: Multiple estimators for all maps

//	#define T39_USE_KALMAN
	#define T39_UNSIGNED_BITS
//	#define T39_APPLY_SPATIAL
	#define T39_USE_ARRAYS
//	#define T39_DISABLE_REC
//	#define T39_PROB_TWEAK	//X

#define T39_NMAPS 15	//14	31		135 HALF HOUR PER IMAGE

#ifndef T39_DISABLE_REC
#define T39_N_REC_ESTIMATORS 6		//15
#define T39_NESTIMATORS ((T39_N_REC_ESTIMATORS+1)*T39_NMAPS)
#else
#define T39_NESTIMATORS T39_NMAPS
#endif
//#define T39_PRINT_ESTIMATOR_CR

//#define T39_BLOCKSIZE_X 2048	//512
//#define T39_BLOCKSIZE_Y 2048	//512
//	#define T39_DISABLE_COUNTER

#ifdef T39_USE_KALMAN
typedef struct KalmanInfoStruct
{
	double R, H, Q, P, Uhat, K;
} KalmanInfo;
static void kalman_init(KalmanInfo *info)
{
	info->R=0.001;//noise covariance, higher values reject more noise
	info->H=1;
	info->Q=10;
	info->P=0;
	info->Uhat=0;
	info->K=0;
}
static void kalman_predict(KalmanInfo *info, double U)
{
	info->K=info->P*info->H/(info->H*info->P*info->H+info->R);//update kalman gain
	info->Uhat+=info->K*(U-info->H*info->Uhat);//update estimated
	info->P=(1-info->K*info->H)*info->P+info->Q;//update error covariance
}
#endif
typedef struct T39NodeStruct
{
#ifndef T39_USE_ARRAYS
	int key;
#endif
	int n[2];
#ifndef T39_DISABLE_REC
	unsigned short rec[T39_N_REC_ESTIMATORS];
#endif
} T39Node;
#ifndef T39_USE_ARRAYS
static CmpRes t39_cmp_node(const void *key, const void *candidate)
{
	int const *k=(int const*)key;
	T35CtxNode const *c=(T35CtxNode const*)candidate;
	return (*k>c->key)-(*k<c->key);
}
//static void t39_debugprinter(RBNodeHandle *node, int depth)
//{
//	T35CtxNode *p;
//	if(node)
//	{
//		p=(T35CtxNode*)node[0]->data;
//		printf("0x%08X %5d %5d\n", p->key, p->n[0], p->n[1]);
//	}
//}
#endif
typedef struct T39CtxStruct
{
	int pred14;
#ifdef T39_USE_KALMAN
	KalmanInfo kalman[T39_NMAPS];
#endif
	int context[T39_NMAPS];
#ifdef T39_USE_ARRAYS
	ArrayHandle maps[24][T39_NMAPS];//3*(256+512+1024+2048+4096+8192+16384+32768)*15*sizeof(T39Node) = 56.03 MB for 14 maps with 6 rec estimators
#else
	Map maps[24][T39_NMAPS];
	int found[T39_NMAPS];
#endif
	T39Node *node[T39_NMAPS];

	int p0arr[T39_NESTIMATORS], p0_0, p0, p0rev;//p0_0 isn't clamped
	int weights[24][T39_NESTIMATORS];
	long long wsum;

#ifdef T39_PROB_TWEAK
	long long proberrors[24];
	int bitsprocessed;
	int hits[24];
	double csizes[24];//needs fixed prec
#endif

	int nnodes;
#ifdef T39_PRINT_ESTIMATOR_CR
	float csizes_est[24*T39_NESTIMATORS];
#endif
} T39Ctx;
void t39_explore(T39Ctx *ctx, const char *src, int iw, int ih)
{
	int res=iw*ih;
	int hist[256]={0};
	int kc=0,//red (orangeness)
		ke=0;//zero predictor
	//double entropy=0;
	for(int sym=0;sym<256;++sym)
	{
		if(sym==128)//
			printf("");

		int prob=0x10000;
		int context=0x80;
		for(int kb=7;kb>=0;--kb)
		{
			ArrayHandle map=ctx->maps[kc<<3|kb][ke];
			int bit=sym>>kb&1;
			T39Node *node=array_at(&map, context);

			//int pb=bit?0x10000-node->rec[3]:node->rec[3];
			int pb=(int)(((long long)node->n[bit]<<16)/(node->n[0]+node->n[1]));

			pb=CLAMP(1, pb, 0xFFFF);
			prob=(int)(((long long)prob*pb+0x8000)>>16);
			context|=bit<<(8+7-kb);
		}
		hist[sym]=prob;

		//printf("%3d  0x%04X\n", sym, prob);
		//if(prob)
		//{
		//	double p=(double)prob/res;
		//	entropy-=p*log2(p);				//X  need to use cross-entropy with image histogram
		//}
	}
	//double invCR=entropy/8;
	//printf("CR %lf\n", 1/invCR);

	double csize=0;
	for(int k=0;k<res;++k)//Zipf's law
	{
		unsigned char sym=src[k<<2|kc]+128;
		int prob=hist[sym];
		if(prob)
		{
			double p=(double)prob/res;
			csize-=log2(p);
		}
	}
	csize/=8;
	printf("C%d  csize %lf  CR %lf\n", kc, csize, res/csize);
}
T39Ctx* t39_ctx_init()
{
	int val=0x8000;
#ifdef T39_USE_ARRAYS
	T39Node node0={{1, 1}};
#ifndef T39_DISABLE_REC
	for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
#endif
	T39Ctx *ctx=(T39Ctx*)malloc(sizeof(T39Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T39Ctx));
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T39_NMAPS;++k2)
		{
#ifdef T39_USE_ARRAYS
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T39Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T39Node), sizeof(T39Node));
			ctx->nnodes+=nnodes;
#else
			MAP_INIT(ctx->maps[k]+k2, T39Node, t39_cmp_node, 0);
#endif
		}
	}
#ifdef T39_USE_KALMAN
	for(int k=0;k<T39_NMAPS;++k)
		kalman_init(ctx->kalman+k);
#endif
	return ctx;
}
void t39_ctx_clear(T39Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T39_NMAPS;++k2)
#ifdef T39_USE_ARRAYS
			array_free(ctx[0]->maps[k]+k2);
#else
			MAP_CLEAR(ctx[0]->maps[k]+k2);
#endif
	}
	free(*ctx);
	*ctx=0;
}
void t39_ctx_reset(T39Ctx *ctx, int hardreset)
{
#ifdef T39_USE_ARRAYS
	T39Node node0={{1, 1}};
#ifndef T39_DISABLE_REC
	for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
#endif
	for(int k=0;k<24;++k)
	{
		if(hardreset)
		{
			for(int k2=0;k2<T39_NMAPS;++k2)
#ifdef T39_USE_ARRAYS
				memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T39Node), sizeof(T39Node));
#else
				MAP_CLEAR(ctx->maps[k]+k2);
#endif
		}
#ifdef T39_PROB_TWEAK
		ctx->proberrors[k]=0;
#endif
	}
	if(hardreset)
		ctx->nnodes=0;
#ifdef T39_PROB_TWEAK
	ctx->bitsprocessed=0;
#endif
}
//void t39_ctx_get_context(T39Ctx *ctx, char *buf, int iw, int ih, int kc, int kx, int ky, int x1, int x2, int y1, int y2)
void t39_ctx_get_context(T39Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
//#define LOAD(C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X)-x1)<(unsigned)(x2-x1)&&(unsigned)(ky-Y-y1)<(unsigned)(y2-y1)?buf[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
//#define LOAD(CO, XO, YO) (unsigned)(kc-CO)<3u&&(unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?buf[(iw*(ky-YO)+kx-(XO))<<2|(kc-CO)]:0
#if 1
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	char
		NNWW =LOAD(buf, 0,  2, 2),
		NNW  =LOAD(buf, 0,  1, 2),
		NN   =LOAD(buf, 0,  0, 2),
		NNE  =LOAD(buf, 0, -1, 2),
		NNEE =LOAD(buf, 0, -2, 2),
		NWW  =LOAD(buf, 0,  2, 1),
		NW   =LOAD(buf, 0,  1, 1),
		N    =LOAD(buf, 0,  0, 1),
		NE   =LOAD(buf, 0, -1, 1),
		NEE  =LOAD(buf, 0, -2, 1),
		WW   =LOAD(buf, 0,  2, 0),
		W    =LOAD(buf, 0,  1, 0),
		eNNWW=LOAD(ebuf, 0,  2, 2),
		eNNW =LOAD(ebuf, 0,  1, 2),
		eNN  =LOAD(ebuf, 0,  0, 2),
		eNNE =LOAD(ebuf, 0, -1, 2),
		eNNEE=LOAD(ebuf, 0, -2, 2),
		eNWW =LOAD(ebuf, 0,  2, 1),
		eNW  =LOAD(ebuf, 0,  1, 1),
		eN   =LOAD(ebuf, 0,  0, 1),
		eNE  =LOAD(ebuf, 0, -1, 1),
		eNEE =LOAD(ebuf, 0, -2, 1),
		eWW  =LOAD(ebuf, 0,  2, 0),
		eW   =LOAD(ebuf, 0,  1, 0),

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);
#if 0
	int W   =LOAD(buf, 0,  1, 0),
		NW  =LOAD(buf, 0,  1, 1),
		N   =LOAD(buf, 0,  0, 1),
		NE  =LOAD(buf, 0, -1, 1),
		NN  =LOAD(buf, 0,  0, 2),

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);
#endif

	int j=-1;

	//bit, channel-bitplane, compressibility			based on kodim13
	//
	//Orangeness:			best pred				worst pred
	// 0	0-0		*		(N+W)/2					0
	// 1	0-1		**		(N+W)/2					0
	// 2	0-2		***		(N+W)/2					NW+NE-NN
	// 3	0-3		****	(N+W)/2					NW+NE-NN
	// 4	0-4		****	W						NW+NE-NN
	// 5	0-5		****	0						NW+NE-NN
	// 6	0-6		****	0						NW+NE-NN
	// 7	0-7		*		NW+NE-NN				0
	//
	//Luma:
	// 8	1-0		*		NW+NE-NN				NW+NE-NN
	// 9	1-1		*		(W+N+m1)/3				NW+NE-NN
	//10	1-2		*		(W+N+m1)/3				NW+NE-NN
	//11	1-3		*		W						0
	//12	1-4		**		W						0
	//13	1-5		***		(N+W-NW + m2)>>1		NW+NE-NN
	//14	1-6		****	(W+N+m1)/3				NW+NE-NN
	//15	1-7		*		NW+NE-NN				0
	//
	//Blueness:
	//16	2-0		*		W						clamp4(N+m1-Nm1, N, m1, Nm1, NW)
	//17	2-1		**		N						clamp4(N+m1-Nm1, N, m1, Nm1, NW)
	//18	2-2		***		W						(N+W-NW + m1)>>1
	//19	2-3		****	(N+W-NW + m2)>>1		(N+W-NW + m1)>>1
	//20	2-4		****	(N+W-NW + m2)>>1		(N+W-NW + m1)>>1
	//21	2-5		****	m2						(N+W-NW + m1)>>1
	//22	2-6		****	0						(N+W-NW + m1)>>1
	//23	2-7		*		m2						0

	ctx->context[++j]=0;//0
	ctx->context[++j]=N;//1
	ctx->context[++j]=W;//2
	ctx->context[++j]=NW;//3
	ctx->context[++j]=m1;//4
	ctx->context[++j]=W+NE-N;//5
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;//6
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);//7
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);//8
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);//9
	ctx->context[++j]=NW+NE-NN;//10
	//ctx->context[++j]=(N+W+NW+NE)>>2;//10_v2
	ctx->context[++j]=(N+W-NW + m1)>>1;//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=(N+W-NW + m2)>>1;//13

	//kodim13
#if 1
	switch(kc)
	{
	case 0:
		ctx->pred14=(
			+0x00C6* NNWW-0x0267* NNW+0x0400* NN-0x04C4* NNE+0x0188* NNEE
			-0x0166* NWW -0x0013* NW +0x034D* N +0x05BF* NE +0x002C* NEE
			+0x03E6* WW  +0x05BA* W
			+0x0076*eNNWW+0x0228*eNNW-0x001D*eNN+0x0241*eNNE-0x0167*eNNEE
			+0x00DA*eNWW +0x0404*eNW +0x0468*eN -0x0008*eNE +0x00C5*eNEE
			+0x0077*eWW  +0x05E7*eW
		)>>12;
		break;
	case 1:
		ctx->pred14=(
			+0x0010* NNWW+0x0020* NNW-0x0163* NN-0x000E* NNE+0x012B* NNEE
			-0x000E* NWW +0x02A6* NW +0x02F1* N +0x0081* NE +0x01D5* NEE
			+0x0215* WW  +0x063A* W
			-0x0074*eNNWW-0x015D*eNNW-0x00CC*eNN-0x00EB*eNNE-0x0156*eNNEE
			-0x0167*eNWW +0x0011*eNW +0x04A7*eN +0x00A5*eNE -0x00B2*eNEE
			-0x01D2*eWW  +0x0650*eW
		)>>12;
		break;
	case 2:
		ctx->pred14=(
			+0x00AC* NNWW-0x02B4* NNW+0x021E* NN-0x0010* NNE+0x0036* NNEE
			-0x0257* NWW +0x0076* NW +0x054D* N +0x00F9* NE +0x00BF* NEE
			+0x02D3* WW  +0x07D5* W
			+0x001D*eNNWW+0x0125*eNNW+0x009D*eNN+0x00EA*eNNE+0x007F*eNNEE
			+0x00AA*eNWW +0x035E*eNW +0x0669*eN +0x03FD*eNE -0x0044*eNEE
			+0x005C*eWW  +0x0526*eW
		)>>12;
		break;
	}
	ctx->pred14=CLAMP(-128, ctx->pred14, 127);
	ctx->context[++j]=ctx->pred14;//14
#endif
	
	//CLIC16
#if 0
	switch(kc)
	{
	case 0:
		ctx->pred14=(
			+0x00C8* NNWW-0x01B9* NNW+0x01CB* NN+0x0170* NNE-0x00E7* NNEE
			-0x01DA* NWW +0x00A8* NW +0x03FD* N +0x01DD* NE +0x00AB* NEE
			+0x00A5* WW  +0x0900* W
			+0x00FF*eNNWW+0x0040*eNNW-0x02B6*eNN+0x000D*eNNE+0x0182*eNNEE
			+0x0031*eNWW +0x00DE*eNW +0x065B*eN +0x0220*eNE +0x0056*eNEE
			-0x02EB*eWW  +0x036E*eW
		)>>12;
		break;
	case 1:
		ctx->pred14=(
			+0x0080* NNWW-0x00DF* NNW-0x00FC* NN+0x0195* NNE+0x0052* NNEE
			+0x0146* NWW -0x021E* NW +0x050F* N +0x0285* NE -0x012A* NEE
			-0x00E5* WW  +0x0AD7* W
			-0x003E*eNNWW+0x0126*eNNW-0x011B*eNN-0x011B*eNNE+0x0049*eNNEE
			-0x008A*eNWW +0x0174*eNW +0x048C*eN +0x007B*eNE +0x025C*eNEE
			-0x0199*eWW  +0x0428*eW
		)>>12;
		break;
	case 2:
		ctx->pred14=(
			-0x00B6* NNWW-0x0040* NNW+0x009D* NN-0x00E4* NNE+0x0088* NNEE
			-0x0263* NWW +0x02C8* NW +0x0420* N +0x0300* NE +0x0034* NEE
			+0x0131* WW  +0x07F9* W
			+0x0139*eNNWW-0x006E*eNNW-0x00E6*eNN+0x0019*eNNE-0x0020*eNNEE
			-0x0099*eNWW -0x002F*eNW +0x0693*eN +0x014F*eNE +0x003C*eNEE
			-0x0124*eWW  +0x05CC*eW
		)>>12;
		break;
	}
	ctx->pred14=CLAMP(-128, ctx->pred14, 127);
	ctx->context[++j]=ctx->pred14;//14
#endif

	//ctx->context[++j]=clamp4((W+NE-N + NW+NE-NN)>>1, N, W, NW, NE);
	//ctx->context[++j]=Nm1+Wm1-NWm1;
	//ctx->context[++j]=Nm1;
	//ctx->context[++j]=Wm1;
	//ctx->context[++j]=NWm1;
#endif
	
#if 0
	//offsets are in NW direction
	int WWWWWW  =LOAD(buf, 0,  6, 0),
		WWWWW   =LOAD(buf, 0,  5, 0),
		WWWW    =LOAD(buf, 0,  4, 0),
		WWW     =LOAD(buf, 0,  3, 0),
		WW      =LOAD(buf, 0,  2, 0),
		W       =LOAD(buf, 0,  1, 0),
		NWWWW   =LOAD(buf, 0,  4, 1),
		NWWW    =LOAD(buf, 0,  3, 1),
		NWW     =LOAD(buf, 0,  2, 1),
		NW      =LOAD(buf, 0,  1, 1),
		N       =LOAD(buf, 0,  0, 1),
		NE      =LOAD(buf, 0, -1, 1),
		NEE     =LOAD(buf, 0, -2, 1),
		NEEE    =LOAD(buf, 0, -3, 1),
		NEEEE   =LOAD(buf, 0, -4, 1),
		NEEEEEE =LOAD(buf, 0, -6, 1),
		NNWWW   =LOAD(buf, 0,  3, 2),
		NNWW    =LOAD(buf, 0,  2, 2),
		NNW     =LOAD(buf, 0,  1, 2),
		NN      =LOAD(buf, 0,  0, 2),
		NNE     =LOAD(buf, 0, -1, 2),
		NNEE    =LOAD(buf, 0, -2, 2),
		NNEEE   =LOAD(buf, 0, -3, 2),
		NNNWWWW =LOAD(buf, 0,  4, 3),
		NNNWWW  =LOAD(buf, 0,  3, 3),
		NNNWW   =LOAD(buf, 0,  2, 3),
		NNNW    =LOAD(buf, 0,  1, 3),
		NNN     =LOAD(buf, 0,  0, 3),
		NNNE    =LOAD(buf, 0, -1, 3),
		NNNEE   =LOAD(buf, 0, -2, 3),
		NNNEEE  =LOAD(buf, 0, -3, 3),
		NNNNW   =LOAD(buf, 0,  1, 4),
		NNNN    =LOAD(buf, 0,  0, 4),
		NNNNE   =LOAD(buf, 0, -1, 4),
		NNNNN   =LOAD(buf, 0,  0, 5),
		NNNNNN  =LOAD(buf, 0,  0, 6),
		WWWWWWp1=LOAD(buf, 1,  6, 0),
		WWWWp1  =LOAD(buf, 1,  4, 0),
		WWWp1   =LOAD(buf, 1,  3, 0),
		WWp1    =LOAD(buf, 1,  2, 0),
		Wp1     =LOAD(buf, 1,  1, 0),
		p1      =LOAD(buf, 1,  0, 0),
		NWWp1   =LOAD(buf, 1,  2, 1),
		NWp1    =LOAD(buf, 1,  1, 1),
		Np1     =LOAD(buf, 1,  0, 1),
		NEp1    =LOAD(buf, 1, -1, 1),
		NEEp1   =LOAD(buf, 1, -2, 1),
		NNWWp1  =LOAD(buf, 1,  2, 2),
		NNp1    =LOAD(buf, 1,  0, 2),
		NNEp1   =LOAD(buf, 1, -1, 2),
		NNEEp1  =LOAD(buf, 1, -2, 2),
		NNNWp1  =LOAD(buf, 1,  1, 3),
		NNNp1   =LOAD(buf, 1,  0, 3),
		NNNEp1  =LOAD(buf, 1, -1, 3),
		NNNNp1  =LOAD(buf, 1,  0, 4),
		NNNNNNp1=LOAD(buf, 1,  0, 6),
		WWWWWWp2=LOAD(buf, 2,  6, 0),
		WWWWp2  =LOAD(buf, 2,  4, 0),
		WWWp2   =LOAD(buf, 2,  3, 0),
		WWp2    =LOAD(buf, 2,  2, 0),
		Wp2     =LOAD(buf, 2,  1, 0),
		p2      =LOAD(buf, 2,  0, 0),
		NWWp2   =LOAD(buf, 2,  2, 1),
		NWp2    =LOAD(buf, 2,  1, 1),
		Np2     =LOAD(buf, 2,  0, 1),
		NEp2    =LOAD(buf, 2, -1, 1),
		NEEp2   =LOAD(buf, 2, -2, 1),
		NNWWp2  =LOAD(buf, 2,  2, 2),
		NNp2    =LOAD(buf, 2,  0, 2),
		NNEp2   =LOAD(buf, 2, -1, 2),
		NNEEp2  =LOAD(buf, 2, -2, 2),
		NNNWp2  =LOAD(buf, 2,  1, 3),
		NNNp2   =LOAD(buf, 2,  0, 3),
		NNNEp2  =LOAD(buf, 2, -1, 3),
		NNNNp2  =LOAD(buf, 2,  0, 4),
		NNNNNNp2=LOAD(buf, 2,  0, 6);
	int j=-1;

	//ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;//#74
	//ctx->context[++j] = (N + W + 1) >> 1;//#71
	//ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;//#70
	//ctx->context[++j] = N;//#75

	ctx->context[++j]=0;
	ctx->context[++j] = clamp4(N + p1 - Np1, W, NW, N, NE);
	ctx->context[++j] = clamp4(N + p2 - Np2, W, NW, N, NE);
	ctx->context[++j] = (W + clamp4(NE * 3 - NNE * 3 + NNNE, W, N, NE, NEE)) / 2;
	ctx->context[++j] = clamp4((W + clip(NE * 2 - NNE)) / 2, W, NW, N, NE);
	ctx->context[++j] = (W + NEE) / 2;
	ctx->context[++j] = ((WWW - 4 * WW + 6 * W + (NE * 4 - NNE * 6 + NNNE * 4 - NNNNE)) / 4);
	ctx->context[++j] = ((-WWWW + 5 * WWW - 10 * WW + 10 * W + clamp4(NE * 4 - NNE * 6 + NNNE * 4 - NNNNE, N, NE, NEE, NEEE)) / 5);
	ctx->context[++j] = ((-4 * WW + 15 * W + 10 * (NE * 3 - NNE * 3 + NNNE) - (NEEE * 3 - NNEEE * 3 + NNNEEE)) / 20);
	ctx->context[++j] = ((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p1 - (Wp1 + (NEp1 * 2 - NNEp1)) / 2);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p2 - (Wp2 + (NEp2 * 2 - NNEp2)) / 2);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p1 -(-3 * WWp1 + 8 * Wp1 + (NEEp1 * 2 - NNEEp1)) / 6);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p2 -(-3 * WWp2 + 8 * Wp2 + (NEEp2 * 2 - NNEEp2)) / 6);
	ctx->context[++j] = ((W + NEE) / 2 + p1 - (Wp1 + NEEp1) / 2);
	ctx->context[++j] = ((W + NEE) / 2 + p2 - (Wp2 + NEEp2) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p1 - (WWp1 + (NEEp1 * 2 - NNEEp1)) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p2 - (WWp2 + (NEEp2 * 2 - NNEEp2)) / 2);
	ctx->context[++j] = (WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1));
	ctx->context[++j] = (WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2));
	ctx->context[++j] = (W + N - NW);
	ctx->context[++j] = (W + N - NW + p1 - (Wp1 + Np1 - NWp1));
	ctx->context[++j] = (W + N - NW + p2 - (Wp2 + Np2 - NWp2));
	ctx->context[++j] = (W + NE - N);
	ctx->context[++j] = (N + NW - NNW);
	ctx->context[++j] = (N + NW - NNW + p1 - (Np1 + NWp1 - NNEp1));
	ctx->context[++j] = (N + NW - NNW + p2 - (Np2 + NWp2 - NNEp2));
	ctx->context[++j] = (N + NE - NNE);
	ctx->context[++j] = (N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1));
	ctx->context[++j] = (N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2));
	ctx->context[++j] = (N + NN - NNN);
	ctx->context[++j] = (N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1));
	ctx->context[++j] = (N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2));
	ctx->context[++j] = (W + WW - WWW);
	ctx->context[++j] = (W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1));
	ctx->context[++j] = (W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2));
	ctx->context[++j] = (W + NEE - NE);
	ctx->context[++j] = (W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1));
	ctx->context[++j] = (W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2));
	ctx->context[++j] = (NN + p1 - NNp1);
	ctx->context[++j] = (NN + p2 - NNp2);
	ctx->context[++j] = (NN + W - NNW);
	ctx->context[++j] = (NN + W - NNW + p1 - (NNp1 + Wp1 - NNEp1));
	ctx->context[++j] = (NN + W - NNW + p2 - (NNp2 + Wp2 - NNEp2));
	ctx->context[++j] = (NN + NW - NNNW);
	ctx->context[++j] = (NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1));
	ctx->context[++j] = (NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2));
	ctx->context[++j] = (NN + NE - NNNE);
	ctx->context[++j] = (NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1));
	ctx->context[++j] = (NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2));
	ctx->context[++j] = (NN + NNNN - NNNNNN);
	ctx->context[++j] = (NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1));
	ctx->context[++j] = (NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2));
	ctx->context[++j] = (WW + p1 - WWp1);
	ctx->context[++j] = (WW + p2 - WWp2);
	ctx->context[++j] = (WW + WWWW - WWWWWW);
	ctx->context[++j] = (WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1));
	ctx->context[++j] = (WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2));
	ctx->context[++j] = (N * 2 - NN + p1 - (Np1 * 2 - NNp1));
	ctx->context[++j] = (N * 2 - NN + p2 - (Np2 * 2 - NNp2));
	ctx->context[++j] = (W * 2 - WW + p1 - (Wp1 * 2 - WWp1));
	ctx->context[++j] = (W * 2 - WW + p2 - (Wp2 * 2 - WWp2));
	ctx->context[++j] = (N * 3 - NN * 3 + NNN);
	ctx->context[++j] = clamp4(N * 3 - NN * 3 + NNN, W, NW, N, NE);
	ctx->context[++j] = clamp4(W * 3 - WW * 3 + WWW, W, NW, N, NE);
	ctx->context[++j] = clamp4(N * 2 - NN, W, NW, N, NE);
	ctx->context[++j] = ((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW, W, NW, N, NN)) / 6);
	ctx->context[++j] = ((NNNEEE - 4 * NNEE + 6 * NE + (W * 4 - NW * 6 + NNW * 4 - NNNW)) / 4);
	ctx->context[++j] = (((N + 3 * NW) / 4) * 3 - ((NNW + NNWW) / 2) * 3 + (NNNWW * 3 + NNNWWW) / 4);
	ctx->context[++j] = ((W * 2 + NW) - (WW + 2 * NWW) + NWWW);
	ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;
	ctx->context[++j] = (N + W + 1) >> 1;
	ctx->context[++j] = (NEEEE + NEEEEEE + 1) >> 1;
	ctx->context[++j] = (WWWWWW + WWWW + 1) >> 1;
	ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;
	ctx->context[++j] = N;
	ctx->context[++j] = NN;
	ctx->context[++j] = N + p1 - Np1;
	ctx->context[++j] = N + p2 - Np2;
	ctx->context[++j] = W + p1 - Wp1;
	ctx->context[++j] = W + p2 - Wp2;
	ctx->context[++j] = NW + p1 - NWp1;
	ctx->context[++j] = NW + p2 - NWp2;
	ctx->context[++j] = NE + p1 - NEp1;
	ctx->context[++j] = NE + p2 - NEp2;
	ctx->context[++j] = NN + p1 - NNp1;
	ctx->context[++j] = NN + p2 - NNp2;
	ctx->context[++j] = WW + p1 - WWp1;
	ctx->context[++j] = WW + p2 - WWp2;
	ctx->context[++j] = W + N - NW;
	ctx->context[++j] = W + N - NW + p1 - Wp1 - Np1 + NWp1;
	ctx->context[++j] = W + N - NW + p2 - Wp2 - Np2 + NWp2;
	ctx->context[++j] = W + NE - N;
	ctx->context[++j] = W + NE - N + p1 - Wp1 - NEp1 + Np1;
	ctx->context[++j] = W + NE - N + p2 - Wp2 - NEp2 + Np2;
	ctx->context[++j] = W + NEE - NE;
	ctx->context[++j] = W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1;
	ctx->context[++j] = W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2;
	ctx->context[++j] = N + NN - NNN;
	ctx->context[++j] = N + NN - NNN + p1 - Np1 - NNp1 + NNNp1;
	ctx->context[++j] = N + NN - NNN + p2 - Np2 - NNp2 + NNNp2;
	ctx->context[++j] = N + NE - NNE;
	ctx->context[++j] = N + NE - NNE + p1 - Np1 - NEp1 + NNEp1;
	ctx->context[++j] = N + NE - NNE + p2 - Np2 - NEp2 + NNEp2;
	ctx->context[++j] = N + NW - NNW;
	ctx->context[++j] = N + NW - NNW + p1 - Np1 - NWp1 + NNEp1;
	ctx->context[++j] = N + NW - NNW + p2 - Np2 - NWp2 + NNEp2;
	ctx->context[++j] = NE + NW - NN;
	ctx->context[++j] = NE + NW - NN + p1 - NEp1 - NWp1 + NNp1;
	ctx->context[++j] = NE + NW - NN + p2 - NEp2 - NWp2 + NNp2;
	ctx->context[++j] = NW + W - NWW;
	ctx->context[++j] = NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1;
	ctx->context[++j] = NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2;
	ctx->context[++j] = W * 2 - WW;
	ctx->context[++j] = W * 2 - WW + p1 - Wp1 * 2 + WWp1;
	ctx->context[++j] = W * 2 - WW + p2 - Wp2 * 2 + WWp2;
	ctx->context[++j] = N * 2 - NN;
	ctx->context[++j] = N * 2 - NN + p1 - Np1 * 2 + NNp1;
	ctx->context[++j] = N * 2 - NN + p2 - Np2 * 2 + NNp2;
	ctx->context[++j] = NW * 2 - NNWW;
	ctx->context[++j] = NW * 2 - NNWW + p1 - NWp1 * 2 + NNWWp1;
	ctx->context[++j] = NW * 2 - NNWW + p2 - NWp2 * 2 + NNWWp2;
	ctx->context[++j] = NE * 2 - NNEE;
	ctx->context[++j] = NE * 2 - NNEE + p1 - NEp1 * 2 + NNEEp1;
	ctx->context[++j] = NE * 2 - NNEE + p2 - NEp2 * 2 + NNEEp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p1 - Np1 * 3 + NNp1 * 3 - NNNp1;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p2 - Np2 * 3 + NNp2 * 3 - NNNp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN;
	ctx->context[++j] = (W + NE * 2 - NNE + 1) >> 1;
	ctx->context[++j] = (W + NE * 3 - NNE * 3 + NNNE+1) >> 1;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1 * 2 - NNEp1) / 2;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2 * 2 - NNEp2) / 2;
	ctx->context[++j] = NNE + NE - NNNE;
	ctx->context[++j] = NNE + W - NN;
	ctx->context[++j] = NNW + W - NNWW;
#endif
#undef LOAD
	for(int k=0;k<T39_NMAPS;++k)
	{
#ifdef T39_USE_KALMAN
		kalman_predict(ctx->kalman+k, ctx->context[k]);
		ctx->context[k]=(int)round(CLAMP(-128, ctx->kalman[k].Uhat, 127));
#endif
		ctx->context[k]+=128;
#ifdef T39_USE_ARRAYS
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
#else
		ctx->context[k]<<=8;
#endif
	}
	
#ifdef T39_PROB_TWEAK
	ctx->bitsprocessed+=!kc;
#endif
}
int t39_ctx_map_context(int *context, int kp, int workidx)//replacement for context[kp]
{
	return context[kp];

	//static const int rep[]={ 0,  0, 10, 10, 10, 10, 10,  0,     10, 10, 10,  0,  0, 10, 10,  0,      8,  8, 11, 11, 11, 11, 11,  0};
	//static const int sub[]={ 6,  6,  6,  6,  2,  0,  0, 10,     10,  6,  6,  2,  2, 13,  6, 10,      2,  1,  2, 13, 13, 12,  0, 12};
	//return context[kp==rep[workidx]?sub[workidx]:kp];
}
void t39_ctx_estimate_p0(T39Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
#ifndef T39_USE_ARRAYS
	RBNodeHandle *hnode;
#endif
	T39Node *node;
	for(int kp=0;kp<T39_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		int context=t39_ctx_map_context(ctx->context, kp, kc);
#ifdef T39_USE_ARRAYS
		ArrayHandle map=ctx->maps[workidx][kp];
		node=ctx->node[kp]=(T39Node*)array_at(&map, context);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T39_DISABLE_REC
		for(;k2<T39_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
#else
		hnode=map_insert(ctx->maps[workidx]+kp, &context, ctx->found+kp);
		node=ctx->node[kp]=(T39Node*)hnode[0]->data;
		if(ctx->found[kp])
		{
			sum=node->n[0]+node->n[1];
			ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
			++k2;
			for(;k2<T39_N_REC_ESTIMATORS+1;++k2)
				ctx->p0arr[p0idx+k2]=node->rec[k2-1];
			p0idx+=k2;
		}
		else
		{
			for(;k2<T39_N_REC_ESTIMATORS+1;++k2)
				ctx->p0arr[p0idx+k2]=0x8000;
			//{
			//	int bestbit=ctx->context[k2]>>kb&1;
			//	ctx->p0arr[p0idx+k2]=0x8000+0x0000*((bestbit<<1)-1);
			//}
			p0idx+=k2;
		}
#endif
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T39_NESTIMATORS;++k)
	{
#ifdef T39_DISABLE_COUNTER
		if(k%(T39_N_REC_ESTIMATORS+1))//
#endif
		{
			sum+=(long long)ctx->p0arr[k]*wk[k];
			ctx->wsum+=wk[k];
		}
	}
	//ctx->p0=ctx->wsum?(int)((sum+(ctx->wsum>>1))/ctx->wsum):0x8000;//same CR
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);


	//if(ctx->p0!=0x8000)//
	//	printf("");
	
#ifdef T39_PROB_TWEAK
	//if(workidx==0)//
	//	printf("");

	int bitsprocessed=ctx->bitsprocessed-1;
	int prevproberror=bitsprocessed?(int)(ctx->proberrors[workidx]/bitsprocessed):0;
	double success=bitsprocessed?(double)ctx->hits[workidx]/bitsprocessed:0;
	if(success>0.5&&ctx->p0!=0x8000)
	//if(prevproberror&&prevproberror<0x8000&&ctx->p0!=0x8000)
	{
		double invCR=ctx->csizes[workidx]/bitsprocessed;
		double success2=pow((success-0.5)*2, invCR)*0.5+0.5;
		int p_MPS=(int)(success2*0x10000);
		//int p_MPS=(int)(((long long)ctx->hits[workidx]<<16)/bitsprocessed);
		//int p_MPS=0xFFFF;
		//int p_MPS=0x10000-prevproberror;
		int p_MPS0=ctx->p0<0x8000?0x10000-ctx->p0:ctx->p0;
		if(p_MPS<p_MPS0)
			p_MPS=p_MPS0;
		//else
		//	p_MPS=0xFFFF;
		ctx->p0rev=ctx->p0<0x8000?0x10000-p_MPS:p_MPS;
		ctx->p0rev=CLAMP(1, ctx->p0rev, 0xFFFF);

		//ctx->p0rev=ctx->p0rev+((ctx->p0-ctx->p0rev)>>0);
		//ctx->p0rev=(ctx->p0rev+ctx->p0+(ctx->p0>0x8000))>>1;
	}
	else
		ctx->p0rev=ctx->p0;

	//int p_MPS=ctx->hits[24]-1?(int)(((long long)ctx->hits[workidx]*3<<16)/(ctx->hits[24]-1)):(ctx->p0<0x8000?0x10000-ctx->p0:ctx->p0);
	//if(p_MPS>0x8000)
	//{
	//	ctx->p0rev=ctx->p0<0x8000?0x10000-p_MPS:p_MPS;
	//	ctx->p0rev=CLAMP(1, ctx->p0rev, 0xFFFF);
	//}
	//else
	//	ctx->p0rev=ctx->p0;
#else
	ctx->p0rev=ctx->p0;
#endif
}
void t39_ctx_update(T39Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
	
#ifdef T39_PROB_TWEAK
	int prob=(bit?0x10000-ctx->p0:ctx->p0);
	double p=(double)prob/0x10000;
	double bitsize=-log2(p);
	ctx->csizes[workidx]+=bitsize;
#endif

#ifdef T39_PRINT_ESTIMATOR_CR
	for(int k=0;k<T39_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T39_NESTIMATORS*workidx+k]+=bitsize;
		}
	}
#endif
	//bwd
	int *wk=ctx->weights[workidx];
	if(ctx->p0_0>=1&&ctx->p0_0<=0xFFFF)
	{
		int p_bit=bit?0x10000-ctx->p0:ctx->p0;
		long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
		dL_dp0^=-bit;
		dL_dp0+=bit;
		for(int k=0;k<T39_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
#ifndef T39_DISABLE_REC
	static const int shifts[]=
	{
		//7, 7, 7, 6, 5, 0, 0, 8,//T39_N_REC_ESTIMATORS 1
		//8, 7, 7, 6, 6, 6, 7, 7,
		//7, 7, 7, 5, 4, 4, 2, 7,

		5, 5, 4, 2, 1, 0, 0, 5,//T39_N_REC_ESTIMATORS 6
		7, 5, 5, 4, 4, 5, 4, 6,
		5, 5, 5, 2, 1, 1, 1, 5,

		//6, 5, 5, 4, 2, 0, 0, 6,
		//9, 7, 6, 5, 4, 5, 4, 6,
		//6, 5, 5, 4, 4, 3, 1, 6,
	};
	//static const int bitrating[]=//higher is more compressible
	//{
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//	1, 1, 1, 1, 2, 3, 4, 1,
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//};
#endif
	T39Node *node;
	for(int kp=0;kp<T39_NMAPS;++kp)
	{
		node=ctx->node[kp];
#ifndef T39_USE_ARRAYS
		if(ctx->found[kp])
		{
#endif
			++node->n[bit];
#ifndef T39_DISABLE_REC
			//static const int lgdens[]={0, 1, 2, 4, 8, 14};
			for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
			{
				//int lgden=k+2;
				//int lgden=k+shifts[workidx];
				//int lgden=k+1;
				//int lgden=(k+1)<<1;//X
				//int lgden=k+((4-bitrating[workidx])<<1);
				//int lgden=k+4-bitrating[workidx];
				//int lgden=k+3;
				int lgden=k;
				//int lgden=lgdens[k];
				//int lgden=((k+1)<<1)-1;
				int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
				node->rec[k]=CLAMP(1, temp, 0xFFFF);
			}
#endif
#ifndef T39_USE_ARRAYS
		}
		else
		{
			int context=t39_ctx_map_context(ctx->context, kp, kc);
			node->key=context;
			node->n[0]=1;
			node->n[1]=1;
			for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
				node->rec[k]=0x8000;
			++ctx->nnodes;
		}
		ctx->context[kp]|=bit<<kb;
#else
		ctx->context[kp]|=bit<<(8+7-kb);//append bits in reverse
		//ctx->context[kp]<<=1;
		//ctx->context[kp]|=bit;
#endif
	}
#ifdef T39_PROB_TWEAK
	ctx->proberrors[workidx]+=abs((bit<<16)-ctx->p0);
	ctx->hits[workidx]+=bit==(ctx->p0<0x8000);
	//if(ctx->hits[workidx]<0)
	//	printf("");
	//if(workidx==15)
	//	printf("");
	//ctx->hits[workidx]+=bit==(ctx->p0<0x8000);
#endif
}
int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);
	T39Ctx *t39_ctx=t39_ctx_init();
	if(!buf2||!ebuf||!t39_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	memset(ebuf, 0, (size_t)res<<2);
#ifdef T39_APPLY_SPATIAL
	apply_transforms_fwd(buf2, iw, ih);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//buffer is signed
#else
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_YCbCr_R_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//X  the buffer is signed
#endif

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	float csizes[24]={0};
	//int hits[24]={0};
	
#if 0
	int nx=(iw+T39_BLOCKSIZE_X-1)/T39_BLOCKSIZE_X,
		ny=(ih+T39_BLOCKSIZE_Y-1)/T39_BLOCKSIZE_Y;
	for(int by=0, blockNo=0;by<ny;++by)
	{
		int y1=by*T39_BLOCKSIZE_Y, y2=y1+T39_BLOCKSIZE_Y;
		if(y2>ih)y2=ih;

		for(int bx=0;bx<nx;++bx, ++blockNo)
		{
			int x1=bx*T39_BLOCKSIZE_X, x2=x1+T39_BLOCKSIZE_X;
			if(x2>iw)x2=iw;
			
#ifdef T39_PROB_TWEAK
			t39_ctx_reset(t39_ctx, 0);
#else
			t39_ctx_reset(t39_ctx, 1);
#endif
			//t39_ctx_reset(t39_ctx, !(blockNo&63));

			for(int ky=y1;ky<y2;++ky)
			{
				for(int kx=x1;kx<x2;++kx)
				{
					//if(kx==5&&ky==5)//
					//	printf("");

					for(int kc=0;kc<3;++kc)
					{
						t39_ctx_get_context(t39_ctx, (char*)buf2, iw, ih, kc, kx, ky, x1, x2, y1, y2);
						for(int kb=7;kb>=0;--kb)//MSB -> LSB
						{
							t39_ctx_estimate_p0(t39_ctx, kc, kb);

							int bit=(buf2[(iw*ky+kx)<<2|kc]+128)>>kb&1;//signed -> unsigned
							abac_enc(&ec, t39_ctx->p0rev, bit);

							int prob_ideal=0x10000&-bit;//
							float prob_error=(float)abs(prob_ideal-t39_ctx->p0rev)*(1.f/0x10000);
							//hits[kc<<3|kb]+=prob_error;
							float bitsize=-log2f(prob_error);
							//if(bitsize<0)
							//	LOG_ERROR("Bitsize negative");
							csizes[kc<<3|kb]+=bitsize;//

							t39_ctx_update(t39_ctx, kc, kb, bit);
						}
					}
				}
			}
		}

		if(loud)
		{
			static float csize1=0;
			if(!by)
				csize1=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k];
			csize/=8;

			//if(iw*(y2-y1)*3/(csize-csize1)<0)
			//	LOG_ERROR("Size estimation error");

			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f%c", by+1, ny, 100.*(by+1)/ny, iw*(y2+1)*3/csize, iw*(y2-y1)*3/(csize-csize1), loud==2?'\n':'\r');
			csize1=csize;
		}
	}
#endif
#if 1
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t39_ctx_get_context(t39_ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t39_ctx_estimate_p0(t39_ctx, kc, kb);
#ifdef T39_UNSIGNED_BITS
					int bit=(buf2[idx]+128)>>kb&1;
#else
					int bit=buf2[idx]>>kb&1;
#endif
					ac_enc_bin(&ec, t39_ctx->p0, bit);
					
					int prob=bit?0x10000-t39_ctx->p0:t39_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t39_ctx_update(t39_ctx, kc, kb, bit);
				}
				ebuf[idx]=buf2[idx]-t39_ctx->pred14;
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f\r", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev));
			//printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f%c", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev), loud==2?'\n':'\r');
			csize_prev=csize;
		}
		//if(!((ky+1)&127))
		//	t39_ctx_reset(&t39_ctx, 0);
	}
#endif
	ac_enc_flush(&ec);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t39_ctx->nnodes*sizeof(T39Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
#if 0
		double csize=0;
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld    CR %lf    WH %d*%d  bitplane %g\n", list.nobj, 3.*iw*ih/list.nobj, iw, ih, iw*ih/8.);
		printf("\n");
#endif
		
		float chsizes[4]={0};
		//printf("\t\tC0\t\t\t\tC1\t\t\t\tC2\n\n");
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=7;kb>=0;--kb)
		{
			printf("B%d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc<<3|kb;
				float size=csizes[idx];
				//printf("       %12.3f %12.2f", iw*ih/size, hits[idx]);
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);

#ifdef T39_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t39_ctx->csizes_est+T39_NESTIMATORS*kb;
				for(int ke=1;ke<T39_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T39_NESTIMATORS;++ke)
			{
				float *sizes=t39_ctx->csizes_est+ke;
#ifndef T39_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T39_N_REC_ESTIMATORS+1), ke%(T39_N_REC_ESTIMATORS+1));
#else
				printf("E%3d ", ke);
#endif
				for(int kb=0;kb<24;++kb)
				{
					char c;
					if(ke==minidx[kb])
						c='*';
					else if(ke==maxidx[kb])
						c='L';
					else
						c=' ';
					printf("%8.2f %c", iw*ih/sizes[T39_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T39_NESTIMATORS*kb]/t39_ctx->csizes_est[T39_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T39_DISABLE_REC
				if(!((ke+1)%(T39_N_REC_ESTIMATORS+1)))
#else
				if(!((ke+1)%8))
#endif
				{
					printf("\n");
					printf("\t\t*         **        ***       ****      ****      ****      ****      *             *         *         *         *         **        ***       ****      *             *         **        ***       ****      ****      ****      ****      *\n");
					printf("\t\t0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7\n");
				}
			}
		}
#endif
		//t39_explore(t39_ctx, buf2, iw, ih);
	}
	t39_ctx_clear(&t39_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	//int debug_index=0;
	char *ebuf=(char*)malloc((size_t)res<<2);
	T39Ctx *t39_ctx=t39_ctx_init();
	if(!ebuf||!t39_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t39_ctx_init(t39_ctx);
	
#if 0
	int nx=(iw+T39_BLOCKSIZE_X-1)/T39_BLOCKSIZE_X,
		ny=(ih+T39_BLOCKSIZE_Y-1)/T39_BLOCKSIZE_Y;
	for(int by=0, blockNo=0;by<ny;++by)
	{
		int y1=by*T39_BLOCKSIZE_Y, y2=y1+T39_BLOCKSIZE_Y;
		if(y2>ih)y2=ih;
		for(int bx=0;bx<nx;++bx, ++blockNo)
		{
			int x1=bx*T39_BLOCKSIZE_X, x2=x1+T39_BLOCKSIZE_X;
			if(x2>iw)x2=iw;
			
#ifdef T39_PROB_TWEAK
			t39_ctx_reset(t39_ctx, 0);
#else
			t39_ctx_reset(t39_ctx, 1);
#endif

			for(int ky=y1;ky<y2;++ky)
			{
				for(int kx=x1;kx<x2;++kx)
				{
					for(int kc=0;kc<3;++kc)
					{
						t39_ctx_get_context(t39_ctx, (char*)buf, iw, ih, kc, kx, ky, x1, x2, y1, y2);
						for(int kb=7;kb>=0;--kb)//MSB -> LSB
						{
							t39_ctx_estimate_p0(t39_ctx, kc, kb);
					
							int bit=abac_dec(&ec, t39_ctx->p0rev);
							buf[(iw*ky+kx)<<2|kc]|=bit<<kb;

							t39_ctx_update(t39_ctx, kc, kb, bit);
						}
						buf[(iw*ky+kx)<<2|kc]+=128;//unsigned -> signed
					}
				}
			}
		}

		if(loud)
			printf("%5d/%5d  %6.2lf%%%c", by+1, ny, 100.*(by+1)/ny, loud==2?'\n':'\r');
	}
#endif
#if 1
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t39_ctx_get_context(t39_ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t39_ctx_estimate_p0(t39_ctx, kc, kb);
					
					int bit=ac_dec_bin(&ec, t39_ctx->p0);
					buf[idx]|=bit<<kb;

					t39_ctx_update(t39_ctx, kc, kb, bit);
				}
#ifdef T39_UNSIGNED_BITS
				buf[idx]+=128;//unsigned -> signed
#endif
				ebuf[idx]=buf[idx]-t39_ctx->pred14;
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
		//if(!((ky+1)&127))
		//	t39_ctx_reset(&t39_ctx, 0);
	}
#endif
	t39_ctx_clear(&t39_ctx);
	
#ifdef T39_APPLY_SPATIAL
	addbuf(buf, iw, ih, 3, 4, 128);//buffer is signed
	apply_transforms_inv(buf, iw, ih);
#else
	//addbuf(buf, iw, ih, 3, 4, 128);//X  the buffer is signed
	colortransform_YCbCr_R_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}


//typedef enum MATestTypeEnum
//{
//	//MA_RCT_LUMA=1,
//	MA_RCT_CHROMA=1,
//	MA_PRED=2,
//	//MA_PRED_LUMA=2,
//	//MA_PRED_CHROMA=4,
//} MATestType;
typedef enum RCTTypeEnum
{
	RCT_SubG,
	RCT_JPEG2000,
	RCT_YCoCg_R,
	RCT_YCbCr_R,
	RCT_YCbCr_R_v2,
	RCT_YCbCr_R_v3,
	RCT_YCbCr_R_v4,
	RCT_YCbCr_R_v5,
	RCT_YCbCr_R_v6,
	RCT_DCT,
} RCTType;
//void print_ma_test(int testtype)
//{
//#define PRINT_TESTTYPE(LABEL) if(testtype&LABEL)printf("[ v] " #LABEL "\n"); else printf("[X ] " #LABEL "\n");
//	//PRINT_TESTTYPE(MA_RCT_LUMA)
//	PRINT_TESTTYPE(MA_RCT_CHROMA)
//	PRINT_TESTTYPE(MA_PRED)
//	//PRINT_TESTTYPE(MA_PRED_LUMA)
//	//PRINT_TESTTYPE(MA_PRED_CHROMA)
//#undef  PRINT_TESTTYPE
//}
size_t ma_test(const unsigned char *src, int iw, int ih, int enable_RCT_MA, int RCTtype, int enable_rounding, int loud)
{
	int res=iw*ih;
	short *buf1=(short*)malloc((size_t)res*sizeof(short[4]));
	short *buf2=(short*)malloc((size_t)res*sizeof(short[4]));
	int *hist=(int*)malloc(1024*sizeof(int));
	if(!buf1||!buf2||!hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}

	if(loud)
	{
		printf("%s\n", enable_RCT_MA?"RCT MA ON":"RCT MA OFF");
		//printf("%d:\n", testtype);
		//print_ma_test(testtype);
	}

	int nbits_ch[]={8, 9, 9};//8 + number of subtractions

	//luma retains range after YCbCr-RCT

	if(RCTtype==RCT_DCT)
		++nbits_ch[0];
	//if(testtype&MA_RCT_CHROMA)
	if(enable_RCT_MA)
	{
		if(RCTtype==RCT_DCT)
			--nbits_ch[0];
		--nbits_ch[1];
		--nbits_ch[2];
	}
	//if(testtype&MA_PRED)
	//{
	//	--nbits_ch[0];
	//	--nbits_ch[1];
	//	--nbits_ch[2];
	//}
	int vmin[6]={0}, vmax[6]={0};
#define UPDATE_MINMAX(IDX, COMP) if(vmin[IDX]>COMP)vmin[IDX]=COMP; if(vmax[IDX]<COMP)vmax[IDX]=COMP;
	for(int k=0;k<res;++k)
	{
		int r=src[k<<2|0]-128, g=src[k<<2|1]-128, b=src[k<<2|2]-128;
		
		//DCT		X  all channels inflate by 1 bit
#if 0
		r-=(g+b+1)>>1;
		b+=g+((r+1)>>1);
		g-=((b+(r>>3)+1)>>1);
		r-=(g+g+g+2)>>2;
#endif

		//YCbCr-R
#if 0
		r-=g;		//diff(r, g)            [ 1      -1      0  ].RGB	Cr
		g+=r>>1;
		b-=g;		//diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB	Cb
		g+=b>>1;	//av(b, av(r, g))       [ 1/4     1/4    1/2].RGB	Y
#endif
		switch(RCTtype)
		{
		case RCT_SubG:
			r-=g;
			b-=g;
			break;
		case RCT_JPEG2000:
			r-=g;
			b-=g;
			g+=(r+b+(enable_rounding<<1))>>2;
			break;
		case RCT_YCoCg_R:
			r-=b;				//co = r-b			diff(r, b)
			b+=(r+enable_rounding)>>1;	//(r+b)/2
			g-=b;				//cg = g-(r+b)/2		diff(g, av(r, b))
			b+=(g+enable_rounding)>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4		av(g, av(r, b))
			{
				char temp;
				SWAPVAR(g, b, temp);//g must contain luma which is always 8-bit
			}
			break;
		case RCT_YCbCr_R:
			r-=g;				//diff(r, g)            [ 1      -1      0  ].RGB	Cr
			g+=(r+enable_rounding)>>1;
			b-=g;				//diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB	Cb
			g+=(b+enable_rounding)>>1;	//av(b, av(r, g))       [ 1/4     1/4    1/2].RGB	Y
			break;
		case RCT_YCbCr_R_v2:
			r-=g;					//Cr =	[1	-1	0].RGB
			g+=(r+enable_rounding)>>1;		//	[1/2	1/2	0]
			b-=g;					//Cb =	[-1/2	-1/2	1]
			g+=(2*b-r+(enable_rounding<<2))>>3;	//Y  =	[1/4	1/2	1/4]	v2
			break;
		case RCT_YCbCr_R_v3:
			r-=g;					//Cr =	[1	-1	0].RGB
			g+=(r+enable_rounding)>>1;		//	[1/2	1/2	0]
			b-=g;					//Cb =	[-1/2	-1/2	1]
			g+=(2*b+r+(enable_rounding<<2))>>3;	//Y  =	[1/2	1/4	1/4]	v3
			break;
		case RCT_YCbCr_R_v4:
			r-=g;				//Cr =	[1	-1	0].RGB
			g+=(r+enable_rounding)>>1;	//	[1/2	1/2	0]
			b-=g;				//Cb =	[-1/2	-1/2	1]
			g+=(b+enable_rounding)/3;	//Y  =	[1/3	1/3	1/3]	v4
			break;
		case RCT_YCbCr_R_v5:
			r-=g;					//Cr =	[1	-1	0].RGB
			g+=(r+enable_rounding)>>1;		//	[1/2	1/2	0]
			b-=g;					//Cb =	[-1/2	-1/2	1]
			g+=(b*6+(enable_rounding<<3))>>4;	//Y  =	[5/16	5/16	6/16]	v5
			break;
		case RCT_YCbCr_R_v6:
			r-=g;					//Cr =	[1	-1	0].RGB
			g+=(r+enable_rounding)>>1;		//	[1/2	1/2	0]
			b-=g;					//Cb =	[-1/2	-1/2	1]
			g+=(b*14+(enable_rounding<<4))>>5;	//Y  =	[9/32	9/32	14/32]	v6
			break;
		case RCT_DCT:
			r-=(g+b+enable_rounding)>>1;
			b+=g+((r+enable_rounding)>>1);
			g-=((b+(r>>3)+enable_rounding)>>1);
			r-=(g+g+g+(enable_rounding<<1))>>2;
			{
				char temp;
				SWAPVAR(g, b, temp);//g must contain luma which is always 8-bit
			}
			break;
		}
		
		//if(g<-128||g>127)//no hit
		//	LOG_ERROR("");

		buf1[k<<2|0]=g;//Y	8-bit because the average (luma) can't spill out of input domain
		//if(testtype&MA_RCT_CHROMA)
		if(enable_RCT_MA)
		{
			buf1[k<<2|1]=((b+128)&0xFF)-128;
			buf1[k<<2|2]=((r+128)&0xFF)-128;
		}
		else
		{
			buf1[k<<2|1]=b;//Cb	9-bit
			buf1[k<<2|2]=r;//Cr	9-bit
		}
		buf1[k<<2|3]=0xFFFF;
		UPDATE_MINMAX(0, g)
		UPDATE_MINMAX(1, b)
		UPDATE_MINMAX(2, r)
	}
	if(loud)
		printf("RCT  YCbCr [%d %d]  [%d %d]  [%d %d]\n", vmin[0], vmax[0], vmin[1], vmax[1], vmin[2], vmax[2]);
	memset(buf2, 0, res*sizeof(short[4]));
	for(int kc=0;kc<3;++kc)//clamped gradient
	{
		int nlevels=1<<nbits_ch[kc];
		int half=nlevels>>1;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?buf1[(iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				short
					N =LOAD( 0, -1),
					W =LOAD(-1,  0),
					NW=LOAD(-1, -1);
#undef  LOAD
				int pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);

				int idx=(iw*ky+kx)<<2|kc;
				buf2[idx]=((buf1[idx]-pred+half)&(nlevels-1))-half;
				UPDATE_MINMAX(kc+3, buf2[idx])
			}
		}
	}
	if(loud)
		printf("pred YCbCr [%d %d]  [%d %d]  [%d %d]\n", vmin[3], vmax[3], vmin[4], vmax[4], vmin[5], vmax[5]);
	double invCR[3];
	size_t csize_bytes[3]={0};
	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<nbits_ch[kc];
		int half=nlevels>>1;
		memset(hist, 0, nlevels*sizeof(int));
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx)
				++hist[buf2[idx<<2|kc]+half];
		}
		double entropy=0;
		for(int sym=0;sym<nlevels;++sym)//measure entropy using Shannon law
		{
			int freq=hist[sym];
			if(freq)
			{
				double p=(double)freq/res;
				entropy-=p*log2(p);
			}
		}
		invCR[kc]=entropy/8;//invCR = entropy/lg(nlevels) * lg(nlevels)/8	invCR gets amplified by lg(nlevels)/8

		for(int k=0, sum=0;k<nlevels;++k)//accumulate CDF
		{
			int freq=hist[k];
			hist[k]=(int)((long long)sum*(0x10000-nlevels)/res)+k;
			sum+=freq;
		}
		unsigned state=0x10000;
		for(int k=0;k<res;++k)//measure compressed size using ANS
		{
			unsigned short sym=buf2[k<<2|kc]+half;
			int cdf=hist[sym], freq=(sym==nlevels-1?0x10000:hist[sym+1])-cdf;
			if(state>((unsigned)freq<<16))//renorm
			{
				csize_bytes[kc]+=2;
				state>>=16;
			}
			state=state/freq<<16|(cdf+state%freq);//update
		}
		csize_bytes[kc]+=4;

	}
	if(loud)
	{
		printf(
			"TYCbCr %lf %lf %lf %lf\n",
			3/(invCR[0]+invCR[1]+invCR[2]),
			1/invCR[0],
			1/invCR[1],
			1/invCR[2]
		);
		printf(
			"TYCbCr %lf %lf %lf %lf\n",
			(double)(res*3)/(csize_bytes[0]+csize_bytes[1]+csize_bytes[2]),
			(double)res/csize_bytes[0],
			(double)res/csize_bytes[1],
			(double)res/csize_bytes[2]
		);
		printf(
			"TYCbCr %lld %lld %lld %lld\n",
			csize_bytes[0]+csize_bytes[1]+csize_bytes[2],
			csize_bytes[0],
			csize_bytes[1],
			csize_bytes[2]
		);
		printf("Done.\n\n");
	}
	
#undef  UPDATE_MINMAX
	free(buf1);
	free(buf2);
	free(hist);
	return csize_bytes[0]+csize_bytes[1]+csize_bytes[2];
}