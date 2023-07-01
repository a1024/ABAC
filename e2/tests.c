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
#define LOAD(PX) ((PX)<<MAGBITS)
#define STORE(VAL) (int)(((VAL)+(1<<(MAGBITS-1))-1)>>MAGBITS)
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
#define LOAD(PX) MUL(PX+0.5f, 1.f/256)
#define STORE(VAL) (char)(MUL(VAL, 256)-0.5f)
//#define LEARNING_RATE(G) ((G)*0.001f)
#define ABS fabsf
#define FROMFLOAT(X) X
#endif
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
	double t_start=time_ms();
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
		timedelta2str(0, 0, time_ms()-t_start);
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

	if(!BETWEEN(1, param->gwidth, t25_limits.gwidth)||!BETWEEN(0, param->mleft, t25_limits.mleft)||!BETWEEN(0, param->mtop, t25_limits.mtop)||!BETWEEN(0, param->mright, t25_limits.mright)||!BETWEEN(0, param->alpha, t25_limits.alpha)||!BETWEEN(0, param->maxinc, t25_limits.maxinc))
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
	if(!BETWEEN(1, param->gwidth, t25_limits.gwidth)||!BETWEEN(0, param->mleft, t25_limits.mleft)||!BETWEEN(0, param->mtop, t25_limits.mtop)||!BETWEEN(0, param->mright, t25_limits.mright)||!BETWEEN(0, param->alpha, t25_limits.alpha)||!BETWEEN(0, param->maxinc, t25_limits.maxinc))
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

	for(int ks=0;ks<COUNTOF(steps);++ks)
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
	for(int ks=0;ks<COUNTOF(steps);++ks)
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
		for(int ks=0;ks<COUNTOF(steps);++ks)
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
	double t_start=time_ms();
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
			ANSEncContext ctx;
			ans_enc_init(&ctx, CDF2, &list);
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
							ans_enc(&ctx, buf2[(iw*ky+kx)<<2|kc], kc);
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
			ans_enc_flush(&ctx);
			//ans_enc_flush(&list, state);
			//dlist_push_back(&list, &state, 4);
		}
		else
		{
			ACEncContext ctx;
			ac_enc_init(&ctx, CDF2, &list);
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
							ac_enc(&ctx, buf2[(iw*ky+kx)<<2|kc]);
						}
					}
				}
			}
			ac_enc_flush(&ctx);
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
		timedelta2str(0, 0, time_ms()-t_start);
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
			ANSDecContext ctx;
			ans_dec_init(&ctx, CDF2, srcstart, srcend);
			//ans_dec_init(&ctx, CDF2, data+(kc?ansbookmarks[kc-1]:overhead), data+ansbookmarks[kc]);
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
							buf[(iw*ky+kx)<<2|kc]=ans_dec(&ctx, kc);
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
			ACDecContext ctx;
			ac_dec_init(&ctx, CDF2, srcstart, srcend);

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
							buf[(iw*ky+kx)<<2|kc]=ac_dec(&ctx);
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
	double t_start=time_ms();
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
			ANSEncContext ctx;
			ans_enc_init(&ctx, CDF2, &list);

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

						ans_enc(&ctx, buf2[(iw*ky+kx)<<2|kc], kc);
					}
				}
			}
			ans_enc_flush(&ctx);
		}
		else
		{
			ACEncContext ctx;
			ac_enc_init(&ctx, CDF2, &list);

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
					
						ac_enc(&ctx, buf2[(iw*ky+kx)<<2|kc]);
					}
				}
			}
			ac_enc_flush(&ctx);
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
		timedelta2str(0, 0, time_ms()-t_start);
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
			ANSDecContext ctx;
			ans_dec_init(&ctx, CDF2, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			//unsigned state_lo=0, state_hi=0xFFFFFFFF, code, cache;
			//int nbits=32;
			//srcptr=kc?data+bookmarks[kc-1]:data+overhead;
			//srcend=data+bookmarks[kc];

			//if(ctx.srcend-ctx.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ctx.code, ctx.srcptr, 4);
			//ctx.srcptr+=4;
			//
			//if(ctx.srcend-ctx.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ctx.cache, ctx.srcptr, 4);
			//ctx.srcptr+=4;

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
						buf[(iw*ky+kx)<<2|kc]=ans_dec(&ctx, kc);
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
			ACDecContext ctx;
			ac_dec_init(&ctx, CDF2, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					int success=t26_prepblock(buf, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);
					if(!success)
						LOG_ERROR("t26_prepblock error");
					for(int kx=x1;kx<x2;++kx)//for each pixel
						buf[(iw*ky+kx)<<2|kc]=ac_dec(&ctx);
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
	double t_start=time_ms();
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
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
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
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);

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
					abac_enc(&ctx, p0, bit);
					
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
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
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
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);

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
							abac_enc(&ctx, p0, bit);
					
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
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
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
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	double *ptree=(double*)malloc(255LL*3*2*sizeof(double));
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
	double inc=1;

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
				double *tree=ptree+(255*kc<<1);
				int step=1;
				unsigned char sym=0;
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kc==1&&kb==5)//
					//	printf("");

					double *c=tree+(sym<<1);
					int p0=(int)(c[0]+c[1]?c[0]*0x10000/(c[0]+c[1]):0x8000);
					//int p0=c[0];
					//int p0=c[0]+c[1]?(int)(((long long)c[0]<<16)/(c[0]+c[1])):0x8000;
					p0=CLAMP(1, p0, 0xFFFF);
					int bit=ptr[kc]>>kb&1;
					abac_enc(&ctx, p0, bit);

					c[bit]+=inc;//CR 2.009992
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
	abac_enc_flush(&ctx);
	
	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		double csize=0;
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
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
		timedelta2str(0, 0, time_ms()-t_start);
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
	double t_start=time_ms();
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
		timedelta2str(0, 0, time_ms()-t_start);
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