#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
#ifdef ALLOW_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include<CL/cl.h>
#pragma comment(lib, "OpenCL.lib")
#endif
//#include<immintrin.h>
static const char file[]=__FILE__;


//	#define FIXEDPREC


//Xoroshiro128+ 1.0 by David Blackman and Sebastiano Vigna		https://prng.di.unimi.it/xoroshiro128plus.c
unsigned long long xoroshiro128_state[2]={0xDF900294D8F554A5, 0x170865DF4B3201FC};
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
//int hamming_weight(unsigned long long x)
//{
//	x=(x&0x5555555555555555)+(x>> 1&0x5555555555555555);
//	x=(x&0x3333333333333333)+(x>> 2&0x3333333333333333);
//	x=(x&0x0F0F0F0F0F0F0F0F)+(x>> 4&0x0F0F0F0F0F0F0F0F);
//	x=(x&0x00FF00FF00FF00FF)+(x>> 8&0x00FF00FF00FF00FF);
//	x=(x&0x0000FFFF0000FFFF)+(x>>16&0x0000FFFF0000FFFF);
//	x=(x&0x00000000FFFFFFFF)+(x>>32&0x00000000FFFFFFFF);
//	return (int)x;
//}

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

#define LEARNING_RATE(G, LR) MUL(G, FROMFLOAT(LR))
#define SGN(X) (DataType)((X)<0?-ONE:((X)>0?ONE:0))
static void initialize(DataType *w, int count, DataType sqrt_fan_in)
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
static void vec_scalar(DataType *dst, DataType *vec, DataType scalar, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=MUL(vec[k], scalar);
}
static void vec_ew(DataType *dst, const DataType *v1, const DataType *v2, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=MUL(v1[k], v2[k]);
}
static void vec_outer(DataType *dst, const DataType *left, const DataType *right, int lh, int rw)
{
	for(int ky=0;ky<lh;++ky)
	{
		for(int kx=0;kx<rw;++kx)
			dst[rw*ky+kx]=MUL(left[ky], right[kx]);
	}
}
static void matmul(DataType *dst, const DataType *m1, const DataType *m2, int h1, int w1h2, int w2)
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


//#define NFEATURES 20
//#define NF0 (NFEATURES<<1)
#define NF0 94
#define NF1 32
#define NF2 32
#define NF3 32
typedef struct TempsStruct
{
	DataType
		nb[NF0],
		net1[NF1], x1[NF1],
		net2[NF2], x2[NF2],
		net3[NF3], x3[NF3],
		pred, expr, loss;
} Temps;
typedef struct BwdTempsStruct
{
	DataType dL_dp, dL_dn3[NF3], dL_dn2[NF2], actdash_n2[NF1], dL_dn1[NF1], actdash_n1[NF1];
} BwdTemps;
typedef struct ParamsStruct
{
	DataType
		weight1[NF1*NF0], bias1[NF1],
		weight2[NF2*NF1], bias2[NF2],
		weight3[NF3*NF2], bias3[NF3],
		weight4[1*NF3];
} Params;
//static DataType gradient(DataType T, DataType L, DataType TL)
//{
//	DataType vmin=L<T?L:T, vmax=L>T?L:T;
//	DataType g=T+L-TL;
//	g=CLAMP(vmin, g, vmax);
//	return g;
//}
static DataType clamp4(DataType p, DataType a, DataType b, DataType c, DataType d)
{
	DataType vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmin>c)vmin=c;
	if(vmin>d)vmin=d;
	if(vmax<b)vmax=b;
	if(vmax<c)vmax=c;
	if(vmax<d)vmax=d;
	p=CLAMP(vmin, p, vmax);
	return p;
}
static DataType clip(DataType x)
{
#ifdef FIXED_PREC
	x=CLAMP(-0x80000, x, 0x7F000);
#else
	x=CLAMP(-1, x, 1);
#endif
	return x;
}
#if 0
static void get_neighbors(const char *buf, int iw, int ih, int x, int y, DataType *ctx)
{
	const int CW=7;
	int idx=iw*y+x, w3=iw*3, w2=iw*2;
	DataType nb[]=
	{
		x-3>=0&&y-3>=0?LOAD(buf[idx-w3-3]):0,
		x-2>=0&&y-3>=0?LOAD(buf[idx-w3-2]):0,
		x-1>=0&&y-3>=0?LOAD(buf[idx-w3-1]):0,
		        y-3>=0?LOAD(buf[idx-w3  ]):0,
		x+1<iw&&y-3>=0?LOAD(buf[idx-w3+1]):0,
		x+2<iw&&y-3>=0?LOAD(buf[idx-w3+2]):0,
		x+3<iw&&y-3>=0?LOAD(buf[idx-w3+3]):0,

		x-3>=0&&y-2>=0?LOAD(buf[idx-w2-3]):0,
		x-2>=0&&y-2>=0?LOAD(buf[idx-w2-2]):0,
		x-1>=0&&y-2>=0?LOAD(buf[idx-w2-1]):0,
		        y-2>=0?LOAD(buf[idx-w2  ]):0,
		x+1<iw&&y-2>=0?LOAD(buf[idx-w2+1]):0,
		x+2<iw&&y-2>=0?LOAD(buf[idx-w2+2]):0,
		x+3<iw&&y-2>=0?LOAD(buf[idx-w2+3]):0,

		x-3>=0&&y-1>=0?LOAD(buf[idx-iw-3]):0,
		x-2>=0&&y-1>=0?LOAD(buf[idx-iw-2]):0,
		x-1>=0&&y-1>=0?LOAD(buf[idx-iw-1]):0,
		        y-1>=0?LOAD(buf[idx-iw  ]):0,
		x+1<iw&&y-1>=0?LOAD(buf[idx-iw+1]):0,
		x+2<iw&&y-1>=0?LOAD(buf[idx-iw+2]):0,
		x+3<iw&&y-1>=0?LOAD(buf[idx-iw+3]):0,

		x-3>=0        ?LOAD(buf[idx   -3]):0,
		x-2>=0        ?LOAD(buf[idx   -2]):0,
		x-1>=0        ?LOAD(buf[idx   -1]):0,
	};
	DataType *ptr=nb+_countof(nb);

	*ctx++ = ptr[-1];//L
	*ctx++ = ptr[CW*-1];//T
	*ctx++ = ptr[CW*-1-1];//TL
	*ctx++ = ptr[CW*-1+1];//TR
	*ctx++ = ptr[-1]*2-ptr[-2];//L*2-L2
	*ctx++ = ptr[CW*-1]*2-ptr[CW*-2];//T*2-T2
	*ctx++ = ptr[-1]*3-ptr[-2]*3+ptr[-3];//L*3-L2*3+L3
	*ctx++ = ptr[CW*-1]*3-ptr[CW*-2]*3+ptr[CW*-3];//T*3-T2*3+T3
	*ctx++ = clamp4(ptr[-1]+ptr[CW*-1]-ptr[CW*-1-1], ptr[-1], ptr[CW*-1], ptr[CW*-1-1], ptr[CW*-1+1]);//T+L-TL (grad)
	*ctx++ = clamp4(ptr[-1]+ptr[CW*-1+1]-ptr[CW*-1], ptr[-1], ptr[CW*-1], ptr[CW*-1-1], ptr[CW*-1+1]);//L+TR-T
	*ctx++ = clamp4(ptr[CW*-1]+ptr[CW*-1-1]-ptr[CW*-2-1], ptr[-1], ptr[CW*-1], ptr[CW*-1-1], ptr[CW*-1+1]);//T+TL-T2L
	*ctx++ = clamp4(ptr[CW*-1]+ptr[CW*-1+1]-ptr[CW*-2+1], ptr[-1], ptr[CW*-1+1], ptr[CW*-1], ptr[CW*-1+2]);//T+TR-T2R
	*ctx++ = (ptr[-1]+ptr[CW*-1+1])/2;
	*ctx++ = ptr[CW*-3];//T3
	*ctx++ = clamp4(ptr[CW*-1]*2-ptr[CW*-2], ptr[-1], ptr[CW*-1], ptr[CW*-1+1], ptr[CW*-1+2]);
	*ctx++ = (ptr[CW*-1]+ptr[CW*-3])/2;
	*ctx++ = ((ptr[-1]+ptr[CW*-1])*3-ptr[CW*-1-1]*2)/2;
	*ctx++ = (ptr[-1]+ptr[CW*-1])/2;//(L+T)/2
	*ctx++ = (-ptr[-2]+ptr[-1]*4+ptr[CW*-1]*4-ptr[CW*-2])/6;//(-L2+L*4+T*4-T2)
	*ctx++ = (ptr[-3]-ptr[-2]*6+ptr[-1]*15+ptr[CW*-1]*15-ptr[CW*-2]*6+ptr[CW*-3])/20;//(L3-L2*6+L*15+T*15-T2*6+T3)
#if 0
	*ctx++ = x-3>=0&&y-3>=0?LOAD(buf[idx-w3-3]):0;
	*ctx++ = x-2>=0&&y-3>=0?LOAD(buf[idx-w3-2]):0;
	*ctx++ = x-1>=0&&y-3>=0?LOAD(buf[idx-w3-1]):0;
	*ctx++ =         y-3>=0?LOAD(buf[idx-w3  ]):0;
	*ctx++ = x+1<iw&&y-3>=0?LOAD(buf[idx-w3+1]):0;
	*ctx++ = x+2<iw&&y-3>=0?LOAD(buf[idx-w3+2]):0;
	*ctx++ = x+3<iw&&y-3>=0?LOAD(buf[idx-w3+3]):0;

	*ctx++ = x-3>=0&&y-2>=0?LOAD(buf[idx-w2-3]):0;
	*ctx++ = x-2>=0&&y-2>=0?LOAD(buf[idx-w2-2]):0;
	*ctx++ = x-1>=0&&y-2>=0?LOAD(buf[idx-w2-1]):0;
	*ctx++ =         y-2>=0?LOAD(buf[idx-w2  ]):0;
	*ctx++ = x+1<iw&&y-2>=0?LOAD(buf[idx-w2+1]):0;
	*ctx++ = x+2<iw&&y-2>=0?LOAD(buf[idx-w2+2]):0;
	*ctx++ = x+3<iw&&y-2>=0?LOAD(buf[idx-w2+3]):0;

	*ctx++ = x-3>=0&&y-1>=0?LOAD(buf[idx-iw-3]):0;
	*ctx++ = x-2>=0&&y-1>=0?LOAD(buf[idx-iw-2]):0;
	*ctx++ = x-1>=0&&y-1>=0?LOAD(buf[idx-iw-1]):0;
	*ctx++ =         y-1>=0?LOAD(buf[idx-iw  ]):0;
	*ctx++ = x+1<iw&&y-1>=0?LOAD(buf[idx-iw+1]):0;
	*ctx++ = x+2<iw&&y-1>=0?LOAD(buf[idx-iw+2]):0;
	*ctx++ = x+3<iw&&y-1>=0?LOAD(buf[idx-iw+3]):0;

	*ctx++ = x-3>=0        ?LOAD(buf[idx   -3]):0;
	*ctx++ = x-2>=0        ?LOAD(buf[idx   -2]):0;
	*ctx++ = x-1>=0        ?LOAD(buf[idx   -1]):0;
#endif
}
#endif
static void get_nb2(const char *buf, const char *errors, int iw, int ih, int kx, int ky, DataType *ctx)
{
	const int CW=7;
	int idx=iw*ky+kx, w3=iw*3, w2=iw*2;
	DataType
		cT6  =         ky-6>=0?LOAD(buf[idx-iw*6  ]):0,
		cT5  =         ky-5>=0?LOAD(buf[idx-iw*5  ]):0,
		cT4L3=kx-3>=0&&ky-4>=0?LOAD(buf[idx-iw*4-3]):0,
		cT4  =         ky-4>=0?LOAD(buf[idx-iw*4  ]):0,
		cT4R3=kx+3<iw&&ky-4>=0?LOAD(buf[idx-iw*4+3]):0,
		cT3L5=kx-5>=0&&ky-3>=0?LOAD(buf[idx-iw*3-5]):0,
		cT3L4=kx-4>=0&&ky-3>=0?LOAD(buf[idx-iw*3-4]):0,
		cT3L2=kx-2>=0&&ky-3>=0?LOAD(buf[idx-iw*3-2]):0,
		cT3L =kx-1>=0&&ky-3>=0?LOAD(buf[idx-iw*3-1]):0,
		cT3  =         ky-3>=0?LOAD(buf[idx-iw*3  ]):0,
		cT3R =kx+1<iw&&ky-3>=0?LOAD(buf[idx-iw*3+1]):0,
		cT3R2=kx+2<iw&&ky-3>=0?LOAD(buf[idx-iw*3+2]):0,
		cT3R3=kx+3<iw&&ky-3>=0?LOAD(buf[idx-iw*3+3]):0,
		cT3R4=kx+4<iw&&ky-3>=0?LOAD(buf[idx-iw*3+4]):0,
		cT2L3=kx-3>=0&&ky-2>=0?LOAD(buf[idx-iw*2-3]):0,
		cT2L2=kx-2>=0&&ky-2>=0?LOAD(buf[idx-iw*2-2]):0,
		cT2L =kx-1>=0&&ky-2>=0?LOAD(buf[idx-iw*2-1]):0,
		cT2  =         ky-2>=0?LOAD(buf[idx-iw*2  ]):0,
		cT2R =kx+1<iw&&ky-2>=0?LOAD(buf[idx-iw*2+1]):0,
		cT2R2=kx+2<iw&&ky-2>=0?LOAD(buf[idx-iw*2+2]):0,
		cT2R3=kx+3<iw&&ky-2>=0?LOAD(buf[idx-iw*2+3]):0,
		cT2R4=kx+4<iw&&ky-2>=0?LOAD(buf[idx-iw*2+4]):0,
		cTL3 =kx-3>=0&&ky-1>=0?LOAD(buf[idx-iw  -3]):0,
		cTL2 =kx-2>=0&&ky-1>=0?LOAD(buf[idx-iw  -2]):0,
		cTL  =kx-1>=0&&ky-1>=0?LOAD(buf[idx-iw  -1]):0,
		cT   =kx  <iw&&ky-1>=0?LOAD(buf[idx-iw    ]):0,
		cTR  =kx+1<iw&&ky-1>=0?LOAD(buf[idx-iw  +1]):0,
		cTR2 =kx+2<iw&&ky-1>=0?LOAD(buf[idx-iw  +2]):0,
		cTR3 =kx+3<iw&&ky-1>=0?LOAD(buf[idx-iw  +3]):0,
		cTR4 =kx+4<iw&&ky-1>=0?LOAD(buf[idx-iw  +4]):0,
		cTR5 =kx+5<iw&&ky-1>=0?LOAD(buf[idx-iw  +5]):0,
		cTR6 =kx+6<iw&&ky-1>=0?LOAD(buf[idx-iw  +6]):0,
		cTR7 =kx+7<iw&&ky-1>=0?LOAD(buf[idx-iw  +7]):0,
		cL6  =kx-6>=0         ?LOAD(buf[idx     -6]):0,
		cL5  =kx-5>=0         ?LOAD(buf[idx     -5]):0,
		cL4  =kx-4>=0         ?LOAD(buf[idx     -4]):0,
		cL3  =kx-3>=0         ?LOAD(buf[idx     -3]):0,
		cL2  =kx-2>=0         ?LOAD(buf[idx     -2]):0,
		cL   =kx-1>=0         ?LOAD(buf[idx     -1]):0;

	*ctx++ = clamp4(cL+cT-cTL, cL, cTL, cT, cTR);//0
	*ctx++ = clip(cL+cT-cTL);//1
	*ctx++ = clamp4(cL+cTR-cT, cL, cTL, cT, cTR);//2
	*ctx++ = clip(cL+cTR-cT);//3
	*ctx++ = clamp4(cT+cTL-cT2L, cL, cTL, cT, cTR);//4
	*ctx++ = clip(cT+cTL-cT2L);//5
	*ctx++ = clamp4(cT+cTR-cT2R, cL, cT, cTR, cTR2);//6
	*ctx++ = clip(cT+cTR-cT2R);//7
	*ctx++ = (cL+cTR2)/2;//8
	*ctx++ = clip(cT*3-cT2*3+cT3);//9
	*ctx++ = clip(cL*3-cL2*3+cL3);//10
	*ctx++ = (cL+clip(cTR*3-cT2R*3+cT3R))/2;//11
	*ctx++ = (cL+clip(cTR2*3-cT2R3*3+cT3R4))/2;//12
	*ctx++ = clip(cT2+cT4-cT6);//13
	*ctx++ = clip(cL2+cL4-cL6);//14
	*ctx++ = clip((cT5-6*cT4+15*cT3-20*cT2+15*cT+clamp4(cL*2-cTL2, cL, cTL, cT, cT2))/6);//15
	*ctx++ = clip((-3*cL2+8*cL+clamp4(3*cTR2-3*cT2R2+cT3R2, cTR, cTR2, cTR3, cTR4))/6);//16
	*ctx++ = clip(cT2+cTL-cT3L);//17
	*ctx++ = clip(cT2+cTR-cT3R);//18
	*ctx++ = clip((cL*2+cTL) - (cL2*2+cTL2) + cL3);//19
	*ctx++ = clip(3*(cTL+cTL2)/2-cT2L3*3+(cT3L4+cT3L5)/2);//20
	*ctx++ = clip(cTR2+cTR-cT2R3);//21
	*ctx++ = clip(cTL2+cL2-cL4);//22
	*ctx++ = clip(((cL+cTL)*3-cTL2*6+cTL3+cT2L3)/2);//23
	*ctx++ = clip((cTR*2+cTR2) - (cT2R2+cT3R2*2) + cT4R3);//24
	*ctx++ = cT6;//25
	*ctx++ = (cTR4+cTR6)/2;//26
	*ctx++ = (cL4+cL6)/2;//27
	*ctx++ = (cL+cT+cTR5+cTR7)/4;//28
	*ctx++ = clip(cTR3+cL-cTR2);//29
	*ctx++ = clip(4*cT3-3*cT4);//30
	*ctx++ = clip(cT+cT2-cT3);//31
	*ctx++ = clip(cL+cL2-cL3);//32
	*ctx++ = clip(cL+cTR2-cTR);//33
	*ctx++ = clip(cL2+cTR2-cT);//34
	*ctx++ = (clip(cL*2-cTL)+clip(cL*2-cTL2)+cT+cTR)/4;//35
	*ctx++ = clamp4(cT*2-cT2, cL, cT, cTR, cTR2);//36
	*ctx++ = (cT+cT3)/2;//37
	*ctx++ = clip(cT2+cL-cT2L);//38
	*ctx++ = clip(cTR2+cT-cT2R2);//39
	*ctx++ = clip((4*cL3-15*cL2+20*cL+clip(cTR2*2-cT2R2))/10);//40
	*ctx++ = clip((cT3R3-4*cT2R2+6*cTR+clip(cL*3-cTL*3+cT2L))/4);//41
	*ctx++ = clip((cT*2+cTR) - (cT2+2*cT2R) + cT3R);//42
	*ctx++ = clip((cTL*2+cT2L) - (cT2L2+cT3L2*2) + cT4L3);//43
	*ctx++ = clip(cT2L2+cL-cT2L3);//44
	*ctx++ = clip((-cT4+5*cT3-10*cT2+10*cT+clip(cL*4-cTL2*6+cT2L3*4-cT3L4))/5);//45
	*ctx++ = clip(cTR2+clip(cTR3*2-cT2R4-cTR4));//46
	*ctx++ = clip(cTL+cL-cTL2);//47
	*ctx++ = clip((cT*2+cTL) - (cT2+2*cT2L) + cT3L);//48
	*ctx++ = clip(cT2+clip(cTR2*2-cT2R3) - cT2R);//49
	*ctx++ = clip((-cL4+5*cL3-10*cL2+10*cL+clip(cTR*2-cT2R))/5);//50
	*ctx++ = clip((-cL5+4*cL4-5*cL3+5*cL+clip(cTR*2-cT2R))/4);//51
	*ctx++ = clip((cL3-4*cL2+6*cL+clip(cTR*3-cT2R*3+cT3R))/4);//52
	*ctx++ = clip((-cT2R2+3*cTR+clip(4*cL-6*cTL+4*cT2L-cT3L))/3);//53
	*ctx++ = ((cL+cT)*3-cTL*2)/4;//54
	*ctx++=         ky-6>=0?LOAD(errors[idx-iw*6  ]):0,
	*ctx++=         ky-5>=0?LOAD(errors[idx-iw*5  ]):0,
	*ctx++=kx-3>=0&&ky-4>=0?LOAD(errors[idx-iw*4-3]):0,
	*ctx++=         ky-4>=0?LOAD(errors[idx-iw*4  ]):0,
	*ctx++=kx+3<iw&&ky-4>=0?LOAD(errors[idx-iw*4+3]):0,
	*ctx++=kx-5>=0&&ky-3>=0?LOAD(errors[idx-iw*3-5]):0,
	*ctx++=kx-4>=0&&ky-3>=0?LOAD(errors[idx-iw*3-4]):0,
	*ctx++=kx-2>=0&&ky-3>=0?LOAD(errors[idx-iw*3-2]):0,
	*ctx++=kx-1>=0&&ky-3>=0?LOAD(errors[idx-iw*3-1]):0,
	*ctx++=         ky-3>=0?LOAD(errors[idx-iw*3  ]):0,
	*ctx++=kx+1<iw&&ky-3>=0?LOAD(errors[idx-iw*3+1]):0,
	*ctx++=kx+2<iw&&ky-3>=0?LOAD(errors[idx-iw*3+2]):0,
	*ctx++=kx+3<iw&&ky-3>=0?LOAD(errors[idx-iw*3+3]):0,
	*ctx++=kx+4<iw&&ky-3>=0?LOAD(errors[idx-iw*3+4]):0,
	*ctx++=kx-3>=0&&ky-2>=0?LOAD(errors[idx-iw*2-3]):0,
	*ctx++=kx-2>=0&&ky-2>=0?LOAD(errors[idx-iw*2-2]):0,
	*ctx++=kx-1>=0&&ky-2>=0?LOAD(errors[idx-iw*2-1]):0,
	*ctx++=         ky-2>=0?LOAD(errors[idx-iw*2  ]):0,
	*ctx++=kx+1<iw&&ky-2>=0?LOAD(errors[idx-iw*2+1]):0,
	*ctx++=kx+2<iw&&ky-2>=0?LOAD(errors[idx-iw*2+2]):0,
	*ctx++=kx+3<iw&&ky-2>=0?LOAD(errors[idx-iw*2+3]):0,
	*ctx++=kx+4<iw&&ky-2>=0?LOAD(errors[idx-iw*2+4]):0,
	*ctx++=kx-3>=0&&ky-1>=0?LOAD(errors[idx-iw  -3]):0,
	*ctx++=kx-2>=0&&ky-1>=0?LOAD(errors[idx-iw  -2]):0,
	*ctx++=kx-1>=0&&ky-1>=0?LOAD(errors[idx-iw  -1]):0,
	*ctx++=kx  <iw&&ky-1>=0?LOAD(errors[idx-iw    ]):0,
	*ctx++=kx+1<iw&&ky-1>=0?LOAD(errors[idx-iw  +1]):0,
	*ctx++=kx+2<iw&&ky-1>=0?LOAD(errors[idx-iw  +2]):0,
	*ctx++=kx+3<iw&&ky-1>=0?LOAD(errors[idx-iw  +3]):0,
	*ctx++=kx+4<iw&&ky-1>=0?LOAD(errors[idx-iw  +4]):0,
	*ctx++=kx+5<iw&&ky-1>=0?LOAD(errors[idx-iw  +5]):0,
	*ctx++=kx+6<iw&&ky-1>=0?LOAD(errors[idx-iw  +6]):0,
	*ctx++=kx+7<iw&&ky-1>=0?LOAD(errors[idx-iw  +7]):0,
	*ctx++=kx-6>=0         ?LOAD(errors[idx     -6]):0,
	*ctx++=kx-5>=0         ?LOAD(errors[idx     -5]):0,
	*ctx++=kx-4>=0         ?LOAD(errors[idx     -4]):0,
	*ctx++=kx-3>=0         ?LOAD(errors[idx     -3]):0,
	*ctx++=kx-2>=0         ?LOAD(errors[idx     -2]):0,
	*ctx++=kx-1>=0         ?LOAD(errors[idx     -1]):0;
}
static void eval_fwd(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params *p)
{
	get_nb2(src, errors, iw, ih, kx, ky, t->nb);

	//get_neighbors(src, iw, ih, kx, ky, t->nb);
	//get_neighbors(errors, iw, ih, kx, ky, t->nb+NFEATURES);

	linear(t->net1, p->weight1, t->nb, p->bias1, NF0, NF1);	//n1 = w1*nb + b1
	act(t->x1, t->net1, NF1);								//x1 = act(n1)

	linear(t->net2, p->weight2, t->x1, p->bias2, NF1, NF2);	//n2 = w2*x1 + b2
	act(t->x2, t->net2, NF2);								//x2 = act(n2)

	linear(t->net3, p->weight3, t->x2, p->bias3, NF2, NF3);	//n3 = w3*x2
	act(t->x3, t->net3, NF3);								//x3 = act(n3)

	linear(&t->pred, p->weight4, t->x3, 0, NF3, 1);			//pred = w4*x3
}
static void train(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params *p, BwdTemps *b, Params *g, float lr)
{
	//fwd
	eval_fwd(src, errors, iw, ih, kx, ky, t, p);

	//loss
	DataType curr=LOAD(src[iw*ky+kx]);
	t->expr=t->pred-curr;
	t->loss=ABS(t->expr);


	//bwd
	b->dL_dp=SGN(t->expr);//dL_dp = sgn(p-x)

	vec_scalar(g->weight4, t->x3, b->dL_dp, NF3);//dL_dw3 = dL_dp * dp_dw3		= sgn(p-x) * x2T

	//dL_dn3 = dL_dp * dp_dx3 * dx3_dn3			= sgn(p-x) * w4T * act_dash(n3)
	act_dash(b->dL_dn3, t->net3, NF3);
	vec_ew(b->dL_dn3, b->dL_dn3, p->weight4, NF3);
	vec_scalar(b->dL_dn3, b->dL_dn3, b->dL_dp, NF3);


	memcpy(g->bias3, b->dL_dn3, NF3*sizeof(DataType));//dL_db3 = dL_dn3
	vec_outer(g->weight3, b->dL_dn3, t->x2, NF3, NF2);//dL_dw3 = dL_dn3 * dn3_dw3		= dL_dn3 * x2T (outer product)

	//dL_dn2 = matmul(dL_dn3T, w3)T .* act_dash(n2)
	matmul(b->dL_dn2, b->dL_dn3, p->weight3, 1, NF3, NF2);
	act_dash(b->actdash_n2, t->net2, NF2);
	vec_ew(b->dL_dn2, b->dL_dn2, b->actdash_n2, NF2);


	memcpy(g->bias2, b->dL_dn2, NF2*sizeof(DataType));//dL_db2 = dL_dn2
	vec_outer(g->weight2, b->dL_dn2, t->x1, NF2, NF1);//dL_dw2 = dL_dn2 * dn2_dw2		= dL_dn2 * x1T (outer product)

	//dL_dn1 = matmul(dL_dn2T, w2)T .* act_dash(n1)
	matmul(b->dL_dn1, b->dL_dn2, p->weight2, 1, NF2, NF1);
	act_dash(b->actdash_n1, t->net1, NF1);
	vec_ew(b->dL_dn1, b->dL_dn1, b->actdash_n1, NF1);
				

	memcpy(g->bias1, b->dL_dn1, NF1*sizeof(DataType));//dL_db1 = dL_dn1
	vec_outer(g->weight1, b->dL_dn1, t->nb, NF1, NF0);//dL_dw1 = dL_dn1 * dn1_dw1


	//update
	DataType *params=(DataType*)p, *gradient=(DataType*)g;
	for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
		params[k]-=LEARNING_RATE(gradient[k], lr);
		//params[k]-=LEARNING_RATE(gradient[k], 0.001f);
}
#if 0
int pred_learned_eval(char *buf, int iw, int ih, int x, int y, int kc, short *params)//32 params (24 actual + 8 dummy)
{
	int idx=iw*y+x, w3=iw*3, w2=iw*2;
	__declspec(align(32)) short nb[32]=
	{
		x-3>=0&&y-3>=0?buf[(idx-w3-3)<<2|kc]<<8:0,
		x-2>=0&&y-3>=0?buf[(idx-w3-2)<<2|kc]<<8:0,
		x-1>=0&&y-3>=0?buf[(idx-w3-1)<<2|kc]<<8:0,
		        y-3>=0?buf[(idx-w3  )<<2|kc]<<8:0,
		x+1<iw&&y-3>=0?buf[(idx-w3+1)<<2|kc]<<8:0,
		x+2<iw&&y-3>=0?buf[(idx-w3+2)<<2|kc]<<8:0,
		x+3<iw&&y-3>=0?buf[(idx-w3+3)<<2|kc]<<8:0,
		
		x-3>=0&&y-2>=0?buf[(idx-w2-3)<<2|kc]<<8:0,
		x-2>=0&&y-2>=0?buf[(idx-w2-2)<<2|kc]<<8:0,
		x-1>=0&&y-2>=0?buf[(idx-w2-1)<<2|kc]<<8:0,
		        y-2>=0?buf[(idx-w2  )<<2|kc]<<8:0,
		x+1<iw&&y-2>=0?buf[(idx-w2+1)<<2|kc]<<8:0,
		x+2<iw&&y-2>=0?buf[(idx-w2+2)<<2|kc]<<8:0,
		x+3<iw&&y-2>=0?buf[(idx-w2+3)<<2|kc]<<8:0,
		
		x-3>=0&&y-1>=0?buf[(idx-iw-3)<<2|kc]<<8:0,
		x-2>=0&&y-1>=0?buf[(idx-iw-2)<<2|kc]<<8:0,
		x-1>=0&&y-1>=0?buf[(idx-iw-1)<<2|kc]<<8:0,
		        y-1>=0?buf[(idx-iw  )<<2|kc]<<8:0,
		x+1<iw&&y-1>=0?buf[(idx-iw+1)<<2|kc]<<8:0,
		x+2<iw&&y-1>=0?buf[(idx-iw+2)<<2|kc]<<8:0,
		x+3<iw&&y-1>=0?buf[(idx-iw+3)<<2|kc]<<8:0,
		
		x-3>=0        ?buf[(idx   -3)<<2|kc]<<8:0,
		x-2>=0        ?buf[(idx   -2)<<2|kc]<<8:0,
		x-1>=0        ?buf[(idx   -1)<<2|kc]<<8:0,
	};
	__m256i vn1=_mm256_load_si256((__m256i*)nb);
	__m256i vn2=_mm256_load_si256((__m256i*)nb+1);
	__m256i vp1=_mm256_load_si256((__m256i*)params);
	__m256i vp2=_mm256_load_si256((__m256i*)params+1);
	_mm256_mulhi_epu16(vn1, vp1);
}
#endif
//void pred_learned_train(char *buf, int iw, int ih, int kx, int ky, short *params)
//{
//}
static void initialize_all(Params *p)
{
	XOROSHIRO128_RESET();
	initialize(p->weight1, _countof(p->weight1), (DataType)sqrtf(NF0));
	initialize(p->bias1, _countof(p->bias1), (DataType)sqrtf(NF0));
	initialize(p->weight2, _countof(p->weight2), (DataType)sqrtf(NF1));
	initialize(p->bias2, _countof(p->bias2), (DataType)sqrtf(NF1));
	initialize(p->weight3, _countof(p->weight3), (DataType)sqrtf(NF2));
	initialize(p->bias3, _countof(p->bias3), (DataType)sqrtf(NF2));
	initialize(p->weight4, _countof(p->weight4), (DataType)sqrtf(NF3));
}
void pred_learned(char *buf, int iw, int ih, int fwd)
{
	double t_start=time_ms();
	int res=iw*ih;
	char *src=(char*)malloc(res), *dst=(char*)malloc(res);
	Temps *t=(Temps*)malloc(sizeof(Temps));
	Params *p=(Params*)malloc(sizeof(Params));
	BwdTemps *b=(BwdTemps*)malloc(sizeof(BwdTemps));
	Params *g=(Params*)malloc(sizeof(Params));
	if(!src||!dst||!t||!p||!b||!g)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	ArrayHandle title;
	STR_ALLOC(title, 1024);
	get_window_title(title->data, 1024);
	//memset(t, 0, sizeof(Temps));
			
	//initialization
	//const int nparams=sizeof(Params)/sizeof(int);
	//DataType *ptr=(DataType*)p;
	//for(int k=0;k<nparams;++k)
	//	ptr[k]=sample;
	//	//ptr[k]=rand()&((ONE-1)>>FRACBITS/2);

	int pred, idx;
	char *buf2=fwd?src:dst, *errors=fwd?dst:src;
	for(int kc=0;kc<3;++kc)
	{
		initialize_all(p);

		for(int k=0;k<res;++k)
			src[k]=buf[k<<2|kc];

		for(int ky=0;ky<ih;++ky)
		{
			float lr=0.001f;
			if(!(ky&15))
			{
				TimeInfo ti;
				parsetimedelta(time_ms()-t_start, &ti);
				set_window_title("%d/3, %d/%d, %.2lf%% - %02d-%02d-%06.3f", kc+1, ky+1, ih, 100.*(ih*kc+ky+1)/(ih*3), ti.hours, ti.mins, ti.secs);
			}

			for(int kx=0;kx<iw;++kx)
			{
				//for(int k=0;k<1;++k)
				//{
				//	//train(buf2, errors, iw, ih, kx, ky-3, t, p, b, g);
				//	//train(buf2, errors, iw, ih, kx-3, ky, t, p, b, g);
				//	//train(buf2, errors, iw, ih, kx, ky-2, t, p, b, g);
				//	//train(buf2, errors, iw, ih, kx-2, ky, t, p, b, g);
				//	//train(buf2, errors, iw, ih, kx, ky-1, t, p, b, g);
				//	train(buf2, errors, iw, ih, kx-1, ky, t, p, b, g, lr);
				//	//lr*=0.9f;
				//}
#if 1
				int niter=kc==1?1:6;
				for(int k=0;k<niter;++k)
				{
					train(buf2, errors, iw, ih, kx-5, ky-5, t, p, b, g, lr/sqrtf(5*5+5*5));
					train(buf2, errors, iw, ih, kx-4, ky-5, t, p, b, g, lr/sqrtf(4*4+5*5));
					train(buf2, errors, iw, ih, kx-3, ky-5, t, p, b, g, lr/sqrtf(3*3+5*5));
					train(buf2, errors, iw, ih, kx-2, ky-5, t, p, b, g, lr/sqrtf(2*2+5*5));
					train(buf2, errors, iw, ih, kx-1, ky-5, t, p, b, g, lr/sqrtf(1*1+5*5));
					train(buf2, errors, iw, ih, kx  , ky-5, t, p, b, g, lr/sqrtf(0*0+5*5));
					train(buf2, errors, iw, ih, kx+1, ky-5, t, p, b, g, lr/sqrtf(1*1+5*5));
					train(buf2, errors, iw, ih, kx+2, ky-5, t, p, b, g, lr/sqrtf(2*2+5*5));
					train(buf2, errors, iw, ih, kx+3, ky-5, t, p, b, g, lr/sqrtf(3*3+5*5));
					train(buf2, errors, iw, ih, kx+4, ky-5, t, p, b, g, lr/sqrtf(4*4+5*5));
					train(buf2, errors, iw, ih, kx+5, ky-5, t, p, b, g, lr/sqrtf(5*5+5*5));
				
					train(buf2, errors, iw, ih, kx-5, ky-4, t, p, b, g, lr/sqrtf(5*5+4*4));
					train(buf2, errors, iw, ih, kx-4, ky-4, t, p, b, g, lr/sqrtf(4*4+4*4));
					train(buf2, errors, iw, ih, kx-3, ky-4, t, p, b, g, lr/sqrtf(3*3+4*4));
					train(buf2, errors, iw, ih, kx-2, ky-4, t, p, b, g, lr/sqrtf(2*2+4*4));
					train(buf2, errors, iw, ih, kx-1, ky-4, t, p, b, g, lr/sqrtf(1*1+4*4));
					train(buf2, errors, iw, ih, kx  , ky-4, t, p, b, g, lr/sqrtf(0*0+4*4));
					train(buf2, errors, iw, ih, kx+1, ky-4, t, p, b, g, lr/sqrtf(1*1+4*4));
					train(buf2, errors, iw, ih, kx+2, ky-4, t, p, b, g, lr/sqrtf(2*2+4*4));
					train(buf2, errors, iw, ih, kx+3, ky-4, t, p, b, g, lr/sqrtf(3*3+4*4));
					train(buf2, errors, iw, ih, kx+4, ky-4, t, p, b, g, lr/sqrtf(4*4+4*4));
					train(buf2, errors, iw, ih, kx+5, ky-4, t, p, b, g, lr/sqrtf(5*5+4*4));
				
					train(buf2, errors, iw, ih, kx-5, ky-3, t, p, b, g, lr/sqrtf(5*5+3*3));
					train(buf2, errors, iw, ih, kx-4, ky-3, t, p, b, g, lr/sqrtf(4*4+3*3));
					train(buf2, errors, iw, ih, kx-3, ky-3, t, p, b, g, lr/sqrtf(3*3+3*3));
					train(buf2, errors, iw, ih, kx-2, ky-3, t, p, b, g, lr/sqrtf(2*2+3*3));
					train(buf2, errors, iw, ih, kx-1, ky-3, t, p, b, g, lr/sqrtf(1*1+3*3));
					train(buf2, errors, iw, ih, kx  , ky-3, t, p, b, g, lr/sqrtf(0*0+3*3));
					train(buf2, errors, iw, ih, kx+1, ky-3, t, p, b, g, lr/sqrtf(1*1+3*3));
					train(buf2, errors, iw, ih, kx+2, ky-3, t, p, b, g, lr/sqrtf(2*2+3*3));
					train(buf2, errors, iw, ih, kx+3, ky-3, t, p, b, g, lr/sqrtf(3*3+3*3));
					train(buf2, errors, iw, ih, kx+4, ky-3, t, p, b, g, lr/sqrtf(4*4+3*3));
					train(buf2, errors, iw, ih, kx+5, ky-3, t, p, b, g, lr/sqrtf(5*5+3*3));
				
					train(buf2, errors, iw, ih, kx-5, ky-2, t, p, b, g, lr/sqrtf(5*5+2*2));
					train(buf2, errors, iw, ih, kx-4, ky-2, t, p, b, g, lr/sqrtf(4*4+2*2));
					train(buf2, errors, iw, ih, kx-3, ky-2, t, p, b, g, lr/sqrtf(3*3+2*2));
					train(buf2, errors, iw, ih, kx-2, ky-2, t, p, b, g, lr/sqrtf(2*2+2*2));
					train(buf2, errors, iw, ih, kx-1, ky-2, t, p, b, g, lr/sqrtf(1*1+2*2));
					train(buf2, errors, iw, ih, kx  , ky-2, t, p, b, g, lr/sqrtf(0*0+2*2));
					train(buf2, errors, iw, ih, kx+1, ky-2, t, p, b, g, lr/sqrtf(1*1+2*2));
					train(buf2, errors, iw, ih, kx+2, ky-2, t, p, b, g, lr/sqrtf(2*2+2*2));
					train(buf2, errors, iw, ih, kx+3, ky-2, t, p, b, g, lr/sqrtf(3*3+2*2));
					train(buf2, errors, iw, ih, kx+4, ky-2, t, p, b, g, lr/sqrtf(4*4+2*2));
					train(buf2, errors, iw, ih, kx+5, ky-2, t, p, b, g, lr/sqrtf(5*5+2*2));
				
					train(buf2, errors, iw, ih, kx-5, ky-1, t, p, b, g, lr/sqrtf(5*5+1*1));
					train(buf2, errors, iw, ih, kx-4, ky-1, t, p, b, g, lr/sqrtf(4*4+1*1));
					train(buf2, errors, iw, ih, kx-3, ky-1, t, p, b, g, lr/sqrtf(3*3+1*1));
					train(buf2, errors, iw, ih, kx-2, ky-1, t, p, b, g, lr/sqrtf(2*2+1*1));
					train(buf2, errors, iw, ih, kx-1, ky-1, t, p, b, g, lr/sqrtf(1*1+1*1));
					train(buf2, errors, iw, ih, kx  , ky-1, t, p, b, g, lr/sqrtf(0*0+1*1));
					train(buf2, errors, iw, ih, kx+1, ky-1, t, p, b, g, lr/sqrtf(1*1+1*1));
					train(buf2, errors, iw, ih, kx+2, ky-1, t, p, b, g, lr/sqrtf(2*2+1*1));
					train(buf2, errors, iw, ih, kx+3, ky-1, t, p, b, g, lr/sqrtf(3*3+1*1));
					train(buf2, errors, iw, ih, kx+4, ky-1, t, p, b, g, lr/sqrtf(4*4+1*1));
					train(buf2, errors, iw, ih, kx+5, ky-1, t, p, b, g, lr/sqrtf(5*5+1*1));
				
					train(buf2, errors, iw, ih, kx-5, ky  , t, p, b, g, lr/5);
					train(buf2, errors, iw, ih, kx-4, ky  , t, p, b, g, lr/4);
					train(buf2, errors, iw, ih, kx-3, ky  , t, p, b, g, lr/3);
					train(buf2, errors, iw, ih, kx-2, ky  , t, p, b, g, lr/2);
					train(buf2, errors, iw, ih, kx-1, ky  , t, p, b, g, lr/1);
				}
#endif

				//for(int k=0;k<1;++k)
				//{
				//	train(buf2, errors, iw, ih, kx-1, ky, t, p, b, g);//CHEAT
				//}

				eval_fwd(buf2, errors, iw, ih, kx, ky, t, p);

				pred=STORE(t->pred);
				pred=CLAMP(-128, pred, 127);

				//if(pred)
				//	printf("");

				idx=iw*ky+kx;
				if(fwd)
					dst[idx]=src[idx]-pred;
				else
					dst[idx]=src[idx]+pred;

				//free(malloc(1));
				//train(buf2, iw, ih, kx, ky, t, p, b, g);
			}
		}

		for(int k=0;k<res;++k)
			buf[k<<2|kc]=dst[k];
	}
	free(src);
	free(dst);
	free(t);
	free(p);
	free(b);
	free(g);
#if 1
	{
		TimeInfo ti;
		parsetimedelta(time_ms()-t_start, &ti);
		set_window_title("3/3, %d/%d, 100%% - %02d-%02d-%06.3f", ih, ih, ti.hours, ti.mins, ti.secs);
	}
#else
	set_window_title("%s", title->data);
#endif
	array_free(&title);
}


#ifdef ALLOW_OPENCL
const char*		clerr2str(int error)
{
	const char *a=0;
#define 		EC(x)		case x:a=(const char*)#x;break;
#define 		EC2(n, x)	case n:a=(const char*)#x;break;
	switch(error)
	{
	EC(CL_SUCCESS)
	EC(CL_DEVICE_NOT_FOUND)
	EC(CL_DEVICE_NOT_AVAILABLE)
	EC(CL_COMPILER_NOT_AVAILABLE)
	EC(CL_MEM_OBJECT_ALLOCATION_FAILURE)
	EC(CL_OUT_OF_RESOURCES)
	EC(CL_OUT_OF_HOST_MEMORY)
	EC(CL_PROFILING_INFO_NOT_AVAILABLE)
	EC(CL_MEM_COPY_OVERLAP)
	EC(CL_IMAGE_FORMAT_MISMATCH)
	EC(CL_IMAGE_FORMAT_NOT_SUPPORTED)
	EC(CL_BUILD_PROGRAM_FAILURE)
	EC(CL_MAP_FAILURE)
//#ifdef CL_VERSION_1_1
	EC(CL_MISALIGNED_SUB_BUFFER_OFFSET)
	EC(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_COMPILE_PROGRAM_FAILURE)
	EC(CL_LINKER_NOT_AVAILABLE)
	EC(CL_LINK_PROGRAM_FAILURE)
	EC(CL_DEVICE_PARTITION_FAILED)
	EC(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
//#endif
	EC(CL_INVALID_VALUE)
	EC(CL_INVALID_DEVICE_TYPE)
	EC(CL_INVALID_PLATFORM)
	EC(CL_INVALID_DEVICE)
	EC(CL_INVALID_CONTEXT)
	EC(CL_INVALID_QUEUE_PROPERTIES)
	EC(CL_INVALID_COMMAND_QUEUE)
	EC(CL_INVALID_HOST_PTR)
	EC(CL_INVALID_MEM_OBJECT)
	EC(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
	EC(CL_INVALID_IMAGE_SIZE)
	EC(CL_INVALID_SAMPLER)
	EC(CL_INVALID_BINARY)
	EC(CL_INVALID_BUILD_OPTIONS)
	EC(CL_INVALID_PROGRAM)
	EC(CL_INVALID_PROGRAM_EXECUTABLE)
	EC(CL_INVALID_KERNEL_NAME)
	EC(CL_INVALID_KERNEL_DEFINITION)
	EC(CL_INVALID_KERNEL)
	EC(CL_INVALID_ARG_INDEX)
	EC(CL_INVALID_ARG_VALUE)
	EC(CL_INVALID_ARG_SIZE)
	EC(CL_INVALID_KERNEL_ARGS)
	EC(CL_INVALID_WORK_DIMENSION)
	EC(CL_INVALID_WORK_GROUP_SIZE)
	EC(CL_INVALID_WORK_ITEM_SIZE)
	EC(CL_INVALID_GLOBAL_OFFSET)
	EC(CL_INVALID_EVENT_WAIT_LIST)
	EC(CL_INVALID_EVENT)
	EC(CL_INVALID_OPERATION)
	EC(CL_INVALID_GL_OBJECT)
	EC(CL_INVALID_BUFFER_SIZE)
	EC(CL_INVALID_MIP_LEVEL)
	EC(CL_INVALID_GLOBAL_WORK_SIZE)
//#ifdef CL_VERSION_1_1
	EC(CL_INVALID_PROPERTY)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_INVALID_IMAGE_DESCRIPTOR)
	EC(CL_INVALID_COMPILER_OPTIONS)
	EC(CL_INVALID_LINKER_OPTIONS)
	EC(CL_INVALID_DEVICE_PARTITION_COUNT)
//#endif
//#ifdef CL_VERSION_2_0
	EC2(-69, CL_INVALID_PIPE_SIZE)
	EC2(-70, CL_INVALID_DEVICE_QUEUE)
//#endif
//#ifdef CL_VERSION_2_2
	EC2(-71, CL_INVALID_SPEC_ID)
	EC2(-72, CL_MAX_SIZE_RESTRICTION_EXCEEDED)
//#endif
//	EC(CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR)
//	EC(CL_PLATFORM_NOT_FOUND_KHR)
	EC2(-1002, CL_INVALID_D3D10_DEVICE_KHR)
	EC2(-1003, CL_INVALID_D3D10_RESOURCE_KHR)
	EC2(-1004, CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1005, CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1006, CL_INVALID_D3D11_DEVICE_KHR)
	EC2(-1007, CL_INVALID_D3D11_RESOURCE_KHR)
	EC2(-1008, CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1009, CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR)
#ifndef __linux__
	EC2(-1010, CL_INVALID_D3D9_DEVICE_NV_or_CL_INVALID_DX9_DEVICE_INTEL)
	EC2(-1011, CL_INVALID_D3D9_RESOURCE_NV_or_CL_INVALID_DX9_RESOURCE_INTEL)
	EC2(-1012, CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV_or_CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL)
	EC2(-1013, CL_D3D9_RESOURCE_NOT_ACQUIRED_NV_or_CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL)
#endif
	EC2(-1092, CL_EGL_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1093, CL_INVALID_EGL_OBJECT_KHR)
	EC2(-1094, CL_INVALID_ACCELERATOR_INTEL)
	EC2(-1095, CL_INVALID_ACCELERATOR_TYPE_INTEL)
	EC2(-1096, CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL)
	EC2(-1097, CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL)
	EC2(-1098, CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL)
	EC2(-1099, CL_INVALID_VA_API_MEDIA_SURFACE_INTEL)
	EC2(-1101, CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL)
	case 1:a="File failure";break;//
	default:
		a="???";
		break;
	}
#undef			EC
#undef			EC2
	return a;
}
#define CHECKCL(E) (!(E)||LOG_ERROR("CL Error %s", clerr2str(E)))

static int clctxinitialized=0;
static cl_platform_id platform=0;
static cl_device_id device=0;
static cl_context context=0;
static cl_command_queue commandqueue=0;
static cl_program program=0;

#define CLKERNELNALELIST CLKERNEL(train)
typedef enum		CLKernelIdxEnum
{
#define				CLKERNEL(LABEL)	OCL_##LABEL,
CLKERNELNALELIST
#undef				CLKERNEL
	OCL_NKERNELS,
} CLKernelIdx;
const char			*kernelnames[]=
{
#define				CLKERNEL(LABEL)	#LABEL,
CLKERNELNALELIST
#undef				CLKERNEL
};
static cl_kernel kernels[OCL_NKERNELS]={0};
int init_ocl()
{
	int error=0;
	unsigned count=0;
	if(clctxinitialized)
		return 1;
	clctxinitialized=1;

	//get platform
	error=clGetPlatformIDs(0, 0, &count);	CHECKCL(error);
	if(!count)
	{
		LOG_ERROR("No OpenCL platforms");
		return 0;
	}
	error=clGetPlatformIDs(1, &platform, 0);	CHECKCL(error);

	//get device
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, 0, &count);	CHECKCL(error);
	if(!count)
	{
		LOG_ERROR("No OpenCL devices");
		return 0;
	}
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, 0);	CHECKCL(error);
	
	//create context
#ifdef CL_GL_INTEROP
	cl_context_properties properties[8]={};
	if(cl_gl_interop)
	{
		auto gl_context=eglGetCurrentContext();//changes when resuming
		auto egl_display=eglGetCurrentDisplay();
		properties[0]=CL_GL_CONTEXT_KHR,	properties[1]=(cl_context_properties)gl_context;//https://stackoverflow.com/questions/26802905/getting-opengl-buffers-using-opencl
		properties[2]=CL_EGL_DISPLAY_KHR,	properties[3]=(cl_context_properties)egl_display;
		properties[4]=CL_CONTEXT_PLATFORM,	properties[5]=(cl_context_properties)platform;
		properties[6]=0, properties[7]=0;
	}
	else
	{
		properties[0]=CL_CONTEXT_PLATFORM, properties[1]=(cl_context_properties)platform;
		properties[2]=0, properties[3]=0;
	}
	context=clCreateContext(properties, 1, &devices[0], nullptr, nullptr, &error);	CHECKCL(error);
#else
	context=clCreateContext(0, 1, &device, 0, 0, &error);	CHECKCL(error);
#endif
	
	//create command queue
#if CL_TARGET_OPENCL_VERSION>=200
	//cl_queue_properties properties[]=
	//{
	//	CL_QUEUE_PROPERTIES, 0,
	//	CL_QUEUE_SIZE, 0,
	//	0,
	//};
	commandqueue=clCreateCommandQueueWithProperties(context, device, 0, &error);	CHECKCL(error);
#else
	commandqueue=clCreateCommandQueue(context, device, 0, &error);	CHECKCL(error);
#endif

	//build kernels
	{
		ArrayHandle srctext=load_file("E:/C/pxView3D/pxView3D/cl_kernels.h", 0, 0);
		const char *k_src=(const char*)srctext->data;
		size_t k_len=srctext->count;
		
		program=clCreateProgramWithSource(context, 1, (const char**)&k_src, &k_len, &error);	CHECKCL(error);
		error=clBuildProgram(program, 1, &device, g_buf, 0, 0);
		if(error)
		{
			size_t retlen=0;
			error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, G_BUF_SIZE, g_buf, &retlen);
			if(retlen>G_BUF_SIZE)
			{
				char *buf=(char*)malloc(retlen+10);
				error=clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, retlen+10, buf, &retlen);	CHECKCL(error);
				messagebox(MBOX_OK, "OpenCL compilation failed", "%s", buf);
				//printf("\nOpenCL compilation failed:\n%s\n", buf);
				free(buf);
				LOG_ERROR("Aborting");
			}
			else
				LOG_ERROR("OpenCL Compilation failed:\n%s\n", g_buf);
		}
		array_free(&srctext);
	}

	//fetch entry points
	for(int k=0;k<OCL_NKERNELS;++k)
	{
		kernels[k]=clCreateKernel(program, kernelnames[k], &error);
		if(error)
		{
			LOG_ERROR("Couldn't find kernel %s", kernelnames[k]);
			return 0;
		}
	}
	return 1;
}
#endif