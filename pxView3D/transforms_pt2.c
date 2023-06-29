#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define DEBUG_NAN
//	#define FIXEDPREC	//not supported yet


//Xoroshiro128+ 1.0 by David Blackman and Sebastiano Vigna		https://prng.di.unimi.it/xoroshiro128plus.c
unsigned long long xoroshiro128_state[2]={0xDF900294D8F554A5, 0x170865DF4B3201FC};
static inline unsigned long long rotl(const unsigned long long x, int k)
{
	return (x<<k)|(x>>(64-k));
}
unsigned long long xoroshiro128_next(void)//not thread-safe, because of global state
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
#ifdef DEBUG_NAN
static float MUL(float A, float B)
{
	float p=A*B;
	if(!isfinite(p))
		LOG_ERROR("NaN");
	//p=CLAMP(-128, p, 127);
	//if(fabsf(p)>1e3)
	//	LOG_ERROR("Overflow");
	return p;
}
#else
#define MUL(A, B) ((A)*(B))
#endif
#define INITIALIZE(FIN, FOUT) (float)sample(FIN, FOUT)/0xFFFF
#define LOAD(PX) MUL(PX+0.5f, 1.f/128)		//[-128, 127] -> [-0.99609375, 0.99609375]
#define STORE(VAL) (char)floorf(MUL(VAL, 128))
//#define LOAD(PX) MUL(PX+0.5f, 1.f/256)
//#define STORE(VAL) (char)(MUL(VAL, 256)-0.5f)
//#define LEARNING_RATE(G) ((G)*0.001f)
#define ABS fabsf
#define FROMFLOAT(X) X
#endif

#define LEARNING_RATE(G, LR) MUL(CLAMP(-100, G, 100), FROMFLOAT(LR))
#define SGN(X) (DataType)((X)<0?-ONE:((X)>0?ONE:0))
static void initialize(DataType *w, int count, DataType sqrt_fan_in, DataType alpha0, DataType alpha1)
{
	for(int k=0;k<count;++k)
	{
#ifdef FIXEDPREC
		w[k]=(DataType)((xoroshiro128_next()&0x1FFFF)-0x10000)/sqrt_fan_in;
#else
		int x=(int)(xoroshiro128_next()&0x1FFFF)-0x10000;//[-2^16, 2^16]
		float val=(DataType)x/(0x10000*sqrt_fan_in);//[-2^-4, 2^-4]/sqrt_fan_in
		w[k]=w[k]*alpha0+val*alpha1;
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


//	#define RESNET
	#define RESNET2

#define NF0_MAIN	55		//55
#define NF0_ERRORS	39		//39
#define NF0_PX		0		//39
#define NF0 (NF0_MAIN+NF0_ERRORS+NF0_PX)
#ifdef RESNET
#define NF1 NF0
#define NF2 NF0
#define NF3 NF0
#else
#define NF1 128	//64
#define NF2 128	//32
#define NF3 94	//32
#endif
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
	DataType dL_dp, dL_dx3[NF3], dL_dx2[NF2], dL_dx1[NF1];
} BwdTemps;
typedef struct ParamsStruct
{
	DataType
		weight1[NF1*NF0], bias1[NF1],
		weight2[NF2*NF1], bias2[NF2],
		weight3[NF3*NF2], bias3[NF3],
		weight4[1*NF3];
} Params;
const int NPARAMS=sizeof(Params)/sizeof(DataType);
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
	
	//*ctx++ = clamp4(cL+cT-cTL, cL, cTL, cT, cTR);
	//*ctx++ = cL;
	//*ctx++ = cT;
	//*ctx++ = cTL;
	//*ctx++ = cTR;

#if 1
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
#endif
#if NF0_ERRORS==39
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
#endif
#if NF0_PX==39
	*ctx++=         ky-6>=0?LOAD(buf[idx-iw*6  ]):0,
	*ctx++=         ky-5>=0?LOAD(buf[idx-iw*5  ]):0,
	*ctx++=kx-3>=0&&ky-4>=0?LOAD(buf[idx-iw*4-3]):0,
	*ctx++=         ky-4>=0?LOAD(buf[idx-iw*4  ]):0,
	*ctx++=kx+3<iw&&ky-4>=0?LOAD(buf[idx-iw*4+3]):0,
	*ctx++=kx-5>=0&&ky-3>=0?LOAD(buf[idx-iw*3-5]):0,
	*ctx++=kx-4>=0&&ky-3>=0?LOAD(buf[idx-iw*3-4]):0,
	*ctx++=kx-2>=0&&ky-3>=0?LOAD(buf[idx-iw*3-2]):0,
	*ctx++=kx-1>=0&&ky-3>=0?LOAD(buf[idx-iw*3-1]):0,
	*ctx++=         ky-3>=0?LOAD(buf[idx-iw*3  ]):0,
	*ctx++=kx+1<iw&&ky-3>=0?LOAD(buf[idx-iw*3+1]):0,
	*ctx++=kx+2<iw&&ky-3>=0?LOAD(buf[idx-iw*3+2]):0,
	*ctx++=kx+3<iw&&ky-3>=0?LOAD(buf[idx-iw*3+3]):0,
	*ctx++=kx+4<iw&&ky-3>=0?LOAD(buf[idx-iw*3+4]):0,
	*ctx++=kx-3>=0&&ky-2>=0?LOAD(buf[idx-iw*2-3]):0,
	*ctx++=kx-2>=0&&ky-2>=0?LOAD(buf[idx-iw*2-2]):0,
	*ctx++=kx-1>=0&&ky-2>=0?LOAD(buf[idx-iw*2-1]):0,
	*ctx++=         ky-2>=0?LOAD(buf[idx-iw*2  ]):0,
	*ctx++=kx+1<iw&&ky-2>=0?LOAD(buf[idx-iw*2+1]):0,
	*ctx++=kx+2<iw&&ky-2>=0?LOAD(buf[idx-iw*2+2]):0,
	*ctx++=kx+3<iw&&ky-2>=0?LOAD(buf[idx-iw*2+3]):0,
	*ctx++=kx+4<iw&&ky-2>=0?LOAD(buf[idx-iw*2+4]):0,
	*ctx++=kx-3>=0&&ky-1>=0?LOAD(buf[idx-iw  -3]):0,
	*ctx++=kx-2>=0&&ky-1>=0?LOAD(buf[idx-iw  -2]):0,
	*ctx++=kx-1>=0&&ky-1>=0?LOAD(buf[idx-iw  -1]):0,
	*ctx++=kx  <iw&&ky-1>=0?LOAD(buf[idx-iw    ]):0,
	*ctx++=kx+1<iw&&ky-1>=0?LOAD(buf[idx-iw  +1]):0,
	*ctx++=kx+2<iw&&ky-1>=0?LOAD(buf[idx-iw  +2]):0,
	*ctx++=kx+3<iw&&ky-1>=0?LOAD(buf[idx-iw  +3]):0,
	*ctx++=kx+4<iw&&ky-1>=0?LOAD(buf[idx-iw  +4]):0,
	*ctx++=kx+5<iw&&ky-1>=0?LOAD(buf[idx-iw  +5]):0,
	*ctx++=kx+6<iw&&ky-1>=0?LOAD(buf[idx-iw  +6]):0,
	*ctx++=kx+7<iw&&ky-1>=0?LOAD(buf[idx-iw  +7]):0,
	*ctx++=kx-6>=0         ?LOAD(buf[idx     -6]):0,
	*ctx++=kx-5>=0         ?LOAD(buf[idx     -5]):0,
	*ctx++=kx-4>=0         ?LOAD(buf[idx     -4]):0,
	*ctx++=kx-3>=0         ?LOAD(buf[idx     -3]):0,
	*ctx++=kx-2>=0         ?LOAD(buf[idx     -2]):0,
	*ctx++=kx-1>=0         ?LOAD(buf[idx     -1]):0;
#endif
}
static void eval_fwd(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params const *p)
{
	get_nb2(src, errors, iw, ih, kx, ky, t->nb);

	//residuals
	linear(t->net1, p->weight1, t->nb, p->bias1, NF0, NF1);	//n1 = w1*nb + b1
	act(t->x1, t->net1, NF1);								//x1 = act(n1)
#ifdef RESNET
	add_vec(t->x1, t->x1, t->nb, NF0);						//x1r = x1 + nb		= act(w1*nb + b1) + nb												dx1r_dw1 = act'(n1)*nb		dx1r_db1 = act'(n1)
#endif

	linear(t->net2, p->weight2, t->x1, p->bias2, NF1, NF2);//n2 = w2*x1r + b2
	act(t->x2, t->net2, NF2);								//x2 = act(n2)
#ifdef RESNET
	add_vec(t->x2, t->x2, t->x1, NF0);						//x2r = x2 + x1r	= act(w2*x1r + b2) + x1r		dx2r_dx1r = act'(n2)*w2 + I			dx2r_dw2 = act'(n2)*x1r		dx2r_db2 = act'(n2)
#endif

	linear(t->net3, p->weight3, t->x2, p->bias3, NF2, NF3);	//n3 = w3*x2r + b3
	act(t->x3, t->net3, NF3);								//x3 = act(n3)
#ifdef RESNET2
	add_vec(t->x3, t->x3, t->nb, NF0);
#endif
#ifdef RESNET
	add_vec(t->x3, t->x3, t->x2, NF0);						//x3r = x3 + x2r	= act(w3*x2r + b3) + x2r		dx3r_dx2r = act'(n3)*w3 + I			dx3r_dw3 = act'(n3)*x2r		dx3r_db3 = act'(n3)		<- read right to left
#endif

	linear(&t->pred, p->weight4, t->x3, 0, NF3, 1);			//pred = w4*x3r				dL_dp = sgn(p-x)		dp_dx3r = w4						dp_dw4 = x3r
}
static void train(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params *p, BwdTemps *b, Params *g, float lr, int evalfwd, double *times)
{
	double t_start=time_ms(), t_now;
#define MARK_TIME(N) if(times)t_now=time_ms(), times[N]+=(t_now-t_start)*0.001, t_start=t_now
	kx=CLAMP(0, kx, iw-1);
	ky=CLAMP(0, ky, ih-1);
	//if((unsigned)kx>=(unsigned)iw||(unsigned)ky>=(unsigned)ih)
	//	return;
	
	//fwd
	//[NF1]	x1 = act(w1*nb + b1) + nb
	//[NF2]	x2 = act(w2*x1 + b2) + x1
	//[NF3]	x3 = act(w3*x2 + b3) + x2
	//[1]	pred = w4*x3
	//[1]	L = abs(pred-val)
	
	//bwd
	//[1]		dL_dp = sgn(p-x)
	//[NF3]		g_w4 = dL_dp * x3
	//
	//[NF3]		dL_dx3 = dL_dp * w4
	//[NF3]		g_b3 = dL_dx3 .* act'(n3)
	//[NF3*NF2]	g_w3 = g_b3 o x2
	//
	//[NF2]		dL_dx2 = dL_dx3 * (act'(n3) * w3 + I)		= g_b3 * w3 + dL_dx3
	//[NF2]		g_b2 = dL_dx2 * act'(n2)
	//[NF2*NF1]	g_w2 = g_b2 o x1
	//
	//[NF1]		dL_dx1 = dL_dx2 * (act'(n2) * w2 + I)		= g_b2 * w2 + dL_dx2
	//[NF1]		g_b1 = dL_dx1 * act'(n1)
	//[NF1*NF0]	g_w1 = g_b1 o nb


	//fwd
	if(evalfwd)
		eval_fwd(src, errors, iw, ih, kx, ky, t, p);

	//loss
	DataType curr=LOAD(src[iw*ky+kx]);
	t->expr=t->pred-curr;
	t->loss=ABS(t->expr);
	
	MARK_TIME(0);//fwd

	//bwd
	b->dL_dp=SGN(t->expr);//dL_dp = sgn(p-x)
	mul_vec_scalar(g->weight4, t->x3, b->dL_dp, NF3);//dL_dw4 = dL_dp * dp_dw4		= sgn(p-x) * x3T
	MARK_TIME(1);//bwd layer 4


	mul_vec_scalar(b->dL_dx3, p->weight4, b->dL_dp, NF3);
	act_dash(g->bias3, t->net3, NF3);
	mul_vec_ew(g->bias3, g->bias3, b->dL_dx3, NF3);
	mul_vec_outer(g->weight3, g->bias3, t->x2, NF3, NF2);
	MARK_TIME(2);//bwd layer 3


	mul_mat(b->dL_dx2, g->bias3, p->weight3, 1, NF3, NF2);
#ifdef RESNET
	add_vec(b->dL_dx2, b->dL_dx2, b->dL_dx3, NF2);
#endif
	act_dash(g->bias2, t->net2, NF2);
	mul_vec_ew(g->bias2, g->bias2, b->dL_dx2, NF2);
	mul_vec_outer(g->weight2, g->bias2, t->x1, NF2, NF1);
	MARK_TIME(3);//bwd layer 2


	mul_mat(b->dL_dx1, g->bias2, p->weight2, 1, NF2, NF1);
#ifdef RESNET
	add_vec(b->dL_dx1, b->dL_dx1, b->dL_dx2, NF1);
#endif
	act_dash(g->bias1, t->net1, NF1);
	mul_vec_ew(g->bias1, g->bias1, b->dL_dx1, NF1);
	mul_vec_outer(g->weight1, g->bias1, t->nb, NF1, NF0);
	MARK_TIME(4);//bwd layer 1

#if 0
	//dL_dn3 = dL_dp * dp_dx3 * dx3_dn3			= sgn(p-x) * w4T * act_dash(n3)
	act_dash(g->bias3, t->net3, NF3);
	mul_vec_ew(g->bias3, g->bias3, p->weight4, NF3);
	mul_vec_scalar(g->bias3, g->bias3, b->dL_dp, NF3);
	
	MARK_TIME(1);//bwd layer 3

	//memcpy(g->bias3, b->dL_dn3, NF3*sizeof(DataType));	//dL_db3 = dL_dn3
	mul_vec_outer(g->weight3, g->bias3, t->x2, NF3, NF2);	//dL_dw3 = dL_dn3 * dn3_dw3		= dL_dn3 * x2T (outer product)

	//dL_dn2 = mul_mat(dL_dn3T, w3)T .* act_dash(n2)
	mul_mat(g->bias2, g->bias3, p->weight3, 1, NF3, NF2);
	act_dash(b->actdash_n2, t->net2, NF2);
	mul_vec_ew(g->bias2, g->bias2, b->actdash_n2, NF2);
	
	MARK_TIME(2);//bwd layer 2

	//memcpy(g->bias2, b->dL_dn2, NF2*sizeof(DataType));	//dL_db2 = dL_dn2
	mul_vec_outer(g->weight2, g->bias2, t->x1, NF2, NF1);	//dL_dw2 = dL_dn2 * dn2_dw2		= dL_dn2 * x1T (outer product)

	//dL_dn1 = mul_mat(dL_dn2T, w2)T .* act_dash(n1)
	mul_mat(g->bias1, g->bias2, p->weight2, 1, NF2, NF1);
	act_dash(b->actdash_n1, t->net1, NF1);
	mul_vec_ew(g->bias1, g->bias1, b->actdash_n1, NF1);
	
	MARK_TIME(3);//bwd layer 1

	//memcpy(g->bias1, b->dL_dn1, NF1*sizeof(DataType));	//dL_db1 = dL_dn1
	mul_vec_outer(g->weight1, g->bias1, t->nb, NF1, NF0);	//dL_dw1 = dL_dn1 * dn1_dw1
	MARK_TIME(4);//bwd finish
#endif

	//update
	if(lr)
	{
		DataType *params=(DataType*)p, *gradient=(DataType*)g;
		for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
			params[k]-=LEARNING_RATE(gradient[k], lr);
			//params[k]-=LEARNING_RATE(gradient[k], 0.001f);
	}
	MARK_TIME(5);//update
#undef	MARK_TIME
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
static void initialize_all(Params *p, int reset_prng, float alpha0, float alpha1)
{
	if(reset_prng)
		XOROSHIRO128_RESET();
	initialize(p->weight1, _countof(p->weight1), (DataType)sqrtf(NF0), alpha0, alpha1);
	initialize(p->bias1  , _countof(p->bias1  ), (DataType)sqrtf(NF0), alpha0, alpha1);
	initialize(p->weight2, _countof(p->weight2), (DataType)sqrtf(NF1), alpha0, alpha1);
	initialize(p->bias2  , _countof(p->bias2  ), (DataType)sqrtf(NF1), alpha0, alpha1);
	initialize(p->weight3, _countof(p->weight3), (DataType)sqrtf(NF2), alpha0, alpha1);
	initialize(p->bias3  , _countof(p->bias3  ), (DataType)sqrtf(NF2), alpha0, alpha1);
	initialize(p->weight4, _countof(p->weight4), (DataType)sqrtf(NF3), alpha0, alpha1);
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
	memset(p, 0, sizeof(Params));
			
	//initialization
	//const int nparams=sizeof(Params)/sizeof(int);
	//DataType *ptr=(DataType*)p;
	//for(int k=0;k<nparams;++k)
	//	ptr[k]=sample;
	//	//ptr[k]=rand()&((ONE-1)>>FRACBITS/2);

	int pred, idx;
	char *buf2=fwd?src:dst, *errors=fwd?dst:src;
	float lr=0.001f, lr2;
	for(int kc=0;kc<3;++kc)
	{
		int reach=0, niter=0;
		switch(kc)
		{
		case 0:reach=5, niter=9;break;
		case 1:reach=5, niter=9;break;
		case 2:reach=5, niter=9;break;
		}
		initialize_all(p, 1, 0, 1);

		for(int k=0;k<res;++k)
			src[k]=buf[k<<2|kc];

		for(int ky=0;ky<ih;++ky)
		{
			//if(!(ky&7))
			{
				TimeInfo ti;
				parsetimedelta(time_ms()-t_start, &ti);
				set_window_title("%d/3, %d/%d, %.2lf%% - %02d-%02d-%06.3f", kc, ky, ih, 100.*(ih*kc+ky)/(ih*3), ti.hours, ti.mins, ti.secs);
			}

			for(int kx=0;kx<iw;++kx)
			{
				int niter2=niter+10/(1+ky);
				//int niter2=niter+(iw*ih*100)/(iw*ih+(iw*ky+kx)*500);
				//int niter2=kx==2&&ky==0?100:niter;
				//int niter2=kx==0&&ky==0?100:niter;
				//int niter2=niter;
				//if(ky>1)
				if(0)
				{
					for(int k=0;k<niter2;++k)//deterministic random coords
					{
						int rx, ry;
						do
						{
							rx=xoroshiro128_next()%(reach<<1|1)-reach;
							ry=xoroshiro128_next()%(reach+1)-reach;
						}while(!ry&&rx>=0);
						lr2=lr*2/((float)(rx*rx+ry*ry));
						train(buf2, errors, iw, ih, kx+rx, ky+ry, t, p, b, g, lr2, 1, 0);
					}
				}
				else
				{
					for(int k=0;k<niter2;++k)
					{
						for(int ky2=-reach;ky2<0;++ky2)
						{
							for(int kx2=-reach;kx2<=reach;++kx2)
							{
								lr2=lr/((float)(kx2*kx2+ky2*ky2));
								train(buf2, errors, iw, ih, kx+kx2, ky+ky2, t, p, b, g, lr2, 1, 0);
							}
						}
						for(int kx2=-reach;kx2<0;++kx2)
						{
							lr2=lr/abs(kx2);
							train(buf2, errors, iw, ih, kx+kx2, ky, t, p, b, g, lr2, 1, 0);
						}
					}
				}
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
#if 0
				//int niter=kc==1?1:6;
				for(int k=0;k<5;++k)
				{
					train(buf2, errors, iw, ih, kx-5, ky-5, t, p, b, g, lr/sqrtf(5*5+5*5), 1);
					train(buf2, errors, iw, ih, kx-4, ky-5, t, p, b, g, lr/sqrtf(4*4+5*5), 1);
					train(buf2, errors, iw, ih, kx-3, ky-5, t, p, b, g, lr/sqrtf(3*3+5*5), 1);
					train(buf2, errors, iw, ih, kx-2, ky-5, t, p, b, g, lr/sqrtf(2*2+5*5), 1);
					train(buf2, errors, iw, ih, kx-1, ky-5, t, p, b, g, lr/sqrtf(1*1+5*5), 1);
					train(buf2, errors, iw, ih, kx  , ky-5, t, p, b, g, lr/sqrtf(0*0+5*5), 1);
					train(buf2, errors, iw, ih, kx+1, ky-5, t, p, b, g, lr/sqrtf(1*1+5*5), 1);
					train(buf2, errors, iw, ih, kx+2, ky-5, t, p, b, g, lr/sqrtf(2*2+5*5), 1);
					train(buf2, errors, iw, ih, kx+3, ky-5, t, p, b, g, lr/sqrtf(3*3+5*5), 1);
					train(buf2, errors, iw, ih, kx+4, ky-5, t, p, b, g, lr/sqrtf(4*4+5*5), 1);
					train(buf2, errors, iw, ih, kx+5, ky-5, t, p, b, g, lr/sqrtf(5*5+5*5), 1);
				
					train(buf2, errors, iw, ih, kx-5, ky-4, t, p, b, g, lr/sqrtf(5*5+4*4), 1);
					train(buf2, errors, iw, ih, kx-4, ky-4, t, p, b, g, lr/sqrtf(4*4+4*4), 1);
					train(buf2, errors, iw, ih, kx-3, ky-4, t, p, b, g, lr/sqrtf(3*3+4*4), 1);
					train(buf2, errors, iw, ih, kx-2, ky-4, t, p, b, g, lr/sqrtf(2*2+4*4), 1);
					train(buf2, errors, iw, ih, kx-1, ky-4, t, p, b, g, lr/sqrtf(1*1+4*4), 1);
					train(buf2, errors, iw, ih, kx  , ky-4, t, p, b, g, lr/sqrtf(0*0+4*4), 1);
					train(buf2, errors, iw, ih, kx+1, ky-4, t, p, b, g, lr/sqrtf(1*1+4*4), 1);
					train(buf2, errors, iw, ih, kx+2, ky-4, t, p, b, g, lr/sqrtf(2*2+4*4), 1);
					train(buf2, errors, iw, ih, kx+3, ky-4, t, p, b, g, lr/sqrtf(3*3+4*4), 1);
					train(buf2, errors, iw, ih, kx+4, ky-4, t, p, b, g, lr/sqrtf(4*4+4*4), 1);
					train(buf2, errors, iw, ih, kx+5, ky-4, t, p, b, g, lr/sqrtf(5*5+4*4), 1);
				
					train(buf2, errors, iw, ih, kx-5, ky-3, t, p, b, g, lr/sqrtf(5*5+3*3), 1);
					train(buf2, errors, iw, ih, kx-4, ky-3, t, p, b, g, lr/sqrtf(4*4+3*3), 1);
					train(buf2, errors, iw, ih, kx-3, ky-3, t, p, b, g, lr/sqrtf(3*3+3*3), 1);
					train(buf2, errors, iw, ih, kx-2, ky-3, t, p, b, g, lr/sqrtf(2*2+3*3), 1);
					train(buf2, errors, iw, ih, kx-1, ky-3, t, p, b, g, lr/sqrtf(1*1+3*3), 1);
					train(buf2, errors, iw, ih, kx  , ky-3, t, p, b, g, lr/sqrtf(0*0+3*3), 1);
					train(buf2, errors, iw, ih, kx+1, ky-3, t, p, b, g, lr/sqrtf(1*1+3*3), 1);
					train(buf2, errors, iw, ih, kx+2, ky-3, t, p, b, g, lr/sqrtf(2*2+3*3), 1);
					train(buf2, errors, iw, ih, kx+3, ky-3, t, p, b, g, lr/sqrtf(3*3+3*3), 1);
					train(buf2, errors, iw, ih, kx+4, ky-3, t, p, b, g, lr/sqrtf(4*4+3*3), 1);
					train(buf2, errors, iw, ih, kx+5, ky-3, t, p, b, g, lr/sqrtf(5*5+3*3), 1);
				
					train(buf2, errors, iw, ih, kx-5, ky-2, t, p, b, g, lr/sqrtf(5*5+2*2), 1);
					train(buf2, errors, iw, ih, kx-4, ky-2, t, p, b, g, lr/sqrtf(4*4+2*2), 1);
					train(buf2, errors, iw, ih, kx-3, ky-2, t, p, b, g, lr/sqrtf(3*3+2*2), 1);
					train(buf2, errors, iw, ih, kx-2, ky-2, t, p, b, g, lr/sqrtf(2*2+2*2), 1);
					train(buf2, errors, iw, ih, kx-1, ky-2, t, p, b, g, lr/sqrtf(1*1+2*2), 1);
					train(buf2, errors, iw, ih, kx  , ky-2, t, p, b, g, lr/sqrtf(0*0+2*2), 1);
					train(buf2, errors, iw, ih, kx+1, ky-2, t, p, b, g, lr/sqrtf(1*1+2*2), 1);
					train(buf2, errors, iw, ih, kx+2, ky-2, t, p, b, g, lr/sqrtf(2*2+2*2), 1);
					train(buf2, errors, iw, ih, kx+3, ky-2, t, p, b, g, lr/sqrtf(3*3+2*2), 1);
					train(buf2, errors, iw, ih, kx+4, ky-2, t, p, b, g, lr/sqrtf(4*4+2*2), 1);
					train(buf2, errors, iw, ih, kx+5, ky-2, t, p, b, g, lr/sqrtf(5*5+2*2), 1);
				
					train(buf2, errors, iw, ih, kx-5, ky-1, t, p, b, g, lr/sqrtf(5*5+1*1), 1);
					train(buf2, errors, iw, ih, kx-4, ky-1, t, p, b, g, lr/sqrtf(4*4+1*1), 1);
					train(buf2, errors, iw, ih, kx-3, ky-1, t, p, b, g, lr/sqrtf(3*3+1*1), 1);
					train(buf2, errors, iw, ih, kx-2, ky-1, t, p, b, g, lr/sqrtf(2*2+1*1), 1);
					train(buf2, errors, iw, ih, kx-1, ky-1, t, p, b, g, lr/sqrtf(1*1+1*1), 1);
					train(buf2, errors, iw, ih, kx  , ky-1, t, p, b, g, lr/sqrtf(0*0+1*1), 1);
					train(buf2, errors, iw, ih, kx+1, ky-1, t, p, b, g, lr/sqrtf(1*1+1*1), 1);
					train(buf2, errors, iw, ih, kx+2, ky-1, t, p, b, g, lr/sqrtf(2*2+1*1), 1);
					train(buf2, errors, iw, ih, kx+3, ky-1, t, p, b, g, lr/sqrtf(3*3+1*1), 1);
					train(buf2, errors, iw, ih, kx+4, ky-1, t, p, b, g, lr/sqrtf(4*4+1*1), 1);
					train(buf2, errors, iw, ih, kx+5, ky-1, t, p, b, g, lr/sqrtf(5*5+1*1), 1);
				
					train(buf2, errors, iw, ih, kx-5, ky  , t, p, b, g, lr/5, 1);
					train(buf2, errors, iw, ih, kx-4, ky  , t, p, b, g, lr/4, 1);
					train(buf2, errors, iw, ih, kx-3, ky  , t, p, b, g, lr/3, 1);
					train(buf2, errors, iw, ih, kx-2, ky  , t, p, b, g, lr/2, 1);
					train(buf2, errors, iw, ih, kx-1, ky  , t, p, b, g, lr/1, 1);
				}
#endif

				//for(int k=0;k<1;++k)
				//{
				//	train(buf2, errors, iw, ih, kx-1, ky, t, p, b, g, 1);//CHEAT
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
				//train(buf2, iw, ih, kx, ky, t, p, b, g, 1);
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


//	#define SHOW_PRED
//	#define DEBUG_V2

static void check_fbuf(const float *buf, int count)
{
	for(int k=0;k<count;++k)
	{
		if(buf[k]!=buf[k])
			LOG_ERROR("NaN");
	}
}
static void update_weights(Params *p, Params *g, DataType lr)
{
	DataType *params=(DataType*)p, *gradient=(DataType*)g;
	
	const int count=sizeof(Params)/sizeof(DataType);
	int k;
	__m256 mlr=_mm256_set1_ps(lr);
	for(k=0;k<count-15;k+=16)
	{
		__m256 g0=_mm256_loadu_ps(gradient+k  );
		__m256 g1=_mm256_loadu_ps(gradient+k+8);
		__m256 p0=_mm256_loadu_ps(params+k  );
		__m256 p1=_mm256_loadu_ps(params+k+8);
		g0=_mm256_mul_ps(g0, mlr);
		g1=_mm256_mul_ps(g1, mlr);
		p0=_mm256_sub_ps(p0, g0);
		p1=_mm256_sub_ps(p1, g1);
		_mm256_storeu_ps(params+k  , p0);
		_mm256_storeu_ps(params+k+8, p0);
	}
	for(;k<count;++k)
		params[k]-=LEARNING_RATE(gradient[k], lr);

	//for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
	//	params[k]-=LEARNING_RATE(gradient[k], lr);
}
static void add_grad(Params *dst, Params const *src, float decay)
{
	const int count=sizeof(Params)/sizeof(DataType);
	DataType *ptr1=(DataType*)dst;
	DataType const *ptr2=(DataType const*)src;
	for(int k=0;k<count;++k)
		ptr1[k]+=ptr2[k]*decay;
}
static void param_weightav(Params *dst, Params const **inputs)
{
	const int count=sizeof(Params)/sizeof(DataType);
	float *fdst=(float*)dst;
	float weights[]=
	{//	N		W		NW		NE
		//1,      0,      0,      0,
		//0.5f,   0.3f,   0.1f,   0.1f,
		//0.45f,  0.45f,  0.05f,  0.05f,
		0.3f,   0.3f,   0.2f,   0.2f,
		//0.5f,   0.5f,   0,      0,
		//0,      1,      0,      0,
	};
	float wsum=0;
	for(int k=0;k<4;++k)
	{
		if(inputs[k])
			wsum+=weights[k];
	}
	if(!wsum)
	{
		initialize_all(dst, 1, 0, 1);
		return;
	}
	memset(fdst, 0, count*sizeof(DataType));
	wsum=0;
	for(int k=0;k<4;++k)
	{
		const float *input=(const float*)inputs[k];
		if(input)
		{
			for(int k2=0;k2<count;++k2)
				fdst[k2]+=input[k2]*weights[k];
		}
		else
			initialize_all(dst, 0, 1, weights[k]);
		wsum+=weights[k];
	}
	wsum=1/wsum;
	for(int k2=0;k2<count;++k2)
		fdst[k2]*=wsum;

	//const float
	//	*fN=(const float*)N,
	//	*fW=(const float*)W,
	//	*fNW=(const float*)NW,
	//	*fNE=(const float*)NE;
	//for(int k=0;k<count;++k)
	//	fdst[k]=0.3f*fN[k]+0.3f*fW[k]+0.2f*fNW[k]+0.2f*fNE[k];
}
void pred_learned_v2(char *buf, int iw, int ih, int fwd)
{
	double t_start=time_ms();
#if 0
	const int reaches[]=
	{
		5,
		5,
		5,
	};
	const int niters[]=
	{
		20,
		20,
		20,
	};
	int maxreach=reaches[0];
	for(int k=0;k<3;++k)
	{
		if(maxreach<reaches[k])
			maxreach=reaches[k];
	}
	maxreach+=5;
#endif
	//int gradbuflen=iw*ih*sizeof(Params);
	//int gradbuflen=2*iw*sizeof(Params);

	int res=iw*ih;
	char *src=(char*)malloc(res), *dst=(char*)malloc(res);
#ifdef SHOW_PRED
	char *residues=(char*)malloc(res);
#endif
	Temps *t=(Temps*)malloc(sizeof(Temps));
	Params *p=(Params*)malloc(sizeof(Params)*iw*2);
	//Params *p0=(Params*)malloc(sizeof(Params));
	BwdTemps *b=(BwdTemps*)malloc(sizeof(BwdTemps));
	Params *g=(Params*)malloc(sizeof(Params));
#ifdef DEBUG_V2
	Params *g_debug=(Params*)malloc(sizeof(Params));
	if(!g_debug)
	{
		LOG_ERROR("Allocation error");
		return;
	}
#endif
	if(!src||!dst
#ifdef SHOW_PRED
		||!residues
#endif
		||!t||!p||!b||!g)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	ArrayHandle title;
	STR_ALLOC(title, 1024);
	get_window_title(title->data, 1024);

	double times[3]={0};
	//int debug_x, debug_y;//

	float px_sum=0, pred_sum=0, error_sum=0;

	int pred, idx;
	char *buf2=fwd?src:dst, *errors=fwd?dst:src;
	float lr=0.005f;
	//float lr2;
	for(int kc=0;kc<3;++kc)
	{
		//int reach=reaches[kc], niter=niters[kc];
		//initialize_all(p0, 1, 0, 1);

		for(int k=0;k<res;++k)
			src[k]=buf[k<<2|kc];

		for(int ky=0;ky<ih;++ky)
		{
			//if(!(ky&4))
			{
				TimeInfo ti;
				parsetimedelta(time_ms()-t_start, &ti);
				set_window_title("%d/3, %d/%d, %.2lf%% - %02d-%02d-%06.3f,  %lf, %lf, %lf", kc, ky, ih, 100.*(ih*kc+ky)/(ih*3), ti.hours, ti.mins, ti.secs, times[0], times[1], times[2]);
			}

			double t1=time_ms(), t2;
			for(int kx=0;kx<iw;++kx)
			{
				//if(ky==(ih>>1)&&kx==(iw>>1))
				//	printf("");

				Params const *paramsets[]=
				{
					(unsigned)(ky-1)<(unsigned)ih?p+iw*((ky-1)&1)+kx  :0,//0 N
					(unsigned)(kx-1)<(unsigned)iw?p+iw*((ky  )&1)+kx-1:0,//1 W
					(unsigned)(kx-1)<(unsigned)iw&&(unsigned)(ky-1)<(unsigned)ih?p+iw*((ky-1)&1)+kx-1:0,//2 NW
					(unsigned)(kx+1)<(unsigned)iw&&(unsigned)(ky-1)<(unsigned)ih?p+iw*((ky-1)&1)+kx+1:0,//3 NE
				};
				Params *p_curr=p+iw*(ky&1)+kx;
				//param_weightav(p_curr, paramsets);

#if 0
				//memcpy(p, p0, sizeof(Params));
				//int niter2=niter+10/(1+ky);
				//int niter2=niter;
				
				if(kx>0)//learn from left grad
				{
					idx=iw*(ky&1)+kx-1;
					update_weights(p, g+idx, 0.01f);
				}
				if(ky>0)//learn from top grad
				{
					//MODVAR(ky3, ky3-1, gbufheight);
					idx=iw*((ky-1)&1)+kx;
					update_weights(p, g+idx, 0.01f);
				}
#endif
#if 0
				int gbufheight=maxreach;
				int ky3;
				for(int k=0;k<niter2;++k)
				{
					for(int ky2=-reach;ky2<0;++ky2)
					{
						if((unsigned)(ky+ky2)>=(unsigned)ih)
							continue;
						//ky3=ky+ky2;
						MODVAR(ky3, ky+ky2, gbufheight);
						for(int kx2=-reach;kx2<=reach;++kx2)
						{
							lr2=lr/((float)(kx2*kx2+ky2*ky2));
							//lr2=lr/sqrt((float)(kx2*kx2+ky2*ky2));
							if((unsigned)(kx+kx2)<(unsigned)iw)
							{
								idx=iw*ky3+kx+kx2;
								update_weights(p, g+idx, lr2);
#ifdef DEBUG_V2
								train(buf2, errors, iw, ih, kx+kx2, ky+ky2, t, p, b, g_debug, lr2, 1);//
								if(ky>reach)//
								{
									DataType *g2=(DataType*)g_debug, *g3=(DataType*)(g+idx);
									for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
									{
										if(fabsf(g2[k]-g3[k])>1e-6)
											LOG_ERROR("");
									}
								}
#endif

								//int LOL_1=(int)g[idx].weight1[0];//
								//debug_x=LOL_1%iw, debug_y=LOL_1/iw;
								//if(LOL_1&&(debug_x!=kx+kx2||debug_y!=ky3))
								//	LOG_ERROR("");
							}
						}
					}
					//ky3=ky;
					MODVAR(ky3, ky, gbufheight);
					for(int kx2=-reach;kx2<0;++kx2)
					{
						lr2=lr/(kx2*kx2);
						//lr2=lr/abs(kx2);
						if((unsigned)(kx+kx2)<(unsigned)iw)
						{
							idx=iw*ky3+kx+kx2;
							update_weights(p, g+idx, lr2);
#ifdef DEBUG_V2
							train(buf2, errors, iw, ih, kx+kx2, ky, t, p, b, g_debug, lr2, 1);//
							if(ky>reach)//
							{
								DataType *g2=(DataType*)g_debug, *g3=(DataType*)(g+idx);
								for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
								{
									if(fabsf(g2[k]-g3[k])>1e-4)
										LOG_ERROR("");
								}
							}
#endif

							//int LOL_1=(int)g[idx].weight1[0];//
							//debug_x=LOL_1%iw, debug_y=LOL_1/iw;
							//if(LOL_1&&(debug_x!=kx+kx2||debug_y!=ky3))
							//	LOG_ERROR("");
						}
					}
				}
#endif
#if 0
				if(ky>1)
				{
					for(int k=0;k<niter2;++k)//deterministic random coords
					{
						int rx, ry;
						do
						{
							rx=xoroshiro128_next()%(reach<<1|1)-reach;
							ry=xoroshiro128_next()%(reach+1)-reach;
						}while(!ry&&rx>=0);
						lr2=lr*2/((float)(rx*rx+ry*ry));
						train(buf2, errors, iw, ih, kx+rx, ky+ry, t, p, b, g, lr2, 1);
					}
				}
				else
				{
					for(int k=0;k<niter2;++k)
					{
						for(int ky2=-reach;ky2<0;++ky2)
						{
							for(int kx2=-reach;kx2<=reach;++kx2)
							{
								lr2=lr/((float)(kx2*kx2+ky2*ky2));
								train(buf2, errors, iw, ih, kx+kx2, ky+ky2, t, p, b, g, lr2, 1);
							}
						}
						for(int kx2=-reach;kx2<0;++kx2)
						{
							lr2=lr/abs(kx2);
							train(buf2, errors, iw, ih, kx+kx2, ky, t, p, b, g, lr2, 1);
						}
					}
				}
#endif
				t2=time_ms();
				times[0]+=(t2-t1)*0.001;
				t1=t2;


				float weights[]=
				{
					0.3f, 0.3f, 0.2f, 0.2f,
				};
				float fpreds[4], wsum=0, fpred;
				for(int k=0;k<4;++k)
				{
					if(paramsets[k])
					{
						eval_fwd(buf2, errors, iw, ih, kx, ky, t, paramsets[k]);
						fpreds[k]=t->pred;
						wsum+=weights[k];
					}
					else
					{
						fpreds[k]=0;
						weights[k]=0;
					}
				}
				fpred=0;
				if(wsum)
				{
					for(int k=0;k<4;++k)
						fpred+=weights[k]*fpreds[k];
					fpred/=wsum;
				}
				//fpred=0.3f*fpreds[0]+0.3f*fpreds[1]+0.2f*fpreds[2]+0.2f*fpreds[3];

				//eval_fwd(buf2, errors, iw, ih, kx, ky, t, p);

				t2=time_ms();
				times[1]+=(t2-t1)*0.001;
				t1=t2;

				pred=STORE(t->pred);
				pred=CLAMP(-128, pred, 127);

				idx=iw*ky+kx;
#ifdef SHOW_PRED
				residues[idx]=pred;
#endif
				if(fwd)
					dst[idx]=src[idx]-pred;
				else
					dst[idx]=src[idx]+pred;
				
				px_sum+=buf2[idx]*buf2[idx];
				pred_sum+=pred*pred;
				error_sum+=errors[idx]*errors[idx];

				float val0=LOAD(src[idx]);
				int bestmodel=0;
				for(int k=1;k<4;++k)
				{
					if(paramsets[bestmodel]&&fabsf(fpreds[bestmodel]-val0)>fabsf(fpreds[k]-val0))
						bestmodel=k;
				}
				if(paramsets[bestmodel])
					memcpy(p_curr, paramsets[bestmodel], sizeof(Params));
				else
					initialize_all(p_curr, 1, 0, 1);

				//if(kx==13)
				//	printf("");
				lr=0.001f;
				//float max_loss=0;
				for(int k=0;k<5;++k)			//<- hyperparameter
				{
					train(buf2, errors, iw, ih, kx, ky, t, p_curr, b, g, lr, 1, 0);
					//if(!max_loss||max_loss<t->loss)
					//	max_loss=t->loss;
					//lr*=t->loss/max_loss;
					lr*=0.95f;
					if(!t->loss)
						break;

					//check_fbuf((float*)t     , sizeof(*t     )/sizeof(float));
					//check_fbuf((float*)b     , sizeof(*b     )/sizeof(float));
					//check_fbuf((float*)g     , sizeof(*g     )/sizeof(float));
					//check_fbuf((float*)p_curr, sizeof(*p_curr)/sizeof(float));
				}
#if 0
				idx=iw*(ky&1)+kx;
				train(buf2, errors, iw, ih, kx, ky, t, p, b, g+idx, 0, 1);
				if(kx>0)//add gradient from left * decay
					add_grad(g+idx, g+idx-1, 0.5);
				if(ky>0)//add gradient from top * decay
				{
					//MODVAR(ky3, ky3-1, gbufheight);
					int idx2=iw*((ky-1)&1)+kx;
					add_grad(g+idx, g+idx2, 0.5);
				}
#endif
				//{
				//	DataType *g2=(DataType*)(g+idx);
				//	g2[0]=idx;
				//}
				t2=time_ms();
				times[2]+=(t2-t1)*0.001;
				t1=t2;
			}
		}

		for(int k=0;k<res;++k)
#ifdef SHOW_PRED
			buf[k<<2|kc]=residues[k];
#else
			buf[k<<2|kc]=dst[k];
#endif
	}
	free(src);
	free(dst);
#ifdef SHOW_PRED
	free(residues);
#endif
	free(t);
	free(p);
	//free(p0);
	free(b);
	free(g);
#ifdef DEBUG_V2
	free(g_debug);
#endif
#if 1
	{
		TimeInfo ti;
		parsetimedelta(time_ms()-t_start, &ti);
		set_window_title("3/3, %d/%d, 100%% - %02d-%02d-%06.3f,  %lf, %lf, %lf,  %f, %f, %f", ih, ih, ti.hours, ti.mins, ti.secs, times[0], times[1], times[2], sqrtf(px_sum/res), sqrtf(pred_sum/res), sqrtf(error_sum/res));
	}
#else
	set_window_title("%s", title->data);
#endif
	array_free(&title);
}



//	#define MARK_THREAD_DONE
	#define PRINT_PROGRESS_FROM_BASE_THREAD
//	#define VALIDATE_IDX
//	#define DISABLE_THREADS

#define NTHREADS 8
#define REACH 7
#define CONC 1		//(50+1000/(idx+1))
typedef struct ThreadArgsStruct
{
	int threadidx, kc;
#ifdef MARK_THREAD_DONE
	int done;
#endif
	//char *title;
	char *src, *dst, *buf2, *errors;
	float *ebuf2;
	int iw, ih;
	Temps *t;
	Params *p;
	BwdTemps *b;
	Params *g;
	int fwd;
	double t_start;
#ifdef VALIDATE_IDX
	int bounds[6];
#endif
	double times[6];
} ThreadArgs;
static unsigned pred_learned_v3_thread(void *args)
{
	double times[3]={0};
	ThreadArgs *a=(ThreadArgs*)args;
	int pred, idx, kx, ky;
	float *pred_errors[]=
	{
		a->ebuf2+a->iw*2*0,
		a->ebuf2+a->iw*2*1,
		a->ebuf2+a->iw*2*2,
		a->ebuf2+a->iw*2*3,
	};
	//double t_start=time_ms();
//#ifdef MARK_THREAD_DONE
//	for(int k=0;k<a->iw*2*sizeof(Params)/sizeof(float);++k)
//	{
//		float val=((float*)a->p)[k];
//		if(val)
//			LOG_ERROR("Buffer error");
//	}
//#endif
	for(ky=0;ky<a->ih;++ky)
	{
#if 0
		int offset_schedule[NTRAINS<<1]=
		{//	x	y
			-2, -2,
			-1, -2,
			 0, -2,
			 1, -2,
			 2, -2,

			 2, -1,
			 1, -1,
			 0, -1,
			-1, -1,
			-2, -1,
			
			 0,  0,

			-2,  0,
			-1,  0,

			 0,  0,
			
			 0,  0,
			 0, -1,
			-1,  0,
			 0,  0,
			 0, -1,
			-1,  0,
			 0,  0,
			 0, -1,
			-1,  0,
			 0,  0,
			 0,  0,
			 0, -1,
			-1,  0,
			 0,  0,
			 0,  0,
			 0,  0,
		};
#endif
		float weights[]=
		{
			0.3f, 0.3f, 0.2f, 0.2f,
		};
		float weights2[4];

		//if(!(ky&4))
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
#ifndef DISABLE_THREADS
		if(!a->threadidx)
#endif
#else
		if(0)
#endif
		{
			TimeInfo ti;
			parsetimedelta(time_ms()-a->t_start, &ti);
			set_window_title("%d/3, %d/%d, %.2lf%% - %02d-%02d-%06.3f,  %lf, %lf, %lf", a->kc, ky, a->ih, 100.*(a->ih*a->kc+ky)/(a->ih*3), ti.hours, ti.mins, ti.secs, times[0], times[1], times[2]);
		}

		double t1=time_ms(), t2;
		for(kx=0;kx<a->iw;++kx)
		{
			//if(ky==(ih>>1)&&kx==(iw>>1))
			//	printf("");
#ifdef VALIDATE_IDX
#define VALIDATE(IDX) ((unsigned)(IDX)>(unsigned)(a->iw*2)?LOG_ERROR("OOB"):IDX)
#else
#define VALIDATE(IDX) IDX
#endif

			Params const *paramsets[]=
			{
				(unsigned)(ky-1)<(unsigned)a->ih?a->p+VALIDATE(a->iw*((ky-1)&1)+kx  ):0,//0 N
				(unsigned)(kx-1)<(unsigned)a->iw?a->p+VALIDATE(a->iw*((ky  )&1)+kx-1):0,//1 W
				(unsigned)(kx-1)<(unsigned)a->iw&&(unsigned)(ky-1)<(unsigned)a->ih?a->p+VALIDATE(a->iw*((ky-1)&1)+kx-1):0,//2 NW
				(unsigned)(kx+1)<(unsigned)a->iw&&(unsigned)(ky-1)<(unsigned)a->ih?a->p+VALIDATE(a->iw*((ky-1)&1)+kx+1):0,//3 NE
			};
#ifdef VALIDATE_IDX
			if(a->iw*(ky&1)+kx>=a->bounds[3])
				LOG_ERROR("OOB");
#endif
			Params *p_curr=a->p+a->iw*(ky&1)+kx;
			t2=time_ms();
			times[0]+=(t2-t1)*0.001;
			t1=t2;


			float fpreds[4], wsum=0, fpred;
			for(int k=0;k<4;++k)
			{
				float *errors=pred_errors[k];
				if(paramsets[k])
				{
					eval_fwd(a->buf2, a->errors, a->iw, a->ih, kx, ky, a->t, paramsets[k]);
					fpreds[k]=a->t->pred;
					float
						NW=kx-1>=0   &&ky-1>=0?errors[a->iw*((ky-1)&1)+kx-1]:0,
						N =ky-1>=0            ?errors[a->iw*((ky-1)&1)+kx  ]:0,
						NE=kx+1<a->iw&&ky-1>=0?errors[a->iw*((ky-1)&1)+kx+1]:0;
					weights2[k]=weights[k]/(NW+N+NE+1.f/128);
					wsum+=weights2[k];
				}
				else
				{
					fpreds[k]=0;
					//weights[k]=0;
					weights2[k]=0;
				}
			}
			fpred=0;
			if(wsum)
			{
				for(int k=0;k<4;++k)
					fpred+=weights2[k]*fpreds[k];
				fpred/=wsum;
			}
			//fpred=0.3f*fpreds[0]+0.3f*fpreds[1]+0.2f*fpreds[2]+0.2f*fpreds[3];

			//eval_fwd(buf2, errors, iw, ih, kx, ky, t, p);

			t2=time_ms();
			times[1]+=(t2-t1)*0.001;
			t1=t2;

			pred=STORE(a->t->pred);
			pred=CLAMP(-128, pred, 127);

			idx=a->iw*ky+kx;
#ifdef SHOW_PRED
			residues[idx]=pred;
#endif
			if(a->fwd)
				a->dst[idx]=a->src[idx]-pred;
			else
				a->dst[idx]=a->src[idx]+pred;
				
			//px_sum+=buf2[idx]*buf2[idx];
			//pred_sum+=pred*pred;
			//error_sum+=errors[idx]*errors[idx];
			//if(ky)
			//	printf("");

			float val0=LOAD(a->src[idx]);
			int bestmodel=-1;
			for(int k=0;k<4;++k)
			{
				if(paramsets[k])
				{
					if((bestmodel==-1||fabsf(fpreds[bestmodel]-val0)>fabsf(fpreds[k]-val0)))
						bestmodel=k;
					weights[k]+=0.001f/fabsf(val0-fpreds[k]);
				}
				float *errors=pred_errors[k], e=fabsf(val0-fpreds[k]);
				errors[a->iw*(ky&1)+kx]=e;
				if(ky-1>=0&&kx+1<a->iw)//add curr to NE, such that N includes W (from jxl)
					errors[a->iw*((ky-1)&1)+kx+1]+=e;
			}

			//Params *p0=p_curr;
			if(bestmodel!=-1&&paramsets[bestmodel])
				memcpy(p_curr, paramsets[bestmodel], sizeof(Params));
			else if(kx||ky)
				LOG_ERROR("Algorithm error");
			//else
			//	initialize_all(p_curr, 0, 0, 1);

			float lr=0.0001f, lr2;
			int rx, ry;
			//int niter=20+(256-20)/(ky+1);
			for(ry=-REACH;ry<0;++ry)
			{
				for(rx=-REACH;rx<=REACH;++rx)
				{
					lr2=lr/(float)(rx*rx+ry*ry);
					train(a->buf2, a->errors, a->iw, a->ih, kx+rx, ky+ry, a->t, p_curr, a->b, a->g, lr2, 1, a->times);
				}
			}
			for(rx=-REACH;rx<0;++rx)
			{
				lr2=lr/(float)rx*rx;
				train(a->buf2, a->errors, a->iw, a->ih, kx+rx, ky, a->t, p_curr, a->b, a->g, lr2, 1, a->times);
			}
			lr2=lr;
			for(int k=0;k<CONC;++k)
				train(a->buf2, a->errors, a->iw, a->ih, kx, ky, a->t, p_curr, a->b, a->g, lr2, 1, a->times);
#if 0
			for(int k=0;k<NTRAINS;++k)			//<- hyperparameter
			{
				rx=offset_schedule[k<<1], ry=offset_schedule[k<<1|1];
				lr2=(float)(rx*rx+ry*ry);
				lr2+=!lr2*0.5f;
				lr2=lr/lr2;
				train(a->buf2, a->errors, a->iw, a->ih, kx+rx, ky+ry, a->t, p_curr, a->b, a->g, lr2, 1);
				//if(lr>0.00001)
				//	lr*=0.95f;
				if(a->t->loss<0.0001)
					break;
				//if(a->t->loss>0.05)
				//{
				//	if(lr<0.1)
				//		lr*=2;
				//	//--k;
				//}
			}
#endif
			t2=time_ms();
			times[2]+=(t2-t1)*0.001;
			t1=t2;
			//break;
		}
	}
#ifdef MARK_THREAD_DONE
	a->done=1;
#endif
#ifndef DISABLE_THREADS
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
	if(a->threadidx)
#endif
	{
		_endthreadex(0);
		//ExitThread(0);
		*(int*)0=0;//hit
	}
#endif
	return 0;
}
void pred_learned_v3(char *buf, int iw, int ih, int fwd)
{
	double t_start=time_ms();
	int res=iw*ih;
	char *src=(char*)malloc(res), *dst=(char*)malloc(res);
#ifdef SHOW_PRED
	char *residues=(char*)malloc(res);
#endif
	Temps *t=(Temps*)malloc(sizeof(Temps)*NTHREADS);
	Params *p=(Params*)malloc(sizeof(Params)*iw*2*NTHREADS);
	BwdTemps *b=(BwdTemps*)malloc(sizeof(BwdTemps)*NTHREADS);
	Params *g=(Params*)malloc(sizeof(Params)*NTHREADS);
	ThreadArgs *args=(ThreadArgs*)malloc(sizeof(ThreadArgs)*NTHREADS), *a;
	float *ebuf2=(float*)malloc(iw*2LL*4*sizeof(float)*NTHREADS);
	if(!src||!dst
#ifdef SHOW_PRED
		||!residues
#endif
		||!t||!p||!b||!g||!args||!ebuf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	ArrayHandle title;
	STR_ALLOC(title, 1024);
	get_window_title(title->data, 1024);

	//float px_sum=0, pred_sum=0, error_sum=0;

	char *buf2=fwd?src:dst, *errors=fwd?dst:src;
	float lr=0.005f;
	int blockh=(ih+NTHREADS-1)/NTHREADS;
	for(int kt=0;kt<NTHREADS;++kt)
	{
		a=args+kt;
		int offset=iw*blockh*kt;
		a->src   =src   +offset;
		a->dst   =dst   +offset;
		a->buf2  =buf2  +offset;
		a->errors=errors+offset;
		a->ebuf2=ebuf2+iw*2*4*kt;

		a->iw=iw;
		a->ih=blockh*(kt+1);
		if(a->ih>ih)
			a->ih=ih;
		a->ih-=blockh*kt;

		a->t=t+kt;
		a->p=p+iw*2*kt;
		a->b=b+kt;
		a->g=g+kt;
		a->fwd=fwd;
		a->t_start=t_start;
#ifdef VALIDATE_IDX
		a->bounds[0]=a->iw*a->ih;
		a->bounds[1]=a->iw*a->ih;
		a->bounds[2]=1;
		a->bounds[3]=iw*2;
		a->bounds[4]=1;
		a->bounds[5]=1;
#endif
		memset(a->times, 0, sizeof(a->times));
	}
	for(int kc=0;kc<3;++kc)
	{
		XOROSHIRO128_RESET();
		for(int k=0;k<res;++k)
			src[k]=buf[k<<2|kc];
		
		memset(ebuf2, 0, iw*2LL*4*sizeof(float)*NTHREADS);
		//memset(p, 0, sizeof(Params)*iw*2*NTHREADS);//

		HANDLE threads[NTHREADS]={0};
		for(int kt=0;kt<NTHREADS;++kt)
		{
			a=args+kt;
			a->threadidx=kt;
			a->kc=kc;
#ifdef MARK_THREAD_DONE
			a->done=0;
#endif
			initialize_all(a->p, 0, 0, 1);
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
			if(kt)
#endif
			{
#ifdef DISABLE_THREADS
				pred_learned_v3_thread(a);
#else
				threads[kt]=(HANDLE)_beginthreadex(0, 0, pred_learned_v3_thread, a, 0, 0);
				if(!threads[kt]||threads[kt]==(HANDLE)-1)
				{
					LOG_ERROR("Thread error %d", GetLastError());
					return;
				}
#endif
			}
		}
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
		pred_learned_v3_thread(args);//thread 0 is here, because SetWindowText doesn't work from another thread, and main thread is blocked
#endif
		
#ifndef DISABLE_THREADS
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
		int ret=WaitForMultipleObjects(NTHREADS-1, threads+1, TRUE, INFINITE);
#else
		int ret=WaitForMultipleObjects(NTHREADS, threads, TRUE, INFINITE);
#endif
		if(ret>=NTHREADS)
			LOG_ERROR("Thread error");
#ifdef PRINT_PROGRESS_FROM_BASE_THREAD
		for(int kt=1;kt<NTHREADS;++kt)
#else
		for(int kt=0;kt<NTHREADS;++kt)
#endif
		{
			int success=CloseHandle(threads[kt]);
			if(!success)
				LOG_ERROR("CloseHandle returned %d", GetLastError());//6 ERROR_INVALID_HANDLE
		}
#endif

		for(int k=0;k<res;++k)
#ifdef SHOW_PRED
			buf[k<<2|kc]=residues[k];
#else
			buf[k<<2|kc]=dst[k];
#endif
	}
#if 1
	{
		TimeInfo ti;
		parsetimedelta(time_ms()-t_start, &ti);
		set_window_title("3/3, %d/%d, 100%% - %02d-%02d-%06.3f,  fwd %lf, bwd4 %lf, bwd3 %lf, bwd2 %lf, bwd1 %lf, update %lf", ih, ih, ti.hours, ti.mins, ti.secs, args->times[0], args->times[1], args->times[2], args->times[3], args->times[4], args->times[5]);
	}
#else
	set_window_title("%s", title->data);
#endif
	array_free(&title);
	free(src);
	free(dst);
#ifdef SHOW_PRED
	free(residues);
#endif
	free(t);
	free(p);
	free(b);
	free(g);
	free(args);
	free(ebuf2);
}