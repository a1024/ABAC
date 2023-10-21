#include"best.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;

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

//T42: T39 with 'custom3' filter

//	#define T42_DISABLE_REC
//	#define T42_DISABLE_COUNTER
//	#define T42_PRINT_ESTIMATOR_CR

#define T42_LR (int)(0.07*0x10000+0.5)
#define T42_NMAPS 15

#define C3_REACH 3	//changing this requires re-training the filter
#define C3_NNB (C3_REACH*(C3_REACH+1)*4)
#define C3_NPARAMS (C3_NNB*9+6)

#ifndef T42_DISABLE_REC
#define T42_N_REC_ESTIMATORS 6		//15
#define T42_NESTIMATORS ((T42_N_REC_ESTIMATORS+1)*T42_NMAPS)
#else
#define T42_NESTIMATORS T42_NMAPS
#endif
typedef struct Custom3Struct//arbitrary size filter
{
	int reach, pred;
	short *params;
	char *errors;
} Custom3;
typedef struct Custom3ParamsStruct
{
	short c00[C3_NNB  ], c01[C3_NNB  ], c02[C3_NNB];//fixed 1.14 bit
	short c10[C3_NNB+2], c11[C3_NNB  ], c12[C3_NNB];
	short c20[C3_NNB+2], c21[C3_NNB+2], c22[C3_NNB];
} Custom3Params;
short filter[]=
{
	//CUSTOM3-r2 CLIC16
#if 0
	 0x00FE, 0x0398, 0x0001,-0x01E1,-0x01E5,-0x0E5E,-0x0019, 0x02FD,-0x000A, 0x0201,
	 0x0070, 0x03C6,-0x1CA5, 0x055A, 0x32FC,-0x041F,-0x0234, 0x0EEA, 0x0059, 0x010A,
	-0x009D,-0x0BB0, 0x2C77, 0x0407,
	
	-0x0001, 0x0028,-0x0001,-0x0014, 0x0003,-0x0043, 0x0001, 0x0001,-0x0005,-0x0001,
	 0x0001, 0x0001,-0x002C,-0x0024, 0x0028,-0x0012,-0x000C, 0x00AF, 0x0001,-0x0004,
	-0x0001, 0x0001, 0x0013, 0x0055,
	
	 0x0001, 0x0023,-0x0065, 0x0001,-0x003E, 0x00BA, 0x005C, 0x0034,-0x0076, 0x0130,
	 0x003B,-0x02A8, 0x00B9, 0x0294,-0x013B, 0x0486, 0x015E, 0x0075, 0x0025, 0x0027,
	-0x004A, 0x01EF,-0x009B, 0x0163,
	
	-0x00F6,-0x0049, 0x0003, 0x00B7, 0x0001, 0x0344,-0x0423,-0x0029, 0x0003, 0x0003,
	 0x0275,-0x00E2, 0x00C3,-0x0022, 0x00D4,-0x10BF, 0x0387,-0x05B6,-0x001F, 0x020F,
	-0x0007, 0x0165, 0x047A,-0x0673,-0x07B7,-0x2563,
	
	 0x0128, 0x0165,-0x01EF,-0x0273,-0x000B,-0x093C,-0x00B8, 0x02DA, 0x013D, 0x0126,
	 0x00AC, 0x010C,-0x0AEE, 0x049A, 0x1FA8, 0x1030, 0x03AE, 0x03E1, 0x00D5,-0x0002,
	 0x0147,-0x040C, 0x2575, 0x150A,
	
	 0x009B, 0x002E, 0x0005, 0x02EE, 0x00D8, 0x00CF, 0x0054, 0x00AC, 0x011C, 0x02EC,
	-0x0174,-0x0178, 0x00DF, 0x01A9,-0x00A3, 0x0C93,-0x0205, 0x083B, 0x0247,-0x03A3,
	 0x00F7, 0x01D7,-0x029B,-0x0239,
	
	-0x0057, 0x008A,-0x0003, 0x017C,-0x0041, 0x0099,-0x0003,-0x00D9,-0x0274, 0x0200,
	 0x0001, 0x00EB,-0x0073,-0x0018,-0x0065,-0x008F, 0x0440,-0x0458,-0x00F1, 0x0147,
	 0x00E4, 0x0299,-0x00D9,-0x0060, 0x005F,-0x0883,
	
	-0x0015,-0x0005, 0x0011, 0x0000, 0x0003,-0x0056, 0x0003,-0x0003,-0x00A3, 0x0023,
	-0x0001, 0x0001,-0x0059, 0x002D,-0x0019, 0x004F, 0x0042,-0x0087,-0x0003,-0x003E,
	-0x0003,-0x0001, 0x002F, 0x010E, 0x0085, 0x02C4,
	
	 0x0036, 0x012B,-0x0083, 0x00E0,-0x002F,-0x0359,-0x00D6, 0x0166, 0x0359,-0x00B4,
	 0x00D7, 0x0073,-0x0C88, 0x0E2C, 0x2021, 0x14BA, 0x02D0, 0x0C3C,-0x0034, 0x005D,
	 0x00B8,-0x04C2, 0x2481, 0x1435,
#endif

	//CUSTOM3-r3 kodim13
#if 1
	-0x0085, 0x0000, 0x000A, 0x0116, 0x0004, 0x01C3, 0x0008, 0x0176, 0x0002, 0x0143, 0x0025, 0x004B, 0x0007, 0x0026,
	 0x0034, 0x000F, 0x0342, 0x024E,-0x0B11, 0x0676, 0x113C,-0x048D,-0x1126, 0x06FC, 0x03F1,-0x03C2,-0x0005, 0x008F,
	-0x00B5,-0x006A,-0x0768, 0x041C,-0x0236, 0x0EDC, 0x0E5C, 0x0E9E, 0x1718,-0x01DB, 0x0117, 0x0210, 0x006C, 0x029A,
	 0x0011, 0x0001, 0x1035,-0x0069, 0x1667, 0x15E1,

	 0x0004, 0x0006, 0x0008, 0x000A,-0x000C,-0x0002,-0x0003, 0x0001,-0x0006,-0x0006, 0x0002, 0x0007, 0x0001, 0x0001,
	-0x0001,-0x0006,-0x0004, 0x0000, 0x0010, 0x0006, 0x0010, 0x0005,-0x000A, 0x0012, 0x0001, 0x0002,-0x0001,-0x0002,
	-0x0004,-0x0008, 0x001C, 0x0000,-0x0008, 0x002A, 0x000A, 0x0009, 0x0005, 0x0002,-0x0002,-0x0018,-0x0001,-0x0009,
	-0x0002, 0x0002,-0x000C, 0x000C, 0x0003,-0x0037,

	-0x0001, 0x0006,-0x0001,-0x000A,-0x0003,-0x0001,-0x0004, 0x0004, 0x0001, 0x000B, 0x0002, 0x002F, 0x0007,-0x0009,
	-0x0006, 0x004B, 0x0003, 0x005D, 0x000F,-0x0008,-0x0006, 0x0031,-0x0002, 0x0037, 0x0002, 0x00E7, 0x0000, 0x0005,
	-0x0023, 0x0004,-0x004C, 0x0016, 0x0056, 0x0031, 0x0004, 0x004F,-0x0003, 0x004A,-0x0007,-0x0026, 0x0005,-0x0004,
	 0x0008,-0x0002, 0x003F,-0x004B,-0x006D,-0x0069,

	-0x0095,-0x01A1,-0x0118,-0x014B, 0x012C,-0x024E, 0x0022, 0x038B, 0x00FC,-0x0030, 0x00B4,-0x0233,-0x0001,-0x0001,
	 0x0003,-0x042B, 0x0873,-0x02EF,-0x005A,-0x00C5, 0x03C6,-0x071E,-0x03A3,-0x00E5, 0x055E,-0x00E8, 0x0002, 0x0006,
	 0x000B, 0x0003,-0x006D,-0x000E,-0x0489,-0x0973,-0x03A3,-0x0F78,-0x0100,-0x01DA,-0x02E0, 0x0297, 0x0002,-0x0020,
	-0x0084, 0x01D0, 0x005D, 0x0171, 0x001C,-0x0AB5,-0x0554,-0x09B2,

	-0x000C, 0x0029,-0x0003,-0x0029, 0x0001,-0x0009, 0x0001, 0x02B2, 0x00A9,-0x0046,-0x0020,-0x018A,-0x006E,-0x0075,
	 0x0014, 0x0005,-0x0197,-0x00FE, 0x0286,-0x060E,-0x048B,-0x019F,-0x00F2,-0x01F4, 0x03EE,-0x0679, 0x0000,-0x0002,
	 0x0007,-0x0035, 0x0088,-0x0675, 0x097E, 0x01A4, 0x0769, 0x165E, 0x0274,-0x000A, 0x0ABF,-0x0642, 0x0001, 0x017B,
	 0x01D4, 0x00E1, 0x07C3,-0x06D2, 0x1730, 0x1838,

	 0x00BD, 0x0009, 0x0005,-0x00F0,-0x0008, 0x00E2,-0x0010, 0x0173, 0x006A,-0x0054,-0x0079,-0x00BC,-0x009D,-0x0029,
	 0x0159,-0x0205, 0x00F2,-0x0111, 0x01AA, 0x0085, 0x0187,-0x0533, 0x0266,-0x0045,-0x0005, 0x01AB,-0x0003, 0x012D,
	-0x0006, 0x0175,-0x00C1,-0x0134,-0x0154,-0x08F1, 0x0002,-0x0C83,-0x03D1,-0x0777, 0x0094,-0x019E, 0x0007,-0x00D0,
	-0x00B7, 0x016D, 0x0234,-0x0397,-0x038F,-0x053D,

	 0x01D1,-0x01B2,-0x02B3, 0x0152, 0x020F,-0x019E,-0x0138, 0x01B9, 0x0089, 0x00D6,-0x0100, 0x003A,-0x0035,-0x0008,
	-0x008B, 0x0098, 0x020E,-0x0066,-0x01B0, 0x0053, 0x02DB, 0x0068,-0x015C, 0x015B, 0x01B3,-0x0063,-0x00E0,-0x0105,
	-0x0006, 0x0000,-0x01E3,-0x01D1, 0x0392,-0x030B,-0x0340, 0x003E,-0x0003,-0x01BB, 0x0003,-0x0369, 0x0081,-0x00B6,
	-0x0069,-0x0079,-0x0012,-0x015A,-0x026E,-0x018A, 0x0216,-0x00A7,

	 0x0000,-0x0002,-0x0004,-0x0001, 0x0001, 0x0006,-0x0008, 0x0027, 0x0006, 0x0027,-0x0006, 0x0010,-0x0007,-0x0008,
	 0x0000, 0x0001, 0x0063,-0x0025, 0x0000,-0x0017, 0x0013,-0x0023, 0x0013,-0x0018, 0x005C,-0x0055,-0x0001, 0x000B,
	-0x0006, 0x0038,-0x0006, 0x0020, 0x0003, 0x0045,-0x001B, 0x0012,-0x0001,-0x0066, 0x0049,-0x003B,-0x0003, 0x0002,
	-0x0001, 0x0005, 0x0103,-0x0061,-0x0022, 0x00DC,-0x01B6, 0x00F3,

	 0x0008, 0x0003, 0x0003, 0x0083, 0x0000, 0x017D,-0x0005, 0x0059,-0x000C, 0x01B7,-0x0002,-0x0069, 0x0001,-0x0001,
	-0x005D, 0x00E9, 0x0251, 0x005D,-0x0B4E, 0x03DA, 0x09BF,-0x004F,-0x0268, 0x0469, 0x01B9,-0x011F, 0x000D, 0x00D4,
	-0x0007, 0x0270,-0x090C, 0x03D9,-0x0082, 0x09EE, 0x1819, 0x1572, 0x0436, 0x0C0E, 0x01EF, 0x009F, 0x00EE, 0x00D4,
	 0x0007, 0x00B4, 0x0B3C, 0x0214, 0x1E9A, 0x1552,
#endif
};
typedef struct T42NodeStruct
{
	int n[2];
#ifndef T42_DISABLE_REC
	unsigned short rec[T42_N_REC_ESTIMATORS];
#endif
} T42Node;
typedef struct T42CtxStruct
{
	int pred14;
	int context[T42_NMAPS];
	ArrayHandle maps[24][T42_NMAPS];//(256+512+1024+2048+4096+8192+16384+32768)*20*14 = 17.43 MB for 14 maps with 6 rec estimators
	T42Node *node[T42_NMAPS];

	int p0arr[T42_NESTIMATORS], p0_0, p0;//p0_0 isn't clamped
	int weights[24][T42_NESTIMATORS];
	long long wsum;

	int nnodes;
#ifdef T42_PRINT_ESTIMATOR_CR
	float csizes_est[24*T42_NESTIMATORS];
#endif
} T42Ctx;
T42Ctx* t42_ctx_init()
{
	int val=0x8000;
	T42Node node0={{1, 1}};
#ifndef T42_DISABLE_REC
	for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	T42Ctx *ctx=(T42Ctx*)malloc(sizeof(T42Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T42Ctx));
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T42_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T42Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T42Node), sizeof(T42Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t42_ctx_clear(T42Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T42_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
static int custom3_loadnb(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-C3_REACH;ky2<0;++ky2)
	{
		for(int kx2=-C3_REACH;kx2<=C3_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-C3_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static int custom3_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-15;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((__m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((__m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
	__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
	s2=_mm_add_epi32(s2, hi);
	s2=_mm_hadd_epi32(s2, s2);
	int s3=_mm_extract_epi32(s2, 0);
	for(;k<count;++k)
		s3+=a[k]*b[k];
	return s3;
}
void t42_ctx_get_context(T42Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	int
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

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);

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
	ctx->context[++j]=(N+W-NW + m1)>>1;//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=(N+W-NW + m2)>>1;//13

	{//CUSTOM3
		Custom3Params *params=(Custom3Params*)filter;
		const short *coeffs[]=
		{
			params->c00, params->c01, params->c02,
			params->c10, params->c11, params->c12,
			params->c20, params->c21, params->c22,

			//0, 0, 0,
			//0, 0, 0,
			//0, 0, 0,
		};
		short nb[3][C3_NNB+2]={0};
		int count[3], idx, idx2;
		for(int kc=0;kc<3;++kc)
			count[kc]=custom3_loadnb(buf, ebuf, iw, ih, kc, kx, ky, nb[kc]);
		idx=(iw*ky+kx)<<2;
		idx2=0;
		switch(kc)
		{
		case 1:
			nb[0][C3_NNB  ]=buf[idx];
			nb[0][C3_NNB+1]=ebuf[idx];
			count[0]+=2;
			//++idx;
			idx2+=3;
			break;
		case 2:
			nb[0][C3_NNB  ]=buf [idx  ];
			nb[0][C3_NNB+1]=ebuf[idx  ];
			nb[1][C3_NNB  ]=buf [idx|1];
			nb[2][C3_NNB+1]=ebuf[idx|1];
			count[0]+=2;
			count[1]+=2;
			//idx+=2;
			idx2+=6;
			break;
		}
		ctx->pred14=0;
		for(int kc=0;kc<3;++kc)
		{
			if(coeffs[idx2+kc])
				ctx->pred14+=custom3_dot(coeffs[idx2+kc], nb[kc], count[kc]);
		}

		ctx->pred14+=1<<13;
		ctx->pred14>>=14;
		ctx->pred14=CLAMP(-128, ctx->pred14, 127);
		ctx->context[++j]=ctx->pred14;
	}
#undef LOAD
	for(int k=0;k<T42_NMAPS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
int t42_ctx_map_context(int *context, int kp, int workidx)//replacement for context[kp]
{
	return context[kp];
}
void t42_ctx_estimate_p0(T42Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	T42Node *node;
	for(int kp=0;kp<T42_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		int context=t42_ctx_map_context(ctx->context, kp, kc);
		ArrayHandle map=ctx->maps[workidx][kp];
		//node=ctx->node[kp]=ARRAY_AT(T42Node, map, context);
		node=ctx->node[kp]=(T42Node*)array_at(&map, context);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T42_DISABLE_REC
		for(;k2<T42_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T42_NESTIMATORS;++k)
	{
#ifdef T42_DISABLE_COUNTER
		if(k%(T42_N_REC_ESTIMATORS+1))//
#endif
		{
			sum+=(long long)ctx->p0arr[k]*wk[k];
			ctx->wsum+=wk[k];
		}
	}
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t42_ctx_update(T42Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
	
#ifdef T42_PRINT_ESTIMATOR_CR
	for(int k=0;k<T42_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T42_NESTIMATORS*workidx+k]+=bitsize;
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
		for(int k=0;k<T42_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=T42_LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
	T42Node *node;
	for(int kp=0;kp<T42_NMAPS;++kp)
	{
		node=ctx->node[kp];
		++node->n[bit];
#ifndef T42_DISABLE_REC
		for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		{
			int lgden=k;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);
	}
}
int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T42 Enc  CUSTOM  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);
	T42Ctx *t42_ctx=t42_ctx_init();
	if(!buf2||!ebuf||!t42_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	memset(ebuf, 0, (size_t)res<<2);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd(buf2, iw, ih);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(t42_ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(t42_ctx, kc, kb);
					int bit=(buf2[idx]+128)>>kb&1;
					abac_enc(&ctx, t42_ctx->p0, bit);
					
					int prob=bit?0x10000-t42_ctx->p0:t42_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t42_ctx_update(t42_ctx, kc, kb, bit);
				}
				ebuf[idx]=buf2[idx]-t42_ctx->pred14;
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
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t42_ctx->nnodes*sizeof(T42Node)/(1024*1024));
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

#ifdef T42_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t42_ctx->csizes_est+T42_NESTIMATORS*kb;
				for(int ke=1;ke<T42_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T42_NESTIMATORS;++ke)
			{
				float *sizes=t42_ctx->csizes_est+ke;
#ifndef T42_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T42_N_REC_ESTIMATORS+1), ke%(T42_N_REC_ESTIMATORS+1));
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
					printf("%8.2f %c", iw*ih/sizes[T42_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T42_NESTIMATORS*kb]/t42_ctx->csizes_est[T42_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T42_DISABLE_REC
				if(!((ke+1)%(T42_N_REC_ESTIMATORS+1)))
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
	}
	t42_ctx_clear(&t42_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	//int debug_index=0;
	char *ebuf=(char*)malloc((size_t)res<<2);
	T42Ctx *t42_ctx=t42_ctx_init();
	if(!ebuf||!t42_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t42_ctx_init(t42_ctx);
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(t42_ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(t42_ctx, kc, kb);
					
					int bit=abac_dec(&ctx, t42_ctx->p0);
					buf[idx]|=bit<<kb;

					t42_ctx_update(t42_ctx, kc, kb, bit);
				}
				buf[idx]+=128;//unsigned -> signed
				ebuf[idx]=buf[idx]-t42_ctx->pred14;
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
		//if(!((ky+1)&127))
		//	t42_ctx_reset(&t42_ctx, 0);
	}
	t42_ctx_clear(&t42_ctx);
	
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}