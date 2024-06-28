#include"best.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<tmmintrin.h>
static const char file[]=__FILE__;

#ifndef assert
#define assert(X) (void)((X)||LOG_ERROR("%s", #X))
#endif


//T44	paq8pxd ported to C

#define ENABLE_SIMD
#define LOUD_UPDATE_PERIOD 16

#define T44_RCT 7

#if T44_RCT==0
	#define T44_RCT_FWD(...)
	#define T44_RCT_INV(...)
#elif T44_RCT==1
	#define T44_RCT_FWD colortransform_YCoCg_R_fwd
	#define T44_RCT_INV colortransform_YCoCg_R_inv
#elif T44_RCT==2
	#define T44_RCT_FWD colortransform_YCbCr_R_fwd
	#define T44_RCT_INV colortransform_YCbCr_R_inv
#elif T44_RCT==3
	#define T44_RCT_FWD colortransform_YCbCr_R_v2_fwd
	#define T44_RCT_INV colortransform_YCbCr_R_v2_inv
#elif T44_RCT==4
	#define T44_RCT_FWD colortransform_YCbCr_R_v3_fwd
	#define T44_RCT_INV colortransform_YCbCr_R_v3_inv
#elif T44_RCT==5
	#define T44_RCT_FWD colortransform_YCbCr_R_v4_fwd
	#define T44_RCT_INV colortransform_YCbCr_R_v4_inv
#elif T44_RCT==6
	#define T44_RCT_FWD colortransform_JPEG2000_fwd
	#define T44_RCT_INV colortransform_JPEG2000_inv
#elif T44_RCT==7
	#define T44_RCT_FWD colortransform_subgreen_fwd
	#define T44_RCT_INV colortransform_subgreen_inv
#endif

#define T44_LEVEL 8
#define T44_MEM (0x10000<<T44_LEVEL)
#define SQ(X) ((X)*(X))

static const unsigned char t44_State_table[256][4]=
{
	{  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7, 10, 2, 0}, // 0-3
	{  8, 12, 1, 1},{  9, 13, 1, 1},{ 11, 14, 0, 2},{ 15, 19, 3, 0}, // 4-7
	{ 16, 23, 2, 1},{ 17, 24, 2, 1},{ 18, 25, 2, 1},{ 20, 27, 1, 2}, // 8-11
	{ 21, 28, 1, 2},{ 22, 29, 1, 2},{ 26, 30, 0, 3},{ 31, 33, 4, 0}, // 12-15
	{ 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1}, // 16-19
	{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2}, // 20-23
	{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 36, 39, 1, 3},{ 36, 39, 1, 3}, // 24-27
	{ 36, 39, 1, 3},{ 36, 39, 1, 3},{ 38, 40, 0, 4},{ 41, 43, 5, 0}, // 28-31
	{ 42, 45, 4, 1},{ 42, 45, 4, 1},{ 44, 47, 3, 2},{ 44, 47, 3, 2}, // 32-35
	{ 46, 49, 2, 3},{ 46, 49, 2, 3},{ 48, 51, 1, 4},{ 48, 51, 1, 4}, // 36-39
	{ 50, 52, 0, 5},{ 53, 43, 6, 0},{ 54, 57, 5, 1},{ 54, 57, 5, 1}, // 40-43
	{ 56, 59, 4, 2},{ 56, 59, 4, 2},{ 58, 61, 3, 3},{ 58, 61, 3, 3}, // 44-47
	{ 60, 63, 2, 4},{ 60, 63, 2, 4},{ 62, 65, 1, 5},{ 62, 65, 1, 5}, // 48-51
	{ 50, 66, 0, 6},{ 67, 55, 7, 0},{ 68, 57, 6, 1},{ 68, 57, 6, 1}, // 52-55
	{ 70, 73, 5, 2},{ 70, 73, 5, 2},{ 72, 75, 4, 3},{ 72, 75, 4, 3}, // 56-59
	{ 74, 77, 3, 4},{ 74, 77, 3, 4},{ 76, 79, 2, 5},{ 76, 79, 2, 5}, // 60-63
	{ 62, 81, 1, 6},{ 62, 81, 1, 6},{ 64, 82, 0, 7},{ 83, 69, 8, 0}, // 64-67
	{ 84, 71, 7, 1},{ 84, 71, 7, 1},{ 86, 73, 6, 2},{ 86, 73, 6, 2}, // 68-71
	{ 44, 59, 5, 3},{ 44, 59, 5, 3},{ 58, 61, 4, 4},{ 58, 61, 4, 4}, // 72-75
	{ 60, 49, 3, 5},{ 60, 49, 3, 5},{ 76, 89, 2, 6},{ 76, 89, 2, 6}, // 76-79
	{ 78, 91, 1, 7},{ 78, 91, 1, 7},{ 80, 92, 0, 8},{ 93, 69, 9, 0}, // 80-83
	{ 94, 87, 8, 1},{ 94, 87, 8, 1},{ 96, 45, 7, 2},{ 96, 45, 7, 2}, // 84-87
	{ 48, 99, 2, 7},{ 48, 99, 2, 7},{ 88,101, 1, 8},{ 88,101, 1, 8}, // 88-91
	{ 80,102, 0, 9},{103, 69,10, 0},{104, 87, 9, 1},{104, 87, 9, 1}, // 92-95
	{106, 57, 8, 2},{106, 57, 8, 2},{ 62,109, 2, 8},{ 62,109, 2, 8}, // 96-99
	{ 88,111, 1, 9},{ 88,111, 1, 9},{ 80,112, 0,10},{113, 85,11, 0}, // 100-103
	{114, 87,10, 1},{114, 87,10, 1},{116, 57, 9, 2},{116, 57, 9, 2}, // 104-107
	{ 62,119, 2, 9},{ 62,119, 2, 9},{ 88,121, 1,10},{ 88,121, 1,10}, // 108-111
	{ 90,122, 0,11},{123, 85,12, 0},{124, 97,11, 1},{124, 97,11, 1}, // 112-115
	{126, 57,10, 2},{126, 57,10, 2},{ 62,129, 2,10},{ 62,129, 2,10}, // 116-119
	{ 98,131, 1,11},{ 98,131, 1,11},{ 90,132, 0,12},{133, 85,13, 0}, // 120-123
	{134, 97,12, 1},{134, 97,12, 1},{136, 57,11, 2},{136, 57,11, 2}, // 124-127
	{ 62,139, 2,11},{ 62,139, 2,11},{ 98,141, 1,12},{ 98,141, 1,12}, // 128-131
	{ 90,142, 0,13},{143, 95,14, 0},{144, 97,13, 1},{144, 97,13, 1}, // 132-135
	{ 68, 57,12, 2},{ 68, 57,12, 2},{ 62, 81, 2,12},{ 62, 81, 2,12}, // 136-139
	{ 98,147, 1,13},{ 98,147, 1,13},{100,148, 0,14},{149, 95,15, 0}, // 140-143
	{150,107,14, 1},{150,107,14, 1},{108,151, 1,14},{108,151, 1,14}, // 144-147
	{100,152, 0,15},{153, 95,16, 0},{154,107,15, 1},{108,155, 1,15}, // 148-151
	{100,156, 0,16},{157, 95,17, 0},{158,107,16, 1},{108,159, 1,16}, // 152-155
	{100,160, 0,17},{161,105,18, 0},{162,107,17, 1},{108,163, 1,17}, // 156-159
	{110,164, 0,18},{165,105,19, 0},{166,117,18, 1},{118,167, 1,18}, // 160-163
	{110,168, 0,19},{169,105,20, 0},{170,117,19, 1},{118,171, 1,19}, // 164-167
	{110,172, 0,20},{173,105,21, 0},{174,117,20, 1},{118,175, 1,20}, // 168-171
	{110,176, 0,21},{177,105,22, 0},{178,117,21, 1},{118,179, 1,21}, // 172-175
	{110,180, 0,22},{181,115,23, 0},{182,117,22, 1},{118,183, 1,22}, // 176-179
	{120,184, 0,23},{185,115,24, 0},{186,127,23, 1},{128,187, 1,23}, // 180-183
	{120,188, 0,24},{189,115,25, 0},{190,127,24, 1},{128,191, 1,24}, // 184-187
	{120,192, 0,25},{193,115,26, 0},{194,127,25, 1},{128,195, 1,25}, // 188-191
	{120,196, 0,26},{197,115,27, 0},{198,127,26, 1},{128,199, 1,26}, // 192-195
	{120,200, 0,27},{201,115,28, 0},{202,127,27, 1},{128,203, 1,27}, // 196-199
	{120,204, 0,28},{205,115,29, 0},{206,127,28, 1},{128,207, 1,28}, // 200-203
	{120,208, 0,29},{209,125,30, 0},{210,127,29, 1},{128,211, 1,29}, // 204-207
	{130,212, 0,30},{213,125,31, 0},{214,137,30, 1},{138,215, 1,30}, // 208-211
	{130,216, 0,31},{217,125,32, 0},{218,137,31, 1},{138,219, 1,31}, // 212-215
	{130,220, 0,32},{221,125,33, 0},{222,137,32, 1},{138,223, 1,32}, // 216-219
	{130,224, 0,33},{225,125,34, 0},{226,137,33, 1},{138,227, 1,33}, // 220-223
	{130,228, 0,34},{229,125,35, 0},{230,137,34, 1},{138,231, 1,34}, // 224-227
	{130,232, 0,35},{233,125,36, 0},{234,137,35, 1},{138,235, 1,35}, // 228-231
	{130,236, 0,36},{237,125,37, 0},{238,137,36, 1},{138,239, 1,36}, // 232-235
	{130,240, 0,37},{241,125,38, 0},{242,137,37, 1},{138,243, 1,37}, // 236-239
	{130,244, 0,38},{245,135,39, 0},{246,137,38, 1},{138,247, 1,38}, // 240-243
	{140,248, 0,39},{249,135,40, 0},{250, 69,39, 1},{ 80,251, 1,39}, // 244-247
	{140,252, 0,40},{249,135,41, 0},{250, 69,40, 1},{ 80,251, 1,40}, // 248-251
	{140,252, 0,41}//252, 253-255 are reserved
};
#define nex(state, sel) t44_State_table[state][sel]

typedef struct T44_PRNGStruct
{
	unsigned table[64];
	int i;
} T44_PRNG;
static void t44_rnd_init(T44_PRNG *rnd)
{
	rnd->table[0]=123456789;
	rnd->table[1]=987654321;
	for(int k=2;k<64;++k)
		rnd->table[k]=rnd->table[k-1]*11+rnd->table[k-2]*23/16;
	rnd->i=0;
}
static unsigned t44_rnd(T44_PRNG *rnd)
{
	++rnd->i;
	rnd->i&=63;
	rnd->table[rnd->i]=rnd->table[(rnd->i-24)&63]^rnd->table[(rnd->i-55)&63];
	return rnd->table[rnd->i];
}

static int t44_squash(int d)//p1 = squash(s = sum: wi*ti) = 1/(1+exp(-s))  (sigmoid)		clamp12(signed) -> uint12
{
	static const int t[33]=
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	if(d>2047)return 4095;
	if(d<-2047)return 0;
	int w=d&127;
	d=(d>>7)+16;
	return (t[d]*(128-w)+t[(d+1)]*w+64) >> 7;
}
static int t44_stretch(int p)//t = stretch(p1) = ln(p1/(1-p1))		uint12 -> signed int12
{
	static short t[4096];
	static int initialized=0;
	if(!initialized)
	{
		initialized=1;
		
		int pi=0;
		for(int x=-2047;x<=2047;++x)//invert squash()
		{
			int i=t44_squash(x);
			for(int j=pi;j<=i;++j)
				t[j]=x;
			pi=i+1;
		}
		t[4095]=2047;
	}
	return t[p];
}
static int t44_ilog(unsigned short x)
{
	static int initialized=0;
	static unsigned char t[65536];
	if(!initialized)
	{
		initialized=1;

		unsigned x=14155776;
		for(int i=2;i<65536;++i)
		{
			x+=774541002/(i*2-1);//numerator is 2^29/ln 2
			t[i]=x>>24;
		}
	}
	return t[x];
}
static int t44_dot_product(short *t, short *w, int n)
{
#ifdef ENABLE_SIMD
	__m128i sum=_mm_setzero_si128();

	//n=(n+7)&-8;
	while(n)
	{
		__m128i mt, mw;
		
		n-=8;
		mt=_mm_loadu_si128((__m128i*)(t+n));
		mw=_mm_loadu_si128((__m128i*)(w+n));
		mt=_mm_madd_epi16(mt, mw);
		mt=_mm_srai_epi32(mt, 8);
		sum=_mm_add_epi32(sum, mt);
	}
	sum=_mm_add_epi32(sum, _mm_srli_si128(sum, 8));
	sum=_mm_add_epi32(sum, _mm_srli_si128(sum, 4));
	return _mm_cvtsi128_si32(sum);
#else
	int sum=0;
	n=(n+7)&-8;
	for(int i=0;i<n;i+=2)
		sum+=(t[i]*w[i]+t[i+1]*w[i+1])>>8;
	return sum;
#endif
}
static unsigned t44_hash(unsigned a, unsigned b, unsigned c, unsigned d, unsigned e)
{
	unsigned h=a*200002979u+b*30005491u+c*50004239u+d*70004807u+e*110002499u;
	return h^h>>9^a>>2^b>>3^c>>4^d>>5^e>>6;
}
static void t44_train(short *t, short *w, int n, int err)
{
#ifdef ENABLE_SIMD
	__m128i merr=_mm_set1_epi16((short)err);
	__m128i one=_mm_set1_epi16(1);
	while(n)
	{
		__m128i mt, mw;

		n-=8;
		mt=_mm_loadu_si128((__m128i*)(t+n));
		mw=_mm_loadu_si128((__m128i*)(w+n));
		mt=_mm_adds_epi16(mt, mt);
		mt=_mm_mulhi_epi16(mt, merr);
		mt=_mm_adds_epi16(mt, one);
		mt=_mm_srai_epi16(mt, 1);
		mw=_mm_adds_epi16(mw, mt);
		_mm_storeu_si128((__m128i*)(w+n), mw);
	}
#else
	n=(n+7)&-8;
	for(int i=0;i<n;++i)
	{
		int wt=w[i]+(((t[i]*err*2>>16)+1)>>1);
		wt=CLAMP(-32768, wt, 32767);
		w[i]=wt;
	}
#endif
}
static int t44_dt(int x)//x -> 16K/(x+3)
{
	static int initialized=0;
	static int dt[1024];
	if(!initialized)
	{
		for(int i=0;i<1024;++i)
			dt[i]=16384/(i+i+3);
		initialized=1;
	}
	return dt[x&1023];
}

typedef struct T44MixerStruct
{
	int N, M, S;//max inputs, max contexts, max context sets
	int ncxt;//number of contexts (0 to S)
	int base;//offset of next context
	int nx;//Number of inputs in tx, 0 to N
	ArrayHandle tx;//<short>	N inputs from add()
	ArrayHandle wx;//<short>	N*M weights
	ArrayHandle cxt;//<int>		S contexts
	ArrayHandle pr;//<int>		last result (scaled 12 bits)
	struct T44MixerStruct *mp;//points to a Mixer to combine results
} T44Mixer;
static T44Mixer* t44_mixer_alloc(int n, int m, int s, short w)
{
	T44Mixer *mixer=(T44Mixer*)malloc(sizeof(T44Mixer));
	if(!mixer)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	mixer->N=(n+7)&-8;
	mixer->M=m;
	mixer->S=s;
	ARRAY_ALLOC(short, mixer->tx, 0, mixer->N, 0, 0);
	int msize=mixer->N*mixer->M;
	mixer->wx=array_construct(&w, sizeof(short), 1, msize, 0, 0);
	ARRAY_ALLOC(int, mixer->cxt, 0, mixer->S, 0, 0);
	mixer->ncxt=0;
	mixer->base=0;
	mixer->nx=0;
	const int val=2048;
	mixer->pr=array_construct(&val, sizeof(int), 1, mixer->S, 0, 0);
	if(mixer->S>1)
		mixer->mp=t44_mixer_alloc(mixer->S, 1, 1, 0);
	else
		mixer->mp=0;
	return mixer;
}
static void t44_mixer_free(T44Mixer *mixer)
{
	if(mixer)
	{
		array_free(&mixer->tx);
		array_free(&mixer->wx);
		array_free(&mixer->cxt);
		array_free(&mixer->pr);
		t44_mixer_free(mixer->mp);
		free(mixer);
	}
}
static void t44_mixer_update(T44Mixer *mixer, int bit)
{
	short *tx=(short*)array_at(&mixer->tx, 0);
	short *wx=(short*)array_at(&mixer->wx, 0);
	int *pr=(int*)array_at(&mixer->pr, 0);
	int *cxt=(int*)array_at(&mixer->cxt, 0);
	for(int i=0;i<mixer->ncxt;++i)
	{
		int err=((bit<<12)-pr[i])*7;
		assert(err>=-32768 && err<32768);
		if(err)
			t44_train(&tx[0], &wx[cxt[i]*mixer->N], mixer->nx, err);
	}
	mixer->nx=mixer->base=mixer->ncxt=0;
}
static void t44_mixer_add(T44Mixer *mixer, int x)
{
	short *p=(short*)array_at(&mixer->tx, mixer->nx);
	*p=x;
	++mixer->nx;
}
static void t44_mixer_set(T44Mixer *mixer, int cx, int range)
{
	assert(range>=0);
	assert(mixer->ncxt<mixer->S);
	assert(cx>=0);
	assert(mixer->base+cx<mixer->M);
	int *cxt=(int*)array_at(&mixer->cxt, 0);
	cxt[mixer->ncxt++]=mixer->base+cx;
	mixer->base+=range;
}
static int t44_mixer_predict(T44Mixer *mixer, int bit)
{
	short *tx=(short*)array_at(&mixer->tx, 0);
	short *wx=(short*)array_at(&mixer->wx, 0);
	int *cxt=(int*)array_at(&mixer->cxt, 0);
	int *pr=(int*)array_at(&mixer->pr, 0);

	while(mixer->nx&7)
		tx[mixer->nx++]=0;//pad
	if(mixer->mp)//combine outputs
	{
		t44_mixer_update(mixer->mp, bit);
		for(int i=0;i<mixer->ncxt;++i)
		{
			pr[i]=t44_squash(t44_dot_product(&tx[0], &wx[cxt[i]*mixer->N], mixer->nx)>>5);
			t44_mixer_add(mixer->mp, t44_stretch(pr[i]));
		}
		t44_mixer_set(mixer->mp, 0, 1);
		return t44_mixer_predict(mixer->mp, bit);
	}
	return pr[0]=t44_squash(t44_dot_product(&tx[0], &wx[0], mixer->nx)>>8);//S=1 context
}
typedef struct T44_APM1Struct
{
	int index;//last p, context
	int N;//number of contexts
	ArrayHandle t;//<unsigned short>	[N][33]:  p, context -> p
} T44_APM1;
static void t44_apm1_init(T44_APM1 *apm, int n)
{
	apm->index=0;
	apm->N=n;
	ARRAY_ALLOC(unsigned short, apm->t, 0, 33*n, 0, 0);
	unsigned short *ptr=(unsigned short*)apm->t->data;
	for(int ky=0;ky<n;++ky)
	{
		for(int kx=0;kx<33;++kx)
			ptr[33*ky+kx]=!ky?t44_squash((kx-16)<<7)<<4:ptr[kx];
	}
}
static void t44_apm1_clear(T44_APM1 *apm)
{
	array_free(&apm->t);
}
static int t44_apm1_predict(T44_APM1 *apm, int pr, int cxt, int rate, int bit)
{
	unsigned short *t=(unsigned short*)array_at(&apm->t, 0);
	assert(pr>=0 && pr<4096 && cxt>=0 && cxt<apm->N && rate>0 && rate<32);
	pr=t44_stretch(pr);
	int g=(bit<<16)+(bit<<rate)-bit-bit;
	t[apm->index]+=(g-t[apm->index])>>rate;
	t[apm->index+1]+=(g-t[apm->index+1])>>rate;
	const int w=pr&127;  // interpolation weight (33 points)
	apm->index=((pr+2048)>>7)+cxt*33;
	return (t[apm->index]*(128-w)+t[apm->index+1]*w) >> 11;
}
typedef struct T44SmallStationaryContextMapStruct
{
	ArrayHandle t;//<unsigned short>
	int cxt;
	unsigned short *cp;
} T44SmallStationaryContextMap;
static void t44_scm_init(T44SmallStationaryContextMap *scm, int m)
{
	const unsigned short val=0x8000;
	scm->t=array_construct(&val, sizeof(unsigned short), 1, m>>1, 0, 0);
	scm->cxt=0;
	scm->cp=(unsigned short*)scm->t->data;
}
static void t44_scm_clear(T44SmallStationaryContextMap *scm)
{
	array_free(&scm->t);
}
static void t44_scm_set(T44SmallStationaryContextMap *scm, unsigned cx)
{
	scm->cxt=cx<<8&(int)(scm->t->count-256);
}
static void t44_scm_mix(T44SmallStationaryContextMap *scm, T44Mixer *m, int rate, int c0, int bit)
{
	unsigned short *t=(unsigned short*)array_at(&scm->t, 0);
	*scm->cp+=((bit<<16)-*scm->cp+(1<<rate>>1))>>rate;
	scm->cp=t+scm->cxt+c0;
	t44_mixer_add(m, t44_stretch((*scm->cp)>>4));
}
typedef struct T44StateMapStruct
{
	int N;//Number of contexts
	int cxt;//Context of last prediction
	ArrayHandle t;//<unsigned>	cxt -> prediction in high 22 bits, count in low 10 bits
} T44StateMap;
static void t44_sm_init(T44StateMap *sm, int n)
{
	sm->N=n;
	sm->cxt=0;
	unsigned val=1<<31;
	sm->t=array_construct(&val, sizeof(unsigned), 1, n, 0, 0);
}
static void t44_sm_clear(T44StateMap *sm)
{
	array_free(&sm->t);
}
static void t44_sm_update(T44StateMap *sm, int limit, int bit)
{
	unsigned *t=(unsigned*)array_at(&sm->t, 0);
	assert(sm->cxt>=0 && sm->cxt<sm->N);
	unsigned *p=&t[sm->cxt], p0=p[0];
	int n=p0&1023, pr=p0>>10;//count, prediction
	if(n<limit)
		++p0;
	else
		p0=(p0&0xfffffc00)|limit;
	p0+=(((bit<<22)-pr)>>3)*t44_dt(n)&0xFFFFFC00;
	p[0]=p0;
}
static int t44_sm_predict(T44StateMap *sm, int cx, int limit, int bit)
{
	assert(cx>=0 && cx<sm->N);
	assert(limit>0 && limit<1024);
	t44_sm_update(sm, limit, bit);
	unsigned *t=(unsigned*)array_at(&sm->t, 0);
	return t[sm->cxt=cx]>>20;
}
typedef struct T44HashElemStruct//hash element, 14+1+49=64 bytes
{
	unsigned short chk[7];//byte context checksums
	unsigned char last;//last 2 accesses (0-6) in low, high nibble
	unsigned char bh[7][7];//byte context, 3-bit context -> bit history state
	//bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
	//bh[][0] is also a replacement priority, 0 = empty
} T44HashElem;
static unsigned char* t44_hash_get(T44HashElem *h, unsigned short ch, int n)//Find element (0-6) matching checksum. If not found, insert or replace lowest priority (not last).
{
	//Find or create hash element matching checksum ch
	int b=0xFFFF, bi=0;
	if(h->chk[h->last&15]==ch)
		return &h->bh[h->last&15][0];
	for(int i=0;i<7;++i)
	{
		if(h->chk[i]==ch)
		{
			h->last=h->last<<4|i;
			return &h->bh[i][0];
		}
		int pri=h->bh[i][0];
		if(pri<b && (h->last&15)!=i && h->last>>4!=i)
		{
			b=pri;
			bi=i;
		}
	}
	h->last=0xf0|bi;
	h->chk[bi]=ch;
	return (unsigned char*)memset(&h->bh[bi][0], 0, 7);
}
typedef struct T44ContextMapStruct
{
	int C;//max number of contexts
	int cn;//Next context to set by set()
	ArrayHandle t;//<T44HashElem>		bit histories for bits 0-1, 2-4, 5-7	For 0-1, also contains a run count in bh[][4] and value in bh[][5] and pending update count in bh[7]
	ArrayHandle cp;//<unsigned char*>	C pointers to current bit history
	ArrayHandle cp0;//<unsigned char*>	First element of 7 element array containing cp[i]
	ArrayHandle cxt;//<unsigned>		C whole byte contexts (hashes)
	ArrayHandle runp;//<unsigned char*>	C [0..3] = count, value, unused, unused
	ArrayHandle sm;//<StateMap>		C maps of state -> p
} T44ContextMap;
static void t44_cm_init(T44ContextMap *cm, int m, int c)//Construct using m bytes of memory for c contexts
{
	cm->C=c;
	cm->cn=0;
	ARRAY_ALLOC(T44HashElem, cm->t, 0, m>>6, 0, 0);
	ARRAY_ALLOC(unsigned char*, cm->cp, 0, c, 0, 0);
	ARRAY_ALLOC(unsigned char*, cm->cp0, 0, c, 0, 0);
	ARRAY_ALLOC(unsigned, cm->cxt, 0, c, 0, 0);
	ARRAY_ALLOC(unsigned char*, cm->runp, 0, c, 0, 0);
	ARRAY_ALLOC(T44StateMap, cm->sm, 0, c, 0, 0);
	T44StateMap *sm=(T44StateMap*)array_at(&cm->sm, 0);
	for(int k=0;k<c;++k)
		t44_sm_init(sm+k, 256);

	unsigned char
		**cp=(unsigned char**)cm->cp->data,
		**cp0=(unsigned char**)cm->cp0->data,
		**runp=(unsigned char**)cm->runp->data;
	T44HashElem *h=(T44HashElem*)cm->t->data;
	for(int k=0;k<c;++k)
	{
		cp0[k]=cp[k]=h->bh[0];
		runp[k]=cp[k]+3;
	}
}
static void t44_cm_clear(T44ContextMap *cm)
{
	array_free(&cm->t);
	array_free(&cm->cp);
	array_free(&cm->cp0);
	array_free(&cm->cxt);
	array_free(&cm->runp);
	T44StateMap *sm=(T44StateMap*)array_at(&cm->sm, 0);
	for(int k=0;k<(int)cm->sm->count;++k)
		t44_sm_clear(sm+k);
	array_free(&cm->sm);
}
static void t44_cm_set(T44ContextMap *cm, unsigned cx)
{
	int i=cm->cn;
	++cm->cn;
	unsigned *p=(unsigned*)array_at(&cm->cxt, i);
	cx=cx*987654323+i;//permute (don't hash) cx to spread the distribution
	cx=cx<<16|cx>>16;
	*p=cx*123456791+i;
}
static int t44_mix2(T44Mixer *m, int s, T44StateMap *sm, int bit)
{
	int p1=t44_sm_predict(sm, s, 1023, bit);
	int n0=-!nex(s, 2);
	int n1=-!nex(s, 3);
	int st=t44_stretch(p1)>>2;
	t44_mixer_add(m, st);
	p1>>=4;
	int p0=255-p1;
	t44_mixer_add(m, p1-p0);
	t44_mixer_add(m, st*(n1-n0));
	t44_mixer_add(m, (p1&n0)-(p0&n1));
	t44_mixer_add(m, (p1&n1)-(p0&n0));
	return s>0;
}
static int t44_cm_mix(T44ContextMap *cm, T44_PRNG *rnd, T44Mixer *m, int bpos, int prev, int c0, int bit, int debug_idx)
{
	int result=0;
	for(int i=0;i<cm->cn;++i)
	{
		unsigned char **cp=(unsigned char**)array_at(&cm->cp, 0);
		unsigned char **cp0=(unsigned char**)array_at(&cm->cp0, 0);
		unsigned *cxt=(unsigned*)array_at(&cm->cxt, 0);
		unsigned char **runp=(unsigned char**)array_at(&cm->runp, 0);
		T44HashElem *t=(T44HashElem*)array_at(&cm->t, 0);
		int tsize_m1=(int)cm->t->count-1;
		T44StateMap *sm=(T44StateMap*)array_at(&cm->sm, 0);

		if(cp[i])
		{
			assert(cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tsize_m1].bh[6][6]);
			//assert(((size_t)cp[i]&63)>=15);
			int ns=nex(*cp[i], bit);
			if(ns>=204 && t44_rnd(rnd)<<((452-ns)>>3))//probabilistic increment
				ns-=4;
			*cp[i]=ns;
		}

		//Update context pointers
		if(bpos>1 && runp[i][0]==0)
			cp[i]=0;
		else
		{
			switch(bpos)
			{
			case 1:case 3:case 6:
				cp[i]=cp0[i]+1+(c0&1);
				break;
			case 4:case 7:
				cp[i]=cp0[i]+3+(c0&3);
				break;
			case 2:case 5:
				cp0[i]=cp[i]=t44_hash_get(t+((cxt[i]+c0)&tsize_m1), cxt[i]>>16, 0);
				break;
			default:
				cp0[i]=cp[i]=t44_hash_get(t+((cxt[i]+c0)&tsize_m1), cxt[i]>>16, sm[i].N);

				//Update pending bit histories for bits 2-7
				if(cp0[i][3]==2)
				{
					const int c=cp0[i][4]+256;
					unsigned char *p=t44_hash_get(t+((cxt[i]+(c>>6))&tsize_m1), cxt[i]>>16, 0);
					p[0]=1+((c>>5)&1);
					p[1+((c>>5)&1)]=1+((c>>4)&1);
					p[3+((c>>4)&3)]=1+((c>>3)&1);
					p=t44_hash_get(t+((cxt[i]+(c>>3))&tsize_m1), cxt[i]>>16, 0);
					p[0]=1+((c>>2)&1);
					p[1+((c>>2)&1)]=1+((c>>1)&1);
					p[3+((c>>1)&3)]=1+(c&1);
					cp0[i][6]=0;
				}
				//Update run count of previous context
				if(runp[i][0]==0)  // new context
					runp[i][0]=2, runp[i][1]=prev;
				else if(runp[i][1]!=prev)  // different byte in context
					runp[i][0]=1, runp[i][1]=prev;
				else if(runp[i][0]<254)  // same byte in context
					runp[i][0]+=2;
				else if(runp[i][0]==255)
					runp[i][0]=128;
				runp[i]=cp0[i]+3;
				break;
			}
		}

		//predict from last byte in context
		if((runp[i][1]+256)>>(8-bpos)==c0)
		{
			int rc=runp[i][0];//count*2, +1 if 2 different bytes seen
			int b=(runp[i][1]>>(7-bpos)&1)*2-1;//predicted bit, + for 1, - for 0
			int c=t44_ilog(rc+1)<<(2+(~rc&1));
			t44_mixer_add(m, b*c);
		}
		else
			t44_mixer_add(m, 0);

		//predict from bit context
		if(cp[i])
			result+=t44_mix2(m, *cp[i], sm+i, bit);
		else
			t44_mix2(m, 0, sm+i, bit);
	}
	if(bpos==7)
		cm->cn=0;
	return result;
}
typedef struct T44StateStruct
{
	//Predictor::update()
	T44_APM1 a[7];

	//contextModel2()
	T44Mixer *m;

	//matchModel()
	ArrayHandle t;//<int>		hash table of pointers to contexts
	int
		h,//hash of last 7 bytes
		ptr,//points to next byte of match if any
		len,//length of match, or 0 if no match
		result;
	T44SmallStationaryContextMap scm1;

	//im24bitModel()
	T44SmallStationaryContextMap scm[10];
	T44ContextMap cm;
	T44_PRNG rnd;

	size_t bitidx;//absolute index of the next bit being predicted
	int bpos;//bits in c0 (0 to 7)
	int col;//bit idx (0 to 23)
	int c0;//partially decoded byte
	int pr;//next prediction

	int decode;
	int w3, iw, ih;
} T44State;
static void t44_state_init(T44State *state, int iw, int ih, int decode)
{
	//Predictor::update()
	t44_apm1_init(state->a+0, 256);
	t44_apm1_init(state->a+1, 0x10000);
	t44_apm1_init(state->a+2, 0x10000);
	t44_apm1_init(state->a+3, 0x10000);
	t44_apm1_init(state->a+4, 0x10000);
	t44_apm1_init(state->a+5, 0x10000);
	t44_apm1_init(state->a+6, 0x10000);
	
	//contextModel2()
	state->m=t44_mixer_alloc(845+80, 3095, 7, 0);
	
	//matchModel()
	ARRAY_ALLOC(int, state->t, 0, T44_MEM, 0, 0);
	state->h=0;
	state->ptr=0;
	state->len=0;
	state->result=0;
	t44_scm_init(&state->scm1, 0x20000);

	//im24bitModel
	const int SC=0x20000;
	t44_scm_init(state->scm+0, SC);
	t44_scm_init(state->scm+1, SC);
	t44_scm_init(state->scm+2, SC);
	t44_scm_init(state->scm+3, SC);
	t44_scm_init(state->scm+4, SC);
	t44_scm_init(state->scm+5, SC);
	t44_scm_init(state->scm+6, SC);
	t44_scm_init(state->scm+7, SC);
	t44_scm_init(state->scm+8, SC<<1);
	t44_scm_init(state->scm+9, 512);
	t44_cm_init(&state->cm, T44_MEM*4, 13);
	t44_rnd_init(&state->rnd);
	state->bitidx=0;
	state->bpos=0;
	state->col=0;
	state->c0=1;
	state->pr=2048;
	state->decode=decode;
	state->w3=3*iw;
	state->iw=iw;
	state->ih=ih;
}
static void t44_state_clear(T44State *state)
{
	for(int k=0;k<_countof(state->a);++k)
		t44_apm1_clear(state->a+k);

	t44_mixer_free(state->m);

	array_free(&state->t);
	t44_scm_clear(&state->scm1);

	for(int k=0;k<_countof(state->scm);++k)
		t44_scm_clear(state->scm+k);
	t44_cm_clear(&state->cm);
}
#define T44_MAXLEN 65534//longest allowed match + 1
#define LOADU(IDX) ((IDX)>=0?buf[IDX]&0xFF:0)
static void t44_matchModel(T44State *state, const char *buf, int c0, int bit)//finds the longest matching context and returns its length
{
	int *t=(int*)array_at(&state->t, 0);
	int tsize_m1=(int)state->t->count-1;
	int idx=(int)(state->bitidx>>3);
	if(!state->bpos)
	{
		state->h=(state->h*997*8+LOADU(idx-1)+1)&tsize_m1;//update context hash
		if(state->len)
			++state->len, ++state->ptr;
		else//find match
		{
			state->ptr=t[state->h];
			if(state->ptr && idx-state->ptr<T44_MEM)//just for compatibility with paq8pxd
			{
				for(;state->ptr-(state->len+1)>=0 && LOADU(idx-(state->len+1))==LOADU(state->ptr-(state->len+1)) && state->len<T44_MAXLEN;++state->len);
			}
		}
		t[state->h]=idx;
		state->result=state->len;
		//if(result>0 && !(result&0xFFF)) printf("pos %d  len %d  ptr %d\n", pos, len, ptr);
		t44_scm_set(&state->scm1, idx);
	}
	if(state->len)
	{
		if(LOADU(idx-1)==LOADU(state->ptr-1) && c0==(LOADU(state->ptr)+256)>>(8-state->bpos))
		{
			if(state->len>T44_MAXLEN)
				state->len=T44_MAXLEN;
			if(LOADU(state->ptr)>>(7-state->bpos)&1)
			{
				t44_mixer_add(state->m, t44_ilog(state->len)<<2);
				t44_mixer_add(state->m, MINVAR(state->len, 32)<<6);
			}
			else
			{
				t44_mixer_add(state->m, -(t44_ilog(state->len)<<2));
				t44_mixer_add(state->m, -(MINVAR(state->len, 32)<<6));
			}
		}
		else
		{
			state->len=0;
			t44_mixer_add(state->m, 0);
			t44_mixer_add(state->m, 0);
		}
	}
	else
	{
		t44_mixer_add(state->m, 0);
		t44_mixer_add(state->m, 0);
	}
	t44_scm_mix(&state->scm1, state->m, 7, c0, bit);
}
static void t44_update(T44State *state, int bit, const char *buf)
{
	//start of Predictor::update()
	state->c0<<=1;
	state->c0|=bit;
	if(state->c0>=256)
	{
		if((buf[state->bitidx>>3]&0xFF)!=(state->c0&0xFF))
			LOG_ERROR("Decode error");
		//if(state->decode)
		//	buf[state->bitidx>>3]=state->c0&0xFF;
		state->c0=1;
	}
	++state->bitidx;
	++state->bpos;
	state->bpos&=7;
	int idx=(int)(state->bitidx>>3), kc=idx%3;

	//start of contextModel2()
	t44_mixer_update(state->m, bit);
	t44_mixer_add(state->m, 256);
	t44_matchModel(state, buf, state->c0, bit);

	//start of im24bitModel()
	if(!state->bpos)
	{
#define LOAD(C, X, Y) LOADU(idx+(Y)*state->w3+3*(X)+(C))
		unsigned char
			NNWW	=LOAD(0, -2, -2),
			NN	=LOAD(0,  0, -2),
			NNE	=LOAD(0,  1, -2),
			NNEE	=LOAD(0,  2, -2),
			NW	=LOAD(0, -1, -1),
			NWp2	=LOAD(2, -1, -1),
			N	=LOAD(0,  0, -1),
		//	Np1	=LOAD(1,  0, -1),
			Np2	=LOAD(2,  0, -1),
			NE	=LOAD(0,  1, -1),
			WW	=LOAD(0, -2,  0),
			WWp2	=LOAD(2, -2,  0),
			W	=LOAD(0, -1,  0),
			Wp1	=LOAD(1, -1,  0),
			Wp2	=LOAD(2, -1,  0);
#undef  LOAD
		int mean=W+NE+N+NW;
		int var=(SQ(W)+SQ(NE)+SQ(N)+SQ(NW)-SQ(mean)/4)>>2;
		mean>>=2;
		int logvar=t44_ilog(var), i=kc%3<<4;
		t44_cm_set(&state->cm, t44_hash(++i, W, -1, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, W, Wp2, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, W, Wp2, Wp1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, N, -1, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, N, Wp2, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, N, Wp2, Wp1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, (W+N)/8, Wp2/16, Wp1/16, -1));
		t44_cm_set(&state->cm, t44_hash(++i, Wp2, Wp1, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, W, Wp2-WWp2, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, W+Wp2-WWp2, -1, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, N, Wp2-NWp2, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, N+Wp2-NWp2, -1, -1, -1));
		t44_cm_set(&state->cm, t44_hash(++i, mean, logvar>>4, -1, -1));
		t44_scm_set(state->scm+0, W+N-NW);
		t44_scm_set(state->scm+1, W+NE-N);
		t44_scm_set(state->scm+2, 2*W-WW);
		t44_scm_set(state->scm+3, 2*N-NN);
		t44_scm_set(state->scm+4, 2*NW-NNWW);
		t44_scm_set(state->scm+5, 2*NE-NNEE);
		t44_scm_set(state->scm+6, NE+Wp2-Np2);
		t44_scm_set(state->scm+7, N+NE-NNE);
		t44_scm_set(state->scm+8, (logvar<<1&0x180)|mean>>1);
	}
	//Predict next bit
	if(++state->col>=24)
		state->col=0;
	t44_scm_mix(state->scm+0, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+1, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+2, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+3, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+4, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+5, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+6, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+7, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+8, state->m, 7, state->c0, bit);
	t44_scm_mix(state->scm+9, state->m, 7, state->c0, bit);
	t44_cm_mix(&state->cm, &state->rnd, state->m, state->bpos, LOADU(idx-1), state->c0, bit, idx);
	t44_mixer_set(state->m, 2, 8);
	t44_mixer_set(state->m, state->col, 24);
	t44_mixer_set(state->m, 3*(LOADU(idx-1)>>4)+kc, 48);
	t44_mixer_set(state->m, state->c0, 256);
	//end of im24bitModel()
	
	int pr0=t44_mixer_predict(state->m, bit);
	//end of contextModel2()
	
	state->pr=t44_apm1_predict(state->a+0, pr0, state->c0, 7, bit);

	int pr1=t44_apm1_predict(state->a+1, pr0, state->c0+256*LOADU(idx-1), 7, bit);
	int pr2=t44_apm1_predict(state->a+2, pr0, (state->c0^t44_hash(LOADU(idx-1), LOADU(idx-2), -1, -1, -1))&0xFFFF, 7, bit);
	int pr3=t44_apm1_predict(state->a+3, pr0, (state->c0^t44_hash(LOADU(idx-1), LOADU(idx-2), LOADU(idx-3), -1, -1))&0xFFFF, 7, bit);
	pr0=(pr0+pr1+pr2+pr3+2)>>2;

	pr1=t44_apm1_predict(state->a+4, state->pr, state->c0+256*LOADU(idx-1), 7, bit);
	pr2=t44_apm1_predict(state->a+5, state->pr, (state->c0^t44_hash(LOADU(idx-1), LOADU(idx-2), -1, -1, -1))&0xFFFF, 7, bit);
	pr3=t44_apm1_predict(state->a+6, state->pr, (state->c0^t44_hash(LOADU(idx-1), LOADU(idx-2), LOADU(idx-3), -1, -1))&0xFFFF, 7, bit);
	state->pr=(state->pr+pr1+pr2+pr3+2)>>2;

	state->pr=(state->pr+pr0+1)>>1;
	//end of Predictor::update()
}
#undef  LOADU
int t44_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	double t_start=time_sec();
	double csize_prev=0;
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T44 paq8pxd  Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int res=iw*ih;
	char *buf=(char*)malloc((size_t)res<<2);
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf, src, (size_t)res<<2);
	addbuf((unsigned char*)buf, iw, ih, 3, 4, 128);
	T44_RCT_FWD(buf, iw, ih);
	pack3_fwd(buf, res);
	
	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
	
	double csizes[24]={0};

	T44State state;
	t44_state_init(&state, iw, ih, 0);
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(short kc=0;kc<3;++kc, ++idx)
			{
				for(short kb=7;kb>=0;--kb)
				{
					int p0=0x10000-(state.pr<<4);
					p0=CLAMP(1, p0, 0xFFFF);

					int bit=buf[idx]>>kb&1;
					ac_enc_bin(&ec, p0, bit);
					
					if(loud)
						csizes[kc<<3|kb]-=log2((double)(bit?0x10000-p0:p0)*(1./0x10000));

					t44_update(&state, bit, buf);
				}
			}
		}
		if(loud&&(ky&(LOUD_UPDATE_PERIOD-1))==LOUD_UPDATE_PERIOD-1)
		{
			double csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f\r", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*LOUD_UPDATE_PERIOD*3/(csize-csize_prev));
			csize_prev=csize;
		}
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		
		double chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=7;kb>=0;--kb)
		{
			printf("B%d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc<<3|kb;
				double size=csizes[idx];
				printf(" %15.6lf", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6lf %15.6lf %15.6lf %15.6lf\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6lf\n", (int)list.nobj, iw*ih*3./list.nobj);
	}
	t44_state_clear(&state);
	dlist_clear(&list);
	free(buf);
	return 1;
}
int t44_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	double t_start=time_sec();
	int res=iw*ih;
	memset(buf, 0, (size_t)res<<2);
	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);
	
	T44State state;
	t44_state_init(&state, iw, ih, 1);
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc, ++idx)
			{
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					int p0=0x10000-(state.pr<<4);
					p0=CLAMP(1, p0, 0xFFFF);
					
					int bit=ac_dec_bin(&ec, p0);
					buf[idx]|=bit<<kb;
					
					t44_update(&state, bit, (char*)buf);
				}
			}
		}
		if(loud&&(ky&63)==63)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
	}
	pack3_inv((char*)buf, res);
	T44_RCT_INV((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	t44_state_clear(&state);
	return 1;
}