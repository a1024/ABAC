#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define DEBUG_AVX2


#define BLOCKW 32
#define BLOCKH 32


#define V16_ZERO(DST)\
	DST[0]=_mm256_setzero_ps(),\
	DST[1]=_mm256_setzero_ps()

#define V16_LOAD(DST, PTR, IDX)\
	DST[0]=_mm256_load_ps(PTR+((IDX)<<4LL  )),\
	DST[1]=_mm256_load_ps(PTR+((IDX)<<4LL|8))

#define V16_STORE(PTR, IDX, SRC)\
	_mm256_store_ps(PTR+((IDX)<<4LL  ), SRC[0]),\
	_mm256_store_ps(PTR+((IDX)<<4LL|8), SRC[1])

#define V16_ASSIGN(DST, X)\
	DST[0]=X[0],\
	DST[1]=X[1]

#define V16_LT_C(DST, A, C)\
	DST[0]=_mm256_cmp_ps(A[0], C, _CMP_LT_OQ),\
	DST[1]=_mm256_cmp_ps(A[1], C, _CMP_LT_OQ)

#define V16_GT_C(DST, A, C)\
	DST[0]=_mm256_cmp_ps(C, A[0], _CMP_LT_OQ),\
	DST[1]=_mm256_cmp_ps(C, A[1], _CMP_LT_OQ)

#define V16_AND(DST, A, B)\
	DST[0]=_mm256_and_ps(A[0], B[0]),\
	DST[1]=_mm256_and_ps(A[1], B[1])

#define V16_AND_C(DST, A, C)\
	DST[0]=_mm256_and_ps(A[0], C),\
	DST[1]=_mm256_and_ps(A[1], C)

#define V16_ANDNOT_C(DST, A, C)\
	DST[0]=_mm256_andnot_ps(A[0], C),\
	DST[1]=_mm256_andnot_ps(A[1], C)

#define V16_OR(DST, A, B)\
	DST[0]=_mm256_or_ps(A[0], B[0]),\
	DST[1]=_mm256_or_ps(A[1], B[1])

#define V16_ADD(DST, A, B)\
	DST[0]=_mm256_add_ps(A[0], B[0]),\
	DST[1]=_mm256_add_ps(A[1], B[1])

#define V16_SUB(DST, A, B)\
	DST[0]=_mm256_sub_ps(A[0], B[0]),\
	DST[1]=_mm256_sub_ps(A[1], B[1])

#define V16_MUL(DST, VEC1, VEC2)\
	DST[0]=_mm256_mul_ps(VEC1[0], VEC2[0]),\
	DST[1]=_mm256_mul_ps(VEC1[1], VEC2[1])

#define V16_MUL_C(DST, VEC, C)\
	DST[0]=_mm256_mul_ps(VEC[0], C),\
	DST[1]=_mm256_mul_ps(VEC[1], C)

#define V16_MIN(DST, A, B)\
	DST[0]=_mm256_min_ps(A[0], B[0]),\
	DST[1]=_mm256_min_ps(A[1], B[1])

#define V16_MAX(DST, A, B)\
	DST[0]=_mm256_max_ps(A[0], B[0]),\
	DST[1]=_mm256_max_ps(A[1], B[1])

#define V16_MIN_C(DST, A, C)\
	DST[0]=_mm256_min_ps(A[0], C),\
	DST[1]=_mm256_min_ps(A[1], C)

#define V16_MAX_C(DST, A, C)\
	DST[0]=_mm256_max_ps(A[0], C),\
	DST[1]=_mm256_max_ps(A[1], C)
	
#define V16_ADD_PPN(DST, POS1, POS2, NEG)\
	V16_ADD(DST, POS1, POS2),\
	V16_SUB(DST, DST, NEG)

#define V16_AVERAGE(DST, A, B)\
	V16_ADD(DST, A, B),\
	V16_MUL_C(DST, DST, mhalf)

#define V16_PARABOLIC(DST, XP3, XN3, XP1)\
	V16_SUB(DST, XP3, XN3),\
	V16_MUL_C(DST, DST, m3),\
	V16_ADD(DST, DST, XP1)

#if 1
static int assert_ptr(const float *p)
{
	size_t val=(size_t)p;
	if(val&31)
		LOG_ERROR("Invalid alignment: %p", p);
	return 0;
}
#else
#define assert_ptr(...)
#endif

#if 0
//#define LOAD(BUF, X, Y) ((unsigned)(X)<iw&&(unsigned)(Y)<ih?BUF[iw*(Y)+(X)]:0)
#define LOAD_UY(BUF, XLO, XHI, X, Y) ((unsigned)(X-XLO)<XHI&&(unsigned)(X)<iw?BUF[iw*(Y)+(X)]:0)		//X
#define LOADVEC8(BUF, X, Y)\
	((unsigned)(Y)<ih?_mm256_set_ps(\
		LOAD_UY(BUF, X, Y),\
		LOAD_UY(BUF, X+BLOCKW, Y),\
		LOAD_UY(BUF, X+BLOCKW*2, Y),\
		LOAD_UY(BUF, X+BLOCKW*3, Y),\
		LOAD_UY(BUF, X+BLOCKW*4, Y),\
		LOAD_UY(BUF, X+BLOCKW*5, Y),\
		LOAD_UY(BUF, X+BLOCKW*6, Y),\
		LOAD_UY(BUF, X+BLOCKW*7, Y))\
	:_mm256_setzero_ps())
#define DECLVEC16(VEC, BUF, X, Y) VEC[2]={LOADVEC8(BUF, X, Y), LOADVEC8(BUF, X+BLOCKW*8, Y)}
#define LOADVEC16(VEC, BUF, X, Y) (VEC)[0]=LOADVEC8(BUF, X, Y), (VEC)[1]=LOADVEC8(BUF, X+BLOCKW*8, Y)
#endif

#ifdef DEBUG_AVX2
const unsigned char *gg_buf;
int g_iw=0, g_ih=0, g_kx0=0, g_kx=0, g_ky;
static float LOAD(int OFFSET)
{
	if(((g_kx0+OFFSET)&~(BLOCKW-1))==((g_kx+OFFSET)&~(BLOCKW-1))&&(unsigned)(g_kx+OFFSET)<(unsigned)g_iw)
	{
		const unsigned char *row=gg_buf+g_iw*g_ky+g_kx;
		if((unsigned)g_ky>=(unsigned)g_ih||(unsigned)(g_kx+OFFSET)>=(unsigned)g_iw)
			LOG_ERROR("ACCESS VIOLATION");
		return (float)row[g_kx+OFFSET];
	}
	return 0;
}
#else
#define LOAD(OFFSET) ((kx0+OFFSET)&~(BLOCKW-1))==((kx+OFFSET)&~(BLOCKW-1))&&(unsigned)(kx+OFFSET)<(unsigned)iw?(float)row[kx+OFFSET]:0
#endif
static void load_vec16(__m256 *dst, const unsigned char *buf, int iw, int ih, int kx, int ky, int ox, int oy)
{
	kx+=ox;
	ky+=oy;
	if((unsigned)ky<(unsigned)ih)
	{
		const unsigned char *row=buf+iw*ky+kx;
		int kx0=kx-ox;
#ifdef DEBUG_AVX2
		gg_buf=buf;
		g_iw=iw;
		g_ih=ih;
		g_kx0=kx0;
		g_kx=kx;
		g_ky=ky;
#endif
		dst[0]=_mm256_set_ps(
			LOAD(0),
			LOAD(BLOCKW  ),
			LOAD(BLOCKW*2),
			LOAD(BLOCKW*3),
			LOAD(BLOCKW*4),
			LOAD(BLOCKW*5),
			LOAD(BLOCKW*6),
			LOAD(BLOCKW*7)
		);
		dst[1]=_mm256_set_ps(
			LOAD(BLOCKW* 8),
			LOAD(BLOCKW* 9),
			LOAD(BLOCKW*10),
			LOAD(BLOCKW*11),
			LOAD(BLOCKW*12),
			LOAD(BLOCKW*13),
			LOAD(BLOCKW*14),
			LOAD(BLOCKW*15)
		);
#undef LOAD
	}
	else
	{
		dst[0]=_mm256_setzero_ps();
		dst[1]=_mm256_setzero_ps();
	}
}
static void store_vec16(__m256 const *src, unsigned char *buf, int iw, int ih, int kx, int ky)
{
	unsigned char *row=buf+iw*ky+kx;
	__declspec(align(32)) int vec[16];
	__m256i isrc[2]=
	{
		_mm256_cvtps_epi32(src[0]),
		_mm256_cvtps_epi32(src[1]),
	};
	_mm256_store_si256((__m256i*)vec  , isrc[0]);
	_mm256_store_si256((__m256i*)vec+1, isrc[1]);
	if((unsigned)ky>=(unsigned)ih)
		return;
	for(int kx2=0;kx2<16;++kx2)
		if((unsigned)(kx+BLOCKW*kx2)<(unsigned)iw)
			row[BLOCKW*kx2]=vec[kx2];
}

static void v16_sgn(float *dst, const float *v16)
{
	__m256 const one=_mm256_castsi256_ps(_mm256_set1_epi32(1));
	__m256 vec[2], neg[2], pos[2];

	V16_LOAD(vec, v16, 0);

	neg[0]=_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_setzero_si256(), _mm256_castps_si256(vec[0])));
	neg[1]=_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_setzero_si256(), _mm256_castps_si256(vec[1])));
	//neg[0]=_mm256_castsi256_ps(_mm256_srai_epi32(_mm256_castps_si256(vec[0]), 31));
	//neg[1]=_mm256_castsi256_ps(_mm256_srai_epi32(_mm256_castps_si256(vec[1]), 31));

	pos[0]=_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_castps_si256(vec[0]), _mm256_setzero_si256()));
	pos[1]=_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_castps_si256(vec[1]), _mm256_setzero_si256()));
	V16_AND_C(pos, pos, one);

	V16_OR(neg, neg, pos);

	neg[0]=_mm256_cvtepi32_ps(_mm256_castps_si256(neg[0]));
	neg[1]=_mm256_cvtepi32_ps(_mm256_castps_si256(neg[1]));

	V16_STORE(dst, 0, neg);
}

//count is number of v16 pseudo-registers
static void act(float *dst, const float *src, int count)
{
	__m256 const
		onep=_mm256_set1_ps(0.01f),
		vmin=_mm256_set1_ps(-10),
		vmax=_mm256_set1_ps(10);
	assert_ptr(src);
	assert_ptr(dst);
	for(int k=0;k<count;++k)
	{
		__m256 val[2], negpart[2];

		V16_LOAD(val, src, k);
		V16_MUL_C(negpart, val, onep);
		V16_MAX(val, val, negpart);
		V16_MIN_C(val, val, vmax);
		V16_MAX_C(val, val, vmin);
		V16_STORE(dst, k, val);

		//DataType val=src[k];
		//DataType negpart=MUL(val, ONE_PERCENT);
		//val=val>negpart?val:negpart;
		//val=CLAMP(-MAXMAG, val, MAXMAG);
		//dst[k]=val;
	}
}
static void act_dash(float *dst, const float *src, int count)
{
	__m256 const
		zero=_mm256_setzero_ps(),
		onep=_mm256_set1_ps(0.01f),
		vmin=_mm256_set1_ps(-10),
		vmax=_mm256_set1_ps(10);
	assert_ptr(src);
	assert_ptr(dst);
	for(int k=0;k<count;++k)
	{
		__m256 val[2], cmp[2];

		V16_LOAD(val, src, k);

		V16_GT_C(cmp, val, vmin);
		V16_AND(val, val, cmp);

		V16_LT_C(cmp, val, vmax);
		V16_AND(val, val, cmp);

		V16_LT_C(cmp, val, zero);
		V16_AND(val, val, cmp);
		V16_ANDNOT_C(cmp, cmp, onep);
		V16_OR(val, val, cmp);
		V16_STORE(dst, k, val);

		//DataType val=src[k];
		//if(val<-MAXMAG||val>MAXMAG)
		//	val=0;
		//else
		//	val=val<0?ONE_PERCENT:ONE;
		//dst[k]=val;
	}
}
static void vec_scalar(float *dst, const float *vec, const float *scalar, int count)
{
	__m256 ms[2];

	assert_ptr(vec);
	assert_ptr(scalar);
	assert_ptr(dst);

	V16_LOAD(ms, scalar, 0);
	for(int k=0;k<count;++k)
	{
		__m256 val[2];

		V16_LOAD(val, vec, k);
		V16_MUL(val, val, ms);
		V16_STORE(dst, k, val);
	}
	//	dst[k]=MUL(vec[k], scalar);
}
static void vec_ew(float *dst, const float *v1, const float *v2, int count)
{
	assert_ptr(v1);
	assert_ptr(v2);
	assert_ptr(dst);
	for(int k=0;k<count;++k)
	{
		__m256 val[2], val2[2];
		V16_LOAD(val, v1, k);
		V16_LOAD(val2, v2, k);
		V16_MUL(val, val, val2);
		V16_STORE(dst, k, val);
	}
	//	dst[k]=MUL(v1[k], v2[k]);
}
static void vec_outer(float *dst, const float *left, const float *right, int lh, int rw)
{
	assert_ptr(left);
	assert_ptr(right);
	assert_ptr(dst);
	for(int ky=0;ky<lh;++ky)
	{
		for(int kx=0;kx<rw;++kx)
		{
			__m256 vl[2], vr[2];

			V16_LOAD(vl, left, ky);
			V16_LOAD(vr, right, kx);
			V16_MUL(vl, vl, vr);
			V16_STORE(dst, rw*ky+kx, vl);
		}
		//	dst[rw*ky+kx]=MUL(left[ky], right[kx]);
	}
}
static void matmul(float *dst, const float *m1, const float *m2, int h1, int w1h2, int w2)
{
	assert_ptr(m1);
	assert_ptr(m2);
	assert_ptr(dst);
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			__m256 sum[2];

			V16_ZERO(sum);
			for(int k=0;k<w1h2;++k)
			{
				__m256 mleft[2], mright[2];

				V16_LOAD(mleft, m1, w1h2*ky+k);
				V16_LOAD(mright, m2, w2*k+kx);
				V16_MUL(mleft, mleft, mright);
				V16_ADD(sum, sum, mleft);
			}
			//	sum+=MUL(m1[w1h2*ky+k], m2[w2*k+kx]);
			V16_STORE(dst, w2*ky+kx, sum);
			//dst[w2*ky+kx]=sum;
		}
	}
}
static void linear(float *dst, const float *mat, const float *vec, const float *bias, int win, int hout)
{
	assert_ptr(mat);
	assert_ptr(vec);
	assert_ptr(bias);
	assert_ptr(dst);
	for(int ko=0;ko<hout;++ko)
	{
		__m256 sum[2];

		if(bias)
			V16_LOAD(sum, bias, 0);
		else
			V16_ZERO(sum);
		//DataType temp=bias?bias[ko]:0;
		for(int ki=0;ki<win;++ki)
		{
			__m256 mleft[2], mright[2];

			V16_LOAD(mleft, mat, win*ko+ki);
			V16_LOAD(mright, vec, ki);
			V16_MUL(mleft, mleft, mright);
			V16_ADD(sum, sum, mleft);
		}
		//	temp+=MUL(mat[win*ko+ki], vec[ki]);
		V16_STORE(dst, ko, sum);
		//dst[ko]=temp;
	}
}


#define NF0 94
#define NF1 32
#define NF2 32
#define NF3 32
typedef struct TempsStruct
{
	__m256
		nb[NF0][2],
		net1[NF1][2], x1[NF1][2],
		net2[NF2][2], x2[NF2][2],
		net3[NF3][2], x3[NF3][2],
		pred[2], expr[2], loss[2];
} Temps;
typedef struct BwdTempsStruct
{
	__m256 dL_dp[2], dL_dn3[NF3][2], dL_dn2[NF2][2], actdash_n2[NF1][2], dL_dn1[NF1][2], actdash_n1[NF1][2];
} BwdTemps;
typedef struct ParamsStruct
{
	__m256
		weight1[NF1*NF0][2], bias1[NF1][2],
		weight2[NF2*NF1][2], bias2[NF2][2],
		weight3[NF3*NF2][2], bias3[NF3][2],
		weight4[1*NF3][2];
} Params;

static void initialize(float *w, int count, float sqrt_fan_in)
{
	for(int k=0;k<count;++k)
	{
		int x=(int)(xoroshiro128_next()&0x1FFFF)-0x10000;
		__m256 wk=_mm256_set1_ps((float)x/(0x10000*sqrt_fan_in));
		_mm256_store_ps(w+(k<<4  ), wk);
		_mm256_store_ps(w+(k<<4|8), wk);
	}
}
static void initialize_all(Params *p)
{
	XOROSHIRO128_RESET();
	initialize((float*)p->weight1	, _countof(p->weight1)	, sqrtf(NF0));
	initialize((float*)p->bias1		, _countof(p->bias1)	, sqrtf(NF0));
	initialize((float*)p->weight2	, _countof(p->weight2)	, sqrtf(NF1));
	initialize((float*)p->bias2		, _countof(p->bias2)	, sqrtf(NF1));
	initialize((float*)p->weight3	, _countof(p->weight3)	, sqrtf(NF2));
	initialize((float*)p->bias3		, _countof(p->bias3)	, sqrtf(NF2));
	initialize((float*)p->weight4	, _countof(p->weight4)	, sqrtf(NF3));
}

static void clamp4(__m256 *dst, __m256 const *p, __m256 const *a, __m256 const *b, __m256 const *c, __m256 const *d)
{
	__m256 vmin[2], vmax[2];

	vmin[0]=_mm256_min_ps(a[0], b[0]);
	vmin[1]=_mm256_min_ps(a[1], b[1]);
	vmax[0]=_mm256_max_ps(a[0], b[0]);
	vmax[1]=_mm256_max_ps(a[1], b[1]);
	vmin[0]=_mm256_min_ps(vmin[0], c[0]);
	vmin[1]=_mm256_min_ps(vmin[1], c[1]);
	vmax[0]=_mm256_max_ps(vmax[0], c[0]);
	vmax[1]=_mm256_max_ps(vmax[1], c[1]);
	vmin[0]=_mm256_min_ps(vmin[0], d[0]);
	vmin[1]=_mm256_min_ps(vmin[1], d[1]);
	vmax[0]=_mm256_max_ps(vmax[0], d[0]);
	vmax[1]=_mm256_max_ps(vmax[1], d[1]);

	dst[0]=_mm256_max_ps(p[0], vmin[0]);
	dst[1]=_mm256_max_ps(p[1], vmin[1]);
	dst[0]=_mm256_min_ps(dst[0], vmax[0]);
	dst[1]=_mm256_min_ps(dst[1], vmax[1]);
}
static void clip(__m256 *dst, __m256 const *x)
{
	__m256 const one=_mm256_set1_ps(1), mone=_mm256_set1_ps(-1);

	dst[0]=_mm256_max_ps(x[0], mone);
	dst[1]=_mm256_max_ps(x[1], mone);
	dst[0]=_mm256_min_ps(dst[0], one);
	dst[1]=_mm256_min_ps(dst[1], one);
}
static void get_nb2(const char *buf, const char *errors, int iw, int ih, int kx, int ky, __m256 *ctx)
{
	const int CW=7;
	int idx=iw*ky+kx, w3=iw*3, w2=iw*2;
	__m256
		NNNNNN[2], NNNNN[2],

		NNNNWWW	[2],
		NNNN	[2],
		NNNNEEE	[2],
		
		NNNWWWWW[2],
		NNNWWWW	[2],
		NNNWW	[2],
		NNNW	[2],
		NNN		[2],
		NNNE	[2],
		NNNEE	[2],
		NNNEEE	[2],
		NNNEEEE	[2],

		NNWWW	[2],
		NNWW	[2],
		NNW		[2],
		NN		[2],
		NNE		[2],
		NNEE	[2],
		NNEEE	[2],
		NNEEEE	[2],

		NWWW	[2],
		NWW		[2],
		NW		[2],
		N		[2],
		NE		[2],
		NEE		[2],
		NEEE	[2],
		NEEEE	[2],
		NEEEEE	[2],
		NEEEEEE	[2],
		NEEEEEEE[2],

		WWWWWW	[2],
		WWWWW	[2],
		WWWW	[2],
		WWW		[2],
		WW		[2],
		W		[2];

	load_vec16(NNNNNN	, buf, iw, ih, kx, ky, 0, -6);
	load_vec16(NNNNN	, buf, iw, ih, kx, ky, 0, -5);

	load_vec16(NNNNWWW	, buf, iw, ih, kx, ky, -3, -4);
	load_vec16(NNNN		, buf, iw, ih, kx, ky,  0, -4);
	load_vec16(NNNNEEE	, buf, iw, ih, kx, ky,  3, -4);

	load_vec16(NNNWWWWW	, buf, iw, ih, kx, ky, -5, -3);
	load_vec16(NNNWWWW	, buf, iw, ih, kx, ky, -4, -3);
	load_vec16(NNNWW	, buf, iw, ih, kx, ky, -2, -3);
	load_vec16(NNNW		, buf, iw, ih, kx, ky, -1, -3);
	load_vec16(NNN		, buf, iw, ih, kx, ky, 	0, -3);
	load_vec16(NNNE		, buf, iw, ih, kx, ky,  1, -3);
	load_vec16(NNNEE	, buf, iw, ih, kx, ky,  2, -3);
	load_vec16(NNNEEE	, buf, iw, ih, kx, ky,  3, -3);
	load_vec16(NNNEEEE	, buf, iw, ih, kx, ky,  4, -3);

	load_vec16(NNWWW	, buf, iw, ih, kx, ky, -3, -2);
	load_vec16(NNWW		, buf, iw, ih, kx, ky, -2, -2);
	load_vec16(NNW		, buf, iw, ih, kx, ky, -1, -2);
	load_vec16(NN		, buf, iw, ih, kx, ky,  0, -2);
	load_vec16(NNE		, buf, iw, ih, kx, ky,  1, -2);
	load_vec16(NNEE		, buf, iw, ih, kx, ky,  2, -2);
	load_vec16(NNEEE	, buf, iw, ih, kx, ky,  3, -2);
	load_vec16(NNEEEE	, buf, iw, ih, kx, ky,  4, -2);

	load_vec16(NWWW		, buf, iw, ih, kx, ky, -3, -1);
	load_vec16(NWW		, buf, iw, ih, kx, ky, -2, -1);
	load_vec16(NW		, buf, iw, ih, kx, ky, -1, -1);
	load_vec16(N		, buf, iw, ih, kx, ky, 	0, -1);
	load_vec16(NE		, buf, iw, ih, kx, ky,  1, -1);
	load_vec16(NEE		, buf, iw, ih, kx, ky,  2, -1);
	load_vec16(NEEE		, buf, iw, ih, kx, ky,  3, -1);
	load_vec16(NEEEE	, buf, iw, ih, kx, ky,  4, -1);
	load_vec16(NEEEEE	, buf, iw, ih, kx, ky,  5, -1);
	load_vec16(NEEEEEE	, buf, iw, ih, kx, ky,  6, -1);
	load_vec16(NEEEEEEE	, buf, iw, ih, kx, ky,  7, -1);
	
	load_vec16(WWWWWW	, buf, iw, ih, kx, ky, -6, 0);
	load_vec16(WWWWW	, buf, iw, ih, kx, ky, -5, 0);
	load_vec16(WWWW		, buf, iw, ih, kx, ky, -4, 0);
	load_vec16(WWW		, buf, iw, ih, kx, ky, -3, 0);
	load_vec16(WW		, buf, iw, ih, kx, ky, -2, 0);
	load_vec16(W		, buf, iw, ih, kx, ky, -1, 0);
#if 0
	__m256
		DECLVEC16(NNNNNN	, buf, kx	, ky-6),

		DECLVEC16(NNNNN		, buf, kx	, ky-5),

		DECLVEC16(NNNNWWW	, buf, kx-3	, ky-4),
		DECLVEC16(NNNN		, buf, kx	, ky-4),
		DECLVEC16(NNNNEEE	, buf, kx+3	, ky-4),

		DECLVEC16(NNNWWWWW	, buf, kx-5	, ky-3),
		DECLVEC16(NNNWWWW	, buf, kx-4	, ky-3),
		DECLVEC16(NNNWW		, buf, kx-2	, ky-3),
		DECLVEC16(NNNW		, buf, kx-1	, ky-3),
		DECLVEC16(NNN		, buf, kx	, ky-3),
		DECLVEC16(NNNE		, buf, kx+1	, ky-3),
		DECLVEC16(NNNEE		, buf, kx+2	, ky-3),
		DECLVEC16(NNNEEE	, buf, kx+3	, ky-3),
		DECLVEC16(NNNEEEE	, buf, kx+4	, ky-3),

		DECLVEC16(NNWWW		, buf, kx-3	, ky-2),
		DECLVEC16(NNWW		, buf, kx-2	, ky-2),
		DECLVEC16(NNW		, buf, kx-1	, ky-2),
		DECLVEC16(NN		, buf, kx	, ky-2),
		DECLVEC16(NNE		, buf, kx+1	, ky-2),
		DECLVEC16(NNEE		, buf, kx+2	, ky-2),
		DECLVEC16(NNEEE		, buf, kx+3	, ky-2),
		DECLVEC16(NNEEEE	, buf, kx+4	, ky-2),

		DECLVEC16(NWWW		, buf, kx-3	, ky-1),
		DECLVEC16(NWW		, buf, kx-2	, ky-1),
		DECLVEC16(NW		, buf, kx-1	, ky-1),
		DECLVEC16(N			, buf, kx	, ky-1),
		DECLVEC16(NE		, buf, kx+1	, ky-1),
		DECLVEC16(NEE		, buf, kx+2	, ky-1),
		DECLVEC16(NEEE		, buf, kx+3	, ky-1),
		DECLVEC16(NEEEE		, buf, kx+4	, ky-1),
		DECLVEC16(NEEEEE	, buf, kx+5	, ky-1),
		DECLVEC16(NEEEEEE	, buf, kx+6	, ky-1),
		DECLVEC16(NEEEEEEE	, buf, kx+7	, ky-1),
		
		DECLVEC16(WWWWWW	, buf, kx-6	, ky  ),
		DECLVEC16(WWWWW		, buf, kx-5	, ky  ),
		DECLVEC16(WWWW		, buf, kx-4	, ky  ),
		DECLVEC16(WWW		, buf, kx-3	, ky  ),
		DECLVEC16(WW		, buf, kx-2	, ky  ),
		DECLVEC16(W			, buf, kx-1	, ky  );
#endif
	__m256 const
		m_10=_mm256_set1_ps(0.1f), m_6=_mm256_set1_ps(1.f/6), m_5=_mm256_set1_ps(0.2f), m_4=_mm256_set1_ps(0.25f), m_3=_mm256_set1_ps(1.f/3), mhalf=_mm256_set1_ps(0.5f),
		m3=_mm256_set1_ps(3), m4=_mm256_set1_ps(4), m5=_mm256_set1_ps(5), m6=_mm256_set1_ps(6), m8=_mm256_set1_ps(8), m10=_mm256_set1_ps(10), m15=_mm256_set1_ps(15), m20=_mm256_set1_ps(20);
	__m256 temp[2], temp2[2];


	//0
	V16_ADD_PPN(temp, W, N, NW);
	clamp4(ctx, temp, W, NW, N, NE);
	ctx+=2;

	//1
	clip(ctx, temp);
	ctx+=2;

	//2
	V16_ADD_PPN(temp, W, NE, N);
	clamp4(ctx, temp, W, NW, N, NE);
	ctx+=2;

	//3
	clip(ctx, temp);
	ctx+=2;

	//4
	V16_ADD_PPN(temp, W, NW, NNW);
	clamp4(ctx, temp, W, NW, N, NE);
	ctx+=2;

	//5
	clip(ctx, temp);
	ctx+=2;

	//6
	V16_ADD_PPN(temp, N, NE, NNE);
	clamp4(ctx, temp, W, NW, N, NE);
	ctx+=2;

	//7
	clip(ctx, temp);
	ctx+=2;

	//8
	V16_AVERAGE(temp, W, NE);
	V16_ASSIGN(ctx, temp);
	ctx+=2;

	//9
	V16_PARABOLIC(temp, N, NN, NNN);
	clip(ctx, temp);
	ctx+=2;

	//10
	V16_PARABOLIC(temp, W, WW, WWW);
	clip(ctx, temp);
	ctx+=2;

	//11
	V16_PARABOLIC(temp, NE, NNE, NNNE);
	clip(temp, temp);
	V16_ADD(ctx, temp, W);
	ctx+=2;

	//12
	V16_PARABOLIC(temp, NEE, NNEEE, NNNEEEE);
	clip(temp, temp);
	V16_ADD(ctx, temp, W);
	ctx+=2;

	//13
	V16_ADD_PPN(temp, NN, NNNN, NNNNNN);
	clip(ctx, temp);
	ctx+=2;

	//14
	V16_ADD_PPN(temp, WW, WWWW, WWWWWW);
	clip(ctx, temp);
	ctx+=2;

	//15
	V16_ADD(temp, W, W);
	//V16_MUL_C(temp, W, m2);
	V16_SUB(temp, temp, NWW);
	clamp4(temp, temp, W, NW, N, NN);
	V16_ADD(temp, temp, NNNNN);
	V16_MUL_C(temp2, NNNN, m6), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, NNN, m15), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, NN, m20), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, N, m15), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp, temp, m_6);
	clip(ctx, temp);
	ctx+=2;

	//16
	V16_PARABOLIC(temp, NEE, NNEE, NNNEE);
	clamp4(temp, temp, NE, NEE, NEEE, NEEEE);
	V16_MUL_C(temp2, W, m8), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, WW, m3), V16_SUB(temp, temp, temp2);
	clip(ctx, temp);
	ctx+=2;

	//17
	V16_ADD_PPN(temp, NN, NW, NNNW);
	clip(ctx, temp);
	ctx+=2;

	//18
	V16_ADD_PPN(temp, NN, NE, NNNE);
	clip(ctx, temp);
	ctx+=2;

	//19
	V16_SUB(temp, W, WW);
	V16_ADD(temp, temp, temp);
	//V16_MUL_C(temp, temp, m2);
	V16_SUB(temp2, NW, NWW);
	V16_ADD(temp, temp, temp2);
	V16_ADD(temp, temp, WWW);
	clip(ctx, temp);
	ctx+=2;

	//20
	V16_AVERAGE(temp, NW, NWW);
	V16_AVERAGE(temp2, NNNWWWW, NNNWWWWW);
	V16_PARABOLIC(temp, temp, NNWWW, temp2);
	clip(ctx, temp);
	ctx+=2;

	//21
	V16_ADD_PPN(temp, NEE, NE, NNEEE);
	clip(ctx, temp);
	ctx+=2;

	//22
	V16_ADD_PPN(temp, NWW, WW, WWWW);
	clip(ctx, temp);
	ctx+=2;

	//23
	V16_AVERAGE(temp, W, NW);
	V16_AVERAGE(temp2, NWWW, NNWWW);
	V16_PARABOLIC(temp, temp, NWW, temp2);
	clip(ctx, temp);
	ctx+=2;

	//24
	V16_SUB(temp, NE, NNNEE);
	V16_ADD(temp, temp, temp);
	//V16_MUL_C(temp, temp, m2);
	V16_SUB(temp2, NEE, NNEE);
	V16_ADD(temp, temp, temp2);
	V16_ADD(temp, temp, NNNNEEE);
	clip(ctx, temp);
	ctx+=2;

	//25
	V16_ASSIGN(ctx, NNNNNN);
	ctx+=2;

	//26
	V16_AVERAGE(ctx, NEEEE, NEEEEEE);
	ctx+=2;

	//27
	V16_AVERAGE(ctx, WWWW, WWWWWW);
	ctx+=2;

	//28
	V16_ADD(temp, W, N);
	V16_ADD(temp, temp, NEEEEE);
	V16_ADD(temp, temp, NEEEEEEE);
	V16_MUL_C(ctx, temp, m_4);
	ctx+=2;

	//29
	V16_ADD_PPN(temp, NEEE, W, NEE);
	clip(ctx, temp);
	ctx+=2;

	//30
	V16_MUL_C(temp, NNN, m4);
	V16_MUL_C(temp2, NNNN, m3);
	V16_SUB(temp, temp, temp2);
	clip(ctx, temp);
	ctx+=2;

	//31
	V16_ADD_PPN(temp, N, NN, NNN);
	clip(ctx, temp);
	ctx+=2;

	//32
	V16_ADD_PPN(temp, W, WW, WWW);
	clip(ctx, temp);
	ctx+=2;

	//33
	V16_ADD_PPN(temp, W, NEE, NE);
	clip(ctx, temp);
	ctx+=2;

	//34
	V16_ADD_PPN(temp, WW, NEE, N);
	clip(ctx, temp);
	ctx+=2;

	//35
	V16_ADD(temp2, W, W);
	V16_SUB(temp, temp2, NW);
	clip(temp, temp);
	V16_SUB(temp2, temp2, NWW);
	clip(temp2, temp2);
	V16_ADD(temp, temp, temp2);
	V16_ADD(temp, temp, N);
	V16_ADD(temp, temp, NE);
	V16_MUL_C(temp, temp, m_4);
	ctx+=2;

	//36
	V16_ADD_PPN(temp, N, N, NN);
	clamp4(ctx, temp, W, N, NE, NEE);
	ctx+=2;

	//37
	V16_AVERAGE(ctx, N, NNN);
	ctx+=2;

	//38
	V16_ADD_PPN(temp, NN, W, NNW);
	clip(ctx, temp);
	ctx+=2;

	//39
	V16_ADD_PPN(temp, NEE, N, NNEE);
	clip(ctx, temp);
	ctx+=2;

	//40
	V16_ADD_PPN(temp, NEE, NEE, NNEE);
	clip(temp, temp);
	V16_MUL_C(temp2, W, m20), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, WW, m15), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, WWW, m4), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp, temp, m_10);
	clip(ctx, temp);
	ctx+=2;

	//41
	V16_PARABOLIC(temp, W, NW, NNW);
	clip(temp, temp);
	V16_MUL_C(temp2, NE, m6), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, NNEE, m4), V16_SUB(temp, temp, temp2);
	V16_ADD(temp, temp, NNNEEE);
	V16_MUL_C(temp, temp, m_4);
	clip(ctx, temp);
	ctx+=2;

	//42
	V16_SUB(temp, N, NNE);
	V16_ADD(temp, temp, temp);
	V16_ADD(temp, temp, NNW);
	V16_SUB(temp, temp, NN);
	V16_ADD(temp, temp, NNNE);
	clip(ctx, temp);
	ctx+=2;

	//43
	V16_SUB(temp, NW, NNNWW);
	V16_ADD(temp, temp, temp);
	V16_ADD(temp, temp, NNW);
	V16_SUB(temp, temp, NNWW);
	V16_ADD(temp, temp, NNNNWWW);
	clip(ctx, temp);
	ctx+=2;

	//44
	V16_ADD_PPN(temp, NNWW, W, NNWWW);
	clip(ctx, temp);
	ctx+=2;

	//45
	V16_MUL_C(temp, W, m4);
	V16_MUL_C(temp2, NWW, m6), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, NNWWW, m4), V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, NNNWWWW);
	clip(temp, temp);
	V16_MUL_C(temp2, N, m10), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, NN, m10), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, NNN, m5), V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, NNNN);
	V16_MUL_C(temp, temp, m_5);
	clip(ctx, temp);
	ctx+=2;

	//46
	V16_ADD_PPN(temp, NEEE, NEEE, NNEEEE);
	V16_SUB(temp, temp, NEEEE);
	clip(temp, temp);
	V16_ADD(temp, temp, NEE);
	clip(ctx, temp);
	ctx+=2;

	//47
	V16_ADD_PPN(temp, NW, W, NWW);
	clip(ctx, temp);
	ctx+=2;

	//48
	V16_SUB(temp, N, NNW);
	V16_ADD(temp, temp, temp);
	V16_ADD(temp, temp, NW);
	V16_SUB(temp, temp, NN);
	V16_ADD(temp, temp, NNNW);
	clip(ctx, temp);
	ctx+=2;

	//49
	V16_ADD_PPN(temp, NEE, NEE, NNEEE);
	clip(temp, temp);
	V16_SUB(temp, temp, NNE);
	V16_ADD(temp, temp, NN);
	clip(ctx, temp);
	ctx+=2;

	//50
	V16_ADD_PPN(temp, NE, NE, NNE);
	clip(temp, temp);
	V16_MUL_C(temp2, W, m10), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, WW, m10), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, WWW, m5), V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, WWWW);
	V16_MUL_C(temp, temp, m_5);
	clip(ctx, temp);
	ctx+=2;

	//51
	V16_ADD_PPN(temp, NE, NE, NNE);
	clip(temp, temp);
	V16_MUL_C(temp2, W, m5), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, WWW, m5), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, WWWW, m4), V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, WWWWW);
	V16_MUL_C(temp, temp, m_4);
	clip(ctx, temp);
	ctx+=2;

	//52
	V16_PARABOLIC(temp, NE, NNE, NNNE);
	clip(temp, temp);
	V16_MUL_C(temp2, W, m6), V16_ADD(temp, temp, temp2);
	V16_MUL_C(temp2, WW, m4), V16_SUB(temp, temp, temp2);
	V16_ADD(temp, temp, WWW);
	V16_MUL_C(temp, temp, m_4);
	clip(ctx, temp);
	ctx+=2;

	//53
	V16_MUL_C(temp, W, m4);
	V16_MUL_C(temp2, NW, m6), V16_SUB(temp, temp, temp2);
	V16_MUL_C(temp2, NNW, m4), V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, NNNW);
	clip(temp, temp);
	V16_MUL_C(temp2, NE, m3);
	V16_ADD(temp, temp, temp2);
	V16_SUB(temp, temp, NNEE);
	V16_MUL_C(temp, temp, m_3);
	clip(ctx, temp);
	ctx+=2;

	//54		((W+N)*3-NW*2)/4 = (W+N)/2*3 - NW/2
	V16_AVERAGE(temp, W, N);
	V16_MUL_C(temp, temp, m3);
	V16_MUL_C(temp2, NW, mhalf);
	V16_SUB(ctx, temp, temp2);
	ctx+=2;

	load_vec16(ctx, errors, iw, ih, kx, ky, 0, -6); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, 0, -5); ctx+=2;
	
	load_vec16(ctx, errors, iw, ih, kx, ky, -3, -4); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  0, -4); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  3, -4); ctx+=2;

	load_vec16(ctx, errors, iw, ih, kx, ky, -5, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -4, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -2, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -1, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  0, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  1, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  2, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  3, -3); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  4, -3); ctx+=2;

	load_vec16(ctx, errors, iw, ih, kx, ky, -3, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -2, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -1, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  0, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  1, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  2, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  3, -2); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  4, -2); ctx+=2;

	load_vec16(ctx, errors, iw, ih, kx, ky, -3, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -2, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -1, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  0, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  1, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  2, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  3, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  4, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  5, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  6, -1); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky,  7, -1); ctx+=2;
	
	load_vec16(ctx, errors, iw, ih, kx, ky, -6, 0); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -5, 0); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -4, 0); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -3, 0); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -2, 0); ctx+=2;
	load_vec16(ctx, errors, iw, ih, kx, ky, -1, 0); ctx+=2;

	
//#undef V16_ASSIGN
//#undef V16_ADD
//#undef V16_SUB
//#undef V16_MUL_C
//#undef V16_ADD_PPN
//#undef V16_AVERAGE
//#undef V16_PARABOLIC


#if 0
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
#endif
#if 0
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
#endif
}

static void eval_fwd(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params *p)
{
	get_nb2(src, errors, iw, ih, kx, ky, (__m256*)t->nb);

	linear((float*)t->net1, (float*)p->weight1, (float*)t->nb, (float*)p->bias1, NF0, NF1);	//n1 = w1*nb + b1
	act((float*)t->x1, (float*)t->net1, NF1);								//x1 = act(n1)

	linear((float*)t->net2, (float*)p->weight2, (float*)t->x1, (float*)p->bias2, NF1, NF2);	//n2 = w2*x1 + b2
	act((float*)t->x2, (float*)t->net2, NF2);								//x2 = act(n2)

	linear((float*)t->net3, (float*)p->weight3, (float*)t->x2, (float*)p->bias3, NF2, NF3);	//n3 = w3*x2
	act((float*)t->x3, (float*)t->net3, NF3);								//x3 = act(n3)

	linear((float*)t->pred, (float*)p->weight4, (float*)t->x3, 0, NF3, 1);			//pred = w4*x3
}
static void train(const char *src, const char *errors, int iw, int ih, int kx, int ky, Temps *t, Params *p, BwdTemps *b, Params *g, float lr)
{
	eval_fwd(src, errors, iw, ih, kx, ky, t, p);

	//loss
	__m256 const signmask=_mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
	__m256 curr[2];
	load_vec16(curr, src, iw, ih, kx, ky, 0, 0);
	V16_SUB(t->expr, t->pred, curr);
	V16_AND_C(t->loss, t->expr, signmask);


	//bwd
	v16_sgn((float*)b->dL_dp, (float*)t->expr);//dL_dp = sgn(p-x)

	vec_scalar((float*)g->weight4, (float*)t->x3, (float*)b->dL_dp, NF3);//dL_dw3 = dL_dp * dp_dw3		= sgn(p-x) * x2T

	//dL_dn3 = dL_dp * dp_dx3 * dx3_dn3			= sgn(p-x) * w4T * act_dash(n3)
	act_dash((float*)b->dL_dn3, (float*)t->net3, NF3);
	vec_ew((float*)b->dL_dn3, (float*)b->dL_dn3, (float*)p->weight4, NF3);
	vec_scalar((float*)b->dL_dn3, (float*)b->dL_dn3, (float*)b->dL_dp, NF3);


	memcpy(g->bias3, b->dL_dn3, sizeof(g->bias3));//dL_db3 = dL_dn3
	vec_outer((float*)g->weight3, (float*)b->dL_dn3, (float*)t->x2, NF3, NF2);//dL_dw3 = dL_dn3 * dn3_dw3		= dL_dn3 * x2T (outer product)

	//dL_dn2 = matmul(dL_dn3T, w3)T .* act_dash(n2)
	matmul((float*)b->dL_dn2, (float*)b->dL_dn3, (float*)p->weight3, 1, NF3, NF2);
	act_dash((float*)b->actdash_n2, (float*)t->net2, NF2);
	vec_ew((float*)b->dL_dn2, (float*)b->dL_dn2, (float*)b->actdash_n2, NF2);


	memcpy(g->bias2, b->dL_dn2, sizeof(g->bias2));//dL_db2 = dL_dn2
	vec_outer((float*)g->weight2, (float*)b->dL_dn2, (float*)t->x1, NF2, NF1);//dL_dw2 = dL_dn2 * dn2_dw2		= dL_dn2 * x1T (outer product)

	//dL_dn1 = matmul(dL_dn2T, w2)T .* act_dash(n1)
	matmul((float*)b->dL_dn1, (float*)b->dL_dn2, (float*)p->weight2, 1, NF2, NF1);
	act_dash((float*)b->actdash_n1, (float*)t->net1, NF1);
	vec_ew((float*)b->dL_dn1, (float*)b->dL_dn1, (float*)b->actdash_n1, NF1);
				

	memcpy(g->bias1, b->dL_dn1, sizeof(g->bias1));//dL_db1 = dL_dn1
	vec_outer((float*)g->weight1, (float*)b->dL_dn1, (float*)t->nb, NF1, NF0);//dL_dw1 = dL_dn1 * dn1_dw1


	//update
	{
		const float *gradient=(const float*)g;
		float *params=(float*)p;
		__m256 vlr=_mm256_set1_ps(lr);
		for(int k=0;k<sizeof(Params)/sizeof(__m256[2]);++k)
		{
			__m256 vp[2], vg[2];

			V16_LOAD(vp, params, k);
			V16_LOAD(vg, gradient, k);
			V16_MUL_C(vg, vg, vlr);
			V16_SUB(vp, vp, vg);
			V16_STORE(params, k, vp);
		}
	}
	//DataType *params=(DataType*)p, *gradient=(DataType*)g;
	//for(int k=0;k<sizeof(Params)/sizeof(DataType);++k)
	//	params[k]-=LEARNING_RATE(gradient[k], lr);
}
void pred_learned_cpu(char *buf, int iw, int ih, int fwd)
{
	double t_start=time_ms();
	int res=iw*ih;
	char *src=(char*)malloc(res), *dst=(char*)malloc(res);
	Temps *t=(Temps*)_mm_malloc(sizeof(Temps), 32);
	Params *p=(Params*)_mm_malloc(sizeof(Params), 32);
	BwdTemps *b=(BwdTemps*)_mm_malloc(sizeof(BwdTemps), 32);
	Params *g=(Params*)_mm_malloc(sizeof(Params), 32);
	if(!src||!dst||!t||!p||!b||!g)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	ArrayHandle title;
	STR_ALLOC(title, 1024);
	get_window_title(title->data, 1024);
	
	const int reach=1;//5

	char *buf2=fwd?src:dst, *errors=fwd?dst:src;
	float lr=0.001f, lr2;
	__m256 curr[2];
	for(int kc=0;kc<3;++kc)
	{
		initialize_all(p);

		for(int k=0;k<res;++k)
			src[k]=buf[k<<2|kc];
		
		for(int ky=0;ky<ih;++ky)
		{
			//if(ky>=(ih>>1))//
			//	printf("");

			if(!(ky&15))
			{
				TimeInfo ti;
				parsetimedelta(time_ms()-t_start, &ti);
				set_window_title("%d/3, %d/%d, %.2lf%% - %02d-%02d-%06.3f", kc+1, ky+1, ih, 100.*(ih*kc+ky+1)/(ih*3), ti.hours, ti.mins, ti.secs);
			}
			
			//for(int kx0=0;kx0<iw-(BLOCKW*16-1);kx0+=BLOCKW*16)
			for(int kx0=0;kx0<iw;kx0+=BLOCKW*16)
			{
				for(int kx=0;kx<BLOCKW;++kx)
				{
					//int niter=kc==1?1:6;
					for(int k=0;k<1;++k)
					{
						for(int ky2=-reach;ky2<0;++ky2)
						{
							for(int kx2=-reach;kx2<=reach;++kx2)
							{
								lr2=lr/sqrtf((float)(kx2*kx2+ky2*ky2));
								train(buf2, errors, iw, ih, kx0+kx+kx2, ky+ky2, t, p, b, g, lr2);
							}
						}
						for(int kx2=-reach;kx2<0;++kx2)
						{
							lr2=lr/abs(kx2);
							train(buf2, errors, iw, ih, kx0+kx+kx2, ky, t, p, b, g, lr2);
						}
					}

					eval_fwd(buf2, errors, iw, ih, kx0+kx, ky, t, p);
					load_vec16(curr, buf2, iw, ih, kx0+kx, ky, 0, 0);
					V16_SUB(curr, curr, t->pred);
					store_vec16(curr, errors, iw, ih, kx0+kx, ky);
				}
			}
		}

		for(int k=0;k<res;++k)
			buf[k<<2|kc]=dst[k];
	}

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
	_mm_free(g);
	_mm_free(b);
	_mm_free(p);
	_mm_free(t);
	free(dst);
	free(src);
}