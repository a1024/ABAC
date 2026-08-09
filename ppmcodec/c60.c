#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>
#if defined _MSC_VER || defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif

#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE

	#define C60TEST2
//	#define C60TEST

//	#define ANSVAL
//	#define FIFOVAL
#endif

	#define ENABLE_ICT		//better
	#define USE_SIMPLE_PRED		//better
//	#define USE_LEGALL53		//worse


enum
{
	DIST_DEFAULT=7,
	DIST_LL=3,

	DWT_SH_HQ=6,	//min dimension = 2<<sh
	DWT_SH=4,
	DWT_SCALE=12,

#ifndef USE_SIMPLE_PRED
	SHIFT=12,
	NPREDS=4,
#endif
	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,
	
	FSE_PROBBITS=12,
	NCTX=6,
};


#ifdef C60TEST
enum
{
	D1=5,
	D2=40,
};
typedef struct _ResultC60
{
	double bps, psnr;
} ResultC60;
static ResultC60 results[D2+1]={0};
#endif


//runtime
#if 1
#define CLAMP2(X, L, H) X=X<(L)?L:X, X=X>(H)?H:X
#ifdef _MSC_VER
#	define	ALIGN(N) __declspec(align(N))
#	define INLINE __forceinline static
#	define LZCNT32 _lzcnt_u32
#	define LZCNT64 _lzcnt_u64
#	define TZCNT32 _tzcnt_u32
#	define TZCNT64 _tzcnt_u64
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#	define LZCNT32 __builtin_clz
#	define LZCNT64 __builtin_clzll
#	define TZCNT32 __builtin_ctz
#	define TZCNT64 __builtin_ctzll
#endif
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
static double time_sec(void)
{
#ifdef _WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static uint8_t *g_image=0;
static double g_sqe[3]={0};
static void guide_save(uint8_t *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static double guide_psnr(double invres)
{
	double rmse[]=
	{
		sqrt(g_sqe[0]*invres),
		sqrt(g_sqe[1]*invres),
		sqrt(g_sqe[2]*invres),
		sqrt((g_sqe[0]+g_sqe[1]+g_sqe[2])*invres*(1./3)),
	};
	double psnr[]=
	{
		-20*log10(rmse[0]*(1./255)),
		-20*log10(rmse[1]*(1./255)),
		-20*log10(rmse[2]*(1./255)),
		-20*log10(rmse[3]*(1./255)),
	};
#if !defined C60TEST2
	printf(
		"RMSE  PSNR\n"
		"T %12.6lf  %12.6lf\n"
		"R %12.6lf  %12.6lf\n"
		"G %12.6lf  %12.6lf\n"
		"B %12.6lf  %12.6lf\n"
		, rmse[3], psnr[3]
		, rmse[0], psnr[0]
		, rmse[1], psnr[1]
		, rmse[2], psnr[2]
	);
#endif
	return psnr[3];
}
#endif
#endif


#ifdef ANSVAL
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	uint8_t data[];
} ANSVALNode;
static ANSVALNode *debugstack=0;
static int ansvalidx=0, ansvalmax=0;
static void ansval_push(const void *data, int esize, int count)
{
	int size=count*esize;
	ANSVALNode *node=(ANSVALNode*)malloc(sizeof(ANSVALNode)+size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, sizeof(ANSVALNode)+size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size);
	debugstack=node;
	++ansvalmax;
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const uint8_t *p=(const uint8_t*)data, *p2=(const uint8_t*)xdata;
	int size=count*esize, k;
	for(k=0;k<size;k+=esize)
	{
		int k2=esize-1;
		printf(" ");
		for(;k2>=0;--k2)
		{
			int val=p[k+k2];
			if(p2)
				val^=p2[k+k2];
			if(p2&&!val)
				printf("--");
			else
				printf("%02X", val);
		}
	}
	printf("\n");
}
static void ansval_check(const void *data, int esize, int count)
{
	--ansvalidx;
	if(!debugstack)
	{
		printf("Debug stack is empty\n");
		ansval_printr(data, esize, count, 0);
		CRASH("");
	}
	else if(debugstack->esize!=esize||debugstack->count!=count||memcmp(data, debugstack->data, esize*count))
	{
		printf("\n\nValidation Error  [enc ^ | v dec]\n");
		if(debugstack->above)
		{
			ANSVALNode *node=debugstack->above;
			if(node->above)
			{
				ANSVALNode *node2=node->above;
				printf("[%10d] Verified:   ", node2->idx);
				ansval_printr(node2->data, node2->esize, node2->count, 0);
			}
			printf("[%10d] Verified:   ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			printf("\n");
		}

		printf("[%10d] Original:   ", debugstack->idx);
		ansval_printr(debugstack->data, esize, count, 0);

		printf("[%10d] Corrupt:    ", debugstack->idx);
		ansval_printr(data, esize, count, 0);
		
		if(debugstack->esize==esize&&debugstack->count==count)
		{
			printf("[%10d] XOR:        ", debugstack->idx);
			ansval_printr(debugstack->data, esize, count, data);
		}
		if(debugstack->below)
		{
			ANSVALNode *node=debugstack->below;
			printf("\n");
			printf("[%10d] Below:      ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			if(node->below)
			{
				node=node->below;
				printf("[%10d] Below:      ", node->idx);
				ansval_printr(node->data, node->esize, node->count, 0);
			}
		}
		printf("\n\n");
		CRASH("");
	}
	if(debugstack->below)
		debugstack=debugstack->below;
}
#else
#define ansval_push(...)
#define ansval_check(...)
#endif
#ifdef FIFOVAL
static ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
static uint32_t *fifoval=0;
static void fifoval_enqueue(uint32_t val)
{
	if(fifoidx+1>=fifocap)
	{
		void *p=0;

		if(!fifocap)
			fifocap=1;
		fifocap<<=1;
		p=realloc(fifoval, fifocap*sizeof(uint32_t));
		if(!p)
		{
			CRASH("Alloc error");
			return;
		}
		fifoval=(uint32_t*)p;
	}
	fifoval[fifoidx++]=val;
}
static void fifoval_check(uint32_t val)
{
	uint32_t val0=fifoval[fifoidx2++];
	if(val!=val0)
	{
		--fifoidx2;
		printf(
			"\n"
			"FIFO Error  at %10lld,  remaining %10lld\n"
			"    0x%08X  !=  original 0x%08X\n"
			"\n"
			, fifoidx2
			, fifoidx-fifoidx2//current element was not decoded successfully
			, val, val0
		);
		for(int k=-32;k<32;++k)
		{
			ptrdiff_t idx=fifoidx2+k;
			if((size_t)idx<(size_t)fifoidx)
			{
				printf("%10td  0x%08X", idx, fifoval[idx]);
				if(idx<fifoidx2)
					printf("  OK");
				if(idx==fifoidx2)
					printf("  !=  corrupt 0x%08X", val);
				printf("\n");
			}
		}
		CRASH("");
	}
}
#endif

static const uint8_t ict_perm[]=
{
	0, 1, 2,
	2, 0, 1,
	1, 2, 0,
	2, 1, 0,
	1, 0, 2,
	0, 2, 1,
};
#define FIX12(X) (int)((X)*4096+0.5)
static int ict_select(uint8_t *image, int iw, int ih)
{
	int64_t sums[6]={0};
	int prev[6][3]={0};
	uint8_t *imptr, *end;
	int bestp=0;

	end=image+(intptr_t)3*iw*ih;
	imptr=image;
	while(imptr<end)
	{
		int c0, c1, c2, y, u, v;
		const uint8_t *pptr=ict_perm;
		for(int k=0;k<6;++k)
		{
			c0=imptr[pptr[0]]<<4;
			c1=imptr[pptr[1]]<<4;
			c2=imptr[pptr[2]]<<4;
			y=(+FIX12(0.299   )*c0+FIX12(0.587   )*c1+FIX12(0.114   )*c2)>>12;
			u=(-FIX12(0.168736)*c0-FIX12(0.331264)*c1+FIX12(0.5     )*c2)>>12;
			v=(+FIX12(0.5     )*c0-FIX12(0.418688)*c1-FIX12(0.081312)*c2)>>12;
			pptr+=3;
			sums[k]+=
				+(int64_t)abs(y-prev[k][0])
				+(int64_t)abs(u-prev[k][1])
				+(int64_t)abs(v-prev[k][2])
			;
			prev[k][0]=y;
			prev[k][1]=u;
			prev[k][2]=v;
		}
		imptr+=3;
	}
	for(int k=1;k<6;++k)
	{
		if(sums[bestp]>sums[k])
			bestp=k;
	}
	return bestp;
}
static void ict_fwd(uint8_t *image, int iw, int ih, int ict)
{
	intptr_t usize=(intptr_t)3*iw*ih;
	int idx0=ict_perm[3*ict+0];
	int idx1=ict_perm[3*ict+1];
	int idx2=ict_perm[3*ict+2];
	int c0, c1, c2, y, u, v;
	uint8_t *imptr=image+usize;
	int16_t *imptr2=(int16_t*)image+usize;

	while(imptr>image)
	{
		imptr-=3;
		imptr2-=3;
#ifdef ENABLE_ICT
		c0=imptr[idx0];
		c1=imptr[idx1];
		c2=imptr[idx2];
		y=(+FIX12(0.299   )*c0+FIX12(0.587   )*c1+FIX12(0.114   )*c2+(1<<12>>1))>>12;
		u=(-FIX12(0.168736)*c0-FIX12(0.331264)*c1+FIX12(0.5     )*c2+(1<<12>>1))>>12;
		v=(+FIX12(0.5     )*c0-FIX12(0.418688)*c1-FIX12(0.081312)*c2+(1<<12>>1))>>12;
		imptr2[0]=y-128;
		imptr2[1]=u;
		imptr2[2]=v;
#else
		c0=imptr[0];
		c1=imptr[1];
		c2=imptr[2];
		imptr2[0]=c0;
		imptr2[1]=c1;
		imptr2[2]=c2;
		(void)y;
		(void)u;
		(void)v;
#endif
	}
}
static void ict_inv(uint8_t *image, int iw, int ih, int ict)
{
	intptr_t usize=(intptr_t)3*iw*ih;
	int idx0=ict_perm[3*ict+0];
	int idx1=ict_perm[3*ict+1];
	int idx2=ict_perm[3*ict+2];
	int c0, c1, c2, y, u, v;
	uint8_t *imptr=image;
	int16_t *imptr2=(int16_t*)image, *imend=imptr2+usize;
	
#ifdef ENABLE_GUIDE
	memset(g_sqe, 0, sizeof(g_sqe));
#endif
	while(imptr2<imend)
	{
#ifdef ENABLE_ICT
		y=(imptr2[0]+128)<<12;
		u=imptr2[1];
		v=imptr2[2];
		c0=(y+FIX12(1.402)*v+(1<<12>>1))>>12;
		c1=(y-FIX12(0.344136)*u-FIX12(0.714136)*v+(1<<12>>1))>>12;
		c2=(y+FIX12(1.772)*u+(1<<12>>1))>>12;
		CLAMP2(c0, 0, 255);
		CLAMP2(c1, 0, 255);
		CLAMP2(c2, 0, 255);
		imptr[idx0]=c0;
		imptr[idx1]=c1;
		imptr[idx2]=c2;
#else
		c0=imptr2[0];
		c1=imptr2[1];
		c2=imptr2[2];
		CLAMP2(c0, 0, 255);
		CLAMP2(c1, 0, 255);
		CLAMP2(c2, 0, 255);
		imptr[0]=c0;
		imptr[1]=c1;
		imptr[2]=c2;
		(void)y;
		(void)u;
		(void)v;
#endif
#ifdef ENABLE_GUIDE
		{
			int diff;

			diff=g_image[imptr-image+0]-imptr[0]; g_sqe[0]+=diff*diff;
			diff=g_image[imptr-image+1]-imptr[1]; g_sqe[1]+=diff*diff;
			diff=g_image[imptr-image+2]-imptr[2]; g_sqe[2]+=diff*diff;
		}
#endif
		imptr+=3;
		imptr2+=3;
	}
}

static void ppm_save(const char *fn, const uint8_t *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
static void ppm_save_16bit(const char *fn, const int16_t *image, int iw, int ih, int vmax)
{
	size_t size=sizeof(uint16_t[3])*iw*ih;
	uint16_t *bigbuf=(uint16_t*)malloc(size);
	int vmax2=0;

	if(!bigbuf)
	{
		CRASH("Alloc error");
		return;
	}
	for(int k=0;k<size/sizeof(uint16_t);++k)
	{
		uint16_t val=image[k]^0x8000;
		bigbuf[k]=val<<8|val>>8;
		if(vmax2<val)vmax2=val;
	}
	if(!vmax)vmax=vmax2;
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n%d\n", iw, ih, vmax);
	//fprintf(fdst, "P6\n%d %d\n65535\n", iw, ih);
	fwrite(bigbuf, 1, size, fdst);
	fclose(fdst);
	free(bigbuf);
}
static void ppm_save_16as8(const char *fn, const int16_t *image, int iw, int ih)
{
	size_t size=sizeof(uint8_t[3])*iw*ih;
	uint8_t *buf=(uint8_t*)malloc(size);
	int vmin=0x7FFF, vmax=-0x8000;
	FILE *fdst=0;

	if(!buf)
	{
		CRASH("Alloc error");
		return;
	}
	for(int k=0;k<size;++k)
	{
		int val=(int16_t)image[k];
		if(vmin>val)vmin=val;
		if(vmax<val)vmax=val;
	}
	if(vmin<0)
	{
		vmin=-vmin;
		if(vmax<vmin)vmax=vmin;
	}
	if(vmax>0)
	{
		for(int k=0;k<size;++k)
			buf[k]=((int16_t)image[k]+vmax)*255/(2*vmax+1);
		//	buf[k]=((int16_t)image[k]-vmin)*255/(vmax-vmin);
		//	buf[k]=image[k]>>8^128;
	}
	else
		memset(buf, 128, size);
	fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(buf, 1, size, fdst);
	fclose(fdst);
	free(buf);
}
#define SCALEUP(X) (int)((X)*(1<<DWT_SCALE)+0.5)
	#define KAPPA 1.1496043988602
static const int32_t cdf97_coeffs[]=
{
	+SCALEUP(1.58613434342059),		//alpha		predict
	-SCALEUP(0.0529801185729),		//beta		update
	-SCALEUP(0.8829110755309),		//gamma		predict
	+SCALEUP(0.4435068520439),		//delta		update
	+SCALEUP(KAPPA), +SCALEUP(1/KAPPA),	//kappa		output gain is 1.89

//	+SCALEUP(1.01), +SCALEUP(1/1.01),
//	+SCALEUP(1.02), +SCALEUP(1/1.02),
//	+SCALEUP(1.03),	+SCALEUP(1/1.03),
//	+SCALEUP(1.04),	+SCALEUP(1/1.04),
//	+SCALEUP(1.05),	+SCALEUP(1/1.05),
//	+SCALEUP(1.06),	+SCALEUP(1/1.06),
//	+SCALEUP(1.07),	+SCALEUP(1/1.07),
//	+SCALEUP(1.08),	+SCALEUP(1/1.08),
};
static int dwt_getNiter(int iw, int ih, int dist)
{
	int sh=dist<=5?DWT_SH_HQ:DWT_SH;
	int niter=31^LZCNT32((iw|ih)>>sh);

	//niter=0;
	//while(iw>mindim&&ih>mindim)
	//{
	//	++niter;
	//	iw>>=1;
	//	ih>>=1;
	//}
	return niter;
}
INLINE void dwt1d_cdf97_predict(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t coeff)
{
	even[0]-=(odd[0]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	even[1]-=(odd[1]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	even[2]-=(odd[2]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	for(int k=3;k<nodd3;k+=3)//predict
	{
#if defined _MSC_VER && 0
		int p=(coeff*(odd[k-3+0]+odd[k+0])+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		int dst=even[k+0]-p;
		if(dst!=(int16_t)dst)
			CRASH("");
#endif
		even[k+0]-=(coeff*(odd[k-3+0]+odd[k+0])+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[k+1]-=(coeff*(odd[k-3+1]+odd[k+1])+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[k+2]-=(coeff*(odd[k-3+2]+odd[k+2])+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
	if(extraeven)
	{
		even[nodd3+0]-=(odd[nodd3-3+0]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[nodd3+1]-=(odd[nodd3-3+1]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[nodd3+2]-=(odd[nodd3-3+2]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
}
INLINE void dwt1d_cdf97_update(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t coeff)
{
	int count=nodd3-3*!extraeven;

	for(int k=0;k<count;k+=3)//update
	{
		odd[k+0]+=((even[k+0]+even[k+3+0])*coeff+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[k+1]+=((even[k+1]+even[k+3+1])*coeff+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[k+2]+=((even[k+2]+even[k+3+2])*coeff+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
	if(!extraeven)
	{
		odd[nodd3-3+0]+=(even[nodd3-3+0]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[nodd3-3+1]+=(even[nodd3-3+1]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[nodd3-3+2]+=(even[nodd3-3+2]*coeff*2+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
}
INLINE void dwt1d_cdf97_scale(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t co, int64_t ce)
{
	for(int k=0;k<nodd3;k+=3)
	{
		odd[k+0]=(odd[k+0]*co+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[k+1]=(odd[k+1]*co+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		odd[k+2]=(odd[k+2]*co+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[k+0]=(even[k+0]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[k+1]=(even[k+1]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[k+2]=(even[k+2]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
	if(extraeven)
	{
		even[nodd3+0]=(even[nodd3+0]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[nodd3+1]=(even[nodd3+1]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
		even[nodd3+2]=(even[nodd3+2]*ce+(1<<DWT_SCALE>>1))>>DWT_SCALE;
	}
}
INLINE void dwt1d_cdf97_shl(int16_t *buf, intptr_t count3, int sh)
{
	intptr_t k;

	for(k=0;k<count3;k+=3)
	{
		int c0=buf[k+0];
		int c1=buf[k+1];
		int c2=buf[k+2];
		buf[k+0]=c0<<sh;
		buf[k+1]=c1<<sh;
		buf[k+2]=c2<<sh;
	}
}
INLINE void dwt1d_cdf97_shr(int16_t *buf, intptr_t count3, int sh)
{
	intptr_t k;
	int bias;

	bias=1<<sh>>1;
	for(k=0;k<count3;k+=3)
	{
		int c0=buf[k+0];
		int c1=buf[k+1];
		int c2=buf[k+2];
		buf[k+0]=(c0+bias)>>sh;
		buf[k+1]=(c1+bias)>>sh;
		buf[k+2]=(c2+bias)>>sh;
	}
}
#if 0
INLINE void dwt1d_cdf97_predict(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t coeff)
{
	even[0]=(even[0]<<DWT_SCALE)-odd[0]*coeff*2;
	even[1]=(even[1]<<DWT_SCALE)-odd[1]*coeff*2;
	even[2]=(even[2]<<DWT_SCALE)-odd[2]*coeff*2;
	for(int k=3;k<nodd3;k+=3)//predict
	{
		even[k+0]=(even[k+0]<<DWT_SCALE)-coeff*(odd[k-3+0]+odd[k+0]);
		even[k+1]=(even[k+1]<<DWT_SCALE)-coeff*(odd[k-3+1]+odd[k+1]);
		even[k+2]=(even[k+2]<<DWT_SCALE)-coeff*(odd[k-3+2]+odd[k+2]);
	}
	if(extraeven)
	{
		even[nodd3+0]=(even[nodd3+0]<<DWT_SCALE)-odd[nodd3-3+0]*coeff*2;
		even[nodd3+1]=(even[nodd3+1]<<DWT_SCALE)-odd[nodd3-3+1]*coeff*2;
		even[nodd3+2]=(even[nodd3+2]<<DWT_SCALE)-odd[nodd3-3+2]*coeff*2;
	}
}
INLINE void dwt1d_cdf97_update(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t coeff)
{
	int count=nodd3-3*!extraeven;

	for(int k=0;k<count;k+=3)//update
	{
		odd[k+0]=(odd[k+0]<<DWT_SCALE)-(even[k+0]+even[k+3+0])*coeff;
		odd[k+1]=(odd[k+1]<<DWT_SCALE)-(even[k+1]+even[k+3+1])*coeff;
		odd[k+2]=(odd[k+2]<<DWT_SCALE)-(even[k+2]+even[k+3+2])*coeff;
	}
	if(!extraeven)
	{
		odd[nodd3-3+0]=(odd[nodd3-3+0]<<DWT_SCALE)-even[nodd3-3+0]*coeff*2;
		odd[nodd3-3+1]=(odd[nodd3-3+1]<<DWT_SCALE)-even[nodd3-3+1]*coeff*2;
		odd[nodd3-3+2]=(odd[nodd3-3+2]<<DWT_SCALE)-even[nodd3-3+2]*coeff*2;
	}
}
INLINE void dwt1d_cdf97_scale(int64_t *odd, int64_t *even, int nodd3, int extraeven, int64_t co, int64_t ce)
{
	for(int k=0;k<nodd3;k+=3)
	{
		odd[k+0]*=co;
		odd[k+1]*=co;
		odd[k+2]*=co;
		even[k+0]*=ce;
		even[k+1]*=ce;
		even[k+2]*=ce;
	}
	if(extraeven)
	{
		even[nodd3+0]*=ce;
		even[nodd3+1]*=ce;
		even[nodd3+2]*=ce;
	}
}
INLINE void dwt1d_cdf97_shr(int64_t *buf, int count3, int sh)
{
	for(int k=0;k<count3;k+=3)
	{
		int64_t c0=buf[k+0];
		int64_t c1=buf[k+1];
		int64_t c2=buf[k+2];
		int64_t mask0=c0>>63;
		int64_t mask1=c1>>63;
		int64_t mask2=c2>>63;
		c0^=mask0;
		c1^=mask1;
		c2^=mask2;
		c0-=mask0;
		c1-=mask1;
		c2-=mask2;
		c0>>=sh;
		c1>>=sh;
		c2>>=sh;
		c0^=mask0;
		c1^=mask1;
		c2^=mask2;
		c0-=mask0;
		c1-=mask1;
		c2-=mask2;
		buf[k+0]=c0;
		buf[k+1]=c1;
		buf[k+2]=c2;
	}
}
INLINE void dwt1d_cdf97_shl(int64_t *buf, int count3, int sh)
{
	for(int k=0;k<count3;k+=3)
	{
		int64_t c0=buf[k+0];
		int64_t c1=buf[k+1];
		int64_t c2=buf[k+2];
		buf[k+0]=c0<<sh;
		buf[k+1]=c1<<sh;
		buf[k+2]=c2<<sh;
	}
}
#endif
#ifdef _MSC_VER
#define DEBUGCHECK(A, B, C) if(A<-0x8000||A>0x7FFF||B<-0x8000||B>0x7FFF||C<-0x8000||C>0x7FFF)CRASH("%d %d %d", A, B, C)
#else
#define DEBUGCHECK(...)
#endif
static void dwt1d_cdf97_fwd(int16_t *buffer, int count, int stride, int64_t *b2, int level, int dist)
{
	int nodd=count>>1, extraeven=count&1, nodd3=nodd*3;
	int64_t *odd=b2, *even=b2+nodd3;
	
	for(int k=0, ks=0;k<nodd3;k+=3, ks+=stride*2)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k+0]=buffer[ks+0];
		even[k+1]=buffer[ks+1];
		even[k+2]=buffer[ks+2];
		odd[k+0]=buffer[ks+stride+0];
		odd[k+1]=buffer[ks+stride+1];
		odd[k+2]=buffer[ks+stride+2];
	}
	if(extraeven)
	{
		even[nodd3+0]=buffer[stride*(count-1)+0];
		even[nodd3+1]=buffer[stride*(count-1)+1];
		even[nodd3+2]=buffer[stride*(count-1)+2];
	}
	
#ifdef USE_LEGALL53
	dwt1d_cdf97_predict(odd, even, nodd3, extraeven, 0x0800);
	dwt1d_cdf97_update (odd, even, nodd3, extraeven, 0x0400);
#else
	//if(level>2)
	//{
	//	dwt1d_cdf97_predict(odd, even, nodd3, extraeven, 0x0800);
	//	dwt1d_cdf97_update (odd, even, nodd3, extraeven, 0x0400);
	//}
	//else
	{
		dwt1d_cdf97_predict(odd, even, nodd3, extraeven, cdf97_coeffs[0]);
		dwt1d_cdf97_update (odd, even, nodd3, extraeven, cdf97_coeffs[1]);
		dwt1d_cdf97_predict(odd, even, nodd3, extraeven, cdf97_coeffs[2]);
		dwt1d_cdf97_update (odd, even, nodd3, extraeven, cdf97_coeffs[3]);
		if(dist>10)
			dwt1d_cdf97_scale  (odd, even, nodd3, extraeven, cdf97_coeffs[5], cdf97_coeffs[4]);//scale: co=1.14, ce=1/1.14 in J2K
	}
#endif
	(void)&dwt1d_cdf97_scale;

	count*=3;
	for(int k=0, ks=0;k<count;k+=3, ks+=stride)
	{
		int64_t v0=b2[k+0];
		int64_t v1=b2[k+1];
		int64_t v2=b2[k+2];
	//	DEBUGCHECK(v0, v1, v2);
		CLAMP2(v0, -0x8000, 0x7FFF);
		CLAMP2(v1, -0x8000, 0x7FFF);
		CLAMP2(v2, -0x8000, 0x7FFF);
		buffer[ks+0]=(int16_t)v0;
		buffer[ks+1]=(int16_t)v1;
		buffer[ks+2]=(int16_t)v2;
	}
}
static void dwt1d_cdf97_inv(int16_t *buffer, int count, int stride, int64_t *b2, int level, int dist)
{
	int nodd=count>>1, extraeven=count&1, nodd3=nodd*3;
	int64_t *odd=b2, *even=b2+nodd3;

	for(int k=0, ks=0;k<3*count;k+=3, ks+=stride)
	{
		b2[k+0]=buffer[ks+0];
		b2[k+1]=buffer[ks+1];
		b2[k+2]=buffer[ks+2];
	}
	
#ifdef USE_LEGALL53
	dwt1d_cdf97_update (odd, even, nodd3, extraeven, -0x0400);
	dwt1d_cdf97_predict(odd, even, nodd3, extraeven, -0x0800);
#else
	//if(level>2)
	//{
	//	dwt1d_cdf97_update (odd, even, nodd3, extraeven, -0x0400);
	//	dwt1d_cdf97_predict(odd, even, nodd3, extraeven, -0x0800);
	//}
	//else
	{
		if(dist>10)
			dwt1d_cdf97_scale  (odd, even, nodd3, extraeven, cdf97_coeffs[4], cdf97_coeffs[5]);//unscale: co=1/1.14, ce=1.14 in J2K
		dwt1d_cdf97_update (odd, even, nodd3, extraeven, -cdf97_coeffs[3]);
		dwt1d_cdf97_predict(odd, even, nodd3, extraeven, -cdf97_coeffs[2]);
		dwt1d_cdf97_update (odd, even, nodd3, extraeven, -cdf97_coeffs[1]);
		dwt1d_cdf97_predict(odd, even, nodd3, extraeven, -cdf97_coeffs[0]);
	}
#endif

	for(int k=0, ks=0;k<nodd3;k+=3, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		int64_t e0=even[k+0];
		int64_t e1=even[k+1];
		int64_t e2=even[k+2];
		int64_t o0=odd[k+0];
		int64_t o1=odd[k+1];
		int64_t o2=odd[k+2];
	//	DEBUGCHECK(e0, e1, e2);
	//	DEBUGCHECK(o0, o1, o2);
		CLAMP2(e0, -0x8000, 0x7FFF);
		CLAMP2(e1, -0x8000, 0x7FFF);
		CLAMP2(e2, -0x8000, 0x7FFF);
		CLAMP2(o0, -0x8000, 0x7FFF);
		CLAMP2(o1, -0x8000, 0x7FFF);
		CLAMP2(o2, -0x8000, 0x7FFF);
		buffer[ks+0]=(int16_t)e0;
		buffer[ks+1]=(int16_t)e1;
		buffer[ks+2]=(int16_t)e2;
		buffer[ks+stride+0]=(int16_t)o0;
		buffer[ks+stride+1]=(int16_t)o1;
		buffer[ks+stride+2]=(int16_t)o2;
	}
	if(extraeven)
	{
		buffer[stride*(count-1)+0]=(int16_t)even[nodd3+0];
		buffer[stride*(count-1)+1]=(int16_t)even[nodd3+1];
		buffer[stride*(count-1)+2]=(int16_t)even[nodd3+2];
	}
}
static void dwt_cdf97(int16_t *image, int iw, int ih, int niter, int dist, int fwd)
{
	int64_t *temp=(int64_t*)malloc((iw>ih?iw:ih)*sizeof(int64_t[3]));
	intptr_t size=(intptr_t)3*iw*ih;
	int rowstride=3*iw;

	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	if(fwd)
	{
		dwt1d_cdf97_shl(image, size, 2);
		for(int it=0;it<niter-1;++it)
		{
			int w2=iw>>it, h2=ih>>it;

			for(int ky=0;ky<h2;++ky)//horizontal DWT
				dwt1d_cdf97_fwd(image+rowstride*ky, w2, 3, temp, it, dist);

			for(int kx=0;kx<w2;++kx)//vertical DWT
				dwt1d_cdf97_fwd(image+3*kx, h2, rowstride, temp, it, dist);
		}
		dwt1d_cdf97_shr(image, size, 2);
	}
	else
	{
		dwt1d_cdf97_shl(image, size, 2);
		for(int it=niter-2;it>=0;--it)
		{
			int w2=iw>>it, h2=ih>>it;

			for(int kx=0;kx<w2;++kx)//vertical IDWT
				dwt1d_cdf97_inv(image+3*kx, h2, rowstride, temp, it, dist);

			for(int ky=0;ky<h2;++ky)//horizontal IDWT
				dwt1d_cdf97_inv(image+rowstride*ky, w2, 3, temp, it, dist);
		}
		dwt1d_cdf97_shr(image, size, 2);
	}
	free(temp);
}

static void getqsteps(int16_t *image, int iw, int ih, int x1, int x2, int y1, int y2, int dist, int *qsteps)
{
	const int minstep=3<<16;
	uint64_t sum[3]={0};
	int rowstride=3*iw, count=(x2-x1)*(y2-y1);
	int vmax[3]={0};

	if(!count)
		return;
	for(int ky=y1;ky<y2;++ky)
	{
		int16_t *imptr=image+rowstride*ky;
		for(int kx=x1;kx<x2;++kx)
		{
			int c0=abs(imptr[0]);
			int c1=abs(imptr[1]);
			int c2=abs(imptr[2]);
			imptr+=3;
			if(vmax[0]<c0)vmax[0]=c0;
			if(vmax[1]<c1)vmax[1]=c1;
			if(vmax[2]<c2)vmax[2]=c2;
			sum[0]+=c0;
			sum[1]+=c1;
			sum[2]+=c2;
		}
	}
	qsteps[0]=(int)((sum[0]*dist<<16)/((uint64_t)count));
	qsteps[1]=(int)((sum[1]*dist<<16)/((uint64_t)count));
	qsteps[2]=(int)((sum[2]*dist<<16)/((uint64_t)count));
	if(qsteps[0]<vmax[0]<<10)qsteps[0]=vmax[0]<<10;//qstep >= 2*vmax for signbit
	if(qsteps[1]<vmax[1]<<10)qsteps[1]=vmax[1]<<10;
	if(qsteps[2]<vmax[2]<<10)qsteps[2]=vmax[2]<<10;
	if(qsteps[0]<minstep)qsteps[0]=minstep;
	if(qsteps[1]<minstep)qsteps[1]=minstep;
	if(qsteps[2]<minstep)qsteps[2]=minstep;
	qsteps[2]=qsteps[1]=qsteps[0]=(int)round(cbrt((double)qsteps[0]*qsteps[1]*qsteps[2]));
//	qsteps[2]=qsteps[1]=qsteps[0]=(qsteps[0]+qsteps[1]+qsteps[2])/3;
	(void)dist;
}
static void quantization(int16_t *image, int iw, int ih, int x1, int x2, int y1, int y2, const int *qsteps, int fwd)
{
	int rowstride=3*iw;
	if(fwd)
	{
		int invsteps[]=
		{
			(int)((0x100000000+qsteps[0]-1)/qsteps[0]),
			(int)((0x100000000+qsteps[1]-1)/qsteps[1]),
			(int)((0x100000000+qsteps[2]-1)/qsteps[2]),
		};
#if defined LOUD && 0
		printf(
			"Qstep/inv  %12.6lf / %12.6lf    %12.6lf / %12.6lf    %12.6lf / %12.6lf\n"
			, qsteps[0]/65536., invsteps[0]/65536.
			, qsteps[1]/65536., invsteps[1]/65536.
			, qsteps[2]/65536., invsteps[2]/65536.
		);
#endif
		for(int ky=y1;ky<y2;++ky)
		{
			int16_t *imptr=image+rowstride*ky+3*x1;
			for(int kx=x1;kx<x2;++kx)
			{
				int c0=imptr[0];
				int c1=imptr[1];
				int c2=imptr[2];
				int mask0=c0>>31;
				int mask1=c1>>31;
				int mask2=c2>>31;

				c0^=mask0;
				c1^=mask1;
				c2^=mask2;
				c0-=mask0;
				c1-=mask1;
				c2-=mask2;
				c0=c0*invsteps[0]>>16;
				c1=c1*invsteps[1]>>16;
				c2=c2*invsteps[2]>>16;
				c0-=c0>0;//deadzone
				c1-=c1>0;
				c2-=c2>0;
				c0^=mask0;
				c1^=mask1;
				c2^=mask2;
				c0-=mask0;
				c1-=mask1;
				c2-=mask2;
				c0=c0<<1^c0>>31;
				c1=c1<<1^c1>>31;
				c2=c2<<1^c2>>31;
				imptr[0]=c0;
				imptr[1]=c1;
				imptr[2]=c2;
				imptr+=3;
			}
		}
	}
	else
	{
		for(int ky=y1;ky<y2;++ky)
		{
			int16_t *imptr=image+rowstride*ky+3*x1;
			for(int kx=x1;kx<x2;++kx)
			{
				int c0=imptr[0];
				int c1=imptr[1];
				int c2=imptr[2];
				c0=c0>>1^c0<<31>>31;
				c1=c1>>1^c1<<31>>31;
				c2=c2>>1^c2<<31>>31;
				c0=(int)((int64_t)(2*c0+3*((c0>0)-(c0<0)))*qsteps[0]>>(16+1));//bias 0.5
				c1=(int)((int64_t)(2*c1+3*((c1>0)-(c1<0)))*qsteps[1]>>(16+1));
				c2=(int)((int64_t)(2*c2+3*((c2>0)-(c2<0)))*qsteps[2]>>(16+1));
				imptr[0]=c0;
				imptr[1]=c1;
				imptr[2]=c2;
				imptr+=3;
			}
		}
	}
}
static void predictLL(int16_t *image, int iw, int ih, int llw, int llh, int16_t *pixels, int psize, int dist, int fwd)
{
#ifndef USE_SIMPLE_PRED
	int32_t coeffs[3][NPREDS]={0}, bias[3]={1<<SHIFT>>1, 1<<SHIFT>>1, 1<<SHIFT>>1};
#endif
	int rowstride=3*iw;
	int invdist=((1<<16)+dist-1)/dist;

	memset(pixels, 0, psize);
	for(int ky=0;ky<llh;++ky)
	{
		int16_t *imptr=image+rowstride*ky;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<llw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					eN	=rows[1][1+0*NCH*NROWS*NVAL],
					eNE	=rows[1][1+1*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int pred, vmin, vmax, curr;
#ifndef USE_SIMPLE_PRED
				int estim[]=
				{
					2*W-WW,
					2*N-NN,
					N+W-NW,
					N+NE-NNE,
				};

				pred=(bias[kc]
					+coeffs[kc][0]*estim[0]
					+coeffs[kc][1]*estim[1]
					+coeffs[kc][2]*estim[2]
					+coeffs[kc][3]*estim[3]
				)>>SHIFT;
				vmax=N; vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
#else
				//pred=W-((eW[kc]+8)>>4);

				pred=N+W-NW+((eW-eN+eNE+4)>>3);
				vmax=N; vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(pred, vmin, vmax);

				//pred=(3*(N+W)-2*NW+2)>>2;

				//pred=(N+W)>>1;
#endif
				curr=*imptr;
				if(fwd)
				{
					curr-=pred;
					curr=(curr*invdist>>16)-(curr>>31);
					*imptr=curr<<1^curr>>31;
				}
				else
					curr=curr>>1^curr<<31>>31;
				curr=curr*dist+pred;
				if(!fwd)
					*imptr=curr;
				++imptr;
#ifndef USE_SIMPLE_PRED
				pred=(curr>pred)-(curr<pred);
				bias[kc]+=pred;
				coeffs[kc][0]+=pred*estim[0];
				coeffs[kc][1]+=pred*estim[1];
				coeffs[kc][2]+=pred*estim[2];
				coeffs[kc][3]+=pred*estim[3];
#endif
				rows[0][0]=curr;
				rows[0][1]=curr-pred;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				(void)NN	;
				(void)NNE	;
				(void)NW	;
				(void)N		;
				(void)NE	;
				(void)NEEE	;
				(void)WWW	;
				(void)WW	;
				(void)W		;
				(void)eN	;
				(void)eNE	;
				(void)eW	;
			}
#if 0
			int y=imptr[0];
			int u=imptr[1];
			int v=imptr[2];
			
			NW[0]=imptr[0-1*3-1*rowstride];
			NW[1]=imptr[1-1*3-1*rowstride];
			NW[2]=imptr[2-1*3-1*rowstride];
			N [0]=imptr[0+0*3-1*rowstride];
			N [1]=imptr[1+0*3-1*rowstride];
			N [2]=imptr[2+0*3-1*rowstride];
			W [0]=imptr[0-1*3-0*rowstride];
			W [1]=imptr[1-1*3-0*rowstride];
			W [2]=imptr[2-1*3-0*rowstride];
			vmax[0]=N[0];
			vmax[1]=N[1];
			vmax[2]=N[2];
			vmin[0]=W[0];
			vmin[1]=W[1];
			vmin[2]=W[2];
			if(N[0]<W[0])vmin[0]=N[0], vmin[0]=W[0];
			if(N[1]<W[1])vmin[1]=N[1], vmin[1]=W[1];
			if(N[2]<W[2])vmin[2]=N[2], vmin[2]=W[2];
			pred[0]=N[0]+W[0]-NW[0];
			pred[1]=N[1]+W[1]-NW[1];
			pred[2]=N[2]+W[2]-NW[2];
			CLAMP2(pred[0], vmin[0], vmax[0]);
			CLAMP2(pred[1], vmin[1], vmax[1]);
			CLAMP2(pred[2], vmin[2], vmax[2]);
			pred[0]^=fwdmask;
			pred[1]^=fwdmask;
			pred[2]^=fwdmask;
			pred[0]-=fwdmask;
			pred[1]-=fwdmask;
			pred[2]-=fwdmask;
			pred[0]+=y;
			pred[1]+=u;
			pred[2]+=v;
#endif
		}
	}
}

//FSE		https://github.com/Cyan4973/FiniteStateEntropy
#if 1
typedef struct _FSE_CState_t
{
	ptrdiff_t value;
	const void *stateTable, *symbolTT;
	uint32_t stateLog;
} FSE_CState_t;
typedef struct _FSE_DState_t
{
	size_t state;
	const void *table;
} FSE_DState_t;
typedef struct _FSE_symbolCompressionTransform
{
	int32_t deltaFindState;
	uint32_t deltaNbBits;
} FSE_symbolCompressionTransform;
typedef uint32_t FSE_CTable;
typedef uint32_t FSE_DTable;
typedef struct _FSE_DTableHeader
{
	uint16_t tableLog, fastMode;
} FSE_DTableHeader;
typedef struct _FSE_decode_t
{
	uint16_t newState;
	uint8_t symbol, nbBits;
} FSE_decode_t;
#define FSE_MAX_SYMBOL_VALUE 255
#define FSE_FUNCTION_TYPE uint8_t
#define FSE_TABLELOG_ABSOLUTE_MAX 15
#define FSE_MAX_MEMORY_USAGE 14
#define FSE_DEFAULT_MEMORY_USAGE 13
#define FSE_MIN_TABLELOG 5
#define FSE_MAX_TABLELOG (FSE_MAX_MEMORY_USAGE-2)
#define FSE_DEFAULT_TABLELOG (FSE_DEFAULT_MEMORY_USAGE-2)
#define FSE_TABLESTEP(tableSize) (((tableSize)>>1)+((tableSize)>>3)+3)
#define FSE_CTABLE_SIZE_U32(maxTableLog, maxSymbolValue) (1+(1<<((maxTableLog)-1))+(((maxSymbolValue)+1)*2))
#define FSE_DTABLE_SIZE_U32(maxTableLog) (1+(1<<(maxTableLog)))
#define FSE_DECODE_TYPE FSE_decode_t
typedef struct _FSEo1Tables
{
	FSE_symbolCompressionTransform symbolTT[NCH*NCTX][256];
	uint16_t stateTable[NCH*NCTX][1<<FSE_PROBBITS];
	FSE_decode_t decodeTable[NCH*NCTX][1<<FSE_PROBBITS];
	
	uint8_t wksp[1<<14>>2];
} FSEo1Tables;
#endif

//LIFO Bit Packer inspired by FPC
#if 1
typedef struct _BitPackerLIFO
{
	int64_t bitidx;
	uint64_t cache;
	uint8_t *start, *ptr, *end;//fwd enc / bwd dec
} BitPackerLIFO;
INLINE void bitpacker_enc_init(BitPackerLIFO *ec, uint8_t *start, uint8_t *end)
{
	ec->bitidx=0;
	ec->cache=0;
	ec->start=start;
	ec->ptr=start+sizeof(uint64_t);
	ec->end=end;
	*(uint64_t*)start=0;
}
INLINE void bitpacker_enc(BitPackerLIFO *ec, int64_t bitcount, uint64_t sym)//1 <= bitcount <= 56
{
	ansval_push(&sym, sizeof(sym), 1);
	ec->cache|=sym<<ec->bitidx;
	ec->bitidx+=bitcount;
	*(uint64_t*)ec->ptr=ec->cache;
	ec->cache>>=ec->bitidx&56;
	ec->ptr+=ec->bitidx>>3;
	ec->bitidx&=7;
}
INLINE uint8_t* bitpacker_enc_flush(BitPackerLIFO *ec)
{
	bitpacker_enc(ec, 1, 1);
	++ec->ptr;
	return ec->ptr;
}
INLINE uint64_t bitpacker_dec(BitPackerLIFO *ec, int64_t bitcount)
{
	uint64_t sym;

	ec->bitidx-=bitcount;
	sym=((uint64_t*)ec->ptr)[-1]>>ec->bitidx&((1ULL<<bitcount)-1);
	ansval_check(&sym, sizeof(sym), 1);
	ec->ptr-=(63-ec->bitidx)>>3;
	ec->bitidx|=56;
	return sym;
}
INLINE void bitpacker_dec_init(BitPackerLIFO *ec, uint8_t *start, uint8_t *end)
{
	ec->start=start;
	ec->ptr=end;
	ec->end=end;
	ec->bitidx=64-LZCNT64(((uint64_t*)ec->ptr)[-1]);
	bitpacker_dec(ec, 1);
}
#endif

static void enc_packhist(BitPackerLIFO *ec, const int *hist, uint64_t bypassmask, int ctxidx, int userans, int probbits)//histogram must be normalized to PROBBITS, with spike at 128
{
	int sum, cdfW, CDFlevels, startsym, ks;
	uint16_t hist2[256], CDF[257];

	if(bypassmask>>ctxidx&1)
		return;
	sum=0;
	if(userans)//32-bit hist
	{
		for(ks=0;ks<256;++ks)
			hist2[ks]=hist[ks];
	}
	else//16-bit hist
		memcpy(hist2, hist, sizeof(hist2));
	for(ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist2[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<probbits;
	
	cdfW=CDF[0];
	CDFlevels=1<<probbits;
	startsym=0;
	for(ks=1;ks<=256;++ks)//push GR.k
	{
		int next=CDF[ks], freq=next-cdfW;
		int nbypass=31-LZCNT32(CDFlevels);
		if(ks>1)
			nbypass-=7;
		if(nbypass<0)
			nbypass=0;
		CDF[ks]=nbypass<<probbits|freq;
		cdfW=next;
		CDFlevels-=freq;
		startsym=ks;
		if(!CDFlevels)
			break;
	}
	for(ks=startsym;ks>0;--ks)//encode GR
	{
		int freq=CDF[ks], nbypass=freq>>probbits;
		freq&=(1<<probbits)-1;
		int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
#ifdef ANSVAL
		//ansval_push(&freq, sizeof(freq), 1);
#endif
		if(nbypass)
			bitpacker_enc(ec, nbypass, bypass);
		bitpacker_enc(ec, 1, 1);
		while(nzeros)
		{
			bitpacker_enc(ec, 1, 0);
			--nzeros;
		}
#ifdef ANSVAL
		//ansval_push(&ks, sizeof(ks), 1);
#endif
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, uint64_t bypassmask, int ctxidx, uint16_t *hist, int probbits)
{
	if(bypassmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<probbits)/256;
	}
	else
	{
		uint16_t CDF[257]={0};
		int CDFlevels=1<<probbits;
		CDF[0]=0;
		for(int ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=31-LZCNT32(CDFlevels);
			int ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANSVAL
			//ansval_check(&ks2, sizeof(ks2), 1);
#endif
			int bit=0;
			do
			{
				bit=(int)bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|(int)bitpacker_dec(ec, nbypass);
#ifdef ANSVAL
			//ansval_check(&freq, sizeof(freq), 1);
#endif

			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					CRASH("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			CRASH("CDF unpack error");
		for(int ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	//int sum=0;
	//for(int ks=0;ks<256;++ks)//integrate
	//{
	//	int freq=hist[ks];
	//	hist[ks]=sum;
	//	sum+=freq;
	//}
	//hist[256]=1<<PROBBITS;
	//for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	//{
	//	int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
	//	int val=(freq<<PROBBITS|0)<<8|ks;
	//	for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
	//		CDF2sym[ks2]=val;
	//}
}

static const uint64_t g_bitmask[]=
{
	0x0,			0x1,			0x3,			0x7,
	0xF,			0x1F,			0x3F,			0x7F,
	0xFF,			0x1FF,			0x3FF,			0x7FF,
	0xFFF,			0x1FFF,			0x3FFF,			0x7FFF,
	0xFFFF,			0x1FFFF,		0x3FFFF,		0x7FFFF,
	0xFFFFF,		0x1FFFFF,		0x3FFFFF,		0x7FFFFF,
	0xFFFFFF,		0x1FFFFFF,		0x3FFFFFF,		0x7FFFFFF,
	0xFFFFFFF,		0x1FFFFFFF,		0x3FFFFFFF,		0x7FFFFFFF,
	0xFFFFFFFF,		0x1FFFFFFFF,		0x3FFFFFFFF,		0x7FFFFFFFF,
	0xFFFFFFFFF,		0x1FFFFFFFFF,		0x3FFFFFFFFF,		0x7FFFFFFFFF,
	0xFFFFFFFFFF,		0x1FFFFFFFFFF,		0x3FFFFFFFFFF,		0x7FFFFFFFFFF,
	0xFFFFFFFFFFF,		0x1FFFFFFFFFFF,		0x3FFFFFFFFFFF,		0x7FFFFFFFFFFF,
	0xFFFFFFFFFFFF,		0x1FFFFFFFFFFFF,	0x3FFFFFFFFFFFF,	0x7FFFFFFFFFFFF,
	0xFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFF,	0x3FFFFFFFFFFFFF,	0x7FFFFFFFFFFFFF,
	0xFFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFFF,	0x3FFFFFFFFFFFFFF,	0x7FFFFFFFFFFFFFF,
	0xFFFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFFFF,	0x3FFFFFFFFFFFFFFF,	0x7FFFFFFFFFFFFFFF,
};
static int32_t hists[3*NCTX*256];
static int16_t nhists[3*NCTX*256];
static void normalizehist(const uint32_t *hist, uint16_t *nhist, int probbits)
{
	int hsum=0, nusedlevels=0;
	for(int ks=0;ks<256;++ks)//faster than maintaining hist sum
	{
		int freq=hist[ks];
		hsum+=freq;
		nusedlevels+=freq!=0;
	}
	int64_t rsum=(((1LL<<probbits)-nusedlevels)<<24)/hsum;//adaptive: allow all symbols
	uint16_t CDF[257]={0};
	for(int ks=0, c=0, c2=0;ks<256;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*rsum>>24)+c2;
		c+=freq;
		c2+=freq!=0;
	}
	CDF[256]=1<<probbits;
	for(int ks=0;ks<256;++ks)
		nhist[ks]=CDF[ks+1]-CDF[ks];
}
static void subband_code_fse(int16_t *image
	, int iw, int ih
	, int x1, int x2, int y1, int y2
	, int16_t *pixels, int psize
	, int fwd
	, uint8_t *ctxbuf	//ctxbuf is size iw*ih bytes
	, uint8_t **pstreamptr, uint8_t *streamend
)
{
	enum
	{
		RBITS=2,
	};
	int kx=0, ky=0, kc=0;
	uint64_t degenmask=0, bypassmask=0;
	FSEo1Tables *fse=0;
	BitPackerLIFO ec={0};
	uint32_t states[3]={0};
	int pixelcount=(y2-y1)*(x2-x1);
	
	fse=(FSEo1Tables*)malloc(sizeof(*fse));//FIXME don't malloc here
	if(!fse)
	{
		CRASH("Alloc error");
		return;
	}
	memset(fse, 0, sizeof(*fse));
	memset(pixels, 0, psize);
	if(fwd)
	{
		uint8_t *ctxptr=ctxbuf;

		memset(hists, 0, sizeof(hists));
		for(ky=y1;ky<y2;++ky)
		{
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			};
			int16_t *imptr=image+3*(iw*ky+x1);
			for(kx=x1;kx<x2;++kx)
			{
				int
					eNW0	=rows[1][0-1*NCH*NROWS*NVAL],
					eNW1	=rows[1][1-1*NCH*NROWS*NVAL],
					eNW2	=rows[1][2-1*NCH*NROWS*NVAL],
					eN0	=rows[1][0+0*NCH*NROWS*NVAL],
					eN1	=rows[1][1+0*NCH*NROWS*NVAL],
					eN2	=rows[1][2+0*NCH*NROWS*NVAL],
					eNE0	=rows[1][0+1*NCH*NROWS*NVAL],
					eNE1	=rows[1][1+1*NCH*NROWS*NVAL],
					eNE2	=rows[1][2+1*NCH*NROWS*NVAL],
					eW0	=rows[0][0-1*NCH*NROWS*NVAL],
					eW1	=rows[0][1-1*NCH*NROWS*NVAL],
					eW2	=rows[0][2-1*NCH*NROWS*NVAL];
				int ctx0=2*((eW0|eN0)!=0)+((eNW0|eNE0)!=0);
				int ctx1=2*((eW1|eN1)!=0)+((eNW1|eNE1)!=0);
				int ctx2=2*((eW2|eN2)!=0)+((eNW2|eNE2)!=0);
				//int ctx0=4*(eN0!=0)+2*(eW0!=0)+((eNW0|eNE0)!=0);
				//int ctx1=4*(eN1!=0)+2*(eW1!=0)+((eNW1|eNE1)!=0);
				//int ctx2=4*(eN2!=0)+2*(eW2!=0)+((eNW2|eNE2)!=0);
				//int ctx0=(eW0|eN0)!=0;
				//int ctx1=(eW1|eN1)!=0;
				//int ctx2=(eW2|eN2)!=0;
				//int ctx0=31^LZCNT32(eW0+1);
				//int ctx1=31^LZCNT32(eW1+1);
				//int ctx2=31^LZCNT32(eW2+1);
				if(ctx0>NCTX-1)ctx0=NCTX-1;
				if(ctx1>NCTX-1)ctx1=NCTX-1;
				if(ctx2>NCTX-1)ctx2=NCTX-1;
				int sym0=imptr[0];
				int sym1=imptr[1];
				int sym2=imptr[2];
				ctx0+=0*NCTX;
				ctx1+=1*NCTX;
				ctx2+=2*NCTX;
				++hists[ctx0*256+sym0];
				++hists[ctx1*256+sym1];
				++hists[ctx2*256+sym2];
				ctxptr[0]=ctx0;
				ctxptr[1]=ctx1;
				ctxptr[2]=ctx2;
				imptr+=3;
				ctxptr+=3;
				rows[0][0+0*NCH*NROWS*NVAL]=sym0;
				rows[0][1+0*NCH*NROWS*NVAL]=sym1;
				rows[0][2+0*NCH*NROWS*NVAL]=sym2;
				//{
				//	int
				//		eNEE0	=rows[1][0+2*NCH*NROWS*NVAL],
				//		eNEE1	=rows[1][1+2*NCH*NROWS*NVAL],
				//		eNEE2	=rows[1][2+2*NCH*NROWS*NVAL],
				//		eNEEE0	=rows[1][0+3*NCH*NROWS*NVAL],
				//		eNEEE1	=rows[1][1+3*NCH*NROWS*NVAL],
				//		eNEEE2	=rows[1][2+3*NCH*NROWS*NVAL];
				//	rows[0][0+0*NCH*NROWS*NVAL]=(2*eW0+(sym0<<RBITS)+(eNEE0>eNEEE0?eNEE0:eNEEE0))>>2;
				//	rows[0][1+0*NCH*NROWS*NVAL]=(2*eW1+(sym1<<RBITS)+(eNEE1>eNEEE1?eNEE1:eNEEE1))>>2;
				//	rows[0][2+0*NCH*NROWS*NVAL]=(2*eW2+(sym2<<RBITS)+(eNEE2>eNEEE2?eNEE2:eNEEE2))>>2;
				//}
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				rows[2]+=NCH*NROWS*NVAL;
				rows[3]+=NCH*NROWS*NVAL;
			}
		}
		for(kc=0;kc<3*NCTX;++kc)
		{
			enum
			{
				TABLESIZE=1<<FSE_PROBBITS,
				TABLEMASK=TABLESIZE-1,
				STEP=FSE_TABLESTEP(TABLESIZE),
				VMAX=255,
			};
			int16_t *currnhist=nhists+256*kc;
			int32_t *hcurr=hists+256*kc;
			int sum=0, nlevels=0, vmax=0, fmax=0;
			FSE_FUNCTION_TYPE *const tableSymbol=(FSE_FUNCTION_TYPE*)fse->wksp;
			uint16_t *tableU16=fse->stateTable[kc];
			FSE_symbolCompressionTransform *const symbolTT=fse->symbolTT[kc];
			uint32_t CDF[FSE_MAX_SYMBOL_VALUE+2];
			int highthreshold=TABLESIZE-1;
			
			//FSE_PROBBITS<=15 required for threshold strategy to work
			if(FSE_PROBBITS>15||((size_t)1<<FSE_PROBBITS)*sizeof(FSE_FUNCTION_TYPE)>sizeof(fse->wksp))
			{
				CRASH("tableLog too large");
				return;
			}
			for(int ks=0;ks<256;++ks)
			{
				int freq=hcurr[ks];
				sum+=freq;
				nlevels+=freq!=0;
				if(fmax<freq)
				{
					fmax=freq;
					vmax=ks;
				}
			}
			if(nlevels==1)//degenerate distribution
			{
				memset(currnhist, 0, sizeof(int16_t[256]));
				currnhist[vmax]=(1<<FSE_PROBBITS)-1;
				currnhist[(uint8_t)(vmax+1)]=1;
				degenmask|=1LL<<kc;
			}
			else if(sum<=256)//bypass distribution
			{
				for(int ks=0;ks<256;++ks)
					currnhist[ks]=1<<FSE_PROBBITS>>8;
				bypassmask|=1LL<<kc;
			}
			else//ordinary distribution
				normalizehist((uint32_t*)hcurr, (uint16_t*)currnhist, FSE_PROBBITS);
			
#if defined LOUD && 0
			if(sum)
			{
				double p0a=(double)hcurr[0]/sum;
				double p0b=(double)currnhist[0]/4096.;
				double e0=0, e1=0;
				if(p0a>0)
					e0-=hcurr[0]*log2(p0a);
				if(p0a<1)
					e0-=(sum-hcurr[0])*log2(1-p0a);
				e0/=8;
				if(p0b>0)
					e1-=hcurr[0]*log2(p0b);
				if(p0b<1)
					e1-=(sum-hcurr[0])*log2(1-p0b);
				e1/=8;
				printf("C%d CTX%2d  %8.4lf%%->%8.4lf%% zeros  %7d  %12.2lf %12.2lf\n"
					, kc/NCTX
					, kc%NCTX
					, 100.*hcurr[0]/sum
					, currnhist[0]*(100./4096)
					, sum
					, e0
					, e1
				);
			}
#endif
			{//Spread symbols
				uint32_t position=0, symbol;
				for(symbol=0;symbol<=VMAX;++symbol)
				{
					int freq=currnhist[symbol], nbOccurrences;
					for(nbOccurrences=0;nbOccurrences<freq;++nbOccurrences)
					{
						tableSymbol[position]=(FSE_FUNCTION_TYPE)symbol;
						position=(position+STEP)&TABLEMASK;
						while((int)position>highthreshold)
							position=(position+STEP)&TABLEMASK;//Low probability area
					}
				}
				if(position)//Must have initialized all positions
				{
					CRASH("Must have initialized all positions");
					return;
				}
			}

#ifdef __clang_analyzer__
			memset(fse->wksp, 0, sizeof(fse->wksp));//useless initialization, just to keep scan-build happy
#endif
			//For explanations on how to distribute symbol values over the table
			//http://fastcompression.blogspot.fr/2014/02/fse-distributing-symbol-values.html
			{//symbol start positions
				uint32_t u;
				CDF[0]=0;
				for(u=1;u<=VMAX+1;++u)
				{
					if(currnhist[u-1]==-1)//Low probability symbol
					{
						CDF[u]=CDF[u-1]+1;
						tableSymbol[highthreshold--]=(FSE_FUNCTION_TYPE)(u-1);
					}
					else
						CDF[u]=CDF[u-1]+currnhist[u-1];
				}
				CDF[VMAX+1]=TABLESIZE+1;
			}
			{//Build table
				uint32_t u;
				for(u=0;u<TABLESIZE;++u)
				{
					FSE_FUNCTION_TYPE s=tableSymbol[u];//note : static analyzer may not understand tableSymbol is properly initialized
					tableU16[CDF[s]++]=(uint16_t)(TABLESIZE+u);//TableU16 : sorted by symbol order; gives next state value
				}
			}
			{//Build Symbol Transformation Table
				uint32_t total=0, s;
				for(s=0;s<=VMAX;++s)
				{
					switch(currnhist[s])
					{
					case  0://filling nonetheless, for compatibility with FSE_getMaxNbBits()
						symbolTT[s].deltaNbBits=((FSE_PROBBITS+1)<<16)-(1<<FSE_PROBBITS);
						break;
					case -1:
					case  1:
						symbolTT[s].deltaNbBits=(FSE_PROBBITS<<16)-(1<<FSE_PROBBITS);
						symbolTT[s].deltaFindState=total-1;
						++total;
						break;
					default:
						{
							uint32_t const maxBitsOut=FSE_PROBBITS-(31-LZCNT32(currnhist[s]-1));
							uint32_t const minStatePlus=currnhist[s]<<maxBitsOut;
							symbolTT[s].deltaNbBits=(maxBitsOut<<16)-minStatePlus;
							symbolTT[s].deltaFindState=total-currnhist[s];
							total+=currnhist[s];
						}
						break;
					}
				}
			}
#if defined _DEBUG && 0		//debug : symbol costs
			printf("table statistics\n");
			{
				uint32_t symbol;
				for(symbol=0;symbol<=maxSymbolValue;++symbol)
				{
					printf("%3u: w=%3i,   maxBits=%u, fracBits=%.2f"
						, symbol
						, normalizedCounter[symbol]
						, FSE_getMaxNbBits(symbolTT, symbol)
						, (double)FSE_bitCost(symbolTT, tableLog, symbol, 8)/256
					);
				}
			}
#endif
			if(nlevels==1)//RLE
				currnhist[0]=vmax;
		}
		bitpacker_enc_init(&ec, *pstreamptr, streamend);
		for(kc=0;kc<3;++kc)
			states[kc]=1<<FSE_PROBBITS;
		ctxptr=ctxbuf+3*(pixelcount-1);
		for(ky=y2-1;ky>=y1;--ky)
		{
			int16_t *imptr=image+3*(iw*ky+x2-1);
			for(kx=x2-1;kx>=x1;--kx)
			{
				const FSE_symbolCompressionTransform symbolTT2=fse->symbolTT[ctxptr[2]][imptr[2]];
				const FSE_symbolCompressionTransform symbolTT1=fse->symbolTT[ctxptr[1]][imptr[1]];
				const FSE_symbolCompressionTransform symbolTT0=fse->symbolTT[ctxptr[0]][imptr[0]];
				const uint32_t nbBitsOut2=(states[2]+symbolTT2.deltaNbBits)>>16;
				const uint32_t nbBitsOut1=(states[1]+symbolTT1.deltaNbBits)>>16;
				const uint32_t nbBitsOut0=(states[0]+symbolTT0.deltaNbBits)>>16;
				uint64_t sym2=states[2]&g_bitmask[nbBitsOut2];
				uint64_t sym1=states[1]&g_bitmask[nbBitsOut1];
				uint64_t sym0=states[0]&g_bitmask[nbBitsOut0];
				states[2]=(states[2]>>nbBitsOut2)+symbolTT2.deltaFindState;
				states[1]=(states[1]>>nbBitsOut1)+symbolTT1.deltaFindState;
				states[0]=(states[0]>>nbBitsOut0)+symbolTT0.deltaFindState;
				sym2|=sym1<<nbBitsOut2;
				sym2|=sym0<<(nbBitsOut2+nbBitsOut1);
				bitpacker_enc(&ec, (int64_t)nbBitsOut2+nbBitsOut1+nbBitsOut0, sym2);
				states[2]=fse->stateTable[ctxptr[2]][states[2]];
				states[1]=fse->stateTable[ctxptr[1]][states[1]];
				states[0]=fse->stateTable[ctxptr[0]][states[0]];
				imptr-=3;
				ctxptr-=3;
			}
		}
		for(kc=3-1;kc>=0;--kc)
			bitpacker_enc(&ec, FSE_PROBBITS, states[kc]&((1<<FSE_PROBBITS)-1));
		for(kc=3*NCTX-1;kc>=0;--kc)
		{
			int16_t *currnhist=nhists+256*kc;
			if(degenmask>>kc&1)//degenerate distribution
				bitpacker_enc(&ec, 8, (uint8_t)currnhist[0]);
			else if(!(bypassmask>>kc&1))//ordinary distribution
				enc_packhist(&ec, (int*)currnhist, bypassmask, kc, 0, FSE_PROBBITS);
		}
		bitpacker_enc(&ec, (int64_t)3*NCTX, bypassmask);
		bitpacker_enc(&ec, (int64_t)3*NCTX, degenmask);
		{
			uint8_t *p0=*pstreamptr;
			uint8_t *p=bitpacker_enc_flush(&ec);
			intptr_t csize=p-p0;
			*(uint64_t*)p0=csize;
			*pstreamptr=p;
		}
	}
	else
	{
		uint8_t *streamptr=*pstreamptr;
		uint64_t csize=*(uint64_t*)streamptr;

		if(csize>(uint64_t)iw*ih)
			CRASH("subband size %lld", csize);
		bitpacker_dec_init(&ec, streamptr, streamptr+csize);
		degenmask=bitpacker_dec(&ec, (int64_t)3*NCTX);
		bypassmask=bitpacker_dec(&ec, (int64_t)3*NCTX);
		for(kc=0;kc<3*NCTX;++kc)
		{
			int16_t nhist[256]={0};
			if(degenmask>>kc&1)//degenerate distribution
			{
				int sym=(int)bitpacker_dec(&ec, 8);
				memset(nhist, 0, sizeof(nhist));
				nhist[sym]=(1<<FSE_PROBBITS)-1;
				nhist[(uint8_t)(sym+1)]=1;
			}
			else if(bypassmask>>kc&1)//bypass distribution
			{
				for(int ks=0;ks<256;++ks)
					nhist[ks]=1<<FSE_PROBBITS>>8;
			}
			else//ordinary distribution
				dec_unpackhist(&ec, bypassmask, kc, (uint16_t*)nhist, FSE_PROBBITS);

			{
				enum
				{
					TABLESIZE=1<<FSE_PROBBITS,
					TABLEMASK=TABLESIZE-1,
					STEP=FSE_TABLESTEP(TABLESIZE),
					VMAX=255,
					NLEVELS=VMAX+1,
				};
				FSE_DECODE_TYPE *const tableDecode=(FSE_DECODE_TYPE*)fse->decodeTable[kc];
				uint16_t symbolNext[FSE_MAX_SYMBOL_VALUE+1];
				uint32_t highThreshold=TABLESIZE-1;

				if(VMAX>FSE_MAX_SYMBOL_VALUE)//Sanity Checks
				{
					CRASH("maxSymbolValue too large");
					return;
				}
				if(FSE_PROBBITS>FSE_MAX_TABLELOG)
				{
					CRASH("tableLog too large");
					return;
				}
				{//Init, lay down lowprob symbols
					uint32_t s;
					for(s=0;s<NLEVELS;++s)
					{
						if(nhist[s]==-1)
						{
							tableDecode[highThreshold--].symbol=(FSE_FUNCTION_TYPE)s;
							symbolNext[s]=1;
						}
						else
							symbolNext[s]=nhist[s];
					}
				}
				{//Spread symbols
					uint32_t s, position=0;
					for(s=0;s<NLEVELS;++s)
					{
						int i;
						for(i=0;i<nhist[s];++i)
						{
							tableDecode[position].symbol=(FSE_FUNCTION_TYPE)s;
							position=(position+STEP)&TABLEMASK;
							while(position>highThreshold)
								position=(position+STEP)&TABLEMASK;//lowprob area
						}
					}
					if(position!=0)//position must reach all cells once, otherwise normalizedCounter is incorrect
					{
						CRASH("position must reach all cells once, otherwise normalizedCounter is incorrect");
						return;
					}
				}
				{//Build Decoding table
					uint32_t u;
					for(u=0;u<TABLESIZE;++u)
					{
						FSE_FUNCTION_TYPE const symbol=(FSE_FUNCTION_TYPE)tableDecode[u].symbol;
						uint32_t const nextState=symbolNext[symbol]++;
						tableDecode[u].nbBits=(uint8_t)(FSE_PROBBITS-31+LZCNT32(nextState));
						tableDecode[u].newState=(uint16_t)((nextState<<tableDecode[u].nbBits)-TABLESIZE);
					}
				}
			}
		}
		for(kc=0;kc<3;++kc)
			states[kc]=(uint32_t)bitpacker_dec(&ec, FSE_PROBBITS);
		for(ky=y1;ky<y2;++ky)
		{
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			};
			int16_t *imptr=image+3*(iw*ky+x1);
			for(kx=x1;kx<x2;++kx)
			{
				int
					eNW0	=rows[1][0-1*NCH*NROWS*NVAL],
					eNW1	=rows[1][1-1*NCH*NROWS*NVAL],
					eNW2	=rows[1][2-1*NCH*NROWS*NVAL],
					eN0	=rows[1][0+0*NCH*NROWS*NVAL],
					eN1	=rows[1][1+0*NCH*NROWS*NVAL],
					eN2	=rows[1][2+0*NCH*NROWS*NVAL],
					eNE0	=rows[1][0+1*NCH*NROWS*NVAL],
					eNE1	=rows[1][1+1*NCH*NROWS*NVAL],
					eNE2	=rows[1][2+1*NCH*NROWS*NVAL],
					eW0	=rows[0][0-1*NCH*NROWS*NVAL],
					eW1	=rows[0][1-1*NCH*NROWS*NVAL],
					eW2	=rows[0][2-1*NCH*NROWS*NVAL];
				int ctx0=2*((eW0|eN0)!=0)+((eNW0|eNE0)!=0);
				int ctx1=2*((eW1|eN1)!=0)+((eNW1|eNE1)!=0);
				int ctx2=2*((eW2|eN2)!=0)+((eNW2|eNE2)!=0);
				//int ctx0=4*(eN0!=0)+2*(eW0!=0)+((eNW0|eNE0)!=0);
				//int ctx1=4*(eN1!=0)+2*(eW1!=0)+((eNW1|eNE1)!=0);
				//int ctx2=4*(eN2!=0)+2*(eW2!=0)+((eNW2|eNE2)!=0);
				//int ctx0=(eW0|eN0)!=0;
				//int ctx1=(eW1|eN1)!=0;
				//int ctx2=(eW2|eN2)!=0;
				//int ctx0=31^LZCNT32(eW0+1);
				//int ctx1=31^LZCNT32(eW1+1);
				//int ctx2=31^LZCNT32(eW2+1);
				if(ctx0>NCTX-1)ctx0=NCTX-1;
				if(ctx1>NCTX-1)ctx1=NCTX-1;
				if(ctx2>NCTX-1)ctx2=NCTX-1;
				ctx0+=0*NCTX;
				ctx1+=1*NCTX;
				ctx2+=2*NCTX;
				uint32_t DInfo0=*(uint32_t*)&fse->decodeTable[ctx0][states[0]];
				uint32_t DInfo1=*(uint32_t*)&fse->decodeTable[ctx1][states[1]];
				uint32_t DInfo2=*(uint32_t*)&fse->decodeTable[ctx2][states[2]];
				int sym0=imptr[0]=(uint8_t)(DInfo0>>16);
				int sym1=imptr[1]=(uint8_t)(DInfo1>>16);
				int sym2=imptr[2]=(uint8_t)(DInfo2>>16);
				uint64_t n0=(uint8_t)(DInfo0>>24);
				uint64_t n1=(uint8_t)(DInfo1>>24);
				uint64_t n2=(uint8_t)(DInfo2>>24);
				states[0]=(uint16_t)DInfo0;
				states[1]=(uint16_t)DInfo1;
				states[2]=(uint16_t)DInfo2;
				int64_t c1=n2;
				int64_t c2=n2+n1;
				size_t lowBits=bitpacker_dec(&ec, n0+n1+n2);
#ifdef LOUD
				if(ec.ptr<ec.start)
					CRASH("stream ended at %12.6lf%% of subband"
						, 100.*((double)(x2-x1)*(ky-y1)+kx-x1)/((double)(x2-x1)*(y2-y1))
					);
#endif
				states[0]+=(uint32_t)(lowBits>>c2&g_bitmask[n0]);
				states[1]+=(uint32_t)(lowBits>>c1&g_bitmask[n1]);
				states[2]+=(uint32_t)(lowBits    &g_bitmask[n2]);
				imptr+=3;
				rows[0][0+0*NCH*NROWS*NVAL]=sym0;
				rows[0][1+0*NCH*NROWS*NVAL]=sym1;
				rows[0][2+0*NCH*NROWS*NVAL]=sym2;
				//{
				//	int
				//		eNEE0	=rows[1][0+2*NCH*NROWS*NVAL],
				//		eNEE1	=rows[1][1+2*NCH*NROWS*NVAL],
				//		eNEE2	=rows[1][2+2*NCH*NROWS*NVAL],
				//		eNEEE0	=rows[1][0+3*NCH*NROWS*NVAL],
				//		eNEEE1	=rows[1][1+3*NCH*NROWS*NVAL],
				//		eNEEE2	=rows[1][2+3*NCH*NROWS*NVAL];
				//	rows[0][0+0*NCH*NROWS*NVAL]=(2*eW0+(sym0<<RBITS)+(eNEE0>eNEEE0?eNEE0:eNEEE0))>>2;
				//	rows[0][1+0*NCH*NROWS*NVAL]=(2*eW1+(sym1<<RBITS)+(eNEE1>eNEEE1?eNEE1:eNEEE1))>>2;
				//	rows[0][2+0*NCH*NROWS*NVAL]=(2*eW2+(sym2<<RBITS)+(eNEE2>eNEEE2?eNEE2:eNEEE2))>>2;
				//}
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				rows[2]+=NCH*NROWS*NVAL;
				rows[3]+=NCH*NROWS*NVAL;
			}
		}
		*pstreamptr+=csize;
	}
	free(fse);
}
static void subband_code_ac(int16_t *image
	, int iw, int ih
	, int x1, int x2, int y1, int y2
	, int16_t *pixels, int psize
	, int fwd
	, uint8_t **pstreamptr, uint8_t *streamend
)
{
	enum
	{
		RBITS=2,
		NLEVELS=256,
	};
	int kx=0, ky=0, kc=0;
	int sym=0;
	uint64_t lo=0, hi=0xFFFFFFFFFFFF, code=0, x=0;
	int den=0, cdf=0, freq=0, tmp=0;
	uint8_t *streamptr=*pstreamptr;

	memset(hists, 0, sizeof(hists));
	//for(int ctx=0;ctx<3*NCTX;++ctx)
	//{
	//	int32_t *currhist=hists+NLEVELS*ctx;
	//	int total=0x8000, sum=0;
	//	for(int ks=0;ks<NLEVELS-1;++ks)
	//		sum+=currhist[ks]=total>>=1;
	//	currhist[NLEVELS-1]=sum;
	//}
	memset(pixels, 0, psize);
	if(!fwd)
	{
		code=*(uint64_t*)streamptr;//load
		code=code<<32|code>>32;
		streamptr+=sizeof(uint64_t);
	}
	for(ky=y1;ky<y2;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		int16_t *imptr=image+3*(iw*ky+x1);
		for(kx=x1;kx<x2;++kx)
		{
			for(kc=0;kc<3;++kc)
			{
				int
					eNW	=rows[1][0-1*NCH*NROWS*NVAL],
					eN	=rows[1][0+0*NCH*NROWS*NVAL],
					eNE	=rows[1][0+1*NCH*NROWS*NVAL],
				//	eNEE	=rows[1][0+2*NCH*NROWS*NVAL],
				//	eNEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					eW	=rows[0][0-1*NCH*NROWS*NVAL];
				int32_t *currhist;
				int ctx;

				ctx=2*((eW|eN)!=0)+((eNW|eNE)!=0);
				//ctx=31^LZCNT32(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
				ctx+=NCTX*kc;
				currhist=hists+NLEVELS*ctx;
				den=currhist[NLEVELS-1]+NLEVELS;
				if(fwd)
				{
					sym=*imptr;
					//if(sym==255)
					//	CRASH("");
					x=hi-lo;
					if(x<=0xFFFF)
					{
						*(uint32_t*)streamptr=(uint32_t)(lo>>32);
						streamptr+=sizeof(uint32_t);
						lo<<=32;
						hi=hi<<32|~0U;
						if(hi<lo)
							hi=~0ULL;
						x=hi-lo;
					}
					for(tmp=0, cdf=0;;++tmp)
					{
						freq=currhist[tmp]+1;
						if(tmp>=sym)
							break;
						cdf+=freq;
					}
#ifdef _MSC_VER
					if(!freq||(uint32_t)freq>(uint32_t)den||(uint32_t)cdf>(uint32_t)den)
						CRASH("");
#endif
#ifdef FIFOVAL
					fifoval_enqueue(freq<<16|cdf);
#endif
					lo+=x*cdf/den;
					hi=lo+x*freq/den-1;
				}
				else
				{
					x=hi-lo;
					if(x<=0xFFFF)//stall: unpredictable branch
					{
#ifdef _MSC_VER
						if(streamptr>streamend)
							CRASH("");
#endif
						code=code<<32|*(uint32_t*)streamptr;
						streamptr+=sizeof(uint32_t);
						lo<<=32;
						hi=hi<<32|~0U;
						if(hi<lo)
							hi=~0ULL;
						x=hi-lo;
					}
					tmp=(int)(((code-lo+1)*den-1)/x);
					for(sym=0, cdf=0;;)
					{
						freq=currhist[sym]+1;
						if(cdf+freq>tmp)
							break;
						cdf+=freq;
						++sym;
					}
#ifdef FIFOVAL
					fifoval_check(freq<<16|cdf);
#endif
					lo+=x*cdf/den;
					hi=lo+x*freq/den-1;

					*imptr=sym;
				}
				++currhist[sym];
				++currhist[NLEVELS-1];
				if(currhist[NLEVELS-1]>=0x8000)//rescale
				{
					den=0;
					for(int k=0;k<NLEVELS-1;++k)
						den+=currhist[k]>>=1;
					currhist[NLEVELS-1]=den;
				}
				rows[0][0]=sym;
				//rows[0][0]=(2*eW+(sym<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
		}
	}
	if(fwd)
	{
		*(uint64_t*)streamptr=lo<<32|lo>>32;//flush
		streamptr+=sizeof(uint64_t);
	}
	//*(uint32_t*)*pstreamptr=streamptr-**pstreamptr;
	*pstreamptr=streamptr;
}
typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end;
} ACState;
INLINE void codebit(ACState *ac, uint32_t *pcell, int32_t *pbit, const int fwd)
{
	enum
	{
		RICE_USEBITS=14,
		RICE_STOREBITS=20,
		RICE_CTRBITS=9,
	};
	uint64_t x;
	int32_t cell, prob, count, p1, sh;
	int bit;

	cell=*pcell;
	x=ac->hi-ac->lo;
	prob=(int32_t)cell>>RICE_CTRBITS;
	count=cell&((1<<RICE_CTRBITS)-1);
	p1=(int32_t)cell>>(RICE_STOREBITS+RICE_CTRBITS-RICE_USEBITS);
	sh=31^LZCNT32(count+2);
	p1+=p1<0;
	count+=count<(1<<RICE_CTRBITS)-1;
	p1+=1<<RICE_USEBITS>>1;
	if(x<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
			CRASH("inflation");
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->lo>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=sizeof(uint32_t);
		ac->lo<<=32;
		ac->hi=ac->hi<<32|~0U;
		if(ac->hi<ac->lo)ac->hi=~0ULL;
		x=ac->hi-ac->lo;
	}
#ifdef _MSC_VER
	if((uint32_t)(p1-1)>(uint32_t)((1<<RICE_USEBITS)-2))
		CRASH("Invalid p1 0x%08X / %d bit", p1, RICE_USEBITS);
#endif
	x=ac->lo+(x*(uint32_t)p1>>RICE_USEBITS);
	bit=*pbit;
	bit=fwd?bit:ac->code<x;
	*pbit=bit;
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue(bit<<24^p1);
	else
		fifoval_check(bit<<24^p1);
#endif
	*(bit?&ac->hi:&ac->lo)=x-bit;
	prob+=(int32_t)((bit<<RICE_STOREBITS)-(1<<RICE_STOREBITS>>1)-prob+(1<<sh>>1))>>sh;
	*pcell=prob<<RICE_CTRBITS|count;
}
static void subband_code_riceac(int16_t *image
	, int iw, int ih
	, int x1, int x2, int y1, int y2
	, int16_t *pixels, int psize
	, int fwd, int LL
	, uint8_t **pstreamptr, uint8_t *streamend
)
{
	enum
	{
		RBITS=2,
		NLEVELS=256,
	};
	int kx=0, ky=0, kc=0;
	int sym=0;
	int nzeros=0, tidx=0;
	ACState ac={0};

	memset(hists, 0, sizeof(hists));
	memset(pixels, 0, psize);
	ac.hi=0xFFFFFFFFFFFF;
	ac.ptr=*pstreamptr;
	ac.end=streamend;
	if(!fwd)
	{
		ac.code=*(uint64_t*)ac.ptr;//load
		ac.code=ac.code<<32|ac.code>>32;
		ac.ptr+=sizeof(uint64_t);
	}
	for(ky=y1;ky<y2;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		int16_t *imptr=image+3*(iw*ky+x1);
		for(kx=x1;kx<x2;++kx)
		{
			for(kc=0;kc<3;++kc)
			{
				int
					eNW	=rows[1][0-1*NCH*NROWS*NVAL],
					eN	=rows[1][0+0*NCH*NROWS*NVAL],
					eNE	=rows[1][0+1*NCH*NROWS*NVAL],
				//	eNEE	=rows[1][0+2*NCH*NROWS*NVAL],
				//	eNEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					eW	=rows[0][0-1*NCH*NROWS*NVAL];
				int32_t *currhist;
				int ctx, nbypass;

				nbypass=31^LZCNT32(eW+1);
				ctx=LL?nbypass:2*((eW|eN)!=0)+((eNW|eNE)!=0);
				//ctx=31^LZCNT32(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
				ctx+=NCTX*kc;
				currhist=hists+NLEVELS*ctx;

				sym=0;
				nzeros=0;
				if(fwd)
				{
					sym=*imptr++;
					nzeros=sym>>nbypass;
				}
				tidx=0;
				for(;;)
				{
					int bit=tidx>=nzeros;
					int c=tidx;
					if(c>255)c=255;
					codebit(&ac, (uint32_t*)currhist+c, &bit, fwd);
					if(bit)
						break;
					++tidx;
				}
				for(int kb=nbypass-1;kb>=0;--kb)
				{
					int bit=sym>>kb&1;
					codebit(&ac, (uint32_t*)currhist+255-kb, &bit, fwd);
					sym|=bit<<kb;
				}
				if(!fwd)
				{
					sym|=tidx<<nbypass;
					*imptr++=sym;
				}
				rows[0][0]=sym;
				//rows[0][0]=(2*eW+(sym<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
		}
	}
	if(fwd)
	{
		*(uint64_t*)ac.ptr=ac.lo<<32|ac.lo>>32;//flush
		ac.ptr+=sizeof(uint64_t);
	}
	//*(uint32_t*)*pstreamptr=streamptr-**pstreamptr;
	*pstreamptr=ac.ptr;
}
static void subband_proc(int16_t *image
	, int iw, int ih
	, int x1, int x2, int y1, int y2
	, int16_t *pixels, int psize
	, int niter, int level, int LL, int fwd, int dist
	, uint8_t *ctxbuf
	, uint8_t **pstreamptr, uint8_t *streamend
)
{
#define UPSCALE16(X) (int)((double)(X)*(1<<16)+0.5)
	static const int coeffs[]=
	{
#if 1
		UPSCALE16(1.0),
		UPSCALE16(1.0),
		UPSCALE16(1.1),
		UPSCALE16(1.2),
		UPSCALE16(1.4),
		UPSCALE16(1.8),
	//	UPSCALE16(2.0),
	//	UPSCALE16(4.0),
	//	UPSCALE16(8.0),
#endif
#if 0
		UPSCALE16(1.0),
		UPSCALE16(1.1),
		UPSCALE16(1.2),
		UPSCALE16(1.3),
		UPSCALE16(1.4),
		UPSCALE16(1.5),
		UPSCALE16(1.6),
		UPSCALE16(1.7),
#endif
#if 0
		UPSCALE16(1.7),
		UPSCALE16(1.6),
		UPSCALE16(1.5),
		UPSCALE16(1.4),
		UPSCALE16(1.3),
		UPSCALE16(1.2),
		UPSCALE16(1.1),
		UPSCALE16(1.0),
#endif
		//UPSCALE16(1),
		//UPSCALE16(1),
		//UPSCALE16(0.75),
		//UPSCALE16(0.75),
		//UPSCALE16(0.5),
		//UPSCALE16(0.5),
		//UPSCALE16(0.25),
		//UPSCALE16(0.25),

		//0x08000,
		//0x10000,
		//0x20000,
		//0x40000,
		//0x80000,
	};
	int d2=level, d3=level-1;
	if(d2>_countof(coeffs)-1)
		d2=_countof(coeffs)-1;
	if(d3<0)
		d3=0;
	if(d3>_countof(coeffs)-1)
		d3=_countof(coeffs)-1;
	int qsteps[]=
	{
		dist*coeffs[d2],//y
		dist*coeffs[d3],//u
		dist*coeffs[d3],//v
		//dist<<16,
		//dist<<16,
		//dist<<16,
		//(dist-level*8)<<16,
		//(dist-level*8)<<16,
		//(dist-level*8)<<16,
		//dist<<(16-level),
		//dist<<(16-level),
		//dist<<(16-level),
		//(int)round((double)dist*pow(1.25, -level*0.5)*(1<<16)),
		//(int)round((double)dist*pow(1.25, -level*0.5)*(1<<16)),
		//(int)round((double)dist*pow(1.25, -level*0.5)*(1<<16)),
	};
	uint8_t *streamptr=*pstreamptr;

	if(qsteps[0]<1<<16)qsteps[0]=1<<16;
	if(qsteps[1]<1<<16)qsteps[1]=1<<16;
	if(qsteps[2]<1<<16)qsteps[2]=1<<16;
	if(fwd)
	{
#if defined _MSC_VER && 0
		printf("%12.6lf %12.6lf %12.6lf  level %d\n"
			, qsteps[0]*(1./(1<<16))
			, qsteps[1]*(1./(1<<16))
			, qsteps[2]*(1./(1<<16))
			, level
		);
#endif
		if(LL)
			predictLL(image, iw, ih, x2-x1, y2-y1, pixels, psize, DIST_LL, fwd);
		else
		{
			//getqsteps(image, iw, ih, x1, x2, y1, y2, dist, qsteps);
			quantization(image, iw, ih, x1, x2, y1, y2, qsteps, fwd);
			//memcpy(streamptr, qsteps, sizeof(qsteps));
			//streamptr+=sizeof(qsteps);
		}
	}
	else
	{
		if(!LL)
		{
			//memcpy(qsteps, streamptr, sizeof(qsteps));
			//streamptr+=sizeof(qsteps);
		}
	}
	subband_code_riceac	(image, iw, ih, x1, x2, y1, y2, pixels, psize, fwd, LL, &streamptr, streamend);
//	subband_code_ac		(image, iw, ih, x1, x2, y1, y2, pixels, psize, fwd, &streamptr, streamend);
//	subband_code_fse	(image, iw, ih, x1, x2, y1, y2, pixels, psize, fwd, ctxbuf, &streamptr, streamend);
	(void)&subband_code_riceac;
	(void)&subband_code_ac;
	(void)&subband_code_fse;
	if(!fwd)
	{
		if(LL)
			predictLL(image, iw, ih, x2-x1, y2-y1, pixels, psize, DIST_LL, fwd);
		else
			quantization(image, iw, ih, x1, x2, y1, y2, qsteps, fwd);
	}
	*pstreamptr=streamptr;
}
int c60_codec(int argc, char **argv)
{
	const uint16_t tag='6'|'0'<<8;
	
	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0, dist=0;
	int64_t usize=0, ccap=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint8_t *image=0, *stream=0, *streamptr=0, *streamend=0;
	int ict=0;
	int niter=0;
	uint8_t *ctxbuf=0;
#ifdef LOUD
	double t=time_sec();
#endif

	if((uint32_t)(argc-3)>(uint32_t)(4-3))
	{
		printf(
			"Usage:  \"%s\"  src  dst  [d]\n"
			"  d: Distortion. Always lossy. Default %d.\n"
			"Only for 24-bit PPM images\n"
			"Built on %s %s\n"
			, argv[0]
			, DIST_DEFAULT
			, __DATE__
			, __TIME__
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	dist=argc<4?DIST_DEFAULT:atoi(argv[3]);
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	fread(&c, 1, 2, fsrc);
	fwd=c==('P'|'6'<<8);
	if(!fwd&&c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(fwd)//parse header
	{
		c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		c=fgetc(fsrc);
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		iw=0;
		while((uint32_t)(c-'0')<10)
		{
			iw=10*iw+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		ih=0;
		while((uint32_t)(c-'0')<10)
		{
			ih=10*ih+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		c|=(int64_t)fgetc(fsrc)<<8*1;
		c|=(int64_t)fgetc(fsrc)<<8*2;
		c|=(int64_t)fgetc(fsrc)<<8*3;
		c|=(int64_t)fgetc(fsrc)<<8*4;
		if(c!=(
			(uint64_t)'\n'<<8*0|
			(uint64_t) '2'<<8*1|
			(uint64_t) '5'<<8*2|
			(uint64_t) '5'<<8*3|
			(uint64_t)'\n'<<8*4
		))
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		ccap=(int64_t)4*iw*ih;
	}
	else
	{
		iw=0;
		ih=0;
		ict=0;
		dist=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
	//	fread(&ict, 1, 1, fsrc);
		fread(&dist, 1, 2, fsrc);
		{
			struct stat info={0};

			stat(srcfn, &info);
			csize=info.st_size;
			ccap=csize-ftell(fsrc);
		}
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	image=(uint8_t*)malloc(usize*2);
	stream=(uint8_t*)malloc(ccap);
	psize=(iw+2*XPAD)*(int)sizeof(int32_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!image||!stream||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	streamptr=stream;
	streamend=stream+ccap;
	if(fwd)
	{
		c=fread(image, 1, usize, fsrc);
#ifdef ENABLE_GUIDE
#ifdef C60TEST2
		if(!g_image)
#endif
			guide_save(image, iw, ih);
#endif
		ict=0;
		//ict=ict_select(image, iw, ih);
#if defined LOUD && 0
		printf("PERM %d: %d %d %d\n"
			, ict
			, ict_perm[ict*3+0]
			, ict_perm[ict*3+1]
			, ict_perm[ict*3+2]
		);
		//exit(0);
#endif
	}
	else
		c=fread(stream, 1, ccap, fsrc);
	fclose(fsrc);

	niter=dwt_getNiter(iw, ih, dist);
	if(fwd)
	{
		//int64_t ctxbufsize=(int64_t)iw*ih;
		//
		//if(niter<=1)
		//	ctxbufsize=3*ctxbufsize;
		ict_fwd(image, iw, ih, ict);
		dwt_cdf97((int16_t*)image, iw, ih, niter, dist, fwd);
#ifdef _MSC_VER
	//	ppm_save_16as8("20260808Sa_0616pm.ppm", (int16_t*)image, iw, ih);//
	//	ppm_save_16bit("20260808Sa_0616pm.ppm", (int16_t*)image, iw, ih, 0xFFFF);//
#endif
		//ctxbuf=(uint8_t*)malloc(ctxbufsize);
		//if(!ctxbuf)
		//{
		//	CRASH("Alloc error");
		//	return 1;
		//}
	//	memset(hists, 0, sizeof(hists));
		subband_proc((int16_t*)image, iw, ih, 0, iw>>(niter-1), 0, ih>>(niter-1), pixels, psize, 2*niter, 2*(niter-1), 1, fwd, dist, ctxbuf, &streamptr, streamend);//LL
		for(int it=0;it<niter-1;)
		{
			int w2, h2;//2x
			int w1, h1;

			w2=iw>>it;
			h2=ih>>it;
			++it;
			w1=iw>>it;
			h1=ih>>it;
		//	memset(hists, 0, sizeof(hists));
			subband_proc((int16_t*)image, iw, ih,  0, w1, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SW
			subband_proc((int16_t*)image, iw, ih, w1, w2,  0, h1, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//NE
			subband_proc((int16_t*)image, iw, ih, w1, w2, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+0, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SE HH
		}
		//free(ctxbuf);
	}
	else
	{
	//	memset(hists, 0, sizeof(hists));
		subband_proc((int16_t*)image, iw, ih, 0, iw>>(niter-1), 0, ih>>(niter-1), pixels, psize, 2*niter, 2*(niter-1), 1, fwd, dist, ctxbuf, &streamptr, streamend);//LL
		for(int it=0;it<niter-1;)
		{
			int w2, h2;//2x
			int w1, h1;

			w2=iw>>it;
			h2=ih>>it;
			++it;
			w1=iw>>it;
			h1=ih>>it;
		//	memset(hists, 0, sizeof(hists));
			subband_proc((int16_t*)image, iw, ih,  0, w1, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SW
			subband_proc((int16_t*)image, iw, ih, w1, w2,  0, h1, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//NE
			subband_proc((int16_t*)image, iw, ih, w1, w2, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+0, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SE HH
		}
		dwt_cdf97((int16_t*)image, iw, ih, niter, dist, fwd);
		ict_inv(image, iw, ih, ict);
	}
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
	//	csize+=fwrite(&ict, 1, 1, fdst);
		csize+=fwrite(&dist, 1, 2, fdst);
		csize+=fwrite(stream, 1, streamptr-stream, fdst);
	}
	else
	{
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(image, 1, usize, fdst);
	}
	fclose(fdst);
	free(pixels);
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
#if !defined C60TEST2
	if(fwd)
	{
		printf("WH %5d*%5d  %8.2lf MB  \"%s\"\n", iw, ih, 3.*iw*ih/(1<<20), srcfn);
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
			, 8.*csize/usize
		);
#ifdef C60TEST
		results[dist].bps=(double)csize/usize*8;
#endif
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
#endif
	if(!fwd)
	{
		static const double jxlcurve[]=
		{//	BPS	PSNR

			//DIV2K
#if 0
			0.0677, 26.617987,//10
			0.0902, 28.435101,//20
			0.1202, 29.648848,//30
			0.1353, 30.156079,//40
			0.1567, 30.831602,//50
			0.1908, 31.723736,//60
			0.2463, 32.956886,//70
			0.2850, 33.663252,//75
			0.3414, 34.609978,//80

#endif

			//DIV2K/0820
#if 1
			0.113, 25.484142,//10	0.340/3
			0.139, 27.991428,//20	0.419/3
			0.184, 29.615995,//30	0.552/3
			0.205, 30.297540,//40	0.617/3
			0.234, 31.133699,//50	0.704/3
			0.278, 32.215792,//60	0.835/3
			0.346, 33.726246,//70	1.040/3
			0.391, 34.617069,//75	1.174/3
			0.453, 35.742132,//80	1.361/3
			0.688, 38.923345,//90	2.065/3
			0.985, 41.586402,//95	2.957/3
			1.761, 45.738214,//99	5.283/3
			3.197, INFINITY,//100	9.592/3
#endif
		};
		extern double g_score;
		double bps=8.*csize/usize;
		double psnr=guide_psnr(1/((double)iw*ih));
		double exbps=0;
		double score=0;
		int bestidx=0;

		for(int k=0;k<_countof(jxlcurve)/2-2;++k)
		{
			double dx=jxlcurve[2*(k+1)+0]-jxlcurve[2*(k+0)+0];
			double dy=jxlcurve[2*(k+1)+1]-jxlcurve[2*(k+0)+1];
			double sd=(dx*(psnr-jxlcurve[2*(k+0)+1])-dy*(bps-jxlcurve[2*(k+0)+0]))/sqrt(dx*dx+dy*dy);
			if(!score||score<sd)
				score=sd, bestidx=k;
			//if(psnr<jxlcurve[2*(k+1)+1])
			//{
			//	exbps=(psnr-jxlcurve[2*(k+0)+1])
			//		/(jxlcurve[2*(k+1)+1]-jxlcurve[2*(k+0)+1])
			//		*(jxlcurve[2*(k+1)+0]-jxlcurve[2*(k+0)+0])
			//		+jxlcurve[2*(k+0)+0];
			//	break;
			//}
		}
	//	double exbps=5*pow(rmse, -1.75);//expected BPS from JXL curve

#ifdef C60TEST2
		printf("dist %3d  BPS %12.6lf  PSNR %12.6lf  score %12.6lf", dist, bps, psnr, score);
		{
			enum
			{
				WIDTH=256,
			};
			int nstars=(int)round(score*WIDTH);
			if(nstars<0)
			{
				for(int k=-WIDTH/4;k<nstars;++k)
					printf(" ");
				for(int k=nstars;k<0;++k)
					printf("-");
			}
			else
			{
				for(int k=-WIDTH/4;k<0;++k)
					printf(" ");
				for(int k=0;k<nstars;++k)
					printf("+");
			}
		}
		printf("\n");
		g_score+=score;
#endif
		//printf(
		//	"JXL1      %12.6lf  %12.6lf\n"
		//	"BPS PSNR  %12.6lf  %12.6lf\n"
		//	"JXL2      %12.6lf  %12.6lf\n"
		//	"sd = %12.8lf\n"
		//	, jxlcurve[2*(bestidx+0)+0], jxlcurve[2*(bestidx+0)+1]
		//	, bps, psnr
		//	, jxlcurve[2*(bestidx+1)+0], jxlcurve[2*(bestidx+1)+1]
		//	, dmax
		//);
	//	printf("RMSE %12.6lf  BPS %12.6lf  exBPS %12.6lf  BPS/exBPS %12.6lf%%\n", psnr, bps, exbps, 100.*bps/exbps);
#ifdef C60TEST
		results[dist].psnr=psnr;
#endif
#ifndef C60TEST2
		free(g_image);
#endif
	}
#endif
	(void)csize;
	(void)&ppm_save;
	(void)&ppm_save_16bit;
	(void)&ppm_save_16as8;
	(void)&time_sec;
	(void)&ict_select;
	(void)&getqsteps;
	return 0;
}
#ifdef C60TEST
void c60_test(const char *argv0, const char *fn, const char *tmpfn, const char *tmpfn2)
{
	char arg3[64]={0};
	const char *argse[]=
	{
		argv0,
		fn,
		tmpfn,
		arg3,
	};
	const char *argsd[]=
	{
		argv0,
		tmpfn,
		tmpfn2,
	};

	for(int d=D1;d<=D2;++d)
	{
		snprintf(arg3, _countof(arg3)-1, "%d", d);
		c60_codec(_countof(argse), (char**)argse);
		c60_codec(_countof(argsd), (char**)argsd);
	}
	for(int d=D1;d<=D2;++d)
	{
		ResultC60 *res=results+d;

		printf("%3d  %12.6lf  %12.6lf\n", d, res->bps, res->psnr);
	}
}
#endif
