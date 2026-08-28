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

//	#define C60TEST2
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
		imptr2[0]=c0-128;
		imptr2[1]=c1-128;
		imptr2[2]=c2-128;
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
		c0=imptr2[0]+128;
		c1=imptr2[1]+128;
		c2=imptr2[2]+128;
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

static int32_t hists[3*NCTX*256];

typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end;
} ACState;
INLINE void codebit(ACState *ac, uint32_t *pcell, int32_t *pbit, const int fwd)
{
	enum
	{//STOREBITS+CTRBITS+HISTBITS <= 32
		CTRBITS=9,
		HISTBITS=4,
		STOREBITS=32-CTRBITS-HISTBITS,
		USEBITS=13,
	};
	uint64_t x;
	int32_t cell, prob, count, p1, sh, bias;
	int bit, prev4;

	cell=*pcell;
	x=ac->hi-ac->lo;
	prev4=cell&((1<<HISTBITS)-1);
	cell>>=HISTBITS;
	prob=(int32_t)cell>>CTRBITS;
	count=cell&((1<<CTRBITS)-1);
	p1=(int32_t)cell>>(STOREBITS-USEBITS+CTRBITS);
	sh=31^LZCNT32(count+3);
	p1+=p1<0;
	count+=count<(1<<CTRBITS)-1;
	p1+=1<<USEBITS>>1;
	bias=(1<<sh>>1)-(1<<STOREBITS>>1)-prob;
	/*
	p1 = prev4/32 + p1*31/32
	p1 = prev3/32 + p1*31/32
	p1 = prev2/32 + p1*31/32
	p1 = prev1/32 + p1*31/32

	p1 = prev1/32 + (prev2/32 + (prev3/32 + (prev4/32 + p1*31/32)*31/32)*31/32)*31/32
	p1 = prev1/32 + prev2*31/32^2 + prev3*31^2/32^3 + prev4*31^3/32^4 + p1*(31/32)^4
*/
#if 0
	{
#define MAKEFRAC(N, D, S) (int)((((int64_t)(N)<<(S))+(D)-1)/(D))
#define MAKECOEFF(X) (X>>3&1)*C1+(X>>2&1)*C2+(X>>1&1)*C3+(X>>0&1)*C4	//prev[1234]
		enum
		{
			C1=MAKEFRAC(32, 1, USEBITS),
			C2=MAKEFRAC(31, 1, USEBITS),
			C3=MAKEFRAC(31*31, 32, USEBITS),
			C4=MAKEFRAC(31*31*31, 32*32, USEBITS),
			CP=MAKEFRAC(31*31*31*31, 32*32, 0),
			A0=MAKECOEFF(0x0),
			A1=MAKECOEFF(0x1),
			A2=MAKECOEFF(0x2),
			A3=MAKECOEFF(0x3),
			A4=MAKECOEFF(0x4),
			A5=MAKECOEFF(0x5),
			A6=MAKECOEFF(0x6),
			A7=MAKECOEFF(0x7),
			A8=MAKECOEFF(0x8),
			A9=MAKECOEFF(0x9),
			AA=MAKECOEFF(0xA),
			AB=MAKECOEFF(0xB),
			AC=MAKECOEFF(0xC),
			AD=MAKECOEFF(0xD),
			AE=MAKECOEFF(0xE),
			AF=MAKECOEFF(0xF),
		};
		static const int coeffs[]=
		{
			A0,
			A1,
			A2,
			A3,
			A4,
			A5,
			A6,
			A7,
			A8,
			A9,
			AA,
			AB,
			AC,
			AD,
			AE,
			AF,
#undef  MAKECOEFF
		};
		p1=(coeffs[prev4]+CP*p1)>>10;
	}
#endif
#if 1
	{
		enum
		{
#define MAKEFRAC(N, D, S) (int)((((int64_t)(N)<<(S))+(D)-1)/(D))
			C1=MAKEFRAC(32, 1, USEBITS),
			C2=MAKEFRAC(31, 1, USEBITS),
			C3=MAKEFRAC(31*31, 32, USEBITS),
			C4=MAKEFRAC(31*31*31, 32*32, USEBITS),
			CP=MAKEFRAC(31*31*31*31, 32*32, 0),
#undef  MAKEFRAC
		};
		static const int coeffs[]=
		{//prev[1234]
#define MAKECOEFF(X) (X>>3&1)*C4+(X>>2&1)*C3+(X>>1&1)*C2+(X>>0&1)*C1
			MAKECOEFF(0x0),
			MAKECOEFF(0x1),
			MAKECOEFF(0x2),
			MAKECOEFF(0x3),
			MAKECOEFF(0x4),
			MAKECOEFF(0x5),
			MAKECOEFF(0x6),
			MAKECOEFF(0x7),
			MAKECOEFF(0x8),
			MAKECOEFF(0x9),
			MAKECOEFF(0xA),
			MAKECOEFF(0xB),
			MAKECOEFF(0xC),
			MAKECOEFF(0xD),
			MAKECOEFF(0xE),
			MAKECOEFF(0xF),
#undef  MAKECOEFF
		};
		p1=(coeffs[prev4]+CP*p1)>>10;
		prev4=prev4<<1&((1<<HISTBITS)-1);
	}
#endif
#if 0
	{
#define MAKEFRAC(N, D, S) (int)((((int64_t)(N)<<(S))+(D)-1)/(D))
		enum
		{
			C1=MAKEFRAC(32, 1, USEBITS),
			C2=MAKEFRAC(31, 1, USEBITS),
			C3=MAKEFRAC(31*31, 32, USEBITS),
			C4=MAKEFRAC(31*31*31, 32*32, USEBITS),
			CP=MAKEFRAC(31*31*31*31, 32*32, 0),
		};
		p1=(
			+(-prev1&C1)
			+(-prev2&C2)
			+(-prev3&C3)
			+(-prev4&C4)
			+p1*CP
		)>>10;
#undef  MAKEFRAC
	}
#endif
	//p1+=((prev4<<USEBITS)-p1+(1<<5>>1))>>5;
	//p1+=((prev3<<USEBITS)-p1+(1<<5>>1))>>5;
	//p1+=((prev2<<USEBITS)-p1+(1<<5>>1))>>5;
	//p1+=((prev1<<USEBITS)-p1+(1<<5>>1))>>5;
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
	x=ac->lo+(x*(uint32_t)p1>>USEBITS);
	bit=*pbit;
	bit=fwd?bit:ac->code<x;
	*pbit=bit;
#ifdef _MSC_VER
	if((uint32_t)(p1-1)>(uint32_t)((1<<USEBITS)-2))
		CRASH("Invalid p1 0x%08X / %d bit", p1, USEBITS);
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue(bit<<24^p1);
	else
		fifoval_check(bit<<24^p1);
#endif
#endif
	*(bit?&ac->hi:&ac->lo)=x-bit;
	prob+=(int32_t)((bit<<STOREBITS)+bias)>>sh;
	cell=prob<<CTRBITS|count;
	cell=cell<<HISTBITS|prev4|bit;
	*pcell=cell;
}
static void subband_code_riceac(int16_t *image
	, int iw, int ih
	, int x1, int x2, int y1, int y2
	, int16_t *pixels, int psize
	, int fwd, int LL, int niter, int level
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
					eW	=rows[0][0-1*NCH*NROWS*NVAL];
				int32_t *currhist;
				int ctx, nbypass;

				nbypass=31^LZCNT32(eW+1);
			//	ctx=LL?nbypass:2*((eW|eN)!=0)+(((niter>>1)==(level>>1)?eNW|eNE:image[3*(iw*(ky>>1)+(kx>>1))+kc])!=0);
			//	ctx=LL?nbypass:2*((eW|eN)!=0)+(image[3*(iw*(ky>>1)+(kx>>1))+kc]!=0);
				ctx=LL?nbypass:2*((eW|eN)!=0)+((eNW|eNE)!=0);
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
		UPSCALE16(1.0),
		UPSCALE16(1.0),
		UPSCALE16(1.1),
		UPSCALE16(1.2),
		UPSCALE16(1.4),
		UPSCALE16(1.8),
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
			quantization(image, iw, ih, x1, x2, y1, y2, qsteps, fwd);
	}
	subband_code_riceac(image, iw, ih, x1, x2, y1, y2, pixels, psize, fwd, LL, niter, level, &streamptr, streamend);
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
		ict_fwd(image, iw, ih, ict);
		dwt_cdf97((int16_t*)image, iw, ih, niter, dist, fwd);
#ifdef _MSC_VER
	//	ppm_save_16as8("20260808Sa_0616pm.ppm", (int16_t*)image, iw, ih);//
	//	ppm_save_16bit("20260808Sa_0616pm.ppm", (int16_t*)image, iw, ih, 0xFFFF);//
#endif
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
			subband_proc((int16_t*)image, iw, ih,  0, w1, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SW
			subband_proc((int16_t*)image, iw, ih, w1, w2,  0, h1, pixels, psize, 2*niter, 2*(niter-1-it)+1, 0, fwd, dist, ctxbuf, &streamptr, streamend);//NE
			subband_proc((int16_t*)image, iw, ih, w1, w2, h1, h2, pixels, psize, 2*niter, 2*(niter-1-it)+0, 0, fwd, dist, ctxbuf, &streamptr, streamend);//SE HH
		}
	}
	else
	{
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
		}

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
