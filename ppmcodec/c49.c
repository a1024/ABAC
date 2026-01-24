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
#include<math.h>
#include<sys/stat.h>
#if defined _MSC_VER || defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD
	#define PROFILE_SIZE
	#define MEASURE_PSNR
	#define FIFOVAL
#endif

	#define RGB_PROCESSING
	#define USE_LEGALL53


enum
{
	MAXSTAGES=5,
	
	NCTX=18,
	NLEVELS=512,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=1,
};

//runtime
#if 1
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#ifdef _MSC_VER
#	define	ALIGN(N) __declspec(align(N))
#	define AWM_INLINE __forceinline static
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define AWM_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
#else
AWM_INLINE int floor_log2_64(uint64_t n)
{
	int	logn=-!n;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
AWM_INLINE int floor_log2_32(uint32_t n)
{
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?floor_log2_64(X):floor_log2_32((uint32_t)(X)))
#endif
#define CVTFP32_I32(X)  _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X)  _mm_cvtsd_si64(_mm_set_sd(X))
#define CVTTFP64_I64(X) _mm_cvttsd_si64(_mm_set_sd(X))
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
#endif

#ifdef FIFOVAL
static ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
static uint32_t *fifoval=0;
static void valfifo_enqueue(uint32_t val)
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
static void valfifo_check(uint32_t val)
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
				printf(
					"%10td  0x%08X"
					, idx
					, fifoval[idx]
				);
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

#ifdef MEASURE_PSNR
static uint8_t *g_image=0;
static void psnr_save(uint8_t *image, int iw, int ih)
{
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void psnr_measure(uint8_t *image, int iw, int ih)
{
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	double rmse[4]={0}, psnr[4]={0}, invres=3./size;
	for(ptrdiff_t k=0;k<size;k+=3)
	{
		int diff;

		diff=image[k+0]-g_image[k+0]; rmse[0]+=diff*diff;
		diff=image[k+1]-g_image[k+1]; rmse[1]+=diff*diff;
		diff=image[k+2]-g_image[k+2]; rmse[2]+=diff*diff;
	}
	rmse[3]=sqrt((rmse[0]+rmse[1]+rmse[2])/size);
	rmse[0]=sqrt(rmse[0]*invres);
	rmse[1]=sqrt(rmse[1]*invres);
	rmse[2]=sqrt(rmse[2]*invres);
	psnr[0]=-20*log10(rmse[0]/255);
	psnr[1]=-20*log10(rmse[1]/255);
	psnr[2]=-20*log10(rmse[2]/255);
	psnr[3]=-20*log10(rmse[3]/255);
	printf("RMSE PSNR\n");
	printf("T %12.6lf %12.6lf\n", rmse[3], psnr[3]);
	printf("Y %12.6lf %12.6lf\n", rmse[0], psnr[0]);
	printf("U %12.6lf %12.6lf\n", rmse[1], psnr[1]);
	printf("V %12.6lf %12.6lf\n", rmse[2], psnr[2]);
}
#endif

#ifdef PROFILE_SIZE
double csizes[3];
#endif
static const float quantdc[3]=
{//	y, u, v
#ifdef RGB_PROCESSING
	4, 4, 4
#else
//	1.5, 1, 1
	1, 1, 1
#endif
};
static const float quant[3]=
{//	y, u, v
#ifdef RGB_PROCESSING
	10, 10, 10
#else
	6, 2, 2
#endif
};
static const float CDF97coeffs[]=
{
	//'factoring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	  1.58613434342059f,	//alpha		predict is negative
	 -0.0529801185729f,	//beta		update is positive
	 -0.8829110755309f,	//gamma		predict
	  0.4435068520439f,	//delta		update
	  1.230174105f,		//K

//	  1.1496043988602f,	//zeta		output gain is 1.89
};
/*
neven >= nodd

0	1	2	3	[4]
e	o	e	o	[e]
     /     \         /     \
    2       \       /      [2]
   /         \     /         \
e2	o	e2	o	[e2]	predict even (HF)
   \         /     \         /
    \       /     ![2]     [/]
     \     /         \     /
e2	o2	e2	o2	[e2]	update odd (LF)
*/

static void dwt1d_predict(float *buf, int count, float predict)
{
	int ecount=count&~1;
	__m128 mpred=_mm_set1_ps(predict);
	{
		__m128 me=_mm_load_ps(buf+0*4);
		__m128 mo=_mm_load_ps(buf+1*4);
		mo=_mm_add_ps(mo, mo);
		mo=_mm_mul_ps(mo, mpred);
		me=_mm_sub_ps(me, mo);
		_mm_store_ps(buf+0*4, me);
	}
	if(ecount>2)
	{
		__m128 mol=_mm_load_ps(buf+4*(2-1));
		for(int k=2;k<ecount;k+=2)
		{
			__m128 me=_mm_load_ps(buf+4*k);
			__m128 mor=_mm_load_ps(buf+4*(k+1));
			mol=_mm_add_ps(mol, mor);
			mol=_mm_mul_ps(mol, mpred);
			me=_mm_sub_ps(me, mol);
			_mm_store_ps(buf+4*k, me);
			mol=mor;
		}
	}
	if(count&1)
	{
		__m128 mo=_mm_load_ps(buf+(count-2)*4);
		__m128 me=_mm_load_ps(buf+(count-1)*4);
		mo=_mm_add_ps(mo, mo);
		mo=_mm_mul_ps(mo, mpred);
		me=_mm_sub_ps(me, mo);
		_mm_store_ps(buf+(count-1)*4, me);
	}
}
static void dwt1d_update(float *buf, int count, float update)
{
	__m128 mupdate=_mm_set1_ps(update);
	if(count>2)
	{
		__m128 mel=_mm_load_ps(buf+0*4);
		for(int k=1;k<count-1;k+=2)
		{
			__m128 mo=_mm_load_ps(buf+k*4);
			__m128 mer=_mm_load_ps(buf+(k+1)*4);
			mel=_mm_add_ps(mel, mer);
			mel=_mm_mul_ps(mel, mupdate);
			mo=_mm_add_ps(mo, mel);
			_mm_store_ps(buf+k*4, mo);
			mel=mer;
		}
	}
	if(!(count&1))//even count
	{
		__m128 me=_mm_load_ps(buf+(count-2)*4);
		__m128 mo=_mm_load_ps(buf+(count-1)*4);
		me=_mm_add_ps(me, me);
		me=_mm_mul_ps(me, mupdate);
		mo=_mm_add_ps(mo, me);
		_mm_store_ps(buf+(count-1)*4, mo);
	}
}
static void dwt1d_scale(float *buf, int count, float higain, float logain)
{
	__m128 mghi=_mm_set1_ps(higain);
	__m128 mglo=_mm_set1_ps(logain);
	for(int k=0;k<count-1;k+=2)
	{
		__m128 me=_mm_load_ps(buf+4*k);
		__m128 mo=_mm_load_ps(buf+4*(k+1));
		me=_mm_mul_ps(me, mghi);
		mo=_mm_mul_ps(mo, mglo);
		_mm_store_ps(buf+4*k, me);
		_mm_store_ps(buf+4*k+4, mo);
	}
	if(count&1)//odd
	{
		__m128 me=_mm_load_ps(buf+4*(count-1));
		me=_mm_mul_ps(me, mghi);
		_mm_store_ps(buf+4*(count-1), me);
	}
}

static uint32_t hweight[3*NCTX];
static uint32_t hists[3*NCTX*NLEVELS];
int c49_codec(int argc, char **argv)
{
	const uint16_t tag='4'|'9'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0, rowstride=0;
	int64_t usize=0, ccap=0, csize=0;
	ptrdiff_t tsize=0;
	int t2size=0;
	float *tbuf=0, *tbuf2=0;
	uint8_t *image=0, *stream=0, *streamptr=0;
	int32_t *imptr=0;
#ifdef _MSC_VER
	uint8_t *streamend=0;
#endif
	int32_t stages[MAXSTAGES][2]={0}, nstages=0;
	int32_t ctx=0, error=0, sym=0, cdf=0, freq=0;
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef LOUD
	double t=0;
#endif
#ifdef _MSC_VER
	float vmin[3][2]={0}, vmax[3][2]={0};
#endif

	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  src  dst\n"
			"Only for 24-bit PPM images\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	srcfn=argv[1];
	dstfn=argv[2];
	
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
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
		{
			struct stat info={0};

			stat(srcfn, &info);
			ccap=(int64_t)info.st_size-ftell(fsrc);
		}
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	rowstride=3*iw;
	usize=(int64_t)3*iw*ih;
	image=(uint8_t*)malloc(usize);
	stream=(uint8_t*)malloc(ccap);
	tsize=(ptrdiff_t)iw*ih*sizeof(float[3]);
	tbuf=(float*)malloc(tsize);
	t2size=(iw>ih?iw:ih)*(int)sizeof(__m128);
	tbuf2=(float*)_mm_malloc(t2size, sizeof(__m128));
	if(!image||!stream||!tbuf||!tbuf2)
	{
		CRASH("Alloc error");
		return 1;
	}
	streamptr=stream;
#ifdef _MSC_VER
	streamend=stream+ccap;
#endif
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
#ifdef MEASURE_PSNR
		psnr_save(image, iw, ih);
#endif
	}
	else
	{
		c=fread(stream, 1, ccap, fsrc);
	}
	fclose(fsrc);
	
#ifdef PROFILE_SIZE
	memset(csizes, 0, sizeof(csizes));
#endif
	nstages=0;
	for(int w2=iw, h2=ih;nstages<MAXSTAGES;++nstages)
	{
		if(w2<=8||ih<=8)
			break;
		stages[nstages][0]=w2;
		stages[nstages][1]=h2;
		w2>>=1;
		h2>>=1;
	}
	memset(hists, 0, sizeof(hists));
	memset(hweight, 0, sizeof(hweight));
	if(fwd)//fwd DWT + quantization
	{
		//extract YUV
		{
			static const float ycbcr[]=
			{
#ifdef RGB_PROCESSING
				1, 0, 0,
				0, 1, 0,
				0, 0, 1,
#else
				0.299f,
				0.587f,
				0.114f,
				-0.168736f,
				-0.331264f,
				0.5f,
				0.5f,
				-0.418688f,
				-0.081312f,
#endif
			};
			float *dstptr=tbuf;
			uint8_t *srcptr=image;
			for(int ky=0;ky<ih;++ky)
			{
				for(int kx=0;kx<iw;++kx, srcptr+=3)
				{
					*dstptr++=
						+ycbcr[0]*srcptr[0]
						+ycbcr[1]*srcptr[1]
						+ycbcr[2]*srcptr[2]
						-128.f
					;
					*dstptr++=
						+ycbcr[3]*srcptr[0]
						+ycbcr[4]*srcptr[1]
						+ycbcr[5]*srcptr[2]
					;
					*dstptr++=
						+ycbcr[6]*srcptr[0]
						+ycbcr[7]*srcptr[1]
						+ycbcr[8]*srcptr[2]
					;
				}
			}
		}

		//apply fwd DWT
		for(int kc=0;kc<3;++kc)
		{
			for(int ks=0;ks<nstages;++ks)
			{
				int w2=stages[ks][0], h2=stages[ks][1];
#if 1
				for(int ky=0;ky<h2;ky+=4)//horizontal
				{
					int nrows=h2-ky;
					if(nrows>4)
						nrows=4;

					//load N rows
					memset(tbuf2, 0, w2*sizeof(float[4]));
					for(int kr=0;kr<nrows;++kr)
					{
						float *srcptr=tbuf+rowstride*(ky+kr)+kc;
						float *dstptr=tbuf2+kr;
						for(int kx=0;kx<w2;++kx, srcptr+=3, dstptr+=4)
							*dstptr=*srcptr;
					}
#ifdef USE_LEGALL53
					dwt1d_predict	(tbuf2, w2, 0.5f);
					dwt1d_update	(tbuf2, w2, 0.25f);
				//	dwt1d_scale	(tbuf2, w2, 1.f/(4+0.125f*(nstages-ks)), 1);
#else
					dwt1d_predict	(tbuf2, w2, CDF97coeffs[0]);
					dwt1d_update	(tbuf2, w2, CDF97coeffs[1]);
					dwt1d_predict	(tbuf2, w2, CDF97coeffs[2]);
					dwt1d_update	(tbuf2, w2, CDF97coeffs[3]);
					dwt1d_scale	(tbuf2, w2, CDF97coeffs[4], 1/CDF97coeffs[4]);
#endif
					//store N rows with lazy wavelet
					for(int kr=0;kr<nrows;++kr)
					{
						float *dstptr=tbuf+rowstride*(ky+kr)+kc;
						float *srcptr=tbuf2+kr;
						for(int kx=1;kx<w2;kx+=2, dstptr+=3)//odd (LF)
							*dstptr=srcptr[4*kx];
						for(int kx=0;kx<w2;kx+=2, dstptr+=3)//even (HF)
							*dstptr=srcptr[4*kx];
					}
				}
#endif
#if 1
				for(int kx=0;kx<w2;kx+=4)//vertical
				{
					int ncols=w2-kx;
					if(ncols>4)
						ncols=4;

					//load N cols
					memset(tbuf2, 0, h2*sizeof(float[4]));
					for(int kr=0;kr<ncols;++kr)
					{
						float *srcptr=tbuf+3*(kx+kr)+kc;
						float *dstptr=tbuf2+kr;
						for(int kx=0;kx<h2;++kx, dstptr+=4, srcptr+=rowstride)
							*dstptr=*srcptr;
					}
#ifdef USE_LEGALL53
					dwt1d_predict	(tbuf2, h2, 0.5f);
					dwt1d_update	(tbuf2, h2, 0.25f);
				//	dwt1d_scale	(tbuf2, w2, 1.f/(4+0.125f*(nstages-ks)), 1);
#else
					dwt1d_predict	(tbuf2, h2, CDF97coeffs[0]);
					dwt1d_update	(tbuf2, h2, CDF97coeffs[1]);
					dwt1d_predict	(tbuf2, h2, CDF97coeffs[2]);
					dwt1d_update	(tbuf2, h2, CDF97coeffs[3]);
					dwt1d_scale	(tbuf2, h2, CDF97coeffs[4], 1/CDF97coeffs[4]);
#endif
					//store N rows with lazy wavelet
					for(int kr=0;kr<ncols;++kr)
					{
						float *dstptr=tbuf+3*(kx+kr)+kc;
						float *srcptr=tbuf2+kr;
						for(int kx=1;kx<h2;kx+=2, dstptr+=rowstride)
							*dstptr=srcptr[4*kx];
						for(int kx=0;kx<h2;kx+=2, dstptr+=rowstride)
							*dstptr=srcptr[4*kx];
					}
				}
#endif
			}
		}

		//quantize
		{
			int dcw=stages[nstages-1][0];
			int dch=stages[nstages-1][1];
			float *srcptr=tbuf;
			for(int ky=0;ky<ih;++ky)
			{
				for(int kx=0;kx<iw;++kx)
				{
					int ac=ky>=dch||kx>=dcw;
					float q[]=
					{
						ac?1.f/quant[0]:1.f/quantdc[0],
						ac?1.f/quant[1]:1.f/quantdc[1],
						ac?1.f/quant[2]:1.f/quantdc[2],
					};
					for(int kc=0;kc<3;++kc)
					{
						float val=*srcptr;
#ifdef _MSC_VER
						if(vmin[kc][ac]>val)vmin[kc][ac]=val;
						if(vmax[kc][ac]<val)vmax[kc][ac]=val;
#endif
						val*=q[kc];
						if(ac)
						{
							if(fabsf(val)<1)//deadzone
								val=0;
						}
//#ifdef _MSC_VER
//						if(fabsf(val)>256)
//							CRASH("");
//#endif
					//	CLAMP2(val, -NLEVELS/2, NLEVELS/2-1);
						*(int32_t*)srcptr=CVTFP32_I32(val);
						++srcptr;
					}
				}
			}
		}
	}
	else
	{
		code=0;
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);//load
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
	}

	//entropy coding
	int psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	uint16_t *pixels=(uint16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	imptr=(int32_t*)tbuf;
	for(int ky=0;ky<ih;++ky)
	{
		uint32_t *currhist=0;
		int den=0;
		uint16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc, ++imptr)
			{
				int
					eNEE	=rows[1][+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][+3*NCH*NROWS*NVAL],
					eW	=rows[0][-1*NCH*NROWS*NVAL];
				ctx=FLOOR_LOG2(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;

				ctx+=NCTX*kc;
				currhist=hists+NLEVELS*ctx;
				den=hweight[ctx]+NLEVELS;
				if(fwd)
				{
					int t;

					error=*imptr;
#ifdef _MSC_VER
					if(abs(error)>NLEVELS/2-1)
						CRASH("");
#endif
					CLAMP2(error, -NLEVELS/2, NLEVELS/2-1);
					sym=error<<1^error>>31;
					if(range<=0xFFFF)
					{
						*(uint32_t*)streamptr=(uint32_t)(low>>32);
						streamptr+=4;
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					//if(ky==7&&kx==946&&kc==0)//
					//	printf("");
					for(t=0, cdf=0;;++t)
					{
						freq=currhist[t]+1;
						if(t>=sym)
							break;
						cdf+=freq;
					}
#ifdef _MSC_VER
					if(!freq||(uint32_t)freq>(uint32_t)den||(uint32_t)cdf>(uint32_t)den)
						CRASH("");
#endif
#ifdef FIFOVAL
					valfifo_enqueue(freq<<16|cdf);
#endif
					low+=range*cdf/den;
					range=range*freq/den-1;
#ifdef PROFILE_SIZE
					csizes[kc]+=-log2((double)freq/den);
#endif
				}
				else
				{
					int c2;

					if(range<=0xFFFF)//stall: unpredictable branch
					{
#ifdef _MSC_VER
						if(streamptr>streamend)
							CRASH("");
#endif
						code=code<<32|*(uint32_t*)streamptr;
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					c2=(int)(((code-low+1)*den-1)/range);
					sym=0;
					cdf=0;
					for(;;)
					{
						freq=currhist[sym]+1;
						if(cdf+freq>c2)
							break;
						cdf+=freq;
						++sym;
					}
#ifdef FIFOVAL
					valfifo_check(freq<<16|cdf);
#endif
					low+=range*cdf/den;
					range=range*freq/den-1;

					error=sym>>1^-(sym&1);
					*imptr=error;
				}
				++currhist[sym];
				++hweight[ctx];
				if(hweight[ctx]>=0x4000)//rescale
				{
					den=0;
					for(int k=0;k<NLEVELS;++k)
						den+=currhist[k]>>=1;
					hweight[ctx]=den;
				}
				rows[0][0]=(2*eW+(sym<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
		}
	}
	if(!fwd)//dequantization + inv DWT
	{
		//dequantize
		{
			int dcw=stages[nstages-1][0];
			int dch=stages[nstages-1][1];
			float *dstptr=tbuf;
			int32_t *srcptr=(int32_t*)tbuf;
			for(int ky=0;ky<ih;++ky)
			{
				for(int kx=0;kx<iw;++kx)
				{
					if(ky<dch&&kx<dcw)
					{
						*dstptr++=quantdc[0]*(float)*srcptr++;
						*dstptr++=quantdc[1]*(float)*srcptr++;
						*dstptr++=quantdc[2]*(float)*srcptr++;
						continue;
					}
					*dstptr++=quant[0]*(float)*srcptr++;
					*dstptr++=quant[1]*(float)*srcptr++;
					*dstptr++=quant[2]*(float)*srcptr++;
				}
			}
		}
		
		//apply inv DWT
		for(int kc=0;kc<3;++kc)
		{
			for(int ks=0;ks<nstages;++ks)
			{
				int w2=stages[nstages-1-ks][0], h2=stages[nstages-1-ks][1];
#if 1
				for(int kx=0;kx<w2;kx+=4)//vertical
				{
					int ncols=w2-kx;
					if(ncols>4)
						ncols=4;
					
					//load N rows with lazy wavelet
					memset(tbuf2, 0, h2*sizeof(float[4]));
					for(int kr=0;kr<ncols;++kr)
					{
						float *srcptr=tbuf+3*(kx+kr)+kc;
						float *dstptr=tbuf2+kr;
						for(int kx=1;kx<h2;kx+=2, srcptr+=rowstride)
							dstptr[4*kx]=*srcptr;
						for(int kx=0;kx<h2;kx+=2, srcptr+=rowstride)
							dstptr[4*kx]=*srcptr;
					}
#ifdef USE_LEGALL53
				//	dwt1d_scale	(tbuf2, w2, (float)(4+0.125f*(nstages-ks)), 1);
					dwt1d_update	(tbuf2, h2, -0.25f);
					dwt1d_predict	(tbuf2, h2, -0.5f);
#else
					dwt1d_scale	(tbuf2, h2, 1/CDF97coeffs[4], CDF97coeffs[4]);
					dwt1d_update	(tbuf2, h2, -CDF97coeffs[3]);
					dwt1d_predict	(tbuf2, h2, -CDF97coeffs[2]);
					dwt1d_update	(tbuf2, h2, -CDF97coeffs[1]);
					dwt1d_predict	(tbuf2, h2, -CDF97coeffs[0]);
#endif
					//store N rows
					for(int kr=0;kr<ncols;++kr)
					{
						float *dstptr=tbuf+3*(kx+kr)+kc;
						float *srcptr=tbuf2+kr;
						for(int kx=0;kx<h2;++kx, dstptr+=rowstride, srcptr+=4)
							*dstptr=*srcptr;
					}
				}
#endif
#if 1
				for(int ky=0;ky<h2;ky+=4)//hoizontal
				{
					int nrows=h2-ky;
					if(nrows>4)
						nrows=4;
					
					//load N rows with lazy wavelet
					memset(tbuf2, 0, w2*sizeof(float[4]));
					for(int kr=0;kr<nrows;++kr)
					{
						float *srcptr=tbuf+rowstride*(ky+kr)+kc;
						float *dstptr=tbuf2+kr;
						for(int kx=1;kx<w2;kx+=2, srcptr+=3)
							dstptr[4*kx]=*srcptr;
						for(int kx=0;kx<w2;kx+=2, srcptr+=3)
							dstptr[4*kx]=*srcptr;
					}
#ifdef USE_LEGALL53
				//	dwt1d_scale	(tbuf2, w2, (float)(4+0.125f*(nstages-ks)), 1);
					dwt1d_update	(tbuf2, w2, -0.25f);
					dwt1d_predict	(tbuf2, w2, -0.5f);
#else
					dwt1d_scale	(tbuf2, w2, 1/CDF97coeffs[4], CDF97coeffs[4]);
					dwt1d_update	(tbuf2, w2, -CDF97coeffs[3]);
					dwt1d_predict	(tbuf2, w2, -CDF97coeffs[2]);
					dwt1d_update	(tbuf2, w2, -CDF97coeffs[1]);
					dwt1d_predict	(tbuf2, w2, -CDF97coeffs[0]);
#endif
					//store N rows
					for(int kr=0;kr<nrows;++kr)
					{
						float *dstptr=tbuf+rowstride*(ky+kr)+kc;
						float *srcptr=tbuf2+kr;
						for(int kx=0;kx<w2;++kx, dstptr+=3, srcptr+=4)
							*dstptr=*srcptr;
					}
				}
#endif
			}
		}

		//inv YCbCr
		{
			static const float invYCbCr[]=
			{
#ifdef RGB_PROCESSING
				1, 0, 0,
				0, 1, 0,
				0, 0, 1,
#else
				1,
				0,
				1.402f,
				1,
				-0.344136f,
				-0.714136f,
				1,
				1.772f,
				0,
#endif
			};
			float *srcptr=tbuf;
			uint8_t *dstptr=image;
			for(int ky=0;ky<ih;++ky)
			{
				for(int kx=0;kx<iw;++kx, srcptr+=3)
				{
					float y=srcptr[0]+128.f, u=srcptr[1], v=srcptr[2];
					float val;

					val=
						+invYCbCr[0]*y
						+invYCbCr[1]*u
						+invYCbCr[2]*v
					;
					CLAMP2(val, 0, 255);
					*dstptr++=CVTFP32_I32(val);

					val=
						+invYCbCr[3]*y
						+invYCbCr[4]*u
						+invYCbCr[5]*v
					;
					CLAMP2(val, 0, 255);
					*dstptr++=CVTFP32_I32(val);

					val=
						+invYCbCr[6]*y
						+invYCbCr[7]*u
						+invYCbCr[8]*v
					;
					CLAMP2(val, 0, 255);
					*dstptr++=CVTFP32_I32(val);
				}
			}
		}
	}
	free(pixels);
	free(tbuf);
	_mm_free(tbuf2);

	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;
		
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		csize+=fwrite(stream, 1, streamptr-stream, fdst);
	}
	else
	{
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(image, 1, usize, fdst);
	}
	fclose(fdst);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef PROFILE_SIZE
		double esize=csizes[0]+csizes[1]+csizes[2];
		printf("T %12.2lf\n", esize/8);
		printf("Y %12.2lf\n", csizes[0]/8);
		printf("U %12.2lf\n", csizes[1]/8);
		printf("V %12.2lf\n", csizes[2]/8);
#endif
		printf("%d stages\n", nstages);
		printf("DC size WH %d*%d  scale %dX\n", stages[nstages-1][0], stages[nstages-1][1], 1<<nstages);

		printf("Y DC %12.2f ~ %12.2f    AC %12.2f ~ %12.2f\n", vmin[0][0], vmax[0][0], vmin[0][1], vmax[0][1]);
		printf("U DC %12.2f ~ %12.2f    AC %12.2f ~ %12.2f\n", vmin[1][0], vmax[1][0], vmin[1][1], vmax[1][1]);
		printf("V DC %12.2f ~ %12.2f    AC %12.2f ~ %12.2f\n", vmin[2][0], vmax[2][0], vmin[2][1], vmax[2][1]);

		printf("WH %5d*%5d  \"%s\"\n", iw, ih, srcfn);
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
			, 8.*csize/usize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
	if(!fwd)
	{
#ifdef MEASURE_PSNR
		psnr_measure(image, iw, ih);
		free(g_image);
#endif
	}
#endif
	free(image);
	(void)&time_sec;
	(void)csize;
	(void)&dwt1d_scale;
	(void)&CDF97coeffs;
	return 0;
}
