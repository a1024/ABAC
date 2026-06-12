#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#endif
#include<stddef.h>
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
#include<assert.h>
#if defined _M_X64 || defined __x86_64__
#include<immintrin.h>
#endif


#ifdef _MSC_VER
	#define LOUD
	#define RICE2ESTIM

//	#define FIFOVAL 3
#endif

	#define USE_AC


enum
{
	DEPTH=9,
	RLIMIT=12,
	BUFSIZE=512<<10,
#ifdef USE_AC
	NCTX=16,//POT
	RSHIFT=4,
	PROBBITS=11,
	NLEVELS=1<<DEPTH,
#else
	NCTX=6,
	RSHIFT=3,
#endif

	XPAD=8,
	NCH=3,
	NROWS=1,
	NVAL=2,
};


//runtime
#if 1
#ifndef MEDIAN3V_CLOB
//clobbers A B C
#define MEDIAN3V_CLOB(M, A, B, C)\
	M=A, A=A>B?B:A, B=M>B?M:B,\
	M=B, B=B>C?C:B, C=M>C?M:C,\
	M=A, A=A>B?B:A, M=M>B?M:B
#endif
#ifndef CLAMP2
#define CLAMP2(X, L, H) X=X>(L)?X:L, X=X<(H)?X:H
#endif
#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define INLINE __forceinline static
#	define LIKELY(C) C
#else
#	define ALIGN(N) __attribute__((aligned(N)))
#	define INLINE __attribute__((always_inline)) inline static
#	define LIKELY(C) __builtin_expect(C, 1)
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#if defined _M_X64 || defined __x86_64__
#define LZCNT32 _lzcnt_u32
#define LZCNT64 _lzcnt_u64
#define TZCNT32 _tzcnt_u32
#define TZCNT64 _tzcnt_u64
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))
#else
#define LZCNT32 __builtin_clz
#define LZCNT64 __builtin_clzll
#define TZCNT32 __builtin_ctz
#define TZCNT64 __builtin_ctzll
#define ROUND32(X) (int32_t)roundf(X)
#define ROUND64(X) (int64_t)round(X)
#define TRUNC32(X) (int32_t)truncf(X)
#define TRUNC64(X) (int64_t)trunc(X)
#endif
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
static double time_sec2(void)
{
#ifdef _WIN32
	static int64_t t0=0;
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
static int fifoval_check(uint32_t val)
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
		//return 1;
	}
	return 0;
}
#endif
#endif


static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t[4])]={0}, wtbuf[BUFSIZE+sizeof(uint64_t[4])]={0};
#if 0
INLINE uint64_t acme_read(uint8_t **pptr, const ptrdiff_t size, FILE *f)
{
	uint8_t *ptr=*pptr;
	uint64_t data;
	
	memcpy(&data, ptr, sizeof(uint64_t));
	if(ptr>=rdbuf+sizeof(uint64_t)+BUFSIZE-size)
	{
		uint64_t d2;

		fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, f);
		ptr-=BUFSIZE;
		memcpy(&d2, ptr, sizeof(uint64_t));
		data|=d2;
	}
	*pptr=ptr+size;
	return data;
}
INLINE void acme_write(uint8_t **pptr, const ptrdiff_t size, FILE *f, uint64_t data)
{
	uint8_t *ptr=*pptr;
	
	memcpy(ptr, &data, sizeof(uint64_t));
	if(ptr>=wtbuf+sizeof(uint64_t)+BUFSIZE-size)
	{
		fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, f);
		ptr-=BUFSIZE;
		memcpy(ptr, &data, sizeof(uint64_t));
	}
	*pptr=ptr+size;
}
#endif

#ifdef USE_AC
static uint8_t ctxtable[1<<(DEPTH+RSHIFT)];
static uint32_t cdftable[NCTX*(NLEVELS+1)];
static void print_hist(int32_t *hist)
{
	int hmax=0;
	for(int ks=0;ks<NLEVELS;++ks)
	{
		int freq=hist[ks];
		if(hmax<freq)
			hmax=freq;
	}
	for(int ks=0;ks<16;++ks)
	{
		int freq=hist[ks], nstars=hmax?freq*128/hmax:0;
		printf("%3d  %5d  ", ks, hist[ks]);
		for(int k2=0;k2<nstars;++k2)
			printf("-");
		printf("\n");
	}
}
#else
static uint8_t log2table[1<<(DEPTH+RSHIFT)];
static uint32_t enctable[8<<DEPTH];
#endif
static int16_t packsign[1024], *const packsignptr=packsign+512;
int c59_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'9'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	int64_t c=0;
	int fwd=0, iw=0, ih=0;
	uint8_t *rdptr=0, *wtptr=0;
	int yuv[3]={0};
	int sym[3]={0};
	int64_t usize=0;
	int psize=0;
	int16_t *pixels=0;
#ifdef USE_AC
	uint64_t lo=0, hi=0xFFFFFFFFFFFF, code=0;
#else
	uint64_t cache=0;
	int nbits=0;
#endif
#ifdef LOUD
	double t=time_sec2();
#endif
#ifdef RICE2ESTIM
	//int zcnt[3]={0}, tcnt[3]={0};
	int zrun[3]={0}, prev[3]={0}, xhat0[3]={1<<RSHIFT, 1<<RSHIFT, 1<<RSHIFT}, xhat[3]={1<<RSHIFT, 1<<RSHIFT, 1<<RSHIFT};
	int64_t ricesize=0, rice2size=0;
#endif

	if(argc!=3)
	{
		printf("Usage:  \"%s\"  src  dst\n", argv[0]);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	c=fgetc(fsrc);
	c|=(int64_t)fgetc(fsrc)<<8;
	if(c==('P'|'6'<<8))
		fwd=1;
	else if(c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(fwd)
	{
		c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported file \"%s\"", srcfn);
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
			iw=10*iw+(int)c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		ih=0;
		while((uint32_t)(c-'0')<10)
		{
			ih=10*ih+(int)c-'0';
			c=fgetc(fsrc);
		}
		if(c!='\n')
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		while(c<=' ')
			c=fgetc(fsrc);
		c|=(int64_t)fgetc(fsrc)<< 8;
		c|=(int64_t)fgetc(fsrc)<<16;
		c|=(int64_t)fgetc(fsrc)<<24;
		if(c!=('2'|'5'<<8|'5'<<16|'\n'<<24))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
	}
	else
	{
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		fclose(fsrc);
		fclose(fdst);
		return 1;
	}
	memset(pixels, 0, psize);
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	memset(rdbuf, 0, sizeof(rdbuf));
	memset(wtbuf, 0, sizeof(wtbuf));
	rdptr=rdbuf+sizeof(uint64_t)+BUFSIZE;
	wtptr=wtbuf+sizeof(uint64_t);
#ifdef USE_AC
#if 1
	static const double xestim4[]=
	{
#define ESTIM(X) (X)*(1<<RSHIFT)
		ESTIM(1./32),
	//	ESTIM(1./16),
	//	ESTIM(1./8),
	//	ESTIM(1./4),
	//	ESTIM(1./2),
	//	ESTIM(1),
		ESTIM(2),
		ESTIM(4),
		ESTIM(8),
		ESTIM(16),
		ESTIM(32),
		ESTIM(64),
		ESTIM(128),
		ESTIM(256),
		ESTIM(512),
		ESTIM(1024),
		ESTIM(2048),
		ESTIM(4096),
		ESTIM(8192),
		ESTIM(16384),
		ESTIM(32768),
		ESTIM(65536),
#undef  ESTIM
	};
#endif
#if 0
	static const double xestim3[]=
	{
#define ESTIM(X) (X)*(1<<RSHIFT)
		ESTIM(1),
		ESTIM(1.414),
		ESTIM(2),
		ESTIM(2.828),
		ESTIM(4),
		ESTIM(5.657),
		ESTIM(8),
		ESTIM(11.314),
		ESTIM(16),
		ESTIM(22.627),
		ESTIM(32),
		ESTIM(45.254),
		ESTIM(64),
		ESTIM(90.509),
		ESTIM(128),
		ESTIM(181.019),
		ESTIM(256),
#undef  ESTIM
	};
#endif
#if 0
	static const double xestim[]=
	{
#define ESTIM(X) (X)*(1<<RSHIFT)
		ESTIM(0),
		ESTIM(0.25),
		ESTIM(0.35),
		ESTIM(0.5),
		ESTIM(0.75),
		ESTIM(1),
		ESTIM(1.5),
		ESTIM(2),
		ESTIM(3),
		ESTIM(5),
		ESTIM(8),
		ESTIM(16),
		ESTIM(32),
		ESTIM(64),
		ESTIM(128),
		ESTIM(256),
		ESTIM(512),
#undef  ESTIM
	};
#endif
	double xestim2[NCTX+1];
	for(int ctx=0;ctx<=NCTX;++ctx)
		xestim2[ctx]=exp2((double)ctx*(9-0.6)/NCTX+0.6)*(1<<RSHIFT);
	//	xestim2[ctx]=exp2((double)ctx*(9-0.58)/NCTX+0.58)*(1<<RSHIFT);

	//int ctx0=-1;
	for(int hx=0, ctx=0;hx<_countof(ctxtable);++hx)
	{
		ctx+=hx>xestim4[ctx+1];
		//if(ctx>NCTX-1)
		//	CRASH("");
		ctxtable[hx]=ctx;
#if 0
		double x, nbits;
		int ctx;

		x=(double)hx;
		if(x<1./(1<<(PROBBITS+RSHIFT)))
			x=1./(1<<(PROBBITS+RSHIFT));
		nbits=log2(x)*0.25;
		ctx=(int)ROUND64(nbits);
		//if(ctx<0)
		//	ctx=0;
		//int ctx=(31-6)-LZCNT32(hx*hx+1);
		CLAMP2(ctx, 0, NCTX-1);
		ctxtable[hx]=ctx;
		//if(!hx||!ctx||hx==_countof(ctxtable)-1)
		if(ctx!=ctx0)
			printf("hx %5d  ctx %2d\n", hx, ctx);
		ctx0=ctx;
#endif
	}
	for(int ctx=0;ctx<NCTX;++ctx)
	{
		int32_t hist[NLEVELS+1]={0}, hsum=0, c=0;
		double invx=(1<<RSHIFT)/xestim4[ctx+1], codelen=0;
		//double invx=(2<<RSHIFT)/(xestim4[ctx]+xestim4[ctx+1]), codelen=0;
		//double invx=exp2((double)-(ctx-6)), codelen=0;
		//double invx=exp2(ctx*0.5), codelen=0;
		//if(invx<1./(1<<PROBBITS))
		//	invx=1./(1<<PROBBITS);
		//invx=1./invx;
		for(int ks=0;ks<NLEVELS;++ks)
		{
			int freq;
			//double scale=ks;
			//
			//CLAMP2(scale, 1, 2);
			codelen+=invx;
			if(codelen>PROBBITS)
				codelen=PROBBITS;
			freq=(int)ROUND64(exp2(-(double)TRUNC64(codelen-PROBBITS+2)+2));
		//	freq=(int)ROUND64(exp2(-codelen)*(1<<PROBBITS));
			hsum+=hist[ks]=freq;
		}
		if(fwd)//
		{
			printf("ctx %3d  x  %12.6lf\n", ctx, 1/invx);
			print_hist(hist);//
		}
		for(int ks=0;ks<=NLEVELS;++ks)
		{
			int32_t freq=hist[ks];
			int32_t cdf=hsum?(int32_t)(((uint64_t)c<<PROBBITS)/hsum):ks<<PROBBITS>>DEPTH;
			CLAMP2(cdf, ks, (1<<PROBBITS)-NLEVELS+ks);
			//if(cdf>(1<<PROBBITS)-NLEVELS+ks)
			//	cdf=(1<<PROBBITS)-NLEVELS+ks;
			cdftable[ks<<RSHIFT|ctx]=cdf;

			if(ks&&cdftable[(ks-1)<<RSHIFT|ctx]>=cdftable[ks<<RSHIFT|ctx])//
				CRASH("Invalid CDF");
			if(ks<NLEVELS&&cdftable[ks<<RSHIFT|ctx]>(1<<PROBBITS)-1)//
				CRASH("Invalid CDF");

			c+=freq;
		}
	}
#if 0
	//for(int k=0;k<NCTX;++k)
	//{
	//	int nbypass=31-LZCNT32((bounds[k]>>RSHIFT)+1);
	//	printf("%2d  %2d  %5d  %3d\n", k, nbypass, bounds[k], bounds[k]>>RSHIFT);
	//}
	for(int ctx=0;ctx<NCTX;++ctx)
	{
		int32_t hist[NLEVELS+1]={0}, hsum=0, c=0;
		double invx=1/sqrt(exp2(ctx)-1+(1./(1<<PROBBITS))), codelen=0;
		//double nbypass=(double)ctx/2-RSHIFT+1, codelen=0;
		//double nbypass=log2((pow(2, ctx*0.5)-1)/(1<<RSHIFT)+1);
		for(int ks=0;ks<NLEVELS;++ks)
		{
			int freq;

			freq=(int)ROUND64(exp2(-codelen)*(1<<PROBBITS));
			//if(freq<2)
			//	freq=2;
			hsum+=hist[ks]=freq;
			codelen+=invx;
			if(codelen>PROBBITS)
				codelen=PROBBITS;
		}
		if(fwd)//
		{
			printf("ctx %3d  x  %12.6lf\n", ctx, 1/invx);
			print_hist(hist);//
		}
		for(int ks=0;ks<=NLEVELS;++ks)
		{
			int32_t freq=hist[ks];
			int32_t cdf=(c<<PROBBITS)/hsum;
			if(cdf>(1<<PROBBITS)-NLEVELS+ks)
				cdf=(1<<PROBBITS)-NLEVELS+ks;
			cdftable[ks<<RSHIFT|ctx]=cdf;
			if(ks&&cdftable[(ks-1)<<RSHIFT|ctx]>=cdftable[ks<<RSHIFT|ctx])//
				CRASH("Invalid CDF");
			//if(cdftable[ks<<RSHIFT|ctx]>(1<<PROBBITS)-1)//
			//	CRASH("Invalid CDF");
			c+=freq;
		}
	}
#endif
#if 1
	if(fwd)
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-3;
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint32_t);//AC renorm

		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<(32-DEPTH)>>(32-DEPTH);
			val=val<<1^val>>31;
			packsign[k]=val<<RSHIFT;
		}
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			for(int kx=0;kx<iw;++kx)
			{
				int lookup[3], cdf[3], freq[3];
				uint32_t reg;

				reg=*(uint32_t*)rdptr;
				if(rdptr>=rdend)
				{
					fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
					rdptr-=BUFSIZE;
					reg|=*(uint32_t*)rdptr;
				}
				rdptr+=3;
				yuv[0]=(uint8_t)(reg>>0*8);
				yuv[1]=(uint8_t)(reg>>1*8);
				yuv[2]=(uint8_t)(reg>>2*8);
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
				yuv[1]+=(yuv[0]+yuv[2])>>2;
				yuv[2]-=yuv[0]>>2;
				sym[0]=packsignptr[yuv[0]-((rptr[0+(0+0*NCH)*NROWS*NVAL]+rptr[0+(0-1*NCH)*NROWS*NVAL])>>1)];
				sym[1]=packsignptr[yuv[1]-((rptr[0+(1+0*NCH)*NROWS*NVAL]+rptr[0+(1-1*NCH)*NROWS*NVAL])>>1)];
				sym[2]=packsignptr[yuv[2]-((rptr[0+(2+0*NCH)*NROWS*NVAL]+rptr[0+(2-1*NCH)*NROWS*NVAL])>>1)];
				lookup[0]=sym[0]+ctxtable[rptr[1+(0-1*NCH)*NROWS*NVAL]];
				lookup[1]=sym[1]+ctxtable[rptr[1+(1-1*NCH)*NROWS*NVAL]];
				lookup[2]=sym[2]+ctxtable[rptr[1+(2-1*NCH)*NROWS*NVAL]];
#ifdef RICE2ESTIM
				{
					int nbypass[3];
					int s2[]=
					{
						sym[0]>>RSHIFT,
						sym[1]>>RSHIFT,
						sym[2]>>RSHIFT,
					};
					
					nbypass[0]=31-LZCNT32((xhat[0]>>RSHIFT)+1);
					nbypass[1]=31-LZCNT32((xhat[1]>>RSHIFT)+1);
					nbypass[2]=31-LZCNT32((xhat[2]>>RSHIFT)+1);
					ricesize+=(int64_t)(s2[0]>>nbypass[0])+1+nbypass[0];
					ricesize+=(int64_t)(s2[1]>>nbypass[1])+1+nbypass[1];
					ricesize+=(int64_t)(s2[2]>>nbypass[2])+1+nbypass[2];
					//nbypass[0]=32-LZCNT32(xhat[0]);
					//nbypass[1]=32-LZCNT32(xhat[1]);
					//nbypass[2]=32-LZCNT32(xhat[2]);
					//nbypass[0]=31^LZCNT32((rptr[1+(0-1*NCH)*NROWS*NVAL]>>RSHIFT)+1);
					//nbypass[1]=31^LZCNT32((rptr[1+(1-1*NCH)*NROWS*NVAL]>>RSHIFT)+1);
					//nbypass[2]=31^LZCNT32((rptr[1+(2-1*NCH)*NROWS*NVAL]>>RSHIFT)+1);

					//Rice 2
					for(int kc=0;kc<3;++kc)
					{
						if(!s2[kc])
						{
							if(prev[kc])
								xhat0[kc]=xhat[kc];
							++zrun[kc];
						}
						else
						{
							int ricek, negk, unarylen;

							ricek=32-RSHIFT-LZCNT32(xhat0[kc]);
							negk=-ricek;
							if(ricek<0)ricek=0;
							if(negk<0)negk=0;
							while(zrun[kc]>=1<<negk)
							{
								++rice2size;
								zrun[kc]-=1<<negk;
							}
							unarylen=0;
							while(negk)
							{
								--negk;
								++unarylen;
								if(zrun[kc]&1<<negk)
									rice2size+=(int64_t)unarylen+1;
								zrun[kc]&=~(1<<negk);
							}
							rice2size+=(int64_t)(s2[kc]>>ricek)+1+ricek;
							xhat0[kc]=xhat[kc];
						}
						xhat[kc]+=((s2[kc]<<RSHIFT)-xhat[kc])>>(RSHIFT+1);
						prev[kc]=s2[kc];
					}
					//if(nbypass[0]<RSHIFT)
					//{
					//	int shift=1<<(RSHIFT-nbypass[0]);
					//	if(sym[0])
					//	{
					//		while(zrun[0])
					//		{
					//			int sym2=zrun[0];
					//			if(sym2>shift)
					//				sym2=shift;
					//			rice2size+=;
					//			zrun[0]-=shift;
					//		}
					//		zrun[0]=0;
					//	}
					//}
					//zrun[0]+=!sym[0];
					//zrun[1]+=!sym[1];
					//zrun[2]+=!sym[2];
					//zcnt[0]+=!sym[0];
					//zcnt[1]+=!sym[1];
					//zcnt[2]+=!sym[2];
					//++tcnt[0];
					//++tcnt[1];
					//++tcnt[2];
				}
#endif
				rptr[0+(0+0*NCH)*NROWS*NVAL]=yuv[0];
				rptr[0+(1+0*NCH)*NROWS*NVAL]=yuv[1];
				rptr[0+(2+0*NCH)*NROWS*NVAL]=yuv[2];
				rptr[1+(0+0*NCH)*NROWS*NVAL]=(2*rptr[1+(0-1*NCH)*NROWS*NVAL]+sym[0]+rptr[1+(0+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=(2*rptr[1+(1-1*NCH)*NROWS*NVAL]+sym[1]+rptr[1+(1+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=(2*rptr[1+(2-1*NCH)*NROWS*NVAL]+sym[2]+rptr[1+(2+3*NCH)*NROWS*NVAL])>>2;
				rptr+=NCH*NROWS*NVAL;
				cdf[0]=cdftable[lookup[0]+(0<<RSHIFT)];
				cdf[1]=cdftable[lookup[1]+(0<<RSHIFT)];
				cdf[2]=cdftable[lookup[2]+(0<<RSHIFT)];
				freq[0]=cdftable[lookup[0]+(1<<RSHIFT)]-cdf[0];
				freq[1]=cdftable[lookup[1]+(1<<RSHIFT)]-cdf[1];
				freq[2]=cdftable[lookup[2]+(1<<RSHIFT)]-cdf[2];
#ifdef FIFOVAL
				sym[0]>>=RSHIFT;
				sym[1]>>=RSHIFT;
				sym[2]>>=RSHIFT;
				fifoval_enqueue(sym[2]<<18^sym[1]<<9^sym[0]);
#endif
				for(int kc=0;kc<3;++kc)
				{
					uint64_t x;

					x=hi-lo;
					if(x<=0xFFFF)
					{
						*(uint32_t*)wtptr=(uint32_t)(lo>>32);
						if(wtptr>=wtend)
						{
							fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
							wtptr-=BUFSIZE;
							*(uint32_t*)wtptr=(uint32_t)(lo>>32);
						}
						wtptr+=sizeof(uint32_t);
						lo<<=32;
						hi=hi<<32|0xFFFFFFFF;
						if(hi<lo)
							hi=~0ULL;
						x=hi-lo;
					}
					//if(freq[kc]>(1<<PROBBITS)-1||cdf[kc]>(1<<PROBBITS)-1)//
					//	CRASH("");

					lo+=x*cdf[kc]>>PROBBITS;
					hi=lo+(x*freq[kc]>>PROBBITS)-1;
				}
			}
		}
		lo=lo<<32|lo>>32;
		*(uint64_t*)wtptr=lo;
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			*(uint64_t*)wtptr=lo;
		}
		wtptr+=sizeof(uint64_t);
		(void)rdend;
		(void)wtend;
	}
	else//dec
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint32_t);//AC renorm
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-3;
		int CDF2symsize=(int)sizeof(uint32_t[NCTX<<PROBBITS]);
		uint32_t *CDF2sym=(uint32_t*)malloc(CDF2symsize);
		//uint64_t cache2=0;

		if(!CDF2sym)
			CRASH("Allloc error");
		for(int ctx=0;ctx<NCTX;++ctx)
		{
			int idx=0;
			for(int ks=0;ks<=NLEVELS;++ks)
			{
				int start=cdftable[ks<<RSHIFT|ctx], end=cdftable[(ks+1)<<RSHIFT|ctx], freq=end-start;
				while(idx<end)
					CDF2sym[ctx<<PROBBITS|idx++]=(freq<<PROBBITS|start)<<DEPTH|ks;
			}
			//if(idx!=1<<PROBBITS)//
			//	CRASH("");
		}
		for(int k=0;k<512;++k)
		{
			int val=k;
			val=val>>1^-(val&1);
		//	val=val<<(32-DEPTH)>>(32-DEPTH);
			packsign[k]=val;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fread(&code, 1, 8, fsrc);
		code=code<<32|code>>32;
		//fread(&cache, 1, 8, fsrc);
		//fread(&cache2, 1, 8, fsrc);
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			for(int kx=0;kx<iw;++kx)
			{
				uint32_t reg;
				int ctx[3];

				ctx[0]=ctxtable[rptr[1+(0-1*NCH)*NROWS*NVAL]];
				ctx[1]=ctxtable[rptr[1+(1-1*NCH)*NROWS*NVAL]];
				ctx[2]=ctxtable[rptr[1+(2-1*NCH)*NROWS*NVAL]];
				for(int kc=0;kc<3;++kc)
				{
					uint64_t x;
					int tmp, cdf, freq;

					x=hi-lo;
					if(x<=0xFFFF)//unpredictable branch
					{
						uint32_t c2;

						c2=*(uint32_t*)rdptr;
						if(rdptr>=rdend)
						{

							fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
							rdptr-=BUFSIZE;
							c2|=*(uint32_t*)rdptr;
						}
						code=code<<32|c2;
						rdptr+=sizeof(uint32_t);
						lo<<=32;
						hi=hi<<32|0xFFFFFFFF;
						if(hi<lo)
							hi=~0ULL;
						x=hi-lo;
					}
					tmp=(int)(((code-lo)<<PROBBITS|((1ULL<<PROBBITS)-1))/x);//128-bit div

					//if((uint32_t)tmp>(1<<PROBBITS)-1)//
					//	CRASH("");
					
					tmp=CDF2sym[ctx[kc]<<PROBBITS|tmp];//cache miss
					sym[kc]=tmp&(NLEVELS-1);
					cdf=(tmp>>DEPTH)&((1<<PROBBITS)-1);
					freq=tmp>>(PROBBITS+DEPTH);
					lo+=x*cdf>>PROBBITS;
					hi=lo+(x*freq>>PROBBITS)-1;

					//if(code<lo||code>=hi)//
					//	CRASH("");
				}
#ifdef FIFOVAL
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
					CRASH("");
#endif
				rptr[1+(0+0*NCH)*NROWS*NVAL]=(2*rptr[1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rptr[1+(0+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=(2*rptr[1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rptr[1+(1+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=(2*rptr[1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rptr[1+(2+3*NCH)*NROWS*NVAL])>>2;
				sym[0]=packsign[sym[0]];
				sym[1]=packsign[sym[1]];
				sym[2]=packsign[sym[2]];
				sym[0]+=(rptr[0+(0+0*NCH)*NROWS*NVAL]+rptr[0+(0-1*NCH)*NROWS*NVAL])>>1;
				sym[1]+=(rptr[0+(1+0*NCH)*NROWS*NVAL]+rptr[0+(1-1*NCH)*NROWS*NVAL])>>1;
				sym[2]+=(rptr[0+(2+0*NCH)*NROWS*NVAL]+rptr[0+(2-1*NCH)*NROWS*NVAL])>>1;
				sym[0]<<=32-DEPTH;
				sym[1]<<=32-DEPTH;
				sym[2]<<=32-DEPTH;
				sym[0]>>=32-DEPTH;
				sym[1]>>=32-DEPTH;
				sym[2]>>=32-DEPTH;
				rptr[0+(0+0*NCH)*NROWS*NVAL]=sym[0];
				rptr[0+(1+0*NCH)*NROWS*NVAL]=sym[1];
				rptr[0+(2+0*NCH)*NROWS*NVAL]=sym[2];
				rptr+=NCH*NROWS*NVAL;
				sym[2]+=sym[0]>>2;
				sym[1]-=(sym[0]+sym[2])>>2;
				sym[2]+=sym[1];
				sym[0]+=sym[1];
				reg=sym[2]<<16|sym[1]<<8|sym[0];
				memcpy(wtptr, &reg, sizeof(reg));
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					memcpy(wtptr, &reg, sizeof(uint64_t));
				}
				wtptr+=3;
			}
		}
		(void)rdend;
		(void)wtend;
	}
#endif
#else
	for(int ks=0;ks<1<<(DEPTH+RSHIFT);++ks)
	{
		int val=(ks>>RSHIFT)+1;
		val=LZCNT32(val);
		val^=31;
		if(val>NCTX-1)
			val=NCTX-1;
		log2table[ks]=val;
	}
	if(fwd)
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-3;
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);

		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<(32-DEPTH)>>(32-DEPTH);
			val=val<<1^val>>31;
			packsign[k]=val<<RSHIFT;
		}
		for(int ks=0;ks<1<<DEPTH;++ks)
		{
			for(int kb=0;kb<8;++kb)
			{
				uint32_t code;
				int nbypass, codelen, nzeros;

				nzeros=ks>>kb;
				codelen=nzeros+1+kb;
				code=ks;
				nbypass=kb^31;
				code<<=nbypass;
				code|=1<<31;
				code>>=nbypass;
				if(nzeros>RLIMIT-1)
					code=ks, codelen=RLIMIT+DEPTH;
				enctable[ks<<3|kb]=code<<8|codelen;
			}
		}
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			for(int kx=0;kx<iw;++kx)
			{
				uint64_t code[3];
				uint32_t codelen[3], reg;

				memcpy(&reg, rdptr, sizeof(reg));
				if(rdptr>=rdend)
				{
					uint64_t d2;

					fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
					rdptr-=BUFSIZE;
					memcpy(&d2, rdptr, sizeof(uint64_t));
					reg|=d2;
				}
				rdptr+=3;
				yuv[0]=(uint8_t)(reg>>0*8);
				yuv[1]=(uint8_t)(reg>>1*8);
				yuv[2]=(uint8_t)(reg>>2*8);
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
				yuv[1]+=(yuv[0]+yuv[2])>>2;
				yuv[2]-=yuv[0]>>2;
				sym[0]=packsignptr[yuv[0]-((rptr[0+(0+0*NCH)*NROWS*NVAL]+rptr[0+(0-1*NCH)*NROWS*NVAL])>>1)];
				sym[1]=packsignptr[yuv[1]-((rptr[0+(1+0*NCH)*NROWS*NVAL]+rptr[0+(1-1*NCH)*NROWS*NVAL])>>1)];
				sym[2]=packsignptr[yuv[2]-((rptr[0+(2+0*NCH)*NROWS*NVAL]+rptr[0+(2-1*NCH)*NROWS*NVAL])>>1)];
				code[0]=enctable[sym[0]+log2table[rptr[1+(0-1*NCH)*NROWS*NVAL]]];
				code[1]=enctable[sym[1]+log2table[rptr[1+(1-1*NCH)*NROWS*NVAL]]];
				code[2]=enctable[sym[2]+log2table[rptr[1+(2-1*NCH)*NROWS*NVAL]]];
				rptr[0+(0+0*NCH)*NROWS*NVAL]=yuv[0];
				rptr[0+(1+0*NCH)*NROWS*NVAL]=yuv[1];
				rptr[0+(2+0*NCH)*NROWS*NVAL]=yuv[2];
				rptr[1+(0+0*NCH)*NROWS*NVAL]=(2*rptr[1+(0-1*NCH)*NROWS*NVAL]+sym[0]+rptr[1+(0+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=(2*rptr[1+(1-1*NCH)*NROWS*NVAL]+sym[1]+rptr[1+(1+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=(2*rptr[1+(2-1*NCH)*NROWS*NVAL]+sym[2]+rptr[1+(2+3*NCH)*NROWS*NVAL])>>2;
#ifdef FIFOVAL
				sym[0]>>=3;
				sym[1]>>=3;
				sym[2]>>=3;
				fifoval_enqueue(sym[2]<<18^sym[1]<<9^sym[0]);
#endif
				rptr+=NCH*NROWS*NVAL;
				codelen[0]=(uint8_t)code[0];
				codelen[1]=(uint8_t)code[1];
				codelen[2]=(uint8_t)code[2];
				code[0]>>=8;
				code[1]>>=8;
				code[2]>>=8;
				{
					uint8_t s0=64-codelen[0];
					uint8_t s1=64-codelen[0]-codelen[1];
					uint8_t s2=64-codelen[0]-codelen[1]-codelen[2];

					code[0]=
						code[0]<<s0
					|	code[1]<<s1
					|	code[2]<<s2
					;
					codelen[2]+=codelen[0]+codelen[1];
				}
				codelen[2]+=nbits;
				cache|=code[0]>>nbits;
				if(codelen[2]>=64)
				{
					memcpy(wtptr, &cache, sizeof(uint64_t));
					if(wtptr>=wtend)
					{
						fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
						wtptr-=BUFSIZE;
						memcpy(wtptr, &cache, sizeof(uint64_t));
					}
					wtptr+=sizeof(uint64_t);
					cache=code[0]<<(64-nbits);
					codelen[2]-=64;
					if(!nbits)
						cache=0;
				}
				nbits=codelen[2];
			}
		}
		memcpy(wtptr, &cache, sizeof(uint64_t));
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			memcpy(wtptr, &cache, sizeof(uint64_t));
		}
		wtptr+=sizeof(uint64_t);
		(void)rdend;
		(void)wtend;
	}
	else//dec
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-3;
		uint64_t cache2=0;

		for(int k=0;k<512;++k)
		{
			int val=k;
			val=val>>1^-(val&1);
		//	val=val<<(32-DEPTH)>>(32-DEPTH);
			packsign[k]=val;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fread(&cache, 1, 8, fsrc);
		fread(&cache2, 1, 8, fsrc);
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			for(int kx=0;kx<iw;++kx)
			{
				uint64_t code;
				uint32_t reg;
				int nzeros, prefix, nbypass[3];

				nbypass[0]=log2table[rptr[1+(0-1*NCH)*NROWS*NVAL]];
				nbypass[1]=log2table[rptr[1+(1-1*NCH)*NROWS*NVAL]];
				nbypass[2]=log2table[rptr[1+(2-1*NCH)*NROWS*NVAL]];
				if(nbits>=64)
				{
					cache=cache2;
					memcpy(&cache2, rdptr, sizeof(uint64_t));
					if(rdptr>=rdend)
					{
						uint64_t d2;

						fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						memcpy(&d2, rdptr, sizeof(uint64_t));
						cache2|=d2;
					}
					rdptr+=8;
					nbits-=64;
				}
				code=cache<<nbits;
				if(nbits)
					code|=cache2>>(64-nbits);

				nzeros=(int)LZCNT64(code);
				prefix=nzeros<<nbypass[0];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[0]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[0]=(int)(code>>(64-nbypass[0]));
				if(!nbypass[0])
					sym[0]=0;
				code<<=nbypass[0];
				sym[0]|=prefix;
				nbits+=nzeros+1+nbypass[0];

				nzeros=(int)LZCNT64(code);
				prefix=nzeros<<nbypass[1];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[1]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[1]=(int)(code>>(64-nbypass[1]));
				if(!nbypass[1])
					sym[1]=0;
				code<<=nbypass[1];
				sym[1]|=prefix;
				nbits+=nzeros+1+nbypass[1];

				nzeros=(int)LZCNT64(code);
				prefix=nzeros<<nbypass[2];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[2]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[2]=(int)(code>>(64-nbypass[2]));
				if(!nbypass[2])
					sym[2]=0;
			//	code<<=nbypass[2];
				sym[2]|=prefix;
#ifdef FIFOVAL
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
				{
					printf("%016llX\n", cache);
					printf("%016llX\n", cache2);
					CRASH("");
				}
#endif
				nbits+=nzeros+1+nbypass[2];
				rptr[1+(0+0*NCH)*NROWS*NVAL]=(2*rptr[1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rptr[1+(0+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=(2*rptr[1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rptr[1+(1+3*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=(2*rptr[1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rptr[1+(2+3*NCH)*NROWS*NVAL])>>2;
				sym[0]=packsign[sym[0]];
				sym[1]=packsign[sym[1]];
				sym[2]=packsign[sym[2]];
				sym[0]+=(rptr[0+(0+0*NCH)*NROWS*NVAL]+rptr[0+(0-1*NCH)*NROWS*NVAL])>>1;
				sym[1]+=(rptr[0+(1+0*NCH)*NROWS*NVAL]+rptr[0+(1-1*NCH)*NROWS*NVAL])>>1;
				sym[2]+=(rptr[0+(2+0*NCH)*NROWS*NVAL]+rptr[0+(2-1*NCH)*NROWS*NVAL])>>1;
				sym[0]<<=32-DEPTH;
				sym[1]<<=32-DEPTH;
				sym[2]<<=32-DEPTH;
				sym[0]>>=32-DEPTH;
				sym[1]>>=32-DEPTH;
				sym[2]>>=32-DEPTH;
				rptr[0+(0+0*NCH)*NROWS*NVAL]=sym[0];
				rptr[0+(1+0*NCH)*NROWS*NVAL]=sym[1];
				rptr[0+(2+0*NCH)*NROWS*NVAL]=sym[2];
				rptr+=NCH*NROWS*NVAL;
				sym[2]+=sym[0]>>2;
				sym[1]-=(sym[0]+sym[2])>>2;
				sym[2]+=sym[1];
				sym[0]+=sym[1];
				reg=sym[2]<<16|sym[1]<<8|sym[0];
				memcpy(wtptr, &reg, sizeof(reg));
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					memcpy(wtptr, &reg, sizeof(uint64_t));
				}
				wtptr+=3;
			}
		}
		(void)rdend;
		(void)wtend;
	}
#endif
	if(wtptr>wtbuf+sizeof(uint64_t))
		fwrite(wtbuf+sizeof(uint64_t), 1, wtptr-(wtbuf+sizeof(uint64_t)), fdst);
	fclose(fsrc);
	fclose(fdst);
	free(pixels);
#ifdef LOUD
	t=time_sec2()-t;
	if(fwd)
	{
		int64_t csize=0;
		struct stat info={0};
		
#ifdef RICE2ESTIM
		printf("%12.2lf\n", ricesize/8.);
		printf("%12.2lf\n", rice2size/8.);
#endif
		stat(dstfn, &info);
		csize=info.st_size;
		printf("3*%7d*%7d  \"%s\"\n", iw, ih, srcfn);
		printf("%12lld->%12lld  %12.6lf:1\n", usize, csize, (double)usize/csize);
	}
	printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
#endif
	(void)usize;
	(void)&time_sec2;
	(void)packsignptr;
	return 0;
}
