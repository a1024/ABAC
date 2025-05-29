#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>
#ifdef _MSC_VER
#include<intrin.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#elif defined __GNUC__
#include<x86intrin.h>
#include<time.h>
#endif
#include<emmintrin.h>


#ifdef _MSC_VER
	#define LOUD
#ifdef _DEBUG
	#define ESTENT
	#define PROBHIST
//	#define PRINT_PROBTABLE
	#define ESTIMATE_BINSIZES
//	#define PRINT_RCT
	#define ENABLE_GUIDE
	#define PRINTBITS1 1000
#endif
#endif


#define NCOEFFS 15
#define COEFF_SH 18

#define NCTX 12

#define NPREDS 8

#define PROBBITS_STORE	12
#define PROBBITS_USE	12	//paq8l stretch/squash requires 12 bit


#ifdef _MSC_VER
#define AWM_INLINE __forceinline static
#else
#define AWM_INLINE __attribute__((always_inline)) inline static
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
AWM_INLINE int32_t floor_log2(uint32_t n)
{
#ifdef _MSC_VER
	unsigned long index=0;
	_BitScanReverse(&index, n);//3 cycles
	return n?index:-1;
#elif defined __GNUC__
	return n?31-__builtin_clz(n):-1;
#else
	if(!n)
		return -1;
	if(n>0x7FFFFFFF)
		return 31;
	{
		//9 cycles excluding CMOVs above
		__m128i t0=_mm_castpd_si128(_mm_cvtsi32_sd(_mm_setzero_pd(), n));
		t0=_mm_srli_epi64(t0, 52);
		return _mm_cvtsi128_si32(t0)-1023;
	}
#endif
}
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
		return;
#endif
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define FILLMEM(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
static double time_sec(void)
{
#ifdef _MSC_VER
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
static void crash(const char *format, ...)
{
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
#ifdef _MSC_VER
	system("pause");
#endif
	exit(1);
}
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		crash("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		crash("Guide error  XY %d %d", kx, ky);
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
static unsigned char* load_file(const char *fn, ptrdiff_t *ret_size)
{
	struct stat info={0};
	int error;
	ptrdiff_t size, nread;
	FILE *fsrc;
	unsigned char *buf;

	error=stat(fn, &info);
	if(error)
	{
		printf("Cannot stat \"%s\"\n", fn);
		return 0;
	}
	size=info.st_size;
	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"\n", fn);
		return 0;
	}
	buf=(unsigned char*)malloc(size+16);
	if(!buf)
	{
		crash("Alloc error\n");
		return 0;
	}
	nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
		printf("Error reading \"%s\"\n", fn);
	fclose(fsrc);
	*ret_size=size;
	return buf;
}
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
//	OCH_R2,
//	OCH_G2,
//	OCH_B2,

	OCH_COUNT,

	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,

} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},// 0
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},// 1
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},// 2
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},// 3
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},// 4
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},// 5
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},// 6
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},// 7
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},// 8
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},// 9
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},//10
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},//11
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},//12
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},//13
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},//14
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},//15
};
static int32_t stats[3][8][256][256];
static int32_t
	mixer11[3][8],		//o0
	mixer12[3][3*3*3*3][8],	//quantize(ctx)
	mixer13[3][256][8],	//[tidx]
	mixer14[3][256][8][8],	//[N+W][kb]
	mixer2[3][4];		//level 2
#ifdef ESTIMATE_BINSIZES
static double shannon_table[1<<PROBBITS_USE];
static double binsizes[6]={0};
#endif
static uint32_t probtable[1<<PROBBITS_USE];
#ifdef PROBHIST
static int32_t probhist[1<<PROBBITS_USE];
#endif
#ifdef ESTENT
static int32_t estent[3][8][256];
#endif
//static int16_t stretchtable[4096];//paq8l
AWM_INLINE int32_t squash2(int32_t p0)
{
	CLAMP2(p0, -(1<<PROBBITS_USE>>1), (1<<PROBBITS_USE>>1)-1);
	p0+=1<<PROBBITS_USE>>1;
	p0=probtable[p0];

	return p0;
}
AWM_INLINE int32_t squash(int32_t p0)
{
	/*
	1/(1+exp(-13(x-0.5)))			//sigmoid
	x>0.5?1-8*sq(sq(1-x)):8*sq(sq(x))	//good approximation
	x>0.5?1-2*sq(1-x):2sq(x)		//smooth
	(x-0.5-re((x-0.5)^3))*4/3+0.5		//almost linear

	x>0?0.5*(1-sq(sq(1-2*x))):-0.5*(1-sq(sq(1+2*x)))		//signed->signed version

	1/(1+exp(-13(x)))-0.5
	0.5*sgn(x)*(1-sq(sq(1-2*abs(x))))	//branchless
	*/
	int32_t negmask;

	CLAMP2(p0, -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1);
	negmask=p0>>31;
	p0<<=1;
	p0^=negmask;
	p0-=negmask;
	p0=(1<<PROBBITS_STORE)-p0;
	p0=(int32_t)((int64_t)p0*p0>>PROBBITS_STORE);
	p0=(int32_t)((int64_t)p0*p0>>PROBBITS_STORE);
	p0=(1<<PROBBITS_STORE)-p0;
	p0^=negmask;
	p0-=negmask;
	p0>>=1;

	return p0;
}
AWM_INLINE int squash_paq8l(int d)
{
	static const int t[33]=
	{
		1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
		1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
		4050,4068,4079,4085,4089,4092,4093,4094,
	};
	int w;

	if(d>+2047)return 4095;
	if(d<-2047)return 1;
	w=d&127;
	d=(d>>7)+16;
	return t[d]+(((t[d+1]-t[d])*w+64)>>7);
}
int s04_codec(const char *command, const char *srcfn, const char *dstfn)
{
	//const int mem=sizeof(stats)+sizeof(mixer11)+sizeof(mixer12)+sizeof(mixer13)+sizeof(mixer14)+sizeof(mixer2);
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *streamptr=0, *streamend=0;
	unsigned char *imptr=0;
	int bestrct=0;
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef PRINTBITS1
	uint32_t bitidx=0;
#endif
	uint64_t low=0, range=0xFFFFFFFFFFFFFFFF, code=0;

	(void)memfill;
	//{
	//	int k, pi;
	//
	//	for(k=-2048, pi=0;k<2047;++k)
	//	{
	//		int i=squash_paq8l(k), k2;
	//		for(k2=pi;k2<=i;++k2)
	//			stretchtable[k2]=k;
	//		pi=i+1;
	//	}
	//	stretchtable[4095]=2047;
	//}
#ifdef ESTIMATE_BINSIZES
	{
		int32_t kp;

		//for(kp=0;kp<(1<<PROBBITS_USE);++kp)
		//{
		//	int prob=kp;
		//	CLAMP2(prob, 1, (1<<PROBBITS_USE)-1);
		//	Shannon[kp]=(int32_t)(-log((double)prob*(1./(1<<PROBBITS_USE)))*((1<<16)/M_LN2));
		//}
		for(kp=0;kp<(1<<PROBBITS_USE);++kp)
		{
			int prob=kp;
			CLAMP2(prob, 1, (1<<PROBBITS_USE)-1);
			shannon_table[kp]=-log((double)prob*(1./(1<<PROBBITS_USE)))*(1./(8*M_LN2));
		}
	}
#endif
	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		crash("Cannot open \"%s\"", srcfn);
		return 1;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'4'<<8))
		{
			crash("Unsupported file  \"%s\"\n", srcfn);
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	if(fwd)
	{
		if(*srcptr++!='\n')
		{
			crash("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		iw=0;
		while((uint32_t)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr<=' ')++srcptr;//skip space
		ih=0;
		while((uint32_t)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(*srcptr++!='\n')
		{
			crash("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		if(memcmp(srcptr, "255\n", 4))
		{
			crash("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		srcptr+=4;
		image=srcptr;
		guide_save(image, iw, ih);
	}
	else
	{
		iw=((uint32_t*)srcptr)[0];
		ih=((uint32_t*)srcptr)[1];
		srcptr+=sizeof(uint32_t[2]);
	}
	if(iw<1||ih<1)
	{
		crash("Invalid file\n");
		free(srcbuf);
		return 1;
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	dstbuf=(unsigned char*)malloc(usize);
	if(!dstbuf)
	{
		crash("Alloc error");
		free(srcbuf);
		return 1;
	}
	if(fwd)
	{
		//analysis
		int64_t counters[OCH_COUNT]={0}, minerr=0;
		int prev[OCH_COUNT]={0};

		imptr=image;
		while(imptr<srcend)
		{
			int r, g, b, rg, gb, br;

			r=imptr[0];
			g=imptr[1];
			b=imptr[2];
			imptr+=3;
			rg=r-g;
			gb=g-b;
			br=b-r;
			counters[0]+=abs(r	-prev[0]);
			counters[1]+=abs(g	-prev[1]);
			counters[2]+=abs(b	-prev[2]);
			counters[3]+=abs(rg	-prev[3]);
			counters[4]+=abs(gb	-prev[4]);
			counters[5]+=abs(br	-prev[5]);
			prev[0]=r;
			prev[1]=g;
			prev[2]=b;
			prev[3]=rg;
			prev[4]=gb;
			prev[5]=br;
		}
		imptr=image;
		{
			int kt;

#ifdef PRINT_RCT
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			for(kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
				long long currerr=
					+counters[rct[0]]
					+counters[rct[1]]
					+counters[rct[2]]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
#ifdef PRINT_RCT
				printf("RCT%02d %16lld%s\n", kt, currerr, kt==bestrct?" <-":"");
#endif
			}
		}
		streamend=dstbuf+usize;
		streamptr=dstbuf;
		*streamptr++=bestrct;
	}
	else
	{
		imptr=dstbuf;
		streamend=srcptr+srcsize;
		streamptr=srcptr;

		bestrct=*streamptr++;

		code=0;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;

		csize=srcsize;
	}
#if 1
	{
		int k;

		for(k=0;k<(1<<PROBBITS_USE);++k)
		{
			int p0e=k;

			p0e+=p0e<(1<<PROBBITS_USE>>1);
			p0e-=1<<PROBBITS_USE>>1;
		//	p0e+=p0e>>2;
		//	p0e+=p0e>>3;
		//	p0e<<=1;
			CLAMP2(p0e, -(1<<PROBBITS_USE>>1)+1, (1<<PROBBITS_USE>>1)-1);
			{
				int32_t negmask=p0e>>31;

				p0e<<=1;
				p0e^=negmask;
				p0e-=negmask;
				p0e=(1<<PROBBITS_USE)-p0e;
				p0e=(int32_t)((int64_t)p0e*p0e>>PROBBITS_USE);
				p0e=(int32_t)((int64_t)p0e*p0e>>PROBBITS_USE);
				p0e=(1<<PROBBITS_USE)-p0e;
				p0e^=negmask;
				p0e-=negmask;
				p0e>>=1;

				CLAMP2(p0e, -(1<<PROBBITS_USE>>1)+1, (1<<PROBBITS_USE>>1)-1);
			//	p0e+=1<<PROBBITS_USE>>1;
			//	CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
			}
			probtable[k]=p0e;
#ifdef PRINT_PROBTABLE
			if(fwd)
			{
				int k2;

				printf("%8d  ", k);
				for(k2=0;k2<(p0e>>3);++k2)
					printf("*");
				printf("\n");
			}
#endif
		}
	}
#endif
	{
		int
			yidx=rct_indices[bestrct][3+0],
			uidx=rct_indices[bestrct][3+1],
			vidx=rct_indices[bestrct][3+2],
			uhelpidx=rct_indices[bestrct][6+0],
			vhelpidx=rct_indices[bestrct][6+1];
		int32_t ky, kx, idx;
		int32_t psize=0;
		int32_t *pixels=0;
		int32_t paddedwidth=iw+16;

		psize=(int32_t)sizeof(int32_t[4*3*2])*paddedwidth;//4 padded rows * 3 channels * {pixels, nbypass}
		pixels=(int32_t*)malloc(psize);
		if(!pixels)
		{
			crash("Alloc error\n");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		memset(stats, 0, sizeof(stats));
		memset(mixer11, 0, sizeof(mixer11));
		memset(mixer12, 0, sizeof(mixer12));
		memset(mixer13, 0, sizeof(mixer13));
		memset(mixer14, 0, sizeof(mixer14));
		memset(mixer2, 0, sizeof(mixer2));
#ifdef PROBHIST
		memset(probhist, 0, sizeof(probhist));
#endif
#ifdef ESTENT
		memset(estent, 0, sizeof(estent));
#endif
		for(ky=0, idx=0;ky<ih;++ky)
		{
			unsigned char yuv[4]={0};
			int32_t *rows[]=
			{
				pixels+(paddedwidth*((ky-0)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-1)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-2)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-3)&3)+8)*3*2,
			};
			for(kx=0;kx<iw;++kx, ++idx)
			{
				int kc;
				int offset;

				if(fwd)
				{
					yuv[0]=imptr[yidx];
					yuv[1]=imptr[uidx];
					yuv[2]=imptr[vidx];
				}
				offset=0;
				for(kc=0;kc<3;++kc)
				{
					int32_t
						NNN	=rows[3][0+0*3*2],
						NNNE	=rows[3][0+1*3*2],
						NN	=rows[2][0+0*3*2],
						NNE	=rows[2][0+1*3*2],
						NWW	=rows[1][0-2*3*2],
						NW	=rows[1][0-1*3*2],
						N	=rows[1][0+0*3*2],
						NE	=rows[1][0+1*3*2],
						NEE	=rows[1][0+2*3*2],
						NEEE	=rows[1][0+3*3*2],
						NEEEE	=rows[1][0+4*3*2],
						WWWWW	=rows[0][0-5*3*2],
						WWWW	=rows[0][0-4*3*2],
						WWW	=rows[0][0-3*3*2],
						WW	=rows[0][0-2*3*2],
						W	=rows[0][0-1*3*2];
					int32_t error;
					int32_t qctx, qctx2;
					int32_t kb, tidx, bit;
					int32_t *statsptr[8];
					int32_t ctxs[]=
					{
						N,
						W,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
						N+W-NW,

					//	(4*(N+W)-NN-WW)/6,
						W+NE-N,

					//	(4*(N+W)+NE-NW)>>3,
					//	(N+W)>>1,
					//	(N+NE-NNE+W+NW-NWW)>>1,
						N+NE-NNE,

						(15*W-20*WW+15*WWW-6*WWWW+WWWWW+3*(NE-NNE)+NNNE)/6,
					//	(WWWW+NNN+NEEE+NEEEE)/4,
					//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
					};
					ctxs[0]+=offset;
					ctxs[1]+=offset;
					ctxs[2]+=offset;
					ctxs[3]+=offset;
					ctxs[4]+=offset;
					ctxs[5]+=offset;
					ctxs[6]+=offset;
					ctxs[7]+=offset;

					CLAMP2(ctxs[0], 0, 255);
					CLAMP2(ctxs[1], 0, 255);
					CLAMP2(ctxs[2], 0, 255);
					CLAMP2(ctxs[3], 0, 255);
					CLAMP2(ctxs[4], 0, 255);
					CLAMP2(ctxs[5], 0, 255);
					CLAMP2(ctxs[6], 0, 255);
					CLAMP2(ctxs[7], 0, 255);

					qctx=(N>W)-(N<W)+1;
					qctx=3*qctx+(N>NW)-(N<NW)+1;
					qctx=3*qctx+(W>NW)-(W<NW)+1;
					qctx=3*qctx+(NE>N)-(NE<N)+1;

					qctx2=((N+W)>>1)+offset;
					CLAMP2(qctx2, 0, 255);
					
					statsptr[0]=stats[kc][0][(unsigned char)ctxs[0]];
					statsptr[1]=stats[kc][1][(unsigned char)ctxs[1]];
					statsptr[2]=stats[kc][2][(unsigned char)ctxs[2]];
					statsptr[3]=stats[kc][3][(unsigned char)ctxs[3]];
					statsptr[4]=stats[kc][4][(unsigned char)ctxs[4]];
					statsptr[5]=stats[kc][5][(unsigned char)ctxs[5]];
					statsptr[6]=stats[kc][6][(unsigned char)ctxs[6]];
					statsptr[7]=stats[kc][7][(unsigned char)ctxs[7]];
					if(fwd)
					{
						error=yuv[kc];
#ifdef ESTENT
						++estent[kc][0][(unsigned char)(error-ctxs[0]+128)];
						++estent[kc][1][(unsigned char)(error-ctxs[1]+128)];
						++estent[kc][2][(unsigned char)(error-ctxs[2]+128)];
						++estent[kc][3][(unsigned char)(error-ctxs[3]+128)];
						++estent[kc][4][(unsigned char)(error-ctxs[4]+128)];
						++estent[kc][5][(unsigned char)(error-ctxs[5]+128)];
						++estent[kc][6][(unsigned char)(error-ctxs[6]+128)];
						++estent[kc][7][(unsigned char)(error-ctxs[7]+128)];
#endif
					}
					else
						error=0;

					//if(ky==13&&kx==785)//
					if(ky==10&&kx==iw/2)//
						printf("");

					for(kb=7, tidx=1;kb>=0;--kb)
					{
						int32_t *pp0[]=
						{
							statsptr[0]+tidx,
							statsptr[1]+tidx,
							statsptr[2]+tidx,
							statsptr[3]+tidx,
							statsptr[4]+tidx,
							statsptr[5]+tidx,
							statsptr[6]+tidx,
							statsptr[7]+tidx,
						};
						int32_t p00[]=
						{
							*pp0[0],
							*pp0[1],
							*pp0[2],
							*pp0[3],
							*pp0[4],
							*pp0[5],
							*pp0[6],
							*pp0[7],
						};
						int32_t *currmixer[]=
						{
							mixer11[kc],
							mixer12[kc][qctx],
							mixer13[kc][tidx-1],
							mixer14[kc][qctx2][kb],
							mixer2[kc],
						};
						int64_t logits[4], p0;
						int32_t p0e;

						logits[0]=(int64_t)currmixer[0][0]*p00[0];
						logits[1]=(int64_t)currmixer[1][0]*p00[0];
						logits[2]=(int64_t)currmixer[2][0]*p00[0];
						logits[3]=(int64_t)currmixer[3][0]*p00[0];
						logits[0]+=(int64_t)currmixer[0][1]*p00[1];
						logits[1]+=(int64_t)currmixer[1][1]*p00[1];
						logits[2]+=(int64_t)currmixer[2][1]*p00[1];
						logits[3]+=(int64_t)currmixer[3][1]*p00[1];
						logits[0]+=(int64_t)currmixer[0][2]*p00[2];
						logits[1]+=(int64_t)currmixer[1][2]*p00[2];
						logits[2]+=(int64_t)currmixer[2][2]*p00[2];
						logits[3]+=(int64_t)currmixer[3][2]*p00[2];
						logits[0]+=(int64_t)currmixer[0][3]*p00[3];
						logits[1]+=(int64_t)currmixer[1][3]*p00[3];
						logits[2]+=(int64_t)currmixer[2][3]*p00[3];
						logits[3]+=(int64_t)currmixer[3][3]*p00[3];
						logits[0]+=(int64_t)currmixer[0][4]*p00[4];
						logits[1]+=(int64_t)currmixer[1][4]*p00[4];
						logits[2]+=(int64_t)currmixer[2][4]*p00[4];
						logits[3]+=(int64_t)currmixer[3][4]*p00[4];
						logits[0]+=(int64_t)currmixer[0][5]*p00[5];
						logits[1]+=(int64_t)currmixer[1][5]*p00[5];
						logits[2]+=(int64_t)currmixer[2][5]*p00[5];
						logits[3]+=(int64_t)currmixer[3][5]*p00[5];
						logits[0]+=(int64_t)currmixer[0][6]*p00[6];
						logits[1]+=(int64_t)currmixer[1][6]*p00[6];
						logits[2]+=(int64_t)currmixer[2][6]*p00[6];
						logits[3]+=(int64_t)currmixer[3][6]*p00[6];
						logits[0]+=(int64_t)currmixer[0][7]*p00[7];
						logits[1]+=(int64_t)currmixer[1][7]*p00[7];
						logits[2]+=(int64_t)currmixer[2][7]*p00[7];
						logits[3]+=(int64_t)currmixer[3][7]*p00[7];
						logits[0]=squash2((int32_t)(logits[0]>>24));
						logits[1]=squash2((int32_t)(logits[1]>>20));
						logits[2]=squash2((int32_t)(logits[2]>>20));
						logits[3]=squash2((int32_t)(logits[3]>>16));
						//logits[0]=squash2((int32_t)(logits[0]>>36));
						//logits[1]=squash2((int32_t)(logits[1]>>35));
						//logits[2]=squash2((int32_t)(logits[2]>>35));
						//logits[3]=squash2((int32_t)(logits[3]>>29));
						p0=(
							+(int64_t)currmixer[4][0]*logits[0]
							+(int64_t)currmixer[4][1]*logits[1]
							+(int64_t)currmixer[4][2]*logits[2]
							+(int64_t)currmixer[4][3]*logits[3]
						)>>15;
						p0=squash2((int32_t)p0);
						p0e=(int32_t)((p0+(1<<PROBBITS_STORE>>1))>>(PROBBITS_STORE-PROBBITS_USE));
						CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
#if 0
#define SH 9
						int32_t p0=(
							+p00[0]
							+p00[1]
							+p00[2]
							+p00[3]
							+p00[4]
							+p00[5]
							+p00[6]
							+p00[7]
						//	+(1<<(PROBBITS_USE+SH)>>1)
							+(1<<SH>>1)
						)>>SH;
						p0=squash_paq8l(p0);
						//p0+=1<<PROBBITS_USE>>1;
						//CLAMP2(p0, 1, (1<<PROBBITS_USE)-1);
						//p0=probtable[p0];
#endif
#ifdef PROBHIST
						++probhist[p0e];
#endif
						if(range<(1<<PROBBITS_USE))
						{
							if(streamptr>=streamend)
							{
								int symidx=3*idx+kc, totalsyms=3*iw*ih;

								crash("ERROR at %d/%d  inflation %8.4lf%%\n"
									, (int32_t)symidx
									, (int32_t)totalsyms
									, 100.*totalsyms/symidx
								);
								return 1;
							}
							if(fwd)
								*(uint32_t*)streamptr=(uint32_t)(low>>32);
							else
								code=code<<32|*(uint32_t*)streamptr;
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						bit=error>>kb&1;
						{
							uint64_t r2=(range>>PROBBITS_USE)*p0e;
							uint64_t mid=low+r2;
							range-=r2;
							if(!fwd)
								bit=code>=mid;
							if(bit)
								low=mid;
							else
								range=r2-1;
						}
						error|=bit<<kb;
						{
							int32_t error=(!bit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1)-(int32_t)p0, k2;
							//error=(error>0)-(error<0);
							for(k2=0;k2<4;++k2)
							{
								int64_t t0, t1, t2, t3, t4, t5, t6, t7;

								t0=currmixer[k2][0]+(((int64_t)error*p00[0]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t1=currmixer[k2][1]+(((int64_t)error*p00[1]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t2=currmixer[k2][2]+(((int64_t)error*p00[2]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t3=currmixer[k2][3]+(((int64_t)error*p00[3]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t4=currmixer[k2][4]+(((int64_t)error*p00[4]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t5=currmixer[k2][5]+(((int64_t)error*p00[5]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t6=currmixer[k2][6]+(((int64_t)error*p00[6]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t7=currmixer[k2][7]+(((int64_t)error*p00[7]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								//CLAMP2(t0, -0x8000, 0x7FFF);
								//CLAMP2(t1, -0x8000, 0x7FFF);
								//CLAMP2(t2, -0x8000, 0x7FFF);
								//CLAMP2(t3, -0x8000, 0x7FFF);
								//CLAMP2(t4, -0x8000, 0x7FFF);
								//CLAMP2(t5, -0x8000, 0x7FFF);
								//CLAMP2(t6, -0x8000, 0x7FFF);
								//CLAMP2(t7, -0x8000, 0x7FFF);
								currmixer[k2][0]=(int32_t)t0;
								currmixer[k2][1]=(int32_t)t1;
								currmixer[k2][2]=(int32_t)t2;
								currmixer[k2][3]=(int32_t)t3;
								currmixer[k2][4]=(int32_t)t4;
								currmixer[k2][5]=(int32_t)t5;
								currmixer[k2][6]=(int32_t)t6;
								currmixer[k2][7]=(int32_t)t7;
							}
							{
								int64_t t0, t1, t2, t3;

								t0=currmixer[4][0]+(((int64_t)error*logits[0]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t1=currmixer[4][1]+(((int64_t)error*logits[1]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t2=currmixer[4][2]+(((int64_t)error*logits[2]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								t3=currmixer[4][3]+(((int64_t)error*logits[3]+(1<<(PROBBITS_STORE+6)>>1))>>(PROBBITS_STORE+6));
								//CLAMP2(t0, -0x8000, 0x7FFF);
								//CLAMP2(t1, -0x8000, 0x7FFF);
								//CLAMP2(t2, -0x8000, 0x7FFF);
								//CLAMP2(t3, -0x8000, 0x7FFF);
								currmixer[4][0]=(int32_t)t0;
								currmixer[4][1]=(int32_t)t1;
								currmixer[4][2]=(int32_t)t2;
								currmixer[4][3]=(int32_t)t3;
							}
						}
#if 1
						{
							int32_t e;

							e=(!bit<<PROBBITS_USE)+(1<<7>>1);
							*pp0[0]=p00[0]+((e-p00[0])>>7);
							*pp0[1]=p00[1]+((e-p00[1])>>7);
							*pp0[2]=p00[2]+((e-p00[2])>>7);
							*pp0[3]=p00[3]+((e-p00[3])>>7);
							*pp0[4]=p00[4]+((e-p00[4])>>7);
							*pp0[5]=p00[5]+((e-p00[5])>>7);
							*pp0[6]=p00[6]+((e-p00[6])>>7);
							*pp0[7]=p00[7]+((e-p00[7])>>7);

							//e=(!bit<<PROBBITS_USE)-p0e;
							//*pp0[0]=p00[0]+e;
							//*pp0[1]=p00[1]+e;
							//*pp0[2]=p00[2]+e;
							//*pp0[3]=p00[3]+e;
							//*pp0[4]=p00[4]+e;
							//*pp0[5]=p00[5]+e;
							//*pp0[6]=p00[6]+e;
							//*pp0[7]=p00[7]+e;
						}
#endif
#ifdef ESTIMATE_BINSIZES
						binsizes[kc]+=shannon_table[bit?(1<<PROBBITS_USE)-p0e:p0e];
#endif
#ifdef PRINTBITS1
						{
							const int checkpoint=200000;
							if(fwd&&(uint32_t)(bitidx-checkpoint)<PRINTBITS1)printf("%c", '0'+bit);
							if(fwd&&bitidx==checkpoint+PRINTBITS1)printf("\n\n");
							++bitidx;
						}
#endif
						tidx=2*tidx+bit;
					}
					if(!fwd)
						yuv[kc]=(char)error;
					rows[0][0]=yuv[kc]-offset;
					offset=kc?yuv[vhelpidx]:yuv[uhelpidx];
					rows[0]+=2;
					rows[1]+=2;
					rows[2]+=2;
					rows[3]+=2;
				}
				if(!fwd)
				{
					imptr[yidx]=yuv[0];
					imptr[uidx]=yuv[1];
					imptr[vidx]=yuv[2];
					guide_check(dstbuf, kx, ky);
				}
				imptr+=3;
			}
		}
		free(pixels);
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			crash("Cannot open \"%s\" for writing\n", dstfn);
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		if(fwd)
		{
			*(uint32_t*)streamptr=(uint32_t)(low>>32);
			streamptr+=4;
			*(uint32_t*)streamptr=(uint32_t)low;
			streamptr+=4;

			csize=streamptr-dstbuf;

			dstsize+=fwrite("04", 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 4, fdst);
			dstsize+=fwrite(&ih, 1, 4, fdst);
			dstsize+=fwrite(dstbuf, 1, csize, fdst);
		}
		else
		{
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(srcbuf);
	free(dstbuf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef PRINTBITS1
		printf("\n%d\n", bitidx);
#endif
#ifdef PROBHIST
		{
			int k;

			k=0;
			if((1<<PROBBITS_USE)>1000)
				k=(1<<PROBBITS_USE)-1000;
			for(;k<(1<<PROBBITS_USE);++k)
				printf("%5d %8d\n", k, probhist[k]);
		}
#endif
#ifdef ESTENT
		{
			double ent[3][8]={0};
			int kc, ke, ks;

			for(kc=0;kc<3;++kc)
			{
				for(ke=0;ke<8;++ke)
				{
					double e=0, norm;
					int32_t *currhist=estent[kc][ke];
					int32_t sum=0;
					for(ks=0;ks<256;++ks)
						sum+=currhist[ks];
					norm=1./sum;
					for(ks=0;ks<256;++ks)
					{
						int freq=currhist[ks];
						if(freq)
							e-=freq*log(freq*norm);
					}
					e*=1./(8*M_LN2);
					ent[kc][ke]=e;
					//printf("C%d E%d %12.2lf\n", kc, ke, e);
				}
				//printf("\n");
			}
			for(ke=0;ke<8;++ke)
				printf("E%d %12.2lf   %12.2lf %12.2lf %12.2lf\n"
					, ke
					, ent[0][ke]+ent[1][ke]+ent[2][ke]
					, ent[0][ke]
					, ent[1][ke]
					, ent[2][ke]
				);
		}
#endif
#ifdef ESTIMATE_BINSIZES
		printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf\n"
			, binsizes[0]+binsizes[1]+binsizes[2]
			, binsizes[0]
			, binsizes[1]
			, binsizes[2]
		);
#endif
		printf("%9d->%9d  %8.4lf%%  %12.6lf\n"
			, srcsize
			, dstsize
			, 100.*dstsize/srcsize
			, (double)srcsize/dstsize
		);
	}
	printf("%c  %12.6lf  %12.6lf MB/s\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
	);
#endif
	(void)time_sec;
	return 0;
}