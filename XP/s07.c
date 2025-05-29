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
//	#define PRINT_RCT
	#define ENABLE_GUIDE
//	#define ENABLE_VAL
#endif
#endif


//from libjxl		sym = packsign(delta) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of FLOOR_LOG2(sym)-E,  bypass = 0bBB...BB  FLOOR_LOG2(sym)-(M+L) bits
#define CONFIG_EXP 4	//revise AVX2 CDF search to change config
#define CONFIG_MSB 1	//410->1+24+1	411->1+32+1	421->1+48+1
#define CONFIG_LSB 1	//		511->1+44+1
#define NLEVELS 33

#define NCTX 12

#define NPREDS 8
#define L1SH 19

#define PROBBITS_STORE	24
#define PROBBITS_USE	16

#define RATE 10


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
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
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
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
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
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("Guide error  XY %d %d", kx, ky);
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef ENABLE_VAL
typedef struct _ValInfo
{
	int cdf, freq;
	unsigned long long low0, range0, low1, range1;
} ValInfo;
static uint32_t valcap=1, valcount=0, validx=0;
static ValInfo *valdata=0;
static void val_enc(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long low1, unsigned long long range1)
{
	{
		void *ptr;

		++valcount;
		if(valcount>=valcap)
		{
			valcap<<=1;
			ptr=realloc(valdata, sizeof(ValInfo)*valcap);
			if(!ptr)
			{
				CRASH("Alloc error");
				return;
			}
			valdata=(ValInfo*)ptr;
		}
	}
	{
		ValInfo info={cdf, freq, low0, range0, low1, range1};
		valdata[valcount-1]=info;
	}
}
static int val_dec(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long code0, unsigned long long low1, unsigned long long range1, unsigned long long code1)
{
	ValInfo *info=valdata+validx;
	if(info->cdf!=cdf||info->freq!=freq||info->low0!=low0||info->range0!=range0||info->low1!=low1||info->range1!=range1)
	{
		printf("OG  %d %d\n", info->cdf, info->freq);
		printf("  0: %016llX %016llX\n", info->low0, info->range0);
		printf("  1: %016llX %016llX\n", info->low1, info->range1);
		printf("X   %d %d\n", cdf, freq);
		printf("  0: %016llX %016llX %016llX\n", low0, range0, code0);
		printf("  1: %016llX %016llX %016llX\n", low1, range1, code1);
		return 1;
	}
	++validx;
	return 0;
}
#define VAL_DEC(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1, ...) if(val_dec(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1))CRASH(__VA_ARGS__)
#else
#define val_enc(...)
#define val_dec(...)
#define VAL_DEC(...)
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
		CRASH("Alloc error\n");
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
typedef struct _QuantInfoFwd
{
	uint8_t token, nbypass;
	uint16_t bypass;
} QuantInfoFwd;
typedef struct _QuantInfoInv
{
	uint8_t sym, nbits;
} QuantInfoInv;
static uint32_t stats[3][NCTX][NLEVELS+1], mixin[NLEVELS*2];
int s07_codec(const char *command, const char *srcfn, const char *dstfn)
{
	//const int mem=sizeof(stats)+sizeof(mixer11)+sizeof(mixer12)+sizeof(mixer13)+sizeof(mixer14)+sizeof(mixer2);
	unsigned char *srcbuf=0, *srcptr=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *streamptr=0, *streamend=0;
	unsigned char *imptr=0;
	int bestrct=0, NWratio=0;
#ifdef LOUD
	double t=time_sec();
#endif
	uint64_t low=0, range=0xFFFFFFFF, code=0;

	(void)memfill;

	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	srcptr=srcbuf;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'4'<<8))
		{
			CRASH("Unsupported file  \"%s\"\n", srcfn);
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	if(fwd)
	{
		if(*srcptr++!='\n')
		{
			CRASH("Unsupported file\n");
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
			CRASH("Unsupported file\n");
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
			CRASH("Unsupported file\n");
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
		CRASH("Invalid file\n");
		free(srcbuf);
		return 1;
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	dstbuf=(unsigned char*)malloc(usize);
	if(!dstbuf)
	{
		CRASH("Alloc error");
		free(srcbuf);
		return 1;
	}
	if(fwd)
	{
		//analysis
		int ky, kx;
		int64_t counters[OCH_COUNT*2]={0}, minerr=0;
		int rowstride=3*iw;

		imptr=image+rowstride;
		for(ky=1;ky<ih;++ky)
		{
			int W[OCH_COUNT]={0};
			for(kx=0;kx<iw;++kx, imptr+=3)
			{
				int N[OCH_COUNT], curr[OCH_COUNT];
				
				N[0]=imptr[0-rowstride];
				N[1]=imptr[1-rowstride];
				N[2]=imptr[2-rowstride];
				N[3]=N[0]-N[1];
				N[4]=N[1]-N[2];
				N[5]=N[2]-N[0];
				curr[0]=imptr[0];
				curr[1]=imptr[1];
				curr[2]=imptr[2];
				curr[3]=curr[0]-curr[1];
				curr[4]=curr[1]-curr[2];
				curr[5]=curr[2]-curr[0];
				counters[0+0*OCH_COUNT]+=abs(curr[0]-N[0]);
				counters[1+0*OCH_COUNT]+=abs(curr[1]-N[1]);
				counters[2+0*OCH_COUNT]+=abs(curr[2]-N[2]);
				counters[3+0*OCH_COUNT]+=abs(curr[3]-N[3]);
				counters[4+0*OCH_COUNT]+=abs(curr[4]-N[4]);
				counters[5+0*OCH_COUNT]+=abs(curr[5]-N[5]);
				counters[0+1*OCH_COUNT]+=abs(curr[0]-W[0]);
				counters[1+1*OCH_COUNT]+=abs(curr[1]-W[1]);
				counters[2+1*OCH_COUNT]+=abs(curr[2]-W[2]);
				counters[3+1*OCH_COUNT]+=abs(curr[3]-W[3]);
				counters[4+1*OCH_COUNT]+=abs(curr[4]-W[4]);
				counters[5+1*OCH_COUNT]+=abs(curr[5]-W[5]);
				W[0]=curr[0];
				W[1]=curr[1];
				W[2]=curr[2];
				W[3]=curr[3];
				W[4]=curr[4];
				W[5]=curr[5];
			}
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
				const unsigned char *rctinfo=rct_indices[kt];
				int64_t currerr=
					+counters[rctinfo[0]]+counters[rctinfo[0]+OCH_COUNT]
					+counters[rctinfo[1]]+counters[rctinfo[1]+OCH_COUNT]
					+counters[rctinfo[2]]+counters[rctinfo[2]+OCH_COUNT]
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
		{
			const unsigned char *rctinfo=rct_indices[bestrct];
			int64_t Nerror=
				+counters[rctinfo[0]]
				+counters[rctinfo[1]]
				+counters[rctinfo[2]]
			;
			int64_t Werror=
				+counters[rctinfo[0]+OCH_COUNT]
				+counters[rctinfo[1]+OCH_COUNT]
				+counters[rctinfo[2]+OCH_COUNT]
			;
			NWratio=(int32_t)(((Werror+1)<<8)/(Nerror+Werror+2));
			CLAMP2(NWratio, 1, 255);
		}

		streamend=dstbuf+usize;
		streamptr=dstbuf;

		*streamptr++=bestrct;
		*streamptr++=NWratio;
	}
	else
	{
		imptr=dstbuf;
		streamend=srcptr+srcsize;
		streamptr=srcptr;

		bestrct=*streamptr++;
		NWratio=*streamptr++;

		code=0;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;

		csize=srcsize;
	}
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
		int32_t coeffs[3][NPREDS]={0};
		QuantInfoFwd qtablefwd[256];
		QuantInfoInv qtableinv[33];
		
		if(fwd)//init qtables
		{
			int ks;

			for(ks=0;ks<256;++ks)
			{
				QuantInfoFwd *info=qtablefwd+ks;

				if(ks<(1<<CONFIG_EXP))
				{
					info->nbypass=0;
					info->token=ks;
					info->bypass=0;
				}
				else
				{
					info->nbypass=floor_log2(ks)-(CONFIG_LSB+CONFIG_MSB);
					info->token=(1<<CONFIG_EXP)-((CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB))<<(CONFIG_MSB+CONFIG_LSB)) + (
						info->nbypass<<(CONFIG_MSB+CONFIG_LSB)|
						(ks>>info->nbypass&((1<<CONFIG_MSB)-1)<<CONFIG_LSB)|
						(ks&((1<<CONFIG_LSB)-1))
					);
					info->bypass=ks>>CONFIG_LSB&((1LL<<info->nbypass)-1);
				}
			}
		}
		else
		{
			int kt;

			for(kt=0;kt<33;++kt)
			{
				QuantInfoInv *info=qtableinv+kt;
				if(kt<(1<<CONFIG_EXP))
				{
					info->sym=kt;
					info->nbits=0;
				}
				else
				{
					info->nbits=(kt>>(CONFIG_MSB+CONFIG_LSB))-((1<<(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)))-(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)));
					info->sym=
						+(((1<<(CONFIG_MSB+CONFIG_LSB))+(kt&((1<<CONFIG_MSB)-1)<<CONFIG_LSB))<<info->nbits)//7 instructions
					//	+(bypass<<CONFIG_LSB)
						+(kt&((1<<CONFIG_LSB)-1))
					;
				}
			}
		}
		psize=(int32_t)sizeof(int32_t[4*3*2])*paddedwidth;//4 padded rows * 3 channels * {pixels, nbypass}
		pixels=(int32_t*)malloc(psize);
		if(!pixels)
		{
			CRASH("Alloc error\n");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		memset(stats, 0, sizeof(stats));
		{
			int kctx, ks;

			for(kctx=0;kctx<NCTX;++kctx)
			{
				for(ks=0;ks<NLEVELS;++ks)//init bypass
					stats[0][kctx][ks]=ks*(((1<<PROBBITS_STORE)-(NLEVELS<<(PROBBITS_STORE-PROBBITS_USE)))/NLEVELS);
				stats[0][kctx][NLEVELS]=1<<PROBBITS_STORE;
			}
		}
		memcpy(stats[1], stats[0], sizeof(stats[0]));
		memcpy(stats[2], stats[0], sizeof(stats[0]));
		memset(mixin, 0, sizeof(mixin));//
		FILLMEM(
			(uint32_t*)mixin,
			(1<<RATE>>1),
			sizeof(int32_t[NLEVELS]),
			sizeof(int32_t)
		);
		FILLMEM(
			(uint32_t*)mixin+NLEVELS,
			(1<<PROBBITS_STORE)-(NLEVELS<<(PROBBITS_STORE-PROBBITS_USE))+(1<<RATE>>1),
			sizeof(int32_t[NLEVELS]),
			sizeof(int32_t)
		);
		FILLMEM((int32_t*)coeffs, (1<<L1SH)/NPREDS, sizeof(coeffs), sizeof(int32_t));
		for(ky=0, idx=0;ky<ih;++ky)
		{
			char yuv[4]={0};
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
					yuv[0]=imptr[yidx]-128;
					yuv[1]=imptr[uidx]-128;
					yuv[2]=imptr[vidx]-128;
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
						W	=rows[0][0-1*3*2],
						eNEE	=rows[1][1+2*3*2],
						eNEEE	=rows[1][1+3*3*2],
						eW	=rows[0][1-1*3*2];
					int32_t p0, pred, ctx, error, k, token, cdf, freq, nbypass, bypass=0;
					uint32_t *currstats, *currstats2;
					int32_t preds[]=
					{
						N,
						W,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
						N+W-NW,
						W+NE-N,
						N+NE-NNE,
						(WWWW+NNN+NEEE+NEEEE)/4,
					};
#ifdef ENABLE_VAL
					uint64_t val_low, val_range, val_code;
#endif
					(void)NNN	;
					(void)NNNE	;
					(void)NN	;
					(void)NNE	;
					(void)NWW	;
					(void)NW	;
					(void)N		;
					(void)NE	;
					(void)NEE	;
					(void)NEEE	;
					(void)NEEEE	;
					(void)WWWWW	;
					(void)WWWW	;
					(void)WWW	;
					(void)WW	;
					(void)W		;
					(void)eNEE	;
					(void)eNEEE	;
					(void)eW	;
					pred=1<<L1SH>>1;
					for(k=0;k<NPREDS;++k)
						pred+=coeffs[kc][k]*preds[k];
					pred>>=L1SH;
					p0=pred;
					pred+=offset;
					CLAMP2(pred, -128, 127);
					ctx=floor_log2(eW*eW+1);
					if(ctx>NCTX-1)
						ctx=NCTX-1;
					currstats=stats[kc][ctx];
					currstats2=stats[kc][ctx+(ctx<NCTX-1)];
					if(fwd)
					{
						QuantInfoFwd *info;

						error=(char)(yuv[kc]-pred);
						error=error<<1^error>>31;

						info=qtablefwd+error;
						token=info->token;
						nbypass=info->nbypass;
						bypass=info->bypass;
#ifdef _DEBUG
						if(bypass>(1<<nbypass))
							CRASH("");
#endif
					}
					else
					{
						error=0;
						token=0;
					}

					//if(ky==13&&kx==785)//
					//if(ky==10&&kx==iw/2)//
					//if(ky==0&&kx==2)//
					//if(ky==0&&kx==17&&kc==1)//
					//if(ky==13&&kx==2039&&kc==2)//
					//if(ky==31&&kx==2037&&kc==0)//
					//	printf("");

					//token
					if(range<(1<<PROBBITS_USE))
					{
						if(streamptr>=streamend)
						{
							int symidx=3*idx+kc, totalsyms=3*iw*ih;
							
							CRASH("inflation %d/%d  %8.4lf%%\n"
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
					if(fwd)
					{
						cdf =((currstats[token+0]+currstats2[token+0])>>(PROBBITS_STORE+1-PROBBITS_USE))+token+0;
						freq=((currstats[token+1]+currstats2[token+1])>>(PROBBITS_STORE+1-PROBBITS_USE))+token+1;
						//cdf=(currstats[token]>>(PROBBITS_STORE-PROBBITS_USE))+token;
						//freq=(currstats[token+1]>>(PROBBITS_STORE-PROBBITS_USE))+token+1;
					}
					else
					{
						int32_t c=(int32_t)(((code-low)<<PROBBITS_USE|((1<<PROBBITS_USE)-1))/range);
						token=0;
						cdf=0;
						for(;;)
						{
							freq=((currstats[token+1]+currstats2[token+1])>>(PROBBITS_STORE+1-PROBBITS_USE))+token+1;
							if(freq>c)
								break;
#ifdef _DEBUG
							if(token>=NLEVELS)
								CRASH("");
#endif
							++token;
							cdf=freq;
						}
						{
							QuantInfoInv *info=qtableinv+token;
							nbypass=info->nbits;
							error=info->sym;
						}
					}
					freq-=cdf;
#ifdef _DEBUG
					if(
						freq<0
						||(uint32_t)cdf>=(unsigned)(1<<PROBBITS_USE)
						||(uint32_t)(cdf+freq)>(uint32_t)(1<<PROBBITS_USE)
					)
						CRASH("Invalid stats");
					{
						uint64_t low0=low, range0=range;
#endif
#ifdef ENABLE_VAL
					val_low=low;
					val_range=range;
					val_code=code;
#endif
					low+=range*cdf>>PROBBITS_USE;
					range=(range*freq>>PROBBITS_USE)-1;
#ifdef ENABLE_VAL
					if(fwd)
						val_enc(cdf, freq, val_low, val_range, low, range);//
					else
						VAL_DEC(cdf, freq, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
#endif
#ifdef _DEBUG
					if(!fwd&&(code<low||code>low+range||low<low0||low+range>low0+range0||range>=range0))
						CRASH("");
					}
#endif
					{
						int ks;
						uint32_t *currmixin=mixin+NLEVELS-1-token;

						for(ks=1;ks<NLEVELS;++ks)
						{
							uint32_t cell=currstats[ks];
							uint32_t cell2=currstats2[ks];
							cell+=(int32_t)(currmixin[ks]-cell)>>RATE;
							cell2+=(int32_t)(currmixin[ks]-cell2)>>RATE;
#ifdef _DEBUG
							if(cell>(1<<PROBBITS_STORE)-NLEVELS||(ks&&cell<currstats[ks-1]))//
								printf("");
#endif
							currstats[ks]=cell;
							currstats2[ks]=cell2;
						}
					}

					if(nbypass)//bypass
					{
						uint32_t nlevels=1<<nbypass;
						if(range<nlevels)
						{
							if(streamptr>=streamend)
							{
								int symidx=3*idx+kc, totalsyms=3*iw*ih;

								CRASH("inflation %d/%d  %8.4lf%%\n"
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
						if(!fwd)
						{
							bypass=(int32_t)((((code-low+1)<<nbypass)-1)/range);
#ifdef _DEBUG
							if(bypass>=(1<<nbypass))
								CRASH("");
#endif
						}
#ifdef _DEBUG
						{
							uint64_t low0=low, range0=range;
#endif
#ifdef ENABLE_VAL
						val_low=low;
						val_range=range;
						val_code=code;
#endif
						low+=range*bypass>>nbypass;
						range=(range>>nbypass)-1;
#ifdef ENABLE_VAL
						if(fwd)
							val_enc(bypass, 1, val_low, val_range, low, range);//
						else
							VAL_DEC(bypass, 1, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
#endif
#ifdef _DEBUG
						if(!fwd&&(code<low||code>low+range||low<low0||low+range>low0+range0||range>=range0))
							CRASH("");
						}
#endif
					}
					else
						bypass=0;

					if(!fwd)
					{
						error+=bypass<<CONFIG_LSB;
						error=error>>1^-(error&1);
						yuv[kc]=(char)(error+pred);
					}
					rows[0][0]=yuv[kc]-offset;
					error=yuv[kc]-pred;
					rows[0][1]=(2*eW+(error<<1^error>>31)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					{
						int32_t e=rows[0][0];
						e=(e>p0)-(e<p0);
						for(k=0;k<NPREDS;++k)
							coeffs[kc][k]+=e*preds[k];
					}
					offset=kc?yuv[vhelpidx]:yuv[uhelpidx];
					rows[0]+=2;
					rows[1]+=2;
					rows[2]+=2;
					rows[3]+=2;
				}
				if(!fwd)
				{
					imptr[yidx]=yuv[0]+128;
					imptr[uidx]=yuv[1]+128;
					imptr[vidx]=yuv[2]+128;
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
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
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