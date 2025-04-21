#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
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


#ifdef _MSC_VER
	#define LOUD
//	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
#endif


	#define USE_MA	//MA is better here


#define NPREDS 8
#define WP_PREDLIST\
	WP_PRED(260, N)\
	WP_PRED(260, W)\
	WP_PRED(150, 3*(N-NN)+NNN)\
	WP_PRED(150, 3*(W-WW)+WWW)\
	WP_PRED(170, W+NE-N)\
	WP_PRED(170, N+W-NW)\
	WP_PRED(160, N+NE-NNE)\
	WP_PRED(100, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)

#define NCTX 8
#define GRLIMIT 16
#define PROBBITS_STORE 24
#define PROBBITS_USE 9
#define PROBBITS_USE2 16
#define MIXBITS 16
#define SHIFT_LEFT_SIGNED(X, SH) ((SH)<0?(X)>>-(SH):(X)<<(SH))

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
static int crash(const char *file, int line, const char *msg, ...)
{
	printf("%s(%d): ", file, line);
	if(msg)
	{
		va_list args;
		va_start(args, msg);
		vprintf(msg, args);
		va_end(args);
		printf("\n");
	}
	exit(1);
}
#define CRASH(MSG, ...) crash(__FILE__, __LINE__, MSG,##__VA_ARGS__)
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
#ifdef _MSC_VER
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
#ifdef LOUD
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
#endif
#define FRAC_BITS 24
unsigned long long exp2_fix24(int x)
{
	/*
	transcendental fractional powers of two
	x					inv(x)
	2^-0x0.000001 = 0x0.FFFFFF4F...		0x1.000000B1... = 2^0x0.000001
	2^-0x0.000002 = 0x0.FFFFFE9D...		0x1.00000163... = 2^0x0.000002
	2^-0x0.000004 = 0x0.FFFFFD3A...		0x1.000002C6... = 2^0x0.000004
	2^-0x0.000008 = 0x0.FFFFFA74...		0x1.0000058C... = 2^0x0.000008
	2^-0x0.000010 = 0x0.FFFFF4E9...		0x1.00000B17... = 2^0x0.000010
	2^-0x0.000020 = 0x0.FFFFE9D2...		0x1.0000162E... = 2^0x0.000020
	2^-0x0.000040 = 0x0.FFFFD3A3...		0x1.00002C5D... = 2^0x0.000040
	2^-0x0.000080 = 0x0.FFFFA747...		0x1.000058B9... = 2^0x0.000080
	2^-0x0.000100 = 0x0.FFFF4E8E...		0x1.0000B172... = 2^0x0.000100
	2^-0x0.000200 = 0x0.FFFE9D1D...		0x1.000162E5... = 2^0x0.000200
	2^-0x0.000400 = 0x0.FFFD3A3B...		0x1.0002C5CD... = 2^0x0.000400
	2^-0x0.000800 = 0x0.FFFA747F...		0x1.00058BA0... = 2^0x0.000800
	2^-0x0.001000 = 0x0.FFF4E91C...		0x1.000B175F... = 2^0x0.001000
	2^-0x0.002000 = 0x0.FFE9D2B3...		0x1.00162F39... = 2^0x0.002000
	2^-0x0.004000 = 0x0.FFD3A752...		0x1.002C605E... = 2^0x0.004000
	2^-0x0.008000 = 0x0.FFA75652...		0x1.0058C86E... = 2^0x0.008000
	2^-0x0.010000 = 0x0.FF4ECB59...		0x1.00B1AFA6... = 2^0x0.010000
	2^-0x0.020000 = 0x0.FE9E115C...		0x1.0163DAA0... = 2^0x0.020000
	2^-0x0.040000 = 0x0.FD3E0C0D...		0x1.02C9A3E7... = 2^0x0.040000
	2^-0x0.080000 = 0x0.FA83B2DB...		0x1.059B0D32... = 2^0x0.080000
	2^-0x0.100000 = 0x0.F5257D15...		0x1.0B5586D0... = 2^0x0.100000
	2^-0x0.200000 = 0x0.EAC0C6E8...		0x1.172B83C8... = 2^0x0.200000
	2^-0x0.400000 = 0x0.D744FCCB...		0x1.306FE0A3... = 2^0x0.400000
	2^-0x0.800000 = 0x0.B504F334...		0x1.6A09E667... = 2^0x0.800000
	*/
	static const unsigned long long frac_pots[]=
	{
		0x1000000B1,//round(2^(0x000001*2^-24)*2^32)		//extra 8 bits of precision
		0x100000163,
		0x1000002C6,
		0x10000058C,
		0x100000B17,
		0x10000162E,
		0x100002C5D,
		0x1000058B9,
		0x10000B172,
		0x1000162E5,
		0x10002C5CC,
		0x100058BA0,
		0x1000B175F,
		0x100162F39,
		0x1002C605E,
		0x10058C86E,
		0x100B1AFA6,
		0x10163DAA0,
		0x102C9A3E7,
		0x1059B0D31,
		0x10B5586D0,
		0x1172B83C8,
		0x1306FE0A3,
		0x16A09E668,//round(2^(0x800000*2^-24)*2^32)
	};
	int neg=x<0;
	x=abs(x);
	if(x>=0x27000000)
		return neg?0:0xFFFFFFFFFFFFFFFF;
	unsigned long long result=0x1000000;
	for(int k=0;k<FRAC_BITS;++k)//up to 24 muls
	{
		if(x&1)
		{
			result*=frac_pots[k];
			result>>=32;
		}
		x>>=1;
	}
	result=SHIFT_LEFT_SIGNED(result, x);
	if(neg)
		result=0x1000000000000/result;
	return result;
}
static unsigned char* file_load(const char *fn, ptrdiff_t *ret_size)
{
	unsigned char *buf=0;
	ptrdiff_t size=0;
	{
		struct stat info={0};
		int error=stat(fn, &info);
		if(error)
		{
			printf("Cannot stat \"%s\"", fn);
			return 0;
		}
		int regularfile=
#if defined _MSC_VER || defined _WIN32
			(info.st_mode&S_IFMT)==S_IFREG
#else
			S_ISREG(info.st_mode)
#endif
		;
		if(!regularfile)
		{
			printf("Not a file  \"%s\"\n", fn);
			return 0;
		}
		size=info.st_size;
	}
	{
		FILE *fsrc=fopen(fn, "rb");
		if(!fsrc)
		{
			printf("Cannot open \"%s\"\n", fn);
			return 0;
		}
		buf=(unsigned char*)malloc(size+16);
		if(!buf)
		{
			CRASH("Alloc error");
			return 0;
		}
		ptrdiff_t nread=fread(buf, 1, size, fsrc);
		if(nread!=size)
			printf("File size %9td  read %9td bytes\n", size, nread);
		fclose(fsrc);
	}
	*ret_size=size;
	return buf;
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
		CRASH("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
//	OCH_R2=6,
//	OCH_G2,
//	OCH_B2,

	OCH_COUNT=6,
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
static int stats1[3][NCTX][256][256], mixer[3][256][NCTX];
//static int stats1[3][NCTX][256][14][GRLIMIT+8];
//static unsigned stats1[3][512][14][GRLIMIT+8];
static unsigned long long g_low, g_range, g_code;
static unsigned char *g_ptr, *g_ptrstart, *g_ptrend;
static int g_fwd;
static unsigned g_tidx, g_bit;
#ifdef _MSC_VER
static int g_ky, g_iw, g_ih;
#endif
static int *statsptr[NCTX]={0}, *mixptr=0;
static int probtable[1<<PROBBITS_USE];
AWM_INLINE int codebit(int rcon, int sh)
{
	int cells[]=
	{
		statsptr[0][g_tidx],
		statsptr[1][g_tidx],
		statsptr[2][g_tidx],
		statsptr[3][g_tidx],
		statsptr[4][g_tidx],
		statsptr[5][g_tidx],
		statsptr[6][g_tidx],
		statsptr[7][g_tidx],
	};
	long long p0raw=
		+(long long)mixptr[0]*cells[0]
		+(long long)mixptr[1]*cells[1]
		+(long long)mixptr[2]*cells[2]
		+(long long)mixptr[3]*cells[3]
		+(long long)mixptr[4]*cells[4]
		+(long long)mixptr[5]*cells[5]
		+(long long)mixptr[6]*cells[6]
		+(long long)mixptr[7]*cells[7]
	;
//	p0raw/=(long long)mixptr[0]+mixptr[1]+mixptr[2]+mixptr[3]+mixptr[4]+mixptr[5]+mixptr[6]+mixptr[7]+1;
	p0raw>>=MIXBITS+2;
//	CLAMP2(p0raw, 1, (1<<PROBBITS_STORE)-1);
	int p0e=(int)((p0raw+(1<<(PROBBITS_STORE-PROBBITS_USE)>>1))>>(PROBBITS_STORE-PROBBITS_USE));
	p0e+=1<<PROBBITS_USE>>1;
//	p0raw>>=PROBBITS_STORE+MIXBITS+3-PROBBITS_USE;
//	int p0e=(int)p0raw;
	CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
	p0e=probtable[p0e];

	//unsigned p0=g_p0;
	//unsigned p0e=p0>>(PROBBITS_STORE-PROBBITS_USE);
	//CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
	unsigned long long r2;
	if(g_fwd)
	{		
		if(g_range<(1ULL<<PROBBITS_USE2))
		{
			*(unsigned*)g_ptr=(unsigned)(g_low>>32);
			g_ptr+=4;
			g_low=g_low<<32;
			g_range=g_range<<32|0xFFFFFFFF;
			if(g_range>~g_low)
				g_range=~g_low;
		}
		r2=g_range*p0e>>PROBBITS_USE2;
		if(g_bit)
		{
			g_low+=r2;
			g_range-=r2;
		}
		else
			g_range=r2-1;
	}
	else
	{
		if(g_range<(1ULL<<PROBBITS_USE2))
		{
			g_code=g_code<<32|*(unsigned*)g_ptr;
			g_ptr+=4;
			g_low=g_low<<32;
			g_range=g_range<<32|0xFFFFFFFF;
			if(g_range>~g_low)
				g_range=~g_low;
		}
		r2=g_range*p0e>>PROBBITS_USE2;
		unsigned long long mid=g_low+r2;
		g_range-=r2;
		g_bit=g_code>=mid;
		if(g_bit)
			g_low=mid;
		else
			g_range=r2-1;
	}
#ifdef _MSC_VER
	if(g_ptr>g_ptrend)
		CRASH("Pointer OOB  %td vs allowed %td  processed %d/%d  inflation %7.2lf%%",
			g_ptr-g_ptrstart,
			g_ptrend-g_ptrstart,
			g_ky,
			g_ih,
			100.*(g_ptr-g_ptrstart)/(3.*g_iw*g_ky)
		);
#endif
	int proberror=(!g_bit<<PROBBITS_STORE)-(p0e<<(PROBBITS_STORE-PROBBITS_USE2));
//	int truth=(!g_bit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1);
//	int proberror=truth-(int)p0raw;
//	int proberror=(!g_bit<<PROBBITS_USE)-p0e;
	int cells0[]=
	{
		cells[0],
		cells[1],
		cells[2],
		cells[3],
		cells[4],
		cells[5],
		cells[6],
		cells[7],
	};
	//cells[0]=truth+((cells[0]-truth+(1<<3>>1))>>3); statsptr[0][g_tidx]=cells[0];
	//cells[1]=truth+((cells[1]-truth+(1<<3>>1))>>3); statsptr[1][g_tidx]=cells[1];
	//cells[2]=truth+((cells[2]-truth+(1<<3>>1))>>3); statsptr[2][g_tidx]=cells[2];
	//cells[3]=truth+((cells[3]-truth+(1<<3>>1))>>3); statsptr[3][g_tidx]=cells[3];
	//cells[4]=truth+((cells[4]-truth+(1<<3>>1))>>3); statsptr[4][g_tidx]=cells[4];
	//cells[5]=truth+((cells[5]-truth+(1<<3>>1))>>3); statsptr[5][g_tidx]=cells[5];
	//cells[6]=truth+((cells[6]-truth+(1<<3>>1))>>3); statsptr[6][g_tidx]=cells[6];
	//cells[7]=truth+((cells[7]-truth+(1<<3>>1))>>3); statsptr[7][g_tidx]=cells[7];
	//cells[0]+=(truth-cells[0]+(1<<3>>1))>>3; statsptr[0][g_tidx]=cells[0];
	//cells[1]+=(truth-cells[1]+(1<<3>>1))>>3; statsptr[1][g_tidx]=cells[1];
	//cells[2]+=(truth-cells[2]+(1<<3>>1))>>3; statsptr[2][g_tidx]=cells[2];
	//cells[3]+=(truth-cells[3]+(1<<3>>1))>>3; statsptr[3][g_tidx]=cells[3];
	//cells[4]+=(truth-cells[4]+(1<<3>>1))>>3; statsptr[4][g_tidx]=cells[4];
	//cells[5]+=(truth-cells[5]+(1<<3>>1))>>3; statsptr[5][g_tidx]=cells[5];
	//cells[6]+=(truth-cells[6]+(1<<3>>1))>>3; statsptr[6][g_tidx]=cells[6];
	//cells[7]+=(truth-cells[7]+(1<<3>>1))>>3; statsptr[7][g_tidx]=cells[7];
	cells[0]+=(int)((long long)proberror*mixptr[0]>>23); CLAMP2(cells[0], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[0][g_tidx]=cells[0];
	cells[1]+=(int)((long long)proberror*mixptr[1]>>23); CLAMP2(cells[1], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[1][g_tidx]=cells[1];
	cells[2]+=(int)((long long)proberror*mixptr[2]>>23); CLAMP2(cells[2], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[2][g_tidx]=cells[2];
	cells[3]+=(int)((long long)proberror*mixptr[3]>>23); CLAMP2(cells[3], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[3][g_tidx]=cells[3];
	cells[4]+=(int)((long long)proberror*mixptr[4]>>23); CLAMP2(cells[4], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[4][g_tidx]=cells[4];
	cells[5]+=(int)((long long)proberror*mixptr[5]>>23); CLAMP2(cells[5], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[5][g_tidx]=cells[5];
	cells[6]+=(int)((long long)proberror*mixptr[6]>>23); CLAMP2(cells[6], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[6][g_tidx]=cells[6];
	cells[7]+=(int)((long long)proberror*mixptr[7]>>23); CLAMP2(cells[7], -(1<<PROBBITS_STORE>>1), (1<<PROBBITS_STORE>>1)-1); statsptr[7][g_tidx]=cells[7];
	mixptr[0]+=(int)((long long)proberror*cells0[0]>>35); CLAMP2(mixptr[0], 1, 1<<MIXBITS);
	mixptr[1]+=(int)((long long)proberror*cells0[1]>>35); CLAMP2(mixptr[1], 1, 1<<MIXBITS);
	mixptr[2]+=(int)((long long)proberror*cells0[2]>>35); CLAMP2(mixptr[2], 1, 1<<MIXBITS);
	mixptr[3]+=(int)((long long)proberror*cells0[3]>>35); CLAMP2(mixptr[3], 1, 1<<MIXBITS);
	mixptr[4]+=(int)((long long)proberror*cells0[4]>>35); CLAMP2(mixptr[4], 1, 1<<MIXBITS);
	mixptr[5]+=(int)((long long)proberror*cells0[5]>>35); CLAMP2(mixptr[5], 1, 1<<MIXBITS);
	mixptr[6]+=(int)((long long)proberror*cells0[6]>>35); CLAMP2(mixptr[6], 1, 1<<MIXBITS);
	mixptr[7]+=(int)((long long)proberror*cells0[7]>>35); CLAMP2(mixptr[7], 1, 1<<MIXBITS);
	//cells[0]+=(int)((long long)proberror*mixptr[0]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[0], 0, 1<<PROBBITS_STORE); statsptr[0][g_tidx]=cells[0];
	//cells[1]+=(int)((long long)proberror*mixptr[1]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[1], 0, 1<<PROBBITS_STORE); statsptr[1][g_tidx]=cells[1];
	//cells[2]+=(int)((long long)proberror*mixptr[2]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[2], 0, 1<<PROBBITS_STORE); statsptr[2][g_tidx]=cells[2];
	//cells[3]+=(int)((long long)proberror*mixptr[3]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[3], 0, 1<<PROBBITS_STORE); statsptr[3][g_tidx]=cells[3];
	//cells[4]+=(int)((long long)proberror*mixptr[4]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[4], 0, 1<<PROBBITS_STORE); statsptr[4][g_tidx]=cells[4];
	//cells[5]+=(int)((long long)proberror*mixptr[5]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[5], 0, 1<<PROBBITS_STORE); statsptr[5][g_tidx]=cells[5];
	//cells[6]+=(int)((long long)proberror*mixptr[6]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[6], 0, 1<<PROBBITS_STORE); statsptr[6][g_tidx]=cells[6];
	//cells[7]+=(int)((long long)proberror*mixptr[7]>>(MIXBITS+PROBBITS_STORE-16)); CLAMP2(cells[7], 0, 1<<PROBBITS_STORE); statsptr[7][g_tidx]=cells[7];
	//mixptr[0]+=(int)((long long)proberror*cells0[0]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[0], 0, 1<<MIXBITS);
	//mixptr[1]+=(int)((long long)proberror*cells0[1]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[1], 0, 1<<MIXBITS);
	//mixptr[2]+=(int)((long long)proberror*cells0[2]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[2], 0, 1<<MIXBITS);
	//mixptr[3]+=(int)((long long)proberror*cells0[3]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[3], 0, 1<<MIXBITS);
	//mixptr[4]+=(int)((long long)proberror*cells0[4]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[4], 0, 1<<MIXBITS);
	//mixptr[5]+=(int)((long long)proberror*cells0[5]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[5], 0, 1<<MIXBITS);
	//mixptr[6]+=(int)((long long)proberror*cells0[6]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[6], 0, 1<<MIXBITS);
	//mixptr[7]+=(int)((long long)proberror*cells0[7]>>(PROBBITS_STORE*2-7)); CLAMP2(mixptr[7], 0, 1<<MIXBITS);
	//int cells0[]=
	//{
	//	(cells[0]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[1]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[2]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[3]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[4]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[5]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[6]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//	(cells[7]>>(PROBBITS_STORE-PROBBITS_USE))-(1<<PROBBITS_USE>>1),
	//};
	//cells[0]+=proberror*mixptr[0]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[0], 0, 1<<PROBBITS_STORE); statsptr[0][g_tidx]=cells[0];
	//cells[1]+=proberror*mixptr[1]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[1], 0, 1<<PROBBITS_STORE); statsptr[1][g_tidx]=cells[1];
	//cells[2]+=proberror*mixptr[2]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[2], 0, 1<<PROBBITS_STORE); statsptr[2][g_tidx]=cells[2];
	//cells[3]+=proberror*mixptr[3]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[3], 0, 1<<PROBBITS_STORE); statsptr[3][g_tidx]=cells[3];
	//cells[4]+=proberror*mixptr[4]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[4], 0, 1<<PROBBITS_STORE); statsptr[4][g_tidx]=cells[4];
	//cells[5]+=proberror*mixptr[5]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[5], 0, 1<<PROBBITS_STORE); statsptr[5][g_tidx]=cells[5];
	//cells[6]+=proberror*mixptr[6]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[6], 0, 1<<PROBBITS_STORE); statsptr[6][g_tidx]=cells[6];
	//cells[7]+=proberror*mixptr[7]>>(MIXBITS+PROBBITS_USE-16); CLAMP2(cells[7], 0, 1<<PROBBITS_STORE); statsptr[7][g_tidx]=cells[7];
	//mixptr[0]+=proberror*cells0[0]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[0], 0, 1<<MIXBITS);
	//mixptr[1]+=proberror*cells0[1]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[1], 0, 1<<MIXBITS);
	//mixptr[2]+=proberror*cells0[2]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[2], 0, 1<<MIXBITS);
	//mixptr[3]+=proberror*cells0[3]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[3], 0, 1<<MIXBITS);
	//mixptr[4]+=proberror*cells0[4]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[4], 0, 1<<MIXBITS);
	//mixptr[5]+=proberror*cells0[5]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[5], 0, 1<<MIXBITS);
	//mixptr[6]+=proberror*cells0[6]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[6], 0, 1<<MIXBITS);
	//mixptr[7]+=proberror*cells0[7]>>(PROBBITS_USE*2-7); CLAMP2(mixptr[7], 0, 1<<MIXBITS);
	//g_p0=p0+((int)((!g_bit<<PROBBITS_STORE)-p0+rcon)>>sh);
	return p0e;
}
int a03_codec(int argc, char **argv)
{
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef ESTIMATE_SIZE
static double csizes[4]={0};
#endif
	if(argc<3)
	{
		printf("Usage:  \"%s\"  input  output\n", argv[0]);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	ptrdiff_t srcsize=0;
	unsigned char *srcbuf=file_load(srcfn, &srcsize);
	if(!srcbuf)
		return 1;
	unsigned char *srcptr=srcbuf, *srcend=srcbuf+srcsize;
	unsigned char *dstbuf=0, *image=0, *dstptr=0;
	{
		int tag=*(unsigned short*)srcptr;
		g_fwd=tag==('P'|'6'<<8);
		if(!g_fwd&&tag!=('0'|'3'<<8))
		{
			printf("Unsupported file  \"%s\"\n", srcfn);
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	int iw=0, ih=0;
	ptrdiff_t res=0, usize=0;
	if(g_fwd)
	{
		if(*srcptr++!='\n')
		{
			CRASH("Unsupported file");
			return 1;
		}
		while(*srcptr=='#')
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr<=' ')++srcptr;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(*srcptr++!='\n')
		{
			CRASH("Unsupported file");
			return 1;
		}
		while(*srcptr=='#')
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		if(memcmp(srcptr, "255\n", 4))
		{
			CRASH("Unsupported file");
			return 1;
		}
		srcptr+=4;

		image=srcptr;
		res=(ptrdiff_t)iw*ih;
		usize=3*res;

		dstbuf=(unsigned char*)malloc(4*res);
		if(!dstbuf)
		{
			CRASH("Alloc error");
			return 1;
		}
	}
	else
	{
		iw=((int*)srcptr)[0];
		ih=((int*)srcptr)[1];
		srcptr+=8;

		res=(ptrdiff_t)iw*ih;
		usize=3*res;
		dstbuf=(unsigned char*)malloc(usize);
		if(!dstbuf)
		{
			CRASH("Alloc error");
			return 1;
		}
		image=dstbuf;
		dstptr=image;
	}
	if(iw<1||ih<1)
	{
		CRASH("Invalid file");
		return 1;
	}
	int bestrct=0;
#ifdef _MSC_VER
	g_iw=iw;
	g_ih=ih;
#endif
	g_low=0;
	g_range=0xFFFFFFFF;
	if(g_fwd)
	{
		guide_save(srcptr, iw, ih);
		{
			long long counters[OCH_COUNT]={0};
			int prev[OCH_COUNT]={0};
			for(srcptr=image;srcptr<srcend;)
			{
				int
					r=srcptr[0],
					g=srcptr[1],
					b=srcptr[2],
					rg=r-g,
					gb=g-b,
					br=b-r;
				srcptr+=3;
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
			long long minerr=0;
			for(int kt=0;kt<_countof(rct_indices);++kt)
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
			}
#ifdef LOUD
			const unsigned char *rct=rct_indices[bestrct];
			long long currerr=
				+counters[rct[0]]
				+counters[rct[1]]
				+counters[rct[2]]
			;
			printf("%s  WH %d*%d\n", srcfn, iw, ih);
			printf("E/U = %12lld / %12td = %12.6lf\n", currerr, usize, (double)currerr/usize);
#endif
		}
		srcptr=image;
		g_ptr=dstbuf;
		g_ptrstart=dstbuf;
		g_ptrend=dstbuf+4*res;
		g_low+=g_range*bestrct>>4;
		g_range=(g_range>>4)-1;
	}
	else
	{
		g_ptr=srcptr;
		g_ptrstart=srcptr;
		g_ptrend=srcend;
		g_code=*(unsigned*)g_ptr;
		g_ptr+=4;
		g_code=g_code<<4|*(unsigned*)g_ptr;
		g_ptr+=4;

		bestrct=(int)((((g_code-g_low+1)<<4)-1)/g_range);
		g_low+=g_range*bestrct>>4;
		g_range=(g_range>>4)-1;
	}
	int wperrors[3][NPREDS]={0};
	int psize=(iw+16LL)*(int)sizeof(short[4*3*3]);//4 padded rows * 3 channels * {pixels, errors, nbypass}
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	int
		yidx=rct_indices[bestrct][3+0],
		uidx=rct_indices[bestrct][3+1],
		vidx=rct_indices[bestrct][3+2],
		uhelpidx=rct_indices[bestrct][6+0],
		vhelpidx=rct_indices[bestrct][6+1];
	for(int k=0;k<(1<<PROBBITS_USE);++k)
	{
		//squash(x) = sigmoid(x) = 1/(1+exp(-x))
		int input=(k-(1<<PROBBITS_USE>>1))<<(24-PROBBITS_USE)|1<<(24-PROBBITS_USE)>>1;
		probtable[k]=(int)((0x1000000ULL<<PROBBITS_USE2)/(0x1000000+exp2_fix24(-(input*55>>1))));

		//if(g_fwd)
		//	printf("%8d\n", probtable[k]);//
	}
	memset(stats1, 0, sizeof(stats1));
	//FILLMEM((unsigned*)stats1, 1<<PROBBITS_STORE>>1, sizeof(stats1), sizeof(int));
	FILLMEM((unsigned*)mixer, (1<<MIXBITS)/NCTX, sizeof(mixer), sizeof(int));
	//memset(mixer, 0, sizeof(mixer));
	for(int ky=0;ky<ih;++ky)
	{
		short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*3*3,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*3*3,
		};
		char yuv[4]={0};
		int errors[3]={0};
		int preds[NPREDS]={0};
		int pred=0;
		//int ctx=0;
#ifdef _MSC_VER
		g_ky=ky;
#endif
		for(int kx=0;kx<iw;++kx)
		{
#ifndef USE_MA
			int offset=0;
#endif
			if(g_fwd)
			{
				yuv[0]=srcptr[yidx]-128;
				yuv[1]=srcptr[uidx]-128;
				yuv[2]=srcptr[vidx]-128;
				srcptr+=3;
#ifdef USE_MA
				yuv[2]-=yuv[vhelpidx];
				yuv[1]-=yuv[uhelpidx];
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
				int//		rows[-Y][E+X*3*3],
					NNN	=rows[3][0+0*3*3],
					NN	=rows[2][0+0*3*3],
					NNE	=rows[2][0+1*3*3],
					NW	=rows[1][0-1*3*3],
					N	=rows[1][0+0*3*3],
					NE	=rows[1][0+1*3*3],
					NEE	=rows[1][0+2*3*3],
					NEEE	=rows[1][0+3*3*3],
					NEEEE	=rows[1][0+4*3*3],
					WWWW	=rows[0][0-4*3*3],
					WWW	=rows[0][0-3*3*3],
					WW	=rows[0][0-2*3*3],
					W	=rows[0][0-1*3*3],
					eNNN	=rows[3][2+0*3*3],
					eNN	=rows[2][2+0*3*3],
					eNNE	=rows[2][2+1*3*3],
					eNW	=rows[1][2-1*3*3],
					eN	=rows[1][2+0*3*3],
					eNE	=rows[1][2+1*3*3],
					eWWW	=rows[0][2-3*3*3],
					eWW	=rows[0][2-2*3*3],
					eW	=rows[0][2-1*3*3],
					bW	=rows[0][1-1*3*3];
				(void)eNNN;
				(void)eNN;
				(void)eNNE;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eWWW;
				(void)eWW;
				(void)eW;
				static const int wg_weights[]=
				{
#define WP_PRED(WEIGHT, EXPR) WEIGHT,
					WP_PREDLIST
					//WG_PREDLIST0
					//WG_PREDLIST1
					//WG_PREDLIST2
#undef  WP_PRED
				};
				int j=0;
#define WP_PRED(WEIGHT, PRED) preds[j++]=PRED;
				WP_PREDLIST
#undef  WP_PRED
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
				static int divlookup[DIVLUTSIZE]={0};
				if(!*divlookup)
				{
					for(int k=0;k<DIVLUTSIZE;++k)
						divlookup[k]=(1<<NUMBITS)/(k+1);
				}
				unsigned coeff[NPREDS], wsum=0;
				int ipred=0;
				for(int kp=0;kp<NPREDS;++kp)
				{
					//			1	1
					//		2	2	1
					//	1	2	?
					coeff[kp]=//maxerror = 18<<depth
						+wperrors[kc][kp]		//+I	*8
					//	+cN	[kp-1*4*NPREDS]*2	//+eNW	*2
					//	+cN	[kp+0*4*NPREDS]*2	//+eN	*2
					//	+ccurr	[kp-1*4*NPREDS]*2	//+eW	*2
					//	+ccurr	[kp+0*4*NPREDS]		//+eNN	*1
					//	+ccurr	[kp+1*4*NPREDS]		//+eNNE	*1
					//	+cN	[kp+1*4*NPREDS]		//+eNE	*1
					//	+ccurr	[kp-2*4*NPREDS]		//+eWW	*1
					;
				}
				int sh[NPREDS];
				sh[0]=31-(DENBITS-1)-_lzcnt_u32(coeff[0]+1);//invert errros
				sh[1]=31-(DENBITS-1)-_lzcnt_u32(coeff[1]+1);
				sh[2]=31-(DENBITS-1)-_lzcnt_u32(coeff[2]+1);
				sh[3]=31-(DENBITS-1)-_lzcnt_u32(coeff[3]+1);
				sh[4]=31-(DENBITS-1)-_lzcnt_u32(coeff[4]+1);
				sh[5]=31-(DENBITS-1)-_lzcnt_u32(coeff[5]+1);
				sh[6]=31-(DENBITS-1)-_lzcnt_u32(coeff[6]+1);
				sh[7]=31-(DENBITS-1)-_lzcnt_u32(coeff[7]+1);
				if(sh[0]<0)sh[0]=0;
				if(sh[1]<0)sh[1]=0;
				if(sh[2]<0)sh[2]=0;
				if(sh[3]<0)sh[3]=0;
				if(sh[4]<0)sh[4]=0;
				if(sh[5]<0)sh[5]=0;
				if(sh[6]<0)sh[6]=0;
				if(sh[7]<0)sh[7]=0;
				coeff[0]=(wg_weights[0]*divlookup[coeff[0]>>sh[0]]>>sh[0])+(1<<DENBITS>>2)/NPREDS;
				coeff[1]=(wg_weights[1]*divlookup[coeff[1]>>sh[1]]>>sh[1])+(1<<DENBITS>>2)/NPREDS;
				coeff[2]=(wg_weights[2]*divlookup[coeff[2]>>sh[2]]>>sh[2])+(1<<DENBITS>>2)/NPREDS;
				coeff[3]=(wg_weights[3]*divlookup[coeff[3]>>sh[3]]>>sh[3])+(1<<DENBITS>>2)/NPREDS;
				coeff[4]=(wg_weights[4]*divlookup[coeff[4]>>sh[4]]>>sh[4])+(1<<DENBITS>>2)/NPREDS;
				coeff[5]=(wg_weights[5]*divlookup[coeff[5]>>sh[5]]>>sh[5])+(1<<DENBITS>>2)/NPREDS;
				coeff[6]=(wg_weights[6]*divlookup[coeff[6]>>sh[6]]>>sh[6])+(1<<DENBITS>>2)/NPREDS;
				coeff[7]=(wg_weights[7]*divlookup[coeff[7]>>sh[7]]>>sh[7])+(1<<DENBITS>>2)/NPREDS;
				wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];//invert coeff sum
				sh[0]=31-(DENBITS-2)-_lzcnt_u32(wsum);
				wsum=0;
				coeff[0]>>=sh[0];
				coeff[1]>>=sh[0];
				coeff[2]>>=sh[0];
				coeff[3]>>=sh[0];
				coeff[4]>>=sh[0];
				coeff[5]>>=sh[0];
				coeff[6]>>=sh[0];
				coeff[7]>>=sh[0];
				wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];
				ipred=(wsum>>1)
					+coeff[0]*preds[0]
					+coeff[1]*preds[1]
					+coeff[2]*preds[2]
					+coeff[3]*preds[3]
					+coeff[4]*preds[4]
					+coeff[5]*preds[5]
					+coeff[6]*preds[6]
					+coeff[7]*preds[7]
				;
				pred=ipred*divlookup[wsum-1]>>NUMBITS;
				{
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
				}
#ifndef USE_MA
				pred+=offset;
				CLAMP2(pred, -128, 127);
#endif
				int nbypass=31-_lzcnt_u32(bW+1);
				//ctx=nbypass;
				nbypass-=6;
				if(nbypass<0)
					nbypass=0;
				errors[kc]=0;
				//int nzeros=-1, bypass=0;
				g_tidx=0;
				statsptr[0]=stats1[kc][0][(preds[0]	+128)&255];
				statsptr[1]=stats1[kc][1][(preds[1]	+128)&255];
				statsptr[2]=stats1[kc][2][(preds[2]	+128)&255];
				statsptr[3]=stats1[kc][3][(preds[3]	+128)&255];
				statsptr[4]=stats1[kc][4][(preds[4]	+128)&255];
				statsptr[5]=stats1[kc][5][(preds[5]	+128)&255];
				statsptr[6]=stats1[kc][6][(preds[6]	+128)&255];
				statsptr[7]=stats1[kc][7][(pred		+128)&255];
				mixptr=mixer[kc][(preds[7]+128)&255];
				//statsptr[0]=stats1[kc][0][preds[0]+128&255][ctx];
				//statsptr[1]=stats1[kc][1][preds[1]+128&255][ctx];
				//statsptr[2]=stats1[kc][2][preds[2]+128&255][ctx];
				//statsptr[3]=stats1[kc][3][preds[3]+128&255][ctx];
				//statsptr[4]=stats1[kc][4][preds[4]+128&255][ctx];
				//statsptr[5]=stats1[kc][5][preds[5]+128&255][ctx];
				//statsptr[6]=stats1[kc][6][preds[6]+128&255][ctx];
				//statsptr[7]=stats1[kc][7][preds[7]+128&255][ctx];
				//mixptr=mixer[kc][(pred+128)&255];
			//	unsigned *statsptr=stats1[kc][(pred*2+eW/8+255)&511][ctx];
#if 1
				if(g_fwd)
					errors[kc]=(char)yuv[kc];
				else
					errors[kc]=0;
				g_tidx=1;
				for(int kb=7;kb>=0;--kb)
				{
					g_bit=errors[kc]>>kb&1;
					int p0=codebit(1<<7>>1, 7);
					(void)p0;
					errors[kc]|=g_bit<<kb;
					g_tidx=g_tidx*2+g_bit;
#ifdef ESTIMATE_SIZE
					if(g_fwd)
					{
						int norm=1<<PROBBITS_USE;
						csizes[kc]-=log2((double)(g_bit?norm-p0:p0)/norm);
					}
#endif
				}
				if(!g_fwd)
					yuv[kc]=(char)errors[kc];
#endif
#if 0
				if(g_fwd)
				{
					errors[kc]=(char)(yuv[kc]-pred);
					errors[kc]=errors[kc]<<1^errors[kc]>>31;
					nzeros=errors[kc]>>nbypass;
					bypass=errors[kc]&((1<<nbypass)-1);
				}
				else
					errors[kc]=0;
				do
				{
#define GRSHLIST\
	GRSH(6)\
	GRSH(6)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)\
	GRSH(7)
					static const int sh[]=
					{
#define GRSH(X) X,
						GRSHLIST
#undef  GRSH
					};
					static const int rcon[]=
					{
#define GRSH(X) 1<<(X)>>1,
						GRSHLIST
#undef  GRSH
					};
					//unsigned *p0ptr=statsptr+g_tidx;
					if(g_fwd)
						g_bit=nzeros--<=0;
					//g_p0=*p0ptr;
					int p0=codebit(rcon[g_tidx], sh[g_tidx]);
					(void)p0;
					//*p0ptr=g_p0;
#ifdef ESTIMATE_SIZE
					if(g_fwd)
					{
						int norm=1<<PROBBITS_USE;
						//unsigned p0e=g_p0>>(PROBBITS_STORE-PROBBITS_USE);
						//CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
						double bitsize=-log2((double)(g_bit?norm-p0:p0)/norm);
						csizes[kc]+=bitsize;

						//int norm=1<<PROBBITS_STORE;
						//double bitsize=-log2((double)(g_bit?norm-g_p0:g_p0)/norm);
						//csizes[kc]+=bitsize;

						if(isinf(bitsize))
							CRASH("");
					}
#endif
					++g_tidx;
				}while(!g_bit&&g_tidx<GRLIMIT);
				if(g_tidx==GRLIMIT)
				{
					nbypass=8;
					if(g_fwd)
						bypass=errors[kc];
				}
				int tidx=g_tidx;
				for(int kb=nbypass-1, tidx2=1;kb>=0;--kb)
				{
					g_tidx=GRLIMIT+8-nbypass+kb;
					//unsigned *p0ptr=statsptr+GRLIMIT+8-nbypass+kb;
					if(g_fwd)
						g_bit=bypass>>kb&1;
					//g_p0=*p0ptr;
					int p0=codebit(1<<9>>1, 9);
					(void)p0;
					//*p0ptr=g_p0;
#ifdef ESTIMATE_SIZE
					if(g_fwd)
					{
						int norm=1<<PROBBITS_USE;
						//unsigned p0e=g_p0>>(PROBBITS_STORE-PROBBITS_USE);
						//CLAMP2(p0e, 1, (1<<PROBBITS_USE)-1);
						double bitsize=-log2((double)(g_bit?norm-p0:p0)/norm);
						csizes[kc]+=bitsize;

						//int norm=1<<PROBBITS_STORE;
						//double bitsize=-log2((double)(g_bit?norm-g_p0:g_p0)/norm);
						//csizes[kc]+=bitsize;

						if(isinf(bitsize))
							CRASH("");
					}
#endif
					if(!g_fwd)
						bypass|=g_bit<<kb;
					tidx2=tidx2*2+g_bit;
				}
				if(!g_fwd)
				{
					errors[kc]=bypass;
					if(tidx<GRLIMIT)
						errors[kc]|=(tidx-1)<<nbypass;
					errors[kc]=errors[kc]>>1^-(errors[kc]&1);
					yuv[kc]=(char)(errors[kc]+pred);
				}
#endif
#ifdef USE_MA
				int curr=yuv[kc];
#else
				int curr=yuv[kc]-offset;
#endif
				int e=curr-pred;
				e=e<<1^e>>31;
				rows[0][0+0*3*3]=curr;
				rows[0][1+0*3*3]=(2*bW+(e<<6)+rows[1][1+3*4*3])>>2;
				rows[0][2+0*3*3]=e;
				int e2[NPREDS], best=0x7FFFFFFF;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int k=0;k<NPREDS;++k)
				{
					e=abs(curr-preds[k]);
					if(best>e)
						best=e;
					e2[k]=e;
				}
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int k=0;k<NPREDS;++k)
				{
					e=e2[k]-best;
					wperrors[kc][k]+=((e<<8)-wperrors[kc][k]+(1<<3>>1))>>3;
				}
#if 0
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int k=0;k<NPREDS;++k)
				{
					e=yuv[kc]-preds[k];
					e=e<<1^e>>31;
					wperrors[kc][k]+=((e<<8)-wperrors[kc][k]+(1<<3>>1))>>3;
				//	wperrors[kc][k]+=((abs(yuv[kc]-preds[k])<<8)-wperrors[kc][k]+(1<<3>>1))>>3;
				}
#endif
#ifndef USE_MA
				switch(kc)
				{
				case 0:
					offset=yuv[uhelpidx];
					break;
				case 1:
					offset=yuv[vhelpidx];
					break;
				}
#endif
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
			}
			if(!g_fwd)
			{
#ifdef USE_MA
				yuv[1]+=yuv[uhelpidx];
				yuv[2]+=yuv[vhelpidx];
#endif
				dstptr[yidx]=(unsigned char)(yuv[0]+128);
				dstptr[uidx]=(unsigned char)(yuv[1]+128);
				dstptr[vidx]=(unsigned char)(yuv[2]+128);
				dstptr+=3;

				guide_check(dstbuf, kx, ky);//
			}
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(g_fwd)
		{
			//flush
			*(unsigned*)g_ptr=(unsigned)(g_low>>32);
			g_ptr+=4;
			*(unsigned*)g_ptr=(unsigned)g_low;
			g_ptr+=4;

			fwrite("03", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(dstbuf, 1, g_ptr-dstbuf, fdst);
		}
		else
		{
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}

	free(pixels);
	free(srcbuf);
	free(dstbuf);

#if defined _MSC_VER && defined LOUD
	t=time_sec()-t;
	if(g_fwd)
	{
#ifdef ESTIMATE_SIZE
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		csizes[3]=csizes[0]+csizes[1]+csizes[2];
		printf("TYUV%14.2lf%14.2lf%14.2lf%14.2lf  "
			, csizes[3]
			, csizes[0]
			, csizes[1]
			, csizes[2]
		);
#endif
		ptrdiff_t csize=g_ptr-dstbuf+10;
		printf("%12td/%12td  %12.6lf  %8.4lf%%\n", csize, usize, (double)usize/csize, (double)csize*100/usize);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+g_fwd, t, usize/(t*1024*1024));
#endif
	return 0;
}
