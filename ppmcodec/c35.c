#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#undef NDEBUG//assert.h
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<assert.h>
#include<sys/stat.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#else
#include<time.h>
#endif


#if defined _MSC_VER && !defined RELEASE
	#define ENABLE_GUIDE
#endif

	#define USE_CARRYLESS


#ifdef _MSC_VER
#define INLINE __forceinline static
#else
#define INLINE __attribute__((always_inline)) inline static
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		printf("\n");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
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
#ifdef ENABLE_GUIDE
static unsigned char *g_buffer;
static ptrdiff_t g_usize;
static void guide_save(unsigned char *buffer, int usize)
{
	assert(g_buffer=(unsigned char*)malloc(usize));
	g_usize=usize;
	memcpy(g_buffer, buffer, usize);
}
static void guide_save_file(FILE *fsrc)
{
	ptrdiff_t idx=ftell(fsrc);
	fseek(fsrc, 0, SEEK_END);
	ptrdiff_t end=ftell(fsrc);
	g_usize=end-idx;
	assert(g_buffer=(unsigned char*)malloc(g_usize));
	fseek(fsrc, (long)idx, SEEK_SET);
	fread(g_buffer, 1, g_usize, fsrc);
	fseek(fsrc, (long)idx, SEEK_SET);
}
static void guide_check(int sym, ptrdiff_t idx)
{
	if(g_buffer[idx]!=sym)
	{
		CRASH("Guide error  IDX %td", idx);
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_save_file(...)
#define guide_check(...)
#endif
static ptrdiff_t filesize(const char *fn)
{
	struct stat info={0};
	assert(!stat(fn, &info));
	return info.st_size;
}
enum
{
	PROBBITS=15,
	NSTATES=1<<15,
};

typedef struct _FSMCell
{
	uint16_t next[2];
	uint32_t p0;
} FSMCell;
static FSMCell fsm[NSTATES];
static uint16_t stats[0x10000][256];
static void codec(const char *srcfn, const char *dstfn, const char *fsmfn, const int fwd)
{
	ptrdiff_t usize;
	FILE *fsrc, *fdst, *fu, *fc;
	double t;

	t=time_sec();
	{
		FILE *ffsm;
		int s0=0, s1=0, p0=0, k;

		assert(ffsm=fopen(fsmfn, "rb"));
		for(k=0;k<NSTATES;++k)
		{
			FSMCell *cell=fsm+k;
			int nread=fscanf(ffsm, " %d,%d,%d", &s0, &s1, &p0);
			if(nread!=3)
				break;
			//assert(nread==3);
			CLAMP2(s0, 0, NSTATES-1);
			CLAMP2(s1, 0, NSTATES-1);
			CLAMP2(p0, 1, (1<<PROBBITS)-1);
			cell->next[0]=s0;
			cell->next[1]=s1;
			cell->p0=p0;
		}
		for(;k<NSTATES;++k)
		{
			FSMCell *cell=fsm+k;
			cell->next[0]=0;
			cell->next[1]=0;
			cell->p0=1<<PROBBITS>>1;
		}
		fclose(ffsm);
	}
	memset(stats, 0, sizeof(stats));

	usize=0;
	if(fwd)
		usize=filesize(srcfn);
	assert(fsrc=fopen(srcfn, "rb"));
	assert(fdst=fopen(dstfn, "wb"));
	if(fwd)
		fu=fsrc, fc=fdst;
	else
		fu=fdst, fc=fsrc;
	if(fwd)
	{
		fwrite(&usize, 1, 4, fc);
		guide_save_file(fu);
	}
	else
		fread(&usize, 1, 4, fc);
	{
		ptrdiff_t k;
		int sym=0;
		uint32_t ctx=0;
#ifdef USE_CARRYLESS
		uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
#else
		uint32_t range=-1, code=0, ffnum=0, cache=-1;
		uint64_t lowc=0;
#endif
		if(!fwd)
		{
#ifdef USE_CARRYLESS
			code=0;
			code|=(uint64_t)fgetc(fc);
			code|=(uint64_t)fgetc(fc)<<8;
			code|=(uint64_t)fgetc(fc)<<16;
			code|=(uint64_t)fgetc(fc)<<24;
			code<<=32;
			code|=(uint64_t)fgetc(fc);
			code|=(uint64_t)fgetc(fc)<<8;
			code|=(uint64_t)fgetc(fc)<<16;
			code|=(uint64_t)fgetc(fc)<<24;
#else
			code=code<<8|fgetc(fc);
			code=code<<8|fgetc(fc);
			code=code<<8|fgetc(fc);
			code=code<<8|fgetc(fc);
#endif
		}
		for(k=0;k<usize;++k)
		{
			uint16_t *statsptr=stats[ctx];
			if(fwd)
				sym=fgetc(fu);
			else
				sym=0;

			for(int kb=7, tidx=1;kb>=0;--kb)
			{
				FSMCell *cell=fsm+statsptr[tidx];
				int bit=sym>>kb&1;
				
#ifdef USE_CARRYLESS
				if(range<0xFFFF)
				{
					if(fwd)
					{
						uint32_t c=(uint32_t)(low>>32);
						putc((uint8_t)c, fc); c>>=8;
						putc((uint8_t)c, fc); c>>=8;
						putc((uint8_t)c, fc); c>>=8;
						putc((uint8_t)c, fc);
					}
					else
					{
						code<<=32;
						code|=(uint64_t)fgetc(fc);
						code|=(uint64_t)fgetc(fc)<<8;
						code|=(uint64_t)fgetc(fc)<<16;
						code|=(uint64_t)fgetc(fc)<<24;
					}
					low<<=32;
					range=range<<32|0xFFFFFFFF;
				}
				{
					uint64_t r2=range*cell->p0>>PROBBITS;
					uint64_t mid=low+r2;
					range-=r2;
					if(!fwd)
						bit=code>=mid;
					if(bit)
						low=mid;
					else
						range=r2-1;
				}
#else
				while(range<=0xFFFFFF)
				{
					range<<=8;
					if(!fwd)
						code=code<<8|fgetc(fc);
					else
					{
						uint32_t carry=(uint32_t)(lowc>>32), low=(uint32_t)lowc;
						if(low<0xFF000000||carry)
						{
							if(cache!=-1)
								putc(cache+carry, fc);
							while(ffnum)
							{
								putc(carry-1, fc);
								--ffnum;
							}
							cache=low>>24;
						}
						else
							++ffnum;
						low<<=8;
						lowc=low;
					}
				}
				{
					uint32_t r2=(uint64_t)range*cell->p0>>PROBBITS;
					if(!fwd)
						bit=code>=r2;
					if(bit)
					{
						range-=r2;
						if(fwd)
							lowc+=r2;
						else
							code-=r2;
					}
					else
						range=r2;
				}
#endif
				sym|=bit<<kb;

				//update
				statsptr[tidx]=cell->next[bit];

				tidx=2*tidx+bit;
			}
			if(!fwd)
			{
				guide_check(sym, k);
				putc(sym, fu);
			}
			ctx=(uint8_t)ctx<<8|sym;
		}
		if(fwd)
		{
#ifdef USE_CARRYLESS
			uint32_t c=(uint32_t)(low>>32);
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc);
			c=(uint32_t)low;
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc); c>>=8;
			putc((uint8_t)c, fc);
#else
			for(k=0;k<5;++k)
			{
				uint32_t carry=(uint32_t)(lowc>>32), low=(uint32_t)lowc;
				if(low<0xFF000000||carry)
				{
					if(cache!=-1)
						putc(cache+carry, fc);
					while(ffnum)
					{
						putc(carry-1, fc);
						--ffnum;
					}
					cache=low>>24;
				}
				else
					++ffnum;
				low<<=8;
				lowc=low;
			}
#endif
		}
	}
	fclose(fc);
	fclose(fu);

	t=time_sec()-t;
	if(fwd)
	{
		ptrdiff_t csize=filesize(dstfn);
		printf("%12lld->%12lld bytes  %12.6lf:1  %12.6lf%%\n", usize, csize, (double)usize/csize, 100.*csize/usize);
	}
	printf("%12.6lf sec  %12.6lf MB/s\n", t, usize/(t*1024*1024));
}
int c35_codec(int argc, char **argv)
{
#ifdef __GNUC__
	const char *command, *srcfn, *dstfn, *fsmfn;
	if(argc!=5)
	{
		printf("Usage  \"%s\"  c|d  input  output  FSM.txt\n", argv[0]);
		return 1;
	}
	command=argv[1];
	srcfn=argv[2];
	dstfn=argv[3];
	fsmfn=argv[4];
	if(command[0]=='c'&&!command[1])
		codec(srcfn, dstfn, fsmfn, 1);
	else if(command[0]=='d'&&!command[1])
		codec(srcfn, dstfn, fsmfn, 0);
#elif defined _MSC_VER
	codec("E:/C/Grapher 2 2010/Grapher 2/g2.cpp",	"E:/C/a/FSM_v1/zzz.m01", "E:/C/a/FSM_v1/FSM0.txt", 1);
	codec("E:/C/a/FSM_v1/zzz.m01",			"E:/C/a/FSM_v1/zzz.zzz", "E:/C/a/FSM_v1/FSM0.txt", 0);

//	codec("E:/C/Grapher 2 2013/Grapher 2 2013/g2.cpp",	"E:/C/a/FSM_v1/zzz.m01", "E:/C/a/FSM_v1/FSM0.txt", 1);
//	codec("E:/C/a/FSM_v1/zzz.m01",				"E:/C/a/FSM_v1/zzz.zzz", "E:/C/a/FSM_v1/FSM0.txt", 0);

//	codec("D:/ML/LPCB001-decorr.ppm",	"E:/C/a/FSM_v1/zzz.m01", "E:/C/a/FSM_v1/FSM0.txt", 1);
//	codec("E:/C/a/FSM_v1/zzz.m01",		"E:/C/a/FSM_v1/zzz.zzz", "E:/C/a/FSM_v1/FSM0.txt", 0);

//	codec("E:/C/a/FSM_v1/book1", "E:/C/a/FSM_v1/zzz.m01", "E:/C/a/FSM_v1/FSM0.txt", 1);
//	codec("E:/C/a/FSM_v1/zzz.m01", "E:/C/a/FSM_v1/zzz.txt", "E:/C/a/FSM_v1/FSM0.txt", 0);

	exit(0);
#endif
	return 0;
}
