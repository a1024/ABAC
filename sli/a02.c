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

#define SRCBUFSIZE (1<<20)
#define DSTBUFSIZE (1<<20)
//#define SRCBUFSIZE (1<<26)
//#define DSTBUFSIZE (1<<26)
static unsigned char srcbuf[SRCBUFSIZE], dstbuf[DSTBUFSIZE];
static int g_fwd=0;
static ptrdiff_t process(const char *srcfn, const char *dstfn, ptrdiff_t srcbufsize, ptrdiff_t dstbufsize)
{
	FILE *fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"\n", srcfn);
		return 0;
	}
	ptrdiff_t nread=fread(srcbuf, 1, 128, fsrc);
	unsigned char *srcptr=srcbuf, *srcend=srcbuf+nread;
	{
		int tag=*(unsigned short*)srcptr;
		g_fwd=tag==('P'|'6'<<8);
		if(!g_fwd&&tag!=('0'|'1'<<8))
		{
			printf("Unsupported file  \"%s\"\n", srcfn);
			return 0;
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
			return 0;
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
			return 0;
		}
		while(*srcptr=='#')
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		if(memcmp(srcptr, "255\n", 4))
		{
			CRASH("Unsupported file");
			return 0;
		}
		srcptr+=4;

		res=(ptrdiff_t)iw*ih;
		usize=3*res;
	}
	else
		CRASH("This is not a codec\n");
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		printf("Cannot open \"%s\" for writing\n", dstfn);
		return 0;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih/3);
	unsigned char *dstptr=dstbuf, *dstend=dstbuf+dstbufsize;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				if(srcptr>=srcend)
				{
					srcptr=srcbuf;
					nread=fread(srcbuf, 1, srcbufsize, fsrc);
					if(!nread)
						goto finish;
					srcend=srcbuf+nread;
				}
				if(!(ky%3))
				{
					if(dstptr>=dstend)
					{
						fwrite(dstbuf, 1, dstbufsize, fdst);
						dstptr=dstbuf;
					}
					int val=*srcptr;
					*dstptr++=val;
				}
				++srcptr;
			}
		}
	}
finish:
	if(dstptr>dstbuf)
		fwrite(dstbuf, 1, dstptr-dstbuf, fdst);
	fclose(fsrc);
	fclose(fdst);
	return usize;
}
int a02_codec(int argc, char **argv)
{
	if(argc<3)
	{
		printf("Usage:  \"%s\"  input  output\n", argv[0]);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	double tbest=0;
	int sbest=0, dbest=0;
	ptrdiff_t usize2=0;
	const int rstart=7, rend=20+1;
	printf("s\\d  ");
	for(int kd=rstart;kd<rend;++kd)
		printf("%7d", kd);
	printf("\n\n");
	for(int ks=rstart, it=0;ks<rend;++ks)
	{
		printf("%3d  ", ks);
		for(int kd=rstart;kd<rend;++kd, ++it)
		{
			double t=time_sec();
			ptrdiff_t usize=process(srcfn, dstfn, 1ULL<<ks, 1ULL<<kd);
			t=time_sec()-t;
			printf("%7d", (int)(usize/(t*1024*1024)));
		//	printf("%3d %3d  %lf sec  %lf MB/s\n", ks, kd, t, usize/(t*1024*1024));
			if(!it||tbest>t)
			{
				tbest=t;
				sbest=ks;
				dbest=kd;
				usize2=usize;
			}
		}
		printf("\n");
	}
	printf("Best: (1<<ks, 1<<kd)\n");
	printf("%3d %3d  %lf sec  %lf MB/s\n", sbest, dbest, tbest, usize2/(tbest*1024*1024));
	return 0;
}
