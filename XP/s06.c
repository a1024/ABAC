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


#ifdef _MSC_VER
	#define LOUD
#ifdef _DEBUG
//	#define PRINT_RANGE
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
	#define ENABLE_GUIDE2
#endif
#endif

	#define USE_L1
//	#define USE_DIVPRED


#define L1NPREDS 8

#define GRBITS 6


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
AWM_INLINE uint32_t floor_log2(uint32_t n)//n >= 1
{
#ifdef M_LZCNT
	return 31-_lzcnt_u32(n);
#else
	if(n>0x7FFFFFFF)
		return 31;
	{
		__m128i t0=_mm_castpd_si128(_mm_cvtsi32_sd(_mm_setzero_pd(), n));
		t0=_mm_srli_epi64(t0, 52);
		return _mm_cvtsi128_si32(t0)-1023;
	}
#if 0
	float x=(float)n;
	size_t addr=(size_t)&x;
	uint32_t bits=*(uint32_t*)addr;
	return (bits>>23)-127;
#endif
#endif
}
#if 0
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
#endif
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
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)

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
		CRASH("Guide error");
}
#else
#define guide_save(...)
#define guide_check(...)
#endif

#ifdef ENABLE_GUIDE2
static int g2_iw=0, g2_ih=0;
static unsigned char *g2_image=0;
static void guide2_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g2_iw=iw;
	g2_ih=ih;
	g2_image=(unsigned char*)malloc(size);
	if(!g2_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g2_image, image, size);
}
static void guide2_check(unsigned char *image, int idx, int size)
{
	int32_t k;

	for(k=idx;k<idx+size;++k)
	{
		if(image[k]!=g2_image[k])
			CRASH("Guide error  XY %d %d", k%g2_iw, k/g2_iw);
	}
}
#else
#define guide2_save(...)
#define guide2_check(...)
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
		printf("Alloc error\n");
		exit(1);
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
#ifdef ESTIMATE_SIZE
static int32_t hist[3][256]={0};
#endif
AWM_INLINE void predict(unsigned char *image, int32_t iw, int32_t ih, int rct, const int fwd)
{
	int32_t kx, ky, kc;
	int32_t paddedwidth=iw+16;
	unsigned char *imptr=image;
	int
		yidx=rct_indices[rct][3+0],
		uidx=rct_indices[rct][3+1],
		vidx=rct_indices[rct][3+2],
		uhelpidx=rct_indices[rct][6+0],
		vhelpidx=rct_indices[rct][6+1];

	int32_t psize=(int32_t)sizeof(int16_t[4*3*2])*(iw+16);//4 padded rows * 3 channels * {pixels, nbypass}
	int16_t *pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0;ky<ih;++ky)
	{
		char yuv[4]={0};
		int16_t *rows[]=
		{
			pixels+(paddedwidth*((ky-0)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-1)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-2)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-3)&3)+8)*3*2,
		};
		for(kx=0;kx<iw;++kx, imptr+=3)
		{
			int offset=0;

			if(fwd)
			{
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
			}
			for(kc=0;kc<3;++kc)
			{
				int
					NW	=rows[1][-1*3*2],
					N	=rows[1][+0*3*2],
					W	=rows[0][-1*3*2];
				int pred;

				pred=abs(N-NW)<abs(W-NW)?W:N;
#if 0
				int vmax=N, vmin=W;
				if(N<W)
					vmin=N, vmax=W;
				pred=N+W-NW;
				CLAMP2(pred, vmin, vmax);
#endif
				if(kc)
				{
					pred+=offset;
					CLAMP2(pred, -128, 127);
				}
				if(fwd)
				{
					int error=(char)(yuv[kc]-pred);
					error=error<<1^error>>31;
					imptr[kc]=error;
				}
				else
				{
					int error=imptr[kc];
					error=error>>1^-(error&1);
					yuv[kc]=(char)(error+pred);
				}
				rows[0][0]=yuv[kc]-offset;
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
			}
		}
	}
	free(pixels);
}
#define BINLIMIT 16
AWM_INLINE void bin_enc(const unsigned char *syms, int32_t count, unsigned char **pstreamptr, const unsigned char *streamend, unsigned char *imageptr)
{
	unsigned char *streamptr=*pstreamptr;
	if(count<BINLIMIT)
	{
		memcpy(streamptr, syms, count);
		streamptr+=count;
	}
	else
	{
		int32_t k, param;

		for(k=0, param=0;k<count;++k)
			param+=syms[k];
		param/=count;
		param=floor_log2(param+1)-1;
		if(param<0)
			param=0;
		*streamptr++=param;
		{
			uint32_t reg=0;
			int32_t nbits=0;//number of written bits
			int step=1<<param, mask=step-1, p1=param+1;

			for(k=0;k<count;++k)
			{
				int sym=syms[k], nzeros=sym>>param, bypass=1<<param|(sym&mask);

				if(nbits+nzeros<32)
				{
					nbits+=nzeros;
					reg<<=nzeros;
				}
				else
				{
					nzeros-=32-nbits;
					reg<<=32-nbits;
					*(uint32_t*)streamptr=reg;
					streamptr+=4;
					while(nzeros>32)
					{
						*(uint32_t*)streamptr=0;
						streamptr+=4;
						nzeros-=32;
					}
					nbits=nzeros;
					reg=0;
				}
				{
					int t0=32-nbits;
					int sh=p1-t0;
					if(sh<0)
					{
						nbits+=p1;
						reg=reg<<p1|bypass;
					}
					else
					{
						reg=reg<<t0|bypass>>sh;
						bypass&=(1<<sh)-1;
						*(uint32_t*)streamptr=reg;
						streamptr+=4;
						nbits=sh;
						reg=bypass;
					}
				}
				//while(val>=step)
				//{
				//	if(nbits>=32)
				//	{
				//		*(uint32_t*)streamptr=reg;
				//		streamptr+=4;
				//		nbits=0;
				//		reg=0;
				//	}
				//	++nbits;
				//	reg<<=1;
				//	val-=step;
				//}
				//if(nbits>=32)
				//{
				//	*(uint32_t*)streamptr=reg;
				//	streamptr+=4;
				//	nbits=0;
				//	reg=0;
				//}
				//reg=reg<<1|1;
				//++nbits;
				//{
				//	int sh=nbits+param-32;
				//	if(sh>=0)
				//	{
				//		reg=reg<<(32-nbits)|val>>sh;
				//		val&=(1<<sh)-1;
				//		*(uint32_t*)streamptr=reg;
				//		streamptr+=4;
				//		nbits=sh;
				//		reg=val;
				//	}
				//	else
				//		reg=reg<<param|val;
				//}
			}
			if(nbits)
			{
				*(uint32_t*)streamptr=reg<<(32-nbits);
				streamptr+=4;
			}
		}
	}
	*pstreamptr=streamptr;
}
AWM_INLINE void bin_dec(unsigned char *syms, int32_t count, const unsigned char **pstreamptr, const unsigned char *streamend, unsigned char *imageptr, int iw)
{
	const unsigned char *streamptr=*pstreamptr;
	if(count<BINLIMIT)
	{
		memcpy(syms, streamptr, count);
		streamptr+=count;
	}
	else
	{
		uint32_t reg;
		int32_t nbits;//number of unread bits
		int32_t k;
		int param=*streamptr++, mask=(1<<param)-1, p1=param+1;
		if(param>8)
		{
			CRASH("Invalid GR param %d", param);
			return;
		}

		reg=*(uint32_t*)streamptr;
		streamptr+=4;
		nbits=32;
		for(k=0;k<count;++k)
		{
			int sym=0, rem;

			while(!reg)
			{
				sym+=nbits;
				nbits=32;
				reg=*(uint32_t*)streamptr;
				streamptr+=4;
			}
			rem=31-floor_log2(reg);
			sym+=rem;
			nbits-=rem;
			reg<<=rem;
			if(nbits<p1)
			{
				int32_t sh=p1-nbits;
				rem=reg>>(32-nbits);
				nbits=32-sh;
				reg=*(uint32_t*)streamptr;
				streamptr+=4;
				rem=rem<<sh|reg>>nbits;
				reg<<=sh;
				sym=sym<<param|(rem&mask);
			}
			else
			{
				sym=sym<<param|(reg>>(32-p1)&mask);
				nbits-=p1;
				reg<<=p1;
			}
			syms[k]=sym;
#ifdef ENABLE_GUIDE2
			if(syms[k]!=g2_image[syms+k-imageptr])
			{
				int32_t idx=(int32_t)(syms+k-imageptr);
				int kc=idx%3;
				idx/=3;
				CRASH("Guide error  XYC %d %d %d", idx%iw, idx/iw, kc);
			}
#else
			(void)imageptr;
			(void)iw;
#endif
		}
	}
	*pstreamptr=streamptr;
}
int s06_codec(const char *command, const char *srcfn, const char *dstfn)
{
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *streamptr=0, *streamstart=0, *streamend=0;
	unsigned char *imptr=0;
	int bestrct=0;
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef ESTIMATE_SIZE
	int32_t csizes_z[3]={0};
#endif

	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		CRASH("");
		return 1;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'6'<<8))
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
#ifdef LOUD
		printf("\"%s\"  WH %d*%d\n", srcfn, iw, ih);
#endif
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

#ifdef ESTIMATE_SIZE
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
#ifdef ESTIMATE_SIZE
				printf("RCT%02d %16lld%s\n", kt, currerr, kt==bestrct?" <-":"");
#endif
			}
		}
		memset(dstbuf, 0, usize);
		streamstart=dstbuf;
		streamend=dstbuf+usize;
		streamptr=dstbuf;
		(void)streamstart;

		*streamptr++=bestrct;
	}
	else
	{
		imptr=image=dstbuf;
		streamstart=srcptr;
		streamend=srcbuf+srcsize;
		streamptr=srcptr;

		bestrct=*streamptr++;

		csize=srcsize;
	}
	if(fwd)
	{
		unsigned char *imend=image+usize, *dataptr=image;
#ifdef LOUD
		int ectr=0;
#endif

		predict(image, iw, ih, bestrct, 1);
		guide2_save(image, iw, ih);
#ifdef ESTIMATE_SIZE
		{
			int32_t k;

			for(k=0;k<usize;k+=3)
			{
				++hist[0][image[k+0]];
				++hist[1][image[k+1]];
				++hist[2][image[k+2]];
			}
		}
#endif
		imptr=image;
		while(imptr<imend)
		{
#define MINRUN 8
			unsigned char *ptr2;

			for(ptr2=imptr;ptr2<imend&&!*ptr2;++ptr2);
			if(ptr2-imptr<MINRUN)//skip small RLE as data
			{
				for(imptr=ptr2;imptr<imend&&*imptr;++imptr);
				continue;
			}
			if(dataptr<imptr)
			{
				int32_t size=(int32_t)(imptr-dataptr);
#ifdef _DEBUG
				if(streamptr+8+size>streamend)
					CRASH("Inflation %8.4lf%%", 100.*usize/(ptr2-imptr));
#endif
				{
					int32_t emit=size<<1|1;
					do
					{
						*streamptr++=(emit>127)<<7|(emit&127);
						emit>>=7;
					}while(emit);
				}
#ifdef LOUD
				++ectr;
#endif
#ifdef PRINT_RANGE
				{
					static int32_t count=0;
					int32_t k, vmin, vmax;

					if(count<1000)
					{
						vmax=vmin=128;
						for(k=0;k<size;++k)
						{
							int val=dataptr[k];
							if(vmax<val)vmax=val;
							if(vmin>val)vmin=val;
						}
						printf("%8d", size);
						printf("%*s", vmin, "");
						for(k=vmin;k<vmax;++k)
							printf("*");
						printf("\n");
						++count;
					}
				}
#endif
				bin_enc(dataptr, size, &streamptr, streamend, image);
			}
			{
				int32_t emit=(int32_t)(ptr2-imptr-MINRUN)<<1;
				do
				{
					*streamptr++=(emit>127)<<7|(emit&127);
					emit>>=7;
				}while(emit);
#ifdef LOUD
				++ectr;
#endif
			}
			dataptr=imptr=ptr2;
		}
		if(dataptr<imptr)
		{
			int32_t size=(int32_t)(imptr-dataptr);
			{
				int32_t emit=size<<1|1;
				do
				{
					*streamptr++=(emit>127)<<7|(emit&127);
					emit>>=7;
				}while(emit);
			}
			bin_enc(dataptr, size, &streamptr, streamend, image);
			//memcpy(streamptr, dataptr, size);
			//streamptr+=size;
		}
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing\n", dstfn);
				return 1;
			}
			csize+=fwrite("06", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 4, fdst);
			csize+=fwrite(&ih, 1, 4, fdst);
			csize+=fwrite(dstbuf, 1, streamptr-dstbuf, fdst);
			fclose(fdst);
		}
		dstsize=csize;
	}
	else
	{
#ifdef LOUD
		int ectr=0;
#endif

		imptr=image;
		while(streamptr<streamend)
		{
			int32_t code;

			code=*streamptr++;//1
			if(code>127)
			{
				int32_t c2=*streamptr++;
				code=(c2&127)<<7|(code&127);//2
				if(c2>127)
				{
					c2=*streamptr++;
					code|=(c2&127)<<14;//3
					if(c2>127)//runs up to 0x1FFFFFFF = 512 MB
						code|=*streamptr++<<21;//4
				}
			}
			if(code&1)//data
			{
				code>>=1;
				if(!code||imptr+code>image+usize)
				{
					CRASH("Data block  size %d = remaining %d", code, image+usize-imptr);
					return 1;
				}
				bin_dec(imptr, code, (const unsigned char**)&streamptr, streamend, image, iw);
				//memcpy(imptr, streamptr, code);
				//streamptr+=code;
				guide2_check(image, (int32_t)(imptr-image), code);
			}
			else//zero
			{
				code=(code>>1)+MINRUN;
				if(!code||imptr+code>image+usize)
				{
					CRASH("Zero block  size %d  remaining %d", code, image+usize-imptr);
					return 1;
				}
				memset(imptr, 0, code);
				guide2_check(image, (int32_t)(imptr-image), code);
			}
			imptr+=code;
#ifdef LOUD
			++ectr;
#endif
		}
		if(streamptr!=srcend)
			CRASH("Bad archive  ptr %08X != end %08X", (int32_t)(size_t)streamptr, (int32_t)(size_t)srcend);
		predict(image, iw, ih, bestrct, 0);
#ifdef ENABLE_GUIDE
		{
			int32_t k;
			for(k=0;k<usize;k+=3)
			{
				if(memcmp(image+k, g_image+k, 3))
					CRASH("Guide error XY %d %d", k%iw, k/iw);
			}
		}
#endif
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing\n", dstfn);
				free(srcbuf);
				free(dstbuf);
				return 1;
			}
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(dstbuf, 1, usize, fdst);
			fclose(fdst);
		}
	}
	free(srcbuf);
	free(dstbuf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef ESTIMATE_SIZE
		{
			double csizes[3]={0};
			int kc, ks;
			for(kc=0;kc<3;++kc)
			{
				double norm;
				int32_t sum=0;
				for(ks=0;ks<256;++ks)
					sum+=hist[kc][ks];
				norm=1./sum;
				for(ks=0;ks<256;++ks)
				{
					int32_t freq=hist[kc][ks];
					if(freq)
						csizes[kc]-=freq*log(freq*norm);
				}
				csizes[kc]*=1./(M_LN2*8);//convert ln->log2
			}
			printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf\n"
				, csizes[0]+csizes[1]+csizes[2]
				, csizes[0]
				, csizes[1]
				, csizes[2]
			);
		}
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
