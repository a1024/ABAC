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
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
	#define PRINT_ESTIMS
//	#define PRINTBITS
#endif
#endif

	#define USE_L1


#define NLEARNERS 8

#define NPREDS 8


#define GRBITS 6
#define GRLIMIT 24
#define PROBBITS_STORE 16
#define PROBBITS_USE 9


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
int floor_log2(int n)
{
	__m128i t0=_mm_castpd_si128(_mm_cvtsi32_sd(_mm_setzero_pd(), n));
	t0=_mm_srli_epi64(t0, 52);
	return _mm_cvtsi128_si32(t0)-1023;
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
#ifdef ESTIMATE_SIZE
static int32_t hist[3][256]={0};
#endif
static uint16_t stats[3][256][256][NLEARNERS];
static uint32_t Shannon[1<<PROBBITS_USE];
static double log2table[1<<PROBBITS_STORE];
int s04_codec(const char *command, const char *srcfn, const char *dstfn)
{
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
	uint32_t low=0, range=0xFFFFFFFF, code=0;
	int64_t totalsizes[3][NLEARNERS]={0};
#ifdef PRINT_ESTIMS
	double allsizes[3][NLEARNERS]={0};
	double esizes[3]={0};
#endif

	{
		int32_t kp;

		for(kp=0;kp<(1<<PROBBITS_USE);++kp)
		{
			int prob=kp;
			CLAMP2(prob, 1, (1<<PROBBITS_USE)-1);
			Shannon[kp]=(int32_t)(-log((double)prob*(1./(1<<PROBBITS_USE)))*((1<<16)/M_LN2));
		}
		for(kp=0;kp<(1<<PROBBITS_STORE);++kp)
		{
			int prob=kp;
			CLAMP2(prob, 1, (1<<PROBBITS_STORE)-1);
			log2table[kp]=-log((double)prob*(1./(1<<PROBBITS_STORE)))*(1./M_LN2);
		}
	}
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

		code=code<<16|*(uint16_t*)streamptr;
		streamptr+=2;
		code=code<<16|*(uint16_t*)streamptr;
		streamptr+=2;

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
		int16_t *pixels=0;
		int32_t ssize=0;
		int32_t *sbuf=0;
		int32_t paddedwidth=iw+16;
#ifdef USE_L1
		int32_t coeffs[NPREDS+1]={0};
#endif

		psize=(int32_t)sizeof(int16_t[4*3*2])*paddedwidth;//4 padded rows * 3 channels * {pixels, nbypass}
		pixels=(int16_t*)malloc(psize);
		ssize=(int32_t)sizeof(int32_t[4*3*NLEARNERS])*paddedwidth;//4 padded rows * 3 channels * NLEARNERS rates
		sbuf=(int32_t*)malloc(ssize);
		if(!pixels)
		{
			crash("Alloc error\n");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		memset(sbuf, 0, ssize);
		FILLMEM((uint16_t*)stats, 0x8000, sizeof(stats), sizeof(int16_t));
		for(ky=0, idx=0;ky<ih;++ky)
		{
			char yuv[4]={0};
			int16_t *rows[]=
			{
				pixels+(paddedwidth*((ky-0)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-1)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-2)&3)+8)*3*2,
				pixels+(paddedwidth*((ky-3)&3)+8)*3*2,
			};
			int32_t *srows[]=
			{
				sbuf+(paddedwidth*((ky-0)&3)+8)*3*NLEARNERS,
				sbuf+(paddedwidth*((ky-1)&3)+8)*3*NLEARNERS,
				sbuf+(paddedwidth*((ky-2)&3)+8)*3*NLEARNERS,
				sbuf+(paddedwidth*((ky-3)&3)+8)*3*NLEARNERS,
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
						NNE	=rows[2][0+1*3*2],
						NN	=rows[2][0+0*3*2],
						NW	=rows[1][0-1*3*2],
						N	=rows[1][0+0*3*2],
						NE	=rows[1][0+1*3*2],
						NEE	=rows[1][0+2*3*2],
						NEEE	=rows[1][0+3*3*2],
						NEEEE	=rows[1][0+4*3*2],
						WWWW	=rows[0][0-4*3*2],
						WWW	=rows[0][0-3*3*2],
						WW	=rows[0][0-2*3*2],
						W	=rows[0][0-1*3*2];
					//	eNEEE	=rows[1][1-3*3*2],
					//	eW	=rows[0][1-1*3*2];
					int32_t
						*sNW	=srows[1]-1*3*NLEARNERS,
						*sN	=srows[1]+0*3*NLEARNERS,
						*sNE	=srows[1]+1*3*NLEARNERS,
						*sNEE	=srows[1]+2*3*NLEARNERS,
						*sNEEE	=srows[1]+3*3*NLEARNERS,
						*sW	=srows[0]-1*3*NLEARNERS,
						*scurr	=srows[0];
					int32_t model;
					int32_t pred, vmax, vmin, pred0;
					int32_t error;
					int32_t kb, tidx, bit;
					int32_t currsizes[NLEARNERS]={0};
#ifdef USE_L1
					int32_t preds[]=
					{
						N,
						W,
						3*(W-WW)+WWW,
						3*(N-NN)+NNN,
						N+W-NW,
						W+NE-N,
						N+NE-NNE,
						(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
					};

					pred=coeffs[NPREDS];
					{
						int j=0;
						for(;j<NPREDS;++j)
							pred+=coeffs[j]*preds[j];
					}
#define L1SH 19
					pred+=1<<L1SH>>1;
					pred>>=L1SH;
					pred0=pred;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
#else
					pred=N+W-NW;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(pred, vmin, vmax);
#endif
					pred+=offset;
					CLAMP2(pred, -128, 127);
					{
						//int32_t sh=floor_log2(idx+1);
						int64_t bestsize=totalsizes[kc][0];
						//int32_t bestsize=sNW[0]+sN[0]+sNE[0]+sNE[0]+sW[0];
						int32_t kp;

						model=0;
						for(kp=1;kp<NLEARNERS;++kp)
						{
							int64_t size=totalsizes[kc][kp];
							//int32_t size=sNW[kp]+sN[kp]+sNE[kp]+sNE[kp]+sW[kp];
							if(bestsize>size)
							{
								bestsize=size;
								model=kp;
							}
						}
					}
					if(fwd)
					{
						error=yuv[kc]-pred;
#ifdef ESTIMATE_SIZE
						++hist[kc][(unsigned char)(error+128)];
#endif
					}
					else
						error=0;
					for(kb=7, tidx=1;kb>=0;--kb)
					{
						uint16_t *pp0=stats[kc][(pred+128)&255][tidx];
						int32_t p0;

#if 0
						p0=(
							+((uint32_t)(totalsizes[kc][0]>>32)+1)*(pp0[0]>>4)
							+((uint32_t)(totalsizes[kc][1]>>32)+1)*(pp0[1]>>4)
							+((uint32_t)(totalsizes[kc][2]>>32)+1)*(pp0[2]>>4)
							+((uint32_t)(totalsizes[kc][3]>>32)+1)*(pp0[3]>>4)
							+((uint32_t)(totalsizes[kc][4]>>32)+1)*(pp0[4]>>4)
							+((uint32_t)(totalsizes[kc][5]>>32)+1)*(pp0[5]>>4)
							+((uint32_t)(totalsizes[kc][6]>>32)+1)*(pp0[6]>>4)
							+((uint32_t)(totalsizes[kc][7]>>32)+1)*(pp0[7]>>4)
						)/(
							+(uint32_t)(totalsizes[kc][0]>>32)
							+(uint32_t)(totalsizes[kc][1]>>32)
							+(uint32_t)(totalsizes[kc][2]>>32)
							+(uint32_t)(totalsizes[kc][3]>>32)
							+(uint32_t)(totalsizes[kc][4]>>32)
							+(uint32_t)(totalsizes[kc][5]>>32)
							+(uint32_t)(totalsizes[kc][6]>>32)
							+(uint32_t)(totalsizes[kc][7]>>32)
							+8
						);
#endif
						p0=pp0[model]>>(PROBBITS_STORE-PROBBITS_USE);
						CLAMP2(p0, 1, (1<<PROBBITS_USE)-1);
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
								*(uint16_t*)streamptr=(uint16_t)(low>>16);
							else
								code=code<<16|*(uint16_t*)streamptr;
							streamptr+=2;
							low<<=16;
							range=range<<16|0xFFFF;
							if(range>~low)
								range=~low;
						}
						bit=error>>kb&1;
						{
							uint32_t r2=(range>>PROBBITS_USE)*p0;
							uint32_t mid=low+r2;
							range-=r2;
							if(!fwd)
								bit=code>=mid;
							if(bit)
								low=mid;
							else
								range=r2-1;
						}
						error|=bit<<kb;
#ifdef PRINT_ESTIMS
						esizes[kc]+=log2table[bit?(1<<PROBBITS_STORE)-pp0[model]:pp0[model]];
#endif
						{
							int32_t sh=3, rcon=1<<sh>>1, kp;

							for(kp=0;kp<NLEARNERS;++kp)
							{
								p0=pp0[kp]>>(PROBBITS_STORE-PROBBITS_USE);
#ifdef PRINT_ESTIMS
								allsizes[kc][kp]+=log2table[bit?(1<<PROBBITS_STORE)-pp0[kp]:pp0[kp]];//
#endif
								CLAMP2(p0, 1, (1<<PROBBITS_USE)-1);
								currsizes[kp]+=Shannon[bit?(1<<PROBBITS_USE)-p0:p0];
								pp0[kp]+=((!bit<<PROBBITS_STORE)-pp0[kp]+rcon)>>sh;
								++sh;
								rcon<<=1;
							}
						}
						tidx=2*tidx+bit;
					}
					if(!fwd)
						yuv[kc]=(char)(error+pred);
					
					//if(ky==10&&kx==10)//
					//	printf("");//

				//	error=yuv[kc]-pred;
				//	error=error<<1^error>>31;
					rows[0][0]=yuv[kc]-offset;
				//	rows[0][1]=(2*eW+(error<<GRBITS)+eNEEE)>>2;
#ifdef USE_L1
					{
						int32_t k, e=yuv[kc]-offset;
						e=(e>pred0)-(e<pred0);
						coeffs[NPREDS]+=e;
						for(k=0;k<NPREDS;++k)
							coeffs[k]+=e*preds[k];
					}
#endif
					{
						int kp;

						for(kp=0;kp<NLEARNERS;++kp)
							scurr[kp]+=(2*sW[kp]+currsizes[kp]+(sNEE[kp]>sNEEE[kp]?sNEE[kp]:sNEEE[kp]))>>2;
						for(kp=0;kp<NLEARNERS;++kp)
							totalsizes[kc][kp]+=currsizes[kp];
						//	totalsizes[kc][kp]+=(currsizes[kp]-totalsizes[kc][kp]+(1<<7>>1))>>7;
					}
					offset=kc?yuv[vhelpidx]:yuv[uhelpidx];
					rows[0]+=2;
					rows[1]+=2;
					rows[2]+=2;
					rows[3]+=2;
					srows[0]+=NLEARNERS;
					srows[1]+=NLEARNERS;
					srows[2]+=NLEARNERS;
					srows[3]+=NLEARNERS;
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
		free(sbuf);
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
			*(uint16_t*)streamptr=(uint16_t)(low>>16);
			streamptr+=2;
			*(uint16_t*)streamptr=(uint16_t)low;
			streamptr+=2;

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
#ifdef PRINTBITS
		printf("\n");
#endif
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
		{
			int kc, kp;

			for(kp=0;kp<NLEARNERS;++kp)
			{
				printf("%2d", kp);
				for(kc=0;kc<3;++kc)
					printf(" %16lld", totalsizes[kc][kp]>>(16+3));
				printf("\n");
			}
			esizes[0]/=8;
			esizes[1]/=8;
			esizes[2]/=8;
			printf("TYUV %16.2lf %16.2lf %16.2lf %16.2lf\n"
				, esizes[0]+esizes[1]+esizes[2]
				, esizes[0]
				, esizes[1]
				, esizes[2]
			);
			for(kp=0;kp<NLEARNERS;++kp)
			{
				printf("%2d", kp);
				for(kc=0;kc<3;++kc)
					printf(" %16.2lf", allsizes[kc][kp]/8);
				printf("\n");
			}
		}
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