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
//#ifdef _DEBUG
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
//	#define PRINTBITS
//	#define PRINTGR
//#endif
#endif


#define L1NPREDS 8


#define GRBITS 6
#define GRLIMIT 24
#define PROBBITS_STORE 16
#define PROBBITS_USE 9

#define MIXBITS 16


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
#ifdef _DEBUG
static void debugcrash(void)
{
	printf("CRASH\n");
	exit(1);
}
#else
#define debugcrash()
#endif
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		debugcrash();
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		debugcrash();
		printf("");
	}
}
static void guide_update(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx), diff;
	diff=g_image[idx+0]-image[idx+0]; g_sqe+=diff*diff;
	diff=g_image[idx+1]-image[idx+1]; g_sqe+=diff*diff;
	diff=g_image[idx+2]-image[idx+2]; g_sqe+=diff*diff;
}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update(...)
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
//static uint16_t stats[3][256][256];
//static uint16_t stats[3][8][256][256];
static uint16_t stats[3][256][8][GRLIMIT+8];
typedef struct _ACState
{
	uint32_t low, range, code, fwd;
	unsigned char *ptr, *end;
	uint32_t symidx, totalsyms;
} ACState;
AWM_INLINE void codebit(ACState *ac, uint16_t *pp0a, int32_t *bit)
{
	uint32_t r2, mid;
	int32_t p0a=*pp0a;
	int32_t p0=p0a>>(PROBBITS_STORE-PROBBITS_USE);
	if(ac->range<(1<<PROBBITS_USE))
	{
		if(ac->ptr>=ac->end)
		{
			printf("ERROR at %d/%d  inflation %8.4lf%%\n"
				, (int32_t)ac->symidx
				, (int32_t)ac->totalsyms
				, 100.*ac->totalsyms/ac->symidx
			);
			debugcrash();
			exit(1);
		}
		if(ac->fwd)
			*(uint16_t*)ac->ptr=(uint16_t)(ac->low>>16);
		else
			ac->code=ac->code<<16|*(uint16_t*)ac->ptr;
		ac->ptr+=2;
		ac->low<<=16;
		ac->range=ac->range<<16|0xFFFF;
		if(ac->range>~ac->low)
			ac->range=~ac->low;
	}
	r2=(ac->range>>PROBBITS_USE)*(p0+(p0<(1<<PROBBITS_USE>>1)));
	mid=ac->low+r2;
	ac->range-=r2;
	if(!ac->fwd)
		*bit=ac->code>=mid;
	if(*bit)
		ac->low=mid;
	else
		ac->range=r2-1;
	int32_t truth=!*bit<<PROBBITS_STORE;
	*pp0a=p0a+((int32_t)(truth-p0a)>>6);
}
int c12_codec(int argc, char **argv)
{
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:  \"%s\"  input  output  [dist]\n"
			"  dist is optional for lossy\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int dist=argc<4?1:atoi(argv[3]);
	if(dist!=1)
		CLAMP2(dist, 4, 16);
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *ptr=0, *imptr=0;
	int bestrct=0;
	ACState ac=
	{
		0, 0xFFFFFFFF, 0, 0,
		0, 0,
	};
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef PRINTGR
	long long gr_symlen=0, gr_bypsum=0;
#endif

	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		debugcrash();
		return 1;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('1'|'2'<<8))
		{
			printf("Unsupported file  \"%s\"\n", srcfn);
			debugcrash();
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	if(fwd)
	{
		if(*srcptr++!='\n')
		{
			printf("Unsupported file\n");
			debugcrash();
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
			printf("Unsupported file\n");
			debugcrash();
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
			printf("Unsupported file\n");
			debugcrash();
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
		printf("Invalid file\n");
		debugcrash();
		free(srcbuf);
		return 1;
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	dstbuf=(unsigned char*)malloc(usize);
	if(!dstbuf)
	{
		printf("Alloc error");
		debugcrash();
		free(srcbuf);
		return 1;
	}
	ac.fwd=fwd;
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
		ptr=dstbuf;
		*ptr++=bestrct;
		*ptr++=dist;

		ac.ptr=ptr;
		ac.end=dstbuf+usize;
	}
	else
	{
		imptr=dstbuf;
		ptr=srcptr;
		bestrct=*ptr++;
		dist=*ptr++;

		ac.code=ac.code<<16|*(uint16_t*)ptr;
		ptr+=2;
		ac.code=ac.code<<16|*(uint16_t*)ptr;
		ptr+=2;
		ac.ptr=ptr;
		ac.end=srcend;

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
		int32_t paddedwidth=iw+16;

		int32_t coeffs[L1NPREDS+1]={0};

		int32_t invdist=((1<<16)+dist-1)/dist;

		psize=(int32_t)sizeof(int16_t[4*3*2])*(iw+16);//4 padded rows * 3 channels * {pixels, nbypass}
		pixels=(int16_t*)malloc(psize);
		if(!pixels)
		{
			printf("Alloc error\n");
			debugcrash();
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
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
						W	=rows[0][0-1*3*2],
						eNEEE	=rows[1][1-3*3*2],
						eW	=rows[0][1-1*3*2];
					int32_t pred, vmax, vmin, pred0;
					int32_t error;
					int32_t nbypass, nzeros=-1, bypass=0;
					int32_t tidx=0;
					uint16_t *statsptr;
					int32_t bit;

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
					pred=coeffs[L1NPREDS];
					{
						int j=0;
						for(;j<L1NPREDS;++j)
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

					pred+=offset;
					CLAMP2(pred, -128, 127);

					{
						float fval=(float)(eW+1);
						size_t addr=(size_t)&fval;
						int32_t bits=*(int32_t*)addr;
						nbypass=(bits>>23)-127-GRBITS;
						if(nbypass<0)
							nbypass=0;
					}
#if 0
					statsptr=stats[kc][eW>>GRBITS];
					//statsptr=stats[kc][nbypass][(pred+128)&255];
					int step=1<<nbypass;
					if(fwd)
					{
						if(dist>1)
						{
							int e2=yuv[kc]-pred;
							if(e2<0)
								e2+=dist-1;
							e2=e2*invdist>>16;
							error=e2<<1^e2>>31;
							yuv[kc]=e2*dist+pred;
							CLAMP2(yuv[kc], -128, 127);
						}
						else
						{
							error=(char)(yuv[kc]-pred);
							error=error<<1^error>>31;
						}
						nzeros=error>>nbypass;
						bypass=error&((1<<nbypass)-1);
#ifdef ESTIMATE_SIZE
						++hist[kc][error];
#endif
					}
					else
						error=0;
					tidx=0;
					do
					{
						int tidx2=tidx+step;
						bit=tidx2<error;
						codebit(&ac, statsptr+tidx, &bit);
						tidx=tidx2;
					}while(bit);
					tidx-=step;
					while(step)
					{
						int floorhalf=step>>1;
						int mid=tidx+floorhalf;
						bit=mid<error;
						codebit(&ac, statsptr+mid, &bit);
						tidx+=(step-floorhalf)&-bit;
						step=floorhalf;
					}
					if(!fwd)
					{
						error=tidx;
						error=error>>1^-(error&1);
						if(dist>1)
						{
							yuv[kc]=error*dist+pred;
							CLAMP2(yuv[kc], -128, 127);
						}
						else
							yuv[kc]=(char)(error+pred);
					}
#ifdef _DEBUG
					else if(error!=tidx)
						debugcrash();
#endif
#endif
#if 1
					statsptr=stats[kc][(pred+128)&255][nbypass];
					if(fwd)
					{
						if(dist>1)
						{
							int e2=yuv[kc]-pred;
							if(e2<0)
								e2+=dist-1;
							e2=e2*invdist>>16;
							error=e2<<1^e2>>31;
							yuv[kc]=e2*dist+pred;
							CLAMP2(yuv[kc], -128, 127);
						}
						else
						{
							error=(char)(yuv[kc]-pred);
							error=error<<1^error>>31;
						}
						nzeros=error>>nbypass;
						bypass=error&((1<<nbypass)-1);
#ifdef ESTIMATE_SIZE
						++hist[kc][error];
#endif
#ifdef PRINTGR
						{
							float fval=(float)(error+1);
							size_t addr=(size_t)&fval;
							int32_t bits=*(int32_t*)addr;
							bits=(bits>>23)-127;
							gr_bypsum+=nbypass;
							gr_symlen+=bits;
						}
#endif
					}
					else
						error=0;
					do
					{
						bit=nzeros--<=0;
						codebit(&ac, statsptr+tidx, &bit);
#ifdef PRINTBITS
						if(fwd&&(unsigned)(idx-(usize>>2))<1000)printf("%c", '0'+bit);//
#endif
						++tidx;
						if(tidx==GRLIMIT)
						{
							tidx=1;
							nbypass=8;
							if(fwd)
								bypass=error;
							break;
						}
					}while(!bit);
					{
						int32_t kb=nbypass-1;

						for(;kb>=0;--kb)
						{
							bit=bypass>>kb&1;
							codebit(&ac, statsptr+GRLIMIT+8-nbypass+kb, &bit);
							bypass|=bit<<kb;
#ifdef PRINTBITS
							if(fwd&&(unsigned)(idx-(usize>>2))<1000)printf("%c", '0'+bit);//
#endif
						}
					}
					if(!fwd)
					{
						error=(tidx-1)<<nbypass|bypass;
						error=error>>1^-(error&1);
						if(dist>1)
						{
							yuv[kc]=error*dist+pred;
							CLAMP2(yuv[kc], -128, 127);
						}
						else
							yuv[kc]=(char)(error+pred);
					}
#endif

					{
						int32_t k, e=yuv[kc]-offset;
						e=(e>pred0)-(e<pred0);
						coeffs[L1NPREDS]+=e;
						for(k=0;k<L1NPREDS;++k)
							coeffs[k]+=e*preds[k];
					}
					error=abs(yuv[kc]-pred);
					rows[0][0]=yuv[kc]-offset;
					rows[0][1]=(2*eW+(error<<GRBITS)+eNEEE)>>2;
					if(kc)
						offset=yuv[vhelpidx];
					else
						offset=yuv[uhelpidx];
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
#ifdef ENABLE_GUIDE
					if(dist>1)
						guide_update(dstbuf, kx, ky);
					else
						guide_check(dstbuf, kx, ky);
#endif
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
			printf("Cannot open \"%s\" for writing\n", dstfn);
			debugcrash();
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		if(fwd)
		{
			*(uint16_t*)ac.ptr=(uint16_t)(ac.low>>16);
			ac.ptr+=2;
			*(uint16_t*)ac.ptr=(uint16_t)ac.low;
			ac.ptr+=2;

			csize=ac.ptr-dstbuf;

			dstsize+=fwrite("12", 1, 2, fdst);
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
#ifdef PRINTGR
		printf("%12lld bypass %12.6lf bits\n", gr_bypsum, (double)gr_bypsum/usize);
		printf("%12lld symlen %12.6lf bits\n", gr_symlen, (double)gr_symlen/usize);
#endif
		printf("%9td->%9td  %8.4lf%%  %12.6lf\n"
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
#ifdef ENABLE_GUIDE
	if(!fwd&&dist>1)
	{
		double rmse=sqrt(g_sqe/((double)3*iw*ih)), psnr=20*log10(255/rmse);
		printf("RMSE %12.6lf PSNR %12.6lf\n", rmse, psnr);
	}
#endif
	(void)time_sec;
	return 0;
}
