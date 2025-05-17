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
	#define PROBHIST
//	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
	#define CHECK_BUFOVERFLOW
#endif
//	#define ENABLE_L1


#define L1NPREDS 8


#define PROBBITS_STORE 20
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
		printf("\n");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
	exit(1);
}
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe=0;
static void guide_save(FILE *fsrc, int iw, int ih)
{
	long fidx=ftell(fsrc);
	ptrdiff_t size=(ptrdiff_t)3*iw*ih, k;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		crash("Alloc error");
		return;
	}
	for(k=0;k<size-2;)
	{
		g_image[k++]=fgetc(fsrc);
		g_image[k++]=fgetc(fsrc);
		g_image[k++]=fgetc(fsrc);
	}
	fseek(fsrc, fidx, SEEK_SET);
}
static void guide_check(unsigned char *rgb, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(rgb, g_image+idx, 3))
	{
		crash("Guide  XY %d %d", kx, ky);
		printf("");
	}
}
static void guide_update(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx), diff;
	diff=g_image[idx+0]-image[0]; g_sqe+=diff*diff;
	diff=g_image[idx+1]-image[1]; g_sqe+=diff*diff;
	diff=g_image[idx+2]-image[2]; g_sqe+=diff*diff;
}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update(...)
#endif
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
#define RCT_COUNT 16
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
#ifdef LOUD
static const char *rctnames[]=
{
	"R G  B ",
	"R G  BG",
	"R G  BR",
	"R GR BR",
	"R GR BG",
	"R BR GB",
	"G B  RG",
	"G B  RB",
	"G BG RG",
	"G BG RB",
	"G RG BR",
	"B R  GR",
	"B R  GB",
	"B RB GB",
	"B RB GR",
	"B GB RG",
};
#endif
#ifdef PROBHIST
static int32_t probhist[1<<PROBBITS_USE]={0};
#endif
#ifdef ESTIMATE_SIZE
#define ELEVELS 256
static int32_t hist[3][ELEVELS]={0};
#endif
static uint32_t low, range, code, fwd;
static int64_t symidx, totalsyms;
static int32_t stats[3][8][256][256];
static int32_t mixer[3][256][8];
//static int32_t mixer[3][256][8*8];
static FILE *fsrc, *fdst;
static ptrdiff_t usize, csize;
static uint16_t buf[512*1024], *bufptr;
int c12_codec(int argc, char **argv)
{
	const char *srcfn, *dstfn;
	int dist;
	int32_t iw, ih, bestrct;
	ptrdiff_t res;
#ifdef LOUD
	double t=time_sec();
#endif

	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:  \"%s\"  input  output  [dist]\n"
			"  dist is optional for lossy\n"
			, argv[0]
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	dist=argc<4?1:atoi(argv[3]);
	if(dist!=1)
		CLAMP2(dist, 4, 16);

	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		crash("Cannot open \"%s\"", srcfn);
		return 1;
	}

	//parse src header
	{
		int32_t nread, c;

		c=0;
		nread=(int32_t)fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(nread!=2||(!fwd&&c!=('1'|'2'<<8)))
		{
			crash("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			c=fgetc(fsrc);
			while((uint32_t)c<=' ')c=fgetc(fsrc);//skip whitespace
			while(c=='#')//strip comments
			{
				c=fgetc(fsrc);
				while(c!='\n')c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			iw=0;
			while((uint32_t)(c-'0')<10)iw=10*iw+c-'0', c=fgetc(fsrc);
			while((uint32_t)c<=' ')c=fgetc(fsrc);//skip whitespace
			ih=0;
			while((uint32_t)(c-'0')<10)ih=10*ih+c-'0', c=fgetc(fsrc);
			while(c=='#')//strip comments
			{
				c=fgetc(fsrc);
				while(c!='\n')c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			nread=(int32_t)fread(&c, 1, 4, fsrc);
			if(c!=('2'|'5'<<8|'5'<<16|'\n'<<24))
			{
				crash("Unsupported file \"%s\"", srcfn);
				return 1;
			}

		}
		else
		{
			int nread=0;

			nread+=(int)fread(&iw, 1, 4, fsrc);
			nread+=(int)fread(&ih, 1, 4, fsrc);
			if(nread!=8)
				crash("Unsupported file \"%s\"", srcfn);
			nread=0;
			bestrct=0;
			nread+=(int)fread(&bestrct, 1, 1, fsrc);
			dist=0;
			nread+=(int)fread(&dist, 1, 1, fsrc);
			if(nread!=2)
				crash("Unsupported file \"%s\"", srcfn);

			fread((unsigned char*)&code+2, 1, 2, fsrc);
			fread((unsigned char*)&code+0, 1, 2, fsrc);
		}
		if(iw<1||ih<1)
		{
			crash("Unsupported file \"%s\"", srcfn);
			return 1;
		}
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	csize=0;

	low=0;
	range=0xFFFFFFFF;
	bufptr=fwd?buf:buf+_countof(buf);

	//write dst header
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		crash("Cannot open \"%s\" for writing", fdst);
		return 1;
	}
	if(fwd)
	{
		csize+=fwrite("12", 1, 2, fdst);
		csize+=fwrite(&iw, 1, 4, fdst);
		csize+=fwrite(&ih, 1, 4, fdst);
	}
	else
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);

	//analysis
	if(fwd)
	{
		int64_t counters[OCH_COUNT]={0}, minerr=0;
		int prev[OCH_COUNT]={0};
		ptrdiff_t k;
		long fidx=ftell(fsrc);
		unsigned char rgb[4]={0};

		for(k=0;k<res;++k)
		{
			int r, g, b, rg, gb, br;

			fread(rgb, 1, 3, fsrc);

			r=rgb[0];
			g=rgb[1];
			b=rgb[2];
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
		fseek(fsrc, fidx, SEEK_SET);
		guide_save(fsrc, iw, ih);
		{
			int kt;

#ifdef ESTIMATE_SIZE
			printf("\"%s\"  WH %d %d  dist %d\n", srcfn, iw, ih, dist);
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			bestrct=0;
			for(kt=0;kt<RCT_COUNT;++kt)
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
				printf("RCT%02d %s %16lld%s\n", kt, rctnames[kt], currerr, kt==bestrct?" <-":"");
#endif
			}
		}
		fwrite(&bestrct, 1, 1, fdst);
		fwrite(&dist, 1, 1, fdst);
	}
	{
		int
			yidx=rct_indices[bestrct][3+0],
			uidx=rct_indices[bestrct][3+1],
			vidx=rct_indices[bestrct][3+2],
			uhelpidx=rct_indices[bestrct][6+0],
			vhelpidx=rct_indices[bestrct][6+1];
		const int32_t psize=(int32_t)sizeof(int16_t[4*3])*(iw+16);//4 padded rows * 3 channels * {pixels}
		int32_t paddedwidth=iw+16;
		int16_t *pixels=0;
	//	int32_t coeffs[L1NPREDS+1]={0};
	//	int32_t invdist=((1<<16)+dist-1)/dist;
		int32_t ky, kx, idx;
		unsigned char rgb[3]={0};

		pixels=(int16_t*)malloc(psize);
		if(!pixels)
		{
			crash("Alloc error\n");
			return 1;
		}
		memset(pixels, 0, psize);

		//FILLMEM((uint32_t*)stats, (1<<PROBBITS_STORE>>1), sizeof(stats), sizeof(int32_t));
		memset(stats, 0, sizeof(stats));
		FILLMEM((int32_t*)mixer, (1<<PROBBITS_STORE)/8, sizeof(mixer), sizeof(int32_t));
		//memset(mixer, 0, sizeof(mixer));
#ifdef ESTIMATE_SIZE
		memset(hist, 0, sizeof(hist));
#endif
		totalsyms=usize;
		symidx=0;
		for(ky=0, idx=0;ky<ih;++ky)
		{
			unsigned char yuv[4]={0};
			int16_t
				*NNNptr		=pixels+(paddedwidth*((ky-3)&3)+8)*3,
				*NNptr		=pixels+(paddedwidth*((ky-2)&3)+8)*3,
				*Nptr		=pixels+(paddedwidth*((ky-1)&3)+8)*3,
				*currptr	=pixels+(paddedwidth*((ky-0)&3)+8)*3;
			for(kx=0;kx<iw;++kx, ++idx)
			{
				int offset=0, kc;

				if(fwd)
				{
					fread(rgb, 1, 3, fsrc);
					yuv[0]=rgb[yidx];
					yuv[1]=rgb[uidx];
					yuv[2]=rgb[vidx];
					yuv[2]-=yuv[vhelpidx]-128;
					yuv[1]-=yuv[uhelpidx]-128;
				}
				for(kc=0;kc<3;++kc)
				{
					int32_t
						NNN	=NNNptr		[+0*3],
						NNE	=NNptr		[+1*3],
						NN	=NNptr		[+0*3],
						NW	=Nptr		[-1*3],
						N	=Nptr		[+0*3],
						NE	=Nptr		[+1*3],
						NEE	=Nptr		[+2*3],
						NEEE	=Nptr		[+3*3],
						NEEEE	=Nptr		[+4*3],
						WWWW	=currptr	[-4*3],
						WWW	=currptr	[-3*3],
						WW	=currptr	[-2*3],
						W	=currptr	[-1*3];
					//int32_t pred, vmax, vmin;
					int32_t error;
					int32_t *statsptrs[8], *mixptr;
					int32_t kb, tidx, k2;
					int32_t preds[]=
					{
						N,
						W,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
						N+W-NW,
						W+NE-N,
						N+NE-NNE,
						(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2,
					};
#ifdef ENABLE_L1
					int32_t pred0;
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
#else
					//pred=N+W-NW; vmax=N; vmin=W;
					//if(N<W)vmin=N, vmax=W;
					//CLAMP2(pred, vmin, vmax);
#endif
					//pred+=offset;
					//CLAMP2(pred, 0, 255);
					//if(ky==10&&kx==10)//
					//	printf("");
					statsptrs[0]=stats[kc][0][preds[0]&255];
					statsptrs[1]=stats[kc][1][preds[1]&255];
					statsptrs[2]=stats[kc][2][preds[2]&255];
					statsptrs[3]=stats[kc][3][preds[3]&255];
					statsptrs[4]=stats[kc][4][preds[4]&255];
					statsptrs[5]=stats[kc][5][preds[5]&255];
					statsptrs[6]=stats[kc][6][preds[6]&255];
					statsptrs[7]=stats[kc][7][preds[7]&255];
					mixptr=mixer[kc][(N+W)>>1];
					if(fwd)
					{
						//if(dist>1)
						//{
						//	error=yuv[kc]-pred;
						//	error=(error*invdist>>16)+((uint32_t)error>>31);
						//	yuv[kc]=error*dist+pred;
						//	CLAMP2(yuv[kc], 0, 255);
						//}
						//else
							error=yuv[kc];
#ifdef ESTIMATE_SIZE
						++hist[kc][(error+128)&255];
#endif
					}
					else
						error=0;
					for(kb=7, tidx=1;kb>=0;--kb)
					{
						int32_t bit=error>>kb&1;
						int32_t *probptrs[]=
						{
							statsptrs[0]+tidx,
							statsptrs[1]+tidx,
							statsptrs[2]+tidx,
							statsptrs[3]+tidx,
							statsptrs[4]+tidx,
							statsptrs[5]+tidx,
							statsptrs[6]+tidx,
							statsptrs[7]+tidx,
						};
						int32_t probs[]=
						{
							*probptrs[0],
							*probptrs[1],
							*probptrs[2],
							*probptrs[3],
							*probptrs[4],
							*probptrs[5],
							*probptrs[6],
							*probptrs[7],
						};
						int64_t mp=0;
						int32_t p0;
						for(k2=0;k2<8;++k2)
							mp+=(int64_t)mixptr[k2]*probs[k2];
						p0=(int32_t)(mp>>29);
						//p0=(int32_t)(mp>>33);
						p0+=1<<PROBBITS_USE>>1;
						CLAMP2(p0, 1, (1<<PROBBITS_USE)-1);
#ifdef PROBHIST
						if(fwd)
							++probhist[p0];
#endif
						if(range<(1<<PROBBITS_USE))
						{
							if(fwd)
							{
#ifdef CHECK_BUFOVERFLOW
								if(csize>usize+1024)
								{
#ifdef PROBHIST
									for(int k=0;k<(1<<PROBBITS_USE);++k)//
									{
										if(probhist[k])
											printf("0x%04X %8d\n", k, probhist[k]);
									}
#endif
									crash("ERROR at %d/%d  inflation %8.4lf%%\n"
										, (int32_t)symidx
										, (int32_t)totalsyms
										, 100.*totalsyms/symidx
									);
									return 1;
								}
#endif
								if(bufptr>=buf+_countof(buf))
								{
									csize+=fwrite(buf, 1, sizeof(buf), fdst);
									bufptr=buf;
								}
								*bufptr++=(uint16_t)(low>>16);
							}
							else
							{
								if(bufptr>=buf+_countof(buf))
								{
									fread(buf, 1, sizeof(buf), fsrc);
									bufptr=buf;
								}
								code=code<<16|*bufptr++;
							}
							low<<=16;
							range=range<<16|0xFFFF;
							if(range>~low)
								range=~low;
						}
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
						{
							//int32_t e=2*bit-1;//L1 = abs(bit-p0)		dL/dwi = -sgn(bit-p0) * pi
							//int32_t e=(bit<<PROBBITS_USE)-p0;//L2 = 1/2 (bit-p0)^2		dL/dwi = -(bit-p0) * pi
							//int32_t e=-(1<<PROBBITS_USE)/((bit<<PROBBITS_USE)-p0);//L = -log2(bit?1-p0:p0)	dL/dwi = -pi / ((bit-p0) * ln2)
							//e=-e;
							for(k2=0;k2<8;++k2)
							{
								int32_t tp0=probs[k2];
							//	int64_t tm0=mixptr[k2];
								int32_t tp1=((!bit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1)-tp0+(1<<7>>1))>>7;
								//int64_t tp1=(int64_t)e*mixptr[k2]>>7;
							//	int64_t tm1=(int64_t)e*probs[k2]>>12;
								tp1+=tp0;
							//	tm1+=tm0;
								//CLAMP2(tp1, -(1<<PROBBITS_STORE), (1<<PROBBITS_STORE));
							//	CLAMP2(tm1, -(1<<PROBBITS_STORE), (1<<PROBBITS_STORE));
								*probptrs[k2]=tp1;
							//	mixptr[k2]=(int32_t)tm1;
							}
						}
						tidx=2*tidx+bit;
						//mixptr+=8;
					}
					if(!fwd)
					{
						//if(dist>1)
						//{
						//	yuv[kc]=error*dist+pred;
						//	CLAMP2(yuv[kc], 0, 255);
						//}
						//else
							yuv[kc]=(unsigned char)error;
					}
					currptr[0]=yuv[kc]-offset;
					{
						int32_t e=currptr[0], k;
						int32_t pred=0;
						for(k=0;k<L1NPREDS;++k)
							pred+=mixptr[k]*preds[k];
						pred+=1<<19>>1;
						pred>>=19;
						e=(e>pred)-(e<pred);
						for(k=0;k<L1NPREDS;++k)
							mixptr[k]+=e*preds[k];
					}
#ifdef ENABLE_L1
					{
						int32_t e=currptr[0], k;
						e=(e>pred0)-(e<pred0);
						coeffs[L1NPREDS]+=e;
						for(k=0;k<L1NPREDS;++k)
							coeffs[k]+=e*preds[k];
					}
#endif
					//offset=kc?yuv[vhelpidx]:yuv[uhelpidx];
					++NNNptr;
					++NNptr;
					++Nptr;
					++currptr;
					++symidx;
				}
				if(!fwd)
				{
					yuv[1]+=yuv[uhelpidx]-128;
					yuv[2]+=yuv[vhelpidx]-128;
					rgb[yidx]=yuv[0];
					rgb[uidx]=yuv[1];
					rgb[vidx]=yuv[2];
					fwrite(rgb, 1, 3, fdst);
#ifdef ENABLE_GUIDE
					if(dist>1)
						guide_update(rgb, kx, ky);
					else
						guide_check(rgb, kx, ky);
#endif
				}
			}
		}
#ifdef PROBHIST
		if(fwd)
		{
			for(int k=0;k<(1<<PROBBITS_USE);++k)//
			{
				if(probhist[k])
					printf("0x%04X %8d\n", k, probhist[k]);
			}
		}
#endif
		free(pixels);
	}
	if(fwd)
	{
		if(bufptr>buf)
			csize+=fwrite(buf, 1, (size_t)bufptr-(size_t)buf, fdst);
		csize+=fwrite((unsigned char*)&low+2, 1, 2, fdst);
		csize+=fwrite((unsigned char*)&low+0, 1, 2, fdst);
	}
	fclose(fsrc);
	fclose(fdst);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		ptrdiff_t srcsize=usize, dstsize=csize;
#ifdef ESTIMATE_SIZE
		{
			double csizes[3]={0};
			int kc, ks;

			for(kc=0;kc<3;++kc)
			{
				int32_t *currhist=hist[kc], sum, count;
				double norm, e;

				for(sum=0, count=0, ks=0;ks<ELEVELS;++ks)
				{
					sum+=currhist[ks];
					count+=currhist[ks]!=0;
				}
				norm=1./sum;
				for(e=0, ks=0;ks<ELEVELS;++ks)
				{
					int32_t freq=currhist[ks];
					if(freq)
					{
						int32_t p=(int32_t)((int64_t)freq*(0x1000LL-count)/sum)+1;
						e-=freq*log(p*(1./0x1000));
						if(!isfinite(csizes[kc]))
							crash("");
					}
				}
				e*=1/(8*M_LN2);//base e -> base 256
				csizes[kc]=e;
			}
			printf("UTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf\n"
				, (double)usize
				, csizes[0]+csizes[1]+csizes[2]
				, csizes[0]
				, csizes[1]
				, csizes[2]
			);
		}
#endif
		printf("%9td->%9td  %8.4lf%%  %12.6lf\n"
			, srcsize
			, dstsize
			, 100.*dstsize/srcsize
			, (double)srcsize/dstsize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
	);
#endif
#ifdef ENABLE_GUIDE
	if(!fwd&&dist>1)
	{
		double rmse=sqrt(g_sqe/((double)3*iw*ih)), psnr=20*log10(255/rmse);
		printf("BPP %8.4lf  RMSE %12.6lf PSNR %12.6lf\n", 8.*csize/usize, rmse, psnr);
	}
#endif
	(void)time_sec;
	return 0;
}
