#include"codec.h"
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
//#include"blist.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define ENABLE_VAL


	#define USE_MIXIN
//	#define ENABLE_RECIPROCAL
//	#define ENABLE_SSE


//	#define AC_VALIDATE
//	#define AC_IMPLEMENTATION
#include"entropy.h"
#ifdef ENABLE_VAL
typedef struct _ValInfo
{
	int cdf, freq;
	unsigned long long low0, range0, low1, range1;
} ValInfo;
static long long validx=0;
static ArrayHandle valdata=0;
static void val_enc(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long low1, unsigned long long range1)
{
	if(!valdata)
	{
		ARRAY_ALLOC(ValInfo, valdata, 0, 0, 0, 0);
		if(!valdata)
		{
			LOG_ERROR("Alloc error");
			return;
		}
	}
	ValInfo info={cdf, freq, low0, range0, low1, range1};
	ARRAY_APPEND(valdata, &info, 1, 1, 0);
}
static int val_dec(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long code0, unsigned long long low1, unsigned long long range1, unsigned long long code1)
{
	ValInfo *info=(ValInfo*)array_at(&valdata, validx);
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
#define VAL_DEC(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1, ...) if(val_dec(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1))LOG_ERROR(__VA_ARGS__)
#else
#define val_enc(...)
#define val_dec(...)
#define VAL_DEC(...)
#endif
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
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("XY %d %d  RGB 0x%08X 0x%08X OG", kx, ky, *(unsigned*)(image+idx), *(unsigned*)(g_image+idx));
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

	OCH_COUNT=6,
} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},
};
#ifdef ENABLE_SSE
#define SSE_FBITS 6
#define SSE_LBITS 8
static int sse[3][64][64];
//static int sse[3][256];
#endif
#ifdef USE_MIXIN
#define RATE 9
static unsigned short hist[3][16][256];
static unsigned short mixin[512];
#else
static int hist[3][16][257];
#ifdef ENABLE_RECIPROCAL
static unsigned long long divtable[256];
#endif
#endif
int c10_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
	const int headersize=3+4+4;
	ptrdiff_t usize=0, overhead=0, csize=0, streamsize=0;
	unsigned char *buf=0, *image=0, *stream=0;
	int iw=0, ih=0, fwd=0;
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		size_t tag=0, nread;
		fread(&tag, 1, 3, fsrc);
		if(tag==('P'|'6'<<8|'\n'<<16))
		{
			fwd=1;
			nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				LOG_ERROR("Unsupported file");
				fclose(fsrc);
				return 1;
			}
			char str[6]={0};
			fread(str, 1, 5, fsrc);
			if(strcmp(str, "\n255\n"))
			{
				LOG_ERROR("Unsupported file");
				fclose(fsrc);
				return 1;
			}
			usize=(ptrdiff_t)3*iw*ih;
			overhead=(ptrdiff_t)iw*ih;
			buf=(unsigned char*)malloc(usize+overhead);
			if(!buf)
			{
				LOG_ERROR("Alloc error");
				fclose(fsrc);
				return 1;
			}
			memset(buf, 128, overhead);
			stream=buf;
			image=buf+overhead;
			nread=fread(image, 1, usize, fsrc);
			if(nread!=usize)
				printf("Warning: source image truncated at %td/%td\n", nread, usize);
			guide_save(image, iw, ih);
		}
		else if(tag==('1'|'0'<<8|'\n'<<16))
		{
			nread=fread(&iw, 1, 4, fsrc);
			nread+=fread(&ih, 1, 4, fsrc);
			if(nread!=8)
			{
				LOG_ERROR("Unsupported file");
				fclose(fsrc);
				return 1;
			}
			csize=get_filesize(srcfn);
			streamsize=csize-headersize;
			usize=(ptrdiff_t)3*iw*ih;
			overhead=(ptrdiff_t)iw*ih;
			buf=(unsigned char*)malloc(usize+overhead);
			if(!buf)
			{
				LOG_ERROR("Alloc error");
				fclose(fsrc);
				return 1;
			}
			memset(buf, 128, overhead);
			image=buf+(ptrdiff_t)9*iw;
			stream=buf+usize+overhead-streamsize;
			nread=fread(stream, 1, streamsize, fsrc);
			if(nread!=streamsize)
				printf("Warning: source stream truncated at %td/%td\n", nread, streamsize);
		}
		else
		{
			LOG_ERROR("Unsupported file");
			fclose(fsrc);
			return 1;
		}
		fclose(fsrc);
	}
	unsigned char *streamptr=stream;
	unsigned long long low=0, range=0xFFFFFFFF, code=0;
	int bestrct=0;
	if(fwd)//analysis
	{
		unsigned char *ptr=image, *end=image+usize;
		unsigned long long counters[6]={0};
		//int prev3[6]={0};
		int prev2[6]={0}, prev[6]={0};
		while(ptr<end)
		{
			int
				r=ptr[0],
				g=ptr[1],
				b=ptr[2],
				rg=r-g,
				gb=g-b,
				br=b-r;
			ptr+=3;
			//counters[0]+=abs(r -3*(prev[0]-prev2[0])-prev3[0]);
			//counters[1]+=abs(g -3*(prev[1]-prev2[1])-prev3[1]);
			//counters[2]+=abs(b -3*(prev[2]-prev2[2])-prev3[2]);
			//counters[3]+=abs(rg-3*(prev[3]-prev2[3])-prev3[3]);
			//counters[4]+=abs(gb-3*(prev[4]-prev2[4])-prev3[4]);
			//counters[5]+=abs(br-3*(prev[5]-prev2[5])-prev3[5]);
			counters[0]+=abs(r -2*prev[0]+prev2[0]);
			counters[1]+=abs(g -2*prev[1]+prev2[1]);
			counters[2]+=abs(b -2*prev[2]+prev2[2]);
			counters[3]+=abs(rg-2*prev[3]+prev2[3]);
			counters[4]+=abs(gb-2*prev[4]+prev2[4]);
			counters[5]+=abs(br-2*prev[5]+prev2[5]);
			//counters[0]+=abs(r -prev[0]);
			//counters[1]+=abs(g -prev[1]);
			//counters[2]+=abs(b -prev[2]);
			//counters[3]+=abs(rg-prev[3]);
			//counters[4]+=abs(gb-prev[4]);
			//counters[5]+=abs(br-prev[5]);
			//prev3[0]=prev2[0];
			//prev3[1]=prev2[1];
			//prev3[2]=prev2[2];
			//prev3[3]=prev2[3];
			//prev3[4]=prev2[4];
			//prev3[5]=prev2[5];
			prev2[0]=prev[0];
			prev2[1]=prev[1];
			prev2[2]=prev[2];
			prev2[3]=prev[3];
			prev2[4]=prev[4];
			prev2[5]=prev[5];
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
		low+=range*bestrct>>4;
		range=(range>>4)-1;
	}
	else
	{
		code=*(unsigned*)streamptr;//big-endian
		streamptr+=4;
		code=code<<4|*(unsigned*)streamptr;
		streamptr+=4;

		bestrct=(int)((((code-low)<<4)+(1LL<<4)-1)/range);
		low+=range*bestrct>>4;
		range=(range>>4)-1;
	}
	short *ebuf=(short*)malloc((iw+16LL)*sizeof(short[3]));
	if(!ebuf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(ebuf, 0, (iw+16LL)*sizeof(short[3]));
	int
		yidx=rct_indices[bestrct][3+0],
		uidx=rct_indices[bestrct][3+1],
		vidx=rct_indices[bestrct][3+2],
		uhelpidx=rct_indices[bestrct][6+0],
		vhelpidx=rct_indices[bestrct][6+1];
	unsigned char *NNptr=image-(ptrdiff_t)6*iw;
	unsigned char *Nptr=image-(ptrdiff_t)3*iw;
	unsigned char *ptr=image;
#ifdef USE_MIXIN
	for(int k=0;k<256;++k)
		hist[0][0][k]=k*((0x10000-256)/256);
	memfill((unsigned short*)hist+256, hist, sizeof(hist)-sizeof(short[256]), sizeof(short[256]));
	memset(mixin, 0, sizeof(mixin));
	FILLMEM(mixin+256, (0x10000-256)>>RATE, sizeof(mixin)-sizeof(short[256]), sizeof(short));
#else
	FILLMEM((int*)hist, 1, sizeof(int[256]), sizeof(int));
	hist[0][0][256]=256;
	memfill((int*)hist+257, hist, sizeof(hist)-sizeof(int[257]), sizeof(int[257]));
#endif
#ifdef ENABLE_RECIPROCAL
	{
		unsigned long long rem=0;
		for(int k=2;k<256;++k)
			UDIV128(divtable[k], rem, 1, k-1ULL, k);
		(void)rem;
	}
#endif
#ifdef ENABLE_SSE
	memset(sse, 0, sizeof(sse));
#endif
	for(int ky=0;ky<ih;++ky)
	{
		short *eptr=ebuf+(ptrdiff_t)8*3;
		int vmax[3]={0}, vmin[3]={0}, preds[3]={0};
		char NW[4]={0}, N[4]={0}, WWW[4]={0}, WW[4]={0}, W[4]={0}, yuv[4]={0};
		//int errors3[3]={0};
		int errors2[3]={0}, errors[3]={0}, deltas[3]={0};
		int sym, den, cdf, freq;
		(void)den;
#ifdef ENABLE_RECIPROCAL
		unsigned long long invden=0, hi0=0, hi1=0;
#endif
#ifdef USE_MIXIN
		unsigned short *curr_hist[3]={0};
#else
		int *curr_hist[3]={0};
#endif
		N[0]=Nptr[yidx]-128;
		N[1]=Nptr[uidx]-128;
		N[2]=Nptr[vidx]-128;
		N[2]-=N[vhelpidx];
		N[1]-=N[uhelpidx];
		for(int kx0=iw;kx0--;)
		{
#ifdef _DEBUG
			int kx=iw-1-kx0;
#endif
			char NN[4]=
			{
				NNptr[yidx+0]-128,
				NNptr[uidx+0]-128,
				NNptr[vidx+0]-128,
			};
			NN[2]-=NN[vhelpidx];
			NN[1]-=NN[uhelpidx];
			char NE[4]=
			{
				Nptr[yidx+3]-128,
				Nptr[uidx+3]-128,
				Nptr[vidx+3]-128,
			};
			NE[2]-=NE[vhelpidx];
			NE[1]-=NE[uhelpidx];

			vmax[0]=N[0], vmin[0]=W[0];
			vmax[1]=N[1], vmin[1]=W[1];
			vmax[2]=N[2], vmin[2]=W[2];
			if(N[0]<W[0])vmin[0]=N[0], vmax[0]=W[0];
			if(N[1]<W[1])vmin[1]=N[1], vmax[1]=W[1];
			if(N[2]<W[2])vmin[2]=N[2], vmax[2]=W[2];
			//if(vmin[0]>NE[0])vmin[0]=NE[0];
			//if(vmax[0]<NE[0])vmax[0]=NE[0];
			//if(vmin[1]>NE[1])vmin[1]=NE[1];
			//if(vmax[1]<NE[1])vmax[1]=NE[1];
			//if(vmin[2]>NE[2])vmin[2]=NE[2];
			//if(vmax[2]<NE[2])vmax[2]=NE[2];
			preds[0]=N[0]+W[0]-NW[0];
			preds[1]=N[1]+W[1]-NW[1];
			preds[2]=N[2]+W[2]-NW[2];
			//preds[0]=(3*(N[0]+W[0])-2*NW[0])>>2;
			//preds[1]=(3*(N[1]+W[1])-2*NW[1])>>2;
			//preds[2]=(3*(N[2]+W[2])-2*NW[2])>>2;
			//preds[0]=(4*(N[0]+W[0])+NE[0]-NW[0])>>3;
			//preds[1]=(4*(N[1]+W[1])+NE[1]-NW[1])>>3;
			//preds[2]=(4*(N[2]+W[2])+NE[2]-NW[2])>>3;
			//errors[0]=abs(errors[0]);
			//errors[1]=abs(errors[1]);
			//errors[2]=abs(errors[2]);
			errors2[0]+=((abs(errors[0])<<3)-errors2[0])>>2;
			errors2[1]+=((abs(errors[1])<<3)-errors2[1])>>2;
			errors2[2]+=((abs(errors[2])<<3)-errors2[2])>>2;
			//if(errors2[0]<0||errors2[1]<0||errors2[2]<0)
			//	LOG_ERROR("");
			int ctx[]=
			{
				//abs(N[0]-W[0])+abs(N[0]-NW[0])+abs(NE[0]-N[0])+abs(W[0]-NW[0])+((abs(N[0]-NN[0])+abs(W[0]-WW[0]))>>1)+1,
				//abs(N[1]-W[1])+abs(N[1]-NW[1])+abs(NE[1]-N[1])+abs(W[1]-NW[1])+((abs(N[1]-NN[1])+abs(W[1]-WW[1]))>>1)+1,
				//abs(N[2]-W[2])+abs(N[2]-NW[2])+abs(NE[2]-N[2])+abs(W[2]-NW[2])+((abs(N[2]-NN[2])+abs(W[2]-WW[2]))>>1)+1,
				abs(N[0]-W[0])+abs(N[0]-NW[0])+abs(NE[0]-N[0])+abs(W[0]-NW[0])+((abs(N[0]-NN[0])+abs(W[0]-WW[0]))>>1)+errors2[0]+eptr[0+9]+1,
				abs(N[1]-W[1])+abs(N[1]-NW[1])+abs(NE[1]-N[1])+abs(W[1]-NW[1])+((abs(N[1]-NN[1])+abs(W[1]-WW[1]))>>1)+errors2[1]+eptr[1+9]+1,
				abs(N[2]-W[2])+abs(N[2]-NW[2])+abs(NE[2]-N[2])+abs(W[2]-NW[2])+((abs(N[2]-NN[2])+abs(W[2]-WW[2]))>>1)+errors2[2]+eptr[2+9]+1,

				//abs(N[0]-W[0])+abs(N[0]-NW[0])+abs(NE[0]-N[0])+abs(W[0]-NW[0])+((abs(N[0]-NN[0])+abs(W[0]-WW[0]))>>1)+2*(errors[0]+errors2[0])+errors3[0]+1,
				//abs(N[1]-W[1])+abs(N[1]-NW[1])+abs(NE[1]-N[1])+abs(W[1]-NW[1])+((abs(N[1]-NN[1])+abs(W[1]-WW[1]))>>1)+2*(errors[1]+errors2[1])+errors3[1]+1,
				//abs(N[2]-W[2])+abs(N[2]-NW[2])+abs(NE[2]-N[2])+abs(W[2]-NW[2])+((abs(N[2]-NN[2])+abs(W[2]-WW[2]))>>1)+2*(errors[2]+errors2[2])+errors3[2]+1,
				//abs(N[0]-W[0])+abs(N[0]-NW[0])+abs(NE[0]-N[0])+abs(W[0]-NW[0])+2*(abs(errors[0])+abs(errors2[0]))+abs(errors3[0]),
				//abs(N[1]-W[1])+abs(N[1]-NW[1])+abs(NE[1]-N[1])+abs(W[1]-NW[1])+2*(abs(errors[1])+abs(errors2[1]))+abs(errors3[1]),
				//abs(N[2]-W[2])+abs(N[2]-NW[2])+abs(NE[2]-N[2])+abs(W[2]-NW[2])+2*(abs(errors[2])+abs(errors2[2]))+abs(errors3[2]),
				//(4*abs(errors[0])+3*abs(errors2[0])+abs(errors3[0]))>>3&127,
				//(4*abs(errors[1])+3*abs(errors2[1])+abs(errors3[1]))>>3&127,
				//(4*abs(errors[2])+3*abs(errors2[2])+abs(errors3[2]))>>3&127,
				//(3*(abs(errors[0])+abs(errors2[0]))+2*abs(errors3[0]))>>3&127,
				//(3*(abs(errors[1])+abs(errors2[1]))+2*abs(errors3[1]))>>3&127,
				//(3*(abs(errors[2])+abs(errors2[2]))+2*abs(errors3[2]))>>3&127,
			};
			eptr[0]=errors2[0];
			eptr[1]=errors2[1];
			eptr[2]=errors2[2];
			eptr+=3;
			//errors3[0]=errors2[0];
			//errors3[1]=errors2[1];
			//errors3[2]=errors2[2];
			//errors2[0]=errors[0];
			//errors2[1]=errors[1];
			//errors2[2]=errors[2];
#ifdef ENABLE_SSE
			int *sse_ptrs[]=
			{
				&sse[0][(N[0]-3*(W[0]-WW[0])-WWW[0])>>4&63][preds[0]>>2&63],//best
				&sse[1][(N[1]-3*(W[1]-WW[1])-WWW[1])>>4&63][preds[1]>>2&63],
				&sse[2][(N[2]-3*(W[2]-WW[2])-WWW[2])>>4&63][preds[2]>>2&63],
				//&sse[0][(N[0]-3*(W[0]-WW[0])-WWW[0])>>4&63][(errors2[0]+eptr[0+9])>>5&63],
				//&sse[1][(N[1]-3*(W[1]-WW[1])-WWW[1])>>4&63][(errors2[1]+eptr[1+9])>>5&63],
				//&sse[2][(N[2]-3*(W[2]-WW[2])-WWW[2])>>4&63][(errors2[2]+eptr[2+9])>>5&63],
				//&sse[0][(errors2[0]+eptr[0+9])>>6&63][preds[0]>>2&63],
				//&sse[1][(errors2[1]+eptr[1+9])>>6&63][preds[1]>>2&63],
				//&sse[2][(errors2[2]+eptr[2+9])>>6&63][preds[2]>>2&63],
				//&sse[0][ctx[0]*3>>8&63][preds[0]>>2&63],
				//&sse[1][ctx[1]*3>>8&63][preds[1]>>2&63],
				//&sse[2][ctx[2]*3>>8&63][preds[2]>>2&63],
				//&sse[0][(N[0]+NE[0]-6*(W[0]-WW[0])-2*WWW[0])>>4&31][preds[0]>>2&63],
				//&sse[1][(N[1]+NE[1]-6*(W[1]-WW[1])-2*WWW[1])>>4&31][preds[1]>>2&63],
				//&sse[2][(N[2]+NE[2]-6*(W[2]-WW[2])-2*WWW[2])>>4&31][preds[2]>>2&63],
				//&sse[0][(N[0]-2*W[0]+WW[0])>>3&31][preds[0]>>2&63],
				//&sse[1][(N[1]-2*W[1]+WW[1])>>3&31][preds[1]>>2&63],
				//&sse[2][(N[2]-2*W[2]+WW[2])>>3&31][preds[2]>>2&63],
				//&sse[0][(N[0]-W[0])>>3&31][preds[0]>>2&63],
				//&sse[1][(N[1]-W[1])>>3&31][preds[1]>>2&63],
				//&sse[2][(N[2]-W[2])>>3&31][preds[2]>>2&63],
				//&sse[0][N[0]>>2&63][W[0]>>2&63],
				//&sse[1][N[1]>>2&63][W[1]>>2&63],
				//&sse[2][N[2]>>2&63][W[2]>>2&63],

				//&sse[0][preds[0]&255],
				//&sse[1][preds[1]&255],
				//&sse[2][preds[2]&255],
			};
			preds[0]+=(*sse_ptrs[0]+(1<<SSE_FBITS>>1))>>SSE_FBITS;
			CLAMP2(preds[0], vmin[0], vmax[0]);
			preds[1]+=(*sse_ptrs[1]+(1<<SSE_FBITS>>1))>>SSE_FBITS;
			CLAMP2(preds[1], vmin[1], vmax[1]);
			preds[2]+=(*sse_ptrs[2]+(1<<SSE_FBITS>>1))>>SSE_FBITS;
			CLAMP2(preds[2], vmin[2], vmax[2]);
#else
			CLAMP2(preds[0], vmin[0], vmax[0]);
			CLAMP2(preds[1], vmin[1], vmax[1]);
			CLAMP2(preds[2], vmin[2], vmax[2]);
			(void)WWW;
#endif
			ctx[0]=FLOOR_LOG2(ctx[0]);
			ctx[1]=FLOOR_LOG2(ctx[1]);
			ctx[2]=FLOOR_LOG2(ctx[2]);
			curr_hist[0]=hist[0][ctx[0]];
			curr_hist[1]=hist[1][ctx[1]];
			curr_hist[2]=hist[2][ctx[2]];
#ifdef AC_VALIDATE
			unsigned long long lo0, r0;
#endif
#ifdef ENABLE_VAL
			unsigned long long val_low, val_range, val_code;
#endif
			if(fwd)
			{
				yuv[0]=ptr[yidx]-128;
				yuv[1]=ptr[uidx]-128;
				yuv[2]=ptr[vidx]-128;
				preds[1]+=yuv[uhelpidx]; CLAMP2(preds[1], -128, 127);
				preds[2]+=yuv[vhelpidx]; CLAMP2(preds[2], -128, 127);
				errors[0]=(char)(yuv[0]-preds[0]);
				errors[1]=(char)(yuv[1]-preds[1]);
				errors[2]=(char)(yuv[2]-preds[2]);
				deltas[0]=errors[0]<<1^errors[0]>>31;
				deltas[1]=errors[1]<<1^errors[1]>>31;
				deltas[2]=errors[2]<<1^errors[2]>>31;

#ifdef USE_MIXIN
				cdf=curr_hist[0][deltas[0]];
				freq=(deltas[0]>=255?0x10000-256:curr_hist[0][deltas[0]+1])-cdf+1;
				cdf+=deltas[0];

#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
#endif
				if(range<0x10000)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_enc(deltas[0], cdf, freq, lo0, lo0+r0, low, low+range, 0, 0);
				val_enc(cdf, freq, val_low, val_range, low, range);

				cdf=curr_hist[1][deltas[1]];
				freq=(deltas[1]>=255?0x10000-256:curr_hist[1][deltas[1]+1])-cdf+1;
				cdf+=deltas[1];

#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
#endif
				if(range<0x10000)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_enc(deltas[1], cdf, freq, lo0, lo0+r0, low, low+range, 0, 0);
				val_enc(cdf, freq, val_low, val_range, low, range);
				
				cdf=curr_hist[2][deltas[2]];
				freq=(deltas[2]>=255?0x10000-256:curr_hist[2][deltas[2]+1])-cdf+1;
				cdf+=deltas[2];

#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
#endif
				if(range<0x10000)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_enc(deltas[2], cdf, freq, lo0, lo0+r0, low, low+range, 0, 0);
				val_enc(cdf, freq, val_low, val_range, low, range);
#else
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[0][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[0][256]+256;
#endif
				cdf=deltas[0];
				freq=curr_hist[0][cdf]+1;
				for(int k=0;k<deltas[0];++k)
					cdf+=curr_hist[0][k];
				if(range<den)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
#ifdef ENABLE_RECIPROCAL
				//unsigned long long low0=low, range0=range;
				//unsigned long long low2=low, range2=range;
				//low2+=(range2*cdf>>8)/(den>>8);
				//range2=(range2*freq>>8)/(den>>8)-1;
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				//MULHI64(hi0, range*cdf>>8, invden);	//these 5 lines take 50% of encode time
				//MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
				//if((low^low2)||(range^range2))
				//	LOG_ERROR("");
#else
				low+=range*cdf/den;//these DIVs take 31.55% E, 31.46% D
				range=range*freq/den-1;
#endif
				
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[1][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[1][256]+256;
#endif
				cdf=deltas[1];
				freq=curr_hist[1][cdf]+1;
				for(int k=0;k<deltas[1];++k)
					cdf+=curr_hist[1][k];
				if(range<den)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
#ifdef ENABLE_RECIPROCAL
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				MULHI64(hi0, range*cdf>>8, invden);
				MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
#else
				low+=range*cdf/den;
				range=range*freq/den-1;
#endif
				
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[2][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[2][256]+256;
#endif
				cdf=deltas[2];
				freq=curr_hist[2][cdf]+1;
				for(int k=0;k<deltas[2];++k)
					cdf+=curr_hist[2][k];
				if(range<den)
				{
					*(unsigned*)streamptr=(unsigned)(low>>32);
					streamptr+=4;
#ifdef _DEBUG
					if(streamptr>NNptr)
						LOG_ERROR("");
#endif
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
#ifdef ENABLE_RECIPROCAL
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				MULHI64(hi0, range*cdf>>8, invden);
				MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
#else
				low+=range*cdf/den;
				range=range*freq/den-1;
#endif
#endif
			}
			else
			{
				int c;
#ifdef USE_MIXIN
#ifdef _DEBUG
				if(low>code||low+range<code)
					LOG_ERROR("XY %d %d", kx, ky);
#endif
#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
				val_code=code;
#endif
				if(range<0x10000)
				{
#ifdef _DEBUG
					if(ptr>=streamptr)
						LOG_ERROR("");
#endif
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)<<16|0xFFFF)/range);
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=sym>=255?0x10000:curr_hist[0][sym+1]+sym+1;
					//if(curr_hist[0][sym]+sym>freq)
					//	LOG_ERROR("Invalid CDF %d, %d", curr_hist[0][sym]+sym, freq);
					if(freq>c)
						break;
#ifdef _DEBUG
					if(sym>=256)
						LOG_ERROR("0XY %d %d sym OOB %d CDF %d C %d", kx, ky, sym, freq, c);
#endif
					++sym;
					cdf=freq;
				}
				freq-=cdf;
				if(freq<0||cdf<0||cdf>=0x10000)
					LOG_ERROR("");
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_dec(sym, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
				VAL_DEC(cdf, freq, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
				
				deltas[0]=sym;
				errors[0]=deltas[0]>>1^-(deltas[0]&1);
				yuv[0]=errors[0]+preds[0];
				preds[1]+=yuv[uhelpidx]; CLAMP2(preds[1], -128, 127);

#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
				val_code=code;
#endif
				if(range<0x10000)
				{
#ifdef _DEBUG
					if(ptr>=streamptr)
						LOG_ERROR("");
#endif
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)<<16|0xFFFF)/range);
#ifdef _DEBUG
				if(low>code||low+range<code)
					LOG_ERROR("");
#endif
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=sym>=255?0x10000:curr_hist[1][sym+1]+sym+1;
					if(freq>c)
						break;
#ifdef _DEBUG
					if(sym>=256)
						LOG_ERROR("1XY %d %d sym OOB %d", kx, ky, sym);
#endif
					++sym;
					cdf=freq;
				}
				freq-=cdf;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_dec(sym, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
				VAL_DEC(cdf, freq, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
				
				deltas[1]=sym;
				errors[1]=deltas[1]>>1^-(deltas[1]&1);
				yuv[1]=errors[1]+preds[1];
				preds[2]+=yuv[vhelpidx]; CLAMP2(preds[2], -128, 127);
				
#ifdef AC_VALIDATE
				lo0=low, r0=range;
#endif
#ifdef ENABLE_VAL
				val_low=low;
				val_range=range;
				val_code=code;
#endif
				if(range<0x10000)
				{
#ifdef _DEBUG
					if(ptr>=streamptr)
						LOG_ERROR("");
#endif
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)<<16|0xFFFF)/range);
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=sym>=255?0x10000:curr_hist[2][sym+1]+sym+1;
					if(freq>c)
						break;
#ifdef _DEBUG
					if(sym>=256)
						LOG_ERROR("2XY %d %d sym OOB %d", kx, ky, sym);
#endif
					++sym;
					cdf=freq;
				}
				freq-=cdf;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
				acval_dec(sym, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
				VAL_DEC(cdf, freq, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);

				deltas[2]=sym;
				errors[2]=deltas[2]>>1^-(deltas[2]&1);
				yuv[2]=errors[2]+preds[2];
#else
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[0][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[0][256]+256;
#endif
				if(range<den)
				{
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)*den+den-1)/range);
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=curr_hist[0][sym]+1;
					int cdf2=cdf+freq;
					if(cdf2>c)
						break;
					++sym;
					cdf=cdf2;
				}
#ifdef ENABLE_RECIPROCAL
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				MULHI64(hi0, range*cdf>>8, invden);
				MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
#else
				low+=range*cdf/den;
				range=range*freq/den-1;
#endif
				deltas[0]=sym;
				errors[0]=deltas[0]>>1^-(deltas[0]&1);
				yuv[0]=errors[0]+preds[0];
				preds[1]+=yuv[uhelpidx]; CLAMP2(preds[1], -128, 127);
				
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[1][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[1][256]+256;
#endif
				if(range<den)
				{
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)*den+den-1)/range);
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=curr_hist[1][sym]+1;
					int cdf2=cdf+freq;
					if(cdf2>c)
						break;
					++sym;
					cdf=cdf2;
				}
#ifdef ENABLE_RECIPROCAL
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				MULHI64(hi0, range*cdf>>8, invden);
				MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
#else
				low+=range*cdf/den;
				range=range*freq/den-1;
#endif
				deltas[1]=sym;
				errors[1]=deltas[1]>>1^-(deltas[1]&1);
				yuv[1]=errors[1]+preds[1];
				preds[2]+=yuv[vhelpidx]; CLAMP2(preds[2], -128, 127);
				
#ifdef ENABLE_RECIPROCAL
				den=(curr_hist[2][256]+(256+(1<<8)-1))&~255;
#else
				den=curr_hist[2][256]+256;
#endif
				if(range<den)
				{
					code=code<<32|*(unsigned*)streamptr;
					streamptr+=4;
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					if(range>~low)
						range=~low;
				}
				c=(int)(((code-low)*den+den-1)/range);
				cdf=0;
				sym=0;
				for(;;)
				{
					freq=curr_hist[2][sym]+1;
					int cdf2=cdf+freq;
					if(cdf2>c)
						break;
					++sym;
					cdf=cdf2;
				}
#ifdef ENABLE_RECIPROCAL
				invden=divtable[den>>8];
				_mulx_u64(range*cdf>>8, invden, &hi0);
				_mulx_u64(range*freq>>8, invden, &hi1);
				MULHI64(hi0, range*cdf>>8, invden);
				MULHI64(hi1, range*freq>>8, invden);
				low+=hi0;
				range=hi1-1;
#else
				low+=range*cdf/den;
				range=range*freq/den-1;
#endif
				deltas[2]=sym;
				errors[2]=deltas[2]>>1^-(deltas[2]&1);
				yuv[2]=errors[2]+preds[2];
#endif
				ptr[yidx]=yuv[0]+128;
				ptr[uidx]=yuv[1]+128;
				ptr[vidx]=yuv[2]+128;

				guide_check(image, kx, ky);
			}
#ifdef _DEBUG
			if(deltas[0]>>8||deltas[1]>>8||deltas[2]>>8)
				LOG_ERROR("");
#endif
#ifdef USE_MIXIN
#if 0
			for(int ks=1;ks<256;++ks)//56.06%
			{
				//curr_hist[0][ks]+=(((0x10000-256)&-(ks>deltas[0]))-curr_hist[0][ks]+(1<<RATE>>1))>>RATE;
				//curr_hist[1][ks]+=(((0x10000-256)&-(ks>deltas[1]))-curr_hist[1][ks]+(1<<RATE>>1))>>RATE;
				//curr_hist[2][ks]+=(((0x10000-256)&-(ks>deltas[2]))-curr_hist[2][ks]+(1<<RATE>>1))>>RATE;
			
				int mc0=curr_hist[0][ks];
				int mc1=curr_hist[1][ks];
				int mc2=curr_hist[2][ks];
				mc0-=mc0>>RATE;
				mc1-=mc1>>RATE;
				mc2-=mc2>>RATE;
				mc0+=(0x10000-256)>>RATE&-(ks>deltas[0]);
				mc1+=(0x10000-256)>>RATE&-(ks>deltas[1]);
				mc2+=(0x10000-256)>>RATE&-(ks>deltas[2]);
				curr_hist[0][ks]=mc0;
				curr_hist[1][ks]=mc1;
				curr_hist[2][ks]=mc2;
			}
#else
			__m256i t0, t1, t2, t3;
			__m256i mc0, mc1, mc2, mc3;
			unsigned short *curr_mixin;

			curr_mixin=mixin+255-deltas[0];
			mc0=_mm256_loadu_si256((__m256i*)curr_hist[0]+0);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[0]+1);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[0]+2);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[0]+3);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+0));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+1));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+2));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+3));
			_mm256_storeu_si256((__m256i*)curr_hist[0]+0, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+1, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+2, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+3, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[0]+4);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[0]+5);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[0]+6);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[0]+7);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+4));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+5));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+6));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+7));
			_mm256_storeu_si256((__m256i*)curr_hist[0]+4, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+5, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+6, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+7, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[0]+ 8);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[0]+ 9);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[0]+10);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[0]+11);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+ 8));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+ 9));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+10));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+11));
			_mm256_storeu_si256((__m256i*)curr_hist[0]+ 8, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+ 9, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+10, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+11, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[0]+12);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[0]+13);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[0]+14);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[0]+15);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+12));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+13));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+14));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+15));
			_mm256_storeu_si256((__m256i*)curr_hist[0]+12, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+13, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+14, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[0]+15, mc3);

			curr_mixin=mixin+255-deltas[1];
			mc0=_mm256_loadu_si256((__m256i*)curr_hist[1]+0);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[1]+1);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[1]+2);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[1]+3);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+0));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+1));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+2));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+3));
			_mm256_storeu_si256((__m256i*)curr_hist[1]+0, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+1, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+2, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+3, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[1]+4);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[1]+5);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[1]+6);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[1]+7);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+4));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+5));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+6));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+7));
			_mm256_storeu_si256((__m256i*)curr_hist[1]+4, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+5, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+6, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+7, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[1]+ 8);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[1]+ 9);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[1]+10);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[1]+11);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+ 8));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+ 9));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+10));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+11));
			_mm256_storeu_si256((__m256i*)curr_hist[1]+ 8, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+ 9, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+10, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+11, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[1]+12);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[1]+13);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[1]+14);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[1]+15);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+12));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+13));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+14));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+15));
			_mm256_storeu_si256((__m256i*)curr_hist[1]+12, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+13, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+14, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[1]+15, mc3);

			curr_mixin=mixin+255-deltas[2];
			mc0=_mm256_loadu_si256((__m256i*)curr_hist[2]+0);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[2]+1);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[2]+2);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[2]+3);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+0));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+1));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+2));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+3));
			_mm256_storeu_si256((__m256i*)curr_hist[2]+0, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+1, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+2, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+3, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[2]+4);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[2]+5);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[2]+6);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[2]+7);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+4));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+5));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+6));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+7));
			_mm256_storeu_si256((__m256i*)curr_hist[2]+4, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+5, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+6, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+7, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[2]+ 8);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[2]+ 9);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[2]+10);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[2]+11);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+ 8));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+ 9));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+10));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+11));
			_mm256_storeu_si256((__m256i*)curr_hist[2]+ 8, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+ 9, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+10, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+11, mc3);

			mc0=_mm256_loadu_si256((__m256i*)curr_hist[2]+12);
			mc1=_mm256_loadu_si256((__m256i*)curr_hist[2]+13);
			mc2=_mm256_loadu_si256((__m256i*)curr_hist[2]+14);
			mc3=_mm256_loadu_si256((__m256i*)curr_hist[2]+15);
			t0=_mm256_srli_epi16(mc0, RATE);
			t1=_mm256_srli_epi16(mc1, RATE);
			t2=_mm256_srli_epi16(mc2, RATE);
			t3=_mm256_srli_epi16(mc3, RATE);
			mc0=_mm256_sub_epi16(mc0, t0);
			mc1=_mm256_sub_epi16(mc1, t1);
			mc2=_mm256_sub_epi16(mc2, t2);
			mc3=_mm256_sub_epi16(mc3, t3);
			mc0=_mm256_add_epi16(mc0, _mm256_loadu_si256((__m256i*)curr_mixin+12));
			mc1=_mm256_add_epi16(mc1, _mm256_loadu_si256((__m256i*)curr_mixin+13));
			mc2=_mm256_add_epi16(mc2, _mm256_loadu_si256((__m256i*)curr_mixin+14));
			mc3=_mm256_add_epi16(mc3, _mm256_loadu_si256((__m256i*)curr_mixin+15));
			_mm256_storeu_si256((__m256i*)curr_hist[2]+12, mc0);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+13, mc1);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+14, mc2);
			_mm256_storeu_si256((__m256i*)curr_hist[2]+15, mc3);
#endif
#else
			++curr_hist[0][deltas[0]];
			++curr_hist[1][deltas[1]];
			++curr_hist[2][deltas[2]];
			++curr_hist[0][256];
			++curr_hist[1][256];
			++curr_hist[2][256];
			if(curr_hist[0][256]>=0x10000-256-256)
			{
				int sum=0;
				for(int ks=0;ks<256;++ks)
					sum+=curr_hist[0][ks]>>=1;
				//	sum+=hist[0][ks]=(hist[0][ks]+1)>>1;
				curr_hist[0][256]=sum;
			}
			if(curr_hist[1][256]>=0x10000-256-256)
			{
				int sum=0;
				for(int ks=0;ks<256;++ks)
					sum+=curr_hist[1][ks]>>=1;
				//	sum+=hist[1][ks]=(hist[1][ks]+1)>>1;
				curr_hist[1][256]=sum;
			}
			if(curr_hist[2][256]>=0x10000-256-256)
			{
				int sum=0;
				for(int ks=0;ks<256;++ks)
					sum+=curr_hist[2][ks]>>=1;
				//	sum+=hist[2][ks]=(hist[2][ks]+1)>>1;
				curr_hist[2][256]=sum;
			}
#endif
			yuv[2]-=yuv[vhelpidx];
			yuv[1]-=yuv[uhelpidx];
			//errors[0]=yuv[0]-preds[0];
			//errors[1]=yuv[1]-preds[1];
			//errors[2]=yuv[2]-preds[2];
#ifdef ENABLE_SSE
			*sse_ptrs[0]+=((errors[0]<<SSE_FBITS)-*sse_ptrs[0])>>SSE_LBITS;
			*sse_ptrs[1]+=((errors[1]<<SSE_FBITS)-*sse_ptrs[1])>>SSE_LBITS;
			*sse_ptrs[2]+=((errors[2]<<SSE_FBITS)-*sse_ptrs[2])>>SSE_LBITS;
#endif
			NW[0]=N[0];
			NW[1]=N[1];
			NW[2]=N[2];
			N[0]=NE[0];
			N[1]=NE[1];
			N[2]=NE[2];
			WWW[0]=WW[0];
			WWW[1]=WW[1];
			WWW[2]=WW[2];
			WW[0]=W[0];
			WW[1]=W[1];
			WW[2]=W[2];
			W[0]=yuv[0];
			W[1]=yuv[1];
			W[2]=yuv[2];
			//*(unsigned*)NW=*(unsigned*)N;//SLOW
			//*(unsigned*)W=*(unsigned*)yuv;

			NNptr+=3;
			Nptr+=3;
			ptr+=3;
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
			*(unsigned*)streamptr=(unsigned)(low>>32);//flush
			streamptr+=4;
			*(unsigned*)streamptr=(unsigned)low;
			streamptr+=4;

			fwrite("10\n", 1, 3, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(stream, 1, streamptr-stream, fdst);
#ifdef _DEBUG
			printf("C %td bytes\n", streamptr-stream+headersize);
#endif
		}
		else
		{
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(buf);
	return 0;
}