#include"fpc/fpc.h"
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


	#define LOUD
	#define PROFILE_TIME

//	#define USE_W
	#define USE_CG
//	#define USE_WGCG

//	#define DISABLE_ANALYSIS


#define ENC_DY 128
#define ENC_DX 512
#define YPAD 2

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
#ifdef PROFILE_TIME
#define PROFLIST\
	PROFLABEL(enc_read)\
	PROFLABEL(enc_interleave)\
	PROFLABEL(enc_analysis)\
	PROFLABEL(enc_predict)\
	PROFLABEL(enc_FPC)\
	PROFLABEL(enc_write)\
	PROFLABEL(enc_count)\
	PROFLABEL(dec_read)\
	PROFLABEL(dec_FPC)\
	PROFLABEL(dec_predict)\
	PROFLABEL(dec_interleave)\
	PROFLABEL(dec_write)
typedef enum _ProfLabel
{
#define PROFLABEL(LABEL) PROF_##LABEL,
	PROFLIST
#undef  PROFLABEL
	PROF_COUNT,
} ProfLabel;
static const char *profnames[]=
{
#define PROFLABEL(LABEL) #LABEL,
	PROFLIST
#undef  PROFLABEL
};
static double g_profinfo[128]={0};
static ptrdiff_t g_volume[128]={0};
static double g_t=0;
static int g_profstart=0, g_profend=0;
static void prof_start()
{
	g_t=time_sec();
	memset(g_profinfo, 0, sizeof(g_profinfo));
	memset(g_volume, 0, sizeof(g_volume));
}
static void prof_checkpoint(int idx, ptrdiff_t volume)
{
	if(idx>=_countof(g_profinfo))
	{
		CRASH("Profiler OOB");
		return;
	}
	{
		double t2=time_sec();
		g_profinfo[idx]+=t2-g_t;
		g_volume[idx]+=volume;
		g_t=t2;
		if(g_profend<idx+1)
			g_profend=idx+1;
	}
}
#define PROF(LABEL, VOLUME) prof_checkpoint(PROF_##LABEL, VOLUME)
static void prof_print(void)
{
	double sum=0;
	for(int k=g_profstart;k<g_profend;++k)
		sum+=g_profinfo[k];
	const int scale=4;//ms
	printf("1 char = %d ms\n", scale);
	printf("|");
	double csum=0;
	int prev=0;
	for(int k=g_profstart;k<g_profend;++k)
	{
		double val=g_profinfo[k];
		csum+=val;
		int next=(int)(csum*1000/scale);
		int nstars=next-prev;
		prev=next;
		for(int k2=0;k2<nstars/2;++k2)
			printf("-");
		printf("%d", k-g_profstart+1);
		for(int k2=nstars/2;k2<nstars;++k2)
			printf("-");
		printf("|");
	}
	printf("\n");
	ptrdiff_t smax=0;
	for(int k=g_profstart;k<g_profend;++k)
	{
		double elapsed=g_profinfo[k];
		printf("%8.4lf%%  %12.6lf sec  %10td bytes  %12.6lf MB/s  %2d  %s\n",
			100.*elapsed/sum,
			elapsed,
			g_volume[k],
			g_volume[k]/(elapsed*1024*1024),
			k-g_profstart+1,
			profnames[k]
		);
		if(smax<g_volume[k])
			smax=g_volume[k];
	}
	printf("\n");
	printf("%lf sec  %12.6lf MB/s\n"
		, sum
		, smax/(sum*1024*1024)
	);
	printf("\n");
}
#else
#define prof_start()
#define PROF(...)
#define prof_print()
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
static char im2[3][(ENC_DY+YPAD)*ENC_DX];
static char cim[3][(ENC_DY+YPAD)*ENC_DX];
#ifndef USE_CG
#define DIVTABLESIZE 2048
static int divtable[DIVTABLESIZE];
#endif
int a04_codec(int argc, char **argv)
{
	prof_start();
	if(argc<3)
	{
		printf("Usage:  \"%s\"  input  output\n", argv[0]);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];

	FILE *fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"\n", srcfn);
		return 0;
	}
	int iw=0, ih=0, fwd=0;
	{
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'4'<<8))
		{
			CRASH("Unsupported file  tag 0x%04X", tag);
			return 1;
		}
	}
	if(fwd)
	{
		char c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported file");
			return 1;
		}
		c=fgetc(fsrc);
		while(c=='#')//skip comments
		{
			c=fgetc(fsrc);
			while(c!='\n')c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		iw=0;
		while((unsigned)(c-'0')<10)
		{
			iw=10*iw+c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')c=fgetc(fsrc);
		ih=0;
		while((unsigned)(c-'0')<10)
		{
			ih=10*ih+c-'0';
			c=fgetc(fsrc);
		}
		while(c!='\n')c=getc(fsrc);
		while(c=='#')//skip comments
		{
			c=fgetc(fsrc);
			while(c!='\n')c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		int tag=0;
		fread(&tag, 1, 4, fsrc);
		if(memcmp(&tag, "255\n", 4))
		{
			CRASH("Unsupported file");
			return 1;
		}
#ifdef LOUD
		printf("\"%s\" %d*%d\n", srcfn, iw, ih);//
#endif
	}
	else
	{
		fread(&iw, 1, 4, fsrc);
		fread(&ih, 1, 4, fsrc);
#ifdef PROFILE_TIME
		g_profstart=PROF_enc_count+1;
#endif
	}
	//int maxsum=0;
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fdst);
		return 1;
	}
	if(fwd)
	{
		fwrite("04", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
	}
	else
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	int imsize=iw*ENC_DY*3;
	unsigned char *image=(unsigned char*)malloc(imsize);
	if(!image)
	{
		CRASH("Alloc error");
		return 1;
	}
#ifndef USE_CG
	for(int k=1;k<DIVTABLESIZE;++k)
		divtable[k]=(0x1000000+k)/(k<<1);
#endif
	for(int y1=0;y1<ih;y1+=ENC_DY)
	{
		int y2=y1+ENC_DY;
		if(y2>ih)
			y2=ih;
		int dy=y2-y1;
		if(fwd)
		{
			int req=iw*dy*3;
			ptrdiff_t nread=fread(image, 1, req, fsrc);
			if(nread!=req)
			{
				CRASH("Truncated image");
				return 1;
			}
			PROF(enc_read, req);
			for(int x1=0;x1<iw;x1+=ENC_DX)//enc block
			{
				int x2=x1+ENC_DX;
				if(x2>iw)
					x2=iw;
				int dx=x2-x1;
				int blocksize=dx*dy, blockbytes=3*blocksize;
				memset(im2, 0, sizeof(im2));

				//deinterleave
				{
					char
						*ydstptr=im2[0]+YPAD*ENC_DX,
						*udstptr=im2[1]+YPAD*ENC_DX,
						*vdstptr=im2[2]+YPAD*ENC_DX;
					for(int ky=0;ky<dy;++ky)
					{
						const unsigned char *row=image+3*(iw*ky+x1);
						for(int kx=0;kx<dx;++kx)
						{
							*ydstptr++=*row++-128;
							*udstptr++=*row++-128;
							*vdstptr++=*row++-128;
						}
					}
				}
				PROF(enc_interleave, blockbytes);

				//analysis
				int bestrct=8;
#ifndef DISABLE_ANALYSIS
				//int Wcoeff=0;
				{
					//const char
					//	*yNptr=im2[0]+YPAD*ENC_DX-dx,
					//	*uNptr=im2[1]+YPAD*ENC_DX-dx,
					//	*vNptr=im2[2]+YPAD*ENC_DX-dx;
					const char
						*yptr=im2[0]+YPAD*ENC_DX,
						*uptr=im2[1]+YPAD*ENC_DX,
						*vptr=im2[2]+YPAD*ENC_DX;
					long long counters[OCH_COUNT]={0};
					//long long counters2[OCH_COUNT]={0};
					for(int ky=0, idx=0;ky<dy;++ky)
					{
						//int preN[OCH_COUNT]={0};
						int prev[OCH_COUNT]={0};
						for(int kx=0;kx<dx;++kx, ++idx)
						{
							int r=yptr[idx];
							int g=uptr[idx];
							int b=vptr[idx];
							int rg=r-g;
							int gb=g-b;
							int br=b-r;
							//int rN=yptr[idx-dx];
							//int gN=uptr[idx-dx];
							//int bN=vptr[idx-dx];
							//int rgN=rN-gN;
							//int gbN=gN-bN;
							//int brN=bN-rN;
							//counters[0]+=abs(r	-(prev[0]+rN)/2);
							//counters[1]+=abs(g	-(prev[1]+gN)/2);
							//counters[2]+=abs(b	-(prev[2]+bN)/2);
							//counters[3]+=abs(rg	-(prev[3]+rgN)/2);
							//counters[4]+=abs(gb	-(prev[4]+gbN)/2);
							//counters[5]+=abs(br	-(prev[5]+brN)/2);
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
							//counters2[0]+=abs(r	-rN);
							//counters2[1]+=abs(g	-gN);
							//counters2[2]+=abs(b	-bN);
							//counters2[3]+=abs(rg	-rgN);
							//counters2[4]+=abs(gb	-gbN);
							//counters2[5]+=abs(br	-brN);
							//preN[0]=rN;
							//preN[1]=gN;
							//preN[2]=bN;
							//preN[3]=rgN;
							//preN[4]=gbN;
							//preN[5]=brN;
						}
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
					//int bestrct2=0;
					//long long minerr2=0;
					//for(int kt=0;kt<_countof(rct_indices);++kt)
					//{
					//	const unsigned char *rct=rct_indices[kt];
					//	long long currerr=
					//		+counters2[rct[0]]
					//		+counters2[rct[1]]
					//		+counters2[rct[2]]
					//	;
					//	if(!kt||minerr>currerr)
					//	{
					//		minerr=currerr;
					//		bestrct2=kt;
					//	}
					//}
					//if(minerr>minerr2)
					//	bestrct=bestrct2;
					//long long minerrsum=minerr+minerr2;
					//Wcoeff=(int)(((minerr<<8)+(minerrsum>>1))/minerrsum);
					//Wcoeff=(Wcoeff-128)/32+128;
					//CLAMP2(Wcoeff, 1, 255);
				}
#endif
				PROF(enc_analysis, blockbytes);
				
				//decorrelation
				{
					int
						yidx=rct_indices[bestrct][3+0],
						uidx=rct_indices[bestrct][3+1],
						vidx=rct_indices[bestrct][3+2],
						uhelpidx=rct_indices[bestrct][6+0],
						vhelpidx=rct_indices[bestrct][6+1];
#ifndef USE_W
#ifndef USE_CG
					const char
						*yNNptr=im2[0]+YPAD*ENC_DX-2*dx,
						*uNNptr=im2[1]+YPAD*ENC_DX-2*dx,
						*vNNptr=im2[2]+YPAD*ENC_DX-2*dx;
#endif
					const char
						*yNptr=im2[0]+YPAD*ENC_DX-dx,
						*uNptr=im2[1]+YPAD*ENC_DX-dx,
						*vNptr=im2[2]+YPAD*ENC_DX-dx;
#endif
					char
						*ycurrptr=im2[0]+YPAD*ENC_DX,//YUV
						*ucurrptr=im2[1]+YPAD*ENC_DX,
						*vcurrptr=im2[2]+YPAD*ENC_DX;
					const char
						*ysrcptr=im2[yidx]+YPAD*ENC_DX,//RGB -> YUV
						*usrcptr=im2[uidx]+YPAD*ENC_DX,
						*vsrcptr=im2[vidx]+YPAD*ENC_DX;
					char
						*ydstptr=im2[0],//residuals
						*udstptr=im2[1],
						*vdstptr=im2[2];
#ifdef USE_W
					int ypred=0;
					int upred=0;
					int vpred=0;
#endif
#ifdef __GNUC__
#pragma GCC unroll 2
#endif
					for(int k=0;k<blocksize;++k)
					{
#ifdef USE_CG
						int yN=yNptr[0];
						int uN=uNptr[0];
						int vN=vNptr[0];
						int yW=ycurrptr[-1];
						int uW=ucurrptr[-1];
						int vW=vcurrptr[-1];
						int ymax=yN, ymin=yW;
						int umax=uN, umin=uW;
						int vmax=vN, vmin=vW;
						if(yN<yW)ymin=yN, ymax=yW;
						if(uN<uW)umin=uN, umax=uW;
						if(vN<vW)vmin=vN, vmax=vW;
					//	int ypred=2*yN+((yW-yN)*Wcoeff>>7)-yNptr[-1];
					//	int upred=2*uN+((uW-uN)*Wcoeff>>7)-uNptr[-1];
					//	int vpred=2*vN+((vW-vN)*Wcoeff>>7)-vNptr[-1];
						int ypred=yN+yW-yNptr[-1];
						int upred=uN+uW-uNptr[-1];
						int vpred=vN+vW-vNptr[-1];
						CLAMP2(ypred, ymin, ymax);
						CLAMP2(upred, umin, umax);
						CLAMP2(vpred, vmin, vmax);
#elif defined USE_WGCG
						int yN=yNptr[0];
						int uN=uNptr[0];
						int vN=vNptr[0];
						int yW=ycurrptr[-1];
						int uW=ucurrptr[-1];
						int vW=vcurrptr[-1];
						int ygx=abs(yW-ycurrptr[-2])+abs(yN-yNptr[-1])+abs(yNptr[1]-yNptr[0])+1;//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
						int ugx=abs(uW-ucurrptr[-2])+abs(uN-uNptr[-1])+abs(uNptr[1]-uNptr[0])+1;
						int vgx=abs(vW-vcurrptr[-2])+abs(vN-vNptr[-1])+abs(vNptr[1]-vNptr[0])+1;
						int ygy=abs(yW-yNptr[-1])+abs(yNptr[0]-yNNptr[0])+abs(yNptr[1]-yNNptr[1])+1;//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
						int ugy=abs(uW-uNptr[-1])+abs(uNptr[0]-uNNptr[0])+abs(uNptr[1]-uNNptr[1])+1;
						int vgy=abs(vW-vNptr[-1])+abs(vNptr[0]-vNNptr[0])+abs(vNptr[1]-vNNptr[1])+1;
						int ysum=ygx+ygy;
						int usum=ugx+ugy;
						int vsum=vgx+vgy;
						int ypred=(ygx+ysum)*yN+(ygy+ysum)*yW-ysum*yNptr[-1];
						int upred=(ugx+usum)*uN+(ugy+usum)*uW-usum*uNptr[-1];
						int vpred=(vgx+vsum)*vN+(vgy+vsum)*vW-vsum*vNptr[-1];
						ypred=(long long)ypred*divtable[ysum]>>24;
						upred=(long long)upred*divtable[usum]>>24;
						vpred=(long long)vpred*divtable[vsum]>>24;
					//	int ypred=((ygx+ysum)*yN+(ygy+ysum)*yW-ysum*yNptr[-1])/(2*ysum);//pred=((gx+sum)*N+(gy+sum)*W-sum*NW)/(2*sum)
					//	int upred=((ugx+usum)*uN+(ugy+usum)*uW-usum*uNptr[-1])/(2*usum);
					//	int vpred=((vgx+vsum)*vN+(vgy+vsum)*vW-vsum*vNptr[-1])/(2*vsum);
						//if(maxsum<ysum)maxsum=ysum;//
						//if(maxsum<usum)maxsum=usum;//
						//if(maxsum<vsum)maxsum=vsum;//

						//int ypred=(ygx*yN+ygy*yW)/ysum;//pred=(gx*N+gy*W)/sum
						//int upred=(ugx*uN+ugy*uW)/usum;
						//int vpred=(vgx*vN+vgy*vW)/vsum;
						int ymax=yN, ymin=yW;
						int umax=uN, umin=uW;
						int vmax=vN, vmin=vW;
						if(yN<yW)ymin=yN, ymax=yW;
						if(uN<uW)umin=uN, umax=uW;
						if(vN<vW)vmin=vN, vmax=vW;
						int yNE=yNptr[1];
						int uNE=uNptr[1];
						int vNE=vNptr[1];
						if(ymin>yNE)ymin=yNE;
						if(umin>uNE)umin=uNE;
						if(vmin>vNE)vmin=vNE;
						if(ymax<yNE)ymax=yNE;
						if(umax<uNE)umax=uNE;
						if(vmax<vNE)vmax=vNE;
						CLAMP2(ypred, ymin, ymax);
						CLAMP2(upred, umin, umax);
						CLAMP2(vpred, vmin, vmax);
						++yNNptr;
						++uNNptr;
						++vNNptr;
#endif
#ifndef USE_W
						++yNptr;
						++uNptr;
						++vNptr;
#endif
						char yuv[4]=
						{
							*ysrcptr++,
							*usrcptr++,
							*vsrcptr++,
							0,
						};
#ifdef DISABLE_ANALYSIS
						yuv[1]-=yuv[0];
						yuv[2]-=yuv[0];
						yuv[0]+=(yuv[1]+yuv[2])>>2;
						yuv[2]-=yuv[1]>>2;
						(void)vhelpidx;
						(void)uhelpidx;
#else
						yuv[2]-=yuv[vhelpidx];//need int16 buffer for cRCT
						yuv[1]-=yuv[uhelpidx];
#endif
#ifdef USE_W
						int yW=yuv[0];
						int uW=yuv[1];
						int vW=yuv[2];
#endif
						*ycurrptr++=yuv[0];//store YUV
						*ucurrptr++=yuv[1];
						*vcurrptr++=yuv[2];
						yuv[0]-=ypred;
						yuv[1]-=upred;
						yuv[2]-=vpred;
						*ydstptr++=yuv[0];//store residuals
						*udstptr++=yuv[1];
						*vdstptr++=yuv[2];
#ifdef USE_W
						ypred=yW;
						upred=uW;
						vpred=vW;
#endif
						//*ydstptr++=(char)(yuv[0]-ypred);
						//upred+=yuv[uhelpidx];
						//CLAMP2(upred, -128, 127);
						//*udstptr++=(char)(yuv[1]-upred);
						//vpred+=yuv[vhelpidx];
						//CLAMP2(vpred, -128, 127);
						//*vdstptr++=(char)(yuv[2]-vpred);
					}
				}
				PROF(enc_predict, blockbytes);

				//write
				int bsize=blocksize;
				if(bsize>0x4000)
					bsize=0x4000;
				size_t csize0=FPC_compress(cim[0], im2[0], blocksize, bsize);
				size_t csize1=FPC_compress(cim[1], im2[1], blocksize, bsize);
				size_t csize2=FPC_compress(cim[2], im2[2], blocksize, bsize);
				if(!csize0||!csize1||!csize2)
				{
					CRASH("FPC_compress %d %d %d", csize0, csize1, csize2);
					return 1;
				}
				PROF(enc_FPC, blockbytes);
				fwrite(&bestrct, 1, 1, fdst);
			//	fwrite(&Wcoeff, 1, 1, fdst);
				fwrite(&csize0, 1, 4, fdst);
				fwrite(&csize1, 1, 4, fdst);
				fwrite(&csize2, 1, 4, fdst);
				fwrite(cim[0], 1, csize0, fdst);
				fwrite(cim[1], 1, csize1, fdst);
				fwrite(cim[2], 1, csize2, fdst);
				PROF(enc_write, csize0+csize1+csize2+1+3*4);
			}
		}
		else
		{
			for(int x1=0;x1<iw;x1+=ENC_DX)//dec block
			{
				int x2=x1+ENC_DX;
				if(x2>iw)
					x2=iw;
				int dx=x2-x1;
				int blocksize=dx*dy, blockbytes=3*blocksize;

				//read
				int bestrct=0, csize0=0, csize1=0, csize2=0;
			//	int Wcoeff=0;
				ptrdiff_t nread=0;
				nread+=fread(&bestrct, 1, 1, fsrc);
			//	nread+=fread(&Wcoeff, 1, 1, fsrc);
				nread+=fread(&csize0, 1, 4, fsrc);
				nread+=fread(&csize1, 1, 4, fsrc);
				nread+=fread(&csize2, 1, 4, fsrc);
				if((int)nread!=1+3*4)
				{
					CRASH("Invalid file");
					return 1;
				}
				nread=0;
				nread+=fread(cim[0], 1, csize0, fsrc);
				nread+=fread(cim[1], 1, csize1, fsrc);
				nread+=fread(cim[2], 1, csize2, fsrc);
				if(nread!=csize0+csize1+csize2)
				{
					CRASH("Invalid file");
					return 1;
				}
				PROF(dec_read, csize0+csize1+csize2+1+3*4);

				//decompress
				size_t usize0=FPC_decompress(im2[0], blocksize, cim[0], csize0);
				size_t usize1=FPC_decompress(im2[1], blocksize, cim[1], csize1);
				size_t usize2=FPC_decompress(im2[2], blocksize, cim[2], csize2);
				if(usize0!=blocksize||usize1!=blocksize||usize2!=blocksize)
				{
					CRASH("FPC_decompress %d %d %d vs %d"
						, usize0
						, usize1
						, usize2
						, blocksize
					);
					return 1;
				}
				PROF(dec_FPC, usize0+usize1+usize2);
				
				memset(cim, 0, sizeof(cim));

				//reconstruct
				{
					int
						yidx=rct_indices[bestrct][3+0],
						uidx=rct_indices[bestrct][3+1],
						vidx=rct_indices[bestrct][3+2],
						uhelpidx=rct_indices[bestrct][6+0],
						vhelpidx=rct_indices[bestrct][6+1];
#ifndef USE_W
#ifndef USE_CG
					const char
						*yNNptr=cim[0]+YPAD*ENC_DX-2*dx,
						*uNNptr=cim[1]+YPAD*ENC_DX-2*dx,
						*vNNptr=cim[2]+YPAD*ENC_DX-2*dx;
#endif
					const char
						*yNptr=cim[0]+YPAD*ENC_DX-dx,
						*uNptr=cim[1]+YPAD*ENC_DX-dx,
						*vNptr=cim[2]+YPAD*ENC_DX-dx;
#endif
					char
						*ycurrptr=cim[0]+YPAD*ENC_DX,//YUV
						*ucurrptr=cim[1]+YPAD*ENC_DX,
						*vcurrptr=cim[2]+YPAD*ENC_DX;
					const char
						*ysrcptr=im2[0],//residuals -> RGB
						*usrcptr=im2[1],
						*vsrcptr=im2[2];
					char
						*ydstptr=im2[yidx],//RGB
						*udstptr=im2[uidx],
						*vdstptr=im2[vidx];
#ifdef USE_W
					int ypred=0;
					int upred=0;
					int vpred=0;
#endif
#ifdef __GNUC__
#pragma GCC unroll 2
#endif
					for(int k=0;k<blocksize;++k)
					{
#ifdef USE_CG
						int yN=yNptr[0];
						int uN=uNptr[0];
						int vN=vNptr[0];
						int yW=ycurrptr[-1];
						int uW=ucurrptr[-1];
						int vW=vcurrptr[-1];
						int ymax=yN, ymin=yW;
						int umax=uN, umin=uW;
						int vmax=vN, vmin=vW;
						if(yN<yW)ymin=yN, ymax=yW;
						if(uN<uW)umin=uN, umax=uW;
						if(vN<vW)vmin=vN, vmax=vW;
					//	int ypred=2*yN+((yW-yN)*Wcoeff>>7)-yNptr[-1];
					//	int upred=2*uN+((uW-uN)*Wcoeff>>7)-uNptr[-1];
					//	int vpred=2*vN+((vW-vN)*Wcoeff>>7)-vNptr[-1];
						int ypred=yN+yW-yNptr[-1];
						int upred=uN+uW-uNptr[-1];
						int vpred=vN+vW-vNptr[-1];
						CLAMP2(ypred, ymin, ymax);
						CLAMP2(upred, umin, umax);
						CLAMP2(vpred, vmin, vmax);
#elif defined USE_WGCG
						int yN=yNptr[0];
						int uN=uNptr[0];
						int vN=vNptr[0];
						int yW=ycurrptr[-1];
						int uW=ucurrptr[-1];
						int vW=vcurrptr[-1];
						int ygx=abs(yW-ycurrptr[-2])+abs(yN-yNptr[-1])+abs(yNptr[1]-yNptr[0])+1;//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
						int ugx=abs(uW-ucurrptr[-2])+abs(uN-uNptr[-1])+abs(uNptr[1]-uNptr[0])+1;
						int vgx=abs(vW-vcurrptr[-2])+abs(vN-vNptr[-1])+abs(vNptr[1]-vNptr[0])+1;
						int ygy=abs(yW-yNptr[-1])+abs(yNptr[0]-yNNptr[0])+abs(yNptr[1]-yNNptr[1])+1;//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
						int ugy=abs(uW-uNptr[-1])+abs(uNptr[0]-uNNptr[0])+abs(uNptr[1]-uNNptr[1])+1;
						int vgy=abs(vW-vNptr[-1])+abs(vNptr[0]-vNNptr[0])+abs(vNptr[1]-vNNptr[1])+1;
						int ysum=ygx+ygy;
						int usum=ugx+ugy;
						int vsum=vgx+vgy;
						int ypred=(ygx+ysum)*yN+(ygy+ysum)*yW-ysum*yNptr[-1];
						int upred=(ugx+usum)*uN+(ugy+usum)*uW-usum*uNptr[-1];
						int vpred=(vgx+vsum)*vN+(vgy+vsum)*vW-vsum*vNptr[-1];
						ypred=(long long)ypred*divtable[ysum]>>24;
						upred=(long long)upred*divtable[usum]>>24;
						vpred=(long long)vpred*divtable[vsum]>>24;
					//	int ypred=((ygx+ysum)*yN+(ygy+ysum)*yW-ysum*yNptr[-1])/(2*ysum);//pred=((gx+sum)*N+(gy+sum)*W-sum*NW)/(2*sum)
					//	int upred=((ugx+usum)*uN+(ugy+usum)*uW-usum*uNptr[-1])/(2*usum);
					//	int vpred=((vgx+vsum)*vN+(vgy+vsum)*vW-vsum*vNptr[-1])/(2*vsum);
						//int ypred=(ygx*yN+ygy*yW)/ysum;//pred=(gx*N+gy*W)/sum
						//int upred=(ugx*uN+ugy*uW)/usum;
						//int vpred=(vgx*vN+vgy*vW)/vsum;
						int ymax=yN, ymin=yW;
						int umax=uN, umin=uW;
						int vmax=vN, vmin=vW;
						if(yN<yW)ymin=yN, ymax=yW;
						if(uN<uW)umin=uN, umax=uW;
						if(vN<vW)vmin=vN, vmax=vW;
						int yNE=yNptr[1];
						int uNE=uNptr[1];
						int vNE=vNptr[1];
						if(ymin>yNE)ymin=yNE;
						if(umin>uNE)umin=uNE;
						if(vmin>vNE)vmin=vNE;
						if(ymax<yNE)ymax=yNE;
						if(umax<uNE)umax=uNE;
						if(vmax<vNE)vmax=vNE;
						CLAMP2(ypred, ymin, ymax);
						CLAMP2(upred, umin, umax);
						CLAMP2(vpred, vmin, vmax);
						++yNNptr;
						++uNNptr;
						++vNNptr;
#endif
						char yuv[4]=
						{
							*ysrcptr++,
							*usrcptr++,
							*vsrcptr++,
							0,
						};
						yuv[0]+=ypred;
						yuv[1]+=upred;
						yuv[2]+=vpred;
#ifndef USE_W
						++yNptr;
						++uNptr;
						++vNptr;
#endif
						*ycurrptr++=yuv[0];//store YUV
						*ucurrptr++=yuv[1];
						*vcurrptr++=yuv[2];
#ifdef USE_W
						int yW=yuv[0];
						int uW=yuv[1];
						int vW=yuv[2];
#endif
#ifdef DISABLE_ANALYSIS
						yuv[2]+=yuv[1]>>2;
						yuv[0]-=(yuv[1]+yuv[2])>>2;
						yuv[2]+=yuv[0];
						yuv[1]+=yuv[0];
						(void)vhelpidx;
						(void)uhelpidx;
#else
						yuv[1]+=yuv[uhelpidx];//need int16 buffer for cRCT
						yuv[2]+=yuv[vhelpidx];
#endif
						*ydstptr++=yuv[0];//store RGB
						*udstptr++=yuv[1];
						*vdstptr++=yuv[2];
						//*ycurrptr++=(char)(yuv[0]+ypred);
						//upred+=yuv[uhelpidx];
						//CLAMP2(upred, -128, 127);
						//*ucurrptr++=(char)(yuv[1]+upred);
						//vpred+=yuv[vhelpidx];
						//CLAMP2(vpred, -128, 127);
						//*vcurrptr++=(char)(yuv[2]+vpred);
#ifdef USE_W
						ypred=yW;
						upred=uW;
						vpred=vW;
#endif
					}
				}
				PROF(dec_predict, blockbytes);

				//interleave
				{
					const char
						*ydstptr=im2[0],
						*udstptr=im2[1],
						*vdstptr=im2[2];
					for(int ky=0;ky<dy;++ky)
					{
						unsigned char *row=image+3*(iw*ky+x1);
						for(int kx=0;kx<dx;++kx)
						{
							*row++=(unsigned char)(*ydstptr++ + 128);
							*row++=(unsigned char)(*udstptr++ + 128);
							*row++=(unsigned char)(*vdstptr++ + 128);
						}
					}
				}
				PROF(dec_interleave, blockbytes);
			}
			fwrite(image, 1, (ptrdiff_t)3*iw*dy, fdst);
			PROF(dec_write, (ptrdiff_t)3*iw*dy);
		}
	}
	//if(fwd)//
	//	printf("%d\n", maxsum);//
	fclose(fsrc);
	fclose(fdst);
	free(image);
	prof_print();
	return 0;
}
