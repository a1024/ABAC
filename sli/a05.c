#include"fse/fse.h"
#include"fpc/fpc.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#elif defined __GNUC__
#include<time.h>
#endif


	#define LOUD
	#define PROFILE_TIME

	#define USE_FPC


#define ENC_DY 256
#define ENC_DX 512
#define YPAD 2

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
	PROFLABEL(enc_predict)\
	PROFLABEL(enc_entropycoding)\
	PROFLABEL(enc_write)\
	PROFLABEL(enc_count)\
	PROFLABEL(dec_read)\
	PROFLABEL(dec_entropycoding)\
	PROFLABEL(dec_predict)\
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
static char im2[3][(ENC_DY+YPAD)*ENC_DX];
static char cim[3][(ENC_DY+YPAD)*ENC_DX];
int a05_codec(int argc, char **argv)
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
						*ydstptr=im2[0],
						*udstptr=im2[1],
						*vdstptr=im2[2];
					char W[3]={0};
					for(int ky=0;ky<dy;++ky)
					{
						const unsigned char *row=image+3*(iw*ky+x1);
						for(int kx=0;kx<dx;++kx)
						{
							char yuv[]=
							{
								row[1]-128,//y = g
								row[2]-128,//u = b
								row[0]-128,//v = r
							};
							row+=3;
							yuv[1]-=yuv[0];
							yuv[2]-=yuv[0];
							//yuv[0]+=(yuv[1]+yuv[2])>>2;
							//yuv[2]-=yuv[1]>>2;
							*ydstptr++=yuv[0]-W[0];
							*udstptr++=yuv[1]-W[1];
							*vdstptr++=yuv[2]-W[2];
							W[0]=yuv[0];
							W[1]=yuv[1];
							W[2]=yuv[2];
						}
					}
				}
				PROF(enc_predict, blockbytes);

				//write
#ifdef USE_FPC
				int bsize=blocksize;
				if(bsize>0x4000)
					bsize=0x4000;
				size_t csize0=FPC_compress(cim[0], im2[0], blocksize, bsize);
				size_t csize1=FPC_compress(cim[1], im2[1], blocksize, bsize);
				size_t csize2=FPC_compress(cim[2], im2[2], blocksize, bsize);
#else
				size_t csize0=FSE_compress(cim[0], sizeof(cim[0]), im2[0], blocksize);
				size_t csize1=FSE_compress(cim[1], sizeof(cim[1]), im2[1], blocksize);
				size_t csize2=FSE_compress(cim[2], sizeof(cim[2]), im2[2], blocksize);
#endif
				if(!csize0||!csize1||!csize2)
				{
					CRASH("FPC_compress %d %d %d", csize0, csize1, csize2);
					return 1;
				}
				PROF(enc_entropycoding, blockbytes);
				(void)blockbytes;
				fwrite(&csize0, 1, 4, fdst);
				fwrite(&csize1, 1, 4, fdst);
				fwrite(&csize2, 1, 4, fdst);
				fwrite(cim[0], 1, csize0, fdst);
				fwrite(cim[1], 1, csize1, fdst);
				fwrite(cim[2], 1, csize2, fdst);
				PROF(enc_write, csize0+csize1+csize2+3*4);
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
				int csize0=0, csize1=0, csize2=0;
				ptrdiff_t nread=0;
				nread+=fread(&csize0, 1, 4, fsrc);
				nread+=fread(&csize1, 1, 4, fsrc);
				nread+=fread(&csize2, 1, 4, fsrc);
				if((int)nread!=3*4)
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
				PROF(dec_read, csize0+csize1+csize2+3*4);

				//decompress
#ifdef USE_FPC
				size_t usize0=FPC_decompress(im2[0], blocksize, cim[0], csize0);
				size_t usize1=FPC_decompress(im2[1], blocksize, cim[1], csize1);
				size_t usize2=FPC_decompress(im2[2], blocksize, cim[2], csize2);
#else
				size_t usize0=FSE_decompress(im2[0], sizeof(im2[0]), cim[0], csize0);
				size_t usize1=FSE_decompress(im2[1], sizeof(im2[1]), cim[1], csize1);
				size_t usize2=FSE_decompress(im2[2], sizeof(im2[2]), cim[2], csize2);
#endif
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
				PROF(dec_entropycoding, usize0+usize1+usize2);

				//interleave
				{
					const char
						*ydstptr=im2[0],
						*udstptr=im2[1],
						*vdstptr=im2[2];
					char W[3]={0};
					for(int ky=0;ky<dy;++ky)
					{
						unsigned char *row=image+3*(iw*ky+x1);
						for(int kx=0;kx<dx;++kx)
						{
							char yuv[]=
							{
								*ydstptr++,
								*udstptr++,
								*vdstptr++,
							};
							yuv[0]+=W[0];
							yuv[1]+=W[1];
							yuv[2]+=W[2];
							W[0]=yuv[0];
							W[1]=yuv[1];
							W[2]=yuv[2];
							//yuv[2]+=yuv[1]>>2;
							//yuv[0]-=(yuv[1]+yuv[2])>>2;
							yuv[2]+=yuv[0];
							yuv[1]+=yuv[0];
							*row++=(unsigned char)(yuv[2]+128);//v = r
							*row++=(unsigned char)(yuv[0]+128);//y = g
							*row++=(unsigned char)(yuv[1]+128);//u = b
						}
					}
				}
				PROF(dec_predict, blockbytes);
				(void)blockbytes;
			}
			fwrite(image, 1, (ptrdiff_t)3*iw*dy, fdst);
			PROF(dec_write, (ptrdiff_t)3*iw*dy);
		}
	}
	fclose(fsrc);
	fclose(fdst);
	free(image);
	prof_print();
	(void)time_sec;
	return 0;
}
