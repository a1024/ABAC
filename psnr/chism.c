#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<sys/stat.h>
#if defined _WIN32 || defined WIN32
#include<Windows.h>
#else
#include<time.h>
#endif

	#define PROFILE_TIME

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
#if defined _WIN32 || defined WIN32
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
typedef struct _SpeedProfilerInfo
{
	double t;
	ptrdiff_t size;
	const char *msg;
} SpeedProfilerInfo;
#define PROF_CAP 128
static double prof_timestamp=0;
static SpeedProfilerInfo prof_data[PROF_CAP]={0};
static int prof_count=0;
static void prof_checkpoint(ptrdiff_t size, const char *msg)
{
	double t2=time_sec();
	if(prof_timestamp)
	{
		SpeedProfilerInfo *info=prof_data+prof_count++;
		if(prof_count>PROF_CAP)
		{
			CRASH("Profiler OOB");
			return;
		}
		info->t=t2-prof_timestamp;
		info->size=size;
		info->msg=msg;
	}
	prof_timestamp=t2;
}
static void prof_print(ptrdiff_t usize)
{
	static char buf[2048]={0};
	double timesum=0, tmax=0;
	for(int k=0;k<prof_count;++k)
	{
		double t=prof_data[k].t;
		timesum+=t;
		if(tmax<t)
			tmax=t;
	}
	int prev=0;
	double csum=0;
	printf("1 char = 4 ms\n");
	printf("|");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		int curr=(int)(csum*250);//fixed scale
		int space=curr-prev;
		int len=0;
		if(info->msg)
			len=(int)strlen(info->msg);
		if(space>2047)//printf("%*s", HUGE, ""); CRASHES
			space=2047;
		if(info->msg&&space>=len)
		{
			int labelstart=(space-len)>>1;
			int labelend=labelstart+len;

			memset(buf, '-', labelstart);
			buf[labelstart]=0;
			printf("%s%s", buf, info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			printf("%s", buf);
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			printf("%s", buf);
		}
		printf("|");
		prev=curr;
	}
	printf("\n");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
		if(info->size)
			printf(" %16.6lf MB/s %12.6lf ms/MB %10td bytes ", info->size/(info->t*1024*1024), info->t*1000*1024*1024/info->size, info->size);
		if(info->msg)
			printf("%-20s", info->msg);
		else// if(nstars)
			printf("%-20s", "");
		printf("\n");
	}
	printf("\n");
	printf("%16.7lf ms %12.6lf MB/s %12.6lf mb/MB Total\n", timesum*1000, usize/(timesum*1024*1024), timesum*1024*1024*1000/usize);
	printf("\n");
	prof_count=0;
	prof_timestamp=0;
}
#else
#define prof_checkpoint(...)
#define prof_print(...)
#endif

static uint8_t* ppm_load(const char *fn, int *ret_iw, int *ret_ih)
{
	FILE *fsrc;
	int iw, ih, c;

	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"\n", fn);
		return 0;
	}
	c=0;
	fread(&c, 1, 2, fsrc);
	if(c!=('P'|'6'<<8))
	{
		CRASH("Invalid PPM \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		CRASH("Invalid PPM file");
		return 0;
	}
	c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	iw=0;
	while((unsigned)(c-'0')<10)
	{
		iw=10*iw+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	ih=0;
	while((unsigned)(c-'0')<10)
	{
		ih=10*ih+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
	{
		CRASH("Invalid PPM file");
		return 0;
	}
	{
		ptrdiff_t srcsize=(ptrdiff_t)3*iw*ih;
		uint8_t *buffer=(uint8_t*)malloc(srcsize);
		if(!buffer)
		{
			CRASH("Alloc error");
			return 0;
		}
		fread(buffer, 1, srcsize, fsrc);
		fclose(fsrc);
		if(ret_iw)*ret_iw=iw;
		if(ret_ih)*ret_ih=ih;
		return buffer;
	}
}
static uint8_t* pgm_load(const char *fn, int *ret_iw, int *ret_ih)
{
	FILE *fsrc;
	int iw, ih, c;

	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"\n", fn);
		return 0;
	}
	c=0;
	fread(&c, 1, 2, fsrc);
	if(c!=('P'|'5'<<8))
	{
		CRASH("Invalid PGM \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		CRASH("Invalid PGM file");
		return 0;
	}
	c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	iw=0;
	while((unsigned)(c-'0')<10)
	{
		iw=10*iw+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	ih=0;
	while((unsigned)(c-'0')<10)
	{
		ih=10*ih+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
	{
		CRASH("Invalid PGM file");
		return 0;
	}
	{
		ptrdiff_t srcsize=(ptrdiff_t)iw*ih;
		uint8_t *buffer=(uint8_t*)malloc(srcsize);
		if(!buffer)
		{
			CRASH("Alloc error");
			return 0;
		}
		fread(buffer, 1, srcsize, fsrc);
		fclose(fsrc);
		if(ret_iw)*ret_iw=iw;
		if(ret_ih)*ret_ih=ih;
		return buffer;
	}
}
static void ppm_save(const char *fn, uint8_t *buffer, int iw, int ih, int overwrite)
{
	struct stat info={0};
	int error=stat(fn, &info);
	if(!overwrite&&error!=-1)
		return;
	{
		FILE *fdst=fopen(fn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", fn);
			return;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(buffer, 1, (ptrdiff_t)3*iw*ih, fdst);
		fclose(fdst);
	}
}
static void pgm_save(const char *fn, uint8_t *buffer, int iw, int ih, int overwrite)
{
	struct stat info={0};
	int error=stat(fn, &info);
	if(!overwrite&&error!=-1)
		return;
	{
		FILE *fdst=fopen(fn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", fn);
			return;
		}
		fprintf(fdst, "P5\n%d %d\n255\n", iw, ih);
		fwrite(buffer, 1, (ptrdiff_t)iw*ih, fdst);
		fclose(fdst);
	}
}
static ptrdiff_t isplit(int overwrite, const char *iname, const char *rname, const char *gname, const char *bname)
{
	uint8_t *ibuf=0, *rbuf=0, *gbuf=0, *bbuf=0;
	ptrdiff_t res, srcsize;
	int rowstride, iw, ih;
	
	{
		struct stat info[4]={0};
		char absent[4]={0};

		absent[0]=stat(iname, info+0)==-1;
		if(!overwrite)
		{
			absent[1]=stat(rname, info+1)==-1;
			absent[2]=stat(gname, info+2)==-1;
			absent[3]=stat(bname, info+3)==-1;
		}

		if(absent[0])
			printf("Cannot stat \"%s\"\n", iname);
		if(!overwrite)
		{
			if(!absent[1])printf("Already exists \"%s\"\n", rname);
			if(!absent[2])printf("Already exists \"%s\"\n", gname);
			if(!absent[3])printf("Already exists \"%s\"\n", bname);
		}
		if(absent[0])
			return 0;
		if(!overwrite)
		{
			if(!absent[1]||!absent[2]||!absent[3])
				return 0;
		}
	}
	ibuf=ppm_load(iname, &iw, &ih);
	res=(ptrdiff_t)iw*ih;
	srcsize=(ptrdiff_t)3*res;
	prof_checkpoint(srcsize, "load");
	
	rowstride=3*iw;
	rbuf=(uint8_t*)malloc(srcsize);
	if(!rbuf)
	{
		CRASH("Alloc error");
		return 0;
	}
	gbuf=rbuf+res;
	bbuf=gbuf+res;
	prof_checkpoint(srcsize, "alloc");
	{
		uint8_t
			*iptr=ibuf,
			*rptr=rbuf,
			*gptr=gbuf,
			*bptr=bbuf;
		int kx, ky;

		for(ky=0;ky<ih;++ky)
		{
			//AVX2
			//FEDCBA9876543210FEDCBA9876543210 FEDCBA9876543210FEDCBA9876543210
			//rbgrbgrbgrbgrbgrbgrbgrbgrbgrbgrb grbgrbgrbgrbgrbgrbgrbgrbgrbgrbgr
			//grbgrbgrbgrbgrbgrbgrbgrbgrbgrbgr bgrbgrbgrbgrbgrbgrbgrbgrbgrbgrbg
			//bgrbgrbgrbgrbgrbgrbgrbgrbgrbgrbg rbgrbgrbgrbgrbgrbgrbgrbgrbgrbgrb
			for(kx=0;kx<=rowstride-32*3;kx+=32*3)
			{
				uint8_t b2[32*3];
				memcpy(b2, iptr, sizeof(char[32*3]));
				rptr[0x00]=iptr[0x00*3+0]; gptr[0x00]=iptr[0x00*3+1]; bptr[0x00]=iptr[0x00*3+2];
				rptr[0x01]=iptr[0x01*3+0]; gptr[0x01]=iptr[0x01*3+1]; bptr[0x01]=iptr[0x01*3+2];
				rptr[0x02]=iptr[0x02*3+0]; gptr[0x02]=iptr[0x02*3+1]; bptr[0x02]=iptr[0x02*3+2];
				rptr[0x03]=iptr[0x03*3+0]; gptr[0x03]=iptr[0x03*3+1]; bptr[0x03]=iptr[0x03*3+2];
				rptr[0x04]=iptr[0x04*3+0]; gptr[0x04]=iptr[0x04*3+1]; bptr[0x04]=iptr[0x04*3+2];
				rptr[0x05]=iptr[0x05*3+0]; gptr[0x05]=iptr[0x05*3+1]; bptr[0x05]=iptr[0x05*3+2];
				rptr[0x06]=iptr[0x06*3+0]; gptr[0x06]=iptr[0x06*3+1]; bptr[0x06]=iptr[0x06*3+2];
				rptr[0x07]=iptr[0x07*3+0]; gptr[0x07]=iptr[0x07*3+1]; bptr[0x07]=iptr[0x07*3+2];
				rptr[0x08]=iptr[0x08*3+0]; gptr[0x08]=iptr[0x08*3+1]; bptr[0x08]=iptr[0x08*3+2];
				rptr[0x09]=iptr[0x09*3+0]; gptr[0x09]=iptr[0x09*3+1]; bptr[0x09]=iptr[0x09*3+2];
				rptr[0x0A]=iptr[0x0A*3+0]; gptr[0x0A]=iptr[0x0A*3+1]; bptr[0x0A]=iptr[0x0A*3+2];
				rptr[0x0B]=iptr[0x0B*3+0]; gptr[0x0B]=iptr[0x0B*3+1]; bptr[0x0B]=iptr[0x0B*3+2];
				rptr[0x0C]=iptr[0x0C*3+0]; gptr[0x0C]=iptr[0x0C*3+1]; bptr[0x0C]=iptr[0x0C*3+2];
				rptr[0x0D]=iptr[0x0D*3+0]; gptr[0x0D]=iptr[0x0D*3+1]; bptr[0x0D]=iptr[0x0D*3+2];
				rptr[0x0E]=iptr[0x0E*3+0]; gptr[0x0E]=iptr[0x0E*3+1]; bptr[0x0E]=iptr[0x0E*3+2];
				rptr[0x0F]=iptr[0x0F*3+0]; gptr[0x0F]=iptr[0x0F*3+1]; bptr[0x0F]=iptr[0x0F*3+2];
				rptr[0x10]=iptr[0x10*3+0]; gptr[0x10]=iptr[0x10*3+1]; bptr[0x10]=iptr[0x10*3+2];
				rptr[0x11]=iptr[0x11*3+0]; gptr[0x11]=iptr[0x11*3+1]; bptr[0x11]=iptr[0x11*3+2];
				rptr[0x12]=iptr[0x12*3+0]; gptr[0x12]=iptr[0x12*3+1]; bptr[0x12]=iptr[0x12*3+2];
				rptr[0x13]=iptr[0x13*3+0]; gptr[0x13]=iptr[0x13*3+1]; bptr[0x13]=iptr[0x13*3+2];
				rptr[0x14]=iptr[0x14*3+0]; gptr[0x14]=iptr[0x14*3+1]; bptr[0x14]=iptr[0x14*3+2];
				rptr[0x15]=iptr[0x15*3+0]; gptr[0x15]=iptr[0x15*3+1]; bptr[0x15]=iptr[0x15*3+2];
				rptr[0x16]=iptr[0x16*3+0]; gptr[0x16]=iptr[0x16*3+1]; bptr[0x16]=iptr[0x16*3+2];
				rptr[0x17]=iptr[0x17*3+0]; gptr[0x17]=iptr[0x17*3+1]; bptr[0x17]=iptr[0x17*3+2];
				rptr[0x18]=iptr[0x18*3+0]; gptr[0x18]=iptr[0x18*3+1]; bptr[0x18]=iptr[0x18*3+2];
				rptr[0x19]=iptr[0x19*3+0]; gptr[0x19]=iptr[0x19*3+1]; bptr[0x19]=iptr[0x19*3+2];
				rptr[0x1A]=iptr[0x1A*3+0]; gptr[0x1A]=iptr[0x1A*3+1]; bptr[0x1A]=iptr[0x1A*3+2];
				rptr[0x1B]=iptr[0x1B*3+0]; gptr[0x1B]=iptr[0x1B*3+1]; bptr[0x1B]=iptr[0x1B*3+2];
				rptr[0x1C]=iptr[0x1C*3+0]; gptr[0x1C]=iptr[0x1C*3+1]; bptr[0x1C]=iptr[0x1C*3+2];
				rptr[0x1D]=iptr[0x1D*3+0]; gptr[0x1D]=iptr[0x1D*3+1]; bptr[0x1D]=iptr[0x1D*3+2];
				rptr[0x1E]=iptr[0x1E*3+0]; gptr[0x1E]=iptr[0x1E*3+1]; bptr[0x1E]=iptr[0x1E*3+2];
				rptr[0x1F]=iptr[0x1F*3+0]; gptr[0x1F]=iptr[0x1F*3+1]; bptr[0x1F]=iptr[0x1F*3+2];
				iptr+=32*3;
				rptr+=32;
				gptr+=32;
				bptr+=32;
			}
			for(;kx<rowstride;kx+=3)
			{
				*rptr++=iptr[0];
				*gptr++=iptr[1];
				*bptr++=iptr[2];
				iptr+=3;
			}
		}
	}
	prof_checkpoint(srcsize, "split");
	free(ibuf);
	pgm_save(rname, rbuf, iw, ih, overwrite);
	pgm_save(gname, gbuf, iw, ih, overwrite);
	pgm_save(bname, bbuf, iw, ih, overwrite);
	free(rbuf);
	prof_checkpoint(srcsize, "save");
	return (ptrdiff_t)3*iw*ih;
}
static ptrdiff_t ijoin(int overwrite, const char *iname, const char *rname, const char *gname, const char *bname)
{
	uint8_t *ibuf=0, *rbuf=0, *gbuf=0, *bbuf=0;
	ptrdiff_t res, srcsize;
	int rowstride, iw, ih;
	
	{
		struct stat info[4]={0};
		char absent[4]={0};

		if(!overwrite)
			absent[0]=stat(iname, info+0)==-1;
		absent[1]=stat(rname, info+1)==-1;
		absent[2]=stat(gname, info+2)==-1;
		absent[3]=stat(bname, info+3)==-1;
		
		if(!overwrite)
		{
			if(!absent[0])
				printf("Already exists \"%s\"\n", iname);
		}
		if(absent[1])printf("Cannot stat \"%s\"\n", rname);
		if(absent[2])printf("Cannot stat \"%s\"\n", gname);
		if(absent[3])printf("Cannot stat \"%s\"\n", bname);
		if(!overwrite)
		{
			if(!absent[0])
				return 0;
		}
		if(absent[1]||absent[2]||absent[3])
			return 0;
	}
	{
		int iw2, ih2, iw3, ih3;

		rbuf=pgm_load(rname, &iw, &ih);
		gbuf=pgm_load(gname, &iw2, &ih2);
		bbuf=pgm_load(bname, &iw3, &ih3);
		if(iw!=iw2||iw!=iw3||ih!=ih2||ih!=ih3)
		{
			CRASH("Dimension mismatch  WH  r %d*%d  g %d*%d  b %d*%d"
				, iw, ih
				, iw2, ih2
				, iw3, ih3
			);
			return 0;
		}
	}
	res=(ptrdiff_t)iw*ih;
	srcsize=(ptrdiff_t)3*res;
	prof_checkpoint(srcsize, "load");
	
	rowstride=3*iw;
	ibuf=(uint8_t*)malloc(srcsize);
	if(!ibuf)
	{
		CRASH("Alloc error");
		return 0;
	}
	prof_checkpoint(srcsize, "alloc");
	{
		uint8_t
			*iptr=ibuf,
			*rptr=rbuf,
			*gptr=gbuf,
			*bptr=bbuf;
		int kx, ky;

		for(ky=0;ky<ih;++ky)
		{
			//AVX2
			//FEDCBA9876543210FEDCBA9876543210 FEDCBA9876543210FEDCBA9876543210
			//rbgrbgrbgrbgrbgrbgrbgrbgrbgrbgrb grbgrbgrbgrbgrbgrbgrbgrbgrbgrbgr
			//grbgrbgrbgrbgrbgrbgrbgrbgrbgrbgr bgrbgrbgrbgrbgrbgrbgrbgrbgrbgrbg
			//bgrbgrbgrbgrbgrbgrbgrbgrbgrbgrbg rbgrbgrbgrbgrbgrbgrbgrbgrbgrbgrb
			for(kx=0;kx<=rowstride-32*3;kx+=32*3)
			{
				uint8_t b2[32*3];
				memcpy(b2, iptr, sizeof(char[32*3]));
				iptr[0x00*3+0]=rptr[0x00]; iptr[0x00*3+1]=gptr[0x00]; iptr[0x00*3+2]=bptr[0x00];
				iptr[0x01*3+0]=rptr[0x01]; iptr[0x01*3+1]=gptr[0x01]; iptr[0x01*3+2]=bptr[0x01];
				iptr[0x02*3+0]=rptr[0x02]; iptr[0x02*3+1]=gptr[0x02]; iptr[0x02*3+2]=bptr[0x02];
				iptr[0x03*3+0]=rptr[0x03]; iptr[0x03*3+1]=gptr[0x03]; iptr[0x03*3+2]=bptr[0x03];
				iptr[0x04*3+0]=rptr[0x04]; iptr[0x04*3+1]=gptr[0x04]; iptr[0x04*3+2]=bptr[0x04];
				iptr[0x05*3+0]=rptr[0x05]; iptr[0x05*3+1]=gptr[0x05]; iptr[0x05*3+2]=bptr[0x05];
				iptr[0x06*3+0]=rptr[0x06]; iptr[0x06*3+1]=gptr[0x06]; iptr[0x06*3+2]=bptr[0x06];
				iptr[0x07*3+0]=rptr[0x07]; iptr[0x07*3+1]=gptr[0x07]; iptr[0x07*3+2]=bptr[0x07];
				iptr[0x08*3+0]=rptr[0x08]; iptr[0x08*3+1]=gptr[0x08]; iptr[0x08*3+2]=bptr[0x08];
				iptr[0x09*3+0]=rptr[0x09]; iptr[0x09*3+1]=gptr[0x09]; iptr[0x09*3+2]=bptr[0x09];
				iptr[0x0A*3+0]=rptr[0x0A]; iptr[0x0A*3+1]=gptr[0x0A]; iptr[0x0A*3+2]=bptr[0x0A];
				iptr[0x0B*3+0]=rptr[0x0B]; iptr[0x0B*3+1]=gptr[0x0B]; iptr[0x0B*3+2]=bptr[0x0B];
				iptr[0x0C*3+0]=rptr[0x0C]; iptr[0x0C*3+1]=gptr[0x0C]; iptr[0x0C*3+2]=bptr[0x0C];
				iptr[0x0D*3+0]=rptr[0x0D]; iptr[0x0D*3+1]=gptr[0x0D]; iptr[0x0D*3+2]=bptr[0x0D];
				iptr[0x0E*3+0]=rptr[0x0E]; iptr[0x0E*3+1]=gptr[0x0E]; iptr[0x0E*3+2]=bptr[0x0E];
				iptr[0x0F*3+0]=rptr[0x0F]; iptr[0x0F*3+1]=gptr[0x0F]; iptr[0x0F*3+2]=bptr[0x0F];
				iptr[0x10*3+0]=rptr[0x10]; iptr[0x10*3+1]=gptr[0x10]; iptr[0x10*3+2]=bptr[0x10];
				iptr[0x11*3+0]=rptr[0x11]; iptr[0x11*3+1]=gptr[0x11]; iptr[0x11*3+2]=bptr[0x11];
				iptr[0x12*3+0]=rptr[0x12]; iptr[0x12*3+1]=gptr[0x12]; iptr[0x12*3+2]=bptr[0x12];
				iptr[0x13*3+0]=rptr[0x13]; iptr[0x13*3+1]=gptr[0x13]; iptr[0x13*3+2]=bptr[0x13];
				iptr[0x14*3+0]=rptr[0x14]; iptr[0x14*3+1]=gptr[0x14]; iptr[0x14*3+2]=bptr[0x14];
				iptr[0x15*3+0]=rptr[0x15]; iptr[0x15*3+1]=gptr[0x15]; iptr[0x15*3+2]=bptr[0x15];
				iptr[0x16*3+0]=rptr[0x16]; iptr[0x16*3+1]=gptr[0x16]; iptr[0x16*3+2]=bptr[0x16];
				iptr[0x17*3+0]=rptr[0x17]; iptr[0x17*3+1]=gptr[0x17]; iptr[0x17*3+2]=bptr[0x17];
				iptr[0x18*3+0]=rptr[0x18]; iptr[0x18*3+1]=gptr[0x18]; iptr[0x18*3+2]=bptr[0x18];
				iptr[0x19*3+0]=rptr[0x19]; iptr[0x19*3+1]=gptr[0x19]; iptr[0x19*3+2]=bptr[0x19];
				iptr[0x1A*3+0]=rptr[0x1A]; iptr[0x1A*3+1]=gptr[0x1A]; iptr[0x1A*3+2]=bptr[0x1A];
				iptr[0x1B*3+0]=rptr[0x1B]; iptr[0x1B*3+1]=gptr[0x1B]; iptr[0x1B*3+2]=bptr[0x1B];
				iptr[0x1C*3+0]=rptr[0x1C]; iptr[0x1C*3+1]=gptr[0x1C]; iptr[0x1C*3+2]=bptr[0x1C];
				iptr[0x1D*3+0]=rptr[0x1D]; iptr[0x1D*3+1]=gptr[0x1D]; iptr[0x1D*3+2]=bptr[0x1D];
				iptr[0x1E*3+0]=rptr[0x1E]; iptr[0x1E*3+1]=gptr[0x1E]; iptr[0x1E*3+2]=bptr[0x1E];
				iptr[0x1F*3+0]=rptr[0x1F]; iptr[0x1F*3+1]=gptr[0x1F]; iptr[0x1F*3+2]=bptr[0x1F];
				iptr+=32*3;
				rptr+=32;
				gptr+=32;
				bptr+=32;
			}
			for(;kx<rowstride;kx+=3)
			{
				iptr[0]=*rptr++;
				iptr[1]=*gptr++;
				iptr[2]=*bptr++;
				iptr+=3;
			}
		}
	}
	prof_checkpoint(srcsize, "join");
	free(rbuf);
	free(gbuf);
	free(bbuf);
	ppm_save(iname, ibuf, iw, ih, overwrite);
	free(ibuf);
	prof_checkpoint(srcsize, "save");
	return (ptrdiff_t)3*iw*ih;
}

static void print_usage(const char *argv0)
{
	printf(
		"Usage:  \"%s\"  s|j  image.ppm  red.pgm  green.pgm  blue.pgm\n"
		"  s     split\n"
		"  S     split and OVERWRITE\n"
		"  j     join\n"
		"  J     join and OVERWRITE\n"
		, argv0
	);
}
int main(int argc, char **argv)
{
	//volatile double elapsed;
	ptrdiff_t processed;
	const char *iname, *rname, *gname, *bname;
	int overwrite;
	char c;

#if defined RELEASE || defined __GNUC__
	if(argc!=6||!argv[1][0]||argv[1][1])
	{
		print_usage(argv[0]);
		return 1;
	}
	c=argv[1][0];
	iname=argv[2];
	rname=argv[3];
	gname=argv[4];
	bname=argv[5];
#else
	c='S';
	iname="D:/ML/serpinski.ppm";
//	iname="D:/ML/zzz.ppm";
	rname="D:/ML/r.pgm";
	gname="D:/ML/g.pgm";
	bname="D:/ML/b.pgm";
#endif
	
	prof_checkpoint(0, 0);
	//elapsed=time_sec();
	overwrite=!(c>>5&1);
	c&=0xDF;
	if(c=='S')
		processed=isplit(overwrite, iname, rname, gname, bname);
	else if(c=='J')
		processed=ijoin(overwrite, iname, rname, gname, bname);
	else
	{
		print_usage(argv[0]);
		return 1;
	}
	prof_print(processed);
	//elapsed=time_sec()-elapsed;
	//printf("%12.6lf sec  %12.6lf MB/s\n", elapsed, processed/(elapsed*1024*1024));
	return 0;
}
