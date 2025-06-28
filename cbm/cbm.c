#include"util.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<time.h>
#include<Windows.h>
#include<Psapi.h>
#include<tlhelp32.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define TRACK_THREADS

//	#define PSNR_MAX_8 128
//	#define PSNR_MAX_16 0x8000
	#define PSNR_MAX_8 255
	#define PSNR_MAX_16 0xFFFF

static char g_buf2[8192]={0};

typedef union _DateTime
{
	struct
	{
		uint8_t ds, second, minute, hour, day, month;
		uint16_t year;
	};
	uint64_t timestamp;
} DateTime;
typedef enum _CmdFlags
{
	CMDFLAG_VERIFY_BITEXACT=0,
	CMDFLAG_NOP=1,
	CMDFLAG_PSNR_8BIT=2,
	CMDFLAG_PSNR_16BIT=3,
	CMDFLAG_SSIM_PPM=4,

	CMDFLAG_PRINT_RIVALS=5,
} CmdFlags;

static int acme_getline(char *buf, int len, FILE *f)
{
	memset(buf, '\n', len);
	fgets(buf, len, f);
	int k=0;
	for(;k<len&&buf[k]!='\n';++k);
	buf[k]=0;
	return k;
}
static int acme_strnimatch(const char *s1, ptrdiff_t len1, const char *s2, ptrdiff_t len2)//return 1: match		FIXME return ASCII order & index of first difference
{
	if(len1!=len2)
		return 0;
	const char *end1=s1+len1;
	while(s1<end1&&tolower(*s1)==tolower(*s2))++s1, ++s2;//check then increment  (blind increment misses the last character at s1==end1)
	return s1==end1;
}
static void verify_files(const char *fn0, const char *fn1)
{
	int broken=0;
	ArrayHandle data0=load_file(fn0, 1, 32, 1);
	ArrayHandle data1=load_file(fn1, 1, 32, 1);
	if(data0->count!=data1->count)
	{
		broken=1;
		printf("  verify_files %zd vs %zd\n", data0->count, data1->count);
	}
	ptrdiff_t size=MINVAR(data0->count, data1->count);
	ptrdiff_t k=0;
	const unsigned char *ptr0=data0->data, *ptr1=data1->data;
	for(;k<size-63;k+=64)
	{
		__m256i a0=_mm256_loadu_si256((__m256i*)ptr0+0);
		__m256i a1=_mm256_loadu_si256((__m256i*)ptr0+1);
		__m256i b0=_mm256_loadu_si256((__m256i*)ptr1+0);
		__m256i b1=_mm256_loadu_si256((__m256i*)ptr1+1);
		a0=_mm256_cmpeq_epi64(a0, b0);
		a1=_mm256_cmpeq_epi64(a1, b1);
		int mask0=_mm256_movemask_pd(_mm256_castsi256_pd(a0));
		int mask1=_mm256_movemask_pd(_mm256_castsi256_pd(a1));
		if((mask0|mask1<<4)!=255)
		{
			broken=1;
			printf("  verify_files [%td]:\n", k);
			printf("Original:  ");
			for(int k2=0;k2<64;++k2)
				printf(" %02X", ptr0[k2]);
			printf("\n");
			printf("Corrupt:   ");
			for(int k2=0;k2<64;++k2)
				printf(" %02X", ptr1[k2]);
			printf("\n");
			printf("XOR:       ");
			for(int k2=0;k2<64;++k2)
			{
				int val=ptr0[k2]^ptr1[k2];
				if(val)
					printf(" %02X", val);
				else
					printf(" --");
			}
			printf("\n");
			break;
		}
		ptr0+=64;
		ptr1+=64;
	}
	if(!broken)
	{
		for(;k<size;++k)
		{
			if(*ptr1!=*ptr0)
			{
				broken=1;
				printf("  verify_files [%td]: Original %02X vs Corrupt %02X  XOR %02X\n", k, *ptr0, *ptr1, *ptr0^*ptr1);
				break;
			}
			++ptr0;
			++ptr1;
		}
	}
	if(broken)
		LOG_ERROR("");
	array_free(&data0);
	array_free(&data1);
}
static unsigned char* load_ppm(const char *fn, int *ret_iw, int *ret_ih)
{
	int iw=0, ih=0;
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		LOG_ERROR("Cannot open \"%s\"", fn);
		return 0;
	}
	int c=0;
	ptrdiff_t nread=fread(&c, 1, 2, fsrc);
	if(nread!=2||c!=('P'|'6'<<8))
	{
		LOG_ERROR("Inalid file \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		LOG_ERROR("Invalid PPM file");
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
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *image=(unsigned char*)malloc(size);
	if(!image)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	nread=fread(image, 1, size, fsrc);
	if(nread!=size)
	{
		LOG_ERROR("Truncated file  expected %td  got %td", size, nread);
		free(image);
		return 0;
	}
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	return image;
}
#define SSIM_SIZE 11
static void measure_ssim_ppm_avx2(const char *fn0, const char *fn1, long long csize, double *ret_ssim, double *ret_weight)//buggy
{
	const int psize=(int)sizeof(int[3][SSIM_SIZE*SSIM_SIZE][8]);
	int
		*patch1=(int*)_mm_malloc(psize, sizeof(__m256i)),
		*patch2=(int*)_mm_malloc(psize, sizeof(__m256i));
	if(!patch1||!patch2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	int iw=0, ih=0;
	unsigned char *im1=0, *im2=0;
	{
		int w2=0, h2=0;
		im1=load_ppm(fn0, &iw, &ih);
		im2=load_ppm(fn1, &w2, &h2);
		if(!im1||!im2||iw!=w2||ih!=h2)
		{
			LOG_ERROR("Dimension mismatch %d*%d vs %d*%d", iw, ih, w2, h2);
			return;
		}
	}
	//printf("WH %d*%d\n", iw, ih);
	ptrdiff_t count=0;
	__m256 mssim[4]={0};
	__m256 c1=_mm256_set1_ps(0.01f*255), c2=_mm256_set1_ps(0.03f*255);
	c1=_mm256_mul_ps(c1, c1);
	c2=_mm256_mul_ps(c2, c2);
	for(int ky=SSIM_SIZE/2;ky<ih-SSIM_SIZE/2-1;ky+=2)
	{
		for(int kx=SSIM_SIZE/2;kx<iw-SSIM_SIZE/2-3;kx+=4)
		{
			int idx=3*(iw*ky+kx);
			unsigned char *src1=im1+idx, *src2=im2+idx;
			__m256i half=_mm256_set1_epi32(128);
			for(int ky3=0;ky3<2;++ky3)
			{
				for(int kx3=0;kx3<4;++kx3)
				{
					for(int ky2=0;ky2<SSIM_SIZE;++ky2)
					{
						for(int kx2=0;kx2<SSIM_SIZE;++kx2)
						{
							int srcidx=3*(iw*(ky3+ky2-SSIM_SIZE/2)+kx3+kx2-SSIM_SIZE/2);
							int dstidx=(SSIM_SIZE*ky2+kx2)*3*8+ky3*4+kx3;
							patch1[8*0+dstidx]=src1[srcidx+0];
							patch2[8*0+dstidx]=src2[srcidx+0];
							patch1[8*1+dstidx]=src1[srcidx+1];
							patch2[8*1+dstidx]=src2[srcidx+1];
							patch1[8*2+dstidx]=src1[srcidx+2];
							patch2[8*2+dstidx]=src2[srcidx+2];
						}
					}
				}
			}
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)//RGB->YUV
			{
				__m256i r1=_mm256_load_si256((__m256i*)patch1+3*k+0);
				__m256i g1=_mm256_load_si256((__m256i*)patch1+3*k+1);
				__m256i b1=_mm256_load_si256((__m256i*)patch1+3*k+2);
				__m256i r2=_mm256_load_si256((__m256i*)patch2+3*k+0);
				__m256i g2=_mm256_load_si256((__m256i*)patch2+3*k+1);
				__m256i b2=_mm256_load_si256((__m256i*)patch2+3*k+2);
				r1=_mm256_sub_epi32(r1, half);
				g1=_mm256_sub_epi32(g1, half);
				b1=_mm256_sub_epi32(b1, half);
				r2=_mm256_sub_epi32(r2, half);
				g2=_mm256_sub_epi32(g2, half);
				b2=_mm256_sub_epi32(b2, half);
				r1=_mm256_sub_epi32(r1, g1);
				b1=_mm256_sub_epi32(b1, g1);
				g1=_mm256_add_epi32(g1, _mm256_srai_epi32(_mm256_add_epi32(r1, g1), 2));
				r2=_mm256_sub_epi32(r2, g2);
				b2=_mm256_sub_epi32(b2, g2);
				g2=_mm256_add_epi32(g2, _mm256_srai_epi32(_mm256_add_epi32(r2, g2), 2));
				_mm256_store_ps((float*)patch1+8*(3*k+0), _mm256_cvtepi32_ps(r1));
				_mm256_store_ps((float*)patch1+8*(3*k+1), _mm256_cvtepi32_ps(g1));
				_mm256_store_ps((float*)patch1+8*(3*k+2), _mm256_cvtepi32_ps(b1));
				_mm256_store_ps((float*)patch2+8*(3*k+0), _mm256_cvtepi32_ps(r2));
				_mm256_store_ps((float*)patch2+8*(3*k+1), _mm256_cvtepi32_ps(g2));
				_mm256_store_ps((float*)patch2+8*(3*k+2), _mm256_cvtepi32_ps(b2));
			}
			__m256 mean1[3], mean2[3];
			memset(mean1, 0, sizeof(mean1));
			memset(mean2, 0, sizeof(mean2));
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
			{
				mean1[0]=_mm256_add_ps(mean1[0], _mm256_load_ps((float*)patch1+8*(3*k+0)));
				mean1[1]=_mm256_add_ps(mean1[1], _mm256_load_ps((float*)patch1+8*(3*k+1)));
				mean1[2]=_mm256_add_ps(mean1[2], _mm256_load_ps((float*)patch1+8*(3*k+2)));
				mean2[0]=_mm256_add_ps(mean2[0], _mm256_load_ps((float*)patch2+8*(3*k+0)));
				mean2[1]=_mm256_add_ps(mean2[1], _mm256_load_ps((float*)patch2+8*(3*k+1)));
				mean2[2]=_mm256_add_ps(mean2[2], _mm256_load_ps((float*)patch2+8*(3*k+2)));
			}
			__m256 mcount=_mm256_set1_ps(1.f/(SSIM_SIZE*SSIM_SIZE));
			mean1[0]=_mm256_mul_ps(mean1[0], mcount);
			mean1[1]=_mm256_mul_ps(mean1[1], mcount);
			mean1[2]=_mm256_mul_ps(mean1[2], mcount);
			mean2[0]=_mm256_mul_ps(mean2[0], mcount);
			mean2[1]=_mm256_mul_ps(mean2[1], mcount);
			mean2[2]=_mm256_mul_ps(mean2[2], mcount);
			__m256 var1[3], var2[3], cov[3];
			memset(var1, 0, sizeof(var1));
			memset(var2, 0, sizeof(var2));
			memset(cov, 0, sizeof(cov));
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
			{
				__m256 d1[]=
				{
					_mm256_sub_ps(mean1[0], _mm256_load_ps((float*)patch1+8*(3*k+0))),
					_mm256_sub_ps(mean1[1], _mm256_load_ps((float*)patch1+8*(3*k+1))),
					_mm256_sub_ps(mean1[2], _mm256_load_ps((float*)patch1+8*(3*k+2))),
				};
				__m256 d2[]=
				{
					_mm256_sub_ps(mean2[0], _mm256_load_ps((float*)patch2+8*(3*k+0))),
					_mm256_sub_ps(mean2[1], _mm256_load_ps((float*)patch2+8*(3*k+1))),
					_mm256_sub_ps(mean2[2], _mm256_load_ps((float*)patch2+8*(3*k+2))),
				};
				var1[0]=_mm256_add_ps(var1[0], _mm256_mul_ps(d1[0], d1[0]));
				var2[0]=_mm256_add_ps(var2[0], _mm256_mul_ps(d2[0], d2[0]));
				cov[0]=_mm256_add_ps(cov[0], _mm256_mul_ps(d1[0], d2[0]));
			}
			var1[0]=_mm256_mul_ps(var1[0], mcount);
			var1[1]=_mm256_mul_ps(var1[1], mcount);
			var1[2]=_mm256_mul_ps(var1[2], mcount);
			var2[0]=_mm256_mul_ps(var2[0], mcount);
			var2[1]=_mm256_mul_ps(var2[1], mcount);
			var2[2]=_mm256_mul_ps(var2[2], mcount);
			cov[0]=_mm256_mul_ps(cov[0], mcount);
			cov[1]=_mm256_mul_ps(cov[1], mcount);
			cov[2]=_mm256_mul_ps(cov[2], mcount);
			__m256 num1[]=
			{
				_mm256_mul_ps(mean1[0], mean2[0]),
				_mm256_mul_ps(mean1[1], mean2[1]),
				_mm256_mul_ps(mean1[2], mean2[2]),
			};
			__m256 num2[]=
			{
				_mm256_add_ps(cov[0], cov[0]),
				_mm256_add_ps(cov[1], cov[1]),
				_mm256_add_ps(cov[2], cov[2]),
			};
			num1[0]=_mm256_add_ps(num1[0], num1[0]);
			num1[1]=_mm256_add_ps(num1[1], num1[1]);
			num1[2]=_mm256_add_ps(num1[2], num1[2]);
			num1[0]=_mm256_add_ps(num1[0], c1);
			num1[1]=_mm256_add_ps(num1[1], c1);
			num1[2]=_mm256_add_ps(num1[2], c1);
			num2[0]=_mm256_add_ps(num2[0], c2);
			num2[1]=_mm256_add_ps(num2[1], c2);
			num2[2]=_mm256_add_ps(num2[2], c2);
			num1[0]=_mm256_mul_ps(num1[0], num2[0]);
			num1[1]=_mm256_mul_ps(num1[1], num2[1]);
			num1[2]=_mm256_mul_ps(num1[2], num2[2]);
			__m256 den1[]=
			{
				_mm256_add_ps(_mm256_mul_ps(mean1[0], mean1[0]), _mm256_mul_ps(mean2[0], mean2[0])),
				_mm256_add_ps(_mm256_mul_ps(mean1[1], mean1[1]), _mm256_mul_ps(mean2[1], mean2[1])),
				_mm256_add_ps(_mm256_mul_ps(mean1[2], mean1[2]), _mm256_mul_ps(mean2[2], mean2[2])),
			};
			den1[0]=_mm256_add_ps(den1[0], c1);
			den1[1]=_mm256_add_ps(den1[1], c1);
			den1[2]=_mm256_add_ps(den1[2], c1);
			__m256 den2[]=
			{
				_mm256_add_ps(_mm256_add_ps(var1[0], var2[0]), c2),
				_mm256_add_ps(_mm256_add_ps(var1[1], var2[1]), c2),
				_mm256_add_ps(_mm256_add_ps(var1[2], var2[2]), c2),
			};
			den1[0]=_mm256_mul_ps(den1[0], den2[0]);
			den1[1]=_mm256_mul_ps(den1[1], den2[1]);
			den1[2]=_mm256_mul_ps(den1[2], den2[2]);
			num1[0]=_mm256_div_ps(num1[0], den1[0]);
			num1[1]=_mm256_div_ps(num1[1], den1[1]);
			num1[2]=_mm256_div_ps(num1[2], den1[2]);
			mssim[0]=_mm256_add_ps(mssim[0], num1[0]);
			mssim[1]=_mm256_add_ps(mssim[1], num1[1]);
			mssim[2]=_mm256_add_ps(mssim[2], num1[2]);
			//double curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));
			//ssim[kc]+=curr_ssim;
			++count;
		}
	}
	free(im1);
	free(im2);
	_mm_free(patch1);
	_mm_free(patch2);
	ALIGN(32) float ssim[4][8]={0};
	if(count)
	{
		__m256 norm=_mm256_set1_ps(1.f/(8*count));
		mssim[0]=_mm256_mul_ps(mssim[0], norm);
		mssim[1]=_mm256_mul_ps(mssim[1], norm);
		mssim[2]=_mm256_mul_ps(mssim[2], norm);
		_mm256_store_ps(ssim[0], mssim[0]);
		_mm256_store_ps(ssim[1], mssim[1]);
		_mm256_store_ps(ssim[2], mssim[2]);
		for(int k=1;k<8;++k)
		{
			ssim[0][0]+=ssim[0][k];
			ssim[1][0]+=ssim[1][k];
			ssim[2][0]+=ssim[2][k];
		}
	}
	ssim[3][0]=(6*ssim[0][0]+ssim[1][0]+ssim[2][0])/8;
	long long res=(long long)iw*ih;
	if(res)
		printf("  BPD%7.4lf SSIM%9.6lf%9.6lf%9.6lf%9.6lf"
			, 8.*csize/(3*res)
			, ssim[3][0]
			, ssim[0][0]
			, ssim[1][0]
			, ssim[2][0]
		);
	double w=(double)res/(1024*1024);
	ret_ssim[0]+=ssim[0][0]*w;
	ret_ssim[1]+=ssim[1][0]*w;
	ret_ssim[2]+=ssim[2][0]*w;
	ret_ssim[3]+=ssim[3][0]*w;
	*ret_weight+=w;
}
static void measure_ssim_ppm(const char *fn0, const char *fn1, long long csize, double *ret_ssim, double *ret_weight)
{
	int iw=0, ih=0;
	unsigned char *im1=0, *im2=0;
	{
		int w2=0, h2=0;
		im1=load_ppm(fn0, &iw, &ih);
		im2=load_ppm(fn1, &w2, &h2);
		if(!im1||!im2||iw!=w2||ih!=h2)
		{
			LOG_ERROR("Dimension mismatch %d*%d vs %d*%d", iw, ih, w2, h2);
			return;
		}
	}
	//printf("WH %d*%d\n", iw, ih);
	ptrdiff_t count=0;
	double ssim[4]={0};
	double c1=0.01*255, c2=0.03*255;
	c1*=c1;
	c2*=c2;
	for(int ky=SSIM_SIZE/2;ky<ih-SSIM_SIZE/2;ky+=5)
	{
		for(int kx=SSIM_SIZE/2;kx<iw-SSIM_SIZE/2;kx+=5)
		{
			int idx=3*(iw*ky+kx);
			unsigned char *p1=im1+idx, *p2=im2+idx;
			short patch1[3][SSIM_SIZE*SSIM_SIZE], patch2[3][SSIM_SIZE*SSIM_SIZE];
			for(int ky2=0;ky2<SSIM_SIZE;++ky2)
			{
				for(int kx2=0;kx2<SSIM_SIZE;++kx2)
				{
					int srcidx=3*(iw*(ky2-SSIM_SIZE/2)+kx2-SSIM_SIZE/2);
					int dstidx=SSIM_SIZE*ky2+kx2;
					patch1[0][dstidx]=p1[srcidx+0]-128;
					patch2[0][dstidx]=p2[srcidx+0]-128;
					patch1[1][dstidx]=p1[srcidx+1]-128;
					patch2[1][dstidx]=p2[srcidx+1]-128;
					patch1[2][dstidx]=p1[srcidx+2]-128;
					patch2[2][dstidx]=p2[srcidx+2]-128;
				}
			}
			for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)//RGB->YUV
			{
				int r1=patch1[0][k], g1=patch1[1][k], b1=patch1[2][k];
				int r2=patch2[0][k], g2=patch2[1][k], b2=patch2[2][k];

				r1-=g1;
				b1-=g1;
				g1+=(r1+b1)>>2;
				r2-=g2;
				b2-=g2;
				g2+=(r2+b2)>>2;

				patch1[0][k]=g1;
				patch1[1][k]=b1;
				patch1[2][k]=r1;
				patch2[0][k]=g2;
				patch2[1][k]=b2;
				patch2[2][k]=r2;
			}
			for(int kc=0;kc<3;++kc)
			{
				double mean1=0, mean2=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[kc][k];
					mean2+=patch2[kc][k];
				}
				mean1/=SSIM_SIZE*SSIM_SIZE;
				mean2/=SSIM_SIZE*SSIM_SIZE;
				double var1=0, var2=0, cov=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					double d1=patch1[kc][k]-mean1;
					double d2=patch2[kc][k]-mean2;
					var1+=d1*d1;
					var2+=d2*d2;
					cov+=d1*d2;
				}
				var1/=SSIM_SIZE*SSIM_SIZE;
				var2/=SSIM_SIZE*SSIM_SIZE;
				cov/=SSIM_SIZE*SSIM_SIZE;
				double curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));
				ssim[kc]+=curr_ssim;
			}
			++count;
		}
	}
	free(im1);
	free(im2);
	if(count)
	{
		ssim[0]/=count;
		ssim[1]/=count;
		ssim[2]/=count;
	}
	ssim[3]=(6*ssim[0]+ssim[1]+ssim[2])/8;
	long long res=(long long)iw*ih;
	if(res)
		printf("  BPD%7.4lf SSIM%9.6lf%9.6lf%9.6lf%9.6lf", 8.*csize/(3*res), ssim[3], ssim[0], ssim[1], ssim[2]);
	double w=(double)res/(1024*1024);
	ret_ssim[0]+=ssim[0]*w;
	ret_ssim[1]+=ssim[1]*w;
	ret_ssim[2]+=ssim[2]*w;
	ret_ssim[3]+=ssim[3]*w;
	*ret_weight+=w;
}
#if 0
static void measure_ssim_ppm(const char *fn0, const char *fn1, long long csize, double *ret_ssim, long long *ret_bmpsize)
{
	int iw=0, ih=0;
	unsigned char *im1=0, *im2=0;
	{
		int w2=0, h2=0;
		im1=load_ppm(fn0, &iw, &ih);
		im2=load_ppm(fn1, &w2, &h2);
		if(iw!=w2||ih!=h2)
		{
			LOG_ERROR("Dimension mismatch %d*%d vs %d*%d", iw, ih, w2, h2);
			return;
		}
	}
	ptrdiff_t count=0;
	double ssim=0;
	double c1=0.01*255, c2=0.03*255;
	c1*=c1;
	c2*=c2;
	for(int ky=SSIM_SIZE/2;ky<ih-SSIM_SIZE/2;++ky)
	{
		for(int kx=SSIM_SIZE/2;kx<iw-SSIM_SIZE/2;++kx)
		{
			int idx=3*(iw*ky+kx);
			unsigned char *p1=im1+idx, *p2=im2+idx;
			for(int kc=0;kc<3;++kc)
			{
				unsigned char patch1[SSIM_SIZE*SSIM_SIZE], patch2[SSIM_SIZE*SSIM_SIZE];
				for(int ky2=0;ky2<SSIM_SIZE;++ky2)
				{
					for(int kx2=0;kx2<SSIM_SIZE;++kx2)
					{
						int idx=3*(iw*(ky2-SSIM_SIZE/2)+kx-SSIM_SIZE/2)+kx;
						patch1[SSIM_SIZE*ky2+kx2]=p1[idx];
						patch2[SSIM_SIZE*ky2+kx2]=p2[idx];
					}
				}
				double mean1=0, mean2=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					mean1+=patch1[k];
					mean2+=patch2[k];
				}
				double var1=0, var2=0, cov=0;
				for(int k=0;k<SSIM_SIZE*SSIM_SIZE;++k)
				{
					double d1=patch1[k]-mean1;
					double d2=patch2[k]-mean2;
					var1+=d1*d1;
					var2+=d2*d2;
					cov+=d1*d2;
				}
				double curr_ssim=(2*mean1*mean2+c1)*(2*cov+c2)/((mean1*mean1+mean2*mean2+c1)*(var1+var2+c2));
				ssim+=curr_ssim;
				++count;
			}
		}
	}
	if(count)
		ssim/=count;
	long long bmpsize=(long long)3*iw*ih;
	if(bmpsize)
		printf("  BPD%7.4lf SSIM%9.6lf", 8.*csize/bmpsize, ssim);
	*ret_bmpsize+=bmpsize;
	*ret_ssim+=ssim*bmpsize;
	free(im1);
	free(im2);
}
#endif
static void measure_psnr_8bit(const char *fn0, const char *fn1, long long csize, double *ret_sqe)
{
	ArrayHandle data0=load_file(fn0, 1, 32, 1);
	ArrayHandle data1=load_file(fn1, 1, 32, 1);
	ptrdiff_t size=MINVAR(data0->count, data1->count), k=0;
	const unsigned char *ptr0=data0->data, *ptr1=data1->data;
	double error=0;
	for(;k<size;++k)
	{
		int diff=*ptr0++-*ptr1++;
		error+=(double)diff*diff;
	}
	if(ret_sqe)*ret_sqe+=error;
	double rmse=sqrt(error/size), psnr=20*log10(PSNR_MAX_8/rmse);
	printf("  BPD%7.4lf rmse %12.6lf psnr %12.6lf", 8.*csize/size, rmse, psnr);
//	printf(" rmse %12.6lf psnr %12.6lf", rmse, psnr);
	if(data0->count!=data1->count)
		printf("different sizes %zd vs %zd:\n", data0->count, data1->count);
	array_free(&data0);
	array_free(&data1);
}
static void measure_psnr_16bit(const char *fn0, const char *fn1, long long csize, double *ret_sqe)
{
	ArrayHandle data0=load_file(fn0, 1, 32, 1);
	ArrayHandle data1=load_file(fn1, 1, 32, 1);
	ptrdiff_t size=MINVAR(data0->count, data1->count)>>1, k=0;
	const unsigned short *ptr0=(const unsigned short*)data0->data, *ptr1=(const unsigned short*)data1->data;
	double error=0;
	for(;k<size;++k)
	{
		int diff=(short)(*ptr0++-*ptr1++);
		error+=(double)diff*diff;
	}
	if(ret_sqe)*ret_sqe+=error;
	double rmse=sqrt(error/size), psnr=20*log10(PSNR_MAX_16/rmse);
	printf("  rmse %12.6lf psnr %12.6lf", rmse, psnr);
	if(data0->count!=data1->count)
		printf("different sizes %zd vs %zd:\n", data0->count, data1->count);
	array_free(&data0);
	array_free(&data1);
}

//GDCC score
#define CALCSCORE(CSIZE, ENC, DEC) ((CSIZE)*(1./(1024*1024))+(ENC)+(DEC)*2)

typedef struct _UInfo
{
	long long usize;
	ArrayHandle filename;
} UInfo;
static void free_uinfo(void *p)
{
	UInfo *info=(UInfo*)p;
	array_free(&info->filename);
}
typedef struct _CellInfo
{
	long long csize;
	double etime, dtime;
	long long emem, dmem;
} CellInfo;
typedef struct _TestInfo
{
	ArrayHandle codecname;
	DateTime datetime;
	CellInfo total;
	ArrayHandle cells;
} TestInfo;
static void free_testinfo(void *p)
{
	TestInfo *info=(TestInfo*)p;
	array_free(&info->codecname);
	array_free(&info->cells);
}
typedef struct _Range
{
	int issrc, start, end;
} Range;
typedef struct _CommandFormat
{
	ArrayHandle format;
	ArrayHandle bounds;//Range
} CommandFormat;
static int getlineno(const char *start, const char *ptr)
{
	int lineno=1;
	while(ptr<start)
		lineno+=*ptr++=='\n';
	return lineno;
}
static void skipspace(const char **ptr)
{
	const char *ptr2=*ptr;
	while(*ptr2&&isspace(*ptr2))++ptr2;
	*ptr=ptr2;
}
static void skip2space(const char **ptr)
{
	const char *ptr2=*ptr;
	while(*ptr2&&!isspace(*ptr2))++ptr2;
	*ptr=ptr2;
}
static int peeklabel(const char *ptr, const char *label)
{
	const char *lptr=label;
	while(*lptr&&*ptr++==*lptr++);
	return !*lptr;
}
static void skiplabel(const char *start, const char **ptr, const char *label)
{
	if(!*ptr)
	{
		LOG_ERROR("");
		return;
	}
	const char *ptr2=*ptr, *lptr=label;
	while(*lptr&&*ptr2++==*lptr++);
	if(*lptr)
	{
		int lineno=getlineno(start, *ptr);
		int len=(int)strlen(*ptr);
		if(len>(int)(lptr-label))
			len=(int)(lptr-label);
		printf("Parser line %d expected \"%s\" got \"%.*s\"", lineno, label, len, *ptr);
		LOG_ERROR("");
	}
	*ptr=ptr2;
}
static ArrayHandle parse_str(const char *start, const char **ptr, char delim, int allow_eof)//delim must be one of the whitespace characters
{
	ArrayHandle str=0;
	const char *ptr2=*ptr;
	for(;*ptr2&&*ptr2!=delim;++ptr2);
	if(allow_eof&&!*ptr2)
		goto skipcheck;
	if(*ptr2!=delim)
	{
		int lineno=getlineno(start, *ptr);
		int len=(int)strlen(*ptr);
		if(len>10)
			len=10;
		printf("Parser line %d expected \"%c\"=0x%02X got \"%.*s\"", lineno, delim, delim, len, *ptr);
		LOG_ERROR("");
	}
skipcheck:
	{
		const char *ptr3=ptr2;
		while(ptr3>*ptr&&isspace(*(ptr3-1)))--ptr3;
		STR_COPY(str, *ptr, ptr3-*ptr);
	}
	//skipspace(&ptr2);
	//while(*ptr2&&isspace(*ptr2++));
	//ptr2-=*ptr2!=0;
	*ptr=ptr2;
	return str;
}
static long long parse_uint(const char **ptr)
{
	const char *ptr2=*ptr;
	long long val=0;
	skipspace(&ptr2);
	while((unsigned)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	*ptr=ptr2;
	return val;
}
static double parse_float(const char **ptr)
{
	const char *ptr2=*ptr;
	double val=0;
	skipspace(&ptr2);
	while((unsigned)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	if(*ptr2=='.')
	{
		++ptr2;
		double p=0.1;
		while((unsigned)(*ptr2-'0')<10)
		{
			val+=(*ptr2++-'0')*p;
			p*=0.1;
		}
	}
	*ptr=ptr2;
	return val;
}
static void parse_cell(const char **ptr, CellInfo *cell)
{
	//{csize  etime  dtime  emem  dmem}
	cell->csize=parse_uint(ptr);
	cell->etime=parse_float(ptr);
	cell->dtime=parse_float(ptr);
	cell->emem=parse_uint(ptr);
	cell->dmem=parse_uint(ptr);
}
static int strimatch(const char *text, const char *label)//return 1: match
{
	while(*label&&tolower(*text)==tolower(*label))++text, ++label;
	return !*label;
}
static ArrayHandle parse_cmdformatext(const char *text, int len, const char *srcext, int requiredst, int *srcdstbounds)
{
	ArrayHandle bounds=0;
	const char *ptr=text, *end=text+len-1;
	memset(srcdstbounds, -1, sizeof(int[4]));
	while(ptr<end)
	{
		short c=*(short*)ptr;
		if(c==('*'|'.'<<8))
		{
			Range *range=(Range*)array_append(&bounds, 0, sizeof(Range), 1, 1, 0, 0);
			int idx=(int)(ptr-text);
			range->issrc=strimatch(text+idx+2, srcext);//add 2 to skip "*."
			range->start=idx;
			skip2space(&ptr);
			range->end=(int)(ptr-text);
			int *ridx=range->issrc?srcdstbounds:srcdstbounds+2;
			if(*ridx==-1)
			{
				ridx[0]=range->start;
				ridx[1]=range->end;
			}
			else if(!acme_strnimatch(
				text+range->start, (ptrdiff_t)range->end-range->start,
				text+ridx[0], (ptrdiff_t)ridx[1]-ridx[0]
			))
			{
				printf(
					"\n"
					"Invalid template:  %s\n"
					"Expected\n"
					"  [%d] \"%.*s\"  to match\n"
					"  [%d] \"%.*s\"  (case-insensitive, no spaces)\n",
					text,
					ridx[0], ridx[1]-ridx[0], text+ridx[0],
					range->start, range->end-range->start, text+range->start
				);
				LOG_ERROR("");
			}
			continue;
		}
		++ptr;
	}
	if(!bounds)
	{
		printf(
			"\n"
			"Invalid template:  %s\n"
			"Missing source%s \"*.extension%s\"\n",
			text,
			requiredst?" and codec":"",
			requiredst?"s":""
		);
		LOG_ERROR("");
	}
	if(requiredst&&srcdstbounds[2]==-1)
	{
		printf(
			"\n"
			"Invalid template:  %s\n"
			"Missing codec \"*.extension\"\n",
			text
		);
		LOG_ERROR("");
	}
	return bounds;
}
static void substitute_cmdplaceholders(char **pptr, const char *end, const CommandFormat *cmd, const char *t1fn, const char *t2fn)
{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wrestrict"
#endif
	char *ptr=*pptr;
	int lastidx=0;
	for(int k=0;k<(int)cmd->bounds->count;++k)
	{
		const Range *range=(const Range*)array_at((ArrayHandle*)(size_t)&cmd->bounds, k);
		ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" ",
			range->start-lastidx, (char*)cmd->format->data+lastidx,
			range->issrc?t1fn:t2fn
		);
		lastidx=range->end;
	}
	ptr+=snprintf(ptr, end-ptr, "%s",
		(char*)cmd->format->data+lastidx
	)+1;
	*pptr=ptr;
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
}
static const char* get_extension(const char *filename, ptrdiff_t len)//excludes the dot
{
	ptrdiff_t idx;

	idx=acme_strrchr(filename, len, '.');
	if(idx==-1)
		return 0;
	return filename+idx+1;
}
static ArrayHandle get_uinfo(const char *path, const char *ext)//extension without '.'
{
	//prepare searchpath
	ArrayHandle searchpath=filter_path(path, -1);
	STR_APPEND(searchpath, "*", 1, 1);

	WIN32_FIND_DATAA data={0};
	void *hSearch=FindFirstFileA((char*)searchpath->data, &data);//skip .
	if(hSearch==INVALID_HANDLE_VALUE)
		return 0;
	FindNextFileA(hSearch, &data);//skip ..

	ArrayHandle uinfo;
	ARRAY_ALLOC(UInfo, uinfo, 0, 0, 0, free_uinfo);
	while(FindNextFileA(hSearch, &data))
	{
		ptrdiff_t len=strlen(data.cFileName);
		const char *curr_ext=get_extension(data.cFileName, len);
		if(!(data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)&&!acme_stricmp(curr_ext, ext))
		{
			UInfo *info=(UInfo*)ARRAY_APPEND(uinfo, 0, 1, 1, 0);
			info->usize=(long long)data.nFileSizeHigh<<32|data.nFileSizeLow;
			STR_COPY(info->filename, searchpath->data, searchpath->count-1);
			STR_APPEND(info->filename, data.cFileName, len, 1);
		}
	}
	int success=FindClose(hSearch);
	if(!success)
	{
		SYSTEMERROR("FindClose");
		return 0;
	}
	array_free(&searchpath);
	return uinfo;
}
static void write_cell(FILE *fdst, CellInfo *cell)//tabs for excel
{
	fprintf(fdst, "\t%10lld\t%12.6lf\t%12.6lf\t%10lld\t%10lld", cell->csize, cell->etime, cell->dtime, cell->emem, cell->dmem);
}
static void qualify_command(ArrayHandle *cmd)
{
	const char *ptr=(char*)cmd[0]->data;
	skip2space(&ptr);//assume no spaces

	//while(*ptr&&!isspace(*ptr++));
	//ptr-=*ptr!=0;

	char fn1[MAX_PATH+1]={0}, fn2[MAX_PATH+2]={0};
	int len1=(int)(ptr-(char*)cmd[0]->data);
	memcpy(fn1, cmd[0]->data, len1);
	if(!strimatch(fn1+len1-4, ".exe"))
	{
		memcpy(fn1+len1, ".exe", 4);
		len1+=4;
	}
	int len2=SearchPathA(0, fn1, 0, sizeof(fn2)-2, fn2+1, 0);
	//int len=GetFullPathNameA(fn1, sizeof(fn2)-2, fn2+1, 0);
	if(!len2)
	{
		SYSTEMERROR("GetFullPathNameA");
		return;
	}
	ptrdiff_t size=get_filesize(fn2+1);
	if(size<FSIZE_EMPTYFILE)
	{
		printf("Program is inaccessible:\n");
		printf("  \"%s\"\n", fn1);
		printf("  \"%s\"\n", fn2+1);
		printf("  size %td bytes", size);
		LOG_ERROR("");
	}
	len2+=2;
	fn2[0]='\"';
	fn2[len2-1]='\"';
//	char *cmd2=(char*)
	array_replace(cmd, 0, ptr-(char*)cmd[0]->data, fn2, len2, 1, 1);
//	printf("  %s\n", cmd2);//
}
static void exec_process2(char *cmd, const char *currdir, int loud, double *elapsed, long long *maxmem, int *threadcount)
{
	int success;
	STARTUPINFOA si={0};
	PROCESS_INFORMATION pi={0};
	si.cb=sizeof(si);
	if(!loud)
	{
		si.dwFlags=STARTF_USESTDHANDLES;
		si.hStdOutput=CreateFileA("NUL", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdOutput==INVALID_HANDLE_VALUE)
		{
			SYSTEMERROR("CreateFileA");
			return;
		}
		si.hStdError=CreateFileA("NUL", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdError==INVALID_HANDLE_VALUE)
		{
			SYSTEMERROR("CreateFileA");
			return;
		}
		si.hStdInput=CreateFileA("NUL", GENERIC_READ, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdInput==INVALID_HANDLE_VALUE)
		{
			SYSTEMERROR("CreateFileA");
			return;
		}
	}
	success=CreateProcessA(0, cmd, 0, 0, 0, CREATE_SUSPENDED, 0, currdir, &si, &pi);
	if(!success)
	{
		SYSTEMERROR("CreateProcessA");
		return;
	}
	ptrdiff_t memusage=0;
	double t=time_sec();
	int suspendcount=ResumeThread(pi.hThread);
	if(suspendcount==(DWORD)-1)
	{
		SYSTEMERROR("CreateProcessA");
		return;
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
#if 0
#ifdef TRACK_THREADS
	int k=0;
#endif
	while(WaitForSingleObject(pi.hProcess, 10)==WAIT_TIMEOUT)
	{
		PROCESS_MEMORY_COUNTERS pmc={0};
		pmc.cb=sizeof(pmc);
		success=GetProcessMemoryInfo(pi.hProcess, &pmc, sizeof(pmc));
		if(!success)
		{
			SYSTEMERROR("GetProcessMemoryInfo");
			return;
		}
		if(memusage<(ptrdiff_t)pmc.WorkingSetSize)
			memusage=pmc.WorkingSetSize;
#ifdef TRACK_THREADS
		if(!k)
		{
			HANDLE snapshot=CreateToolhelp32Snapshot(TH32CS_SNAPALL, pi.dwProcessId);
			if(snapshot==INVALID_HANDLE_VALUE)
			{
				continue;
				//SYSTEMERROR("CreateToolhelp32Snapshot");
				//return;
			}
			THREADENTRY32 threadentry={sizeof(THREADENTRY32)};
			int nthreads=0;
			for(int k2=0;;k2|=1)
			{
				if(!k2)
				{
					if(!Thread32First(snapshot, &threadentry))
						break;
				}
				else
				{
					if(!Thread32Next(snapshot, &threadentry))
						break;
				}
				nthreads+=threadentry.th32OwnerProcessID==pi.dwProcessId;
			}
			CloseHandle(snapshot);
			if(*threadcount<nthreads)
				*threadcount=nthreads;
		}
		k|=1;
#endif
	}
#endif
	if(elapsed)*elapsed=time_sec()-t;
	FILETIME tstart={0}, tfinish={0}, tkernel={0}, tuser={0};
	//timer7
	GetProcessTimes(pi.hProcess, &tstart, &tfinish, &tkernel, &tuser);
	PROCESS_MEMORY_COUNTERS counters={sizeof(PROCESS_MEMORY_COUNTERS)};
	GetProcessMemoryInfo(pi.hProcess, &counters, sizeof(counters));
	memusage=counters.PeakWorkingSetSize;
	if(maxmem)*maxmem=memusage;
	success=CloseHandle(pi.hThread);
	if(!success)
	{
		SYSTEMERROR("CloseHandle");
		return;
	}
	success=CloseHandle(pi.hProcess);
	if(!success)
	{
		SYSTEMERROR("CloseHandle");
		return;
	}
	if(!loud)
	{
		success=CloseHandle(si.hStdOutput);
		if(!success)
		{
			SYSTEMERROR("CloseHandle");
			return;
		}
		success=CloseHandle(si.hStdError);
		if(!success)
		{
			SYSTEMERROR("CloseHandle");
			return;
		}
		success=CloseHandle(si.hStdInput);
		if(!success)
		{
			SYSTEMERROR("CloseHandle");
			return;
		}
	}
}
static int print_scicolor(double xprint, double xcolor, int is_ratio)//is_ratio = 0: diff  1: ratio
{
	int dB=0, red=0, green=0;
	/*
	diff:
	x		log10(x)	round((log10(x)+2)*256)
	200KB/8MB	-1.60		101
	300KB/10MB	-1.52		121
	400KB/10MB			178

	ratio:
	x		log10(x)	round(log10(x)*256)
	0.01x		-2		-256
	0.2x		-0.699		-89
	1x		0		0
	2x		0.301		38
	90x		1.954		250
	*/
	if(is_ratio)
	{
		dB=(int)(log10(xcolor)*256);
		CLAMP2(dB, -255, 255);
		if(xcolor<1)//small ratio is good for the rival	eg: 0.5x time -> 2x faster
			red=-dB;
		else if(xcolor>1)
			green=dB;
	}
	else if(xcolor)
	{
		dB=(int)round(log10(fabs(xcolor))*256+2*256);
		CLAMP2(dB, 0, 255);
		if(xcolor<0)//negative diff is good for the rival	eg: -300 KB
			red=dB;
		else if(xcolor>0)
			green=dB;
	}
	static char str[64]={0};
	int printed=snprintf(str, sizeof(str)-1, "%+5.0e", xprint);
	if(printed!=6)
	{
		printed=printf("\nprint_scicolor:  %s\n", str);
		LOG_ERROR("");
		return printed;
	}
	// 0123456
	//"+2e+05"
	//"+25\0"	diff
	//"2+5\0"	ratio
	static const char expdigit[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	if(is_ratio)
	{
		str[0]=str[1];
		str[1]=str[3];
	}
	str[2]=expdigit[(str[4]-'0')<<4|(str[5]-'0')];
	str[3]=0;
	return colorprintf(green>128?COLORPRINTF_TXT_DEFAULT^0xFFFFFF:COLORPRINTF_TXT_DEFAULT, green<<8|red, "%s", str);
}
static void print_rivals_v2(ArrayHandle besttestidxs, ArrayHandle testinfo, int sampleidx, const CellInfo *currcell, ptrdiff_t usize)
{
	for(int k2=0;k2<(int)besttestidxs->count;++k2)
	{
		int *idx=(int*)array_at(&besttestidxs, k2);
		TestInfo *test=(TestInfo*)array_at(&testinfo, *idx);
		CellInfo *cell=&test->total;
		if(sampleidx>=0)
			cell=(CellInfo*)array_at(&test->cells, sampleidx);
		printf(" %X", (k2+1)&15);

		double x=(double)cell->csize-currcell->csize;
		print_scicolor(x, x/usize, 0);

		x=cell->etime/currcell->etime;
		print_scicolor(x, x, 1);

		x=cell->dtime/currcell->dtime;
		print_scicolor(x, x, 1);
	}
}
static void print_summary(ArrayHandle besttestidxs, ArrayHandle testinfo, ptrdiff_t usize, int special, int printnotation)
{
	if(printnotation)
	{
		const double usize=800000;
		const double csize1=450000, etime1=2.5, dtime1=0.5;
		const double csize2=360000, etime2=1.3, dtime2=1.2;
		double score[]=
		{
			csize1-csize2,
			etime1/etime2,
			dtime1/dtime2,
		};
		printf("Table notation \"+SSE+ED+D\":    For example {%g-%g, %g/%g, %g/%g} = {%+.2lf KB, %gx, %gx} -> ",
			csize1, csize2,
			etime1, etime2,
			dtime1, dtime2,
			score[0]/1024, score[1], score[2]
		);
		print_scicolor(score[0], score[0]/usize, 0);
		print_scicolor(score[1], score[1], 1);
		print_scicolor(score[2], score[2], 1);
		printf(".\n");
		printf("  +SS      = size diff   = {SIGN MSDIGIT * 10 ^      EXP} = listsize - yoursize\n");
		printf("  E+E, D+D = time ratios = {     MSDIGIT * 10 ^ SIGN EXP} = listtime / yourtime\n");
	}
	for(int k2=0;k2<(int)besttestidxs->count;++k2)
	{
		int *idx=(int*)array_at(&besttestidxs, k2);
		TestInfo *test=(TestInfo*)array_at(&testinfo, *idx);
		if(k2==special)
			printf("\n");
		printf("%10lld B  %12.6lf %12.6lf sec  %12.6lf %12.6lf MB/s %8.2lf %8.2lf MB  ",
			test->total.csize,
			test->total.etime,
			test->total.dtime,
			usize/(test->total.etime*1024*1024),
			usize/(test->total.dtime*1024*1024),
			(double)test->total.emem/(1024*1024),
			(double)test->total.dmem/(1024*1024)
		);
		printf("%04d%02d%02d_%02d%02d%02d"
			, test->datetime.year
			, test->datetime.month
			, test->datetime.day
			, test->datetime.hour
			, test->datetime.minute
			, test->datetime.second
		);
		if((unsigned)special<(unsigned)besttestidxs->count)
		{
			if(k2==special)
				printf("  %X <- *", (k2+1)&15);
			else if(k2<special)
				printf("  %X     ", (k2+1)&15);
			else
				printf("  %X <- %X", (k2+1)&15, k2&15);
		}
		else
			printf("  %X", (k2+1)&15);
		printf(" %s\n", (char*)test->codecname->data);
		if(k2==special)
			printf("\n");
	}
}
static void ascii_deletefile(const char *fn)
{
	ptrdiff_t size=get_filesize(fn);
	if(size>=FSIZE_EMPTYFILE)
	{
		int success=DeleteFileA(fn);
		if(!success)
			SYSTEMERROR("DeleteFileA");
	}
	else
	{
		if(size==-2)
			printf("File already gone:\n");
		else
			printf("This isn't a file:\n");
		printf("  \"%s\"\n", fn);
		LOG_ERROR("");
	}
}
int main(int argc, char **argv)
{
	const char *datasetname=0, *codecname=0;
	int flags=CMDFLAG_VERIFY_BITEXACT;
//#ifndef _DEBUG
#ifdef __GNUC__
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:    %s  DATASET  CODEC  [FLAGS]\n"
			"You will be prompted to define DATASET and CODEC.\n"
			"[FLAGS] (optional):\n"
			"  [0]  Verify files (default).\n"
			"  [1]  Don't verify bit-exact decodes.\n"
			"  [2]  Measure PSNR (8-bit).\n"
			"  [3]  Measure PSNR (16-bit).\n"
			"  [4]  Measure SSIM (8-bit PPM).\n"
			"  [5]  Print rivals. It\'s recommended to zoom out the terminal. Can be added to other flags.\n"
			"Examples:\n"
			"  %s  div2k  jxl7        Verifies that the decoded files are bit-exact.\n"
			"  %s  div2k  j2k    %d    Doesn't verify files.\n"
			"  %s  div2k  c32n   %d    Measures 8-bit PSNR.\n"
			"  %s  div2k  c32n   %d    Measures 8-bit SSIM.\n"
		//	"  %s  div2k  qlic2  %d    Prints rivals.\n"
			, argv[0]
			, argv[0]
			, argv[0], CMDFLAG_NOP
			, argv[0], CMDFLAG_PSNR_8BIT
			, argv[0], CMDFLAG_SSIM_PPM
		//	, argv[0], CMDFLAG_PRINT_RIVALS
		);
		return 0;
	}
	datasetname=argv[1];
	codecname=argv[2];
	flags=argc==4?atoi(argv[3]):CMDFLAG_VERIFY_BITEXACT;
#else
	datasetname="div2k";
	codecname="c32";
#endif
	const char placeholdertag[]="*.";
	const int placeholderlen=sizeof(placeholdertag)-1;
	char programpath[MAX_PATH+1]={0};
	ArrayHandle currdir=0, tmpfn1=0, tmpfn2=0;
	ArrayHandle srcpath=0, ext=0, uinfo=0, testinfo=0;
	CommandFormat enccmd={0}, deccmd={0};
	int srcdstbounds[4]={0};
	char *srctitle=0, *dsttitle=0;

	ArrayHandle besttestidxs=0;

	//1. get program path
	{
		int len=GetModuleFileNameA(0, programpath, sizeof(programpath)-1);
		if(!len||len==sizeof(programpath)-1)
		{
			SYSTEMERROR("GetModuleFileNameA");
			return 0;
		}
		int k=len-1;
		for(;k>=0;--k)
		{
			if(programpath[k]=='/'||programpath[k]=='\\')
			{
				++k;
				break;
			}
		}
		programpath[k]=0;
	//	printf("%s\n", programpath);//
	}

	//2. get temp filenames
	{
		int len, val;
#ifdef __GNUC__
		len=GetTempPathA(sizeof(g_buf)-1, g_buf);
#elif defined _MSC_VER
		len=GetTempPath2A(sizeof(g_buf)-1, g_buf);
#endif
		if(!len)
		{
			SYSTEMERROR("GetTempPath2A");
			return 0;
		}
		STR_COPY(currdir, g_buf, len);
		STR_ALLOC(tmpfn1, MAX_PATH+1);
		val=GetTempFileNameA(g_buf, "t1_", 0, (char*)tmpfn1->data);
		if(!val)
		{
			SYSTEMERROR("GetTempFileNameA");
			return 0;
		}
		tmpfn1->count=strlen((char*)tmpfn1->data);
		STR_ALLOC(tmpfn2, MAX_PATH+1);
		val=GetTempFileNameA(g_buf, "t2_", 0, (char*)tmpfn2->data);
		if(!val)
		{
			SYSTEMERROR("GetTempFileNameA");
			return 0;
		}
		tmpfn2->count=strlen((char*)tmpfn2->data);
		ascii_deletefile((char*)tmpfn1->data);//these were generated by GetTempFileNameA()
		ascii_deletefile((char*)tmpfn2->data);
	}

	//3. get dataset & previous tests
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
	{
		ArrayHandle text=load_file(g_buf, 0, 16, 0);//null terminated
		if(text)
		{
/*
dataset_DATASET.TXT
path	extension
"files:"
{usize  filetitle}[NFILES]
"tests:"
{codec  timestamp  {ctotal  etotal  dtotal  emax  dmax}  {csize  etime  dtime  emem  dmem}[NFILES]}

codec_CODEC.TXT
enc command template
dec command template
*/
			//int dst=0;
			const char *start=(char*)text->data, *ptr=start, *end=start+text->count;
			srcpath=parse_str(start, &ptr, '\t', 0);
			skipspace(&ptr);
			ext=parse_str(start, &ptr, '\n', 0);
			skipspace(&ptr);

			ARRAY_ALLOC(UInfo, uinfo, 0, 0, 0, free_uinfo);
			skiplabel(start, &ptr, "files:");
			while(!peeklabel(ptr, "tests:"))
			{
				UInfo *info=(UInfo*)ARRAY_APPEND(uinfo, 0, 1, 1, 0);
				info->usize=parse_uint(&ptr);
				skipspace(&ptr);
				info->filename=parse_str(start, &ptr, '\n', 0);
				skipspace(&ptr);
			}
			if(!uinfo->count)
			{
				LOG_ERROR("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			//{
			//	time_t current=time(0);
			//	struct tm *local=localtime(&current);
			//	dst=local->tm_isdst;
			//}
			ARRAY_ALLOC(TestInfo, testinfo, 0, 0, 0, free_testinfo);
			skiplabel(start, &ptr, "tests:");
			skipspace(&ptr);
			while(ptr<end)
			{
				TestInfo *info=(TestInfo*)ARRAY_APPEND(testinfo, 0, 1, 1, 0);
				info->codecname=parse_str(start, &ptr, '\t', 0);

				//YYYYmmdd_HHMMSS
#if 1
				skipspace(&ptr);
				info->datetime.year=10*(10*(10*(ptr[0]-'0')+ptr[1]-'0')+ptr[2]-'0')+ptr[3]-'0';
				ptr+=4;
				info->datetime.month=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.day=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=3;//skip '_'
				info->datetime.hour=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.minute=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.second=10*(ptr[0]-'0')+ptr[1]-'0';
				info->datetime.ds=0;
				ptr+=2;
#endif
#if 0
				int d=(int)parse_uint(&ptr);
				ptr+=*ptr=='_';
				int t=(int)parse_uint(&ptr);
#ifdef _DEBUG
				int d0=d, t0=t;
				(void)d0;
				(void)t0;
#endif
				struct tm date={0};
				date.tm_sec=t%100;
				t/=100;
				date.tm_min=t%100;
				t/=100;
				date.tm_hour=t;
				date.tm_mday=d%100;
				d/=100;
				date.tm_mon=d%100-1;
				d/=100;
				date.tm_year=d-1900;
				if((unsigned)date.tm_sec>=60
					||(unsigned)date.tm_min>=60
					||(unsigned)date.tm_hour>=24
					||(unsigned)(date.tm_mday-1)>=31
					||(unsigned)date.tm_mon>=12
					||date.tm_year<2025-1900
				)
					LOG_ERROR("Invalid timestamp %04d-%02d-%02d_%02d-%02d-%02d",
						date.tm_year+1900, date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec
					);
				date.tm_isdst=dst;
				info->timestamp=mktime(&date);
#endif
				parse_cell(&ptr, &info->total);
				ARRAY_ALLOC(CellInfo, info->cells, 0, 0, 0, 0);
				for(int k=0;k<(int)uinfo->count;++k)
				{
					if(ptr>=end)
						LOG_ERROR("Unexpected EOF");
					CellInfo *cell=(CellInfo*)ARRAY_APPEND(info->cells, 0, 1, 1, 0);
					parse_cell(&ptr, cell);
				}

				skipspace(&ptr);
			}
			array_free(&text);
		}
		else
		{
			ptrdiff_t size=0;
			int len=0;
			for(;;)
			{
				printf("Define %s path:  ", datasetname);
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				srcpath=filter_path(g_buf, len);
				size=get_filesize((char*)srcpath->data);
				if(size==FSIZE_FOLDER)
					break;
				array_free(&srcpath);
			}
			for(;;)
			{
				printf("Extension:  ");
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				if(!len)
					continue;
				int valid=1;
				for(int k=0;k<len;++k)
				{
					if(g_buf[k]<=' ')
					{
						valid=0;
						break;
					}
				}
				if(valid)
					break;
			}
			STR_COPY(ext, g_buf, len);

			uinfo=get_uinfo((char*)srcpath->data, (char*)ext->data);
			if(!uinfo||!uinfo->count)
			{
				LOG_ERROR("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			ARRAY_ALLOC(TestInfo, testinfo, 0, 0, 0, free_testinfo);
		}
	}
	int titlecolwidth=(int)strlen(codecname);
	for(int k=0;k<(int)uinfo->count;++k)//get filetitle column width
	{
		UInfo *info=(UInfo*)array_at(&uinfo, k);
		int start=0, end=0;
		get_filetitle((char*)info->filename->data, (int)info->filename->count, &start, &end);
		int width=end-start;
		if(titlecolwidth<width)
			titlecolwidth=width;
	}
	ARRAY_ALLOC(int, besttestidxs, 0, 0, testinfo->count, 0);
	{
		int *idx=(int*)ARRAY_APPEND(besttestidxs, 0, 1, 1, 0);
		*idx=0;
	}
	for(int k=1;k<(int)testinfo->count;++k)//select best prev tests
	{
		TestInfo *test1=(TestInfo*)array_at(&testinfo, k);
		double score1=CALCSCORE(test1->total.csize, test1->total.etime, test1->total.dtime);
		int encountered=0, bestsofar=1;
		for(int k2=0;k2<k;++k2)//compare with previous tests
		{
			TestInfo *test2=(TestInfo*)array_at(&testinfo, k2);
			if(acme_strnimatch((char*)test2->codecname->data, test2->codecname->count, (char*)test1->codecname->data, test1->codecname->count))
			{
				double score2=CALCSCORE(test2->total.csize, test2->total.etime, test2->total.dtime);
				encountered=1;
				if(score2<score1)//new test loses
				{
					bestsofar=0;
					break;
				}
				for(int k3=0;k3<(int)besttestidxs->count;++k3)//new test surpasses previous test
				{
					int *idx=(int*)array_at(&besttestidxs, k3);
					if(*idx==k2)
					{
						array_erase(&besttestidxs, k3, 1);
						break;
					}
				}
			}
		}
		if(!encountered||bestsofar)
		{
			int *idx=(int*)ARRAY_APPEND(besttestidxs, 0, 1, 1, 0);
			*idx=k;
		}
	}
	for(int k=0;k<(int)besttestidxs->count-1;++k)//rank besttestidxs by csize (insertion sort for simplicity)
	{
		int *idx1=(int*)array_at(&besttestidxs, k);
		TestInfo *test1=(TestInfo*)array_at(&testinfo, *idx1);
		for(int k2=k+1;k2<(int)besttestidxs->count;++k2)
		{
			int *idx2=(int*)array_at(&besttestidxs, k2);
			TestInfo *test2=(TestInfo*)array_at(&testinfo, *idx2);
			if(test1->total.csize>test2->total.csize)//test2 is tighter
			{
				int temp;
				SWAPVAR(*idx1, *idx2, temp);
				test1=test2;
			}
		}
	}

	//4. get command templates
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzcode_%s.txt", programpath, codecname);
	{
		ArrayHandle text=load_file(g_buf, 0, 16, 0);//null terminated
		ArrayHandle enc0=0, dec0=0;
		int is_new=!text;
		if(text)
		{
			const char *start=(char*)text->data, *ptr=start;
			enccmd.format=parse_str(start, &ptr, '\n', 0);
			skipspace(&ptr);
			deccmd.format=parse_str(start, &ptr, '\n', 1);
		}
		else//new
		{
			int len=0;
			printf("\n");
			printf("  Use \"*.extension\" as filename placeholders.\n");
			printf("Define %s encode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
			STR_COPY(enccmd.format, g_buf, len);
			printf("Define %s decode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
			STR_COPY(deccmd.format, g_buf, len);
			STR_COPY(enc0, enccmd.format->data, enccmd.format->count);
			STR_COPY(dec0, deccmd.format->data, deccmd.format->count);
		}
		qualify_command(&enccmd.format);
		qualify_command(&deccmd.format);
		int srcdstbounds2[4]={0};
		enccmd.bounds=parse_cmdformatext((char*)enccmd.format->data, (int)enccmd.format->count, (char*)ext->data, 1, srcdstbounds);
		deccmd.bounds=parse_cmdformatext((char*)deccmd.format->data, (int)deccmd.format->count, (char*)ext->data, 0, srcdstbounds2);
		int valid=acme_strnimatch(
			(char*)enccmd.format->data+srcdstbounds[2], (ptrdiff_t)srcdstbounds[3]-srcdstbounds[2],
			(char*)deccmd.format->data+srcdstbounds2[2], (ptrdiff_t)srcdstbounds2[3]-srcdstbounds2[2]
		);
		if(srcdstbounds2[0]!=-1)
		{
			valid+=acme_strnimatch(
				(char*)enccmd.format->data+srcdstbounds[0], (ptrdiff_t)srcdstbounds[1]-srcdstbounds[0],
				(char*)deccmd.format->data+srcdstbounds2[0], (ptrdiff_t)srcdstbounds2[1]-srcdstbounds2[0]
			)*2;
		}
		//int srcbounds[2]={-1}, dstbounds[2]={-1};
		//for(int k=0;k<(int)enccmd.bounds->count;++k)
		//{
		//	const Range *range=(const Range*)array_at(&enccmd.bounds, k);
		//	int *dst=range->issrc?srcbounds:dstbounds;
		//	if(dst[0]==-1)
		//	{
		//		dst[0]=range->start;
		//		dst[1]=range->end;
		//	}
		//	else if(!acme_strnimatch(
		//		(char*)enccmd.format->data+range->start, (ptrdiff_t)range->end-range->start,
		//		(char*)enccmd.format->data+dst[0], (ptrdiff_t)dst[1]-dst[0]
		//	))
		//	{
		//		valid=0;
		//		break;
		//	}
		//}
		//if(!enccmd.bounds||!deccmd.bounds||!valid)
		if(!(valid&1)||(srcdstbounds2[0]!=-1&&!(valid&2)))//dec dst can be omitted, eg 7zip
		{
			printf("Template \"*.extension\" mismatch\n");
			printf("Enc:  %s\n", (char*)enccmd.format->data);
			printf("Dec:  %s\n", (char*)deccmd.format->data);
			printf("\n");
			printf("Enc template requires source and destination placeholders\n");
			printf("Dec template requires source placeholder at least\n");
			printf("Source      extensions must match.\n");
			printf("Destination extensions must match.\n");
			if(!is_new)
				printf("\nCheck \"%s\"\n", g_buf);
			LOG_ERROR("");
		}
		if(is_new)
		{
			int printed=snprintf(g_buf, sizeof(g_buf)-1, "%s\n%s\n", enc0->data, dec0->data);
			char fn[128]={0};
			snprintf(fn, sizeof(fn)-1, "%szzzcode_%s.txt", programpath, codecname);
			save_file(fn, g_buf, printed, 0);
			array_free(&enc0);
			array_free(&dec0);
		}
		else
			array_free(&text);
	}
	char *t1fn=0, *t2fn, *encline=0, *decline=0;
	//int t1len=0, t2len=0, enclen=0, declen=0;
	{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wrestrict"
#endif
		char *ptr=g_buf2, *end=g_buf2+sizeof(g_buf2)-1;
		t1fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn1->data, srcdstbounds[1]-(srcdstbounds[0]+placeholderlen), (char*)enccmd.format->data+srcdstbounds[0]+placeholderlen
		)+1;//skip null terminator
		srctitle=ptr;
		while(srctitle>t1fn&&srctitle[-1]!='/'&&srctitle[-1]!='\\')--srctitle;

		t2fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn2->data, srcdstbounds[3]-(srcdstbounds[2]+placeholderlen), (char*)enccmd.format->data+srcdstbounds[2]+placeholderlen
		)+1;
		dsttitle=ptr;
		while(dsttitle>t2fn&&dsttitle[-1]!='/'&&dsttitle[-1]!='\\')--dsttitle;
		
		encline=ptr;
		substitute_cmdplaceholders(&ptr, end, &enccmd, srctitle, dsttitle);
		decline=ptr;
		substitute_cmdplaceholders(&ptr, end, &deccmd, srctitle, dsttitle);
#if 0
		char *ptr=g_buf2, *end=g_buf2+sizeof(g_buf2)-1;
		t1fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn1->data, enccmd.srcbounds[1]-(enccmd.srcbounds[0]+3), (char*)enccmd.format->data+enccmd.srcbounds[0]+3
		)+1;//skip null terminator
		t1len=(int)(ptr-t1fn-1);

		t2fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn2->data, enccmd.dstbounds[1]-(enccmd.dstbounds[0]+3), (char*)enccmd.format->data+enccmd.dstbounds[0]+3
		)+1;
		t2len=(int)(ptr-t2fn-1);

		encline=ptr;
		if(enccmd.srcbounds[0]<enccmd.dstbounds[0])
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.srcbounds[0], (char*)enccmd.format->data,
				t1fn,
				enccmd.dstbounds[0]-enccmd.srcbounds[1], (char*)enccmd.format->data+enccmd.srcbounds[1],
				t2fn,
				(char*)enccmd.format->data+enccmd.dstbounds[1]
			)+1;
		}
		else
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.dstbounds[0], (char*)enccmd.format->data,
				t2fn,
				enccmd.srcbounds[0]-enccmd.dstbounds[1], (char*)enccmd.format->data+enccmd.dstbounds[1],
				t1fn,
				(char*)enccmd.format->data+enccmd.srcbounds[1]
			)+1;
		}
		enclen=(int)(ptr-encline-1);
		
		decline=ptr;
		if(deccmd.dstbounds[0]==-1)
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %s",
				deccmd.srcbounds[0], (char*)deccmd.format->data,
				t2fn,
				(char*)deccmd.format->data+deccmd.srcbounds[1]
			)+1;
		}
		else if(deccmd.srcbounds[0]<deccmd.dstbounds[0])
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.srcbounds[0], (char*)deccmd.format->data,
				t2fn,
				deccmd.dstbounds[0]-deccmd.srcbounds[1], (char*)deccmd.format->data+deccmd.srcbounds[1],
				t1fn,
				(char*)deccmd.format->data+deccmd.dstbounds[1]
			)+1;
		}
		else
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.dstbounds[0], (char*)deccmd.format->data,
				t1fn,
				deccmd.srcbounds[0]-deccmd.dstbounds[1], (char*)deccmd.format->data+deccmd.dstbounds[1],
				t2fn,
				(char*)deccmd.format->data+deccmd.srcbounds[1]
			)+1;
		}
		declen=(int)(ptr-decline-1);
		(void)t1len;
		(void)t2len;
		(void)enclen;
		(void)declen;
#endif
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
		printf("Temp filenames:\n");//
		printf("  \"%s\"\n", t1fn);
		printf("  \"%s\"\n", t2fn);
		printf("Wirking directory:\n");
		printf("  \"%s\"\n", (char*)currdir->data);
		printf("Commands:\n");
		printf("  %s\n", encline);//
		printf("  %s\n", decline);//
		printf("\n");
		if(ptr>end)
		{
			printf("\n\nsnprintf OOB  ptr %016zX > %016zX\n", (size_t)ptr, (size_t)end);
			LOG_ERROR("");
		}
	}

	//5. test		DON'T MODIFY g_buf2 BELOW THIS POINT
	double sqe=0;
	ptrdiff_t usize=0;
	double ssim[4]={0}, ssim_weight=0;
	int maxencthreads=0, maxdecthreads=0;
	for(int k=0;k<(int)uinfo->count;++k)
	{
		UInfo *info=(UInfo*)array_at(&uinfo, k);
		usize+=info->usize;
	}
	TestInfo *currtest=(TestInfo*)ARRAY_APPEND(testinfo, 0, 1, 1, 0);
	STR_COPY(currtest->codecname, codecname, strlen(codecname));
	ARRAY_ALLOC(CellInfo, currtest->cells, 0, uinfo->count, 0, 0);
	print_summary(besttestidxs, testinfo, usize, -1, flags/CMDFLAG_PRINT_RIVALS);
	printf("\n");
	print_currtimestamp("%Y-%m-%d_%H%M%S");
	printf("  ");
	colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xC060FF, "%s", datasetname);
	printf(" ");
	colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xFFC060, "%s", codecname);
	printf("  %d images\n", (int)uinfo->count);
//	printf("  %s %d  %s\n", datasetname, (int)uinfo->count, codecname);
	{
		time_t t=time(0);
		struct tm *date=localtime(&t);
		currtest->datetime.year		=date->tm_year+1900;
		currtest->datetime.month	=date->tm_mon+1;
		currtest->datetime.day		=date->tm_mday;
		currtest->datetime.hour		=date->tm_hour;
		currtest->datetime.minute	=date->tm_min;
		currtest->datetime.second	=date->tm_sec;
	}
	for(int k=0;k<(int)uinfo->count;++k)
	{
		UInfo *info=(UInfo*)array_at(&uinfo, k);
		CellInfo *currcell=(CellInfo*)array_at(&currtest->cells, k);
//#ifdef _DEBUG
//		printf("\"%s\"\n", (char*)info->filename->data);//
//#endif

		int titlestart=0, titleend=0;
		get_filetitle((char*)info->filename->data, (int)info->filename->count, &titlestart, &titleend);

		//Print 1:  idx filetitle usize			//print filetitle first in case of CRASH
		printf("%5d %-*.*s %10lld",
			k+1,
			titlecolwidth, titleend-titlestart, (char*)info->filename->data+titlestart,
			info->usize
		);

		int success=CopyFileA((char*)info->filename->data, t1fn, 1);
		if(!success)
		{
			printf("Source:       %s\n", (char*)info->filename->data);
			printf("Destination:  %s\n", t1fn);
			SYSTEMERROR("CopyFileA");
			return 0;
		}

		int encthreads=0, decthreads=0;
		exec_process2(encline, (char*)currdir->data, 0, &currcell->etime, &currcell->emem, &encthreads);
		if(maxencthreads<encthreads)
			maxencthreads=encthreads;
		currcell->csize=get_filesize(t2fn);

		//Print 2:  -> csize etime
		printf(" -> %10lld B  %12.6lf", currcell->csize, currcell->etime);
		if(currcell->csize<=FSIZE_EMPTYFILE)
		{
			const char *msg="malicious";
			if(currcell->csize==FSIZE_EMPTYFILE)
				msg="empty";
			else if(currcell->csize==FSIZE_FOLDER)
				msg="a folder";
			else if(currcell->csize==FSIZE_INACCESSIBLE)
				msg="not found";
			printf("\n");
			printf("Temp. file #1 is %s:\n", msg);
			printf("  \"%s\"\n", t1fn);
			printf("The encode command was:\n");
			printf("  %s\n", encline);
			LOG_ERROR("");
		}
		
		ascii_deletefile(t1fn);
		exec_process2(decline, (char*)currdir->data, 0, &currcell->dtime, &currcell->dmem, &decthreads);
		if(maxdecthreads<decthreads)
			maxdecthreads=decthreads;
		
		//Print 3:  dtime espeed dspeed emem dmem  rivals
		printf(
			" %12.6lf sec  %12.6lf %12.6lf MB/s %8.2lf %8.2lf MB "
#ifdef TRACK_THREADS
			"%2d %2d "
#endif
			, currcell->dtime
			, info->usize/(currcell->etime*1024*1024)
			, info->usize/(currcell->dtime*1024*1024)
			, (double)currcell->emem/(1024*1024)
			, (double)currcell->dmem/(1024*1024)
#ifdef TRACK_THREADS
			, encthreads
			, decthreads
#endif
		);
		if(flags/CMDFLAG_PRINT_RIVALS)
			print_rivals_v2(besttestidxs, testinfo, k, currcell, info->usize);

		switch(flags%CMDFLAG_PRINT_RIVALS)
		{
		case CMDFLAG_VERIFY_BITEXACT:
			verify_files((char*)info->filename->data, t1fn);
			break;
		case CMDFLAG_PSNR_8BIT:
			measure_psnr_8bit((char*)info->filename->data, t1fn, currcell->csize, &sqe);
			break;
		case CMDFLAG_PSNR_16BIT:
			measure_psnr_16bit((char*)info->filename->data, t1fn, currcell->csize, &sqe);
			break;
		case CMDFLAG_SSIM_PPM:
		//	measure_ssim_ppm_avx2((char*)info->filename->data, t1fn, currcell->csize, ssim, &ssim_weight);
			measure_ssim_ppm((char*)info->filename->data, t1fn, currcell->csize, ssim, &ssim_weight);
			(void)measure_ssim_ppm_avx2;
			break;
		}
		ascii_deletefile(t1fn);
		ascii_deletefile(t2fn);
		printf("\n");

		currtest->total.csize+=currcell->csize;
		currtest->total.etime+=currcell->etime;
		if(currtest->total.emem<currcell->emem)
			currtest->total.emem=currcell->emem;
		currtest->total.dtime+=currcell->dtime;
		if(currtest->total.dmem<currcell->dmem)
			currtest->total.dmem=currcell->dmem;
	}
	printf("\n");
	//print summary
	{
		printf(
			"%5d %*s %10lld -> %10lld B  %12.6lf %12.6lf sec  %12.6lf %12.6lf MB/s %8.2lf %8.2lf MB "
#ifdef TRACK_THREADS
			"%2d %2d "
#endif
			, (int)uinfo->count
			, titlecolwidth, ""
			, usize
			, currtest->total.csize
			, currtest->total.etime
			, currtest->total.dtime
			, usize/(currtest->total.etime*1024*1024)
			, usize/(currtest->total.dtime*1024*1024)
			, (double)currtest->total.emem/(1024*1024)
			, (double)currtest->total.dmem/(1024*1024)
#ifdef TRACK_THREADS
			, maxencthreads
			, maxdecthreads
#endif
		);
		switch(flags%CMDFLAG_PRINT_RIVALS)
		{
		case CMDFLAG_PSNR_8BIT:
		case CMDFLAG_PSNR_16BIT:
			{
				int vmax=(flags%CMDFLAG_PRINT_RIVALS)==CMDFLAG_PSNR_8BIT?PSNR_MAX_8:PSNR_MAX_16;
				double rmse=sqrt(sqe/usize), psnr=20*log10(vmax/rmse);
				printf("  BPD%7.4lf rmse %12.6lf psnr %12.6lf", 8.*currtest->total.csize/usize, rmse, psnr);
			}
			break;
		case CMDFLAG_SSIM_PPM:
			printf("  BPD%7.4lf SSIM%9.6lf%9.6lf%9.6lf%9.6lf"
				, 8.*currtest->total.csize/usize
				, ssim[3]/ssim_weight
				, ssim[0]/ssim_weight
				, ssim[1]/ssim_weight
				, ssim[2]/ssim_weight
			);
		//	printf("  BPD%7.4lf SSIM%9.6lf", 8.*currtest->total.csize/usize, ssim/ssimden);
			break;
		}
#ifndef TRACK_THREADS
		(void)maxencthreads;
		(void)maxdecthreads;
#endif
		if(flags/CMDFLAG_PRINT_RIVALS)
			print_rivals_v2(besttestidxs, testinfo, -1, &currtest->total, usize);
		printf("\n");

		print_currtimestamp("%Y-%m-%d_%H%M%S");
		printf("  ");
		colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xC060FF, "%s", datasetname);
		printf(" ");
		colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xFFC060, "%s", codecname);
		printf("  %lf sec\n", currtest->total.etime+currtest->total.dtime);
		
		printf("\n");
		int rank=0, *idx=0;
		for(;rank<(int)besttestidxs->count;++rank)
		{
			idx=(int*)array_at(&besttestidxs, rank);
			TestInfo *test=(TestInfo*)array_at(&testinfo, *idx);
			int won=0;
			if(currtest->total.csize==test->total.csize)
				won=currtest->total.etime+currtest->total.dtime*2<test->total.etime+test->total.dtime*2;
			else
				won=currtest->total.csize<test->total.csize;
			if(won)
				break;
		}
		idx=(int*)array_insert(&besttestidxs, rank, 0, 1, 1, 0);
		*idx=(int)testinfo->count-1;
		print_summary(besttestidxs, testinfo, usize, rank, flags/CMDFLAG_PRINT_RIVALS);
	}

	//6. save		g_buf2 can be modified now
	{
		snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
		FILE *fdst=fopen(g_buf, "w");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", g_buf);
			return 0;
		}
		fprintf(fdst, "%s\t%s\nfiles:\n", srcpath->data, ext->data);
		for(int k=0;k<(int)uinfo->count;++k)
		{
			UInfo *info=(UInfo*)array_at(&uinfo, k);
			fprintf(fdst, "%10lld\t%s\n", info->usize, info->filename->data);
		}
		fprintf(fdst, "tests:\n");
		for(int k=0;k<(int)testinfo->count;++k)
		{
			TestInfo *info=(TestInfo*)array_at(&testinfo, k);
			fprintf(fdst, "%-20s\t%04d%02d%02d_%02d%02d%02d",
				(char*)info->codecname->data,
				info->datetime.year,
				info->datetime.month,
				info->datetime.day,
				info->datetime.hour,
				info->datetime.minute,
				info->datetime.second
			);
			//struct tm *t=localtime(&info->timestamp);
			//fprintf(fdst, "%-20s\t%04d%02d%02d_%02d%02d%02d",
			//	(char*)info->codecname->data,
			//	t->tm_year+1900,
			//	t->tm_mon+1,
			//	t->tm_mday,
			//	t->tm_hour,
			//	t->tm_min,
			//	t->tm_sec
			//);
			write_cell(fdst, &info->total);
			for(int k2=0;k2<(int)info->cells->count;++k2)
			{
				CellInfo *cell=(CellInfo*)array_at(&info->cells, k2);
				write_cell(fdst, cell);
			}
			fprintf(fdst, "\n");
		}
		fclose(fdst);
	}

	//7. free
	array_free(&enccmd.format);
	array_free(&deccmd.format);
	array_free(&srcpath);
	array_free(&ext);
	array_free(&uinfo);
	array_free(&testinfo);
	array_free(&tmpfn1);
	array_free(&tmpfn2);
	array_free(&besttestidxs);
	array_free(&currdir);
	return 0;
}