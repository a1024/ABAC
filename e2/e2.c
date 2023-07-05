#include "e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

int hist[256*4];

const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
};
void batch_test(const char *path)
{
	double t_start=time_ms();
	ArrayHandle filenames=get_filenames(path, g_extensions, COUNTOF(g_extensions), 1);
	if(!filenames)
	{
		printf("No images in \"%s\"\n", path);
		return;
	}
	long long
		count_PNG=0, count_JPEG=0,
		sum_cPNGsize=0, sum_cJPEGsize=0,
		sum_uPNGsize=0, sum_uJPEGsize=0,
		sum_testsize=0;
		//sum_test3size[2]={0};
#if 0	
	long long *hist=e10_start();//
	if(!hist)
	{
		exit(0);
		return;
	}
#endif

	for(ptrdiff_t k=0;k<(ptrdiff_t)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);

		if(!fn)
		{
			LOG_ERROR("filename read error");
			continue;
		}

		ptrdiff_t formatsize=get_filesize(fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images
			continue;

		int iw=0, ih=0, nch0=3, stride=4;
		long long cycles=__rdtsc();
		unsigned char *buf=image_load(fn[0]->data, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}

		ptrdiff_t res=(ptrdiff_t)iw*ih, len=res*stride, usize=res*nch0;
		double ratio=(double)usize/formatsize;
#ifdef BATCHTEST_NO_B2
		printf("%3lld/%3lld  %.2lf%%\r", k+1, filenames->count, (k+1)*100./filenames->count);
#else
		printf("%3lld/%3lld  \"%s\"\tCR %lf (%lf BPP) Dec %lf CPB", k+1, filenames->count, fn[0]->data, ratio, 8/ratio, (double)cycles/usize);
#endif
		if(!acme_stricmp(fn[0]->data+fn[0]->count-3, "PNG"))
		{
			sum_cPNGsize+=formatsize;
			sum_uPNGsize+=usize;
			++count_PNG;
		}
		else//assumed
		{
			sum_cJPEGsize+=formatsize;
			sum_uJPEGsize+=usize;
			++count_JPEG;
		}
#ifndef BATCHTEST_NO_B2
		unsigned char *b2=(unsigned char*)malloc(len);
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memset(b2, 0, len);
#endif
		
		//T34 ABAC + adaptive Bayesian inference
#if 1
		{
			ArrayHandle cdata=0;
			printf("\n");
#if 1
			memcpy(b2, buf, len);
			addbuf(b2, iw, ih, 3, 4, 128);
			colortransform_ycocb_fwd((char*)b2, iw, ih);
			pred_opt_opt_v6((char*)b2, iw, ih, 1);
			memset(b2, 0, len);
			pred_opt_printparam();
#endif

			printf("\nT35 (ABAC + context tree)\n");
			t35_encode(buf, iw, ih, &cdata, 1);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");

			t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T35", 0);

			//printf("\nT34 (ABAC + adaptive Bayesian inference)\n");
			//t34_encode(buf, iw, ih, &cdata, 1);

			printf("\n");
		}
#endif

		//T29
#if 0
		{
			ArrayHandle cdata=0;
			double elapsed;

			elapsed=time_ms();
			cycles=__rdtsc();
			t29_encode(buf, iw, ih, &cdata, 0);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
			
			array_free(&cdata);
		}
#endif

		//test26: T16 with range coder
#if 0
		{
			int use_ans=0;
			ArrayHandle cdata=0;
			double elapsed;
			T26Params params[]=
			{
				{ 8, 26, 26, 26, 0xD3, 52},
				{23, 37, 37, 37, 0xD3, 52},
				{ 8, 26, 26, 26, 0xD3, 52},
			};
			printf("\nT26 (%s)\n", use_ans?"ANS":"AC");
			//printf("\nT26\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t26_encode(buf, iw, ih, params, use_ans, &cdata, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t26_decode(cdata->data, cdata->count, iw, ih, params, use_ans, b2, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Dec %11lf CPB  ", (double)cycles/usize);
			timedelta2str(0, 0, elapsed);
			printf("\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T26", 0);
			memset(b2, 0, len);
			printf("\n");
		}
#endif

		//T25: T16 optimizer
#if 0
		{
			int use_ans=0;
			ArrayHandle cdata=0;
			//int blockw[]={96, 96, 96}, blockh[]={96, 96, 96};
			int blockw[]={128, 128, 128}, blockh[]={128, 128, 128};
			double elapsed;
			printf("\nT25 (%s)\n", use_ans?"ANS":"AC");
			elapsed=time_ms();
			cycles=__rdtsc();
			t25_encode(buf, iw, ih, blockw, blockh, use_ans, &cdata, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);
			
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t25_decode(cdata->data, cdata->count, iw, ih, blockw, blockh, use_ans, b2, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Dec %11lf CPB  ", (double)cycles/usize);
			timedelta2str(0, 0, elapsed);
			printf("\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T25", 0);
			memset(b2, 0, len);
			printf("\n");
		}
#endif
		
		//test16 - THE BEST
#if 0
		{
			printf("\nT16\n");
#if 1
			memcpy(b2, buf, len);
			addbuf(b2, iw, ih, 3, 4, 128);
			colortransform_ycocb_fwd((char*)b2, iw, ih);
			pred_opt_opt_v3((char*)b2, iw, ih, 1);
			memset(b2, 0, len);
			pred_opt_printparam();
#endif

			ArrayHandle cdata=0;
			int alpha=0xD3E7,
				blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
				blockh[]={ 1,  1,  1},
				margin[]={26, 37, 26};

			cycles=__rdtsc();
			test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, 0);
			cycles=__rdtsc()-cycles;
			printf("Enc %lf CPB  CR %lf  csize %lld", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");

			cycles=__rdtsc();
			test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
			cycles=__rdtsc()-cycles;
			printf("Dec %lf CPB\n", (double)cycles/usize);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, 3, 4, "T16", 0);
			memset(b2, 0, len);

			printf("\n");
		}
#endif

		//test16 estimate
#if 0
		apply_transforms_fwd(buf, iw, ih);
		double csize=test16_estimate_csize(buf, iw, ih, 32, 0.6f, 0);
		sum_testsize+=(long long)ceil(csize);
		printf("\tCR2 %f", usize/csize);
		if(csize<formatsize)
			printf(" !!!");
		printf("\n");
#endif

		//printf("\n");
		free(buf);
#ifndef BATCHTEST_NO_B2
		free(b2);
#endif
	}
	printf("Batch elapsed ");
	timedelta2str(0, 0, time_ms()-t_start);
	printf("\n");
#if 0
	e10_print(hist);
	free(hist);
#else
	ptrdiff_t totalusize=sum_uPNGsize+sum_uJPEGsize;
	if(totalusize)
	{
		printf("\nOn average:\n");
		printf("BMP     csize %9lld\n", totalusize);
		if(sum_cPNGsize)
			printf("PNG     csize %9lld  CR %lf  (%lld images)\n", sum_cPNGsize, (double)sum_uPNGsize/sum_cPNGsize, count_PNG);
		if(sum_cJPEGsize)
			printf("JPEG    csize %9lld  CR %lf  (%lld images)\n", sum_cJPEGsize, (double)sum_uJPEGsize/sum_cJPEGsize, count_JPEG);
		printf("test    csize %9lld  CR %lf\n", sum_testsize, (double)totalusize/sum_testsize);
		//printf("test3s  CR %lf\n", (double)totalusize/sum_test3size[0]);
		//printf("test3sd CR %lf\n", (double)totalusize/sum_test3size[1]);
	}
	else
		printf("\nNo valid images found\n");
#endif

	array_free(&filenames);

	printf("\nDone.\n");
	pause();
}
int main(int argc, char **argv)
{
	//const int LOL_1=(-1)/2;//0

	//int width=10,
	//	n=pyramid_getsize(width);
	//for(int k=0;k<n;++k)
	//	printf("%c", pyramid_getchar(width, k));
	//DCTtest();
	//DCTtest2();
	//test4();
	//test5();
	//test6();
	//test7();
	//test8();
	//test9();
	//test_swar();

	printf("EntropyBattle\n");
#if 1
	long long cycles;
	int iw=0, ih=0, nch0=3,
		nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	const char *fn=0;
#ifdef _DEBUG
	fn="D:/ML/dataset-kodak";
	//fn="D:/ML/dataset-kodak/kodim13.png";
#endif
	if(fn||argc==2)
	{
		if(!fn)
			fn=argv[1];
		ptrdiff_t formatsize=get_filesize(fn);
		if(formatsize==-1)
		{
			LOG_ERROR("Cannot open \"%s\"", fn);
			return 0;
		}
		if(!formatsize)//path
		{
			batch_test(fn);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		cycles=__rdtsc();
		buf=image_load(fn, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			LOG_ERROR("Couldn't open \"%s\"", fn);
			return 0;
		}
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		printf("Format Dec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", (double)cycles/(resolution*nch0), iw, ih, nch0, formatsize, (double)resolution*nch0/formatsize);
	}
	else if(argc==3)
	{
		const char *fn1=argv[1], *fn2=argv[2];
		int w2, h2;
		buf=image_load(fn1, &iw, &ih);
		b2 =image_load(fn2, &w2, &h2);
		if(!buf)
		{
			printf("Couldn't open %s\n", fn1);
			return 1;
		}
		if(!b2)
		{
			printf("Couldn't open %s\n", fn2);
			return 1;
		}
		if(iw!=w2||ih!=h2)
		{
			printf("Expected two images of SAME RESOLUTION. %dx%d != %dx%d\n", iw, ih, w2, h2);
			return 1;
		}
		ptrdiff_t formatsize=get_filesize(fn2);
		int res=iw*ih;
		long long sum[3]={0};
		for(int k=0;k<res;++k)
		{
			int dr=buf[k<<2  ]-b2[k<<2  ],
				dg=buf[k<<2|1]-b2[k<<2|1],
				db=buf[k<<2|2]-b2[k<<2|2];
			sum[0]+=dr*dr;
			sum[1]+=dg*dg;
			sum[2]+=db*db;
		}
		double rmse[]=
		{
			sqrt((double)sum[0]/res),
			sqrt((double)sum[1]/res),
			sqrt((double)sum[2]/res),
			sqrt((double)(sum[0]+sum[1]+sum[2])/(res*3)),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		double CR=res*3./formatsize;
		printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)formatsize, CR, 8/CR);
		printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
		printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
		printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
		return 0;
	}
	else
	{
		printf("Usage: e2.exe  file_or_path\n");
		pause();
		return 0;
#if 0
		iw=1920, ih=1080, nch0=3,//1080*1920*3	640*480		50		4*4*1
			nch=4;
		resolution=(size_t)iw*ih, len=resolution*nch;
		buf=(unsigned char*)malloc(len);
		if(!buf)
			return 0;
		//srand((unsigned)__rdtsc());
	
#ifdef UNIFORM
		printf("Generating test data (uniform)...\n");
		fill_uniform(buf, len);
#else
		int unibits=256;
		printf("Generating test data (%d bit binomial)...\n", unibits);
		fill_halfbinomial(buf, len, unibits);
#endif
#endif
	}

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<len;k+=nch)
			buf[k]=0xFF;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	size_t usize=len*nch0>>2;

	printf("\n");
	
	ArrayHandle cdata=0;
	//const void *ptr, *end;
	
	int loud=0;

	//test16
#if 0
	{
		//debug_ptr=buf;//
		int besta=0, bestb=0, bestm=0, bestc=0;
		int it=0;
		for(int m=30;m<=30;++m)
		{
			for(int b=15;b<=15;++b)
			{
				for(int a=0xD3E7;a<=0xD3E7;++a, ++it)
				{
					int alpha=a,//(60<<16)/100;	0xA51F	0xD3DA	0xD3D6	0xD3EB
						bsize=b,//32 24 26 20
						margin=m;
					printf("T16 a 0x%04X b %3d m %3d ", alpha, bsize, margin);
					cycles=__rdtsc();
					test16_encode(buf, iw, ih, alpha, bsize, margin, &cdata, 1);
					cycles=__rdtsc()-cycles;
					printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

					if(!it||bestc>(int)cdata->count)
						besta=alpha, bestb=bsize, bestm=m, bestc=(int)cdata->count;

					cycles=__rdtsc();
					test16_decode(cdata->data, cdata->count, iw, ih, alpha, bsize, margin, b2);
					cycles=__rdtsc()-cycles;
					printf("Dec %11lf CPB ", (double)cycles/usize);

					array_free(&cdata);
					compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
					memset(b2, 0, len);
				}
				//printf("\n");
			}
		}
		if(it>1)
			printf("T16 best ABMS 0x%04X %2d %2d %d CR %lf\n", besta, bestb, bestm, bestc, (double)usize/bestc);
		else
			printf("\n");
	}
#endif


	//test16 codec with jxl predictor optimizer
#if 0
	{
		int alpha=0xD3E7,
			blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
			blockh[]={ 1,  1,  1},
			margin[]={26, 37, 26};

#if 0
		int res=iw*ih;
		double step=0.001, CR0=0, CR, csize[3]={0};
		estimate_csize_from_transforms(buf, b2, iw, ih, csize);
		CR=res*3/(csize[0]+csize[1]+csize[2]);
		printf("%4d TRGB %lf [%lf %lf %lf]\n", 0, CR, res/csize[0], res/csize[1], res/csize[2]);

		for(int k=0;k<256;++k)
		{
			int idx=k%33;
			if(!(k+1)%33)
				step*=0.9;
			do
			{
				CR0=CR;
				jxlpred_params[idx]+=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);

			do
			{
				CR0=CR;
				jxlpred_params[idx]-=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);
		}
#endif

		//for(int k=0;k<3;++k)
		//	printf("%g\t%g\t%g\t%g\n", jxlpred_params[k<<2], jxlpred_params[k<<2|1], jxlpred_params[k<<2|2], jxlpred_params[k<<2|3]);
		//for(int k=12;k<33;k+=7)
		//	printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", jxlpred_params[k], jxlpred_params[k+1], jxlpred_params[k+2], jxlpred_params[k+3], jxlpred_params[k+4], jxlpred_params[k+5], jxlpred_params[k+6]);
		//for(int k=0;k<33;++k)
		//	printf("%3d  %lf\n", k, jxlpred_params[k]);
		//printf("\n");
		
		printf("T16\n");
		int bestcsizes[3]={0}, bestw[3]={0}, besth[3]={0}, bestm[3]={0};
		int it=0;
		//for(int bw=19;bw<25;++bw)//
		{
			//for(int m=31;m<48;++m)
			{
				int csizes[3];
				//blockw[1]=bw, blockh[1]=1, margin[1]=m;
				//blockw[1]=1+k, blockh[1]=1;
				//blockw[1]=4+k%10, blockh[1]=1+k/10;
			
				//blockw[0]=blockw[2]=blockw[1];
				//blockh[0]=blockh[2]=blockh[1];

				cycles=__rdtsc();
				test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, csizes);
				cycles=__rdtsc()-cycles;
				printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

				for(int kc=0;kc<3;++kc)
				{
					if(!it||bestcsizes[kc]>csizes[kc])
						bestcsizes[kc]=csizes[kc], bestw[kc]=blockw[kc], besth[kc]=blockh[kc], bestm[kc]=margin[kc];
				}

				cycles=__rdtsc();
				test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
				cycles=__rdtsc()-cycles;
				printf("Dec %11lf CPB ", (double)cycles/usize);

				array_free(&cdata);
				compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
				memset(b2, 0, len);
				printf("\n");
				++it;
			}
		}
		int res=iw*ih;
		printf("R %7d %lf %2dx%d  M %d\n", bestcsizes[0], (double)res/bestcsizes[0], bestw[0], besth[0], bestm[0]);
		printf("G %7d %lf %2dx%d  M %d\n", bestcsizes[1], (double)res/bestcsizes[1], bestw[1], besth[1], bestm[1]);
		printf("B %7d %lf %2dx%d  M %d\n", bestcsizes[2], (double)res/bestcsizes[2], bestw[2], besth[2], bestm[2]);
		printf("\n");
	}
#endif

	//test25: T16 optimizer
#if 0
	{
		int use_ans=0;
		printf("T25 (%s)\n", use_ans?"ANS":"AC");
		double elapsed;
		int blockw[]={128, 128, 128}, blockh[]={128, 128, 128};
		//int lbsizes[]=
		//{
		//	32, 32,
		//	32, 32,
		//	32, 32,
		//};
		//int sbsizes[]=//unused
		//{
		//	16, 16,
		//	16, 16,
		//	16, 16,
		//};
		elapsed=time_ms();
		cycles=__rdtsc();
		t25_encode(buf, iw, ih, blockw, blockh, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t25_decode(cdata->data, cdata->count, iw, ih, blockw, blockh, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T25", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//test26: T16 with range coder
#if 0
	{
		int use_ans=0;
		double elapsed;
		T26Params params[]=
		{
			{ 8, 26, 26, 26, 0xD3, 52},
			{23, 37, 37, 37, 0xD3, 52},
			{ 8, 26, 26, 26, 0xD3, 52},
		};
		printf("T26 (%s)\n", use_ans?"ANS":"AC");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t26_encode(buf, iw, ih, params, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t26_decode(cdata->data, cdata->count, iw, ih, params, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T26", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//T27: ABAC
#if 1
	//printf("T29 (ABAC + Bayesian inference)\n");
	//t29_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T30 (ABAC + Bayesian inference + predictor)\n");	//X
	//t30_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T31 (ABAC + adaptive Bayesian inference)\n");
	//t31_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T33 (ABAC + adaptive Bayesian inference with circular buffer)\n");	//X
	//t33_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T34 (ABAC + adaptive Bayesian inference)\n");
	//t34_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);
	
	printf("T35 Entropy coding with context tree\n");
	//printf("T35 Combines spatial transform with entropy coding\n");
	t35_encode(buf, iw, ih, &cdata, 1);
	t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T35", 0);
	memset(b2, 0, len);
	printf("\n");
#endif

	//predict image
	apply_transforms_fwd(buf, iw, ih);
	//lodepng_encode_file("kodim21-YCoCgT-unplane.PNG", buf, iw, ih, LCT_RGBA, 8);//
	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//

	//colortransform_ycocb_fwd(buf, iw, ih);
	//save_channel(buf+1, iw, ih, 4, 0, "kodim13-YCoCb-jxl-luma.PNG");
#if 0
	printf("Predict image...\n");
	{
		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		char *temp=(char*)malloc(MAXVAR(iw, ih));
		
		addbuf(buf, iw, ih, nch0, nch, 128);//unsigned char -> signed char
		
		//colortransform_ycocg_fwd((char*)buf, iw, ih);
		//colortransform_xgz_fwd((char*)buf, iw, ih);
		//colortransform_xyz_fwd((char*)buf, iw, ih);

		//char *b3=(char*)malloc(iw), *b4=(char*)malloc(iw);
		//if(!b3||!b4)
		//	return 0;
		//memcpy(b3, buf, iw);
		//dwt1d_squeeze_fwd(b3, iw, 1, b4);
		//dwt1d_squeeze_inv(b3, iw, 1, b4);
		//compare_bufs_uint8((unsigned char*)b3, buf, iw, 1, 1, 1, "squeeze row", 0);
		//free(b3);
		//free(b4);

#if 1
		memcpy(b2, buf, len);

		colortransform_ycocb_fwd((char*)b2, iw, ih);
		float jxlparams[33]=
		{
			 0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
			 0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
			 0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
		};
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 1);
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 0);

		//colortransform_xyz_fwd(b2, iw, ih);
		//colortransform_xyz_inv(b2, iw, ih);
		//for(int kc=0;kc<3;++kc)
		//{
		//	dwt2d_squeeze_fwd((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//	dwt2d_squeeze_inv((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//}
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "transform", 0);
		printf("\n");
#endif
		
		//for(int kc=0;kc<3;++kc)
		//	//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	dwt2d_squeeze_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, 2, 4, (char*)temp);
		//	//dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

		addbuf(buf, iw, ih, nch0, nch, 128);

		//save_DWT_int8("kodim21-squeeze-stage", buf, (DWTSize*)sizes->data, 2, 4);//
		//lodepng_encode_file("kodim21-cubic.PNG", buf, iw, ih, LCT_RGBA, 8);//

		array_free(&sizes);
		free(temp);
	}
	//squeeze_8bit_lossy(buf, iw, ih, nch0, nch);
//	image_pred(buf, iw, ih, nch0, nch);

	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//
#endif

	for(int kc=0;kc<nch0;++kc)
		calc_histogram(buf+kc, len, nch, hist+((size_t)kc<<8));
	//print_histogram(hist, 1);

	double entropy[6]={0};
	for(int kc=0;kc<nch0;++kc)
	{
		int freq;
		double p;
		for(int k=0;k<256;++k)
		{
			freq=hist[kc<<8|k];
			if(freq)
			{
				p=(double)freq/(len>>2);
				p*=0x10000-255;
				++p;
				p/=0x10000;
				entropy[kc]+=-p*log2(p);
			}

			//if(!kc)//
			//	printf("%3d %6d %lf\n", k, freq, p);//
		}
		printf("ch %d E = %lf / 8, optimal ratio = %lf\n", kc, entropy[kc], 8/entropy[kc]);
		entropy[4]+=entropy[kc];
	}
	entropy[4]/=nch0;
	printf("Av. E = %lf / 8, optimal ratio = %lf <- true limit\n", entropy[4], 8/entropy[4]);


	free(buf);
	free(b2);
#endif
	printf("Done.\n");
	pause();
	return 0;
}