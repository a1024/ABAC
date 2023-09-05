#include"lossy.h"
#include<stdio.h>
#include<math.h>
#include"lodepng.h"
#include"turbojpeg.h"
#pragma comment(lib, "turbojpeg.lib")
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
static const char file[]=__FILE__;
#define CHECKTJ(E, CTX)	(!(E)||LOG_ERROR("%s", tjGetErrorStr2(CTX)))


const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
};
void batch_test(const char *path)
{
	int known_dataset=0;
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Start %s\n", g_buf);
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
	{
		ArrayHandle path2=filter_path(path);
		if(!acme_stricmp(path2->data, "D:/ML/dataset-CLIC30/"))
			known_dataset=1;
		else if(!acme_stricmp(path2->data, "D:/ML/dataset-Kodak/"))
			known_dataset=2;
		array_free(&path2);
	}

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
		unsigned char *buf=stbi_load(fn[0]->data, &iw, &ih, 0, 4);
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
		
		//T34+: ABAC + adaptive Bayesian inference
#if 0
		{
			ArrayHandle cdata=0;
			printf("\n");
#if 0
			if(known_dataset==1&&k<30)//pre-trained
				g_param_ptr=jxl_CLIC30_params+33*k;
			else if(known_dataset==2&&k<24)
				g_param_ptr=jxl_Kodak_params+33*k;
			else
			{
				memcpy(b2, buf, len);
				addbuf(b2, iw, ih, 3, 4, 128);
				colortransform_ycocb_fwd((char*)b2, iw, ih);
				pred_opt_opt_v6((char*)b2, iw, ih, 1);
				memset(b2, 0, len);
				pred_opt_printparam();
			}
#endif

			//printf("\nT35 (ABAC + context tree)\n");
			//printf("\nT39 Multiple estimators for all maps  WH %d*%d\n", iw, ih);

			//t35_encode(buf, iw, ih, &cdata, 1);
			//t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			//t39_encode(buf, iw, ih, &cdata, 1);//prev record
			//t39_decode(cdata->data, cdata->count, iw, ih, b2, 1);//prev record

			//t40_encode(buf, iw, ih, &cdata, 1);
			//t40_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			//t42_encode(buf, iw, ih, &cdata, 1);
			//t42_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			//t43_encode(buf, iw, ih, &cdata, 1);
			//t43_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "Test", 0);

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
	if(totalusize&&sum_testsize)
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
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);

	printf("\nDone.\n");
	pause();
}
int main(int argc, char **argv)
{
	//void test_mq();
	//test_mq();

	//for(int k=0;k<16;++k)
	//{
	//	int x=k&3, y=k>>2;
	//	int d=(y-x)&3, a=(y+(d>>1))&3;
	//	//d&=3, a&=3;
	//	printf("%d%d%d%d%c", d>>1, d&1, a>>1, a&1, x==3?'\n':'\t');
	//	//printf("%d %d%c", d, a, x==3?'\n':'\t');
	//}
	//pause();
	//exit(0);

	//void DCT_test();
	//DCT_test();

	printf("Lossy Compression Benchmark\n");
	
	long long cycles;
	int iw=0, ih=0, nch0=3,
		nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	const char *fn=0;
	double rmse[4], psnr[4];
#ifdef _DEBUG
	//fn="C:/Projects/datasets/CLIC11-crop4-2.PNG";
	//fn="C:/Projects/datasets/CLIC11-small4.PNG";
	//fn="C:/Projects/datasets/dataset-CLIC30/11.png";
	//fn="C:/Projects/datasets/dataset-kodak";
	fn="C:/Projects/datasets/dataset-kodak/kodim13.png";

	//fn="D:/ML/dataset-CLIC30";
	//fn="D:/ML/dataset-kodak";
	//fn="D:/ML/dataset-CLIC30/16.png";//hardest noiseless CLIC30 image
	//fn="D:/ML/dataset-CLIC30/17.png";
	//fn="D:/ML/dataset-kodak/kodim13.png";
	//fn="D:/ML/dataset-kodak-small/13.PNG";
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
			LOG_ERROR("");//
			batch_test(fn);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		cycles=__rdtsc();
		buf=stbi_load(fn, &iw, &ih, 0, 4);
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
		buf=stbi_load(fn1, &iw, &ih, 0, 4);
		b2 =stbi_load(fn2, &w2, &h2, 0, 4);
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

		measure_distortion(buf, b2, iw, ih, rmse, psnr);
		
		int res=iw*ih;
		double CR=res*3./formatsize;
		printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)formatsize, CR, 8/CR);
		printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
		printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
		printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
		pause();
		return 0;
	}
	else
	{
		printf("Usage: lossy.exe  file_or_path\n");
		pause();
		return 0;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	size_t usize=len*nch0>>2;
	printf("\n");


	ArrayHandle data=0;
	int res=iw*ih;
#if 0
	t46_encode(buf, iw, ih, &data, 1);
#endif

	//JPEG		R-D graph
#if 0
	{
		int error=0;
		unsigned char *jpegbuf=0;
		unsigned long jpegsize=0;
		tjhandle encoder=tjInitCompress();
		tjhandle decoder=tjInitDecompress();
		for(int q=0;q<=100;++q)
		{
			double t_enc=time_ms();
			error=tjCompress2(encoder, buf, iw, 0, ih, TJPF_RGBA, &jpegbuf, &jpegsize, TJSAMP_420, q, TJFLAG_FASTDCT);	CHECKTJ(error, encoder);
			t_enc=time_ms()-t_enc;
			
			int subsamp=0;
			double t_dec=time_ms();
			error=tjDecompressHeader2(decoder, jpegbuf, jpegsize, &iw, &ih, &subsamp);	CHECKTJ(error, decoder);
			error=tjDecompress2(decoder, jpegbuf, jpegsize, b2, iw, 0, ih, TJPF_RGBA, TJFLAG_FASTDCT);	CHECKTJ(error, decoder);
			t_dec=time_ms()-t_dec;

			double CR=iw*ih*3./jpegsize;
			measure_distortion(buf, b2, iw, ih, rmse, psnr);
			printf("%3d%% Q  %15.6lf BPP  %15.6lf PSNR  Enc %15.6lf ms  Dec %15.6lf ms\n", q, 8/CR, psnr[3], t_enc, t_dec);
			
			tjFree(jpegbuf);
			jpegsize=0;
			jpegbuf=0;
		}
		tjDestroy(encoder);
		tjDestroy(decoder);
	}
#endif

	//T47 pseudo-J2K	R-D graph
#if 0
	extern float t47_q;
	//psnr[3]=0;//
	for(t47_q=1;t47_q<=1024;t47_q*=2)
	{
		t47_encode(buf, iw, ih, &data, 0);
		t47_decode(data->data, data->count, iw, ih, b2, 0);
		measure_distortion(buf, b2, iw, ih, rmse, psnr);

		double CR=res*3./data->count;

		printf("%15.6lf Q  %15.6lf BPP  %15.6lf PSNR\n", t47_q, 8/CR, psnr[3]);
		array_free(&data);
	}
#endif

	//T45 pseudo-JPEG	R-D graph
#if 0
	extern float t45_q;
	//psnr[3]=0;//
	for(t45_q=0.0625f;t45_q<=128;t45_q*=2)
	{
		t45_encode(buf, iw, ih, &data, 0);
		t45_decode(data->data, data->count, iw, ih, b2, 0);
		measure_distortion(buf, b2, iw, ih, rmse, psnr);

		double CR=res*3./data->count;

		printf("%15.6lf Q  %15.6lf BPP  %15.6lf PSNR\n", t45_q, 8/CR, psnr[3]);
		array_free(&data);
	}
#endif

	//testbench
#if 1
	t47_encode(buf, iw, ih, &data, 1);
	t47_decode(data->data, data->count, iw, ih, b2, 1);

	//t46_encode(buf, iw, ih, &data, 1);
	//t46_decode(data->data, data->count, iw, ih, b2, 1);

	//t45_encode(buf, iw, ih, &data, 1);
	//t45_decode(data->data, data->count, iw, ih, b2, 1);

	measure_distortion(buf, b2, iw, ih, rmse, psnr);

	double CR=res*3./data->count;
	printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)data->count, CR, 8/CR);
	printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
	printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
	printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
	
	array_free(&data);
	//system("cd");
	lodepng_encode_file("dump.PNG", b2, iw, ih, LCT_RGBA, 8);
#endif


	free(buf);
	free(b2);

	printf("Done.\n");
	pause();
	return 0;
}
