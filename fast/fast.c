#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#ifdef _MSC_VER
#include<intrin.h>
#include<Windows.h>
#include<process.h>
#define THREAD_CALL __stdcall
typedef unsigned THREAD_RET;
#else
#include<x86intrin.h>
#include<pthread.h>
#define THREAD_CALL
typedef void *THREAD_RET;
#endif
#include"huffman.h"
static const char file[]=__FILE__;


#define CODECID     19
#define CODECNAME "F19"
#define ENCODE     f19_encode
#define DECODE     f19_decode


static const char *g_extensions[]=
{
	"png",
	"jpg", "jpeg",
	"ppm", "pgm",
	"bmp",
	"tif", "tiff",
};
typedef struct ThreadArgsStruct
{
	ArrayHandle title;
	Image src, dst;
	size_t usize, csize1, csize2;
	double fdec, enc, dec;//time in secs
	int error, unused;//whether the image was recovered successfully
	ptrdiff_t idx;
} ThreadArgs;
static THREAD_RET THREAD_CALL sample_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	ArrayHandle cdata=0;

	args->usize=image_getBMPsize(&args->src);

	double t=time_sec();
	ENCODE(&args->src, &cdata, 0);
	t=time_sec()-t;
	args->enc=t;

	args->csize2=cdata->count;
	
	t=time_sec();
	DECODE(cdata->data, cdata->count, &args->dst, 0);
	t=time_sec()-t;
	args->dec=t;

	array_free(&cdata);
	args->error=compare_bufs_16(args->dst.data, args->src.data, args->src.iw, args->src.ih, args->src.nch, args->src.nch, CODECNAME, 0, 0);
	image_clear(&args->src);
	image_clear(&args->dst);
	return 0;
}
typedef struct ResultStruct
{
	//ptrdiff_t idx;
	size_t usize, csize1, csize2, error;
	double fdec, enc, dec;
} Result;
typedef struct ProcessCtxStruct
{
	int nstarted, nfinished;
	ArrayHandle threadargs;//<ThreadArgs>	*thread_count
	ArrayHandle results;//<Result>		*nsamples
} ProcessCtx;
static double start_time=0, check_time=0;
static double g_total_usize=0, g_total_csize=0;
static void print_result(Result *res, const char *title, int width, int print_timestamp)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	g_total_usize+=res->usize;
	g_total_csize+=res->csize2;
	printf("%-*s  %10zd  format %10zd %10.6lf%% D %12lf sec  test %10zd %10.6lf%% E %12lf D %12lf sec %s  all %10.6lf%%",
		width, title, res->usize,
		res->csize1, 100./CR1, res->fdec,
		res->csize2, 100./CR2, res->enc, res->dec, res->error?"ERROR":"SUCCESS",
		100.*g_total_csize/g_total_usize
	);
	if(print_timestamp)
	{
		printf(" ");
		double t=time_sec();
		if(print_timestamp==2)
			timedelta2str(0, 0, t-start_time);
		else
		{
			timedelta2str(0, 0, t-check_time);
			check_time=t;
		}
	}
	printf("\n");
}
static void process_file(ProcessCtx *ctx, ArrayHandle title, int maxlen, Image *image, size_t csize1, double fdec, ptrdiff_t idx, int nthreads)
{
	if(!ctx->nstarted)
	{
		ARRAY_ALLOC(ThreadArgs, ctx->threadargs, 0, 0, nthreads, 0);
		ARRAY_ALLOC(Result, ctx->results, 0, 0, 0, 0);
	}

	if(image)
	{
		ThreadArgs *threadargs=(ThreadArgs*)ARRAY_APPEND(ctx->threadargs, 0, 1, 1, 0);
		size_t nvals=(size_t)image->iw*image->ih*image->nch;
		threadargs->src=*image;
		memset(&threadargs->dst, 0, sizeof(threadargs->dst));
		int success=image_copy_nodata(&threadargs->dst, image);
		if(!success)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(threadargs->dst.data, 0, nvals*sizeof(short));
		threadargs->title=title;
		threadargs->usize=0;
		threadargs->csize1=csize1;
		threadargs->csize2=0;
		threadargs->fdec=fdec;
		threadargs->enc=0;
		threadargs->dec=0;
		threadargs->error=0;
		threadargs->idx=idx;
		++ctx->nstarted;
	}

	if(ctx->nstarted==1)//first sample initializes memory, can't parallelize yet
	{
		ThreadArgs *threadargs=(ThreadArgs*)array_at(&ctx->threadargs, ctx->threadargs->count-1);
		sample_thread(threadargs);
		Result result=
		{
			threadargs->usize,
			threadargs->csize1,
			threadargs->csize2,
			(size_t)threadargs->error,
			threadargs->fdec,
			threadargs->enc,
			threadargs->dec,
		};
		ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
		print_result(&result, (char*)threadargs->title->data, maxlen, 1);

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
	else if(ctx->nstarted-ctx->nfinished>=nthreads)
	{
		int n=ctx->nstarted-ctx->nfinished;
		ArrayHandle handles;
#ifdef _MSC_VER
		ARRAY_ALLOC(HANDLE, handles, 0, n, 0, 0);
		for(int k=0;k<n;++k)
		{
			HANDLE *h=(HANDLE*)array_at(&handles, k);
			ThreadArgs *threadargs=(ThreadArgs*)array_at(&ctx->threadargs, k);
			*h=(void*)_beginthreadex(0, 0, sample_thread, threadargs, 0, 0);
			if(!*h)
			{
				LOG_ERROR("Alloc error");
				return;
			}
		}
		WaitForMultipleObjects(n, (HANDLE*)handles->data, TRUE, INFINITE);
		for(int k=0;k<n;++k)
		{
			HANDLE *h=(HANDLE*)array_at(&handles, k);
			CloseHandle(*h);
		}
#else
		ARRAY_ALLOC(pthread_t, handles, 0, n, 0, 0);
		for(int k=0;k<n;++k)
		{
			pthread_t *h=(pthread_t*)array_at(&handles, k);
			ThreadArgs *threadargs=(ThreadArgs*)array_at(&ctx->threadargs, k);
			int error=pthread_create(h, 0, sample_thread, threadargs);
			if(error)
			{
				LOG_ERROR("Alloc error");
				return;
			}
		}
		for(int k=0;k<n;++k)
		{
			pthread_t *h=(pthread_t*)array_at(&handles, k);
			pthread_join(*h, 0);
		}
#endif
		
		for(int k=0;k<n;++k)
		{
			ThreadArgs *threadargs=(ThreadArgs*)array_at(&ctx->threadargs, k);
			Result result=
			{
				threadargs->usize,
				threadargs->csize1,
				threadargs->csize2,
				(size_t)threadargs->error,
				threadargs->fdec,
				threadargs->enc,
				threadargs->dec,
			};
			ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
			print_result(&result, (char*)threadargs->title->data, maxlen, k>=n-1);
		}

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
}
static void batch_test_mt(const char *path, int nthreads)
{
	g_total_usize=0;
	g_total_csize=0;
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("%s\n", g_buf);
	printf("Multithreaded Batch Test %s\n", CODECNAME);
	double t_start=time_sec();
	check_time=start_time=t_start;
	ArrayHandle filenames=get_filenames(path, g_extensions, _countof(g_extensions), 1);
	if(!filenames)
	{
		printf("No supported images in \"%s\"\n", path);
		return;
	}
	ArrayHandle titles;
	ARRAY_ALLOC(ArrayHandle, titles, 0, 0, filenames->count, (void(*)(void*))array_free);
	int width=6;//"Total:"
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		int start=0, end=0;
		get_filetitle((char*)fn[0]->data, (int)fn[0]->count, &start, &end);
		ArrayHandle title;
		STR_COPY(title, (char*)fn[0]->data+start, end-start);
		if(width<(int)title->count)
			width=(int)title->count;
		ARRAY_APPEND(titles, &title, 1, 1, 0);
	}
	ProcessCtx processctx={0};
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		ArrayHandle *title=(ArrayHandle*)array_at(&titles, k);

		ptrdiff_t formatsize=get_filesize((char*)fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images, this check is useless because get_filenames() has already filtered the list
			continue;

		Image image={0};
		double t=time_sec();
		int success=image_load((char*)fn[0]->data, &image);
		t=time_sec()-t;
		if(!success)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}
		process_file(&processctx, *title, width, &image, formatsize, t, k, nthreads);
		//printf("");
	}
	process_file(&processctx, 0, width, 0, 0, 0, 0, 0);//set nthreads=0 to flush queued images
	if(processctx.results)
	{
		Result total={0};
		for(int k=0;k<(int)processctx.results->count;++k)
		{
			Result *result=(Result*)array_at(&processctx.results, k);

			total.usize+=result->usize;
			total.csize1+=result->csize1;
			total.csize2+=result->csize2;
			total.error+=result->error;
			total.fdec+=result->fdec;
			total.enc+=result->enc;
			total.dec+=result->dec;
		}
		printf("\n");
		print_result(&total, "Total:", width, 2);
		array_free(&processctx.results);
	}
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);
}

int main(int argc, char **argv)
{
#if 0
	{
#ifdef _DEBUG
		const char *path=
		//	"D:/ML/dataset-CLIC30-ppm"
		//	"D:/ML/dataset-CLIC"
			"D:/ML/dataset-kodak-ppm"
			;
#else
		if(argc!=2)
			return 0;
		const char *path=argv[1];
#endif
		f12_statstest(path);
		//f10_mptest("D:/ML/dataset-sintel-ppm");
		pause();
		return 0;
	}
#endif
	printf("FastEntropy\n\n");
#ifndef _DEBUG
//#if 0
	if((unsigned)(argc-2)>1)
	{
		printf(
			"Usage:\n"
			"  %s  filename\n"
			"  %s  folder [nthreads]\n",
			argv[0], argv[0]
		);
		pause();
		return 1;
	}
	const char *fn=argv[1];
#else
	const char *fn=
		"D:/ML/dataset-kodak/kodim13.png"
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
	//	"D:/ML/big_building.PPM"
	//	"C:/dataset-LPCB-ppm/PIA13785.ppm"
	//	"C:/dataset-LPCB-ppm/STA13456.ppm"	//uncorrelated channels
	//	"C:/dataset-LPCB-ppm/PIA13799.ppm"
	//	"D:/ML/dataset-RAW/a0001-jmac_DSC1459.dng"

	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim13.ppm"
	//	"C:/Projects/datasets/dataset-kodak/kodim13.png"
	//	"C:/Projects/datasets/kodim13-small4.PNG"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-ic-rgb16bit/deer.png"
	//	"C:/Projects/datasets/jupiter.PNG"
	//	"C:/Projects/datasets/space-8k-CROP.PPM"
	//	"C:/Projects/datasets/dataset-CLIC30"
	//	"C:/Projects/datasets/dataset-LPCB-ppm"
		;
#endif
	ptrdiff_t formatsize=get_filesize(fn);
	if(formatsize<0)
	{
		printf("Not a file or directory:  \'%s\'\n", fn);
		pause();
		return 1;
	}
	if(!formatsize)
	{
		int nthreads=4;
		if(argc==3)
			nthreads=atoi(argv[2]);
		batch_test_mt(fn, nthreads);
	}
	else
	{
		printf("File: \"%s\"\n", fn);
		Image src={0};
		double t=time_sec();
		int success=image_load(fn, &src);
		t=time_sec()-t;
		if(!success)
		{
			printf("Cannot open \'%s\'\n", fn);
			return 1;
		}
		size_t usize=image_getBMPsize(&src);
		printf("CWHD %d*%d*%d*%d/8 = %zd bytes  format %td bytes %10.6lf%% CR %lf  D %lf sec\n",
			src.nch, src.iw, src.ih, src.depth,
			usize,
			formatsize, 100.*formatsize/usize, (double)usize/formatsize,
			t
		);
		printf("%s:\n", CODECNAME);
		Image dst={0};
		image_copy_nodata(&dst, &src);

		ArrayHandle cdata=0;
#if 0
		extern int f09_disable_ctx;
		ptrdiff_t basesize=0;
		for(f09_disable_ctx=-1;f09_disable_ctx<F09_NCTX;++f09_disable_ctx)
		{
			ENCODE(&src, &cdata, 0);
			DECODE(cdata->data, cdata->count, &dst, 0);
			if(f09_disable_ctx==-1)
				basesize=cdata->count;
			int error=compare_bufs_16(dst.data, src.data, src.iw, src.ih, src.nch, src.nch, CODECNAME, 0, 0);
			ptrdiff_t dispsize=f09_disable_ctx==-1?basesize:cdata->count-basesize;
			printf("Disable CTX %2d  %12lld  %s\n", f09_disable_ctx, dispsize, f09_disable_ctx>=0?f09_prednames[f09_disable_ctx]:"");
			if(error)
				LOG_ERROR("ERROR");
			array_free(&cdata);
		}
#endif

		//single test
#if 1
		ENCODE(&src, &cdata, 1);
		
		DECODE(cdata->data, cdata->count, &dst, 1);

		//pred_wp_deferred(&dst, 1);
		//pred_wp_deferred(&dst, 0);
		//char depths[]={dst.depth, dst.depth, dst.depth};
		//pred_clampgrad(&dst, 1, depths);
		//pred_clampgrad(&dst, 0, depths);
		//rct_JPEG2000_32(&dst, 1);
		//rct_JPEG2000_32(&dst, 0);
		//pred_avx2(&dst, 1, depths);
		//pred_avx2(&dst, 0, depths);
		
		compare_bufs_16(dst.data, src.data, src.iw, src.ih, src.nch, src.nch, CODECNAME, 0, 1);
#endif
#if 0
		ptrdiff_t nvals=(ptrdiff_t)src.nch*src.iw*src.ih;
		success=1;
		for(ptrdiff_t k=0;k<nvals;++k)
		{
			if(dst.data[k]!=src.data[k])
			{
				ptrdiff_t k2=k;
				int kc=(int)(k2%src.nch);
				k2/=src.nch;
				int kx=(int)(k2%src.iw);
				k2/=src.iw;
				int ky=(int)k2;
				printf(
					"ERROR\n"
					"  decoded != original\n"
					"  %d != %d\n"
					"  0x%04X != 0x%04X\n"
					"  CXY %d %d %d\n",
					dst.data[k], src.data[k],
					(unsigned short)dst.data[k], (unsigned short)src.data[k],
					kc, kx, ky
				);
				success=0;
				LOG_ERROR("");
				break;
			}
		}
		if(success)
			printf("%s: SUCCESS\n", CODECNAME);
#endif

		//huff_test(&src);

#if 0
		char depths[]=
		{
			src.depth,
			src.depth+1,
			src.depth+1,
			0,
		};
		double csizes[5]={0}, csizes_vlc[5]={0}, csizes_shannon[5]={0}, csizes_bin[5]={0};
		double times[6]={0};

		times[0]=time_sec();
		rct_JPEG2000_32(&src, 1);
		times[0]=time_sec()-times[0];
		
		times[1]=time_sec();
		pred_simd(&src, 1, depths);
		//pred_clampgrad_fast(&src, 1, depths);
		//pred_clampgrad(&src, 1, depths);
		times[1]=time_sec()-times[1];
		
		times[2]=time_sec();
		calc_csize(&src, depths, csizes_shannon, 0);
		times[2]=time_sec()-times[2];
		
		times[3]=time_sec();
		calc_csize_vlc(&src, depths, csizes, csizes_vlc);
		times[3]=time_sec()-times[3];
		
		times[4]=time_sec();
		//calc_csize_bin(&src, depths, csizes_bin);
		times[4]=time_sec()-times[4];
		
		size_t csize_abac=0;
		times[5]=time_sec();
		csize_abac=calc_csize_ABAC(&src, depths);
		times[5]=time_sec()-times[5];

		printf("\nC       Shannon       Ada-Zipf        Ada-VLC        Ada-Bin\tusize%12lld\n", image_getBMPsize(&src));
		for(int k=0;k<src.nch+1;++k)
			printf("%c%14.2lf %14.2lf %14.2lf %14.2lf\n", "TYUVA"[k], csizes_shannon[k], csizes[k], csizes_vlc[k], csizes_bin[k]);
		printf("csize_ABAC %14lld\n", csize_abac);
		const char *oplabels[]=
		{
			"RCT_JPEG2000",
			"clampgrad",
			"Shannon",
			"Ada-Zipf/VLC",
			"Ada-Bin",
			"ABAC",
		};
		for(int k=0;k<_countof(times);++k)
			printf("%-16s %14.3lf sec\n", oplabels[k], times[k]);
		//pred_avx2(&src, 0, depths);
		//rct_JPEG2000_32(&src, 0);
		//compare_bufs_16(dst.data, src.data, src.iw, src.ih, src.nch, src.nch, CODECNAME, 0, 1);
#endif

		array_free(&cdata);
		image_clear(&src);
		image_clear(&dst);
	}
	printf("\nDone.\n");
	pause();
	return 0;
}