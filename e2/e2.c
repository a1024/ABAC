#include"e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
static const char file[]=__FILE__;


//override program
//	#define DIFF_AV
//	#define BENCHMARK_DIV
//	#define SORT_TEST
//	#define DSP_TEST


#define CODECID     47
#define CODECNAME "T47"
#define ENCODE     t47_encode
#define DECODE     t47_decode

#if CODECID==47
	#define PRINT_RCT//comment when not applicable
#endif


static const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
	"ppm",
	"pgm",
};
typedef struct ThreadArgsStruct
{
	ArrayHandle title;
	Image *src, *dst;
	size_t usize, csize1, csize2;
	double fdec, enc, dec;//time in secs
#ifdef PRINT_RCT
	SLIC5Curiosity curiosity;
#endif
	int error;//whether the image was recovered successfully
	ptrdiff_t idx;
} ThreadArgs;
static void sample_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	ArrayHandle cdata=0;

	args->usize=(int)ceil(image_getBMPsize(args->src));

	double t=time_sec();
#ifdef PRINT_RCT
	ENCODE(args->src, &cdata, &args->curiosity, 0);
#else
	ENCODE(args->src, &cdata, 0);
#endif
	t=time_sec()-t;
	args->enc=t;

	args->csize2=cdata->count;
	
	//rct_JPEG2000_32(args->src, 1);
	t=time_sec();
	DECODE(cdata->data, cdata->count, args->dst, 0);
	t=time_sec()-t;
	args->dec=t;

	array_free(&cdata);
	args->error=compare_bufs_32(args->dst->data, args->src->data, args->src->iw, args->src->ih, args->src->nch, 4, CODECNAME, 0, 0);
	free(args->src);
	free(args->dst);
	args->dst=0;
}
typedef struct ResultStruct
{
	//ptrdiff_t idx;
	size_t usize, csize1, csize2, error;
	double fdec, enc, dec;
#ifdef PRINT_RCT
	SLIC5Curiosity curiosity;
#endif
} Result;
typedef struct ProcessCtxStruct
{
	int nstarted, nfinished;
	ArrayHandle threadargs;//<ThreadArgs>	*thread_count
	ArrayHandle results;//<Result>		*nsamples
} ProcessCtx;
static double start_time=0, check_time=0;
static double g_total_usize=0, g_total_csize=0;
#ifdef PRINT_RCT
static void curiosity_add(SLIC5Curiosity *dst, SLIC5Curiosity const *src)
{
#ifndef SLIC5_OPTIMIZE_RCT
	for(int k=0;k<RCT_COUNT;++k)
		dst->rct_sizes[k]+=src->rct_sizes[k];
#endif
	for(int k=0;k<SLIC5_NPREDS;++k)
		dst->pred_errors[k]+=src->pred_errors[k];
}
static void curiosity_print(SLIC5Curiosity const *src)
{
#ifndef SLIC5_OPTIMIZE_RCT
	ptrdiff_t usize=0;
	for(int k=0;k<RCT_COUNT;++k)
		usize+=src->coverage[k];
	printf("RCT csizes & coverage:\n");
	for(int k=0;k<RCT_COUNT;++k)
		printf("%20lf %8.4lf%% %s\n", src->rct_sizes[k], 100.*src->coverage[k]/usize, rct_names[k]);
#endif

	printf("Pred errors:\n");
	long long sum=0;
	for(int k=0;k<SLIC5_NPREDS;++k)
		sum+=src->pred_errors[k];
	for(int k=0;k<SLIC5_NPREDS;++k)
		printf("  %2d  %17lld  %6.2lf%%  %s\n", k, src->pred_errors[k], 100.*src->pred_errors[k]/sum, slic5_prednames[k]);
}
#endif
static void print_result(Result *res, const char *title, int width, int print_timestamp, int print_rct)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	g_total_usize+=res->usize;
	g_total_csize+=res->csize2;
	printf("%-*s  %10zd  format %10zd %10.6lf%% D %12lf sec %8.3lf MB/s  test %10zd %10.6lf%% E %12lf D %12lf sec %8.3lf MB/s %s",
		width, title, res->usize,
		res->csize1, 100./CR1, res->fdec, res->usize/(res->fdec*1024*1024),
		res->csize2, 100./CR2, res->enc, res->dec, res->usize/(res->dec*1024*1024),
		res->error?"ERROR":"OK"
	);
	//printf("%-*s  %10zd  format %10zd %10.6lf%% D %12lf sec  test %10zd %10.6lf%% E %12lf D %12lf sec %s  all %10.6lf%%",
	//	width, title, res->usize,
	//	res->csize1, 100./CR1, res->fdec,
	//	res->csize2, 100./CR2, res->enc, res->dec, res->error?"ERROR":"SUCCESS",
	//	100.*g_total_csize/g_total_usize
	//);
#ifdef PRINT_RCT
#ifdef SLIC5_OPTIMIZE_RCT
	if(print_rct)
	{
		printf(" ");
		orct_print_compact(res->curiosity.rct_params);
		//printf(" [%s", slic5_orct_permutationnames[res->curiosity.rct_params[ORCT_NPARAMS]]);
		//for(int k=0;k<ORCT_NPARAMS;++k)
		//{
		//	int val=res->curiosity.rct_params[k];
		//	printf("%c%02X", val<0?'-':'+', abs(val));
		//}
		//	//printf("%02X%c", res->curiosity.rct_params[k]&0xFF, k<_countof(res->curiosity.rct_params)-1?' ':']');
		//printf("]");
	}
#else
	if(print_rct)
	{
		const char *a;
		if((unsigned)res->curiosity.rct<_countof(rct_names))
			a=rct_names[res->curiosity.rct];
		else
			a="ERROR";
		printf(" %-11s", a);
	}
#endif
#else
	(void)print_rct;
#endif
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
static void process_file(ProcessCtx *ctx, ArrayHandle title, int maxlen, Image *image, size_t csize1, double fdec, ptrdiff_t idx, int nthreads, SLIC5Curiosity *curiosity)
{
	if(!ctx->nstarted)
	{
		ARRAY_ALLOC(ThreadArgs, ctx->threadargs, 0, 0, nthreads, 0);
		ARRAY_ALLOC(Result, ctx->results, 0, 0, 0, 0);
	}

	if(image)
	{
		ThreadArgs *threadargs=(ThreadArgs*)ARRAY_APPEND(ctx->threadargs, 0, 1, 1, 0);
		size_t res=(size_t)image->iw*image->ih;
		threadargs->src=image;
		image_copy_nodata(&threadargs->dst, image);
		if(!threadargs->dst)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(threadargs->dst->data, 0, res*sizeof(int[4]));
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
		ThreadArgs *threadargs=array_at(&ctx->threadargs, ctx->threadargs->count-1);
		sample_thread(threadargs);
		Result result=
		{
			threadargs->usize,
			threadargs->csize1,
			threadargs->csize2,
			threadargs->error,
			threadargs->fdec,
			threadargs->enc,
			threadargs->dec,
#ifdef PRINT_RCT
			threadargs->curiosity,
#endif
		};
		ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
		print_result(&result, (char*)threadargs->title->data, maxlen, 1, 1);
#ifdef PRINT_RCT
		curiosity_add(curiosity, &threadargs->curiosity);
		curiosity->coverage[threadargs->curiosity.rct]+=threadargs->usize;
#else
		(void)curiosity;
#endif

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
	else if(ctx->nstarted-ctx->nfinished>=nthreads)
	{
		int n=ctx->nstarted-ctx->nfinished;
		void *hthreads=mt_exec(sample_thread, ctx->threadargs->data, (int)ctx->threadargs->esize, (int)ctx->threadargs->count);
		mt_finish(hthreads);
#if 0
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
#endif
		for(int k=0;k<n;++k)
		{
			ThreadArgs *threadargs=(ThreadArgs*)array_at(&ctx->threadargs, k);
			Result result=
			{
				threadargs->usize,
				threadargs->csize1,
				threadargs->csize2,
				threadargs->error,
				threadargs->fdec,
				threadargs->enc,
				threadargs->dec,
#ifdef PRINT_RCT
				threadargs->curiosity,
#endif
			};
			ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
			print_result(&result, (char*)threadargs->title->data, maxlen, k>=n-1, 1);
			//print_result(&result, ctx->nfinished+k+1);
#ifdef PRINT_RCT
			curiosity_add(curiosity, &threadargs->curiosity);
			curiosity->coverage[threadargs->curiosity.rct]+=threadargs->usize;
#endif
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
	SLIC5Curiosity curiosity={0};
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		ArrayHandle *title=(ArrayHandle*)array_at(&titles, k);

		ptrdiff_t formatsize=get_filesize((char*)fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images, this check is useless because get_filenames() has already filtered the list
			continue;

		double t=time_sec();
		Image *image=image_load((char*)fn[0]->data);
		t=time_sec()-t;
		if(!image)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}
		process_file(&processctx, *title, width, image, formatsize, t, k, nthreads, &curiosity);
	}
	process_file(&processctx, 0, width, 0, 0, 0, 0, 0, &curiosity);//set nthreads=0 to flush queued images
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
		print_result(&total, "Total:", width, 2, 0);
		array_free(&processctx.results);
	}
#ifdef PRINT_RCT
	curiosity_print(&curiosity);
#endif
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);
}

static void print_usage(const char *argv0)
{
	//skip the full path if present, to print only the program title
	int len=(int)strlen(argv0), ks;
	for(ks=len-1;ks>=0&&argv0[ks]!='/'&&argv0[ks]!='\\';--ks);
	++ks;
	len-=ks;
	argv0+=ks;

	printf(
		"Usage:\n"
		" %s  <srcfn>                       Test the current codec.\n"
		" %s  <folder>  [<nthreads>]        Test the current codec on all images in the folder.\n"
		" %*s                                nthreads defaults to 1.\n"
	//	" %*s                                nthreads defaults to the available number of CPU cores.\n"
		" %s  c  <srcfn>  <dstfn.LSIM>      Losslessly encode src image to dst.\n"
		" %s  d  <srcfn.LSIM>  <dstfn.PNG>  Losslessly decode src image to dst as a PNG file.\n"
		" %s  t  <im1>  <im2>               Measure MSE & PSNR between two images.\n",
		argv0, argv0, len, "", argv0, argv0, argv0
	);
	//printf(
	//	"Usage:\n"
	//	" %s  <file_or_path>\n"
	//	" %s  <path> <nthreads>\n"
	//	" %s  \"mse\" <file1> <file2>\n",
	//	argv0, argv0, argv0
	//);
}
typedef enum ProgOpEnum
{
	OP_INVALID,
	OP_TESTFILE,
	OP_TESTFOLDER,
	OP_COMPRESS,
	OP_DECOMPRESS,
	OP_COMPARE,
} ProgOp;
typedef struct ProgArgsStruct
{
	ProgOp op;
	int nthreads;
	ptrdiff_t formatsize;
	const char *fn1, *fn2;
} ProgArgs;
static int parse_cmdargs(int argc, char **argv, ProgArgs *args)
{
	switch(argc)
	{
	case 2://file or folder
		args->formatsize=get_filesize(argv[1]);
		if(args->formatsize<0)//not file nor folder
		{
			args->op=OP_INVALID;
			return 0;
		}
		args->op=args->formatsize?OP_TESTFILE:OP_TESTFOLDER;
		args->nthreads=1;
		//args->nthreads=args->formatsize?1:query_cpu_cores();//default to number of available cores
		args->fn1=argv[1];
		args->fn2=0;
		break;
	case 3://folder nthreads
		args->formatsize=get_filesize(argv[1]);
		if(args->formatsize)//not a folder
		{
			args->op=OP_INVALID;
			return 0;
		}
		args->op=OP_TESTFOLDER;
		args->nthreads=atoi(argv[2]);
		if(args->nthreads<1||args->nthreads>25)
		{
			args->op=OP_INVALID;
			return 0;
		}
		args->fn1=argv[1];
		args->fn2=0;
		break;
	case 4://encode/decode/compare
		if(strlen(argv[1])!=1)
		{
			args->op=OP_INVALID;
			return 0;
		}
		switch(argv[1][0]&0xDF)
		{
		case 'C':args->op=OP_COMPRESS;break;
		case 'D':args->op=OP_DECOMPRESS;break;
		case 'T':args->op=OP_COMPARE;break;
		}
		args->nthreads=1;
		args->fn1=argv[2];
		args->fn2=argv[3];
		break;
	default:
		args->op=OP_INVALID;
		return 0;
	}
	return 1;//valid
}
static void test_one(const char *fn, ptrdiff_t formatsize)
{
	if(!formatsize)
		formatsize=get_filesize(fn);
	printf("Testing \"%s\"\n", fn);
	double t=time_sec();
	Image *src=image_load(fn);
	t=time_sec()-t;
	if(!src)
	{
		printf("Cannot open \"%s\"\n", fn);
		return;
	}
	double usize=image_getBMPsize(src);
	printf("Opened in %lf sec  csize %zd  invCR %lf%%\n", t, formatsize, 100.*formatsize/usize);

	//test_alphaVSbin(src);
#if 1
	//encode
	ArrayHandle cdata=0;
#ifdef PRINT_RCT
	ENCODE(src, &cdata, 0, 1);
#else
	ENCODE(src, &cdata, 1);
#endif
	if(!cdata)
	{
		printf("Encode error\n");
		return;
	}
#if 1
	//decode
	Image *dst=0;
	image_copy_nodata(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	size_t res=(size_t)dst->iw*dst->ih;
	memset(dst->data, 0, res*sizeof(int[4]));
	
	//rct_JPEG2000_32(src, 1);
	DECODE(cdata->data, cdata->count, dst, 1);

	//check & cleanup
	compare_bufs_32(dst->data, src->data, src->iw, src->ih, 3, 4, CODECNAME, 0, 1);
	//memset(dst->data, 0, res*sizeof(int[4]));
#endif
	array_free(&cdata);
	free(dst);
#endif
	printf("\n");

	free(src);
}
static void encode(const char *srcfn, const char *dstfn)
{
	printf("Encoding \"%s\"\n", srcfn);
	int start=0, end=0;
	get_filetitle(srcfn, -1, &start, &end);
	//if(!_stricmp(srcfn+end, ".PPM"))
	//{
	//	t47_from_ppm(srcfn, dstfn);
	//	return;
	//}
	Image *src=image_load(srcfn);
	if(!src)
	{
		printf("Cannot open \"%s\"\n", srcfn);
		return;
	}
	ArrayHandle cdata=0;
	lsim_writeheader(&cdata, src->iw, src->ih, src->nch, src->src_depth, CODECID);
#ifdef PRINT_RCT
	ENCODE(src, &cdata, 0, 1);
#else
	ENCODE(src, &cdata, 1);
#endif

	int success=save_file(dstfn, cdata->data, cdata->count, 1);
	printf("%s\n", success?"Saved.":"Failed to save.");

	array_free(&cdata);
	free(src);
}
static void decode(const char *srcfn, const char *dstfn)
{
	printf("Decoding \"%s\"\n", srcfn);
	int start=0, end=0;
	get_filetitle(dstfn, -1, &start, &end);
	//if(!_stricmp(dstfn+end, ".PPM"))
	//{
	//	t47_to_ppm(srcfn, dstfn);
	//	return;
	//}
	//else
	if(!_stricmp(dstfn+end, ".PNG"))
	{
		ArrayHandle cdata=load_file(srcfn, 1, 16, 0);
		if(!cdata)
		{
			printf("Cannot open \'%s\'", srcfn);
			return;
		}
		LSIMHeader header;
		size_t idx=lsim_readheader(cdata->data, cdata->count, &header);
		Image *image=0;
		image_from_lsimheader(&image, &header);

		int success=DECODE(cdata->data+idx, cdata->count-idx, image, 1);
		array_free(&cdata);
		if(!success)
		{
			printf("Failed to decode \'%s\'", srcfn);
			return;
		}
		success=image_save_native(dstfn, image, !image->depth[3]);
		printf("%s\n", success?"Saved.":"Failed to save.");
		free(image);
	}
}
static void compare(const char *fn1, const char *fn2)
{
	Image *im1=image_load(fn1), *im2=image_load(fn2);
	if(!im1)
	{
		printf("Cannot open %s\n", fn1);
		return;
	}
	if(!im2)
	{
		printf("Cannot open %s\n", fn2);
		return;
	}
	if(im2->iw!=im1->iw||im2->ih!=im1->ih)
	{
		printf("Dimension mismatch  %dx%d != %dx%d\n", im1->iw, im1->ih, im2->iw, im2->ih);
		return;
	}
	ptrdiff_t formatsize1=get_filesize(fn1);
	ptrdiff_t formatsize2=get_filesize(fn2);
	double usize=image_getBMPsize(im1);
	int res=im1->iw*im1->ih;
	long long sum[4]={0};
	for(int k=0;k<res;++k)
	{
		int
			dr=im1->data[k<<2|0]-im2->data[k<<2|0],
			dg=im1->data[k<<2|1]-im2->data[k<<2|1],
			db=im1->data[k<<2|2]-im2->data[k<<2|2],
			da=im1->data[k<<2|3]-im2->data[k<<2|3];
		sum[0]+=dr*dr;
		sum[1]+=dg*dg;
		sum[2]+=db*db;
		sum[3]+=da*da;
	}
	int nch=get_nch32(im1->data, res);
	double rmse[]=
	{
		sqrt((double)sum[0]/res),
		sqrt((double)sum[1]/res),
		sqrt((double)sum[2]/res),
		sqrt((double)sum[3]/res),
		sqrt((double)(sum[0]+sum[1]+sum[2]+sum[3])/(res*nch)),
	};
	double psnr[]=
	{
		20*log10(255/rmse[0]),
		20*log10(255/rmse[1]),
		20*log10(255/rmse[2]),
		20*log10(255/rmse[3]),
		20*log10(255/rmse[4]),
	};
	double CR1=usize/formatsize1;
	double CR2=usize/formatsize2;
	printf("T RMSE %lf PSNR %lf  %lf  L %d %lf%%  R %d %lf%%\n", rmse[4], psnr[4], usize, (int)formatsize1, 100./CR1, (int)formatsize2, 100./CR2);
	printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
	printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
	printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
	if(nch==4)
		printf("A RMSE %lf PSNR %lf\n", rmse[3], psnr[3]);
	return;
}
ProgArgs args=
{
#if 1
	OP_TESTFILE, 1, 0,//op, nthreads, formatsize

//	"D:/ML/dataset-kodak/kodim13.png",
//	"D:/ML/dataset-kodak-ppm/kodim13.ppm",
//	"C:/dataset-LPCB-ppm/canon_eos_1100d_02.ppm",
//	"D:/ML/20240407 blank.PNG",
//	"D:/ML/dataset-ic-rgb16bit/artificial.png",
//	"D:/ML/dataset-ic-rgb16bit/big_building.png",
//	"D:/ML/dataset-ic-rgb16bit/big_tree.png",
//	"D:/ML/dataset-ic-rgb16bit/bridge.png",
//	"D:/ML/dataset-ic-rgb16bit/cathedral.png",
//	"D:/ML/dataset-ic-rgb16bit/deer.png",
//	"D:/ML/dataset-CLIC30/11.png",

//	"C:/Projects/datasets/dataset-kodak/kodim02.png",
//	"C:/Projects/datasets/dataset-kodak/kodim05.png",
	"C:/Projects/datasets/dataset-kodak/kodim13.png",
//	"C:/Projects/datasets/dataset-kodak/kodim19.png",
//	"C:/Projects/datasets/dataset-kodak/kodim21.png",
//	"C:/Projects/datasets/dataset-kodak/kodim22.png",
//	"C:/Projects/datasets/dataset-kodak/kodim23.png",
//	"C:/Projects/datasets/kodim13-small4.PNG",
//	"C:/Projects/datasets/xplogo.jpg",
//	"C:/Projects/datasets/rubberduck.jpg",
//	"C:/Projects/datasets/dataset-train/AWM/world.jpg",//odd dims
//	"C:/Projects/datasets/dataset-kodak-pgm/kodim02.pgm",
//	"C:/Projects/datasets/dataset-kodak-CLIC30/01.png",
//	"C:/Projects/datasets/dataset-kodak-CLIC30/02.png",
//	"C:/Projects/datasets/dataset-kodak-CLIC30/03.png",//prefers RCT_NONE
//	"C:/Projects/datasets/dataset-kodak-CLIC30/05.png",
//	"C:/Projects/datasets/space-HUGE.ppm",
//	"C:/Projects/datasets/Screenshots/Screenshot 2023-04-10 153155.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/artificial.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/big_building.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/cathedral.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/deer.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/hdr.png",
//	"C:/Projects/datasets/dataset-LPCB/canon_eos_1100d_01.PNG",
//	"C:/Projects/datasets/dataset-LPCB/canon_eos_1100d_02.PNG",
//	"C:/Projects/datasets/dataset-LPCB/canon_eos_1100d_03.PNG",
//	"C:/Projects/datasets/dataset-LPCB/PIA13815.PNG",
//	"C:/Projects/datasets/dataset-LPCB/PIA13833.PNG",
//	"C:/Projects/datasets/dataset-LPCB/STA13456.PNG",//prefers RCT_NONE
//	"C:/Projects/datasets/dataset-LPCB/STA13782.PNG",//prefers RCT_NONE
//	"C:/Projects/datasets/dataset-LPCB/STA13942.PNG",
//	"C:/Projects/datasets/Screenshots/Screenshot 2023-03-12 181054.png",

	0,
#else
	OP_TESTFOLDER, 1, 0,
	"/media/awm/Toshiba/ML/dataset-kodak",
	0,

//	OP_COMPRESS, 1, 0,
//	"C:/Projects/datasets/space.ppm",
//	"C:/Projects/datasets/space.lsim",

//	OP_DECOMPRESS, 1, 0,
//	"C:/Projects/datasets/dataset-kodak-ppm/kodim13.lsim",
//	"C:/Projects/datasets/dataset-kodak-ppm/kodim13-dec.ppm",

//	OP_TESTFOLDER, 1, 0,
//	//"C:/Projects/datasets/dataset-LPCB",
//	"C:/Projects/datasets/Screenshots",
#endif
};
//{
//	OP_DECOMPRESS, 1, 0,
//	"C:/Projects/datasets/kodim13.lsim",
//	"C:/Projects/datasets/kodim13-dec.png",
//};

#ifdef DSP_TEST
#include"lodepng.h"
typedef struct ACPStruct
{
	double alpha;
	int x, y;
} ACP;
#define DSP_NNWW	-2, -2
#define DSP_NNW		-1, -2
#define DSP_NN		 0, -2
#define DSP_NNE		 1, -2
#define DSP_NNEE	 2, -2
#define DSP_NWW		-2, -1
#define DSP_NW		-1, -1
#define DSP_N		 0, -1
#define DSP_NE		 1, -1
#define DSP_NEE		 2, -1
#define DSP_NEEE	 3, -1
#define DSP_WW		-2,  0
#define DSP_W		-1,  0
ACP dsp_info[]=
{//	alpha, x, y
//	{0.125,		DSP_NN	},
//	{0.998,		DSP_NW	},
//	{0.7499,	DSP_N	},
//	{0.499,		DSP_NE	},
//	{0.3333,	DSP_NEE	},
//	{0.25,		DSP_NEEE},
//	{0.125,		DSP_WW	},
//	{0.001,		DSP_W	},

	{0.25,		DSP_NEEE},
	{0.749999,	DSP_W	},
	
//	{0.3333,	DSP_NEE	},
//	{0.6665,	DSP_W	},

//	{0.499,		DSP_NE	},
//	{0.5,		DSP_W	},
	
//	{0.8749,	DSP_N	},//X
//	{0.125,		DSP_W	},
};
#define DSP_GAMMA 1
//#define DSP_ALPHA1 0		//W
//#define DSP_ALPHA2 0.331	//N
//#define DSP_ALPHA3 0.331	//NE
//#define DSP_ALPHA4 0.331	//NEE
#define DSP_REACH 256//controls buffer dimensions
#define DSP_WIDTH (DSP_REACH<<1|1)
#define DSP_HEIGHT (DSP_REACH+1)
//#define DSP_NITER 30
double g_image[DSP_HEIGHT][DSP_WIDTH];
//char g_mask[DSP_HEIGHT][DSP_WIDTH];
#endif
#ifdef SORT_TEST
void sort_bitonic(short *data, int lgcount)
{
	int count=1<<lgcount;
	for(int it=0;it<lgcount;++it)
	{
		for(int k2=it>>1;k2>0;k2>>=1)
		{
			for(int k3=0;k3<count;++k3)
			{
				int k4=k3^k2;
				short temp;
				if(k4>k3&&!(k3&it)!=(data[k3]<data[k4]))
					SWAPVAR(data[k3], data[k4], temp);
			}
		}
	}
}
static int threeway_short(const void *p1, const void *p2)
{
	const short
		*s1=(const short*)p1,
		*s2=(const short*)p2;
	return (*s1>*s2)-(*s1<*s2);
}
void init_rand(short *data, int count)
{
	for(int k=0;k<count;++k)
		data[k]=rand();
}
void sort_radix_nlsb(short *data, int count, int nsteps, short *temp)
{
	if(nsteps<1||nsteps>16||(nsteps&(nsteps-1)))//nsteps in {1, 2, 4, 8, 16}
		return;
	int nbits=sizeof(short[8])/nsteps, histsize=1<<nbits, mask=histsize-1;
	int *hist=(int*)malloc(histsize*sizeof(int));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kb=0;kb<nsteps;++kb)
	{
		int shift=nbits*kb, val;
		memset(hist, 0, histsize*sizeof(int));
		for(int k=0;k<count;++k)
		{
			val=data[k]>>shift&mask;
			++hist[val];
		}
		for(int k=0, sum=0;k<histsize;++k)
		{
			val=hist[k];
			hist[k]=sum;
			sum+=val;
		}
		for(int k=0;k<count;++k)
		{
			val=data[k]>>shift&mask;
			temp[hist[val]]=data[k];
			++hist[val];
		}
		memcpy(data, temp, count*sizeof(short));
	}
	free(hist);
}
void sort_radix_lsd(short *data, int count, short *temp)
{
	int hist[16];
	for(int kb=0;kb<16;kb+=4)
	{
		memset(hist, 0, sizeof(hist));
		for(int k=0;k<count;++k)
		{
			int val=data[k]>>kb&15;
			++hist[val];
		}
		for(int k=0, sum=0;k<16;++k)
		{
			int freq=hist[k];
			hist[k]=sum;
			sum+=freq;
		}
		for(int k=0;k<count;++k)
		{
			int val=data[k]>>kb&15;
			temp[hist[val]]=data[k];
			++hist[val];
		}
		memcpy(data, temp, count*sizeof(short));
	}
}
void sort_radix_msd(short *data, int count, short *temp, int kb)
{
	int hist[16]={0}, hist2[16];
	for(int k=0;k<count;++k)
	{
		int val=data[k]>>(12-kb)&15;
		++hist[val];
	}
	for(int k=0, sum=0;k<16;++k)
	{
		int freq=hist[k];
		hist[k]=sum;
		sum+=freq;
	}
	memcpy(hist2, hist, sizeof(hist2));
	for(int k=0;k<count;++k)
	{
		int val=data[k]>>(12-kb)&15;
		temp[hist[val]]=data[k];
		++hist[val];
	}
	memcpy(data, temp, count*sizeof(short));
	if(kb>=12)
		return;
	for(int k=0;k<16;++k)
	{
		int pos=hist2[k], count2=(k<15?hist2[k+1]:count)-pos;
		if(count2>0)
			sort_radix_msd(data+pos, count2, temp, kb+4);
	}
}
#endif
#ifdef DIFF_AV
#include"lodepng.h"
//unsigned char buf[32*32*4];
#endif
int main(int argc, char **argv)
{
#ifdef DIFF_AV
	printf("Diff-Av\n");
	Image *src=image_load("C:/Projects/e2/20240301 6 src.PNG");
	int srcdim=MINVAR(src->iw, src->ih), dstdim=srcdim<<1, half=srcdim>>1;
	Image *dst=image_alloc(dstdim, dstdim, 3, 8, 8, 8, 0, 1, -1);
	int a, b, temp;
	for(int ky=0;ky<srcdim;++ky)
	{
		for(int kx=0;kx<srcdim;++kx)
		{
			a=kx-half, b=ky-half;

			SWAPVAR(a, b, temp);
			b=-b;

			a-=b;
			b+=a>>1;

			a=-a;

			memcpy(dst->data+((dst->iw*((size_t)b+srcdim)+a+srcdim)<<2), src->data+(((size_t)src->iw*ky+kx)<<2), sizeof(int[4]));
		}
	}
	image_snapshot(dst);
#if 0
	memset(buf, -1, sizeof(buf));
	for(int ky=0;ky<16;++ky)
	{
		for(int kx=0;kx<16;++kx)
		{
			int a=kx, b=ky;
			a-=8;
			b-=8;

			a-=b;
			b+=a>>1;

			//a-=b;
			//a=((a+8)&15)-8;
			//b+=a>>1;
			//b=((b+8)&15)-8;
			//
			//a+=8;
			//b+=8;
			//a&=15;
			//b&=15;

			//a-=8;
			//b-=8;
			//
			//b-=a>>1;
			//b=((b+8)&15)-8;
			//a+=b;
			//a=((a+8)&15)-8;
			//
			//a+=8;
			//b+=8;
			//a&=15;
			//b&=15;
			printf(" %X%X", (a+8)&15, (b+8)&15);
			buf[((b+16)<<5|(a+16))<<2|0]=kx<<4;
			buf[((b+16)<<5|(a+16))<<2|1]=ky<<4;
			buf[((b+16)<<5|(a+16))<<2|2]=0;
			buf[((b+16)<<5|(a+16))<<2|3]=0xFF;
		}
		printf("\n");
	}
	lodepng_encode_file("free-dst.PNG", buf, 32, 32, LCT_RGBA, 8);
#endif
	printf("Done.\n");
	pause();
	return 0;
#endif
#ifdef BENCHMARK_DIV
	{
#define BENCHMARK_SIZE 0x1000
#define BENCHMARK_COUNT 0x1000
		printf("benchmark\n");
		double *ddata=(double*)malloc(sizeof(double[BENCHMARK_SIZE]));
		float *fdata=(float*)malloc(sizeof(float[BENCHMARK_SIZE]));
		long long *ldata=(long long*)malloc(sizeof(long long[BENCHMARK_SIZE]));
		int *idata=(int*)malloc(sizeof(int[BENCHMARK_SIZE]));
		if(!ddata||!fdata||!ldata||!idata)
		{
			LOG_ERROR("Buy more RAM");
			return 0;
		}
		for(int k=0;k<BENCHMARK_SIZE;++k)
		{
			ddata[k]=(double)rand()/256;
			fdata[k]=(float)rand()/256;
			ldata[k]=(long long)rand()<<8;
			idata[k]=rand();
		}
		volatile double t[16]={0}, ttemp;
		volatile long long c[16]={0}, ctemp;

#define BENCHLOOP(IDX, DATA, OP)\
	ttemp=time_sec();\
	ctemp=__rdtsc();\
	for(int k=0;k<BENCHMARK_SIZE-1;++k)\
		DATA[k] OP (DATA[k+1]);\
	t[IDX]+=time_sec()-ttemp;\
	c[IDX]+=__rdtsc()-ctemp;
		
		for(int k=0;k<BENCHMARK_COUNT;++k)
		{
			BENCHLOOP( 0, idata, =abs)
			BENCHLOOP( 1, fdata, =fabsf)
			BENCHLOOP( 2, ldata, =llabs)
			BENCHLOOP( 3, ddata, =fabs)
			BENCHLOOP( 4, idata, =1+)
			BENCHLOOP( 5, fdata, =1+)
			BENCHLOOP( 6, ldata, =1+)
			BENCHLOOP( 7, ddata, =1+)
			BENCHLOOP( 8, idata, /=)
			BENCHLOOP( 9, fdata, /=)
			BENCHLOOP(10, ldata, /=)
			BENCHLOOP(11, ddata, /=)
			BENCHLOOP(12, idata, *=)
			BENCHLOOP(13, fdata, *=)
			BENCHLOOP(14, ldata, *=)
			BENCHLOOP(15, ddata, *=)
			printf("%6.2lf%%\r", 100.*(k+1)/BENCHMARK_COUNT);
		}
		printf("\n");

		double dsum=0;
		float fsum=0;
		int isum=0;
		long long lsum=0;
		for(int k=0;k<BENCHMARK_SIZE;++k)
		{
			dsum+=ddata[k];
			fsum+=fdata[k];
			isum+=idata[k];
			lsum+=ldata[k];
		}
		const char *data_labels[]=
		{
			"i32",
			"f32",
			"i64",
			"f64",
		};
		const char *op_labels[]=
		{
			"abs",
			"inc",
			"div",
			"mul",
		};
		double maxtime=0;
		for(int k=0;k<16;++k)
		{
			UPDATE_MAX(maxtime, t[k]);
		}
		for(int k=0;k<16;++k)
		{
			printf("%s %s %10lf ms  %15lld cycles  ", data_labels[k&3], op_labels[k>>2], t[k]*1000, c[k]);
			for(int k2=0, nstars=(int)(t[k]*60/maxtime);k2<nstars;++k2)
				printf("*");
			printf("\n");
		}
		printf("Result  %d %f %lld %lf\n", isum, fsum, lsum, dsum);

		free(ddata);
		free(fdata);
		free(ldata);
		free(idata);

		printf("Done.\n");
		pause();
		return 0;
	}
#endif
#ifdef SORT_TEST
	{
#define SORT_COUNT 0x10000
		for(int k=0;k<32;++k)
			printf("%d\n", rand()&0xFF);
		short *orig=(short*)malloc(sizeof(short[SORT_COUNT]));
		short *sort1=(short*)malloc(sizeof(short[SORT_COUNT]));
		//short *sort2=(short*)malloc(sizeof(short[SORT_COUNT]));
		short *sort3=(short*)malloc(sizeof(short[SORT_COUNT]));
		short *temp=(short*)malloc(sizeof(short[SORT_COUNT]));
		if(!orig||!sort1||!sort3)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		init_rand(orig, SORT_COUNT);
		memcpy(sort1, orig, sizeof(short[SORT_COUNT]));
		//memcpy(sort2, orig, sizeof(short[SORT_COUNT]));
		memcpy(sort3, orig, sizeof(short[SORT_COUNT]));

		printf("%d elements\n", SORT_COUNT);

		double t=time_sec();
		long long c=__rdtsc();
		qsort(sort1, SORT_COUNT, sizeof(short), threeway_short);
		c=__rdtsc()-c;
		t=time_sec()-t;
		printf("qsort  %lf sec  %10lld cycles\n", t, c);
		
		//t=time_sec();
		//sort_radix_lsd(sort2, SORT_COUNT, temp);
		//t=time_sec()-t;
		//int success=!memcmp(sort2, sort1, sizeof(short[SORT_COUNT]));
		//printf("sort_radix_lsd      %lf  %s\n", t, success?"OK":"ERROR");

		//t=time_sec();
		//sort_radix_msd(sort3, SORT_COUNT, temp, 0);
		//t=time_sec()-t;
		//success=!memcmp(sort3, sort1, sizeof(short[SORT_COUNT]));
		//printf("sort_radix_msd      %lf  %s\n", t, success?"OK":"ERROR");

		for(int k=0;k<=4;++k)
		{
			memcpy(sort3, orig, sizeof(short[SORT_COUNT]));

			t=time_sec();
			c=__rdtsc();
			sort_radix_nlsb(sort3, SORT_COUNT, 1<<k, temp);
			c=__rdtsc()-c;
			t=time_sec()-t;
			int success=!memcmp(sort3, sort1, sizeof(short[SORT_COUNT]));
			printf("sort_radix_nlsb(%2d-bit)  %lf sec  %10lld cycles  %s\n", (int)(sizeof(short[8])/(1LL<<k)), t, c, success?"OK":"ERROR");
			
			//for(int k=0;k<SORT_COUNT;++k)
			//	printf("[%4d] 0x%08X 0x%08X\n", k, sort1[k], sort3[k]);
			//printf("\n");
		}
		//for(int k=0;k<SORT_COUNT;++k)
		//	printf("[%4d] 0x%08X 0x%08X 0x%08X\n", k, sort1[k], sort2[k], sort3[k]);

		free(orig);
		free(sort1);
		//free(sort2);
		free(sort3);
		free(temp);
		pause();
		return 0;
	}
#endif
#if 0
	int sum=0;
	double t=time_sec();
	for(int k=0;k<1000000000;++k)
	{
		unsigned n=rand()<<15|rand();
		int lgn=FLOOR_LOG2(n);
		sum+=lgn;
	}
	printf("%lf  %d\n", time_sec()-t, sum);
	pause();
	return 0;
#endif
#ifdef DSP_TEST
	double asum=0;
	for(int k=0;k<_countof(dsp_info);++k)
		asum+=dsp_info[k].alpha;
	if(asum>1)
		printf("ALPHA SUM %lf\n", asum);
	//for(int k=0;k<DSP_HEIGHT*DSP_WIDTH;++k)
	//	((double*)g_image)[k]=-1;
	memset(g_image, 0, sizeof(g_image));
	//memset(g_mask, 0, sizeof(g_mask));
	g_image[DSP_HEIGHT-1][DSP_REACH]=0x10000;
	//g_mask[DSP_HEIGHT-1][DSP_REACH]=1;
	//double sum=0;
	for(int ky=DSP_HEIGHT-1;ky>=0;--ky)//the kernel is causal, so 1 iteration in reverse is enough
	{
		for(int kx=DSP_WIDTH-1;kx>=0;--kx)
		{
			double *curr=g_image[ky]+kx;
//#define LOAD(X, Y) ((unsigned)(kx+(X))<(unsigned)DSP_W&&(unsigned)(ky+(Y))<(unsigned)DSP_H?g_image[ky-(Y)]+kx-(X):0)
//#define UPDATE_NB(NB, ALPHA) if(NB)*NB=(*NB>=0?*NB:0)+*curr*ALPHA
			for(int k=0;k<_countof(dsp_info);++k)
			{
				int kx2=kx+dsp_info[k].x, ky2=ky+dsp_info[k].y;
				if((unsigned)kx2<(unsigned)DSP_WIDTH&&(unsigned)ky2<(unsigned)DSP_HEIGHT)
					g_image[ky2][kx2]+=*curr*dsp_info[k].alpha;
			}
			*curr*=1-asum;
			//double *nb1=LOAD( 1,  0);//W
			//double *nb2=LOAD( 0,  1);//N
			//double *nb3=LOAD(-1,  1);//NE
			//double *nb4=LOAD(-2,  1);//NEE
			//if(*curr>=0&&(nb1&&*nb1<0||nb2&&*nb2<0||nb3&&*nb3<0||nb4&&*nb4<0))
			//{
			//	UPDATE_NB(nb1, DSP_ALPHA1);
			//	UPDATE_NB(nb2, DSP_ALPHA2);
			//	UPDATE_NB(nb3, DSP_ALPHA3);
			//	UPDATE_NB(nb4, DSP_ALPHA4);
			//	*curr*=1-asum;
			//	updated=1;
			//}
//#undef  LOAD
//#undef  UPDATE_NB
		}
	}
	//sum=0;
	//for(int ky=0;ky<DSP_H;++ky)
	//{
	//	for(int kx=0;kx<DSP_W;++kx)
	//	{
	//		double val=g_image[ky][kx];
	//		char c=kx<DSP_W-1?' ':'\n';
	//		if(val<0)
	//			printf("%s%c", ky==DSP_H-1&&kx==DSP_REACH+1?" [-] ":"  -  ", c);
	//		else
	//		{
	//			printf("%05d%c", (int)round(val*100000/0x10000), c);
	//			//printf("%5.1lf%c", val, c);
	//			sum+=val;
	//		}
	//	}
	//}
	//printf("[NE  W  curr].[%lf  %lf  %lf]  wsum %lf\n\n",
	//	(double)DSP_ALPHA1, (double)DSP_ALPHA2, 1.-DSP_ALPHA1-DSP_ALPHA2, sum/0x10000
	//);

	unsigned char *result=(unsigned char*)malloc(DSP_WIDTH*DSP_HEIGHT*sizeof(char[4]));
	if(!result)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(result, 0, DSP_WIDTH*DSP_HEIGHT*sizeof(char[4]));
	double vmin=INFINITY, vmax=-INFINITY, wsum=0;
	for(int ky=0;ky<DSP_HEIGHT;++ky)
	{
		for(int kx=0;kx<DSP_WIDTH;++kx)
		{
			double val=g_image[ky][kx];
			if(val)
			{
				//if(val<1e-6)//this clamp makes sum > 1
				//	val=1e-6;
				//double lg_val=log2(val/0x10000);
				UPDATE_MIN(vmin, val);
				UPDATE_MAX(vmax, val);
				wsum+=val;
			}
		}
	}
	for(int ky=0;ky<DSP_HEIGHT;++ky)
	{
		for(int kx=0;kx<DSP_WIDTH;++kx)
		{
			double val=g_image[ky][kx];
			int idx=(DSP_WIDTH*ky+kx)<<2;
			//if(ky==DSP_HEIGHT-2&&kx==DSP_REACH+5)//
			//	printf("");
			if(fabs(val)<1e-9)
			{
				if(ky==DSP_HEIGHT-1&&kx==DSP_REACH+1)//mark the target blue
				{
					result[idx|0]=0;
					result[idx|1]=0;
					result[idx|2]=0xFF;
				}
				else//unused is yellow
				{
					result[idx|0]=255;
					result[idx|1]=255;
					result[idx|2]=0;
				}
			}
			else
			{
				val=(val-vmin)/(vmax-vmin);
				val=pow(val, DSP_GAMMA);
				val*=255;
				//if(val<1e-6)//this clamp makes sum > 1
				//	val=1e-6;
				//val=255*(log2(val/0x10000)-vmin)/(vmax-vmin);
				unsigned char level=(unsigned char)val;
				result[idx|0]=level;
				result[idx|1]=level;
				result[idx|2]=level;
			}
			result[idx|3]=0xFF;
		}
	}
	int printed=snprintf(g_buf, G_BUF_SIZE, "%s", argv[0]);
	for(;printed>=0&&g_buf[printed-1]!='/'&&g_buf[printed-1]!='\\';--printed);
	double west=g_image[DSP_HEIGHT-1][DSP_REACH];
	printed+=acme_strftime(g_buf+printed, G_BUF_SIZE-printed, "%Y%m%d_%H-%M-%S.PNG");
	printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "_W%.2lfpercent.PNG", 100.*west/wsum);
	//printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "-[%g %g %g %g].PNG",
	//	(double)DSP_ALPHA1, (double)DSP_ALPHA2, (double)DSP_ALPHA3, (double)DSP_ALPHA4
	//);
	printf("wsum  %lf (west %6.2lf%%  rest %6.2lf%%)  About to save \"%s\"\n", wsum/0x10000, 100.*west/wsum, 100.*(wsum-west)/wsum, g_buf);
	pause();
	lodepng_encode_file(g_buf, result, DSP_WIDTH, DSP_HEIGHT, LCT_RGBA, 8);
	return 0;
#endif
#if 0
	printf("FIXED PREC MATH TEST\n");
	printf("x\tlgx\texact\tdiff\n");
	for(int k=1;k<=256;++k)
	{
		unsigned long long x=(unsigned long long)k<<24;
		
		if(x==0x3000000)//
			printf("");

		int lgx=log2_fix24(x);
		double trial=(double)lgx/0x1000000;
		double truth=log2((double)x/0x1000000);
		printf("0x%016llX>>24  %16lf  %16lf  %16lf\n", x, trial, truth, truth-trial);
		//printf("0x%016llX>>24  %c0x%08X>>24  %c0x%08X>>24\n", x, lgx<0?'-':'+', abs(lgx), true_lgx<0?'-':'+', abs(true_lgx));
	}
	pause();
	return 0;
#endif
#if 0
	printf("FIXED PREC MATH TEST\n");
	for(;;)
	{
#if 0
#define POW_FIX24(B, E) exp2_fix24((int)((long long)E*log2_fix24(B)>>24))
		int base=0, e=0;

		printf("Enter base (hex fix24 uint32 no prefix): ");
		while(!scanf("%X", &base));
		base=abs(base);
		printf("Enter exponent (hex fix24 int32 no prefix): ");
		while(!scanf("%X", &e));
		//base=0x800000, e=0x2000000;

		unsigned long long trial=POW_FIX24(base, e);
		//int temp=log2_fix24(base);
		//unsigned long long trial=exp2_fix24((int)((long long)e*temp>>24));
		unsigned long long truth=(unsigned long long)(pow((double)base/0x1000000, (double)e/0x1000000)*0x1000000);
		printf("  TRIAL 0x%016llX>>24\n  TRUTH 0x%016llX>>24\n", trial, truth);
#endif

		//int k=0;
		//printf("Enter an int32 in fix24 (hex without prefix): ");
		//while(!scanf("%X", &k));
		//unsigned long long trial=exp2_fix24(k);
		//unsigned long long truth=(unsigned long long)(exp2((double)k/0x1000000)*0x1000000);
		//printf("exp2(%c0x%08X>>24)  TRIAL 0x%016llX>>24  TRUTH 0x%016llX>>24\n", k<0?'-':'+', abs(k), trial, truth);


		//unsigned long long k=0;
		//printf("Enter a uint64 in fix24 (hex without prefix): ");
		//while(!scanf("%llX", &k));
		//int trial=log2_fix24(k);
		//int truth=(int)(log2((double)k/0x1000000)*0x1000000);
		//printf("log2(0x%016llX>>24)  TRIAL %c0x%08X>>24  TRUTH %c0x%08X>>24\n", k, trial<0?'-':'+', abs(trial), truth<0?'-':'+', abs(truth));

		unsigned long long k=0;
		printf("Enter a uint63: ");
		while(!scanf("%llu", &k));
		unsigned s=floor_sqrt(k);
		printf("SQRT(%llu) = %u\n", k, s);
	}
#endif
#if 0
	//t47_analyze_preds("C:/Projects/datasets/dataset-LPCB");
	t47_analyze_preds("C:/Projects/datasets/dataset-CLIC30");
	//t47_analyze_preds("C:/Projects/datasets/temp");
	return 0;
#endif
#if 0
	{
		const char
			*src="C:/Projects/datasets/dataset-kodak-ppm/kodim13.ppm",
			*lsim="C:/Projects/datasets/dataset-kodak-ppm/kodim13.lsim",
			*dec="C:/Projects/datasets/dataset-kodak-ppm/kodim13-dec.ppm";
		t47_from_ppm(src, lsim);
		t47_to_ppm(lsim, dec);

		printf("Done.\n");
		pause();
		return 0;
	}
#endif
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
	//system("cd");
	//init_vk();//
	//test1();
	//void ac_vs_ans();
	//ac_vs_ans();
	//test_calendar();

	printf("Entropy2\n\n");
#ifndef _DEBUG
	parse_cmdargs(argc, argv, &args);
#endif
	switch(args.op)
	{
	case OP_INVALID:
		print_usage(argv[0]);
		break;
	case OP_TESTFILE:
		test_one(args.fn1, args.formatsize);
		break;
	case OP_TESTFOLDER:
		batch_test_mt(args.fn1, args.nthreads);
		break;
	case OP_COMPRESS:
		encode(args.fn1, args.fn2);
		break;
	case OP_DECOMPRESS:
		decode(args.fn1, args.fn2);
		break;
	case OP_COMPARE:
		compare(args.fn1, args.fn2);
		break;
	}
	printf("Done.\n");
	pause();
	return 0;
}