#include"e2.h"
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
#include"lodepng.h"
static const char file[]=__FILE__;


//override program
//	#define DSP_TEST


#define CODECID     47
#define CODECNAME "T47"
#define ENCODE     t47_encode
#define DECODE     t47_decode

	#define PRINT_RCT//comment when not applicable


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
	SLIC5Curiosity curiosity;
	int error;//whether the image was recovered successfully
	ptrdiff_t idx;
} ThreadArgs;
static THREAD_RET THREAD_CALL sample_thread(void *param)
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
	return 0;
}
typedef struct ResultStruct
{
	//ptrdiff_t idx;
	size_t usize, csize1, csize2, error;
	double fdec, enc, dec;
	SLIC5Curiosity curiosity;
} Result;
typedef struct ProcessCtxStruct
{
	int nstarted, nfinished;
	ArrayHandle threadargs;//<ThreadArgs>	*thread_count
	ArrayHandle results;//<Result>		*nsamples
} ProcessCtx;
static double start_time=0;
static double g_total_usize=0, g_total_csize=0;
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
	printf("RCT csizes:\n");
	for(int k=0;k<RCT_COUNT;++k)
		printf("%20lf %s\n", src->rct_sizes[k], rct_names[k]);
#endif

	printf("Pred errors:\n");
	long long sum=0;
	for(int k=0;k<SLIC5_NPREDS;++k)
		sum+=src->pred_errors[k];
	for(int k=0;k<SLIC5_NPREDS;++k)
		printf("  %2d  %17lld  %6.2lf%%  %s\n", k, src->pred_errors[k], 100.*src->pred_errors[k]/sum, slic5_prednames[k]);
}
static void print_result(Result *res, const char *title, int width, int print_timestamp, int print_rct)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	g_total_usize+=res->usize;
	g_total_csize+=res->csize2;
	printf("%-*s  %10lld  format %10lld %10.6lf D %12lf sec  test %10lld %10.6lf E %12lf D %12lf sec %s  CCR %10.6lf",
		width, title, res->usize,
		res->csize1, CR1, res->fdec,
		res->csize2, CR2, res->enc, res->dec, res->error?"ERROR":"SUCCESS",
		g_total_usize/g_total_csize
	);
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
#endif
	if(print_timestamp)
	{
		printf(" ");
		timedelta2str(0, 0, time_sec()-start_time);
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
			threadargs->curiosity,
		};
		ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
		print_result(&result, (char*)threadargs->title->data, maxlen, 1, 1);
		curiosity_add(curiosity, &threadargs->curiosity);

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
				threadargs->error,
				threadargs->fdec,
				threadargs->enc,
				threadargs->dec,
				threadargs->curiosity,
			};
			ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
			print_result(&result, (char*)threadargs->title->data, maxlen, k>=n-1, 1);
			//print_result(&result, ctx->nfinished+k+1);
			curiosity_add(curiosity, &threadargs->curiosity);
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
	start_time=t_start;
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
		size_t start=0, end=0;
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
		for(int k=0;k<processctx.results->count;++k)
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
		print_result(&total, "Total:", width, 1, 0);
		array_free(&processctx.results);
	}
	curiosity_print(&curiosity);
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
	printf("Opened in %lf sec  csize %lld  CR %lf\n", t, formatsize, usize/formatsize);

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
	printf("\n");

	array_free(&cdata);
	free(src);
	free(dst);
}
static void encode(const char *srcfn, const char *dstfn)
{
	printf("Encoding \"%s\"\n", srcfn);
	size_t start=0, end=0;
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
	size_t start=0, end=0;
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
	ptrdiff_t formatsize=get_filesize(fn2);
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
	double CR=res*3./formatsize;
	printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[4], psnr[4], res*nch, (int)formatsize, CR, 8/CR);
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
//	"D:/ML/dataset-ic-rgb16bit/artificial.png",
//	"D:/ML/dataset-ic-rgb16bit/big_building.png",
//	"D:/ML/dataset-ic-rgb16bit/big_tree.png",
//	"D:/ML/dataset-ic-rgb16bit/bridge.png",
//	"D:/ML/dataset-ic-rgb16bit/cathedral.png",
//	"D:/ML/dataset-ic-rgb16bit/deer.png",

//	"C:/Projects/datasets/dataset-kodak/kodim02.png",
//	"C:/Projects/datasets/dataset-kodak/kodim05.png",
	"C:/Projects/datasets/dataset-kodak/kodim13.png",
//	"C:/Projects/datasets/dataset-kodak/kodim22.png",
//	"C:/Projects/datasets/kodim13-small4.PNG",
//	"C:/Projects/datasets/dataset-kodak/kodim19.png",
//	"C:/Projects/datasets/dataset-kodak/kodim21.png",
//	"C:/Projects/datasets/dataset-kodak-pgm/kodim02.pgm",
//	"C:/Projects/datasets/dataset-kodak-CLIC30/01.png",
//	"C:/Projects/datasets/dataset-kodak-CLIC30/03.png",//prefers RCT_NONE
//	"C:/Projects/datasets/dataset-kodak-CLIC30/05.png",
//	"C:/Projects/datasets/Screenshots/Screenshot 2023-04-10 153155.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/artificial.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/big_building.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/cathedral.png",
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
#else
	OP_COMPRESS, 1, 0,
	"C:/Projects/datasets/space.ppm",
	"C:/Projects/datasets/space.lsim",

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
int main(int argc, char **argv)
{
#if 0
	int sum=0;
	double t=time_sec();
	for(int k=0;k<1000000000;++k)
	{
		unsigned n=rand()<<15|rand();
		int lgn=floor_log2_32(n);
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


		//printf("Enter a uint64: ");
		//while(!scanf("%llu", &k));
		//unsigned s=floor_sqrt(k);
		//printf("SQRT(%llu) = %u\n", k, s);
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

#if 0
	long long cycles;
	int iw=0, ih=0, nch0=3, nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	const char *fn=0;
	int nthreads=0;
#ifdef _DEBUG
	nthreads=10;
	//fn="C:/Projects/datasets/CLIC11-crop4-2.PNG";
	//fn="C:/Projects/datasets/CLIC11-small4.PNG";
	//fn="C:/Projects/datasets/dataset-CLIC30/11.png";
	//fn="C:/Projects/datasets/dataset-kodak";
	//fn="C:/Projects/datasets/dataset-kodak-pgm/kodim01.pgm";
	//fn="C:/Projects/datasets/kodim13-small4.PNG";
	//fn="C:/Projects/datasets/dataset-ic-rgb8bit/big_building.PNG";
	fn="C:/Projects/datasets/dataset-kodak/kodim13.png";

//	fn="D:/ML/dataset-LPCB";
	//fn="D:/ML/dataset-CLIC30";
	//fn="D:/ML/dataset-kodak";
	//fn="D:/ML/dataset-CLIC30/16.png";//hardest noiseless CLIC30 image
	//fn="D:/ML/dataset-CLIC30/17.png";
	//fn="D:/ML/dataset-kodak/kodim13.png";
//	fn="D:/ML/dataset-kodak/kodim18.png";
	//fn="D:/ML/dataset-kodak-small/13.PNG";
#endif
	if(fn||argc==2||argc==3)
	{
		if(!fn)
			fn=argv[1];
		ptrdiff_t formatsize=get_filesize(fn);
		if(formatsize==-1)
		{
			LOG_ERROR("Cannot open \"%s\"", fn);
			print_usage(argv[0]);
			return 0;
		}
		if(!formatsize)//path
		{
			if(!nthreads)
			{
				nthreads=1;
				if(argc==3)
					nthreads=atoi(argv[2]);
			}
			if(nthreads==1)
				batch_test(fn);
			else
				batch_test_mt(fn, nthreads);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		double t_dec=time_sec();
		cycles=__rdtsc();
		buf=image_load(fn, &iw, &ih);
		cycles=__rdtsc()-cycles;
		t_dec=time_sec()-t_dec;
		if(!buf)
		{
			LOG_ERROR("Couldn't open \"%s\"", fn);
			return 0;
		}
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		printf("Format Dec %lfsec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", t_dec, (double)cycles/(resolution*nch0), iw, ih, nch0, (long long)formatsize, (double)resolution*nch0/formatsize);
	}
	else if(argc==4)
	{
		if(!_stricmp(argv[1], "MSE"))
		{
			const char *fn1=argv[2], *fn2=argv[3];
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
			print_usage(argv[0]);
			pause();
			return 0;
		}
	}
	else
	{
		print_usage(argv[0]);
		pause();
		return 0;
	}

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<(int)len;k+=nch)
			buf[k]=0xFF;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	//size_t usize=len*nch0>>2;

	printf("\n");
	
	ArrayHandle cdata=0;
	//const void *ptr, *end;

	//{
	//	size_t bestsize=0, lastsize=0;
	//	int besttest=0;
	//	for(int k=0;k<4;++k)
	//	{
	//		size_t size=ma_test(buf, iw, ih, k, 1, 1);
	//		if(!k||bestsize>size)
	//			bestsize=size, besttest=k;
	//		if(k==7)
	//			lastsize=size;
	//	}
	//	printf("\nBest combination  %d  %lld  %lf\n", besttest, bestsize, (double)lastsize/bestsize);
	//	print_ma_test(besttest);
	//	printf("\n");
	//	//pause();
	//	exit(0);
	//}
	//grad_explore(buf, iw, ih);

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
					compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0, 1);
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
				compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0, 1);
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
		elapsed=time_sec();
		cycles=__rdtsc();
		t25_encode(buf, iw, ih, blockw, blockh, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_sec()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_sec();
		cycles=__rdtsc();
		t25_decode(cdata->data, cdata->count, iw, ih, blockw, blockh, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_sec()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T25", 0, 1);
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
		
		elapsed=time_sec();
		cycles=__rdtsc();
		t26_encode(buf, iw, ih, params, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_sec()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_sec();
		cycles=__rdtsc();
		t26_decode(cdata->data, cdata->count, iw, ih, params, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_sec()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T26", 0, 1);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//T27+: ABAC
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

	//old record
#if 0
	printf("T35 Entropy coding with context tree\n");
	//printf("T35 Combines spatial transform with entropy coding\n");
	t35_encode(buf, iw, ih, &cdata, 1);
	t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T35", 0, 1);
	memset(b2, 0, len);
	printf("\n");
#endif

	//printf("T36 stretch & squish\n");
	//t36_encode(buf, iw, ih, &cdata, 1);
	//t36_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T36", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

	//printf("T37 Fixed array as binary tree predictor\n");
	//t37_encode(buf, iw, ih, &cdata, 1);
	//t37_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T37", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

	//printf("T38 Single simple bit predictor\n");
	//t38_encode(buf, iw, ih, &cdata, 1);
	//t38_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T38", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

	//old record
#if 0
	printf("T39 Multiple estimators for all maps\n");
	t39_encode(buf, iw, ih, &cdata, 1);
	t39_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T39", 0, 1);
	memset(b2, 0, len);
	printf("\n");
#endif

	//t40_encode(buf, iw, ih, &cdata, 2);	//X
	//t40_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T40", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

	//t41_encode(buf, iw, ih, &cdata, 2);	//X
	//t41_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T41", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

#if 0
	t42_encode(buf, iw, ih, &cdata, 1);		//prev record
	t42_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T42", 0, 1);
	memset(b2, 0, len);
	printf("\n");
#endif

	//t43_encode(buf, iw, ih, &cdata, 2);
	//t43_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T43", 0, 1);
	//memset(b2, 0, len);
	//printf("\n");

#if 0
	t44_encode(buf, iw, ih, &cdata, 1);		//current record: from paq8pxd
	t44_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T44", 0, 1);
	memset(b2, 0, len);
	printf("\n");
#endif

#if 1
	t45_encode(buf, iw, ih, &cdata, 1);
	t45_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, CODECNAME, 0, 1);
	memset(b2, 0, len);
	printf("\n");
#endif
#endif

	//SLIC
#if 0
	{
		double t_enc=0, t_dec=0;

		t_enc=time_sec();
		slic_encode(iw, ih, 4, 8, buf, &cdata, 1);
		t_enc=time_sec()-t_enc;

		int iw2=0, ih2=0, nch2=0, depth2=0;
		t_dec=time_sec();
		unsigned char *ret=(unsigned char*)slic_decode(cdata->data, (int)cdata->count, &iw2, &ih2, &nch2, &depth2);
		t_dec=time_sec()-t_dec;

		printf("\nSLI %8d  CR %lf    Enc %lfsec  Dec %lfsec\n", (int)cdata->count, iw*ih*3./cdata->count, t_enc, t_dec);
		array_free(&cdata);

		compare_bufs_uint8(ret, buf, iw, ih, 3, 4, "SLIC", 0, 1);
		free(ret);
	}
#endif

	//SLIC2
#if 0
	{
		double t_enc=0, t_dec=0;

		t_enc=time_sec();
		int retlen=0;
		unsigned char *data=slic2_encode(iw, ih, 4, 8, buf, &retlen);
		t_enc=time_sec()-t_enc;

		int iw2=0, ih2=0, nch2=0, depth2=0;
		t_dec=time_sec();
		unsigned char *ret=(unsigned char*)slic2_decode(data, retlen, &iw2, &ih2, &nch2, &depth2, 0, 1);
		t_dec=time_sec()-t_dec;

		printf("\nSLI %8d  CR %lf    Enc %lfsec  Dec %lfsec\n", (int)retlen, iw*ih*3./retlen, t_enc, t_dec);
		compare_bufs_uint8(ret, buf, iw, ih, 3, 4, "SLI2", 0, 1);
		free(data);
		free(ret);
	}
#endif

	//SIF by Marcio Pais
#if 0
	{
		SIF_content_descriptor_t info=
		{
			iw, ih, 3,//dims, nch
			0,//flags

			//2<<2,
		};
		size_t csize=0;
		memcpy(b2, buf, (size_t)iw*ih<<2);
		pack3_fwd(b2, iw*ih);

		double t_enc=time_sec();
		void *cdata=SIF_compressImage(&info, b2, (size_t)iw*ih<<2, &csize);
		t_enc=time_sec()-t_enc;

		uint64_t usize=0;
		double t_dec=time_sec();
		void *udata=SIF_decompressImage(&info, cdata, csize, &usize);
		t_dec=time_sec()-t_dec;

		printf("SIF  %lld  CR %lf  enc %lf dec %lf  ", csize, (double)iw*ih*4/csize, t_enc, t_dec);

		compare_bufs_uint8(b2, (unsigned char*)udata, iw, ih, 3, 3, "SIF", 0, 1);

		printf("\n");

		SIF_FREE(cdata);
		SIF_FREE(udata);
	}
#endif

	//QOI
#if 0
	{
		for(int k=0, res=iw*ih;k<res;++k)
		{
			b2[k*3+0]=buf[k<<2|0];
			b2[k*3+1]=buf[k<<2|1];
			b2[k*3+2]=buf[k<<2|2];
		}
		qoi_desc desc={iw, ih, 3, QOI_LINEAR};
		int csize=0;
		double t0=time_sec();
		void *cdata_qoi=qoi_encode(b2, &desc, &csize);
		t0=time_sec()-t0;
		printf("\nQOI %8d  CR %lf    Enc %lfsec\n\n", csize, iw*ih*3./csize, t0);
		free(cdata_qoi);
	}
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
		//compare_bufs_uint8((unsigned char*)b3, buf, iw, 1, 1, 1, "squeeze row", 0, 1);
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
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "transform", 0, 1);
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