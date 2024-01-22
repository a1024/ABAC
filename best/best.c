#include"best.h"
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
static const char file[]=__FILE__;


#define CODECID     47
#define CODECNAME "T47"
#define ENCODE     t47_encode
#define DECODE     t47_decode


	#define BATCHTEST_PRINTTABLE


const char *g_extensions[]=
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
	//int iw, ih;
	//unsigned char *src, *dst;
	size_t usize, csize1, csize2;
	double fdec, enc, dec;//time in secs
	int error;//whether the image was recovered successfully
	ptrdiff_t idx;
} ThreadArgs;
//static void free_threadargs(void *p)
//{
//	ThreadArgs *args=(ThreadArgs*)p;
//	free(args->src);
//	free(args->dst);
//	args->dst=0;
//}
static THREAD_RET THREAD_CALL sample_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	ArrayHandle cdata=0;

	args->usize=(int)ceil(image_getBMPsize(args->src));

	double t=time_sec();
	ENCODE(args->src, &cdata, 0);
	t=time_sec()-t;
	args->enc=t;

	args->csize2=cdata->count;

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
} Result;
typedef struct ProcessCtxStruct
{
	int nstarted, nfinished;
	ArrayHandle threadargs;//<ThreadArgs>	*thread_count
	ArrayHandle results;//<Result>		*nsamples
} ProcessCtx;
static void print_result(Result *res, const char *title, int width)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	printf("%-*s  %10lld  format %10lld %10.6lf D %12lf sec  test %10lld %10.6lf E %12lf D %12lf sec %s\n",
		width, title, res->usize,
		res->csize1, CR1, res->fdec,
		res->csize2, CR2, res->enc, res->dec, res->error?"ERROR":"SUCCESS"
	);
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
		};
		ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
		print_result(&result, (char*)threadargs->title->data, maxlen);

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
			};
			ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
			print_result(&result, (char*)threadargs->title->data, maxlen);
			//print_result(&result, ctx->nfinished+k+1);
		}

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
}
void batch_test_mt(const char *path, int nthreads)
{
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("%s\n", g_buf);
	printf("Multithreaded Batch Test\n");
	double t_start=time_sec();
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
		ArrayHandle title=get_filetitle((char*)fn[0]->data, (int)fn[0]->count);
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

		double t=time_sec();
		Image *image=image_load((char*)fn[0]->data);
		t=time_sec()-t;
		if(!image)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}
		process_file(&processctx, *title, width, image, formatsize, t, k, nthreads);
	}
	process_file(&processctx, 0, width, 0, 0, 0, 0, 0);//set nthreads=0 to flush queued images
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
		print_result(&total, "Total:", width);
		array_free(&processctx.results);
	}
	printf("Batch elapsed ");
	timedelta2str(0, 0, time_sec()-t_start);
	printf("\n");
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);
}

void print_usage(const char *argv0)
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
		" %s  c  <srcfn>  <dstfn.LSIM>      Losslessly encode src image to dst.\n"
		" %s  d  <srcfn.LSIM>  <dstfn.PNG>  Losslessly decode src image to dst as a PNG file.\n"
		" %s  t  <im1>  <im2>               Measure MSE & PSNR between two images.\n",
		argv0, argv0, len, "", argv0, argv0, argv0
	);
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
	ENCODE(src, &cdata, 1);
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
	Image *src=image_load(srcfn);
	if(!src)
	{
		printf("Cannot open \"%s\"\n", srcfn);
		return;
	}
	ArrayHandle cdata=0;
	lsim_writeheader(&cdata, src->iw, src->ih, src->nch, src->src_depth, CODECID);
	ENCODE(src, &cdata, 1);

	int success=save_file(dstfn, cdata->data, cdata->count, 1);
	printf("%s\n", success?"Saved.":"Failed to save.");

	array_free(&cdata);
	free(src);
}
static void decode(const char *srcfn, const char *dstfn)
{
	printf("Decoding \"%s\"\n", srcfn);
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
	OP_TESTFILE, 1, 0,//op, nthreads, formatsize

	"D:/ML/dataset-kodak/kodim13.png",

//	"C:/Projects/datasets/dataset-kodak/kodim02.png",
//	"C:/Projects/datasets/dataset-kodak/kodim13.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/artificial.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/big_building.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/cathedral.png",
//	"C:/Projects/datasets/dataset-ic-rgb16bit/hdr.png",
};
int main(int argc, char **argv)
{
	printf("Best  %s %s\n\n", __DATE__, __TIME__);
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