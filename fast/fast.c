﻿#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
static const char file[]=__FILE__;


//	#define _DEBUG
//	#define DNGCODEC
//	#define BENCH_QOI

#define CODECID     23
#define CODECNAME "F23"
#define ENCODE     f23_encode
#define DECODE     f23_decode


typedef struct ThreadArgsStruct
{
	ArrayHandle title, ext;
	Image src, dst;
	size_t usize, csize1, csize2;
	double fdec, enc, dec;//time in secs
	int error, unused;//whether the image was recovered successfully
	ptrdiff_t idx;
} ThreadArgs;
static void sample_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
#ifdef BENCH_QOI
	args->usize=image_getBMPsize(&args->src);
	qoi_test(&args->src, &args->csize2, &args->enc, &args->dec, &args->error, 0);
#else
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
#endif
	image_clear(&args->src);
	image_clear(&args->dst);
}
typedef struct ResultStruct
{
	size_t usize, csize1, csize2, error;
	double fdec, enc, dec;
	ptrdiff_t idx;
} Result;
typedef struct ProcessCtxStruct
{
	int nstarted, nfinished;
	ArrayHandle threadargs;//<ThreadArgs>	*thread_count
	ArrayHandle results;//<Result>		*nsamples
} ProcessCtx;
static double start_time=0, check_time=0;
static double g_total_usize=0, g_total_csize=0;
static void print_result(Result *res, const char *title, const char *ext, int twidth, int ewidth, int print_timestamp)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	g_total_usize+=res->usize;
	g_total_csize+=res->csize2;
	printf("%-*s  %10zd  %-*s %10zd %10.6lf%% D %12lf sec %9.4lf MB/s %5td %10zd %10.6lf%% E %12lf D %12lf sec %9.4lf MB/s %s",
		twidth, title,
		res->usize,
		ewidth, ext,
		res->csize1, 100./CR1, res->fdec, res->usize/(res->fdec*1024*1024),
		res->idx+1,
		res->csize2, 100./CR2, res->enc, res->dec, res->usize/(res->dec*1024*1024),
		res->error?"ERROR":"OK"
	);
	//printf("%-*s  %10zd  format %10zd %10.6lf%% D %12lf sec  test %10zd %10.6lf%% E %12lf D %12lf sec %s  all %10.6lf%%",
	//	width, title, res->usize,
	//	res->csize1, 100./CR1, res->fdec,
	//	res->csize2, 100./CR2, res->enc, res->dec, res->error?"ERROR":"OK",
	//	100.*g_total_csize/g_total_usize
	//);
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
static void process_file(ProcessCtx *ctx, ArrayHandle title, ArrayHandle ext, int twidth, int ewidth, Image *image, size_t csize1, double fdec, ptrdiff_t idx, int nthreads)
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
		threadargs->ext=ext;
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
			threadargs->idx,
		};
		ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
		print_result(&result, (char*)threadargs->title->data, (char*)threadargs->ext->data, twidth, ewidth, 1);

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
	else if(ctx->nstarted-ctx->nfinished>=nthreads)
	{
		int n=ctx->nstarted-ctx->nfinished;

		void *hthreads=mt_exec(sample_thread, ctx->threadargs->data, (int)ctx->threadargs->esize, (int)ctx->threadargs->count);
		mt_finish(hthreads);

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
				threadargs->idx,
			};
			ARRAY_APPEND(ctx->results, &result, 1, 1, 0);
			print_result(&result, (char*)threadargs->title->data, (char*)threadargs->ext->data, twidth, ewidth, k>=n-1);
		}

		array_clear(&ctx->threadargs);
		ctx->nfinished=ctx->nstarted;
	}
}
static void batch_test_mt(const char *path, int nthreads)
{
	static const char *ext[]=
	{
#ifdef DNGCODEC
		"dng",
#else
		"png",
		"jpg", "jpeg",
		"ppm", "pgm",
		"bmp",
		"tif", "tiff",
#endif
	};
	
	ArrayHandle filenames, titles, exts;
	ProcessCtx processctx={0};
	double t_start;
	int twidth=6;//"Total:"
	int ewidth=0;

	if(CODECID>=23)//codecs starting from F23 are multithreaded in shaa Allah
		nthreads=1;
	g_total_usize=0;
	g_total_csize=0;
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("%s\n", g_buf);
	printf("Multithreaded Batch Test %s  \"%s\"\n", CODECNAME, path);
	filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames)
	{
		printf("No supported images in \"%s\"\n", path);
		return;
	}
	ARRAY_ALLOC(ArrayHandle, titles, 0, 0, filenames->count, (void(*)(void*))array_free);
	ARRAY_ALLOC(ArrayHandle, exts, 0, 0, filenames->count, (void(*)(void*))array_free);
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle title;
		ArrayHandle ext;
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		int start=0, end=0;
		get_filetitle((char*)fn[0]->data, (int)fn[0]->count, &start, &end);
		STR_COPY(title, (char*)fn[0]->data+start, end-start);
		STR_COPY(ext, (char*)fn[0]->data+end, fn[0]->count-end);
		ARRAY_APPEND(titles, &title, 1, 1, 0);
		ARRAY_APPEND(exts, &ext, 1, 1, 0);
		if(twidth<(int)title->count)
			twidth=(int)title->count;
		if(ewidth<(int)ext->count)
			ewidth=(int)ext->count;
	}
	check_time=start_time=t_start=time_sec();
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		ArrayHandle *title=(ArrayHandle*)array_at(&titles, k);
		ArrayHandle *ext2=(ArrayHandle*)array_at(&exts, k);

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
		process_file(&processctx, *title, *ext2, twidth, ewidth, &image, formatsize, t, k, nthreads);
		//printf("");
	}
	process_file(&processctx, 0, 0, twidth, ewidth, 0, 0, 0, 0, 0);//set nthreads=0 to flush queued images
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
			total.idx=result->idx;
		}
		printf("\n");
		print_result(&total, "Total:", "", twidth, ewidth, 2);
		array_free(&processctx.results);
	}
#if CODECID==24
	f24_curiosity();
#elif CODECID==26
	f26_curiosity();
#elif CODECID==27
	f27_curiosity();
#endif
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);
	array_free(&titles);
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
	//printf("FastEntropy\n\n");
#ifndef _DEBUG
//#if 0
	if((unsigned)(argc-2)>1)
	{
		printf(
			"Usage:\n"
			"  %s  filename                        Test codec without saving\n"
			"  %s  folder [nthreads]               Batch test without saving\n"
#ifdef DNGCODEC
			"  %s  input.DNG  output.LSIM          Encode image\n"
#else
			"  %s  input.PPM/PGM/PNG  output.LSIM  Encode image\n"
#endif
			"  %s  input.LSIM  output.PPM/PGM/PNG  Decode image\n",
			argv[0], argv[0], argv[0], argv[0]
		);
#ifndef __GNUC__
		pause();
#endif
		return 1;
	}
	const char *fn=argv[1], *arg2=argc==3?argv[2]:0;
#else
	const char *arg2=
		0
	//	"D:/ML/kodim13.lsim.ppm"
	//	"D:/ML/big_building.LSIM.PPM"
		;
	const char *fn=
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
	//	"D:/ML/kodim13.lsim"
	//	"D:/ML/dataset-kodak-ppm/kodim20.ppm"
	//	"D:/ML/dataset-kodak/kodim13.png"
		"D:/ML/dataset-CLIC30-ppm/03.ppm"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/big_building.LSIM"
	//	"D:/ML/20240407 blank.ppm"
	//	"D:/ML/PIA13757-crop.PNG"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_03.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_04.ppm"
	//	"C:/dataset-LPCB-ppm/fujifilm_finepix_x100_01.ppm"
	//	"C:/dataset-LPCB-ppm/olympus_xz1_14.ppm"
	//	"C:/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/dataset-LPCB-ppm/PIA12813.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13757.ppm"	//MER-3D
	//	"C:/dataset-LPCB-ppm/PIA13785.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13943.ppm"	//mostly blank
	//	"C:/dataset-LPCB-ppm/PIA13799.ppm"	//mars crater, needs palette?
	//	"C:/dataset-LPCB-ppm/PIA13882.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13912.ppm"
	//	"C:/dataset-LPCB-ppm/STA13456.ppm"	//uncorrelated channels
	//	"C:/dataset-LPCB-ppm/STA13844.ppm"	//space gas clouds, 6800^2
	//	"C:/dataset-LPCB-ppm"
	//	"D:/ML/dataset-RAW/a0001-jmac_DSC1459.dng"
	//	"D:/ML/dataset-CID22-ppm/pexels-photo-2802032.PPM"
	//	"D:/ML/dataset-CID22-ppm"
	//	"D:/ML/dataset-kodak-small"

	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim13.ppm"
	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim24.ppm"	//borderless
	//	"C:/Projects/datasets/dataset-kodak-pgm/kodim13.pgm"
	//	"C:/Projects/datasets/dataset-kodak/kodim13.png"
	//	"C:/Projects/datasets/kodim13-small4.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-CLIC16-ppm/2048x1320_adam-willoughby-knox-56406.ppm"
	//	"C:/Projects/datasets/dataset-CLIC16-ppm/2048x1320_eric-huang-35182.ppm"
	//	"C:/Projects/datasets/kodim13-small4.PNG"
	//	"C:/Projects/datasets/PNG_transparency_demonstration_1.png"
	//	"C:/Projects/datasets/dataset-DNG/L1020006.DNG"
	//	"C:/Projects/datasets/dataset-DNG"
	//	"C:/Projects/datasets/dataset-CID22-ppm/3637739.PPM"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_02.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13891.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13943.ppm"	//highly skewed
	//	"C:/Projects/datasets/dataset-ic-rgb16bit/flower_foveon.png"	//smallest 16-bit image
	//	"C:/Projects/datasets/dataset-ic-rgb16bit/deer.png"
	//	"C:/Projects/datasets/jupiter.PNG"		//actually 8-bit
	//	"C:/Projects/datasets/space-8k-CROP.PPM"
	//	"C:/Projects/datasets/dataset-CID22-ppm"
	//	"C:/Projects/datasets/dataset-CLIC30"
	//	"C:/Projects/datasets/dataset-LPCB-ppm"
		;
#endif
	ptrdiff_t formatsize=get_filesize(fn);
	if(arg2)
	{
		const char *ext1=strrchr(fn, '.'), *ext2=strrchr(arg2, '.');
		if(ext2)
		{
			if((!strcmp_ci(ext1, ".PPM")||!strcmp_ci(ext1, ".PGM")||!strcmp_ci(ext1, ".PNG"))&&!strcmp_ci(ext2, ".LSIM"))
			{
				Image src={0};
				image_load(fn, &src);
				if(!src.data)
				{
					printf("Canot open \"%s\"\n", fn);
					return 1;
				}
				
				ArrayHandle cdata=0;
				lsim_writeheader(&cdata, src.iw, src.ih, src.nch, src.depth, CODECID);
				ENCODE(&src, &cdata, 0);
				{
					int success=save_file(arg2, cdata->data, cdata->count, 1);
					if(!success)
						printf("Failed to save\n");
					//printf("%s\n", success?"Saved.":"Failed to save.");
				}

				array_free(&cdata);
				image_clear(&src);
				return 0;
			}
			if(!strcmp_ci(ext1, ".LSIM")&&(!strcmp_ci(ext2, ".PPM")||!strcmp_ci(ext2, ".PGM")||!strcmp_ci(ext2, ".PNG")))
			{
				int e, success;
				size_t idx;
				LSIMHeader header;
				Image dst={0};
				ArrayHandle cdata=load_file(fn, 1, 16, 0);
				if(!cdata)
				{
					printf("Cannot open \'%s\'", fn);
					return 1;
				}
				idx=lsim_readheader(cdata->data, cdata->count, &header);
				image_from_lsimheader(&dst, &header);

				e=DECODE(cdata->data+idx, cdata->count-idx, &dst, 0);
				array_free(&cdata);
				if(e)
				{
					printf("Failed to decode \'%s\'", fn);
					return 1;
				}
				if(!strcmp_ci(ext2, ".PNG"))
					success=image_save_native(arg2, &dst);
				else
					success=image_save_ppm(arg2, &dst);
				if(!success)
					printf("Failed to save\n");
				//printf("%s\n", success?"Saved.":"Failed to save.");
				image_clear(&dst);
				return 0;
			}
		}
	}
	if(formatsize<0)
	{
		printf("Not a file nor directory:  \'%s\'\n", fn);
#ifndef __GNUC__
		pause();
#endif
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

		//qoi_test(&src, 0, 0, 0, 0, 1);
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
#ifndef __GNUC__
	pause();
#endif
	return 0;
}