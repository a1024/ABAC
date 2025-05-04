#include"best.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;


#define CODECLIST\
	CODEC(   39, T39, "T39: 8-bit only, very efficient/slow, ABAC")\
	CODEC(   42, T42, "T42: 8-bit only, very efficient/slow, ABAC")\
	CODEC(   44, T44, "T44: 8-bit only, most efficient/slowest, paq8px ported to C")\
	CODEC(   45, T45, "T45: 8-bit only, kind of efficient/fast, CALIC clone, supports color, fast, AC")\
	CODEC(   46, T46, "T46: [8~16] bit, moderate efficiency/speed, ANS")\
	CODEC(   47, T47, "T47: [8~16] bit, moderate efficiency/speed, AC")\
	CODEC(   54, T54, "T54: [8~16] bit, fast/inefficient, AC")\
	CODEC( 1000, F02, "F02: [8~16] bypass")\
	CODEC( 1023, F23, "F23: [8~16] bit, extremely fast/inefficient, multithreaded, GR")
typedef enum CodecChoiceEnum
{
#define CODEC(ID, NAME, DESC) CODEC_##NAME,
	CODECLIST
#undef  CODEC
	CODEC_COUNT,
} CodecChoice;
static const char *codecnames[]=
{
#define CODEC(ID, NAME, DESC) #NAME,
	CODECLIST
#undef  CODEC
};
static const int codecIDs[]=
{
#define CODEC(ID, NAME, DESC) ID,
	CODECLIST
#undef  CODEC
};
static const char *codecdesc[]=
{
#define CODEC(ID, NAME, DESC) DESC,
	CODECLIST
#undef  CODEC
};
static int cli_select_codec(void)
{
	int choice=0;
	printf("Input   Codec\n");
	for(int k=0;k<_countof(codecIDs);++k)
		printf("%3d\t%s\n", k, codecdesc[k]);
	//printf(
	//	"Input  Codec\n"
	//	"  0     MEMCPY\n"
	//	"  1     T39: 8-bit only, very efficient/slow, ABAC\n"
	//	"  2     T42: 8-bit only, most efficient/slow, ABAC\n"
	//	"  3     T45: 8-bit only, kind of efficient/slow, CALIC clone, supports color, fast, AC\n"
	//	"  4     T46: [8~16] bit, moderate efficiency/slow, ANS\n"
	//	"  5     T47: [8~16] bit, moderate efficiency/slow, AC\n"
	//	"  6     T54: [8~16] bit, fast/inefficient, AC\n"
	//	"  7     F23: [8~16] bit, extremely fast/inefficient, multithreaded, GR\n"
	//);
	for(;;)
	{
		printf("Choose a codec: ");
		//printf("Enter 0 for T42 (better/slower, 8-bit only), other keys for T46 (faster, supports high bit depth): ");
		while(!scanf("%d", &choice));
		if((unsigned)choice<_countof(codecIDs))
			break;
		printf("Invalid choice \"%d\"\n\n", choice);
	}
	printf("\n");
	return choice;
}
static int encode(const Image *src, const unsigned char *src8, int iw, int ih, ArrayHandle *cdata, int codecid, int loud)
{
	switch(codecid)
	{
	case CODEC_T39:return t39_encode(src8, iw, ih, cdata, loud);
	case CODEC_T42:return t42_encode(src8, iw, ih, cdata, loud);
	case CODEC_T44:return t44_encode(src8, iw, ih, cdata, loud);
	case CODEC_T45:return t45_encode(src8, iw, ih, cdata, loud);
	case CODEC_T46:return t46_encode(src, cdata, loud);
	case CODEC_T47:return t47_encode(src, cdata, loud);
	case CODEC_T54:return t54_encode(src, cdata, loud);
	case CODEC_F02:return f02_encode(src, cdata, loud);
	case CODEC_F23:return f23_encode(src, cdata, loud);
	}
	LOG_ERROR("Unsupported codec ID");
	return 1;
}
static int decode(const unsigned char *cdata, size_t clen, int iw, int ih, unsigned char *dst8, Image *dst, int codecid, int loud)
{
	switch(codecid)
	{
	case CODEC_T39:return t39_decode(cdata, clen, iw, ih, dst8, loud);
	case CODEC_T42:return t42_decode(cdata, clen, iw, ih, dst8, loud);
	case CODEC_T44:return t44_decode(cdata, clen, iw, ih, dst8, loud);
	case CODEC_T45:return t45_decode(cdata, clen, iw, ih, dst8, loud);
	case CODEC_T46:return t46_decode(cdata, clen, dst, loud);
	case CODEC_T47:return t47_decode(cdata, clen, dst, loud);
	case CODEC_T54:return t54_decode(cdata, clen, dst, loud);
	case CODEC_F02:return f02_decode(cdata, clen, dst, loud);
	case CODEC_F23:return f23_decode(cdata, clen, dst, loud);
	}
	LOG_ERROR("Unsupported codec ID");
	return 1;
}
static int codec_is_8bit(int codecid)
{
	return codecid==CODEC_T39||codecid==CODEC_T42||codecid==CODEC_T44||codecid==CODEC_T45;
}

typedef enum CodecModeEnum
{
	MODE_TEST,//enc, dec, compare
	MODE_ENC,//PNG->LSIM
	MODE_DEC,//LSIM->PNG
} CodecMode;
typedef struct ThreadArgsStruct
{
	Image *src, *dst;
	
	unsigned char *src8, *dst8;
	int iw, ih, nch, depth;

	ArrayHandle cdata;
	const char *dstfn;

	int mode, loud, codecid, error;
	ptrdiff_t idx, usize, csize1, csize2;
	double fdec, enc, dec;
} ThreadArgs;
static int open_fancy(const char *fn, int _8bit_codec, ThreadArgs *arg, int point_idx, ptrdiff_t formatsize)
{
	volatile double t;

	t=time_sec();
	if(_8bit_codec)
	{
		extern unsigned char* stbi_load(char const *filename, int *x, int *y, int *channels_in_file, int desired_channels);

		arg->src8=stbi_load(fn, &arg->iw, &arg->ih, &arg->nch, 4);
		if(!arg->src8)
		{
			if(arg->loud)
				printf("Cannot open \"%s\"\n", fn);
			return 1;
		}
		arg->depth=8;
		arg->usize=(ptrdiff_t)arg->iw*arg->ih*arg->nch;
	}
	else
	{
		arg->src=image_load(fn);
		if(!arg->src)
		{
			if(arg->loud)
				printf("Cannot open \"%s\"\n", fn);
			return 1;
		}
		arg->iw=arg->src->iw;
		arg->ih=arg->src->ih;
		arg->nch=arg->src->nch;
		arg->depth=arg->src->src_depth[0];
		arg->usize=image_getBMPsize(arg->src);
	}
	t=time_sec()-t;
	arg->fdec=t;
	
	if(arg->loud)
	{
		printf("Opened %s in %lf sec  %lf MB/s  CWHD %d*%d*%d*%d  %td/%td bytes  %lf%%  CR %lf\n",
			fn+point_idx,
			t, arg->usize/(t*1024*1024),
			arg->nch, arg->iw, arg->ih, arg->depth,
			formatsize, arg->usize,
			100.*formatsize/arg->usize, (double)arg->usize/formatsize
		);
	}
	return 0;
}
static void process_sample(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	volatile double t;
	ptrdiff_t datastart=0;
	int _8bit_codec;

	if(args->mode==MODE_TEST||args->mode==MODE_DEC)
	{
		LSIMHeader header={0};
		if(args->mode==MODE_DEC)
		{
			datastart=lsim_readheader(args->cdata->data, args->cdata->count, &header);
			args->iw=header.iw;
			args->ih=header.ih;
			args->nch=header.nch;
			args->depth=header.depth[0];

			args->codecid=-1;
			for(int k=0;k<_countof(codecIDs);++k)
			{
				if(header.codec_id==codecIDs[k])
				{
					args->codecid=k;
					break;
				}
			}
			if(args->codecid<0)
			{
				LOG_ERROR("Unsupported codec ID  %d", header.codec_id);
				return;
			}
		}
		_8bit_codec=codec_is_8bit(args->codecid);
		if(_8bit_codec)
		{
			args->usize=sizeof(char[4])*args->iw*args->ih;
			args->dst8=(unsigned char*)malloc(args->usize);
			if(!args->dst8)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			memset(args->dst8, 0, args->usize);
			args->usize=(ptrdiff_t)args->nch*args->iw*args->ih;
		}
		else
		{
			if(args->mode==MODE_DEC)
				image_from_lsimheader(&args->dst, &header);
			else
			{
				args->usize=image_getBMPsize(args->src);
				image_copy_nodata(&args->dst, args->src);
				if(!args->dst)
				{
					LOG_ERROR("Alloc error");
					return;
				}
			}
		}
	}
	else
	{
		char depths[4]={0};
		memfill(depths, &args->depth, sizeof(char)*args->nch, sizeof(char));
		lsim_writeheader(&args->cdata, args->iw, args->ih, args->nch, depths, codecIDs[args->codecid]);
		_8bit_codec=codec_is_8bit(args->codecid);
	}

	if(args->mode==MODE_TEST||args->mode==MODE_ENC)
	{
		t=time_sec();
		encode(args->src, args->src8, args->iw, args->ih, &args->cdata, args->codecid, args->loud);
		t=time_sec()-t;
		args->enc=t;
	}

	if(!args->cdata)
	{
		LOG_ERROR("Encode error");
		return;
	}
	args->csize2=args->cdata->count;
	
	if(args->mode==MODE_TEST||args->mode==MODE_DEC)
	{
		t=time_sec();
		decode(args->cdata->data+datastart, args->cdata->count-datastart, args->iw, args->ih, args->dst8, args->dst, args->codecid, args->loud);
		t=time_sec()-t;
		args->dec=t;
	}
	
	switch(args->mode)
	{
	case MODE_TEST:
		if(_8bit_codec)
			args->error=compare_bufs_uint8(args->dst8, args->src8, args->iw, args->ih, args->nch, 4, codecnames[args->codecid], 0, args->loud);
		else
			args->error=compare_bufs_32(args->dst->data, args->src->data, args->src->iw, args->src->ih, args->src->nch, 4, codecnames[args->codecid], 0, args->loud);
		//args->error=compare(args->image, dst, args->image8, dst8, args->iw, args->ih, args->nch, args->codecid, args->loud);
		break;
	case MODE_ENC:
		{
			int success=save_file(args->dstfn, args->cdata->data, args->cdata->count, 1);
			if(args->loud)
				printf("%s  \"%s\".\n", success?"Saved":"Failed to save", args->dstfn);
		}
		break;
	case MODE_DEC:
		{
			int error=0;
			if(_8bit_codec)
				error=image_save_buf8(args->dstfn, args->dst8, args->iw, args->ih, args->nch);
			else
				error=!image_save_native(args->dstfn, args->dst, !args->dst->depth[3]);
			if(args->loud)
				printf("%s  \"%s\".\n", error?"Failed to save":"Saved", args->dstfn);
		}
		break;
	}

	array_free(&args->cdata);
	if(_8bit_codec)
	{
		free(args->src8);
		free(args->dst8);
	}
	else
	{
		free(args->src);
		free(args->dst);
	}
}
static void print_result(ThreadArgs *res, const char *title, int titlecolumn, int print_timestamp, double elapsed)
{
	double
		CR1=(double)res->usize/res->csize1,
		CR2=(double)res->usize/res->csize2;
	printf(
		"%-*s  %10zd  format %10zd %10.6lf%% %12lf sec %8.3lf MB/s  test %10zd %10.6lf%% %12lf %12lf sec  %8.3lf %8.3lf MB/s %s",
		titlecolumn, title, res->usize,
		res->csize1, 100./CR1, res->fdec,
		res->usize/(res->fdec*1024*1024),
		res->csize2, 100./CR2, res->enc, res->dec,
		res->usize/(res->enc*1024*1024),
		res->usize/(res->dec*1024*1024),
		res->error?"ERROR":"OK"
	);
	if(print_timestamp)
	{
		printf(" ");
		timedelta2str(0, 0, elapsed);
	}
	printf("\n");
}
static void batch_test(const char *path, int nthreads)
{
	static const char *ext[]=
	{
		"PNG",
		"JPG", "JPEG",
		"PPM", "PGM",
		"BMP",
		"TIF", "TIFF",
	};

	int codecid, _8bit_codec;
	ArrayHandle filenames, titles, threadargs;
	int titlecolumn;
	volatile double t_start;
	ThreadArgs total={0};

	//print date & time
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("%s\n", g_buf);

	//read directory
	filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames)
	{
		printf("No supported images in \"%s\"\n", path);
		return;
	}
	ARRAY_ALLOC(ArrayHandle, titles, 0, 0, filenames->count, (void(*)(void*))array_free);
	titlecolumn=6;//"Total:"
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle title;
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		int start=0, end=0;
		get_filetitle((char*)fn[0]->data, (int)fn[0]->count, &start, &end);
		STR_COPY(title, (char*)fn[0]->data+start, end-start);
		if(titlecolumn<(int)title->count)
			titlecolumn=(int)title->count;
		ARRAY_APPEND(titles, &title, 1, 1, 0);
	}

	//select codec
	codecid=cli_select_codec();
	_8bit_codec=codec_is_8bit(codecid);
	printf("Multithreaded Batch Test  %s\n", codecnames[codecid]);
	if(codecid==CODEC_F23&&nthreads!=1)
	{
		printf("F23 is already multithreaded, defaulting to 1 thread\n");
		nthreads=1;
	}
	ARRAY_ALLOC(ThreadArgs, threadargs, 0, 0, nthreads, 0);

	t_start=time_sec();
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn0=(ArrayHandle*)array_at(&filenames, k);
		const char *fn=(char*)fn0[0]->data;
		ThreadArgs arg={0};
		
		arg.csize1=get_filesize(fn);

		arg.error=open_fancy(fn, _8bit_codec, &arg, 0, arg.csize1);
		if(arg.error)
			continue;
		
		arg.mode=MODE_TEST;
		arg.loud=0;
		arg.codecid=codecid;
		arg.idx=k;
		ARRAY_APPEND(threadargs, &arg, 1, 1, 0);
		if((int)threadargs->count>=nthreads||k+1>=(int)filenames->count)
		{
			volatile double t;

			t=time_sec();
			if(threadargs->count==1)
				process_sample(threadargs->data);
			else
			{
				void *hthreads=mt_exec(process_sample, threadargs->data, (int)threadargs->esize, (int)threadargs->count);
				mt_finish(hthreads);
			}
			t=time_sec()-t;
			for(int k2=0;k2<(int)threadargs->count;++k2)
			{
				ThreadArgs *res=(ThreadArgs*)array_at(&threadargs, k2);
				ArrayHandle *title=(ArrayHandle*)array_at(&titles, res->idx);
				print_result(res, (char*)title[0]->data, titlecolumn, k2>=(int)threadargs->count-1, t);

				total.usize+=res->usize;
				total.csize1+=res->csize1;
				total.csize2+=res->csize2;
				total.fdec+=res->fdec;
				total.enc+=res->enc;
				total.dec+=res->dec;
			}
			array_clear(&threadargs);
		}
	}
	t_start=time_sec()-t_start;
	printf("\n");
	print_result(&total, "Total:", titlecolumn, 1, t_start);
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("%s elapsed %lf sec\n", g_buf, t_start);
	array_free(&filenames);
	array_free(&titles);
	array_free(&threadargs);
}

static void print_usage(const char *argv0)
{
	//skip the full path if present, to print only the program title
	int start=0, end=0;
	get_filetitle(argv0, -1, &start, &end);
	argv0+=start;
	end-=start;
	printf("Usage:\n");
	printf(" %.*s  srcfn                     Test the current codec.\n",					end, argv0);
	printf(" %.*s  folder  [nthreads]        Test the current codec on all images in the folder.\n",	end, argv0);
	printf(" %.*s  c/e  srcfn  dstfn.LSIM    Losslessly encode src image to dst.\n",			end, argv0);
	printf(" %.*s  d  srcfn.LSIM  dstfn.PNG  Losslessly decode src image to dst as a PNG file.\n",		end, argv0);
	printf(" %.*s  t  im1  im2               Measure MSE & PSNR between two images.\n",			end, argv0);

#if 0
	printf(
		"Usage:\n"
		" %.*s  srcfn                     Test the current codec.\n"
		" %.*s  folder  [nthreads]        Test the current codec on all images in the folder.\n"
		 " %*s                                nthreads defaults to 1.\n"
		" %.*s  c  srcfn  dstfn.LSIM      Losslessly encode src image to dst.\n"
		" %.*s  d  srcfn.LSIM  dstfn.PNG  Losslessly decode src image to dst as a PNG file.\n"
		" %.*s  t  im1  im2               Measure MSE & PSNR between two images.\n",
		end, argv0,
		end, argv0,
		end, "",
		end, argv0,
		end, argv0,
		end, argv0
	);
#endif
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
		case 'E':args->op=OP_COMPRESS;break;
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

static void test_one_file(const char *fn, ptrdiff_t formatsize)
{
	ThreadArgs arg={0};
	int codecid, _8bit_codec;
	int start=0, end=0;

	get_filetitle(fn, -1, &start, &end);

	if(!formatsize)
		formatsize=get_filesize(fn);
	if(formatsize<1)
	{
		printf("Not a valid file  \"%s\"", fn);
		return;
	}
#ifdef _DEBUG
	codecid=CODEC_T44;
#else
	codecid=cli_select_codec();
#endif
	_8bit_codec=codec_is_8bit(codecid);
	
	printf("Testing  \"%s\"\n", fn);

	arg.mode=MODE_TEST;
	arg.loud=1;
	arg.codecid=codecid;
	arg.csize1=formatsize;
	open_fancy(fn, _8bit_codec, &arg, end, formatsize);
	
	process_sample(&arg);
}
static void encode_one_file(const char *srcfn, const char *dstfn)
{
	ThreadArgs arg={0};
	int codecid, _8bit_codec;
	int start=0, end=0;
	ptrdiff_t formatsize;

	get_filetitle(srcfn, -1, &start, &end);

	formatsize=get_filesize(srcfn);
	if(formatsize<1)
	{
		printf("Not a valid file  \"%s\"\n", srcfn);
		return;
	}
	
	codecid=cli_select_codec();
	_8bit_codec=codec_is_8bit(codecid);
	
	printf("Testing  \"%s\"\n", srcfn);

	arg.dstfn=dstfn;
	arg.mode=MODE_ENC;
	arg.loud=1;
	arg.codecid=codecid;
	arg.csize1=formatsize;
	open_fancy(srcfn, _8bit_codec, &arg, end, formatsize);
	
	process_sample(&arg);
}
static void decode_one_file(const char *srcfn, const char *dstfn)
{
	ThreadArgs arg={0};
	
	printf("Decoding  \"%s\"\n", srcfn);
	arg.cdata=load_file(srcfn, 1, 16, 0);
	if(!arg.cdata)
	{
		printf("Cannot open \'%s\'", srcfn);
		return;
	}
	arg.mode=MODE_DEC;
	arg.loud=1;
	process_sample(&arg);
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
		sum[0]+=(long long)dr*dr;
		sum[1]+=(long long)dg*dg;
		sum[2]+=(long long)db*db;
		sum[3]+=(long long)da*da;
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
	OP_TESTFILE, 1, 0,//operation, nthreads, formatsize

	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
//	"D:/ML/big_building.PPM"

//	"C:/Projects/datasets/kodim13.ppm"
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
		test_one_file(args.fn1, args.formatsize);
		break;
	case OP_TESTFOLDER:
		batch_test(args.fn1, args.nthreads);
		break;
	case OP_COMPRESS:
		encode_one_file(args.fn1, args.fn2);
		break;
	case OP_DECOMPRESS:
		decode_one_file(args.fn1, args.fn2);
		break;
	case OP_COMPARE:
		compare(args.fn1, args.fn2);
		break;
	}
	printf("Done.\n");
	pause();
	return 0;
}