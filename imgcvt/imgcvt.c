#include "imgcvt.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#include<Windows.h>
//#include<process.h>
//#define THREAD_CALL __stdcall
//typedef unsigned THREAD_RET;
//#else
//#include<x86intrin.h>
//#include<pthread.h>
//#define THREAD_CALL
//typedef void *THREAD_RET;
//#endif
static const char file[]=__FILE__;

	#define EXECUTE

#define BUFLEN 8192
#define MAXEXTLEN 128
static const char *g_extensions[]=
{
	"ppm", "pgm",
	"png",
	"jxl",
	"webp",
	"avif",
	"gra111",
	"flic",
	"qlic2",
	"qlic",
	"qic",
	"lea",
	"bmp",
	"tif", "tiff",
	"jpg", "jpeg",
};
typedef enum FormatEnum
{
	FMT_NONE,
	FMT_PNM,
	FMT_PGM,
	FMT_PPM,
	FMT_PNG,
	FMT_JXL,
	FMT_WEBP,
	FMT_AVIF,
	FMT_GRA111,
	FMT_FLIC,
	FMT_QLIC2,
	FMT_QLIC,
	FMT_QIC,
	FMT_LEA,
	FMT_BMP,
	FMT_TIF,
	FMT_JPEG,
	FMT_COUNT,
} Format;
static const char* get_codecname(Format dstfmt)
{
	switch(dstfmt)
	{
	case FMT_PNM:		return "FFmpeg";
	case FMT_PGM:		return "FFmpeg";
	case FMT_PPM:		return "FFmpeg";
	case FMT_PNG:		return "FFmpeg";
	case FMT_JXL:		return "libjxl";
	case FMT_WEBP:		return "libwebp";
	case FMT_AVIF:		return "libavif";
	case FMT_GRA111:	return "Gralic by Alexander Rhatushnyak";
	case FMT_FLIC:		return "Flic by Alexander Rhatushnyak";
	case FMT_QLIC2:		return "Qlic2 by Alexander Rhatushnyak";
	case FMT_QLIC:		return "Qlic by Alexander Rhatushnyak";
	case FMT_QIC:		return "Qic by Alexander Rhatushnyak";
	case FMT_LEA:		return "LEA by Marcio Pais";
	case FMT_BMP:		return "FFmpeg";
	case FMT_TIF:		return "FFmpeg";
	case FMT_JPEG:		return "FFmpeg";
	default:
		break;
	}
	return "UNKNOWN CODEC";
}
static Format classify_extension(const char *ext)
{
	if(!_stricmp(ext, ".PNM"))	return FMT_PNM;
	if(!_stricmp(ext, ".PPM"))	return FMT_PPM;
	if(!_stricmp(ext, ".PGM"))	return FMT_PGM;
	if(!_stricmp(ext, ".PNG"))	return FMT_PNG;
	if(!_stricmp(ext, ".JXL"))	return FMT_JXL;
	if(!_stricmp(ext, ".WEBP"))	return FMT_WEBP;
	if(!_stricmp(ext, ".AVIF"))	return FMT_AVIF;
	if(!_stricmp(ext, ".GRA111"))	return FMT_GRA111;
	if(!_stricmp(ext, ".FLIC"))	return FMT_FLIC;
	if(!_stricmp(ext, ".QLIC2"))	return FMT_QLIC2;
	if(!_stricmp(ext, ".QLIC"))	return FMT_QLIC;
	if(!_stricmp(ext, ".QIC"))	return FMT_QIC;
	if(!_stricmp(ext, ".LEA"))	return FMT_LEA;
	if(!_stricmp(ext, ".BMP"))	return FMT_BMP;
	if(!_stricmp(ext, ".TIF"))	return FMT_TIF;
	if(!_stricmp(ext, ".TIFF"))	return FMT_TIF;
	if(!_stricmp(ext, ".JPG"))	return FMT_JPEG;
	if(!_stricmp(ext, ".JPEG"))	return FMT_JPEG;
	LOG_ERROR("Unknown extension \"%s\"", ext);
	return FMT_NONE;
}
static int is_codec_pnm(Format codec)
{
	switch(codec)
	{
	case FMT_PNM:		return 0;
	case FMT_PGM:		return 0;
	case FMT_PPM:		return 0;
	case FMT_PNG:		return 0;
	case FMT_JXL:		return 0;
	case FMT_WEBP:		return 0;
	case FMT_AVIF:		return 0;
	case FMT_GRA111:	return 1;
	case FMT_FLIC:		return 1;
	case FMT_QLIC2:		return 1;
	case FMT_QLIC:		return 1;
	case FMT_QIC:		return 1;
	case FMT_LEA:		return 1;
	case FMT_BMP:		return 0;
	case FMT_TIF:		return 0;
	case FMT_JPEG:		return 0;
	default:
		break;
	}
	LOG_ERROR("Unknown format");
	return 0;
}
static const char* cmd_enc_pnm(Format srcfmt)
{
	switch(srcfmt)
	{
	case FMT_PNM:		return 0;
	case FMT_PGM:		return 0;
	case FMT_PPM:		return 0;
	case FMT_PNG:		return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_JXL:		return "djxl";
	case FMT_WEBP:		return "dwebp";
	case FMT_AVIF:		return "avifdec";
	case FMT_GRA111:	return "gralic111d d";
	case FMT_FLIC:		return "flic d";
	case FMT_QLIC2:		return "qlic2 d";
	case FMT_QLIC:		return "qlic d";
	case FMT_QIC:		return "qic d";
	case FMT_LEA:		return "dlea";
	case FMT_BMP:		return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_TIF:		return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_JPEG:		return "ffmpeg -hide_banner -loglevel error -i";
	default:
		break;
	}
	LOG_ERROR("Unknown format");
	return 0;
}
static const char* cmd_enc_direct(Format dstfmt, int *swap_args)
{
	switch(dstfmt)
	{
	case FMT_PNM:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_PGM:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_PPM:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_PNG:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_JXL:		*swap_args=0; return "cjxl -d 0 -e 9";
	case FMT_WEBP:		*swap_args=0; return "cwebp -lossless";
	case FMT_AVIF:		*swap_args=0; return "avifenc -l -s 0";
	case FMT_GRA111:	*swap_args=1; return "gralic111d c";
	case FMT_FLIC:		*swap_args=1; return "flic c";
	case FMT_QLIC2:		*swap_args=1; return "qlic2 c";
	case FMT_QLIC:		*swap_args=1; return "qlic c";
	case FMT_QIC:		*swap_args=1; return "qic c";
	case FMT_LEA:		*swap_args=0; return "clea";
	case FMT_BMP:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_TIF:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	case FMT_JPEG:		*swap_args=0; return "ffmpeg -hide_banner -loglevel error -i";
	default:
		break;
	}
	LOG_ERROR("Unknown format");
	return 0;
}
static double va_system(char *tmp, size_t tmplen, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	vsnprintf(tmp, tmplen-1, fmt, args);
	va_end(args);
	
	double t=time_sec();
#ifdef EXECUTE
	system(tmp);
#else
	printf("%s\n", tmp);
#endif
	return time_sec()-t;
}
typedef struct CvtResultStruct
{
	size_t srcsize, dstsize;
	double elapsed[2];
} CvtResult;
static int convert1(const char *srcfn, const char *dstfn, CvtResult *result)
{
	ptrdiff_t srcsize=get_filesize(srcfn);
	ptrdiff_t dstsize=get_filesize(dstfn);
	if(srcsize<1)//source is not accessible or not a file
	{
		if(srcsize)
			printf("Src doesn't exist:\t\"%s\"\n", srcfn);
		else
			printf("Src is a folder:\t\"%s\"\n", srcfn);
		return 1;
	}
	if(dstsize>=0)//destination already exists
	{
		printf("Dst already exists:\t\"%s\"\n", dstfn);
		return 1;
	}
	int srcstart=0, srcend=0;
	int dststart=0, dstend=0;
	get_filetitle(srcfn, -1, &srcstart, &srcend);
	get_filetitle(dstfn, -1, &dststart, &dstend);
	Format srcfmt=classify_extension(srcfn+srcend);
	Format dstfmt=classify_extension(dstfn+dstend);
	int srcpnm=is_codec_pnm(srcfmt);
	int dstpnm=is_codec_pnm(dstfmt);

	char *buf=(char*)malloc(BUFLEN);
	if(!buf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	const char *cmd=0;
	char *tmpfn=0;
	double elapsed[2]={0};
	if(
		(srcpnm&&(unsigned)(dstfmt-FMT_PNM)>(unsigned)(FMT_PPM-FMT_PNM))||
		(dstpnm&&(unsigned)(srcfmt-FMT_PNM)>(unsigned)(FMT_PPM-FMT_PNM))
	)//create temp PPM when you need it but didn't have it
	{
		cmd=cmd_enc_pnm(srcfmt);
		if(!cmd)
			return 1;
		tmpfn=(char*)malloc(BUFLEN);
		snprintf(tmpfn, BUFLEN-1LL, "%.*s.PPM", dstend, dstfn);

		elapsed[1]=va_system(buf, BUFLEN, "\"%s \"%s\" \"%s\"\"", cmd, srcfn, tmpfn);
	}
	int swapargs=0;
	cmd=cmd_enc_direct(dstfmt, &swapargs);
	if(!cmd)
		return 1;
	if(swapargs)
		elapsed[0]=va_system(buf, BUFLEN, "\"%s \"%s\" \"%s\"\"", cmd, dstfn, tmpfn?tmpfn:srcfn);
	else
		elapsed[0]=va_system(buf, BUFLEN, "\"%s \"%s\" \"%s\"\"", cmd, tmpfn?tmpfn:srcfn, dstfn);
	dstsize=get_filesize(dstfn);
	if(result)
	{
		result->srcsize+=srcsize;
		result->dstsize+=dstsize;
		result->elapsed[0]+=elapsed[0];
		result->elapsed[1]+=elapsed[1];
	}
	if(tmpfn)
		free(tmpfn);
	free(buf);
	return 0;
}
static void print_result(CvtResult const *result)
{
	printf("%12zd -> %12zd bytes", result->srcsize, result->dstsize);
	if(result->elapsed[1])
		printf("  A %12lf  B %12lf sec\n", result->elapsed[1], result->elapsed[0]);
	else
		printf("  E %12lf sec\n", result->elapsed[0]);
}

static void print_usage(const char *argv0)
{
	int start=0, end=0;
	get_filetitle(argv0, -1, &start, &end);
	printf(
		"Usage:\n"
		"  %.*s  srcpath  dstpath  format\n"
		"    Converts all applicable images to the selected format.\n"
		"  %.*s  srcimg  dstimg\n"
		"    Converts one image to format from the extension.\n"
		"Input formats:  PPM/PGM/PNG/TIFF/BMP/JPEG\n"
		"Output formats:\n"
		"  \"pnm\"       FFmpeg (PPM/PGM)\n"
		"  \"png\"       FFmpeg\n"
		"  \"jxl\"       cjxl -d 0 -e 9 / djxl\n"
		"  \"webp\"      cwebp -lossless\n"
		"  \"avif\"      avifenc -l -s 0 / avifdec\n"
		"  \"gra111\"    Gralic111d  by Alexander Rhatushnyak\n"
		"  \"flic\"      Flic        by Alexander Rhatushnyak\n"
		"  \"qlic2\"     Qlic2       by Alexander Rhatushnyak\n"
		"  \"qlic\"      Qlic        by Alexander Rhatushnyak\n"
		"  \"qic\"       Qic         by Alexander Rhatushnyak\n"
		"  \"lea\"       clea/dlea   by Marcio Pais\n"
		"\n",
		end-start, argv0+start,
		end-start, argv0+start
	);
}
int main(int argc, char **argv)
{
	double t=time_sec();
	printf("ImageConverter  %s, %s\n\n", __DATE__, __TIME__);
#ifdef _DEBUG
	const char
		*srcarg="D:/ML/dataset-kodak-ppm",
		*dstarg="D:/ML/temp",
		*fmtarg="flic";

	//	*srcarg="D:/ML/dataset-sintel-ppm",
	//	*dstarg="D:/ML/dataset-sintel-flic",
	//	*fmtarg="flic";

	//	*srcarg="D:/ML/dataset-sintel",
	//	*dstarg="D:/ML/dataset-sintel-gralic",
	//	*fmtarg="gra111";
	
	//	*srcarg="C:/Projects/datasets/dataset-kodak",
	//	*dstarg="C:/Projects/datasets/dataset-temp",
	//	*fmtarg="lea";

	//	*srcarg="C:/Projects/datasets/dataset-nasa",
	//	*dstarg="C:/Projects/datasets/dataset-temp",
	//	*fmtarg="gra111";

	//	*srcarg="C:/Projects/datasets/dataset-kodak-grey-gralic",
	//	*dstarg="C:/Projects/datasets/dataset-temp",
	//	*fmtarg="PNG";
#else
	if(argc!=3&&argc!=4)
	{
		print_usage(argv[0]);
		pause();
		return 1;
	}
	const char *srcarg=argv[1], *dstarg=argv[2], *fmtarg=argc==4?argv[3]:0;
#endif
	size_t size0=get_filesize(srcarg);
	if(size0<0)//inaccessible
	{
		printf("Cannot open \"%s\"\n", srcarg);
		pause();
		return 1;
	}
	CvtResult total={0};
	if(size0)//file
	{
#ifndef _DEBUG
		if(argc!=3)
		{
			print_usage(argv[0]);
			pause();
			return 1;
		}
#endif
		convert1(srcarg, dstarg, &total);
		print_result(&total);
	}
	else//folder
	{
#ifndef _DEBUG
		if(argc!=4)
		{
			print_usage(argv[0]);
			pause();
			return 1;
		}
#endif
		char ext[MAXEXTLEN]={'.'};
		int extlen=(int)strlen(fmtarg);
		if(extlen>_countof(ext)-2)
		{
			printf("Invalid codec  \"%s\"\n", fmtarg);
			return 1;
		}
		strcpy(ext+1, fmtarg);
		Format dstfmt=classify_extension(ext);
		const char *codecname=get_codecname(dstfmt);
		if(codecname)
			printf("%s:\n", codecname);

		ArrayHandle filenames=get_filenames(srcarg, g_extensions, _countof(g_extensions), 1);
		if(!filenames)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		if(!filenames->count)
		{
			printf("No applicable images found in \"%s\"\n", srcarg);
			return 1;
		}
		ArrayHandle dstpath=filter_path(dstarg, -1);
		char *dstfn=(char*)malloc(BUFLEN);
		if(!dstfn)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		int ignoreerrors=0, errorpass=rand()%10;
		for(int k=0;k<(int)filenames->count;++k)
		{
			ArrayHandle *srcfn=(ArrayHandle*)array_at(&filenames, k);

			int tidx=0, pidx=0;
			get_filetitle((char*)srcfn[0]->data, (int)srcfn[0]->count, &tidx, &pidx);

			snprintf(dstfn, BUFLEN-1, "%s%.*s.%s", dstpath->data, pidx-tidx, srcfn[0]->data+tidx, fmtarg);

			CvtResult result={0};
			int error=convert1((char*)srcfn[0]->data, dstfn, &result);
			print_result(&result);

			if(error)
			{
				printf(
					"Conversion ERROR.\n"
					"Enter %d to continue, ignoring further errors: ",
					errorpass
				);
				while(!scanf(" %d", &ignoreerrors));
				if(ignoreerrors!=errorpass)
					break;
			}
			else
			{
				total.srcsize+=result.srcsize;
				total.dstsize+=result.dstsize;
				total.elapsed[0]+=result.elapsed[0];
				total.elapsed[1]+=result.elapsed[1];
			}
		}
		printf("Total:\n");
		print_result(&total);
		free(dstfn);
	}

	printf("\nDone.  %12lf sec\n", time_sec()-t);
	pause();
	return 0;
}
