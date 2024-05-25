#include "imgcvt.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#endif
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

	#define CMD_EXEC

#define CMD_MAX 4
#define BUFLEN 65536
#define MAXEXTLEN 128
static int effort=-1;
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
	"emma",
	"lea",
	"halic",
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
	FMT_EMMA,
	FMT_LEA,
	FMT_HALIC,
	FMT_BMP,
	FMT_TIF,
	FMT_JPEG,
	FMT_COUNT,
} Format;
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
	if(!_stricmp(ext, ".GRALIC"))	return FMT_GRA111;
	if(!_stricmp(ext, ".FLIC"))	return FMT_FLIC;
	if(!_stricmp(ext, ".QLIC2"))	return FMT_QLIC2;
	if(!_stricmp(ext, ".QLIC"))	return FMT_QLIC;
	if(!_stricmp(ext, ".QIC"))	return FMT_QIC;
	if(!_stricmp(ext, ".EMMA"))	return FMT_EMMA;
	if(!_stricmp(ext, ".LEA"))	return FMT_LEA;
	if(!_stricmp(ext, ".HALIC"))	return FMT_HALIC;
	if(!_stricmp(ext, ".BMP"))	return FMT_BMP;
	if(!_stricmp(ext, ".TIF"))	return FMT_TIF;
	if(!_stricmp(ext, ".TIFF"))	return FMT_TIF;
	if(!_stricmp(ext, ".JPG"))	return FMT_JPEG;
	if(!_stricmp(ext, ".JPEG"))	return FMT_JPEG;
	LOG_ERROR("Unknown extension \"%s\"", ext);
	return FMT_NONE;
}
#if 0
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
	case FMT_EMMA:		return "EMMA by Marcio Pais";
	case FMT_LEA:		return "LEA by Marcio Pais";
	case FMT_HALIC:		return "HALIC by Hakan Abbas";
	case FMT_BMP:		return "FFmpeg";
	case FMT_TIF:		return "FFmpeg";
	case FMT_JPEG:		return "FFmpeg";
	default:
		break;
	}
	return "UNKNOWN CODEC";
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
#endif
typedef struct CvtResultStruct
{
	size_t srcsize, dstsize[CMD_MAX];
	double elapsed[CMD_MAX];
} CvtResult;
static int convert1(int loud, const char *srcfn, const char *dstfn, CvtResult *result)
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

	char *buf=(char*)malloc(BUFLEN);
	if(!buf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	char *tmpfn=0;
#if 1
	if(dstfmt==FMT_JXL&&effort<0)
	{
		printf("Enter effort [0~10]: ");
		while(!scanf(" %d", &effort)&&(unsigned)effort>10);
	}
	const char *cmd[CMD_MAX]={0}, *dst[CMD_MAX]={0};
	int ncmd=0, printed=0;
#ifdef _WIN32
	int justcopy=0;
#define CMD_APPEND(STR, ...) cmd[ncmd]=buf+printed, printed+=snprintf(buf+printed, BUFLEN-printed-1, "\"" STR "\"", __VA_ARGS__)+1
#else
#define CMD_APPEND(STR, ...) cmd[ncmd]=buf+printed, printed+=snprintf(buf+printed, BUFLEN-printed-1, STR, __VA_ARGS__)+1
#endif
#define DST_APPEND(STR, ...) dst[ncmd]=buf+printed, printed+=snprintf(buf+printed, BUFLEN-printed-1, STR, __VA_ARGS__)+1
	if(srcfmt==dstfmt)
	{
#ifdef _WIN32
		justcopy=1;
		//CMD_APPEND("copy /b \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
#else
		CMD_APPEND("cp \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
#endif
	}
	else
	{
		switch(srcfmt)
		{
		case FMT_PNM:case FMT_PGM:case FMT_PPM:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
			case FMT_PNG:case FMT_BMP:case FMT_TIF:case FMT_JPEG:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("cjxl -d 0 -e %d \"%s\" \"%s\"", effort, srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("cwebp -lossless \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("gralic111d c \"%s\" \"%s\"", dstfn, srcfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("flic c \"%s\" \"%s\"", dstfn, srcfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("qlic2 c \"%s\" \"%s\"", dstfn, srcfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("qlic c \"%s\" \"%s\"", dstfn, srcfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("qic c \"%s\" \"%s\"", dstfn, srcfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("emma_c \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("clea \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%s\" \"%s\" -mt", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_PNG:case FMT_BMP:case FMT_TIF:case FMT_JPEG:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
			case FMT_PNG:case FMT_BMP:case FMT_TIF:case FMT_JPEG:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("cjxl -d 0 -e %d \"%s\" \"%s\"", effort, srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("cwebp -lossless \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				if(srcfmt==FMT_PNG)
					CMD_APPEND("avifenc -l -s 4 \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				else
				{
					CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
					CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				}
				break;
			case FMT_GRA111:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_JXL:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
			case FMT_PNG:case FMT_BMP:case FMT_TIF:case FMT_JPEG:
				CMD_APPEND("djxl \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("cjxl -d 0 -e %d \"%s\" \"%s\"", effort, srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("djxl \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_WEBP:
			switch(dstfmt)
			{
			case FMT_PNM:
			case FMT_PPM:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PGM:
				CMD_APPEND("dwebp \"%s\" -pgm -o \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
				CMD_APPEND("dwebp \"%s\" -o \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_BMP:
				CMD_APPEND("dwebp \"%s\" -bmp -o \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_TIF:
				CMD_APPEND("dwebp \"%s\" -tiff -o \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JPEG:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				break;
			case FMT_AVIF:
				CMD_APPEND("dwebp \"%s\" -o \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("djxl \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("dwebp \"%s\" -ppm -o \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_AVIF:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
				CMD_APPEND("avifdec \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PNG\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				break;
			case FMT_GRA111:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("avifdec \"%s\" \"%.*s.PNG\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PNG\" \"%.*s.PPM\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_GRA111:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("gralic111d d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				break;
			case FMT_FLIC:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("gralic111d d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_FLIC:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("flic d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				break;
			case FMT_QLIC2:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("flic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_QLIC2:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("qlic2 d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				break;
			case FMT_QLIC:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("qlic2 d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_QLIC:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("qlic d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				break;
			case FMT_QIC:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("qlic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_QIC:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("qic d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				break;
			case FMT_EMMA:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("qic d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_EMMA:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("emma_d \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				break;
			case FMT_LEA:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				CMD_APPEND("emma_d \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_LEA:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("dlea \"%s\" \"%s\"", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				break;
			case FMT_HALIC:
				CMD_APPEND("dlea \"%s\" \"%.*s.PPM\"", srcfn, dstend, dstfn), DST_APPEND("%.*s.PPM", dstend, dstfn), ++ncmd;
				CMD_APPEND("HALIC_ENCODE_V.0.7.1_MT \"%.*s.PPM\" \"%s\" -mt", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			default:
				break;
			}
			break;
		case FMT_HALIC:
			switch(dstfmt)
			{
			case FMT_PNM:case FMT_PGM:case FMT_PPM:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%s\" -mt", srcfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_PNG:
			case FMT_BMP:
			case FMT_TIF:
			case FMT_JPEG:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_JXL:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("cjxl -d 0 -e %d \"%.*s.PPM\" \"%s\"", effort, dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_WEBP:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("cwebp -lossless \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_AVIF:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("ffmpeg -hide_banner -loglevel error -i \"%.*s.PPM\" \"%.*s.PNG\"", dstend, dstfn, dstend, dstfn), DST_APPEND("%.*s.PNG", dstend, dstfn), ++ncmd;
				CMD_APPEND("avifenc -l -s 4 \"%.*s.PNG\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_GRA111:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("gralic111d c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_FLIC:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("flic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC2:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("qlic2 c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QLIC:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("qlic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_QIC:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("qic c \"%s\" \"%.*s.PPM\"", dstfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_EMMA:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("emma_c \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_LEA:
				CMD_APPEND("HALIC_DECODE_V.0.7.1_MT \"%s\" \"%.*s.PPM\" -mt", srcfn, dstend, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				CMD_APPEND("clea \"%.*s.PPM\" \"%s\"", dstend, dstfn, dstfn), DST_APPEND("%s", dstfn), ++ncmd;
				break;
			case FMT_HALIC:
				break;
			default:
				break;
			}
			break;
		default:
			break;
		}
	}

	result->srcsize+=srcsize;
#ifdef _WIN32
	if(justcopy)
	{
		int success;
		double t;
		
		t=time_sec();
		success=CopyFileA(srcfn, dstfn, 1);
		t=time_sec()-t;

		if(success)
		{
			ptrdiff_t size2=get_filesize(dstfn);
			result->dstsize[0]+=size2;
			result->elapsed[0]+=t;
		}
		else
		{
			unsigned e=GetLastError();
			printf("Error %d copying  \"%s\"\n", e, srcfn);
		}
	}
	else
#endif
	for(int k=0;k<ncmd;++k)
	{
#ifdef CMD_EXEC
		double t;
		ptrdiff_t size2;
#endif
		if(loud)
			printf("%s\n", cmd[k]);
#ifdef CMD_EXEC
		t=time_sec();
		system(cmd[k]);
		t=time_sec()-t;

		size2=get_filesize(dst[k]);
		if(size2>0)
		{
			result->dstsize[k]+=size2;
			result->elapsed[k]+=t;
		}
		else
		{
			printf("Failed to convert  \"%s\"\n", srcfn);
			system("echo %errorlevel%");
		}
#endif
	}
#endif
#if 0
	double elapsed[2]={0};
	const char *cmd=0;
	int srcpnm=is_codec_pnm(srcfmt);
	int dstpnm=is_codec_pnm(dstfmt);
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
#endif
	//dstsize=get_filesize(dstfn);
	//if(result)
	//{
	//	result->srcsize+=srcsize;
	//	result->dstsize+=dstsize;
	//	result->elapsed[0]+=elapsed[0];
	//	result->elapsed[1]+=elapsed[1];
	//}
	if(tmpfn)
		free(tmpfn);
	free(buf);
	return 0;
}
static void print_result(CvtResult const *result)
{
	printf("%12zd", result->srcsize);
	for(int k=0;k<CMD_MAX;++k)
	{
		if(!result->elapsed[k])
			break;
		printf(" -> %12zd", result->dstsize[k]);
	}
	printf(" bytes  E");
	for(int k=0;k<CMD_MAX;++k)
	{
		if(!result->elapsed[k])
			break;
		printf(" %12lf", result->elapsed[k]);
	}
	printf(" sec\n");
	//printf("%12zd -> %12zd bytes", result->srcsize, result->dstsize);
	//if(result->elapsed[1])
	//	printf("  A %12lf  B %12lf sec\n", result->elapsed[1], result->elapsed[0]);
	//else
	//	printf("  E %12lf sec\n", result->elapsed[0]);
}

static void print_usage(const char *argv0)
{
	int start=0, end=0;
	get_filetitle(argv0, -1, &start, &end);
	printf(
		"Usage:\n"
		"  %.*s  srcpath  dstpath  format  [q]\n"
		"    Converts all applicable images to the selected format.\n"
		"  %.*s  srcimg  dstimg  [q]\n"
		"    Converts one image to format from the extension.\n"
		"  \"q\"         Silences commands.\n"
		"Formats:\n"
		"  \"ppm\"       FFmpeg\n"
		"  \"pgm\"       FFmpeg\n"
		"  \"bmp\"       FFmpeg\n"
		"  \"tiff\"      FFmpeg\n"
		"  \"jpeg\"      FFmpeg\n"
		"  \"jxl\"       libjxl -d 0 -e 9\n"
		"  \"webp\"      libwebp -lossless\n"
		"  \"avif\"      libavif -l -s 0\n"
		"  \"gra111\"    Gralic111d      by Alexander Rhatushnyak\n"
		"  \"flic\"      Flic            by Alexander Rhatushnyak\n"
		"  \"qlic2\"     Qlic2           by Alexander Rhatushnyak\n"
		"  \"qlic\"      Qlic            by Alexander Rhatushnyak\n"
		"  \"qic\"       Qic             by Alexander Rhatushnyak\n"
		"  \"emma\"      emma_c/emma_d   by Marcio Pais\n"
		"  \"lea\"       clea/dlea       by Marcio Pais\n"
		"  \"halic\"     HALIC_ENCODE_V.0.7.1_MT / HALIC_DECODE_V.0.7.1_MT    by Hakan Abbas\n"
		"These codecs must be present in PATH.\n"
		"\n",
		end-start, argv0+start,
		end-start, argv0+start
	);
}
int main(int argc, char **argv)
{
	double t=time_sec();
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H:%M:%S");
	printf("ImageConverter  %s, %s  started on %s\n\n", __DATE__, __TIME__, g_buf);
	int loud=1;
#ifdef _DEBUG
	const char
		*srcarg="D:/ML/dataset-CLIC30-ppm",
		*dstarg="D:/ML/dataset-CLIC30-temp",
		*fmtarg="ppm";

	//	*srcarg="D:/ML/dataset-LPCB-halic",
	//	*dstarg="D:/ML/dataset-LPCB-jxl",
	//	*fmtarg="ppm";

	//	*srcarg="D:/ML/dataset-LPCB-lea",
	//	*dstarg="D:/ML/dataset-LPCB-ppm",
	//	*fmtarg="ppm";

	//	*srcarg="D:/ML/dataset-kodak-ppm",
	//	*dstarg="D:/ML/temp",
	//	*fmtarg="flic";

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
	if(argc!=3&&argc!=4&&argc!=5)
	{
		print_usage(argv[0]);
		pause();
		return 1;
	}
	const char *srcarg=argv[1], *dstarg=argv[2], *fmtarg=argc>=4?argv[3]:0;
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
		if(argc!=3&&argc!=4)
		{
			print_usage(argv[0]);
			pause();
			return 1;
		}
		if(argc==4&&strlen(argv[3])==1&&(argv[3][0]&0xDF)=='Q')
			loud=0;
#endif
		convert1(loud, srcarg, dstarg, &total);
		print_result(&total);
	}
	else//folder
	{
#ifndef _DEBUG
		if(argc!=4&&argc!=5)
		{
			print_usage(argv[0]);
			pause();
			return 1;
		}
		if(argc==5&&strlen(argv[4])==1&&(argv[4][0]&0xDF)=='Q')
			loud=0;
#endif
		char ext[MAXEXTLEN]={'.'};
		int extlen=(int)strlen(fmtarg);
		if(extlen>_countof(ext)-2)
		{
			printf("Invalid codec  \"%s\"\n", fmtarg);
			return 1;
		}
		strcpy(ext+1, fmtarg);
		//Format dstfmt=classify_extension(ext);
		//const char *codecname=get_codecname(dstfmt);
		//if(codecname)
		//	printf("%s:\n", codecname);

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
			int error=convert1(loud, (char*)srcfn[0]->data, dstfn, &result);
			printf("%5d/%5d  ", k+1, (int)filenames->count);
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
				for(int k=0;k<CMD_MAX;++k)
				{
					total.dstsize[k]+=result.dstsize[k];
					total.elapsed[k]+=result.elapsed[k];
				}
				//total.dstsize+=result.dstsize;
				//total.elapsed[0]+=result.elapsed[0];
				//total.elapsed[1]+=result.elapsed[1];
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
