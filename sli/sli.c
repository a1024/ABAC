﻿
int a01_codec(int argc, char **argv);//best
int a02_codec(int argc, char **argv);//R/W experiment
int a03_codec(int argc, char **argv);//experimental
int a04_codec(int argc, char **argv);//fast
int a05_codec(int argc, char **argv);//extremely fast, compromise


#ifndef CODEC_EXT
//	#define CODEC_EXT a01
//	#define CODEC_EXT a02
//	#define CODEC_EXT a03
	#define CODEC_EXT a04
//	#define CODEC_EXT a05
#endif
#define STR_EXPAND(X) #X
#define STRINGIFY(X) STR_EXPAND(X)
#define GLUE_EXPAND(A, B) A##B
#define GLUE(A, B) GLUE_EXPAND(A, B)
#define CODEC_FUNC GLUE(CODEC_EXT, _codec)
#ifndef _countof
#define _countof(A) (sizeof(A)/sizeof(*(A)))
#endif

int main(int argc, char **argv)
{
#if defined __GNUC__ || (!defined _DEBUG && !defined RELWITHDEBINFO)
	return CODEC_FUNC(argc, argv);
#else
	const char *srcfn=

		"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"

	//	"C:/Projects/datasets/dataset-DSLR-ppm/DSC_0825.ppm"
	//	"C:/dataset-NEF-ppm/DSC_0185.ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/juno.ppm"
	//	"D:/ML/tobruk.ppm"

	;
	const char *encargs[]=
	{
		argv[0],
		srcfn,
		"C:/Projects/datasets/mystery.alic",
	//	"C:/dataset-HUGE-ppm/mystery.alic",
	//	"D:/ML/mystery.alic",
	};
	const char *decargs[]=
	{
		argv[0],
		"C:/Projects/datasets/mystery.alic",
		"C:/Projects/datasets/mystery.ppm",
	//	"C:/dataset-HUGE-ppm/mystery.alic",
	//	"C:/dataset-HUGE-ppm/mystery.ppm",
	//	"D:/ML/mystery.alic",
	//	"D:/ML/mystery.ppm",
	};
	if(CODEC_FUNC(_countof(encargs), (char**)(long long)encargs))
		return 1;
	if(CODEC_FUNC(_countof(decargs), (char**)(long long)decargs))
		return 1;
	return 0;
#endif
}
