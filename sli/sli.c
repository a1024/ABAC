
int a01_codec(int argc, char **argv);
int a02_codec(int argc, char **argv);

#ifndef CODEC_EXT
//	#define CODEC_EXT a01
	#define CODEC_EXT a02
//	#define CODEC_EXT a03
//	#define CODEC_EXT a04
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
#ifdef __GNUC__
	return CODEC_FUNC(argc, argv);
#else
	const char *srcfn=
		
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
		"C:/dataset-HUGE-ppm/juno.ppm"
	//	"D:/ML/tobruk.ppm"

	;
	const char *encargs[]=
	{
		argv[0],
		srcfn,
		"C:/dataset-HUGE-ppm/mystery.alic",
	//	"D:/ML/mystery.alic",
	};
	const char *decargs[]=
	{
		argv[0],
		"C:/dataset-HUGE-ppm/mystery.alic",
		"C:/dataset-HUGE-ppm/mystery.ppm",
	//	"D:/ML/mystery.alic",
	//	"D:/ML/mystery.ppm",
	};
	if(CODEC_FUNC(_countof(encargs), (char**)(long long)encargs))
		return 1;
	//if(CODEC_FUNC(_countof(decargs), (char**)(long long)decargs))
	//	return 1;
	return 0;
#endif
}
