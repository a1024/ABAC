#include<stdio.h>
#include<stdlib.h>

int s01_codec(const char *command, const char *srcfn, const char *dstfn);//binary
int s02_codec(const char *command, const char *srcfn, const char *dstfn);//ANS
int s03_codec(const char *command, const char *srcfn, const char *dstfn);//GR
int s04_codec(const char *command, const char *srcfn, const char *dstfn);//binary
int s05_codec(const char *command, const char *srcfn, const char *dstfn);//FSE  (WIP)
int s06_codec(const char *command, const char *srcfn, const char *dstfn);//RLE+GR
int s07_codec(const char *command, const char *srcfn, const char *dstfn);//AC

#ifndef CODECFUNC
//	#define CODECFUNC s01_codec
//	#define CODECFUNC s02_codec
//	#define CODECFUNC s03_codec
//	#define CODECFUNC s04_codec
//	#define CODECFUNC s05_codec
//	#define CODECFUNC s06_codec
	#define CODECFUNC s07_codec
#endif
int main(int argc, char **argv)
{
#ifdef _DEBUG
	{
		const char *srcfn=0, *tmpfn=0, *tmpf2=0;

		srcfn=

		//	"C:/Projects/datasets/dataset-synth-ppm/20231214 1 screenshot.ppm"
		//	"C:/Projects/datasets/dataset-synth-ppm/20231214 2 screenshot LZ.ppm"
		//	"C:/Projects/datasets/20240513 screenshot.PPM"

			"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
		//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
		//	"F:/Projects/dataset-GDCC2020-ppm/astro-01.ppm"

		;
		tmpfn=
			"C:/Projects/datasets/zzz.sli"
		//	"F:/Projects/zzz_tmp.s01"
		;
		tmpf2=
			"C:/Projects/datasets/zzz.ppm"
		//	"F:/Projects/zzz_tmp.ppm"
		;
		if(CODECFUNC("e", srcfn, tmpfn))
			return 1;
		if(CODECFUNC("d", tmpfn, tmpf2))
			return 1;
#ifdef _MSC_VER
		system("pause");
#endif
		return 0;
	}
#else
	if(argc!=4)
	{
		printf(
			"Usage: \"%s\"  e|d  input  output\n"
			"Only 24-bit PPM is supported\n"
			, argv[0]
		);
		return 1;
	}
	CODECFUNC(argv[1], argv[2], argv[3]);
#endif
	return 0;
}