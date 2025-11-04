#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

//	#define RELEASE
//	#define PROFILER

#ifdef PROFILER
#include"util.h"
#endif
#include<stdlib.h>
#if 1
int c01_codec(int argc, char **argv);//MT AC (mix4)
int c02_codec(int argc, char **argv);//MT binary
int c03_codec(int argc, char **argv);//MT binary (A2/WG4)
int c04_codec(int argc, char **argv);//MT Golomb-Rice WP3
int c05_codec(int argc, char **argv);//MT Golomb-Rice 3/7 simple preds
int c06_codec(int argc, char **argv);//MT binary
int c07_codec(int argc, char **argv);//MT o0 (binary) AC speed test
int c08_codec(int argc, char **argv);//ST disk AC5 test
int c09_codec(int argc, char **argv);//ST J2K CG o0 disk symbol ANS
int c10_codec(int argc, char **argv);//ST J2K CG o0 disk symbol AC
int c11_codec(int argc, char **argv);//ST SubG CG o0 disk nibble AC
int c12_codec(int argc, char **argv);//ST RCT o1 binary AC
int c13_codec(int argc, char **argv);//MT WG4_12
int c14_codec(int argc, char **argv);//MT WG4_8 binary AC
int c15_codec(int argc, char **argv);//MT WG4_8 binary AC
int c16_codec(int argc, char **argv);//ST SubG CG AC o0 separate loops
int c17_codec(int argc, char **argv);//ST SubPrev GR (speed test)
int c18_codec(int argc, char **argv);//MT SubG AC/GR
int c19_codec(int argc, char **argv);//ST CG3D static o0 AC		C10 is better
int c20_codec(int argc, char **argv);//ST/3T AVX2 3A3 static o1 rANS
int c21_codec(int argc, char **argv);//ST  deferred WP  snapshot-CDF  AVX2 AC/RC 12 prob bits
int c22_codec(int argc, char **argv);//ST/3T AVX2 484 static o1 rANS
int c23_codec(int argc, char **argv);//ST/3T AVX2 blockwise Shannon, static o1 rANS
int c24_codec(int argc, char **argv);//MT AC mix4
int c25_codec(int argc, char **argv);//MT AC mix8
int c26_codec(int argc, char **argv);//ST ABAC
int c27_codec(int argc, char **argv);//ST GR
int c28_codec(int argc, char **argv);//ST AC with GR context
int c29_codec(int argc, char **argv);//C29 ST interleaved AVX2 cRCT WG4/CG GRctx static-rANS
int c30_codec(int argc, char **argv);//ST OLS
int c31_codec(int argc, char **argv);//cross-platform C29, slower	X
int c32_codec(int argc, char **argv);//C32: like C29 but 16 coders
int c33_codec(int argc, char **argv);//C33: speed priority CG-only
int c34_codec(int argc, char **argv);//MT rANS
int c35_codec(int argc, char **argv);//binary FSM
int c36_codec(int argc, char **argv);//video test
int c37_codec(int argc, char **argv);//WP vs L1 test
int c38_codec(int argc, char **argv);//Synth-Natural codec
#endif


#ifndef CODEC_EXT
//	#define CODEC_EXT c01
//	#define CODEC_EXT c02
//	#define CODEC_EXT c03
//	#define CODEC_EXT c04
//	#define CODEC_EXT c05
//	#define CODEC_EXT c06
//	#define CODEC_EXT c07
//	#define CODEC_EXT c08
//	#define CODEC_EXT c09
//	#define CODEC_EXT c10
//	#define CODEC_EXT c11
//	#define CODEC_EXT c12
//	#define CODEC_EXT c13
//	#define CODEC_EXT c14
//	#define CODEC_EXT c15
//	#define CODEC_EXT c16
//	#define CODEC_EXT c17
//	#define CODEC_EXT c18
//	#define CODEC_EXT c19
//	#define CODEC_EXT c20
//	#define CODEC_EXT c21
//	#define CODEC_EXT c22
//	#define CODEC_EXT c23
//	#define CODEC_EXT c24
//	#define CODEC_EXT c25
//	#define CODEC_EXT c26
//	#define CODEC_EXT c27
//	#define CODEC_EXT c28
//	#define CODEC_EXT c29
//	#define CODEC_EXT c30
//	#define CODEC_EXT c31
//	#define CODEC_EXT c32
//	#define CODEC_EXT c33
//	#define CODEC_EXT c34
//	#define CODEC_EXT c35
//	#define CODEC_EXT c36
//	#define CODEC_EXT c37
	#define CODEC_EXT c38
#endif
#define STR_EXPAND(X) #X
#define STRINGIFY(X) STR_EXPAND(X)
#define GLUE_EXPAND(A, B) A##B
#define GLUE(A, B) GLUE_EXPAND(A, B)
#define CODEC_FUNC GLUE(CODEC_EXT, _codec)


int main(int argc, char **argv)
{
	int retcode=0;
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif
#if defined __GNUC__ || defined RELEASE
	retcode=CODEC_FUNC(argc, argv);
//#elif 1
//	const char *args[]=
//	{
//		argv[0],
//		"C:/dataset-DIV2K-ppm/0801.ppm"
//	};
//	if(CODEC_FUNC(_countof(args), (char**)args))
//		return 1;
//#elif 1
//	const char *args[]=
//	{
//		argv[0],
//
//		"D:/ML/zzz_deletethis.lsim",
//	//	"C:/dataset-panasonic-ppm/zzz.c34",
//		
//		"D:/ML/zzz_deletethis.ppm"
//	//	"C:/dataset-panasonic-ppm/zzz.ppm",
//	};
//	return c34_codec(_countof(args), (char**)args);
#else
	const char dstfn[]=//OVERWRITTEN
	//	"C:/Projects/datasets/zzz.ppm"

		"C:/dataset-a-temp/zzz.ppm"
	//	"D:/ML/zzz_deletethis.ppm"
	;
	const char tmpfn[]=//OVERWRITTEN
	//	"C:/Projects/datasets/zzz.lsim"
		
		"C:/dataset-a-temp/zzz.lsim"
	//	"D:/ML/zzz_deletethis.lsim"
	;
	const char srcfn[]=
	//	"C:/Projects/datasets/0868-ecrop.ppm"
	//	"C:/Projects/datasets/20240806 6 why me.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/Projects/datasets/dataset-CLIC30-ppm/03.ppm"
	//	"C:/Projects/datasets/dataset-DIV2K-ppm/0843.ppm"
	//	"C:/Projects/datasets/0801-cg.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim23.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13882.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"	//large
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim13-small16.PPM"
	//	"C:/Projects/datasets/kodim13-small4.PPM"
	//	"C:/Projects/datasets/temp.c18"
	//	"C:/Projects/datasets/20240414-noise.LSIM"
	//	"C:/Projects/datasets/20240414-noise.PPM"
	//	"C:/Projects/datasets/20240513 screenshot.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/kodim13-small16.ppm"
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim24.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
		

	//	"C:/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_cosmic-timetraveler-29758.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_rosan-harmens-18703.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-panasonic-ppm/P1000058.ppm"
	//	"C:/dataset-panasonic-ppm/P1000169.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13914.ppm"
	//	"C:/dataset-sintel-ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DSLR2x4-ppm/DSC_0133.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-HUGE-ppm/kodak.PPM"
	//	"C:/dataset-HUGE-ppm/space_huge.ppm"
	//	"C:/dataset-DSLR2-ppm/DSC_0320.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/dataset-sony-ppm/DSC00315.ppm"
	//	"D:/ML/zzz_halfbright.PPM"
	//	"C:/dataset-LPCB-ppm/STA13845.ppm"
	//	"C:/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13803.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13833.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13912.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13915.ppm"	//false color terrain
	//	"C:/dataset-LPCB-ppm/STA13843.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13844.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13845.ppm"	//space clouds
	//	"C:/dataset-synthetic-ppm/20240409 1 LPCB.ppm"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/dataset-CLIC303-ppm/2048x1320_lucas-lof-388.ppm"
	//	"D:/ML/dataset-kodak-ppm/kodim13.c01"
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_alberto-restifo-4549.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-DIV2K-ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DIV2K-ppm/0805.ppm"
	//	"C:/dataset-DIV2K-ppm/0807.ppm"
	//	"C:/dataset-DIV2K-ppm/0823.ppm"
	//	"C:/dataset-DIV2K-ppm/0843.ppm"
	//	"C:/dataset-DIV2K-ppm/0859.ppm"
	//	"C:/dataset-DIV2K-ppm/0864.ppm"
	//	"C:/dataset-DIV2K-ppm/0880.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-02.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-06.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-14.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-20.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-30.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-43.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-03.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-05.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-49.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-52.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-67.ppm"
	//	"C:/dataset-HUGE-ppm/blackmarble.ppm"
	//	"C:/dataset-HUGE-ppm/chaos1.ppm"
	//	"C:/dataset-HUGE-ppm/diagram.ppm"
	//	"C:/dataset-HUGE-ppm/gaia.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE2-ppm/andromeda.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13757.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_02.ppm"
	//	"C:/dataset-RAW-ppm/a0014-WP_CRW_6320.ppm"
	//	"C:/dataset-a70-ppm/20240816_113656_966.ppm"
	//	"C:/dataset-meme-ppm/emoji_u1f628.ppm"
	//	"C:/dataset-memes-ppm/usa.ppm"
	//	"C:/dataset-synth-ppm/20240421 1 the front.ppm"
	//	"C:/dataset-synth-ppm/20240516 4 DSC_0054.ppm"
		"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240407 blank.ppm"
	//	"C:/dataset-synth2-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-synth2-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-synth2-ppm/20240412 2 gralic enc.ppm"
	//	"C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm"
	//	"C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm"
	//	"C:/dataset-synth2-ppm/20240419 3.ppm"
	//	"C:/dataset-synth2-ppm/20240422 1.PPM"
	//	"C:/dataset-synth2-ppm/20240524 numbers.ppm"
	//	"C:/dataset-synth2-ppm/20241006 linux is cursed.ppm"
	//	"C:/dataset-synth2-ppm/art.ppm"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/checkboard.PPM"
	//	"D:/ML/dataset-CID22-ppm/pexels-photo-1933873.PPM"
	//	"D:/ML/kodim13.ppm"
	//	"D:/ML/kodim24.ppm"
	//	"D:/ML/nice_clock_face.ppm"
	//	"D:/ML/zzz_halfbright.PPM"
	//	"D:/Programs/c29/song.ppm"
	;
	const char *encargs[]=
	{
		argv[0],
		srcfn,
		tmpfn,
	//	"9",

	//	"-e", "0",

	//	"0",	//param1
	//	"7",	//near
	};
	const char *decargs[]=
	{
		argv[0],
		tmpfn,
		dstfn,
	};
	if(CODEC_FUNC(_countof(encargs), (char**)encargs))
		return 1;
	if(CODEC_FUNC(_countof(decargs), (char**)decargs))
		return 1;
#endif
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return retcode;
}
