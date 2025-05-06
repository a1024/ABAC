#include"codec.h"
//	#define PROFILER
#ifdef PROFILER
#include"util.h"
#endif
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//static const char file[]=__FILE__;


	#define SINGLE_THREAD_CLI
//	#define _DEBUG


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
	#define CODEC_EXT c32
//	#define CODEC_EXT c33
#endif
#define STR_EXPAND(X) #X
#define STRINGIFY(X) STR_EXPAND(X)
#define GLUE_EXPAND(A, B) A##B
#define GLUE(A, B) GLUE_EXPAND(A, B)
#define CODEC_FUNC GLUE(CODEC_EXT, _codec)


int main(int argc, char **argv)
{
	int nthreads=0;
	const char *srcfn, *dstfn;
	
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif
#ifdef __GNUC__
//#ifndef _DEBUG
//#if 0
	if(argc<2||argc>4)
	{
		printf("Usage:\n");
#ifdef SINGLE_THREAD_CLI
		printf("  %s  input.ppm  output.%s  [N]    Encode file\n", argv[0], STRINGIFY(CODEC_EXT));
		printf("  %s  input.%s  output.ppm  [4]    Decode file\n", argv[0], STRINGIFY(CODEC_EXT));
		printf("N  =  1 Force CG | 2 Force WG4 | 4 Profile\n");
#else
		printf("  %s  input.ppm  output.%s  [N]    Encode file\n", argv[0], STRINGIFY(CODEC_EXT));
		printf("  %s  input.%s  output.ppm  [N]    Decode file\n", argv[0], STRINGIFY(CODEC_EXT));
		printf("\n");
		printf("Where N is an optional number of threads.\n");
		printf("  0: Use as many threads as there are cores (default).\n");
		printf("  1: Single-threaded.\n");
	//	printf("  %s  input.PPM          Test without saving\n", argv[0]);
#endif
		return 0;
	}
	srcfn=argv[1];
	dstfn=argc>=3?argv[2]:0;
	if(argc==4)
		nthreads=atoi(argv[3]);
	CODEC_FUNC(srcfn, dstfn, nthreads);
#else
	srcfn=
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13882.ppm"
	//	"C:/Projects/datasets/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/Projects/datasets/big_building.PPM"
		"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim23.ppm"
	//	"C:/Projects/datasets/kodim13-small4.PPM"
	//	"C:/Projects/datasets/kodim13-small16.PPM"
	//	"C:/Projects/datasets/dataset-CLIC30-ppm/03.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"	//large
	//	"C:/Projects/datasets/20240806 6 why me.PPM"
	//	"C:/Projects/datasets/temp.c18"

	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13912.ppm"
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
	//	"D:/ML/dataset-kodak-ppm/kodim13.c01"
	//	"C:/dataset-synthetic-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-HUGE-ppm/space_huge.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-DSLR2x4-ppm/DSC_0133.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"D:/ML/dataset-CLIC303-ppm/2048x1320_lucas-lof-388.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_cosmic-timetraveler-29758.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_rosan-harmens-18703.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-HUGE-ppm/kodak.PPM"
	//	"D:/ML/big_building.PPM"
	//	"C:/dataset-LPCB-ppm/PIA13803.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13833.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13915.ppm"	//false color terrain
	//	"C:/dataset-LPCB-ppm/STA13843.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13844.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13845.ppm"	//space clouds
		;
	dstfn=
		0
	//	"C:/Projects/datasets/temp.ppm"
	//	"C:/Projects/datasets/temp.c18"
	//	"C:/Projects/datasets/kodim13.c18"
	//	"D:/ML/dataset-kodak-ppm/kodim13.c18"
		;
//	CODEC_FUNC("D:/ML/kodim13.ppm", "D:/ML/kodim13.lsim", nthreads);
//	CODEC_FUNC("D:/ML/kodim13.lsim", "D:/ML/kodim13_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/kodim24.ppm", "D:/ML/kodim24.lsim", nthreads);
//	CODEC_FUNC("D:/ML/kodim24.lsim", "D:/ML/kodim24_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240407 blank.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240524 numbers.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240409 1 LPCB.ppm",		"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",					"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240412 2 gralic enc.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20241006 linux is cursed.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",					"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/checkboard.PPM",		"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-01.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-02.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-06.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-14.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-20.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-30.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/photo-03.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/photo-05.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/photo-49.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/zzz_halfbright.PPM",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",	"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240422 1.PPM",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/nice_clock_face.ppm",		"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm", 0, 0);//

//	CODEC_FUNC("C:/dataset-HUGE2-ppm/andromeda.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/space_huge.ppm",	"C:/Projects/datasets/mystery.lsim",	nthreads);
//	CODEC_FUNC("C:/Projects/datasets/mystery.lsim",		"C:/Projects/datasets/mystery_dec.ppm",	nthreads);

//	CODEC_FUNC("C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm",	"C:/Projects/datasets/mystery.lsim",	nthreads);
//	CODEC_FUNC("C:/Projects/datasets/mystery.lsim",				"C:/Projects/datasets/mystery_dec.ppm",	nthreads);

//	CODEC_FUNC("C:/Projects/datasets/20240513 screenshot.PPM",	"C:/Projects/datasets/mystery.lsim",	nthreads);
//	CODEC_FUNC("C:/Projects/datasets/mystery.lsim",			"C:/Projects/datasets/mystery_dec.ppm",	nthreads);

	CODEC_FUNC("C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm",	"C:/Projects/datasets/mystery.lsim",	nthreads);
	CODEC_FUNC("C:/Projects/datasets/mystery.lsim",			"C:/Projects/datasets/mystery_dec.ppm",	nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0801.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0864.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("D:/Programs/c29/song.ppm",		"D:/Programs/c29/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/Programs/c29/mystery.lsim",	"D:/Programs/c29/mystery_dec.ppm", nthreads);
//	CODEC_FUNC("D:/Programs/c29/0801.c29",		"D:/Programs/c29/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0805.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0807.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0823.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0843.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0859.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-DIV2K-ppm/0880.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/photo-52.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/photo-67.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/art.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240409 1 LPCB.ppm",		"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-meme-ppm/emoji_u1f628.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-LPCB-ppm/canon_eos_1100d_02.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-LPCB-ppm/PIA13757.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-RAW-ppm/a0014-WP_CRW_6320.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",					"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth2-ppm/20240419 3.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/jwst.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/gaia.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-01.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-GDCC2020-ppm/astro-43.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/dataset-CID22-ppm/pexels-photo-1933873.PPM",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/chaos1.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/diagram.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/blackmarble.ppm",	"D:/ML/mystery.lsim", nthreads);		//HUGE
//	CODEC_FUNC("D:/ML/mystery.lsim",			"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-memes-ppm/usa.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",		"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth-ppm/20240421 1 the front.ppm",	"D:/ML/mystery.lsim", nthreads);//synth
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-synth-ppm/20240516 4 DSC_0054.ppm",	"D:/ML/mystery.lsim", nthreads);//RCT
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-CLIC303-ppm/2048x1320_alberto-restifo-4549.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",					"D:/ML/mystery.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-a70-ppm/20240816_113656_966.ppm",	"D:/ML/mystery.lsim", nthreads);
//	CODEC_FUNC("D:/ML/mystery.lsim",				"D:/ML/mystery.ppm", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/jwst.ppm", "C:/dataset-HUGE-ppm/jwst.lsim", nthreads);
//	CODEC_FUNC("C:/dataset-HUGE-ppm/jwst.lsim", "C:/dataset-HUGE-ppm/jwst_dec.ppm", nthreads);

//	CODEC_FUNC("D:/ML/big_building.PPM", "D:/ML/big_building.LSIM", nthreads);
//	CODEC_FUNC("D:/ML/big_building.LSIM", "D:/ML/big_building_dec.PPM", nthreads);

//	CODEC_FUNC("C:/dataset-HUGE-ppm/jwst.ppm", "C:/dataset-HUGE-ppm/jwst.LSIM", nthreads);
//	CODEC_FUNC("C:/dataset-HUGE-ppm/jwst.LSIM", "C:/dataset-HUGE-ppm/jwst_dec.PPM", nthreads);




//	CODEC_FUNC("C:/Projects/datasets/kodim13.ppm", "C:/Projects/datasets/kodim13.lsim", nthreads);
//	CODEC_FUNC("C:/Projects/datasets/kodim13.lsim", "C:/Projects/datasets/kodim13_dec.ppm", nthreads);
	
//	CODEC_FUNC("C:/Projects/datasets/kodim24.ppm", "C:/Projects/datasets/kodim24.lsim", nthreads);
//	CODEC_FUNC("C:/Projects/datasets/kodim24.lsim", "C:/Projects/datasets/kodim24_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/kodim13-small16.ppm", "C:/Projects/datasets/kodim13-small16.lsim", nthreads);
//	CODEC_FUNC("C:/Projects/datasets/kodim13-small16.lsim", "C:/Projects/datasets/kodim13-small16_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm", "C:/Projects/datasets/mystery.lsim", nthreads);
//	CODEC_FUNC("C:/Projects/datasets/mystery.lsim", "C:/Projects/datasets/mystery.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm", "C:/Projects/datasets/PIA12811.lsim", nthreads);
//	CODEC_FUNC("C:/Projects/datasets/PIA12811.lsim", "C:/Projects/datasets/PIA12811_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/big_building.PPM", "C:/Projects/datasets/big_building.lsim", nthreads);		//large image
//	CODEC_FUNC("C:/Projects/datasets/big_building.lsim", "C:/Projects/datasets/big_building_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/space_huge.ppm", "C:/Projects/datasets/space_huge.lsim", nthreads);		//very large image
//	CODEC_FUNC("C:/Projects/datasets/space_huge.lsim", "C:/Projects/datasets/space_huge_dec.ppm", nthreads);

//	CODEC_FUNC("C:/Projects/datasets/20240414-noise.PPM", "C:/Projects/datasets/20240414-noise.LSIM", nthreads);	//bypass
//	CODEC_FUNC("C:/Projects/datasets/20240414-noise.LSIM", "C:/Projects/datasets/20240414-noise_dec.PPM", nthreads);

//	CODEC_FUNC(srcfn, dstfn);

//	CODEC_FUNC(srcfn, dstfn);
#endif
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
