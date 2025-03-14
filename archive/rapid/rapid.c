﻿#include"rapid.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
static const char file[]=__FILE__;

#define CODEC_FUNC r01_codec

const char *progname=0;
static void print_usage(const char *progname)
{
	printf(
		"Rapid Lossless Image Codec (RaLIC) by AWM (2024)\n"
		"Usage:\n"
		"  %s  c  input.PPM   output.LSIM  Encode silently\n"
		"  %s  d  input.LSIM  output.PPM   Decode silently\n"
		"  %s  C  input.PPM   output.LSIM  Encode and print summary\n"
		"  %s  D  input.LSIM  output.PPM   Decode and print summary\n"
		, progname, progname, progname, progname
	);
}
int main(int argc, char **argv)
{
	progname=argv[0];
//#ifdef __GNUC__
#ifndef _DEBUG
	if(argc!=4||strlen(argv[1])!=1)
	{
		print_usage(argv[0]);
		return 0;
	}
	switch(argv[1][0])
	{
	case 'c':CODEC_FUNC(argv[2], argv[3], 1, 0);break;
	case 'd':CODEC_FUNC(argv[2], argv[3], 0, 0);break;
	case 'C':CODEC_FUNC(argv[2], argv[3], 1, 1);break;
	case 'D':CODEC_FUNC(argv[2], argv[3], 0, 1);break;
	default:
		print_usage(argv[0]);
		break;
	}
#else
	CODEC_FUNC("C:/Projects/datasets/kodim13.ppm", "C:/Projects/datasets/kodim13.lsim", 1, 1);
	CODEC_FUNC("C:/Projects/datasets/kodim13.lsim", "C:/Projects/datasets/kodim13_dec.ppm", 0, 1);

//	CODEC_FUNC("C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm", "C:/Projects/datasets/canon_eos_1100d_01.lsim", 1, 1);//X  out of resources
//	CODEC_FUNC("C:/Projects/datasets/canon_eos_1100d_01.lsim", "C:/Projects/datasets/canon_eos_1100d_01-dec.ppm", 0, 1);

//	CODEC_FUNC("C:/Projects/datasets/space_huge.ppm", "C:/Projects/datasets/space_huge.lsim", 1, 1);		//X  out of resources	screen flicker
//	CODEC_FUNC("C:/Projects/datasets/space_huge.lsim", "C:/Projects/datasets/space_huge_dec.ppm", 0, 1);
#endif
	return 0;
}