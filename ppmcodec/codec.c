#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
static const char file[]=__FILE__;


//	#define _DEBUG


#ifndef CODEC_FUNC
#define CODEC_FUNC c22_codec
#endif


static void ppm_skip(const unsigned char **ptr, const unsigned char *end)
{
	for(;*ptr<end;)
	{
		if(**ptr=='#')
		{
			for(;*ptr<end&&**ptr!='\n';++*ptr);
		}
		else if(**ptr<'!')
			++*ptr;
		else
			return;
	}
}
int header_read(const unsigned char *src, int len, int *iw, int *ih, CodecID *codec)
{
	const unsigned char *ptr, *end;
	int vmax;

	ptr=src;
	end=src+len;
	if(!memcmp(ptr, "P6", 2))
	{
		*codec=CODEC_PPM;
		ptr+=2;
	}
	else if(!memcmp(ptr, "P5", 2))
	{
		*codec=CODEC_PGM;
		ptr+=2;
	}
	else if(!memcmp(ptr, "C01", 3))
	{
		*codec=CODEC_C01;
		ptr+=3;
	}
	ppm_skip(&ptr, end);
	*iw=strtol((char*)ptr, (char**)&ptr, 10);
	ppm_skip(&ptr, end);
	*ih=strtol((char*)ptr, (char**)&ptr, 10);
	if(*codec==CODEC_PPM||*codec==CODEC_PGM)
	{
		ppm_skip(&ptr, end);
		vmax=strtol((char*)ptr, (char**)&ptr, 10);//255
		if(vmax!=255)
			*codec=CODEC_INVALID;
	}
	ptr+=*ptr=='\n';//ppm_skip skips binary data < 33 as well
	return (int)(ptr-src);
}
int compare_bufs_8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud)
{
	ptrdiff_t len=(ptrdiff_t)chstride*iw*ih;
	int inc=chstride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-chstride:0;k>=0&&k<len;k+=inc)
	{
		if(memcmp(b1+k, b0+k, nch*sizeof(char)))
		{
			if(loud)
			{
				ptrdiff_t idx=k/chstride, kx=idx%iw, ky=idx/iw;
				printf("\n%s error IDX %td  XY (%5td, %5td) / %5d x %5d  b1 != b0\n", name, k, kx, ky, iw, ih);
				for(int kc=0;kc<nch;++kc)
				{
					char c=(unsigned char)b1[k+kc]==(unsigned char)b0[k+kc]?'=':'!';
					printf("C%d  0x%04X %c= 0x%04X    %d %c= %d\n",
						kc,
						(unsigned char)b1[k+kc], c, (unsigned char)b0[k+kc],
						b1[k+kc], c, b0[k+kc]
					);
				}
				LOG_ERROR("");
			}
			return 1;
		}
	}
	if(loud)
		printf("%s:\tSUCCESS\n", name);
	return 0;
}
int main(int argc, char **argv)
{
	const char *srcfn, *dstfn;
	
#ifdef __GNUC__
//#ifndef _DEBUG
//#if 0
	if((unsigned)(argc-2)>1)
	{
		printf("Usage:\n");
	//	printf("  %s  input.PPM          Test without saving\n", argv[0]);
		printf("  %s  input.PPM  output  Encode file\n", argv[0]);
		printf("  %s  input  output.PPM  Decode file\n", argv[0]);
		return 0;
	}
	srcfn=argv[1];
	dstfn=argc==3?argv[2]:0;
	CODEC_FUNC(srcfn, dstfn);
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
//	c23_codec("D:/ML/kodim13.ppm", "D:/ML/kodim13.lsim");
//	c23_codec("D:/ML/kodim13.lsim", "D:/ML/kodim13_dec.ppm");

//	c22_codec("D:/ML/kodim24.ppm", "D:/ML/kodim24.lsim");
//	c22_codec("D:/ML/kodim24.lsim", "D:/ML/kodim24_dec.ppm");

	c24_codec("C:/dataset-memes-ppm/1725754154794738.ppm",	"D:/ML/mystery.lsim");//synth
	c24_codec("D:/ML/mystery.lsim",				"D:/ML/mystery_dec.ppm");

//	c24_codec("C:/dataset-synth-ppm/20240421 1 the front.ppm",	"D:/ML/mystery.lsim");//synth
//	c24_codec("D:/ML/mystery.lsim",					"D:/ML/mystery_dec.ppm");

//	c24_codec("C:/dataset-synth-ppm/20240516 4 DSC_0054.ppm",	"D:/ML/mystery.lsim");//RCT
//	c24_codec("D:/ML/mystery.lsim",					"D:/ML/mystery_dec.ppm");

//	c22_codec("C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm",	"D:/ML/mystery.lsim");
//	c22_codec("D:/ML/mystery.lsim",					"D:/ML/mystery.ppm");

//	c21_codec("C:/dataset-CLIC303-ppm/2048x1320_alberto-restifo-4549.ppm",	"D:/ML/mystery.lsim");
//	c21_codec("D:/ML/mystery.lsim",						"D:/ML/mystery.ppm");

//	c22_codec("C:/dataset-a70-ppm/20240816_113656_966.ppm",	"D:/ML/mystery.lsim");
//	c22_codec("D:/ML/mystery.lsim",				"D:/ML/mystery.ppm");

//	c24_codec("C:/dataset-HUGE-ppm/jwst.ppm", "C:/dataset-HUGE-ppm/jwst.lsim");
//	c24_codec("C:/dataset-HUGE-ppm/jwst.lsim", "C:/dataset-HUGE-ppm/jwst_dec.ppm");

//	c10_codec("D:/ML/big_building.PPM", "D:/ML/big_building.LSIM");
//	c10_codec("D:/ML/big_building.LSIM", "D:/ML/big_building_dec.PPM");

//	c01_codec("C:/dataset-HUGE-ppm/jwst.ppm", "C:/dataset-HUGE-ppm/jwst.LSIM");
//	c01_codec("C:/dataset-HUGE-ppm/jwst.LSIM", "C:/dataset-HUGE-ppm/jwst_dec.PPM");




//	c20_codec("C:/Projects/datasets/kodim13.ppm", "C:/Projects/datasets/kodim13.lsim");
//	c20_codec("C:/Projects/datasets/kodim13.lsim", "C:/Projects/datasets/kodim13_dec.ppm");
	
//	c20_codec("C:/Projects/datasets/kodim24.ppm", "C:/Projects/datasets/kodim24.lsim");
//	c20_codec("C:/Projects/datasets/kodim24.lsim", "C:/Projects/datasets/kodim24_dec.ppm");

//	c08_codec("C:/Projects/datasets/kodim13-small16.ppm", "C:/Projects/datasets/kodim13-small16.lsim");
//	c08_codec("C:/Projects/datasets/kodim13-small16.lsim", "C:/Projects/datasets/kodim13-small16_dec.ppm");

//	c21_codec("C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm", "C:/Projects/datasets/mystery.lsim");
//	c21_codec("C:/Projects/datasets/mystery.lsim", "C:/Projects/datasets/mystery.ppm");

//	c20_codec("C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm", "C:/Projects/datasets/PIA12811.lsim");
//	c20_codec("C:/Projects/datasets/PIA12811.lsim", "C:/Projects/datasets/PIA12811_dec.ppm");

//	c20_codec("C:/Projects/datasets/big_building.PPM", "C:/Projects/datasets/big_building.lsim");		//large image
//	c20_codec("C:/Projects/datasets/big_building.lsim", "C:/Projects/datasets/big_building_dec.ppm");

//	c20_codec("C:/Projects/datasets/space_huge.ppm", "C:/Projects/datasets/space_huge.lsim");		//very large image
//	c20_codec("C:/Projects/datasets/space_huge.lsim", "C:/Projects/datasets/space_huge_dec.ppm");

//	c08_codec("C:/Projects/datasets/20240414-noise.PPM", "C:/Projects/datasets/20240414-noise.LSIM");	//bypass
//	c08_codec("C:/Projects/datasets/20240414-noise.LSIM", "C:/Projects/datasets/20240414-noise_dec.PPM");

//	c07_codec(srcfn, dstfn);

//	CODEC_FUNC(srcfn, dstfn);
#endif
	return 0;
}
