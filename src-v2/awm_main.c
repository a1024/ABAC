//awm_main.c - AAC test
//Copyright (C) 2022  Ayman Wagih Mohsen
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include"awm_ac.h"
#include<stdio.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"

	#define	FILE_OPERATION
	#define	USE_SSE2

#define	LOUD		1

#ifndef FILE_OPERATION
#define	DATA_SIZE	1		//can be 1, 2, 4 or 8
#define IW 16
#define IH 16


#if DATA_SIZE==1
typedef unsigned char DataType;
#elif DATA_SIZE==2
typedef unsigned short DataType;
#elif DATA_SIZE==4
typedef unsigned DataType;
#elif DATA_SIZE==8
typedef unsigned long long DataType;
#else
#error Invalid data type
#endif

#define image_size	(IW*IH)
#define		ODDV		0x80
DataType buf[image_size]=
{
	//0,

	ODDV, 0x00, 0x00, 0x00, ODDV, 0x00, 0x00, 0x00,
	0x00, ODDV, 0x00, 0x00, 0x00, ODDV, 0x00, 0x00,
	0x00, 0x00, ODDV, 0x00, 0x00, 0x00, ODDV, 0x00,
	0x00, 0x00, 0x00, ODDV, 0x00, 0x00, 0x00, ODDV,
	ODDV, 0x00, 0x00, 0x00, ODDV, 0x00, 0x00, 0x00,
	0x00, ODDV, 0x00, 0x00, 0x00, ODDV, 0x00, 0x00,
	0x00, 0x00, ODDV, 0x00, 0x00, 0x00, ODDV, 0x00,
	0x00, 0x00, 0x00, ODDV, 0x00, 0x00, 0x00, ODDV,

	//0, 255,
	//255, 0,

	// 7,  7,  7,  7,  8,  8,  8,  7,
	// 7,  7,  7,  7,  7,  8,  8,  8,
	// 7,  8,  7,  7,  7,  7,  7,  8,
	// 7,  7,  8,  7,  7,  7,  7,  7,
	// 7,  7,  7,  7,  6,  7,  7,  7,
	// 7,  7,  7,  7,  6,  6,  6,  6,
	// 6,  6,  7,  7,  7,  7,  7,  7,
	// 7,  7,  7,  8,  8,  8,  8,  8,

	//0, 0, 0,
	//0, 0, 0xFF,
	//0, 0, 0xFF,

	//10,   8,  10,   9,   0,
	//22,  19,  16,  15,  -6,
	//24,  16,  22,  21,  -6,
	//24,  15,  28,  16,  -5,
	//13,   0,   5,  10,  -7,

	//182, 181, 183, 183, 185, 186, 182, 184, 191, 189, 189, 188, 190, 191, 192, 192,//iw=ih=16
	//187, 181, 184, 188, 191, 191, 191, 191, 192, 192, 191, 190, 192, 192, 193, 195,
	//191, 185, 187, 190, 191, 191, 191, 192, 193, 192, 192, 191, 191, 193, 193, 195,
	//194, 190, 189, 191, 190, 192, 192, 193, 194, 194, 194, 194, 194, 193, 193, 195,
	//193, 191, 192, 193, 192, 192, 194, 195, 196, 196, 197, 197, 197, 195, 196, 198,
	//192, 192, 192, 193, 193, 193, 194, 194, 195, 195, 195, 196, 196, 196, 198, 200,
	//193, 193, 193, 194, 195, 196, 195, 195, 195, 195, 195, 196, 196, 197, 198, 199,
	//197, 195, 194, 195, 195, 196, 196, 195, 196, 197, 196, 195, 197, 199, 200, 199,
	//199, 197, 196, 196, 196, 196, 196, 196, 196, 201, 199, 197, 199, 201, 202, 201,
	//200, 198, 197, 198, 198, 199, 199, 198, 197, 199, 198, 198, 198, 200, 201, 201,
	//200, 198, 198, 199, 199, 200, 200, 199, 198, 201, 200, 199, 200, 200, 201, 199,
	//199, 199, 199, 200, 200, 201, 201, 198, 197, 201, 200, 199, 199, 199, 200, 199,
	//200, 200, 200, 200, 201, 202, 202, 199, 199, 202, 201, 200, 200, 200, 200, 200,
	//199, 200, 200, 200, 201, 203, 203, 202, 201, 202, 201, 200, 200, 201, 201, 201,
	//200, 200, 200, 201, 202, 203, 204, 203, 203, 202, 201, 201, 201, 202, 203, 203,
	//201, 201, 202, 202, 203, 205, 205, 205, 204, 203, 202, 202, 202, 203, 204, 204,

	//182, 181, 183, 183, 185, 186, 182, 184, 191, 189, 189, 188, 190, 191, 192, 192, 191, 193, 189, 186, 186, 188, 191, 191, 191, 192, 192, 192, 193, 195, 195,//iw=31, ih=16
	//187, 181, 184, 188, 191, 191, 191, 191, 192, 192, 191, 190, 192, 192, 193, 195, 196, 191, 190, 191, 193, 196, 197, 197, 196, 200, 201, 202, 203, 203, 200,
	//191, 185, 187, 190, 191, 191, 191, 192, 193, 192, 192, 191, 191, 193, 193, 195, 194, 194, 191, 192, 193, 194, 195, 196, 195, 198, 198, 198, 199, 199, 199,
	//194, 190, 189, 191, 190, 192, 192, 193, 194, 194, 194, 194, 194, 193, 193, 195, 196, 196, 195, 195, 195, 195, 195, 196, 197, 198, 198, 198, 198, 199, 200,
	//193, 191, 192, 193, 192, 192, 194, 195, 196, 196, 197, 197, 197, 195, 196, 198, 199, 196, 196, 197, 196, 195, 195, 197, 199, 200, 200, 200, 200, 201, 201,
	//192, 192, 192, 193, 193, 193, 194, 194, 195, 195, 195, 196, 196, 196, 198, 200, 199, 196, 196, 196, 195, 195, 195, 197, 199, 199, 200, 201, 202, 202, 202,
	//193, 193, 193, 194, 195, 196, 195, 195, 195, 195, 195, 196, 196, 197, 198, 199, 199, 196, 197, 197, 197, 197, 198, 199, 200, 199, 200, 201, 202, 202, 202,
	//197, 195, 194, 195, 195, 196, 196, 195, 196, 197, 196, 195, 197, 199, 200, 199, 198, 200, 199, 199, 200, 201, 201, 201, 201, 201, 201, 202, 203, 203, 203,
	//199, 197, 196, 196, 196, 196, 196, 196, 196, 201, 199, 197, 199, 201, 202, 201, 199, 202, 202, 201, 202, 203, 204, 203, 202, 202, 202, 203, 203, 204, 205,
	//200, 198, 197, 198, 198, 199, 199, 198, 197, 199, 198, 198, 198, 200, 201, 201, 201, 202, 202, 201, 201, 201, 201, 201, 202, 199, 201, 202, 203, 204, 204,
	//200, 198, 198, 199, 199, 200, 200, 199, 198, 201, 200, 199, 200, 200, 201, 199, 199, 201, 201, 201, 201, 201, 201, 201, 201, 200, 201, 202, 203, 204, 204,
	//199, 199, 199, 200, 200, 201, 201, 198, 197, 201, 200, 199, 199, 199, 200, 199, 199, 201, 201, 202, 202, 202, 202, 202, 201, 201, 201, 201, 202, 204, 205,
	//200, 200, 200, 200, 201, 202, 202, 199, 199, 202, 201, 200, 200, 200, 200, 200, 200, 202, 203, 204, 205, 205, 204, 204, 203, 203, 202, 201, 202, 205, 206,
	//199, 200, 200, 200, 201, 203, 203, 202, 201, 202, 201, 200, 200, 201, 201, 201, 201, 201, 202, 203, 204, 204, 203, 203, 202, 203, 201, 200, 202, 204, 205,
	//200, 200, 200, 201, 202, 203, 204, 203, 203, 202, 201, 201, 201, 202, 203, 203, 202, 202, 202, 203, 203, 204, 203, 203, 203, 205, 204, 203, 204, 205, 206,
	//201, 201, 202, 202, 203, 205, 205, 205, 204, 203, 202, 202, 202, 203, 204, 204, 204, 203, 203, 203, 204, 204, 204, 204, 204, 206, 206, 206, 206, 206, 207,
	//200, 202, 202, 203, 204, 206, 206, 206, 205, 204, 204, 203, 203, 204, 204, 204, 204, 205, 204, 204, 205, 205, 205, 206, 206, 207, 207, 207, 207, 207, 208,

	//0, 0, 0, 0, 0, 1, 3, 4,
	//1, 0, 0, 0, 0, 0, 0, 0,
	//1, 0, 1, 0, 0, 0, 0, 0,
	//3, 2, 1, 0, 0, 0, 0, 0,
	//6, 5, 3, 0, 0, 0, 0, 0,
	//2, 3, 1, 0, 0, 1, 1, 0,
	//0, 1, 1, 0, 2, 5, 5, 2,
	//0, 2, 3, 1, 2, 4, 3, 0,

	//0, 0, 0, 0, 0, 1, 3, 4, 12, 7, 1, 0, 0, 1, 2, 2,
	//1, 0, 0, 0, 0, 0, 0, 0,  8, 4, 1, 0, 0, 2, 2, 0,
	//1, 0, 1, 0, 0, 0, 0, 0,  3, 1, 0, 0, 0, 0, 0, 0,
	//3, 2, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
	//6, 5, 3, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,
	//2, 3, 1, 0, 0, 1, 1, 0,  1, 0, 0, 1, 0, 0, 0, 2,
	//0, 1, 1, 0, 2, 5, 5, 2,  1, 0, 2, 4, 2, 0, 0, 4,
	//0, 2, 3, 1, 2, 4, 3, 0,  0, 0, 3, 6, 4, 0, 1, 5,
	//0, 3, 6, 4, 2, 3, 2, 0,  0, 2, 3, 1, 1, 2, 0, 0,
	//0, 2, 4, 3, 1, 1, 2, 3,  2, 3, 3, 1, 1, 3, 1, 1,
	//3, 3, 4, 4, 2, 1, 2, 6,  2, 3, 2, 0, 1, 3, 4, 4,
	//8, 4, 3, 4, 4, 2, 2, 5,  0, 0, 0, 0, 0, 2, 3, 3,
	//6, 1, 0, 2, 4, 2, 2, 3,  0, 0, 0, 0, 1, 1, 2, 2,
	//1, 0, 0, 0, 3, 3, 1, 0,  0, 0, 0, 0, 1, 0, 1, 2,
	//0, 0, 0, 0, 0, 1, 1, 0,  1, 0, 0, 0, 1, 1, 3, 4,
	//0, 2, 4, 0, 0, 0, 0, 0,  4, 1, 0, 1, 1, 1, 3, 6,

	//0x55555555, 0x55555555,
	//0x55555555, 0x55555555,

	//0xBAADF00D, 0xBAADF00D, 0xBAADF00D, 0xBAADF00D,
	//0xBAADF00D, 0xBAADF00D, 0xBAADF00D, 0xBAADF00D,
	//0xBAADF00D, 0xBAADF00D, 0xBAADF00D, 0xBAADF00D,
	//0xBAADF00D, 0xBAADF00D, 0xBAADF00D, 0xBAADF00D,

	//0xFF020100, 0xFF050403, 0xFF080706, 0xFF0B0A09,
	//0xFF0E0D0C, 0xFF11100F, 0xFF141312, 0xFF171615,
	//0xFF1A1918, 0xFF1D1C1B, 0xFF201F1E, 0xFF232221,
	//0xFF262524, 0xFF292827, 0xFF2C2B2A, 0xFF2F2E2D,

	//0x00010203, 0x04050607,
	//0x08090A0B, 0x0C0D0E0F,
};
DataType buf2[image_size]={0};
void print_vbuffer(DataType *buffer, int bw, int bh)
{
	int symbolchars=sizeof(DataType)*2+1, unspecw=bw==-1;
	if(bw*symbolchars>80)
	{
		bw=80/symbolchars;
		if(unspecw)
			bh/=bw;
	}
	if(bh*symbolchars>160)
		bh=160/symbolchars;

	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
#if DATA_SIZE==1
			printf("%02X ", (int)buffer[bw*ky+kx]);
#elif DATA_SIZE==2
			printf("%04X ", (int)buffer[bw*ky+kx]);
#elif DATA_SIZE==4
			printf("%08X ", buffer[bw*ky+kx]);
#elif DATA_SIZE==8
			printf("%016llX ", buffer[bw*ky+kx]);
#endif
		printf("\n");
	}
	printf("\n");
}
#endif
void		print_buffer(unsigned char *buffer, int bw, int bh)
{
	int symbolchars=2+1, unspecw=bh&-(bw==-1);
	if(unspecw||bw*symbolchars>80)
	{
		bw=80/symbolchars;
		if(unspecw)
			bh/=bw;
	}
	if(bh*symbolchars>160)
		bh=160/symbolchars;

	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%02X ", (int)buffer[bw*ky+kx]);
		printf("\n");
	}
	if(unspecw)
	{
		for(int k=bw*bh;k<unspecw;++k)
			printf("%02X ", (int)buffer[k]);
		printf("\n");
	}
	printf("\n");
}
void		set_console_buffer_size(short x, short y)
{
	COORD coord={x, y};
	SMALL_RECT Rect={0, 0, x-1, y-1};

    HANDLE Handle=GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleScreenBufferSize(Handle, coord);
    SetConsoleWindowInfo(Handle, TRUE, &Rect);
}


int pred2_init=0, pred2_val=0;
unsigned short pred2_prob[32]={0};
void pred2(const int *data, int iw, int ih, int bpp, int x, int y, unsigned short *hit_hist, unsigned short *ret_prob)
{
	if(!pred2_init)
	{
		pred2_init=1;
		for(int kb=0;kb<bpp;++kb)
			pred2_prob[kb]=0x8000;
	}
	else
	{
		for(int kb=0;kb<bpp;++kb)
		{
			int bit=pred2_val>>kb&1;
			pred2_prob[kb]=!bit<<15|pred2_prob[kb]>>1;//kodim23.png 0.982491	slide-19.jpg 1.462140
		}
	}
	memcpy(ret_prob, pred2_prob, bpp*sizeof(short));
	pred2_val=data[iw*y+x];
}


//const int weights[]=
//{
//	0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0,
//	1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,//kodim23.png 1.605809
//};
//const int weight_sum=0x27FFD, rx=16, ry=1;

const int weights[]=
{
	1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32767,//kodim23.png 0.966653 ???
};
const int weight_sum=0xFFFF, rx=16, ry=0;

//const int weights[]=
//{
//	0, 0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0, 0,
//	0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0,
//	1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,//kodim23.png 0.895339 ?
//};
//const int weight_sum=0x33FFB, rx=16, ry=1;

//const int weights[]=
//{
//	1, 2, 4, 8, 16, 8, 4, 2, 1,//kodim23.png 0.176958
//	2, 4, 8, 16, 32, 16, 8, 4, 2,
//	4, 8, 16, 32, 64, 32, 16, 8, 4,
//	8, 16, 32, 64, 128, 64, 32, 16, 8,
//	16, 32, 64, 128,
//};
//const int weight_sum=930, rx=4, ry=4;

//const int weights[]=
//{
//	1, 2, 3, 4, 5, 4, 3, 2, 1,//kodim23.png 0.176391
//	2, 3, 4, 5, 6, 5, 4, 3, 2,
//	3, 4, 5, 6, 7, 6, 5, 4, 3,
//	4, 5, 6, 7, 8, 7, 6, 5, 4,
//	5, 6, 7, 8,
//};
//const int weight_sum=180, rx=4, ry=4;

//const int weights[]=
//{
//	0, 0, 1, 2,  3, 2, 1, 0, 0,//kodim23.png 0.172076
//	0, 1, 2, 3,  4, 3, 2, 1, 0,
//	1, 2, 3, 4,  6, 4, 3, 2, 1,
//	2, 3, 4, 6, 12, 6, 4, 3, 2,
//	3, 4, 6, 12,
//};
//const int weight_sum=118, rx=4, ry=4;

//const int weights[]=
//{
//	0, 0, 0, 0, 0, 0, 0, 0, 0,//kodim23.png 0.108342
//	0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0,
//	1, 2, 4, 8,
//};
//const int weight_sum=15, rx=4, ry=4;

void pred1(const int *data, int iw, int ih, int bpp, int x, int y, unsigned short *hit_hist, unsigned short *ret_prob)
{
	memset(ret_prob, 0, bpp*sizeof(short));
	int ww=rx<<1|1;
	for(int kb=0;kb<bpp;++kb)
	{
		int prob=0;
		for(int ky=-ry;ky<=0;++ky)
		{
			int ky2=y+ky;
			if(ky2>=0)
			{
				for(int kx=-rx;kx<=rx;++kx)
				{
					int kx2=kx+x;
					if(kx2>=0&&kx2<iw&&(ky||kx<0))//!(!ky&&kx>=0) causal neighbor condition
					{
						int bit=data[iw*ky2+kx2]>>kb&1;
						prob+=weights[ww*(ry+ky)+rx+kx]&-!bit;
					}
				}
			}
		}
		if(prob<0)
			printf("OVERFLOW: %08X\n", prob);
		//else if(prob)
		//	prob=prob;
		unsigned long long p2=(unsigned long long)prob*0xFFFF/weight_sum;
		//printf("\t%llX\n", p2);
		if(p2>0xFFFF)
			printf("OVERFLOW: %016llX > 0xFFFF\n", p2);
		ret_prob[kb]=(unsigned short)p2;
	}
	//for(int kb=0;kb<bpp;++kb)
	//	ret_prob[kb]=(int)((unsigned long long)ret_prob[kb]*0x10000/weight_sum);
}


int			main(int argc, char **argv)
{
	//int sum=0;
	//for(int k=0;k<_countof(weights);++k)
	//	sum+=weights[k];
	set_console_buffer_size(120, 4096);
#if 0
	{
		for(int k=0;k<32;++k)
		{
			int val=(int)(10*sin((double)k/10));
			for(int k2=0;k2<6;++k2)
				printf("%c", '0'+(val>>(5-k2)&1));
			printf(" %3d\n", val);
			//printf("\n");
		}
		pause();
	}
#endif

#ifdef FILE_OPERATION
	if(argc!=2)
	{
		printf("Usage: program input_file\n\n");
		return 1;
	}
	int iw=0, ih=0, nch=0;
	unsigned char *image=stbi_load(argv[1], &iw, &ih, &nch, 4);
	if(!image)
	{
		printf("Failed to open \'%s\'\n\n", argv[1]);
		return 1;
	}
	int imsize=iw*ih;

#if 0
	ArrayHandle cbuf=0;
	//YCoCg_fwd((int*)image, iw*ih);
	test_encode((int*)image, iw, ih, nch<<3, pred1, &cbuf, 1);
	
	array_free(&cbuf);
#ifdef _DEBUG
	pause();
#endif
#endif

#if 1
	double t1, t2;
	int success=1;
	ArrayHandle cbuf;
	t1=time_ms();
	cbuf=0;
	for(int k=0;k<nch;++k)
#ifdef USE_SSE2
		success=abac0a_encode((unsigned char*)image+k, imsize, nch, &cbuf, 1);
#else
		success=abac4_encode(image, image_size, k<<3, 8, nch, &cbuf, 1);
#endif
	t2=time_ms();
	t2-=t1;
	printf("Encode: %lf ms, consumption rate: %lf MB/s\n", t2, imsize*nch/(1.024*1024*t2));
	
	unsigned char *buf2=(unsigned char*)malloc(imsize*nch);
	memset(buf2, 0, imsize*nch);
	t1=time_ms();
	const void *ptr=cbuf->data, *end=cbuf->data+cbuf->count;
	for(int k=0;k<nch;++k)
#ifdef USE_SSE2
		ptr=abac0a_decode(ptr, end, (unsigned char*)buf2+k, imsize, nch, 1);
#else
		ptr=abac4_decode(ptr, end, buf2, imsize, k<<3, 8, nch, 1);
#endif
	t2=time_ms();
	t2-=t1;
	printf("Decode: %lf ms, production rate: %lf MB/s\n", t2, imsize*nch/(1.024*1024*t2));
	
	int clean=1;
	for(int k=0;k<imsize;++k)
	{
		buf2[k]^=image[k];
		if(buf2[k])
		{
			clean=0;
			printf("Error at (%d, %d)\n", k%iw, k/iw);
			pause();
		}
	}
	if(clean)
		printf("SUCCESS\n");
	free(buf2);
	array_free(&cbuf);
	free(image);
#endif
#ifdef _MSC_VER
	//pause();
#endif
#endif
#ifndef FILE_OPERATION
	ArrayHandle cbuf;
	int success;
	const void *ptr, *end;
	double t1, t2;

	//for(int k=0;k<image_size;++k)
	//	buf[k]=rand()<RAND_MAX/8?1:0;

	//print_vbuffer(buf, IW, IH);

#if 1
	cbuf=0;
	t1=time_ms();
	for(int k=0;k<sizeof(DataType);++k)
#ifdef USE_SSE2
		success=abac0a_encode((unsigned char*)buf+k, image_size, sizeof(DataType), &cbuf, 1);
#else
		success=abac4_encode(buf, image_size, k<<3, 8, sizeof(DataType), &cbuf, 1);
#endif
	t2=time_ms();
	t2-=t1;
	printf("Encode: %lf ms, consumption rate: %lf MB/s\n", t2, image_size*sizeof(DataType)/(1.024*1024*t2));

	ptr=cbuf->data, end=cbuf->data+cbuf->count;
	t1=time_ms();
	for(int k=0;k<sizeof(DataType);++k)
#ifdef USE_SSE2
		ptr=abac0a_decode(ptr, end, (unsigned char*)buf2+k, image_size, sizeof(DataType), 1);
#else
		ptr=abac4_decode(ptr, end, buf2, image_size, k<<3, 8, sizeof(DataType), 1);
#endif
	t2=time_ms();
	t2-=t1;
	printf("Decode: %lf ms, production rate: %lf MB/s\n", t2, image_size*sizeof(DataType)/(1.024*1024*t2));
#endif

#if 0
	cbuf=0;
	for(int k=0;k<sizeof(DataType);++k)
		success=abac4_encode(buf, image_size, k<<3, 8, sizeof(DataType), &cbuf, 1);

	ptr=cbuf->data, end=cbuf->data+cbuf->count;
	for(int k=0;k<sizeof(DataType);++k)
		ptr=abac4_decode(ptr, end, buf2, image_size, k<<3, 8, sizeof(DataType), 1);
#endif

	//print_buffer(cbuf->data, -1, (int)cbuf->count);
	//print_vbuffer(buf2, IW, IH);
	
	//printf("Diff:\n");
	int clean=1;
	for(int k=0;k<image_size;++k)
	{
		buf2[k]^=buf[k];
		if(buf2[k])
		{
			clean=0;
			printf("Error at (%d, %d)\n", k%IW, k/IW);
			pause();
		}
	}
	if(clean)
		printf("SUCCESS\n");
	//print_vbuffer(buf2, IW, IH);
	
	array_free(&cbuf);
	pause();
#endif
	return 0;
}