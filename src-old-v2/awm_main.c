//awm_main.c - Compression efficiency benchmark
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
#include"stb_image.h"//for loading PNGs
#include"lodepng.h"//for saving PNGs

#define _USE_MATH_DEFINES
#include<math.h>

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


#if 1

	#define		ESTIMATE_MAGNITUDE
//	#define		CAP_NBLOCKS	2
//	#define		ONEBLOCKGUIDE

typedef float InternalType;
InternalType b2[65536];
#if defined CAP_NBLOCKS && defined ONEBLOCKGUIDE
short b4[65536];
#endif
short *b3=(short*)b2;
float mDCT8[64], mDCT16[256], mDCT32[1024];
#ifdef CAP_NBLOCKS
int nblocks=0;
#endif
#ifdef ESTIMATE_MAGNITUDE
float pmin[9]={0}, pmax[9]={0};
int peaksrccoord_x[9]={0}, peaksrccoord_y[9]={0};
#endif
void gen_name()
{
	static int callcount=0;
	sprintf_s(g_buf, G_BUF_SIZE, "out-%03d.PNG", callcount);
	++callcount;
}
void save_short(const short *buffer, int iw, int ih)
{
	gen_name();
	lodepng_encode_file(g_buf, (unsigned char*)buffer, iw, ih, LCT_GREY, 16);
}
void save_float(const float *buffer, int iw, int ih)
{
	int size=iw*ih;
	unsigned char *b2=(unsigned char*)malloc(size);
	float vmin=buffer[0], vmax=buffer[0];
	for(int k=0;k<size;++k)
	{
		if(vmin>buffer[k])
			vmin=buffer[k];
		if(vmax<buffer[k])
			vmax=buffer[k];
	}
	float gain=vmax-vmin;
	if(gain)
	{
		gain=255/gain;
		for(int k=0;k<size;++k)
			b2[k]=(unsigned char)((buffer[k]-vmin)*gain);
	}
	else
		memset(b2, 0, size);
	gen_name();
	lodepng_encode_file(g_buf, b2, iw, ih, LCT_GREY, 8);
	free(b2);
}
int try_encode(ArrayHandle *out, InternalType *src, int iw, int ih, int x, int y, int lgdx, int lgdy, float vmin, float vmax)//this assumes the src is padded
{
	int bw=1<<lgdx, bh=1<<lgdy, nstages=lgdx<lgdy?lgdx:lgdy, size=bw*bh;
	unsigned long long hist[64]={0};
	int nbits=0;
	size_t outsize;
	int ky, ylim=bh<ih-y?bh:ih-y, xlim=bw<iw-x?bw:iw-x;
	unsigned char *flag;

	for(ky=0;ky<ylim;++ky)//read block (zero-padded)
	{
		memcpy(b2+bw*ky, src+iw*(y+ky)+x, xlim*sizeof(InternalType));
		if(xlim<bw)
			memset(b2+bw*ky+xlim, 0, (bw-xlim)*sizeof(InternalType));
	}
	if(ky<bh)
		memset(b2+bw*ky, 0, (bh-ky)*bw*sizeof(InternalType));
	outsize=out[0]->count;
	switch(nstages)
	{
	case 2://bypass 8x8
		{
			flag=(unsigned char*)ARRAY_APPEND(*out, 0, 65, 1, 0);
			flag[0]=2;
			for(int k=0;k<64;++k)
			{
				float val=b2[k];
				if(val<0)
					val=0;
				if(val>255)
					val=255;
				flag[k+1]=(unsigned char)val;
			}
		}
		return (int)(out[0]->count-outsize);

		//try DCT		FIXME support for nonsquare blocks
	case 3:
		apply_DCT_2D(mDCT8, mDCT8, b2, 8, 8, b2+32486);
		break;
	case 4:
		apply_DCT_2D(mDCT16, mDCT16, b2, 16, 16, b2+32486);
		break;
	case 5:
		apply_DCT_2D(mDCT32, mDCT32, b2, 32, 32, b2+32486);
		break;

	case 6:case 7:case 8://try DWT
		dwt2_2d_fwd(b2, bw, bh, nstages-2);
		break;
	}
	vmin=vmax=b2[0];
	for(int k=0;k<size;++k)//quantize & analyze
	{
		float fval;
		short sval;
		unsigned short code, neg;

		fval=b2[k];
		//fval=sqrtf(b2[k]*10);
#ifdef ESTIMATE_MAGNITUDE
		if(vmin>fval)//
			vmin=fval;
		if(vmax<fval)
			vmax=fval;
#endif
		if(fval>40000)
			fval=fval;
		
		//sval=(short)fval;
		sval=(short)(fval*0.1f);

		neg=sval<0;
		sval^=-neg;
		sval+=neg;
		sval<<=1;
		sval|=neg;

		code=b3[k]=sval;
		if(k)
			hist[code>>10]|=1LL<<(code>>4&0x3F);
	}
#ifdef ESTIMATE_MAGNITUDE
	if(pmin[nstages]>vmin)//
		pmin[nstages]=vmin;
	if(pmax[nstages]<vmax)
		pmax[nstages]=vmax, peaksrccoord_x[nstages]=x, peaksrccoord_y[nstages]=y;
#endif

	for(int k=0;k<64*64;++k)
		nbits+=hist[k>>6]>>(k&63)&1;
	if(nbits>32)//test fails if more than [quarter] levels was used
		return 0;

	//save_short(b3, bw, bh);//
#ifdef ONEBLOCKGUIDE
	if(!nblocks)//
		memcpy(b4, b3, size*sizeof(short));//
#endif

	flag=(unsigned char*)ARRAY_APPEND(*out, 0, 1, 1, 0);
	*flag=nstages;
	int csize=abac4_encode(b3, size, 0, 12, 2, out, 0);//entropy coding
	
	printf("Encode %dx%d, csize=%d, ratio=%lf, %f~%f\n", bw, bh, csize, (double)size/csize, vmin, vmax);//
	//printf("Encode %dx%d, csize=%d, ratio=%lf\n", bw, bh, csize, (double)size/csize);//
	return (int)(out[0]->count-outsize);
}
void encode_buffer(ArrayHandle *out, InternalType *buffer, int bw, int bh, float vmin, float vmax)
{
	int ds;//change in size
	for(int ky=0;ky<bh;ky+=256)
	{
		for(int kx=0;kx<bw;kx+=256)
		{
			//printf("\r");
			ds=try_encode(out, buffer, bw, bh, kx, ky, 8, 8, vmin, vmax);
#ifdef CAP_NBLOCKS
			if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
			if(!ds)
			{
				for(int k7=0;k7<4;++k7)
				{
					int x7=kx+((k7&1)<<7), y7=ky+((k7>>1)<<7);
					ds=try_encode(out, buffer, bw, bh, x7, y7, 7, 7, vmin, vmax);
#ifdef CAP_NBLOCKS
					if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
					if(!ds)
					{
						for(int k6=0;k6<4;++k6)
						{
							int x6=x7+((k6&1)<<6), y6=y7+((k6>>1)<<6);
							ds=try_encode(out, buffer, bw, bh, x6, y6, 6, 6, vmin, vmax);
#ifdef CAP_NBLOCKS
							if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
							if(!ds)
							{
								for(int k5=0;k5<4;++k5)
								{
									int x5=x6+((k5&1)<<5), y5=y6+((k5>>1)<<5);
									ds=try_encode(out, buffer, bw, bh, x5, y5, 5, 5, vmin, vmax);
#ifdef CAP_NBLOCKS
									if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
									if(!ds)
									{
										for(int k4=0;k4<4;++k4)
										{
											int x4=x5+((k4&1)<<4), y4=y5+((k4>>1)<<4);
											ds=try_encode(out, buffer, bw, bh, x4, y4, 4, 4, vmin, vmax);
#ifdef CAP_NBLOCKS
											if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
											if(!ds)
											{
												for(int k3=0;k3<4;++k3)
												{
													int x3=x4+((k3&1)<<3), y3=y4+((k3>>1)<<3);
													ds=try_encode(out, buffer, bw, bh, x3, y3, 3, 3, vmin, vmax);
#ifdef CAP_NBLOCKS
													if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
													if(!ds)//bypass
													{
														ds=try_encode(out, buffer, bw, bh, x3, y3, 2, 2, vmin, vmax);
#ifdef CAP_NBLOCKS
														if(ds){++nblocks; if(nblocks>=CAP_NBLOCKS)return;}
#endif
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

#if 0
			for(int ky2=0;ky2<256;++ky2)
				memcpy(b2+(ky2<<8), buffer+bw*ky+kx, 256*sizeof(InternalType));

			//analyze block
			memcpy(b3, b2, 65536*sizeof(InternalType));
			dwt2_2d_fwd(b3, 256, 256, 8);
			vmin=b3[1], vmax=b3[1];
			for(int k=2;k<65536;++k)
			{
				if(vmin>b3[k])
					vmin=b3[k];
				if(vmax<b3[k])
					vmax=b3[k];
			}
			if(vmax-vmin<4)//pass the test
#endif
		}
	}
	//printf("\n");
}
const unsigned char *srcstart=0;
int decode_block(const unsigned char **srcptr, const unsigned char *srcend, float *dst, int iw, int ih, int x, int y, int lgdim)//recursive
{
	if(**srcptr<lgdim)
	{
		for(int k=0;k<4;++k)
		{
			if(!decode_block(srcptr, srcend, dst, iw, ih, x+(1<<(lgdim-1)&-(k&1)), y+(1<<(lgdim-1)&-(k>>1)), lgdim-1))
				return 0;
#ifdef CAP_NBLOCKS
			if(!nblocks)
				return 1;
#endif
		}
	}
	else//decode dim x dim block at (x, y)
	{
		++*srcptr;
		if(lgdim==2)//bypass 8x8
		{
			for(int ky=0;ky<8;++ky)
			{
				for(int kx=0;kx<8&&*srcptr<srcend;++kx, ++*srcptr)
					dst[iw*(y+ky)+x+kx]=**srcptr;
			}
		}
		else
		{
			int dim=1<<lgdim, size=dim*dim;
			memset(b3, 0, size*sizeof(short));
			*srcptr=(unsigned char*)abac4_decode(*srcptr, srcend, b3, size, 0, 12, 2, 0);
#ifdef ONEBLOCKGUIDE
			if(nblocks==CAP_NBLOCKS)
			{
				for(int k=0;k<size;++k)
				{
					if(b3[k]!=b4[k])
					{
						printf("GUIDE: Decode error at %d: 0x%04X != original 0x%04X\n", k, b3[k], b4[k]);
						return 0;
					}
				}
			}
#endif
			if(!*srcptr)
			{
				printf("Failed to decode %dx%d block at (%d, %d)\n", dim, dim, x, y);
				return 0;
			}
			//if(*srcptr<srcend&&(**srcptr<2||**srcptr>8))
			//{
			//	printf("Invalid next block flag at %d\n", (int)(*srcptr-srcstart));
			//	return 0;
			//}
			//if(*srcptr<srcend&&!(srcptr[0][0]=='A'&&srcptr[0][0]=='C'&&srcptr[0][0]=='0'&&srcptr[0][0]=='4'))//abac4
			//{
			//	printf("Invalid next block tag at %d\n", (int)(*srcptr-srcstart));
			//	return 0;
			//}
			for(int k=size-1;k>=0;--k)//dequantize
			{
				short sval;
				unsigned short neg;

				sval=b3[k];

				neg=sval&1;
				sval>>=1;
				sval^=-neg;
				sval+=neg;

				b2[k]=(float)(10*sval);
			}
			switch(lgdim)
			{
			case 3://DCT 8x8
				apply_DCT_2D(mDCT8, mDCT8, b2, 8, 8, b2+32486);
				break;
			case 4://DCT 16x16
				apply_DCT_2D(mDCT16, mDCT16, b2, 16, 16, b2+32486);
				break;
			case 5://DCT 32x32
				apply_DCT_2D(mDCT32, mDCT32, b2, 32, 32, b2+32486);
				break;
			case 6://DWT 64x64
			case 7://DWT 128x128
			case 8://DWT 256x256
				dwt2_2d_inv(b2, dim, dim, lgdim-2);
				break;
			}
			int ky, ylim=dim<ih-y?dim:ih-y, xlim=dim<iw-x?dim:iw-x;
			for(ky=0;ky<ylim;++ky)//read block (zero-padded)
				memcpy(dst+iw*(y+ky)+x, b2+dim*ky, xlim*sizeof(InternalType));
		}
	}
#ifdef CAP_NBLOCKS
	--nblocks;
#endif
	return 1;
}
int decode_buffer(const unsigned char **srcptr, const unsigned char *srcend, float *dst, int iw, int ih)
{
	for(int ky=0;ky<ih;ky+=256)
	{
		for(int kx=0;kx<iw;kx+=256)
		{
			if(!decode_block(srcptr, srcend, dst, iw, ih, kx, ky, 8))
				return 0;
#ifdef CAP_NBLOCKS
			if(!nblocks)
				return 1;
#endif
		}
	}
	return 1;
}
void test2(int argc, char **argv)
{
	int iw, ih, nch;
	int *image;
	int size;
	InternalType *bufY, *bufCo, *bufCg;
	ArrayHandle out;

	if(argc!=2)
	{
		printf("Usage: program input_file\n\n");
		exit(1);
	}
	image=(int*)stbi_load(argv[1], &iw, &ih, &nch, 4);
	if(!image)
	{
		printf("Cannot open \'%s\'\n\n", argv[1]);
		exit(1);
	}

	printf("Encoding \'%s\'...\n", argv[1]);//encode
	init_fwdDCT(mDCT8, 3);
	init_fwdDCT(mDCT16, 4);
	init_fwdDCT(mDCT32, 5);
	size=iw*ih;
	bufY=(InternalType*)malloc(size*sizeof(InternalType));
	bufCo=(InternalType*)malloc(size*sizeof(InternalType));
	bufCg=(InternalType*)malloc(size*sizeof(InternalType));
	ARRAY_ALLOC(unsigned char, out, 0, 0, 0, 0);
	
	//color transform
	YCoCg_f32_fwd(image, size, bufY, bufCo, bufCg);
	encode_buffer(&out, bufY, iw, ih, 0, 255);
#ifndef CAP_NBLOCKS
	encode_buffer(&out, bufCo, iw, ih, -255, 255);
	encode_buffer(&out, bufCg, iw, ih, -255, 255);
#endif

#ifdef ESTIMATE_MAGNITUDE
	printf("Ranges:\n");
	for(int k=0;k<9;++k)
		printf("%d\t%f ~ %f at (%d, %d)\n", k, pmin[k], pmax[k], peaksrccoord_x[k], peaksrccoord_y[k]);
	printf("\n");
#endif
#if 0
	const char hello[]="helloooooooooooooooooooo woooooooooooooooooooooorld!!!!!!!!!!!!";
	int blocksize=sizeof(hello);
	for(int k=0;k<4;++k)
	{
		int csize=abac4_encode(hello, blocksize, 0, 8, 1, &out, 0);
		printf("Encoded %d -> %d\n", blocksize, csize);
	}
#endif
#if 0
	//long long tail=0xEFCDAB8967452301;
	//ARRAY_APPEND(out, &tail, 8, 1, 0);
	for(int k=0;k<(int)out->count;++k)
		printf("%02X", out->data[k]&0xFF);
	printf("\n");
#endif

	printf("%dx%dx3 = %d -> %d, ratio = %lf\n", iw, ih, size*3, (int)out->count, (double)size*3/out->count);
	
	printf("Decoding...\n");//decode
	init_invDCT(mDCT8, 3);
	init_invDCT(mDCT16, 4);
	init_invDCT(mDCT32, 5);

	const unsigned char *srcptr=out->data, *srcend=out->data+out->count;
	srcstart=out->data;
	int success;
	success=decode_buffer(&srcptr, srcend, bufY, iw, ih);
	if(!success)
		goto quit;
#ifndef CAP_NBLOCKS
	success=decode_buffer(&srcptr, srcend, bufCo, iw, ih);
	if(!success)
		goto quit;
	success=decode_buffer(&srcptr, srcend, bufCg, iw, ih);
	if(!success)
		goto quit;
#endif
	printf("Inverse color transform...\n");
	YCoCg_f32_inv(image, size, bufY, bufCo, bufCg);

	printf("Saving...\n");
	lodepng_encode_file("out.PNG", (unsigned char*)image, iw, ih, LCT_RGBA, 8);
	printf("SUCCESS?\n");
	//spatial decorrelating transforms
	//quantization
	//entropy coding
	//dequantization & inverse spatial transforms
	//inverse channel transform
	//save out.PNG
quit:
	pause();
	exit(0);
}
void test3()
{
	const char hello[]="helloooooooooooooooooooo woooooooooooooooooooooorld!!!!!!!!!!!!";
	int blocksize=sizeof(hello);
	ArrayHandle out;
	int csize;
	unsigned char *buf;
	const unsigned char *srcptr, *srcend;
	int ntimes=4;

	ARRAY_ALLOC(unsigned char, out, 0, 0, 0, 0);
	for(int k=0;k<ntimes;++k)
	{
		csize=abac4_encode(hello, blocksize, 0, 8, 1, &out, 0);
		printf("Encoded %d -> %d\n", blocksize, csize);
	}

	buf=(unsigned char*)malloc(blocksize);
	memset(buf, 0, blocksize);
	srcptr=out->data;
	srcend=out->data+out->count;
	for(int k=0;k<ntimes;++k)
	{
		srcptr=(unsigned char*)abac4_decode(srcptr, srcend, buf, blocksize, 0, 8, 1, 0);
		if(!srcptr)
		{
			printf("Decode failure at block %d\n", k);
			goto failure;
		}
		int clean=1;
		for(int k2=0;k2<blocksize;++k2)
		{
			if(buf[k2]!=hello[k2])
			{
				printf("Decode failure at block %d at sym %d: 0x%02X != original 0x%02X\n", k, k2, buf[k2]&0xFF, hello[k2]&0xFF);
				goto failure;
			}
		}
	}
	printf("SUCCESS\n");
failure:
	free(buf);
	pause();
	exit(0);
}
void test4()
{
	memset(b2, 0, sizeof(b2));
	b2[0]=1000;
	//dwt2_2d_fwd(b2, 256, 256, 6);
	dwt2_2d_inv(b2, 256, 256, 6);
	save_float(b2, 256, 256);

	pause();
	exit(0);
}
#endif
#if 0
void print_buffer_double(double *data, int count)
{
	for(int k=0;k<count;++k)
		printf("%g ", data[k]);
	printf("\n");
}
void DCT8f_fwd(double *data);
void DCT8f_inv(double *data);
void DCT8_ref(double *data, int inv);
void DCT_test()
{
	double data[]=
	{
		1, 2, 3, 4, 3, 2, 2, 2,
		//0.21, 0.36, 0.51, 0.14, 0.99, 0.01, 0.73, 0.28,
	};
	double d2[8];
	long long t1, t2;
	printf("Original:\n"), print_buffer_double(data, 8);
	printf("\n");

	t1=__rdtsc();
	DCT8_ref(data, 0);
	t2=__rdtsc();
	printf("Forward Reference: %lldc\n", t2-t1), print_buffer_double(data, 8);
	memcpy(d2, data, 8*sizeof(double));
	t1=__rdtsc();
	DCT8_ref(data, 1);
	t2=__rdtsc();
	printf("Inverse Reference: %lldc\n", t2-t1), print_buffer_double(data, 8);
	printf("\n");

	t1=__rdtsc();
	DCT8f_fwd(data);
	t2=__rdtsc();
	printf("Forward Optimized: %lldc\n", t2-t1), print_buffer_double(data, 8);
	for(int k=0;k<8;++k)
		d2[k]/=data[k];
	printf("Reference/Optimized:\n"), print_buffer_double(d2, 8);
	t1=__rdtsc();
	DCT8f_inv(data);
	t2=__rdtsc();
	printf("Inverse Optimized: %lldc\n", t2-t1), print_buffer_double(data, 8);
	printf("\n");
	
	t1=__rdtsc();
	DCT8f_inv(data);
	t2=__rdtsc();
	printf("Inverse Optimized: %lldc\n", t2-t1), print_buffer_double(data, 8);
	t1=__rdtsc();
	DCT8f_fwd(data);
	t2=__rdtsc();
	printf("Forward Optimized: %lldc\n", t2-t1), print_buffer_double(data, 8);

	pause();
	exit(0);
}
#endif


int			main(int argc, char **argv)
{
	set_console_buffer_size(120, 4096);

	//DCT_test();//
//	test4();
//	test3();
	test2(argc, argv);

	//int sum=0;
	//for(int k=0;k<_countof(weights);++k)
	//	sum+=weights[k];
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
	//YCoCg_i32_fwd((int*)image, iw*ih);
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
#ifdef _DEBUG
	pause();
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