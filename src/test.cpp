#include"ac.h"
#include"huffman.h"
#include"lz77.h"

#include"lodepng.h"
#ifndef __linux__
#include<Windows.h>
#endif
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<string>
#include<fstream>
#ifdef __linux__
#define	__rdtsc	__builtin_ia32_rdtsc
#define scanf_s scanf
#define sprintf_s snprintf
#else
#include<intrin.h>
#include<conio.h>
#endif
#include<time.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"

//	#define		DEBUG_DWT
//	#define		ENABLE_RVL

#define	SIZEOF(STATIC_ARRAY)	(sizeof(STATIC_ARRAY)/sizeof(*STATIC_ARRAY))

#define			G_BUF_SIZE	1024
char			g_buf[G_BUF_SIZE]={};
#ifdef __linux__
typedef unsigned char byte;
#define			set_console_buffer_size(...)
#else
int				set_console_buffer_size(short w, short h)
{
	COORD coords={w, h};
	HANDLE handle=GetStdHandle(STD_OUTPUT_HANDLE);
	int success=SetConsoleScreenBufferSize(handle, coords);
	if(!success)
		printf("Failed to resize console buffer: %d\n\n", GetLastError());
	return success;
}
#endif
bool			open_text(const char *filename, std::string &data)
{
	std::ifstream input(filename);
	if(!input.is_open())
		return false;
	data.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
	input.close();
	return true;
}

void			print_hex(short *buffer, int imsize, int depth)
{
	for(int k=0;k<imsize;++k)
		printf("%0*X-", (depth+3)>>2, buffer[k]);
	printf("\n");
}
void			print_rgba8(unsigned *data, int npixels)
{
	printf("AABBGGRR:\n");
	for(int k=0;k<npixels;++k)
		printf("%08X-", data[k]);
	printf("\n");
}
void			print_bin(int x, int nbits)
{
	//printf("0b");
	for(int k=nbits-1;k>=0;--k)
	{
		printf("%d", x>>k&1);
		if(k&&!(k&3))
			printf("_");
	}
}
void			print_bitplane(const void *src, int imsize, int bytestride, int bitplane)
{
	auto buffer=(byte*)src;
	int bit_offset=bitplane>>3, bit_shift=bitplane&7;
	for(int k=0, k2=0;k<imsize;++k, k2+=bytestride)
		printf("%d", buffer[k2+bit_offset]>>bit_shift&1);
	printf("\n");
}
void			print_image(const int *buffer, int bw, int bh, int x1, int x2, int y1, int y2)
{
	if(x1<0)
		x1=0;
	if(x2>bw)
		x2=bw;
	if(y1<0)
		y1=0;
	if(y2>bh)
		y2=bh;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
			printf("%08X-", buffer[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
}
void			print_fdata(const float *data, int bw, int bh, int x1, int x2, int y1, int y2)
{
	if(x1<0)
		x1=0;
	if(x2>bw)
		x2=bw;
	if(y1<0)
		y1=0;
	if(y2>bh)
		y2=bh;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
			printf("%f ", data[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
}
void			print_data(const float *data, int count, int stride)
{
	for(int k=0, end=count*stride;k<end;k+=stride)
		printf("%f ", data[k]);
	printf("\n");
}
void			print_diff(const int *buffer, const int *b2, int bw, int bh, int x1, int x2, int y1, int y2, int mask)
{
	if(x1<0)
		x1=0;
	if(x2>bw)
		x2=bw;
	if(y1<0)
		y1=0;
	if(y2>bh)
		y2=bh;
	auto m=(const byte*)&mask;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			auto p=(const byte*)(buffer+bw*ky+kx), p2=(const byte*)(b2+bw*ky+kx);
			printf("%02X%02X%02X%02X-", abs((p[3]&m[3])-(p2[3]&m[3])), abs((p[2]&m[2])-(p2[2]&m[2])), abs((p[1]&m[1])-(p2[1]&m[1])), abs((p[0]&m[0])-(p2[0]&m[0])));
		}
			//printf("%08X-", abs(buffer[bw*ky+kx]-b2[bw*ky+kx]));
		printf("\n");
	}
	printf("\n");
}

//old transforms
void			differentiate_rows(const short *src, int bw, int bh, short *dst)//can be the same buffer
{
	for(int ky=0;ky<bh;++ky)
	{
		dst[bw*ky]=src[bw*ky];
		for(int kx=bw-1;kx>=1;--kx)
		{
			int delta=src[bw*ky+kx]-src[bw*ky+kx-1];
			int neg=delta<0;
			delta^=-neg, delta+=neg;
			delta=delta<<1|neg;//RVL
			dst[bw*ky+kx]=delta;
		}
	}
}
void			integrate_rows(const short *src, int bw, int bh, short *dst)
{
	for(int ky=0;ky<bh;++ky)
	{
		dst[bw*ky]=src[bw*ky];
		for(int kx=1;kx<bw;++kx)
		{
			int delta=src[bw*ky+kx];
			int neg=delta&1;
			delta>>=1;
			delta^=-neg, delta+=neg;
			delta+=src[bw*ky+kx-1];
			dst[bw*ky+kx]=delta;
		}
	}
}

void			RVL_rows(const int *src, int bw, int bh, int *dst)//can be the same buffer		24 -> 27bit
{
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=bw-1;kx>=1;--kx)
		{
			auto p=(const byte*)(src+bw*ky+kx), p0=(const byte*)(src+bw*ky+kx-1);
			int R=p[0]-p0[0], G=p[1]-p0[1], B=p[2]-p0[2], neg;
			neg=R<0, R^=-neg, R+=neg, R=R<<1|neg;
			neg=G<0, G^=-neg, G+=neg, G=G<<1|neg;
			neg=B<0, B^=-neg, B+=neg, B=B<<1|neg;
			dst[bw*ky+kx]=B<<18|G<<9|R;
		}
		auto ps=(const byte*)(src+bw*ky);
		dst[bw*ky]=ps[2]<<18|ps[1]<<9|ps[0];
	}
}
void			invRVL_rows(const int *src, int bw, int bh, int *dst)//27 -> 24bit
{
	for(int ky=0;ky<bh;++ky)
	{
		int ps=src[bw*ky];
		dst[bw*ky]=(ps>>18&0xFF)<<16|(ps>>9&0xFF)<<8|ps&0xFF;
		if(src==dst)
		{
			for(int kx=1;kx<bw;++kx)
			{
				int px=src[bw*ky+kx];
				auto p0=(const byte*)(src+bw*ky+kx-1);
				int R=px&0x1FF, G=px>>9&0x1FF, B=px>>18&0x1FF, neg;
				neg=R&1, R>>=1, R^=-neg, R+=neg, R+=p0[0];
				neg=G&1, G>>=1, G^=-neg, G+=neg, G+=p0[1];
				neg=B&1, B>>=1, B^=-neg, B+=neg, B+=p0[2];
				dst[bw*ky+kx]=B<<16|G<<8|R;
			}
		}
		else
		{
			for(int kx=1;kx<bw;++kx)
			{
				int px=src[bw*ky+kx], px0=src[bw*ky+kx-1];
				int R=px&0x1FF, G=px>>9&0x1FF, B=px>>18&0x1FF, neg;
				neg=R&1, R>>=1, R^=-neg, R+=neg, R+=px0&0x1FF;
				neg=G&1, G>>=1, G^=-neg, G+=neg, G+=px0>>9&0x1FF;
				neg=B&1, B>>=1, B^=-neg, B+=neg, B+=px0>>18&0x1FF;
				dst[bw*ky+kx]=B<<16|G<<8|R;
			}
		}
	}
}
void			RVL_rows_gray(const int *src, int bw, int bh, int depth0, int *dst)//can be the same buffer		8 -> 9bit
{
	const int mask=(1<<depth0)-1;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=bw-1;kx>=1;--kx)
		{
			auto p=(const byte*)(src+bw*ky+kx), p0=(const byte*)(src+bw*ky+kx-1);
			int val=(src[bw*ky+kx]&mask)-(src[bw*ky+kx-1]&mask), neg;
			neg=val<0, val^=-neg, val+=neg, val=val<<1|neg;
			dst[bw*ky+kx]=val;
		}
		dst[bw*ky]=src[bw*ky]&mask;
	}
}
void			invRVL_rows_gray(const int *src, int bw, int bh, int depth0, int *dst)//8 -> 9bit
{
	const int mask=(1<<depth0)-1, m2=(1<<(depth0+1))-1;
	for(int ky=0;ky<bh;++ky)
	{
		int ps=src[bw*ky];
		dst[bw*ky]=src[bw*ky]&mask;
		if(src==dst)
		{
			for(int kx=1;kx<bw;++kx)
			{
				int val=src[bw*ky+kx]&m2, neg;
				neg=val&1, val>>=1, val^=-neg, val+=neg, val+=src[bw*ky+kx-1];
				dst[bw*ky+kx]=val;
			}
		}
		else
		{
			for(int kx=1;kx<bw;++kx)
			{
				int val=src[bw*ky+kx]&m2, px0=src[bw*ky+kx-1], neg;
				neg=val&1, val>>=1, val^=-neg, val+=neg, val+=dst[bw*ky+kx-1];
				dst[bw*ky+kx]=val;
			}
		}
	}
}
void			apply_RCT27_29(int *dst, const int *src, int imsize)
{
	for(int k=0;k<imsize;++k)
	{
		int px=src[k];
		int R=px&0x1FF, G=px>>9&0x1FF, B=px>>18&0x1FF;
		int Y=(R+(G<<1)+B)>>2, Cb=(B-G)&0x3FF, Cr=(R-G)&0x3FF;
		dst[k]=Cr<<19|Cb<<9|Y;
	}
}
void			apply_invRCT29_27(int *dst, const int *src, int imsize)
{
	int mask=0xFFFFFC00;
	for(int k=0;k<imsize;++k)
	{
		int px=src[k];
		int Y=px&0x1FF, Cb=px>>9&0x3FF, Cr=px>>19&0x3FF;
		Cb|=mask&-(Cb>>9&1);
		Cr|=mask&-(Cr>>9&1);
		int G=Y-((Cb+Cr)>>2), R=Cr+G, B=Cb+G;
		dst[k]=B<<18|G<<9|R;
	}
}

void			apply_RCT24_26(int *dst, const int *src, int imsize)
{
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(src+k);
		byte R=p[0], G=p[1], B=p[2];
		int Y=(R+(G<<1)+B)>>2, Cb=(B-G)&0x1FF, Cr=(R-G)&0x1FF;
		dst[k]=Cr<<17|Cb<<8|Y;
	}
}
void			apply_invRCT26_24(int *dst, const int *src, int imsize)
{
	int mask=0xFFFFFE00;
	for(int k=0;k<imsize;++k)
	{
		int px=src[k];
		int Y=px&0xFF, Cb=px>>8&0x1FF, Cr=px>>17&0x1FF;
		Cb|=mask&-(Cb>>8&1);
		Cr|=mask&-(Cr>>8&1);
		byte G=Y-((Cb+Cr)>>2), R=Cr+G, B=Cb+G;
		dst[k]=B<<16|G<<8|R;
	}
}

void			apply_DWT24(int *buffer, int bw, int bh)
{
	//int imsize=bw*bh;
	//for(int k=0;k<imsize;++k)//preprocess
	//	buffer[k]^=0x808080;
	//{
	//	auto p=(byte*)(buffer+k);
	//	p[0]+=128;
	//	p[1]+=128;
	//	p[2]+=128;
	//}
	int maxdim=bw;
	if(maxdim<bh)
		maxdim=bh;
	auto temp=new int[maxdim];
	for(int ky=0;ky<bh;++ky)
	{
		auto row=buffer+bw*ky;
		for(int size=bw;size>2;)
		{
			int nodd=size>>1, neven=nodd+(size&1);
			for(int kx=0;kx+1<size;kx+=2)//predict: hi analysis {-0.5, 1, -0.5}
			{
				auto p0=(byte*)(row+kx), p1=(byte*)(row+kx+1), p2=(byte*)(row+kx+(2&-(kx+2<size)));
				p1[0]-=(p0[0]+p2[0])>>1;
				p1[1]-=(p0[1]+p2[1])>>1;
				p1[2]-=(p0[2]+p2[2])>>1;
			}
			for(int kx=0;kx+1<size;kx+=2)//update: lo analysis {0.25, 1, 0.25}
			{
				auto p0=(byte*)(row+kx+(kx==0)-(kx>0)), p1=(byte*)(row+kx), p2=(byte*)(row+kx+1);
				p1[0]+=(p0[0]+p2[0])>>2;
				p1[1]+=(p0[1]+p2[1])>>2;
				p1[2]+=(p0[2]+p2[2])>>2;
			}
			for(int kx=0, kx2=0;kx<neven;++kx, kx2+=2)//split to even & odd samples
				temp[kx]=row[kx2];
			for(int kx=neven, kx2=1;kx<size;++kx, kx2+=2)
				temp[kx]=row[kx2];
			memcpy(row, temp, size*sizeof(*row));
			size=neven;
		}
	}
	delete[] temp;
}
void			apply_invDWT24(int *buffer, int bw, int bh)
{
	int maxdim=bw;
	if(maxdim<bh)
		maxdim=bh;
	auto temp=new int[maxdim];
	std::vector<int> sizes;
	for(int size=bw;size>2;size=(size>>1)+(size&1))
		sizes.push_back(size);
	for(int ky=0;ky<bh;++ky)
	{
		auto row=buffer+bw*ky;
		for(int ks=sizes.size()-1;ks>=0;--ks)
		{
			int size=sizes[ks];
			int nodd=size>>1, neven=nodd+(size&1);
			for(int kx=0, kx2=0;kx<neven;++kx, kx2+=2)//mix even & odd samples
				temp[kx2]=row[kx];
			for(int kx=neven, kx2=1;kx<size;++kx, kx2+=2)
				temp[kx2]=row[kx];
			memcpy(row, temp, size*sizeof(*row));
			for(int kx=0;kx+1<size;kx+=2)//lo synthesis {-0.25, 1, -0.25}
			{
				auto p0=(byte*)(row+kx+(kx==0)-(kx>0)), p1=(byte*)(row+kx), p2=(byte*)(row+kx+1);
				p1[0]-=(p0[0]+p2[0])>>2;
				p1[1]-=(p0[1]+p2[1])>>2;
				p1[2]-=(p0[2]+p2[2])>>2;
			}
			for(int kx=0;kx+1<size;kx+=2)//hi synthesis {0.5, 1, 0.5}
			{
				auto p0=(byte*)(row+kx), p1=(byte*)(row+kx+1), p2=(byte*)(row+kx+(2&-(kx+2<size)));
				p1[0]+=(p0[0]+p2[0])>>1;
				p1[1]+=(p0[1]+p2[1])>>1;
				p1[2]+=(p0[2]+p2[2])>>1;
			}
			size=neven;
		}
	}
	delete[] temp;
	//int imsize=bw*bh;
	//for(int k=0;k<imsize;++k)//preprocess
	//	buffer[k]^=0x808080;
}


//irreversible color transform
void			apply_ICT_BT709(const int *buffer, int imsize, float *Y, float *Cb, float *Cr)
{
	const float gain=1.f/128;
	const float m[9]=
	{
		//gain,		0,		0,
		//0,		gain,	0,
		//0,		0,		gain,
		gain* 0.2126f,		gain* 0.7152f,		gain* 0.0722f,
		gain*-0.114572f,	gain*-0.385428f,	gain* 0.5f,
		gain* 0.5f,			gain* 0.454153f,	gain*-0.045847f,
	};
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		float R=float(p[0]-128), G=float(p[1]-128), B=float(p[2]-128);
		Y [k]=m[0]*R+m[1]*G+m[2]*B;
		Cb[k]=m[3]*R+m[4]*G+m[5]*B;
		Cr[k]=m[6]*R+m[7]*G+m[8]*B;
	}
}
void			apply_invICT_BT709(const float *Y, const float *Cb, const float *Cr, int imsize, int *buffer)
{
	const float gain=128;
	const float m[9]=
	{
		//gain,		0,		0,
		//0,		gain,	0,
		//0,		0,		gain,
		gain*-1.48852f,	gain*0.466159f,		gain*2.73974f,
		gain* 1.73974f,	gain*-0.325895f,	gain*-0.814412f,
		gain,			gain*1.8556f,		gain*-6.99456e-007f,
	};
	double Ramp=0, Gamp=0, Bamp=0;
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		float R=m[0]*Y[k]+m[1]*Cb[k]+m[2]*Cr[k]+gain;
		float G=m[3]*Y[k]+m[4]*Cb[k]+m[5]*Cr[k]+gain;
		float B=m[6]*Y[k]+m[7]*Cb[k]+m[8]*Cr[k]+gain;
		if(R<0)
			R=0;
		if(R>255)
			R=255;
		if(G<0)
			G=0;
		if(G>255)
			G=255;
		if(B<0)
			B=0;
		if(B>255)
			B=255;
		p[0]=(byte)R;
		p[1]=(byte)G;
		p[2]=(byte)B;
	}
}

//ICT subtracts mean
#if 0
void			apply_ICT_BT709(const int *buffer, int imsize, float *Y, float *Cb, float *Cr, float *means, float *amps)
{
	double Rmean=0, Gmean=0, Bmean=0;
	memset(amps, 0, 3*sizeof(float));
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		Rmean+=p[0];
		Gmean+=p[1];
		Bmean+=p[2];
		if(amps[0]<p[0])
			amps[0]=p[0];
		if(amps[1]<p[1])
			amps[1]=p[1];
		if(amps[2]<p[2])
			amps[2]=p[2];
	}
	double inv_imsize=1./imsize;
	means[0]=float(Rmean*inv_imsize);
	means[1]=float(Gmean*inv_imsize);
	means[2]=float(Bmean*inv_imsize);

	const float gain=2.f/255;
	const float m[9]=
	{
		gain* 0.2126,	gain* 0.7152,	gain* 0.0722,
		gain*-0.114572,	gain*-0.385428,	gain* 0.5,
		gain* 0.5,		gain* 0.454153,	gain*-0.045847,
	};
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		float R=p[0]-means[0], G=p[1]-means[1], B=p[2]-means[2];
		Y [k]=m[0]*R+m[1]*G+m[2]*B;//-1;
		Cb[k]=m[3]*R+m[4]*G+m[5]*B;//-1;
		Cr[k]=m[6]*R+m[7]*G+m[8]*B;//-1;
	}
}
void			apply_invICT_BT709(const float *means, const float *amps, float *Y, float *Cb, float *Cr, int imsize, int *buffer)//YCbCr are overwritten
{
	const float gain=255.f/2;
	const float m[9]=
	{
		gain*-1.48852,	gain*0.466159,	gain*2.73974,
		gain* 1.73974,	gain*-0.325895,	gain*-0.814412,
		gain,			gain*1.8556,	gain*-6.99456e-007,
	};
	double Ramp=0, Gamp=0, Bamp=0;
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		float R=m[0]*Y[k]+m[1]*Cb[k]+m[2]*Cr[k]+means[0];
		float G=m[3]*Y[k]+m[4]*Cb[k]+m[5]*Cr[k]+means[1];
		float B=m[6]*Y[k]+m[7]*Cb[k]+m[8]*Cr[k]+means[2];
		Y[k]=R, Cb[k]=G, Cr[k]=B;
		if(Ramp<R)
			Ramp=R;
		if(Gamp<G)
			Gamp=G;
		if(Bamp<B)
			Bamp=B;
	}
	Ramp=amps[0]/Ramp;
	Gamp=amps[1]/Gamp;
	Bamp=amps[2]/Bamp;
	printf("InvICT gains: %lf, %lf, %lf\n", Ramp, Gamp, Bamp);//
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
	//	auto R=int(Y[k]*Ramp), G=int(Cb[k]*Gamp), B=int(Cr[k]*Bamp);
		auto R=int(Y[k]), G=int(Cb[k]), B=int(Cr[k]);
		p[0]=byte(R), p[1]=byte(G), p[2]=byte(B);
	}
}
#endif

#if 0
void			apply_ICT_BT709_0(const int *buffer, int imsize, float *Y, float *Cb, float *Cr, float *means, float *amps)
{
	double Rmean=0, Gmean=0, Bmean=0;
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		Rmean+=p[0];
		Gmean+=p[1];
		Bmean+=p[2];
	}
	double inv_imsize=1./imsize;
	means[0]=float(Rmean*inv_imsize);
	means[1]=float(Gmean*inv_imsize);
	means[2]=float(Bmean*inv_imsize);

	const float gain=2.f/255;
	//const float gain=1;
	//const float Kr=0.2126, Kb=0.0722, Kg=1-Kr-Kb;
	const float m[9]=
	{
		gain* 0.2126,	gain* 0.7152,	gain* 0.0722,
		gain*-0.114572,	gain*-0.385428,	gain* 0.5,
		gain* 0.5,		gain* 0.454153,	gain*-0.045847,
		//Kr, Kg, Kb,
		//1, -Kg/(1-Kr), -Kb/(1-Kr),
		//-Kr/(1-Kb), -Kg/(1-Kb), 1,
	};
	memset(amps, 0, 3*sizeof(float));
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		float R=p[0]-means[0], G=p[1]-means[1], B=p[2]-means[2];
		Y [k]=m[0]*R+m[1]*G+m[2]*B;//-1;
		Cb[k]=m[3]*R+m[4]*G+m[5]*B;//-1;
		Cr[k]=m[6]*R+m[7]*G+m[8]*B;//-1;
		R=abs(R);
		G=abs(G);
		B=abs(B);
		if(amps[0]<R)
			amps[0]=R;
		if(amps[1]<G)
			amps[1]=G;
		if(amps[2]<B)
			amps[2]=B;
	}
}
void			apply_invICT_BT709_0(const float *means, const float *amps, float *Y, float *Cb, float *Cr, int imsize, int *buffer)//YCbCr are overwritten
{
	const float gain=255.f/2;
	//const float gain=1;
	const float m[9]=
	{
		gain*-1.48852,	gain*0.466159,	gain*2.73974,
		gain* 1.73974,	gain*-0.325895,	gain*-0.814412,
		gain,			gain*1.8556,	gain*-6.99456e-007,
		//gain, 0,				gain*1.5748,
		//gain, gain*-0.187324,	gain*-0.468124,
		//gain, gain*1.8556,	0,
	};
	double Ramp=0, Gamp=0, Bamp=0;
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		float R=m[0]*Y[k]+m[1]*Cb[k]+m[2]*Cr[k];
		float G=m[3]*Y[k]+m[4]*Cb[k]+m[5]*Cr[k];
		float B=m[6]*Y[k]+m[7]*Cb[k]+m[8]*Cr[k];
		Y[k]=R, Cb[k]=G, Cr[k]=B;
		R=abs(R);
		G=abs(G);
		B=abs(B);
		if(Ramp<R)
			Ramp=R;
		if(Gamp<G)
			Gamp=G;
		if(Bamp<B)
			Bamp=B;
	}
	Ramp=amps[0]/Ramp;
	Gamp=amps[1]/Gamp;
	Bamp=amps[2]/Bamp;
	printf("InvICT gains: %lf, %lf, %lf\n", Ramp, Gamp, Bamp);//
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		auto R=int(means[0]+Y[k]*Ramp), G=int(means[1]+Cb[k]*Gamp), B=int(means[2]+Cr[k]*Bamp);
		p[0]=byte(R), p[1]=byte(G), p[2]=byte(B);
	}
}
#endif


//	#define		QUANTIZER_NONLINEAR

//quantizer
float			get_amplitude(const float *data, int imsize)
{
	float amplitude=0;
	for(int k=0;k<imsize;++k)
	{
		float val=abs(data[k]);
		if(amplitude<val)
			amplitude=val;
	}
	return amplitude;
}
void			DWT_quantize(const float *Y, const float *Cb, const float *Cr, int imsize, int *buffer, float *amps)
{
	amps[0]=get_amplitude(Y, imsize);
	amps[1]=get_amplitude(Cb, imsize);
	amps[2]=get_amplitude(Cr, imsize);
#ifdef QUANTIZER_NONLINEAR
	float invY=1/amps[0], invCb=1/amps[1], invCr=1/amps[2];
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		float Yk=Y[k]*invY, Cbk=Cb[k]*invCb, Crk=Cr[k]*invCr;
		p[0]=(int)(127*sqrt(Yk))<<1|(Yk<0);
		p[1]=(int)(127*sqrt(Cbk))<<1|(Cbk<0);
		p[2]=(int)(127*sqrt(Crk))<<1|(Crk<0);
	}
#else
	float f_Y=127/amps[0], f_Cb=127/amps[1], f_Cr=127/amps[2];
	for(int k=0;k<imsize;++k)
	{
		auto p=(byte*)(buffer+k);
		int Yk=int(Y[k]*f_Y), Cbk=int(Cb[k]*f_Cb), Crk=int(Cr[k]*f_Cr), neg;
		neg=Yk <0, Yk ^=-neg, Yk +=neg, p[0]=Yk <<1|neg;
		neg=Cbk<0, Cbk^=-neg, Cbk+=neg, p[1]=Cbk<<1|neg;
		neg=Crk<0, Crk^=-neg, Crk+=neg, p[2]=Crk<<1|neg;
	}
#endif
}
void			DWT_dequantize(const int *buffer, int imsize, const float *amps, float *Y, float *Cb, float *Cr)
{
#ifdef QUANTIZER_NONLINEAR
	float ampY=amps[0]/127, ampCb=amps[1]/127, ampCr=amps[2]/127;
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		int neg;
		neg=(p[0]&1)<<31, Y [k]=ampY *(p[0]>>1), Y [k]*=Y [k], (int&)Y [k]=(int&)Y [k]|neg;
		neg=(p[1]&1)<<31, Cb[k]=ampCb*(p[1]>>1), Cb[k]*=Cb[k], (int&)Cb[k]=(int&)Cb[k]|neg;
		neg=(p[2]&1)<<31, Cr[k]=ampCr*(p[2]>>1), Cr[k]*=Cr[k], (int&)Cr[k]=(int&)Cr[k]|neg;
	}
#else
	float f_Y=amps[0]/127, f_Cb=amps[1]/127, f_Cr=amps[2]/127;
	for(int k=0;k<imsize;++k)
	{
		auto p=(const byte*)(buffer+k);
		int Yk=p[0], Cbk=p[1], Crk=p[2], neg;
		neg=Yk &1, Yk >>=1, Yk ^=-neg, Yk +=neg, Y [k]=Yk *f_Y ;
		neg=Cbk&1, Cbk>>=1, Cbk^=-neg, Cbk+=neg, Cb[k]=Cbk*f_Cb;
		neg=Crk&1, Crk>>=1, Crk^=-neg, Crk+=neg, Cr[k]=Crk*f_Cr;
	}
#endif
}
void			DWT_quantize_ch(const float *data, int imsize, int channel, int *buffer, float *amps)
{
	amps[channel]=get_amplitude(data, imsize);
	float gain=127/amps[0];
	for(int k=0;k<imsize;++k)
	{
		int val=int(data[k]*gain), neg;
		neg=val<0, val^=-neg, val+=neg, buffer[k]|=(val<<1|neg)<<(channel<<3);
	}
}
void			DWT_dequantize_ch(const int *buffer, int imsize, int channel, const float *amps, float *data)
{
	float gain=amps[channel]/127;
	for(int k=0;k<imsize;++k)
	{
		int val=buffer[k]>>(channel<<3)&0xFF, neg;
		neg=val&1, val>>=1, val^=-neg, val+=neg, data[k]=val*gain;
	}
}

//old transforms
void			extract_channel(const int *buffer, int imsize, int channel, float *data)//[-1, 1[
{
	auto p=(const byte*)buffer+channel;
	const double inv128=1./128;
	for(int k=0;k<imsize;++k, p+=4)
		data[k]=float((*p-128)*inv128);
}
void			assign_channel(const float *data, int imsize, int channel, int *buffer)
{
	auto p=(byte*)buffer+channel;
	for(int k=0;k<imsize;++k, p+=4)
		*p=int(data[k]*128)+128;
}
void			quantize_DWT(const float *data, int imsize, int channel, float amplitude, int *buffer)
{
	auto p=(byte*)buffer+channel;
	float factor=127/amplitude;
	for(int k=0;k<imsize;++k, p+=4)
	{
		int coeff=int(data[k]*factor);
		int neg=coeff<0;
		coeff^=-neg, coeff+=neg;
		*p=coeff<<1|neg;
	}
}
void			extract_DWT(const int *buffer, int imsize, int channel, float amplitude, float *data)
{
	auto p=(const byte*)buffer+channel;
	float factor=amplitude/127;
	for(int k=0;k<imsize;++k, p+=4)
	{
		int coeff=*p;
		int neg=coeff&1;
		coeff>>=1, coeff^=-neg, coeff+=neg;
		data[k]=coeff*factor;
	}
}


//DWT2 via lifting scheme
void			apply_DWT_1D(float *data, int size, int stride, float *temp, double *filt, int filtsize, double *norms)
{
#ifdef DEBUG_DWT
	printf("DWT1D input:\n"), print_data(data, size, stride);
#endif
	int nodd=size>>1, neven=nodd+(size&1), s2=stride<<1, bsize=size*stride;
	for(int k=0;k<filtsize;k+=2)
	{
		double predict=filt[k], update=filt[k+1];
		int kx;
		//predict: odd (high pass) += coeff * even neighbors
		for(kx=0;kx+s2<bsize;kx+=s2)
			data[kx+stride]+=float(predict*(data[kx]+data[kx+s2]));
		if(kx+stride<bsize)
			data[kx+stride]+=float(predict*(data[kx]+data[kx]));
#ifdef DEBUG_DWT
	printf("DWT1D predict %d:\n", k), print_data(data, size, stride);
#endif
			
		//update: even (low pass) += coeff * odd neighbors
		data[0]+=float(update*(data[stride]+data[stride]));
		for(kx=s2;kx+stride<bsize;kx+=s2)
			data[kx]+=float(update*(data[kx-stride]+data[kx+stride]));
#ifdef DEBUG_DWT
	printf("DWT1D update %d:\n", k+1), print_data(data, size, stride);
#endif
	}
	for(int kx=0, kx2=0;kx<neven;++kx, kx2+=s2)//split even & odd samples
		temp[kx]=float(norms[0]*data[kx2]);
	//	temp[kx]=data[kx2];
	for(int kx=neven, kx2=stride;kx<size;++kx, kx2+=s2)
		temp[kx]=float(norms[1]*data[kx2]);
	//	temp[kx]=data[kx2];
	if(stride==1)
		memcpy(data, temp, size*sizeof(*data));
	else
	{
		for(int ks=0, kd=0;ks<size;kd+=stride, ++ks)
			data[kd]=temp[ks];
	}
#ifdef DEBUG_DWT
	printf("DWT1D even-odd permutation:\n"), print_data(data, size, stride);
#endif
}
void			apply_invDWT_1D(float *data, int size, int stride, float *temp, double *filt, int filtsize, double *norms)
{
#ifdef DEBUG_DWT
	printf("InvDWT1D input:\n"), print_data(data, size, stride);
#endif
	int nodd=size>>1, neven=nodd+(size&1), s2=stride<<1, bsize=size*stride;
	if(stride==1)
		memcpy(temp, data, size*sizeof(*data));
	else
	{
		for(int ks=0, kd=0;ks<size;kd+=stride, ++ks)
			temp[ks]=data[kd];
	}
	double evennorm=1/norms[0], oddnorm=1/norms[1];
	for(int kx=0, kx2=0;kx<neven;++kx, kx2+=s2)//split even & odd samples
		data[kx2]=float(evennorm*temp[kx]);
	//	data[kx2]=temp[kx];
	for(int kx=neven, kx2=stride;kx<size;++kx, kx2+=s2)
		data[kx2]=float(oddnorm*temp[kx]);
	//	data[kx2]=temp[kx];
#ifdef DEBUG_DWT
	printf("InvDWT1D even-odd permutation:\n"), print_data(data, size, stride);
#endif
	for(int k=filtsize-2;k>=0;k-=2)
	{
		double predict=filt[k], update=filt[k+1];
		int kx;
		//update
		data[0]-=float(update*(data[stride]+data[stride]));
		for(kx=s2;kx+stride<bsize;kx+=s2)
			data[kx]-=float(update*(data[kx-stride]+data[kx+stride]));
#ifdef DEBUG_DWT
	printf("InvDWT1D update %d:\n", k), print_data(data, size, stride);
#endif

		//predict
		for(kx=0;kx+s2<bsize;kx+=s2)
			data[kx+stride]-=float(predict*(data[kx]+data[kx+s2]));
		if(kx+stride<bsize)
			data[kx+stride]-=float(predict*(data[kx]+data[kx]));
#ifdef DEBUG_DWT
	printf("InvDWT1D predict %d:\n", k), print_data(data, size, stride);
#endif
	}
}
void			apply_DWT_2D(float *buffer, int bw, int bh, double *filt, int filtsize, double *norms)//lifting scheme, filtsize is even, norms={lonorm, hinorm}
{
	int maxdim=bw;
	if(maxdim<bh)
		maxdim=bh;
	auto temp=new float[maxdim];
	for(int xsize=bw, ysize=bh;xsize>2&&ysize>2;)
	{
		//int xodd=xsize>>1, xeven=xodd+(xsize&1);
		//int yodd=ysize>>1, yeven=yodd+(ysize&1);
		for(int ky=0;ky<ysize;++ky)
			apply_DWT_1D(buffer+bw*ky, xsize, 1, temp, filt, filtsize, norms);
	/*	{
			auto row=buffer+bw*ky;
			for(int k=0;k<filtsize;k+=2)
			{
				double predict=filt[k], update=filt[k+1];
				int kx=0;
				for(;kx+2<xsize;kx+=2)
					row[kx+1]+=float(predict*(row[kx]+row[kx+2]));
				if(kx+1<xsize)
					row[kx+1]+=float(predict*(row[kx]+row[kx]));
			
				row[0]+=float(update*(row[1]+row[1]));
				for(kx=2;kx+2<xsize;kx+=2)
					row[kx]+=float(update*(row[kx-1]+row[kx+1]));
			}
			for(int kx=0, kx2=0;kx<xeven;++kx, kx2+=2)//split to even & odd samples
				temp[kx]=row[kx2];
			for(int kx=xeven, kx2=1;kx<xsize;++kx, kx2+=2)
				temp[kx]=row[kx2];
			memcpy(row, temp, xsize*sizeof(*row));
		}//*/
		for(int kx=0;kx<xsize;++kx)
			apply_DWT_1D(buffer+kx, ysize, bw, temp, filt, filtsize, norms);

		xsize-=xsize>>1, ysize-=ysize>>1;
	}
	delete[] temp;
}
void			DWT_sizes(int size, std::vector<int> &sizes)
{
	for(int ks=size;ks>2;ks-=ks>>1)
		sizes.push_back(ks);
}
void			apply_invDWT_2D(float *buffer, int bw, int bh, double *filt, int filtsize, double *norms)//lifting scheme, filtsize is even, norms={lonorm, hinorm}
{
	int maxdim=bw;
	if(maxdim<bh)
		maxdim=bh;
	auto temp=new float[maxdim];
	std::vector<int> xsizes, ysizes;
	DWT_sizes(bw, xsizes);
	DWT_sizes(bh, ysizes);
	int ks=xsizes.size()-1;
	if(ks>(int)ysizes.size()-1)
		ks=ysizes.size()-1;
	for(;ks>=0;--ks)
	{
		for(int kx=0;kx<xsizes[ks];++kx)
			apply_invDWT_1D(buffer+kx, ysizes[ks], bw, temp, filt, filtsize, norms);
		for(int ky=0;ky<ysizes[ks];++ky)
			apply_invDWT_1D(buffer+bw*ky, xsizes[ks], 1, temp, filt, filtsize, norms);
	}
	delete[] temp;
}


#if 0
const int	depth=8,
			nlevels=1<<depth;
int			histogram[nlevels]={};
#endif
void			load_image(const char *filename, int *&buffer, int &iw, int &ih)
{
	int nch=0;
	byte *original_image=stbi_load(filename, &iw, &ih, &nch, 4);
	if(original_image)
		printf("Opened \'%s\'\n", filename);
	else
	{
		printf("Failed to open \'%s\'\n", filename);
#ifdef _MSC_VER
		scanf_s("%d", &nch);
#endif
		exit(1);
	}
	buffer=(int*)original_image;
}
void			save_image(const char *filename, const int *buffer, int iw, int ih)
{
	printf("About to save result to:\n\n\t\'%s\'\n\nEnter ZERO to save: ", filename);
//	printf("Enter ZERO to save result image: ");
	int x=0;
	scanf_s("%d", &x);
	if(x)
		printf("Didn't save.\n");
	else
		lodepng::encode(filename, (const byte*)buffer, iw, ih);
}
void			gen_filename(double compression_ratio=0)
{
	time_t t=time(nullptr);
#ifdef __linux__
	auto &now=*localtime(&t);
#else
	tm now={};
	localtime_s(&now, &t);
#endif
	if(compression_ratio)
		sprintf_s(g_buf, G_BUF_SIZE, "%04d%02d%02d_%02d%02d%02d_%g.PNG", 1900+now.tm_year, now.tm_mon+1, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec, compression_ratio);
	else
		sprintf_s(g_buf, G_BUF_SIZE, "%04d%02d%02d_%02d%02d%02d.PNG", 1900+now.tm_year, now.tm_mon+1, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec);
}
void			match_buffers(int *buffer, int *b2, int imsize, int depth)
{
	int nerrors=0, kp=0, kb=0;
	int depthmask=(1<<depth)-1;
	for(int k=0;k<imsize;++k)
	{
		if((b2[k]^buffer[k])&depthmask)
		{
			if(!nerrors)
			{
				kb=k;
				for(kp=0;kp<depth&&!((b2[k]^buffer[k])>>kp&1);++kp);
			}
			++nerrors;
			printf("%d: %0*X != %0*X\n", k, (depth+3)>>2, buffer[k]&depthmask, (depth+3)>>2, b2[k]&depthmask);
			//if(nerrors>=100)
				break;
		}
	}
	if(nerrors)
	{
		const int width=64;
		printf("Error in bit plane %d, pixel %d:\n", kp, kb);
		int start=kb-(width>>1);
		if(start<0)
			start=0;
		int end=start+width;
		if(end>imsize)
			end=imsize;
		for(int k=start;k<end;++k)
			printf("%d", k%10);
		printf("\n");
		print_bitplane(buffer+start, end-start, sizeof(*buffer), kp);//
		print_bitplane(b2+start, end-start, sizeof(*buffer), kp);//
		//print_image(buffer, iw, ih, 0, iw, 0, ih);
		//print_image(b2, iw, ih, 0, iw, 0, ih);
	}
}
enum			Transform
{
	TR_97DWT,
	TR_RAW,
	TR_COLOR_97DWT,
	TR_RVL,
};
enum			Coder
{
	CODER_RANS,
	CODER_RANS32,
	CODER_HUFFMAN,
	CODER_ABAC2,
	CODER_ZPAQ0_AC,
	CODER_LZ77_RANS,
	CODER_LZ77_RANS32,
	CODER_LZ77_HUFFMAN,
	CODER_LZ77_ABAC2,
	CODER_LZ77_ZPAQ0_AC,
};
int				main(int argc, char **argv)
{
	set_console_buffer_size(120, 4000);
	if(argc<2)
	{
		printf("Pass filename as command argument\n");
		return 1;
	}

	//interactive
#if 1
	int *buffer=nullptr, iw=0, ih=0, imsize=0, nch=0, depths[4]={}, total_depth=0;
	int transform=0;

	int len=strlen(argv[1]);
	if(	!strcmp(argv[1]+len-4, ".txt")||!strcmp(argv[1]+len-4, ".TXT")||
		!strcmp(argv[1]+len-4, ".cpp")||!strcmp(argv[1]+len-4, ".CPP")||
		!strcmp(argv[1]+len-2, ".c")||!strcmp(argv[1]+len-2, ".C"))
	{
		std::string text;
		if(open_text(argv[1], text))
			printf("Opened \'%s\'\n", argv[1]);
		else
		{
			printf("Failed to open \'%s\'\n", argv[1]);
			return 1;
		}
		iw=text.size(), ih=1, imsize=iw*ih, nch=1, depths[1]=8, total_depth=depths[1];
		buffer=new int[imsize];
		for(int k=0;k<imsize;++k)
			buffer[k]=text[k];
		transform=TR_RAW;
	}
	else
	{
		auto original_image=stbi_load(argv[1], &iw, &ih, &nch, 4);
		if(original_image)
			printf("Opened \'%s\'\n", argv[1]);
		else
		{
			printf("Failed to open \'%s\'\n", argv[1]);
			return 1;
		}
		printf("Image has %d channels\n", nch);
		if(nch==4)
			printf(
				"\n"
				"Select channels to take.\n"
				"    0   Only the red channel\n"
				"   [1]  RGB\n"
				"    2   RGBA\n"
				);
		else
			printf(
				"\n"
				"Select channels to take.\n"
				"    0   Only the red channel\n"
				"   [1]  RGB\n"
				);
		int mask=_getch();
		//int mask=0;
		//scanf_s("%d", &mask);
		switch(mask)
		{
		case '0':
			printf("%c: Only the red channel\n", (char)mask);
			mask=0xFF;
			depths[0]=8, total_depth=8;
			break;
		case '2':
			if(nch==4)
			{
				printf("%c: RGBA\n", (char)mask);
				mask=0xFFFFFFFF;
				depths[0]=depths[1]=depths[2]=depths[3]=8, total_depth=32;
				break;
			}
			//no break
		default:
		case '1':
			printf("%c: RGB\n", (char)mask);
			mask=0xFFFFFF;
			depths[0]=depths[1]=depths[2]=8, total_depth=24;
			break;
		}
	/*	printf("Enter a mask to select channels (ABGR, 1 or 0 per digit, eg: BGR=111, R=1, ...etc): ");
		int mask=0;
		for(;;)
		{
			scanf_s("%d", &mask);
			if(mask)
				break;
			printf("No channels selected. Enter a channel mask: ");
		}
		for(int k=0;k<nch;++k, mask/=10)
			if(mask%10)
				depths[k]=8, total_depth+=depths[k];
		mask=((1<<depths[3])-1)<<24|((1<<depths[2])-1)<<16|((1<<depths[1])-1)<<8|((1<<depths[0])-1);//*/

		imsize=iw*ih;
		buffer=new int[imsize];
		auto b2=(int*)original_image;
		for(int k=0;k<imsize;++k)
			buffer[k]=b2[k]&mask;
		STBI_FREE(original_image);
		printf(
			"\n"
			"Select transform:\n"
			"   [0]  9-7 DWT\n"
			"    1   Raw\n"
			"    2   Color transform + 9-7 DWT\n"
			"    3   RVL (+1 depth bit per channel)\n"
			);
		transform=_getch();
		switch(transform)
		{
		default:
		case '0':printf("%c: 9-7 DWT\n", (char)transform);	transform=TR_97DWT;break;
		case '1':printf("%c: RAW\n", (char)transform);		transform=TR_RAW;break;
		case '2':
			if(depths[0]||depths[1]&&depths[2])
			{
				printf("%c: Color transform + 9-7 DWT\n", (char)transform);
				transform=TR_COLOR_97DWT;
			}
			else
			{
				printf("Can't apply color transform because not all color channels were selected\n");
				transform=TR_RAW;
			}
			break;
		case '3':
			if(depths[3])
			{
				printf("%c: Can't select RVL because alpha is selected. Deferring to RAW.\n", (char)transform);
				transform=TR_RAW;
			}
			else
			{
				printf("%c: RVL\n", (char)transform);
				transform=TR_RVL;
			}
			break;
		}
	}
	printf(
		"\n"
		"Select coder:\n"
		"    Raw   LZ77\n"
		"   [0]    A     rANS\n"
		"    1     B     rANS32\n"
		"    2     C     Huffman\n"
		"    3     D     ABAC2\n"
		"    4     E     ZPAQ0 AC\n"
		);
	int coder=_getch();
	switch(coder)
	{
	default:
	case '0':printf("%c: CODER_RANS\n",				(char)coder);	coder=CODER_RANS;break;
	case '1':printf("%c: CODER_RANS32\n",			(char)coder);	coder=CODER_RANS32;break;
	case '2':printf("%c: CODER_HUFFMAN\n",			(char)coder);	coder=CODER_HUFFMAN;break;
	case '3':printf("%c: CODER_ABAC2\n",			(char)coder);	coder=CODER_ABAC2;break;
	case '4':printf("%c: CODER_ZPAQ0_AC\n",			(char)coder);	coder=CODER_ZPAQ0_AC;break;
	case 'A':printf("%c: CODER_LZ77_RANS\n",		(char)coder);	coder=CODER_LZ77_RANS;break;
	case 'B':printf("%c: CODER_LZ77_RANS32\n",		(char)coder);	coder=CODER_LZ77_RANS32;break;
	case 'C':printf("%c: CODER_LZ77_HUFFMAN\n",		(char)coder);	coder=CODER_LZ77_HUFFMAN;break;
	case 'D':printf("%c: CODER_LZ77_ABAC2\n",		(char)coder);	coder=CODER_LZ77_ABAC2;break;
	case 'E':printf("%c: CODER_LZ77_ZPAQ0_AC\n",	(char)coder);	coder=CODER_LZ77_ZPAQ0_AC;break;
	}
	printf("\n");
	
	double filt[]=
	{
		// 1,
		// 1,
		 1.1496043988602962433651033986614476,
		 0.86986445162473959153241758552174375,

		-1.5861343420594238292020515785937243,
		-0.052980118573376671019344897514278511,
		 0.882911075528503100647487289764588328,
		 0.443506852044983007635941975182636061,
	};
	float *f_channels[4]={};
	auto b2=new int[imsize];
	memset(b2, 0, imsize*sizeof(*b2));
	float amps[4]={};
	float gain=1.f/128;
	switch(transform)//decorrelating transform		buffer -> b2
	{
	case TR_97DWT:
		{
			auto t1=__rdtsc();
			for(int kc=0;kc<4;++kc)
			{
				if(depths[kc])
				{
					f_channels[kc]=new float[imsize];
					for(int k=0;k<imsize;++k)
						f_channels[kc][k]=gain*((buffer[k]>>(kc<<3)&0xFF)-128);
					apply_DWT_2D(f_channels[kc], iw, ih, filt+2, 4, filt);
					DWT_quantize_ch(f_channels[kc], imsize, kc, b2, amps);
				}
			}
			auto t2=__rdtsc();
			printf("DWT + Quantization: %lld cycles, %lf CPB\n", t2-t1, (double)(t2-t1)/(imsize*total_depth>>3));
		}
		break;
	case TR_RAW:
		memcpy(b2, buffer, imsize*sizeof(*buffer));
		break;
	case TR_COLOR_97DWT:
		f_channels[0]=new float[imsize];
		f_channels[1]=new float[imsize];
		f_channels[2]=new float[imsize];
		if(depths[3])
		{
			f_channels[3]=new float[imsize];
			for(int k=0;k<imsize;++k)
				f_channels[3][k]=gain*((buffer[k]>>24&0xFF)-128);
			apply_DWT_2D(f_channels[3], iw, ih, filt+2, 4, filt);
		}
		apply_ICT_BT709(buffer, imsize, f_channels[0], f_channels[1], f_channels[2]);
		apply_DWT_2D(f_channels[0], iw, ih, filt+2, 4, filt);
		apply_DWT_2D(f_channels[1], iw, ih, filt+2, 4, filt);
		apply_DWT_2D(f_channels[2], iw, ih, filt+2, 4, filt);
		DWT_quantize(f_channels[0], f_channels[1], f_channels[2], imsize, b2, amps);
		if(depths[3])
			DWT_quantize_ch(f_channels[3], imsize, 3, b2, amps);
		break;
	case TR_RVL:
		if(total_depth==8)
		{
			RVL_rows_gray(buffer, iw, ih, 8, b2);
			depths[0]=9, total_depth=9;
		}
		else if(total_depth==24)
		{
			RVL_rows(buffer, iw, ih, b2);
			depths[0]=depths[1]=depths[2]=9, total_depth=27;
		}
		break;
	}
	switch(coder)
	{
	case CODER_LZ77_RANS:
	case CODER_LZ77_HUFFMAN:
	case CODER_LZ77_ABAC2:
	case CODER_LZ77_ZPAQ0_AC:
		printf("LZ77. TODO.\n");
		break;
	}
	switch(coder)//entropy coding		b2 -> buffer
	{
	case CODER_RANS:
	case CODER_LZ77_RANS:
		{
			std::vector<unsigned short> data;
			if(total_depth==8)
			{
				auto freqs=rans_gray_start<int, 8>(b2, imsize, true);
				rans_gray_encode<int, 8>(b2, imsize, freqs, data, 2);
				rans_gray_decode<int, 8>(data.data(), data.size(), imsize, freqs, buffer, true);
				rans_gray_finish(freqs);
			}
			if(total_depth==9)
			{
				auto freqs=rans_gray_start<int, 9>(b2, imsize, true);
				rans_gray_encode<int, 9>(b2, imsize, freqs, data, 2);
				rans_gray_decode<int, 9>(data.data(), data.size(), imsize, freqs, buffer, true);
				rans_gray_finish(freqs);
			}
			else if(total_depth==24)
			{
				auto freqs=rans_rgb_start<int, 8, 8, 8>(b2, imsize, true);
				rans_rgb_encode<int, 8, 8, 8>(b2, imsize, freqs, data, 2);
				rans_rgb_decode<int, 8, 8, 8>(data.data(), data.size(), imsize, freqs, buffer, true);
				rans_rgb_finish(freqs);
			}
			else if(total_depth==27)
			{
				auto freqs=rans_rgb_start<int, 9, 9, 9>(b2, imsize, true);
				rans_rgb_encode<int, 9, 9, 9>(b2, imsize, freqs, data, 2);
				rans_rgb_decode<int, 9, 9, 9>(data.data(), data.size(), imsize, freqs, buffer, true);
				rans_rgb_finish(freqs);
			}
		}
		break;
	case CODER_RANS32:
	case CODER_LZ77_RANS32:
		{
			std::vector<unsigned> data;
			unsigned short freqs[1024];
			int s2=imsize;
			rans3_encode(b2, s2, data, freqs, RANS_SERIAL, 2);
			rans3_decode(data.data(), data.size(), freqs, buffer, s2, RANS_SERIAL, true);
			//std::string data;
			//unsigned char freqs[128]={};
			//int s2=imsize;
			//rans2_encode(b2, s2, data, freqs, RANS_SERIAL, 2);
			//rans2_decode(data.data(), data.size(), freqs, buffer, s2, RANS_SSE2, true);
		}
		break;
	case CODER_HUFFMAN:
	case CODER_LZ77_HUFFMAN:
		{
			printf("Huffman expects 16bit buffer, not 32bit. TODO.\n");
			memcpy(buffer, b2, imsize*sizeof(*buffer));//
		}
		break;
	case CODER_ABAC2:
	case CODER_LZ77_ABAC2:
		{
			std::string data;
			auto sizes=new int[total_depth];
			auto conf=new int[total_depth];
			abac2_encode(b2, imsize, total_depth, sizeof(*b2), data, sizes, conf, true);
			abac2_decode(data.data(), sizes, conf, buffer, imsize, total_depth, sizeof(*buffer), true);
			delete[] sizes, conf;
		}
		break;
	case CODER_ZPAQ0_AC:
	case CODER_LZ77_ZPAQ0_AC:
		{
			printf("ABAC3 expects 16bit buffer, not 32bit. TODO.\n");
			memcpy(buffer, b2, imsize*sizeof(*buffer));//
		}
		break;
	}
/*	switch(coder)
	{
	case CODER_LZ77_RANS:
	case CODER_LZ77_HUFFMAN:
	case CODER_LZ77_ABAC2:
	case CODER_LZ77_ZPAQ0_AC:
		{
			printf("LZ77. TODO.\n");
		//	memcpy(buffer, b2, imsize*sizeof(*buffer));//
		}
		break;
	}//*/
	match_buffers(b2, buffer, imsize, total_depth);
	switch(transform)//inverse transform	buffer -> b2
	{
	case TR_97DWT:
		{
			auto t1=__rdtsc();
			for(int kc=0;kc<4;++kc)
			{
				if(depths[kc])
				{
					DWT_dequantize_ch(buffer, imsize, kc, amps, f_channels[kc]);
					apply_invDWT_2D(f_channels[kc], iw, ih, filt+2, 4, filt);
					gain=127/amps[kc];
					for(int k=0;k<imsize;++k)
						b2[k]|=((128+(int)(gain*f_channels[kc][k]))&0xFF)<<(kc<<3);
				}
			}
			auto t2=__rdtsc();
			printf("Dequantization + InvDWT: %lld cycles, %lf CPB\n", t2-t1, (double)(t2-t1)/(imsize*total_depth>>3));
		}
		break;
	case TR_RAW:
		memcpy(b2, buffer, imsize*sizeof(*buffer));
		break;
	case TR_COLOR_97DWT:
		DWT_dequantize(buffer, imsize, amps, f_channels[0], f_channels[1], f_channels[2]);
		apply_invDWT_2D(f_channels[0], iw, ih, filt+2, 4, filt);
		apply_invDWT_2D(f_channels[1], iw, ih, filt+2, 4, filt);
		apply_invDWT_2D(f_channels[2], iw, ih, filt+2, 4, filt);
		apply_invICT_BT709(f_channels[0], f_channels[1], f_channels[2], imsize, b2);
		delete[] f_channels[0], f_channels[1], f_channels[2];
		if(depths[3])
		{
			DWT_dequantize_ch(buffer, imsize, 3, amps, f_channels[3]);
			apply_invDWT_2D(f_channels[3], iw, ih, filt+2, 4, filt);
			for(int k=0;k<imsize;++k)
				b2[k]|=((128+(int)(gain*f_channels[3][k]))&0xFF)<<24;
			delete[] f_channels[3];
		}
		break;
	case TR_RVL:
		if(total_depth==9)
		{
			invRVL_rows_gray(buffer, iw, ih, 8, b2);
			depths[0]=8, total_depth=8;
		}
		else if(total_depth==27)
		{
			invRVL_rows(buffer, iw, ih, b2);
			depths[0]=depths[1]=depths[2]=8, total_depth=24;
		}
		break;
	}
	delete[] buffer, b2;
#endif

	//LZ77
#if 0
	int iw=0, ih=0, nch=0;
	//byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	byte *original_image=stbi_load("20211129 1 confidence.PNG", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Couldn't open image\n");
		return 1;
	}
	auto buffer=(char*)original_image;
	int imsize=iw*ih;//*/

/*	std::string str;
//	open_text("D:/C/Compression Sandbox/Compression Sandbox/ac.cpp", str);
	open_text("D:/C/Compression Sandbox/g2.cpp", str);
	auto buffer=str.c_str();
	int imsize=str.size();//*/

/*	const char buffer[]="01234555555555555555555555555555555555555556789";
	//const char buffer[]="Hello World! Sample Text. What is going on??!";
	//const char buffer[]="0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_";
	printf("Original:\n%s\n", buffer);
	int imsize=sizeof(buffer);//*/

	//const int depth=8;
	char *b2=new char[imsize];
	memset(b2, 0, imsize);

	std::string data;
	lz77_encode(buffer, imsize, data, true);
	lz77_decode(data.data(), data.size(), imsize, b2, true);
	bool error=false;
	for(int k=0;k<imsize;++k)
	{
		if(b2[k]!=buffer[k])
		{
			error=true;
			printf("Error %d: %d -> %d, %c -> %c\n", k, b2[k]&0xFF, buffer[k]&0xFF, b2[k], buffer[k]);
			int start=k-100;
			if(start<0)
				start=0;
			int end=start+200;
			if(end>imsize)
				end=imsize;
			printf("Before:\n%.*s\nAfter:\n%.*s\n", end-start, buffer+start, end-start, b2+start);
			break;
		}
	}
	if(!error)
	{
		int preview=imsize;
		if(preview>200)
			preview=200;
		printf("Decoded:\n%.*s\n", preview, b2);
	}
	
	//STBI_FREE(original_image);
	delete[] b2;
#endif

	//Huffman - text
#if 0
	std::string text;
	open_text("abac.cpp", text);
	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];
	
	vector_bool bits;
	huff_encode(buffer, imsize, 8, bits, true);
/*	auto t1=__rdtsc();
	calculate_histogram(buffer, imsize, histogram, nlevels);
	auto t_hist=__rdtsc();
	build_tree(histogram, nlevels);
	auto t_tree=__rdtsc();
	std::vector<vector_bool> alphabet;
	make_alphabet(alphabet);
	auto t_alpha=__rdtsc();
	int size_estimated=huff_rough_size_estimate(alphabet, imsize);
	vector_bool bits;
	bits.data.reserve(size_estimated);
	for(int k=0;k<imsize;++k)
		bits.push_back(alphabet[buffer[k]]);
	bits.clear_tail();
	auto t2=__rdtsc();

	int compressed_bytesize=bits.data.size()*sizeof(int);
	printf("Huffman:  %lld cycles\n", t2-t1);
	printf("Size: %d -> %d bytes,  ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
	printf("Histogram:\t%lld\n", t_hist-t1);
	printf("Build tree:\t%lld\n", t_tree-t_hist);
	printf("Alphabet:\t%lld\n", t_alpha-t_tree);
	printf("Pack:\t%lld\n", t2-t_alpha);//*/
#endif

	//Huffman - image
#if 0
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Failed to open image\n");
		return 1;
	}
	auto image=(int*)original_image;
	//unsigned iw=0, ih=0;
	//std::vector<unsigned char> vimage;
	//int error=lodepng::decode(vimage, iw, ih, "example.png");
	//if(error)
	//{
	//	printf("%s\n", lodepng_error_text(error));
	//	return 1;
	//}
	//auto image=(int*)vimage.data();

	int imsize=iw*ih;
//	int imsize=800;
	short *buffer=new short[imsize];
	const int depth=8;
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;
	short *b2=new short[imsize];
	STBI_FREE(original_image);


	int histogram[256]={};
	vector_bool data;
	huff_encode(buffer, imsize, 8, histogram, data, true);

	huff_decode(data.data.data(), data.bitSize, imsize, 8, histogram, b2, true);


	delete[] buffer;
	delete[] b2;
#endif
	
	//AC - quantized float
#if 0
#if 1
#define FREE_ORIGINAL_IM
	int iw=0, ih=0;
	int *buffer=nullptr;
	load_image(argv[1], buffer, iw, ih);
#else
	int iw=8, ih=8;
	int buffer[]=
	{
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255,
	};

#endif
	int imsize=iw*ih;
	auto b2=new int[imsize];
	memset(b2, 0xFF, imsize*sizeof(*b2));

	double filt[]=
	{
		 1,
		 1,
		 //1.1496043988602962433651033986614476,
		 //0.86986445162473959153241758552174375,

		-1.5861343420594238292020515785937243,
		-0.052980118573376671019344897514278511,
		 0.882911075528503100647487289764588328,
		 0.443506852044983007635941975182636061,
	};
	
/*	const float gain=1;
	const float m1[9]=
	{
		gain* 0.2126,	gain* 0.7152,	gain* 0.0722,
		gain*-0.114572,	gain*-0.385428,	gain* 0.5,
		gain* 0.5,		gain* 0.454153,	gain*-0.045847,
	};
	const float m2[9]=
	{
		gain*-1.48852,	gain*0.466159,	gain*2.73974,
		gain* 1.73974,	gain*-0.325895,	gain*-0.814412,
		gain,			gain*1.8556,	gain*-6.99456e-007,
	};
	float m3[9]={};
	for(int k=0;k<3;++k)
		for(int k2=0;k2<3;++k2)
			for(int k3=0;k3<3;++k3)
				m3[k*3+k2]+=m1[k*3+k3]*m2[k3*3+k2];
	print_fdata(m3, 3, 3, 0, 3, 0, 3);
	_getch();
	exit(0);//*/
	auto Y=new float[imsize], Cb=new float[imsize], Cr=new float[imsize];
	apply_ICT_BT709(buffer, imsize, Y, Cb, Cr);
	apply_DWT_2D(Y, iw, ih, filt+2, 4, filt);
	apply_DWT_2D(Cb, iw, ih, filt+2, 4, filt);
	apply_DWT_2D(Cr, iw, ih, filt+2, 4, filt);

	float amps[3]={};
	DWT_quantize(Y, Cb, Cr, imsize, b2, amps);
	
	const int depth=24;
#if 1
#define USE_RANS
	std::vector<unsigned short> data;
	auto freqs=rans_rgb888_start(b2, imsize, true);
	rans_rgb888_encode(b2, imsize, freqs, data, 2);
	rans_rgb888_decode(data.data(), data.size(), imsize, freqs, buffer, true);
	//memcpy(buffer, b2, imsize*sizeof(*buffer));//
	rans_rgb888_finish(freqs);
#else
	std::string cdata;
	int sizes[depth]={};
	int conf[depth]={};
	abac2_encode(b2, imsize, depth, sizeof(*b2), cdata, sizes, conf, true);
	abac2_decode(cdata.data(), sizes, conf, buffer, imsize, depth, sizeof(*b2), true);
#endif
	match_buffers(b2, buffer, imsize, depth);
	
	gen_filename();//TODO: mention amplitudes in filename
	save_image(g_buf, b2, iw, ih);

	DWT_dequantize(b2, imsize, amps, Y, Cb, Cr);
	apply_invDWT_2D(Y, iw, ih, filt+2, 4, filt);
	apply_invDWT_2D(Cb, iw, ih, filt+2, 4, filt);
	apply_invDWT_2D(Cr, iw, ih, filt+2, 4, filt);//*/
	apply_invICT_BT709(Y, Cb, Cr, imsize, b2);
	
#ifdef USE_RANS
	gen_filename((double)imsize*3/(data.size()*sizeof(short)));
#else
	gen_filename((double)imsize*3/cdata.size());
#endif
	//gen_filename();
	save_image(g_buf, b2, iw, ih);


//	auto data=new float[imsize];
	
/*	for(int ky=0;ky<ih;++ky)
		for(int kx=0;kx<iw;++kx)
			data[iw*ky+kx]=(kx^ky)&1;
	int x1=0, x2=8, y1=0, y2=8;
	print_fdata(data, iw, ih, x1, x2, y1, y2);
	apply_DWT_2D	(data, iw, ih, filt, 4, filt+4);
	print_fdata(data, iw, ih, x1, x2, y1, y2);
	apply_invDWT_2D	(data, iw, ih, filt, 4, filt+4);
	print_fdata(data, iw, ih, x1, x2, y1, y2);//*/

/*	float temp[8]={};
	print_data(data, 8, 1);
	apply_DWT_1D(data+1, 4, 2, temp, filt, 4, filt+4);
	print_data(data, 8, 1);
	apply_invDWT_1D(data+1, 4, 2, temp, filt, 4, filt+4);
	print_data(data, 8, 1);//*/

/*	for(;;)
	{
		for(int k=0;k<imsize;++k)//
			buffer[k]=rand()&0xFF;
		extract_channel(buffer, imsize, 0, data);
		float temp[8]={};
		apply_DWT_1D(data, 8, 1, temp, filt, 4, filt+4);
		apply_invDWT_1D(data, 8, 1, temp, filt, 4, filt+4);
		assign_channel(data, imsize, 0, b2);

		print_image(buffer, iw, ih, 0, iw, 0, ih);
		print_image(b2, iw, ih, 0, iw, 0, ih);
		print_diff(buffer, b2, iw, ih, 0, 8, 0, 8, 0x000000FF);
		_getch();
	}//*/

//	memcpy(b2, buffer, imsize*sizeof(*b2));

/*	float amplitude[3]={};
	int x1=0, x2=8, y1=0, y2=8;
	for(int kc=0;kc<3;++kc)//encode channels
	{
		printf("Encoding channel %d...\n", kc);
		extract_channel(buffer, imsize, kc, data);
		//for(int ky=0;ky<ih;++ky)
		//	for(int kx=0;kx<iw;++kx)
		//		data[iw*ky+kx]=(kx^ky)&1;
		//	printf("Before DWT:\n"), print_fdata(data, iw, ih, x1, x2, y1, y2);
		apply_DWT_2D(data, iw, ih, filt, 4, filt+4);
		amplitude[kc]=get_amplitude(data, imsize);
		printf("Channel %d amplitude: %f\n", kc, amplitude);
		quantize_DWT(data, imsize, kc, amplitude[kc], b2);

	//	extract_DWT(b2, imsize, kc, amplitude[kc], data);//
	//	apply_invDWT_2D(data, iw, ih, filt, 4, filt+4);//
	//	assign_channel(data, imsize, kc, b2);//

		//	printf("DWT:\n"), print_fdata(data, iw, ih, x1, x2, y1, y2);
		//apply_invDWT_2D(data, iw, ih, filt, 4, filt+4);
		//	printf("Inverse DWT:\n"), print_fdata(data, iw, ih, x1, x2, y1, y2), printf("\n");
		//assign_channel(data, imsize, kc, b2);
	}

	const int depth=24;
	std::string cdata;
	int sizes[depth]={};
	int conf[depth]={};
	abac2_encode(b2, imsize, depth, sizeof(*b2), cdata, sizes, conf, true);
	abac2_decode(cdata.data(), sizes, conf, buffer, imsize, depth, sizeof(*b2), true);
	match_buffers(b2, buffer, imsize, depth);

	gen_filename();
	save_image(g_buf, b2, iw, ih);

	for(int kc=0;kc<3;++kc)//decode channels
	{
		printf("Decoding channel %d...\n", kc);
		extract_DWT(b2, imsize, kc, amplitude[kc], data);//
		apply_invDWT_2D(data, iw, ih, filt, 4, filt+4);//
		assign_channel(data, imsize, kc, b2);//
	}
//	float amplitude[]={101.232773, 82.474693, 165.159637};//kodim23.png
//	for(int kc=0;kc<3;++kc)//decode channels
//	{
//		printf("Processing channel %d...\n", kc);
//		extract_DWT(buffer, imsize, kc, amplitude[kc], data);
//		apply_invDWT_2D(data, iw, ih, filt, 4, filt+4);
//		assign_channel(data, imsize, kc, b2);
//	}

	gen_filename((double)imsize*3/cdata.size());
	save_image(g_buf, b2, iw, ih);//*/

	delete[] Y, Cb, Cr;
	//delete[] data;
	delete[] b2;
#ifdef FREE_ORIGINAL_IM
	stbi_image_free(buffer);
#endif
#endif

	//AC - color image
#if 0
#if 1
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load(argv[1], &iw, &ih, &nch, 4);
#define FREE_ORIGINAL_IM
	if(original_image)
		printf("Opened \'%s\'\n", argv[1]);
	else
	{
		printf("Failed to open \'%s\'\n", argv[1]);
		return 1;
	}
	auto buffer=(int*)original_image;
#else
	int iw=8, ih=8;
	int buffer[]=
	{
		1, 2, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
	};
#endif

	int imsize=iw*ih;
	auto b2=new int[imsize];
	const int depth0=24;
	//const int depth=29;
	const int depth=26;
	
	//for(int k=0;k<imsize;++k)//
	//	buffer[k]=rand();
	memcpy(b2, buffer, imsize*sizeof(*buffer));
	//	printf("Before:\n"), print_image(b2, iw, ih, iw>>1, (iw>>1)+8, ih>>1, (ih>>1)+8);
	//apply_DWT24(b2, iw, ih);
	//	printf("DWT:\n"), print_image(b2, iw, ih, iw>>1, (iw>>1)+8, ih>>1, (ih>>1)+8);

	//RVL_rows(buffer, iw, ih, b2);
	//apply_RCT27_29(b2, b2, imsize);
	//apply_RCT24_26(b2, b2, imsize);
		//printf("RCT26:\n"), print_image(b2, iw, ih, 0, iw, 0, ih);
	
#if 0
#define	USE_RANS
	std::vector<unsigned short> data;
	int subtract=0;//742	801
	if(subtract)
		printf("imsize: %d\nleaving last %d pixels\n", imsize, subtract);

	auto freqs=rans_rgb888_start(buffer, imsize-subtract, true);
	rans_rgb888_encode(buffer, imsize-subtract, freqs, data, 2);
	rans_rgb888_decode(data.data(), data.size(), imsize-subtract, freqs, b2, true);

	rans_rgb888_finish(freqs);
#else
	std::string data;
	int sizes[depth]={};
	int conf[depth]={};
	abac2_encode(b2, imsize, depth, sizeof(*b2), data, sizes, conf, true);
	abac2_decode(data.data(), sizes, conf, b2, imsize, depth, sizeof(*b2), true);
#endif

	//apply_invRCT29_27(b2, b2, imsize);
	//invRVL_rows(b2, iw, ih, b2);
	//apply_invRCT26_24(b2, b2, imsize);//*/
		//printf("IRCT26:\n"), print_image(b2, iw, ih, 0, iw, 0, ih);

	//apply_invDWT24(b2, iw, ih);
		//printf("IDWT:\n"), print_image(b2, iw, ih, 0, iw, 0, ih);
	
	int nerrors=0, kp=0, kb=0;
	int depthmask=(1<<depth0)-1;
	//printf("Depthmask: %04X\n", depthmask);
#ifdef USE_RANS
	for(int k=imsize-1;k>=0;--k)
#else
	for(int k=0;k<imsize;++k)
#endif
	{
		if((b2[k]^buffer[k])&depthmask)
		{
			if(!nerrors)
			{
				kb=k;
				for(kp=0;kp<depth0&&!((b2[k]^buffer[k])>>kp&1);++kp);
			}
			++nerrors;
			printf("%d: %0*X != %0*X\n", k, (depth0+3)>>2, buffer[k]&depthmask, (depth0+3)>>2, b2[k]&depthmask);
			//if(nerrors>=100)
				break;
		}
	}
	if(nerrors)
	{
		const int width=64;
		printf("Error in bit plane %d, pixel %d:\n", kp, kb);
		int start=kb-(width>>1);
		if(start<0)
			start=0;
		int end=start+width;
		if(end>imsize)
			end=imsize;
		for(int k=start;k<end;++k)
			printf("%d", k%10);
		printf("\n");
		print_bitplane(buffer+start, end-start, sizeof(*buffer), kp);//
		print_bitplane(b2+start, end-start, sizeof(*buffer), kp);//
		//print_image(buffer, iw, ih, 0, iw, 0, ih);
		//print_image(b2, iw, ih, 0, iw, 0, ih);
	}
#ifdef FREE_ORIGINAL_IM
	stbi_image_free(original_image);
#endif
	delete[] b2;
#endif

	//AC - image
#if 0
/*	int iw=0, ih=0, nch=0;
#define	FREE_IMAGE
	byte *original_image=stbi_load(argv[1], &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("20211129 1 confidence.PNG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(original_image)
		printf("Opened \'%s\'\n", argv[1]);
	else
	{
		printf("Failed to open \'%s\'\n", argv[1]);
		return 1;
	}
	auto image=(int*)original_image;
#ifdef ENABLE_RVL
	const int depth=9;
#else
	const int depth=8;
#endif
	int imsize=iw*ih;//*/

	const char text[]="Sample text yes";
	//const int depth=8;
	const int depth=1;
	const int imsize=sizeof(text);
	const char *image=text;//*/

	short *buffer=new short[imsize];
	for(int k=0;k<imsize;++k)//extract red channel
	//	buffer[k]=image[k]&0xFF;
		buffer[k]=image[k]>>6&1;
#ifdef FREE_IMAGE
	stbi_image_free(original_image);
#endif//*/

/*	printf("Opening \'%s\'\n", argv[1]);
	std::string text;
	open_text(argv[1], text);
	//open_text("D:/C/ABAC Linux/Makefile", text);
	const int depth=8;
	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];//*/


/*	int imsize=1024;
	const int depth=1;
	short *buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=rand()&1;
		//buffer[k]=!(k&1);
		//buffer[k]=0;//*/


/*	unsigned iw=0, ih=0;
	std::vector<unsigned char> vimage;
	int error=lodepng::decode(vimage, iw, ih, "example.png");
	if(error)
	{
		printf("%s\n", lodepng_error_text(error));
		return 1;
	}
	auto image=(int*)vimage.data();
	int imsize=iw*ih;
//	int imsize=800;
	short *buffer=new short[imsize];
	const int depth=8;
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;//*/
	//memset(buffer, 0, imsize*sizeof(short));
	//for(int k=0;k<100000;++k)
	//	buffer[rand()%imsize]=rand();
	auto b2=new short[imsize];
	memset(b2, 0, imsize*sizeof(*b2));

	//print_rgba8((unsigned*)image.data(), imsize);
	//print_hex(buffer, imsize, depth);
	
#if 0
	vector_bool bits;
	int histogram[256];
	huff_encode(buffer, imsize, 8, histogram, bits, true);
	huff_decode(bits.data.data(), bits.bitSize, imsize, 8, histogram, b2, true);
#endif
#if 1
	std::string data;
	int sizes[depth]={};
	int conf[depth]={};
//	int probs[depth]={};

	//abac_estimate(buffer, imsize, depth, 2, true);
	
	//memset(buffer, 0, imsize*sizeof(*buffer));//
	//rans_encode(buffer, imsize, depth, sizeof(*buffer), data, sizes, conf, true);
	//rans_decode(data.data(), sizes, conf, b2, imsize, depth, sizeof(*b2), true);
	
//#ifdef ENABLE_RVL
//	differentiate_rows(buffer, iw, ih, buffer);
//#endif
//	abac3_encode(buffer, imsize, depth, data, sizes, conf, true);
//	abac3_decode(data.data(), sizes, conf, b2, imsize, depth, true);
//#ifdef ENABLE_RVL
//	integrate_rows(buffer, iw, ih, buffer);
//	integrate_rows(b2, iw, ih, b2);
//#endif

	//abac3_encode(buffer, imsize, depth, data, sizes, conf, true);
	//abac3_decode(data.data(), sizes, conf, b2, imsize, depth, true);

	abac2_encode(buffer, imsize, depth, data, sizes, conf, true);
	abac2_decode(data.data(), sizes, conf, b2, imsize, depth, true);

	//abac_encode(buffer, imsize, depth, data, sizes, true);
	//abac_decode(data.data(), sizes, b2, imsize, depth, true);

	//abac_encode_sse2(buffer, imsize, depth, data, sizes, true);
	//abac_decode_sse2(data.data(), sizes, b2, imsize, depth, true);

	//abac_encode_avx2(buffer, imsize, depth, data, sizes, true);
	//abac_decode_avx2(data.data(), sizes, b2, imsize, depth, true);

	//ac_encode(buffer, imsize, depth, data, sizes, probs, true);
	//ac_decode(data.data(), sizes, probs, b2, imsize, depth, true);
	
	//ac_debug(buffer, imsize, depth, data, sizes, probs, b2, true);
#endif


	//print_hex(b2, imsize, depth);

	int nerrors=0, kp=0, kb=0;
	int depthmask=(1<<depth)-1;
	//printf("Depthmask: %04X\n", depthmask);
	for(int k=0;k<imsize;++k)
	{
		if((b2[k]^buffer[k])&depthmask)
		{
			if(!nerrors)
			{
				kb=k;
				for(kp=0;kp<depth&&!((b2[k]^buffer[k])>>kp&1);++kp);
			}
			++nerrors;
			printf("%d: %04X != %04X\n", k, (int)buffer[k], (int)b2[k]);
			//if(nerrors>=100)
				break;
		}
	}
	//printf("original, recovered, xor:\n");
	//for(int k=0;k<imsize;++k)
	//{
	//	printf("%3d ", k);
	//	print_bin(buffer[k], depth);
	//	printf(", ");
	//	print_bin(b2[k], depth);
	//	printf(", ");
	//	print_bin(buffer[k]^b2[k], depth);
	//	printf("\n");
	//}

	//for(int k=0;k<depth;++k)//print all bitplanes
	//{
	//	printf("Plane %d:\n", k);
	//	print_bitplane(buffer, imsize, k);//
	//	print_bitplane(b2, imsize, k);//
	//}
	//printf("\n");

	if(nerrors)
	{
		const int width=64;
		printf("Error in bit plane %d, pixel %d:\n", kp, kb);
		int start=kb-(width>>1);
		if(start<0)
			start=0;
		int end=start+width;
		if(end>imsize)
			end=imsize;
		for(int k=start;k<end;++k)
			printf("%d", k%10);
		printf("\n");
		print_bitplane(buffer+start, end-start, sizeof(*buffer), kp);//
		print_bitplane(b2+start, end-start, sizeof(*b2), kp);//
		for(int k=0;k<imsize;++k)
		{
			printf("%04X ", buffer[k]);
			for(int k2=depth-1;k2>=0;--k2)
				printf("%d", buffer[k]>>k2&1);
			printf(" ");
			printf("%04X ", b2[k]);
			for(int k2=depth-1;k2>=0;--k2)
				printf("%d", b2[k]>>k2&1);
			printf("\n");
		}
	}

	//print_bitplane(buffer, imsize, 0);//
	//print_bitplane(b2, imsize, 0);//

	//print_bitplane(buffer+2800, 100, 7);//
	//print_bitplane(b2+2800, 100, 7);//
	//print_bitplane(buffer+700, 100, 0);//
	//print_bitplane(b2+700, 100, 0);//
	//print_bitplane(buffer+215000, 100, 0);//
	//print_bitplane(b2+215000, 100, 0);//
	
	delete[] b2;
	delete[] buffer;
#endif

	//AC text
#if 0
	const char str[]="111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100000111111";
//	const char str[]="000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="111100000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="11111111111111110100000001";
//	const char str[]="00111111100000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="011110100000100011101110000011010010110100110100101011000110111101111011000111001110110111001100111";
//	const char str[]="1011110011010110000010110001111000111010111101001010100100011101";
//	const char str[]="0000000100000001000000010000000100000001000000010000000100000001";
	const int imsize=sizeof(str);
	short buffer[imsize]={};
	for(int k=0;k<imsize;++k)
		buffer[k]=str[k]-'0';
	//	buffer[k]=rand()&1;

	//int dmask=0;
	//ac_test_bitplane_differentiation(buffer, imsize, 1, dmask);
	//if(dmask)
	//	ac_differentiate_bitplanes(buffer, imsize, 1, dmask);
	
	const int depth=1;
	std::string data;
	int sizes[depth], probs[depth];
	short b2[imsize]={};

	//for(int k=0;k<imsize-1;++k)
	//{
	//	printf("k=%d\n", k);
	//	memset(b2, 0, sizeof(b2));
	//	ac_debug(buffer+k, imsize-k, 1, data, sizes, probs, b2, true);
	//}
	ac_debug(buffer, imsize, 1, data, sizes, probs, b2, true);
//	ac_encode(buffer, imsize, 1, data, sizes, true);
//	ac_decode(data.c_str(), sizes.data(), b2, imsize, 1, true);

	//if(dmask)
	//{
	//	ac_integrate_bitplanes(buffer, imsize, 8, dmask);
	//	ac_integrate_bitplanes(b2, imsize, 8, dmask);
	//}

	print_bitplane(buffer, imsize, 0);
	print_bitplane(b2, imsize, 0);//*/
#endif

	//AC
#if 0
	const int imsize=100;
	short buffer[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=1;
	buffer[62]&=0;
	std::string data;
	int sizes[imsize], probs[imsize];
	short b2[imsize]={};
	ac_debug(buffer, imsize, 1, data, sizes, probs, b2, true);
	print_bitplane(buffer, imsize, 0);
	print_bitplane(b2, imsize, 0);
#endif

	//AC - text
#if 0
	//std::string text=R"(void			ac_print_summary())";//*/
	
	std::string text;
	for(int k=0;k<100;++k)//100
		text.push_back(0xFF);
	text[62]&=0;
//	text[62]&=0xDF;//62
	
	//std::string text;
	//if(!open_text("Makefile", text))
	//{
	//	printf("Failed to open file\n");
	//	return 1;
	//}

	//text.resize(100);
	//for(int k=0;k<text.size();++k)
	//	text[k]=rand()&1;

	//memset(&text[0], 'A', text.size());

	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];
	
//	const int depth=8;
	const int depth=1;

	//printf("Bitplane 2:\n");
	//print_bitplane(buffer, imsize, 2);
	int dmask=0;
	ac_test_bitplane_differentiation(buffer, imsize, depth, dmask);
	dmask=0;//
	if(dmask)
	{
		printf("dmask: ");
		print_bin(dmask, depth);
		printf("\n");
	//	ac_differentiate_bitplanes(buffer, imsize, depth, dmask);
	}//*/
/*	int bitplane=0;
	int p1[depth]={};
	for(int k=0;k<imsize;++k)
		for(int k2=0;k2<depth;++k2)
			p1[k2]+=buffer[k]>>k2&1;
	for(int k=0, latch=0;k<imsize;++k)//differentiator
	{
		if(k<imsize-1)
			latch=buffer[k]^buffer[k+1];
		else
			latch=buffer[k];
		buffer[k]=latch;
	}
	int p2[depth]={};
	for(int k=0;k<imsize;++k)
		for(int k2=0;k2<depth;++k2)
			p2[k2]+=buffer[k]>>k2&1;
		printf("One\'s count: (size=%d)\n", imsize);
	for(int k2=0;k2<depth;++k2)
		printf("%d: %d -> %d (50 + %lf%% -> 50 + %lf%%)\n", k2, p1[k2], p2[k2], 50-100.*p1[k2]/imsize, 50-100.*p2[k2]/imsize);
	print_bitplane(buffer, imsize, bitplane);//*/

	std::string data;
	int sizes[depth]={}, prob[depth]={};
	ac_encode(buffer, imsize, depth, data, sizes, prob, true);

	auto b2=new short[imsize];
	ac_decode(data.c_str(), sizes, prob, b2, imsize, depth, true);
	if(dmask)
	{
	//	ac_integrate_bitplanes(buffer, imsize, depth, dmask);
	//	ac_integrate_bitplanes(b2, imsize, depth, dmask);
	}
	for(int k=0, error_count=0;k<imsize;++k)
	{
		int error=(buffer[k]^b2[k])&((1<<depth)-1);
		if(error)
		{
			printf("Error %d: ", k);
			print_bin(buffer[k], depth);
			printf("!=");
			print_bin(b2[k], depth);
			printf(", xor=");
			print_bin(error, depth);
			printf(",%c!=%c\n", (char)buffer[k], (char)b2[k]);
			//printf("Error %d: %04X!=%04X,%3d!=%3d,%c!=%c\n", k, (int)buffer[k], (int)b2[k], (int)buffer[k], (int)b2[k], (char)buffer[k], (char)b2[k]);
			++error_count;
			if(error_count>=100)
				break;
		}
	}
	int printsize=imsize;
	if(printsize>100)
		printsize=100;
	printf("Original: ");
	print_bitplane(buffer, printsize, 0);//5
	printf("Result:   ");
	print_bitplane(b2, printsize, 0);//5

	//std::vector<Symbol> data;
	//ac_encode(buffer, imsize, depth, data, true);
	//
	//ac_decode(data.data(), buffer, imsize, depth, true);
	std::string str(printsize, '\0');
	for(int k=0;k<(int)str.size();++k)
		str[k]=b2[k];
	printf("Result:\n%.*s\n", printsize, str.c_str());
/*	for(int k=0;k<head;++k)
	{
		printf("text:\t");
		print_bin(text[k]&((1<<depth)-1), depth);
		printf(" %d %c\n", text[k], text[k]);

		printf("dec:\t");
		print_bin(str[k]&((1<<depth)-1), depth);
		printf(" %d %c\n\n", str[k], str[k]);
	}
	//	printf("%02X: %c\n", str[k]&255, str[k]);//*/

	delete[] buffer;
#endif

	//AC v1 - raw PNG
#if 0
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Failed to open image\n");
		return 1;
	}
	int imsize=iw*ih;
	short *red=new short[imsize];
	for(int k=0;k<imsize;++k)//extract red channel
		red[k]=original_image[k]&0xFF;
	stbi_image_free(original_image);
	
	std::vector<Symbol> data;
	
	//auto t1=__rdtsc();
	ac_encode(red, imsize, 8, data, true);
	//auto t2=__rdtsc();

/*	int compressed_bytesize=(int)data.size()*sizeof(Symbol);
	printf("u%d/u%d:\n", SYMBOL_BITS, SYMBOL_BITS*2);
	printf("Size: %d -> %d, ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
	//printf("Size: %dx%d=%d -> %d bytes\n", iw, ih, imsize, compressed_bytesize);
	//printf("Compresison ratio: %lf\n", (double)imsize/compressed_bytesize);
	printf("Elapsed: %lld cycles\n", t2-t1);//*/
	delete[] red;
#endif

#ifdef _MSC_VER
	printf("\nDone.\n");
	int i=0;
	scanf_s("%d", &i);
#endif
	return 0;
}