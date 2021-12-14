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
#else
#include<intrin.h>
#endif
#include<time.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"

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
void			print_image(int *buffer, int bw, int bh, int x1, int x2, int y1, int y2)
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
void			print_fdata(float *data, int bw, int bh, int x1, int x2, int y1, int y2)
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

void			extract_channel(const int *buffer, int imsize, int channel, float *data)//[1, -1[
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
		*p=128+int(data[k]*128);
}
void			apply_DWT_1D(float *data, int size, int stride, float *temp, double *filt, int filtsize, double *norms)
{
	int nodd=size>>1, neven=nodd+(size&1), s2=stride<<1, bsize=size*stride;
	for(int k=0;k<filtsize;k+=2)
	{
		double predict=filt[k], update=filt[k+1];
		int kx=0;
		//predict
		for(;kx+s2<bsize;kx+=s2)
			data[kx+stride]+=float(predict*(data[kx]+data[kx+s2]));
		if(kx+stride<bsize)
			data[kx+stride]+=float(predict*(data[kx]+data[kx]));
			
		//update
		data[0]+=float(update*(data[stride]+data[stride]));
		for(kx=s2;kx+s2<bsize;kx+=s2)
			data[kx]+=float(update*(data[kx-stride]+data[kx+stride]));
	}
	for(int kx=0, kx2=0;kx<neven;++kx, kx2+=s2)//split even & odd samples
		temp[kx]=float(norms[0]*data[kx2]);
	for(int kx=neven, kx2=stride;kx<size;++kx, kx2+=s2)
		temp[kx]=float(norms[1]*data[kx2]);
	if(stride==1)
		memcpy(data, temp, size*sizeof(*data));
	else
	{
		for(int ks=0, kd=0;ks<size;kd+=stride, ++ks)
			data[kd]=temp[ks];
	}
}
void			apply_invDWT_1D(float *data, int size, int stride, float *temp, double *filt, int filtsize, double *norms)
{
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
	for(int kx=neven, kx2=stride;kx<size;++kx, kx2+=s2)
		data[kx2]=float(oddnorm*temp[kx]);
	for(int k=filtsize-2;k>=0;--k)
	{
		double predict=filt[k], update=filt[k+1];
		int kx=0;
		//update
		data[0]-=float(update*(data[stride]+data[stride]));
		for(kx=s2;kx+s2<bsize;kx+=s2)
			data[kx]-=float(update*(data[kx-stride]+data[kx+stride]));

		//predict
		for(;kx+s2<bsize;kx+=s2)
			data[kx+stride]-=float(predict*(data[kx]+data[kx+s2]));
		if(kx+stride<bsize)
			data[kx+stride]-=float(predict*(data[kx]+data[kx]));
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
void			apply_invDWT_2D(float *buffer, int bw, int bh, double *filt, int filtsize, double *norms)//lifting scheme, filtsize is even, norms={lonorm, hinorm}
{
	int maxdim=bw;
	if(maxdim<bh)
		maxdim=bh;
	auto temp=new float[maxdim];
	for(int xsize=bw, ysize=bh;xsize>2&&ysize>2;)
	{
		for(int ky=0;ky<ysize;++ky)
			apply_invDWT_1D(buffer+bw*ky, xsize, 1, temp, filt, filtsize, norms);
		for(int kx=0;kx<xsize;++kx)
			apply_invDWT_1D(buffer+kx, ysize, bw, temp, filt, filtsize, norms);

		xsize-=xsize>>1, ysize-=ysize>>1;
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
	printf("Enter ZERO to save result image: ");
	int x=0;
	scanf_s("%d", &x);
	if(!x)
		lodepng::encode(filename, (const byte*)buffer, iw, ih);
}
void			gen_filename()
{
	time_t t=time(nullptr);
	tm now={};
	localtime_s(&now, &t);
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
int				main(int argc, char **argv)
{
	set_console_buffer_size(120, 4000);
	if(argc<2)
	{
		printf("Pass filename as command argument\n");
		return 1;
	}

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
#if 1
#if 0
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
		-1.5861343420594238292020515785937243,
		-0.052980118573376671019344897514278511,
		 0.882911075528503100647487289764588328,
		 0.443506852044983007635941975182636061,

		 1.1496043988602962433651033986614476,
		 0.86986445162473959153241758552174375,
	};

	auto data=new float[imsize];
	for(int kc=0;kc<3;++kc)
	{
		printf("Processing channel %d...\n", kc);
		extract_channel(buffer, imsize, kc, data);
			printf("Before DWT:\n"), print_fdata(data, iw, ih, 0, 8, 0, 8);
		apply_DWT_2D	(data, iw, ih, filt, 4, filt+4);
			printf("DWT:\n"), print_fdata(data, iw, ih, 0, 8, 0, 8);
		apply_invDWT_2D	(data, iw, ih, filt, 4, filt+4);
			printf("Inverse DWT:\n"), print_fdata(data, iw, ih, 0, 8, 0, 8);
		assign_channel(data, imsize, kc, b2);
	}

	//gen_filename();
	//save_image(g_buf, b2, iw, ih);
	print_image(buffer, iw, ih, 0, iw, 0, ih);
	print_image(b2, iw, ih, 0, iw, 0, ih);

	delete[] data;
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
#endif
#if 0
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

	std::string data;
	int sizes[depth]={};
	int conf[depth]={};
	abac2_encode(b2, imsize, depth, sizeof(*b2), data, sizes, conf, true);
	abac2_decode(data.data(), sizes, conf, b2, imsize, depth, sizeof(*b2), true);

	//apply_invRCT29_27(b2, b2, imsize);
	//invRVL_rows(b2, iw, ih, b2);
	//apply_invRCT26_24(b2, b2, imsize);//*/
		//printf("IRCT26:\n"), print_image(b2, iw, ih, 0, iw, 0, ih);

	//apply_invDWT24(b2, iw, ih);
		//printf("IDWT:\n"), print_image(b2, iw, ih, 0, iw, 0, ih);
	
	int nerrors=0, kp=0, kb=0;
	int depthmask=(1<<depth0)-1;
	//printf("Depthmask: %04X\n", depthmask);
	for(int k=0;k<imsize;++k)
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
	int iw=0, ih=0, nch=0;
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

/*	const char text[]="Sample text";
	const int depth=8;
	int imsize=sizeof(text);
	const char *image=text;//*/

	short *buffer=new short[imsize];
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;
	stbi_image_free(original_image);//*/

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
	
#ifdef ENABLE_RVL
	differentiate_rows(buffer, iw, ih, buffer);
#endif
	abac3_encode(buffer, imsize, depth, data, sizes, conf, true);
	abac3_decode(data.data(), sizes, conf, b2, imsize, depth, true);
#ifdef ENABLE_RVL
	integrate_rows(buffer, iw, ih, buffer);
	integrate_rows(b2, iw, ih, b2);
#endif

	//abac3_encode(buffer, imsize, depth, data, sizes, conf, true);
	//abac3_decode(data.data(), sizes, conf, b2, imsize, depth, true);

	//abac2_encode(buffer, imsize, depth, data, sizes, conf, true);
	//abac2_decode(data.data(), sizes, conf, b2, imsize, depth, true);

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
		print_bitplane(buffer+start, end-start, kp);//
		print_bitplane(b2+start, end-start, kp);//
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