//ac.cpp - Arithmetic Coder implementation
//Copyright (C) 2021  Ayman Wagih Mohsen, unless source link provided.
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

#include"ac.h"
#include<stdio.h>
#include<string.h>
#include<vector>

#ifdef __linux__
#define	__rdtsc	__builtin_ia32_rdtsc
#else
#include<intrin.h>
#endif

	#define	LOG_WINDOW_SIZE		16	//[2, 16]

	#define	PRINT_ERROR//should be enabled

//	#define	PRINT_VISUAL	//disable with release & large files
//	#define	PRINT_BITPLANES	//ditto
//	#define	DEBUG_PREDICTOR	//ditto

#ifdef PRINT_VISUAL
const int visual_plane=0, visual_start=0, visual_end=100;
//const int visual_plane=1, visual_start=215000, visual_end=215100;
//const int visual_plane=0, visual_start=700, visual_end=800;
#endif
void			print_bitplane(const short *buffer, int imsize, int bitplane);
typedef unsigned long long u64;

void			breakpoint()
{
	int x=0;
	x=1;
}
//void			ac_print_summary(int nsymbols, int original_size, long long elapsed, bool encode)
//{
//	int compressed_bytesize=nsymbols*sizeof(Symbol);
//	printf("u%d/u%d %s:\n", SYMBOL_BITS, SYMBOL_BITS*2, encode?"encode":"decode");
//	printf("Size: %d -> %d, ratio: %lf\n", original_size, compressed_bytesize, (double)original_size/compressed_bytesize);
//	printf("Elapsed: %lld cycles\n", elapsed);
//}
//static char	history[0x2000];
inline int		mod(int x, int n)
{
	x%=n;
	x+=n&-(x<0);
	return x;
}
int				floor_log2(unsigned long long n)
{
	int logn=0;
	int sh=(n>=1ULL<<32)<<5;logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
int				clamp(int lo, int x, int hi)
{
	if(x<lo)
		x=lo;
	if(x>hi)
		x=hi;
	return x;
}

#ifdef PRINT_VISUAL
struct	Range
{
	unsigned start, end, middle;
	int bit;
	Range():start(0), end(0), middle(0), bit(0){}
	Range(unsigned start, unsigned end, unsigned middle, int bit):start(start), end(end), middle(middle), bit(bit){}
};
std::vector<Range> ranges;
#endif
static void		print_range(unsigned start, unsigned end, unsigned middle, unsigned code, int kb, int bit)
{
	const int width=48;
	//const int width=64;
	//const int width=128;

	if(start>end)
		printf("Fatal error: start>end\n");

	if(middle<start)
		printf("middle < start\n");
	if(middle>=end)
		printf("middle >= end\n");

	if(code<start)
		printf("code < start\n");
	if(code>end)
		printf("code > end\n");

	int vstart=(u64)start*width/0xFFFFFFFF, vend=(u64)end*width/0xFFFFFFFF, vmiddle=(u64)middle*width/0xFFFFFFFF;
	vstart=clamp(0, vstart, width);
	vend=clamp(0, vend, width);
	vmiddle=clamp(0, vmiddle, width);

	for(int k=0;k<vstart;++k)
		printf("-");
	for(int k=vstart;k<vmiddle;++k)
		printf("0");
	for(int k=vmiddle;k<vend;++k)
		printf("1");
	for(int k=vend;k<width;++k)
		printf("-");
	printf(" %08X~%08X", start, end);
	printf("\n");
//	printf("\t");
	if(start<end)//zoomed view
	{
		int ratio=(u64)(middle-start)*width/(end-start);
		ratio=clamp(0, ratio, width);
		for(int k=0;k<ratio;++k)
			printf("0");
		for(int k=ratio;k<width;++k)
			printf("1");
		printf("\n");
		for(int k=0;k<ratio;++k)
			printf(" ");
		printf("^ %08X\n", middle);
		//if(code!=middle)
		{
			ratio=(u64)(code-start)*width/(end-start);
			ratio=clamp(0, ratio, width);
			for(int k=0;k<ratio;++k)
				printf(" ");
			printf("^ %08X code, bit %d = %d\n", code, kb, bit);
		}
	}
//	printf("\n%08X\n%08X\n%08X\n", start, middle, end);//
	printf("\n");
}
inline void		update_histogram(short value, int depth, int *h)
{
	for(int k2=0;k2<depth;++k2)
		h[k2]+=value>>k2&1;
}
void			ac_test_bitplane_differentiation(short *buffer, int imsize, int depth, int &dmask)
{
	int p1[16], p2[16];
	for(int k=0, latch=0;k<imsize;++k)//differentiator
	{
		if(k<imsize-1)
			latch=buffer[k]^buffer[k+1];
		else
			latch=buffer[k];
		update_histogram(buffer[k], depth, p1);
		update_histogram(latch, depth, p2);
	}
	dmask=0;
	for(int k2=0, halfsize=imsize>>1;k2<depth;++k2)
		dmask|=(abs(halfsize-p2[k2])>abs(halfsize-p1[k2]))<<k2;
}
void			ac_differentiate_bitplanes(short *buffer, int imsize, int depth, int dmask)
{
	for(int k=0;k<imsize-1;++k)//differentiator
	{
		int latch=buffer[k]^buffer[k+1];
		buffer[k]=latch&dmask|buffer[k]&~dmask;
	}
}
void			ac_integrate_bitplanes(short *buffer, int imsize, int depth, int dmask)
{
	for(int k=imsize-1;k>0;--k)
	{
		int latch=buffer[k]^buffer[k-1];
		buffer[k-1]=latch&dmask|buffer[k-1]&~dmask;
	}
}
const int		window_size=1<<LOG_WINDOW_SIZE, prob_mask=window_size-1, prob_max=window_size-2, prob_init=(1<<(LOG_WINDOW_SIZE-1))-1;
void			ac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, bool loud)
{
#ifdef PRINT_VISUAL
	ranges.clear();
	int vkc_start=0, vkc_end=0;
#endif
	if(!imsize)
		return;
	auto t1=__rdtsc();

	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		auto &plane=planes[depth-1-kp];
		plane.reserve(imsize>>8);

		//int prob=prob_init;
		int prob=0;
		for(int k=0;k<imsize;++k)//calculate probability of one
			prob+=buffer[k]>>kp&1;
		prob=(u64)prob*(1<<LOG_WINDOW_SIZE)/imsize;
		prob=(1<<LOG_WINDOW_SIZE)-prob;//probability of zero
		prob=clamp(1, prob, prob_max);//avoid exact 0 and 100% probability
		out_probabilities[depth-1-kp]=prob;
		
		//u64 start=0, end=0x100000000;
		unsigned start=0, end=0xFFFFFFFF;
		for(int kb=0;kb<imsize;)//bit-pixel loop		http://mattmahoney.net/dc/dce.html
		{
			u64 range=end-start;
			if(range==1)
			{
				plane.push_back(start>>24);
				plane.push_back(start>>16&0xFF);
				plane.push_back(start>>8&0xFF);
				plane.push_back(start&0xFF);
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			unsigned middle=start+(unsigned)(range*prob>>LOG_WINDOW_SIZE);
			int bit=buffer[kb]>>kp&1;
			//int prevbit;
			//if(kb>=window_size)
			//	prob+=bit-(buffer[kb-window_size]>>kp&1);
			//if(kb<window_size)
			//	prevbit=kb&1;
			//else
			//	prevbit=buffer[kb-window_size]>>kp&1;
			//prob+=bit-prevbit;
			//prob=clamp(1, prob, prob_max);//probability shouldn't be exactly 0 or 100%

			//if(kp==5)
			//	if(kb>55&&kb<70)//error at plane 5, pixel 62
			//		print_range(start, end, middle, middle, kb, bit);//
#ifdef PRINT_VISUAL
			if(kp==visual_plane)
			{
				if(kb==visual_start)
					vkc_start=plane.size();
				if(kb==visual_end)
					vkc_end=plane.size();
				ranges.push_back(Range(start, end, middle, bit));
			}
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
				print_range(start, end, middle, middle, kb, bit);
#endif
			
			if(bit)
				start=middle;
			else
				end=middle-1;
			++kb;
			
			if((start^end)<0x1000000)//most significant byte has stabilized			zpaq 1.10
			{
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: ");//
#endif
				do
				{
#ifdef PRINT_VISUAL
					if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
						printf("%08X~%08X, ", start, end);//
#endif
					plane.push_back(start>>24);
					start<<=8;
					end=end<<8|0xFF;
				}while((start^end)<0x1000000);
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\n");//
#endif
			}
		}
		plane.push_back(start>>24&0xFF);//big-endian
		plane.push_back(start>>16&0xFF);
		plane.push_back(start>>8&0xFF);
		plane.push_back(start&0xFF);
#ifdef PRINT_BITPLANES
		if(kp==visual_plane)
		{
			//for(int k=vkc_start;k<vkc_end;++k)
			for(int k=vkc_start-1;k<vkc_end+10;++k)//
			//for(int k=0;k<(int)plane.size();++k)
				printf("%02X-", plane[k]&0xFF);
			printf("\n");
		}
#endif
	}
	auto t_enc=__rdtsc();
	//out_sizes.resize(depth);
	out_data.clear();
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
/*	for(int k=0, size=depth;k<depth;++k)
	{
		size+=planes[k].size();
		auto p=(unsigned char*)&size;
		data.push_back(p[0]);
		data.push_back(p[1]);
		data.push_back(p[2]);
		data.push_back(p[3]);
	}//*/
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)out_data.size()<<3;
		printf("AC encode:  %lld cycles (Enc: %lld cycles)\n", t2-t1, t_enc-t1);
		printf("Size: %d -> %d, ratio: %lf\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize);
		printf("Bit\tkp\tSymbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%2d\t%5d\n", depth-1-k, k, out_sizes[k]);
		
		printf("Preview:\n");
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
}
int				load32_big(const unsigned char *data)
{
	return data[0]<<24|data[1]<<16|data[2]<<8|data[3];
}
void			ac_decode(const char *data, const int *sizes, const int *probabilities, short *buffer, int imsize, int depth, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));
	
	for(int kp=depth-1, cusize=0;kp>=0;--kp)//bit-plane loop
	{
		int ncodes=sizes[depth-1-kp];
		auto plane=data+cusize;

		int prob=probabilities[depth-1-kp];
		//int prob=prob_init;

		//u64 start=0, end=0x100000000;
		unsigned start=0, end=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane);
		for(int kc=4, kb=0;kb<imsize;)//bit-pixel loop
		{
			u64 range=end-start;
			if(range==1)
			{
				code=load32_big((unsigned char*)plane+kc);
				kc+=4;
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			unsigned middle=start+(unsigned)(range*prob>>LOG_WINDOW_SIZE);
			int bit=code>=middle;
			//int prevbit;
		/*	if(code==middle)
			{
#ifdef PRINT_VISUAL
				printf("code==middle at kp=%d,kb=%d\n", kp, kb);
#endif
				int start2=start, end2=end, middle2=middle, code2=code, kc2=kc;
				for(;kc2<ncodes;)//look ahead on match
				{
					code2=code2<<8|(unsigned char)plane[kc2];
					++kc2;
					start2<<=8;
					end2=end2<<8|0xFF;
					middle2=((u64)start2+end2)>>1;
					//middle2=start2+((u64)(end2-start2)*(prob_mask-prob)>>LOG_WINDOW_SIZE);
					if(middle2!=code2)
						break;
				}
				bit=code2>=middle2;
			}//*/
			//if(kb>=window_size)
			//	prob+=bit-(buffer[kb-window_size]>>kp&1);
			//if(kb<window_size)
			//	prevbit=kb&1;
			//else
			//	prevbit=buffer[kb-window_size]>>kp&1;
			//prob+=bit-prevbit;
			//prob=clamp(1, prob, prob_max);//probability shouldn't be exactly 0 or 100%
			//if(code==middle)
			//	printf("(kp=%d,kb=%d)\n", kp, kb);
			//if(kp==5)
			//	if(kb>55&&kb<70)//error at plane 5, pixel 62
			//		print_range(start, end, middle, code, kb, bit);//BUG: sosmetimes wrong when code == middle
#ifdef PRINT_VISUAL
			//if(kp==0&&kb>700&&kb<800)//
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
			{
				if(kb<(int)ranges.size())
					print_range(ranges[kb].start, ranges[kb].end, ranges[kb].middle, code, kb, ranges[kb].bit);
				print_range(start, end, middle, code, kb, bit);
				printf("\n\n");
			}
#endif
			
			if(bit)
				start=middle;
			else
				end=middle-1;
			
			buffer[kb]|=bit<<kp;
			++kb;

			while(kc<ncodes&&(start^end)<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: %08X~%08X, code %08X\n", start, end, code);//
#endif
				code=code<<8|(unsigned char)plane[kc];
				++kc;
				start<<=8;
				end=end<<8|0xFF;
			}
		}
		cusize+=ncodes;
	}

	auto t2=__rdtsc();
	if(loud)
	{
		printf("AC decode:  %lld cycles\n", t2-t1);
	}
}
#if 0
const int
	ac_logwindowbits=5,
	ac_windowbits=1<<ac_logwindowbits,
	ac_windowbitsmask=ac_windowbits-1;
const Symbol
	ac_den=ac_windowbits+2,
	ac_invden=((Symbol)1<<SYMBOL_BITS/2)/ac_den,
	ac_range_end=(Symbol)1<<(sizeof(Symbol)*8-1)|1;
void			ac_encode(const short *buffer, int imsize, int depth, std::vector<Symbol> &data, bool loud)
{
	auto t1=__rdtsc();

	std::vector<std::vector<Symbol>> planes(depth);
	for(int kp=0;kp<depth;++kp)
	{
		auto &plane=planes[kp];

		int history[ac_windowbits];
		for(int k=0;k<ac_windowbits;++k)
			history[k]=k&1;
		Symbol num=(ac_windowbits>>1)+1;
		for(int kb=0;kb<imsize;)
		{
			Symbol start=0, end=ac_range_end;
			//Symbol start=0, end=0x100000000;
			for(;kb<imsize;)
			{
				auto middle=start+((MulType)(end-start)*(ac_den-num)*ac_invden>>SYMBOL_BITS/2);
			//	long long middle=start+((end-start)*(ac_den-num)*ac_invden>>SYMBOL_BITS/2);
				if(end-middle<=1||middle-start<=1)
					break;

				int bit=buffer[kb]>>kp&1;
				if(bit)
					start=middle;
				else
					end=middle;

				int &sieve=history[kb&ac_windowbitsmask];
				num+=bit-sieve;
				sieve=bit;

				++kb;
				//if(kp==2&&kb==60)
				//	int LOL_1=0;
			}
			plane.push_back(start);
		}
	}
	data.clear();
	for(int k=0, size=depth;k<depth;++k)
	{
		size+=planes[k].size();
		data.push_back(size);
	}
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		data.insert(data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int compressed_bytesize=(int)data.size()*sizeof(Symbol);
		printf("u%d/u%d encode:  %lld cycles\n", SYMBOL_BITS, SYMBOL_BITS*2, t2-t1);
		printf("Size: %d -> %d, ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
		printf("Plane\tCululative symbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%5d\n", k, data[k]);
	}
}
void			ac_decode(const Symbol *data, short *buffer, int imsize, int depth, bool loud)
{
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));
	
	for(int kp=0;kp<depth;++kp)
	{
		int ncodes=data[kp];
		auto plane=data+(kp?data[kp-1]:depth);

		int history[ac_windowbits];
		for(int k=0;k<ac_windowbits;++k)
			history[k]=k&1;
		Symbol num=(ac_windowbits>>1)+1;
		for(int kc=0, kb=0;kc<ncodes&&kb<imsize;++kc)
		{
			Symbol code=plane[kc];
			Symbol start=0, end=ac_range_end;
			//Symbol start=0, end=0x100000000;
			for(;kb<imsize;)
			{
				auto middle=start+((MulType)(end-start)*(ac_den-num)*ac_invden>>SYMBOL_BITS/2);
			//	Symbol middle=start+((end-start)*(ac_den-num)*ac_invden>>SYMBOL_BITS/2);
				if(end-middle<=1||middle-start<=1)
					break;

				int bit=code>=middle;
				if(bit)
					start=middle;
				else
					end=middle;

				int &sieve=history[kb&ac_windowbitsmask];
				num+=bit-sieve;
				sieve=bit;

				buffer[kb]|=bit<<kp;
				++kb;
				//if(kp==2&&kb==60)
				//	int LOL_1=0;
			}
		}
	}

	auto t2=__rdtsc();
	if(loud)
	{
		printf("u%d/u%d decode:  %lld cycles\n", SYMBOL_BITS, SYMBOL_BITS*2, t2-t1);
	}
}
#endif
void			ac_debug(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, short *out, bool loud)
{
	memset(out, 0, imsize*sizeof(short));
#ifdef PRINT_VISUAL
	ranges.clear();
	int vkc_start=0, vkc_end=0;
#endif
	if(!imsize)
		return;
	auto t1=__rdtsc();

#ifdef PRINT_ERROR
	int nerrors=0;
#endif
	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		auto &plane=planes[depth-1-kp];
		plane.reserve(imsize>>8);

		//int prob=prob_init;
		int prob=0;
		for(int k=0;k<imsize;++k)//calculate probability of one
			prob+=buffer[k]>>kp&1;
		prob=(u64)prob*(1<<LOG_WINDOW_SIZE)/imsize;
		prob=(1<<LOG_WINDOW_SIZE)-prob;//probability of zero
		prob=clamp(1, prob, prob_max);//avoid exact 0 and 100% probability
		out_probabilities[depth-1-kp]=prob;
		
		unsigned start=0, end=0xFFFFFFFF, middle=0;
		for(int kb=0;kb<imsize;)//bit-pixel loop		http://mattmahoney.net/dc/dce.html
		{
			u64 range=end-start;
			if(range==1)
			{
				plane.push_back(start>>24);
				plane.push_back(start>>16&0xFF);
				plane.push_back(start>>8&0xFF);
				plane.push_back(start&0xFF);
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			middle=start+(unsigned)(range*prob>>LOG_WINDOW_SIZE);
			//middle=start+(range*(prob_mask-prob)>>LOG_WINDOW_SIZE);
			int bit=buffer[kb]>>kp&1;
			//int prevbit;

#ifdef PRINT_VISUAL
			if(kp==visual_plane)
			{
				if(kb==visual_start)
					vkc_start=plane.size();
				if(kb==visual_end)
					vkc_end=plane.size();
				ranges.push_back(Range(start, end, middle, bit));
			}
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
				print_range(start, end, middle, middle, kb, bit);
#endif
			
			if(bit)
				start=middle;
				//start=middle+1;
			else
				end=middle-1;
				//end=middle;
			++kb;
			
			if((start^end)<0x1000000)//most significant byte has stabilized			zpaq 1.10
			{
				unsigned char codebyte=start>>24;
				if(codebyte==0)
					breakpoint();
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: ");//
#endif
				do
				{
#ifdef PRINT_VISUAL
					if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
						printf("%08X~%08X, ", start, end);//
#endif
					plane.push_back(start>>24);
					start<<=8;
					end=end<<8|0xFF;
				}while((start^end)<0x1000000);
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\n");//
#endif
			}
		}
		middle=start+((u64)(end-start)*prob>>LOG_WINDOW_SIZE);
		plane.push_back(middle>>24&0xFF);//big-endian
		plane.push_back(middle>>16&0xFF);
		plane.push_back(middle>>8&0xFF);
		plane.push_back(middle&0xFF);
		
#ifdef PRINT_BITPLANES
		if(kp==visual_plane)
		{
			//for(int k=vkc_start;k<vkc_end;++k)
			for(int k=vkc_start-1;k<vkc_end+10;++k)
			//for(int k=0;k<(int)plane.size();++k)
				printf("%02X-", plane[k]&0xFF);
			printf("\n");
		}
#endif

		int ncodes=plane.size();
		start=0, end=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane.c_str());
		for(int kc=4, kb=0;kb<imsize;)//bit-pixel loop
		{
			u64 range=end-start;
			if(range==1)
			{
				code=load32_big((unsigned char*)plane.data()+kc);
				kc+=4;
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			unsigned middle=start+(unsigned)(range*prob>>LOG_WINDOW_SIZE);
			//if(middle>=end)
			//	printf("Check: %08X~%08X, mid=%08X, prob=%d\n", start, end, middle, prob);
			int bit=code>=middle;
			//int prevbit;
		/*	if(code==middle)
			{
#ifdef PRINT_VISUAL
				printf("code==middle at kp=%d,kb=%d, %08X~%08X, mid %08X, code %08X\n", kp, kb, start, end, middle, code);
#endif
				int start2=start, end2=end, middle2=middle, code2=code, kc2=kc;
				do
				{
					code2=code2<<8|(unsigned char)plane[kc2];
					++kc2;
					start2<<=8;
					end2=end2<<8|0xFF;
					//middle2=((u64)start2+end2)>>1;
					middle2=start2+((u64)(end2-start2)*(prob_mask-prob)>>LOG_WINDOW_SIZE);
#ifdef PRINT_VISUAL
					printf("Look-ahead %d:\n", kc2-kc);
					print_range(start2, end2, middle2, code2, kb, code2>=middle2);
#endif
				}while(code2==middle2);
				bit=code2>=middle2;
			}//*/
#ifdef PRINT_VISUAL
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
			{
				if(kb<(int)ranges.size())
					print_range(ranges[kb].start, ranges[kb].end, ranges[kb].middle, code, kb, ranges[kb].bit);
				print_range(start, end, middle, code, kb, bit);
				printf("\n\n");
			}
#endif
			out[kb]|=bit<<kp;
#ifdef PRINT_ERROR
			int bit0=buffer[kb]>>kp&1;
#ifdef PRINT_VISUAL
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
#endif
			if(loud&&bit!=bit0&&nerrors<100)
			{
				printf("Error pixel %d, plane %d: original %d != decoded %d, prob=%d\n", kb, kp, bit0, bit, prob);
				int estart=kb-50, eend=kb+50;
				//int estart=kb-100, eend=kb+100;
				if(estart<0)
					estart=0;
				if(eend>imsize)
					eend=imsize;
				printf("Before:\t");
				print_bitplane(buffer+estart, eend-estart, kp);
				printf("After:\t");
				print_bitplane(out+estart, kb+1-estart, kp);
				print_range(start, end, middle, code, kb, bit);
				++nerrors;
			}
#endif
			
			if(bit)
				start=middle;
				//start=middle+1;
			else
				//end=middle;
				end=middle-1;
			
#ifdef PRINT_VISUAL
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
				printf("bit %d: %d\n", kb, bit);
#endif
			++kb;

			while(kc<ncodes&&(start^end)<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: %08X~%08X, code %08X\n", start, end, code);//
#endif
				code=code<<8|(unsigned char)plane[kc];
				++kc;
				start<<=8;
				end=end<<8|0xFF;
			}
		}
	}
	auto t_enc=__rdtsc();
	//out_sizes.resize(depth);
	out_data.clear();
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)out_data.size()<<3;
		printf("AC encode:  %lld cycles (Enc: %lld cycles)\n", t2-t1, t_enc-t1);
		printf("Size: %d -> %d, ratio: %lf\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize);
		printf("Bit\tkp\tSymbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%2d\t%5d\n", depth-1-k, k, out_sizes[k]);
		
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
}

#ifdef DEBUG_PREDICTOR
struct			ProbInfo
{
	int kb, kc, acc, p0;
	ProbInfo():kb(0), kc(0), acc(0), p0(0){}
	ProbInfo(int kb, int kc, int acc, int p0):kb(kb), kc(kc), acc(acc), p0(p0){}
};
static std::vector<ProbInfo> probs;
#endif
void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud)
{
#ifdef PRINT_VISUAL
	ranges.clear();
	int vkc_start=0, vkc_end=0;
#endif
	if(!imsize)
		return;
	auto t1=__rdtsc();

	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		auto &plane=planes[depth-1-kp];
		plane.reserve(imsize>>8);

		//int prob=prob_init;
		
		int prob=1<<(LOG_WINDOW_SIZE-1);//cheap weighted average predictor
		//int prob=0x55555555&prob_init;

		//int prob=0;
		//for(int k=0;k<imsize;++k)//calculate probability of one
		//	prob+=buffer[k]>>kp&1;
		//prob=(u64)prob*(1<<LOG_WINDOW_SIZE)/imsize;
		//prob=(1<<LOG_WINDOW_SIZE)-prob;//probability of zero

		//prob=clamp(1, prob, prob_max);//avoid exact 0 and 100% probability
		//out_probabilities[depth-1-kp]=prob;
		
		//u64 start=0, end=0x100000000;
		unsigned start=0, end=0xFFFFFFFF;
		for(int kb=0;kb<imsize;)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
		{
			//if(kp==0&&kb==6)
			//	int LOL_1=0;
			int bit=buffer[kb]>>kp&1;
			//if(bit)//impossible non-causal predictor
			//	prob=65536/8;
			//else
			//	prob=65536*7/8;

			u64 range=end-start;
			if(range==1)
			{
				plane.push_back(start>>24);
				plane.push_back(start>>16&0xFF);
				plane.push_back(start>>8&0xFF);
				plane.push_back(start&0xFF);
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}

			unsigned p0=prob;
			if(plane.size())
			{
				double ratio=(double)kb/(plane.size()<<3);
				p0=0x7FFF+(int)(((double)p0-0x7FFF)*ratio/(1+ratio));
#ifdef DEBUG_PREDICTOR
				if(kp==0)
					printf("r=%lf, acc=%04X, p0=%d\n", ratio, prob, p0);
#endif
				//if(!(kb%10000))
				//	printf("kp=%d, kb=%d, psize=%d, r=%lf, c=%lf, acc=%04X, P(0)=%d\n", kp, kb, plane.size(), ratio, ratio/(1+ratio), prob, p0);

				//int ratio=(kb<<3)/plane.size();
				//p0=0x7FFF+(int)((u64)((p0-0x7FFF)*ratio)/(1+ratio));//taking confidence into account
			}
			p0=clamp(1, p0, prob_max);
#ifdef DEBUG_PREDICTOR
			if(kp==0)
				probs.push_back(ProbInfo(kb, plane.size(), prob, p0));
#endif
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);

			//unsigned middle=start+(unsigned)(range*clamp(1, prob, prob_max)>>LOG_WINDOW_SIZE);
			//unsigned middle=start+(unsigned)(range*prob>>LOG_WINDOW_SIZE);
			//int prevbit;
			//if(kb>=window_size)
			//	prob+=bit-(buffer[kb-window_size]>>kp&1);
			//if(kb<window_size)
			//	prevbit=kb&1;
			//else
			//	prevbit=buffer[kb-window_size]>>kp&1;
			//prob+=bit-prevbit;

			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16

			//prob=clamp(1, prob, prob_max);//probability shouldn't be exactly 0 or 100%

			//if(kp==5)
			//	if(kb>55&&kb<70)//error at plane 5, pixel 62
			//		print_range(start, end, middle, middle, kb, bit);//
#ifdef PRINT_VISUAL
			if(kp==visual_plane)
			{
				if(kb==visual_start)
					vkc_start=plane.size();
				if(kb==visual_end)
					vkc_end=plane.size();
				ranges.push_back(Range(start, end, middle, bit));
			}
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
				print_range(start, end, middle, middle, kb, bit);
#endif
			
			if(bit)
				start=middle;
			else
				end=middle-1;
			++kb;
			
			if((start^end)<0x1000000)//most significant byte has stabilized			zpaq 1.10
			{
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: ");//
#endif
				do
				{
#ifdef PRINT_VISUAL
					if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
						printf("%08X~%08X, ", start, end);//
#endif
					plane.push_back(start>>24);
					start<<=8;
					end=end<<8|0xFF;
				}while((start^end)<0x1000000);
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\n");//
#endif
			}
		}
		plane.push_back(start>>24&0xFF);//big-endian
		plane.push_back(start>>16&0xFF);
		plane.push_back(start>>8&0xFF);
		plane.push_back(start&0xFF);
#ifdef PRINT_BITPLANES
		if(kp==visual_plane)
		{
			//for(int k=vkc_start;k<vkc_end;++k)
			for(int k=vkc_start-1;k<vkc_end+10;++k)//
			//for(int k=0;k<(int)plane.size();++k)
				printf("%02X-", plane[k]&0xFF);
			printf("\n");
		}
#endif
	}
	auto t_enc=__rdtsc();
	//out_sizes.resize(depth);
	out_data.clear();
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
/*	for(int k=0, size=depth;k<depth;++k)
	{
		size+=planes[k].size();
		auto p=(unsigned char*)&size;
		data.push_back(p[0]);
		data.push_back(p[1]);
		data.push_back(p[2]);
		data.push_back(p[3]);
	}//*/
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)out_data.size()<<3;
		printf("AC encode:  %lld cycles (Enc: %lld cycles)\n", t2-t1, t_enc-t1);
		printf("Size: %d -> %d, ratio: %lf\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize);
		printf("Bit\tkp\tSymbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%2d\t%5d\n", depth-1-k, k, out_sizes[k]);
		
		printf("Preview:\n");
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
}
void			abac_decode(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));
	
	for(int kp=depth-1, cusize=0;kp>=0;--kp)//bit-plane loop
	{
		int ncodes=sizes[depth-1-kp];
		auto plane=data+cusize;
		
		int prob=1<<(LOG_WINDOW_SIZE-1);//cheap weighted average predictor
		//int prob=0x55555555&prob_init;
		//int prob=probabilities[depth-1-kp];
		//int prob=prob_init;

		//u64 start=0, end=0x100000000;
		unsigned start=0, end=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane);
		for(int kc=4, kb=0;kb<imsize;)//bit-pixel loop
		{
			//if(kp==0&&kb==6)
			//	int LOL_1=0;
			u64 range=end-start;
			if(range==1)
			{
				code=load32_big((unsigned char*)plane+kc);
				kc+=4;
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			unsigned p0=prob;
			if(kc-4>0)
			{
				double ratio=(double)kb/((kc-4)<<3);
				p0=0x7FFF+(int)(((double)p0-0x7FFF)*ratio/(1+ratio));
#ifdef DEBUG_PREDICTOR
				if(kp==0)
					printf("r=%lf, acc=%04X, p0=%d\n", ratio, prob, p0);
#endif
			}
			p0=clamp(1, p0, prob_max);
#ifdef DEBUG_PREDICTOR
			if(kp==0)
			{
				if(p0!=probs[kb].p0)
					printf("Error kp=%d, kb=%d: kb %d->%d, kc %d->%d, acc %04X->%04X, p0%d->%d\n", kp, kb, probs[kb].kb, kb, probs[kb].kc, kc-4, probs[kb].acc, prob, probs[kb].p0, p0);
				//	printf("Error kp=%d, kb=%d: original p0 %d != calc p0 %d\n", kp, kb, probs[kb], p0);
			}
#endif
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
			int bit=code>=middle;
			//int prevbit;
		/*	if(code==middle)
			{
#ifdef PRINT_VISUAL
				printf("code==middle at kp=%d,kb=%d\n", kp, kb);
#endif
				int start2=start, end2=end, middle2=middle, code2=code, kc2=kc;
				for(;kc2<ncodes;)//look ahead on match
				{
					code2=code2<<8|(unsigned char)plane[kc2];
					++kc2;
					start2<<=8;
					end2=end2<<8|0xFF;
					middle2=((u64)start2+end2)>>1;
					//middle2=start2+((u64)(end2-start2)*(prob_mask-prob)>>LOG_WINDOW_SIZE);
					if(middle2!=code2)
						break;
				}
				bit=code2>=middle2;
			}//*/
			//if(kb>=window_size)
			//	prob+=bit-(buffer[kb-window_size]>>kp&1);
			//if(kb<window_size)
			//	prevbit=kb&1;
			//else
			//	prevbit=buffer[kb-window_size]>>kp&1;
			//prob+=bit-prevbit;

			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16

			//prob=clamp(1, prob, prob_max);//probability shouldn't be exactly 0 or 100%
			//if(code==middle)
			//	printf("(kp=%d,kb=%d)\n", kp, kb);
			//if(kp==5)
			//	if(kb>55&&kb<70)//error at plane 5, pixel 62
			//		print_range(start, end, middle, code, kb, bit);//BUG: sosmetimes wrong when code == middle
#ifdef PRINT_VISUAL
			//if(kp==0&&kb>700&&kb<800)//
			if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
			{
				if(kb<(int)ranges.size())
					print_range(ranges[kb].start, ranges[kb].end, ranges[kb].middle, code, kb, ranges[kb].bit);
				print_range(start, end, middle, code, kb, bit);
				printf("\n\n");
			}
#endif
			
			if(bit)
				start=middle;
			else
				end=middle-1;
			
			buffer[kb]|=bit<<kp;
			++kb;

			while(kc<ncodes&&(start^end)<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef PRINT_VISUAL
				if(kp==visual_plane&&kb>=visual_start&&kb<visual_end)
					printf("\nShift-out: %08X~%08X, code %08X\n", start, end, code);//
#endif
				code=code<<8|(unsigned char)plane[kc];
				++kc;
				start<<=8;
				end=end<<8|0xFF;
			}
		}
		cusize+=ncodes;
	}

	auto t2=__rdtsc();
	if(loud)
	{
		printf("AC decode:  %lld cycles\n", t2-t1);
	}
}