//ac.cpp - Arithmetic Coder implementation
//Arithmetic Coder from ZPAQ0
//modifications by Ayman Wagih Mohsen, unless source link provided.
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
#include<math.h>
#include<vector>

#ifdef __linux__
#define	__rdtsc	__builtin_ia32_rdtsc
#define	scanf_s	scanf
#else
#include<immintrin.h>
//#include<tmmintrin.h>
//#include<intrin.h>
#define	scanf	scanf_s
#endif

	#define	LOG_WINDOW_SIZE		16	//[2, 16]	do not change

	#define	PRINT_ERROR//should be enabled
	#define	MEASURE_PREDICTION//disable on benchmark/release

//choose one:
	#define	PREDICTOR_AWM_CONFBIT//best, but slow
//	#define	PREDICTOR_AWM_DIVFREE
//	#define	PREDICTOR_AWM_CONFSIZE
//	#define	PREDICTOR_ZPAQ0

//	#define	PREDICTOR_AWM_CONFBIT_WND//bad
//	#define PREDICTOR_AWM_PREVBIT//bad
//	#define	PREDICTOR_AWM_FULLCONF//bad


//disable all of these at release:
//	#define	PRINT_VISUAL
//	#define	PRINT_BITPLANES
//	#define	DEBUG_PREDICTOR
//	#define	DEBUG_SIMD_ENC
//	#define	DEBUG_SIMD_DEC
//	#define	HARD_AVX2_PROFILE
//	#define	DEBUG_PREDICTOR2
//	#define	DEBUG_PORTABILITY
//	#define	DEBUG_CONFBIT
//	#define	DEBUG_CONFBIT_WND
	#define DEBUG_DIVFREE

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
	scanf("%d", &x);
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

#ifdef DEBUG_PREDICTOR
struct			ProbInfo
{
	int kb, kc, acc, p0;
	ProbInfo():kb(0), kc(0), acc(0), p0(0){}
	ProbInfo(int kb, int kc, int acc, int p0):kb(kb), kc(kc), acc(acc), p0(p0){}
};
static std::vector<ProbInfo> probs;
#endif
#ifdef DEBUG_SIMD_ENC
struct			ProbInfo
{
	unsigned start, end, acc, p0;
	ProbInfo():start(0), end(0), acc(0), p0(0){}
	ProbInfo(int start, int end, int acc, int p0):start(start), end(end), acc(acc), p0(p0){}
};
static std::vector<ProbInfo> probs;
const int examined_plane=4;
#endif
#ifdef DEBUG_PREDICTOR2
const int examined_plane=6, examined_bit=7677;
#endif
inline int		weighted_bitsum(int x)
{
	int sum=0;
	for(int k=0;k<16;++k)
	{
		int bit=x>>k&1;
		sum+=(k+1)&-bit;
	}
	return sum;
}
inline int		hamming_weight(int x)
{
	x-=x>>1&0x55555555;
	x=(x&0x33333333)+(x>>2&0x33333333);
	return ((x+(x>>4))&0x0F0F0F0F)*0x01010101>>24;
}
void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud)
{
#ifdef PRINT_VISUAL
	ranges.clear();
	int vkc_start=0, vkc_end=0;
#endif
	if(!imsize)
		return;
	auto t1=__rdtsc();
#ifdef MEASURE_PREDICTION
	u64 hitnum=0, hitden=0;//prediction efficiency
#ifdef PREDICTOR_AWM_DIVFREE
	u64 hitnum2=0, hitden2=0;
#endif
#endif
	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		auto &plane=planes[depth-1-kp];
		plane.reserve(imsize>>8);
#ifdef PREDICTOR_AWM_CONFBIT
		int prob=1<<(LOG_WINDOW_SIZE-1), hitcount=1;//cheap weighted average predictor
#endif
#ifdef PREDICTOR_AWM_DIVFREE
		int prob=1<<(LOG_WINDOW_SIZE-1);
		int logden=floor_log2(imsize)+1, hitcount=1<<(logden-1);
	//	int logden=0, hitcount=1;
#endif
#if defined PREDICTOR_AWM_CONFSIZE || defined PREDICTOR_AWM_FULLCONF
		int prob=1<<(LOG_WINDOW_SIZE-1);//cheap weighted average predictor
#endif
#ifdef PREDICTOR_AWM_PREVBIT
		int prevbit=0;
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
		int prob=1<<(LOG_WINDOW_SIZE-1), hitcount=0xAAAAAAAA&((1<<LOG_WINDOW_SIZE)-1);//cheap weighted average predictor
#endif
#ifdef PREDICTOR_ZPAQ0
		int n0=1, n1=1;
#endif

		unsigned start=0, end=0xFFFFFFFF;
		for(int kb=0;kb<imsize;)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
		{
			int bit=buffer[kb]>>kp&1;
			//if(bit)//impossible non-causal predictor
			//	prob=65536/16;
			//else
			//	prob=65536*15/16;
			
			//if(kp==2&&kb==71831)
			//	int LOL_1=0;
			u64 range=end-start;
			if(range<3)
			//if(range==1)
			{
				plane.push_back(start>>24);
				plane.push_back(start>>16&0xFF);
				plane.push_back(start>>8&0xFF);
				plane.push_back(start&0xFF);
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
			
#ifdef PREDICTOR_AWM_CONFBIT
			//if(kp==0&&kb==1416)
		//	if(kp==2&&kb==982)
		//		int LOL_1=0;
			//int p0=prob;
			int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
			//printf("\rconfidence=%lf", (double)hitcount/(kb+1));
#ifdef DEBUG_CONFBIT
			if(kp==6&&kb>=6056900&&kb<6057000)
				printf("%d %08X~%08X %08X acc=%04X hit=%d p0=%04X %08X %d\n", kb, start, end, (unsigned)range, prob, hitcount, p0, middle, bit);
#endif
#ifdef DEBUG_PORTABILITY
			//if(kp==2&&kb>=1100&&kb<1300)
			if(kp==2&&kb>=700&&kb<1000)
			//if(kp==0&&kb>=1400&&kb<1600)
			//if(kp==0&&kb>=0&&kb<200)
			//if(kp==0&&kb>=700&&kb<800)
				printf("%d %d|%08X~%08X|%04X|%d/%d|%d|%08X\n", kb, bit, start, end, prob, hitcount, kb+1, p0, middle);
			//	printf("%d %d|%08X~%08X|acc %04X|%d/%d|p0=%d|mid=%08X\n", kb, bit, start, end, prob, hitcount, kb+1, p0, middle);
			//	printf("%d: %08X~%08X, acc=%04X, %d/%d, p0=%d, mid=%08X\n", kb, start, end, prob, hitcount, kb+1, p0, middle);
			//	printf("%d: %d/%d = %lf, p0 = %d\n", kb, hitcount, kb+1, (double)hitcount/(kb+1), p0);
#endif
#endif
#ifdef PREDICTOR_AWM_DIVFREE
			int p0=0x8000+((long long)(prob-0x8000)*hitcount>>logden);
//#ifdef DEBUG_DIVFREE
//			if(p0<1||p0>0xFFFF)
//			{
//				printf("kp %d, kb %d, p0 %04X\n", kp, kb, p0);
//				int x=0;
//				scanf_s("%d", &x);
//			}
//#endif
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
			int p0=0x8000+((long long)(prob-0x8000)*hitcount>>LOG_WINDOW_SIZE);
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);

			middle+=(middle==start)-(middle==end);
			//middle-=middle==end;//half-open
			
#ifdef DEBUG_CONFBIT_WND
			//if(kp==0&&kb>700&&kb<800)
			//if(kp==2&&kb>71700&&kb<71900)
			//if(kp==4&&kb>8100&&kb<8300)
			//if(kp==3&&kb>8800&&kb<9000)
			if(kp==7&&kb>2800&&kb<3000)
			{
				if(kb==8880)
					int LOL_1=0;
				printf("%d %08X~%08X %08X acc=%04X hit=%04X p0=%04X %08X %d\n", kb, start, end, (unsigned)range, prob, hitcount, p0, middle, bit);
			}
#endif
#endif
#ifdef PREDICTOR_AWM_FULLCONF
			int p0=clamp(1, prob, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#ifdef DEBUG_PREDICTOR2
			if(kp==examined_plane&&kb==examined_bit)
				printf("%08X~%08X, acc=%04X, p0=%d, middle=%08X, bit=%d\n", start, end, prob, p0, middle, bit);
#endif
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			int p0=1+(0xFFFD&-!prevbit);
			//p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#endif
			middle+=(middle==start)-(middle==end);
#ifdef PREDICTOR_AWM_CONFSIZE
			//int p0=prob;
			//if(conf_den)
			//	p0=0x8000+((long long)(p0-0x8000)*kb*conf_invden>>16);
#if 1
			int p0=prob;		//1.471149

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*hamming_weight(p_alt)>>6);		//1.450343

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*p_alt>>16);						//1.361238
			//printf("new=%d, acc=%04X, alt=%04X, P(0)=%d=%lf\n", bit, prob, ((prob^prob<<1)&0xFFFF), p0, (double)p0/0x10000);

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*weighted_bitsum(p_alt)>>7);		//1.240445

			int conf_den=(plane.size()<<3)+kb;
#ifdef DEBUG_SIMD_ENC
			if(kp==examined_plane&&kb==40)
			//if(kp==examined_plane&&kb==7)
				int LOL_1=0;
#endif
			if(conf_den)
			//if(plane.size())
			{
				p0=0x8000+(long long)(p0-0x8000)*kb/conf_den;
				//double conf=(double)kb/conf_den;
				//p0=0x8000+(int)floor(((int)p0-0x8000)*conf);
				//p0=(unsigned)(0x8000+((int)p0-0x8000)*conf);
				//p0=(unsigned)(0x7FFF+((int)p0-0x7FFF)*conf);

			//	double ratio=(double)kb/(plane.size()<<3);
			//	p0=0x7FFF+(int)(((double)p0-0x7FFF)*ratio/(1+ratio));
#ifdef DEBUG_PREDICTOR
				if(kp==0)
					printf("r=%lf, acc=%04X, p0=%d\n", ratio, prob, p0);
#endif
				//if(!(kb%10000))
				//	printf("kp=%d, kb=%d, psize=%d, r=%lf, c=%lf, acc=%04X, P(0)=%d\n", kp, kb, plane.size(), ratio, ratio/(1+ratio), prob, p0);

				//int ratio=(kb<<3)/plane.size();
				//p0=0x7FFF+(int)((u64)((p0-0x7FFF)*ratio)/(1+ratio));//taking confidence into account
			}
			//else
			//	p0=0x8000;
#endif
			p0=clamp(1, p0, prob_max);
#ifdef DEBUG_PREDICTOR
			if(kp==0)
				probs.push_back(ProbInfo(kb, plane.size(), prob, p0));
#endif
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#endif
#ifdef PREDICTOR_ZPAQ0
			int sh=(n0>=0x10000)|(n1>=0x10000);
			n0>>=sh, n1>>=sh;
			n0+=n0<1, n1+=n1<1;
			unsigned middle=start+(unsigned)(range*n0/(n0+n1));
#endif
#ifdef DEBUG_SIMD_ENC
			if(kp==examined_plane&&kb==40)
			//if(kp==examined_plane&&kb==7)
			//if(kp==1&&(kb==904||kb==905))//
			//if(kp==5&&kb==1)//
				printf("kp=%d, kb=%d: %08X~%08X, mid=%08X, p0=%d\n", kp, kb, start, end, middle, p0);//
			if(kp==examined_plane)
				probs.push_back(ProbInfo(start, end, prob, p0));
#endif

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
			
#ifdef PREDICTOR_AWM_CONFBIT
			hitcount+=bit^(p0>=0x8000);
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_DIVFREE
		//	int sh=kb+1>=(1<<logden);
		//	hitcount=hitcount<<sh;//|sh;
		//	logden+=sh;
			int correct=bit^(p0>=0x8000);
			hitcount+=correct-!correct;
			hitcount=clamp(1, hitcount, (1<<logden)-1);
		//	hitcount-=!correct;
		//	hitcount+=hitcount==0;
			
#ifdef DEBUG_DIVFREE
		//	printf("\r%6d: %6d / %6d = %lf%%\t", kb, hitcount, 1<<logden, 100.*hitcount/(1<<logden));
			if(hitcount<=0)
			{
				printf("kp %d, kb %d, hit %d\n", kp, kb, hitcount);
				int x=0;
				scanf_s("%d", &x);
			}
#endif
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#if defined PREDICTOR_AWM_CONFSIZE || defined PREDICTOR_AWM_FULLCONF
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
			hitcount=(bit^(p0>=0x8000))<<15|hitcount>>1;
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			prevbit=bit;
#endif
#ifdef PREDICTOR_ZPAQ0
			n0+=!bit, n1+=bit;
#endif
#ifdef MEASURE_PREDICTION
			hitnum+=bit^(p0>=0x8000), ++hitden;
#endif

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
			//	start=middle;
				start=middle+1;
			else
			//	end=middle;
				end=middle-1;
			//if(start>end)
			//	int LOL_1=0;
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
#if defined MEASURE_PREDICTION && defined PREDICTOR_AWM_DIVFREE
		hitnum2+=hitcount, hitden2+=1LL<<logden;
#endif
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
#ifdef MEASURE_PREDICTION
		printf("Predicted (actual): %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#ifdef PREDICTOR_AWM_DIVFREE
		printf("Predicted (used):   %6lld / %6lld = %lf%%\n", hitnum2, hitden2, 100.*hitnum2/hitden2);
#endif
#endif
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
		
#ifdef PREDICTOR_AWM_CONFBIT
		int prob=1<<(LOG_WINDOW_SIZE-1), hitcount=1;//cheap weighted average predictor
#endif
#ifdef PREDICTOR_AWM_DIVFREE
		int prob=1<<(LOG_WINDOW_SIZE-1);
		int logden=floor_log2(imsize)+1, hitcount=1<<(logden-1);
	//	int logden=0, hitcount=1;
#endif
#if defined PREDICTOR_AWM_CONFSIZE || defined PREDICTOR_AWM_FULLCONF
		int prob=1<<(LOG_WINDOW_SIZE-1);//cheap weighted average predictor
#endif
#ifdef PREDICTOR_AWM_PREVBIT
		int prevbit=0;
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
		int prob=1<<(LOG_WINDOW_SIZE-1), hitcount=0xAAAAAAAA&((1<<LOG_WINDOW_SIZE)-1);//cheap weighted average predictor
#endif
#ifdef PREDICTOR_ZPAQ0
		int n0=1, n1=1;
#endif

		unsigned start=0, end=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane);
		for(int kc=4, kb=0;kb<imsize;)//bit-pixel loop
		{
			u64 range=end-start;
			if(range<3)
			//if(range==1)
			{
				code=load32_big((unsigned char*)plane+kc);
				kc+=4;
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*weighted_bitsum(p_alt)>>7);

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*hamming_weight(p_alt)>>6);

			//int p_alt=(prob^prob<<1)&0xFFFF;
			//int p0=prob+((~prob-prob)*p_alt>>16);
			
#ifdef PREDICTOR_AWM_CONFBIT
			//int p0=prob;
			int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#ifdef DEBUG_CONFBIT
			if(kp==6&&kb>=6056900&&kb<6057000)
				printf("%d %08X~%08X %08X acc=%04X hit=%08X p0=%04X %08X %08X %d\n", kb, start, end, (unsigned)range, prob, hitcount, p0, middle, code, code>=middle);
#endif
#endif
#ifdef PREDICTOR_AWM_DIVFREE
			int p0=0x8000+((long long)(prob-0x8000)*hitcount>>logden);
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
			int p0=0x8000+((long long)(prob-0x8000)*hitcount>>LOG_WINDOW_SIZE);
			p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);

			//middle-=middle==end;//half-open
#ifdef DEBUG_CONFBIT_WND
			//if(kp==0&&kb>700&&kb<800)
			//if(kp==2&&kb>71700&&kb<71900)
			//if(kp==4&&kb>8100&&kb<8300)
			//if(kp==3&&kb>8800&&kb<9000)
			if(kp==7&&kb>2800&&kb<3000)
			{
				//if(kb==71821)
				//if(kb==8217)
				//if(kb==8880)
				if(kb==2915)
					int LOL_1=0;
				printf("%d %08X~%08X %08X acc=%04X hit=%04X p0=%04X %08X %08X %d\n", kb, start, end, (unsigned)range, prob, hitcount, p0, middle, code, code>=middle);
			}
#endif
#endif
#ifdef PREDICTOR_AWM_FULLCONF
			int p0=clamp(1, prob, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#ifdef DEBUG_PREDICTOR2
			if(kp==examined_plane&&kb==examined_bit)
				printf("%08X~%08X, acc=%04X, p0=%d, middle=%08X, code=%08X, bit=%d\n", start, end, prob, p0, middle, code, code>=middle);
#endif
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			int p0=1+(0xFFFD&-!prevbit);
			//p0=clamp(1, p0, prob_max);
			unsigned middle=start+(unsigned)(range*p0>>LOG_WINDOW_SIZE);
#endif
			middle+=(middle==start)-(middle==end);
#ifdef PREDICTOR_AWM_CONFSIZE
			int p0=prob;

			int conf_den=((kc-4)<<3)+kb;
			if(conf_den)
			//if(kc-4>0)
			{
				p0=0x8000+(long long)(p0-0x8000)*kb/conf_den;
				//double conf=(double)kb/conf_den;
				//p0=0x8000+(int)floor(((int)p0-0x8000)*conf);
				//p0=(unsigned)(0x8000+((int)p0-0x8000)*conf);
				//p0=(unsigned)(0x7FFF+((int)p0-0x7FFF)*conf);

			//	double ratio=(double)kb/((kc-4)<<3);
			//	p0=0x7FFF+(int)(((double)p0-0x7FFF)*ratio/(1+ratio));
#ifdef DEBUG_PREDICTOR
				//if(kp==0)
				//	printf("r=%lf, acc=%04X, p0=%d\n", ratio, prob, p0);
#endif
			}
			//else
			//	p0=0x8000;
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
#endif
#ifdef PREDICTOR_ZPAQ0
			int sh=(n0>=0x10000)|(n1>=0x10000);
			n0>>=sh, n1>>=sh;
			n0+=n0<1, n1+=n1<1;
			unsigned middle=start+(unsigned)(range*n0/(n0+n1));
#endif
			int bit=code>middle;
		//	int bit=code>=middle;
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
			
#ifdef PREDICTOR_AWM_CONFBIT
			hitcount+=bit^(p0>=0x8000);
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_DIVFREE
		//	int sh=kb+1>=(1<<logden);
		//	hitcount=hitcount<<sh;//|sh;
		//	logden+=sh;
			int correct=bit^(p0>=0x8000);
			hitcount+=correct-!correct;
			hitcount=clamp(1, hitcount, (1<<logden)-1);
		//	hitcount-=!correct;
		//	hitcount+=hitcount==0;
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#if defined PREDICTOR_AWM_CONFSIZE || defined PREDICTOR_AWM_FULLCONF
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_CONFBIT_WND
			hitcount=(bit^(p0>=0x8000))<<15|hitcount>>1;
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			prevbit=bit;
#endif
#ifdef PREDICTOR_ZPAQ0
			n0+=!bit, n1+=bit;
#endif

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
			//	start=middle;
				start=middle+1;
			else
			//	end=middle;
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

#if !defined __GNUC__ && !defined NO_SIMD
#ifdef DEBUG_SIMD_DEC
const int		examined_plane=0, examined_pixel=776;
//const int		examined_plane=0, examined_pixel=772;
//const int		examined_plane=0, examined_pixel=8;
//const int		examined_plane=4, examined_pixel=1227515;
int examine_flag=0;
#endif
inline __m128i	mul_sr16(__m128i const &a32, __m128i const &b16)
{
	__m128i lo=_mm_mul_epu32(a32, b16);
	lo=_mm_srli_epi64(lo, 16);
	__m128i a2=_mm_shuffle_epi32(a32, _MM_SHUFFLE(2, 3, 0, 1));
	__m128i b2=_mm_shuffle_epi32(b16, _MM_SHUFFLE(2, 3, 0, 1));
	__m128i hi=_mm_mul_epu32(a2, b2);
	hi=_mm_srli_epi64(hi, 16);
	__m128i product=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(lo), _mm_castsi128_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
	product=_mm_shuffle_epi32(product, _MM_SHUFFLE(3, 1, 2, 0));
	return product;
}
void			abac_encode_sse2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();

	std::vector<std::string> planes(depth);
	for(int k=0;k<depth;++k)
		planes[k].reserve(imsize>>8);

	__m128i bitmask=_mm_set_epi32(8, 4, 2, 1);
	__m128i m_ones=_mm_set1_epi32(-1);
	__m128i probhalf=_mm_set1_epi32(0x8000);
	__m128i msbyte_limit=_mm_set1_epi32(0x1000000);
	__m128i m_one=_mm_set1_epi32(1), m_probmax=_mm_set1_epi32(prob_max);
	__m128i mask_lo16=_mm_set1_epi32(0x0000FFFF);
	__m128 f_one=_mm_set1_ps(1);
	__m128 f_half=_mm_set1_ps(0.5f), f_minushalf=_mm_set1_ps(-0.5f);
	for(int kp=0;kp+3<depth;kp+=4)//bit-plane loop: encode 4 planes simultaneously
	{
		__m128i i_outbitcounts=_mm_setzero_si128();
		__m128 f_outbitcounts=_mm_cvtepi32_ps(i_outbitcounts);
		__m128i prob=_mm_set1_epi32(1<<(LOG_WINDOW_SIZE-1));//estimated P(next bit = 0): weighted sum (sum i=1 to 16: !bit[-i]*2^i)
		__m128i hitcount=_mm_set1_epi32(0xAAAAAAAA&((1<<LOG_WINDOW_SIZE)-1));
		//__m128i hitcount=_mm_set1_epi32(1);
#ifdef PREDICTOR_AWM_PREVBIT
		const __m128i pb_base=_mm_set1_epi32(2), pb_offset=_mm_set1_epi32(0xFFFF-2*2);
		__m128i prevbit=_mm_setzero_si128();
#endif
		__m128i m_kb=_mm_set1_epi32(1);
		__m128i start=_mm_setzero_si128(), end=_mm_set1_epi32(-1);//inclusive range
		for(int kb=0;kb<imsize;)
		{
			__m128i range=_mm_sub_epi32(end, start);
			__m128i eq_one=_mm_cmpeq_epi32(range, _mm_set1_epi32(1));
			int condition_eq_one=_mm_movemask_ps(_mm_castsi128_ps(eq_one));
			if(condition_eq_one)//because 1=0.9999...
			{
				for(int kch=0;kch<4;++kch)
				{
					if(condition_eq_one>>kch&1)
					{
						unsigned istart=start.m128i_u32[kch];
						auto &plane=planes[depth-1-(kp+kch)];
						plane.push_back(istart>>24);
						plane.push_back(istart>>16&0xFF);
						plane.push_back(istart>>8&0xFF);
						plane.push_back(istart&0xFF);
						start.m128i_u32[kch]=0, end.m128i_u32[kch]=0xFFFFFFFF;
					}
				}
				range=_mm_sub_epi32(end, start);
			}

		/*	//calculate p0=0x8000+(long long)(p0-0x8000)*kb/conf_den;
			__m128i p0=_mm_sub_epi32(prob, probhalf);
			__m128i lo=_mm_mul_epi32(p0, m_kb);//SSE4.1
			__m128i hi=_mm_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			hi=_mm_mul_epi32(hi, m_kb);//SSE4.1
			__m128i den=_mm_add_epi32(i_outbitcounts, m_kb);
			if(den.m128i_u32[0])
				lo.m128i_i64[0]/=den.m128i_u32[0];
			if(den.m128i_u32[1])
				hi.m128i_i64[0]/=den.m128i_u32[1];
			if(den.m128i_u32[2])
				lo.m128i_i64[1]/=den.m128i_u32[2];
			if(den.m128i_u32[3])
				hi.m128i_i64[1]/=den.m128i_u32[3];
			p0=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(lo), _mm_castsi128_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			p0=_mm_shuffle_epi32(p0, _MM_SHUFFLE(3, 1, 2, 0));
			p0=_mm_add_epi32(p0, probhalf);

			//if confidence denominator==0 then confidence=100%
			den=_mm_cmpeq_epi32(den, _mm_setzero_si128());
			__m128i mask=_mm_xor_si128(den, m_ones);
			p0=_mm_and_si128(p0, mask);
			den=_mm_and_si128(den, prob);
			p0=_mm_or_si128(p0, den);//*/
		
#ifdef PREDICTOR_AWM_CONFBIT_WND
			//calculate p0 = 0x8000 + ((acc-0x8000)*hitcount>>16)
			__m128i p0=_mm_sub_epi32(prob, probhalf);
			p0=_mm_mulhi_epu16(p0, hitcount);
			p0=_mm_and_si128(p0, mask_lo16);
			p0=_mm_add_epi32(p0, probhalf);
#endif
#ifdef PREDICTOR_AWM_CONFBIT
			//calculate confidence = (hitcount+1)/(kb+1)
			__m128 conf_num=_mm_cvtepi32_ps(hitcount);
			__m128 conf_den=_mm_cvtepi32_ps(m_kb);
			__m128 conf=_mm_div_ps(conf_num, conf_den);
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			__m128i p0=_mm_and_si128(pb_offset, prevbit);
			p0=_mm_add_epi32(p0, pb_base);
#endif
#ifdef PREDICTOR_AWM_CONFSIZE
			//calculate confidence = kb / (8kc+kb)		where kb: in bit idx, kc: out byte idx
			__m128 conf_num=_mm_set1_ps((float)kb);
			__m128 conf_den=_mm_add_ps(conf_num, f_outbitcounts);
			__m128 conf=_mm_div_ps(conf_num, conf_den);

			//if confidence denominator==0 then confidence=100%
			__m128 den_nz=_mm_cmpneq_ps(conf_den, _mm_setzero_ps());
			conf=_mm_and_ps(conf, den_nz);
			den_nz=_mm_xor_ps(den_nz, _mm_castsi128_ps(m_ones));
			den_nz=_mm_and_ps(den_nz, f_one);
			conf=_mm_or_ps(conf, den_nz);
#endif
#if !defined PREDICTOR_AWM_CONFBIT_WND && !defined PREDICTOR_AWM_PREVBIT
			//predict P(0) = 1/2 + (W(0)-1/2)*confidence
			__m128i p0=_mm_sub_epi32(prob, probhalf);
			__m128 f_p0=_mm_cvtepi32_ps(p0);
			f_p0=_mm_mul_ps(f_p0, conf);
			
			f_p0=_mm_round_ps(f_p0, _MM_FROUND_TRUNC);//SSE4.1
			//__m128 f_mask=_mm_cmplt_ps(f_p0, _mm_setzero_ps());//correction to cast product to integer
			//__m128 f_correction=_mm_and_ps(f_mask, f_half);
			//f_mask=_mm_xor_ps(f_mask, _mm_castsi128_ps(m_ones));
			//f_mask=_mm_and_ps(f_mask, f_minushalf);
			//f_correction=_mm_or_ps(f_correction, f_mask);
			//f_p0=_mm_add_ps(f_p0, f_correction);

			p0=_mm_cvtps_epi32(f_p0);			//FIX THIS
			p0=_mm_add_epi32(p0, probhalf);//*/
#endif

			//clamp probability: P(0) = clamp(2^-L, P(0), 1-2^-L)	where L=16		//probability shouldn't be exactly 0 or 100%
			__m128i c_min=_mm_cmpgt_epi32(p0, m_one);
			__m128i c_max=_mm_cmplt_epi32(p0, m_probmax);
			p0=_mm_and_si128(p0, c_min);
			p0=_mm_and_si128(p0, c_max);
			c_min=_mm_xor_si128(c_min, m_ones);
			c_max=_mm_xor_si128(c_max, m_ones);
			c_min=_mm_and_si128(c_min, m_one);
			c_max=_mm_and_si128(c_max, m_probmax);
			p0=_mm_or_si128(p0, c_min);
			p0=_mm_or_si128(p0, c_max);

#ifdef DEBUG_SIMD_ENC
			if(examined_plane>=kp&&examined_plane<kp+4)
			{
				auto &p=probs[kb];
				if(p0.m128i_u32[0]!=p.p0)
					int LOL_1=0;
			}
			//if(kp==0&&kb==5)
			//if(1>=kp&&1<kp+4&&(kb==904||kb==905))//
			//if(5>=kp&&5<kp+4&&kb==1)//
			//	int LOL_1=0;//
#endif

			//calculate middle = start+((u64)(end-start)*p0>>L)
			__m128i middle=mul_sr16(range, p0);
		/*	__m128i lo=_mm_mul_epu32(range, p0);
			lo=_mm_srli_epi64(lo, 16);
			range=_mm_shuffle_epi32(range, _MM_SHUFFLE(2, 3, 0, 1));
			eq_one=_mm_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			__m128i hi=_mm_mul_epu32(range, eq_one);
			hi=_mm_srli_epi64(hi, 16);
			__m128i middle=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(lo), _mm_castsi128_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			middle=_mm_shuffle_epi32(middle, _MM_SHUFFLE(3, 1, 2, 0));//*/
			middle=_mm_add_epi32(middle, start);

			//get next bit
			__m128i bit_c=_mm_set1_epi32(buffer[kb]);
			bit_c=_mm_srli_epi32(bit_c, kp);
			bit_c=_mm_and_si128(bit_c, bitmask);
			bit_c=_mm_cmpeq_epi32(bit_c, _mm_setzero_si128());
			__m128i bit=_mm_xor_si128(bit_c, m_ones);
			
#ifdef PREDICTOR_AWM_PREVBIT
			prevbit=bit_c;
#endif
#ifdef DEBUG_SIMD_DEC
			if(examined_plane>=kp&&examined_plane<kp+4&&kb==examined_pixel)
			//if(5>=kp&&5<kp+4&&kb==2326)
			//if(2>=kp&&2<kp+4&&kb==1283)
			{
				printf("buffer[%d]=%04X, p0=%d, %08X~%08X, mid=%08X\n", kb, buffer[kb], p0.m128i_i32[examined_plane&3], start.m128i_u32[examined_plane&3], end.m128i_u32[examined_plane&3], middle.m128i_u32[examined_plane&3]);
				examine_flag=1;
			}
			//	int LOL_1=0;
#endif

#ifdef PREDICTOR_AWM_CONFBIT_WND
			//update hitcount: hit = !bit != (p0<0x8000),	hitcount = hit<<15|hitcount>>1
			__m128i b2=_mm_cmplt_epi32(p0, probhalf);
			b2=_mm_xor_si128(b2, bit_c);
			b2=_mm_and_si128(b2, probhalf);
			hitcount=_mm_srli_epi32(hitcount, 1);
			hitcount=_mm_or_si128(hitcount, b2);
#endif
#ifdef PREDICTOR_AWM_CONFBIT
			//update hitcount = !bit != (p0<0x8000)
			__m128i b2=_mm_cmplt_epi32(p0, probhalf);
			b2=_mm_xor_si128(b2, bit_c);
			hitcount=_mm_sub_epi32(hitcount, b2);
#endif

			//update weighted sum W(0) = !bit_new | W_prev(0)>>1	 = (sum i=1 to 16: !bit[-i]*2^i) / 2^16
			prob=_mm_srli_epi32(prob, 1);
			__m128i bit15=_mm_and_si128(bit_c, probhalf);
			prob=_mm_or_si128(prob, bit15);

			//select range: if(bit)start=middle; else end=middle-1;		same as: start=bit?middle:start, end=bit?end:middle-1;
			__m128i midm1=_mm_sub_epi32(middle, m_one);
			__m128i midp1=_mm_add_epi32(middle, m_one);
			//__m128i midm1=_mm_sub_epi32(middle, _mm_set1_epi32(2));//X
			
			midp1=_mm_and_si128(midp1, bit);
			start=_mm_and_si128(start, bit_c);
			start=_mm_or_si128(start, midp1);
			//middle=_mm_and_si128(middle, bit);
			//start=_mm_and_si128(start, bit_c);
			//start=_mm_or_si128(start, middle);

			midm1=_mm_and_si128(midm1, bit_c);
			end=_mm_and_si128(end, bit);
			end=_mm_or_si128(end, midm1);

			++kb;
			m_kb=_mm_add_epi32(m_kb, m_one);

			//check if most-significant byte is ready
			__m128i diff=_mm_xor_si128(start, end);
			diff=_mm_srli_epi32(diff, 24);
			diff=_mm_cmpeq_epi32(diff, _mm_setzero_si128());
			int ready_mask=_mm_movemask_ps(_mm_castsi128_ps(diff));
			if(ready_mask)
			{
				for(int kch=0;kch<4;++kch)
				{
					if(ready_mask>>kch&1)
					{
						unsigned istart=start.m128i_u32[kch], iend=end.m128i_u32[kch];
						auto &plane=planes[depth-1-(kp+kch)];
						do
						{
							plane.push_back(istart>>24);
#ifdef DEBUG_SIMD_DEC
							if(examine_flag&&kp+kch==examined_plane)
								printf("Shift-out: %02X\n", istart>>24);
#endif
							istart<<=8, iend=iend<<8|0xFF;
							i_outbitcounts.m128i_u32[kch]+=8;
							f_outbitcounts=_mm_cvtepi32_ps(i_outbitcounts);
						}while((istart^iend)<0x1000000);
						start.m128i_u32[kch]=istart, end.m128i_u32[kch]=iend;
#ifdef DEBUG_SIMD_DEC
						if(examine_flag&&kp+kch==examined_plane)
						{
							printf("New bounds: %08X~%08X\n", start.m128i_u32[kch], end.m128i_u32[kch]);
							examine_flag=0;
						}
#endif
					}
				}
			}
		}
		for(int kch=0;kch<4;++kch)
		{
			auto &plane=planes[depth-1-(kp+kch)];
			unsigned istart=start.m128i_u32[kch];
			plane.push_back(istart>>24&0xFF);//big-endian
			plane.push_back(istart>>16&0xFF);
			plane.push_back(istart>>8&0xFF);
			plane.push_back(istart&0xFF);
		}
	}
	auto t_enc=__rdtsc();
	//out_sizes.resize(depth);
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
	out_data.clear();
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)out_data.size()<<3;
		printf("AC encode SSE2:  %lld cycles (Enc: %lld cycles)\n", t2-t1, t_enc-t1);
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
void			abac_decode_sse2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));

	int *plane_offsets=new int[depth];
	for(int kp=depth-1, cusize=0;kp>=0;--kp)
	{
		plane_offsets[kp]=cusize;
		cusize+=sizes[depth-1-kp];
	}
	
	__m128i bitmask=_mm_set_epi32(8, 4, 2, 1);
	__m128i m_ones=_mm_set1_epi32(-1);
	__m128i probhalf=_mm_set1_epi32(0x8000);
	__m128i msbyte_limit=_mm_set1_epi32(0x1000000);
	__m128i m_one=_mm_set1_epi32(1), m_probmax=_mm_set1_epi32(prob_max);
//	__m128i shuffle_bigendian=_mm_set_epi8(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3);
	__m128i mask_lo16=_mm_set1_epi32(0x0000FFFF);
	__m128 f_one=_mm_set1_ps(1), f_half=_mm_set1_ps(0.5f), f_minushalf=_mm_set1_ps(-0.5f);
	for(int kp=0, cusize=0;kp+3<depth;kp+=4)//bit-plane loop
	{
		__m128i i_outbitcounts=_mm_setzero_si128();
		__m128 f_outbitcounts=_mm_cvtepi32_ps(i_outbitcounts);
		__m128i prob=_mm_set1_epi32(1<<(LOG_WINDOW_SIZE-1));//estimated P(next bit = 0): weighted sum (sum i=1 to 16: !bit[-i]*2^i)
		__m128i hitcount=_mm_set1_epi32(0xAAAAAAAA&((1<<LOG_WINDOW_SIZE)-1));
		//__m128i hitcount=_mm_set1_epi32(1);
#ifdef PREDICTOR_AWM_PREVBIT
		const __m128i pb_base=_mm_set1_epi32(2), pb_offset=_mm_set1_epi32(0xFFFF-2*2);
		__m128i prevbit=_mm_setzero_si128();
#endif
		__m128i m_kb=_mm_set1_epi32(1);
		__m128i start=_mm_setzero_si128(), end=_mm_set1_epi32(-1);//inclusive range
		__m128i code=_mm_set_epi32(
			load32_big((unsigned char*)data+plane_offsets[kp+3]),
			load32_big((unsigned char*)data+plane_offsets[kp+2]),
			load32_big((unsigned char*)data+plane_offsets[kp+1]),
			load32_big((unsigned char*)data+plane_offsets[kp]));
	//	code=_mm_shuffle_epi8(code, shuffle_bigendian);
		__m128i kc=_mm_set1_epi32(4);
		for(int kb=0;kb<imsize;)//bit-pixel loop
		{
			__m128i range=_mm_sub_epi32(end, start);
			__m128i eq_one=_mm_cmpeq_epi32(range, _mm_set1_epi32(1));
			int condition_eq_one=_mm_movemask_ps(_mm_castsi128_ps(eq_one));
			if(condition_eq_one)//because 1=0.9999...
			{
				for(int kch=0;kch<4;++kch)
				{
					if(condition_eq_one>>kch&1)
					{
						auto plane=data+plane_offsets[kp+kch];
						code.m128i_u32[kch]=load32_big((unsigned char*)plane+kc.m128i_i32[kch]);
						kc.m128i_i32[kch]+=4;
						start.m128i_u32[kch]=0, end.m128i_u32[kch]=0xFFFFFFFF;
					}
				}
				range=_mm_sub_epi32(end, start);
			}
			
#ifdef PREDICTOR_AWM_CONFBIT_WND
			//calculate p0 = 0x8000 + ((acc-0x8000)*hitcount>>16)
			__m128i p0=_mm_sub_epi32(prob, probhalf);
			p0=_mm_mulhi_epu16(p0, hitcount);
			p0=_mm_and_si128(p0, mask_lo16);
			p0=_mm_add_epi32(p0, probhalf);
#endif
#ifdef PREDICTOR_AWM_CONFBIT
			//if(kp==0&&kb==1416)
			if(kp==0&&kb==982)
				int LOL_1=0;
			//calculate confidence = (hitcount+1)/(kb+1)
			__m128 conf_num=_mm_cvtepi32_ps(hitcount);
			__m128 conf_den=_mm_cvtepi32_ps(m_kb);
			__m128 conf=_mm_div_ps(conf_num, conf_den);
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			__m128i p0=_mm_and_si128(pb_offset, prevbit);
			p0=_mm_add_epi32(p0, pb_base);
#endif
#ifdef PREDICTOR_AWM_CONFSIZE
			//calculate confidence = kb / (8kc+kb)		where kb: in bit idx, kc: out byte idx
			__m128 conf_num=_mm_set1_ps((float)kb);
			__m128 conf_den=_mm_add_ps(conf_num, f_outbitcounts);
			__m128 conf=_mm_div_ps(conf_num, conf_den);

			//if confidence denominator==0 then confidence=100%
			__m128 den_nz=_mm_cmpneq_ps(conf_den, _mm_setzero_ps());
			conf=_mm_and_ps(conf, den_nz);
			den_nz=_mm_xor_ps(den_nz, _mm_castsi128_ps(m_ones));
			den_nz=_mm_and_ps(den_nz, f_one);
			conf=_mm_or_ps(conf, den_nz);
#endif
#if !defined PREDICTOR_AWM_CONFBIT_WND && !defined PREDICTOR_AWM_PREVBIT
			//predict P(0) = 1/2 + (W(0)-1/2)*confidence
			__m128i p0=_mm_sub_epi32(prob, probhalf);
			__m128 f_p0=_mm_cvtepi32_ps(p0);
			f_p0=_mm_mul_ps(f_p0, conf);
			
			f_p0=_mm_round_ps(f_p0, _MM_FROUND_TRUNC);//SSE4.1
			//__m128 f_mask=_mm_cmplt_ps(f_p0, _mm_setzero_ps());//correction to cast product to integer
			//__m128 f_correction=_mm_and_ps(f_mask, f_half);
			//f_mask=_mm_xor_ps(f_mask, _mm_castsi128_ps(m_ones));
			//f_mask=_mm_and_ps(f_mask, f_minushalf);
			//f_correction=_mm_or_ps(f_correction, f_mask);
			//f_p0=_mm_add_ps(f_p0, f_correction);

			p0=_mm_cvtps_epi32(f_p0);
			p0=_mm_add_epi32(probhalf, p0);//*/
#endif

			//clamp probability: P(0) = clamp(2^-L, P(0), 1-2^-L)	where L=16		//probability shouldn't be exactly 0 or 100%
			__m128i c_min=_mm_cmpgt_epi32(p0, m_one);
			__m128i c_max=_mm_cmplt_epi32(p0, m_probmax);
			p0=_mm_and_si128(p0, c_min);
			p0=_mm_and_si128(p0, c_max);
			c_min=_mm_xor_si128(c_min, m_ones);
			c_max=_mm_xor_si128(c_max, m_ones);
			c_min=_mm_and_si128(c_min, m_one);
			c_max=_mm_and_si128(c_max, m_probmax);
			p0=_mm_or_si128(p0, c_min);
			p0=_mm_or_si128(p0, c_max);
			
			//calculate middle = start+((u64)(end-start)*p0>>L)
			__m128i middle=mul_sr16(range, p0);
		/*	__m128i lo=_mm_mul_epu32(range, p0);
			lo=_mm_srli_epi64(lo, 16);
			range=_mm_shuffle_epi32(range, _MM_SHUFFLE(2, 3, 0, 1));
			eq_one=_mm_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			__m128i hi=_mm_mul_epu32(range, eq_one);
			hi=_mm_srli_epi64(hi, 16);
			__m128i middle=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(lo), _mm_castsi128_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			middle=_mm_shuffle_epi32(middle, _MM_SHUFFLE(3, 1, 2, 0));//*/
			middle=_mm_add_epi32(middle, start);

			//unpack next bit	bit=code>=middle	unsigned compare
			__m128i bit=_mm_max_epu32(code, middle);//SSE4.1
			bit=_mm_cmpeq_epi32(bit, code);
			__m128i bit_c=_mm_xor_si128(bit, m_ones);
			//__m128i diff=_mm_sub_epi32(code, middle);
			//__m128i bit_c=_mm_cmplt_epi32(diff, _mm_setzero_si128());
			//__m128i bit=_mm_xor_si128(bit_c, m_ones);
#ifdef DEBUG_PORTABILITY
			//if(kp==0&&kb>=1100&&kb<1300)
			if(kp==0&&kb>=700&&kb<1000)
				printf("%d %d|%08X~%08X|%04X|%d/%d=%f|%d|%08X c%08X\n", kb, bit.m128i_i32[2]==-1, start.m128i_u32[2], end.m128i_u32[2], prob.m128i_u32[2], hitcount.m128i_i32[2], kb+1, conf.m128_f32[2], p0.m128i_i32[2], middle.m128i_u32[2], code.m128i_u32[2]);
			//if(kp==0&&kb>=1400&&kb<1600)
			//if(kp==0&&kb>=0&&kb<200)
			//if(kp==0&&kb>=700&&kb<800)
			//	printf("%d %d|%08X~%08X|%04X|%d/%d=%f|%d|%08X c%08X\n", kb, bit.m128i_i32[0]==-1, start.m128i_u32[0], end.m128i_u32[0], prob.m128i_u32[0], hitcount.m128i_i32[0], kb+1, conf.m128_f32[0], p0.m128i_i32[0], middle.m128i_u32[0], code.m128i_u32[0]);
			//	printf("%d: %08X~%08X, acc=%04X, %d/%d=%f, p0=%d, mid=%08X\n", kb, start.m128i_u32[0], end.m128i_u32[0], prob.m128i_u32[0], hitcount.m128i_i32[0], kb+1, conf.m128_f32[0], p0.m128i_i32[0], middle.m128i_u32[0]);
			//	printf("%d: %d/%d = %lf = %f, p0 = %d\n", kb, hitcount.m128i_i32[0], kb+1, (double)hitcount.m128i_i32[0]/(kb+1), conf.m128_f32[0], p0.m128i_i32[0]);
#endif
#ifdef PREDICTOR_AWM_PREVBIT
			prevbit=bit_c;
#endif

#ifdef DEBUG_SIMD_DEC
			if(examined_plane>=kp&&examined_plane<kp+4&&kb==examined_pixel)
			//if(5>=kp&&5<kp+4&&kb==2326)
			//if(2>=kp&&2<kp+4&&kb==1283)
				printf("buffer[%d]=%04X, p0=%d, code=%08X, %08X~%08X, mid=%08X, bit %d = %d\n", kb, buffer[kb], p0.m128i_i32[examined_plane&3], code.m128i_i32[examined_plane&3], start.m128i_u32[examined_plane&3], end.m128i_u32[examined_plane&3], middle.m128i_u32[examined_plane&3], examined_plane, bit.m128i_i32[examined_plane&3]==-1);
			//	int LOL_1=0;
#endif
			
#ifdef PREDICTOR_AWM_CONFBIT_WND
			//update hitcount: hit = !bit != (p0<0x8000),	hitcount = hit<<15|hitcount>>1
			__m128i b2=_mm_cmplt_epi32(p0, probhalf);
			b2=_mm_xor_si128(b2, bit_c);
			b2=_mm_and_si128(b2, probhalf);
			hitcount=_mm_srli_epi32(hitcount, 1);
			hitcount=_mm_or_si128(hitcount, b2);
#endif
#ifdef PREDICTOR_AWM_CONFBIT
			//update hitcount = !bit != (p0<0x8000)
			__m128i b2=_mm_cmplt_epi32(p0, probhalf);
			b2=_mm_xor_si128(b2, bit_c);
			hitcount=_mm_sub_epi32(hitcount, b2);
#endif

			//update weighted sum W(0) = !bit_new | W_prev(0)>>1	 = (sum i=1 to 16: !bit[-i]*2^i) / 2^16
			prob=_mm_srli_epi32(prob, 1);
			__m128i bit15=_mm_and_si128(bit_c, probhalf);
			prob=_mm_or_si128(prob, bit15);

			//select range: if(bit)start=middle; else end=middle-1;		same as: start=bit?middle:start, end=bit?end:middle-1;
		/*	__m128i temp=_mm_and_si128(middle, bit);
			start=_mm_and_si128(start, bit_c);
			start=_mm_or_si128(start, temp);

			temp=_mm_and_si128(middle, bit_c);
			end=_mm_and_si128(end, bit);
			end=_mm_or_si128(end, temp);//*/

			__m128i midm1=_mm_sub_epi32(middle, m_one);
			__m128i midp1=_mm_add_epi32(middle, m_one);
			
			midp1=_mm_and_si128(midp1, bit);
			start=_mm_and_si128(start, bit_c);
			start=_mm_or_si128(start, midp1);
			//middle=_mm_and_si128(middle, bit);
			//start=_mm_and_si128(start, bit_c);
			//start=_mm_or_si128(start, middle);

			midm1=_mm_and_si128(midm1, bit_c);
			end=_mm_and_si128(end, bit);
			end=_mm_or_si128(end, midm1);//*/

			//extract bits
			int bits=_mm_movemask_ps(_mm_castsi128_ps(bit));
			buffer[kb]|=bits<<kp;

			++kb;
			m_kb=_mm_add_epi32(m_kb, m_one);
			
			//check if most-significant byte is ready
			__m128i diff=_mm_xor_si128(start, end);
			diff=_mm_srli_epi32(diff, 24);
			diff=_mm_cmpeq_epi32(diff, _mm_setzero_si128());
			int ready_mask=_mm_movemask_ps(_mm_castsi128_ps(diff));
			if(ready_mask)
			{
				for(int kch=0;kch<4;++kch)
				{
					if(ready_mask>>kch&1)
					{
						unsigned istart=start.m128i_u32[kch], iend=end.m128i_u32[kch];
						auto plane=data+plane_offsets[kp+kch];
						do
						{
							code.m128i_u32[kch]=code.m128i_u32[kch]<<8|(unsigned char)plane[kc.m128i_u32[kch]];
							++kc.m128i_u32[kch];
							istart<<=8, iend=iend<<8|0xFF;
							i_outbitcounts.m128i_u32[kch]+=8;
							f_outbitcounts=_mm_cvtepi32_ps(i_outbitcounts);
						}while((istart^iend)<0x1000000);
						start.m128i_u32[kch]=istart, end.m128i_u32[kch]=iend;
					}
				}
			}
		}//end pixel loop
	}//end bit-plane loop
	delete[] plane_offsets, plane_offsets=nullptr;

	auto t2=__rdtsc();
	if(loud)
	{
		printf("AC decode:  %lld cycles\n", t2-t1);
	}
}

void			abac_encode_avx2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud)
{
	if(!imsize)
		return;
#ifdef HARD_AVX2_PROFILE
	long long ti[8]={}, ti1, ti2;
#endif
	auto t1=__rdtsc();

	std::vector<std::string> planes(depth);
	__m256i bitmask=_mm256_set_epi32(128, 64, 32, 16, 8, 4, 2, 1);
	__m256i m_ones=_mm256_set1_epi32(-1);
	__m256i probhalf=_mm256_set1_epi32(0x8000);
	__m256i msbyte_limit=_mm256_set1_epi32(0x1000000);
	__m256i m_one=_mm256_set1_epi32(1), m_probmax=_mm256_set1_epi32(prob_max);
	__m256 f_one=_mm256_set1_ps(1);
	for(int k=0;k<depth;++k)
		planes[k].reserve(imsize>>8);
	for(int kp=0;kp+7<depth;kp+=8)//bit-plane loop: encode 4 planes simultaneously
	{
		__m256i i_outbitcounts=_mm256_setzero_si256();
		__m256 f_outbitcounts=_mm256_cvtepi32_ps(i_outbitcounts);
		__m256i prob=_mm256_set1_epi32(1<<(LOG_WINDOW_SIZE-1));//estimated P(next bit = 0): weighted sum (sum i=1 to 16: !bit[-i]*2^i)
		__m256i m_kb=_mm256_setzero_si256();
		__m256i start=_mm256_setzero_si256(), end=_mm256_set1_epi32(-1);//inclusive range
		for(int kb=0;kb<imsize;)
		{
#ifdef HARD_AVX2_PROFILE
			int part=0;
			ti1=__rdtsc();
#endif
			__m256i range=_mm256_sub_epi32(end, start);
			__m256i eq_one=_mm256_cmpeq_epi32(range, _mm256_set1_epi32(1));
			int condition_eq_one=_mm256_movemask_ps(_mm256_castsi256_ps(eq_one));
			if(condition_eq_one)//because 1=0.9999...
			{
				for(int kch=0;kch<8;++kch)
				{
					if(condition_eq_one>>kch&1)
					{
						unsigned istart=start.m256i_u32[kch];
						auto &plane=planes[depth-1-(kp+kch)];
						plane.push_back(istart>>24);
						plane.push_back(istart>>16&0xFF);
						plane.push_back(istart>>8&0xFF);
						plane.push_back(istart&0xFF);
						start.m256i_u32[kch]=0, end.m256i_u32[kch]=0xFFFFFFFF;
					}
				}
				range=_mm256_sub_epi32(end, start);
			}
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

		/*	//calculate p0=0x8000+(long long)(p0-0x8000)*kb/conf_den;
			__m256i p0=_mm256_sub_epi32(prob, probhalf);
			__m256i lo=_mm256_mul_epi32(p0, m_kb);//SSE4.1 equivalent
			__m256i hi=_mm256_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			hi=_mm256_mul_epi32(hi, m_kb);//SSE4.1 equivalent
			__m256i den=_mm256_add_epi32(i_outbitcounts, m_kb);
			if(den.m256i_u32[0])lo.m256i_i64[0]/=den.m256i_u32[0];
			if(den.m256i_u32[1])hi.m256i_i64[0]/=den.m256i_u32[1];
			if(den.m256i_u32[2])lo.m256i_i64[1]/=den.m256i_u32[2];
			if(den.m256i_u32[3])hi.m256i_i64[1]/=den.m256i_u32[3];
			if(den.m256i_u32[4])lo.m256i_i64[2]/=den.m256i_u32[4];
			if(den.m256i_u32[5])hi.m256i_i64[2]/=den.m256i_u32[5];
			if(den.m256i_u32[6])lo.m256i_i64[3]/=den.m256i_u32[6];
			if(den.m256i_u32[7])hi.m256i_i64[3]/=den.m256i_u32[7];
			p0=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(lo), _mm256_castsi256_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			p0=_mm256_shuffle_epi32(p0, _MM_SHUFFLE(3, 1, 2, 0));
			p0=_mm256_add_epi32(p0, probhalf);
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

			//if confidence denominator==0 then confidence=100%
			den=_mm256_cmpeq_epi32(den, _mm256_setzero_si256());
			__m256i mask=_mm256_xor_si256(den, m_ones);
			p0=_mm256_and_si256(p0, mask);
			den=_mm256_and_si256(den, prob);
			p0=_mm256_or_si256(p0, den);//*/
			
			//calculate confidence = kb / (8kc+kb)		where kb: in bit idx, kc: out byte idx
			__m256 conf_num=_mm256_set1_ps((float)kb);
			__m256 conf_den=_mm256_add_ps(conf_num, f_outbitcounts);
			__m256 conf=_mm256_div_ps(conf_num, conf_den);

			//if confidence denominator==0 then confidence=100%
			__m256 den_nz=_mm256_cmp_ps(conf_den, _mm256_setzero_ps(), _CMP_NEQ_UQ);
			conf=_mm256_and_ps(conf, den_nz);
			den_nz=_mm256_xor_ps(den_nz, _mm256_castsi256_ps(m_ones));
			den_nz=_mm256_and_ps(den_nz, f_one);
			conf=_mm256_or_ps(conf, den_nz);

			//predict P(0) = 1/2 + (W(0)-1/2)*confidence
			__m256i p0=_mm256_sub_epi32(prob, probhalf);
			__m256 f_p0=_mm256_cvtepi32_ps(p0);
			f_p0=_mm256_mul_ps(f_p0, conf);
			f_p0=_mm256_floor_ps(f_p0);//SSE4.1
			p0=_mm256_cvtps_epi32(f_p0);
			p0=_mm256_add_epi32(probhalf, p0);//*/

#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

			//clamp probability: P(0) = clamp(2^-L, P(0), 1-2^-L)	where L=16		//probability shouldn't be exactly 0 or 100%
			__m256i c_min=_mm256_cmpgt_epi32(p0, m_one);
			__m256i c_max=_mm256_cmpgt_epi32(m_probmax, p0);
			p0=_mm256_and_si256(p0, c_min);
			p0=_mm256_and_si256(p0, c_max);
			c_min=_mm256_xor_si256(c_min, m_ones);
			c_max=_mm256_xor_si256(c_max, m_ones);
			c_min=_mm256_and_si256(c_min, m_one);
			c_max=_mm256_and_si256(c_max, m_probmax);
			p0=_mm256_or_si256(p0, c_min);
			p0=_mm256_or_si256(p0, c_max);
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

#ifdef DEBUG_SIMD_ENC
			if(examined_plane>=kp&&examined_plane<kp+8)
			{
				auto &p=probs[kb];
				if(p0.m256i_u32[examined_plane&7]!=p.p0)
					int LOL_1=0;
			}
			//if(kp==0&&kb==5)
			//if(1>=kp&&1<kp+4&&(kb==904||kb==905))//
			//if(5>=kp&&5<kp+4&&kb==1)//
			//	int LOL_1=0;//
#endif

			//calculate middle = start+((u64)(end-start)*p0>>L)
			__m256i lo=_mm256_mul_epu32(range, p0);
			lo=_mm256_srli_epi64(lo, 16);
			range=_mm256_shuffle_epi32(range, _MM_SHUFFLE(2, 3, 0, 1));
			p0=_mm256_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			__m256i hi=_mm256_mul_epu32(range, p0);
			hi=_mm256_srli_epi64(hi, 16);
			__m256i middle=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(lo), _mm256_castsi256_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			middle=_mm256_shuffle_epi32(middle, _MM_SHUFFLE(3, 1, 2, 0));
			middle=_mm256_add_epi32(middle, start);
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

			//get next bit
			__m256i bit_c=_mm256_set1_epi32(buffer[kb]);
			bit_c=_mm256_srli_epi32(bit_c, kp);
			bit_c=_mm256_and_si256(bit_c, bitmask);
			bit_c=_mm256_cmpeq_epi32(bit_c, _mm256_setzero_si256());
			__m256i bit=_mm256_xor_si256(bit_c, m_ones);
#ifdef DEBUG_SIMD_DEC
			if(examined_plane>=kp&&examined_plane<kp+4&&kb==examined_pixel)
			{
				printf("buffer[%d]=%04X, p0=%d, %08X~%08X, mid=%08X\n", kb, buffer[kb], p0.m256i_i32[examined_plane&3], start.m256i_u32[examined_plane&3], end.m256i_u32[examined_plane&3], middle.m256i_u32[examined_plane&3]);
				examine_flag=1;
			}
#endif
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

			//update weighted sum W(0) = !bit_new | W_prev(0)>>1	 = (sum i=1 to 16: !bit[-i]*2^i) / 2^16
			prob=_mm256_srli_epi32(prob, 1);
			__m256i bit15=_mm256_and_si256(bit_c, probhalf);
			prob=_mm256_or_si256(prob, bit15);

			//select range: if(bit)start=middle; else end=middle-1;		same as: start=bit?middle:start, end=bit?end:middle-1;
			__m256i midm1=_mm256_sub_epi32(middle, m_one);
			__m256i midp1=_mm256_add_epi32(middle, m_one);

			midp1=_mm256_and_si256(midp1, bit);
			start=_mm256_and_si256(start, bit_c);
			start=_mm256_or_si256(start, midp1);

			midm1=_mm256_and_si256(midm1, bit_c);
			end=_mm256_and_si256(end, bit);
			end=_mm256_or_si256(end, midm1);

			++kb;
			m_kb=_mm256_add_epi32(m_kb, m_one);
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif

			//check if most-significant byte is ready
			__m256i diff=_mm256_xor_si256(start, end);
			diff=_mm256_srli_epi32(diff, 24);
			diff=_mm256_cmpeq_epi32(diff, _mm256_setzero_si256());
			int ready_mask=_mm256_movemask_ps(_mm256_castsi256_ps(diff));
			if(ready_mask)
			{
				for(int kch=0;kch<8;++kch)
				{
					if(ready_mask>>kch&1)
					{
						unsigned istart=start.m256i_u32[kch], iend=end.m256i_u32[kch];
						auto &plane=planes[depth-1-(kp+kch)];
						do
						{
							plane.push_back(istart>>24);
#ifdef DEBUG_SIMD_DEC
							if(examine_flag&&kp+kch==examined_plane)
								printf("Shift-out: %02X\n", istart>>24);
#endif
							istart<<=8, iend=iend<<8|0xFF;
							i_outbitcounts.m256i_u32[kch]+=8;
							f_outbitcounts=_mm256_cvtepi32_ps(i_outbitcounts);
						}while((istart^iend)<0x1000000);
						start.m256i_u32[kch]=istart, end.m256i_u32[kch]=iend;
#ifdef DEBUG_SIMD_DEC
						if(examine_flag&&kp+kch==examined_plane)
						{
							printf("New bounds: %08X~%08X\n", start.m256i_u32[kch], end.m256i_u32[kch]);
							examine_flag=0;
						}
#endif
					}
				}
			}
#ifdef HARD_AVX2_PROFILE
			ti2=__rdtsc();
			ti[part]+=ti2-ti1, ++part;
			ti1=__rdtsc();
#endif
		}
		for(int kch=0;kch<8;++kch)
		{
			auto &plane=planes[depth-1-(kp+kch)];
			unsigned istart=start.m256i_u32[kch];
			plane.push_back(istart>>24&0xFF);//big-endian
			plane.push_back(istart>>16&0xFF);
			plane.push_back(istart>>8&0xFF);
			plane.push_back(istart&0xFF);
		}
	}
	auto t_enc=__rdtsc();
	//out_sizes.resize(depth);
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
	out_data.clear();
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}

	auto t2=__rdtsc();
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)out_data.size()<<3;
		printf("AC encode AVX2:  %lld cycles (Enc: %lld cycles)\n", t2-t1, t_enc-t1);
		printf("Size: %d -> %d, ratio: %lf\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize);
		printf("Bit\tkp\tSymbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%2d\t%5d\n", depth-1-k, k, out_sizes[k]);
		
		printf("Preview:\n");
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
#ifdef HARD_AVX2_PROFILE
		long long sum=0;
#define SIZEOF(STATIC_ARRAY)	sizeof(STATIC_ARRAY)/sizeof(*STATIC_ARRAY)
		for(int k=0;k<SIZEOF(ti);++k)
			sum+=ti[k];
		printf("AVX2 Profile:\n");
		for(int k=0;k<SIZEOF(ti);++k)
			printf("%d: %10lld, %9lf%%\n", k, ti[k], 100.*((double)ti[k]/sum));
#endif
	}
}
void			abac_decode_avx2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));

	int *plane_offsets=new int[depth];
	for(int kp=depth-1, cusize=0;kp>=0;--kp)
	{
		plane_offsets[kp]=cusize;
		cusize+=sizes[depth-1-kp];
	}
	
	__m256i bitmask=_mm256_set_epi32(128, 64, 32, 16, 8, 4, 2, 1);
	__m256i m_ones=_mm256_set1_epi32(-1);
	__m256i probhalf=_mm256_set1_epi32(0x8000);
	__m256i msbyte_limit=_mm256_set1_epi32(0x1000000);
	__m256i m_one=_mm256_set1_epi32(1), m_probmax=_mm256_set1_epi32(prob_max);
	__m256 f_one=_mm256_set1_ps(1);
	for(int kp=0, cusize=0;kp+7<depth;kp+=8)//bit-plane loop
	{
		__m256i i_outbitcounts=_mm256_setzero_si256();
		__m256 f_outbitcounts=_mm256_cvtepi32_ps(i_outbitcounts);
		__m256i prob=_mm256_set1_epi32(1<<(LOG_WINDOW_SIZE-1));//estimated P(next bit = 0): weighted sum (sum i=1 to 16: !bit[-i]*2^i)
		__m256i start=_mm256_setzero_si256(), end=_mm256_set1_epi32(-1);//inclusive range
		__m256i code=_mm256_set_epi32(
			load32_big((unsigned char*)data+plane_offsets[kp+7]),
			load32_big((unsigned char*)data+plane_offsets[kp+6]),
			load32_big((unsigned char*)data+plane_offsets[kp+5]),
			load32_big((unsigned char*)data+plane_offsets[kp+4]),
			load32_big((unsigned char*)data+plane_offsets[kp+3]),
			load32_big((unsigned char*)data+plane_offsets[kp+2]),
			load32_big((unsigned char*)data+plane_offsets[kp+1]),
			load32_big((unsigned char*)data+plane_offsets[kp]));
		__m256i kc=_mm256_set1_epi32(4);
		for(int kb=0;kb<imsize;)//bit-pixel loop
		{
			__m256i range=_mm256_sub_epi32(end, start);
			__m256i eq_one=_mm256_cmpeq_epi32(range, _mm256_set1_epi32(1));
			int condition_eq_one=_mm256_movemask_ps(_mm256_castsi256_ps(eq_one));
			if(condition_eq_one)//because 1=0.9999...
			{
				for(int kch=0;kch<8;++kch)
				{
					if(condition_eq_one>>kch&1)
					{
						auto plane=data+plane_offsets[kp+kch];
						code.m256i_u32[kch]=load32_big((unsigned char*)plane+kc.m256i_i32[kch]);
						kc.m256i_i32[kch]+=4;
						start.m256i_u32[kch]=0, end.m256i_u32[kch]=0xFFFFFFFF;
					}
				}
				range=_mm256_sub_epi32(end, start);
			}
			
			//calculate confidence = kb / (8kc+kb)		where kb: in bit idx, kc: out byte idx
			__m256 conf_num=_mm256_set1_ps((float)kb);
			__m256 conf_den=_mm256_add_ps(conf_num, f_outbitcounts);
			__m256 conf=_mm256_div_ps(conf_num, conf_den);

			//if confidence denominator==0 then confidence=100%
			__m256 den_nz=_mm256_cmp_ps(conf_den, _mm256_setzero_ps(), _CMP_NEQ_UQ);
			conf=_mm256_and_ps(conf, den_nz);
			den_nz=_mm256_xor_ps(den_nz, _mm256_castsi256_ps(m_ones));
			den_nz=_mm256_and_ps(den_nz, f_one);
			conf=_mm256_or_ps(conf, den_nz);

			//predict P(0) = 1/2 + (W(0)-1/2)*confidence
			__m256i p0=_mm256_sub_epi32(prob, probhalf);
			__m256 f_p0=_mm256_cvtepi32_ps(p0);
			f_p0=_mm256_mul_ps(f_p0, conf);
			f_p0=_mm256_floor_ps(f_p0);//SSE4.1 equivalent
			p0=_mm256_cvtps_epi32(f_p0);
			p0=_mm256_add_epi32(probhalf, p0);//*/

			//clamp probability: P(0) = clamp(2^-L, P(0), 1-2^-L)	where L=16		//probability shouldn't be exactly 0 or 100%
			__m256i c_min=_mm256_cmpgt_epi32(p0, m_one);
			__m256i c_max=_mm256_cmpgt_epi32(m_probmax, p0);
			p0=_mm256_and_si256(p0, c_min);
			p0=_mm256_and_si256(p0, c_max);
			c_min=_mm256_xor_si256(c_min, m_ones);
			c_max=_mm256_xor_si256(c_max, m_ones);
			c_min=_mm256_and_si256(c_min, m_one);
			c_max=_mm256_and_si256(c_max, m_probmax);
			p0=_mm256_or_si256(p0, c_min);
			p0=_mm256_or_si256(p0, c_max);
			
			//calculate middle = start+((u64)(end-start)*p0>>L)
			__m256i lo=_mm256_mul_epu32(range, p0);
			lo=_mm256_srli_epi64(lo, 16);
			range=_mm256_shuffle_epi32(range, _MM_SHUFFLE(2, 3, 0, 1));
			p0=_mm256_shuffle_epi32(p0, _MM_SHUFFLE(2, 3, 0, 1));
			__m256i hi=_mm256_mul_epu32(range, p0);
			hi=_mm256_srli_epi64(hi, 16);
			__m256i middle=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(lo), _mm256_castsi256_ps(hi), _MM_SHUFFLE(2, 0, 2, 0)));
			middle=_mm256_shuffle_epi32(middle, _MM_SHUFFLE(3, 1, 2, 0));
			middle=_mm256_add_epi32(middle, start);

			//unpack next bit	bit=code>=middle	unsigned compare
			__m256i bit=_mm256_max_epu32(code, middle);//SSE4.1 equivalent
			bit=_mm256_cmpeq_epi32(bit, code);
			__m256i bit_c=_mm256_xor_si256(bit, m_ones);

#ifdef DEBUG_SIMD_DEC
			if(examined_plane>=kp&&examined_plane<kp+4&&kb==examined_pixel)
				printf("buffer[%d]=%04X, p0=%d, code=%08X, %08X~%08X, mid=%08X, bit %d = %d\n", kb, buffer[kb], p0.m256i_i32[examined_plane&3], code.m256i_i32[examined_plane&3], start.m256i_u32[examined_plane&3], end.m256i_u32[examined_plane&3], middle.m256i_u32[examined_plane&3], examined_plane, bit.m256i_i32[examined_plane&3]==-1);
#endif

			//update weighted sum W(0) = !bit_new | W_prev(0)>>1	 = (sum i=1 to 16: !bit[-i]*2^i) / 2^16
			prob=_mm256_srli_epi32(prob, 1);
			__m256i bit15=_mm256_and_si256(bit_c, probhalf);
			prob=_mm256_or_si256(prob, bit15);

			//select range: if(bit)start=middle; else end=middle-1;		same as: start=bit?middle:start, end=bit?end:middle-1;
			__m256i midm1=_mm256_sub_epi32(middle, m_one);
			__m256i midp1=_mm256_add_epi32(middle, m_one);
			
			midp1=_mm256_and_si256(midp1, bit);
			start=_mm256_and_si256(start, bit_c);
			start=_mm256_or_si256(start, midp1);

			midm1=_mm256_and_si256(midm1, bit_c);
			end=_mm256_and_si256(end, bit);
			end=_mm256_or_si256(end, midm1);

			int bits=_mm256_movemask_ps(_mm256_castsi256_ps(bit));
			buffer[kb]|=bits<<kp;

			++kb;
			
			//check if most-significant byte is ready
			__m256i diff=_mm256_xor_si256(start, end);
			diff=_mm256_srli_epi32(diff, 24);
			diff=_mm256_cmpeq_epi32(diff, _mm256_setzero_si256());
			int ready_mask=_mm256_movemask_ps(_mm256_castsi256_ps(diff));
			if(ready_mask)
			{
				for(int kch=0;kch<8;++kch)
				{
					if(ready_mask>>kch&1)
					{
						unsigned istart=start.m256i_u32[kch], iend=end.m256i_u32[kch];
						auto plane=data+plane_offsets[kp+kch];
						do
						{
							code.m256i_u32[kch]=code.m256i_u32[kch]<<8|(unsigned char)plane[kc.m256i_u32[kch]];
							++kc.m256i_u32[kch];
							istart<<=8, iend=iend<<8|0xFF;
							i_outbitcounts.m256i_u32[kch]+=8;
							f_outbitcounts=_mm256_cvtepi32_ps(i_outbitcounts);
						}while((istart^iend)<0x1000000);
						start.m256i_u32[kch]=istart, end.m256i_u32[kch]=iend;
					}
				}
			}
		}//end pixel loop
	}//end bit-plane loop
	delete[] plane_offsets, plane_offsets=nullptr;

	auto t2=__rdtsc();
	if(loud)
	{
		printf("AC decode:  %lld cycles\n", t2-t1);
	}
}
#endif

	#define			ESTIMATE_CONFBIT
//	#define			ESTIMATE_PREVBIT

int				abac_estimate(const void *src, int imsize, int depth, int bytestride, bool loud, int *sizes)
{
	auto buffer=(const unsigned char*)src;
	int compr_size=0;
	if(!imsize)
		return compr_size;
	auto t1=__rdtsc();
	
	auto planesizes=new int[depth];
	memset(planesizes, 0, depth*sizeof(int));
//	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		int bit_offset=kp>>3, bit_shift=kp&7;
	//	int bit_offset=kp/(bytestride<<3), bit_shift=kp%(bytestride<<3);
	//	auto &plane=planes[depth-1-kp];
	//	plane.reserve(imsize>>8);
#ifdef ESTIMATE_CONFBIT
		int prob=1<<(LOG_WINDOW_SIZE-1), hitcount=1;//cheap weighted average predictor
#endif
#ifdef ESTIMATE_PREVBIT
		int prevbit=0;
#endif

		unsigned start=0, end=0xFFFFFFFF;
		for(int kb=0, kb2=0;kb<imsize;kb2+=bytestride)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
		{
			int bit=buffer[kb2+bit_offset]>>bit_shift&1;
			u64 range=end-start;
			if(range<3)
			{
				planesizes[depth-1-kp]+=4;
			//	plane.push_back(start>>24);
			//	plane.push_back(start>>16&0xFF);
			//	plane.push_back(start>>8&0xFF);
			//	plane.push_back(start&0xFF);
				start=0, end=0xFFFFFFFF;//because 1=0.9999...
				range=end-start;
			}
#ifdef ESTIMATE_CONFBIT
			int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
			p0=clamp(1, p0, prob_max);
#endif
#ifdef ESTIMATE_PREVBIT
			int p0=1+(0xFFFD&-!prevbit);
#endif
			unsigned middle=start+(unsigned)(range*p0>>16);
			middle+=(middle==start)-(middle==end);
#ifdef ESTIMATE_CONFBIT
			hitcount+=bit^(p0>=0x8000);
			prob=!bit<<15|prob>>1;//cheap weighted average predictor		P(0) = (sum i=1 to 16: !b[-i]*2^i) / 2^16
#endif
			if(bit)
				start=middle+1;
			else
				end=middle-1;
#ifdef ESTIMATE_PREVBIT
			prevbit=bit;
#endif
			++kb;
			while((start^end)<0x1000000)//most significant byte has stabilized			zpaq 1.10
			{
				++planesizes[depth-1-kp];
			//	plane.push_back(start>>24);
				start<<=8;
				end=end<<8|0xFF;
			}
		}
		planesizes[depth-1-kp]+=4;
	//	plane.push_back(start>>24&0xFF);//big-endian
	//	plane.push_back(start>>16&0xFF);
	//	plane.push_back(start>>8&0xFF);
	//	plane.push_back(start&0xFF);
	}
/*	auto t_enc=__rdtsc();
	out_data.clear();
	for(int k=0;k<depth;++k)
		out_sizes[k]=planes[k].size();
	for(int k=0;k<depth;++k)
	{
		auto &plane=planes[k];
		out_data.insert(out_data.end(), plane.begin(), plane.end());
	}//*/
	auto t2=__rdtsc();

	for(int kp=0;kp<depth;++kp)
		compr_size+=planesizes[depth-1-kp];
	if(sizes)
		memcpy(sizes, planesizes, depth*sizeof(int));
	if(loud)
	{
		int original_bitsize=imsize*depth, compressed_bitsize=(int)compr_size<<3;
		printf("AC estimate:  %lld cycles\n", t2-t1);
		printf("Size: %d -> %d, ratio: %lf\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize);
		printf("Bit\tkp\tSymbol count\n");
		for(int k=0;k<depth;++k)
			printf("%2d\t%2d\t%5d\n", depth-1-k, k, planesizes[k]);
		
		//printf("Preview:\n");
		//int kprint=out_data.size()<200?out_data.size():200;
		//for(int k=0;k<kprint;++k)
		//	printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
	delete[] planesizes;
	return compr_size;
}

#if 1
	#define		LOG_CONFBOOST	14
	#define		ABAC2_CONF_MSB_RELATION

//	#define		DEBUG_ABAC2

const double	boost_power=4, min_conf=0.55;
#ifdef DEBUG_ABAC2
const int		examined_plane=4, examined_start=0, examined_end=100;
#endif
/*int			test_conf[]=
{
	213316,
	211642,
	227643,
	272892,
	318084,
	349064,
	369961,
	383717,
};//*/
static int		abac2_normalize16(int x, int logx)
{
	if(logx>15)
		x>>=logx-15;
	else if(logx<15)
		x<<=15-logx;
	return x;
}
static void		abac2_invimsize(int imsize, int &logimsize, int &invimsize)
{
	logimsize=floor_log2(imsize);
	invimsize=abac2_normalize16(imsize, logimsize);
	invimsize=(1<<31)/invimsize;
}
void			abac2_encode(const void *src, int imsize, int depth, int bytestride, std::string &out_data, int *out_sizes, int *out_conf, bool loud)
{
	auto buffer=(const unsigned char*)src;
	if(!imsize)
		return;
	auto t1=__rdtsc();
	
#ifdef MEASURE_PREDICTION
	u64 hitnum=0, hitden=0;//prediction efficiency
#endif
	//int logimsize=0, invimsize=0;
	//abac2_invimsize(imsize, logimsize, invimsize);

	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		int bit_offset=kp>>3, bit_shift=kp&7;
		int bit_offset2=(kp+1)>>3, bit_shift2=(kp+1)&7;
	//	int bit_offset=kp/(bytestride<<3), bit_shift=kp%(bytestride<<3);
	//	int bit_offset2=(kp+1)/(bytestride<<3), bit_shift2=(kp+1)%(bytestride<<3);
		auto &plane=planes[depth-1-kp];
		int prob=0x8000, prob_correct=0x8000;//cheap weighted average predictor

		u64 hitcount=1;

		for(int kb=0, kb2=0;kb<imsize;++kb, kb2+=bytestride)//analyze bitplane
		{
			int bit=buffer[kb2+bit_offset]>>bit_shift&1;
		//	int bit=buffer[kb]>>kp&1;
			int p0=((long long)(prob-0x8000)*prob_correct>>16);
			p0+=0x8000;
			//int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
			p0=clamp(1, p0, prob_max);
			int correct=bit^(p0>=0x8000);
			//if(kp==1)
			//	printf("%d", bit);//actual bits
			//	printf("%d", p0<0x8000);//predicted bits
			//	printf("%d", !correct);//prediction error
			hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
		}
		out_conf[depth-1-kp]=(int)hitcount;

		if(hitcount<imsize*min_conf)
		{
			plane.resize((imsize+7)>>3, 0);
			for(int kb=0, kb2=0, b=0;kb<imsize;++kb, kb2+=bytestride)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				int bit=buffer[kb2+bit_offset]>>bit_shift&1;
			//	int bit=buffer[kb]>>kp&1;
				plane[byte_idx]|=bit<<bit_idx;
			}
		}
		else
		{
			int hitratio_sure=int(0x10000*pow((double)hitcount/imsize, 1/boost_power)), hitratio_notsure=int(0x10000*pow((double)hitcount/imsize, boost_power));
			int hitratio_delta=hitratio_sure-hitratio_notsure;
		//	int hitratio_sure=int(0x10000*cbrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)imsize*imsize*imsize));
		//	int hitratio_sure=int(0x10000*sqrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)imsize*imsize));
			hitcount=(hitcount<<16)/imsize;

			//hitcount=unsigned(((u64)hitcount<<16)/imsize);
			//hitcount=abac2_normalize16(hitcount, logimsize);
			//hitcount*=invimsize;

			prob_correct=prob=0x8000;

#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit0=0;
#endif
			
			plane.reserve(imsize>>8);
			unsigned start=0;
			u64 range=0xFFFFFFFF;
			for(int kb=0, kb2=0;kb<imsize;kb2+=bytestride)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
			{
				int bit=buffer[kb2+bit_offset]>>bit_shift&1;
			//	int bit=buffer[kb]>>kp&1;
#ifdef ABAC2_CONF_MSB_RELATION
				int prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
			//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
				
				if(range<3)
				{
					plane.push_back(start>>24);
					plane.push_back(start>>16&0xFF);
					plane.push_back(start>>8&0xFF);
					plane.push_back(start&0xFF);
					start=0, range=0xFFFFFFFF;//because 1=0.9999...
				}
				
				int p0=prob-0x8000;
				p0=p0*prob_correct>>16;
				p0=p0*prob_correct>>16;
				int sure=-(prevbit==prevbit0);
				p0=p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
				//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
				//p0=(long long)p0*hitcount>>16;
				p0+=0x8000;
				//if(prevbit!=prevbit0)
				//	p0=0x8000;
				//	p0=0xFFFF-p0;

				//int p0=0x8000+((long long)(prob-0x8000)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

				//int p0=(long long)(prob-0x8000)*sqrthitcount>>16;
				//if(prevbit==prevbit0)
				//	p0=(long long)p0*hitcount>>16;
				//p0+=0x8000;

				//int confboost=prevbit==prevbit0;
				//confboost-=!confboost;
				//confboost<<=LOG_CONFBOOST;
				//int p0=0x8000+((long long)(prob-0x8000)*(hitcount+confboost)>>16);

			//	int p0=0x8000+(int)((prob-0x8000)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/imsize):(double)test_conf[kp]*test_conf[kp]/((double)imsize*imsize)));
			//	int p0=prevbit==prevbit0?prob:0x8000;
			//	int p0=0x8000+(long long)(prob-0x8000)*test_conf[kp]/imsize;
			//	int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
				p0=clamp(1, p0, prob_max);
				unsigned r2=(unsigned)(range*p0>>16);
				r2+=(r2==0)-(r2==range);
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("%6d %6d %d %08X+%08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, start+r2);
#endif

				int correct=bit^(p0>=0x8000);
			//	hitcount+=correct;
				prob=!bit<<15|prob>>1;
				prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
				prevbit0=prevbit;
#endif
#ifdef MEASURE_PREDICTION
				hitnum+=correct, ++hitden;
#endif
				auto start0=start;
				if(bit)
				{
					++r2;
					start+=r2, range-=r2;
				}
				//	start=middle+1;
				else
					range=r2-1;
				//	end=middle-1;
				if(start<start0)//
				{
					printf("OVERFLOW\nstart = %08X -> %08X, r2 = %08X", start0, start, r2);
					int k=0;
					scanf_s("%d", &k);
				}
				++kb;
				
				while((start^(start+(unsigned)range))<0x1000000)//most significant byte has stabilized			zpaq 1.10
				{
#ifdef DEBUG_ABAC2
					if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
						printf("range %08X byte-out %02X\n", (int)range, start>>24);
#endif
					plane.push_back(start>>24);
					start<<=8;
					range=range<<8|0xFF;
				}
			}
			plane.push_back(start>>24&0xFF);//big-endian
			plane.push_back(start>>16&0xFF);
			plane.push_back(start>>8&0xFF);
			plane.push_back(start&0xFF);
		}
		if(loud)
			printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, out_conf[depth-1-kp], imsize, 100.*out_conf[depth-1-kp]/imsize);
		//	printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, hitcount, imsize, 100.*hitcount/imsize);
	}
	auto t_enc=__rdtsc();
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
#ifdef MEASURE_PREDICTION
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", imsize>>3);
		for(int k=0;k<depth;++k)
			printf("%2d\t%5d\t%lf\n", depth-1-k, out_sizes[k], (double)imsize/(out_sizes[k]<<3));
		
		printf("Preview:\n");
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
}
void			abac2_decode(const char *data, const int *sizes, const int *conf, void *dst, int imsize, int depth, int bytestride, bool loud)
{
	auto buffer=(unsigned char*)dst;
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*bytestride);
	
	//int logimsize=0, invimsize=0;
	//abac2_invimsize(imsize, logimsize, invimsize);
	
	for(int kp=depth-1, cusize=0;kp>=0;--kp)//bit-plane loop
	{
		int bit_offset=kp>>3, bit_shift=kp&7;
		int bit_offset2=(kp+1)>>3, bit_shift2=(kp+1)&7;
	//	int bit_offset=kp/(bytestride<<3), bit_shift=kp%(bytestride<<3);
	//	int bit_offset2=(kp+1)/(bytestride<<3), bit_shift2=(kp+1)%(bytestride<<3);
		int ncodes=sizes[depth-1-kp];
		auto plane=data+cusize;
		
		int prob=0x8000, prob_correct=0x8000;
#if 1
		u64 hitcount=conf[depth-1-kp];
		if(hitcount<imsize*min_conf)
		{
			for(int kb=0, kb2=0, b=0;kb<imsize;++kb, kb2+=bytestride)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				int bit=plane[byte_idx]>>bit_idx&1;
				buffer[kb2+bit_offset]|=bit<<bit_shift;
			//	buffer[kb]|=bit<<kp;
			}
			cusize+=ncodes;
			continue;
		}
#ifdef ABAC2_CONF_MSB_RELATION
		int prevbit0=0;
#endif
		int hitratio_sure=int(0x10000*pow((double)hitcount/imsize, 1/boost_power)), hitratio_notsure=int(0x10000*pow((double)hitcount/imsize, boost_power));
		int hitratio_delta=hitratio_sure-hitratio_notsure;
		//int hitratio_sure=int(0x10000*cbrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)imsize*imsize*imsize));
		//int hitratio_sure=int(0x10000*sqrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)imsize*imsize));
		hitcount=(hitcount<<16)/imsize;
		//hitcount=unsigned(((u64)hitcount<<16)/imsize);
		//hitcount=abac2_normalize16(hitcount, logimsize);
		//hitcount*=invimsize;
#endif

		unsigned start=0;
		u64 range=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane);
		for(int kc=4, kb=0, kb2=0;kb<imsize;kb2+=bytestride)//bit-pixel loop
		{
			if(range<3)
			{
				code=load32_big((unsigned char*)plane+kc);
				kc+=4;
				start=0, range=0xFFFFFFFF;//because 1=0.9999...
			}
#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
		//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
			int p0=prob-0x8000;
			p0=p0*prob_correct>>16;
			p0=p0*prob_correct>>16;
			int sure=-(prevbit==prevbit0);
			p0=p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
			//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
			//p0=(long long)p0*hitcount>>16;
			p0+=0x8000;
			//if(prevbit!=prevbit0)
			//	p0=0x8000;
			//	p0=0xFFFF-p0;

			//int p0=0x8000+((long long)(prob-0x8000)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

			//int p0=(long long)(prob-0x8000)*sqrthitcount>>16;
			//if(prevbit==prevbit0)
			//	p0=(long long)p0*hitcount>>16;
			//p0+=0x8000;

			//int confboost=prevbit==prevbit0;
			//confboost-=!confboost;
			//confboost<<=LOG_CONFBOOST;
			//int p0=0x8000+((long long)(prob-0x8000)*(hitcount+confboost)>>16);

		//	int p0=0x8000+(int)((prob-0x8000)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/imsize):(double)test_conf[kp]*test_conf[kp]/((double)imsize*imsize)));
		//	int p0=prevbit==prevbit0?prob:0x8000;
		//	int p0=0x8000+(long long)(prob-0x8000)*test_conf[kp]/imsize;
		//	int p0=0x8000+(long long)(prob-0x8000)*hitcount/(kb+1);
			p0=clamp(1, p0, prob_max);
			unsigned r2=(unsigned)(range*p0>>16);
			r2+=(r2==0)-(r2==range);
			unsigned middle=start+r2;
			int bit=code>middle;
#ifdef DEBUG_ABAC2
			if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
				printf("%6d %6d %d %08X+%08X %08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, middle, code);
#endif
			
			int correct=bit^(p0>=0x8000);
		//	hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
			prevbit0=prevbit;
#endif
			
			if(bit)
			{
				++r2;
				start+=r2, range-=r2;
			}
			//	start=middle+1;
			else
				range=r2-1;
			//	end=middle-1;
			
			buffer[kb2+bit_offset]|=bit<<bit_shift;
		//	buffer[kb]|=bit<<kp;
			++kb;
			
			while((start^(start+(unsigned)range))<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range %08X byte-out %02X\n", (int)range, code>>24);
#endif
				code=code<<8|(unsigned char)plane[kc];
				++kc;
				start<<=8;
				range=range<<8|0xFF;
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
#endif

#if 1

	#define		ABAC3_ZPAQ0

void			abac3_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_conf, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	
#ifdef MEASURE_PREDICTION
	u64 hitnum=0, hitden=0;//prediction efficiency
#endif

	std::vector<std::string> planes(depth);
	for(int kp=depth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		auto &plane=planes[depth-1-kp];
#ifdef ABAC3_ZPAQ0
		int zpaq0_num=1, zpaq0_den=2;
#else
		int prob=0x8000, prob2=0x8000, prob2_parity=1, prob_correct=0x8000;
	//	int p0_bins[4]={1, 1, 1, 1}, countof0=4;
	//	int p0_bins[16]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, countof0=16;
	//	int p0_bins[16]={8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
		u64 hitcount=1;

		for(int kb=0;kb<imsize;++kb)//analyze bitplane
		{
			int bit=buffer[kb]>>kp&1, bit2=buffer[kb]>>(kp+1)&1;
			prob2_parity^=bit2^prob2&1;
			prob2=!bit2<<15|prob2>>1;

			//int p0=((long long)(p0_bins[prob]-(countof0>>2))*prob_correct>>16);
			//p0+=countof0>>2;
			//p0=clamp(1, p0, countof0>>2);
			//int correct=bit^(p0>=2);

			//int p0=((long long)(p0_bins[prob]-(countof0>>4))*prob_correct>>16);
			//p0+=countof0>>4;
			//p0=clamp(1, p0, countof0>>4);
			//int correct=bit^(p0>=8);

			//int p0=((long long)(p0_bins[prob]-8)*prob_correct>>16);
			//p0+=8;
			//p0=clamp(1, p0, 14);
			//int correct=bit^(p0>=8);

			int p0=((long long)(prob-0x8000)*prob_correct>>16);
			p0+=0x8000;
			p0=clamp(1, p0, 0xFFFE);
			int correct=bit^(p0>=0x8000);
			
			if(kb<120)
			//if(kp==0)
				printf("%d", bit);//actual bits
			//	printf("%d", p0<0x8000);//predicted bits
			//	printf("%d", !correct);//prediction error
			hitcount+=correct;
			//p0_bins[prob]+=!bit, countof0+=bit;
			//p0_bins[prob]=!bit<<3|p0_bins[prob]>>1;
			prob=!bit<<15|prob>>1;
			//prob=!bit<<3|prob>>1;
			//prob=!bit<<1|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
		}
		out_conf[depth-1-kp]=hitcount;

		if(hitcount<imsize*min_conf)
		{
			plane.resize((imsize+7)>>3, 0);
			for(int kb=0, b=0;kb<imsize;++kb)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				int bit=buffer[kb]>>kp&1;
				plane[byte_idx]|=bit<<bit_idx;
			}
			goto done;
		}
		
		int hitratio_sure=int(0x10000*pow((double)hitcount/imsize, 1/boost_power)), hitratio_notsure=int(0x10000*pow((double)hitcount/imsize, boost_power));
		int hitratio_delta=hitratio_sure-hitratio_notsure;
		hitcount=(hitcount<<16)/imsize;

		prob_correct=prob=0x8000;
#ifdef ABAC2_CONF_MSB_RELATION
		int prevbit0=0;
#endif
#endif
		
		plane.reserve(imsize>>8);
		unsigned start=0;
		u64 range=0xFFFFFFFF;
		for(int kb=0;kb<imsize;)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
		{
			int bit=buffer[kb]>>kp&1;
#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit=buffer[kb]>>(kp+1)&1;
#endif
			
			if(range<3)
			{
				plane.push_back(start>>24);
				plane.push_back(start>>16&0xFF);
				plane.push_back(start>>8&0xFF);
				plane.push_back(start&0xFF);
				start=0, range=0xFFFFFFFF;//because 1=0.9999...
			}
#ifdef ABAC3_ZPAQ0
			unsigned r2=(unsigned)(range*zpaq0_num/zpaq0_den);
			if(r2<=0)
				r2=1;
			else if(r2>=range)
				r2=(unsigned)range-1;
#else
			int p0=prob-0x8000;
			p0=p0*prob_correct>>16;
			p0=p0*prob_correct>>16;
			int sure=-(prevbit==prevbit0);
			p0=p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
			p0+=0x8000;

			p0=clamp(1, p0, prob_max);
			unsigned r2=(unsigned)(range*p0>>16);
			r2+=(r2==0)-(r2==range);
#endif
#ifdef DEBUG_ABAC2
			if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
				printf("%6d %6d %d %08X+%08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, start+r2);
#endif
			
#ifdef ABAC3_ZPAQ0
			int correct=bit^(zpaq0_num>=(zpaq0_den>>1));
			zpaq0_num+=!bit, ++zpaq0_den;
			int sh=zpaq0_den>0x10000;
			zpaq0_num>>=sh, zpaq0_den>>=sh;
#else
			int correct=bit^(p0>=0x8000);
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
			prevbit0=prevbit;
#endif
#endif
#ifdef MEASURE_PREDICTION
			hitnum+=correct, ++hitden;
#endif
			auto start0=start;
			if(bit)
			{
				++r2;
				start+=r2, range-=r2;
			}
			else
				range=r2-1;
			if(start<start0)//
			{
				printf("OVERFLOW\nstart = %08X -> %08X, r2 = %08X", start0, start, r2);
				int k=0;
				scanf_s("%d", &k);
			}
			++kb;
			
			while((start^(start+(unsigned)range))<0x1000000)//most significant byte has stabilized			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range %08X byte-out %02X\n", (int)range, start>>24);
#endif
				plane.push_back(start>>24);
				start<<=8;
				range=range<<8|0xFF;
			}
		}
		plane.push_back(start>>24&0xFF);//big-endian
		plane.push_back(start>>16&0xFF);
		plane.push_back(start>>8&0xFF);
		plane.push_back(start&0xFF);
#ifndef ABAC3_ZPAQ0
	done:
#endif
		if(loud)
			printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, out_conf[depth-1-kp], imsize, 100.*out_conf[depth-1-kp]/imsize);
		//	printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, hitcount, imsize, 100.*hitcount/imsize);
	}
	auto t_enc=__rdtsc();
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
#ifdef MEASURE_PREDICTION
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", imsize>>3);
		for(int k=0;k<depth;++k)
			printf("%2d\t%5d\t%lf\n", depth-1-k, out_sizes[k], (double)imsize/(out_sizes[k]<<3));
		
		printf("Preview:\n");
		int kprint=out_data.size()<200?out_data.size():200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out_data[k]&0xFF);
		printf("\n");
	}
}
void			abac3_decode(const char *data, const int *sizes, const int *conf, short *buffer, int imsize, int depth, bool loud)
{
	if(!imsize)
		return;
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(short));
	
	for(int kp=depth-1, cusize=0;kp>=0;--kp)//bit-plane loop
	{
		int ncodes=sizes[depth-1-kp];
		auto plane=data+cusize;
		
#ifdef ABAC3_ZPAQ0
		int zpaq0_num=1, zpaq0_den=2;
#else
		int prob=0x8000, prob_correct=0x8000;
		u64 hitcount=conf[depth-1-kp];
		if(hitcount<imsize*min_conf)
		{
			for(int kb=0, b=0;kb<imsize;++kb)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				int bit=plane[byte_idx]>>bit_idx&1;
				buffer[kb]|=bit<<kp;
			}
			cusize+=ncodes;
			continue;
		}
		int prevbit0=0;
		int hitratio_sure=int(0x10000*pow((double)hitcount/imsize, 1/boost_power)), hitratio_notsure=int(0x10000*pow((double)hitcount/imsize, boost_power));
		int hitratio_delta=hitratio_sure-hitratio_notsure;
		hitcount=(hitcount<<16)/imsize;
#endif
		unsigned start=0;
		u64 range=0xFFFFFFFF;
		unsigned code=load32_big((unsigned char*)plane);
		for(int kc=4, kb=0;kb<imsize;)//bit-pixel loop
		{
			if(range<3)
			{
				code=load32_big((unsigned char*)plane+kc);
				kc+=4;
				start=0, range=0xFFFFFFFF;//because 1=0.9999...
			}
#ifdef ABAC3_ZPAQ0
			unsigned r2=(unsigned)(range*zpaq0_num/zpaq0_den);
			if(r2<=0)
				r2=1;
			else if(r2>=range)
				r2=(unsigned)range-1;
#else
			int prevbit=buffer[kb]>>(kp+1)&1;
			int p0=prob-0x8000;
			p0=p0*prob_correct>>16;
			p0=p0*prob_correct>>16;
			int sure=-(prevbit==prevbit0);
			p0=p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
			p0+=0x8000;

			p0=clamp(1, p0, prob_max);
			unsigned r2=(unsigned)(range*p0>>16);
			r2+=(r2==0)-(r2==range);
#endif
			unsigned middle=start+r2;
			int bit=code>middle;
#ifdef DEBUG_ABAC2
			if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
				printf("%6d %6d %d %08X+%08X %08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, middle, code);
#endif
			
#ifdef ABAC3_ZPAQ0
			int correct=bit^(zpaq0_num>=(zpaq0_den>>1));
			zpaq0_num+=!bit, ++zpaq0_den;
			int sh=zpaq0_den>0x10000;
			zpaq0_num>>=sh, zpaq0_den>>=sh;
#else
			int correct=bit^(p0>=0x8000);
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
			prevbit0=prevbit;
#endif
#endif
			if(bit)
			{
				++r2;
				start+=r2, range-=r2;
			}
			else
				range=r2-1;
			
			buffer[kb]|=bit<<kp;
			++kb;
			
			while((start^(start+(unsigned)range))<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range %08X byte-out %02X\n", (int)range, code>>24);
#endif
				code=code<<8|(unsigned char)plane[kc];
				++kc;
				start<<=8;
				range=range<<8|0xFF;
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
#endif