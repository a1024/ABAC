#ifndef ABAC_H
#define ABAC_H
#include<stdlib.h>
#include<string>

#include<stdio.h>
#include<stdarg.h>
#include<vector>
#include<algorithm>

#if defined __GNUC__ ||1
#define			NO_SIMD
#endif

	#define		ENABLE_ASSERT

typedef unsigned long long u64;
int				floor_log2(unsigned long long n);

static bool		crash(const char *file, int line, const char *condition, const char *format, ...)
{
	printf("\nCRASH\n");
	printf("%s(%d):\n\t( %s ) == 0\n\n", file, line, condition);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	int k=0;
	scanf_s("%d", &k);
	exit(1);
	return false;
}
#ifdef ENABLE_ASSERT
#define			MY_ASSERT(SUCCESS, FORMAT, ...)		((SUCCESS)!=0||crash(__FILE__, __LINE__, #SUCCESS, FORMAT,##__VA_ARGS__))
#else
#define			MY_ASSERT(...)
#endif

//Asymmetric Numeral Systems coder		J Duda 2009			NEED TO STORE HISTOGRAM FOR EACH CHANNEL
//rANS	https://github.com/rygorous/ryg_rans
//tANS	https://github.com/JarekDuda/AsymmetricNumeralSystemsToolkit

//	#define		CHECK_ANS_STATE		//probably useless
//	#define		ANS_PRINT_RGB		//only for debugging
//	#define		ANS_PRINT_HISTOGRAM	//annoying

#define			ANS_PROB_BITS	16
#define			ANS_L			(1<<ANS_PROB_BITS)
struct			SortedHistInfo
{
	int idx, freq, qfreq;
	SortedHistInfo():idx(0), freq(0), qfreq(0){}
};
//template because same code is repeated with different compile-time bytestride depending on type of buffer
template<typename T>inline void ans_calc_histogram(const T *buffer, int imsize, int bit0, int depth, int *histogram, int *CDF, int *CDF2sym, bool loud)
{
	int nlevels=1<<depth, mask=nlevels-1;
	MY_ASSERT(nlevels<ANS_L, "Too many levels");//what if nlevels = 2^N-1 ?
	if(!imsize)
	{
		memset(histogram, 0, nlevels*sizeof(*histogram));
		return;
	}
	std::vector<SortedHistInfo> h(nlevels);
	for(int k=0;k<nlevels;++k)
		h[k].idx=k;
	for(int k=0;k<imsize;++k)
		++h[buffer[k]>>bit0&mask].freq;
	if(imsize!=ANS_L)
	{
		for(int k=0;k<nlevels;++k)
			h[k].qfreq=((long long)h[k].freq<<ANS_PROB_BITS)/imsize;

		std::sort(h.begin(), h.end(), [](SortedHistInfo const &a, SortedHistInfo const &b)
		{
			return a.freq<b.freq;
		});
		int frac_start=0, one_start=0;
		for(frac_start=0;frac_start<nlevels&&!h[frac_start].freq;++frac_start);
		for(one_start=frac_start;one_start<nlevels&&!h[one_start].qfreq;++one_start);
		if(frac_start<one_start)
		{
			for(int k=frac_start;k<one_start;++k)
				++h[k].qfreq;
			//int correction=one_start-frac_start, normal_start=0;
			//for(;correction;)
			//{
			//	for(normal_start=one_start;normal_start<nlevels&&h[normal_start].qfreq<2;++normal_start);
			//	for(int k=normal_start;k<nlevels&&correction>0;++k)
			//		--h[k].qfreq, --correction;
			//}
		}

		int error=-ANS_L;//too much -> +ve error & vice versa
		for(int k=0;k<nlevels;++k)
			error+=h[k].qfreq;
		if(error>0)
		{
			while(error)
			{
				for(int k=0;k<nlevels&&error;++k)
				{
					int dec=h[k].qfreq>1;
					h[k].qfreq-=dec, error-=dec;
				}
			}
		}
		else
		{
			while(error)
			{
				for(int k=nlevels-1;k>=0&&error;--k)
					++h[k].qfreq, ++error;
			}
		}
		//int sign=1-((error<0)<<1);
		//for(int k=0;k<nlevels&&error;++k)
		//{
		//	if(sign<0||h[k].qfreq>1)
		//		h[k].qfreq-=sign, error-=sign;
		//}

		std::sort(h.begin(), h.end(), [](SortedHistInfo const &a, SortedHistInfo const &b)
		{
			return a.idx<b.idx;
		});
	}
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		histogram[k]=h[k].qfreq;
		CDF[k]=sum;
		sum+=histogram[k];
		if(k)
		{
			for(int k2=CDF[k-1];k2<CDF[k];++k2)
				CDF2sym[k2]=k-1;
		}
	}
	for(int k2=CDF[nlevels-1];k2<ANS_L;++k2)
		CDF2sym[k2]=nlevels-1;
	if(loud)
	{
#ifdef ANS_PRINT_HISTOGRAM
		printf("s\tf0\tfq\tCDF\n");
		for(int k=0;k<nlevels;++k)
		{
			if(histogram[k])
				printf("%3d\t%5d\t%04X\t%04X=%5d\n", k, h[k].freq, histogram[k], CDF[k], CDF[k]);
		}
#endif
		if(sum==ANS_L)
			printf("qfreq sum: %08X\n", sum);
		else
			printf("qfreq sum: %08X != %08X CDF ERROR\n", sum, ANS_L);
	}

	//int nlevels=1<<depth, mask=nlevels-1;
	//memset(CDF, 0, nlevels*sizeof(*CDF));//CDF is used as temp buffer
	//for(int k=0;k<imsize;++k)
	//	++CDF[buffer[k]>>bit0&mask];

	//quantize frequencies
	//for(int k=0;k<nlevels;++k)//round(freq*2^N/imsize)
	//	histogram[k]=(((long long)CDF[k]<<ANS_PROB_BITS)+(imsize>>1))/imsize;

	//normalize (quantize) frequencies
/*	int minfreq=floor_log2(imsize);
	if(minfreq>ANS_PROB_BITS)
		minfreq=(1<<(minfreq-ANS_PROB_BITS))-1;
	else
		minfreq=0;
	for(int k=0;k<nlevels;++k)
	{
		if(histogram[k])
		{
			if(histogram[k]<minfreq)
				histogram[k]=1, ++debt;
			else
				histogram[k]=(int)(((long long)histogram[k]<<ANS_PROB_BITS)/imsize);
		}
	}//*/

/*	//check that summation is still unity
	int debt=-(1<<ANS_PROB_BITS);
	for(int k=0;k<nlevels;++k)
		debt+=histogram[k];
	if(debt<0)
	{
		for(int k=0;k<nlevels&&debt;++k)
		{
			int small=histogram[k]==1;
			histogram[k]+=small;
			debt-=small;
		}
	}
	else if(debt>0)
	{
		for(int k=0;k<nlevels&&debt;++k)
	}//*/
/*	int diff=0;
	for(int k=0;k<nlevels;++k)
	{
		if(histogram[k])
		{
			int temp=((u64)histogram[k]<<ANS_PROB_BITS)/imsize;
			if(!temp)
				histogram[k]=-1, ++diff;
		}
	}
	for(int k=0;k<nlevels;++k)
	{
		if(histogram[k]==-1)
			histogram[k]=1;
		else if(histogram[k])
		{
			int dec=diff!=0;
			histogram[k]+=dec;
			diff-=dec;
			((u64)histogram[k]<<ANS_PROB_BITS)/imsize;
		}
	}//*/

	//int shift=floor_log2(imsize);
	//shift<<=(1<<shift)<imsize;//ceil_log2
	//int mag=((1LL<<(shift+31))+imsize-1)/imsize;
	//for(int k=0;k<nlevels;++k)
	//	histogram[k]=(long long)histogram[k]*mag<<shift;
}
/*inline void		tans_calc_histogram(const void *src, int imsize, int bit0, int depth, int bytestride, int *histogram)
{
	int nlevels=1<<depth, mask=nlevels-1;
	memset(histogram, 0, nlevels*sizeof(*histogram));
	if(bytestride==1)
	{
		auto buffer=(const unsigned char*)src;
		for(int k=0, ks=0;k<imsize;++k)
			++histogram[buffer[k]>>bit0&mask];
	}
	else if(bytestride==2)
	{
		auto buffer=(const unsigned short*)src;
		for(int k=0, ks=0;k<imsize;++k)
			++histogram[buffer[k]>>bit0&mask];
	}
	else if(bytestride==4)
	{
		auto buffer=(const unsigned*)src;
		for(int k=0, ks=0;k<imsize;++k)
			++histogram[buffer[k]>>bit0&mask];
	}
	else if(bytestride==8)
	{
		auto buffer=(const unsigned long long*)src;
		for(int k=0, ks=0;k<imsize;++k)
			++histogram[buffer[k]>>bit0&mask];
	}
}//*/
inline void rans_encode_start(u64 &state){state=ANS_L;}
template<typename T, int BIT0, int DEPTHMASK>inline void rans_encode(T const *srcptr, std::vector<unsigned short> &dst, u64 &state, const int *histogram, const int *CDF)//goes forward in src & dst
{
	int s=*srcptr>>BIT0&DEPTHMASK;
	//++srcptr;		//POST-INCREMENT SRCPTR IN YOUR LOOP

	MY_ASSERT(histogram[s], "Zero frequency symbol");

	if(state>=(unsigned)(histogram[s]<<(32-ANS_PROB_BITS)))//renormalize
		dst.push_back((unsigned short)state), state>>=16;
	//if(!histogram[s])//
	//{
	//	printf("freq: %d\n", histogram[s]);
	//	int k=0;
	//	scanf_s("%d", &k);
	//}
#ifdef CHECK_ANS_STATE
	long long s2=((long long)state/histogram[s]<<ANS_PROB_BITS)+state%histogram[s]+CDF[s];
	MY_ASSERT(s2>=0&&s2<0xFFFFFFFF, "s2 = %016llX\n= %08X/%04X<<16+%08X%%%04X+%04X", s2, state, histogram[s], state, histogram[s], CDF[s]);
	if(s2>0xFFFFFFFF)
		dst.push_back((unsigned short)s2), s2>>=16;
	state=(unsigned)s2;
#else
	state=(state/histogram[s]<<ANS_PROB_BITS)+state%histogram[s]+CDF[s];
#endif
}
inline void rans_encode_finish(std::vector<unsigned short> &dst, u64 &state)
{
	dst.push_back((unsigned short)state);
	dst.push_back((unsigned short)(state>>16));
	dst.push_back((unsigned short)(state>>32));
}

template<typename T>inline void rans_decode_start(const unsigned short *srcbuf, const unsigned short *&srcptr, T *dstbuf, T *&dstptr, int csize, int imsize, u64 &state)
{
	if(srcbuf)
		srcptr=srcbuf+csize;
	if(dstbuf)
		dstptr=dstbuf+imsize;

	--srcptr, state=*srcptr;
	--srcptr, state=state<<16|*srcptr;
	--srcptr, state=state<<16|*srcptr;
}
template<typename T, int BIT0, int DEPTHMASK>inline void rans_decode(const unsigned short *&srcptr, T *dstptr, u64 &state, const int *histogram, const int *CDF, const int *CDF2sym)//goes backwards in src & dst
{
	int c=state&(ANS_L-1);
	int s=CDF2sym[c];
	MY_ASSERT(s>=0&&s<=DEPTHMASK, "s=%08X", s);
	
	//--dstptr;		//PRE-DECREMENT DSTPTR IN YOUR LOOP
	*dstptr|=s<<BIT0;
	
#ifdef CHECK_ANS_STATE
	long long s2=(long long)histogram[s]*(state>>ANS_PROB_BITS)+(state&(ANS_L-1))-CDF[s];
	MY_ASSERT(s2>=0&&s2<0xFFFFFFFF, "s2=%016llX", s2);
	state=(unsigned)s2;
#else
	state=histogram[s]*(state>>ANS_PROB_BITS)+(state&(ANS_L-1))-CDF[s];
#endif

	if(state<ANS_L)
		--srcptr, state=state<<16|*srcptr;
}

//rANS usage example
#ifdef ANS_PRINT_RGB
const int print_rgb_pixel=193673, print_rgb_total=10;
#endif
inline int*		rans_rgb888_start(int *buffer, int imsize, bool loud)
{
	auto t1=__rdtsc();
	const int nlevels=1<<8, bufsize=nlevels*6+ANS_L*3;
	auto freqs=new int[bufsize];//193.5KB total,	first 3KB are histograms stored with compressed data
	auto
		h_R=freqs,				CDF_R=freqs+nlevels*3,	CDF2sym_R=freqs+nlevels*6,
		h_G=freqs+nlevels*1,	CDF_G=freqs+nlevels*4,	CDF2sym_G=freqs+nlevels*6+ANS_L,
		h_B=freqs+nlevels*2,	CDF_B=freqs+nlevels*5,	CDF2sym_B=freqs+nlevels*6+ANS_L*2;
	ans_calc_histogram(buffer, imsize,  0, 8, h_R, CDF_R, CDF2sym_R, true);
	ans_calc_histogram(buffer, imsize,  8, 8, h_G, CDF_G, CDF2sym_G, true);
	ans_calc_histogram(buffer, imsize, 16, 8, h_B, CDF_B, CDF2sym_B, true);
	auto t2=__rdtsc();
	if(loud)
	{
		//CPB_encode += CPB_start
		//CPB_decode += CPB_start
		printf("rANS freqs: %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(imsize*3));
	}
	return freqs;
}
inline void		rans_rgb888_finish(int *&freqs){delete[] freqs, freqs=nullptr;}
inline void		rans_rgb888_encode(const int *buffer, int imsize, const int *freqs, std::vector<unsigned short> &dst, int loud)
{
#ifdef ANS_PRINT_RGB
	if(loud)//
	{
		int start=print_rgb_pixel-(print_rgb_total>>1);
		if(start<0)
			start=0;
		int end=start+print_rgb_total;
		if(end>imsize)
			end=imsize;
		printf("\nk\tsize-k\t00BBGGRR:\n");
		for(int k=start;k<end;++k)
			printf("%6d\t%6d\t%08X\n", k, imsize-1-k, buffer[k]);
		printf("\n");
	}
#endif
	auto t1=__rdtsc();
	const int nlevels=1<<8;
	auto
		h_R=freqs,				CDF_R=freqs+nlevels*3,	CDF2sym_R=freqs+nlevels*6,
		h_G=freqs+nlevels*1,	CDF_G=freqs+nlevels*4,	CDF2sym_G=freqs+nlevels*6+ANS_L,
		h_B=freqs+nlevels*2,	CDF_B=freqs+nlevels*5,	CDF2sym_B=freqs+nlevels*6+ANS_L*2;

	dst.clear();
	dst.reserve(imsize*sizeof(int)/sizeof(short));
	u64 state_R=0, state_G=0, state_B=0;
	rans_encode_start(state_R);
	state_B=state_G=state_R;
	for(int k=0;k<imsize;++k)//encode loop: read & write forward (little-endian)
	{
#ifdef ANS_PRINT_RGB
		int debugpos=k>=print_rgb_pixel-(print_rgb_total>>1)&&k<print_rgb_pixel+(print_rgb_total>>1);
		if(debugpos)
		{
			printf("kb %d %08X  vvv %d vvv\n", k, buffer[k], print_rgb_pixel+(print_rgb_total>>1)-1-k);
			printf(" kc %6d R %016llX %02X", dst.size(), state_R, buffer[k]&0xFF);
		}
		rans_encode<int,  0, 0xFF>(buffer+k, dst, state_R, h_R, CDF_R);
		if(debugpos)
		{
			printf(" -> %016llX\n", state_R);
			printf(" kc %6d G %016llX %02X", dst.size(), state_G, buffer[k]>>8&0xFF);
		}
		rans_encode<int,  8, 0xFF>(buffer+k, dst, state_G, h_G, CDF_G);
		if(debugpos)
		{
			printf(" -> %016llX\n", state_G);
			printf(" kc %6d B %016llX %02X", dst.size(), state_B, buffer[k]>>16&0xFF);
		}
		rans_encode<int, 16, 0xFF>(buffer+k, dst, state_B, h_B, CDF_B);
		if(debugpos)
			printf(" -> %016llX\n\n", state_B);
#else
		rans_encode<int,  0, 0xFF>(buffer+k, dst, state_R, h_R, CDF_R);
		rans_encode<int,  8, 0xFF>(buffer+k, dst, state_G, h_G, CDF_G);
		rans_encode<int, 16, 0xFF>(buffer+k, dst, state_B, h_B, CDF_B);
#endif
	}
	rans_encode_finish(dst, state_R);
	rans_encode_finish(dst, state_G);
	rans_encode_finish(dst, state_B);
	auto t2=__rdtsc();
	if(loud)
	{
		printf("rANS encode: %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(imsize*3));
		int s0=imsize*3, sc=dst.size()*sizeof(short);
		printf("size: %d -> %d, ratio: %lf, %lf bpp\n", s0, sc, (double)s0/sc, 8.*sc/imsize);

		if(loud==2)
		{
			printf("Preview:\n");
			int kprint=dst.size()<100?dst.size():100;
			for(int k=0;k<kprint;++k)
				printf("%04X-", dst[k]);
			printf("\n");
		}
	}
}
inline void		rans_rgb888_decode(const unsigned short *src, int csize, int imsize, const int *freqs, int *buffer, bool loud)
{
	auto t1=__rdtsc();
	memset(buffer, 0, imsize*sizeof(*buffer));
	const int nlevels=1<<8;
	auto
		h_R=freqs,				CDF_R=freqs+nlevels*3,	CDF2sym_R=freqs+nlevels*6,
		h_G=freqs+nlevels*1,	CDF_G=freqs+nlevels*4,	CDF2sym_G=freqs+nlevels*6+ANS_L,
		h_B=freqs+nlevels*2,	CDF_B=freqs+nlevels*5,	CDF2sym_B=freqs+nlevels*6+ANS_L*2;
	const unsigned short *srcptr=nullptr;
	int *dstptr=nullptr;
	u64 state_R=0, state_G=0, state_B=0;
	rans_decode_start(src, srcptr, buffer, dstptr, csize, imsize, state_B);
	rans_decode_start(nullptr, srcptr, (int*)nullptr, dstptr, csize, imsize, state_G);
	rans_decode_start(nullptr, srcptr, (int*)nullptr, dstptr, csize, imsize, state_R);
	for(int k=imsize-1;k>=0&&srcptr>=src;--k)//decode loop: read & write backwards
	{
		--dstptr;
		MY_ASSERT(dstptr>=buffer&&dstptr<buffer+imsize, "kd=%d, size=", dstptr-buffer, imsize);
#ifdef ANS_PRINT_RGB
		int debugpos=k>=print_rgb_pixel-(print_rgb_total>>1)&&k<print_rgb_pixel+(print_rgb_total>>1);
		if(debugpos)
			printf(" kc%6d B %016llX", srcptr-src, state_B);
		rans_decode<int, 16, 0xFF>(srcptr, dstptr, state_B, h_B, CDF_B, CDF2sym_B);
		if(debugpos)
		{
			printf(" -> %016llX %02X\n", state_B, buffer[k]>>16&0xFF);
			printf(" kc%6d G %016llX", srcptr-src, state_G);
		}
		rans_decode<int,  8, 0xFF>(srcptr, dstptr, state_G, h_G, CDF_G, CDF2sym_G);
		if(debugpos)
		{
			printf(" -> %016llX %02X\n", state_G, buffer[k]>>8&0xFF);
			printf(" kc%6d R %016llX", srcptr-src, state_R);
		}
		rans_decode<int,  0, 0xFF>(srcptr, dstptr, state_R, h_R, CDF_R, CDF2sym_R);
		if(debugpos)
		{
			printf(" -> %016llX %02X\n", state_R, buffer[k]&0xFF);
			printf("kb %6d %08X  \t\t^^^ %d ^^^\n\n", k, buffer[k], print_rgb_pixel+(print_rgb_total>>1)-1-k);
		}
#else
		rans_decode<int, 16, 0xFF>(srcptr, dstptr, state_B, h_B, CDF_B, CDF2sym_B);
		rans_decode<int,  8, 0xFF>(srcptr, dstptr, state_G, h_G, CDF_G, CDF2sym_G);
		rans_decode<int,  0, 0xFF>(srcptr, dstptr, state_R, h_R, CDF_R, CDF2sym_R);
#endif
	}
	auto t2=__rdtsc();
	if(loud)
	{
		printf("rANS decode: %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(imsize*3));
#ifdef ANS_PRINT_RGB
		int start=print_rgb_pixel-(print_rgb_total>>1);
		if(start<0)
			start=0;
		int end=start+print_rgb_total;
		if(end>imsize)
			end=imsize;
		printf("\nk\tsize-k\t00BBGGRR:\n");
		for(int k=start;k<end;++k)
			printf("%6d\t%6d\t%08X\n", k, imsize-1-k, buffer[k]);
		printf("\n");
#endif
	}
}
/*template<int BIT0, int DEPTH>inline void tans_encode(const void *src, int imsize, const int *depths, int nchannels, int bytestride, std::string &out_data, bool loud)
{
	auto buffer=(const unsigned char*)src;
	auto t1=__rdtsc();
	auto t2=__rdtsc();
}
void			tans_decode(const void *data, void *dst, int imsize, const int *depths, int nchannels, int bytestride, bool loud)
{
}//*/

//void			rans_encode(const void *src, int imsize, int depth, int bytestride, std::string &out_data, int *out_sizes, int *out_conf, bool loud);
//void			rans_decode(const void *data, const int *sizes, const int *conf, void *dst, int imsize, int depth, int bytestride, bool loud);


void			abac3_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_conf, bool loud);
void			abac3_decode(const char *data, const int *sizes, const int *conf, short *buffer, int imsize, int depth, bool loud);

void			abac2_encode(const void *src, int imsize, int depth, int bytestride, std::string &out_data, int *out_sizes, int *out_conf, bool loud=false);
void			abac2_decode(const char *data, const int *sizes, const int *conf, void *dst, int imsize, int depth, int bytestride, bool loud=false);


#ifndef NO_SIMD
void			abac_encode_sse2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_sse2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);

void			abac_encode_avx2(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode_avx2(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);
#endif
void			abac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, bool loud=false);
void			abac_decode(const char *data, const int *sizes, short *buffer, int imsize, int depth, bool loud=false);
int				abac_estimate(const void *src, int imsize, int bitdepth, int bytestride, bool loud=false, int *sizes=nullptr);


void			ac_encode(const short *buffer, int imsize, int depth, std::string &out_data, int *out_sizes, int *out_probabilities, bool loud=false);
void			ac_decode(const char *data, const int *sizes, const int *probabilities, short *buffer, int imsize, int depth, bool loud=false);


void			ac_test_bitplane_differentiation(short *buffer, int imsize, int depth, int &dmask);
void			ac_differentiate_bitplanes(short *buffer, int imsize, int depth, int dmask);
void			ac_integrate_bitplanes(short *buffer, int imsize, int depth, int dmask);

#endif