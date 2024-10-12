#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#include"blist.h"
static const char file[]=__FILE__;


	#define ESTIMATE_SIZE
	#define DEBUG_ANS


#ifdef DEBUG_ANS
typedef struct DebugANSInfoStruct
{
	unsigned s0, state;
	int den, cdf, freq, sym, id;
} DebugANSInfo;
static void debug_enc_update(unsigned state, int den, int cdf, int freq, int sym);
static void debug_dec_update(unsigned state, int den, int cdf, int freq, int sym);
//#ifdef AC_IMPLEMENTATION
static SList states={0};
static void debug_ans_print(DebugANSInfo *info, int dec)
{
	printf("%6d state 0x%08X%s0x%08X den%8d cdf%8d freq%8d sym%8d\n",
		info->id,
		info->s0,
		dec?"<-":"->",
		info->state,
		info->den,
		info->cdf,
		info->freq,
		info->sym
	);
}
static void debug_enc_dump(DebugANSInfo *i2)
{
	debug_ans_print(i2, 1);
	printf("\n");
	for(int k=0;k<20&&states.count;++k)
	{
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states);
		debug_ans_print(i0, 0);
		STACK_POP(&states);
	}
}
static void debug_enc_update(unsigned state, int den, int cdf, int freq, int sym)
{
	unsigned s0=state;
	if(freq<=0||cdf+freq>den)
		LOG_ERROR2("ANS: invalid frequency %d", freq);
	if(!states.count)
		slist_init(&states, sizeof(DebugANSInfo), 0);
	state=state/freq*den+cdf+state%freq;//enc update
	{
		DebugANSInfo info={s0, state, den, cdf, freq, sym, (unsigned)states.count};
		STACK_PUSH(&states, &info);
	}
}
static void debug_dec_update(unsigned state, int den, int cdf, int freq, int sym)
{
	if(freq<=0||cdf+freq>den)
		LOG_ERROR2("ANS: invalid frequency %d", freq);
	if(!states.count)
		LOG_ERROR2("Nothing to decode");
	{
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states), info;
		memcpy(&info, i0, sizeof(info));

		unsigned s0=state/den*freq+state%den-cdf;//dec update

		if(info.s0!=s0||info.state!=state||info.den!=den||info.cdf!=cdf||info.freq!=freq||info.sym!=sym)
		{
			DebugANSInfo i2={s0, state, den, cdf, freq, sym, (unsigned)states.count-1};
			debug_enc_dump(&i2);
			LOG_ERROR2("Decode error  (%d decodes remaining)", info.id);
		}

		STACK_POP(&states);
	}
}
//#endif
#else
#define debug_enc_update(...)
#define debug_dec_update(...)
#endif

#if 1
#define CLAMPGRAD(PRED, N, W, GRAD)\
	do\
	{\
		int _vmin=N, _vmax=W, _G;\
		if(_vmax<_vmin)\
		{\
			_vmax=N;\
			_vmin=W;\
		}\
		_G=GRAD;\
		if(_G<_vmin)\
			_G=_vmin;\
		if(_G>_vmax)\
			_G=_vmax;\
		PRED=_G;\
	}while(0)
#else
#define CLAMPGRAD(PRED, N, W, GRAD) MEDIAN3_32(PRED, N, W, GRAD)
#endif
int c20_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
	double csizes[3]={0};
	ptrdiff_t csize_actual=0;
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize;
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<1)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+16);
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
		srcbuf[srcsize]=0;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	if(srcsize<=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	unsigned char *image=0;
	if(fwd)//encode
	{
		if(*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++ - '0';
		if(iw<1||*srcptr++ != ' ')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++ - '0';
		if(ih<1||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		int temp=0;
		while((unsigned)(*srcptr-'0')<10)
			temp=10*temp+*srcptr++ - '0';
		if(temp!=255||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		image=srcptr;


		int hperiod=0x8000/iw, hcount=(ih+hperiod-1)/hperiod+1;
		int hsize=hcount*sizeof(short[3][256]);
		unsigned short *hist=(unsigned short*)malloc(hsize);
		int denssize=hcount*sizeof(short[3]);
		unsigned short *dens=(unsigned short*)malloc(denssize);
		if(!hist||!dens)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(hist, 0, hsize);
		memset(dens, 0, denssize);
		unsigned short *curr_dens[]=
		{
			dens+hcount*0LL+0,//Y
			dens+hcount*0LL+1,
			dens+hcount*1LL+0,//U
			dens+hcount*1LL+1,
			dens+hcount*2LL+0,//V
			dens+hcount*2LL+1,
		};
		unsigned short *curr_hists[]=
		{
			hist+256*(hcount*0LL+0),//Y
			hist+256*(hcount*0LL+1),
			hist+256*(hcount*1LL+0),//U
			hist+256*(hcount*1LL+1),
			hist+256*(hcount*2LL+0),//V
			hist+256*(hcount*2LL+1),
		};

		ALIGN(16) short
			errorsY[8]={0},
			errorsU[8]={0},
			errorsV[8]={0};
		__m128i getY=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,    8,  7,    6,  5,    4,  3,    2,  1,    0
			-1, -1, -1, -1, -1, -1, -1, 12+1, -1,  9+1, -1,  6+1, -1,  3+1, -1,  0+1
		);
		__m128i getU=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,    8,  7,    6,  5,    4,  3,    2,  1,    0
			-1, -1, -1, -1, -1, -1, -1, 12+0, -1,  9+2, -1,  6+2, -1,  3+2, -1,  0+2
		);
		__m128i getV=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,    8,  7,    6,  5,    4,  3,    2,  1,    0
			-1, -1, -1, -1, -1, -1, -1, 12+0, -1,  9+0, -1,  6+0, -1,  3+0, -1,  0+0
		);
		__m128i mmin=_mm_set1_epi16(0);
		__m128i mmax=_mm_set1_epi16(255);
		for(int ky=ih-1, ky3=0;ky>=100;--ky)		//predict backwards (inplace)
		{
			int kx=iw-1;
			unsigned char *Nptr	=image+3*(iw*(ky-1)+kx);
			unsigned char *ptr	=image+3*(iw*(ky+0)+kx);
			for(;kx>iw-5;--kx)
			{
				int NW, N, W, curr;
				int offset, pred, sym;

				NW	=Nptr	[1-1*3]-128;
				N	=Nptr	[1+0*3]-128;
				W	=ptr	[1-1*3]-128;
				curr	=ptr	[1+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				sym=(curr-pred)<<(32-8);//signed modular arithmetic
				sym=sym>>(31-8)^sym>>31;//pack sign
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;

				NW	=Nptr	[2-1*3]-Nptr	[1-1*3];
				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				W	=ptr	[2-1*3]-ptr	[1-1*3];
				curr	=ptr	[2+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];

				NW	=Nptr	[0-1*3]-Nptr	[1-1*3];
				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				W	=ptr	[0-1*3]-ptr	[1-1*3];
				curr	=ptr	[0+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];

				Nptr	-=3;
				ptr	-=3;
			}
			for(;kx>4;kx-=5)	//fast predict
			{
				__m128i mNW	=_mm_loadu_si128((__m128i*)(Nptr-1*3));
				__m128i mN	=_mm_loadu_si128((__m128i*)(Nptr+0*3));
				__m128i mW	=_mm_loadu_si128((__m128i*)(ptr	-1*3));
				__m128i mcurr	=_mm_loadu_si128((__m128i*)(ptr	+0*3));
				__m128i mNW0	=_mm_shuffle_epi8(mNW, getY);
				__m128i mNW1	=_mm_shuffle_epi8(mNW, getU);
				__m128i mNW2	=_mm_shuffle_epi8(mNW, getV);
				__m128i mN0	=_mm_shuffle_epi8(mN, getY);
				__m128i mN1	=_mm_shuffle_epi8(mN, getU);
				__m128i mN2	=_mm_shuffle_epi8(mN, getV);
				__m128i mW0	=_mm_shuffle_epi8(mW, getY);
				__m128i mW1	=_mm_shuffle_epi8(mW, getU);
				__m128i mW2	=_mm_shuffle_epi8(mW, getV);
				__m128i mcurr0	=_mm_shuffle_epi8(mcurr, getY);
				__m128i mcurr1	=_mm_shuffle_epi8(mcurr, getU);
				__m128i mcurr2	=_mm_shuffle_epi8(mcurr, getV);
				
				mNW1=_mm_sub_epi16(mNW1, mNW0);
				mNW2=_mm_sub_epi16(mNW2, mNW0);
				mN1=_mm_sub_epi16(mN1, mN0);
				mN2=_mm_sub_epi16(mN2, mN0);
				mW1=_mm_sub_epi16(mW1, mW0);
				mW2=_mm_sub_epi16(mW2, mW0);
				__m128i vmin0=_mm_min_epi16(mN0, mW0);
				__m128i vmin1=_mm_min_epi16(mN1, mW1);
				__m128i vmin2=_mm_min_epi16(mN2, mW2);
				__m128i vmax0=_mm_max_epi16(mN0, mW0);
				__m128i vmax1=_mm_max_epi16(mN1, mW1);
				__m128i vmax2=_mm_max_epi16(mN2, mW2);
				__m128i mp0=_mm_add_epi16(mN0, mW0);
				__m128i mp1=_mm_add_epi16(mN1, mW1);
				__m128i mp2=_mm_add_epi16(mN2, mW2);
				mp0=_mm_sub_epi16(mp0, mNW0);
				mp1=_mm_sub_epi16(mp1, mNW1);
				mp2=_mm_sub_epi16(mp2, mNW2);
				mp0=_mm_max_epi16(mp0, vmin0);
				mp1=_mm_max_epi16(mp1, vmin1);
				mp2=_mm_max_epi16(mp2, vmin2);
				mp0=_mm_min_epi16(mp0, vmax0);
				mp1=_mm_min_epi16(mp1, vmax1);
				mp2=_mm_min_epi16(mp2, vmax2);
				mp1=_mm_add_epi16(mp1, mcurr0);
				mp2=_mm_add_epi16(mp2, mcurr0);
				mp1=_mm_max_epi16(mp1, mmin);
				mp2=_mm_max_epi16(mp2, mmin);
				mp1=_mm_min_epi16(mp1, mmax);
				mp2=_mm_min_epi16(mp2, mmax);
				mp0=_mm_sub_epi16(mcurr0, mp0);
				mp1=_mm_sub_epi16(mcurr1, mp1);
				mp2=_mm_sub_epi16(mcurr2, mp2);
				mp0=_mm_slli_epi16(mp0, 16-8);
				mp1=_mm_slli_epi16(mp1, 16-8);
				mp2=_mm_slli_epi16(mp2, 16-8);
				mp0=_mm_xor_si128(_mm_srli_epi16(mp0, 15), _mm_srli_epi16(mp0, 7));//pack sign
				mp1=_mm_xor_si128(_mm_srli_epi16(mp1, 15), _mm_srli_epi16(mp1, 7));
				mp2=_mm_xor_si128(_mm_srli_epi16(mp2, 15), _mm_srli_epi16(mp2, 7));
				_mm_store_si128((__m128i*)errorsY, mp0);
				_mm_store_si128((__m128i*)errorsU, mp1);
				_mm_store_si128((__m128i*)errorsV, mp2);
				
				ptr[1+0*3]=(unsigned char)errorsY[0];
				ptr[2+0*3]=(unsigned char)errorsU[0];
				ptr[0+0*3]=(unsigned char)errorsV[0];
				ptr[1+1*3]=(unsigned char)errorsY[1];
				ptr[2+1*3]=(unsigned char)errorsU[1];
				ptr[0+1*3]=(unsigned char)errorsV[1];
				ptr[1+2*3]=(unsigned char)errorsY[2];
				ptr[2+2*3]=(unsigned char)errorsU[2];
				ptr[0+2*3]=(unsigned char)errorsV[2];
				ptr[1+3*3]=(unsigned char)errorsY[3];
				ptr[2+3*3]=(unsigned char)errorsU[3];
				ptr[0+3*3]=(unsigned char)errorsV[3];
				ptr[1+4*3]=(unsigned char)errorsY[4];
				ptr[2+4*3]=(unsigned char)errorsU[4];
				ptr[0+4*3]=(unsigned char)errorsV[4];
				++curr_hists	[2*0+0][errorsY[0]];
				++curr_hists	[2*0+1][errorsY[0]];
				++curr_hists	[2*1+0][errorsU[0]];
				++curr_hists	[2*1+1][errorsU[0]];
				++curr_hists	[2*2+0][errorsV[0]];
				++curr_hists	[2*2+1][errorsV[0]];
				++curr_hists	[2*0+0][errorsY[1]];
				++curr_hists	[2*0+1][errorsY[1]];
				++curr_hists	[2*1+0][errorsU[1]];
				++curr_hists	[2*1+1][errorsU[1]];
				++curr_hists	[2*2+0][errorsV[1]];
				++curr_hists	[2*2+1][errorsV[1]];
				++curr_hists	[2*0+0][errorsY[2]];
				++curr_hists	[2*0+1][errorsY[2]];
				++curr_hists	[2*1+0][errorsU[2]];
				++curr_hists	[2*1+1][errorsU[2]];
				++curr_hists	[2*2+0][errorsV[2]];
				++curr_hists	[2*2+1][errorsV[2]];
				++curr_hists	[2*0+0][errorsY[3]];
				++curr_hists	[2*0+1][errorsY[3]];
				++curr_hists	[2*1+0][errorsU[3]];
				++curr_hists	[2*1+1][errorsU[3]];
				++curr_hists	[2*2+0][errorsV[3]];
				++curr_hists	[2*2+1][errorsV[3]];
				++curr_hists	[2*0+0][errorsY[4]];
				++curr_hists	[2*0+1][errorsY[4]];
				++curr_hists	[2*1+0][errorsU[4]];
				++curr_hists	[2*1+1][errorsU[4]];
				++curr_hists	[2*2+0][errorsV[4]];
				++curr_hists	[2*2+1][errorsV[4]];
				curr_dens	[2*0+0][0]+=4;
				curr_dens	[2*0+1][0]+=4;
				curr_dens	[2*1+0][0]+=4;
				curr_dens	[2*1+1][0]+=4;
				curr_dens	[2*2+0][0]+=4;
				curr_dens	[2*2+1][0]+=4;

				Nptr	-=3*5;
				ptr	-=3*5;
			}
			for(;kx>0;--kx)
			{
				int NW, N, W, curr;
				int offset, pred, sym;

				NW	=Nptr	[1-1*3]-128;
				N	=Nptr	[1+0*3]-128;
				W	=ptr	[1-1*3]-128;
				curr	=ptr	[1+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				sym=(curr-pred)<<(32-8);//signed modular arithmetic
				sym=sym>>(31-8)^sym>>31;//pack sign
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;

				NW	=Nptr	[2-1*3]-Nptr	[1-1*3];
				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				W	=ptr	[2-1*3]-ptr	[1-1*3];
				curr	=ptr	[2+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];

				NW	=Nptr	[0-1*3]-Nptr	[1-1*3];
				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				W	=ptr	[0-1*3]-ptr	[1-1*3];
				curr	=ptr	[0+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];

				Nptr	-=3;
				ptr	-=3;
			}
			{
				int N, curr;
				int offset, pred, sym;

				N	=Nptr	[1+0*3]-128;
				curr	=ptr	[1+0*3]-128;
				pred=N;
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;

				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				curr	=ptr	[2+0*3]-128;
				pred=N;
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];

				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				curr	=ptr	[0+0*3]-128;
				pred=N;
				pred+=offset;
				CLAMP2(pred, -128, 127);
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];

				Nptr	-=3;
				ptr	-=3;
			}
			if(!(ky%hperiod))
			{
				curr_hists	[2*0+(ky3&1)]+=256*2;//skip 2 histograms (current and next)
				curr_hists	[2*1+(ky3&1)]+=256*2;
				curr_hists	[2*2+(ky3&1)]+=256*2;
				curr_dens	[2*0+(ky3&1)]+=2;
				curr_dens	[2*1+(ky3&1)]+=2;
				curr_dens	[2*2+(ky3&1)]+=2;
				++ky3;
			}
		}
		{
			int kx=iw-1;
		//	unsigned char *Nptr	=image+3*(iw*(0-1)+kx);
			unsigned char *ptr	=image+3*(iw*(0+0)+kx);
			for(;kx>=0;--kx)
			{
				int W, curr;
				int offset, pred;

				W	=kx?ptr	[1-1*3]-128:0;
				curr	=ptr	[1+0*3]-128;
				pred=W;
				++curr_hists	[2*0+0][(curr-pred+128)&255];
				++curr_hists	[2*0+1][(curr-pred+128)&255];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;

				W	=kx?ptr	[2-1*3]-ptr	[1-1*3]:0;
				curr	=ptr	[2+0*3]-128;
				pred=W;
				pred+=offset;
				CLAMP2(pred, -128, 127);
				++curr_hists	[2*1+0][(curr-pred+128)&255];
				++curr_hists	[2*1+1][(curr-pred+128)&255];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];

				W	=kx?ptr	[0-1*3]-ptr	[1-1*3]:0;
				curr	=ptr	[0+0*3]-128;
				pred=W;
				pred+=offset;
				CLAMP2(pred, -128, 127);
				++curr_hists	[2*2+0][(curr-pred+128)&255];
				++curr_hists	[2*2+1][(curr-pred+128)&255];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];

			//	--Nptr;
				--ptr;
			}
		}
		ptrdiff_t dstbufsize=4LL*iw*ih;
		unsigned short *dstbuf=(unsigned short*)malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		unsigned short *dstptr=dstbuf;
		curr_dens[0]=dens+hcount*0LL+0;//Y
		curr_dens[1]=dens+hcount*0LL+1;
		curr_dens[2]=dens+hcount*1LL+0;//U
		curr_dens[3]=dens+hcount*1LL+1;
		curr_dens[4]=dens+hcount*2LL+0;//V
		curr_dens[5]=dens+hcount*2LL+1;
		curr_hists[0]=hist+256*(hcount*0LL+0);//Y
		curr_hists[1]=hist+256*(hcount*0LL+1);
		curr_hists[2]=hist+256*(hcount*1LL+0);//U
		curr_hists[3]=hist+256*(hcount*1LL+1);
		curr_hists[4]=hist+256*(hcount*2LL+0);//V
		curr_hists[5]=hist+256*(hcount*2LL+1);
		unsigned state=0x10000;
		for(int ky=ih-1, ky3=0;ky>=0;--ky)	//encode
		{
			int kx=iw-1;
			unsigned char *ptr=image+3*((size_t)iw*ky+kx);
			for(;kx>=0;--kx)
			{
				int yuv[]=
				{
					ptr[1],//Y
					ptr[2],//U
					ptr[0],//V
				};
				int sym, den, cdf, freq;

				//enc V
#define KC 2
				sym=yuv[KC];
				if(!curr_hists[2*KC+0][sym])LOG_ERROR("");//
				if(!curr_hists[2*KC+1][sym])LOG_ERROR("");//
				--curr_hists	[2*KC+0][sym];
				--curr_hists	[2*KC+1][sym];
				--curr_dens	[2*KC+0][0];
				--curr_dens	[2*KC+1][0];
				den=curr_dens[2*KC+0][0]+curr_dens[2*KC+1][0]+256;
				cdf=0;
				for(int k=0;;++k)
				{
					freq=curr_hists[2*KC+0][k]+curr_hists[2*KC+1][k]+1;
					if(k>=sym)
						break;
					cdf+=freq;
				}
				if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
				{
					*dstptr++=(unsigned short)state;
					state>>=16;
				}
				debug_enc_update(state, den, cdf, freq, sym);
				state=state/freq*den+cdf+state%freq;
#ifdef ESTIMATE_SIZE
				csizes[KC]-=log2((double)freq/den);
#endif
#undef  KC
				
				//enc U
#define KC 1
				sym=yuv[KC];
				if(!curr_hists[2*KC+0][sym])LOG_ERROR("");//
				if(!curr_hists[2*KC+1][sym])LOG_ERROR("");//
				--curr_hists	[2*KC+0][sym];
				--curr_hists	[2*KC+1][sym];
				--curr_dens	[2*KC+0][0];
				--curr_dens	[2*KC+1][0];
				den=curr_dens[2*KC+0][0]+curr_dens[2*KC+1][0]+256;
				cdf=0;
				for(int k=0;;++k)
				{
					freq=curr_hists[2*KC+0][k]+curr_hists[2*KC+1][k]+1;
					if(k>=sym)
						break;
					cdf+=freq;
				}
				if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
				{
					*dstptr++=(unsigned short)state;
					state>>=16;
				}
				debug_enc_update(state, den, cdf, freq, sym);
				state=state/freq*den+cdf+state%freq;
#ifdef ESTIMATE_SIZE
				csizes[KC]-=log2((double)freq/den);
#endif
#undef  KC
				
				//enc Y
#define KC 0
				sym=yuv[KC];
				if(!curr_hists[2*KC+0][sym])LOG_ERROR("");//
				if(!curr_hists[2*KC+1][sym])LOG_ERROR("");//
				--curr_hists	[2*KC+0][sym];
				--curr_hists	[2*KC+1][sym];
				--curr_dens	[2*KC+0][0];
				--curr_dens	[2*KC+1][0];
				den=curr_dens[2*KC+0][0]+curr_dens[2*KC+1][0]+256;
				cdf=0;
				for(int k=0;;++k)
				{
					freq=curr_hists[2*KC+0][k]+curr_hists[2*KC+1][k]+1;
					if(k>=sym)
						break;
					cdf+=freq;
				}
				if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
				{
					*dstptr++=(unsigned short)state;
					state>>=16;
				}
				debug_enc_update(state, den, cdf, freq, sym);
				state=state/freq*den+cdf+state%freq;
#ifdef ESTIMATE_SIZE
				csizes[KC]-=log2((double)freq/den);
#endif
#undef  KC

				ptr-=3;
			}
			if(!(ky%hperiod))
			{
				curr_hists	[2*0+(ky3&1)]+=256*2;//skip 2 histograms (current and next)
				curr_hists	[2*1+(ky3&1)]+=256*2;
				curr_hists	[2*2+(ky3&1)]+=256*2;
				curr_dens	[2*0+(ky3&1)]+=2;
				curr_dens	[2*1+(ky3&1)]+=2;
				curr_dens	[2*2+(ky3&1)]+=2;
				++ky3;
			}
		}
		*dstptr++=(unsigned short)state;
		*dstptr++=(unsigned short)(state>>16);

		{//save dstptr-dstbuf
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Canno open \"%s\" for writing", dstfn);
				return 1;
			}
			ptrdiff_t nemitts=dstptr-dstbuf;
			fwrite("CH", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(&nemitts, 1, 4, fdst);
			fwrite(dstbuf, 2, nemitts, fdst);
			fclose(fdst);
#ifdef ESTIMATE_SIZE
			csize_actual=nemitts*2LL+4LL*3+2;
#endif
		}
		free(dstbuf);
		free(hist);
		free(dens);
	}
	else//decode
	{
		int nemitts=0;
		if(srcptr+4*3>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
		memcpy(&nemitts, srcptr, 4); srcptr+=4;
		if(iw<1||ih<1||srcend-srcptr<nemitts*2LL)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		unsigned short *loptr=(unsigned short*)srcptr, *ptr=loptr+nemitts-1;

		ptrdiff_t dstbufsize=3LL*iw*ih;
		unsigned char *dstbuf=(unsigned char*)malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		image=dstbuf;


		int hperiod=0x8000/iw;
		int hsize=sizeof(short[3][256*2]);
		unsigned short *hist=(unsigned short*)malloc(hsize);
		int psize=(iw+16LL)*sizeof(short[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
		short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
		if(!hist||!pixels)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(hist, 0, hsize);
		memset(pixels, 0, psize);
		unsigned short *hists[]=
		{
			hist+256*0,
			hist+256*1,
			hist+256*2,
			hist+256*3,
			hist+256*4,
			hist+256*5,
		};
		unsigned short dens[6]={0};
		unsigned state=ptr[0]<<16|ptr[-1];
		ptr-=2;
		for(int ky=0, idx=0, ky2=hperiod, ky3=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&1)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&1)+8LL)*4*2,
			};
			ALIGN(16) short preds[8]={0};
			for(int kx=0;kx<iw;++kx, idx+=3)//decode
			{
				short *curr=rows[0];
				__m128i mNW	=_mm_load_si128((__m128i*)(rows[1]-1*4*2));
				__m128i mN	=_mm_load_si128((__m128i*)(rows[1]+0*4*2));
				__m128i mW	=_mm_load_si128((__m128i*)(rows[0]-1*4*2));
				__m128i vmin=_mm_min_epi16(mN, mW);
				__m128i vmax=_mm_max_epi16(mN, mW);
				__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
				mp=_mm_max_epi16(mp, vmin);
				mp=_mm_min_epi16(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);

				int den, quotient, remainder, cdf, freq, sym;
				
				//dec Y
#define KC 0
				den=dens[2*KC+0]+dens[2*KC+1]+256;
				quotient=state/den;
				remainder=state%den;
				sym=0;
				cdf=0;
				for(;;)
				{
					unsigned cdf2;

					freq=hists[2*KC+0][sym]+hists[2*KC+1][sym]+1;
					cdf2=cdf+freq;
					if(cdf2>(unsigned)remainder)
						break;
					cdf=cdf2;
					++sym;
				}
				debug_dec_update(state, den, cdf, freq, sym);
				state=quotient*freq+remainder-cdf;
				if(state<(unsigned)den)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				sym=sym>>1^-(sym&1);
				curr[KC]=(sym+preds[KC])<<24>>24;
#undef  KC

				preds[1]+=curr[0];
				CLAMP2(preds[1], -128, 127);
				preds[2]+=curr[0];
				CLAMP2(preds[2], -128, 127);
				
				//dec U
#define KC 1
				den=dens[2*KC+0]+dens[2*KC+1]+256;
				quotient=state/den;
				remainder=state%den;
				sym=0;
				cdf=0;
				for(;;)
				{
					unsigned cdf2;

					freq=hists[2*KC+0][sym]+hists[2*KC+1][sym]+1;
					cdf2=cdf+freq;
					if(cdf2>(unsigned)remainder)
						break;
					cdf=cdf2;
					++sym;
				}
				debug_dec_update(state, den, cdf, freq, sym);
				state=quotient*freq+remainder-cdf;
				if(state<(unsigned)den)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				sym=sym>>1^-(sym&1);
				curr[KC]=(sym+preds[KC])<<24>>24;
#undef  KC
				
				//dec V
#define KC 2
				den=dens[2*KC+0]+dens[2*KC+1]+256;
				quotient=state/den;
				remainder=state%den;
				sym=0;
				cdf=0;
				for(;;)
				{
					unsigned cdf2;

					freq=hists[2*KC+0][sym]+hists[2*KC+1][sym]+1;
					cdf2=cdf+freq;
					if(cdf2>(unsigned)remainder)
						break;
					cdf=cdf2;
					++sym;
				}
				debug_dec_update(state, den, cdf, freq, sym);
				state=quotient*freq+remainder-cdf;
				if(state<(unsigned)den)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				curr[KC]=(sym+preds[KC])<<24>>24;
#undef  KC

			//	curr[4+0]=rows[0]-preds[0];
			//	curr[4+1]=rows[1]-preds[1];
			//	curr[4+2]=rows[2]-preds[2];
				image[idx+1]=(unsigned char)(curr[0]+128);
				image[idx+2]=(unsigned char)(curr[1]+128);
				image[idx+0]=(unsigned char)(curr[2]+128);
				curr[1]-=curr[0];
				curr[2]-=curr[0];

				rows[0]+=4*2;
				rows[1]+=4*2;
			}
			if(ky==ky2)
			{
				memset(hists+256*(2LL*0+(ky3&1)), 0, sizeof(short[256]));
				memset(hists+256*(2LL*1+(ky3&1)), 0, sizeof(short[256]));
				memset(hists+256*(2LL*2+(ky3&1)), 0, sizeof(short[256]));
				dens[2LL*0+(ky3&1)]=0;
				dens[2LL*1+(ky3&1)]=0;
				dens[2LL*2+(ky3&1)]=0;
				ky2+=hperiod;
				++ky3;
			}
		}
		free(hist);
		_mm_free(pixels);
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, 3LL*iw*ih, fdst);
			fclose(fdst);
		}
	}
	free(srcbuf);
#ifdef ESTIMATE_SIZE
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
	if(fwd)
	{
		for(int kc=0;kc<3;++kc)
			printf("%c%12.2lf\n", "YUV"[kc], csizes[kc]);
		printf(" %12td   /%12td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	}
	printf("%c%12lf sec %12lf MB/s  %12lld cycles %12lf C/B\n",
		'D'+fwd,
		elapsed,
		usize/(elapsed*1024*1024),
		cycles,
		(double)cycles/usize
	);
#endif
	return 0;
}