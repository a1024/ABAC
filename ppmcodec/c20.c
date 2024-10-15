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
	#define ENABLE_GUIDE
//	#define DISABLE_PRED


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
	printf("%6d state 0x%08X%s0x%08X   den%8d   cdf%8d   freq%8d   sym%8d\n",
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

		if(info.s0!=s0||info.state!=state||
			info.den!=den||info.cdf!=cdf||info.freq!=freq||info.sym!=sym
		)
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

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
		LOG_ERROR("");
}
#else
#define guide_save(...)
#define guide_check(...)
#endif

#if 0
static void deb(unsigned short **curr_hists, unsigned short **curr_dens)
{
	for(int kc=0;kc<6;++kc)
	{
		int sum=0;
		for(int ks=0;ks<256;++ks)
			sum+=curr_hists[kc][ks];
		if(sum!=curr_dens[kc][0])
			LOG_ERROR("");
	}
}
static void deb2(unsigned short **curr_hists, unsigned short *curr_dens)
{
	for(int kc=0;kc<6;++kc)
	{
		int sum=0;
		for(int ks=0;ks<256;++ks)
			sum+=curr_hists[kc][ks];
		if(sum!=curr_dens[kc])
			LOG_ERROR("");
	}
}
#define DEB() deb(curr_hists, curr_dens)
#define DEB2() deb2(hists, dens)
#else
#define DEB()
#define DEB2()
#endif
#if 0
#define DEB\
	do\
	{\
		for(int kc=0;kc<6;++kc)\
		{\
			int sum=0;\
			for(int ks=0;ks<256;++ks)\
				sum+=curr_hists[kc][ks];\
			if(sum!=curr_dens[kc][0])\
				LOG_ERROR("");\
		}\
	}while(0)
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
	//unsigned char *image=0;
	//ptrdiff_t dstbufsize=0;
	//unsigned char *dstbuf;
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
		//image=srcptr;

		//guide_save(image, iw, ih);

		//dstbufsize=4LL*iw*ih;
		//dstbuf=(unsigned char*)malloc(dstbufsize);
		//if(!dstbuf)
		//{
		//	LOG_ERROR("Alloc error");
		//	return 1;
		//}
		//image=dstbuf;
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

		//dstbufsize=3LL*iw*ih;
		//dstbuf=(unsigned char*)malloc(dstbufsize);
		//if(!dstbuf)
		//{
		//	LOG_ERROR("Alloc error");
		//	return 1;
		//}
		//image=dstbuf;
	}
	
	int psize=(iw+32LL)*sizeof(short[4*4]);//4 padded rows * 4 channels max
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m256i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	
	unsigned short stats[3][256]={0};
	FILLMEM((unsigned short*)stats, 0x8000, sizeof(stats), sizeof(short));

	ALIGN(32) unsigned short low[48]={0}, range[48]={0};
	for(int k=0;k<48;++k)
		range[k]=0xFFFF;
	__m256i mlow0=_mm256_load_si256((__m256i*)low+0);
	__m256i mlow1=_mm256_load_si256((__m256i*)low+1);
	__m256i mlow2=_mm256_load_si256((__m256i*)low+2);
	__m256i mrange0=_mm256_load_si256((__m256i*)range+0);
	__m256i mrange1=_mm256_load_si256((__m256i*)range+1);
	__m256i mrange2=_mm256_load_si256((__m256i*)range+2);
	if(fwd)
	{
		int dstbufsize=iw*ih;
		unsigned char *dstbufs[48]={0}, *dstptrs[48]={0};
		for(int k=0;k<48;++k)
		{
			dstbufs[k]=(unsigned char*)malloc(dstbufsize);
			if(!dstbufs[k])
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			dstptrs[k]=dstbufs[k];
		}
		unsigned char *image=srcptr;
		guide_save(image, iw, ih);
		for(int ky=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
			ALIGN(32) unsigned char tidx[48]={0};
			ALIGN(32) short syms[48]={0};
			ALIGN(32) unsigned short bits[48]={0};
			int idx=3*iw*ky;
			for(int kx=0;kx<iw;kx+=16, idx+=48)
			{
				short
					*NW	=rows[1]-1*4,
					*N	=rows[1]+0*4,
					*W	=rows[0]-1*4,
					*curr	=rows[0]+0*4;

				int remainder=iw-kx;
				if(remainder>16)
					remainder=16;
				remainder*=3;
				for(int k=0;k<remainder;k+=3)//TODO
				{
					curr[k+0]=image[idx+k+1];//Y
					curr[k+1]=image[idx+k+2];//U
					curr[k+2]=image[idx+k+0];//V
				}

				__m256i himask=_mm256_set1_epi16(0xFF00);
				__m256i lomask=_mm256_set1_epi16(0x00FF);
				__m256i mask8=_mm256_set1_epi16(8);
				__m256i probmin=_mm256_set1_epi16(0x0100);
				__m256i probmax=himask;
				__m256i ones=_mm256_set1_epi16(-1);
				__m256i one=_mm256_set1_epi16(1);
				__m256i bitmask=_mm256_set1_epi16(0x80);

				__m256i msym0=_mm256_load_si256((__m256i*)syms+0);
				__m256i msym1=_mm256_load_si256((__m256i*)syms+1);
				__m256i msym2=_mm256_load_si256((__m256i*)syms+2);
				memset(tidx, 1, 48);
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
				for(int kb=7;kb>=0;--kb)
				{
					__m256i mprob0=_mm256_set_epi16(
						stats[0][tidx[16*0+15]], stats[0][tidx[16*0+14]], stats[0][tidx[16*0+13]], stats[0][tidx[16*0+12]],
						stats[0][tidx[16*0+11]], stats[0][tidx[16*0+10]], stats[0][tidx[16*0+ 9]], stats[0][tidx[16*0+ 8]],
						stats[0][tidx[16*0+ 7]], stats[0][tidx[16*0+ 6]], stats[0][tidx[16*0+ 5]], stats[0][tidx[16*0+ 4]],
						stats[0][tidx[16*0+ 3]], stats[0][tidx[16*0+ 2]], stats[0][tidx[16*0+ 1]], stats[0][tidx[16*0+ 0]]
					);
					__m256i mprob1=_mm256_set_epi16(
						stats[1][tidx[16*1+15]], stats[1][tidx[16*1+14]], stats[1][tidx[16*1+13]], stats[1][tidx[16*1+12]],
						stats[1][tidx[16*1+11]], stats[1][tidx[16*1+10]], stats[1][tidx[16*1+ 9]], stats[1][tidx[16*1+ 8]],
						stats[1][tidx[16*1+ 7]], stats[1][tidx[16*1+ 6]], stats[1][tidx[16*1+ 5]], stats[1][tidx[16*1+ 4]],
						stats[1][tidx[16*1+ 3]], stats[1][tidx[16*1+ 2]], stats[1][tidx[16*1+ 1]], stats[1][tidx[16*1+ 0]]
					);
					__m256i mprob2=_mm256_set_epi16(
						stats[2][tidx[16*2+15]], stats[2][tidx[16*2+14]], stats[2][tidx[16*2+13]], stats[2][tidx[16*2+12]],
						stats[2][tidx[16*2+11]], stats[2][tidx[16*2+10]], stats[2][tidx[16*2+ 9]], stats[2][tidx[16*2+ 8]],
						stats[2][tidx[16*2+ 7]], stats[2][tidx[16*2+ 6]], stats[2][tidx[16*2+ 5]], stats[2][tidx[16*2+ 4]],
						stats[2][tidx[16*2+ 3]], stats[2][tidx[16*2+ 2]], stats[2][tidx[16*2+ 1]], stats[2][tidx[16*2+ 0]]
					);
					mprob0=_mm256_max_epu16(mprob0, probmin);
					mprob1=_mm256_max_epu16(mprob1, probmin);
					mprob2=_mm256_max_epu16(mprob2, probmin);
					mprob0=_mm256_min_epu16(mprob0, probmax);
					mprob1=_mm256_min_epu16(mprob1, probmax);
					mprob2=_mm256_min_epu16(mprob2, probmax);

					__m256i mcond0=_mm256_and_si256(mrange0, himask);
					__m256i mcond1=_mm256_and_si256(mrange1, himask);
					__m256i mcond2=_mm256_and_si256(mrange2, himask);
					mcond0=_mm256_cmpeq_epi16(mcond0, _mm256_setzero_si256());
					mcond1=_mm256_cmpeq_epi16(mcond1, _mm256_setzero_si256());
					mcond2=_mm256_cmpeq_epi16(mcond2, _mm256_setzero_si256());
					int cond0=_mm256_movemask_epi8(mcond0);
					int cond1=_mm256_movemask_epi8(mcond1);
					int cond2=_mm256_movemask_epi8(mcond2);
					if(cond0&3<<2* 0)*dstptrs[16*0+ 0]++=((unsigned char*)low)[(16*0+ 0)<<1|1];//i hope it's {TEST, CMOV, INC}
					if(cond0&3<<2* 1)*dstptrs[16*0+ 1]++=((unsigned char*)low)[(16*0+ 1)<<1|1];
					if(cond0&3<<2* 2)*dstptrs[16*0+ 2]++=((unsigned char*)low)[(16*0+ 2)<<1|1];
					if(cond0&3<<2* 3)*dstptrs[16*0+ 3]++=((unsigned char*)low)[(16*0+ 3)<<1|1];
					if(cond0&3<<2* 4)*dstptrs[16*0+ 4]++=((unsigned char*)low)[(16*0+ 4)<<1|1];
					if(cond0&3<<2* 5)*dstptrs[16*0+ 5]++=((unsigned char*)low)[(16*0+ 5)<<1|1];
					if(cond0&3<<2* 6)*dstptrs[16*0+ 6]++=((unsigned char*)low)[(16*0+ 6)<<1|1];
					if(cond0&3<<2* 7)*dstptrs[16*0+ 7]++=((unsigned char*)low)[(16*0+ 7)<<1|1];
					if(cond0&3<<2* 8)*dstptrs[16*0+ 8]++=((unsigned char*)low)[(16*0+ 8)<<1|1];
					if(cond0&3<<2* 9)*dstptrs[16*0+ 9]++=((unsigned char*)low)[(16*0+ 9)<<1|1];
					if(cond0&3<<2*10)*dstptrs[16*0+10]++=((unsigned char*)low)[(16*0+10)<<1|1];
					if(cond0&3<<2*11)*dstptrs[16*0+11]++=((unsigned char*)low)[(16*0+11)<<1|1];
					if(cond0&3<<2*12)*dstptrs[16*0+12]++=((unsigned char*)low)[(16*0+12)<<1|1];
					if(cond0&3<<2*13)*dstptrs[16*0+13]++=((unsigned char*)low)[(16*0+13)<<1|1];
					if(cond0&3<<2*14)*dstptrs[16*0+14]++=((unsigned char*)low)[(16*0+14)<<1|1];
					if(cond0&3<<2*15)*dstptrs[16*0+15]++=((unsigned char*)low)[(16*0+15)<<1|1];
					if(cond1&3<<2* 0)*dstptrs[16*1+ 0]++=((unsigned char*)low)[(16*1+ 0)<<1|1];
					if(cond1&3<<2* 1)*dstptrs[16*1+ 1]++=((unsigned char*)low)[(16*1+ 1)<<1|1];
					if(cond1&3<<2* 2)*dstptrs[16*1+ 2]++=((unsigned char*)low)[(16*1+ 2)<<1|1];
					if(cond1&3<<2* 3)*dstptrs[16*1+ 3]++=((unsigned char*)low)[(16*1+ 3)<<1|1];
					if(cond1&3<<2* 4)*dstptrs[16*1+ 4]++=((unsigned char*)low)[(16*1+ 4)<<1|1];
					if(cond1&3<<2* 5)*dstptrs[16*1+ 5]++=((unsigned char*)low)[(16*1+ 5)<<1|1];
					if(cond1&3<<2* 6)*dstptrs[16*1+ 6]++=((unsigned char*)low)[(16*1+ 6)<<1|1];
					if(cond1&3<<2* 7)*dstptrs[16*1+ 7]++=((unsigned char*)low)[(16*1+ 7)<<1|1];
					if(cond1&3<<2* 8)*dstptrs[16*1+ 8]++=((unsigned char*)low)[(16*1+ 8)<<1|1];
					if(cond1&3<<2* 9)*dstptrs[16*1+ 9]++=((unsigned char*)low)[(16*1+ 9)<<1|1];
					if(cond1&3<<2*10)*dstptrs[16*1+10]++=((unsigned char*)low)[(16*1+10)<<1|1];
					if(cond1&3<<2*11)*dstptrs[16*1+11]++=((unsigned char*)low)[(16*1+11)<<1|1];
					if(cond1&3<<2*12)*dstptrs[16*1+12]++=((unsigned char*)low)[(16*1+12)<<1|1];
					if(cond1&3<<2*13)*dstptrs[16*1+13]++=((unsigned char*)low)[(16*1+13)<<1|1];
					if(cond1&3<<2*14)*dstptrs[16*1+14]++=((unsigned char*)low)[(16*1+14)<<1|1];
					if(cond1&3<<2*15)*dstptrs[16*1+15]++=((unsigned char*)low)[(16*1+15)<<1|1];
					if(cond2&3<<2* 0)*dstptrs[16*2+ 0]++=((unsigned char*)low)[(16*2+ 0)<<1|1];
					if(cond2&3<<2* 1)*dstptrs[16*2+ 1]++=((unsigned char*)low)[(16*2+ 1)<<1|1];
					if(cond2&3<<2* 2)*dstptrs[16*2+ 2]++=((unsigned char*)low)[(16*2+ 2)<<1|1];
					if(cond2&3<<2* 3)*dstptrs[16*2+ 3]++=((unsigned char*)low)[(16*2+ 3)<<1|1];
					if(cond2&3<<2* 4)*dstptrs[16*2+ 4]++=((unsigned char*)low)[(16*2+ 4)<<1|1];
					if(cond2&3<<2* 5)*dstptrs[16*2+ 5]++=((unsigned char*)low)[(16*2+ 5)<<1|1];
					if(cond2&3<<2* 6)*dstptrs[16*2+ 6]++=((unsigned char*)low)[(16*2+ 6)<<1|1];
					if(cond2&3<<2* 7)*dstptrs[16*2+ 7]++=((unsigned char*)low)[(16*2+ 7)<<1|1];
					if(cond2&3<<2* 8)*dstptrs[16*2+ 8]++=((unsigned char*)low)[(16*2+ 8)<<1|1];
					if(cond2&3<<2* 9)*dstptrs[16*2+ 9]++=((unsigned char*)low)[(16*2+ 9)<<1|1];
					if(cond2&3<<2*10)*dstptrs[16*2+10]++=((unsigned char*)low)[(16*2+10)<<1|1];
					if(cond2&3<<2*11)*dstptrs[16*2+11]++=((unsigned char*)low)[(16*2+11)<<1|1];
					if(cond2&3<<2*12)*dstptrs[16*2+12]++=((unsigned char*)low)[(16*2+12)<<1|1];
					if(cond2&3<<2*13)*dstptrs[16*2+13]++=((unsigned char*)low)[(16*2+13)<<1|1];
					if(cond2&3<<2*14)*dstptrs[16*2+14]++=((unsigned char*)low)[(16*2+14)<<1|1];
					if(cond2&3<<2*15)*dstptrs[16*2+15]++=((unsigned char*)low)[(16*2+15)<<1|1];
					__m256i t0=_mm256_srli_epi16(mlow0, 8);
					__m256i t1=_mm256_srli_epi16(mlow1, 8);
					__m256i t2=_mm256_srli_epi16(mlow2, 8);
					__m256i r0=_mm256_srli_epi16(mrange0, 8);
					__m256i r1=_mm256_srli_epi16(mrange1, 8);
					__m256i r2=_mm256_srli_epi16(mrange2, 8);
					r0=_mm256_or_si256(r0, lomask);
					r1=_mm256_or_si256(r1, lomask);
					r2=_mm256_or_si256(r2, lomask);
					mlow0=_mm256_blendv_epi8(mlow0, t0, mcond0);
					mlow1=_mm256_blendv_epi8(mlow1, t1, mcond1);
					mlow2=_mm256_blendv_epi8(mlow2, t2, mcond2);
					mrange0=_mm256_blendv_epi8(mrange0, r0, mcond0);
					mrange1=_mm256_blendv_epi8(mrange1, r1, mcond1);
					mrange2=_mm256_blendv_epi8(mrange2, r2, mcond2);
					__m256i rmax0=_mm256_xor_si256(mlow0, ones);
					__m256i rmax1=_mm256_xor_si256(mlow1, ones);
					__m256i rmax2=_mm256_xor_si256(mlow2, ones);
					mcond0=_mm256_min_epu16(mrange0, rmax0);
					mcond1=_mm256_min_epu16(mrange1, rmax1);
					mcond2=_mm256_min_epu16(mrange2, rmax2);
					mcond0=_mm256_cmpeq_epi16(mcond0, rmax0);
					mcond1=_mm256_cmpeq_epi16(mcond1, rmax1);
					mcond2=_mm256_cmpeq_epi16(mcond2, rmax2);
					mrange0=_mm256_blendv_epi8(mrange0, rmax0, mcond0);
					mrange1=_mm256_blendv_epi8(mrange1, rmax1, mcond1);
					mrange2=_mm256_blendv_epi8(mrange2, rmax2, mcond2);

					__m256i mid0=_mm256_mulhi_epu16(mrange0, mprob0);//get mids
					__m256i mid1=_mm256_mulhi_epu16(mrange1, mprob1);
					__m256i mid2=_mm256_mulhi_epu16(mrange2, mprob2);
					t0=_mm256_add_epi16(mlow0, mid0);
					t1=_mm256_add_epi16(mlow1, mid1);
					t2=_mm256_add_epi16(mlow2, mid2);
					mrange0=_mm256_sub_epi16(mrange0, mid0);
					mrange1=_mm256_sub_epi16(mrange1, mid1);
					mrange2=_mm256_sub_epi16(mrange2, mid2);

					//get bits
					__m256i mbit0=_mm256_and_si256(msym0, bitmask);
					__m256i mbit1=_mm256_and_si256(msym1, bitmask);
					__m256i mbit2=_mm256_and_si256(msym2, bitmask);
					mbit0=_mm256_cmpeq_epi16(mbit0, _mm256_setzero_si256());
					mbit1=_mm256_cmpeq_epi16(mbit1, _mm256_setzero_si256());
					mbit2=_mm256_cmpeq_epi16(mbit2, _mm256_setzero_si256());

					//update ranges
					mlow0=_mm256_blendv_epi8(mlow0, t0, mbit0);
					mlow1=_mm256_blendv_epi8(mlow1, t1, mbit1);
					mlow2=_mm256_blendv_epi8(mlow2, t2, mbit2);
					mrange0=_mm256_blendv_epi8(mid0, mrange0, mbit0);
					mrange1=_mm256_blendv_epi8(mid1, mrange1, mbit1);
					mrange2=_mm256_blendv_epi8(mid2, mrange2, mbit2);

					mbit0=_mm256_and_si256(mbit0, one);
					mbit1=_mm256_and_si256(mbit1, one);
					mbit2=_mm256_and_si256(mbit2, one);
					mbit0=_mm256_xor_si256(mbit0, one);
					mbit1=_mm256_xor_si256(mbit1, one);
					mbit2=_mm256_xor_si256(mbit2, one);
					_mm256_store_si256((__m256i*)bits+0, mbit0);
					_mm256_store_si256((__m256i*)bits+1, mbit1);
					_mm256_store_si256((__m256i*)bits+2, mbit2);
					for(int k=0;k<16;++k)
					{
						stats[0][tidx[16*0+k]]+=((bits[16*0+k]<<16)-stats[0][tidx[16*0+k]]+(1<<5>>1))>>5;
						stats[1][tidx[16*1+k]]+=((bits[16*1+k]<<16)-stats[1][tidx[16*1+k]]+(1<<5>>1))>>5;
						stats[2][tidx[16*2+k]]+=((bits[16*2+k]<<16)-stats[2][tidx[16*2+k]]+(1<<5>>1))>>5;
						tidx[16*0+k]=tidx[16*0+k]+bits[16*0+k];
						tidx[16*1+k]=tidx[16*1+k]+bits[16*1+k];
						tidx[16*2+k]=tidx[16*2+k]+bits[16*2+k];
					}
				}
#if 0
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
				for(int k=0;k<48;++k)
				{
					unsigned short rmax;
#if 0
					int cond=range[k]<256;
					int cmask=-cond;
					int sh=cmask&8;
					dstptrs[k]+=cond;
					*dstptrs[k]=low[k]>>8;
					low[k]<<=sh;//blendv
					range[k]<<=sh;
					range[k]|=(unsigned char)cmask;
					rmax=~low[k];
					if(range[k]>rmax)
						range[k]=rmax;
#endif
#if 1
					if(range[k]<256)
					{
						*dstptrs[k]++=low[k]>>8;
						low[k]<<=8;
						range[k]=range[k]<<8|255;
						rmax=~low[k];
						if(range[k]>rmax)
							range[k]=rmax;
					}
#endif
				}
#endif
				rows[0]+=4;
				rows[1]+=4;
			}
		}
		for(int k=0;k<48;++k)
			free(dstbufs[k]);
	}
	else
	{
	}
	_mm_free(pixels);
#if 0
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

		guide_save(image, iw, ih);

		int hperiod=0x7FFF/iw, hcount=(ih+hperiod-1)/hperiod+2;
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
			-1, -1, -1, -1, -1, -1, -1, 12+2, -1,  9+2, -1,  6+2, -1,  3+2, -1,  0+2
		);
		__m128i getV=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,    8,  7,    6,  5,    4,  3,    2,  1,    0
			-1, -1, -1, -1, -1, -1, -1, 12+0, -1,  9+0, -1,  6+0, -1,  3+0, -1,  0+0
		);
		__m128i mmin=_mm_set1_epi16(0);
		__m128i mmax=_mm_set1_epi16(255);
		for(int ky=ih-1, ky3=0;ky>=1;--ky)		//predict backwards (inplace)
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
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);//signed modular arithmetic
				sym=sym>>(31-8)^sym>>31;//pack sign
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;
				DEB();

				NW	=Nptr	[2-1*3]-Nptr	[1-1*3];
				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				W	=ptr	[2-1*3]-ptr	[1-1*3];
				curr	=ptr	[2+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];
				DEB();

				NW	=Nptr	[0-1*3]-Nptr	[1-1*3];
				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				W	=ptr	[0-1*3]-ptr	[1-1*3];
				curr	=ptr	[0+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];
				DEB();

				Nptr	-=3;
				ptr	-=3;
			}
			kx-=5-1;
			Nptr	-=3*(5-1);
			ptr	-=3*(5-1);
			for(;;)			//fast predict
			{
				//if(ky==1&&kx>=38-5&&kx<38)//
				//	printf("");

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
#ifdef DISABLE_PRED
				mp0=_mm_setzero_si128();
				mp1=_mm_setzero_si128();
				mp2=_mm_setzero_si128();
#endif
				mp0=_mm_sub_epi16(mcurr0, mp0);
				mp1=_mm_sub_epi16(mcurr1, mp1);
				mp2=_mm_sub_epi16(mcurr2, mp2);
				mp0=_mm_slli_epi16(mp0, 16-8);
				mp1=_mm_slli_epi16(mp1, 16-8);
				mp2=_mm_slli_epi16(mp2, 16-8);
				mp0=_mm_xor_si128(_mm_srai_epi16(mp0, 15), _mm_srai_epi16(mp0, 7));//pack sign
				mp1=_mm_xor_si128(_mm_srai_epi16(mp1, 15), _mm_srai_epi16(mp1, 7));
				mp2=_mm_xor_si128(_mm_srai_epi16(mp2, 15), _mm_srai_epi16(mp2, 7));
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
				curr_dens	[2*0+0][0]+=5;
				curr_dens	[2*0+1][0]+=5;
				curr_dens	[2*1+0][0]+=5;
				curr_dens	[2*1+1][0]+=5;
				curr_dens	[2*2+0][0]+=5;
				curr_dens	[2*2+1][0]+=5;
				DEB();

				if(kx<5+1)
					break;
				kx-=5;
				Nptr	-=3*5;
				ptr	-=3*5;
			}
			--kx;
			Nptr	-=3;
			ptr	-=3;
			for(;kx>0;--kx)
			{
				int NW, N, W, curr;
				int offset, pred, sym;

				NW	=Nptr	[1-1*3]-128;
				N	=Nptr	[1+0*3]-128;
				W	=ptr	[1-1*3]-128;
				curr	=ptr	[1+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);//signed modular arithmetic
				sym=sym>>(31-8)^sym>>31;//pack sign
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;
				DEB();

				NW	=Nptr	[2-1*3]-Nptr	[1-1*3];
				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				W	=ptr	[2-1*3]-ptr	[1-1*3];
				curr	=ptr	[2+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];
				DEB();

				NW	=Nptr	[0-1*3]-Nptr	[1-1*3];
				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				W	=ptr	[0-1*3]-ptr	[1-1*3];
				curr	=ptr	[0+0*3]-128;
				CLAMPGRAD(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];
				DEB();

				Nptr	-=3;
				ptr	-=3;
			}
			{
				int N, curr;
				int offset, pred, sym;

				N	=Nptr	[1+0*3]-128;
				curr	=ptr	[1+0*3]-128;
				pred=N;
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;
				DEB();

				N	=Nptr	[2+0*3]-Nptr	[1+0*3];
				curr	=ptr	[2+0*3]-128;
				pred=N;
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];
				DEB();

				N	=Nptr	[0+0*3]-Nptr	[1+0*3];
				curr	=ptr	[0+0*3]-128;
				pred=N;
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];
				DEB();

				Nptr	-=3;
				ptr	-=3;
			}
			if(!(ky%hperiod))
			{
				DEB();
				curr_hists	[2*0+(ky3&1)]+=256*2;//skip 2 histograms (current and next)
				curr_hists	[2*1+(ky3&1)]+=256*2;
				curr_hists	[2*2+(ky3&1)]+=256*2;
				curr_dens	[2*0+(ky3&1)]+=2;
				curr_dens	[2*1+(ky3&1)]+=2;
				curr_dens	[2*2+(ky3&1)]+=2;
				DEB();
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
				int offset, pred, sym;

				W	=kx?ptr	[1-1*3]-128:0;
				curr	=ptr	[1+0*3]-128;
				pred=W;
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[1+0*3]=sym;
				++curr_hists	[2*0+0][sym];
				++curr_hists	[2*0+1][sym];
				++curr_dens	[2*0+0][0];
				++curr_dens	[2*0+1][0];
				offset=curr;
				DEB();

				W	=kx?ptr	[2-1*3]-ptr	[1-1*3]:0;
				curr	=ptr	[2+0*3]-128;
				pred=W;
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[2+0*3]=sym;
				++curr_hists	[2*1+0][sym];
				++curr_hists	[2*1+1][sym];
				++curr_dens	[2*1+0][0];
				++curr_dens	[2*1+1][0];
				DEB();

				W	=kx?ptr	[0-1*3]-ptr	[1-1*3]:0;
				curr	=ptr	[0+0*3]-128;
				pred=W;
				pred+=offset;
				CLAMP2(pred, -128, 127);
#ifdef DISABLE_PRED
				pred=0;
#endif
				sym=(curr-pred)<<(32-8);
				sym=sym>>(31-8)^sym>>31;
				ptr[0+0*3]=sym;
				++curr_hists	[2*2+0][sym];
				++curr_hists	[2*2+1][sym];
				++curr_dens	[2*2+0][0];
				++curr_dens	[2*2+1][0];
				DEB();

			//	Nptr	-=3;
				ptr	-=3;
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

				//if(ky==1&&kx==38)//
				//if(ky==2&&kx==158)//
				//if(ky==2&&kx==157)//
				//	printf("");

				//enc V
#define KC 2
				sym=yuv[KC];
				DEB();
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
				if((state>>16)>=freq)
				//if(state>=((unsigned long long)freq<<32)/den)
				//if((unsigned long long)state*den>=(unsigned long long)freq<<32)
				//if(state>=((0x10000/den)<<16)*freq)
				//if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
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
				DEB();
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
				if((state>>16)>=freq)
				//if(state>=((unsigned long long)freq<<32)/den)
				//if((unsigned long long)state*den>=(unsigned long long)freq<<32)
				//if(state>=((0x10000/den)<<16)*freq)
				//if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
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
				DEB();
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
				if((state>>16)>=freq)
				//if(state>=((unsigned long long)freq<<32)/den)
				//if((unsigned long long)state*den>=(unsigned long long)freq<<32)
				//if(state>=((0x10000/den)<<16)*freq)
				//if((unsigned)((unsigned long long)state*(unsigned)den>>32)>=(unsigned)freq)
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
				DEB();
				curr_hists	[2*0+(ky3&1)]+=256*2;//skip 2 histograms (current and next)
				curr_hists	[2*1+(ky3&1)]+=256*2;
				curr_hists	[2*2+(ky3&1)]+=256*2;
				curr_dens	[2*0+(ky3&1)]+=2;
				curr_dens	[2*1+(ky3&1)]+=2;
				curr_dens	[2*2+(ky3&1)]+=2;
				DEB();
				++ky3;
			}
		}
		*dstptr++=(unsigned short)state;
		*dstptr++=(unsigned short)(state>>16);

		{//save dstptr-dstbuf
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
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


		int hperiod=0x7FFF/iw;
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
			if(ky==ky2)
			{
				DEB2();
				memset(hist+256*(2LL*0+(ky3&1)), 0, sizeof(short[256]));
				memset(hist+256*(2LL*1+(ky3&1)), 0, sizeof(short[256]));
				memset(hist+256*(2LL*2+(ky3&1)), 0, sizeof(short[256]));
				dens[2LL*0+(ky3&1)]=0;
				dens[2LL*1+(ky3&1)]=0;
				dens[2LL*2+(ky3&1)]=0;
				ky2+=hperiod;
				++ky3;
			}
			for(int kx=0;kx<iw;++kx, idx+=3)//decode
			{
				//if(ky==1&&kx==37)//
				//if(ky==2&&kx==158)//
				if(ky==2&&kx==157)//
					printf("");
				
				short *curr=rows[0];
				__m128i mNW	=_mm_load_si128((__m128i*)(rows[1]-1*4*2));
				__m128i mN	=_mm_load_si128((__m128i*)(rows[1]+0*4*2));
				__m128i mW	=_mm_load_si128((__m128i*)(rows[0]-1*4*2));
				if(!kx)
					mNW=mW=mN;
				if(!ky)
					mNW=mN=mW;
				__m128i vmin=_mm_min_epi16(mN, mW);
				__m128i vmax=_mm_max_epi16(mN, mW);
				__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
				mp=_mm_max_epi16(mp, vmin);
				mp=_mm_min_epi16(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);
#ifdef DISABLE_PRED
				_mm_store_si128((__m128i*)preds, _mm_setzero_si128());
#endif

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
				if(state<0x10000)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				DEB2();
				sym=sym>>1^-(sym&1);
				curr[KC]=(sym+preds[KC])<<24>>24;
#undef  KC

				preds[1]+=curr[0];
				CLAMP2(preds[1], -128, 127);
				preds[2]+=curr[0];
				CLAMP2(preds[2], -128, 127);
#ifdef DISABLE_PRED
				_mm_store_si128((__m128i*)preds, _mm_setzero_si128());
#endif

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
				if(state<0x10000)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				DEB2();
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
				if(state<0x10000)
					state=state<<16|*ptr--;
				++hists	[2*KC+0][sym];
				++hists	[2*KC+1][sym];
				++dens	[2*KC+0];
				++dens	[2*KC+1];
				DEB2();
				sym=sym>>1^-(sym&1);
				curr[KC]=(sym+preds[KC])<<24>>24;
#undef  KC

			//	curr[4+0]=curr[0]-preds[0];
			//	curr[4+1]=curr[1]-preds[1];
			//	curr[4+2]=curr[2]-preds[2];
				image[idx+1]=(unsigned char)(curr[0]+128);
				image[idx+2]=(unsigned char)(curr[1]+128);
				image[idx+0]=(unsigned char)(curr[2]+128);
				guide_check(image, kx, ky);
				curr[1]-=curr[0];
				curr[2]-=curr[0];

				rows[0]+=4*2;
				rows[1]+=4*2;
			}
			//if(ptr<loptr)
			//	LOG_ERROR("");
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
#endif
	free(srcbuf);
#ifdef ESTIMATE_SIZE
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
	if(fwd)
	{
		for(int kc=0;kc<3;++kc)
			printf("%c%12.2lf\n", "YUV"[kc], csizes[kc]/8);
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
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