#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#include"blist.h"
static const char file[]=__FILE__;


//	#define ESTIMATE_SIZE
//	#define ENABLE_GUIDE
//	#define DISABLE_PRED
//	#define ENABLE_VALIDATION


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

#ifdef ENABLE_VALIDATION
static ArrayHandle g_probs=0;
static ptrdiff_t g_probidx=0;
static void val_append(unsigned short *probs, int nprobs)
{
	if(!g_probs)
		ARRAY_ALLOC(short, g_probs, 0, 0, 0, 0);
	ARRAY_APPEND(g_probs, probs, nprobs, 1, 0);
}
static void val_check(unsigned short *probs, int nprobs)
{
	unsigned short *correctprobs=(unsigned short*)array_at(&g_probs, g_probidx);
	if(memcmp(probs, correctprobs, nprobs*sizeof(short)))
	{
		for(int k=0;k<nprobs;++k)
			printf("%8td %04X %04X\n", g_probidx+k, probs[k], correctprobs[k]);
		LOG_ERROR("Val Error");
	}
	g_probidx+=nprobs;
}
#else
#define val_append(...)
#define val_check(...)
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
#ifdef _MSC_VER
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
#endif
#ifdef ESTIMATE_SIZE
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
	}
	else//decode
	{
		if(srcptr+4*3>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
	}
	
	int psize=(iw+32LL)*sizeof(short[4*3]);//4 padded rows * 3 channels
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
	if(fwd)//encode
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
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*0)+16LL),//Y
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*1)+16LL),//U
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*2)+16LL),//V
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*0)+16LL),
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*1)+16LL),
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*2)+16LL),
			};
			ALIGN(32) unsigned short bits[48]={0};
			ALIGN(32) unsigned char tidx[48]={0};
			__m256i himask=_mm256_set1_epi16(0xFF00);
			__m256i lomask=_mm256_set1_epi16(0x00FF);
			__m256i probmin=_mm256_set1_epi16(0x0102);//0xFF*0x0102>>16 = 1
			__m256i probmax=_mm256_set1_epi16(0xFEFE);//0xFF*0xFEFE>>16 = 253
			__m256i ones=_mm256_set1_epi16(-1);
			__m256i one=_mm256_set1_epi16(1);
			__m256i vmin=_mm256_set1_epi16(-128);
			__m256i vmax=_mm256_set1_epi16(+127);

			int idx=3*iw*ky;
			for(int kx=0;kx<iw;kx+=16, idx+=48)
			{
				short
					*NW0	=rows[2*0+1]-1,
					*NW1	=rows[2*1+1]-1,
					*NW2	=rows[2*2+1]-1,
					*N0	=rows[2*0+1]+0,
					*N1	=rows[2*1+1]+0,
					*N2	=rows[2*2+1]+0,
					*W0	=rows[2*0+0]-1,
					*W1	=rows[2*1+0]-1,
					*W2	=rows[2*2+0]-1,
					*curr0	=rows[2*0+0]+0,
					*curr1	=rows[2*1+0]+0,
					*curr2	=rows[2*2+0]+0;

				int remainder=iw-kx;
				if(remainder>16)
					remainder=16;
				remainder*=3;
				for(int kd=0, k=0;k<remainder;++kd, k+=3)
				{
					curr0[kd]=image[idx+k+1]-128;//Y
					curr1[kd]=image[idx+k+2]-128;//U
					curr2[kd]=image[idx+k+0]-128;//V
				}
				__m256i mNW0=	_mm256_loadu_si256((__m256i*)NW0);
				__m256i mNW1=	_mm256_loadu_si256((__m256i*)NW1);
				__m256i mNW2=	_mm256_loadu_si256((__m256i*)NW2);
				__m256i mN0=	_mm256_loadu_si256((__m256i*)N0);
				__m256i mN1=	_mm256_loadu_si256((__m256i*)N1);
				__m256i mN2=	_mm256_loadu_si256((__m256i*)N2);
				__m256i mW0=	_mm256_loadu_si256((__m256i*)W0);
				__m256i mW1=	_mm256_loadu_si256((__m256i*)W1);
				__m256i mW2=	_mm256_loadu_si256((__m256i*)W2);
				__m256i mcurr0=	_mm256_loadu_si256((__m256i*)curr0);
				__m256i mcurr1=	_mm256_loadu_si256((__m256i*)curr1);
				__m256i mcurr2=	_mm256_loadu_si256((__m256i*)curr2);
				
				mNW1=_mm256_sub_epi16(mNW1, mNW0);
				mNW2=_mm256_sub_epi16(mNW2, mNW0);
				mN1=_mm256_sub_epi16(mN1, mN0);
				mN2=_mm256_sub_epi16(mN2, mN0);
				mW1=_mm256_sub_epi16(mW1, mW0);
				mW2=_mm256_sub_epi16(mW2, mW0);
				__m256i vmin0=_mm256_min_epi16(mN0, mW0);
				__m256i vmin1=_mm256_min_epi16(mN1, mW1);
				__m256i vmin2=_mm256_min_epi16(mN2, mW2);
				__m256i vmax0=_mm256_max_epi16(mN0, mW0);
				__m256i vmax1=_mm256_max_epi16(mN1, mW1);
				__m256i vmax2=_mm256_max_epi16(mN2, mW2);
				__m256i mp0=_mm256_add_epi16(mN0, mW0);
				__m256i mp1=_mm256_add_epi16(mN1, mW1);
				__m256i mp2=_mm256_add_epi16(mN2, mW2);
				mp0=_mm256_sub_epi16(mp0, mNW0);
				mp1=_mm256_sub_epi16(mp1, mNW1);
				mp2=_mm256_sub_epi16(mp2, mNW2);
				mp0=_mm256_max_epi16(mp0, vmin0);
				mp1=_mm256_max_epi16(mp1, vmin1);
				mp2=_mm256_max_epi16(mp2, vmin2);
				mp0=_mm256_min_epi16(mp0, vmax0);
				mp1=_mm256_min_epi16(mp1, vmax1);
				mp2=_mm256_min_epi16(mp2, vmax2);

				mp1=_mm256_add_epi16(mp1, mcurr0);
				mp2=_mm256_add_epi16(mp2, mcurr0);
				mp1=_mm256_max_epi16(mp1, vmin);
				mp2=_mm256_max_epi16(mp2, vmin);
				mp1=_mm256_min_epi16(mp1, vmax);
				mp2=_mm256_min_epi16(mp2, vmax);

#ifdef DISABLE_PRED
				mp0=_mm256_setzero_si256();
				mp1=_mm256_setzero_si256();
				mp2=_mm256_setzero_si256();
#endif
				//if(ky==2&&kx==16)//
				//	printf("");

				__m256i bitmask=_mm256_set1_epi16(0x80);
				__m256i msym0=_mm256_sub_epi16(mcurr0, mp0);
				__m256i msym1=_mm256_sub_epi16(mcurr1, mp1);
				__m256i msym2=_mm256_sub_epi16(mcurr2, mp2);
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

					_mm256_store_si256((__m256i*)low+0, mlow0);
					_mm256_store_si256((__m256i*)low+1, mlow1);
					_mm256_store_si256((__m256i*)low+2, mlow2);
					__m256i mcond0=_mm256_and_si256(mrange0, himask);
					__m256i mcond1=_mm256_and_si256(mrange1, himask);
					__m256i mcond2=_mm256_and_si256(mrange2, himask);
					mcond0=_mm256_cmpeq_epi16(mcond0, _mm256_setzero_si256());
					mcond1=_mm256_cmpeq_epi16(mcond1, _mm256_setzero_si256());
					mcond2=_mm256_cmpeq_epi16(mcond2, _mm256_setzero_si256());
					int cond0=_mm256_movemask_epi8(mcond0);
					int cond1=_mm256_movemask_epi8(mcond1);
					int cond2=_mm256_movemask_epi8(mcond2);

					//if(cond0)//
					//	printf("0");
					//if(cond1)//
					//	printf("1");
					//if(cond2)//
					//	printf("2");
					//{
					//	_mm256_store_si256((__m256i*)range+0, mrange0);
					//	_mm256_store_si256((__m256i*)range+1, mrange1);
					//	_mm256_store_si256((__m256i*)range+2, mrange2);
					//	int cond[]={cond0, cond1, cond2};
					//	for(int k=0;k<48;++k)//
					//	{
					//		if(cond[k>>4]&3<<2*(k&15)&&((unsigned char*)range)[k<<1|1])
					//			LOG_ERROR("");
					//	}
					//}

					if(cond0)
					{
						if(cond0&3<<2* 0)*dstptrs[16*0+ 0]++=((unsigned char*)low)[(16*0+ 0)<<1|1];//should be {TEST, CMOV, INC}, but it's not
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
					}
					if(cond1)
					{
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
					}
					if(cond2)
					{
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
					}
					__m256i t0=_mm256_slli_epi16(mlow0, 8);
					__m256i t1=_mm256_slli_epi16(mlow1, 8);
					__m256i t2=_mm256_slli_epi16(mlow2, 8);
					__m256i r0=_mm256_slli_epi16(mrange0, 8);
					__m256i r1=_mm256_slli_epi16(mrange1, 8);
					__m256i r2=_mm256_slli_epi16(mrange2, 8);
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
					
					//if(ky==0&&kx==16&&kb==6)//
					//if(ky==0&&kx==752&&kb==0)//
					//	printf("");

					mrange0=_mm256_sub_epi16(mrange0, mid0);
					mrange1=_mm256_sub_epi16(mrange1, mid1);
					mrange2=_mm256_sub_epi16(mrange2, mid2);

					//get bits
					__m256i mbit0=_mm256_and_si256(msym0, bitmask);
					__m256i mbit1=_mm256_and_si256(msym1, bitmask);
					__m256i mbit2=_mm256_and_si256(msym2, bitmask);
					mbit0=_mm256_cmpeq_epi16(mbit0, bitmask);
					mbit1=_mm256_cmpeq_epi16(mbit1, bitmask);
					mbit2=_mm256_cmpeq_epi16(mbit2, bitmask);
					bitmask=_mm256_srli_epi16(bitmask, 1);
#ifdef ESTIMATE_SIZE
					{
						ALIGN(32) unsigned short probs[48];
						_mm256_store_si256((__m256i*)probs+0, mprob0);
						_mm256_store_si256((__m256i*)probs+1, mprob1);
						_mm256_store_si256((__m256i*)probs+2, mprob2);
						_mm256_store_si256((__m256i*)bits+0, mbit0);
						_mm256_store_si256((__m256i*)bits+1, mbit1);
						_mm256_store_si256((__m256i*)bits+2, mbit2);
						for(int k=0;k<16;++k)
						{
							csizes[0]-=log2((double)(bits[16*0+k]?0x10000-probs[16*0+k]:probs[16*0+k])/0x10000);
							csizes[1]-=log2((double)(bits[16*1+k]?0x10000-probs[16*1+k]:probs[16*1+k])/0x10000);
							csizes[2]-=log2((double)(bits[16*2+k]?0x10000-probs[16*2+k]:probs[16*2+k])/0x10000);
						}
					}
#endif
#ifdef ENABLE_VALIDATION
					{
						ALIGN(32) unsigned short probs[48];
						_mm256_store_si256((__m256i*)probs+0, mprob0);
						_mm256_store_si256((__m256i*)probs+1, mprob1);
						_mm256_store_si256((__m256i*)probs+2, mprob2);
						_mm256_store_si256((__m256i*)bits+0, mbit0);
						_mm256_store_si256((__m256i*)bits+1, mbit1);
						_mm256_store_si256((__m256i*)bits+2, mbit2);
						for(int k=0;k<48;++k)
						{
							probs[k]^=bits[k];
							probs[k]-=bits[k];
						}
						val_append(probs, 48);
					}
#endif
					//update ranges
					mlow0=_mm256_blendv_epi8(mlow0, t0, mbit0);
					mlow1=_mm256_blendv_epi8(mlow1, t1, mbit1);
					mlow2=_mm256_blendv_epi8(mlow2, t2, mbit2);
					mid0=_mm256_sub_epi16(mid0, one);
					mid1=_mm256_sub_epi16(mid1, one);
					mid2=_mm256_sub_epi16(mid2, one);
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
#if 1
#ifdef __GNUC__
#pragma GCC unroll 16
#endif
					for(int k=0;k<16;++k)
					{
						unsigned short
							*s0=stats[0]+tidx[16*0+k],
							*s1=stats[1]+tidx[16*1+k],
							*s2=stats[2]+tidx[16*2+k];
						*s0+=((bits[16*0+k]<<16)-*s0+(1<<5>>1))>>5;
						*s1+=((bits[16*1+k]<<16)-*s1+(1<<5>>1))>>5;
						*s2+=((bits[16*2+k]<<16)-*s2+(1<<5>>1))>>5;
						//stats[0][tidx[16*0+k]]+=((bits[16*0+k]<<16)-stats[0][tidx[16*0+k]]+(1<<5>>1))>>5;
						//stats[1][tidx[16*1+k]]+=((bits[16*1+k]<<16)-stats[1][tidx[16*1+k]]+(1<<5>>1))>>5;
						//stats[2][tidx[16*2+k]]+=((bits[16*2+k]<<16)-stats[2][tidx[16*2+k]]+(1<<5>>1))>>5;
						tidx[16*0+k]=2*tidx[16*0+k]+bits[16*0+k];
						tidx[16*1+k]=2*tidx[16*1+k]+bits[16*1+k];
						tidx[16*2+k]=2*tidx[16*2+k]+bits[16*2+k];
					}
#endif
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
				rows[0]+=16;
				rows[1]+=16;
				rows[2]+=16;
				rows[3]+=16;
				rows[4]+=16;
				rows[5]+=16;
			}
		}
		_mm256_store_si256((__m256i*)low+0, mlow0);
		_mm256_store_si256((__m256i*)low+1, mlow1);
		_mm256_store_si256((__m256i*)low+2, mlow2);
		for(int k=0;k<48;++k)//flush
		{
			*dstptrs[k]++=low[k]>>8;
			*dstptrs[k]++=low[k]&255;
		}
		{
			ptrdiff_t csize=0;
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			csize+=fwrite("CH", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 4, fdst);
			csize+=fwrite(&ih, 1, 4, fdst);
			for(int k=0;k<48;++k)
			{
				ptrdiff_t size=dstptrs[k]-dstbufs[k];
				csize+=fwrite(&size, 1, 4, fdst);
#ifdef ESTIMATE_SIZE
				printf("%2d: %10td\n", k, size);
#endif
			}
			for(int k=0;k<48;++k)
			{
				ptrdiff_t size=dstptrs[k]-dstbufs[k];
				csize+=fwrite(dstbufs[k], 1, size, fdst);
			}
#ifdef ESTIMATE_SIZE
			csize_actual=csize;
#endif
			fclose(fdst);
		}
		for(int k=0;k<48;++k)
			free(dstbufs[k]);
	}
	else//decode
	{
		int imsize=3*iw*ih;
		unsigned char *image=(unsigned char*)malloc(imsize);
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}

		int streamsizes[48]={0};
		unsigned char *srcptrs[48]={0};
		memcpy(streamsizes, srcptr, 4LL*48); srcptr+=4LL*48;
		int csize=0;
		for(int k=0;k<48;++k)
		{
			int size=streamsizes[k];
			srcptrs[k]=srcptr;
			srcptr+=size;
			csize+=size;
		}
		if(csize+2LL+4*(2LL+48)!=srcptr-srcbuf)//"CH", iw, ih, 48 stream sizes
		{
			LOG_ERROR("Error: expected csize %d != actual %td", csize+4*2+4*48, srcptr-srcbuf);
			return 1;
		}
		ALIGN(32) unsigned short code[48]={0};
		for(int k=0;k<48;++k)
		{
			code[k]=srcptrs[k][0]<<8|srcptrs[k][1];
			srcptrs[k]+=2;
		}
		//__m256i mcode0=_mm256_load_si256((__m256i*)code+0);
		//__m256i mcode1=_mm256_load_si256((__m256i*)code+1);
		//__m256i mcode2=_mm256_load_si256((__m256i*)code+2);
		for(int ky=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*0)+16LL),//Y
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*1)+16LL),//U
				pixels+((iw+32LL)*(((ky-0LL)&1)+2*2)+16LL),//V
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*0)+16LL),
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*1)+16LL),
				pixels+((iw+32LL)*(((ky-1LL)&1)+2*2)+16LL),
			};
			ALIGN(32) unsigned short syms[48]={0};
			ALIGN(32) unsigned short bits[48]={0};
			ALIGN(32) unsigned char tidx[48]={0};
			__m256i himask=_mm256_set1_epi16(0xFF00);
			__m256i lomask=_mm256_set1_epi16(0x00FF);
			__m256i probmin=_mm256_set1_epi16(0x0102);
			__m256i probmax=_mm256_set1_epi16(0xFEFE);
			__m256i ones=_mm256_set1_epi16(-1);
			__m256i one=_mm256_set1_epi16(1);

			int idx=3*iw*ky;
			for(int kx=0;kx<iw;kx+=16)
			{
				int remainder=iw-kx;
				if(remainder>16)
					remainder=16;
				
				__m256i bitmask=_mm256_set1_epi16(0x80);
				__m256i msym0=_mm256_setzero_si256();
				__m256i msym1=_mm256_setzero_si256();
				__m256i msym2=_mm256_setzero_si256();
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
					
					//_mm256_store_si256((__m256i*)code+0, mcode0);
					//_mm256_store_si256((__m256i*)code+1, mcode1);
					//_mm256_store_si256((__m256i*)code+2, mcode2);
					__m256i mcond0=_mm256_and_si256(mrange0, himask);
					__m256i mcond1=_mm256_and_si256(mrange1, himask);
					__m256i mcond2=_mm256_and_si256(mrange2, himask);
					mcond0=_mm256_cmpeq_epi16(mcond0, _mm256_setzero_si256());
					mcond1=_mm256_cmpeq_epi16(mcond1, _mm256_setzero_si256());
					mcond2=_mm256_cmpeq_epi16(mcond2, _mm256_setzero_si256());
					int cond0=_mm256_movemask_epi8(mcond0);
					int cond1=_mm256_movemask_epi8(mcond1);
					int cond2=_mm256_movemask_epi8(mcond2);
					if(cond0)
					{
						if(cond0&3<<2* 0)((unsigned char*)code)[(16*0+ 0)<<1|1]=((unsigned char*)code)[(16*0+ 0)<<1|0], ((unsigned char*)code)[(16*0+ 0)<<1|0]=*srcptrs[16*0+ 0]++;
						if(cond0&3<<2* 1)((unsigned char*)code)[(16*0+ 1)<<1|1]=((unsigned char*)code)[(16*0+ 1)<<1|0], ((unsigned char*)code)[(16*0+ 1)<<1|0]=*srcptrs[16*0+ 1]++;
						if(cond0&3<<2* 2)((unsigned char*)code)[(16*0+ 2)<<1|1]=((unsigned char*)code)[(16*0+ 2)<<1|0], ((unsigned char*)code)[(16*0+ 2)<<1|0]=*srcptrs[16*0+ 2]++;
						if(cond0&3<<2* 3)((unsigned char*)code)[(16*0+ 3)<<1|1]=((unsigned char*)code)[(16*0+ 3)<<1|0], ((unsigned char*)code)[(16*0+ 3)<<1|0]=*srcptrs[16*0+ 3]++;
						if(cond0&3<<2* 4)((unsigned char*)code)[(16*0+ 4)<<1|1]=((unsigned char*)code)[(16*0+ 4)<<1|0], ((unsigned char*)code)[(16*0+ 4)<<1|0]=*srcptrs[16*0+ 4]++;
						if(cond0&3<<2* 5)((unsigned char*)code)[(16*0+ 5)<<1|1]=((unsigned char*)code)[(16*0+ 5)<<1|0], ((unsigned char*)code)[(16*0+ 5)<<1|0]=*srcptrs[16*0+ 5]++;
						if(cond0&3<<2* 6)((unsigned char*)code)[(16*0+ 6)<<1|1]=((unsigned char*)code)[(16*0+ 6)<<1|0], ((unsigned char*)code)[(16*0+ 6)<<1|0]=*srcptrs[16*0+ 6]++;
						if(cond0&3<<2* 7)((unsigned char*)code)[(16*0+ 7)<<1|1]=((unsigned char*)code)[(16*0+ 7)<<1|0], ((unsigned char*)code)[(16*0+ 7)<<1|0]=*srcptrs[16*0+ 7]++;
						if(cond0&3<<2* 8)((unsigned char*)code)[(16*0+ 8)<<1|1]=((unsigned char*)code)[(16*0+ 8)<<1|0], ((unsigned char*)code)[(16*0+ 8)<<1|0]=*srcptrs[16*0+ 8]++;
						if(cond0&3<<2* 9)((unsigned char*)code)[(16*0+ 9)<<1|1]=((unsigned char*)code)[(16*0+ 9)<<1|0], ((unsigned char*)code)[(16*0+ 9)<<1|0]=*srcptrs[16*0+ 9]++;
						if(cond0&3<<2*10)((unsigned char*)code)[(16*0+10)<<1|1]=((unsigned char*)code)[(16*0+10)<<1|0], ((unsigned char*)code)[(16*0+10)<<1|0]=*srcptrs[16*0+10]++;
						if(cond0&3<<2*11)((unsigned char*)code)[(16*0+11)<<1|1]=((unsigned char*)code)[(16*0+11)<<1|0], ((unsigned char*)code)[(16*0+11)<<1|0]=*srcptrs[16*0+11]++;
						if(cond0&3<<2*12)((unsigned char*)code)[(16*0+12)<<1|1]=((unsigned char*)code)[(16*0+12)<<1|0], ((unsigned char*)code)[(16*0+12)<<1|0]=*srcptrs[16*0+12]++;
						if(cond0&3<<2*13)((unsigned char*)code)[(16*0+13)<<1|1]=((unsigned char*)code)[(16*0+13)<<1|0], ((unsigned char*)code)[(16*0+13)<<1|0]=*srcptrs[16*0+13]++;
						if(cond0&3<<2*14)((unsigned char*)code)[(16*0+14)<<1|1]=((unsigned char*)code)[(16*0+14)<<1|0], ((unsigned char*)code)[(16*0+14)<<1|0]=*srcptrs[16*0+14]++;
						if(cond0&3<<2*15)((unsigned char*)code)[(16*0+15)<<1|1]=((unsigned char*)code)[(16*0+15)<<1|0], ((unsigned char*)code)[(16*0+15)<<1|0]=*srcptrs[16*0+15]++;
					}
					if(cond1)
					{
						if(cond1&3<<2* 0)((unsigned char*)code)[(16*1+ 0)<<1|1]=((unsigned char*)code)[(16*1+ 0)<<1|0], ((unsigned char*)code)[(16*1+ 0)<<1|0]=*srcptrs[16*1+ 0]++;
						if(cond1&3<<2* 1)((unsigned char*)code)[(16*1+ 1)<<1|1]=((unsigned char*)code)[(16*1+ 1)<<1|0], ((unsigned char*)code)[(16*1+ 1)<<1|0]=*srcptrs[16*1+ 1]++;
						if(cond1&3<<2* 2)((unsigned char*)code)[(16*1+ 2)<<1|1]=((unsigned char*)code)[(16*1+ 2)<<1|0], ((unsigned char*)code)[(16*1+ 2)<<1|0]=*srcptrs[16*1+ 2]++;
						if(cond1&3<<2* 3)((unsigned char*)code)[(16*1+ 3)<<1|1]=((unsigned char*)code)[(16*1+ 3)<<1|0], ((unsigned char*)code)[(16*1+ 3)<<1|0]=*srcptrs[16*1+ 3]++;
						if(cond1&3<<2* 4)((unsigned char*)code)[(16*1+ 4)<<1|1]=((unsigned char*)code)[(16*1+ 4)<<1|0], ((unsigned char*)code)[(16*1+ 4)<<1|0]=*srcptrs[16*1+ 4]++;
						if(cond1&3<<2* 5)((unsigned char*)code)[(16*1+ 5)<<1|1]=((unsigned char*)code)[(16*1+ 5)<<1|0], ((unsigned char*)code)[(16*1+ 5)<<1|0]=*srcptrs[16*1+ 5]++;
						if(cond1&3<<2* 6)((unsigned char*)code)[(16*1+ 6)<<1|1]=((unsigned char*)code)[(16*1+ 6)<<1|0], ((unsigned char*)code)[(16*1+ 6)<<1|0]=*srcptrs[16*1+ 6]++;
						if(cond1&3<<2* 7)((unsigned char*)code)[(16*1+ 7)<<1|1]=((unsigned char*)code)[(16*1+ 7)<<1|0], ((unsigned char*)code)[(16*1+ 7)<<1|0]=*srcptrs[16*1+ 7]++;
						if(cond1&3<<2* 8)((unsigned char*)code)[(16*1+ 8)<<1|1]=((unsigned char*)code)[(16*1+ 8)<<1|0], ((unsigned char*)code)[(16*1+ 8)<<1|0]=*srcptrs[16*1+ 8]++;
						if(cond1&3<<2* 9)((unsigned char*)code)[(16*1+ 9)<<1|1]=((unsigned char*)code)[(16*1+ 9)<<1|0], ((unsigned char*)code)[(16*1+ 9)<<1|0]=*srcptrs[16*1+ 9]++;
						if(cond1&3<<2*10)((unsigned char*)code)[(16*1+10)<<1|1]=((unsigned char*)code)[(16*1+10)<<1|0], ((unsigned char*)code)[(16*1+10)<<1|0]=*srcptrs[16*1+10]++;
						if(cond1&3<<2*11)((unsigned char*)code)[(16*1+11)<<1|1]=((unsigned char*)code)[(16*1+11)<<1|0], ((unsigned char*)code)[(16*1+11)<<1|0]=*srcptrs[16*1+11]++;
						if(cond1&3<<2*12)((unsigned char*)code)[(16*1+12)<<1|1]=((unsigned char*)code)[(16*1+12)<<1|0], ((unsigned char*)code)[(16*1+12)<<1|0]=*srcptrs[16*1+12]++;
						if(cond1&3<<2*13)((unsigned char*)code)[(16*1+13)<<1|1]=((unsigned char*)code)[(16*1+13)<<1|0], ((unsigned char*)code)[(16*1+13)<<1|0]=*srcptrs[16*1+13]++;
						if(cond1&3<<2*14)((unsigned char*)code)[(16*1+14)<<1|1]=((unsigned char*)code)[(16*1+14)<<1|0], ((unsigned char*)code)[(16*1+14)<<1|0]=*srcptrs[16*1+14]++;
						if(cond1&3<<2*15)((unsigned char*)code)[(16*1+15)<<1|1]=((unsigned char*)code)[(16*1+15)<<1|0], ((unsigned char*)code)[(16*1+15)<<1|0]=*srcptrs[16*1+15]++;
					}
					if(cond2)
					{
						if(cond2&3<<2* 0)((unsigned char*)code)[(16*2+ 0)<<1|1]=((unsigned char*)code)[(16*2+ 0)<<1|0], ((unsigned char*)code)[(16*2+ 0)<<1|0]=*srcptrs[16*2+ 0]++;
						if(cond2&3<<2* 1)((unsigned char*)code)[(16*2+ 1)<<1|1]=((unsigned char*)code)[(16*2+ 1)<<1|0], ((unsigned char*)code)[(16*2+ 1)<<1|0]=*srcptrs[16*2+ 1]++;
						if(cond2&3<<2* 2)((unsigned char*)code)[(16*2+ 2)<<1|1]=((unsigned char*)code)[(16*2+ 2)<<1|0], ((unsigned char*)code)[(16*2+ 2)<<1|0]=*srcptrs[16*2+ 2]++;
						if(cond2&3<<2* 3)((unsigned char*)code)[(16*2+ 3)<<1|1]=((unsigned char*)code)[(16*2+ 3)<<1|0], ((unsigned char*)code)[(16*2+ 3)<<1|0]=*srcptrs[16*2+ 3]++;
						if(cond2&3<<2* 4)((unsigned char*)code)[(16*2+ 4)<<1|1]=((unsigned char*)code)[(16*2+ 4)<<1|0], ((unsigned char*)code)[(16*2+ 4)<<1|0]=*srcptrs[16*2+ 4]++;
						if(cond2&3<<2* 5)((unsigned char*)code)[(16*2+ 5)<<1|1]=((unsigned char*)code)[(16*2+ 5)<<1|0], ((unsigned char*)code)[(16*2+ 5)<<1|0]=*srcptrs[16*2+ 5]++;
						if(cond2&3<<2* 6)((unsigned char*)code)[(16*2+ 6)<<1|1]=((unsigned char*)code)[(16*2+ 6)<<1|0], ((unsigned char*)code)[(16*2+ 6)<<1|0]=*srcptrs[16*2+ 6]++;
						if(cond2&3<<2* 7)((unsigned char*)code)[(16*2+ 7)<<1|1]=((unsigned char*)code)[(16*2+ 7)<<1|0], ((unsigned char*)code)[(16*2+ 7)<<1|0]=*srcptrs[16*2+ 7]++;
						if(cond2&3<<2* 8)((unsigned char*)code)[(16*2+ 8)<<1|1]=((unsigned char*)code)[(16*2+ 8)<<1|0], ((unsigned char*)code)[(16*2+ 8)<<1|0]=*srcptrs[16*2+ 8]++;
						if(cond2&3<<2* 9)((unsigned char*)code)[(16*2+ 9)<<1|1]=((unsigned char*)code)[(16*2+ 9)<<1|0], ((unsigned char*)code)[(16*2+ 9)<<1|0]=*srcptrs[16*2+ 9]++;
						if(cond2&3<<2*10)((unsigned char*)code)[(16*2+10)<<1|1]=((unsigned char*)code)[(16*2+10)<<1|0], ((unsigned char*)code)[(16*2+10)<<1|0]=*srcptrs[16*2+10]++;
						if(cond2&3<<2*11)((unsigned char*)code)[(16*2+11)<<1|1]=((unsigned char*)code)[(16*2+11)<<1|0], ((unsigned char*)code)[(16*2+11)<<1|0]=*srcptrs[16*2+11]++;
						if(cond2&3<<2*12)((unsigned char*)code)[(16*2+12)<<1|1]=((unsigned char*)code)[(16*2+12)<<1|0], ((unsigned char*)code)[(16*2+12)<<1|0]=*srcptrs[16*2+12]++;
						if(cond2&3<<2*13)((unsigned char*)code)[(16*2+13)<<1|1]=((unsigned char*)code)[(16*2+13)<<1|0], ((unsigned char*)code)[(16*2+13)<<1|0]=*srcptrs[16*2+13]++;
						if(cond2&3<<2*14)((unsigned char*)code)[(16*2+14)<<1|1]=((unsigned char*)code)[(16*2+14)<<1|0], ((unsigned char*)code)[(16*2+14)<<1|0]=*srcptrs[16*2+14]++;
						if(cond2&3<<2*15)((unsigned char*)code)[(16*2+15)<<1|1]=((unsigned char*)code)[(16*2+15)<<1|0], ((unsigned char*)code)[(16*2+15)<<1|0]=*srcptrs[16*2+15]++;
					}
					__m256i mcode0=_mm256_load_si256((__m256i*)code+0);
					__m256i mcode1=_mm256_load_si256((__m256i*)code+1);
					__m256i mcode2=_mm256_load_si256((__m256i*)code+2);
					__m256i t0=_mm256_slli_epi16(mlow0, 8);
					__m256i t1=_mm256_slli_epi16(mlow1, 8);
					__m256i t2=_mm256_slli_epi16(mlow2, 8);
					__m256i r0=_mm256_slli_epi16(mrange0, 8);
					__m256i r1=_mm256_slli_epi16(mrange1, 8);
					__m256i r2=_mm256_slli_epi16(mrange2, 8);
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

					//if(ky==0&&kx==16&&kb==6)//
					//if(ky==0&&kx==752&&kb==0)//
					//	printf("");

					//decode bits	bit = (code >= low+mid)
					__m256i mbit0=_mm256_min_epu16(mcode0, t0);
					__m256i mbit1=_mm256_min_epu16(mcode1, t1);
					__m256i mbit2=_mm256_min_epu16(mcode2, t2);
					mbit0=_mm256_cmpeq_epi16(mbit0, t0);
					mbit1=_mm256_cmpeq_epi16(mbit1, t1);
					mbit2=_mm256_cmpeq_epi16(mbit2, t2);

					//update ranges
					mlow0=_mm256_blendv_epi8(mlow0, t0, mbit0);
					mlow1=_mm256_blendv_epi8(mlow1, t1, mbit1);
					mlow2=_mm256_blendv_epi8(mlow2, t2, mbit2);
					mid0=_mm256_sub_epi16(mid0, one);
					mid1=_mm256_sub_epi16(mid1, one);
					mid2=_mm256_sub_epi16(mid2, one);
					mrange0=_mm256_blendv_epi8(mid0, mrange0, mbit0);
					mrange1=_mm256_blendv_epi8(mid1, mrange1, mbit1);
					mrange2=_mm256_blendv_epi8(mid2, mrange2, mbit2);
					
					//set bits
					t0=_mm256_and_si256(mbit0, bitmask);
					t1=_mm256_and_si256(mbit1, bitmask);
					t2=_mm256_and_si256(mbit2, bitmask);
					msym0=_mm256_or_si256(msym0, t0);
					msym1=_mm256_or_si256(msym1, t1);
					msym2=_mm256_or_si256(msym2, t2);
					bitmask=_mm256_srli_epi16(bitmask, 1);
#ifdef ENABLE_VALIDATION
					{
						ALIGN(32) unsigned short probs[48];
						_mm256_store_si256((__m256i*)probs+0, mprob0);
						_mm256_store_si256((__m256i*)probs+1, mprob1);
						_mm256_store_si256((__m256i*)probs+2, mprob2);
						_mm256_store_si256((__m256i*)bits+0, mbit0);
						_mm256_store_si256((__m256i*)bits+1, mbit1);
						_mm256_store_si256((__m256i*)bits+2, mbit2);
						for(int k=0;k<48;++k)
						{
							probs[k]^=bits[k];
							probs[k]-=bits[k];
						}
						val_check(probs, 48);
					}
#endif

					mbit0=_mm256_and_si256(mbit0, one);
					mbit1=_mm256_and_si256(mbit1, one);
					mbit2=_mm256_and_si256(mbit2, one);
					mbit0=_mm256_xor_si256(mbit0, one);
					mbit1=_mm256_xor_si256(mbit1, one);
					mbit2=_mm256_xor_si256(mbit2, one);
					_mm256_store_si256((__m256i*)bits+0, mbit0);
					_mm256_store_si256((__m256i*)bits+1, mbit1);
					_mm256_store_si256((__m256i*)bits+2, mbit2);
#if 1
#ifdef __GNUC__
#pragma GCC unroll 16
#endif
					for(int k=0;k<16;++k)
					{
						unsigned short
							*s0=stats[0]+tidx[16*0+k],
							*s1=stats[1]+tidx[16*1+k],
							*s2=stats[2]+tidx[16*2+k];
						*s0+=((bits[16*0+k]<<16)-*s0+(1<<5>>1))>>5;
						*s1+=((bits[16*1+k]<<16)-*s1+(1<<5>>1))>>5;
						*s2+=((bits[16*2+k]<<16)-*s2+(1<<5>>1))>>5;
						tidx[16*0+k]=2*tidx[16*0+k]+bits[16*0+k];
						tidx[16*1+k]=2*tidx[16*1+k]+bits[16*1+k];
						tidx[16*2+k]=2*tidx[16*2+k]+bits[16*2+k];
					}
#endif
				}
				//if(ky==2&&kx==16)//
				//	printf("");

				//msym0=_mm256_slli_epi16(msym0, 8);
				//msym1=_mm256_slli_epi16(msym1, 8);
				//msym2=_mm256_slli_epi16(msym2, 8);
				//msym0=_mm256_srai_epi16(msym0, 8);
				//msym1=_mm256_srai_epi16(msym1, 8);
				//msym2=_mm256_srai_epi16(msym2, 8);
				_mm256_store_si256((__m256i*)syms+0, msym0);
				_mm256_store_si256((__m256i*)syms+1, msym1);
				_mm256_store_si256((__m256i*)syms+2, msym2);
				for(int k=0;k<remainder;++k, idx+=3)
				{
					short
						NW0	=rows[2*0+1][-1],
						NW1	=rows[2*1+1][-1],
						NW2	=rows[2*2+1][-1],
						N0	=rows[2*0+1][+0],
						N1	=rows[2*1+1][+0],
						N2	=rows[2*2+1][+0],
						W0	=rows[2*0+0][-1],
						W1	=rows[2*1+0][-1],
						W2	=rows[2*2+0][-1],
						*curr0	=rows[2*0+0] +0,
						*curr1	=rows[2*1+0] +0,
						*curr2	=rows[2*2+0] +0;
					
					NW1-=NW0;
					NW2-=NW0;
					N1-=N0;
					N2-=N0;
					W1-=W0;
					W2-=W0;

					int vmin0=N0, vmax0=W0; if(W0<N0)vmin0=W0, vmax0=N0;
					int vmin1=N1, vmax1=W1; if(W1<N1)vmin1=W1, vmax1=N1;
					int vmin2=N2, vmax2=W2; if(W2<N2)vmin2=W2, vmax2=N2;
					int pred0=N0+W0-NW0; CLAMP2(pred0, vmin0, vmax0);
					int pred1=N1+W1-NW1; CLAMP2(pred1, vmin1, vmax1);
					int pred2=N2+W2-NW2; CLAMP2(pred2, vmin2, vmax2);
					
#ifdef DISABLE_PRED
					pred0=0;
#endif

					int c0=(syms[16*0+k]+pred0)<<24>>24;

					pred1+=c0;
					CLAMP2(pred1, -128, 127);
					pred2+=c0;
					CLAMP2(pred2, -128, 127);

#ifdef DISABLE_PRED
					pred1=0;
					pred2=0;
#endif

					int c1=(syms[16*1+k]+pred1)<<24>>24;
					int c2=(syms[16*2+k]+pred2)<<24>>24;
					*curr0=c0;
					*curr1=c1;
					*curr2=c2;

					//if(idx+2>=3*iw*ih)//
					//	LOG_ERROR("");

					image[idx+1]=c0+128;
					image[idx+2]=c1+128;
					image[idx+0]=c2+128;

					guide_check(image, idx/3%iw, idx/3/iw);

					++rows[0];
					++rows[1];
					++rows[2];
					++rows[3];
					++rows[4];
					++rows[5];
				}
			}
		}
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, imsize, fdst);
			fclose(fdst);
		}
		free(image);
	}
	_mm_free(pixels);
	free(srcbuf);
#ifdef _MSC_VER
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
#endif
#ifdef ESTIMATE_SIZE
	if(fwd)
	{
		printf("T%12.2lf\n", (csizes[0]+csizes[1]+csizes[2])/8);
		for(int kc=0;kc<3;++kc)
			printf("%c%12.2lf\n", "YUV"[kc], csizes[kc]/8);
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	}
#endif
#ifdef _MSC_VER
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