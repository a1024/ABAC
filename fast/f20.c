#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


//	#define PROFILE_CSIZE
//	#define ENABLE_GUIDE


#include"ac.h"
#ifdef ENABLE_GUIDE
static Image *guide=0;
#endif
#define CTX_PRED (eN+eW)>>1
#define CTX_UPDATE (eN+eW+eNEEE+val)>>2
int f20_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif

	size_t ebufsize=sizeof(short[4*8])*(image->iw+8LL);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)_mm_malloc(ebufsize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	int minmag=1<<image->depth>>7;
	pixels[0]=0;
	pixels[1]=0;
	pixels[2]=0;
	pixels[3]=0;
	pixels[4]=minmag;
	pixels[5]=minmag;
	pixels[6]=minmag;
	pixels[7]=minmag;
	memfill(pixels+8, pixels, ebufsize-sizeof(short[8]), sizeof(short[8]));
	__m128i mone=_mm_set1_epi16(1);
	__m128i mminmag=_mm_load_si128((__m128i*)pixels);
	if(image->depth==8)
	{
		mminmag=_mm_srli_epi16(mminmag, 2);
		minmag>>=2;
	}
	//minmag>>=image->depth==16;//X
	//memset(pixels, 0, ebufsize);
	int perm[]={image->nch!=1, 2, 0, 3};
	ptrdiff_t usize=(ptrdiff_t)image->iw*image->ih*image->nch*image->depth>>3;
	ALIGN(16) short pred[8]={0}, mag[8]={0};
	if(fwd)
	{
#ifdef PROFILE_CSIZE
		ptrdiff_t csizes[4*3]={0};//{unary, bypass}
#endif
		GolombRiceCoder ec;
		DList list;
		dlist_init(&list, 1, 1024, 0);
		gr_enc_init(&ec, &list);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				short *curr=rows[0];
				memcpy(curr, image->data+idx, sizeof(short)*image->nch);
				__m128i mNW	=_mm_load_si128((__m128i*)rows[1]-1);
				__m128i mN	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i mW	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i mp=_mm_add_epi16(mN, mW);
				__m128i mmin=_mm_min_epi16(mN, mW);
				__m128i mmax=_mm_max_epi16(mN, mW);
				__m128i mpc=_mm_srai_epi16(mp, 1);
				mp=_mm_sub_epi16(mp, mNW);
				mpc=_mm_max_epi16(mpc, mminmag);
				mp=_mm_max_epi16(mp, mmin);
				mp=_mm_min_epi16(mp, mmax);
				mpc=_mm_add_epi16(mpc, mone);
				_mm_store_si128((__m128i*)pred, mp);
				_mm_store_si128((__m128i*)mag, mpc);
				if(image->nch>=3)
				{
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
				}
				for(int kc0=0;kc0<image->nch;++kc0)
				{
					int
						kc	=perm[kc0],
					//	NW	=rows[1][kc-1*8+0],
					//	N	=rows[1][kc+0*8+0],
					//	NE	=rows[1][kc+1*8+0],
					//	WW	=rows[0][kc-2*8+0],
					//	W	=rows[0][kc-1*8+0],
					//	eNW	=rows[1][kc-1*8+4],
						eN	=rows[1][kc+0*8+4],
					//	eNE	=rows[1][kc+1*8+4],
					//	eNEE	=rows[1][kc+2*8+4],
						eNEEE	=rows[1][kc+3*8+4],
					//	eNEEEE	=rows[1][kc+4*8+4],
					//	eWW	=rows[0][kc-2*8+4],
						eW	=rows[0][kc-1*8+4];
					//int pred;
					//MEDIAN3_32(pred, N, W, N+W-NW);
#if 0
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					//UPDATE_MIN(vmin, NE);
					//UPDATE_MAX(vmax, NE);
					//int pred=W+((5*(N-NW)+NE-WW)>>3);
					//int pred=(N+W)>>1;
					int pred=N+W-NW;
					pred=CLAMP(vmin, pred, vmax);
					//pred=(62*pred+N+W)>>6;
#endif
					//int magnitude=CTX_PRED;
					//int magnitude=!ky?eW:(!kx?eN:CTX_PRED);
					//if(kc0>1)
					//	magnitude=(magnitude+rows[0][0+0*8+4])>>1;
					//UPDATE_MAX(magnitude, minmag);
					//int magnitude=!kx||!ky?depth-2:(eN+eW)>>1;
					//int magnitude=MAXVAR(eN, eW);
					//magnitude-=magnitude>0;//X
					//magnitude+=magnitude<depth;//X
					//if(idx==1155516&&kc0==2)//
					//if(ky==10&&kx==10)//
					//	printf("");

					int val=curr[kc]-pred[kc];
					val=val<<1^-(val<0);
					curr[kc+4]=CTX_UPDATE;
#ifdef PROFILE_CSIZE
					++magnitude;
					csizes[0<<2|kc]+=val/magnitude;
					++csizes[1<<2|kc];
					int bypass=val%magnitude;
					int nbypass=floor_log2_32(magnitude)+1;
					csizes[2<<2|kc]+=nbypass-(bypass<(1<<nbypass)-magnitude);
					--magnitude;
#endif
					gr_enc(&ec, val, mag[kc+4]);
				}
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
			}
		}
		gr_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		if(loud)
		{
			ptrdiff_t csize=list.nobj;
			t0=time_sec()-t0;
			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
#ifdef PROFILE_CSIZE
			printf("C         unary         stop_bit          bypass\n");
			double total_unary=0, total_stop=0, total_bypass=0;
			for(int kc=0;kc<image->nch;++kc)
			{
				double size_unary=csizes[0<<2|kc]/8., size_stop=csizes[1<<2|kc]/8., size_bypass=csizes[2<<2|kc]/8.;
				printf("%d %16.2lf %16.2lf %16.2lf\n", kc, size_unary, size_stop, size_bypass);
				total_unary+=size_unary;
				total_stop+=size_stop;
				total_bypass+=size_bypass;
			}
			printf("T %16.2lf %16.2lf %16.2lf\n", total_unary, total_stop, total_bypass);
#endif
			printf("E  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
		dlist_clear(&list);
	}
	else
	{
		GolombRiceCoder ec;
		gr_dec_init(&ec, cbuf, cbuf+clen);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				short *curr=rows[0];
				__m128i mNW	=_mm_load_si128((__m128i*)rows[1]-1);
				__m128i mN	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i mW	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i mp=_mm_add_epi16(mN, mW);
				__m128i mmin=_mm_min_epi16(mN, mW);
				__m128i mmax=_mm_max_epi16(mN, mW);
				__m128i mpc=_mm_srai_epi16(mp, 1);
				mp=_mm_sub_epi16(mp, mNW);
				mpc=_mm_max_epi16(mpc, mminmag);
				mp=_mm_max_epi16(mp, mmin);
				mp=_mm_min_epi16(mp, mmax);
				mpc=_mm_add_epi16(mpc, mone);
				_mm_store_si128((__m128i*)pred, mp);
				_mm_store_si128((__m128i*)mag, mpc);
				//__m128i mNW	=_mm_load_si128((__m128i*)rows[1]-1);
				//__m128i mN	=_mm_load_si128((__m128i*)rows[1]+0);
				//__m128i mW	=_mm_load_si128((__m128i*)rows[0]-1);
				//__m128i mp=_mm_add_epi16(mN, mW);
				//__m128i mmin=_mm_min_epi16(mN, mW);
				//__m128i mmax=_mm_max_epi16(mN, mW);
				//mp=_mm_sub_epi16(mp, mNW);
				//mp=_mm_max_epi16(mp, mmin);
				//mp=_mm_min_epi16(mp, mmax);
				//_mm_store_si128((__m128i*)pred, mp);
				for(int kc0=0;kc0<image->nch;++kc0)
				{
					int
						kc	=perm[kc0],
					//	NW	=rows[1][kc-1*8+0],
					//	N	=rows[1][kc+0*8+0],
					//	NE	=rows[1][kc+1*8+0],
					//	WW	=rows[0][kc-2*8+0],
					//	W	=rows[0][kc-1*8+0],
					//	eNW	=rows[1][kc-1*8+4],
						eN	=rows[1][kc+0*8+4],
					//	eNE	=rows[1][kc+1*8+4],
					//	eNEE	=rows[1][kc+2*8+4],
						eNEEE	=rows[1][kc+3*8+4],
					//	eNEEEE	=rows[1][kc+4*8+4],
					//	eWW	=rows[0][kc-2*8+4],
						eW	=rows[0][kc-1*8+4];
					//int pred;
					//MEDIAN3_32(pred, N, W, N+W-NW);
#if 0
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					//UPDATE_MIN(vmin, NE);
					//UPDATE_MAX(vmax, NE);
					//int pred=W+((5*(N-NW)+NE-WW)>>3);
					//int pred=(N+W)>>1;
					int pred=N+W-NW;
					pred=CLAMP(vmin, pred, vmax);
					//pred=(62*pred+N+W)>>6;
#endif
					//int magnitude=CTX_PRED;
					//int magnitude=!ky?eW:(!kx?eN:CTX_PRED);
					//int magnitude=!kx||!ky?depth-2:(eN+eW)>>1;
					//int magnitude=MAXVAR(eN, eW);
					//if(kc0>1)
					//	magnitude=(magnitude+rows[0][0+0*8+4])>>1;
					//UPDATE_MAX(magnitude, minmag);

					int val=gr_dec(&ec, mag[kc+4]);
					curr[kc+4]=CTX_UPDATE;
					val=val>>1^-(val&1);
					curr[kc]=val+pred[kc];
				}
				short *c2=image->data+idx;
				memcpy(c2, curr, sizeof(short)*image->nch);
				if(image->nch>=3)
				{
					c2[1]-=(c2[0]+c2[2])>>2;
					c2[2]+=c2[1];
					c2[0]+=c2[1];
				}
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(c2, guide->data+idx, sizeof(short)*image->nch))
				{
					short c0[4]={0};
					memcpy(c0, guide->data+idx, sizeof(short)*image->nch);
					c0[0]-=c0[1];
					c0[2]-=c0[1];
					c0[1]+=(c0[0]+c0[2])>>2;
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
					LOG_ERROR("Guide error IDX %d/%d", idx, image->nch*image->iw*image->ih);
					printf("");//
				}
#endif
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
			}
		}
		if(loud)
		{
			t0=time_sec()-t0;
			printf("D  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	_mm_free(pixels);
	return 0;
}