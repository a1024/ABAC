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


//	#define ENABLE_GUIDE
//	#define PRINT_SELECTIONS


#include"ac.h"
#ifdef ENABLE_GUIDE
static Image *guide=0;
#endif
#define UPDATE_MASK 3
#define PACK_SIGN(X) ((X)<<1^-((X)<0))
#define UNPACK_SIGN(X) ((X)>>1^-((X)&1))
#define MAG_UPDATE (eN+eW+eNEEE+curr)>>2
static int CG(int N, int W, int NW)
{
	__m128i mN	=_mm_set_epi32(0, 0, 0, N);
	__m128i mW	=_mm_set_epi32(0, 0, 0, W);
	__m128i mNW	=_mm_set_epi32(0, 0, 0, NW);
	__m128i pred=_mm_add_epi32(mN, mW);
	__m128i mmin=_mm_min_epi32(mN, mW);
	pred=_mm_sub_epi32(pred, mNW);
	__m128i mmax=_mm_max_epi32(mN, mW);
	pred=_mm_max_epi32(pred, mmin);
	pred=_mm_min_epi32(pred, mmax);
	ALIGN(16) int apred[4];
	_mm_store_si128((__m128i*)apred, pred);
	return apred[0];
}
static void CGx4(int *nb)
{
	__m128i mN	=_mm_load_si128((__m128i*)nb+0);
	__m128i mW	=_mm_load_si128((__m128i*)nb+1);
	__m128i mNW	=_mm_load_si128((__m128i*)nb+2);
	__m128i pred=_mm_add_epi32(mN, mW);
	__m128i mmin=_mm_min_epi32(mN, mW);
	pred=_mm_sub_epi32(pred, mNW);
	__m128i mmax=_mm_max_epi32(mN, mW);
	pred=_mm_max_epi32(pred, mmin);
	pred=_mm_min_epi32(pred, mmax);
	_mm_store_si128((__m128i*)nb, pred);
}
static int clamp(int vmin, int x, int vmax)
{
	__m128i mx	=_mm_set_epi32(0, 0, 0, x);
	__m128i mmin	=_mm_set_epi32(0, 0, 0, vmin);
	__m128i mmax	=_mm_set_epi32(0, 0, 0, vmax);
	mx=_mm_max_epi32(mx, mmin);
	mx=_mm_min_epi32(mx, mmax);
	ALIGN(16) int ax[4];
	_mm_store_si128((__m128i*)ax, mx);
	return ax[0];
}
static void sort2x4(int *a, int *b)
{
	__m128i ma=_mm_load_si128((__m128i*)a);
	__m128i mb=_mm_load_si128((__m128i*)b);
	__m128i mmin=_mm_min_epi32(ma, mb);
	__m128i mmax=_mm_max_epi32(ma, mb);
	_mm_store_si128((__m128i*)a, mmin);
	_mm_store_si128((__m128i*)b, mmax);
}
static void rct0_fwd(short *comp)
{
	comp[0]-=comp[1];
	comp[2]-=comp[1];
	comp[1]+=(comp[0]+comp[2])>>2;
}
static void rct0_inv(short *comp)
{
	comp[1]-=(comp[0]+comp[2])>>2;
	comp[2]+=comp[1];
	comp[0]+=comp[1];
}
static void rct1_fwd(short *comp)
{
	comp[0]-=comp[1];
	comp[1]+=comp[0]>>1;
	comp[2]-=comp[1];
	comp[1]+=comp[2]*0x5555>>16;//division by 3
}
static void rct1_inv(short *comp)
{
	comp[1]-=comp[2]*0x5555>>16;
	comp[2]+=comp[1];
	comp[1]-=comp[0]>>1;
	comp[0]+=comp[1];
}
static void rct2_fwd(short *comp)
{
	comp[2]-=(87*comp[0]+169*comp[1]+128)>>8;
	comp[0]-=comp[1];
	comp[1]+=(86*comp[0]+29*comp[2]+128)>>8;
}
static void rct2_inv(short *comp)
{
	comp[1]-=(86*comp[0]+29*comp[2]+128)>>8;
	comp[0]+=comp[1];
	comp[2]+=(87*comp[0]+169*comp[1]+128)>>8;
}
#define NPREDS 23
#define NPREDS_G 8
#define NPREDS_B 6
#define NPREDS_R 9
#define NPREDS_RCT0 6
#define NPREDS_RCT1 6
#define NPREDS_RCT2 6
static const short g_pools[]=
{
	//each 3 consequent pools make a channel group
	8, 6, 9,
	2, 2, 2,
	2, 2, 2,
	2, 2, 2,
	0,
};
#define PREDLIST\
	PRED(g)\
	PRED(g-((gN+gW)>>1))\
	PRED(g-gCG)\
	PRED(g-((16*gCG+gN+gW+gNW+gNE-2*(gNN+gWW))>>4))\
	PRED(g-clamp(gmin, (3*(gN+2*gW-gNW)+2*gNE)>>3, gmax))\
	PRED(g-clamp(gmin, (4*(gN+gW)+gW-2*gNW+gNE)>>3, gmax))\
	PRED(g-clamp(gmin, (3*(gN+gW)-2*gNW)>>2, gmax))\
	PRED(g-clamp(gmin, (4*(gN+gW)+gNE-gNW)>>3, gmax))\
	PRED(b)\
	PRED(b-((bN+bW)>>1))\
	PRED(b-bCG)\
	PRED(b-g)\
	PRED(b-clamp(-half, ((bN-gN+bW-gW)>>1)+g, half-1))\
	PRED(b-clamp(-half, bgCG+g, half-1))\
	PRED(r)\
	PRED(r-((rN+rW)>>1))\
	PRED(r-CG(rN, rW, rNW))\
	PRED(r-g)\
	PRED(r-clamp(-half, ((rN-gN+rW-gW)>>1)+g, half-1))\
	PRED(r-clamp(-half, rgCG+g, half-1))\
	PRED(r-((g+b)>>1))\
	PRED(r-clamp(-half, ((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half-1))\
	PRED(r-clamp(-half, (rgbCG+g+b)>>1, half-1))\
	PRED(y0)\
	PRED(y0-CG(y0N, y0W, y0NW))\
	PRED(cb0)\
	PRED(cb0-CG(cb0N, cb0W, cb0NW))\
	PRED(cr0)\
	PRED(cr0-CG(cr0N, cr0W, cr0NW))\
	PRED(y1)\
	PRED(y1-CG(y1N, y1W, y1NW))\
	PRED(cb1)\
	PRED(cb1-CG(cb1N, cb1W, cb1NW))\
	PRED(cr1)\
	PRED(cr1-CG(cr1N, cr1W, cr1NW))\
	PRED(y2)\
	PRED(y2-CG(y2N, y2W, y2NW))\
	PRED(cb2)\
	PRED(cb2-CG(cb2N, cb2W, cb2NW))\
	PRED(cr2)\
	PRED(cr2-CG(cr2N, cr2W, cr2NW))
static const char *prednames[]=
{
#define PRED(X) #X,
	PREDLIST
#undef  PRED
};
//#define CTX_PRED (eN+eW)>>1
#define CTX_UPDATE (eN+eW+eNEEE+val)>>2
static int argmin(int *data, int start, int end)
{
	int kmin=start;
	for(int k=start+1;k<end;++k)
	{
		if(data[kmin]>data[k])
			kmin=k;
	}
	return kmin;
}
static void update_selection(int *csizes, const short *pools, int *selection)
{
	int currsize=csizes[selection[0]]+csizes[selection[1]]+csizes[selection[2]];
	for(int kp=0;kp<_countof(g_pools)-1;kp+=3)
	{
		int idx0=argmin(csizes, pools[kp+0], pools[kp+1]);
		int idx1=argmin(csizes, pools[kp+1], pools[kp+2]);
		int idx2=argmin(csizes, pools[kp+2], pools[kp+3]);
		int csize=csizes[idx0]+csizes[idx1]+csizes[idx2];
		if(currsize>csize)
		{
			selection[0]=idx0;
			selection[1]=idx1;
			selection[2]=idx2;
			currsize=csize;
		}
	}
	memset(csizes, 0, sizeof(int)*pools[12]);
	//for(int kc=0;kc<pools[12];++kc)
	//	csizes[kc]>>=1;
}
int f21_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	
	short pools[_countof(g_pools)]={0};
	int totalcount=0;
	for(int k=0;k<_countof(g_pools);++k)
	{
		int count=g_pools[k];
		pools[k]=totalcount;
		totalcount+=count;
	}
	int depth=image->depth, nlevels=1<<depth, half=nlevels>>1;

	size_t ebufsize=sizeof(short[4*2])*(image->iw+8LL)*totalcount;//4 padded rows * total channel count * {pixels, errors}
	short *pixels=(short*)_mm_malloc(ebufsize, sizeof(__m128i));
	short *mags=(short*)_mm_malloc(ebufsize, sizeof(__m128i));
	int *csizes=(int*)malloc(totalcount*sizeof(int));
	if(!pixels||!mags||!csizes)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, ebufsize);

	int minmag=1<<image->depth>>7;
	*mags=minmag;
	memfill(mags+1, mags, ebufsize-sizeof(short), sizeof(short));
	minmag>>=(image->depth==8)<<1;

	memset(csizes, 0, totalcount*sizeof(int));

	int perm[]={image->nch!=1, 2, 0, 3};
	ptrdiff_t usize=(ptrdiff_t)image->iw*image->ih*image->nch*image->depth>>3;
	int selection[3]=//must be in ascending order
	{
		pools[3]+1,
		pools[4]+1,
		pools[5]+1,
	};
	GolombRiceCoder ec;
	if(fwd)
	{
		DList list;
		dlist_init(&list, 1, 1024, 0);
		gr_enc_init(&ec, &list);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+((image->iw+8LL)*((ky-0LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-1LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-2LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-3LL)&3)+4)*totalcount,
			};
			short *magrows[]=
			{
				mags+((image->iw+8LL)*((ky-0LL)&3)+4)*totalcount,
				mags+((image->iw+8LL)*((ky-1LL)&3)+4)*totalcount,
				mags+((image->iw+8LL)*((ky-2LL)&3)+4)*totalcount,
				mags+((image->iw+8LL)*((ky-3LL)&3)+4)*totalcount,
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int
					gNN	=rows[2][pools[0]+0*totalcount],
					gNW	=rows[1][pools[0]-1*totalcount],
					gN	=rows[1][pools[0]+0*totalcount],
					gNE	=rows[1][pools[0]+1*totalcount],
					gWW	=rows[0][pools[0]-2*totalcount],
					gW	=rows[0][pools[0]-1*totalcount],
					g	=image->data[idx+1];
				int
					bNW	=rows[1][pools[1]-1*totalcount],
					bN	=rows[1][pools[1]+0*totalcount],
					bW	=rows[0][pools[1]-1*totalcount],
					b	=image->data[idx+2];
				int
					rNW	=rows[1][pools[2]-1*totalcount],
					rN	=rows[1][pools[2]+0*totalcount],
					rW	=rows[0][pools[2]-1*totalcount],
					r	=image->data[idx+0];
				short comp[]=
				{
					r, g, b,
					r, g, b,
					r, g, b,
				};
				rct0_fwd(comp+0*3);
				rct1_fwd(comp+1*3);
				rct2_fwd(comp+2*3);
				int
					y0NW	=rows[1][pools[3]-1*totalcount],
					y0N	=rows[1][pools[3]+0*totalcount],
					y0W	=rows[0][pools[3]-1*totalcount],
					y0	=comp[1+0*3],
					cb0NW	=rows[1][pools[4]-1*totalcount],
					cb0N	=rows[1][pools[4]+0*totalcount],
					cb0W	=rows[0][pools[4]-1*totalcount],
					cb0	=comp[2+0*3],
					cr0NW	=rows[1][pools[5]-1*totalcount],
					cr0N	=rows[1][pools[5]+0*totalcount],
					cr0W	=rows[0][pools[5]-1*totalcount],
					cr0	=comp[0+0*3];
				int
					y1NW	=rows[1][pools[6]-1*totalcount],
					y1N	=rows[1][pools[6]+0*totalcount],
					y1W	=rows[0][pools[6]-1*totalcount],
					y1	=comp[1+1*3],
					cb1NW	=rows[1][pools[7]-1*totalcount],
					cb1N	=rows[1][pools[7]+0*totalcount],
					cb1W	=rows[0][pools[7]-1*totalcount],
					cb1	=comp[2+1*3],
					cr1NW	=rows[1][pools[8]-1*totalcount],
					cr1N	=rows[1][pools[8]+0*totalcount],
					cr1W	=rows[0][pools[8]-1*totalcount],
					cr1	=comp[0+1*3];
				int
					y2NW	=rows[1][pools[9]-1*totalcount],
					y2N	=rows[1][pools[9]+0*totalcount],
					y2W	=rows[0][pools[9]-1*totalcount],
					y2	=comp[1+2*3],
					cb2NW	=rows[1][pools[10]-1*totalcount],
					cb2N	=rows[1][pools[10]+0*totalcount],
					cb2W	=rows[0][pools[10]-1*totalcount],
					cb2	=comp[2+2*3],
					cr2NW	=rows[1][pools[11]-1*totalcount],
					cr2N	=rows[1][pools[11]+0*totalcount],
					cr2W	=rows[0][pools[11]-1*totalcount],
					cr2	=comp[0+2*3];
				ALIGN(16) int cgdata[]=
				{
					gN,	bN,	bN-gN,		rN-gN,
					gW,	bW,	bW-gW,		rW-gW,
					gNW,	bNW,	bNW-gNW,	rNW-gNW,

					y0N,	cb0N,	cr0N,		2*rN-(gN+bN),
					y0W,	cb0W,	cr0W,		2*rW-(gW+bW),
					y0NW,	cb0NW,	cr0NW,		2*rNW-(gNW+bNW),

					y1N,	cb1N,	cr1N,		0,
					y1W,	cb1W,	cr1W,		0,
					y1NW,	cb1NW,	cr1NW,		0,

					y2N,	cb2N,	cr2N,		0,
					y2W,	cb2W,	cr2W,		0,
					y2NW,	cb2NW,	cr2NW,		0,
				};
				CGx4(cgdata+0*12);
				CGx4(cgdata+1*12);
				CGx4(cgdata+2*12);
				CGx4(cgdata+3*12);
				int
					gCG	=cgdata[0+0*12],
					bCG	=cgdata[1+0*12],
					bgCG	=cgdata[2+0*12],
					rgCG	=cgdata[3+0*12],
					y0CG	=cgdata[0+1*12],
					cb0CG	=cgdata[1+1*12],
					cr0CG	=cgdata[2+1*12],
					rgbCG	=cgdata[3+1*12],
					y1CG	=cgdata[0+2*12],
					cb1CG	=cgdata[1+2*12],
					cr1CG	=cgdata[2+2*12],

					y2CG	=cgdata[0+3*12],
					cb2CG	=cgdata[1+3*12],
					cr2CG	=cgdata[2+3*12];
				ALIGN(16) int agmin[]={gN, 0, 0, 0}, agmax[]={gW, 0, 0, 0};
				sort2x4(agmin, agmax);
				int
					gmin	=agmin[0],
					gmax	=agmax[0];
				short preds[]=
				{
#define PRED(X) X,
					PREDLIST
#undef  PRED
				};
				memcpy(rows[0], preds, sizeof(short)*totalcount);
				for(int kc=0, kkc=0;kc<pools[12];++kc)
				{
					int
						eN	=magrows[1][kc+0*totalcount],
						eNEEE	=magrows[1][kc+3*totalcount],
						eW	=magrows[0][kc-1*totalcount],
						curr	=rows[0][kc+0*totalcount];
					int mag=((eN+eW)>>1)+1;
					curr=PACK_SIGN(curr);
					
					int nzeros=0, nbypass=0;
					if(kc==selection[kkc])
					{
						gr_enc_track(&ec, curr, mag, &nzeros, &nbypass);
						kkc+=kkc<2;
					}
					else
					{
						//copied from ac.h/gr_enc_track()
						nbypass=FLOOR_LOG2_P1(mag);
						unsigned bypass=curr%mag;
						nzeros=curr/mag;
						unsigned nunused=(1<<nbypass)-mag;//truncated binary code
						if(bypass<nunused)	//emit(bypass, nbypass-1)
							--nbypass;
						else			//emit(bypass+nunused, nbypass)
							bypass+=nunused;

						bypass|=1<<nbypass;//append 1 stop bit
						++nbypass;
					}
					magrows[0][kc+0*totalcount]=MAG_UPDATE;
					csizes[kc]+=nzeros+nbypass;
				}
				rows[0]+=totalcount;
				rows[1]+=totalcount;
				rows[2]+=totalcount;
				rows[3]+=totalcount;
				magrows[0]+=totalcount;
				magrows[1]+=totalcount;
				magrows[2]+=totalcount;
				magrows[3]+=totalcount;
			}
			if((ky&UPDATE_MASK)==UPDATE_MASK)
				update_selection(csizes, pools, selection);
#ifdef PRINT_SELECTIONS
			if(loud)
				printf("%5d  %s  %s  %s\n",
					ky+1,
					prednames[selection[0]],
					prednames[selection[1]],
					prednames[selection[2]]
				);
				//printf("%5d  %2d %2d %2d\n", ky+1, selection[0], selection[1], selection[2]);
#endif
		}
		gr_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		if(loud)
		{
			ptrdiff_t csize=list.nobj;
			t0=time_sec()-t0;
			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("E  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
		dlist_clear(&list);
	}
	else
	{
		LOG_ERROR("Protorype");
		gr_dec_init(&ec, cbuf, cbuf+clen);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+((image->iw+8LL)*((ky-0LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-1LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-2LL)&3)+4)*totalcount,
				pixels+((image->iw+8LL)*((ky-3LL)&3)+4)*totalcount,
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				for(int kkc=0;kkc<3;++kkc)
				{
					int kc=selection[kkc];
					int
						N	=rows[1][kc+0*totalcount],
						W	=rows[0][kc-1*totalcount];
					int mag=((PACK_SIGN(N)+PACK_SIGN(W))>>1)+1;
					int nzeros=0, nbypass=0;
					int curr=gr_dec_track(&ec, mag, &nzeros, &nbypass);
					curr=UNPACK_SIGN(curr);
					rows[0][kc+0*totalcount]=curr;
					csizes[kc]+=nzeros+nbypass;
				}
				rows[0]+=totalcount;
				rows[1]+=totalcount;
				rows[2]+=totalcount;
				rows[3]+=totalcount;
			}
			update_selection(csizes, pools, selection);
		}
		if(loud)
		{
			t0=time_sec()-t0;
			printf("D  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	_mm_free(pixels);
	free(csizes);
	return 0;
}