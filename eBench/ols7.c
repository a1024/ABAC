#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


#define OLS7_NPARAMS 3
#define OLS7_PADX 16	//64
#define OLS7_PADY 16	//multiple of 4
static void learn(const int *ctx, int target, long long *acc, long long *params)
{
#if OLS7_NPARAMS==3
	//8 bit only
	acc[0]+=((long long)ctx[0]*ctx[0]-acc[0])*7>>3;//sqa
	acc[1]+=((long long)ctx[1]*ctx[1]-acc[1])*7>>3;//sqb
	acc[2]+=((long long)ctx[2]*ctx[2]-acc[2])*7>>3;//sqc
	acc[3]+=((long long)ctx[0]*ctx[1]-acc[3])*7>>3;//ab
	acc[4]+=((long long)ctx[1]*ctx[2]-acc[4])*7>>3;//bc
	acc[5]+=((long long)ctx[2]*ctx[0]-acc[5])*7>>3;//ca
	acc[6]+=((long long)ctx[0]*target-acc[6])*7>>3;//at
	acc[7]+=((long long)ctx[1]*target-acc[7])*7>>3;//bt
	acc[8]+=((long long)ctx[2]*target-acc[8])*7>>3;//ct
	params[0]=(acc[1]*acc[2]-acc[4]*acc[4])*acc[6]+(acc[5]*acc[4]-acc[3]*acc[2])*acc[7]+(acc[3]*acc[4]-acc[1]*acc[5])*acc[8];//param0 = [sqb*sqc-bc^2  ac*bc-ab*sqc  ab*bc-sqb*ac] . [at bt ct]
	params[1]=(acc[5]*acc[4]-acc[3]*acc[2])*acc[6]+(acc[0]*acc[2]-acc[5]*acc[5])*acc[7]+(acc[3]*acc[5]-acc[0]*acc[4])*acc[8];//param1 = [ac*bc-ab*sqc  sqa*sqc-ac^2  ab*ac-sqa*bc] . [at bt ct]
	params[2]=(acc[3]*acc[4]-acc[5]*acc[1])*acc[6]+(acc[3]*acc[5]-acc[0]*acc[4])*acc[7]+(acc[0]*acc[1]-acc[3]*acc[3])*acc[8];//param2 = [ab*bc-ac*sqb  ab*ac-sqa*bc  sqa*sqb-ab^2] . [at bt ct]
	params[3]=acc[0]*acc[1]*acc[2]+2*acc[3]*acc[4]*acc[5]-acc[0]*acc[4]*acc[4]-acc[1]*acc[5]*acc[5]-acc[2]*acc[3]*acc[3];//den = sqa*sqb*sqc + 2*ab*bc*ac - sqa*bc^2 - sqb*ac^2 - sqc*ab^2
#else
	acc[0]+=((long long)ctx[0]*ctx[0]-acc[0])*7>>3;//sqa
	acc[1]+=((long long)ctx[1]*ctx[1]-acc[1])*7>>3;//sqb
	acc[2]+=((long long)ctx[0]*ctx[1]-acc[1])*7>>3;//sab
	acc[3]+=((long long)ctx[0]*target-acc[1])*7>>3;//at
	acc[4]+=((long long)ctx[1]*target-acc[1])*7>>3;//bt
	params[0]=acc[1]*acc[3]-acc[2]*acc[4];//param0 = sqb*at-sab*bt
	params[1]=acc[0]*acc[4]-acc[2]*acc[3];//param1 = sqa*bt-sab*at
	params[2]=acc[0]*acc[1]-acc[2]*acc[2];//den = sqa*sqb-sab^2
#endif
}
void pred_ols7(Image *src, int fwd)
{
	long long acc[3][OLS7_NPARAMS*(OLS7_NPARAMS+1)/2+OLS7_NPARAMS]={0};//{sqa, sqb, sab, at, bt}
	long long params[3][OLS7_NPARAMS+1]={0};//{param0, param1, den}

	int bufsize=(src->iw+OLS7_PADX*2)*(int)sizeof(short[OLS7_PADY*4]);//16 padded rows * 4 channels max
	short *pixels=(short*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	memset(pixels, 0, bufsize);
	for(int ky=0;ky<src->ih;++ky)
	{
		__m256i rowstride=_mm256_set1_epi64x(sizeof(short[4]));
		ALIGN(32) short *rows[OLS7_PADY]={0};
		for(int k=0;k<_countof(rows);++k)
			rows[k]=pixels+((src->iw+OLS7_PADX*2LL)*((ky-(size_t)k)%OLS7_PADY)+OLS7_PADX)*4;
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int
					NN	=rows[2][kc+0*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4];
				int ctx[]=
				{
					N,
					W,
					NW,
				};
				//int ctx2[]=
				//{
				//	N,
				//	W,
				//	NW,
				//	//4*(N+W-NW),
				//	//2*(N+W),
				//};
				//int ctx3[2];
				int pred;
				//if(ky==10&&kx==10&&!kc)//
				//	printf("");

				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				int half=1<<src->depth[kc]>>1;
				vmin+=vmin<vmax;
				vmax-=vmax>vmin;
				//vmin=MAXVAR(vmin-1, -half);
				//vmax=MINVAR(vmax+1, half-1);
				pred=(N+NE)/2+W-NW;
				UPDATE_MIN(pred, vmax);
				UPDATE_MAX(pred, vmin);
				//MEDIAN3_32(pred, N, W, (N+W-NW+2*N-NN+2*W-WW+2+NE)/4);
				//CLAMP3_32(pred, (N+W-NW+2*N-NN+2*W-WW+2+NE)/4, N, W, NE);
				//MEDIAN3_32(pred, N, W, N+W-NW);
				//pred=(2*pred+2*N-NN+2*W-WW+2)>>2;
#if 0
#if OLS7_NPARAMS==3
				if(params[kc][3])
				{
					pred=(int)((params[kc][0]*ctx[0]+params[kc][1]*ctx[1]+params[kc][2]*ctx[2])/params[kc][3]);
					CLAMP3_32(pred, pred, N, W, NE);
				}
#else
				if(params[kc][2])
				{
					pred=(int)((params[kc][0]*ctx[0]+params[kc][1]*ctx[1])/params[kc][2]);
					CLAMP3_32(pred, pred, N, W, NE);
				}
#endif
				else
					//pred=0;
					MEDIAN3_32(pred, N, W, N+W-NW);
#endif
				int idx=(src->iw*ky+kx)<<2|kc;
				int curr=0;
				if(fwd)
				{
					curr=src->data[idx];
					pred=curr-pred;
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
				}
				else
				{
					pred+=src->data[idx];
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
					curr=pred;
				}
				src->data[idx]=pred;
				rows[0][kc]=curr;

				//learn(ctx, curr, acc[kc], params[kc]);
				//acc[kc][0]+=((long long)ctx[0]*ctx[0]-acc[kc][0])*7>>3;//sqa
				//acc[kc][1]+=((long long)ctx[1]*ctx[1]-acc[kc][1])*7>>3;//sqb
				//acc[kc][2]+=((long long)ctx[0]*ctx[1]-acc[kc][1])*7>>3;//sab
				//acc[kc][3]+=((long long)ctx[0]*curr-acc[kc][1])*7>>3;//at
				//acc[kc][4]+=((long long)ctx[1]*curr-acc[kc][1])*7>>3;//bt
				//params[kc][0]=acc[kc][1]*acc[kc][3]-acc[kc][2]*acc[kc][4];//param0 = sqb*at-sab*bt
				//params[kc][1]=acc[kc][0]*acc[kc][4]-acc[kc][2]*acc[kc][3];//param1 = sqa*bt-sab*at
				//params[kc][2]=acc[kc][0]*acc[kc][1]-acc[kc][2]*acc[kc][2];//den = sqa*sqb-sab^2
			}
			for(int k=0;k<sizeof(rows)/sizeof(__m256i);++k)
			{
				__m256i mr=_mm256_load_si256((__m256i*)rows+k);
				mr=_mm256_add_epi64(mr, rowstride);
				_mm256_store_si256((__m256i*)rows+k, mr);
			}
		}
	}
	free(pixels);
}