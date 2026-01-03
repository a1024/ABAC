#include"ebench.h"
#include<stdint.h>
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


#define CTXBITS1 0
#define CTXBITS2 4


#define BOOSTTRAIN 3


#define L1SH 19

#if 1
#define BIAS0 0
#define NPREDS 17
#define PREDLIST\
	PRED(100000, N+W-NW)\
	PRED(200000, N+W-NW)\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED(100000, 3*(N-NN)+NNN)\
	PRED(100000, 3*(W-WW)+WWW)\
	PRED(100000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED(100000, W+((NEEE+NEEEEE-N-W)>>3))\
	PRED( 50000, W+NW-NWW)\
	PRED( 50000, N+NW-NNW)\
	PRED( 50000, NE+NEE-NNEEE)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 40000, NW)\
	PRED( 40000, NE)\
	PRED( 40000, NN)\
	PRED( 40000, WW)
#endif

#if 1
#define NPREDS2 12
#define PREDLIST2\
	PRED(100000, eN+eW-eNW)\
	PRED(100000, eN)\
	PRED(100000, eW)\
	PRED(100000, 3*(eN-eNN)+eNNN)\
	PRED(100000, 3*(eW-eWW)+eWWW)\
	PRED(100000, eW+eNE-eN)\
	PRED(100000, eN+eNE-eNNE)\
	PRED(100000, eW+((eNEEE+eNEEEEE-eN-eW)>>3))\
	PRED( 50000, eW+eNW-eNWW)\
	PRED( 50000, eN+eNW-eNNW)\
	PRED( 50000, eNE+eNEE-eNNEEE)\
	PRED( 50000, (eWWWWW+eWW-eW+eNNN+eN+eNEEEEE)>>2)
#endif

void pred_ols8(Image *src, int fwd)
{
	int amin[]=
	{
		-(1<<src->depth[0]>>1),
		-(1<<src->depth[1]>>1),
		-(1<<src->depth[2]>>1),
		-(1<<src->depth[3]>>1),
	};
	int amax[]=
	{
		(1<<src->depth[0]>>1)-1,
		(1<<src->depth[1]>>1)-1,
		(1<<src->depth[2]>>1)-1,
		(1<<src->depth[3]>>1)-1,
	};
	int invdist=((1<<16)+g_dist-1)/g_dist;
	const int wsize=sizeof(long long[4][1<<CTXBITS1][NPREDS+1]);
	long long *weights=(long long*)malloc(wsize);
	const int w2size=sizeof(long long[4][1<<CTXBITS2][NPREDS2+1]);
	long long *weights2=(long long*)malloc(w2size);
	int bufsize=(src->iw+8*2)*(int)sizeof(short[6*4*3]);//6 padded rows * 4 channels max * {pixels, residuals1, residuals2}
	short *pixels=(short*)malloc(bufsize);
	if(!pixels||!weights||!weights2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
//	{
//		int j=0;
//#define PRED(W0, EXPR) weights[0*(NPREDS+1)+j]=weights[1*(NPREDS+1)+j]=weights[2*(NPREDS+1)+j]=weights[3*(NPREDS+1)+j]=W0; ++j;
//		PREDLIST
//#undef  PRED
//		weights[0*(NPREDS+1)+NPREDS]=weights[1*(NPREDS+1)+NPREDS]=weights[2*(NPREDS+1)+NPREDS]=weights[3*(NPREDS+1)+NPREDS]=BIAS0;
//	}
	//FILLMEM(weights, (1<<L1SH)/NPREDS, wsize, sizeof(long long));
	memset(weights, 0, wsize);
	memset(weights2, 0, w2size);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-1LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-2LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-3LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-4LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-5LL+6)%6)+8)*4-1)*3,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				rows[4]+=3;
				rows[5]+=3;
				if(!src->depth[kc])
					continue;
				int
					NNNNN		=rows[5][+0*4*3],
					NNNNWW		=rows[4][-2*4*3],
					NNNNW		=rows[4][-1*4*3],
					NNNN		=rows[4][+0*4*3],
					NNNNE		=rows[4][+1*4*3],
					NNNNEE		=rows[4][+2*4*3],
					NNNNEEEE	=rows[4][+4*4*3],
					NNNWWW		=rows[3][-3*4*3],
					NNNW		=rows[3][-1*4*3],
					NNN		=rows[3][+0*4*3],
					NNNE		=rows[3][+1*4*3],
					NNNEE		=rows[3][+2*4*3],
					NNNEEE		=rows[3][+3*4*3],
					NNNEEEE		=rows[3][+4*4*3],
					NNWWWW		=rows[2][-4*4*3],
					NNWWW		=rows[2][-3*4*3],
					NNWW		=rows[2][-2*4*3],
					NNW		=rows[2][-1*4*3],
					NN		=rows[2][+0*4*3],
					NNE		=rows[2][+1*4*3],
					NNEE		=rows[2][+2*4*3],
					NNEEE		=rows[2][+3*4*3],
					NNEEEE		=rows[2][+4*4*3],
					NWWWW		=rows[1][-4*4*3],
					NWWW		=rows[1][-3*4*3],
					NWW		=rows[1][-2*4*3],
					NW		=rows[1][-1*4*3],
					N		=rows[1][+0*4*3],
					NE		=rows[1][+1*4*3],
					NEE		=rows[1][+2*4*3],
					NEEE		=rows[1][+3*4*3],
					NEEEE		=rows[1][+4*4*3],
					NEEEEE		=rows[1][+5*4*3],
					NEEEEEE		=rows[1][+6*4*3],
					NEEEEEEE	=rows[1][+7*4*3],
					NEEEEEEEE	=rows[1][+8*4*3],
					WWWWWWWWW	=rows[0][-9*4*3],
					WWWWWWWW	=rows[0][-8*4*3],
					WWWWWWW		=rows[0][-7*4*3],
					WWWWWW		=rows[0][-6*4*3],
					WWWWW		=rows[0][-5*4*3],
					WWWW		=rows[0][-4*4*3],
					WWW		=rows[0][-3*4*3],
					WW		=rows[0][-2*4*3],
					W		=rows[0][-1*4*3];
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(preds[0], vmin, vmax);
				//long long *currw=weights+(NPREDS+1)*((1<<CTXBITS1)*kc+((N+W)/2<<CTXBITS1>>src->depth[kc]&((1<<CTXBITS1)-1)));
				long long *currw=weights+(NPREDS+1)*kc;
				long long pred1=currw[NPREDS];
				for(int k=0;k<NPREDS;++k)
					pred1+=currw[k]*preds[k];
				int
					eNNNNN		=rows[5][1+0*4*3],
					eNNNNWW		=rows[4][1-2*4*3],
					eNNNNW		=rows[4][1-1*4*3],
					eNNNN		=rows[4][1+0*4*3],
					eNNNNE		=rows[4][1+1*4*3],
					eNNNNEE		=rows[4][1+2*4*3],
					eNNNNEEEE	=rows[4][1+4*4*3],
					eNNNWWW		=rows[3][1-3*4*3],
					eNNNW		=rows[3][1-1*4*3],
					eNNN		=rows[3][1+0*4*3],
					eNNNE		=rows[3][1+1*4*3],
					eNNNEE		=rows[3][1+2*4*3],
					eNNNEEE		=rows[3][1+3*4*3],
					eNNNEEEE	=rows[3][1+4*4*3],
					eNNWWWW		=rows[2][1-4*4*3],
					eNNWWW		=rows[2][1-3*4*3],
					eNNWW		=rows[2][1-2*4*3],
					eNNW		=rows[2][1-1*4*3],
					eNN		=rows[2][1+0*4*3],
					eNNE		=rows[2][1+1*4*3],
					eNNEE		=rows[2][1+2*4*3],
					eNNEEE		=rows[2][1+3*4*3],
					eNNEEEE		=rows[2][1+4*4*3],
					eNWWWW		=rows[1][1-4*4*3],
					eNWWW		=rows[1][1-3*4*3],
					eNWW		=rows[1][1-2*4*3],
					eNW		=rows[1][1-1*4*3],
					eN		=rows[1][1+0*4*3],
					eNE		=rows[1][1+1*4*3],
					eNEE		=rows[1][1+2*4*3],
					eNEEE		=rows[1][1+3*4*3],
					eNEEEE		=rows[1][1+4*4*3],
					eNEEEEE		=rows[1][1+5*4*3],
					eNEEEEEE	=rows[1][1+6*4*3],
					eNEEEEEEE	=rows[1][1+7*4*3],
					eNEEEEEEEE	=rows[1][1+8*4*3],
					eWWWWWWWWW	=rows[0][1-9*4*3],
					eWWWWWWWW	=rows[0][1-8*4*3],
					eWWWWWWW	=rows[0][1-7*4*3],
					eWWWWWW		=rows[0][1-6*4*3],
					eWWWWW		=rows[0][1-5*4*3],
					eWWWW		=rows[0][1-4*4*3],
					eWWW		=rows[0][1-3*4*3],
					eWW		=rows[0][1-2*4*3],
					eW		=rows[0][1-1*4*3];
				int preds2[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST2
#undef  PRED
				};
				long long *currw2=weights2+(NPREDS2+1)*((1<<CTXBITS2)*kc+((N+W)/2<<CTXBITS2>>src->depth[kc]&((1<<CTXBITS2)-1)));
				//long long *currw2=weights2+(NPREDS2+1)*kc;
				long long pred2=currw2[NPREDS2];
				for(int k=0;k<NPREDS2;++k)
					pred2+=currw2[k]*preds2[k];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

				pred1+=pred2<<3;//[sic]
				pred1+=1<<L1SH>>1;
				pred1>>=L1SH;
				pred2+=1<<L1SH>>1;
				pred2>>=L1SH;

				long long pred2s=pred1+pred2, predc=pred2s;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(predc, vmin, vmax);
				
				int curr=src->data[idx];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)predc;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx]=curr;

						curr=g_dist*curr+(int)predc;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+(int)predc;
						CLAMP2(curr, amin[kc], amax[kc]);

						src->data[idx]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						int error=curr-(int)predc;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx]=error;
					}
					else
					{
						curr+=(int)predc;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx]=curr;
					}
				}
				rows[0][0]=curr;
				rows[0][1]=curr-(int)pred1;
				rows[0][2]=curr-(int)pred2s;

				//update
				int e=(curr>pred1)-(curr<pred1);
				currw[NPREDS]+=e;//bias
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];//coeffs

				e=(curr>pred2s)-(curr<pred2s);
				currw2[NPREDS2]+=e;//bias
				for(int k=0;k<NPREDS2;++k)
					currw2[k]+=e*preds2[k];//coeffs
			}
		}
	}
	free(pixels);
	free(weights);
	free(weights2);
}
void pred_ols8_crct(Image *src, int fwd)
{
	int amin[]=
	{
		-(1<<src->depth[0]>>1),
		-(1<<src->depth[1]>>1),
		-(1<<src->depth[2]>>1),
		-(1<<src->depth[3]>>1),
	};
	int amax[]=
	{
		(1<<src->depth[0]>>1)-1,
		(1<<src->depth[1]>>1)-1,
		(1<<src->depth[2]>>1)-1,
		(1<<src->depth[3]>>1)-1,
	};
	int invdist=((1<<16)+g_dist-1)/g_dist;
	const int wsize=sizeof(int[4][1<<CTXBITS1][NPREDS+1]);
	int *weights=(int*)malloc(wsize);
	const int w2size=sizeof(int[4][1<<CTXBITS2][NPREDS2+1]);
	int *weights2=(int*)malloc(w2size);
	int bufsize=(src->iw+8*2)*(int)sizeof(short[6*4*3]);//6 padded rows * 4 channels max * {pixels, residuals1, residuals2}
	short *pixels=(short*)malloc(bufsize);
	if(fwd)
		src->rct=crct_analysis(src);
	const unsigned char *combination=rct_combinations[src->rct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	if(!pixels||!weights||!weights2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
//	{
//		int j=0;
//#define PRED(W0, EXPR) weights[0*(NPREDS+1)+j]=weights[1*(NPREDS+1)+j]=weights[2*(NPREDS+1)+j]=weights[3*(NPREDS+1)+j]=W0; ++j;
//		PREDLIST
//#undef  PRED
//		weights[0*(NPREDS+1)+NPREDS]=weights[1*(NPREDS+1)+NPREDS]=weights[2*(NPREDS+1)+NPREDS]=weights[3*(NPREDS+1)+NPREDS]=BIAS0;
//	}
	//FILLMEM(weights, (1<<L1SH)/NPREDS, wsize, sizeof(long long));
	memset(weights, 0, wsize);
	memset(weights2, 0, w2size);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-1LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-2LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-3LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-4LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-5LL+6)%6)+8)*4-1)*3,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0;
			int yuv[]=
			{
				src->data[idx+yidx],
				src->data[idx+uidx],
				src->data[idx+vidx],
			};
			for(int kc=0;kc<4;++kc)
			{
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				rows[4]+=3;
				rows[5]+=3;
				if(!src->depth[kc])
					continue;
				int
					NNNNN		=rows[5][+0*4*3],
					NNNNWW		=rows[4][-2*4*3],
					NNNNW		=rows[4][-1*4*3],
					NNNN		=rows[4][+0*4*3],
					NNNNE		=rows[4][+1*4*3],
					NNNNEE		=rows[4][+2*4*3],
					NNNNEEEE	=rows[4][+4*4*3],
					NNNWWW		=rows[3][-3*4*3],
					NNNW		=rows[3][-1*4*3],
					NNN		=rows[3][+0*4*3],
					NNNE		=rows[3][+1*4*3],
					NNNEE		=rows[3][+2*4*3],
					NNNEEE		=rows[3][+3*4*3],
					NNNEEEE		=rows[3][+4*4*3],
					NNWWWW		=rows[2][-4*4*3],
					NNWWW		=rows[2][-3*4*3],
					NNWW		=rows[2][-2*4*3],
					NNW		=rows[2][-1*4*3],
					NN		=rows[2][+0*4*3],
					NNE		=rows[2][+1*4*3],
					NNEE		=rows[2][+2*4*3],
					NNEEE		=rows[2][+3*4*3],
					NNEEEE		=rows[2][+4*4*3],
					NWWWW		=rows[1][-4*4*3],
					NWWW		=rows[1][-3*4*3],
					NWW		=rows[1][-2*4*3],
					NW		=rows[1][-1*4*3],
					N		=rows[1][+0*4*3],
					NE		=rows[1][+1*4*3],
					NEE		=rows[1][+2*4*3],
					NEEE		=rows[1][+3*4*3],
					NEEEE		=rows[1][+4*4*3],
					NEEEEE		=rows[1][+5*4*3],
					NEEEEEE		=rows[1][+6*4*3],
					NEEEEEEE	=rows[1][+7*4*3],
					NEEEEEEEE	=rows[1][+8*4*3],
					WWWWWWWWW	=rows[0][-9*4*3],
					WWWWWWWW	=rows[0][-8*4*3],
					WWWWWWW		=rows[0][-7*4*3],
					WWWWWW		=rows[0][-6*4*3],
					WWWWW		=rows[0][-5*4*3],
					WWWW		=rows[0][-4*4*3],
					WWW		=rows[0][-3*4*3],
					WW		=rows[0][-2*4*3],
					W		=rows[0][-1*4*3];
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(preds[0], vmin, vmax);
				//int *currw=weights+(NPREDS+1)*((1<<CTXBITS1)*kc+((N+W)/2<<CTXBITS1>>src->depth[kc]&((1<<CTXBITS1)-1)));
				int *currw=weights+(NPREDS+1)*kc;
				int pred1=currw[NPREDS];
				for(int k=0;k<NPREDS;++k)
					pred1+=currw[k]*preds[k];
				int
					eNNNNN		=rows[5][1+0*4*3],
					eNNNNWW		=rows[4][1-2*4*3],
					eNNNNW		=rows[4][1-1*4*3],
					eNNNN		=rows[4][1+0*4*3],
					eNNNNE		=rows[4][1+1*4*3],
					eNNNNEE		=rows[4][1+2*4*3],
					eNNNNEEEE	=rows[4][1+4*4*3],
					eNNNWWW		=rows[3][1-3*4*3],
					eNNNW		=rows[3][1-1*4*3],
					eNNN		=rows[3][1+0*4*3],
					eNNNE		=rows[3][1+1*4*3],
					eNNNEE		=rows[3][1+2*4*3],
					eNNNEEE		=rows[3][1+3*4*3],
					eNNNEEEE	=rows[3][1+4*4*3],
					eNNWWWW		=rows[2][1-4*4*3],
					eNNWWW		=rows[2][1-3*4*3],
					eNNWW		=rows[2][1-2*4*3],
					eNNW		=rows[2][1-1*4*3],
					eNN		=rows[2][1+0*4*3],
					eNNE		=rows[2][1+1*4*3],
					eNNEE		=rows[2][1+2*4*3],
					eNNEEE		=rows[2][1+3*4*3],
					eNNEEEE		=rows[2][1+4*4*3],
					eNWWWW		=rows[1][1-4*4*3],
					eNWWW		=rows[1][1-3*4*3],
					eNWW		=rows[1][1-2*4*3],
					eNW		=rows[1][1-1*4*3],
					eN		=rows[1][1+0*4*3],
					eNE		=rows[1][1+1*4*3],
					eNEE		=rows[1][1+2*4*3],
					eNEEE		=rows[1][1+3*4*3],
					eNEEEE		=rows[1][1+4*4*3],
					eNEEEEE		=rows[1][1+5*4*3],
					eNEEEEEE	=rows[1][1+6*4*3],
					eNEEEEEEE	=rows[1][1+7*4*3],
					eNEEEEEEEE	=rows[1][1+8*4*3],
					eWWWWWWWWW	=rows[0][1-9*4*3],
					eWWWWWWWW	=rows[0][1-8*4*3],
					eWWWWWWW	=rows[0][1-7*4*3],
					eWWWWWW		=rows[0][1-6*4*3],
					eWWWWW		=rows[0][1-5*4*3],
					eWWWW		=rows[0][1-4*4*3],
					eWWW		=rows[0][1-3*4*3],
					eWW		=rows[0][1-2*4*3],
					eW		=rows[0][1-1*4*3];
				int preds2[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST2
#undef  PRED
				};
				int *currw2=weights2+(NPREDS2+1)*((1<<CTXBITS2)*kc+((N+W)/2<<CTXBITS2>>src->depth[kc]&((1<<CTXBITS2)-1)));
				//int *currw2=weights2+(NPREDS2+1)*kc;
				int pred2=currw2[NPREDS2];
				for(int k=0;k<NPREDS2;++k)
					pred2+=currw2[k]*preds2[k];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

				//if((int64_t)pred2<<4!=(int32_t)(pred2<<4))
				//	LOG_ERROR("");
				pred1+=pred2<<3;//[sic]
				pred1+=1<<L1SH>>1;
				pred1>>=L1SH;
				pred2+=1<<L1SH>>1;
				pred2>>=L1SH;

				int pred2s=pred1+pred2, predc=pred2s;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(predc, vmin, vmax);
				if(kc)
				{
					predc+=offset;
					CLAMP2(predc, amin[kc], amax[kc]);
				}
				
				int curr=yuv[kc];
				if(g_dist>1)
				{
#if 0
					int dist=g_dist+src->depth[kc]-FLOOR_LOG2(abs(eW)+abs(eN)+1);
					//int dist=g_dist/(abs(eW)+abs(eN)+1)+1;
					if(fwd)
					{
						curr-=predc;
						curr/=dist;
						src->data[idx+kc]=curr;
					}
					else
						curr=src->data[idx+kc];
					curr=dist*curr+predc;
#endif
#if 1
					if(fwd)
					{
						curr-=predc;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx+kc]=curr;
					}
					else
						curr=src->data[idx+kc];
					curr=g_dist*curr+predc;
#endif

					CLAMP2(curr, amin[kc], amax[kc]);
					yuv[kc]=curr;
				}
				else
				{
					if(fwd)
					{
						int error=curr-predc;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx+kc]=error;
					}
					else
					{
						curr=src->data[idx+kc]+predc;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						yuv[kc]=curr;
					}
				}
				curr-=offset;
				rows[0][0]=curr;
				rows[0][1]=curr-(int)pred1;
				rows[0][2]=curr-(int)pred2s;

				//update
				int e=(curr>pred1)-(curr<pred1);
				currw[NPREDS]+=e;//bias
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];//coeffs

				e=(curr>pred2s)-(curr<pred2s);
				currw2[NPREDS2]+=e;//bias
				for(int k=0;k<NPREDS2;++k)
					currw2[k]+=e*preds2[k];//coeffs

				offset=kc?(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2:yuv[0]&vfromy;
			}
			if(!fwd)
			{
				src->data[idx+yidx]=yuv[0];
				src->data[idx+uidx]=yuv[1];
				src->data[idx+vidx]=yuv[2];
			}
		}
	}
	free(pixels);
	free(weights);
	free(weights2);
}
