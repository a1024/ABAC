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


#if 1
#define NPREDS 14
#define PREDLIST\
	PRED(N+W-NW)\
	PRED(N+W-NW)\
	PRED(N)\
	PRED(W)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED(N+NE-NNE)\
	PRED(W+((NEEE+NEEEEE-N-W)>>3))\
	PRED(W+NW-NWW)\
	PRED(N+NW-NNW)\
	PRED(NE+NEE-NNEEE)\
	PRED((WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED(1)
#endif

#if 0
#define NPREDS2 12
#define PREDLIST2\
	PRED(eN+eW-eNW)\
	PRED(eN)\
	PRED(eW)\
	PRED(3*(eN-eNN)+eNNN)\
	PRED(3*(eW-eWW)+eWWW)\
	PRED(eW+eNE-eN)\
	PRED(eN+eNE-eNNE)\
	PRED(eW+((eNEEE+eNEEEEE-eN-eW)>>3))\
	PRED(eW+eNW-eNWW)\
	PRED(eN+eNW-eNNW)\
	PRED(eNE+eNEE-eNNEEE)\
	PRED((eWWWWW+eWW-eW+eNNN+eN+eNEEEEE)>>2)
#endif
void pred_ols8(Image *src, int fwd)
{
	int vecsize=sizeof(double[3][NPREDS]);
	double *vec=(double*)malloc(vecsize);
	int covsize=sizeof(double[4][NPREDS*NPREDS]);
	double *cov=(double*)malloc(covsize);
	const int wsize=sizeof(double[4][NPREDS+1]);
	double *weights=(double*)malloc(wsize);
	//const int w2size=sizeof(long long[4][NPREDS2+1]);
	//long long *weights2=(long long*)malloc(w2size);
	//const int w3size=sizeof(long long[4][NPREDS2+1]);
	//long long *weights3=(long long*)malloc(w2size);
	int bufsize=(src->iw+8*2)*(int)sizeof(int[6*4*3]);//6 padded rows * 4 channels max * {pixels, residuals1, residuals2}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels||!weights||!vec||!cov)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(cov, 0, covsize);
	memset(vec, 0, vecsize);
	memset(weights, 0, wsize);
	//memset(weights2, 0, w2size);
	//memset(weights3, 0, w3size);
	memset(pixels, 0, bufsize);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
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
#define PRED(EXPR) EXPR,
					PREDLIST
#undef  PRED

				//	N,
				//	W,
				//	3*(N-NN)+NNN,
				//	3*(W-WW)+WWW,
				//	W+NE-N,
				//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
				//	N+W-NW,
				//	N+NE-NNE,
				};
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(preds[0], vmin, vmax);
				double *curr_params=weights+NPREDS*kc;
				double pred1=0;
				for(int k=0;k<NPREDS;++k)
					pred1+=curr_params[k]*preds[k];
#if 0
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
#define PRED(EXPR) EXPR,
					PREDLIST2
#undef  PRED
				};
				long long *currw2=weights2+(NPREDS2+1)*kc;
				long long pred2=currw2[NPREDS2];
				for(int k=0;k<NPREDS2;++k)
					pred2+=currw2[k]*preds2[k];
#endif
#if 0
				int
					e3NNNNN		=rows[5][2+0*4*3],
					e3NNNNWW	=rows[4][2-2*4*3],
					e3NNNNW		=rows[4][2-1*4*3],
					e3NNNN		=rows[4][2+0*4*3],
					e3NNNNE		=rows[4][2+1*4*3],
					e3NNNNEE	=rows[4][2+2*4*3],
					e3NNNNEEEE	=rows[4][2+4*4*3],
					e3NNNWWW	=rows[3][2-3*4*3],
					e3NNNW		=rows[3][2-1*4*3],
					e3NNN		=rows[3][2+0*4*3],
					e3NNNE		=rows[3][2+1*4*3],
					e3NNNEE		=rows[3][2+2*4*3],
					e3NNNEEE	=rows[3][2+3*4*3],
					e3NNNEEEE	=rows[3][2+4*4*3],
					e3NNWWWW	=rows[2][2-4*4*3],
					e3NNWWW		=rows[2][2-3*4*3],
					e3NNWW		=rows[2][2-2*4*3],
					e3NNW		=rows[2][2-1*4*3],
					e3NN		=rows[2][2+0*4*3],
					e3NNE		=rows[2][2+1*4*3],
					e3NNEE		=rows[2][2+2*4*3],
					e3NNEEE		=rows[2][2+3*4*3],
					e3NNEEEE	=rows[2][2+4*4*3],
					e3NWWWW		=rows[1][2-4*4*3],
					e3NWWW		=rows[1][2-3*4*3],
					e3NWW		=rows[1][2-2*4*3],
					e3NW		=rows[1][2-1*4*3],
					e3N		=rows[1][2+0*4*3],
					e3NE		=rows[1][2+1*4*3],
					e3NEE		=rows[1][2+2*4*3],
					e3NEEE		=rows[1][2+3*4*3],
					e3NEEEE		=rows[1][2+4*4*3],
					e3NEEEEE	=rows[1][2+5*4*3],
					e3NEEEEEE	=rows[1][2+6*4*3],
					e3NEEEEEEE	=rows[1][2+7*4*3],
					e3NEEEEEEEE	=rows[1][2+8*4*3],
					e3WWWWWWWWW	=rows[0][2-9*4*3],
					e3WWWWWWWW	=rows[0][2-8*4*3],
					e3WWWWWWW	=rows[0][2-7*4*3],
					e3WWWWWW	=rows[0][2-6*4*3],
					e3WWWWW		=rows[0][2-5*4*3],
					e3WWWW		=rows[0][2-4*4*3],
					e3WWW		=rows[0][2-3*4*3],
					e3WW		=rows[0][2-2*4*3],
					e3W		=rows[0][2-1*4*3];
				int preds3[]=
				{
#define PRED(EXPR) EXPR,
					PREDLIST3
#undef  PRED
				};
				long long *currw3=weights3+(NPREDS2+1)*kc;
				long long pred3=currw3[NPREDS2];
#ifdef PRINT_L1_BOUNDS
				if(b3min>currw2[NPREDS2])b3min=currw2[NPREDS2];
				if(b3max<currw2[NPREDS2])b3max=currw2[NPREDS2];
#endif
				for(int k=0;k<NPREDS3;++k)
				{
					pred3+=currw3[k]*preds3[k];
#ifdef PRINT_L1_BOUNDS
					if(c3min>currw2[k])c3min=currw2[k];
					if(c3max<currw2[k])c3max=currw2[k];
#endif
				}
#endif
				long long predc=CVTFP64_I64(pred1);
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				//if(vmin>NEEEE)vmin=NEEEE;
				//if(vmax<NEEEE)vmax=NEEEE;
				//if(vmin>NW)vmin=NW;
				//if(vmax<NW)vmax=NW;
				//if(vmin>WWW)vmin=WWW;
				//if(vmax<WWW)vmax=WWW;
				CLAMP2(predc, vmin, vmax);
				
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//if(kc==1)//
				//	printf("");

				int curr=src->data[idx];
				int val=(int)predc;
				if(!keyboard[KEY_SHIFT])
				{
					val^=-fwd;
					val+=fwd;
					val+=curr;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
				}
				src->data[idx]=val;
				if(!fwd)
					curr=val;
				rows[0][0]=curr;
				rows[0][1]=curr-(int)pred1;
				//rows[0][2]=curr-(int)pred2s;

				//update
				double *curr_cov=cov+NPREDS*NPREDS*kc, *curr_vec=vec+NPREDS*kc;
				double curr_cholesky[NPREDS*NPREDS];
				double lr=0.003;
				for(int ky2=0, midx=0;ky2<NPREDS;++ky2)
				{
					for(int kx2=0;kx2<NPREDS;++kx2, ++midx)
						curr_cov[midx]+=(preds[kx2]*preds[ky2]-curr_cov[midx])*lr;
				}
				double lval=curr*lr, lr_comp=1-lr;
				for(int k=0;k<NPREDS;++k)
					curr_vec[k]=lval*preds[k]+lr_comp*curr_vec[k];
				int success=1;
				memcpy(curr_cholesky, curr_cov, sizeof(double[NPREDS*NPREDS]));
				for(int k=0;k<NPREDS*NPREDS;k+=NPREDS+1)
					curr_cholesky[k]+=0.0075;
				double sum;
				for(int i=0;i<NPREDS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=curr_cholesky[i*NPREDS+j];
						for(int k=0;k<j;++k)
							sum-=curr_cholesky[i*NPREDS+k]*curr_cholesky[j*NPREDS+k];
						curr_cholesky[i*NPREDS+j]=sum/curr_cholesky[j*NPREDS+j];
					}
					sum=curr_cholesky[i*NPREDS+i];
					for(int k=0;k<i;++k)
						sum-=curr_cholesky[i*NPREDS+k]*curr_cholesky[i*NPREDS+k];
					if(sum<=1e-8)
					{
						success=0;
						break;
					}
					curr_cholesky[i*NPREDS+i]=sqrt(sum);
				}
				if(success)
				{
					for(int i=0;i<NPREDS;++i)
					{
						sum=curr_vec[i];
						for(int j=0;j<i;++j)
							sum-=curr_cholesky[i*NPREDS+j]*curr_params[j];
						curr_params[i]=sum/curr_cholesky[i*NPREDS+i];
					}
					for(int i=NPREDS-1;i>=0;--i)
					{
						sum=curr_params[i];
						for(int j=i+1;j<NPREDS;++j)
							sum-=curr_cholesky[j*NPREDS+i]*curr_params[j];
						curr_params[i]=sum/curr_cholesky[i*NPREDS+i];
					}
				}
			}
		}
	}
	free(pixels);
	free(weights);
	//free(weights2);
	//free(weights3);
}
