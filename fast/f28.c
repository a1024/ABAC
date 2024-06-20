#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


#include"ac.h"
#define BLOCKSIZE 64


#define ORCT_NPARAMS 8//not counting permutation

#define ORCT_NITER 100
#define ORCT_DELTAGROUP 1
#define ORCT_WATCHDOGTIMEOUT (ORCT_NPARAMS*16)

#define OCRT_PARAMBITS 7
#define ORCT_ONE (1<<OCRT_PARAMBITS)


#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9
#if 1
#define WG_NPREDS	11
#define WG_PREDLIST\
	WG_PRED(132, N+W-NW)\
	WG_PRED(176, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(176, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(176, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(176, (W+NEEE)/2)
#endif
#if 0
#define WG_NPREDS 33
#define WG_PREDLIST\
	WG_PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	WG_PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>8))\
	WG_PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>8))\
	WG_PRED(N+(int)((-eNN*0x0DFLL-eN*0x051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x102)>>8))\
	WG_PRED(3*(N-NN)+NNN)\
	WG_PRED((N+W)>>1)\
	WG_PRED(N+W-NW)\
	WG_PRED((W+NEE)>>1)\
	WG_PRED((3*W+NEEE)>>2)\
	WG_PRED((3*(3*W+NE+NEE)-10*N+2)/5)\
	WG_PRED((3*(3*W+NE+NEE)-10*N)/5)\
	WG_PRED((4*N-2*NN+NW+NE)>>2)\
	WG_PRED(N+NE-NNE-eNNE)\
	WG_PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	WG_PRED(W+((eW-eWW)>>1))\
	WG_PRED(paper_GAP)\
	WG_PRED(calic_GAP)\
	WG_PRED(N+W-((NW+NN+WW+NE)>>2))\
	WG_PRED(((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12)\
	WG_PRED(3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW)\
	WG_PRED(2*(W+NE-N)-(WW+NNEE-NN))\
	WG_PRED((2*W+NEE-N)>>1)\
	WG_PRED(NW+NWW-NNWWW)\
	WG_PRED((14*NE-(NNEE+NNNEE+NNEEE))/11)\
	WG_PRED((NEEE+NEEEE)>>1)\
	WG_PRED((NNNEEEE+NNEEE)>>1)\
	WG_PRED(NNEEEE)\
	WG_PRED((NNWWWW+NNNWWWW)>>1)\
	WG_PRED((WWW+WWWW)>>1)\
	WG_PRED((N+NN)>>1)\
	WG_PRED((NE+NNEE)>>1)\
	WG_PRED((NE+NNE+NEE+NNEE)>>2)\
	WG_PRED(ols)
#endif

#define OCH_COUNT 4
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, x1, x2, y1, y2;
	
	int *pixels, bufsize;
	int *hist, histsize;

	//AC
	int tlevels, clevels, statssize;
	unsigned *stats;
	DList list;
	const unsigned char *decstart, *decend;

	//aux
	double bestsize;

	//WG
	int wg_weights[4*WG_NPREDS];
	int wg_errors[4*WG_NPREDS];
} ThreadArgs;
static const char *slic5_orct_permutationnames[]=
{
	"rgb",//0
	"grb",//1
	"gbr",//2
	"rbg",//3
	"bgr",//4
	"brg",//5
};
static void orct_print_compact(const char *params)
{
	printf("[%s", slic5_orct_permutationnames[(unsigned char)params[ORCT_NPARAMS]]);
	for(int k=0;k<ORCT_NPARAMS;++k)
	{
		int val=params[k];
		printf("%c%02X", val<0?'-':'+', abs(val));
	}
	printf("]");
}
static void orct_unpack_permutation(unsigned char p, unsigned char *permutation)
{
	int temp=0;
	switch(p)
	{
	case 0:temp=0x020100;break;//rgb
	case 1:temp=0x020001;break;//grb
	case 2:temp=0x000201;break;//gbr
	case 3:temp=0x010200;break;//rbg
	case 4:temp=0x000102;break;//bgr
	case 5:temp=0x010002;break;//brg
	default:
		LOG_ERROR("Invalid RGB permutation");
		return;
	}
	memcpy(permutation, &temp, sizeof(char[3]));
}
static void orct_custom_getmatrix(const char *params, double *matrix, int fwd)
{
	unsigned char per[4]={0};
	orct_unpack_permutation(params[ORCT_NPARAMS], per);
	memset(matrix, 0, sizeof(double[9]));
	matrix[0]=1;
	matrix[4]=1;
	matrix[8]=1;
	for(int k=0;k<3;++k)
	{
		double vrtx[]=
		{
			matrix[k+3*(fwd?per[0]:0)],
			matrix[k+3*(fwd?per[1]:1)],
			matrix[k+3*(fwd?per[2]:2)],
		};
		if(fwd)
		{
			vrtx[0]+=(params[0]*vrtx[1]+params[1]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
			vrtx[1]+=(params[2]*vrtx[0]+params[3]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
			vrtx[2]+=(params[4]*vrtx[0]+params[5]*vrtx[1])*(1./(1<<OCRT_PARAMBITS));
			vrtx[1]+=(params[6]*vrtx[0]+params[7]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
		}
		else
		{
			vrtx[1]-=(params[6]*vrtx[0]+params[7]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
			vrtx[2]-=(params[4]*vrtx[0]+params[5]*vrtx[1])*(1./(1<<OCRT_PARAMBITS));
			vrtx[1]-=(params[2]*vrtx[0]+params[3]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
			vrtx[0]-=(params[0]*vrtx[1]+params[1]*vrtx[2])*(1./(1<<OCRT_PARAMBITS));
		}
		matrix[k+3*(fwd?0:per[0])]=vrtx[0];
		matrix[k+3*(fwd?1:per[1])]=vrtx[1];
		matrix[k+3*(fwd?2:per[2])]=vrtx[2];
	}
}
static void orct_fwd(int *comp, const char *params, const unsigned char *permutation)
{
	int c2[3]=
	{
		comp[permutation[0]],
		comp[permutation[1]],
		comp[permutation[2]],
	};
	c2[0]+=(params[0]*c2[1]+params[1]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[1]+=(params[2]*c2[0]+params[3]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[2]+=(params[4]*c2[0]+params[5]*c2[1]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[1]+=(params[6]*c2[0]+params[7]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	memcpy(comp, c2, sizeof(c2));
}
static void orct_inv(int *comp, const char *params, const unsigned char *permutation)
{
	int c2[3];
	memcpy(c2, comp, sizeof(c2));
	c2[1]-=(params[6]*c2[0]+params[7]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[2]-=(params[4]*c2[0]+params[5]*c2[1]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[1]-=(params[2]*c2[0]+params[3]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	c2[0]-=(params[0]*c2[1]+params[1]*c2[2]+((1<<OCRT_PARAMBITS)>>1))>>OCRT_PARAMBITS;
	comp[permutation[0]]=c2[0];
	comp[permutation[1]]=c2[1];
	comp[permutation[2]]=c2[2];
}
static double orct_calcloss(Image const *src, int x1, int x2, int y1, int y2, int *pixels, int bufsize, int *hist, const char *params)
{
	int depth=src->depth+2, nlevels=1<<depth, half=nlevels>>1, mask=nlevels-1;
	__m128i mhalf=_mm_set1_epi32(half);
	__m128i mmask=_mm_set1_epi32(mask);
	unsigned char permutation[3]={0};

	orct_unpack_permutation(params[8], permutation);
	memset(hist, 0, nlevels*sizeof(int[3]));
	memset(pixels, 0, bufsize);
	for(int ky=y1;ky<y2;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			pixels+((BLOCKSIZE+16LL)*((ky-0LL)&1)+8LL)*4,
			pixels+((BLOCKSIZE+16LL)*((ky-1LL)&1)+8LL)*4,
		};
		ALIGN(16) int result[4]={0};
		for(int kx=x1;kx<x2;++kx)
		{
			int
				idx	=src->nch*(src->iw*ky+kx),
				*pNW	=rows[1]-1*4,
				*pN	=rows[1]+0*4,
				*pW	=rows[0]-1*4,
				*pcurr	=rows[0]+0*4;
			pcurr[0]=src->data[idx+0];
			pcurr[1]=src->data[idx+1];
			pcurr[2]=src->data[idx+2];
			orct_fwd(pcurr, params, permutation);
			{
				__m128i NW	=_mm_load_si128((__m128i*)pNW);
				__m128i N	=_mm_load_si128((__m128i*)pN);
				__m128i W	=_mm_load_si128((__m128i*)pW);
				__m128i curr	=_mm_load_si128((__m128i*)pcurr);
				__m128i vmin=_mm_min_epi32(N, W);
				__m128i vmax=_mm_max_epi32(N, W);
				__m128i pred=_mm_add_epi32(_mm_sub_epi32(N, NW), W);
				pred=_mm_max_epi32(pred, vmin);
				pred=_mm_min_epi32(pred, vmax);

				curr=_mm_sub_epi32(curr, pred);
				curr=_mm_add_epi32(curr, mhalf);
				curr=_mm_and_si128(curr, mmask);
				_mm_store_si128((__m128i*)result, curr);
				++hist[0<<depth|result[0]];
				++hist[1<<depth|result[1]];
				++hist[2<<depth|result[2]];
			}
			rows[0]+=4;
			rows[1]+=4;
		}
	}
	{
		double csize;
		int res=(x2-x1)*(y2-y1);
		double csizes[3]={0};
		for(int kc=0;kc<3;++kc)
		{
			int *curr_hist=hist+((size_t)kc<<depth);
			for(int ks=0;ks<nlevels;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
		}
		csize=(csizes[0]+csizes[1]+csizes[2])/8;
		return csize;
	}
}
static double orct_optimize(Image const *src, int x1, int x2, int y1, int y2, char *params, int loud)
{
	static const char *rctnames[]=
	{
		"NONE",
		"SubGreen",
		"JPEG2000",
		"YCbCr_R_v1",
		"A710",
		"YCbCr_R_v3",
		"YCbCr_R_v4",
		"YCbCr_R_v5",
		"YCbCr_R_v6",
		"YCbCr_R_v7",
		"Pei09",
	};
	static const char paramsets[]=//permutations are brute forced
	{
		0,	0,//0	RCT_NONE defeats the purpose of optimization?
		0,	0,
		0,	0,
		0,	0,
		0,//rgb
		
		-ORCT_ONE,	0,//1	RCT_SubG
		0,		0,
		0,		-ORCT_ONE,
		0,		0,
		0,//rgb

		-ORCT_ONE,	0,//2	RCT_JPEG2000
		0,		0,
		0,		-ORCT_ONE,
		ORCT_ONE>>1,	ORCT_ONE>>1,
		0,//rgb

		-ORCT_ONE,	0,//3	RCT_YCbCr-R_v1
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE>>1,
		0,//rgb

		-ORCT_ONE,	0,//4	RCT_A710
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		-(ORCT_ONE>>3),	ORCT_ONE>>4,
		0,//rgb

		-ORCT_ONE,	0,//5	RCT_YCbCr-R_v3
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		ORCT_ONE>>3,	ORCT_ONE>>4,
		0,//rgb

		-ORCT_ONE,	0,//6	RCT_YCbCr-R_v4 "fair luma"
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		(ORCT_ONE+1)/3,
		0,//rgb

		-ORCT_ONE,	0,//7	RCT_YCbCr-R_v5
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE*3>>3,
		0,//rgb

		-ORCT_ONE,	0,//8	RCT_YCbCr-R_v6
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE*7>>4,
		0,//rgb

		-ORCT_ONE,	0,//9	RCT_YCbCr-R_v7
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		-(ORCT_ONE>>5),	ORCT_ONE*10>>5,
		0,//rgb

		-(ORCT_ONE*87>>8),	-(ORCT_ONE*169>>8),	//10	RCT_Pei09
		0,			0,
		0,			-ORCT_ONE,
		ORCT_ONE*29>>8,		ORCT_ONE*86>>8,
		0,//rgb
	};
	double loss_init, loss_bestsofar, loss_prev, loss_curr=0;
	double t_start=time_sec();
	int maxdepth=src->depth+2, maxlevels=1<<maxdepth;
	int bufsize=(BLOCKSIZE+16LL)*sizeof(int[2*4]);//2 rows * 4 channels
	int *pixels=(int*)_mm_malloc(bufsize, sizeof(__m128i));
	int *hist=(int*)malloc(maxlevels*sizeof(int[3]));
	if(!pixels||!hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	char params2[]=
	{
		-ORCT_ONE,	0,		//-1	0		JPEG2000 RCT
		0,		0,		//0	0
		0,		-ORCT_ONE,	//0	-1
		ORCT_ONE>>1,	ORCT_ONE>>1,	//1/2	1/2
		0,//rgb
	};
#define CALC_LOSS(L) L=orct_calcloss(src, x1, x2, y1, y2, pixels, bufsize, hist, params2)

	int best_init=0, bestp=0;
	loss_bestsofar=INFINITY;
	for(int k=0;k<(_countof(paramsets)/(ORCT_NPARAMS+1));++k)//brute force initial RCT
	{
		memcpy(params2, paramsets+(ORCT_NPARAMS+1)*k, sizeof(params2));
		for(int k2=1;k2<6;++k2)//brute force the permutations
		{
			params2[ORCT_NPARAMS]=k2;
			CALC_LOSS(loss_prev);
			if(loss_bestsofar>loss_prev)
				loss_bestsofar=loss_prev, best_init=k, bestp=k2;
		}
	}
	memcpy(params2, paramsets+(ORCT_NPARAMS+1)*best_init, sizeof(params2));
	params2[ORCT_NPARAMS]=bestp;
	loss_init=loss_bestsofar;
	memcpy(params, params2, sizeof(params2));//save
	if(loud)
		printf("init  rgb->%s(%d)  %s(%d)  UC %d->%lf bytes\n",
			slic5_orct_permutationnames[bestp], bestp,
			rctnames[best_init], best_init,
			((x2-x1)*(y2-y1)*src->nch*src->depth+7)>>3,
			loss_bestsofar
		);

	for(int it=0, checkpoint=-1;it<ORCT_NITER;++it)
	{
		double losses[3]={0};
		char params3[ORCT_NPARAMS+1];
		const int deviation2=ORCT_ONE>>4;

		memcpy(params3, params2, sizeof(params3));

		//slider
		int vals[3]=
		{
			CLAMP(-ORCT_ONE, params2[0]-deviation2, ORCT_ONE-1),
			CLAMP(-ORCT_ONE, params2[0]+deviation2, ORCT_ONE-1),
		};
		params2[0]=vals[0], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[0]);
		params2[0]=vals[1], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[1]);
		for(;vals[0]<vals[1]-1;)
		{
			vals[2]=(vals[0]+vals[1])>>1;
			params2[0]=vals[2], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[2]);
			if(losses[2]>=losses[0]&&losses[2]>=losses[1])
			{
				if(losses[1]>losses[0])
					loss_curr=losses[0], params2[0]=vals[0];
				else
					loss_curr=losses[1], params2[0]=vals[1];
				params2[1]=-ORCT_ONE-params2[0];
				break;
			}
			loss_curr=losses[2];
			if(losses[1]>losses[0])
				losses[1]=losses[2], vals[1]=vals[2];
			else
				losses[0]=losses[2], vals[0]=vals[2];
		}

		//luma coefficient 1 and 2 (can take any value)
		const int deviation=ORCT_ONE>>4;
		for(int k=2;k<ORCT_NPARAMS;++k)
		{
			char *p=params2+k;
			switch(k>>1)
			{
			case 2://curr+[-1, 1/2]
				vals[0]=*p-deviation;
				vals[1]=*p+(deviation>>1);
				break;
			case 1://curr+[-1/2, 0.99]
			case 3:
				vals[0]=*p-(deviation>>1);
				vals[1]=*p+deviation;
				break;
			}

			vals[0]=CLAMP(-ORCT_ONE, vals[0], ORCT_ONE-1);
			vals[1]=CLAMP(-ORCT_ONE, vals[1], ORCT_ONE-1);
			*p=vals[0]; CALC_LOSS(losses[0]);
			*p=vals[1]; CALC_LOSS(losses[1]);
			for(;vals[0]<vals[1];)
			{
				vals[2]=(vals[0]+vals[1])>>1;
				*p=vals[2]; CALC_LOSS(losses[2]);
				if(losses[2]>=losses[0]&&losses[2]>=losses[1])
				{
					if(losses[1]>losses[0])
						loss_curr=losses[0], *p=vals[0];
					else
						loss_curr=losses[1], *p=vals[1];
					break;
				}
				loss_curr=losses[2];
				if(losses[1]>losses[0])
					losses[1]=losses[2], vals[1]=vals[2];
				else
					losses[0]=losses[2], vals[0]=vals[2];
			}
		}

		//brute force the permutations
		bestp=params2[ORCT_NPARAMS];
		for(int k=0;k<6;++k)
		{
			if(k==bestp)
				continue;
			params2[ORCT_NPARAMS]=k;
			CALC_LOSS(loss_prev);
			if(loss_bestsofar>loss_prev)
				loss_bestsofar=loss_prev, bestp=k;
		}
		params2[ORCT_NPARAMS]=bestp;
		
		if(loss_bestsofar>loss_curr)
		{
			memcpy(params, params2, sizeof(params2));
			loss_bestsofar=loss_curr;
			checkpoint=it;
		}
		if(loud)
			printf("chk %4d  %4d/%4d  curr %16lf  best %16lf  elapsed %12lf sec\r",
				checkpoint+1, it+1, ORCT_NITER, loss_curr, loss_bestsofar, time_sec()-t_start
			);
		if(!memcmp(params2, params3, sizeof(params2))||it-checkpoint>=5)
			break;
	}
#undef  CALC_LOSS
	if(loud)
	{
		const char chnames[]="rgb";
		unsigned char p[3]={0};
		double matrix[9]={0};

		printf("\n");
		printf("CR improvement %lf -> %lf  %+6.2lf%%\n", loss_init, loss_bestsofar, 100.*(1-loss_init/loss_bestsofar));
		printf("Final RCT params:\n");
		//for(int k=0;k<ORCT_NPARAMS;++k)
		//	printf("  %12g,%c", (double)params[k]/ORCT_ONE, k&1?'\n':' ');
		//printf("  p%d  rgb->%s\n", params[ORCT_NPARAMS], slic5_orct_permutationnames[params[ORCT_NPARAMS]]);
		orct_unpack_permutation(params[ORCT_NPARAMS], p);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[0]], (double)params[0]/ORCT_ONE, chnames[p[1]], (double)params[1]/ORCT_ONE, chnames[p[2]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[1]], (double)params[2]/ORCT_ONE, chnames[p[0]], (double)params[3]/ORCT_ONE, chnames[p[2]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[2]], (double)params[4]/ORCT_ONE, chnames[p[0]], (double)params[5]/ORCT_ONE, chnames[p[1]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[1]], (double)params[6]/ORCT_ONE, chnames[p[0]], (double)params[7]/ORCT_ONE, chnames[p[2]]);

		printf("RCT matrces:\n");
		orct_custom_getmatrix(params, matrix, 1);
		for(int k=0;k<9;++k)
		{
			if(!(k%3))
				printf("%.3s", "V  Y  U  "+k);
			printf("%15g, ", matrix[k]);
			if(!((k+1)%3))
				printf("%c\n", "RGB"[k/3]);
		}
		printf("\n");
		orct_custom_getmatrix(params, matrix, 0);
		for(int k=0;k<9;++k)
		{
			if(!(k%3))
				printf("%.3s", "R  G  B  "+k);
			printf("%15g, ", matrix[k]);
			if(!((k+1)%3))
				printf("%c\n", "YUV"[k/3]);
		}

		printf("Raw data: ");
		orct_print_compact(params);
		printf("\n");
	}

	_mm_free(pixels);
	free(hist);
	return loss_bestsofar;
}

static void wg_init(int *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const int *weights,
	int **rows, const int stride, int kc2, int depth,
	const int *perrors, int *preds
)
{
	double pred2=0, wsum=0;
	int j=0, pred;
#ifdef WG_T47
	short
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNNWWW	=rows[3][kc2-3*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NNNEE	=rows[3][kc2+2*stride+0],
		NNNEEEE	=rows[3][kc2+4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNWWW	=rows[2][kc2-3*stride+0],
		NNWW	=rows[2][kc2-2*stride+0],
		NNW	=rows[2][kc2-1*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEE	=rows[2][kc2+2*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NNEEEE	=rows[2][kc2+4*stride+0],
		NWW	=rows[1][kc2-2*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		NEEEE	=rows[1][kc2+4*stride+0],
		WWWW	=rows[0][kc2-4*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNN	=rows[2][kc2+0*stride+1],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eWW	=rows[0][kc2-2*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int diff=(dy-dx)<<8>>depth, diff2=(d45-d135)<<8>>depth, diff3=NE-NW;
	int paper_GAP, calic_GAP;

	if(dy+dx>32)//[sic]
		paper_GAP=(int)(((long long)dx*N+(long long)dy*W)/((long long)dy+dx));
	else if(diff>12)
		paper_GAP=(N+2*W)/3;
	else if(diff<-12)
		paper_GAP=(2*N+W)/3;
	else
		paper_GAP=(N+W)>>1;

	if(diff2>32)
		paper_GAP+=diff3>>2;
	else if(diff2>16)
		paper_GAP+=diff3*3>>4;
	else if(diff2>=-16)
		paper_GAP+=diff3>>3;
	else if(diff2>=-32)
		paper_GAP+=diff3>>4;

	if(diff>80)
		calic_GAP=W;
	else if(diff<-80)
		calic_GAP=N;
	else if(diff>32)
		calic_GAP=(2*N+6*W+NE-NW)>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
	else if(diff>8)
		calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	//c2	[3/8  5/8  3/16  -3/16]
	else if(diff<-32)
		calic_GAP=(6*N+2*W+NE-NW)>>3;		//c3	[3/4  1/4  1/8  -1/8]
	else if(diff<-8)
		calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	//c4	[5/8  3/8  3/16  -3/16]
	else
		calic_GAP=(((N+W)<<1)+NE-NW)>>2;	//c5	[1/2  1/2  1/4  -1/4]
#else
	int
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eW	=rows[0][kc2-1*stride+1];
#endif

#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
	}
	pred2/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
static void wg_update(int curr, const int *preds, int *perrors, int *weights)
{
#ifdef WG_UPDATE
	double wsum=0;
	int kbest=0, ebest=0;
	//int kworst=0, eworst=0;
#endif
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
#ifdef WG_UPDATE
		if(!k||ebest>e2)
			kbest=k, ebest=e2;
		//if(!k||eworst<e2)
		//	kworst=k, eworst=e2;
		
		//{
		//	//https://stackoverflow.com/questions/11644441/fast-inverse-square-root-on-x64
		//	double t1=e2+1, t2=t1*0.5;
		//	long long t3=*(long long*)&t1;
		//	t3=0x5FE6EB50C7B537A9-(t3>>1);
		//	t1=*(double*)&t3;
		//	t1=t1*(1.5-(t2*t1*t1));
		//	wsum+=weights[k]+=t1;
		//}
		//wsum+=weights[k]+=1./sqrt(e2+1);
		//weights[k]+=1./(e2+1);
#endif
	}
#ifdef WG_UPDATE
	//wsum=1/wsum;
	//for(int k=0;k<WG_NPREDS;++k)
	//	weights[k]*=wsum;

	//weights[kbest]+=1./256;
	++weights[kbest];
	//if(weights[kbest]>22)//352/16.
	if(weights[kbest]>352)
	{
		for(int k=0;k<WG_NPREDS;++k)
			weights[k]>>=1;
			//weights[k]*=0.5;
	}
#endif
}

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
static void quantize_pixel(int val, int *token, int *bypass, int *nbits)
{
	if(val<(1<<CONFIG_EXP))
	{
		*token=val;//token
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
#define CDFSTRIDE 16
static void update_CDF(const int *hist, unsigned *CDF, int tlevels)
{
	int sum=hist[tlevels], c=0;
	for(int ks=0;ks<tlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<PROB_BITS)-tlevels)/sum)+ks;
		c+=freq;
	}
	CDF[tlevels]=1<<PROB_BITS;
}
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void block_thread(void *param)
{
	ArithmeticCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	char rct_params[ORCT_NPARAMS+1]={0}, rct_per[4]={0};
	int depth=image->depth+2, nlevels=1<<depth, half=nlevels>>1, mask=nlevels-1;
	int res=(args->x2-args->x1)*(args->y2-args->y1);
	int nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1, chsize=nctx*cdfstride;

	if(args->fwd)//encode
	{
		args->bestsize=orct_optimize(image, args->x1, args->x2, args->y1, args->y2, rct_params, args->loud);

		dlist_init(&args->list, 1, 1024, 0);
		ac_enc_init(&ec, &args->list);

		for(int k=0;k<ORCT_NPARAMS;++k)
			ac_enc_bypass(&ec, rct_params[k]+128, 8);
		ac_enc_bypass_NPOT(&ec, rct_params[ORCT_NPARAMS], 6);
	}
	else//decode
	{
		ac_dec_init(&ec, args->decstart, args->decend);
		
		for(int k=0;k<ORCT_NPARAMS;++k)
			rct_params[k]=ac_dec_bypass(&ec, 8)-128;
		rct_params[ORCT_NPARAMS]=ac_dec_bypass_NPOT(&ec, 6);
	}
	orct_unpack_permutation(rct_params[ORCT_NPARAMS], rct_per);
	for(int kc=0;kc<image->nch;++kc)
	{
		int *curr_hist=args->hist+chsize*kc;
		unsigned *curr_CDF=args->stats+chsize*kc;
		
		*curr_hist=1;//init bypass
		memfill(curr_hist+1, curr_hist, sizeof(int)*(args->tlevels-1LL), sizeof(int));
		curr_hist[args->tlevels]=args->tlevels;
		update_CDF(curr_hist, curr_CDF, args->tlevels);
		memfill(curr_hist+cdfstride, curr_hist, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
		memfill(curr_CDF+cdfstride, curr_CDF, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	}
	for(int kc=0;kc<image->nch;++kc)
		wg_init(args->wg_weights+WG_NPREDS*kc);
	memset(args->wg_errors, 0, sizeof(args->wg_errors));
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int wg_preds[WG_NPREDS]={0};
		int token=0, bypass=0, nbits=0;
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int
				idx	=image->nch*(image->iw*ky+kx),
				*NNN	=rows[3]+0*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NW=N=W;
					NN=N;
					NNE=NE;
				}
				NNN=NN;
			}
			if(kx<=2)
			{
				if(kx<=1)
				{
					if(!kx)
						NW=W=N;
					WW=W;
				}
				WWW=WW;
			}
			if(kx>=image->iw-3)
			{
				if(kx>=image->iw-2)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEE=NE;
				}
				NEEE=NEE;
			}
			if(args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					yuv[0]=image->data[idx+0];
					yuv[1]=image->data[idx+1];
					yuv[2]=image->data[idx+2];
					orct_fwd(yuv, rct_params, rct_per);
					kc=3;
				}
				if(kc<image->nch)
					yuv[kc]=image->data[kc];
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				int kc2=kc<<1;
				int pred=0, error, sym;
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<8>>depth,
					vy=(abs(W[kc2]-NW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<8>>depth;
				int qeN=FLOOR_LOG2(vy+1);
				int qeW=FLOOR_LOG2(vx+1);
				int cidx=cdfstride*(nctx*kc+args->clevels*MINVAR(qeN, 8)+MINVAR(qeW, 8));
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;

				pred=wg_predict(
					args->wg_weights+WG_NPREDS*kc,
					rows, 4*2, kc2, depth,
					args->wg_errors+WG_NPREDS*kc, wg_preds
				);

				//if(ky==13&&kx==623)//
				//if(ky==80&&kx==510&&kc==2)//
				//if(ky==16&&kx==640&&kc==1)//
				//if(ky==0&&kx==258&&kc==1)//
				//if(ky==17&&kx==458&&kc==2)//
				//if(ky==255&&kx==767&&kc==0)//
				//if(ky==0&&kx==125&&kc==0)//
				//if(ky==0&&kx==102&&kc==1)//
				//if(ky==1&&kx==111&&kc==1)//
				//	printf("");

				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
					{
						int upred=half-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[ch]));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^(sym>>31);//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);

					ac_enc(&ec, token, curr_CDF);
					if(nbits)
						ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits
				}
				else
				{
					token=ac_dec(&ec, curr_CDF, args->tlevels);//try ac_dec_packedsign()
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						int lsb, msb;

						sym-=1<<CONFIG_EXP;
						lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
					{
						int upred=half-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-half));
#endif
						}
						else
						{
							error=sym-upred;
							negmask=-(pred>0);
						}
						error^=negmask;
						error-=negmask;
					}
					curr[kc2+0]=yuv[kc]=error+pred;
					curr[kc2+1]=error;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if(curr_hist[args->tlevels]>=args->tlevels&&!(curr_hist[args->tlevels]&(CDFSTRIDE-1)))
					update_CDF(curr_hist, curr_CDF, args->tlevels);
				if(curr_hist[args->tlevels]>=6144)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist[ks]>>=1;
					curr_hist[args->tlevels]=sum;
				}
				wg_update(curr[kc2], wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
			}
			if(!args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					orct_inv(yuv, rct_params, rct_per);
					image->data[idx+0]=yuv[0];
					image->data[idx+1]=yuv[1];
					image->data[idx+2]=yuv[2];
					kc=3;
				}
				if(kc<image->nch)
					image->data[kc]=yuv[kc];
#ifdef ENABLE_GUIDE
				if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
				{
					short orig[4]={0};
					memcpy(orig, guide->data+idx, image->nch*sizeof(short));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac_enc_flush(&ec);
}
int f28_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int ncores=query_cpu_cores();
	int xblocks=(image->iw+BLOCKSIZE-1)/BLOCKSIZE;
	int yblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE;
	int nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	int coffset=sizeof(int)*nblocks;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	int tlevels=0, clevels=0, statssize=0;
	double esize=0;
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
		start=array_append(data, 0, 1, coffset, 1, 0, 0);
	}
	else//integrity check
	{
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(start!=(ptrdiff_t)clen)
			LOG_ERROR("Corrupt file");
		start=coffset;
	}
	{
		int nlevels=2<<image->depth;//chroma-inflated
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=9;
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
	}
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(int[4*OCH_COUNT*2])*(BLOCKSIZE+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=statssize*sizeof(int);
		arg->hist=(int*)malloc(arg->histsize);

		arg->statssize=statssize;
		arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist||!arg->stats)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		memusage+=arg->statssize;
		
		arg->tlevels=tlevels;
		arg->clevels=clevels;
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
	for(int kt=0;kt<nblocks;kt+=nthreads)
	{
		int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
		for(int kt2=0;kt2<nthreads2;++kt2)
		{
			ThreadArgs *arg=args+kt2;
			int kx, ky;

			kx=kt+kt2;
			ky=kx/xblocks;
			kx%=xblocks;
			arg->x1=BLOCKSIZE*kx;
			arg->y1=BLOCKSIZE*ky;
			arg->x2=MINVAR(arg->x1+BLOCKSIZE, image->iw);
			arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			if(!fwd)
			{
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->decstart=cbuf+start;
				start+=size;
				arg->decend=cbuf+start;
			}
		}
#ifdef DISABLE_MT
		for(int k=0;k<nthreads2;++k)
			block_thread(args+k);
#else
		void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
		mt_finish(ctx);
#endif
		if(fwd)
		{
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				if(loud)
				{
					int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*image->nch*image->depth+7)>>3;
					int kx, ky;

					kx=kt+kt2;
					ky=kx/xblocks;
					kx%=xblocks;
					if(nblocks<2000)
					{
						//if(!(kt+kt2))
						//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
						printf(
							"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes  %10.6lf%%  CR %10lf\n",
							kt+kt2+1, nblocks,
							kx, ky,
							arg->y2-arg->y1,
							arg->x2-arg->x1,
							blocksize,
							arg->bestsize,
							arg->list.nobj,
							100.*arg->list.nobj/blocksize,
							(double)blocksize/arg->list.nobj
						);
					}
					esize+=arg->bestsize;
				}
				memcpy(data[0]->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
				dlist_appendtoarray(&arg->list, data);
				dlist_clear(&arg->list);
			}
		}
	}
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t csize=data[0]->count-start;
			printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
			printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		free(arg->stats);
	}
	free(args);
	return 0;
}