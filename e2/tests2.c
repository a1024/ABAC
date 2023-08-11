#include"e2.h"
#include"ac.h"
#include"ac_simd.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

static int clamp4(int x, int a, int b, int c, int d)
{
	int vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmax<b)vmax=b;
	if(vmin>c)vmin=c;
	if(vmax<c)vmax=c;
	if(vmin>d)vmin=d;
	if(vmax<d)vmax=d;
	x=CLAMP(vmin, x, vmax);
	return x;
}


//T40: random generated predictors
#if 0
	#define T40_PRETRAINED
	#define T40_PRINTPARAMS
	#define T40_UNSIGNED_BITS
//	#define T40_APPLY_SPATIAL
//	#define T40_DISABLE_COUNTER
//	#define T40_DISABLE_REC
//	#define T40_DISABLE_RCT
	#define T40_PRINT_ESTIMATOR_CR
//	#define T40_DISABLE_LEARNING
	#define T40_DISABLE_ADAI
//	#define T40_ADAI_PD

#define T40_NPRED 14
//#define T40_NRANDPRED 1
#define T40_NITER 32
#define T40_NITER2 32
#define T40_LR (int)(0.101*0x10000+0.5)
//#define T40_LEARNINGITER 16

#ifdef T40_NRANDPRED
#define T40_NMAPS (T40_NPRED+T40_NRANDPRED)
#else
#define T40_NMAPS T40_NPRED
#endif
#ifndef T40_DISABLE_REC
#define T40_N_REC_ESTIMATORS 6		//15
#define T40_NESTIMATORS ((T40_N_REC_ESTIMATORS+1)*T40_NMAPS)
#else
#define T40_NESTIMATORS T40_NMAPS
#endif
typedef struct T40NodeStruct
{
	int n[2];
#ifndef T40_DISABLE_REC
	unsigned short rec[T40_N_REC_ESTIMATORS];
#endif
} T40Node;
typedef struct T40CtxStruct
{
#if T40_NRANDPRED
	int pred[T40_NRANDPRED];//original predictions
#endif
	int context[T40_NMAPS];
	ArrayHandle maps[24][T40_NMAPS];
	T40Node *node[T40_NMAPS];

#ifndef T40_DISABLE_ADAI
#ifdef T40_ADAI_PD
	double
#else
	int
#endif
		lr,
		beta[3], eps,//beta[1] is beta[2]^iter
		grad[T40_NESTIMATORS],
		velocity[T40_NESTIMATORS], v_hat[T40_NESTIMATORS], beta1[T40_NESTIMATORS], momentum[T40_NESTIMATORS], beta_acc[T40_NESTIMATORS];//Adai optimizer
#else
	int lr;
#endif

	int p0arr[T40_NESTIMATORS], p0_0, p0;//p0_0 isn't clamped
	int weights[24][T40_NESTIMATORS];
	long long wsum;

	int nnodes;
#ifdef T40_PRINT_ESTIMATOR_CR
	float csizes_est[24*T40_NESTIMATORS];
#endif
} T40Ctx;
T40Ctx* t40_ctx_init()
{
	int val=0x8000;
	T40Node node0={{1, 1}};
#ifndef T40_DISABLE_REC
	for(int k=0;k<T40_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	T40Ctx *ctx=(T40Ctx*)malloc(sizeof(T40Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T40Ctx));
	
#ifndef T40_DISABLE_ADAI
#ifdef T40_ADAI_PD
	ctx->lr=0.01;//Adai optimizer
	ctx->beta[0]=0.1;
	ctx->beta[1]=1;
	ctx->beta[2]=0.99;
	ctx->eps=0.001;
	for(int k=0;k<T40_NESTIMATORS;++k)
		ctx->beta_acc[k]=1;
#else
	ctx->lr=T40_LR;//Adai optimizer
	ctx->beta[0]=(int)(0.1*0x10000+0.5);
	ctx->beta[1]=0x10000;
	ctx->beta[2]=(int)(0.99*0x10000+0.5);
	ctx->eps=(int)(0.001*0x10000+0.5);
	for(int k=0;k<T40_NESTIMATORS;++k)
		ctx->beta_acc[k]=0x10000;
#endif
#else
	ctx->lr=T40_LR;
#endif

	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T40_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T40Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T40Node), sizeof(T40Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t40_ctx_clear(T40Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T40_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
void t40_ctx_get_context(T40Ctx *ctx, const char *pixels, const char **errors, const short *params, int iw, int ih, int kc, int kx, int ky)
{
	char nb[CUSTOM_NPARAMS]={0};
	int idx=0;
	for(int ky2=-CUSTOM_REACH;ky2<0;++ky2)
	{
		for(int kx2=-CUSTOM_REACH;kx2<=CUSTOM_REACH;++kx2, ++idx)
		{
			if((unsigned)(ky+ky2)<(unsigned)ih&&(unsigned)(kx+kx2)<(unsigned)iw)
				nb[idx]=pixels[(iw*(ky+ky2)+kx+kx2)<<2|kc];
		}
	}
	for(int kx2=-CUSTOM_REACH;kx2<0;++kx2, ++idx)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
			nb[idx]=pixels[(iw*ky+kx+kx2)<<2|kc];
	}
	
	int j=-1;
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	char
#if CUSTOM_REACH==2
		NNWW =nb[ 0],
		NNW  =nb[ 1],
		NN   =nb[ 2],
		NNE  =nb[ 3],
		NNEE =nb[ 4],
		NWW  =nb[ 5],
		NW   =nb[ 6],
		N    =nb[ 7],
		NE   =nb[ 8],
		NEE  =nb[ 9],
		WW   =nb[10],
		W    =nb[11],
#else
		NNWW =LOAD(pixels, 0,  2, 2),
		NNW  =LOAD(pixels, 0,  1, 2),
		NN   =LOAD(pixels, 0,  0, 2),
		NNE  =LOAD(pixels, 0, -1, 2),
		NNEE =LOAD(pixels, 0, -2, 2),
		NWW  =LOAD(pixels, 0,  2, 1),
		NW   =LOAD(pixels, 0,  1, 1),
		N    =LOAD(pixels, 0,  0, 1),
		NE   =LOAD(pixels, 0, -1, 1),
		NEE  =LOAD(pixels, 0, -2, 1),
		WW   =LOAD(pixels, 0,  2, 0),
		W    =LOAD(pixels, 0,  1, 0),
#endif

		m1  =LOAD(pixels, 1, 0, 0),
		Nm1 =LOAD(pixels, 1, 0, 1),
		Wm1 =LOAD(pixels, 1, 1, 0),
		NWm1=LOAD(pixels, 1, 1, 1),

		m2  =LOAD(pixels, 2, 0, 0),
		Nm2 =LOAD(pixels, 2, 0, 1),
		Wm2 =LOAD(pixels, 2, 1, 0),
		NWm2=LOAD(pixels, 2, 1, 1);
#undef LOAD
	ctx->context[++j]=0;//0
	ctx->context[++j]=N;//1
	ctx->context[++j]=W;//2
	ctx->context[++j]=NW;//3
	ctx->context[++j]=m1;//4
	ctx->context[++j]=W+NE-N;//5
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;//6
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);//7		best
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);//8
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);//9
	ctx->context[++j]=NW+NE-NN;//10		worst
	ctx->context[++j]=(N+W-NW + m1)>>1;//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=(N+W-NW + m2)>>1;//13
	
#if T40_NRANDPRED
	for(int kp=0;kp<T40_NRANDPRED;++kp)
	{
		idx=0;
		for(int ky2=-CUSTOM_REACH_E;ky2<0;++ky2)
		{
			for(int kx2=-CUSTOM_REACH_E;kx2<=CUSTOM_REACH_E;++kx2, ++idx)
			{
				if((unsigned)(ky+ky2)<(unsigned)ih&&(unsigned)(kx+kx2)<(unsigned)iw)
					nb[CUSTOM_REACH+idx]=errors[kp][(iw*(ky+ky2)+kx+kx2)<<2|kc];
			}
		}
		for(int kx2=-CUSTOM_REACH_E;kx2<0;++kx2, ++idx)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw)
				nb[CUSTOM_REACH+idx]=errors[kp][(iw*ky+kx+kx2)<<2|kc];
		}

		const short *p=params+CUSTOM_NPARAMS*kp;
		int pred=0;
		for(int k2=0;k2<CUSTOM_NPARAMS;++k2)
			pred+=nb[k2]*p[k2];
		//for(int k2=0;k2<CUSTOM_NNB;++k2)
		//	pred+=nb[k2]*p[k2];
		//for(int k2=0;k2<CUSTOM_NNB_E;++k2)
		//	pred+=nb[CUSTOM_REACH+k2]*p[CUSTOM_NNB+k2];
#ifdef CUSTOM_USE_MULHRS
		pred>>=14;//_mm256_mulhrs_epi16
#else
		pred>>=12;//_mm256_mulhi_epi16
#endif
		pred=CLAMP(-128, pred, 127);
		ctx->pred[kp]=pred;
		ctx->context[++j]=ctx->pred[kp];
	}
#endif
	++j;
	for(int k=0;k<j;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
void t40_ctx_estimate_p0(T40Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	T40Node *node;
	for(int kp=0;kp<T40_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		node=ctx->node[kp]=(T40Node*)array_at(ctx->maps[workidx]+kp, ctx->context[kp]);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T40_DISABLE_REC
		for(;k2<T40_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T40_NESTIMATORS;++k)
	{
#ifdef T40_DISABLE_COUNTER
		if(k%(T40_N_REC_ESTIMATORS+1))//
#endif
		{
			sum+=(long long)ctx->p0arr[k]*wk[k];
			ctx->wsum+=wk[k];
		}
	}
	//ctx->p0=ctx->wsum?(int)((sum+(ctx->wsum>>1))/ctx->wsum):0x8000;//same CR
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t40_ctx_update(T40Ctx *ctx, int kc, int kb, int bit, int dummy)
{
	int workidx=kc<<3|kb;

#ifdef T40_PRINT_ESTIMATOR_CR
	for(int k=0;k<T40_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T40_NESTIMATORS*workidx+k]+=bitsize;
		}
	}
#endif
	//bwd
#ifndef T40_DISABLE_LEARNING
	int *wk=ctx->weights[workidx];
	if(ctx->p0_0>=1&&ctx->p0_0<=0xFFFF)
	{
		int p_bit=bit?0x10000-ctx->p0:ctx->p0;
		long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
		dL_dp0^=-bit;
		dL_dp0+=bit;
		
#ifndef T40_DISABLE_ADAI
#ifdef T40_ADAI_PD
		double v_sum=0;//Adai optimizer
		ctx->beta[1]*=ctx->beta[2];
		for(int k=0;k<T40_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			ctx->grad[k] = (double)dL_dp0*diff/(ctx->wsum<<16);

			double g2=ctx->grad[k]*ctx->grad[k];
			ctx->velocity[k]=g2+(ctx->velocity[k]-g2)*ctx->beta[2];
			ctx->v_hat[k]=ctx->velocity[k]/(1-ctx->beta[1]);
			v_sum+=ctx->v_hat[k];
		}
		v_sum/=T40_NESTIMATORS;
		for(int k=0;k<T40_NESTIMATORS;++k)
		{
			double vmax=1-ctx->eps;
			if(v_sum)
			{
				ctx->beta1[k]=1-ctx->beta[0]*ctx->v_hat[k]/v_sum;
				ctx->beta1[k]=CLAMP(0, ctx->beta1[k], vmax);
			}
			else
				ctx->beta1[k]=vmax;
			ctx->momentum[k]=ctx->p0arr[k]+(ctx->momentum[k]-ctx->p0arr[k])*ctx->beta1[k];
			ctx->beta_acc[k]*=ctx->beta1[k];
			double m_hat=ctx->momentum[k]/(1-ctx->beta_acc[k]);
			double delta=ctx->lr*m_hat;
			wk[k]+=(int)round(delta*0x10000);
			wk[k]=CLAMP(1, wk[k], 0xFFFF);
		}
#else
		long long v_sum=0;//Adai optimizer
		ctx->beta[1]=(int)((long long)ctx->beta[1]*ctx->beta[2]>>16);
		for(int k=0;k<T40_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			grad=CLAMP(-0x7FFFFFFF, grad, 0x7FFFFFFF);
			ctx->p0arr[k]=(int)grad;

			int g2=(int)(grad*grad>>16);
			ctx->velocity[k]=g2+(int)((long long)(ctx->velocity[k]-g2)*ctx->beta[2]>>16);
			ctx->v_hat[k]=((long long)ctx->velocity[k]<<16)/(0x10000-ctx->beta[1]);
			v_sum+=ctx->v_hat[k];
		}
		v_sum/=T40_NESTIMATORS;
		for(int k=0;k<T40_NESTIMATORS;++k)
		{
			int vmax=0x10000-ctx->eps;
			if(v_sum)
			{
				ctx->beta1[k]=0x10000-(int)(ctx->beta[0]*ctx->v_hat[k]/v_sum);
				ctx->beta1[k]=CLAMP(0, ctx->beta1[k], vmax);
			}
			else
				ctx->beta1[k]=vmax;
			ctx->momentum[k]=ctx->p0arr[k]+(int)((long long)(ctx->momentum[k]-ctx->p0arr[k])*ctx->beta1[k]>>16);
			ctx->beta_acc[k]=(int)((long long)ctx->beta_acc[k]*ctx->beta1[k]>>16);
			int m_hat=(int)(((long long)ctx->momentum[k]<<16)/(0x10000-ctx->beta_acc[k]));
			int delta=ctx->lr*m_hat>>16;
			wk[k]+=delta;
			wk[k]=CLAMP(1, wk[k], 0xFFFF);
		}
#endif
#else

		//long long delta=0;
		for(int k=0;k<T40_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=ctx->lr*grad>>16;
			//delta+=wnew*wnew;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
#endif
	}
#endif

	if(dummy)
		return;

	//update
#ifndef T40_DISABLE_REC
	static const int shifts[]=
	{
		//7, 7, 7, 6, 5, 0, 0, 8,//T40_N_REC_ESTIMATORS 1
		//8, 7, 7, 6, 6, 6, 7, 7,
		//7, 7, 7, 5, 4, 4, 2, 7,

		5, 5, 4, 2, 1, 0, 0, 5,//T40_N_REC_ESTIMATORS 6
		7, 5, 5, 4, 4, 5, 4, 6,
		5, 5, 5, 2, 1, 1, 1, 5,

		//6, 5, 5, 4, 2, 0, 0, 6,
		//9, 7, 6, 5, 4, 5, 4, 6,
		//6, 5, 5, 4, 4, 3, 1, 6,
	};
	//static const int bitrating[]=//higher is more compressible
	//{
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//	1, 1, 1, 1, 2, 3, 4, 1,
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//};
#endif
	T40Node *node;
	for(int kp=0;kp<T40_NMAPS;++kp)
	{
		node=ctx->node[kp];
		++node->n[bit];
#ifndef T40_DISABLE_REC
		//static const int lgdens[]={0, 1, 2, 4, 8, 14};
		for(int k=0;k<T40_N_REC_ESTIMATORS;++k)
		{
			//int lgden=k<<1;
			//int lgden=k*5>>1;//x2.5
			int lgden=k;
			//int lgden=k*3;//X
			//int lgden=k<<2;//X
			//int lgden=k+1;
			//int lgden=k+2;
			//int lgden=k+shifts[workidx];
			//int lgden=(k+1)<<1;//X
			//int lgden=k+((4-bitrating[workidx])<<1);
			//int lgden=k+4-bitrating[workidx];
			//int lgden=k+3;
			//int lgden=lgdens[k];
			//int lgden=((k+1)<<1)-1;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);
	}
}
#if T40_NRANDPRED
short t40_params[3][T40_NRANDPRED][CUSTOM_NPARAMS];//TODO insert params in codestream
const short t40_params_kodim13[]=
{
	//CLIC11-crop4-2
#if 1
	 0x0020,-0x0045, 0x0005, 0x007C, 0x0006,
	 0x00C1, 0x0181,-0x0001, 0x0124, 0x0169,
	 0x006A, 0x026B,

	 0x0064, 0x014F, 0x0015, 0x0161, 0x003A,
	-0x0002, 0x0041, 0x0106,-0x000D, 0x0044,
	 0x00B0, 0x008C,

	 0x0083, 0x0021,-0x0062,-0x004D,-0x0007,
	 0x00EF, 0x01C7, 0x0118, 0x01BD, 0x007B,
	 0x0109, 0x00FD,

	 0x0048, 0x0090,-0x0020, 0x00E2,-0x0069,
	 0x0105, 0x0026, 0x01F3, 0x0154, 0x00A1,
	-0x0041, 0x0171,

	 0x015C,-0x0047, 0x00D1, 0x00E0, 0x005E,
	 0x0009, 0x004C, 0x00B1, 0x002B, 0x0101,
	 0x0254, 0x006C,

	-0x001C, 0x0117,-0x0066, 0x0075, 0x00EA,
	-0x005C,-0x0008, 0x01C5, 0x01A9, 0x0122,
	 0x0006, 0x0002,
#endif

	//kodim13
#if 0
	-0x0034,-0x001B,-0x000E,-0x0001, 0x000A,
	 0x00B4, 0x00A9,-0x0013, 0x0227, 0x0109,
	 0x0005, 0x0398,

	-0x0070, 0x00F7, 0x0033, 0x00A5, 0x00A6,
	-0x0085,-0x0015, 0x0166, 0x0044, 0x003B,
	 0x0117, 0x0082,

	 0x0088, 0x0009,-0x0088, 0x0031, 0x0087,
	-0x0024, 0x013F, 0x0154, 0x00BB, 0x0025,
	 0x0117, 0x02D5,

	-0x006B,-0x0064,-0x0058,-0x0043,-0x0049,
	-0x009B,-0x0009, 0x0252, 0x0072, 0x006C,
	-0x00B6, 0x0322,

	 0x0052, 0x0002, 0x000C, 0x00BE,-0x001B,
	 0x00A5, 0x00C5, 0x013D, 0x004E, 0x00FD,
	 0x01E0, 0x013D,

	 0x0052, 0x0105, 0x0063,-0x0055, 0x0056,
	-0x0050, 0x0108, 0x018D, 0x013D, 0x0099,
	 0x0082, 0x003B,
#endif
};
#endif
static void print_causalkernel(const short *params, int reach)
{
	int width=reach<<1|1, count=reach*(reach+1)*2;
	for(int ky=-reach, idx=0;ky<=0&&idx<count;++ky)
	{
		for(int kx=-reach;kx<=reach&&idx<count;++kx, ++idx)
		{
			short val=params[idx];
			printf("%c0x%04X,", val<0?'-':' ', abs(val));
			//printf("%6d,", val);
		}
		printf("\n");
	}
	printf("\n");
}
int t40_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T40 Enc  Random generated predictors  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifdef T40_APPLY_SPATIAL
	apply_transforms_fwd(buf2, iw, ih);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//buffer is signed
#else
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
#ifndef T40_DISABLE_RCT
	colortransform_ycocb_fwd(buf2, iw, ih);
#endif
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//X  the buffer is signed
#endif

#if T40_NRANDPRED
#ifdef ALLOW_SRAND
	srand((unsigned)__rdtsc());
#endif
#ifdef T40_PRETRAINED
	memcpy(t40_params, t40_params_kodim13, sizeof(t40_params));
#else
	for(int kc=0;kc<3;++kc)
	{
		for(int k=0;k<T40_NRANDPRED;++k)
		{
			for(int k2=0;k2<CUSTOM_NPARAMS;++k2)
				t40_params[kc][k][k2]=(rand()&0x1FF)-0x100;

			//double best=0;
			float curr=0;
			for(int k2=0;k2<T40_NITER;++k2)
			{
				curr=opt_custom_v2(buf2, iw, ih, kc, T40_NITER2, t40_params[kc][k], curr, loud==2);//120 optimization sequences (encode only)
				//curr=opt_custom(buf2, iw, ih, kc, T40_NITER2, t40_params[kc][k], loud==2);
				if(loud)
				{
					TimeInfo ti;
					parsetimedelta(time_ms()-t_start, &ti);
					printf("[%d/%d  %2d/%2d  %3d/%3d]  %6.2lf%%  CR%11lf  %02d-%02d-%06.3f%c", kc+1, 3, k+1, T40_NRANDPRED, k2+1, T40_NITER, 100.*(T40_NITER*(T40_NRANDPRED*kc+k)+k2+1)/(3*T40_NRANDPRED*T40_NITER), 1/curr, ti.hours, ti.mins, ti.secs, loud==2?'\n':'\r');
				}
				//if(!best)
				//	best=curr;
				//else if(best<curr)
				//	LOG_ERROR("CR decreased from %lf to %lf", 1/best, 1/curr);
				//else
				//	best=curr;
			}
			if(loud)
				printf("\n");
		}
	}
#ifdef T40_PRINTPARAMS
	for(int k=0;k<3*T40_NRANDPRED*CUSTOM_NPARAMS;)
	{
		print_causalkernel((short*)t40_params+k, CUSTOM_REACH);
		k+=CUSTOM_NNB;
		print_causalkernel((short*)t40_params+k, CUSTOM_REACH_E);
		k+=CUSTOM_NNB_E;
	}
#endif
#endif

	char *ebuf=(char*)malloc((size_t)T40_NRANDPRED*res<<2);
	if(!ebuf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ebuf, 0, (size_t)T40_NRANDPRED*res<<2);
	char *ebufs[T40_NRANDPRED];
	for(int k=0;k<T40_NRANDPRED;++k)
		ebufs[k]=ebuf+(res<<2)*k;
#endif

	T40Ctx *t40_ctx=t40_ctx_init();
	if(!t40_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
#if T40_NRANDPRED
				t40_ctx_get_context(t40_ctx, buf2, ebufs, (short*)t40_params[kc], iw, ih, kc, kx, ky);
#else
				t40_ctx_get_context(t40_ctx, buf2, 0, 0, iw, ih, kc, kx, ky);
#endif
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t40_ctx_estimate_p0(t40_ctx, kc, kb);
#ifdef T40_UNSIGNED_BITS
					int bit=(buf2[idx]+128)>>kb&1;
#else
					int bit=buf2[idx]>>kb&1;
#endif
					abac_enc(&ctx, t40_ctx->p0, bit);
					
					int prob=bit?0x10000-t40_ctx->p0:t40_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//
					
#ifdef T40_LEARNINGITER
					for(int k=0;k<T40_LEARNINGITER;++k)
					{
						t40_ctx_update(t40_ctx, kc, kb, bit, 1);
						t40_ctx_estimate_p0(t40_ctx, kc, kb);
					}
#endif
					t40_ctx_update(t40_ctx, kc, kb, bit, 0);
				}
#if T40_NRANDPRED
				for(int k=0;k<T40_NRANDPRED;++k)
					ebufs[k][idx]=buf2[idx]-t40_ctx->pred[k];
#endif
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f%c", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev), loud==2?'\n':'\r');
			csize_prev=csize;
		}
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t40_ctx->nnodes*sizeof(T40Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
		float chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=7;kb>=0;--kb)
		{
			printf("B%d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc<<3|kb;
				float size=csizes[idx];
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);

#ifdef T40_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t40_ctx->csizes_est+T40_NESTIMATORS*kb;
				for(int ke=1;ke<T40_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T40_NESTIMATORS;++ke)
			{
				float *sizes=t40_ctx->csizes_est+ke;
#ifndef T40_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T40_N_REC_ESTIMATORS+1), ke%(T40_N_REC_ESTIMATORS+1));
#else
				printf("E%3d ", ke);
#endif
				for(int kb=0;kb<24;++kb)
				{
					char c;
					if(ke==minidx[kb])
						c='*';
					else if(ke==maxidx[kb])
						c='L';
					else
						c=' ';
					printf("%8.2f %c", iw*ih/sizes[T40_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T40_NESTIMATORS*kb]/t40_ctx->csizes_est[T40_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T40_DISABLE_REC
				if(!((ke+1)%(T40_N_REC_ESTIMATORS+1)))
#else
				if(!((ke+1)%8))
#endif
				{
					printf("\n");
					printf("\t\t*         **        ***       ****      ****      ****      ****      *             *         *         *         *         **        ***       ****      *             *         **        ***       ****      ****      ****      ****      *\n");
					printf("\t\t0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7\n");
				}
			}
		}
#endif
	}
	t40_ctx_clear(&t40_ctx);
	dlist_clear(&list);
	free(buf2);
#if T40_NRANDPRED
	free(ebuf);
#endif
	return 1;
}
int t40_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T40 Dec  %s\n", g_buf);
	}
	
#if T40_NRANDPRED
	char *ebuf=(char*)malloc((size_t)T40_NRANDPRED*res<<2);
	if(!ebuf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ebuf, 0, (size_t)T40_NRANDPRED*res<<2);
	char *ebufs[T40_NRANDPRED];
	for(int k=0;k<T40_NRANDPRED;++k)
		ebufs[k]=ebuf+(res<<2)*k;
#endif

	T40Ctx *t40_ctx=t40_ctx_init();
	if(!t40_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t40_ctx_init(t40_ctx);
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
#if T40_NRANDPRED
				t40_ctx_get_context(t40_ctx, (char*)buf, ebufs, (short*)t40_params[kc], iw, ih, kc, kx, ky);
#else
				t40_ctx_get_context(t40_ctx, (char*)buf, 0, 0, iw, ih, kc, kx, ky);
#endif
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t40_ctx_estimate_p0(t40_ctx, kc, kb);
					
					int bit=abac_dec(&ctx, t40_ctx->p0);
					buf[idx]|=bit<<kb;

#ifdef T40_LEARNINGITER
					for(int k=0;k<T40_LEARNINGITER;++k)
					{
						t40_ctx_update(t40_ctx, kc, kb, bit, 1);
						t40_ctx_estimate_p0(t40_ctx, kc, kb);
					}
#endif
					t40_ctx_update(t40_ctx, kc, kb, bit, 0);
				}
#ifdef T40_UNSIGNED_BITS
				buf[idx]+=128;//unsigned -> signed
#endif
#if T40_NRANDPRED
				for(int k=0;k<T40_NRANDPRED;++k)
					ebufs[k][idx]=buf[idx]-t40_ctx->pred[k];
#endif
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
	}
	
#ifdef T40_APPLY_SPATIAL
	addbuf(buf, iw, ih, 3, 4, 128);//buffer is signed
	apply_transforms_inv(buf, iw, ih);
#else
	//addbuf(buf, iw, ih, 3, 4, 128);//X  the buffer is signed
#ifndef T40_DISABLE_RCT
	colortransform_ycocb_inv((char*)buf, iw, ih);
#endif
	addbuf(buf, iw, ih, 3, 4, 128);
#endif
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	t40_ctx_clear(&t40_ctx);
#if T40_NRANDPRED
	free(ebuf);
#endif
	return 1;
}
#endif


//T41 SIMD ABAC
#if 0
//	#define T41_DISABLE_COUNTER
//	#define T41_DISABLE_REC
//	#define T41_PRINT_ESTIMATOR_CR
//	#define T41_DISABLE_LEARNING

#define T41_NPRED 14
#define T41_LR (int)(0.001*0x100+0.5)
#define T41_NMAPS T41_NPRED
#ifndef T41_DISABLE_REC
#define T41_N_REC_ESTIMATORS 4		//8
#define T41_NESTIMATORS ((T41_N_REC_ESTIMATORS+1)*T41_NMAPS)
#else
#define T41_NESTIMATORS T41_NMAPS
#endif
typedef struct T41NodeStruct
{
	int n[2];
#ifndef T41_DISABLE_REC
	unsigned char rec[T41_N_REC_ESTIMATORS];
#endif
} T41Node;
typedef struct T41CtxStruct
{
	__m256i context[T41_NMAPS];
	ArrayHandle maps[24][T41_NMAPS];
	T41Node *node[T41_NMAPS][16];

	short lr;
	short p0arr[16][T41_NESTIMATORS], p0_0[16];//p0_0 isn't clamped
	__m256i p0;
	short weights[24][16][T41_NESTIMATORS];
	int wsums[16];

	int nnodes;
#ifdef T41_PRINT_ESTIMATOR_CR
	float csizes_est[24*T41_NESTIMATORS];
#endif
} T41Ctx;
T41Ctx* t41_ctx_init()
{
	short val=0x80;
	T41Node node0={{1, 1}};
#ifndef T41_DISABLE_REC
	for(int k=0;k<T41_N_REC_ESTIMATORS;++k)
		node0.rec[k]=(unsigned char)val;
#endif
	T41Ctx *ctx=(T41Ctx*)malloc(sizeof(T41Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T41Ctx));
	ctx->lr=T41_LR;

	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(val));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T41_NMAPS;++k2)
		{
			int nnodes=16*256<<(7-kb);
			ARRAY_ALLOC(T41Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, nnodes*sizeof(T41Node), sizeof(T41Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t41_ctx_clear(T41Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T41_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
void t41_ctx_get_context(T41Ctx *ctx, const char *pixels, int iw, int ih, int kc, int kx, int ky, int x1, int x2, int y1, int y2)
{
	int j=-1;
	int dx=x2-x1, dy=y2-y1;
	//_mm256_set_epi16 has MSB in first arg
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X)-x1)<(unsigned)dx&&(unsigned)(ky-Y-y1)<(unsigned)dy?_mm256_set_epi16(\
	BUF[(iw*(ky-Y+dy*3)+kx-(X)+dx*3)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*3)+kx-(X)+dx*2)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*3)+kx-(X)+dx  )<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*3)+kx-(X)     )<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*2)+kx-(X)+dx*3)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*2)+kx-(X)+dx*2)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*2)+kx-(X)+dx  )<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy*2)+kx-(X)     )<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy  )+kx-(X)+dx*3)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy  )+kx-(X)+dx*2)<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy  )+kx-(X)+dx  )<<2|(kc-C)],\
	BUF[(iw*(ky-Y+dy  )+kx-(X)     )<<2|(kc-C)],\
	BUF[(iw*(ky-Y     )+kx-(X)+dx*3)<<2|(kc-C)],\
	BUF[(iw*(ky-Y     )+kx-(X)+dx*2)<<2|(kc-C)],\
	BUF[(iw*(ky-Y     )+kx-(X)+dx  )<<2|(kc-C)],\
	BUF[(iw*(ky-Y     )+kx-(X)     )<<2|(kc-C)]\
):_mm256_setzero_si256()
	//int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	__m256i
	//	NNWW =LOAD(pixels, 0,  2, 2),
	//	NNW  =LOAD(pixels, 0,  1, 2),
		NN   =LOAD(pixels, 0,  0, 2),
	//	NNE  =LOAD(pixels, 0, -1, 2),
	//	NNEE =LOAD(pixels, 0, -2, 2),
	//	NWW  =LOAD(pixels, 0,  2, 1),
		NW   =LOAD(pixels, 0,  1, 1),
		N    =LOAD(pixels, 0,  0, 1),
		NE   =LOAD(pixels, 0, -1, 1),
	//	NEE  =LOAD(pixels, 0, -2, 1),
	//	WW   =LOAD(pixels, 0,  2, 0),
		W    =LOAD(pixels, 0,  1, 0),

		m1  =LOAD(pixels, 1, 0, 0),
		Nm1 =LOAD(pixels, 1, 0, 1),
		Wm1 =LOAD(pixels, 1, 1, 0),
	//	NWm1=LOAD(pixels, 1, 1, 1),

		m2  =LOAD(pixels, 2, 0, 0);
	//	Nm2 =LOAD(pixels, 2, 0, 1),
	//	Wm2 =LOAD(pixels, 2, 1, 0),
	//	NWm2=LOAD(pixels, 2, 1, 1);
#undef LOAD
	ctx->context[++j]=_mm256_setzero_si256();//0
	ctx->context[++j]=N;//1
	ctx->context[++j]=W;//2
	ctx->context[++j]=NW;//3
	ctx->context[++j]=m1;//4
	ctx->context[++j]=_mm256_sub_epi16(_mm256_add_epi16(W, NE), N);//5
	ctx->context[++j]=_mm256_add_epi16(_mm256_add_epi16(W, N), m1);//6
	ctx->context[++j]=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);//7		best
	ctx->context[++j]=_mm256_sub_epi16(_mm256_add_epi16(N, m1), Nm1);//8
	ctx->context[++j]=_mm256_sub_epi16(_mm256_add_epi16(W, m1), Wm1);//9
	ctx->context[++j]=_mm256_sub_epi16(_mm256_add_epi16(NW, NE), NN);//10		worst
	ctx->context[++j]=_mm256_srli_epi16(_mm256_add_epi16(_mm256_add_epi16(N, W), _mm256_sub_epi16(m1, NW)), 1);//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=_mm256_srli_epi16(_mm256_add_epi16(_mm256_add_epi16(N, W), _mm256_sub_epi16(m2, NW)), 1);//13
	
	++j;
	__m256i half=_mm256_set1_epi16(128), full=_mm256_set1_epi16(255);
	for(int k=0;k<j;++k)
	{
		ctx->context[k]=_mm256_add_epi16(ctx->context[k], half);
		ctx->context[k]=_mm256_max_epi16(ctx->context[k], _mm256_setzero_si256());
		ctx->context[k]=_mm256_min_epi16(ctx->context[k], full);
	}
}
void t41_ctx_estimate_p0(T41Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;

	int p0idx;
	int sum;
	T41Node *node;
	for(int kl=0;kl<16;++kl)
	{
		p0idx=0;
		for(int kp=0;kp<T41_NMAPS;++kp)//for each predictor
		{
			node=ctx->node[kp][kl]=(T41Node*)array_at(ctx->maps[workidx]+kp, ctx->context[kp].m256i_u16[kl]<<4|kl);
		
			sum=node->n[0]+node->n[1];
			ctx->p0arr[kl][p0idx]=sum?(int)(((long long)node->n[0]<<8)/sum):0x80;
			++p0idx;
#ifndef T41_DISABLE_REC
			for(int k2=0;k2<T41_N_REC_ESTIMATORS;++k2, ++p0idx)
				ctx->p0arr[kl][p0idx]=node->rec[k2];
#endif
		}
	}

	for(int kl=0;kl<16;++kl)
	{
		short *wk=ctx->weights[workidx][kl];
		sum=0;
		ctx->wsums[kl]=0;
		for(int kp=0;kp<T41_NESTIMATORS;++kp)
		{
#ifdef T41_DISABLE_COUNTER
			if(kp%(T41_N_REC_ESTIMATORS+1))//
#endif
			{
				sum+=ctx->p0arr[kl][kp]*wk[kp];
				ctx->wsums[kl]+=wk[kp];
			}
		}
		ctx->p0_0[kl]=ctx->wsums[kl]?(int)(sum/ctx->wsums[kl]):0x80;

		ctx->p0.m256i_u16[kl]=CLAMP(1, ctx->p0_0[kl], 0xFF);
	}
}
void t41_ctx_update(T41Ctx *ctx, int kc, int kb, __m256i const *bits)
{
	int workidx=kc<<3|kb;

#ifdef T41_PRINT_ESTIMATOR_CR
	for(int kp=0;kp<T41_NESTIMATORS;++kp)
	{
		int prob=(bit?0x100-ctx->p0arr[kp]:ctx->p0arr[kp]);
		if(prob)
		{
			float p=(float)prob/0x100;
			float bitsize=-log2f(p);
			ctx->csizes_est[T41_NESTIMATORS*workidx+kp]+=bitsize;
		}
	}
#endif
	//bwd
#ifndef T41_DISABLE_LEARNING
	for(int kl=0;kl<16;++kl)
	{
		short *wk=ctx->weights[workidx][kl];
		if(ctx->p0_0[kl]>=1&&ctx->p0_0[kl]<=0xFFFF)
		{
			int bit=bits->m256i_u16[kl]&1;
			int p_bit=bit?0x100-ctx->p0.m256i_u16[kl]:ctx->p0.m256i_u16[kl];
			long long dL_dp0=-(1LL<<32)/p_bit;//fixed 7.8 bit
			dL_dp0^=-bit;
			dL_dp0+=bit;

			for(int kp=0;kp<T41_NESTIMATORS;++kp)
			{
				int diff=ctx->p0arr[kl][kp]-ctx->p0.m256i_u16[kl];//fixed 7.9 bit
				long long grad = dL_dp0*diff/ctx->wsums[kl];
				long long wnew=ctx->lr*grad>>8;
				wnew=wk[kp]-wnew;
				wnew=CLAMP(1, wnew, 0xFF);
				wk[kp]=(int)wnew;
			}
		}
	}
#endif

	//update
	T41Node *node;
	//if(kb==7)//
	//	printf("");
	for(int kp=0;kp<T41_NMAPS;++kp)
	{
		for(int kl=0;kl<16;++kl)
		{
			int bit=bits->m256i_u16[kl];
			node=ctx->node[kp][kl];
			++node->n[bit];
#ifndef T41_DISABLE_REC
			for(int kp=0;kp<T41_N_REC_ESTIMATORS;++kp)
			{
				//int lgden=kp<<1;
				//int lgden=kp*5>>1;//x2.5
				int lgden=kp;
				//int lgden=kp*3;//X
				//int lgden=kp<<2;//X
				//int lgden=kp+1;
				//int lgden=kp+2;
				int temp=node->rec[kp]+(((!bit<<8)-node->rec[kp])>>lgden);
				node->rec[kp]=CLAMP(1, temp, 0xFF);
			}
#endif
			//if(bit)
			//	printf("");
			//unsigned short c0=ctx->context[kp].m256i_u16[kl];
			ctx->context[kp].m256i_u16[kl]|=bit<<(8+7-kb);

			//if((ctx->context[kp>>4].m256i_u16[kp&15]<<4|15)>(16*256<<(7-kb)))
			//	LOG_ERROR("Context overflow");
		}
	}
}
int t41_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T41 Enc  SIMD ABAC  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifdef T41_APPLY_SPATIAL
	apply_transforms_fwd(buf2, iw, ih);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//buffer is signed
#else
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd(buf2, iw, ih);
#endif

#if T41_NRANDPRED
#ifdef ALLOW_SRAND
	srand((unsigned)__rdtsc());
#endif
#ifdef T41_PRETRAINED
	memcpy(t41_params, t41_params_kodim13, sizeof(t41_params));
#else
	for(int kc=0;kc<3;++kc)
	{
		for(int kp=0;kp<T41_NRANDPRED;++kp)
		{
			for(int k2=0;k2<CUSTOM_NPARAMS;++k2)
				t41_params[kc][kp][k2]=(rand()&0x1FF)-0x100;

			//double best=0;
			float curr=0;
			for(int k2=0;k2<T41_NITER;++k2)
			{
				curr=opt_custom_v2(buf2, iw, ih, kc, T41_NITER2, t41_params[kc][kp], curr, loud==2);//120 optimization sequences (encode only)
				//curr=opt_custom(buf2, iw, ih, kc, T41_NITER2, t41_params[kc][kp], loud==2);
				if(loud)
				{
					TimeInfo ti;
					parsetimedelta(time_ms()-t_start, &ti);
					printf("[%d/%d  %2d/%2d  %3d/%3d]  %6.2lf%%  CR%11lf  %02d-%02d-%06.3f%c", kc+1, 3, kp+1, T41_NRANDPRED, k2+1, T41_NITER, 100.*(T41_NITER*(T41_NRANDPRED*kc+kp)+k2+1)/(3*T41_NRANDPRED*T41_NITER), 1/curr, ti.hours, ti.mins, ti.secs, loud==2?'\n':'\r');
				}
				//if(!best)
				//	best=curr;
				//else if(best<curr)
				//	LOG_ERROR("CR decreased from %lf to %lf", 1/best, 1/curr);
				//else
				//	best=curr;
			}
			if(loud)
				printf("\n");
		}
	}
#ifdef T41_PRINTPARAMS
	for(int k=0;k<3*T41_NRANDPRED*CUSTOM_NPARAMS;)
	{
		print_causalkernel((short*)t41_params+k, CUSTOM_REACH);
		k+=CUSTOM_NNB;
		print_causalkernel((short*)t41_params+k, CUSTOM_REACH_E);
		k+=CUSTOM_NNB_E;
	}
#endif
#endif

	char *ebuf=(char*)malloc((size_t)T41_NRANDPRED*res<<2);
	if(!ebuf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ebuf, 0, (size_t)T41_NRANDPRED*res<<2);
	char *ebufs[T41_NRANDPRED];
	for(int k=0;k<T41_NRANDPRED;++k)
		ebufs[k]=ebuf+(res<<2)*k;
#endif

	T41Ctx *t41_ctx=t41_ctx_init();
	if(!t41_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	DList lists[16];
	for(int kl=0;kl<16;++kl)
	{
		dlist_init(lists+kl, 1, 1024, 0);
		dlist_push_back(lists+kl, 0, 4);
	}
	
	ABACSIMDEnc ctx;
	abacsimd_enc_init(&ctx, lists);
	
	float csizes[16][24]={0};
	
	int bw=iw>>2, bh=ih>>2, count=bw*bh, idx2=0;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx, ++idx2)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx[]=
				{
					(iw* ky      +kx     )<<2|kc,
					(iw* ky      +kx+bw  )<<2|kc,
					(iw* ky      +kx+bw*2)<<2|kc,
					(iw* ky      +kx+bw*3)<<2|kc,
					(iw*(ky+bh  )+kx     )<<2|kc,
					(iw*(ky+bh  )+kx+bw  )<<2|kc,
					(iw*(ky+bh  )+kx+bw*2)<<2|kc,
					(iw*(ky+bh  )+kx+bw*3)<<2|kc,
					(iw*(ky+bh*2)+kx     )<<2|kc,
					(iw*(ky+bh*2)+kx+bw  )<<2|kc,
					(iw*(ky+bh*2)+kx+bw*2)<<2|kc,
					(iw*(ky+bh*2)+kx+bw*3)<<2|kc,
					(iw*(ky+bh*3)+kx     )<<2|kc,
					(iw*(ky+bh*3)+kx+bw  )<<2|kc,
					(iw*(ky+bh*3)+kx+bw*2)<<2|kc,
					(iw*(ky+bh*3)+kx+bw*3)<<2|kc,
				};
				t41_ctx_get_context(t41_ctx, buf2, iw, ih, kc, kx, ky, 0, bw, 0, bh);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					//if(kx==1&&ky==0&&kc==0&&kb==7)
					//	printf("");

					t41_ctx_estimate_p0(t41_ctx, kc, kb);

					//if(t41_ctx->p0.m256i_u16[1]!=0x80)//
					//	printf("");

					__m256i bits=_mm256_set_epi16(//first arg is MSB
						(buf2[idx[15]]+128)>>kb&1,
						(buf2[idx[14]]+128)>>kb&1,
						(buf2[idx[13]]+128)>>kb&1,
						(buf2[idx[12]]+128)>>kb&1,
						(buf2[idx[11]]+128)>>kb&1,
						(buf2[idx[10]]+128)>>kb&1,
						(buf2[idx[ 9]]+128)>>kb&1,
						(buf2[idx[ 8]]+128)>>kb&1,
						(buf2[idx[ 7]]+128)>>kb&1,
						(buf2[idx[ 6]]+128)>>kb&1,
						(buf2[idx[ 5]]+128)>>kb&1,
						(buf2[idx[ 4]]+128)>>kb&1,
						(buf2[idx[ 3]]+128)>>kb&1,
						(buf2[idx[ 2]]+128)>>kb&1,
						(buf2[idx[ 1]]+128)>>kb&1,
						(buf2[idx[ 0]]+128)>>kb&1
					);
					abacsimd_enc(&ctx, &t41_ctx->p0, &bits);
					
					for(int kl=0;kl<16;++kl)
					{
						int prob=bits.m256i_u16[kl]?0x100-t41_ctx->p0.m256i_u16[kl]:t41_ctx->p0.m256i_u16[kl];//
						float bitsize=-log2f((float)prob*(1.f/0x100));
						csizes[kl][kc<<3|kb]+=bitsize;//
					}
					
					t41_ctx_update(t41_ctx, kc, kb, &bits);
				}
			}
			//if(loud)
			//	printf("%6.2lf%%  %d/%d\r", 100.*(idx2+1)/count, idx2, count);
		}
		if(loud)
		{
			if(ky>=bh)
				printf("");
			static float csize_prev=0;
			float csize=0;
			if(!idx2)
				csize_prev=0;
			for(int k=0;k<24;++k)
				csize+=csizes[0][k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f%c", ky+1, bh, 100.*(ky+1)/bh, bw*(ky+1)*3/csize, bw*3/(csize-csize_prev), loud==2?'\n':'\r');
			csize_prev=csize;
		}
	}
	abacsimd_enc_flush(&ctx);

	ptrdiff_t dststart0=-1;
	for(int kl=0;kl<16;++kl)
	{
		size_t dststart=dlist_appendtoarray(lists+kl, data);
		if(dststart0==-1)
			dststart0=dststart;
		memcpy(data[0]->data+dststart, &lists[kl].nobj, 4);
	}
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t41_ctx->nnodes*sizeof(T41Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
		int res2=bw*bh;
		for(int kl=0;kl<16;++kl)
		{
			float chsizes[4]={0};
			printf("lane %d:\n", kl);
			printf("\tC0\t\tC1\t\tC2\n");
			for(int kb=7;kb>=0;--kb)
			{
				printf("B%d  ", kb);
				for(int kc=0;kc<3;++kc)
				{
					int idx=kc<<3|kb;
					float size=csizes[kl][idx];
					printf(" %15.6f", res2/size);
					chsizes[kc]+=size;
				}
				printf("\n");
			}
			printf("\n");
			chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
			printf("Total%15.6f %15.6f %15.6f %15.6f\n", res2*8/chsizes[0], res2*8/chsizes[1], res2*8/chsizes[2], res2*24/chsizes[3]);
			printf("\n");
		}
		ptrdiff_t csize=data[0]->count-dststart0;
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)csize, iw*ih*3./csize);

#ifdef T41_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t41_ctx->csizes_est+T41_NESTIMATORS*kb;
				for(int ke=1;ke<T41_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T41_NESTIMATORS;++ke)
			{
				float *sizes=t41_ctx->csizes_est+ke;
#ifndef T41_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T41_N_REC_ESTIMATORS+1), ke%(T41_N_REC_ESTIMATORS+1));
#else
				printf("E%3d ", ke);
#endif
				for(int kb=0;kb<24;++kb)
				{
					char c;
					if(ke==minidx[kb])
						c='*';
					else if(ke==maxidx[kb])
						c='L';
					else
						c=' ';
					printf("%8.2f %c", iw*ih/sizes[T41_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T41_NESTIMATORS*kb]/t41_ctx->csizes_est[T41_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T41_DISABLE_REC
				if(!((ke+1)%(T41_N_REC_ESTIMATORS+1)))
#else
				if(!((ke+1)%8))
#endif
				{
					printf("\n");
					printf("\t\t*         **        ***       ****      ****      ****      ****      *             *         *         *         *         **        ***       ****      *             *         **        ***       ****      ****      ****      ****      *\n");
					printf("\t\t0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7\n");
				}
			}
		}
#endif
	}
	t41_ctx_clear(&t41_ctx);
	for(int kl=0;kl<16;++kl)
		dlist_clear(lists+kl);
	free(buf2);
#if T41_NRANDPRED
	free(ebuf);
#endif
	return 1;
}
int t41_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T41 Dec  %s\n", g_buf);
	}
	T41Ctx *t41_ctx=t41_ctx_init();
	//char *gbuf=(char*)malloc((size_t)res<<2);//
	if(!t41_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	//memcpy(gbuf, guard, (size_t)res<<2);//
	//addbuf((unsigned char*)gbuf, iw, ih, 3, 4, 128);//
	//colortransform_ycocb_fwd(gbuf, iw, ih);//

	ABACSIMDDec ctx;
	ptrdiff_t csizes[16]={0}, csize_sum=0;
	const unsigned char *starts[16], *ends[16];
	for(int kl=0;kl<16;++kl)
	{
		starts[kl]=data+csize_sum;
		memcpy(csizes+kl, starts[kl], 4);
		csize_sum+=csizes[kl];
		starts[kl]+=4;//skip size
		ends[kl]=data+csize_sum;
	}
	abacsimd_dec_init(&ctx, starts, ends);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t41_ctx_init(t41_ctx);
	
	int bw=iw>>2, bh=ih>>2;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx[]=
				{
					(iw* ky      +kx     )<<2|kc,
					(iw* ky      +kx+bw  )<<2|kc,
					(iw* ky      +kx+bw*2)<<2|kc,
					(iw* ky      +kx+bw*3)<<2|kc,
					(iw*(ky+bh  )+kx     )<<2|kc,
					(iw*(ky+bh  )+kx+bw  )<<2|kc,
					(iw*(ky+bh  )+kx+bw*2)<<2|kc,
					(iw*(ky+bh  )+kx+bw*3)<<2|kc,
					(iw*(ky+bh*2)+kx     )<<2|kc,
					(iw*(ky+bh*2)+kx+bw  )<<2|kc,
					(iw*(ky+bh*2)+kx+bw*2)<<2|kc,
					(iw*(ky+bh*2)+kx+bw*3)<<2|kc,
					(iw*(ky+bh*3)+kx     )<<2|kc,
					(iw*(ky+bh*3)+kx+bw  )<<2|kc,
					(iw*(ky+bh*3)+kx+bw*2)<<2|kc,
					(iw*(ky+bh*3)+kx+bw*3)<<2|kc,
				};
				t41_ctx_get_context(t41_ctx, (char*)buf, iw, ih, kc, kx, ky, 0, bw, 0, bh);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t41_ctx_estimate_p0(t41_ctx, kc, kb);
					
					__m256i bits;
					abacsimd_dec(&ctx, &t41_ctx->p0, &bits);
					for(int kl=0;kl<16;++kl)
						buf[idx[kl]]|=bits.m256i_u16[kl]<<kb;

					t41_ctx_update(t41_ctx, kc, kb, &bits);
				}
				for(int kl=0;kl<16;++kl)
				{
					buf[idx[kl]]+=128;//unsigned -> signed

					//if((char)buf[idx[kl]]!=gbuf[idx[kl]])//
					//	LOG_ERROR("Error");
				}
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, bh, 100.*(ky+1)/bh);
	}
	
#ifdef T41_APPLY_SPATIAL
	addbuf(buf, iw, ih, 3, 4, 128);//buffer is signed
	apply_transforms_inv(buf, iw, ih);
#else
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	t41_ctx_clear(&t41_ctx);
#if T41_NRANDPRED
	free(ebuf);
#endif
	//free(gbuf);//
	return 1;
}
#endif


//T42: Bold and brash (T39 with 'custom2' filter)

//	#define T42_APPLY_SPATIAL
//	#define T42_DISABLE_REC
//	#define T42_DISABLE_COUNTER
//	#define T42_PRINT_ESTIMATOR_CR

#define T42_LR (int)(0.07*0x10000+0.5)
#define T42_NMAPS 15	//14	31		135 HALF HOUR PER IMAGE

#define C2_REACH 2	//changing this requires re-training the filter
#define C2_NNB (C2_REACH*(C2_REACH+1)*2*6)
typedef struct Custom2ParamsStruct
{
	short c0[C2_NNB];//fixed 1.14 bit
	short c1[C2_NNB+2];
	short c2[C2_NNB+4];
} Custom2Params;
Custom2Params filter=
{
	{
		 0x0659,-0x0C78, 0x0ABB,-0x0BC8, 0x09FE,
		-0x1249, 0x04C8, 0x1A98, 0x08DB, 0x0042,
		 0x1515, 0x11F9,
		-0x0019,-0x000C, 0x000F, 0x0030,-0x001D,
		-0x0058,-0x00D4, 0x0077,-0x0044, 0x0088,
		-0x0005, 0x00C8,
		-0x0078,-0x0104, 0x02FA, 0x0215,-0x028F,
		 0x0310, 0x038F,-0x0784,-0x02A7, 0x045E,
		-0x0050,-0x01EC,
		-0x016D, 0x048E,-0x0358, 0x0530,-0x051C,
		 0x078F, 0x07EC, 0x0DD8, 0x0C5B, 0x01F6,
		-0x0236, 0x1AD8,
		 0x0017, 0x0004,-0x005A,-0x0017,-0x0007,
		 0x002F,-0x0036,-0x0070,-0x000E,-0x0075,
		-0x00A8,-0x0143,
		 0x0097, 0x022D, 0x0295, 0x00A1, 0x009C,
		 0x000E, 0x0554, 0x0936, 0x006F,-0x03F0,
		 0x01B3, 0x01AD,
	},
	{
		 0x0094,-0x07E8, 0x148A, 0x0BC2, 0x00FD,
		 0x0F7B, 0x0898,-0x1534,-0x013C,-0x0651,
		-0x065A,-0x0CF4,
		 0x021A, 0x0094,-0x00CC, 0x015A, 0x03D8,
		 0x032A, 0x0574, 0x1389,-0x01B4, 0x07F3,
		 0x0AF7, 0x09C5,
		 0x0486,-0x0688,-0x03B1, 0x080E, 0x0028,
		 0x057D, 0x019C,-0x063F,-0x02C7, 0x04CA,
		-0x0002,-0x05B9,
		 0x002C, 0x00BB,-0x0131, 0x025E,-0x03A7,
		 0x00D2, 0x0432, 0x046B,-0x0006, 0x03EA,
		-0x0124,-0x0337,
		 0x0061,-0x043B,-0x06F5,-0x0173,-0x04BD,
		-0x004C, 0x0777, 0x0E22, 0x0563,-0x038F,
		-0x01C7, 0x24CF,
		 0x0037, 0x0299,-0x06D2,-0x0380,-0x0230,
		-0x04DB,-0x05C4,-0x07A6,-0x05F3,-0x050F,
		-0x0344,-0x0341,
		-0x099F,-0x03E9,
	},
	{
		 0x0223,-0x01CE, 0x051C, 0x0094,-0x0005,
		-0x00C3, 0x0324,-0x011D, 0x025C,-0x00F6,
		 0x0644,-0x0A9D,
		 0x006C,-0x0026, 0x00BB, 0x009D,-0x002E,
		 0x00C9, 0x0013,-0x0137, 0x0033,-0x0006,
		 0x007E,-0x0039,
		-0x02F5,-0x01BC, 0x05CB, 0x0108,-0x033C,
		-0x05CE, 0x0254, 0x129B, 0x07DE, 0x03F8,
		 0x0AB4, 0x1A62,
		-0x00FB, 0x00D7,-0x0164, 0x00E0, 0x02A4,
		-0x00B2, 0x0196, 0x024C,-0x0037,-0x001F,
		 0x0149, 0x08D3,
		-0x0065, 0x006E,-0x002A,-0x002A, 0x0068,
		-0x0040, 0x0126, 0x0120,-0x006F, 0x003B,
		-0x0048, 0x00A6,
		 0x0156, 0x0089, 0x0304, 0x016F, 0x01D9,
		 0x01F9, 0x0C27, 0x1670, 0x074D,-0x0198,
		 0x02EF, 0x18E0,
		-0x03D3, 0x0385,-0x0176, 0x009C,
	},
};
#ifndef T42_DISABLE_REC
#define T42_N_REC_ESTIMATORS 6		//15
#define T42_NESTIMATORS ((T42_N_REC_ESTIMATORS+1)*T42_NMAPS)
#else
#define T42_NESTIMATORS T42_NMAPS
#endif
typedef struct T42NodeStruct
{
	int n[2];
#ifndef T42_DISABLE_REC
	unsigned short rec[T42_N_REC_ESTIMATORS];
#endif
} T42Node;
typedef struct T42CtxStruct
{
	int pred14;
	int context[T42_NMAPS];
	ArrayHandle maps[24][T42_NMAPS];//(256+512+1024+2048+4096+8192+16384+32768)*20*14 = 17.43 MB for 14 maps with 6 rec estimators
	T42Node *node[T42_NMAPS];

	int p0arr[T42_NESTIMATORS], p0_0, p0, p0rev;//p0_0 isn't clamped
	int weights[24][T42_NESTIMATORS];
	long long wsum;

	int nnodes;
#ifdef T42_PRINT_ESTIMATOR_CR
	float csizes_est[24*T42_NESTIMATORS];
#endif
} T42Ctx;
T42Ctx* t42_ctx_init()
{
	int val=0x8000;
	T42Node node0={{1, 1}};
#ifndef T42_DISABLE_REC
	for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	T42Ctx *ctx=(T42Ctx*)malloc(sizeof(T42Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T42Ctx));
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T42_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T42Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T42Node), sizeof(T42Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t42_ctx_clear(T42Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T42_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
#if 0
void t42_ctx_reset(T42Ctx *ctx, int hardreset)
{
	T42Node node0={{1, 1}};
#ifndef T42_DISABLE_REC
	for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	for(int k=0;k<24;++k)
	{
		if(hardreset)
		{
			for(int k2=0;k2<T42_NMAPS;++k2)
				memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T42Node), sizeof(T42Node));
		}
	}
	if(hardreset)
		ctx->nnodes=0;
}
#endif
void t42_ctx_get_context(T42Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
//#define LOAD(C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X)-x1)<(unsigned)(x2-x1)&&(unsigned)(ky-Y-y1)<(unsigned)(y2-y1)?buf[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
//#define LOAD(CO, XO, YO) (unsigned)(kc-CO)<3u&&(unsigned)(kx-(XO))<(unsigned)iw&&(unsigned)(ky-YO)<(unsigned)ih?buf[(iw*(ky-YO)+kx-(XO))<<2|(kc-CO)]:0
#if 1
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	int
		NNWW =LOAD(buf, 0,  2, 2),
		NNW  =LOAD(buf, 0,  1, 2),
		NN   =LOAD(buf, 0,  0, 2),
		NNE  =LOAD(buf, 0, -1, 2),
		NNEE =LOAD(buf, 0, -2, 2),
		NWW  =LOAD(buf, 0,  2, 1),
		NW   =LOAD(buf, 0,  1, 1),
		N    =LOAD(buf, 0,  0, 1),
		NE   =LOAD(buf, 0, -1, 1),
		NEE  =LOAD(buf, 0, -2, 1),
		WW   =LOAD(buf, 0,  2, 0),
		W    =LOAD(buf, 0,  1, 0),
		eNNWW=LOAD(ebuf, 0,  2, 2),
		eNNW =LOAD(ebuf, 0,  1, 2),
		eNN  =LOAD(ebuf, 0,  0, 2),
		eNNE =LOAD(ebuf, 0, -1, 2),
		eNNEE=LOAD(ebuf, 0, -2, 2),
		eNWW =LOAD(ebuf, 0,  2, 1),
		eNW  =LOAD(ebuf, 0,  1, 1),
		eN   =LOAD(ebuf, 0,  0, 1),
		eNE  =LOAD(ebuf, 0, -1, 1),
		eNEE =LOAD(ebuf, 0, -2, 1),
		eWW  =LOAD(ebuf, 0,  2, 0),
		eW   =LOAD(ebuf, 0,  1, 0),

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);
#if 0
	int W   =LOAD(buf, 0,  1, 0),
		NW  =LOAD(buf, 0,  1, 1),
		N   =LOAD(buf, 0,  0, 1),
		NE  =LOAD(buf, 0, -1, 1),
		NN  =LOAD(buf, 0,  0, 2),

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);
#endif

	int j=-1;

	//bit, channel-bitplane, compressibility			based on kodim13
	//
	//Orangeness:			best pred				worst pred
	// 0	0-0		*		(N+W)/2					0
	// 1	0-1		**		(N+W)/2					0
	// 2	0-2		***		(N+W)/2					NW+NE-NN
	// 3	0-3		****	(N+W)/2					NW+NE-NN
	// 4	0-4		****	W						NW+NE-NN
	// 5	0-5		****	0						NW+NE-NN
	// 6	0-6		****	0						NW+NE-NN
	// 7	0-7		*		NW+NE-NN				0
	//
	//Luma:
	// 8	1-0		*		NW+NE-NN				NW+NE-NN
	// 9	1-1		*		(W+N+m1)/3				NW+NE-NN
	//10	1-2		*		(W+N+m1)/3				NW+NE-NN
	//11	1-3		*		W						0
	//12	1-4		**		W						0
	//13	1-5		***		(N+W-NW + m2)>>1		NW+NE-NN
	//14	1-6		****	(W+N+m1)/3				NW+NE-NN
	//15	1-7		*		NW+NE-NN				0
	//
	//Blueness:
	//16	2-0		*		W						clamp4(N+m1-Nm1, N, m1, Nm1, NW)
	//17	2-1		**		N						clamp4(N+m1-Nm1, N, m1, Nm1, NW)
	//18	2-2		***		W						(N+W-NW + m1)>>1
	//19	2-3		****	(N+W-NW + m2)>>1		(N+W-NW + m1)>>1
	//20	2-4		****	(N+W-NW + m2)>>1		(N+W-NW + m1)>>1
	//21	2-5		****	m2						(N+W-NW + m1)>>1
	//22	2-6		****	0						(N+W-NW + m1)>>1
	//23	2-7		*		m2						0

	ctx->context[++j]=0;//0
	ctx->context[++j]=N;//1
	ctx->context[++j]=W;//2
	ctx->context[++j]=NW;//3
	ctx->context[++j]=m1;//4
	ctx->context[++j]=W+NE-N;//5
	ctx->context[++j]=count_W_N_m1?(W+N+m1)/count_W_N_m1:0;//6
	ctx->context[++j]=clamp4(N+W-NW, N, W, NW, NE);//7
	ctx->context[++j]=clamp4(N+m1-Nm1, N, m1, Nm1, NW);//8
	ctx->context[++j]=clamp4(W+m1-Wm1, W, m1, Wm1, NW);//9
	ctx->context[++j]=NW+NE-NN;//10
	ctx->context[++j]=(N+W-NW + m1)>>1;//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=(N+W-NW + m2)>>1;//13

#undef  LOAD
#define LOAD(BUF, C, X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|C]:0
	char nb[C2_NNB+4]=
	{
		LOAD(buf, 0, -2, -2),
		LOAD(buf, 0, -1, -2),
		LOAD(buf, 0,  0, -2),
		LOAD(buf, 0,  1, -2),
		LOAD(buf, 0,  2, -2),
		LOAD(buf, 0, -2, -1),
		LOAD(buf, 0, -1, -1),
		LOAD(buf, 0,  0, -1),
		LOAD(buf, 0,  1, -1),
		LOAD(buf, 0,  2, -1),
		LOAD(buf, 0, -2,  0),
		LOAD(buf, 0, -1,  0),
		LOAD(buf, 1, -2, -2),
		LOAD(buf, 1, -1, -2),
		LOAD(buf, 1,  0, -2),
		LOAD(buf, 1,  1, -2),
		LOAD(buf, 1,  2, -2),
		LOAD(buf, 1, -2, -1),
		LOAD(buf, 1, -1, -1),
		LOAD(buf, 1,  0, -1),
		LOAD(buf, 1,  1, -1),
		LOAD(buf, 1,  2, -1),
		LOAD(buf, 1, -2,  0),
		LOAD(buf, 1, -1,  0),
		LOAD(buf, 2, -2, -2),
		LOAD(buf, 2, -1, -2),
		LOAD(buf, 2,  0, -2),
		LOAD(buf, 2,  1, -2),
		LOAD(buf, 2,  2, -2),
		LOAD(buf, 2, -2, -1),
		LOAD(buf, 2, -1, -1),
		LOAD(buf, 2,  0, -1),
		LOAD(buf, 2,  1, -1),
		LOAD(buf, 2,  2, -1),
		LOAD(buf, 2, -2,  0),
		LOAD(buf, 2, -1,  0),
		LOAD(ebuf, 0, -2, -2),
		LOAD(ebuf, 0, -1, -2),
		LOAD(ebuf, 0,  0, -2),
		LOAD(ebuf, 0,  1, -2),
		LOAD(ebuf, 0,  2, -2),
		LOAD(ebuf, 0, -2, -1),
		LOAD(ebuf, 0, -1, -1),
		LOAD(ebuf, 0,  0, -1),
		LOAD(ebuf, 0,  1, -1),
		LOAD(ebuf, 0,  2, -1),
		LOAD(ebuf, 0, -2,  0),
		LOAD(ebuf, 0, -1,  0),
		LOAD(ebuf, 1, -2, -2),
		LOAD(ebuf, 1, -1, -2),
		LOAD(ebuf, 1,  0, -2),
		LOAD(ebuf, 1,  1, -2),
		LOAD(ebuf, 1,  2, -2),
		LOAD(ebuf, 1, -2, -1),
		LOAD(ebuf, 1, -1, -1),
		LOAD(ebuf, 1,  0, -1),
		LOAD(ebuf, 1,  1, -1),
		LOAD(ebuf, 1,  2, -1),
		LOAD(ebuf, 1, -2,  0),
		LOAD(ebuf, 1, -1,  0),
		LOAD(ebuf, 2, -2, -2),
		LOAD(ebuf, 2, -1, -2),
		LOAD(ebuf, 2,  0, -2),
		LOAD(ebuf, 2,  1, -2),
		LOAD(ebuf, 2,  2, -2),
		LOAD(ebuf, 2, -2, -1),
		LOAD(ebuf, 2, -1, -1),
		LOAD(ebuf, 2,  0, -1),
		LOAD(ebuf, 2,  1, -1),
		LOAD(ebuf, 2,  2, -1),
		LOAD(ebuf, 2, -2,  0),
		LOAD(ebuf, 2, -1,  0),
	};
	ctx->pred14=0;
	switch(kc)
	{
	case 0:
		for(int k=0;k<C2_NNB;++k)
			ctx->pred14+=filter.c0[k]*nb[k];
		break;
	case 1:
		nb[72]=LOAD(buf, 0, 0, 0);
		nb[73]=LOAD(ebuf, 0, 0, 0);
		for(int k=0;k<C2_NNB+2;++k)
			ctx->pred14+=filter.c1[k]*nb[k];
		break;
	case 2:
		nb[72]=LOAD(buf, 0, 0, 0);
		nb[73]=LOAD(ebuf, 0, 0, 0);
		nb[72]=LOAD(buf, 1, 0, 0);
		nb[73]=LOAD(ebuf, 1, 0, 0);
		for(int k=0;k<C2_NNB+4;++k)
			ctx->pred14+=filter.c2[k]*nb[k];
		break;
	}
	ctx->pred14+=1<<13;
	ctx->pred14>>=14;
	ctx->context[++j]=ctx->pred14;

	//kodim13
#if 0
	switch(kc)
	{
	case 0:
		ctx->pred14=(
			+0x00C6* NNWW-0x0267* NNW+0x0400* NN-0x04C4* NNE+0x0188* NNEE
			-0x0166* NWW -0x0013* NW +0x034D* N +0x05BF* NE +0x002C* NEE
			+0x03E6* WW  +0x05BA* W
			+0x0076*eNNWW+0x0228*eNNW-0x001D*eNN+0x0241*eNNE-0x0167*eNNEE
			+0x00DA*eNWW +0x0404*eNW +0x0468*eN -0x0008*eNE +0x00C5*eNEE
			+0x0077*eWW  +0x05E7*eW
		)>>12;
		break;
	case 1:
		ctx->pred14=(
			+0x0010* NNWW+0x0020* NNW-0x0163* NN-0x000E* NNE+0x012B* NNEE
			-0x000E* NWW +0x02A6* NW +0x02F1* N +0x0081* NE +0x01D5* NEE
			+0x0215* WW  +0x063A* W
			-0x0074*eNNWW-0x015D*eNNW-0x00CC*eNN-0x00EB*eNNE-0x0156*eNNEE
			-0x0167*eNWW +0x0011*eNW +0x04A7*eN +0x00A5*eNE -0x00B2*eNEE
			-0x01D2*eWW  +0x0650*eW
		)>>12;
		break;
	case 2:
		ctx->pred14=(
			+0x00AC* NNWW-0x02B4* NNW+0x021E* NN-0x0010* NNE+0x0036* NNEE
			-0x0257* NWW +0x0076* NW +0x054D* N +0x00F9* NE +0x00BF* NEE
			+0x02D3* WW  +0x07D5* W
			+0x001D*eNNWW+0x0125*eNNW+0x009D*eNN+0x00EA*eNNE+0x007F*eNNEE
			+0x00AA*eNWW +0x035E*eNW +0x0669*eN +0x03FD*eNE -0x0044*eNEE
			+0x005C*eWW  +0x0526*eW
		)>>12;
		break;
	}
	ctx->pred14=CLAMP(-128, ctx->pred14, 127);
	ctx->context[++j]=ctx->pred14;//14
#endif
	
	//CLIC16
#if 0
	switch(kc)
	{
	case 0:
		ctx->pred14=(
			+0x00C8* NNWW-0x01B9* NNW+0x01CB* NN+0x0170* NNE-0x00E7* NNEE
			-0x01DA* NWW +0x00A8* NW +0x03FD* N +0x01DD* NE +0x00AB* NEE
			+0x00A5* WW  +0x0900* W
			+0x00FF*eNNWW+0x0040*eNNW-0x02B6*eNN+0x000D*eNNE+0x0182*eNNEE
			+0x0031*eNWW +0x00DE*eNW +0x065B*eN +0x0220*eNE +0x0056*eNEE
			-0x02EB*eWW  +0x036E*eW
		)>>12;
		break;
	case 1:
		ctx->pred14=(
			+0x0080* NNWW-0x00DF* NNW-0x00FC* NN+0x0195* NNE+0x0052* NNEE
			+0x0146* NWW -0x021E* NW +0x050F* N +0x0285* NE -0x012A* NEE
			-0x00E5* WW  +0x0AD7* W
			-0x003E*eNNWW+0x0126*eNNW-0x011B*eNN-0x011B*eNNE+0x0049*eNNEE
			-0x008A*eNWW +0x0174*eNW +0x048C*eN +0x007B*eNE +0x025C*eNEE
			-0x0199*eWW  +0x0428*eW
		)>>12;
		break;
	case 2:
		ctx->pred14=(
			-0x00B6* NNWW-0x0040* NNW+0x009D* NN-0x00E4* NNE+0x0088* NNEE
			-0x0263* NWW +0x02C8* NW +0x0420* N +0x0300* NE +0x0034* NEE
			+0x0131* WW  +0x07F9* W
			+0x0139*eNNWW-0x006E*eNNW-0x00E6*eNN+0x0019*eNNE-0x0020*eNNEE
			-0x0099*eNWW -0x002F*eNW +0x0693*eN +0x014F*eNE +0x003C*eNEE
			-0x0124*eWW  +0x05CC*eW
		)>>12;
		break;
	}
	ctx->pred14=CLAMP(-128, ctx->pred14, 127);
	ctx->context[++j]=ctx->pred14;//14
#endif

	//ctx->context[++j]=clamp4((W+NE-N + NW+NE-NN)>>1, N, W, NW, NE);
	//ctx->context[++j]=Nm1+Wm1-NWm1;
	//ctx->context[++j]=Nm1;
	//ctx->context[++j]=Wm1;
	//ctx->context[++j]=NWm1;
#endif
	
#if 0
	//offsets are in NW direction
	int WWWWWW  =LOAD(buf, 0,  6, 0),
		WWWWW   =LOAD(buf, 0,  5, 0),
		WWWW    =LOAD(buf, 0,  4, 0),
		WWW     =LOAD(buf, 0,  3, 0),
		WW      =LOAD(buf, 0,  2, 0),
		W       =LOAD(buf, 0,  1, 0),
		NWWWW   =LOAD(buf, 0,  4, 1),
		NWWW    =LOAD(buf, 0,  3, 1),
		NWW     =LOAD(buf, 0,  2, 1),
		NW      =LOAD(buf, 0,  1, 1),
		N       =LOAD(buf, 0,  0, 1),
		NE      =LOAD(buf, 0, -1, 1),
		NEE     =LOAD(buf, 0, -2, 1),
		NEEE    =LOAD(buf, 0, -3, 1),
		NEEEE   =LOAD(buf, 0, -4, 1),
		NEEEEEE =LOAD(buf, 0, -6, 1),
		NNWWW   =LOAD(buf, 0,  3, 2),
		NNWW    =LOAD(buf, 0,  2, 2),
		NNW     =LOAD(buf, 0,  1, 2),
		NN      =LOAD(buf, 0,  0, 2),
		NNE     =LOAD(buf, 0, -1, 2),
		NNEE    =LOAD(buf, 0, -2, 2),
		NNEEE   =LOAD(buf, 0, -3, 2),
		NNNWWWW =LOAD(buf, 0,  4, 3),
		NNNWWW  =LOAD(buf, 0,  3, 3),
		NNNWW   =LOAD(buf, 0,  2, 3),
		NNNW    =LOAD(buf, 0,  1, 3),
		NNN     =LOAD(buf, 0,  0, 3),
		NNNE    =LOAD(buf, 0, -1, 3),
		NNNEE   =LOAD(buf, 0, -2, 3),
		NNNEEE  =LOAD(buf, 0, -3, 3),
		NNNNW   =LOAD(buf, 0,  1, 4),
		NNNN    =LOAD(buf, 0,  0, 4),
		NNNNE   =LOAD(buf, 0, -1, 4),
		NNNNN   =LOAD(buf, 0,  0, 5),
		NNNNNN  =LOAD(buf, 0,  0, 6),
		WWWWWWp1=LOAD(buf, 1,  6, 0),
		WWWWp1  =LOAD(buf, 1,  4, 0),
		WWWp1   =LOAD(buf, 1,  3, 0),
		WWp1    =LOAD(buf, 1,  2, 0),
		Wp1     =LOAD(buf, 1,  1, 0),
		p1      =LOAD(buf, 1,  0, 0),
		NWWp1   =LOAD(buf, 1,  2, 1),
		NWp1    =LOAD(buf, 1,  1, 1),
		Np1     =LOAD(buf, 1,  0, 1),
		NEp1    =LOAD(buf, 1, -1, 1),
		NEEp1   =LOAD(buf, 1, -2, 1),
		NNWWp1  =LOAD(buf, 1,  2, 2),
		NNp1    =LOAD(buf, 1,  0, 2),
		NNEp1   =LOAD(buf, 1, -1, 2),
		NNEEp1  =LOAD(buf, 1, -2, 2),
		NNNWp1  =LOAD(buf, 1,  1, 3),
		NNNp1   =LOAD(buf, 1,  0, 3),
		NNNEp1  =LOAD(buf, 1, -1, 3),
		NNNNp1  =LOAD(buf, 1,  0, 4),
		NNNNNNp1=LOAD(buf, 1,  0, 6),
		WWWWWWp2=LOAD(buf, 2,  6, 0),
		WWWWp2  =LOAD(buf, 2,  4, 0),
		WWWp2   =LOAD(buf, 2,  3, 0),
		WWp2    =LOAD(buf, 2,  2, 0),
		Wp2     =LOAD(buf, 2,  1, 0),
		p2      =LOAD(buf, 2,  0, 0),
		NWWp2   =LOAD(buf, 2,  2, 1),
		NWp2    =LOAD(buf, 2,  1, 1),
		Np2     =LOAD(buf, 2,  0, 1),
		NEp2    =LOAD(buf, 2, -1, 1),
		NEEp2   =LOAD(buf, 2, -2, 1),
		NNWWp2  =LOAD(buf, 2,  2, 2),
		NNp2    =LOAD(buf, 2,  0, 2),
		NNEp2   =LOAD(buf, 2, -1, 2),
		NNEEp2  =LOAD(buf, 2, -2, 2),
		NNNWp2  =LOAD(buf, 2,  1, 3),
		NNNp2   =LOAD(buf, 2,  0, 3),
		NNNEp2  =LOAD(buf, 2, -1, 3),
		NNNNp2  =LOAD(buf, 2,  0, 4),
		NNNNNNp2=LOAD(buf, 2,  0, 6);
	int j=-1;

	//ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;//#74
	//ctx->context[++j] = (N + W + 1) >> 1;//#71
	//ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;//#70
	//ctx->context[++j] = N;//#75

	ctx->context[++j]=0;
	ctx->context[++j] = clamp4(N + p1 - Np1, W, NW, N, NE);
	ctx->context[++j] = clamp4(N + p2 - Np2, W, NW, N, NE);
	ctx->context[++j] = (W + clamp4(NE * 3 - NNE * 3 + NNNE, W, N, NE, NEE)) / 2;
	ctx->context[++j] = clamp4((W + clip(NE * 2 - NNE)) / 2, W, NW, N, NE);
	ctx->context[++j] = (W + NEE) / 2;
	ctx->context[++j] = ((WWW - 4 * WW + 6 * W + (NE * 4 - NNE * 6 + NNNE * 4 - NNNNE)) / 4);
	ctx->context[++j] = ((-WWWW + 5 * WWW - 10 * WW + 10 * W + clamp4(NE * 4 - NNE * 6 + NNNE * 4 - NNNNE, N, NE, NEE, NEEE)) / 5);
	ctx->context[++j] = ((-4 * WW + 15 * W + 10 * (NE * 3 - NNE * 3 + NNNE) - (NEEE * 3 - NNEEE * 3 + NNNEEE)) / 20);
	ctx->context[++j] = ((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p1 - (Wp1 + (NEp1 * 2 - NNEp1)) / 2);
	ctx->context[++j] = ((W + (NE * 2 - NNE)) / 2 + p2 - (Wp2 + (NEp2 * 2 - NNEp2)) / 2);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p1 -(-3 * WWp1 + 8 * Wp1 + (NEEp1 * 2 - NNEEp1)) / 6);
	ctx->context[++j] = ((-3 * WW + 8 * W + (NEE * 2 - NNEE)) / 6 + p2 -(-3 * WWp2 + 8 * Wp2 + (NEEp2 * 2 - NNEEp2)) / 6);
	ctx->context[++j] = ((W + NEE) / 2 + p1 - (Wp1 + NEEp1) / 2);
	ctx->context[++j] = ((W + NEE) / 2 + p2 - (Wp2 + NEEp2) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p1 - (WWp1 + (NEEp1 * 2 - NNEEp1)) / 2);
	ctx->context[++j] = ((WW + (NEE * 2 - NNEE)) / 2 + p2 - (WWp2 + (NEEp2 * 2 - NNEEp2)) / 2);
	ctx->context[++j] = (WW + NEE - N + p1 - (WWp1 + NEEp1 - Np1));
	ctx->context[++j] = (WW + NEE - N + p2 - (WWp2 + NEEp2 - Np2));
	ctx->context[++j] = (W + N - NW);
	ctx->context[++j] = (W + N - NW + p1 - (Wp1 + Np1 - NWp1));
	ctx->context[++j] = (W + N - NW + p2 - (Wp2 + Np2 - NWp2));
	ctx->context[++j] = (W + NE - N);
	ctx->context[++j] = (N + NW - NNW);
	ctx->context[++j] = (N + NW - NNW + p1 - (Np1 + NWp1 - NNEp1));
	ctx->context[++j] = (N + NW - NNW + p2 - (Np2 + NWp2 - NNEp2));
	ctx->context[++j] = (N + NE - NNE);
	ctx->context[++j] = (N + NE - NNE + p1 - (Np1 + NEp1 - NNEp1));
	ctx->context[++j] = (N + NE - NNE + p2 - (Np2 + NEp2 - NNEp2));
	ctx->context[++j] = (N + NN - NNN);
	ctx->context[++j] = (N + NN - NNN + p1 - (Np1 + NNp1 - NNNp1));
	ctx->context[++j] = (N + NN - NNN + p2 - (Np2 + NNp2 - NNNp2));
	ctx->context[++j] = (W + WW - WWW);
	ctx->context[++j] = (W + WW - WWW + p1 - (Wp1 + WWp1 - WWWp1));
	ctx->context[++j] = (W + WW - WWW + p2 - (Wp2 + WWp2 - WWWp2));
	ctx->context[++j] = (W + NEE - NE);
	ctx->context[++j] = (W + NEE - NE + p1 - (Wp1 + NEEp1 - NEp1));
	ctx->context[++j] = (W + NEE - NE + p2 - (Wp2 + NEEp2 - NEp2));
	ctx->context[++j] = (NN + p1 - NNp1);
	ctx->context[++j] = (NN + p2 - NNp2);
	ctx->context[++j] = (NN + W - NNW);
	ctx->context[++j] = (NN + W - NNW + p1 - (NNp1 + Wp1 - NNEp1));
	ctx->context[++j] = (NN + W - NNW + p2 - (NNp2 + Wp2 - NNEp2));
	ctx->context[++j] = (NN + NW - NNNW);
	ctx->context[++j] = (NN + NW - NNNW + p1 - (NNp1 + NWp1 - NNNWp1));
	ctx->context[++j] = (NN + NW - NNNW + p2 - (NNp2 + NWp2 - NNNWp2));
	ctx->context[++j] = (NN + NE - NNNE);
	ctx->context[++j] = (NN + NE - NNNE + p1 - (NNp1 + NEp1 - NNNEp1));
	ctx->context[++j] = (NN + NE - NNNE + p2 - (NNp2 + NEp2 - NNNEp2));
	ctx->context[++j] = (NN + NNNN - NNNNNN);
	ctx->context[++j] = (NN + NNNN - NNNNNN + p1 - (NNp1 + NNNNp1 - NNNNNNp1));
	ctx->context[++j] = (NN + NNNN - NNNNNN + p2 - (NNp2 + NNNNp2 - NNNNNNp2));
	ctx->context[++j] = (WW + p1 - WWp1);
	ctx->context[++j] = (WW + p2 - WWp2);
	ctx->context[++j] = (WW + WWWW - WWWWWW);
	ctx->context[++j] = (WW + WWWW - WWWWWW + p1 - (WWp1 + WWWWp1 - WWWWWWp1));
	ctx->context[++j] = (WW + WWWW - WWWWWW + p2 - (WWp2 + WWWWp2 - WWWWWWp2));
	ctx->context[++j] = (N * 2 - NN + p1 - (Np1 * 2 - NNp1));
	ctx->context[++j] = (N * 2 - NN + p2 - (Np2 * 2 - NNp2));
	ctx->context[++j] = (W * 2 - WW + p1 - (Wp1 * 2 - WWp1));
	ctx->context[++j] = (W * 2 - WW + p2 - (Wp2 * 2 - WWp2));
	ctx->context[++j] = (N * 3 - NN * 3 + NNN);
	ctx->context[++j] = clamp4(N * 3 - NN * 3 + NNN, W, NW, N, NE);
	ctx->context[++j] = clamp4(W * 3 - WW * 3 + WWW, W, NW, N, NE);
	ctx->context[++j] = clamp4(N * 2 - NN, W, NW, N, NE);
	ctx->context[++j] = ((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW, W, NW, N, NN)) / 6);
	ctx->context[++j] = ((NNNEEE - 4 * NNEE + 6 * NE + (W * 4 - NW * 6 + NNW * 4 - NNNW)) / 4);
	ctx->context[++j] = (((N + 3 * NW) / 4) * 3 - ((NNW + NNWW) / 2) * 3 + (NNNWW * 3 + NNNWWW) / 4);
	ctx->context[++j] = ((W * 2 + NW) - (WW + 2 * NWW) + NWWW);
	ctx->context[++j] = ((W * 2 - NW) + (W * 2 - NWW) + N + NE) / 4;
	ctx->context[++j] = (N + W + 1) >> 1;
	ctx->context[++j] = (NEEEE + NEEEEEE + 1) >> 1;
	ctx->context[++j] = (WWWWWW + WWWW + 1) >> 1;
	ctx->context[++j] = ((W + N) * 3 - NW * 2) >> 2;
	ctx->context[++j] = N;
	ctx->context[++j] = NN;
	ctx->context[++j] = N + p1 - Np1;
	ctx->context[++j] = N + p2 - Np2;
	ctx->context[++j] = W + p1 - Wp1;
	ctx->context[++j] = W + p2 - Wp2;
	ctx->context[++j] = NW + p1 - NWp1;
	ctx->context[++j] = NW + p2 - NWp2;
	ctx->context[++j] = NE + p1 - NEp1;
	ctx->context[++j] = NE + p2 - NEp2;
	ctx->context[++j] = NN + p1 - NNp1;
	ctx->context[++j] = NN + p2 - NNp2;
	ctx->context[++j] = WW + p1 - WWp1;
	ctx->context[++j] = WW + p2 - WWp2;
	ctx->context[++j] = W + N - NW;
	ctx->context[++j] = W + N - NW + p1 - Wp1 - Np1 + NWp1;
	ctx->context[++j] = W + N - NW + p2 - Wp2 - Np2 + NWp2;
	ctx->context[++j] = W + NE - N;
	ctx->context[++j] = W + NE - N + p1 - Wp1 - NEp1 + Np1;
	ctx->context[++j] = W + NE - N + p2 - Wp2 - NEp2 + Np2;
	ctx->context[++j] = W + NEE - NE;
	ctx->context[++j] = W + NEE - NE + p1 - Wp1 - NEEp1 + NEp1;
	ctx->context[++j] = W + NEE - NE + p2 - Wp2 - NEEp2 + NEp2;
	ctx->context[++j] = N + NN - NNN;
	ctx->context[++j] = N + NN - NNN + p1 - Np1 - NNp1 + NNNp1;
	ctx->context[++j] = N + NN - NNN + p2 - Np2 - NNp2 + NNNp2;
	ctx->context[++j] = N + NE - NNE;
	ctx->context[++j] = N + NE - NNE + p1 - Np1 - NEp1 + NNEp1;
	ctx->context[++j] = N + NE - NNE + p2 - Np2 - NEp2 + NNEp2;
	ctx->context[++j] = N + NW - NNW;
	ctx->context[++j] = N + NW - NNW + p1 - Np1 - NWp1 + NNEp1;
	ctx->context[++j] = N + NW - NNW + p2 - Np2 - NWp2 + NNEp2;
	ctx->context[++j] = NE + NW - NN;
	ctx->context[++j] = NE + NW - NN + p1 - NEp1 - NWp1 + NNp1;
	ctx->context[++j] = NE + NW - NN + p2 - NEp2 - NWp2 + NNp2;
	ctx->context[++j] = NW + W - NWW;
	ctx->context[++j] = NW + W - NWW + p1 - NWp1 - Wp1 + NWWp1;
	ctx->context[++j] = NW + W - NWW + p2 - NWp2 - Wp2 + NWWp2;
	ctx->context[++j] = W * 2 - WW;
	ctx->context[++j] = W * 2 - WW + p1 - Wp1 * 2 + WWp1;
	ctx->context[++j] = W * 2 - WW + p2 - Wp2 * 2 + WWp2;
	ctx->context[++j] = N * 2 - NN;
	ctx->context[++j] = N * 2 - NN + p1 - Np1 * 2 + NNp1;
	ctx->context[++j] = N * 2 - NN + p2 - Np2 * 2 + NNp2;
	ctx->context[++j] = NW * 2 - NNWW;
	ctx->context[++j] = NW * 2 - NNWW + p1 - NWp1 * 2 + NNWWp1;
	ctx->context[++j] = NW * 2 - NNWW + p2 - NWp2 * 2 + NNWWp2;
	ctx->context[++j] = NE * 2 - NNEE;
	ctx->context[++j] = NE * 2 - NNEE + p1 - NEp1 * 2 + NNEEp1;
	ctx->context[++j] = NE * 2 - NNEE + p2 - NEp2 * 2 + NNEEp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p1 - Np1 * 3 + NNp1 * 3 - NNNp1;
	ctx->context[++j] = N * 3 - NN * 3 + NNN + p2 - Np2 * 3 + NNp2 * 3 - NNNp2;
	ctx->context[++j] = N * 3 - NN * 3 + NNN;
	ctx->context[++j] = (W + NE * 2 - NNE + 1) >> 1;
	ctx->context[++j] = (W + NE * 3 - NNE * 3 + NNNE+1) >> 1;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p1 - (Wp1 + NEp1 * 2 - NNEp1) / 2;
	ctx->context[++j] = (W + NE * 2 - NNE) / 2 + p2 - (Wp2 + NEp2 * 2 - NNEp2) / 2;
	ctx->context[++j] = NNE + NE - NNNE;
	ctx->context[++j] = NNE + W - NN;
	ctx->context[++j] = NNW + W - NNWW;
#endif
#undef LOAD
	for(int k=0;k<T42_NMAPS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
int t42_ctx_map_context(int *context, int kp, int workidx)//replacement for context[kp]
{
	return context[kp];

	//static const int rep[]={ 0,  0, 10, 10, 10, 10, 10,  0,     10, 10, 10,  0,  0, 10, 10,  0,      8,  8, 11, 11, 11, 11, 11,  0};
	//static const int sub[]={ 6,  6,  6,  6,  2,  0,  0, 10,     10,  6,  6,  2,  2, 13,  6, 10,      2,  1,  2, 13, 13, 12,  0, 12};
	//return context[kp==rep[workidx]?sub[workidx]:kp];
}
void t42_ctx_estimate_p0(T42Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	T42Node *node;
	for(int kp=0;kp<T42_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		int context=t42_ctx_map_context(ctx->context, kp, kc);
		ArrayHandle map=ctx->maps[workidx][kp];
		node=ctx->node[kp]=(T42Node*)array_at(&map, context);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T42_DISABLE_REC
		for(;k2<T42_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T42_NESTIMATORS;++k)
	{
#ifdef T42_DISABLE_COUNTER
		if(k%(T42_N_REC_ESTIMATORS+1))//
#endif
		{
			sum+=(long long)ctx->p0arr[k]*wk[k];
			ctx->wsum+=wk[k];
		}
	}
	//ctx->p0=ctx->wsum?(int)((sum+(ctx->wsum>>1))/ctx->wsum):0x8000;//same CR
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
	ctx->p0rev=ctx->p0;
}
void t42_ctx_update(T42Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
	
#ifdef T42_PRINT_ESTIMATOR_CR
	for(int k=0;k<T42_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T42_NESTIMATORS*workidx+k]+=bitsize;
		}
	}
#endif
	//bwd
	int *wk=ctx->weights[workidx];
	if(ctx->p0_0>=1&&ctx->p0_0<=0xFFFF)
	{
		int p_bit=bit?0x10000-ctx->p0:ctx->p0;
		long long dL_dp0=-(1LL<<32)/p_bit;//fixed 47.16 bit
		dL_dp0^=-bit;
		dL_dp0+=bit;
		for(int k=0;k<T42_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=T42_LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
#ifndef T42_DISABLE_REC
	static const int shifts[]=
	{
		//7, 7, 7, 6, 5, 0, 0, 8,//T42_N_REC_ESTIMATORS 1
		//8, 7, 7, 6, 6, 6, 7, 7,
		//7, 7, 7, 5, 4, 4, 2, 7,

		5, 5, 4, 2, 1, 0, 0, 5,//T42_N_REC_ESTIMATORS 6
		7, 5, 5, 4, 4, 5, 4, 6,
		5, 5, 5, 2, 1, 1, 1, 5,

		//6, 5, 5, 4, 2, 0, 0, 6,
		//9, 7, 6, 5, 4, 5, 4, 6,
		//6, 5, 5, 4, 4, 3, 1, 6,
	};
	//static const int bitrating[]=//higher is more compressible
	//{
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//	1, 1, 1, 1, 2, 3, 4, 1,
	//	1, 2, 3, 4, 4, 4, 4, 1,
	//};
#endif
	T42Node *node;
	for(int kp=0;kp<T42_NMAPS;++kp)
	{
		node=ctx->node[kp];
		++node->n[bit];
#ifndef T42_DISABLE_REC
		//static const int lgdens[]={0, 1, 2, 4, 8, 14};
		for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		{
			//int lgden=k+2;
			//int lgden=k+shifts[workidx];
			//int lgden=k+1;
			//int lgden=(k+1)<<1;//X
			//int lgden=k+((4-bitrating[workidx])<<1);
			//int lgden=k+4-bitrating[workidx];
			//int lgden=k+3;
			int lgden=k;
			//int lgden=lgdens[k];
			//int lgden=((k+1)<<1)-1;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);
		//ctx->context[kp]<<=1;
		//ctx->context[kp]|=bit;
	}
#ifdef T42_PROB_TWEAK
	ctx->proberrors[workidx]+=abs((bit<<16)-ctx->p0);
	ctx->hits[workidx]+=bit==(ctx->p0<0x8000);
	//if(ctx->hits[workidx]<0)
	//	printf("");
	//if(workidx==15)
	//	printf("");
	//ctx->hits[workidx]+=bit==(ctx->p0<0x8000);
#endif
}
int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T42 Enc  CUSTOM2  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);
	T42Ctx *t42_ctx=t42_ctx_init();
	if(!buf2||!ebuf||!t42_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	memset(ebuf, 0, (size_t)res<<2);
#ifdef T42_APPLY_SPATIAL
	apply_transforms_fwd(buf2, iw, ih);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//buffer is signed
#else
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd(buf2, iw, ih);
	//addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);//X  the buffer is signed
#endif

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	//int hits[24]={0};
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(t42_ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(t42_ctx, kc, kb);
					int bit=(buf2[idx]+128)>>kb&1;
					abac_enc(&ctx, t42_ctx->p0, bit);
					
					int prob=bit?0x10000-t42_ctx->p0:t42_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t42_ctx_update(t42_ctx, kc, kb, bit);
				}
				ebuf[idx]=buf2[idx]-t42_ctx->pred14;
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f\r", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev));
			//printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f%c", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev), loud==2?'\n':'\r');
			csize_prev=csize;
		}
		//if(!((ky+1)&127))
		//	t42_ctx_reset(&t42_ctx, 0);
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t42_ctx->nnodes*sizeof(T42Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
#if 0
		double csize=0;
		for(int k=0;k<24;++k)
		{
			if(!(k&7))
			{
				printf("C%d\n", k>>3);
				csize=0;
			}
			printf("bit %2d  size %14f  CR %14f  H %7d %10lf%%\n", k&7, csizes[k]/8, iw*ih/csizes[k], hits[k], 100.*hits[k]/(iw*ih));
			csize+=csizes[k]/8;
			if(!((k+1)&7))
				printf("C%d  size %14lf  CR %14lf\n\n", k>>3, csize, iw*ih/csize);
		}
		printf("Total %lld    CR %lf    WH %d*%d  bitplane %g\n", list.nobj, 3.*iw*ih/list.nobj, iw, ih, iw*ih/8.);
		printf("\n");
#endif
		
		float chsizes[4]={0};
		//printf("\t\tC0\t\t\t\tC1\t\t\t\tC2\n\n");
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=7;kb>=0;--kb)
		{
			printf("B%d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc<<3|kb;
				float size=csizes[idx];
				//printf("       %12.3f %12.2f", iw*ih/size, hits[idx]);
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);

#ifdef T42_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t42_ctx->csizes_est+T42_NESTIMATORS*kb;
				for(int ke=1;ke<T42_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T42_NESTIMATORS;++ke)
			{
				float *sizes=t42_ctx->csizes_est+ke;
#ifndef T42_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T42_N_REC_ESTIMATORS+1), ke%(T42_N_REC_ESTIMATORS+1));
#else
				printf("E%3d ", ke);
#endif
				for(int kb=0;kb<24;++kb)
				{
					char c;
					if(ke==minidx[kb])
						c='*';
					else if(ke==maxidx[kb])
						c='L';
					else
						c=' ';
					printf("%8.2f %c", iw*ih/sizes[T42_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T42_NESTIMATORS*kb]/t42_ctx->csizes_est[T42_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T42_DISABLE_REC
				if(!((ke+1)%(T42_N_REC_ESTIMATORS+1)))
#else
				if(!((ke+1)%8))
#endif
				{
					printf("\n");
					printf("\t\t*         **        ***       ****      ****      ****      ****      *             *         *         *         *         **        ***       ****      *             *         **        ***       ****      ****      ****      ****      *\n");
					printf("\t\t0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7             0         1         2         3         4         5         6         7\n");
				}
			}
		}
#endif
	}
	t42_ctx_clear(&t42_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();

	//int debug_index=0;
	char *ebuf=(char*)malloc((size_t)res<<2);
	T42Ctx *t42_ctx=t42_ctx_init();
	if(!ebuf||!t42_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t42_ctx_init(t42_ctx);
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(t42_ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(t42_ctx, kc, kb);
					
					int bit=abac_dec(&ctx, t42_ctx->p0);
					buf[idx]|=bit<<kb;

					t42_ctx_update(t42_ctx, kc, kb, bit);
				}
				buf[idx]+=128;//unsigned -> signed
				ebuf[idx]=buf[idx]-t42_ctx->pred14;
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
		//if(!((ky+1)&127))
		//	t42_ctx_reset(&t42_ctx, 0);
	}
	t42_ctx_clear(&t42_ctx);
	
#ifdef T42_APPLY_SPATIAL
	addbuf(buf, iw, ih, 3, 4, 128);//buffer is signed
	apply_transforms_inv(buf, iw, ih);
#else
	//addbuf(buf, iw, ih, 3, 4, 128);//X  the buffer is signed
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	return 1;
}