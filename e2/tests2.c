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


//T42: The answer (T39 with 'custom3' filter)
#if 1
//	#define T42_DISABLE_REC
//	#define T42_DISABLE_COUNTER
//	#define T42_PRINT_ESTIMATOR_CR

#define T42_LR (int)(0.07*0x10000+0.5)
#define T42_NMAPS 15

#define C3_REACH 3	//changing this requires re-training the filter
#define C3_NNB (C3_REACH*(C3_REACH+1)*4)
#define C3_NPARAMS (C3_NNB*9+6)

#ifndef T42_DISABLE_REC
#define T42_N_REC_ESTIMATORS 6		//15
#define T42_NESTIMATORS ((T42_N_REC_ESTIMATORS+1)*T42_NMAPS)
#else
#define T42_NESTIMATORS T42_NMAPS
#endif
typedef struct Custom3Struct//arbitrary size filter
{
	int reach, pred;
	short *params;
	char *errors;
} Custom3;
typedef struct Custom3ParamsStruct
{
	short c00[C3_NNB  ], c01[C3_NNB  ], c02[C3_NNB];//fixed 1.14 bit
	short c10[C3_NNB+2], c11[C3_NNB  ], c12[C3_NNB];
	short c20[C3_NNB+2], c21[C3_NNB+2], c22[C3_NNB];
} Custom3Params;
short filter[]=
{
	//CUSTOM3-r2 CLIC16
#if 0
	 0x00FE, 0x0398, 0x0001,-0x01E1,-0x01E5,-0x0E5E,-0x0019, 0x02FD,-0x000A, 0x0201,
	 0x0070, 0x03C6,-0x1CA5, 0x055A, 0x32FC,-0x041F,-0x0234, 0x0EEA, 0x0059, 0x010A,
	-0x009D,-0x0BB0, 0x2C77, 0x0407,
	
	-0x0001, 0x0028,-0x0001,-0x0014, 0x0003,-0x0043, 0x0001, 0x0001,-0x0005,-0x0001,
	 0x0001, 0x0001,-0x002C,-0x0024, 0x0028,-0x0012,-0x000C, 0x00AF, 0x0001,-0x0004,
	-0x0001, 0x0001, 0x0013, 0x0055,
	
	 0x0001, 0x0023,-0x0065, 0x0001,-0x003E, 0x00BA, 0x005C, 0x0034,-0x0076, 0x0130,
	 0x003B,-0x02A8, 0x00B9, 0x0294,-0x013B, 0x0486, 0x015E, 0x0075, 0x0025, 0x0027,
	-0x004A, 0x01EF,-0x009B, 0x0163,
	
	-0x00F6,-0x0049, 0x0003, 0x00B7, 0x0001, 0x0344,-0x0423,-0x0029, 0x0003, 0x0003,
	 0x0275,-0x00E2, 0x00C3,-0x0022, 0x00D4,-0x10BF, 0x0387,-0x05B6,-0x001F, 0x020F,
	-0x0007, 0x0165, 0x047A,-0x0673,-0x07B7,-0x2563,
	
	 0x0128, 0x0165,-0x01EF,-0x0273,-0x000B,-0x093C,-0x00B8, 0x02DA, 0x013D, 0x0126,
	 0x00AC, 0x010C,-0x0AEE, 0x049A, 0x1FA8, 0x1030, 0x03AE, 0x03E1, 0x00D5,-0x0002,
	 0x0147,-0x040C, 0x2575, 0x150A,
	
	 0x009B, 0x002E, 0x0005, 0x02EE, 0x00D8, 0x00CF, 0x0054, 0x00AC, 0x011C, 0x02EC,
	-0x0174,-0x0178, 0x00DF, 0x01A9,-0x00A3, 0x0C93,-0x0205, 0x083B, 0x0247,-0x03A3,
	 0x00F7, 0x01D7,-0x029B,-0x0239,
	
	-0x0057, 0x008A,-0x0003, 0x017C,-0x0041, 0x0099,-0x0003,-0x00D9,-0x0274, 0x0200,
	 0x0001, 0x00EB,-0x0073,-0x0018,-0x0065,-0x008F, 0x0440,-0x0458,-0x00F1, 0x0147,
	 0x00E4, 0x0299,-0x00D9,-0x0060, 0x005F,-0x0883,
	
	-0x0015,-0x0005, 0x0011, 0x0000, 0x0003,-0x0056, 0x0003,-0x0003,-0x00A3, 0x0023,
	-0x0001, 0x0001,-0x0059, 0x002D,-0x0019, 0x004F, 0x0042,-0x0087,-0x0003,-0x003E,
	-0x0003,-0x0001, 0x002F, 0x010E, 0x0085, 0x02C4,
	
	 0x0036, 0x012B,-0x0083, 0x00E0,-0x002F,-0x0359,-0x00D6, 0x0166, 0x0359,-0x00B4,
	 0x00D7, 0x0073,-0x0C88, 0x0E2C, 0x2021, 0x14BA, 0x02D0, 0x0C3C,-0x0034, 0x005D,
	 0x00B8,-0x04C2, 0x2481, 0x1435,
#endif

	//CUSTOM3-r3 kodim13
#if 1
	-0x0085, 0x0000, 0x000A, 0x0116, 0x0004, 0x01C3, 0x0008, 0x0176, 0x0002, 0x0143, 0x0025, 0x004B, 0x0007, 0x0026,
	 0x0034, 0x000F, 0x0342, 0x024E,-0x0B11, 0x0676, 0x113C,-0x048D,-0x1126, 0x06FC, 0x03F1,-0x03C2,-0x0005, 0x008F,
	-0x00B5,-0x006A,-0x0768, 0x041C,-0x0236, 0x0EDC, 0x0E5C, 0x0E9E, 0x1718,-0x01DB, 0x0117, 0x0210, 0x006C, 0x029A,
	 0x0011, 0x0001, 0x1035,-0x0069, 0x1667, 0x15E1,

	 0x0004, 0x0006, 0x0008, 0x000A,-0x000C,-0x0002,-0x0003, 0x0001,-0x0006,-0x0006, 0x0002, 0x0007, 0x0001, 0x0001,
	-0x0001,-0x0006,-0x0004, 0x0000, 0x0010, 0x0006, 0x0010, 0x0005,-0x000A, 0x0012, 0x0001, 0x0002,-0x0001,-0x0002,
	-0x0004,-0x0008, 0x001C, 0x0000,-0x0008, 0x002A, 0x000A, 0x0009, 0x0005, 0x0002,-0x0002,-0x0018,-0x0001,-0x0009,
	-0x0002, 0x0002,-0x000C, 0x000C, 0x0003,-0x0037,

	-0x0001, 0x0006,-0x0001,-0x000A,-0x0003,-0x0001,-0x0004, 0x0004, 0x0001, 0x000B, 0x0002, 0x002F, 0x0007,-0x0009,
	-0x0006, 0x004B, 0x0003, 0x005D, 0x000F,-0x0008,-0x0006, 0x0031,-0x0002, 0x0037, 0x0002, 0x00E7, 0x0000, 0x0005,
	-0x0023, 0x0004,-0x004C, 0x0016, 0x0056, 0x0031, 0x0004, 0x004F,-0x0003, 0x004A,-0x0007,-0x0026, 0x0005,-0x0004,
	 0x0008,-0x0002, 0x003F,-0x004B,-0x006D,-0x0069,

	-0x0095,-0x01A1,-0x0118,-0x014B, 0x012C,-0x024E, 0x0022, 0x038B, 0x00FC,-0x0030, 0x00B4,-0x0233,-0x0001,-0x0001,
	 0x0003,-0x042B, 0x0873,-0x02EF,-0x005A,-0x00C5, 0x03C6,-0x071E,-0x03A3,-0x00E5, 0x055E,-0x00E8, 0x0002, 0x0006,
	 0x000B, 0x0003,-0x006D,-0x000E,-0x0489,-0x0973,-0x03A3,-0x0F78,-0x0100,-0x01DA,-0x02E0, 0x0297, 0x0002,-0x0020,
	-0x0084, 0x01D0, 0x005D, 0x0171, 0x001C,-0x0AB5,-0x0554,-0x09B2,

	-0x000C, 0x0029,-0x0003,-0x0029, 0x0001,-0x0009, 0x0001, 0x02B2, 0x00A9,-0x0046,-0x0020,-0x018A,-0x006E,-0x0075,
	 0x0014, 0x0005,-0x0197,-0x00FE, 0x0286,-0x060E,-0x048B,-0x019F,-0x00F2,-0x01F4, 0x03EE,-0x0679, 0x0000,-0x0002,
	 0x0007,-0x0035, 0x0088,-0x0675, 0x097E, 0x01A4, 0x0769, 0x165E, 0x0274,-0x000A, 0x0ABF,-0x0642, 0x0001, 0x017B,
	 0x01D4, 0x00E1, 0x07C3,-0x06D2, 0x1730, 0x1838,

	 0x00BD, 0x0009, 0x0005,-0x00F0,-0x0008, 0x00E2,-0x0010, 0x0173, 0x006A,-0x0054,-0x0079,-0x00BC,-0x009D,-0x0029,
	 0x0159,-0x0205, 0x00F2,-0x0111, 0x01AA, 0x0085, 0x0187,-0x0533, 0x0266,-0x0045,-0x0005, 0x01AB,-0x0003, 0x012D,
	-0x0006, 0x0175,-0x00C1,-0x0134,-0x0154,-0x08F1, 0x0002,-0x0C83,-0x03D1,-0x0777, 0x0094,-0x019E, 0x0007,-0x00D0,
	-0x00B7, 0x016D, 0x0234,-0x0397,-0x038F,-0x053D,

	 0x01D1,-0x01B2,-0x02B3, 0x0152, 0x020F,-0x019E,-0x0138, 0x01B9, 0x0089, 0x00D6,-0x0100, 0x003A,-0x0035,-0x0008,
	-0x008B, 0x0098, 0x020E,-0x0066,-0x01B0, 0x0053, 0x02DB, 0x0068,-0x015C, 0x015B, 0x01B3,-0x0063,-0x00E0,-0x0105,
	-0x0006, 0x0000,-0x01E3,-0x01D1, 0x0392,-0x030B,-0x0340, 0x003E,-0x0003,-0x01BB, 0x0003,-0x0369, 0x0081,-0x00B6,
	-0x0069,-0x0079,-0x0012,-0x015A,-0x026E,-0x018A, 0x0216,-0x00A7,

	 0x0000,-0x0002,-0x0004,-0x0001, 0x0001, 0x0006,-0x0008, 0x0027, 0x0006, 0x0027,-0x0006, 0x0010,-0x0007,-0x0008,
	 0x0000, 0x0001, 0x0063,-0x0025, 0x0000,-0x0017, 0x0013,-0x0023, 0x0013,-0x0018, 0x005C,-0x0055,-0x0001, 0x000B,
	-0x0006, 0x0038,-0x0006, 0x0020, 0x0003, 0x0045,-0x001B, 0x0012,-0x0001,-0x0066, 0x0049,-0x003B,-0x0003, 0x0002,
	-0x0001, 0x0005, 0x0103,-0x0061,-0x0022, 0x00DC,-0x01B6, 0x00F3,

	 0x0008, 0x0003, 0x0003, 0x0083, 0x0000, 0x017D,-0x0005, 0x0059,-0x000C, 0x01B7,-0x0002,-0x0069, 0x0001,-0x0001,
	-0x005D, 0x00E9, 0x0251, 0x005D,-0x0B4E, 0x03DA, 0x09BF,-0x004F,-0x0268, 0x0469, 0x01B9,-0x011F, 0x000D, 0x00D4,
	-0x0007, 0x0270,-0x090C, 0x03D9,-0x0082, 0x09EE, 0x1819, 0x1572, 0x0436, 0x0C0E, 0x01EF, 0x009F, 0x00EE, 0x00D4,
	 0x0007, 0x00B4, 0x0B3C, 0x0214, 0x1E9A, 0x1552,
#endif
};
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
	ArrayHandle maps[24][T42_NMAPS];//3*(256+512+1024+2048+4096+8192+16384+32768)*15*sizeof(T42Node) = 56.03 MB for 15 maps with 6 rec estimators
	T42Node *node[T42_NMAPS];

	int p0arr[T42_NESTIMATORS], p0_0, p0;//p0_0 isn't clamped
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
T42Ctx* t42_ctx_copy(T42Ctx *ctx)
{
	T42Ctx *ctx2=(T42Ctx*)malloc(sizeof(T42Ctx));
	if(!ctx2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx2, 0, sizeof(T42Ctx));
	memcpy(ctx2->weights, ctx->weights, sizeof(ctx2->weights));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T42_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T42Node, ctx2->maps[k][k2], 0, nnodes, 0, 0);
			memcpy(ctx2->maps[k][k2]->data, ctx->maps[k][k2]->data, ctx->maps[k][k2]->count*ctx->maps[k][k2]->esize);
			ctx2->nnodes+=nnodes;
		}
	}
	return ctx2;
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
static int custom3_loadnb(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-C3_REACH;ky2<0;++ky2)
	{
		for(int kx2=-C3_REACH;kx2<=C3_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-C3_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static int custom3_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-15;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((__m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((__m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
	__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
	s2=_mm_add_epi32(s2, hi);
	s2=_mm_hadd_epi32(s2, s2);
	int s3=_mm_extract_epi32(s2, 0);
	for(;k<count;++k)
		s3+=a[k]*b[k];
	return s3;
}
void t42_ctx_get_context(T42Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
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

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);

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

	{//CUSTOM3
		Custom3Params *params=(Custom3Params*)filter;
		const short *coeffs[]=
		{
			params->c00, params->c01, params->c02,
			params->c10, params->c11, params->c12,
			params->c20, params->c21, params->c22,

			//0, 0, 0,
			//0, 0, 0,
			//0, 0, 0,
		};
		short nb[3][C3_NNB+2]={0};
		int count[3], idx, idx2;
		for(int kc=0;kc<3;++kc)
			count[kc]=custom3_loadnb(buf, ebuf, iw, ih, kc, kx, ky, nb[kc]);
		idx=(iw*ky+kx)<<2;
		idx2=0;
		switch(kc)
		{
		case 1:
			nb[0][C3_NNB  ]=buf[idx];
			nb[0][C3_NNB+1]=ebuf[idx];
			count[0]+=2;
			//++idx;
			idx2+=3;
			break;
		case 2:
			nb[0][C3_NNB  ]=buf [idx  ];
			nb[0][C3_NNB+1]=ebuf[idx  ];
			nb[1][C3_NNB  ]=buf [idx|1];
			nb[2][C3_NNB+1]=ebuf[idx|1];
			count[0]+=2;
			count[1]+=2;
			//idx+=2;
			idx2+=6;
			break;
		}
		ctx->pred14=0;
		for(int kc=0;kc<3;++kc)
		{
			if(coeffs[idx2+kc])
				ctx->pred14+=custom3_dot(coeffs[idx2+kc], nb[kc], count[kc]);
		}

		ctx->pred14+=1<<13;
		ctx->pred14>>=14;
		ctx->pred14=CLAMP(-128, ctx->pred14, 127);
		ctx->context[++j]=ctx->pred14;
	}
#undef LOAD
	for(int k=0;k<T42_NMAPS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
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
		int context=ctx->context[kp];
		ArrayHandle map=ctx->maps[workidx][kp];
		//node=ctx->node[kp]=ARRAY_AT(T42Node, map, context);
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
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
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
	T42Node *node;
	for(int kp=0;kp<T42_NMAPS;++kp)
	{
		node=ctx->node[kp];
		++node->n[bit];
#ifndef T42_DISABLE_REC
		for(int k=0;k<T42_N_REC_ESTIMATORS;++k)
		{
			int lgden=k;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);
	}
}
void t42_explore(void *ctx0)
{
#if 1
	const char *prednames[]=
	{
		"0               ",
		"N               ",
		"W               ",
		"NW              ",
		"m1              ",
		"W+NE-N          ",
		"av(W,N,m1)      ",
		"clamp(N+W-NW)   ",
		"clamp(N+m1-Nm1) ",
		"clamp(W+m1-Wm1) ",
		"clamp(NW+NE-NN) ",
		"(N+W-NW + m1)>>1",
		"m2              ",
		"(N+W-NW + m2)>>1",
		"CUSTOM3-r3      ",
	};
	printf("About to overwrite dump.txt at:\n");
	system("cd");
	pause();
	FILE *f=fopen("dump.txt", "w");
	T42Ctx *ctx=(T42Ctx*)ctx0;
	for(int kw=0;kw<24;++kw)
	{
		//printf("C%d B%d %lf\n", kw>>3, kw&7, csizes[kw]);
		fprintf(f, "C%d B%d\n", kw>>3, kw&7);
		for(int kp=0;kp<T42_NMAPS;++kp)
		{
			ArrayHandle tree=ctx->maps[kw][kp];
			//printf("C%d  B%d  %s\t", kw>>3, kw&7, prednames[kp]);
			fprintf(f, "%s\t", prednames[kp]);
			for(int kv=0;kv<(int)tree->count;++kv)
			{
				T42Node *node=(T42Node*)array_at(&tree, kv);
				int sum=node->n[0]+node->n[1];
				if(sum>2)
				{
					int c=(node->n[1]*15+(sum>>1))/sum;
					fprintf(f, "%X", c);
				}
				else
					fprintf(f, "-");
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
	fclose(f);
	printf("Done.\n");
#endif
#if 0
	int total_count=0, total_use=0;
	for(int kw=0;kw<24;++kw)
	{
		for(int kp=0;kp<T42_NMAPS;++kp)
		{
			ArrayHandle tree=ctx->maps[kw][kp];
			int count=(int)tree->count;
			int unused=0;
			for(int kv=0;kv<count;++kv)
			{
				T42Node *node=(T42Node*)array_at(&tree, kv);
				unused+=node->n[0]==1&&node->n[1]==1;
			}
			printf("C%d  B%d  P%3d  %7d/%7d\n", kw>>3, kw&7, kp, unused, count);

			total_count+=count;
			total_use+=count-unused;
		}
	}
	printf("Usage %d/%d\n", total_use, total_count);
#endif
}
void t42_freectx(void **ctx)
{
	t42_ctx_clear((T42Ctx**)ctx);
}
int t42_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T42 Enc  CUSTOM  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);

	//T42Ctx *ctx=(T42Ctx*)*ctx0;
	//if(!ctx)
	//	ctx=t42_ctx_init();
	T42Ctx *ctx=t42_ctx_init();

	if(!buf2||!ebuf||!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	//*ctx0=t42_ctx_copy(ctx);
	memcpy(buf2, src, (size_t)res<<2);
	memset(ebuf, 0, (size_t)res<<2);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd(buf2, iw, ih);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ac;
	abac_enc_init(&ac, &list);
	
	float csizes[24]={0};
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(ctx, kc, kb);
					unsigned short p0=ctx->p0;

					int bit=(buf2[idx]+128)>>kb&1;
					abac_enc(&ac, p0, bit);
					
					int prob=bit?0x10000-p0:p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t42_ctx_update(ctx, kc, kb, bit);
				}
				ebuf[idx]=buf2[idx]-ctx->pred14;
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f\r", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev));
			csize_prev=csize;
		}
	}
	abac_enc_flush(&ac);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)ctx->nnodes*sizeof(T42Node)/(1024*1024));
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

		//t42_explore(t42_ctx, csizes);

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
	t42_ctx_clear(&ctx);
	//t42_ctx_clear(&t42_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t42_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();

	char *ebuf=(char*)malloc((size_t)res<<2);
	//if(!*ctx0)
	//	*ctx0=t42_ctx_init();
	T42Ctx *ctx=t42_ctx_init();
	if(!ebuf||!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	//T42Ctx *ctx=(T42Ctx*)*ctx0;

	ABACDecContext ac;
	abac_dec_init(&ac, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t42_ctx_get_context(ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t42_ctx_estimate_p0(ctx, kc, kb);
					
					int bit=abac_dec(&ac, ctx->p0);
					buf[idx]|=bit<<kb;

					t42_ctx_update(ctx, kc, kb, bit);
				}
				buf[idx]+=128;//unsigned -> signed
				ebuf[idx]=buf[idx]-ctx->pred14;
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
	}
	t42_ctx_clear(&ctx);
	
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	return 1;
}
#endif


//T43: Wisdom of the crowd
#if 1
//	#define T43_DISABLE_REC
//	#define T43_DISABLE_COUNTER
//	#define T43_PRINT_ESTIMATOR_CR

#define T43_REACH 3	//changing this requires re-training the filter
#define T43_NNB (T43_REACH*(T43_REACH+1)*4)
#define T43_NPARAMS (T43_NNB*9+6)
#define T43_NFILTERS	1

#define T43_LR (int)(0.07*0x10000+0.5)
#define T43_NMAPS (14+T43_NFILTERS)

#ifndef T43_DISABLE_REC
#define T43_N_REC_ESTIMATORS 6		//15
#define T43_NESTIMATORS ((T43_N_REC_ESTIMATORS+1)*T43_NMAPS)
#else
#define T43_NESTIMATORS T43_NMAPS
#endif
typedef struct T43C3ParamsStruct
{
	short c00[T43_NNB  ], c01[T43_NNB  ], c02[T43_NNB];//fixed 1.14 bit
	short c10[T43_NNB+2], c11[T43_NNB  ], c12[T43_NNB];
	short c20[T43_NNB+2], c21[T43_NNB+2], c22[T43_NNB];
} T43C3Params;
short t43_allparams[]=
{
	//kodim13	r3	overall good filter for natural content
#if 1
	-0x0088, 0x0003, 0x0007, 0x0113, 0x0004, 0x01C1, 0x000A, 0x017A,-0x0002, 0x013F, 0x0023, 0x0046, 0x0009, 0x0027,
	 0x0033, 0x000B, 0x0342, 0x024B,-0x0B12, 0x067A, 0x113B,-0x048A,-0x1121, 0x06FF, 0x03F2,-0x03C7,-0x0006, 0x0093,
	-0x00B9,-0x0066,-0x0764, 0x041D,-0x0236, 0x0ED6, 0x0E5D, 0x0E9F, 0x1717,-0x01DA, 0x0116, 0x020B, 0x0071, 0x029E,
	 0x0009,-0x0001, 0x1035,-0x0065, 0x1666, 0x15E2,
	
	 0x0008, 0x0004, 0x0009, 0x000A,-0x0011, 0x0004,-0x0001,-0x0001,-0x0007,-0x0001, 0x0001, 0x000B,-0x0007, 0x0001,
	 0x0004,-0x0007,-0x0008,-0x0005, 0x0012, 0x0007, 0x0012, 0x0008,-0x0004, 0x000F, 0x0001,-0x0003, 0x0000, 0x0000,
	-0x0002,-0x000A, 0x0019, 0x0001,-0x000F, 0x0029, 0x000A, 0x0008, 0x0004, 0x0001,-0x0004,-0x0020, 0x0005,-0x0007,
	-0x0001, 0x0002,-0x000F, 0x0009, 0x0006,-0x0035,
	
	-0x0005, 0x0007, 0x0001,-0x000D,-0x0005, 0x0002,-0x0004, 0x0003, 0x0004, 0x000B, 0x0007, 0x002E, 0x0007,-0x0008,
	-0x0006, 0x0049, 0x0003, 0x005A, 0x000F,-0x000D,-0x0009, 0x002D,-0x0009, 0x003A,-0x0001, 0x00EA, 0x0001, 0x0009,
	-0x0028, 0x0005,-0x004E, 0x001A, 0x005D, 0x0033, 0x0009, 0x004B,-0x0005, 0x004A,-0x0007,-0x002A, 0x0008,-0x0004,
	 0x0008,-0x0006, 0x0041,-0x0044,-0x006E,-0x006E,
	
	-0x0093,-0x01A5,-0x0117,-0x0147, 0x012B,-0x0246, 0x0023, 0x038B, 0x00FA,-0x0031, 0x00B3,-0x0230,-0x0004, 0x0000,
	 0x0002,-0x0424, 0x086D,-0x02F0,-0x0058,-0x00C3, 0x03CB,-0x0723,-0x03A2,-0x00E5, 0x055D,-0x00E7, 0x0003, 0x0004,
	 0x000D, 0x0000,-0x006D,-0x0011,-0x0486,-0x0974,-0x03A4,-0x0F6D,-0x00FE,-0x01DA,-0x02DB, 0x029A, 0x0006,-0x001E,
	-0x0087, 0x01CD, 0x005C, 0x016E, 0x0020,-0x0AB5,-0x0557,-0x09AF,
	
	-0x000E, 0x002B,-0x0003,-0x002D, 0x0001,-0x000A, 0x0002, 0x02B2, 0x00AA,-0x0049,-0x0024,-0x018F,-0x006B,-0x0079,
	 0x0015, 0x0004,-0x019B,-0x0102, 0x0283,-0x060D,-0x0486,-0x019B,-0x00EF,-0x01F0, 0x03EC,-0x067A, 0x0006, 0x0000,
	 0x000B,-0x003A, 0x0085,-0x0678, 0x0978, 0x01A4, 0x0769, 0x1661, 0x0272,-0x0009, 0x0AB8,-0x0642, 0x0002, 0x017F,
	 0x01D1, 0x00E1, 0x07C2,-0x06CF, 0x1732, 0x1839,
	
	 0x00BA, 0x000D, 0x000A,-0x00F0,-0x0004, 0x00E1,-0x0010, 0x017C, 0x006D,-0x0054,-0x007C,-0x00BD,-0x0099,-0x0028,
	 0x0152,-0x0203, 0x00EE,-0x010D, 0x01A6, 0x007D, 0x0189,-0x0533, 0x0264,-0x004C,-0x0008, 0x01A8,-0x0003, 0x012D,
	 0x0002, 0x0174,-0x00C1,-0x0134,-0x0151,-0x08F6, 0x0001,-0x0C82,-0x03D1,-0x077A, 0x0098,-0x01A5, 0x000B,-0x00CF,
	-0x00AE, 0x016C, 0x0235,-0x0394,-0x038D,-0x053D,
	
	 0x01D5,-0x01B3,-0x02B4, 0x0150, 0x020E,-0x01A0,-0x0137, 0x01B4, 0x008F, 0x00D8,-0x00FE, 0x0035,-0x0037,-0x0006,
	-0x0092, 0x0097, 0x0211,-0x0067,-0x01AE, 0x0054, 0x02DB, 0x0068,-0x015F, 0x015F, 0x01B5,-0x0064,-0x00E2,-0x0104,
	-0x0002,-0x0003,-0x01DF,-0x01D2, 0x0396,-0x030E,-0x0341, 0x0043, 0x0001,-0x01B9, 0x0004,-0x0368, 0x0082,-0x00BB,
	-0x006F,-0x0079,-0x0010,-0x0157,-0x026D,-0x0186, 0x0217,-0x00A7,
	
	-0x0004, 0x0000,-0x0003, 0x0001,-0x0004, 0x0005,-0x0003, 0x002D, 0x0008, 0x0026,-0x0009, 0x0011,-0x0008, 0x0000,
	-0x0002,-0x0004, 0x0067,-0x0028,-0x0002,-0x0016, 0x0016,-0x0025, 0x001B,-0x001A, 0x005B,-0x0057,-0x0004, 0x000A,
	-0x0007, 0x0037,-0x0006, 0x0023,-0x0004, 0x0041,-0x001C, 0x000F,-0x0002,-0x0067, 0x004A,-0x003A,-0x0006, 0x0004,
	-0x0001,-0x0002, 0x0105,-0x0065,-0x0021, 0x00D5,-0x01B2, 0x00F6,
	
	 0x0005, 0x0005, 0x0005, 0x0084,-0x0003, 0x017B,-0x0009, 0x005A,-0x000A, 0x01B1,-0x0003,-0x0068, 0x0001,-0x0002,
	-0x005F, 0x00E9, 0x0255, 0x005F,-0x0B4E, 0x03DF, 0x09BE,-0x004B,-0x026A, 0x046A, 0x01BE,-0x011A, 0x0010, 0x00D4,
	-0x0004, 0x0271,-0x0911, 0x03DC,-0x0084, 0x09F3, 0x1817, 0x156F, 0x0435, 0x0C0F, 0x01EF, 0x0099, 0x00E8, 0x00D8,
	 0x000B, 0x00B1, 0x0B3E, 0x0212, 0x1EA0, 0x1551,
#endif

	//kodim13	r2
#if 0
	 0x02A1, 0x016D,-0x0B30, 0x06A8, 0x1193,-0x0466,-0x10CD, 0x07A6, 0x03FF,-0x0403,
	-0x07A6, 0x040B,-0x01AB, 0x0C07, 0x0E4E, 0x0E0E, 0x172F,-0x0173, 0x0099, 0x020E,
	 0x102B,-0x0056, 0x1622, 0x1506,
	
	-0x0026, 0x0020,-0x0009, 0x0029, 0x000E,-0x0009, 0x001B, 0x0005, 0x0003, 0x0004,
	 0x0011, 0x000C,-0x004D, 0x0012, 0x0028,-0x0039, 0x002F, 0x0009,-0x001C, 0x0011,
	 0x0016,-0x001D, 0x0000,-0x003B,
	
	 0x000B,-0x000D, 0x0042,-0x0069,-0x0010, 0x0035,-0x001F, 0x0039, 0x003E, 0x008A,
	-0x001C,-0x006F, 0x005E,-0x0009, 0x00D8, 0x0045,-0x0098, 0x002A,-0x0027,-0x0043,
	 0x0009,-0x004A,-0x00F4,-0x004E,
	
	 0x0A93,-0x0596,-0x0122, 0x00D9, 0x0033,-0x05BD,-0x0220,-0x01F5, 0x0494,-0x0119,
	-0x002E,-0x0013,-0x03C2,-0x0A7C,-0x02DE,-0x1262,-0x00AE,-0x00DC,-0x01FD, 0x005B,
	 0x0037, 0x01FE, 0x0099,-0x0BE7,-0x0531,-0x08F1,
	
	-0x009C,-0x00F4, 0x02A1,-0x04E6,-0x05EA,-0x0272, 0x00C4,-0x0452, 0x03FE,-0x0579,
	 0x004C,-0x069D, 0x0A7A, 0x0155, 0x0737, 0x1648, 0x029E,-0x009D, 0x0AC8,-0x070F,
	 0x0790,-0x06CC, 0x1773, 0x18F5,
	
	 0x0100, 0x014D, 0x0273,-0x0071, 0x0064,-0x0584, 0x01D2,-0x0049, 0x001B, 0x00AF,
	 0x014C,-0x01C1,-0x027E,-0x0809,-0x012C,-0x0AAF,-0x0305,-0x05EA, 0x00D2,-0x0021,
	 0x0280,-0x025F,-0x0508,-0x06E8,
	
	 0x021F,-0x0035,-0x0206, 0x0085, 0x03D8, 0x014E,-0x0275, 0x0094, 0x0274,-0x00CB,
	-0x01CB,-0x00BC, 0x03FA,-0x0324,-0x04B4, 0x01A2, 0x001B,-0x01B5,-0x0047,-0x01CD,
	-0x002D,-0x0079,-0x037E,-0x01C1, 0x0275,-0x005D,
	
	-0x0011, 0x000E, 0x0019,-0x0032, 0x0020,-0x004C, 0x0048,-0x003C, 0x002F,-0x002D,
	 0x002A,-0x002E, 0x0032, 0x0014, 0x000D,-0x0018,-0x0035,-0x006B, 0x003E,-0x003A,
	 0x00E7,-0x0052,-0x000F, 0x00B1,-0x01CA, 0x00EA,
	
	 0x0248, 0x0081,-0x0B21, 0x041F, 0x09B4,-0x0030,-0x0274, 0x03F7, 0x01BB,-0x00DB,
	-0x0A8D, 0x038C, 0x00D8, 0x0910, 0x1850, 0x154D, 0x042B, 0x0C04, 0x01E7, 0x0078,
	 0x0BC8, 0x00B2, 0x1EC3, 0x14F2,
#endif

	//24 mediocre filters
#if 0
	//01
	 0x1101, 0x07CD, 0x1432, 0x1735, 0x089C, 0x0B84,
	 0x1231, 0x11D7,
	 0x0138,-0x00E6,-0x00B0, 0x000C, 0x0017,-0x0091,
	-0x00A9, 0x0011,
	 0x00FF,-0x03F9,-0x0334,-0x0491,-0x0114, 0x0015,
	 0x03CB,-0x0636,

	 0x0314,-0x07BD,-0x0673,-0x1449,-0x0094, 0x0435,
	-0x03A1,-0x1A78,-0x0A5A,-0x0CAF,
	-0x0577, 0x0ADE, 0x12EA, 0x14DA, 0x02F3,-0x049A,
	 0x282B, 0x129D,
	 0x034C,-0x0702,-0x0C54,-0x06B7,-0x05B7, 0x0C30,
	-0x0A90,-0x0BB3,

	-0x003E,-0x02F6,-0x01FE,-0x09D9,-0x0258, 0x0071,
	 0x02AD,-0x09BD, 0x001A,-0x083A,
	 0x00E6,-0x00AD,-0x0119, 0x0029, 0x000B,-0x0071,
	-0x0001,-0x00A2,-0x0034,-0x00C1,
	 0x06FA, 0x0585, 0x11AC, 0x1311, 0x099A, 0x0AC2,
	 0x1A3D, 0x0B09,

	 
	//02
	-0x3F75, 0x18BD, 0x3FF9,-0x0FB2,-0x0037, 0x151D,
	 0x3F67,-0x0E8B,
	-0x0004, 0x0005,-0x0026,-0x0009, 0x004A,-0x00C1,
	-0x0004, 0x00BC,
	-0x0015, 0x0560,-0x0005,-0x147B,-0x0092,-0x0363,
	 0x0025,-0x13CA,
	
	 0x02AA, 0x03F5,-0x0048, 0x03D2,-0x023C, 0x075A,
	-0x0110, 0x0200, 0x0013, 0x0C4D,
	-0x3761, 0x0F2F, 0x3B0D,-0x1DA1,-0x0153,-0x0189,
	 0x3D63,-0x14D7,
	 0x022D,-0x0390,-0x0265,-0x072C,-0x0128, 0x03CF,
	-0x001D,-0x082E,
	
	-0x0010,-0x026D,-0x0013,-0x07EC,-0x0003,-0x075D,
	 0x0011,-0x0689, 0x002C,-0x1050,
	-0x002F, 0x0014,-0x000C, 0x0086,-0x0001,-0x0021,
	-0x0002,-0x000F, 0x001E,-0x009E,
	-0x3F98, 0x1A1C, 0x3FF9,-0x19AC, 0x0064, 0x0DC1,
	 0x3F72,-0x17B6,
	
	 
	//03
	 0x0C67, 0x046D, 0x12A0, 0x0666, 0x10C3,-0x0064,
	 0x102A, 0x0AD0,
	-0x0042,-0x0032,-0x005C,-0x0085, 0x006D,-0x005A,
	 0x003A,-0x0148,
	 0x00DD, 0x0128,-0x00F0,-0x005F, 0x02BC,-0x005F,
	-0x02BD, 0x01B6,
	
	-0x00DF,-0x0059,-0x09C8,-0x002A,-0x017D,-0x0124,
	-0x02F0, 0x0039, 0x0F16, 0x007F,
	 0x119A, 0x042B, 0x160B, 0x05F1, 0x0457,-0x032F,
	 0x1407, 0x19EA,
	-0x01BC, 0x000B,-0x03BF, 0x000F,-0x000B, 0x013A,
	 0x0574,-0x0288,
	
	 0x0078,-0x0028,-0x010F,-0x048F, 0x0470,-0x029E,
	 0x0060,-0x049E,-0x0449,-0x03DC,
	-0x03DB, 0x01C0, 0x01E5,-0x01EB, 0x00E9,-0x0193,
	 0x0392,-0x025A,-0x0288, 0x01DC,
	 0x190D,-0x0043, 0x0CA8, 0x0E5D, 0x119E, 0x0003,
	 0x089F, 0x119E,
	
	 
	//04
	 0x0ECA, 0x08C3, 0x0D7D, 0x1844, 0x1145, 0x00FB,
	 0x12E1, 0x1281,
	 0x03F4,-0x00DF,-0x03B7, 0x0263, 0x01DB,-0x00F7,
	-0x01BE, 0x00C2,
	-0x02E3,-0x038A, 0x0819,-0x152C,-0x02DA, 0x02FC,
	-0x0261,-0x094B,
	
	-0x0317,-0x0109,-0x0306,-0x024D,-0x052A, 0x00D0,
	-0x0389,-0x011F, 0x0E31,-0x0074,
	 0x0FE4, 0x07FD, 0x0DAE, 0x1865, 0x10A3,-0x0472,
	 0x11A4, 0x1A7A,
	-0x002C, 0x0167,-0x0242, 0x02E6,-0x01E4,-0x0010,
	 0x0423,-0x0156,
	
	 0x0211,-0x0181, 0x042D,-0x0425,-0x000E,-0x001E,
	-0x0215,-0x03EC,-0x0432,-0x0E1D,
	 0x005B, 0x0022,-0x0001,-0x004D, 0x003F,-0x003C,
	-0x01B9, 0x0152, 0x0112,-0x0031,
	 0x0CC9,-0x0164, 0x0F95, 0x0387, 0x133D, 0x0179,
	 0x104C, 0x043B,
	
	 
	//05
	 0x0D8C, 0x092C, 0x0809, 0x17B8, 0x1541, 0x03D6,
	 0x1527, 0x117C,
	 0x01BE,-0x006B,-0x01AC, 0x0049, 0x0160,-0x0139,
	-0x014F, 0x00B1,
	 0x0594, 0x0079,-0x055B, 0x00BD, 0x054C,-0x0250,
	-0x0586, 0x02B1,
	
	-0x0134,-0x017F, 0x04C3,-0x08B7,-0x082E, 0x0600,
	-0x01B3,-0x078A, 0x0602,-0x01D2,
	 0x0A30, 0x1047, 0x0F63, 0x2243, 0x0ACD,-0x00E2,
	 0x1A3E, 0x1DEF,
	-0x03DC,-0x07AC, 0x032B,-0x066E, 0x001E, 0x0345,
	 0x006D,-0x0C60,
	
	 0x036A,-0x0061,-0x030D,-0x074A, 0x026B,-0x0157,
	-0x00B4,-0x05C3,-0x0238,-0x02C2,
	-0x001E,-0x0141, 0x00B3,-0x0239,-0x0146, 0x003F,
	-0x00B7,-0x0172, 0x0172,-0x0264,
	 0x0440, 0x05EC, 0x1585, 0x1324, 0x08AC, 0x0BC9,
	 0x1D52, 0x1069,
	
	 
	//06
	 0x09EF, 0x09F9, 0x15B9, 0x0CE0, 0x0B9A, 0x0981,
	 0x1492, 0x1155,
	 0x00F7,-0x00B7,-0x0022, 0x00FF,-0x012F, 0x00E1,
	 0x0052, 0x0058,
	 0x031C,-0x0132,-0x07F4,-0x002C, 0x0277,-0x02ED,
	 0x0268,-0x0537,
	
	 0x014C, 0x0589,-0x07AE, 0x0685,-0x0373, 0x0BDF,
	 0x0D0F,-0x08B1,-0x0526, 0x0692,
	 0x0A35,-0x01D9, 0x0D3C, 0x067A, 0x0AA6,-0x047C,
	 0x1C94, 0x21FB,
	-0x02B2,-0x08C8, 0x0157,-0x10D4,-0x0014,-0x0426,
	-0x02B4,-0x106D,
	
	-0x03B0, 0x016C,-0x0072,-0x0A97, 0x071F,-0x052F,
	 0x0001,-0x0776,-0x04D0,-0x062E,
	 0x01D4,-0x009F,-0x01AC, 0x006F, 0x007A,-0x013D,
	-0x00A2, 0x0018,-0x00D8,-0x00B0,
	 0x07AB, 0x0780, 0x0A1D, 0x1254, 0x0F54, 0x04FB,
	 0x1B77, 0x0F7C,
	
	 
	//07
	 0x1279, 0x05AD, 0x0BCD, 0x13CD, 0x0EEF, 0x06F1,
	 0x12E8, 0x1346,
	-0x01ED,-0x0059, 0x04BF,-0x036F,-0x01D9, 0x0193,
	-0x00EA, 0x0090,
	-0x0132,-0x0030,-0x00C6,-0x0351,-0x0119, 0x02AE,
	 0x031C,-0x0335,
	
	 0x021F,-0x0496,-0x0295,-0x0614, 0x042F,-0x05E7,
	-0x0462,-0x070C, 0x016E,-0x0679,
	 0x0AE1, 0x0C8D, 0x0F00, 0x1B7D, 0x0D4F, 0x0601,
	 0x181E, 0x1DB9,
	-0x0043,-0x0424, 0x01F0,-0x049B, 0x006A,-0x00F3,
	-0x00D6,-0x02BA,
	
	 0x03FB,-0x0488, 0x0368,-0x04DE,-0x050F, 0x02A3,
	-0x0460,-0x058C, 0x0237,-0x0474,
	-0x0138,-0x0064, 0x0413, 0x000D,-0x0207, 0x0292,
	 0x0372, 0x009C,-0x0471, 0x0532,
	 0x0F15, 0x0943, 0x11D4, 0x149E, 0x12F0, 0x03A1,
	 0x0C0C, 0x17E2,
	
	 
	//08
	 0x0BA2, 0x073D, 0x121E, 0x1966, 0x157A, 0x0400,
	 0x0C23, 0x1356,
	 0x02E7,-0x006E,-0x0307, 0x0236, 0x01FE,-0x022D,
	-0x0219, 0x01DC,
	-0x0266,-0x03C1,-0x0035,-0x089F, 0x02ED,-0x0353,
	-0x00AB,-0x0500,
	
	 0x096E,-0x00FD,-0x055A,-0x1026, 0x0249, 0x0B66,
	-0x0741,-0x066D,-0x0576, 0x04B7,
	-0x01B6, 0x007A, 0x1572, 0x1549, 0x0D8B,-0x0210,
	 0x1CE8, 0x1169,
	-0x0781, 0x050C, 0x060A,-0x099F, 0x05C3, 0x0689,
	-0x0A41,-0x03BE,
	
	 0x0580,-0x0795,-0x05B7,-0x0969, 0x0184,-0x0346,
	-0x02D1,-0x0A6D, 0x0271,-0x0E1A,
	 0x00B0,-0x015F,-0x0005,-0x00B3,-0x0018,-0x007E,
	 0x003C,-0x00CA,-0x008A,-0x0078,
	 0x0F3F, 0x0902, 0x0F00, 0x1CA8, 0x1394, 0x00C6,
	 0x0EEE, 0x1500,
	
	 
	//09
	 0x0DA9, 0x0389, 0x10EB, 0x0B4E, 0x0D44, 0x054A,
	 0x139B, 0x0AF8,
	 0x0237,-0x01E0,-0x0293, 0x016A, 0x0066,-0x017D,
	-0x0005,-0x00CE,
	 0x05B0,-0x01C8,-0x005D,-0x0131,-0x01DB, 0x010C,
	-0x03AD,-0x00DC,
	
	-0x013B, 0x03FB, 0x0656,-0x0ADE,-0x006F, 0x05E4,
	-0x0463,-0x0720,-0x0017,-0x0440,
	 0x1032,-0x00F1, 0x0D88, 0x10C0, 0x0CFB,-0x03C2,
	 0x14CE, 0x1276,
	 0x001A, 0x009F,-0x04FD,-0x0349, 0x021B,-0x00BC,
	 0x02C2,-0x08D1,
	
	 0x01FD,-0x00A3,-0x0061,-0x0205,-0x00E7, 0x0133,
	-0x011C,-0x0269, 0x00A5, 0x00BC,
	 0x03C5,-0x0274,-0x0154,-0x0038, 0x0102,-0x026D,
	-0x0412, 0x0249, 0x00AF,-0x00F0,
	 0x0ECB, 0x0305, 0x125B, 0x0F5D, 0x0F73, 0x05C7,
	 0x0FB0, 0x11D1,
	
	 
	//10
	 0x112D, 0x0402, 0x0C0C, 0x124B, 0x088B, 0x0B22,
	 0x1374, 0x0D33,
	 0x00F4,-0x0187,-0x00A3,-0x0092,-0x01BE, 0x000C,
	 0x0021,-0x018C,
	 0x04C1,-0x04A0,-0x0319,-0x01C3,-0x05DB, 0x02D8,
	-0x0087,-0x0460,
	
	 0x0147,-0x04E5, 0x069E,-0x0945,-0x04B3, 0x0631,
	 0x01B9,-0x0A7B,-0x0588, 0x00A2,
	 0x080E, 0x01F3, 0x1164, 0x18DE, 0x0ABE, 0x0272,
	 0x1B4A, 0x0F29,
	 0x0327,-0x02A6, 0x04CA,-0x06E3,-0x074C, 0x02CF,
	-0x01D3,-0x0404,
	
	 0x0662,-0x020D,-0x0313,-0x01BE, 0x003A, 0x0124,
	-0x033A,-0x0007,-0x0427, 0x033E,
	-0x002E,-0x0126,-0x00F6,-0x01BF, 0x001D,-0x011A,
	-0x00DE,-0x011E, 0x0128,-0x0186,
	 0x0B5D, 0x0817, 0x112B, 0x12D1, 0x0DDF, 0x073E,
	 0x12BA, 0x1076,
	
	 
	//11
	 0x0A91, 0x0BF8, 0x1281, 0x1598, 0x10D7, 0x089C,
	 0x1186, 0x16E5,
	 0x0119,-0x006D, 0x00B0,-0x0030,-0x0040, 0x0060,
	-0x018C, 0x018A,
	-0x0007,-0x0565,-0x0011,-0x0ABD,-0x009E,-0x014A,
	-0x0084,-0x08ED,
	
	-0x00C0, 0x000E,-0x0350,-0x0206,-0x0321, 0x049F,
	 0x049B,-0x049D, 0x025A,-0x008D,
	 0x0E31,-0x0118, 0x1178, 0x0965, 0x0EDA,-0x08A1,
	 0x110D, 0x2059,
	 0x04EB, 0x00F2,-0x0256, 0x072B, 0x027B, 0x055D,
	-0x03CA,-0x07FC,
	
	-0x0290,-0x02F5, 0x058B,-0x0D5F, 0x02B3, 0x0044,
	-0x027C,-0x034F,-0x02FF,-0x09A5,
	 0x00A2,-0x01A4,-0x0040,-0x0099,-0x018C, 0x00C0,
	-0x002B,-0x01F5, 0x014F,-0x021B,
	 0x0BD9, 0x082A, 0x0B91, 0x1371, 0x1334, 0x020C,
	 0x154B, 0x1046,
	
	 
	//12
	 0x0F35, 0x0436, 0x0EA0, 0x0D74, 0x0DA1, 0x0556,
	 0x125B, 0x0DA9,
	-0x0192,-0x015D, 0x042C,-0x03C8,-0x0271, 0x01D7,
	-0x00A2,-0x00D7,
	 0x004A,-0x002D,-0x00E5,-0x02A8,-0x01AA,-0x0024,
	-0x0052,-0x009F,
	
	 0x04E0,-0x0140, 0x00B1, 0x03A6,-0x0156,-0x004E,
	-0x016E,-0x025E,-0x035A, 0x08BB,
	 0x0612, 0x056E, 0x0A38, 0x062F, 0x173B,-0x0D92,
	 0x1828, 0x1BBD,
	-0x0096, 0x0039, 0x0102,-0x003B,-0x000E, 0x015F,
	-0x003C,-0x03A8,
	
	 0x06CE,-0x031D,-0x0550, 0x0288,-0x00C8, 0x0156,
	-0x00DE,-0x00B2,-0x0028, 0x0295,
	 0x02B2, 0x008B,-0x043E, 0x01F6, 0x02BB,-0x031E,
	-0x00DB, 0x0115,-0x0074, 0x017A,
	 0x098B, 0x0624, 0x0FCA, 0x0C1E, 0x1218, 0x0291,
	 0x1482, 0x0CD5,
	
	 
	//13
	 0x0CC5, 0x08F3, 0x0AD3, 0x139E, 0x0E94, 0x08A7,
	 0x193F, 0x10BE,
	-0x007E,-0x0067, 0x01C9,-0x0137,-0x0106, 0x0101,
	-0x001B,-0x002E,
	-0x0009, 0x0131,-0x0444, 0x0101,-0x0003,-0x001F,
	 0x03EF,-0x036D,
	
	 0x048C,-0x0BD4,-0x012B, 0x0258, 0x055D,-0x02D6,
	-0x0769, 0x04E5,-0x0733, 0x0662,
	 0x0E0B, 0x0301, 0x0501, 0x1152, 0x10D9,-0x070B,
	 0x19F3, 0x1642,
	-0x050F,-0x0436, 0x0273,-0x085E, 0x03B7,-0x0A0C,
	-0x0164,-0x0726,
	
	-0x0192,-0x01BA,-0x034C,-0x05BF, 0x0172,-0x02EC,
	-0x0333,-0x065B, 0x068A,-0x0435,
	 0x007A, 0x0028,-0x00DE, 0x0007, 0x023A,-0x020C,
	 0x00D6, 0x0055,-0x02BA, 0x01DD,
	 0x06AA, 0x088A, 0x0B8F, 0x167E, 0x0F4F, 0x075A,
	 0x1DF6, 0x11DE,
	
	 
	//14
	 0x1085, 0x04B5, 0x0AC1, 0x0DF4, 0x10D3, 0x0368,
	 0x13C3, 0x0E46,
	 0x0142,-0x0046, 0x00E3,-0x00C4, 0x009B,-0x0083,
	-0x02B8, 0x0256,
	 0x01DA,-0x0059,-0x00D0,-0x00F3, 0x01BA,-0x00C2,
	-0x02A9,-0x001A,
	
	 0x0286, 0x001C,-0x01DF,-0x00F9,-0x061F, 0x0022,
	-0x08A8,-0x02EC, 0x0E9D,-0x0656,
	 0x0E9D, 0x062F, 0x0DAF, 0x129F, 0x0FE1,-0x03D9,
	 0x1349, 0x2047,
	 0x016D,-0x013C,-0x006A, 0x003B,-0x00DD, 0x00B4,
	 0x00F9,-0x0572,
	
	 0x0803,-0x02A1,-0x01BD,-0x0523, 0x048D,-0x02E7,
	-0x0509,-0x088A,-0x0561,-0x0D1B,
	-0x0066, 0x0026, 0x00C0,-0x0034,-0x0021,-0x0024,
	 0x007D, 0x002F,-0x00AE, 0x00DA,
	 0x1113, 0x0356, 0x0E00, 0x0CA0, 0x1224,-0x0053,
	 0x0E9B, 0x11AE,
	
	 
	//15
	-0x3FFF, 0x0DEA, 0x401A,-0x0C1B,-0x0007, 0x0F5F,
	 0x401C,-0x1395,
	 0x0037,-0x0086,-0x0001,-0x00B0,-0x002A,-0x009E,
	 0x0001,-0x0093,
	 0x0089, 0x00DF,-0x005C, 0x0018,-0x004E, 0x0082,
	 0x0041,-0x0258,
	
	 0x0002, 0x01F4,-0x0137,-0x03F2,-0x006A, 0x0174,
	-0x00B8,-0x0ADB, 0x02A0, 0x1732,
	-0x4003, 0x0E5A, 0x3FF0,-0x13EF, 0x002B, 0x07B3,
	 0x3FFC,-0x1F82,
	-0x0005, 0x0207, 0x0013,-0x03C7, 0x001A, 0x00EB,
	 0x0023,-0x06E8,
	
	 0x0003, 0x0141,-0x002B,-0x048A, 0x0000,-0x0184,
	 0x0029,-0x018F,-0x0006,-0x05F5,
	 0x0000,-0x0026, 0x0004, 0x0051,-0x0006,-0x0045,
	 0x0004,-0x001E, 0x0000,-0x002C,
	-0x3FD5,-0x0010, 0x4026,-0x003F,-0x000C,-0x0022,
	 0x3FC7,-0x0035,
	
	 
	//16
	-0x3DED, 0x18CB, 0x3DED,-0x20D7, 0x007A, 0x0E11,
	 0x3F83,-0x16DA,
	 0x000C, 0x0055, 0x0004, 0x0045, 0x0003, 0x0060,
	-0x0007, 0x002A,
	-0x0001, 0x025D,-0x0186,-0x03DC, 0x0068,-0x01C8,
	 0x013D,-0x0847,
	
	 0x008C, 0x06FC,-0x0200,-0x0200, 0x01E7,-0x033F,
	 0x0038,-0x0532, 0x0049, 0x0357,
	-0x3F7E,-0x0651, 0x3FFE,-0x3027,-0x00D1, 0x0324,
	 0x4067, 0x0B00,
	-0x0147, 0x01E7, 0x01DE,-0x0972,-0x0057,-0x03E7,
	 0x0075,-0x0C60,
	
	 0x001B, 0x003E,-0x0104,-0x037C, 0x00BE,-0x035D,
	-0x00EC,-0x04DE, 0x00F9,-0x0600,
	-0x0031, 0x002B,-0x000A,-0x0085, 0x0034,-0x00B1,
	-0x0001,-0x012D,-0x0007,-0x00E7,
	-0x3E9A, 0x16EA, 0x3DF2,-0x1B8B, 0x01B5, 0x0F38,
	 0x3EB1,-0x0E3E,
	
	 
	//17
	 0x106F, 0x03B7, 0x0AB2, 0x0D15, 0x0F7F, 0x027E,
	 0x1371, 0x08D3,
	 0x0127,-0x0147,-0x027C, 0x00E2, 0x0084,-0x012D,
	 0x00AE,-0x021D,
	 0x06D9,-0x022F,-0x0398,-0x0026, 0x000A,-0x0026,
	-0x0448,-0x021B,
	
	 0x017A,-0x037B,-0x0533,-0x0AE4, 0x01D2,-0x0322,
	 0x00ED,-0x04C6, 0x01D7,-0x063E,
	 0x0D8C, 0x073D, 0x10CD, 0x156E, 0x1257,-0x02F1,
	 0x0F4C, 0x18C1,
	-0x0188,-0x00E9,-0x008C,-0x03B6, 0x013B, 0x01AC,
	 0x018C,-0x0332,
	
	 0x099A,-0x03D1, 0x001D,-0x01DF, 0x0081, 0x007A,
	-0x06D0, 0x02BA,-0x04AF, 0x01D4,
	-0x0112, 0x0044,-0x0100,-0x0056, 0x0094,-0x004E,
	-0x0005,-0x00FC, 0x016E,-0x00E8,
	 0x05DF, 0x08A3, 0x1190, 0x0F35, 0x0A5E, 0x09DD,
	 0x1CF0, 0x0A22,
	
	 
	//18
	 0x0A3B, 0x0CBA, 0x0CE5, 0x17FF, 0x1110, 0x06A5,
	 0x171A, 0x11BB,
	-0x0122, 0x0078,-0x0173, 0x0032, 0x0199,-0x0135,
	 0x019F,-0x016B,
	 0x0131,-0x00E8, 0x00D9,-0x02F9,-0x01C7, 0x0178,
	-0x0067,-0x0165,
	
	-0x05A5,-0x0256,-0x0213,-0x09D8, 0x03D5, 0x0853,
	-0x015D,-0x02BD, 0x06BA,-0x0633,
	 0x0D94, 0x0581, 0x14CA, 0x1619, 0x0B8F, 0x0181,
	 0x0F0C, 0x1749,
	 0x00F2,-0x0A18, 0x0121,-0x0A8C, 0x0221,-0x0304,
	-0x04FD,-0x0EE1,
	
	 0x031F,-0x0190,-0x0071,-0x056F, 0x017A, 0x02B7,
	-0x0345,-0x01C7,-0x00CA,-0x059E,
	 0x0450,-0x0214, 0x0046, 0x00DF,-0x0041,-0x0059,
	-0x0139, 0x01F2,-0x0332, 0x0273,
	 0x11A1, 0x09E2, 0x0D6E, 0x1C40, 0x12A3, 0x0548,
	 0x0E23, 0x17CF,
	
	 
	//19
	 0x0B53, 0x0798, 0x1761, 0x094A, 0x1077, 0x03C7,
	 0x0C0E, 0x13D9,
	-0x0055, 0x00AB, 0x0100,-0x00DF, 0x0052, 0x0042,
	-0x00EB, 0x00B4,
	 0x02FD,-0x01FB,-0x02BB,-0x06D3, 0x0243,-0x0167,
	-0x02BB,-0x04B7,
	
	-0x0285, 0x0130,-0x03C4, 0x02B4, 0x007A,-0x012B,
	 0x0112,-0x048B, 0x02F7, 0x0358,
	 0x072E,-0x01C7, 0x11CB, 0x1387, 0x04ED,-0x007E,
	 0x217A, 0x09FF,
	-0x0342,-0x04BE, 0x06AE,-0x0830,-0x0400, 0x0211,
	-0x0154,-0x0A06,
	
	 0x04F5,-0x0301,-0x05B4,-0x0428, 0x0569,-0x027A,
	 0x00FF,-0x0962,-0x0608,-0x0920,
	 0x0067,-0x01A0, 0x03F7,-0x0010,-0x01B0, 0x01DF,
	 0x0402,-0x00A7,-0x06D4, 0x05ED,
	 0x0A29, 0x07B9, 0x1461, 0x1378, 0x0D51, 0x0819,
	 0x13B6, 0x1076,
	
	 
	//20
	 0x1168, 0x0029, 0x0B51, 0x0D02, 0x0CB7, 0x02EF,
	 0x1655, 0x0B70,
	 0x001F, 0x0098,-0x005D, 0x0013, 0x0047, 0x002D,
	-0x000A,-0x0045,
	-0x01A1, 0x008D, 0x0036,-0x018F, 0x00A4, 0x00E7,
	 0x00FF,-0x0191,
	
	 0x06BA, 0x0435,-0x0311, 0x0446, 0x0663, 0x0091,
	 0x04F9,-0x0466,-0x1015,-0x14CB,
	 0x0FFA, 0x04D0, 0x0E35, 0x148A, 0x137F,-0x05D2,
	 0x0E5A, 0x1578,
	 0x00A2,-0x03B8,-0x046B,-0x03A6, 0x00F6,-0x0170,
	 0x02E6,-0x0814,
	
	-0x00EA,-0x03D4, 0x01D4,-0x08CB, 0x0193, 0x022B,
	 0x008F,-0x0808,-0x01E6,-0x12F2,
	-0x0195, 0x0020,-0x00F5,-0x01E1, 0x008E,-0x0167,
	 0x003F,-0x02B3, 0x01CF,-0x0190,
	 0x0C41,-0x0186, 0x0F12, 0x0975, 0x0CE6, 0x0475,
	 0x1793, 0x08E9,
	
	 
	//21
	 0x1645,-0x0135, 0x102A, 0x0C49, 0x076C, 0x0B9C,
	 0x10C4, 0x1000,
	-0x00DE, 0x00B9, 0x00C3,-0x00DA, 0x00F4,-0x0074,
	-0x00D7, 0x00CD,
	-0x0373, 0x00B7, 0x0510,-0x061E,-0x02B3, 0x01C5,
	 0x0024,-0x0359,
	
	 0x033C,-0x0185, 0x04B3,-0x04D1,-0x01A0, 0x04F1,
	-0x0466, 0x0135,-0x03D2, 0x0185,
	 0x080D, 0x00B2, 0x1149, 0x0987, 0x0DAE,-0x00DF,
	 0x1839, 0x18D9,
	 0x0103,-0x03D6, 0x0134,-0x03D9, 0x03F4,-0x0498,
	-0x0718,-0x0386,
	
	-0x0550, 0x032F,-0x00CC,-0x080E, 0x06B7,-0x02BC,
	 0x041F,-0x112B,-0x062B,-0x0ACA,
	 0x0099,-0x0162, 0x03ED,-0x0229,-0x0174, 0x00EC,
	-0x0186, 0x0194,-0x01D4, 0x017C,
	 0x08DA, 0x0BDC, 0x0F60, 0x11F8, 0x1522, 0x0127,
	 0x1177, 0x17EC,
	
	 
	//22
	 0x1214, 0x0780, 0x04F3, 0x1C4D, 0x13C2, 0x0300,
	 0x1510, 0x105C,
	 0x01B5,-0x0045,-0x024C, 0x01B1, 0x012E,-0x00EA,
	-0x00AA, 0x00F5,
	 0x02AE, 0x014D, 0x036E,-0x031C, 0x005E, 0x0333,
	-0x0657, 0x05C2,
	
	 0x01D6,-0x00A7, 0x0256,-0x05B1, 0x00C4, 0x04BD,
	-0x0840,-0x0752, 0x03DC,-0x0D15,
	 0x139A, 0x03FD, 0x0D7C, 0x1ACC, 0x1001,-0x03F7,
	 0x0EA7, 0x17F1,
	 0x01DE,-0x015C,-0x0394,-0x0445, 0x0520,-0x00CB,
	-0x029B,-0x0100,
	
	 0x00B6,-0x0127,-0x017C,-0x0347,-0x01AE, 0x0460,
	 0x02B9,-0x06FC,-0x005F,-0x07C6,
	 0x0182,-0x0147, 0x0035, 0x0037,-0x014B, 0x003A,
	-0x009D, 0x00B5, 0x0013, 0x00BA,
	 0x08DF, 0x09E6, 0x128C, 0x194B, 0x0BE9, 0x0842,
	 0x1861, 0x1128,
	
	 
	//23
	 0x106B, 0x0455, 0x126F, 0x1127, 0x1056, 0x02CE,
	 0x0CE1, 0x0EB9,
	-0x0180, 0x0112,-0x00AB,-0x0027, 0x015D,-0x00CA,
	 0x00DC,-0x01CD,
	 0x03D1,-0x006C,-0x0460,-0x015C,-0x0154, 0x00DC,
	 0x01D6,-0x03A5,
	
	-0x00B5,-0x002E, 0x00D0, 0x000F,-0x008A, 0x0051,
	 0x0006,-0x00E7, 0x006F, 0x000E,
	 0x0EC7, 0x044F, 0x0B7A, 0x1474, 0x14AA,-0x03B1,
	 0x1119, 0x0F2C,
	-0x04F1, 0x015D,-0x0396, 0x0231, 0x0458, 0x0153,
	 0x0411,-0x01AD,
	
	-0x0381, 0x00F1, 0x02D0,-0x07D1, 0x0337,-0x01B1,
	-0x01C4,-0x02AE,-0x00CE,-0x0737,
	 0x00BC, 0x0033,-0x0122, 0x014B, 0x011B,-0x00C9,
	 0x0082, 0x0075,-0x0139, 0x0351,
	 0x0C68, 0x0179, 0x103F, 0x1109, 0x10FB, 0x03CF,
	 0x1255, 0x0A57,
	
	 
	//24
	 0x11A4, 0x0626, 0x0491, 0x15AB, 0x17C2,-0x02F5,
	 0x11C7, 0x0D86,
	 0x02DF,-0x00DA, 0x003E, 0x0020,-0x004E,-0x0003,
	-0x02D8, 0x01FF,
	 0x0536,-0x021E, 0x00B5, 0x0084,-0x01E5, 0x01D1,
	-0x0440, 0x0253,
	
	 0x08C3,-0x07B7,-0x02BE,-0x0661, 0x003A, 0x040A,
	 0x024D,-0x102E,-0x059F,-0x0C81,
	 0x0CB1, 0x057E, 0x0FB3, 0x19D0, 0x104C,-0x045D,
	 0x137B, 0x16B5,
	 0x08A2, 0x02CA,-0x0440, 0x0099,-0x039A,-0x020D,
	-0x0048,-0x036B,
	
	 0x009F,-0x0304, 0x02F1,-0x086E,-0x01FD, 0x02CE,
	-0x023D,-0x05BF, 0x0016,-0x0A5E,
	 0x0125,-0x0111, 0x006F, 0x005E,-0x003E,-0x0017,
	 0x014C,-0x007C,-0x028F, 0x01FE,
	 0x0B04, 0x065F, 0x0E25, 0x138B, 0x0F4B, 0x02FF,
	 0x16C1, 0x0DF6,
#endif
};
typedef struct T43NodeStruct
{
	int n[2];
#ifndef T43_DISABLE_REC
	unsigned short rec[T43_N_REC_ESTIMATORS];
#endif
} T43Node;
typedef struct T43CtxStruct
{
	int pred[T43_NFILTERS];
	int context[T43_NMAPS];
	ArrayHandle maps[24][T43_NMAPS];//(256+512+1024+2048+4096+8192+16384+32768)*20*14 = 17.43 MB for 14 maps with 6 rec estimators
	T43Node *node[T43_NMAPS];

	int p0arr[T43_NESTIMATORS], p0_0, p0;//p0_0 isn't clamped
	int weights[24][T43_NESTIMATORS];
	long long wsum;

	int nnodes;
#ifdef T43_PRINT_ESTIMATOR_CR
	float csizes_est[24*T43_NESTIMATORS];
#endif
} T43Ctx;
T43Ctx* t43_ctx_init()
{
	int val=0x8000;
	T43Node node0={{1, 1}};
#ifndef T43_DISABLE_REC
	for(int k=0;k<T43_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	T43Ctx *ctx=(T43Ctx*)malloc(sizeof(T43Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T43Ctx));
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T43_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T43Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T43Node), sizeof(T43Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t43_ctx_clear(T43Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T43_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
#if 0
static int custom3_loadnb(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-T43_REACH;ky2<0;++ky2)
	{
		for(int kx2=-T43_REACH;kx2<=T43_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-T43_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static int custom3_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-15;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((__m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((__m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
	__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
	s2=_mm_add_epi32(s2, hi);
	s2=_mm_hadd_epi32(s2, s2);
	int s3=_mm_extract_epi32(s2, 0);
	for(;k<count;++k)
		s3+=a[k]*b[k];
	return s3;
}
#endif
void t43_ctx_get_context(T43Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
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

		m1  =LOAD(buf, 1, 0, 0),
		Nm1 =LOAD(buf, 1, 0, 1),
		Wm1 =LOAD(buf, 1, 1, 0),
		NWm1=LOAD(buf, 1, 1, 1),

		m2  =LOAD(buf, 2, 0, 0),
		Nm2 =LOAD(buf, 2, 0, 1),
		Wm2 =LOAD(buf, 2, 1, 0),
		NWm2=LOAD(buf, 2, 1, 1);

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

	for(int k=0;k<T43_NFILTERS;++k)
	{//CUSTOM3
		T43C3Params *params=(T43C3Params*)t43_allparams+k;
		const short *coeffs[]=
		{
			params->c00, params->c01, params->c02,
			params->c10, params->c11, params->c12,
			params->c20, params->c21, params->c22,

			//0, 0, 0,
			//0, 0, 0,
			//0, 0, 0,
		};
		short nb[3][T43_NNB+2]={0};
		int count[3], idx, idx2;
		for(int kc=0;kc<3;++kc)
			count[kc]=custom3_loadnb(buf, ebuf, iw, ih, kc, kx, ky, nb[kc]);
		idx=(iw*ky+kx)<<2;
		idx2=0;
		switch(kc)
		{
		case 1:
			nb[0][T43_NNB  ]=buf[idx];
			nb[0][T43_NNB+1]=ebuf[idx];
			count[0]+=2;
			//++idx;
			idx2+=3;
			break;
		case 2:
			nb[0][T43_NNB  ]=buf [idx  ];
			nb[0][T43_NNB+1]=ebuf[idx  ];
			nb[1][T43_NNB  ]=buf [idx|1];
			nb[2][T43_NNB+1]=ebuf[idx|1];
			count[0]+=2;
			count[1]+=2;
			//idx+=2;
			idx2+=6;
			break;
		}
		ctx->pred[k]=0;
		for(int kc=0;kc<3;++kc)
		{
			if(coeffs[idx2+kc])
				ctx->pred[k]+=custom3_dot(coeffs[idx2+kc], nb[kc], count[kc]);
		}

		ctx->pred[k]+=1<<13;
		ctx->pred[k]>>=14;
		ctx->pred[k]=CLAMP(-128, ctx->pred[k], 127);
		ctx->context[++j]=ctx->pred[k];
	}
#undef LOAD
	for(int k=0;k<T43_NMAPS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
void t43_ctx_estimate_p0(T43Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	T43Node *node;
	for(int kp=0;kp<T43_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		int context=ctx->context[kp];
		ArrayHandle map=ctx->maps[workidx][kp];
		node=ctx->node[kp]=(T43Node*)array_at(&map, context);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T43_DISABLE_REC
		for(;k2<T43_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T43_NESTIMATORS;++k)
	{
#ifdef T43_DISABLE_COUNTER
		if(k%(T43_N_REC_ESTIMATORS+1))//
#endif
		{
			sum+=(long long)ctx->p0arr[k]*wk[k];
			ctx->wsum+=wk[k];
		}
	}
	ctx->p0=ctx->wsum?(int)(sum/ctx->wsum):0x8000;
	ctx->p0_0=ctx->p0;

	ctx->p0=CLAMP(1, ctx->p0, 0xFFFF);
}
void t43_ctx_update(T43Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;
	
#ifdef T43_PRINT_ESTIMATOR_CR
	for(int k=0;k<T43_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T43_NESTIMATORS*workidx+k]+=bitsize;
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
		for(int k=0;k<T43_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=T43_LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
	T43Node *node;
	for(int kp=0;kp<T43_NMAPS;++kp)
	{
		node=ctx->node[kp];
		++node->n[bit];
#ifndef T43_DISABLE_REC
		for(int k=0;k<T43_N_REC_ESTIMATORS;++k)
		{
			int lgden=k;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);
	}
}
void t43_update_error(T43Ctx *ctx, char pixel, char *error)
{
	int pred=0;
	for(int k=0;k<T43_NFILTERS;++k)
		pred+=ctx->pred[k];
	pred+=T43_NFILTERS/2;
	pred/=T43_NFILTERS;
	*error=pixel-pred;
}
int t43_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T43: Wisdom of the crowd   Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);
	T43Ctx *t43_ctx=t43_ctx_init();
	if(!buf2||!ebuf||!t43_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	memset(ebuf, 0, (size_t)res<<2);
	addbuf((unsigned char*)buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd(buf2, iw, ih);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	ABACEncContext ctx;
	abac_enc_init(&ctx, &list);
	
	float csizes[24]={0};
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t43_ctx_get_context(t43_ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				
				if(kx==10&&ky==10)//
					printf("");

				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t43_ctx_estimate_p0(t43_ctx, kc, kb);
					int bit=(buf2[idx]+128)>>kb&1;
					abac_enc(&ctx, t43_ctx->p0, bit);
					
					int prob=bit?0x10000-t43_ctx->p0:t43_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t43_ctx_update(t43_ctx, kc, kb, bit);
				}
				t43_update_error(t43_ctx, buf2[idx], ebuf+idx);
			}
		}
		if(loud)
		{
			static float csize_prev=0;
			float csize=0;
			for(int k=0;k<24;++k)
				csize+=csizes[k]/8;
			printf("%5d/%5d  %6.2lf%%  CR%11f  CR_delta%11f\r", ky+1, ih, 100.*(ky+1)/ih, iw*(ky+1)*3/csize, iw*3/(csize-csize_prev));
			csize_prev=csize;
		}
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t43_ctx->nnodes*sizeof(T43Node)/(1024*1024));
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

#ifdef T43_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t43_ctx->csizes_est+T43_NESTIMATORS*kb;
				for(int ke=1;ke<T43_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T43_NESTIMATORS;++ke)
			{
				float *sizes=t43_ctx->csizes_est+ke;
#ifndef T43_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T43_N_REC_ESTIMATORS+1), ke%(T43_N_REC_ESTIMATORS+1));
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
					printf("%8.2f %c", iw*ih/sizes[T43_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T43_NESTIMATORS*kb]/t43_ctx->csizes_est[T43_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T43_DISABLE_REC
				if(!((ke+1)%(T43_N_REC_ESTIMATORS+1)))
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
	t43_ctx_clear(&t43_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t43_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();

	char *ebuf=(char*)malloc((size_t)res<<2);
	T43Ctx *t43_ctx=t43_ctx_init();
	if(!ebuf||!t43_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t43_ctx_init(t43_ctx);
	
	for(int ky=0, idx;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				idx=(iw*ky+kx)<<2|kc;
				t43_ctx_get_context(t43_ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t43_ctx_estimate_p0(t43_ctx, kc, kb);
					
					int bit=abac_dec(&ctx, t43_ctx->p0);
					buf[idx]|=bit<<kb;

					t43_ctx_update(t43_ctx, kc, kb, bit);
				}
				buf[idx]+=128;//unsigned -> signed
				t43_update_error(t43_ctx, buf[idx], ebuf+idx);
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
	}
	t43_ctx_clear(&t43_ctx);
	
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	return 1;
}
#endif


void ac_vs_ans()
{
	DList list_abac, list_bans;
	dlist_init(&list_abac, 1, 1024, 0);
	dlist_init(&list_bans, 1, 1024, 0);

	ABACEncContext abac;
	abac_enc_init(&abac, &list_abac);

	BANSEncContext bans;
	bans_enc_init(&bans, &list_bans);

	int p0=0xFFFF;
	
	for(int k=0;k<1000000;++k)
	{
		p0=rand()|1;
		abac_enc(&abac, p0, 0);
		bans_enc(&bans, p0, 0);
	}

	//double t_abac=time_ms();
	//for(int k=0;k<1000000000;++k)
	//	abac_enc(&abac, p0, 0);
	//abac_enc_flush(&abac);
	//t_abac=time_ms()-t_abac;
	//
	//double t_bans=time_ms();
	//for(int k=0;k<1000000000;++k)
	//	bans_enc(&bans, p0, 0);
	//bans_enc_flush(&bans);
	//t_bans=time_ms()-t_bans;

	//printf("ABAC  %8lld bytes  %lf ms\n", list_abac.nobj, t_abac);//5724 bytes  3160.941600 ms
	//printf("BANS  %8lld bytes  %lf ms\n", list_bans.nobj, t_bans);//2618 bytes  8928.841800 ms
	printf("ABAC  %8lld bytes\n", list_abac.nobj);
	printf("BANS  %8lld bytes\n", list_bans.nobj);
	dlist_clear(&list_abac);
	dlist_clear(&list_bans);
	pause();
	exit(0);
}