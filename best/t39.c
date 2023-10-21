#include"best.h"
#define AC_IMPLEMENTATION
#include"ac.h"
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


//T39: Multiple estimators for all maps

//	#define T39_DISABLE_COUNTER
//	#define T39_DISABLE_REC
//	#define T39_PRINT_ESTIMATOR_CR

#define T39_LR (int)(0.07*0x10000+0.5)
#define T39_NMAPS 15	//14	31		135 HALF HOUR PER IMAGE

#ifndef T39_DISABLE_REC
#define T39_N_REC_ESTIMATORS 6		//15
#define T39_NESTIMATORS ((T39_N_REC_ESTIMATORS+1)*T39_NMAPS)
#else
#define T39_NESTIMATORS T39_NMAPS
#endif
typedef struct T39NodeStruct
{
	int n[2];
#ifndef T39_DISABLE_REC
	unsigned short rec[T39_N_REC_ESTIMATORS];
#endif
} T39Node;
typedef struct T39CtxStruct
{
	int pred14;
	int context[T39_NMAPS];
	ArrayHandle maps[24][T39_NMAPS];//3*(256+512+1024+2048+4096+8192+16384+32768)*15*sizeof(T39Node) = 56.03 MB for 14 maps with 6 rec estimators
	T39Node *node[T39_NMAPS];

	int p0arr[T39_NESTIMATORS], p0_0, p0;//p0_0 isn't clamped
	int weights[24][T39_NESTIMATORS];
	long long wsum;

	int nnodes;
#ifdef T39_PRINT_ESTIMATOR_CR
	float csizes_est[24*T39_NESTIMATORS];
#endif
} T39Ctx;
#if 0
void t39_explore(T39Ctx *ctx, const char *src, int iw, int ih)
{
	int res=iw*ih;
	int hist[256]={0};
	int kc=0,//red (orangeness)
		ke=0;//zero predictor
	//double entropy=0;
	for(int sym=0;sym<256;++sym)
	{
		if(sym==128)//
			printf("");

		int prob=0x10000;
		int context=0x80;
		for(int kb=7;kb>=0;--kb)
		{
			ArrayHandle map=ctx->maps[kc<<3|kb][ke];
			int bit=sym>>kb&1;
			T39Node *node=array_at(&map, context);

			//int pb=bit?0x10000-node->rec[3]:node->rec[3];
			int pb=(int)(((long long)node->n[bit]<<16)/(node->n[0]+node->n[1]));

			pb=CLAMP(1, pb, 0xFFFF);
			prob=(int)(((long long)prob*pb+0x8000)>>16);
			context|=bit<<(8+7-kb);
		}
		hist[sym]=prob;

		//printf("%3d  0x%04X\n", sym, prob);
		//if(prob)
		//{
		//	double p=(double)prob/res;
		//	entropy-=p*log2(p);				//X  need to use cross-entropy with image histogram
		//}
	}
	//double invCR=entropy/8;
	//printf("CR %lf\n", 1/invCR);

	double csize=0;
	for(int k=0;k<res;++k)//Zipf's law
	{
		unsigned char sym=src[k<<2|kc]+128;
		int prob=hist[sym];
		if(prob)
		{
			double p=(double)prob/res;
			csize-=log2(p);
		}
	}
	csize/=8;
	printf("C%d  csize %lf  CR %lf\n", kc, csize, res/csize);
}
#endif
T39Ctx* t39_ctx_init()
{
	int val=0x8000;
	T39Node node0={{1, 1}};
#ifndef T39_DISABLE_REC
	for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	T39Ctx *ctx=(T39Ctx*)malloc(sizeof(T39Ctx));
	if(!ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(ctx, 0, sizeof(T39Ctx));
	memfill(ctx->weights, &val, sizeof(ctx->weights), sizeof(int));
	for(int k=0;k<24;++k)
	{
		int kb=k&7;
		for(int k2=0;k2<T39_NMAPS;++k2)
		{
			int nnodes=256<<(7-kb);
			ARRAY_ALLOC(T39Node, ctx->maps[k][k2], 0, nnodes, 0, 0);
			memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T39Node), sizeof(T39Node));
			ctx->nnodes+=nnodes;
		}
	}
	return ctx;
}
void t39_ctx_clear(T39Ctx **ctx)
{
	for(int k=0;k<24;++k)
	{
		for(int k2=0;k2<T39_NMAPS;++k2)
			array_free(ctx[0]->maps[k]+k2);
	}
	free(*ctx);
	*ctx=0;
}
#if 0
void t39_ctx_reset(T39Ctx *ctx, int hardreset)
{
	T39Node node0={{1, 1}};
#ifndef T39_DISABLE_REC
	for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
		node0.rec[k]=0x8000;
#endif
	for(int k=0;k<24;++k)
	{
		if(hardreset)
		{
			for(int k2=0;k2<T39_NMAPS;++k2)
				memfill(ctx->maps[k][k2]->data, &node0, ctx->maps[k][k2]->count*sizeof(T39Node), sizeof(T39Node));
		}
	}
	if(hardreset)
		ctx->nnodes=0;
}
#endif
void t39_ctx_get_context(T39Ctx *ctx, const char *buf, const char *ebuf, int iw, int ih, int kc, int kx, int ky)
{
#define LOAD(BUF, C, X, Y) (unsigned)(kc-C)<3&&(unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)]:0
	int count_W_N_m1=(kx-1>=0)+(ky-1>=0)+(kc-1>=0);
	char
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
	//ctx->context[++j]=(N+W+NW+NE)>>2;//10_v2
	ctx->context[++j]=(N+W-NW + m1)>>1;//11
	ctx->context[++j]=m2;//12
	ctx->context[++j]=(N+W-NW + m2)>>1;//13

	//kodim13
#if 1
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
#endif//kodim13
	
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
#endif//CLIC16
	
#undef LOAD
	for(int k=0;k<T39_NMAPS;++k)
	{
		ctx->context[k]+=128;
		ctx->context[k]=CLAMP(0, ctx->context[k], 255);
	}
}
int t39_ctx_map_context(int *context, int kp, int workidx)//replacement for context[kp]
{
	return context[kp];

	//static const int rep[]={ 0,  0, 10, 10, 10, 10, 10,  0,     10, 10, 10,  0,  0, 10, 10,  0,      8,  8, 11, 11, 11, 11, 11,  0};
	//static const int sub[]={ 6,  6,  6,  6,  2,  0,  0, 10,     10,  6,  6,  2,  2, 13,  6, 10,      2,  1,  2, 13, 13, 12,  0, 12};
	//return context[kp==rep[workidx]?sub[workidx]:kp];
}
void t39_ctx_estimate_p0(T39Ctx *ctx, int kc, int kb)
{
	int workidx=kc<<3|kb;
	int *wk=ctx->weights[workidx];

	int p0idx=0;
	long long sum;
	T39Node *node;
	for(int kp=0;kp<T39_NMAPS;++kp)//for each predictor
	{
		int k2=0;
		int context=t39_ctx_map_context(ctx->context, kp, kc);
		ArrayHandle map=ctx->maps[workidx][kp];
		node=ctx->node[kp]=(T39Node*)array_at(&map, context);
		
		sum=node->n[0]+node->n[1];
		ctx->p0arr[p0idx+k2]=sum?(int)(((long long)node->n[0]<<16)/sum):0x8000;
		++k2;
#ifndef T39_DISABLE_REC
		for(;k2<T39_N_REC_ESTIMATORS+1;++k2)
			ctx->p0arr[p0idx+k2]=node->rec[k2-1];
#endif
		p0idx+=k2;
	}

	sum=0;
	ctx->wsum=0;
	for(int k=0;k<T39_NESTIMATORS;++k)
	{
#ifdef T39_DISABLE_COUNTER
		if(k%(T39_N_REC_ESTIMATORS+1))//
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
void t39_ctx_update(T39Ctx *ctx, int kc, int kb, int bit)
{
	int workidx=kc<<3|kb;

#ifdef T39_PRINT_ESTIMATOR_CR
	for(int k=0;k<T39_NESTIMATORS;++k)
	{
		int prob=(bit?0x10000-ctx->p0arr[k]:ctx->p0arr[k]);
		if(prob)
		{
			float p=(float)prob/0x10000;
			float bitsize=-log2f(p);
			ctx->csizes_est[T39_NESTIMATORS*workidx+k]+=bitsize;
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
		for(int k=0;k<T39_NESTIMATORS;++k)
		{
			int diff=ctx->p0arr[k]-ctx->p0;//fixed 15.16 bit
			long long grad = dL_dp0*diff/ctx->wsum;
			long long wnew=T39_LR*grad>>16;
			wnew=wk[k]-wnew;
			wnew=CLAMP(1, wnew, 0xFFFF);
			wk[k]=(int)wnew;
		}
	}

	//update
	T39Node *node;
	for(int kp=0;kp<T39_NMAPS;++kp)
	{
		node=ctx->node[kp];
			++node->n[bit];
#ifndef T39_DISABLE_REC
		for(int k=0;k<T39_N_REC_ESTIMATORS;++k)
		{
			int lgden=k;
			int temp=node->rec[k]+(((!bit<<16)-node->rec[k])>>lgden);
			node->rec[k]=CLAMP(1, temp, 0xFFFF);
		}
#endif
		ctx->context[kp]|=bit<<(8+7-kb);//append bits in reverse
	}
}
int t39_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T39 Enc  Multiple estimators for all maps with CUSTOM1 filter  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	char *buf2=(char*)malloc((size_t)res<<2);
	char *ebuf=(char*)malloc((size_t)res<<2);
	T39Ctx *t39_ctx=t39_ctx_init();
	if(!buf2||!ebuf||!t39_ctx)
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
	//int hits[24]={0};
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t39_ctx_get_context(t39_ctx, (char*)buf2, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t39_ctx_estimate_p0(t39_ctx, kc, kb);
					int bit=(buf2[idx]+128)>>kb&1;
					abac_enc(&ctx, t39_ctx->p0, bit);
					
					int prob=bit?0x10000-t39_ctx->p0:t39_ctx->p0;//
					float bitsize=-log2f((float)prob*(1.f/0x10000));
					csizes[kc<<3|kb]+=bitsize;//

					t39_ctx_update(t39_ctx, kc, kb, bit);
				}
				ebuf[idx]=buf2[idx]-t39_ctx->pred14;
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
	}
	abac_enc_flush(&ctx);

	size_t dststart=dlist_appendtoarray(&list, data);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Used %f MB of memory\n", (float)t39_ctx->nnodes*sizeof(T39Node)/(1024*1024));
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
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

#ifdef T39_PRINT_ESTIMATOR_CR
		if(loud==2)
		{
			printf("Estimator efficiencies:\n");
			int minidx[24]={0}, maxidx[24]={0};
			for(int kb=0;kb<24;++kb)
			{
				float *sizes=t39_ctx->csizes_est+T39_NESTIMATORS*kb;
				for(int ke=1;ke<T39_NESTIMATORS;++ke)
				{
					if(sizes[minidx[kb]]>sizes[ke])
						minidx[kb]=ke;
					if(sizes[maxidx[kb]]<sizes[ke])
						maxidx[kb]=ke;
				}
			}
			for(int ke=0;ke<T39_NESTIMATORS;++ke)
			{
				float *sizes=t39_ctx->csizes_est+ke;
#ifndef T39_DISABLE_REC
				printf("E%3d-%02d-%02d ", ke, ke/(T39_N_REC_ESTIMATORS+1), ke%(T39_N_REC_ESTIMATORS+1));
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
					printf("%8.2f %c", iw*ih/sizes[T39_NESTIMATORS*kb], c);
					//printf(" %7.2f%c", sizes[T39_NESTIMATORS*kb]/t39_ctx->csizes_est[T39_NESTIMATORS*kb+minidx[kb]], c);
					if(kb+1<24&&!((kb+1)&7))
						printf("    ");
				}
				printf("\n");
#ifndef T39_DISABLE_REC
				if(!((ke+1)%(T39_N_REC_ESTIMATORS+1)))
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
		//t39_explore(t39_ctx, buf2, iw, ih);
	}
	t39_ctx_clear(&t39_ctx);
	dlist_clear(&list);
	free(buf2);
	return 1;
}
int t39_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T39 Dec  %s\n", g_buf);
	}

	//int debug_index=0;
	char *ebuf=(char*)malloc((size_t)res<<2);
	T39Ctx *t39_ctx=t39_ctx_init();
	if(!ebuf||!t39_ctx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	ABACDecContext ctx;
	abac_dec_init(&ctx, data, data+srclen);

	int black=0xFF000000;
	memfill(buf, &black, res*sizeof(int), sizeof(int));
	t39_ctx_init(t39_ctx);
	
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=(iw*ky+kx)<<2|kc;
				t39_ctx_get_context(t39_ctx, (char*)buf, ebuf, iw, ih, kc, kx, ky);
				for(int kb=7;kb>=0;--kb)//MSB -> LSB
				{
					t39_ctx_estimate_p0(t39_ctx, kc, kb);
					
					int bit=abac_dec(&ctx, t39_ctx->p0);
					buf[idx]|=bit<<kb;

					t39_ctx_update(t39_ctx, kc, kb, bit);
				}
				buf[idx]+=128;//unsigned -> signed
				ebuf[idx]=buf[idx]-t39_ctx->pred14;
			}
		}
		if(loud)
			printf("%5d/%5d  %6.2lf%%\r", ky+1, ih, 100.*(ky+1)/ih);
	}
	t39_ctx_clear(&t39_ctx);
	
	//addbuf(buf, iw, ih, 3, 4, 128);//X  the buffer is signed
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}