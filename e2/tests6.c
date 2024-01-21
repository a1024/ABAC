#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

	#define SLIC5_DISABLE_PALETTE

#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0

static const int qlevels[]=
{
	1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47, 55, 63, 79, 95, 111, 127
	//1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392, 500//from libjxl
};

typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
static void hybriduint_encode(unsigned val, HybridUint *hu)
{
	int token, bypass, nbits;
	if(val<(1<<SLIC5_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC5_CONFIG_EXP) + (
				(lgv-SLIC5_CONFIG_EXP)<<(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC5_CONFIG_MSB))<<SLIC5_CONFIG_LSB|
				(mantissa&((1<<SLIC5_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB);
		bypass=val>>SLIC5_CONFIG_LSB&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}
static int quantize(int x)
{
	int q;
	if(x<qlevels[8])
	{
		if(x<qlevels[4])
		{
			if(x<qlevels[2])
				q=x<qlevels[1]?x>=qlevels[0]:2;
			else
				q=x<qlevels[3]?3:4;
		}
		else
		{
			if(x<qlevels[6])
				q=x<qlevels[5]?5:6;
			else
				q=x<qlevels[7]?7:8;
		}
	}
	else
	{
		if(x<qlevels[12])
		{
			if(x<qlevels[10])
				q=x<qlevels[9]?9:10;
			else
				q=x<qlevels[11]?11:12;
		}
		else
		{
			if(x<qlevels[14])
				q=x<qlevels[13]?13:14;
			else
				q=x<qlevels[15]?15:16;
		}
	}
	return q;
}
static int quantize_signed(int x)
{
	int x2=quantize(abs(x));
	//int x2=quantize(abs(x))>>1;
	int neg=-(x<0);
	x2^=neg;
	x2-=neg;
	x2<<=1;
	x2^=neg;
	//x=quantize(abs(x))<<1|(x<0);
	//if(x2>_countof(qlevels))
	//	x2=_countof(qlevels);
	return x2;
}
#define PAD_SIZE 2
#define NPREDS 10
#define PRED_PREC 8
#define PARAM_PREC 8
#define SSE_STAGES 8
#define NHIST (_countof(qlevels)+1)
#define SSE_SIZE (NHIST*NHIST<<2)
//#define HASH_FUNC(A, B) ((A+B+1)*(A+B)/2+B)
static const int wp_params[]=
{
	// 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x0004,
	// 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,


	0x0004AB65,//wp1
	0x0002841C,//wp2
	0x00034147,//wp3
	0x00016056,//wp4
	0x00016056,//wp5
	0x00016056,//wp6
	0x00016056,//wp7
	0x00016056,//clamp grad
	0x00028F8A,//journal GAP
	0x00022F58,//CALIC GAP


//	0xD000, 0xC000, 0xC000, 0xB000,//wp weights
//	0xC000,//journal GAP
//	0xD000,//CALIC GAP

//	 0x0DB8,  0x0E22,  0x181F,  0x0BF3,//wp weights
//	//0x181F,//clamp grad
//	 0x181F,//journal GAP
//	 0x1A1F,//CALIC GAP

	-0x005C,
	-0x005B,
	 0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,

	//0xd, 0xc, 0xc, 0xb,
	//8,
	//8,
	//4, 0, 3, 23, 2,
};
typedef struct SLIC5CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int nlevels[4], min_allowed[4], max_allowed[4];
	int cdfsize;
	int *pred_errors, *errors, *pixels;
	long long *hist, histsums[NHIST];
	unsigned *CDFs;
	long long esum[4];//sum of abs errors so far in current row
	int shift[4];//depth compensation
	long long *sse, sse_idx[SSE_STAGES], sse_sum[SSE_STAGES];
	int sse_count[SSE_STAGES];
	int preds[NPREDS], params[_countof(wp_params)];
	long long pred;
	int pred_final;
	int hist_idx;
	int sse_corr;
	int kc, kx, ky, kym0, kym1, kym2;
	ArithmeticCoder *ec;
} SLIC5Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx-(X)+PAD_SIZE)<<2|kc]
static int slic5_init(SLIC5Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)
{
	if(iw<1||ih<1||nch<1)
	{
		LOG_ERROR("Invalid image");
		return 0;
	}
	memset(pr, 0, sizeof(*pr));
	pr->iw=iw;
	pr->ih=ih;
	pr->nch=nch;
	if(nch<3)
		memcpy(pr->depths, depths, nch);
	else
	{
		pr->depths[0]=depths[1];  //Y
		pr->depths[1]=depths[2]+1;//Cb
		pr->depths[2]=depths[0]+1;//Cr
		if(nch==4)
			pr->depths[3]=depths[3];//a
	}
	int maxdepth=0;
	for(int kc=0;kc<nch;++kc)
	{
		pr->nlevels[kc]=1<<pr->depths[kc];
		pr->min_allowed[kc]=-(pr->nlevels[kc]>>1);
		pr->max_allowed[kc]=(pr->nlevels[kc]>>1)-1;
		if(maxdepth<pr->depths[kc])
			maxdepth=pr->depths[kc];
	}
	HybridUint hu;
	hybriduint_encode(-((1<<maxdepth)>>1), &hu);//encode -half
	pr->cdfsize=hu.token+1;

	pr->pred_errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));//NPREDS * 4 rows * 4 comps
	pr->errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->pixels=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->hist=(long long*)malloc(pr->cdfsize*sizeof(long long[NHIST]));//WH: cdfsize * NHIST
	pr->CDFs=(unsigned*)malloc((pr->cdfsize+1LL)*sizeof(unsigned[NHIST]));
	pr->sse=(long long*)malloc(sizeof(long long[SSE_STAGES*SSE_SIZE]));
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->CDFs||!pr->sse)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pred_errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));
	memset(pr->errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	memset(pr->pixels, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	long long fillval=1;
	memfill(pr->hist, &fillval, pr->cdfsize*sizeof(long long[NHIST]), sizeof(fillval));
	fillval=pr->cdfsize;
	memfill(pr->histsums, &fillval, sizeof(pr->histsums), sizeof(fillval));
	memset(pr->CDFs, 0, (pr->cdfsize+1LL)*sizeof(unsigned[NHIST]));
	memset(pr->sse, 0, sizeof(long long[SSE_STAGES*SSE_SIZE]));
	pr->ec=ec;
	for(int k=0;k<_countof(wp_params);++k)
		pr->params[k]=wp_params[k];
	return 1;
}
static void slic5_free(SLIC5Ctx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->pixels);
	free(pr->sse);
	free(pr->hist);
	free(pr->CDFs);
}
static void slic5_update_CDFs(SLIC5Ctx *pr)
{
	for(int kt=0;kt<NHIST;++kt)//update CDFs
	{
		long long *curr_hist=pr->hist+pr->cdfsize*kt;
		unsigned *curr_CDF=pr->CDFs+(pr->cdfsize+1)*kt;
		long long sum=pr->histsums[kt];
		long long c=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			long long freq=curr_hist[ks];
			curr_CDF[ks]=(int)(c*(0x10000LL-pr->cdfsize)/sum)+ks;
			c+=freq;
		}
		curr_CDF[pr->cdfsize]=0x10000;
	}
}
static void slic5_nextrow(SLIC5Ctx *pr, int ky)
{
	if(ky)
	{
		//slic5_update_CDFs(pr);
		for(int kc=0;kc<pr->nch;++kc)
		{
			int x=(int)(pr->esum[kc]/pr->iw);
			x=floor_log2(x)+1-PRED_PREC;
			pr->shift[kc]=MAXVAR(8, x)-8+PRED_PREC;
			pr->esum[kc]=0;
		}
	}
	else
	{
		for(int kc=0;kc<pr->nch;++kc)
			pr->shift[kc]=MAXVAR(8, pr->depths[kc])-8+PRED_PREC;
	}
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic5_predict(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	int idx=(pr->iw*ky+kx)<<2|kc;
	if(idx<0xFFF||!(idx&0xFFF))
	//if(!(idx&idx>>1))
		slic5_update_CDFs(pr);
	//XY are flipped, no need to check if indices OOB due to padding
	int
		NNWW    =LOAD(pr->pixels,  2, 2)<<PRED_PREC,
		NNW     =LOAD(pr->pixels,  1, 2)<<PRED_PREC,
		NN      =LOAD(pr->pixels,  0, 2)<<PRED_PREC,
		NNE     =LOAD(pr->pixels, -1, 2)<<PRED_PREC,
		NNEE    =LOAD(pr->pixels, -2, 2)<<PRED_PREC,
		NWW     =LOAD(pr->pixels,  2, 1)<<PRED_PREC,
		NW      =LOAD(pr->pixels,  1, 1)<<PRED_PREC,
		N       =LOAD(pr->pixels,  0, 1)<<PRED_PREC,
		NE      =LOAD(pr->pixels, -1, 1)<<PRED_PREC,
		NEE     =LOAD(pr->pixels, -2, 1)<<PRED_PREC,
		WW      =LOAD(pr->pixels,  2, 0)<<PRED_PREC,
		W       =LOAD(pr->pixels,  1, 0)<<PRED_PREC;
	
	long long
		eN =LOAD(pr->errors,  0, 1),//error = (curr<<8) - pred
		eW =LOAD(pr->errors,  1, 0),
		eNW=LOAD(pr->errors,  1, 1),
		eNE=LOAD(pr->errors, -1, 0);

	int j=-1;
	pr->preds[++j]=(int)(W+NE-N);
	pr->preds[++j]=(int)(N-((eN+eW+eNE)*pr->params[NPREDS+0]>>PARAM_PREC));
	pr->preds[++j]=(int)(W-((eN+eW+eNW)*pr->params[NPREDS+1]>>PARAM_PREC));
	pr->preds[++j]=(int)(N-((
		eNW*pr->params[NPREDS+2]+
		eN*pr->params[NPREDS+3]+
		eNE*pr->params[NPREDS+4]+
		((long long)NN-N)*pr->params[NPREDS+5]+
		((long long)NW-W)*pr->params[NPREDS+6]
	)>>PARAM_PREC));
	
#if 0
	//pr->preds[++j]=(((N+W)<<1)+NW+NE+NN+WW+4)>>3;//X

	pr->preds[++j]=(int)(W -(eW *pr->params[NPREDS+ 7]>>8));
	pr->preds[++j]=(int)(N -(eN *pr->params[NPREDS+ 8]>>8));
	pr->preds[++j]=(int)(NW-(eNW*pr->params[NPREDS+ 9]>>8));
	pr->preds[++j]=(int)(NE-(eNE*pr->params[NPREDS+10]>>8));

	int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
	int vmin2=MINVAR(vmin, NW), vmax2=MAXVAR(vmax, NW);
	vmin2=MINVAR(vmin2, NE), vmax2=MAXVAR(vmax, NE);
	++j, pr->preds[j]=N+W-NW, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=W+NE-N, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=N+NW-NNW, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	vmin2=MINVAR(vmin, NE), vmax2=MAXVAR(vmax, NE);
	vmin2=MINVAR(vmin2, NEE), vmax2=MAXVAR(vmax, NEE);
	++j, pr->preds[j]=N+NE-NNE, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	
	++j, pr->preds[j]=(W+NEE)/2;
	++j, pr->preds[j]=NNNNNN;
	++j, pr->preds[j]=(NEEEE+NEEEEEE)/2;
	++j, pr->preds[j]=(WWWW+WWWWWW)/2;
	++j, pr->preds[j]=(N+W+NEEEEE+NEEEEEEE)/4;
	++j, pr->preds[j]=N*2-NN, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=(N+NNN)/2;
	++j, pr->preds[j]=((N+W)*3-NW*2)/4;
	++j, pr->preds[j]=(int)(N+W-NW), pr->preds[j]=(int)CLAMP(vmin, pr->preds[j], vmax);
#endif
#if 1
	++j, pr->preds[j]=N>>(abs((int)eN)>=abs(N));//v
	++j, pr->preds[j]=W>>(abs((int)eW)>=abs(W));//v
	++j, pr->preds[j]=(int)(N+W+NW+NE+((eN+eW+eNW+eNE)>>2))>>2;//v
	//++j, pr->preds[j]=(NNWW+NNW+NN+NNE+NNEE+NWW+NW+N+NE+NEE+WW+W)/12;//~
	//++j, pr->preds[j]=(4*(N+W)+2*(NW+NE)+eN+eW+eNW+eNE)>>3;//X
	//++j, pr->preds[j]=2*W-WW+(eW>>1), pr->preds[j]=CLAMP(pr->min_allowed, pr->preds[j], pr->max_allowed);//X

	++j, pr->preds[j]=(int)(N+W-NW);//, pr->preds[j]=(int)MEDIAN3(N, W, pr->preds[j]);
	
	++j;
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int sum2=(dy+dx)>>pr->shift[kc], diff=(dy-dx)>>pr->shift[kc], diff2=NE-NW;
	if(sum2>(32<<(PRED_PREC+1)))
		pr->preds[j]=(dx*N+dy*W)/sum2+diff2/8;
	else if(diff>(12<<(PRED_PREC+1)))
		pr->preds[j]=(N+2*W)/3+diff2/8;
	else if(diff<-(12<<(PRED_PREC+1)))
		pr->preds[j]=(2*N+W)/3+diff2/8;
	else
		pr->preds[j]=(N+W)/2+diff2/8;
	diff=d45-d135;
	if(diff>(32<<(PRED_PREC+1)))
		pr->preds[j]+=diff2/8;
	else if(diff>(16<<(PRED_PREC+1)))
		pr->preds[j]+=diff2/16;
	else if(diff<-(32<<(PRED_PREC+1)))
		pr->preds[j]-=diff2/8;
	else if(diff<-(16<<(PRED_PREC+1)))
		pr->preds[j]-=diff2/16;
	
	++j;
	int disc=(dy-dx)>>pr->shift[kc];
	if(disc>80)
		pr->preds[j]=W;
	else if(disc<-80)
		pr->preds[j]=N;
	else
	{
		pr->preds[j]=(N+W)/2+(NE-NW)/4;
		if(disc>32)
			pr->preds[j]=(pr->preds[j]+W)/2;
		else if(disc>8)
			pr->preds[j]=(3*pr->preds[j]+W)/4;
		else if(disc<-32)
			pr->preds[j]=(pr->preds[j]+N)/2;
		else if(disc<-8)
			pr->preds[j]=(3*pr->preds[j]+N)/4;
	}
#endif

	long long weights[NPREDS]={0}, wsum=0;
	for(int k=0;k<NPREDS;++k)
	{
		weights[k]=
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+1+PAD_SIZE)+k)<<2|kc]+//peNE
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx  +PAD_SIZE)+k)<<2|kc]+//peN
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx-1+PAD_SIZE)+k)<<2|kc];//peNW
		weights[k]=((long long)pr->params[k]<<8)/(weights[k]+1);
		wsum+=weights[k];
	}
	if(wsum)
	{
		pr->pred=0;
		for(int k=0;k<NPREDS;++k)
			pr->pred+=weights[k]*pr->preds[k];
		pr->pred+=(wsum>>1)-1;
		pr->pred/=wsum;
	}
	else
		pr->pred=pr->preds[0];
	
	int
		qx  =quantize_signed(dx  >>(pr->shift[kc]+2)),
		qy  =quantize_signed(dy  >>(pr->shift[kc]+2)),
		q45 =quantize_signed(d45 >>(pr->shift[kc]+2)),
		q135=quantize_signed(d135>>(pr->shift[kc]+2));
	//int
	//	gx=(int)(W-WW+NE-NW+NN-NNE)>>(pr->shift[kc]+2), qx=quantize_signed(gx),//X
	//	gy=(int)(W-NW+N-NN+NE-NNE)>>(pr->shift[kc]+2), qy=quantize_signed(gy),
	//	g45=(int)(((N-W)<<1)+NN-WW)>>(pr->shift[kc]+2), q45=quantize_signed(g45),
	//	g135=(int)(N-NNW+NE-NN)>>(pr->shift[kc]+2), q135=quantize_signed(g135);
	pr->hist_idx=MAXVAR(qx, qy);
	pr->hist_idx=MAXVAR(pr->hist_idx, q45);
	pr->hist_idx=MAXVAR(pr->hist_idx, q135);
	pr->hist_idx=MINVAR(pr->hist_idx, _countof(qlevels));
	
	const int g[]=
	{
		qx,
		qy,
		q45,
		q135,
		quantize_signed((int)eW >>pr->shift[kc]),
		quantize_signed((int)eNW>>pr->shift[kc]),
		quantize_signed((int)eN >>pr->shift[kc]),
		quantize_signed((int)eNE>>pr->shift[kc]),
	};
	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES;++k)
	{
		pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
		int qp=quantize_signed(pr->pred_final>>pr->shift[kc]);
		pr->sse_idx[k]=((_countof(qlevels)+1)<<1)*qp+g[k];
		//pr->sse_idx[k]=CLAMP(0, pr->sse_idx[k], SSE_SIZE-1);
		long long sse_val=pr->sse[SSE_SIZE*k+pr->sse_idx[k]];
		pr->sse_count[k]=(int)(sse_val&0xFFF);
		pr->sse_sum[k]=(int)(sse_val>>12);
		int sse_corr=pr->sse_count[k]?(int)(pr->sse_sum[k]/pr->sse_count[k]):0;
		pr->pred+=sse_corr;
		pr->sse_corr+=sse_corr;
	}
	
	long long pred=pr->pred;
	//if(((eN^eW)|(eN^eNW))<0)//clamp only if signs are different
	//{
		int clamp_lo=(int)N, clamp_hi=(int)N;
		clamp_lo=(int)MINVAR(clamp_lo, W);
		clamp_hi=(int)MAXVAR(clamp_hi, W);
		clamp_lo=(int)MINVAR(clamp_lo, NE);
		clamp_hi=(int)MAXVAR(clamp_hi, NE);
		pred=CLAMP(clamp_lo, pred, clamp_hi);
	//}
	pr->pred=pred;
	pr->pred_final=(int)((pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
	pr->kc=kc;
	pr->kx=kx;
}
static void slic5_update(SLIC5Ctx *pr, int curr, int token)
{
	int kc=pr->kc, kx=pr->kx;

	int error=(curr<<PRED_PREC)-(int)pr->pred;
	pr->esum[kc]+=abs(error);
	LOAD(pr->errors, 0, 0)=error;
	int errors[NPREDS]={0}, kbest=0;
	for(int k=0;k<NPREDS;++k)
	{
		errors[k]=abs((curr<<PRED_PREC)-pr->preds[k]);
		//if(e<0)
		//	LOG_ERROR("");
		pr->pred_errors[(NPREDS*(pr->kym0+pr->kx+PAD_SIZE)+k)<<2|kc]=errors[k];
		if(pr->ky&&pr->kx+1<pr->iw)
			pr->pred_errors[(NPREDS*(pr->kym1+pr->kx+1+PAD_SIZE)+k)<<2|kc]+=errors[k];
		if(errors[kbest]>errors[k])
			kbest=k;
	}
	++pr->params[kbest];

	//if(pr->hist_idx==23&&token==32)//
	//	printf("");

	LOAD(pr->pixels, 0, 0)=curr;
	++pr->hist[pr->cdfsize*pr->hist_idx+token];
	++pr->histsums[pr->hist_idx];
	
	for(int k=0;k<SSE_STAGES;++k)
	{
		++pr->sse_count[k];
		pr->sse_sum[k]+=error;
		if(pr->sse_count[k]>640)
		{
			pr->sse_count[k]>>=1;
			pr->sse_sum[k]>>=1;
		}
		pr->sse[SSE_SIZE*k+pr->sse_idx[k]]=pr->sse_sum[k]<<12|pr->sse_count[k];
	}
}
#undef  LOAD
static void slic5_enc(SLIC5Ctx *pr, int curr, int kc, int kx, int ky)
{
	slic5_predict(pr, kc, kx, ky);

	int error=curr-pr->pred_final;
	int negmask=-(pr->sse_corr<0);
	error^=negmask;
	error-=negmask;
	if(error)//pack sign from CALIC
	{
		int cond=(pr->pred_final<0)==(pr->sse_corr<0);
		int upred=pr->pred_final;
		upred^=-cond;
		upred+=cond;
		upred+=pr->nlevels[pr->kc]>>1;//upred = half + ((pred<0)==(corr<0) ? -pred : pred)
		if(abs(error)<=upred)
			error=(error<<1)^-(error<0);//pack sign
		else
			error=upred+abs(error);//error sign is known
	}
#if 0
	if(error)
	{
		if((ctx[0]<0)!=(ctx[2]<0))
		{
			int upred=(nlevels>>1)+ctx[0];
			if(abs(error)<=upred)
				error=(error<<1)^-(error<0);//pack sign
			else
				error=upred+abs(error);//error sign is known
		}
		else
		{
			int upred=(nlevels>>1)-ctx[0];
			if(abs(error)<=upred)
				error=(error<<1)^-(error<0);//pack sign
			else
				error=upred+abs(error);
		}
	}
#endif

	HybridUint hu;
	hybriduint_encode(error, &hu);
	if(hu.token>=pr->cdfsize)
		LOG_ERROR("Token OOB %d/%d", hu.token, pr->cdfsize);

	ac_enc(pr->ec, hu.token, pr->CDFs+(pr->cdfsize+1)*pr->hist_idx, pr->cdfsize, 1);
	if(hu.nbits)
	{
		int bypass=hu.bypass, nbits=hu.nbits;
		while(nbits>8)
		{
			ac_enc(pr->ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 0x10000>>8);
			nbits-=8;
		}
		ac_enc(pr->ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 0x10000>>nbits);
	}

	slic5_update(pr, curr, hu.token);
}
static int slic5_dec(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	slic5_predict(pr, kc, kx, ky);

	int token=ac_dec(pr->ec, pr->CDFs+(pr->cdfsize+1)*pr->hist_idx, pr->cdfsize, 1);
	int error=token;
	if(error>=(1<<SLIC5_CONFIG_EXP))
	{
		error-=1<<SLIC5_CONFIG_EXP;
		int lsb=error&((1<<SLIC5_CONFIG_LSB)-1);
		error>>=SLIC5_CONFIG_LSB;
		int msb=error&((1<<SLIC5_CONFIG_MSB)-1);
		error>>=SLIC5_CONFIG_MSB;
		int nbits=error+SLIC5_CONFIG_EXP-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB), n=nbits;
		int bypass=0;
		while(n>8)
		{
			n-=8;
			bypass|=ac_dec(pr->ec, 0, 1<<8, 0x10000>>8)<<n;
		}
		bypass|=ac_dec(pr->ec, 0, 1<<n, 0x10000>>n);
		//int bypass=nbits&7?ac_dec(pr->ec, 0, 1<<(nbits&7), 0x10000>>nbits):0;
		//for(int k=0, nb=nbits>>3;k<nb;++k)
		//	bypass=bypass<<8|ac_dec(pr->ec, 0, 256, 0x10000>>8);
		error=1;
		error<<=SLIC5_CONFIG_MSB;
		error|=msb;
		error<<=nbits;
		error|=bypass;
		error<<=SLIC5_CONFIG_LSB;
		error|=lsb;
	}
	int cond=(pr->pred_final<0)==(pr->sse_corr<0);
	int upred=pr->pred_final;
	upred^=-cond;
	upred+=cond;
	upred+=pr->nlevels[kc]>>1;
	if(error<=(upred<<1))//sign pack from CALIC
		error=error>>1^-(error&1);//unpack sign
	else
	{
		error=error-upred;//error sign is known
		error^=-cond;
		error+=cond;
	}
#if 0
	if((ctx[0]<0)!=(ctx[2]<0))
	{
		int upred=(nlevels>>1)+ctx[0];
		if(sym<=(upred<<1))
			error=sym>>1^-(sym&1);
		else
			error=sym-upred;
	}
	else
	{
		int upred=(nlevels>>1)-ctx[0];
		if(sym<=(upred<<1))
			error=sym>>1^-(sym&1);
		else
			error=upred-sym;
	}
#endif
	int negmask=-(pr->sse_corr<0);
	error^=negmask;
	error-=negmask;
	int curr=error+pr->pred_final;
	slic5_update(pr, curr, token);
	return curr;
}
static const Image *guide=0;
int t47_encode(Image const *src, ArrayHandle *data, int loud)
{
	guide=src;
	double t_start=time_sec();
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T47 SLIC5-AC  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, src->nch, src->iw, src->ih, maxdepth);
	}
	DList list;
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	int success=slic5_init(&pr, src->iw, src->ih, src->nch, src->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	
#ifndef SLIC5_DISABLE_PALETTE
	unsigned short *palettes[4]={0};
	unsigned short pal_sizes[4]={0};
	//memset(header.palettesizes, 0, sizeof(header.palettesizes));
	int *pal_hist=(int*)malloc(maxlevels*sizeof(int));
	if(!pal_hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int kc=0;kc<nch;++kc)
	{
		int depth=im2->depth[kc], nlevels=1<<depth;
		memset(hist, 0, nlevels*sizeof(int));
		for(ptrdiff_t k=0;k<res;++k)
		{
			int val=im2->data[k<<2|kc]+(nlevels>>1);
			val=CLAMP(0, val, nlevels-1);
			++hist[val];
		}
		int palettesize=0;
		for(int k=0;k<nlevels;++k)
			palettesize+=hist[k]!=0;
		if(palettesize>((nlevels+15)>>4))
			pal_sizes[kc]=0;
		else
		{
			pal_sizes[kc]=palettesize;
			im2->depth[kc]=floor_log2_32(palettesize-1)+1;
			palettes[kc]=(unsigned short*)malloc(palettesize*sizeof(short));
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			for(int k=0, idx=0;k<nlevels;++k)
			{
				if(hist[k])
				{
					palettes[kc][idx]=k;
					++idx;
				}
			}
			int half=(palettesize+1)>>1;//ceil_half
			for(int k=0;k<res;++k)
			{
				unsigned short val=im2->data[k<<2|kc];
				int idx=0;
				int L=0, R=palettesize-1;
				while(L<=R)
				{
					idx=(L+R)>>1;
					if(palettes[kc][idx]<val)
						L=idx+1;
					else if(palettes[kc][idx]>val)
						R=idx-1;
					else
						break;
				}
				im2->data[k<<2|kc]=idx-half;
			}
		}
	}
	free(pal_hist);
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			//if(kx==7&&ky==3)//
			//if(kx==674&&ky==2)//
			//	printf("");

			int kc=0;
			if(src->nch>=3)
			{
				int r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];

				r-=g;
				b-=g;
				g+=(r+b)>>2;

				slic5_enc(&pr, g, 0, kx, ky);//Y
				slic5_enc(&pr, b, 1, kx, ky);//Cb
				slic5_enc(&pr, r, 2, kx, ky);//Cr
				kc+=3;
			}
			while(kc<src->nch)//gray/alpha
			{
				int pixel=src->data[idx|kc];
				slic5_enc(&pr, pixel, kc, kx, ky);
				++kc;
			}
		}
	}
	ac_enc_flush(&ec);
	slic5_free(&pr);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);
	}
	dlist_clear(&list);
	return 1;
}
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
#ifndef SLIC5_DISABLE_PALETTE
	unsigned short pal_sizes[4];
	//size_t emithistsize=cdfsize*sizeof(short[_countof(qlevels)+1]);
	if(srclen<sizeof(pal_sizes))
	{
		LOG_ERROR("Corrupt file\n");
		//printf("File smaller than header size\n");
		return 0;
	}
	memcpy(pal_sizes, data, sizeof(pal_sizes));
	data+=sizeof(pal_sizes);
	srclen-=sizeof(pal_sizes);
	unsigned short *palettes[4]={0};
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(pal_sizes[kc]>0)
		{
			int bytesize=pal_sizes[kc]*sizeof(short);
			palettes[kc]=(unsigned short*)malloc(bytesize);
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			if(bytesize>srclen)
			{
				LOG_ERROR("Invalid file\n");
				return 0;
			}
			memcpy(palettes[kc], data, bytesize);
			data+=bytesize;
			srclen-=bytesize;
		}
	}
#endif
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	int success=slic5_init(&pr, dst->iw, dst->ih, dst->nch, dst->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=4)
		{
			//if(kx==7&&ky==3)//
			//if(kx==674&&ky==2)//
			//if(kx==1173&&ky==11)//
			//if(kx==0&&ky==0)//
			//	printf("");

			int kc=0;
			if(dst->nch>=3)
			{
				int g=slic5_dec(&pr, 0, kx, ky);//Y
				int b=slic5_dec(&pr, 1, kx, ky);//Cb
				int r=slic5_dec(&pr, 2, kx, ky);//Cr
					
				g-=(r+b)>>2;
				b+=g;
				r+=g;
				//if(guide&&(r!=guide->data[(guide->iw*ky+kx)<<2|0]||g!=guide->data[(guide->iw*ky+kx)<<2|1]||b!=guide->data[(guide->iw*ky+kx)<<2|2]))
				//	LOG_ERROR("Guide error XY %d %d", kx, ky);

				dst->data[idx|0]=r;
				dst->data[idx|1]=g;
				dst->data[idx|2]=b;

				kc+=3;
			}
			while(kc<dst->nch)//gray/alpha
			{
				int pixel=slic5_dec(&pr, kc, kx, ky);
				dst->data[idx|kc]=pixel;
				++kc;
			}
		}
	}
#ifndef SLIC5_DISABLE_PALETTE
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(palettes[kc])
		{
			for(int k=0;k<res;++k)
			{
				int val=dst->data[k<<2|kc];
				val+=(pal_sizes[kc]+1)>>1;
				if((unsigned)val>=pal_sizes[kc])
				{
					LOG_ERROR("Palette error");
					return 0;
				}
				dst->data[k<<2|kc]=palettes[kc][val];
			}
			free(palettes[kc]);
		}
	}
#endif
	if(loud)
	{
		printf("\n");
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}