#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

	#define SLIC4_DISABLE_PALETTE
//	#define SLIC4_USE_SIMPLE_HYBRID//4 1 0

#define SLIC4_CONFIG_EXP 5
#define SLIC4_CONFIG_MSB 2
#define SLIC4_CONFIG_LSB 0

static const int qlevels[]=
{
	1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47, 55, 63, 79, 95, 111, 127
	//1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392, 500
};

typedef struct TempHybridStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} TempHybrid;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
static void hybriduint_encode(unsigned val, TempHybrid *hu)
{
	int token, bypass, nbits;
#ifdef SLIC4_USE_SIMPLE_HYBRID
	if(val<16)
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		nbits=floor_log2_32(val>>2)+1;
		token=16+((nbits-3)<<1)+(val>>nbits)-2;
		bypass=val&((1<<nbits)-1);
	}
#else
	if(val<(1<<SLIC4_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC4_CONFIG_EXP) + (
				(lgv-SLIC4_CONFIG_EXP)<<(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC4_CONFIG_MSB))<<SLIC4_CONFIG_LSB|
				(mantissa&((1<<SLIC4_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
		bypass=val>>SLIC4_CONFIG_LSB&((1LL<<nbits)-1);
	}
#endif
	hu->token=(unsigned short)token;
	hu->nbits=(unsigned short)nbits;
	hu->bypass=bypass;
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
#define HASH_FUNC(A, B) ((A+B+1)*(A+B)/2+B)
#define NPREDS 6
#define PRED_PREC 8
#define PARAM_PREC 8
#define SSE_STAGES 6
//#define SSE_SIZE 0x100000
//#define SSE_SIZE ((_countof(qlevels)+1)*(_countof(qlevels)+1)*(_countof(qlevels)+1)*(_countof(qlevels)+1)<<1)
//#define SSE_SIZE 1024
#define SSE_SIZE ((_countof(qlevels)+1)*(_countof(qlevels)+1)<<2)
//#define SSE_SIZE (_countof(qlevels)+1)

#define TAG_HIST_START 7
#define TAG_HIST_LEN 4
typedef struct PredictorCtxStruct
{
	int iw;
	int *pred_errors, *errors;
	long long esum;//sum of abs errors so far in current row
	int sigma;//average of abs errors from prev row
	int shift;//depth compensation
	long long *sse, sse_idx[SSE_STAGES], sse_sum[SSE_STAGES];
	int sse_count[SSE_STAGES];
	int preds[NPREDS];
	long long pred;
	int pred_final;
	//long long r_weights[NPREDS];
	int hist_idx;
	int sse_corr;
	int depth, nlevels, min_allowed, max_allowed;
	int kx, ky, ky1;
} PredictorCtx;
static int slic4_pred_init(PredictorCtx *pr, int iw)
{
	pr->iw=iw;
	pr->pred_errors=(int*)malloc(iw*sizeof(int[NPREDS*2]));
	pr->errors=(int*)malloc(iw*sizeof(int[2]));
	pr->sse=(long long*)malloc(sizeof(long long[SSE_STAGES*SSE_SIZE]));
	if(!pr->pred_errors||!pr->errors||!pr->sse)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	return 1;
}
static void slic4_pred_free(PredictorCtx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->sse);
}
static void slic4_pred_nextchannel(PredictorCtx *pr, int depth)
{
	memset(pr->pred_errors, 0, pr->iw*sizeof(int[NPREDS*2]));
	memset(pr->errors, 0, pr->iw*sizeof(int[2]));
	memset(pr->sse, 0, sizeof(long long[SSE_STAGES*SSE_SIZE]));
	pr->sigma=0;
	pr->esum=0;
	pr->depth=depth;
	pr->nlevels=1<<depth;
	pr->min_allowed=-(pr->nlevels>>1);
	pr->max_allowed=(pr->nlevels>>1)-1;
}
static void slic4_pred_nextrow(PredictorCtx *pr, int ky)
{
	pr->sigma=(int)(pr->esum/pr->iw);
	pr->esum=0;
	//pr->sse_idx=0;
	if(ky)
	{
		pr->shift=FLOOR_LOG2_P1(pr->sigma)-PRED_PREC;
		pr->shift=MAXVAR(8, pr->shift)-8+PRED_PREC;
	}
	else
		pr->shift=MAXVAR(8, pr->depth)-8+PRED_PREC;
	pr->ky=ky;
	pr->ky1=ky&1;
}
static const short wp_params[]=
{
	// 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x0004,
	// 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,

	 0x0DB8,  0x0E22,  0x181F,  0x0BF3,//wp weights
	//0x181F,//clamp grad
	 0x181F,//journal GAP
	 0x181F,//CALIC GAP

	-0x005C,
	-0x005B,
	 0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,

	//0xd, 0xc, 0xc, 0xb,
	//8,
	//8,
	//4, 0, 3, 23, 2,
};
#if 0
void pred_jxl(Image *srcim, int fwd, int enable_ma)
{
	const short *params=wp_params;
	ptrdiff_t res=(ptrdiff_t)srcim->iw*srcim->ih;
	int *temp_w10=(int*)malloc((size_t)srcim->iw*10*sizeof(int));
	int *dst=(int*)malloc(res*sizeof(int[4]));
	if(!temp_w10||!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int errorbuflen=srcim->iw<<1, rowlen=srcim->iw<<2;
	int *error=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int iw=srcim->iw;
	int *src=srcim->data;
	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<srcim->depth[kc];
		int idx=kc;
		const int *pixels=fwd?src:dst, *errors=fwd?dst:src;
		for(int ky=0;ky<srcim->ih;++ky)
		{
			int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
			for(int kx=0;kx<iw;++kx, idx+=4)
			{
				int pred, curr;
			
				int
					ctt      =         ky-2>=0?pixels[idx-rowlen*2]:0,
					ctopleft =kx-1>=0&&ky-1>=0?pixels[idx-rowlen-4]:0,
					ctop     =kx  <iw&&ky-1>=0?pixels[idx-rowlen  ]:0,
					ctopright=kx+1<iw&&ky-1>=0?pixels[idx-rowlen+4]:0,
					cleft    =kx-1>=0         ?pixels[idx       -4]:0;

				//if(kx==(iw>>1)&&ky==(ih>>1))
				//	kx=iw>>1;

				//w0   w1   w2   w3
				//p1C  p2c
				//p3Ca p3Cb p3Cc p3Cd p3Ce
			
				int weights[4];//fixed 23.8 bit
				for(int k=0;k<4;++k)
				{
					int w=
						(ky-1>=0         ?pred_errors[k][prevrow+kx  ]:0)+//eN
						(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+//eNE
						(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);//eNW
					weights[k]=(params[k]<<8)/(w+1);
				}

				int
					etop=ky-1>=0?error[prevrow+kx]:0,
					eleft=kx-1>=0?error[currrow+kx-1]:0,
					etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
					etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
					etopplusleft=etop+eleft;
				long long predictions[]=//fixed 23.8 bit
				{
					(cleft+ctopright-ctop)<<8,
					(ctop<<8)-((etopplusleft+etopright)*params[4]>>8),
					(cleft<<8)-((etopplusleft+etopleft)*params[5]>>8),
					(ctop<<8)-((etopleft*params[6]+etop*params[7]+etopright*params[8]+((ctt-ctop)<<8)*params[9]+((ctopleft-cleft)<<8)*params[10])>>8),
					//(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
				};

				int sum=weights[0]+weights[1]+weights[2]+weights[3];
				if(sum)
					pred=(int)((predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum);
				else
					pred=(int)predictions[0];

				int vmin=cleft, vmax=cleft;
				if(vmin>ctopright)vmin=ctopright;
				if(vmax<ctopright)vmax=ctopright;

				if(vmin>ctop)vmin=ctop;
				if(vmax<ctop)vmax=ctop;

				vmin<<=8;
				vmax<<=8;

				//if(kc==0&&kx==0&&ky==1)//
				//if(kc==0&&kx==256&&ky==256)//
				//if(kc==0&&kx==1&&ky==0)//
				//if(kc==0&&kx==4&&ky==2)//
				//	printf("");

				pred=CLAMP(vmin, pred, vmax);

				int pred_final=(pred+127)>>8;
				pred_final^=-fwd;
				pred_final+=fwd;
				pred_final+=src[idx];
				if(enable_ma)
				{
					pred_final+=nlevels>>1;
					pred_final&=nlevels-1;
					pred_final-=nlevels>>1;
				}
				dst[idx]=pred_final;
				curr=pixels[idx]<<8;

				error[currrow+kx]=curr-pred;
				for(int k=0;k<4;++k)
				{
					int e=abs(curr-(int)predictions[k]);
					pred_errors[k][currrow+kx]=e;
					if(ky&&kx+1<iw)
						pred_errors[k][prevrow+kx+1]+=e;
				}
			}
		}
	}
	memcpy(src, dst, res*sizeof(int[4]));
	free(temp_w10);
	free(dst);
}
#endif
static void slic4_predict(Image const *im, int kc, int kx, int ky, PredictorCtx *pr)
{
	const int *pixels=im->data;
	pr->kx=kx;
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)im->iw&&(unsigned)(ky+(Y))<(unsigned)im->ih?BUF[(im->iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0
	int
	//	NNNNNN  =LOAD(pixels,  0, -6),
	//	NNNNWWWW=LOAD(pixels, -4, -4),
	//	NNNN    =LOAD(pixels,  0, -4),
	//	NNNNEEEE=LOAD(pixels,  4, -4),
	//	NNNWWW  =LOAD(pixels, -3, -3),
	//	NNN     =LOAD(pixels,  0, -3),
	//	NNNEEE  =LOAD(pixels,  3, -3),
		NNWW    =LOAD(pixels, -2, -2),
		NNW     =LOAD(pixels, -1, -2),
		NN      =LOAD(pixels,  0, -2),
		NNE     =LOAD(pixels,  1, -2),
		NNEE    =LOAD(pixels,  2, -2),
		NWW     =LOAD(pixels, -2, -1),
		NW      =LOAD(pixels, -1, -1),
		N       =LOAD(pixels,  0, -1),
		NE      =LOAD(pixels,  1, -1),
	//	NEE     =LOAD(pixels,  2, -1),
	//	NEEEE   =LOAD(pixels,  4, -1),
	//	NEEEEE  =LOAD(pixels,  5, -1),
	//	NEEEEEE =LOAD(pixels,  6, -1),
	//	NEEEEEEE=LOAD(pixels,  7, -1),
	//	WWWWWW  =LOAD(pixels, -6,  0),
	//	WWWW    =LOAD(pixels, -4,  0),
	//	WWW     =LOAD(pixels, -3,  0),
		WW      =LOAD(pixels, -2,  0),
		W       =LOAD(pixels, -1,  0);
	//long long
	//	idx=im->iw*ky+kx,
	//	NN=ky>=2 ?im->data[(idx-im->iw*2)<<2|kc]<<PRED_PREC:0,
	//	WW=kx>=2 ?im->data[(idx       -2)<<2|kc]<<PRED_PREC:0,
	//	N =ky    ?im->data[(idx-im->iw  )<<2|kc]<<PRED_PREC:0,
	//	W =kx    ?im->data[(idx       -1)<<2|kc]<<PRED_PREC:0,
	//	NW=kx&&ky?im->data[(idx-im->iw-1)<<2|kc]<<PRED_PREC:0,
	//	NE=kx+1<im->iw  &&ky   ?im->data[(idx-im->iw  +1)<<2|kc]<<PRED_PREC:0,
	//	NWW =kx>=2      &&ky>=1?im->data[(idx-im->iw  -2)<<2|kc]<<PRED_PREC:0,
	//	NEE =kx+2<im->iw&&ky>=1?im->data[(idx-im->iw  +2)<<2|kc]<<PRED_PREC:0,
	//	NNW =kx>=1      &&ky>=2?im->data[(idx-im->iw*2-1)<<2|kc]<<PRED_PREC:0,
	//	NNE =kx+1<im->iw&&ky>=2?im->data[(idx-im->iw*2+1)<<2|kc]<<PRED_PREC:0,
	//	NNWW=kx>=2      &&ky>=2?im->data[(idx-im->iw*2-2)<<2|kc]<<PRED_PREC:0,
	//	NNEE=kx+2<im->iw&&ky>=2?im->data[(idx-im->iw*2+2)<<2|kc]<<PRED_PREC:0;
	
	long long
		eN =ky             ?pr->errors[pr->iw*!pr->ky1+kx  ]:0,//error = curr - pred
		eW =kx             ?pr->errors[pr->iw* pr->ky1+kx-1]:0,
		eNW=kx&&ky         ?pr->errors[pr->iw*!pr->ky1+kx-1]:0,
		eNE=kx+1<pr->iw&&ky?pr->errors[pr->iw*!pr->ky1+kx+1]:0;

	int j=-1;
	pr->preds[++j]=(int)(W+NE-N);
	pr->preds[++j]=(int)(N-((eN+eW+eNE)*wp_params[NPREDS+0]>>PARAM_PREC));
	pr->preds[++j]=(int)(W-((eN+eW+eNW)*wp_params[NPREDS+1]>>PARAM_PREC));
	pr->preds[++j]=(int)(N-((eNW*wp_params[NPREDS+2]+eN*wp_params[NPREDS+3]+eNE*wp_params[NPREDS+4]+(NN-N)*wp_params[NPREDS+5]+(NW-W)*wp_params[NPREDS+6])>>PARAM_PREC));
	
#if 0
	//pr->preds[++j]=(((N+W)<<1)+NW+NE+NN+WW+4)>>3;//X

	pr->preds[++j]=(int)(W -(eW *wp_params[NPREDS+ 7]>>8));
	pr->preds[++j]=(int)(N -(eN *wp_params[NPREDS+ 8]>>8));
	pr->preds[++j]=(int)(NW-(eNW*wp_params[NPREDS+ 9]>>8));
	pr->preds[++j]=(int)(NE-(eNE*wp_params[NPREDS+10]>>8));

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
	//++j, pr->preds[j]=(int)(N+W-NW), pr->preds[j]=(int)MEDIAN3(N, W, pr->preds[j]);
	
	++j;
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int sum2=(dy+dx)>>pr->shift, diff=(dy-dx)>>pr->shift, diff2=NE-NW;
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
	int disc=(dy-dx)>>pr->shift;
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
#if 0
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int sum2=dy+dx, diff=dy-dx, diff2=NE-NW;
	pr->preds[0]=0;
	if(sum2>(32<<(PRED_PREC+1)))
		pr->preds[0]=(dx*N+dy*W)/sum2+diff2/8;
	else if(diff>(12<<(PRED_PREC+1)))
		pr->preds[0]=(N+2*W)/3+diff2/8;
	else if(diff<-(12<<(PRED_PREC+1)))
		pr->preds[0]=(2*N+W)/3+diff2/8;
	else
		pr->preds[0]=(N+W)/2+diff2/8;
	diff=d45-d135;
	if(diff>(32<<(PRED_PREC+1)))
		pr->preds[0]+=diff2/8;
	else if(diff>(16<<(PRED_PREC+1)))
		pr->preds[0]+=diff2/16;
	else if(diff<-(32<<(PRED_PREC+1)))
		pr->preds[0]-=diff2/8;
	else if(diff<-(16<<(PRED_PREC+1)))
		pr->preds[0]-=diff2/16;
	pr->preds[0]+=(eN+eW)>>4;

	pr->preds[1]=N+W-NW;
	pr->preds[1]=MEDIAN3(N, W, pr->preds[1]);

	pr->preds[2]=W+NE-N+((eN+eW)>>4);

	pr->preds[3]=N+(eN>>2);

	pr->preds[4]=W+(eW>>2);
#endif

	long long weights[NPREDS], wsum=0;
	for(int k=0;k<NPREDS;++k)
	{
		weights[k]=
			(ky             ?pr->pred_errors[NPREDS*(pr->iw*!pr->ky1+kx  )+k]:0)+//peN
			(ky&&kx+1<pr->iw?pr->pred_errors[NPREDS*(pr->iw*!pr->ky1+kx+1)+k]:0)+//peNE
			(ky&&kx>=1      ?pr->pred_errors[NPREDS*(pr->iw*!pr->ky1+kx-1)+k]:0);//peNW
		weights[k]=((long long)wp_params[k]<<8)/(weights[k]+1);
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
	
	//int p=abs((int)eN);
	//if(p>abs((int)eW))p=abs((int)eW);
	//if(p>abs((int)eNW))p=abs((int)eNW);
	//if(p>abs((int)eNE))p=abs((int)eNE);
	//pr->sse_idx=pr->hist_idx=quantize(p);
	int
		qx  =quantize_signed(dx  >>(pr->shift+2)),
		qy  =quantize_signed(dy  >>(pr->shift+2)),
		q45 =quantize_signed(d45 >>(pr->shift+2)),
		q135=quantize_signed(d135>>(pr->shift+2));
	//int
	//	gx=(int)(W-WW+NE-NW+NN-NNE)>>(pr->shift+2), qx=quantize_signed(gx),
	//	gy=(int)(W-NW+N-NN+NE-NNE)>>(pr->shift+2), qy=quantize_signed(gy),
	//	g45=(int)(((N-W)<<1)+NN-WW)>>(pr->shift+2), q45=quantize_signed(g45),
	//	g135=(int)(N-NNW+NE-NN)>>(pr->shift+2), q135=quantize_signed(g135);
	pr->hist_idx=MAXVAR(qx, qy);
	pr->hist_idx=MAXVAR(pr->hist_idx, q45);
	pr->hist_idx=MAXVAR(pr->hist_idx, q135);
	pr->hist_idx=MINVAR(pr->hist_idx, (int)_countof(qlevels));
	
	qx=(int)(W-WW+NE-NW+NN-NNE)>>(pr->shift+2), qx=quantize_signed(qx);
	qy=(int)(W-NW+N -NN+NE-NNE)>>(pr->shift+2), qy=quantize_signed(qy);
	const int g[]=
	{
		qx,
	//	q45,
		qy,
	//	q135,
		quantize_signed((int)eW>>pr->shift),
		quantize_signed((int)eNW>>pr->shift),
		quantize_signed((int)eN>>pr->shift),
		quantize_signed((int)eNE>>pr->shift),
	};
	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES;++k)
	{
		pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
		int qp=quantize_signed(pr->pred_final>>pr->shift);
		pr->sse_idx[k]=((_countof(qlevels)+1)<<1)*qp+g[k];
		//pr->sse_idx[k]=CLAMP(0, pr->sse_idx[k], SSE_SIZE-1);
		long long sse_val=pr->sse[SSE_SIZE*k+pr->sse_idx[k]];
		pr->sse_count[k]=(int)(sse_val&0xFFF);
		pr->sse_sum[k]=(int)(sse_val>>12);
		int sse_corr=pr->sse_count[k]?(int)(pr->sse_sum[k]/pr->sse_count[k]):0;
		pr->pred+=sse_corr;
		pr->sse_corr+=sse_corr;
	}

	//pr->sse_idx=(_countof(qlevels)+1)*(_countof(qlevels)+1)*qy*qx+q45*q135;
	//pr->sse_idx%=SSE_SIZE;

	//pr->hist_idx=pr->sse_idx=0;

	//GAP
#if 0
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	//int sum2=dy+dx, diff=dy-dx, diff2=NE-NW;
	//int pred=0;
	//if(sum2>32)
	//	pred=(dx*N+dy*W)/sum2+diff2/8;
	//else if(diff>12)
	//	pred=(N+2*W)/3+diff2/8;
	//else if(diff<-12)
	//	pred=(2*N+W)/3+diff2/8;
	//else
	//	pred=(N+W)/2+diff2/8;
	//diff=d45-d135;
	//if(diff>32)
	//	pred+=diff2/8;
	//else if(diff>16)
	//	pred+=diff2/16;
	//else if(diff<-32)
	//	pred-=diff2/8;
	//else if(diff<-16)
	//	pred-=diff2/16;

	//int disc=(dy-dx)>>shift;
	//int pred;
	//if(disc>80)
	//	pred=W;
	//else if(disc<-80)
	//	pred=N;
	//else
	//{
	//	pred=(N+W)/2+(NE-NW)/4;
	//	if(disc>32)
	//		pred=(pred+W)/2;
	//	else if(disc>8)
	//		pred=(3*pred+W)/4;
	//	else if(disc<-32)
	//		pred=(pred+N)/2;
	//	else if(disc<-8)
	//		pred=(3*pred+N)/4;
	//}

	//int pred=N+W-NW;
	//pred=MEDIAN3(N, W, pred);
	//int pred=dy<dx?N:W;
	int sum2=dx+dy+d45+d135, pred;
	if(sum2)
		pred=(dx*N+dy*W+d45*NE+d135*NW)/sum2;
	else
	{
		pred=N+W-NW;
		pred=MEDIAN3(N, W, pred);
	}

	int ctx[]=
	{
		quantize_signed(prev_error>>shift),
		quantize_signed(dx>>(shift+2)),
		quantize_signed(dy>>(shift+2)),
		quantize_signed((((N-W)<<1)+NN-WW)>>(shift+2)),
		quantize_signed((N-NNW+NE-NN)>>(shift+2)),
	};
	//int hash=HASH_FUNC(ctx[0], ctx[1]);
	//hash=HASH_FUNC(hash, ctx[2]);
	//hash=HASH_FUNC(hash, ctx[3]);
	//hash=HASH_FUNC(hash, ctx[4]);
	//hash&=SSE_SIZE-1;
	//int hash=(ctx[4]&15)<<16|(ctx[3]&15)<<12|(ctx[2]&15)<<8|(ctx[1]&15)<<4|(ctx[0]&15);
	//int hash=ctx[4]+ctx[3]+ctx[2]+ctx[1]+ctx[0];
	int hist_ctx=ctx[0]>>1;
	int sse_ctx=0;

	//int qx=quantize_signed((W-WW+NE-NW+NN-NNE)>>shift), qy=quantize_signed((W-NW+N-NN+NE-NNE)>>shift),
	//	q45=quantize_signed((((N-W)<<1)+NN-WW)>>shift), q135=quantize_signed((N-NNW+NE-NN)>>shift);
	//int hist_ctx=MAXVAR(qx, qy);
	//hist_ctx=MAXVAR(hist_ctx, q45);
	//hist_ctx=MAXVAR(hist_ctx, q135);
	//hist_ctx=MINVAR(hist_ctx, _countof(qlevels));
	//int sse_ctx=(_countof(qlevels)+1)*(_countof(qlevels)+1)*qy*qx+q45*q135;
	//sse_ctx%=SSE_SIZE;
#endif
#if 0
	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);
	int qx=quantize_signed((W-WW+NE-NW+NN-NNE)>>shift), qy=quantize_signed((W-NW+N-NN+NE-NNE)>>shift),
		q45=quantize_signed((((N-W)<<1)+NN-WW)>>shift), q135=quantize_signed((N-NNW+NE-NN)>>shift);
	//int qx=quantize_signed((W-WW+NE-NW)>>shift), qy=quantize_signed((((N-NN)<<1)+W-NW)>>shift), q45=quantize_signed((((N-W)<<1)+NN-WW)>>shift);
	//int qx=quantize_signed((W-WW+NE-NW)>>shift), qy=quantize_signed((N-NN+W-NW)>>shift), q45=quantize_signed((N-W+NN-WW)>>shift);
	//int qx=quantize_signed(W-WW+N-((NE+NW)>>1)), qy=quantize_signed(((N-NN)<<1)+W-NW), q45=quantize_signed(((N-W)<<1)+NN-WW);
	//int gradx=abs(W-WW+NE-NW), grady=abs(((N-NN)<<1)+W-NW), grad45=abs(((N-W)<<1)+NN-WW);
	//int gradx=abs(W-WW)+abs(N-NW)+abs(NE-N), grady=(abs(N-NN)<<1)+abs(W-NW), grad45=(abs(N-W)<<1)+abs(NN-WW);
	//int qx=quantize(gradx>>(shift+4)), qy=quantize(grady>>(shift+4)), q45=quantize(grad45>>(shift+4));
	int hist_ctx=MAXVAR(qx, qy);
	//int hist_ctx=qx+qy;
	//int hist_ctx=(qx+qy)>>1;
	hist_ctx=MAXVAR(hist_ctx, q45);
	hist_ctx=MAXVAR(hist_ctx, q135);
	//hist_ctx>>=1;
	hist_ctx=MINVAR(hist_ctx, _countof(qlevels));
	int sse_ctx=(_countof(qlevels)+1)*(_countof(qlevels)+1)*qy*qx+q45*q135;
	//int sse_ctx=qy*qx+q45*q135;
	//int sse_ctx=(_countof(qlevels)+1)*((_countof(qlevels)+1)*qy+qx)+q45;
	//int sse_ctx=(_countof(qlevels)+1)*((_countof(qlevels)+1)*qy+qx)+((q45+q135)>>1);
	//int sse_ctx=((_countof(qlevels)+1)*qy*qx+q45);
	sse_ctx%=SSE_SIZE;
#endif
#if 0
	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);
	int
		qx=quantize_signed((W-WW)>>shift)*quantize_signed((NE-NW)>>shift),
		qy=quantize_signed((N-NN)>>shift)*quantize_signed((W-NW)>>shift),
		q45=quantize_signed((N-W)>>shift)*quantize_signed((NN-WW)>>shift);
	if(qx>_countof(qlevels))qx=_countof(qlevels);
	if(qy>_countof(qlevels))qy=_countof(qlevels);
	if(q45>_countof(qlevels))q45=_countof(qlevels);
	int hist_ctx=MAXVAR(qx, qy);
	hist_ctx=MAXVAR(hist_ctx, q45);
	int sse_ctx=((_countof(qlevels)+1)*((_countof(qlevels)+1)*qy+qx)+q45);
#endif
#if 0
	int vmax=MAXVAR(N, W), vmin=MINVAR(N, W), energy=vmax-vmin;
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	//int pred=(N+W)/2;
	int
		hist_ctx=quantize((energy+abs(W-WW)+abs(N-NN)+abs(prev_error)*2)>>(shift+2)),
		//hist_ctx=quantize((abs(W-WW)+abs(N-NN)+abs(prev_error))>>2),
		//hist_ctx=quantize(abs(prev_error)),
		//hist_ctx=quantize(energy>>shift),
		sse_ctx=(_countof(qlevels)+1)*(_countof(qlevels)+1)*quantize((kc?im->data[idx<<2]:abs(N-W))>>(shift+4))+quantize((abs(W-WW)+abs(N-NN))>>(shift+2))*quantize(abs(prev_error));//added info must be new
		//sse_ctx=(_countof(qlevels)+1)*(_countof(qlevels)+1)*quantize((kc?im->data[idx<<2]:abs(N-W))>>(shift+4))+quantize((abs(W-WW)+abs(N-NN)+abs(NE-NW))>>(shift+2))*quantize(abs(prev_error));
		//sse_ctx=(_countof(qlevels)+1)*quantize((abs(W-WW)+abs(N-NN))>>(shift+2))+quantize(abs(prev_error));
		//sse_ctx=hist_ctx;
		//sse_ctx=(_countof(qlevels)+1)*quantize(abs(W-WW)>>(shift+1))+quantize(abs(N-NN)>>(shift+1));
		//sse_ctx=quantize(abs(W-WW)>>(shift+1))*quantize(abs(N-NN)>>(shift+1));
		//sse_ctx=0;//X
	//sse_ctx+=(_countof(qlevels)+1)*(_countof(qlevels)+1)*quantize((kc?im->data[idx<<2]:(N+W)>>1)>>(shift+4));//chroma channel: luma is always available, otherwise use neighbor luma
	//if(kx&&kc)//luma affects context
	//	sse_ctx+=(_countof(qlevels)+1)*(_countof(qlevels)+1)*quantize(im->data[idx<<2]>>(shift+2));
		//sse_ctx*=quantize((im->data[(idx-1)<<2|kc]+im->data[(idx-2)<<2|kc]+im->data[(idx-3)<<2|kc]+W)>>2);
#endif

	//long long sse_cell=pr->sse[pr->sse_idx];
	//pr->sse_count=(int)(sse_cell&0xFFF);
	//pr->sse_sum=sse_cell>>12;
	//pr->sse_corr=pr->sse_count?(int)(pr->sse_sum/pr->sse_count):0;
	//pr->pred+=pr->sse_corr;

	//if(kc==0&&kx==0&&ky==0)//
	//	printf("");
	//if(kc==0&&kx==1&&ky==1)//
	//if(kc==0&&kx==4&&ky==2)//
	//if(kc==0&&kx==256&&ky==256)//
	//	printf("");
	
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

	//if(guide)
	//{
	//	int nlevels=1<<guide->depth[kc];
	//	int error=pixels->data[idx<<2|kc]-pr->pred_final;
	//	error+=nlevels>>1;
	//	error&=nlevels-1;
	//	error-=nlevels>>1;
	//	if(error!=guide->data[idx<<2|kc])
	//		LOG_ERROR("Pred error");
	//}
	//pred+=corr;
	//ret[0]=pred;
	//ret[1]=hist_ctx;
	//ret[2]=corr;
}
static void slic4_pred_update(PredictorCtx *pr, int curr)
{
	int error=(curr<<PRED_PREC)-(int)pr->pred;
	pr->esum+=abs(error);
	pr->errors[pr->iw*pr->ky1+pr->kx]=error;
	for(int k=0;k<NPREDS;++k)
	{
		int e=abs((curr<<PRED_PREC)-pr->preds[k]);
		pr->pred_errors[NPREDS*(pr->iw*pr->ky1+pr->kx)+k]=e;
		if(pr->ky&&pr->kx+1<pr->iw)
			pr->pred_errors[NPREDS*(pr->iw*!pr->ky1+pr->kx+1)+k]+=e;
	}
	
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
	//++pr->sse_count;
	//pr->sse_sum+=error;
	//if(pr->sse_count>640)
	//{
	//	pr->sse_count>>=1;
	//	pr->sse_sum>>=1;
	//}
	//pr->sse[pr->sse_idx]=pr->sse_sum<<12|pr->sse_count;
}
static const Image *guide=0;
int t46_encode(Image const *src, ArrayHandle *data, int loud)
{
	guide=src;
	double t_start=time_sec();
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	//double bpp=
	//	((double)src->src_depth[0]+src->src_depth[1]+src->src_depth[2]+src->src_depth[3])/
	//	(8*((src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0)));
	double usize=image_getBMPsize(src);
	int nch=src->nch;
	//int nch=get_nch32(src->data, res);//FIXME: differentiate between just gray and just alpha
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T46 SLIC4  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	Image *im2=0, *im3=0, *im4=0;
//#ifdef _DEBUG
//	Image *im5=0, *im6=0, *im7=0;//
//	image_copy(&im5, src);
//	image_copy(&im6, src);
//#endif
	image_copy(&im2, src);
	image_copy(&im3, src);
	image_copy(&im4, src);
	TempHybrid hu;
	hybriduint_encode(-maxlevels, &hu);//-(half<<1) to make way for chroma inflation
	int cdfsize=hu.token+1;
	unsigned *hist=(unsigned*)malloc((cdfsize+1)*sizeof(int[_countof(qlevels)+1]));
	unsigned short *emit_hist=(unsigned short*)malloc(cdfsize*sizeof(short));
	//unsigned short *emit_hist=(unsigned short*)malloc(cdfsize*sizeof(short[_countof(qlevels)+1]));
	//long long *sse=(long long*)malloc(SSE_SIZE*sizeof(long long));
	PredictorCtx pr;
	int success=slic4_pred_init(&pr, src->iw);
	if(!im2||!im3||!im4||!hist||!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	
	unsigned short *palettes[4]={0};
	unsigned short pal_sizes[4]={0};
#ifndef SLIC4_DISABLE_PALETTE
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
	rct_JPEG2000_32(im2, 1);
	unsigned char rct_depths[]=
	{
		src->depth[1],	//Y
		src->depth[2]+1,//Cb
		src->depth[0]+1,//Cr
		src->depth[3],	//a
	};
//#ifdef _DEBUG
//	image_copy(&im7, im2);
//	memcpy(im7->depth, rct_depths, sizeof(rct_depths));
//	pred_jxl(im7, 1, 1);
//#endif
	//memset(header.hist, 0, sizeof(header.hist));
	memset(hist, 0, (cdfsize+1)*sizeof(int[_countof(qlevels)+1]));
	//memset(emit_hist, 0, cdfsize*sizeof(short[_countof(qlevels)+1]));
	for(int kc=0;kc<nch;++kc)//for each channel
	{
		if(pal_sizes[kc]==1)
			continue;
		int nlevels=1<<rct_depths[kc];
		//int shift=(MAXVAR(8, rct_depths[kc])-8)>>2;
		slic4_pred_nextchannel(&pr, rct_depths[kc]);
		//memset(sse, 0, SSE_SIZE*sizeof(long long));
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			slic4_pred_nextrow(&pr, ky);
			//int prev_error=0;
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				//if(kc==0&&kx==114&&ky==2)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==4&&ky==0)//
				//if(kc==0&&kx==509&&ky==0)//
				//if(kc==0&&kx==23&&ky==0)//
				//if(kc==0&&kx==0&&ky==0)//
				//if(kc==1&&kx==648&&ky==354)//
				//if(kc==0&&kx==237&&ky==205)//
				//if(kc==0&&kx==108&&ky==2)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==2&&ky==0)//
				//if(kc==0&&kx==0&&ky==511)//
				//if(kc==0&&kx==1&&ky==0)//
				//if(kc==0&&kx==405&&ky==39)//
				//	printf("");

				slic4_predict(im2, kc, kx, ky, &pr);
				int curr=im2->data[idx<<2|kc];
				int error=curr-pr.pred_final;
				int negmask=-(pr.sse_corr<0);
				error^=negmask;
				error-=negmask;
				//error+=nlevels>>1;//no need for modular arithmetic
				//error&=nlevels-1;
				//error-=nlevels>>1;
				if(error)//sign pack from CALIC
				{
					int cond=(pr.pred_final<0)==(pr.sse_corr<0);
					int upred=pr.pred_final;
					upred^=-cond;
					upred+=cond;
					upred+=nlevels>>1;//upred = half + ((pred<0)==(corr<0) ? -pred : pred)
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
				hybriduint_encode(error, &hu);
				if(hu.token>=cdfsize)
					LOG_ERROR("Token OOB %d/%d", hu.token, cdfsize);
				++hist[(cdfsize+1)*pr.hist_idx+hu.token];
				im3->data[idx<<2|kc]=hu.token<<16|pr.hist_idx;
				im4->data[idx<<2|kc]=hu.bypass;
				slic4_pred_update(&pr, curr);
				//prev_error=error0;
//#ifdef _DEBUG
//				im5->data[idx<<2|kc]=((curr-pr.pred_final+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
//				im6->data[idx<<2|kc]=curr-((pr.preds[4]+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
//#endif
			}
		}
	}
	free(im2);
	if(loud)
	{
		printf("Hist WH %d*%d:\n", cdfsize, (int)_countof(qlevels)+1);
		for(int kt=0;kt<(int)_countof(qlevels)+1;++kt)
		{
			for(int ks=0;ks<cdfsize;++ks)
				printf("%d ", hist[(cdfsize+1)*kt+ks]);
			printf("\n");
		}
//#ifdef _DEBUG
//		double csizes[4]={0};
//		calc_csize(im7, csizes);
//		printf("WP0 TYUVA %lf %lf %lf %lf %lf  CR %lf\n",
//			csizes[0]+csizes[1]+csizes[2]+csizes[3],
//			csizes[0], csizes[1], csizes[2], csizes[3],
//			usize/(csizes[0]+csizes[1]+csizes[2]+csizes[3])
//		);
//		calc_csize(im5, csizes);
//		printf("WP  TYUVA %lf %lf %lf %lf %lf  CR %lf\n",
//			csizes[0]+csizes[1]+csizes[2]+csizes[3],
//			csizes[0], csizes[1], csizes[2], csizes[3],
//			usize/(csizes[0]+csizes[1]+csizes[2]+csizes[3])
//		);
//		calc_csize(im6, csizes);
//		printf("CG  TYUVA %lf %lf %lf %lf %lf  CR %lf\n",
//			csizes[0]+csizes[1]+csizes[2]+csizes[3],
//			csizes[0], csizes[1], csizes[2], csizes[3],
//			usize/(csizes[0]+csizes[1]+csizes[2]+csizes[3])
//		);
//#endif
		//printf("SSE:\n");
		//for(int k1=0;k1<_countof(qlevels)+1;++k1)
		//{
		//	for(int k2=0;k2<_countof(qlevels)+1;++k2)
		//	{
		//		long long *cell=sse+17*k1+k2;
		//		int count=(int)(*cell&0xFFF);
		//		long long sum=*cell>>12;
		//		int corr=count?(int)(sum/count):0;
		//		printf("%3d ", corr);
		//	}
		//	printf("\n");
		//}
	}
//#ifdef _DEBUG
//	free(im5);
//	free(im6);
//	free(im7);
//#endif
	slic4_pred_free(&pr);
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, pal_sizes, sizeof(pal_sizes));
	//free(sse);
	size_t overhead=0;
	for(int kt=0;kt<(int)_countof(qlevels)+1;++kt)
	{
		//if(kt==11)
		//	printf("");
		unsigned *curr_hist=hist+(cdfsize+1)*kt;
		//unsigned short *curr_emit_hist=emit_hist+cdfsize*kt;
		int weight=0, bins_used=0, last_symbol=0;
		for(int ks=0;ks<cdfsize;++ks)
		{
			unsigned freq=curr_hist[ks];
			weight+=freq;
			bins_used+=freq!=0;
			if(freq)
				last_symbol=ks;
		}
		unsigned short CDFrange[2];
		if(bins_used>1)
		{
			for(int ks=0, ks2=0, c=0;ks<cdfsize;++ks)
			{
				unsigned freq=curr_hist[ks];
				curr_hist[ks]=(unsigned)((unsigned long long)c*(0x10000LL-bins_used)/weight)+ks2;
				ks2+=freq!=0;
				c+=freq;
			}
			curr_hist[cdfsize]=0x10000;
			int start, end;
			for(start=0;start<cdfsize+1&&!curr_hist[start];++start);
			for(end=cdfsize;end>=start&&curr_hist[end]==0x10000;--end);
			++end;
			int n=end-start;
			CDFrange[0]=TAG_HIST_START<<12|(start&0xFFF);
			CDFrange[1]=TAG_HIST_LEN<<12|(n&0xFFF);
			dlist_push_back(&list, CDFrange, sizeof(CDFrange));
			if(n>0)
			{
				for(int k=0;k<n;++k)
					emit_hist[k]=(unsigned short)curr_hist[start+k];
				dlist_push_back(&list, emit_hist, n*sizeof(short));
			}
			overhead+=sizeof(CDFrange)+n*sizeof(short);
		}
		else if(bins_used)//degenerate histogram
		{
			//memset(curr_emit_hist, 0, (last_symbol+1)*sizeof(short));
			int c=0x10000;
			//memset(curr_emit_hist+last_symbol+1, -1, (cdfsize-(last_symbol+1))*sizeof(short));
			curr_hist[last_symbol]=0;
			memfill(curr_hist+last_symbol+1, &c, (cdfsize+1-(last_symbol+1))*sizeof(int), sizeof(int));
			CDFrange[0]=TAG_HIST_START<<12|(last_symbol&0xFFF);
			CDFrange[1]=TAG_HIST_LEN<<12|0;
			dlist_push_back(&list, CDFrange, sizeof(CDFrange));
			//for(int ks=last_symbol+1;ks<cdfsize;++ks)
			//	curr_emit_hist[ks]=curr_hist[ks]=0xFFFF;
			overhead+=sizeof(CDFrange);
		}
		else//null histogram, no encounters
		{
			CDFrange[0]=TAG_HIST_START<<12|0;
			CDFrange[1]=TAG_HIST_LEN<<12|0;
			dlist_push_back(&list, CDFrange, sizeof(CDFrange));
			//memset(curr_emit_hist, 0, cdfsize*sizeof(short));
			curr_hist[cdfsize]=0x10000;
			overhead+=sizeof(CDFrange);
		}
	}
	free(emit_hist);
	//dlist_push_back(&list, emit_hist, cdfsize*sizeof(short[_countof(qlevels)+1]));//TODO: histogram end marked with zero
	for(int kc=0;kc<4;++kc)//insert palettes, if any
	{
		if(palettes[kc])
		{
			dlist_push_back(&list, palettes[kc], pal_sizes[kc]*sizeof(short));
			free(palettes[kc]);
			palettes[kc]=0;
		}
	}
	ANSCoder ec;
	ans_enc_init(&ec, &list);
	for(int kc=nch-1;kc>=0;--kc)//for each channel
	{
		if(pal_sizes[kc]==1)//degenerate channel
			continue;
		for(ptrdiff_t k=res-1;k>=0;--k)
		{
			//if((k<<2|kc)==6600)//
			//if((k<<2|kc)==4)//
			//if((k<<2|kc)==16)//
			//if((k<<2|kc)==2306)//
			//if((k<<2|kc)==92)//
			//if((k<<2|kc)==0)//
			//if((k<<2|kc)==215340)//
			//	printf("");

			int token=im3->data[k<<2|kc]>>16, ctx=im3->data[k<<2|kc]&0xFFFF, bypass=im4->data[k<<2|kc];
#ifdef SLIC4_USE_SIMPLE_HYBRID
			if(token>=16)
#else
			if(token>=(1<<SLIC4_CONFIG_EXP))
#endif
			{
#ifdef SLIC4_USE_SIMPLE_HYBRID
				int nbits=((token-16)>>1)+3;
#else
				int nbits=((token-(1<<SLIC4_CONFIG_EXP))>>(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB))+SLIC4_CONFIG_EXP-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
#endif
				if(!nbits)
					LOG_ERROR("Bypass nbits 0");
				//unsigned char *ptr=(unsigned char*)&bypass;
				while(nbits>8)//encode bypass LSB-first
				{
					ans_enc(&ec, bypass&0xFF, 0, 256);
					nbits-=8;
					bypass>>=8;
				}
				//while(nbits>8)//encode bypass MSB-first
				//{
				//	nbits-=8;
				//	ans_enc(&ec, bypass>>nbits&0xFF, 0, 256);
				//}
				ans_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits);
			}
			ans_enc(&ec, token, hist+(cdfsize+1)*ctx, cdfsize);
		}
	}
	ans_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);
		printf("new overhead   %5d\n", (int)(sizeof(pal_sizes)+overhead));
		printf("naive overhead %5d\n", (int)(sizeof(pal_sizes)+cdfsize*sizeof(short[_countof(qlevels)+1])));
	}
	dlist_clear(&list);
	free(hist);
	free(im3);
	free(im4);
	return 1;
}
int t46_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	int res=dst->iw*dst->ih;
	int maxdepth=calc_maxdepth(dst, 0), maxlevels=1<<maxdepth;
	TempHybrid hu;
	hybriduint_encode(-maxlevels, &hu);//-(half<<1) to make way for chroma inflation
	int cdfsize=hu.token+1;
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

	unsigned *hist=(unsigned*)malloc((cdfsize+1LL)*sizeof(int[_countof(qlevels)+1]));
	//unsigned short *emit_hist=(unsigned short*)malloc(emithistsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int kt=0;kt<(int)_countof(qlevels)+1;++kt)
	{
		//if(kt==11)
		//	printf("");
		unsigned short CDFrange[2];
		if(srclen<sizeof(CDFrange))
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		memcpy(CDFrange, data, sizeof(CDFrange));
		data+=sizeof(CDFrange);
		srclen-=sizeof(CDFrange);

		if((CDFrange[0]>>12)!=TAG_HIST_START||(CDFrange[1]>>12)!=TAG_HIST_LEN)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		CDFrange[0]&=0xFFF;
		CDFrange[1]&=0xFFF;

		CDFrange[1]+=CDFrange[0];//{start, len} -> {start, end}

		if((unsigned)CDFrange[0]>(unsigned)cdfsize+1||(unsigned)CDFrange[1]>(unsigned)cdfsize+1||CDFrange[0]>CDFrange[1])
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}

		unsigned *curr_hist=hist+(cdfsize+1)*kt;
		//if(CDFrange[0]==CDFrange[1])//null histogram
		//{
		//	memset(curr_hist, 0, cdfsize*sizeof(int));
		//	curr_hist[cdfsize]=0x10000;
		//}
		if(CDFrange[0]==CDFrange[1])//degenerate histogram
		{
			memset(curr_hist, 0, (CDFrange[0]+1)*sizeof(int));
			if(CDFrange[0]+1<cdfsize+1)
			{
				int fillval=0x10000;
				memfill(curr_hist+CDFrange[0]+1, &fillval, (cdfsize+1-(CDFrange[0]+1))*sizeof(int), sizeof(int));
			}
		}
		else
		{
			if(CDFrange[0])
				memset(curr_hist, 0, CDFrange[0]*sizeof(int));
			for(int ks=CDFrange[0];ks<CDFrange[1];++ks)
			{
				unsigned short CDFval;
				if(srclen<sizeof(CDFval))
				{
					LOG_ERROR("Corrupt file");
					return 0;
				}
				memcpy(&CDFval, data, sizeof(CDFval));
				data+=sizeof(CDFval);
				srclen-=sizeof(CDFval);
				curr_hist[ks]=CDFval;
			}
			if(CDFrange[1]<cdfsize+1)
			{
				int fillval=0x10000;
				memfill(curr_hist+CDFrange[1], &fillval, (cdfsize+1-CDFrange[1])*sizeof(int), sizeof(int));
			}
		}
	}
	//memcpy(emit_hist, data, emithistsize);
	//data+=emithistsize;
	//srclen-=emithistsize;

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
			if(bytesize>(int)srclen)
			{
				LOG_ERROR("Invalid file\n");
				return 0;
			}
			memcpy(palettes[kc], data, bytesize);
			data+=bytesize;
			srclen-=bytesize;
		}
	}
	PredictorCtx pr;
	int success=slic4_pred_init(&pr, dst->iw);
	//long long *sse=(long long*)malloc(SSE_SIZE*sizeof(long long));
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	//for(int kt=0;kt<_countof(qlevels)+1;++kt)
	//{
	//	for(int ks=0, overflow=0;ks<cdfsize;++ks)
	//	{
	//		if(overflow)
	//			hist[(cdfsize+1)*kt+ks]=0x10000;
	//		else
	//		{
	//			unsigned cdf=emit_hist[cdfsize*kt+ks];
	//			hist[(cdfsize+1)*kt+ks]=cdf;
	//			if(ks<cdfsize-1)
	//				overflow|=cdf>emit_hist[cdfsize*kt+ks+1];
	//		}
	//	}
	//	hist[(cdfsize+1)*kt+cdfsize]=0x10000;
	//}
	//free(emit_hist);
	ANSCoder ec;
	ans_dec_init(&ec, data, data+srclen);
	int rct_depths[]=
	{
		dst->depth[1],//Y
		dst->depth[2]+1,//Cb
		dst->depth[0]+1,//Cr
		dst->depth[3],//a
	};
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(pal_sizes[kc]==1)
		{
			for(int k=0;k<res;++k)
				dst->data[k<<2|kc]=palettes[kc][0];
			continue;
		}
		int nlevels=1<<rct_depths[kc];
		//int shift=(MAXVAR(8, rct_depths[kc])-8)>>2;
		slic4_pred_nextchannel(&pr, rct_depths[kc]);
		//memset(sse, 0, SSE_SIZE*sizeof(long long));
		for(int ky=0, idx=0;ky<dst->ih;++ky)
		{
			slic4_pred_nextrow(&pr, ky);
			//int prev_error=0;
			for(int kx=0;kx<dst->iw;++kx, ++idx)
			{
				//if(kc==0&&kx==114&&ky==2)//
				//if(kc==0&&kx==1&&ky==0)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==4&&ky==0)//
				//if(kc==0&&kx==509&&ky==0)//
				//if(kc==0&&kx==0&&ky==0)//
				//if(kc==0&&kx==23&&ky==0)//
				//if(kc==0&&kx==0&&ky==0)//
				//if(kc==1&&kx==648&&ky==354)//
				//if(kc==0&&kx==108&&ky==2)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==2&&ky==0)//
				//if(kc==0&&kx==0&&ky==511)//
				//if(kc==0&&kx==405&&ky==39)//
				//	printf("");

				slic4_predict(dst, kc, kx, ky, &pr);
				unsigned *CDF=hist+(cdfsize+1)*pr.hist_idx;
				int sym=ans_dec(&ec, CDF, cdfsize);
#ifdef SLIC4_USE_SIMPLE_HYBRID
				if(sym>=16)
				{
					int nbits=((sym-16)>>1)+3;
					int bypass=ans_dec(&ec, 0, 1<<nbits);
					sym=1<<(nbits+1)|(sym&1)<<nbits|bypass;
				}
#else
				if(sym>=(1<<SLIC4_CONFIG_EXP))
				{
					sym-=1<<SLIC4_CONFIG_EXP;
					int lsb=sym&((1<<SLIC4_CONFIG_LSB)-1);
					sym>>=SLIC4_CONFIG_LSB;
					int msb=sym&((1<<SLIC4_CONFIG_MSB)-1);
					sym>>=SLIC4_CONFIG_MSB;
					int nbits=sym+SLIC4_CONFIG_EXP-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
					int bypass=nbits&7?ans_dec(&ec, 0, 1<<(nbits&7)):0;
					for(int k=0, nb=nbits>>3;k<nb;++k)
						bypass=bypass<<8|ans_dec(&ec, 0, 256);
					sym=1;
					sym<<=SLIC4_CONFIG_MSB;
					sym|=msb;
					sym<<=nbits;
					sym|=bypass;
					sym<<=SLIC4_CONFIG_LSB;
					sym|=lsb;
				}
#endif
				int error;
				int cond=(pr.pred_final<0)==(pr.sse_corr<0);
				int upred=pr.pred_final;
				upred^=-cond;
				upred+=cond;
				upred+=nlevels>>1;
				if(sym<=(upred<<1))//sign pack from CALIC
					error=sym>>1^-(sym&1);//unpack sign
				else
				{
					error=sym-upred;//error sign is known
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
				int negmask=-(pr.sse_corr<0);
				error^=negmask;
				error-=negmask;
				int curr=error+pr.pred_final;
				//curr+=nlevels>>1;//no need for modular arithmetic
				//curr&=nlevels-1;
				//curr-=nlevels>>1;
				//if(guide&&curr!=guide->data[idx<<2|kc])
				//	LOG_ERROR("Guide error CXY %d %d %d", kc, kx, ky);
				dst->data[idx<<2|kc]=curr;
				slic4_pred_update(&pr, curr);
				//prev_error=error;
			}
		}
	}
	//free(sse);
	free(hist);
	slic4_pred_free(&pr);
	rct_JPEG2000_32(dst, 0);
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
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}