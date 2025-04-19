#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define PRINT_MAXERROR

#if 0
#define WP_NPREDS 8
#define WP_PREDLIST\
	WP_PRED(W+NE-N)\
	WP_PRED(N-((eN+eW+eNE)*params[4]>>8))\
	WP_PRED(W-((eN+eW+eNW)*params[5]>>8))\
	WP_PRED(N-((eNW*params[6]+eN*params[7]+eNE*params[8]+(NN-N)*params[9]+(NW-W)*params[10])>>8))\
	WP_PRED(N)\
	WP_PRED(W)\
	WP_PRED(3*(N-NN)+NNN)\
	WP_PRED(3*(W-WW)+WWW)
#endif
#if 1
#define WP_NPREDS 4
#define WP_PREDLIST\
	WP_PRED(W+NE-N)\
	WP_PRED(N-((eN+eW+eNE)*params[4]>>8))\
	WP_PRED(W-((eN+eW+eNW)*params[5]>>8))\
	WP_PRED(N-((eNW*params[6]+eN*params[7]+eNE*params[8]+(NN-N)*params[9]+(NW-W)*params[10])>>8))
#endif

static const short allparams[33]=//signed fixed 7.8 bit
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,//Y
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,//Cb
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,//Cr
};

//JXL WP exact
#if 0
#define WP_NPREDS 4
#define WP_PREDLIST\
	WP_PRED(W+NE-N)\
	WP_PRED(N-((eN+eW+eNE)*p1C>>5))\
	WP_PRED(W-((eN+eW+eNW)*p2C>>5))\
	WP_PRED(N-((eNW*p3Ca+eN*p3Cb+eNE*p3Cc+(NN-N)*p3Cd+(NW-W)*p3Ce)>>5))
#endif

#if 0
#define WP_NPREDS 8
#define WP_PREDLIST\
	WP_PRED(230,	N+eN/4)\
	WP_PRED(230,	W+eW/4)\
	WP_PRED(190,	3*(N-NN)+NNN)\
	WP_PRED(190,	3*(W-WW)+WWW)\
	WP_PRED(170,	W+NE-N)\
	WP_PRED(170,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WP_PRED(120,	N+W-NW)\
	WP_PRED(120,	N+NE-NNE)
#endif
#if 0
#define WP_NPREDS 8
#define WP_PREDLIST\
	WP_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WP_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WP_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WP_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WP_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WP_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WP_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WP_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#if 0
	WP_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-eN-eNN+eNNN)/3)\
	WP_PRED(97,	(6*(W+WWW)+20*WW+(eW-eWW+eWWW)/3)/32)\
	WP_PRED(65,	(W+3*NW-NWW-NNWW)/2+eNW/4+eW/6)\
	WP_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#endif
#endif

#define NBITS 3

#ifdef PRINT_MAXERROR
static unsigned maxerror=0;
#endif
FORCE_INLINE int wp_predict(
	const int *errors, const unsigned short *cN, const unsigned short *ccurr,
	const short *NNNptr, const short *NNptr, const short *Nptr, const short *currptr,
	const short *params, int *preds
)
{
	int
		NNNN	=currptr[0+0*4*2],
		NNNWWWW	=NNNptr	[0-4*4*2],
		NNNWWW	=NNNptr	[0-3*4*2],
		NNNW	=NNNptr	[0-1*4*2],
		NNN	=NNNptr	[0+0*4*2],
		NNNE	=NNNptr	[0+1*4*2],
		NNNEEE	=NNNptr	[0+3*4*2],
		NNWWWW	=NNptr	[0-4*4*2],
		NNWW	=NNptr	[0-2*4*2],
		NNW	=NNptr	[0-1*4*2],
		NN	=NNptr	[0+0*4*2],
		NNE	=NNptr	[0+1*4*2],
		NNEE	=NNptr	[0+2*4*2],
		NNEEE	=NNptr	[0+3*4*2],
		NNEEEE	=NNptr	[0+4*4*2],
		NWWW	=Nptr	[0-3*4*2],
		NWW	=Nptr	[0-2*4*2],
		NW	=Nptr	[0-1*4*2],
		N	=Nptr	[0+0*4*2],
		NE	=Nptr	[0+1*4*2],
		NEE	=Nptr	[0+2*4*2],
		NEEE	=Nptr	[0+3*4*2],
		NEEEE	=Nptr	[0+4*4*2],
		WWWWW	=currptr[0-5*4*2],
		WWWW	=currptr[0-4*4*2],
		WWW	=currptr[0-3*4*2],
		WW	=currptr[0-2*4*2],
		W	=currptr[0-1*4*2],
		eNNN	=NNNptr	[4+0*4*2],
		eNN	=NNptr	[4+0*4*2],
		eNNE	=NNptr	[4+1*4*2],
		eNW	=Nptr	[4-1*4*2],
		eN	=Nptr	[4+0*4*2],
		eNE	=Nptr	[4+1*4*2],
		eNEE	=Nptr	[4+2*4*2],
		eNEEE	=Nptr	[4+3*4*2],
		eWWWW	=currptr[4-4*4*2],
		eWWW	=currptr[4-3*4*2],
		eWW	=currptr[4-2*4*2],
		eW	=currptr[4-1*4*2];

	int j=0;
	
#define WP_PRED(EXPR) preds[j++]=EXPR;
	WP_PREDLIST
#undef  WP_PRED
#if 0
	unsigned weights[WP_NPREDS];
	//		[1]	[1]	[1]
	//	 1	 1	 ?		+I*4		9*nlevels
	for(int k=0;k<WP_NPREDS;++k)
		weights[k]=(errors[k]>>4)+cN[k-1*4*WP_NPREDS]+cN[k+0*4*WP_NPREDS]+cN[k+1*4*WP_NPREDS]+2;
	for(int k=0;k<4;++k)
		weights[k]=(params[k]<<8)/weights[k];
	for(int k=4;k<WP_NPREDS;++k)
		weights[k]=(256<<8)/weights[k];
	int wsum=0;
	for(int k=0;k<WP_NPREDS;++k)
		wsum+=weights[k];
	long long ipred=0;
	for(int k=0;k<WP_NPREDS;++k)
		ipred+=(long long)weights[k]*preds[k];
	int pred=(int)(ipred/wsum);
#endif
#if 1
	unsigned weights[WP_NPREDS];
	//		[1]	[1]	[1]
	//	 1	 1	 ?		+I*4		9*nlevels
	weights[0]=(errors[0]>>4)+cN[0-1*4*WP_NPREDS]+cN[0+0*4*WP_NPREDS]+cN[0+1*4*WP_NPREDS]+2;
	weights[1]=(errors[1]>>4)+cN[1-1*4*WP_NPREDS]+cN[1+0*4*WP_NPREDS]+cN[1+1*4*WP_NPREDS]+2;
	weights[2]=(errors[2]>>4)+cN[2-1*4*WP_NPREDS]+cN[2+0*4*WP_NPREDS]+cN[2+1*4*WP_NPREDS]+2;
	weights[3]=(errors[3]>>4)+cN[3-1*4*WP_NPREDS]+cN[3+0*4*WP_NPREDS]+cN[3+1*4*WP_NPREDS]+2;
	//unsigned best=weights[0];
	//if(best>weights[1])best=weights[1];
	//if(best>weights[2])best=weights[2];
	//if(best>weights[3])best=weights[3];
	//best>>=4;
	weights[0]=(params[0]<<8)/weights[0];
	weights[1]=(params[1]<<8)/weights[1];
	weights[2]=(params[2]<<8)/weights[2];
	weights[3]=(params[3]<<8)/weights[3];
	int wsum=weights[0]+weights[1]+weights[2]+weights[3];
	int ipred=weights[0]*preds[0]+weights[1]*preds[1]+weights[2]*preds[2]+weights[3]*preds[3];
	int pred=ipred/wsum;
#endif
#if 0
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[DIVLUTSIZE]={0};
	if(!*divlookup)
	{
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	unsigned coeff[WP_NPREDS], wsum=0;
	int ipred=0;
	int pred=0;
	for(int kp=0;kp<WP_NPREDS;++kp)
	{
		//		[1]	[1]	[1]
		//	 1	 1	 ?
		coeff[kp]=//maxerror = 18<<depth
			+(errors[kp]>>4)		//+I	*8
			+cN	[kp-1*4*WP_NPREDS]	//+eNW+eWW
			+cN	[kp+0*4*WP_NPREDS]	//+eN +eW
			+cN	[kp+1*4*WP_NPREDS]	//+eNE
		;

		//			1	1
		//		2	2	1
		//	1	2	?
		//coeff[kp]=//maxerror = 18<<depth
		//	+(errors[kp]>>3)		//+I	*8
		//	+cN	[kp-1*4*WP_NPREDS]*2	//+eNW	*2
		//	+cN	[kp+0*4*WP_NPREDS]*2	//+eN	*2
		//	+ccurr	[kp-1*4*WP_NPREDS]*2	//+eW	*2
		//	+ccurr	[kp+0*4*WP_NPREDS]	//+eNN	*1
		//	+ccurr	[kp+1*4*WP_NPREDS]	//+eNNE	*1
		//	+cN	[kp+1*4*WP_NPREDS]	//+eNE	*1
		//	+ccurr	[kp-2*4*WP_NPREDS]	//+eWW	*1
		//;
	}
	int sh[WP_NPREDS];
	sh[0]=31-(DENBITS-1)-_lzcnt_u32(coeff[0]+1);//invert errros
	sh[1]=31-(DENBITS-1)-_lzcnt_u32(coeff[1]+1);
	sh[2]=31-(DENBITS-1)-_lzcnt_u32(coeff[2]+1);
	sh[3]=31-(DENBITS-1)-_lzcnt_u32(coeff[3]+1);
	if(sh[0]<0)sh[0]=0;
	if(sh[1]<0)sh[1]=0;
	if(sh[2]<0)sh[2]=0;
	if(sh[3]<0)sh[3]=0;
	//libjxl/context_predict.h
//	const int w0=13, w1=12, w2=12, w3=12,	p1C=16, p2C=10, p3Ca=7, p3Cb=7, p3Cc=7, p3Cd=0, p3Ce=0;//lossless16
//	const int w0=13, w1=12, w2=12, w3=11,	p1C=8, p2C=8, p3Ca=4, p3Cb=0, p3Cc=3, p3Cd=23, p3Ce=2;//default lossless8
//	const int w0=13, w1=12, w2=13, w3=12,	p1C=10, p2C=9, p3Ca=7, p3Cb=0, p3Cc=0, p3Cd=16, p3Ce=9;//west lossless8
//	const int w0=13, w1=13, w2=12, w3=12,	p1C=16, p2C=8, p3Ca=0, p3Cb=16, p3Cc=0, p3Cd=23, p3Ce=0;//north lossless8
//	const int w0=13, w1=12, w2=12, w3=12,	p1C=10, p2C=10, p3Ca=5, p3Cb=5, p3Cc=5, p3Cd=12, p3Ce=4;//something else
	coeff[0]=(params[0]*divlookup[coeff[0]>>sh[0]]>>sh[0])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[1]=(params[1]*divlookup[coeff[1]>>sh[1]]>>sh[1])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[2]=(params[2]*divlookup[coeff[2]>>sh[2]]>>sh[2])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[3]=(params[3]*divlookup[coeff[3]>>sh[3]]>>sh[3])+(1<<DENBITS>>2)/WP_NPREDS;
	wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3];//invert coeff sum
	sh[0]=FLOOR_LOG2(wsum)-(DENBITS-2);
	wsum=0;
	coeff[0]>>=sh[0];
	coeff[1]>>=sh[0];
	coeff[2]>>=sh[0];
	coeff[3]>>=sh[0];
	wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3];
	ipred=(wsum>>1)
		+coeff[0]*preds[0]
		+coeff[1]*preds[1]
		+coeff[2]*preds[2]
		+coeff[3]*preds[3]
	;
	pred=ipred*divlookup[wsum-1]>>NUMBITS;
#endif
#if 0
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[DIVLUTSIZE]={0};
	if(!*divlookup)
	{
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	unsigned coeff[WP_NPREDS], wsum=0;
	int ipred=0;
	int pred=0;
	for(int kp=0;kp<WP_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		coeff[kp]=//maxerror = 18<<depth
			+(errors[kp]>>3)		//+I	*8
			+cN	[kp-1*4*WP_NPREDS]*2	//+eNW	*2
			+cN	[kp+0*4*WP_NPREDS]*2	//+eN	*2
			+ccurr	[kp-1*4*WP_NPREDS]*2	//+eW	*2
			+ccurr	[kp+0*4*WP_NPREDS]	//+eNN	*1
			+ccurr	[kp+1*4*WP_NPREDS]	//+eNNE	*1
			+cN	[kp+1*4*WP_NPREDS]	//+eNE	*1
			+ccurr	[kp-2*4*WP_NPREDS]	//+eWW	*1
		;
	}
	int sh[WP_NPREDS];
	sh[0]=31-(DENBITS-1)-_lzcnt_u32(coeff[0]+1);//invert errros
	sh[1]=31-(DENBITS-1)-_lzcnt_u32(coeff[1]+1);
	sh[2]=31-(DENBITS-1)-_lzcnt_u32(coeff[2]+1);
	sh[3]=31-(DENBITS-1)-_lzcnt_u32(coeff[3]+1);
	sh[4]=31-(DENBITS-1)-_lzcnt_u32(coeff[4]+1);
	sh[5]=31-(DENBITS-1)-_lzcnt_u32(coeff[5]+1);
	sh[6]=31-(DENBITS-1)-_lzcnt_u32(coeff[6]+1);
	sh[7]=31-(DENBITS-1)-_lzcnt_u32(coeff[7]+1);
	if(sh[0]<0)sh[0]=0;
	if(sh[1]<0)sh[1]=0;
	if(sh[2]<0)sh[2]=0;
	if(sh[3]<0)sh[3]=0;
	if(sh[4]<0)sh[4]=0;
	if(sh[5]<0)sh[5]=0;
	if(sh[6]<0)sh[6]=0;
	if(sh[7]<0)sh[7]=0;
	coeff[0]=(wp_weights[0]*divlookup[coeff[0]>>sh[0]]>>sh[0])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[1]=(wp_weights[1]*divlookup[coeff[1]>>sh[1]]>>sh[1])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[2]=(wp_weights[2]*divlookup[coeff[2]>>sh[2]]>>sh[2])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[3]=(wp_weights[3]*divlookup[coeff[3]>>sh[3]]>>sh[3])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[4]=(wp_weights[4]*divlookup[coeff[4]>>sh[4]]>>sh[4])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[5]=(wp_weights[5]*divlookup[coeff[5]>>sh[5]]>>sh[5])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[6]=(wp_weights[6]*divlookup[coeff[6]>>sh[6]]>>sh[6])+(1<<DENBITS>>2)/WP_NPREDS;
	coeff[7]=(wp_weights[7]*divlookup[coeff[7]>>sh[7]]>>sh[7])+(1<<DENBITS>>2)/WP_NPREDS;
	wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];//invert coeff sum
	sh[0]=FLOOR_LOG2(wsum)-(DENBITS-2);
	wsum=0;
	coeff[0]>>=sh[0];
	coeff[1]>>=sh[0];
	coeff[2]>>=sh[0];
	coeff[3]>>=sh[0];
	coeff[4]>>=sh[0];
	coeff[5]>>=sh[0];
	coeff[6]>>=sh[0];
	coeff[7]>>=sh[0];
	wsum=coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]+coeff[6]+coeff[7];
	ipred=(wsum>>1)
		+coeff[0]*preds[0]
		+coeff[1]*preds[1]
		+coeff[2]*preds[2]
		+coeff[3]*preds[3]
		+coeff[4]*preds[4]
		+coeff[5]*preds[5]
		+coeff[6]*preds[6]
		+coeff[7]*preds[7]
	;
	pred=ipred*divlookup[wsum-1]>>NUMBITS;
#endif
	//if(((eN^eW)|(eN^eNW))<=0)
	{
		int vmax=N, vmin=W;
		if(N<W)vmin=N, vmax=W;
		if(vmin>NE)vmin=NE;
		if(vmax<NE)vmax=NE;
		if(vmin>NEEE)vmin=NEEE;
		if(vmax<NEEE)vmax=NEEE;
		CLAMP2(pred, vmin, vmax);
	}
	return pred;
}
FORCE_INLINE void wp_update(
	int curr, const int *preds,
	unsigned short *cN,
	unsigned short *ccurr,
	int *errors
)
{
	int e2[WP_NPREDS], best=0x7FFFFFFF;
	for(int kp=0;kp<WP_NPREDS;++kp)
	{
		int e=(abs(curr-preds[kp])+(1<<NBITS>>1))>>NBITS;
		if(best>e)
			best=e;
		e2[kp]=e;
	}
	for(int kp=0;kp<WP_NPREDS;++kp)
	{
		int e=e2[kp]-best;
		errors[kp]+=((e<<6)-errors[kp]+(1<<3>>1))>>3;
		ccurr	[kp+0*4*WP_NPREDS]=e;
		cN	[kp+1*4*WP_NPREDS]+=e;
	}
}
void pred_wpred7(Image *src, int fwd)
{
	ALIGN(32) int wp_errors[4][WP_NPREDS]={0}, wp_preds[WP_NPREDS]={0};

	int fwdmask=-fwd;
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	int pesize=(src->iw+16)*(int)sizeof(short[2][4][WP_NPREDS]);//2 padded rows * 4 channels * WG_NPREDS
	unsigned short *ebuf=(unsigned short*)_mm_malloc(pesize, sizeof(__m256i));
	int bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	short *pixels=(short*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int cpred=0, curr=0;
		short *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		unsigned short *erows[]=
		{
			ebuf+((src->iw+16LL)*((ky-0LL)&1)+8LL)*4*WP_NPREDS,
			ebuf+((src->iw+16LL)*((ky-1LL)&1)+8LL)*4*WP_NPREDS,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			//if(ky==src->ih/2&&kx==src->iw/2)//
			//	printf("");//

			for(int kc=0;kc<src->nch;++kc)
			{
				int pred=wp_predict(
					wp_errors[kc],
					erows[1]+(ptrdiff_t)WP_NPREDS*kc,
					erows[0]+(ptrdiff_t)WP_NPREDS*kc,
					rows[3]+kc,
					rows[2]+kc,
					rows[1]+kc,
					rows[0]+kc,
					allparams+11*kc,
					wp_preds
				);
				curr=src->data[idx+kc];

				cpred=(pred+(1<<NBITS>>1))>>NBITS;
				cpred^=fwdmask;
				cpred-=fwdmask;
				cpred+=curr;
				cpred<<=32-src->depth[kc];
				cpred>>=32-src->depth[kc];

				src->data[idx+kc]=cpred;
				if(!fwd)
					curr=cpred;
				rows[0][kc+0]=curr<<NBITS;
				rows[0][kc+4]=rows[0][kc+0]-pred;
				wp_update(
					rows[0][kc+0],
					wp_preds,
					erows[1]+(ptrdiff_t)WP_NPREDS*kc,
					erows[0]+(ptrdiff_t)WP_NPREDS*kc,
					wp_errors[kc]
				);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
			erows[0]+=4*WP_NPREDS;
			erows[1]+=4*WP_NPREDS;
		}
	}
#ifdef PRINT_MAXERROR
	if(loud_transforms)
		messagebox(MBOX_OK, "Info", "max error %d", maxerror);
#endif
	free(pixels);
	_mm_free(ebuf);
}