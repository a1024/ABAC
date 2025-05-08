#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define PRINT_MAXCOEFF
//	#define PRINT_MAXERROR


#if 0
#define WG_NPREDS 8
#define WG_PREDLIST\
	WG_PRED(230,	N+eN/4)\
	WG_PRED(230,	W+eW/4)\
	WG_PRED(190,	3*(N-NN)+NNN)\
	WG_PRED(190,	3*(W-WW)+WWW)\
	WG_PRED(170,	W+NE-N)\
	WG_PRED(170,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#endif
#if 0
#define WG_NPREDS 7
#define WG_PREDLIST\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(180,	3*(N-NN)+NNN)\
	WG_PRED(180,	3*(W-WW)+WWW)\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(120,	N+NE-NNE)
//	WG_PRED(120,	(2*W-WW+2*NE-(NEE+NNE)/2)/2)
#endif
#if 0
#define WG_NPREDS 6
#define WG_PREDLIST\
	WG_PRED(120,	N)\
	WG_PRED(180,	W)\
	WG_PRED(180,	NW)\
	WG_PRED(240,	NE)\
	WG_PRED(240,	NN)\
	WG_PRED(140,	WW)
#endif
#if 0
#define WG_NPREDS 11
#define WG_PREDLIST\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	3*(N-NN)+NNN)\
	WG_PRED(180,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(120,	N+NE-NNE)\
	WG_PRED(180,	N+NW-NNW)\
	WG_PRED(180,	W+NW-NWW)\
	WG_PRED(160,	NE+NEE-NNEEE)\
	WG_PRED(160,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)
#endif
#if 0
#define WG_NPREDS 2
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)
#endif
#if 0
#define WG_NPREDS 4
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	(3*(N-NN)+NNN+3*(W-WW)+WWW)>>1)\
	WG_PRED(180,	N+W-NW)
#endif
#if 0
#define WG_NPREDS 8
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	3*(N-NN)+NNN)\
	WG_PRED(180,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(160,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#endif
#if 0
#define WG_NPREDS 12
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	2*N-NN)\
	WG_PRED(180,	2*W-WW)\
	WG_PRED(180,	3*(N-NN)+NNN)\
	WG_PRED(180,	3*(W-WW)+WWW)\
	WG_PRED(180,	4*(N+NNN)-6*NN-NNNN)\
	WG_PRED(180,	4*(W+WWW)-6*WW-WWWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(160,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#endif
#if 0
#define WG_NPREDS 8
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	3*(N-NN)-NNN)\
	WG_PRED(180,	3*(W-WW)-WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(160,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#endif
#if 1
#define WG_NPREDS 8
#define WG_PREDLIST\
	WG_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WG_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WG_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WG_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WG_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WG_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WG_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#if 0
	WG_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-eN-eNN+eNNN)/3)\
	WG_PRED(97,	(6*(W+WWW)+20*WW+(eW-eWW+eWWW)/3)/32)\
	WG_PRED(65,	(W+3*NW-NWW-NNWW)/2+eNW/4+eW/6)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#endif
#endif
#if 0
#define WG_NPREDS	12	//multiple of 4
#define WG_PREDLIST0\
	WG_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WG_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WG_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WG_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WG_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WG_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WG_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))\
	WG_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-eN-eNN+eNNN)/3)\
	WG_PRED(97,	(6*(W+WWW)+20*WW+(eW-eWW+eWWW)/3)/32)\
	WG_PRED(65,	(W+3*NW-NWW-NNWW)/2+eNW/4+eW/6)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#define WG_PREDLIST1\
	WG_PRED(250,	N+(3*eN+eNW+eW)/6)\
	WG_PRED(250,	W+(3*eW+eNW+eN)/6)\
	WG_PRED(175,	3*(N-NN)+NNN+(eN/2+eNN/2-eWW)/2)\
	WG_PRED(175,	3*(W-WW)+WWW+(eW/2+eWW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+((eN+eNE+eNNE+8)>>4))\
	WG_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-2*eN-eNN+eNNN)/3)\
	WG_PRED(57,	(W+WW+(eW-eWW+eWWW)/4)/2)\
	WG_PRED(35,	(W+3*NW-NWW-NNWW+eNW)/2+eW/3)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#define WG_PREDLIST2\
	WG_PRED(270,	N+(3*eN+eW)/6)\
	WG_PRED(270,	W+(3*eW+eN)/6)\
	WG_PRED(200,	3*(N-NN)+NNN+(eN/2-eWW)/2)\
	WG_PRED(200,	3*(W-WW)+WWW+(eW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((2*eN+eW+31)>>5))\
	WG_PRED(165,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+(eN+eNE+eNNE)/16)\
	WG_PRED(55,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW)/3-eN-eNN+eNNN)\
	WG_PRED(47,	(W+WW+(eW+eWWW)/3)/2)\
	WG_PRED(22,	(W+3*NW-NWW-NNWW+eNW)/2+eW/3)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#endif

//#define WG_NCTX 1	//power of 2
#define WG_NBITS 1	//why 0 is best?

#ifdef PRINT_MAXCOEFF
static long long maxcoeff=0;
#endif
#ifdef PRINT_MAXERROR
static unsigned maxerror=0;
#endif
static const int wg_weights[]=
{
#define WG_PRED(WEIGHT, EXPR) WEIGHT,
	WG_PREDLIST
	//WG_PREDLIST0
	//WG_PREDLIST1
	//WG_PREDLIST2
#undef  WG_PRED
};
//static void wg_init(int *weights)
//{
//	int j=0;
//#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
//	WG_PREDLIST
//#undef  WG_PRED
//}
static int wg_predict(
	const int *errors,
	const int *coeffs, const unsigned short *cN, const unsigned short *ccurr,
	const short *NNNptr, const short *NNptr, const short *Nptr, const short *currptr,
	int *subpreds, const int *mixer, const int *sse,
	int *preds, int depth, int kc
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
	
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
	//if(!kc)
	//{
	//	WG_PREDLIST0
	//}
	//else if(kc==1)
	//{
	//	WG_PREDLIST1
	//}
	//else if(kc==2)
	//{
	//	WG_PREDLIST2
	//}
	//else
	//	return W;
#undef  WG_PRED
#if 1
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[DIVLUTSIZE]={0};
	if(!*divlookup)
	{
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	unsigned coeff[WG_NPREDS], wsum=0;
	int ipred=0;
	int pred=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		coeff[kp]=//maxerror = 18<<depth
			+(errors[kp]>>3)		//+I	*8
			+cN	[kp-1*4*WG_NPREDS]*2	//+eNW	*2
			+cN	[kp+0*4*WG_NPREDS]*2	//+eN	*2
			+ccurr	[kp-1*4*WG_NPREDS]*2	//+eW	*2
			+ccurr	[kp+0*4*WG_NPREDS]	//+eNN	*1
			+ccurr	[kp+1*4*WG_NPREDS]	//+eNNE	*1
			+cN	[kp+1*4*WG_NPREDS]	//+eNE	*1
			+ccurr	[kp-2*4*WG_NPREDS]	//+eWW	*1
		;
	}
	int sh[WG_NPREDS];
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
	coeff[0]=(wg_weights[0]*divlookup[coeff[0]>>sh[0]]>>sh[0])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[1]=(wg_weights[1]*divlookup[coeff[1]>>sh[1]]>>sh[1])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[2]=(wg_weights[2]*divlookup[coeff[2]>>sh[2]]>>sh[2])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[3]=(wg_weights[3]*divlookup[coeff[3]>>sh[3]]>>sh[3])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[4]=(wg_weights[4]*divlookup[coeff[4]>>sh[4]]>>sh[4])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[5]=(wg_weights[5]*divlookup[coeff[5]>>sh[5]]>>sh[5])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[6]=(wg_weights[6]*divlookup[coeff[6]>>sh[6]]>>sh[6])+(1<<DENBITS>>2)/WG_NPREDS;
	coeff[7]=(wg_weights[7]*divlookup[coeff[7]>>sh[7]]>>sh[7])+(1<<DENBITS>>2)/WG_NPREDS;
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
	//static const int divlookup[]=
	//{
	//	//(1<<24)/(i+1)		libjxl/lib/jxl/modular/encode/context_predict.h
	//	16777216,  8388608, 5592405, 4194304, 3355443, 2796202, 2396745, 2097152,
	//	 1864135,  1677721, 1525201, 1398101, 1290555, 1198372, 1118481, 1048576,
	//	  986895,   932067,  883011,  838860,  798915,  762600,  729444,  699050,
	//	  671088,   645277,  621378,  599186,  578524,  559240,  541200,  524288,
	//	  508400,   493447,  479349,  466033,  453438,  441505,  430185,  419430,
	//	  409200,   399457,  390167,  381300,  372827,  364722,  356962,  349525,
	//	  342392,   335544,  328965,  322638,  316551,  310689,  305040,  299593,
	//	  294337,   289262,  284359,  279620,  275036,  270600,  266305,  262144,
	//};
	unsigned coeff[WG_NPREDS], wsum=0, sh=0;
	int ipred=0;
	int pred=0;
	//const int *weights=wg_weights+WG_NPREDS*kc;
	const int *weights=wg_weights;
//again:
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		unsigned e=//maxerror = 18<<depth
			+(errors[kp]>>2)		//+I	*8
			+cN	[kp-1*4*WG_NPREDS]*2	//+eNW	*2
			+cN	[kp+0*4*WG_NPREDS]*2	//+eN	*2
			+ccurr	[kp-1*4*WG_NPREDS]*2	//+eW	*2
			+ccurr	[kp+0*4*WG_NPREDS]	//+eNN	*1
			+ccurr	[kp+1*4*WG_NPREDS]	//+eNNE	*1
			+cN	[kp+1*4*WG_NPREDS]	//+eNE	*1
			+ccurr	[kp-2*4*WG_NPREDS]	//+eWW	*1
		;
		//unsigned e=
		//	+(errors[kp]>>9)
		//	+cN	[kp-1*4*WG_NPREDS]	//+eNW+eWW
		//	+cN	[kp+0*4*WG_NPREDS]	//+eN+eW
		//	+cN	[kp+1*4*WG_NPREDS]	//+eNE
		//	+ccurr	[kp-1*4*WG_NPREDS]	//+eW
		//	+ccurr	[kp+0*4*WG_NPREDS]	//+eNN+eNW
		//	+ccurr	[kp+1*4*WG_NPREDS]	//+eNNE+eN
		//;
		//if(e>5000)
		//	printf("");
#ifdef PRINT_MAXERROR
		if(maxerror<e)
			maxerror=e;
#endif
		sh=FLOOR_LOG2(e+1);
		coeff[kp]=(wg_weights[kp]*divlookup[e<<(DENBITS-1)>>sh]<<(DENBITS-1)>>sh)+(1<<DENBITS>>2)/WG_NPREDS;
		//sh=FLOOR_LOG2(e+1)-5;
		//if(sh<0)sh=0;
		//coeff[kp]=(weights[kp]*divlookup[e>>sh]>>sh)+(1<<DENBITS>>2)/WG_NPREDS;
		wsum+=coeff[kp];
	}
	//if(wsum<16)
	//	LOG_ERROR("wsum %d < 16", wsum);
	sh=FLOOR_LOG2(wsum)-(DENBITS-2);
	wsum=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int c=(int)(coeff[kp]>>sh);
		//if((unsigned)c>64)
		//	goto again;
		wsum+=c;
		coeff[kp]=c;
	}
	ipred=wsum>>1;
	//ipred=(int)((wsum>>1)-1LL);
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		ipred+=(int)(coeff[kp]*preds[kp]);
#ifdef PRINT_MAXCOEFF
		if(maxcoeff<coeff[kp])
			maxcoeff=coeff[kp];
#endif
	}
	//if(wsum<0)
	//	goto again;
	pred=ipred*divlookup[wsum-1]>>NUMBITS;
#endif
#if 0
#define NUMBITS 20
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[DIVLUTSIZE]={0};
	if(!*divlookup)
	{
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	//static const int divlookup[]=
	//{
	//	//(1<<24)/(i+1)		libjxl/lib/jxl/modular/encode/context_predict.h
	//	16777216,  8388608, 5592405, 4194304, 3355443, 2796202, 2396745, 2097152,
	//	 1864135,  1677721, 1525201, 1398101, 1290555, 1198372, 1118481, 1048576,
	//	  986895,   932067,  883011,  838860,  798915,  762600,  729444,  699050,
	//	  671088,   645277,  621378,  599186,  578524,  559240,  541200,  524288,
	//	  508400,   493447,  479349,  466033,  453438,  441505,  430185,  419430,
	//	  409200,   399457,  390167,  381300,  372827,  364722,  356962,  349525,
	//	  342392,   335544,  328965,  322638,  316551,  310689,  305040,  299593,
	//	  294337,   289262,  284359,  279620,  275036,  270600,  266305,  262144,
	//};
	long long coeff[WG_NPREDS], wsum=0;
	int ipred=0;
	int sh=0;
	int pred=0;
//again:
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		unsigned e=
			+(errors[kp]>>9)
			+cN	[kp-1*4*WG_NPREDS]	//+eNW+eWW
			+cN	[kp+0*4*WG_NPREDS]	//+eN+eW
			+cN	[kp+1*4*WG_NPREDS]	//+eNE
			+ccurr	[kp-1*4*WG_NPREDS]	//+eW
			+ccurr	[kp+0*4*WG_NPREDS]	//+eNN+eNW
			+ccurr	[kp+1*4*WG_NPREDS]	//+eNNE+eN
		;
		sh=FLOOR_LOG2(e+1);
		coeff[kp]=((long long)wg_weights[kp]*divlookup[e<<(DENBITS-1)>>sh]<<(DENBITS-1)>>sh)+(1<<DENBITS>>2)/WG_NPREDS;
		//sh=FLOOR_LOG2(e+1)-5;
		//if(sh<0)sh=0;
		//coeff[kp]=((long long)wg_weights[kp]*divlookup[e>>sh]>>sh)+4;
		wsum+=coeff[kp];
	}
	//if(wsum<16)
	//	LOG_ERROR("wsum %d < 16", wsum);
	sh=FLOOR_LOG2(wsum)-(DENBITS-2);
	wsum=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int c=(int)(coeff[kp]>>sh);
		//if((unsigned)c>64)
		//	goto again;
		wsum+=c;
		coeff[kp]=c;
	}
	ipred=(int)((wsum>>1)-1LL);
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		ipred+=(int)(coeff[kp]*preds[kp]);
#ifdef PRINT_MAXCOEFF
		if(maxcoeff<coeff[kp])
			maxcoeff=coeff[kp];
#endif
	}
	//if(wsum<0)
	//	goto again;
	pred=(int)(ipred*divlookup[wsum-1]>>NUMBITS);
#endif
#if 0
	long long ipred=0, wsum=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		unsigned e=
			+(errors[kp]>>9)
			+cN	[kp-1*4*WG_NPREDS]	//+eNW+eWW
			+cN	[kp+0*4*WG_NPREDS]	//+eN+eW
			+cN	[kp+1*4*WG_NPREDS]	//+eNE
			+ccurr	[kp-1*4*WG_NPREDS]	//+eW
			+ccurr	[kp+0*4*WG_NPREDS]	//+eNN+eNW
			+ccurr	[kp+1*4*WG_NPREDS]	//+eNNE+eN
		;
		unsigned weight=(wg_weights[kp]<<10)/(e+1LL);
		ipred+=(long long)weight*preds[kp];
		wsum+=weight;
	}
	ipred+=wsum>>1;
	ipred/=wsum;
	int pred=(int)ipred;
#endif
#if 0
#define PREC 8
	preds[ 9]=preds[ 1]+((preds[ 2]-preds[ 1])*errors[0]>>PREC);
	preds[ 7]=preds[ 3]+((preds[ 4]-preds[ 3])*errors[1]>>PREC);
	preds[ 8]=preds[ 5]+((preds[ 6]-preds[ 5])*errors[2]>>PREC);
	preds[10]=preds[ 7]+((preds[ 8]-preds[ 7])*errors[3]>>PREC);
	preds[11]=preds[ 9]+((preds[10]-preds[ 9])*errors[4]>>PREC);
	preds[12]=preds[11]+((preds[ 0]-preds[11])*errors[5]>>PREC);
	int pred=preds[12];
	//int pred=preds[12]+(errors[6]>>12);
#endif
#if 0
#define PREC 8
	preds[ 6]=preds[ 0]+((preds[ 1]-preds[ 0])*errors[0]>>PREC);
	preds[ 7]=preds[ 2]+((preds[ 3]-preds[ 2])*errors[1]>>PREC);
	preds[ 8]=preds[ 4]+((preds[ 5]-preds[ 4])*errors[2]>>PREC);
	preds[ 9]=preds[ 7]+((preds[ 8]-preds[ 7])*errors[3]>>PREC);
	preds[10]=preds[ 9]+((preds[ 6]-preds[ 9])*errors[4]>>PREC);
	int pred=preds[10];
#endif
#if 0
#define PREC 8
	preds[15]=preds[ 1]+((preds[ 2]-preds[ 1])*errors[0]>>PREC);
	preds[11]=preds[ 3]+((preds[ 4]-preds[ 3])*errors[1]>>PREC);
	preds[12]=preds[ 5]+((preds[ 6]-preds[ 5])*errors[2]>>PREC);
	preds[13]=preds[ 7]+((preds[ 8]-preds[ 7])*errors[3]>>PREC);
	preds[14]=preds[ 9]+((preds[10]-preds[ 9])*errors[4]>>PREC);
	preds[16]=preds[11]+((preds[12]-preds[11])*errors[5]>>PREC);
	preds[17]=preds[13]+((preds[14]-preds[13])*errors[6]>>PREC);
	preds[18]=preds[16]+((preds[17]-preds[16])*errors[7]>>PREC);
	preds[19]=preds[18]+((preds[15]-preds[18])*errors[8]>>PREC);
	preds[20]=preds[19]+((preds[ 0]-preds[19])*errors[9]>>PREC);
	int pred=preds[20];
#endif
#if 0
#define PREC 8
	preds[0x8]=preds[0x0]+((preds[0x1]-preds[0x0])*errors[0]>>PREC);
	preds[0x9]=preds[0x2]+((preds[0x3]-preds[0x2])*errors[1]>>PREC);
	preds[0xA]=preds[0x4]+((preds[0x5]-preds[0x4])*errors[2]>>PREC);
	preds[0xB]=preds[0x6]+((preds[0x7]-preds[0x6])*errors[3]>>PREC);
	preds[0xC]=preds[0x8]+((preds[0x9]-preds[0x8])*errors[4]>>PREC);
	preds[0xD]=preds[0xA]+((preds[0xB]-preds[0xA])*errors[5]>>PREC);
	preds[0xE]=preds[0xC]+((preds[0xD]-preds[0xC])*errors[6]>>PREC);
	int pred=preds[0xE];
#endif
#if 0
#define PREC 8
	preds[0x3]=preds[0x0]+((preds[0x1]-preds[0x0])*errors[0]>>PREC);
	int pred=preds[0x3];
#endif
#if 0
#define PREC 8
	preds[0x4]=preds[0x0]+((preds[0x1]-preds[0x0])*errors[0]>>PREC);
	preds[0x5]=preds[0x2]+((preds[0x3]-preds[0x2])*errors[1]>>PREC);
	preds[0x6]=preds[0x4]+((preds[0x5]-preds[0x4])*errors[2]>>PREC);
	int pred=preds[0x6];
#endif
#if 0
#define PREC 8
	int mix[7]={0};
	//			1	1
	//		2	2	1
	//	1	2	?+I
	//mix[0]=errors[0]+cN[0-1*4*WG_NPREDS]+cN[0+0*4*WG_NPREDS]+cN[0+1*4*WG_NPREDS]+ccurr[0-1*4*WG_NPREDS]+ccurr[0+0*4*WG_NPREDS]+ccurr[0+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[1]=errors[1]+cN[1-1*4*WG_NPREDS]+cN[1+0*4*WG_NPREDS]+cN[1+1*4*WG_NPREDS]+ccurr[1-1*4*WG_NPREDS]+ccurr[1+0*4*WG_NPREDS]+ccurr[1+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[2]=errors[2]+cN[2-1*4*WG_NPREDS]+cN[2+0*4*WG_NPREDS]+cN[2+1*4*WG_NPREDS]+ccurr[2-1*4*WG_NPREDS]+ccurr[2+0*4*WG_NPREDS]+ccurr[2+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[3]=errors[3]+cN[3-1*4*WG_NPREDS]+cN[3+0*4*WG_NPREDS]+cN[3+1*4*WG_NPREDS]+ccurr[3-1*4*WG_NPREDS]+ccurr[3+0*4*WG_NPREDS]+ccurr[3+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[4]=errors[4]+cN[4-1*4*WG_NPREDS]+cN[4+0*4*WG_NPREDS]+cN[4+1*4*WG_NPREDS]+ccurr[4-1*4*WG_NPREDS]+ccurr[4+0*4*WG_NPREDS]+ccurr[4+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[5]=errors[5]+cN[5-1*4*WG_NPREDS]+cN[5+0*4*WG_NPREDS]+cN[5+1*4*WG_NPREDS]+ccurr[5-1*4*WG_NPREDS]+ccurr[5+0*4*WG_NPREDS]+ccurr[5+1*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[6]=errors[6]+cN[6-1*4*WG_NPREDS]+cN[6+0*4*WG_NPREDS]+cN[6+1*4*WG_NPREDS]+ccurr[6-1*4*WG_NPREDS]+ccurr[6+0*4*WG_NPREDS]+ccurr[6+1*4*WG_NPREDS]+(1<<PREC>>1);
	mix[0]=errors[0];
	mix[1]=errors[1];
	mix[2]=errors[2];
	mix[3]=errors[3];
	mix[4]=errors[4];
	mix[5]=errors[5];
	mix[6]=errors[6];
	//mix[0]=cN[0+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[1]=cN[1+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[2]=cN[2+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[3]=cN[3+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[4]=cN[4+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[5]=cN[5+0*4*WG_NPREDS]+(1<<PREC>>1);
	//mix[6]=cN[6+0*4*WG_NPREDS]+(1<<PREC>>1);
	//CLAMP2(mix[0], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[1], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[2], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[3], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[4], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[5], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[6], -(1<<PREC>>3), (1<<PREC)+(1<<PREC>>3));
	//CLAMP2(mix[0], 0, (1<<PREC));
	//CLAMP2(mix[1], 0, (1<<PREC));
	//CLAMP2(mix[2], 0, (1<<PREC));
	//CLAMP2(mix[3], 0, (1<<PREC));
	//CLAMP2(mix[4], 0, (1<<PREC));
	//CLAMP2(mix[5], 0, (1<<PREC));
	//CLAMP2(mix[6], 0, (1<<PREC));
	preds[0x8]=preds[0x0]+((preds[0x1]-preds[0x0])*mix[0]>>PREC);
	preds[0x9]=preds[0x2]+((preds[0x3]-preds[0x2])*mix[1]>>PREC);
	preds[0xA]=preds[0x4]+((preds[0x5]-preds[0x4])*mix[2]>>PREC);
	preds[0xB]=preds[0x6]+((preds[0x7]-preds[0x6])*mix[3]>>PREC);
	preds[0xC]=preds[0x8]+((preds[0x9]-preds[0x8])*mix[4]>>PREC);
	preds[0xD]=preds[0xA]+((preds[0xB]-preds[0xA])*mix[5]>>PREC);
	preds[0xE]=preds[0xC]+((preds[0xD]-preds[0xC])*mix[6]>>PREC);
	//preds[0x8]=preds[0x0]+((preds[0x1]-preds[0x0])*mix[0]>>PREC)+(sse[0]>>9);
	//preds[0x9]=preds[0x2]+((preds[0x3]-preds[0x2])*mix[1]>>PREC)+(sse[1]>>9);
	//preds[0xA]=preds[0x4]+((preds[0x5]-preds[0x4])*mix[2]>>PREC)+(sse[2]>>9);
	//preds[0xB]=preds[0x6]+((preds[0x7]-preds[0x6])*mix[3]>>PREC)+(sse[3]>>9);
	//preds[0xC]=preds[0x8]+((preds[0x9]-preds[0x8])*mix[4]>>PREC)+(sse[4]>>9);
	//preds[0xD]=preds[0xA]+((preds[0xB]-preds[0xA])*mix[5]>>PREC)+(sse[5]>>9);
	//preds[0xE]=preds[0xC]+((preds[0xD]-preds[0xC])*mix[6]>>PREC)+(sse[6]>>9);
	int pred=preds[0xE];
#endif
#if 0
	preds[0x8]=preds[0x0]+((preds[0x1]-preds[0x0])*mixer[0]>>16);
	preds[0x9]=preds[0x2]+((preds[0x3]-preds[0x2])*mixer[1]>>16);
	preds[0xA]=preds[0x4]+((preds[0x5]-preds[0x4])*mixer[2]>>16);
	preds[0xB]=preds[0x6]+((preds[0x7]-preds[0x6])*mixer[3]>>16);
	preds[0xC]=preds[0x8]+((preds[0x9]-preds[0x8])*mixer[4]>>16);
	preds[0xD]=preds[0xA]+((preds[0xB]-preds[0xA])*mixer[5]>>16);
	preds[0xE]=preds[0xC]+((preds[0xD]-preds[0xC])*mixer[6]>>16);
	int pred=preds[0xE];
#endif
#if 0
	float fpred=0;
//	float wsum=0;
	unsigned e2[WG_NPREDS], best=~0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		unsigned e=
			+(errors[kp]>>9)
			+cN[kp-1*4*WG_NPREDS]		//+eNW+eWW
			+cN[kp+0*4*WG_NPREDS]		//+eN+eW
			+cN[kp+1*4*WG_NPREDS]		//+eNE
			+ccurr[kp-1*4*WG_NPREDS]	//+eW
			+ccurr[kp+0*4*WG_NPREDS]	//+eNN+eNW
			+ccurr[kp+1*4*WG_NPREDS]	//+eNNE+eN
		;
		if(best>e)
			best=e;
		e2[kp]=e;
	}
	float softmax[WG_NPREDS], norm=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		unsigned e=(e2[kp]-(best>>1))*2+1;
		float x=(float)wg_weights[kp]/e;
	//	float x=(float)cbrt(exp(-0.003*e)/((double)e*e));
	//	x=roundf(x*0x400);
	//	if(x>0xFFFF)x=0xFFFF;
	//	if(x>0xFFFF)
	//		printf("");
		norm+=x;
		softmax[kp]=x;
	}
	norm=1/norm;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
	//	float coeff=1/(float)e;
	//	float coeff=RCPSS((float)e);
		float coeff=softmax[kp]*norm;
	//	coeff*=wg_weights[kp]*100;
	//	coeff=roundf(coeff);
		fpred+=coeff*preds[kp];
	//	wsum+=coeff;
	}
//	fpred/=wsum;
	int pred=CVTFP32_I32(fpred);
#endif
#if 0
	long long pred=0, wsum=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		//			1	1
		//		2	2	1
		//	1	2	?
		unsigned e=(
			+(errors[kp]>>9)
			+cN[kp-1*4*WG_NPREDS]		//+eNW+eWW
			+cN[kp+0*4*WG_NPREDS]		//+eN+eW
			+cN[kp+1*4*WG_NPREDS]		//+eNE
			+ccurr[kp-1*4*WG_NPREDS]	//+eW
			+ccurr[kp+0*4*WG_NPREDS]	//+eNN+eNW
			+ccurr[kp+1*4*WG_NPREDS]	//+eNNE+eN
		)*8+1;
#ifdef PRINT_MAXERROR
		if(maxerror<e)
			maxerror=e;
#endif
		unsigned coeff=(wg_weights[kp]<<12)/e;
		pred+=(long long)coeff*preds[kp];
		wsum+=coeff;
	}
	pred/=wsum;
#endif
	int vmax=N, vmin=W;
	if(N<W)vmin=N, vmax=W;
	if(vmin>NE)vmin=NE;
	if(vmax<NE)vmax=NE;
	CLAMP2(pred, vmin, vmax);
	return (int)pred;
}
static void wg_update(
	int curr, int pred,
	int *subpreds, int *mixer, int *sse,
	const int *preds,
	int *coeffs, unsigned short *cN, unsigned short *ccurr,
	int *errors, int depth
)
{
#if 0
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int e=abs(curr-preds[kp]);
		int e0=errors[kp];
		errors[kp]+=((e<<12)-errors[kp]+(1<<3>>1))>>3;
		//if(errors[kp]<0)
		//	LOG_ERROR("");
		ccurr[kp]=e;
		cN[kp+1*4*WG_NPREDS]+=e;//eNE+=e
	}
#endif
	int e2[WG_NPREDS], best=0x7FFFFFFF;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int e=abs(curr-preds[kp]);
		if(best>e)
			best=e;
		e2[kp]=e;
	}
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int e=e2[kp]-best;
		errors[kp]+=((e<<6)-errors[kp]+(1<<3>>1))>>3;
		ccurr[kp]=e;
	//	ccurr[kp]=ccurr[kp-1*4*WG_NPREDS]+((e-ccurr[kp-1*4*WG_NPREDS]+(1<<2>>1))>>2);//X
	//	ccurr[kp]=(2*ccurr[kp-1*4*WG_NPREDS]+e+cN[kp+3*4*WG_NPREDS])>>2;//X
	//	cN[kp+1*4*WG_NPREDS]+=e;//eNE+=e
	}
#if 0
	int update[]=
	{
		abs(curr-preds[ 1])-abs(curr-preds[ 2]),
		abs(curr-preds[ 3])-abs(curr-preds[ 4]),
		abs(curr-preds[ 5])-abs(curr-preds[ 6]),
		abs(curr-preds[ 7])-abs(curr-preds[ 8]),
		abs(curr-preds[ 9])-abs(curr-preds[10]),
		abs(curr-preds[11])-abs(curr-preds[ 0]),
	};
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<4>>1))>>4;
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<4>>1))>>4;
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[4]+(1<<4>>1))>>4;
	errors[5]+=(((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[5]+(1<<4>>1))>>4;
	//errors[6]+=(((curr-preds[12])<<12)-errors[6]+(1<<9>>1))>>9;
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[ 1])-abs(curr-preds[ 2]),
		abs(curr-preds[ 3])-abs(curr-preds[ 4]),
		abs(curr-preds[ 5])-abs(curr-preds[ 6]),
		abs(curr-preds[ 7])-abs(curr-preds[ 8]),
		abs(curr-preds[ 9])-abs(curr-preds[10]),
		abs(curr-preds[11])-abs(curr-preds[12]),
		abs(curr-preds[13])-abs(curr-preds[14]),
		abs(curr-preds[16])-abs(curr-preds[17]),
		abs(curr-preds[18])-abs(curr-preds[15]),
		abs(curr-preds[19])-abs(curr-preds[ 0]),
	};
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<4>>1))>>4;
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<4>>1))>>4;
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[4]+(1<<4>>1))>>4;
	errors[5]+=(((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[5]+(1<<4>>1))>>4;
	errors[6]+=(((update[6]>0)-(update[6]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[6]+(1<<4>>1))>>4;
	errors[7]+=(((update[7]>0)-(update[7]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[7]+(1<<4>>1))>>4;
	errors[8]+=(((update[8]>0)-(update[8]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[8]+(1<<4>>1))>>4;
	errors[9]+=(((update[9]>0)-(update[9]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[9]+(1<<4>>1))>>4;
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[ 1])-abs(curr-preds[ 0]),
		abs(curr-preds[ 3])-abs(curr-preds[ 2]),
		abs(curr-preds[ 5])-abs(curr-preds[ 4]),
		abs(curr-preds[ 8])-abs(curr-preds[ 7]),
		abs(curr-preds[ 6])-abs(curr-preds[ 9]),
	};
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<4>>1))>>4;
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<4>>1))>>4;
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[4]+(1<<4>>1))>>4;
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[0x0])-abs(curr-preds[0x1]),
		abs(curr-preds[0x2])-abs(curr-preds[0x3]),
		abs(curr-preds[0x4])-abs(curr-preds[0x5]),
		abs(curr-preds[0x6])-abs(curr-preds[0x7]),
		abs(curr-preds[0x8])-abs(curr-preds[0x9]),
		abs(curr-preds[0xA])-abs(curr-preds[0xB]),
		abs(curr-preds[0xC])-abs(curr-preds[0xD]),
	};
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<3>>1))>>3;
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<3>>1))>>3;
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[4]+(1<<5>>1))>>5;
	errors[5]+=(((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[5]+(1<<4>>1))>>4;
	errors[6]+=(((update[6]>0)-(update[6]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[6]+(1<<5>>1))>>5;
#endif
#if 0
	int update=abs(curr-preds[0x0])-abs(curr-preds[0x1]);
	errors[0]+=(((update>0)-(update<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<3>>1))>>3;
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[0x0])-abs(curr-preds[0x1]),
		abs(curr-preds[0x2])-abs(curr-preds[0x3]),
		abs(curr-preds[0x4])-abs(curr-preds[0x5]),
	};
	int sh[]=
	{
		FLOOR_LOG2(abs(curr-preds[0x4]))+3,
		FLOOR_LOG2(abs(curr-preds[0x5]))+3,
		FLOOR_LOG2(abs(curr-preds[0x6]))+3,
	};
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<sh[0]>>1))>>sh[0];
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<sh[1]>>1))>>sh[1];
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<sh[2]>>1))>>sh[2];
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[0x0])-abs(curr-preds[0x1]),
		abs(curr-preds[0x2])-abs(curr-preds[0x3]),
		abs(curr-preds[0x4])-abs(curr-preds[0x5]),
		abs(curr-preds[0x6])-abs(curr-preds[0x7]),
		abs(curr-preds[0x8])-abs(curr-preds[0x9]),
		abs(curr-preds[0xA])-abs(curr-preds[0xB]),
		abs(curr-preds[0xC])-abs(curr-preds[0xD]),
	};
	//errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[0]+(1<<3>>1))>>3;
	//errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[1]+(1<<3>>1))>>3;
	//errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[2]+(1<<4>>1))>>4;
	//errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[3]+(1<<4>>1))>>4;
	//errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[4]+(1<<5>>1))>>5;
	//errors[5]+=(((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[5]+(1<<4>>1))>>4;
	//errors[6]+=(((update[6]>0)-(update[6]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1)-errors[6]+(1<<5>>1))>>5;
	int sh[]=
	{
		abs(curr-preds[0x8]),
		abs(curr-preds[0x9]),
		abs(curr-preds[0xA]),
		abs(curr-preds[0xB]),
		abs(curr-preds[0xC]),
		abs(curr-preds[0xD]),
		abs(curr-preds[0xE]),
	};
	sh[0]=(FLOOR_LOG2(sh[0]+1)>>1)+3;
	sh[1]=(FLOOR_LOG2(sh[1]+1)>>1)+3;
	sh[2]=(FLOOR_LOG2(sh[2]+1)>>1)+3;
	sh[3]=(FLOOR_LOG2(sh[3]+1)>>1)+3;
	sh[4]=(FLOOR_LOG2(sh[4]+1)>>1)+3;
	sh[5]=(FLOOR_LOG2(sh[5]+1)>>1)+3;
	sh[6]=(FLOOR_LOG2(sh[6]+1)>>1)+3;
	errors[0]+=(((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[0]+(1<<sh[0]>>1))>>sh[0];
	errors[1]+=(((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[1]+(1<<sh[1]>>1))>>sh[1];
	errors[2]+=(((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[2]+(1<<sh[2]>>1))>>sh[2];
	errors[3]+=(((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[3]+(1<<sh[3]>>1))>>sh[3];
	errors[4]+=(((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[4]+(1<<sh[4]>>1))>>sh[4];
	errors[5]+=(((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[5]+(1<<sh[5]>>1))>>sh[5];
	errors[6]+=(((update[6]>0)-(update[6]<0))*((1<<PREC>>1)+(1<<PREC>>5))+(1<<PREC>>1)-errors[6]+(1<<sh[6]>>1))>>sh[6];
	//errors[0]+=((((update[0]>0)-(update[0]<0)+1)<<PREC>>1)-errors[0]+(1<<sh[0]>>1))>>sh[0];
	//errors[1]+=((((update[1]>0)-(update[1]<0)+1)<<PREC>>1)-errors[1]+(1<<sh[1]>>1))>>sh[1];
	//errors[2]+=((((update[2]>0)-(update[2]<0)+1)<<PREC>>1)-errors[2]+(1<<sh[2]>>1))>>sh[2];
	//errors[3]+=((((update[3]>0)-(update[3]<0)+1)<<PREC>>1)-errors[3]+(1<<sh[3]>>1))>>sh[3];
	//errors[4]+=((((update[4]>0)-(update[4]<0)+1)<<PREC>>1)-errors[4]+(1<<sh[4]>>1))>>sh[4];
	//errors[5]+=((((update[5]>0)-(update[5]<0)+1)<<PREC>>1)-errors[5]+(1<<sh[5]>>1))>>sh[5];
	//errors[6]+=((((update[6]>0)-(update[6]<0)+1)<<PREC>>1)-errors[6]+(1<<sh[6]>>1))>>sh[6];
#endif
#if 0
	int update[]=
	{
		((abs(curr-preds[0x0])-abs(curr-preds[0x1]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0x2])-abs(curr-preds[0x3]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0x4])-abs(curr-preds[0x5]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0x6])-abs(curr-preds[0x7]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0x8])-abs(curr-preds[0x9]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0xA])-abs(curr-preds[0xB]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
		((abs(curr-preds[0xC])-abs(curr-preds[0xD]))*((1<<PREC>>1)+(1<<PREC>>4))*2>>depth)+(1<<PREC>>1),
	};
	if((preds[0x0]<curr&&curr<preds[0x1])||(preds[0x0]>curr&&curr>preds[0x1]))update[0]=preds[0x1]==preds[0x0]?0:((curr-preds[0x0])<<PREC)/(preds[0x1]-preds[0x0]);
	if((preds[0x2]<curr&&curr<preds[0x3])||(preds[0x2]>curr&&curr>preds[0x3]))update[1]=preds[0x3]==preds[0x2]?0:((curr-preds[0x2])<<PREC)/(preds[0x3]-preds[0x2]);
	if((preds[0x4]<curr&&curr<preds[0x5])||(preds[0x4]>curr&&curr>preds[0x5]))update[2]=preds[0x5]==preds[0x4]?0:((curr-preds[0x4])<<PREC)/(preds[0x5]-preds[0x4]);
	if((preds[0x6]<curr&&curr<preds[0x7])||(preds[0x6]>curr&&curr>preds[0x7]))update[3]=preds[0x7]==preds[0x6]?0:((curr-preds[0x6])<<PREC)/(preds[0x7]-preds[0x6]);
	if((preds[0x8]<curr&&curr<preds[0x9])||(preds[0x8]>curr&&curr>preds[0x9]))update[4]=preds[0x9]==preds[0x8]?0:((curr-preds[0x8])<<PREC)/(preds[0x9]-preds[0x8]);
	if((preds[0xA]<curr&&curr<preds[0xB])||(preds[0xA]>curr&&curr>preds[0xB]))update[5]=preds[0xB]==preds[0xA]?0:((curr-preds[0xA])<<PREC)/(preds[0xB]-preds[0xA]);
	if((preds[0xC]<curr&&curr<preds[0xD])||(preds[0xC]>curr&&curr>preds[0xD]))update[6]=preds[0xD]==preds[0xC]?0:((curr-preds[0xC])<<PREC)/(preds[0xD]-preds[0xC]);
	//if((preds[0x0]<curr)==(curr<preds[0x1]))
	//volatile unsigned int flags=0;
	//__asm volatile("pushf; pop %0" : "=r" (flags));
	//int update[]=
	//{
	//	abs(curr-preds[0x0])-abs(curr-preds[0x1]),
	//	abs(curr-preds[0x2])-abs(curr-preds[0x3]),
	//	abs(curr-preds[0x4])-abs(curr-preds[0x5]),
	//	abs(curr-preds[0x6])-abs(curr-preds[0x7]),
	//	abs(curr-preds[0x8])-abs(curr-preds[0x9]),
	//	abs(curr-preds[0xA])-abs(curr-preds[0xB]),
	//	abs(curr-preds[0xC])-abs(curr-preds[0xD]),
	//};
	//update[0]=((update[0]>0)-(update[0]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[1]=((update[1]>0)-(update[1]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[2]=((update[2]>0)-(update[2]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[3]=((update[3]>0)-(update[3]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[4]=((update[4]>0)-(update[4]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[5]=((update[5]>0)-(update[5]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	//update[6]=((update[6]>0)-(update[6]<0))*((1<<PREC>>1)+(1<<PREC>>4))+(1<<PREC>>1);
	errors[0]+=(update[0]-errors[0]+(1<<3>>1))>>3;
	errors[1]+=(update[1]-errors[1]+(1<<3>>1))>>3;
	errors[2]+=(update[2]-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(update[3]-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(update[4]-errors[4]+(1<<5>>1))>>5;
	errors[5]+=(update[5]-errors[5]+(1<<4>>1))>>4;
	errors[6]+=(update[6]-errors[6]+(1<<5>>1))>>5;
#endif
#if 0
	int update[]=
	{
		preds[0x1]==preds[0x0]?0:((curr-preds[0x0])<<PREC)/(preds[0x1]-preds[0x0]),//table[0]==0
		preds[0x3]==preds[0x2]?0:((curr-preds[0x2])<<PREC)/(preds[0x3]-preds[0x2]),
		preds[0x5]==preds[0x4]?0:((curr-preds[0x4])<<PREC)/(preds[0x5]-preds[0x4]),
		preds[0x7]==preds[0x6]?0:((curr-preds[0x6])<<PREC)/(preds[0x7]-preds[0x6]),
		preds[0x9]==preds[0x8]?0:((curr-preds[0x8])<<PREC)/(preds[0x9]-preds[0x8]),
		preds[0xB]==preds[0xA]?0:((curr-preds[0xA])<<PREC)/(preds[0xB]-preds[0xA]),
		preds[0xD]==preds[0xC]?0:((curr-preds[0xC])<<PREC)/(preds[0xD]-preds[0xC]),
	};
	CLAMP2(update[0], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[1], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[2], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[3], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[4], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[5], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	CLAMP2(update[6], -(1<<PREC>>1), (1<<PREC)+(1<<PREC>>1));
	errors[0]+=(update[0]-errors[0]+(1<<7>>1))>>7;
	errors[1]+=(update[1]-errors[1]+(1<<7>>1))>>7;
	errors[2]+=(update[2]-errors[2]+(1<<7>>1))>>7;
	errors[3]+=(update[3]-errors[3]+(1<<7>>1))>>7;
	errors[4]+=(update[4]-errors[4]+(1<<7>>1))>>7;
	errors[5]+=(update[5]-errors[5]+(1<<7>>1))>>7;
	errors[6]+=(update[6]-errors[6]+(1<<7>>1))>>7;
#endif
#if 0
	int update[]=
	{
		abs(curr-preds[0x0])-abs(curr-preds[0x1]),
		abs(curr-preds[0x2])-abs(curr-preds[0x3]),
		abs(curr-preds[0x4])-abs(curr-preds[0x5]),
		abs(curr-preds[0x6])-abs(curr-preds[0x7]),
		abs(curr-preds[0x8])-abs(curr-preds[0x9]),
		abs(curr-preds[0xA])-abs(curr-preds[0xB]),
		abs(curr-preds[0xC])-abs(curr-preds[0xD]),
		//(preds[0x0]<curr)==(preds[0x1]<curr)?abs(curr-preds[0x0])<abs(curr-preds[0x1])?0:1<<PREC:1<<PREC>>1,
		//(preds[0x2]<curr)==(preds[0x3]<curr)?abs(curr-preds[0x2])<abs(curr-preds[0x3])?0:1<<PREC:1<<PREC>>1,
		//(preds[0x4]<curr)==(preds[0x5]<curr)?abs(curr-preds[0x4])<abs(curr-preds[0x5])?0:1<<PREC:1<<PREC>>1,
		//(preds[0x6]<curr)==(preds[0x7]<curr)?abs(curr-preds[0x6])<abs(curr-preds[0x7])?0:1<<PREC:1<<PREC>>1,
		//(preds[0x8]<curr)==(preds[0x9]<curr)?abs(curr-preds[0x8])<abs(curr-preds[0x9])?0:1<<PREC:1<<PREC>>1,
		//(preds[0xA]<curr)==(preds[0xB]<curr)?abs(curr-preds[0xA])<abs(curr-preds[0xB])?0:1<<PREC:1<<PREC>>1,
		//(preds[0xC]<curr)==(preds[0xD]<curr)?abs(curr-preds[0xC])<abs(curr-preds[0xD])?0:1<<PREC:1<<PREC>>1,
	};
	//errors[0]+=update[0];
	//errors[1]+=update[1];
	//errors[2]+=update[2];
	//errors[3]+=update[3];
	//errors[4]+=update[4];
	//errors[5]+=update[5];
	//errors[6]+=update[6];
	//CLAMP2(errors[0], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[1], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[2], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[3], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[4], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[5], -(1<<PREC>>1), (1<<PREC>>1));
	//CLAMP2(errors[6], -(1<<PREC>>1), (1<<PREC>>1));
	//ccurr[0+0*4*WG_NPREDS]=update[0];
	//ccurr[1+0*4*WG_NPREDS]=update[1];
	//ccurr[2+0*4*WG_NPREDS]=update[2];
	//ccurr[3+0*4*WG_NPREDS]=update[3];
	//ccurr[4+0*4*WG_NPREDS]=update[4];
	//ccurr[5+0*4*WG_NPREDS]=update[5];
	//ccurr[6+0*4*WG_NPREDS]=update[6];
	//cN[0+1*4*WG_NPREDS]+=update[0];
	//cN[1+1*4*WG_NPREDS]+=update[1];
	//cN[2+1*4*WG_NPREDS]+=update[2];
	//cN[3+1*4*WG_NPREDS]+=update[3];
	//cN[4+1*4*WG_NPREDS]+=update[4];
	//cN[5+1*4*WG_NPREDS]+=update[5];
	//cN[6+1*4*WG_NPREDS]+=update[6];
	//update[0]=((update[0]>>31)-((-update[0])>>31)+1)<<PREC>>1;
	//update[1]=((update[1]>>31)-((-update[1])>>31)+1)<<PREC>>1;
	//update[2]=((update[2]>>31)-((-update[2])>>31)+1)<<PREC>>1;
	//update[3]=((update[3]>>31)-((-update[3])>>31)+1)<<PREC>>1;
	//update[4]=((update[4]>>31)-((-update[4])>>31)+1)<<PREC>>1;
	//update[5]=((update[5]>>31)-((-update[5])>>31)+1)<<PREC>>1;
	//update[6]=((update[6]>>31)-((-update[6])>>31)+1)<<PREC>>1;
	update[0]=((update[0]>0)-(update[0]<0)+1)<<PREC>>1;
	update[1]=((update[1]>0)-(update[1]<0)+1)<<PREC>>1;
	update[2]=((update[2]>0)-(update[2]<0)+1)<<PREC>>1;
	update[3]=((update[3]>0)-(update[3]<0)+1)<<PREC>>1;
	update[4]=((update[4]>0)-(update[4]<0)+1)<<PREC>>1;
	update[5]=((update[5]>0)-(update[5]<0)+1)<<PREC>>1;
	update[6]=((update[6]>0)-(update[6]<0)+1)<<PREC>>1;
	//update[0]+=1<<PREC>>1;
	//update[1]+=1<<PREC>>1;
	//update[2]+=1<<PREC>>1;
	//update[3]+=1<<PREC>>1;
	//update[4]+=1<<PREC>>1;
	//update[5]+=1<<PREC>>1;
	//update[6]+=1<<PREC>>1;
	errors[0]+=(update[0]-errors[0]+(1<<4>>1))>>4;
	errors[1]+=(update[1]-errors[1]+(1<<4>>1))>>4;
	errors[2]+=(update[2]-errors[2]+(1<<4>>1))>>4;
	errors[3]+=(update[3]-errors[3]+(1<<4>>1))>>4;
	errors[4]+=(update[4]-errors[4]+(1<<4>>1))>>4;
	errors[5]+=(update[5]-errors[5]+(1<<4>>1))>>4;
	errors[6]+=(update[6]-errors[6]+(1<<4>>1))>>4;
	//sse[0]+=(((curr-preds[0x8])<<8)-sse[0]+(1<<7>>1))>>7;
	//sse[1]+=(((curr-preds[0x9])<<8)-sse[1]+(1<<7>>1))>>7;
	//sse[2]+=(((curr-preds[0xA])<<8)-sse[2]+(1<<7>>1))>>7;
	//sse[3]+=(((curr-preds[0xB])<<8)-sse[3]+(1<<7>>1))>>7;
	//sse[4]+=(((curr-preds[0xC])<<8)-sse[4]+(1<<7>>1))>>7;
	//sse[5]+=(((curr-preds[0xD])<<8)-sse[5]+(1<<7>>1))>>7;
	//sse[6]+=(((curr-preds[0xE])<<8)-sse[6]+(1<<7>>1))>>7;
#endif
#if 0
	mixer[0]+=(abs(curr-preds[0x0])-abs(curr-preds[0x1]))<<10;
	mixer[1]+=(abs(curr-preds[0x2])-abs(curr-preds[0x3]))<<10;
	mixer[2]+=(abs(curr-preds[0x4])-abs(curr-preds[0x5]))<<10;
	mixer[3]+=(abs(curr-preds[0x6])-abs(curr-preds[0x7]))<<10;
	mixer[4]+=(abs(curr-preds[0x8])-abs(curr-preds[0x9]))<<10;
	mixer[5]+=(abs(curr-preds[0xA])-abs(curr-preds[0xB]))<<10;
	mixer[6]+=(abs(curr-preds[0xC])-abs(curr-preds[0xD]))<<10;
	CLAMP2(mixer[0], 0, 0x10000);
	CLAMP2(mixer[1], 0, 0x10000);
	CLAMP2(mixer[2], 0, 0x10000);
	CLAMP2(mixer[3], 0, 0x10000);
	CLAMP2(mixer[4], 0, 0x10000);
	CLAMP2(mixer[5], 0, 0x10000);
	CLAMP2(mixer[6], 0, 0x10000);
#endif
#if 0
	int abserrors[WG_NPREDS];
	int best=0x7FFFFFFF, worst=0;
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int e=abs(curr-preds[kp]);
		if(best>e)
			best=e;
		if(worst<e)
			worst=e;
		abserrors[kp]=e;
	}
	for(int kp=0;kp<WG_NPREDS;++kp)
	{
		int e=abserrors[kp]-best;
	//	int e=((abserrors[kp]+1-best)<<10)/(worst-best+1);
	//	int e=abs(curr-preds[kp]);
		errors[kp]+=((e<<12)-errors[kp]+(1<<3>>1))>>3;
		ccurr[kp]=e;
		cN[kp+1*4*WG_NPREDS]+=e;//eNE+=e
	}
#endif
}
void pred_wgrad6(Image *src, int fwd)
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
	ALIGN(32) int wg_coeffs[4][WG_NPREDS]={0}, wg_errors[4][WG_NPREDS]={0}, wg_preds[WG_NPREDS*2]={0};
	int subpreds[2]={0}, mixer[4][WG_NPREDS]={0}, sse[4][WG_NPREDS]={0};

	int fwdmask=-fwd;
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	int pesize=(src->iw+16)*(int)sizeof(short[2][4][WG_NPREDS]);//2 padded rows * 4 channels * WG_NPREDS
	unsigned short *ebuf=(unsigned short*)_mm_malloc(pesize, sizeof(__m256i));
	int bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	short *pixels=(short*)malloc(bufsize);
	//int wgerrorsize=sizeof(int[4][WG_NCTX][WG_NPREDS]);
	//int *wg_errors=(int*)malloc(wgerrorsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	//FILLMEM(ebuf, 0x10000/WG_NPREDS, pesize, sizeof(int));
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	FILLMEM((int*)wg_coeffs, 0x10000/WG_NPREDS, sizeof(wg_coeffs), sizeof(int));
#ifdef PREC
	FILLMEM((int*)wg_errors, (1<<PREC>>1), sizeof(wg_errors), sizeof(int));
#endif
	//FILLMEM(wg_errors, 0x10000/WG_NPREDS, wgerrorsize, sizeof(int));
	//memset(wg_errors, 0, wgerrorsize);
	UPDATE_MAX(nch, src->nch);
	//wg_init(wg_weights);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		//if(!(ky%(src->ih>>4)))
		//{
		//	memset(wg_errors, 0, sizeof(wg_errors));
		//	memset(ebuf, 0, pesize);
		//}
	//	int cpred2=0;
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
			ebuf+((src->iw+16LL)*((ky-0LL)&1)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-1LL)&1)+8LL)*4*WG_NPREDS,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				//if(kx==src->iw/2&&ky==src->ih/2)//
				//	printf("");

				//int
				//	kc2=kc<<1,
				//	kc3=kc*WG_NPREDS,
				//	*eNNE	=erows[2]+kc3+1*4*WG_NPREDS,
				//	*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
				//	*eN	=erows[1]+kc3+0*4*WG_NPREDS,
				//	*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
				//	*eW	=erows[0]+kc3-1*4*WG_NPREDS,
				//	*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
			//	int ctx=0;
			//	int ctx=cpred*WG_NCTX>>src->depth[kc]&(WG_NCTX-1);
			//	int ctx=curr*WG_NCTX>>src->depth[kc]&(WG_NCTX-1);
			//	int *currerrors=wg_errors+WG_NPREDS*((ptrdiff_t)WG_NCTX*kc+ctx);
				int pred=wg_predict(
					wg_errors[kc],
					wg_coeffs[kc],
					erows[1]+(ptrdiff_t)WG_NPREDS*kc,
					erows[0]+(ptrdiff_t)WG_NPREDS*kc,
					rows[3]+kc,
					rows[2]+kc,
					rows[1]+kc,
					rows[0]+kc,
					subpreds,
					mixer[kc], sse[kc],
					wg_preds,
					src->depth[kc],
					kc
				);
				curr=src->data[idx+kc];

				int val;
				if(fwd)
				{
					val=(curr-(int)pred+g_dist/2)/g_dist;
					curr=g_dist*val+(int)pred;
				}
				else
				{
					val=g_dist*curr+(int)pred;
					curr=val;
					CLAMP2(val, amin[kc], amax[kc]);
				}
				src->data[idx+kc]=val;
				rows[0][kc+0]=curr;
				rows[0][kc+4]=curr-pred;
#if 0
			//	cpred2=cpred;
				cpred=pred;
			//	cpred+=((1<<custom_params[0])-1)&cpred>>31;//rounding to zero
			//	cpred>>=custom_params[0];
			// 
			//	cpred=(pred+(1<<custom_params[0]>>1))>>custom_params[0];//rounding to nearest
				cpred^=fwdmask;
				cpred-=fwdmask;
				cpred+=curr;

				cpred<<=32-src->depth[kc];
				cpred>>=32-src->depth[kc];

				src->data[idx+kc]=cpred;
				rows[0][kc+0]=(fwd?curr:cpred)<<custom_params[0];
				rows[0][kc+4]=rows[0][kc+0]-pred;
#endif
				//if(kx==1)//
				//	printf("");
				wg_update(
					rows[0][kc+0], pred,
					subpreds, mixer[kc], sse[kc],
					wg_preds,
					wg_coeffs[kc],
					erows[1]+(ptrdiff_t)WG_NPREDS*kc,
					erows[0]+(ptrdiff_t)WG_NPREDS*kc,
					wg_errors[kc],
					src->depth[kc]
				);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
			erows[0]+=4*WG_NPREDS;
			erows[1]+=4*WG_NPREDS;
			//erows[2]+=3*WG_NPREDS;
			//erows[3]+=3*WG_NPREDS;
		}
	}
#ifdef PRINT_MAXCOEFF
	if(loud_transforms)
		messagebox(MBOX_OK, "Info", "max coeff %d", maxcoeff);
#endif
#ifdef PRINT_MAXERROR
	if(loud_transforms)
		messagebox(MBOX_OK, "Info", "max error %d", maxerror);
#endif
	free(pixels);
	//free(wg_errors);
	_mm_free(ebuf);
}