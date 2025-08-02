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


	#define ENABLE_EXTENDED_RCT


#if 0
#define L1SH 19
#define NPREDS 8
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 40000, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)\
	PRED( 50000, W+NE-N)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)
#endif
#if 1
//#define L1SH 17
#define L1SH 19
#define NPREDS 11		//up to 11, otherwise slow
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 40000, NNN)\
	PRED( 40000, WWW)\
	PRED( 40000, NEEE)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 50000, W+NE-N)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)\
	PRED( 40000, (WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2)
//	PRED( 40000, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)
#endif
#if 0
#define L1SH 19
#define NPREDS 15
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 50000, W+NE-N)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)\
	PRED( 40000, N+NW-NNW)\
	PRED( 40000, W+NW-NWW)\
	PRED( 40000, NEEE)\
	PRED( 40000, NW)\
	PRED( 40000, NE)\
	PRED( 40000, NN)\
	PRED( 40000, WW)
#endif
#if 0
#define NPREDS 10
#define PREDLIST\
	PRED( 40000, N)\
	PRED( 40000, W)\
	PRED( 40000, 3*(N-NN)+NNN)\
	PRED( 40000, 3*(W-WW)+WWW)\
	PRED( 40000, W+NE-N)\
	PRED(160000, N+W-NW)\
	PRED( 40000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 40000, N+NE-NNE)\
	PRED( 40000, W+NW-NWW)\
	PRED( 40000, NEEE)
#endif
#if 0
#define NPREDS 10
#define PREDLIST\
	PRED( 38000, N)\
	PRED( 69000, W)\
	PRED( 41000, 3*(N-NN)+NNN)\
	PRED( 72000, 3*(W-WW)+WWW)\
	PRED( 70000, W+NE-N)\
	PRED( 83000, N+W-NW)\
	PRED(-10000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 61000, N+NE-NNE)\
	PRED( 81000, W+NW-NWW)\
	PRED( 18000, NEEE)
#endif

void pred_ols7(Image *src, int fwd)
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
	int weights[4][NPREDS]={0};
	int invdist=((1<<16)+g_dist-1)/g_dist;
	int bufsize=(src->iw+8*2)*(int)sizeof(short[4*4*1]);//4 padded rows * 4 channels max * {pixels}
	short *pixels=(short*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
//	static const int w0[]=
//	{
//#define PRED(W0, EXPR) W0,
//		PREDLIST
//#undef  PRED
//	};
//	for(int kc=0;kc<4;++kc)
//	{
//		for(int kp=0;kp<NPREDS;++kp)
//			weights[kc][kp]=w0[kp];
//	}
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*1,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
				int
					NNNWWW		=rows[3][-3*4*1],
					NNNW		=rows[3][-1*4*1],
					NNN		=rows[3][+0*4*1],
					NNNE		=rows[3][+1*4*1],
					NNNEE		=rows[3][+2*4*1],
					NNNEEE		=rows[3][+3*4*1],
					NNNEEEE		=rows[3][+4*4*1],
					NNWWWW		=rows[2][-4*4*1],
					NNWWW		=rows[2][-3*4*1],
					NNWW		=rows[2][-2*4*1],
					NNW		=rows[2][-1*4*1],
					NN		=rows[2][+0*4*1],
					NNE		=rows[2][+1*4*1],
					NNEE		=rows[2][+2*4*1],
					NNEEE		=rows[2][+3*4*1],
					NNEEEE		=rows[2][+4*4*1],
					NWWWW		=rows[1][-4*4*1],
					NWWW		=rows[1][-3*4*1],
					NWW		=rows[1][-2*4*1],
					NW		=rows[1][-1*4*1],
					N		=rows[1][+0*4*1],
					NE		=rows[1][+1*4*1],
					NEE		=rows[1][+2*4*1],
					NEEE		=rows[1][+3*4*1],
					NEEEE		=rows[1][+4*4*1],
					NEEEEE		=rows[1][+5*4*1],
					NEEEEEE		=rows[1][+6*4*1],
					NEEEEEEE	=rows[1][+7*4*1],
					NEEEEEEEE	=rows[1][+8*4*1],
					WWWWWWWWW	=rows[0][-9*4*1],
					WWWWWWWW	=rows[0][-8*4*1],
					WWWWWWW		=rows[0][-7*4*1],
					WWWWWW		=rows[0][-6*4*1],
					WWWWW		=rows[0][-5*4*1],
					WWWW		=rows[0][-4*4*1],
					WWW		=rows[0][-3*4*1],
					WW		=rows[0][-2*4*1],
					W		=rows[0][-1*4*1];
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int *currw=weights[kc];
				int p0=1LL<<L1SH>>1;

				//int idx2=custom_params[1];
				//if(idx2>=0&&idx2<NPREDS)
				//	currw[idx2]=0;
				
				for(int k=0;k<NPREDS;++k)
					p0+=currw[k]*preds[k];
				p0>>=L1SH;

				//for(int k=0;k<NPREDS;++k)
				//	p0+=(currw[k]>>8)*preds[k];
				//p0+=1LL<<(L1SH-8)>>1;
				//p0>>=L1SH-8;
			//	p0-=p0>>31;//X  deadzone bad with advanced pred
				int predc=p0;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(predc, vmin, vmax);

				int curr=src->data[idx];
				//if(ky==src->ih/2&&kx==src->iw/2)
				//	printf("");
				
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=predc;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx]=curr;

						curr=g_dist*curr+predc;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+predc;
						CLAMP2(curr, amin[kc], amax[kc]);

						src->data[idx]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						int error=curr-predc;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx]=error;
					}
					else
					{
						curr+=predc;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx]=curr;
					}
				}
			//	//near lossless
			//	if(fwd)
			//	{
			//		curr-=predc;
			//		curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
			//		src->data[idx]=curr;
			//
			//		curr=g_dist*curr+predc;
			//		CLAMP2(curr, amin[kc], amax[kc]);
			//	}
			//	else
			//	{
			//		curr=g_dist*curr+predc;
			//		CLAMP2(curr, amin[kc], amax[kc]);
			//
			//		src->data[idx]=curr;
			//	}
				rows[0][0]=curr;

				//update
			//	int e=curr-p0;//L2
				int e=(curr>p0)-(curr<p0);//L1
			//	currw[NPREDS]+=e;
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];
			}
		}
	}
	free(pixels);
}
void pred_mixN(Image *src, int fwd)
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
	int nch;
	int fwdmask=-fwd;
	
	int invdist=((1<<16)+g_dist-1)/g_dist;
	//int coeffsA[4][4]={0};
	//int coeffsB[4][4]={0};
	//int coeffsC[4][4]={0};
	//int coeffsD[4][4]={0};
	int coeffs[4][14]={0};
	//int coeffs[4][9]={0};

	int bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max * {pixel, error}
	int *pixels=(int*)malloc(bufsize);

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	//FILLMEM((int*)coeffs, (1<<L1SH)/NPREDS, sizeof(coeffs), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4*2,
		};
		int asum=0;
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<3;++kc)
			{
				int
					NNN	=rows[3][kc+0*4*2+0],
					NNWW	=rows[2][kc-2*4*2+0],
					NNW	=rows[2][kc-1*4*2+0],
					NN	=rows[2][kc+0*4*2+0],
					NNE	=rows[2][kc+1*4*2+0],
					NNEE	=rows[2][kc+2*4*2+0],
					NWW	=rows[1][kc-2*4*2+0],
					NW	=rows[1][kc-1*4*2+0],
					N	=rows[1][kc+0*4*2+0],
					NE	=rows[1][kc+1*4*2+0],
					NEE	=rows[1][kc+2*4*2+0],
					NEEE	=rows[1][kc+3*4*2+0],
					NEEEE	=rows[1][kc+4*4*2+0],
					WWWW	=rows[0][kc-4*4*2+0],
					WWW	=rows[0][kc-3*4*2+0],
					WW	=rows[0][kc-2*4*2+0],
					W	=rows[0][kc-1*4*2+0],
					eNW	=rows[1][kc-1*4*2+1],
					eN	=rows[1][kc+0*4*2+1],
					eNE	=rows[1][kc+1*4*2+1],
					eWW	=rows[0][kc-2*4*2+1],
					eW	=rows[0][kc-1*4*2+1];
#if 0
				int *currw=coeffs[kc];

				//             NNN
				//             NN
				//         NW  N  NE  NEE
				//WWW  WW  W   ?

				//int paeth=abs(N-NW)>abs(W-NW)?N:W;
				//if(abs(NW-((N+W)>>1))<?)
				//	paeth=NW;
				int px=3*(W-WW)+WWW, py=3*(N-NN)+NNN;
				int g=N+W-NW;
				int pred=
					+currw[0]*N
					+currw[1]*W
					+currw[2]*NE
					+currw[3]*g
					+currw[4]*py
					+currw[5]*px

				//	+currw[6]*paeth//X
				//	+currw[6]*(N+NE-NNE)
				;
				pred+=1<<18>>1;
				pred>>=18;
			//	pred-=pred>>31;
				int p0=pred;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx+kc];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)pred;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx+kc]=curr;

						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);

						src->data[idx+kc]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						int error=curr-pred;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx+kc]=error;
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx+kc]=curr;
					}
				}
				rows[0][kc+0]=curr;
			//	rows[0][kc+1]=curr-pred;

				int t0;
				t0=N;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[0]+=t0;
				t0=W;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[1]+=t0;
				t0=NE;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[2]+=t0;
				t0=g;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[3]+=t0;
				t0=py;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[4]+=t0;
				t0=px;		if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[5]+=t0;

			//	t0=paeth;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[6]+=t0;//X
			//	t0=N+NE-NNE;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; currw[6]+=t0;
#if 0
				int e=(curr>p0)-(curr<p0);//L1
				currw[0]+=e*N;
				currw[1]+=e*W;
				currw[2]+=e*NE;
				currw[3]+=e*(3*(N-NN)+NNN);
				currw[4]+=e*(3*(W-WW)+WWW);
				currw[5]+=e*(N+W-NW);
			//	currw[6]+=e*NW;
#endif
#endif
#if 0
				int *currw=coeffs[kc];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

				//            NNN
				//            NN
				//    NWW NW  N   NE
				//WWW WW  W   ?
				int dW=W-WW, dWW=WW-WWW, dN=N-NN, dNN=NN-NNN, g=N+W-NW, dNE=N-NE, dNWW=NW-NWW;
				int pred=
					+currw[0]*dNN
					+currw[1]*dN
					+currw[2]*g
					+currw[3]*N
					+currw[4]*dNE
					+currw[5]*dWW
					+currw[6]*dW
					+currw[7]*W
					+currw[8]*dNWW
				;
				pred+=1<<17>>1;
				pred>>=17;
				int p0=pred;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx+kc];
				if(fwd)
				{
					curr-=(int)pred;
					curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
					src->data[idx+kc]=curr;

					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);
				}
				else
				{
					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);

					src->data[idx+kc]=curr;
				}
				rows[0][kc+0]=curr;
				rows[0][kc+1]=curr-pred;
				int e=(curr>p0)-(curr<p0);//L1
				currw[0]+=e*dNN;
				currw[1]+=e*dN;
				currw[2]+=e*g;
				currw[3]+=e*N;
				currw[4]+=e*dNE;
				currw[5]+=e*dWW;
				currw[6]+=e*dW;
				currw[7]+=e*W;
				currw[8]+=e*dNWW;
#endif
#if 1
				int *currw=coeffs[kc];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

				//NW  N   NE
				//W   ?
				int pred=
					+currw[0]*W
					+currw[1]*(N-NW)
					+currw[2]*(NE-N)
					+currw[3]*(N-W)
				;
			//	int pred=
			//		+currw[0]*N
			//		+currw[1]*W
			//		+currw[2]*NW
			//		+currw[3]*NE
			//	;
			//	pred+=((1<<17)-1)&pred>>31;//rounding to zero	X
				pred+=1<<17>>1;//rounding to nearest
				pred>>=17;
				int p0=pred;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				//if(vmin>NEEE)vmin=NEEE;
				//if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx+kc];
				if(fwd)
				{
					curr-=(int)pred;
					curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
					src->data[idx+kc]=curr;

					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);
				}
				else
				{
					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);

					src->data[idx+kc]=curr;
				}
				rows[0][kc+0]=curr;
			//	rows[0][kc+1]=curr-pred;
				int e=(curr>p0)-(curr<p0);//L1
			//	int e=curr-p0;//L2
				currw[0]+=e*W;
				currw[1]+=e*(N-NW);
				currw[2]+=e*(NE-N);
				currw[3]+=e*(N-W);
			//	CLAMP2(currw[0], (1<<17)/4, (1<<17));
			//	CLAMP2(currw[1], -(1<<17)/4, (1<<17)*3/4);
			//	currw[2]=(1<<17)/4;
			//	currw[3]=(1<<17)-currw[0];

				//currw[0]+=e*N;
				//currw[1]+=e*W;
				//currw[2]+=e*NW;
				//currw[3]+=e*NE;
#endif
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
#if 0
		for(int it=0;it<1;++it)
		{
			rows[0]=pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4*2;
			rows[1]=pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4*2;
			rows[2]=pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4*2;
			rows[3]=pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4*2;
			for(int kx=0;kx<src->iw;++kx)
			{
				for(int kc=0;kc<3;++kc)
				{
					int
						NNN	=rows[3][kc+0*4*2+0],
						NNWW	=rows[2][kc-2*4*2+0],
						NNW	=rows[2][kc-1*4*2+0],
						NN	=rows[2][kc+0*4*2+0],
						NNE	=rows[2][kc+1*4*2+0],
						NNEE	=rows[2][kc+2*4*2+0],
						NWW	=rows[1][kc-2*4*2+0],
						NW	=rows[1][kc-1*4*2+0],
						N	=rows[1][kc+0*4*2+0],
						NE	=rows[1][kc+1*4*2+0],
						NEE	=rows[1][kc+2*4*2+0],
						NEEE	=rows[1][kc+3*4*2+0],
						NEEEE	=rows[1][kc+4*4*2+0],
						WWWW	=rows[0][kc-4*4*2+0],
						WWW	=rows[0][kc-3*4*2+0],
						WW	=rows[0][kc-2*4*2+0],
						W	=rows[0][kc-1*4*2+0],
						eNW	=rows[1][kc-1*4*2+1],
						eN	=rows[1][kc+0*4*2+1],
						eNE	=rows[1][kc+1*4*2+1],
						eWW	=rows[0][kc-2*4*2+1],
						eW	=rows[0][kc-1*4*2+1];
					int curr=rows[0][kc];
					int px=3*(W-WW)+WWW, py=3*(N-NN)+NNN;
					int g=N+W-NW;
					int p0=
						+coeffs[kc][0]*N
						+coeffs[kc][1]*W
						+coeffs[kc][2]*NE
						+coeffs[kc][3]*g
						+coeffs[kc][4]*py
						+coeffs[kc][5]*px
					;
					p0+=1<<18>>1;
					p0>>=18;
					int t0;
					t0=N;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][0]+=t0;
					t0=W;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][1]+=t0;
					t0=NE;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][2]+=t0;
					t0=g;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][3]+=t0;
					t0=py;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][4]+=t0;
					t0=px;	if(curr<p0)t0=-t0; if(curr==p0)t0=0; coeffs[kc][5]+=t0;
				}
				rows[0]+=4*2;
				rows[1]+=4*2;
				rows[2]+=4*2;
				rows[3]+=4*2;
			}
		}
#endif
	}
	free(pixels);
}


int crct_analysis(Image *src)
{
	long long counters[OCH_COUNT]={0};
	int prev[OCH_COUNT]={0};
	for(ptrdiff_t k=0, len=(ptrdiff_t)src->iw*src->ih*4;k<len;k+=4)
	{
		int
			r=src->data[k+0]<<2,
			g=src->data[k+1]<<2,
			b=src->data[k+2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
		counters[0]+=abs(r -prev[0]);
		counters[1]+=abs(g -prev[1]);
		counters[2]+=abs(b -prev[2]);
		counters[3]+=abs(rg-prev[3]);
		counters[4]+=abs(gb-prev[4]);
		counters[5]+=abs(br-prev[5]);
		prev[0]=r;
		prev[1]=g;
		prev[2]=b;
		prev[3]=rg;
		prev[4]=gb;
		prev[5]=br;
#ifdef ENABLE_EXTENDED_RCT
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0)\
	do\
	{\
		int a0=A0, b0=B0, c0=C0;\
		counters[IDXA]+=abs(a0-prev[IDXA]);\
		counters[IDXB]+=abs(b0-prev[IDXB]);\
		counters[IDXC]+=abs(c0-prev[IDXC]);\
		prev[IDXA]=a0;\
		prev[IDXB]=b0;\
		prev[IDXC]=c0;\
	}while(0)
		//r-(3*g+b)/4 = r-g-(b-g)/4
		//g-(3*r+b)/4 = g-r-(b-r)/4
		//b-(3*r+g)/4 = b-r-(g-r)/4
		UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, rg+(gb>>2), rg+(br>>2), br+(rg>>2));

		//r-(g+3*b)/4 = r-b-(g-b)/4
		//g-(r+3*b)/4 = g-b-(r-b)/4
		//b-(r+3*g)/4 = b-g-(r-g)/4
		UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, br+(gb>>2), gb+(br>>2), gb+(rg>>2));

		//r-(g+b)/2 = (r-g + r-b)/2
		//g-(r+b)/2 = (g-r + g-b)/2
		//b-(r+g)/2 = (b-r + b-g)/2
		UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, (rg-br)>>1, (gb-rg)>>1, (br-gb)>>1);
#undef  UPDATE
#endif
	}
	int bestrct=0;
	long long minerr=0;
	for(int kt=0;kt<RCT_COUNT;++kt)
	{
		const unsigned char *rct=rct_combinations[kt];
		long long currerr=
			+counters[rct[0]]
			+counters[rct[1]]
			+counters[rct[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			bestrct=kt;
		}
	}
	return bestrct;
}
void pred_l1crct(Image *src, int fwd)
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
	int weights[4][NPREDS]={0};
	int bufsize=(src->iw+8*2)*(int)sizeof(short[4*4*1]);//4 padded rows * 4 channels max * {pixels}
	short *pixels=(short*)malloc(bufsize);
	int invdist=((1<<16)+g_dist-1)/g_dist;
	if(fwd)
		src->rct=crct_analysis(src);
	const unsigned char *combination=rct_combinations[src->rct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*1,
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
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
				int
					NNNWWW		=rows[3][-3*4*1],
					NNNW		=rows[3][-1*4*1],
					NNN		=rows[3][+0*4*1],
					NNNE		=rows[3][+1*4*1],
					NNNEE		=rows[3][+2*4*1],
					NNNEEE		=rows[3][+3*4*1],
					NNNEEEE		=rows[3][+4*4*1],
					NNWWWW		=rows[2][-4*4*1],
					NNWWW		=rows[2][-3*4*1],
					NNWW		=rows[2][-2*4*1],
					NNW		=rows[2][-1*4*1],
					NN		=rows[2][+0*4*1],
					NNE		=rows[2][+1*4*1],
					NNEE		=rows[2][+2*4*1],
					NNEEE		=rows[2][+3*4*1],
					NNEEEE		=rows[2][+4*4*1],
					NWWWW		=rows[1][-4*4*1],
					NWWW		=rows[1][-3*4*1],
					NWW		=rows[1][-2*4*1],
					NW		=rows[1][-1*4*1],
					N		=rows[1][+0*4*1],
					NE		=rows[1][+1*4*1],
					NEE		=rows[1][+2*4*1],
					NEEE		=rows[1][+3*4*1],
					NEEEE		=rows[1][+4*4*1],
					NEEEEE		=rows[1][+5*4*1],
					NEEEEEE		=rows[1][+6*4*1],
					NEEEEEEE	=rows[1][+7*4*1],
					NEEEEEEEE	=rows[1][+8*4*1],
					WWWWWWWWW	=rows[0][-9*4*1],
					WWWWWWWW	=rows[0][-8*4*1],
					WWWWWWW		=rows[0][-7*4*1],
					WWWWWW		=rows[0][-6*4*1],
					WWWWW		=rows[0][-5*4*1],
					WWWW		=rows[0][-4*4*1],
					WWW		=rows[0][-3*4*1],
					WW		=rows[0][-2*4*1],
					W		=rows[0][-1*4*1];
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int *currw=weights[kc];
				int p0=0;
				for(int k=0;k<NPREDS;++k)
					p0+=currw[k]*preds[k];
				p0+=1LL<<L1SH>>1;
				p0>>=L1SH;
				int predc=p0;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
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
					if(fwd)
					{
						curr-=predc;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx+kc]=curr;
						
						curr=g_dist*curr+predc;
					}
					else
						curr=g_dist*src->data[idx+kc]+predc;
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

				//update
				int e=(curr>p0)-(curr<p0);//L1
			//	int e=curr-p0;//L2 (faster rise, worse steady state)
			//	currw[NPREDS]+=e;
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];

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
}

void pred_grfilt(Image *src, int fwd)
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
	int weights[4][NPREDS]={0};
	int invdist=((1<<16)+g_dist-1)/g_dist;
	int bufsize=(src->iw+8*2)*(int)sizeof(short[4*4*3]);//4 padded rows * 4 channels max * {pixels, error}
	short *pixels=(short*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL+4)%4)+8)*4*3,
			pixels+((src->iw+16LL)*((ky-1LL+4)%4)+8)*4*3,
			pixels+((src->iw+16LL)*((ky-2LL+4)%4)+8)*4*3,
			pixels+((src->iw+16LL)*((ky-3LL+4)%4)+8)*4*3,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				if(!src->depth[kc])
					continue;
				int
					NNNWWW		=rows[3][0-3*4*3],
					NNNW		=rows[3][0-1*4*3],
					NNN		=rows[3][0+0*4*3],
					NNNE		=rows[3][0+1*4*3],
					NNNEE		=rows[3][0+2*4*3],
					NNNEEE		=rows[3][0+3*4*3],
					NNNEEEE		=rows[3][0+4*4*3],
					NNWWWW		=rows[2][0-4*4*3],
					NNWWW		=rows[2][0-3*4*3],
					NNWW		=rows[2][0-2*4*3],
					NNW		=rows[2][0-1*4*3],
					NN		=rows[2][0+0*4*3],
					NNE		=rows[2][0+1*4*3],
					NNEE		=rows[2][0+2*4*3],
					NNEEE		=rows[2][0+3*4*3],
					NNEEEE		=rows[2][0+4*4*3],
					NWWWW		=rows[1][0-4*4*3],
					NWWW		=rows[1][0-3*4*3],
					NWW		=rows[1][0-2*4*3],
					NW		=rows[1][0-1*4*3],
					N		=rows[1][0+0*4*3],
					NE		=rows[1][0+1*4*3],
					NEE		=rows[1][0+2*4*3],
					NEEE		=rows[1][0+3*4*3],
					NEEEE		=rows[1][0+4*4*3],
					NEEEEE		=rows[1][0+5*4*3],
					NEEEEEE		=rows[1][0+6*4*3],
					NEEEEEEE	=rows[1][0+7*4*3],
					NEEEEEEEE	=rows[1][0+8*4*3],
					WWWWWWWWW	=rows[0][0-9*4*3],
					WWWWWWWW	=rows[0][0-8*4*3],
					WWWWWWW		=rows[0][0-7*4*3],
					WWWWWW		=rows[0][0-6*4*3],
					WWWWW		=rows[0][0-5*4*3],
					WWWW		=rows[0][0-4*4*3],
					WWW		=rows[0][0-3*4*3],
					WW		=rows[0][0-2*4*3],
					W		=rows[0][0-1*4*3],
					eNE		=rows[1][1+1*4*3],
					eNEE		=rows[1][1+2*4*3],
					eNEEE		=rows[1][1+3*4*3],
					eW		=rows[0][1-1*4*3],
					e2NE		=rows[1][2+1*4*3],
					e2NEE		=rows[1][2+2*4*3],
					e2NEEE		=rows[1][2+3*4*3],
					e2W		=rows[0][2-1*4*3];
				
				int pred=0;
				int curr=src->data[idx];

				src->data[idx]=eW>>6;
				//if(fwd)
				//{
				//	int error=curr-pred;
				//	error<<=32-src->depth[kc];
				//	error>>=32-src->depth[kc];
				//	src->data[idx]=error;
				//}
				//else
				//{
				//	curr+=pred;
				//	curr<<=32-src->depth[kc];
				//	curr>>=32-src->depth[kc];
				//	src->data[idx]=curr;
				//}
				rows[0][0]=curr;
				{
					int error=abs(curr);
					//int error=curr<<1^curr>>31;
					rows[0][1]=(eW+(eW<eNE?eW:eNE)+(error<<6)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				//	rows[0][1]=(2*eW+error+eNEEE)>>2;

					rows[0][2]=(e2W+(e2W<e2NE?e2W:e2NE)+eW+(e2NEE>e2NEEE?e2NEE:e2NEEE))>>2;
				}
			}
		}
	}
	free(pixels);
}
