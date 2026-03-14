#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>//log2
#include<immintrin.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#elif defined __GNUC__
//#include<x86intrin.h>
//#endif
static const char file[]=__FILE__;


//	#define ESTIMATE_SIZE
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
#define L1SH 21
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
enum
{
#define PRED(WEIGHT, EXPR) +1
	NPREDS=PREDLIST,
#undef  PRED
};
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
	int psize=(src->iw+16*2)*(int)sizeof(short[4*4*1]);//4 padded rows * 4 channels max * {pixels}
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
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
			pixels+(((src->iw+16LL*2)*((ky-0LL+4)%4)+16)*4-1)*1,
			pixels+(((src->iw+16LL*2)*((ky-1LL+4)%4)+16)*4-1)*1,
			pixels+(((src->iw+16LL*2)*((ky-2LL+4)%4)+16)*4-1)*1,
			pixels+(((src->iw+16LL*2)*((ky-3LL+4)%4)+16)*4-1)*1,
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
				int predc=1<<L1SH>>1;
				for(int k=0;k<NPREDS;++k)
					predc+=currw[k]*preds[k];
				predc>>=L1SH;
			//	p0-=p0>>31;//X  deadzone bad with advanced pred
				int p0=predc;
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
						//curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
						curr=(curr*invdist>>16)-(curr>>31);//curr/=g_dist
						src->data[idx]=curr;
					}
					curr=g_dist*curr+predc;
					CLAMP2(curr, amin[kc], amax[kc]);
					if(!fwd)
						src->data[idx]=curr;
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

static void noise_analysis(Image *src)
{
#define LAGLIST\
	LAG(leak[kc][0]>>16)\
	LAG(leak[kc][1]>>16)\
	LAG(leak[kc][2]>>16)\
	LAG(leak[kc][3]>>16)\
	LAG(leak[kc][4]>>16)\
	LAG(leak[kc][5]>>16)\
	LAG(leak[kc][6]>>16)\
	LAG(leak[kc][7]>>16)\

	enum
	{
#define LAG(...) +1
		NLAGS=LAGLIST
#undef  LAG
	};
	ptrdiff_t k=0, size=0;
	int *ptr=0;
	uint64_t ctr[4][NLAGS]={0};
	int32_t leak[4][NLAGS]={0};
	volatile double t=0;

	t=time_sec();
	for(k=0, size=(ptrdiff_t)src->iw*src->ih, ptr=src->data;k<size;k+=4, ptr+=4)
	{
		int kc=0;

		for(kc=0;kc<4;++kc)
		{
			if(!src->depth[kc])
				continue;
			int
				WWWW	=ptr[kc-4*4],
				WWW	=ptr[kc-3*4],
				WW	=ptr[kc-2*4],
				W	=ptr[kc-1*4],
				curr	=ptr[kc+0*4];
			int j=0;
#define LAG(E) ctr[kc][j++]+=abs(curr-(E));
			j=0;
			LAGLIST
#undef  LAG
			leak[kc][0]+=((curr<<16)-leak[kc][0])>>0;
			leak[kc][1]+=((curr<<16)-leak[kc][1])>>1;
			leak[kc][2]+=((curr<<16)-leak[kc][2])>>2;
			leak[kc][3]+=((curr<<16)-leak[kc][3])>>3;
			leak[kc][4]+=((curr<<16)-leak[kc][4])>>4;
			leak[kc][5]+=((curr<<16)-leak[kc][5])>>5;
			leak[kc][6]+=((curr<<16)-leak[kc][6])>>6;
			leak[kc][7]+=((curr<<16)-leak[kc][7])>>7;
		}
	}
	t=time_sec()-t;
	{
		char buf[4096]={0};
		int printed=0;
		static const char *lagnames[]=
		{
#define LAG(E) #E,
			LAGLIST
#undef  LAG
		};

		for(int k=0;k<NLAGS;++k)
			printed+=snprintf(buf+printed, sizeof(buf)-1-printed, "%14lld %14lld %14lld  %s\n"
				, ctr[0][k]
				, ctr[1][k]
				, ctr[2][k]
				, lagnames[k]
			);
		printed+=snprintf(buf+printed, sizeof(buf)-1-printed, "%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
			, t
			, (3.*src->iw*src->ih)/(t*1024*1024)
			, (t*1024*1024*1000)/(3.*src->iw*src->ih)
		);
		copy_to_clipboard(buf, printed);
		messagebox(MBOX_OK, "Copied", "%s", buf);
	}
}
static void crct2_analysis(Image *src, int32_t *alphas)
{
	int maxdepth=src->depth[0];
	if(maxdepth<src->depth[1])maxdepth=src->depth[1];
	if(maxdepth<src->depth[2])maxdepth=src->depth[2];
	if(maxdepth<src->depth[3])maxdepth=src->depth[3];
	int hsize=(int)sizeof(int32_t[6])<<maxdepth;
	int32_t *hists=(int32_t*)malloc(hsize);
	if(!hists)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hists, 0, hsize);
	int prev[6]={0};
	int half[6]=
	{
		1<<src->depth[0]>>1,
		1<<src->depth[1]>>1,
		1<<src->depth[2]>>1,
		1<<maxdepth>>1,
		1<<maxdepth>>1,
		1<<maxdepth>>1,
	};
	int mask[6]=
	{
		(1<<src->depth[0])-1,
		(1<<src->depth[1])-1,
		(1<<src->depth[2])-1,
		(1<<maxdepth)-1,
		(1<<maxdepth)-1,
		(1<<maxdepth)-1,
	};
	for(ptrdiff_t k=0, size=(ptrdiff_t)4*src->iw*src->ih;k<size;k+=4)
	{
		int r=src->data[k+0];
		int g=src->data[k+1];
		int b=src->data[k+2];
		int rg=r-g, gb=g-b, br=b-r;
		++hists[0<<maxdepth|((r -prev[0]+half[0])&mask[0])];
		++hists[1<<maxdepth|((g -prev[1]+half[1])&mask[1])];
		++hists[2<<maxdepth|((b -prev[2]+half[2])&mask[2])];
		++hists[3<<maxdepth|((rg-prev[3]+half[3])&mask[3])];
		++hists[4<<maxdepth|((gb-prev[4]+half[4])&mask[4])];
		++hists[5<<maxdepth|((br-prev[5]+half[5])&mask[5])];
		prev[0]=r;
		prev[1]=g;
		prev[2]=b;
		prev[3]=rg;
		prev[4]=gb;
		prev[5]=br;
	}
	double csizes[6]={0};
	for(int kc=0;kc<6;++kc)
	{
		int32_t *currhist=hists+((ptrdiff_t)kc<<maxdepth);
		int32_t sum=0;
		for(int ks=0;ks<=mask[kc];++ks)
			sum+=currhist[ks];
		if(!sum)
			continue;
		double invsum=1./sum, e=0;
		for(int ks=0;ks<=mask[kc];++ks)
		{
			int32_t freq=currhist[ks];
			if(freq)
				e-=freq*log2(freq*invsum);
		}
		csizes[kc]=e/8;
	}
	free(hists);
	
	//alphas decrease with correlation
	double a;
	
//	a=(csizes[3]+csizes[0])/(csizes[1]+csizes[0]);//(rg+r)/(g+r)
	a=csizes[3]/csizes[1];//rg/g
	a*=a;
	a*=a;
	a*=a;
	a*=a;
	CLAMP2(a, 0, 1);
	alphas[0]=(int32_t)CVTFP64_I64(a*0x10000);
	
//	a=(csizes[5]+csizes[0])/(csizes[2]+csizes[0]);//(br+r)/(b+r)
	a=csizes[5]/csizes[2];//br/b
	a*=a;
	a*=a;
	a*=a;
	a*=a;
	CLAMP2(a, 0, 1);
	alphas[1]=(int32_t)CVTFP64_I64(a*0x10000);
	
//	a=(csizes[4]+csizes[1])/(csizes[2]+csizes[1]);//(gb+g)/(b+g)
	a=csizes[4]/csizes[2];//gb/b
	a*=a;
	a*=a;
	a*=a;
	a*=a;
	CLAMP2(a, 0, 1);
	alphas[2]=(int32_t)CVTFP64_I64(a*0x10000);
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
	int rmin[]=
	{
		amin[0]/g_dist,
		amin[1]/g_dist,
		amin[2]/g_dist,
		amin[3]/g_dist,
	};
	int rmax[]=
	{
		amax[0]/g_dist,
		amax[1]/g_dist,
		amax[2]/g_dist,
		amax[3]/g_dist,
	};
	int invdist=((1<<16)+g_dist-1)/g_dist;

	enum
	{
		MIXPREDS=4,

		SHIFT=18+6,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=1,
	};
	ALIGN(16) int32_t coeffs[4][MIXPREDS]={0}, bias[4]={1<<SHIFT>>1}, estims[MIXPREDS]={0};

	/*
	cache-friendly layout:
	...
	A(NNN NN N C)WW

	Y(NNN NN N C)W
	U(NNN NN N C)W
	V(NNN NN N C)W
	A(NNN NN N C)W

	Y(NNN NN N C)curr
	...
	*/
	int psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	int32_t *pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	FILLMEM((int32_t*)coeffs, (1<<SHIFT)/MIXPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[3]=bias[2]=bias[1]=bias[0];
	memset(pixels, 0, psize);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				int32_t
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL];
				//	eN	=rows[1][1+0*NCH*NROWS*NVAL],
				//	eNE	=rows[1][1+1*NCH*NROWS*NVAL],
				//	eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int curr=src->data[idx];
				int j=0;

				//mix 2
#if 0
				estims[j++]=N;
				estims[j++]=W;
#endif

				//mix 4
#if 1
				//		NN
				//	NW	N	NE
				//	W	?		216 MB/s  4.49 ms/MB  i7-13700KF

				estims[j++]=(W+NE)>>1;
				estims[j++]=(WWWW+WWW+NNN+NEEE)>>2;
				estims[j++]=NW-((NN+WW)>>1);
				estims[j++]=N+W-NW;

				//estims[j++]=4*(W+NE);
				//estims[j++]=2*(WWWW+WWW+NNN+NEEE);
				//estims[j++]=4*(2*NW-(NN+WW));
				//estims[j++]=8*(N+W-NW);

				//estims[j++]=W;
				//estims[j++]=NE;
				//estims[j++]=2*N-NN;
				//estims[j++]=N+W-NW;
#endif

				//mix 8 - c32
#if 0
				//					NNN
				//					NN	NNE
				//				NW	N	NE	NEE	NEEE	NEEEE
				//	WWWW	WWW	WW	W	?
				estims[j++]=N;
				estims[j++]=W;
				estims[j++]=3*(N-NN)+NNN;
				estims[j++]=3*(W-WW)+WWW;
				estims[j++]=W+NE-N;
				estims[j++]=(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4;
				estims[j++]=N+W-NW;
				estims[j++]=N+NE-NNE;
#endif
#if 0
				for(int k=0;k<10;++k)//cheat
				{
					int p1=(int)((bias[kc]
						+(int64_t)coeffs[kc][0]*estims[0]
						+(int64_t)coeffs[kc][1]*estims[1]
						+(int64_t)coeffs[kc][2]*estims[2]
						+(int64_t)coeffs[kc][3]*estims[3]
						+(int64_t)coeffs[kc][4]*estims[4]
					)>>SHIFT);
					int e=(curr>p1)-(curr<p1);//L1
					bias[kc]+=e;
					coeffs[kc][0]+=e*estims[0];
					coeffs[kc][1]+=e*estims[1];
					coeffs[kc][2]+=e*estims[2];
					coeffs[kc][3]+=e*estims[3];
					coeffs[kc][4]+=e*estims[4];
				}
#endif
				int p1=(int)((bias[kc]
					+(int64_t)coeffs[kc][0]*estims[0]
					+(int64_t)coeffs[kc][1]*estims[1]
					+(int64_t)coeffs[kc][2]*estims[2]
					+(int64_t)coeffs[kc][3]*estims[3]
				)>>SHIFT);
				int pred=p1;

				//cheat
#if 0
				//int E=kx<src->iw-1?src->data[idx+4]:0, EE=kx<src->iw-2?src->data[idx+2*4]:0;
				//int S=ky<src->ih-1?src->data[idx+4*src->iw]:0, SS=ky<src->ih-2?src->data[idx+2*4*src->iw]:0;
				//int e1=(W+E)>>1, e2=(N+S)>>1;
				//int e1=W+E-((WW+EE)>>1);
				//int e2=N+S-((NN+SS)>>1);
				//pred=(W+E)>>1;
				//int e1=W;
				//int e2=E;
				//int e1=(N+W)>>1;
				//int e2=N+W-NW;
				//pred=32*(abs(curr-e1)<abs(curr-e2));
				//pred=abs(curr-e1)<abs(curr-e2)?e1:e2;
#endif
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
				
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)pred;
						//curr=curr<0?-(-curr>>1):curr>>1;//IMG0008 d3  17.9% 33.3 dB  19.19% 34.0 dB

						//curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
						curr=(curr*invdist>>16)-(curr>>31);//curr/=g_dist
						CLAMP2(curr, rmin[kc], rmax[kc]);
						src->data[idx]=curr;
					}
					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);
					if(!fwd)
						src->data[idx]=curr;
				}
				else
				{
					if(fwd)
					{
						int error=curr-pred;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx]=error;
						//src->data[idx]=32*(abs(curr-e1)<abs(curr-e2));
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx]=curr;
					}
				}
				rows[0][0]=curr;
				//rows[0][1]=curr-p1;

				//150 MB/s  6.65 ms/MB  i5-1145G7
#if 0
				{
					__m128i p=_mm_load_si128((__m128i*)estims);
					__m128i c=_mm_load_si128((__m128i*)coeffs[kc]);
					p=_mm_sign_epi32(p, _mm_set1_epi32(curr-p1));
					c=_mm_add_epi32(c, p);
					_mm_store_si128((__m128i*)coeffs[kc], c);
				}
#endif

				//147 MB/s  6.77 ms/MB  i5-1145G7
#if 1
				//int e=(curr>p1)-(curr<p1);//L1
				int e=((curr-p1)>>31)-((p1-curr)>>31);
				//int e=curr-p1; CLAMP2(e, -1, 1);//jump?
				bias[kc]+=e;
				coeffs[kc][0]+=(int16_t)((int16_t)e*(int16_t)estims[0]);//casts prevent pmulld
				coeffs[kc][1]+=(int16_t)((int16_t)e*(int16_t)estims[1]);
				coeffs[kc][2]+=(int16_t)((int16_t)e*(int16_t)estims[2]);
				coeffs[kc][3]+=(int16_t)((int16_t)e*(int16_t)estims[3]);
#endif
			}
		}
	}
	_mm_free(pixels);
}
void pred_mixN_crct2(Image *src, int fwd)
{
	enum
	{
		MIXPREDS=4,
		SHIFT=18,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,
	};
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
	int invdist=((1<<16)+g_dist-1)/g_dist;

	ALIGN(16) int coeffs[4][MIXPREDS]={0}, estims[MIXPREDS]={0}, coeff2[4][MIXPREDS]={0}, estim2[MIXPREDS]={0};

	int psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	int16_t *pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));

#ifdef _DEBUG
	static//either  reversible  or  multithreaded batch test
#endif
	int32_t alphas[3]={0};

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	FILLMEM((int32_t*)coeffs, (1<<SHIFT)/MIXPREDS, sizeof(coeffs), sizeof(int32_t));
	memset(pixels, 0, psize);
	if(fwd)
		crct2_analysis(src, alphas);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			int offset=0;
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				int16_t
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL];
				int16_t
					aNNN	=rows[3][1+0*NCH*NROWS*NVAL],
					aNNWW	=rows[2][1-2*NCH*NROWS*NVAL],
					aNNW	=rows[2][1-1*NCH*NROWS*NVAL],
					aNN	=rows[2][1+0*NCH*NROWS*NVAL],
					aNNE	=rows[2][1+1*NCH*NROWS*NVAL],
					aNNEE	=rows[2][1+2*NCH*NROWS*NVAL],
					aNWW	=rows[1][1-2*NCH*NROWS*NVAL],
					aNW	=rows[1][1-1*NCH*NROWS*NVAL],
					aN	=rows[1][1+0*NCH*NROWS*NVAL],
					aNE	=rows[1][1+1*NCH*NROWS*NVAL],
					aNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					aNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					aNEEEE	=rows[1][1+4*NCH*NROWS*NVAL],
					aWWWW	=rows[0][1-4*NCH*NROWS*NVAL],
					aWWW	=rows[0][1-3*NCH*NROWS*NVAL],
					aWW	=rows[0][1-2*NCH*NROWS*NVAL],
					aW	=rows[0][1-1*NCH*NROWS*NVAL];
				int curr=src->data[idx];
				int *weights=coeffs[kc];
				int vmax=N, vmin=W;

				int j=0;
				//mix 4
				//		NN
				//	NW	N	NE
				//	W	?		133 MB/s  7.47 ms/MB  4T
				j=0;
				estims[j++]=W;
				estims[j++]=N+W-NW;
				estims[j++]=2*N-NN;
				estims[j++]=NE;
				int pred=((1<<SHIFT>>1)
					+weights[0]*estims[0]
					+weights[1]*estims[1]
					+weights[2]*estims[2]
					+weights[3]*estims[3]
				)>>SHIFT, pred2=0, p2=0;
				int p1=pred;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
				if(kc)
				{
					//pred2=aN+aW-aNW;
					//vmax=aN, vmin=aW;
					//if(aN<aW)vmin=aN, vmax=aW;
					//CLAMP2(pred2, vmin, vmax);

				//	pred2=0;
#if 1
					j=0;
					estim2[j++]=aW;
					estim2[j++]=aN+aW-aNW;
					estim2[j++]=2*aN-aNN;
					estim2[j++]=aNE;
					pred2=((1<<SHIFT>>1)
						+coeff2[kc][0]*estim2[0]
						+coeff2[kc][1]*estim2[1]
						+coeff2[kc][2]*estim2[2]
						+coeff2[kc][3]*estim2[3]
					)>>SHIFT;
					p2=pred2;
					vmax=aN, vmin=aW;
					if(aN<aW)vmin=aN, vmax=aW;
					if(vmin>aNE)vmin=aNE;
					if(vmax<aNE)vmax=aNE;
					if(vmin>aNEEE)vmin=aNEEE;
					if(vmax<aNEEE)vmax=aNEEE;
					CLAMP2(pred2, vmin, vmax);
#endif

					pred+=offset;
				//	int p=(pred2-pred)*(0x10000-alphas[kc-1]);
					int p=(pred2-pred)*alphas[kc-1];
					pred+=p<0?-(-p>>16):p>>16;
					CLAMP2(pred, amin[kc], amax[kc]);
				}
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)pred;
						//curr=curr<0?-(-curr>>1):curr>>1;//IMG0008 d3  17.9% 33.3 dB  19.19% 34.0 dB

						//curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
						curr=(curr*invdist>>16)-(curr>>31);//curr/=g_dist
						src->data[idx]=curr;
					}
					curr=g_dist*curr+(int)pred;
					CLAMP2(curr, amin[kc], amax[kc]);
					if(!fwd)
						src->data[idx]=curr;
				}
				else
				{
					if(fwd)
					{
						int error=curr-pred;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx]=error;
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx]=curr;
					}
				}
				rows[0][1]=curr;
				curr-=offset;
				rows[0][0]=curr;
				int e=(curr>p1)-(curr<p1);//L1
				weights[0]+=e*estims[0];
				weights[1]+=e*estims[1];
				weights[2]+=e*estims[2];
				weights[3]+=e*estims[3];
				if(kc)
				{
#if 1
					int acurr=rows[0][1];
					int e=(acurr>p2)-(acurr<p2);
					coeff2[kc][0]+=e*estim2[0];
					coeff2[kc][1]+=e*estim2[1];
					coeff2[kc][2]+=e*estim2[2];
					coeff2[kc][3]+=e*estim2[3];
#endif
				}
				else
					offset=curr;
			}
		}
	}
	_mm_free(pixels);
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
	enum
	{
		MIXPREDS=4,
		SHIFT=18,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,
	};
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
#ifdef ESTIMATE_SIZE
#define PREDBITS 1
#define NCTX 16
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int hstart[]=
	{
		0,
		nlevels[0],
		nlevels[0]+nlevels[1],
		nlevels[0]+nlevels[1]+nlevels[2],
		nlevels[0]+nlevels[1]+nlevels[2]+nlevels[3],
	};
	int hsize=sizeof(int[NCTX])*hstart[4]<<PREDBITS;
	int *hist=(int*)malloc(hsize);
	int esize=(src->iw+8*2)*(int)sizeof(short[4*4*1]);//4 padded rows * 4 channels max * {errors}
	short *ebuf=(short*)malloc(esize);
	if(!hist||!ebuf)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, hsize);
	memset(ebuf, 0, esize);
#endif
	int32_t weights[4][NPREDS]={0};
	int psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	int16_t *pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
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
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			//pixels+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*1,
			//pixels+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*1,
			//pixels+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*1,
			//pixels+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*1,
		};
#ifdef ESTIMATE_SIZE
		short *erows[]=
		{
			ebuf+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*1,
			ebuf+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*1,
			ebuf+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*1,
			ebuf+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*1,
		};
#endif
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
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
#ifdef ESTIMATE_SIZE
				++erows[0];
				++erows[1];
				++erows[2];
				++erows[3];
#endif
				if(!src->depth[kc])
					continue;
				int
					NNNWWW		=rows[3][0-3*NCH*NROWS*NVAL],
					NNNW		=rows[3][0-1*NCH*NROWS*NVAL],
					NNN		=rows[3][0+0*NCH*NROWS*NVAL],
					NNNE		=rows[3][0+1*NCH*NROWS*NVAL],
					NNNEE		=rows[3][0+2*NCH*NROWS*NVAL],
					NNNEEE		=rows[3][0+3*NCH*NROWS*NVAL],
					NNNEEEE		=rows[3][0+4*NCH*NROWS*NVAL],
					NNWWWW		=rows[2][0-4*NCH*NROWS*NVAL],
					NNWWW		=rows[2][0-3*NCH*NROWS*NVAL],
					NNWW		=rows[2][0-2*NCH*NROWS*NVAL],
					NNW		=rows[2][0-1*NCH*NROWS*NVAL],
					NN		=rows[2][0+0*NCH*NROWS*NVAL],
					NNE		=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE		=rows[2][0+2*NCH*NROWS*NVAL],
					NNEEE		=rows[2][0+3*NCH*NROWS*NVAL],
					NNEEEE		=rows[2][0+4*NCH*NROWS*NVAL],
					NWWWW		=rows[1][0-4*NCH*NROWS*NVAL],
					NWWW		=rows[1][0-3*NCH*NROWS*NVAL],
					NWW		=rows[1][0-2*NCH*NROWS*NVAL],
					NW		=rows[1][0-1*NCH*NROWS*NVAL],
					N		=rows[1][0+0*NCH*NROWS*NVAL],
					NE		=rows[1][0+1*NCH*NROWS*NVAL],
					NEE		=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE		=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE		=rows[1][0+4*NCH*NROWS*NVAL],
					NEEEEE		=rows[1][0+5*NCH*NROWS*NVAL],
					NEEEEEE		=rows[1][0+6*NCH*NROWS*NVAL],
					NEEEEEEE	=rows[1][0+7*NCH*NROWS*NVAL],
					NEEEEEEEE	=rows[1][0+8*NCH*NROWS*NVAL],
					WWWWWWWWW	=rows[0][0-9*NCH*NROWS*NVAL],
					WWWWWWWW	=rows[0][0-8*NCH*NROWS*NVAL],
					WWWWWWW		=rows[0][0-7*NCH*NROWS*NVAL],
					WWWWWW		=rows[0][0-6*NCH*NROWS*NVAL],
					WWWWW		=rows[0][0-5*NCH*NROWS*NVAL],
					WWWW		=rows[0][0-4*NCH*NROWS*NVAL],
					WWW		=rows[0][0-3*NCH*NROWS*NVAL],
					WW		=rows[0][0-2*NCH*NROWS*NVAL],
					W		=rows[0][0-1*NCH*NROWS*NVAL],
					aNW		=rows[1][1-1*NCH*NROWS*NVAL],
					aN		=rows[1][1+0*NCH*NROWS*NVAL],
					aNE		=rows[1][1+1*NCH*NROWS*NVAL],
					aNEEE		=rows[1][1+3*NCH*NROWS*NVAL],
					aW		=rows[0][1-1*NCH*NROWS*NVAL];
#ifdef ESTIMATE_SIZE
				int
					eNEE		=erows[1][+2*4*1],
					eNEEE		=erows[1][+3*4*1],
					eW		=erows[0][-1*4*1];
				int ctx=FLOOR_LOG2(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
#endif
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int *currw=weights[kc];
				int p0=1LL<<L1SH>>1;
				for(int k=0;k<NPREDS;++k)
					p0+=currw[k]*preds[k];
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
#if 0
					vmax=aN, vmin=aW;
					if(aN<aW)vmin=aN, vmax=aW;
					if(vmin>aNE)vmin=aNE;
					if(vmax<aNE)vmax=aNE;
					if(vmin>aNEEE)vmin=aNEEE;
					if(vmax<aNEEE)vmax=aNEEE;
					if(vmin>aNW)vmin=aNW;
					if(vmax<aNW)vmax=aNW;
					CLAMP2(predc, vmin, vmax);
#endif
					CLAMP2(predc, amin[kc], amax[kc]);
				}

				int curr=yuv[kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=predc;
					//	curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						curr=(curr*invdist>>16)-(curr>>31);
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
#ifdef ESTIMATE_SIZE
				int e2=curr-predc;
				e2<<=32-src->depth[kc];//MA
				e2>>=32-src->depth[kc];
				//if((unsigned)(hstart[4]*((ctx<<PREDBITS)+((predc+(1<<src->depth[kc]>>1))>>(src->depth[kc]-PREDBITS)))+hstart[kc]+e2+(1<<src->depth[kc]>>1))>=(unsigned)hsize)
				//	LOG_ERROR("");
				++hist[hstart[4]*((ctx<<PREDBITS)+((predc+(1<<src->depth[kc]>>1))>>(src->depth[kc]-PREDBITS)))+hstart[kc]+e2+(1<<src->depth[kc]>>1)];
				e2=e2<<1^e2>>31;
				erows[0][0]=(2*eW+(e2<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
#endif
				rows[0][1]=curr;
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
	_mm_free(pixels);
#ifdef ESTIMATE_SIZE
	if(loud_transforms)
	{
		double csize=0, overhead=0;
		for(int kctx=0;kctx<NCTX<<PREDBITS;++kctx)
		{
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				int *hcurr=hist+hstart[4]*kctx+hstart[kc];
				int sum=0;
				for(int ks=0;ks<nlevels[kc];++ks)
					sum+=hcurr[ks];
				if(!sum)
					continue;

				//simulate nonzero-bin tracking guard
				int count=0;
				for(int ks=0;ks<nlevels[kc];++ks)
					count+=hcurr[ks]!=0;
				for(int ks=0;ks<nlevels[kc];++ks)
				{
					int freq=hcurr[ks];
					if(freq)
					{
						int prob=(int)((long long)freq*(0x1000LL-count)/sum)+1;
						csize-=freq*log2(prob*(1./0x1000));
					}
				}
			
				//overhead size
				const int probbits=12;
				int cdfW=0;
				int sum2=0;
				int codelen=probbits+1, CDFlevels=1<<probbits;
				int nlevels2=1<<src->depth[kc], half2=nlevels2>>1, mask2=nlevels2-1;
				for(int ks=0, ks2=0;ks<nlevels2;++ks)//calc overhead size
				{
					int sym=((ks>>1^-(ks&1))+half2)&mask2;
					int freq=hcurr[sym];
					int cdf=sum2*((1ULL<<probbits)-count)/sum+ks2;
					ks2+=freq!=0;
					int csym=cdf-cdfW;
					if(ks&&CDFlevels)//CDF[0] is always zero
					{
						//GR
						int nbypass=FLOOR_LOG2(CDFlevels);
						if(ks>1)
							nbypass-=7;
						if(nbypass<0)
							nbypass=0;
						overhead+=(csym>>nbypass)+1+nbypass;
					}
					CDFlevels-=csym;
					cdfW=cdf;
					sum2+=freq;
				}

				//double norm=1./sum;
				//for(int ks=0;ks<nlevels[kc];++ks)
				//{
				//	int freq=hcurr[ks];
				//	if(freq)
				//		csize-=freq*log2(freq*norm);
				//}
			}
		}
		csize/=8;
		overhead/=8;
		set_window_title("%10.2lf + %10.2lf = %10.2lf"
			, csize
			, overhead
			, csize+overhead
		);
	}
	free(ebuf);
	free(hist);
#endif
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
	int psize=(src->iw+8*2)*(int)sizeof(short[4*4*3]);//4 padded rows * 4 channels max * {pixels, error}
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
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

void pred_adaquant(Image *src, int fwd)
{
	enum
	{
		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,

		NPREDS=4,
		SHIFT=18,
	};

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
	int32_t weights[NCH][NPREDS]={0}, estims[NPREDS]={0};
	int32_t LPF[4]={0};
	int psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	int16_t *pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		int drift[4]={0};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				int16_t
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					eN	=rows[1][1+0*NCH*NROWS*NVAL],
					eNE	=rows[1][1+1*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
			//	int sh=0;
				int pred=1<<SHIFT>>1, p1, j=0, curr, error;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				j=0;
				estims[j++]=W;
				estims[j++]=N+W-NW;
				estims[j++]=2*N-NN;
				estims[j++]=NE;
				pred=(
					+weights[kc][0]*estims[0]
					+weights[kc][1]*estims[1]
					+weights[kc][2]*estims[2]
					+weights[kc][3]*estims[3]
				)>>SHIFT;
				p1=pred;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

			//	int qden=eW+eNE+1;
				int qden=eW+1;
				if(qden>(1<<g_dist))//just to decrease penalty for easy regions
					qden=(1<<g_dist);
			//	LPF[kc]+=(eW*eW-LPF[kc])>>3;
			//	sh=FLOOR_LOG2(eW*g_dist+1);
			//	if(sh>g_dist)
			//		sh=g_dist;
				//if(kc)//lighter quantization on chroma
				//	CLAMP2(sh, g_dist>>3, g_dist>>1);
				//else
				//	CLAMP2(sh, g_dist>>2, g_dist);

				curr=src->data[idx];
				if(fwd)
				{
					curr-=pred;
					//curr<<=32-src->depth[kc];
					//curr>>=32-src->depth[kc];
					//curr=curr<0?-(-curr>>sh):curr>>sh;
					curr/=qden;
					error=curr;
					//curr<<=sh;
					curr*=qden;
					curr+=pred;
					CLAMP2(curr, amin[kc], amax[kc]);
					src->data[idx]=error;
				}
				else
				{
					error=curr;
					//curr<<=sh;
					curr*=qden;
					curr+=pred;
					CLAMP2(curr, amin[kc], amax[kc]);
					src->data[idx]=curr;
				}
				rows[0][0]=curr;
			//	rows[0][1]=error<<1^error>>31;
			//	rows[0][1]=eW+((((error<<1^error>>31)<<g_dist)+eN-2*eW)>>3);
				rows[0][1]=(2*eW+((error<<1^error>>31)<<g_dist)+eNEEE)>>2;

				int e=(curr>p1)-(curr<p1);//L1
				weights[kc][0]+=e*estims[0];
				weights[kc][1]+=e*estims[1];
				weights[kc][2]+=e*estims[2];
				weights[kc][3]+=e*estims[3];
			}
		}
	}
	_mm_free(pixels);
}

void pred_gray(Image *src, int fwd)
{
	int half[]=
	{
		1<<src->depth[0]>>1,
		1<<src->depth[1]>>1,
		1<<src->depth[2]>>1,
		1<<src->depth[3]>>1,
	};
	if(fwd)
	{
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				for(int kc=0;kc<4;++kc, ++idx)
				{
					if(!src->depth[kc])
						continue;
					int val=src->data[idx]+half[kc];
					val^=val>>1;
					src->data[idx]=val-half[kc];
				}
			}
		}
	}
	else
	{
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				for(int kc=0;kc<4;++kc, ++idx)
				{
					if(!src->depth[kc])
						continue;
					int val=src->data[idx]+half[kc];
					val^=val>>16;
					val^=val>> 8;
					val^=val>> 4;
					val^=val>> 2;
					val^=val>> 1;
					src->data[idx]=val-half[kc];
				}
			}
		}
	}
}

static void init_weights(int32_t *weights, int npreds, int sh)
{
	for(int k=0;k<npreds;++k)
		weights[k]=(1<<sh)/npreds;
}
#define LOAD(X, Y) (uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?src->data[4*(src->iw*(ky+(Y))+kx+(X))+kc]:0
static void awav_queen(Image *src, int fwd, int kc, int iw, int ih)
{
	enum
	{
	//	NPRED=4,
	//	SHIFT=15,
	//	NPRED=12,
	//	SHIFT=19,
		NPRED=14,
		SHIFT=20,
	};
	int32_t weights[NPRED]={0};
	int amin=-(1<<src->depth[kc]>>1), amax=(1<<src->depth[kc]>>1)-1;

	init_weights(weights, NPRED, SHIFT);
	for(int ky=0;ky<ih;ky+=2)
	{
		for(int kx=1;kx<iw;kx+=2)
		{
			int
				NNN	=LOAD(+0, -3),
				WWW	=LOAD(-3, +0),
				SSS	=LOAD(+0, +3),
				EEE	=LOAD(+3, +0),
				NW	=LOAD(-1, -1),
				NE	=LOAD(+1, -1),
				SW	=LOAD(-1, +1),
				SE	=LOAD(+1, +1),
				N	=LOAD(+0, -1),
				W	=LOAD(-1, +0),
				S	=LOAD(+0, +1),
				E	=LOAD(+1, +0);
			int estim[]=
			{
				N,
				W,
				S,
				E,
				N+W-NW,
				N+E-NE,
				S+W-SW,
				S+E-SE,
				(3*N-NNN)>>1,
				(3*W-WWW)>>1,
				(3*E-EEE)>>1,
				(3*S-SSS)>>1,
				(N+S)>>1,
				(W+E)>>1,
			};
			int vmin, vmax;
			int64_t p1=1LL<<SHIFT>>1;
			int32_t pred;
			for(int k=0;k<NPRED;++k)
				p1+=(int64_t)weights[k]*estim[k];
			p1>>=SHIFT;
			pred=(int)p1;
			vmax=N; vmin=W;
			if(N<W)vmin=N, vmax=W;
			if(vmin>S)vmin=S;
			if(vmax<S)vmax=S;
			if(vmin>E)vmin=E;
			if(vmax<E)vmax=E;
			CLAMP2(pred, vmin, vmax);

			int idx=4*(src->iw*ky+kx)+kc;
			int curr, val=src->data[idx];
			if(fwd)
			{
				pred=-pred;
				curr=val;
			}
			val+=pred;
			if(g_dist<=1)
			{
				val<<=32-src->depth[kc];
				val>>=32-src->depth[kc];
			}
			else if(!fwd)
				CLAMP2(val, amin, amax);
			if(!fwd)
				curr=val;
			src->data[idx]=val;

			int e=(curr>(int32_t)p1)-(curr<(int32_t)p1);
			for(int k=0;k<NPRED;++k)
				weights[k]+=e*estim[k];
		}
	}
#if 0
	if(iw==src->iw)
		messagebox(MBOX_OK, "Info",
			" 0 %8d\n"
			" 1 %8d\n"
			" 2 %8d\n"
			" 3 %8d\n"
			" 4 %8d\n"
			" 5 %8d\n"
			" 6 %8d\n"
			" 7 %8d\n"
			" 8 %8d\n"
			" 9 %8d\n"
			" A %8d\n"
			" B %8d\n"
			, weights[0x0]
			, weights[0x1]
			, weights[0x2]
			, weights[0x3]
			, weights[0x4]
			, weights[0x5]
			, weights[0x6]
			, weights[0x7]
			, weights[0x8]
			, weights[0x9]
			, weights[0xA]
			, weights[0xB]
		);
#endif
}
static void awav_rooks(Image *src, int fwd, int kc, int iw, int ih)
{
	enum
	{
		NPRED=14,
		SHIFT=19,
	};
	int32_t weights[NPRED]={0};
	int amin=-(1<<src->depth[kc]>>1), amax=(1<<src->depth[kc]>>1)-1;

	init_weights(weights, NPRED, SHIFT);
	for(int ky=1;ky<ih;ky+=2)
	{
		for(int kx=0;kx<iw;kx+=2)
		{
			int
				NNN	=LOAD(+0, -3),
				WWW	=LOAD(-3, +0),
				SSS	=LOAD(+0, +3),
				EEE	=LOAD(+3, +0),
				NEE	=LOAD(+2, -1),
				SEE	=LOAD(+2, +1),
				NWW	=LOAD(-2, -1),
				SWW	=LOAD(-2, +1),
				NNE	=LOAD(+1, -2),
				SSE	=LOAD(+1, +2),
				NNW	=LOAD(-1, -2),
				SSW	=LOAD(-1, +2),
				N	=LOAD(+0, -1),
				W	=LOAD(-1, +0),
				S	=LOAD(+0, +1),
				E	=LOAD(+1, +0);
			int estim[]=
			{
				N,
				W,
				S,
				E,
				(4*N-NNE-NNW)>>1,
				(4*S-SSE-SSW)>>1,
				(4*E-NEE-SEE)>>1,
				(4*W-NWW-SWW)>>1,
				(3*N-NNN)>>1,
				(3*W-WWW)>>1,
				(3*E-EEE)>>1,
				(3*S-SSS)>>1,
				(N+S)>>1,
				(W+E)>>1,
			};
			int vmin, vmax;
			int64_t p1=1LL<<SHIFT>>1;
			int32_t pred;
			for(int k=0;k<NPRED;++k)
				p1+=(int64_t)weights[k]*estim[k];
			p1>>=SHIFT;
			pred=(int)p1;
			vmax=N; vmin=W;
			if(N<W)vmin=N, vmax=W;
			if(vmin>S)vmin=S;
			if(vmax<S)vmax=S;
			if(vmin>E)vmin=E;
			if(vmax<E)vmax=E;
			CLAMP2(pred, vmin, vmax);

			int idx=4*(src->iw*ky+kx)+kc;
			int curr, val=src->data[idx];
			if(fwd)
			{
				pred=-pred;
				curr=val;
			}
			val+=pred;
			if(g_dist<=1)
			{
				val<<=32-src->depth[kc];
				val>>=32-src->depth[kc];
			}
			else if(!fwd)
				CLAMP2(val, amin, amax);
			if(!fwd)
				curr=val;
			src->data[idx]=val;

			int e=(curr>(int32_t)p1)-(curr<(int32_t)p1);
			for(int k=0;k<NPRED;++k)
				weights[k]+=e*estim[k];
		}
	}
}
static void awav_bishop(Image *src, int fwd, int kc, int iw, int ih)
{
	enum
	{
		NPRED=14,
		SHIFT=19,
	};
	int32_t weights[NPRED]={0};
	int amin=-(1<<src->depth[kc]>>1), amax=(1<<src->depth[kc]>>1)-1;

	init_weights(weights, NPRED, SHIFT);
	for(int ky=0;ky<ih;ky+=2)
	{
		for(int kx=0;kx<iw;kx+=2)
		{
			int
				NNNWWW	=LOAD(-1, -1),
				NNNEEE	=LOAD(+1, -1),
				SSSWWW	=LOAD(-1, +1),
				SSSEEE	=LOAD(+1, +1),
				NEEE	=LOAD(+3, -1),
				NWWW	=LOAD(-3, -1),
				SEEE	=LOAD(+3, +1),
				SWWW	=LOAD(-3, +1),
				NNNE	=LOAD(+1, -3),
				NNNW	=LOAD(-1, -3),
				SSSE	=LOAD(+1, +3),
				SSSW	=LOAD(-1, +3),
				NW	=LOAD(-1, -1),
				NE	=LOAD(+1, -1),
				SW	=LOAD(-1, +1),
				SE	=LOAD(+1, +1);
			int estim[]=
			{
				NW,
				NE,
				SW,
				SE,
				(3*NW-NNNWWW)>>1,
				(3*NE-NNNEEE)>>1,
				(3*SW-SSSWWW)>>1,
				(3*SE-SSSEEE)>>1,
				(4*NW-NWWW-NNNW)>>1,
				(4*NE-NEEE-NNNE)>>1,
				(4*SW-SWWW-SSSW)>>1,
				(4*SE-SEEE-SSSE)>>1,
				(NW+SE)>>1,
				(SW+NE)>>1,
			};
			int vmin, vmax;
			int64_t p1=1LL<<SHIFT>>1;
			int32_t pred;
			for(int k=0;k<NPRED;++k)
				p1+=(int64_t)weights[k]*estim[k];
			p1>>=SHIFT;
			pred=(int)p1;
			vmax=NW; vmin=NE;
			if(NW<NE)vmin=NW, vmax=NE;
			if(vmin>SW)vmin=SW;
			if(vmax<SW)vmax=SW;
			if(vmin>SE)vmin=SE;
			if(vmax<SE)vmax=SE;
			CLAMP2(pred, vmin, vmax);

			int idx=4*(src->iw*ky+kx)+kc;
			int curr, val=src->data[idx];
			if(fwd)
			{
				pred=-pred;
				curr=val;
			}
			val+=pred;
			if(g_dist<=1)
			{
				val<<=32-src->depth[kc];
				val>>=32-src->depth[kc];
			}
			else if(!fwd)
				CLAMP2(val, amin, amax);
			if(!fwd)
				curr=val;
			src->data[idx]=val;

			int e=(curr>(int32_t)p1)-(curr<(int32_t)p1);
			for(int k=0;k<NPRED;++k)
				weights[k]+=e*estim[k];
		}
	}
}
static void awav_king(Image *src, int fwd, int kc, int iw, int ih)
{
	enum
	{
		NPRED=14,
		SHIFT=15,
		
		XPAD=8,
		NROWS=4,
		NCH=1,
		NVAL=1,
	};
	int32_t weights[NPRED]={0};
	int amin=-(1<<src->depth[kc]>>1), amax=(1<<src->depth[kc]>>1)-1;
	int psize=(iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	int16_t *pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);

	init_weights(weights, NPRED, SHIFT);
	for(int ky=1;ky<ih;ky+=2)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=1;kx<iw;kx+=2)
		{
			rows[0]+=NROWS*NVAL;
			rows[1]+=NROWS*NVAL;
			rows[2]+=NROWS*NVAL;
			rows[3]+=NROWS*NVAL;
			int
				NNN	=LOAD(+0, -3),
				WWW	=LOAD(-3, +0),
				SSS	=LOAD(+0, +3),
				EEE	=LOAD(+3, +0),
				NW	=LOAD(-1, -1),
				NE	=LOAD(+1, -1),
				SW	=LOAD(-1, +1),
				SE	=LOAD(+1, +1),
				N	=LOAD(+0, -1),
				W	=LOAD(-1, +0),
				S	=LOAD(+0, +1),
				E	=LOAD(+1, +0),
				NNWW	=rows[1][-1*NCH*NROWS*NVAL],
				NN	=rows[1][+0*NCH*NROWS*NVAL],
				WW	=rows[0][-1*NCH*NROWS*NVAL];
			int p0=NN+WW-NNWW;
			int estim[]=
			{
				N,
				W,
				S,
				E,
				N+W-NW,
				N+E-NE,
				S+W-SW,
				S+E-SE,
				(3*N-NNN)>>1,
				(3*W-WWW)>>1,
				(3*E-EEE)>>1,
				(3*S-SSS)>>1,
				(N+S)>>1,
				(W+E)>>1,
			};
			int vmin, vmax;
			int64_t p1=1LL<<SHIFT>>1;
			int32_t pred;
			for(int k=0;k<NPRED;++k)
				p1+=(int64_t)weights[k]*estim[k];
			p1>>=SHIFT;
			pred=(int)p1;
			vmax=N; vmin=W;
			if(N<W)vmin=N, vmax=W;
			if(vmin>S)vmin=S;
			if(vmax<S)vmax=S;
			if(vmin>E)vmin=E;
			if(vmax<E)vmax=E;
			CLAMP2(pred, vmin, vmax);
		//	int pred=(3*(N+W+S+E)+(NW+NE+SW+SE))>>4;
			pred=-pred;

			int idx=4*(src->iw*ky+kx)+kc;
			int curr, val=src->data[idx];
			if(fwd)
			{
				pred=-pred;
				curr=val;
			}
			val+=pred;
		//	if(g_dist<=1)
			{
				val<<=32-src->depth[kc];
				val>>=32-src->depth[kc];
			}
		//	else if(!fwd)
		//		CLAMP2(val, amin, amax);
			if(!fwd)
				curr=val;
			src->data[idx]=val;
			rows[0][0]=curr;

			int e=(curr-p0>(int32_t)p1)-(curr-p0<(int32_t)p1);
			for(int k=0;k<NPRED;++k)
				weights[k]-=e*estim[k];
		}
	}
	_mm_free(pixels);
}
static void awav_pawn(Image *src, int fwd, int kc, int iw, int ih)
{
	if(g_dist<=1)
		return;
	int invdist=((1<<16)+g_dist-1)/g_dist;
	//	bishop	queen
	//	rooks	DC
	for(int ky=0;ky<ih;ky+=2)
	{
		for(int kx=0;kx<iw;kx+=2)
		{
			int idx=4*(src->iw*ky+kx)+kc;
			{
				int val=src->data[idx];
				if(fwd)
					val=(val*invdist>>16)-(val>>31);
				else
					val*=g_dist;
				src->data[idx]=val;
			}
			if(kx+1<iw)
			{
				int val=src->data[idx+4];
				if(fwd)
					val=(val*invdist>>16)-(val>>31);
				else
					val*=g_dist;
				src->data[idx+4]=val;
			}
			if(ky+1<ih)
			{
				int val=src->data[idx+4*src->iw];
				if(fwd)
					val=(val*invdist>>16)-(val>>31);
				else
					val*=g_dist;
				src->data[idx+4*src->iw]=val;
			}
		}
	}
}
static void awav_horse(Image *src, int fwd, int kc, int iw, int ih, int *tmp)
{
	if(fwd)
	{
		for(int ky=0;ky<ih;++ky)
		{
			int half=iw>>1;
			for(int kx=0;kx<iw;++kx)
				tmp[kx]=src->data[4*(src->iw*ky+kx)+kc];
			for(int kx=0;kx<half;++kx)
			{
				src->data[4*(src->iw*ky+kx)+kc]=tmp[2*kx+1];
				src->data[4*(src->iw*ky+kx+half)+kc]=tmp[2*kx+0];
			}
			if(iw&1)
				src->data[4*(src->iw*ky+iw-1)+kc]=tmp[iw-1];
		}
		for(int kx=0;kx<iw;++kx)
		{
			int half=ih>>1;
			for(int ky=0;ky<ih;++ky)
				tmp[ky]=src->data[4*(src->iw*ky+kx)+kc];
			for(int ky=0;ky<half;++ky)
			{
				src->data[4*(src->iw*ky+kx)+kc]=tmp[2*ky+1];
				src->data[4*(src->iw*(ky+half)+kx)+kc]=tmp[2*ky+0];
			}
			if(ih&1)
				src->data[4*(src->iw*(ih-1)+kx)+kc]=tmp[ih-1];
		}
	}
	else
	{
		for(int kx=0;kx<iw;++kx)
		{
			int half=ih>>1;
			for(int ky=0;ky<ih;++ky)
				tmp[ky]=src->data[4*(src->iw*ky+kx)+kc];
			for(int ky=0;ky<half;++ky)
			{
				src->data[4*(src->iw*(2*ky+1)+kx)+kc]=tmp[ky];
				src->data[4*(src->iw*(2*ky+0)+kx)+kc]=tmp[ky+half];
			}
			if(ih&1)
				src->data[4*(src->iw*(ih-1)+kx)+kc]=tmp[ih-1];
		}
		for(int ky=0;ky<ih;++ky)
		{
			int half=iw>>1;
			for(int kx=0;kx<iw;++kx)
				tmp[kx]=src->data[4*(src->iw*ky+kx)+kc];
			for(int kx=0;kx<half;++kx)
			{
				src->data[4*(src->iw*ky+2*kx+1)+kc]=tmp[kx];
				src->data[4*(src->iw*ky+2*kx+0)+kc]=tmp[kx+half];
			}
			if(iw&1)
				src->data[4*(src->iw*ky+iw-1)+kc]=tmp[iw-1];
		}
	}
}
#undef  LOAD
void pred_awav(Image *src, int fwd)
{
	enum
	{
		MINDIM=64,
		MAXIT=3,
	};
	int dim=src->iw>src->ih?src->iw:src->ih;
	int *tmp=(int*)malloc(sizeof(int)*dim);
	if(!tmp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
		int it=0;
		if(!src->depth[kc])
			continue;
		if(fwd)
		{
			for(int h2=src->ih, w2=src->iw;h2>MINDIM&&w2>MINDIM&&it<MAXIT;h2>>=1, w2>>=1, ++it)
			{
				awav_queen	(src, fwd, kc, w2, h2);
				awav_rooks	(src, fwd, kc, w2, h2);
				awav_bishop	(src, fwd, kc, w2, h2);
			//	awav_king	(src, fwd, kc, w2, h2);

				awav_pawn	(src, fwd, kc, w2, h2);
				awav_horse	(src, fwd, kc, w2, h2, tmp);
			}
		}
		else
		{
			int sizes[16][2]={0}, nsizes=0;
			for(int h2=src->ih, w2=src->iw;h2>MINDIM&&w2>MINDIM&&it<MAXIT;h2>>=1, w2>>=1, ++it)
			{
				sizes[nsizes][0]=w2;
				sizes[nsizes][1]=h2;
				++nsizes;
			}
			for(int k=nsizes-1;k>=0;--k)
			{
				int w2=sizes[k][0], h2=sizes[k][1];
				
				awav_horse	(src, fwd, kc, w2, h2, tmp);
				awav_pawn	(src, fwd, kc, w2, h2);
				
			//	awav_king	(src, fwd, kc, w2, h2);
				awav_bishop	(src, fwd, kc, w2, h2);
				awav_rooks	(src, fwd, kc, w2, h2);
				awav_queen	(src, fwd, kc, w2, h2);
			}
		}
	}
	free(tmp);
}
