#include"e2.h"
//#define EC_USE_ARRAY
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


	#define CLASSIFY_IMAGE
//	#define CHECK_OOB
//	#define PROFILER 1		//0: use __rdtsc()	1: use time_sec()
//	#define ESTIMATE_CSIZES
//	#define ENABLE_GUIDE



#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(GETCTX)\
	CHECKPOINT(QUANTIZE)\
	CHECKPOINT(ENTROPYCODER)\
	CHECKPOINT(RESCALE)\
	CHECKPOINT(UPDATE_CDF)\
	CHECKPOINT(TOARRAY)
#endif
#include"profiler.h"


//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 5
#define CONFIG_MSB 2
#define CONFIG_LSB 0
	#define CLEVELS 9
//	#define CLEVELS 33
#define NCTX (CLEVELS*CLEVELS)
static void update_CDF(const int *hist, unsigned *CDF, int qlevels)
{
	int sum=hist[qlevels], c=0;
	for(int ks=0;ks<qlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<16)-qlevels)/sum)+ks;
		c+=freq;
	}
	CDF[qlevels]=1<<16;
}
int quantize_ctx(int val)
{
	int negmask=-(val<0);
	val=abs(val);
	val=floor_log2_32(val)+1;//[0~0x10000] -> [0~0x10]
	val=floor_log2_32(val)+1;//[0~0x10] -> [0~4]
	val^=negmask;
	val-=negmask;
	val+=CLEVELS>>1;
#ifdef CHECK_OOB
	if((unsigned)val>=CLEVELS)
		LOG_ERROR("Context OOB");
#endif
	return val;
}
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#ifdef CLASSIFY_IMAGE
typedef enum PredTypeEnum
{
	PRED_ZERO,
	//PRED_W,
	//PRED_W_NW,
	//PRED_NW,
	//PRED_N_NW,
	//PRED_N,
	//PRED_N_NE,
	//PRED_NE,
	//PRED_N_W_2,
	//PRED_W_NE_N,
	//PRED_Nbias,
	//PRED_Wbias,
	//PRED_AV4,
	//PRED_N_W,
	PRED_N_W_adj,
	PRED_CGRAD,

	PRED_COUNT,
} PredType;
#endif
static const int compidx[]=
{
	0, 1, 2, 3,
	0, 1, 2, 3,
	1, 2, 0, 3,
	2, 0, 1, 3,
};
static void rct_fwd(int *comp, int krct)
{
	if(krct)
	{
		const int *idx=compidx+((size_t)krct<<2);
		int Y=comp[idx[0]], U=comp[idx[1]], V=comp[idx[2]];
		U-=Y;
		V-=Y;
		Y+=(U+V)>>2;
		comp[0]=Y;
		comp[1]=U;
		comp[2]=V;
	}
}
static void rct_inv(int *comp, int krct)
{
	if(krct)
	{
		const int *idx=compidx+((size_t)krct<<2);
		int Y=comp[0], U=comp[1], V=comp[2];
		Y-=(U+V)>>2;
		V+=Y;
		U+=Y;
		comp[idx[0]]=Y;
		comp[idx[1]]=U;
		comp[idx[2]]=V;
	}
}
static int predict(int N, int W, int NW, int NE, PredType kp)
{
	int pred=0;
	switch(kp)
	{
	//case PRED_W:
	//	pred=W;
	//	break;
	//case PRED_W_NW:
	//	pred=(W+NW)>>1;
	//	break;
	//case PRED_NW:
	//	pred=NW;
	//	break;
	//case PRED_N_NW:
	//	pred=(N+NW)>>1;
	//	break;
	//case PRED_N:
	//	pred=N;
	//	break;
	//case PRED_N_NE:
	//	pred=(N+NE)>>1;
	//	break;
	//case PRED_NE:
	//	pred=NE;
	//	break;
	//case PRED_N_W_2:
	//	pred=(N+W)>>2;
	//	break;
	//case PRED_W_NE_N:
	//	pred=W+NE-N;
	//	break;
	//case PRED_Nbias:
	//	pred=(3*N+W)>>2;
	//	break;
	//case PRED_Wbias:
	//	pred=(N+3*W)>>2;
	//	break;
	//case PRED_AV4:
	//	pred=(N+W+NW+NE)>>2;
	//	break;
	//case PRED_N_W:
	//	pred=(N+W)>>1;
	//	break;
	case PRED_N_W_adj:
		pred=(4*(N+W)+NE-NW)>>3;
		break;
	case PRED_CGRAD:
		pred=N+W-NW;
		pred=MEDIAN3(N, W, pred);
		break;
	}
	return pred;
}
int t54_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	PROF_START();
	double t_start=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;

	//Image *debugbuf=0;//
	//if(fwd)//
	//	image_copy(&debugbuf, image);//

#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	int nch=(image->depth[0]!=0)+(image->depth[1]!=0)+(image->depth[2]!=0)+(image->depth[3]!=0);
	UPDATE_MIN(nch, image->nch);
	if(loud)
	{
		int maxdepth=calc_maxdepth(image, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T54  %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, image->iw, image->ih, maxdepth);
	}
	int use_rct=1, use_pred=1;
	int *hist=(int*)malloc(sizeof(int[NCTX*82*4]));
	unsigned *CDF=(unsigned*)malloc(sizeof(int[NCTX*82*4]));
	int *pixels=(int*)malloc((image->iw+4LL)*sizeof(int[4*16]));//padded 4 rows * ({4 pixels, 4 errors} OR {RCT0, RCT1, RCT2, RCT3})
	if(!hist||!CDF||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	char depths[4];
	memcpy(depths, image->depth, sizeof(char[4]));
	if(nch>=3)
	{
		++depths[0];
		++depths[1];
		++depths[2];
		//char temp;
		//ROTATE3(depths[0], depths[1], depths[2], temp);
	}
	int nlevels[4], qlevels[4];
	for(int kc=0;kc<4;++kc)
	{
		nlevels[kc]=1<<depths[kc];
		switch(depths[kc])
		{
		case 0:		qlevels[kc]=1;break;
		case 8:		qlevels[kc]=45;break;
		case 9:		qlevels[kc]=49;break;
		case 16:	qlevels[kc]=77;break;
		case 17:	qlevels[kc]=81;break;
		default:
			LOG_ERROR("Unsupported bit depth %d", depths[kc]);
			return 0;
		}
	}
#ifdef CLASSIFY_IMAGE
	PredType config=PRED_ZERO;
	if(fwd)
	{
		int maxdepth=depths[0];
		UPDATE_MAX(maxdepth, depths[1]);
		UPDATE_MAX(maxdepth, depths[2]);
		UPDATE_MAX(maxdepth, depths[3]);
		int maxlevels=1<<maxdepth;
		int *hist2=(int*)malloc(sizeof(int[4*PRED_COUNT*4])<<maxdepth);//4 channels * 10 predictors * 4 RCTs
		if(!hist2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(hist2, 0, sizeof(int[4*PRED_COUNT*4])<<maxdepth);
		memset(pixels, 0, (image->iw+4LL)*sizeof(int[4*16]));
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int kym[]=
			{
				(image->iw+4)*((ky-0)&3),
				(image->iw+4)*((ky-1)&3),
				(image->iw+4)*((ky-2)&3),
				(image->iw+4)*((ky-3)&3),
			};
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				int *comp=pixels+((kym[0]+2LL+kx+0)<<4);
				memfill(comp, image->data+((size_t)idx<<2), sizeof(int[4*4]), sizeof(int[4]));
				for(int krct=0;krct<=(3&-(nch>=3));++krct)
				{
					rct_fwd(comp, krct);
					for(int kc=0;kc<nch;++kc)
					{
						int
							NW	=pixels[(kym[1]+kx+2-1)<<4|krct<<2|kc],
							N	=pixels[(kym[1]+kx+2+0)<<4|krct<<2|kc],
							NE	=pixels[(kym[1]+kx+2+1)<<4|krct<<2|kc],
							W	=pixels[(kym[0]+kx+2-1)<<4|krct<<2|kc];
						int v0=comp[kc];
						for(int kp=0;kp<PRED_COUNT;++kp)
						{
							int pred=predict(N, W, NW, NE, (PredType)kp);
							int val=v0-pred;
							val+=maxlevels>>1;
							val&=maxlevels-1;
							++hist2[((kp<<2|krct)<<2|kc)<<maxdepth|val];
						}
					}
					comp+=4;
				}
			}
		}
		double entropy[4*PRED_COUNT*4]={0};
		int res=image->iw*image->ih;
		double bestsize=0;
		for(int krct=0;krct<4;++krct)
		{
			for(int kp=0;kp<PRED_COUNT;++kp)
			{
				int idx=(kp<<2|krct)<<2;
				//if(loud)
				//	printf("%-8s %-10s RGBAT:", RCTnames[krct], prednames[kp]);
				for(int kc=0;kc<4;++kc)
				{
					for(int ks=0;ks<maxlevels;++ks)
					{
						int freq=hist2[(idx|kc)<<maxdepth|ks];
						if(freq)
							entropy[idx|kc]-=freq*log2((double)freq/res);
					}
					entropy[idx|kc]/=8;
					//if(loud)
					//	printf("%16.3lf ", entropy[idx|kc]);
				}
				double csize=
					entropy[idx|0]+
					entropy[idx|1]+
					entropy[idx|2]+
					entropy[idx|3];
				//if(loud)
				//	printf("%16.3lf\n", csize);
				if(!idx||bestsize>csize)
					config=idx>>2, bestsize=csize;
			}
		}
		//if(loud)
		//	printf("Selected: %-8s %-10s %16.3lf  %lf sec\n", RCTnames[config&3], prednames[config>>2], bestsize, time_sec()-t_start);
		if(loud)
		{
			const char *RCTnames[]=
			{
				" R G B ",
				"[R]G B ",
				" R[G]B ",
				" R G[B]",
			};
			const char *prednames[]=
			{
				"zero",
				//"W",
				//"(W+NW)/2",
				//"NW",
				//"(N+NW)/2",
				//"N",
				//"(N+NE)/2",
				//"NE",
				//"(N+W)/4",
				//"W+NE-N",
				//"(3*N+W)/4",
				//"(N+3*W)/4",
				//"(N+W+NW+NE)/4",
				//"(N+W)/2",
				"(4*(N+W)+NE-NW)/8",
				"cgrad",
			};
			for(int krct=0;krct<4;++krct)
			{
				for(int kp=0;kp<PRED_COUNT;++kp)
				{
					int idx=(kp<<2|krct)<<2;
					double csize=
						entropy[idx|0]+
						entropy[idx|1]+
						entropy[idx|2]+
						entropy[idx|3];
					printf("%-8s %-20s RGBAT:", RCTnames[krct], prednames[kp]);
					for(int kc=0;kc<4;++kc)
						printf("%16.3lf ", entropy[idx|kc]);
					printf("%16.3lf", csize);
					if(config==(idx>>2))
						printf(" <- %lf sec", time_sec()-t_start);
					printf("\n");
				}
			}
		}
		free(hist2);
	}
#endif
	memset(pixels, 0, (image->iw+4LL)*sizeof(int[4*16]));
#ifdef ESTIMATE_CSIZES
	double csizes[2*4]={0};
#endif
#ifndef EC_USE_ARRAY
	DList list;
	dlist_init(&list, 1, 256, 0);
#endif
#ifdef CLASSIFY_IMAGE
	if(fwd)
#ifdef EC_USE_ARRAY
		array_append(data, &config, 1, 1, 1, 0, 0);
#else
		dlist_push_back1(&list, &config);
#endif
	else
	{
		config=*cbuf;
		++cbuf;
		--clen;
	}
#endif
	ArithmeticCoder ec;
	if(fwd)
	{
#ifdef EC_USE_ARRAY
		ac_enc_init(&ec, data);
#else
		ac_enc_init(&ec, &list);
#endif
	}
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	for(int kc=0;kc<4;++kc)
	{
		int *curr_hist=hist+NCTX*82*kc;
		unsigned *curr_CDF=CDF+NCTX*82*kc;
		memset(curr_hist, 0, sizeof(int[82]));
		memset(curr_CDF, 0, sizeof(int[82]));
		int sum=0;
		for(int k=0;k<qlevels[kc];++k)
		{
			int val=qlevels[kc]-k;
			val*=val;
			sum+=curr_hist[k]=val*val;
		}
		curr_hist[qlevels[kc]]=sum;
		update_CDF(curr_hist, curr_CDF, qlevels[kc]);
		memfill(curr_hist+82, curr_hist, sizeof(int[82*(NCTX-1)]), sizeof(int[82]));
		memfill(curr_CDF+82, curr_CDF, sizeof(int[82*(NCTX-1)]), sizeof(int[82]));
	}
	PROF(INIT);
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		int kym[]=
		{
			(image->iw+4)*((ky-0)&3),
			(image->iw+4)*((ky-1)&3),
			(image->iw+4)*((ky-2)&3),
			(image->iw+4)*((ky-3)&3),
		};
		for(int kx=0;kx<image->iw;++kx, ++idx)
		{
			int *comp=pixels+(((size_t)kym[0]+kx+2+0)<<3);
			if(fwd)
			{
				int *comp=pixels+(((size_t)kym[0]+kx+2+0)<<3);
				memcpy(comp, image->data+((size_t)idx<<2), sizeof(int[4]));
#ifdef CLASSIFY_IMAGE
				rct_fwd(comp, config&3);
#else
				if(nch>=3)
				{
					comp[0]-=comp[1];
					comp[2]-=comp[1];
					comp[1]+=(comp[0]+comp[2])>>2;
					int temp;
					ROTATE3(comp[0], comp[1], comp[2], temp);
				}
#endif
			}
			for(int kc=0;kc<nch;++kc)
			{
				int nlevels_kc=nlevels[kc], qlevels_kc=qlevels[kc];
				int
					NW	=pixels[(kym[1]+kx+2-1)<<3|0|kc],
					N	=pixels[(kym[1]+kx+2+0)<<3|0|kc],
					NE	=pixels[(kym[1]+kx+2+1)<<3|0|kc],
					W	=pixels[(kym[0]+kx+2-1)<<3|0|kc],
					eN	=pixels[(kym[1]+kx+2+0)<<3|4|kc],
					eW	=pixels[(kym[0]+kx+2-1)<<3|4|kc];
#ifndef CLASSIFY_IMAGE
				int cgrad=N+W-NW-((eW+eN)>>7);
				cgrad=MEDIAN3(N, W, cgrad);
#endif
				int ctx=CLEVELS*quantize_ctx(eN)+quantize_ctx(eW);
				ctx=82*(NCTX*kc+ctx);
				int *curr_hist=hist+ctx;
				unsigned *curr_CDF=CDF+ctx;
				PROF(GETCTX);
				int token=0, delta=0, curr=0;
				if(fwd)
				{
					curr=pixels[(kym[0]+kx+2+0)<<3|0|kc];
					delta=curr;
#ifdef CLASSIFY_IMAGE
					delta-=predict(N, W, NW, NE, (PredType)(config>>2));
#else
					delta-=cgrad;
#endif
					delta+=nlevels_kc>>1;
					delta&=nlevels_kc-1;
					delta-=nlevels_kc>>1;

					int val=delta<<1^-(delta<0);

					int bypass=0, nbits=0;
					if(val<(1<<CONFIG_EXP))
					{
						token=val;//token
						nbits=0;
						bypass=0;
					}
					else
					{
						int lgv=floor_log2_32((unsigned)val);
						int mantissa=val-(1<<lgv);
						token = (1<<CONFIG_EXP) + (
								(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
								(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
								(mantissa&((1<<CONFIG_LSB)-1))
							);
						nbits=lgv-CONFIG_MSB+CONFIG_LSB;
						bypass=val>>CONFIG_LSB&((1LL<<nbits)-1);
					}
					PROF(QUANTIZE);
					
#ifdef ESTIMATE_CSIZES
					if(loud&&fwd)
					{
						double p=(double)(curr_CDF[token+1]-curr_CDF[token])/0x10000;
						csizes[kc<<1|0]-=log2(p);
						csizes[kc<<1|1]+=nbits;
					}
#endif

					ac_enc(&ec, token, curr_CDF, qlevels_kc, 0);
					if(nbits)
					{
						while(nbits>8)
						{
							ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
							nbits-=8;
						}
						ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
					}
				}
				else
				{
					token=ac_dec(&ec, curr_CDF, qlevels_kc, 0);
					delta=token;
					if(delta>=(1<<CONFIG_EXP))
					{
						delta-=1<<CONFIG_EXP;
						int lsb=delta&((1<<CONFIG_LSB)-1);
						delta>>=CONFIG_LSB;
						int msb=delta&((1<<CONFIG_MSB)-1);
						delta>>=CONFIG_MSB;
						int nbits=delta+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB), n=nbits;
						int bypass=0;
						while(n>8)
						{
							n-=8;
							bypass|=ac_dec(&ec, 0, 1<<8, 16-8)<<n;
						}
						bypass|=ac_dec(&ec, 0, 1<<n, 16-n);
						delta=1;
						delta<<=CONFIG_MSB;
						delta|=msb;
						delta<<=nbits;
						delta|=bypass;
						delta<<=CONFIG_LSB;
						delta|=lsb;
					}
					delta=delta>>1^-(delta&1);
					curr=delta;
#ifdef CLASSIFY_IMAGE
					curr+=predict(N, W, NW, NE, (PredType)(config>>2));
#else
					curr+=cgrad;
#endif
					curr+=nlevels_kc>>1;
					curr&=nlevels_kc-1;
					curr-=nlevels_kc>>1;
					pixels[(kym[0]+kx+2+0)<<3|0|kc]=curr;
				}
				PROF(ENTROPYCODER);
				
				//if(fwd)//
				//	debugbuf->data[idx<<2|kc]=delta;
				comp[4|kc]=delta;
				if(curr_hist[qlevels_kc]+1>0x1000)
				{
					int sum=0;
					for(int k=0;k<qlevels_kc;++k)
						sum+=curr_hist[k]=(curr_hist[k]+1)>>1;
					curr_hist[qlevels_kc]=sum;
				}
				++curr_hist[token];
				++curr_hist[qlevels_kc];
				PROF(RESCALE);

				if(!(idx&63))
					update_CDF(curr_hist, curr_CDF, qlevels_kc);
				PROF(UPDATE_CDF);
			}
			if(!fwd)
			{
				int *comp=pixels+(((size_t)kym[0]+kx+2+0)<<3);
				int *comp2=dst->data+((size_t)idx<<2);
				memcpy(comp2, comp, sizeof(int[4]));
#ifdef CLASSIFY_IMAGE
				rct_inv(comp2, config&3);
#else
				if(nch>=3)
				{
					int temp;
					ROTATE3(comp2[2], comp2[1], comp2[0], temp);
					comp2[1]-=(comp2[0]+comp2[2])>>2;
					comp2[2]+=comp2[1];
					comp2[0]+=comp2[1];
				}
#endif
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(comp2, guide->data+((size_t)idx<<2), sizeof(int[4])))
				{
					int comp0[4];
					memcpy(comp0, guide->data+((size_t)idx<<2), sizeof(int[4]));
					comp2[0]-=comp2[1];
					comp2[2]-=comp2[1];
					comp2[1]+=(comp2[0]+comp2[2])>>2;
					comp0[0]-=comp0[1];
					comp0[2]-=comp0[1];
					comp0[1]+=(comp0[0]+comp0[2])>>2;
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
#ifndef EC_USE_ARRAY
		dlist_appendtoarray(&list, data);
#endif
	}
	if(loud)
	{
		t_start=time_sec()-t_start;
		if(fwd)
		{
#ifdef ESTIMATE_CSIZES
			for(int kc=0;kc<nch;++kc)
				printf("C%d  %12.2lf  %12.2lf\n", kc, csizes[kc<<1|0]/8, csizes[kc<<1|1]/8);
#endif
			printf("csize %12lld  %7.3lf%%\n",
#ifdef EC_USE_ARRAY
				data[0]->count,
#else
				list.nobj,
#endif
				100.*list.nobj/image_getBMPsize(image)
			);
		}
		printf("%c %12.3lf sec\n", 'D'+fwd, t_start);
		prof_print();
	}
	//if(fwd)//
	//{
	//	image_snapshot(debugbuf);
	//	free(debugbuf);
	//}

#ifndef EC_USE_ARRAY
	dlist_clear(&list);
#endif
	free(hist);
	free(CDF);
	free(pixels);
	return 1;
}