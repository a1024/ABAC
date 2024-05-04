#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


#include"ac.h"
//#include"profiler.h"

static int quantize_ctx(int x)
{
	x=abs(x);
	x=floor_log2_32(x)+1;
	//x=floor_log2_32(x)+1;
	return x;
}
static int CG(int N, int W, int NW)
{
	int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	return pred;
}
static int clamp(int x, int half)
{
	return CLAMP(-half, x, half-1);
}
#define NCTX 2

#define NGCH 5
#define NBCH 6
#define NRCH 9
#define NCH (NRCH+NGCH+NBCH)
#define GPREDLIST\
	PRED(g)\
	PRED(g-((gN+gW)>>1))\
	PRED(g-gCG)\
	PRED(g-((16*gCG+gN+gW+gNW+gNE-2*(gNN+gWW))>>4))\
	PRED(postfilter)
#define BPREDLIST\
	PRED(b)\
	PRED(b-((bN+bW)>>1))\
	PRED(b-CG(bN, bW, bNW))\
	PRED(b-g)\
	PRED(b-clamp(((bN-gN+bW-gW)>>1)+g, half))\
	PRED(b-clamp(CG(bN-gN, bW-gW, bNW-gNW)+g, half))
#define RPREDLIST\
	PRED(r)\
	PRED(r-((rN+rW)>>1))\
	PRED(r-CG(rN, rW, rNW))\
	PRED(r-g)\
	PRED(r-clamp(((rN-gN+rW-gW)>>1)+g, half))\
	PRED(r-clamp(CG(rN-gN, rW-gW, rNW-gNW)+g, half))\
	PRED(r-((g+b)>>1))\
	PRED(r-clamp(((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half))\
	PRED(r-clamp((CG(2*rN-(gN+bN), 2*rW-(gW+bW), 2*rNW-(gNW+bNW))+g+b)>>1, half))

static const char *chnames[]=
{
#define PRED(X) #X,
	GPREDLIST
	BPREDLIST
	RPREDLIST
#undef  PRED
};
static void predict(Image *src, int *predidx, int fwd)
{
	int *pixels=(int*)malloc((src->iw+4LL)*sizeof(int[4*4]));//2 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+4LL)*sizeof(int[4*4]));
	int depth=src->depth, nlevels=1<<depth, half=nlevels>>1;
	int perm[]={1, 2, 0, 3};
	int fwdmask=-fwd;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+4LL)*((ky-0LL)&3)+2)<<2),
			pixels+(((src->iw+4LL)*((ky-1LL)&3)+2)<<2),
			pixels+(((src->iw+4LL)*((ky-2LL)&3)+2)<<2),
			pixels+(((src->iw+4LL)*((ky-3LL)&3)+2)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=src->nch)
		{
			if(!fwd&&src->nch>=3&&predidx[1]==4)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				int luma=src->data[idx+1]-update;//subtract update
				luma+=half;
				luma&=nlevels-1;
				luma-=half;
				src->data[idx+1]=luma;
			}
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pidx=predidx[kc];
				int
					NN	=rows[2][kc+0*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4],
					offset	=0;
				//if(kc0>0)
				if((unsigned)(pidx-8)<3||(unsigned)(pidx-14)<6)
				{
					offset+=rows[0][1];
					if((unsigned)(pidx-18)<3)
						offset=(offset+rows[0][0])>>1;
				}
				int pred=0, vmin, vmax;
				switch(pidx)
				{
				case 0:
				case 5:
				case 8:
				case 11:
				case 14:
				case 17:
					pred=0;
					break;
				case 1:
				case 6:
				case 9:
				case 12:
				case 15:
				case 18:
					pred=(N+W)>>1;
					break;
				case 2:
				case 7:
				case 10:
				case 13:
				case 16:
				case 19:
					pred=N+W-NW;
					vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					pred=CLAMP(vmin, pred, vmax);
					break;
				case 3:
				case 4:
					pred=N+W-NW;
					vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					pred=CLAMP(vmin, pred, vmax);
					pred=(16*pred+N+W+NW+NE-2*(NN+WW))>>4;
					break;
				}
				pred+=offset;
				pred=CLAMP(-half, pred, half-1);

				int curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;

				src->data[idx+kc]=pred;
				rows[0][kc]=(fwd?curr:pred)-offset;
			}
			if(fwd&&src->nch>=3&&predidx[1]==4)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				int luma=src->data[idx+1]+update;//add update
				luma+=half;
				luma&=nlevels-1;
				luma-=half;
				src->data[idx+1]=luma;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}
static int find_min(double *arr, int count)
{
	int idx=0;
	for(int k=1;k<count;++k)
	{
		if(arr[idx]>arr[k])
			idx=k;
	}
	return idx;
}
static void update_CDF(int sym, unsigned *CDF)
{
#if 0
	for(int ks=0;ks<16;++ks)
		CDF[ks]+=(int)(((0x10000-16)&-(ks>sym))+ks-CDF[ks])>>7;
#else
	__m256i ramp0=_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
	__m256i ramp1=_mm256_set_epi32(15, 14, 13, 12, 11, 10, 9, 8);
	__m256i mamp=_mm256_set1_epi32(0x10000-16);

	__m256i msym=_mm256_set1_epi32(sym);
	__m256i mcdf0=_mm256_load_si256((__m256i*)CDF+0);
	__m256i mcdf1=_mm256_load_si256((__m256i*)CDF+1);
	__m256i update0=_mm256_cmpgt_epi32(ramp0, msym);
	__m256i update1=_mm256_cmpgt_epi32(ramp1, msym);
	update0=_mm256_and_si256(update0, mamp);
	update1=_mm256_and_si256(update1, mamp);
	update0=_mm256_add_epi32(update0, ramp0);
	update1=_mm256_add_epi32(update1, ramp1);
	update0=_mm256_sub_epi32(update0, mcdf0);
	update1=_mm256_sub_epi32(update1, mcdf1);
	update0=_mm256_srai_epi32(update0, 7);
	update1=_mm256_srai_epi32(update1, 7);
	mcdf0=_mm256_add_epi32(mcdf0, update0);
	mcdf1=_mm256_add_epi32(mcdf1, update1);
	_mm256_store_si256((__m256i*)CDF+0, mcdf0);
	_mm256_store_si256((__m256i*)CDF+1, mcdf1);
#endif

#if 0
	for(int ks=0;ks<15;++ks)
	{
		if(CDF[ks]>=CDF[ks+1])
			LOG_ERROR("CDF[%d]=0x%08X\nCDF[%d]=0x%08X", ks, CDF[ks], ks+1, CDF[ks+1]);
	}
#endif
}
int f14_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	//PROF_START();
	double t0=time_sec();
	int fwd=src!=0;
	Image tempimage={0}, *errorimage=&tempimage;
	Image const *image=fwd?src:dst;
	if(image->depth&3)
		LOG_ERROR("Bit depth must be divisible by 4, got %d", image->depth);
	int ridx=0, gidx=0, bidx=0;
	int depth=image->depth, nlevels=1<<depth, half=nlevels>>1;
	int rowstride=image->iw*3;
	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 256, 0);
	if(fwd)
	{
		size_t histsize=sizeof(int[NCH])<<image->depth;
		int *hist=(int*)malloc(histsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(hist, 0, histsize);
		for(int ky=2;ky<image->ih;++ky)
		{
			for(int kx=2;kx<image->iw-1;++kx)
			{
#define LOAD(C, X, Y) image->data[idx+rowstride*(Y)+3*(X)+C]
				int
					idx=(image->iw*ky+kx)*3,
					rNW	=LOAD(0, -1, -1),
					rN	=LOAD(0,  0, -1),
					rW	=LOAD(0, -1,  0),
					r	=LOAD(0,  0,  0),
					gNN	=LOAD(1,  0, -2),
					gNW	=LOAD(1, -1, -1),
					gN	=LOAD(1,  0, -1),
					gNE	=LOAD(1,  1, -1),
					gWW	=LOAD(1, -2,  0),
					gW	=LOAD(1, -1,  0),
					g	=LOAD(1,  0,  0),
					bNW	=LOAD(2, -1, -1),
					bN	=LOAD(2,  0, -1),
					bW	=LOAD(2, -1,  0),
					b	=LOAD(2,  0,  0);
#undef  LOAD
				const int postfilter=0;
				int gCG=CG(gN, gW, gNW);
				int preds[]=
				{
#define PRED(X) X,
					GPREDLIST
					BPREDLIST
					RPREDLIST
#undef  PRED
				};
				preds[4]=preds[3]+((preds[10]+preds[16])>>2);

				for(int kt=0;kt<NCH;++kt)
				{
					preds[kt]+=half;
					preds[kt]&=nlevels-1;
					++hist[kt<<depth|preds[kt]];
				}
			}
		}
		double t1=time_sec();
		ptrdiff_t field=(image->iw-3LL)*(image->ih-2LL);
		double esizes[NRCH+NGCH+NBCH]={0};
		for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
		{
			int *curr_hist=hist+((size_t)kt<<depth);
			for(int ks=0;ks<nlevels;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					esizes[kt]-=freq*log2((double)freq/field);
			}
			esizes[kt]/=8;
		}
		double t2=time_sec();
		gidx=find_min(esizes, NGCH-1);
		bidx=find_min(esizes+NGCH, NBCH)+NGCH;
		ridx=find_min(esizes+NGCH+NBCH, NRCH)+NGCH+NBCH;
		double csize0=esizes[ridx]+esizes[gidx]+esizes[bidx];
		double csize1=esizes[16]+esizes[NGCH-1]+esizes[10];
		if(csize1<csize0)
			ridx=16, gidx=NGCH-1, bidx=10;
		if(loud)
		{
			ptrdiff_t usize=field*depth>>3;
			printf("U   %11td\n", usize);
			for(int kt=0;kt<NCH;++kt)
			{
				char c=' ';
				if(kt==ridx)
					c='r';
				if(kt==gidx)
					c='g';
				if(kt==bidx)
					c='b';
				printf("%2d  %14.2lf  %6.2lf%%  %14.6lf  %c  %s\n", kt, esizes[kt], 100.*esizes[kt]/usize, usize/esizes[kt], c, chnames[kt]);
			}
			double esize=esizes[ridx]+esizes[gidx]+esizes[bidx];
			printf("Estimated size %.2lf bytes  CR %lf\n", esize, 3*usize/esize);
			printf("Sampling\t%16.6lf sec  %lf B/s\n", t1-t0, usize/(t1-t0));
			printf("Analysis\t%16.6lf sec\n", t2-t1);
		}
		free(hist);
		int predidx=gidx;			//Y
		predidx=NBCH*predidx+bidx-NGCH;		//Cb
		predidx=NRCH*predidx+ridx-(NGCH+NBCH);	//Cr
		dlist_push_back(&list, &predidx, 2);
		ac_enc_init(&ec, &list);
		
		errorimage=&tempimage;
		image_copy(errorimage, image);
		int predictors[]={ridx, gidx, bidx};
		predict(errorimage, predictors, 1);

		if(loud)
		{
			double t3=time_sec();
			printf("Predict \t%16.6lf sec\n", t3-t2);
		}
	}
	else
	{
		int predidx=0;
		memcpy(&predidx, cbuf, sizeof(short));
		ridx=predidx%NRCH+NGCH+NBCH,	predidx/=NRCH;//Cr
		bidx=predidx%NBCH+NGCH,		predidx/=NBCH;//Cb
		gidx=predidx%NGCH,		predidx/=NGCH;//Y
		if(predidx)
		{
			LOG_ERROR("Corrupt file");
			return 1;
		}
		cbuf+=2;
		clen-=2;
		ac_dec_init(&ec, cbuf, cbuf+clen);

		errorimage=dst;
	}
	int ctxsize=quantize_ctx(nlevels>>1)+1;
	int treesize=0;
	for(int k=0, p=1;k<depth;k+=4, p<<=4)
		treesize+=p;
	size_t statssize=sizeof(int[NCTX*16])*treesize*ctxsize*image->nch;
	unsigned *stats=(unsigned*)_mm_malloc(statssize, sizeof(__m256i));
	size_t bufsize=sizeof(short[4*4*1])*(image->iw+4LL);//4 padded rows * 4 channels max * {errors}
	short *pixels=(short*)malloc(bufsize);
	if(!stats||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	//prepare stats
	for(int ks=0;ks<16;++ks)
		stats[ks]=ks<<12;
	memfill(stats+16, stats, sizeof(int[16])*(treesize-1LL), sizeof(int[16]));
	for(int ks=0;ks<nlevels;++ks)
	{
		int tidx=0;
		for(int kb=depth-4;kb>=0;kb-=4)
		{
			int sym=ks>>kb&15;
			unsigned *CDF=stats+((size_t)tidx<<4);
			update_CDF(sym, CDF);
			tidx=(tidx<<4)+sym+1;
		}
	}
	memfill(stats+((size_t)treesize<<4), stats, statssize-sizeof(int[16])*treesize, sizeof(int[16])*treesize);

	memset(pixels, 0, bufsize);
	int perm[]={1, 2, 0, 3};
	double t3=time_sec();
	//PROF(PREP);
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((image->iw+4LL)*((ky-0LL)&3)+2)<<2),
			pixels+(((image->iw+4LL)*((ky-1LL)&3)+2)<<2),
			pixels+(((image->iw+4LL)*((ky-2LL)&3)+2)<<2),
			pixels+(((image->iw+4LL)*((ky-3LL)&3)+2)<<2),
		};
		for(int kx=0;kx<image->iw;++kx, idx+=3)
		{
			short *N=rows[1];
			short *W=rows[0]-4;
			short *curr=rows[0];

			int ctx1, ctx2, tidx, tidx2, token, sym;
			unsigned *tree1, *CDF1;
			unsigned *tree2, *CDF2;

			for(int kc0=0;kc0<3;++kc0)
			{
				//if(ky==0&&kx==517)//
				//if(ky==2&&kx==649&&kc0==1)//
				//	printf("");

				int kc=perm[kc0];
				ctx1=quantize_ctx(N[kc]);
				ctx2=quantize_ctx(W[kc]);
				tree1=stats+(treesize*(ctxsize*2LL*kc+ctx1)<<4);
				tree2=stats+(treesize*(ctxsize*2LL*kc+ctx2+ctxsize)<<4);
				tidx=0;
				if(fwd)
				{
					token=curr[kc]=errorimage->data[idx+kc];
					token=token<<1^-(token<0);//pack sign for exponential distribution
					for(int kb=depth-4;kb>=0;kb-=4)
					{
						tidx2=tidx<<4;
						CDF1=tree1+tidx2;
						CDF2=tree2+tidx2;
						sym=token>>kb&15;
						int cdf=(CDF1[sym]+CDF2[sym])>>1;
						int freq=(sym>=15?0x10000:((CDF1[sym+1]+CDF2[sym+1])>>1))-cdf;
						ac_enc_update(&ec, cdf, freq);
						update_CDF(sym, CDF1);
						update_CDF(sym, CDF2);
						tidx=(tidx<<4)+sym+1;
					}
				}
				else
				{
					token=0;
					for(int kb=depth-4;kb>=0;kb-=4)
					{
						tidx2=tidx<<4;
						CDF1=tree1+tidx2;
						CDF2=tree2+tidx2;
						int cdf=ac_dec_getcdf(&ec);

						__m256i mlevel=_mm256_set1_epi32(cdf);
						__m256i mcdfA0=_mm256_load_si256((__m256i*)CDF1+0);
						__m256i mcdfA1=_mm256_load_si256((__m256i*)CDF1+1);
						__m256i mcdfB0=_mm256_load_si256((__m256i*)CDF2+0);
						__m256i mcdfB1=_mm256_load_si256((__m256i*)CDF2+1);
						mcdfA0=_mm256_add_epi32(mcdfA0, mcdfB0);
						mcdfA1=_mm256_add_epi32(mcdfA1, mcdfB1);
						mcdfA0=_mm256_srai_epi32(mcdfA0, 1);
						mcdfA1=_mm256_srai_epi32(mcdfA1, 1);
						mcdfA0=_mm256_cmpgt_epi32(mcdfA0, mlevel);
						mcdfA1=_mm256_cmpgt_epi32(mcdfA1, mlevel);
						unsigned short mask=_mm256_movemask_ps(_mm256_castsi256_ps(mcdfA1))<<8|_mm256_movemask_ps(_mm256_castsi256_ps(mcdfA0));
						sym=get_lsb_index16(mask)-1;

						cdf=(CDF1[sym]+CDF2[sym])>>1;
						int freq=(sym>=15?0x10000:((CDF1[sym+1]+CDF2[sym+1])>>1))-cdf;
						ac_dec_update(&ec, cdf, freq);
						token|=sym<<kb;
						update_CDF(sym, CDF1);
						update_CDF(sym, CDF2);
						tidx=(tidx<<4)+sym+1;
					}
					token=token>>1^-(token&1);
					//if((unsigned)(token+half)>=nlevels)//
					//	LOG_ERROR("");

					curr[kc]=dst->data[idx+kc]=token;
				}
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	double t4=time_sec();
	if(fwd)
	{
		image_clear(&tempimage);
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	else
	{
		int predictors[]={ridx, gidx, bidx};
		predict(dst, predictors, 0);
	}
	if(loud)
	{
		double t5=time_sec();
		if(fwd)
		{
			ptrdiff_t usize=((ptrdiff_t)src->iw*src->ih*src->nch*src->depth+7)>>3;
			ptrdiff_t csize=list.nobj;
			printf("Memory usage:      %17.2lf MB\n", (statssize+bufsize+usize)/(1024.*1024));
			printf("  Stats:           %14zd bytes\n", statssize);
			printf("  Circular buffer: %14zd bytes\n", bufsize);
			printf("  Temporary image: %14zd bytes\n", usize);

			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("  Encode\t%16.6lf sec\n", t4-t3);
			printf("  Append\t%16.6lf sec\n", t5-t4);
		}
		else
		{
			printf("  Decode\t%16.6lf sec\n", t4-t3);
			printf("  Unpred\t%16.6lf sec\n", t5-t4);
		}
		printf("%c       \t%16.6lf sec\n", 'D'+fwd, t5-t0);
		if(fwd)
			printf("\n");
		//prof_print();
	}
	if(fwd)
		dlist_clear(&list);
	_mm_free(stats);
	free(pixels);
	return 0;
}