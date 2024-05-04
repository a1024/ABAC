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
static int clamp(int vmin, int x, int vmax)
{
	return CLAMP(vmin, x, vmax);
}
#define NBITS 4		//{1, 2, 4, 8}
#define NCTX 4

#define NGCH 8
#define NBCH 6
#define NRCH 9
#define NCH (NRCH+NGCH+NBCH)
#define GPREDLIST\
	PRED(g)\
	PRED(g-((gN+gW)>>1))\
	PRED(g-gCG)\
	PRED(g-((16*gCG+gN+gW+gNW+gNE-2*(gNN+gWW))>>4))\
	PRED(g-clamp(vmin, (3*(gN+2*gW-gNW)+2*gNE)>>3, vmax))\
	PRED(g-clamp(vmin, (4*(gN+gW)+gW-2*gNW+gNE)>>3, vmax))\
	PRED(g-clamp(vmin, (3*(gN+gW)-2*gNW)>>2, vmax))\
	PRED(g-clamp(vmin, (4*(gN+gW)+gNE-gNW)>>3, vmax))
#define BPREDLIST\
	PRED(b)\
	PRED(b-((bN+bW)>>1))\
	PRED(b-CG(bN, bW, bNW))\
	PRED(b-g)\
	PRED(b-clamp(-half, ((bN-gN+bW-gW)>>1)+g, half-1))\
	PRED(b-clamp(-half, CG(bN-gN, bW-gW, bNW-gNW)+g, half-1))
#define RPREDLIST\
	PRED(r)\
	PRED(r-((rN+rW)>>1))\
	PRED(r-CG(rN, rW, rNW))\
	PRED(r-g)\
	PRED(r-clamp(-half, ((rN-gN+rW-gW)>>1)+g, half-1))\
	PRED(r-clamp(-half, CG(rN-gN, rW-gW, rNW-gNW)+g, half-1))\
	PRED(r-((g+b)>>1))\
	PRED(r-clamp(-half, ((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half-1))\
	PRED(r-clamp(-half, (CG(2*rN-(gN+bN), 2*rW-(gW+bW), 2*rNW-(gNW+bNW))+g+b)>>1, half-1))

static const char *chnames[]=
{
#define PRED(X) #X,
	GPREDLIST
	BPREDLIST
	RPREDLIST
#undef  PRED
};
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
#if NBITS==8 || NBITS==4
ALIGN(32) static const unsigned short gramp[]=
{
	  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
	 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
	 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
	 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
	 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
	 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
	112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
	128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
	144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
	160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
	176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
	192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
	208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
	224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
	240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,
};
#endif
static void update_CDF(int sym, unsigned short *CDF)
{
#if 0
	for(int ks=1;ks<(1<<NBITS);++ks)
		CDF[ks]+=(int)(((0x10000-(1<<NBITS))&-(ks>sym))+ks-CDF[ks])>>7;
#elif NBITS==8
	__m256i mamp=_mm256_set1_epi16((1<<14)-(1<<NBITS)/4);

	__m256i msym=_mm256_set1_epi16(sym);
	for(int k=0;k<16;k+=2)
	{
		__m256i ramp0=_mm256_load_si256((__m256i*)gramp+k+0);
		__m256i ramp1=_mm256_load_si256((__m256i*)gramp+k+1);
		__m256i mcdf0=_mm256_load_si256((__m256i*)CDF+k+0);
		__m256i mcdf1=_mm256_load_si256((__m256i*)CDF+k+1);
		__m256i update0=_mm256_cmpgt_epi16(ramp0, msym);
		__m256i update1=_mm256_cmpgt_epi16(ramp1, msym);
		update0=_mm256_and_si256(update0, mamp);
		update1=_mm256_and_si256(update1, mamp);
		//update0=_mm256_add_epi16(update0, ramp0);
		//update1=_mm256_add_epi16(update1, ramp1);
		update0=_mm256_sub_epi16(update0, mcdf0);
		update1=_mm256_sub_epi16(update1, mcdf1);
		update0=_mm256_srai_epi16(update0, 8);
		update1=_mm256_srai_epi16(update1, 8);
		mcdf0=_mm256_add_epi16(mcdf0, update0);
		mcdf1=_mm256_add_epi16(mcdf1, update1);
		_mm256_store_si256((__m256i*)CDF+k+0, mcdf0);
		_mm256_store_si256((__m256i*)CDF+k+1, mcdf1);
	}
#elif NBITS==4
	__m256i ramp=_mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	__m256i mamp=_mm256_set1_epi16((1<<14)-(1<<NBITS)/4);

	__m256i msym=_mm256_set1_epi16(sym);
	__m256i mcdf=_mm256_load_si256((__m256i*)CDF);
	__m256i update=_mm256_cmpgt_epi16(ramp, msym);
	update=_mm256_and_si256(update, mamp);
	//update=_mm256_add_epi16(update, ramp);
	update=_mm256_sub_epi16(update, mcdf);
	update=_mm256_srai_epi16(update, 7);
	mcdf=_mm256_add_epi16(mcdf, update);
	_mm256_store_si256((__m256i*)CDF, mcdf);
#elif NBITS==2
	CDF[1]+=(int)((((1<<14)-(1<<NBITS)/4)&-(1>sym))-CDF[1])>>8;
	CDF[2]+=(int)((((1<<14)-(1<<NBITS)/4)&-(2>sym))-CDF[1])>>8;
	CDF[3]+=(int)((((1<<14)-(1<<NBITS)/4)&-(3>sym))-CDF[1])>>8;
	//__m128i ramp=_mm_set_epi32(3, 2, 1, 0);
	//__m128i mamp=_mm_set1_epi32(0x10000-(1<<NBITS));
	//
	//__m128i msym=_mm_set1_epi32(sym);
	//__m128i mcdf=_mm_load_si128((__m128i*)CDF);
	//__m128i update=_mm_cmpgt_epi32(ramp, msym);
	//update=_mm_and_si128(update, mamp);
	//update=_mm_add_epi32(update, ramp);
	//update=_mm_sub_epi32(update, mcdf);
	//update=_mm_srai_epi32(update, 8);
	//mcdf=_mm_add_epi32(mcdf, update);
	//_mm_store_si128((__m128i*)CDF, mcdf);
#elif NBITS==1
	CDF[1]+=(int)((((1<<14)-(1<<NBITS)/4)&-(1>sym))-CDF[1])>>8;
#endif

#if 1
	for(int ks=0;ks<(1<<NBITS)-1;++ks)
	{
		if(CDF[ks]>CDF[ks+1])
			LOG_ERROR("CDF[%d]=0x%08X\nCDF[%d]=0x%08X", ks, CDF[ks], ks+1, CDF[ks+1]);
	}
#endif
}
int f16_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	//PROF_START();
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	if(image->depth&(NBITS-1))
	{
		LOG_ERROR("Bit depth must be divisible by %d, got %d", NBITS, image->depth);
		return 2;
	}
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
				int gCG=CG(gN, gW, gNW);
				int vmin=MINVAR(gN, gW);
				int vmax=MAXVAR(gN, gW);
				UPDATE_MIN(vmin, gNE);
				UPDATE_MAX(vmax, gNE);
				int preds[]=
				{
#define PRED(X) X,
					GPREDLIST
					BPREDLIST
					RPREDLIST
#undef  PRED
				};
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
	}
	int ctxsize=quantize_ctx(nlevels>>1)+1;
	int nnodes=0;
	for(int k=0, p=1;k<depth;k+=NBITS, p<<=NBITS)
		nnodes+=p;
	size_t statssize=sizeof(short[NCTX<<NBITS])*nnodes*ctxsize*image->nch;
	unsigned short *stats=(unsigned short*)_mm_malloc(statssize, sizeof(__m256i));
	size_t bufsize=sizeof(short[4*4*2])*(image->iw+4LL);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)malloc(bufsize);
	if(!stats||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	//prepare stats
	for(int ks=0;ks<(1<<NBITS);++ks)
		stats[ks]=ks<<(14-NBITS);
	memfill(stats+(1LL<<NBITS), stats, sizeof(short[1LL<<NBITS])*(nnodes-1LL), sizeof(short[1LL<<NBITS]));
	for(int ks=0;ks<nlevels;++ks)
	{
		int tidx=0;
		for(int kb=depth-NBITS;kb>=0;kb-=NBITS)
		{
			int sym=ks>>kb&((1<<NBITS)-1);
			unsigned short *CDF=stats+((size_t)tidx<<NBITS);
			update_CDF(sym, CDF);
			tidx=(tidx<<NBITS)+sym+1;
		}
	}
	memfill(stats+((size_t)nnodes<<NBITS), stats, statssize-sizeof(short[1LL<<NBITS])*nnodes, sizeof(short[1LL<<NBITS])*nnodes);

	memset(pixels, 0, bufsize);
	int perm[]={1, 2, 0, 3};
	double t3=time_sec();
	//PROF(PREP);
	int predictors[]={ridx, gidx, bidx};
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((image->iw+4LL)*((ky-0LL)&3)+2)<<3),
			pixels+(((image->iw+4LL)*((ky-1LL)&3)+2)<<3),
			pixels+(((image->iw+4LL)*((ky-2LL)&3)+2)<<3),
			pixels+(((image->iw+4LL)*((ky-3LL)&3)+2)<<3),
		};
		for(int kx=0;kx<image->iw;++kx, idx+=3)
		{
			short
				*NN	=rows[2]+0*8,
				*NW	=rows[1]-1*8,
				*N	=rows[1]+0*8,
				*NE	=rows[1]+1*8,
				*WW	=rows[0]-2*8,
				*W	=rows[0]-1*8,
				*curr	=rows[0]+0*8;
			int offset=0;
			int ctx1, ctx2, ctx3, ctx4;
			int tidx, token, sym, cdf, freq;
			unsigned short *tree1, *CDF1;
			unsigned short *tree2, *CDF2;
			unsigned short *tree3, *CDF3;
			unsigned short *tree4, *CDF4;

			for(int kc0=0;kc0<3;++kc0)
			{
				//if(ky==0&&kx==517)//
				//if(ky==2&&kx==649&&kc0==1)//
				//if(ky==1&&kx==1)//
				if(ky==0&&kx==17)//
					printf("");

				int kc=perm[kc0], pidx=predictors[kc];
				if((unsigned)(pidx-(NGCH+3))<3||(unsigned)(pidx-(NGCH+NBCH+3))<6)
				{
					offset+=rows[0][1];
					if((unsigned)(pidx-(NGCH+NBCH+6))<3)
						offset=(offset+rows[0][0])>>1;
				}
				int pred=0, vmin, vmax;
				switch(pidx)
				{
				case 0:
				case NGCH+0:case NGCH+3:
				case NGCH+NBCH+0:case NGCH+NBCH+3:case NGCH+NBCH+6:
					pred=0;
					break;
				case 1:
				case NGCH+1:case NGCH+4:
				case NGCH+NBCH+1:case NGCH+NBCH+4:case NGCH+NBCH+7:
					pred=(N[kc]+W[kc])>>1;
					break;
				case 2:
				case NGCH+2:case NGCH+5:
				case NGCH+NBCH+2:case NGCH+NBCH+5:case NGCH+NBCH+8:
					pred=N[kc]+W[kc]-NW[kc];
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					pred=CLAMP(vmin, pred, vmax);
					break;
				case 3:
					pred=N[kc]+W[kc]-NW[kc];
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					pred=CLAMP(vmin, pred, vmax);
					pred=(16*pred+N[kc]+W[kc]+NW[kc]+NE[kc]-2*(NN[kc]+WW[kc]))>>4;
					break;
				case 4:
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					UPDATE_MIN(vmin, NE[kc]);
					UPDATE_MAX(vmax, NE[kc]);
					pred=clamp(vmin, (3*(N[kc]+2*W[kc]-NW[kc])+2*NE[kc])>>3, vmax);
					break;
				case 5:
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					UPDATE_MIN(vmin, NE[kc]);
					UPDATE_MAX(vmax, NE[kc]);
					pred=clamp(vmin, (4*(N[kc]+W[kc])+W[kc]-2*NW[kc]+NE[kc])>>3, vmax);
					break;
				case 6:
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					UPDATE_MIN(vmin, NE[kc]);
					UPDATE_MAX(vmax, NE[kc]);
					pred=clamp(vmin, (3*(N[kc]+W[kc])-2*NW[kc])>>2, vmax);
					break;
				case 7:
					vmin=MINVAR(N[kc], W[kc]);
					vmax=MAXVAR(N[kc], W[kc]);
					UPDATE_MIN(vmin, NE[kc]);
					UPDATE_MAX(vmax, NE[kc]);
					pred=clamp(vmin, (4*(N[kc]+W[kc])+NE[kc]-NW[kc])>>3, vmax);
					break;
				}
				pred+=offset;
				pred=CLAMP(-half, pred, half-1);

				ctx1=quantize_ctx(N [kc+4]);
				ctx2=quantize_ctx(W [kc+4]);
				ctx3=quantize_ctx(NW[kc+4]);
				ctx4=quantize_ctx(NE[kc+4]);
				tree1=stats+(nnodes*((size_t)ctxsize*NCTX*kc+ctx1+ctxsize*0LL)<<NBITS);
				tree2=stats+(nnodes*((size_t)ctxsize*NCTX*kc+ctx2+ctxsize*1LL)<<NBITS);
				tree3=stats+(nnodes*((size_t)ctxsize*NCTX*kc+ctx3+ctxsize*2LL)<<NBITS);
				tree4=stats+(nnodes*((size_t)ctxsize*NCTX*kc+ctx4+ctxsize*3LL)<<NBITS);
				tidx=0;
				if(fwd)
				{
					curr[kc]=image->data[idx+kc];
					token=curr[kc+4]=curr[kc]-pred;
					token+=half;
					token&=nlevels-1;//error distribution is irrelevant due to movemask/tzcnt
					for(int kb=depth-NBITS;kb>=0;kb-=NBITS)
					{
						CDF1=tree1+((size_t)tidx<<NBITS);
						CDF2=tree2+((size_t)tidx<<NBITS);
						CDF3=tree3+((size_t)tidx<<NBITS);
						CDF4=tree4+((size_t)tidx<<NBITS);
						sym=token>>kb&((1<<NBITS)-1);
						cdf=CDF1[sym]+CDF2[sym]+CDF3[sym]+CDF4[sym]+sym;
						freq=(sym>=(1<<NBITS)-1?0x10000:CDF1[sym+1]+CDF2[sym+1]+CDF3[sym+1]+CDF4[sym+1]+sym+1)-cdf;
						ac_enc_update(&ec, cdf, freq);
						update_CDF(sym, CDF1);
						update_CDF(sym, CDF2);
						update_CDF(sym, CDF3);
						update_CDF(sym, CDF4);
						tidx=(tidx<<NBITS)+sym+1;
					}
				}
				else
				{
					token=0;
					for(int kb=depth-NBITS;kb>=0;kb-=NBITS)
					{
						CDF1=tree1+((size_t)tidx<<NBITS);
						CDF2=tree2+((size_t)tidx<<NBITS);
						CDF3=tree3+((size_t)tidx<<NBITS);
						CDF4=tree4+((size_t)tidx<<NBITS);
						int cdf=ac_dec_getcdf(&ec);

#if NBITS==8
						__m256i mlevel=_mm256_set1_epi16(cdf);
						sym=-2;
						for(int k=0;k<16;k+=2)
						{
							__m256i ramp0=_mm256_load_si256((__m256i*)gramp+k+0);
							__m256i ramp1=_mm256_load_si256((__m256i*)gramp+k+1);
							__m256i mcdfA0=_mm256_load_si256((__m256i*)CDF1+k+0);
							__m256i mcdfA1=_mm256_load_si256((__m256i*)CDF1+k+1);
							__m256i mcdfB0=_mm256_load_si256((__m256i*)CDF2+k+0);
							__m256i mcdfB1=_mm256_load_si256((__m256i*)CDF2+k+1);
							__m256i mcdfC0=_mm256_load_si256((__m256i*)CDF3+k+0);
							__m256i mcdfC1=_mm256_load_si256((__m256i*)CDF3+k+1);
							__m256i mcdfD0=_mm256_load_si256((__m256i*)CDF4+k+0);
							__m256i mcdfD1=_mm256_load_si256((__m256i*)CDF4+k+1);
							mcdfA0=_mm256_add_epi16(mcdfA0, mcdfB0);
							mcdfA1=_mm256_add_epi16(mcdfA1, mcdfB1);
							mcdfA0=_mm256_add_epi16(mcdfA0, mcdfC0);
							mcdfA1=_mm256_add_epi16(mcdfA1, mcdfC1);
							mcdfA0=_mm256_add_epi16(mcdfA0, mcdfD0);
							mcdfA1=_mm256_add_epi16(mcdfA1, mcdfD1);
							//mcdfA0=_mm256_srai_epi32(mcdfA0, 2);
							//mcdfA1=_mm256_srai_epi32(mcdfA1, 2);
							mcdfA0=_mm256_add_epi16(mcdfA0, ramp0);
							mcdfA1=_mm256_add_epi16(mcdfA1, ramp1);
							mcdfA0=_mm256_cmpgt_epi32(mcdfA0, mlevel);
							mcdfA1=_mm256_cmpgt_epi32(mcdfA1, mlevel);
							unsigned long long mask=(unsigned long long)_mm256_movemask_epi8(mcdfA1)<<32|_mm256_movemask_epi8(mcdfA0);
							sym+=get_lsb_index64(mask);
						}
						sym>>=1;
#elif NBITS==4
						__m256i ramp=_mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
						__m256i mlevel=_mm256_set1_epi16(cdf);
						__m256i mcdfA=_mm256_load_si256((__m256i*)CDF1);
						__m256i mcdfB=_mm256_load_si256((__m256i*)CDF2);
						__m256i mcdfC=_mm256_load_si256((__m256i*)CDF3);
						__m256i mcdfD=_mm256_load_si256((__m256i*)CDF4);
						mcdfA=_mm256_add_epi16(mcdfA, mcdfB);
						mcdfA=_mm256_add_epi16(mcdfA, mcdfC);
						mcdfA=_mm256_add_epi16(mcdfA, mcdfD);
						//mcdfA=_mm256_srai_epi16(mcdfA, 2);
						mcdfA=_mm256_add_epi16(mcdfA, ramp);
						mcdfA=_mm256_cmpgt_epi16(mcdfA, mlevel);
						unsigned mask=_mm256_movemask_epi8(mcdfA);
						sym=(get_lsb_index32(mask)>>1)-1;
#elif NBITS==2
						sym =cdf>=CDF1[1]+CDF2[1]+CDF3[1]+CDF4[1]+1;
						sym+=cdf>=CDF1[2]+CDF2[2]+CDF3[2]+CDF4[2]+2;
						sym+=cdf>=CDF1[3]+CDF2[3]+CDF3[3]+CDF4[3]+3;
						//__m128i mlevel=_mm_set1_epi32(cdf);
						//__m128i mcdfA=_mm_load_si128((__m128i*)CDF1);
						//__m128i mcdfB=_mm_load_si128((__m128i*)CDF2);
						//__m128i mcdfC=_mm_load_si128((__m128i*)CDF3);
						//__m128i mcdfD=_mm_load_si128((__m128i*)CDF4);
						//mcdfA=_mm_add_epi32(mcdfA, mcdfB);
						//mcdfA=_mm_add_epi32(mcdfA, mcdfC);
						//mcdfA=_mm_add_epi32(mcdfA, mcdfD);
						//mcdfA=_mm_srai_epi32(mcdfA, 2);
						//mcdfA=_mm_cmpgt_epi32(mcdfA, mlevel);
						//unsigned short mask=_mm_movemask_ps(_mm_castsi128_ps(mcdfA));
						//sym=_tzcnt_u16(mask|0xFFF0)-1;
#elif NBITS==1
						sym=cdf>=CDF1[1]+CDF2[1]+CDF3[1]+CDF4[1];
#endif
						
						cdf=CDF1[sym]+CDF2[sym]+CDF3[sym]+CDF4[sym]+sym;
						freq=(sym>=(1<<NBITS)-1?0x10000:CDF1[sym+1]+CDF2[sym+1]+CDF3[sym+1]+CDF4[sym+1]+sym+1)-cdf;
						ac_dec_update(&ec, cdf, freq);
						token|=sym<<kb;
						update_CDF(sym, CDF1);
						update_CDF(sym, CDF2);
						update_CDF(sym, CDF3);
						update_CDF(sym, CDF4);
						tidx=(tidx<<NBITS)+sym+1;
					}
					token+=pred;
					token&=nlevels-1;
					token-=half;
					dst->data[idx+kc]=curr[kc]=token;
					curr[kc+4]=token-pred;
				}
				curr[kc]-=offset;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	double t4=time_sec();
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
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
			printf("  Decode\t%16.6lf sec\n", t4-t3);
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