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
	x=FLOOR_LOG2_P1(x);
	//x=FLOOR_LOG2_P1(x);
	return x;
}
#define NBITS 4
#define NCTX 4

static void update_CDF(int sym, unsigned short *CDF)
{
#if 0
	for(int ks=1;ks<(1<<NBITS);++ks)
		CDF[ks]+=(int)(((0x10000-(1<<NBITS))&-(ks>sym))+ks-CDF[ks])>>7;
#else
	__m256i ramp=_mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	__m256i mamp=_mm256_set1_epi16((1<<14)-(1<<NBITS)/4);

	__m256i msym=_mm256_set1_epi16(sym);
	__m256i mcdf=_mm256_load_si256((__m256i*)CDF);
	__m256i update=_mm256_cmpgt_epi16(ramp, msym);
	update=_mm256_and_si256(update, mamp);
	update=_mm256_sub_epi16(update, mcdf);
	update=_mm256_srai_epi16(update, 7);
	mcdf=_mm256_add_epi16(mcdf, update);
	_mm256_store_si256((__m256i*)CDF, mcdf);
#endif

#if 0
	for(int ks=0;ks<(1<<NBITS)-1;++ks)
	{
		if(CDF[ks]>CDF[ks+1])
			LOG_ERROR("\nCDF[%d]=0x%08X\nCDF[%d]=0x%08X", ks, CDF[ks], ks+1, CDF[ks+1]);
	}
#endif
}
static void predict(short **rows, short *dst, int half)
{
	static const int perm[]={1, 2, 0, 3};
	short
		*NW	=rows[1]-1*8,
		*N	=rows[1]+0*8,
		*W	=rows[0]-1*8,
		*curr	=rows[0]+0*8;
	for(int kc0=0;kc0<3;++kc0)
	{
		int kc=perm[kc0];
		int offset=0;
		if(kc0>0)
		{
			offset+=curr[1];
			if(kc0>1)
				offset=(2*offset+curr[2])>>1;//(2*g+[b-g])/2 = (g+b)/2
		}
		int pred=0, vmin, vmax, val;
		pred=N[kc]+W[kc]-NW[kc];
		vmin=MINVAR(N[kc], W[kc]);
		vmax=MAXVAR(N[kc], W[kc]);
		pred=CLAMP(vmin, pred, vmax);
		pred+=offset;
		pred=CLAMP(-half, pred, half-1);

		if(dst)//decoding
		{
			val=dst[kc]+pred;
			val+=half;
			val&=(half<<1)-1;
			val-=half;
			dst[kc]=val;
			curr[kc]=val-offset;
		}
		else//encoding
		{
			val=curr[kc]-pred;
			val+=half;
			val&=(half<<1)-1;
			val-=half;
			curr[kc+4]=val;
			curr[kc]-=offset;
		}
	}
}
int f17_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
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
	int depth=image->depth, nlevels=1<<depth, half=nlevels>>1;
	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 65536, 0);
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
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
				*NW	=rows[1]-1*8,
				*N	=rows[1]+0*8,
				*NE	=rows[1]+1*8,
				*W	=rows[0]-1*8,
				*curr	=rows[0]+0*8;
			int ctx1, ctx2, ctx3, ctx4;
			int tidx, token, sym, cdf, freq;
			unsigned short *tree1, *CDF1;
			unsigned short *tree2, *CDF2;
			unsigned short *tree3, *CDF3;
			unsigned short *tree4, *CDF4;

			if(fwd&&image->nch>=3)
			{
				memcpy(curr, image->data+idx, sizeof(short[3]));
				predict(rows, 0, half);
				int luma=curr[4+1]+((curr[4+0]+curr[4+2])>>2);//update
				luma+=half;
				luma&=nlevels-1;
				luma-=half;
				curr[4+1]=luma;
			}
			for(int kc0=0;kc0<3;++kc0)
			{
				int kc=perm[kc0];

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
					token=curr[kc+4];
					token=token<<1^-(token<0);
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

						__m256i ramp=_mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
						__m256i mlevel=_mm256_set1_epi16(cdf);
						__m256i mcdfA=_mm256_load_si256((__m256i*)CDF1);
						__m256i mcdfB=_mm256_load_si256((__m256i*)CDF2);
						__m256i mcdfC=_mm256_load_si256((__m256i*)CDF3);
						__m256i mcdfD=_mm256_load_si256((__m256i*)CDF4);
						__m256i mcdf=_mm256_add_epi16(mcdfA, mcdfB);
						mcdf=_mm256_add_epi16(mcdf, mcdfC);
						mcdf=_mm256_add_epi16(mcdf, mcdfD);
						mcdf=_mm256_add_epi16(mcdf, ramp);
						mcdf=_mm256_subs_epu16(mcdf, mlevel);
						mcdf=_mm256_cmpeq_epi16(mcdf, _mm256_setzero_si256());
						unsigned mask=_mm256_movemask_epi8(mcdf);
						sym=FLOOR_LOG2(mask)>>1;
						
						cdf=CDF1[sym]+CDF2[sym]+CDF3[sym]+CDF4[sym]+sym;
						freq=(sym>=(1<<NBITS)-1?0x10000:CDF1[sym+1]+CDF2[sym+1]+CDF3[sym+1]+CDF4[sym+1]+sym+1)-cdf;
						ac_dec_update(&ec, cdf, freq);
						token|=sym<<kb;

						__m256i mamp=_mm256_set1_epi16((1<<14)-(1<<NBITS)/4);

						__m256i msym=_mm256_set1_epi16(sym);
						__m256i update=_mm256_cmpgt_epi16(ramp, msym);
						update=_mm256_and_si256(update, mamp);
						__m256i updateA=_mm256_sub_epi16(update, mcdfA);
						__m256i updateB=_mm256_sub_epi16(update, mcdfB);
						__m256i updateC=_mm256_sub_epi16(update, mcdfC);
						__m256i updateD=_mm256_sub_epi16(update, mcdfD);
						updateA=_mm256_srai_epi16(updateA, 7);
						updateB=_mm256_srai_epi16(updateB, 7);
						updateC=_mm256_srai_epi16(updateC, 7);
						updateD=_mm256_srai_epi16(updateD, 7);
						mcdfA=_mm256_add_epi16(mcdfA, updateA);
						mcdfB=_mm256_add_epi16(mcdfB, updateB);
						mcdfC=_mm256_add_epi16(mcdfC, updateC);
						mcdfD=_mm256_add_epi16(mcdfD, updateD);
						_mm256_store_si256((__m256i*)CDF1, mcdfA);
						_mm256_store_si256((__m256i*)CDF2, mcdfB);
						_mm256_store_si256((__m256i*)CDF3, mcdfC);
						_mm256_store_si256((__m256i*)CDF4, mcdfD);

#if 0
						for(int ks=0;ks<(1<<NBITS)-1;++ks)
						{
							if(CDF[ks]>CDF[ks+1])
								LOG_ERROR("\nCDF[%d]=0x%08X\nCDF[%d]=0x%08X", ks, CDF[ks], ks+1, CDF[ks+1]);
						}
#endif
						//update_CDF(sym, CDF1);
						//update_CDF(sym, CDF2);
						//update_CDF(sym, CDF3);
						//update_CDF(sym, CDF4);
						tidx=(tidx<<NBITS)+sym+1;
					}
					token=token>>1^-(token&1);
					curr[kc+4]=token;
				}
			}
			if(!fwd&&image->nch>=3)
			{
				memcpy(dst->data+idx, curr+4, sizeof(short[3]));
				int luma=dst->data[idx+1]-((dst->data[idx+0]+dst->data[idx+2])>>2);//unupdate
				luma+=half;
				luma&=nlevels-1;
				luma-=half;
				dst->data[idx+1]=luma;
				predict(rows, dst->data+idx, half);
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