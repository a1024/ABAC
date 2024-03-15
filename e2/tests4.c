#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

//auxiliary functions
double calc_bitsize(unsigned *CDF, int nlevels, int sym)
{
	int freq=CDF?CDF[sym+1]-CDF[sym]:0x10000/nlevels;
	//if(!freq)
	//	LOG_ERROR("ZPS");
	double prob=(double)freq/0x10000;
	double bitsize=-log2(prob);

	return bitsize;
}
#if 0
double calc_csize_from_hist(int *hist, int nlevels, double *ret_usize)//works only with unit-increment histograms initialized with ones
{
	double csize=0;
	int sum=0;
	for(int ks=0;ks<nlevels;++ks)
		sum+=hist[ks];
	if(!sum)
	{
		if(ret_usize)
			*ret_usize=0;
		return 0;
	}
	for(int ks=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		if(freq)
		{
			double p=(double)freq/sum;
			double bitsize=-freq*log2(p);
			csize+=bitsize;
		}
	}
	if(ret_usize)
		*ret_usize=(sum-nlevels)*log2(nlevels)/8;
	return csize/8;
}
#endif


//T45 CALIC


#define CALIC_SSE_BITS (2+8)
#define CALIC_REACH 2
typedef enum NBIDEnum
{
	NB_NNWW,
	NB_NNW,
	NB_NN,
	NB_NNE,
	NB_NNEE,

	NB_NWW,
	NB_NW,
	NB_N,
	NB_NE,
	NB_NEE,

	NB_WW,
	NB_W,

	NB_eN,
	NB_eW,

	NB_COUNT,
} NBID;
typedef struct EntropyTablesStruct
{
	unsigned
		ctr_ct0[ 18],//continuous-tone mode histograms
		ctr_ct1[ 26],
		ctr_ct2[ 34],
		ctr_ct3[ 50],
		ctr_ct4[ 66],
		ctr_ct5[ 82],
		ctr_ct6[114],
		ctr_ct7[256],

		ctr_bin[32][3];//binary mode histograms
} EntropyTables;
typedef struct CalicStateStruct
{
	Image const *image;
	Image *errors;
	int nb[NB_COUNT], nb2[6], unb[6];
	int nunique, dh, dv, d45, d135, energy, delta, beta, sse_correction, pred, error, e2;
	long long currrow_sum_absE;
	int prevrow_av_absE, lg_lambda;
	int half, clamp_lo, clamp_hi;
	unsigned CDF[257];
	int *sse;
	EntropyTables *tables;
	unsigned *ct_tables[8], ct_tsizes[8], nbypass[8], ct_inc[8], bin_inc[32];
	ArithmeticCoder ec;

	int loud, kc;
	double csizes[4];
} CalicState;
static void calic_init(CalicState *state, Image const *image)
{
	int res=image->iw*image->ih;
	memset(state, 0, sizeof(*state));
	state->image=image;

	image_copy_nodata(&state->errors, image);
	state->sse=(int*)malloc((1LL<<CALIC_SSE_BITS)*sizeof(int));
	state->tables=(EntropyTables*)malloc(sizeof(EntropyTables));
	if(!state->errors||!state->sse||!state->tables)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(state->errors->data, 0, (size_t)res*sizeof(int[4]));

	state->ct_tables[0]=state->tables->ctr_ct0, state->ct_tsizes[0]=_countof(state->tables->ctr_ct0);
	state->ct_tables[1]=state->tables->ctr_ct1, state->ct_tsizes[1]=_countof(state->tables->ctr_ct1);
	state->ct_tables[2]=state->tables->ctr_ct2, state->ct_tsizes[2]=_countof(state->tables->ctr_ct2);
	state->ct_tables[3]=state->tables->ctr_ct3, state->ct_tsizes[3]=_countof(state->tables->ctr_ct3);
	state->ct_tables[4]=state->tables->ctr_ct4, state->ct_tsizes[4]=_countof(state->tables->ctr_ct4);
	state->ct_tables[5]=state->tables->ctr_ct5, state->ct_tsizes[5]=_countof(state->tables->ctr_ct5);
	state->ct_tables[6]=state->tables->ctr_ct6, state->ct_tsizes[6]=_countof(state->tables->ctr_ct6);
	state->ct_tables[7]=state->tables->ctr_ct7, state->ct_tsizes[7]=_countof(state->tables->ctr_ct7);
}
static void calic_free(CalicState *state)
{
	free(state->errors);
	free(state->sse);
	free(state->tables);
}
static void calic_nextchannel(CalicState *state, int kc)
{
	int fillval=1;
	memfill(state->tables, &fillval, sizeof(*state->tables), sizeof(int));
	fillval=0x2000;
	memfill(state->ct_inc, &fillval, sizeof(state->ct_inc), sizeof(int));
	memfill(state->bin_inc, &fillval, sizeof(state->bin_inc), sizeof(int));
	memset(state->sse, 0, (1LL<<CALIC_SSE_BITS)*sizeof(int));
	state->kc=kc;
	int nlevels=1<<state->image->depth[kc];
	state->half=nlevels>>1;
	state->clamp_lo=-state->half;
	state->clamp_hi=state->half-1;
}
static void calic_nextrow(CalicState *state)
{
	state->prevrow_av_absE=(int)((state->currrow_sum_absE+(state->image->iw>>1))/state->image->iw);

	state->lg_lambda=ceil_log2(state->prevrow_av_absE)-5;//lambda = 2^(-(z-8)/2-max(0, ceil_log2(sigma)-5))
	if(state->lg_lambda<0)
		state->lg_lambda=0;
	state->lg_lambda-=(state->image->depth[state->kc]-8)>>1;

	state->currrow_sum_absE=0;

	for(int kt=0;kt<_countof(state->ct_tsizes);++kt)
	{
		int nlevels=state->ct_tsizes[kt];
		if(state->prevrow_av_absE>(nlevels>>1))
			state->nbypass[kt]=floor_log2(state->prevrow_av_absE/(nlevels>>1))+1;
		else
			state->nbypass[kt]=0;
	}
}
static void calic_prepctx(CalicState *state, int kc, int kx, int ky)
{
	int idx=0;
	const int *pixels=state->image->data, *errors=state->errors->data;
	for(int ky2=-CALIC_REACH;ky2<=0;++ky2)
	{
		for(int kx2=-CALIC_REACH;kx2<=CALIC_REACH;++kx2, ++idx)
		{
			if(!ky2&&!kx2)
				break;
			if((unsigned)(ky+ky2)<(unsigned)state->image->ih&&(unsigned)(kx+kx2)<(unsigned)state->image->iw)
				state->nb[idx]=pixels[(state->image->iw*(ky+ky2)+kx+kx2)<<2|kc];
			else
				state->nb[idx]=0;
		}
	}
	state->nb[idx++]=ky?errors[(state->image->iw*(ky-1)+kx)<<2|kc]:0;
	state->nb[idx++]=kx?errors[(state->image->iw*ky+kx-1)<<2|kc]:0;

	state->nb2[0]=state->nb[NB_W];
	state->nb2[1]=state->nb[NB_N];
	state->nb2[2]=state->nb[NB_NW];
	state->nb2[3]=state->nb[NB_NE];
	state->nb2[4]=state->nb[NB_WW];
	state->nb2[5]=state->nb[NB_NN];

	state->unb[0]=state->nb2[0];
	state->nunique=1;
	for(int k=1;k<_countof(state->nb2);++k)
	{
		int val=state->nb2[k];
		int found=0;
		for(int k2=0;k2<state->nunique;++k2)
		{
			if(val==state->unb[k2])
			{
				found=1;
				break;
			}
		}
		if(!found)
			state->unb[state->nunique++]=val;
	}
}
static const int calic_qlevels[]={5, 15, 25, 42, 60, 85, 140};
static int calic_ct(CalicState *state, int curr, int enc)//continuous-tone mode
{
	state->dh=abs(state->nb[NB_W]-state->nb[NB_WW])+abs(state->nb[NB_N]-state->nb[NB_NW])+abs(state->nb[NB_NE]-state->nb[NB_N]);//dh
	state->dv=abs(state->nb[NB_W]-state->nb[NB_NW])+abs(state->nb[NB_N]-state->nb[NB_NN])+abs(state->nb[NB_NE]-state->nb[NB_NNE]);//dv
	state->d45=abs(state->nb[NB_W]-state->nb[NB_NWW])+abs(state->nb[NB_NW]-state->nb[NB_NNWW])+abs(state->nb[NB_N]-state->nb[NB_NNW]);//d45
	state->d135=abs(state->nb[NB_NE]-state->nb[NB_NNEE])+abs(state->nb[NB_N]-state->nb[NB_NNE])+abs(state->nb[NB_W]-state->nb[NB_N]);//d135

	int scaled_dh, scaled_dv, scaled_d45, scaled_d135, scaled_eN, scaled_eW;
	if(state->lg_lambda<0)//depth probably > 8
	{
		scaled_dh  =state->dh  >>-state->lg_lambda;
		scaled_dv  =state->dv  >>-state->lg_lambda;
		scaled_d45 =state->d45 >>-state->lg_lambda;
		scaled_d135=state->d135>>-state->lg_lambda;
		scaled_eN  =state->nb[NB_eN]>>-state->lg_lambda;
		scaled_eW  =state->nb[NB_eW]>>-state->lg_lambda;
	}
	else if(state->lg_lambda>0)//depth probably < 8
	{
		scaled_dh  =state->dh  <<state->lg_lambda;
		scaled_dv  =state->dv  <<state->lg_lambda;
		scaled_d45 =state->d45 <<state->lg_lambda;
		scaled_d135=state->d135<<state->lg_lambda;
		scaled_eN  =state->nb[NB_eN]<<state->lg_lambda;
		scaled_eW  =state->nb[NB_eW]<<state->lg_lambda;
	}
	else
	{
		scaled_dh  =state->dh;
		scaled_dv  =state->dv;
		scaled_d45 =state->d45;
		scaled_d135=state->d135;
		scaled_eN  =state->nb[NB_eN];
		scaled_eW  =state->nb[NB_eW];
	}

	int sum=scaled_dv+state->dh, diff=scaled_dv-scaled_dh;//dv-dh
	if(sum>32)//sharp edge
		state->pred=(scaled_dv*state->nb[NB_W]+scaled_dh*state->nb[NB_N])/sum+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff>12)//horizontal edge
		state->pred=(2*state->nb[NB_W]+state->nb[NB_N])/3+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff<-12)//vertical edge
		state->pred=(state->nb[NB_W]+2*state->nb[NB_N])/3+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else//smooth area
		state->pred=(state->nb[NB_W]+state->nb[NB_N])/2+(state->nb[NB_NE]-state->nb[NB_NW])/8;

	diff=scaled_d45-scaled_d135;//d45-d135
	if(diff>32)//sharp 135-deg diagonal edge
		state->pred+=(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff>16)//135-deg diagonal edge
		state->pred+=(state->nb[NB_NE]-state->nb[NB_NW])/16;
	else if(diff<-32)//sharp 45-deg diagonal edge
		state->pred-=(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff<-16)//45-deg diagonal edge
		state->pred-=(state->nb[NB_NE]-state->nb[NB_NW])/16;

	state->energy=scaled_dh+scaled_dv+abs(scaled_eN)+abs(scaled_eW);
	//state->energy=scaled_dh+scaled_dv+abs(state->nb[NB_eN])+abs(state->nb[NB_eW]);
	if(state->energy<calic_qlevels[3])
	{
		if(state->energy<calic_qlevels[1])
			state->delta=state->energy<calic_qlevels[0]?0:1;
		else
			state->delta=state->energy<calic_qlevels[2]?2:3;
	}
	else
	{
		if(state->energy<calic_qlevels[5])
			state->delta=state->energy<calic_qlevels[4]?4:5;
		else
			state->delta=state->energy<calic_qlevels[6]?6:7;
	}
	state->beta=
		(state->nb[NB_W ]<state->pred)<<7|
		(state->nb[NB_NW]<state->pred)<<6|
		(state->nb[NB_N ]<state->pred)<<5|
		(state->nb[NB_NE]<state->pred)<<4|
		(state->nb[NB_WW]<state->pred)<<3|
		(state->nb[NB_NN]<state->pred)<<2|
		(2*state->nb[NB_N]-state->nb[NB_NN]<state->pred)<<1|
		(2*state->nb[NB_W]-state->nb[NB_WW]<state->pred);
	int ctx=(state->delta>>(11-CALIC_SSE_BITS))<<8|state->beta;
	int *cell=state->sse+ctx;
	int sse_count=*cell&0xFFF;
	int sse_sum=*cell>>12;
	state->sse_correction=sse_count?sse_sum/sse_count:0;
	state->pred+=state->sse_correction;
	state->pred=CLAMP(state->clamp_lo, state->pred, state->clamp_hi);
	//state->pred=MEDIAN3(state->nb[NB_N], state->nb[NB_W], state->pred);//clamp prediction

	int upred=state->pred+state->half;
	if(enc)
	{
		state->error=curr-state->pred;//curr in [-128, 127], error in [-128-pred, 127-pred]
		state->e2=state->sse_correction<0?-state->error:state->error;//to skew the histogram (predict the sign of error from SSE correction)
		if(state->e2)
		{
			//L/JXL permutation:		{0, -1,  1, -2,  2, ...}	(e<<1)^-(e<0)
			//CALIC/FLIF permutation:	{0,  1, -1,  2, -2, ...}	(e<0?-2*e:2*e-1) == (abs(e)<<1)-(e>0)
			if((upred<state->half)!=(state->sse_correction<0))
			{
				if(abs(state->e2)<=upred)
					state->e2=state->e2<<1^-(state->e2<0);
				else
					state->e2=upred+abs(state->e2);
			}
			else
			{
				upred=(state->half<<1)-upred;
				if(abs(state->e2)<=upred)
					state->e2=state->e2<<1^-(state->e2<0);
				else
					state->e2=upred+abs(state->e2);
			}
		}
	}
	
	int bypass=0, e2=enc?state->e2:0, delta=state->delta;
	unsigned *hist, nlevels, *CDF, fmin, f2;
	int sym;
	if(state->nbypass[delta])
	{
		int nbits=state->nbypass[delta];
		if(enc)
		{
			if(state->loud==2)
				state->csizes[state->kc]+=nbits;
			bypass=e2&((1<<nbits)-1);
			e2>>=nbits;
			while(nbits>8)
			{
				ac_enc(&state->ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 0x10000>>8);
				nbits-=8;
			}
			ac_enc(&state->ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 0x10000>>nbits);
		}
		else
		{
			while(nbits>8)
			{
				nbits-=8;
				bypass|=ac_dec(&state->ec, 0, 1<<8, 0x10000>>8)<<nbits;
			}
			bypass|=ac_dec(&state->ec, 0, 1<<nbits, 0x10000>>nbits);
		}
	}
	do//encode e2
	{
		if(delta<8)
		{
			CDF=state->CDF;
			hist=state->ct_tables[delta];
			nlevels=state->ct_tsizes[delta];

			//integrate & quantize CDF		for adaptive non-emitted CDFs
			sum=0;
			for(int ks=0;ks<(int)nlevels;++ks)
				sum+=hist[ks];
			int c=hist[0];
			CDF[0]=0;
			fmin=0;
			for(int ks=1;ks<(int)nlevels;++ks)
			{
				int freq=hist[ks];
				CDF[ks]=(unsigned)(c*(0x10000LL-nlevels)/sum)+ks;
				f2=CDF[ks]-CDF[ks-1];
				if(!fmin||fmin>f2)
					fmin=f2;
				c+=freq;
			}
			CDF[nlevels]=0x10000;
			f2=CDF[nlevels]-CDF[nlevels-1];
			if(!fmin||fmin>f2)
				fmin=f2;
		}
		else//bypass
		{
			CDF=0;
			nlevels=256;//tune this
			fmin=0x10000/256;
		}

		if(enc)
		{
			if(e2>=(int)nlevels-1)
			{
				sym=nlevels-1;//escape symbol
				e2-=sym;
			}
			else
				sym=e2;

			if((unsigned)sym>(unsigned)nlevels)//
				LOG_ERROR("Symbol OOB");

			ac_enc(&state->ec, sym, CDF, nlevels, fmin);
			if(state->loud==2)
				state->csizes[state->kc]+=calc_bitsize(CDF, nlevels, sym);
		}
		else
		{
			sym=ac_dec(&state->ec, CDF, nlevels, fmin);
			e2+=sym;
		}

		if(CDF)
		{
			hist[sym]+=state->ct_inc[delta];//update hist
			sum+=state->ct_inc[delta];
			if(sum>0x4000)//rescale
			{
				for(int ks=0;ks<(int)nlevels;++ks)
					hist[ks]=(hist[ks]+1)>>1;//ceil_half
				state->ct_inc[delta]>>=state->ct_inc[delta]>1;
			}
		}

		++delta;
	}while(sym==nlevels-1);

	if(!enc)
	{
		int nbits=state->nbypass[state->delta];
		e2<<=nbits;
		e2|=bypass;
		if((upred<state->half)!=(state->sse_correction<0))
		{
			if(e2<=(upred<<1))
				state->e2=e2>>1^-(e2&1);
			else
				state->e2=e2-upred;
		}
		else
		{
			upred=(state->half<<1)-upred;
			if(e2<=(upred<<1))
				state->e2=e2>>1^-(e2&1);
			else
				state->e2=upred-e2;
		}
		state->error=state->sse_correction<0?-state->e2:state->e2;
		curr=state->error+state->pred;
	}

	sse_sum+=state->error;//update SSE table
	++sse_count;
	if(sse_count>640)
	{
		sse_count>>=1;
		sse_sum>>=1;
	}
	*cell=sse_sum<<12|sse_count;
	state->currrow_sum_absE+=abs(state->error);
	return curr;
}
static int calic_bin(CalicState *state, int sym, int enc)//binary mode: 2 symbols in context (sparse context)
{
	state->beta=
		(state->nb2[5]!=state->nb2[0])<<4|
		(state->nb2[4]!=state->nb2[0])<<3|
		(state->nb2[3]!=state->nb2[0])<<2|
		(state->nb2[2]!=state->nb2[0])<<1|
		(state->nb2[1]!=state->nb2[0]);
	unsigned *hist=state->tables->ctr_bin[state->beta];
	unsigned sum=hist[0]+hist[1]+hist[2], c=0, fmin, f2;
	state->CDF[0]=0, c+=hist[0];
	state->CDF[1]=(int)(c*(0x10000LL-4)/sum)+1, c+=hist[1];		fmin=state->CDF[1]-state->CDF[0];
	state->CDF[2]=(int)(c*(0x10000LL-4)/sum)+2, c+=hist[2];		f2=state->CDF[2]-state->CDF[1];		if(fmin>f2)fmin=f2;
	state->CDF[3]=0x10000;						f2=state->CDF[3]-state->CDF[2];		if(fmin>f2)fmin=f2;
	
	if(enc)
	{
		ac_enc(&state->ec, sym, state->CDF, 3, fmin);
		if(state->loud==2)
			state->csizes[state->kc]+=calc_bitsize(state->CDF, 3, sym);
	}
	else
		sym=ac_dec(&state->ec, state->CDF, 3, fmin);
	
	hist[sym]+=state->bin_inc[state->beta];
	sum+=state->bin_inc[state->beta];
	if(sum>0x4000)//rescale
	{
		hist[0]=(hist[0]+1)>>1;//ceil_half
		hist[1]=(hist[1]+1)>>1;
		hist[2]=(hist[2]+1)>>1;
		state->bin_inc[state->beta]>>=state->bin_inc[state->beta]>1;
	}
	return sym;
}
int t45_encode(Image const *src, ArrayHandle *data, int loud)
{
	int res=src->iw*src->ih;
	double t_start=time_sec();
	int maxdepth=calc_maxdepth(src, 0);
	double bpp=
		((double)src->src_depth[0]+src->src_depth[1]+src->src_depth[2]+src->src_depth[3])/
		(8*((src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0)));
	double usize=image_getBMPsize(src);
	int nch=src->nch;
	//int nch=get_nch32(src->data, res);//FIXME: differentiate between just gray and just alpha
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T45 CALIC  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	if(!nch)
	{
		array_append(data, src->data, 1, 4, 1, 0, 0);
		return 1;
	}
	Image *im2=0;
	image_copy(&im2, src);
	if(!im2)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	rct_JPEG2000_32(im2, 1);
	{
		char temp;
		ROTATE3(im2->depth[0], im2->depth[1], im2->depth[2], temp);
		im2->depth[1]+=im2->depth[1]<24;
		im2->depth[2]+=im2->depth[2]<24;
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back1(&list, &nch);
	
	CalicState state;
	calic_init(&state, im2);
	ac_enc_init(&state.ec, &list);
	state.loud=loud;

	for(int kc=0;kc<nch;++kc)
	{
		calic_nextchannel(&state, kc);
		for(int ky=0;ky<src->ih;++ky)
		{
			calic_nextrow(&state);
			for(int kx=0;kx<src->iw;++kx)
			{
				//if(kc==0&&kx==1678&&ky==326)//
				//if(kc==0&&kx==1835&&ky==326)//
				//if(kc==0&&kx==1678&&ky==326)//
				//	printf("");

				int idx=(src->iw*ky+kx)<<2|kc;
				int curr=im2->data[idx];
				calic_prepctx(&state, kc, kx, ky);
				if(state.nunique<=2)
				{
					if(curr==state.unb[0])
					{
						calic_bin(&state, 0, 1);
						state.errors->data[idx]=0;
					}
					else if(state.nunique==2&&curr==state.unb[1])
					{
						calic_bin(&state, 1, 1);
						state.errors->data[idx]=0;
					}
					else
					{
						calic_bin(&state, 2, 1);
						calic_ct(&state, curr, 1);
						state.errors->data[idx]=state.error;
					}
				}
				else
				{
					calic_ct(&state, curr, 1);
					state.errors->data[idx]=state.error;
				}
			}
			if(loud)
			{
				int prog=kc*src->ih+ky+1;
				printf("%6.2lf%%  CR %10.6lf\r", (double)prog*100/(nch*src->ih), (double)prog*src->iw*bpp/list.nobj);
			}
		}
	}
	ac_enc_flush(&state.ec);
	dlist_appendtoarray(&list, data);

	if(loud)
	{
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		
		if(loud==2)
		{
			double csize=state.csizes[0]+state.csizes[1]+state.csizes[2]+state.csizes[3];
			printf("%-*.*s csize %12.2lf", 5, nch+1, "TYUVA", csize/8);
			for(int kc=0;kc<nch;++kc)
				printf(" %12.2lf", state.csizes[kc]/8);
			printf("\n");
			printf("%-*.*s CR    %12.6lf", 5, nch+1, "TYUVA", usize*8/csize);
			for(int kc=0;kc<nch;++kc)
				printf(" %12.6lf", res*8/state.csizes[kc]);
			printf("\n");
		}

		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);

		double proper_csizes[4]={0};
		calc_csize(state.errors, proper_csizes);
		printf("Proper csizes TYUVA %lf %lf %lf %lf %lf\n",
			proper_csizes[0]+proper_csizes[1]+proper_csizes[2]+proper_csizes[3],
			proper_csizes[0], proper_csizes[1], proper_csizes[2], proper_csizes[3]
		);
	}

	dlist_clear(&list);
	calic_free(&state);
	free(im2);
	return 1;
}
int t45_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	int res=dst->iw*dst->ih;
	double t_start=time_sec();

	if(srclen<4)
		return 0;
	int nch=data[0], comps[4]={0, 0, 0, 0};
	if(srclen==4)
	{
		memcpy(&nch, data, sizeof(int));
		comps[0]=nch;
		comps[1]=nch;
		comps[2]=nch;
		comps[3]=0;
		memfill(dst->data, comps, res*sizeof(int[4]), sizeof(comps));
		return 1;
	}
	++data;//skip tag
	--srclen;

	memfill(dst->data, comps, res*sizeof(int[4]), sizeof(comps));
	{
		char temp;
		ROTATE3(dst->depth[0], dst->depth[1], dst->depth[2], temp);
		dst->depth[1]+=dst->depth[1]<24;
		dst->depth[2]+=dst->depth[2]<24;
	}
	
	CalicState state;
	calic_init(&state, dst);
	ac_dec_init(&state.ec, data, data+srclen);

	for(int kc=0;kc<nch;++kc)
	{
		calic_nextchannel(&state, kc);
		for(int ky=0;ky<dst->ih;++ky)
		{
			calic_nextrow(&state);
			for(int kx=0;kx<dst->iw;++kx)
			{
				//if(kc==0&&kx==1678&&ky==326)//
				//if(kc==0&&kx==1835&&ky==326)//
				//if(kc==0&&kx==1678&&ky==326)//
				//	printf("");

				int idx=(dst->iw*ky+kx)<<2|kc;
				calic_prepctx(&state, kc, kx, ky);
				if(state.nunique<=2)
				{
					int T=calic_bin(&state, 0, 0);
					switch(T)
					{
					case 0:
						dst->data[idx]=state.unb[0];
						state.errors->data[idx]=0;
						break;
					case 1:
						if(state.nunique<2)
						{
							LOG_ERROR("Decode error CXY %d %d %d", kc, kx, ky);
							return 0;
						}
						dst->data[idx]=state.unb[1];
						state.errors->data[idx]=0;
						break;
					case 2:
						dst->data[idx]=calic_ct(&state, 0, 0);
						state.errors->data[idx]=state.error;
						break;
					}
				}
				else
				{
					dst->data[idx]=calic_ct(&state, 0, 0);
					state.errors->data[idx]=state.error;
				}
				//if(guide&&dst->data[idx]!=guide->data[idx])
				//	LOG_ERROR("Decode error CXY %d %d %d", kc, kx, ky);
			}
		}
	}
	{
		char temp;
		dst->depth[1]-=dst->depth[1]>dst->src_depth[1];
		dst->depth[2]-=dst->depth[2]>dst->src_depth[2];
		ROTATE3(dst->depth[2], dst->depth[1], dst->depth[0], temp);
	}
	rct_JPEG2000_32(dst, 0);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	calic_free(&state);
	return 1;
}