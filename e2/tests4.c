#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;


//static double p_sum=0;
//static int p_count=0;
double calc_bitsize(unsigned *CDF, int nlevels, int sym)
{
	int freq=CDF?CDF[sym+1]-CDF[sym]:0x10000/nlevels;
	//if(!freq)
	//	LOG_ERROR("ZPS");
	double prob=(double)freq/0x10000;
	double bitsize=-log2(prob);

	//p_sum+=prob;
	//++p_count;
	//if(!(p_count&0xFF))
	//	printf("%lf\n", p_sum/p_count);

	return bitsize;
}


//T45 CALIC

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
	char *pixels, *errors;
	int iw, ih;
	char nb[NB_COUNT], nb2[6], unb[6];
	int nunique, dh, dv, d45, d135, energy, delta, beta, sse_correction, pred, error, e2;
	unsigned CDF[257];
	int *sse;
	EntropyTables *tables;
	unsigned *ct_tables[8], ct_tsizes[8], ct_inc[8], bin_inc[32];
	ArithmeticCoder ec;

	int loud, kc;
	double csizes[4];
} CalicState;
static void calic_init(CalicState *state, char *pixels, int iw, int ih)
{
	int res=iw*ih;
	memset(state, 0, sizeof(*state));
	state->pixels=pixels;
	state->iw=iw;
	state->ih=ih;

	state->errors=(char*)malloc((size_t)res<<2);
	state->sse=(int*)malloc(sizeof(int[1024]));
	state->tables=(EntropyTables*)malloc(sizeof(EntropyTables));
	if(!state->errors||!state->sse||!state->tables)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(state->errors, 0, (size_t)res<<2);

	state->ct_tables[0]=state->tables->ctr_ct0, state->ct_tsizes[0]=_countof(state->tables->ctr_ct0);
	state->ct_tables[1]=state->tables->ctr_ct1, state->ct_tsizes[1]=_countof(state->tables->ctr_ct1);
	state->ct_tables[2]=state->tables->ctr_ct2, state->ct_tsizes[2]=_countof(state->tables->ctr_ct2);
	state->ct_tables[3]=state->tables->ctr_ct3, state->ct_tsizes[3]=_countof(state->tables->ctr_ct3);
	state->ct_tables[4]=state->tables->ctr_ct4, state->ct_tsizes[4]=_countof(state->tables->ctr_ct4);
	state->ct_tables[5]=state->tables->ctr_ct5, state->ct_tsizes[5]=_countof(state->tables->ctr_ct5);
	state->ct_tables[6]=state->tables->ctr_ct6, state->ct_tsizes[6]=_countof(state->tables->ctr_ct6);
	state->ct_tables[7]=state->tables->ctr_ct7, state->ct_tsizes[7]=_countof(state->tables->ctr_ct7);
}
static void calic_resetfornextchannel(CalicState *state, int kc)
{
	int fillval=1;
	memfill(state->tables, &fillval, sizeof(*state->tables), sizeof(int));
	fillval=0x2000;
	memfill(state->ct_inc, &fillval, sizeof(state->ct_inc), sizeof(int));
	memfill(state->bin_inc, &fillval, sizeof(state->bin_inc), sizeof(int));
	memset(state->sse, 0, sizeof(int[1024]));
	state->kc=kc;
}
static void calic_free(CalicState *state)
{
	free(state->errors);
	free(state->sse);
	free(state->tables);
}
static void calic_prepctx(CalicState *state, int kc, int kx, int ky)
{
	int idx=0;
	for(int ky2=-CALIC_REACH;ky2<=0;++ky2)
	{
		for(int kx2=-CALIC_REACH;kx2<=CALIC_REACH;++kx2, ++idx)
		{
			if(!ky2&&!kx2)
				break;
			if((unsigned)(ky+ky2)<(unsigned)state->ih&&(unsigned)(kx+kx2)<(unsigned)state->iw)
				state->nb[idx]=state->pixels[(state->iw*(ky+ky2)+kx+kx2)<<2|kc];
			else
				state->nb[idx]=0;
		}
	}
	state->nb[idx++]=ky?state->errors[(state->iw*(ky-1)+kx)<<2|kc]:0;
	state->nb[idx++]=kx?state->errors[(state->iw*ky+kx-1)<<2|kc]:0;

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
		char val=state->nb2[k];
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

	int sum=state->dv+state->dh, diff=state->dv-state->dh;//dv-dh
	if(sum>32)//sharp edge
		state->pred=(state->dv*state->nb[NB_W]+state->dh*state->nb[NB_N])/sum+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff>12)//horizontal edge
		state->pred=(2*state->nb[NB_W]+state->nb[NB_N])/3+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff<-12)//vertical edge
		state->pred=(state->nb[NB_W]+2*state->nb[NB_N])/3+(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else//smooth area
		state->pred=(state->nb[NB_W]+state->nb[NB_N])/2+(state->nb[NB_NE]-state->nb[NB_NW])/8;

	diff=state->d45-state->d135;//d45-d135
	if(diff>32)//sharp 135-deg diagonal edge
		state->pred+=(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff>16)//135-deg diagonal edge
		state->pred+=(state->nb[NB_NE]-state->nb[NB_NW])/16;
	else if(diff<-32)//sharp 45-deg diagonal edge
		state->pred-=(state->nb[NB_NE]-state->nb[NB_NW])/8;
	else if(diff<-16)//45-deg diagonal edge
		state->pred-=(state->nb[NB_NE]-state->nb[NB_NW])/16;

	state->energy=state->dh+state->dv+abs(state->nb[NB_eW])+abs(state->nb[NB_eN]);
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
	int ctx=(state->delta>>1)<<8|state->beta;
	int *cell=state->sse+ctx;
	int count=*cell&0xFFF;
	sum=*cell>>12;
	state->sse_correction=count?sum/count:0;
	state->pred+=state->sse_correction;
	state->pred=MEDIAN3(state->nb[NB_N], state->nb[NB_W], state->pred);//clamp prediction

	if(enc)
	{
		state->error=curr-state->pred;//curr in [-128, 127], error in [-128-pred, 127-pred]
	
		state->e2=state->sse_correction<0?-state->error:state->error;//to skew the histogram (predict the sign of error from SSE correction)

		if(state->e2)
		{
			//L/JXL permutation:		{0, -1,  1, -2,  2, ...}	(e<<1)^-(e<0)
			//CALIC/FLIF permutation:	{0,  1, -1,  2, -2, ...}	(e<0?-2*e:2*e-1) == (abs(e)<<1)-(e>0)
			int upred=state->pred+128;
			if(upred<128)
			{
				if(abs(state->e2)<=upred)
					state->e2=state->e2<<1^-(state->e2<0);
				else
					state->e2+=upred;
			}
			else
			{
				upred=256-upred;
				if(abs(state->e2)<=upred)
					state->e2=state->e2<<1^-(state->e2<0);
				else
					state->e2=upred-state->e2;
			}
		}
	}
	
	int e2=enc?state->e2:0, delta=state->delta;
	unsigned *hist, nlevels, *CDF;
	int sym;
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
			for(int ks=0, c=0;ks<(int)nlevels;++ks)
			{
				int freq=hist[ks];
				CDF[ks]=(unsigned)(c*(0x10000LL-nlevels)/sum)+ks;
				c+=freq;
			}
			CDF[nlevels]=0x10000;

			//int nused=0, sum=0;		//X  this is designed for static-emitted CDFs
			//for(int ks=0;ks<tsize;++ks)
			//{
			//	sum+=table[ks];
			//	nused=table[ks]!=0;
			//}
			//for(int ks=0, ks2=0, c=0;ks<tsize;++ks)
			//{
			//	int freq=table[ks];
			//	CDF[ks]=(unsigned)(c*(0x10000LL-nused)/sum)+ks2;
			//	c+=freq;
			//	ks2+=freq!=0;
			//}
			//CDF[nlevels]=0x10000;
		}
		else//bypass
		{
			CDF=0;
			nlevels=32;//tune this
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

			ac_enc(&state->ec, sym, CDF, nlevels);
			if(state->loud)
				state->csizes[state->kc]+=calc_bitsize(CDF, nlevels, sym);
		}
		else
		{
			sym=ac_dec(&state->ec, CDF, nlevels);
			e2+=sym;
		}

		if(CDF)
		{
#if 1
			hist[sym]+=state->ct_inc[delta];//update hist
			sum+=state->ct_inc[delta];
			if(sum>0x4000)//rescale
			{
				for(int ks=0;ks<(int)nlevels;++ks)
					hist[ks]=(hist[ks]+1)>>1;//ceil_half
				state->ct_inc[delta]>>=state->ct_inc[delta]>1;
			}
#endif
			//++hist[sym];
		}

		++delta;
	}while(sym==nlevels-1);

	if(!enc)
	{
		int upred=state->pred+128;
		if(upred<128)
		{
			if(e2<=(upred<<1))
				state->e2=e2>>1^-(e2&1);
			else
				state->e2=e2-upred;
		}
		else
		{
			upred=256-upred;
			if(e2<=(upred<<1))
				state->e2=e2>>1^-(e2&1);
			else
				state->e2=upred-e2;
		}
		state->error=state->sse_correction<0?-state->e2:state->e2;
		curr=state->error+state->pred;
	}

	sum+=state->error;//update SSE table
	++count;
	if(count>640)
	{
		count>>=1;
		sum>>=1;
	}
	*cell=sum<<12|count;
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
	unsigned sum=hist[0]+hist[1]+hist[2], c=0;
	state->CDF[0]=0, c+=hist[0];
	state->CDF[1]=(int)(c*(0x10000LL-4)/sum)+1, c+=hist[1];
	state->CDF[2]=(int)(c*(0x10000LL-4)/sum)+2, c+=hist[2];
	state->CDF[3]=0x10000;
	
	if(enc)
	{
		ac_enc(&state->ec, sym, state->CDF, 3);
		if(state->loud)
			state->csizes[state->kc]+=calc_bitsize(state->CDF, 3, sym);
	}
	else
		sym=ac_dec(&state->ec, state->CDF, 3);

#if 1
	hist[sym]+=state->bin_inc[state->beta];
	sum+=state->bin_inc[state->beta];
	if(sum>0x4000)//rescale
	{
		hist[0]=(hist[0]+1)>>1;//ceil_half
		hist[1]=(hist[1]+1)>>1;
		hist[2]=(hist[2]+1)>>1;
		state->bin_inc[state->beta]>>=state->bin_inc[state->beta]>1;
	}
#endif
	//++hist[sym];
	return sym;
}
int t45_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T45 CALIC  Enc %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int nch=get_nch(src, iw*ih);//FIXME: differentiate between just gray and just alpha
	if(!nch)
	{
		//printf("Skipping image filled with solid 0x%02X\n", src[0]);
		if(*data)
		{
			if(data[0]->esize!=1)
				return 0;
			ARRAY_APPEND(*data, src, 1, 1, 0);
		}
		else
			ARRAY_ALLOC(char, *data, src, 1, 0, 0);
		return 1;
	}
	char *pixels=(char*)malloc((size_t)res<<2);
	if(!pixels)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(pixels, src, (size_t)res<<2);
	addbuf((unsigned char*)pixels, iw, ih, nch==1?3:nch, 4, 128);
	colortransform_JPEG2000_fwd(pixels, iw, ih);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back1(&list, &nch);
	
	CalicState state;
	calic_init(&state, pixels, iw, ih);
	ac_enc_init(&state.ec, &list);
	state.loud=loud;

	for(int kc=0;kc<nch;++kc)
	{
		calic_resetfornextchannel(&state, kc);
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				//if(kc==0&&kx==647&&ky==355)//
				//if(kc==0&&kx==648&&ky==354)//
				//if(kc==0&&kx==128&&ky==1)//
				//if(kc==0&&kx==507&&ky==98)//
				//	printf("");

				int idx=(state.iw*ky+kx)<<2|kc;
				char curr=state.pixels[idx];
				calic_prepctx(&state, kc, kx, ky);
				if(state.nunique<=2)
				{
					if(curr==state.unb[0])
					{
						calic_bin(&state, 0, 1);
						state.errors[idx]=0;
					}
					else if(state.nunique==2&&curr==state.unb[1])
					{
						calic_bin(&state, 1, 1);
						state.errors[idx]=0;
					}
					else
					{
						calic_bin(&state, 2, 1);
						calic_ct(&state, curr, 1);
						state.errors[idx]=state.error;
					}
				}
				else
				{
					calic_ct(&state, curr, 1);
					state.errors[idx]=state.error;
				}
			}
			if(loud)
				printf("%5.2lf%%  CR %10.6lf\r", (double)(ky+1)*100/ih, (double)(ky+1)*iw*nch/list.nobj);
		}
	}
	ac_enc_flush(&state.ec);
	size_t dststart=dlist_appendtoarray(&list, data);

	if(loud)
	{
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		
		size_t usize=(size_t)iw*ih*nch;
		double csize=state.csizes[0]+state.csizes[1]+state.csizes[2]+state.csizes[3];
		printf("%-*.*s %10.6lf", 5, nch+1, "TYUVA", usize*8/csize);
		for(int kc=0;kc<nch;++kc)
			printf(" %10.6lf", res*8/state.csizes[kc]);
		printf("\n");

		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, (double)usize/list.nobj);
	}

	dlist_clear(&list);
	calic_free(&state);
	free(pixels);
	return 1;
}
int t45_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *pixels, int loud)
{
	int res=iw*ih;
	double t_start=time_sec();

	if(!srclen)
		return 0;
	unsigned char nch=data[0];
	if(srclen==1)
	{
		int color=0xFF000000|nch<<16|nch<<8|nch;
		memfill(pixels, &color, sizeof(int), (size_t)iw*ih<<2);
		return 1;
	}
	++data;//skip tag
	--srclen;

	char *errors=(char*)malloc((size_t)res<<2);
	if(!errors)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	int black=0xFF000000;
	memfill(pixels, &black, res*sizeof(int), sizeof(int));
	
	CalicState state;
	calic_init(&state, pixels, iw, ih);
	ac_dec_init(&state.ec, data, data+srclen);

	for(int kc=0;kc<nch;++kc)
	{
		calic_resetfornextchannel(&state, kc);
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				//if(kc==0&&kx==647&&ky==355)//
				//if(kc==0&&kx==2087&&ky==1570)//
				//if(kc==0&&kx==507&&ky==98)//
				//	printf("");

				int idx=(state.iw*ky+kx)<<2|kc;
				calic_prepctx(&state, kc, kx, ky);
				if(state.nunique<=2)
				{
					int T=calic_bin(&state, 0, 0);
					switch(T)
					{
					case 0:
						pixels[idx]=state.unb[0];
						state.errors[idx]=0;
						break;
					case 1:
						if(state.nunique<2)
							LOG_ERROR("Decode error");
						pixels[idx]=state.unb[1];
						state.errors[idx]=0;
						break;
					case 2:
						pixels[idx]=calic_ct(&state, 0, 0);
						state.errors[idx]=state.error;
						break;
					}
				}
				else
				{
					pixels[idx]=calic_ct(&state, 0, 0);
					state.errors[idx]=state.error;
				}
			}
		}
	}
	colortransform_JPEG2000_inv(pixels, iw, ih);
	addbuf(pixels, iw, ih, nch==1?3:nch, 4, 128);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	calic_free(&state);
	free(errors);
	return 1;
}