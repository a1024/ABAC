#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
//#include"blist.h"
static const char file[]=__FILE__;


	#define BYPASS_ON_INFLATION
//	#define ENABLE_FILEGUARD	//makes using scripts harder

	#define USE_RCTJ2K
	#define USE_WG

//	#define ESTIMATE_SIZE		//FOR DEBUGGING
//	#define AC_VALIDATE		//FOR DEBUGGING

#ifdef AC_VALIDATE
#define AC_IMPLEMENTATION
#include"entropy.h"
#endif

#define PADX 8
#define PADY 4

#ifdef USE_WG
#define WG_NPREDS	8	//multiple of 4
#define WG_PREDLIST0\
	WG_PRED(340,	N-eN/3)\
	WG_PRED(340,	W-eW/3)\
	WG_PRED(205,	3*(N-NN)+NNN-eN/6-eNN/6+eNNN*2/3)\
	WG_PRED(205,	3*(W-WW)+WWW-eW/6-eWW/6+eWWW*2/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW/4-(eW>>7))\
	WG_PRED(240,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW-(4*(eN+eW)+eNN+eWW)/2)/3)\
	WG_PRED(120,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(120,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#define WG_PREDLIST1\
	WG_PRED(330,	N+(2*eN+eW)/6)\
	WG_PRED(330,	W+(2*eW+eN)/6)\
	WG_PRED(175,	3*(N-NN)+NNN+eN/6+eNN/6-eWW*2/3)\
	WG_PRED(175,	3*(W-WW)+WWW+eW/6+eWW/6-eNN*2/3)\
	WG_PRED(180,	W+NE-N-((eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+(2*eN+eNE)/10)
#define WG_PREDLIST2\
	WG_PRED(330,	N+(2*eN+eW)/6)\
	WG_PRED(330,	W+(2*eW+eN)/6)\
	WG_PRED(200,	3*(N-NN)+NNN+(eN-eWW)/3)\
	WG_PRED(200,	3*(W-WW)+WWW+(eW-eNN)/3)\
	WG_PRED(180,	W+NE-N-((5*eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/8)\
	WG_PRED(150,	N+NE-NNE+(2*eN+eNE)/10)
#endif
static void update_CDF(const unsigned *hist, int nlevels, unsigned *CDF)
{
	const int hsum=hist[nlevels];
	for(int ks=0, c=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<16)-nlevels)/hsum)+ks;
		c+=freq;
	}
	CDF[nlevels]=1<<16;
}
static void rescale_hist(unsigned *hist, int nlevels)
{
	int hsum=0;
	for(int ks=0;ks<nlevels;++ks)
		hsum+=hist[ks]=(hist[ks]+1)>>1;
	hist[nlevels]=hsum;
}
int c15_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t dstsize;
#ifdef ENABLE_FILEGUARD
	dstsize=get_filesize(dstfn);
	if(dstsize>=0)
	{
		LOG_ERROR("Destination file already exists");
		return 1;
	}
#endif
	FILE *fsrc=fopen(srcfn, "rb");
	int tag=0;
	size_t nread=fread(&tag, 1, 2, fsrc), nwritten=0;
	if(nread!=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int fwd=tag==('P'|'6'<<8);
	int bypass=tag==('B'|'P'<<8);
	if(!fwd&&tag!=('C'|'H'<<8)&&!bypass)
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	int usize=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
	}
	else
	{
		nread=fread(&iw, 1, 4, fsrc);
		nread+=fread(&ih, 1, 4, fsrc);
		if(nread!=4LL+4)
		{
			LOG_ERROR("Unsupported archive");
			return 1;
		}
	}
	usize=iw*ih*3;
	unsigned char *buffer=(unsigned char*)malloc(usize);
	if(!buffer)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	FILE *fdst=fopen(dstfn, "wb");
	nwritten=0;
	if(fwd)
	{
		nread=fread(buffer, 1, usize, fsrc);
		if(nread!=usize)
		{
			LOG_ERROR("Corrupt PPM");
			return 1;
		}

		nwritten+=fwrite("CH", 1, 2, fdst);
		nwritten+=fwrite(&iw, 1, 4, fdst);
		nwritten+=fwrite(&ih, 1, 4, fdst);
	}
	else
	{
		nwritten+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		if(bypass)
		{
			printf("Bypass decode\n");
			nread=fread(buffer, 1, usize, fsrc);
			if(nread!=usize)
				printf("Error  read %d/%d bytes\n", (int)nread, usize);
			fwrite(buffer, 1, usize, fdst);
			fclose(fsrc);
			fclose(fdst);
			free(buffer);
			return 0;
		}
		memset(buffer, 0, usize);
	}
#ifdef USE_WG
	ALIGN(32) double wg_weights[3][WG_NPREDS]={0};
	ALIGN(32) int wg_perrors[3][WG_NPREDS]={0}, wg_preds[3][WG_NPREDS]={0};
	{
		int j;

		j=0;
#define WG_PRED(WEIGHT, EXPR) wg_weights[0][j++]=WEIGHT;
		WG_PREDLIST0
#undef  WG_PRED
		j=0;
#define WG_PRED(WEIGHT, EXPR) wg_weights[1][j++]=WEIGHT;
		WG_PREDLIST1
#undef  WG_PRED
		j=0;
#define WG_PRED(WEIGHT, EXPR) wg_weights[2][j++]=WEIGHT;
		WG_PREDLIST2
#undef  WG_PRED
	}
#endif
	int psize=(iw+16LL)*(int)sizeof(int[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
	int *pixels=(int*)_mm_malloc(psize, sizeof(__m128i));
#ifdef USE_WG
	int ebufsize=(iw+PADX*2)*(int)sizeof(int[PADY*4*WG_NPREDS]);//PADY padded rows * 4 channels max * WG_NPREDS pred errors
	int *ebuf=(int*)_mm_malloc(ebufsize, sizeof(__m128i));
#endif
	int hsize=(int)sizeof(int[257+513+513]);
	unsigned *hist=(unsigned*)malloc(hsize);
	unsigned *CDF=(unsigned*)malloc(hsize);
	int sbufsize=iw*16;//rowsize = iw*3, allocate iw*16 just in case
	unsigned short *sbuf=(unsigned short*)malloc(sbufsize);//stream buffer
	if(!pixels||!hist||!CDF||!sbuf
#ifdef USE_WG
		||!ebuf
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(sbuf, 0, sbufsize);
	memset(hist, 0, hsize);
	for(int k=0;k<257;++k)
		CDF[k]=k<<8;
	for(int ks=0;ks<513;++ks)
	{
		CDF[257+ks]=ks<<7;
		CDF[257+513+ks]=ks<<7;
	}
#ifdef USE_WG
	memset(ebuf, 0, ebufsize);
#endif
	memset(pixels, 0, psize);
	int rowsize=iw*3;
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) int *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(32) int *erows[]=
		{
			ebuf+((iw+16LL)*((ky-0LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((iw+16LL)*((ky-1LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((iw+16LL)*((ky-2LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((iw+16LL)*((ky-3LL)&3)+8LL)*4*WG_NPREDS,
		};
		ALIGN(16) int preds[4]={0};
		unsigned *curr_CDF=0;
		unsigned cdfs[3]={0}, freqs[3]={0};
		//unsigned c2;
		int yuv[3]={0}, errors[3]={0};
		unsigned long long low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef AC_VALIDATE
		unsigned long long lo0, r0;
#endif
		unsigned short *sptr=sbuf;
		int streamsize=0;
		if(!fwd)
		{
			nread=fread(&streamsize, 1, 4, fsrc);
			if(streamsize)
			{
				nread=fread(sbuf, sizeof(short), streamsize, fsrc);
				code=*sptr++;
				code=*sptr++|code<<16;
				code=*sptr++|code<<16;
			}
			else
				nread=fread(buffer+rowsize*ky, 1, rowsize, fsrc);
		}
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
#ifdef USE_WG
			int j;
			{
				int
					NNNN	=rows[0][0+0*4*2+0],
					NNNWWWW	=rows[3][0-4*4*2+0],
					NNNWWW	=rows[3][0-3*4*2+0],
					NNNW	=rows[3][0-1*4*2+0],
					NNN	=rows[3][0+0*4*2+0],
					NNNE	=rows[3][0+1*4*2+0],
					NNNEEE	=rows[3][0+3*4*2+0],
					NNWWWW	=rows[2][0-4*4*2+0],
					NNWW	=rows[2][0-2*4*2+0],
					NNW	=rows[2][0-1*4*2+0],
					NN	=rows[2][0+0*4*2+0],
					NNE	=rows[2][0+1*4*2+0],
					NNEE	=rows[2][0+2*4*2+0],
					NNEEE	=rows[2][0+3*4*2+0],
					NNEEEE	=rows[2][0+4*4*2+0],
					NWWW	=rows[1][0-3*4*2+0],
					NWW	=rows[1][0-2*4*2+0],
					NW	=rows[1][0-1*4*2+0],
					N	=rows[1][0+0*4*2+0],
					NE	=rows[1][0+1*4*2+0],
					NEE	=rows[1][0+2*4*2+0],
					NEEE	=rows[1][0+3*4*2+0],
					NEEEE	=rows[1][0+4*4*2+0],
					WWWWW	=rows[0][0-5*4*2+0],
					WWWW	=rows[0][0-4*4*2+0],
					WWW	=rows[0][0-3*4*2+0],
					WW	=rows[0][0-2*4*2+0],
					W	=rows[0][0-1*4*2+0],
					eNNN	=rows[3][0+0*4*2+4],
					eNN	=rows[2][0+0*4*2+4],
					eNNE	=rows[2][0+1*4*2+4],
					eNW	=rows[1][0-1*4*2+4],
					eN	=rows[1][0+0*4*2+4],
					eNE	=rows[1][0+1*4*2+4],
					eNEE	=rows[1][0+2*4*2+4],
					eNEEE	=rows[1][0+3*4*2+4],
					eWWWW	=rows[0][0-4*4*2+4],
					eWWW	=rows[0][0-3*4*2+4],
					eWW	=rows[0][0-2*4*2+4],
					eW	=rows[0][0-1*4*2+4];
				(void)NNNN;
				(void)NNNWWWW;
				(void)NNNWWW;
				(void)NNNW;
				(void)NNN;
				(void)NNNE;
				(void)NNNEEE;
				(void)NNWWWW;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NNEEE;
				(void)NNEEEE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNNN;
				(void)eNN;
				(void)eNNE;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eNEEE;
				(void)eWWWW;
				(void)eWWW;
				(void)eWW;
				(void)eW;
				j=0;
#define WG_PRED(WEIGHT, EXPR) wg_preds[0][j++]=EXPR;
				WG_PREDLIST0
#undef  WG_PRED
			}
			{
				int
					NNNN	=rows[0][1+0*4*2+0],
					NNNWWWW	=rows[3][1-4*4*2+0],
					NNNWWW	=rows[3][1-3*4*2+0],
					NNNW	=rows[3][1-1*4*2+0],
					NNN	=rows[3][1+0*4*2+0],
					NNNE	=rows[3][1+1*4*2+0],
					NNNEEE	=rows[3][1+3*4*2+0],
					NNWWWW	=rows[2][1-4*4*2+0],
					NNWW	=rows[2][1-2*4*2+0],
					NNW	=rows[2][1-1*4*2+0],
					NN	=rows[2][1+0*4*2+0],
					NNE	=rows[2][1+1*4*2+0],
					NNEE	=rows[2][1+2*4*2+0],
					NNEEE	=rows[2][1+3*4*2+0],
					NNEEEE	=rows[2][1+4*4*2+0],
					NWWW	=rows[1][1-3*4*2+0],
					NWW	=rows[1][1-2*4*2+0],
					NW	=rows[1][1-1*4*2+0],
					N	=rows[1][1+0*4*2+0],
					NE	=rows[1][1+1*4*2+0],
					NEE	=rows[1][1+2*4*2+0],
					NEEE	=rows[1][1+3*4*2+0],
					NEEEE	=rows[1][1+4*4*2+0],
					WWWWW	=rows[0][1-5*4*2+0],
					WWWW	=rows[0][1-4*4*2+0],
					WWW	=rows[0][1-3*4*2+0],
					WW	=rows[0][1-2*4*2+0],
					W	=rows[0][1-1*4*2+0],
					eNNN	=rows[3][1+0*4*2+4],
					eNN	=rows[2][1+0*4*2+4],
					eNNE	=rows[2][1+1*4*2+4],
					eNW	=rows[1][1-1*4*2+4],
					eN	=rows[1][1+0*4*2+4],
					eNE	=rows[1][1+1*4*2+4],
					eNEE	=rows[1][1+2*4*2+4],
					eNEEE	=rows[1][1+3*4*2+4],
					eWWWW	=rows[0][1-4*4*2+4],
					eWWW	=rows[0][1-3*4*2+4],
					eWW	=rows[0][1-2*4*2+4],
					eW	=rows[0][1-1*4*2+4];
				(void)NNNN;
				(void)NNNWWWW;
				(void)NNNWWW;
				(void)NNNW;
				(void)NNN;
				(void)NNNE;
				(void)NNNEEE;
				(void)NNWWWW;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NNEEE;
				(void)NNEEEE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNNN;
				(void)eNN;
				(void)eNNE;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eNEEE;
				(void)eWWWW;
				(void)eWWW;
				(void)eWW;
				(void)eW;
				j=0;
#define WG_PRED(WEIGHT, EXPR) wg_preds[1][j++]=EXPR;
				WG_PREDLIST1
#undef  WG_PRED
			}
			{
				int
					NNNN	=rows[0][2+0*4*2+0],
					NNNWWWW	=rows[3][2-4*4*2+0],
					NNNWWW	=rows[3][2-3*4*2+0],
					NNNW	=rows[3][2-1*4*2+0],
					NNN	=rows[3][2+0*4*2+0],
					NNNE	=rows[3][2+1*4*2+0],
					NNNEEE	=rows[3][2+3*4*2+0],
					NNWWWW	=rows[2][2-4*4*2+0],
					NNWW	=rows[2][2-2*4*2+0],
					NNW	=rows[2][2-1*4*2+0],
					NN	=rows[2][2+0*4*2+0],
					NNE	=rows[2][2+1*4*2+0],
					NNEE	=rows[2][2+2*4*2+0],
					NNEEE	=rows[2][2+3*4*2+0],
					NNEEEE	=rows[2][2+4*4*2+0],
					NWWW	=rows[1][2-3*4*2+0],
					NWW	=rows[1][2-2*4*2+0],
					NW	=rows[1][2-1*4*2+0],
					N	=rows[1][2+0*4*2+0],
					NE	=rows[1][2+1*4*2+0],
					NEE	=rows[1][2+2*4*2+0],
					NEEE	=rows[1][2+3*4*2+0],
					NEEEE	=rows[1][2+4*4*2+0],
					WWWWW	=rows[0][2-5*4*2+0],
					WWWW	=rows[0][2-4*4*2+0],
					WWW	=rows[0][2-3*4*2+0],
					WW	=rows[0][2-2*4*2+0],
					W	=rows[0][2-1*4*2+0],
					eNNN	=rows[3][2+0*4*2+4],
					eNN	=rows[2][2+0*4*2+4],
					eNNE	=rows[2][2+1*4*2+4],
					eNW	=rows[1][2-1*4*2+4],
					eN	=rows[1][2+0*4*2+4],
					eNE	=rows[1][2+1*4*2+4],
					eNEE	=rows[1][2+2*4*2+4],
					eNEEE	=rows[1][2+3*4*2+4],
					eWWWW	=rows[0][2-4*4*2+4],
					eWWW	=rows[0][2-3*4*2+4],
					eWW	=rows[0][2-2*4*2+4],
					eW	=rows[0][2-1*4*2+4];
				(void)NNNN;
				(void)NNNWWWW;
				(void)NNNWWW;
				(void)NNNW;
				(void)NNN;
				(void)NNNE;
				(void)NNNEEE;
				(void)NNWWWW;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NNEEE;
				(void)NNEEEE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNNN;
				(void)eNN;
				(void)eNNE;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eNEEE;
				(void)eWWWW;
				(void)eWWW;
				(void)eWW;
				(void)eW;
				j=0;
#define WG_PRED(WEIGHT, EXPR) wg_preds[2][j++]=EXPR;
				WG_PREDLIST2
#undef  WG_PRED
			}
			int
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			int
				*eNNE	=erows[2]+1*4*WG_NPREDS,
				*eNW	=erows[1]-1*4*WG_NPREDS,
				*eN	=erows[1]+0*4*WG_NPREDS,
				*eNE	=erows[1]+1*4*WG_NPREDS,
			//	*eNEE	=erows[1]+2*4*WG_NPREDS,
			//	*eNEEE	=erows[1]+3*4*WG_NPREDS,
				*eW	=erows[0]-1*4*WG_NPREDS,
				*ecurr	=erows[0]+0*4*WG_NPREDS;
			{
#define WG_KC 0
				__m128i one=_mm_set1_epi32(1);
				__m128i me0=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+0);
				__m128i me1=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+1);
			//	__m128i me2=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+2);
			//	__m128i me3=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+3);
				__m256d w0=_mm256_load_pd(wg_weights[WG_KC]+0*4);
				__m256d w1=_mm256_load_pd(wg_weights[WG_KC]+1*4);
			//	__m256d w2=_mm256_load_pd(wg_weights[WG_KC]+2*4);
			//	__m256d w3=_mm256_load_pd(wg_weights[WG_KC]+3*4);
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+0), 1));
				me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+1), 1));
			//	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+2), 1));
			//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+3), 1));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, one);
				me1=_mm_add_epi32(me1, one);
			//	me2=_mm_add_epi32(me2, one);
			//	me3=_mm_add_epi32(me3, one);
				__m256d de0=_mm256_cvtepi32_pd(me0);
				__m256d de1=_mm256_cvtepi32_pd(me1);
			//	__m256d de2=_mm256_cvtepi32_pd(me2);
			//	__m256d de3=_mm256_cvtepi32_pd(me3);
				w0=_mm256_div_pd(w0, de0);
				w1=_mm256_div_pd(w1, de1);
			//	w2=_mm256_div_pd(w2, de2);
			//	w3=_mm256_div_pd(w3, de3);
				de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+0));
				de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+1));
			//	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+2));
			//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+3));
				de0=_mm256_mul_pd(de0, w0);
				de1=_mm256_mul_pd(de1, w1);
			//	de2=_mm256_mul_pd(de2, w2);
			//	de3=_mm256_mul_pd(de3, w3);
				w0=_mm256_add_pd(w0, w1);
			//	w0=_mm256_add_pd(w0, w2);
			//	w0=_mm256_add_pd(w0, w3);
				de0=_mm256_add_pd(de0, de1);
			//	de0=_mm256_add_pd(de0, de2);
			//	de0=_mm256_add_pd(de0, de3);
				//[num3 num2 num1 num0]
				//[den3 den2 den1 den0]
				//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
				//lo=_mm256_extractf128_pd(r, 0)
				//hi=_mm256_extractf128_pd(r, 1)
				//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
				w0=_mm256_hadd_pd(de0, w0);
				__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
				dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
				__m128i mp=_mm_cvtpd_epi32(dp);
				_mm_store_ss((float*)preds+WG_KC, _mm_castsi128_ps(mp));
#undef  WG_KC
			}
			{
#define WG_KC 1
				__m128i one=_mm_set1_epi32(1);
				__m128i me0=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+0);
				__m128i me1=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+1);
			//	__m128i me2=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+2);
			//	__m128i me3=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+3);
				__m256d w0=_mm256_load_pd(wg_weights[WG_KC]+0*4);
				__m256d w1=_mm256_load_pd(wg_weights[WG_KC]+1*4);
			//	__m256d w2=_mm256_load_pd(wg_weights[WG_KC]+2*4);
			//	__m256d w3=_mm256_load_pd(wg_weights[WG_KC]+3*4);
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+0), 1));
				me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+1), 1));
			//	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+2), 1));
			//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+3), 1));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, one);
				me1=_mm_add_epi32(me1, one);
			//	me2=_mm_add_epi32(me2, one);
			//	me3=_mm_add_epi32(me3, one);
				__m256d de0=_mm256_cvtepi32_pd(me0);
				__m256d de1=_mm256_cvtepi32_pd(me1);
			//	__m256d de2=_mm256_cvtepi32_pd(me2);
			//	__m256d de3=_mm256_cvtepi32_pd(me3);
				w0=_mm256_div_pd(w0, de0);
				w1=_mm256_div_pd(w1, de1);
			//	w2=_mm256_div_pd(w2, de2);
			//	w3=_mm256_div_pd(w3, de3);
				de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+0));
				de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+1));
			//	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+2));
			//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+3));
				de0=_mm256_mul_pd(de0, w0);
				de1=_mm256_mul_pd(de1, w1);
			//	de2=_mm256_mul_pd(de2, w2);
			//	de3=_mm256_mul_pd(de3, w3);
				w0=_mm256_add_pd(w0, w1);
			//	w0=_mm256_add_pd(w0, w2);
			//	w0=_mm256_add_pd(w0, w3);
				de0=_mm256_add_pd(de0, de1);
			//	de0=_mm256_add_pd(de0, de2);
			//	de0=_mm256_add_pd(de0, de3);
				//[num3 num2 num1 num0]
				//[den3 den2 den1 den0]
				//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
				//lo=_mm256_extractf128_pd(r, 0)
				//hi=_mm256_extractf128_pd(r, 1)
				//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
				w0=_mm256_hadd_pd(de0, w0);
				__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
				dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
				__m128i mp=_mm_cvtpd_epi32(dp);
				_mm_store_ss((float*)preds+WG_KC, _mm_castsi128_ps(mp));
#undef  WG_KC
			}
			{
#define WG_KC 2
				__m128i one=_mm_set1_epi32(1);
				__m128i me0=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+0);
				__m128i me1=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+1);
			//	__m128i me2=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+2);
			//	__m128i me3=_mm_load_si128((__m128i*)wg_perrors[WG_KC]+3);
				__m256d w0=_mm256_load_pd(wg_weights[WG_KC]+0*4);
				__m256d w1=_mm256_load_pd(wg_weights[WG_KC]+1*4);
			//	__m256d w2=_mm256_load_pd(wg_weights[WG_KC]+2*4);
			//	__m256d w3=_mm256_load_pd(wg_weights[WG_KC]+3*4);
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNW+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+0), 1));
				me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+1), 1));
			//	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+2), 1));
			//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)(eN+WG_KC*WG_NPREDS)+3), 1));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+0));
				me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+1));
			//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+2));
			//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)(eNNE+WG_KC*WG_NPREDS)+3));
				me0=_mm_add_epi32(me0, one);
				me1=_mm_add_epi32(me1, one);
			//	me2=_mm_add_epi32(me2, one);
			//	me3=_mm_add_epi32(me3, one);
				__m256d de0=_mm256_cvtepi32_pd(me0);
				__m256d de1=_mm256_cvtepi32_pd(me1);
			//	__m256d de2=_mm256_cvtepi32_pd(me2);
			//	__m256d de3=_mm256_cvtepi32_pd(me3);
				w0=_mm256_div_pd(w0, de0);
				w1=_mm256_div_pd(w1, de1);
			//	w2=_mm256_div_pd(w2, de2);
			//	w3=_mm256_div_pd(w3, de3);
				de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+0));
				de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+1));
			//	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+2));
			//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)wg_preds[WG_KC]+3));
				de0=_mm256_mul_pd(de0, w0);
				de1=_mm256_mul_pd(de1, w1);
			//	de2=_mm256_mul_pd(de2, w2);
			//	de3=_mm256_mul_pd(de3, w3);
				w0=_mm256_add_pd(w0, w1);
			//	w0=_mm256_add_pd(w0, w2);
			//	w0=_mm256_add_pd(w0, w3);
				de0=_mm256_add_pd(de0, de1);
			//	de0=_mm256_add_pd(de0, de2);
			//	de0=_mm256_add_pd(de0, de3);
				//[num3 num2 num1 num0]
				//[den3 den2 den1 den0]
				//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
				//lo=_mm256_extractf128_pd(r, 0)
				//hi=_mm256_extractf128_pd(r, 1)
				//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
				w0=_mm256_hadd_pd(de0, w0);
				__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
				dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
				__m128i mp=_mm_cvtpd_epi32(dp);
				_mm_store_ss((float*)preds+WG_KC, _mm_castsi128_ps(mp));
#undef  WG_KC
			}
			{
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i mNE	=_mm_load_si128((__m128i*)NE);
				__m128i mp=_mm_load_si128((__m128i*)preds);
				__m128i vmin=_mm_min_epi32(mN, mW);
				__m128i vmax=_mm_max_epi32(mN, mW);
				vmin=_mm_min_epi32(vmin, mNE);
				vmax=_mm_max_epi32(vmax, mNE);
				mp=_mm_max_epi32(mp, vmin);
				mp=_mm_min_epi32(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);
			}
#else
			int
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i vmin=_mm_min_epi16(mN, mW);
			__m128i vmax=_mm_max_epi16(mN, mW);
			__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
			mp=_mm_max_epi16(mp, vmin);
			mp=_mm_min_epi16(mp, vmax);
			_mm_store_si128((__m128i*)preds, mp);
#endif
			
			//if(ky==1&&kx==39)//
			//if(ky==1&&kx==5)//
			//if(ky==1&&kx==743)//
			//	printf("");
			if(fwd||!streamsize)
			{
				yuv[0]=buffer[idx+1]-128;//g
				yuv[1]=buffer[idx+2]-128;//b
				yuv[2]=buffer[idx+0]-128;//r
#ifdef USE_RCTJ2K
				yuv[1]-=yuv[0];
				yuv[2]-=yuv[0];
				yuv[0]+=(yuv[1]+yuv[2])>>2;
#else
				//Pei09 RCT		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
				yuv[2]-=yuv[0];
				yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
#endif
				curr[0]=yuv[0];
				curr[1]=yuv[1];
				curr[2]=yuv[2];

				errors[0]=yuv[0]-preds[0];
				errors[1]=yuv[1]-preds[1];
				errors[2]=yuv[2]-preds[2];
				curr[4]=errors[0];
				curr[5]=errors[1];
				curr[6]=errors[2];
				errors[0]=errors[0]<<(32-8)>>(32-8);
				errors[1]=errors[1]<<(32-9)>>(32-9);
				errors[2]=errors[2]<<(32-9)>>(32-9);
				errors[0]=errors[0]<<1^errors[0]>>31;//pack sign
				errors[1]=errors[1]<<1^errors[1]>>31;
				errors[2]=errors[2]<<1^errors[2]>>31;

				if(fwd)
				{
					curr_CDF=CDF;
					cdfs[0]=curr_CDF[errors[0]];
					freqs[0]=curr_CDF[errors[0]+1]-cdfs[0];

					curr_CDF=CDF+257;
					cdfs[1]=curr_CDF[errors[1]];
					freqs[1]=curr_CDF[errors[1]+1]-cdfs[1];
			
					curr_CDF=CDF+257+513;
					cdfs[2]=curr_CDF[errors[2]];
					freqs[2]=curr_CDF[errors[2]+1]-cdfs[2];
#ifdef ESTIMATE_SIZE
					csizes[0]-=log2((double)freqs[0]/0x10000);
					csizes[1]-=log2((double)freqs[1]/0x10000);
					csizes[2]-=log2((double)freqs[2]/0x10000);
#endif
					while(range<0x10000)
					{
						*sptr++=(unsigned short)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[0]>>16;
					range=(range*freqs[0]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[0], freqs[0], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

					while(range<0x10000)
					{
						*sptr++=(unsigned short)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[1]>>16;
					range=(range*freqs[1]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[1], freqs[1], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

					while(range<0x10000)
					{
						*sptr++=(unsigned short)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[2]>>16;
					range=(range*freqs[2]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[2], freqs[2], lo0, lo0+r0, low, low+range, 0, 0);//
#endif
				}
			}
			else
			{
				unsigned long long code2;
				unsigned cdf, freq, sym;
				//if(!ky&&kx==513)//
				//	printf("");
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[0]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF+257;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[1]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF+257+513;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[2]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				yuv[0]=errors[0]>>1^-(errors[0]&1);//unpack sign
				yuv[1]=errors[1]>>1^-(errors[1]&1);
				yuv[2]=errors[2]>>1^-(errors[2]&1);
				//if(!ky&&kx&&yuv[0])
				//	printf("");
				yuv[0]+=preds[0];
				yuv[1]+=preds[1];
				yuv[2]+=preds[2];
				yuv[0]=yuv[0]<<(32-8)>>(32-8);
				yuv[1]=yuv[1]<<(32-9)>>(32-9);
				yuv[2]=yuv[2]<<(32-9)>>(32-9);
				curr[0]=yuv[0];
				curr[1]=yuv[1];
				curr[2]=yuv[2];
				curr[4]=yuv[0]-preds[0];
				curr[5]=yuv[1]-preds[1];
				curr[6]=yuv[2]-preds[2];
#ifdef USE_RCTJ2K
				yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[2]+=yuv[0];
				yuv[1]+=yuv[0];
#else
				//inverse Pei09 RCT
				yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
				yuv[2]+=yuv[0];
				yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
#endif
				buffer[idx+1]=(unsigned char)(yuv[0]+128);
				buffer[idx+2]=(unsigned char)(yuv[1]+128);
				buffer[idx+0]=(unsigned char)(yuv[2]+128);
			}
#ifdef USE_WG
			static const int factors[]={97, 99, 99};
#define WG_KC 0
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
			for(int k=0;k<WG_NPREDS;++k)
			{
				int pr=wg_preds[WG_KC][k];
				int e2=abs(curr[WG_KC]-pr)<<1;
				wg_perrors[WG_KC][k]=(wg_perrors[WG_KC][k]+e2)*factors[WG_KC]>>7;
				ecurr[k]=(eW[k]+e2+eN[k]+eNE[k])>>2;
				eNE[k]+=e2;
			}
#undef  WG_KC

#define WG_KC 1
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
			for(int k=0;k<WG_NPREDS;++k)
			{
				int pr=wg_preds[WG_KC][k];
				int e2=abs(curr[WG_KC]-pr)<<1;
				wg_perrors[WG_KC][k]=(wg_perrors[WG_KC][k]+e2)*factors[WG_KC]>>7;
				ecurr[k]=(eW[k]+e2+eN[k]+eNE[k])>>2;
				eNE[k]+=e2;
			}
#undef  WG_KC

#define WG_KC 2
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
			for(int k=0;k<WG_NPREDS;++k)
			{
				int pr=wg_preds[WG_KC][k];
				int e2=abs(curr[WG_KC]-pr)<<1;
				wg_perrors[WG_KC][k]=(wg_perrors[WG_KC][k]+e2)*factors[WG_KC]>>7;
				ecurr[k]=(eW[k]+e2+eN[k]+eNE[k])>>2;
				eNE[k]+=e2;
			}
#undef  WG_KC
#endif
			++hist[errors[0]];
			++hist[256];
			++hist[257+errors[1]];
			++hist[257+512];
			++hist[257+513+errors[2]];
			++hist[257+513+512];
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
		if(fwd)
		{
			*sptr++=(unsigned short)(low>>32);
			*sptr++=(unsigned short)(low>>16);
			*sptr++=(unsigned short)low;

			streamsize=(int)(sptr-sbuf);//number of emits
			if(streamsize<<1>rowsize+4&&0)//bypass row
			{
				streamsize=0;
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(buffer+rowsize*ky, 1, rowsize, fdst);
			}
			else
			{
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(sbuf, sizeof(short), streamsize, fdst);
			}
		}

		update_CDF(hist, 256, CDF);
		update_CDF(hist+257, 512, CDF+257);
		update_CDF(hist+257+513, 512, CDF+257+513);
		rescale_hist(hist, 256);
		rescale_hist(hist+257, 512);
		rescale_hist(hist+257+513, 512);
	}
	if(!fwd)
		fwrite(buffer, 1, usize, fdst);
	fclose(fsrc);
	fclose(fdst);
	_mm_free(pixels);
	_mm_free(ebuf);
	free(hist);
	free(CDF);
	free(sbuf);
	dstsize=get_filesize(dstfn);
#ifdef ESTIMATE_SIZE
	t=time_sec()-t;
	if(fwd)
	{
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		printf("%.2lf + %.2lf + %.2lf = %.2lf  ->  %td\n",
			csizes[0],
			csizes[1],
			csizes[2],
			csizes[0]+csizes[1]+csizes[2],
			dstsize
		);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
#ifdef BYPASS_ON_INFLATION
	if(fwd&&dstsize>usize)
	{
		printf("Bypass encode\n");
		fdst=fopen(dstfn, "wb");
		fwrite("BP", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
		fwrite(buffer, 1, usize, fdst);
		fclose(fdst);
	}
#endif
	free(buffer);
	(void)nwritten;
	return 0;
}