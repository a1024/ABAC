#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"lodepng.h"
static const char file[]=__FILE__;

//	#define RELATION
	#define DISABLE_DECORRELATION
//	#define BITCTR
//	#define ABACHIST	//X

#define HISTBITS 8
#define HISTSIZE (1<<HISTBITS)
#define HISTHALF (1<<HISTBITS>>1)
#define HISTMASK ((1<<HISTBITS)-1)

static const char *ext[]=
{
	"PPM",
	"PNG",
	"BMP",
	"JPG",
	"JPEG",
};
static void hist_snapshot(size_t *hist, unsigned char *result, int idx)
{
	//static const char *rnames[]=
	//{
	//	"xW-yN",
	//	"xW-yNW",
	//	"xN-yNE",
	//	"xNW-yNE",
	//};
	//for(int ks=0;ks<HISTSIZE;++ks)
	//	printf("%3d  %8zd\n", ks, hist[((0<<HISTBITS|HISTHALF)<<HISTBITS|HISTHALF)<<HISTBITS|ks]);
	for(int kc=0;kc<4;++kc)
	{
		int x=(kc&1)<<HISTBITS, y=(kc>>1&1)<<HISTBITS;
		for(int kN=0;kN<HISTSIZE;++kN)
		{
			for(int kW=0;kW<HISTSIZE;++kW)
			{
				double mean=0, wsum=0;
				for(int ks=0;ks<HISTSIZE;++ks)
				{
					double freq=(double)hist[((kc<<HISTBITS|kN)<<HISTBITS|kW)<<HISTBITS|ks];
					mean+=freq*ks;
					wsum+=freq;
				}
				if(wsum)
					mean/=wsum;
				else
					mean=HISTHALF;
				int kmean=(int)round(255*mean/(HISTSIZE-1));
#if 1
				kmean-=(kN+kW)>>1;
				kmean+=128;
				kmean&=255;
				if(!wsum)
					kmean=128;
#endif
				result[(y|kN)<<(HISTBITS+1)|x|kW]=(unsigned char)kmean;
				//result[kN<<HISTBITS|kW]=(int)round(255*mean/(HISTSIZE-1));
				//printf(" %X", (int)round(15*mean/(HISTSIZE-1)));
			}
			//printf("\n");
		}
	}
	snprintf(g_buf, G_BUF_SIZE, "sample_%03d.PNG", idx+1);
	//snprintf(g_buf, G_BUF_SIZE, "%s_%03d.PNG", rnames[kc], idx+1);
	lodepng_encode_file(g_buf, result, HISTSIZE<<1, HISTSIZE<<1, LCT_GREY, 8);
}
int f12_statstest(const char *path)
{
	printf("F12 stats\n");
	ArrayHandle filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames||!filenames->count)
	{
		LOG_ERROR("No media in \"%s\"", path);
		return 1;
	}
#ifdef RELATION
	size_t histsize=sizeof(size_t[4<<HISTBITS*3]);
#else
	size_t histsize=sizeof(size_t[4<<HISTBITS]);
#endif
	size_t *hist=(size_t*)malloc(histsize);
#ifdef BITCTR
	size_t *ctr=(size_t*)malloc(sizeof(size_t[4*2<<HISTBITS]));
	double *hist2=(double*)malloc(sizeof(double[4<<HISTBITS]));
#endif
	if(!hist
#ifdef BITCTR
		||!ctr||!hist2
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(hist, 0, histsize);
#ifdef BITCTR
	memset(ctr, 0, sizeof(size_t[4*2<<HISTBITS]));
	memset(hist2, 0, sizeof(double[4<<HISTBITS]));
#endif
#ifdef RELATION
	size_t rsize=sizeof(char[4*HISTSIZE*HISTSIZE]);
	unsigned char *result=(unsigned char*)malloc(rsize);
	if(!result)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(result, 0, rsize);
#endif
#ifdef ABACHIST
	size_t abachistsize=sizeof(size_t[HISTSIZE]);
	size_t *abachist=(size_t*)malloc(abachistsize);
	if(!abachist)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(abachist, 0, abachistsize);
#endif

	Image image={0};
	size_t ctr_total[4]={0}, ctr_hit[4]={0};
	int maxlen=0;
	for(int kf=0;kf<(int)filenames->count;++kf)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, kf);
		int len=(int)fn[0]->count;
		UPDATE_MAX(maxlen, len);
	}
	for(int kf=0;kf<(int)filenames->count;++kf)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, kf);
		printf("%5d/%5d  %-*s\r", kf+1, (int)filenames->count, maxlen, (char*)fn[0]->data);
		image_load((char*)fn[0]->data, &image);
		if(!image.data)
		{
			LOG_ERROR("Cannot load \"%s\"", (char*)fn[0]->data);
			return 1;
		}
		int depths[4]={0};
		memfill(depths, &image.depth, sizeof(depths), sizeof(image.depth));
#ifndef DISABLE_DECORRELATION
		if(image.nch>=3)
		{
			rct_JPEG2000_32(&image, 1);
			++depths[1];
			++depths[2];
		}
#endif
		char cdepths[]=
		{
			(char)depths[0],
			(char)depths[1],
			(char)depths[2],
			(char)depths[3],
		};
#ifndef DISABLE_DECORRELATION
		pred_simd(&image, 1, cdepths);
#endif
		for(int ky=0, idx=0;ky<image.ih;++ky)
		{
			for(int kx=0;kx<image.iw;++kx, idx+=image.nch)
			{
				int kc=1;
				//for(int kc=0;kc<image.nch;++kc)
				{
#define LOAD(X, Y) (((unsigned)(ky+(Y))<(unsigned)image.ih&&(unsigned)(kx+(X))<(unsigned)image.iw?image.data[(image.iw*(ky+(Y))+kx+(X))*image.nch+kc]>>(depths[kc]-HISTBITS):0)+HISTHALF)&HISTMASK
					int
						NW	=LOAD(-1, -1),
						N	=LOAD( 0, -1),
						NE	=LOAD( 1, -1),
						W	=LOAD(-1,  0),
						curr	=LOAD( 0,  0);
					(void)NW;
					(void)N;
					(void)NE;
					(void)W;
					(void)curr;
					if(!kx)
						NW=N, W=N;
					if(!ky)
						NW=W, N=W;
					if(kx==image.iw-1)
						NE=N;
					
#ifdef ABACHIST
					if(kc==1)
					{
						int low=0, range=HISTSIZE;
						for(int kb=HISTBITS-1;kb>=0;--kb)
						{
							int bit=curr>>kb&1;
						
							int floorhalf=range>>1;
							if(bit)
								low+=floorhalf;
							range=floorhalf;
							for(int k2=0;k2<range;++k2)	//X  should increment once at the end, leading to same multisymbol histogram
								++abachist[low+k2];
						}
					}
#endif
#ifdef RELATION
					int cgrad=N+W-NW;
					int vmin=W, vmax=N;
					if(N<W)
						vmin=N, vmax=W;
					cgrad=CLAMP(vmin, cgrad, vmax);
					int av2=(N+W)>>1;
					//++hist[((kc<<HISTBITS|N)<<HISTBITS|W)<<HISTBITS|curr];

					//X-Y:
					//- |	- \
					//| /	\ /	loss?

					//W-N	W-NW
					//N-NE	NW-NE
					++hist[((0<<HISTBITS|N )<<HISTBITS|W )<<HISTBITS|curr];//XY 0 0    X-Y W-N
					++hist[((1<<HISTBITS|av2)<<HISTBITS|cgrad)<<HISTBITS|curr];
				//	++hist[((1<<HISTBITS|NW)<<HISTBITS|W )<<HISTBITS|curr];//XY 1 0    X-Y W-NW
					++hist[((2<<HISTBITS|NE)<<HISTBITS|N )<<HISTBITS|curr];//XY 0 1    X-Y N-NE
					++hist[((3<<HISTBITS|NE)<<HISTBITS|NW)<<HISTBITS|curr];//XY 1 1    X-Y NW-NE
#else
					++hist[kc<<HISTBITS|curr];
#endif
					++ctr_hit[kc];
					++ctr_total[kc];
#ifdef BITCTR
					size_t *currctr=ctr+((size_t)kc<<(HISTBITS+1));
					int idx2=1;
					for(int kb=HISTBITS-1;kb>=0;--kb)
					{
						int bit=curr>>kb&1;
						idx2=idx2<<1|bit;
						++currctr[idx2];
					}
#endif
				}
			}
		}
#ifdef RELATION
		hist_snapshot(hist, result, kf);
		memset(hist, 0, histsize);
#endif
		image_clear(&image);
	}
	printf("\n");

#ifdef RELATION
	//hist_snapshot(hist, result, 0);
#endif

#ifndef RELATION
#ifdef BITCTR
	for(int kc=0;kc<4;++kc)
	{
		size_t *currctr=ctr+((size_t)kc<<(HISTBITS+1));
		if(!ctr_total[kc])
			continue;
		for(int ks=0;ks<HISTSIZE;++ks)
		{
			double prob=1;
			int idx2=1;
			for(int kb=HISTBITS-1;kb>=0;--kb)
			{
				int bit=ks>>kb&1;
				idx2<<=1;
				double sum=(double)(currctr[idx2|0]+currctr[idx2|1]);
				idx2|=bit;
				if(sum)
					prob*=currctr[idx2]/sum;
			}
			hist2[kc<<HISTBITS|ks]=prob;
		}
	}
#endif
	printf("ks  freq  freq%%  binctr->freq%%\n");
	for(int kc=0;kc<4;++kc)
	{
		if(!ctr_total[kc])
			continue;
		printf("C%d  %zd/%zd = %8.4lf%%\n", kc, ctr_hit[kc], ctr_total[kc], 100.*ctr_hit[kc]/ctr_total[kc]);
		if(!ctr_hit[kc])
			continue;
		for(int ks=0;ks<HISTSIZE;++ks)
		{
			size_t freq=hist[kc<<HISTBITS|ks];
			printf("%3d  %8zd  %8.4lf%%", ks, freq, 100.*freq/ctr_hit[kc]);
#ifdef ABACHIST
			if(kc==1)
			{
				freq=abachist[ks];
				printf("  %16lld%%", freq);
			}
#endif
#ifdef BITCTR
			printf("  %8.4lf%%", 100.*hist2[kc<<HISTBITS|ks]);
#endif
			printf("\n");
		}
		break;
	}
#endif
	printf("\nDone.\n");
#ifdef RELATION
	free(result);
#endif
#ifdef BITCTR
	free(ctr);
	free(hist2);
#endif
	free(hist);
	array_free(&filenames);

	LOG_ERROR("This isn't a codec.");
	return 0;
}