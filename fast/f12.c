#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

	#define BITCTR

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
int f12_statstest(const char *path)
{
	printf("F12 stats\n");
	double t0=time_sec();
	ArrayHandle filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames||!filenames->count)
	{
		LOG_ERROR("No media in \"%s\"", path);
		return 0;
	}
	size_t *hist=(size_t*)malloc(sizeof(size_t[4<<HISTBITS]));
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
		return 0;
	}
	memset(hist, 0, sizeof(size_t[4<<HISTBITS]));
#ifdef BITCTR
	memset(ctr, 0, sizeof(size_t[4*2<<HISTBITS]));
	memset(hist2, 0, sizeof(double[4<<HISTBITS]));
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
			return 0;
		}
		int depths[4]={0};
		memfill(depths, &image.depth, sizeof(depths), sizeof(image.depth));
		if(image.nch>=3)
		{
			rct_JPEG2000_32(&image, 1);
			++depths[1];
			++depths[2];
		}
		char cdepths[]=
		{
			depths[0],
			depths[1],
			depths[2],
			depths[3],
		};
		pred_simd(&image, 1, cdepths);
		//int nlevels[]=
		//{
		//	1<<depths[0],
		//	1<<depths[1],
		//	1<<depths[2],
		//	1<<depths[3],
		//};
		//int halfs[]=
		//{
		//	nlevels[0]>>1,
		//	nlevels[1]>>1,
		//	nlevels[2]>>1,
		//	nlevels[3]>>1,
		//};
		int rowstride=image.iw*image.nch;
		for(int ky=0, idx=0;ky<image.ih;++ky)
		{
			for(int kx=0;kx<image.iw;++kx, idx+=image.nch)
			{
				for(int kc=0;kc<image.nch;++kc)
				{
#define LOAD(X, Y) (((unsigned)(ky+(Y))<(unsigned)image.ih&&(unsigned)(kx+(X))<(unsigned)image.iw?image.data[(image.iw*(ky+(Y))+kx+(X))*image.nch+kc]>>(depths[kc]-HISTBITS):0)+HISTHALF)&HISTMASK
					int
						NW	=LOAD(-1, -1),
						N	=LOAD( 0, -1),
						NE	=LOAD( 1, -1),
						W	=LOAD(-1,  0),
						curr	=LOAD( 0,  0);

					//if(kx==100&&ky==100)//
					//	printf("");

					//int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					//if(BETWEEN_EXC(64, vmin, 65)&&BETWEEN_EXC(192, vmax, 193))
					{
						++hist[kc<<HISTBITS|curr];
						++ctr_hit[kc];
					}
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
		image_clear(&image);
	}
	printf("\n");
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
		printf("C%d  %lld/%lld = %8.4lf%%\n", kc, ctr_hit[kc], ctr_total[kc], 100.*ctr_hit[kc]/ctr_total[kc]);
		if(!ctr_hit[kc])
			continue;
		for(int ks=0;ks<HISTSIZE;++ks)
		{
			size_t freq=hist[kc<<HISTBITS|ks];
			printf("%3d  %8lld  %8.4lf%%", ks, freq, 100.*freq/ctr_hit[kc]);
#ifdef BITCTR
			printf("  %8.4lf%%", 100.*hist2[kc<<HISTBITS|ks]);
#endif
			printf("\n");
		}
		break;
	}
	printf("\nDone.\n");
#ifdef BITCTR
	free(ctr);
	free(hist2);
#endif
	free(hist);
	array_free(&filenames);

	LOG_ERROR("This isn't a codec.");
	return 0;
}