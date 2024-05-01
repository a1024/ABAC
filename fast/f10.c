#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;

//#include"ac.h"
//#include"profiler.h"

//	#define REUSE_HIST

static const char *ext[]=
{
	"PPM",
	"PNG",
	"BMP",
	"JPG",
	"JPEG",
};
int f10_mptest(const char *path)
{
	printf("F10 video test\n");
	double t0=time_sec();
	ArrayHandle filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames||!filenames->count)
	{
		LOG_ERROR("No media in \"%s\"", path);
		return 0;
	}
	double elapsed=0;
	size_t usize_total=0;
	double csizes[4]={0};
	Image prev={0}, curr={0};
	short *pixels=0;
#ifdef REUSE_HIST
	int *hist=0;//TODO size_t
#endif
	for(int kf=0;kf<(int)filenames->count;++kf)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, kf);
		image_load((char*)fn[0]->data, &curr);
		if(!curr.data)
		{
			LOG_ERROR("Cannot load \"%s\"", (char*)fn[0]->data);
			return 0;
		}
		if(prev.data&&(curr.depth!=prev.depth||curr.nch!=prev.nch||curr.iw!=prev.iw||curr.ih!=prev.ih))
		{
			LOG_ERROR("Change at frame %d: CWHD %d*%d*%d*%d -> %d*%d*%d*%d", kf,
				prev.nch, prev.iw, prev.ih, prev.depth,
				curr.nch, curr.iw, curr.ih, curr.depth
			);
			return 0;
		}
		int nlevels=1<<curr.depth, half=nlevels>>1;
		if(!pixels)
		{
			size_t bufsize=(curr.iw+2LL)*sizeof(short[2*2*4]);
			pixels=(short*)malloc(bufsize);
			if(!pixels)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(pixels, 0, bufsize);
			usize_total=((size_t)curr.iw*curr.ih*curr.nch*filenames->count*curr.depth)>>3;
		}
#ifdef REUSE_HIST
		if(!hist)
		{
			size_t histsize=sizeof(int[4])<<curr.depth;
			hist=(int*)malloc(histsize);
			if(!hist)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(hist, 0, histsize);
		}
#else
		size_t histsize=sizeof(int[4])<<curr.depth;
		int *hist=(int*)malloc(histsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(hist, 0, histsize);
#endif
		double t=time_sec();
		for(int ky=0, idx=0;ky<curr.ih;++ky)
		{
			short *rows[]=
			{
				pixels+(((curr.iw+2LL)*(((kf-0LL)&1)<<1|((ky-0LL)&1))+1LL)<<2),
				pixels+(((curr.iw+2LL)*(((kf-0LL)&1)<<1|((ky-1LL)&1))+1LL)<<2),
				pixels+(((curr.iw+2LL)*(((kf-1LL)&1)<<1|((ky-0LL)&1))+1LL)<<2),
				pixels+(((curr.iw+2LL)*(((kf-1LL)&1)<<1|((ky-1LL)&1))+1LL)<<2),
			};
			for(int kx=0;kx<curr.iw;++kx, idx+=curr.nch)
			{
				for(int kc=0;kc<curr.nch;++kc)
				{
					int
						NWp1	=rows[3][kc-4],
						Np1	=rows[3][kc+0],
						Wp1	=rows[2][kc-4],
						p1	=rows[2][kc+0],
						NW	=rows[1][kc-4],
						N	=rows[1][kc+0],
						W	=rows[0][kc-4],
						offset	=0;
					if(kc)
						offset+=rows[0][0];
					int pred=N+W+p1-(NW+Np1+Wp1)+NWp1;
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					UPDATE_MIN(vmin, p1);
					UPDATE_MAX(vmax, p1);
					pred=CLAMP(vmin, pred, vmax);
					pred+=offset;
					pred=CLAMP(-half, pred, half-1);

					pred=curr.data[idx+kc]-pred;
					pred+=half;
					pred&=nlevels-1;
					//pred-=half;
					++hist[kc<<curr.depth|pred];
					rows[0][kc]=(short)(curr.data[idx+kc]-offset);
				}
				for(int k=0;k<4;++k)
					rows[k]+=4;
			}
		}
#ifdef REUSE_HIST
		int res=curr.iw*curr.ih*(kf+1);
#else
		int res=curr.iw*curr.ih;
#endif
		for(int kc=0;kc<curr.nch;++kc)
		{
			double csize=0;
			for(int ks=0;ks<nlevels;++ks)
			{
				int freq=hist[kc<<curr.depth|ks];
				if(freq)
					csize-=freq*log2((double)freq/res);
			}
#ifdef REUSE_HIST
			csizes[kc]=csize;
#else
			csizes[kc]+=csize;
#endif
		}
		elapsed+=time_sec()-t;
		size_t usizesofar=((size_t)curr.iw*curr.ih*curr.nch*(kf+1LL)*curr.depth)>>3;
		double csizesofar=(csizes[0]+csizes[1]+csizes[2]+csizes[3])/8;
		printf("Progress %7d/%7d = %6.2lf%%  %16.2lf/%14zd = %10.6lf  CR %16.6lf\r",
			kf+1, (int)filenames->count, 100.*(kf+1)/filenames->count,
			csizesofar, usizesofar, csizesofar/usizesofar, usizesofar/csizesofar
		);
#ifndef REUSE_HIST
		free(hist);
#endif
		image_clear(&prev);
		prev=curr;
	}
#ifdef REUSE_HIST
	free(hist);
#endif
	image_clear(&prev);
	double csize_total=0;
	t0=time_sec()-t0;
	printf("\nYUVA ");
	for(int kc=0;kc<4;++kc)
	{
		double csize=csizes[kc]/8;
		printf(" %14.2lf", csize);
		csize_total+=csize;
	}
	printf("\n");
	printf("Total %16.2lf/%14zd = %10.6lf  CR %16.6lf\n", csize_total, usize_total, csize_total/usize_total, usize_total/csize_total);
	printf("Loop    %14lf sec\n", elapsed);
	printf("Elapsed %14lf sec\n", t0);
	free(pixels);

	LOG_ERROR("This isn't a codec.");
	return 1;
}