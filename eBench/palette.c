#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


static ArrayHandle levels[4]={0};
void pred_palette(Image *src, int fwd)
{
	ptrdiff_t res;
	int maxdepth, histsize, *hist;

	maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	histsize=sizeof(int)<<maxdepth;
	hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	res=(ptrdiff_t)src->iw*src->ih;
	for(int kc=0;kc<4;++kc)
	{
		int *vals, nvals;
		int depth, nlevels, half, mask;

		depth=src->depth[kc];
		if(!depth)
			continue;
		nlevels=1<<depth;
		half=nlevels>>1;
		mask=nlevels-1;
		memset(hist, 0, histsize);
		if(fwd)
		{
			array_append(levels+kc, 0, sizeof(int), 1, 1, 0, 0);
			array_clear(levels+kc);
			for(ptrdiff_t k=0;k<res;++k)
			{
				int val=src->data[k<<2|kc];
				val+=half;
				val&=mask;
				++hist[val];
			}
			for(int ks=0, ks2=0;ks<nlevels;++ks)
			{
				int freq=hist[ks];
				hist[ks]=ks2;
				if(freq)
				{
					array_append(levels+kc, &ks, sizeof(int), 1, 1, 0, 0);
					++ks2;
				}
			}
			vals=(int*)levels[kc]->data;
			nvals=(int)levels[kc]->count;
			for(ptrdiff_t k=0;k<res;++k)
			{
				int val=src->data[k<<2|kc];
				val+=half;
				val&=mask;
				src->data[k<<2|kc]=hist[val]-(nvals>>1);
			}
		}
		else
		{
			if(!levels[kc])
				continue;
			vals=(int*)levels[kc]->data;
			nvals=(int)levels[kc]->count;
			for(ptrdiff_t k=0;k<res;++k)
			{
				int idx=src->data[k<<2|kc];
				idx+=nvals>>1;
				idx%=nvals;
				src->data[k<<2|kc]=vals[idx]+half;
			}
		}
	}
}