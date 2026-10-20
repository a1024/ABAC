#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#if defined _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif


static ArrayHandle levels[4]={0};
void pred_palette(Image *src, int fwd)
{
	intptr_t k, res;
	int maxdepth, histsize, *hist, kc, ks, ks2;

	maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	histsize=sizeof(int)<<maxdepth;
	hist=(int*)malloc(histsize);
	if(!hist)
	{
		CRASH("Alloc error");
		return;
	}
	res=(intptr_t)src->iw*src->ih;
	for(kc=0;kc<4;++kc)
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
			for(k=0;k<res;++k)
			{
				int val=((int32_t*)HDR_DATA(src))[k<<2|kc];
				val+=half;
				val&=mask;
				++hist[val];
			}
			for(ks=0, ks2=0;ks<nlevels;++ks)
			{
				int freq=hist[ks];
				hist[ks]=ks2;
				if(freq)
				{
					array_append(levels+kc, &ks, sizeof(int), 1, 1, 0, 0);
					++ks2;
				}
			}
			vals=(int*)HDR_DATA(levels[kc]);
			nvals=(int)levels[kc]->count;
			for(k=0;k<res;++k)
			{
				int val=((int32_t*)HDR_DATA(src))[k<<2|kc];
				val+=half;
				val&=mask;
				((int32_t*)HDR_DATA(src))[k<<2|kc]=hist[val]-(nvals>>1);
			}
		}
		else
		{
			if(!levels[kc])
				continue;
			vals=(int*)HDR_DATA(levels[kc]);
			nvals=(int)levels[kc]->count;
			for(k=0;k<res;++k)
			{
				int idx=((int32_t*)HDR_DATA(src))[k<<2|kc];
				idx+=nvals>>1;
				idx%=nvals;
				((int32_t*)HDR_DATA(src))[k<<2|kc]=vals[idx]-half;
			}
		}
	}
}