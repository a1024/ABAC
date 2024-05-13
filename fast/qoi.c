#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

#define QOI_IMPLEMENTATION
#include"qoi.h"

void qoi_test(Image const *src, size_t *csize2, double *enc, double *dec, int *error, int loud)
{
	unsigned char *buf=image_export8(src);
	if(!buf)
	{
		LOG_ERROR("Alloc error");
		return;
	}

	qoi_desc desc={(unsigned)src->iw, (unsigned)src->ih, (unsigned char)src->nch, QOI_SRGB};
	int csize=0;
	volatile double t0=time_sec();
	void *data=qoi_encode(buf, &desc, &csize);
	t0=time_sec()-t0;
	if(enc)
		*enc=t0;
	if(csize2)
		*csize2=csize;
	
	ptrdiff_t usize=(ptrdiff_t)src->iw*src->ih*src->nch*src->depth>>3;
	if(loud)
	{
		printf("%14d/%14zd  %8.4lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		printf("E  %16.6lf  %lf MB/s\n", t0, usize/(t0*1024*1024));
	}
	
	t0=time_sec();
	unsigned char *im2=(unsigned char*)qoi_decode(data, csize, &desc, src->nch);
	t0=time_sec()-t0;
	if(dec)
		*dec=t0;
	
	if(loud)
		printf("D  %16.6lf  %lf MB/s\n", t0, usize/(t0*1024*1024));
	int e2=compare_bufs_8(im2, buf, src->iw, src->ih, src->nch, src->nch, "QOI", 0, loud);
	if(error)
		*error=e2;

	free(buf);
	free(data);
	free(im2);
}