#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>

/*
"A Block-sorting Lossless Data Compression Algorithm" by M. Burrows and D.J. Wheeler
https://github.com/MichaelDipperstein/bwt
*/

#define WRAP(VALUE, LIMIT) ((VALUE)-((LIMIT)&-((VALUE)>=(LIMIT))))	/*wraps array index within array bounds (assumes value < 2 * limit)*/
typedef struct _BWTSortInfo
{
	const uint8_t *ptr;
	int count;
} BWTSortInfo;
/*
https://stackoverflow.com/questions/4210689/pass-extra-parameter-to-comparator-for-qsort
*/
#if defined _MSC_VER || __STDC__ < 199012
static int bwt_cmp_sorted(void *param, const void *s1, const void *s2)/*for qsort_s (MSVC)*/
#elif defined __GNUC__
static int bwt_cmp_sorted(const void *s1, const void *s2, void *param)/*for qsort_r (GNU C, C11)*/
#endif
{
	const BWTSortInfo *ctx=(const BWTSortInfo*)param;
	int offset1, offset2;
	int i;

	/*
	Compare 1 character at a time until there's difference or the end of the block is reached.
	Since we're only sorting strings that already match at the first two characters, start with the third character.
	*/
	offset1=*(const int*)s1+2;
	offset2=*(const int*)s2+2;
	for(i=2;i<ctx->count;++i)
	{
		uint8_t c1, c2;

		if(offset1>=ctx->count)offset1-=ctx->count;/*ensure that offsets are properly bounded*/
		if(offset2>=ctx->count)offset2-=ctx->count;

		c1=ctx->ptr[offset1];
		c2=ctx->ptr[offset2];

		if (c1>c2)
			return 1;
		if (c2>c1)
			return -1;

		offset1++;/*strings match to here, try next character*/
		offset2++;
	}
	return 0;/*strings are identical*/
}

/*
Forward BWT:
src: byte[count]
hist:		int[nlevels]
CDF:		int[nlevels]
rotationidx:	int[count]
buf:		int[count]
dst: byte[count]
returns: s0idx for reversibility
*/
static int bwt_fwd(uint8_t *ptr, int count, int nlevels, int *hist, int *CDF, int *rotationidx, int *buf, uint8_t *dst)
{
	BWTSortInfo ctx={0};
	int i, j, k, sym, sum, s0idx;

	ctx.ptr=ptr;
	ctx.count=count;
/*
Sort the rotated strings in the block.
A radix sort is performed on the first to characters of all the rotated strings (2nd character then 1st).
All rotated strings with matching initial characters are then quicksorted. - Q4..Q7
*/

	/*radix sort by second character in rotation*/
	memset(hist, 0, sizeof(int)*nlevels);
	for(k=0;k<count;++k)/*count number of characters for radix sort*/
		++hist[ptr[k]];
	for(sym=0, sum=0;sym<nlevels;++sym)/*accumulate CDF*/
	{
		int freq=hist[sym];
		CDF[sym]=sum;
		sum+=freq;
	}

	for(k=0;k<count-1;++k)/*sort by 2nd character*/
	{
		j=ptr[k+1];
		buf[CDF[j]]=k;
		++CDF[j];
	}
	j=ptr[0];/*handle wrap around for string starting at end of block*/
	buf[CDF[j]]=k;
	CDF[0]=0;

	/*radix sort by first character in rotation*/
	for(sym=0, sum=0;sym<nlevels;++sym)/*accumulate CDF*/
	{
		int freq=hist[sym];
		CDF[sym]=sum;
		sum+=freq;
	}
	for(k=0;k<count;++k)
	{
		j=ptr[buf[k]];
		rotationidx[CDF[j]]=buf[k];
		++CDF[j];
	}

	/*
	now rotationIdx contains the sort order of all strings sorted by their first 2 characters.
	Use qsort to sort the strings that have their first two characters matching.
	*/
	for(i=0, k=0;i<nlevels&&k<count-1;++i)
	{
		for(j=0;j<nlevels&&k<count-1;++j)
		{
			int first=k;
			while(i==ptr[rotationidx[k]]&&j==ptr[WRAP(rotationidx[k]+1, count)])
			{
				++k;
				if(k==count)/*we've searched the whole block*/
					break;
			}
			if(k-first>1)/*there are at least 2 strings staring with ij, sort them*/
				qsort_s(rotationidx+first, (size_t)k-first, sizeof(int), bwt_cmp_sorted, &ctx);
		}
	}

	/*find last characters of rotations (L) - C2*/
	s0idx=0;
	for(i=0;i<count;++i)
	{
		if(rotationidx[i])
			dst[i]=ptr[rotationidx[i]-1];
		else
		{
			s0idx=i;/*unrotated string 1st character is end of string*/
			dst[i]=ptr[count-1];
		}
	}
	return s0idx;
}

/*
Inverse BWT:
ptr:	byte[count]
s0idx:	value returned by bwt_fwd
hist:	int[nlevels]
pred:	int[count]
dst:	byte[count]
*/
static void bwt_inv(uint8_t *ptr, int s0idx, int count, int nlevels, int *hist, int *pred, uint8_t *dst)
{
	int i, j, k, sum;

	/*based on pseudo code from section 4.2 (D1 and D2)*/
	memset(hist, 0, sizeof(int)*nlevels);
	for(k=0;k<count;++k)
	{
		pred[k]=hist[ptr[k]];		/*Set pred[i] to the number of times block[i] appears in the substring block[0 .. i - 1].*/
		++hist[ptr[k]];			/*As a useful side effect count[i] will be the number of times character i appears in block.*/
	}
	for(k=0, sum=0;k<nlevels;++k)	/*Finally, set count[i] to the number of characters in block lexicographically less than i.*/
	{
		int freq=hist[k];
		hist[k]=sum;
		sum+=freq;
	}
	for(i=s0idx, j=count;j>0;--j)	/*construct the initial unrotated string (S[0])*/
	{
		dst[j-1]=ptr[i];
		i=pred[i]+hist[ptr[i]];
	}
}

void prep_BWT_x(Image **psrc, int fwd)
{
	Image *src=*psrc;
	int maxdepth, maxlevels, kc, ky, kx;

	maxdepth=calc_maxdepth(src, 0);
	maxlevels=1<<maxdepth;
#if 0
	if(src->depth[0]!=8)
	{
		CRASH("BWT: Expected 8-bit depth");
		return;
	}
#endif
	if(src->iw>0x10000||src->ih>0x10000)
	{
		CRASH("BWT: Image dimensions must be up to 65536");
		return;
	}
	if(fwd)
	{
		int h2=src->ih, w2=src->iw+2, size2=h2*w2;
		uint8_t *sbuf=(uint8_t*)malloc(w2);
		uint8_t *dbuf=(uint8_t*)malloc(w2);
#if 0
		uint8_t *debugbuf=(uint8_t*)malloc(w2);
#endif
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*maxlevels);
		int *t3=(int*)malloc(sizeof(int)*w2);
		int *t4=(int*)malloc(sizeof(int)*w2);
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));/*don't access Image *src until it is reassigned!*/
			if(!ptr||!sbuf||!dbuf||!t1||!t2||!t3||!t4)
			{
				CRASH("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}
		for(kc=0;kc<src->nch;++kc)
		{
			for(ky=src->ih-1;ky>=0;--ky)
			{
				int bwtidx;
				int *ptr=(int32_t*)HDR_DATA(src)+((size_t)src->iw*ky<<2|kc);
				for(kx=0;kx<src->iw;++kx)
					sbuf[kx]=(ptr[kx<<2]+128)&255;

				bwtidx=bwt_fwd(sbuf, src->iw, 1<<src->depth[kc], t1, t2, t3, t4, dbuf);

#if 0
				bwt_inv(dbuf, bwtidx, src->iw, 1<<src->depth[kc], t1, t3, debugbuf);/**/
				if(memcmp(debugbuf, sbuf, src->iw))
					CRASH("BWT Error ky=%d", ky);
#endif

				ptr=(int32_t*)HDR_DATA(src)+((size_t)w2*ky<<2|kc);
				for(kx=0;kx<src->iw;++kx)
					ptr[kx<<2]=dbuf[kx]-128;
				ptr[(src->iw+(intptr_t)0)<<2]=(bwtidx>>0&255)-128;/*little-endian*/
				ptr[(src->iw+(intptr_t)1)<<2]=(bwtidx>>8&255)-128;
			}
		}
		src->ih=h2;
		src->iw=w2;

		free(sbuf);
		free(dbuf);
		free(t1);
		free(t2);
		free(t3);
		free(t4);
	/*	free(debugbuf);*/
	}
	else
	{
		int h2=src->ih, w2=src->iw-2, size2=h2*w2;
		uint8_t *sbuf=(uint8_t*)malloc(w2);
		uint8_t *dbuf=(uint8_t*)malloc(w2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*w2);
		if(!sbuf||!dbuf||!t1||!t2)
		{
			CRASH("Alloc error");
			return;
		}

		for(kc=0;kc<src->nch;++kc)
		{
			for(ky=0;ky<src->ih;++ky)
			{
				int *ptr=(int32_t*)HDR_DATA(src)+((size_t)src->iw*ky<<2|kc);
				int bwtidx=
					((ptr[(w2+(intptr_t)0)<<2]+128)&255)<<0|/*little-endian*/
					((ptr[(w2+(intptr_t)1)<<2]+128)&255)<<8;
				for(kx=0;kx<w2;++kx)
					sbuf[kx]=(ptr[kx<<2]+128)&255;

#if 0
				if((unsigned)bwtidx>=(unsigned)w2)
				{
					CRASH("Don't click Inv BWT without Fwd BWT!");
					return;
				}
#endif
				bwtidx%=w2;
				bwt_inv(sbuf, bwtidx, w2, 1<<src->depth[kc], t1, t2, dbuf);
				
				ptr=(int32_t*)HDR_DATA(src)+((size_t)w2*ky<<2|kc);
				for(kx=0;kx<w2;++kx)
					ptr[kx<<2]=dbuf[kx]-128;
			}
		}
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));/*don't access Image *src until it is reassigned!*/
			if(!ptr)
			{
				CRASH("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}
		src->ih=h2;
		src->iw=w2;

		free(sbuf);
		free(dbuf);
		free(t1);
		free(t2);
	}
	*psrc=src;
}
void prep_BWT_y(Image **psrc, int fwd)
{
	Image *src=*psrc;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	Image *tempimg=0;
	int kc, kx, ky;

	if(src->iw>0x10000||src->ih>0x10000)
	{
		CRASH("BWT: Image dimensions must be up to 65536");
		return;
	}
	image_copy(&tempimg, src);
	if(!tempimg)
	{
		CRASH("Alloc error");
		return;
	}
	if(fwd)
	{
		int h2=src->ih+2, w2=src->iw, size2=h2*w2;
		uint8_t *sbuf=(uint8_t*)malloc(h2);
		uint8_t *dbuf=(uint8_t*)malloc(h2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*maxlevels);
		int *t3=(int*)malloc(sizeof(int)*h2);
		int *t4=(int*)malloc(sizeof(int)*h2);
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));/*don't access Image *src until it is reassigned!*/
			if(!ptr||!sbuf||!dbuf||!t1||!t2||!t3||!t4)
			{
				CRASH("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}

		for(kc=0;kc<src->nch;++kc)
		{
			for(kx=0;kx<src->iw;++kx)
			{
				int bwtidx;
				int *ptr=(int32_t*)HDR_DATA(tempimg)+((size_t)kx<<2|kc);
				for(ky=0;ky<src->ih;++ky)
					sbuf[ky]=(ptr[(size_t)src->iw*ky<<2]+128)&255;

				bwtidx=bwt_fwd(sbuf, src->ih, 1<<src->depth[kc], t1, t2, t3, t4, dbuf);

				ptr=(int32_t*)HDR_DATA(src)+((size_t)kx<<2|kc);
				for(ky=0;ky<src->ih;++ky)
					ptr[(size_t)w2*ky<<2]=dbuf[ky]-128;
				ptr[((size_t)w2*(src->ih+(intptr_t)0))<<2]=(bwtidx>>0&255)-128;/*little-endian*/
				ptr[((size_t)w2*(src->ih+(intptr_t)1))<<2]=(bwtidx>>8&255)-128;
			}
		}
		src->ih=h2;
		src->iw=w2;

		free(sbuf);
		free(dbuf);
		free(t1);
		free(t2);
		free(t3);
		free(t4);
	}
	else
	{
		int h2=src->ih-2, w2=src->iw, size2=h2*w2;
		uint8_t *sbuf=(uint8_t*)malloc(h2);
		uint8_t *dbuf=(uint8_t*)malloc(h2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*h2);
		if(!sbuf||!dbuf||!t1||!t2)
		{
			CRASH("Alloc error");
			return;
		}

		for(kc=0;kc<src->nch;++kc)
		{
			for(kx=0;kx<src->iw;++kx)
			{
				int *ptr=(int32_t*)HDR_DATA(tempimg)+((size_t)kx<<2|kc);
				int bwtidx=
					((ptr[((size_t)src->iw*(h2+(intptr_t)0))<<2]+128)&255)<<0|/*little-endian*/
					((ptr[((size_t)src->iw*(h2+(intptr_t)1))<<2]+128)&255)<<8;
				for(ky=0;ky<h2;++ky)
					sbuf[ky]=(ptr[(size_t)src->iw*ky<<2]+128)&255;
				
				bwtidx%=h2;
				bwt_inv(sbuf, bwtidx, h2, 1<<src->depth[kc], t1, t2, dbuf);
				
				ptr=(int32_t*)HDR_DATA(src)+((size_t)kx<<2|kc);
				for(ky=0;ky<h2;++ky)
					ptr[(size_t)src->iw*ky<<2]=dbuf[ky]-128;
			}
		}
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));/*don't access Image *src until it is reassigned!*/
			if(!ptr)
			{
				CRASH("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}
		src->ih=h2;
		src->iw=w2;

		free(sbuf);
		free(dbuf);
		free(t1);
		free(t2);
	}
	free(tempimg);
	*psrc=src;
}