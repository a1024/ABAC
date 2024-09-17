#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;

//"A Block-sorting Lossless Data Compression Algorithm" by M. Burrows and D.J. Wheeler
//https://github.com/MichaelDipperstein/bwt

#define WRAP(VALUE, LIMIT) ((VALUE)-((LIMIT)&-((VALUE)>=(LIMIT))))	//wraps array index within array bounds (assumes value < 2 * limit)
typedef struct _BWTSortInfo
{
	const unsigned char *ptr;
	int count;
} BWTSortInfo;
//https://stackoverflow.com/questions/4210689/pass-extra-parameter-to-comparator-for-qsort
//#ifdef __GNUC__
//static int bwt_cmp_sorted(const void *s1, const void *s2, void *param)//for qsort_r (GNU C, C11)
//#else
static int bwt_cmp_sorted(void *param, const void *s1, const void *s2)//for qsort_s (MSVC)
//#endif
{
	const BWTSortInfo *ctx=(const BWTSortInfo*)param;
	int offset1, offset2;
	int i;

	//Compare 1 character at a time until there's difference or the end of the block is reached.
	//Since we're only sorting strings that already match at the first two characters, start with the third character.
	offset1=*(const int*)s1+2;
	offset2=*(const int*)s2+2;
	for(i=2;i<ctx->count;++i)
	{
		unsigned char c1, c2;

		if(offset1>=ctx->count)offset1-=ctx->count;//ensure that offsets are properly bounded
		if(offset2>=ctx->count)offset2-=ctx->count;

		c1=ctx->ptr[offset1];
		c2=ctx->ptr[offset2];

		if (c1>c2)
			return 1;
		if (c2>c1)
			return -1;

		offset1++;//strings match to here, try next character
		offset2++;
	}
	return 0;//strings are identical
}

//Forward BWT:
//src: byte[count]
//hist:		int[nlevels]
//CDF:		int[nlevels]
//rotationidx:	int[count]
//buf:		int[count]
//dst: byte[count]
//returns: s0idx for reversibility
static int bwt_fwd(unsigned char *ptr, int count, int nlevels, int *hist, int *CDF, int *rotationidx, int *buf, unsigned char *dst)
{
	BWTSortInfo ctx={ptr, count};
	//Sort the rotated strings in the block.
	//A radix sort is performed on the first to characters of all the rotated strings (2nd character then 1st).
	//All rotated strings with matching initial characters are then quicksorted. - Q4..Q7

	//radix sort by second character in rotation
	memset(hist, 0, sizeof(int)*nlevels);
	for(int k=0;k<count;++k)//count number of characters for radix sort
		++hist[ptr[k]];
	for(int sym=0, sum=0;sym<nlevels;++sym)//accumulate CDF
	{
		int freq=hist[sym];
		CDF[sym]=sum;
		sum+=freq;
	}

	{
		int j, k;
		for(k=0;k<count-1;++k)//sort by 2nd character
		{
			j=ptr[k+1];
			buf[CDF[j]]=k;
			++CDF[j];
		}
		j=ptr[0];//handle wrap around for string starting at end of block
		buf[CDF[j]]=k;
		CDF[0]=0;
	}

	//radix sort by first character in rotation
	for(int sym=0, sum=0;sym<nlevels;++sym)//accumulate CDF
	{
		int freq=hist[sym];
		CDF[sym]=sum;
		sum+=freq;
	}
	for(int k=0;k<count;++k)
	{
		int j=ptr[buf[k]];
		rotationidx[CDF[j]]=buf[k];
		++CDF[j];
	}

	//now rotationIdx contains the sort order of all strings sorted by their first 2 characters.
	//Use qsort to sort the strings that have their first two characters matching.
	for(int i=0, k=0;i<nlevels&&k<count-1;++i)
	{
		for(int j=0;j<nlevels&&k<count-1;++j)
		{
			int first=k;
			while(i==ptr[rotationidx[k]]&&j==ptr[WRAP(rotationidx[k]+1, count)])
			{
				++k;
				if(k==count)//we've searched the whole block
					break;
			}
			if(k-first>1)//there are at least 2 strings staring with ij, sort them
				qsort_s(rotationidx+first, (size_t)k-first, sizeof(int), bwt_cmp_sorted, &ctx);
		}
	}

	//find last characters of rotations (L) - C2
	{
		int s0idx=0;
		for(int i=0;i<count;++i)
		{
			if(rotationidx[i])
				dst[i]=ptr[rotationidx[i]-1];
			else
			{
				s0idx=i;//unrotated string 1st character is end of string
				dst[i]=ptr[count-1];
			}
		}
		return s0idx;
	}
}

//Inverse BWT:
//ptr:		byte[count]
//s0idx:	value returned by bwt_fwd
//hist:	int[nlevels]
//pred:	int[count]
//dst:		byte[count]
static void bwt_inv(unsigned char *ptr, int s0idx, int count, int nlevels, int *hist, int *pred, unsigned char *dst)
{
	//based on pseudo code from section 4.2 (D1 and D2)
	memset(hist, 0, sizeof(int)*nlevels);
	for(int k=0;k<count;++k)
	{
		pred[k]=hist[ptr[k]];		//Set pred[i] to the number of times block[i] appears in the substring block[0 .. i - 1].
		++hist[ptr[k]];			//As a useful side effect count[i] will be the number of times character i appears in block.
	}
	for(int k=0, sum=0;k<nlevels;++k)	//Finally, set count[i] to the number of characters in block lexicographically less than i.
	{
		int freq=hist[k];
		hist[k]=sum;
		sum+=freq;
	}
	for(int i=s0idx, j=count;j>0;--j)	//construct the initial unrotated string (S[0])
	{
		dst[j-1]=ptr[i];
		i=pred[i]+hist[ptr[i]];
	}
}

//static void image_add_offset(Image *src, int offset, int mask)
//{
//	int res=src->iw*src->ih;
//	for(int kc=0;kc<src->nch;++kc)
//	{
//		for(int k=0;k<res;++k)
//		{
//			int val=src->data[k<<2|kc];
//			val+=offset;
//			val&=mask;
//			src->data[k<<2|kc]=val;
//		}
//	}
//}
void prep_BWT_x(Image **psrc, int fwd)
{
	Image *src=*psrc;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	//if(src->depth[0]!=8)
	//{
	//	LOG_ERROR("BWT: Expected 8-bit depth");
	//	return;
	//}
	if(src->iw>0x10000||src->ih>0x10000)
	{
		LOG_ERROR("BWT: Image dimensions must be up to 65536");
		return;
	}
	if(fwd)
	{
		int h2=src->ih, w2=src->iw+2, size2=h2*w2;
		unsigned char *sbuf=(unsigned char*)malloc(w2);
		unsigned char *dbuf=(unsigned char*)malloc(w2);
	//	unsigned char *debugbuf=(unsigned char*)malloc(w2);//
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*maxlevels);
		int *t3=(int*)malloc(sizeof(int)*w2);
		int *t4=(int*)malloc(sizeof(int)*w2);
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));//don't access Image *src until it is reassigned!
			if(!ptr||!sbuf||!dbuf||!t1||!t2||!t3||!t4)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}

#if 0
		{//
			const char msg[]=
			//	"BANANA"
				"The Burrows–Wheeler transform (BWT, also called block-sorting compression) rearranges a character string into runs of similar characters."
				;
			const int size=(int)sizeof(msg)-1;
			memcpy(sbuf, msg, size+1LL);

			int bwtidx=bwt_fwd(sbuf, size, 256, t1, t2, t3, t4, dbuf);
			bwt_inv(dbuf, bwtidx, size, 256, t1, t3, debugbuf);//
			if(memcmp(debugbuf, sbuf, size))
				LOG_ERROR("BWT Test Error");
			printf("");
		}
#endif
		for(int kc=0;kc<src->nch;++kc)
		{
			for(int ky=src->ih-1;ky>=0;--ky)
			{
				int bwtidx;
				int *ptr=src->data+((size_t)src->iw*ky<<2|kc);
				for(int kx=0;kx<src->iw;++kx)
					sbuf[kx]=(ptr[kx<<2]+128)&255;

				bwtidx=bwt_fwd(sbuf, src->iw, 1<<src->depth[kc], t1, t2, t3, t4, dbuf);

#if 0
				bwt_inv(dbuf, bwtidx, src->iw, 1<<src->depth[kc], t1, t3, debugbuf);//
				if(memcmp(debugbuf, sbuf, src->iw))
					LOG_ERROR("BWT Error ky=%d", ky);
#endif

				ptr=src->data+((size_t)w2*ky<<2|kc);
				for(int kx=0;kx<src->iw;++kx)
					ptr[kx<<2]=dbuf[kx]-128;
				ptr[(src->iw+0LL)<<2]=(bwtidx>>0&255)-128;//little-endian
				ptr[(src->iw+1LL)<<2]=(bwtidx>>8&255)-128;
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
	//	free(debugbuf);//
	}
	else
	{
		int h2=src->ih, w2=src->iw-2, size2=h2*w2;
		unsigned char *sbuf=(unsigned char*)malloc(w2);
		unsigned char *dbuf=(unsigned char*)malloc(w2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*w2);
		if(!sbuf||!dbuf||!t1||!t2)
		{
			LOG_ERROR("Alloc error");
			return;
		}

		for(int kc=0;kc<src->nch;++kc)
		{
			for(int ky=0;ky<src->ih;++ky)
			{
				int *ptr=src->data+((size_t)src->iw*ky<<2|kc);
				int bwtidx=
					((ptr[(w2+0LL)<<2]+128)&255)<<0|//little-endian
					((ptr[(w2+1LL)<<2]+128)&255)<<8;
				for(int kx=0;kx<w2;++kx)
					sbuf[kx]=(ptr[kx<<2]+128)&255;

				//if((unsigned)bwtidx>=(unsigned)w2)
				//{
				//	LOG_ERROR("Don't click Inv BWT without Fwd BWT!");
				//	return;
				//}
				bwtidx%=w2;
				bwt_inv(sbuf, bwtidx, w2, 1<<src->depth[kc], t1, t2, dbuf);
				
				ptr=src->data+((size_t)w2*ky<<2|kc);
				for(int kx=0;kx<w2;++kx)
					ptr[kx<<2]=dbuf[kx]-128;
			}
		}
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));//don't access Image *src until it is reassigned!
			if(!ptr)
			{
				LOG_ERROR("Alloc error");
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
	if(src->iw>0x10000||src->ih>0x10000)
	{
		LOG_ERROR("BWT: Image dimensions must be up to 65536");
		return;
	}
	image_copy(&tempimg, src);
	if(!tempimg)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	if(fwd)
	{
		int h2=src->ih+2, w2=src->iw, size2=h2*w2;
		unsigned char *sbuf=(unsigned char*)malloc(h2);
		unsigned char *dbuf=(unsigned char*)malloc(h2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*maxlevels);
		int *t3=(int*)malloc(sizeof(int)*h2);
		int *t4=(int*)malloc(sizeof(int)*h2);
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));//don't access Image *src until it is reassigned!
			if(!ptr||!sbuf||!dbuf||!t1||!t2||!t3||!t4)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			src=(Image*)ptr;
		}

		for(int kc=0;kc<src->nch;++kc)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				int bwtidx;
				int *ptr=tempimg->data+((size_t)kx<<2|kc);
				for(int ky=0;ky<src->ih;++ky)
					sbuf[ky]=(ptr[(size_t)src->iw*ky<<2]+128)&255;

				bwtidx=bwt_fwd(sbuf, src->ih, 1<<src->depth[kc], t1, t2, t3, t4, dbuf);

				ptr=src->data+((size_t)kx<<2|kc);
				for(int ky=0;ky<src->ih;++ky)
					ptr[(size_t)w2*ky<<2]=dbuf[ky]-128;
				ptr[((size_t)w2*(src->ih+0LL))<<2]=(bwtidx>>0&255)-128;//little-endian
				ptr[((size_t)w2*(src->ih+1LL))<<2]=(bwtidx>>8&255)-128;
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
		unsigned char *sbuf=(unsigned char*)malloc(h2);
		unsigned char *dbuf=(unsigned char*)malloc(h2);
		int *t1=(int*)malloc(sizeof(int)*maxlevels);
		int *t2=(int*)malloc(sizeof(int)*h2);
		if(!sbuf||!dbuf||!t1||!t2)
		{
			LOG_ERROR("Alloc error");
			return;
		}

		for(int kc=0;kc<src->nch;++kc)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				int *ptr=tempimg->data+((size_t)kx<<2|kc);
				int bwtidx=
					((ptr[((size_t)src->iw*(h2+0LL))<<2]+128)&255)<<0|//little-endian
					((ptr[((size_t)src->iw*(h2+1LL))<<2]+128)&255)<<8;
				for(int ky=0;ky<h2;++ky)
					sbuf[ky]=(ptr[(size_t)src->iw*ky<<2]+128)&255;
				
				bwtidx%=h2;
				bwt_inv(sbuf, bwtidx, h2, 1<<src->depth[kc], t1, t2, dbuf);
				
				ptr=src->data+((size_t)kx<<2|kc);
				for(int ky=0;ky<h2;++ky)
					ptr[(size_t)src->iw*ky<<2]=dbuf[ky]-128;
			}
		}
		{
			void *ptr=realloc(src, sizeof(Image)+size2*sizeof(int[4]));//don't access Image *src until it is reassigned!
			if(!ptr)
			{
				LOG_ERROR("Alloc error");
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