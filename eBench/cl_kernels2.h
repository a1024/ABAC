#ifdef _MSC_VER
#include<stdio.h>
#include<math.h>
#define __kernel
#define __global
#define __constant
#define get_global_id(X) X
#define ABS fabsf
#define LOG2 log2f
#else
#define ABS fabs
#define LOG2 log2
#endif
#define CLAMP(LO, X, HI)    ((X)>(LO)?(X)<(HI)?(X):(HI):(LO))


//CUSTOM4 (CUSTOM3 with GPU training)

//[0]RO	int   indices[3] = {reach, iw, ih}
//[1]RO	short params[NTHREADS*NPARAMS]
//[2]RO	char  pixels[IW*IH]
//[3]RW	char  allerrors[NTHREADS*(REACH+1)*IW]
//[4]RW	short neighbors[NTHREADS*(KSIZE+2)*3]
//[5]RW	int   histograms[NTHREADS*768]
//[6]WO	float invCRs[NTHREADS*3]
__kernel void custom3_eval(__constant int *indices, __constant short *params, __constant char *pixels, __global char *allerrors, __global short *neighbors, __global int *histograms, __global float *invCRs)
{
	int threadidx=get_global_id(0);
	int reach=indices[0], iw=indices[1], ih=indices[2], res=iw*ih;
	int nnb=reach*(reach+1)*2, ksize=nnb*2, nparams=ksize*9+6, eh=reach+1;
	__constant short *filter=params+nparams*threadidx;
	__constant short *coeffs[]=
	{
		filter+ksize*0+0, filter+ksize*1+0, filter+ksize*2+0,
		filter+ksize*3+0, filter+ksize*4+2, filter+ksize*5+2,
		filter+ksize*6+2, filter+ksize*7+4, filter+ksize*8+6,
	};
	__global char *errors=allerrors+eh*iw*threadidx;
	__global short *nb[]=
	{
		neighbors+(ksize+2)*3*threadidx,
		neighbors+(ksize+2)*3*threadidx+ksize+2,
		neighbors+(ksize+2)*3*threadidx+(ksize+2)*2,
	};
	__global int *hist=histograms+768*threadidx;
	int count[3];

	for(int k=0, end=eh*iw;k<end;++k)
		errors[k]=0;
	for(int k=0;k<768;++k)//clear histogram
		hist[k]=0;

	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx, pred;
			for(int kc=0;kc<3;++kc)//load neighbors
			{
				__global short *nbk=nb[kc];
				idx=-1;
				for(int ky2=-reach;ky2<0;++ky2)
				{
					int y=ky+ky2;
					for(int kx2=-reach;kx2<=reach;++kx2)
					{
						int x=kx+kx2;
						if((unsigned)x<(unsigned)iw&&(unsigned)y<(unsigned)ih)
						{
							nbk[++idx]=pixels[(iw*y+x)<<2|kc];
							nbk[++idx]=errors[(iw*(y%eh)+x)<<2|kc];
						}
						else
						{
							nbk[++idx]=0;
							nbk[++idx]=0;
						}
					}
				}
				for(int kx2=-reach;kx2<0;++kx2)
				{
					int x=kx+kx2;
					if((unsigned)x<(unsigned)iw)
					{
						nbk[++idx]=pixels[(iw*ky+x)<<2|kc];
						nbk[++idx]=errors[(iw*(ky%eh)+x)<<2|kc];
					}
					else
					{
						nbk[++idx]=0;
						nbk[++idx]=0;
					}
				}
				count[kc]=++idx;
			}

			idx=0;
			for(int kdst=0;kdst<3;++kdst)
			{
				int srcidx=(iw*ky+kx)<<2|kdst,
					dstidx=(iw*(ky%eh)+kx)<<2|kdst;

				pred=0;
				for(int kc=0;kc<3;++kc)
				{
					__constant short *a=coeffs[idx+kc];
					__global short *b=nb[kc];
					for(int k=0;k<count[kc];++k)
						pred+=a[k]*b[k];
				}
				pred+=1<<13;
				pred>>=14;
				pred=CLAMP(-128, pred, 127);

				errors[dstidx]=pixels[srcidx]-pred;

				++hist[kdst<<8|(errors[dstidx]+128)];
				
				nb[kdst][ksize  ]=pixels[srcidx];
				nb[kdst][ksize+1]=errors[dstidx];
				count[kdst]+=2;
				idx+=3;
			}
		}
	}
	for(int kc=0;kc<3;++kc)
	{
		float entropy=0;
		for(int sym=0;sym<256;++sym)
		{
			int freq=hist[kc<<8|sym];
			if(freq)
			{
				float p=(float)freq/res;
				entropy-=p*LOG2(p);
			}
		}
		invCRs[threadidx*3+kc]=entropy/8;
	}
}