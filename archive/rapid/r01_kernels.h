#ifdef __OPEN_CL__
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable
#else
#include<stdio.h>
#include<math.h>
#define __kernel
#define __global
#define __constant
//#define half float
#define get_global_id(X) X
#endif

#define CLAMP2(X, LO, HI)\
	do\
	{\
		if(X<LO)X=LO;\
		if(X>HI)X=HI;\
	}while(0)

//worksize: iw*ih			indices: {iw, ih, c0, c1}	{c0, c1} = {00, 40, 31, 22, 13, 04}
__kernel void prep_planes(__global const int *indices, __global unsigned char *image, __global short *planes)
{
	int iw=indices[0], ih=indices[1], c0=indices[2], c1=indices[3];
	int res=iw*ih;
	int idx=get_global_id(0);
	int idx2=idx*3;
	int r=image[idx2+0]-128;
	int g=image[idx2+1]-128;
	int b=image[idx2+2]-128;
	planes[res*0+idx]=r-((c0*g+c1*b)>>2);
	planes[res*1+idx]=g-((c0*b+c1*r)>>2);
	planes[res*2+idx]=b-((c0*r+c1*g)>>2);
}

#define BLOCKSIZE 64
#define NPREDS 11

//worksize: nblocks * 3			indices: {iw, ih, xblocks, yblocks}
__kernel void pred_planes(__global const int *indices, __global const short *planes, __global int *hist)
{
	int iw=indices[0], ih=indices[1], xblocks=indices[2];
	int res=iw*ih;
	int idx=get_global_id(0);
	int coffset, x1, y1;
	int kx, ky;
	__global int *curr_hist=hist+idx*(256*NPREDS);

	for(kx=0;kx<256*NPREDS;++kx)//clear hist
		curr_hist[kx]=0;
	coffset=idx%3*res;
	x1=idx/3%xblocks*BLOCKSIZE;
	y1=idx/(xblocks*3)*BLOCKSIZE;
	for(ky=2;ky<BLOCKSIZE;++ky)
	{
		int y=y1+ky;
		__global const short *NNptr, *Nptr, *currptr;
		if(y>=ih)
			continue;
		NNptr	=planes+coffset+iw*(y-2),
		Nptr	=planes+coffset+iw*(y-1),
		currptr	=planes+coffset+iw*(y+0);
		for(kx=2;kx<BLOCKSIZE-2;++kx)
		{
			int x=x1+kx;
			if(x>iw-2)
				continue;
			int
				NNWW	=NNptr	[x-2],
				NNW	=NNptr	[x-1],
				NN	=NNptr	[x+0],
				NNE	=NNptr	[x+1],
				NNEE	=NNptr	[x+2],
				NWW	=Nptr	[x-2],
				NW	=Nptr	[x-1],
				N	=Nptr	[x+0],
				NE	=Nptr	[x+1],
				NEE	=Nptr	[x+2],
				WW	=currptr[x-2],
				W	=currptr[x-1],
				curr	=currptr[x+0];

			int gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1;
			int gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1;
			int preds[NPREDS]=
			{
				N,
				W,
				(N+W+1)>>1,		//AV2
				(gx*N+gy*W)/(gx+gy),	//WG
				N+W-NW,			//CG
				(3*(N+W)-2*NW+2)>>2,	//AV3
				(4*(N+W)+NE-NW+4)>>3,	//AV4
				W+((5*(N-NW)+NE-WW+4)>>3),	//AV5
				W+((6*N-5*NW+NE-NN-WW+4)>>3),	//AV6
				W+((10*N-9*NW+4*NE-2*(NN+WW)-NNE+NNW-NWW+8)>>4),	//AV9
				(4*NNWW+3*NNW-31*NN-38*NNE + 7*NWW-158*NW+219*N+30*NE+19*NEE - 42*WW+243*W+128)>>8,//AVB
			};
			int vmax=N, vmin=W;

			if(N<W)vmin=N, vmax=W;
			CLAMP2(preds[ 4], vmin, vmax);
			if(vmin<NE)vmin=NE;
			if(vmax>NE)vmax=NE;
			CLAMP2(preds[ 5], vmin, vmax);
			CLAMP2(preds[ 6], vmin, vmax);
			CLAMP2(preds[ 7], vmin, vmax);
			CLAMP2(preds[ 8], vmin, vmax);
			CLAMP2(preds[ 9], vmin, vmax);
			CLAMP2(preds[10], vmin, vmax);

			//sizeof(hist) = 6 * yblocks*xblocks*3 * sizeof(int[NPREDS*256])
			//sizeof(curr_hist) = sizeof(int[NPREDS*256])
			++curr_hist[ 0<<8|((curr-preds[ 0]+128)&255)];
			++curr_hist[ 1<<8|((curr-preds[ 1]+128)&255)];
			++curr_hist[ 2<<8|((curr-preds[ 2]+128)&255)];
			++curr_hist[ 3<<8|((curr-preds[ 3]+128)&255)];
			++curr_hist[ 4<<8|((curr-preds[ 4]+128)&255)];
			++curr_hist[ 5<<8|((curr-preds[ 5]+128)&255)];
			++curr_hist[ 6<<8|((curr-preds[ 6]+128)&255)];
			++curr_hist[ 7<<8|((curr-preds[ 7]+128)&255)];
			++curr_hist[ 8<<8|((curr-preds[ 8]+128)&255)];
			++curr_hist[ 9<<8|((curr-preds[ 9]+128)&255)];
			++curr_hist[10<<8|((curr-preds[10]+128)&255)];
		}
	}
}

//workidx: nblocks * 3 * NPREDS		indices: {dstoffset}
__kernel void calc_entropy(__global const int *indices, __global const int *hist, __global float *csizes)
{
	int dstoffset=indices[0];
	int idx=get_global_id(0);
	__global const int *curr_hist=hist+idx*256;
	int hsum;
	float entropy;
	int ks;

	hsum=0;
	for(ks=0;ks<256;++ks)
		hsum+=curr_hist[ks];
	if(!hsum)
		return;
	entropy=0;
	for(ks=0;ks<256;++ks)
	{
		int freq=curr_hist[ks];
		if(freq)
			entropy-=freq*log2((float)freq/hsum);
	}
	csizes[dstoffset+idx]=entropy*0.125f;
}
