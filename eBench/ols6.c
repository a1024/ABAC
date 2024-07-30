#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


	#define ENABLE_V3
//	#define ENABLE_V2
//	#define ENABLE_LOGGING


#define V3_NPARAMS 16

#define OLS6_REACH 1
#define OLS6_STEP 1
#define OLS6_SAMPLEREACH 8
#define OLS6_NPARAMS0 (2*(OLS6_REACH+1)*OLS6_REACH*3+0+1)
#define OLS6_NPARAMS1 (2*(OLS6_REACH+1)*OLS6_REACH*3+1+1)
#define OLS6_NPARAMS2 (2*(OLS6_REACH+1)*OLS6_REACH*3+2+1)
//#define OLS6_NPARAMS_TOTAL (2*(OLS6_REACH+1)*OLS6_REACH*9+3+3)
#define OLS6_NSAMPLES ((2*OLS6_SAMPLEREACH+OLS6_STEP+1)*OLS6_SAMPLEREACH)
static int solve_Mx_v_cholesky(double *matrix, const double *vec, int n, double *solution)
{
	int success=1;
	double sum;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<i;++j)
		{
			sum=matrix[i*n+j];
			for(int k=0;k<j;++k)
				sum-=matrix[i*n+k]*matrix[j*n+k];
			matrix[i*n+j]=sum/matrix[j*n+j];
		}
		sum=matrix[i*n+i];
		for(int k=0;k<i;++k)
			sum-=matrix[i*n+k]*matrix[i*n+k];
		if(sum<=1e-8)
		{
			success=0;
			break;
		}
		matrix[i*n+i]=sqrt(sum);
	}
	if(success)
	{
		for(int ky=0;ky<n;++ky)
		{
			//solution[0] = vec[0]/matrix[0][0]
			//solution[1] = (vec[1] - solution[0]*matrix[1][0])/matrix[1][1]
			//solution[2] = (vec[2] - solution[0]*matrix[2][0] - solution[1]*matrix[2][1])/matrix[2][2]
			//...
			sum=vec[ky];
			for(int kx=0;kx<ky;++kx)//lower triangle downwards	_*\ |
				sum-=matrix[ky*n+kx]*solution[kx];
			solution[ky]=sum/matrix[(n+1)*ky];
		}
		for(int kx=n-1;kx>=0;--kx)
		{
			sum=solution[kx];
			for(int ky=kx+1;ky<n;++ky)//upper triangle upwards	_ \*|
				sum-=matrix[ky*n+kx]*solution[ky];
			solution[kx]=sum/matrix[(n+1)*kx];
		}
	}
	return success;
}
static int solve_Mx_v(double *matrix, double *vec, int size)
{
	int success=1;
	for(int it=0;it<size;++it)
	{
		int kp;
		double pivot;

		kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[size*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			double temp;
			for(int k=it;k<size;++k)
				SWAPVAR(matrix[size*it+k], matrix[size*kp+k], temp);
			SWAPVAR(vec[it], vec[kp], temp);
		}
		pivot=matrix[size*it+it];
		for(int kx=it;kx<size;++kx)
			matrix[size*it+kx]/=pivot;
		vec[it]/=pivot;
		for(int ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[size*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(int kx=it;kx<size;++kx)
					matrix[size*ky+kx]-=matrix[size*it+kx]*factor;
				vec[ky]-=vec[it]*factor;
			}
		}
	}
	return success;
}
static int invert_matrix(double *matrix, int size, double *temprow)
{
	int success;

	success=1;
	for(int it=0;it<size;++it)
	{
		int kp;
		double pivot;

		kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[((size_t)size<<1)*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			memcpy(temprow, matrix+((size_t)size<<1)*it, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*it, matrix+((size_t)size<<1)*kp, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*kp, temprow, size*sizeof(double[2]));
		}
		pivot=matrix[((size_t)size<<1)*it+it];
		for(int kx=it;kx<(size<<1);++kx)
			matrix[((size_t)size<<1)*it+kx]/=pivot;
		for(int ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[((size_t)size<<1)*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(int kx=it;kx<(size<<1);++kx)
					matrix[((size_t)size<<1)*ky+kx]-=matrix[((size_t)size<<1)*it+kx]*factor;
			}
		}
	}
	//_freea(temp);
	return success;
}
static int ols6_fallbackpred(const int *pixels, int iw, int kc, int kx, int ky, int idx)
{
	int
		N =ky?pixels[(idx-iw)<<2|kc]:0,
		W =kx?pixels[(idx-1)<<2|kc]:0,
		NW=ky&&kx?pixels[(idx-iw-1)<<2|kc]:0;
	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);
	return pred;
}
static void ols6_add_sample(const int *pixels, double gain, int iw, int kc, int kx, int ky, double *sample, double *matrix, int remove)
{
	int nparams=OLS6_NPARAMS0+kc;
	if(remove)
	{
		for(int ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]-=vy*vx;
			}
		}
	}
	sample[0]=1;
	for(int ky3=-OLS6_REACH, idx2=1;ky3<=0;++ky3)
	{
		for(int kx3=-OLS6_REACH;kx3<=OLS6_REACH;++kx3)
		{
			for(int kc2=0;kc2<3;++kc2, ++idx2)
			{
				//if(idx2>=nparams+1)
				//	LOG_ERROR("");
				sample[idx2]=(double)pixels[(iw*(ky+ky3)+kx+kx3)<<2|kc2]*gain;
				if(!ky3&&!kx3&&kc2==kc)//last element is target
					goto finish;
			}
		}
	}
finish:
	if(!remove)
	{
		for(int ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]+=vy*vx;
			}
		}
	}
}
void pred_ols6(Image *src, int fwd)
{
#ifdef ENABLE_V3
	static const int xindices[]=
	{
		0,
		0-1*8, 1-1*8, 2-1*8,
		0+0*8, 1+0*8, 2+0*8,
		0+1*8, 1+1*8, 2+1*8,
		0-1*8, 1-1*8, 2-1*8,
		0+0*8, 1+0*8, 2+0*8,
	};
	static const int yindices[]=
	{
		0,
		1, 1, 1,
		1, 1, 1,
		1, 1, 1,
		0, 0, 0,
		0, 0, 0,
	};
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
	};
	double gains[]=
	{
		1./half[0],
		1./half[1],
		1./half[2],
		//1,
		//1,
		//1,
	};
	int failcount=0;
	double matrix[V3_NPARAMS*V3_NPARAMS];
	double vec[V3_NPARAMS];
	double params[V3_NPARAMS*3]={0};
	int bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);//16 padded rows * 4 channels max * {pixels, preds/errors}
	int *pixels=(int*)malloc(bufsize);
	int samplessize=(src->iw+6)*(int)sizeof(double[(V3_NPARAMS+1)*3]);
	double *samples=(double*)malloc(samplessize);
	double t_start=time_sec();
	if(!pixels||!samples)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	//memset(samples, 0, sizeof(double[V3_NPARAMS+1]));
	//*samples=1;
	//memfill(samples+V3_NPARAMS+1, samples, samplessize-sizeof(double[V3_NPARAMS+1]), sizeof(double[V3_NPARAMS+1]));
	memset(samples, 0, samplessize);
	memset(pixels, 0, bufsize);
	for(int ky=0, idx=0;ky<src->ih;++ky)//O(H * W * C * NPARAMS * NPARAMS)
	{
		ALIGN(16) int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			//int
			//	*NNWW	=rows[2]-2*8,
			//	*NNW	=rows[2]-1*8,
			//	*NN	=rows[2]+0*8,
			//	*NNE	=rows[2]+1*8,
			//	*NNEE	=rows[2]+2*8,
			//	*NWW	=rows[1]-2*8,
			//	*NW	=rows[1]-1*8,
			//	*N	=rows[1]+0*8,
			//	*NE	=rows[1]+1*8,
			//	*NEE	=rows[1]+2*8,
			//	*NEEE	=rows[1]+3*8,
			//	*WW	=rows[0]-2*8,
			//	*W	=rows[0]-1*8,
			//	*curr	=rows[0]+0*8;
			for(int kc=0;kc<3;++kc)
			{
				int
					NNWW	=rows[2][kc-2*8+0],
					NNW	=rows[2][kc-1*8+0],
					NN	=rows[2][kc+0*8+0],
					NNE	=rows[2][kc+1*8+0],
					NNEE	=rows[2][kc+2*8+0],
					NWW	=rows[1][kc-2*8+0],
					NW	=rows[1][kc-1*8+0],
					N	=rows[1][kc+0*8+0],
					NE	=rows[1][kc+1*8+0],
					NEE	=rows[1][kc+2*8+0],
					NEEE	=rows[1][kc+3*8+0],
					WW	=rows[0][kc-2*8+0],
					W	=rows[0][kc-1*8+0],
					eNNWW	=rows[2][kc-2*8+4],
					eNNW	=rows[2][kc-1*8+4],
					eNN	=rows[2][kc+0*8+4],
					eNNE	=rows[2][kc+1*8+4],
					eNNEE	=rows[2][kc+2*8+4],
					eNWW	=rows[1][kc-2*8+4],
					eNW	=rows[1][kc-1*8+4],
					eN	=rows[1][kc+0*8+4],
					eNE	=rows[1][kc+1*8+4],
					eNEE	=rows[1][kc+2*8+4],
					eNEEE	=rows[1][kc+3*8+4],
					eWW	=rows[0][kc-2*8+4],
					eW	=rows[0][kc-1*8+4];
				int pred, success;
				const int nparams=12;
				//int nparams=OLS6_NPARAMS0+kc;
				double *curr_params=params+kc*V3_NPARAMS;
				double *samplesWW	=samples+(V3_NPARAMS+1)*((kx+0)*3+kc);
				double *samplesW	=samples+(V3_NPARAMS+1)*((kx+1)*3+kc);
				double *samplesN	=samples+(V3_NPARAMS+1)*((kx+2)*3+kc);
				double *samplesNE	=samples+(V3_NPARAMS+1)*((kx+3)*3+kc);
				double *samplesNEE	=samples+(V3_NPARAMS+1)*((kx+4)*3+kc);
				double *samplesNEEE	=samples+(V3_NPARAMS+1)*((kx+5)*3+kc);
			//	double *samplesNEEEE	=samples+(V3_NPARAMS+1)*((kx+6)*3+kc);
				//samplesW[0]=1;
				//samplesN[0]=1;
				//samplesNEEE[0]=1;
				//params = inv(sum i: sample[i] o sample[i]) * (sum i: sample[i] .* <target[i], ...>)
				memset(matrix, 0, sizeof(matrix));
				memset(vec, 0, sizeof(vec));
				for(int my=0;my<nparams;++my)
				{
					for(int mx=0;mx<nparams;++mx)
						matrix[my*nparams+mx]+=samplesW[my]*samplesW[mx];
					vec[my]+=samplesW[my]*samplesW[nparams];
				}
				//for(int my=0;my<nparams;++my)
				//{
				//	for(int mx=0;mx<nparams;++mx)
				//		matrix[my*nparams+mx]+=samplesNE[my]*samplesNE[mx];
				//	vec[my]+=samplesNE[my]*samplesNE[nparams];
				//}
				//if(ky==50&&kx==50)
				//	printf("");
				for(int k=0;k<nparams;++k)
					matrix[(nparams+1)*k]+=0.005;
				memset(curr_params, 0, sizeof(double[V3_NPARAMS]));
				if(1)
					success=solve_Mx_v_cholesky(matrix, vec, nparams, curr_params);
				else
				{
					success=solve_Mx_v(matrix, vec, nparams);
					memcpy(curr_params, vec, nparams*sizeof(double));
					//for(int k=0;k<nparams;++k)
					//	curr_params[k]+=(vec[k]-curr_params[k])*0.2;
				}
				double preds[V3_NPARAMS+1]=
				{
					//v3.4
#if 1
					N+W-NW,				//+0x2E00000
					N+NE-NNE,			//+0x1800000
					(3*(3*W+NE+NEE)-10*N)/5.,	//+0x1400000
					2*N-NN,				//+0x1400000
					2*W-WW,				//+0x1400000
					(4*W+N+NE+NEE+NEEE)/8.,		//+0x1000000
					N+eN/2.,			//+0x0F00000
					W+eW/2.,			//+0x0F00000
					(3*W+NEEE)/4.,			//+0x0F00000
					(4*(N+W)-NE+NW)/8.,		//+0x0800000
					(2*W+NEE-N)/2.,			//+0x0800000
					(3*W+NEEE-N)/3.,		//+0x0800000
					W+NE-N,				//+0x0500000
					W+NE-N+NW-(NWW+NN-NNW),		//+0x0500000
					(2*(W+NE-N)-(WW+NNEE-NN))/2.,	//+0x0200000
					(N+W)/2.,			//-0x0800000
				//	(W+NE)/2.,			//-0x0E00000
				//	(2*W+NEE)/3.,			//-0x0E00000
				//	rows[0][0],
				//	rows[0][1],
				//	rows[0][2],
					0,//target
#endif

					//v3.3
#if 0
					half[kc],
					NW	[kc],
					N[kc]+W[kc]-NW[kc],
					(N[kc]+W[kc])/2,
					N	[0],	N	[1],	N	[2],
					NE	[kc],
					(W[kc]+NE[kc])/2,
					W[kc]+NE[kc]-N[kc],
					W	[0],	W	[1],	W	[2],
					rows[0]	[0],	rows[0]	[1],	rows[0]	[2],
#endif

					//v3.2
#if 0
					half[kc],
					NW	[0],	NW	[1],	NW	[2],
					N	[0],	N	[1],	N	[2],
					NE	[0],	NE	[1],	NE	[2],
					W	[0],	W	[1],	W	[2],
					rows[0]	[0],	rows[0]	[1],	rows[0]	[2],
#endif
				};
				if(success)
				{
					double fpred=0;
					for(int ks=0;ks<nparams;++ks)
						fpred+=curr_params[ks]*preds[ks];
					//double fpred=curr_params[0]/gains[kc];
					//for(int ks=1;ks<nparams;++ks)
					//	fpred+=curr_params[ks]*rows[yindices[ks]][xindices[ks]];
					CLAMP2_32(pred, (int)fpred, -half[kc], half[kc]-1);
				}
				else
				{
					MEDIAN3_32(pred, rows[1][kc], rows[0][kc-8], rows[1][kc]+rows[0][kc-8]-rows[1][kc-8]);
					++failcount;
				}
				pred^=-fwd;
				pred+=fwd;
				pred+=rows[0][kc+!fwd*4]=src->data[idx<<2|kc];
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				rows[0][kc+ fwd*4]=src->data[idx<<2|kc]=pred;

				preds[nparams]=rows[0][kc];
				//	[0] 1		2 3		4 5		6 7		8 9
				//	10 11		12 13		[14] 15		[16] 17		[18] 19
				//	[20] 21		[22] 23		?

				//int wsum=
				//	custom_params[0]+
				//	custom_params[14]+
				//	custom_params[16]+
				//	custom_params[18]+
				//	custom_params[20]+
				//	custom_params[22];
				//if(wsum)
				for(int ks=0;ks<nparams+1;++ks)
					samplesN[ks]=(2*samplesW[ks]+preds[ks]*gains[kc]+samplesNEEE[ks])*0.25;
				//	samplesN[ks]=(
				//		custom_params[0]*preds[ks]*gains[kc]+
				//		custom_params[14]*samplesN[ks]+
				//		custom_params[16]*samplesNE[ks]+
				//		custom_params[18]*samplesNEE[ks]+
				//		custom_params[20]*samplesWW[ks]+
				//		custom_params[22]*samplesW[ks]
				//	)/65536;
				//	samplesN[ks]=(2*samplesN[ks]+samplesNE[ks]+4*samplesNEE[ks]+2*samplesWW[ks]+6*samplesW[ks]+preds[ks]*gains[kc])*(1./16);
				//	samplesN[ks]=preds[ks];
				//for(int ks=1;ks<nparams+1;++ks)
				//	samplesN[ks]=(2*samplesW[ks]+rows[yindices[ks]][xindices[ks]]+samplesNEEE[ks])*0.25;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
		if(loud_transforms&&(!((ky+1)&63)||ky==src->ih-1))
			set_window_title("%d/%d = %7.3lf%%  failrate %lf%%  %lf sec",
				ky+1,
				src->ih,
				100.*(ky+1)/src->ih,
				100.*failcount/(src->nch*src->iw*(ky+1)),
				time_sec()-t_start
			);
	}
	free(pixels);
	free(samples);
#elif defined ENABLE_V2
	static const int xindices[]=
	{
		0,
		0-1*8, 1-1*8, 2-1*8,
		0+0*8, 1+0*8, 2+0*8,
		0+1*8, 1+1*8, 2+1*8,
		0-1*8, 1-1*8, 2-1*8,
		0+0*8, 1+0*8, 2+0*8,
	};
	static const int yindices[]=
	{
		0,
		1, 1, 1,
		1, 1, 1,
		1, 1, 1,
		0, 0, 0,
		0, 0, 0,
	};
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
	};
	double gains[]=
	{
		1./half[0],
		1./half[1],
		1./half[2],
		//1,
		//1,
		//1,
	};
	int failcount=0;
	double params[OLS6_NPARAMS2*3]={0};
	double nb[OLS6_NPARAMS2+1];
	double vec[OLS6_NPARAMS2];
	double matrix[OLS6_NPARAMS2*OLS6_NPARAMS2];
	int bufsize=(src->iw+32)*(int)sizeof(int[16*4*2]);//16 padded rows * 4 channels max * {pixels, preds/errors}
	int *pixels=(int*)malloc(bufsize);
	double t_start=time_sec();
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	memset(pixels, 0, bufsize);
	for(int ky=0, idx=0;ky<src->ih;++ky)//O(H * W * C * DY * DX * NPARAMS * NPARAMS)
	{
		ALIGN(16) int *rows[]=
		{
			pixels+((src->iw+32LL)*((ky- 0LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 1LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 2LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 3LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 4LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 5LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 6LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 7LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 8LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky- 9LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-10LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-11LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-12LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-13LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-14LL)&15)+16LL)*4*2,
			pixels+((src->iw+32LL)*((ky-15LL)&15)+16LL)*4*2,
		};
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int pred, success;
				int nparams=OLS6_NPARAMS0+kc;
				double *curr_params=params+kc*OLS6_NPARAMS2;
				memset(matrix, 0, sizeof(matrix));
				memset(vec, 0, sizeof(vec));
				//params = inv(sum i: sample[i] o sample[i]) * (sum i: sample[i] .* <target[i], ...>)
				for(int dy=-OLS6_SAMPLEREACH;dy<=0;++dy)
				{
					for(int dx=-OLS6_SAMPLEREACH;dx<=OLS6_SAMPLEREACH;++dx)
					{
						if(!dy&&!dx)
							break;
						nb[0]=1;
						for(int ks=1;ks<nparams+1;++ks)
							nb[ks]=rows[yindices[ks]-dy][xindices[ks]+dx]*gains[kc];
						for(int my=0;my<nparams;++my)
						{
							for(int mx=0;mx<nparams;++mx)
								matrix[my*nparams+mx]+=nb[my]*nb[mx];
							vec[my]+=nb[my]*nb[nparams];
						}
					}
				}
				for(int k=0;k<nparams;++k)
					matrix[(nparams+1)*k]+=0.005;
				if(1)
					success=solve_Mx_v_cholesky(matrix, vec, nparams, curr_params);
				else
				{
					success=solve_Mx_v(matrix, vec, nparams);
					for(int k=0;k<nparams;++k)
						curr_params[k]+=(vec[k]-curr_params[k])*0.2;
					//memcpy(curr_params, vec, nparams*sizeof(double));
				}
				if(success)
				{
					double fpred=curr_params[0]/gains[kc];
					for(int ks=1;ks<nparams;++ks)
						fpred+=curr_params[ks]*rows[yindices[ks]][xindices[ks]];
					CLAMP2_32(pred, (int)(fpred+0.5), -half[kc], half[kc]-1);
				}
				else
				{
					MEDIAN3_32(pred, rows[1][kc], rows[0][kc-8], rows[1][kc]+rows[0][kc-8]-rows[1][kc-8]);
					++failcount;
				}
				int curr=src->data[idx<<2|kc];
				rows[0][kc+!fwd*4]=curr;
				pred^=-fwd;
				pred+=fwd;
				pred+=curr;
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				src->data[idx<<2|kc]=pred;
				rows[0][kc+ fwd*4]=pred;
			}
			rows[ 0]+=8;
			rows[ 1]+=8;
			rows[ 2]+=8;
			rows[ 3]+=8;
			rows[ 4]+=8;
			rows[ 5]+=8;
			rows[ 6]+=8;
			rows[ 7]+=8;
			rows[ 8]+=8;
			rows[ 9]+=8;
			rows[10]+=8;
			rows[11]+=8;
			rows[12]+=8;
			rows[13]+=8;
			rows[14]+=8;
			rows[15]+=8;
		}
		if(loud_transforms)
			set_window_title("%d/%d = %7.3lf%%  failrate %lf%%  %lf sec",
				ky+1,
				src->ih,
				100.*(ky+1)/src->ih,
				100.*failcount/(src->nch*src->iw*(ky+1)),
				time_sec()-t_start
			);
	}
	free(pixels);
#else
	double t_start;
	int successcount;
	Image *dst;
	double *params, *samples, *matrix, *matrix2, *temp;
	const int *pixels;
//	const int *errors;
//	double avparams[OLS6_NPARAMS2*3]={0};
//	double params[OLS6_NPARAMS2*3]={0};
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
	};
	double gains[]=
	{
		1./half[0],
		1./half[1],
		1./half[2],
		//1,
		//1,
		//1,
	};
	size_t paramssize=(src->iw+1LL)*sizeof(double[OLS6_NPARAMS2*3]);
	int paramidx=0;

	t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
#ifdef ENABLE_LOGGING
	console_start();
#endif
	successcount=0;
	dst=0;
	image_copy(&dst, src);
	params=(double*)malloc(paramssize);
	samples=(double*)malloc(sizeof(double[OLS6_NSAMPLES*(OLS6_NPARAMS2+1)*3]));
	matrix=(double*)malloc(sizeof(double[OLS6_NPARAMS2*OLS6_NPARAMS2*3<<1]));
	matrix2=(double*)malloc(sizeof(double[OLS6_NPARAMS2*OLS6_NPARAMS2*3<<1]));
	temp=(double*)malloc(sizeof(double[OLS6_NPARAMS2<<1]));
	if(!dst||!params||!samples||!matrix||!matrix2||!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	double *pend=params+(src->iw+1LL)*OLS6_NPARAMS2*3;
	memset(params, 0, paramssize);
	memset(samples, 0, sizeof(double[OLS6_NSAMPLES*(OLS6_NPARAMS2+1)*3]));
	pixels=fwd?dst->data:src->data;
//	errors=fwd?src->data:dst->data;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int initialized[3]={0};
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			for(int kc=0;kc<3;++kc)
			{
				double *cparamsW	=params+OLS6_NPARAMS2*((paramidx+0)%(src->iw+1)*3+kc);
				double *cparamsNW	=params+OLS6_NPARAMS2*((paramidx+1)%(src->iw+1)*3+kc);
				double *cparamsN	=params+OLS6_NPARAMS2*((paramidx+2)%(src->iw+1)*3+kc);
				double *cparamsNE	=params+OLS6_NPARAMS2*((paramidx+3)%(src->iw+1)*3+kc);
				double *cparamsNEE	=params+OLS6_NPARAMS2*((paramidx+4)%(src->iw+1)*3+kc);
				double *cparamsNEEE	=params+OLS6_NPARAMS2*((paramidx+5)%(src->iw+1)*3+kc);
				//if(
				//	paramidx<0||
				//	cparamsW	+OLS6_NPARAMS2>pend||
				//	cparamsNW	+OLS6_NPARAMS2>pend||
				//	cparamsN	+OLS6_NPARAMS2>pend||
				//	cparamsNE	+OLS6_NPARAMS2>pend||
				//	cparamsNEE	+OLS6_NPARAMS2>pend||
				//	cparamsNEEE	+OLS6_NPARAMS2>pend
				//)
				//	LOG_ERROR("");
			//	double *cparams=params+OLS6_NPARAMS2*kc;
				double *csamples=samples+(OLS6_NPARAMS2+1)*OLS6_NSAMPLES*kc;
				int nparams=OLS6_NPARAMS0+kc;
				ALIGN(16) int pred[4];

				if(kx<OLS6_SAMPLEREACH+OLS6_REACH||kx>src->iw-OLS6_SAMPLEREACH-OLS6_REACH-OLS6_STEP||ky<OLS6_SAMPLEREACH+OLS6_REACH)
					pred[0]=ols6_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
				else
				{
					int success=0;
					double
						*mat1=matrix+(OLS6_NPARAMS2*OLS6_NPARAMS2<<1)*kc,
						*mat2=matrix2+(OLS6_NPARAMS2*OLS6_NPARAMS2<<1)*kc;
					//if(kx==OLS6_SAMPLEREACH+OLS6_REACH)//initialize row
					if(!(kx%OLS6_STEP))
					{
#if 0
						//params = ((inv(sum i: sample[i][x] o sample[i][y]) * (sample[i] .* <target[i], ...>)
						for(int ky2=0;ky2<nparams;++ky2)
							memcpy(mat2+nparams*ky2, mat1+nparams*2*ky2, sizeof(double)*nparams);
						success=solve_Mx_v_cholesky(mat2, vec2, nparams, params);
#endif
#if 1
						memset(mat1, 0, sizeof(double[2])*nparams*nparams);
						for(int k=0;k<nparams;++k)
						{
							mat1[(nparams<<1|1)*k]=0.0005;//add small lambda for regularization
							mat1[(nparams<<1|1)*k+nparams]=1;//identity
						}
						for(int ky2=-OLS6_SAMPLEREACH, idx2=0;ky2<=0;++ky2)
						{
							for(int kx2=-OLS6_SAMPLEREACH;kx2<=OLS6_SAMPLEREACH+OLS6_STEP-1&&(ky2||kx2);++kx2, ++idx2)
								ols6_add_sample(pixels, gains[kc], src->iw, kc, kx+kx2, ky+ky2, csamples+(OLS6_NPARAMS2+1)*idx2, mat1, 0);
						}

						memcpy(mat2, mat1, sizeof(double[2])*nparams*nparams);
						success=invert_matrix(mat2, nparams, temp);
#endif
					}
#if 0
					else if(!initialized[kc]||!(kx%OLS6_STEP))//INCOMPATIBLE WITH STEP>1
					{
						int idx2;
						for(int ky2=0;ky2<OLS6_SAMPLEREACH-1;++ky2)
						{
							idx2=(2*OLS6_SAMPLEREACH+OLS6_STEP)*ky2+(kx-OLS6_SAMPLEREACH-1)%(2*OLS6_SAMPLEREACH+OLS6_STEP);
							ols6_add_sample(pixels, gains[kc], src->iw, kc, kx-OLS6_SAMPLEREACH-1, ky-OLS6_SAMPLEREACH+ky2, csamples+(OLS6_NPARAMS2+1)*idx2, mat1, 1);
							ols6_add_sample(pixels, gains[kc], src->iw, kc, kx+OLS6_SAMPLEREACH-1, ky-OLS6_SAMPLEREACH+ky2, csamples+(OLS6_NPARAMS2+1)*idx2, mat1, 0);
						}
						idx2=(2*OLS6_SAMPLEREACH+OLS6_STEP)*(OLS6_SAMPLEREACH-1)+(kx-OLS6_SAMPLEREACH-1)%OLS6_SAMPLEREACH;
						ols6_add_sample(pixels, gains[kc], src->iw, kc, kx-OLS6_SAMPLEREACH, ky, csamples+(OLS6_NPARAMS2+1)*idx2, mat1, 1);
						ols6_add_sample(pixels, gains[kc], src->iw, kc, kx-1, ky, csamples+(OLS6_NPARAMS2+1)*idx2, mat1, 0);

						memcpy(mat2, mat1, sizeof(double[2])*nparams*nparams);
						success=invert_matrix(mat2, nparams, temp);
					}
#endif
					//NNWW NNW NN NNE NNEE
					//NWW  NW  N  NE  NEE
					//WW   W   ?
					if(!success&&!initialized[kc])
						pred[0]=ols6_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
					else
					{
						if(success)
						{
							//params = ((inv(NB * NBT) * NB * Targets) - params)*lr
							for(int ky2=0;ky2<nparams;++ky2)//temp = NB * Targets
							{
								double sum=0;
								for(int kx2=0;kx2<OLS6_NSAMPLES;++kx2)
									sum+=csamples[(OLS6_NPARAMS2+1)*kx2+ky2]*csamples[(OLS6_NPARAMS2+1)*kx2+nparams];
								temp[ky2]=sum;
							}
							for(int ky2=0;ky2<nparams;++ky2)//params = matrix * temp
							{
								double sum=0;
								for(int j=0;j<nparams;++j)
									sum+=mat2[(nparams<<1)*ky2+j+nparams]*temp[j];
								cparamsN[ky2]=(9*cparamsW[ky2]+5*sum-cparamsNW[ky2]+cparamsNE[ky2]+cparamsNEE[ky2]+cparamsNEEE[ky2])/16;
								//cparamsN[ky2]=(6*cparamsW[ky2]+4*sum-cparamsN[ky2]+cparamsNEEE[ky2])/10;
								//cparamsN[ky2]=(6*cparamsW[ky2]+4*sum-cparamsN[ky2]+cparamsNE[ky2])/10;
								//cparamsN[ky2]=(6*cparamsW[ky2]+4*sum)/10;
								//cparamsN[ky2]=(cparamsW[ky2]+2*sum+cparamsNEEE[ky2])*0.25;
								//cparamsN[ky2]=(cparamsW[ky2]+2*sum+cparamsNEE[ky2])*0.25;
								//if(initialized[kc])
								//	cparams[ky2]+=(sum-cparams[ky2])*0.4;
								//else
								//	cparams[ky2]=sum;
							//	avparams[OLS6_NPARAMS2*kc+ky2]+=cparams[ky2];
							}
							initialized[kc]=1;
						}
						{
							double *cparams=cparamsN;
							int nb[14], j=0;
							//pred = nb * params
							int
								N =pixels[(idx-src->iw)<<2|kc],
								W =pixels[(idx-1)<<2|kc],
								NE=pixels[(idx-src->iw+1)<<2|kc];
							double fpred=cparams[0]/gains[kc];
							for(int ky3=-OLS6_REACH, idx2=1;ky3<=0;++ky3)
							{
								for(int kx3=-OLS6_REACH;kx3<=OLS6_REACH;++kx3)
								{
									for(int kc2=0;kc2<3;++kc2, ++idx2)
									{
										if(!ky3&&!kx3&&kc2==kc)//exclude current pixel
											goto finish;
										fpred+=pixels[(src->iw*(ky+ky3)+kx+kx3)<<2|kc2]*cparams[idx2];
										nb[j++]=pixels[(src->iw*(ky+ky3)+kx+kx3)<<2|kc2];
									}
								}
							}
						finish:
							{
								__m128i mpred=_mm_cvtpd_epi32(_mm_set_sd(fpred));
								__m128i mN	=_mm_set_epi32(0, 0, 0, N);
								__m128i mW	=_mm_set_epi32(0, 0, 0, W);
								__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
								__m128i vmin=_mm_min_epi32(mN, mW);
								__m128i vmax=_mm_max_epi32(mN, mW);
								vmin=_mm_min_epi32(vmin, mNE);
								vmax=_mm_max_epi32(vmax, mNE);
								mpred=_mm_max_epi32(mpred, vmin);
								mpred=_mm_min_epi32(mpred, vmax);
								_mm_store_si128((__m128i*)pred, mpred);
							}
							++successcount;
#ifdef ENABLE_LOGGING
							if(!kc)
							{
								console_log("YXC (%d, %d, %d): pred %d  curr %d  error %d\n", ky, kx, kc, pred[0], src->data[idx<<2|kc], src->data[idx<<2|kc]-pred[0]);
								//for(int k=0;k<nparams;++k)
								//	console_log("%lf,%c", cparams[k], k&3?' ':'\n');
								//console_log("\n");
								console_log("params:\n");
								for(int k=0;k<nparams;++k)
								{
									int val=(int)round(cparams[k]*256);
									console_log("%3d,%c", val, k%3?' ':'\n');
									//console_log("%c0x%02X,%c", val<0?'-':'+', abs(val), k%3?' ':'\n');
								}
								console_log("\n");
								console_log("NB:\n");
								for(int k=0;k<nparams-1;++k)
								{
									int val=nb[k];
									console_log("%3d,%c", val, (k+1)%3?' ':'\n');
									//console_log("%c0x%02X,%c", val<0?'-':'+', abs(val), (k+1)%3?' ':'\n');
								}
								console_log("\n");
							}
							if(kc==2)
							{
								console_pause();
								console_log("\n");
							}
#endif
						}
					}
				}
				{
					int curr=pred[0];
					curr^=-fwd;
					curr+=fwd;
					curr+=src->data[idx<<2|kc];
					curr<<=32-src->depth[kc];
					curr>>=32-src->depth[kc];
					src->data[idx<<2|kc]=curr;
				}
			}
			paramidx=(paramidx+1)%(src->iw+1);
		}
		if(loud_transforms)
			set_window_title("%d/%d = %7.3lf%%  OLS3 rate %lf%%  %lf sec", ky+1, src->ih, 100.*(ky+1)/src->ih, 100.*successcount/(src->nch*src->iw*(ky+1)), time_sec()-t_start);
	}
	//memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	//if(loud_transforms)
	//{
	//	set_window_title("OLS3 %lf%%  %lf sec", 100.*successcount/(src->nch*src->iw*src->ih), time_sec()-t_start);
	//	for(int k=0;k<(int)_countof(avparams);++k)
	//		avparams[k]*=(double)src->nch/successcount;
	//	int printed=0;
	//	printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1, "NW N NE W curr\n");
	//	for(int kg=0;kg<3;++kg)
	//	{
	//		for(int ky=0;ky<3;++ky)
	//		{
	//			for(int kx=0;kx<5;++kx)
	//				printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1, " %+12.8lf", avparams[(kg*3+ky)*5+kx]);
	//			printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1, "\n");
	//		}
	//		printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1, "\n");
	//	}
	//	copy_to_clipboard(g_buf, printed);
	//	messagebox(MBOX_OK, "Copied to clipboard", g_buf);
	//}
	free(dst);
	free(samples);
	free(matrix);
	free(matrix2);
	free(temp);
	free(params);
#endif
}