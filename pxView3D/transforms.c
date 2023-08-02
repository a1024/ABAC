#include"pxview3d.h"
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include"intercept_malloc.h"
static const char file[]=__FILE__;

void addhalf(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=128;
	}
}

void colortransform_xgz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xgz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xyz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		//r-=g;
		//g+=r>>1;
		//b-=g;
		//g+=b>>1;

		r-=b;
		g-=b;
		b+=(r+g)>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;

		//buf[k  ]=r-g+128;//XGZ
		//buf[k|1]=g;
		//buf[k|2]=b-g+128;

		//buf[k  ]=r;//RYZ
		//buf[k|1]=g-r+128;
		//buf[k|2]=b-r+128;

		//buf[k  ]=r-b+128;//XYBdash
		//buf[k|1]=g-b+128;
		//buf[k|2]=b;
	}
}
void colortransform_xyz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		//g-=b>>1;
		//b+=g;
		//g-=r>>1;
		//r+=g;

		b-=(r+g)>>1;
		g+=b;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=b;		//co = r-b			diff(r, b)
		b+=r>>1;	//(r+b)/2
		g-=b;		//cg = g-(r+b)/2	diff(g, av(r, b))
		b+=g>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		b-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocb_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;		//diff(r, g)
		g+=r>>1;
		b-=g;		//diff(b, av(r, g))
		g+=b>>1;	//av(b, av(r, g))

		buf[k  ]=r;//Co
		buf[k|1]=g;//Y
		buf[k|2]=b;//Cb
	}
}
void colortransform_ycocb_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}

void c_mul_add(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]+=r, dst[1]+=i;
}
void c_mul_sub(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]-=r, dst[1]-=i;
}
double c_abs2(const double *z)
{
	return z[0]*z[0]+z[1]*z[1];
}
void c_div(double *dst, double const *a, double const *b)
{
	double invabsb2=1/c_abs2(b);
	double
		r=(a[0]*b[0]+a[1]*b[1])*invabsb2,
		i=(a[1]*b[0]-a[0]*b[1])*invabsb2;
	dst[0]=r, dst[1]=i;
}
void c_exp(double *dst, double const *x)
{
	double m=exp(x[0]);
	dst[0]=m*cos(x[1]);
	dst[1]=m*sin(x[1]);
}
void c_ln(double *dst, double const *x)
{
	double
		r=log(sqrt(c_abs2(x))),
		i=atan2(x[1], x[0]);
	dst[0]=r, dst[1]=i;
}
void c_cbrt(double *dst, const double *x)//sqrt(x)=exp(0.5lnx)
{
	double temp[2];
	if(x[0]||x[1])
	{
		c_ln(temp, x);
		temp[0]*=1./3, temp[1]*=1./3;
		c_exp(dst, temp);
	}
	else
		dst[0]=x[0], dst[1]=x[1];
}
void impl_solve_cubic(const double *coeffs, double *roots)//finds r[0], r[1], & r[2], the solutions of x^3 + c[2]x^2 + c[1]x + c[0] = 0
{
	//https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
	double p=coeffs[2], q=coeffs[1], r=coeffs[0],
		p2, p3, q2, pq, A[2], B[2];
	double sqrt27=3*sqrt(3), inv3cbrt2=1./(3*cbrt(2)), ninth=1./9;
	double cm[]={-0.5, -sqrt(3)*0.5}, cp[]={-0.5, sqrt(3)*0.5};

	if(!p&&!q)
	{
		roots[4]=roots[2]=roots[0]=-cbrt(r);
		roots[5]=roots[3]=roots[1]=0;
		return;
	}

	p2=p*p;
	p3=p2*p;
	q2=q*q;
	pq=p*q;

	A[0]=(4*q-p2)*q2+(4*p3-18*pq+27*r)*r;
	if(A[0]<0)
		A[1]=sqrt(fabs(A[0])), A[0]=0;
	else
		A[0]=sqrt(A[0]), A[1]=0;
	A[0]*=sqrt27;
	A[1]*=sqrt27;
	A[0]-=27*r;
	A[0]+=9*pq-2*p3;
	c_cbrt(A, A);
	A[0]*=inv3cbrt2;
	A[1]*=inv3cbrt2;
	B[0]=3*q-p2;
	c_div(B, B, A);
	B[0]*=ninth;
	B[1]*=ninth;
	roots[4]=roots[2]=roots[0]=p*(-1./3);
	roots[5]=roots[3]=roots[1]=0;
	roots[0]+=A[0]-B[0];
	roots[1]+=A[1]-B[1];
	c_mul_add(roots+2, A, cm);
	c_mul_sub(roots+2, B, cp);
	c_mul_add(roots+4, A, cp);
	c_mul_sub(roots+4, B, cm);
}
void impl_rref(double *m, short dx, short dy)
{
#ifdef _DEBUG
	double pivot;
#endif
	double coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(fabs(m[dx*ky+it])>1e-10)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot row
			{
				coeff=1/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*kpivot+kx]*=coeff;
			}
			else//subtract pivot row from all other rows
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			//print_matrix_debug(m, dx, dy);//
		}
	}
}
int impl_nullspace(double *M, int dx, int dy, double *solution, int sstart, char *dep_flags, short *row_idx)
{
	//M is rref'ed
	//solution allocated size dy*dy, pre-memset solution to zero
	//p_flags & row_idx both size dx
	//returns number of vectors in solution
	int kx, kxdst, ky, keq, kfree, idx, idx2, nvec;
#ifdef DEBUG_NULLSPACE
	printf("Before RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	impl_rref(M, dx, dy);//first, get rref(M)
#ifdef DEBUG_NULLSPACE
	printf("After RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	memset(dep_flags, 0, dx);
	memset(row_idx, 0, dx*sizeof(short));
	for(ky=0;ky<dy;++ky)//find pivots (dependent variables)
	{
		for(kx=ky;kx<dx;++kx)//find first nonzero element
		{
			idx=dx*ky+kx;
			if(fabs(M[idx])>1e-10)
				break;
		}
		if(kx<dx)//if found
			dep_flags[kx]=1, row_idx[ky]=kx;
		else
			break;
	}
	nvec=dx-ky;
	for(ky=0, keq=0, kfree=0;ky<dx;++ky)
	{
		if(dep_flags[ky])//pivot, dependent variable
		{
			for(kx=0, kxdst=0;kx<dx;++kx)
			{
				if(dep_flags[kx])
					continue;
				idx=dx*ky+kxdst, idx2=dx*keq+kx;
				if(sstart+idx>=9)
					LOG_ERROR("");
				solution[sstart+idx]=-M[idx2];
				++kxdst;
			}
			++keq;
		}
		else//free variable
		{
			idx=dx*ky;
			if(sstart+idx+nvec>9)
				LOG_ERROR("");
			memset(solution+sstart+idx, 0, nvec*sizeof(double));
			if(sstart+idx+kfree>=9)
				LOG_ERROR("");
			solution[sstart+idx+kfree]=1;
			++kfree;
		}
#ifdef DEBUG_NULLSPACE
		//printf("Nullspace row %d:\n", ky);
		print_matrix_debug(solution+dx*ky, nvec, 1);//
#endif
	}
	return nvec;
}
int impl_egvec(double const *M, int n, const double *lambdas, double *S)
{
	int kv, kx, nvec, size=n*n, again;
	double *temp=(double*)malloc(size*sizeof(double));
	char *dep_flags=(char*)malloc(n);
	short *row_idx=(short*)malloc(n*sizeof(short));
	if(!temp||!dep_flags||!row_idx)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(S, 0, size*sizeof(double));
	for(kv=0, nvec=0;kv<n&&nvec<n;++kv)//get eigenvectors
	{
		again=0;
		for(kx=0;kx<kv;++kx)//check for repeated eigenvalues
		{
			if(fabs(lambdas[kx]-lambdas[kv])<1e-10)
			{
				again=1;
				break;
			}
		}
		if(again)
			continue;
		memcpy(temp, M, size*sizeof(double));
		for(kx=0;kx<n;++kx)
			temp[(n+1)*kx]-=lambdas[kv];

		nvec+=impl_nullspace(temp, n, n, S, nvec, dep_flags, row_idx);
		//if(nvec>6)
		//	LOG_ERROR("OOB");
		//print_matrix_debug(S, n, n);//
		//for(ky=0;ky<n;++ky)
		//{
		//	for(kx=ky;kx<n;++kx)
		//		if(fabs(temp[n*ky+kx].r)>1e-10||fabs(temp[n*ky+kx].i)>1e-10)
		//			break;
		//}
	}
	free(dep_flags), free(row_idx);
	free(temp);
	return nvec;
}
void colortransform_adaptive(char *buf, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int black=0xFF000000;
	memfill(b2, &black, (size_t)res<<2, 4);
	char *src=fwd?buf:b2;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
#define ACT_N 12
			char r[ACT_N]={0}, g[ACT_N]={0}, b[ACT_N]={0};//topleft, top, topright, left
			if(ky-2>=0)
			{
				if(kx-2>=0)	r[0]=src[(idx-iw*2-2)<<2], g[0]=src[(idx-iw*2-2)<<2|1], b[0]=src[(idx-iw*2-2)<<2|2];
				if(kx-1>=0)	r[1]=src[(idx-iw*2-1)<<2], g[1]=src[(idx-iw*2-1)<<2|1], b[1]=src[(idx-iw*2-1)<<2|2];
							r[2]=src[(idx-iw*2  )<<2], g[2]=src[(idx-iw*2  )<<2|1], b[2]=src[(idx-iw*2  )<<2|2];
				if(kx+1<iw)	r[3]=src[(idx-iw*2+1)<<2], g[3]=src[(idx-iw*2+1)<<2|1], b[3]=src[(idx-iw*2+1)<<2|2];
				if(kx+2<iw)	r[4]=src[(idx-iw*2+2)<<2], g[4]=src[(idx-iw*2+2)<<2|1], b[4]=src[(idx-iw*2+2)<<2|2];
			}
			if(ky-1>=0)
			{
				if(kx-2>=0)	r[5]=src[(idx-iw-2)<<2], g[5]=src[(idx-iw-2)<<2|1], b[5]=src[(idx-iw-2)<<2|2];
				if(kx-1>=0)	r[6]=src[(idx-iw-1)<<2], g[6]=src[(idx-iw-1)<<2|1], b[6]=src[(idx-iw-1)<<2|2];
							r[7]=src[(idx-iw  )<<2], g[7]=src[(idx-iw  )<<2|1], b[7]=src[(idx-iw  )<<2|2];
				if(kx+1<iw)	r[8]=src[(idx-iw+1)<<2], g[8]=src[(idx-iw+1)<<2|1], b[8]=src[(idx-iw+1)<<2|2];
				if(kx+2<iw)	r[9]=src[(idx-iw+2)<<2], g[9]=src[(idx-iw+2)<<2|1], b[9]=src[(idx-iw+2)<<2|2];
			}
			if(kx-2>=0)	r[10]=src[(idx-2)<<2], g[10]=src[(idx-2)<<2|1], b[10]=src[(idx-2)<<2|2];
			if(kx-1>=0)	r[11]=src[(idx-1)<<2], g[11]=src[(idx-1)<<2|1], b[11]=src[(idx-1)<<2|2];

			//if(kx-1>=0&&ky-1>=0)
			//	r[0]=src[(idx-iw-1)<<2], g[0]=src[(idx-iw-1)<<2|1], b[0]=src[(idx-iw-1)<<2|2];
			//if(ky-1>=0)
			//	r[1]=src[(idx-iw)<<2], g[1]=src[(idx-iw)<<2|1], b[1]=src[(idx-iw)<<2|2];
			//if(kx+1<iw&&ky-1>=0)
			//	r[2]=src[(idx-iw+1)<<2], g[2]=src[(idx-iw+1)<<2|1], b[2]=src[(idx-iw+1)<<2|2];
			//if(kx-1>=0)
			//	r[3]=src[(idx-1)<<2], g[3]=src[(idx-1)<<2|1], b[3]=src[(idx-1)<<2|2];
			
			//get prediction range
			int rmin=r[0], bmin=b[0];
			int rmax=r[0], bmax=b[0];

			for(int k=1;k<ACT_N;++k)
			{
				if(rmin>r[k]) rmin=r[k];
				if(rmax<r[k]) rmax=r[k];
				//if(gmin>g[k]) gmin=g[k];
				//if(gmax<g[k]) gmax=g[k];
				if(bmin>b[k]) bmin=b[k];
				if(bmax<b[k]) bmax=b[k];
			}

			char rcurr=buf[idx<<2], gcurr=buf[idx<<2|1], bcurr=buf[idx<<2|2];

#if 0
			int pred=((g[3]?(r[3]*gcurr+(g[3]>>1))/g[3]:r[3])+(b[3]?(r[3]*bcurr+(b[3]>>1))/b[3]:r[3])+1)>>1;
			pred=CLAMP(rmin, pred, rmax);
			rcurr-=pred;

			pred=g[3]?(b[3]*gcurr+(g[3]>>1))/g[3]:b[3];
			pred=CLAMP(bmin, pred, bmax);
			bcurr-=pred;
#endif

#if 1
			//if(kx==405&&ky==170)
			//if(kx==4&&ky==2)
			//	kx=4;

			double mean[3]={0};
			double points[ACT_N*3];
			for(int k=0;k<ACT_N;++k)
			{
				mean[0]+=r[k];
				mean[1]+=g[k];
				mean[2]+=b[k];
			}
			mean[0]/=ACT_N;
			mean[1]/=ACT_N;
			mean[2]/=ACT_N;
			for(int k=0;k<ACT_N;++k)
			{
				points[3*k  ]=r[k]-mean[0];
				points[3*k+1]=g[k]-mean[1];
				points[3*k+2]=b[k]-mean[2];
			}
			double cov[9]={0};
			for(int k=0;k<4;++k)
			{
				double *point=points+k*3, xy=point[0]*point[1], yz=point[1]*point[2], zx=point[2]*point[0];
				cov[0]+=point[0]*point[0]; cov[1]+=xy; cov[2]+=zx;//outer product
				cov[3]+=xy; cov[4]+=point[1]*point[1]; cov[5]+=yz;
				cov[6]+=zx; cov[7]+=yz; cov[8]+=point[2]*point[2];
			}
			for(int k=0;k<9;++k)//divide matrix by N-1 because the mean was from these samples
				cov[k]*=1./3;

			//get eigenvalues
			double evalues[6]={0};
			double C[3];
				
			//C[2] = -tr(M) = -(M[0]+M[4]+M[8]);
			C[2]=-(cov[0]+cov[4]+cov[8]);

			//C[1] = cof0+cof4+cof8 = M[4]*M[8]-M[5]*M[7] + M[0]*M[8]-M[2]*M[6] + M[0]*M[4]-M[1]*M[3];
			//C[1]=cov[4]*cov[8]-cov[5]*cov[7];
			//C[0]=C[1];
			//C[1]+=cov[0]*cov[8];
			//C[1]-=cov[2]*cov[6];
			//C[1]+=cov[0]*cov[4];
			//C[1]-=cov[1]*cov[3];
			C[1]=cov[4]*cov[8]-cov[5]*cov[7] + cov[0]*cov[8]-cov[2]*cov[6] + cov[0]*cov[4]-cov[1]*cov[3];

			//C[0] = -det(M) = -(M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6]));
			//C[0]*=cov[0];
			//C[0]-=cov[1]*(cov[3]*cov[8]-cov[5]*cov[6]);
			//C[0]+=cov[2]*(cov[3]*cov[7]-cov[4]*cov[6]);
			//C[0]=-C[0];
			C[0]=-(cov[0]*(cov[4]*cov[8]-cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8]-cov[5]*cov[6]) + cov[2]*(cov[3]*cov[7]-cov[4]*cov[6]));

			//if(!C[0]&&!C[1]&&!C[2])
			//	goto skip;

			impl_solve_cubic(C, evalues);

			evalues[1]=evalues[2];//take real parts
			evalues[2]=evalues[4];
			//evalues[0]=sqrt(evalues[0]*evalues[0]+evalues[1]*evalues[1]);
			//evalues[1]=sqrt(evalues[2]*evalues[2]+evalues[3]*evalues[3]);
			//evalues[2]=sqrt(evalues[4]*evalues[4]+evalues[5]*evalues[5]);

			if(fabs(evalues[0])<fabs(evalues[1]))SWAPVAR(evalues[0], evalues[1], evalues[3]);//sort in descending order
			if(fabs(evalues[1])<fabs(evalues[2]))SWAPVAR(evalues[1], evalues[2], evalues[3]);
			if(fabs(evalues[0])<fabs(evalues[1]))SWAPVAR(evalues[0], evalues[1], evalues[3]);

			//get eigenvectors
			double evecs[9];
			impl_egvec(cov, 3, evalues, evecs);

			//normalize v0
			evalues[3]=1/sqrt(evecs[0]*evecs[0]+evecs[1]*evecs[1]+evecs[2]*evecs[2]);
			evecs[0]*=evalues[3];
			evecs[1]*=evalues[3];
			evecs[2]*=evalues[3];

			//predict components
			int pred;
			int count=0;
			double t=0;
			if(fwd)
			{
				if(evecs[1])
					t+=(gcurr-mean[1])/evecs[1], ++count;
				if(evecs[2])
					t+=(bcurr-mean[2])/evecs[2], ++count;
				if(count)
				{
					t/=count;
					pred=(int)(t*evecs[0]+mean[0]);
					pred=CLAMP(rmin, pred, rmax);
				}
				else
					pred=(gcurr+bcurr)>>1;
				rcurr-=pred;

				if(evecs[1])
				{
					t=(gcurr-mean[1])/evecs[1];
					pred=(int)(t*evecs[2]+mean[2]);
					pred=CLAMP(bmin, pred, bmax);
				}
				else
					pred=gcurr;
				bcurr-=pred;
			}
			else
			{
				if(evecs[1])
				{
					t=(gcurr-mean[1])/evecs[1];
					pred=(int)(t*evecs[2]+mean[2]);
					pred=CLAMP(bmin, pred, bmax);
				}
				else
					pred=gcurr;
				bcurr+=pred;
				
				if(evecs[1])
					t+=(gcurr-mean[1])/evecs[1], ++count;
				if(evecs[2])
					t+=(bcurr-mean[2])/evecs[2], ++count;
				if(count)
				{
					t/=count;
					pred=(int)(t*evecs[0]+mean[0]);
				}
				else
					pred=(gcurr+bcurr)>>1;
				pred=CLAMP(rmin, pred, rmax);
				rcurr+=pred;
			}
#endif
#undef ACT_N

#if 0
			//if(red&&green&&blue)
			//{
			//	int LOL_1=0;
			//}
			double coeff[3]={0};
			
			coeff[0]=green+blue?(double)red/(green+blue):1;
			coeff[0]=CLAMP(-128, coeff[0], 127);
			red-=(int)(coeff[0]*(green+blue));
			red=CLAMP(-128, red, 127);

			coeff[1]=red+green?(double)blue/(red+green):1;
			coeff[1]=CLAMP(-128, coeff[1], 127);
			blue-=(int)(coeff[1]*(red+green));
			blue=CLAMP(-128, blue, 127);

			coeff[2]=red+blue?(double)green/(red+blue):1;
			coeff[2]=CLAMP(-128, coeff[2], 127);
			//green+=(int)(coeff[2]*(red+blue));
			//green=CLAMP(-128, green, 127);

			int temp;
			if(fwd)
			{
				temp=(int)(coeff[0]*(cg+cb));
				cr-=CLAMP(-128, temp, 127);

				temp=(int)(coeff[1]*(cr+cg));
				cb-=CLAMP(-128, temp, 127);

				temp=(int)(coeff[2]*(cr+cb));
				cg-=CLAMP(-128, temp, 127);
			}
			else
			{
				temp=(int)(coeff[2]*(cr+cb));
				cg+=CLAMP(-128, temp, 127);

				temp=(int)(coeff[1]*(cr+cg));
				cb+=CLAMP(-128, temp, 127);

				temp=(int)(coeff[0]*(cg+cb));
				cr+=CLAMP(-128, temp, 127);
			}
#endif

			b2[idx<<2  ]=rcurr;
			b2[idx<<2|1]=gcurr;
			b2[idx<<2|2]=bcurr;
		}
	}
	memcpy(buf, b2, (size_t)res<<2);
	free(b2);
}

void colortransform_exp_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char
				r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2], r2=0, g2=0, b2=0, r3=0, g3=0, b3=0;
			if(kx+1<iw)//right
				r2=buf[(idx+1)<<2], g2=buf[(idx+1)<<2|1], b2=buf[(idx+1)<<2|2];
			if(ky+1<ih)//below
				r3=buf[(idx+iw)<<2], g3=buf[(idx+iw)<<2|1], b3=buf[(idx+iw)<<2|2];

			if(abs(r2-g2)<abs(r2-b2)&&abs(r3-g3)<abs(r3-b3))
			{
				r-=g;
				g+=r>>1;
			}
			else if(abs(r2-g2)>abs(r2-b2)&&abs(r3-g3)>abs(r3-b3))
			{
				r-=b;
				b+=r>>1;
			}

			//if(abs(r2-g2-(r3-g3))<abs(r2-r3))
			//	r-=g;
			//else
			//	r-=(r2+r3)>>1;

			//r-=(g+g+r2+r3)>>2;

			//if(abs(r2-r3)<abs(g-b))
			//	r-=(r2+r3)>>1;
			//else
			//	r-=g;

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}
void colortransform_exp_inv(char *buf, int iw, int ih)
{
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx)
		{
			int idx=iw*ky+kx;
			char
				r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2], r2=0, g2=0, b2=0, r3=0, g3=0, b3=0;
			if(kx+1<iw)//right
				r2=buf[(idx+1)<<2], g2=buf[(idx+1)<<2|1], b2=buf[(idx+1)<<2|2];
			if(ky+1<ih)//below
				r3=buf[(idx+iw)<<2], g3=buf[(idx+iw)<<2|1], b3=buf[(idx+iw)<<2|2];
			
			if(abs(r2-g2)<abs(r2-b2)&&abs(r3-g3)<abs(r3-b3))
			{
				r+=g;
				g-=r>>1;
			}
			else if(abs(r2-g2)>abs(r2-b2)&&abs(r3-g3)>abs(r3-b3))
			{
				r+=b;
				b-=r>>1;
			}

			//if(abs(r2-g2-(r3-g3))<abs(r2-r3))
			//	r+=g;
			//else
			//	r+=(r2+r3)>>1;
			
			//r+=(g+g+r2+r3)>>2;

			//if(abs(r2-r3)<abs(g-b))
			//	r-=(r2+r3)>>1;
			//else
			//	r-=g;

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}

#if 0
float lrct[]=
{
	-0.235968455672264100f, -0.005472748540341854f,  0.546232640743255600f,
	 0.042025841772556305f, -0.281145840883255000f, -0.041713438928127290f,
	-0.270333796739578250f,  0.382541775703430200f,  0.122593350708484650f,
	-0.654174447059631300f, -0.111721098423004150f,  0.049731694161891940f,
	-0.180651888251304630f, -0.665180742740631100f, -0.449275195598602300f,
	-0.634068489074707000f, -0.438444226980209350f, -0.079102315008640290f,
	-0.492376208305358900f, -0.011669249273836613f,  0.612039744853973400f,
	-0.140227079391479500f,  0.255208671092987060f, -0.405863404273986800f,
	-0.133800372481346130f, -0.049969632178545000f, -0.006972887087613344f,
	-0.498686790466308600f, -0.417643666267395000f,  0.129942402243614200f,
	-0.584400236606597900f, -0.273861020803451540f,  0.104215130209922790f,
	 0.297003507614135740f, -0.688771069049835200f,  0.330031186342239400f,
};
char lift(char v1, char v2, float *coeff)
{
	return (char)(127*(v1*coeff[0]/127+v2*coeff[1]/127+coeff[2]));
}
void colortransform_learned_fwd(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		r-=lift(g, b, lrct+3* 0);
		g-=lift(r, b, lrct+3* 1);
		b-=lift(r, g, lrct+3* 2);
		r+=lift(g, b, lrct+3* 3);
		g+=lift(r, b, lrct+3* 4);
		b+=lift(r, g, lrct+3* 5);
		r-=lift(g, b, lrct+3* 6);
		g-=lift(r, b, lrct+3* 7);
		b-=lift(r, g, lrct+3* 8);
		r+=lift(g, b, lrct+3* 9);
		g+=lift(r, b, lrct+3*10);
		b+=lift(r, g, lrct+3*11);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
void colortransform_learned_inv(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		b-=lift(r, g, lrct+3*11);
		g-=lift(r, b, lrct+3*10);
		r-=lift(g, b, lrct+3* 9);
		b+=lift(r, g, lrct+3* 8);
		g+=lift(r, b, lrct+3* 7);
		r+=lift(g, b, lrct+3* 6);
		b-=lift(r, g, lrct+3* 5);
		g-=lift(r, b, lrct+3* 4);
		r-=lift(g, b, lrct+3* 3);
		b+=lift(r, g, lrct+3* 2);
		g+=lift(r, b, lrct+3* 1);
		r+=lift(g, b, lrct+3* 0);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
#endif
#if 0
#if 1
float lrt_biases[]=
{
	-0.007593977730721235f,
	-0.20463228225708008f,
	0.048144787549972534f,
	-0.052841395139694214f,
	-0.229848250746727f,
	-0.10854046046733856f,
	-0.27017149329185486f,
	-0.28025728464126587f,
	0.14076033234596252f,
	0.11415114998817444f,
	-0.24117304384708405f,
	-0.22909477353096008f,
};
float lrt_c01[]=
{
	 0.21996481716632843f,
	 0.3447851538658142f ,
	 0.4309721887111664f ,
	-0.3914514482021332f ,
	 0.1072743684053421f ,
	 0.09292115271091461f,
	 0.09527700394392014f,
	-0.7720031142234802f ,
	-0.03891817852854729f,
	 0.18422403931617737f,
	 0.6045548319816589f ,
};
float lrt_c02[]=
{
	0.16338559985160828f,
	-0.21537575125694275f,
	-0.01003316044807434f,
	0.0806010365486145f,
	0.1516939401626587f,
	0.22873657941818237f,
	-0.13795271515846252f,
	0.13615399599075317f,
	-0.07778285443782806f,
	0.291059672832489f,
	0.22079020738601685f,
};
float lrt_c03[]=
{
	0.010352730751037598f,
	0.1270458996295929f,
	-0.18215587735176086f,
	-0.164699524641037f,
	0.23722952604293823f,
	0.27438366413116455f,
	-0.05265103280544281f,
	0.18550443649291992f,
	-0.21034197509288788f,
	-0.13812671601772308f,
	-0.015580296516418457f,
};
float lrt_c04[]=
{
	-0.06830897927284241f,
	-0.0477299690246582f,
	-0.2479163110256195f,
	-0.04657846689224243f,
	-0.0802721232175827f,
	-0.004736065864562988f,
	-0.07743765413761139f,
	-0.26771998405456543f,
	0.013397663831710815f,
	-0.23979492485523224f,
	0.12400856614112854f,
};
float lrt_c05[]=
{
	-0.04351922869682312f,
	-0.11087171733379364f,
	-0.11081892251968384f,
	0.301426351070404050f,
	0.090425968170166020f,
	0.171547293663024900f,
	0.166505038738250730f,
	0.086036622524261470f,
	-0.25268214941024780f,
	0.168854206800460820f,
	0.221817135810852050f,
};
float lrt_c06[]=
{
	0.216941952705383300f,
	0.145165383815765380f,
	-0.30066022276878357f,
	0.062808305025100710f,
	-0.11871317028999329f,
	-0.03821510076522827f,
	0.057682424783706665f,
	-0.15345700085163116f,
	-0.10404434800148010f,
	0.159864872694015500f,
	0.030692487955093384f,
};
float lrt_c07[]=
{
	-0.11159592866897583f,
	0.267521142959594700f,
	-0.04237860441207886f,
	-0.12383817136287689f,
	0.221899330615997310f,
	0.031740218400955200f,
	-0.11875738203525543f,
	0.109959483146667480f,
	-0.17683275043964386f,
	-0.12659740447998047f,
	-0.21873277425765990f,
};
float lrt_c08[]=
{
	-0.18569713830947876f,
	0.224896728992462160f,
	-0.28754058480262756f,
	0.111558407545089720f,
	-0.18325410783290863f,
	-0.27818745374679565f,
	0.241845905780792240f,
	0.179279416799545300f,
	-0.05782620608806610f,
	-0.10470718145370483f,
	0.128715872764587400f,
};
float lrt_c09[]=
{
	0.0736332833766937300f,
	-0.278686881065368650f,
	-0.062273472547531130f,
	0.2923494577407837000f,
	0.1254649162292480500f,
	0.0480349659919738800f,
	0.2714382410049438500f,
	-0.016680717468261720f,
	0.0569746196269989000f,
	-0.021338820457458496f,
	0.0394396781921386700f,
};
float lrt_c10[]=
{
	0.068611204624176030f,
	-0.12615610659122467f,
	-0.07715198397636414f,
	0.252162516117095950f,
	-0.02322089672088623f,
	0.062368571758270264f,
	0.138979107141494750f,
	0.039735138416290280f,
	0.171300053596496580f,
	-0.13656683266162872f,
	-0.07565732300281525f,
};
float lrt_c11[]=
{
	0.1174532473087310800f,
	0.0668951570987701400f,
	-0.071375116705894470f,
	-0.240586921572685240f,
	0.1582727432250976600f,
	0.1294519901275634800f,
	0.2438962459564209000f,
	-0.014519184827804565f,
	0.1907848119735717800f,
	-0.215079009532928470f,
	0.0485148429870605500f,
};
float lrt_c12[]=
{
	0.086784452199935910f,
	0.276734590530395500f,
	-0.24387010931968690f,
	0.102736681699752810f,
	-0.28351089358329773f,
	0.072889894247055050f,
	-0.18153974413871765f,
	0.022794693708419800f,
	-0.27390283346176150f,
	-0.27018001675605774f,
	0.102353841066360470f,
};
#endif
void learnedtransform_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			//if(!kx&&!ky)
			//	kx=0;

			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x01-=(char)(127*(lrt_c01[0]*x02/127+lrt_c01[1]*x03/127+lrt_c01[2]*x04/127+lrt_c01[3]*x05/127+lrt_c01[4]*x06/127+lrt_c01[5]*x07/127+lrt_c01[6]*x08/127+lrt_c01[7]*x09/127+lrt_c01[8]*x10/127+lrt_c01[9]*x11/127+lrt_c01[10]*x12/127+lrt_biases[ 0]));
			x02-=(char)(127*(lrt_c02[0]*x01/127+lrt_c02[1]*x03/127+lrt_c02[2]*x04/127+lrt_c02[3]*x05/127+lrt_c02[4]*x06/127+lrt_c02[5]*x07/127+lrt_c02[6]*x08/127+lrt_c02[7]*x09/127+lrt_c02[8]*x10/127+lrt_c02[9]*x11/127+lrt_c02[10]*x12/127+lrt_biases[ 1]));
			x03-=(char)(127*(lrt_c03[0]*x01/127+lrt_c03[1]*x02/127+lrt_c03[2]*x04/127+lrt_c03[3]*x05/127+lrt_c03[4]*x06/127+lrt_c03[5]*x07/127+lrt_c03[6]*x08/127+lrt_c03[7]*x09/127+lrt_c03[8]*x10/127+lrt_c03[9]*x11/127+lrt_c03[10]*x12/127+lrt_biases[ 2]));
			x04-=(char)(127*(lrt_c04[0]*x01/127+lrt_c04[1]*x02/127+lrt_c04[2]*x03/127+lrt_c04[3]*x05/127+lrt_c04[4]*x06/127+lrt_c04[5]*x07/127+lrt_c04[6]*x08/127+lrt_c04[7]*x09/127+lrt_c04[8]*x10/127+lrt_c04[9]*x11/127+lrt_c04[10]*x12/127+lrt_biases[ 3]));
			x05-=(char)(127*(lrt_c05[0]*x01/127+lrt_c05[1]*x02/127+lrt_c05[2]*x03/127+lrt_c05[3]*x04/127+lrt_c05[4]*x06/127+lrt_c05[5]*x07/127+lrt_c05[6]*x08/127+lrt_c05[7]*x09/127+lrt_c05[8]*x10/127+lrt_c05[9]*x11/127+lrt_c05[10]*x12/127+lrt_biases[ 4]));
			x06-=(char)(127*(lrt_c06[0]*x01/127+lrt_c06[1]*x02/127+lrt_c06[2]*x03/127+lrt_c06[3]*x04/127+lrt_c06[4]*x05/127+lrt_c06[5]*x07/127+lrt_c06[6]*x08/127+lrt_c06[7]*x09/127+lrt_c06[8]*x10/127+lrt_c06[9]*x11/127+lrt_c06[10]*x12/127+lrt_biases[ 5]));
			x07-=(char)(127*(lrt_c07[0]*x01/127+lrt_c07[1]*x02/127+lrt_c07[2]*x03/127+lrt_c07[3]*x04/127+lrt_c07[4]*x05/127+lrt_c07[5]*x06/127+lrt_c07[6]*x08/127+lrt_c07[7]*x09/127+lrt_c07[8]*x10/127+lrt_c07[9]*x11/127+lrt_c07[10]*x12/127+lrt_biases[ 6]));
			x08-=(char)(127*(lrt_c08[0]*x01/127+lrt_c08[1]*x02/127+lrt_c08[2]*x03/127+lrt_c08[3]*x04/127+lrt_c08[4]*x05/127+lrt_c08[5]*x06/127+lrt_c08[6]*x07/127+lrt_c08[7]*x09/127+lrt_c08[8]*x10/127+lrt_c08[9]*x11/127+lrt_c08[10]*x12/127+lrt_biases[ 7]));
			x09-=(char)(127*(lrt_c09[0]*x01/127+lrt_c09[1]*x02/127+lrt_c09[2]*x03/127+lrt_c09[3]*x04/127+lrt_c09[4]*x05/127+lrt_c09[5]*x06/127+lrt_c09[6]*x07/127+lrt_c09[7]*x08/127+lrt_c09[8]*x10/127+lrt_c09[9]*x11/127+lrt_c09[10]*x12/127+lrt_biases[ 8]));
			x10-=(char)(127*(lrt_c10[0]*x01/127+lrt_c10[1]*x02/127+lrt_c10[2]*x03/127+lrt_c10[3]*x04/127+lrt_c10[4]*x05/127+lrt_c10[5]*x06/127+lrt_c10[6]*x07/127+lrt_c10[7]*x08/127+lrt_c10[8]*x09/127+lrt_c10[9]*x11/127+lrt_c10[10]*x12/127+lrt_biases[ 9]));
			x11-=(char)(127*(lrt_c11[0]*x01/127+lrt_c11[1]*x02/127+lrt_c11[2]*x03/127+lrt_c11[3]*x04/127+lrt_c11[4]*x05/127+lrt_c11[5]*x06/127+lrt_c11[6]*x07/127+lrt_c11[7]*x08/127+lrt_c11[8]*x09/127+lrt_c11[9]*x10/127+lrt_c11[10]*x12/127+lrt_biases[10]));
			x12-=(char)(127*(lrt_c12[0]*x01/127+lrt_c12[1]*x02/127+lrt_c12[2]*x03/127+lrt_c12[3]*x04/127+lrt_c12[4]*x05/127+lrt_c12[5]*x06/127+lrt_c12[6]*x07/127+lrt_c12[7]*x08/127+lrt_c12[8]*x09/127+lrt_c12[9]*x10/127+lrt_c12[10]*x11/127+lrt_biases[11]));
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04;
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08;
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_fwd(buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);
}
void learnedtransform_inv(char *buf, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_inv(buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);

	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			//if(!kx&&!ky)
			//	kx=0;

			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x12+=(char)(127*(lrt_c12[0]*x01/127+lrt_c12[1]*x02/127+lrt_c12[2]*x03/127+lrt_c12[3]*x04/127+lrt_c12[4]*x05/127+lrt_c12[5]*x06/127+lrt_c12[6]*x07/127+lrt_c12[7]*x08/127+lrt_c12[8]*x09/127+lrt_c12[9]*x10/127+lrt_c12[10]*x11/127+lrt_biases[11]));
			x11+=(char)(127*(lrt_c11[0]*x01/127+lrt_c11[1]*x02/127+lrt_c11[2]*x03/127+lrt_c11[3]*x04/127+lrt_c11[4]*x05/127+lrt_c11[5]*x06/127+lrt_c11[6]*x07/127+lrt_c11[7]*x08/127+lrt_c11[8]*x09/127+lrt_c11[9]*x10/127+lrt_c11[10]*x12/127+lrt_biases[10]));
			x10+=(char)(127*(lrt_c10[0]*x01/127+lrt_c10[1]*x02/127+lrt_c10[2]*x03/127+lrt_c10[3]*x04/127+lrt_c10[4]*x05/127+lrt_c10[5]*x06/127+lrt_c10[6]*x07/127+lrt_c10[7]*x08/127+lrt_c10[8]*x09/127+lrt_c10[9]*x11/127+lrt_c10[10]*x12/127+lrt_biases[ 9]));
			x09+=(char)(127*(lrt_c09[0]*x01/127+lrt_c09[1]*x02/127+lrt_c09[2]*x03/127+lrt_c09[3]*x04/127+lrt_c09[4]*x05/127+lrt_c09[5]*x06/127+lrt_c09[6]*x07/127+lrt_c09[7]*x08/127+lrt_c09[8]*x10/127+lrt_c09[9]*x11/127+lrt_c09[10]*x12/127+lrt_biases[ 8]));
			x08+=(char)(127*(lrt_c08[0]*x01/127+lrt_c08[1]*x02/127+lrt_c08[2]*x03/127+lrt_c08[3]*x04/127+lrt_c08[4]*x05/127+lrt_c08[5]*x06/127+lrt_c08[6]*x07/127+lrt_c08[7]*x09/127+lrt_c08[8]*x10/127+lrt_c08[9]*x11/127+lrt_c08[10]*x12/127+lrt_biases[ 7]));
			x07+=(char)(127*(lrt_c07[0]*x01/127+lrt_c07[1]*x02/127+lrt_c07[2]*x03/127+lrt_c07[3]*x04/127+lrt_c07[4]*x05/127+lrt_c07[5]*x06/127+lrt_c07[6]*x08/127+lrt_c07[7]*x09/127+lrt_c07[8]*x10/127+lrt_c07[9]*x11/127+lrt_c07[10]*x12/127+lrt_biases[ 6]));
			x06+=(char)(127*(lrt_c06[0]*x01/127+lrt_c06[1]*x02/127+lrt_c06[2]*x03/127+lrt_c06[3]*x04/127+lrt_c06[4]*x05/127+lrt_c06[5]*x07/127+lrt_c06[6]*x08/127+lrt_c06[7]*x09/127+lrt_c06[8]*x10/127+lrt_c06[9]*x11/127+lrt_c06[10]*x12/127+lrt_biases[ 5]));
			x05+=(char)(127*(lrt_c05[0]*x01/127+lrt_c05[1]*x02/127+lrt_c05[2]*x03/127+lrt_c05[3]*x04/127+lrt_c05[4]*x06/127+lrt_c05[5]*x07/127+lrt_c05[6]*x08/127+lrt_c05[7]*x09/127+lrt_c05[8]*x10/127+lrt_c05[9]*x11/127+lrt_c05[10]*x12/127+lrt_biases[ 4]));
			x04+=(char)(127*(lrt_c04[0]*x01/127+lrt_c04[1]*x02/127+lrt_c04[2]*x03/127+lrt_c04[3]*x05/127+lrt_c04[4]*x06/127+lrt_c04[5]*x07/127+lrt_c04[6]*x08/127+lrt_c04[7]*x09/127+lrt_c04[8]*x10/127+lrt_c04[9]*x11/127+lrt_c04[10]*x12/127+lrt_biases[ 3]));
			x03+=(char)(127*(lrt_c03[0]*x01/127+lrt_c03[1]*x02/127+lrt_c03[2]*x04/127+lrt_c03[3]*x05/127+lrt_c03[4]*x06/127+lrt_c03[5]*x07/127+lrt_c03[6]*x08/127+lrt_c03[7]*x09/127+lrt_c03[8]*x10/127+lrt_c03[9]*x11/127+lrt_c03[10]*x12/127+lrt_biases[ 2]));
			x02+=(char)(127*(lrt_c02[0]*x01/127+lrt_c02[1]*x03/127+lrt_c02[2]*x04/127+lrt_c02[3]*x05/127+lrt_c02[4]*x06/127+lrt_c02[5]*x07/127+lrt_c02[6]*x08/127+lrt_c02[7]*x09/127+lrt_c02[8]*x10/127+lrt_c02[9]*x11/127+lrt_c02[10]*x12/127+lrt_biases[ 1]));
			x01+=(char)(127*(lrt_c01[0]*x02/127+lrt_c01[1]*x03/127+lrt_c01[2]*x04/127+lrt_c01[3]*x05/127+lrt_c01[4]*x06/127+lrt_c01[5]*x07/127+lrt_c01[6]*x08/127+lrt_c01[7]*x09/127+lrt_c01[8]*x10/127+lrt_c01[9]*x11/127+lrt_c01[10]*x12/127+lrt_biases[ 0]));
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04;
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08;
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
}
#endif

const int customparam_ct_w=2, customparam_ct_h=6, customparam_st_reach=2;
/*const double customparam_ct0[]=
{
	 1,   0,
	-0.5, 0,
	 0,   1,
	 0,   0,
	 0,   1,
	 0,   0,
};
const double customparam_st0[]=
{
	0,  0, 0, 0, 0,
	0, -1, 1, 0, 0,
	0,  1,
};*/
int customparam_sel=12;
double customparam_ct[12]={0}, customparam_st[12*6]={0};
int customparam_ch_idx=0;
int customparam_clamp[2]={-128, 127};
void customtransforms_resetparams()
{
	memset(customparam_ct, 0, sizeof(customparam_ct));
	memset(customparam_st+12*customparam_ch_idx, 0, sizeof(customparam_st)/6);
	//memcpy(customparam_ct, customparam_ct0, sizeof(customparam_ct));
	//memcpy(customparam_st, customparam_st0, sizeof(customparam_st));
}
void colortransform_custom_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2];
			
			r-=(char)(customparam_ct[ 0]*g+customparam_ct[ 1]*b);
			g-=(char)(customparam_ct[ 2]*r+customparam_ct[ 3]*b);
			b-=(char)(customparam_ct[ 4]*r+customparam_ct[ 5]*g);
			r+=(char)(customparam_ct[ 6]*g+customparam_ct[ 7]*b);
			g+=(char)(customparam_ct[ 8]*r+customparam_ct[ 9]*b);
			b+=(char)(customparam_ct[10]*r+customparam_ct[11]*g);

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}
void colortransform_custom_inv(char *buf, int iw, int ih)
{
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx)
		{
			int idx=iw*ky+kx;
			char r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2];
			
			b-=(char)(customparam_ct[10]*r+customparam_ct[11]*g);
			g-=(char)(customparam_ct[ 8]*r+customparam_ct[ 9]*b);
			r-=(char)(customparam_ct[ 6]*g+customparam_ct[ 7]*b);
			b+=(char)(customparam_ct[ 4]*r+customparam_ct[ 5]*g);
			g+=(char)(customparam_ct[ 2]*r+customparam_ct[ 3]*b);
			r+=(char)(customparam_ct[ 0]*g+customparam_ct[ 1]*b);

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}


//chroma from luma (CfL)
short cfl_params[4]={0};//{alpha, beta} x2 chroma channels
void pred_cfl_calcparams(char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int kl, int kc, short *params)
{
	double alpha, beta, term[4]={0};
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			char L=buf[idx|kl], C=buf[idx|kc];
			term[0]+=L*C;
			term[1]+=L;
			term[2]+=C;
			term[3]+=L*L;
		}
	}
	int bw=x2-x1, bh=y2-y1, count=bw*bh;
	if(!count||!term[1]&&!term[3])
	{
		params[0]=0, params[1]=0;
		return;
	}
	alpha=(count*term[0]-term[1]*term[2])/(count*term[3]-term[1]*term[1]);
	beta=(term[2]-alpha*term[1])/count;
	params[0]=(short)(alpha*256);
	params[1]=(short)(beta*256);
}
void pred_cfl(char *buf, int iw, int ih, int fwd)
{
#if 0
	if(fwd)
	{
		pred_cfl_calcparams(buf, iw, ih, 0, iw, 0, ih, 1, 0, cfl_params);
		pred_cfl_calcparams(buf, iw, ih, 0, iw, 0, ih, 1, 2, cfl_params+2);
	}
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			char Co=buf[idx], Y=buf[idx|1], Cb=buf[idx|2];
			int pred_o=(cfl_params[0]*Y+cfl_params[1]+127)>>8;
			int pred_b=(cfl_params[2]*Y+cfl_params[3]+127)>>8;
			if(fwd)
			{
				Co-=pred_o;
				Cb-=pred_b;
			}
			else
			{
				Co+=pred_o;
				Cb+=pred_b;
			}
			buf[idx  ]=Co;
			buf[idx|2]=Cb;
		}
	}
#endif
}


//spatial transforms


const int permute_idx[]=
{
	136, 229, 247, 217,  82, 226,  99, 203, 139,  29, 238, 193, 141, 134,   7, 160,
	 32, 253, 165,  51,  37, 113, 115,  60,  94, 140,  78,  70, 202,  24, 218, 151,
	219, 190, 246, 161,  81, 163, 135, 248, 198, 220, 155, 102,  65,  35, 213, 122,
	242, 171,  57, 249,  12, 148, 119,  97, 232, 221, 132,  46, 117,  41,  84, 183,
	188,  93, 129, 157, 145,  90,  40, 243, 182, 154, 194,  91, 177,  72, 252, 121,
	 52, 101, 167,  92,  34,  38, 168, 235,  88, 215, 196, 211,  22,  63, 175, 199,
	 75,  16, 192,  87, 186,   2,  56, 245,  58, 137, 201,  61,  71, 234, 111,  76,
	170, 240,   1, 166,   4,  68,  43, 228, 208, 244,  15, 153,  53, 237, 123, 108,
	 95,  18,   8,  14,  23,  67,  62, 231,  74, 227,  27, 254,  33, 250, 207,  96,
	  0, 120, 173, 181,  39, 176, 225,  55, 107,  26, 251, 179, 239,  25, 112,  79,
	128,  64, 130, 162, 143,   9,  21,   6, 223, 216,  86, 191,  45, 144,   3, 106,
	 69, 147,  83, 114, 118, 169,  77,  31,  47,  42, 100, 174, 172, 138, 159, 164,
	 19, 209,  80, 131, 210,  89, 150, 236, 133, 241,  20, 224, 212,  36, 104, 184,
	 10,  49, 124, 255,  44, 110,   5,  50, 197, 222,  54, 109,  85, 146, 206, 158,
	126, 125, 105, 187, 185, 200,  30,  66, 152, 204,  48,  11,  28, 149,  17, 230,
	116, 233, 214,  98, 127, 103, 156, 142, 178, 205, 180, 195, 189,  59,  73,  13,
};
char permute_temp[256];
void shuffle(char *buf, int iw, int ih, int fwd)
{
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih-1;ky+=2)
		{
			int y=fwd?ky:((ih-1)&~1)-ky;
			y=(y+permute_idx[y&0xFF])%ih;
			int ycount=16;
			if(ycount>ih-y)
				ycount=ih-y;
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int x=fwd?kx:((iw-1)&~1)-kx;
				x=(x+permute_idx[0xFF-(x&0xFF)])%iw;
				int xcount=16;
				if(xcount>iw-x)
					xcount=iw-x;
				int blockcount=xcount*ycount;
				for(int ky2=0;ky2<ycount;++ky2)
				{
					for(int kx2=0;kx2<xcount;++kx2)
						permute_temp[xcount*ky2+kx2]=buf[(iw*(y+ky2)+x+kx2)<<2|kc];
				}
				int idx, idx2;
				char temp;
				if(fwd)
				{
					for(int ky2=0;ky2<ycount;++ky2)
					{
						for(int kx2=0;kx2<xcount;++kx2)
						{
							idx=xcount*ky2+kx2, idx2=permute_idx[idx]%blockcount;
							if(idx!=idx2)
								SWAPVAR(permute_temp[idx], permute_temp[idx2], temp);
						}
					}
				}
				else
				{
					for(int ky2=ycount-1;ky2>=0;--ky2)
					{
						for(int kx2=xcount-1;kx2>=0;--kx2)
						{
							idx=xcount*ky2+kx2, idx2=permute_idx[idx]%blockcount;
							if(idx!=idx2)
								SWAPVAR(permute_temp[idx], permute_temp[idx2], temp);
						}
					}
				}
				for(int ky2=0;ky2<ycount;++ky2)
				{
					for(int kx2=0;kx2<xcount;++kx2)
						buf[(iw*(y+ky2)+x+kx2)<<2|kc]=permute_temp[xcount*ky2+kx2];
				}
			}
		}
	}
}

//learned predictor
#if 1
void mulmatvec_pd(double *dst, const double *matrix, const double *vector, int mw, int mh)
{
	for(int ky=0;ky<mh;++ky)
	{
		double sum=0;
		for(int kx=0;kx<mw;++kx)
			sum+=matrix[mw*ky+kx]*vector[kx];
		dst[ky]=sum;
	}
}
void mulvTmat_pd(double *dst, const double *vT, const double *matrix, int mw, int mh)
{
	for(int kx=0;kx<mw;++kx)
	{
		double sum=0;
		for(int ky=0;ky<mh;++ky)
			sum+=vT[ky]*matrix[mw*ky+kx];
		dst[kx]=sum;
	}
}
void mulvecs_pd(double *dst, const double *a, const double *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]*b[k];
}
void scalevec_pd(double *vec, int count, double factor)
{
	for(int k=0;k<count;++k)
		vec[k]*=factor;
}
void filternan_pd(double *vec, int count)
{
	for(int k=0;k<count;++k)
	{
		if(!isfinite(vec[k])||fabs(vec[k])>1e6)
			vec[k]=0;
	}
}
void addvecs_pd(double *dst, const double *a, const double *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b[k];
}
void addvec1_pd(double *dst, const double *a, const double b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b;
}
double vecsum_pd(const double *data, int count)
{
	double sum=0;
	for(int k=0;k<count;++k)
		sum+=data[k];
	return sum;
}
void negbuffer_pd(double *data, int count)
{
	for(int k=0;k<count;++k)
		data[k]=-data[k];
}
void subvecs_pd(double *dst, const double *pos, const double *neg, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=pos[k]-neg[k];
}
void leakyReLU_pd(double *dst, const double *src, int count)
{
	for(int k=0;k<count;++k)
	{
		if(src[k]<0)
			dst[k]=src[k]*0.01;
		else
			dst[k]=src[k];
	}
}
void leakyReLUdash_pd(double *data, int count)
{
	for(int k=0;k<count;++k)
	{
		if(data[k]<0)
			data[k]=0.01;
		else
			data[k]=1;
	}
}
double calc_rmse_pd(const double *error, int count)
{
	double sum=0;
	for(int k=0;k<count;++k)
		sum+=error[k]*error[k];
	sum/=count;
	sum=sqrt(sum);
	return sum;
}
void fill_matrow_unchecked(const unsigned char *buf, int iw, int ih, int kc, int kx, int ky, double *mrow, double *target)
{
	int idx=iw*ky+kx;
	mrow[0]=buf[(idx-iw*2-2)<<2|kc];
	mrow[1]=buf[(idx-iw*2-1)<<2|kc];
	mrow[2]=buf[(idx-iw*2)<<2|kc];
	mrow[3]=buf[(idx-iw*2+1)<<2|kc];
	mrow[4]=buf[(idx-iw*2+2)<<2|kc];

	mrow[5]=buf[(idx-iw-2)<<2|kc];
	mrow[6]=buf[(idx-iw-1)<<2|kc];
	mrow[7]=buf[(idx-iw)<<2|kc];
	mrow[8]=buf[(idx-iw+1)<<2|kc];
	mrow[9]=buf[(idx-iw+2)<<2|kc];

	mrow[10]=buf[(idx-2)<<2|kc];
	mrow[11]=buf[(idx-1)<<2|kc];

	*target=buf[idx<<2|kc];
}
double leakyReLU1(double x)
{
	if(x<0)
		x*=0.01;
	return x;
}
#define OPT_N 12
#define OPT_B 12
static double g_mat[OPT_N*OPT_B], g_y[OPT_N], g_v[OPT_B], g_v2[OPT_B], g_grad[OPT_N+1];
double opt_causal_reach2(unsigned char *buf, int iw, int ih, int kc, double *x, double *bias, double lr, int test)
{
	addhalf(buf, iw, ih, 3, 4);
	int nupdates=0;
	double rmse=0;
	for(int ky=2, row=0;ky<ih-2;ky+=3)
	{
		for(int kx=2;kx<iw-4;kx+=5)
		{
			fill_matrow_unchecked(buf, iw, ih, kc, kx, ky, g_mat+row*12, g_y+row);
			row=(row+1)%12;
			if(!row)
			{
#if 0
				//forward		pred = leakyReLU(M x + b)		x & b are params, M & y are pixels
				mulmatvec_pd(g_v,  g_mat, x,  OPT_N, OPT_B);	//     M x
				addvec1_pd(g_v,  g_v, *bias,  OPT_N);			//v  = M x + b
				leakyReLU_pd(g_v2,  g_v,  OPT_N);				//     leakyReLU(M x + b)
				subvecs_pd(g_v2,  g_y, g_v2,  OPT_N);			//v2 = y - leakyReLU(M x + b)
				
				rmse+=calc_rmse_pd(g_v2, OPT_N);				//L  = sum i: (y[i] - leakyReLU(M x[i] + b))^2

				if(!test)
				{
					//backward
					leakyReLUdash_pd(g_v, OPT_N);				//v  = leakyReLU'(M x + b)
					negbuffer_pd(g_v, OPT_N);					//v  = -leakyReLU'(M x + b)

					mulvecs_pd(g_v2, g_v2, g_v, OPT_N);			//v2 = 2 * (y - leakyReLU(M x + b)) * -leakyReLU'(M x + b)
					mulvTmat_pd(g_grad, g_v2, g_mat, OPT_N, OPT_N);	//gx = v2T M
					g_grad[OPT_N]=vecsum_pd(g_v2, OPT_N);			//gb = sum(v2)
				
					scalevec_pd(g_grad, OPT_N+1, lr);

					subvecs_pd(x, x, g_grad, OPT_N);
					*bias-=g_grad[OPT_N];
				}
				++nupdates;
#endif
#if 1
				mulmatvec_pd(g_v, g_mat, x, OPT_N, OPT_N);
				subvecs_pd(g_v, g_y, g_v, OPT_N);
				rmse+=calc_rmse_pd(g_v, OPT_N);
				mulvTmat_pd(g_grad, g_v, g_mat, OPT_N, OPT_N);
				scalevec_pd(g_grad, OPT_N, lr);
				addvecs_pd(x, x, g_grad, OPT_N);
				filternan_pd(x, OPT_N);
				++nupdates;
#endif
			}
		}
	}
	addhalf(buf, iw, ih, 3, 4);
	rmse/=nupdates;
	return rmse;
}


#if 0
typedef struct JointContextStruct
{
	int key;
	int hist[256];
} JointContext;
Map j_map;
CmpRes cmp_jctx(const void *key, const void *candidate)
{
	const int *k=(const int*)key, *c=(const int*)candidate;
	return (*k>*c)-(*k<*c);
}
void pred_bitwise(char *buf, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *buf2=(char*)malloc((size_t)res<<2);
	//int *hist=(int*)malloc(0x1000000*sizeof(int));
	if(!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int pred=0xFF000000;
	memfill(buf2, &pred, (size_t)res<<2, 4);
	//int grads[256]={0};
	MAP_INIT(&j_map, JointContext, cmp_jctx, 0);
	//memset(hist, 0, 0x1000000*sizeof(int));
	//for(int k=0;k<res;++k)
	//	((int*)buf2)[k]=0xFF00000;
	//memcpy(buf2, buf, (size_t)res<<2);
	const char *pixels=fwd?buf:buf2, *errors=fwd?buf2:buf;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==(iw>>1)&&ky==(ih>>1))
			//if(kx==4&&ky==3)
			//	printf("");
			for(int kc=0;kc<3;++kc)
			{
#define LOAD(BUF, C, X, Y) (unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|(kc-C)%3]:0
				char
					W =LOAD(pixels, 0,  1, 0),
					NW=LOAD(pixels, 0,  1, 1),
					N =LOAD(pixels, 0,  0, 1),
					NE=LOAD(pixels, 0, -1, 1);
#undef LOAD
				int idx=(iw*ky+kx)<<2|kc;
				int grad=N+W-NW;
				grad=CLAMP(-128, grad, 127);

				//int found0=grads[grad+128];
				//if(found0)
				//	printf("");

				int found=0;
				int key=kc<<16|(N+128)<<8|(W+128);
				RBNodeHandle *hnode=map_insert(&j_map, &key, &found);
				JointContext *node=(JointContext*)hnode[0]->data;
				int pred;
				if(found)
				{
					int best=0, sum=0;
					for(int sym=0;sym<256;++sym)
					{
						if(node->hist[best]<node->hist[sym])
							best=sym;
						sum+=node->hist[sym];
					}
					if(sum)
						pred=best-128;
					else
						pred=grad;
				}
				else
				{
					node->key=key;
					pred=grad;
				}

				if(fwd)
					buf2[idx]=buf[idx]-pred;
				else
					buf2[idx]=buf[idx]+pred;
				int curr=pixels[idx]+128;
				++node->hist[curr];
				//grads[grad+128]=1;

#if 0
				//int vmin, vmax;
				//int emin, emax;
				//if(N<W)
				//	vmin=N, vmax=W;
				//else
				//	vmin=W, vmax=N;
				//grad=CLAMP(vmin, grad, vmax);
#if 1
				grad+=128;
				int pixel=0, error=0;
				int kb=7, bit0, pred, diff=0;
				for(;kb>=0;--kb)
				{
					int vmax=pixel+(1<<(kb+1))-1;
					grad=CLAMP(pixel, grad, vmax);
					pred=grad>=pixel+(1<<kb);
					if(fwd)
					{
						bit0=(pixels[idx]+128)>>kb&1;//original bit
						diff=bit0-pred;
					}
					else
					{
						diff=errors[idx]>>kb&1;//"zero" bit
						bit0=diff+pred;
					}
					pixel^=bit0<<kb;
					error^=(diff&1)<<kb;
					if(diff)
					{
						int mask, rem;

						mask=(1<<kb)-1;
						//pred=-(diff>0)&mask;//X
						pred=0;
						//if(mask)
						//	printf("");
						if(fwd)
						{
							rem=(pixels[idx]+128)&mask;
							error+=rem-pred;
						}
						else
						{
							rem=errors[idx]&mask;
							pixel+=rem+pred;
						}
#if 0
						--kb;
						pred=diff>0;
						int pixel2=0, error2=0, kb0=kb;
						for(;kb>=0;--kb)
						{
							if(fwd)
							{
								bit0=pixels[idx]>>kb&1;//original bit
								diff=bit0-pred;
							}
							else
							{
								diff=errors[idx]>>kb&1;//"zero" bit
								bit0=diff+pred;
							}
							pixel2+=bit0<<kb;
							error2+=(diff&1)<<kb;
						}
						pixel|=pixel2&(1<<(kb0+1))-1;
						error|=error2&(1<<(kb0+1))-1;
#endif
						break;
					}
				}
				buf2[idx]=fwd?error:pixel;
#endif
#if 0
				buf2[idx]=buf[idx];
				int diff=0;
				for(int kb=7, mask=0xFFFFFF00;kb>=0;--kb, mask>>=1)
				{
					vmin=pixels[idx]&mask;
					vmax=vmin+(1<<(kb+1))-1;
					//emin=errors[idx]&mask;
					//emax=emin+(1<<(kb+1))-1;
					grad=CLAMP(vmin, grad, vmax);
					//if(kb<7)
					//{
					//}
					int vmid=(vmin+vmax+1)>>1;
					pred=grad>=vmid;
					if(fwd)
						buf2[idx]-=pred<<kb;
					else
						buf2[idx]+=pred<<kb;
					int bit0=pixels[idx]>>kb&1;
					diff=bit0-pred;
					if(diff)
					{
						--kb;
						for(;kb>=0;--kb)
						{
							pred=diff>0;
							if(fwd)
								buf2[idx]-=pred<<kb;
							else
								buf2[idx]+=pred<<kb;
						}
					}
				}
#endif
#endif
			}
		}
	}
	set_window_title("%d nodes, %lf MB", j_map.nnodes, (double)j_map.nnodes*sizeof(JointContext)/(1024*1024));
	memcpy(buf, buf2, (size_t)res<<2);
	free(buf2);
	MAP_CLEAR(&j_map);
}
#endif


//	#define LOGIC_WEIGHTED_AVERAGE1
//	#define LOGIC_WEIGHTED_AVERAGE2

int logic_opt_timeron=0;
short logic_params[LOGIC_TOTALPARAMS]={0};
void pred_logic_prealloc(const char *src, int iw, int ih, int kc, int fwd, char *dst, const short *params)//params are fixed 3.12 bit
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int idx;
	long long pred=0;
#if 1
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//fixed 19.12 bit
			int nb[LOGIC_NNB*2]=
			{
#define LOAD(BUF, XO, YO) (unsigned)(kx+XO)<(unsigned)iw&&(unsigned)(ky+YO)<(unsigned)ih?src[(iw*(ky+YO)+kx+XO)<<2|kc]<<12:0
#if LOGIC_REACH>=2
				LOAD(pixels, -2, -2),
				LOAD(pixels, -1, -2),
				LOAD(pixels,  0, -2),
				LOAD(pixels,  1, -2),
				LOAD(pixels,  2, -2),
				LOAD(pixels, -2, -1),
#endif
				LOAD(pixels, -1, -1),
				LOAD(pixels,  0, -1),
				LOAD(pixels,  1, -1),
#if LOGIC_REACH>=2
				LOAD(pixels,  2, -1),
				LOAD(pixels, -2,  0),
#endif
				LOAD(pixels, -1,  0),
				
#if LOGIC_REACH>=2
				LOAD(errors, -2, -2),
				LOAD(errors, -1, -2),
				LOAD(errors,  0, -2),
				LOAD(errors,  1, -2),
				LOAD(errors,  2, -2),
				LOAD(errors, -2, -1),
#endif
				LOAD(errors, -1, -1),
				LOAD(errors,  0, -1),
				LOAD(errors,  1, -1),
#if LOGIC_REACH>=2
				LOAD(errors,  2, -1),
				LOAD(errors, -2,  0),
#endif
				LOAD(errors, -1,  0),
#undef  LOAD
			};
			int temp[LOGIC_NF0]={0};
			const short *lastrow=params+LOGIC_ROWPARAMS*LOGIC_NF0;
			for(int kc=0;kc<LOGIC_NF0;++kc)
			{
				const short *row=params+LOGIC_ROWPARAMS*kc;
				long long sum=0;

				//weighted average
				long long wsum=0;
				for(int kx=0;kx<LOGIC_NNB*2;++kx)
				{
					sum+=(long long)nb[kx]*row[kx];
					wsum+=row[kx];
				}
				if(wsum)
					temp[kc]=(int)(sum/wsum);
				else
					temp[kc]=(int)((sum+(1LL<<11))>>12);

				//kernel
				//for(int kx=0;kx<LOGIC_NNB*2;++kx)
				//	sum+=(long long)nb[kx]*row[kx];
				//temp[kc]=(int)((sum+(1LL<<11))>>12);
			}

			pred=0;
			long long wsum=0;
			for(int kc=0;kc<LOGIC_NF0;++kc)
			{
				pred+=(long long)temp[kc]*lastrow[kc];
				wsum+=lastrow[kc];
			}
			if(wsum)
				pred/=wsum;
			else
				pred>>=12;
			pred+=1LL<<11;//rounding
			pred>>=12;
			pred=CLAMP(-128, pred, 128);

			idx=(iw*ky+kx)<<2|kc;
			if(fwd)
				dst[idx]=(char)(src[idx]-pred);
			else
				dst[idx]=(char)(src[idx]+pred);
		}
	}
#endif
#if 0
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//fixed 19.12 bit
			int nb[LOGIC_NNB*2]=
			{
#define LOAD(BUF, XO, YO) (unsigned)(kx+XO)<(unsigned)iw&&(unsigned)(ky+YO)<(unsigned)ih?src[(iw*(ky+YO)+kx+XO)<<2|kc]<<12:0
#if LOGIC_REACH>=2
				LOAD(pixels, -2, -2),
				LOAD(pixels, -1, -2),
				LOAD(pixels,  0, -2),
				LOAD(pixels,  1, -2),
				LOAD(pixels,  2, -2),
				LOAD(pixels, -2, -1),
#endif
				LOAD(pixels, -1, -1),
				LOAD(pixels,  0, -1),
				LOAD(pixels,  1, -1),
#if LOGIC_REACH>=2
				LOAD(pixels,  2, -1),
				LOAD(pixels, -2,  0),
#endif
				LOAD(pixels, -1,  0),
				
#if LOGIC_REACH>=2
				LOAD(errors, -2, -2),
				LOAD(errors, -1, -2),
				LOAD(errors,  0, -2),
				LOAD(errors,  1, -2),
				LOAD(errors,  2, -2),
				LOAD(errors, -2, -1),
#endif
				LOAD(errors, -1, -1),
				LOAD(errors,  0, -1),
				LOAD(errors,  1, -1),
#if LOGIC_REACH>=2
				LOAD(errors,  2, -1),
				LOAD(errors, -2,  0),
#endif
				LOAD(errors, -1,  0),
#undef  LOAD
			};
			int temp[LOGIC_NF0]={0};
			const short *lastrow=params+LOGIC_ROWPARAMS*LOGIC_NF0;
#ifdef LOGIC_WEIGHTED_AVERAGE1
			for(int kc=0;kc<LOGIC_NF0;++kc)
			{
				const short *row=params+LOGIC_ROWPARAMS*kc;
				long long
					sum=row[LOGIC_NNB*2],//bias
					wsum=0;
				for(int kx=0;kx<LOGIC_NNB*2;++kx)
				{
					sum+=(long long)nb[kx]*row[kx];
					wsum+=row[kx];
				}
				if(wsum)
					temp[kc]=(int)(sum/wsum);
				else
					temp[kc]=(int)((sum+(1LL<<11))>>12);
			}
#else
			for(int kc=0;kc<LOGIC_NF0;++kc)
			{
				const short *row=params+LOGIC_ROWPARAMS*kc;
				long long sum=row[LOGIC_NNB*2];//bias
				for(int kx=0;kx<LOGIC_NNB*2;++kx)
					sum+=(long long)nb[kx]*row[kx];
				temp[kc]=(int)((sum+(1LL<<11))>>12);
			}
#endif
			for(int kc=0;kc<LOGIC_NF1;++kc)
			{
				int cond=temp[kc+LOGIC_NF1*2LL];
				cond=CLAMP(-0x1000, cond, 0x1000);
				temp[kc]=temp[kc]+(int)((long long)(temp[kc+LOGIC_NF1]-temp[kc])*cond>>12);
			}

			pred=0;
#ifdef LOGIC_WEIGHTED_AVERAGE2
			long long wsum=0;
			for(int kc=0;kc<LOGIC_NF1;++kc)
			{
				pred+=(long long)temp[kc]*lastrow[kc];
				wsum+=lastrow[kc];
			}
			if(wsum)
				pred/=wsum;
			else
				pred>>=12;
			pred+=1LL<<11;//rounding
			pred>>=12;
#else
			for(int kc=0;kc<LOGIC_NF0;++kc)
				pred+=(long long)temp[kc]*lastrow[kc];
			pred+=1LL<<27;//rounding
			pred>>=28;
#endif
			pred=CLAMP(-128, pred, 128);

			idx=(iw*ky+kx)<<2|kc;
			if(fwd)
				dst[idx]=(char)(src[idx]-pred);
			else
				dst[idx]=(char)(src[idx]+pred);
		}
	}
#endif
}
void pred_logic_apply(char *buf, int iw, int ih, const short *allparams, int fwd)
{
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(temp, buf, (size_t)res<<2);//copy alpha
	pred_logic_prealloc(buf, iw, ih, 0, fwd, temp, allparams);
	pred_logic_prealloc(buf, iw, ih, 1, fwd, temp, allparams+LOGIC_PARAMS_PER_CH);
	pred_logic_prealloc(buf, iw, ih, 2, fwd, temp, allparams+LOGIC_PARAMS_PER_CH*2);
	memcpy(buf, temp, (size_t)res<<2);
	free(temp);
}
typedef struct LogicInfoStruct
{
	double loss;
	short params[LOGIC_PARAMS_PER_CH];
} LogicInfo;
int cmp_logicinfo(const void *left, const void *right)
{
	LogicInfo const *a, *b;

	a=(LogicInfo const*)left;
	b=(LogicInfo const*)right;
	return (a->loss>b->loss)-(a->loss<b->loss);//ascending order
}
typedef struct LogicThreadInfoStruct
{
	char *buf, *temp;//not const for free()
	int hist[256];
	LogicInfo *params;
	int iw, ih, kc;
	short *channel_params;
	float progressbar;//goes from 0 to 100
	float t_elapsed;//seconds
	float loss;
	int bitidx;
	int toggle_histogram[LOGIC_PARAMS_PER_CH<<4];
} LogicThreadInfo;
LogicThreadInfo logic_info={0};
HANDLE logic_hthread=0, ghMutex=0;
void logic_thread_update(LogicThreadInfo *info, const short *params, float progressbar, float t_start, float loss, const char *buf, int bitidx)
{
	int waitstatus=WaitForSingleObject(ghMutex, INFINITE);
	switch(waitstatus)
	{
	case WAIT_OBJECT_0:
		if(params)
			memcpy(info->channel_params, params, LOGIC_PARAMS_PER_CH*sizeof(short));
		if(buf)
		{
			int res=info->iw*info->ih;
			for(int k=0;k<res;++k)
			{
				unsigned char sym=buf[k<<2|info->kc]+128;
				image[k<<2  ]=sym;
				image[k<<2|1]=sym;
				image[k<<2|2]=sym;
			}
		}
		//if(buf)
		//	memcpy(image, buf, (size_t)info->iw*info->ih<<2);
		info->progressbar=progressbar;
		info->t_elapsed=(float)((time_ms()-t_start)*0.001);
		info->loss=loss;
		info->bitidx=bitidx;
		++info->toggle_histogram[bitidx];
		ReleaseMutex(ghMutex);
		break;
	case WAIT_ABANDONED:
		break;
	}
}
double logic_calcloss(short *params, LogicThreadInfo *info, char *temp, int *hist, float progressbar, float t_start, int bitidx)
{
	int res=info->iw*info->ih;
	pred_logic_prealloc(info->buf, info->iw, info->ih, info->kc, 1, temp, params);

	memset(hist, 0, 256*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char sym=temp[k<<2|info->kc];
		++hist[sym];
	}
	double entropy=0;
	for(int sym=0;sym<256;++sym)
	{
		int freq=hist[sym];
		if(freq)
		{
			double prob=(double)freq/res;
			double bitsize=-log2(prob);
			if(isinf(bitsize))
				LOG_ERROR("loss error");
			entropy+=prob*bitsize;
		}
	}
	double invCR=entropy/8;
	logic_thread_update(info, params, progressbar, t_start, (float)invCR, temp, bitidx);
	return invCR;
}
#if 0
unsigned __stdcall logic_opt_thread(LogicThreadInfo *info)//Nelder-Mead optimizer
{
	double t_start=time_ms();
	int res=info->iw*info->ih;
	const int nv=LOGIC_PARAMS_PER_CH, np=LOGIC_PARAMS_PER_CH+1;
	//char *temp=(char*)malloc((size_t)res<<2);
	//int *hist=malloc(256*sizeof(int));
	//LogicInfo *params=(LogicInfo*)malloc((np+3LL)*sizeof(LogicInfo));
	//if(!temp||!hist||!params)
	//{
	//	LOG_ERROR("Allocation error");
	//	return 1;
	//}
	LogicInfo
		*best=info->params,
		*worst=info->params+np-1,
		*x0=info->params+np,
		*xr=x0+1,
		*x2=xr+1;
	
#define CALC_LOSS(X, P) (X)->loss=logic_calcloss((X)->params, info, info->temp, info->hist, P, (float)t_start)
	
	//initialize N+1 param sets
	srand((unsigned)__rdtsc());
	for(int kp=0;kp<np;++kp)
	{
		LogicInfo *x=best+kp;
		//info->progressbar=(float)(kp+1)/np;
		//set_window_title("C%d it 0/100: %lf, elapsed %lf, %d/%d", kc, params->loss, (time_ms()-t_start)*0.001, kp, np);
		if(kp)
		{
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=info->channel_params[k2]+(rand()&0x1FFF)-0x1000;
		}
		else
			memcpy(x->params, info->channel_params, sizeof(x->params));
		CALC_LOSS(x, (float)(kp+1)/np);
		//logic_thread_update(info, kp?0:params->params, (float)params->loss, (float)(kp+1)/np, temp);
	}
	double loss0=best->loss;
	const int alpha=0x10000, gamma=0x20000, rho=0x8000, sigma=0x8000;
	for(int ki=0;ki<100;++ki)
	{
		//1  order
		isort(best, np, sizeof(LogicInfo), cmp_logicinfo);
		
		//logic_thread_update(info, params->params, (float)params->loss, (float)(ki+1), temp);
		//set_window_title("C%d it %d/100: %lf, elapsed %lf", kc, ki+1, params->loss, (time_ms()-t_start)*0.001);
		//memcpy(channel_params, params->params, sizeof(params->params));
		//io_render();

		//2  get the centroid of all points except worst
		memset(x0->params, 0, sizeof(x0->params));
		for(int k2=0;k2<np-1;++k2)//exclude the worst point
		{
			LogicInfo *x=best+k2;
			for(int k3=0;k3<nv;++k3)
				x0->params[k3]+=x->params[k3];
		}
		for(int k2=0;k2<nv;++k2)
			x0->params[k2]/=nv;

		//3  reflection
		for(int k2=0;k2<nv;++k2)
			xr->params[k2]=x0->params[k2]+(int)((long long)(x0->params[k2]-worst->params[k2])*alpha>>16);
		CALC_LOSS(xr, (float)(ki+1));
		if(xr->loss>best->loss&&xr->loss<worst[-1].loss)//if xr is between best and 2nd worst, replace worst with xr
		{
			memcpy(worst, xr, sizeof(LogicInfo));
			continue;
		}

		//4  expansion
		if(xr->loss<best->loss)//if xr is best so far
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(int)((long long)(xr->params[k2]-x0->params[k2])*gamma>>16);
			CALC_LOSS(x2, (float)(ki+1));
			if(x2->loss<xr->loss)
				memcpy(worst, x2, sizeof(LogicInfo));
			else
				memcpy(worst, xr, sizeof(LogicInfo));
			continue;
		}

		//5  contraction
		if(xr->loss<worst->loss)//if xr is between 2nd worst and worst
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(int)((long long)(xr->params[k2]-x0->params[k2])*rho>>16);
			CALC_LOSS(x2, (float)(ki+1));
			if(x2->loss<xr->loss)//if contracted point is better than xr
			{
				memcpy(worst, x2, sizeof(LogicInfo));
				continue;
			}
		}
		else
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(int)((long long)(worst->params[k2]-x0->params[k2])*rho>>16);
			CALC_LOSS(x2, (float)(ki+1));
			if(x2->loss<worst->loss)//if contracted point is better than xr
			{
				memcpy(worst, x2, sizeof(LogicInfo));
				continue;
			}
		}

		//6  shrink
		for(int kp=1;kp<np;++kp)
		{
			//logic_thread_update(info, 0, (float)params->loss, ki+1+(float)(kp+1)/np, temp);
			//info->progressbar=ki+1+(float)(kp+1)/np;
			//set_window_title("C%d it %d/100: %lf, elapsed %lf, %d/%d", kc, ki+1, params->loss, (time_ms()-t_start)*0.001, kp, np);
			LogicInfo *x=best+kp;
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=best->params[k2]+(int)((long long)(x->params[k2]-best->params[k2])*sigma>>16);
			CALC_LOSS(x, ki+1+(float)(kp+1)/np);
		}
	}
	CALC_LOSS(best, 100);
	//info->progressbar=101;
#undef CALC_LOSS
	//if(params->loss>loss0)
	//{
	//	messagebox(MBOX_OK, "Error", "Loss has increased");
	//}
	//logic_thread_update(info, 0, (float)params->loss, 100, temp);
	//memcpy(info->channel_params, params->params, sizeof(params->params));

	//free(params);
	//free(temp);
	//free(hist);
	info->progressbar=-1;
	return 0;
}
#endif
unsigned __stdcall logic_opt_thread_v2(LogicThreadInfo *info)
{
	double t_start=time_ms();
	const int nv=LOGIC_PARAMS_PER_CH;
	int bitidx=0;
	double bestloss;
	LogicInfo *set1=info->params, *set2=info->params+1;
	memcpy(set1->params, info->channel_params, sizeof(set1->params));
	memcpy(set2->params, info->channel_params, sizeof(set2->params));
	
#define GET_RAND_IDX() (xoroshiro128_next()%(LOGIC_PARAMS_PER_CH<<4))
#define TOGGLE_BIT(X, IDX) (X)->params[IDX>>4]^=1<<(IDX&15)
#define CALC_LOSS(X, P) (X)->loss=logic_calcloss((X)->params, info, info->temp, info->hist, P, (float)t_start, bitidx)
	
	srand((unsigned)__rdtsc());
	xoroshiro128_state[0]^=rand();
	xoroshiro128_state[1]^=rand();

	const int niter=4800;
	//const int niter=250;
	CALC_LOSS(set1, 0);
	bestloss=set1->loss;
	for(int ki=0, set2divergence=0;ki<niter;++ki, ++set2divergence)
	{
		bitidx=GET_RAND_IDX();
		//bitidx=(bitidx+rand())%(LOGIC_PARAMS_PER_CH<<4);
		TOGGLE_BIT(set1, bitidx);
		//set1->params[bitidx>>4]^=1<<(bitidx&15);
		CALC_LOSS(set1, ki*100.f/niter);
		if(set1->loss>bestloss)//revert set1 if worse
		{
			set1->loss=bestloss;
			TOGGLE_BIT(set1, bitidx);
			//set1->params[bitidx>>4]^=1<<(bitidx&15);
		}
		else if(set1->loss<bestloss)
		{
			bestloss=set1->loss;
			memcpy(set2, set1, sizeof(*set2));//set1 won, overwrite set2
		}
		
		bitidx=GET_RAND_IDX();
		if(set2divergence>25)//set2 went too far, overwrite set2
		{
			set2divergence=0;
			memcpy(set2, set1, sizeof(*set2));
		}
		TOGGLE_BIT(set2, bitidx);
		CALC_LOSS(set2, ki*100.f/niter);//don't revert set2 if worse
		if(set2->loss<bestloss)
		{
			bestloss=set2->loss;
			memcpy(set1, set2, sizeof(*set1));//set2 won, overwrite set1
			set2divergence=0;
		}
	}
	CALC_LOSS(set1, 100);
#undef GET_RAND_IDX
#undef TOGGLE_BIT
#undef CALC_LOSS
	info->progressbar=-1;
	return 0;
}
static void logic_opt_freeresources()
{
	CloseHandle(logic_hthread);
	logic_hthread=0;
	CloseHandle(ghMutex);
	ghMutex=0;
	free(logic_info.buf);
	free(logic_info.temp);
	free(logic_info.params);
	memset(&logic_info, 0, sizeof(logic_info));
	timer_stop(TIMER_ID_MONITOR);
}
void logic_opt_checkonthread(float *info)
{
	float progressbar=-1, t_elapsed=-1, loss=-1;
	int bitidx=-1;
	if(logic_hthread)
	{
		int waitstatus=WaitForSingleObject(ghMutex, INFINITE);
		switch(waitstatus)
		{
		case WAIT_OBJECT_0:
			progressbar=logic_info.progressbar;
			t_elapsed=logic_info.t_elapsed;
			loss=logic_info.loss;
			bitidx=logic_info.bitidx;
			if(progressbar<0)
			{
				logic_opt_freeresources();
				//CloseHandle(logic_hthread);
				//logic_hthread=0;
				//CloseHandle(ghMutex);
				//ghMutex=0;
				//free(logic_info.buf);
				//free(logic_info.temp);
				//free(logic_info.params);
				//memset(&logic_info, 0, sizeof(logic_info));
				update_image();
			}
			else
				ReleaseMutex(ghMutex);
			break;
		case WAIT_ABANDONED:
			break;
		}
	}
	if(info)
	{
		info[0]=progressbar;
		info[1]=t_elapsed;
		info[2]=loss;
		((int*)info)[3]=bitidx;
	}
}
void logic_opt_forceclosethread()
{
	if(logic_hthread)
	{
		TerminateThread(logic_hthread, 0);
		logic_opt_freeresources();
		//CloseHandle(logic_hthread);
		//logic_hthread=0;
		//CloseHandle(ghMutex);
		//ghMutex=0;
	}
}
void logic_opt(char *buf, int iw, int ih, int kc, short *channel_params)
{
	int res=iw*ih;
    ghMutex=CreateMutexA(0, 0, 0);//default security parameters, initially not owned, unnamed
	if(!ghMutex)
	{
		LOG_ERROR("Mutex allocation error");
		return;
	}
	logic_info.buf=buf;
	logic_info.temp=(char*)malloc((size_t)res<<2);
	logic_info.params=(LogicInfo*)malloc(2*sizeof(LogicInfo));
	//logic_info.params=(LogicInfo*)malloc((LOGIC_PARAMS_PER_CH+1+3LL)*sizeof(LogicInfo));
	if(!logic_info.temp||!logic_info.params)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	logic_info.iw=iw;
	logic_info.ih=ih;
	logic_info.kc=kc;
	logic_info.channel_params=channel_params;
	logic_info.progressbar=0;
	memset(logic_info.toggle_histogram, 0, sizeof(logic_info.toggle_histogram));
	logic_hthread=(HANDLE)_beginthreadex(0, 0, logic_opt_thread_v2, &logic_info, 0, 0);
	//logic_hthread=(HANDLE)_beginthreadex(0, 0, logic_opt_thread, &logic_info, 0, 0);
	if(!logic_hthread)
		LOG_ERROR("Thread error");
	timer_start(TIMER_MONITOR_MS, TIMER_ID_MONITOR);
}

#define O2_NPARAMS (_countof(customparam_st)/3)
typedef struct OptCR2InfoStruct
{
	double loss, params[O2_NPARAMS];
} OptCR2Info;
int cmp_optcr2info(const void *left, const void *right)
{
	OptCR2Info const *a, *b;

	a=(OptCR2Info const*)left;
	b=(OptCR2Info const*)right;
	return (a->loss>b->loss)-(a->loss<b->loss);//ascending order
}
double opt_cr2_calcloss(double *params, const char *buf, int iw, int ih, int kc, char *temp, int *hist)
{
	int res=iw*ih;
	//memcpy(customparam_st+O2_NPARAMS*customparam_st_ch, params, sizeof(customparam_st)/6);
	//memcpy(temp, buf, (size_t)res<<2);
	pred_custom_prealloc(buf, iw, ih, kc, 1, params, temp);
	//pred_custom_fwd(temp+kc, iw, ih, 1, 4, params);

	memset(hist, 0, 256*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char sym=temp[k<<2|kc];
		++hist[sym];
	}
	double entropy=0;
	for(int sym=0;sym<256;++sym)//Shannon law
	{
		int freq=hist[sym];
		if(freq)
		{
			double prob=(double)freq/res;
			entropy-=prob*log2(prob);
		}
	}
	double invCR=entropy/8;
#if 0
	double csize=0;
	for(int k=0;k<res;++k)//Zipf's law
	{
		unsigned char sym=temp[k<<2|kc];
		int freq=hist[sym];
		if(freq)
		{
			double prob=(double)freq/res;
			csize-=log2(prob);
		}
	}
	csize/=8;
	double csize2=res*invCR;
	if(fabs(csize-csize2)>1e-3)
		LOG_ERROR("BROKEN");
#endif
	//for(int k=0;k<res;++k)
	//{
	//	image[k<<2  ]=temp[k<<2|kc]+128;
	//	image[k<<2|1]=temp[k<<2|kc]+128;
	//	image[k<<2|2]=temp[k<<2|kc]+128;
	//}
	//memcpy(image, temp, (size_t)res<<2);
	//io_render();

	return invCR;
}
void opt_cr2_v2(const char *buf, int iw, int ih, int kc)
{
	//double params[_countof(customparam_st)];
	//memcpy(params, customparam_st, sizeof(params));//save params
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	int *hist=malloc(256*sizeof(int));
	const int nv=O2_NPARAMS, np=O2_NPARAMS+1;
	OptCR2Info *params=(OptCR2Info*)malloc((np+3LL)*sizeof(OptCR2Info));
	if(!temp||!hist||!params)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, (size_t)res<<2);
	OptCR2Info
		*best=params,
		*worst=params+np-1,
		*x0=params+np,
		*xr=x0+1,
		*x2=xr+1;
	
	double initial_step=0.1;
#define CALC_LOSS(X) (X)->loss=opt_cr2_calcloss((X)->params, buf, iw, ih, kc, temp, hist)
	
	//initialize N+1 param sets
	srand((unsigned)__rdtsc());
	for(int kp=0;kp<np;++kp)
	{
		OptCR2Info *x=best+kp;
		//if(kp)
		//{
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=customparam_st[O2_NPARAMS*kc+k2]+((double)rand()-((RAND_MAX>>1)+1))*(initial_step/RAND_MAX);
		//}
		//else
		//	memcpy(x->params, customparam_st+O2_NPARAMS*kc, sizeof(x->params));
		CALC_LOSS(x);
		//x->loss=opt_cr2_calcloss(x->params, buf, iw, ih, kc, temp, hist);
	}
	memcpy(x0->params, customparam_st+O2_NPARAMS*kc, sizeof(x0->params));
	CALC_LOSS(x0);
	double loss0=x0->loss;
	//extern float ch_cr[4];
	//if(fabs(1/loss0-ch_cr[kc])>1e-3)
	//	LOG_ERROR("CR: loss %lf ch_cr %lf", 1/loss0, ch_cr[kc]);
	const double alpha=1, gamma=2, rho=0.5, sigma=0.5;
	for(int ki=0;ki<100;++ki)
	{
		//1  order
		isort(best, np, sizeof(OptCR2Info), cmp_optcr2info);
		
		set_window_title("it %d/100: %lf", ki+1, best->loss);

		//2  get the centroid of all points except worst
		memset(x0->params, 0, sizeof(x0->params));
		for(int k2=0;k2<nv;++k2)//exclude the worst point
		{
			OptCR2Info *x=best+k2;
			for(int k3=0;k3<nv;++k3)
				x0->params[k3]+=x->params[k3];
		}
		for(int k2=0;k2<nv;++k2)
			x0->params[k2]/=nv;

		//3  reflection
		for(int k2=0;k2<nv;++k2)
			xr->params[k2]=x0->params[k2]+(x0->params[k2]-worst->params[k2])*alpha;
		CALC_LOSS(xr);
		if(xr->loss>params->loss&&xr->loss<worst[-1].loss)//if xr is between best and 2nd worst, replace worst with xr
		{
			memcpy(worst, xr, sizeof(OptCR2Info));
			continue;
		}

		//4  expansion
		if(xr->loss<best->loss)//if xr is best so far
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(xr->params[k2]-x0->params[k2])*gamma;
			CALC_LOSS(x2);
			if(x2->loss<xr->loss)
				memcpy(worst, x2, sizeof(OptCR2Info));
			else
				memcpy(worst, xr, sizeof(OptCR2Info));
			continue;
		}

		//5  contraction
		if(xr->loss<worst->loss)//if xr is between 2nd worst and worst
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(xr->params[k2]-x0->params[k2])*rho;
			CALC_LOSS(x2);
			if(x2->loss<xr->loss)//if contracted point is better than xr
			{
				memcpy(worst, x2, sizeof(OptCR2Info));
				continue;
			}
		}
		else
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(worst->params[k2]-x0->params[k2])*rho;
			CALC_LOSS(x2);
			if(x2->loss<worst->loss)//if contracted point is better than xr
			{
				memcpy(worst, x2, sizeof(OptCR2Info));
				continue;
			}
		}

		//6  shrink
		for(int kp=1;kp<np;++kp)
		{
			OptCR2Info *x=best+kp;
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=best->params[k2]+(x->params[k2]-best->params[k2])*sigma;
			CALC_LOSS(x2);
		}
	}
#undef CALC_LOSS
	//if(best->loss>loss0)
	//{
	//	messagebox(MBOX_OK, "Error", "Loss has increased");
	//}
	if(best->loss<loss0)
		memcpy(customparam_st+O2_NPARAMS*kc, best->params, sizeof(customparam_st)/3);

#if 0
	double l1, losses[_countof(customparam_st)<<1];
	double steps[]={2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625};
	
	srand((unsigned)__rdtsc());
	l1=opt_cr2_calcloss(buf, iw, ih, kc, temp, hist);
	for(int ks=0;ks<_countof(steps);++ks)
	{
		double step=steps[rand()%_countof(steps)];
		for(int kp=0;kp<_countof(losses);++kp)
		{
			customparam_st[kp>>1]+=kp&1?-step:step;
			losses[kp]=opt_cr2_calcloss(buf, iw, ih, kc, temp, hist);
			customparam_st[kp>>1]-=kp&1?-step:step;
		}
	}
#endif
	free(params);
	free(temp);
	free(hist);
}

double opt_cr2_calcloss_v3(double *params, const char *buf, int iw, int ih, int kc, char *temp, int *hist, float t_start, float progress)
{
	int res=iw*ih;
	pred_custom_prealloc(buf, iw, ih, kc, 1, params, temp);

	memset(hist, 0, 256*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char sym=temp[k<<2|kc];
		++hist[sym];
	}
	double entropy=0;
	for(int sym=0;sym<256;++sym)
	{
		int freq=hist[sym];
		if(freq)
		{
			double prob=(double)freq/res;
			entropy-=prob*log2(prob);
		}
	}
	double invCR=entropy/8;
	//x->loss=invCR;

	TimeInfo ti;
	parsetimedelta(time_ms()-t_start, &ti);
	set_window_title("%6.2f%%  %02d-%02d-%06.3f  CR %lf", progress, ti.hours, ti.mins, ti.secs, 1/invCR);

	return invCR;
}
typedef struct OptCR2Info3Struct
{
	double loss;
	short params[O2_NPARAMS];
} OptCR2Info3;
void set_params_pd(double *params_pd, const short *params_i16, int count)
{
	for(int k=0;k<count;++k)
		params_pd[k]=params_i16[k]/4096.;
}
void get_params_pd(short *params_i16, const double *params_pd, int count)
{
	for(int k=0;k<count;++k)
		params_i16[k]=(int)(params_pd[k]*4096);
}
void opt_cr2_v3(const char *buf, int iw, int ih, int kc)//X  inefficient
{
	double t_start=time_ms();
	double *params0=customparam_st+O2_NPARAMS*kc;
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	int *hist=malloc(256*sizeof(int));
	const int nv=O2_NPARAMS;
	OptCR2Info3 sets[2];
	double params_pd[O2_NPARAMS];
	int bitidx=0;
	double bestloss;
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	get_params_pd(sets->params, params0, O2_NPARAMS);//copy params
	memcpy(sets[1].params, sets->params, sizeof(sets->params));
	
#define GET_RAND_IDX() (xoroshiro128_next()%(O2_NPARAMS<<4))
#define TOGGLE_BIT(X, IDX) (X)->params[IDX>>4]^=1<<(IDX&15)
#define CALC_LOSS(X, P) set_params_pd(params_pd, (X)->params, O2_NPARAMS), (X)->loss=opt_cr2_calcloss_v3(params_pd, buf, iw, ih, kc, temp, hist, (float)t_start, P)
	
	const int niter=2*(O2_NPARAMS<<4);
	srand((unsigned)__rdtsc());
	xoroshiro128_state[0]^=rand();
	xoroshiro128_state[1]^=rand();
	CALC_LOSS(sets, 0);
	bestloss=sets[0].loss;
	for(int ki=0, set2divergence=0;ki<niter;++ki, ++set2divergence)
	{
		bitidx=GET_RAND_IDX();
		TOGGLE_BIT(sets, bitidx);
		CALC_LOSS(sets, ki*100.f/niter);
		if(sets[0].loss>bestloss)//revert set1 if worse
		{
			sets[0].loss=bestloss;
			TOGGLE_BIT(sets, bitidx);
		}
		else if(sets[0].loss<bestloss)
		{
			bestloss=sets[0].loss;
			memcpy(sets+1, sets, sizeof(*sets));//set1 won, overwrite set2
		}
		
		bitidx=GET_RAND_IDX();
		if(set2divergence>25)//set2 went too far, overwrite set2
		{
			set2divergence=0;
			memcpy(sets+1, sets, sizeof(*sets));
		}
		TOGGLE_BIT(sets+1, bitidx);
		CALC_LOSS(sets+1, ki*100.f/niter);//don't revert set2 if worse
		if(sets[1].loss<bestloss)
		{
			bestloss=sets[1].loss;
			memcpy(sets, sets+1, sizeof(*sets));//set2 won, overwrite set1
			set2divergence=0;
		}
	}
	CALC_LOSS(sets, 100);
#undef GET_RAND_IDX
#undef TOGGLE_BIT
#undef CALC_LOSS
	if(sets->loss>bestloss)
	{
		messagebox(MBOX_OK, "Error", "Loss has increased");
	}
	set_params_pd(params0, sets->params, O2_NPARAMS);

	free(temp);
	free(hist);
}

#define O2_N 72
double customparam_hybrid[(7*3+3)*3*3]={0};//shape [3, 72]		maps 24 causal neighbors (72 values) -> 1 pixel (3 values)
double opt_causal_hybrid_r3(unsigned char *buf, int iw, int ih, double lr)
{
	double *x=customparam_hybrid;
	ArrayHandle tensors;
	ARRAY_ALLOC(double, tensors, 0, 0, 0, 0);
	int off_src=(int)ARRAY_APPEND_OFFSET(tensors, 0, O2_N, 1, 0);
	int off_grad=(int)ARRAY_APPEND_OFFSET(tensors, 0, O2_N*3, 1, 0);

	double *t=(double*)tensors->data;
	
	addhalf(buf, iw, ih, 3, 4);
	int nupdates=0;
	double rmse=0;
	for(int ky=3, row=0;ky<ih-3;ky+=4)
	{
		for(int kx=3;kx<iw-6;kx+=7)
		{
			int idx=iw*ky+kx;
			for(int ky2=0, kd=0;ky2<4;++ky2)//fetch causal neighbors
			{
				for(int kx2=0;kx2<7&&kd<O2_N;++kx2)
				{
					for(int kc=0;kc<3;++kc, ++kd)
						t[off_src+kd]=(char)buf[(idx+iw*(ky2-3)+kx2-3)<<2|kc];
				}
			}
			double v[3]={(char)buf[idx<<2], (char)buf[idx<<2|1], (char)buf[idx<<2|2]};

			for(int kc=0;kc<3;++kc)//v[3] = y[3] - x[3, 72] src.flatten[72]
			{
				for(int k=0;k<O2_N;++k)
					v[kc]-=x[kc*O2_N+k]*t[off_src+k];
			}

			rmse+=calc_rmse_pd(v, 3);

			for(int kc=0;kc<3;++kc)//grad = v[3] * -src.flatten[72]T		outer product
			{
				for(int k=0;k<O2_N;++k)
					t[off_grad+kc*O2_N+k]=v[kc]*-t[off_src+k];
			}
			
			scalevec_pd(t+off_grad, O2_N*3, lr);
			subvecs_pd(x, x, t+off_grad, O2_N*3);
			filternan_pd(x, O2_N*3);

			++nupdates;
		}
	}
	array_free(&tensors);
	rmse/=nupdates;
	addhalf(buf, iw, ih, 3, 4);
	return rmse;
}
void predict_hybrid(const char *buf, int iw, int kx, int ky, int idx, int rowlen, int *ret_pred)
{
	char cn[72];
	for(int ky2=0, kd=0;ky2<4;++ky2)//fetch causal neighbors
	{
		for(int kx2=0;kx2<7&&kd<O2_N;++kx2)
		{
			for(int kc=0;kc<3;++kc, ++kd)
				cn[kd]=ky2-3>=0&&BETWEEN_EXC(3, kx2, iw+3)>=0?buf[(idx+((iw*(ky2-3)+kx2-3)<<2))|kc]:0;
		}
	}
	for(int kc=0;kc<3;++kc)//v[3] = y[3] - x[3, 72] src.flatten[72]
	{
		double sum=0;
		for(int k=0;k<O2_N;++k)
			sum+=customparam_hybrid[kc*O2_N+k]*cn[k];
		sum=CLAMP(customparam_clamp[0], sum, customparam_clamp[1]);
		ret_pred[kc]=(int)round(sum);
	}
}
void pred_hybrid_fwd(char *buf, int iw, int ih)
{
	int rowlen=iw<<2;
	int idx=(iw*ih-1)<<2;
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx, idx-=4)
		{
			int pred[3];
			predict_hybrid(buf, iw, kx, ky, idx, rowlen, pred);
			
			buf[idx  ]-=pred[0];
			buf[idx|1]-=pred[1];
			buf[idx|2]-=pred[2];
		}
	}
}
void pred_hybrid_inv(char *buf, int iw, int ih)
{
	int rowlen=iw<<2;
	int idx=0;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred[3];
			predict_hybrid(buf, iw, kx, ky, idx, rowlen, pred);

			buf[idx  ]+=pred[0];
			buf[idx|1]+=pred[1];
			buf[idx|2]+=pred[2];
		}
	}
}
#endif

#define PREDICT_LINEAR(A, B)		(((B)<<1)-(A))
#define PREDICT_QUADRATIC(A, B, C)	((A)-3*(B)+3*(C))
#if 0
static int predict_simple(char topleft, char top, char left)
{
	int xdelta=top-topleft, ydelta=left-topleft, pred;
	if((xdelta>0)==(ydelta>0))
		pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
	else
		pred=topleft+xdelta+ydelta;//average slope
	return pred;
}
#endif
//int testhist[3]={0};//
int  predict_diff2d(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char top[]=
	{
		ky-4>=0?buf[idx-rowlen*4]:0,
		ky-3>=0?buf[idx-rowlen*3]:0,
		ky-2>=0?buf[idx-(rowlen<<1)]:0,
		ky-1>=0?buf[idx-rowlen]:0,
	};
	char left[]=
	{
		kx-4>=0?buf[idx- bytestride*4  ]:0,
		kx-3>=0?buf[idx- bytestride*3  ]:0,
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	char topleft[]=
	{
		kx-4>=0&&ky-4>=0?buf[idx- rowlen*4  - bytestride*4  ]:0,
		kx-3>=0&&ky-3>=0?buf[idx- rowlen*3  - bytestride*3  ]:0,
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx- rowlen    - bytestride    ]:0,
	};
	char topright[]=
	{
		kx+4<iw&&ky-4>=0?buf[idx- rowlen*4  + bytestride*4  ]:0,
		kx+3<iw&&ky-3>=0?buf[idx- rowlen*3  + bytestride*3  ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
		kx+1<iw&&ky-1>=0?buf[idx- rowlen    + bytestride    ]:0,
	};
	int errors[]=
	{
		abs(left[3]-PREDICT_LINEAR(left[1], left[2])),
		abs(top[3]-PREDICT_LINEAR(top[1], top[2])),
		abs(topleft[3]-PREDICT_LINEAR(topleft[1], topleft[2])),
		abs(topright[3]-PREDICT_LINEAR(topright[1], topright[2])),
	};

	int eidx=0;
	if(errors[eidx]>errors[1])
		eidx=1;
	if(errors[eidx]>errors[2])
		eidx=2;
	if(errors[eidx]>errors[3])
		eidx=3;
	char *ptr[]={left, top, topleft, topright}, *ptr2=ptr[eidx];

	int weights[]=
	{
		0x10000/(errors[0]+1),
		0x10000/(errors[1]+1),
		0x10000/(errors[2]+1),
		0x10000/(errors[3]+1),
	};
	int sum=weights[0]+weights[1]+weights[2]+weights[3], pred;
	int p[]=
	{
		PREDICT_LINEAR(left[2], left[3]),
		PREDICT_LINEAR(top[2], top[3]),
		PREDICT_LINEAR(topleft[2], topleft[3]),
		PREDICT_LINEAR(topright[2], topright[3]),
	};
	pred=(p[0]*weights[0]+p[1]*weights[1]+p[2]*weights[2]+p[3]*weights[3]+(sum>>1))/sum;
	//int pred=PREDICT_LINEAR(p2[1], p2[2]);

	int vmin=ptr2[0], vmax=ptr2[0];
	if(vmin>ptr2[1]) vmin=ptr2[1];
	if(vmax<ptr2[1]) vmax=ptr2[1];
	if(vmin>ptr2[2]) vmin=ptr2[2];
	if(vmax<ptr2[2]) vmax=ptr2[2];
	pred=CLAMP(vmin, pred, vmax);
	return pred;
#if 0
	char cn[]=
	{
		kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride]:0,
		         ky-1>=0?buf[idx-rowlen           ]:0,
		
		kx-1>=0?buf[idx-bytestride]:0,
	};
	int pred=cn[1]+cn[2]-cn[0], vmin=cn[0], vmax=cn[0];
	if(vmin>cn[1])
		vmin=cn[1];
	if(vmin>cn[2])
		vmin=cn[2];
	if(vmax<cn[1])
		vmax=cn[1];
	if(vmax<cn[2])
		vmax=cn[2];
	pred=CLAMP(vmin, pred, vmax);//equivalent to gradient

	//pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
#endif

#if 0
	char cn[]=
	{
		kx-3>=0&&ky-3>=0?buf[idx-rowlen*3- bytestride*3  ]:0,
		kx-2>=0&&ky-3>=0?buf[idx-rowlen*3-(bytestride<<1)]:0,
		kx-1>=0&&ky-3>=0?buf[idx-rowlen*3- bytestride    ]:0,
		         ky-3>=0?buf[idx-rowlen*3                ]:0,
		kx+1<iw&&ky-3>=0?buf[idx-rowlen*3+ bytestride    ]:0,
		kx+2<iw&&ky-3>=0?buf[idx-rowlen*3+(bytestride<<1)]:0,
		kx+3<iw&&ky-3>=0?buf[idx-rowlen*3+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
		kx+3<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-1>=0?buf[idx-rowlen- bytestride*3  ]:0,
		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,
		kx+3<iw&&ky-1>=0?buf[idx-rowlen+ bytestride*3  ]:0,
		
		kx-3>=0?buf[idx- bytestride*3  ]:0,
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};

	int predx0=cn[23], predy0=cn[17];
	int predx1=PREDICT_LINEAR(cn[22], cn[23]), predy1=PREDICT_LINEAR(cn[10], cn[17]);
	int predx2=PREDICT_QUADRATIC(cn[21], cn[22], cn[23]), predy2=PREDICT_QUADRATIC(cn[3], cn[10], cn[17]);
	int pred=(predx0+predy0)>>1;
	int best=abs(predx0-predy0),
		test=abs(predx1-predy1);
	if(best>test)
		best=test, pred=(predx1+predy1)>>1;
	test=abs(predx2-predy2);
	if(best>test)
		pred=(predx2+predy2)>>1;

	//if(pred==(predx0+predy0)>>1)
	//	++testhist[0];
	//else if(pred==(predx1+predy1)>>1)
	//	++testhist[1];
	//else if(pred==(predx2+predy2)>>1)
	//	++testhist[2];

	if(pred<-128)
		pred=-128;
	if(pred>127)
		pred=127;

	return pred;
#endif
#if 0
	char cn[]=
	{
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,

		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,

		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	int error=
		predict_simple(cn[0], cn[1], cn[5])-cn[6]+
		predict_simple(cn[1], cn[2], cn[6])-cn[7]+
		predict_simple(cn[4], cn[3], cn[9])-cn[8]+
		predict_simple(cn[5], cn[6], cn[10])-cn[11];
	int pred=predict_simple(cn[7], cn[8], cn[11]);
	pred=pred-(pred*error>>9);//10

	return pred;
#endif

	//if(pred<-64)
	//	pred=-64;
	//if(pred>64)
	//	pred=64;
	//return pred;
}
void pred_diff2d_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				if(kx==(iw>>1)&&ky==(ih>>1))//
					kx=iw>>1;

				char pred=predict_diff2d(buf, iw, kx, ky, idx, bytestride, rowlen);
#if 0
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=left+top-topleft;
				//if(kx||ky)
				//	pred-=128;
#endif
				buf[idx]-=pred;
			}
		}
	}
}
void pred_diff2d_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_diff2d(buf, iw, kx, ky, idx, bytestride, rowlen);
#if 0
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=left+top-topleft;
				//if(kx||ky)
				//	pred-=128;
#endif
				buf[idx]+=pred;
			}
		}
	}
}

void pred_hpf_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=(left+top)>>1;

				pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				buf[idx]-=pred;
			}
		}
	}
}
void pred_hpf_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=(left+top)>>1;

				pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				buf[idx]+=pred;
			}
		}
	}
}

int  predict_median(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char cn[]=
	{
	//	kx-3>=0&&ky-3>=0?buf[idx-(rowlen*3)- bytestride*3  ]:0,
	//	kx-2>=0&&ky-3>=0?buf[idx-(rowlen*3)-(bytestride<<1)]:0,
	//	kx-1>=0&&ky-3>=0?buf[idx-(rowlen*3)- bytestride    ]:0,
	//	         ky-3>=0?buf[idx-(rowlen*3)                ]:0,
	//	kx+1<iw&&ky-3>=0?buf[idx-(rowlen*3)+ bytestride    ]:0,
	//	kx+2<iw&&ky-3>=0?buf[idx-(rowlen*3)+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-3>=0?buf[idx-(rowlen*3)+ bytestride*3  ]:0,
		
	//	kx-3>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
	//	kx-3>=0&&ky-1>=0?buf[idx-rowlen- bytestride*3  ]:0,
		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-1>=0?buf[idx-rowlen+ bytestride*3  ]:0,
		
	//	kx-3>=0?buf[idx- bytestride*3  ]:0,
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	for(int k=1;k<COUNTOF(cn);++k)//insertion sort
	{
		char val=cn[k];
		int L=0, R=k-1;
		int idx=-1;
		while(L<=R)
		{
			int mid=(L+R)>>1;
			if(cn[mid]<val)
				L=mid+1;
			else if(cn[mid]>val)
				R=mid-1;
			else
			{
				idx=mid;
				break;
			}
		}
		if(idx==-1)
			idx=L+(L<k&&cn[L]<val);
		if(idx<k)
		{
			char temp=cn[k];
			memmove(cn+idx+1, cn+idx, (size_t)k-idx);
			cn[idx]=temp;
		}
	}
	int pred=(cn[(COUNTOF(cn)>>1)-1]+cn[COUNTOF(cn)>>1])>>1;

	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_median_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char pred=predict_median(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_median_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_median(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}

//double bestslope=0;
int get_mean(char *buf, int iw, int x1, int y1, int x2, int y2)
{
	int sum=0, count=-1;

	int dx=abs(x2-x1), dy=abs(y2-y1), xa, ya, xb, yb;
	int error, inc, cmp;
	if(dx>=dy)//horizontal & 45 degrees
	{
		error=dx>>1;//initialize yerror fraction with half the denominator
		if(x1<x2)
			xa=x1, ya=y1, xb=x2, yb=y2;
		else
			xa=x2, ya=y2, xb=x1, yb=y1;
		inc=ya<=yb?1:-1;
		for(int kx=xa, ky=ya;kx<=xb;++kx)
		{
			//if(kx<0||kx>=bw||ky<0||ky>=bh)//
			//	return;//
			if(kx<0||kx>=iw||ky<0)
				break;
			if(count>=0)
				sum+=buf[iw*ky+kx];
			++count;
			//buffer[bw*ky+kx]=color;

			error+=dy;//add slope to fraction 'yerror'
			cmp=-(error>=dx);	//if fraction >1 then:
			ky+=inc&cmp;		//	...increment y
			error-=dx&cmp;		//	...and subtract 1 from yerror fraction
		}
	}
	else//vertical & 45 degrees (steep)
	{
		error=dy>>1;//initialize xerror fraction with half the denominator
		if(y1<y2)
			xa=x1, ya=y1, xb=x2, yb=y2;
		else
			xa=x2, ya=y2, xb=x1, yb=y1;
		inc=xa<=xb?1:-1;
		for(int ky=ya, kx=xa;ky<=yb;++ky)
		{
			//if(kx<0||kx>=bw||ky<0||ky>=bh)//
			//	return;//
			if(kx<0||kx>=iw||ky<0)
				break;
			if(count>=0)
				sum+=buf[iw*ky+kx];
			++count;
			//buffer[bw*ky+kx]=color;

			error+=dx;//add invslope to fraction 'xerror'
			cmp=-(error>=dy);	//if fraction >1 then:
			kx+=inc&cmp;		//	...increment x
			error-=dy&cmp;		//	...and subtract 1 from xerror fraction
		}
	}
	if(count<=0)
		return 0;
	sum=(int)(((long long)sum<<16)/count);
	return sum;
#if 0
	int error=2*sx-sy;
	int x=kx, y=ky;
	if(error>0)
	{
		--x;
		error-=sy<<1;
	}
	error+=sx<<1;
	--y;
	for(;;--y)//Bressenham's line algorithm
	{
		if(x<0||y<0||x>=iw||count>=distance)
			break;

		sum+=buf[iw*y+x];
		++count;

		if(error>0)
		{
			--x;
			error-=sy<<1;
		}
		error+=sx<<1;
	}
	if(!count)
		return 0;
	sum=(int)(((long long)sum<<16)/count);
	return sum;
#endif
}
#if 0
int get_rmse(char *a1, char *a2, int count)
{
	long long sum=0;
	for(int k=0;k<count;++k)
		sum+=a1[k]*a2[k];
	sum=(long long)(sqrt((double)sum/count)*0x10000);
	return (int)sum;
}
#endif
char predict_slope(char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	if(!kx||!ky)
		return 0;
	int distance=10;
	int besterror=0x7FFFFFFF, bestslope1=0, bestslope2=0;
	for(int angle=1;angle<16-1;++angle)
	{
		double th=M_PI*angle/16;
		double rx=distance*cos(th), ry=distance*sin(th);
		int mean=get_mean(buf, iw, kx-1, ky, kx-1+(int)round(rx), ky-(int)round(ry));
		int error=abs((buf[iw*ky+kx-1]<<16)-mean);
		if(besterror>error)
			besterror=error, bestslope1=angle;
	}
	besterror=0x7FFFFFFF;
	for(int angle=1;angle<16-1;++angle)
	{
		double th=M_PI*angle/16;
		double rx=distance*cos(th), ry=distance*sin(th);
		int mean=get_mean(buf, iw, kx, ky-1, kx+(int)round(rx), ky-1-(int)round(ry));
		int error=abs((buf[iw*ky+kx-1]<<16)-mean);
		if(besterror>error)
			besterror=error, bestslope2=angle;
	}
	double th=M_PI*(bestslope1+bestslope2)/32;
	double rx=distance*cos(th), ry=distance*sin(th);
	int pred=get_mean(buf, iw, kx, ky, kx+(int)round(rx), ky-(int)round(ry));
	pred=(pred+(1<<15))>>16;
	return pred;

#if 0
	char cn[]=
	{
		kx-3>=0&&ky-3>=0?buf[idx-(rowlen*3)- bytestride*3  ]:0,//0
		kx-2>=0&&ky-3>=0?buf[idx-(rowlen*3)-(bytestride<<1)]:0,
		kx-1>=0&&ky-3>=0?buf[idx-(rowlen*3)- bytestride    ]:0,
		         ky-3>=0?buf[idx-(rowlen*3)                ]:0,
		kx+1<iw&&ky-3>=0?buf[idx-(rowlen*3)+ bytestride    ]:0,
		kx+2<iw&&ky-3>=0?buf[idx-(rowlen*3)+(bytestride<<1)]:0,
		kx+3<iw&&ky-3>=0?buf[idx-(rowlen*3)+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,//7
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
		kx+3<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-1>=0?buf[idx-rowlen- bytestride*3  ]:0,//14
		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,
		kx+3<iw&&ky-1>=0?buf[idx-rowlen+ bytestride*3  ]:0,
		
		kx-3>=0?buf[idx- bytestride*3  ]:0,//21
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	int bestslope12=0, bestrmse12=0x7FFFFFFF,
		bestslope23=0, bestrmse23=0x7FFFFFFF, rmse;
	rmse=get_rmse(cn+7  , cn+14  , 7  ); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12= 0;
	rmse=get_rmse(cn+7+1, cn+14  , 7-1); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12= 1;
	rmse=get_rmse(cn+7  , cn+14+1, 7-1); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12=-1;
	rmse=get_rmse(cn+7+2, cn+14  , 7-2); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12= 2;
	rmse=get_rmse(cn+7  , cn+14+2, 7-2); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12=-2;
	rmse=get_rmse(cn+7+3, cn+14  , 7-3); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12= 3;
	rmse=get_rmse(cn+7  , cn+14+3, 7-3); if(bestrmse12>rmse)bestrmse12=rmse, bestslope12=-3;
	
	rmse=get_rmse(cn    , cn+ 7  , 7  ); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23= 0;
	rmse=get_rmse(cn  +1, cn+ 7  , 7-1); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23= 1;
	rmse=get_rmse(cn    , cn+ 7+1, 7-1); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23=-1;
	rmse=get_rmse(cn  +2, cn+ 7  , 7-2); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23= 2;
	rmse=get_rmse(cn    , cn+ 7+2, 7-2); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23=-2;
	rmse=get_rmse(cn  +3, cn+ 7  , 7-3); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23= 3;
	rmse=get_rmse(cn    , cn+ 7+3, 7-3); if(bestrmse23>rmse)bestrmse23=rmse, bestslope23=-3;

	bestslope+=(bestslope12*3+bestslope23)*0.25;

	int pred=(cn[14+3+bestslope12]+cn[7+3+bestslope23])>>1;
	//int pred=cn[14+3+bestslope12];

	return pred;
#endif

#if 0
	if(x0>=x&&y0>=y)
		return 0;
	if(y0>=y)
		return buf[iw*y+x0];
	if(x0>=x)
		return buf[iw*y0+x];
	return (buf[iw*y+x0]+buf[iw*y0+x])>>1;
#endif
#if 0
	if(x0>=x&&y0>=y)
		return 0;
	int sum=0;
	if(y0>=y)//horizontal prediction
	{
		for(int kx=x0;kx<x;++kx)
			sum+=buf[iw*y+kx];
		sum/=x-x0;
		return sum;
	}
	if(x0>=x)//vertical prediction
	{
		for(int ky=y0;ky<y;++ky)
			sum+=buf[iw*ky+x];
		sum/=y-y0;
		return sum;
	}
	for(int k=0, kend=x-x0+y-y0;k<kend;++k)
	{
		if(k<y)
			sum+=buf[iw*k+x-k];
		else if(k<x)
			sum+=buf[iw*k+x-k];
		else
			sum+=buf[iw*y+k];
	}
#endif
}
void pred_slope_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//bestslope=0;
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				//if(kx==iw>>1&&ky==ih>>1)
				//	kx=iw>>1;
				char pred=predict_slope(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
	//bestslope/=iw*ih;
}
void pred_slope_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//bestslope=0;
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_slope(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
	//bestslope/=iw*ih;
}

double grad2_csize[GRAD2PREDCOUNT];
int grad2_hits[GRAD2PREDCOUNT];
int grad2_hist[GRAD2PREDCOUNT*256];
static double calc_csize(int *hist)
{
	int sum=0;
	for(int sym=0;sym<256;++sym)
		sum+=hist[sym];
	if(!sum)
		return 0;
	double csize=0;
	for(int sym=0;sym<256;++sym)
	{
		int freq=hist[sym];
		if(freq)
		{
			double p=(double)freq/sum;
			csize-=p*log2(p);
		}
	}
	return 8/csize;
}
int  pred_grad(char top, char left, char topleft)
{
	char vmin, vmax;
	int pred;
	if(top<left)
		vmin=top, vmax=left;
	else
		vmin=left, vmax=top;
	if(topleft<vmin)
		pred=vmax;
	else if(topleft>vmax)
		pred=vmin;
	else
		pred=top+left-topleft;
	return pred;
}
int	 pred_paeth(char top, char left, char topleft)
{
	int p=top+left-topleft;
	int closest=top;
	if(abs(closest-p)>abs(left-p))
		closest=left;
	if(abs(closest-p)>abs(topleft-p))
		closest=topleft;
	return closest;
}
int  predict_grad2(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen, char *curr, int is_fwd)
{
	char
		ctl3     =kx-3>=0&&ky-3>=0?buf[idx-rowlen*3-bytestride*3]:0,
		ct3      =         ky-3>=0?buf[idx-rowlen*3             ]:0,
		ctr3     =kx+3<iw&&ky-3>=0?buf[idx-rowlen*3+bytestride*3]:0,

		ctltl    =kx-2>=0&&ky-2>=0?buf[idx-rowlen*2-bytestride*2]:0,
		cttl     =kx-2>=0&&ky-2>=0?buf[idx-rowlen*2-bytestride  ]:0,
		ctt      =         ky-2>=0?buf[idx-rowlen*2             ]:0,
		ctrtr    =kx+2<iw&&ky-2>=0?buf[idx-rowlen*2+bytestride*2]:0,

		ctll     =kx-2>=0&&ky-1>=0?buf[idx-rowlen-bytestride*2]:0,
		ctopleft =kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride  ]:0,
		ctop     =kx  <iw&&ky-1>=0?buf[idx-rowlen             ]:0,
		ctopright=kx+1<iw&&ky-1>=0?buf[idx-rowlen+bytestride  ]:0,
		
		cl3      =kx-3>=0?buf[idx-bytestride*3]:0,
		cll      =kx-2>=0?buf[idx-bytestride*2]:0,
		cleft    =        buf[idx-bytestride  ]  ,
		ccurr    =kx  <iw?buf[idx             ]:0;
	int
		g45tl=ctopleft-ctltl,
		gxtl=ctop-ctopleft,
		gyt=ctop-ctt,
		gxtr=ctop-ctopright,
		g45tr=ctopright-ctrtr,
		gxl=cleft-cll,
		gyl=cleft-ctopleft,
		g270tl=ctopleft-ctt,
		g270tr=ctopright-ctt;

	int pred[GRAD2PREDCOUNT]={0};

	pred[0]=pred_grad(ctop, cleft, ctopleft);//0: grad predictor
	//if((gyl<0)!=(gxtl<0))
	//	pred[0]=ctt+gyl+gxtl;
	//else if(g270tl<0)
	//	pred[0]=ctt+MINVAR(gyl, gxtl);
	//else
	//	pred[0]=ctt+MAXVAR(gyl, gxtl);

	pred[1]=pred_grad(ctopleft, ctopright, ctt);//1: grad45 predictor
	//if((g270tl<0)!=(g270tr<0))
	//	pred[1]=ctt+g270tl+g270tr;
	//else if(g270tl<0)
	//	pred[1]=ctt+MINVAR(g270tl, g270tr);
	//else
	//	pred[1]=ctt+MAXVAR(g270tl, g270tr);
						
	if((gxl<0)!=(gyt<0))//2: path predictor
		pred[2]=(cleft+gxl+ctop+gyt)>>1;
	else if(gxl<0)
		pred[2]=MINVAR(ctop, cleft);
	else
		pred[2]=MAXVAR(ctop, cleft);

	if((g45tl<0)!=(g45tr<0))//3: path45 predictor
		pred[3]=(ctopleft+g45tl+ctopright+g45tr)>>1;
	else if(g45tl<0)
		pred[3]=MINVAR(ctopleft, ctopright);
	else
		pred[3]=MAXVAR(ctopleft, ctopright);

	int temp;
	if(gyl<0)//4: gamma predictor
	{
		if(gxtl<0)
		{
			if(gxtr<0)	//hole
				pred[4]=MINVAR(ctop, cleft)-((gxtl+gxtr)>>1);
			else		//bottom-right descends
				temp=ctopleft+ctopright, pred[4]=MINVAR(temp, cleft<<1)>>1;
		}
		else
		{
			if(gxtr<0)	//bottom-left descends
				pred[4]=(ctopleft+ctopright+ctop*2+cleft*2)/6;
			else		//roof, bottom descends
				pred[4]=((ctopleft+ctopright)>>1)+gyl;
		}
	}
	else
	{
		if(gxtl<0)
		{
			if(gxtr<0)	//valley, top descends
				pred[4]=((ctopleft+ctopright)>>1)+gyl;
			else		//top-right descends
				pred[4]=(ctopleft+ctopright+ctop*2+cleft*2)/6;
		}
		else
		{
			if(gxtr<0)	//top-left descends
				temp=ctopleft+ctopright, pred[4]=MAXVAR(temp, cleft<<1)>>1;
			else		//peak
				pred[4]=MAXVAR(ctop, cleft)+((gxtl+gxtr)>>1);
		}
	}

	int d1, d2;//5: select predictor
	d1=abs(cll-((cl3+cleft)>>1));                        pred[5]=(cleft+(cleft<<1)-cll)>>1;//av(const, lin) is the best among the 4 select versions
	d2=abs(ctltl-((ctl3+ctopleft)>>1));  if(d2<d1)d1=d2, pred[5]=(ctopleft+(ctopleft<<1)-ctltl)>>1;
	d2=abs(ctt-((ct3+ctop)>>1));         if(d2<d1)d1=d2, pred[5]=(ctop+(ctop<<1)-ctt)>>1;
	d2=abs(ctrtr-((ctr3+ctopright)>>1)); if(d2<d1)d1=d2, pred[5]=(ctopright+(ctopright<<1)-ctrtr)>>1;

	//d1=abs(cll-((cl3+cleft)>>1));                        pred[5]=(cleft<<1)-cll;
	//d2=abs(ctltl-((ctl3+ctopleft)>>1));  if(d2<d1)d1=d2, pred[5]=(ctopleft<<1)-ctltl;
	//d2=abs(ctt-((ct3+ctop)>>1));         if(d2<d1)d1=d2, pred[5]=(ctop<<1)-ctt;
	//d2=abs(ctrtr-((ctr3+ctopright)>>1)); if(d2<d1)d1=d2, pred[5]=(ctopright<<1)-ctrtr;
	
	//d1=abs(cll-((cl3+cleft)>>1));                        pred[5]=cleft;
	//d2=abs(ctltl-((ctl3+ctopleft)>>1));  if(d2<d1)d1=d2, pred[5]=ctopleft;
	//d2=abs(ctt-((ct3+ctop)>>1));         if(d2<d1)d1=d2, pred[5]=ctop;
	//d2=abs(ctrtr-((ctr3+ctopright)>>1)); if(d2<d1)d1=d2, pred[5]=ctopright;
#if 0
	int param[]=//select predictor v2
	{
		abs(cll  -((cl3 +cleft    )>>1)),//0
		abs(ctt  -((ct3 +ctop     )>>1)),//1
		abs(ctltl-((ctl3+ctopleft )>>1)),//2
		abs(ctrtr-((ctr3+ctopright)>>1)),//3
	};
	int den[3]=
	{
		param[0]+param[1],
		param[2]+param[3],
	};
	den[2]=den[0]+den[1];
	int pred90=den[0]?(((cleft<<1)-cll)*param[1]+((ctop<<1)-ctt)*param[0])/den[0]:(cleft+ctop)>>1,
		pred45=den[1]?(((ctopleft<<1)-ctltl)*param[3]+((ctopright<<1)-ctrtr)*param[2])/den[1]:(ctopleft+ctopright)>>1;
	pred[5]=den[2]?(pred90*den[0]+pred45*den[1])/den[2]:(pred90*5+pred45*3)>>3;
#endif


	int g2left=pred_grad(ctopleft, cll, ctll), g2top=pred_grad(ctt, ctopleft, cttl);//6: grad2
	int g2error=(abs(cleft-g2left)+abs(ctop-g2top))>>1;
	pred[6]=pred[0]*(255-g2error)/255;

	//pred[7]=ctop+cleft-ctopleft;//youchat
	//if(pred[7]<0)
	//	pred[7]=ctop;
	//else
	//{
	//	pred[7]=cleft;
	//	if(abs(ctop-cleft)<=abs(ctop-ctopleft))
	//		pred[7]=ctop;
	//}

	pred[7]=(pred[2]+pred[3]+pred[4])/3;
	
	int pred2=0;
#if 1
	int weights[]=
	{
		405214, 238662, 113409, 76552, 107912, 105129, 94864, 37906
		//489270, 244357, 116299, 79014, 109588, 107769, 33351
	}, sum=0;
	for(int k=0;k<GRAD2PREDCOUNT;++k)
	{
		pred2+=pred[k]*weights[k];
		sum+=weights[k];
	}
	pred2/=sum;
#else
	for(int k=0;k<GRAD2PREDCOUNT;++k)
		pred2+=pred[k];
	pred2/=GRAD2PREDCOUNT;
#endif
	//int pred2=(pred[0]*247756+pred[1]*279346+pred[2]*315392+pred[3]*152132+pred[4]*185022)/1179648;

	pred2=CLAMP(customparam_clamp[0], pred2, customparam_clamp[1]);

	char val;
	if(is_fwd)
		val=*curr, *curr-=pred2;
	else
		*curr+=pred2, val=*curr;

	int delta, winner=0, winner_se=0;
	for(int k=0;k<GRAD2PREDCOUNT;++k)
	{
		delta=val-pred[k];
		if(!((idx+3)&3))
			++grad2_hist[k<<8|(delta+128)];
		delta*=delta;
		//grad2_rmse[k]+=delta;
		if(!k||winner_se>delta)
			winner_se=delta, winner=k;
	}
	++grad2_hits[winner];

	//*curr=winner<<5|val>>3;

	//if(!((idx+3)&3))//
	//	*curr=winner*255/GRAD2PREDCOUNT;
	//else
	//	*curr=val;

	//*curr=abs(val-pred[4]);//
	//*curr=abs(val-pred[kx*5/iw]);//

	return pred2;
}
void pred_grad2_finish(int res)
{
	for(int k=0;k<GRAD2PREDCOUNT;++k)
		grad2_csize[k]=calc_csize(grad2_hist+(k<<8));
	//	grad2_rmse[k]=sqrt(grad2_rmse[k]/res);
}
void pred_grad2_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	memset(grad2_hist, 0, sizeof(grad2_hist));
	memset(grad2_hits, 0, sizeof(grad2_hits));
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_grad2(buf, iw, kx, ky, idx, bytestride, rowlen, buf+idx, 1);

				//buf[idx]-=pred;
			}
		}
	}
	pred_grad2_finish(iw*ih);
}
void pred_grad2_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	memset(grad2_hist, 0, sizeof(grad2_hist));
	memset(grad2_hits, 0, sizeof(grad2_hits));
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_grad2(buf, iw, kx, ky, idx, bytestride, rowlen, buf+idx, 0);

				//buf[idx]+=pred;
			}
		}
	}
	pred_grad2_finish(iw*ih);
}

int  predict_path(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char
		ctl3     =kx-3>=0&&ky-3>=0?buf[idx-rowlen*3-bytestride*3]:0,
		ct3      =         ky-3>=0?buf[idx-rowlen*3             ]:0,
		ctr3     =kx+3<iw&&ky-3>=0?buf[idx-rowlen*3+bytestride*3]:0,

		ctltl    =kx-2>=0&&ky-2>=0?buf[idx-rowlen*2-bytestride*2]:0,
		ctt      =         ky-2>=0?buf[idx-rowlen*2             ]:0,
		ctrtr    =kx+2<iw&&ky-2>=0?buf[idx-rowlen*2+bytestride*2]:0,

		ctopleft =kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride]:0,
		ctop     =         ky-1>=0?buf[idx-rowlen           ]:0,
		ctopright=kx+1<iw&&ky-1>=0?buf[idx-rowlen+bytestride]:0,
		
		cl3      =kx-3>=0?buf[idx-bytestride*3]:0,
		cll      =kx-2>=0?buf[idx-bytestride*2]:0,
		cleft    =kx-1>=0?buf[idx-bytestride  ]:0,
		ccurr    =        buf[idx             ]  ;
	int
		g45tl=ctopleft-ctltl,
		gxtl=ctop-ctopleft,
		gyt=ctop-ctt,
		gxtr=ctop-ctopright,
		g45tr=ctopright-ctrtr,
		gxl=cleft-cll,
		gyl=cleft-ctopleft,
		g270tl=ctopleft-ctt,
		g270tr=ctopright-ctt;
	
	int pred=0;
	if((gxl<0)!=(gyt<0))//2: path predictor
		pred+=(cleft+gxl+ctop+gyt)>>1;
	else if(gxl<0)
		pred+=MINVAR(ctop, cleft);
	else
		pred+=MAXVAR(ctop, cleft);

	//pred*=3;

	//if((g45tl<0)!=(g45tr<0))//3: path45 predictor
	//	pred+=(ctopleft+g45tl+ctopright+g45tr)>>1;
	//else if(g45tl<0)
	//	pred+=MINVAR(ctopleft, ctopright);
	//else
	//	pred+=MAXVAR(ctopleft, ctopright);

#if 0
	int param[8]=
	{
		cll  -((cl3 +cleft    )>>1),//0
		ctt  -((ct3 +ctop     )>>1),//1
		ctltl-((ctl3+ctopleft )>>1),//2
		ctrtr-((ctr3+ctopright)>>1),//3
	};
	param[4]=param[0]<0, param[0]=abs(param[0]);
	param[5]=param[1]<0, param[1]=abs(param[1]);
	param[6]=param[2]<0, param[2]=abs(param[2]);
	param[7]=param[3]<0, param[3]=abs(param[3]);
	int den[3]=
	{
		param[0]+param[1],
		param[2]+param[3],
	};
	den[2]=den[0]+den[1];
	int pred90=den[0]?(((cleft<<1)-cll)*param[1]+((ctop<<1)-ctt)*param[0])/den[0]:(cleft+ctop)>>1,
		pred45=den[1]?(((ctopleft<<1)-ctltl)*param[3]+((ctopright<<1)-ctrtr)*param[2])/den[1]:(ctopleft+ctopright)>>1;
	pred=den[2]?(pred90*den[0]+pred45*den[1])/den[2]:(pred90*5+pred45*3)>>3;
#endif

	//int d1, d2;//select predictor
	//d1=cll-((cl3+cleft)>>1);                 pred=cleft;
	//d2=ctltl-((ctl3+ctopleft)>>1);  if(d2<d1)pred=ctopleft;
	//d2=ctt-((ct3+ctop)>>1);         if(d2<d1)pred=ctop;
	//d2=ctrtr-((ctr3+ctopright)>>1); if(d2<d1)pred=ctopright;

	//d1=cll-((cl3+cleft)>>1);                 pred=(cleft+(cleft<<1)-cll)>>1;
	//d2=ctltl-((ctl3+ctopleft)>>1);  if(d2<d1)pred=(ctopleft+(ctopleft<<1)-ctltl)>>1;
	//d2=ctt-((ct3+ctop)>>1);         if(d2<d1)pred=(ctop+(ctop<<1)-ctt)>>1;
	//d2=ctrtr-((ctr3+ctopright)>>1); if(d2<d1)pred=(ctopright+(ctopright<<1)-ctrtr)>>1;

	//pred>>=2;
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_path_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_path(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_path_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_path(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}

int adagrad_hits[ADAGRADCOUNT];
double adagrad_rmse[ADAGRADCOUNT], adagrad_csize[ADAGRADCOUNT];
double adagrad_abserror[ADAGRADCOUNT], adagrad_signederror[ADAGRADCOUNT];
void pred_adaptive(char *buf, int iw, int ih, int nch, int bytestride, int fwd)
{
	memset(adagrad_hits, 0, sizeof(adagrad_hits));
	memset(grad2_hist, 0, sizeof(grad2_hist));
	memset(adagrad_abserror, 0, sizeof(adagrad_abserror));
	memset(adagrad_signederror, 0, sizeof(adagrad_signederror));
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res*bytestride);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(b2, 0, (size_t)res*bytestride);
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			int error[8]={0};
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred;
				int current_val=buf[idx];
				
				char *src=fwd?buf:b2;
				char
					ctl3     =kx-3>=0&&ky-3>=0?src[idx-rowlen*3-bytestride*3]:0,
					ct3      =         ky-3>=0?src[idx-rowlen*3             ]:0,
					ctr3     =kx+3<iw&&ky-3>=0?src[idx-rowlen*3+bytestride*3]:0,

					ctltl    =kx-2>=0&&ky-2>=0?src[idx-rowlen*2-bytestride*2]:0,
					ctt      =         ky-2>=0?src[idx-rowlen*2             ]:0,
					ctrtr    =kx+2<iw&&ky-2>=0?src[idx-rowlen*2+bytestride*2]:0,

					ctopleft =kx-1>=0&&ky-1>=0?src[idx-rowlen-bytestride  ]:0,
					ctop     =kx  <iw&&ky-1>=0?src[idx-rowlen             ]:0,
					ctopright=kx+1<iw&&ky-1>=0?src[idx-rowlen+bytestride  ]:0,
					ctrr     =kx+2<iw&&ky-1>=0?src[idx-rowlen+bytestride*2]:0,
		
					cl3      =kx-3>=0?src[idx-bytestride*3]:0,
					cll      =kx-2>=0?src[idx-bytestride*2]:0,
					cleft    =kx-1>=0?src[idx-bytestride  ]:0;
				int
					g45tl=ctopleft-ctltl,
					gxtl=ctop-ctopleft,
					gyt=ctop-ctt,
					gxtr=ctop-ctopright,
					g45tr=ctopright-ctrtr,
					gxl=cleft-cll,
					gyl=cleft-ctopleft;
				
				//int gx=gxl+gxtl, gy=gyl+gyt,
				//	T=customparam_st[0];//44
				//if(gy>T+gx)
				//	pred=cleft;
				//else if(gy+T<gx)
				//	pred=ctop;
				//else
				//	pred=ctop+cleft-ctopleft;

				int type;

				//grad
#if 0
				char vmin, vmax;
				if(ctop<cleft)
					vmin=ctop, vmax=cleft;
				else
					vmin=cleft, vmax=ctop;

				//CHEAT
#if 0
				pred=vmin, type=0;
				if(abs(buf[idx]-pred)>=abs(buf[idx]-vmax))
					pred=vmax, type=1;
				int p2=ctop+cleft-ctopleft;
				if(abs(buf[idx]-pred)>=abs(buf[idx]-p2))
					pred=p2, type=2;
#endif
				if(ctopleft<vmin)
					pred=vmax, type=0;
				else if(ctopleft>vmax)
					pred=vmin, type=1;
				else
				{
					pred=ctop+cleft-ctopleft, type=2;
					//pred=pred*(255-abs(error))/255;
				}
#endif

				//linear extrapolation
#if 0
				int predleft=(cleft<<1)-cll, predtop=(ctop<<1)-ctt;
				if(!error[0]&&!error[1])
					pred=predleft;
				else if((error[0]<0)==(error[1]<0))
					pred=(predleft*abs(error[1])+predtop*abs(error[0]))/(abs(error[0])+abs(error[1])), type=0;
				else
					pred=(predleft+predtop)>>1, type=1;
#endif

				//multigrad
#if 0
				type=0;
				int prednear=pred_grad(ctop, cleft, ctopleft),
					predtop=pred_grad(ctopleft, ctopright, ctt),
					predfar=pred_grad(ctt, cll, ctltl);
				int den=abs(error[1])+abs(error[2]);
				if(den)
					pred=(predfar*abs(error[1])+predtop*abs(error[2]))/den;
				else
					pred=predtop;
				int den2=den+abs(error[0]);
				if(den2)
					pred=(pred*abs(error[0])+prednear*den)/den2;
				else
					pred=prednear;
#endif

				//gamma v2
#if 0
				int temp;
				if(gyl<0)//gamma predictor
				{
					if(gxtl<0)
					{
						if(gxtr<0)	//hole
							pred=MINVAR(ctop, cleft), type=0;
							//pred=MINVAR(ctop, cleft)+(int)(gxtl*(-0.16+customparam_st[0])+gxtr*(-0.06+customparam_st[5])), type=0;
							//pred=MINVAR(ctop, cleft)+(int)(gxtl*customparam_st[0]+gxtr*customparam_st[5]), type=0;//params: -0.16 -0.06
							//pred=MINVAR(ctop, cleft)+(int)((gxtl+gxtr)*customparam_st[0]), type=0;
							//pred=MINVAR(ctop, cleft), type=0;
							//pred=MINVAR(ctop, cleft)-(gxtl+gxtr)/4, type=0;
						else		//bottom-right descends
							pred=MINVAR(ctop, cleft), type=1;
							//pred=(int)(gyl*(0.5+customparam_st[1])+(ctopleft+gxtl*(0.51+customparam_st[6])+ctopright+gxtr)*0.5), type=1;
							//pred=(int)(gyl*(1+customparam_st[1])+(ctopleft+gxtl*(1+customparam_st[6])+ctopright+gxtr)*0.5), type=1;//params: -0.5 -0.49
							//temp=ctop*2+ctopleft+ctopright, pred=MINVAR(temp, cleft<<2)/4, type=1;
							//pred2+=(ctopleft+ctopright+ctop+cleft)*0.25f;
							//pred+=ctopright;
					}
					else
					{
						if(gxtr<0)	//bottom-left descends
							pred=ctop+cleft-ctopleft, type=2;
							//pred=(ctopleft+ctopright+ctop*2+cleft*2)/6, type=2;
							//pred=ctop+gyl, type=2;
						else		//roof, bottom descends
							pred=ctop+cleft-ctopleft, type=3;
							//pred=(ctopleft+ctopright)/2+gyl, type=3;
					}
				}
				else
				{
					if(gxtl<0)
					{
						if(gxtr<0)	//valley, top descends
							pred=ctop+cleft-ctopleft, type=4;
							//pred=(ctopleft+ctopright)/2+gyl, type=4;
							//pred=(ctopleft+ctopright)>>1, type=4;
						else		//top-right descends
							pred=ctop+cleft-ctopleft, type=5;
							//pred=(ctopleft+ctopright+ctop*2+cleft*2)/6, type=5;
							//pred2+=ctop+gyl;
					}
					else
					{
						if(gxtr<0)	//top-left descends
							pred=MAXVAR(ctop, cleft), type=6;
							//temp=ctopleft+ctopright, pred=MAXVAR(temp, cleft<<1)/2, type=6;
							//pred2+=(ctopleft+ctopright+ctop+cleft)*0.25f;
						else		//peak
							pred=MAXVAR(ctop, cleft), type=7;
							//pred=MAXVAR(ctop, cleft)+(gxtl+gxtr)/16, type=7;
					}
				}
#endif

#if 1
				int predictors[]=
				{
					pred_grad(ctop, cleft, ctopleft),
					(-2*ctt+6*ctop+3*ctopright+ctrr+cll+7*cleft+8)>>4,
					cleft,
					ctop,
					ctopleft,
					ctopright,
					(cleft<<1)-cll,
					(ctop<<1)-ctt,
					//pred_paeth(ctop, cleft, ctopleft),
				};

				type=(kx+ky+kc)&7;
				pred=predictors[type];
				//type=0;
				//pred=predictors[0];
				//for(int k=1;k<COUNTOF(predictors);++k)
				//{
				//	if(error[type]>error[k])
				//		type=k, pred=predictors[k];
				//}
#endif

				//pred=pred_grad(ctop, cleft, ctopleft)+((error[0]-error[1])>>2);

				//pred=pred*(255-abs(error))/255;
				
				//if(kc==1&&kx==(iw>>1)&&ky==(ih>>1))//
				//	kx=iw>>1;

				pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				if(fwd)
				{
					current_val=buf[idx];
					b2[idx]=buf[idx]-pred;

					//b2[idx]=type*255/7-128;
					//b2[idx]=type<<6|current_error>>2&0x3F;//
				}
				else
				{
					b2[idx]=buf[idx]+pred;
					current_val=b2[idx];
				}
				
				++adagrad_hits[type];
				int delta=current_val-pred;
				adagrad_rmse[type]+=delta*delta;
				++grad2_hist[type<<8|(delta+128)];//

				for(int k=0;k<COUNTOF(error);++k)
				{
					error[k]=current_val-predictors[k];
					adagrad_abserror[k]+=abs(error[k]);
					adagrad_signederror[k]+=error[k];
				}

				//++grad2_hist[type<<8|(current_val+128)];//

				//error[0]=current_val-prednear;
				//error[1]=current_val-predtop;
				//error[2]=current_val-predfar;

				//error[0]=current_val-predleft;
				//error[1]=current_val-predtop;

				//error[2]=error[1];
				//error[1]=error[0];
				//error[0]=current_error;

				//error=current_error+(error-current_error)*3/8;
			}
		}
	}
	for(int kc=nch;kc<bytestride;++kc)
	{
		for(int k=0;k<res;++k)
			b2[k*bytestride+kc]=buf[k*bytestride+kc];
	}
	memcpy(buf, b2, (size_t)res*bytestride);
	free(b2);
	for(int k=0;k<ADAGRADCOUNT;++k)
	{
		adagrad_rmse[k]=sqrt(adagrad_rmse[k]/adagrad_hits[k]);
		double invCR=calc_entropy(grad2_hist+(k<<8), -1)/8;
		adagrad_csize[k]=adagrad_hits[k]*invCR;
		//adagrad_abserror[k]/=adagrad_hits[k];
	}
}


int jointpredparams[]=//fixed 19.12 bit
{
	//left						top					topleft				topright
	//r      g  b   errors
	0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,//r
	0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,//g
	0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,//b

	//wtop, wleft, wtopleft, wtopright
	0x1000, 0x1000, 0x1000, 0x1000,//r
	0x1000, 0x1000, 0x1000, 0x1000,//g
	0x1000, 0x1000, 0x1000, 0x1000,//b
};
void   pred_joint(const char *src, int iw, int ih, int *params, int fwd, char *dst, int *temp_w24)
{
	const char *src2=fwd?src:dst, *error=fwd?dst:src;
	for(int ky=0;ky<ih;++ky)
	{
		int *currrow=temp_w24+(ky&1?0:iw*12), *prevrow=temp_w24+(ky&1?iw*12:0);
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			int in[]=
			{
				kx-1>=0         ?src2 [(idx   -1)<<2  ]<<12:0,//Rleft
				kx-1>=0         ?src2 [(idx   -1)<<2|1]<<12:0,//Gleft
				kx-1>=0         ?src2 [(idx   -1)<<2|2]<<12:0,//Bleft

				kx-1>=0         ?error[(idx   -1)<<2  ]<<12:0,//REleft
				kx-1>=0         ?error[(idx   -1)<<2|1]<<12:0,//GEleft
				kx-1>=0         ?error[(idx   -1)<<2|2]<<12:0,//BEleft

				kx  <iw&&ky-1>=0?src2 [(idx-iw  )<<2  ]<<12:0,//Rtop
				kx  <iw&&ky-1>=0?src2 [(idx-iw  )<<2|1]<<12:0,//Gtop
				kx  <iw&&ky-1>=0?src2 [(idx-iw  )<<2|2]<<12:0,//Btop

				kx  <iw&&ky-1>=0?error[(idx-iw  )<<2  ]<<12:0,//REtop
				kx  <iw&&ky-1>=0?error[(idx-iw  )<<2|1]<<12:0,//GEtop
				kx  <iw&&ky-1>=0?error[(idx-iw  )<<2|2]<<12:0,//BEtop

				kx-1>=0&&ky-1>=0?src2 [(idx-iw-1)<<2  ]<<12:0,//Rtopleft
				kx-1>=0&&ky-1>=0?src2 [(idx-iw-1)<<2|1]<<12:0,//Gtopleft
				kx-1>=0&&ky-1>=0?src2 [(idx-iw-1)<<2|2]<<12:0,//Btopleft

				kx-1>=0&&ky-1>=0?error[(idx-iw-1)<<2  ]<<12:0,//REtopleft
				kx-1>=0&&ky-1>=0?error[(idx-iw-1)<<2|1]<<12:0,//GEtopleft
				kx-1>=0&&ky-1>=0?error[(idx-iw-1)<<2|2]<<12:0,//BEtopleft

				kx+1<iw&&ky-1>=0?src2 [(idx-iw+1)<<2  ]<<12:0,//Rtopright
				kx+1<iw&&ky-1>=0?src2 [(idx-iw+1)<<2|1]<<12:0,//Gtopright
				kx+1<iw&&ky-1>=0?src2 [(idx-iw+1)<<2|2]<<12:0,//Btopright

				kx+1<iw&&ky-1>=0?error[(idx-iw+1)<<2  ]<<12:0,//REtopright
				kx+1<iw&&ky-1>=0?error[(idx-iw+1)<<2|1]<<12:0,//GEtopright
				kx+1<iw&&ky-1>=0?error[(idx-iw+1)<<2|2]<<12:0,//BEtopright
			};

			//if(kx==(iw>>1)&&ky==(ih>>1))//
			//	kx=iw>>1;

			int weights[12], w;
			for(int k=0;k<12;++k)
			{
				w=(ky-1>=0?prevrow[12*kx+k]:0)+(ky-1>=0&&kx+1<iw?prevrow[12*(kx+1)+k]:0)+(ky-1>=0&&kx-1>=0?prevrow[12*(kx-1)+k]:0);
				weights[k]=(int)(((long long)params[72+k]<<24)/(w+1));
			}

			int vmin[]={in[0], in[1], in[2]}, vmax[]={in[0], in[1], in[2]};
			for(int k=6;k<24;k+=6)
			{
				if(vmin[0]>in[k  ])vmin[0]=in[k  ];
				if(vmax[0]<in[k  ])vmax[0]=in[k  ];
				if(vmin[1]>in[k+1])vmin[1]=in[k+1];
				if(vmax[1]<in[k+1])vmax[1]=in[k+1];
				if(vmin[2]>in[k+2])vmin[2]=in[k+2];
				if(vmax[2]<in[k+2])vmax[2]=in[k+2];
			}

			long long cand[12]={0};//12 candidate predictions
			for(int k=0;k<6;++k)
			{
				cand[ 0]+=params[k   ]*in[k   ];
				cand[ 1]+=params[k+ 6]*in[k+ 6];
				cand[ 2]+=params[k+12]*in[k+12];
				cand[ 3]+=params[k+18]*in[k+18];

				cand[ 4]+=params[k+24]*in[k   ];
				cand[ 5]+=params[k+30]*in[k+ 6];
				cand[ 6]+=params[k+36]*in[k+12];
				cand[ 7]+=params[k+42]*in[k+18];

				cand[ 8]+=params[k+48]*in[k   ];
				cand[ 9]+=params[k+54]*in[k+ 6];
				cand[10]+=params[k+60]*in[k+12];
				cand[11]+=params[k+66]*in[k+18];
			}

			int fpred[3];//3 final predictions
			int sum[]=
			{
				weights[ 0]+weights[ 1]+weights[ 2]+weights[ 3],
				weights[ 4]+weights[ 5]+weights[ 6]+weights[ 7],
				weights[ 8]+weights[ 9]+weights[10]+weights[11],
			};
			for(int k=0;k<3;++k)
			{
				if(sum[k])
				{
					w=(int)((weights[4*k]*cand[4*k]+weights[4*k+1]*cand[4*k+1]+weights[4*k+2]*cand[4*k+2]+weights[4*k+3]*cand[4*k+3])>>12);
					fpred[k]=(int)((w+(sum[k]>>1)-1)/sum[k]);
				}
				else
					fpred[k]=(int)(cand[4*k]>>12);
			}
			//for(int k=0;k<3;++k)
			//	fpred[k]=CLAMP(vmin[k], fpred[k], vmax[k]);
			int curr[3];
			if(fwd)
			{
				curr[0]=src[idx<<2  ]<<12, dst[idx<<2  ]=src[idx<<2  ]-((fpred[0]+0x7FF)>>12);
				curr[1]=src[idx<<2|1]<<12, dst[idx<<2|1]=src[idx<<2|1]-((fpred[1]+0x7FF)>>12);
				curr[2]=src[idx<<2|2]<<12, dst[idx<<2|2]=src[idx<<2|2]-((fpred[2]+0x7FF)>>12);
			}
			else
			{
				dst[idx<<2  ]=src[idx<<2  ]+((fpred[0]+0x7FF)>>12), curr[0]=dst[idx<<2  ]<<12;
				dst[idx<<2|1]=src[idx<<2|1]+((fpred[1]+0x7FF)>>12), curr[1]=dst[idx<<2|1]<<12;
				dst[idx<<2|2]=src[idx<<2|2]+((fpred[2]+0x7FF)>>12), curr[2]=dst[idx<<2|2]<<12;
			}
			
			for(int k=0;k<12;++k)
			{
				int e=abs(curr[k>>2]-(int)(cand[k]>>12));
				currrow[12*kx+k]=e;
				if(kx+1<iw)
					prevrow[12*(kx+1)+k]+=e;
			}
		}
	}
}
double pred_joint_calcsize(const char *src, int iw, int ih, int *params, int *temp, char *dst, int *hist, int loud, int iter, int itcount)
{
	int res=iw*ih;
	pred_joint(src, iw, ih, params, 1, dst, temp);
	if(loud)
		addhalf((unsigned char*)dst, iw, ih, 3, 4);
	//addhalf((unsigned char*)dst+kc, iw, ih, 1, 4);

	double csize=0;
	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(dst+0, (ptrdiff_t)res<<2, 4, hist);

		double entropy=0;
		for(int k=0;k<256;++k)
		{
			int freq=hist[k];
			if(freq)
			{
				double p=(double)freq/res;
				p*=0x10000-255;
				++p;
				p/=0x10000;
				entropy-=p*log2(p);
			}
		}
		double invCR=entropy/8;
		csize+=res*invCR;
	}
	if(loud)
	{
		double CR=(3*iw*ih)/csize;
		set_window_title("[%d/%d] CR %lf", iter+1, itcount, CR);
		io_render();
	}
	return csize;
}
void vec_add(int *a, const int *b, int count)
{
	for(int k=0;k<count;++k)
		a[k]+=b[k];
}
void vec_sub(int *a, const int *b, int count)
{
	for(int k=0;k<count;++k)
		a[k]-=b[k];
}
void vec_initrand(int *vec, int count, int vmin, int vmax)
{
	for(int k=0;k<count;++k)
		vec[k]=vmin+rand()%(vmax+1-vmin);
}
void   pred_joint_optimize(const char *src, int iw, int ih, int *params, int step, char *dst, int loud)
{
	int *temp=(int*)malloc(iw*24LL*sizeof(int));
	int *hist=(int*)malloc(256*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}

	ArrayHandle title;
	STR_ALLOC(title, 1024);
	get_window_title(title->data, 1024);
	unsigned char *image0=image;
	image=(unsigned char*)dst;
	
	int res=iw*ih;
	double csize0, csize, csize00;

	srand((unsigned)__rdtsc());
	//int idx=rand()%PREDJOINT_NPARAMS;
	//static int idx=0;
	//for(int idx=0;idx<PREDJOINT_NPARAMS;++idx)
	//{
		int p0[PREDJOINT_NPARAMS], vec[PREDJOINT_NPARAMS], pos, count, start;
		memcpy(p0, params, sizeof(p0));
		//int p0=params[idx];
		csize=pred_joint_calcsize(src, iw, ih, params, temp, dst, hist, loud, 0, PREDJOINT_NPARAMS);
		csize00=csize;

		int subit;
		for(subit=0;subit<20;++subit)
		{
			pos=rand()%PREDJOINT_NPARAMS;
			count=rand()%(PREDJOINT_NPARAMS>>2);
			if(pos-(count>>1)<0)
			{
				if(pos)
					count=pos<<1;
				else
					pos=1, count=2;
			}
			if(pos-(count>>1)+count>PREDJOINT_NPARAMS)
			{
				if(pos<PREDJOINT_NPARAMS)
					count=(PREDJOINT_NPARAMS-pos)<<1;
				else
					pos=PREDJOINT_NPARAMS-1, count=2;
			}
			start=pos-(count>>1);
			vec_initrand(vec+start, count, -0xFF, 0xFF);
			//vec_initrand(vec, PREDJOINT_NPARAMS, -0xF, 0xF);
			//step=rand()&0xFF;
			//step+=!step;

			csize0=csize;
			vec_add(params+start, vec+start, count);
			//vec_add(params, vec, PREDJOINT_NPARAMS);
			//params[idx]+=step;
			csize=pred_joint_calcsize(src, iw, ih, params, temp, dst, hist, loud, 0, PREDJOINT_NPARAMS);
			if(csize>csize0)//cancel last change and break
			{
				vec_sub(params+start, vec+start, count);
				//params[idx]-=step;
				csize=csize0;
				//break;
			}
		}
		
		for(subit=0;subit<20;++subit)
		{
			pos=rand()%PREDJOINT_NPARAMS;
			count=rand()%(PREDJOINT_NPARAMS>>2);
			if(pos-(count>>1)<0)
			{
				if(pos)
					count=pos<<1;
				else
					pos=1, count=2;
			}
			if(pos-(count>>1)+count>PREDJOINT_NPARAMS)
			{
				if(pos<PREDJOINT_NPARAMS)
					count=(PREDJOINT_NPARAMS-pos)<<1;
				else
					pos=PREDJOINT_NPARAMS-1, count=2;
			}
			start=pos-(count>>1);
			vec_initrand(vec+start, count, -0xFF, 0xFF);

			csize0=csize;
			vec_sub(params+start, vec+start, count);
			//params[idx]-=step;
			csize=pred_joint_calcsize(src, iw, ih, params, temp, dst, hist, loud, 0, PREDJOINT_NPARAMS);
			if(csize>csize0)
			{
				vec_add(params+start, vec+start, count);
				//params[idx]+=step;
				csize=csize0;
				//break;
			}
		}
		if(csize>csize00)//prevent CR from worsening
		{
			memcpy(params, p0, sizeof(p0));
			//params[idx]=p0;
			csize=csize00;
		}
	//}
	//idx=(idx+1)%PREDJOINT_NPARAMS;
	
	set_window_title("%s", title->data);
	image=(unsigned char*)image0;
	free(hist);
	free(temp);
}
void   pred_joint_apply(char *buf, int iw, int ih, int *allparams, int fwd)
{
	int res=iw*ih;
	int *temp=(int*)malloc(iw*24LL*sizeof(int));
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_joint(buf, iw, ih, allparams, fwd, buf2, temp);

	for(int k=0;k<res;++k)
	{
		buf[k<<2  ]=buf2[k<<2  ];
		buf[k<<2|1]=buf2[k<<2|1];
		buf[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}


double pw2_errors[PW2_NPRED]={0};//
short pw2_params[PW2_NPARAM*3]=
{
	//0
	
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,-0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C, 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,

	// 0x007D, 0x0040, 0x0039, 0x004A, 0x0007, 0x0003, 0x001D, 0x0007, 0x0000, 0x0007, 0x0002, 0x000F, 0x001B, 0x0018, 0x000A, 0x0008, 0x001C,-0x0008, 0x0004, 0x0005, 0x0006,-0x0022,-0x003B,-0x0041,-0x00C8,-0x0040,-0x0085,-0x0050, 0x0060,
	// 0x007A, 0x00F9, 0x0165, 0x00C2, 0x0036, 0x0100, 0x0054, 0x0000,-0x0081,-0x0078, 0x0020, 0x004D,-0x0010, 0x0028, 0x00BD, 0x009D, 0x0020,-0x0082,-0x003F, 0x0060, 0x002A, 0x0161,-0x004E,-0x001D, 0x0123,-0x0008,-0x0080, 0x0020, 0x003C,
	// 0x0010, 0x0039, 0x002E, 0x0037, 0x000E,-0x0010, 0x0014, 0x0008,-0x0007,-0x001C, 0x0074, 0x0019, 0x0010, 0x001B, 0x000D, 0x0047, 0x000A, 0x001C, 0x0008, 0x0004, 0x0023,-0x0012,-0x0156,-0x0074,-0x00A0,-0x0002,-0x0088,-0x0060, 0x0102,
};
static int pred_w2_paeth2(int T, int L, int TL, int TR)
{
	int p=T+L-TL, closest=T;
	if(abs(closest-p)>abs(L-p))
		closest=L;
	if(abs(closest-p)>abs(TL-p))
		closest=TL;
	if(abs(closest-p)>abs(TR-p))
		closest=TR;
	return closest;
}
static int pred_w2_select(int T, int L, int TL)
{
	int p=T+L-TL, pT=abs(p-T), pL=abs(p-L);
	return pT<pL?L:T;
}
static pred_w2_cgrad(int T, int L, int TL)
{
	int vmin, vmax, grad;

	if(T<L)
		vmin=T, vmax=L;
	else
		vmin=L, vmax=T;
	grad=T+L-TL;
	grad=CLAMP(vmin, grad, vmax);
	return grad;
}
static int clamp4(int p, int a, int b, int c, int d)
{
	int vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmin>c)vmin=c;
	if(vmin>d)vmin=d;
	if(vmax<b)vmax=b;
	if(vmax<c)vmax=c;
	if(vmax<d)vmax=d;
	p=CLAMP(vmin, p, vmax);
	return p;
}
static int clip(int x)
{
	x=CLAMP(-128, x, 127);
	return x;
}
void pred_w2_prealloc(const char *src, int iw, int ih, int kc, short *params, int fwd, char *dst, int *temp)//temp is (PW2_NPRED+1)*2w
{
	int errorbuflen=iw<<1, rowlen=iw<<2;
	int *error=temp, *pred_errors[PW2_NPRED];
	for(int k=0;k<PW2_NPRED;++k)
		pred_errors[k]=temp+errorbuflen*(k+1);
	int idx=kc;

	const char *src2=fwd?src:dst;
	memset(pw2_errors, 0, sizeof(pw2_errors));
	for(int ky=0;ky<ih;++ky)
	{
		//int pred_left=0, pred_left2=0;

		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			//           T3
			//   T2L2    T2    T2R2
			//        TL T  TR TR2
			//L3 L2   L  X
			int
				cT6  =         ky-6>=0?src2[idx-rowlen*6   ]<<8:0,
				cT5  =         ky-5>=0?src2[idx-rowlen*5   ]<<8:0,

				cT4L3=kx-3>=0&&ky-4>=0?src2[idx-rowlen*4-12]<<8:0,
				cT4  =         ky-4>=0?src2[idx-rowlen*4   ]<<8:0,
				cT4R3=kx+3<iw&&ky-4>=0?src2[idx-rowlen*4+12]<<8:0,
				
				cT3L5=kx-5>=0&&ky-3>=0?src2[idx-rowlen*3-20]<<8:0,
				cT3L4=kx-4>=0&&ky-3>=0?src2[idx-rowlen*3-16]<<8:0,
				cT3L2=kx-2>=0&&ky-3>=0?src2[idx-rowlen*3- 8]<<8:0,
				cT3L =kx-1>=0&&ky-3>=0?src2[idx-rowlen*3- 4]<<8:0,
				cT3  =         ky-3>=0?src2[idx-rowlen*3   ]<<8:0,
				cT3R =kx+1<iw&&ky-3>=0?src2[idx-rowlen*3+ 4]<<8:0,
				cT3R2=kx+2<iw&&ky-3>=0?src2[idx-rowlen*3+ 8]<<8:0,
				cT3R3=kx+3<iw&&ky-3>=0?src2[idx-rowlen*3+12]<<8:0,
				cT3R4=kx+4<iw&&ky-3>=0?src2[idx-rowlen*3+16]<<8:0,
				
				cT2L3=kx-3>=0&&ky-2>=0?src2[idx-rowlen*2-12]<<8:0,
				cT2L2=kx-2>=0&&ky-2>=0?src2[idx-rowlen*2- 8]<<8:0,
				cT2L =kx-1>=0&&ky-2>=0?src2[idx-rowlen*2- 4]<<8:0,
				cT2  =         ky-2>=0?src2[idx-rowlen*2   ]<<8:0,
				cT2R =kx+1<iw&&ky-2>=0?src2[idx-rowlen*2+ 4]<<8:0,
				cT2R2=kx+2<iw&&ky-2>=0?src2[idx-rowlen*2+ 8]<<8:0,
				cT2R3=kx+3<iw&&ky-2>=0?src2[idx-rowlen*2+12]<<8:0,
				cT2R4=kx+4<iw&&ky-2>=0?src2[idx-rowlen*2+16]<<8:0,
				
				cTL3 =kx-3>=0&&ky-1>=0?src2[idx-rowlen  -12]<<8:0,
				cTL2 =kx-2>=0&&ky-1>=0?src2[idx-rowlen  - 8]<<8:0,
				cTL  =kx-1>=0&&ky-1>=0?src2[idx-rowlen  - 4]<<8:0,
				cT   =kx  <iw&&ky-1>=0?src2[idx-rowlen     ]<<8:0,
				cTR  =kx+1<iw&&ky-1>=0?src2[idx-rowlen  + 4]<<8:0,
				cTR2 =kx+2<iw&&ky-1>=0?src2[idx-rowlen  + 8]<<8:0,
				cTR3 =kx+3<iw&&ky-1>=0?src2[idx-rowlen  +12]<<8:0,
				cTR4 =kx+4<iw&&ky-1>=0?src2[idx-rowlen  +16]<<8:0,
				cTR5 =kx+5<iw&&ky-1>=0?src2[idx-rowlen  +20]<<8:0,
				cTR6 =kx+6<iw&&ky-1>=0?src2[idx-rowlen  +24]<<8:0,
				cTR7 =kx+7<iw&&ky-1>=0?src2[idx-rowlen  +28]<<8:0,

				cL6  =kx-6>=0         ?src2[idx         -24]<<8:0,
				cL5  =kx-5>=0         ?src2[idx         -20]<<8:0,
				cL4  =kx-4>=0         ?src2[idx         -16]<<8:0,
				cL3  =kx-2>=0         ?src2[idx         -12]<<8:0,
				cL2  =kx-2>=0         ?src2[idx         - 8]<<8:0,
				cL   =kx-1>=0         ?src2[idx         - 4]<<8:0;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c

			//if(kx==(iw>>1)&&ky==(ih>>1))//
			//	kx=iw>>1;
			
			int weights[PW2_NPRED];//fixed 23.8 bit
			for(int k=0;k<PW2_NPRED;++k)
			{
				int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
				weights[k]=(params[k]<<8)/(w+1);
			}

			//TL T TR
			//L  X
			int
				eT=ky-1>=0?error[prevrow+kx]:0,
				eL=kx-1>=0?error[currrow+kx-1]:0,
				eTL=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				eTR=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
				eT_L=eT+eL;

			//pred_left=pred_left+((cL-pred_left)*params[PW2_NPRED+11]>>8);
			//pred_left2=pred_left2+((cL-pred_left2)*abs(eL)>>16);

			//int pred_w=(cT*(0x10000-eT)+cL*(0x10000-eL)+cTL*(0x10000-eTL)+cTR*(0x10000-eTR))>>16;

			int predictions[PW2_NPRED]=//fixed 23.8 bit
			{
				//from jxl
				cT-((eTL*params[PW2_NPRED]+eT*params[PW2_NPRED+1]+eTR*params[PW2_NPRED+2]+(cT2-cT)*params[PW2_NPRED+3]+(cTL-cL)*params[PW2_NPRED+4])>>8),//k13: 1.998458 how many optimizations?
				//k13: 1.737220, 1.837736, 1.842377, 1.847278, 1.865093, 1.861586, 1.866601, 1.872176, 1.878888, 1.883146, 1.883862, 1.883863, <

				cL-((eT_L+eTL)*params[PW2_NPRED+5]>>8),
				//k13: 1.737220, 1.932766, 1.960469, 1.966698, 1.970106, 1.971381, 1.971891, 1.972187, 1.972368, 1.972492, 1.972565, 1.972580, 1.972594, 1.972602, <

				cT-((eT_L+eTR)*params[PW2_NPRED+6]>>8),
				//k13: 1.737220, 1.934911, 1.951553, 1.962075, 1.973068, 1.979630, 1.983403, 1.985205, 1.986109, 1.986488, 1.986609, 1.986677, 1.986688, 1.986689, <

				cL+cTR-cT,
				//k13: 1.909206, 1.918954, 1.946385, 1.963439, 1.970334, 1.981971, 1.983898, 1.984527, 1.985102, 1.985773, 1.986008, 1.986331, 1.986678, 1.986722, 1.986756			...1.998458
#if 1
				cL-(eL*params[PW2_NPRED+7]>>8),
				cT-(eT*params[PW2_NPRED+8]>>8),
				//k13: 1.909206, 1.922112, 1.931268, 1.954690, 1.964123, 1.978742, 1.981578, 1.984138, 1.985535, 1.986711, 1.987659, 1.988190, 1.988474, 1.988532, 1.988550
				cTL-(eTL*params[PW2_NPRED+9]>>8),
				cTR-(eTR*params[PW2_NPRED+10]>>8),
				//k13: 1.909206, 1.921490, 1.932631, 1.949766, 1.950930, 1.951645, 1.951977, 1.960758, 1.967595, 1.969669, 1.972408, 1.973050, 1.973506, 1.974268, 1.975184			...1.977183
#endif
#if 0
				pred_left,
				pred_left2,//k13: 1.977869 how many optimizations?
				(cL*params[PW2_NPRED+11]+cT*params[PW2_NPRED+12]+cTL*params[PW2_NPRED+13]+cTR*params[PW2_NPRED+14])>>8,//1.909206 -> 
				(cL*params[PW2_NPRED+11]+cT*params[PW2_NPRED+12]+cTL*params[PW2_NPRED+13]+cTR*params[PW2_NPRED+14])>>8,//1.909206 -> 
#endif
				//pred_w,//k13: 1.909206
#if 1
				//paq8px by Matt Mahoney

				//k13: 1.737220, 1.958513, 1.973267, 1.979685, 1.983374, 1.985860, 1.987622, 1.989731, 1.991147, 1.992018, 1.992707, 1.993444, 1.994374, 1.995238, 1.996056,   1.996876, 1.997423, 1.997708, 1.997946, 1.998162, 1.998320, 1.998364, 1.998611, 1.998815, 1.998948, 1.999125, 1.999207, 1.999222, 1.999229, 1.999235, 1.999241, 1.999242, 1.999247, 1.999248, 1.999250, 1.999251, <

				clamp4(cL+cT-cTL, cL, cTL, cT, cTR),//0
				//clip(cL+cT-cTL),//1
				clamp4(cL+cTR-cT, cL, cTL, cT, cTR),//2
				//clip(cL+cTR-cT),//3
				clamp4(cT+cTL-cT2L, cL, cTL, cT, cTR),//4
				//clip(cT+cTL-cT2L),//5
				clamp4(cT+cTR-cT2R, cL, cT, cTR, cTR2),//6
				//clip(cT+cTR-cT2R),//7
				(cL+cTR2)>>1,//8
				//clip(cT*3-cT2*3+cT3),//9
				//clip(cL*3-cL2*3+cL3),//10
				//(cL+clip(cTR*3-cT2R*3+cT3R))>>1,//11
				//(cL+clip(cTR2*3-cT2R3*3+cT3R4))>>1,//12
				//clip(cT2+cT4-cT6),//13
				//clip(cL2+cL4-cL6),//14
				//clip((cT5-6*cT4+15*cT3-20*cT2+15*cT+clamp4(cL*2-cTL2, cL, cTL, cT, cT2))/6),//15
				//clip((-3*cL2+8*cL+clamp4(3*cTR2-3*cT2R2+cT3R2, cTR, cTR2, cTR3, cTR4))/6),//16
				//clip(cT2+cTL-cT3L),//17
				//clip(cT2+cTR-cT3R),//18
				//clip((cL*2+cTL) - (cL2*2+cTL2) + cL3),//19
				//clip(3*(cTL+cTL2)/2-cT2L3*3+(cT3L4+cT3L5)/2),//20
				//clip(cTR2+cTR-cT2R3),//21
				//clip(cTL2+cL2-cL4),//22
				//clip(((cL+cTL)*3-cTL2*6+cTL3+cT2L3)/2),//23
				//clip((cTR*2+cTR2) - (cT2R2+cT3R2*2) + cT4R3),//24
				cT6,//25
				(cTR4+cTR6)>>1,//26
				(cL4+cL6)>>1,//27
				(cL+cT+cTR5+cTR7)>>2,//28
				//clip(cTR3+cL-cTR2),//29
				//clip(4*cT3-3*cT4),//30
				//clip(cT+cT2-cT3),//31
				//clip(cL+cL2-cL3),//32
				//clip(cL+cTR2-cTR),//33
				//clip(cL2+cTR2-cT),//34
				//(clip(cL*2-cTL)+clip(cL*2-cTL2)+cT+cTR)>>2,//35
				clamp4(cT*2-cT2, cL, cT, cTR, cTR2),//36
				(cT+cT3)>>1,//37
				//clip(cT2+cL-cT2L),//38
				//clip(cTR2+cT-cT2R2),//39
				//clip((4*cL3-15*cL2+20*cL+clip(cTR2*2-cT2R2))/10),//40
				//clip((cT3R3-4*cT2R2+6*cTR+clip(cL*3-cTL*3+cT2L))/4),//41
				//clip((cT*2+cTR) - (cT2+2*cT2R) + cT3R),//42
				//clip((cTL*2+cT2L) - (cT2L2+cT3L2*2) + cT4L3),//43
				//clip(cT2L2+cL-cT2L3),//44
				//clip((-cT4+5*cT3-10*cT2+10*cT+clip(cL*4-cTL2*6+cT2L3*4-cT3L4))/5),//45
				//clip(cTR2+clip(cTR3*2-cT2R4-cTR4)),//46
				//clip(cTL+cL-cTL2),//47
				//clip((cT*2+cTL) - (cT2+2*cT2L) + cT3L),//48
				//clip(cT2+clip(cTR2*2-cT2R3) - cT2R),//49
				//clip((-cL4+5*cL3-10*cL2+10*cL+clip(cTR*2-cT2R))/5),//50
				//clip((-cL5+4*cL4-5*cL3+5*cL+clip(cTR*2-cT2R))>>2),//51
				//clip((cL3-4*cL2+6*cL+clip(cTR*3-cT2R*3+cT3R))>>2),//52
				//clip((-cT2R2+3*cTR+clip(4*cL-6*cTL+4*cT2L-cT3L))/3),//53
				((cL+cT)*3-cTL*2)>>2,//54
#endif
#if 0

				pred_w2_paeth2(cT, cL, cTL, cTR),
				(cL<<1)-cL2 + eL*params[PW2_NPRED+7],
				(cT<<1)-cT2,
				(cTL<<1)-cT2L2,
				(cTR<<1)-cT2R2,

				0,
				cL,
				cT,
				pred_w2_select(cT, cL, cTL),
				pred_w2_cgrad(cT, cL, cTL),
				cTL,
				cTR,
				cL2,
				(cL+cT)>>1,
				(cL+cTL)>>1,
				(cT+cTL)>>1,
				(cT+cTR)>>1,
				(6*cT-2*cT2+7*cL+cL2+cTR2+3*cTR+8)>>4,
#endif
			};

			long long pred, sum=0;
			for(int k=0;k<PW2_NPRED;++k)
				sum+=weights[k];
			if(sum)
			{
				pred=(sum>>1)-1;
				for(int k=0;k<PW2_NPRED;++k)
					pred+=predictions[k]*weights[k];
				pred/=sum;
			}
			else
				pred=predictions[0];

			int vmin=cL, vmax=cL, curr;
			//if(vmin>cTR)
			//	vmin=cTR;
			if(vmin>cT)
				vmin=cT;

			//if(vmax<cTR)
			//	vmax=cTR;
			if(vmax<cT)
				vmax=cT;

			pred=CLAMP(vmin, pred, vmax);
			if(fwd)
			{
				curr=src[idx]<<8;
				dst[idx]=src[idx]-(int)((pred+127)>>8);
			}
			else
			{
				dst[idx]=src[idx]+(int)((pred+127)>>8);
				curr=dst[idx]<<8;
			}

			error[currrow+kx]=curr-(int)pred;
			for(int k=0;k<PW2_NPRED;++k)
			{
				int e=abs(curr-predictions[k]);
				pw2_errors[k]+=e;//
				pred_errors[k][currrow+kx]=e;
				if(kx+1<iw)
					pred_errors[k][prevrow+kx+1]+=e;
			}
		}
	}
	for(int k=0;k<PW2_NPRED;++k)
		pw2_errors[k]/=iw*ih*256;
}
void pred_w2_apply(char *buf, int iw, int ih, short *allparams, int fwd)
{
	int res=iw*ih;
	int *temp=(int*)malloc((size_t)iw*(PW2_NPRED+1)*2*sizeof(int));
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_w2_prealloc(buf, iw, ih, 0, allparams             , fwd, buf2, temp);
	pred_w2_prealloc(buf, iw, ih, 1, allparams+PW2_NPARAM  , fwd, buf2, temp);
	pred_w2_prealloc(buf, iw, ih, 2, allparams+PW2_NPARAM*2, fwd, buf2, temp);

	for(int k=0;k<res;++k)
	{
		buf[k<<2  ]=buf2[k<<2  ];
		buf[k<<2|1]=buf2[k<<2|1];
		buf[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}
double pred_w2_calcloss(const char *src, int iw, int ih, int kc, short *params, int *temp, char *dst, int *hist)
{
	int res=iw*ih;
	pred_w2_prealloc(src, iw, ih, kc, params, 1, dst, temp);
	//addhalf((unsigned char*)dst+kc, iw, ih, 1, 4);
	calc_histogram(dst+kc, (ptrdiff_t)res<<2, 4, hist);

	double entropy=0;
	for(int k=0;k<256;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/res;
			entropy-=p*log2(p);
		}
	}
	double invCR=entropy/8, csize=res*invCR;
	return csize;
}
void pred_w2_opt_v2(char *buf2, int iw, int ih, short *params, int loud)
{
	int res=iw*ih;
	char *buf3=(char*)malloc((size_t)res<<2);
	int *temp=(int*)malloc((size_t)iw*(PW2_NPRED+1)*2*sizeof(int));
	int *hist=(int*)malloc(256*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	char title0[256];
	get_window_title(title0, 256);
	int steps[]={128, 64, 32, 16, 8, 4, 2, 1};
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*PW2_NPARAM;
		double csize0=pred_w2_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
		for(int ks=0;ks<COUNTOF(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<PW2_NPARAM;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=step;
				csize=pred_w2_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_w2_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;

				set_window_title("Ch%d csize %lf [%d/%d %d/%d]...", kc, csize0, kc*COUNTOF(steps)+ks+1, COUNTOF(steps)*3, idx+1, PW2_NPARAM);//
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
			//set_window_title("Ch%d csize %lf [%d/%d]...", kc, csize0, kc*COUNTOF(steps)+ks+1, COUNTOF(steps)*3);//
		}
	}
	free(hist);
	free(temp);
	free(buf3);
	set_window_title("%s", title0);//
}
void pred_jmj_apply(char *buf, int iw, int ih, int fwd)
{
	int res=iw*ih;
	int *temp=(int*)malloc((size_t)iw*(PW2_NPRED+1)*2*sizeof(int));
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_jxl_prealloc(buf, iw, ih, 0, jxlparams_i16        , fwd, buf2, temp);
	pred_w2_prealloc (buf, iw, ih, 1, pw2_params+PW2_NPARAM, fwd, buf2, temp);
	pred_jxl_prealloc(buf, iw, ih, 2, jxlparams_i16+22     , fwd, buf2, temp);

	for(int k=0;k<res;++k)
	{
		buf[k<<2  ]=buf2[k<<2  ];
		buf[k<<2|1]=buf2[k<<2|1];
		buf[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}

short jxlparams_i16[33]=//signed fixed 7.8 bit
{
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,

	//0x0DD3, 0x1002, 0x11F9, 0x0BEC, 0x0002,-0x0029,-0x010E,-0x00E8,-0x00A2,-0x0089, 0x0064,
	//0x0FAF, 0x0C9E, 0x139F, 0x0E13,-0x009D,-0x005C, 0x0041, 0x0013, 0x00DD,-0x0004,-0x00C3,
	//0x09A8, 0x0F00, 0x1191, 0x0B9E,-0x0026,-0x0027,-0x0086,-0x0084,-0x005F,-0x0038, 0x00B2,

	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,  -0x0148,-0x0029, 0, 0, 0, 0, 0,
	//0x0A14, 0x0829, 0x1548, 0x0CA4,   0x047B, 0x0052, 0, 0, 0, 0, 0,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,  -0x0148,-0x00F6, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,-0x0148,-0x0029, 0x0971, 0x01EC,-0x01C3, 0x047B, 0x0AB8,
	//0x0A14, 0x0829, 0x1548, 0x0CA4, 0x047B, 0x0052,-0x011F, 0x0000, 0x0029, 0x063D, 0x0266,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,-0x0148,-0x00F6, 0x0614, 0x00A4,-0x007B, 0x019A, 0x0E8F,

	//0x63, 0x5A, 0x50, 0x59,    -0x0A, -0x01,  0x4B, 0x0F, -0x0E, -0x24, 0x55,
	//0x50, 0x41, 0xA9, 0x64,     0x24,  0x03, -0x09,    0,  0x01,  0x32, 0x13,
	//0x59, 0x61, 0x6D, 0x8C,    -0x0A, -0x08,  0x30, 0x05, -0x04,  0x0D, 0x74,
};
//float jxlparams_ps[33]=
//{
//	0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
//	0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
//	0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
//};
void pred_jxl_prealloc(const char *src, int iw, int ih, int kc, short *params, int fwd, char *dst, int *temp_w10)
{
#if 0
	params[0]=0;
	params[1]=0;
	params[2]=0;
	params[3]=0x100;
#endif
	int errorbuflen=iw<<1, rowlen=iw<<2;
	int *error=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int idx=kc;
	const char *src2=fwd?src:dst;
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred, curr;
			
			char
				ctt      =         ky-2>=0?src2[idx-rowlen*2]:0,
				ctopleft =kx-1>=0&&ky-1>=0?src2[idx-rowlen-4]:0,
				ctop     =kx  <iw&&ky-1>=0?src2[idx-rowlen  ]:0,
				ctopright=kx+1<iw&&ky-1>=0?src2[idx-rowlen+4]:0,
				cleft    =kx-1>=0         ?src2[idx       -4]:0;

			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	kx=iw>>1;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c
			
			int weights[4];//fixed 23.8 bit
			for(int k=0;k<4;++k)
			{
				int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
				weights[k]=(params[k]<<8)/(w+1);
			}

			int
				etop=ky-1>=0?error[prevrow+kx]:0,
				eleft=kx-1>=0?error[currrow+kx-1]:0,
				etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
				etopplusleft=etop+eleft;
			int predictions[]=//fixed 23.8 bit
			{
				(cleft+ctopright-ctop)<<8,
				(ctop<<8)-((etopplusleft+etopright)*params[4]>>8),
				(cleft<<8)-((etopplusleft+etopleft)*params[5]>>8),
				(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
			};

			int sum=weights[0]+weights[1]+weights[2]+weights[3];
			if(sum)
				pred=(predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum;
			else
				pred=predictions[0];

			int vmin=cleft, vmax=cleft;
			if(vmin>ctopright)
				vmin=ctopright;
			if(vmin>ctop)
				vmin=ctop;

			if(vmax<ctopright)
				vmax=ctopright;
			if(vmax<ctop)
				vmax=ctop;

			vmin<<=8;
			vmax<<=8;

			pred=CLAMP(vmin, pred, vmax);
			if(fwd)
			{
				curr=src[idx]<<8;
				dst[idx]=src[idx]-((pred+127)>>8);
			}
			else
			{
				dst[idx]=src[idx]+((pred+127)>>8);
				curr=dst[idx]<<8;
			}

			error[currrow+kx]=curr-pred;
			for(int k=0;k<4;++k)
			{
				int e=abs(curr-predictions[k]);
				pred_errors[k][currrow+kx]=e;
				if(kx+1<iw)
					pred_errors[k][prevrow+kx+1]+=e;
			}
		}
	}
}
void calc_histogram(unsigned char *buf, ptrdiff_t len, ptrdiff_t stride, int *hist)
{
	memset(hist, 0, 256*sizeof(int));
	for(ptrdiff_t k=0, end=len-(stride-1);k<end;k+=stride)
		++hist[buf[k]];
}
double pred_jxl_calcloss(const char *src, int iw, int ih, int kc, short *params, int *temp, char *dst, int *hist)
{
	int res=iw*ih;
	pred_jxl_prealloc(src, iw, ih, kc, params, 1, dst, temp);
	//addhalf((unsigned char*)dst+kc, iw, ih, 1, 4);
	calc_histogram(dst+kc, (ptrdiff_t)res<<2, 4, hist);

	double entropy=0;
	for(int k=0;k<256;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/res;
			p*=0x10000-256;
			++p;
			p/=0x10000;
			entropy-=p*log2(p);
		}
	}
	double invCR=entropy/8, csize=res*invCR;
	return csize;
}
void pred_jxl_opt_v2(char *buf2, int iw, int ih, short *params, int loud)
{
	int res=iw*ih;
	char *buf3=(char*)malloc((size_t)res<<2);
	int *temp=(int*)malloc((size_t)iw*10*sizeof(int));
	int *hist=(int*)malloc(256*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int steps[]={256, 128, 64, 32, 16, 8, 4, 2, 1};
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*11;
		double csize0=pred_jxl_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
		for(int ks=0;ks<COUNTOF(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<11;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=step;
				csize=pred_jxl_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_jxl_calcloss(buf2, iw, ih, kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
			//set_window_title("%lf", csize0);//
		}
	}
	free(hist);
	free(temp);
	free(buf3);
}
#if 0
void pred_jxl_optimize(const char *src, int iw, int ih, int kc, short *params, int step, int pidx, char *dst, int loud)
{
	int *temp=(int*)malloc((size_t)iw*10*sizeof(int));
	int *hist=(int*)malloc(256*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	int res=iw*ih;
	double csize0, csize, csize00;
	short p0=params[pidx];
	csize=pred_jxl_calcloss(src, iw, ih, kc, params, temp, dst, hist);
	csize00=csize;

	//static int it=0;
	//for(int k=it, end=it+niter;k<end;++k, ++it)
	//{
		//int idx=it%11;
	int subit;
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		params[pidx]+=step;
		csize=pred_jxl_calcloss(src, iw, ih, kc, params, temp, dst, hist);
		if(csize>csize0)//cancel last change and break
		{
			params[pidx]-=step;
			csize=csize0;
			break;
		}
	}
		
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		params[pidx]-=step;
		csize=pred_jxl_calcloss(src, iw, ih, kc, params, temp, dst, hist);
		if(csize>csize0)
		{
			params[pidx]+=step;
			csize=csize0;
			break;
		}
	}
	if(csize>csize00)//prevent CR from worsening
	{
		params[pidx]=p0;
		csize=csize00;
	}
	//}
	//it%=33;
	//if(loud)
	//	printf("\n");
	
	free(hist);
	free(temp);
}
#endif
void pred_jxl_apply(char *buf, int iw, int ih, short *allparams, int fwd)
{
	int res=iw*ih;
	int *temp=(int*)malloc((size_t)iw*10*sizeof(int));
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_jxl_prealloc(buf, iw, ih, 0, allparams   , fwd, buf2, temp);
	pred_jxl_prealloc(buf, iw, ih, 1, allparams+11, fwd, buf2, temp);
	pred_jxl_prealloc(buf, iw, ih, 2, allparams+22, fwd, buf2, temp);

	for(int k=0;k<res;++k)
	{
		buf[k<<2  ]=buf2[k<<2  ];
		buf[k<<2|1]=buf2[k<<2|1];
		buf[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}

#if 0
void pred_jxl(char *buf, int iw, int ih, int nch, int bytestride, int fwd)
{
	if(!(customparam_st[0]+customparam_st[1]+customparam_st[2]+customparam_st[3]))//reset params if all weights are zero
	{
		customparam_st[0]=1, customparam_st[1]=0.5, customparam_st[2]=1.1, customparam_st[3]=1.1;
		customparam_st[5]=0.5, customparam_st[6]=customparam_st[7]=0, customparam_st[8]=0.5, customparam_st[9]=0.15;
		customparam_st[10]=0.2, customparam_st[11]=0.002;
	}
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res*bytestride);
	int errorbuflen=iw*2;
	char *pred_errors[]=
	{
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
	};
	char *error=(char*)malloc((size_t)iw*2);
	if(!b2||!pred_errors[0]||!pred_errors[1]||!pred_errors[2]||!pred_errors[3]||!error)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(pred_errors[0], 0, errorbuflen);
	memset(pred_errors[1], 0, errorbuflen);
	memset(pred_errors[2], 0, errorbuflen);
	memset(pred_errors[3], 0, errorbuflen);
	memset(error, 0, errorbuflen);
	memset(b2, 0, (size_t)res*bytestride);
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred, curr;
				
				char *src=fwd?buf:b2;
				char
					ctl3     =kx-3>=0&&ky-3>=0?src[idx-rowlen*3-bytestride*3]:0,
					ct3      =         ky-3>=0?src[idx-rowlen*3             ]:0,
					ctr3     =kx+3<iw&&ky-3>=0?src[idx-rowlen*3+bytestride*3]:0,

					ctltl    =kx-2>=0&&ky-2>=0?src[idx-rowlen*2-bytestride*2]:0,
					ctt      =         ky-2>=0?src[idx-rowlen*2             ]:0,
					ctrtr    =kx+2<iw&&ky-2>=0?src[idx-rowlen*2+bytestride*2]:0,

					ctopleft =kx-1>=0&&ky-1>=0?src[idx-rowlen-bytestride  ]:0,
					ctop     =kx  <iw&&ky-1>=0?src[idx-rowlen             ]:0,
					ctopright=kx+1<iw&&ky-1>=0?src[idx-rowlen+bytestride  ]:0,
					ctrr     =kx+2<iw&&ky-1>=0?src[idx-rowlen+bytestride*2]:0,
		
					cl3      =kx-3>=0?src[idx-bytestride*3]:0,
					cll      =kx-2>=0?src[idx-bytestride*2]:0,
					cleft    =kx-1>=0?src[idx-bytestride  ]:0;


				//w0   w1   w2   w3
				//p3Ca p3Cb p3Cc p3Cd p3Ce
				//p1C  p2c

				double weights[4];
				for(int k=0;k<4;++k)
				{
					int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
					weights[k]=customparam_st[k]/(w+1);
				}

				char
					etop=ky-1>=0?error[prevrow+kx]:0,
					eleft=kx-1>=0?error[currrow+kx-1]:0,
					etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
					etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0;
				double predictions[]=
				{
					cleft+ctopright-ctop,
					//ctop+cleft-ctopleft,
					ctop-(int)((etop+eleft+etopright)*customparam_st[10]),
					cleft-(int)((etop+eleft+etopleft)*customparam_st[11]),
					ctop-(int)(etopleft*customparam_st[5]+ctop*customparam_st[6]+ctopright*customparam_st[7]+(ctt-ctop)*customparam_st[8]+(ctopleft-cleft)*customparam_st[9]),
				};

				double wsum=weights[0]+weights[1]+weights[2]+weights[3];
				if(wsum)
					pred=(int)round((predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3])/wsum);
				else
					pred=(int)round(predictions[0]);

				//if(((etop<0)!=(eleft<0))||((etop<0)!=(etopleft<0)))
				//if(((etop^eleft)|(etop^etopleft))<=0)
				{
					int vmin=cleft, vmax=cleft;
					if(vmin>ctopright)
						vmin=ctopright;
					if(vmin>ctop)
						vmin=ctop;

					if(vmax<ctopright)
						vmax=ctopright;
					if(vmax<ctop)
						vmax=ctop;

					pred=CLAMP(vmin, pred, vmax);
				}
				//pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				if(fwd)
				{
					curr=buf[idx];
					b2[idx]=buf[idx]-pred;
				}
				else
				{
					b2[idx]=buf[idx]+pred;
					curr=b2[idx];
				}

				error[currrow+kx]=pred-curr;
				for(int k=0;k<4;++k)
				{
					int e=(int)round(fabs(curr-predictions[k]));
					pred_errors[k][currrow+kx]=e;
					if(kx+1<iw)
						pred_errors[k][prevrow+kx+1]+=e;
				}
			}
		}
	}
	for(int kc=nch;kc<bytestride;++kc)
	{
		for(int k=0;k<res;++k)
			b2[k*bytestride+kc]=buf[k*bytestride+kc];
	}
	memcpy(buf, b2, (size_t)res*bytestride);
	free(b2);
	free(pred_errors[0]);
	free(pred_errors[1]);
	free(pred_errors[2]);
	free(pred_errors[3]);
	free(error);
}
#endif


int sortnb_cases[SORTNBCASES];
double sortnb_rmse[SORTNBCASES];
void pred_sortnb(char *buf, int iw, int ih, int nch, int bytestride, int fwd)
{
	memset(sortnb_cases, 0, sizeof(sortnb_cases));
	memset(sortnb_rmse, 0, sizeof(sortnb_rmse));
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res*bytestride);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(b2, 0, (size_t)res*bytestride);
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			int error[3]={0};
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred;

				int curr=buf[idx];
				char *src=fwd?buf:b2;
				char
					top     =         ky-1>=0?src[idx-rowlen           ]:0,
					left    =kx-1>=0         ?src[idx       -bytestride]:0,
					topleft =kx-1>=0&&ky-1>=0?src[idx-rowlen-bytestride]:0,
					topright=kx+1<iw&&ky-1>=0?src[idx-rowlen+bytestride]:0;
				char nb[]={topleft, top, topright, left};
				int permutation=0;
				char temp;

#define SORT_STEP(A, B)\
				if(nb[A]<nb[B])\
					permutation+=0;\
				else if(nb[A]>nb[B])\
					permutation+=1, temp=nb[A], nb[A]=nb[B], nb[B]=temp;\
				else\
					permutation+=2;

				SORT_STEP(0, 1);
				permutation*=3;
				SORT_STEP(0, 2);
				permutation*=3;
				SORT_STEP(0, 3);

				permutation*=3;
				SORT_STEP(1, 2);
				permutation*=3;
				SORT_STEP(1, 3);

				permutation*=3;
				SORT_STEP(2, 3);
#undef  SORT_STEP

				int c=0;
				switch(permutation)
				{
				case 232:case 253:case 259:case 256:case 283:case 293:case 301:case 310:case 320:case 364:
				case 391:case 415:case 416:case 448:case 475:case 547:case 617:case 644:case 701:case 728:
					pred=nb[0];//A
					c=0;
					break;
				case 50:case 205:case 250:case 254:case 266:case 269:case 274:case 337:case 340:case 520:
					pred=(nb[0]+nb[1])>>1;//(A+B)/2
					c=1;
					break;
				case 26:case 31:case 245:case 334:case 335:case 590:case 593:
					pred=nb[1];//B
					c=2;
					break;
				case 94:case 148:case 244:
					pred=(nb[1]+nb[2])>>1;//(B+C)/2
					c=3;
					break;
				case 4:case 11:case 58:case 77:case 92:case 173:case 333:case 414:case 487:case 488:
					pred=nb[2];//C
					c=4;
					break;
				case 7:case 13:case 16:case 67:case 97:case 171:case 172:case 261:case 585:
					pred=(nb[2]+nb[3])>>1;//(C+D)/2
					c=5;
					break;
				case 0:case 2:case 9:case 10:case 18:case 23:case 486:case 666:
					pred=nb[3];//D
					c=6;
					break;
				case 1:case 40:case 90:case 91:case 121:case 243:case 247:case 252:default:
					pred=pred_grad(top, left, topleft);
					c=7;
					break;
				}

				if(fwd)
					b2[idx]=buf[idx]-pred;
				else
					b2[idx]=buf[idx]+pred;

				++sortnb_cases[c];
				if(fwd)
					sortnb_rmse[c]+=b2[idx]*b2[idx];
				else
					sortnb_rmse[c]+=buf[idx]*buf[idx];
			}
		}
	}
	for(int kc=nch;kc<bytestride;++kc)
	{
		for(int k=0;k<res;++k)
			b2[k*bytestride+kc]=buf[k*bytestride+kc];
	}
	memcpy(buf, b2, (size_t)res*bytestride);
	free(b2);
	for(int k=0;k<SORTNBCASES;++k)
		sortnb_rmse[k]=sqrt(sortnb_rmse[k]/sortnb_cases[k]);
}

int  predict_grad(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char
		left=kx?buf[idx-bytestride]:0,
		top=ky?buf[idx-rowlen]:0,
		topleft=kx&&ky?buf[idx-rowlen-bytestride]:0;

	int pred;

	char vmax, vmin;
	if(top<left)
		vmin=top, vmax=left;
	else
		vmin=left, vmax=top;

	if(topleft>vmax)//choose steepest slope if both downward or upward
		pred=vmin;
	else if(topleft<vmin)
		pred=vmax;
	else
		pred=left+top-topleft;//planar prediction (unplane)

	//char xdelta=top-topleft, ydelta=left-topleft;
	//if((xdelta>0)==(ydelta>0))
	//	pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
	//else
	//	pred=topleft+xdelta+ydelta;//average slope
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_grad_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_grad(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_grad_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_grad(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}

int  predict_custom(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, const double *params, double sum)
{
#if 0
#define LOAD(BUF, X, Y) (unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|kc]:0
	char
		NNNWWW =LOAD(pixels,  2, 3),
		NNWW   =LOAD(pixels,  2, 2),//1
		NNW    =LOAD(pixels,  1, 2),//4
		NN     =LOAD(pixels,  0, 2),//3
		NNE    =LOAD(pixels, -1, 2),//3
		NNEE   =LOAD(pixels, -2, 2),//2
		NWW    =LOAD(pixels,  2, 1),//2
		NW     =LOAD(pixels,  1, 1),//7
		N      =LOAD(pixels,  0, 1),//13
		NE     =LOAD(pixels, -1, 1),//8
		NEE    =LOAD(pixels, -2, 1),//2
		WW     =LOAD(pixels,  2, 0),//3
		W      =LOAD(pixels,  1, 0),//11
		eNNNWWW=LOAD(errors,  2, 2),
		eNNWW  =LOAD(errors,  2, 2),
		eNNW   =LOAD(errors,  1, 2),
		eNN    =LOAD(errors,  0, 2),
		eNNE   =LOAD(errors, -1, 2),
		eNNEE  =LOAD(errors, -2, 2),
		eNWW   =LOAD(errors,  2, 1),
		eNW    =LOAD(errors,  1, 1),
		eN     =LOAD(errors,  0, 1),
		eNE    =LOAD(errors, -1, 1),
		eNEE   =LOAD(errors, -2, 1),
		eWW    =LOAD(errors,  2, 0),
		eW     =LOAD(errors,  1, 0);

	double features[]=
	{
		NNWW, NNW, NN, NNE, NNEE, NWW, NW, N, NE, NEE, WW, W,
		eNNWW, eNNW, eNN, eNNE, eNNEE, eNWW, eNW, eN, eNE, eNEE, eWW, eW,
		//NNWW+eNNWW, NNW+eNNW, NN+eNN, NNE+eNNE, NNEE+eNNEE, NWW+eNWW, NW+eNW, N+eN, NE+eNE, NEE+eNEE, WW+eWW, W+eW
	};
#if 0
	double features[]=
	{
		(N+W-NW) + (eN+eW-eNW)*params[23],
		N + eN*params[23],
		W + eW*params[23],
		(N+W)*0.5*params[23] + (eN+eW)*0.5*params[23],
		(N*2-NN) + (eN*2-eNN)*params[23],
		(W*2-WW) + (eW*2-eWW)*params[23],
		(NW*2-NNWW) + (eNW*2-eNNWW)*params[23],
		(NE+NW-NN) + (eNE+eNW-eNN)*params[23],
		(NE*2-NNEE) + (eNE*2-eNNEE)*params[23],
		(W+NE-N) + (eW+eNE-eN)*params[23],
		(N*3-NNW-NE) + (eN*3-eNNW-eNE)*params[23],
		(W*4+NEE*2-(WW+NE*2))/3 + (eW*4+eNEE*2-(eWW+eNE*2))/3*params[23],
		(N*3-NNE-NW) + (eN*3-eNNE-eNW)*params[23],
		(W+NW-NWW) + (eW+eNW-eNWW)*params[23],
		(N+NE-NNE) + (eN+eNE-eNNE)*params[23],
		(N+NW-NNW) + (eN+eNW-eNNW)*params[23],
		(N-(NE*2+WW)+(W*2+NEE)) + (eN-(eNE*2+eWW)+(eW*2+eNEE))*params[23],
		((N+NW+W+NN)/2-NNW) + ((eN+eNW+eW+eNN)/2-eNNW)*params[23],
		((N+NE+W+NNEE)/2-NNE) + ((eN+eNE+eW+eNNEE)/2-eNNE)*params[23],
		((N+W)-(NNW+NWW)/2) + ((eN+eW)-(eNNW+eNWW)/2)*params[23],
		((NNW+W+NW*2)/2-NNWW) + ((eNNW+eW+eNW*2)/2-eNNWW)*params[23],
		((W+NNE+N*2)/2-NN) + ((eW+eNNE+eN*2)/2-eNN)*params[23],
		((N+WW+NW+W)/2-NWW) + ((eN+eWW+eNW+eW)/2-eNWW)*params[23],
	};
#endif
	double fpred=0;
	for(int k=0;k<_countof(features);++k)
		fpred+=features[k]*params[k];
	if(sum)
		fpred/=sum;
	int pred=(int)round(fpred);
#endif

#if 0
	int pred;
	if(!sum)
		pred=N+W-NW;
	else
		pred=(int)((
			params[ 0]*(N+W-NW)+
			params[ 1]*N+
			params[ 2]*W+
			params[ 3]*(N+W)*0.5+
			params[ 4]*(N*2-NN)+
			params[ 5]*(W*2-WW)+
			params[ 6]*(NW*2-NNWW)+
			params[ 7]*(NE+NW-NN)+
			params[ 8]*(NE*2-NNEE)+
			params[ 9]*(W+NE-N)+
			params[10]*(N*3-NNW-NE)+
			params[11]*(W*4+NEE*2-(WW+NE*2))/3+
			params[12]*(eN+eW+eNW)+
			params[13]*eN+
			params[14]*eW+
			params[15]*(eN+eW)*0.5+
			params[16]*(eN*2-eNN)+
			params[17]*(eW*2-eWW)+
			params[18]*(eNW*2-eNNWW)+
			params[19]*(eNE+eNW-eNN)+
			params[20]*(eNE*2-eNNEE)+
			params[21]*(eW+eNE-eN)+
			params[22]*(eN*3-eNNW-eNE)+
			params[23]*(eW*4+eNEE*2-(eWW+eNE*2))/3
		)/sum);
#endif

#if 1
#define LOAD(BUF, X, Y) (unsigned)(kx-(X))<(unsigned)iw&&(unsigned)(ky-Y)<(unsigned)ih?BUF[(iw*(ky-Y)+kx-(X))<<2|kc]:0
	char comp[]=
	{
#if 0
		LOAD(pixels,  3, 3),
		LOAD(pixels,  2, 3),
		LOAD(pixels,  1, 3),
		LOAD(pixels,  0, 3),
		LOAD(pixels, -1, 3),
		LOAD(pixels, -2, 3),
		LOAD(pixels, -3, 3),
		
		LOAD(pixels,  3, 2),
		LOAD(pixels,  2, 2),
		LOAD(pixels,  1, 2),
		LOAD(pixels,  0, 2),
		LOAD(pixels, -1, 2),
		LOAD(pixels, -2, 2),
		LOAD(pixels, -3, 2),
		
		LOAD(pixels,  3, 1),
		LOAD(pixels,  2, 1),
		LOAD(pixels,  1, 1),
		LOAD(pixels,  0, 1),
		LOAD(pixels, -1, 1),
		LOAD(pixels, -2, 1),
		LOAD(pixels, -3, 1),
		
		LOAD(pixels,  3, 0),
		LOAD(pixels,  2, 0),
		LOAD(pixels,  1, 0),
#endif
#if 1
		LOAD(pixels,  2, 2),
		LOAD(pixels,  1, 2),
		LOAD(pixels,  0, 2),
		LOAD(pixels, -1, 2),
		LOAD(pixels, -2, 2),
		
		LOAD(pixels,  2, 1),
		LOAD(pixels,  1, 1),
		LOAD(pixels,  0, 1),
		LOAD(pixels, -1, 1),
		LOAD(pixels, -2, 1),
		
		LOAD(pixels,  2, 0),
		LOAD(pixels,  1, 0),
		
		LOAD(errors,  2, 2),
		LOAD(errors,  1, 2),
		LOAD(errors,  0, 2),
		LOAD(errors, -1, 2),
		LOAD(errors, -2, 2),
		
		LOAD(errors,  2, 1),
		LOAD(errors,  1, 1),
		LOAD(errors,  0, 1),
		LOAD(errors, -1, 1),
		LOAD(errors, -2, 1),
		
		LOAD(errors,  2, 0),
		LOAD(errors,  1, 0),
#endif
	};
#if 0
	char comp[]=
	{
		kx-2>=0&&ky-2>=0?pixels[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?pixels[idx-(rowlen<<1)-bytestride]:0,
		         ky-2>=0?pixels[idx-(rowlen<<1)]:0,
		kx+1<iw&&ky-2>=0?pixels[idx-(rowlen<<1)+bytestride]:0,
		kx+2<iw&&ky-2>=0?pixels[idx-(rowlen<<1)+(bytestride<<1)]:0,

		kx-2>=0&&ky-1>=0?pixels[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?pixels[idx-rowlen-bytestride]:0,
		         ky-1>=0?pixels[idx-rowlen]:0,
		kx+1<iw&&ky-1>=0?pixels[idx-rowlen+bytestride]:0,
		kx+2<iw&&ky-1>=0?pixels[idx-rowlen+(bytestride<<1)]:0,

		kx-2>=0?pixels[idx-(bytestride<<1)]:0,
		kx-1>=0?pixels[idx-bytestride]:0,

		kx-2>=0&&ky-2>=0?errors[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?errors[idx-(rowlen<<1)-bytestride]:0,
		         ky-2>=0?errors[idx-(rowlen<<1)]:0,
		kx+1<iw&&ky-2>=0?errors[idx-(rowlen<<1)+bytestride]:0,
		kx+2<iw&&ky-2>=0?errors[idx-(rowlen<<1)+(bytestride<<1)]:0,

		kx-2>=0&&ky-1>=0?errors[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?errors[idx-rowlen-bytestride]:0,
		         ky-1>=0?errors[idx-rowlen]:0,
		kx+1<iw&&ky-1>=0?errors[idx-rowlen+bytestride]:0,
		kx+2<iw&&ky-1>=0?errors[idx-rowlen+(bytestride<<1)]:0,

		kx-2>=0?errors[idx-(bytestride<<1)]:0,
		kx-1>=0?errors[idx-bytestride]:0,
	};
#endif

	int pred=(int)(
		params[ 0]*comp[ 0]+
		params[ 1]*comp[ 1]+
		params[ 2]*comp[ 2]+
		params[ 3]*comp[ 3]+
		params[ 4]*comp[ 4]+
		params[ 5]*comp[ 5]+
		params[ 6]*comp[ 6]+
		params[ 7]*comp[ 7]+
		params[ 8]*comp[ 8]+
		params[ 9]*comp[ 9]+
		params[10]*comp[10]+
		params[11]*comp[11]+
		params[12]*comp[12]+
		params[13]*comp[13]+
		params[14]*comp[14]+
		params[15]*comp[15]+
		params[16]*comp[16]+
		params[17]*comp[17]+
		params[18]*comp[18]+
		params[19]*comp[19]+
		params[20]*comp[20]+
		params[21]*comp[21]+
		params[22]*comp[22]+
		params[23]*comp[23]
	);
#endif
	//int vmin, vmax;
	//if(comp[7]<comp[11])
	//	vmin=comp[7], vmax=comp[11];
	//else
	//	vmin=comp[11], vmax=comp[7];
	//pred=CLAMP(vmin, pred, vmax);
	pred=CLAMP(-128, pred, 127);
	return pred;
}
void pred_custom_prealloc(const char *src, int iw, int ih, int kc, int fwd, const double *ch_params, char *dst)
{
	int idx=0, rowlen=iw*4;
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	double sum=0;
	for(int k=0;k<23;++k)
		sum+=ch_params[k];
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			idx=(iw*ky+kx)<<2|kc;
			int pred=predict_custom(pixels, errors, iw, ih, kc, kx, ky, ch_params, sum);

			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
		}
	}
}
void pred_custom_apply(char *src, int iw, int ih, int fwd, const double *allparams)
{
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	pred_custom_prealloc(src, iw, ih, 0, fwd, allparams, dst);
	pred_custom_prealloc(src, iw, ih, 1, fwd, allparams+24, dst);
	pred_custom_prealloc(src, iw, ih, 2, fwd, allparams+24*2, dst);
	for(int k=0;k<res;++k)//set alpha
		dst[k<<2|3]=0xFF;
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
}
#if 0
void pred_custom_fwd(char *buf, int iw, int ih, int nch, int bytestride, const double *params)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		const double *p2=params+12*kc;
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char pred=predict_custom(buf, iw, kx, ky, idx, bytestride, rowlen, p2);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_custom_inv(char *buf, int iw, int ih, int nch, int bytestride, const double *params)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		const double *p2=params+12*kc;
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_custom(buf, iw, kx, ky, idx, bytestride, rowlen, p2);

				buf[idx]+=pred;
			}
		}
	}
}
#endif


//DWTs
void dwt2d_grad_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	//for(int it=sizes_start;it<1;++it)
	{
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				char v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				//if((unsigned char)v[0]==0xFF||(unsigned char)v[1]==0xFF||(unsigned char)v[2]==0xFF||(unsigned char)v[3]==0xFF)
				//	v[0]=0xFF;

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]-=vmax;
				else if(v[0]>vmax)
					v[3]-=vmin;
				else
					v[3]-=v[1]+v[2]-v[0];

				v[2]-=v[0];
				v[1]-=v[0];
				v[0]+=(v[1]+v[2])>>2;

				buffer[ idx            *stride]=v[3];//grad
				buffer[(idx         +1)*stride]=v[2];//diffy
				buffer[(idx+sizes->w  )*stride]=v[1];//diffx
				buffer[(idx+sizes->w+1)*stride]=v[0];//av

				//buffer[ idx            <<2|kc]=v[0];//av
				//buffer[(idx         +1)<<2|kc]=v[1];//diffx
				//buffer[(idx+sizes->w  )<<2|kc]=v[2];//diffy
				//buffer[(idx+sizes->w+1)<<2|kc]=v[3];//grad
			}
		}
		dwt2d_lazy_fwd(buffer, sizes, it, it+2, stride, temp);
	}
}
void dwt2d_grad_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		dwt2d_lazy_inv(buffer, sizes, it, it+2, stride, temp);
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				char v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				v[0]-=(v[1]+v[2])>>2;
					
				v[1]+=v[0];
				v[2]+=v[0];

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]+=vmax;
				else if(v[0]>vmax)
					v[3]+=vmin;
				else
					v[3]+=v[1]+v[2]-v[0];

				buffer[ idx            *stride]=v[3];
				buffer[(idx         +1)*stride]=v[2];
				buffer[(idx+sizes->w  )*stride]=v[1];
				buffer[(idx+sizes->w+1)*stride]=v[0];
			}
		}
	}
}

double customdwtparams[12]={0};
void dwt1d_custom_fwd(char *buffer, int count, int stride, char *b2, const double *params)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	even[0]-=(char)floor(odd[0]*(params[0]+params[5]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[0]+odd[k]*params[5]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[0]+params[5]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*params[1]+even[k+1]*params[6]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(params[1]+params[6]));


	even[0]-=(char)(odd[0]*(params[2]+params[7]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[2]+odd[k]*params[7]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[2]+params[7]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*params[3]+even[k+1]*params[8]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(params[3]+params[8]));


	even[0]-=(char)floor(odd[0]*(params[4]+params[9]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[4]+odd[k]*params[9]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[4]+params[9]));


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_custom_inv(char *buffer, int count, int stride, char *b2, const double *params)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	even[0]+=(char)floor(odd[0]*(params[4]+params[9]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[4]+odd[k]*params[9]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[4]+params[9]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*params[3]+even[k+1]*params[8]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(params[3]+params[8]));
	
	even[0]+=(char)floor(odd[0]*(params[2]+params[7]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[2]+odd[k]*params[7]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[2]+params[7]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*params[1]+even[k+1]*params[6]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(params[1]+params[6]));
	
	even[0]+=(char)floor(odd[0]*(params[0]+params[5]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[0]+odd[k]*params[5]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[0]+params[5]));


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_custom_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_custom_fwd(buffer+rowlen*ky, w2, stride, temp, params);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_custom_fwd(buffer+stride*kx, h2, rowlen, temp, params);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_custom_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_custom_inv(buffer+stride*kx, h2, rowlen, temp, params);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_custom_inv(buffer+rowlen*ky, w2, stride, temp, params);
	}
}

void dwt1d_exp_fwd(char *buffer, int count, int stride, char *b2, const double *paramsx12)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd+extraeven;++k)//predict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			paramsx12[0]*(prev+next)+
			paramsx12[1]*(prev2+next2)+
			paramsx12[2]*(prev3+next3)+
			paramsx12[3]*(prev4+next4)+
			paramsx12[4]*(prev5+next5)
		));
		even[k]-=pred;
		//even[k]-=(9*(prev+next)+prev2+next2)>>4;
	}
	for(int k=0;k<nodd;++k)//update
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			paramsx12[5]*(prev+next)+
			paramsx12[6]*(prev2+next2)+
			paramsx12[7]*(prev3+next3)+
			paramsx12[8]*(prev4+next4)+
			paramsx12[9]*(prev5+next5)
		));
		odd[k]+=update;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_exp_inv(char *buffer, int count, int stride, char *b2, const double *paramsx12)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			paramsx12[5]*(prev+next)+
			paramsx12[6]*(prev2+next2)+
			paramsx12[7]*(prev3+next3)+
			paramsx12[8]*(prev4+next4)+
			paramsx12[9]*(prev5+next5)
		));
		odd[k]-=update;
	}
	for(int k=0;k<nodd+extraeven;++k)//unpredict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			paramsx12[0]*(prev+next)+
			paramsx12[1]*(prev2+next2)+
			paramsx12[2]*(prev3+next3)+
			paramsx12[3]*(prev4+next4)+
			paramsx12[4]*(prev5+next5)
		));
		even[k]+=pred;
		//even[k]+=(9*(prev+next)+prev2+next2)>>4;
	}


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_exp_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_exp_fwd(buffer+rowlen*ky, w2, stride, temp, params);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_exp_fwd(buffer+stride*kx, h2, rowlen, temp, params);
	}
}
void dwt2d_exp_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_exp_inv(buffer+stride*kx, h2, rowlen, temp, params);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_exp_inv(buffer+rowlen*ky, w2, stride, temp, params);
	}
}

ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	ArrayHandle sizes;
	int lw=floor_log2(iw), lh=floor_log2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	int nstages=0;
	DWTSize *p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=iw;
	p->h=ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=w2;
		p->h=h2;
	}
	return sizes;
}

//lifting-based 8bit lazy DWT
void dwt1d_lazy_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_lazy_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_lazy_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit Haar DWT
void dwt1d_haar_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd;++k)//const predictor (difference)
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd;++k)//update (O+(E-O)/2 = average)
		odd[k]+=even[k]>>1;


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_haar_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
		odd[k]-=even[k]>>1;
	
	for(int k=0;k<nodd;++k)//unpredict
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_haar_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit squeeze DWT
char smoothtendency(char B, char a, char n)
{
	char diff=0;
	if(B>=a&&a>=n)
	{
		diff=(4*B-3*n-a+6)/12;
		//      2C = a<<1 + diff - diff&1 <= 2B  so diff - diff&1 <= 2B - 2a
		//      2D = a<<1 - diff - diff&1 >= 2n  so diff + diff&1 <= 2a - 2n
		if(diff-(diff&1)>2*(B-a))
			diff=2*(B-a)+1;
		if(diff+(diff&1)>2*(a-n))
			diff=2*(a-n);
	}
	else if(B<=a&&a<=n)
	{
		diff=(4*B-3*n-a-6)/12;
		//      2C = a<<1 + diff + diff&1 >= 2B  so diff + diff&1 >= 2B - 2a
		//      2D = a<<1 - diff + diff&1 <= 2n  so diff - diff&1 >= 2a - 2n
		if(diff+(diff&1)<2*(B-a))
			diff=2*(B-a)-1;
		if(diff-(diff&1)<2*(a-n))
			diff=2*(a-n);
	}
	return diff;
}
void dwt1d_squeeze_fwd(char *buffer, int count, int stride, char *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		char
			o1=kd>0?buffer[ks-stride]:0,
			e =buffer[ks],
			o =(kd<<1)+1<count?buffer[ks+stride]:0,
			e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:0,//n-1-(idx-(n-1)) = ((n-1)<<1)-idx
			o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:0;

		e-=o;		//diff
		o+=e>>1;	//av
		e2-=o2;
		o2+=e2>>1;
		e-=smoothtendency(o1, o, o2);
		
		if(kd<nodd)
			odd[kd]=o;
		even[kd]=e;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_squeeze_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(int ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
	{
		char
			o1=ks>0?buffer[kd-stride]:0,
			e =even[ks],
			o=ks<nodd?odd[ks]:0,
			o2=ks+1<nodd?odd[ks+1]:0;

		e+=smoothtendency(o1, o, o2);
		o-=e>>1;
		e+=o;

		buffer[kd]=e;
		if(ks<nodd)
			buffer[kd+stride]=o;
	}
}
void dwt2d_squeeze_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_squeeze_fwd(buffer+stride*kx, h2, rowlen, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dB.PNG", it);
	}
}
void dwt2d_squeeze_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_squeeze_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 5/3 DWT
void dwt1d_cdf53_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

#if 0
	if(fabs(customparam_st[0])<1.5)//linear
	{
		even[0]-=odd[0];
		for(int k=1;k<nodd;++k)//linear predictor (deviation from linearity)
			even[k]-=(odd[k-1]+odd[k]+1)>>1;
		if(extraeven)
			even[nodd]-=odd[nodd-1];
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<3.5)//cubic
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev3=k-2>=0?odd[k-2]:0, prev1=k-1>=0?odd[k-1]:0, next1=k<nodd?odd[k]:0, next3=k+1<nodd?odd[k+1]:0;
			even[k]-=(-(prev3+next3)+(prev1+next1)*9+8)>>4;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;

		//for(int k=0;k<nodd;++k)
		//{
		//	int prev3=k-1>=0?even[k-1]:0, prev1=even[k], next1=k+1<nodd+extraeven?even[k+1]:0, next3=k+2<nodd+extraeven?even[k+2]:0;
		//	odd[k]+=((prev3+next3)+(prev1+next1)*3+4)>>3;
		//}
	}
	else if(fabs(customparam_st[0])<5.5)//power 5
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0;
			even[k]-=(3*(prev5+next5)-25*(prev3+next3)+150*(prev1+next1)+128)>>8;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<7.5)//power 7
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev7=k-4>=0?odd[k-4]:0,
				prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0,
				next7=k+3<nodd?odd[k+3]:0;
			even[k]-=(-5*(prev7+next7)+49*(prev5+next5)-245*(prev3+next3)+1225*(prev1+next1)+1024)>>11;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<9.5)//power 9
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev9=k-5>=0?odd[k-5]:0,
				prev7=k-4>=0?odd[k-4]:0,
				prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0,
				next7=k+3<nodd?odd[k+3]:0,
				next9=k+4<nodd?odd[k+4]:0;
			even[k]-=(35*(prev9+next9)-405*(prev7+next7)+2268*(prev5+next5)-8820*(prev3+next3)+39690*(prev1+next1)+0x8000)>>16;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
#endif

	//orginal CDF 5/3 DWT
#if 1
	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//linear predictor (deviation from linearity)
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update (smoothing)
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;
#endif


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf53_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf53_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf53_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	-0x1960C,	//-1.58613434342059f,	//alpha
	-0x00D90,	//-0.0529801185729f,	//beta
	 0x0E206,	// 0.8829110755309f,	//gamma
	 0x07189,	// 0.4435068520439f,	//delta
	 0x1264C,	// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]-=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_update(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unupdate(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_scale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=buf[k]*coeff>>16;
}
static void dwt1d_u8_unscale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=(buf[k]<<16)/coeff;
}
void dwt1d_cdf97_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[3]);
	//dwt1d_u8_scale(b2, count, cdf97_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	//dwt1d_u8_unscale(b2, count, cdf97_coeff[4]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

void dwt2d_dec_fwd(char *buffer, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int w2=iw, h2=ih;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	free(temp);
}
void dwt2d_dec_inv(char *buffer, int iw, int ih)
{
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 3);
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[2].w;++kx)//vertical invDWT
			dwt1d_haar_inv(buffer+4*kx+kc, psizes[2].h, 4*iw, temp);
		for(int ky=0;ky<psizes[2].h;++ky)//horizontal invDWT
			dwt1d_haar_inv(buffer+4*iw*ky+kc, psizes[2].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[1].w;++kx)//vertical invDWT
			dwt1d_cdf53_inv(buffer+4*kx+kc, psizes[1].h, 4*iw, temp);
		for(int ky=0;ky<psizes[1].h;++ky)//horizontal invDWT
			dwt1d_cdf53_inv(buffer+4*iw*ky+kc, psizes[1].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[0].w;++kx)//vertical invDWT
			dwt1d_cdf97_inv(buffer+4*kx+kc, psizes[0].h, 4*iw, temp);
		for(int ky=0;ky<psizes[0].h;++ky)//horizontal invDWT
			dwt1d_cdf97_inv(buffer+4*iw*ky+kc, psizes[0].w, 4, temp);
	}
	free(temp);
	array_free(&sizes);
}


//DCTs

static void dct4_fwd_i8(char *x)
{
	x[3]=x[0]-x[3];
	x[0]-=x[3]>>1;
	x[2]=x[1]-x[2];
	x[1]-=x[2]>>1;
	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;
	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;
}
static void dct4_inv_i8(char *x)
{
	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];
	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];
	x[1]+=x[2]>>1;
	x[2]=x[1]-x[2];
	x[0]+=x[3]>>1;
	x[3]=x[0]-x[3];
}
void image_dct4_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[ kx>>2           ]=x[0];
				temp[(kx>>2)+(iw>>2)  ]=x[1];
				temp[(kx>>2)+(iw>>2)*2]=x[2];
				temp[(kx>>2)+(iw>>2)*3]=x[3];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[(ky>>2)          ]=x[0];
				temp[(ky>>2)+(ih>>2)  ]=x[1];
				temp[(ky>>2)+(ih>>2)*2]=x[2];
				temp[(ky>>2)+(ih>>2)*3]=x[3];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct4_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>2)          ],
					temp[(ky>>2)+(ih>>2)  ],
					temp[(ky>>2)+(ih>>2)*2],
					temp[(ky>>2)+(ih>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>2)          ],
					temp[(kx>>2)+(iw>>2)  ],
					temp[(kx>>2)+(iw>>2)*2],
					temp[(kx>>2)+(iw>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
}

static void dct8_fwd_i8(char *x)
{
	//binDCT-C7
	x[7]=x[0]-x[7];
	x[6]=x[1]-x[6];
	x[5]=x[2]-x[5];
	x[4]=x[3]-x[4];
	x[0]-=x[7]>>1;
	x[1]-=x[6]>>1;
	x[2]-=x[5]>>1;
	x[3]-=x[4]>>1;

	x[3]=x[0]-x[3];
	x[2]=x[1]-x[2];
	x[0]-=x[3]>>1;
	x[1]-=x[2]>>1;

	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;

	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;

	x[5]-=x[6]*13>>5;
	x[6]+=x[5]*11>>4;
	x[5]=(x[6]*15>>5)-x[5];

	x[5]=x[4]-x[5];
	x[6]=x[7]-x[6];
	x[4]-=x[5]>>1;
	x[7]-=x[6]>>1;

	x[4]=(x[7]*3>>4)-x[4];
	x[7]-=x[4]*3>>4;

	x[5]+=x[6]*11>>4;
	x[6]-=x[5]*15>>5;
}
static void dct8_inv_i8(char *x)
{
	//invBinDCT-C7
	x[6]+=x[5]*15>>5;
	x[5]-=x[6]*11>>4;

	x[7]+=x[4]*3>>4;
	x[4]=(x[7]*3>>4)-x[4];

	x[7]+=x[6]>>1;
	x[4]+=x[5]>>1;
	x[6]=x[7]-x[6];
	x[5]=x[4]-x[5];

	x[5]=(x[6]*15>>5)-x[5];
	x[6]-=x[5]*11>>4;
	x[5]+=x[6]*13>>5;

	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];

	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];

	x[1]+=x[2]>>1;
	x[0]+=x[3]>>1;
	x[2]=x[1]-x[2];
	x[3]=x[0]-x[3];

	x[3]+=x[4]>>1;
	x[2]+=x[5]>>1;
	x[1]+=x[6]>>1;
	x[0]+=x[7]>>1;
	x[4]=x[3]-x[4];
	x[5]=x[2]-x[5];
	x[6]=x[1]-x[6];
	x[7]=x[0]-x[7];
}
void image_dct8_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
					image[(idx+4)<<2|kc],
					image[(idx+5)<<2|kc],
					image[(idx+6)<<2|kc],
					image[(idx+7)<<2|kc],
				};

				//char y[8];
				//memcpy(y, x, 8);
				//dct8_fwd_i8(y);
				//dct8_inv_i8(y);
				//if(memcmp(x, y, 8))
				//	x[0]=y[0];

				dct8_fwd_i8(x);

				temp[ kx>>3           ]=x[0];
				temp[(kx>>3)+(iw>>3)  ]=x[1];
				temp[(kx>>3)+(iw>>3)*2]=x[2];
				temp[(kx>>3)+(iw>>3)*3]=x[3];
				temp[(kx>>3)+(iw>>3)*4]=x[4];
				temp[(kx>>3)+(iw>>3)*5]=x[5];
				temp[(kx>>3)+(iw>>3)*6]=x[6];
				temp[(kx>>3)+(iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
					image[(idx+iw*4)<<2|kc],
					image[(idx+iw*5)<<2|kc],
					image[(idx+iw*6)<<2|kc],
					image[(idx+iw*7)<<2|kc],
				};

				dct8_fwd_i8(x);

				temp[(ky>>3)          ]=x[0];
				temp[(ky>>3)+(ih>>3)  ]=x[1];
				temp[(ky>>3)+(ih>>3)*2]=x[2];
				temp[(ky>>3)+(ih>>3)*3]=x[3];
				temp[(ky>>3)+(ih>>3)*4]=x[4];
				temp[(ky>>3)+(ih>>3)*5]=x[5];
				temp[(ky>>3)+(ih>>3)*6]=x[6];
				temp[(ky>>3)+(ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct8_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>3)          ],
					temp[(ky>>3)+(ih>>3)  ],
					temp[(ky>>3)+(ih>>3)*2],
					temp[(ky>>3)+(ih>>3)*3],
					temp[(ky>>3)+(ih>>3)*4],
					temp[(ky>>3)+(ih>>3)*5],
					temp[(ky>>3)+(ih>>3)*6],
					temp[(ky>>3)+(ih>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
				image[(idx+iw*4)<<2|kc]=x[4];
				image[(idx+iw*5)<<2|kc]=x[5];
				image[(idx+iw*6)<<2|kc]=x[6];
				image[(idx+iw*7)<<2|kc]=x[7];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>3)          ],
					temp[(kx>>3)+(iw>>3)  ],
					temp[(kx>>3)+(iw>>3)*2],
					temp[(kx>>3)+(iw>>3)*3],
					temp[(kx>>3)+(iw>>3)*4],
					temp[(kx>>3)+(iw>>3)*5],
					temp[(kx>>3)+(iw>>3)*6],
					temp[(kx>>3)+(iw>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
				image[(idx+4)<<2|kc]=x[4];
				image[(idx+5)<<2|kc]=x[5];
				image[(idx+6)<<2|kc]=x[6];
				image[(idx+7)<<2|kc]=x[7];
			}
		}
#endif
	}
	free(temp);
}


void predict_dct3_prep(float *dct3, float *dct4)
{
	dct3[0]=2,          dct3[1]= 2, dct3[2]=2;
	dct3[3]=sqrtf(3.f), dct3[4]= 0, dct3[5]=-sqrtf(3.f);
	dct3[6]=1,          dct3[7]=-2, dct3[8]=1;

	for(int k=0;k<9;++k)
		dct3[k]/=3.f;

	float a=(float)cos(M_PI/8), b=(float)cos(M_PI/4), c=(float)cos(M_PI*3/8);
	dct4[ 0]=0.5f, dct4[ 1]= a, dct4[ 2]= b, dct4[ 3]= c;
	dct4[ 4]=0.5f, dct4[ 5]= a, dct4[ 6]=-b, dct4[ 7]=-c;
	dct4[ 8]=0.5f, dct4[ 9]=-a, dct4[10]=-b, dct4[11]= c;
	dct4[12]=0.5f, dct4[13]=-a, dct4[14]= b, dct4[15]=-c;
}
int  predict_dct3(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen, const float *dct3, const float *dct4)
{
	float
		left[]=
	{
		kx-3>=0?(float)buf[idx-bytestride*3]:0,
		kx-2>=0?(float)buf[idx-bytestride*2]:0,
		kx-1>=0?(float)buf[idx-bytestride  ]:0,
	},
		top []=
	{
		ky-3>=0?(float)buf[idx-rowlen*3]:0,
		ky-2>=0?(float)buf[idx-rowlen*2]:0,
		ky-1>=0?(float)buf[idx-rowlen  ]:0,
	},
		topleft[]=
	{
		kx-3>=0&&ky-3>=0?(float)buf[idx-(rowlen+bytestride)*3]:0,
		kx-2>=0&&ky-2>=0?(float)buf[idx-(rowlen+bytestride)*2]:0,
		kx-1>=0&&ky-1>=0?(float)buf[idx-(rowlen+bytestride)  ]:0,
	},
		topright[]=
	{
		kx+3<iw&&ky-3>=0?(float)buf[idx-(rowlen-bytestride)*3]:0,
		kx+2<iw&&ky-2>=0?(float)buf[idx-(rowlen-bytestride)*2]:0,
		kx+1<iw&&ky-1>=0?(float)buf[idx-(rowlen-bytestride)  ]:0,
	};
	float x[]=
	{
		(dct3[0]*left[0]+dct3[1]*left[1]+dct3[2]*left[2] + dct3[0]*top[0]+dct3[1]*top[1]+dct3[2]*top[2] + dct3[0]*topleft[0]+dct3[1]*topleft[1]+dct3[2]*topleft[2] + dct3[0]*topright[0]+dct3[1]*topright[1]+dct3[2]*topright[2])*0.25f,
		(dct3[3]*left[0]+dct3[4]*left[1]+dct3[5]*left[2] + dct3[3]*top[0]+dct3[4]*top[1]+dct3[5]*top[2] + dct3[3]*topleft[0]+dct3[4]*topleft[1]+dct3[5]*topleft[2] + dct3[3]*topright[0]+dct3[4]*topright[1]+dct3[5]*topright[2])*0.25f,
		0,
		(dct3[6]*left[0]+dct3[7]*left[1]+dct3[8]*left[2] + dct3[6]*top[0]+dct3[7]*top[1]+dct3[8]*top[2] + dct3[6]*topleft[0]+dct3[7]*topleft[1]+dct3[8]*topleft[2] + dct3[6]*topright[0]+dct3[7]*topright[1]+dct3[8]*topright[2])*0.25f,
	};
	x[2]=x[1];
	float pred=dct4[12]*x[0]+dct4[13]*x[1]+dct4[14]*x[2]+dct4[15]*x[3];
	pred=roundf(pred);
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return (int)pred;
}
void pred_dct3_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	float coeff[25];
	predict_dct3_prep(coeff, coeff+9);

	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_dct3(buf, iw, kx, ky, idx, bytestride, rowlen, coeff, coeff+9);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_dct3_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	float coeff[25];
	predict_dct3_prep(coeff, coeff+9);

	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_dct3(buf, iw, kx, ky, idx, bytestride, rowlen, coeff, coeff+9);

				buf[idx]+=pred;
			}
		}
	}
}


//other
void image_split_fwd(char *image, int iw, int ih)
{
	char *b2=(char*)malloc((size_t)iw*ih<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, image, (size_t)iw*ih<<2);
	//int maxdim=MAXVAR(iw, ih);
	//char *temp=(char*)malloc(maxdim);
	//if(!temp)
	//{
	//	LOG_ERROR("Allocation error");
	//	return;
	//}
	//memset(temp, 0, maxdim);
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih-1;ky+=2)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					b2[ idx      <<2|kc],
					b2[(idx+1   )<<2|kc],
					b2[(idx  +iw)<<2|kc],
					b2[(idx+1+iw)<<2|kc],
				};
#if 0
				for(int k=1;k<4;++k)//mini-CDF (insertion sort)
				{
					int L=0, R=k-1, mid, found=0;
					while(L<=R)
					{
						mid=(L+R)>>1;
						if(x[mid]<x[k])
							L=mid+1;
						else if(x[mid]>x[k])
							R=mid-1;
						else
						{
							found=1;
							break;
						}
					}
					if(!found)
						mid=L+(L<k&&x[L]<x[k]);
				}
#endif
#if 1
				char temp;
				if(x[0]>x[1])temp=x[0], x[0]=x[1], x[1]=temp;//mini-CDF (dedicated sort)
				if(x[2]>x[3])temp=x[2], x[2]=x[3], x[3]=temp;//https://stackoverflow.com/questions/6145364/sort-4-number-with-few-comparisons
				if(x[0]>x[2])temp=x[0], x[0]=x[2], x[2]=temp;
				if(x[1]>x[3])temp=x[1], x[1]=x[3], x[3]=temp;
				if(x[1]>x[2])temp=x[1], x[1]=x[2], x[2]=temp;
#endif
				int idx2=iw*(ky>>1)+(kx>>1);
				image[ idx2                    <<2|kc]=x[0];
				image[(idx2           +(iw>>1))<<2|kc]=x[1];
				image[(idx2+iw*(ih>>1)        )<<2|kc]=x[2];
				image[(idx2+iw*(ih>>1)+(iw>>1))<<2|kc]=x[3];
			}
		}
#if 0
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
				};

				temp[ kx>>1           ]=x[0];
				temp[(kx>>1)+(iw>>1)  ]=x[1];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 0
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
				};

				temp[(ky>>1)          ]=x[0];
				temp[(ky>>1)+(ih>>1)  ]=x[1];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(b2);
	//free(temp);
}
void image_split_inv(char *image, int iw, int ih)
{
	char *b2=(char*)malloc((size_t)iw*ih<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, image, (size_t)iw*ih<<2);
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih-1;ky+=2)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx, idx2=iw*(ky>>1)+(kx>>1);
				char x[]=
				{
					b2[ idx2                    <<2|kc],
					b2[(idx2           +(iw>>1))<<2|kc],
					b2[(idx2+iw*(ih>>1)        )<<2|kc],
					b2[(idx2+iw*(ih>>1)+(iw>>1))<<2|kc],
				};
				
#if 0
				char temp;
				if(x[0]>x[1])temp=x[0], x[0]=x[1], x[1]=temp;//mini-CDF (dedicated sort)
				if(x[2]>x[3])temp=x[2], x[2]=x[3], x[3]=temp;//https://stackoverflow.com/questions/6145364/sort-4-number-with-few-comparisons
				if(x[0]>x[2])temp=x[0], x[0]=x[2], x[2]=temp;
				if(x[1]>x[3])temp=x[1], x[1]=x[3], x[3]=temp;
				if(x[1]>x[2])temp=x[1], x[1]=x[2], x[2]=temp;
#endif
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+1   )<<2|kc]=x[1];
				image[(idx  +iw)<<2|kc]=x[2];
				image[(idx+1+iw)<<2|kc]=x[3];
			}
		}
	}
	free(b2);
}

static unsigned qhist[256]={0};
void channel_entropy(unsigned char *buf, int resolution, int nch, int bytestride, float *cr, int *usage)
{
	double entropy[4]={0};
	memset(usage, 0, 4*sizeof(int));
	for(int kc=0;kc<nch;++kc)
	{
		memset(qhist, 0, 256*sizeof(unsigned));
		for(int k=0, end=resolution*bytestride;k<end;k+=bytestride)
		{
			unsigned char val=buf[k+kc];
			++qhist[val];
		}
		for(int ks=0;ks<256;++ks)
		{
			unsigned freq=qhist[ks];
			if(freq)
			{
				double p=(double)freq/resolution;
				p*=0xFF00;
				++p;
				p/=0x10000;
				entropy[kc]-=p*log2(p);
				++usage[kc];
			}
		}
		cr[kc]=(float)(8/entropy[kc]);
	}
	
	//calculate csize with joint histogram
	unsigned *h2=(unsigned*)malloc(0x1000000*sizeof(unsigned));
	if(!h2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(h2, 0, 0x1000000*sizeof(unsigned));
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		++h2[color];
	}
	for(int k=0;k<0x1000000;++k)
		usage[3]+=h2[k]!=0;
	double csize=0;
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		double p=(double)h2[color]/resolution;
		p*=0xFFFFFFFF;
		++p;
		p/=0x100000000;
		double bitsize=-log2(p);
		csize+=bitsize;
	}
	free(h2);
	csize/=8;
	cr[3]=(float)(resolution*3/csize);
}
void jointhistogram(unsigned char *buf, int iw, int ih, int nbits, ArrayHandle *hist, int space_not_color)
{
	int nlevels=1<<nbits, hsize=nlevels*nlevels*nlevels;
	if(*hist)
		array_free(hist);
	ARRAY_ALLOC(unsigned, *hist, 0, hsize, 0, 0);
	//unsigned *hist=(unsigned*)malloc(hsize*sizeof(unsigned));
	unsigned *h=(unsigned*)hist[0]->data;

	int res=iw*ih;
	switch(space_not_color)
	{
	case 0://show correlation in color
		for(int k=0;k<res;++k)
		{
			unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
			int color=b<<(nbits<<1)|g<<nbits|r;

			++h[color];
		}
		break;
	case 1://show correlation in space x (CURR, W, WW)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=kx-2>=0?buf[(iw*ky+kx-2)<<2|1]>>(8-nbits):0,//WW
					v1=kx-1>=0?buf[(iw*ky+kx-1)<<2|1]>>(8-nbits):0,//W
					v0=kx-0>=0?buf[(iw*ky+kx-0)<<2|1]>>(8-nbits):0;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++h[color];
			}
		}
		break;
	case 2://show correlation in space x (CURR, N, NN)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=ky-2>=0?buf[(iw*(ky-2)+kx)<<2|1]>>(8-nbits):0,//NN
					v1=ky-1>=0?buf[(iw*(ky-1)+kx)<<2|1]>>(8-nbits):0,//N
					v0=ky-0>=0?buf[(iw*(ky-0)+kx)<<2|1]>>(8-nbits):0;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++h[color];
			}
		}
		break;
	case 3://show correlation in space x (CURR, N, W)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=kx-1>=0?buf[(iw* ky   +kx-1)<<2|1]>>(8-nbits):0,//W
					v1=ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|1]>>(8-nbits):0,//N
					v0=        buf[(iw* ky   +kx  )<<2|1]>>(8-nbits)  ;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++h[color];
			}
		}
		break;
	}

	//X  don't calculate csize from downsampled histogram
#if 0
	double csize=0;
	for(int k=0;k<res;++k)//calculate csize using Zipf's law
	{
		unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
		int color=b<<(nbits<<1)|g<<nbits|r;

		unsigned freq=h[color];
		double p=(double)freq/res;

		p*=0xFFFFFFFF;//guard against quantized zero
		++p;
		p/=0x100000000;

		csize-=log2(p);
	}
	csize/=nbits;
	float CR=(float)(res*3/csize);
#endif

#if 0
	float entropy=0;//also BPP
	for(int k=0;k<hsize;++k)//calculate entropy using Shannon's law
	{
		int freq=h[k];
		if(freq)
		{
			float p=(float)freq/res;
			entropy-=p*log2f(p);
		}
	}
	float CR=3*nbits/entropy;
#endif

	unsigned histmin=0;
	unsigned histmax=0;
	for(int k=0;k<hsize;++k)//get min & max
	{
		if(histmin>h[k])
			histmin=h[k];
		if(histmax<h[k])
			histmax=h[k];
	}

	unsigned char *h2=(unsigned char*)h;
	for(int k=0;k<hsize;++k)//normalize
		h2[k]=(h[k]<<8)/histmax;

	//for(int k=0;k<hsize;++k)
	//	h[k]=(unsigned)(((long long)h[k]<<32)/histmax);

	//for(int kz=0;kz<nlevels-1;++kz)
	//{
	//	for(int ky=0;ky<nlevels-1;++ky)
	//	{
	//		for(int kx=0;kx<nlevels-1;++kx)
	//		{
	//		}
	//	}
	//}

	//return CR;
}


E24Params e24_params[3]=
{
	{ 8, 26, 26, 26, 0xD4, 71},
	{23, 37, 37, 37, 0xD4, 86},
	{ 8, 26, 26, 26, 0xD4, 77},
};
double e24_cr[3]={0};
#define E24_SETPARAM(P, IDX, VAL) ((short*)(P))[IDX]=VAL
int e24_incparam(E24Params *p, int pidx, int step)
{
	int prevval=0;
	switch(pidx)
	{
	case 0:prevval=p->gwidth, p->gwidth+=step; if(p->gwidth<1)p->gwidth=1; break;
	case 1:prevval=p->mleft , p->mleft +=step; if(p->mleft <0)p->mleft =0; break;
	case 2:prevval=p->mtop  , p->mtop  +=step; if(p->mtop  <0)p->mtop  =0; break;
	case 3:prevval=p->mright, p->mright+=step; if(p->mright<0)p->mright=0; break;
	case 4:prevval=p->alpha , p->alpha +=step; if(p->alpha <0)p->alpha =0; else if(p->alpha>0xFF)p->alpha=0xFF; break;
	case 5:prevval=p->maxinc, p->maxinc+=step; if(p->maxinc<1)p->maxinc=1; break;
	}
	return prevval;
}
void e24_normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	if(!nsymbols)//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
		return;
	}
	unsigned sum=0, qfreq;
	for(int sym=0;sym<nlevels;++sym)
	{
		qfreq=((long long)srchist[sym]<<16)/nsymbols;
		CDF[sym]=sum;
		sum+=qfreq;
	}
}
void e24_addhist(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0a, int x0b, int y0, int maxinc, unsigned *CDF2)
{
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(iw*ky+kx)<<2|kc];
			int dist=abs(ky-y0);
			if(kx<x0a)
				dist+=abs(kx-x0a);
			else if(kx>x0b)
				dist+=abs(kx-x0b);
			int inc=maxinc-dist;
			if(inc>0)
			{
				CDF2[sym]+=inc;
				CDF2[256]+=inc;
			}
		}
	}
}
double e24_calcloss(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, E24Params const *p, const unsigned short *CDF0, unsigned *CDF2)
{
	double chsize=0;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;)
		{
			int xend=kx+p->gwidth;
			if(xend>x2)
				xend=x2;
			memset(CDF2, 0, 257*sizeof(unsigned));
			if(p->mtop)
				e24_addhist(buf, iw, ih, kc, kx-p->mleft, xend+p->mright, ky-p->mtop, ky, kx, xend, ky, p->maxinc, CDF2);
			if(p->mleft)
				e24_addhist(buf, iw, ih, kc, kx-p->mleft, kx, ky, ky+1, kx, xend, ky, p->maxinc, CDF2);

			int overflow=0;
			int sum, cdf1, f1, f2, freq;
			if(CDF2[256])
			{
				sum=0;
				for(int sym=0;sym<256;++sym)
				{
					cdf1=!overflow?CDF0[sym]:0x10000;
					if(sym<255)
						overflow|=cdf1>CDF0[sym+1];
					f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

					f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize

					freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

					freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
					if(freq<0||freq>0xFF01)
						LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
					CDF2[sym]=sum;
					sum+=freq;
					if(sum>0x10000)
						LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
				}
			}
			else
			{
				for(int sym=0;sym<256;++sym)
				{
					if(overflow)
						CDF2[sym]=0xFF00|sym;
					else
					{
						int cdf=CDF0[sym];
						CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
						if(sym<255)
							overflow|=cdf>CDF0[sym+1];
					}
				}
			}
			CDF2[256]=0x10000;
				
			int kx2=kx;
			for(;kx2<kx+p->gwidth;++kx2)
			{
				unsigned char sym=buf[(iw*ky+kx2)<<2|kc];
				int freq=CDF2[sym+1]-CDF2[sym];
				double p=(double)freq/0x10000, bitsize=-log2(p);//Zipf's law
				chsize+=bitsize;
			}
			kx+=p->gwidth;
		}
	}
	return chsize;
}
double e24_optimize(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, E24Params *p, int pidx, int step, const unsigned short *CDF0, unsigned *CDF2)
{
	double csize00, csize0, csize;
	int prevval0, prevval;

	prevval0=((short*)p)[pidx];
	csize00=csize=e24_calcloss(buf, iw, ih, kc, x1, x2, y1, y2, p, CDF0, CDF2);
		
	int subit;
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		prevval=e24_incparam(p, pidx, step);
		csize=e24_calcloss(buf, iw, ih, kc, x1, x2, y1, y2, p, CDF0, CDF2);
		if(csize>=csize0)//cancel last change and break
		{
			E24_SETPARAM(p, pidx, prevval);
			csize=csize0;
			break;
		}
	}
		
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		prevval=e24_incparam(p, pidx, -step);
		csize=e24_calcloss(buf, iw, ih, kc, x1, x2, y1, y2, p, CDF0, CDF2);
		if(csize>=csize0)
		{
			E24_SETPARAM(p, pidx, prevval);
			csize=csize0;
			break;
		}
	}
	if(csize>=csize00)//prevent CR from worsening
	{
		E24_SETPARAM(p, pidx, prevval0);
		csize=csize00;
	}
	return csize;
}
int e24_optimizeall(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int loud)
{
	int res=iw*ih;
	unsigned short *CDF0=(unsigned short*)malloc(256LL*3*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++CDF2[sym];
		}
		e24_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}
	
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	
	double csizes[3]={0};
	int steps[]={4, 2, 1};
	int usize=(x2-x1)*(y2-y1);
	for(int kc=0;kc<3;++kc)
	{
		for(int ks=0;ks<COUNTOF(steps);++ks)
		{
			for(int it=0, improve=1;it<64&&improve;++it)
			{
				improve=0;
				for(int pidx=0;pidx<sizeof(e24_params)/sizeof(e24_params->gwidth);++pidx)
				{
					csizes[kc]=e24_optimize(buf, iw, ih, kc, x1, x2, y1, y2, e24_params+kc, pidx, steps[ks], CDF0+((size_t)kc<<8), CDF2);
					if(loud)
						io_render();
				}
			}
		}
		csizes[kc]/=8;
		e24_cr[kc]=usize>0&&csizes[kc]?(x2-x1)*(y2-y1)/csizes[kc]:0;
	}

	free(CDF0);
	free(CDF2);
	return 1;
}
void e24_estimate(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2)
{
	int res=iw*ih;
	unsigned short *CDF0=(unsigned short*)malloc(256LL*3*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++CDF2[sym];
		}
		e24_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}
	
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;

	double csizes[3];
	int usize=(x2-x1)*(y2-y1);
	csizes[0]=e24_calcloss(buf, iw, ih, 0, x1, x2, y1, y2, e24_params  , CDF0, CDF2);
	csizes[1]=e24_calcloss(buf, iw, ih, 1, x1, x2, y1, y2, e24_params+1, CDF0, CDF2);
	csizes[2]=e24_calcloss(buf, iw, ih, 2, x1, x2, y1, y2, e24_params+2, CDF0, CDF2);

	csizes[0]/=8;
	csizes[1]/=8;
	csizes[2]/=8;

	e24_cr[0]=usize>0&&csizes[0]?(x2-x1)*(y2-y1)/csizes[0]:0;
	e24_cr[1]=usize>0&&csizes[1]?(x2-x1)*(y2-y1)/csizes[1]:0;
	e24_cr[2]=usize>0&&csizes[2]?(x2-x1)*(y2-y1)/csizes[2]:0;

	free(CDF0);
	free(CDF2);
}