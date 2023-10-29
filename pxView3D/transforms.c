#include"pxview3d.h"
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include"intercept_malloc.h"
#include"lodepng.h"//
#include<immintrin.h>
static const char file[]=__FILE__;


unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);

int fast_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-31;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((__m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((__m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
	__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
	s2=_mm_add_epi32(s2, hi);
	s2=_mm_hadd_epi32(s2, s2);
	int s3=_mm_extract_epi32(s2, 0);
	for(;k<count;++k)
		s3+=a[k]*b[k];
	return s3;
}
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
void colortransform_ycmcb_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;       //diff(r, g)            [ 1      -1      0  ].RGB
		g+=r>>1;
		b-=g;       //diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB
		g+=b>>1;    //av(b, av(r, g))       [ 1/4     1/4    1/2].RGB

		buf[k  ]=r;//Cm
		buf[k|1]=g;//Y
		buf[k|2]=b;//Cb
	}
}
void colortransform_ycmcb_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];//Cm Y Cb
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_subg_fwd(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;       //r-g							[1     -1     0  ].RGB
		b-=g;       //b-g							[0     -1     1  ].RGB
		g+=(r+b)>>2;//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB

		buf[k  ]=r;//C?
		buf[k|1]=g;//Y
		buf[k|2]=b;//C?
	}
}
void colortransform_subg_inv(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];//Cm Y Cb
		
		g-=(r+b)>>2;
		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void lossy_colortransform_ycbcr(char *buf, int iw, int ih, int fwd)
{
	for(int ky=0;ky<ih;++ky)//for demonstration purposes only
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			
			if(fwd)
			{
				double
					r=(double)(buf[idx|0]+128),
					g=(double)(buf[idx|1]+128),
					b=(double)(buf[idx|2]+128);

				double
					Y=0.299*r+0.587*g+0.114*b,
					Cb=128-0.168736*r-0.331264*g+0.5*b,
					Cr=128+0.5*r-0.418688*g-0.081312*b;

				Y=CLAMP(0, Y, 255);
				Cb=CLAMP(0, Cb, 255);
				Cr=CLAMP(0, Cr, 255);

				buf[idx|0]=(char)Y-128;
				buf[idx|1]=(char)Cb-128;
				buf[idx|2]=(char)Cr-128;
			}
			else
			{
				double
					Y=(double)buf[idx|0],
					Cb=(double)buf[idx|1],
					Cr=(double)buf[idx|2];

				double
					r=Y+1.402*Cr,
					g=Y-0.344136*Cb-0.714136*Cr,
					b=Y+1.772*Cb;

				buf[idx|0]=(char)r;
				buf[idx|1]=(char)g;
				buf[idx|2]=(char)b;
			}
		}
	}
}
void lossy_colortransform_xyb(char *buf, int iw, int ih, int fwd)
{
	const double bias=0.00379307325527544933;
	double cbrt_bias=cbrt(bias);
	for(int ky=0;ky<ih;++ky)//for demonstration purposes only
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;

			if(fwd)
			{
				double
					r=(double)(buf[idx|0]+128)/255,
					g=(double)(buf[idx|1]+128)/255,
					b=(double)(buf[idx|2]+128)/255;

				double
					L=0.3*r+0.622*g+0.078*b,
					M=0.23*r+0.692*g+0.078*b,
					S=0.24342268924547819*r+0.20476744424496821*g+0.55180986650955360*b;
				L=cbrt(L+bias)-cbrt_bias;
				M=cbrt(M+bias)-cbrt_bias;
				S=cbrt(S+bias)-cbrt_bias;
				r=(L-M)*0.5;//X
				g=(L+M)*0.5;//Y
				b=S-g;		//B-Y

				r*=4096;//customparam_ct[0]*1000
				g*=144;
				b*=1024;
				r=CLAMP(-128, r, 127);
				g=CLAMP(-128, g, 127);
				b=CLAMP(-128, b, 127);

				buf[idx|0]=(char)r;
				buf[idx|1]=(char)g;
				buf[idx|2]=(char)b;
			}
			else
			{
				double
					r=(double)buf[idx|0]/4096,
					g=(double)buf[idx|1]/144,
					b=(double)buf[idx|2]/1024;

				double
					L=g+r,
					M=g-r,//g is the sum
					S=b+g;
				L+=cbrt_bias, L*=L*L, L-=bias;
				M+=cbrt_bias, M*=M*M, M-=bias;
				S+=cbrt_bias, S*=S*S, S-=bias;
				r=11.03160*L-9.86694*M-0.164623*S;
				g=-3.25415*L+4.41877*M-0.164623*S;
				b=-3.65885*L+2.71292*M+1.945930*S;

				buf[idx|0]=(char)(255*r-128);
				buf[idx|1]=(char)(255*g-128);
				buf[idx|2]=(char)(255*b-128);
			}
		}
	}
}

void colortransform_quad(char *src, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);
	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			int idx=(iw*ky+kx)<<2;
			char comp[]=
			{
				src[idx  ], src[idx+4  ], src[idx+(iw<<2)  ], src[idx+(iw<<2)+4  ],//r
				src[idx+1], src[idx+4+1], src[idx+(iw<<2)+1], src[idx+(iw<<2)+4+1],//g
				src[idx+2], src[idx+4+2], src[idx+(iw<<2)+2], src[idx+(iw<<2)+4+2],//b
			};

			if(fwd)
			{
				for(int k=0;k<4;++k)
				{
					comp[k+0]-=comp[k+4];
					comp[k+4]+=comp[k+0]>>1;
					comp[k+8]-=comp[k+4];
					comp[k+4]+=comp[k+8]>>1;
				}
				for(int k=0;k<12;k+=4)
				{
					comp[k|1]-=comp[k|0];
					comp[k|0]+=comp[k|1]>>1;
					comp[k|3]-=comp[k|2];
					comp[k|2]+=comp[k|3]>>1;
					comp[k|2]-=comp[k|0];
					comp[k|0]+=comp[k|2]>>1;
				}
			}
			else
			{
				for(int k=0;k<12;k+=4)
				{
					comp[k|0]-=comp[k|2]>>1;
					comp[k|2]+=comp[k|0];
					comp[k|2]-=comp[k|3]>>1;
					comp[k|3]+=comp[k|2];
					comp[k|0]-=comp[k|1]>>1;
					comp[k|1]+=comp[k|0];
				}
				for(int k=0;k<4;++k)
				{
					comp[k+4]-=comp[k+8]>>1;
					comp[k+8]+=comp[k+4];
					comp[k+4]-=comp[k+0]>>1;
					comp[k+0]+=comp[k+4];
				}
			}

			dst[(iw*(ky>>1)+(kx>>1))<<2|0]=comp[0|0];
			dst[(iw*(ky>>1)+(kx>>1))<<2|1]=comp[4|0];
			dst[(iw*(ky>>1)+(kx>>1))<<2|2]=comp[8|0];

			dst[(iw*(ky>>1)+(iw>>1)+(kx>>1))<<2|0]=comp[0|1];
			dst[(iw*(ky>>1)+(iw>>1)+(kx>>1))<<2|1]=comp[4|1];
			dst[(iw*(ky>>1)+(iw>>1)+(kx>>1))<<2|2]=comp[8|1];

			dst[(iw*((ih>>1)+(ky>>1))+(kx>>1))<<2|0]=comp[0|2];
			dst[(iw*((ih>>1)+(ky>>1))+(kx>>1))<<2|1]=comp[4|2];
			dst[(iw*((ih>>1)+(ky>>1))+(kx>>1))<<2|2]=comp[8|2];

			dst[(iw*((ih>>1)+(ky>>1))+(iw>>1)+(kx>>1))<<2|0]=comp[0|3];
			dst[(iw*((ih>>1)+(ky>>1))+(iw>>1)+(kx>>1))<<2|1]=comp[4|3];
			dst[(iw*((ih>>1)+(ky>>1))+(iw>>1)+(kx>>1))<<2|2]=comp[8|3];

			//char r[]={src[idx  ], src[idx+4  ], src[idx+(iw<<2)  ], src[idx+(iw<<2)+4  ]};
			//char g[]={src[idx+1], src[idx+4+1], src[idx+(iw<<2)+1], src[idx+(iw<<2)+4+1]};
			//char b[]={src[idx+2], src[idx+4+2], src[idx+(iw<<2)+2], src[idx+(iw<<2)+4+2]};
			//
			//for(int k=0;k<4;++k)
			//{
			//	r[k]-=g[k];       //diff(r, g)            [ 1      -1      0  ].RGB
			//	g[k]+=r[k]>>1;
			//	b[k]-=g[k];       //diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB
			//	g[k]+=b[k]>>1;    //av(b, av(r, g))       [ 1/4     1/4    1/2].RGB
			//}
			//r[0]-=r[1];
			//r[1]+=r[0]>>1;
			//r[2]-=
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
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

int pythagoras(int a, int b)
{
	a*=a;
	b*=b;
	a+=b;
	return a;
}
void act_fwd(short *coeff, char *r, char *g, char *b)
{
	*r+=(coeff[ 0]**g+coeff[ 1]**b+128)>>8;
	*g+=(coeff[ 2]**r+coeff[ 3]**b+128)>>8;
	*b+=(coeff[ 4]**r+coeff[ 5]**g+128)>>8;
	*r+=(coeff[ 6]**g+coeff[ 7]**b+128)>>8;
	*g+=(coeff[ 8]**r+coeff[ 9]**b+128)>>8;
	*b+=(coeff[10]**r+coeff[11]**g+128)>>8;
}
void act_train(short *coeff, char r0, char g0, char b0, int e0)
{
	int idx[]=
	{
		xoroshiro128_next()%12,
		xoroshiro128_next()%12,
		xoroshiro128_next()%12,
	};
	int inc[]=
	{
		((xoroshiro128_next()&1)<<1)-1,
		((xoroshiro128_next()&1)<<1)-1,
		((xoroshiro128_next()&1)<<1)-1,
		//(xoroshiro128_next()&63)-32,
		//(xoroshiro128_next()&63)-32,
		//(xoroshiro128_next()&63)-32,
	};

	coeff[idx[0]]+=inc[0];
	coeff[idx[1]]+=inc[1];
	coeff[idx[2]]+=inc[2];

	char r=r0, g=g0, b=b0;
	//act_fwd(coeff, &r, &g, &b);
	r+=(coeff[ 0]*g+coeff[ 1]*b+128)>>8;
	g+=(coeff[ 2]*r+coeff[ 3]*b+128)>>8;
	b+=(coeff[ 4]*r+coeff[ 5]*g+128)>>8;
	r+=(coeff[ 6]*g+coeff[ 7]*b+128)>>8;
	g+=(coeff[ 8]*r+coeff[ 9]*b+128)>>8;
	b+=(coeff[10]*r+coeff[11]*g+128)>>8;

	int error=pythagoras(r, b);
	if(error>=e0)
	{
		coeff[idx[0]]-=inc[0];
		coeff[idx[1]]-=inc[1];
		coeff[idx[2]]-=inc[2];
	}
}
void colortransform_adaptive(char *src, int iw, int ih, int fwd)//does mediocre job at spatially predicting chroma, doesn't spatially predict luma
{
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	int *blurrer=(int*)malloc(iw*4*sizeof(int));
	if(!dst||!blurrer)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	const char *pixels=fwd?src:dst;
	int idx;
	memcpy(dst, src, (size_t)res<<2);
	memset(blurrer, 0, iw*3*sizeof(int));
	for(int ky=0;ky<ih;++ky)
	{
		int predr=0, predg=0, predb=0;
		for(int kx=0;kx<iw;++kx)
		{
#if 0
			int predr=0, predg=0, predb=0, c=0;
			if(kx)
			{
				idx=(iw*ky+kx-1)<<2;
				++c;
				predr+=pixels[idx|0];
				predg+=pixels[idx|1];
				predb+=pixels[idx|2];
			}
			if(ky)
			{
				idx=(iw*(ky-1)+kx)<<2;
				++c;
				predr+=pixels[idx|0];
				predg+=pixels[idx|1];
				predb+=pixels[idx|2];
				if(kx)
				{
					idx=(iw*(ky-1)+kx-1)<<2;
					++c;
					predr+=pixels[idx|0];
					predg+=pixels[idx|1];
					predb+=pixels[idx|2];
				}
				if(kx<iw-1)
				{
					idx=(iw*(ky-1)+kx+1)<<2;
					++c;
					predr+=pixels[idx|0];
					predg+=pixels[idx|1];
					predb+=pixels[idx|2];
				}
			}
#endif
			idx=(iw*ky+kx)<<2;
			char r=src[idx|0], g=src[idx|1], b=src[idx|2];
#if 0
			if(fwd)
			{
				if(predg)
					r-=g*predr/predg;
				else
					r-=g;
				g+=r/2;
				b-=g;
				g+=b>>1;
			}
			else
			{
				g-=b>>1;
				b+=g;
				g-=r/2;
				if(predg)
					r+=g*predr/predg;
				else
					r+=g;
			}
#endif
#if 1
			if(ky)//N
			{
				predr=blurrer[kx<<2|0];
				predg=blurrer[kx<<2|1];
				predb=blurrer[kx<<2|2];
			}
			if(kx)//W
			{
				predr=blurrer[(kx-1)<<2|0];
				predg=blurrer[(kx-1)<<2|1];
				predb=blurrer[(kx-1)<<2|2];
			}
			char r0, g0, b0;
			if(fwd)
			{
				r0=r;
				g0=g;
				b0=b;

				if(predg&&abs(predr/predg)<127)
					r-=(int)((long long)g*predr/predg);
				else
					r-=g;

				if(predg&&abs(predb/predg)<127)
					b-=(int)((long long)g*predb/predg);
				else
					b-=g;

				g+=(r+b)>>2;
			}
			else
			{

				g-=(r+b)>>2;

				if(predg&&abs(predb/predg)<127)
					b+=(int)((long long)g*predb/predg);
				else
					b+=g;

				if(predg&&abs(predr/predg)<127)
					r+=(int)((long long)g*predr/predg);
				else
					r+=g;

				r0=r;
				g0=g;
				b0=b;
			}
			blurrer[kx<<2|0]=((r0<<16)-predr)>>3;
			blurrer[kx<<2|1]=((g0<<16)-predg)>>3;
			blurrer[kx<<2|2]=((b0<<16)-predb)>>3;
			//predr+=((r0<<16)-predr)>>4;
			//predg+=((g0<<16)-predg)>>4;
			//predb+=((b0<<16)-predb)>>4;
#endif
			dst[idx|0]=r;
			dst[idx|1]=g;
			dst[idx|2]=b;
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);

#if 0
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	short *coeff=(short*)malloc(iw*sizeof(short[9]));
	if(!dst||!coeff)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);
	//int black=0xFF000000;
	//memfill(dst, &black, (size_t)res<<2, 4);
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	short coeff0[]=
	{
		 0x0100,  0x0000, -0x0100,
		 0x0040,  0x0080,  0x0040,
		-0x0080,  0x0100, -0x0080,
	};
	//short coeff0[12]=
	//{
	//	-0x0100,  0x0000,
	//	 0x0080,  0x0000,
	//	 0x0000, -0x0100,
	//	 0x0000,  0x0000,
	//	 0x0000,  0x0080,
	//	 0x0000,  0x0000,
	//};
	memfill(coeff, coeff0, iw*sizeof(short[9]), sizeof(short[9]));
	XOROSHIRO128_RESET();
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			char
				rN=0, gN=0, bN=0,
				rW=0, gW=0, bW=0;
			if(kx)
			{
				rW=pixels[(idx-(1<<2))|0];
				gW=pixels[(idx-(1<<2))|1];
				bW=pixels[(idx-(1<<2))|2];
			}
			if(ky)
			{
				rN=pixels[(idx-(iw<<2))|0];
				gN=pixels[(idx-(iw<<2))|1];
				bN=pixels[(idx-(iw<<2))|2];
			}
			char
				pr=(rN+rW)>>1,
				pg=(gN+gW)>>1,
				pb=(bN+bW)>>1;
			char r=src[idx], g=src[idx|1], b=src[idx|2];

			if(fwd)
			{
				r-=g;
				//r-=pg?(int)(g*sqrt((double)pr/pg)):g;//X  noisy: inefficient with spatial predictor
				g+=r>>1;
				b-=g;
				g+=b>>1;
			}
			else
			{
				g-=b>>1;
				b+=g;
				g-=r>>1;
				r+=g;
				//r+=pg?(int)(g*sqrt((double)pr/pg)):g;//X
			}

			dst[idx|0]=r;
			dst[idx|1]=g;
			dst[idx|2]=b;
#if 0
			char
				rN=0, gN=0, bN=0,
				rW=0, gW=0, bW=0;
			int eN=0,
				eW=0;
			if(kx)
			{
				rW=pixels[(idx-(1<<2))|0];
				gW=pixels[(idx-(1<<2))|1];
				bW=pixels[(idx-(1<<2))|2];
				eW=pythagoras(errors[(idx-(1<<2))|0], errors[(idx-(1<<2))|2]);
			}
			if(ky)
			{
				rN=pixels[(idx-(iw<<2))|0];
				gN=pixels[(idx-(iw<<2))|1];
				bN=pixels[(idx-(iw<<2))|2];
				eN=pythagoras(errors[(idx-(iw<<2))|0], errors[(idx-(iw<<2))|2]);
			}
			for(int k=0;k<1;++k)
			{
				if(kx)
					act_train(coeff, rW, gW, bW, eW);
				if(ky)
					act_train(coeff, rN, gN, bN, eN);
			}

			char r=src[idx], g=src[idx|1], b=src[idx|2];
			char y, c0, c1;
			if(fwd)
			{
				y =(r*coeff0[0]+g*coeff0[1]+b*coeff0[2]+128)>>8;
				c0=(r*coeff0[3]+g*coeff0[4]+b*coeff0[5]+128)>>8;
				c1=(r*coeff0[6]+g*coeff0[7]+b*coeff0[8]+128)>>8;
				//r+=(coeff[ 0]*g+coeff[ 1]*b+128)>>8;
				//g+=(coeff[ 2]*r+coeff[ 3]*b+128)>>8;
				//b+=(coeff[ 4]*r+coeff[ 5]*g+128)>>8;
				//r+=(coeff[ 6]*g+coeff[ 7]*b+128)>>8;
				//g+=(coeff[ 8]*r+coeff[ 9]*b+128)>>8;
				//b+=(coeff[10]*r+coeff[11]*g+128)>>8;
			}
			else
			{
				b-=(coeff[10]*r+coeff[11]*g+128)>>8;
				g-=(coeff[ 8]*r+coeff[ 9]*b+128)>>8;
				r-=(coeff[ 6]*g+coeff[ 7]*b+128)>>8;
				b-=(coeff[ 4]*r+coeff[ 5]*g+128)>>8;
				g-=(coeff[ 2]*r+coeff[ 3]*b+128)>>8;
				r-=(coeff[ 0]*g+coeff[ 1]*b+128)>>8;
			}
			dst[idx|0]=r;
			dst[idx|1]=g;
			dst[idx|2]=b;

			//int e0=pythagoras(r, b);
			//const char r0=pixels[idx], g0=pixels[idx|1], b0=pixels[idx|2];
			//for(int k=0;k<2;++k)
			//{
			//	//train
			//}
#endif

#if 0
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

			char rcurr=src[idx<<2], gcurr=src[idx<<2|1], bcurr=src[idx<<2|2];

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
#endif
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
#endif
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


double lossygrad_rmse[3]={0}, lossygrad_psnr[3]={0};
void preproc_grad(char *src, int iw, int ih)
{
	memset(lossygrad_rmse, 0, sizeof(lossygrad_rmse));
	memset(lossygrad_psnr, 0, sizeof(lossygrad_psnr));
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);//copy alpha
	long long MSE[3]={0};
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|kc]:0

				//calc sdev
				int sdev=0;
				int idx=(iw*ky+kx)<<2|kc;
				int curr=src[idx];
				char vals[]=
				{
					LOAD(src,  0, -1),
					LOAD(src, -1,  0),
					LOAD(src, -1, -1),
					curr
				};
				int mean=0;
				for(int k=0;k<_countof(vals);++k)
					mean+=vals[k];
				mean>>=2;
				for(int k=0;k<_countof(vals);++k)
				{
					int dev=vals[k]-mean;
					sdev+=dev*dev;
				}
				sdev/=3;
				sdev=(int)round(sqrt(sdev));

				int tolerance=sdev;//tolerance is proportional to sdev of neighborhood
				if(tolerance>4)
					tolerance=4;
				
				if(kx==(iw>>1)&&ky==(ih>>1))//
					printf("");
				
				//calc pred
				char
					N =LOAD(dst,  0, -1),
					W =LOAD(dst, -1,  0),
					NW=LOAD(dst, -1, -1);
				char vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				int pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);

				//select quantization
				int delta=curr-pred, qdelta, recon;
				int kb=8;
				do//{7, ...1} because bit 0 can't be quantized
				{
					--kb;

					//quantize
					qdelta=delta+(1<<(kb-1));//add half
					qdelta>>=kb;
					qdelta<<=kb;
					//qdelta&=~((1<<kb)-1);

					recon=pred+qdelta;
					if(recon<-128||recon>127)
						continue;
				}while(abs(curr-recon)>tolerance&&kb>=1);
				qdelta=delta;
				if(kb>0)
				{
					qdelta+=(1<<(kb-1));//add half
					qdelta>>=kb;
					qdelta<<=kb;
				}
				recon=pred+qdelta;
				recon=CLAMP(-128, recon, 127);
				//if(recon<-128||recon>127)
				//	LOG_ERROR("Overflow");

				int error=curr-recon;
				MSE[kc]+=error*error;

				dst[idx]=recon;
#undef  LOAD
			}
		}
		lossygrad_rmse[kc]=sqrt((double)MSE[kc]/res);//RMSE is the Euclidean distance (l2 norm)
		lossygrad_psnr[kc]=20*log10(255/lossygrad_rmse[kc]);
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
}

#define X_BLOCKSIZE 3
void preproc_x(char *src, int iw, int ih)//https://github.com/Equationist/xpng/blob/master/main.c
{
	int res=iw*ih;
	unsigned char *noise=(unsigned char*)malloc((size_t)res<<2);
	int *freq=(int*)malloc(768*sizeof(int));
	char *dst=(char*)malloc((size_t)res<<2);
	char *diff=(char*)malloc((size_t)res<<2);
	if(!noise||!freq||!dst||!diff)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int kc=0;kc<3;++kc)//calc noise
	{
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				int n=0;
				long long x=0, x2=0;
				for(int ky2=ky-X_BLOCKSIZE/2;ky2<=ky+X_BLOCKSIZE/2;++ky2)
				{
					if((unsigned)ky2<(unsigned)ih)
					{
						for(int kx2=ky-X_BLOCKSIZE/2;kx2<=kx+X_BLOCKSIZE/2;++kx2)
						{
							if((unsigned)kx2<(unsigned)iw)
							{
								char px=src[(iw*ky2+kx2)<<2|kc];
								x+=px;
								x2+=px*px;
								++n;
							}
						}
					}
				}
				noise[(iw*ky+kx)<<2|kc]=(int)round(sqrt((double)(x2-x*x)/n));
			}
		}
	}

	for(int step=1;step<=128;step<<=1)//set preference levels
	{
		for(int k=-128;k<128;k+=step)
			freq[(unsigned char)k]=step;
	}
	memfill(freq+256, freq, 512*sizeof(int), 256*sizeof(int));
	
	memcpy(dst, src, (size_t)res<<2);//copy alpha

	const int clevel=6;//0~40
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih;++ky)
		{
			//long long refdist=1;
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|kc]:0
				char
					N =LOAD(dst,  0, -1),
					W =LOAD(dst, -1,  0),
					NW=LOAD(dst, -1, -1);
				char vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				int pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
				pred+=128;
#undef  LOAD
				int idx=(iw*ky+kx)<<2|kc;
				int curr=src[idx];
				curr+=128;

				int tolerance=clevel+noise[idx]*clevel/(10+clevel);

				int ltarget=0, rtarget=255;
				if(tolerance<curr)
					ltarget=curr-tolerance;
				if(curr<255-tolerance)
					rtarget=curr+tolerance;

				//find best quantization within allowable tolerance
				int approx=curr;
				long long f=0;
				int a=ltarget;
				do
				{
					unsigned char d=a-pred;
					//if((unsigned)(kc<<8|d)>=768)
					//	LOG_ERROR("OOB");
					if(freq[kc<<8|d]>f)
					{
						approx=a;
						f=freq[kc<<8|d];
					}
				}while(a++!=rtarget);
				diff[idx]=approx-pred;
				dst[idx]=approx-128;

				//if(abs(curr-(pred+(unsigned char)diff[idx-refdist]))<tolerance)//RLE-related
				//{
				//	diff[idx]=diff[idx-refdist];
				//	dst[idx]=pred+diff[idx];
				//	continue;
				//}
				//for(size_t rd=1;rd<(kx<<2|kc)&&rd<=4;++rd)
				//{
				//	if(abs(curr-(pred+(unsigned char)diff[idx-rd]))<tolerance)
				//	{
				//		refdist=rd;
				//		diff[idx]=diff[idx-rd];
				//		dst[idx]=pred+diff[idx];
				//	}
				//}
			}
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(noise);
	free(freq);
	free(dst);
	free(diff);
}

static void x2_ycocb_fwd(char *comp)
{
	comp[0]-=comp[1];
	comp[1]+=comp[0]>>1;
	comp[2]-=comp[1];
	comp[1]+=comp[2]>>1;
}
static void x2_ycocb_inv(char *comp)
{
	comp[1]-=comp[2]>>1;
	comp[2]+=comp[1];
	comp[1]-=comp[0]>>1;
	comp[0]+=comp[1];
}
static int x2_grad(int N, int W, int NW)
{
	char vmin, vmax;
	if(N<W)
		vmin=N, vmax=W;
	else
		vmin=W, vmax=N;
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	return pred;
}
void preproc_x2(char *src, int iw, int ih)//https://github.com/Equationist/xpng/blob/master/main.c
{
	long long MSE[3]={0};
	memset(lossygrad_rmse, 0, sizeof(lossygrad_rmse));
	memset(lossygrad_psnr, 0, sizeof(lossygrad_psnr));

	int res=iw*ih;
	unsigned char *noise=(unsigned char*)malloc((size_t)res<<2);
	int *pref=(int*)malloc(256*sizeof(int));
	char *dst=(char*)malloc((size_t)res<<2);
	//char *diff=(char*)malloc((size_t)res<<2);
	if(!noise||!pref||!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int kc=0;kc<3;++kc)//calc noise
	{
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				int n=0;
				long long x=0, x2=0;
				for(int ky2=ky-X_BLOCKSIZE/2;ky2<=ky+X_BLOCKSIZE/2;++ky2)
				{
					if((unsigned)ky2<(unsigned)ih)
					{
						for(int kx2=ky-X_BLOCKSIZE/2;kx2<=kx+X_BLOCKSIZE/2;++kx2)
						{
							if((unsigned)kx2<(unsigned)iw)
							{
								char px=src[(iw*ky2+kx2)<<2|kc];
								x+=px;
								x2+=px*px;
								++n;
							}
						}
					}
				}
				noise[(iw*ky+kx)<<2|kc]=(int)round(sqrt((double)(x2-x*x)/n));
			}
		}
	}

	for(int step=1;step<=128;step<<=1)//set preference levels
	{
		for(int k=-128;k<128;k+=step)
			pref[(unsigned char)k]=step;
	}
	
	memcpy(dst, src, (size_t)res<<2);//copy alpha

	const int clevel=6;//0~40, default 6
	for(int ky=0;ky<ih;++ky)
	{
		set_window_title("%d/%d", ky+1, ih);
		for(int kx=0;kx<iw;++kx)
		{
			int idx;
			char N[3]={0}, W[3]={0}, NW[3]={0};
			if(ky>0)
			{
				idx=(iw*(ky-1)+kx)<<2;
				N[0]=dst[idx];
				N[1]=dst[idx|1];
				N[2]=dst[idx|2];
				if(kx>0)
				{
					idx=(iw*(ky-1)+kx-1)<<2;
					NW[0]=dst[idx];
					NW[1]=dst[idx|1];
					NW[2]=dst[idx|2];
				}
			}
			if(kx>0)
			{
				idx=(iw*ky+kx-1)<<2;
				W[0]=dst[idx];
				W[1]=dst[idx|1];
				W[2]=dst[idx|2];
			}
			x2_ycocb_fwd(N);
			x2_ycocb_fwd(W);
			x2_ycocb_fwd(NW);
			int pred[]=
			{
				x2_grad(N[0], W[0], NW[0]),
				x2_grad(N[1], W[1], NW[1]),
				x2_grad(N[2], W[2], NW[2]),
			};
			idx=(iw*ky+kx)<<2;
			unsigned char curr[]=
			{
				src[idx],
				src[idx|1],
				src[idx|2],
			};
			//unsigned char c2[]=
			//{
			//	src[idx],
			//	src[idx|1],
			//	src[idx|2],
			//};
			//x2_ycocb_fwd((char*)c2);
			
			pred[0]+=128;
			pred[1]+=128;
			pred[2]+=128;
			curr[0]+=128;
			curr[1]+=128;
			curr[2]+=128;
			//c2[0]+=128;
			//c2[1]+=128;
			//c2[2]+=128;

			if(kx==(iw>>1)&&ky==(ih>>1))
				printf("");

			int ltarget[3], rtarget[3];
			for(int kc=0;kc<3;++kc)
			{
				int tolerance=clevel+noise[idx|kc]*clevel/(10+clevel);

				ltarget[kc]=curr[kc]-tolerance,
				rtarget[kc]=curr[kc]+tolerance;
				if(ltarget[kc]<  0)ltarget[kc]=  0;
				if(rtarget[kc]>255)rtarget[kc]=255;
			}

			//find best quantization within allowable tolerance
			unsigned char approx[]=
			{
				curr[0],
				curr[1],
				curr[2],
			};
			int bestscore=0;
			for(int kr=ltarget[0];kr<=rtarget[0];++kr)
			{
				for(int kg=ltarget[1];kg<=rtarget[1];++kg)
				{
					for(int kb=ltarget[2];kb<=rtarget[2];++kb)
					{
						char c2[]=
						{
							kr-128,
							kg-128,
							kb-128,
						};
						x2_ycocb_fwd(c2);
						unsigned char delta[]=
						{
							c2[0]-pred[0],
							c2[1]-pred[1],
							c2[2]-pred[2],
						};
						int score=pref[delta[0]]+pref[delta[1]]+pref[delta[2]];
						if(bestscore<score)
						{
							bestscore=score;
							approx[0]=kr;
							approx[1]=kg;
							approx[2]=kb;
						}
					}
				}
			}
			curr[0]-=128;
			curr[1]-=128;
			curr[2]-=128;
			approx[0]-=128;
			approx[1]-=128;
			approx[2]-=128;

			dst[idx]=approx[0];
			dst[idx|1]=approx[1];
			dst[idx|2]=approx[2];

			int error;
			error=curr[0]-approx[0], MSE[0]+=error*error;
			error=curr[1]-approx[1], MSE[1]+=error*error;
			error=curr[2]-approx[2], MSE[2]+=error*error;
		}
	}
	for(int kc=0;kc<3;++kc)
	{
		lossygrad_rmse[kc]=sqrt((double)MSE[kc]/res);//RMSE is the Euclidean distance (l2 norm)
		lossygrad_psnr[kc]=20*log10(255/lossygrad_rmse[kc]);
	}

	ArrayHandle str;
	STR_ALLOC(str, 0);
	str_append(&str, "RMSE RGB %14lf %14lf %14lf\n", lossygrad_rmse[0], lossygrad_rmse[1], lossygrad_rmse[2]);
	str_append(&str, "PSNR RGB %14lf %14lf %14lf\n", lossygrad_psnr[0], lossygrad_psnr[1], lossygrad_psnr[2]);
	copy_to_clipboard((char*)str->data, (int)str->count);
	array_free(&str);

	memcpy(src, dst, (size_t)res<<2);
	free(noise);
	free(pref);
	free(dst);
	//free(diff);
}

#if 0
typedef struct KalmanInfoStruct
{
	double R, H, Q, P, Uhat, K;
} KalmanInfo;
static void kalman_init(KalmanInfo *info)
{
	info->R=-0.0001;//noise covariance, higher values reject more noise
	info->H=1;
	info->Q=10;
	info->P=0;
	info->Uhat=0;
	info->K=0;
}
static void kalman_predict(KalmanInfo *info, double U)
{
	info->K=info->P*info->H/(info->H*info->P*info->H+info->R);//update kalman gain
	info->Uhat+=info->K*(U-info->H*info->Uhat);//update estimated
	info->P=(1-info->K*info->H)*info->P+info->Q;//update error covariance
}
void kalman_apply(char *src, int iw, int ih, int fwd)
{
	//predict from grad
#if 1
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);//copy alpha

	KalmanInfo kalman[3];
	kalman_init(kalman+0);
	kalman_init(kalman+1);
	kalman_init(kalman+2);
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?pixels[(iw*(ky+(Y))+kx+(X))<<2|kc]:0
				char
					N =LOAD( 0, -1),
					W =LOAD(-1,  0),
					NW=LOAD(-1, -1);
				char vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				int pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
				
				kalman_predict(kalman+kc, pred);
				pred=(int)round(CLAMP(vmin, kalman[kc].Uhat, vmax));
				//pred=(int)round(CLAMP(-128, kalman[kc].Uhat, 127));
				int idx=(iw*ky+kx)<<2|kc;
				if(fwd)
					dst[idx]=src[idx]-pred;
				else
					dst[idx]=src[idx]+pred;
#undef  LOAD
			}
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
#endif

	//predict from left
#if 0
	KalmanInfo c[3];
	kalman_init(c);
	kalman_init(c+1);
	kalman_init(c+2);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;

			char curr[3];
			char pred[]=
			{
				(char)round(CLAMP(-128, c[0].Uhat, 127)),
				(char)round(CLAMP(-128, c[1].Uhat, 127)),
				(char)round(CLAMP(-128, c[2].Uhat, 127)),
			};
			if(fwd)
			{
				curr[0]=src[idx+0];
				curr[1]=src[idx+1];
				curr[2]=src[idx+2];
				src[idx+0]-=pred[0];
				src[idx+1]-=pred[1];
				src[idx+2]-=pred[2];
			}
			else
			{
				src[idx+0]+=pred[0];
				src[idx+1]+=pred[1];
				src[idx+2]+=pred[2];
				curr[0]=src[idx+0];
				curr[1]=src[idx+1];
				curr[2]=src[idx+2];
			}
			kalman_predict(c+0, curr[0]);
			kalman_predict(c+1, curr[1]);
			kalman_predict(c+2, curr[2]);
		}
	}
#endif
}
#endif


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
	//int *hist=(int*)malloc(256*sizeof(int));
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
	int *hist=(int*)malloc(256*sizeof(int));
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
		
		set_window_title("it %d/100: %lf", ki+1, 1/best->loss);

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
	int *hist=(int*)malloc(256*sizeof(int));
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
				//eNW + eN + eNE
				int w=
					 (ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0)
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0);
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
			if(vmin>cT)vmin=cT;
			if(vmax<cT)vmax=cT;
			//if(vmin>cTR)vmin=cTR;
			//if(vmax<cTR)vmax=cTR;
			pred=CLAMP(vmin, pred, vmax);

			//if(kc==1&&kx==384&&ky==256)//
			//if(pred)
			//	printf("");

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
				if(kx+1<iw)//add current error to eNE, such that eN = eN0 + eW
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
			if(vmin>ctopright)vmin=ctopright;
			if(vmax<ctopright)vmax=ctopright;

			if(vmin>ctop)vmin=ctop;
			if(vmax<ctop)vmax=ctop;

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

//clamped gradient predictor, aka LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS
int  predict_grad(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char
		W=kx?buf[idx-bytestride]:0,
		N=ky?buf[idx-rowlen]:0,
		NW=kx&&ky?buf[idx-rowlen-bytestride]:0;

	int pred;

	char vmax, vmin;
	if(N<W)
		vmin=N, vmax=W;
	else
		vmin=W, vmax=N;

	if(NW>vmax)//choose steepest slope if both downward or upward
		pred=vmin;
	else if(NW<vmin)
		pred=vmax;
	else
		pred=N+W-NW;

	//char xdelta=top-topleft, ydelta=left-topleft;
	//if((xdelta>0)==(ydelta>0))
	//	pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
	//else
	//	pred=topleft+xdelta+ydelta;//average slope
	
	//pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
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


#define CTX_NPRED 4
static int median3(int a, int b, int c)
{
	if(a<b)
	{
		if(b<c)
			return b;
		return a<c?c:a;
	}
	if(a<c)
		return a;
	return b<c?c:b;
}
//typedef struct Pred3InfoStruct
//{
//	//char key[2];
//	int val, count;
//} Pred3Info;
static double sigmoid(double x)
{
	if(x<-6)
		return 0;
	if(x>6)
		return 1;
	return 1/(1+exp(-x));
}
void pred_ctx(char *src, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	int *hist=(int*)malloc(0x10000*sizeof(int));
	//int *temp=(int*)malloc(iw*(CTX_NPRED<<1)*sizeof(int));
	//Pred3Info *context=(Pred3Info*)malloc(0x10000*sizeof(Pred3Info));
	if(!dst||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);//copy alpha
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	for(int kc=0;kc<3;++kc)
	{
		memset(hist, 0, 0x10000*sizeof(int));
		//memset(temp, 0, iw*(CTX_NPRED<<1)*sizeof(int));
		//memset(context, 0, 0x10000*sizeof(Pred3Info));
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?pixels[(iw*(ky+(Y))+kx+(X))<<2|kc]:0)

				char
					NNW=LOAD(-1, -2),
					NW =LOAD(-1, -1),
					N  =LOAD( 0, -1),
					NE =LOAD( 1, -1),
					W  =LOAD(-1,  0);

				//softclamp grad
#if 0
				char vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				int range=vmax-vmin;
				double mid=(vmax+vmin)*0.5;
				int pred;
				if(range)
				{
					pred=vmin+(int)(range/(1+exp(-((double)(N+W-NW)-mid)*10/range)));
					//pred=mid+(int)(range*(sigmoid(((double)(N+W-NW)-mid)*10/range)-0.5));
					//pred=mid+(int)(range*0.5*tanh(((double)(N+W-NW)-mid)*6/range));
				}
				else
					pred=W;
				if(pred<-128||pred>127)
					pred=CLAMP(-128, pred, 127);
#endif

				//collinearity
#if 0
				int pred=N+W-NW;
				//if(abs(N-(NW+NE)/2)>2&&abs(NW-(W+NNW)/2)>2)
				{
					char vmin, vmax;
					if(N<W)
						vmin=N, vmax=W;
					else
						vmin=W, vmax=N;
					pred=CLAMP(vmin, pred, vmax);
				}
#endif

				//blur
#if 0
				int preds[CTX_NPRED]={0};
				long long num=0;
				int den=0;
				int pred=0;
				for(int kp=0;kp<CTX_NPRED;++kp)
				{
					int weight=0;
					if(kx>0)//W
					{
						preds[kp]+=temp[iw*kp+kx-1];
						weight+=temp[iw*(CTX_NPRED+kp)+kx-1];
					}
					if(ky>0)//N
					{
						preds[kp]+=temp[iw*kp+kx];
						weight+=temp[iw*(CTX_NPRED+kp)+kx];
					}
					int sh=kx>0&&ky>0;
					preds[kp]>>=sh;
					weight>>=sh;
					weight=0x1000000/(weight+1);
					num+=(long long)preds[kp]*weight;
					den+=weight;
				}
				if(den)
					pred=(int)((num+(den>>1))/den);
#endif

				//clamped grad
#if 1
				char vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				int pred0=N+W-NW;
				//int pred=NW+(N-NW)+(W-NW);//same

				//int pred=NW+((N-NW)+(NE-N))/2+(W-NW);
				//int pred=NW+(NE-NW)/2+W-NW;
				//int pred=W+(NE-NW)/2;

				pred0=CLAMP(vmin, pred0, vmax);

				pred0+=128;
				pred0&=0xFF;
				int *dist=hist+(pred0<<8), pred=0, vmax2=0;
				//int start=pred0-32, end=pred0+32;
				//if(start<0)
				//	start=0;
				//if(end>256)
				//	end=256;
				for(int k=0;k<256;++k)//pred = argmax(distribution)
				{
					if(vmax2<dist[k])
						vmax2=dist[k], pred=k;
				}
				if(vmax2)
					pred-=128;
				else
					pred=pred0-128;
				if(pred<vmin||pred>vmax)
					pred=CLAMP(vmin, pred, vmax);

				//int pred=CLAMP(vmin, pred0, vmax);
#endif
				//int pred=median3(N, W, N+W-NW);//equivalent to clamped gradient

				//context
#if 0
				int ctx=(unsigned char)N<<8|(unsigned char)W;
				//int ctx=(unsigned char)(N-NW)<<8|(unsigned char)(W-NW);
				//int ctx=((unsigned char)(W-NW)>>6)<<6|((unsigned char)(N-NW)>>3)<<3|((unsigned char)(NE-N)>>3);
				//int ctx=((unsigned char)NE>>5)<<9|((unsigned char)NW>>5)<<6|((unsigned char)W>>5)<<3|((unsigned char)N>>5);
				ctx&=0xFFFF;
				Pred3Info *p=context+ctx;
				int pred;
				if(p->count)
					pred=(p->val+(p->count>>1))/p->count;
				else
				{
					char vmin, vmax;
					if(N<W)
						vmin=N, vmax=W;
					else
						vmin=W, vmax=N;
					pred=N+W-NW;
					pred=CLAMP(vmin, pred, vmax);
				}
#endif

				//median
#if 0
				char nb[5]=
				{
					LOAD( 0, -1),//N
					LOAD(-1,  0),//W
					LOAD(-1, -1),//NW
					//LOAD( 1, -1),//NE
				};
				nb[_countof(nb)-1]=nb[0]+nb[1]-nb[2];
				//char vmin, vmax;
				//if(nb[0]<nb[1])
				//	vmin=nb[0], vmax=nb[1];
				//else
				//	vmin=nb[1], vmax=nb[0];
				//nb[_countof(nb)-1]=CLAMP(vmin, nb[_countof(nb)-1], vmax);
				for(int i=0;i<_countof(nb)-1;++i)
				{
					for(int j=i+1;j<_countof(nb);++j)
					{
						if(nb[j]<nb[i])
						{
							char temp=nb[i];
							nb[i]=nb[j];
							nb[j]=temp;
						}
					}
				}
				char pred=(nb[1]+nb[2])>>1;
#endif

				int idx=(iw*ky+kx)<<2|kc;

				//if(kx==384&&ky==256)
				//	printf("");

				//pred+=128;
				//pred>>=8;
				if(fwd)
					dst[idx]=src[idx]-pred;
				else
					dst[idx]=src[idx]+pred;

				int curr=pixels[idx]+128;
				++hist[pred0<<8|curr];

			//	pred0+=128;
			//	pred0&=0xFF;
			//	int curr=pixels[idx]+128;
			//	++hist[curr<<8|pred0];

			//	int curr=pixels[idx];
			//	//curr<<=8;
			//	for(int kp=0;kp<CTX_NPRED;++kp)
			//	{
			//		temp[iw*kp+kx]+=(curr-temp[iw*kp+kx])>>kp;
			//		temp[iw*(CTX_NPRED+kp)+kx]=abs(curr-pred);
			//	}

				//p->val+=pixels[idx];
				//++p->count;
				//if(p->count>=0x10000)
				//{
				//	p->count>>=1;
				//	p->val/=2;
				//}
#undef  LOAD
			}
		}
#if 0
		int vmax=0;
		for(int k=0;k<0x10000;++k)
		{
			if(vmax<hist[k])
				vmax=hist[k];
		}
		const char *fn=0;
		switch(kc)
		{
		case 0:fn="dump_c0.PNG";break;
		case 1:fn="dump_c1.PNG";break;
		case 2:fn="dump_c2.PNG";break;
		}
		for(int k=0;k<0x10000;++k)
			((unsigned char*)hist)[k]=hist[k]*0xFF/vmax;
		lodepng_encode_file(fn, (unsigned char*)hist, 256, 256, LCT_GREY, 8);
#endif
#if 0
		FILE *f=fopen(fn, "w");
		for(int ky=0;ky<256;++ky)
		{
			for(int kx=0;kx<256;++kx)
				fprintf(f, " %02X", hist[ky<<8|kx]*0xFF/vmax);
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
		fclose(f);
#endif
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
	free(hist);
	//free(temp);
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
#undef  LOAD
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


//	#define C2_WSUM
//	#define C2_DISABLE_LEARNING

#define LOAD(BUF, C, X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|C]:0
static int custom2_loadnb(const char *pixels, const char *errors, int iw, int ih, int kx, int ky, int kc, char *nb)
{
	int j=-1;
	
	nb[++j]=LOAD(pixels, 0, -2, -2);
	nb[++j]=LOAD(pixels, 0, -1, -2);
	nb[++j]=LOAD(pixels, 0,  0, -2);
	nb[++j]=LOAD(pixels, 0,  1, -2);
	nb[++j]=LOAD(pixels, 0,  2, -2);
	nb[++j]=LOAD(pixels, 0, -2, -1);
	nb[++j]=LOAD(pixels, 0, -1, -1);
	nb[++j]=LOAD(pixels, 0,  0, -1);//7: c0.top
	nb[++j]=LOAD(pixels, 0,  1, -1);
	nb[++j]=LOAD(pixels, 0,  2, -1);
	nb[++j]=LOAD(pixels, 0, -2,  0);
	nb[++j]=LOAD(pixels, 0, -1,  0);//11: c0.left
	nb[++j]=LOAD(pixels, 1, -2, -2);
	nb[++j]=LOAD(pixels, 1, -1, -2);
	nb[++j]=LOAD(pixels, 1,  0, -2);
	nb[++j]=LOAD(pixels, 1,  1, -2);
	nb[++j]=LOAD(pixels, 1,  2, -2);
	nb[++j]=LOAD(pixels, 1, -2, -1);
	nb[++j]=LOAD(pixels, 1, -1, -1);
	nb[++j]=LOAD(pixels, 1,  0, -1);//19: c1.top
	nb[++j]=LOAD(pixels, 1,  1, -1);
	nb[++j]=LOAD(pixels, 1,  2, -1);
	nb[++j]=LOAD(pixels, 1, -2,  0);
	nb[++j]=LOAD(pixels, 1, -1,  0);//23: c1.left
	nb[++j]=LOAD(pixels, 2, -2, -2);
	nb[++j]=LOAD(pixels, 2, -1, -2);
	nb[++j]=LOAD(pixels, 2,  0, -2);
	nb[++j]=LOAD(pixels, 2,  1, -2);
	nb[++j]=LOAD(pixels, 2,  2, -2);
	nb[++j]=LOAD(pixels, 2, -2, -1);
	nb[++j]=LOAD(pixels, 2, -1, -1);
	nb[++j]=LOAD(pixels, 2,  0, -1);//31: c2.top
	nb[++j]=LOAD(pixels, 2,  1, -1);
	nb[++j]=LOAD(pixels, 2,  2, -1);
	nb[++j]=LOAD(pixels, 2, -2,  0);
	nb[++j]=LOAD(pixels, 2, -1,  0);//35: c2.left
	nb[++j]=LOAD(errors, 0, -2, -2);
	nb[++j]=LOAD(errors, 0, -1, -2);
	nb[++j]=LOAD(errors, 0,  0, -2);
	nb[++j]=LOAD(errors, 0,  1, -2);
	nb[++j]=LOAD(errors, 0,  2, -2);
	nb[++j]=LOAD(errors, 0, -2, -1);
	nb[++j]=LOAD(errors, 0, -1, -1);
	nb[++j]=LOAD(errors, 0,  0, -1);
	nb[++j]=LOAD(errors, 0,  1, -1);
	nb[++j]=LOAD(errors, 0,  2, -1);
	nb[++j]=LOAD(errors, 0, -2,  0);
	nb[++j]=LOAD(errors, 0, -1,  0);
	nb[++j]=LOAD(errors, 1, -2, -2);
	nb[++j]=LOAD(errors, 1, -1, -2);
	nb[++j]=LOAD(errors, 1,  0, -2);
	nb[++j]=LOAD(errors, 1,  1, -2);
	nb[++j]=LOAD(errors, 1,  2, -2);
	nb[++j]=LOAD(errors, 1, -2, -1);
	nb[++j]=LOAD(errors, 1, -1, -1);
	nb[++j]=LOAD(errors, 1,  0, -1);
	nb[++j]=LOAD(errors, 1,  1, -1);
	nb[++j]=LOAD(errors, 1,  2, -1);
	nb[++j]=LOAD(errors, 1, -2,  0);
	nb[++j]=LOAD(errors, 1, -1,  0);
	nb[++j]=LOAD(errors, 2, -2, -2);
	nb[++j]=LOAD(errors, 2, -1, -2);
	nb[++j]=LOAD(errors, 2,  0, -2);
	nb[++j]=LOAD(errors, 2,  1, -2);
	nb[++j]=LOAD(errors, 2,  2, -2);
	nb[++j]=LOAD(errors, 2, -2, -1);
	nb[++j]=LOAD(errors, 2, -1, -1);
	nb[++j]=LOAD(errors, 2,  0, -1);
	nb[++j]=LOAD(errors, 2,  1, -1);
	nb[++j]=LOAD(errors, 2,  2, -1);
	nb[++j]=LOAD(errors, 2, -2,  0);
	nb[++j]=LOAD(errors, 2, -1,  0);

#if 0
#if C2_REACH>=2
	nb[++j]=LOAD(pixels, 0, -2,  0);
	nb[++j]=LOAD(pixels, 0, -2, -1);
	nb[++j]=LOAD(pixels, 0, -2, -2);
	nb[++j]=LOAD(pixels, 0, -1, -2);
	nb[++j]=LOAD(pixels, 0,  0, -2);
	nb[++j]=LOAD(pixels, 0,  1, -2);
	nb[++j]=LOAD(pixels, 0,  2, -2);
	nb[++j]=LOAD(pixels, 0,  2, -1);
	nb[++j]=LOAD(pixels, 1, -2,  0);
	nb[++j]=LOAD(pixels, 1, -2, -1);
	nb[++j]=LOAD(pixels, 1, -2, -2);
	nb[++j]=LOAD(pixels, 1, -1, -2);
	nb[++j]=LOAD(pixels, 1,  0, -2);
	nb[++j]=LOAD(pixels, 1,  1, -2);
	nb[++j]=LOAD(pixels, 1,  2, -2);
	nb[++j]=LOAD(pixels, 1,  2, -1);
	nb[++j]=LOAD(pixels, 2, -2,  0);
	nb[++j]=LOAD(pixels, 2, -2, -1);
	nb[++j]=LOAD(pixels, 2, -2, -2);
	nb[++j]=LOAD(pixels, 2, -1, -2);
	nb[++j]=LOAD(pixels, 2,  0, -2);
	nb[++j]=LOAD(pixels, 2,  1, -2);
	nb[++j]=LOAD(pixels, 2,  2, -2);
	nb[++j]=LOAD(pixels, 2,  2, -1);
	nb[++j]=LOAD(errors, 0, -2,  0);
	nb[++j]=LOAD(errors, 0, -2, -1);
	nb[++j]=LOAD(errors, 0, -2, -2);
	nb[++j]=LOAD(errors, 0, -1, -2);
	nb[++j]=LOAD(errors, 0,  0, -2);
	nb[++j]=LOAD(errors, 0,  1, -2);
	nb[++j]=LOAD(errors, 0,  2, -2);
	nb[++j]=LOAD(errors, 0,  2, -1);
	nb[++j]=LOAD(errors, 1, -2,  0);
	nb[++j]=LOAD(errors, 1, -2, -1);
	nb[++j]=LOAD(errors, 1, -2, -2);
	nb[++j]=LOAD(errors, 1, -1, -2);
	nb[++j]=LOAD(errors, 1,  0, -2);
	nb[++j]=LOAD(errors, 1,  1, -2);
	nb[++j]=LOAD(errors, 1,  2, -2);
	nb[++j]=LOAD(errors, 1,  2, -1);
	nb[++j]=LOAD(errors, 2, -2,  0);
	nb[++j]=LOAD(errors, 2, -2, -1);
	nb[++j]=LOAD(errors, 2, -2, -2);
	nb[++j]=LOAD(errors, 2, -1, -2);
	nb[++j]=LOAD(errors, 2,  0, -2);
	nb[++j]=LOAD(errors, 2,  1, -2);
	nb[++j]=LOAD(errors, 2,  2, -2);
	nb[++j]=LOAD(errors, 2,  2, -1);
#endif
	nb[++j]=LOAD(pixels, 0, -1, -1);
	nb[++j]=LOAD(pixels, 0,  0, -1);
	nb[++j]=LOAD(pixels, 0,  1, -1);
	nb[++j]=LOAD(pixels, 0, -1,  0);
	nb[++j]=LOAD(pixels, 1, -1, -1);
	nb[++j]=LOAD(pixels, 1,  0, -1);
	nb[++j]=LOAD(pixels, 1,  1, -1);
	nb[++j]=LOAD(pixels, 1, -1,  0);
	nb[++j]=LOAD(pixels, 2, -1, -1);
	nb[++j]=LOAD(pixels, 2,  0, -1);
	nb[++j]=LOAD(pixels, 2,  1, -1);
	nb[++j]=LOAD(pixels, 2, -1,  0);
	nb[++j]=LOAD(errors, 0, -1, -1);
	nb[++j]=LOAD(errors, 0,  0, -1);
	nb[++j]=LOAD(errors, 0,  1, -1);
	nb[++j]=LOAD(errors, 0, -1,  0);
	nb[++j]=LOAD(errors, 1, -1, -1);
	nb[++j]=LOAD(errors, 1,  0, -1);
	nb[++j]=LOAD(errors, 1,  1, -1);
	nb[++j]=LOAD(errors, 1, -1,  0);
	nb[++j]=LOAD(errors, 2, -1, -1);
	nb[++j]=LOAD(errors, 2,  0, -1);
	nb[++j]=LOAD(errors, 2,  1, -1);
	nb[++j]=LOAD(errors, 2, -1,  0);
#endif

	if(kc>=1)
	{
		nb[++j]=LOAD(pixels, 0, 0, 0);
		nb[++j]=LOAD(errors, 0, 0, 0);
		if(kc==2)
		{
			nb[++j]=LOAD(pixels, 1, 0, 0);
			nb[++j]=LOAD(errors, 1, 0, 0);
		}
	}
	return ++j;
}
static int custom2_eval_fwd(const char *nb, const short *params, int count, int *ret_sum)
{
	int pred=0;
#ifdef C2_WSUM
	int sum=0;
#endif
	for(int k=0;k<count;++k)
	{
		pred+=nb[k]*params[k];
#ifdef C2_WSUM
		sum+=params[k];
#endif
	}
#ifdef C2_WSUM
	pred=sum?pred/sum:0;//pred = (params . nb)/(sum params)
	if(ret_sum)
		*ret_sum=sum;
#else
	pred+=1<<13;
	pred>>=14;
	//pred=CLAMP(-128, pred, 127);
#endif
	return pred;
}
static char clamp2(char pred, char v1, char v2)
{
	char vmin, vmax;
	if(v1<v2)
		vmin=v1, vmax=v2;
	else
		vmin=v2, vmax=v1;
	pred=CLAMP(vmin, pred, vmax);
	return pred;
}
static void custom2_train(const char *pixels, const char *errors, int iw, int ih, int kx, int ky, short *params, int kc)
{
	char nb[C2_REACH*(C2_REACH+1)*2*6+4];
	int nparams=custom2_loadnb(pixels, errors, iw, ih, kx, ky, kc, nb);
	int curr=LOAD(pixels, kc, 0, 0);

	//fwd pass
#ifdef C2_WSUM
	int sum=0;
	int pred=custom2_eval_fwd(nb, params, nparams, &sum);
#else
	int pred=custom2_eval_fwd(nb, params, nparams, 0);
#endif
	//int loss=((curr<<14)-pred);//loss = 1/2 (curr-pred)^2
	//loss*=loss;
	//loss>>=1;

	//bwd pass
	const int lr=(int)(0.02*0x10000+0.5);
	//const int lr=(int)(36*0x10000+0.5);
	int dL_dpred=pred-curr;//dL_dpred = (curr-pred)(-1) = pred-curr
#ifdef C2_WSUM
	for(int k=0;k<nparams;++k)//dpred_dparam[i] = (nb[i] - pred)/(sum params)
	{
		int grad=(int)((long long)lr*dL_dpred*(nb[k]-pred)/sum);
		
		int newval=params[k]-grad;
		newval=CLAMP(1, newval, 0x7FFF);
		params[k]=newval;
	}
#else
	for(int k=0;k<nparams;++k)//dpred_dparam[i] = nb[i]
	{
		int grad=dL_dpred*nb[k];
		grad=(int)(((long long)lr*grad+0x8000)>>16);

		//if(grad)//
		//	printf("");

		int newval=params[k]-grad;
		newval=CLAMP(-0x7FFF, newval, 0x7FFF);
		params[k]=newval;
	}
#endif
}
static void custom2_train_v2(const char *pixels, const char *errors, int iw, int ih, int kx, int ky, short *params, int kc)
{
	if((unsigned)kx>=(unsigned)iw||(unsigned)ky>=(unsigned)ih)
		return;
	char nb[C2_REACH*(C2_REACH+1)*2*6+4];
	int nparams=custom2_loadnb(pixels, errors, iw, ih, kx, ky, kc, nb);

	int curr=LOAD(pixels, kc, 0, 0);
	
#ifdef C2_WSUM
	int sum=0;
	int pred=custom2_eval_fwd(nb, params, nparams, &sum);
#else
	int pred=custom2_eval_fwd(nb, params, nparams, 0);
#endif

	int loss0=curr-pred;
	for(int it=0;it<(nparams<<4);++it)
	{
		int bitidx=(unsigned)xoroshiro128_next()%(nparams<<4);

		//if((bitidx>>4)>=nparams)
		//	LOG_ERROR("");

		params[bitidx>>4]^=1<<(bitidx&7);

		pred=0;
		for(int k=0;k<nparams;++k)
		{
			pred+=nb[k]*params[k];
#ifdef C2_WSUM
			sum+=params[k];
#endif
		}
#ifdef C2_WSUM
		pred/=sum;
#else
		pred+=1<<13;
		pred>>=14;
		pred=CLAMP(-128, pred, 127);
#endif
		int loss2=curr-pred;

		if(abs(loss2)>abs(loss0))//revert
			params[bitidx>>4]^=1<<(bitidx&7);
	}
}
void pred_custom2_apply(char *src, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *buf=(char*)malloc((size_t)res<<2);
	//Custom2Params *params=(Custom2Params*)malloc(iw*sizeof(Custom2Params));
	Custom2Params params=
	{
		{
#if C2_REACH>=2
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
#endif
			-0x2000, 0x4000, -0x2000,
			 0x4000,
			0, 0, 0, 0,
			0, 0, 0, 0,

			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
		},
		{
#if C2_REACH>=2
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
#endif
			 0x2000, -0x4000,  0x2000,
			-0x4000,
			-0x2000, 0x4000, -0x2000,
			 0x4000,
			0, 0, 0, 0,

			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,

			0x4000, 0,
		},
		{
#if C2_REACH>=2
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
#endif
			 0x1000, -0x2000,  0x1000,
			-0x2000,
			 0x1000, -0x2000,  0x1000,
			-0x2000,
			-0x2000, 0x4000, -0x2000,
			 0x4000,

			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,

			0x2000, 0,
			0x2000, 0,
		},
	};
	if(!buf)
	{
		LOG_ERROR("Alllocation error");
		return;
	}
	XOROSHIRO128_RESET();
	//for(int k=0;k<iw;++k)
	//{
	//	short *p=(short*)(params+k);
	//	for(int k2=0;k2<sizeof(Custom2Params)/sizeof(short);++k2)
	//		p[k2]=(short)(xoroshiro128_next()&0xFFFF);
	//}
#ifdef C2_WSUM
	{
		short *p=(short*)&params;
		for(int k=0;k<sizeof(Custom2Params)/sizeof(short);++k)
		{
			short val=(short)(xoroshiro128_next()&0x1FF)+1;
			val=CLAMP(1, val, 0x7FFF);
			p[k]=val;
		}
	}
#else
	//{
	//	short *p=(short*)&params;
	//	for(int k=0;k<sizeof(Custom2Params)/sizeof(short);++k)
	//		p[k]=(short)(xoroshiro128_next()&0x1FF)-0x100;
	//}
	//memset(&params, 0, sizeof(params));
#endif

	for(int k=0;k<res;++k)//copy alpha
		buf[k<<2|3]=src[k<<2|3];
	
	const int niter=32, train_xmask=0;
	const int schedule[]=
	{
		-2,  0,
		-2, -1,
		-2, -2,
		-1, -2,
		 0, -2,
		 1, -2,
		 2, -2,
		 2, -1,

		 1, -1,//dx, dy
		 0, -1,
		-1, -1,
		-1,  0,
	};
	const char *pixels=fwd?src:buf, *errors=fwd?buf:src;
	int pred, idx;
	long long errorsum=0, predsum=0;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==(iw>>1)&&ky==(ih>>1))//
			if(kx==1&&ky==1)//
				printf("");

			char nb[C2_REACH*(C2_REACH+1)*2*6+4];
			int count=custom2_loadnb(pixels, errors, iw, ih, kx, ky, 0, nb);
			//Custom2Params *p=params+kx;
			
			idx=(iw*ky+kx)<<2;

#ifndef C2_DISABLE_LEARNING
			int xoffset, yoffset;
			if(kx)
				xoffset=-1, yoffset=0;
			else
				xoffset=0, yoffset=-1;
			if(!(kx&train_xmask))
			{
				//custom2_train_v2(pixels, errors, iw, ih, kx+xoffset, ky+yoffset, params.c0, 0);
				for(int k2=0;k2<niter;++k2)
				{
					for(int k=0;k<_countof(schedule);k+=2)
						custom2_train(pixels, errors, iw, ih, kx+schedule[k], ky+schedule[k+1], params.c0, 0);
				}
			}
#endif
			pred=custom2_eval_fwd(nb, params.c0, C2_REACH*(C2_REACH+1)*2*6, 0);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				buf[idx]=src[idx]-pred;
			else
				buf[idx]=src[idx]+pred;

			errorsum+=abs(errors[idx]);
			predsum+=pred;
			nb[count++]=pixels[idx];
			nb[count++]=errors[idx];
			++idx;

			
#ifndef C2_DISABLE_LEARNING
			if(!(kx&train_xmask))
			{
				//custom2_train_v2(pixels, errors, iw, ih, kx+xoffset, ky+yoffset, params.c1, 1);
				for(int k2=0;k2<niter;++k2)
				{
					for(int k=0;k<_countof(schedule);k+=2)
						custom2_train(pixels, errors, iw, ih, kx+schedule[k], ky+schedule[k+1], params.c1, 1);
				}
			}
#endif
			pred=custom2_eval_fwd(nb, params.c1, C2_REACH*(C2_REACH+1)*2*6+2, 0);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				buf[idx]=src[idx]-pred;
			else
				buf[idx]=src[idx]+pred;
			
			errorsum+=abs(errors[idx]);
			predsum+=pred;
			nb[count++]=pixels[idx];
			nb[count++]=errors[idx];
			++idx;

			
#ifndef C2_DISABLE_LEARNING
			if(!(kx&train_xmask))
			{
				//custom2_train_v2(pixels, errors, iw, ih, kx+xoffset, ky+yoffset, params.c2, 2);
				for(int k2=0;k2<niter;++k2)
				{
					for(int k=0;k<_countof(schedule);k+=2)
						custom2_train(pixels, errors, iw, ih, kx+schedule[k], ky+schedule[k+1], params.c2, 2);
				}
			}
#endif
			pred=custom2_eval_fwd(nb, params.c2, C2_REACH*(C2_REACH+1)*2*6+4, 0);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				buf[idx]=src[idx]-pred;
			else
				buf[idx]=src[idx]+pred;
			errorsum+=abs(errors[idx]);
			predsum+=pred;
		}
		set_window_title("%d/%d: %6.2lf%%  P%lf  E%lf", ky+1, ih, 100.*(ky+1)/ih, (double)predsum/((ky+1)*iw), (double)errorsum/((ky+1)*iw));
	}
	memcpy(src, buf, (size_t)res<<2);

	//free(params);
	free(buf);
}
#undef  LOAD


	#define C2_USE_GA
//	#define C2_VALIDATE

#ifdef C2_VALIDATE
#define C2_BUFSIZE 0x1000
char c2_pixel[C2_BUFSIZE]={0}, c2_error[C2_BUFSIZE]={0}, c2_pred[C2_BUFSIZE]={0};
int c2_idx=0, c2_debugstate=0;
void c2_start_record()
{
	c2_idx=0;
	c2_debugstate=1;
}
void c2_start_check()
{
	c2_idx=0;
	c2_debugstate=2;
}
void c2_checkpoint(char pixel, char error, char pred)
{
	if(c2_debugstate==1)
	{
		if(c2_idx<C2_BUFSIZE)
		{
			c2_pixel[c2_idx]=pixel;
			c2_error[c2_idx]=error;
			c2_pred[c2_idx]=pred;
			++c2_idx;
		}
	}
	else if(c2_debugstate==2)
	{
		if(c2_idx<C2_BUFSIZE)
		{
			if(pixel!=c2_pixel[c2_idx]||error!=c2_error[c2_idx]||pred!=c2_pred[c2_idx])
				LOG_ERROR("Validation error");
			++c2_idx;
		}
	}
}
#endif

#if 0
double kalman(double U)
{
	static const double
		R=40,//noise covariance (actually 10)
		H=1;//measurement map scalar
	static double
		Q=10,//initial estimated covariance
		P=0,//initial error covariance (must be 0)
		Uhat=0,//initial estimated state (assume we don't know)
		K=0;//initial Kalman gain

	//begin
	K=P*H/(H*P*H+R);//update kalman gain
	Uhat+=K*(U-H*Uhat);//update estimated
	P=(1-K*H)*P+Q;//update error covariance
	return Uhat;
}
#endif
#if 0
void c2_test(char *src, int count)
{
	static const double
		R=40,//noise covariance (actually 10)
		H=1;//measurement map scalar
	static double
		Q=10,//initial estimated covariance
		P=0,//initial error covariance (must be 0)
		Uhat=0,//initial estimated state (assume we don't know)
		K=0;//initial Kalman gain

	console_start();
	console_log("idx  src  pred  diff\n");
	for(int k=0;k<count;++k)
	{
		double U=src[k<<2];//red channel
		K=P*H/(H*P*H+R);//update kalman gain
		Uhat+=K*(U-H*Uhat);//update estimated
		P=(1-K*H)*P+Q;//update error covariance

		console_log("%4d %4d %11lf  %11lf\n", k, src[k<<2], Uhat, src[k<<2]-Uhat);
	}
	console_log("Done.\n");
	pause();
}
#endif

#define C2_NPARAMSTOTAL (C2_REACH*(C2_REACH+1)*2*6*3+6)
//unsigned char *debug_buf=0;
Custom2Params c2_params=
{
	{//for c0
		0, 0, 0, 0, 0,//pixels from c0
		0, -0x4000, 0x4000, 0, 0,
		0, 0x4000,
	
		0, 0, 0, 0, 0,//pixels from c1
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,//pixels from c2
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	},
	{//for c1
		0, 0, 0, 0, 0,//pixels from c0
		0, 0x4000, -0x4000, 0, 0,
		0, -0x4000,
	
		0, 0, 0, 0, 0,//pixels from c1
		0, -0x4000, 0x4000, 0, 0,
		0, 0x4000,
	
		0, 0, 0, 0, 0,//pixels from c2
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,

		0x4000, 0,
	},
	{//for c2
		0, 0, 0, 0, 0,//pixels from c0
		0, 0x4000, -0x4000, 0, 0,
		0, -0x4000,
	
		0, 0, 0, 0, 0,//pixels from c1
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,//pixels from c2
		0, -0x4000, 0x4000, 0, 0,
		0, 0x4000,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,
	
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0,

		0x4000, 0,
		0, 0,
	},
};
#ifdef C2_USE_GA
typedef short C2OptParam_t;
#else
typedef double C2OptParam_t;
#endif
typedef struct Custom2OptInfoStruct
{
	double loss, invCR[3];//loss is the average
	C2OptParam_t params[C2_NPARAMSTOTAL];
	//Custom2Params params;//params are double to calculate centroid and have smoother training
} Custom2OptInfo;
//void print_histogram(int *hist)
//{
//	console_start();
//	for(int k=0;k<256;++k)
//	{
//		console_log("%3d  0x%04X", k, hist[k]);
//	}
//}
static void custom2_prealloc(const char *src, int iw, int ih, int fwd, Custom2Params const *params, char *dst)
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int pred, idx;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==(iw>>1)&&ky==(ih>>1))//
			//if(kx==1&&ky==1)//
			//	printf("");

			char nb[C2_REACH*(C2_REACH+1)*2*6+4];
			int count=custom2_loadnb(pixels, errors, iw, ih, kx, ky, 0, nb);
			
			idx=(iw*ky+kx)<<2;
			
			pred=custom2_eval_fwd(nb, params->c0, C2_REACH*(C2_REACH+1)*2*6, 0);
			//pred=clamp2(pred, nb[7], nb[11]);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
#ifdef C2_VALIDATE
			c2_checkpoint(pixels[idx], errors[idx], pred);
#endif

			nb[count++]=pixels[idx];
			nb[count++]=errors[idx];
			++idx;

			
			pred=custom2_eval_fwd(nb, params->c1, C2_REACH*(C2_REACH+1)*2*6+2, 0);
			//pred=clamp2(pred, nb[19], nb[23]);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
#ifdef C2_VALIDATE
			c2_checkpoint(pixels[idx], errors[idx], pred);
#endif
			
			nb[count++]=pixels[idx];
			nb[count++]=errors[idx];
			++idx;

			
#ifdef C2_VALIDATE
			if((c2_debugstate==1||c2_debugstate==2)&&c2_idx==2306)//record
				printf("");
#endif
			pred=custom2_eval_fwd(nb, params->c2, C2_REACH*(C2_REACH+1)*2*6+4, 0);
			//pred=clamp2(pred, nb[31], nb[35]);
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
#ifdef C2_VALIDATE
			c2_checkpoint(pixels[idx], errors[idx], pred);
			++idx;
#endif
		}
	}
}
void custom2_apply(char *src, int iw, int ih, int fwd, Custom2Params const *params)
{
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	custom2_prealloc(src, iw, ih, fwd, params, temp);

	for(int k=0;k<res;++k)//copy alpha
		temp[k<<2|3]=src[k<<2|3];
	
	memcpy(src, temp, (size_t)res<<2);
	free(temp);
}
static void custom2_cvtparams(const C2OptParam_t *src, short *dst, int count)
{
#ifdef C2_USE_GA
	memcpy(dst, src, count*sizeof(short));
#else
	for(int k=0;k<count;++k)
		dst[k]=(short)CLAMP(-0x7FFF, src[k], 0x7FFF);
#endif
}
static void custom2_opt2_calcloss(const char *src, int iw, int ih, Custom2OptInfo *info, char *temp, int *hist)
{
	int res=iw*ih;
	short params[_countof(info->params)];
	custom2_cvtparams(info->params, params, _countof(info->params));
	custom2_prealloc(src, iw, ih, 1, (Custom2Params*)params, temp);

	memset(hist, 0, 768*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char *p=temp+(k<<2);
		++hist[p[0]];
		++hist[p[1]|256];
		++hist[p[2]|512];
	}
	info->loss=0;
	for(int kc=0;kc<3;++kc)
	{
		double entropy[3]={0};
		for(int sym=0;sym<256;++sym)//Shannon law
		{
			int freq=hist[kc<<8|sym];
			if(freq)
			{
				double prob=(double)freq/res;
				entropy[kc]-=prob*log2(prob);
			}
		}
		info->invCR[kc]=entropy[kc]/8;
		info->loss+=info->invCR[kc];
	}
	info->loss/=3;
}
void custom2_opt(const char *src, int iw, int ih, Custom2Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;
	++call_idx;
	//c2_test(src+iw*(ih>>1), 128);

#ifdef C2_USE_GA
	int res=iw*ih;
	unsigned char *temp=(unsigned char*)malloc((size_t)res<<2);
	int *hist=(int*)malloc(768*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int k=0;k<res;++k)//set alpha for preview
		temp[k<<2|3]=0xFF;
	Custom2OptInfo info;
#define CALC_LOSS() custom2_opt2_calcloss(src, iw, ih, &info, (char*)temp, hist)
	const int nd=C2_NPARAMSTOTAL;
	//const int nd[]=
	//{
	//	C2_REACH*(C2_REACH+1)*2*6,
	//	C2_REACH*(C2_REACH+1)*2*6+2,
	//	C2_REACH*(C2_REACH+1)*2*6+4,
	//};
	short *p0=(short*)srcparams;
	memcpy(info.params, srcparams, sizeof(info.params));
	if(loud)
		srand((unsigned)__rdtsc());//
	
	//random starting point
#if 1
	if(call_idx==1)
	{
		int sum=0;
		for(int k=0;k<nd;++k)
			sum+=info.params[k];
		if(!sum)
		{
			for(int k=0;k<nd;++k)
				//info.params[k]=(rand()&0x1FFF)-0x1000;
				info.params[k]=(rand()&0x1FF)-0x100;
				//info.params[k]=(rand()&0x1F)-0x10;
		}
	}
#endif

	CALC_LOSS();
	double invCR0[3];
	memcpy(invCR0, info.invCR, sizeof(invCR0));
	//const int niter=C2_NPARAMSTOTAL*10;
	if(!niter)
		niter=2220;
	for(int it=0;it<niter;++it)
	{
		if(loud)
			set_window_title("%d %4d/%4d: TRGB %lf  %lf %lf %lf  %d", loud, it+1, niter, 1/info.loss, 1/info.invCR[0], 1/info.invCR[1], 1/info.invCR[2], call_idx);
		//short delta[C2_NPARAMSTOTAL];//X
		//for(int k=0;k<_countof(delta);++k)
		//	delta[k]=(rand()%3)-1;
		int idx=it%nd;
		//int idx=it*nd/niter;
		//int idx=rand()%nd,
		int inc=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1));
		//int inc=(rand()&0x1FFF)-0x1000;
		//int inc=(rand()&0x1FF)-0x100;
		//int inc=(rand()&0x1F)-0x10;
		//int inc=((rand()&1)<<1)-1;
		//int inc=(int)(((rand()<<1)-0x8000)*pow(info.loss, 20));
		//for(int k=0;k<_countof(delta);++k)
		//	info.params[k]+=delta[k];
		info.params[idx]+=inc;
		CALC_LOSS();
		double loss0=(invCR0[0]+invCR0[1]+invCR0[2])/3;
		if(info.loss>loss0)
		{
			info.loss=loss0;
			memcpy(info.invCR, invCR0, sizeof(info.invCR));
			//for(int k=0;k<_countof(delta);++k)
			//	info.params[k]-=delta[k];
			info.params[idx]-=inc;
		}
		else
			memcpy(invCR0, info.invCR, sizeof(invCR0));

		//X
#if 0
		int idx[]=
		{
			rand()%nd[0],
			nd[0]+rand()%nd[1],
			nd[0]+nd[1]+rand()%nd[2],
		};
		int inc[]=
		{
			(rand()&0x1FFF)-0x1000,
			(rand()&0x1FFF)-0x1000,
			(rand()&0x1FFF)-0x1000,
		};
		info.params[idx[0]]+=inc[0];
		info.params[idx[1]]+=inc[1];
		info.params[idx[2]]+=inc[2];
		//int bitidx=rand()%(nd<<4);
		//info.params[bitidx>>4]^=1<<(bitidx&15);
		
#ifdef C2_VALIDATE
		if(it==niter-1)
			c2_start_record();
		//	printf("");
#endif
		CALC_LOSS();

		for(int k=0;k<3;++k)
		{
			if(info.invCR[k]>invCR0[k])
				info.invCR[k]=invCR0[k], info.params[idx[k]]-=inc[k];
			else
				invCR0[k]=info.invCR[k];
		}
		info.loss=(info.invCR[0]+info.invCR[1]+info.invCR[2])/3;
#endif

		//preview
#if 1
		if(loud)
		{
			memcpy(srcparams, info.params, sizeof(info.params));
			ch_cr[0]=(float)(1/info.invCR[0]);
			ch_cr[1]=(float)(1/info.invCR[1]);
			ch_cr[2]=(float)(1/info.invCR[2]);
			unsigned char *ptr;
			addhalf(temp, iw, ih, 3, 4);
			SWAPVAR(image, temp, ptr);
			io_render();
			SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	if(!loud)
		memcpy(srcparams, info.params, sizeof(info.params));
	if(loss)
		memcpy(loss, invCR0, sizeof(invCR0));
	
#ifdef C2_VALIDATE
	c2_start_check();
#endif
#if 0
	void *ptr=realloc(debug_buf, (size_t)res<<2);
	if(!ptr)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	debug_buf=(unsigned char*)ptr;
	memcpy(debug_buf, temp, (size_t)res<<2);
#endif
	free(temp);
	free(hist);
#else
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	int *hist=(int*)malloc(768*sizeof(int));
	const int
		nv=C2_NPARAMSTOTAL,
		np=C2_NPARAMSTOTAL+1;
	Custom2OptInfo *params=(Custom2OptInfo*)malloc((np+3LL)*sizeof(Custom2OptInfo));
	if(!temp||!hist||!params)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, (size_t)res<<2);
	Custom2OptInfo
		*best=params,
		*worst=params+np-1,
		*x0=params+np,
		*xr=x0+1,
		*x2=xr+1;
	
#define CALC_LOSS(X) custom2_opt2_calcloss(src, iw, ih, X, temp, hist)
	
	//initialize N+1 param sets
	srand((unsigned)__rdtsc());
	short *p0=(short*)srcparams;
	for(int kp=0;kp<np;++kp)
	{
		Custom2OptInfo *x=best+kp;
		//short *p=(short*)&x->params;
		//if(kp)
		//{
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=p0[k2]+(rand()&0x1FF)-0x100;
				//x->params[k2]=p0[k2]+rand()&0xFFFF;
		//}
		//else
		//	memcpy(x->params, customparam_st+O2_NPARAMS*kc, sizeof(x->params));
		CALC_LOSS(x);
		//x->loss=opt_cr2_calcloss(x->params, buf, iw, ih, kc, temp, hist);
		set_window_title("init %d/%d", kp+1, np);
	}
	for(int k2=0;k2<nv;++k2)
		x0->params[k2]=p0[k2];
	//memcpy(x0->params, srcparams, sizeof(x0->params));
	CALC_LOSS(x0);
	double loss0=x0->loss;
	const double alpha=1, gamma=2, rho=0.5, sigma=0.5;
	for(int ki=0;ki<100;++ki)
	{
		//1  order
		isort(best, np, sizeof(Custom2OptInfo), cmp_optcr2info);
		
		set_window_title("it %d/100: %lf", ki+1, 1/best->loss);

		//2  get the centroid of all points except worst
		memset(x0->params, 0, sizeof(x0->params));
		for(int k2=0;k2<nv;++k2)//exclude the worst point
		{
			Custom2OptInfo *x=best+k2;
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
			memcpy(worst, xr, sizeof(*worst));
			continue;
		}

		//4  expansion
		if(xr->loss<best->loss)//if xr is best so far
		{
			for(int k2=0;k2<nv;++k2)
				x2->params[k2]=x0->params[k2]+(xr->params[k2]-x0->params[k2])*gamma;
			CALC_LOSS(x2);
			if(x2->loss<xr->loss)
				memcpy(worst, x2, sizeof(*worst));
			else
				memcpy(worst, xr, sizeof(*worst));
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
				memcpy(worst, x2, sizeof(*worst));
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
				memcpy(worst, x2, sizeof(*worst));
				continue;
			}
		}

		//6  shrink
		for(int kp=1;kp<np;++kp)
		{
			Custom2OptInfo *x=best+kp;
			for(int k2=0;k2<nv;++k2)
				x->params[k2]=best->params[k2]+(x->params[k2]-best->params[k2])*sigma;
			CALC_LOSS(x2);
		}
	}
#undef CALC_LOSS
	if(best->loss<loss0)
		custom2_cvtparams(best->params, p0, _countof(best->params));
	//	memcpy(srcparams, best->params, sizeof(best->params));

	free(params);
	free(temp);
	free(hist);
#endif
}
void custom2_opt_blocks(const char *src, int iw, int ih, Custom2Params *srcparams)
{
	const int blocksize=128, niter=512;
	static int call_idx=0;
	++call_idx;

	char *buf=(char*)malloc((size_t)blocksize*blocksize<<2);
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	double invCR[4];
	for(int it=0;it<niter;++it)
	{
		int x=rand()%(iw-(blocksize-1)), y=rand()%(ih-(blocksize-1));
		for(int ky=0;ky<blocksize;++ky)
			memcpy(buf+(((size_t)blocksize*ky)<<2), src+((iw*((size_t)y+ky)+x)<<2), (size_t)blocksize<<2);

		custom2_opt(buf, blocksize, blocksize, srcparams, 100, 10, 0, invCR);

		invCR[3]=(invCR[0]+invCR[1]+invCR[2])/3;
		set_window_title("%4d/%4d: TRGB %lf  %lf %lf %lf  %d", it+1, niter, 1/invCR[3], 1/invCR[0], 1/invCR[1], 1/invCR[2], call_idx);
	}
	free(buf);
}
unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
void custom2_opt_batch(Custom2Params *srcparams)
{
	ArrayHandle folder=dialog_open_folder();
	if(!folder)
		return;
	const char *extensions[]=
	{
		"PNG",
		"JPG",
		"JPEG",
	};
	ArrayHandle filenames=get_filenames((char*)folder->data, extensions, _countof(extensions), 1);
	array_free(&folder);
	if(!filenames)
		return;
	for(int ks=0;ks<(int)filenames->count;++ks)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, ks);
		int iw=0, ih=0;
		char *buf=(char*)stbi_load(fn[0]->data, &iw, &ih, 0, 4);
		if(!buf)
			continue;
		addhalf((unsigned char*)buf, iw, ih, 3, 4);
		colortransform_ycmcb_fwd(buf, iw, ih);//
		custom2_opt(buf, iw, ih, srcparams, 0, 10, ks+1, 0);
		free(buf);
	}
	array_free(&filenames);
}


//CUSTOM3
//#define C3_GRIDSEARCH 1
//#define C3_OPT_UNISPHERE//X  bad
#ifdef C3_GRIDSEARCH
#define C3_OPT_NCOMP 3
#else
#define C3_OPT_NCOMP (C3_NPARAMS>>5)
#endif
Custom3Params c3_params={0};
static int custom3_loadnb(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-C3_REACH;ky2<0;++ky2)
	{
		for(int kx2=-C3_REACH;kx2<=C3_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-C3_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static void custom3_prealloc(const char *src, int iw, int ih, int fwd, Custom3Params const *params, char *dst)
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int pred, idx;
	const short *coeffs[]=
	{
		params->c00, params->c01, params->c02,
		params->c10, params->c11, params->c12,
		params->c20, params->c21, params->c22,
	};
	short nb[3][C3_NNB+2]={0};
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int count[3], idx2;
			for(int kc=0;kc<3;++kc)
				count[kc]=custom3_loadnb(pixels, errors, iw, ih, kc, kx, ky, nb[kc]);
			
			idx=(iw*ky+kx)<<2;
			idx2=0;
			
			for(int kdst=0;kdst<3;++kdst)
			{
				pred=0;
				for(int kc=0;kc<3;++kc)
					pred+=fast_dot(coeffs[idx2+kc], nb[kc], count[kc]);
				pred+=1<<13;
				pred>>=14;
				pred=CLAMP(-128, pred, 127);
				if(fwd)
					dst[idx]=src[idx]-pred;
				else
					dst[idx]=src[idx]+pred;

				nb[kdst][C3_NNB  ]=pixels[idx];
				nb[kdst][C3_NNB+1]=errors[idx];
				count[kdst]+=2;
				++idx;
				idx2+=3;
			}
		}
	}
}
void custom3_apply(char *src, int iw, int ih, int fwd, Custom3Params const *params)
{
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(temp, src, (size_t)res<<2);//copy alpha

	custom3_prealloc(src, iw, ih, fwd, params, temp);

	memcpy(src, temp, (size_t)res<<2);
	free(temp);
}
typedef struct Custom3OptInfoStruct
{
	double invCR[4];//loss == invCR[3]=(invCR[0]+invCR[1]+invCR[2])/3
	short params[C3_NPARAMS];
} Custom3OptInfo;
static void custom3_calcloss(const char *src, int iw, int ih, Custom3OptInfo *info, char *temp, int *hist)
{
	int res=iw*ih;
	custom3_prealloc(src, iw, ih, 1, (Custom3Params*)info->params, temp);

	memset(hist, 0, 768*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char *p=temp+(k<<2);
		++hist[p[0]];
		++hist[p[1]|256];
		++hist[p[2]|512];
	}
	info->invCR[3]=0;
	for(int kc=0;kc<3;++kc)
	{
		double entropy=0;
		for(int sym=0;sym<256;++sym)//Shannon's law
		{
			int freq=hist[kc<<8|sym];
			if(freq)
			{
				double prob=(double)freq/res;
				entropy-=prob*log2(prob);
			}
		}
		info->invCR[kc]=entropy/8;
		info->invCR[3]+=info->invCR[kc];
	}
	info->invCR[3]/=3;
}
static float gen_normal()
{
	float val=0;
	for(int k2=0;k2<12;++k2)
		val+=(float)rand()/RAND_MAX;
	val-=6;
	//val*=10;
	return val;
	//return (int)roundf(val);
}
void custom3_opt(const char *src, int iw, int ih, Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;
	++call_idx;

	int res=iw*ih;
	unsigned char *temp=(unsigned char*)malloc((size_t)res<<2);
	int *hist=(int*)malloc(768*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int k=0;k<res;++k)//set alpha for preview
		temp[k<<2|3]=0xFF;
	Custom3OptInfo info;
#define CALC_LOSS() custom3_calcloss(src, iw, ih, &info, (char*)temp, hist)
	const int nd=C3_NPARAMS;
	memcpy(info.params, srcparams, sizeof(info.params));
	if(loud)
		srand((unsigned)__rdtsc());//
	
	//random starting point		X  use ctrl N to add noise periodically
#if 0
	if(call_idx==1)
	{
		int sum=0;
		for(int k=0;k<nd;++k)
			sum+=info.params[k];
		if(!sum)//randomize if starting from all zeros
		{
			for(int k=0;k<nd;++k)
				//info.params[k]=(rand()&0x1FFF)-0x1000;
				//info.params[k]=(rand()&0x1FF)-0x100;
				info.params[k]=(rand()&0x1F)-0x10;
		}
	}
#endif

	CALC_LOSS();
	double invCR[4], loss0;
	memcpy(invCR, info.invCR, sizeof(invCR));
	loss0=info.invCR[3];
	if(!niter)
//#if C3_REACH==1
//		niter=C3_NPARAMS*100;
//#else
		niter=C3_NPARAMS*10;
//#endif
	//if(niter>1000)
	//	niter=1000;
	int shakethreshold=C3_NPARAMS;
	for(int it=0, watchdog=0;it<niter;++it)
	{
		if(loud)
			set_window_title("%d-%d %4d/%4d,%d/%d: %lf RGB %lf %lf %lf", loud, call_idx, it+1, niter, watchdog, shakethreshold, 1/info.invCR[3], 1/info.invCR[0], 1/info.invCR[1], 1/info.invCR[2]);
		int idx[C3_OPT_NCOMP]={0}, inc[C3_OPT_NCOMP]={0}, stuck=0;
		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(info.params, srcparams, sizeof(info.params));
			for(int k=0;k<C3_NPARAMS;++k)
				//info.params[k]+=gen_normal();//X  bump is too hard
				info.params[k]+=((rand()&1)<<1)-1;
				//info.params[k]+=rand()%3-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
#ifdef C3_OPT_UNISPHERE
			//const float power=1.f/C3_OPT_NCOMP;
			float
				fdelta[C3_OPT_NCOMP],
				radius=powf((float)rand()*(1<<maskbits)/RAND_MAX, 1.f/C3_OPT_NCOMP),
				//radius=4096,
				norm=0;
			for(int k=0;k<C3_OPT_NCOMP;++k)
			{
				fdelta[k]=gen_normal()*0.25f;
				norm+=fdelta[k]*fdelta[k];
			}
			radius*=C3_OPT_NCOMP/norm;
			for(int k=0;k<C3_OPT_NCOMP;++k)
			{
				fdelta[k]*=radius;
				fdelta[k]=CLAMP(-0x7FFF, fdelta[k], 0x7FFF);
			}
#endif
			for(int k=0;k<C3_OPT_NCOMP;++k)
			{
				idx[k]=rand()%nd;
				//idx[k]=(it+k)%nd;
#ifdef C3_GRIDSEARCH
				inc[k]=1;
				//inc[k]=((rand()&15)<<1)+1;
#else
#ifdef C3_OPT_UNISPHERE
				inc[k]=(int)roundf(fdelta[k]);
				//inc[k]=gen_normal();
#else
				while(!(inc[k]=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1))));//reject zero delta
#endif

				//int inc=(rand()&0x1FFF)-0x1000;
				//int inc=(rand()&0x1FF)-0x100;
				//int inc=(rand()&0x1F)-0x10;
				//int inc=((rand()&1)<<1)-1;
				//int inc=(int)(((rand()<<1)-0x8000)*pow(info.loss, 20));
		
				info.params[idx[k]]+=inc[k];
#endif
			}
		}

#ifdef C3_GRIDSEARCH
		short
			pz=info.params[idx[2]],
			py=info.params[idx[1]],
			px=info.params[idx[0]];
		short bestx=0, besty=0, bestz=0;
		double bestloss=0;
		int it2=0;
		for(int kz=-C3_GRIDSEARCH*inc[2];kz<=C3_GRIDSEARCH*inc[2];kz+=inc[2])
		{
			info.params[idx[2]]=pz+kz;
			for(int ky=-C3_GRIDSEARCH*inc[1];ky<=C3_GRIDSEARCH*inc[1];ky+=inc[1])
			{
				info.params[idx[1]]=py+ky;
				for(int kx=-C3_GRIDSEARCH*inc[0];kx<=C3_GRIDSEARCH*inc[0];kx+=inc[0], ++it2)
				{
					info.params[idx[0]]=px+kx;
					CALC_LOSS();
					if(!it2||bestloss>info.invCR[3])
					{
						bestloss=info.invCR[3];
						//memcpy(invCR, info.invCR, sizeof(info.invCR));
						bestx=info.params[idx[0]];
						besty=info.params[idx[1]];
						bestz=info.params[idx[2]];
					}
				}
			}
		}
		info.params[idx[0]]=bestx;
		info.params[idx[1]]=besty;
		info.params[idx[2]]=bestz;
#else
		CALC_LOSS();
#endif
		
		if(info.invCR[3]>invCR[3])//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				memcpy(invCR, info.invCR, sizeof(info.invCR));
			else
			{
				memcpy(info.invCR, invCR, sizeof(info.invCR));
				for(int k=0;k<C3_OPT_NCOMP;++k)
					info.params[idx[k]]-=inc[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss0>info.invCR[3])
			{
				memcpy(srcparams, info.params, sizeof(info.params));
				loss0=info.invCR[3];
				--it;//bis
			}
			memcpy(invCR, info.invCR, sizeof(invCR));
			watchdog=0;
		}

		//preview
#if 1
		if(loud)
		{
			ch_cr[0]=(float)(1/info.invCR[0]);
			ch_cr[1]=(float)(1/info.invCR[1]);
			ch_cr[2]=(float)(1/info.invCR[2]);
			unsigned char *ptr;
			addhalf(temp, iw, ih, 3, 4);
			SWAPVAR(image, temp, ptr);
			io_render();
			SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	//if(!loud)
	//	memcpy(srcparams, info.params, sizeof(info.params));
	if(loss)
		memcpy(loss, invCR, sizeof(invCR));
	
	free(temp);
	free(hist);
}
void custom3_opt_batch(Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	ArrayHandle folder=dialog_open_folder();
	if(!folder)
		return;
	const char *extensions[]=
	{
		"PNG",
		"JPG",
		"JPEG",
	};
	ArrayHandle filenames=get_filenames((char*)folder->data, extensions, _countof(extensions), 1);
	array_free(&folder);
	if(!filenames)
		return;

	int iw=0, ih=0;
	ArrayHandle bmp;
	ARRAY_ALLOC(int, bmp, 0, 0, 0, 0);
	for(int ks=0;ks<(int)filenames->count;++ks)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, ks);
		int iw2=0, ih2=0;
		char *buf=(char*)stbi_load(fn[0]->data, &iw2, &ih2, 0, 4);
		if(!buf)
			continue;
		if(!iw||iw==iw2)
		{
			iw=iw2;
			ih+=ih2;
			ARRAY_APPEND(bmp, buf, iw2*ih2, 1, 0);
		}
		free(buf);
	}
	array_free(&filenames);
	if(bmp->count)
	{
		addhalf(bmp->data, iw, ih, 3, 4);
		colortransform_ycmcb_fwd((char*)bmp->data, iw, ih);//
		custom3_opt((char*)bmp->data, iw, ih, srcparams, niter, maskbits, loud, 0);
	}
	array_free(&bmp);
}
typedef struct ImageStruct
{
	int iw, ih;
	char *data;
} Image;
void free_image(void *p)
{
	Image *im=(Image*)p;
	free(im->data);
	im->data=0;
}
void custom3_opt_batch2(Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	ArrayHandle folder=dialog_open_folder();
	if(!folder)
		return;
	const char *extensions[]=
	{
		"PNG",
		"JPG",
		"JPEG",
	};
	ArrayHandle filenames=get_filenames((char*)folder->data, extensions, _countof(extensions), 1);
	array_free(&folder);
	if(!filenames)
		return;

	ArrayHandle images;
	ARRAY_ALLOC(Image, images, 0, 0, filenames->count, free_image);
	for(int ks=0;ks<(int)filenames->count;++ks)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, ks);
		Image im2;
		im2.data=(char*)stbi_load(fn[0]->data, &im2.iw, &im2.ih, 0, 4);
		if(!im2.data)
			continue;
		addhalf(im2.data, im2.iw, im2.ih, 3, 4);
		colortransform_ycmcb_fwd((char*)im2.data, im2.iw, im2.ih);
		ARRAY_APPEND(images, &im2, 1, 1, 0);
	}
	array_free(&filenames);
	for(int it=0;it<10;++it)
	{
		int idx=rand()%images->count;
		Image *im=(Image*)array_at(&images, idx);
		custom3_opt(im->data, im->iw, im->ih, srcparams, C3_NPARAMS, maskbits, idx+1, 0);
	}
	array_free(&images);
}




//CUSTOM4
#define C4_OPT_NCOMP 4
Custom4Params c4_params={0};
static int custom4_loadnb(const char *pixels, const char *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-C4_REACH;ky2<0;++ky2)
	{
		for(int kx2=-C4_REACH;kx2<=C4_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=pixels[idx2];
				nb[++idx]=errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-C4_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=pixels[idx2];
			nb[++idx]=errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static void custom4_prealloc(const char *src, int iw, int ih, int fwd, Custom4Params const *params, char *dst)
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int pred, idx;
	short nb[3][C3_NNB+2]={0};
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int count[3];
			//int idx2;
			for(int kc=0;kc<3;++kc)
				count[kc]=custom4_loadnb(pixels, errors, iw, ih, kc, kx, ky, nb[kc]);
			
			idx=(iw*ky+kx)<<2;
			
			pred=0;
			for(int kh=0;kh<C4_HIDDENNODES;++kh)
			{
				int pred2=0;
				pred2+=fast_dot(params->c00[kh]         , nb[0], C4_NNB);
				pred2+=fast_dot(params->c00[kh]+C4_NNB  , nb[1], C4_NNB);
				pred2+=fast_dot(params->c00[kh]+C4_NNB*2, nb[2], C4_NNB);
				pred2+=1<<13;
				pred2>>=14;
				pred+=(int)(params->c01[kh]*tanh(pred2)*(128/16384.));
				//pred2=(int)(128*tanh(pred2));
				//pred+=(int)((long long)params->c01[kh]*pred2>>14);
			}
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
			nb[0][C3_NNB  ]=pixels[idx];
			nb[0][C3_NNB+1]=errors[idx];
			++idx;
			
			pred=0;
			for(int kh=0;kh<C4_HIDDENNODES;++kh)
			{
				int pred2=0;
				pred2+=fast_dot(params->c10[kh]           , nb[0], C4_NNB+2);
				pred2+=fast_dot(params->c10[kh]+C4_NNB  +2, nb[1], C4_NNB);
				pred2+=fast_dot(params->c10[kh]+C4_NNB*2+2, nb[2], C4_NNB);
				pred2+=1<<13;
				pred2>>=14;
				pred+=(int)(params->c11[kh]*tanh(pred2)*(128/16384.));
				//pred2=(int)(128*tanh(pred2));
				//pred+=(int)((long long)params->c01[kh]*pred2>>14);
			}
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
			nb[1][C3_NNB  ]=pixels[idx];
			nb[1][C3_NNB+1]=errors[idx];
			++idx;
			
			pred=0;
			for(int kh=0;kh<C4_HIDDENNODES;++kh)
			{
				int pred2=0;
				pred2+=fast_dot(params->c20[kh]           , nb[0], C4_NNB+2);
				pred2+=fast_dot(params->c20[kh]+C4_NNB  +2, nb[1], C4_NNB+2);
				pred2+=fast_dot(params->c20[kh]+C4_NNB*2+4, nb[2], C4_NNB);
				pred2+=1<<13;
				pred2>>=14;
				pred+=(int)(params->c21[kh]*tanh(pred2)*(128/16384.));
				//pred2=(int)(128*tanh(pred2));
				//pred+=(int)((long long)params->c01[kh]*pred2>>14);
			}
			pred=CLAMP(-128, pred, 127);
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
		}
	}
}
void custom4_apply(char *src, int iw, int ih, int fwd, Custom4Params const *params)
{
	int res=iw*ih;
	char *temp=(char*)malloc((size_t)res<<2);
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(temp, src, (size_t)res<<2);//copy alpha

	custom4_prealloc(src, iw, ih, fwd, params, temp);

	memcpy(src, temp, (size_t)res<<2);
	free(temp);
}
typedef struct Custom4OptInfoStruct
{
	double invCR[4];//loss == invCR[3]=(invCR[0]+invCR[1]+invCR[2])/3
	short params[C4_NPARAMS];
} Custom4OptInfo;
static void custom4_calcloss(const char *src, int iw, int ih, Custom4OptInfo *info, char *temp, int *hist)
{
	int res=iw*ih;
	custom4_prealloc(src, iw, ih, 1, (Custom4Params*)info->params, temp);

	memset(hist, 0, 768*sizeof(int));
	for(int k=0;k<res;++k)
	{
		unsigned char *p=temp+(k<<2);
		++hist[p[0]];
		++hist[p[1]|256];
		++hist[p[2]|512];
	}
	info->invCR[3]=0;
	for(int kc=0;kc<3;++kc)
	{
		double entropy=0;
		for(int sym=0;sym<256;++sym)//Shannon's law
		{
			int freq=hist[kc<<8|sym];
			if(freq)
			{
				double prob=(double)freq/res;
				entropy-=prob*log2(prob);
			}
		}
		info->invCR[kc]=entropy/8;
		info->invCR[3]+=info->invCR[kc];
	}
	info->invCR[3]/=3;
}
void custom4_opt(const char *src, int iw, int ih, Custom4Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;
	++call_idx;

	int res=iw*ih;
	unsigned char *temp=(unsigned char*)malloc((size_t)res<<2);
	int *hist=(int*)malloc(768*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	for(int k=0;k<res;++k)//set alpha for preview
		temp[k<<2|3]=0xFF;
	Custom4OptInfo info;
#define CALC_LOSS() custom4_calcloss(src, iw, ih, &info, (char*)temp, hist)
	const int nd=C4_NPARAMS;
	memcpy(info.params, srcparams, sizeof(info.params));
	if(loud)
		srand((unsigned)__rdtsc());//

	CALC_LOSS();
	double invCR[4], loss0;
	memcpy(invCR, info.invCR, sizeof(invCR));
	loss0=info.invCR[3];
	if(!niter)
		niter=C4_NPARAMS*10;
	if(niter>1000)
		niter=1000;
	int shakethreshold=C3_NPARAMS;
	for(int it=0, watchdog=0;it<niter;++it)
	{
		if(loud)
			set_window_title("%d-%d %4d/%4d,%d/%d: %lf RGB %lf %lf %lf", loud, call_idx, it+1, niter, watchdog, shakethreshold, 1/info.invCR[3], 1/info.invCR[0], 1/info.invCR[1], 1/info.invCR[2]);
		int idx[C4_OPT_NCOMP]={0}, inc[C4_OPT_NCOMP]={0}, stuck=0;
		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(info.params, srcparams, sizeof(info.params));
			for(int k=0;k<C4_NPARAMS;++k)
				info.params[k]+=((rand()&1)<<1)-1;
				//info.params[k]+=rand()%3-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<C4_OPT_NCOMP;++k)
			{
				idx[k]=rand()%nd;
				//idx[k]=(it+k)%nd;

				while(!(inc[k]=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1))));//reject zero delta

				//int inc=(rand()&0x1FFF)-0x1000;
				//int inc=(rand()&0x1FF)-0x100;
				//int inc=(rand()&0x1F)-0x10;
				//int inc=((rand()&1)<<1)-1;
				//int inc=(int)(((rand()<<1)-0x8000)*pow(info.loss, 20));
		
				info.params[idx[k]]+=inc[k];
			}
		}

		CALC_LOSS();
		
		if(info.invCR[3]>invCR[3])//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				memcpy(invCR, info.invCR, sizeof(info.invCR));
			else
			{
				memcpy(info.invCR, invCR, sizeof(info.invCR));
				for(int k=0;k<C4_OPT_NCOMP;++k)
					info.params[idx[k]]-=inc[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss0>info.invCR[3])
			{
				memcpy(srcparams, info.params, sizeof(info.params));
				loss0=info.invCR[3];
				--it;//bis
			}
			memcpy(invCR, info.invCR, sizeof(invCR));
			watchdog=0;
		}

		//preview
#if 1
		if(loud)
		{
			ch_cr[0]=(float)(1/info.invCR[0]);
			ch_cr[1]=(float)(1/info.invCR[1]);
			ch_cr[2]=(float)(1/info.invCR[2]);
			unsigned char *ptr;
			addhalf(temp, iw, ih, 3, 4);
			SWAPVAR(image, temp, ptr);
			io_render();
			SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	//if(!loud)
	//	memcpy(srcparams, info.params, sizeof(info.params));
	if(loss)
		memcpy(loss, invCR, sizeof(invCR));
	
	free(temp);
	free(hist);
}


void pred_calic(char *buf, int iw, int ih, int fwd)//https://github.com/siddharths2710/CALIC/blob/master/calic.ipynb
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	int *arrN=(int*)malloc(1024*sizeof(int));
	int *arrS=(int*)malloc(1024*sizeof(int));
	if(!b2||!arrN||!arrS)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, buf, (size_t)res<<2);//copy alpha
	const int thresholds[]={5, 15, 25, 42, 60, 85, 140};
	const char *pixels=fwd?buf:b2, *errors=fwd?b2:buf;
	for(int kc=0;kc<3;++kc)//process each channel separately
	{
		int prev_error=0;
		memset(arrN, 0, 1024*sizeof(int));
		memset(arrS, 0, 1024*sizeof(int));
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|kc]:0
				char
					NNWW=LOAD(pixels, -2, -2),
					NNW =LOAD(pixels, -1, -2),
					NN  =LOAD(pixels,  0, -2),
					NNE =LOAD(pixels,  1, -2),
					NNEE=LOAD(pixels,  2, -2),
					NWW =LOAD(pixels, -2, -1),
					NW  =LOAD(pixels, -1, -1),
					N   =LOAD(pixels,  0, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					W   =LOAD(pixels, -1,  0);
#undef  LOAD
				//NNWW NNW NN NNE NNEE
				//NWW  NW  N  NE
				//WW   W   ?
				int dx=abs(W-WW)+abs(N-NW)+abs(N-NE);
				int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
				int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
				int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
				int pred;

				//'A context-based adaptive lossless/nearly-lossless coding scheme for continuous-tone images'
#if 0
				if(dy+dx>32)//sharp edge
					pred=(dy*W+dx*N)/(dy+dx)+(NE-NW)/8;
				else if(dy-dx>12)//horizontal edge
					pred=(2*W+N)/3+(NE-NW)/8;
				else if(dy-dx<-12)//vertical edge
					pred=(W+2*N)/3+(NE-NW)/8;
				else//shallow area
					pred=(W+N)/2+(NE-NW)/8;

				if(d45-d135>32)//sharp 135-deg diagonal edge
					pred+=(NE-NW)/8;
				else if(d45-d135>16)//135-deg diagonal edge
					pred+=(NE-NW)/16;
				else if(d45-d135<-32)//sharp 45-deg diagonal edge
					pred-=(NE-NW)/8;
				else if(d45-d135<-16)//45-deg diagonal edge
					pred-=(NE-NW)/16;
#endif

				//'CALIC - A context-based adaptive lossless image codec'
#if 1
				if(dy-dx>80)
					pred=W;
				else if(dy-dx<-80)
					pred=N;
				else
				{
					//pred=((W+N)*5+(NE+NW)*3)/16;
					//pred=((W+N)*3+(NE+NW))/8;
					pred=(W+N)/2+(NE-NW)/4;
					if(dy-dx>32)
						pred=(pred+W)/2;
					else if(dy-dx>8)
						pred=(pred*3+W)/4;
					else if(dy-dx<-32)
						pred=(pred+N)/2;
					else if(dy-dx<-8)
						pred=(pred*3+N)/4;
				}
#endif
#if 0
				int vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
#endif

				int B=((2*W-WW)<pred)<<7|((2*N-NN)<pred)<<6|(WW<pred)<<5|(NN<pred)<<4|(NE<pred)<<3|(NW<pred)<<2|(W<pred)<<1|(N<pred);

				int delta=dx+dy+2*abs(prev_error);
				int Qdelta=7;
				for(int k=0;k<_countof(thresholds);++k)
				{
					if(delta<=thresholds[k])
					{
						Qdelta=k;
						break;
					}
				}
				int C=B*Qdelta/2;//context (texture identifier)
				//int C=B<<2|Qdelta>>1;//why not like this?

				++arrN[C];
				arrS[C]+=prev_error;
				if(arrN[C]>=255)
				{
					arrN[C]/=2;
					arrS[C]/=2;
				}
				int pred2=pred+arrS[C]/arrN[C];//sum of encountered errors / number of occurrences

				//if(kx==(iw>>1)&&ky==(ih>>1))//
				//	printf("");//

				int idx=(iw*ky+kx)<<2|kc;
				if(fwd)
					b2[idx]=buf[idx]-pred2;
				else
					b2[idx]=buf[idx]+pred2;

				//context[idx]=C;//store context for entropy coder

				prev_error=errors[idx];
			}
		}
	}
	memcpy(buf, b2, (size_t)res<<2);
	free(b2);
	free(arrN);
	free(arrS);
}


//	#define G2_WAVELET
	#define G2_MM

#ifdef G2_WAVELET
#define G2_REACH 8
#define G2_NNB (G2_REACH*(G2_REACH+1)*2)
static const int g2_sums[]=
{
	0x0001122A,
	0x0002A440,
	0x00043654,
	0x0005C47E,
	0x0002A440,
	0x0005C858,
	0x0008EC6A,
	0x000C08AA,
	0x00043654,
	0x0008EC6A,
	0x000DA278,
	0x00124CCC,
	0x0005C47E,
	0x000C08AA,
	0x00124CCC,
	0x00188152,
};
static const short g2_kernels[]=
{
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x0008, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0015, 0x01B9, 0x04B0, 0x01B9, 0x0015, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x01B9, 0x22A5, 0x5E2D, 0x22A5, 0x01B9, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0008, 0x04B0, 0x5E2D,

	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x0006, 0x0008, 0x0006, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0002, 0x0015, 0x007E, 0x01B9, 0x03A6, 0x04B0, 0x03A6, 0x01B9, 0x007E, 0x0015, 0x0002, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0002, 0x002E, 0x01B9, 0x09ED, 0x22A5, 0x4958, 0x5E2D, 0x4958, 0x22A5, 0x09ED, 0x01B9, 0x002E, 0x0002, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0008, 0x007E, 0x04B0, 0x1AFB, 0x5E2D, 0xC75F,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0001, 0x0002, 0x0005, 0x0007, 0x0008, 0x0007, 0x0005, 0x0002, 0x0001, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0005, 0x0015, 0x004A, 0x00CA, 0x01B9, 0x0301, 0x0432, 0x04B0, 0x0432, 0x0301, 0x01B9, 0x00CA, 0x004A, 0x0015, 0x0005, 0x0000,
	 0x0013, 0x0068, 0x01B9, 0x05DB, 0x0FEA, 0x22A5, 0x3C62, 0x5445, 0x5E2D, 0x5445, 0x3C62, 0x22A5, 0x0FEA, 0x05DB, 0x01B9, 0x0068, 0x0013,
	 0x0035, 0x011B, 0x04B0, 0x0FEA, 0x2B44, 0x5E2D, 0xA424, 0xE514,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0001, 0x0002, 0x0004, 0x0006, 0x0007, 0x0008, 0x0007, 0x0006, 0x0004, 0x0002, 0x0001, 0x0000, 0x0000, 0x0000,
	 0x0015, 0x0038, 0x007E, 0x00FB, 0x01B9, 0x02AB, 0x03A6, 0x0467, 0x04B0, 0x0467, 0x03A6, 0x02AB, 0x01B9, 0x00FB, 0x007E, 0x0038, 0x0015,
	 0x01B9, 0x0467, 0x09ED, 0x13BD, 0x22A5, 0x35A9, 0x4958, 0x5878, 0x5E2D, 0x5878, 0x4958, 0x35A9, 0x22A5, 0x13BD, 0x09ED, 0x0467, 0x01B9,
	 0x04B0, 0x0BF9, 0x1AFB, 0x35A9, 0x5E2D, 0x91DD, 0xC75F, 0xF07D,

	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x0008, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x002E, 0x007E, 0x002E, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0015, 0x01B9, 0x04B0, 0x01B9, 0x0015, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x007E, 0x09ED, 0x1AFB, 0x09ED, 0x007E, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x01B9, 0x22A5, 0x5E2D, 0x22A5, 0x01B9, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0006, 0x03A6, 0x4958, 0xC75F, 0x4958, 0x03A6, 0x0006, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0008, 0x04B0, 0x5E2D,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x0006, 0x0008, 0x0006, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x000D, 0x002E, 0x0062, 0x007E, 0x0062, 0x002E, 0x000D, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0002, 0x0015, 0x007E, 0x01B9, 0x03A6, 0x04B0, 0x03A6, 0x01B9, 0x007E, 0x0015, 0x0002, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x000D, 0x007E, 0x02D8, 0x09ED, 0x1503, 0x1AFB, 0x1503, 0x09ED, 0x02D8, 0x007E, 0x000D, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0002, 0x002E, 0x01B9, 0x09ED, 0x22A5, 0x4958, 0x5E2D, 0x4958, 0x22A5, 0x09ED, 0x01B9, 0x002E, 0x0002, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0006, 0x0062, 0x03A6, 0x1503, 0x4958, 0x9B45, 0xC75F, 0x9B45, 0x4958, 0x1503, 0x03A6, 0x0062, 0x0006, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0008, 0x007E, 0x04B0, 0x1AFB, 0x5E2D, 0xC75F,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0001, 0x0002, 0x0005, 0x0007, 0x0008, 0x0007, 0x0005, 0x0002, 0x0001, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0002, 0x0007, 0x0015, 0x002E, 0x0051, 0x0071, 0x007E, 0x0071, 0x0051, 0x002E, 0x0015, 0x0007, 0x0002, 0x0000, 0x0000,
	 0x0000, 0x0005, 0x0015, 0x004A, 0x00CA, 0x01B9, 0x0301, 0x0432, 0x04B0, 0x0432, 0x0301, 0x01B9, 0x00CA, 0x004A, 0x0015, 0x0005, 0x0000,
	 0x0005, 0x001D, 0x007E, 0x01AD, 0x048F, 0x09ED, 0x114C, 0x1825, 0x1AFB, 0x1825, 0x114C, 0x09ED, 0x048F, 0x01AD, 0x007E, 0x001D, 0x0005,
	 0x0013, 0x0068, 0x01B9, 0x05DB, 0x0FEA, 0x22A5, 0x3C62, 0x5445, 0x5E2D, 0x5445, 0x3C62, 0x22A5, 0x0FEA, 0x05DB, 0x01B9, 0x0068, 0x0013,
	 0x0029, 0x00DC, 0x03A6, 0x0C65, 0x21B2, 0x4958, 0x7FD5, 0xB268, 0xC75F, 0xB268, 0x7FD5, 0x4958, 0x21B2, 0x0C65, 0x03A6, 0x00DC, 0x0029,
	 0x0035, 0x011B, 0x04B0, 0x0FEA, 0x2B44, 0x5E2D, 0xA424, 0xE514,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0001, 0x0002, 0x0004, 0x0006, 0x0007, 0x0008, 0x0007, 0x0006, 0x0004, 0x0002, 0x0001, 0x0000, 0x0000, 0x0000,
	 0x0002, 0x0005, 0x000D, 0x001A, 0x002E, 0x0048, 0x0062, 0x0076, 0x007E, 0x0076, 0x0062, 0x0048, 0x002E, 0x001A, 0x000D, 0x0005, 0x0002,
	 0x0015, 0x0038, 0x007E, 0x00FB, 0x01B9, 0x02AB, 0x03A6, 0x0467, 0x04B0, 0x0467, 0x03A6, 0x02AB, 0x01B9, 0x00FB, 0x007E, 0x0038, 0x0015,
	 0x007E, 0x0143, 0x02D8, 0x05A7, 0x09ED, 0x0F5F, 0x1503, 0x1958, 0x1AFB, 0x1958, 0x1503, 0x0F5F, 0x09ED, 0x05A7, 0x02D8, 0x0143, 0x007E,
	 0x01B9, 0x0467, 0x09ED, 0x13BD, 0x22A5, 0x35A9, 0x4958, 0x5878, 0x5E2D, 0x5878, 0x4958, 0x35A9, 0x22A5, 0x13BD, 0x09ED, 0x0467, 0x01B9,
	 0x03A6, 0x0953, 0x1503, 0x29CA, 0x4958, 0x7199, 0x9B45, 0xBB4B, 0xC75F, 0xBB4B, 0x9B45, 0x7199, 0x4958, 0x29CA, 0x1503, 0x0953, 0x03A6,
	 0x04B0, 0x0BF9, 0x1AFB, 0x35A9, 0x5E2D, 0x91DD, 0xC75F, 0xF07D,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0013, 0x0035, 0x0013, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0005, 0x0068, 0x011B, 0x0068, 0x0005, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0015, 0x01B9, 0x04B0, 0x01B9, 0x0015, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x004A, 0x05DB, 0x0FEA, 0x05DB, 0x004A, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001, 0x00CA, 0x0FEA, 0x2B44, 0x0FEA, 0x00CA, 0x0001, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x01B9, 0x22A5, 0x5E2D, 0x22A5, 0x01B9, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0005, 0x0301, 0x3C62, 0xA424, 0x3C62, 0x0301, 0x0005, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0007, 0x0432, 0x5445, 0xE514, 0x5445, 0x0432, 0x0007, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0008, 0x04B0, 0x5E2D,
	 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0005, 0x0013, 0x0029, 0x0035, 0x0029, 0x0013, 0x0005, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0005, 0x001D, 0x0068, 0x00DC, 0x011B, 0x00DC, 0x0068, 0x001D, 0x0005, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0002, 0x0015, 0x007E, 0x01B9, 0x03A6, 0x04B0, 0x03A6, 0x01B9, 0x007E, 0x0015, 0x0002, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0007, 0x004A, 0x01AD, 0x05DB, 0x0C65, 0x0FEA, 0x0C65, 0x05DB, 0x01AD, 0x004A, 0x0007, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0001, 0x0015, 0x00CA, 0x048F, 0x0FEA, 0x21B2, 0x2B44, 0x21B2, 0x0FEA, 0x048F, 0x00CA, 0x0015, 0x0001, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0002, 0x002E, 0x01B9, 0x09ED, 0x22A5, 0x4958, 0x5E2D, 0x4958, 0x22A5, 0x09ED, 0x01B9, 0x002E, 0x0002, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0005, 0x0051, 0x0301, 0x114C, 0x3C62, 0x7FD5, 0xA424, 0x7FD5, 0x3C62, 0x114C, 0x0301, 0x0051, 0x0005, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0007, 0x0071, 0x0432, 0x1825, 0x5445, 0xB268, 0xE514, 0xB268, 0x5445, 0x1825, 0x0432, 0x0071, 0x0007, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0008, 0x007E, 0x04B0, 0x1AFB, 0x5E2D, 0xC75F,
	 
	 0x0000, 0x0000, 0x0000, 0x0003, 0x0009, 0x0013, 0x0022, 0x002F, 0x0035, 0x002F, 0x0022, 0x0013, 0x0009, 0x0003, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0001, 0x0005, 0x0011, 0x002F, 0x0068, 0x00B5, 0x00FD, 0x011B, 0x00FD, 0x00B5, 0x0068, 0x002F, 0x0011, 0x0005, 0x0001, 0x0000,
	 0x0000, 0x0005, 0x0015, 0x004A, 0x00CA, 0x01B9, 0x0301, 0x0432, 0x04B0, 0x0432, 0x0301, 0x01B9, 0x00CA, 0x004A, 0x0015, 0x0005, 0x0000,
	 0x0003, 0x0011, 0x004A, 0x00FD, 0x02B0, 0x05DB, 0x0A34, 0x0E3E, 0x0FEA, 0x0E3E, 0x0A34, 0x05DB, 0x02B0, 0x00FD, 0x004A, 0x0011, 0x0003,
	 0x0009, 0x002F, 0x00CA, 0x02B0, 0x0750, 0x0FEA, 0x1BBE, 0x26B7, 0x2B44, 0x26B7, 0x1BBE, 0x0FEA, 0x0750, 0x02B0, 0x00CA, 0x002F, 0x0009,
	 0x0013, 0x0068, 0x01B9, 0x05DB, 0x0FEA, 0x22A5, 0x3C62, 0x5445, 0x5E2D, 0x5445, 0x3C62, 0x22A5, 0x0FEA, 0x05DB, 0x01B9, 0x0068, 0x0013,
	 0x0022, 0x00B5, 0x0301, 0x0A34, 0x1BBE, 0x3C62, 0x693E, 0x92E1, 0xA424, 0x92E1, 0x693E, 0x3C62, 0x1BBE, 0x0A34, 0x0301, 0x00B5, 0x0022,
	 0x002F, 0x00FD, 0x0432, 0x0E3E, 0x26B7, 0x5445, 0x92E1, 0xCCFD, 0xE514, 0xCCFD, 0x92E1, 0x5445, 0x26B7, 0x0E3E, 0x0432, 0x00FD, 0x002F,
	 0x0035, 0x011B, 0x04B0, 0x0FEA, 0x2B44, 0x5E2D, 0xA424, 0xE514,
		 
	 0x0000, 0x0002, 0x0005, 0x000B, 0x0013, 0x001E, 0x0029, 0x0032, 0x0035, 0x0032, 0x0029, 0x001E, 0x0013, 0x000B, 0x0005, 0x0002, 0x0000,
	 0x0005, 0x000D, 0x001D, 0x003B, 0x0068, 0x00A1, 0x00DC, 0x0109, 0x011B, 0x0109, 0x00DC, 0x00A1, 0x0068, 0x003B, 0x001D, 0x000D, 0x0005,
	 0x0015, 0x0038, 0x007E, 0x00FB, 0x01B9, 0x02AB, 0x03A6, 0x0467, 0x04B0, 0x0467, 0x03A6, 0x02AB, 0x01B9, 0x00FB, 0x007E, 0x0038, 0x0015,
	 0x004A, 0x00BE, 0x01AD, 0x0356, 0x05DB, 0x0911, 0x0C65, 0x0EF3, 0x0FEA, 0x0EF3, 0x0C65, 0x0911, 0x05DB, 0x0356, 0x01AD, 0x00BE, 0x004A,
	 0x00CA, 0x0206, 0x048F, 0x0911, 0x0FEA, 0x18A7, 0x21B2, 0x28A5, 0x2B44, 0x28A5, 0x21B2, 0x18A7, 0x0FEA, 0x0911, 0x048F, 0x0206, 0x00CA,
	 0x01B9, 0x0467, 0x09ED, 0x13BD, 0x22A5, 0x35A9, 0x4958, 0x5878, 0x5E2D, 0x5878, 0x4958, 0x35A9, 0x22A5, 0x13BD, 0x09ED, 0x0467, 0x01B9,
	 0x0301, 0x07AD, 0x114C, 0x2267, 0x3C62, 0x5D86, 0x7FD5, 0x9A32, 0xA424, 0x9A32, 0x7FD5, 0x5D86, 0x3C62, 0x2267, 0x114C, 0x07AD, 0x0301,
	 0x0432, 0x0AB6, 0x1825, 0x3004, 0x5445, 0x8286, 0xB268, 0xD733, 0xE514, 0xD733, 0xB268, 0x8286, 0x5445, 0x3004, 0x1825, 0x0AB6, 0x0432,
	 0x04B0, 0x0BF9, 0x1AFB, 0x35A9, 0x5E2D, 0x91DD, 0xC75F, 0xF07D,
		 
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0015, 0x01B9, 0x04B0, 0x01B9, 0x0015, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0038, 0x0467, 0x0BF9, 0x0467, 0x0038, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x007E, 0x09ED, 0x1AFB, 0x09ED, 0x007E, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001, 0x00FB, 0x13BD, 0x35A9, 0x13BD, 0x00FB, 0x0001, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002, 0x01B9, 0x22A5, 0x5E2D, 0x22A5, 0x01B9, 0x0002, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0004, 0x02AB, 0x35A9, 0x91DD, 0x35A9, 0x02AB, 0x0004, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0006, 0x03A6, 0x4958, 0xC75F, 0x4958, 0x03A6, 0x0006, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0007, 0x0467, 0x5878, 0xF07D, 0x5878, 0x0467, 0x0007, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0008, 0x04B0, 0x5E2D,
		 
	 0x0000, 0x0000, 0x0000, 0x0002, 0x0015, 0x007E, 0x01B9, 0x03A6, 0x04B0, 0x03A6, 0x01B9, 0x007E, 0x0015, 0x0002, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x0005, 0x0038, 0x0143, 0x0467, 0x0953, 0x0BF9, 0x0953, 0x0467, 0x0143, 0x0038, 0x0005, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0000, 0x000D, 0x007E, 0x02D8, 0x09ED, 0x1503, 0x1AFB, 0x1503, 0x09ED, 0x02D8, 0x007E, 0x000D, 0x0000, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0001, 0x001A, 0x00FB, 0x05A7, 0x13BD, 0x29CA, 0x35A9, 0x29CA, 0x13BD, 0x05A7, 0x00FB, 0x001A, 0x0001, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0002, 0x002E, 0x01B9, 0x09ED, 0x22A5, 0x4958, 0x5E2D, 0x4958, 0x22A5, 0x09ED, 0x01B9, 0x002E, 0x0002, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0004, 0x0048, 0x02AB, 0x0F5F, 0x35A9, 0x7199, 0x91DD, 0x7199, 0x35A9, 0x0F5F, 0x02AB, 0x0048, 0x0004, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0006, 0x0062, 0x03A6, 0x1503, 0x4958, 0x9B45, 0xC75F, 0x9B45, 0x4958, 0x1503, 0x03A6, 0x0062, 0x0006, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0007, 0x0076, 0x0467, 0x1958, 0x5878, 0xBB4B, 0xF07D, 0xBB4B, 0x5878, 0x1958, 0x0467, 0x0076, 0x0007, 0x0000, 0x0000,
	 0x0000, 0x0000, 0x0008, 0x007E, 0x04B0, 0x1AFB, 0x5E2D, 0xC75F,
		 
	 0x0000, 0x0005, 0x0015, 0x004A, 0x00CA, 0x01B9, 0x0301, 0x0432, 0x04B0, 0x0432, 0x0301, 0x01B9, 0x00CA, 0x004A, 0x0015, 0x0005, 0x0000,
	 0x0002, 0x000D, 0x0038, 0x00BE, 0x0206, 0x0467, 0x07AD, 0x0AB6, 0x0BF9, 0x0AB6, 0x07AD, 0x0467, 0x0206, 0x00BE, 0x0038, 0x000D, 0x0002,
	 0x0005, 0x001D, 0x007E, 0x01AD, 0x048F, 0x09ED, 0x114C, 0x1825, 0x1AFB, 0x1825, 0x114C, 0x09ED, 0x048F, 0x01AD, 0x007E, 0x001D, 0x0005,
	 0x000B, 0x003B, 0x00FB, 0x0356, 0x0911, 0x13BD, 0x2267, 0x3004, 0x35A9, 0x3004, 0x2267, 0x13BD, 0x0911, 0x0356, 0x00FB, 0x003B, 0x000B,
	 0x0013, 0x0068, 0x01B9, 0x05DB, 0x0FEA, 0x22A5, 0x3C62, 0x5445, 0x5E2D, 0x5445, 0x3C62, 0x22A5, 0x0FEA, 0x05DB, 0x01B9, 0x0068, 0x0013,
	 0x001E, 0x00A1, 0x02AB, 0x0911, 0x18A7, 0x35A9, 0x5D86, 0x8286, 0x91DD, 0x8286, 0x5D86, 0x35A9, 0x18A7, 0x0911, 0x02AB, 0x00A1, 0x001E,
	 0x0029, 0x00DC, 0x03A6, 0x0C65, 0x21B2, 0x4958, 0x7FD5, 0xB268, 0xC75F, 0xB268, 0x7FD5, 0x4958, 0x21B2, 0x0C65, 0x03A6, 0x00DC, 0x0029,
	 0x0032, 0x0109, 0x0467, 0x0EF3, 0x28A5, 0x5878, 0x9A32, 0xD733, 0xF07D, 0xD733, 0x9A32, 0x5878, 0x28A5, 0x0EF3, 0x0467, 0x0109, 0x0032,
	 0x0035, 0x011B, 0x04B0, 0x0FEA, 0x2B44, 0x5E2D, 0xA424, 0xE514,
		 
	 0x0015, 0x0038, 0x007E, 0x00FB, 0x01B9, 0x02AB, 0x03A6, 0x0467, 0x04B0, 0x0467, 0x03A6, 0x02AB, 0x01B9, 0x00FB, 0x007E, 0x0038, 0x0015,
	 0x0038, 0x008F, 0x0143, 0x0282, 0x0467, 0x06D2, 0x0953, 0x0B3F, 0x0BF9, 0x0B3F, 0x0953, 0x06D2, 0x0467, 0x0282, 0x0143, 0x008F, 0x0038,
	 0x007E, 0x0143, 0x02D8, 0x05A7, 0x09ED, 0x0F5F, 0x1503, 0x1958, 0x1AFB, 0x1958, 0x1503, 0x0F5F, 0x09ED, 0x05A7, 0x02D8, 0x0143, 0x007E,
	 0x00FB, 0x0282, 0x05A7, 0x0B3F, 0x13BD, 0x1E93, 0x29CA, 0x3268, 0x35A9, 0x3268, 0x29CA, 0x1E93, 0x13BD, 0x0B3F, 0x05A7, 0x0282, 0x00FB,
	 0x01B9, 0x0467, 0x09ED, 0x13BD, 0x22A5, 0x35A9, 0x4958, 0x5878, 0x5E2D, 0x5878, 0x4958, 0x35A9, 0x22A5, 0x13BD, 0x09ED, 0x0467, 0x01B9,
	 0x02AB, 0x06D2, 0x0F5F, 0x1E93, 0x35A9, 0x531C, 0x7199, 0x8906, 0x91DD, 0x8906, 0x7199, 0x531C, 0x35A9, 0x1E93, 0x0F5F, 0x06D2, 0x02AB,
	 0x03A6, 0x0953, 0x1503, 0x29CA, 0x4958, 0x7199, 0x9B45, 0xBB4B, 0xC75F, 0xBB4B, 0x9B45, 0x7199, 0x4958, 0x29CA, 0x1503, 0x0953, 0x03A6,
	 0x0467, 0x0B3F, 0x1958, 0x3268, 0x5878, 0x8906, 0xBB4B, 0xE1EB, 0xF07D, 0xE1EB, 0xBB4B, 0x8906, 0x5878, 0x3268, 0x1958, 0x0B3F, 0x0467,
	 0x04B0, 0x0BF9, 0x1AFB, 0x35A9, 0x5E2D, 0x91DD, 0xC75F, 0xF07D,
};
void pred_grad2(char *buf, int iw, int ih, int fwd)
{
#if 1
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	int *prederrors=(int*)malloc(iw*_countof(g2_sums)*sizeof(int));
	if(!b2||!prederrors)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, buf, (size_t)res<<2);//copy alpha
	const char *pixels=fwd?buf:b2, *errors=fwd?b2:buf;

	//int vmin2=0, vmax2=0;
	for(int kc=0;kc<3;++kc)
	{
		memset(prederrors, 0, iw*_countof(g2_sums)*sizeof(int));
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				short nb[G2_NNB];
				int idx=0;
				for(int ky2=-G2_REACH;ky2<0;++ky2)
				{
					for(int kx2=-G2_REACH;kx2<=G2_REACH;++kx2)
					{
						nb[idx]=(unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih?pixels[(iw*(ky+ky2)+kx+kx2)<<2|kc]:0;
						++idx;
					}
				}
				for(int kx2=-G2_REACH;kx2<0;++kx2)
				{
					nb[idx]=(unsigned)(kx+kx2)<(unsigned)iw?pixels[(iw*ky+kx+kx2)<<2|kc]:0;
					++idx;
				}

				int pred;
				long long num;
				int preds[_countof(g2_sums)], den;//floor_log2(G2_REACH)^2

				num=0, den=0;
				for(int k=0;k<_countof(g2_sums);++k)
				{
					pred=fast_dot(nb, g2_kernels+G2_NNB*k, G2_NNB);
					preds[k]=(int)(((long long)pred<<8)/g2_sums[k]);
					int e=(unsigned)(kx-1)<(unsigned)iw?prederrors[iw*k+kx-1]:0;
					e+=(unsigned)(ky-1)<(unsigned)ih?prederrors[iw*k+kx]:0;
					int w=0x1000000/(e+1);
					num+=(long long)preds[k]*w;
					den+=w;
				}
				pred=den?(int)(num/den):preds[0];

				//if(vmin2>pred)vmin2=pred;//
				//if(vmax2<pred)vmax2=pred;
				//pred=CLAMP(-(128<<8), pred, 127<<8);
				idx=(iw*ky+kx)<<2|kc;


				//if(kx==(iw>>1)&&ky==(ih>>1))//
				//	printf("");
				
				if(fwd)
					b2[idx]=buf[idx]-((pred+128)>>8);
				else
					b2[idx]=buf[idx]+((pred+128)>>8);

				int curr=pixels[idx]<<8;
				for(int k=0;k<_countof(g2_sums);++k)
				{
					int e=abs(curr-preds[k]);
					prederrors[iw*k+kx]=e;
					if(kx<iw&&ky>0)
						prederrors[iw*k+kx+1]+=e;
				}
			}
		}
	}
	//console_start();
	//console_log("[%d %d]\n", vmin2, vmax2);
	memcpy(buf, b2, (size_t)res<<2);
	free(b2);
#endif
#if 0
#define KERNEL(X, FX, Y, FY) exp(-((X)*(X)/(double)((FX)*(FX))+(Y)*(Y)/(double)((FY)*(FY))))
	console_start();
	int nf=floor_log2(G2_REACH)+1;
	for(int fy=1;fy<=nf;++fy)
	{
		for(int fx=1;fx<=nf;++fx)
		{
			console_log("XY %d %d:\n", fx, fy);
			int sum=0;
			for(int ky=-G2_REACH;ky<0;++ky)
			{
				for(int kx=-G2_REACH;kx<=G2_REACH;++kx)
				{
					double coeff=KERNEL(kx, fx, ky, fy);
					if(!isfinite(coeff))
						LOG_ERROR("");
					int val=(int)(0x10000*coeff);
					console_log("%c0x%04X,", val<0?'-':' ', abs(val));
					sum+=val;
				}
				console_log("\n");
			}
			for(int kx=-G2_REACH;kx<0;++kx)
			{
				double coeff=KERNEL(kx, fx, 0, 1);
				if(!isfinite(coeff))
					LOG_ERROR("");
				int val=(int)(0x10000*coeff);
				console_log("%c0x%04X,", val<0?'-':' ', abs(val));
				sum+=val;
			}
			console_log("\n");
			console_log("Sum 0x%08X\n", sum);
			console_log("\n");
		}
	}
	console_log("Done.\n");
	console_pause();
#endif
}
#elif defined G2_MM
#define G2_NPRED 21
//#define G2_NPRED (21+12)
short g2_weights[]=
{
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 -0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,

	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,

	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,
};
void pred_grad2(char *buf, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	int *perrors=(int*)malloc(iw*(G2_NPRED+1)*2*sizeof(int));
	if(!b2||!perrors)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memcpy(b2, buf, (size_t)res<<2);//copy alpha
	const char *pixels=fwd?buf:b2, *errors=fwd?b2:buf;
	for(int kc=0;kc<3;++kc)
	{
		int maxerror=0;
		short *params=g2_weights+(_countof(g2_weights)/3)*kc;
		int *hireserror=perrors+iw*2*G2_NPRED;
		memset(perrors, 0, 2*iw*(G2_NPRED+1)*sizeof(int));
		//int weights[3]=
		//{
		//	0x8000,
		//	0x8000,
		//	0x8000,
		//};
		//int certainty=512;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?pixels[(iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0
				int
					NNNNNN  =LOAD( 0, -6),
					NNNNWWWW=LOAD(-4, -4),
					NNNN    =LOAD( 0, -4),
					NNNNEEEE=LOAD( 4, -4),
					NNNWWW  =LOAD(-3, -3),
					NNN     =LOAD( 0, -3),
					NNNEEE  =LOAD( 3, -3),
					NNWW    =LOAD(-2, -2),
					NNW     =LOAD(-1, -2),
					NN      =LOAD( 0, -2),
					NNE     =LOAD( 1, -2),
					NNEE    =LOAD( 2, -2),
					NW      =LOAD(-1, -1),
					N       =LOAD( 0, -1),
					NE      =LOAD( 1, -1),
					NEEEE   =LOAD( 4, -1),
					NEEEEE  =LOAD( 5, -1),
					NEEEEEE =LOAD( 6, -1),
					NEEEEEEE=LOAD( 7, -1),
					NEE     =LOAD( 2, -1),
					WWWWWW  =LOAD(-6,  0),
					WWWW    =LOAD(-4,  0),
					WWW     =LOAD(-3,  0),
					WW      =LOAD(-2,  0),
					W       =LOAD(-1,  0);
#undef  LOAD
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?hireserror[iw*((ky+(Y))&1)+kx+(X)]:0
				int
					eNW=LOAD(-1, -1),
					eN =LOAD( 0, -1),
					eNE=LOAD( 1, -1),
					eW =LOAD(-1,  0);
#undef  LOAD
				int preds[G2_NPRED], pred;
				int prederrs[G2_NPRED];

				for(int k=0;k<G2_NPRED;++k)
				{
					//eNW + eN + eNE
					prederrs[k]=
						((unsigned)(kx-1)<(unsigned)iw?perrors[iw*2*k+iw*((ky-1)&1)+kx-1]:0)+
						perrors[iw*2*k+iw*((ky-1)&1)+kx]+
						((unsigned)(kx+1)<(unsigned)iw?perrors[iw*2*k+iw*((ky-1)&1)+kx+1]:0);
				}
				
				//const int den2=2;
				const int correction=0;
				//int correction=((eN+eW)*5+(eNW+eNE)*3)/(16*4);
				int j=-1;
				int vmin, vmax;
#if 1
#define GRAD(pred, N, W, NW, vmin, vmax)\
				do\
				{\
					if(N<W)\
						vmin=N, vmax=W;\
					else\
						vmin=W, vmax=N;\
					pred=N+W-NW;\
					pred=CLAMP(vmin, pred, vmax);\
				}while(0)
#endif
				vmin=N, vmax=N;
				if(vmin>W)vmin=W;
				if(vmax<W)vmax=W;
#if 1
				//the 4 predictors from JPEG XL:
				++j, preds[j]=N-((eNW*params[G2_NPRED]+eN*params[G2_NPRED+1]+eNE*params[G2_NPRED+2]+(NN-N)*params[G2_NPRED+3]+(NW-W)*params[G2_NPRED+4])>>8);
				++j, preds[j]=W-((eN+eW+eNW)*params[G2_NPRED+5]>>8);
				++j, preds[j]=N-((eN+eW+eNE)*params[G2_NPRED+6]>>8);
				++j, preds[j]=W+NE-N;
				//++j, preds[j]=N+((eN+eW+eNE)*10>>5);
				//++j, preds[j]=W+((eN+eW+eNW)*10>>5);
				//++j, preds[j]=N+((eNW*5+eN*5+eNE*5+(NN-N)*12+(NW-W)*4)>>5);
				
				++j, preds[j]=W -(eW *params[G2_NPRED+ 7]>>8);
				++j, preds[j]=N -(eN *params[G2_NPRED+ 8]>>8);
				++j, preds[j]=NW-(eNW*params[G2_NPRED+ 9]>>8);
				++j, preds[j]=NE-(eNE*params[G2_NPRED+10]>>8);
				//++j, preds[j]=W+(eW*10>>5);
				//++j, preds[j]=N+(eN*10>>5);
				//++j, preds[j]=NW+(eNW*8>>5);
				//++j, preds[j]=NE+(eNE*8>>5);
				//++j, preds[j]=W+(prederrs[j]*10>>5);
				//++j, preds[j]=N+(prederrs[j]*10>>5);
				//++j, preds[j]=NW+(prederrs[j]*8>>5);
				//++j, preds[j]=NE+(prederrs[j]*8>>5);
				//++j, preds[j]=W+correction;
				//++j, preds[j]=N+correction;
				//++j, preds[j]=NW+correction;
				//++j, preds[j]=NE+correction;
				
				++j, preds[j]=clamp4(N+W -NW +correction, N, W, NW, NE);
				++j, preds[j]=clamp4(W+NE-N  +correction, N, W, NW, NE);
				++j, preds[j]=clamp4(N+NW-NNW+correction, N, W, NW, NE);
				++j, preds[j]=clamp4(N+NE-NNE+correction, N, W, NE, NEE);
				
				++j, preds[j]=(W+NEE)/2+correction;
				++j, preds[j]=NNNNNN+correction;
				++j, preds[j]=(NEEEE+NEEEEEE)/2+correction;
				++j, preds[j]=(WWWW+WWWWWW)/2+correction;
				//++j, preds[j]=W*3-WW*3+WWW;//parabolic
				//++j, preds[j]=N*3-NN*3+NNN;
				//++j, preds[j]=NW*3-NNWW*3+NNNWWW;
				//++j, preds[j]=NE*3-NNEE*3+NNNEEE;
				//++j, preds[j]=W*2-WW;//linear
				//++j, preds[j]=N*2-NN;
				//++j, preds[j]=NW*2-NNWW;
				//++j, preds[j]=NE*2-NNEE;
				//++j, preds[j]=4*W-6*WW+4*WWW-WWWW, preds[j]=CLAMP(vmin, preds[j], vmax);//cubic
				//++j, preds[j]=4*N-6*NN+4*NNN-NNNN, preds[j]=CLAMP(vmin, preds[j], vmax);
				//++j, preds[j]=4*NW-6*NNWW+4*NNNWWW-NNNNWWWW, preds[j]=CLAMP(vmin, preds[j], vmax);
				//++j, preds[j]=4*NE-6*NNEE+4*NNNEEE-NNNNEEEE, preds[j]=CLAMP(vmin, preds[j], vmax);

				++j, preds[j]=(N+W+NEEEEE+NEEEEEEE)/4+correction;
				++j, preds[j]=clamp4(N*2-NN+correction, N, W, NE, NEE);
				++j, preds[j]=(N+NNN)/2+correction;
				++j, preds[j]=((N+W)*3-NW*2)/4+correction;
#endif

				++j; GRAD(preds[j], N, W, NW, vmin, vmax);
				//++j, preds[j]=W+NE-N, preds[j]=CLAMP(vmin, preds[j], vmax);

				//++j;
				//int gx1=N-NW, gx2=W-WW, gy1=W-NW, gy2=N-NN;
				//if(abs(gx1-gx2)<abs(gy1-gy2))
				//	preds[j]=W+(gx1*11+gx2*5)/16;
				//else
				//	preds[j]=N+(gy1*11+gy2*5)/16;

				//++j, preds[j]=(W+((N-NW)*11+(W-WW)*5)/16 + N+((W-NW)*11+(N-NN)*5)/16)/2;
				//++j, preds[j]=NW+(N-NW)+(W-NW);//==N+W-NW
				//++j, preds[j]=((N+W)*5+(NE+NW)*3)/16;
				//++j, preds[j]=(N+W)/2+(NE-NW)/4;
				//++j, preds[j]=N;
				//++j, preds[j]=W;
				//++j; GRAD(preds[j], NE, NW, NN, vmin, vmax);
				//++j; GRAD(preds[j], NN, WW, NNWW, vmin, vmax);

#if 0
				++j, preds[j]=2*W -WW  , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//linear
				++j, preds[j]=2*N -NN  , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=2*NW-NNWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=2*NE-NNEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);

				++j, preds[j]=3*W -3*WW  +WWW   , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//parabolic
				++j, preds[j]=3*N -3*NN  +NNN   , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=3*NW-3*NNWW+NNNWWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=3*NE-3*NNEE+NNNEEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);

				++j, preds[j]=4*W -6*WW  -4*WWW   +WWWW    , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//cubic
				++j, preds[j]=4*N -6*NN  -4*NNN   +NNNN    , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=4*NW-6*NNWW-4*NNNWWW+NNNNWWWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=4*NE-6*NNEE-4*NNNEEE+NNNNEEEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
#endif

				long long num=0;
				int weights[G2_NPRED], den=0;
				for(int k=0;k<G2_NPRED;++k)
				{
					weights[k]=(params[k]<<8)/(prederrs[k]+1);
					//if(maxerror<prederrs[k])
					//	maxerror=prederrs[k];
					//weights[k]=maxerror-prederrs[k]+1;

					num+=(long long)preds[k]*weights[k];
					den+=weights[k];
				}
				pred=den?(int)(num/den):preds[0];
				
				pred=CLAMP(vmin, pred, vmax);
				
				
				//if(kc==1&&kx==384&&ky==256)//
				//if(pred)
				//	printf("");
				
				//int pred0=pred;
				//int energy=abs(eN)+abs(eW)+abs(eNW)+abs(eNE);
				//pred-=pred*energy/certainty;//pred=pred0*(1-energy/certainty)

				//pred-=pred*(abs(eN)+abs(eW)+abs(eNW))/1108;//improves only with kodim13

				//pred*=(1108-(abs(eN)+abs(eW)+abs(eNW)));
				//pred/=1108;

				//pred/=(abs(eN)+abs(eW)+abs(eNW))/128+1;

				//if(abs(eN)+abs(eW)+abs(eNW)>512)pred/=2;

				//pred+=(eN>0&&eW>0&&eNW>0)-(eN<0&&eW<0&&eNW<0);

				//if(eN>0&&eW>0&&eNW>0)pred=MAXVAR(vmax, NW);
				//else if(eN<0&&eW<0&&eNW<0)pred=MINVAR(vmin, NW);

				//pred=CLAMP(-128, pred, 127);

				int idx=(iw*ky+kx)<<2|kc;
				if(fwd)
					b2[idx]=buf[idx]-((pred+128)>>8);
				else
					b2[idx]=buf[idx]+((pred+128)>>8);

				//update correction

				int curr=pixels[idx]<<8;
				hireserror[iw*(ky&1)+kx]=curr-pred;
				for(int k=0;k<G2_NPRED;++k)
				{
					int e=abs(curr-preds[k]);
					perrors[iw*2*k+iw*(ky&1)+kx]=e;
					if(kx<iw&&ky>0)
						perrors[iw*2*k+iw*((ky-1)&1)+kx+1]+=e;
				}
#if 0
				int curr=pixels[idx];
				//curr=pred0*(1-energy/temperature)	->	temperature = energy/(1 - abs(curr/pred0))

				if(abs(curr)<abs(pred))
				{
					certainty-=10;
					if(certainty<512)
						certainty=512;
				}
				else
					++certainty;

				//certainty+=abs(curr)<abs(pred)?-1:1;

				//int den=pred0?0x10000-abs((curr<<16)/pred0):1;
				//if(den>0)
				//{
				//	int t2=(energy<<16)/den;
				//	if(t2>512)
				//		temperature=t2;
				//}

				//int num=pred0*energy, den=pred0-curr;
				//if(num&&den&&abs(num)>abs(den))
				//	temperature=abs(num/den);
#endif
			}
		}
	}
	memcpy(buf, b2, (size_t)res<<2);
	free(perrors);
	free(b2);
}
#endif

#define LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|kc]:0)
static void pred_wu97_pass3_row(const char *src, char *dst, int iw, int ih, int fwd, int kc, int kx0, int ky)
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	for(int kx=kx0;kx<iw;kx+=2)
	{
		char
			NW=LOAD(pixels, -1, -1),
			N =LOAD(pixels,  0, -1),
			NE=LOAD(pixels,  1, -1),
			W =LOAD(pixels, -1,  0),
			E =LOAD(pixels,  1,  0),
			S =LOAD(pixels,  0,  1);
		int pred=((N+W+S+E)*3-(NW+NE)*2+4)>>3;
		pred=CLAMP(-128, pred, 127);

		int idx=(iw*ky+kx)<<2|kc;
		if(fwd)
			dst[idx]=src[idx]-pred;
		else
			dst[idx]=src[idx]+pred;
	}
}
static void pred_wu97_pass3(const char *src, char *dst, int iw, int ih, int fwd)
{
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih;ky+=2)
		{
			pred_wu97_pass3_row(src, dst, iw, ih, fwd, kc, 1, ky);
			pred_wu97_pass3_row(src, dst, iw, ih, fwd, kc, 0, ky+1);
		}
	}
}
static void pred_wu97_pass2_fwd(const char *pixels, char *errors, int iw, int ih)//diamonds to sticks
{
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass2
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NN  =LOAD(pixels,  0, -2),
					NW  =LOAD(pixels, -1, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					curr=LOAD(pixels,  0,  0),
					EE  =LOAD(pixels,  2,  0),
					SW  =LOAD(pixels, -1,  1),
					SE  =LOAD(pixels,  1,  1),
					SEEE=LOAD(pixels,  3,  1),
					SS  =LOAD(pixels,  0,  2),
					SSSE=LOAD(pixels,  1,  3);
				
				char
					dcurr=curr-SE,
					acurr=SE+(dcurr>>1),
					aE=SEEE+((char)(EE-SEEE)>>1),
					aS=SSSE+((char)(SS-SSSE)>>1);

				//if(kx==(iw>>1)&&ky==(ih>>1))//
				//if(kc==0&&kx==0&&ky==0)//
				//	printf("");

				int pred=(acurr*0xE666+(NE+NW+SW)*0x2AAB-(NN+WW)*0x0CCD-(aE+aS)*0x2666+0x8000)>>16;
				pred=CLAMP(-128, pred, 127);
				pred-=acurr;
				pred<<=1;
				//pred=CLAMP(-128, pred, 127);

				int idx=(iw*ky+kx)<<2|kc;
				errors[idx]=acurr;
				errors[idx+((iw+1)<<2)]=dcurr-pred;
				//errors[idx+((iw+1)<<2)]=dcurr;
			}
		}
#endif
#if 1
		for(int ky=(ih-2)&(-2);ky>=0;ky-=2)//pass1
		{
			for(int kx=(iw-2)&(-2);kx>=0;kx-=2)
			{
				char
					NNWW=LOAD(errors, -2, -2),
					NN  =LOAD(errors,  0, -2),
					NNEE=LOAD(errors,  2, -2),
					WW  =LOAD(errors, -2,  0);
				//int pred=(NN+WW)/2+(NNEE-NNWW)/4;

				int vmin, vmax;
				if(NN<WW)
					vmin=NN, vmax=WW;
				else
					vmin=WW, vmax=NN;
				int pred=NN+WW-NNWW;
				pred=CLAMP(vmin, pred, vmax);

				//if(kc==0&&kx==2&&ky==0)//
				//	printf("");
				//if(kc==0&&kx==0&&ky==2)//
				//	printf("");

				int idx=(iw*ky+kx)<<2|kc;
				errors[idx]-=pred;
			}
		}
#endif
	}
}
static void pred_wu97_pass2_inv(const char *errors, char *pixels, int iw, int ih)//sticks to diamonds
{
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass1
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NNWW=LOAD(pixels, -2, -2),
					NN  =LOAD(pixels,  0, -2),
					NNEE=LOAD(pixels,  2, -2),
					WW  =LOAD(pixels, -2,  0);
				//int pred=(NN+WW)/2+(NNEE-NNWW)/4;
				
				int vmin, vmax;
				if(NN<WW)
					vmin=NN, vmax=WW;
				else
					vmin=WW, vmax=NN;
				int pred=NN+WW-NNWW;
				pred=CLAMP(vmin, pred, vmax);

				//if(kc==0&&kx==2&&ky==0)//
				//	printf("");
				//if(kc==0&&kx==0&&ky==2)//
				//	printf("");

				int idx=(iw*ky+kx)<<2|kc;
				pixels[idx]=errors[idx]+pred;

				//idx+=(iw+1)<<2;
				//pixels[idx]=errors[idx];//copy diff
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass2
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NN  =LOAD(pixels,  0, -2),
					NW  =LOAD(pixels, -1, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					SW  =LOAD(pixels, -1,  1);
#if 1
				char//original
					acurr=LOAD(pixels, 0, 0),
					dcurr=LOAD(errors, 1, 1),
					aE   =LOAD(pixels, 2, 0),
					aS   =LOAD(pixels, 0, 2);
#endif
#if 0
				char//pass2 only
					acurr=LOAD(errors, 0, 0),
					dcurr=LOAD(errors, 1, 1),
					aE   =LOAD(errors, 2, 0),
					aS   =LOAD(errors, 0, 2);
#endif

				//if(kc==0&&kx==0&&ky==0)//
				//	printf("");

				int pred=(acurr*0xE666+(NE+NW+SW)*0x2AAB-(NN+WW)*0x0CCD-(aE+aS)*0x2666+0x8000)>>16;
				pred=CLAMP(-128, pred, 127);
				pred-=acurr;
				pred<<=1;
				//pred=CLAMP(-128, pred, 127);

				dcurr+=pred;

				char curr=dcurr, SE=acurr;
				SE-=curr>>1;	//SE = av - floor(diff/2)
				curr+=SE;		//curr = diff + SE

				int idx=(iw*ky+kx)<<2|kc;
				pixels[idx]=curr;
				pixels[idx+((iw+1)<<2)]=SE;
			}
		}
#endif
	}
}
#undef  LOAD
void pred_wu97(char *buf, int iw, int ih, int fwd)//'Lossless Compression of Continuous-Tone Images via Context Selection, Quantization, and Modeling'
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!b2||!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	//int black=0xFF000000;
	//memfill(b2, &black, (size_t)res<<2, sizeof(int));
	memcpy(b2, buf, (size_t)res<<2);//copy alpha
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 0, 0, 1);
	if(fwd)
	{
		pred_wu97_pass3(buf, b2, iw, ih, fwd);
		pred_wu97_pass2_fwd(buf, b2, iw, ih);
		dwt2d_lazy_fwd(b2  , (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_fwd(b2+1, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_fwd(b2+2, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);

	}
	else
	{
		dwt2d_lazy_inv(buf  , (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_inv(buf+1, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_inv(buf+2, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		pred_wu97_pass2_inv(buf, b2, iw, ih);
		pred_wu97_pass3(buf, b2, iw, ih, fwd);
	}
	array_free(&sizes);
	memcpy(buf, b2, (size_t)res<<2);
	free(temp);
	free(b2);
}




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
#if 0
	if(debug_buf)
	{
		lodepng_encode_file("buf_debug.PNG", debug_buf, 768, resolution/768, LCT_RGBA, 8);
		lodepng_encode_file("buf.PNG", buf, 768, resolution/768, LCT_RGBA, 8);
		for(int k=0;k<resolution;++k)
		{
			debug_buf[k<<2  ]-=buf[k<<2  ]-128;
			debug_buf[k<<2|1]-=buf[k<<2|1]-128;
			debug_buf[k<<2|2]-=buf[k<<2|2]-128;
		}
		lodepng_encode_file("buf_diff.PNG", debug_buf, 768, resolution/768, LCT_RGBA, 8);
		for(int k=0;k<resolution;++k)
		{
			debug_buf[k<<2  ]+=buf[k<<2  ]-128;
			debug_buf[k<<2|1]+=buf[k<<2|1]-128;
			debug_buf[k<<2|2]+=buf[k<<2|2]-128;
		}
		//for(int k=0;k<resolution;++k)
		//{
		//	if(buf[k]!=(unsigned char)debug_buf[k])
		//	{
		//		LOG_ERROR("Error");
		//	}
		//}
	}
#endif

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
				//p*=0xFF00;
				//++p;
				//p/=0x10000;
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
		//p*=0xFFFFFFFF;
		//++p;
		//p/=0x100000000;
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
	unsigned *htemp=(unsigned*)malloc(hsize*sizeof(unsigned));
	if(!htemp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(htemp, 0, hsize*sizeof(unsigned));

	int res=iw*ih;
	switch(space_not_color)
	{
	case 0://show correlation in color
		for(int k=0;k<res;++k)
		{
			unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
			int color=b<<(nbits<<1)|g<<nbits|r;

			++htemp[color];
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

				++htemp[color];
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

				++htemp[color];
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

				++htemp[color];
			}
		}
		break;
	}

	//don't calculate csize from downsampled histogram

	unsigned histmin=0;
	unsigned histmax=0;
	for(int k=0;k<hsize;++k)//get min & max
	{
		if(histmin>htemp[k])
			histmin=htemp[k];
		if(histmax<htemp[k])
			histmax=htemp[k];
	}

	ARRAY_ALLOC(char, *hist, 0, hsize, 0, 0);
	unsigned char *h2=hist[0]->data;
	for(int k=0;k<hsize;++k)//normalize
		h2[k]=htemp[k]*255/histmax;
	free(htemp);
}


#if 0
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
#endif

ArrayHandle bayes_mem[8]={0};
void bayes_estimate(unsigned char *src, int iw, int ih, int x1, int x2, int y1, int y2, int kc)
{
	if(!*bayes_mem)
	{
		for(int kb=7;kb>=0;--kb)
			ARRAY_ALLOC(BayesCounter, bayes_mem[kb], 0, 256LL<<(7-kb), 0, 0);
	}
	int val=1;
	for(int kb=7;kb>=0;--kb)//reset memory
		memfill(bayes_mem[kb]->data, &val, bayes_mem[kb]->count*bayes_mem[kb]->esize, sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			char nb[]=
			{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?src[(iw*(ky+(Y))+kx+(X))<<2|kc]-128:0
				LOAD(-2, -2), LOAD(-1, -2), LOAD( 0, -2), LOAD( 1, -2), LOAD( 2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD( 0, -1), LOAD( 1, -1), LOAD( 2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
#undef  LOAD
			};
			int pred=0;
			{
				double fpred=0;
				for(int k=0;k<12;++k)
					fpred+=(double)nb[k]/12;
					//fpred+=customparam_st[12*(kc<<1)+k]*nb[k];
				pred=(int)CLAMP(-128, fpred, 127)+128;
			}
			int context=pred;
			int idx=(iw*ky+kx)<<2|kc;
			unsigned char sym=src[idx];
			for(int kb=7;kb>=0;--kb)
			{
				int bit=sym>>kb&1;
				ArrayHandle mem=bayes_mem[kb];
				BayesCounter *ctr=(BayesCounter*)array_at(&mem, context);
				++ctr->n[bit];
				context|=bit<<(8+7-kb);
			}
		}
	}
}


#define C03_USE_FLOAT

#ifdef C03_USE_FLOAT
#define C03_DATATYPE double
#define C03_CVT_COEFF(X) (C03_DATATYPE)(X)
#define C03_CVT_PIXEL(X) (C03_DATATYPE)((X)/128.)
#define C03_CVT_PRED(X) (char)(CLAMP(-1, X, 1)*128.+0.5)
#else
#define C03_DATATYPE short
#define C03_CVT_COEFF(X) (C03_DATATYPE)((X)*0x4000)
#define C03_CVT_PIXEL(X) ((C03_DATATYPE)(X)<<8)
#define C03_CVT_PRED(X) (char)((X+128)>>8)
#endif

#define C03_REACH 3
#define C03_NNB (2*(C03_REACH+1)*C03_REACH)
#define C03_CIN (4*(C03_REACH+1)*C03_REACH)
static void c03_dense(const C03_DATATYPE *mat, const C03_DATATYPE *vec, const C03_DATATYPE *bias, C03_DATATYPE *res, int mw, int mh)
{
#ifdef C03_USE_FLOAT
#if 1
	for(int ky=0;ky<mh;++ky, mat+=mw)
	{
		__m256d sum=_mm256_setzero_pd();
		int kx;
		for(kx=0;kx<mw-(sizeof(__m256d)/sizeof(C03_DATATYPE)-1);kx+=sizeof(__m256)/sizeof(C03_DATATYPE))
		{
			__m256d ra=_mm256_loadu_pd(mat+kx);
			__m256d rb=_mm256_loadu_pd(vec+kx);
			ra=_mm256_mul_pd(ra, rb);
			sum=_mm256_add_pd(sum, ra);
		}
		__m128d hi=_mm256_extractf128_pd(sum, 1);
		__m128d lo=_mm256_extractf128_pd(sum, 0);
		lo=_mm_add_pd(hi, lo);
		double s2=lo.m128d_f64[0]+lo.m128d_f64[1];
		for(;kx<mw;++kx)
			s2+=mat[kx]*vec[kx];

		s2+=bias[ky];
		s2=s2<0?s2*0.01:s2;//LeakyReLU
		res[ky]=s2;
	}
#else
	//memset(res, 0, mh*sizeof(float));
	for(int ky=0;ky<mh;++ky, mat+=mw)
	{
		__m256 sum=_mm256_setzero_ps();
		int kx;
		for(kx=0;kx<mw-(sizeof(__m256)/sizeof(C03_DATATYPE)-1);kx+=sizeof(__m256)/sizeof(C03_DATATYPE))
		{
			__m256 ra=_mm256_loadu_ps(mat+kx);
			__m256 rb=_mm256_loadu_ps(vec+kx);
			ra=_mm256_mul_ps(ra, rb);
			sum=_mm256_add_ps(sum, ra);
		}
		__m128 hi=_mm256_extractf128_ps(sum, 1);
		__m128 lo=_mm256_extractf128_ps(sum, 0);
		lo=_mm_add_ps(hi, lo);
		hi=_mm_movehl_ps(hi, lo);
		lo=_mm_add_ps(hi, lo);
		hi=_mm_shuffle_ps(lo, lo, _MM_SHUFFLE(0, 0, 0, 1));
		lo=_mm_add_ss(hi, lo);
		float s2=_mm_cvtss_f32(lo);
		for(;kx<mw;++kx)
			s2+=mat[kx]*vec[kx];
		s2+=bias[ky];
		s2=s2<0?s2*0.01:s2;//LeakyReLU
		res[ky]=s2;
	}
#endif
#if 0
	int k;
	__m256 sum=_mm256_setzero_ps();
	for(k=0;k<count-7;k+=8)
	{
		__m256 ra=_mm256_loadu_ps(va+k);
		__m256 rb=_mm256_loadu_ps(vb+k);
		ra=_mm256_mul_ps(ra, rb);
		sum=_mm256_add_ps(sum, ra);
	}
	__m128 hi=_mm256_extractf128_ps(sum, 1);
	__m128 lo=_mm256_extractf128_ps(sum, 0);
	lo=_mm_add_ps(hi, lo);
	hi=_mm_movehl_ps(hi, lo);
	lo=_mm_add_ps(hi, lo);
	hi=_mm_shuffle_ps(lo, lo, _MM_SHUFFLE(0, 0, 0, 1));
	lo=_mm_add_ss(hi, lo);
	float s2=_mm_cvtss_f32(lo);
	for(;k<count;++k)
		s2+=va[k]*vb[k];
	return s2;
#endif
#else
#error TODO
#endif
}
void pred_c03(char *src, int iw, int ih, int fwd)
{
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return;
	}
#if 1
	static const double weight01[]=
	{
		  -0.085016950964927673,   0.121044263243675232,  -0.133760437369346619,  -0.149269431829452515,   0.029993344098329544,  -0.111950308084487915,  -0.114358991384506226,   0.019766116514801979,  -0.033142376691102982,   0.125196069478988647,  -0.096173897385597229,   0.054200652986764908,  -0.078628942370414734,  -0.099403962492942810,   0.088059775531291962,  -0.015082032419741154,   0.087087437510490417,  -0.033826321363449097,  -0.117229007184505463,  -0.045268345624208450,  -0.021088935434818268,   0.054183855652809143,   0.086013458669185638,   0.023304458707571030,   0.042773321270942688,  -0.055332981050014496,  -0.099078334867954254,   0.097148187458515167,  -0.069732740521430969,  -0.013463201001286507,   0.044459264725446701,   0.119321361184120178,   0.029356343671679497,  -0.122919887304306030,  -0.081503748893737793,  -0.160202667117118835,  -0.159139975905418396,  -0.120996922254562378,   0.000396791321691126,  -0.173649534583091736,  -0.195554465055465698,  -0.054109890013933182,   0.054450053721666336,  -0.130159571766853333,  -0.152187123894691467,   0.025569202378392220,   0.052201811224222183,  -0.145975619554519653,
		  -0.070544220507144928,   0.067921519279479980,  -0.058349817991256714,   0.047339584678411484,  -0.136621221899986267,  -0.055856414139270782,   0.114766813814640045,   0.089406348764896393,   0.133205562829971313,   0.041353125125169754,   0.076319023966789246,  -0.099271468818187714,  -0.071967035531997681,   0.053637418895959854,  -0.122146658599376678,   0.076023675501346588,  -0.072132311761379242,  -0.141777962446212769,  -0.060118224471807480,  -0.120739974081516266,   0.037572570145130157,   0.074924141168594360,   0.025188259780406952,  -0.026306115090847015,   0.065113835036754608,  -0.063755452632904053,   0.052151340991258621,  -0.068313486874103546,  -0.104620590806007385,   0.014844322577118874,  -0.077069848775863647,  -0.023491414263844490,   0.071562491357326508,  -0.045528728514909744,  -0.137169063091278076,  -0.135868832468986511,  -0.075038738548755646,  -0.036117315292358398,   0.098155535757541656,   0.065946735441684723,  -0.049135144799947739,  -0.234922707080841064,  -0.116365075111389160,   0.102719344198703766,  -0.080737866461277008,  -0.119058825075626373,  -0.022510958835482597,  -0.077321343123912811,
		  -0.003103345399722457,   0.087842144072055817,  -0.143391817808151245,  -0.049535807222127914,   0.100916579365730286,   0.003241522237658501,  -0.125092282891273499,   0.065612517297267914,   0.077299825847148895,  -0.075409896671772003,  -0.021854894235730171,  -0.103678204119205475,  -0.024461947381496429,   0.107226617634296417,  -0.136118486523628235,  -0.103559069335460663,  -0.093400083482265472,  -0.149753943085670471,  -0.080239929258823395,   0.101879097521305084,   0.007158371154218912,  -0.110818549990653992,   0.086287096142768860,  -0.028305737301707268,  -0.077712222933769226,  -0.030011299997568130,  -0.019944077357649803,   0.010586144402623177,  -0.028433635830879211,  -0.040263935923576355,  -0.112146206200122833,   0.039863865822553635,  -0.051370400935411453,   0.028520189225673676,  -0.151258409023284912,  -0.136650726199150085,   0.101653337478637695,   0.024024941027164459,   0.019601115956902504,  -0.023499563336372375,   0.050871554762125015,  -0.032546862959861755,  -0.108481608331203461,   0.061433676630258560,   0.070175118744373322,   0.028004638850688934,  -0.035762723535299301,  -0.132525518536567688,
		  -0.081284932792186737,   0.102296106517314911,   0.136263206601142883,   0.042855516076087952,  -0.076870538294315338,  -0.014772039838135242,  -0.048362351953983307,   0.125929355621337891,  -0.098532602190971375,  -0.070248529314994812,   0.070250503718852997,  -0.011117610149085522,   0.133501008152961731,  -0.028302805498242378,  -0.077848784625530243,   0.108451895415782928,  -0.135522156953811646,  -0.103397183120250702,  -0.171887651085853577,  -0.031735572963953018,  -0.025562506169080734,   0.082442678511142731,  -0.070272892713546753,  -0.087941601872444153,   0.027388237416744232,   0.050037272274494171,  -0.128448158502578735,   0.002477742731571198,   0.109615489840507507,  -0.005332966335117817,   0.097578696906566620,   0.043626759201288223,   0.024898041039705276,  -0.063138388097286224,   0.053812161087989807,  -0.034797847270965576,  -0.121874518692493439,   0.033258672803640366,  -0.158950537443161011,  -0.151366427540779114,  -0.183146297931671143,  -0.107115499675273895,  -0.176797613501548767,   0.114639632403850555,   0.054441403597593307,   0.126974672079086304,  -0.069740928709506989,  -0.159324824810028076,
		   0.099402405321598053,   0.146145790815353394,  -0.027423566207289696,   0.000607049267273396,  -0.017154222354292870,  -0.092115625739097595,   0.006149639841169119,  -0.081624135375022888,   0.071452766656875610,   0.004991088062524796,  -0.000677471747621894,   0.083376631140708923,   0.032803330570459366,   0.070147737860679626,   0.009470382705330849,   0.098682753741741180,   0.160697668790817261,   0.017104793339967728,  -0.076067067682743073,   0.106322690844535828,   0.059006389230489731,   0.103637523949146271,   0.146079868078231812,  -0.024362046271562576,   0.156848713755607605,   0.157300174236297607,   0.126126140356063843,   0.081740356981754303,  -0.088592447340488434,   0.001433958532288671,   0.143514141440391541,   0.042920328676700592,  -0.017379255965352058,   0.082152731716632843,   0.033449817448854446,   0.009909911081194878,   0.036733902990818024,  -0.003538192948326468,   0.155675649642944336,   0.010609120130538940,  -0.064754903316497803,   0.195229500532150269,   0.161670878529548645,   0.093289308249950409,   0.017553731799125671,  -0.081894651055335999,   0.090730331838130951,   0.164513304829597473,
		  -0.063882119953632355,  -0.091956682503223419,   0.086917363107204437,   0.082861229777336121,   0.087790712714195251,  -0.125367909669876099,  -0.108314462006092072,  -0.039305582642555237,  -0.025327887386083603,   0.146737098693847656,   0.021569145843386650,  -0.086808905005455017,  -0.077294662594795227,  -0.046061273664236069,  -0.097282566130161285,  -0.061429593712091446,  -0.152288302779197693,  -0.051012482494115829,   0.010970050469040871,  -0.076586306095123291,  -0.166325882077217102,  -0.141952410340309143,   0.129793688654899597,  -0.078501962125301361,  -0.022807762026786804,  -0.050103109329938889,  -0.039266265928745270,  -0.097879961133003235,   0.004987086635082960,  -0.079231068491935730,   0.042797803878784180,   0.108788490295410156,  -0.020038815215229988,  -0.068141654133796692,  -0.038331829011440277,  -0.110749237239360809,  -0.026170479133725166,   0.083834819495677948,  -0.018639404326677322,   0.083407416939735413,   0.046514645218849182,  -0.045570157468318939,   0.038987040519714355,  -0.038952659815549850,   0.106338083744049072,  -0.107322089374065399,  -0.064841672778129578,   0.041111662983894348,
		   0.002033447148278356,  -0.128494173288345337,  -0.172179520130157471,  -0.053132057189941406,  -0.019294237717986107,  -0.045939184725284576,   0.018042691051959991,  -0.043208256363868713,  -0.010939886793494225,  -0.027981603518128395,  -0.136838674545288086,  -0.188074409961700439,   0.053569190204143524,  -0.123980537056922913,  -0.121539510786533356,   0.004599544219672680,  -0.020019624382257462,  -0.118599392473697662,  -0.030526248738169670,  -0.180937007069587708,  -0.168239608407020569,   0.079411752521991730,   0.049109958112239838,  -0.091573059558868408,   0.091307409107685089,  -0.082027934491634369,   0.043096251785755157,  -0.087800055742263794,   0.000782740942668170,   0.003169046016409993,   0.049454480409622192,   0.080821685492992401,   0.032503962516784668,  -0.110203109681606293,  -0.057815067470073700,  -0.090289413928985596,  -0.148942664265632629,   0.093532666563987732,  -0.133212208747863770,  -0.103967472910881042,  -0.059278331696987152,   0.088755927979946136,   0.031388975679874420,   0.044864699244499207,  -0.013409201987087727,  -0.176358163356781006,  -0.182055354118347168,   0.014290302060544491,
		  -0.057777777314186096,  -0.032318294048309326,  -0.038893707096576691,   0.057151772081851959,   0.057432837784290314,  -0.140937954187393188,   0.126361578702926636,  -0.095305517315864563,  -0.138567224144935608,  -0.010915915481746197,  -0.021244525909423828,   0.098328858613967896,   0.128807440400123596,  -0.148049324750900269,   0.025565491989254951,   0.023042056709527969,  -0.113152697682380676,  -0.200296431779861450,   0.020550964400172234,   0.004968507681041956,  -0.146816566586494446,  -0.093152441084384918,  -0.096546903252601624,  -0.195468977093696594,   0.051256861537694931,  -0.106791034340858459,  -0.019607054069638252,   0.092737011611461639,   0.091646373271942139,  -0.123117089271545410,   0.034567881375551224,  -0.039830233901739120,   0.035666797310113907,  -0.110136777162551880,  -0.042352497577667236,  -0.144597247242927551,  -0.146304920315742493,   0.045528605580329895,   0.103850223124027252,   0.055614888668060303,   0.059520963579416275,   0.036169953644275665,   0.012445895932614803,   0.070902988314628601,   0.008376710116863251,   0.072714813053607941,  -0.113990463316440582,  -0.195442378520965576,
		   0.092876553535461426,  -0.156149715185165405,   0.010124005377292633,  -0.122321225702762604,  -0.065328277647495270,   0.037710908800363541,  -0.076757691800594330,  -0.011676868423819542,  -0.002939089899882674,   0.000419732386944816,  -0.152508780360221863,  -0.028132554143667221,   0.055674977600574493,   0.058497197926044464,  -0.136200591921806335,   0.015515110455453396,   0.036845419555902481,  -0.045972015708684921,   0.042176708579063416,   0.047472957521677017,  -0.046107776463031769,  -0.147095441818237305,   0.086734578013420105,  -0.111007690429687500,  -0.002082373714074492,  -0.002194649074226618,  -0.035402659326791763,   0.038862768560647964,  -0.133049324154853821,   0.074419192969799042,  -0.051712334156036377,  -0.091801568865776062,   0.035807773470878601,   0.004815363325178623,   0.076139286160469055,  -0.036499552428722382,  -0.018675101920962334,   0.040134187787771225,   0.086868181824684143,   0.016901032999157906,  -0.178645014762878418,  -0.215882360935211182,  -0.023003162816166878,  -0.087036751210689545,   0.026474572718143463,  -0.058824878185987473,  -0.034671653062105179,  -0.133515626192092896,
		   0.056611146777868271,  -0.076185539364814758,   0.124027848243713379,   0.052959226071834564,  -0.074911259114742279,   0.146648526191711426,   0.021868595853447914,   0.135149449110031128,  -0.073533140122890472,   0.039820961654186249,   0.111644379794597626,   0.077395223081111908,   0.031744964420795441,   0.089622624218463898,   0.020464735105633736,  -0.029063247144222260,   0.044699259102344513,   0.063411332666873932,  -0.149838536977767944,   0.008170351386070251,  -0.004105423577129841,  -0.143348664045333862,  -0.130609884858131409,   0.075912982225418091,   0.120366834104061127,   0.126062557101249695,  -0.121714830398559570,   0.108053825795650482,   0.073219761252403259,  -0.088121831417083740,   0.062297992408275604,   0.092266209423542023,  -0.131320536136627197,   0.102139100432395935,  -0.017563384026288986,  -0.022436898201704025,  -0.004064847249537706,  -0.041674848645925522,  -0.088223531842231750,   0.077535882592201233,  -0.124438688158988953,  -0.169159993529319763,  -0.148288086056709290,   0.121834479272365570,  -0.031951542943716049,   0.006822742987424135,   0.072856627404689789,  -0.089802160859107971,
		  -0.146903201937675476,  -0.040131244808435440,  -0.099174045026302338,   0.013182979077100754,  -0.109871134161949158,  -0.031630981713533401,  -0.057680044323205948,  -0.084400445222854614,  -0.084155112504959106,   0.145312175154685974,  -0.000717414193786681,  -0.103745661675930023,   0.054593231528997421,  -0.045979507267475128,   0.088664829730987549,  -0.059698935598134995,  -0.034144744277000427,  -0.175423398613929749,   0.010338572785258293,   0.082940705120563507,  -0.109391145408153534,   0.093482248485088348,   0.151648059487342834,  -0.102206416428089142,  -0.026267668232321739,  -0.058197036385536194,   0.008579411543905735,   0.035530939698219299,   0.018778985366225243,   0.092865221202373505,   0.007113330531865358,   0.161546543240547180,  -0.044132485985755920,   0.167394861578941345,  -0.064383335411548615,   0.114948339760303497,  -0.089157395064830780,  -0.028752701357007027,   0.025285758078098297,  -0.047698352485895157,   0.059664186090230942,  -0.170510545372962952,   0.099765628576278687,  -0.047102034091949463,  -0.013601510785520077,   0.033057238906621933,  -0.114280216395854950,  -0.197346389293670654,
		   0.122280992567539215,   0.006578224711120129,  -0.010647910647094250,  -0.013140119612216949,  -0.031685244292020798,  -0.063175074756145477,   0.024063615128397942,  -0.104206509888172150,   0.036879643797874451,  -0.059504412114620209,  -0.083848319947719574,  -0.023607965558767319,  -0.101393051445484161,   0.106740579009056091,  -0.106633402407169342,   0.016722068190574646,  -0.100571937859058380,  -0.110155314207077026,  -0.128676757216453552,  -0.002612124662846327,  -0.099474772810935974,  -0.035400588065385818,   0.107379622757434845,  -0.071668624877929688,  -0.050543438643217087,   0.086247362196445465,  -0.066026680171489716,   0.034625098109245300,   0.022724360227584839,  -0.106511831283569336,  -0.059838518500328064,  -0.085102543234825134,  -0.044030088931322098,  -0.082121811807155609,   0.021745026111602783,   0.012783295474946499,  -0.013608242385089397,  -0.081509508192539215,   0.108420133590698242,  -0.061685901135206223,  -0.160383433103561401,   0.044490527361631393,  -0.168767586350440979,   0.065318055450916290,   0.118280492722988129,  -0.113208316266536713,   0.010412113741040230,  -0.151566281914710999,
		  -0.047865863889455795,  -0.119142316281795502,  -0.143994912505149841,   0.028014445677399635,  -0.150423362851142883,  -0.083872124552726746,   0.129987463355064392,  -0.140762001276016235,  -0.069420374929904938,  -0.048370856791734695,   0.001117731095291674,  -0.047850009053945541,  -0.032351382076740265,   0.064959429204463959,  -0.043331716209650040,  -0.128461539745330811,  -0.145479321479797363,  -0.017337262630462646,  -0.094742111861705780,  -0.113058075308799744,   0.054383229464292526,   0.089307673275470734,   0.031934399157762527,  -0.195004299283027649,   0.049059059470891953,  -0.075490325689315796,   0.049327053129673004,   0.065011315047740936,   0.058220028877258301,  -0.002966593485325575,   0.096255362033843994,  -0.112147651612758636,  -0.117532953619956970,  -0.125722542405128479,  -0.051490399986505508,  -0.022029908373951912,  -0.145169407129287720,  -0.015085764229297638,   0.050206892192363739,  -0.145531386137008667,  -0.113423563539981842,  -0.229539483785629272,  -0.148413181304931641,  -0.055879376828670502,  -0.046686772257089615,  -0.084506504237651825,   0.076575644314289093,  -0.196140319108963013,
		   0.125909686088562012,   0.015048954635858536,  -0.081818327307701111,  -0.054355308413505554,  -0.044745888561010361,  -0.082716181874275208,   0.037313655018806458,  -0.036788076162338257,  -0.023752423003315926,   0.105426035821437836,  -0.143757611513137817,  -0.038675185292959213,  -0.073206961154937744,  -0.079028762876987457,   0.102559551596641541,  -0.133772894740104675,   0.107954040169715881,   0.180678382515907288,  -0.109996281564235687,  -0.018091225996613503,  -0.037321072071790695,  -0.041486896574497223,   0.019658079370856285,   0.073165051639080048,  -0.035658702254295349,  -0.123551025986671448,  -0.066068291664123535,   0.008931527845561504,   0.007751515135169029,  -0.053052425384521484,  -0.009962101466953754,   0.056861247867345810,   0.087234832346439362,  -0.156439274549484253,   0.049866262823343277,   0.119536451995372772,   0.038248874247074127,  -0.076201364398002625,  -0.097633324563503265,   0.064466074109077454,  -0.065742813050746918,   0.199949800968170166,  -0.090841770172119141,   0.038836091756820679,  -0.078279815614223480,   0.032206475734710693,   0.008943357504904270,   0.124266847968101501,
		  -0.064751572906970978,   0.106872133910655975,   0.040711306035518646,  -0.106023333966732025,   0.033653695136308670,   0.004909013863652945,   0.016303321346640587,   0.145494952797889709,   0.047181069850921631,  -0.107435978949069977,   0.135541275143623352,   0.099944896996021271,  -0.000561352760996670,  -0.091335244476795197,   0.055749565362930298,  -0.122136846184730530,   0.069689333438873291,  -0.216271206736564636,   0.061914209276437759,  -0.119145311415195465,   0.028486385941505432,   0.000633118208497763,  -0.154085010290145874,   0.047936819493770599,  -0.007048823870718479,   0.124076649546623230,   0.097882851958274841,  -0.037076979875564575,  -0.091679766774177551,   0.021232459694147110,   0.023366892710328102,  -0.038987722247838974,  -0.058122310787439346,   0.126173973083496094,   0.119354851543903351,   0.116714395582675934,  -0.113921411335468292,  -0.106340363621711731,   0.063595466315746307,  -0.134952917695045471,   0.034955024719238281,  -0.073001168668270111,   0.086819335818290710,  -0.080398581922054291,   0.112317144870758057,  -0.148474171757698059,  -0.008703635074198246,  -0.119694173336029053,
		   0.161711305379867554,  -0.036967080086469650,   0.033142086118459702,  -0.109966203570365906,   0.028041601181030273,   0.023791501298546791,  -0.064820617437362671,   0.125121638178825378,   0.161680877208709717,  -0.065716497600078583,  -0.089353591203689575,   0.083489857614040375,  -0.020305320620536804,   0.130189716815948486,  -0.097002796828746796,  -0.000110047170892358,   0.193340256810188293,   0.091567009687423706,   0.205757111310958862,   0.088066853582859039,   0.101215690374374390,   0.156758591532707214,  -0.065510019659996033,   0.157660961151123047,  -0.018743440508842468,  -0.070092044770717621,   0.164896354079246521,   0.012971272692084312,  -0.001261450466699898,  -0.004926422145217657,  -0.100203000009059906,  -0.050844904035329819,  -0.012702214531600475,   0.171908751130104065,   0.129082113504409790,   0.187741667032241821,   0.172071412205696106,   0.170131146907806396,   0.183487862348556519,  -0.044358812272548676,  -0.026981424540281296,  -0.008509210310876369,   0.029123324900865555,   0.111187390983104706,   0.045340258628129959,  -0.029337810352444649,  -0.037236548960208893,   0.206741690635681152,
		   0.024946853518486023,  -0.026946753263473511,  -0.031997438520193100,   0.040496576577425003,  -0.019008524715900421,   0.134063467383384705,   0.088938057422637939,  -0.053713727742433548,   0.153917327523231506,   0.033673085272312164,   0.034633748233318329,  -0.062873028218746185,   0.108121320605278015,   0.143676981329917908,   0.006848237942904234,  -0.066014766693115234,   0.125137388706207275,   0.013555021956562996,   0.177480459213256836,  -0.036737844347953796,  -0.103113815188407898,  -0.076397769153118134,   0.173974782228469849,   0.145414501428604126,   0.117065235972404480,   0.025638790801167488,  -0.049251124262809753,   0.096933573484420776,   0.059845492243766785,  -0.071799613535404205,  -0.001731216441839933,   0.030783379450440407,   0.090298675000667572,   0.000970111053902656,   0.164149522781372070,   0.063842371106147766,  -0.087369710206985474,   0.095455899834632874,  -0.062601044774055481,  -0.030136995017528534,   0.081745110452175140,   0.221583917737007141,   0.038597509264945984,  -0.010786258615553379,   0.072966180741786957,  -0.089446723461151123,   0.090575657784938812,   0.182809531688690186,
		  -0.013525537215173244,   0.133297488093376160,   0.026839636266231537,  -0.040435694158077240,   0.129999145865440369,  -0.126733109354972839,   0.105366922914981842,   0.008027857169508934,   0.072700060904026031,   0.032069955021142960,   0.074582599103450775,   0.020333698019385338,   0.100053727626800537,   0.069907903671264648,   0.115797378122806549,   0.062981620430946350,   0.076469197869300842,  -0.164255574345588684,  -0.052407927811145782,  -0.004377917852252722,   0.019347390159964561,   0.037077799439430237,  -0.111894667148590088,  -0.185458496212959290,   0.007193222176283598,   0.118014298379421234,  -0.110175363719463348,  -0.164532959461212158,  -0.069908872246742249,   0.065195828676223755,   0.033256616443395615,   0.120091669261455536,   0.066516295075416565,  -0.117417320609092712,   0.035993818193674088,  -0.065288208425045013,  -0.111695721745491028,   0.008658925071358681,   0.018438901752233505,  -0.031534980982542038,   0.028438314795494080,   0.045329723507165909,  -0.049899779260158539,  -0.155650123953819275,   0.085530936717987061,  -0.095037154853343964,  -0.161516010761260986,   0.016804387792944908,
		   0.000377165008103475,   0.023061921820044518,  -0.071946807205677032,  -0.095292955636978149,  -0.011255994439125061,   0.051704566925764084,   0.143522635102272034,   0.107004523277282715,   0.138775780797004700,  -0.050760798156261444,  -0.111171826720237732,   0.053630702197551727,   0.053242672234773636,   0.088784530758857727,  -0.064813174307346344,  -0.081711307168006897,   0.071107737720012665,   0.189060553908348083,   0.012266500853002071,   0.025283431634306908,   0.005327078979462385,   0.114698298275470734,   0.081821531057357788,   0.162966877222061157,   0.036295037716627121,   0.025842795148491859,   0.141605719923973083,   0.065142557024955750,   0.146158859133720398,  -0.046347983181476593,  -0.034420192241668701,   0.151525527238845825,  -0.034621126949787140,   0.127542257308959961,   0.038586825132369995,  -0.025873666629195213,  -0.002451492240652442,  -0.016076004132628441,   0.046427506953477859,  -0.066805824637413025,   0.064580410718917847,   0.211376771330833435,  -0.060594774782657623,   0.128011316061019897,  -0.102126814424991608,  -0.043940562754869461,   0.030562672764062881,   0.177347615361213684,
		  -0.035860050469636917,  -0.051614373922348022,   0.077248163521289825,  -0.102939277887344360,  -0.098211005330085754,   0.029827527701854706,  -0.031213572248816490,  -0.130246356129646301,  -0.095771476626396179,   0.022688768804073334,  -0.118954859673976898,  -0.008812169544398785,  -0.078328445553779602,   0.107189886271953583,  -0.007521524094045162,   0.123947277665138245,   0.097529090940952301,   0.070766046643257141,  -0.002885331399738789,   0.007666113320738077,  -0.161234155297279358,  -0.133382603526115417,  -0.112081207334995270,  -0.170659244060516357,  -0.019820285961031914,  -0.018541587516665459,   0.098202109336853027,  -0.129811644554138184,   0.088831894099712372,   0.105888485908508301,   0.101742595434188843,  -0.102774843573570251,   0.044804058969020844,   0.000120085613161791,  -0.057961948215961456,  -0.066869519650936127,  -0.123979508876800537,  -0.032885260879993439,   0.069961071014404297,  -0.007316689006984234,  -0.040984634310007095,  -0.043069094419479370,  -0.178636178374290466,  -0.126615732908248901,   0.030450925230979919,  -0.043647818267345428,  -0.125721618533134460,  -0.023562548682093620,
		  -0.067410044372081757,   0.074607007205486298,   0.117840692400932312,   0.134992077946662903,   0.052972033619880676,   0.073600657284259796,  -0.031376406550407410,   0.140268608927726746,   0.162841901183128357,   0.130389824509620667,   0.093137331306934357,  -0.117407783865928650,  -0.124446608126163483,   0.056633193045854568,  -0.048112206161022186,   0.007101788651198149,  -0.011679258197546005,   0.212209045886993408,  -0.071879372000694275,   0.002048897324129939,   0.093058116734027863,   0.125368356704711914,   0.024892466142773628,   0.059484388679265976,   0.040056899189949036,  -0.036337923258543015,   0.145951420068740845,   0.004772706888616085,  -0.022255167365074158,   0.034051265567541122,   0.145216315984725952,   0.166265264153480530,   0.089070416986942291,   0.090910963714122772,   0.106695048511028290,   0.061202351003885269,   0.028815902769565582,   0.084369443356990814,  -0.004303896334022284,  -0.056840661913156509,   0.013698461465537548,   0.042149461805820465,   0.033954918384552002,   0.160321041941642761,   0.127263307571411133,   0.134500026702880859,   0.053679242730140686,   0.151246324181556702,
		   0.053254216909408569,   0.052117262035608292,   0.147415608167648315,  -0.133413404226303101,   0.091530919075012207,   0.151589274406433105,  -0.016953740268945694,   0.016231948509812355,  -0.035864900797605515,  -0.074136503040790558,   0.073382630944252014,  -0.012733979150652885,   0.131069108843803406,  -0.021613290533423424,  -0.093420051038265228,   0.062464151531457901,  -0.063044764101505280,   0.083620391786098480,   0.029503025114536285,  -0.061954781413078308,   0.063800036907196045,  -0.137306973338127136,  -0.057793162763118744,  -0.203894168138504028,   0.030841108411550522,   0.077617980539798737,   0.109707832336425781,   0.010924233123660088,  -0.008582485839724541,  -0.034068029373884201,  -0.069321289658546448,   0.064039506018161774,  -0.076288230717182159,   0.133126303553581238,   0.011549264192581177,  -0.088986650109291077,   0.088439315557479858,  -0.043931309133768082,   0.076343730092048645,  -0.084136158227920532,  -0.005061340518295765,  -0.149679809808731079,   0.085352182388305664,  -0.099418722093105316,   0.012353688478469849,   0.066492393612861633,  -0.086776949465274811,  -0.100664451718330383,
		  -0.151032894849777222,  -0.011341165751218796,   0.109295859932899475,   0.103995338082313538,  -0.074525393545627594,  -0.059091504663228989,  -0.066026091575622559,  -0.060423213988542557,  -0.041975691914558411,  -0.022525168955326080,   0.135714843869209290,  -0.087333932518959045,   0.084017410874366760,   0.110917173326015472,  -0.036070190370082855,   0.046705186367034912,  -0.134539797902107239,  -0.171390503644943237,  -0.171984836459159851,   0.067617289721965790,   0.043865658342838287,  -0.067688472568988800,  -0.014605383388698101,  -0.100201845169067383,  -0.076780468225479126,  -0.111859984695911407,   0.063333638012409210,  -0.004839523695409298,   0.120617955923080444,  -0.020709704607725143,   0.070117995142936707,  -0.103284344077110291,  -0.128286212682723999,  -0.119394503533840179,  -0.125553712248802185,   0.087977446615695953,   0.098710939288139343,  -0.158015519380569458,   0.014430178329348564,  -0.042508702725172043,  -0.114591501653194427,  -0.088727317750453949,   0.054377775639295578,   0.103317290544509888,   0.049224454909563065,  -0.114509709179401398,   0.049453627318143845,  -0.192175507545471191,
		   0.146672189235687256,   0.098605789244174957,   0.128183037042617798,  -0.014116719365119934,   0.085591658949851990,  -0.073239542543888092,  -0.051354017108678818,   0.051271244883537292,   0.145766332745552063,  -0.012142680585384369,   0.161011621356010437,  -0.014377077110111713,   0.154245302081108093,   0.046543654054403305,   0.113480664789676666,  -0.039192479103803635,   0.060779638588428497,  -0.063511282205581665,   0.186497613787651062,   0.019120318815112114,   0.101494625210762024,   0.019387505948543549,   0.060840535908937454,   0.084813982248306274,   0.083772607147693634,   0.077649734914302826,  -0.040634930133819580,   0.000187015510164201,  -0.054149616509675980,   0.096416600048542023,  -0.006190130487084389,  -0.035185743123292923,  -0.075335033237934113,   0.093449361622333527,   0.046566683799028397,  -0.015555187128484249,   0.130357518792152405,   0.034689407795667648,   0.192162379622459412,   0.064713872969150543,   0.200294882059097290,   0.000616977456957102,   0.125989988446235657,   0.054411750286817551,  -0.006564210169017315,   0.083803355693817139,   0.007808590307831764,   0.058405131101608276,
		  -0.076337546110153198,   0.086204491555690765,  -0.093516632914543152,  -0.076072648167610168,   0.066723957657814026,  -0.030315360054373741,   0.020470106974244118,  -0.085515104234218597,   0.024793395772576332,  -0.038676206022500992,  -0.056648638099431992,   0.100433237850666046,   0.101447917520999908,  -0.107251770794391632,  -0.128152102231979370,  -0.028509987518191338,  -0.064596295356750488,  -0.066627748310565948,  -0.041361331939697266,  -0.144522652029991150,  -0.139153212308883667,  -0.123855829238891602,   0.054756894707679749,  -0.175798505544662476,  -0.135453373193740845,  -0.092796102166175842,  -0.012452377006411552,   0.060101009905338287,   0.063204959034919739,   0.034304879605770111,  -0.093415096402168274,   0.072877541184425354,   0.011929382570087910,   0.097005687654018402,  -0.043283343315124512,  -0.173329219222068787,   0.051216945052146912,  -0.059177208691835403,   0.015295757912099361,   0.057419214397668839,  -0.112448222935199738,  -0.187113970518112183,  -0.080608136951923370,  -0.090233229100704193,  -0.070058949291706085,   0.053560763597488403,  -0.148598894476890564,  -0.214715480804443359,
		  -0.027228219434618950,  -0.079507254064083099,   0.036686155945062637,  -0.067929312586784363,  -0.113506630063056946,   0.103822395205497742,  -0.073097899556159973,   0.100237496197223663,   0.091993056237697601,  -0.081930130720138550,  -0.085693553090095520,  -0.143626496195793152,   0.047946948558092117,  -0.076618589460849762,   0.003991567064076662,  -0.105476923286914825,   0.053823977708816528,  -0.122275881469249725,   0.021736653521656990,   0.096318095922470093,   0.066854491829872131,  -0.050357494503259659,  -0.087765775620937347,  -0.121192231774330139,   0.014284683391451836,   0.073694162070751190,   0.066725566983222961,   0.009336148388683796,   0.049534525722265244,   0.109590604901313782,   0.055702999234199524,   0.014730784110724926,   0.116543568670749664,  -0.001747315283864737,   0.033862814307212830,   0.039974756538867950,  -0.030919373035430908,   0.007665533106774092,   0.042873255908489227,  -0.068688347935676575,  -0.080269768834114075,  -0.146421283483505249,   0.004264321178197861,  -0.013927187770605087,  -0.155780673027038574,  -0.083506226539611816,  -0.135189175605773926,  -0.088855788111686707,
		   0.013134742155671120,   0.141419425606727600,  -0.075794987380504608,   0.037549685686826706,   0.004395066760480404,   0.046200588345527649,   0.028515690937638283,   0.082817308604717255,  -0.014110256917774677,   0.037832960486412048,   0.049588281661272049,   0.050380568951368332,  -0.011342986486852169,   0.038529954850673676,  -0.103250786662101746,   0.035564258694648743,  -0.053722288459539413,   0.202308669686317444,   0.026313060894608498,   0.159621715545654297,  -0.073499627411365509,   0.163663983345031738,   0.036754846572875977,   0.060011137276887894,   0.023057099431753159,   0.068417936563491821,   0.140645012259483337,   0.137983456254005432,  -0.070079848170280457,   0.149607673287391663,   0.047109887003898621,   0.121378399431705475,   0.050741240382194519,   0.057377718389034271,   0.046510446816682816,   0.072229146957397461,   0.004505455493927002,   0.005167886614799500,   0.138310253620147705,   0.108797155320644379,   0.015060511417686939,   0.014133279211819172,  -0.009903846308588982,   0.075530774891376495,   0.078502632677555084,   0.110270485281944275,  -0.089026339352130890,   0.219007179141044617,
		  -0.042828757315874100,   0.107460424304008484,   0.159615620970726013,   0.114373482763767242,   0.055125348269939423,   0.025062082335352898,   0.147788450121879578,  -0.065429598093032837,  -0.082529902458190918,   0.087297432124614716,   0.164037942886352539,  -0.085479930043220520,   0.065591871738433838,  -0.051735732704401016,  -0.051233332604169846,  -0.000805013231001794,   0.108169250190258026,   0.138789221644401550,   0.015298142097890377,   0.149545788764953613,  -0.089882083237171173,  -0.071342393755912781,  -0.045732777565717697,   0.213022515177726746,  -0.094835728406906128,   0.162220165133476257,   0.068320922553539276,   0.132266178727149963,  -0.089018933475017548,   0.112514555454254150,   0.083604656159877777,  -0.092097714543342590,   0.127695128321647644,   0.035805262625217438,   0.102607436478137970,   0.136481076478958130,   0.033865716308355331,  -0.026029422879219055,   0.155570000410079956,  -0.068804100155830383,  -0.000370521709555760,  -0.015922747552394867,   0.099891491234302521,   0.089216142892837524,  -0.015288206748664379,   0.086525961756706238,   0.164334982633590698,   0.028449844568967819,
		   0.000676319410558790,  -0.085510708391666412,  -0.116319276392459869,  -0.140364781022071838,  -0.073834717273712158,  -0.013926716521382332,  -0.134266972541809082,  -0.188165605068206787,  -0.078850060701370239,  -0.012042072601616383,  -0.097279123961925507,   0.074158720672130585,  -0.176739156246185303,  -0.183841079473495483,   0.011137131601572037,  -0.165803253650665283,  -0.159638971090316772,   0.114692948758602142,  -0.013527673669159412,   0.072382979094982147,  -0.123965859413146973,   0.040797889232635498,  -0.077779747545719147,   0.108914926648139954,  -0.006011036224663258,   0.047425117343664169,   0.091881342232227325,  -0.067378386855125427,   0.037277750670909882,  -0.107209861278533936,   0.085284471511840820,  -0.134506434202194214,  -0.060468889772891998,  -0.014356588013470173,  -0.042457133531570435,  -0.066462233662605286,   0.002272151410579681,   0.058392606675624847,  -0.046391367912292480,  -0.117862828075885773,   0.095075488090515137,   0.044567730277776718,   0.073747158050537109,  -0.146550670266151428,   0.059757381677627563,   0.092440061271190643,   0.022657288238406181,  -0.078031741082668304,
		   0.129929453134536743,   0.105967111885547638,  -0.017921417951583862,  -0.141040891408920288,   0.036280289292335510,  -0.067863516509532928,   0.095321863889694214,  -0.040093399584293365,  -0.152468875050544739,  -0.027463756501674652,  -0.154160782694816589,  -0.079746402800083160,   0.037971209734678268,  -0.147948578000068665,  -0.127055943012237549,   0.055877003818750381,  -0.172802075743675232,   0.008753903210163116,   0.010951705276966095,  -0.144756034016609192,  -0.024388432502746582,   0.016232231631875038,   0.004785770084708929,   0.006880389060825109,   0.044572077691555023,  -0.023670012131333351,   0.032971404492855072,  -0.186587274074554443,   0.030035518109798431,   0.057264201343059540,  -0.067846737802028656,   0.015015381388366222,  -0.147152975201606750,   0.077764354646205902,  -0.133039474487304688,  -0.123848237097263336,  -0.074422940611839294,  -0.037833455950021744,   0.035521194338798523,  -0.058409977704286575,  -0.184761673212051392,  -0.025431279093027115,  -0.016609728336334229,  -0.163855448365211487,  -0.107258431613445282,   0.018298471346497536,  -0.123234234750270844,  -0.204733759164810181,
		  -0.064223065972328186,   0.105817899107933044,  -0.029799528419971466,  -0.019217083230614662,  -0.079213343560695648,  -0.140152215957641602,   0.038595512509346008,   0.076672308146953583,  -0.060109235346317291,  -0.151554629206657410,  -0.040622729808092117,  -0.141295656561851501,   0.059982981532812119,   0.082665868103504181,  -0.124415136873722076,  -0.062308646738529205,   0.057360257953405380,   0.075932167470455170,   0.052466798573732376,  -0.118764415383338928,   0.024383721873164177,   0.050106335431337357,  -0.127365931868553162,   0.010414958931505680,   0.100248903036117554,  -0.013666702434420586,  -0.109410092234611511,   0.054759345948696136,   0.038199424743652344,   0.065440624952316284,  -0.093261040747165680,   0.041386045515537262,  -0.005250708200037479,  -0.030425827950239182,   0.070403546094894409,  -0.100217461585998535,   0.134753137826919556,  -0.089541807770729065,  -0.008141573518514633,  -0.028318407014012337,  -0.184256911277770996,  -0.130011469125747681,  -0.010160948149859905,  -0.030888790264725685,  -0.127529859542846680,  -0.126789256930351257,  -0.094544567167758942,   0.028028177097439766,
		   0.086569368839263916,   0.093550100922584534,  -0.110205054283142090,  -0.163228407502174377,   0.087347909808158875,  -0.132297262549400330,  -0.050841096788644791,  -0.125842049717903137,  -0.117724984884262085,  -0.000828967837151140,  -0.144128352403640747,  -0.039011694490909576,  -0.115014068782329559,  -0.020612500607967377,  -0.061965413391590118,   0.125891879200935364,  -0.153212517499923706,  -0.211043238639831543,  -0.136654943227767944,  -0.167121291160583496,   0.051532208919525146,   0.060598544776439667,   0.026112908497452736,  -0.200519263744354248,  -0.033239703625440598,  -0.143843963742256165,  -0.120756655931472778,  -0.137849509716033936,  -0.106055736541748047,  -0.100454948842525482,  -0.069201022386550903,  -0.075339235365390778,   0.004882707260549068,   0.028056455776095390,  -0.094315879046916962,  -0.110156193375587463,  -0.133089214563369751,  -0.075735829770565033,  -0.099075891077518463,   0.113511964678764343,  -0.031699020415544510,  -0.043968673795461655,  -0.034527476876974106,   0.034589163959026337,  -0.027507683262228966,  -0.008586908690631390,  -0.147589281201362610,  -0.163272276520729065,
		   0.095246382057666779,  -0.163208827376365662,   0.023180700838565826,  -0.163730859756469727,  -0.044898074120283127,  -0.140796110033988953,   0.122891373932361603,   0.024093963205814362,   0.119662724435329437,   0.148876190185546875,  -0.107373915612697601,  -0.004234636202454567,  -0.105164915323257446,   0.000256750732660294,  -0.045256603509187698,  -0.041911564767360687,  -0.019073979929089546,  -0.090327933430671692,  -0.009661230258643627,   0.095875956118106842,  -0.164914458990097046,  -0.067618496716022491,   0.044922206550836563,  -0.098950676620006561,   0.022149294614791870,   0.121425412595272064,   0.094611123204231262,  -0.030610308051109314,   0.034584976732730865,   0.126345843076705933,  -0.024690529331564903,  -0.023175364360213280,  -0.047695059329271317,   0.116971023380756378,  -0.006966816727072001,   0.084276251494884491,  -0.140891984105110168,  -0.155145928263664246,   0.034328605979681015,  -0.091374710202217102,  -0.120099470019340515,  -0.063886970281600952,   0.103388711810112000,   0.098979495465755463,  -0.078713349997997284,  -0.019404863938689232,  -0.034621406346559525,   0.011097948066890240,
		  -0.047748956829309464,  -0.086338534951210022,  -0.040574889630079269,  -0.036573536694049835,  -0.126425608992576599,  -0.076353780925273895,  -0.123390398919582367,  -0.153950586915016174,  -0.151130288839340210,  -0.055570796132087708,  -0.146131724119186401,   0.087979696691036224,  -0.155955195426940918,   0.012593761086463928,   0.029360443353652954,  -0.015007952228188515,  -0.043277256190776825,   0.046003077179193497,  -0.116323903203010559,   0.026343636214733124,   0.089651338756084442,  -0.087869659066200256,  -0.054888531565666199,   0.028460068628191948,  -0.143286421895027161,   0.122474670410156250,  -0.086688831448554993,   0.036901362240314484,  -0.148560240864753723,  -0.019464494660496712,  -0.092983350157737732,  -0.132029846310615540,   0.022068247199058533,   0.035307597368955612,   0.103334553539752960,   0.065401829779148102,   0.008213583379983902,   0.017653791233897209,  -0.090065099298954010,   0.033848468214273453,  -0.106555886566638947,  -0.025387261062860489,  -0.032640337944030762,  -0.133828505873680115,  -0.013211183249950409,  -0.077037043869495392,  -0.073522411286830902,   0.013243699446320534,
		   0.046499498188495636,   0.037742022424936295,   0.107056699693202972,   0.046696990728378296,   0.058405440300703049,   0.034709200263023376,  -0.090450130403041840,   0.031776823103427887,  -0.028545938432216644,  -0.047285057604312897,  -0.121534474194049835,   0.077813617885112762,   0.009456032887101173,  -0.102535784244537354,  -0.019670329988002777,  -0.100312083959579468,  -0.025639168918132782,   0.026362011209130287,   0.081884279847145081,   0.083399385213851929,   0.142711073160171509,   0.090553209185600281,   0.035950556397438049,  -0.039331320673227310,   0.020627370104193687,   0.110627509653568268,   0.011675180867314339,  -0.002382035832852125,   0.058626115322113037,  -0.036351218819618225,  -0.032471030950546265,  -0.040170643478631973,   0.091971814632415771,  -0.094908475875854492,  -0.075566895306110382,  -0.082335546612739563,  -0.094404883682727814,   0.099583730101585388,   0.099556736648082733,  -0.139753937721252441,  -0.183966845273971558,  -0.204351931810379028,  -0.176684215664863586,   0.129282146692276001,  -0.125737249851226807,  -0.097000472247600555,  -0.042868111282587051,  -0.215185955166816711,
		   0.037640977650880814,  -0.038436070084571838,  -0.027124460786581039,   0.051706109195947647,   0.023734875023365021,  -0.094711251556873322,  -0.013807144947350025,  -0.048360146582126617,  -0.006956972647458315,   0.118609704077243805,  -0.042834468185901642,  -0.107166267931461334,   0.106096960604190826,  -0.049541238695383072,   0.116104155778884888,   0.146810695528984070,  -0.064234636723995209,  -0.079604290425777435,   0.042712423950433731,   0.026443470269441605,   0.045397184789180756,   0.016141779720783234,   0.015308535657823086,  -0.054579790681600571,   0.052161835134029388,  -0.088792406022548676,   0.103612676262855530,   0.055265903472900391,   0.096590153872966766,   0.034880403429269791,  -0.045763362199068069,  -0.036577921360731125,   0.126188814640045166,  -0.049900140613317490,  -0.038109645247459412,  -0.026388727128505707,  -0.123240858316421509,   0.035405475646257401,  -0.117287814617156982,   0.057999711483716965,   0.049335848540067673,  -0.079310148954391479,  -0.131104886531829834,  -0.139374867081642151,  -0.021381516009569168,   0.124263718724250793,  -0.137337520718574524,  -0.101884685456752777,
		   0.027835603803396225,  -0.095283322036266327,   0.033460486680269241,  -0.024806493893265724,   0.092959240078926086,   0.008354718796908855,   0.012077149003744125,   0.050287090241909027,   0.021077664569020271,   0.092495314776897430,   0.011391565203666687,   0.010961403138935566,  -0.061599008738994598,   0.138240247964859009,   0.136318385601043701,   0.166069656610488892,   0.123488210141658783,   0.211740761995315552,   0.135724842548370361,   0.045738361775875092,   0.042304102331399918,  -0.058700047433376312,   0.023098342120647430,   0.164125561714172363,   0.080566897988319397,   0.085972867906093597,  -0.006635268684476614,  -0.033815309405326843,  -0.087492100894451141,   0.043733287602663040,   0.043564569205045700,   0.029886048287153244,   0.081233993172645569,   0.077126860618591309,   0.147640988230705261,  -0.059499450027942657,   0.029244165867567062,  -0.028227703645825386,   0.042180404067039490,   0.136041536927223206,   0.148606032133102417,   0.049101881682872772,   0.032726898789405823,   0.126801103353500366,   0.046514574438333511,   0.011220919899642467,   0.169431611895561218,   0.098879240453243256,
		   0.074039012193679810,  -0.130188822746276855,  -0.028905790299177170,  -0.079112835228443146,  -0.129015102982521057,   0.046009726822376251,   0.121058978140354156,   0.121510028839111328,  -0.110447242856025696,  -0.115923672914505005,  -0.159202337265014648,  -0.042148925364017487,   0.059359919279813766,  -0.027464788407087326,   0.038767877966165543,  -0.105279870331287384,  -0.093930505216121674,  -0.039405684918165207,   0.041804619133472443,   0.082154467701911926,  -0.004381546750664711,  -0.046390295028686523,  -0.150589331984519958,   0.079659536480903625,  -0.109938375651836395,   0.119043439626693726,   0.022585753351449966,  -0.089118190109729767,  -0.145805791020393372,   0.122780948877334595,  -0.035450775176286697,  -0.094503492116928101,  -0.151632010936737061,  -0.108316630125045776,  -0.127453520894050598,  -0.179500177502632141,   0.067372113466262817,  -0.058496385812759399,  -0.096798278391361237,   0.002442687517032027,  -0.076358444988727570,   0.029809135943651199,  -0.009540378116071224,  -0.001372943748719990,  -0.035917796194553375,  -0.012926465831696987,  -0.043593555688858032,  -0.004957829136401415,
		   0.049145936965942383,   0.157150208950042725,  -0.066094748675823212,  -0.040830660611391068,  -0.039102327078580856,  -0.004550330806523561,   0.094093516469001770,   0.032574456185102463,   0.069707229733467102,  -0.061765246093273163,   0.101591430604457855,   0.072846807539463043,   0.031205767765641212,   0.107656516134738922,  -0.057959288358688354,  -0.117480494081974030,   0.164910987019538879,   0.105376794934272766,   0.109332121908664703,  -0.096481248736381531,  -0.105790793895721436,   0.167993873357772827,  -0.028977353125810623,   0.169474571943283081,   0.029208524152636528,  -0.046261399984359741,   0.148047521710395813,  -0.077854730188846588,   0.039667416363954544,  -0.096585132181644440,   0.127580583095550537,   0.127302989363670349,   0.138883695006370544,  -0.027924399822950363,  -0.064222231507301331,   0.046775896102190018,   0.007625832222402096,  -0.062650240957736969,   0.014238627627491951,   0.054787680506706238,   0.182090520858764648,  -0.019291607663035393,   0.186371877789497375,   0.014631564728915691,   0.104749672114849091,   0.030191797763109207,   0.160923615097999573,   0.088815934956073761,
		   0.126072302460670471,   0.123127892613410950,   0.048576384782791138,  -0.032558135688304901,  -0.114627033472061157,   0.120101332664489746,  -0.107651472091674805,  -0.036826603114604950,  -0.118823520839214325,   0.000382982543669641,   0.012781832367181778,  -0.144826263189315796,  -0.021506715565919876,  -0.083542995154857635,   0.026589615270495415,  -0.120876170694828033,  -0.074524447321891785,  -0.068375974893569946,  -0.144476071000099182,  -0.087073192000389099,   0.067339383065700531,  -0.110781766474246979,  -0.147234618663787842,  -0.002313370350748301,  -0.054690010845661163,  -0.029183775186538696,   0.017072986811399460,  -0.019212510436773300,   0.099016517400741577,  -0.151473358273506165,   0.047955196350812912,  -0.105839647352695465,   0.012174546718597412,  -0.023512937128543854,  -0.121587134897708893,   0.016024140641093254,   0.077021896839141846,  -0.111451216042041779,  -0.074013791978359222,  -0.130977973341941833,   0.020076390355825424,  -0.116957142949104309,  -0.010407884605228901,   0.092777639627456665,   0.029452336952090263,  -0.060476217418909073,   0.020025651901960373,  -0.031454794108867645,
		   0.120768681168556213,   0.000718572875484824,  -0.070353232324123383,   0.050116021186113358,  -0.079647503793239594,   0.013008507899940014,  -0.038263674825429916,   0.008538600057363510,   0.116355277597904205,  -0.103902243077754974,  -0.127355784177780151,  -0.118764951825141907,  -0.051339454948902130,   0.113616250455379486,  -0.027904145419597626,   0.076250746846199036,   0.028300717473030090,  -0.177770838141441345,  -0.048729587346315384,   0.025360573083162308,   0.086697317659854889,   0.037120025604963303,   0.032039858400821686,  -0.189633756875991821,   0.058085441589355469,  -0.041077520698308945,  -0.066415809094905853,  -0.151290476322174072,   0.026299102231860161,  -0.121439330279827118,  -0.035629108548164368,   0.106405235826969147,  -0.006092450581490993,  -0.013792941346764565,   0.022374233230948448,   0.064655542373657227,  -0.016240991652011871,  -0.085499413311481476,  -0.040209222584962845,   0.067220337688922882,   0.068874388933181763,   0.004012745339423418,  -0.027702143415808678,   0.010833236388862133,  -0.132610082626342773,  -0.057887885719537735,  -0.083918623626232147,  -0.221613481640815735,
		  -0.019639575853943825,   0.097329854965209961,  -0.102737978100776672,  -0.104181133210659027,  -0.084600485861301422,   0.047790497541427612,   0.063985481858253479,   0.069516532123088837,  -0.044248931109905243,   0.043706338852643967,   0.005224118009209633,  -0.039544712752103806,   0.143241226673126221,  -0.131393015384674072,   0.078047983348369598,   0.068782247602939606,   0.099513053894042969,   0.211590260267257690,   0.029159938916563988,  -0.092473044991493225,   0.064039610326290131,   0.055964477360248566,   0.095205940306186676,   0.194217860698699951,   0.062176872044801712,   0.022191999480128288,   0.029748931527137756,   0.110934890806674957,  -0.084004767239093781,   0.056059438735246658,  -0.052039925009012222,   0.055620536208152771,  -0.091641604900360107,   0.057467177510261536,   0.060787238180637360,   0.025900511071085930,  -0.032570846378803253,   0.170282393693923950,   0.045919887721538544,  -0.073007568717002869,   0.095488451421260834,  -0.012043260037899017,   0.180170521140098572,   0.015375065617263317,   0.024109877645969391,   0.117864102125167847,  -0.033715523779392242,   0.218541100621223450,
		   0.098706953227519989,  -0.075807780027389526,  -0.039339333772659302,   0.120866619050502777,  -0.122214533388614655,  -0.046351674944162369,  -0.060904759913682938,   0.093150429427623749,   0.052712410688400269,  -0.088584683835506439,   0.133784487843513489,   0.035693954676389694,   0.132750838994979858,   0.063042283058166504,   0.074801728129386902,  -0.119003608822822571,   0.058677624911069870,  -0.069883286952972412,  -0.072932697832584381,   0.026325594633817673,  -0.022836254909634590,  -0.120166853070259094,  -0.025660518556833267,  -0.135781764984130859,   0.033783655613660812,   0.076213084161281586,   0.012972851283848286,   0.086706802248954773,  -0.146069779992103577,  -0.115743957459926605,  -0.092381246387958527,  -0.001091693178750575,  -0.070278041064739227,  -0.114386461675167084,   0.006091843359172344,  -0.012271012179553509,  -0.012298453599214554,   0.064397737383842468,  -0.024315891787409782,  -0.093049883842468262,  -0.012274019420146942,  -0.185638457536697388,  -0.171812519431114197,  -0.041898213326931000,  -0.081414036452770233,  -0.020394826307892799,  -0.144756913185119629,  -0.231727346777915955,
		  -0.068135112524032593,  -0.154385164380073547,   0.031679131090641022,  -0.058155968785285950,  -0.125703185796737671,   0.009666500613093376,  -0.052275616675615311,   0.000437130016507581,   0.005321977194398642,  -0.106495521962642670,  -0.183597683906555176,  -0.149389624595642090,  -0.049342382699251175,  -0.018207279965281487,   0.075721666216850281,   0.044385831803083420,  -0.164314836263656616,   0.020053006708621979,   0.079120844602584839,  -0.003239954821765423,   0.053367786109447479,  -0.076823279261589050,  -0.062116291373968124,   0.020445656031370163,  -0.085600137710571289,   0.097155578434467316,  -0.077037349343299866,  -0.040669035166501999,  -0.023806724697351456,   0.045737531036138535,   0.025571977719664574,   0.072298757731914520,   0.087862238287925720,  -0.164345383644104004,  -0.051130965352058411,   0.014957599341869354,   0.129311919212341309,  -0.042797062546014786,  -0.072765566408634186,  -0.122207857668399811,   0.054202746599912643,  -0.056434649974107742,  -0.091840773820877075,   0.093289285898208618,  -0.140554890036582947,  -0.069226838648319244,  -0.025148216634988785,   0.090776346623897552,
		  -0.032821830362081528,   0.103224806487560272,  -0.094587698578834534,   0.156625360250473022,   0.062292471528053284,   0.129149496555328369,  -0.071340449154376984,  -0.098212860524654388,   0.052407175302505493,   0.034968946129083633,   0.012703401036560535,   0.110496520996093750,   0.121259391307830811,   0.090809814631938934,  -0.090228945016860962,   0.119329407811164856,   0.073152944445610046,   0.106873624026775360,   0.112267129123210907,   0.087888546288013458,   0.156070217490196228,  -0.037717074155807495,   0.047231461852788925,   0.189183503389358521,  -0.028978303074836731,  -0.015194400213658810,  -0.098094977438449860,   0.158999741077423096,   0.002698512049391866,  -0.092593804001808167,   0.020604671910405159,   0.087833501398563385,  -0.045845545828342438,  -0.035146679729223251,   0.040721084922552109,   0.038799528032541275,  -0.058591030538082123,  -0.071666993200778961,   0.048715461045503616,   0.045342773199081421,   0.205716013908386230,   0.188528880476951599,  -0.072928860783576965,  -0.000307910871924832,   0.156634107232093811,   0.048884902149438858,   0.064654424786567688,   0.090113766491413116,
		   0.066499441862106323,  -0.009680038318037987,  -0.000740819261409342,   0.076120205223560333,  -0.003885419107973576,  -0.053393017500638962,  -0.017729839310050011,  -0.051081549376249313,   0.100127495825290680,  -0.101659603416919708,  -0.086794570088386536,  -0.045329831540584564,   0.044860750436782837,   0.111485481262207031,   0.112011320888996124,  -0.061313334852457047,  -0.086745753884315491,   0.099350310862064362,   0.109997160732746124,  -0.053470585495233536,  -0.091510318219661713,   0.101848661899566650,   0.139611899852752686,   0.082634516060352325,  -0.070147663354873657,  -0.017901994287967682,   0.106651201844215393,   0.125676527619361877,   0.029335439205169678,  -0.093411386013031006,  -0.043923292309045792,  -0.032095026224851608,   0.152988970279693604,  -0.099073804914951324,  -0.004373109899461269,  -0.031288515776395798,  -0.031214036047458649,   0.019520720466971397,   0.155119478702545166,  -0.081156164407730103,   0.180501058697700500,   0.129919782280921936,   0.088928483426570892,   0.076648771762847900,   0.015430541709065437,   0.134887889027595520,   0.003278828924521804,   0.114164941012859344,
		  -0.086293965578079224,   0.119624741375446320,   0.008244000375270844,   0.159088596701622009,  -0.068329446017742157,   0.106588892638683319,   0.177362620830535889,   0.079597637057304382,   0.078835204243659973,  -0.060263507068157196,   0.078545339405536652,   0.134617120027542114,  -0.065959751605987549,   0.124900117516517639,  -0.072269022464752197,  -0.078037157654762268,  -0.067647188901901245,   0.022924676537513733,  -0.018228696659207344,   0.000812610378488898,   0.040418989956378937,   0.062638238072395325,   0.098094783723354340,   0.171971052885055542,   0.057360976934432983,  -0.045646637678146362,  -0.094542980194091797,   0.009176282212138176,  -0.082924939692020416,   0.055218767374753952,  -0.106163509190082550,  -0.112972475588321686,   0.056607622653245926,   0.019058668985962868,   0.109896905720233917,   0.125512957572937012,  -0.102543972432613373,   0.042997755110263824,   0.012635045684874058,  -0.049482360482215881,   0.004376615863293409,  -0.095338195562362671,  -0.010711866430938244,   0.167276799678802490,   0.043381735682487488,  -0.060115795582532883,  -0.084237061440944672,  -0.004162692464888096,
		  -0.054800182580947876,   0.052576825022697449,  -0.105417527258396149,   0.162675797939300537,   0.056601151823997498,   0.158786043524742126,   0.035719409584999084,   0.151088297367095947,  -0.046730849891901016,   0.025125732645392418,  -0.098689265549182892,   0.090225480496883392,  -0.016204178333282471,  -0.011277893558144569,  -0.004288212396204472,  -0.106375247240066528,   0.091792300343513489,   0.021168099716305733,   0.080974444746971130,   0.108282364904880524,   0.182037681341171265,   0.115767918527126312,   0.075922861695289612,   0.213375762104988098,   0.187276363372802734,   0.046226423233747482,   0.157502263784408569,   0.095375336706638336,   0.159428656101226807,   0.152365043759346008,   0.027073174715042114,   0.084024861454963684,   0.160055771470069885,  -0.009118296205997467,  -0.015525629743933678,   0.043245963752269745,  -0.045701920986175537,   0.004000536631792784,  -0.077700942754745483,  -0.060583084821701050,   0.079452730715274811,   0.104660324752330780,  -0.040037564933300018,  -0.017935806885361671,  -0.063678078353404999,   0.169098421931266785,   0.102114014327526093,  -0.034208342432975769,
	};
	static const double bias01[]=
	{
		   0.063564635813236237,   0.067409232258796692,   0.037352185696363449,  -0.010961883701384068,   0.139947131276130676,  -0.028114786371588707,   0.175188809633255005,   0.058169372379779816,  -0.046174276620149612,  -0.016275495290756226,  -0.028049904853105545,   0.003518699901178479,   0.010505917482078075,   0.093193583190441132,   0.083317227661609650,   0.020609689876437187,   0.126502379775047302,   0.063168078660964966,  -0.070523172616958618,   0.037394430488348007,   0.023194927722215652,  -0.036656770855188370,  -0.021291155368089676,   0.143945172429084778,   0.061271712183952332,   0.097141407430171967,   0.109715841710567474,   0.135792911052703857,   0.056323986500501633,  -0.110115796327590942,   0.079645358026027679,   0.111144594848155975,   0.006995297968387604,   0.090590417385101318,   0.161971926689147949,   0.092760242521762848,   0.072442762553691864,  -0.000755263434257358,  -0.011050751432776451,   0.094387151300907135,   0.094672694802284241,   0.094633229076862335,   0.160927250981330872,   0.044505130499601364,   0.057330582290887833,   0.008585708215832710,  -0.129347398877143860,   0.061406962573528290,
	};
	static const double weight02[]=
	{
		   0.021736858412623405,   0.099122770130634308,   0.091600112617015839,  -0.066293612122535706,   0.059440180659294128,   0.140994295477867126,   0.122788876295089722,   0.119795992970466614,  -0.159827694296836853,  -0.053958065807819366,  -0.141712397336959839,   0.070911750197410583,  -0.084270164370536804,  -0.021070554852485657,  -0.089736156165599823,  -0.084831081330776215,  -0.040014162659645081,  -0.103010028600692749,   0.136448696255683899,  -0.103033848106861115,  -0.003484718268737197,   0.051602471619844437,   0.083232894539833069,   0.011546194553375244,   0.005491736344993114,  -0.024919172748923302,   0.050235055387020111,  -0.042740199714899063,   0.110517062246799469,  -0.031764842569828033,   0.152497947216033936,  -0.005642423406243324,  -0.096169121563434601,   0.145846471190452576,  -0.048227138817310333,   0.056418154388666153,   0.033536255359649658,   0.094945773482322693,  -0.141848236322402954,  -0.101813562214374542,  -0.129187360405921936,  -0.100657686591148376,  -0.078268878161907196,  -0.012939292006194592,  -0.070663340389728546,  -0.002544535091146827,   0.003054152242839336,  -0.014059534296393394,
		   0.137066170573234558,   0.049393486231565475,  -0.106899760663509369,   0.137942016124725342,  -0.007097891066223383,   0.078240342438220978,   0.059466496109962463,   0.045502103865146637,   0.107722729444503784,  -0.033216979354619980,   0.043970953673124313,   0.104374118149280548,  -0.036255463957786560,   0.085871331393718719,   0.070392973721027374,   0.016900867223739624,   0.104659616947174072,   0.132927969098091125,   0.166561588644981384,   0.122240245342254639,  -0.067935444414615631,  -0.107977725565433502,  -0.033960413187742233,   0.046757660806179047,   0.027569623664021492,   0.054759357124567032,   0.057645019143819809,   0.117653720080852509,  -0.061572652310132980,   0.134963005781173706,  -0.118705272674560547,  -0.049178034067153931,   0.149928048253059387,   0.054908741265535355,   0.097327977418899536,   0.056169103831052780,   0.138392746448516846,  -0.048902895301580429,   0.073859997093677521,   0.082207389175891876,  -0.006018421147018671,   0.147343367338180542,  -0.047824520617723465,  -0.091583423316478729,   0.143404006958007812,  -0.041709441691637039,  -0.078271642327308655,   0.136077359318733215,
		   0.094821885228157043,  -0.054315838962793350,   0.069504663348197937,   0.137264221906661987,  -0.054970055818557739,   0.049673341214656830,   0.053639166057109833,  -0.051069006323814392,  -0.028713367879390717,   0.090959437191486359,   0.123489946126937866,  -0.081509783864021301,   0.126421615481376648,  -0.009016871452331543,   0.106591515243053436,   0.045395947992801666,  -0.064209587872028351,   0.058951023966073990,  -0.118036895990371704,   0.001126339891925454,  -0.095294699072837830,  -0.045338582247495651,   0.079394847154617310,   0.129984766244888306,  -0.007913288660347462,  -0.106923684477806091,   0.144719466567039490,   0.037423174828290939,  -0.078873500227928162,   0.032507233321666718,  -0.146299496293067932,   0.101247645914554596,  -0.018668640404939651,   0.040445487946271896,   0.074284799396991730,   0.088599570095539093,   0.065338172018527985,  -0.164217069745063782,  -0.037528112530708313,   0.089478380978107452,   0.031592082232236862,  -0.028700072318315506,  -0.006672603543847799,   0.091786235570907593,  -0.085942648351192474,  -0.140913367271423340,   0.116532497107982635,   0.118136428296566010,
		  -0.014862730167806149,  -0.037370800971984863,  -0.107672825455665588,   0.047724016010761261,   0.092713899910449982,   0.123125225305557251,   0.017029097303748131,   0.161576494574546814,   0.142378598451614380,   0.012018200941383839,  -0.010675545781850815,  -0.046050403267145157,   0.112062603235244751,  -0.132509008049964905,   0.153919830918312073,  -0.032525166869163513,   0.069992467761039734,   0.039566691964864731,   0.012539654038846493,   0.072832010686397552,  -0.034259732812643051,   0.061486836522817612,   0.166807755827903748,  -0.064685642719268799,   0.164109632372856140,  -0.090480946004390717,   0.036404643207788467,  -0.010263928212225437,  -0.060681965202093124,   0.089942574501037598,   0.115746960043907166,  -0.045776173472404480,  -0.000886630325112492,  -0.120686054229736328,   0.047963518649339676,  -0.013811235316097736,   0.097048468887805939,  -0.030192337930202484,   0.017036832869052887,   0.125185802578926086,   0.157987922430038452,  -0.111605033278465271,   0.036443695425987244,  -0.090135529637336731,   0.013930181041359901,  -0.095079571008682251,  -0.075333483517169952,  -0.078248672187328339,
		   0.163302361965179443,   0.011766158044338226,  -0.007406334858387709,   0.050766255706548691,  -0.070416048169136047,  -0.105082735419273376,   0.053200729191303253,  -0.084274634718894958,   0.131781071424484253,   0.008676655590534210,   0.102942049503326416,   0.139800906181335449,  -0.077830031514167786,  -0.049994204193353653,   0.112222746014595032,   0.139932215213775635,  -0.068107984960079193,  -0.059964887797832489,  -0.099080629646778107,   0.011003239080309868,   0.093102045357227325,   0.069414973258972168,   0.156696811318397522,  -0.142745032906532288,   0.184411883354187012,  -0.002436455804854631,   0.055786877870559692,  -0.083896785974502563,   0.111685194075107574,  -0.071269325911998749,   0.019704090431332588,   0.071582444012165070,   0.147093921899795532,   0.143638953566551208,   0.051520731300115585,   0.171527966856956482,   0.035916421562433243,   0.081358224153518677,   0.032632850110530853,   0.121372118592262268,   0.030238978564739227,  -0.154507353901863098,   0.041088648140430450,   0.126183301210403442,   0.028215117752552032,   0.068425036966800690,   0.116900876164436340,  -0.118539817631244659,
		  -0.056299801915884018,  -0.135868459939956665,   0.081707283854484558,  -0.075093805789947510,  -0.042675696313381195,   0.078780226409435272,   0.074539408087730408,  -0.074542231857776642,   0.084046870470046997,  -0.054149217903614044,   0.030073797330260277,   0.069056957960128784,  -0.122704803943634033,   0.124154031276702881,   0.080444894731044769,   0.030246859416365623,  -0.015670776367187500,  -0.078393839299678802,  -0.056270491331815720,  -0.108328394591808319,   0.130032777786254883,  -0.058690611273050308,  -0.030074877664446831,   0.131638452410697937,  -0.121398203074932098,   0.036448106169700623,   0.117305278778076172,   0.006626063026487827,   0.179943382740020752,   0.126004114747047424,   0.047379385679960251,   0.016854567453265190,   0.139969810843467712,   0.152886092662811279,  -0.035420771688222885,  -0.038952492177486420,   0.108661293983459473,  -0.067486502230167389,   0.056569810956716537,   0.162694677710533142,   0.089099533855915070,   0.097819522023200989,  -0.155239373445510864,   0.063174150884151459,   0.143071532249450684,  -0.064486429095268250,   0.057867418974637985,   0.142764881253242493,
		   0.015291833318769932,   0.129370585083961487,  -0.055602822452783585,  -0.038694232702255249,   0.165045261383056641,   0.116751238703727722,   0.143919542431831360,  -0.017844730988144875,  -0.041775643825531006,   0.014558739028871059,   0.126119539141654968,  -0.000024052831577137,   0.183618769049644470,   0.049567013978958130,  -0.055630940943956375,   0.055371228605508804,   0.107824601233005524,  -0.002184971002861857,   0.133749812841415405,   0.079020775854587555,   0.064959943294525146,   0.137308999896049500,   0.128214344382286072,   0.149249747395515442,  -0.106892831623554230,   0.069350309669971466,   0.117618151009082794,   0.118115052580833435,   0.195724189281463623,  -0.077577404677867889,   0.052987847477197647,   0.156387791037559509,  -0.058499950915575027,  -0.047684270888566971,  -0.138476967811584473,  -0.109493158757686615,   0.133612215518951416,  -0.065886862576007843,   0.131658792495727539,   0.152210056781768799,  -0.014894876629114151,   0.104321464896202087,  -0.084945321083068848,  -0.089855268597602844,  -0.080909281969070435,   0.104087598621845245,   0.133644744753837585,   0.048834722489118576,
		  -0.012759377248585224,  -0.154336050152778625,  -0.107653066515922546,   0.056800547987222672,   0.039438672363758087,  -0.048583842813968658,   0.073466189205646515,  -0.060276333242654800,   0.068053498864173889,  -0.091455996036529541,  -0.001971819903701544,  -0.127049922943115234,  -0.004878222942352295,  -0.065205462276935577,  -0.139822185039520264,  -0.027308931574225426,   0.102415904402732849,   0.006520491093397141,   0.171551316976547241,  -0.095768675208091736,   0.074442990124225616,  -0.001856723101809621,   0.030965283513069153,   0.111556179821491241,  -0.110567606985569000,  -0.131476208567619324,  -0.114671394228935242,   0.020641157403588295,  -0.073296800255775452,  -0.101814344525337219,   0.047761715948581696,  -0.138399004936218262,  -0.046883337199687958,  -0.096633799374103546,   0.075431622564792633,  -0.149727120995521545,   0.102398961782455444,  -0.096159338951110840,   0.116219319403171539,  -0.070040524005889893,  -0.142445728182792664,   0.172705277800559998,   0.055098444223403931,  -0.089037276804447174,   0.115184769034385681,   0.022561021149158478,   0.101618692278862000,  -0.128165006637573242,
		   0.060133621096611023,  -0.060921542346477509,   0.042587108910083771,  -0.020601056516170502,   0.032310456037521362,  -0.052196651697158813,   0.101273342967033386,  -0.021021256223320961,   0.166512966156005859,   0.009912014938890934,   0.009115694090723991,  -0.030002178624272346,   0.048162728548049927,  -0.033754829317331314,   0.130743354558944702,  -0.072610504925251007,  -0.155842855572700500,   0.087774515151977539,   0.064792558550834656,   0.004998633638024330,  -0.002019597683101892,   0.172350496053695679,   0.065219387412071228,  -0.022587710991501808,   0.049212720245122910,   0.088169068098068237,  -0.110973834991455078,   0.054328989237546921,  -0.119856886565685272,   0.139371037483215332,  -0.138576969504356384,   0.126107662916183472,   0.070251740515232086,  -0.021665122359991074,   0.096963360905647278,   0.139033719897270203,   0.056391414254903793,  -0.137937232851982117,  -0.024906020611524582,   0.027116673067212105,   0.062739498913288116,   0.040912043303251266,  -0.058244407176971436,  -0.115213766694068909,  -0.052671346813440323,   0.102708615362644196,  -0.070201188325881958,  -0.103796645998954773,
		  -0.048620183020830154,  -0.062928207218647003,  -0.032753579318523407,  -0.084496699273586273,   0.082019142806529999,   0.011949583888053894,   0.007156055886298418,   0.043100334703922272,   0.110938280820846558,   0.011686492711305618,  -0.052241310477256775,  -0.112727485597133636,  -0.033236533403396606,   0.214095801115036011,   0.071574591100215912,  -0.008756723254919052,  -0.080016709864139557,  -0.172299146652221680,  -0.064850978553295135,   0.098684206604957581,   0.117007948458194733,  -0.021530106663703918,  -0.124763913452625275,   0.165184587240219116,   0.056097824126482010,  -0.008841441012918949,   0.047671958804130554,  -0.104347787797451019,   0.048262421041727066,   0.003341061063110828,   0.112210035324096680,  -0.074423737823963165,  -0.145721748471260071,   0.058072153478860855,  -0.092178240418434143,  -0.076019026339054108,   0.172390535473823547,   0.027014110237360001,   0.156345337629318237,  -0.007422929629683495,   0.106953240931034088,   0.051656551659107208,  -0.155956298112869263,   0.020222431048750877,   0.183195233345031738,   0.115791074931621552,   0.064063444733619690,   0.041810415685176849,
		   0.062003389000892639,   0.099131524562835693,  -0.110108576714992523,  -0.068062327802181244,  -0.148219451308250427,  -0.026705823838710785,   0.119021020829677582,  -0.087286248803138733,  -0.051697961986064911,   0.126532286405563354,   0.083297826349735260,  -0.134586349129676819,  -0.121004670858383179,  -0.141142904758453369,  -0.072084575891494751,  -0.080385133624076843,  -0.017362510785460472,   0.080718919634819031,  -0.004254353232681751,  -0.131540760397911072,   0.094487294554710388,  -0.050190765410661697,  -0.139979198575019836,  -0.038119889795780182,  -0.053164850920438766,   0.013800510205328465,   0.033939696848392487,   0.104595959186553955,  -0.104777358472347260,  -0.070943146944046021,  -0.122872941195964813,  -0.103981673717498779,   0.024890726432204247,  -0.056443069130182266,   0.103778049349784851,  -0.043730575591325760,  -0.163710430264472961,  -0.098256602883338928,  -0.134536013007164001,  -0.120016574859619141,  -0.144732952117919922,  -0.115990847349166870,   0.133798345923423767,   0.025903346017003059,  -0.003199639031663537,  -0.169351294636726379,   0.044719018042087555,  -0.053888052701950073,
		   0.059427130967378616,   0.045746408402919769,   0.003577230032533407,  -0.065747678279876709,  -0.067269191145896912,  -0.068942993879318237,  -0.095777131617069244,  -0.115941211581230164,  -0.156061351299285889,   0.028765199705958366,  -0.015110683627426624,  -0.023185610771179199,  -0.152603641152381897,   0.102781862020492554,   0.015788771212100983,   0.194727465510368347,   0.182041659951210022,   0.023087004199624062,   0.177714869379997253,  -0.051141075789928436,   0.133152306079864502,  -0.135258600115776062,  -0.059493839740753174,   0.061799660325050354,  -0.168928220868110657,   0.062789879739284515,   0.036528468132019043,   0.116136737167835236,   0.080123379826545715,  -0.177334055304527283,  -0.109484016895294189,  -0.031509242951869965,  -0.065542064607143402,   0.134624242782592773,  -0.077601678669452667,   0.093587033450603485,   0.072497509419918060,  -0.124483399093151093,   0.052625719457864761,   0.102142632007598877,   0.014628246426582336,   0.210924625396728516,  -0.008456284180283546,   0.056899037212133408,   0.060771323740482330,   0.181751281023025513,  -0.032074663788080215,   0.076039671897888184,
		  -0.025933627039194107,   0.116258658468723297,   0.013795732520520687,  -0.048399195075035095,   0.165785580873489380,   0.015185052528977394,  -0.065583720803260803,  -0.048970524221658707,  -0.050569992512464523,   0.010271225124597549,   0.038362286984920502,   0.068862192332744598,   0.000397152500227094,   0.029918150976300240,  -0.006935693789273500,   0.012481871992349625,  -0.088999755680561066,   0.017280014231801033,  -0.109969288110733032,   0.074780844151973724,  -0.028609642758965492,   0.125038504600524902,  -0.076321572065353394,   0.030503647401928902,  -0.014201076701283455,  -0.091936208307743073,  -0.024151336401700974,   0.057754222303628922,  -0.035304952412843704,   0.044319178909063339,  -0.123816810548305511,  -0.124320514500141144,   0.001718317158520222,  -0.138063356280326843,  -0.056827723979949951,   0.077538736164569855,   0.070050932466983795,   0.002354839351028204,  -0.009285135194659233,   0.065011397004127502,  -0.037329964339733124,   0.020608482882380486,   0.081092908978462219,  -0.074064165353775024,  -0.038462523370981216,  -0.039211232215166092,   0.118791393935680389,   0.074337974190711975,
		  -0.056693434715270996,  -0.058460038155317307,  -0.087241910398006439,   0.070990972220897675,  -0.098309524357318878,   0.127556383609771729,  -0.019629040732979774,   0.127610743045806885,   0.111463516950607300,   0.065686956048011780,  -0.007429029792547226,  -0.077469274401664734,   0.054639324545860291,  -0.153378829360008240,   0.125815093517303467,  -0.112813711166381836,  -0.110933788120746613,   0.111965946853160858,  -0.048157799988985062,   0.082796484231948853,  -0.009418226778507233,  -0.031963191926479340,   0.130048945546150208,  -0.121316909790039062,   0.148455634713172913,   0.183157533407211304,  -0.096543185412883759,  -0.077993132174015045,   0.111817516386508942,   0.062903247773647308,  -0.020228894427418709,   0.155342891812324524,   0.172830507159233093,  -0.047021079808473587,   0.068567998707294464,   0.056761864572763443,   0.098930098116397858,   0.110688291490077972,   0.028593793511390686,   0.076028168201446533,   0.024306135252118111,  -0.054902877658605576,   0.166841208934783936,  -0.023188401013612747,  -0.041257385164499283,  -0.083843611180782318,  -0.079129211604595184,  -0.125289395451545715,
		   0.071542620658874512,   0.202242434024810791,   0.081044569611549377,   0.168119162321090698,   0.092103742063045502,   0.075089737772941589,  -0.059463188052177429,   0.134065866470336914,   0.082992449402809143,   0.037826962769031525,   0.000436290982179344,   0.128097087144851685,  -0.091000005602836609,  -0.172929331660270691,   0.118752948939800262,   0.103994697332382202,   0.013764837756752968,   0.084407836198806763,  -0.138668462634086609,   0.025751540437340736,  -0.016440363600850105,   0.038603086024522781,   0.166704148054122925,   0.124262645840644836,   0.145485028624534607,  -0.042103823274374008,   0.053815692663192749,   0.046378813683986664,   0.030477587133646011,   0.138656362891197205,   0.126349374651908875,   0.165593370795249939,  -0.055670775473117828,   0.045564465224742889,   0.028260419145226479,   0.136232078075408936,   0.000582432374358177,   0.003486717119812965,   0.007671022787690163,  -0.103183217346668243,   0.178885504603385925,  -0.141120538115501404,  -0.021236104890704155,  -0.105006322264671326,  -0.081925347447395325,  -0.005166931543499231,   0.021170649677515030,   0.159768611192703247,
		  -0.062160823494195938,  -0.021675392985343933,   0.054798830300569534,  -0.095060147345066071,   0.071062594652175903,  -0.135003969073295593,  -0.136811509728431702,   0.050428405404090881,   0.086867645382881165,  -0.114066600799560547,   0.124956041574478149,  -0.129713147878646851,  -0.090797178447246552,   0.091305032372474670,   0.003469971707090735,  -0.022663651034235954,  -0.104043722152709961,  -0.109328694641590118,  -0.113107621669769287,  -0.038257181644439697,  -0.023280400782823563,  -0.121398553252220154,   0.023119803518056870,   0.024320488795638084,   0.107449844479560852,  -0.008016701787710190,  -0.162500128149986267,  -0.146857365965843201,  -0.045538339763879776,  -0.067834012210369110,  -0.040400128811597824,  -0.094793722033500671,   0.003171788761392236,   0.026185762137174606,   0.053077749907970428,  -0.011257176287472248,  -0.173931002616882324,   0.031672257930040359,   0.117275536060333252,  -0.044351052492856979,   0.084851093590259552,  -0.071985736489295959,  -0.050386585295200348,  -0.145760744810104370,  -0.061990119516849518,   0.059284601360559464,   0.127843588590621948,  -0.103565387427806854,
		   0.045478161424398422,   0.133395433425903320,   0.131643325090408325,   0.128677397966384888,   0.037338070571422577,   0.035246875137090683,   0.035389006137847900,  -0.090009115636348724,   0.110401980578899384,  -0.030345747247338295,  -0.067336127161979675,   0.092436738312244415,  -0.006812208797782660,   0.031857095658779144,  -0.077871426939964294,  -0.137012526392936707,   0.098258286714553833,   0.182205319404602051,  -0.158688977360725403,  -0.040248058736324310,  -0.021843707188963890,  -0.046281240880489349,   0.051599405705928802,  -0.044907331466674805,   0.148147836327552795,   0.043471563607454300,   0.107473656535148621,   0.047000445425510406,   0.095051176846027374,  -0.121725477278232574,   0.115233898162841797,  -0.085323564708232880,   0.064650751650333405,  -0.119153961539268494,   0.097833216190338135,   0.120406918227672577,  -0.022025274112820625,  -0.084126651287078857,   0.100474871695041656,   0.159296452999114990,   0.161484286189079285,  -0.087368980050086975,   0.050729941576719284,  -0.091232985258102417,   0.139082938432693481,   0.082069702446460724,   0.017654851078987122,   0.062034092843532562,
		   0.164072290062904358,  -0.008815459907054901,   0.128574252128601074,   0.177377372980117798,   0.004708716180175543,   0.129659101366996765,   0.119507133960723877,   0.055194400250911713,  -0.068955138325691223,   0.039343185722827911,   0.127504378557205200,   0.155539542436599731,   0.089003808796405792,  -0.166803896427154541,  -0.000379174249246716,   0.027570735663175583,  -0.112735204398632050,   0.182497143745422363,   0.091704584658145905,   0.107272297143936157,  -0.015428734943270683,   0.108409360051155090,   0.108378097414970398,   0.072736971080303192,   0.127162948250770569,  -0.065742678940296173,  -0.107049778103828430,   0.089900717139244080,  -0.062906607985496521,  -0.004915283992886543,  -0.004815822932869196,   0.145956411957740784,  -0.084708042442798615,  -0.013614735566079617,  -0.027658855542540550,   0.022546239197254181,  -0.000592968252021819,   0.097363114356994629,   0.075388789176940918,   0.146658390760421753,   0.011756231077015400,  -0.085758745670318604,   0.157789185643196106,  -0.027882168069481850,  -0.052918888628482819,  -0.076114594936370850,   0.131088674068450928,  -0.020866064354777336,
		   0.141912087798118591,  -0.032320600003004074,   0.099015235900878906,  -0.068654038012027740,  -0.137599870562553406,  -0.053221940994262695,   0.082687176764011383,   0.054333485662937164,   0.170765578746795654,   0.031628407537937164,  -0.036189395934343338,  -0.078375384211540222,   0.153924897313117981,   0.038226295262575150,   0.087533310055732727,   0.091531477868556976,  -0.175449967384338379,  -0.107493139803409576,  -0.092299230396747589,   0.104830391705036163,  -0.083048708736896515,  -0.089773371815681458,  -0.062784790992736816,  -0.145179703831672668,   0.042965188622474670,   0.015512755140662193,  -0.023811249062418938,  -0.058117516338825226,  -0.052663479000329971,  -0.032353520393371582,  -0.050377860665321350,  -0.090494446456432343,   0.159007638692855835,   0.130102500319480896,  -0.074801571667194366,   0.136861905455589294,  -0.100483439862728119,   0.119999781250953674,  -0.019117644056677818,  -0.091277040541172028,   0.134392157196998596,  -0.018761273473501205,   0.139137804508209229,   0.057879570871591568,  -0.186391636729240417,  -0.060595713555812836,  -0.022811064496636391,  -0.056962188333272934,
		  -0.087416350841522217,   0.101140372455120087,  -0.055890701711177826,  -0.097803987562656403,   0.056795205920934677,  -0.112910404801368713,   0.029709473252296448,   0.018904708325862885,  -0.014630813151597977,  -0.131564989686012268,   0.096192739903926849,  -0.032375030219554901,  -0.051357310265302658,   0.078783132135868073,  -0.009693135507404804,  -0.063106037676334381,   0.023460835218429565,   0.002750257262960076,  -0.006449729204177856,  -0.119360685348510742,  -0.066334538161754608,   0.034687105566263199,   0.153566181659698486,  -0.044119596481323242,   0.110382154583930969,   0.011021154001355171,  -0.090780459344387054,   0.096795953810214996,   0.030451636761426926,   0.104565590620040894,  -0.095272123813629150,  -0.127224043011665344,  -0.039382554590702057,  -0.103097863495349884,  -0.041261415928602219,  -0.073526874184608459,  -0.041264608502388000,   0.070950649678707123,  -0.033030033111572266,   0.003873419249430299,   0.021590044721961021,   0.005294175818562508,  -0.046111989766359329,   0.092173762619495392,   0.055721268057823181,  -0.141965240240097046,   0.083419069647789001,  -0.120701447129249573,
		   0.074662961065769196,   0.033156689256429672,   0.084677211940288544,  -0.026673914864659309,  -0.058412648737430573,   0.081240184605121613,  -0.030425034463405609,  -0.007514216471463442,   0.142708659172058105,   0.029295250773429871,  -0.017276050522923470,  -0.030027922242879868,   0.092331767082214355,   0.010845246724784374,   0.031328361481428146,  -0.135379940271377563,  -0.036379031836986542,   0.132522553205490112,  -0.090743452310562134,   0.052638586610555649,   0.015877095982432365,   0.059559024870395660,   0.088000789284706116,   0.060590513050556183,   0.187503129243850708,  -0.088187910616397858,   0.045543391257524490,  -0.059985708445310593,   0.032830998301506042,   0.127657920122146606,   0.124404780566692352,   0.093065135180950165,   0.149862334132194519,   0.055591374635696411,   0.059082940220832825,   0.018484879285097122,   0.046746410429477692,   0.103104598820209503,  -0.048991605639457703,   0.146165758371353149,  -0.055103756487369537,  -0.113374449312686920,   0.108824186027050018,   0.111642673611640930,  -0.084510393440723419,  -0.062543608248233795,  -0.103980742394924164,   0.023966468870639801,
		  -0.040518764406442642,   0.045767389237880707,  -0.056939605623483658,   0.094880305230617523,   0.120520152151584625,  -0.042483441531658173,  -0.003524857107549906,   0.156065076589584351,   0.110414527356624603,  -0.037505716085433960,   0.192176029086112976,   0.096308037638664246,   0.181718885898590088,  -0.073318019509315491,  -0.026675181463360786,  -0.092522069811820984,  -0.076706804335117340,   0.073572680354118347,  -0.049350261688232422,   0.141149520874023438,   0.091094292700290680,   0.034242991358041763,   0.127393350005149841,   0.059798389673233032,   0.128875970840454102,   0.157671242952346802,   0.007816335186362267,   0.110007479786872864,   0.048778083175420761,  -0.024850247427821159,   0.042542088776826859,   0.106631755828857422,  -0.080858886241912842,  -0.135810866951942444,   0.018807357177138329,  -0.034796580672264099,  -0.068681605160236359,  -0.105054229497909546,  -0.085907809436321259,   0.078854314982891083,   0.022091591730713844,  -0.130779087543487549,  -0.076775237917900085,   0.119005955755710602,   0.130145087838172913,   0.031391873955726624,  -0.103573694825172424,  -0.018589373677968979,
		  -0.051552630960941315,  -0.146078884601593018,  -0.001195372664369643,  -0.073991961777210236,  -0.095248654484748840,  -0.073258168995380402,   0.151116698980331421,   0.180804714560508728,   0.057476148009300232,  -0.101906932890415192,   0.052410822361707687,   0.002184206154197454,  -0.080490149557590485,   0.218619197607040405,  -0.001672114245593548,  -0.074026443064212799,   0.187642201781272888,  -0.133371159434318542,   0.043884556740522385,   0.111144572496414185,   0.081695221364498138,  -0.043563328683376312,   0.029200579971075058,   0.152276948094367981,   0.008395408280193806,  -0.071583859622478485,   0.134030342102050781,   0.153605833649635315,  -0.011455978266894817,   0.074091434478759766,   0.200839340686798096,   0.041782610118389130,   0.061656706035137177,   0.170455992221832275,  -0.065861821174621582,  -0.142740353941917419,  -0.082127958536148071,   0.049900989979505539,  -0.080100536346435547,   0.024483889341354370,   0.071610808372497559,   0.120842657983303070,  -0.064878448843955994,   0.176285341382026672,   0.078416772186756134,   0.183029413223266602,  -0.055597282946109772,   0.065193630754947662,
		  -0.023453719913959503,   0.126305952668190002,   0.009835464879870415,   0.046441815793514252,   0.117618337273597717,   0.090531244874000549,   0.076373212039470673,  -0.027326485142111778,   0.132646113634109497,  -0.050626080483198166,   0.020996864885091782,   0.128351688385009766,   0.077016845345497131,  -0.015966065227985382,  -0.125091820955276489,   0.036031160503625870,   0.170289590954780579,  -0.032508049160242081,   0.024701755493879318,   0.135642051696777344,   0.148936912417411804,  -0.070990465581417084,  -0.136652901768684387,  -0.071566328406333923,   0.094564162194728851,  -0.017492661252617836,  -0.074556030333042145,   0.144021317362785339,   0.162017464637756348,  -0.008379638195037842,  -0.011421005241572857,  -0.109514102339744568,  -0.115900799632072449,   0.137612715363502502,   0.048606604337692261,   0.023850891739130020,  -0.010898274369537830,  -0.073020279407501221,  -0.049246665090322495,  -0.024639464914798737,  -0.127315655350685120,   0.114675022661685944,  -0.125785529613494873,   0.184167474508285522,   0.171641454100608826,   0.195288807153701782,  -0.028574122115969658,  -0.022379206493496895,
		  -0.035662554204463959,  -0.126788020133972168,  -0.092218071222305298,   0.114989593625068665,  -0.114526376128196716,  -0.168425604701042175,  -0.016637898981571198,  -0.033006388694047928,   0.108141638338565826,   0.032302405685186386,   0.142584547400474548,   0.039882015436887741,  -0.104837670922279358,  -0.033760655671358109,   0.057738792151212692,   0.042757816612720490,  -0.085136242210865021,  -0.137978583574295044,  -0.114381060004234314,   0.079578965902328491,   0.031431712210178375,   0.055932473391294479,   0.043222162872552872,   0.012578821741044521,  -0.067455887794494629,  -0.110020391643047333,  -0.088988542556762695,  -0.004311929922550917,  -0.122756361961364746,  -0.115770965814590454,  -0.152950510382652283,  -0.028423065319657326,   0.026599610224366188,  -0.066894501447677612,  -0.092722967267036438,   0.051670067012310028,  -0.095988087356090546,  -0.151659950613975525,  -0.005310312844812870,   0.085246466100215912,   0.043060444295406342,  -0.105312854051589966,   0.133727744221687317,   0.102603323757648468,   0.077961735427379608,  -0.103317230939865112,  -0.039435651153326035,  -0.003128608688712120,
		   0.147500529885292053,   0.181191056966781616,   0.062931686639785767,  -0.010370008647441864,  -0.043944261968135834,   0.116993799805641174,   0.109909005463123322,   0.118591092526912689,   0.038976900279521942,   0.041001237928867340,   0.005830244626849890,   0.084168799221515656,   0.112415291368961334,  -0.083256892859935760,   0.065575279295444489,   0.000821010908111930,  -0.081083409488201141,  -0.088449403643608093,   0.102253906428813934,   0.042499911040067673,  -0.096286900341510773,   0.112178437411785126,   0.113885216414928436,   0.129769131541252136,  -0.074453279376029968,   0.128340929746627808,  -0.054548162966966629,  -0.028226474300026894,  -0.051938738673925400,  -0.096836544573307037,  -0.102401264011859894,   0.103528290987014771,   0.076072320342063904,   0.124059773981571198,   0.004569165874272585,  -0.013857616111636162,  -0.090115897357463837,   0.075642794370651245,  -0.086326874792575836,   0.036579508334398270,  -0.049478944391012192,  -0.090894825756549835,   0.145599439740180969,   0.094828233122825623,  -0.114835746586322784,  -0.067983895540237427,   0.118564568459987640,   0.088024839758872986,
		   0.153094351291656494,  -0.066977076232433319,  -0.060904338955879211,  -0.063280038535594940,  -0.071629159152507782,  -0.097116932272911072,  -0.084626637399196625,   0.096851743757724762,  -0.153953075408935547,  -0.085827901959419250,  -0.088700070977210999,   0.015866488218307495,  -0.102317601442337036,   0.114224255084991455,  -0.104966185986995697,  -0.079865612089633942,   0.013280519284307957,   0.032503988593816757,   0.101258583366870880,   0.022447295486927032,   0.102127321064472198,  -0.129847407341003418,  -0.007612398825585842,  -0.080321460962295532,  -0.043026860803365707,  -0.032157827168703079,   0.110072843730449677,   0.050729509443044662,   0.127897307276725769,   0.060304399579763412,   0.088672734797000885,   0.003384658368304372,  -0.117672599852085114,  -0.028643390163779259,   0.088403657078742981,  -0.130891144275665283,  -0.020173480734229088,   0.056773360818624496,   0.117302715778350830,   0.158823028206825256,   0.003511594841256738,   0.037477545440196991,  -0.073633648455142975,   0.185416519641876221,  -0.081272833049297333,  -0.031788680702447891,   0.014046907424926758,   0.073349602520465851,
		   0.046568240970373154,   0.031867921352386475,   0.006894602440297604,  -0.029981477186083794,   0.100841395556926727,  -0.130321443080902100,   0.044376112520694733,  -0.088844433426856995,  -0.093400709331035614,  -0.045669011771678925,   0.088864192366600037,  -0.056085169315338135,  -0.035548642277717590,  -0.015949638560414314,  -0.116534419357776642,  -0.035577490925788879,  -0.029868286103010178,  -0.027195883914828300,   0.091716445982456207,  -0.057457052171230316,  -0.067032590508460999,  -0.107736483216285706,   0.004622668027877808,  -0.126347422599792480,  -0.089126288890838623,   0.120933219790458679,   0.117822632193565369,   0.034354273229837418,  -0.065986357629299164,   0.111199617385864258,  -0.112715288996696472,  -0.050309158861637115,  -0.008469704538583755,   0.036644425243139267,   0.049394443631172180,   0.031536903232336044,  -0.073091074824333191,   0.083771802484989166,   0.045396052300930023,  -0.018174052238464355,   0.102583348751068115,  -0.122102864086627960,  -0.057049218565225601,   0.100116126239299774,   0.063650466501712799,   0.031844738870859146,   0.157470792531967163,  -0.083323061466217041,
		   0.061648316681385040,   0.203427657485008240,   0.011119161732494831,   0.082663170993328094,  -0.023786867037415504,   0.130255520343780518,   0.089052177965641022,  -0.078992821276187897,   0.117664903402328491,  -0.061184458434581757,   0.072977431118488312,  -0.032273024320602417,   0.123666979372501373,  -0.157196819782257080,   0.127068385481834412,   0.103506349027156830,  -0.089373856782913208,   0.126256406307220459,  -0.023198749870061874,   0.154675647616386414,   0.081183671951293945,   0.079865187406539917,   0.008269230835139751,   0.082234233617782593,   0.192755892872810364,   0.107005558907985687,  -0.031885895878076553,  -0.022954724729061127,   0.115282721817493439,  -0.018985189497470856,   0.104947462677955627,  -0.054844547063112259,  -0.017083937302231789,   0.073509678244590759,   0.206424698233604431,  -0.077309750020503998,  -0.113259926438331604,   0.004401080310344696,  -0.121716156601905823,  -0.025146512314677238,   0.151612162590026855,  -0.061315655708312988,   0.156330451369285583,   0.129342824220657349,  -0.035394292324781418,  -0.025085724890232086,  -0.005144170951098204,   0.090926766395568848,
		   0.158200263977050781,  -0.071653828024864197,   0.036819599568843842,  -0.005277588032186031,   0.012907611206173897,  -0.052832860499620438,  -0.026066703721880913,  -0.042962394654750824,   0.075265318155288696,   0.033126298338174820,   0.100643955171108246,   0.110289640724658966,   0.028570484369993210,   0.152416110038757324,   0.055086661130189896,   0.049304872751235962,   0.155355617403984070,   0.133501648902893066,   0.014099300839006901,   0.106751412153244019,  -0.067163892090320587,   0.013618800789117813,   0.035199724137783051,   0.163617268204689026,  -0.108683943748474121,  -0.094478175044059753,  -0.047936294227838516,  -0.009566943161189556,   0.139960333704948425,  -0.059263873845338821,   0.006939884275197983,   0.081227399408817291,  -0.031562972813844681,   0.173727527260780334,   0.097455859184265137,  -0.120962128043174744,   0.171823173761367798,   0.005179543048143387,   0.187133654952049255,   0.143841981887817383,  -0.035489898175001144,  -0.065619729459285736,  -0.023249663412570953,   0.130491495132446289,  -0.067749731242656708,   0.111658953130245209,   0.062469858676195145,  -0.105936430394649506,
		   0.042724888771772385,   0.130960062146186829,   0.052014734596014023,  -0.061142489314079285,   0.005231459159404039,   0.123193196952342987,  -0.044409021735191345,   0.076355792582035065,   0.103721320629119873,   0.050722360610961914,   0.077370062470436096,  -0.129381179809570312,   0.106964774429798126,  -0.172438651323318481,   0.016363516449928284,  -0.049359351396560669,   0.072333663702011108,   0.114245161414146423,  -0.114489741623401642,   0.146459951996803284,  -0.087905995547771454,   0.080769903957843781,   0.024202818050980568,   0.127159133553504944,  -0.006405795924365520,   0.018073514103889465,  -0.083795182406902313,   0.027117623016238213,   0.110848896205425262,  -0.036576416343450546,   0.007108157500624657,   0.002951460424810648,   0.115049101412296295,  -0.021206025034189224,   0.118914179503917694,   0.037986494600772858,   0.049533195793628693,  -0.076202407479286194,  -0.124475568532943726,  -0.024593656882643700,   0.160151913762092590,   0.071405015885829926,   0.148673161864280701,   0.106981441378593445,  -0.124246798455715179,   0.074020698666572571,   0.108641088008880615,   0.041801255196332932,
		   0.072373285889625549,  -0.093241892755031586,   0.085053108632564545,   0.049395091831684113,  -0.031787842512130737,   0.062020771205425262,   0.108449198305606842,   0.120352536439895630,  -0.054813776165246964,   0.066801391541957855,   0.119267530739307404,  -0.048662573099136353,  -0.021769717335700989,  -0.043244577944278717,  -0.129766508936882019,   0.161182940006256104,   0.088876411318778992,   0.020867059007287025,   0.082660898566246033,   0.169390618801116943,   0.135088741779327393,  -0.119555242359638214,   0.002269951626658440,   0.123970329761505127,   0.133938297629356384,   0.082345344126224518,   0.188799858093261719,   0.071286581456661224,   0.062389042228460312,  -0.052860941737890244,   0.178622007369995117,  -0.070494748651981354,  -0.023528901860117912,   0.183161839842796326,   0.048308681696653366,  -0.086205199360847473,   0.089599289000034332,  -0.058202371001243591,   0.184488475322723389,   0.013291357085108757,  -0.090012602508068085,  -0.000739888113457710,  -0.129000797867774963,  -0.054675262421369553,  -0.075250104069709778,   0.108020305633544922,   0.119824141263961792,   0.099529482424259186,
		  -0.012535247020423412,   0.163118690252304077,   0.076753646135330200,   0.065292045474052429,  -0.094739340245723724,  -0.020589126273989677,  -0.093994311988353729,  -0.023359928280115128,  -0.010621094144880772,   0.118728198111057281,   0.070399478077888489,   0.027476852759718895,  -0.050687845796346664,   0.082461215555667877,   0.163438081741333008,  -0.109270840883255005,   0.012185396626591682,  -0.073984883725643158,   0.125425323843955994,   0.149408504366874695,  -0.098905175924301147,   0.053889553993940353,   0.110849983990192413,  -0.051860392093658447,   0.032549984753131866,   0.165526747703552246,   0.021338030695915222,  -0.013420595787465572,  -0.091020084917545319,   0.086427398025989532,  -0.057497821748256683,  -0.009037693031132221,   0.106808044016361237,   0.065983884036540985,   0.188168883323669434,   0.070588536560535431,   0.027254920452833176,  -0.108062975108623505,  -0.112302392721176147,   0.041684158146381378,   0.029838608577847481,   0.011391744017601013,   0.199594095349311829,   0.042093366384506226,  -0.140964508056640625,  -0.137049883604049683,  -0.047793544828891754,  -0.031161980703473091,
		   0.088282600045204163,   0.048867430537939072,  -0.103223003447055817,   0.108214974403381348,  -0.065456815063953400,   0.123345628380775452,  -0.047068171203136444,   0.127330973744392395,   0.181487321853637695,   0.054736398160457611,   0.076205223798751831,   0.022886993363499641,   0.078551694750785828,  -0.080745011568069458,   0.015980428084731102,  -0.042441818863153458,   0.064692310988903046,   0.135001033544540405,  -0.060458973050117493,   0.052939500659704208,   0.148882642388343811,   0.105010911822319031,  -0.035251051187515259,   0.015897762030363083,   0.062380775809288025,   0.164551421999931335,  -0.038536820560693741,  -0.066151916980743408,  -0.095490820705890656,   0.105880632996559143,  -0.061386473476886749,   0.000555851380340755,   0.111997500061988831,   0.119926609098911285,   0.086479075253009796,   0.168631851673126221,   0.004658945370465517,  -0.113958448171615601,  -0.020859757438302040,   0.116146594285964966,   0.140952825546264648,   0.119816772639751434,   0.178103968501091003,   0.129062354564666748,  -0.070757254958152771,  -0.130644679069519043,   0.000069082758272998,   0.071047715842723846,
		  -0.080125890672206879,   0.110303364694118500,   0.041769362986087799,   0.002564815571531653,   0.073350913822650909,   0.150279864668846130,   0.082050137221813202,   0.130592733621597290,  -0.017627449706196785,  -0.038270674645900726,   0.092327684164047241,   0.070035077631473541,   0.161027595400810242,  -0.026983750984072685,  -0.067744307219982147,   0.205760464072227478,   0.079022891819477081,  -0.035047542303800583,   0.062200669199228287,   0.162947639822959900,  -0.004936055745929480,  -0.046383976936340332,   0.096971139311790466,   0.093547374010086060,   0.094830073416233063,  -0.089301392436027527,  -0.008862929418683052,   0.099711664021015167,   0.071421295404434204,   0.033231422305107117,  -0.022194212302565575,   0.119881957769393921,   0.147413074970245361,   0.159736007452011108,  -0.138275668025016785,  -0.097639746963977814,   0.165217816829681396,   0.191247344017028809,   0.066075950860977173,   0.162475854158401489,  -0.094867862761020660,   0.153526201844215393,   0.098112970590591431,   0.004369939211755991,   0.174378171563148499,  -0.007171478588134050,  -0.133114844560623169,   0.153407216072082520,
		   0.014395096339285374,  -0.091582961380481720,   0.088591523468494415,  -0.174389079213142395,   0.121269933879375458,  -0.068314030766487122,   0.061457592993974686,  -0.050320859998464584,  -0.082397773861885071,  -0.051629349589347839,  -0.076110549271106720,  -0.014077816158533096,   0.077397726476192474,  -0.018251320347189903,  -0.011181509122252464,   0.115865454077720642,  -0.067922592163085938,  -0.005828076507896185,   0.201960414648056030,   0.041279695928096771,   0.007895087823271751,  -0.043106283992528915,  -0.116871416568756104,   0.067805424332618713,  -0.156529560685157776,   0.025174979120492935,   0.154727891087532043,   0.102192461490631104,   0.151989519596099854,  -0.112209156155586243,  -0.031884800642728806,  -0.117846868932247162,  -0.092644877731800079,  -0.049151003360748291,  -0.020071199163794518,  -0.030543150380253792,   0.060006119310855865,  -0.016938762739300728,   0.160796552896499634,   0.150791451334953308,   0.028348293155431747,   0.037686452269554138,   0.072878919541835785,  -0.083135850727558136,  -0.088645711541175842,  -0.070278950035572052,   0.109819613397121429,   0.144254535436630249,
		  -0.100588560104370117,   0.068142659962177277,   0.099822245538234711,   0.143313050270080566,  -0.093267984688282013,  -0.036330584436655045,  -0.007494651246815920,   0.138019829988479614,  -0.101717658340930939,   0.072683319449424744,   0.189138606190681458,  -0.082516260445117950,  -0.032811067998409271,  -0.153333902359008789,  -0.010524689219892025,  -0.081517145037651062,  -0.028066165745258331,  -0.011158691719174385,  -0.144707605242729187,  -0.045008301734924316,   0.008957443758845329,   0.038951165974140167,  -0.040852900594472885,   0.067477762699127197,   0.141359776258468628,   0.054413076490163803,  -0.063068948686122894,   0.034895658493041992,   0.049123350530862808,   0.013153210282325745,  -0.041591130197048187,   0.053298093378543854,  -0.115304984152317047,  -0.122602142393589020,   0.067560620605945587,   0.153892621397972107,  -0.078833930194377899,   0.105537220835685730,  -0.046549607068300247,  -0.033288680016994476,  -0.038216196000576019,  -0.063297219574451447,  -0.094297118484973907,   0.075539559125900269,  -0.140930205583572388,   0.017034111544489861,   0.026852803304791451,   0.155406504869461060,
		   0.038701903074979782,   0.122903749346733093,  -0.060909725725650787,  -0.123151198029518127,   0.002407208783552051,   0.136816173791885376,  -0.018100498244166374,  -0.024392364546656609,  -0.043390154838562012,  -0.131947934627532959,   0.044681377708911896,   0.097372658550739288,  -0.038306761533021927,  -0.007096272427588701,  -0.016504799947142601,  -0.111661046743392944,   0.005129255354404449,  -0.091946505010128021,   0.075018450617790222,   0.001299773342907429,  -0.149494886398315430,   0.010095108300447464,   0.023025153204798698,   0.065557412803173065,   0.099489793181419373,   0.074121952056884766,   0.071379281580448151,  -0.156576648354530334,   0.079051017761230469,  -0.115639314055442810,  -0.084843263030052185,  -0.115877613425254822,  -0.095634534955024719,  -0.082234680652618408,   0.050433643162250519,  -0.070192851126194000,   0.033093675971031189,   0.049247410148382187,   0.087953798472881317,   0.015255929902195930,  -0.040925610810518265,  -0.015674969181418419,   0.061591200530529022,  -0.040080275386571884,   0.068074509501457214,  -0.051921315491199493,  -0.049818880856037140,  -0.143834814429283142,
		  -0.083593241870403290,   0.012753333896398544,  -0.070677436888217926,  -0.072970800101757050,   0.155478864908218384,   0.135199353098869324,  -0.063375666737556458,  -0.116297587752342224,  -0.116274505853652954,  -0.148414507508277893,   0.003234316362068057,  -0.071203634142875671,   0.100913889706134796,  -0.047330331057310104,   0.107963860034942627,   0.116824500262737274,   0.118880763649940491,  -0.150686144828796387,   0.063782505691051483,   0.115007787942886353,   0.025299986824393272,  -0.000203554925974458,  -0.083410479128360748,   0.069445699453353882,   0.019240867346525192,   0.004589785821735859,   0.068782664835453033,   0.128958821296691895,   0.002282957080751657,   0.021392129361629486,   0.092101044952869415,   0.090411007404327393,   0.056899871677160263,  -0.087267592549324036,  -0.113557226955890656,   0.013460087589919567,  -0.012552570551633835,   0.163499101996421814,  -0.054914332926273346,  -0.053650047630071640,  -0.055484171956777573,  -0.058470357209444046,   0.015184605494141579,   0.148823082447052002,  -0.007182673551142216,   0.074190661311149597,  -0.156029298901557922,   0.024104747921228409,
		  -0.020581152290105820,   0.197213351726531982,   0.096786797046661377,  -0.071517966687679291,  -0.076632648706436157,   0.074929043650627136,  -0.119035907089710236,   0.111876688897609711,   0.071817182004451752,   0.127647295594215393,   0.102202445268630981,   0.140098005533218384,  -0.021060109138488770,  -0.086305499076843262,   0.181056126952171326,   0.074585452675819397,   0.125781610608100891,   0.167257368564605713,   0.003010526997968554,  -0.109924472868442535,   0.085977070033550262,   0.157796144485473633,   0.171338409185409546,   0.113239780068397522,   0.192757815122604370,   0.060339801013469696,   0.120260164141654968,   0.037216369062662125,   0.073570661246776581,  -0.076063014566898346,  -0.041103541851043701,  -0.014305714517831802,   0.136780917644500732,   0.010069795884191990,   0.125046640634536743,   0.188368678092956543,  -0.040767654776573181,   0.071296729147434235,   0.144027724862098694,   0.093959949910640717,   0.159011080861091614,  -0.039513479918241501,   0.120297700166702271,  -0.015154596418142319,  -0.051088966429233551,  -0.046445764601230621,  -0.084696806967258453,   0.000378126802388579,
		   0.150133028626441956,   0.054367341101169586,   0.105081804096698761,   0.117385089397430420,  -0.015743149444460869,  -0.007442950736731291,   0.112875029444694519,  -0.054641727358102798,   0.157739624381065369,  -0.089375369250774384,   0.006348524708300829,   0.038352545350790024,   0.139029949903488159,  -0.183320775628089905,  -0.013894585892558098,   0.121444106101989746,   0.055053584277629852,   0.107755750417709351,  -0.089988753199577332,  -0.026628777384757996,  -0.067351795732975006,   0.131420254707336426,   0.178958177566528320,  -0.054553169757127762,   0.122729331254959106,   0.110843278467655182,   0.015578997321426868,  -0.091252148151397705,   0.007443272974342108,   0.104836009442806244,  -0.119921430945396423,   0.124633193016052246,   0.149306789040565491,  -0.019961260259151459,   0.025390326976776123,   0.041604705154895782,   0.086882233619689941,   0.050660707056522369,  -0.005567304324358702,   0.056674566119909286,  -0.040365096181631088,  -0.080682158470153809,   0.041439559310674667,  -0.028561387211084366,  -0.018927453085780144,  -0.057370621711015701,  -0.011810083873569965,  -0.098370052874088287,
		  -0.133057013154029846,   0.022322831675410271,  -0.014476254582405090,  -0.162974789738655090,  -0.064005434513092041,   0.123865790665149689,   0.146997570991516113,  -0.053380742669105530,  -0.072316810488700867,  -0.022002872079610825,   0.074308574199676514,   0.133767619729042053,  -0.112926274538040161,   0.176847383379936218,   0.024952933192253113,   0.001055724103935063,   0.159753844141960144,  -0.037546649575233459,  -0.010796800255775452,  -0.025284612551331520,   0.089557208120822906,   0.092144861817359924,  -0.022511349990963936,  -0.114573933184146881,  -0.102260962128639221,   0.023443628102540970,  -0.019299641251564026,   0.055089041590690613,   0.099688939750194550,  -0.119307316839694977,   0.140595540404319763,   0.026326248422265053,   0.101178698241710663,  -0.058507304638624191,  -0.050240878015756607,   0.033793590962886810,   0.003285780083388090,   0.064588345587253571,  -0.077078342437744141,   0.163477078080177307,   0.004341075662523508,   0.202077537775039673,  -0.160235077142715454,   0.144233644008636475,   0.153261154890060425,  -0.054556094110012054,  -0.013645349070429802,   0.121513463556766510,
		   0.119640722870826721,  -0.016910502687096596,   0.155832841992378235,   0.017963888123631477,  -0.016051352024078369,   0.028636127710342407,  -0.002629640977829695,   0.026515439152717590,   0.065735869109630585,  -0.113474927842617035,  -0.073414914309978485,  -0.069345071911811829,  -0.133569046854972839,  -0.046603199094533920,  -0.081343442201614380,   0.126236826181411743,   0.034069936722517014,   0.057480782270431519,   0.054583903402090073,  -0.035121325403451920,  -0.053039185702800751,  -0.030667137354612350,  -0.072132535278797150,   0.025839034467935562,  -0.163972824811935425,   0.045333247631788254,  -0.004456780385226011,   0.119742512702941895,  -0.035684656351804733,   0.101292721927165985,  -0.068680666387081146,  -0.002202250296249986,  -0.012751637957990170,  -0.022734578698873520,  -0.141715511679649353,   0.101842172443866730,   0.090227983891963959,  -0.063707210123538971,   0.133568897843360901,  -0.003265073290094733,  -0.126628413796424866,   0.132028400897979736,   0.026658302173018456,   0.012531632557511330,   0.129011765122413635,   0.029679909348487854,  -0.158622652292251587,   0.117998808622360229,
		   0.128729268908500671,   0.185470417141914368,   0.053598504513502121,   0.037042472511529922,   0.028623618185520172,   0.016684025526046753,  -0.106663219630718231,   0.119513556361198425,  -0.059177778661251068,   0.068852990865707397,   0.099865779280662537,   0.155960783362388611,   0.013960661366581917,  -0.147011816501617432,   0.088386267423629761,   0.079904295504093170,  -0.020712364464998245,   0.084332011640071869,   0.097774408757686615,   0.024111857637763023,  -0.037644572556018829,  -0.038984205573797226,  -0.017310999333858490,   0.003735647303983569,  -0.050256147980690002,   0.116926543414592743,  -0.031845878809690475,   0.085697628557682037,  -0.071267917752265930,   0.069142967462539673,  -0.063991375267505646,   0.160395547747612000,  -0.077007092535495758,   0.016043242067098618,   0.183649241924285889,   0.119019977748394012,   0.073871038854122162,   0.004676378332078457,  -0.051928319036960602,   0.079019539058208466,   0.141107037663459778,   0.013402942568063736,  -0.014182006940245628,  -0.045671269297599792,  -0.025540415197610855,   0.006136605050414801,  -0.000554195721633732,   0.118988581001758575,
		   0.004525781609117985,   0.127192169427871704,   0.155890643596649170,   0.019405564293265343,   0.086723566055297852,   0.054832484573125839,  -0.005292054265737534,  -0.025157814845442772,   0.101369701325893402,  -0.086717881262302399,   0.038626287132501602,  -0.074691161513328552,  -0.096380881965160370,   0.100165493786334991,   0.140811085700988770,   0.051341190934181213,  -0.034941032528877258,   0.000659587036352605,   0.088581643998622894,  -0.016980128362774849,   0.135848373174667358,  -0.084513917565345764,   0.197331279516220093,   0.165879175066947937,   0.082604974508285522,   0.048505596816539764,   0.134819775819778442,   0.010213169269263744,  -0.098969079554080963,   0.137353897094726562,   0.083966046571731567,   0.129172012209892273,   0.178320497274398804,   0.112982340157032013,   0.057593472301959991,   0.118699230253696442,   0.031678039580583572,   0.014271513558924198,   0.104899227619171143,   0.095849789679050446,   0.056689690798521042,  -0.067455202341079712,   0.125876784324645996,  -0.015083364211022854,  -0.006959040183573961,   0.001925845281220973,  -0.019195344299077988,   0.088564649224281311,
		  -0.056279674172401428,  -0.072540648281574249,   0.094035588204860687,  -0.077233560383319855,   0.115749195218086243,  -0.100925028324127197,  -0.062894582748413086,   0.037267055362462997,   0.012712524272501469,  -0.011719022877514362,   0.090725295245647430,  -0.086033977568149567,  -0.124453604221343994,  -0.039866745471954346,   0.001402800902724266,   0.197184100747108459,  -0.015743723139166832,   0.000180082890437916,   0.143450438976287842,  -0.096462950110435486,  -0.058039024472236633,  -0.114658221602439880,  -0.017603537067770958,  -0.089935950934886932,  -0.179478913545608521,  -0.094848319888114929,   0.181576550006866455,  -0.014333460479974747,   0.007982471026480198,  -0.134813934564590454,   0.136681735515594482,  -0.048387177288532257,  -0.001953432103618979,  -0.087839655578136444,  -0.037348426878452301,  -0.132411539554595947,   0.149497881531715393,  -0.087109237909317017,   0.149267956614494324,  -0.066819995641708374,  -0.065558582544326782,   0.073566004633903503,  -0.041244249790906906,  -0.014168784022331238,  -0.016530416905879974,  -0.012740063481032848,  -0.112514965236186981,   0.070930317044258118,
		  -0.063961304724216461,  -0.054155718535184860,   0.150990098714828491,  -0.041567556560039520,   0.129075795412063599,   0.023418406024575233,  -0.039496783167123795,  -0.099623568356037140,  -0.024024115875363350,  -0.031699746847152710,  -0.031792178750038147,  -0.026539953425526619,   0.018954543396830559,   0.008632137440145016,  -0.102532044053077698,   0.151881843805313110,   0.132733508944511414,   0.019379403442144394,   0.075408264994621277,   0.121372677385807037,  -0.060424581170082092,  -0.001693759695626795,  -0.038057215511798859,   0.138995736837387085,   0.084418870508670807,   0.159429311752319336,   0.048734981566667557,   0.126483529806137085,  -0.068526551127433777,   0.130399450659751892,   0.052215881645679474,   0.126754119992256165,   0.162492528557777405,  -0.016447037458419800,   0.106459632515907288,   0.042268507182598114,  -0.043900970369577408,   0.135932311415672302,   0.097676977515220642,  -0.076327532529830933,   0.163047790527343750,   0.083698056638240814,   0.123256951570510864,  -0.087053030729293823,   0.134639754891395569,   0.126548156142234802,   0.064611241221427917,   0.021332323551177979,
		  -0.025451542809605598,   0.157218039035797119,   0.131727576255798340,   0.064501158893108368,  -0.012583516538143158,   0.137575909495353699,   0.094247750937938690,   0.119529113173484802,   0.197894349694252014,   0.159850671887397766,   0.124997235834598541,  -0.077490918338298798,  -0.079444117844104767,   0.068290889263153076,   0.066626101732254028,  -0.059429794549942017,   0.010027109645307064,   0.037692453712224960,  -0.002973406109958887,   0.125901550054550171,   0.149203851819038391,   0.021727742627263069,  -0.043946854770183563,  -0.025664724409580231,   0.013498177751898766,   0.037385720759630203,  -0.019555635750293732,   0.118356786668300629,  -0.105633832514286041,  -0.081143997609615326,   0.081220246851444244,   0.034163799136877060,   0.123800367116928101,   0.039260473102331161,   0.033530890941619873,  -0.083510652184486389,  -0.141808584332466125,  -0.045554462820291519,   0.013933402486145496,   0.093241542577743530,   0.115672692656517029,  -0.067914314568042755,   0.101173028349876404,   0.063926510512828827,  -0.084209516644477844,   0.031826008111238480,   0.125986054539680481,  -0.044656410813331604,
	};
	static const double bias02[]=
	{
		   0.092503815889358521,   0.110817253589630127,  -0.094807319343090057,   0.038542266935110092,   0.058064926415681839,   0.012316204607486725,   0.169875353574752808,  -0.125764042139053345,  -0.053020026534795761,  -0.074508287012577057,  -0.111305549740791321,  -0.027911540120840073,   0.119888074696063995,   0.123668506741523743,   0.001274883281439543,  -0.003981476649641991,  -0.010160743258893490,   0.122578248381614685,  -0.027237903326749802,  -0.034415848553180695,   0.049136605113744736,  -0.103235833346843719,   0.060335390269756317,  -0.092366285622119904,  -0.051461189985275269,   0.106756590306758881,   0.173072502017021179,  -0.007278684992343187,   0.091738186776638031,   0.045597802847623825,   0.061411429196596146,  -0.017699711024761200,   0.033835228532552719,   0.066294029355049133,   0.185743764042854309,   0.174056276679039001,  -0.037593793123960495,  -0.126262918114662170,   0.079915106296539307,  -0.090734057128429413,   0.056512758135795593,   0.037018384784460068,   0.137123003602027893,  -0.018218450248241425,  -0.009204948320984840,   0.055353101342916489,   0.183282494544982910,  -0.107207253575325012,
	};
	static const double weight03[]=
	{
		  -0.118858225643634796,  -0.098265610635280609,   0.118601344525814056,   0.138904839754104614,  -0.089409917593002319,  -0.017497487366199493,  -0.122421100735664368,   0.008285209536552429,   0.110105529427528381,   0.000713838962838054,   0.126317501068115234,   0.090169988572597504,  -0.096290729939937592,   0.121892012655735016,  -0.093194499611854553,   0.135112538933753967,   0.083800613880157471,  -0.112078420817852020,   0.116521500051021576,   0.004884067457169294,   0.044702388346195221,   0.083450064063072205,  -0.166212186217308044,   0.058734107762575150,  -0.109541796147823334,  -0.036527443677186966,  -0.116845250129699707,   0.017412714660167694,  -0.011931382119655609,   0.003254310227930546,   0.011804530397057533,  -0.079602815210819244,   0.034991066902875900,  -0.119245871901512146,  -0.142009958624839783,  -0.066534072160720825,   0.141070991754531860,   0.044223770499229431,  -0.093710727989673615,  -0.106339201331138611,   0.000824069953523576,   0.054906684905290604,  -0.159939989447593689,  -0.017373081296682358,   0.080630101263523102,   0.056171245872974396,  -0.033967033028602600,  -0.023293275386095047,
		  -0.024138962849974632,   0.112243935465812683,   0.115342274308204651,  -0.091998927295207977,   0.138315856456756592,  -0.065822094678878784,  -0.016026401892304420,  -0.086770191788673401,   0.132442474365234375,  -0.083554781973361969,   0.014656553976237774,  -0.111081913113594055,   0.034195445477962494,  -0.000577764702029526,   0.100774407386779785,   0.086062848567962646,  -0.025773873552680016,   0.164487257599830627,  -0.049312058836221695,   0.102796487510204315,   0.059802196919918060,   0.180273070931434631,   0.005619590170681477,  -0.091531276702880859,  -0.002827435964718461,   0.153472632169723511,  -0.028023775666952133,  -0.015150585211813450,   0.083718642592430115,   0.143965110182762146,   0.044634845107793808,   0.047432053834199905,   0.075999952852725983,   0.002546216128394008,   0.085810631513595581,   0.093690223991870880,   0.126038104295730591,  -0.130521968007087708,   0.008166913874447346,  -0.048206459730863571,   0.107730850577354431,  -0.159148246049880981,   0.059303566813468933,   0.103092469274997711,  -0.045469157397747040,   0.082245968282222748,  -0.043541904538869858,   0.123454935848712921,
		   0.034282013773918152,  -0.037553492933511734,  -0.049058265984058380,  -0.014701751060783863,   0.052571892738342285,  -0.126043513417243958,  -0.010858167894184589,  -0.021302562206983566,   0.044201344251632690,   0.007265300955623388,  -0.051879316568374634,   0.054728727787733078,   0.019839484244585037,   0.087489292025566101,  -0.044823911041021347,   0.026934007182717323,   0.095976702868938446,   0.150132820010185242,  -0.001195951714180410,  -0.009520934894680977,   0.137877404689788818,   0.000581682892516255,   0.043174996972084045,   0.016280353069305420,   0.129716873168945312,  -0.036475952714681625,  -0.087961979210376740,  -0.056838590651750565,   0.150415807962417603,  -0.012022993527352810,  -0.072765856981277466,   0.076036781072616577,  -0.062489017844200134,   0.119448728859424591,   0.020617000758647919,  -0.105453923344612122,   0.054168447852134705,  -0.137852534651756287,   0.019927566871047020,   0.170007452368736267,   0.073674865067005157,   0.008125003427267075,  -0.057642549276351929,   0.136062368750572205,  -0.015709638595581055,  -0.152957648038864136,  -0.030330115929245949,   0.150423660874366760,
		   0.066872775554656982,   0.051392406225204468,  -0.073799878358840942,   0.124825023114681244,   0.125677078962326050,  -0.029352260753512383,   0.187772512435913086,   0.018017623573541641,  -0.082041583955287933,   0.045134395360946655,  -0.053120788186788559,   0.045048493891954422,   0.038534447550773621,   0.063583455979824066,  -0.105314247310161591,  -0.156507194042205811,  -0.096436306834220886,   0.027850275859236717,  -0.013792547397315502,  -0.120603524148464203,   0.099172756075859070,   0.076927125453948975,   0.186185136437416077,   0.144390583038330078,  -0.183971405029296875,   0.149609580636024475,  -0.027105966582894325,   0.115843392908573151,   0.070963181555271149,   0.120557866990566254,   0.109299667179584503,   0.189370438456535339,   0.086077116429805756,   0.103346109390258789,  -0.077534198760986328,   0.167075455188751221,   0.058885734528303146,   0.021050399169325829,  -0.086592055857181549,  -0.072816461324691772,   0.082528389990329742,  -0.034813247621059418,   0.105773143470287323,  -0.084139958024024963,   0.148024573922157288,   0.164236664772033691,   0.150891304016113281,  -0.032826315611600876,
		   0.019051212817430496,   0.027180284261703491,   0.031201843172311783,  -0.052629373967647552,  -0.020689178258180618,  -0.127154186367988586,   0.037673227488994598,  -0.110625118017196655,  -0.041105572134256363,   0.087534032762050629,   0.019190780818462372,  -0.030417162925004959,  -0.033706311136484146,   0.112522967159748077,   0.101677849888801575,  -0.145556524395942688,  -0.016875252127647400,   0.136806294322013855,   0.096606247127056122,   0.125091150403022766,   0.129836946725845337,   0.121497653424739838,  -0.061616744846105576,  -0.040326137095689774,  -0.083433471620082855,   0.169374927878379822,  -0.087398529052734375,   0.047202635556459427,   0.166487097740173340,   0.033854141831398010,   0.042763561010360718,   0.160186350345611572,   0.048770267516374588,   0.126647263765335083,   0.049366652965545654,  -0.145790621638298035,   0.074639245867729187,  -0.158909156918525696,   0.009204921312630177,   0.182516247034072876,   0.041949529200792313,  -0.014607148244976997,  -0.144755616784095764,   0.156196177005767822,  -0.042019382119178772,  -0.025588089600205421,   0.004631350748240948,   0.136134669184684753,
		  -0.047890659421682358,   0.018513318151235580,  -0.099158175289630890,   0.021729527041316032,  -0.050798136740922928,   0.176480710506439209,   0.170989394187927246,   0.050381321460008621,   0.014034494757652283,   0.001999766798689961,  -0.175665006041526794,   0.153168261051177979,  -0.102728486061096191,  -0.075217641890048981,  -0.128837555646896362,   0.048441991209983826,   0.004491808824241161,  -0.079932607710361481,   0.137249991297721863,  -0.094683587551116943,   0.101422354578971863,  -0.064846321940422058,   0.103028275072574615,   0.140298470854759216,  -0.043964728713035583,  -0.115608617663383484,   0.115536451339721680,  -0.111625693738460541,  -0.016362046822905540,   0.032860998064279556,   0.070879988372325897,  -0.080478161573410034,  -0.001425647642463446,   0.018286615610122681,   0.108274713158607483,  -0.000189325306564569,   0.107480101287364960,   0.042951453477144241,   0.155231997370719910,   0.032117888331413269,   0.022920425981283188,   0.124518342316150665,   0.035474948585033417,   0.030167669057846069,  -0.086965024471282959,   0.129837110638618469,   0.000686835730448365,   0.133376628160476685,
		   0.118063442409038544,  -0.038546510040760040,   0.084053777158260345,  -0.037679970264434814,   0.066283620893955231,   0.063299015164375305,  -0.125366985797882080,  -0.080825276672840118,   0.119126014411449432,  -0.010214258916676044,   0.116148836910724640,  -0.000679858610965312,  -0.130274847149848938,  -0.075626470148563385,   0.107650786638259888,  -0.124123543500900269,  -0.026616659015417099,   0.018822973594069481,  -0.101021699607372284,  -0.039151847362518311,   0.022896038368344307,   0.139462664723396301,  -0.091721937060356140,  -0.022019702941179276,  -0.103265099227428436,   0.155381008982658386,   0.077770322561264038,  -0.056584335863590240,  -0.115021698176860809,  -0.120065733790397644,  -0.029753964394330978,   0.066844776272773743,  -0.016948485746979713,   0.064986877143383026,  -0.008560709655284882,  -0.140499293804168701,   0.120648436248302460,   0.063453473150730133,  -0.039562586694955826,   0.051290590316057205,  -0.028835205361247063,  -0.121223509311676025,   0.110830709338188171,   0.105200059711933136,  -0.072786010801792145,  -0.132696002721786499,   0.055435720831155777,  -0.004004133865237236,
		   0.066244333982467651,  -0.114477582275867462,  -0.045724153518676758,   0.021624663844704628,  -0.036170911043882370,   0.025357764214277267,   0.082445830106735229,   0.051771510392427444,  -0.054031684994697571,  -0.058460790663957596,   0.121334686875343323,   0.024844385683536530,   0.096239902079105377,  -0.084109410643577576,  -0.149162784218788147,   0.033805958926677704,   0.018149815499782562,   0.047023802995681763,   0.047023944556713104,  -0.082520641386508942,  -0.095497027039527893,   0.100782327353954315,  -0.019976379349827766,   0.039904568344354630,  -0.041457239538431168,   0.066846728324890137,  -0.021358858793973923,  -0.068305961787700653,   0.038044080138206482,   0.023092493414878845,  -0.141713351011276245,   0.068973556160926819,   0.037632372230291367,   0.106208302080631256,  -0.137713834643363953,  -0.068276561796665192,   0.024539928883314133,   0.020133186131715775,   0.054678779095411301,  -0.090063877403736115,  -0.107050135731697083,  -0.044503193348646164,  -0.091663070023059845,  -0.071652427315711975,  -0.033489763736724854,  -0.066936694085597992,   0.053552806377410889,   0.031816925853490829,
		   0.100773021578788757,  -0.096353888511657715,  -0.048028901219367981,  -0.128248542547225952,  -0.024486901238560677,  -0.091730944812297821,  -0.036557294428348541,   0.082625374197959900,  -0.088649123907089233,  -0.050206549465656281,  -0.032498214393854141,   0.058253146708011627,   0.106807410717010498,  -0.004817463457584381,  -0.145222544670104980,   0.019499324262142181,  -0.051089283078908920,  -0.147170737385749817,  -0.089513525366783142,  -0.069569088518619537,  -0.031987778842449188,   0.003803249215707183,  -0.014445699751377106,  -0.134231165051460266,   0.017502356320619583,   0.124773345887660980,  -0.046274591237306595,   0.128497675061225891,   0.073034569621086121,   0.006329496856778860,  -0.118594028055667877,   0.107575006783008575,  -0.002633274067193270,  -0.139393150806427002,   0.085104480385780334,  -0.092635110020637512,   0.099944055080413818,   0.074924595654010773,  -0.015933409333229065,  -0.023525852710008621,  -0.155728265643119812,  -0.094587646424770355,   0.041472997516393661,   0.033292379230260849,   0.102777540683746338,  -0.044549908488988876,  -0.138693034648895264,  -0.038134183734655380,
		  -0.106869719922542572,   0.062925554811954498,   0.028888598084449768,   0.158820912241935730,   0.156934112310409546,   0.094450168311595917,   0.048622060567140579,   0.028296159580349922,   0.062803477048873901,  -0.027599859982728958,   0.116768524050712585,   0.075693018734455109,   0.160700887441635132,   0.128295242786407471,   0.142522007226943970,  -0.006626693997532129,   0.109130702912807465,   0.127499237656593323,   0.052875034511089325,   0.107626095414161682,   0.090249516069889069,  -0.037782419472932816,   0.113366432487964630,   0.084373548626899719,  -0.108954071998596191,   0.033985871821641922,  -0.195165425539016724,  -0.015267977491021156,   0.186009958386421204,   0.071925669908523560,   0.082417778670787811,   0.009806028567254543,  -0.056462299078702927,   0.067805632948875427,  -0.107505552470684052,  -0.052714034914970398,   0.047589834779500961,  -0.020552955567836761,   0.088856518268585205,   0.188068687915802002,   0.171878442168235779,   0.007237679325044155,  -0.087934754788875580,   0.128292441368103027,   0.046646337956190109,  -0.145932793617248535,   0.127077594399452209,   0.035192646086215973,
		  -0.066357254981994629,  -0.062568649649620056,  -0.102712251245975494,  -0.005964828655123711,  -0.151634484529495239,  -0.127746850252151489,  -0.101653575897216797,   0.065644294023513794,  -0.110306583344936371,  -0.051160205155611038,  -0.013468121178448200,  -0.095408596098423004,  -0.042549438774585724,   0.041974350810050964,  -0.034914199262857437,   0.000247208459768444,  -0.004265398252755404,  -0.121507883071899414,  -0.028623301535844803,   0.015613347291946411,  -0.003606733866035938,  -0.118736848235130310,  -0.021412694826722145,  -0.151537507772445679,  -0.009999140165746212,   0.100276350975036621,   0.036229085177183151,   0.116265371441841125,  -0.017323099076747894,   0.022562710568308830,   0.112148128449916840,  -0.124710157513618469,  -0.078489422798156738,  -0.063084691762924194,  -0.110943496227264404,  -0.031575549393892288,  -0.049363452941179276,  -0.039308000355958939,   0.019880460575222969,  -0.107067637145519257,  -0.089358113706111908,  -0.061671655625104904,   0.089996039867401123,  -0.108643844723701477,   0.095539823174476624,  -0.123921371996402740,  -0.144279897212982178,  -0.027565311640501022,
		  -0.018977543339133263,  -0.064438112080097198,   0.145277515053749084,   0.046122118830680847,  -0.083585008978843689,   0.115668341517448425,  -0.114793755114078522,   0.000603314256295562,  -0.021747745573520660,  -0.047958586364984512,  -0.054566673934459686,  -0.036574255675077438,   0.000101561876363121,  -0.002595787635073066,   0.140945494174957275,   0.127372786402702332,   0.093182198703289032,   0.063397713005542755,  -0.057138305157423019,  -0.008045251481235027,   0.066344060003757477,   0.105671130120754242,   0.126445814967155457,   0.083189822733402252,   0.060151968151330948,   0.012749845162034035,  -0.038400281220674515,   0.001121015986427665,  -0.052540950477123260,  -0.076829493045806885,  -0.018803317099809647,   0.035918839275836945,   0.089700534939765930,  -0.026975670829415321,   0.139984771609306335,   0.061756987124681473,   0.034232813864946365,   0.004289026372134686,  -0.111758947372436523,   0.092654190957546234,  -0.065031021833419800,  -0.171387255191802979,   0.066954247653484344,   0.077298112213611603,  -0.016379861161112785,  -0.101324923336505890,   0.039088182151317596,   0.086074352264404297,
		   0.116859123110771179,   0.171180203557014465,   0.033493898808956146,  -0.134735777974128723,   0.089433513581752777,   0.000145310375955887,   0.174442097544670105,   0.027852715924382210,  -0.104048088192939758,   0.190816983580589294,  -0.071820430457592010,  -0.034318938851356506,   0.076328709721565247,  -0.063232697546482086,   0.046326208859682083,   0.063974805176258087,   0.123139634728431702,  -0.102702595293521881,  -0.031357694417238235,  -0.110305242240428925,   0.014815013855695724,   0.032191388309001923,   0.211460843682289124,   0.127757936716079712,   0.008113534189760685,   0.090361654758453369,   0.056932259351015091,  -0.106215551495552063,  -0.106462247669696808,  -0.044176049530506134,  -0.144765391945838928,   0.080586723983287811,  -0.035509496927261353,  -0.101777002215385437,   0.173257708549499512,   0.090640559792518616,  -0.116469249129295349,  -0.054759345948696136,   0.161163583397865295,   0.031464800238609314,   0.098654896020889282,   0.168899789452552795,  -0.059572864323854446,   0.081391550600528717,  -0.029527004808187485,   0.185082346200942993,   0.079639256000518799,   0.037894871085882187,
		  -0.007364722434431314,   0.157088339328765869,  -0.063459992408752441,   0.183153361082077026,   0.125885546207427979,  -0.063919648528099060,   0.156940430402755737,   0.081178165972232819,  -0.007654115557670593,  -0.159124970436096191,  -0.003872174769639969,  -0.037196401506662369,  -0.084713056683540344,   0.127530485391616821,   0.146899774670600891,  -0.047810934484004974,   0.125047490000724792,   0.171675086021423340,   0.066081658005714417,   0.036772161722183228,   0.012228908017277718,   0.104334808886051178,  -0.117767974734306335,   0.009632426314055920,   0.023825010284781456,   0.094252802431583405,   0.040615659207105637,   0.068250566720962524,   0.176876664161682129,   0.076827697455883026,   0.071327477693557739,  -0.046419683843851089,  -0.060568209737539291,   0.193082943558692932,  -0.073288775980472565,   0.006420057266950607,   0.113047055900096893,   0.098151497542858124,   0.052650999277830124,  -0.012243683449923992,   0.065239563584327698,  -0.156080454587936401,  -0.102303370833396912,   0.159814938902854919,   0.166564255952835083,   0.107613682746887207,  -0.027613971382379532,   0.134230792522430420,
		  -0.086141303181648254,  -0.051386203616857529,  -0.094429560005664825,  -0.017085233703255653,   0.052153579890727997,   0.054713390767574310,   0.074216231703758240,   0.067934937775135040,  -0.166501581668853760,   0.092803180217742920,   0.011086901649832726,  -0.029286805540323257,   0.086148075759410858,   0.116278961300849915,  -0.065913461148738861,   0.116078369319438934,  -0.170307233929634094,  -0.137559667229652405,  -0.152792468667030334,  -0.072487190365791321,   0.039761673659086227,   0.080445237457752228,  -0.071310073137283325,   0.156058907508850098,  -0.156939119100570679,  -0.060341510921716690,  -0.066570155322551727,   0.032732535153627396,  -0.171054720878601074,   0.107346236705780029,  -0.061759155243635178,   0.036800462752580643,   0.084910020232200623,  -0.186991378664970398,  -0.046383339911699295,   0.157453656196594238,  -0.064475044608116150,  -0.125018283724784851,   0.000115729730168823,  -0.038570918142795563,   0.021536320447921753,   0.047884121537208557,   0.110674597322940826,  -0.062134459614753723,   0.050593063235282898,   0.011493787169456482,   0.149654671549797058,  -0.168133467435836792,
		  -0.024279981851577759,   0.088414043188095093,  -0.003136228304356337,   0.190452709794044495,   0.131105780601501465,   0.059290543198585510,   0.111708410084247589,   0.002019237494096160,  -0.016293322667479515,  -0.048926427960395813,   0.136568412184715271,   0.063690215349197388,   0.057455576956272125,  -0.005128564778715372,   0.184973672032356262,   0.040374442934989929,  -0.106007829308509827,   0.153705045580863953,   0.131588503718376160,   0.090416893362998962,   0.003437038743868470,  -0.032574612647294998,   0.078855067491531372,  -0.159618213772773743,   0.061576072126626968,   0.074874557554721832,   0.029931420460343361,   0.108958177268505096,   0.152059152722358704,  -0.130664810538291931,   0.122357875108718872,   0.002105223247781396,   0.183230400085449219,   0.194548174738883972,  -0.028913425281643867,  -0.133118435740470886,   0.024044513702392578,   0.035089481621980667,  -0.156936958432197571,   0.189439877867698669,  -0.019076971337199211,  -0.058715458959341049,   0.121492981910705566,   0.005713727790862322,  -0.053252212703227997,  -0.144314721226692200,   0.005609277635812759,   0.026798030361533165,
		   0.042884506285190582,  -0.051736235618591309,  -0.131861537694931030,   0.117530681192874908,   0.132288277149200439,   0.096564404666423798,   0.126207038760185242,   0.092120178043842316,  -0.183779150247573853,   0.031989853829145432,   0.015771584585309029,   0.128296703100204468,  -0.078426681458950043,   0.152152195572853088,  -0.100305803120136261,  -0.026661846786737442,  -0.094096824526786804,  -0.022553959861397743,   0.039834130555391312,  -0.018022662028670311,   0.153359755873680115,   0.164696574211120605,   0.027749769389629364,   0.111140020191669464,  -0.134758695960044861,   0.140849053859710693,   0.008621393702924252,   0.102516554296016693,   0.061159558594226837,   0.072118058800697327,  -0.118830032646656036,   0.152923092246055603,  -0.023703867569565773,  -0.033088605850934982,   0.186914369463920593,   0.064779937267303467,   0.085845135152339935,  -0.092735767364501953,  -0.038919702172279358,  -0.011123368516564369,  -0.108497254550457001,   0.029375307261943817,  -0.022481361404061317,   0.065455213189125061,  -0.004523326177150011,  -0.041701395064592361,  -0.083414673805236816,   0.046448014676570892,
		   0.031077010557055473,  -0.048195872455835342,  -0.089825570583343506,   0.031388420611619949,   0.120927825570106506,   0.079884439706802368,  -0.087532557547092438,   0.055056288838386536,  -0.089955300092697144,   0.020666753873229027,  -0.150525301694869995,   0.110912218689918518,  -0.003715147497132421,  -0.028732672333717346,  -0.016543265432119370,  -0.107845202088356018,  -0.110225327312946320,  -0.009324378333985806,  -0.001752485171891749,  -0.024165168404579163,   0.035486474633216858,  -0.128260254859924316,   0.147404372692108154,   0.075725689530372620,  -0.000513895938638598,   0.011611755006015301,   0.140489801764488220,  -0.125759094953536987,   0.067066550254821777,   0.131197988986968994,  -0.026953635737299919,   0.130897074937820435,  -0.149961978197097778,  -0.100967325270175934,   0.148202404379844666,   0.128394067287445068,  -0.002156746108084917,   0.020158424973487854,   0.078947201371192932,  -0.019899573177099228,  -0.055972788482904434,   0.070852011442184448,   0.109135225415229797,   0.054158203303813934,  -0.082648105919361115,   0.005536227021366358,   0.005208370275795460,  -0.140715375542640686,
		  -0.135215982794761658,   0.147204145789146423,   0.114606380462646484,   0.191517889499664307,  -0.025546630844473839,  -0.079602353274822235,   0.165924012660980225,   0.063199989497661591,  -0.079952932894229889,   0.035843722522258759,  -0.158803910017013550,  -0.004833552055060863,   0.130818113684654236,   0.020995395258069038,   0.076660782098770142,   0.076871968805789948,  -0.049295950680971146,   0.085536584258079529,   0.154516950249671936,   0.021512955427169800,  -0.084697782993316650,  -0.096522204577922821,  -0.011227158829569817,  -0.032604731619358063,  -0.150437310338020325,  -0.043489765375852585,  -0.015348915010690689,   0.002204267540946603,   0.084091842174530029,   0.055294260382652283,  -0.002490265062078834,   0.019322754815220833,  -0.002951988251879811,   0.070636078715324402,  -0.000564477813895792,  -0.068158768117427826,  -0.132869258522987366,  -0.014364124275743961,   0.129910618066787720,   0.120143353939056396,  -0.008501426316797733,  -0.071268871426582336,   0.095261827111244202,   0.010434887371957302,   0.097026906907558441,   0.145956560969352722,   0.131959289312362671,   0.150464653968811035,
		  -0.031815405935049057,  -0.129107803106307983,  -0.122790724039077759,  -0.054628215730190277,  -0.068782478570938110,  -0.173679903149604797,  -0.013880729675292969,   0.085736900568008423,  -0.101783253252506256,   0.078201852738857269,   0.148112505674362183,  -0.164396166801452637,  -0.131315082311630249,  -0.061166778206825256,  -0.116641186177730560,   0.134443804621696472,  -0.034348450601100922,  -0.055051099509000778,   0.092657990753650665,  -0.060288760811090469,  -0.075981475412845612,  -0.089186705648899078,   0.025331728160381317,  -0.013787901960313320,  -0.062238171696662903,  -0.013793321326375008,  -0.100515142083168030,  -0.086066454648971558,  -0.111350625753402710,  -0.059855416417121887,  -0.066215232014656067,  -0.159343585371971130,  -0.074473910033702850,   0.082770273089408875,   0.100894726812839508,   0.095265828073024750,  -0.008172573521733284,   0.117047183215618134,  -0.029169268906116486,  -0.003448265371844172,  -0.079788334667682648,   0.095126122236251831,  -0.129984304308891296,   0.071252509951591492,   0.102403670549392700,   0.017688389867544174,  -0.061609517782926559,   0.064599588513374329,
		  -0.114915601909160614,   0.151102617383003235,   0.074179247021675110,   0.080593533813953400,   0.032000843435525894,  -0.016018938273191452,   0.054669573903083801,   0.164667412638664246,   0.152750059962272644,   0.158128574490547180,   0.088001623749732971,   0.120389498770236969,  -0.084583900868892670,   0.005551029928028584,  -0.044234074652194977,  -0.071714565157890320,   0.137727648019790649,   0.030764561146497726,  -0.035653866827487946,   0.089959621429443359,  -0.013472629711031914,  -0.024877807125449181,   0.042259145528078079,   0.167027965188026428,  -0.119325332343578339,   0.113485418260097504,   0.047574896365404129,  -0.082703061401844025,   0.171613723039627075,   0.051218755543231964,  -0.033768408000469208,   0.034904386848211288,  -0.085420548915863037,  -0.017542798072099686,   0.163781806826591492,   0.083211429417133331,  -0.055681686848402023,   0.023784374818205833,  -0.001596280490048230,   0.097468465566635132,   0.154869124293327332,   0.111817359924316406,   0.041581023484468460,   0.151476904749870300,   0.146859511733055115,  -0.025321498513221741,   0.011013976298272610,  -0.010788596235215664,
		   0.061342559754848480,   0.073894806206226349,  -0.197857379913330078,  -0.068186953663825989,  -0.134600132703781128,   0.185392171144485474,   0.042681358754634857,   0.020454268902540207,  -0.192612618207931519,   0.107074260711669922,   0.061808861792087555,   0.149016708135604858,  -0.126177802681922913,  -0.086960665881633759,  -0.043848887085914612,  -0.039976175874471664,  -0.105268515646457672,   0.131378039717674255,  -0.099043175578117371,  -0.094174750149250031,  -0.098739966750144958,  -0.049046419560909271,   0.071389988064765930,   0.049246497452259064,  -0.100791163742542267,  -0.060325268656015396,   0.157057553529739380,   0.102930903434753418,  -0.060259662568569183,   0.124467916786670685,  -0.052157331258058548,   0.136247992515563965,  -0.023173978552222252,   0.082305811345577240,   0.134131297469139099,   0.046676155179738998,   0.042055059224367142,   0.000031038558518048,   0.137496605515480042,  -0.024207172915339470,   0.001407817588187754,   0.031215915456414223,   0.125360608100891113,  -0.136033654212951660,   0.095500186085700989,   0.067848086357116699,   0.021612687036395073,   0.001372527331113815,
		  -0.068936750292778015,   0.136141389608383179,   0.007433269172906876,   0.000840918393805623,   0.026752371340990067,   0.100884579122066498,  -0.005992290563881397,   0.066440224647521973,   0.050964411348104477,   0.060291934758424759,  -0.149137169122695923,   0.021599484607577324,  -0.029754422605037689,   0.026525460183620453,  -0.010601003654301167,  -0.139084205031394958,   0.060276940464973450,  -0.134770169854164124,   0.000530285120476037,  -0.082121297717094421,  -0.079279772937297821,   0.065613731741905212,  -0.074928931891918182,   0.005118009634315968,   0.055513657629489899,   0.019203225150704384,   0.137205451726913452,  -0.107867494225502014,  -0.044067591428756714,   0.093952409923076630,  -0.069061316549777985,   0.118689268827438354,   0.088165313005447388,   0.023821279406547546,   0.076734535396099091,  -0.088020890951156616,  -0.161894783377647400,  -0.087938211858272552,   0.160686358809471130,  -0.089712925255298615,   0.004217817913740873,   0.071950331330299377,   0.013722599484026432,   0.026239383965730667,  -0.078715726733207703,  -0.008145166561007500,   0.133195236325263977,  -0.055107515305280685,
		  -0.100412108004093170,   0.164978668093681335,  -0.005153121426701546,   0.183496445417404175,  -0.066174305975437164,  -0.084514483809471130,  -0.022982904687523842,  -0.069249622523784637,   0.027427081018686295,   0.115000128746032715,  -0.100271031260490417,  -0.144585058093070984,   0.095019832253456116,   0.093193754553794861,   0.104854375123977661,   0.028779093176126480,   0.104293830692768097,   0.124234952032566071,  -0.019007049500942230,  -0.147869154810905457,   0.119800686836242676,   0.001239208970218897,   0.096996225416660309,  -0.017833603546023369,  -0.043519955128431320,   0.018398977816104889,  -0.172138929367065430,   0.006528103258460760,   0.117187999188899994,  -0.115606650710105896,  -0.007175992242991924,   0.027075964957475662,   0.175182223320007324,   0.072154298424720764,   0.050333723425865173,  -0.142471835017204285,   0.102370895445346832,  -0.075402952730655670,   0.074287146329879761,   0.169484272599220276,   0.104224406182765961,  -0.066001199185848236,  -0.122388295829296112,   0.197447195649147034,   0.003017096081748605,  -0.005039311014115810,   0.146142572164535522,   0.168268665671348572,
		   0.022188771516084671,  -0.159123331308364868,  -0.021585121750831604,  -0.074966393411159515,  -0.012547320686280727,   0.104364186525344849,  -0.120906203985214233,  -0.004589081276208162,   0.104731306433677673,   0.126174911856651306,  -0.105066642165184021,  -0.028906766325235367,   0.031344823539257050,  -0.085509121417999268,  -0.100253120064735413,  -0.119342796504497528,  -0.010904333554208279,  -0.127310127019882202,   0.034581478685140610,   0.123478733003139496,  -0.075911849737167358,  -0.164864793419837952,  -0.059487946331501007,   0.059403195977210999,  -0.046143092215061188,  -0.159846752882003784,  -0.118610560894012451,  -0.068505622446537018,   0.118765406310558319,  -0.040895342826843262,   0.098730184137821198,  -0.027859374880790710,   0.116923302412033081,   0.063121266663074493,   0.019788831472396851,   0.026799920946359634,  -0.001296248752623796,  -0.001614106586202979,  -0.003878734074532986,  -0.133999511599540710,   0.005679793655872345,   0.064976058900356293,  -0.057166919112205505,  -0.132132887840270996,  -0.001488990732468665,   0.086779825389385223,  -0.087316162884235382,  -0.090932719409465790,
		  -0.164624497294425964,  -0.107551477849483490,   0.090749479830265045,   0.097202219069004059,  -0.090730048716068268,  -0.071312978863716125,  -0.191093981266021729,  -0.066540598869323730,   0.155139729380607605,  -0.184074342250823975,   0.101736441254615784,  -0.123874001204967499,  -0.036861356347799301,  -0.062700979411602020,  -0.122135482728481293,  -0.129630386829376221,   0.030006283894181252,  -0.151828780770301819,   0.038079410791397095,  -0.006290134042501450,  -0.112292513251304626,   0.132710576057434082,   0.022683670744299889,  -0.096812322735786438,   0.133928999304771423,  -0.049105904996395111,  -0.084724947810173035,  -0.058666538447141647,  -0.111867636442184448,   0.004915479570627213,   0.026126120239496231,  -0.033043753355741501,   0.015102328732609749,   0.090117022395133972,  -0.045904599130153656,  -0.148381963372230530,  -0.008145692758262157,   0.088726446032524109,  -0.029209155589342117,  -0.086604654788970947,  -0.119847275316715240,  -0.181248858571052551,   0.040285266935825348,   0.038882892578840256,  -0.068134367465972900,  -0.066829867660999298,  -0.179444998502731323,   0.027660062536597252,
		   0.065998770296573639,   0.071587838232517242,   0.162911668419837952,   0.108717255294322968,  -0.054821390658617020,   0.104780666530132294,  -0.016947165131568909,  -0.048903983086347580,   0.003046893514692783,  -0.026949310675263405,  -0.119941003620624542,  -0.124184057116508484,  -0.094867639243602753,   0.146862506866455078,  -0.060325972735881805,  -0.134926483035087585,  -0.009931481443345547,  -0.048644158989191055,  -0.101832918822765350,  -0.124637372791767120,   0.103632561862468719,  -0.074992515146732330,  -0.116876237094402313,   0.100957356393337250,  -0.007147786673158407,  -0.079794079065322876,  -0.089686937630176544,  -0.126671686768531799,   0.139286503195762634,  -0.059359587728977203,   0.171086519956588745,  -0.039237312972545624,  -0.073994629085063934,   0.113082483410835266,  -0.087818130850791931,   0.058733507990837097,   0.077797651290893555,  -0.030675141140818596,  -0.056706909090280533,   0.179296657443046570,  -0.025687996298074722,  -0.070341140031814575,   0.070691712200641632,   0.168154984712600708,  -0.038793906569480896,  -0.107801213860511780,   0.071720764040946960,   0.059696942567825317,
		   0.129318907856941223,   0.059395134449005127,  -0.083249464631080627,  -0.093792580068111420,  -0.089064873754978180,  -0.072277739644050598,  -0.021958665922284126,   0.110528081655502319,   0.053851269185543060,   0.026243953034281731,  -0.067534737288951874,   0.168567523360252380,  -0.078772731125354767,   0.046037863940000534,   0.103676915168762207,   0.060671873390674591,   0.087829157710075378,  -0.091347031295299530,   0.117341026663780212,   0.091809920966625214,   0.167871654033660889,  -0.108733355998992920,   0.088237352669239044,   0.138907715678215027,  -0.112828910350799561,   0.156687572598457336,   0.180666267871856689,  -0.061599437147378922,  -0.010558802634477615,   0.103394396603107452,  -0.071386575698852539,  -0.090396150946617126,   0.038962483406066895,   0.018466971814632416,   0.024660306051373482,   0.151118114590644836,  -0.074942834675312042,  -0.077802456915378571,  -0.003516761120408773,  -0.001765047549270093,   0.082974627614021301,   0.101859629154205322,   0.076812334358692169,  -0.138870403170585632,  -0.008388678543269634,  -0.030048256739974022,  -0.084524191915988922,  -0.040793698281049728,
		  -0.046811640262603760,   0.000826931674964726,   0.065179213881492615,  -0.101196050643920898,  -0.012802887707948685,  -0.115112110972404480,   0.004302073735743761,   0.060467429459095001,  -0.027167981490492821,   0.046351265162229538,   0.062486384063959122,   0.014243475161492825,  -0.076616361737251282,   0.089617021381855011,   0.120259188115596771,  -0.055896509438753128,   0.023789841681718826,  -0.027755275368690491,  -0.081538796424865723,   0.121721409261226654,   0.112378202378749847,   0.087422877550125122,  -0.078002907335758209,  -0.120139651000499725,  -0.054116245359182358,  -0.034979362040758133,   0.060900025069713593,   0.075268931686878204,  -0.097984112799167633,  -0.033208653330802917,  -0.045758634805679321,   0.060922365635633469,   0.050309840589761734,   0.144503936171531677,   0.173430055379867554,   0.135402560234069824,   0.098998188972473145,   0.031574632972478867,   0.110514268279075623,   0.166989773511886597,   0.045667473226785660,  -0.061000723391771317,   0.118820875883102417,  -0.044481150805950165,   0.190895080566406250,   0.037395827472209930,   0.189460933208465576,  -0.086051702499389648,
		  -0.015885651111602783,  -0.131197005510330200,  -0.037247974425554276,   0.100131526589393616,  -0.011162091046571732,  -0.048647031188011169,  -0.086564563214778900,   0.085022136569023132,  -0.021591525524854660,  -0.048282101750373840,  -0.028090860694646835,  -0.067265003919601440,   0.131809979677200317,   0.029879821464419365,  -0.178099200129508972,  -0.011894293129444122,   0.052817996591329575,  -0.054493431001901627,   0.077375888824462891,   0.073581114411354065,  -0.028240555897355080,  -0.073044970631599426,  -0.093726262450218201,   0.069628231227397919,   0.079353399574756622,   0.072296924889087677,   0.133222967386245728,  -0.087285213172435760,   0.090678952634334564,   0.015534141100943089,  -0.079961329698562622,  -0.008471424691379070,  -0.096582897007465363,   0.100998915731906891,  -0.165057778358459473,  -0.044585268944501877,   0.101070068776607513,  -0.098304636776447296,   0.115071386098861694,  -0.071351356804370880,  -0.132419243454933167,  -0.050152517855167389,   0.060793507844209671,  -0.128133535385131836,   0.000921675760764629,   0.047210589051246643,   0.067340500652790070,  -0.039710447192192078,
		   0.075135946273803711,   0.035726539790630341,  -0.022088166326284409,   0.003442282322794199,  -0.002745133824646473,   0.053521677851676941,   0.094892203807830811,  -0.093320757150650024,  -0.088853418827056885,  -0.116073369979858398,   0.034482367336750031,  -0.131073668599128723,  -0.038837160915136337,   0.082407556474208832,   0.112934850156307220,  -0.003238211851567030,  -0.058227289468050003,   0.058297548443078995,  -0.001431413693353534,  -0.009775532409548759,  -0.047339212149381638,   0.090956658124923706,   0.089418828487396240,  -0.061324436217546463,  -0.003999263979494572,   0.058032661676406860,   0.100766137242317200,  -0.089726693928241730,   0.126191452145576477,  -0.136642083525657654,  -0.076829142868518829,  -0.040299575775861740,   0.059732384979724884,  -0.027479119598865509,  -0.036449067294597626,   0.046259701251983643,  -0.021579938009381294,  -0.012783559039235115,  -0.119057439267635345,   0.024567048996686935,  -0.067597806453704834,   0.075535401701927185,  -0.131501808762550354,  -0.039766442030668259,   0.120734490454196930,  -0.138554751873016357,   0.047693479806184769,   0.050681956112384796,
		  -0.011924248188734055,  -0.086622923612594604,   0.147269621491432190,   0.086010582745075226,   0.165012210607528687,   0.026620877906680107,   0.110209211707115173,  -0.097863174974918365,   0.098202049732208252,   0.004890603013336658,  -0.123181894421577454,  -0.056320473551750183,  -0.024599045515060425,  -0.057441163808107376,   0.073068082332611084,   0.112926855683326721,   0.110242009162902832,   0.154754996299743652,   0.002030853647738695,  -0.121303968131542206,   0.116209886968135834,   0.158540651202201843,  -0.043064992874860764,   0.121521055698394775,  -0.056157913058996201,   0.059528443962335587,   0.046882532536983490,   0.101166002452373505,   0.038436170667409897,   0.122601665556430817,   0.051345031708478928,  -0.011716201901435852,   0.059327431023120880,   0.077379688620567322,   0.116632662713527679,  -0.117451749742031097,   0.085617229342460632,  -0.080946721136569977,  -0.096054241061210632,   0.178116559982299805,   0.143133431673049927,   0.018951373174786568,   0.034181877970695496,   0.048238579183816910,   0.030032753944396973,  -0.085591495037078857,  -0.009808089584112167,   0.101430445909500122,
		   0.029238490387797356,  -0.018853349611163139,  -0.068047456443309784,   0.124578364193439484,   0.116096653044223785,   0.069457575678825378,   0.104552686214447021,  -0.096414320170879364,  -0.074204057455062866,  -0.105988457798957825,  -0.047535728663206100,  -0.058649551123380661,   0.113330677151679993,   0.157616510987281799,   0.018600363284349442,  -0.109153397381305695,   0.094043694436550140,  -0.112752564251422882,  -0.092656925320625305,   0.135487556457519531,   0.071890071034431458,   0.186605677008628845,   0.072861209511756897,  -0.090557664632797241,   0.155263185501098633,  -0.059672284871339798,  -0.179849430918693542,  -0.131459921598434448,   0.077876709401607513,  -0.100787594914436340,  -0.032348722219467163,   0.030501646921038628,  -0.017801370471715927,   0.105641938745975494,  -0.084186471998691559,  -0.107885748147964478,  -0.107363417744636536,   0.034182094037532806,  -0.088074304163455963,   0.017435112968087196,   0.008316743187606335,   0.065649755299091339,  -0.087305158376693726,   0.019773114472627640,  -0.081914752721786499,  -0.059896834194660187,  -0.100607886910438538,   0.147654742002487183,
		   0.180521711707115173,  -0.063777044415473938,   0.082971356809139252,   0.114962935447692871,   0.007427326403558254,   0.126614198088645935,   0.020572301000356674,   0.079581290483474731,   0.058690007776021957,   0.124580308794975281,  -0.009767518378794193,  -0.144130095839500427,  -0.011363750323653221,   0.044413276016712189,   0.014155262149870396,   0.106367260217666626,  -0.080261133611202240,  -0.127342224121093750,   0.089084804058074951,  -0.031316280364990234,   0.037081778049468994,  -0.012684042565524578,   0.064001768827438354,   0.118733502924442291,  -0.133446976542472839,   0.024676071479916573,   0.028379939496517181,   0.041553292423486710,  -0.138775765895843506,   0.064378768205642700,  -0.034552197903394699,  -0.106674991548061371,   0.027859359979629517,   0.077295199036598206,   0.080795191228389740,  -0.004701354540884495,   0.043134681880474091,  -0.053751159459352493,  -0.009276463650166988,  -0.024891039356589317,  -0.030209030956029892,  -0.051856644451618195,   0.005510922987014055,  -0.097591087222099304,   0.039175216108560562,   0.052102264016866684,   0.064901985228061676,  -0.093065127730369568,
		  -0.158421710133552551,  -0.090993411839008331,   0.007917124778032303,  -0.058341082185506821,   0.080763645470142365,  -0.013057100586593151,  -0.007427508477121592,   0.094715014100074768,  -0.038383644074201584,  -0.026546407490968704,   0.090184375643730164,  -0.119215600192546844,   0.086416587233543396,   0.150637909770011902,  -0.009405431337654591,  -0.097062431275844574,  -0.038111586123704910,  -0.061190009117126465,   0.126544803380966187,   0.015977974981069565,   0.167340129613876343,  -0.093161866068840027,  -0.126956358551979065,   0.050884716212749481,   0.059593785554170609,   0.174308419227600098,  -0.133399888873100281,  -0.019434301182627678,   0.020492333918809891,  -0.051850467920303345,  -0.046024661511182785,  -0.062485288828611374,   0.168412610888481140,   0.109670147299766541,   0.080451287329196930,  -0.094135567545890808,  -0.048385661095380783,   0.045091029256582260,  -0.135076776146888733,   0.010228274390101433,  -0.060250494629144669,  -0.026783036068081856,  -0.122775800526142120,   0.039558008313179016,  -0.063988298177719116,  -0.077985122799873352,   0.152183681726455688,   0.155222475528717041,
		  -0.054273623973131180,   0.176783546805381775,  -0.093418486416339874,   0.172422692179679871,  -0.079290226101875305,   0.096404910087585449,   0.006925460416823626,  -0.006358842831104994,  -0.100108511745929718,  -0.007163254078477621,   0.031732957810163498,   0.068844996392726898,   0.119298934936523438,   0.119355179369449615,   0.066676564514636993,   0.089751720428466797,   0.076040826737880707,   0.004075190983712673,   0.049128144979476929,   0.076634600758552551,   0.093328058719635010,   0.006798856426030397,   0.055224720388650894,   0.108811333775520325,  -0.097199849784374237,   0.026391293853521347,  -0.105446740984916687,   0.014256768859922886,  -0.041907690465450287,   0.170658677816390991,   0.093611985445022583,   0.012823427096009254,   0.102499790489673615,   0.035581652075052261,   0.118441112339496613,   0.014015648514032364,   0.095070526003837585,   0.041843056678771973,  -0.108163870871067047,  -0.012722668237984180,   0.143963679671287537,  -0.090294316411018372,  -0.093221396207809448,   0.006609574891626835,   0.056219857186079025,   0.072737649083137512,  -0.065499551594257355,   0.044151451438665390,
		  -0.149076923727989197,   0.067799858748912811,  -0.079515852034091949,   0.068670034408569336,   0.157791927456855774,  -0.125952646136283875,  -0.048872161656618118,  -0.091110795736312866,  -0.099622562527656555,  -0.033501461148262024,   0.107585176825523376,  -0.132139384746551514,  -0.033213581889867783,   0.113260388374328613,   0.128043830394744873,   0.098411120474338531,   0.116232126951217651,   0.081576362252235413,  -0.074666120111942291,   0.061707463115453720,  -0.026864653453230858,   0.181499287486076355,   0.084201000630855560,  -0.018417971208691597,  -0.053526915609836578,   0.069810204207897186,  -0.052518419921398163,   0.008734368719160557,   0.031374320387840271,  -0.098048307001590729,   0.000778178102336824,   0.054094649851322174,   0.166076272726058960,  -0.029456282034516335,   0.039428539574146271,   0.039595942944288254,  -0.011612620204687119,   0.046028289943933487,  -0.135842859745025635,  -0.067719936370849609,   0.091704882681369781,   0.048361539840698242,   0.131688430905342102,  -0.038869202136993408,   0.051757939159870148,  -0.004269655328243971,  -0.057584397494792938,  -0.082267135381698608,
		  -0.099162280559539795,   0.079197540879249573,   0.029054779559373856,   0.060947727411985397,   0.034205447882413864,  -0.024388104677200317,   0.103774085640907288,   0.074915498495101929,   0.108789816498756409,   0.059869918972253799,  -0.030101092532277107,  -0.121377632021903992,  -0.063937321305274963,   0.069535173475742340,   0.167495116591453552,  -0.050089344382286072,   0.124858744442462921,   0.126774564385414124,   0.150365769863128662,   0.107535004615783691,   0.152171716094017029,   0.016197171062231064,  -0.097000472247600555,  -0.103922717273235321,   0.024408107623457909,   0.049295961856842041,  -0.043224696069955826,   0.100952073931694031,   0.166513219475746155,   0.122361607849597931,   0.078320965170860291,   0.076739378273487091,   0.007597636431455612,  -0.047421026974916458,   0.047722667455673218,   0.011763826943933964,  -0.023527545854449272,  -0.015104408375918865,  -0.103659719228744507,   0.123738199472427368,   0.014085118658840656,  -0.188047692179679871,   0.020824914798140526,   0.127643674612045288,  -0.068017229437828064,   0.046567216515541077,   0.127315983176231384,   0.042887203395366669,
		  -0.058730911463499069,  -0.063175648450851440,  -0.064359404146671295,   0.176601842045783997,   0.161614105105400085,  -0.188426539301872253,   0.081756740808486938,   0.065574631094932556,  -0.055136844515800476,  -0.010116691701114178,  -0.083298258483409882,  -0.004135891329497099,   0.067304633557796478,   0.163341656327247620,   0.140150293707847595,  -0.145397856831550598,   0.071851745247840881,   0.175411090254783630,   0.084572531282901764,  -0.136685803532600403,  -0.078311033546924591,   0.043670766055583954,  -0.035963442176580429,   0.057761065661907196,   0.041182208806276321,   0.021801330149173737,  -0.137013241648674011,  -0.006426113657653332,  -0.029162591323256493,   0.135249108076095581,   0.097662329673767090,   0.052388258278369904,   0.150383532047271729,   0.052057463675737381,   0.011992097832262516,  -0.078590668737888336,  -0.094547048211097717,   0.018252030014991760,  -0.034115411341190338,   0.041753325611352921,   0.074397929012775421,  -0.174409985542297363,  -0.097301758825778961,   0.088398940861225128,   0.137654304504394531,  -0.031297788023948669,   0.130821421742439270,   0.047466389834880829,
		   0.036888923496007919,  -0.031679138541221619,  -0.133823677897453308,   0.115958303213119507,  -0.053239457309246063,   0.070893928408622742,   0.117271877825260162,  -0.082013063132762909,  -0.161555305123329163,   0.143213257193565369,  -0.149572625756263733,  -0.016493033617734909,   0.093147069215774536,  -0.005353782325983047,  -0.024942032992839813,   0.013650320470333099,  -0.085899956524372101,   0.010471648536622524,   0.076273344457149506,  -0.126730129122734070,   0.128119394183158875,   0.069742619991302490,   0.203203096985816956,   0.155895367264747620,  -0.107169605791568756,   0.151639878749847412,   0.080237679183483124,   0.014461833983659744,  -0.090417809784412384,   0.067607827484607697,   0.021243337541818619,   0.100615136325359344,  -0.049539700150489807,   0.053535867482423782,   0.080389052629470825,  -0.022990835830569267,  -0.015180614776909351,   0.044081866741180420,  -0.043898567557334900,   0.103535279631614685,  -0.020035311579704285,   0.017843840643763542,  -0.047713316977024078,  -0.112088508903980255,  -0.019990317523479462,   0.133214145898818970,   0.104376472532749176,   0.179507091641426086,
		   0.023589562624692917,   0.104902073740959167,  -0.176978632807731628,   0.046594925224781036,   0.059358321130275726,   0.135704934597015381,   0.186822831630706787,  -0.027208806946873665,  -0.173215508460998535,  -0.029839072376489639,  -0.045082952827215195,   0.083355300128459930,  -0.018231112509965897,   0.141959592700004578,   0.009124260395765305,  -0.022832868620753288,  -0.130701616406440735,  -0.039084769785404205,  -0.109201639890670776,   0.062250267714262009,  -0.091671027243137360,   0.010616517625749111,  -0.004173902329057455,   0.074361637234687805,   0.050496634095907211,  -0.057835135608911514,   0.059665318578481674,  -0.084298826754093170,   0.099372975528240204,   0.078880570828914642,  -0.075812600553035736,  -0.038930486887693405,  -0.095059797167778015,  -0.058172281831502914,   0.092057496309280396,  -0.072430275380611420,  -0.096455223858356476,   0.052879132330417633,   0.173765331506729126,  -0.144485518336296082,  -0.131784498691558838,   0.126805439591407776,   0.130588144063949585,  -0.032034020870923996,  -0.059218194335699081,  -0.060626357793807983,   0.088836975395679474,   0.039976276457309723,
		   0.000514866318553686,   0.032273951917886734,  -0.008896890096366405,  -0.010476980358362198,  -0.003552405862137675,  -0.145863071084022522,  -0.168210744857788086,  -0.158436864614486694,   0.047057338058948517,  -0.111749894917011261,  -0.120823644101619720,  -0.024237217381596565,   0.083584263920783997,   0.048812828958034515,  -0.058714732527732849,   0.127528771758079529,  -0.127088099718093872,   0.106374450027942657,  -0.033188708126544952,   0.055181913077831268,   0.082326963543891907,  -0.080125026404857635,  -0.003842835314571857,  -0.172716811299324036,   0.042726609855890274,   0.095685735344886780,  -0.150812000036239624,  -0.100189305841922760,  -0.073529317975044250,  -0.139020740985870361,  -0.096088767051696777,  -0.050497502088546753,   0.117579728364944458,  -0.066283263266086578,  -0.096454679965972900,  -0.052085988223552704,  -0.121990054845809937,   0.063108481466770172,  -0.037297233939170837,   0.086854480206966400,  -0.073316924273967743,   0.024924224242568016,  -0.182348355650901794,   0.111724995076656342,  -0.011157966218888760,   0.015886669978499413,  -0.060297548770904541,  -0.060715455561876297,
		   0.030192375183105469,   0.058155592530965805,  -0.004185057710856199,  -0.100413583219051361,  -0.048437979072332382,   0.024975610896945000,   0.145790815353393555,  -0.036298751831054688,  -0.009314691647887230,   0.004332957323640585,  -0.108257971704006195,   0.055569630116224289,  -0.044102244079113007,   0.163054183125495911,   0.074553154408931732,  -0.106016293168067932,   0.072672009468078613,  -0.043863322585821152,   0.082968078553676605,  -0.175763279199600220,   0.130330935120582581,  -0.083650358021259308,   0.154955461621284485,   0.178152069449424744,   0.014690114185214043,   0.127830252051353455,   0.120054982602596283,   0.090109877288341522,  -0.033966470509767532,   0.198763817548751831,   0.005720805376768112,   0.070901297032833099,   0.157839402556419373,  -0.067669011652469635,   0.102847710251808167,   0.123821407556533813,  -0.027020309120416641,   0.089066103100776672,  -0.042564447969198227,  -0.108920514583587646,   0.151871338486671448,   0.184426382184028625,   0.068297639489173889,  -0.005159481428563595,   0.163695991039276123,  -0.043149948120117188,   0.147074759006500244,   0.166753992438316345,
		  -0.030530877411365509,   0.066378556191921234,   0.171318948268890381,  -0.085301451385021210,  -0.088491618633270264,  -0.051844432950019836,  -0.078672558069229126,  -0.035346399992704391,   0.137038648128509521,   0.065104112029075623,  -0.085053406655788422,  -0.081598550081253052,   0.059595968574285507,  -0.087771721184253693,   0.122590765357017517,  -0.118273712694644928,   0.080369763076305389,   0.148322463035583496,  -0.022896308451890945,   0.082914888858795166,   0.014147988520562649,   0.113934807479381561,   0.101456835865974426,  -0.135195612907409668,   0.149791449308395386,   0.171043023467063904,   0.067719988524913788,   0.086187921464443207,   0.105035521090030670,   0.073754683136940002,  -0.055680524557828903,  -0.002303312532603741,   0.063130326569080353,  -0.072160355746746063,   0.164129614830017090,   0.105393201112747192,  -0.052271813154220581,  -0.056216645985841751,   0.018966676667332649,   0.196270078420639038,   0.016186706721782684,  -0.074244014918804169,  -0.124747477471828461,   0.080068528652191162,  -0.045438200235366821,  -0.016674725338816643,  -0.029256716370582581,  -0.084523178637027740,
		  -0.004239727742969990,  -0.034468587487936020,   0.056333106011152267,  -0.068847015500068665,  -0.013274086639285088,   0.081705071032047272,   0.004972043447196484,  -0.064042121171951294,   0.089812666177749634,  -0.001821957877837121,   0.026443937793374062,   0.048632349818944931,  -0.037891209125518799,  -0.109305523335933685,  -0.063415236771106720,  -0.129854217171669006,   0.103414468467235565,   0.096573241055011749,  -0.120062477886676788,   0.087144069373607635,  -0.020763488486409187,  -0.028950816020369530,   0.147428199648857117,   0.117490373551845551,  -0.023095076903700829,   0.167187005281448364,   0.160071402788162231,  -0.112083256244659424,  -0.105079844594001770,   0.159361630678176880,   0.006338099483400583,   0.061822488903999329,  -0.082816183567047119,  -0.094483613967895508,   0.110726341605186462,  -0.004483164288103580,   0.112130090594291687,  -0.133397445082664490,   0.140395596623420715,  -0.131256237626075745,   0.073886036872863770,   0.070608176290988922,  -0.064337700605392456,  -0.031042950227856636,  -0.060956500470638275,  -0.122203357517719269,  -0.058767221868038177,  -0.082605496048927307,
		  -0.163092166185379028,   0.029787920415401459,  -0.046258013695478439,  -0.113963574171066284,  -0.046139918267726898,  -0.091147087514400482,  -0.015317440964281559,   0.040493644773960114,   0.158480420708656311,   0.077937319874763489,   0.035130336880683899,   0.108311802148818970,   0.143344521522521973,   0.012968736700713634,  -0.028501868247985840,   0.103340424597263336,   0.016270363703370094,  -0.065653540194034576,   0.109813891351222992,  -0.094260372221469879,  -0.048784140497446060,   0.093767173588275909,  -0.047197733074426651,  -0.031502958387136459,   0.055516850203275681,   0.054127536714076996,  -0.013247902505099773,  -0.039201594889163971,  -0.052246533334255219,  -0.049215290695428848,  -0.050822481513023376,  -0.106918886303901672,   0.012062275782227516,   0.063478380441665649,   0.053940508514642715,   0.046720951795578003,  -0.043217871338129044,  -0.029759714379906654,   0.002072174334898591,   0.161631986498832703,   0.071145795285701752,  -0.057390239089727402,  -0.115119040012359619,   0.106034636497497559,   0.046211004257202148,  -0.095334902405738831,   0.064393244683742523,  -0.042558223009109497,
		  -0.161914095282554626,  -0.156817644834518433,   0.025590077042579651,  -0.037054006010293961,   0.058672849088907242,   0.006984765175729990,   0.085058026015758514,  -0.015943745151162148,   0.148453161120414734,   0.004281400702893734,  -0.037655625492334366,   0.023657126352190971,  -0.011025970801711082,   0.052386041730642319,  -0.071958594024181366,  -0.038461849093437195,  -0.057385720312595367,   0.094174370169639587,  -0.010713899508118629,  -0.139230489730834961,  -0.162797972559928894,   0.026291284710168839,  -0.186592698097229004,   0.076299607753753662,   0.017276210710406303,  -0.008933557197451591,   0.043600998818874359,  -0.060759671032428741,   0.062895596027374268,  -0.059975408017635345,  -0.034714102745056152,  -0.119432888925075531,   0.023209702223539352,  -0.058715958148241043,  -0.137481287121772766,  -0.058343235403299332,   0.071581736207008362,   0.131860643625259399,  -0.100807078182697296,  -0.042493041604757309,  -0.156633406877517700,  -0.011083472520112991,  -0.022414891049265862,  -0.071550361812114716,  -0.081109844148159027,   0.006509505212306976,  -0.108233585953712463,   0.123346753418445587,
		  -0.138591021299362183,   0.092247590422630310,  -0.093087077140808105,   0.144676238298416138,   0.009457685053348541,   0.106564633548259735,   0.013701031915843487,  -0.076945170760154724,   0.086299985647201538,   0.060363005846738815,   0.014678422361612320,  -0.066817671060562134,  -0.090086907148361206,   0.087190404534339905,  -0.045289870351552963,   0.045426644384860992,   0.005161944311112165,   0.063697084784507751,   0.013244400732219219,  -0.013722164556384087,   0.066273733973503113,  -0.065346956253051758,  -0.071086354553699493,  -0.037504218518733978,   0.022640146315097809,  -0.074542455375194550,   0.042206667363643646,   0.039582449942827225,  -0.034490082412958145,  -0.017266105860471725,   0.163132429122924805,  -0.051868379116058350,   0.122886948287487030,  -0.048948407173156738,   0.043010219931602478,   0.156206160783767700,   0.116997063159942627,   0.089439883828163147,  -0.081334002315998077,   0.095641262829303741,   0.023233465850353241,  -0.079089999198913574,  -0.121956072747707367,  -0.030082322657108307,   0.160107836127281189,  -0.129954323172569275,   0.114731721580028534,   0.056039933115243912,
	};
	static const double bias03[]=
	{
		   0.117858000099658966,  -0.021966768428683281,   0.085116982460021973,   0.146176815032958984,   0.073609493672847748,  -0.055600911378860474,  -0.149970039725303650,  -0.164250820875167847,  -0.048568796366453171,  -0.018937187269330025,  -0.087151154875755310,  -0.012655377388000488,   0.116997905075550079,   0.038090087473392487,  -0.025833740830421448,   0.123697400093078613,   0.085039138793945312,   0.108884178102016449,   0.173833385109901428,   0.006189849693328142,  -0.033453099429607391,   0.170433133840560913,   0.002527193399146199,   0.106386780738830566,  -0.147308751940727234,  -0.138053715229034424,  -0.095164686441421509,  -0.067167013883590698,   0.112427107989788055,  -0.140160486102104187,  -0.072729974985122681,   0.018991494551301003,   0.048257879912853241,   0.017959244549274445,   0.140728086233139038,  -0.045987978577613831,  -0.050362676382064819,  -0.059924654662609100,  -0.006610093638300896,   0.102964013814926147,   0.050578687340021133,  -0.142631784081459045,  -0.028627648949623108,  -0.031932905316352844,   0.114220455288887024,  -0.007224338129162788,  -0.174019291996955872,  -0.016361726447939873,
	};
	static const double weight04[]=
	{
		  -0.115897327661514282,  -0.076449178159236908,   0.122952193021774292,   0.035970181226730347,   0.119324326515197754,   0.029750876128673553,   0.121927477419376373,  -0.049250431358814240,  -0.019577169790863991,   0.055801078677177429,   0.072643205523490906,   0.127298057079315186,  -0.041810695081949234,   0.158112704753875732,  -0.150354593992233276,   0.024383315816521645,   0.152396917343139648,  -0.062766149640083313,   0.124248221516609192,   0.092508427798748016,   0.024174779653549194,  -0.147844552993774414,  -0.001023009303025901,   0.127308934926986694,   0.038233123719692230,   0.084913864731788635,  -0.091292440891265869,  -0.156813278794288635,   0.088310658931732178,  -0.153887167572975159,   0.014635102823376656,   0.179253220558166504,   0.051027569919824600,  -0.052884750068187714,   0.006569168064743280,   0.078556664288043976,  -0.032395269721746445,   0.045484654605388641,   0.037591490894556046,   0.092428848147392273,  -0.141419455409049988,   0.056091248989105225,   0.119990423321723938,   0.111466370522975922,  -0.012550047598779202,   0.119087174534797668,  -0.004579711239784956,   0.120215140283107758,
		  -0.033248003572225571,   0.192599892616271973,  -0.097318559885025024,  -0.005625420715659857,   0.033782143145799637,  -0.184931352734565735,   0.120203316211700439,   0.022799748927354813,   0.109015643596649170,   0.099286675453186035,  -0.010435753501951694,   0.168903559446334839,  -0.127513349056243896,   0.154504209756851196,  -0.014021827839314938,   0.186079844832420349,   0.137739211320877075,  -0.016993790864944458,   0.052839647978544235,   0.079946465790271759,   0.143316969275474548,   0.063031367957592010,  -0.032633341848850250,   0.165103256702423096,  -0.054597247391939163,   0.019580595195293427,  -0.025633834302425385,   0.035938635468482971,   0.040787801146507263,  -0.069199204444885254,  -0.060616482049226761,  -0.030188577249646187,  -0.021457156166434288,  -0.059257384389638901,   0.134516030550003052,   0.112668074667453766,  -0.091659530997276306,   0.088345810770988464,   0.049837186932563782,  -0.028540318831801414,  -0.132265895605087280,  -0.057434849441051483,  -0.033419147133827209,  -0.028261179104447365,   0.021897127851843834,   0.079368658363819122,  -0.121878668665885925,   0.152517706155776978,
		   0.060233082622289658,   0.169061437249183655,  -0.041324462741613388,   0.041140716522932053,   0.064192652702331543,  -0.115057311952114105,  -0.130237713456153870,  -0.055170737206935883,  -0.021550206467509270,   0.127763450145721436,  -0.170284241437911987,   0.024945957586169243,   0.001507605193182826,  -0.084106959402561188,  -0.139286279678344727,   0.026255013421177864,   0.146119400858879089,  -0.153279721736907959,   0.189586535096168518,  -0.078564718365669250,   0.149271354079246521,  -0.073809392750263214,  -0.088977210223674774,   0.161885425448417664,  -0.010130674578249454,  -0.088567145168781281,   0.046493649482727051,  -0.046989545226097107,   0.139106839895248413,  -0.045776218175888062,   0.103995837271213531,   0.080370567739009857,   0.098997965455055237,   0.008442334830760956,   0.176790252327919006,   0.114124163985252380,  -0.006751014385372400,   0.144256904721260071,   0.183764308691024780,  -0.081390514969825745,  -0.160852074623107910,  -0.062707163393497467,   0.062275223433971405,  -0.089716508984565735,   0.068421289324760437,   0.025765938684344292,  -0.078093841671943665,   0.037308827042579651,
		  -0.102223224937915802,   0.051609102636575699,   0.069766156375408173,  -0.054562114179134369,   0.085654221475124359,   0.123642764985561371,   0.103597149252891541,  -0.011279548518359661,  -0.110084109008312225,  -0.033774543553590775,  -0.063396893441677094,  -0.027685044333338737,   0.122054785490036011,   0.106193318963050842,   0.154147520661354065,  -0.002523908391594887,  -0.070101767778396606,   0.151810139417648315,  -0.125202417373657227,  -0.083209209144115448,   0.068457968533039093,   0.157462656497955322,  -0.048034451901912689,   0.010555457323789597,  -0.165188610553741455,  -0.050362356007099152,   0.020381627604365349,   0.061842381954193115,   0.074516281485557556,   0.081853784620761871,   0.047002136707305908,   0.134735673666000366,   0.093097880482673645,  -0.088076665997505188,   0.065891563892364502,   0.098089531064033508,  -0.103725723922252655,   0.140494465827941895,   0.161886066198348999,  -0.036317728459835052,   0.140241906046867371,  -0.159271374344825745,   0.159787327051162720,   0.155082643032073975,  -0.076943233609199524,   0.051458604633808136,  -0.107256531715393066,  -0.122601576149463654,
		   0.071862719953060150,  -0.076918691396713257,   0.004135248251259327,   0.037340022623538971,   0.020613787695765495,  -0.043200612068176270,   0.092524901032447815,  -0.120759978890419006,   0.035642255097627640,   0.056404568254947662,  -0.149012967944145203,   0.104235783219337463,  -0.007984924130141735,  -0.088361278176307678,   0.156812623143196106,   0.028284216299653053,   0.045157544314861298,  -0.034893743693828583,   0.119053423404693604,   0.097449049353599548,   0.083742402493953705,   0.200999245047569275,   0.165614083409309387,   0.121337331831455231,   0.089079141616821289,   0.020273601636290550,  -0.144820824265480042,   0.165004223585128784,  -0.075624629855155945,   0.060912173241376877,   0.086221583187580109,  -0.115109398961067200,   0.107032984495162964,  -0.040678866207599640,  -0.088175214827060699,   0.056028150022029877,  -0.080130666494369507,  -0.145799309015274048,  -0.159735053777694702,   0.075603522360324860,   0.009060886688530445,   0.053962159901857376,   0.070135518908500671,  -0.058049309998750687,   0.132520347833633423,   0.072561539709568024,  -0.073774076998233795,   0.080541923642158508,
		   0.026023687794804573,  -0.112075328826904297,   0.097196325659751892,   0.149884194135665894,  -0.003367675701156259,   0.086583986878395081,   0.041722860187292099,  -0.047359835356473923,  -0.040502808988094330,   0.054597120732069016,  -0.062303178012371063,  -0.004814057610929012,   0.173128798604011536,   0.120132386684417725,   0.171613201498985291,  -0.090044781565666199,  -0.014292991720139980,   0.153984919190406799,   0.121537543833255768,   0.003212999319657683,   0.081470079720020294,   0.173924222588539124,   0.133544877171516418,  -0.108672074973583221,   0.043304517865180969,  -0.167885795235633850,  -0.106385007500648499,   0.173386678099632263,   0.074950762093067169,  -0.121090635657310486,   0.047233521938323975,   0.039421297609806061,  -0.123547807335853577,  -0.039510987699031830,   0.002700397279113531,   0.069836579263210297,   0.112911581993103027,  -0.122108131647109985,  -0.019350919872522354,   0.073115803301334381,  -0.039155852049589157,   0.089195616543292999,  -0.021262686699628830,   0.088122472167015076,  -0.080067589879035950,   0.032546170055866241,  -0.064399018883705139,  -0.036741442978382111,
		  -0.112739771604537964,  -0.051669601351022720,  -0.034433219581842422,   0.015325043350458145,   0.067443020641803741,   0.088457010686397552,   0.133549541234970093,  -0.137348979711532593,  -0.127181157469749451,   0.086291819810867310,  -0.087680220603942871,  -0.088869199156761169,  -0.021688811480998993,  -0.072967015206813812,   0.019988780841231346,  -0.110309913754463196,  -0.058334223926067352,   0.039190884679555893,  -0.107061974704265594,   0.133871018886566162,   0.031396158039569855,  -0.046346362680196762,  -0.127073824405670166,  -0.074240364134311676,   0.122961245477199554,   0.036800086498260498,   0.127206176519393921,   0.110101483762264252,   0.072848662734031677,   0.089791126549243927,   0.095986120402812958,  -0.076691471040248871,  -0.076052606105804443,   0.135962232947349548,  -0.108978100121021271,  -0.093199901282787323,  -0.061092261224985123,   0.031221436336636543,   0.131446391344070435,   0.120282799005508423,  -0.008702601306140423,   0.074821338057518005,  -0.115905694663524628,   0.015427982434630394,  -0.011578019708395004,  -0.020751899108290672,   0.002883298555389047,  -0.130428001284599304,
		   0.091370567679405212,   0.058920428156852722,  -0.001470951712690294,  -0.073113307356834412,  -0.020520323887467384,   0.120597392320632935,   0.093128621578216553,   0.012916733510792255,   0.025029042735695839,   0.017635088413953781,  -0.080726638436317444,  -0.135861769318580627,  -0.050946079194545746,  -0.137261450290679932,  -0.002379801357164979,  -0.144935160875320435,   0.144617825746536255,   0.107428491115570068,   0.049791179597377777,   0.106798790395259857,  -0.007489241659641266,   0.136527389287948608,  -0.055594239383935928,   0.047058727592229843,  -0.133896619081497192,   0.039164345711469650,  -0.152971744537353516,   0.000962738413363695,  -0.132975593209266663,   0.114698760211467743,   0.133526608347892761,   0.038842484354972839,  -0.044739238917827606,  -0.056265827268362045,   0.112119071185588837,   0.043007235974073410,   0.075506344437599182,   0.090179115533828735,   0.087660327553749084,   0.027968004345893860,   0.039030749350786209,   0.077298283576965332,   0.159257814288139343,   0.061025653034448624,   0.081539906561374664,   0.041652392596006393,   0.120657771825790405,  -0.116095408797264099,
		  -0.015476961620151997,   0.032809138298034668,   0.092701025307178497,  -0.048301231116056442,   0.098115488886833191,   0.015865841880440712,   0.013763030059635639,  -0.079558819532394409,  -0.115589007735252380,   0.168113127350807190,  -0.087145954370498657,  -0.091844178736209869,  -0.028987273573875427,   0.115142516791820526,   0.051323007792234421,   0.051957957446575165,   0.044442456215620041,  -0.097471542656421661,   0.042114533483982086,  -0.032796408981084824,   0.004458025097846985,  -0.125774577260017395,  -0.088185057044029236,   0.186186775565147400,  -0.138335034251213074,  -0.084442652761936188,   0.081798195838928223,   0.049920186400413513,   0.089816741645336151,  -0.146981939673423767,   0.027353253215551376,   0.155803456902503967,   0.149531677365303040,  -0.048802535980939865,  -0.058201093226671219,  -0.024658778682351112,  -0.084419578313827515,   0.112490423023700714,   0.127202212810516357,   0.093020237982273102,  -0.017948979511857033,  -0.011257722973823547,   0.153035402297973633,  -0.026588905602693558,  -0.120627097785472870,   0.130315348505973816,   0.095285214483737946,   0.098333261907100677,
		   0.093699693679809570,  -0.057061672210693359,  -0.113118439912796021,   0.128873676061630249,   0.038408536463975906,  -0.102259077131748199,  -0.013989473693072796,  -0.005840197205543518,  -0.010771326720714569,  -0.031679801642894745,  -0.027610210701823235,  -0.126728937029838562,   0.016768500208854675,   0.074309363961219788,   0.135560199618339539,  -0.059518858790397644,  -0.110513851046562195,  -0.007753319572657347,   0.062257129698991776,   0.125410109758377075,  -0.057935535907745361,  -0.022274963557720184,  -0.039647337049245834,  -0.063285842537879944,  -0.121572501957416534,   0.018123516812920570,   0.069806113839149475,   0.121805839240550995,  -0.115655153989791870,  -0.059167094528675079,  -0.037599656730890274,  -0.012399018742144108,  -0.130150958895683289,   0.034819152206182480,  -0.087220154702663422,  -0.028454385697841644,   0.018846411257982254,  -0.101841539144515991,  -0.021380431950092316,  -0.029521796852350235,  -0.066227912902832031,   0.117695823311805725,   0.008024424314498901,   0.056866105645895004,  -0.015772502869367599,  -0.128843098878860474,  -0.049595661461353302,  -0.127649262547492981,
		   0.018125640228390694,  -0.080300301313400269,   0.016118105500936508,  -0.031796723604202271,  -0.059313561767339706,   0.059424530714750290,   0.080695450305938721,   0.053168833255767822,   0.034908551722764969,   0.131101235747337341,  -0.112114883959293365,  -0.060677137225866318,   0.040499538183212280,  -0.127459898591041565,   0.140848010778427124,   0.047972280532121658,   0.012810993939638138,   0.002596860285848379,   0.103750512003898621,   0.037112671881914139,   0.083554409444332123,  -0.032543092966079712,  -0.000574679579585791,  -0.020224999636411667,   0.135528132319450378,   0.071270175278186798,  -0.117028333246707916,   0.131382033228874207,   0.067789621651172638,   0.018136434257030487,   0.081012479960918427,   0.110534116625785828,   0.011084916070103645,   0.021486222743988037,  -0.007736470550298691,   0.080871745944023132,   0.104773305356502533,   0.120264820754528046,  -0.112661845982074738,   0.092197038233280182,   0.016049211844801903,   0.119199790060520172,   0.111406981945037842,   0.048446882516145706,  -0.010578973218798637,   0.107825532555580139,  -0.017620369791984558,  -0.133428499102592468,
		   0.047308202832937241,  -0.089198626577854156,   0.085245169699192047,   0.027961878105998039,   0.067429155111312866,   0.036174155771732330,   0.112255707383155823,   0.032088689506053925,   0.058839399367570877,  -0.048670668154954910,   0.117125116288661957,   0.066468179225921631,   0.095124021172523499,   0.128440588712692261,  -0.041930656880140305,  -0.036435287445783615,   0.011343053542077541,  -0.076917968690395355,   0.064583674073219299,  -0.163687899708747864,   0.125051960349082947,   0.004974269773811102,   0.019204333424568176,   0.030048420652747154,  -0.111094586551189423,   0.153606534004211426,   0.076158858835697174,  -0.020967429503798485,   0.034629806876182556,   0.001292623812332749,  -0.001076119835488498,   0.108971774578094482,  -0.107895448803901672,  -0.142471298575401306,  -0.026442892849445343,   0.035145577043294907,   0.050733000040054321,   0.118622042238712311,   0.062679558992385864,  -0.123305588960647583,  -0.171773478388786316,   0.154865130782127380,  -0.068421825766563416,  -0.002517958870157599,   0.018263787031173706,  -0.075624398887157440,   0.094499059021472931,   0.152292922139167786,
		  -0.102922290563583374,  -0.153248369693756104,  -0.013333288952708244,   0.141065791249275208,   0.039021313190460205,   0.151186957955360413,  -0.049371819943189621,  -0.060123335570096970,  -0.039800498634576797,  -0.115741863846778870,  -0.152427181601524353,   0.008214961737394333,   0.091649882495403290,  -0.118952058255672455,   0.011117864400148392,   0.110083296895027161,   0.125910833477973938,   0.019605187699198723,  -0.075520724058151245,  -0.129156038165092468,   0.160144731402397156,   0.170434370636940002,   0.158920541405677795,  -0.133989825844764709,  -0.120735213160514832,   0.016589555889368057,   0.085311584174633026,  -0.076220400631427765,   0.066231660544872284,   0.108768902719020844,   0.095161996781826019,  -0.030406160280108452,  -0.116862639784812927,   0.175305873155593872,  -0.004890990443527699,  -0.026736903935670853,  -0.107158288359642029,  -0.002197685884311795,  -0.108445622026920319,   0.073452800512313843,  -0.081315927207469940,   0.011127237230539322,   0.060863930732011795,  -0.164265245199203491,  -0.053349107503890991,  -0.109788425266742706,  -0.074398554861545563,  -0.146414950489997864,
		  -0.133264780044555664,   0.084658883512020111,   0.051698666065931320,  -0.055357597768306732,  -0.073066137731075287,   0.117560960352420807,   0.008305262774229050,   0.120161853730678558,  -0.015988534316420555,  -0.047910831868648529,   0.062843292951583862,   0.115926705300807953,   0.093627788126468658,  -0.082817725837230682,   0.143604934215545654,  -0.138411268591880798,  -0.024279642850160599,   0.117436379194259644,   0.016839943826198578,   0.086246132850646973,  -0.044072881340980530,  -0.051346484571695328,  -0.097513906657695770,  -0.029582735151052475,  -0.098913542926311493,  -0.129786819219589233,  -0.002031886950135231,  -0.116188667714595795,   0.077121764421463013,  -0.086411036550998688,  -0.088030144572257996,  -0.085891783237457275,  -0.076701872050762177,   0.030568517744541168,  -0.055061165243387222,   0.007508701644837856,   0.082964986562728882,  -0.124927967786788940,  -0.028454920276999474,  -0.098895564675331116,   0.002583089983090758,   0.048367194831371307,  -0.100972883403301239,   0.123190045356750488,  -0.007833020761609077,  -0.111970543861389160,  -0.139850810170173645,  -0.094563685357570648,
		  -0.121171876788139343,  -0.037189833819866180,   0.026651177555322647,  -0.034802537411451340,   0.182901635766029358,  -0.056791074573993683,  -0.016004092991352081,  -0.160881981253623962,   0.023967061191797256,   0.160633429884910583,  -0.101085036993026733,   0.022122146561741829,  -0.064122438430786133,  -0.087631255388259888,   0.040945537388324738,   0.184918940067291260,   0.087157644331455231,  -0.005168618168681860,   0.128478765487670898,  -0.017008503898978233,   0.110423400998115540,  -0.038938250392675400,  -0.013936867006123066,   0.168925076723098755,   0.059854958206415176,  -0.058586750179529190,   0.021891808137297630,   0.058584291487932205,   0.023249717429280281,  -0.161650523543357849,   0.068173445761203766,   0.133008867502212524,  -0.038309648633003235,  -0.030411921441555023,   0.069175213575363159,  -0.037427097558975220,   0.149240180850028992,   0.133535131812095642,   0.185695052146911621,  -0.094009011983871460,  -0.104562081396579742,  -0.023894645273685455,   0.108399204909801483,   0.147720113396644592,  -0.115923486649990082,   0.068414248526096344,  -0.046470444649457932,  -0.047421108931303024,
		   0.094527810811996460,   0.083842761814594269,   0.174163758754730225,   0.084734879434108734,   0.010203368961811066,  -0.050686739385128021,   0.090158373117446899,  -0.062741763889789581,   0.083122521638870239,   0.180457204580307007,  -0.071135915815830231,   0.077816948294639587,   0.012098865583539009,   0.187763050198554993,  -0.093744248151779175,   0.163758158683776855,  -0.044772908091545105,  -0.133211866021156311,  -0.073556631803512573,  -0.023010434582829475,   0.063654735684394836,   0.061681244522333145,   0.069346904754638672,  -0.068575397133827209,  -0.159654855728149414,  -0.135025694966316223,   0.042228043079376221,   0.043657910078763962,  -0.025761133059859276,  -0.128313496708869934,   0.012926561757922173,  -0.022360596805810928,   0.048938121646642685,  -0.190365269780158997,  -0.036372147500514984,  -0.011213173158466816,   0.127451360225677490,   0.172095179557800293,   0.026388363912701607,   0.149886742234230042,  -0.061489660292863846,   0.060375038534402847,   0.148522958159446716,   0.118460282683372498,  -0.112052567303180695,  -0.010875254869461060,  -0.105675965547561646,  -0.097991950809955597,
		  -0.146886572241783142,   0.076851710677146912,  -0.079672567546367645,  -0.062327042222023010,  -0.075969092547893524,   0.038027595728635788,  -0.175813764333724976,   0.002058747690171003,  -0.138417929410934448,  -0.015636209398508072,  -0.127603799104690552,  -0.087899506092071533,   0.190037548542022705,   0.012697163037955761,   0.179395407438278198,  -0.096075579524040222,   0.184138134121894836,   0.019881658256053925,   0.053025647997856140,  -0.111041493713855743,  -0.055915128439664841,   0.181853070855140686,   0.175998181104660034,  -0.016691766679286957,  -0.062514021992683411,  -0.193330571055412292,   0.002988288644701242,   0.028932621702551842,   0.009053614921867847,   0.095798552036285400,  -0.057074490934610367,  -0.066437475383281708,   0.086289480328559875,  -0.027194209396839142,  -0.081443890929222107,   0.034690953791141510,  -0.078463509678840637,   0.042897652834653854,  -0.032827988266944885,  -0.027051638811826706,   0.054470662027597427,  -0.172119349241256714,   0.071526341140270233,   0.113435722887516022,   0.216300383210182190,  -0.014672688208520412,   0.058274935930967331,   0.062977708876132965,
		   0.048404973000288010,   0.116761840879917145,  -0.107982791960239410,  -0.042928539216518402,   0.134401425719261169,  -0.124060161411762238,   0.000495770073030144,  -0.097223795950412750,   0.020714385434985161,   0.068185433745384216,  -0.045696124434471130,  -0.011983656324446201,   0.081900343298912048,  -0.024610966444015503,   0.018450716510415077,   0.128373116254806519,  -0.015591611154377460,  -0.147288918495178223,   0.025285264477133751,  -0.097695194184780121,   0.079007603228092194,  -0.031172247603535652,   0.076756700873374939,   0.192802295088768005,   0.070901937782764435,   0.108744442462921143,   0.082541070878505707,  -0.051940623670816422,   0.025965970009565353,  -0.096373684704303741,  -0.023593220859766006,   0.020442454144358635,   0.147765889763832092,  -0.082595728337764740,  -0.032679613679647446,   0.109687343239784241,   0.171980321407318115,   0.055198974907398224,   0.138614580035209656,   0.022985044866800308,   0.002000879030674696,   0.065651901066303253,   0.038297552615404129,   0.117829173803329468,  -0.079904295504093170,   0.103200078010559082,   0.025687651708722115,   0.067368760704994202,
		   0.116383001208305359,  -0.009628774598240852,  -0.064001835882663727,  -0.067788168787956238,   0.162780448794364929,  -0.028821354731917381,  -0.022460399195551872,  -0.004059081431478262,  -0.031190674751996994,   0.065474115312099457,  -0.066727407276630402,   0.036615259945392609,  -0.028091602027416229,   0.167003065347671509,  -0.086922079324722290,   0.039739415049552917,  -0.101258777081966400,  -0.105993092060089111,   0.101974934339523315,  -0.041548140347003937,   0.173739016056060791,  -0.144154414534568787,   0.103788748383522034,  -0.075084842741489410,  -0.019968381151556969,  -0.081017136573791504,  -0.051686719059944153,   0.055023659020662308,   0.016761587932705879,   0.025659549981355667,  -0.001206888118758798,   0.172578319907188416,  -0.036610491573810577,  -0.075571849942207336,   0.120019085705280304,   0.068482510745525360,   0.174191281199455261,   0.191031858325004578,   0.017508579418063164,   0.006156698334962130,   0.054763529449701309,  -0.129070475697517395,  -0.065459117293357849,   0.068553574383258820,   0.035263214260339737,   0.030982740223407745,  -0.164433702826499939,   0.062509752810001373,
		   0.084143564105033875,  -0.030807886272668839,  -0.027564946562051773,   0.047794938087463379,  -0.102448761463165283,  -0.015145293436944485,  -0.004922516178339720,   0.066594325006008148,   0.034675449132919312,   0.075159192085266113,  -0.118840500712394714,   0.044604770839214325,   0.089537039399147034,   0.060202419757843018,   0.144439458847045898,  -0.072922475636005402,   0.108523413538932800,   0.148951143026351929,   0.155415177345275879,   0.111942708492279053,   0.072235502302646637,  -0.010284500196576118,   0.091501556336879730,   0.063910648226737976,  -0.056471802294254303,   0.002964147599413991,  -0.154923453927040100,   0.102029852569103241,  -0.083019733428955078,  -0.172215700149536133,  -0.002462001983076334,  -0.073840133845806122,  -0.031805623322725296,  -0.016113290563225746,  -0.074652954936027527,   0.140409395098686218,  -0.156474485993385315,  -0.107819743454456329,  -0.171900153160095215,   0.053344801068305969,   0.094077453017234802,  -0.074737906455993652,   0.017078779637813568,  -0.108215019106864929,   0.020378971472382545,  -0.097965471446514130,   0.007481216918677092,   0.109157048165798187,
		   0.070853397250175476,   0.140621945261955261,  -0.100699484348297119,  -0.046976201236248016,   0.138093695044517517,  -0.115802533924579620,  -0.006944343913346529,  -0.004734279122203588,  -0.094197243452072144,  -0.028494026511907578,   0.095375813543796539,  -0.081375837326049805,   0.116423197090625763,   0.130601719021797180,  -0.073243677616119385,   0.169381991028785706,  -0.022245209664106369,   0.049901649355888367,  -0.039242532104253769,   0.102051250636577606,  -0.042182430624961853,  -0.147691294550895691,   0.060466509312391281,   0.121936164796352386,  -0.149243533611297607,  -0.096213757991790771,   0.180283129215240479,  -0.137654945254325867,   0.189996778964996338,   0.003093535779044032,  -0.072948910295963287,   0.029097840189933777,   0.013934483751654625,   0.020324636250734329,  -0.072232216596603394,  -0.008852294646203518,   0.037290450185537338,   0.047482941299676895,  -0.023040549829602242,  -0.036524321883916855,   0.116728112101554871,  -0.045967806130647659,  -0.052073869854211807,  -0.085708536207675934,  -0.070097863674163818,   0.028653282672166824,   0.027153087779879570,   0.114366419613361359,
		   0.082481920719146729,   0.096789456903934479,   0.108487218618392944,  -0.124955318868160248,   0.074008703231811523,  -0.033024828881025314,  -0.017247920855879784,  -0.088488250970840454,   0.084472991526126862,   0.044799964874982834,   0.132604837417602539,   0.024566493928432465,  -0.088669605553150177,  -0.115407802164554596,   0.062715426087379456,  -0.061199922114610672,  -0.014828287065029144,   0.127606570720672607,  -0.144846707582473755,  -0.087086200714111328,   0.048626605421304703,  -0.022835912182927132,  -0.022825300693511963,   0.122914627194404602,  -0.098182134330272675,  -0.126919239759445190,   0.108895286917686462,   0.033770754933357239,  -0.053225792944431305,  -0.100371748208999634,   0.044813219457864761,  -0.051240161061286926,   0.085793137550354004,  -0.099657267332077026,  -0.148419320583343506,   0.005516600329428911,  -0.071107931435108185,   0.099565677344799042,   0.082651823759078979,  -0.034444127231836319,  -0.100831568241119385,   0.099417112767696381,   0.101918347179889679,  -0.012522447854280472,   0.024018751457333565,  -0.087557904422283173,   0.053299188613891602,   0.042219538241624832,
		   0.084958009421825409,   0.020708745345473289,   0.079909533262252808,   0.046485457569360733,  -0.110039979219436646,   0.056132230907678604,   0.091900490224361420,  -0.027578057721257210,  -0.125247895717620850,  -0.060887504369020462,   0.058405525982379913,  -0.100989110767841339,   0.023657577112317085,  -0.048503275960683823,   0.029211480170488358,  -0.035389948636293411,   0.053076345473527908,  -0.083369605243206024,  -0.073102787137031555,  -0.022461777552962303,   0.059426974505186081,  -0.021091628819704056,   0.146242290735244751,   0.064096957445144653,  -0.073271572589874268,   0.100242935121059418,  -0.017913272604346275,  -0.040602926164865494,   0.040213786065578461,  -0.102219462394714355,   0.056791614741086960,  -0.102636046707630157,  -0.094991661608219147,   0.121284186840057373,   0.050067279487848282,   0.181431859731674194,  -0.041994120925664902,   0.035689018666744232,  -0.086055777966976166,  -0.001627427758648992,  -0.048642296344041824,   0.013602316379547119,   0.068309031426906586,  -0.035827603191137314,   0.149372383952140808,  -0.031928703188896179,   0.039889272302389145,  -0.003027830272912979,
		   0.041217941790819168,  -0.077539995312690735,  -0.034047540277242661,   0.019842267036437988,  -0.005615831818431616,  -0.118092328310012817,  -0.074331954121589661,  -0.059110786765813828,  -0.003473140299320221,  -0.088402494788169861,  -0.101070530712604523,  -0.058962427079677582,  -0.105411864817142487,  -0.011233205907046795,   0.039096858352422714,   0.173697188496589661,  -0.108483083546161652,  -0.120102666318416595,  -0.130510389804840088,  -0.118647634983062744,   0.076955877244472504,  -0.139560103416442871,   0.082538060843944550,   0.036456454545259476,  -0.031379364430904388,  -0.003740052226930857,  -0.121308706700801849,   0.045174907892942429,   0.122128792107105255,  -0.134284451603889465,  -0.097074151039123535,  -0.112630747258663177,   0.035408370196819305,  -0.099881134927272797,  -0.000775563414208591,   0.132861539721488953,   0.059059705585241318,   0.154054313898086548,   0.058968175202608109,  -0.143284603953361511,  -0.149794176220893860,   0.134086653590202332,  -0.041945718228816986,  -0.022914415225386620,  -0.132523611187934875,  -0.106047876179218292,  -0.101665735244750977,  -0.091767765581607819,
		   0.063055664300918579,   0.093418613076210022,   0.004576072562485933,   0.097365401685237885,   0.024680862203240395,  -0.141984030604362488,   0.103011213243007660,  -0.092507623136043549,  -0.158321082592010498,   0.129203602671623230,  -0.095187775790691376,  -0.083763159811496735,  -0.065109848976135254,  -0.091448299586772919,  -0.102588102221488953,  -0.053839836269617081,  -0.054443839937448502,  -0.125875800848007202,   0.121050894260406494,  -0.033402107656002045,  -0.050805181264877319,  -0.050220992416143417,  -0.139351829886436462,   0.176247894763946533,  -0.022127998992800713,  -0.013404021970927715,   0.040213279426097870,   0.099841192364692688,  -0.083732686936855316,  -0.167195409536361694,  -0.122091852128505707,   0.126223281025886536,   0.022567428648471832,  -0.185243874788284302,   0.093531958758831024,   0.168038815259933472,  -0.075609862804412842,  -0.058350063860416412,   0.176563456654548645,   0.081179112195968628,  -0.139567777514457703,  -0.095020733773708344,   0.116835847496986389,  -0.039696544408798218,   0.087253622710704803,  -0.075678035616874695,   0.089686505496501923,   0.169472068548202515,
		  -0.104012303054332733,   0.132246091961860657,  -0.068616159260272980,   0.105320073664188385,  -0.103146158158779144,  -0.026883946731686592,   0.008185348473489285,  -0.153626143932342529,  -0.112847402691841125,   0.187710851430892944,   0.080510646104812622,   0.060743052512407303,   0.026260098442435265,   0.136366173624992371,   0.070244565606117249,   0.177115917205810547,   0.020746164023876190,   0.044781167060136795,   0.106704100966453552,  -0.130872905254364014,   0.106864839792251587,   0.086567819118499756,   0.058242671191692352,   0.160442501306533813,  -0.003317281138151884,  -0.053060635924339294,   0.161562830209732056,   0.093857750296592712,   0.050920050591230392,   0.073934316635131836,   0.138890907168388367,  -0.097681708633899689,   0.124376520514488220,  -0.128226459026336670,  -0.082453303039073944,   0.029168792068958282,   0.094274446368217468,   0.140444099903106689,   0.092064052820205688,   0.008654256351292133,  -0.137673065066337585,  -0.007518874481320381,   0.092483684420585632,  -0.026628054678440094,   0.104027740657329559,  -0.058162145316600800,  -0.117850236594676971,   0.154303714632987976,
		   0.009927367791533470,   0.081906966865062714,   0.046074043959379196,  -0.063124440610408783,   0.049252979457378387,   0.164742201566696167,   0.028009103611111641,  -0.097235210239887238,   0.136941686272621155,  -0.035694681107997894,  -0.103224687278270721,   0.011305827647447586,  -0.095499679446220398,  -0.046731006354093552,   0.127206563949584961,  -0.116426765918731689,   0.135068967938423157,   0.178187042474746704,   0.104184590280056000,  -0.072495728731155396,  -0.007486274931579828,  -0.018966259434819221,   0.036187797784805298,  -0.131531402468681335,   0.010798939503729343,  -0.134031102061271667,   0.018050255253911018,  -0.085963547229766846,  -0.108843840658664703,  -0.148124754428863525,  -0.117923937737941742,  -0.016764249652624130,  -0.062106166034936905,  -0.049035910516977310,   0.010711697861552238,   0.036065131425857544,  -0.033370360732078552,   0.114113010466098785,   0.004094246309250593,   0.120288334786891937,   0.013395000249147415,   0.047097846865653992,   0.025649221614003181,  -0.066180869936943054,   0.185463517904281616,  -0.148666203022003174,  -0.072798408567905426,  -0.015591585077345371,
		  -0.026682564988732338,  -0.085821136832237244,   0.166764065623283386,   0.193998202681541443,  -0.040587898343801498,   0.181879401206970215,  -0.120780028402805328,   0.075711309909820557,   0.058186639100313187,  -0.022060936316847801,  -0.023700907826423645,   0.002623883076012135,   0.145080789923667908,   0.084951825439929962,   0.037766195833683014,  -0.121060460805892944,   0.092744141817092896,   0.132162541151046753,   0.107935220003128052,  -0.079405665397644043,   0.096850305795669556,  -0.093823760747909546,  -0.085364863276481628,  -0.038019236177206039,  -0.027546426281332970,  -0.163230314850807190,   0.056756291538476944,  -0.027273656800389290,   0.031523909419775009,  -0.129030555486679077,  -0.032807029783725739,   0.100208252668380737,   0.007071627303957939,   0.166822671890258789,   0.060506723821163177,  -0.048061132431030273,   0.038687624037265778,  -0.091577783226966858,   0.013729005120694637,   0.158486902713775635,   0.018781198188662529,  -0.124057009816169739,   0.091370552778244019,   0.029967125505208969,  -0.023437540978193283,  -0.150929197669029236,  -0.170436292886734009,   0.085768990218639374,
		  -0.167216137051582336,  -0.009648022241890430,   0.094926692545413971,   0.061398040503263474,  -0.082208983600139618,  -0.067784413695335388,  -0.142463788390159607,  -0.120807029306888580,  -0.129167228937149048,   0.035449896007776260,  -0.180910110473632812,  -0.044830266386270523,   0.049199655652046204,   0.155916944146156311,   0.016401065513491631,   0.060696508735418320,   0.083322331309318542,  -0.000542042427696288,  -0.081618517637252808,  -0.158883303403854370,  -0.057573277503252029,   0.032116528600454330,   0.094051480293273926,  -0.018183790147304535,   0.057490956038236618,  -0.061004694551229477,  -0.010081727989017963,   0.072174854576587677,   0.068931706249713898,  -0.113640725612640381,   0.000217674998566508,   0.128104999661445618,   0.039620384573936462,  -0.108896911144256592,   0.017433691769838333,  -0.021114934235811234,   0.081426493823528290,   0.075653925538063049,  -0.024056993424892426,   0.163834124803543091,  -0.037711448967456818,   0.036828290671110153,   0.198593392968177795,   0.013729969039559364,  -0.095433481037616730,  -0.150771975517272949,   0.016624623909592628,  -0.104964502155780792,
		  -0.045261386781930923,  -0.127639442682266235,   0.040784645825624466,   0.088864400982856750,  -0.114350296556949615,   0.008333771489560604,  -0.077067434787750244,  -0.108581408858299255,   0.086528226733207703,   0.082039207220077515,  -0.090950258076190948,  -0.010815798304975033,   0.200805410742759705,  -0.059198904782533646,   0.074956268072128296,  -0.041679419577121735,  -0.026708334684371948,   0.008108433336019516,   0.074008859694004059,  -0.063718192279338837,   0.010753368958830833,   0.188597187399864197,   0.169457346200942993,   0.065891154110431671,   0.092801399528980255,  -0.046470455825328827,  -0.109499804675579071,   0.168430089950561523,  -0.004504296462982893,   0.053532056510448456,  -0.127519562840461731,  -0.065038368105888367,  -0.013174024410545826,   0.051154185086488724,  -0.081555359065532684,  -0.081303693354129791,   0.009713215753436089,  -0.090290382504463196,   0.019615609198808670,   0.131062850356101990,   0.108824260532855988,  -0.117202349007129669,   0.118525333702564240,  -0.131621122360229492,   0.007249442860484123,  -0.120211459696292877,  -0.135658711194992065,   0.057682000100612640,
		  -0.107777230441570282,   0.092904426157474518,  -0.058357264846563339,   0.114650234580039978,   0.119529120624065399,   0.143259987235069275,   0.016238383948802948,  -0.037781883031129837,   0.068561509251594543,  -0.129460036754608154,  -0.098748959600925446,  -0.058631584048271179,   0.034413494169712067,  -0.126205593347549438,   0.074394427239894867,   0.089090459048748016,   0.013313215225934982,   0.128783196210861206,  -0.097104422748088837,   0.072802759706974030,   0.144467398524284363,   0.051539875566959381,   0.038014635443687439,  -0.090531736612319946,  -0.152028486132621765,  -0.190623641014099121,  -0.102597080171108246,   0.153319314122200012,   0.062591120600700378,   0.069383628666400909,   0.121818684041500092,  -0.151284247636795044,  -0.056189097464084625,   0.067262984812259674,   0.026870265603065491,  -0.013107406906783581,   0.084771744906902313,   0.036821987479925156,  -0.137809947133064270,   0.124232493340969086,  -0.006465752609074116,  -0.034978024661540985,   0.195123478770256042,  -0.149280026555061340,   0.072704017162322998,  -0.109990186989307404,  -0.165269643068313599,   0.006068970542401075,
		  -0.059516526758670807,  -0.011807438917458057,  -0.128857031464576721,  -0.005131216719746590,   0.008008795790374279,  -0.053679186850786209,   0.027787731960415840,   0.144730523228645325,   0.133722037076950073,   0.083428397774696350,  -0.000125465798191726,  -0.088011518120765686,  -0.008891222998499870,  -0.066899314522743225,   0.099725246429443359,   0.070788159966468811,  -0.007806164678186178,  -0.059619396924972534,  -0.126831442117691040,   0.143088519573211670,   0.021540723741054535,   0.101359546184539795,   0.056007768958806992,  -0.095759615302085876,   0.023996233940124512,  -0.092107981443405151,   0.092801287770271301,  -0.089603282511234283,   0.072041861712932587,   0.080419033765792847,  -0.140492811799049377,  -0.115977823734283447,   0.021813483908772469,  -0.107547625899314880,  -0.130753800272941589,  -0.044669546186923981,  -0.086716637015342712,   0.027419794350862503,  -0.120626479387283325,   0.081244744360446930,  -0.071723744273185730,   0.083130948245525360,   0.022109676152467728,  -0.043856300413608551,   0.153579905629158020,  -0.117361903190612793,   0.068158291280269623,  -0.041222404688596725,
		   0.092483952641487122,  -0.066239729523658752,  -0.059233829379081726,  -0.141864746809005737,   0.065126277506351471,  -0.155846148729324341,   0.094671405851840973,  -0.138198539614677429,   0.036935646086931229,  -0.122915282845497131,   0.048747092485427856,  -0.102469787001609802,  -0.003933570813387632,  -0.075764179229736328,   0.031685825437307358,  -0.112088441848754883,  -0.129331782460212708,  -0.174569845199584961,  -0.113260813057422638,  -0.126738369464874268,  -0.045285139232873917,  -0.058227982372045517,  -0.150347515940666199,   0.122846953570842743,   0.070291347801685333,  -0.048078835010528564,  -0.039972290396690369,  -0.151261195540428162,   0.084025204181671143,  -0.069595642387866974,   0.056200262159109116,  -0.069189362227916718,   0.037046268582344055,  -0.090419359505176544,   0.113783806562423706,   0.081715814769268036,  -0.073547415435314178,   0.018785921856760979,  -0.072498120367527008,   0.039294857531785965,  -0.032751839607954025,   0.100009009242057800,  -0.093166291713714600,   0.137015715241432190,  -0.008254189044237137,  -0.043578486889600754,  -0.091023668646812439,   0.039947900921106339,
		   0.069562718272209167,   0.146123260259628296,   0.026207584887742996,   0.160829901695251465,   0.055178422480821609,  -0.064562119543552399,  -0.107801362872123718,  -0.026669068261981010,  -0.025184582918882370,  -0.030869428068399429,   0.029548613354563713,   0.132644340395927429,  -0.077327795326709747,   0.069046929478645325,   0.078930944204330444,   0.082299344241619110,  -0.041227735579013824,   0.154201552271842957,   0.049275182187557220,  -0.038262974470853806,   0.104214690625667572,  -0.083579391241073608,   0.046282611787319183,   0.037659645080566406,   0.052311971783638000,  -0.171531036496162415,   0.024836625903844833,  -0.082720525562763214,   0.192041963338851929,  -0.067609548568725586,   0.055197864770889282,   0.139220654964447021,  -0.092564187943935394,  -0.108710132539272308,   0.059413518756628036,   0.141365647315979004,   0.098898500204086304,  -0.096918463706970215,   0.160932362079620361,   0.022867003455758095,   0.111898288130760193,  -0.155686408281326294,   0.100630491971969604,   0.074435688555240631,  -0.085908249020576477,   0.141786605119705200,  -0.004860396496951580,  -0.025219688192009926,
		   0.025672568008303642,   0.100388228893280029,  -0.106381043791770935,  -0.095908686518669128,   0.094906248152256012,   0.013981085270643234,  -0.023400176316499710,   0.101397924125194550,  -0.079603098332881927,  -0.110545650124549866,  -0.005160337314009666,  -0.018452797085046768,  -0.055574167519807816,  -0.102804079651832581,   0.069106265902519226,  -0.078451067209243774,  -0.134522214531898499,  -0.121779859066009521,  -0.063175410032272339,  -0.084492243826389313,  -0.153201356530189514,   0.118746601045131683,   0.031467497348785400,  -0.108211375772953033,  -0.014800901524722576,   0.013306817039847374,  -0.115636460483074188,  -0.052033331245183945,  -0.056827042251825333,   0.115124143660068512,  -0.152046561241149902,   0.022988906130194664,  -0.105627827346324921,  -0.047224897891283035,   0.090173728764057159,  -0.014895408414304256,  -0.032647382467985153,  -0.142541006207466125,  -0.161549627780914307,  -0.139873281121253967,   0.096379891037940979,   0.076544292271137238,  -0.155445665121078491,  -0.109913751482963562,  -0.150639519095420837,  -0.060806378722190857,  -0.122168436646461487,   0.118828177452087402,
		   0.115447908639907837,  -0.009874944575130939,  -0.098617449402809143,  -0.022835331037640572,  -0.140049129724502563,  -0.095201954245567322,  -0.063580282032489777,  -0.066505074501037598,   0.015824604779481888,   0.136187136173248291,  -0.101248450577259064,  -0.095737874507904053,  -0.003086180193349719,  -0.096863433718681335,   0.083774350583553314,  -0.058454114943742752,  -0.007551420945674181,   0.096650160849094391,   0.027829581871628761,   0.030514296144247055,   0.147395655512809753,   0.080505020916461945,  -0.092228628695011139,  -0.139576107263565063,  -0.044232781976461411,  -0.091042548418045044,  -0.041823036968708038,   0.064094267785549164,   0.148409441113471985,   0.042635314166545868,  -0.095317512750625610,  -0.103676036000251770,   0.131187960505485535,  -0.015810025855898857,  -0.103081449866294861,   0.123063087463378906,  -0.111978985369205475,  -0.100550033152103424,   0.034072738140821457,  -0.006021082866936922,   0.020852494984865189,  -0.024753354489803314,   0.118783794343471527,  -0.137683361768722534,   0.158147588372230530,   0.052900299429893494,   0.005355429835617542,  -0.099617615342140198,
		   0.014846288599073887,  -0.122308067977428436,   0.109539188444614410,  -0.058318387717008591,   0.052688620984554291,  -0.040701266378164291,  -0.031523063778877258,  -0.029539391398429871,   0.140616506338119507,   0.029481671750545502,  -0.051213901489973068,  -0.095996916294097900,  -0.135076612234115601,  -0.037466466426849365,  -0.043740194290876389,  -0.005062276031821966,   0.023182943463325500,  -0.105780348181724548,   0.106798812747001648,   0.045549865812063217,  -0.123310655355453491,  -0.051462795585393906,   0.041066631674766541,  -0.098633341491222382,  -0.048735912889242172,  -0.101384527981281281,   0.044545292854309082,  -0.107866905629634857,   0.028486480936408043,  -0.102693676948547363,   0.083659291267395020,  -0.008920831605792046,  -0.090851701796054840,   0.025799522176384926,  -0.084804669022560120,   0.030255835503339767,   0.129433736205101013,   0.094747871160507202,  -0.014083375222980976,  -0.127770707011222839,   0.016065305098891258,   0.000136317670694552,  -0.092568002641201019,   0.041484236717224121,  -0.050781819969415665,  -0.057667158544063568,   0.136603474617004395,   0.108818218111991882,
		  -0.106491535902023315,  -0.022230708971619606,   0.091793343424797058,  -0.003480253275483847,   0.015406625345349312,   0.029342589899897575,  -0.136387705802917480,  -0.013909798115491867,   0.064389340579509735,  -0.150394260883331299,   0.140973612666130066,  -0.131531566381454468,  -0.088885322213172913,  -0.068316668272018433,  -0.031072096899151802,   0.005992752034217119,  -0.090620853006839752,   0.051235940307378769,  -0.129338264465332031,  -0.067058816552162170,  -0.159267216920852661,   0.011351029388606548,  -0.156659185886383057,   0.030870599672198296,   0.134041383862495422,  -0.114931359887123108,  -0.091569513082504272,   0.048289611935615540,  -0.141006290912628174,   0.040194101631641388,   0.122136242687702179,  -0.118551790714263916,   0.101472534239292145,  -0.024658069014549255,  -0.136262223124504089,  -0.104815550148487091,   0.077670149505138397,  -0.164778083562850952,  -0.040244441479444504,  -0.133141905069351196,  -0.118192389607429504,  -0.045021396130323410,  -0.146016344428062439,   0.075308673083782196,  -0.057356067001819611,   0.087296031415462494,   0.136963427066802979,  -0.015070493333041668,
		   0.133087307214736938,   0.109075501561164856,  -0.095141381025314331,  -0.037469714879989624,  -0.090238429605960846,   0.036219239234924316,   0.021721385419368744,   0.048222556710243225,   0.097662575542926788,   0.115678384900093079,   0.120904661715030670,  -0.160061284899711609,  -0.105997912585735321,  -0.091801770031452179,   0.138533398509025574,  -0.098721548914909363,  -0.011544829234480858,   0.038169234991073608,  -0.065767683088779449,  -0.091436505317687988,   0.110548838973045349,   0.128608018159866333,   0.054044518619775772,  -0.125568985939025879,   0.029768493026494980,  -0.019733970984816551,  -0.122274897992610931,   0.112569831311702728,   0.102608241140842438,   0.149015441536903381,   0.148025706410408020,  -0.059008423238992691,   0.013379237614572048,  -0.094038948416709900,   0.038209211081266403,  -0.021870281547307968,   0.083789288997650146,   0.033431731164455414,   0.115313462913036346,   0.084936574101448059,  -0.039705917239189148,   0.042851205915212631,   0.108376033604145050,  -0.096674688160419464,   0.136303916573524475,  -0.015172314830124378,   0.104250699281692505,  -0.038457125425338745,
		   0.119225494563579559,   0.128146559000015259,   0.144596710801124573,  -0.001114384969696403,   0.044206526130437851,  -0.051834836602210999,  -0.032426055520772934,  -0.013674642890691757,  -0.094405144453048706,  -0.036636125296354294,  -0.139341831207275391,   0.126853808760643005,  -0.097508467733860016,   0.074110753834247589,  -0.162398532032966614,   0.145925998687744141,   0.020201969891786575,  -0.058887209743261337,  -0.015882128849625587,  -0.024117974564433098,   0.056851640343666077,   0.083940722048282623,  -0.001900909352116287,   0.079136840999126434,  -0.070375673472881317,   0.080588236451148987,   0.126731321215629578,   0.031932976096868515,  -0.097485743463039398,  -0.086045950651168823,   0.110526017844676971,  -0.043521534651517868,   0.058948826044797897,  -0.029068250209093094,   0.000556153652723879,  -0.087012231349945068,  -0.082312121987342834,  -0.059700198471546173,   0.113576956093311310,  -0.002357303630560637,  -0.140696033835411072,  -0.043380618095397949,  -0.034299887716770172,  -0.076173238456249237,  -0.002027867594733834,   0.024220854043960571,  -0.063810206949710846,   0.091934502124786377,
		  -0.045751091092824936,  -0.137417450547218323,  -0.013327755033969879,   0.128260597586631775,   0.050722911953926086,  -0.060223221778869629,  -0.078703217208385468,   0.021691340953111649,  -0.056297115981578827,  -0.062755778431892395,  -0.107106499373912811,   0.012124123983085155,   0.193444937467575073,  -0.127108186483383179,   0.146754771471023560,  -0.070276215672492981,   0.192051291465759277,   0.196121677756309509,  -0.099952250719070435,  -0.086244694888591766,   0.166198447346687317,   0.013861330226063728,   0.023073159158229828,  -0.121377244591712952,   0.102305635809898376,  -0.158019557595252991,   0.068907767534255981,   0.083478651940822601,  -0.011770532466471195,  -0.133801087737083435,  -0.001199123798869550,  -0.032201837748289108,  -0.039015538990497589,  -0.063310749828815460,   0.128658771514892578,  -0.035480223596096039,  -0.022920068353414536,   0.064571313560009003,  -0.043340224772691727,   0.119744881987571716,   0.067885212600231171,  -0.134047463536262512,   0.133485600352287292,  -0.107728116214275360,   0.084646493196487427,   0.054487984627485275,  -0.031751655042171478,  -0.066568911075592041,
		   0.085351251065731049,   0.054077487438917160,  -0.107083670794963837,   0.077900856733322144,  -0.014518550597131252,  -0.061514850705862045,   0.064791440963745117,  -0.018587747588753700,   0.012096717953681946,   0.092748083174228668,   0.074322372674942017,  -0.089092113077640533,   0.120233058929443359,   0.045261558145284653,  -0.088982239365577698,  -0.006470632739365101,   0.181896001100540161,   0.015411888249218464,   0.082011476159095764,  -0.143170475959777832,   0.133115947246551514,   0.096357911825180054,  -0.037780106067657471,  -0.113837443292140961,   0.078279137611389160,   0.070700541138648987,   0.060203585773706436,   0.004361857194453478,  -0.077331438660621643,  -0.012429177761077881,  -0.030398437753319740,   0.002431664848700166,   0.158698767423629761,  -0.058842767030000687,  -0.072146050631999969,  -0.028678648173809052,  -0.010789034888148308,   0.084686450660228729,   0.139008715748786926,  -0.062221698462963104,   0.012384538538753986,  -0.114082954823970795,  -0.077139429748058319,   0.111069545149803162,  -0.058602306991815567,   0.117694132030010223,   0.066379509866237640,  -0.014758192002773285,
		   0.035255320370197296,   0.120876416563987732,   0.066450618207454681,   0.030949775129556656,  -0.060530554503202438,   0.077155798673629761,   0.128263920545578003,   0.109551772475242615,   0.076407566666603088,   0.010425882413983345,  -0.056949339807033539,   0.120985604822635651,  -0.104701928794384003,   0.160949602723121643,  -0.127892255783081055,   0.131845250725746155,   0.041441466659307480,  -0.089928381145000458,   0.104109212756156921,   0.060443550348281860,   0.002076450502499938,   0.011563984677195549,   0.032090317457914352,   0.154990807175636292,   0.023835493251681328,  -0.039517942816019058,   0.082348175346851349,  -0.020510958507657051,  -0.107471510767936707,   0.066706053912639618,  -0.137514889240264893,   0.077969923615455627,   0.131907984614372253,  -0.171367779374122620,   0.014435485005378723,   0.062023825943470001,  -0.011166173033416271,  -0.004514154046773911,  -0.014715563505887985,  -0.056589189916849136,  -0.028131255879998207,  -0.055995482951402664,   0.033782619982957840,   0.004649381153285503,  -0.163411110639572144,   0.078444823622703552,  -0.069284968078136444,  -0.039758823812007904,
		  -0.062537327408790588,   0.023462459444999695,  -0.104513928294181824,   0.101867832243442535,  -0.104541927576065063,   0.047449447214603424,  -0.095469810068607330,  -0.051879402250051498,   0.085151843726634979,  -0.033213011920452118,  -0.002088214969262481,   0.093304529786109924,   0.094048924744129181,  -0.078873552381992340,  -0.033573165535926819,  -0.031371824443340302,  -0.096121966838836670,  -0.136636421084403992,  -0.075898043811321259,   0.015483176335692406,   0.020253468304872513,  -0.033748824149370193,  -0.066920705139636993,   0.176718607544898987,  -0.139070317149162292,  -0.090571768581867218,  -0.061825811862945557,   0.139606043696403503,   0.000435803347500041,  -0.116766892373561859,  -0.132710680365562439,  -0.005429098848253489,   0.045310840010643005,   0.127565979957580566,   0.055424444377422333,   0.038921464234590530,   0.011964336968958378,   0.128761515021324158,  -0.017753738909959793,   0.006035054568201303,  -0.067838303744792938,   0.041006751358509064,  -0.045369572937488556,   0.057537395507097244,  -0.126669034361839294,   0.074273020029067993,  -0.112638920545578003,  -0.052727308124303818,
		  -0.098180226981639862,  -0.033197499811649323,  -0.062385182827711105,   0.089746147394180298,  -0.032193850725889206,   0.190204069018363953,   0.075918979942798615,   0.025555586442351341,   0.130808293819427490,   0.059739232063293457,  -0.004324795212596655,   0.132778480648994446,   0.053426638245582581,  -0.123067371547222137,   0.016140831634402275,   0.029527112841606140,   0.031056286767125130,   0.150436967611312866,   0.026001811027526855,  -0.142643839120864868,   0.040805142372846603,  -0.077378161251544952,  -0.085342250764369965,   0.065570123493671417,   0.072691336274147034,  -0.139306083321571350,   0.096225924789905548,  -0.063359282910823822,  -0.002517956076189876,   0.016357427462935448,  -0.029435172677040100,   0.167352989315986633,  -0.059439070522785187,   0.072047740221023560,  -0.068787880241870880,   0.035715844482183456,   0.134516850113868713,   0.015121790580451488,  -0.068559035658836365,   0.156030073761940002,   0.071173191070556641,  -0.173237249255180359,  -0.023067208006978035,   0.006570605561137199,   0.122362397611141205,  -0.071346692740917206,   0.066419184207916260,   0.022103218361735344,
		  -0.058229513466358185,   0.106708161532878876,  -0.003022040938958526,   0.072435729205608368,   0.045036211609840393,  -0.017494847998023033,  -0.114418029785156250,   0.133070126175880432,   0.088692791759967804,   0.110667847096920013,  -0.105392619967460632,  -0.119324237108230591,  -0.102814845740795135,   0.076806053519248962,  -0.106441669166088104,   0.050739485770463943,  -0.123773917555809021,   0.076652467250823975,  -0.129118800163269043,   0.111330710351467133,  -0.070019379258155823,  -0.118426255881786346,   0.154313996434211731,   0.025485340505838394,   0.036047391593456268,   0.109226696193218231,   0.104048438370227814,   0.032233968377113342,  -0.013893358409404755,   0.046530727297067642,   0.110703565180301666,  -0.043014720082283020,  -0.130199134349822998,   0.035401239991188049,   0.079668991267681122,   0.042732097208499908,  -0.069873653352260590,   0.017481463029980659,  -0.025724284350872040,   0.037715360522270203,   0.068530917167663574,  -0.103861212730407715,  -0.083806000649929047,   0.038140799850225449,  -0.003362316172569990,   0.152904719114303589,   0.126416295766830444,   0.066258594393730164,
		   0.025858469307422638,   0.021112881600856781,   0.061051104217767715,   0.093559302389621735,   0.093483984470367432,  -0.143618822097778320,   0.049312964081764221,   0.032491046935319901,  -0.052038554102182388,   0.169821918010711670,  -0.024467660114169121,  -0.064873978495597839,  -0.094593070447444916,  -0.072561010718345642,   0.063870660960674286,   0.147146269679069519,  -0.058218535035848618,  -0.005869075190275908,   0.007010714150965214,  -0.062979072332382202,   0.172093197703361511,  -0.078436449170112610,  -0.136261120438575745,   0.034084998071193695,  -0.140474870800971985,   0.056798025965690613,   0.005793415009975433,   0.036277789622545242,   0.160394623875617981,   0.044662062078714371,  -0.018062876537442207,   0.142273560166358948,  -0.085320949554443359,   0.046434190124273300,   0.146174237132072449,  -0.049030940979719162,   0.085957311093807220,   0.014638955704867840,   0.137180283665657043,  -0.020252877846360207,  -0.158089071512222290,   0.072391144931316376,   0.032570529729127884,   0.162118881940841675,  -0.006152836605906487,   0.082272008061408997,   0.098522089421749115,   0.084164328873157501,
		   0.090668886899948120,  -0.069899834692478180,  -0.093753784894943237,  -0.120756588876247406,  -0.039574503898620605,   0.040698587894439697,   0.041634723544120789,  -0.053325526416301727,  -0.032073996961116791,  -0.036412287503480911,   0.123282894492149353,  -0.141224831342697144,  -0.087568722665309906,  -0.062638536095619202,  -0.063662067055702209,  -0.003277165116742253,   0.120988130569458008,  -0.009595822542905807,  -0.137741178274154663,  -0.064697615802288055,  -0.081954665482044220,   0.153355792164802551,   0.022333480417728424,   0.024566348642110825,   0.014404564164578915,  -0.020331621170043945,   0.117896258831024170,   0.068073935806751251,  -0.105408482253551483,   0.090246804058551788,   0.068338125944137573,  -0.025961220264434814,   0.047522582113742828,   0.079249501228332520,   0.033727876842021942,  -0.155820772051811218,   0.041782621294260025,  -0.020235329866409302,  -0.069948539137840271,  -0.142826437950134277,  -0.025803528726100922,  -0.008602724410593510,   0.084629334509372711,  -0.152066960930824280,   0.054086305201053619,   0.074732057750225067,   0.124688141047954559,   0.116056405007839203,
	};
	static const double bias04[]=
	{
		  -0.034444153308868408,   0.027404276654124260,  -0.036175731569528580,   0.153692752122879028,   0.001608376624062657,  -0.068360283970832825,  -0.044887848198413849,  -0.039495918899774551,   0.040557734668254852,  -0.047557760030031204,   0.019699038937687874,  -0.004210551735013723,  -0.066258639097213745,  -0.057032875716686249,  -0.034439329057931900,  -0.109724089503288269,   0.144170597195625305,  -0.012520448304712772,   0.018250081688165665,   0.128400400280952454,   0.057754173874855042,  -0.039242010563611984,  -0.066304542124271393,  -0.051178019493818283,   0.025654457509517670,  -0.031091177836060524,   0.032859433442354202,   0.045359842479228973,   0.024232454597949982,   0.109065249562263489,   0.094312593340873718,  -0.109796755015850067,   0.021626975387334824,   0.135847553610801697,  -0.136965423822402954,  -0.098282918334007263,  -0.110864102840423584,  -0.098820827901363373,  -0.066214069724082947,  -0.082993738353252411,   0.031074538826942444,   0.002271912992000580,   0.102069042623043060,   0.028509162366390228,   0.050790410488843918,  -0.016996974125504494,   0.036369219422340393,   0.028901457786560059,
	};
	static const double weight05[]=
	{
		   0.125585615634918213,   0.057178370654582977,  -0.037143524736166000,  -0.148183599114418030,   0.108163826167583466,   0.126040458679199219,  -0.032907474786043167,   0.101936422288417816,   0.147323220968246460,  -0.081362538039684296,  -0.071125455200672150,   0.065193526446819305,  -0.142328992486000061,   0.023663360625505447,   0.077915385365486145,   0.086963593959808350,   0.002623175038024783,   0.082515910267829895,  -0.023572651669383049,  -0.009552929550409317,   0.028940835967659950,  -0.140132412314414978,  -0.095330171287059784,  -0.056020982563495636,  -0.016832133755087852,  -0.078386075794696808,  -0.134443417191505432,  -0.090023018419742584,   0.151046961545944214,  -0.018124068155884743,  -0.064721018075942993,   0.136211350560188293,   0.007347978185862303,   0.011923289857804775,  -0.037623871117830276,   0.000898086116649210,  -0.122598059475421906,  -0.099067762494087219,   0.012380475178360939,  -0.114957720041275024,  -0.008117969147861004,   0.149953275918960571,   0.138880386948585510,   0.080454312264919281,   0.024211315438151360,  -0.048561319708824158,  -0.025662949308753014,  -0.032096616923809052,
		  -0.011726432479918003,  -0.011031988076865673,  -0.152882456779479980,   0.059770118445158005,   0.001279058633372188,   0.148346662521362305,  -0.023380858823657036,   0.153980165719985962,  -0.031399250030517578,  -0.108930654823780060,  -0.082779407501220703,  -0.111783243715763092,  -0.002514134859666228,   0.080042250454425812,  -0.130431026220321655,   0.128613486886024475,   0.086342304944992065,   0.016077786684036255,  -0.009358510375022888,   0.078530192375183105,  -0.069386564195156097,   0.043038286268711090,   0.132813438773155212,   0.042305208742618561,   0.000457951333373785,   0.146394416689872742,   0.000356174365151674,   0.034771714359521866,   0.077467009425163269,   0.168037056922912598,   0.097704641520977020,   0.030617792159318924,  -0.113397739827632904,   0.157467395067214966,  -0.115581572055816650,   0.111254975199699402,  -0.018462637439370155,  -0.141461640596389771,   0.135403797030448914,   0.083950281143188477,  -0.037774864584207535,   0.050788190215826035,  -0.086005732417106628,  -0.178316891193389893,   0.099144786596298218,   0.082890182733535767,   0.050035625696182251,   0.014037004671990871,
		   0.157460510730743408,   0.114752158522605896,   0.116324380040168762,  -0.037613540887832642,  -0.034449435770511627,  -0.063563555479049683,  -0.149055242538452148,   0.032171677798032761,   0.033044721931219101,  -0.039136875420808792,  -0.096675485372543335,   0.038317229598760605,  -0.152722284197807312,  -0.030523175373673439,   0.084018401801586151,  -0.030766880139708519,  -0.046113029122352600,  -0.114737607538700104,   0.163884043693542480,  -0.086627855896949768,  -0.048763025552034378,  -0.073071323335170746,  -0.023329099640250206,   0.053827520459890366,   0.044491719454526901,   0.043085832148790359,  -0.095412991940975189,  -0.058477431535720825,  -0.076983034610748291,  -0.040327537804841995,  -0.001097491360269487,  -0.135698795318603516,  -0.043009210377931595,  -0.045651610940694809,   0.087853588163852692,  -0.000559232546947896,  -0.133495405316352844,   0.095446579158306122,  -0.107100322842597961,   0.080326721072196960,  -0.130634680390357971,   0.157731458544731140,   0.082700751721858978,   0.132038101553916931,  -0.040633007884025574,   0.027311602607369423,  -0.022705480456352234,   0.104510232806205750,
		  -0.051885403692722321,   0.031760018318891525,   0.034369993954896927,  -0.075386248528957367,   0.083242028951644897,   0.058015409857034683,  -0.018395565450191498,  -0.158811464905738831,   0.013190739788115025,   0.015626961365342140,   0.023896865546703339,   0.078631155192852020,  -0.049640018492937088,   0.017426507547497749,  -0.109715804457664490,   0.068953745067119598,  -0.008571649901568890,  -0.106354936957359314,   0.081364989280700684,  -0.033580984920263290,  -0.137048110365867615,   0.008960584178566933,  -0.125552371144294739,   0.082128986716270447,  -0.095054574310779572,  -0.053056962788105011,   0.090242698788642883,  -0.181419923901557922,  -0.025024149566888809,  -0.188564017415046692,   0.014320423826575279,   0.129214346408843994,  -0.117582440376281738,  -0.022322507575154305,   0.059720069169998169,   0.083478093147277832,   0.123359076678752899,   0.081391520798206329,  -0.100920483469963074,  -0.131456300616264343,  -0.058086667209863663,  -0.006701818201690912,   0.109357535839080811,   0.030467681586742401,   0.091840736567974091,  -0.047328736633062363,  -0.022139305248856544,  -0.107122525572776794,
		   0.180136278271675110,   0.106489092111587524,   0.071051999926567078,  -0.035387594252824783,   0.051471095532178879,   0.061691395938396454,   0.085223212838172913,  -0.104306332767009735,   0.092263273894786835,  -0.087298557162284851,   0.068414598703384399,  -0.085087694227695465,   0.054826624691486359,   0.039031956344842911,   0.066685721278190613,   0.117552764713764191,  -0.072210215032100677,   0.135759502649307251,   0.122349262237548828,  -0.147376596927642822,   0.031090680509805679,  -0.099669843912124634,  -0.016187086701393127,   0.104530043900012970,  -0.007636271417140961,   0.141284793615341187,   0.068155728280544281,   0.092174537479877472,   0.057369910180568695,  -0.167590916156768799,   0.014103075489401817,  -0.159979462623596191,  -0.072859913110733032,   0.155906781554222107,   0.078104540705680847,  -0.048601154237985611,   0.130407884716987610,  -0.096529617905616760,  -0.013614736497402191,  -0.075788497924804688,  -0.109098337590694427,  -0.125534847378730774,  -0.036236669868230820,   0.023259848356246948,  -0.081587061285972595,  -0.034969199448823929,  -0.007564601488411427,  -0.092126086354255676,
		  -0.085677474737167358,  -0.079364843666553497,   0.099442072212696075,   0.063343018293380737,   0.013053916394710541,  -0.147721841931343079,   0.016461668536067009,  -0.058715533465147018,  -0.067880623042583466,   0.016990834847092628,   0.095580615103244781,  -0.130455106496810913,  -0.065777786076068878,  -0.044212359935045242,   0.106342457234859467,  -0.128977045416831970,  -0.089042134582996368,   0.062873527407646179,  -0.042823821306228638,  -0.060361586511135101,  -0.101782374083995819,   0.075499005615711212,   0.081589609384536743,   0.147838249802589417,  -0.113480456173419952,  -0.171182006597518921,  -0.151556834578514099,  -0.140721365809440613,   0.115819200873374939,  -0.092399694025516510,  -0.138850241899490356,  -0.054383426904678345,   0.150535956025123596,   0.006294653750956059,  -0.032665420323610306,   0.028097620233893394,  -0.005342203658074141,  -0.031511731445789337,   0.010827429592609406,  -0.022390810772776604,  -0.125465348362922668,  -0.003404526971280575,   0.079232320189476013,   0.074865408241748810,  -0.048697762191295624,   0.005498645827174187,  -0.107967898249626160,  -0.015328740701079369,
		   0.140798717737197876,  -0.125318318605422974,   0.120835900306701660,   0.048341050744056702,   0.030755762010812759,   0.060429342091083527,   0.065472155809402466,   0.047089319676160812,  -0.063600368797779083,   0.063363954424858093,   0.134347632527351379,   0.094810225069522858,   0.158702522516250610,  -0.029304377734661102,  -0.077617563307285309,   0.001265107071958482,  -0.037771172821521759,  -0.115276433527469635,  -0.100223056972026825,   0.085243523120880127,  -0.002567852614447474,  -0.080177031457424164,   0.137790888547897339,  -0.106131814420223236,  -0.018006106838583946,   0.146753281354904175,  -0.012174411676824093,   0.197719722986221313,  -0.001195163582451642,   0.190924093127250671,   0.016947852447628975,  -0.146310567855834961,   0.060288175940513611,  -0.038849312812089920,  -0.076616130769252777,   0.168737068772315979,  -0.036749944090843201,   0.030200280249118805,  -0.025606188923120499,   0.114340268075466156,   0.022441776469349861,   0.024192858487367630,  -0.082002565264701843,   0.087533205747604370,  -0.069247692823410034,   0.018784057348966599,  -0.049646306782960892,  -0.044621992856264114,
		  -0.052781950682401657,  -0.142469629645347595,  -0.001223539700731635,   0.136766508221626282,  -0.010191168636083603,   0.170816406607627869,   0.108629822731018066,   0.056798379868268967,   0.096079505980014801,   0.097789824008941650,   0.055975403636693954,  -0.147377103567123413,  -0.061025999486446381,   0.110389113426208496,   0.102214261889457703,   0.063362859189510345,  -0.011397991329431534,   0.042301621288061142,   0.070637919008731842,   0.081695117056369781,  -0.093959808349609375,  -0.005928518250584602,   0.077269583940505981,  -0.117894016206264496,   0.104884438216686249,   0.009306121617555618,   0.089529491961002350,   0.061505489051342010,   0.149606063961982727,   0.075216606259346008,   0.095311179757118225,   0.106502041220664978,  -0.148449957370758057,  -0.109624087810516357,  -0.071064770221710205,  -0.066360905766487122,  -0.137534856796264648,   0.069384038448333740,  -0.069762080907821655,   0.070257082581520081,   0.124887876212596893,   0.074650652706623077,  -0.009291936643421650,   0.005106806289404631,   0.031786330044269562,   0.076378978788852692,   0.026365382596850395,   0.042275112122297287,
		   0.050769492983818054,   0.026666173711419106,   0.076389081776142120,  -0.104722790420055389,  -0.066614694893360138,  -0.176694005727767944,  -0.068073518574237823,  -0.021895920857787132,  -0.089909180998802185,   0.024854471907019615,  -0.004368834663182497,   0.137205600738525391,  -0.087532624602317810,   0.044192451983690262,   0.023852841928601265,   0.106274515390396118,   0.033178783953189850,  -0.128207653760910034,  -0.132931426167488098,  -0.018092162907123566,  -0.056266441941261292,   0.013913610018789768,  -0.091511897742748260,  -0.078701458871364594,   0.099280171096324921,  -0.008302526548504829,   0.006826364900916815,  -0.186266869306564331,  -0.136227190494537354,   0.020594913512468338,  -0.145748823881149292,  -0.113904818892478943,  -0.037832204252481461,   0.003448658855631948,   0.042440589517354965,  -0.121404618024826050,  -0.078886322677135468,   0.106888368725776672,  -0.114423289895057678,   0.122511334717273712,  -0.179403796792030334,  -0.026080932468175888,  -0.103745996952056885,  -0.061912279576063156,   0.000405208702431992,  -0.045489162206649780,  -0.059087458997964859,   0.052682667970657349,
		  -0.018315551802515984,   0.004777090623974800,   0.072523713111877441,  -0.016434332355856895,   0.014853293076157570,  -0.010327761992812157,   0.055925160646438599,   0.001269589993171394,   0.124392680823802948,   0.044071882963180542,  -0.017728354781866074,   0.091559715569019318,  -0.136294096708297729,  -0.128766700625419617,  -0.022383691743016243,  -0.016563188284635544,   0.005702503025531769,   0.026498714461922646,   0.103892378509044647,  -0.134876847267150879,  -0.031288761645555496,  -0.137049853801727295,   0.018654158338904381,   0.017147120088338852,  -0.149398043751716614,   0.095543593168258667,  -0.091795288026332855,  -0.176921606063842773,  -0.102475725114345551,  -0.197823718190193176,  -0.175586059689521790,  -0.082991123199462891,  -0.129542499780654907,  -0.085900284349918365,   0.062240786850452423,  -0.004285248462110758,   0.111034981906414032,   0.024723233655095100,  -0.154487743973731995,   0.052522201091051102,  -0.084124609827995300,   0.060684002935886383,  -0.128109529614448547,  -0.104313954710960388,  -0.096147909760475159,   0.037222567945718765,  -0.137707009911537170,   0.145736709237098694,
		   0.031129864975810051,  -0.134106457233428955,   0.051083471626043320,   0.066992692649364471,  -0.062123391777276993,  -0.054451875388622284,  -0.090152181684970856,   0.079097852110862732,  -0.088758327066898346,  -0.009056517854332924,  -0.137096390128135681,  -0.074819050729274750,  -0.112162284553050995,   0.131816834211349487,  -0.042809918522834778,  -0.008949873037636280,  -0.198130324482917786,  -0.077808342874050140,  -0.128335654735565186,   0.000221190915908664,  -0.050490260124206543,  -0.017559882253408432,  -0.072464741766452789,   0.007993056438863277,  -0.001242097700014710,  -0.020127225667238235,   0.103882759809494019,  -0.116760544478893280,   0.045207679271697998,  -0.007762932684272528,  -0.189331322908401489,   0.052429147064685822,  -0.027017833665013313,  -0.039524547755718231,  -0.128396749496459961,  -0.119844749569892883,   0.099114187061786652,  -0.047914605587720871,  -0.071649752557277679,   0.137338623404502869,   0.073267281055450439,   0.027866955846548080,  -0.010866163298487663,   0.067051991820335388,   0.078261837363243103,  -0.079988203942775726,   0.105784490704536438,   0.091369897127151489,
		  -0.147379189729690552,  -0.082653708755970001,   0.000457040732726455,   0.117540150880813599,   0.133002310991287231,   0.185643553733825684,  -0.069689385592937469,   0.013320460915565491,   0.022239426150918007,   0.090392783284187317,   0.141093641519546509,   0.036773730069398880,   0.127559959888458252,   0.039208285510540009,  -0.107687585055828094,   0.013243934139609337,   0.180552527308464050,   0.055848430842161179,  -0.110865272581577301,   0.143103316426277161,  -0.006948480382561684,   0.097705751657485962,   0.027010148391127586,  -0.022845443338155746,  -0.031804468482732773,  -0.088140219449996948,   0.123923361301422119,   0.075106933712959290,  -0.055791288614273071,   0.053255144506692886,   0.184874102473258972,   0.085185781121253967,  -0.019542325288057327,   0.170642808079719543,  -0.067702159285545349,   0.001671308069489896,  -0.099225461483001709,  -0.069849081337451935,  -0.037644740194082260,  -0.196530297398567200,   0.186816930770874023,   0.040651522576808929,   0.040727093815803528,   0.036705404520034790,   0.165190935134887695,  -0.122881971299648285,  -0.105775639414787292,  -0.063607372343540192,
		  -0.037275917828083038,  -0.147275775671005249,  -0.040902737528085709,   0.000068553599703591,  -0.056194737553596497,  -0.019670424982905388,   0.094844996929168701,  -0.096824131906032562,  -0.044123653322458267,  -0.136558532714843750,  -0.029781762510538101,  -0.171907752752304077,   0.158201664686203003,   0.025726409628987312,  -0.175784930586814880,  -0.117710456252098083,   0.119101867079734802,  -0.023536689579486847,  -0.190884083509445190,  -0.015540624968707561,  -0.016143217682838440,  -0.129226982593536377,  -0.003008216852322221,  -0.135960832238197327,  -0.113590769469738007,   0.087761290371417999,   0.151346370577812195,   0.095455057919025421,   0.019793421030044556,  -0.026986699551343918,   0.153886333107948303,   0.133918032050132751,  -0.103632025420665741,   0.128703668713569641,  -0.090551108121871948,   0.042114749550819397,  -0.111462980508804321,  -0.070558428764343262,   0.168185651302337646,   0.002256345935165882,   0.138178274035453796,   0.026924790814518929,   0.002932841889560223,  -0.017832640558481216,   0.145930454134941101,   0.086732380092144012,  -0.192618578672409058,  -0.031939554959535599,
		   0.024576568976044655,  -0.039037946611642838,  -0.022449407726526260,   0.015543202869594097,  -0.098798476159572601,   0.078827537596225739,   0.061153635382652283,  -0.036167491227388382,   0.003288861131295562,   0.128205239772796631,  -0.057876564562320709,   0.116046100854873657,  -0.155870720744132996,   0.011057233437895775,  -0.099483393132686615,   0.088887393474578857,   0.005075735505670309,   0.090609960258007050,  -0.017001686617732048,  -0.031376093626022339,  -0.118295937776565552,  -0.056818135082721710,   0.009434351697564125,   0.117416232824325562,   0.087280213832855225,  -0.042506255209445953,   0.121546275913715363,  -0.047694616019725800,  -0.024908587336540222,   0.013750393874943256,  -0.171443998813629150,   0.143221735954284668,  -0.067749790847301483,  -0.111116737127304077,   0.056019004434347153,   0.107477433979511261,  -0.087879568338394165,  -0.026214735582470894,   0.081911675631999969,  -0.018999267369508743,  -0.035581879317760468,  -0.042735490947961807,  -0.026101550087332726,  -0.016778413206338882,   0.025750976055860519,   0.090596921741962433,   0.008143696002662182,  -0.053078301250934601,
		  -0.007843360304832458,   0.105244934558868408,   0.057933412492275238,   0.010431133210659027,   0.044634480029344559,   0.101908557116985321,   0.059566043317317963,  -0.019571902230381966,   0.156432166695594788,   0.014519696123898029,  -0.068793259561061859,  -0.011767229065299034,  -0.128868564963340759,  -0.072348177433013916,   0.061591453850269318,   0.175337910652160645,   0.035597998648881912,   0.057201143354177475,  -0.082182087004184723,  -0.060109231621026993,   0.001811491092666984,   0.014187700115144253,   0.058455564081668854,   0.031980086117982864,   0.074092954397201538,  -0.059816177934408188,   0.024769796058535576,   0.126899808645248413,   0.141450688242912292,  -0.074622757732868195,  -0.196244403719902039,   0.119866475462913513,  -0.056177776306867599,  -0.086161993443965912,  -0.091028951108455658,   0.084072597324848175,  -0.078275851905345917,  -0.090245954692363739,  -0.106112405657768250,  -0.050225969403982162,  -0.179506316781044006,   0.030765468254685402,   0.141865193843841553,  -0.143162980675697327,   0.153510704636573792,  -0.089181356132030487,  -0.078952945768833160,  -0.163129478693008423,
		   0.167726084589958191,   0.146039277315139771,   0.071779780089855194,   0.126647949218750000,  -0.180811405181884766,   0.014293114654719830,   0.131778135895729065,  -0.128327205777168274,   0.186530828475952148,   0.038405872881412506,   0.147976458072662354,   0.064770810306072235,  -0.163307279348373413,   0.069366410374641418,  -0.083280459046363831,   0.176655128598213196,   0.029994184151291847,   0.185988843441009521,  -0.068398833274841309,  -0.028606496751308441,   0.151079639792442322,  -0.120636656880378723,   0.013445116579532623,   0.003472589887678623,   0.056845698505640030,   0.113752558827400208,  -0.166576102375984192,  -0.004627510905265808,  -0.003078051377087831,  -0.167119458317756653,   0.085989885032176971,  -0.019443823024630547,   0.057102661579847336,   0.143656373023986816,  -0.058482918888330460,   0.071562260389328003,  -0.010430218651890755,  -0.049588367342948914,  -0.015666471794247627,   0.160340949892997742,  -0.185723215341567993,   0.060708213597536087,   0.107729107141494751,   0.083875626325607300,   0.049315936863422394,  -0.009077825583517551,   0.028295684605836868,  -0.047458603978157043,
		   0.075125545263290405,   0.174168810248374939,   0.124034605920314789,   0.133078038692474365,  -0.157269075512886047,   0.011284904554486275,  -0.089267984032630920,  -0.085852533578872681,   0.133910968899726868,  -0.116166986525058746,  -0.034119453281164169,   0.072883136570453644,  -0.018400406464934349,  -0.077420130372047424,   0.177414238452911377,  -0.036348689347505569,  -0.107406169176101685,   0.115729123353958130,   0.174830406904220581,  -0.144320756196975708,   0.092751830816268921,   0.026744430884718895,   0.106209576129913330,   0.121496051549911499,   0.156547546386718750,  -0.020093498751521111,   0.058282516896724701,   0.140775069594383240,   0.089018926024436951,  -0.102430455386638641,  -0.157596543431282043,  -0.053214225918054581,   0.097996607422828674,   0.134332194924354553,  -0.123576894402503967,   0.060237351804971695,   0.032801952213048935,  -0.075334534049034119,  -0.115099467337131500,   0.178195387125015259,  -0.168771550059318542,   0.161000788211822510,   0.049799177795648575,   0.037968069314956665,   0.018474418669939041,  -0.048160649836063385,   0.145098313689231873,  -0.120759062469005585,
		  -0.030278928577899933,  -0.077625781297683716,   0.059721477329730988,   0.139779493212699890,  -0.097917824983596802,  -0.043708648532629013,  -0.080174162983894348,   0.016209470108151436,   0.062406182289123535,   0.114931993186473846,  -0.020470432937145233,  -0.026451729238033295,  -0.010722359642386436,  -0.027440203353762627,   0.066700324416160583,   0.149239748716354370,  -0.195851624011993408,   0.088742740452289581,   0.052374269813299179,  -0.031582642346620560,  -0.022890524938702583,   0.011143542826175690,  -0.109791480004787445,   0.059891376644372940,   0.125088497996330261,   0.076851703226566315,   0.094208620488643646,   0.072448648512363434,   0.068248197436332703,  -0.164244830608367920,  -0.141039431095123291,  -0.063564330339431763,  -0.084390304982662201,  -0.097468085587024689,  -0.134001657366752625,   0.109135806560516357,   0.051207870244979858,   0.032764766365289688,  -0.137867629528045654,   0.095901146531105042,   0.004417451098561287,  -0.063501589000225067,  -0.027212014421820641,   0.115847013890743256,  -0.069795139133930206,   0.160657450556755066,   0.040965292602777481,   0.058068808168172836,
		  -0.039777785539627075,   0.059236053377389908,  -0.012095023877918720,  -0.023212485015392303,  -0.081164777278900146,  -0.158575639128684998,  -0.139935940504074097,  -0.143529206514358521,   0.122408919036388397,   0.080884918570518494,  -0.019503287971019745,   0.118150390684604645,   0.110093824565410614,   0.048775635659694672,  -0.029976654797792435,  -0.015015087090432644,   0.045646544545888901,   0.073851056396961212,   0.077569946646690369,  -0.065559156239032745,   0.006806406192481518,  -0.011772188358008862,  -0.145194649696350098,  -0.139664962887763977,   0.009834774769842625,  -0.067699238657951355,  -0.085639506578445435,   0.030930496752262115,  -0.116444438695907593,   0.064728550612926483,   0.092299453914165497,  -0.061753802001476288,  -0.128886684775352478,   0.116834402084350586,   0.100857317447662354,   0.089925289154052734,  -0.045693054795265198,  -0.031993243843317032,  -0.111370667815208435,   0.127068638801574707,   0.010054912418127060,  -0.071528106927871704,   0.112031720578670502,  -0.038076531141996384,   0.066555485129356384,  -0.063103131949901581,   0.080676153302192688,   0.106086991727352142,
		   0.093865387141704559,   0.137715965509414673,  -0.118534833192825317,  -0.000802837894298136,   0.103553012013435364,   0.000738535716664046,   0.098521165549755096,   0.024914437904953957,  -0.062182348221540451,  -0.046887621283531189,  -0.032683562487363815,   0.028068717569112778,   0.122313387691974640,   0.078802391886711121,  -0.043671116232872009,  -0.047608479857444763,   0.169516533613204956,  -0.009708138182759285,   0.022663451731204987,   0.189101994037628174,  -0.117349773645401001,  -0.066926799714565277,  -0.095870509743690491,  -0.061821471899747849,   0.065472759306430817,  -0.050125777721405029,   0.053703725337982178,   0.044096145778894424,   0.135556712746620178,   0.140275374054908752,  -0.030146038159728050,  -0.013471880927681923,   0.060800541192293167,  -0.033098012208938599,  -0.069348447024822235,   0.002986596664413810,   0.026281744241714478,  -0.051050685346126556,   0.186794683337211609,  -0.152489215135574341,   0.180582791566848755,   0.017146281898021698,   0.067037433385848999,   0.044303700327873230,   0.049189332872629166,   0.072139695286750793,   0.052216336131095886,   0.125103592872619629,
		   0.040130637586116791,  -0.084394618868827820,  -0.059074882417917252,  -0.003200456267222762,   0.039499927312135696,   0.053741220384836197,  -0.123342752456665039,  -0.147207587957382202,   0.067485988140106201,  -0.122934311628341675,   0.120572425425052643,   0.183704316616058350,   0.047454908490180969,  -0.074172511696815491,  -0.071705989539623260,   0.091236084699630737,   0.006894421298056841,   0.020631533116102219,   0.161829799413681030,  -0.137319415807723999,   0.168973952531814575,   0.153807833790779114,  -0.109980843961238861,   0.020898731425404549,   0.065386891365051270,   0.045843731611967087,   0.052050471305847168,  -0.121021598577499390,  -0.043362271040678024,  -0.075443051755428314,  -0.057407006621360779,  -0.046041384339332581,   0.049959212541580200,   0.121292620897293091,   0.098061069846153259,  -0.095160990953445435,   0.008570769801735878,  -0.054884679615497589,  -0.160716682672500610,  -0.027380900457501411,  -0.063503161072731018,  -0.079279445111751556,  -0.093378446996212006,   0.002815669635310769,   0.054034091532230377,   0.059944476932287216,   0.175430476665496826,  -0.102335728704929352,
		  -0.125589087605476379,   0.040918417274951935,  -0.001734985155053437,   0.161606714129447937,   0.063707098364830017,   0.188694879412651062,  -0.027966460213065147,   0.127806618809700012,   0.042828779667615891,   0.068675592541694641,  -0.079512417316436768,  -0.084598258137702942,   0.062600068747997284,  -0.097823262214660645,   0.070102311670780182,  -0.099327892065048218,   0.188443183898925781,  -0.071115039288997650,  -0.061846774071455002,  -0.059606216847896576,   0.024685515090823174,   0.043632093816995621,   0.080104842782020569,  -0.130897313356399536,   0.088612474501132965,   0.015866564586758614,  -0.038368985056877136,   0.099862165749073029,   0.141209810972213745,   0.169755473732948303,   0.049224533140659332,  -0.063170507550239563,  -0.143490761518478394,  -0.028323210775852203,  -0.077841252088546753,   0.161718770861625671,  -0.149380400776863098,  -0.092040643095970154,  -0.062939524650573730,  -0.181306257843971252,   0.175913855433464050,  -0.044613301753997803,   0.010033825412392616,  -0.011124053969979286,  -0.024736154824495316,  -0.094662003219127655,  -0.088706880807876587,  -0.055320922285318375,
		   0.045211374759674072,   0.106290742754936218,   0.058591336011886597,   0.054282035678625107,   0.017175409942865372,  -0.060199525207281113,  -0.017676820978522301,   0.100274167954921722,  -0.055274538695812225,   0.057004492729902267,   0.068351894617080688,   0.001770490547642112,   0.050896625965833664,  -0.031327843666076660,  -0.045463573187589645,   0.047707226127386093,  -0.178897917270660400,   0.024209061637520790,   0.044218771159648895,  -0.127810031175613403,  -0.142206475138664246,  -0.045950710773468018,  -0.148815184831619263,   0.100681155920028687,  -0.025868387892842293,   0.034498773515224457,  -0.104610711336135864,  -0.127137020230293274,  -0.032807115465402603,  -0.077244400978088379,   0.039492592215538025,  -0.062562592327594757,   0.092362016439437866,  -0.067501083016395569,   0.066720895469188690,   0.007995228283107281,   0.152508929371833801,  -0.058396924287080765,  -0.033104714006185532,  -0.098699860274791718,  -0.139303565025329590,  -0.030336247757077217,  -0.150966912508010864,  -0.115717105567455292,  -0.095952175557613373,  -0.046546351164579391,  -0.150882661342620850,   0.072545349597930908,
		  -0.177483826875686646,  -0.149706333875656128,  -0.024498887360095978,  -0.053497135639190674,  -0.022111227735877037,  -0.141301691532135010,  -0.042828213423490524,   0.105855651199817657,   0.106276609003543854,   0.130869343876838684,   0.071855604648590088,   0.064433179795742035,  -0.084933616220951080,  -0.136181473731994629,   0.062372654676437378,  -0.126466199755668640,   0.051368366926908493,   0.031505584716796875,  -0.101040273904800415,  -0.122183464467525482,  -0.080076120793819427,  -0.071878768503665924,   0.084030948579311371,  -0.076208636164665222,  -0.051107488572597504,  -0.149220064282417297,   0.066554032266139984,  -0.138709962368011475,  -0.125826239585876465,   0.013028004206717014,  -0.133254304528236389,  -0.067002266645431519,  -0.057784833014011383,  -0.155475348234176636,  -0.078713044524192810,  -0.114619821310043335,  -0.084381200373172760,   0.066217936575412750,   0.030591411516070366,   0.118639558553695679,  -0.148413091897964478,   0.080951601266860962,   0.004795453511178493,   0.114998944103717804,   0.040118161588907242,   0.030168319121003151,   0.016566803678870201,  -0.002171367639675736,
	};
	static const double bias05[]=
	{
		   0.090784311294555664,   0.019710289314389229,  -0.043645255267620087,  -0.164307028055191040,   0.051259096711874008,   0.010963270440697670,  -0.036646477878093719,  -0.069213874638080597,  -0.036129213869571686,  -0.108621194958686829,  -0.188859254121780396,   0.167898833751678467,  -0.097859546542167664,   0.080491960048675537,  -0.065500006079673767,   0.068490609526634216,   0.041737735271453857,   0.099689178168773651,  -0.103179670870304108,   0.163395911455154419,  -0.073220990598201752,  -0.014758344739675522,  -0.068251609802246094,  -0.167497798800468445,
	};
	static const double weight06[]=
	{
		  -0.061573594808578491,   0.154305636882781982,   0.121408551931381226,  -0.113109566271305084,   0.010432111099362373,  -0.165441378951072693,   0.076593741774559021,   0.099624186754226685,   0.193615153431892395,  -0.012113876640796661,   0.127756491303443909,  -0.053349439054727554,  -0.163795903325080872,   0.149565964937210083,   0.090389527380466461,   0.003588597290217876,   0.112133182585239410,   0.223607361316680908,  -0.041132628917694092,   0.077351234853267670,   0.168130472302436829,  -0.069547861814498901,  -0.090242177248001099,  -0.120232529938220978,
		  -0.001261305413208902,   0.127436086535453796,  -0.028255771845579147,  -0.059230152517557144,  -0.204355612397193909,  -0.049131426960229874,   0.226905897259712219,   0.186723485589027405,  -0.182638645172119141,   0.161035433411598206,   0.126456573605537415,   0.060945607721805573,   0.210241749882698059,  -0.167275220155715942,  -0.199401646852493286,  -0.134118527173995972,   0.111528173089027405,   0.144385293126106262,   0.166834518313407898,   0.046669110655784607,  -0.048894722014665604,  -0.061065085232257843,  -0.176716178655624390,  -0.222040265798568726,
		  -0.011088541708886623,  -0.128197029232978821,  -0.181330069899559021,   0.106388464570045471,  -0.022323770448565483,  -0.220353379845619202,  -0.142507642507553101,  -0.182143509387969971,   0.184143394231796265,  -0.058546148240566254,   0.085986010730266571,   0.222737237811088562,   0.040017131716012955,  -0.004326886497437954,  -0.124426476657390594,  -0.117237947881221771,  -0.074008978903293610,   0.053895652294158936,   0.002546626143157482,   0.197631672024726868,  -0.047611098736524582,  -0.017905825749039650,  -0.135146006941795349,   0.030470019206404686,
		   0.079169809818267822,  -0.143215939402580261,  -0.096176505088806152,  -0.045727506279945374,   0.091287471354007721,  -0.046112634241580963,   0.065584920346736908,  -0.102312967181205750,  -0.118180081248283386,  -0.063614085316658020,  -0.225180074572563171,   0.015646765008568764,  -0.188164770603179932,  -0.065719410777091980,   0.073604919016361237,   0.159789487719535828,   0.229025542736053467,   0.010915385559201241,   0.009396781213581562,   0.038110788911581039,  -0.017489137127995491,  -0.174215018749237061,  -0.123276114463806152,  -0.155011832714080811,
		   0.225129172205924988,  -0.035636004060506821,   0.021994648501276970,  -0.107074521481990814,  -0.163158923387527466,   0.063356012105941772,  -0.033571500331163406,   0.113842211663722992,  -0.053124886006116867,  -0.152895718812942505,  -0.174849838018417358,   0.118203073740005493,  -0.111074380576610565,   0.031576290726661682,  -0.102412737905979156,  -0.022436678409576416,  -0.168892741203308105,  -0.124616771936416626,  -0.174093276262283325,   0.266223430633544922,  -0.015367407351732254,  -0.061702407896518707,  -0.048704113811254501,   0.031966742128133774,
		  -0.172517806291580200,   0.219090446829795837,  -0.082572102546691895,  -0.142388716340065002,  -0.207209318876266479,  -0.102501623332500458,   0.031172433868050575,   0.129124924540519714,  -0.180039197206497192,   0.012257676571607590,  -0.236455649137496948,   0.083452902734279633,  -0.123277880251407623,   0.065318211913108826,  -0.156456649303436279,  -0.175824195146560669,  -0.082165554165840149,   0.142818346619606018,  -0.082591779530048370,  -0.008980558253824711,   0.066251680254936218,   0.138596981763839722,   0.031421314924955368,  -0.197153344750404358,
		  -0.118509300053119659,  -0.085765585303306580,   0.176450416445732117,   0.110088042914867401,   0.225236952304840088,  -0.093877412378787994,   0.156265631318092346,   0.069452188909053802,  -0.078694954514503479,   0.045207347720861435,   0.047848060727119446,  -0.242250844836235046,   0.058742754161357880,  -0.169319227337837219,   0.142582148313522339,   0.017487978562712669,   0.195289954543113708,   0.244031891226768494,   0.232528954744338989,   0.110138863325119019,   0.162623301148414612,  -0.120391406118869781,   0.202706277370452881,   0.133640512824058533,
		   0.108704194426536560,  -0.004796571563929319,   0.189284205436706543,   0.087469041347503662,   0.210335344076156616,  -0.184346809983253479,  -0.153087124228477478,   0.054921828210353851,   0.146172195672988892,  -0.132020518183708191,  -0.118569366633892059,  -0.085994690656661987,   0.082717575132846832,  -0.202392295002937317,   0.197953611612319946,   0.001000450924038887,   0.225327953696250916,   0.083320327103137970,   0.017833765596151352,  -0.050387963652610779,   0.162010401487350464,   0.104472883045673370,  -0.120948672294616699,  -0.189997375011444092,
		   0.131891876459121704,  -0.072596110403537750,  -0.214286595582962036,  -0.079226918518543243,  -0.013667362742125988,  -0.018526514992117882,   0.152137622237205505,   0.226147770881652832,  -0.171163335442543030,   0.125708982348442078,   0.054990120232105255,  -0.048831708729267120,   0.227375522255897522,  -0.031224688515067101,  -0.234086379408836365,  -0.093196123838424683,   0.161037817597389221,  -0.011975143104791641,   0.147086724638938904,  -0.007076272740960121,  -0.161182314157485962,   0.084048561751842499,  -0.003949416335672140,   0.047417137771844864,
		  -0.174016609787940979,   0.171842157840728760,  -0.196859613060951233,  -0.053426038473844528,   0.069617383182048798,  -0.185905560851097107,   0.173014268279075623,   0.021485030651092529,  -0.167390525341033936,  -0.185185879468917847,  -0.235575959086418152,   0.231451973319053650,   0.203157141804695129,  -0.179199099540710449,  -0.100721597671508789,  -0.064948678016662598,   0.030972789973020554,  -0.120899088680744171,  -0.182158842682838440,   0.166038349270820618,  -0.134558707475662231,   0.193049535155296326,  -0.206814408302307129,  -0.130298718810081482,
		  -0.066888287663459778,   0.018188409507274628,   0.108236052095890045,   0.099524453282356262,  -0.011017865501344204,  -0.010902509093284607,   0.097856417298316956,   0.039022967219352722,  -0.053080666810274124,  -0.057871762663125992,  -0.067993588745594025,  -0.059920411556959152,  -0.135927334427833557,   0.143263816833496094,   0.057382527738809586,  -0.027897454798221588,  -0.211796730756759644,   0.031941901892423630,   0.048023808747529984,   0.154774263501167297,  -0.011418030597269535,  -0.111684918403625488,  -0.151802510023117065,   0.168926060199737549,
		  -0.081654123961925507,  -0.121594585478305817,  -0.108362630009651184,  -0.163093328475952148,  -0.235036447644233704,  -0.208788827061653137,   0.242658466100692749,   0.211433082818984985,  -0.118852913379669189,  -0.176829412579536438,   0.084573768079280853,   0.239658132195472717,   0.057825725525617599,   0.050240430980920792,  -0.116622284054756165,  -0.219730392098426819,   0.019183564931154251,   0.049292944371700287,   0.059042990207672119,   0.013476813212037086,  -0.135060533881187439,   0.043805442750453949,   0.157244369387626648,  -0.018672885373234749,
	};
	static const double bias06[]=
	{
		  -0.207236602902412415,  -0.069941431283950806,  -0.119261406362056732,  -0.039139635860919952,   0.121392160654067993,   0.097777545452117920,  -0.089829504489898682,   0.092655986547470093,  -0.151401072740554810,   0.167809292674064636,  -0.031289990991353989,  -0.042157560586929321,
	};
	static const double weight07[]=
	{
		  -0.228611856698989868,   0.277465879917144775,   0.294378429651260376,  -0.250823587179183960,   0.073585435748100281,   0.169272065162658691,  -0.081664375960826874,  -0.078288242220878601,   0.089558362960815430,   0.288079887628555298,   0.123927362263202667,   0.147792458534240723,
	};
	static const double bias07[]=
	{
		  -0.137596815824508667,
	};
#endif
	const double *src_weights[]=
	{
		weight01, bias01,
		weight02, bias02,
		weight03, bias03,
		weight04, bias04,
		weight05, bias05,
		weight06, bias06,
		weight07, bias07,
	};
	size_t src_sizes[]=
	{
		_countof(weight01), _countof(bias01),
		_countof(weight02), _countof(bias02),
		_countof(weight03), _countof(bias03),
		_countof(weight04), _countof(bias04),
		_countof(weight05), _countof(bias05),
		_countof(weight06), _countof(bias06),
		_countof(weight07), _countof(bias07),
	};
	C03_DATATYPE *coeffs[_countof(src_weights)];
	for(int k=0;k<_countof(src_weights);++k)
	{
		coeffs[k]=(C03_DATATYPE*)_mm_malloc(src_sizes[k]*sizeof(C03_DATATYPE), 32);
		for(int k2=0;k2<src_sizes[k];++k2)
			coeffs[k][k2]=C03_CVT_COEFF(src_weights[k][k2]);
	}
	int nnch[]=
	{
		C03_CIN,
		C03_CIN,
		C03_CIN,
		C03_CIN,
		C03_CIN,
		C03_NNB,
		C03_NNB/2,
		1,
	};

	double avpred=0;//

	int idx, idx2, pred;
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	memcpy(dst, src, (size_t)res<<2);
	for(int kc=0;kc<3;++kc)//channel loop
	{
		for(int ky=0;ky<ih;++ky)//y loop
		{
			for(int kx=0;kx<iw;++kx)//x loop
			{
				//if(kc==0&&kx==10&&ky==10)
				//	printf("");

				ALIGN(32) C03_DATATYPE nb[C03_CIN], temp[C03_CIN];
				idx=0;
				for(int ky2=-C03_REACH;ky2<0;++ky2)//load neighbors
				{
					for(int kx2=-C03_REACH;kx2<=C03_REACH;++kx2, ++idx)
					{
						if((unsigned)(ky+ky2)<(unsigned)ih&&(unsigned)(kx+kx2)<(unsigned)iw)
						{
							idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
							nb[idx]=C03_CVT_PIXEL(pixels[idx2]);
							nb[idx+C03_NNB]=C03_CVT_PIXEL(errors[idx2]);
						}
						else
						{
							nb[idx]=0;
							nb[idx+C03_NNB]=0;
						}
					}
				}
				for(int kx2=-C03_REACH;kx2<0;++kx2, ++idx)
				{
					if((unsigned)(kx+kx2)<(unsigned)iw)
					{
						idx2=(iw*ky+kx+kx2)<<2|kc;
						nb[idx]=C03_CVT_PIXEL(pixels[idx2]);
						nb[idx+C03_NNB]=C03_CVT_PIXEL(errors[idx2]);
					}
					else
					{
						nb[idx]=0;
						nb[idx+C03_NNB]=0;
					}
				}
				
				C03_DATATYPE *vi=nb, *vo=temp;
				for(int k=1;k<_countof(nnch);++k)
				{
					int ci=nnch[k-1], co=nnch[k];
					c03_dense(coeffs[(k-1)<<1], vi, coeffs[(k-1)<<1|1], vo, ci, co);
					
					C03_DATATYPE *t2;
					SWAPVAR(vi, vo, t2);
				}
				pred=C03_CVT_PRED(*vi);
				//int pred2=(int)(((*vi)>(-1) ? (*vi)<(1) ? (*vi) : (1) : (-1)) * 128.f + 0.5f);

				avpred+=fabs(*vi*128);//

				idx=(iw*ky+kx)<<2|kc;
				pred^=-fwd;
				pred+=fwd;
				dst[idx]=src[idx]+pred;
			}//x loop
		}//y loop
	}//channel loop
	memcpy(dst, src, (size_t)res<<2);
	free(dst);
	for(int k=0;k<_countof(src_weights);++k)
		_mm_free(coeffs[k]);

	set_window_title("%lf", avpred/(res*3));//
}