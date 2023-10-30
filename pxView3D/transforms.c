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
#define C03_DATATYPE float
#define C03_CVT_COEFF(X) (C03_DATATYPE)(X)
#define C03_CVT_PIXEL(X) (C03_DATATYPE)(((X)+128)/255.*2-1)
#define C03_CVT_PRED(X) (char)((CLAMP(-1, X, 1)+1)/2*255-128)
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
#if 0
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
		//s2=s2<0?s2*0.01:s2;//LeakyReLU
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
static void c03_activation(C03_DATATYPE *vec, int count)
{
	for(int k=0;k<count;++k)
	{
		C03_DATATYPE val=vec[k];
		if(val<0)//LeakyReLU
			val*=0.01;
		vec[k]=val;
	}
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
		  -0.085383899509906769,   0.118390068411827087,  -0.136094152927398682,  -0.151561245322227478,   0.030368980020284653,  -0.109790407121181488,  -0.112390592694282532,   0.017681185156106949,  -0.034899804741144180,   0.125526353716850281,  -0.093470938503742218,   0.056151844561100006,  -0.078147679567337036,  -0.098730660974979401,   0.089048556983470917,  -0.017100425437092781,   0.087147012352943420,  -0.034936301410198212,  -0.116087533533573151,  -0.044706199318170547,  -0.021290479227900505,   0.050434906035661697,   0.086150385439395905,   0.023349676281213760,   0.031740311533212662,  -0.060442496091127396,  -0.104596666991710663,   0.094257295131683350,  -0.072534427046775818,  -0.024137133732438087,   0.049870446324348450,   0.111887060105800629,   0.029084386304020882,  -0.120553560554981232,  -0.072650514543056488,  -0.166243672370910645,  -0.151480972766876221,  -0.123768836259841919,  -0.018935125321149826,  -0.183374628424644470,  -0.202874198555946350,  -0.061107669025659561,   0.051188673824071884,  -0.136303275823593140,  -0.145607233047485352,   0.025814196094870567,   0.050803288817405701,  -0.152526676654815674,
		  -0.071527726948261261,   0.065523825585842133,  -0.059865076094865799,   0.046400289982557297,  -0.133711576461791992,  -0.050455771386623383,   0.118131630122661591,   0.087230928242206573,   0.131875738501548767,   0.043416678905487061,   0.081556148827075958,  -0.094485603272914886,  -0.069874010980129242,   0.055102590471506119,  -0.120514325797557831,   0.075725302100181580,  -0.070086240768432617,  -0.141116380691528320,  -0.055445991456508636,  -0.119016416370868683,   0.038672074675559998,   0.072363339364528656,   0.028751540929079056,  -0.024036537855863571,   0.056025721132755280,  -0.072045236825942993,   0.049701206386089325,  -0.068358190357685089,  -0.107125960290431976,   0.010339006781578064,  -0.075832843780517578,  -0.028757285326719284,   0.073171719908714294,  -0.042389459908008575,  -0.132986009120941162,  -0.140099644660949707,  -0.075206048786640167,  -0.040282137691974640,   0.083268575370311737,   0.054923374205827713,  -0.058886416256427765,  -0.244638919830322266,  -0.116247072815895081,   0.098645165562629700,  -0.074600100517272949,  -0.116877354681491852,  -0.017066322267055511,  -0.084240019321441650,
		  -0.000460180483059958,   0.088042475283145905,  -0.143661141395568848,  -0.049939282238483429,   0.103497743606567383,   0.006299057509750128,  -0.122538551688194275,   0.067048400640487671,   0.078384600579738617,  -0.074274614453315735,  -0.018591793254017830,  -0.101122058928012848,  -0.023696430027484894,   0.108617730438709259,  -0.132320970296859741,  -0.102961286902427673,  -0.091966181993484497,  -0.149186298251152039,  -0.077539250254631042,   0.103994175791740417,   0.009383607655763626,  -0.111287653446197510,   0.088940382003784180,  -0.026342965662479401,  -0.087295055389404297,  -0.034473147243261337,  -0.032307550311088562,  -0.000277756480500102,  -0.030988633632659912,  -0.050812579691410065,  -0.115324303507804871,   0.033828005194664001,  -0.052957274019718170,   0.023988816887140274,  -0.148886874318122864,  -0.145227611064910889,   0.104951038956642151,   0.021466294303536415,  -0.000895963516086340,  -0.029801705852150917,   0.039252858608961105,  -0.042209863662719727,  -0.111893996596336365,   0.057233475148677826,   0.074059784412384033,   0.025607859715819359,  -0.037837732583284378,  -0.137911602854728699,
		  -0.081440590322017670,   0.099579438567161560,   0.133599966764450073,   0.039643697440624237,  -0.075847022235393524,  -0.012557107023894787,  -0.047102395445108414,   0.125250756740570068,  -0.099792629480361938,  -0.070460617542266846,   0.072240456938743591,  -0.009167276322841644,   0.133299350738525391,  -0.028349567204713821,  -0.076107613742351532,   0.107622779905796051,  -0.135000899434089661,  -0.104312233626842499,  -0.168799430131912231,  -0.030804362148046494,  -0.025212852284312248,   0.080616965889930725,  -0.068797633051872253,  -0.087972685694694519,   0.018724100664258003,   0.040828954428434372,  -0.138900443911552429,  -0.007962954230606556,   0.111153475940227509,  -0.014933977276086807,   0.098583668470382690,   0.037372566759586334,   0.024395940825343132,  -0.064573213458061218,   0.057440441101789474,  -0.040233027189970016,  -0.121238864958286285,   0.031350690871477127,  -0.179775997996330261,  -0.156538471579551697,  -0.193493917584419250,  -0.117110751569271088,  -0.175564706325531006,   0.112588763236999512,   0.055798664689064026,   0.123570010066032410,  -0.068234652280807495,  -0.165356084704399109,
		   0.104637801647186279,   0.148196414113044739,  -0.020716156810522079,   0.005288233980536461,  -0.018761646002531052,  -0.102115064859390259,  -0.002105508698150516,  -0.080324724316596985,   0.069677628576755524,  -0.000516234722454101,  -0.009354592300951481,   0.076439246535301208,   0.027084140107035637,   0.060428678989410400,   0.010036380961537361,   0.095874913036823273,   0.158389315009117126,   0.023905467242002487,  -0.076257593929767609,   0.105195678770542145,   0.059953961521387100,   0.109640948474407196,   0.135035559535026550,  -0.014746472239494324,   0.158362075686454773,   0.155717283487319946,   0.121458016335964203,   0.068729199469089508,  -0.091008387506008148,  -0.004664309322834015,   0.131759718060493469,   0.032833769917488098,  -0.026249077171087265,   0.069176152348518372,   0.025511393323540688,   0.008257143199443817,   0.037578482180833817,  -0.005619410425424576,   0.148932516574859619,   0.014767665416002274,  -0.063705205917358398,   0.206272944808006287,   0.158683761954307556,   0.089372806251049042,   0.011828553862869740,  -0.094597652554512024,   0.072318784892559052,   0.179552197456359863,
		  -0.061777010560035706,  -0.091917425394058228,   0.086529947817325592,   0.081289127469062805,   0.088728778064250946,  -0.124017372727394104,  -0.107058539986610413,  -0.038083422929048538,  -0.024296859279274940,   0.148504033684730530,   0.025014996528625488,  -0.085601560771465302,  -0.078179903328418732,  -0.046188946813344955,  -0.093151360750198364,  -0.059956971555948257,  -0.150475591421127319,  -0.050594512373209000,   0.012189546599984169,  -0.076149836182594299,  -0.165922820568084717,  -0.142131581902503967,   0.132767796516418457,  -0.076410338282585144,  -0.032188866287469864,  -0.047796074301004410,  -0.051434446126222610,  -0.111136920750141144,   0.001513922237791121,  -0.091183200478553772,   0.043219875544309616,   0.102557845413684845,  -0.020098837092518806,  -0.072258152067661285,  -0.034818518906831741,  -0.120237261056900024,  -0.023816434666514397,   0.082423493266105652,  -0.039690218865871429,   0.080430075526237488,   0.037499886006116867,  -0.056136220693588257,   0.035570107400417328,  -0.044422082602977753,   0.110225751996040344,  -0.109800100326538086,  -0.065702229738235474,   0.035133920609951019,
		   0.003971314057707787,  -0.124053329229354858,  -0.166314229369163513,  -0.046547651290893555,  -0.015269814990460873,  -0.044760916382074356,   0.018757205456495285,  -0.040277268737554550,  -0.009004036895930767,  -0.026397855952382088,  -0.137230783700942993,  -0.187192529439926147,   0.056376513093709946,  -0.122627556324005127,  -0.123183749616146088,   0.005771602503955364,  -0.019722364842891693,  -0.115626744925975800,  -0.028967505320906639,  -0.179362371563911438,  -0.166920974850654602,   0.084189370274543762,   0.048071920871734619,  -0.091526627540588379,   0.098732702434062958,  -0.081359602510929108,   0.053600013256072998,  -0.084155917167663574,   0.009269944392144680,   0.010053627192974091,   0.046161841601133347,   0.082315020263195038,   0.026454715058207512,  -0.109696030616760254,  -0.063879594206809998,  -0.084170997142791748,  -0.150692358613014221,   0.095110677182674408,  -0.118722885847091675,  -0.098871752619743347,  -0.052572969347238541,   0.098578870296478271,   0.032013483345508575,   0.049071099609136581,  -0.021895291283726692,  -0.181625500321388245,  -0.183359056711196899,   0.017621304839849472,
		  -0.054699894040822983,  -0.031395144760608673,  -0.038961037993431091,   0.056650619953870773,   0.059446215629577637,  -0.138094663619995117,   0.129716321825981140,  -0.093494586646556854,  -0.136633306741714478,  -0.008618881925940514,  -0.017182992771267891,   0.100979007780551910,   0.129212051630020142,  -0.146698102355003357,   0.031104572117328644,   0.025895765051245689,  -0.110763303935527802,  -0.198732122778892517,   0.023918451741337776,   0.007097783964127302,  -0.144679412245750427,  -0.093990772962570190,  -0.094157539308071136,  -0.196323499083518982,   0.044247727841138840,  -0.106850191950798035,  -0.031494863331317902,   0.082145340740680695,   0.086683653295040131,  -0.132577568292617798,   0.039713010191917419,  -0.045589815825223923,   0.037134047597646713,  -0.114438533782958984,  -0.039849348366260529,  -0.153153091669082642,  -0.146019056439399719,   0.045557022094726562,   0.084922023117542267,   0.052754767239093781,   0.050126828253269196,   0.026554897427558899,   0.010965585708618164,   0.067032761871814728,   0.013651339337229729,   0.070731066167354584,  -0.115135632455348969,  -0.204158157110214233,
		   0.094495259225368500,  -0.156629279255867004,   0.009818433783948421,  -0.123450137674808502,  -0.063529193401336670,   0.039534877985715866,  -0.075046434998512268,  -0.011512297205626965,  -0.002577364444732666,   0.001253015827387571,  -0.148577600717544556,  -0.025639623403549194,   0.055597528815269470,   0.059060152620077133,  -0.133626401424407959,   0.015817033126950264,   0.037960804998874664,  -0.045439559966325760,   0.044360183179378510,   0.048758968710899353,  -0.044631525874137878,  -0.149510771036148071,   0.087645187973976135,  -0.111140452325344086,  -0.011749321594834328,  -0.000076564923801925,  -0.047150168567895889,   0.026970919221639633,  -0.133548215031623840,   0.061079785227775574,  -0.052195634692907333,  -0.099733240902423859,   0.032390564680099487,   0.001305662794038653,   0.084385856986045837,  -0.047990031540393829,  -0.016449188813567162,   0.038690187036991119,   0.065887540578842163,   0.012283770367503166,  -0.189933359622955322,  -0.224218890070915222,  -0.026209453120827675,  -0.089597441256046295,   0.030947018414735794,  -0.066027954220771790,  -0.036922544240951538,  -0.140550091862678528,
		   0.057462040334939957,  -0.074758999049663544,   0.120023742318153381,   0.047004625201225281,  -0.075138740241527557,   0.154395803809165955,   0.027145924046635628,   0.140760496258735657,  -0.066787995398044586,   0.049561083316802979,   0.118814341723918915,   0.082001008093357086,   0.032288517802953720,   0.096602953970432281,   0.029245864599943161,  -0.018244057893753052,   0.053044829517602921,   0.060417421162128448,  -0.144363626837730408,   0.010362757369875908,  -0.003107965225353837,  -0.135443657636642456,  -0.110549405217170715,   0.071278728544712067,   0.117598503828048706,   0.114605680108070374,  -0.119898632168769836,   0.104630753397941589,   0.072770416736602783,  -0.080859571695327759,   0.071434430778026581,   0.103117063641548157,  -0.121725335717201233,   0.111563459038734436,  -0.023631276562809944,  -0.020944511517882347,  -0.011944940313696861,  -0.039479851722717285,  -0.083110868930816650,   0.076675146818161011,  -0.128836557269096375,  -0.184868931770324707,  -0.143163174390792847,   0.125471293926239014,  -0.032582931220531464,   0.021134151145815849,   0.093273557722568512,  -0.111196376383304596,
		  -0.145380899310112000,  -0.040387462824583054,  -0.099400579929351807,   0.012294747866690159,  -0.108507566154003143,  -0.029944021254777908,  -0.055760532617568970,  -0.083575643599033356,  -0.083833381533622742,   0.146625295281410217,   0.002742210635915399,  -0.101951554417610168,   0.054710004478693008,  -0.045062921941280365,   0.093035846948623657,  -0.058828905224800110,  -0.033113390207290649,  -0.176600277423858643,   0.011958690360188484,   0.084141820669174194,  -0.108642481267452240,   0.093165323138237000,   0.154055163264274597,  -0.101055383682250977,  -0.038491599261760712,  -0.055925697088241577,  -0.000856936909258366,   0.025388823822140694,   0.015262376517057419,   0.081074163317680359,   0.007300006691366434,   0.155556708574295044,  -0.048447929322719574,   0.164258673787117004,  -0.056159261614084244,   0.105457454919815063,  -0.087048120796680450,  -0.029539849609136581,   0.005721843801438808,  -0.052963994443416595,   0.049772176891565323,  -0.181449338793754578,   0.096073672175407410,  -0.049697693437337875,  -0.011371552944183350,   0.029216563329100609,  -0.116089500486850739,  -0.204270392656326294,
		   0.123481638729572296,   0.003702368587255478,  -0.013035267591476440,  -0.016295766457915306,  -0.031110117211937904,  -0.059825535863637924,   0.025256291031837463,  -0.104023613035678864,   0.036594763398170471,  -0.058286163955926895,  -0.082161404192447662,  -0.022360613569617271,  -0.102263003587722778,   0.105560906231403351,  -0.101202107965946198,   0.020144926384091377,  -0.096624195575714111,  -0.110916256904602051,  -0.126184895634651184,  -0.002041050465777516,  -0.099139727652072906,  -0.035245355218648911,   0.112516358494758606,  -0.068924784660339355,  -0.059468094259500504,   0.078005261719226837,  -0.074139177799224854,   0.030428400263190269,   0.016806803643703461,  -0.108100414276123047,  -0.065617233514785767,  -0.089701674878597260,  -0.044010449200868607,  -0.081475853919982910,   0.015696996822953224,   0.009806212969124317,  -0.012606726959347725,  -0.085679784417152405,   0.094896763563156128,  -0.060334939509630203,  -0.172253847122192383,   0.034038424491882324,  -0.167589157819747925,   0.061056818813085556,   0.124520674347877502,  -0.107221886515617371,   0.010584514588117599,  -0.158202469348907471,
		  -0.046123575419187546,  -0.119731500744819641,  -0.144239038228988647,   0.027118692174553871,  -0.148928269743919373,  -0.081912070512771606,   0.131712809205055237,  -0.140366449952125549,  -0.069119125604629517,  -0.047299496829509735,   0.004944362677633762,  -0.045831110328435898,  -0.032627008855342865,   0.065528221428394318,  -0.039507649838924408,  -0.127604618668556213,  -0.144811064004898071,  -0.017395449802279472,  -0.092665500938892365,  -0.111924834549427032,   0.055034592747688293,   0.088399760425090790,   0.034282192587852478,  -0.193687990307807922,   0.038271289318799973,  -0.076638065278530121,   0.036467719823122025,   0.053995989263057709,   0.052569720894098282,  -0.013348110951483250,   0.095646470785140991,  -0.118549622595310211,  -0.120414309203624725,  -0.130507960915565491,  -0.047617513686418533,  -0.031848024576902390,  -0.142555087804794312,  -0.015452180057764053,   0.030988095328211784,  -0.150718480348587036,  -0.124587737023830414,  -0.238715916872024536,  -0.150257751345634460,  -0.059581670910120010,  -0.043321110308170319,  -0.086259998381137848,   0.073475956916809082,  -0.203551769256591797,
		   0.126236617565155029,   0.017348209396004677,  -0.077523328363895416,  -0.048946883529424667,  -0.042074777185916901,  -0.082561962306499481,   0.038146499544382095,  -0.035845011472702026,  -0.023664193227887154,   0.102402172982692719,  -0.145481795072555542,  -0.037540432065725327,  -0.069719299674034119,  -0.077300295233726501,   0.099983312189579010,  -0.134226694703102112,   0.103832565248012543,   0.181414261460304260,  -0.108275778591632843,  -0.015857079997658730,  -0.034719537943601608,  -0.038670122623443604,   0.016417708247900009,   0.073741160333156586,  -0.029645137488842010,  -0.123104028403759003,  -0.060167830437421799,   0.013591964729130268,   0.011823143810033798,  -0.042650785297155380,  -0.014328787103295326,   0.054187867790460587,   0.079583883285522461,  -0.162460103631019592,   0.044654928147792816,   0.128098636865615845,   0.039586134254932404,  -0.075309194624423981,  -0.085678227245807648,   0.071695312857627869,  -0.060256410390138626,   0.209568366408348083,  -0.089812248945236206,   0.039207845926284790,  -0.085574142634868622,   0.027354927733540535,   0.002839191118255258,   0.132661610841751099,
		  -0.065481267869472504,   0.104507453739643097,   0.036746755242347717,  -0.110377676784992218,   0.032013192772865295,   0.005341478623449802,   0.016887286677956581,   0.144994974136352539,   0.046478953212499619,  -0.106404304504394531,   0.136638581752777100,   0.099010616540908813,  -0.003494796343147755,  -0.090608879923820496,   0.059943016618490219,  -0.119636774063110352,   0.072003997862339020,  -0.218301624059677124,   0.061041042208671570,  -0.119245782494544983,   0.028556486591696739,   0.001623120740987360,  -0.146166488528251648,   0.048047464340925217,  -0.014046982862055302,   0.115134656429290771,   0.090220607817173004,  -0.040973309427499771,  -0.096829794347286224,   0.013741422444581985,   0.030967814847826958,  -0.041157934814691544,  -0.054501492530107498,   0.128602996468544006,   0.120351754128932953,   0.108165301382541656,  -0.114939890801906586,  -0.103484518826007843,   0.052882049232721329,  -0.136144459247589111,   0.024641353636980057,  -0.087184548377990723,   0.084630720317363739,  -0.080296166241168976,   0.119390174746513367,  -0.138373956084251404,   0.006401062943041325,  -0.127474218606948853,
		   0.159982994198799133,  -0.039654430001974106,   0.037776324898004532,  -0.104602642357349396,   0.029255596920847893,   0.013968098908662796,  -0.070992887020111084,   0.119602181017398834,   0.156138256192207336,  -0.073183871805667877,  -0.094905838370323181,   0.078841298818588257,  -0.021764099597930908,   0.120298117399215698,  -0.101098299026489258,  -0.006418752484023571,   0.187447547912597656,   0.100403569638729095,   0.208419129252433777,   0.089893817901611328,   0.103678047657012939,   0.155674487352371216,  -0.082810722291469574,   0.166682451963424683,  -0.018329033628106117,  -0.063173331320285797,   0.159977540373802185,   0.006093886680901051,  -0.001285710372030735,  -0.010763422586023808,  -0.110453248023986816,  -0.062893696129322052,  -0.019687671214342117,   0.159299463033676147,   0.128956094384193420,   0.184320077300071716,   0.177857577800750732,   0.167636960744857788,   0.177569925785064697,  -0.039107546210289001,  -0.027859702706336975,   0.005532428622245789,   0.026718821376562119,   0.108235456049442291,   0.043186735361814499,  -0.046650499105453491,  -0.056692879647016525,   0.222587287425994873,
		   0.027199441567063332,  -0.026372510939836502,  -0.023735387250781059,   0.050056330859661102,  -0.012562972493469715,   0.127638772130012512,   0.084896363317966461,  -0.055114336311817169,   0.150366902351379395,   0.027789391577243805,   0.032858904451131821,  -0.062489736825227737,   0.109780192375183105,   0.136005923151969910,   0.002706402679905295,  -0.071784257888793945,   0.120296835899353027,   0.022880714386701584,   0.180990248918533325,  -0.034149281680583954,  -0.101041629910469055,  -0.074721723794937134,   0.155351877212524414,   0.153235226869583130,   0.118899308145046234,   0.029837753623723984,  -0.052449617534875870,   0.093112766742706299,   0.061377678066492081,  -0.077707283198833466,  -0.011463488452136517,   0.018250141292810440,   0.082404136657714844,  -0.010343374684453011,   0.165098622441291809,   0.061647802591323853,  -0.082319542765617371,   0.091013409197330475,  -0.066924795508384705,  -0.025907546281814575,   0.083383612334728241,   0.236092478036880493,   0.036129117012023926,  -0.013761483132839203,   0.069094970822334290,  -0.105984285473823547,   0.070724204182624817,   0.198194921016693115,
		  -0.018498105928301811,   0.129379451274871826,   0.021107306703925133,  -0.043593119829893112,   0.131307899951934814,  -0.113841816782951355,   0.112168096005916595,   0.005409646779298782,   0.071227356791496277,   0.041113417595624924,   0.081673026084899902,   0.028220914304256439,   0.102591194212436676,   0.079025611281394958,   0.121862784028053284,   0.068214900791645050,   0.085954181849956512,  -0.168265059590339661,  -0.048064507544040680,   0.001599310780875385,   0.023506186902523041,   0.037558834999799728,  -0.093195334076881409,  -0.191561773419380188,   0.002082026097923517,   0.102932274341583252,  -0.104209043085575104,  -0.158244669437408447,  -0.072606131434440613,   0.071390211582183838,   0.041253633797168732,   0.125942140817642212,   0.069925032556056976,  -0.104536794126033783,   0.031786728650331497,  -0.061721410602331161,  -0.115580759942531586,   0.011220466345548630,   0.017311727628111839,  -0.035453800112009048,   0.023115202784538269,   0.032537333667278290,  -0.046072036027908325,  -0.152752652764320374,   0.089386790990829468,  -0.080860413610935211,  -0.142630249261856079,  -0.001038932241499424,
		   0.001406693481840193,   0.021644307300448418,  -0.067309215664863586,  -0.089147344231605530,  -0.010290267877280712,   0.044577438384294510,   0.139594435691833496,   0.102902159094810486,   0.135523214936256409,  -0.055151604115962982,  -0.114234745502471924,   0.051285643130540848,   0.053325068205595016,   0.081919059157371521,  -0.067791841924190521,  -0.086002275347709656,   0.068579904735088348,   0.197592213749885559,   0.014988054521381855,   0.028286151587963104,   0.009554819203913212,   0.113063991069793701,   0.064504839479923248,   0.171415165066719055,   0.037828680127859116,   0.033138811588287354,   0.137212485074996948,   0.059932209551334381,   0.144529536366462708,  -0.048664174973964691,  -0.043243996798992157,   0.136671736836433411,  -0.040452901273965836,   0.115775190293788910,   0.038716875016689301,  -0.027933306992053986,   0.002070277929306030,  -0.019039034843444824,   0.042402692139148712,  -0.060574151575565338,   0.065344862639904022,   0.225320950150489807,  -0.062362343072891235,   0.122671373188495636,  -0.097947821021080017,  -0.057901177555322647,   0.011568034999072552,   0.195627838373184204,
		  -0.033411994576454163,  -0.050995856523513794,   0.078464232385158539,  -0.101810060441493988,  -0.095142498612403870,   0.033464491367340088,  -0.027710564434528351,  -0.128978833556175232,  -0.094510987401008606,   0.026043772697448730,  -0.113737747073173523,  -0.005044265184551477,  -0.076937817037105560,   0.109447047114372253,  -0.003281770972535014,   0.126161411404609680,   0.101085603237152100,   0.073531024158000946,   0.001679701381362975,   0.010792400687932968,  -0.158280596137046814,  -0.134971246123313904,  -0.111192010343074799,  -0.171177044510841370,  -0.031566258519887924,  -0.022283829748630524,   0.086058542132377625,  -0.138188198208808899,   0.081218138337135315,   0.094262219965457916,   0.101349346339702606,  -0.112574085593223572,   0.039721686393022537,  -0.000757239817176014,  -0.056353881955146790,  -0.075262844562530518,  -0.121607825160026550,  -0.033587176352739334,   0.047018647193908691,  -0.011816503480076790,  -0.052647724747657776,  -0.049485690891742706,  -0.179851919412612915,  -0.130277127027511597,   0.035924922674894333,  -0.046114366501569748,  -0.130109012126922607,  -0.030048668384552002,
		  -0.063065424561500549,   0.076897956430912018,   0.125922232866287231,   0.142887428402900696,   0.055454846471548080,   0.064359463751316071,  -0.036578405648469925,   0.141475960612297058,   0.161684304475784302,   0.124584242701530457,   0.089657247066497803,  -0.120367579162120819,  -0.125542610883712769,   0.048438183963298798,  -0.049593869596719742,   0.002998002804815769,  -0.015433161519467831,   0.219317302107810974,  -0.070575132966041565,   0.002118948614224792,   0.093703866004943848,   0.124790534377098083,   0.007086551282554865,   0.064918108284473419,   0.041714690625667572,  -0.028909888118505478,   0.139935374259948730,  -0.002281527500599623,  -0.022257043048739433,   0.026211580261588097,   0.135250017046928406,   0.154009833931922913,   0.081442780792713165,   0.078069850802421570,   0.107815735042095184,   0.058262962847948074,   0.032703276723623276,   0.082122936844825745,  -0.010565158911049366,  -0.051188942044973373,   0.013905329629778862,   0.055260781198740005,   0.030376406386494637,   0.157513961195945740,   0.125146135687828064,   0.116062432527542114,   0.033365599811077118,   0.167080521583557129,
		   0.056163582950830460,   0.056018512696027756,   0.145543351769447327,  -0.133256763219833374,   0.090808413922786713,   0.153475642204284668,  -0.012435927055776119,   0.021025510504841805,  -0.030122011899948120,  -0.068296626210212708,   0.077424161136150360,  -0.010918001644313335,   0.127322211861610413,  -0.015416597016155720,  -0.087915852665901184,   0.070178665220737457,  -0.061756163835525513,   0.078290767967700958,   0.028862152248620987,  -0.060816079378128052,   0.063631422817707062,  -0.134810432791709900,  -0.043648779392242432,  -0.220649048686027527,   0.028519671410322189,   0.068579368293285370,   0.103371821343898773,   0.012174431234598160,  -0.012954905629158020,  -0.033368863165378571,  -0.058564674109220505,   0.064771808683872223,  -0.069338053464889526,   0.133591085672378540,   0.006234460044652224,  -0.088472209870815277,   0.077863201498985291,  -0.031599760055541992,   0.073074541985988617,  -0.082751870155334473,  -0.018866140395402908,  -0.165093824267387390,   0.085599787533283234,  -0.093317084014415741,   0.007419424597173929,   0.076638579368591309,  -0.065989531576633453,  -0.123912401497364044,
		  -0.148807898163795471,  -0.011277826502919197,   0.109036818146705627,   0.103230379521846771,  -0.072483934462070465,  -0.056470274925231934,  -0.063802517950534821,  -0.059424672275781631,  -0.040287118405103683,  -0.020706417039036751,   0.138982340693473816,  -0.085304334759712219,   0.083156898617744446,   0.111107885837554932,  -0.032169979065656662,   0.049033846706151962,  -0.132682770490646362,  -0.171106413006782532,  -0.170952349901199341,   0.067389555275440216,   0.043692439794540405,  -0.068939834833145142,  -0.012141202576458454,  -0.100393526256084442,  -0.087099574506282806,  -0.113696932792663574,   0.052084214985370636,  -0.017417646944522858,   0.117865793406963348,  -0.029850432649254799,   0.068982705473899841,  -0.106855742633342743,  -0.126065790653228760,  -0.123039938509464264,  -0.125350385904312134,   0.078316114842891693,   0.093413032591342926,  -0.159815683960914612,  -0.004678430967032909,  -0.043737057596445084,  -0.124070033431053162,  -0.100157193839550018,   0.048796866089105606,   0.096512801945209503,   0.051699887961149216,  -0.119115293025970459,   0.047677971422672272,  -0.200944691896438599,
		   0.148881390690803528,   0.103880181908607483,   0.133407965302467346,  -0.010026194155216217,   0.092970944941043854,  -0.060511283576488495,  -0.042560771107673645,   0.058506283909082413,   0.149892508983612061,  -0.003199306782335043,   0.171158060431480408,  -0.004030040465295315,   0.157829001545906067,   0.056859798729419708,   0.118631958961486816,  -0.032955288887023926,   0.069648489356040955,  -0.066150106489658356,   0.192105621099472046,   0.024379475042223930,   0.106751807034015656,   0.022233095020055771,   0.077804848551750183,   0.074176609516143799,   0.080595210194587708,   0.065850771963596344,  -0.035503018647432327,  -0.001116397907026112,  -0.052830465137958527,   0.102618366479873657,  -0.001213572104461491,  -0.032175339758396149,  -0.074263468384742737,   0.098962359130382538,   0.041706595569849014,  -0.014146627858281136,   0.123785108327865601,   0.038207929581403732,   0.189530953764915466,   0.063820287585258484,   0.194566741585731506,  -0.014060521498322487,   0.126179888844490051,   0.059373497962951660,  -0.008687364868819714,   0.086488798260688782,   0.022113243117928505,   0.037351865321397781,
		  -0.074359275400638580,   0.086445011198520660,  -0.094294048845767975,  -0.077706202864646912,   0.067088946700096130,  -0.029161194339394569,   0.021865114569664001,  -0.084570378065109253,   0.026192162185907364,  -0.036861401051282883,  -0.052937597036361694,   0.101707264780998230,   0.100361622869968414,  -0.106582976877689362,  -0.123860061168670654,  -0.026746163144707680,  -0.062404088675975800,  -0.066030986607074738,  -0.039498899132013321,  -0.143150806427001953,  -0.137559056282043457,  -0.124760366976261139,   0.058183420449495316,  -0.175130888819694519,  -0.144152685999870300,  -0.092047654092311859,  -0.024642927572131157,   0.050670668482780457,   0.053584620356559753,   0.024991365149617195,  -0.092182666063308716,   0.066823482513427734,   0.013952892273664474,   0.093287065625190735,  -0.040205292403697968,  -0.182285532355308533,   0.052824921905994415,  -0.058581303805112839,  -0.004833093378692865,   0.051007140427827835,  -0.123070113360881805,  -0.196633815765380859,  -0.082564637064933777,  -0.093668147921562195,  -0.064612530171871185,   0.053760040551424026,  -0.148463949561119080,  -0.222114726901054382,
		  -0.025720268487930298,  -0.080864243209362030,   0.035264756530523300,  -0.069560341536998749,  -0.112952530384063721,   0.105178028345108032,  -0.071631416678428650,   0.100425325334072113,   0.092481099069118500,  -0.080381229519844055,  -0.082225389778614044,  -0.141872212290763855,   0.047559868544340134,  -0.075846828520298004,   0.007983224466443062,  -0.104306019842624664,   0.056552570313215256,  -0.122517012059688568,   0.024031223729252815,   0.097792118787765503,   0.067906431853771210,  -0.051261506974697113,  -0.084711335599422455,  -0.121343836188316345,   0.003699012333527207,   0.067866705358028412,   0.055240381509065628,   0.003384921466931701,   0.043473914265632629,   0.100596733391284943,   0.056438382714986801,   0.007329539395868778,   0.119021713733673096,  -0.001868864521384239,   0.038485266268253326,   0.031822022050619125,  -0.031170064583420753,   0.007045614533126354,   0.023974740877747536,  -0.074239045381546021,  -0.090874969959259033,  -0.156214594841003418,   0.002753272186964750,  -0.017505597323179245,  -0.149852097034454346,  -0.081084080040454865,  -0.135362967848777771,  -0.095162838697433472,
		   0.015126985497772694,   0.142451912164688110,  -0.066881231963634491,   0.046810522675514221,   0.010044508613646030,   0.038766369223594666,   0.022997304797172546,   0.082137085497379303,  -0.016551822423934937,   0.031470675021409988,   0.046155042946338654,   0.048432093113660812,  -0.011860406957566738,   0.029550692066550255,  -0.105984315276145935,   0.031666669994592667,  -0.056238945573568344,   0.212706640362739563,   0.028911633417010307,   0.159871026873588562,  -0.073138184845447540,   0.166613876819610596,   0.021343052387237549,   0.069831863045692444,   0.023372603580355644,   0.073445722460746765,   0.137832671403884888,   0.130798399448394775,  -0.065902590751647949,   0.139296799898147583,   0.036470320075750351,   0.108698114752769470,   0.041501596570014954,   0.044430598616600037,   0.048141874372959137,   0.068347305059432983,   0.009047470055520535,   0.001130225020460784,   0.134019136428833008,   0.113475091755390167,   0.017592938616871834,   0.027246955782175064,  -0.013159720227122307,   0.073106676340103149,   0.072250030934810638,   0.092054314911365509,  -0.106946982443332672,   0.235244184732437134,
		  -0.039489399641752243,   0.109030768275260925,   0.170187667012214661,   0.127079114317893982,   0.064323768019676208,   0.021541897207498550,   0.144965782761573792,  -0.066386438906192780,  -0.086508013308048248,   0.081390656530857086,   0.164445340633392334,  -0.081568509340286255,   0.068769492208957672,  -0.057179942727088928,  -0.054709151387214661,  -0.007221158593893051,   0.103233039379119873,   0.148557826876640320,   0.020797090604901314,   0.154264092445373535,  -0.086937300860881805,  -0.072785109281539917,  -0.065286360681056976,   0.217318773269653320,  -0.094300761818885803,   0.163857042789459229,   0.065486304461956024,   0.127480834722518921,  -0.085542201995849609,   0.104067027568817139,   0.072647221386432648,  -0.105287194252014160,   0.117911823093891144,   0.024320831522345543,   0.105448134243488312,   0.133816868066787720,   0.038801271468400955,  -0.030194446444511414,   0.148214399814605713,  -0.065241105854511261,  -0.001432084594853222,  -0.001553925685584545,   0.097528532147407532,   0.085706859827041626,  -0.020733878016471863,   0.066475845873355865,   0.143035903573036194,   0.041625019162893295,
		  -0.000589703966397792,  -0.084779992699623108,  -0.113460861146450043,  -0.135552734136581421,  -0.071622632443904877,  -0.014191093854606152,  -0.134865731000900269,  -0.188675373792648315,  -0.079492308199405670,  -0.011978305876255035,  -0.099608436226844788,   0.074688851833343506,  -0.174382388591766357,  -0.184387430548667908,   0.007568473462015390,  -0.165643885731697083,  -0.159323990345001221,   0.118388622999191284,  -0.011558629572391510,   0.073865644633769989,  -0.123468682169914246,   0.044690720736980438,  -0.078344389796257019,   0.110737226903438568,   0.000006799600669183,   0.038965441286563873,   0.103511802852153778,  -0.058125738054513931,   0.042326189577579498,  -0.101766422390937805,   0.084453165531158447,  -0.132174625992774963,  -0.065959349274635315,  -0.011344745755195618,  -0.051930278539657593,  -0.057372737675905228,   0.000822855101432651,   0.056954689323902130,  -0.031227758154273033,  -0.110493503510951996,   0.101624891161918640,   0.053863070905208588,   0.077032156288623810,  -0.143612816929817200,   0.053013511002063751,   0.093842014670372009,   0.021593688055872917,  -0.072897836565971375,
		   0.131735444068908691,   0.105589337646961212,  -0.017905896529555321,  -0.140851691365242004,   0.038528889417648315,  -0.065885990858078003,   0.097455851733684540,  -0.039019342511892319,  -0.152242496609687805,  -0.026749715209007263,  -0.150332197546958923,  -0.077456548810005188,   0.038513965904712677,  -0.146350786089897156,  -0.124035716056823730,   0.054961565881967545,  -0.173261746764183044,   0.007709068711847067,   0.011979185044765472,  -0.142891138792037964,  -0.023919777944684029,   0.015106420964002609,   0.006042185705155134,   0.007130179088562727,   0.035658277571201324,  -0.025358345359563828,   0.023725498467683792,  -0.194580778479576111,   0.029841845855116844,   0.040845829993486404,  -0.067711532115936279,   0.007721786387264729,  -0.149652883410453796,   0.070953987538814545,  -0.122897095978260040,  -0.135258823633193970,  -0.070056580007076263,  -0.037602774798870087,   0.015555912628769875,  -0.067117042839527130,  -0.196101501584053040,  -0.035698324441909790,  -0.021376702934503555,  -0.165358528494834900,  -0.108296990394592285,   0.010923988185822964,  -0.126114189624786377,  -0.214565351605415344,
		  -0.062037725001573563,   0.110707819461822510,  -0.024127919226884842,  -0.013349287211894989,  -0.076568670570850372,  -0.138785764575004578,   0.040041435509920120,   0.080137379467487335,  -0.057429961860179901,  -0.149432301521301270,  -0.041398216038942337,  -0.140231102705001831,   0.063308596611022949,   0.084765821695327759,  -0.123975545167922974,  -0.059358716011047363,   0.059455808252096176,   0.080573633313179016,   0.053965181112289429,  -0.117127604782581329,   0.026953851804137230,   0.055075876414775848,  -0.128504484891891479,   0.013973372057080269,   0.108464904129505157,  -0.007577444426715374,  -0.100271128118038177,   0.059810172766447067,   0.040839329361915588,   0.073687948286533356,  -0.097461447119712830,   0.046743232756853104,  -0.012041872367262840,  -0.032028499990701675,   0.064615756273269653,  -0.093097306787967682,   0.132350102066993713,  -0.089369483292102814,   0.009877254255115986,  -0.022342443466186523,  -0.173663452267646790,  -0.119409620761871338,  -0.009764960035681725,  -0.028529895469546318,  -0.133462980389595032,  -0.131263151764869690,  -0.099587596952915192,   0.034563891589641571,
		   0.087405093014240265,   0.092133603990077972,  -0.111609213054180145,  -0.164814919233322144,   0.088473781943321228,  -0.130188286304473877,  -0.048548936843872070,  -0.125765308737754822,  -0.117986641824245453,   0.000176972753251903,  -0.140783473849296570,  -0.036957290023565292,  -0.115204080939292908,  -0.020084980875253677,  -0.058157756924629211,   0.126726284623146057,  -0.151580378413200378,  -0.212010905146598816,  -0.134643718600273132,  -0.165812686085700989,   0.052695710211992264,   0.059055678546428680,   0.028178596869111061,  -0.199139699339866638,  -0.043204881250858307,  -0.145254284143447876,  -0.130483269691467285,  -0.144737735390663147,  -0.110362626612186432,  -0.110894337296485901,  -0.069098204374313354,  -0.082017563283443451,   0.003040237817913294,   0.025658901780843735,  -0.087423197925090790,  -0.119721740484237671,  -0.129657089710235596,  -0.079836785793304443,  -0.118276432156562805,   0.108220681548118591,  -0.041318655014038086,  -0.053434532135725021,  -0.035951130092144012,   0.029879258945584297,  -0.021680727601051331,  -0.013814794830977917,  -0.149678170680999756,  -0.170052826404571533,
		   0.096573263406753540,  -0.164411619305610657,   0.022207848727703094,  -0.164972215890884399,  -0.044025279581546783,  -0.139501556754112244,   0.124339640140533447,   0.024233873933553696,   0.119651086628437042,   0.149972990155220032,  -0.103764541447162628,  -0.002723023528233171,  -0.105943590402603149,   0.000468124606413767,  -0.041658546775579453,  -0.041096586734056473,  -0.017557449638843536,  -0.091317497193813324,  -0.008445464074611664,   0.096865855157375336,  -0.164589911699295044,  -0.069156371057033539,   0.046694055199623108,  -0.098581083118915558,   0.011780865490436554,   0.122432254254817963,   0.087192662060260773,  -0.035533368587493896,   0.030580295249819756,   0.115650177001953125,  -0.023499285802245140,  -0.028978485614061356,  -0.049588296562433243,   0.115576691925525665,   0.003283326979726553,   0.074863806366920471,  -0.137070089578628540,  -0.156865000724792480,   0.013065312989056110,  -0.095383197069168091,  -0.131064504384994507,  -0.072588048875331879,   0.100596599280834198,   0.095862545073032379,  -0.074661709368228912,  -0.020291393622756004,  -0.035056710243225098,   0.005202531814575195,
		  -0.051045194268226624,  -0.087639756500720978,  -0.037517011165618896,  -0.029141353443264961,  -0.120756164193153381,  -0.073390811681747437,  -0.120943769812583923,  -0.155922889709472656,  -0.154601931571960449,  -0.055687349289655685,  -0.146108195185661316,   0.091769702732563019,  -0.150148183107376099,   0.015378664247691631,   0.024526173248887062,  -0.016850439831614494,  -0.043749686330556870,   0.052198953926563263,  -0.110363714396953583,   0.031536746770143509,   0.094064444303512573,  -0.085365988314151764,  -0.057432260364294052,   0.031282879412174225,  -0.143697604537010193,   0.110571995377540588,  -0.076085142791271210,   0.047124568372964859,  -0.143035694956779480,  -0.013723351992666721,  -0.094217337667942047,  -0.131230458617210388,   0.010959403589367867,   0.040800232440233231,   0.097250260412693024,   0.074253439903259277,   0.009252917952835560,   0.017232073470950127,  -0.079715110361576080,   0.039410870522260666,  -0.101773604750633240,  -0.011963455006480217,  -0.029520086944103241,  -0.132230475544929504,  -0.019328454509377480,  -0.080524735152721405,  -0.077110923826694489,   0.019696922972798347,
		   0.042687393724918365,   0.030711427330970764,   0.096324987709522247,   0.039270028471946716,   0.053783234208822250,   0.043157238513231277,  -0.083800211548805237,   0.031131593510508537,  -0.028684815391898155,  -0.040875751525163651,  -0.115702375769615173,   0.080452047288417816,   0.008579503744840622,  -0.095422700047492981,  -0.014595778658986092,  -0.095873922109603882,  -0.020384095609188080,   0.020736392587423325,   0.084865659475326538,   0.082822032272815704,   0.142634838819503784,   0.089088506996631622,   0.053223758935928345,  -0.046415898948907852,   0.014568937942385674,   0.100182913243770599,   0.007696365471929312,  -0.000061638384067919,   0.052349306643009186,  -0.031300701200962067,  -0.021985094994306564,  -0.031850315630435944,   0.096748717129230499,  -0.082597382366657257,  -0.080755971372127533,  -0.084585711359977722,  -0.098626255989074707,   0.100687578320503235,   0.092828884720802307,  -0.140371009707450867,  -0.190543189644813538,  -0.216401636600494385,  -0.174500107765197754,   0.129098385572433472,  -0.118342526257038116,  -0.085055187344551086,  -0.028977170586585999,  -0.230354741215705872,
		   0.033622063696384430,  -0.043065577745437622,  -0.032428529113531113,   0.046737860888242722,   0.023506931960582733,  -0.088688030838966370,  -0.009298408403992653,  -0.052577391266822815,  -0.011587084271013737,   0.120526224374771118,  -0.039239630103111267,  -0.103738576173782349,   0.106800466775894165,  -0.044626604765653610,   0.117852956056594849,   0.145697012543678284,  -0.063388168811798096,  -0.082744039595127106,   0.047879364341497421,   0.029693260788917542,   0.047336049377918243,   0.012357225641608238,   0.023730117827653885,  -0.058230929076671600,   0.044589783996343613,  -0.098658181726932526,   0.102160573005676270,   0.055225960910320282,   0.093320935964584351,   0.033489171415567398,  -0.037228155881166458,  -0.036462653428316116,   0.126652881503105164,  -0.040801230818033218,  -0.035490535199642181,  -0.027976324781775475,  -0.124198220670223236,   0.036834940314292908,  -0.126236811280250549,   0.051388420164585114,   0.041657663881778717,  -0.088704645633697510,  -0.127076089382171631,  -0.139252334833145142,  -0.013100622221827507,   0.136158242821693420,  -0.123456895351409912,  -0.112431578338146210,
		   0.030390594154596329,  -0.095777198672294617,   0.040643636137247086,  -0.015693100169301033,   0.096131794154644012,  -0.000066875953052659,   0.007644267287105322,   0.048779223114252090,   0.019295243546366692,   0.088617712259292603,   0.009376393631100655,   0.008606882765889168,  -0.062419231981039047,   0.129567205905914307,   0.134540960192680359,   0.163106605410575867,   0.122678935527801514,   0.222318336367607117,   0.139113172888755798,   0.048052012920379639,   0.044059216976165771,  -0.057117760181427002,   0.006804092321544886,   0.173806056380271912,   0.081409148871898651,   0.088503278791904449,  -0.012902216985821724,  -0.038638573139905930,  -0.088986568152904510,   0.036538742482662201,   0.033576149493455887,   0.017322028055787086,   0.073113940656185150,   0.065085627138614655,   0.146686851978302002,  -0.064890004694461823,   0.032589010894298553,  -0.031113427132368088,   0.034855019301176071,   0.139965787529945374,   0.147163093090057373,   0.062804080545902252,   0.028837665915489197,   0.122060902416706085,   0.042082738131284714,  -0.005482014268636703,   0.148437038064002991,   0.113536186516284943,
		   0.073621332645416260,  -0.127980679273605347,  -0.026584761217236519,  -0.076330490410327911,  -0.129158422350883484,   0.046237703412771225,   0.121108956634998322,   0.122397705912590027,  -0.109163232147693634,  -0.114652551710605621,  -0.160853758454322815,  -0.043331895023584366,   0.060629464685916901,  -0.026804722845554352,   0.036659657955169678,  -0.104238636791706085,  -0.092638961970806122,  -0.039555765688419342,   0.040489099919795990,   0.081985339522361755,  -0.004275815095752478,  -0.043551150709390640,  -0.150317892432212830,   0.080026954412460327,  -0.101714700460433960,   0.119875289499759674,   0.033591225743293762,  -0.078979931771755219,  -0.148390352725982666,   0.138796344399452209,  -0.036080352962017059,  -0.089313112199306488,  -0.151187360286712646,  -0.104779161512851715,  -0.135690867900848389,  -0.172016814351081848,   0.063152089715003967,  -0.057506330311298370,  -0.076305560767650604,   0.008271171711385250,  -0.065165728330612183,   0.035847358405590057,  -0.006376751232892275,   0.002610589843243361,  -0.039145715534687042,  -0.006366608198732138,  -0.040701050311326981,  -0.000937500561121851,
		   0.052194759249687195,   0.159192755818367004,  -0.061447516083717346,  -0.036587953567504883,  -0.039514970034360886,  -0.014403978362679482,   0.088719196617603302,   0.030712798237800598,   0.067898176610469818,  -0.067678950726985931,   0.096141256392002106,   0.067671865224838257,   0.028617082163691521,   0.097887478768825531,  -0.059889856725931168,  -0.120611406862735748,   0.160904735326766968,   0.113701149821281433,   0.111562252044677734,  -0.096492841839790344,  -0.104930959641933441,   0.169079706072807312,  -0.044972799718379974,   0.178824260830879211,   0.030400065705180168,  -0.037886511534452438,   0.142578557133674622,  -0.086253114044666290,   0.037307124584913254,  -0.100666649639606476,   0.117065563797950745,   0.115150533616542816,   0.132082998752593994,  -0.038877360522747040,  -0.065815389156341553,   0.042313285171985626,   0.010544959455728531,  -0.065847121179103851,   0.009079040028154850,   0.061674389988183975,   0.181839078664779663,  -0.005560427904129028,   0.184549495577812195,   0.010073952376842499,   0.106056474149227142,   0.015092768706381321,   0.141270324587821960,   0.105699539184570312,
		   0.128898099064826965,   0.124238327145576477,   0.050462163984775543,  -0.030259065330028534,  -0.110411860048770905,   0.123552225530147552,  -0.104246646165847778,  -0.036104541271924973,  -0.118612110614776611,   0.002345963846892118,   0.017022864893078804,  -0.141929522156715393,  -0.020373057574033737,  -0.082000032067298889,   0.029267329722642899,  -0.120459564030170441,  -0.073053449392318726,  -0.065243035554885864,  -0.140639901161193848,  -0.084501989185810089,   0.070236966013908386,  -0.110854208469390869,  -0.145543664693832397,  -0.000736591638997197,  -0.066593483090400696,  -0.037267044186592102,   0.005167576018720865,  -0.032733879983425140,   0.096495173871517181,  -0.167109951376914978,   0.046177793294191360,  -0.117238432168960571,   0.002961903112009168,  -0.029408218339085579,  -0.119717128574848175,   0.006551254540681839,   0.083507694303989410,  -0.113434068858623505,  -0.097433097660541534,  -0.132867112755775452,   0.006123237777501345,  -0.124464377760887146,  -0.010065615177154541,   0.090073868632316589,   0.031142009422183037,  -0.066588543355464935,   0.016314439475536346,  -0.039011213928461075,
		   0.120293177664279938,  -0.002141603734344244,  -0.072570346295833588,   0.048356480896472931,  -0.078126892447471619,   0.016701985150575638,  -0.035704832524061203,   0.006933860015124083,   0.115474447607994080,  -0.102332174777984619,  -0.123543389141559601,  -0.115159049630165100,  -0.050539046525955200,   0.114713788032531738,  -0.025206796824932098,   0.076494187116622925,   0.030388297513127327,  -0.175754338502883911,  -0.045652098953723907,   0.026883758604526520,   0.087543271481990814,   0.034141108393669128,   0.034351691603660583,  -0.190530136227607727,   0.047381889075040817,  -0.048913586884737015,  -0.070324838161468506,  -0.152916848659515381,   0.024577591568231583,  -0.128713458776473999,  -0.034567978233098984,   0.102154888212680817,  -0.006030079443007708,  -0.010874248109757900,   0.027125759050250053,   0.059681676328182220,  -0.017526458948850632,  -0.088303007185459137,  -0.055114276707172394,   0.058118354529142380,   0.058812458068132401,  -0.002357033314183354,  -0.029667522758245468,   0.005915432237088680,  -0.123727075755596161,  -0.057681750506162643,  -0.081371136009693146,  -0.230896666646003723,
		  -0.010940249077975750,   0.102701120078563690,  -0.091777190566062927,  -0.093189977109432220,  -0.076448604464530945,   0.042494408786296844,   0.061926953494548798,   0.072122834622859955,  -0.042916599661111832,   0.040159467607736588,   0.003948521334677935,  -0.038303494453430176,   0.144594699144363403,  -0.135234326124191284,   0.074067637324333191,   0.064710043370723724,   0.096129126846790314,   0.223238021135330200,   0.034246362745761871,  -0.088947609066963196,   0.069960989058017731,   0.059509243816137314,   0.077478229999542236,   0.202306091785430908,   0.064793959259986877,   0.026784932240843773,   0.024330811575055122,   0.103965744376182556,  -0.081856206059455872,   0.048960607498884201,  -0.060131989419460297,   0.044010978192090988,  -0.097675606608390808,   0.043541885912418365,   0.059316124767065048,   0.023811757564544678,  -0.029534677043557167,   0.166268259286880493,   0.039482410997152328,  -0.068709522485733032,   0.096354544162750244,   0.003284036414697766,   0.177390098571777344,   0.011949460953474045,   0.020460603758692741,   0.102359600365161896,  -0.052287597209215164,   0.229863315820693970,
		   0.096316531300544739,  -0.080851703882217407,  -0.046239107847213745,   0.114093959331512451,  -0.124952830374240875,  -0.045328881591558456,  -0.060570400208234787,   0.090152353048324585,   0.050348863005638123,  -0.090240649878978729,   0.133677884936332703,   0.035155177116394043,   0.129061669111251831,   0.061982978135347366,   0.077379725873470306,  -0.118462115526199341,   0.059284754097461700,  -0.071322061121463776,  -0.071809753775596619,   0.025576014071702957,  -0.024180768057703972,  -0.120439149439334869,  -0.017617400735616684,  -0.137085303664207458,   0.025066077709197998,   0.064644925296306610,   0.004413934890180826,   0.078900054097175598,  -0.150382161140441895,  -0.118795752525329590,  -0.086484551429748535,   0.001725527108646929,  -0.064720645546913147,  -0.112521171569824219,   0.004284078255295753,  -0.016612621024250984,  -0.017587099224328995,   0.064029507339000702,  -0.034285094588994980,  -0.097471684217453003,  -0.021087396889925003,  -0.198774397373199463,  -0.170317232608795166,  -0.042397290468215942,  -0.075236074626445770,  -0.007897144183516502,  -0.129411697387695312,  -0.240784183144569397,
		  -0.068693175911903381,  -0.152719885110855103,   0.034742776304483414,  -0.054102681577205658,  -0.123917698860168457,   0.009902857244014740,  -0.051435213536024094,   0.000823445734567940,   0.005546332802623510,  -0.106277257204055786,  -0.185755908489227295,  -0.149505257606506348,  -0.047047160565853119,  -0.016621906310319901,   0.072853676974773407,   0.044487684965133667,  -0.164438009262084961,   0.023337850347161293,   0.080594524741172791,  -0.001179278711788356,   0.056010734289884567,  -0.073221802711486816,  -0.062964893877506256,   0.020386690273880959,  -0.077775940299034119,   0.093178883194923401,  -0.067702710628509521,  -0.035308111459016800,  -0.019902685657143593,   0.054094854742288589,   0.025895839557051659,   0.075145907700061798,   0.084723487496376038,  -0.163303002715110779,  -0.060194134712219238,   0.024162424728274345,   0.128555119037628174,  -0.041835252195596695,  -0.056393541395664215,  -0.116613805294036865,   0.061820093542337418,  -0.048893291503190994,  -0.088800728321075439,   0.098410882055759430,  -0.145915076136589050,  -0.069819889962673187,  -0.025107713416218758,   0.094053342938423157,
		  -0.034060001373291016,   0.100248888134956360,  -0.088871091604232788,   0.162616610527038574,   0.065026909112930298,   0.120017208158969879,  -0.077731117606163025,  -0.102662153542041779,   0.046176232397556305,   0.026897311210632324,   0.008736389689147472,   0.107412695884704590,   0.120279371738433838,   0.080312676727771759,  -0.095418445765972137,   0.111696943640708923,   0.067572236061096191,   0.115175068378448486,   0.113634318113327026,   0.087253503501415253,   0.156915679574012756,  -0.041231889277696609,   0.027391903102397919,   0.197545692324638367,  -0.028183085843920708,  -0.008424596861004829,  -0.101120389997959137,   0.156745016574859619,   0.003654502099379897,  -0.099220432341098785,   0.009527455084025860,   0.074139356613159180,  -0.054379023611545563,  -0.045038964599370956,   0.044168733060359955,   0.035389523953199387,  -0.053527869284152985,  -0.075324006378650665,   0.043088927865028381,   0.049118231981992722,   0.206689566373825073,   0.204870909452438354,  -0.075471952557563782,  -0.004739336203783751,   0.152917891740798950,   0.033107556402683258,   0.045159753412008286,   0.109645910561084747,
		   0.067439131438732147,  -0.012188101187348366,   0.000127011764561757,   0.078674696385860443,  -0.006813048385083675,  -0.065696693956851959,  -0.026315506547689438,  -0.053866602480411530,   0.095779046416282654,  -0.109925717115402222,  -0.094104118645191193,  -0.052458446472883224,   0.040162842720746994,   0.099848486483097076,   0.106022208929061890,  -0.067158438265323639,  -0.092028401792049408,   0.105054505169391632,   0.111222229897975922,  -0.054378502070903778,  -0.090924449265003204,   0.103816635906696320,   0.124502055346965790,   0.086241744458675385,  -0.067425891757011414,  -0.014094199053943157,   0.100845962762832642,   0.119730465114116669,   0.029485283419489861,  -0.098568670451641083,  -0.053786370903253555,  -0.042238548398017883,   0.144543856382369995,  -0.108184859156608582,  -0.005366817582398653,  -0.033822044730186462,  -0.028690213337540627,   0.018213886767625809,   0.147563084959983826,  -0.072078876197338104,   0.180306077003479004,   0.143634334206581116,   0.088150672614574432,   0.076768346130847931,   0.013765043579041958,   0.125514790415763855,  -0.012692230753600597,   0.132630601525306702,
		  -0.086589038372039795,   0.121931649744510651,   0.006931612268090248,   0.155580222606658936,  -0.066675372421741486,   0.117214292287826538,   0.182318136096000671,   0.084984920918941498,   0.084016427397727966,  -0.051503822207450867,   0.084675312042236328,   0.139528840780258179,  -0.063823327422142029,   0.134409502148628235,  -0.067180819809436798,  -0.071578912436962128,  -0.062890738248825073,   0.016149390488862991,  -0.015326544642448425,   0.002560165477916598,   0.039007104933261871,   0.068015903234481812,   0.118125915527343750,   0.164498552680015564,   0.055445555597543716,  -0.054497223347425461,  -0.088181190192699432,   0.008194446563720703,  -0.081263095140457153,   0.057834673672914505,  -0.099366217851638794,  -0.101285323500633240,   0.061270937323570251,   0.027885051444172859,   0.104570224881172180,   0.125834017992019653,  -0.109009884297847748,   0.045078132301568985,   0.015835622325539589,  -0.053071394562721252,   0.001602918258868158,  -0.110986731946468353,  -0.007422635797411203,   0.171807423233985901,   0.036849308758974075,  -0.048693105578422546,  -0.065864466130733490,  -0.021894836798310280,
		  -0.062627427279949188,   0.046400040388107300,  -0.099631354212760925,   0.170854181051254272,   0.065346270799636841,   0.153509527444839478,   0.031337745487689972,   0.144349336624145508,  -0.056635800749063492,   0.014625864103436470,  -0.100581362843513489,   0.090206891298294067,  -0.014655879698693752,  -0.020209589973092079,  -0.009355862624943256,  -0.113494470715522766,   0.083413422107696533,   0.030791224911808968,   0.087018646299839020,   0.111804388463497162,   0.187462002038955688,   0.117292769253253937,   0.057854387909173965,   0.217651143670082092,   0.185749933123588562,   0.052928056567907333,   0.156303808093070984,   0.090856224298477173,   0.166534647345542908,   0.145625337958335876,   0.017688117921352386,   0.071570746600627899,   0.150549262762069702,  -0.021870594471693039,  -0.011971394531428814,   0.041129853576421738,  -0.038733337074518204,  -0.000429674430051818,  -0.080038353800773621,  -0.055491648614406586,   0.077745102345943451,   0.116864688694477081,  -0.040603388100862503,  -0.020051939412951469,  -0.068375989794731140,   0.149047121405601501,   0.085977703332901001,  -0.022261887788772583,
	};
	static const double bias01[]=
	{
		   0.067900858819484711,   0.071879133582115173,   0.036661140620708466,  -0.013055714778602123,   0.130241245031356812,  -0.028642373159527779,   0.172458752989768982,   0.058011192828416824,  -0.042439505457878113,  -0.024122206494212151,  -0.022777885198593140,   0.002446191851049662,   0.011639439500868320,   0.097554810345172882,   0.086937569081783295,   0.021169373765587807,   0.124953605234622955,   0.059645477682352066,  -0.064095877110958099,   0.035594906657934189,   0.023782385513186455,  -0.045811995863914490,  -0.023131631314754486,   0.142123222351074219,   0.062571741640567780,   0.096727892756462097,   0.112990103662014008,   0.133822530508041382,   0.053090732544660568,  -0.110916793346405029,   0.089179120957851410,   0.112084470689296722,   0.010125493630766869,   0.089283987879753113,   0.158927887678146362,   0.094365321099758148,   0.074035458266735077,  -0.001835534465499222,  -0.010584095492959023,   0.091251291334629059,   0.098062194883823395,   0.091176688671112061,   0.164215654134750366,   0.041398145258426666,   0.059095900505781174,   0.008346536196768284,  -0.125438436865806580,   0.068651378154754639,
	};
	static const double weight02[]=
	{
		   0.021126762032508850,   0.105141557753086090,   0.091657668352127075,  -0.064708121120929718,   0.062814801931381226,   0.141434580087661743,   0.122931987047195435,   0.120731331408023834,  -0.160610899329185486,  -0.049528561532497406,  -0.140563532710075378,   0.071359202265739441,  -0.084239564836025238,  -0.020661864429712296,  -0.088772565126419067,  -0.069725140929222107,  -0.035800810903310776,  -0.084355771541595459,   0.145795464515686035,  -0.102514550089836121,   0.005344205070286989,   0.051986239850521088,   0.085222780704498291,   0.019230816513299942,   0.005756110884249210,  -0.024351216852664948,   0.057740744203329086,  -0.038776624947786331,   0.110151886940002441,  -0.030698662623763084,   0.153917074203491211,  -0.005896988790482283,  -0.095950923860073090,   0.145026504993438721,  -0.041324835270643234,   0.058092519640922546,   0.041094429790973663,   0.094418190419673920,  -0.140496715903282166,  -0.102310746908187866,  -0.126207470893859863,  -0.098030760884284973,  -0.076825127005577087,  -0.013179335743188858,  -0.062558740377426147,   0.002648606197908521,  -0.002344511449337006,  -0.012416548095643520,
		   0.136936590075492859,   0.047367002815008163,  -0.107752881944179535,   0.138286650180816650,  -0.004999152850359678,   0.076485618948936462,   0.057749480009078979,   0.044343519955873489,   0.108011677861213684,  -0.026931680738925934,   0.043040428310632706,   0.102730534970760345,  -0.037042807787656784,   0.082963295280933380,   0.069674849510192871,   0.009576144628226757,   0.103436931967735291,   0.129953548312187195,   0.162510290741920471,   0.120922461152076721,  -0.067539013922214508,  -0.094700314104557037,  -0.034635722637176514,   0.046048820018768311,   0.026840355247259140,   0.053822215646505356,   0.052835475653409958,   0.115600489079952240,  -0.062951669096946716,   0.133084565401077271,  -0.122072651982307434,  -0.050127968192100525,   0.148575469851493835,   0.053791344165802002,   0.100940510630607605,   0.058729834854602814,   0.137353971600532532,  -0.049330350011587143,   0.075033113360404968,   0.081495329737663269,  -0.006786639802157879,   0.146346360445022583,  -0.045134872198104858,  -0.093136884272098541,   0.137642756104469299,  -0.046657294034957886,  -0.072750203311443329,   0.138431966304779053,
		   0.094759382307529449,  -0.057061418890953064,   0.069258622825145721,   0.136646345257759094,  -0.049993071705102921,   0.049710504710674286,   0.052540943026542664,  -0.050816420465707779,  -0.026849683374166489,   0.096889935433864594,   0.123131915926933289,  -0.082308828830718994,   0.126766949892044067,  -0.011112492531538010,   0.105015501379966736,   0.043285466730594635,  -0.059358924627304077,   0.055021207779645920,  -0.116816632449626923,   0.000977100571617484,  -0.089731588959693909,  -0.033560443669557571,   0.079758994281291962,   0.132531404495239258,  -0.007542323321104050,  -0.106994457542896271,   0.144283011555671692,   0.039884567260742188,  -0.078582547605037689,   0.032327942550182343,  -0.147712066769599915,   0.100545585155487061,  -0.018974747508764267,   0.040744058787822723,   0.075207687914371490,   0.087101347744464874,   0.069971881806850433,  -0.161904320120811462,  -0.032285723835229874,   0.089539453387260437,   0.029954902827739716,  -0.027351161465048790,  -0.005523289088159800,   0.091511964797973633,  -0.085726395249366760,  -0.143258899450302124,   0.122713252902030945,   0.124019682407379150,
		  -0.011938591487705708,  -0.039746928960084915,  -0.107510246336460114,   0.047443114221096039,   0.097815312445163727,   0.123302884399890900,   0.017450118437409401,   0.160956159234046936,   0.143682256340980530,   0.018552552908658981,  -0.010967384092509747,  -0.046303186565637589,   0.112350732088088989,  -0.130038484930992126,   0.153891503810882568,  -0.036674965173006058,   0.071140468120574951,   0.034252185374498367,   0.010188030079007149,   0.073024362325668335,  -0.030567249283194542,   0.070561952888965607,   0.165562093257904053,  -0.062911182641983032,   0.164435163140296936,  -0.090523101389408112,   0.033579140901565552,  -0.010949698276817799,  -0.059541590511798859,   0.089651554822921753,   0.114863276481628418,  -0.045101027935743332,  -0.000713783432729542,  -0.119571477174758911,   0.049943733960390091,  -0.010033607482910156,   0.099968954920768738,  -0.029153704643249512,   0.022186413407325745,   0.125959202647209167,   0.157067716121673584,  -0.103813596069812775,   0.038828879594802856,  -0.089359045028686523,   0.009112430736422539,  -0.092416450381278992,  -0.068803712725639343,  -0.073023036122322083,
		   0.165446430444717407,   0.008917693980038166,  -0.007428215816617012,   0.051093000918626785,  -0.068848900496959686,  -0.105205044150352478,   0.053133074194192886,  -0.084732435643672943,   0.132871508598327637,   0.017289133742451668,   0.102503836154937744,   0.139352366328239441,  -0.077797010540962219,  -0.048905521631240845,   0.114056453108787537,   0.135934010148048401,  -0.071942105889320374,  -0.061616390943527222,  -0.108703248202800751,   0.011152753606438637,   0.091040119528770447,   0.081088364124298096,   0.155563533306121826,  -0.141225323081016541,   0.184667974710464478,  -0.002441769465804100,   0.048117168247699738,  -0.086600363254547119,   0.112159192562103271,  -0.071875810623168945,   0.018692076206207275,   0.071988135576248169,   0.147207498550415039,   0.144361615180969238,   0.052818071097135544,   0.173867389559745789,   0.032382279634475708,   0.082369156181812286,   0.037670861929655075,   0.121987275779247284,   0.029465913772583008,  -0.152347505092620850,   0.043861400336027145,   0.126382112503051758,   0.022111961618065834,   0.067528091371059418,   0.126785442233085632,  -0.116488739848136902,
		  -0.056999448686838150,  -0.131214812397956848,   0.081511050462722778,  -0.073164962232112885,  -0.046026468276977539,   0.079361386597156525,   0.074496574699878693,  -0.072999939322471619,   0.083019010722637177,  -0.056428894400596619,   0.030438864603638649,   0.069286070764064789,  -0.122295066714286804,   0.125378504395484924,   0.081926241517066956,   0.035301141440868378,  -0.016461640596389771,  -0.068770639598369598,  -0.054687321186065674,  -0.107938498258590698,   0.129655003547668457,  -0.059303313493728638,  -0.027271566912531853,   0.132445022463798523,  -0.121222309768199921,   0.036349553614854813,   0.120426669716835022,   0.008868885226547718,   0.179841205477714539,   0.126967683434486389,   0.049750939011573792,   0.016387755051255226,   0.139361158013343811,   0.152363821864128113,  -0.032674361020326614,  -0.040258754044771194,   0.108649395406246185,  -0.067041210830211639,   0.053917877376079559,   0.161966651678085327,   0.090787157416343689,   0.097666539251804352,  -0.154983103275299072,   0.063159696757793427,   0.147221818566322327,  -0.059063781052827835,   0.053801499307155609,   0.140154555439949036,
		   0.014434848912060261,   0.133495911955833435,  -0.059460457414388657,  -0.031433828175067902,   0.162631541490554810,   0.115572363138198853,   0.137456312775611877,  -0.017799979075789452,  -0.044803243130445480,   0.017987348139286041,   0.122928671538829803,  -0.004263755865395069,   0.181908354163169861,   0.051296945661306381,  -0.052918300032615662,   0.061815969645977020,   0.106528185307979584,   0.008653544820845127,   0.135162413120269775,   0.075017735362052917,   0.066204622387886047,   0.157050222158432007,   0.134208679199218750,   0.150737255811691284,  -0.110770538449287415,   0.062868334352970123,   0.119821585714817047,   0.120751120150089264,   0.193644478917121887,  -0.077504917979240417,   0.053725924342870712,   0.151254400610923767,  -0.065066657960414886,  -0.051237482577562332,  -0.134956538677215576,  -0.110256992280483246,   0.134531348943710327,  -0.061325747519731522,   0.131171211600303650,   0.145234808325767517,  -0.013392973691225052,   0.105524100363254547,  -0.081176102161407471,  -0.092343486845493317,  -0.076465450227260590,   0.111178971827030182,   0.135042563080787659,   0.048335853964090347,
		   0.000281111977528781,  -0.151032179594039917,  -0.106746904551982880,   0.063855230808258057,   0.037879124283790588,  -0.041821125894784927,   0.082996025681495667,  -0.055958181619644165,   0.064098119735717773,  -0.093030236661434174,  -0.004013733007013798,  -0.124129898846149445,  -0.003149773227050900,  -0.052897375077009201,  -0.133817404508590698,  -0.021661005914211273,   0.102101676166057587,   0.009844144806265831,   0.173735931515693665,  -0.092855371534824371,   0.073647029697895050,  -0.003583163488656282,   0.030852725729346275,   0.112078517675399780,  -0.115987636148929596,  -0.134791612625122070,  -0.112571142613887787,   0.021895457059144974,  -0.061123486608266830,  -0.094795562326908112,   0.066675715148448944,  -0.132910460233688354,  -0.044499170035123825,  -0.084898836910724640,   0.073078624904155731,  -0.150702595710754395,   0.103192858397960663,  -0.082945972681045532,   0.117068417370319366,  -0.066573478281497955,  -0.140665128827095032,   0.175717756152153015,   0.054170202463865280,  -0.082279562950134277,   0.118960045278072357,   0.029341911897063255,   0.101658925414085388,  -0.129166349768638611,
		   0.059193495661020279,  -0.067966878414154053,   0.041426554322242737,  -0.020781179890036583,   0.033309370279312134,  -0.051610648632049561,   0.099361196160316467,  -0.020929213613271713,   0.168529555201530457,   0.025112139061093330,   0.008237293921411037,  -0.030881628394126892,   0.048484392464160919,  -0.036183021962642670,   0.129986032843589783,  -0.067913427948951721,  -0.151106461882591248,   0.087905243039131165,   0.058644663542509079,   0.004384169820696115,   0.003534861607477069,   0.190271139144897461,   0.065579265356063843,  -0.017405927181243896,   0.049022857099771500,   0.086780853569507599,  -0.112585015594959259,   0.062596902251243591,  -0.119606994092464447,   0.139335244894027710,  -0.139928624033927917,   0.124780833721160889,   0.069375202059745789,  -0.021184705197811127,   0.094691075384616852,   0.134561091661453247,   0.057644657790660858,  -0.135468855500221252,  -0.022586198523640633,   0.026095557957887650,   0.060215376317501068,   0.036811795085668564,  -0.056788656860589981,  -0.116147346794605255,  -0.039675433188676834,   0.100654631853103638,  -0.054034899920225143,  -0.091113045811653137,
		  -0.047583971172571182,  -0.063555300235748291,  -0.033786129206418991,  -0.087441526353359222,   0.078737325966358185,   0.013096633367240429,   0.007654680870473385,   0.044080123305320740,   0.108331508934497833,   0.010842467658221722,  -0.053082630038261414,  -0.113976635038852692,  -0.032701831310987473,   0.218452081084251404,   0.071572959423065186,  -0.002060924656689167,  -0.079950548708438873,  -0.167368888854980469,  -0.062344014644622803,   0.097904317080974579,   0.117283433675765991,  -0.015021083876490593,  -0.122948445379734039,   0.164963349699974060,   0.054692223668098450,  -0.011403086595237255,   0.050677496939897537,  -0.101114012300968170,   0.050822105258703232,   0.004617163911461830,   0.114709101617336273,  -0.075195036828517914,  -0.146873235702514648,   0.059887051582336426,  -0.094472564756870270,  -0.079630389809608459,   0.172853589057922363,   0.032074283808469772,   0.155297622084617615,  -0.008955208584666252,   0.104203984141349792,   0.055183507502079010,  -0.160230994224548340,   0.022422013804316521,   0.188094228506088257,   0.125173732638359070,   0.063564889132976532,   0.040683802217245102,
		   0.067204780876636505,   0.104626037180423737,  -0.105626478791236877,  -0.062730371952056885,  -0.139213532209396362,  -0.026403443887829781,   0.119834870100021362,  -0.088980205357074738,  -0.048606958240270615,   0.134430438280105591,   0.085680842399597168,  -0.133606985211372375,  -0.119703590869903564,  -0.141122788190841675,  -0.064986944198608398,  -0.070357330143451691,  -0.010982462204992771,   0.087173990905284882,   0.001987806288525462,  -0.130599126219749451,   0.104714095592498779,  -0.041898187249898911,  -0.139152318239212036,  -0.026761798188090324,  -0.052268266677856445,   0.020026022568345070,   0.042055908590555191,   0.110527880489826202,  -0.104160316288471222,  -0.069329574704170227,  -0.120063908398151398,  -0.103133760392665863,   0.028511919081211090,  -0.056104041635990143,   0.116285458207130432,  -0.037009317427873611,  -0.158348008990287781,  -0.097092084586620331,  -0.125459432601928711,  -0.117373310029506683,  -0.143585458397865295,  -0.110433518886566162,   0.143926039338111877,   0.026623910292983055,   0.000307589332805946,  -0.166156813502311707,   0.055279098451137543,  -0.044268634170293808,
		   0.056761227548122406,   0.038553379476070404,  -0.001723671914078295,  -0.061825569719076157,  -0.071452178061008453,  -0.063178382813930511,  -0.097118265926837921,  -0.107954733073711395,  -0.168664202094078064,   0.025331478565931320,  -0.032434485852718353,  -0.026951344683766365,  -0.153892219066619873,   0.110425308346748352,   0.019017795100808144,   0.200045317411422729,   0.180630818009376526,   0.028238689526915550,   0.179766312241554260,  -0.051762085407972336,   0.132045149803161621,  -0.136945962905883789,  -0.058949921280145645,   0.061433285474777222,  -0.167379125952720642,   0.056959580630064011,   0.039806243032217026,   0.117463603615760803,   0.084294036030769348,  -0.168843105435371399,  -0.097407653927803040,  -0.034530986100435257,  -0.076040826737880707,   0.136640593409538269,  -0.079856134951114655,   0.088552899658679962,   0.072477489709854126,  -0.113131694495677948,   0.049648817628622055,   0.097412280738353729,   0.013438469730317593,   0.211457550525665283,  -0.010256843641400337,   0.059229176491498947,   0.064717441797256470,   0.187682464718818665,  -0.035922933369874954,   0.073251925408840179,
		  -0.015101339668035507,   0.122329726815223694,   0.020290756598114967,  -0.052222795784473419,   0.169143348932266235,   0.020123364403843880,  -0.058215979486703873,  -0.051749616861343384,  -0.036108344793319702,   0.015037938021123409,   0.056923750787973404,   0.075048655271530151,   0.008903102949261665,   0.027830144390463829,  -0.013069927692413330,   0.008130421862006187,  -0.088207930326461792,   0.011317782104015350,  -0.111221157014369965,   0.081827290356159210,  -0.027587983757257462,   0.131400063633918762,  -0.075957432389259338,   0.029323969036340714,  -0.010635511949658394,  -0.088915944099426270,  -0.027363434433937073,   0.056816432625055313,  -0.026301296427845955,   0.050037872046232224,  -0.129798963665962219,  -0.114729195833206177,   0.013215465471148491,  -0.127962768077850342,  -0.055810205638408661,   0.082791492342948914,   0.070667475461959839,   0.006518759764730930,  -0.005757082719355822,   0.067328028380870819,  -0.038136754184961319,   0.022692216560244560,   0.081555485725402832,  -0.065291665494441986,  -0.041524428874254227,  -0.042461782693862915,   0.124411314725875854,   0.077943749725818634,
		  -0.054339554160833359,  -0.062144454568624496,  -0.087657734751701355,   0.071278758347034454,  -0.093709468841552734,   0.127312913537025452,  -0.020097827538847923,   0.126996919512748718,   0.112581200897693634,   0.076913535594940186,  -0.007988418452441692,  -0.078113146126270294,   0.054611638188362122,  -0.151647925376892090,   0.127083986997604370,  -0.109376884996891022,  -0.111405059695243835,   0.111807040870189667,  -0.048292025923728943,   0.082656726241111755,  -0.004160312470048666,  -0.018455008044838905,   0.129083797335624695,  -0.116962544620037079,   0.148480907082557678,   0.182608038187026978,  -0.097949318587779999,  -0.075379930436611176,   0.112133741378784180,   0.062460947781801224,  -0.021299436688423157,   0.155414775013923645,   0.172593310475349426,  -0.046590313315391541,   0.070802070200443268,   0.058703757822513580,   0.105425938963890076,   0.112005814909934998,   0.036931484937667847,   0.076242759823799133,   0.023237563669681549,  -0.045972291380167007,   0.169654607772827148,  -0.023201597854495049,  -0.039741702377796173,  -0.078367151319980621,  -0.065184481441974640,  -0.118729084730148315,
		   0.073013268411159515,   0.198791354894638062,   0.081127822399139404,   0.167508825659751892,   0.096329323947429657,   0.074434846639633179,  -0.059431534260511398,   0.133307620882987976,   0.084181614220142365,   0.042327549308538437,   0.000332375173456967,   0.127597078680992126,  -0.091078259050846100,  -0.173805400729179382,   0.118738859891891479,   0.098884798586368561,   0.014450115151703358,   0.077346295118331909,  -0.140922620892524719,   0.025763483718037605,  -0.015591699630022049,   0.044803190976381302,   0.165047198534011841,   0.123639710247516632,   0.145601585507392883,  -0.041679050773382187,   0.051072057336568832,   0.045004330575466156,   0.030522810295224190,   0.137594938278198242,   0.124514564871788025,   0.166050538420677185,  -0.055268742144107819,   0.046070657670497894,   0.028474375605583191,   0.139143541455268860,   0.001027143676765263,   0.003245532279834151,   0.010511703789234161,  -0.102521792054176331,   0.177908331155776978,  -0.139952287077903748,  -0.019229540601372719,  -0.105063289403915405,  -0.086219556629657745,  -0.010262533091008663,   0.026206709444522858,   0.162782251834869385,
		  -0.061893377453088760,  -0.023104904219508171,   0.054932761937379837,  -0.095607392489910126,   0.070395067334175110,  -0.134690135717391968,  -0.134852349758148193,   0.049533858895301819,   0.086200781166553497,  -0.114373005926609039,   0.124075330793857574,  -0.129262566566467285,  -0.090420141816139221,   0.091003209352493286,   0.003164146328344941,  -0.024523831903934479,  -0.103574745357036591,  -0.112682759761810303,  -0.113554254174232483,  -0.037543933838605881,  -0.023795776069164276,  -0.122974298894405365,   0.021778166294097900,   0.023030182346701622,   0.107404530048370361,  -0.007483730558305979,  -0.164071276783943176,  -0.148599237203598022,  -0.044307593256235123,  -0.068313032388687134,  -0.040181118994951248,  -0.093661218881607056,   0.003880741074681282,   0.028119025751948357,   0.051244262605905533,  -0.012296628206968307,  -0.174771755933761597,   0.031242774799466133,   0.118348106741905212,  -0.042200028896331787,   0.083098739385604858,  -0.072111003100872040,  -0.051505334675312042,  -0.144641369581222534,  -0.063987098634243011,   0.058577839285135269,   0.128952160477638245,  -0.103607371449470520,
		   0.046772558242082596,   0.129957884550094604,   0.131594806909561157,   0.127831280231475830,   0.041729688644409180,   0.034496150910854340,   0.035174041986465454,  -0.090811699628829956,   0.111599490046501160,  -0.025459581986069679,  -0.067828647792339325,   0.091848619282245636,  -0.006900950800627470,   0.031066562980413437,  -0.078052237629890442,  -0.141246676445007324,   0.099435664713382721,   0.175640508532524109,  -0.159988880157470703,  -0.040315300226211548,  -0.020724333822727203,  -0.038281764835119247,   0.050126869231462479,  -0.045085765421390533,   0.148250475525856018,   0.043676923960447311,   0.104818366467952728,   0.045850444585084915,   0.095097035169601440,  -0.122737318277359009,   0.113208815455436707,  -0.085085086524486542,   0.064924784004688263,  -0.118695415556430817,   0.098674297332763672,   0.122958801686763763,  -0.021148752421140671,  -0.083815209567546844,   0.103483162820339203,   0.159750923514366150,   0.160382136702537537,  -0.085844688117504120,   0.052576202899217606,  -0.091213323175907135,   0.135399743914604187,   0.077726125717163086,   0.022997906431555748,   0.065590985119342804,
		   0.165036842226982117,  -0.012653155252337456,   0.128686413168907166,   0.176642894744873047,   0.008412947878241539,   0.129169717431068420,   0.119442194700241089,   0.054500024765729904,  -0.067759767174720764,   0.043649289757013321,   0.127202734351158142,   0.155186280608177185,   0.088992565870285034,  -0.167911306023597717,  -0.000439200783148408,   0.021935284137725830,  -0.111970618367195129,   0.175138652324676514,   0.089561887085437775,   0.107295773923397064,  -0.014946031384170055,   0.114229820668697357,   0.106724850833415985,   0.071718975901603699,   0.127452522516250610,  -0.065437063574790955,  -0.110152907669544220,   0.088285170495510101,  -0.062851823866367340,  -0.005866240710020065,  -0.006833439227193594,   0.146337404847145081,  -0.084321826696395874,  -0.013171047903597355,  -0.028002945706248283,   0.024787217378616333,  -0.000651251990348101,   0.097160801291465759,   0.077967897057533264,   0.147366434335708618,   0.010603664442896843,  -0.085584878921508789,   0.159551680088043213,  -0.028058893978595734,  -0.057280693203210831,  -0.081776194274425507,   0.135853931307792664,  -0.018022352829575539,
		   0.142380639910697937,  -0.036283265799283981,   0.098404787480831146,  -0.069388240575790405,  -0.129475533962249756,  -0.053368043154478073,   0.082063972949981689,   0.054409880191087723,   0.171199247241020203,   0.039782620966434479,  -0.037842251360416412,  -0.079157665371894836,   0.153778761625289917,   0.039347562938928604,   0.089771777391433716,   0.091890662908554077,  -0.167645037174224854,  -0.113526366651058197,  -0.096077658236026764,   0.104531601071357727,  -0.081178084015846252,  -0.077187821269035339,  -0.064332909882068634,  -0.138829395174980164,   0.043065629899501801,   0.015153140760958195,  -0.023971334099769592,  -0.055298209190368652,  -0.052887376397848129,  -0.033278308808803558,  -0.051230583339929581,  -0.090484023094177246,   0.158628836274147034,   0.130207926034927368,  -0.078144729137420654,   0.135258272290229797,  -0.098763749003410339,   0.121495708823204041,  -0.013966713100671768,  -0.090620070695877075,   0.133268073201179504,  -0.013640221208333969,   0.141054585576057434,   0.057175826281309128,  -0.187511414289474487,  -0.063111849129199982,  -0.014918257482349873,  -0.054020021110773087,
		  -0.084653213620185852,   0.112020052969455719,  -0.043368458747863770,  -0.086353778839111328,   0.058444313704967499,  -0.109095305204391479,   0.029969720169901848,   0.022931367158889771,  -0.008085991255939007,  -0.130122065544128418,   0.105311803519725800,  -0.024387300014495850,  -0.043507389724254608,   0.071496844291687012,  -0.013871303759515285,  -0.068409271538257599,   0.015939930453896523,   0.003200642997398973,  -0.010852163657546043,  -0.120563320815563202,  -0.066257715225219727,   0.040923785418272018,   0.169333428144454956,  -0.043533734977245331,   0.113220006227493286,   0.015246217139065266,  -0.094397649168968201,   0.093453034758567810,   0.026570351794362068,   0.103528335690498352,  -0.097594536840915680,  -0.120926611125469208,  -0.039173115044832230,  -0.104850001633167267,  -0.025969408452510834,  -0.065952718257904053,  -0.046798214316368103,   0.070945411920547485,  -0.038552418351173401,   0.009146320633590221,   0.029533080756664276,  -0.001505716471001506,  -0.034393761307001114,   0.093491844832897186,   0.050242852419614792,  -0.142080098390579224,   0.087838113307952881,  -0.121878638863563538,
		   0.076496012508869171,   0.031053753569722176,   0.084126435220241547,  -0.025673732161521912,  -0.051282204687595367,   0.080835305154323578,  -0.031552113592624664,  -0.008138731122016907,   0.143581464886665344,   0.043149545788764954,  -0.017944840714335442,  -0.031064478680491447,   0.092100381851196289,   0.011768574826419353,   0.033454403281211853,  -0.127526342868804932,  -0.031021462753415108,   0.138301402330398560,  -0.086370356380939484,   0.052128236740827560,   0.024643911048769951,   0.076497353613376617,   0.087728649377822876,   0.069067776203155518,   0.187269940972328186,  -0.089102491736412048,   0.048629283905029297,  -0.052247233688831329,   0.032722368836402893,   0.127136290073394775,   0.122770771384239197,   0.092644147574901581,   0.149144351482391357,   0.055456828325986862,   0.066222503781318665,   0.022897949442267418,   0.055022444576025009,   0.104472249746322632,  -0.038938812911510468,   0.145838052034378052,  -0.055754594504833221,  -0.104972854256629944,   0.112889237701892853,   0.111133717000484467,  -0.078221134841442108,  -0.055935464799404144,  -0.088828720152378082,   0.034781225025653839,
		  -0.038360498845577240,   0.041709665209054947,  -0.057589646428823471,   0.094784446060657501,   0.121868446469306946,  -0.042572520673274994,  -0.003900421084836125,   0.155465468764305115,   0.111682571470737457,  -0.030660625547170639,   0.191552519798278809,   0.095766291022300720,   0.181686311960220337,  -0.071780554950237274,  -0.024157697334885597,  -0.091105990111827850,  -0.077468052506446838,   0.072757191956043243,  -0.052855208516120911,   0.141224965453147888,   0.091966494917869568,   0.045330334454774857,   0.126017630100250244,   0.062274582684040070,   0.129114776849746704,   0.157572239637374878,   0.003777644829824567,   0.107810713350772858,   0.048850003629922867,  -0.025376565754413605,   0.041640475392341614,   0.106905587017536163,  -0.080980665981769562,  -0.135685533285140991,   0.020923178642988205,  -0.033745937049388885,  -0.064587146043777466,  -0.104304730892181396,  -0.083441980183124542,   0.079127341508865356,   0.021559871733188629,  -0.127527534961700439,  -0.073634654283523560,   0.118956044316291809,   0.127471566200256348,   0.030993802472949028,  -0.096162915229797363,  -0.016001420095562935,
		  -0.052890431135892868,  -0.141381353139877319,  -0.001745490357279778,  -0.072405517101287842,  -0.098485134541988373,  -0.073282197117805481,   0.150608345866203308,   0.181425064802169800,   0.055701699107885361,  -0.104302026331424713,   0.052236396819353104,   0.001947515876963735,  -0.080938115715980530,   0.219138577580451965,  -0.000629272661171854,  -0.068503223359584808,   0.186536133289337158,  -0.123588941991329193,   0.045427534729242325,   0.110740475356578827,   0.081649936735630035,  -0.045014500617980957,   0.031106108799576759,   0.152657911181449890,   0.007719034329056740,  -0.072204150259494781,   0.137016803026199341,   0.155199080705642700,  -0.012140562757849693,   0.074349150061607361,   0.202413648366928101,   0.040917448699474335,   0.060557805001735687,   0.169384181499481201,  -0.063566416501998901,  -0.144017457962036133,  -0.081892266869544983,   0.049811135977506638,  -0.082564398646354675,   0.023356556892395020,   0.073268502950668335,   0.120753273367881775,  -0.065392658114433289,   0.175828382372856140,   0.082423657178878784,   0.188227444887161255,  -0.059424333274364471,   0.062501631677150726,
		  -0.024820268154144287,   0.126492232084274292,   0.008184054866433144,   0.046060331165790558,   0.111563280224800110,   0.090502604842185974,   0.074818164110183716,  -0.027106391265988350,   0.131168738007545471,  -0.051943939179182053,   0.021107509732246399,   0.126992970705032349,   0.076313599944114685,  -0.015264013782143593,  -0.128231331706047058,   0.041294958442449570,   0.167183443903923035,  -0.028569672256708145,   0.026342192664742470,   0.134618684649467468,   0.147955223917961121,  -0.067328147590160370,  -0.134775027632713318,  -0.073933131992816925,   0.092812620103359222,  -0.020366808399558067,  -0.073994621634483337,   0.144493594765663147,   0.161871671676635742,  -0.007441086694598198,  -0.011429868638515472,  -0.111265189945697784,  -0.117347717285156250,   0.136710599064826965,   0.044084288179874420,   0.020019765943288803,  -0.011799694038927555,  -0.072010874748229980,  -0.050933599472045898,  -0.027620006352663040,  -0.128831386566162109,   0.114750832319259644,  -0.130279451608657837,   0.183685705065727234,   0.175186008214950562,   0.201059892773628235,  -0.029922287911176682,  -0.024804353713989258,
		  -0.035292241722345352,  -0.128577739000320435,  -0.090985536575317383,   0.113880693912506104,  -0.110036633908748627,  -0.168769448995590210,  -0.016134025529026985,  -0.031624618917703629,   0.111657038331031799,   0.035529427230358124,   0.143818110227584839,   0.040121294558048248,  -0.104289814829826355,  -0.035964813083410263,   0.056920353323221207,   0.043618239462375641,  -0.082906328141689301,  -0.136497035622596741,  -0.116645194590091705,   0.083020418882369995,   0.032809752970933914,   0.063179999589920044,   0.044887736439704895,   0.016925670206546783,  -0.065440267324447632,  -0.107965566217899323,  -0.087118595838546753,  -0.001640065689571202,  -0.122828528285026550,  -0.116733580827713013,  -0.153887063264846802,  -0.027475316077470779,   0.027898980304598808,  -0.066846378147602081,  -0.087596438825130463,   0.055003687739372253,  -0.094347499310970306,  -0.152322858572006226,  -0.003861777484416962,   0.086444310843944550,   0.045707996934652328,  -0.106863036751747131,   0.137837752699851990,   0.102254249155521393,   0.079965189099311829,  -0.106720604002475739,  -0.036064621061086655,  -0.000586110923904926,
		   0.147881418466567993,   0.180398494005203247,   0.062273468822240829,  -0.009067299775779247,  -0.034909423440694809,   0.115783907473087311,   0.108045585453510284,   0.118064217269420624,   0.039340306073427200,   0.055613752454519272,   0.004866470582783222,   0.082479260861873627,   0.111776396632194519,  -0.084770627319812775,   0.066301837563514709,   0.003934013191610575,  -0.076808929443359375,  -0.080529183149337769,   0.099833019077777863,   0.041477199643850327,  -0.087606430053710938,   0.131270915269851685,   0.114163681864738464,   0.138863623142242432,  -0.074987508356571198,   0.127278253436088562,  -0.052746370434761047,  -0.020936813205480576,  -0.053135585039854050,  -0.098287321627140045,  -0.105246089398860931,   0.102493911981582642,   0.074526168406009674,   0.122987873852252960,   0.013753926381468773,  -0.006242053117603064,  -0.082973219454288483,   0.076608024537563324,  -0.077696248888969421,   0.035697467625141144,  -0.049366533756256104,  -0.085799224674701691,   0.150203168392181396,   0.093371234834194183,  -0.111532710492610931,  -0.065176188945770264,   0.131956294178962708,   0.098120667040348053,
		   0.151576235890388489,  -0.062213197350502014,  -0.061025716364383698,  -0.062025293707847595,  -0.073856584727764130,  -0.096865430474281311,  -0.084758982062339783,   0.097624361515045166,  -0.154935270547866821,  -0.088034957647323608,  -0.088220581412315369,   0.016073195263743401,  -0.102404125034809113,   0.114415153861045837,  -0.104257941246032715,  -0.073380306363105774,   0.013210783712565899,   0.043990377336740494,   0.102712303400039673,   0.022564323619008064,   0.102500282227993011,  -0.129904001951217651,  -0.005659949965775013,  -0.076989531517028809,  -0.043275661766529083,  -0.032234318554401398,   0.114101581275463104,   0.053304441273212433,   0.127452373504638672,   0.060900919139385223,   0.090296037495136261,   0.002807332901284099,  -0.118247129023075104,  -0.029393067583441734,   0.091455534100532532,  -0.130830332636833191,  -0.019370885565876961,   0.056505311280488968,   0.114304065704345703,   0.157949686050415039,   0.005622474942356348,   0.036405730992555618,  -0.073421478271484375,   0.185169756412506104,  -0.076001711189746857,  -0.027069356292486191,   0.009472424164414406,   0.070814527571201324,
		   0.055960919708013535,   0.044282056391239166,   0.014328937046229839,  -0.021368375048041344,   0.102198973298072815,  -0.129615187644958496,   0.046089593321084976,  -0.087275646626949310,  -0.082754030823707581,  -0.047896053642034531,   0.102951057255268097,  -0.045158218592405319,  -0.027941731736063957,  -0.019420985132455826,  -0.106095127761363983,  -0.038567058742046356,  -0.029330667108297348,  -0.033184126019477844,   0.087764509022235870,  -0.062903866171836853,  -0.068723231554031372,  -0.106259867548942566,   0.010321511887013912,  -0.122744090855121613,  -0.083374954760074615,   0.124154999852180481,   0.112250477075576782,   0.031576398760080338,  -0.069182038307189941,   0.118298441171646118,  -0.119329802691936493,  -0.044441498816013336,  -0.000292470591375604,   0.034438971430063248,   0.054139234125614166,   0.037465676665306091,  -0.075478263199329376,   0.083343639969825745,   0.047149218618869781,  -0.014522632583975792,   0.110032722353935242,  -0.122962974011898041,  -0.052213169634342194,   0.096827410161495209,   0.058143172413110733,   0.031741652637720108,   0.160144686698913574,  -0.084364719688892365,
		   0.062668122351169586,   0.199686661362648010,   0.011047962121665478,   0.082063019275665283,  -0.021825518459081650,   0.129818484187126160,   0.088760904967784882,  -0.079701296985149384,   0.118872053921222687,  -0.057757489383220673,   0.072709791362285614,  -0.032696984708309174,   0.123648464679718018,  -0.158347755670547485,   0.126414984464645386,   0.096133649349212646,  -0.090726576745510101,   0.118823535740375519,  -0.028370391577482224,   0.154566287994384766,   0.079705782234668732,   0.086648128926753998,   0.006909514777362347,   0.079841189086437225,   0.193022310733795166,   0.107042215764522552,  -0.037376169115304947,  -0.025773085653781891,   0.115265361964702606,  -0.019753446802496910,   0.102815054357051849,  -0.054645098745822906,  -0.016894161701202393,   0.073792524635791779,   0.206018701195716858,  -0.075173869729042053,  -0.115027748048305511,   0.004309875424951315,  -0.121584884822368622,  -0.024639787152409554,   0.150417029857635498,  -0.063126116991043091,   0.157914221286773682,   0.129044890403747559,  -0.041684526950120926,  -0.031083408743143082,  -0.002071554772555828,   0.092184834182262421,
		   0.157354041934013367,  -0.066653609275817871,   0.035986278206110001,  -0.002376051852479577,   0.010453705675899982,  -0.052893001586198807,  -0.027124531567096710,  -0.042281318455934525,   0.073279447853565216,   0.033499088138341904,   0.100209549069404602,   0.109539695084095001,   0.027917798608541489,   0.153706789016723633,   0.057289510965347290,   0.055924396961927414,   0.154763609170913696,   0.143466472625732422,   0.015909383073449135,   0.106000445783138275,  -0.066349819302558899,   0.019700020551681519,   0.038173388689756393,   0.165044039487838745,  -0.109678350389003754,  -0.095704555511474609,  -0.044665709137916565,  -0.006822347175329924,   0.139291614294052124,  -0.058779440820217133,   0.008856075815856457,   0.079920969903469086,  -0.033336840569972992,   0.172449022531509399,   0.100256159901618958,  -0.121816910803318024,   0.172758713364601135,   0.006178404670208693,   0.185833632946014404,   0.142130032181739807,  -0.033267397433519363,  -0.064792193472385406,  -0.022444974631071091,   0.129965126514434814,  -0.063026957213878632,   0.117928519845008850,   0.061050117015838623,  -0.107148729264736176,
		   0.043923422694206238,   0.127284720540046692,   0.051899410784244537,  -0.061899702996015549,   0.008839807473123074,   0.122877888381481171,  -0.044629111886024475,   0.075561523437500000,   0.105006910860538483,   0.055507160723209381,   0.077129535377025604,  -0.129812806844711304,   0.106979109346866608,  -0.173474773764610291,   0.015967750921845436,  -0.055041432380676270,   0.072619535028934479,   0.107218854129314423,  -0.117650501430034637,   0.146323278546333313,  -0.087309390306472778,   0.087148904800415039,   0.022734615951776505,   0.126377329230308533,  -0.006189441308379173,   0.018061917275190353,  -0.087545327842235565,   0.025703476741909981,   0.110980778932571411,  -0.037112154066562653,   0.005156950093805790,   0.003203829983249307,   0.115381143987178802,  -0.020887007936835289,   0.118568584322929382,   0.039615955203771591,   0.049289505928754807,  -0.076105646789073944,  -0.121926873922348022,  -0.024147830903530121,   0.158794209361076355,   0.071673169732093811,   0.150127977132797241,   0.106825120747089386,  -0.129219576716423035,   0.068709574639797211,   0.114107951521873474,   0.044585928320884705,
		   0.071552425622940063,  -0.089962318539619446,   0.083249576389789581,   0.052046153694391251,  -0.034677989780902863,   0.061776321381330490,   0.107022784650325775,   0.120899014174938202,  -0.057911526411771774,   0.065788917243480682,   0.117967270314693451,  -0.050187535583972931,  -0.023049501702189445,  -0.040846548974514008,  -0.127170920372009277,   0.166529282927513123,   0.087354861199855804,   0.029206372797489166,   0.083609320223331451,   0.168167352676391602,   0.134832456707954407,  -0.115808762609958649,   0.005332456901669502,   0.124148622155189514,   0.132450237870216370,   0.080241508781909943,   0.191255122423171997,   0.072858117520809174,   0.061643160879611969,  -0.052868779748678207,   0.180969312787055969,  -0.072338260710239410,  -0.026329744607210159,   0.181560322642326355,   0.048927363008260727,  -0.087905533611774445,   0.089771546423435211,  -0.056827429682016373,   0.182437956333160400,   0.010924120433628559,  -0.087300918996334076,   0.000135023161419667,  -0.128888025879859924,  -0.055413272231817245,  -0.071924768388271332,   0.114114493131637573,   0.116965718567371368,   0.097069941461086273,
		  -0.009732894599437714,   0.160289391875267029,   0.077447414398193359,   0.065144442021846771,  -0.087207064032554626,  -0.020310722291469574,  -0.093331880867481232,  -0.023581877350807190,  -0.009305795654654503,   0.125962555408477783,   0.070557691156864166,   0.027500947937369347,  -0.050309889018535614,   0.084388293325901031,   0.165507361292839050,  -0.105032473802566528,   0.015559711493551731,  -0.075339101254940033,   0.122047081589698792,   0.149620413780212402,  -0.097561575472354889,   0.064893640577793121,   0.109719045460224152,  -0.045779179781675339,   0.032978661358356476,   0.165981888771057129,   0.021163213998079300,  -0.011358038522303104,  -0.089674174785614014,   0.085921131074428558,  -0.058148942887783051,  -0.008086361922323704,   0.107352323830127716,   0.067545324563980103,   0.190349414944648743,   0.074882753193378448,   0.029761962592601776,  -0.106583036482334137,  -0.100721977651119232,   0.042711406946182251,   0.029405497014522552,   0.018949462100863457,   0.202353268861770630,   0.042961210012435913,  -0.141481205821037292,  -0.129819959402084351,  -0.035826243460178375,  -0.023135930299758911,
		   0.089134074747562408,   0.045334435999393463,  -0.103192038834095001,   0.107585147023200989,  -0.061561003327369690,   0.122935339808464050,  -0.047423377633094788,   0.126667723059654236,   0.182808563113212585,   0.059725619852542877,   0.076026007533073425,   0.022463226690888405,   0.078646771609783173,  -0.082183450460433960,   0.015275269746780396,  -0.047534365206956863,   0.065695643424987793,   0.128177687525749207,  -0.062552936375141144,   0.052747838199138641,   0.149889945983886719,   0.112419918179512024,  -0.036379873752593994,   0.015373573638498783,   0.062624804675579071,   0.164480209350585938,  -0.041569951921701431,  -0.067116379737854004,  -0.095484398305416107,   0.105228230357170105,  -0.063617311418056488,   0.000723772740457207,   0.112216979265213013,   0.120209552347660065,   0.086531378328800201,   0.170335710048675537,   0.004944547079503536,  -0.113794915378093719,  -0.018041074275970459,   0.116537563502788544,   0.139679327607154846,   0.119953684508800507,   0.179563611745834351,   0.128804281353950500,  -0.074578322470188141,  -0.135858625173568726,   0.005639062728732824,   0.074494093656539917,
		  -0.081562377512454987,   0.115249395370483398,   0.041221063584089279,   0.004770532716065645,   0.070204056799411774,   0.150039106607437134,   0.081257067620754242,   0.131303101778030396,  -0.019412335008382797,  -0.039352715015411377,   0.092250205576419830,   0.069540672004222870,   0.160537973046302795,  -0.026872720569372177,  -0.067011028528213501,   0.211230248212814331,   0.077796928584575653,  -0.025160167366266251,   0.063564032316207886,   0.162296399474143982,  -0.004835125058889389,  -0.042432107031345367,   0.099580936133861542,   0.093925610184669495,   0.094096861779689789,  -0.090163514018058777,  -0.006052281241863966,   0.101780362427234650,   0.070575088262557983,   0.033394414931535721,  -0.020968051627278328,   0.118809297680854797,   0.145966768264770508,   0.158492982387542725,  -0.135949179530143738,  -0.098265275359153748,   0.165361478924751282,   0.191478028893470764,   0.063742421567440033,   0.161036074161529541,  -0.092891827225685120,   0.153463438153266907,   0.098110027611255646,   0.003794936928898096,   0.178372949361801147,  -0.001412432640790939,  -0.136232495307922363,   0.151356607675552368,
		   0.011501128785312176,  -0.088421694934368134,   0.087263762950897217,  -0.170876264572143555,   0.118149839341640472,  -0.070171006023883820,   0.058515638113021851,  -0.051116380840539932,  -0.087576337158679962,  -0.053295761346817017,  -0.079209342598915100,  -0.015213547274470329,   0.076385401189327240,  -0.017076529562473297,  -0.010257523506879807,   0.122014023363590240,  -0.068689890205860138,   0.002702749799937010,   0.204058974981307983,   0.036869756877422333,   0.008103613741695881,  -0.041467837989330292,  -0.113428696990013123,   0.068051166832447052,  -0.160695359110832214,   0.022271899506449699,   0.157689720392227173,   0.104456298053264618,   0.150590285658836365,  -0.110997810959815979,  -0.029122594743967056,  -0.121793352067470551,  -0.096865154802799225,  -0.052728202193975449,  -0.018677817657589912,  -0.033295441418886185,   0.060476485639810562,  -0.015309754759073257,   0.158995389938354492,   0.146397069096565247,   0.030122598633170128,   0.038288775831460953,   0.071062214672565460,  -0.083353884518146515,  -0.084134086966514587,  -0.063262946903705597,   0.107195973396301270,   0.142298728227615356,
		  -0.102990843355655670,   0.064053662121295929,   0.096664249897003174,   0.143325820565223694,  -0.080103017389774323,  -0.038860108703374863,  -0.011470720171928406,   0.135098397731781006,  -0.104187242686748505,   0.087015613913536072,   0.185876533389091492,  -0.085378214716911316,  -0.035516839474439621,  -0.156538724899291992,  -0.010688350535929203,  -0.069769635796546936,  -0.017071781679987907,  -0.002045209752395749,  -0.140291333198547363,  -0.048780202865600586,   0.020840903744101524,   0.051832225173711777,  -0.042883619666099548,   0.082108862698078156,   0.138348326086997986,   0.050616484135389328,  -0.054995775222778320,   0.046614620834589005,   0.045859422534704208,   0.011104121804237366,  -0.045862600207328796,   0.050110518932342529,  -0.119286403059959412,  -0.126097857952117920,   0.076193667948246002,   0.155364632606506348,  -0.063266627490520477,   0.103608697652816772,  -0.035087700933218002,  -0.036605801433324814,  -0.040920615196228027,  -0.054795920848846436,  -0.093263722956180573,   0.072184599936008453,  -0.131142660975456238,   0.025315832346677780,   0.041683863848447800,   0.168611407279968262,
		   0.035882823169231415,   0.121240861713886261,  -0.061824075877666473,  -0.123769015073776245,   0.001360365655273199,   0.133546233177185059,  -0.019027940928936005,  -0.026586096733808517,  -0.044163577258586884,  -0.133313342928886414,   0.043991044163703918,   0.096741214394569397,  -0.040213268250226974,  -0.008646049536764622,  -0.018126677721738815,  -0.114402987062931061,   0.004117472097277641,  -0.094600096344947815,   0.074333406984806061,   0.000207155433599837,  -0.150662690401077271,   0.007446964271366596,   0.021441256627440453,   0.063119247555732727,   0.097456052899360657,   0.074097074568271637,   0.070358321070671082,  -0.157962977886199951,   0.076490454375743866,  -0.116641491651535034,  -0.085950836539268494,  -0.117059677839279175,  -0.096236459910869598,  -0.084281988441944122,   0.047879658639431000,  -0.071710564196109772,   0.030268501490354538,   0.047730725258588791,   0.086797930300235748,   0.015154525637626648,  -0.041055776178836823,  -0.016985256224870682,   0.058453783392906189,  -0.041072215884923935,   0.064880542457103729,  -0.053317937999963760,  -0.050540231168270111,  -0.145232364535331726,
		  -0.084383688867092133,   0.017673328518867493,  -0.070780500769615173,  -0.071066752076148987,   0.151588380336761475,   0.135401338338851929,  -0.063380122184753418,  -0.115532405674457550,  -0.117619387805461884,  -0.151807591319084167,   0.003522538812831044,  -0.070902891457080841,   0.100798465311527252,  -0.046753782778978348,   0.108691617846488953,   0.121453627943992615,   0.117837682366371155,  -0.141748920083045959,   0.065426297485828400,   0.114951454102993011,   0.024726618081331253,  -0.004959193989634514,  -0.081457920372486115,   0.070237152278423309,   0.018965892493724823,   0.004440129268914461,   0.071888394653797150,   0.130545705556869507,   0.001809174544177949,   0.022005110979080200,   0.093898423016071320,   0.089993879199028015,   0.056390617042779922,  -0.088026598095893860,  -0.110974781215190887,   0.012428075075149536,  -0.013078861869871616,   0.163142040371894836,  -0.058207441121339798,  -0.054251525551080704,  -0.053801558911800385,  -0.059470009058713913,   0.014626762829720974,   0.148551046848297119,  -0.003093667095527053,   0.078825712203979492,  -0.160404026508331299,   0.021076714619994164,
		  -0.018287897109985352,   0.194080218672752380,   0.097025156021118164,  -0.072045922279357910,  -0.072261169552803040,   0.074863478541374207,  -0.118907593190670013,   0.111198090016841888,   0.073330491781234741,   0.133537888526916504,   0.102735333144664764,   0.139724612236022949,  -0.020773911848664284,  -0.086202017962932587,   0.180637702345848083,   0.070351846516132355,   0.127003297209739685,   0.160805225372314453,   0.001644078409299254,  -0.110009901225566864,   0.087370269000530243,   0.166119858622550964,   0.170301795005798340,   0.112741954624652863,   0.193044275045394897,   0.060287445783615112,   0.117421962320804596,   0.036375198513269424,   0.074370361864566803,  -0.076484426856040955,  -0.042923722416162491,  -0.013788142241537571,   0.137096017599105835,   0.010973718017339706,   0.125525996088981628,   0.190858319401741028,  -0.039684876799583435,   0.071996390819549561,   0.147940039634704590,   0.094402387738227844,   0.157976344227790833,  -0.036245554685592651,   0.121907860040664673,  -0.014653845690190792,  -0.054600860923528671,  -0.047947403043508530,  -0.078479476273059845,   0.004701138474047184,
		   0.151941612362861633,   0.052395526319742203,   0.104829460382461548,   0.117960825562477112,  -0.011135498993098736,  -0.007868431508541107,   0.112074419856071472,  -0.055166307836771011,   0.158830523490905762,  -0.077237017452716827,   0.005872785579413176,   0.037471406161785126,   0.138904094696044922,  -0.182620450854301453,  -0.012953515164554119,   0.122471913695335388,   0.056240245699882507,   0.108497194945812225,  -0.094099208712577820,  -0.026950355619192123,  -0.062788248062133789,   0.145817890763282776,   0.178514674305915833,  -0.050072502344846725,   0.122637823224067688,   0.110392779111862183,   0.014152538962662220,  -0.086623333394527435,   0.007431741338223219,   0.104028433561325073,  -0.121572099626064301,   0.124507762491703033,   0.148927450180053711,  -0.019723778590559959,   0.030340313911437988,   0.046331681311130524,   0.088639616966247559,   0.051908824592828751,   0.000699226628057659,   0.056694027036428452,  -0.040916942059993744,  -0.077579528093338013,   0.044959753751754761,  -0.028875671327114105,  -0.017784515395760536,  -0.055338736623525620,   0.001265105674974620,  -0.089577771723270416,
		  -0.133824914693832397,   0.027022048830986023,  -0.014657697640359402,  -0.160886391997337341,  -0.067003369331359863,   0.124320805072784424,   0.146929919719696045,  -0.052190121263265610,  -0.073468945920467377,  -0.024244006723165512,   0.074702709913253784,   0.134035676717758179,  -0.112703226506710052,   0.177792355418205261,   0.026372190564870834,   0.006771192885935307,   0.159044355154037476,  -0.027490328997373581,  -0.009028428234159946,  -0.025053551420569420,   0.089536771178245544,   0.090539053082466125,  -0.020047230646014214,  -0.112245239317417145,  -0.102313898503780365,   0.023351196199655533,  -0.015627596527338028,   0.057451691478490829,   0.099517799913883209,  -0.118435621261596680,   0.142769604921340942,   0.025893040001392365,   0.100616984069347382,  -0.059016685932874680,  -0.046687006950378418,   0.033150941133499146,   0.003716752165928483,   0.064789280295372009,  -0.079712018370628357,   0.162744909524917603,   0.006181270815432072,   0.201930403709411621,  -0.159912973642349243,   0.144176408648490906,   0.157791092991828918,  -0.049051508307456970,  -0.017290716990828514,   0.119160674512386322,
		   0.114966310560703278,  -0.019333411008119583,   0.151349484920501709,   0.026481134817004204,  -0.019192712381482124,   0.026388250291347504,  -0.008394285105168819,   0.029563354328274727,   0.054899588227272034,  -0.115679360926151276,  -0.086945161223411560,  -0.074696227908134460,  -0.133837416768074036,  -0.042437359690666199,  -0.075313337147235870,   0.131888762116432190,   0.033216133713722229,   0.065554507076740265,   0.056581012904644012,  -0.040244676172733307,  -0.053209010511636734,  -0.030786836519837379,  -0.065727762877941132,   0.026365347206592560,  -0.168332993984222412,   0.042961385101079941,  -0.001235696370713413,   0.121582292020320892,  -0.036523979157209396,   0.101706914603710175,  -0.060486044734716415,  -0.008733754977583885,  -0.022121809422969818,  -0.027869850397109985,  -0.140968292951583862,   0.097887873649597168,   0.090351052582263947,  -0.059780359268188477,   0.131525889039039612,  -0.006855376996099949,  -0.124487012624740601,   0.132398143410682678,   0.027642790228128433,   0.010576063767075539,   0.133145377039909363,   0.036357104778289795,  -0.161537989974021912,   0.115752555429935455,
		   0.130432143807411194,   0.182165935635566711,   0.053920004516839981,   0.036607995629310608,   0.032623413950204849,   0.017019528895616531,  -0.106222301721572876,   0.119481801986694336,  -0.057387154549360275,   0.074386969208717346,   0.100122332572937012,   0.155649244785308838,   0.014529003761708736,  -0.146855890750885010,   0.088654905557632446,   0.074927888810634613,  -0.019271448254585266,   0.077708728611469269,   0.095767796039581299,   0.024780973792076111,  -0.036582913249731064,  -0.030635006725788116,  -0.018226213753223419,   0.002969988156110048,  -0.049419801682233810,   0.117372952401638031,  -0.034797891974449158,   0.085285887122154236,  -0.069996781647205353,   0.068935208022594452,  -0.065068788826465607,   0.160931542515754700,  -0.076292634010314941,   0.017488500103354454,   0.183873832225799561,   0.121167257428169250,   0.074263542890548706,   0.005968797486275434,  -0.048711106181144714,   0.079972296953201294,   0.140035420656204224,   0.014935596846044064,  -0.012175948359072208,  -0.045013677328824997,  -0.029238319024443626,   0.002228075172752142,   0.005197777412831783,   0.122544601559638977,
		   0.005742023698985577,   0.124303996562957764,   0.155780151486396790,   0.019251562654972076,   0.089446574449539185,   0.054464254528284073,  -0.005946204531937838,  -0.025761893019080162,   0.102391138672828674,  -0.077546253800392151,   0.038255978375673294,  -0.075501516461372375,  -0.096425503492355347,   0.099670767784118652,   0.140741646289825439,   0.048996560275554657,  -0.033671900629997253,  -0.001806949148885906,   0.087465867400169373,  -0.017371904104948044,   0.139062300324440002,  -0.070993065834045410,   0.196681633591651917,   0.167736276984214783,   0.082532152533531189,   0.047922357916831970,   0.131269469857215881,   0.010390488430857658,  -0.098765306174755096,   0.136715382337570190,   0.082014180719852448,   0.129127308726310730,   0.177950516343116760,   0.113326355814933777,   0.060070019215345383,   0.120684757828712463,   0.033991601318120956,   0.015532207675278187,   0.110147692263126373,   0.095780804753303528,   0.055613733828067780,  -0.064869135618209839,   0.127928391098976135,  -0.015183790586888790,  -0.008309615775942802,  -0.000360955513315275,  -0.009562213905155659,   0.094613984227180481,
		  -0.046190440654754639,  -0.070822358131408691,   0.088000901043415070,  -0.065330423414707184,   0.113600648939609528,  -0.089560449123382568,  -0.062902510166168213,   0.047399539500474930,   0.003699186956509948,  -0.012837319634854794,   0.079905904829502106,  -0.084921978414058685,  -0.122541882097721100,  -0.030664063990116119,   0.008932010270655155,   0.203159585595130920,  -0.016783058643341064,   0.007549040019512177,   0.145281016826629639,  -0.095021098852157593,  -0.057547725737094879,  -0.113746248185634613,  -0.011151040904223919,  -0.089298918843269348,  -0.177569344639778137,  -0.098720148205757141,   0.184566214680671692,  -0.012427373789250851,   0.010794247500598431,  -0.124149993062019348,   0.149799436330795288,  -0.045911516994237900,  -0.005030141677707434,  -0.087839528918266296,  -0.036973245441913605,  -0.134191602468490601,   0.150200784206390381,  -0.074964217841625214,   0.148208662867546082,  -0.065350145101547241,  -0.062893837690353394,   0.075280182063579559,  -0.039353251457214355,  -0.010397708974778652,  -0.012590177357196808,  -0.006237020250409842,  -0.113478623330593109,   0.069016024470329285,
		  -0.062142226845026016,  -0.054561331868171692,   0.150359004735946655,  -0.039802610874176025,   0.131615012884140015,   0.023358004167675972,  -0.041403912007808685,  -0.100035443902015686,  -0.023023357614874840,  -0.019179159775376320,  -0.032168284058570862,  -0.027948904782533646,   0.018837442621588707,   0.009761069901287556,  -0.101774692535400391,   0.161079823970794678,   0.135100394487380981,   0.028998708352446556,   0.078635185956954956,   0.120305322110652924,  -0.053998265415430069,   0.019205944612622261,  -0.036903046071529388,   0.145551830530166626,   0.083840467035770416,   0.157534286379814148,   0.051613021641969681,   0.133531391620635986,  -0.068558961153030396,   0.130647107958793640,   0.050514671951532364,   0.125815093517303467,   0.160990223288536072,  -0.016987510025501251,   0.113639257848262787,   0.046283625066280365,  -0.038291100412607193,   0.139125540852546692,   0.105052925646305084,  -0.077666006982326508,   0.162703558802604675,   0.089053176343441010,   0.126072198152542114,  -0.087644539773464203,   0.141941592097282410,   0.137411445379257202,   0.075506635010242462,   0.028155714273452759,
		  -0.023767644539475441,   0.155541390180587769,   0.131290689110755920,   0.064533643424510956,  -0.006433788686990738,   0.137213408946990967,   0.093546114861965179,   0.119142748415470123,   0.199010014533996582,   0.168231785297393799,   0.124374113976955414,  -0.078272797167301178,  -0.079504430294036865,   0.067949362099170685,   0.068916954100131989,  -0.055191647261381149,   0.008694757707417011,   0.039631612598896027,  -0.011423854157328606,   0.125622794032096863,   0.148269146680831909,   0.036284428089857101,  -0.045013032853603363,  -0.016881391406059265,   0.013728101737797260,   0.037190683186054230,  -0.025440569967031479,   0.113652557134628296,  -0.105742499232292175,  -0.081802882254123688,   0.079821966588497162,   0.034157861024141312,   0.123670071363449097,   0.039351932704448700,   0.038956705480813980,  -0.077976353466510773,  -0.138872623443603516,  -0.044538162648677826,   0.020399764180183411,   0.093324601650238037,   0.115327470004558563,  -0.062958300113677979,   0.104786895215511322,   0.063600257039070129,  -0.089963227510452271,   0.030103726312518120,   0.136841163039207458,  -0.037838749587535858,
	};
	static const double bias02[]=
	{
		   0.093823127448558807,   0.109387964010238647,  -0.097560226917266846,   0.039532817900180817,   0.057389810681343079,   0.013363858684897423,   0.167607486248016357,  -0.124363027513027191,  -0.057174172252416611,  -0.075074113905429840,  -0.103791110217571259,  -0.027664791792631149,   0.118837483227252960,   0.122957393527030945,   0.001524808350950480,  -0.004125546663999557,  -0.010130700655281544,   0.122465036809444427,  -0.026400109753012657,  -0.034022189676761627,   0.048645179718732834,  -0.104092888534069061,   0.060734432190656662,  -0.097601532936096191,  -0.049020566046237946,   0.106649041175842285,   0.173958092927932739,  -0.008930903859436512,   0.090960755944252014,   0.045887835323810577,   0.060475241392850876,  -0.017820913344621658,   0.035059768706560135,   0.065335132181644440,   0.185803592205047607,   0.172580048441886902,  -0.042516946792602539,  -0.128700375556945801,   0.080458380281925201,  -0.090911224484443665,   0.055905487388372421,   0.038145232945680618,   0.137114003300666809,  -0.017633847892284393,  -0.010651409626007080,   0.057077422738075256,   0.178988516330718994,  -0.107578471302986145,
	};
	static const double weight03[]=
	{
		  -0.118278846144676208,  -0.097933888435363770,   0.119583502411842346,   0.137995600700378418,  -0.088732123374938965,  -0.016220621764659882,  -0.121559172868728638,   0.008845733478665352,   0.109602846205234528,  -0.000027619727916317,   0.126392379403114319,   0.090313337743282318,  -0.095152877271175385,   0.121819071471691132,  -0.093065708875656128,   0.135098978877067566,   0.083492383360862732,  -0.111689619719982147,   0.117100566625595093,   0.005005770362913609,   0.045126751065254211,   0.083050191402435303,  -0.165881261229515076,   0.058926586061716080,  -0.109526321291923523,  -0.035868842154741287,  -0.116470873355865479,   0.017770679667592049,  -0.011766391806304455,   0.004337958991527557,   0.012390610761940479,  -0.078720092773437500,   0.034816116094589233,  -0.118627972900867462,  -0.141486421227455139,  -0.065681718289852142,   0.141184747219085693,   0.044220108538866043,  -0.093427374958992004,  -0.106188908219337463,   0.000980915385298431,   0.055596593767404556,  -0.159953951835632324,  -0.017305873334407806,   0.081072270870208740,   0.055750619620084763,  -0.033484097570180893,  -0.022938951849937439,
		  -0.023936463519930840,   0.111939825117588043,   0.116861619055271149,  -0.093804948031902313,   0.139084234833717346,  -0.065306738018989563,  -0.015639757737517357,  -0.086103364825248718,   0.131285056471824646,  -0.087804690003395081,   0.017997493967413902,  -0.112744309008121490,   0.037141233682632446,  -0.000490158738102764,   0.100720115005970001,   0.086731873452663422,  -0.026218712329864502,   0.164598673582077026,  -0.048848267644643784,   0.109192796051502228,   0.060191180557012558,   0.179175779223442078,   0.004829663317650557,  -0.093154944479465485,  -0.000923965068068355,   0.154155150055885315,  -0.028059605509042740,  -0.008399651385843754,   0.083518281579017639,   0.144484668970108032,   0.045618098229169846,   0.047165900468826294,   0.076450020074844360,   0.003571374109014869,   0.085409231483936310,   0.094985865056514740,   0.127545475959777832,  -0.130307972431182861,   0.007360233459621668,  -0.048273738473653793,   0.107868574559688568,  -0.159237399697303772,   0.057821687310934067,   0.103063143789768219,  -0.045669924467802048,   0.079765848815441132,  -0.043184529989957809,   0.123464860022068024,
		   0.033492676913738251,  -0.037899743765592575,  -0.046716172248125076,  -0.014997808262705803,   0.052796877920627594,  -0.125483587384223938,  -0.010831456631422043,  -0.020955670624971390,   0.044201325625181198,   0.004854105878621340,  -0.048915755003690720,   0.056876126676797867,   0.024745261296629906,   0.087067976593971252,  -0.044096026569604874,   0.026644797995686531,   0.095686212182044983,   0.150250092148780823,  -0.000997182214632630,  -0.001571703935042024,   0.138172358274459839,   0.000381400343030691,   0.042425047606229782,   0.016819581389427185,   0.131276085972785950,  -0.036377359181642532,  -0.090091496706008911,  -0.045915760099887848,   0.150525569915771484,  -0.010967858135700226,  -0.072050958871841431,   0.076774083077907562,  -0.063027791678905487,   0.120028272271156311,   0.020571926608681679,  -0.109782353043556213,   0.056816510856151581,  -0.137158215045928955,   0.018059439957141876,   0.172431603074073792,   0.073615871369838715,   0.007738347165286541,  -0.061747044324874878,   0.136776357889175415,  -0.015320381149649620,  -0.161873683333396912,  -0.031263604760169983,   0.150975599884986877,
		   0.065034911036491394,   0.051330152899026871,  -0.073359787464141846,   0.127504006028175354,   0.123401030898094177,  -0.030439319089055061,   0.186786264181137085,   0.017362104728817940,  -0.078945435583591461,   0.048920135945081711,  -0.053870584815740585,   0.046933282166719437,   0.038221344351768494,   0.061709593981504440,  -0.104071162641048431,  -0.157385334372520447,  -0.095671772956848145,   0.027252905070781708,  -0.016518583521246910,  -0.117989391088485718,   0.097879648208618164,   0.077476583421230316,   0.185937255620956421,   0.146580055356025696,  -0.184215158224105835,   0.147490829229354858,  -0.028791870921850204,   0.118633508682250977,   0.070722125470638275,   0.120226174592971802,   0.108365871012210846,   0.189625978469848633,   0.084455676376819611,   0.101925589144229889,  -0.077598474919795990,   0.165327504277229309,   0.063071452081203461,   0.021263934671878815,  -0.086845919489860535,  -0.068968959152698517,   0.081059888005256653,  -0.035562705248594284,   0.106502510607242584,  -0.082510091364383698,   0.148093044757843018,   0.165099769830703735,   0.151231139898300171,  -0.033363956958055496,
		   0.019163966178894043,   0.027348875999450684,   0.033056389540433884,  -0.053568713366985321,  -0.019897686317563057,  -0.126294136047363281,   0.038324296474456787,  -0.110454991459846497,  -0.041664265096187592,   0.083528079092502594,   0.021922163665294647,  -0.031458776444196701,  -0.029283292591571808,   0.112761944532394409,   0.102052062749862671,  -0.145279973745346069,  -0.017045106738805771,   0.137098491191864014,   0.097451426088809967,   0.131397277116775513,   0.130510300397872925,   0.121143810451030731,  -0.062070596963167191,  -0.041125033050775528,  -0.081846028566360474,   0.170318886637687683,  -0.087455607950687408,   0.053935755044221878,   0.166770309209823608,   0.034755807369947433,   0.043934524059295654,   0.160623535513877869,   0.048979550600051880,   0.127878457307815552,   0.049355242401361465,  -0.144556686282157898,   0.076489128172397614,  -0.158706933259963989,   0.008858083747327328,   0.182908535003662109,   0.042311381548643112,  -0.014540129341185093,  -0.145514965057373047,   0.156423196196556091,  -0.041897144168615341,  -0.026855813339352608,   0.004966314882040024,   0.136501491069793701,
		  -0.048349589109420776,   0.018936745822429657,  -0.099719591438770294,   0.023787347599864006,  -0.051665719598531723,   0.175911083817481995,   0.170631960034370422,   0.049760919064283371,   0.015983931720256805,   0.005865744780749083,  -0.178030088543891907,   0.154775142669677734,  -0.104591734707355499,  -0.075525797903537750,  -0.128225266933441162,   0.047564733773469925,   0.005336041096597910,  -0.079908162355422974,   0.136469498276710510,  -0.098194628953933716,   0.100922137498855591,  -0.063772946596145630,   0.103556729853153229,   0.142080545425415039,  -0.045398436486721039,  -0.116565570235252380,   0.115128926932811737,  -0.115126244723796844,  -0.016026286408305168,   0.032701380550861359,   0.070189967751502991,  -0.080070942640304565,  -0.001669300254434347,   0.017357498407363892,   0.108611583709716797,  -0.001586747006513178,   0.107874207198619843,   0.042723618447780609,   0.155525729060173035,   0.033562600612640381,   0.022687483578920364,   0.124471925199031830,   0.036458794027566910,   0.031122280284762383,  -0.086505927145481110,   0.131125539541244507,   0.000682352518197149,   0.133562088012695312,
		   0.122418411076068878,  -0.033923983573913574,   0.089375555515289307,  -0.034204628318548203,   0.069790221750736237,   0.066077426075935364,  -0.121062345802783966,  -0.080673478543758392,   0.121095560491085052,  -0.011435152031481266,   0.116714820265769958,   0.000122260404168628,  -0.123403206467628479,  -0.072485059499740601,   0.113480642437934875,  -0.124265722930431366,  -0.022063367068767548,   0.024008058011531830,  -0.101289205253124237,  -0.038444839417934418,   0.027084721252322197,   0.141105204820632935,  -0.090142287313938141,  -0.021308587864041328,  -0.103609435260295868,   0.160065680742263794,   0.082978352904319763,  -0.053807895630598068,  -0.110840573906898499,  -0.116736523807048798,  -0.025388808920979500,   0.069593816995620728,  -0.012679832987487316,   0.069156393408775330,  -0.004219876602292061,  -0.136524513363838196,   0.128227859735488892,   0.063360445201396942,  -0.037116020917892456,   0.057233445346355438,  -0.025068374350667000,  -0.119448989629745483,   0.113845877349376678,   0.111853592097759247,  -0.067344658076763153,  -0.133742392063140869,   0.061514366418123245,  -0.002065131673589349,
		   0.066915944218635559,  -0.115232527256011963,  -0.046431142836809158,   0.021256376057863235,  -0.035799678415060043,   0.025445841252803802,   0.082205168902873993,   0.051593560725450516,  -0.054353039711713791,  -0.058927401900291443,   0.121273137629032135,   0.023700313642621040,   0.094988025724887848,  -0.083533167839050293,  -0.150065138936042786,   0.033814970403909683,   0.017534058541059494,   0.046744618564844131,   0.047522269189357758,  -0.082634374499320984,  -0.095298551023006439,   0.101073503494262695,  -0.019295932725071907,   0.039784733206033707,  -0.041478954255580902,   0.067019544541835785,  -0.020757909864187241,  -0.068628802895545959,   0.037851706147193909,   0.022666074335575104,  -0.142403200268745422,   0.068604819476604462,   0.037266559898853302,   0.105653509497642517,  -0.137607812881469727,  -0.068640679121017456,   0.023895736783742905,   0.020134394988417625,   0.054913092404603958,  -0.092013642191886902,  -0.106762997806072235,  -0.044100947678089142,  -0.092159211635589600,  -0.072750531136989594,  -0.034104079008102417,  -0.066976509988307953,   0.052344433963298798,   0.031567841768264771,
		   0.100746005773544312,  -0.097242735326290131,  -0.048952415585517883,  -0.128395825624465942,  -0.024738168343901634,  -0.092541888356208801,  -0.037352506071329117,   0.082165583968162537,  -0.088692739605903625,  -0.050423555076122284,  -0.032550677657127380,   0.057105362415313721,   0.105355836451053619,  -0.004716169089078903,  -0.146050065755844116,   0.019517509266734123,  -0.051549185067415237,  -0.147745773196220398,  -0.089605949819087982,  -0.069689676165580750,  -0.032329119741916656,   0.003831124631687999,  -0.014598189853131771,  -0.134810417890548706,   0.017490826547145844,   0.124429136514663696,  -0.046265576034784317,   0.128137260675430298,   0.072637908160686493,   0.005352519452571869,  -0.119343310594558716,   0.106579057872295380,  -0.002764948643743992,  -0.140108779072761536,   0.084518775343894958,  -0.093394398689270020,   0.099473744630813599,   0.074927344918251038,  -0.016082886606454849,  -0.025008808821439743,  -0.155825734138488770,  -0.094839558005332947,   0.040931198745965958,   0.032433494925498962,   0.101998470723628998,  -0.044610250741243362,  -0.139881268143653870,  -0.038561992347240448,
		  -0.107240900397300720,   0.062654301524162292,   0.030572514981031418,   0.157542467117309570,   0.157305508852005005,   0.094776719808578491,   0.048861183226108551,   0.028990563005208969,   0.062192749232053757,  -0.031672000885009766,   0.120344839990139008,   0.074384428560733795,   0.164166361093521118,   0.128101497888565063,   0.142726182937622070,  -0.006140618119388819,   0.108817942440509796,   0.127485707402229309,   0.053033433854579926,   0.115693420171737671,   0.090480260550975800,  -0.038675308227539062,   0.112340822815895081,   0.083126112818717957,  -0.106770113110542297,   0.034425105899572372,  -0.195740014314651489,  -0.006323298439383507,   0.185778066515922546,   0.072518408298492432,   0.083309665322303772,   0.009793137200176716,  -0.056219790130853653,   0.068654894828796387,  -0.108014397323131561,  -0.051905341446399689,   0.049984589219093323,  -0.020075332373380661,   0.087725542485713959,   0.188850075006484985,   0.171831071376800537,   0.006841105874627829,  -0.089820310473442078,   0.128676295280456543,   0.046507369726896286,  -0.148752793669700623,   0.127298533916473389,   0.035175129771232605,
		  -0.065589047968387604,  -0.063548035919666290,  -0.102881714701652527,  -0.008062394335865974,  -0.150725886225700378,  -0.127507090568542480,  -0.101723782718181610,   0.066095113754272461,  -0.111565425992012024,  -0.053239207714796066,  -0.013390200212597847,  -0.097044706344604492,  -0.042968418449163437,   0.042394306510686874,  -0.036317583173513412,   0.000288199342321604,  -0.005575648508965969,  -0.121645793318748474,  -0.027616636827588081,   0.015618751756846905,  -0.003304327838122845,  -0.119861081242561340,  -0.021698344498872757,  -0.152988776564598083,  -0.009973092004656792,   0.101275682449340820,   0.036733500659465790,   0.116227805614471436,  -0.017797762528061867,   0.022344816476106644,   0.112036183476448059,  -0.125114977359771729,  -0.078316375613212585,  -0.062818974256515503,  -0.111367009580135345,  -0.031183900311589241,  -0.049930837005376816,  -0.039298947900533676,   0.019556973129510880,  -0.109825454652309418,  -0.089160136878490448,  -0.061592344194650650,   0.088814474642276764,  -0.110167793929576874,   0.094828784465789795,  -0.124841719865798950,  -0.145261764526367188,  -0.027717540040612221,
		  -0.018495036289095879,  -0.064414903521537781,   0.145950928330421448,   0.044750072062015533,  -0.082431867718696594,   0.116144165396690369,  -0.114102467894554138,  -0.000088420863903593,  -0.022334801033139229,  -0.052615366876125336,  -0.052243188023567200,  -0.038254357874393463,   0.002555046696215868,  -0.002136084716767073,   0.140832185745239258,   0.128006964921951294,   0.092897735536098480,   0.063711337745189667,  -0.056299213320016861,  -0.003458212362602353,   0.066768825054168701,   0.104856453835964203,   0.126103490591049194,   0.081919997930526733,   0.061294227838516235,   0.013747520744800568,  -0.037641424685716629,   0.007767657283693552,  -0.052412152290344238,  -0.076140642166137695,  -0.017552940174937248,   0.035621773451566696,   0.090802684426307678,  -0.025697996839880943,   0.139859944581985474,   0.063156418502330780,   0.035871636122465134,   0.004413718823343515,  -0.112139828503131866,   0.092299297451972961,  -0.064554937183856964,  -0.171252444386482239,   0.066023178398609161,   0.077167257666587830,  -0.016401601955294609,  -0.103590101003646851,   0.039377119392156601,   0.086381085216999054,
		   0.116320729255676270,   0.171149671077728271,   0.031885989010334015,  -0.133140966296195984,   0.088473677635192871,  -0.001019073533825576,   0.173728823661804199,   0.026703583076596260,  -0.102575816214084625,   0.194258585572242737,  -0.074301920831203461,  -0.033332005143165588,   0.073495522141456604,  -0.063486725091934204,   0.046406231820583344,   0.063626900315284729,   0.123695820569992065,  -0.102877460420131683,  -0.032443974167108536,  -0.115789845585823059,   0.014144905842840672,   0.032816682010889053,   0.211679890751838684,   0.128661394119262695,   0.006875921972095966,   0.089277811348438263,   0.056604925543069839,  -0.111794658005237579,  -0.106618747115135193,  -0.045099485665559769,  -0.145672217011451721,   0.080076076090335846,  -0.035458654165267944,  -0.102841258049011230,   0.173246785998344421,   0.089077524840831757,  -0.116853252053260803,  -0.054886110126972198,   0.161310702562332153,   0.031955122947692871,   0.098307080566883087,   0.168457373976707458,  -0.058806061744689941,   0.081700995564460754,  -0.029529072344303131,   0.186333581805229187,   0.079218454658985138,   0.037755358964204788,
		  -0.007333486806601286,   0.156482324004173279,  -0.062356118112802505,   0.181696534156799316,   0.126579284667968750,  -0.063712410628795624,   0.157091677188873291,   0.081312648952007294,  -0.008506027981638908,  -0.163915410637855530,  -0.000638265046291053,  -0.039408408105373383,  -0.082192219793796539,   0.127705752849578857,   0.146813184022903442,  -0.046891327947378159,   0.124493584036827087,   0.171721771359443665,   0.066612474620342255,   0.043989725410938263,   0.012657683342695236,   0.103506267070770264,  -0.118614263832569122,   0.008140498772263527,   0.025832861661911011,   0.094823829829692841,   0.040382295846939087,   0.075602538883686066,   0.176655694842338562,   0.077269695699214935,   0.072169333696365356,  -0.046861406415700912,  -0.060094624757766724,   0.194001823663711548,  -0.073792546987533569,   0.007016263436526060,   0.114821650087833405,   0.098803400993347168,   0.051615405827760696,  -0.012306138873100281,   0.065419897437095642,  -0.156410083174705505,  -0.104447200894355774,   0.159622848033905029,   0.166186958551406860,   0.104873739182949066,  -0.027893120422959328,   0.134272754192352295,
		  -0.079549826681613922,  -0.048438914120197296,  -0.091072551906108856,  -0.010612821206450462,   0.064442299306392670,   0.058431763201951981,   0.076793029904365540,   0.068183593451976776,  -0.155636116862297058,   0.096963368356227875,   0.007491325959563255,  -0.025899920612573624,   0.085992008447647095,   0.120909638702869415,  -0.061649937182664871,   0.112586110830307007,  -0.167276844382286072,  -0.131653830409049988,  -0.140781566500663757,  -0.068882569670677185,   0.049883827567100525,   0.089394517242908478,  -0.066415369510650635,   0.161411806941032410,  -0.159716844558715820,  -0.056147363036870956,  -0.064062915742397308,   0.044992305338382721,  -0.161693245172500610,   0.113369986414909363,  -0.057748701423406601,   0.040003422647714615,   0.086269713938236237,  -0.184205129742622375,  -0.042563468217849731,   0.158557042479515076,  -0.062903061509132385,  -0.126548528671264648,   0.002675119088962674,  -0.031726866960525513,   0.028613669797778130,   0.052236508578062057,   0.113318935036659241,  -0.058749049901962280,   0.054541323333978653,   0.012801205739378929,   0.152703285217285156,  -0.153122976422309875,
		  -0.024061812087893486,   0.088104479014873505,  -0.001473479554988444,   0.189222767949104309,   0.131869733333587646,   0.060047067701816559,   0.112139031291007996,   0.002661661012098193,  -0.017130997031927109,  -0.052856706082820892,   0.139452084898948669,   0.062293183058500290,   0.060585793107748032,  -0.004959530197083950,   0.185052797198295593,   0.040883865207433701,  -0.106581650674343109,   0.153828069567680359,   0.132410496473312378,   0.097592137753963470,   0.004023285117000341,  -0.033041749149560928,   0.078333467245101929,  -0.160478278994560242,   0.063326112926006317,   0.075771860778331757,   0.029764754697680473,   0.116112239658832550,   0.152185201644897461,  -0.129896268248558044,   0.123233251273632050,   0.002424981445074081,   0.183231368660926819,   0.195512920618057251,  -0.029115760698914528,  -0.132365301251411438,   0.025750547647476196,   0.035617709159851074,  -0.157616108655929565,   0.189493015408515930,  -0.018804110586643219,  -0.058664988726377487,   0.120002977550029755,   0.005570434499531984,  -0.053406145423650742,  -0.146136015653610229,   0.005495319142937660,   0.027002505958080292,
		   0.043035481125116348,  -0.051608931273221970,  -0.133071720600128174,   0.119622074067592621,   0.131746858358383179,   0.095771253108978271,   0.125821202993392944,   0.090535528957843781,  -0.181883350014686584,   0.034679230302572250,   0.013357944786548615,   0.128978490829467773,  -0.080972857773303986,   0.152139976620674133,  -0.099884018301963806,  -0.027243388816714287,  -0.093376114964485168,  -0.022380985319614410,   0.039099019020795822,  -0.021646136417984962,   0.153102368116378784,   0.165645688772201538,   0.028444401919841766,   0.112616673111915588,  -0.135848253965377808,   0.140300169587135315,   0.008532257750630379,   0.099011830985546112,   0.061706654727458954,   0.072103142738342285,  -0.119406141340732574,   0.152945473790168762,  -0.023647317662835121,  -0.033905286341905594,   0.187254399061203003,   0.063117042183876038,   0.086293898522853851,  -0.092991322278976440,  -0.038680627942085266,  -0.010109166614711285,  -0.108587808907032013,   0.029475543648004532,  -0.021960394456982613,   0.065945588052272797,  -0.004335090052336454,  -0.041374459862709045,  -0.083583272993564606,   0.046542931348085403,
		   0.030736664310097694,  -0.047731883823871613,  -0.091265797615051270,   0.033177562057971954,   0.120186075568199158,   0.078932181000709534,  -0.087943293154239655,   0.053851578384637833,  -0.088454604148864746,   0.024089975282549858,  -0.153161823749542236,   0.112058915197849274,  -0.005812857765704393,  -0.028743321076035500,  -0.016239237040281296,  -0.108414933085441589,  -0.109223119914531708,  -0.009222346358001232,  -0.002714463742449880,  -0.031012143939733505,   0.035027243196964264,  -0.127569675445556641,   0.147839531302452087,   0.076614521443843842,  -0.002002984983846545,   0.010777817107737064,   0.140563562512397766,  -0.132241651415824890,   0.067147798836231232,   0.130460694432258606,  -0.027546247467398643,   0.130515828728675842,  -0.149504229426383972,  -0.101673476397991180,   0.148441612720489502,   0.127416908740997314,  -0.002800193149596453,   0.019559202715754509,   0.079487249255180359,  -0.019352171570062637,  -0.056118708103895187,   0.070603221654891968,   0.110526405274868011,   0.054726809263229370,  -0.082315705716609955,   0.007037003990262747,   0.005400472786277533,  -0.140690594911575317,
		  -0.135549351572990417,   0.146457418799400330,   0.114788584411144257,   0.190648421645164490,  -0.024983085691928864,  -0.080368980765342712,   0.165717571973800659,   0.059555120766162872,  -0.080234266817569733,   0.027841832488775253,  -0.155900657176971436,  -0.010383958928287029,   0.132596939802169800,   0.021183751523494720,   0.076429769396781921,   0.077664174139499664,  -0.049473755061626434,   0.085666276514530182,   0.154493153095245361,   0.027046037837862968,  -0.084671422839164734,  -0.097768478095531464,  -0.012463135644793510,  -0.034655544906854630,  -0.148722305893898010,  -0.042999681085348129,  -0.015486748889088631,   0.010297165252268314,   0.084125623106956482,   0.055517960339784622,  -0.001485029817558825,   0.018162664026021957,  -0.001551303779706359,   0.071438163518905640,  -0.001360890804789960,  -0.068773649632930756,  -0.130112498998641968,  -0.014001997187733650,   0.128205582499504089,   0.120741598308086395,  -0.008444543927907944,  -0.072182781994342804,   0.090715698897838593,   0.010539171285927296,   0.096565544605255127,   0.138178274035453796,   0.131752982735633850,   0.150427281856536865,
		  -0.031249290332198143,  -0.129080832004547119,  -0.123803302645683289,  -0.054695725440979004,  -0.068084396421909332,  -0.174111232161521912,  -0.013459668494760990,   0.085912257432937622,  -0.102126955986022949,   0.076707728207111359,   0.148058354854583740,  -0.165849789977073669,  -0.131385773420333862,  -0.059799492359161377,  -0.117272622883319855,   0.134472578763961792,  -0.034306816756725311,  -0.054569918662309647,   0.093064658343791962,  -0.060514789074659348,  -0.075512632727622986,  -0.089501038193702698,   0.025463208556175232,  -0.015849115327000618,  -0.062253624200820923,  -0.012760827317833900,  -0.099076434969902039,  -0.086633674800395966,  -0.110890626907348633,  -0.060550350695848465,  -0.065847888588905334,  -0.160142064094543457,  -0.073176957666873932,   0.083554871380329132,   0.101046644151210785,   0.096628092229366302,  -0.008961745537817478,   0.117030151188373566,  -0.028779016807675362,  -0.005993683822453022,  -0.079105839133262634,   0.094990298151969910,  -0.130163043737411499,   0.070608161389827728,   0.102244503796100616,   0.017726918682456017,  -0.061564937233924866,   0.064135581254959106,
		  -0.116121187806129456,   0.152400463819503784,   0.078465022146701813,   0.081251025199890137,   0.031020116060972214,  -0.014244426973164082,   0.054418832063674927,   0.165486782789230347,   0.154007747769355774,   0.163920819759368896,   0.090056568384170532,   0.125076919794082642,  -0.079690486192703247,   0.003837083699181676,  -0.042205084115266800,  -0.073788866400718689,   0.137679696083068848,   0.030758641660213470,  -0.037133418023586273,   0.100031711161136627,  -0.013961951248347759,  -0.024590993300080299,   0.042356949299573898,   0.170927807688713074,  -0.118614338338375092,   0.112316481769084930,   0.044250182807445526,  -0.070142574608325958,   0.171545520424842834,   0.053086679428815842,  -0.033152449876070023,   0.037539433687925339,  -0.087848134338855743,  -0.017721939831972122,   0.164606198668479919,   0.082784466445446014,  -0.049788407981395721,   0.023472813889384270,  -0.001569750602357090,   0.103292658925056458,   0.153973206877708435,   0.112262852489948273,   0.044183142483234406,   0.153435319662094116,   0.148071333765983582,  -0.023559361696243286,   0.013133953325450420,  -0.010254160501062870,
		   0.060957163572311401,   0.074381649494171143,  -0.199179992079734802,  -0.066357977688312531,  -0.135549485683441162,   0.184679582715034485,   0.042315568774938583,   0.019644156098365784,  -0.191476598381996155,   0.111070737242698669,   0.058854792267084122,   0.150668516755104065,  -0.127937823534011841,  -0.087226673960685730,  -0.043637245893478394,  -0.040708672255277634,  -0.104559801518917084,   0.131288900971412659,  -0.100008368492126465,  -0.100718289613723755,  -0.099322356283664703,  -0.048311050981283188,   0.071906134486198425,   0.050413884222507477,  -0.102452404797077179,  -0.061214938759803772,   0.157051250338554382,   0.096343666315078735,  -0.060051847249269485,   0.123813509941101074,  -0.052947007119655609,   0.136357009410858154,  -0.023371731862425804,   0.081363938748836517,   0.134415894746780396,   0.045840002596378326,   0.041070368140935898,  -0.000534824328497052,   0.138230070471763611,  -0.023784777149558067,   0.001098547130823135,   0.031070789322257042,   0.127135396003723145,  -0.135617971420288086,   0.095790795981884003,   0.069926850497722626,   0.021725270897150040,   0.001176417805254459,
		  -0.069048561155796051,   0.136178955435752869,   0.005762428510934114,   0.002613402903079987,   0.025776328518986702,   0.099983230233192444,  -0.006485856138169765,   0.065071821212768555,   0.051960296928882599,   0.063792228698730469,  -0.152100175619125366,   0.022733625024557114,  -0.032089926302433014,   0.026263659819960594,  -0.010579736903309822,  -0.139541938900947571,   0.060578353703022003,  -0.134952083230018616,  -0.000371048925444484,  -0.086681887507438660,  -0.079836457967758179,   0.066436290740966797,  -0.074399918317794800,   0.006253677420318127,   0.054000534117221832,   0.018428299576044083,   0.137211799621582031,  -0.114588387310504913,  -0.043970577418804169,   0.093196123838424683,  -0.069998785853385925,   0.118598081171512604,   0.087714016437530518,   0.022683417424559593,   0.076950363814830780,  -0.089346498250961304,  -0.162416413426399231,  -0.088136598467826843,   0.161292627453804016,  -0.089677266776561737,   0.003894092515110970,   0.071776837110519409,   0.015202985145151615,   0.026195457205176353,  -0.078775785863399506,  -0.006277656648308039,   0.132827386260032654,  -0.055413533002138138,
		  -0.100354954600334167,   0.164575487375259399,  -0.003649603575468063,   0.181931748986244202,  -0.065348148345947266,  -0.084313593804836273,  -0.022707894444465637,  -0.068749777972698212,   0.026687711477279663,   0.110716506838798523,  -0.096745751798152924,  -0.146494820713996887,   0.097903326153755188,   0.093464411795139313,   0.104906328022480011,   0.029703412204980850,   0.104085795581340790,   0.124449655413627625,  -0.018576543778181076,  -0.141366988420486450,   0.120139427483081818,   0.000253397185588256,   0.096029140055179596,  -0.019717702642083168,  -0.041350033134222031,   0.019156720489263535,  -0.172372013330459595,   0.013245392590761185,   0.117069087922573090,  -0.115430161356925964,  -0.005996365565806627,   0.026549756526947021,   0.176042363047599792,   0.073361314833164215,   0.049831300973892212,  -0.141453042626380920,   0.104471147060394287,  -0.075108632445335388,   0.073291920125484467,   0.169612497091293335,   0.104497142136096954,  -0.066367574036121368,  -0.124205872416496277,   0.197562366724014282,   0.002777614165097475,  -0.007320856209844351,   0.146393641829490662,   0.168409615755081177,
		   0.022590760141611099,  -0.159479945898056030,  -0.023074343800544739,  -0.074645005166530609,  -0.012482203543186188,   0.103511601686477661,  -0.121317766606807709,  -0.005281392484903336,   0.104712024331092834,   0.126159310340881348,  -0.105178520083427429,  -0.030129941180348396,   0.029319232329726219,  -0.084878422319889069,  -0.100911021232604980,  -0.119322948157787323,  -0.010904983617365360,  -0.127464696764945984,   0.034551866352558136,   0.123265251517295837,  -0.075970232486724854,  -0.164787083864212036,  -0.059119336307048798,   0.058816719800233841,  -0.046181242913007736,  -0.159834071993827820,  -0.117632657289505005,  -0.069061353802680969,   0.118833251297473907,  -0.041710682213306427,   0.098130509257316589,  -0.028825853019952774,   0.117304816842079163,   0.062640994787216187,   0.019658282399177551,   0.026605403050780296,  -0.001997504848986864,  -0.001618339214473963,  -0.003637236542999744,  -0.135527595877647400,   0.005856968928128481,   0.065054506063461304,  -0.057263098657131195,  -0.132830128073692322,  -0.002008778275921941,   0.086927197873592377,  -0.087865501642227173,  -0.091377109289169312,
		  -0.164019063115119934,  -0.107604116201400757,   0.091741755604743958,   0.096039257943630219,  -0.090105839073657990,  -0.070093885064125061,  -0.190466776490211487,  -0.065731577575206757,   0.154235467314720154,  -0.185500144958496094,   0.101824671030044556,  -0.124338075518608093,  -0.035638350993394852,  -0.062678061425685883,  -0.122327037155628204,  -0.129636138677597046,   0.029292270541191101,  -0.151610940694808960,   0.038903687149286270,  -0.006124733947217464,  -0.111804157495498657,   0.132488876581192017,   0.022961063310503960,  -0.096749499440193176,   0.133950397372245789,  -0.048448100686073303,  -0.084427215158939362,  -0.058263674378395081,  -0.111962348222732544,   0.005718168336898088,   0.026523763313889503,  -0.032276511192321777,   0.014575387351214886,   0.090650290250778198,  -0.045495767146348953,  -0.147657990455627441,  -0.008285324089229107,   0.088731080293655396,  -0.028971718624234200,  -0.087145149707794189,  -0.119644291698932648,  -0.180771186947822571,   0.039897486567497253,   0.038412764668464661,  -0.067955277860164642,  -0.067218810319900513,  -0.179483786225318909,   0.027887742966413498,
		   0.064360961318016052,   0.070648759603500366,   0.164948299527168274,   0.107620827853679657,  -0.055459193885326385,   0.104744762182235718,  -0.017587402835488319,  -0.051862709224224091,   0.002904613036662340,  -0.031204544007778168,  -0.115925982594490051,  -0.122288048267364502,  -0.087109297513961792,   0.145353347063064575,  -0.060228686779737473,  -0.134652167558670044,  -0.010785207152366638,  -0.049258008599281311,  -0.102440252900123596,  -0.113624408841133118,   0.103240936994552612,  -0.075312137603759766,  -0.118008755147457123,   0.102001309394836426,  -0.004928790498524904,  -0.080484546720981598,  -0.093072280287742615,  -0.111722201108932495,   0.138613566756248474,  -0.058953851461410522,   0.170946151018142700,  -0.039293635636568069,  -0.075708307325839996,   0.112710483372211456,  -0.088586129248142242,   0.053319938480854034,   0.079908639192581177,  -0.028592308983206749,  -0.058342427015304565,   0.180689990520477295,  -0.026370279490947723,  -0.071581341326236725,   0.068171046674251556,   0.168150722980499268,  -0.038925521075725555,  -0.115637868642807007,   0.070433646440505981,   0.059967041015625000,
		   0.128429517149925232,   0.060972191393375397,  -0.083539336919784546,  -0.091476358473300934,  -0.090183936059474945,  -0.072969302535057068,  -0.021806512027978897,   0.110624127089977264,   0.055561255663633347,   0.031533949077129364,  -0.069305412471294403,   0.171706527471542358,  -0.076733112335205078,   0.045630216598510742,   0.104735203087329865,   0.059680763632059097,   0.089801199734210968,  -0.091116867959499359,   0.115830972790718079,   0.087605252861976624,   0.167179301381111145,  -0.108207382261753082,   0.088508017361164093,   0.139637693762779236,  -0.113734163343906403,   0.155707389116287231,   0.180677369236946106,  -0.066242747008800507,  -0.010249536484479904,   0.102915771305561066,  -0.071363516151905060,  -0.090201385319232941,   0.039548154920339584,   0.018105000257492065,   0.025120677426457405,   0.152211338281631470,  -0.074312739074230194,  -0.079061917960643768,  -0.002476138295605779,   0.000097704847576097,   0.082620657980442047,   0.101502835750579834,   0.079948574304580688,  -0.136939048767089844,  -0.007452093530446291,  -0.026703292503952980,  -0.082657523453235626,  -0.040982149541378021,
		  -0.046573027968406677,  -0.000661114696413279,   0.066134773194789886,  -0.103196077048778534,  -0.012426370754837990,  -0.115006074309349060,   0.004052583128213882,   0.060912050306797028,  -0.028253063559532166,   0.042100839316844940,   0.065614320337772369,   0.011745055206120014,  -0.077017351984977722,   0.089338608086109161,   0.119627192616462708,  -0.055149629712104797,   0.022304266691207886,  -0.028128992766141891,  -0.081367447972297668,   0.128874316811561584,   0.112437725067138672,   0.086372360587120056,  -0.078995667397975922,  -0.121387146413326263,  -0.052218936383724213,  -0.034579783678054810,   0.060324411839246750,   0.082587584853172302,  -0.098643109202384949,  -0.033064343035221100,  -0.045613981783390045,   0.060371007770299911,   0.050011083483695984,   0.144681632518768311,   0.172605663537979126,   0.134602710604667664,   0.101431116461753845,   0.032372459769248962,   0.108786419034004211,   0.166551932692527771,   0.045442570000886917,  -0.061177156865596771,   0.115678660571575165,  -0.045296773314476013,   0.189950570464134216,   0.033633407205343246,   0.188182741403579712,  -0.086265258491039276,
		  -0.015713600441813469,  -0.130874156951904297,  -0.038029201328754425,   0.100716680288314819,  -0.011405429802834988,  -0.048448160290718079,  -0.086507424712181091,   0.085076145827770233,  -0.021534925326704979,  -0.046579681336879730,  -0.028293935582041740,  -0.066450372338294983,   0.130179405212402344,   0.030146922916173935,  -0.178262203931808472,  -0.011905678547918797,   0.052898544818162918,  -0.054705534130334854,   0.077376291155815125,   0.073263466358184814,  -0.028268033638596535,  -0.072128631174564362,  -0.092875570058822632,   0.070512786507606506,   0.079280436038970947,   0.071887768805027008,   0.133574604988098145,  -0.088144756853580475,   0.090664178133010864,   0.015164108015596867,  -0.080655165016651154,  -0.008227327838540077,  -0.097166836261749268,   0.100360900163650513,  -0.164608091115951538,  -0.044677626341581345,   0.100045591592788696,  -0.098311960697174072,   0.115782789885997772,  -0.072168052196502686,  -0.132260590791702271,  -0.049667630344629288,   0.061927601695060730,  -0.128507822751998901,   0.000890955561771989,   0.048896390944719315,   0.067223690450191498,  -0.039962645620107651,
		   0.075172804296016693,   0.032483700662851334,  -0.020425124093890190,   0.003807694651186466,  -0.004556184168905020,   0.053917706012725830,   0.093123182654380798,  -0.095882624387741089,  -0.088477753102779388,  -0.109499745070934296,   0.033612534403800964,  -0.131895259022712708,  -0.037011735141277313,   0.080545559525489807,   0.113863140344619751,  -0.003377483226358891,  -0.060866404324769974,   0.056599635630846024,  -0.002503897761926055,  -0.008708786219358444,  -0.047158118337392807,   0.093941442668437958,   0.089756429195404053,  -0.057606995105743408,  -0.004059526138007641,   0.055884897708892822,   0.097355186939239502,  -0.093322329223155975,   0.124456167221069336,  -0.137379422783851624,  -0.078955888748168945,  -0.039158117026090622,   0.053558338433504105,  -0.030209956690669060,  -0.036855768412351608,   0.038994964212179184,  -0.020461987704038620,  -0.012167196720838547,  -0.120016932487487793,   0.025777045637369156,  -0.068036317825317383,   0.075476370751857758,  -0.134004473686218262,  -0.041504956781864166,   0.119275704026222229,  -0.134479001164436340,   0.044284481555223465,   0.050772808492183685,
		  -0.012107213959097862,  -0.087101228535175323,   0.148904204368591309,   0.084832437336444855,   0.165572360157966614,   0.027010357007384300,   0.110379964113235474,  -0.097858712077140808,   0.097772166132926941,  -0.000023236072593136,  -0.119717545807361603,  -0.058183882385492325,  -0.022127887234091759,  -0.057474866509437561,   0.073261775076389313,   0.113506741821765900,   0.109808564186096191,   0.154790490865707397,   0.002566998358815908,  -0.112501554191112518,   0.116623491048812866,   0.157957628369331360,  -0.043917849659919739,   0.120594874024391174,  -0.053922928869724274,   0.060194611549377441,   0.046409741044044495,   0.110348939895629883,   0.038485217839479446,   0.123291127383708954,   0.052200578153133392,  -0.011598174460232258,   0.059433467686176300,   0.078233689069747925,   0.116198025643825531,  -0.117546714842319489,   0.088056802749633789,  -0.080115161836147308,  -0.097229555249214172,   0.178822055459022522,   0.143245100975036621,   0.018683183938264847,   0.031705062836408615,   0.048409800976514816,   0.029886575415730476,  -0.089016221463680267,  -0.009964362718164921,   0.101661182940006256,
		   0.029938770458102226,  -0.017120836302638054,  -0.067898511886596680,   0.123678319156169891,   0.117579080164432526,   0.071179956197738647,   0.106020018458366394,  -0.092538975179195404,  -0.075381286442279816,  -0.108226373791694641,  -0.047031484544277191,  -0.049654986709356308,   0.118217088282108307,   0.158323869109153748,   0.018219152465462685,  -0.109723813831806183,   0.094082348048686981,  -0.112091191112995148,  -0.091733060777187347,   0.140908271074295044,   0.072534084320068359,   0.185388505458831787,   0.073597915470600128,  -0.091140598058700562,   0.156098082661628723,  -0.058770623058080673,  -0.177505537867546082,  -0.123071730136871338,   0.078414663672447205,  -0.098748460412025452,  -0.031247975304722786,   0.031367670744657516,  -0.016423670575022697,   0.107058584690093994,  -0.083451695740222931,  -0.102780140936374664,  -0.107363484799861908,   0.033542577177286148,  -0.086558893322944641,   0.017721297219395638,   0.008757404983043671,   0.067113064229488373,  -0.082953952252864838,   0.019923498854041100,  -0.081024847924709320,  -0.066406749188899994,  -0.098437123000621796,   0.147763192653656006,
		   0.180556952953338623,  -0.061237495392560959,   0.085199780762195587,   0.117296800017356873,   0.007213824428617954,   0.127490952610969543,   0.021793708205223083,   0.081233568489551544,   0.060544859617948532,   0.130943089723587036,  -0.010269952006638050,  -0.139430418610572815,  -0.009216574020683765,   0.044378202408552170,   0.016364317387342453,   0.104302123188972473,  -0.077200122177600861,  -0.126262813806533813,   0.088742837309837341,  -0.036928638815879822,   0.037378881126642227,  -0.011240055784583092,   0.065195076167583466,   0.121249631047248840,  -0.135134175419807434,   0.024744536727666855,   0.028750430792570114,   0.038158167153596878,  -0.137397378683090210,   0.065492145717144012,  -0.033792801201343536,  -0.104411810636520386,   0.028303828090429306,   0.078007094562053680,   0.082224830985069275,  -0.003155188867822289,   0.043568465858697891,  -0.055296599864959717,  -0.007774777710437775,  -0.021487802267074585,  -0.029889903962612152,  -0.050948560237884521,   0.009020578116178513,  -0.094446927309036255,   0.041456900537014008,   0.055449031293392181,   0.067668244242668152,  -0.091885484755039215,
		  -0.158817648887634277,  -0.087950341403484344,   0.010914236307144165,  -0.058453548699617386,   0.080804891884326935,  -0.010931540280580521,  -0.005328736267983913,   0.104519560933113098,  -0.038502126932144165,  -0.020557597279548645,   0.092865429818630219,  -0.104556977748870850,   0.101554259657859802,   0.150053232908248901,  -0.007777234073728323,  -0.098346583545207977,  -0.036305639892816544,  -0.060367412865161896,   0.126134082674980164,   0.024014957249164581,   0.167551353573799133,  -0.093622580170631409,  -0.126491263508796692,   0.051462940871715546,   0.060575950890779495,   0.174675256013870239,  -0.132507890462875366,  -0.008276314474642277,   0.021088907495141029,  -0.049679894000291824,  -0.044304117560386658,  -0.060281004756689072,   0.168209373950958252,   0.111018233001232147,   0.081657811999320984,  -0.085834309458732605,  -0.045134171843528748,   0.044123817235231400,  -0.132630899548530579,   0.013349724933505058,  -0.060364961624145508,  -0.025789881125092506,  -0.110906951129436493,   0.041853159666061401,  -0.062088582664728165,  -0.067736819386482239,   0.156280517578125000,   0.155485793948173523,
		  -0.054262701421976089,   0.174032971262931824,  -0.093226835131645203,   0.170845970511436462,  -0.079192452132701874,   0.094980381429195404,   0.006295657716691494,  -0.010516392998397350,  -0.100484102964401245,  -0.014697342179715633,   0.035069003701210022,   0.063602112233638763,   0.120715327560901642,   0.118930995464324951,   0.065960116684436798,   0.090477786958217621,   0.075047403573989868,   0.003874421119689941,   0.048104412853717804,   0.082953482866287231,   0.092832840979099274,   0.004863406531512737,   0.053500395268201828,   0.106294184923171997,  -0.095308668911457062,   0.026636216789484024,  -0.106108129024505615,   0.022971417754888535,  -0.042275745421648026,   0.170427531003952026,   0.094198122620582581,   0.011052616871893406,   0.103654496371746063,   0.035714786499738693,   0.117074854671955109,   0.011949403211474419,   0.098859854042530060,   0.042192567139863968,  -0.111927069723606110,  -0.012379326857626438,   0.143534570932388306,  -0.091323055326938629,  -0.097895532846450806,   0.005574245005846024,   0.055158060044050217,   0.062499135732650757,  -0.067741408944129944,   0.043777059763669968,
		  -0.149010837078094482,   0.068960197269916534,  -0.078088574111461639,   0.067368917167186737,   0.158566176891326904,  -0.124939754605293274,  -0.047646019607782364,  -0.084260284900665283,  -0.100367799401283264,  -0.039225671440362930,   0.110261484980583191,  -0.126124545931816101,  -0.025226930156350136,   0.113310456275939941,   0.128279253840446472,   0.098246850073337555,   0.116711765527725220,   0.082090765237808228,  -0.074728295207023621,   0.067045979201793671,  -0.026588333770632744,   0.179995730519294739,   0.083840467035770416,  -0.020437160506844521,  -0.052184816449880600,   0.070424072444438934,  -0.051443643867969513,   0.016512727364897728,   0.031694117933511734,  -0.096405677497386932,   0.002093825256451964,   0.054496768862009048,   0.167232573032379150,  -0.028129847720265388,   0.039495054632425308,   0.045342598110437393,  -0.009885314851999283,   0.045534282922744751,  -0.135326266288757324,  -0.066587589681148529,   0.091787770390510559,   0.048788104206323624,   0.135198295116424561,  -0.038151025772094727,   0.052819237112998962,  -0.010250424034893513,  -0.055468443781137466,  -0.082267843186855316,
		  -0.099488712847232819,   0.079542547464370728,   0.031456232070922852,   0.059918120503425598,   0.034608233720064163,  -0.023434529080986977,   0.104441709816455841,   0.076610416173934937,   0.108310543000698090,   0.057314764708280563,  -0.026593886315822601,  -0.121020287275314331,  -0.058418776839971542,   0.069259084761142731,   0.168018117547035217,  -0.050092004239559174,   0.124844402074813843,   0.126867145299911499,   0.150653839111328125,   0.115911476314067841,   0.152479842305183411,   0.015559153631329536,  -0.097730152308940887,  -0.104502402245998383,   0.026422264054417610,   0.049872167408466339,  -0.043778587132692337,   0.110252439975738525,   0.166611313819885254,   0.123348899185657501,   0.079402789473533630,   0.077479176223278046,   0.007519797887653112,  -0.046430695801973343,   0.047568649053573608,   0.013415556401014328,  -0.021000415086746216,  -0.014896129257977009,  -0.104306392371654510,   0.124995581805706024,   0.014080710709095001,  -0.188040688633918762,   0.020323237404227257,   0.128424122929573059,  -0.067744977772235870,   0.045201789587736130,   0.128239914774894714,   0.043083608150482178,
		  -0.058756362646818161,  -0.062583185732364655,  -0.062101196497678757,   0.175570860505104065,   0.162200614809989929,  -0.187177792191505432,   0.082654818892478943,   0.067242331802845001,  -0.055717866867780685,  -0.012685503810644150,  -0.079955011606216431,  -0.003429951379075646,   0.072742760181427002,   0.163212120532989502,   0.140668526291847229,  -0.145573765039443970,   0.071820817887783051,   0.175594836473464966,   0.085066631436347961,  -0.127830341458320618,  -0.077881649136543274,   0.043065454810857773,  -0.036360371857881546,   0.057392016053199768,   0.043020796030759811,   0.022505775094032288,  -0.137165039777755737,   0.003199290717020631,  -0.028960350900888443,   0.136466830968856812,   0.098767518997192383,   0.053316220641136169,   0.150348439812660217,   0.053093481808900833,   0.012067598290741444,  -0.076520457863807678,  -0.092277534306049347,   0.018363712355494499,  -0.034338705241680145,   0.042883828282356262,   0.074501499533653259,  -0.174008622765541077,  -0.097362563014030457,   0.089128732681274414,   0.137979686260223389,  -0.032501671463251114,   0.131797358393669128,   0.047697339206933975,
		   0.035790912806987762,  -0.031803257763385773,  -0.134954556822776794,   0.118158794939517975,  -0.054803367704153061,   0.069560348987579346,   0.116449475288391113,  -0.083184324204921722,  -0.159407898783683777,   0.146340996026992798,  -0.151007726788520813,  -0.015220516361296177,   0.092037878930568695,  -0.006269811186939478,  -0.024434830993413925,   0.013156382367014885,  -0.085224367678165436,   0.010085416957736015,   0.074483752250671387,  -0.127991348505020142,   0.127088904380798340,   0.070299230515956879,   0.203097015619277954,   0.156998828053474426,  -0.107705302536487579,   0.150076329708099365,   0.079136721789836884,   0.012684185989201069,  -0.090462960302829742,   0.067052923142910004,   0.020306380465626717,   0.100426055490970612,  -0.049991209059953690,   0.052192609757184982,   0.080229699611663818,  -0.024537332355976105,  -0.013578329235315323,   0.044142425060272217,  -0.044069338589906693,   0.105370469391345978,  -0.020819962024688721,   0.017159717157483101,  -0.046933405101299286,  -0.111189179122447968,  -0.019978400319814682,   0.134039759635925293,   0.104395523667335510,   0.179164364933967590,
		   0.023451888933777809,   0.105188675224781036,  -0.178642630577087402,   0.048080254346132278,   0.058905795216560364,   0.134412467479705811,   0.186316221952438354,  -0.028956947848200798,  -0.171862840652465820,  -0.027005638927221298,  -0.047459855675697327,   0.083392225205898285,  -0.021188534796237946,   0.142301186919212341,   0.009265929460525513,  -0.022916369140148163,  -0.129630893468856812,  -0.038875285536050797,  -0.109885752201080322,   0.053634874522686005,  -0.091951228678226471,   0.011161853559315205,  -0.003858848707750440,   0.074514374136924744,   0.048946086317300797,  -0.058319941163063049,   0.059917088598012924,  -0.092772454023361206,   0.099567569792270660,   0.077880948781967163,  -0.076246008276939392,  -0.039766222238540649,  -0.093989983201026917,  -0.058617781847715378,   0.092123478651046753,  -0.073204800486564636,  -0.097070768475532532,   0.052099570631980896,   0.174145117402076721,  -0.144616007804870605,  -0.131718009710311890,   0.126363709568977356,   0.131597086787223816,  -0.031680952757596970,  -0.059286255389451981,  -0.059017222374677658,   0.088913515210151672,   0.040008515119552612,
		   0.000847409071866423,   0.031702745705842972,  -0.008713041432201862,  -0.011832538060843945,  -0.003115277737379074,  -0.145624712109565735,  -0.168132096529006958,  -0.157965660095214844,   0.046220254153013229,  -0.113409422338008881,  -0.120755068957805634,  -0.025332873687148094,   0.084084019064903259,   0.048930097371339798,  -0.059483010321855545,   0.127551332116127014,  -0.127928301692008972,   0.106300845742225647,  -0.032614067196846008,   0.055226366966962814,   0.082512855529785156,  -0.080875717103481293,  -0.004219670314341784,  -0.173854216933250427,   0.042753614485263824,   0.096349462866783142,  -0.150588646531105042,  -0.100159652531147003,  -0.073943011462688446,  -0.139126926660537720,  -0.095902033150196075,  -0.050668939948081970,   0.117627345025539398,  -0.065861664712429047,  -0.096738487482070923,  -0.051474794745445251,  -0.122210867702960968,   0.063112303614616394,  -0.037447184324264526,   0.085282340645790100,  -0.073252387344837189,   0.024791248142719269,  -0.183024093508720398,   0.110847003757953644,  -0.011542009189724922,   0.015406317077577114,  -0.060777764767408371,  -0.060865897685289383,
		   0.029298903420567513,   0.058798737823963165,  -0.004204698838293552,  -0.098013952374458313,  -0.049718525260686874,   0.023982562124729156,   0.145442947745323181,  -0.036896653473377228,  -0.006620343308895826,   0.008267990313470364,  -0.109206423163414001,   0.057473760098218918,  -0.043990168720483780,   0.162376880645751953,   0.075712151825428009,  -0.106874123215675354,   0.074131518602371216,  -0.043758150190114975,   0.081136181950569153,  -0.176452800631523132,   0.129542142152786255,  -0.083103239536285400,   0.154992565512657166,   0.179362103343009949,   0.014350671321153641,   0.126573935151100159,   0.119142152369022369,   0.089153349399566650,  -0.033796023577451706,   0.198384523391723633,   0.005459884647279978,   0.071005977690219879,   0.158111631870269775,  -0.068393617868423462,   0.103024259209632874,   0.122890226542949677,  -0.024122085422277451,   0.088642887771129608,  -0.042601197957992554,  -0.106117263436317444,   0.151332646608352661,   0.184005305171012878,   0.069505497813224792,  -0.003297823714092374,   0.164259836077690125,  -0.041991148144006729,   0.147724330425262451,   0.166688069701194763,
		  -0.030865857377648354,   0.065807864069938660,   0.172167316079139709,  -0.086763836443424225,  -0.087988905608654022,  -0.052083741873502731,  -0.078674651682376862,  -0.034996517002582550,   0.136192426085472107,   0.060746595263481140,  -0.081618018448352814,  -0.083585254848003387,   0.062590159475803375,  -0.087618574500083923,   0.122439235448837280,  -0.117163918912410736,   0.080136053264141083,   0.148325145244598389,  -0.022788034752011299,   0.088897943496704102,   0.014131646603345871,   0.112757034599781036,   0.100201398134231567,  -0.137471407651901245,   0.151905089616775513,   0.171579301357269287,   0.067380264401435852,   0.092197008430957794,   0.104729384183883667,   0.073513872921466827,  -0.054694652557373047,  -0.003173931268975139,   0.064143411815166473,  -0.071243003010749817,   0.163374096155166626,   0.106543473899364471,  -0.050521630793809891,  -0.055792983621358871,   0.017854634672403336,   0.196186214685440063,   0.016273008659482002,  -0.075041130185127258,  -0.126437589526176453,   0.080014482140541077,  -0.045863654464483261,  -0.018707022070884705,  -0.029200684279203415,  -0.084656320512294769,
		  -0.004193664994090796,  -0.033420421183109283,   0.055260725319385529,  -0.066975392401218414,  -0.013783404603600502,   0.081404902040958405,   0.005045595578849316,  -0.064630702137947083,   0.091129958629608154,   0.002206970704719424,   0.023538837209343910,   0.050513684749603271,  -0.039416234940290451,  -0.109249196946620941,  -0.062924154102802277,  -0.130825638771057129,   0.104672916233539581,   0.096855171024799347,  -0.120816990733146667,   0.080595195293426514,  -0.021058315411210060,  -0.028256976976990700,   0.148309201002120972,   0.118816085159778595,  -0.024721357971429825,   0.166531652212142944,   0.160618573427200317,  -0.117872484028339386,  -0.104602694511413574,   0.159178256988525391,   0.005967935081571341,   0.062128644436597824,  -0.082441218197345734,  -0.094991095364093781,   0.111358247697353363,  -0.004597725346684456,   0.110868752002716064,  -0.134353816509246826,   0.141480937600135803,  -0.130473434925079346,   0.073839113116264343,   0.070970073342323303,  -0.062231957912445068,  -0.030190505087375641,  -0.060220830142498016,  -0.120352327823638916,  -0.057796951383352280,  -0.082508176565170288,
		  -0.167377024888992310,   0.027338134124875069,  -0.038144037127494812,  -0.112825207412242889,  -0.043009437620639801,  -0.092029102146625519,  -0.016549995169043541,   0.041709836572408676,   0.168153136968612671,   0.077477119863033295,   0.042103189975023270,   0.106657035648822784,   0.143777295947074890,   0.013420603238046169,  -0.026373831555247307,   0.104090854525566101,   0.015935730189085007,  -0.065277174115180969,   0.115747086703777313,  -0.081313557922840118,  -0.045980706810951233,   0.101937770843505859,  -0.050930481404066086,  -0.029714321717619896,   0.058986295014619827,   0.057711161673069000,  -0.020270645618438721,  -0.027074113488197327,  -0.051884967833757401,  -0.050967916846275330,  -0.048385407775640488,  -0.106757290661334991,   0.013406234793365002,   0.064721964299678802,   0.052550401538610458,   0.043243758380413055,  -0.024384420365095139,  -0.028843730688095093,  -0.002064639469608665,   0.164817184209823608,   0.075051955878734589,  -0.059749651700258255,  -0.118866585195064545,   0.108321018517017365,   0.045927740633487701,  -0.097316972911357880,   0.062273457646369934,  -0.033483710139989853,
		  -0.160970225930213928,  -0.158086404204368591,   0.025414498522877693,  -0.038460154086351395,   0.059709306806325912,   0.007289425469934940,   0.084778122603893280,  -0.016065835952758789,   0.147713392972946167,   0.002527845092117786,  -0.037647578865289688,   0.021768419072031975,  -0.011881888844072819,   0.052971091121435165,  -0.073132425546646118,  -0.038442291319370270,  -0.058577746152877808,   0.094137825071811676,  -0.009734725579619408,  -0.139263361692428589,  -0.162188872694969177,   0.025992795825004578,  -0.186487093567848206,   0.075632743537425995,   0.017281113192439079,  -0.008000165224075317,   0.043937806040048599,  -0.060864310711622238,   0.062669686973094940,  -0.059997051954269409,  -0.034965597093105316,  -0.119887940585613251,   0.022964226081967354,  -0.058691769838333130,  -0.137779206037521362,  -0.058848638087511063,   0.071050629019737244,   0.131869584321975708,  -0.101083569228649139,  -0.044825144112110138,  -0.156198814511299133,  -0.010850234888494015,  -0.023744229227304459,  -0.073056548833847046,  -0.081911914050579071,   0.005449410527944565,  -0.109770156443119049,   0.123387664556503296,
		  -0.138477981090545654,   0.091959580779075623,  -0.091879747807979584,   0.143182754516601562,   0.010257311165332794,   0.106930352747440338,   0.014039428904652596,  -0.076857507228851318,   0.085334427654743195,   0.055684600025415421,   0.017649766057729721,  -0.068759448826313019,  -0.087058261036872864,   0.087420195341110229,  -0.045285303145647049,   0.046114075928926468,   0.004828961100429296,   0.063878171145915985,   0.013771968893706799,  -0.008172614499926567,   0.066706992685794830,  -0.066264465451240540,  -0.071846663951873779,  -0.039092116057872772,   0.024290215224027634,  -0.073848567903041840,   0.042210787534713745,   0.046099230647087097,  -0.034587748348712921,  -0.016803063452243805,   0.164147526025772095,  -0.052210055291652679,   0.123496413230895996,  -0.047873020172119141,   0.042642068117856979,   0.157247185707092285,   0.118476763367652893,   0.089698679745197296,  -0.082094416022300720,   0.095560535788536072,   0.023478329181671143,  -0.079288922250270844,  -0.123544767498970032,  -0.030116993933916092,   0.159932434558868408,  -0.132297158241271973,   0.114902578294277191,   0.056172661483287811,
	};
	static const double bias03[]=
	{
		   0.117992058396339417,  -0.021675340831279755,   0.081426568329334259,   0.144085705280303955,   0.073787607252597809,  -0.056295040994882584,  -0.143345355987548828,  -0.163472712039947510,  -0.048303071409463882,  -0.019421707838773727,  -0.086498595774173737,  -0.011546599678695202,   0.116687148809432983,   0.038041420280933380,  -0.024983735755085945,   0.123415492475032806,   0.084915086627006531,   0.109289430081844330,   0.174260228872299194,   0.010355518199503422,  -0.037407785654067993,   0.170528978109359741,   0.002562857232987881,   0.106780156493186951,  -0.145391047000885010,  -0.138135671615600586,  -0.100766077637672424,  -0.066095001995563507,   0.111217468976974487,  -0.139868333935737610,  -0.078394822776317596,   0.018147381022572517,   0.050694938749074936,   0.018498230725526810,   0.142205238342285156,  -0.046493671834468842,  -0.048338353633880615,  -0.060520585626363754,  -0.006827924400568008,   0.101972058415412903,   0.051784161478281021,  -0.142036110162734985,  -0.028962479904294014,  -0.031460519880056381,   0.115054324269294739,  -0.015350323170423508,  -0.174021601676940918,  -0.015958715230226517,
	};
	static const double weight04[]=
	{
		  -0.116490170359611511,  -0.075814835727214813,   0.123199604451656342,   0.036562323570251465,   0.119824267923831940,   0.029895056039094925,   0.123145602643489838,  -0.049231890588998795,  -0.019608728587627411,   0.056330244988203049,   0.072515070438385010,   0.126994654536247253,  -0.041901033371686935,   0.158658042550086975,  -0.151265323162078857,   0.025092508643865585,   0.152497276663780212,  -0.062683530151844025,   0.124414332211017609,   0.092306785285472870,   0.024376813322305679,  -0.147817745804786682,  -0.001067185308784246,   0.127542540431022644,   0.038073096424341202,   0.084590405225753784,  -0.088489413261413574,  -0.156651303172111511,   0.088838554918766022,  -0.153891369700431824,   0.014792335219681263,   0.179740831255912781,   0.050496425479650497,  -0.052223023027181625,   0.007098650559782982,   0.078903116285800934,  -0.032752323895692825,   0.046020671725273132,   0.037956058979034424,   0.092540293931961060,  -0.141648203134536743,   0.055935498327016830,   0.120549008250236511,   0.111882746219635010,  -0.011940369382500648,   0.130961179733276367,  -0.004800870548933744,   0.120775520801544189,
		  -0.032667413353919983,   0.193225026130676270,  -0.097303345799446106,  -0.005752740427851677,   0.033859807997941971,  -0.186696708202362061,   0.120994009077548981,   0.022744258865714073,   0.109234668314456940,   0.099114827811717987,  -0.009772947058081627,   0.167108207941055298,  -0.128521129488945007,   0.155073374509811401,  -0.015868233516812325,   0.186609387397766113,   0.137342661619186401,  -0.017719866707921028,   0.052552796900272369,   0.080500446259975433,   0.141995683312416077,   0.062236625701189041,  -0.034019231796264648,   0.164845630526542664,  -0.054636053740978241,   0.019836004823446274,  -0.022876460105180740,   0.035150215029716492,   0.040661089122295380,  -0.068824939429759979,  -0.060356143862009048,  -0.030170982703566551,  -0.021709363907575607,  -0.059128660708665848,   0.136237844824790955,   0.112023249268531799,  -0.091932341456413269,   0.088217973709106445,   0.049697484821081161,  -0.029078379273414612,  -0.133282184600830078,  -0.056889764964580536,  -0.033980596810579300,  -0.028873870149254799,   0.021869303658604622,   0.084193013608455658,  -0.121767275035381317,   0.152842238545417786,
		   0.059482548385858536,   0.169650182127952576,  -0.041027095168828964,   0.041559815406799316,   0.064709171652793884,  -0.114185608923435211,  -0.128648489713668823,  -0.054753504693508148,  -0.021475626155734062,   0.128430768847465515,  -0.170469924807548523,   0.024935435503721237,   0.001519494340755045,  -0.083692528307437897,  -0.140248388051986694,   0.026697585359215736,   0.146479815244674683,  -0.153269127011299133,   0.189512491226196289,  -0.078823916614055634,   0.149749785661697388,  -0.074064493179321289,  -0.088910557329654694,   0.161944970488548279,  -0.010018184781074524,  -0.088758751749992371,   0.049652166664600372,  -0.046394709497690201,   0.139224365353584290,  -0.045691229403018951,   0.104288741946220398,   0.080990038812160492,   0.098496437072753906,   0.009350318461656570,   0.176945090293884277,   0.114670261740684509,  -0.007246469147503376,   0.144916012883186340,   0.184134557843208313,  -0.081248916685581207,  -0.161082684993743896,  -0.062748983502388000,   0.062948159873485565,  -0.089229069650173187,   0.068906418979167938,   0.040621910244226456,  -0.078003950417041779,   0.037701785564422607,
		  -0.102575831115245819,   0.049883875995874405,   0.068904176354408264,  -0.055360104888677597,   0.084608897566795349,   0.124072775244712830,   0.104785054922103882,  -0.011353758163750172,  -0.110556237399578094,  -0.033253289759159088,  -0.064012221992015839,  -0.023910978808999062,   0.121656857430934906,   0.105249002575874329,   0.154695183038711548,  -0.003052570391446352,  -0.071690008044242859,   0.150803610682487488,  -0.124365478754043579,  -0.084069885313510895,   0.069622382521629333,   0.157035171985626221,  -0.048130143433809280,   0.011416711844503880,  -0.165334776043891907,  -0.050329789519309998,   0.024938696995377541,   0.060176704078912735,   0.075020626187324524,   0.081600017845630646,   0.043408345431089401,   0.134035855531692505,   0.090778499841690063,  -0.090946882963180542,   0.062581926584243774,   0.098631992936134338,  -0.105514928698539734,   0.141428217291831970,   0.162189766764640808,  -0.036992460489273071,   0.139244720339775085,  -0.159237667918205261,   0.158949717879295349,   0.156312465667724609,  -0.079491831362247467,   0.061799649149179459,  -0.107124686241149902,  -0.121086053550243378,
		   0.072813108563423157,  -0.077213868498802185,   0.004198337905108929,   0.037088114768266678,   0.020454468205571175,  -0.044473685324192047,   0.090910144150257111,  -0.121532320976257324,   0.035343848168849945,   0.055990464985370636,  -0.148868128657341003,   0.104155905544757843,  -0.008255879394710064,  -0.088397726416587830,   0.156919702887535095,   0.028325546532869339,   0.044665422290563583,  -0.035071566700935364,   0.119297519326210022,   0.097634688019752502,   0.083284050226211548,   0.200998768210411072,   0.165063753724098206,   0.121690399944782257,   0.088741131126880646,   0.020344821736216545,  -0.146792367100715637,   0.164279565215110779,  -0.075452670454978943,   0.060685973614454269,   0.086022853851318359,  -0.115475393831729889,   0.107776634395122528,  -0.041832581162452698,  -0.087799586355686188,   0.055531993508338928,  -0.079400554299354553,  -0.146164417266845703,  -0.159786805510520935,   0.075372569262981415,   0.009102120064198971,   0.054043937474489212,   0.069494999945163727,  -0.058350678533315659,   0.131988868117332458,   0.060036860406398773,  -0.073963418602943420,   0.080694787204265594,
		   0.027020594105124474,  -0.112843811511993408,   0.096964545547962189,   0.149049833416938782,  -0.003802043618634343,   0.085195012390613556,   0.041228707879781723,  -0.047650400549173355,  -0.040508724749088287,   0.054105777293443680,  -0.061921369284391403,  -0.004617616534233093,   0.172271236777305603,   0.119730919599533081,   0.171920046210289001,  -0.090337574481964111,  -0.015457617118954659,   0.153068363666534424,   0.121517002582550049,   0.003500704653561115,   0.080975621938705444,   0.173100292682647705,   0.132420286536216736,  -0.108376778662204742,   0.043379459530115128,  -0.167388141155242920,  -0.107630044221878052,   0.172049120068550110,   0.074891619384288788,  -0.120959654450416565,   0.046221997588872910,   0.038779489696025848,  -0.123374767601490021,  -0.041426885873079300,   0.002483336487784982,   0.069189615547657013,   0.113210327923297882,  -0.122404143214225769,  -0.019429123029112816,   0.072381086647510529,  -0.039893504232168198,   0.089729189872741699,  -0.022405939176678658,   0.087887741625308990,  -0.081665419042110443,   0.027325330302119255,  -0.064055204391479492,  -0.036418743431568146,
		  -0.112699031829833984,  -0.052253447473049164,  -0.034710455685853958,   0.014833118766546249,   0.066943041980266571,   0.086434245109558105,   0.133520200848579407,  -0.137355729937553406,  -0.127168387174606323,   0.085095271468162537,  -0.087635099887847900,  -0.090040072798728943,  -0.022643793374300003,  -0.073360025882720947,   0.020209280773997307,  -0.110501028597354889,  -0.058908145874738693,   0.038596738129854202,  -0.107714608311653137,   0.133903324604034424,   0.029938001185655594,  -0.047132976353168488,  -0.128071591258049011,  -0.074555858969688416,   0.122969239950180054,   0.036831546574831009,   0.125972971320152283,   0.109423957765102386,   0.072441220283508301,   0.089803822338581085,   0.095898419618606567,  -0.077609360218048096,  -0.075755231082439423,   0.135581716895103455,  -0.108129590749740601,  -0.094418130815029144,  -0.060709998011589050,   0.030179092660546303,   0.130904451012611389,   0.119439333677291870,  -0.009265373460948467,   0.074850454926490784,  -0.116975940763950348,   0.014501593075692654,  -0.011753107421100140,  -0.022534083575010300,   0.002897704718634486,  -0.130744650959968567,
		   0.091147258877754211,   0.059071261435747147,  -0.000823101669084281,  -0.072626866400241852,  -0.020194344222545624,   0.120831497013568878,   0.091965109109878540,   0.011555108241736889,   0.023952694609761238,   0.018243687227368355,  -0.081648908555507660,  -0.134059578180313110,  -0.049950305372476578,  -0.136740684509277344,   0.000117549643618986,  -0.144318357110023499,   0.144815221428871155,   0.108062617480754852,   0.051599066704511642,   0.105920992791652679,  -0.006483493372797966,   0.137901857495307922,  -0.054301939904689789,   0.048283636569976807,  -0.134890973567962646,   0.038368485867977142,  -0.153717726469039917,   0.000599048915319145,  -0.131315767765045166,   0.113699004054069519,   0.133317038416862488,   0.039090137928724289,  -0.044030793011188507,  -0.057010002434253693,   0.112402319908142090,   0.043618500232696533,   0.076289802789688110,   0.090928278863430023,   0.088427104055881500,   0.028573913499712944,   0.039804570376873016,   0.076377183198928833,   0.159384235739707947,   0.061941225081682205,   0.081018216907978058,   0.036187153309583664,   0.119649998843669891,  -0.114846043288707733,
		  -0.015695815905928612,   0.033299610018730164,   0.092924125492572784,  -0.048223033547401428,   0.098506674170494080,   0.015754079446196556,   0.015864251181483269,  -0.078862108290195465,  -0.115248948335647583,   0.168432742357254028,  -0.086883127689361572,  -0.092373773455619812,  -0.029619541019201279,   0.115431636571884155,   0.049120128154754639,   0.052264221012592316,   0.044428270310163498,  -0.098159179091453552,   0.041669934988021851,  -0.032676149159669876,   0.004273733124136925,  -0.126927867531776428,  -0.089175790548324585,   0.186066642403602600,  -0.137989386916160583,  -0.084198892116546631,   0.084794677793979645,   0.049989983439445496,   0.089417122304439545,  -0.146548762917518616,   0.027735894545912743,   0.156202659010887146,   0.149154916405677795,  -0.048277176916599274,  -0.057870287448167801,  -0.024491410702466965,  -0.084761694073677063,   0.112811625003814697,   0.127435386180877686,   0.092761926352977753,  -0.018836855888366699,  -0.010778582654893398,   0.153275564312934875,  -0.026429394260048866,  -0.120428614318370819,   0.142559170722961426,   0.095740430057048798,   0.098419241607189178,
		   0.093699865043163300,  -0.057308487594127655,  -0.113280273973941803,   0.128796979784965515,   0.038249872624874115,  -0.101528018712997437,  -0.014005501754581928,  -0.005837518256157637,  -0.010771778412163258,  -0.031599882990121841,  -0.027619084343314171,  -0.126268550753593445,   0.017221620306372643,   0.074103541672229767,   0.136359661817550659,  -0.059929534792900085,  -0.110419243574142456,  -0.007590251043438911,   0.062371358275413513,   0.125407561659812927,  -0.057340241968631744,  -0.021930208429694176,  -0.039104241877794266,  -0.063240714371204376,  -0.121565543115139008,   0.018127834424376488,   0.068877167999744415,   0.121840827167034149,  -0.115599863231182098,  -0.059169746935367584,  -0.037653379142284393,  -0.012420400045812130,  -0.129802271723747253,   0.034680113196372986,  -0.087988466024398804,  -0.028145104646682739,   0.018977645784616470,  -0.101769611239433289,  -0.021443000063300133,  -0.029369959607720375,  -0.065927371382713318,   0.117689698934555054,   0.008101751096546650,   0.056931070983409882,  -0.016017304733395576,  -0.129915565252304077,  -0.049591612070798874,  -0.127822592854499817,
		   0.018718745559453964,  -0.083401314914226532,   0.014033253304660320,  -0.035199452191591263,  -0.061525940895080566,   0.056274872273206711,   0.080831110477447510,   0.053398139774799347,   0.035297751426696777,   0.127916157245635986,  -0.111534148454666138,  -0.063053682446479797,   0.037120889872312546,  -0.130455866456031799,   0.144003450870513916,   0.045259598642587662,   0.009408480487763882,  -0.000264301081188023,   0.100964546203613281,   0.037587918341159821,   0.080499753355979919,  -0.034718003123998642,  -0.003240532707422972,  -0.022591408342123032,   0.135996848344802856,   0.071941442787647247,  -0.117998763918876648,   0.128353789448738098,   0.065401040017604828,   0.018535040318965912,   0.078300118446350098,   0.107631340622901917,   0.009191724471747875,   0.018426612019538879,  -0.008736078627407551,   0.077560007572174072,   0.103014685213565826,   0.117650404572486877,  -0.114762030541896820,   0.088819950819015503,   0.013051425106823444,   0.119859285652637482,   0.107703156769275665,   0.045893106609582901,  -0.014407273381948471,   0.111572951078414917,  -0.016945509240031242,  -0.135214343667030334,
		   0.046814601868391037,  -0.088745638728141785,   0.085434041917324066,   0.028050465509295464,   0.067775413393974304,   0.036916024982929230,   0.114122264087200165,   0.032865781337022781,   0.059141878038644791,  -0.048287786543369293,   0.117212809622287750,   0.066351525485515594,   0.094968624413013458,   0.128614023327827454,  -0.044480007141828537,  -0.036338400095701218,   0.011745489202439785,  -0.077115230262279510,   0.064110659062862396,  -0.163722842931747437,   0.125119522213935852,   0.004065401852130890,   0.018641447648406029,   0.029794495552778244,  -0.110711187124252319,   0.153810679912567139,   0.078833587467670441,  -0.020422045141458511,   0.034034516662359238,   0.001645830925554037,  -0.000615169061347842,   0.109461106359958649,  -0.108258858323097229,  -0.141426548361778259,  -0.026410663500428200,   0.035520255565643311,   0.050311602652072906,   0.118929222226142883,   0.062843285501003265,  -0.123012363910675049,  -0.172079682350158691,   0.155131161212921143,  -0.067997030913829803,  -0.002318186452612281,   0.018700301647186279,  -0.063571982085704803,   0.094921305775642395,   0.152086704969406128,
		  -0.103946298360824585,  -0.153744578361511230,  -0.013906908221542835,   0.140009075403213501,   0.038731221109628677,   0.153983712196350098,  -0.047244779765605927,  -0.058535389602184296,  -0.039268635213375092,  -0.114800639450550079,  -0.152916535735130310,   0.014752100221812725,   0.091781087219715118,  -0.120530843734741211,   0.011738582514226437,   0.107800304889678955,   0.126018315553665161,   0.019374180585145950,  -0.076821692287921906,  -0.129951372742652893,   0.161556914448738098,   0.169984713196754456,   0.159378454089164734,  -0.134926557540893555,  -0.119697339832782745,   0.017271846532821655,   0.087364092469215393,  -0.075334094464778900,   0.064785294234752655,   0.109185926616191864,   0.095594823360443115,  -0.029113408178091049,  -0.115481182932853699,   0.175729751586914062,  -0.011604594066739082,  -0.024956738576292992,  -0.108164049685001373,  -0.001632318715564907,  -0.108879387378692627,   0.073841549456119537,  -0.081205964088439941,   0.011188737116754055,   0.061036568135023117,  -0.163054898381233215,  -0.054062515497207642,  -0.108182393014431000,  -0.073489956557750702,  -0.148785039782524109,
		  -0.133266702294349670,   0.084438763558864594,   0.051568605005741119,  -0.054991018027067184,  -0.073270112276077271,   0.118025809526443481,   0.008283281698822975,   0.120151810348033905,  -0.015993297100067139,  -0.048022814095020294,   0.062832899391651154,   0.116068653762340546,   0.094449765980243683,  -0.082937963306903839,   0.144220933318138123,  -0.138495877385139465,  -0.024042803794145584,   0.118147037923336029,   0.017367342486977577,   0.086243078112602234,  -0.043651957064867020,  -0.050299938768148422,  -0.096894860267639160,  -0.029437700286507607,  -0.098921194672584534,  -0.129801034927368164,  -0.002949489513412118,  -0.116103745996952057,   0.077660113573074341,  -0.086419805884361267,  -0.088165588676929474,  -0.086031600832939148,  -0.076527319848537445,   0.030641091987490654,  -0.055413108319044113,   0.007569057866930962,   0.083086349070072174,  -0.125086709856987000,  -0.028619851917028427,  -0.098530061542987823,   0.003365292679518461,   0.048351123929023743,  -0.100781418383121490,   0.123161263763904572,  -0.007506756577640772,  -0.112921185791492462,  -0.139865145087242126,  -0.094653539359569550,
		  -0.121870361268520355,  -0.036546602845191956,   0.026900989934802055,  -0.034140646457672119,   0.183370321989059448,  -0.056372731924057007,  -0.014997038058936596,  -0.160957232117652893,   0.023829674348235130,   0.161206141114234924,  -0.101339042186737061,   0.021903902292251587,  -0.064001969993114471,  -0.087088979780673981,   0.040335141122341156,   0.185619369149208069,   0.087315328419208527,  -0.004879951477050781,   0.128868222236633301,  -0.017313985154032707,   0.110722459852695465,  -0.038544692099094391,  -0.013628238812088966,   0.169157788157463074,   0.059626150876283646,  -0.059005700051784515,   0.024667844176292419,   0.058723874390125275,   0.023907648399472237,  -0.161761835217475891,   0.068364880979061127,   0.133491665124893188,  -0.038916125893592834,  -0.029642732813954353,   0.069663219153881073,  -0.036997664719820023,   0.148695722222328186,   0.134090870618820190,   0.186047524213790894,  -0.093770422041416168,  -0.104541115462779999,  -0.024230726063251495,   0.109030075371265411,   0.148184359073638916,  -0.115263901650905609,   0.079600796103477478,  -0.046833332628011703,  -0.046918679028749466,
		   0.094714231789112091,   0.084159336984157562,   0.174029648303985596,   0.084504254162311554,   0.010411582887172699,  -0.051549110561609268,   0.092220813035964966,  -0.061994489282369614,   0.083572722971439362,   0.180610790848731995,  -0.070595860481262207,   0.077313110232353210,   0.010885315947234631,   0.187970504164695740,  -0.095654509961605072,   0.163990378379821777,  -0.045342564582824707,  -0.134572133421897888,  -0.074215635657310486,  -0.022701907902956009,   0.063143365085124969,   0.060020480304956436,   0.067846216261386871,  -0.068684808909893036,  -0.159202456474304199,  -0.134451761841773987,   0.045189615339040756,   0.042965944856405258,  -0.026164812967181206,  -0.127674028277397156,   0.012853863649070263,  -0.022218482568860054,   0.048416726291179657,  -0.190573796629905701,  -0.036078110337257385,  -0.011245571076869965,   0.127042591571807861,   0.172350451350212097,   0.026516294106841087,   0.149230226874351501,  -0.063096076250076294,   0.061187505722045898,   0.148315832018852234,   0.118484243750572205,  -0.112641461193561554,   0.000746729958336800,  -0.105022065341472626,  -0.097805880010128021,
		  -0.145854219794273376,   0.076287709176540375,  -0.079859010875225067,  -0.062889821827411652,  -0.076370552182197571,   0.037038680166006088,  -0.177260279655456543,   0.001708839903585613,  -0.138422131538391113,  -0.016264613717794418,  -0.127266973257064819,  -0.088103614747524261,   0.189778596162796021,   0.012309768237173557,   0.179567009210586548,  -0.096472747623920441,   0.183605134487152100,   0.019621586427092552,   0.052683886140584946,  -0.110652923583984375,  -0.056467495858669281,   0.181567952036857605,   0.175373390316963196,  -0.016760094091296196,  -0.062489394098520279,  -0.192982986569404602,   0.000569162715692073,   0.028475416824221611,   0.008634324185550213,   0.095811530947685242,  -0.057359952479600906,  -0.066959537565708160,   0.086872182786464691,  -0.028212008997797966,  -0.081530392169952393,   0.034058265388011932,  -0.077949099242687225,   0.042283061891794205,  -0.033169563859701157,  -0.027330217882990837,   0.054457746446132660,  -0.171792283654212952,   0.070786796510219574,   0.112879365682601929,   0.215680092573165894,  -0.027184799313545227,   0.058457709848880768,   0.062677949666976929,
		   0.048716716468334198,   0.117109365761280060,  -0.108183346688747406,  -0.043298758566379547,   0.134361341595649719,  -0.125147327780723572,   0.001939066336490214,  -0.096767574548721313,   0.021166553720831871,   0.068058714270591736,  -0.045073300600051880,  -0.013258654624223709,   0.080832690000534058,  -0.024417880922555923,   0.016665808856487274,   0.128684610128402710,  -0.015977310016751289,  -0.148144304752349854,   0.024448510259389877,  -0.097210779786109924,   0.078121490776538849,  -0.032351404428482056,   0.075299173593521118,   0.192261010408401489,   0.071217134594917297,   0.109162263572216034,   0.085173822939395905,  -0.052410233765840530,   0.025404255837202072,  -0.095833182334899902,  -0.023485384881496429,   0.020497562363743782,   0.147337332367897034,  -0.082319438457489014,  -0.031928062438964844,   0.109291635453701019,   0.171445310115814209,   0.055088374763727188,   0.138322591781616211,   0.022439399734139442,   0.000901147781405598,   0.066313564777374268,   0.038097262382507324,   0.117337457835674286,  -0.079953983426094055,   0.111639097332954407,   0.026136839762330055,   0.067271627485752106,
		   0.114897280931472778,  -0.008911497890949249,  -0.063649877905845642,  -0.066757462918758392,   0.163333863019943237,  -0.027044076472520828,  -0.021273922175168991,  -0.004351694602519274,  -0.031613659113645554,   0.066418126225471497,  -0.067516423761844635,   0.037303946912288666,  -0.026910079643130302,   0.167646512389183044,  -0.086460344493389130,   0.040516119450330734,  -0.100604265928268433,  -0.104935862123966217,   0.103343002498149872,  -0.042382065206766129,   0.174807429313659668,  -0.142645731568336487,   0.105491541326045990,  -0.074553981423377991,  -0.020434405654668808,  -0.081862054765224457,  -0.048598583787679672,   0.055505901575088501,   0.018089847639203072,   0.025225367397069931,  -0.001183933229185641,   0.173281982541084290,  -0.037285950034856796,  -0.074489973485469818,   0.120076954364776611,   0.069431997835636139,   0.173671126365661621,   0.191958740353584290,   0.018104758113622665,   0.007054185960441828,   0.055850386619567871,  -0.130057901144027710,  -0.064408533275127411,   0.069468505680561066,   0.036204598844051361,   0.046158868819475174,  -0.165186703205108643,   0.063306227326393127,
		   0.084516972303390503,  -0.031367722898721695,  -0.027570126578211784,   0.047518905252218246,  -0.102803841233253479,  -0.015161073766648769,  -0.006094983313232660,   0.066233962774276733,   0.034494500607252121,   0.074881121516227722,  -0.118942804634571075,   0.045356258749961853,   0.089869543910026550,   0.059827778488397598,   0.145456582307815552,  -0.073342382907867432,   0.108436815440654755,   0.149080514907836914,   0.155820325016975403,   0.111884146928787231,   0.072342231869697571,  -0.010050022974610329,   0.091722473502159119,   0.064040973782539368,  -0.056561030447483063,   0.003032146021723747,  -0.157317399978637695,   0.101745948195457458,  -0.082979306578636169,  -0.172428578138351440,  -0.002866296563297510,  -0.074185997247695923,  -0.031385105103254318,  -0.016798043623566628,  -0.075134448707103729,   0.140239968895912170,  -0.156145885586738586,  -0.108074598014354706,  -0.172058075666427612,   0.053393706679344177,   0.094504117965698242,  -0.074825555086135864,   0.016622489318251610,  -0.108303055167198181,   0.020009106025099754,  -0.106822863221168518,   0.007432736456394196,   0.109017543494701385,
		   0.071541339159011841,   0.141041025519371033,  -0.100858971476554871,  -0.047672081738710403,   0.138043642044067383,  -0.117689870297908783,  -0.005607750732451677,  -0.004094983451068401,  -0.093624740839004517,  -0.028816651552915573,   0.096309550106525421,  -0.083104103803634644,   0.114837385714054108,   0.130815669894218445,  -0.076302252709865570,   0.169594794511795044,  -0.022824665531516075,   0.048523832112550735,  -0.040585532784461975,   0.102790005505084991,  -0.043613485991954803,  -0.149470955133438110,   0.058316290378570557,   0.121227696537971497,  -0.148759469389915466,  -0.095473349094390869,   0.182784587144851685,  -0.138453766703605652,   0.188696399331092834,   0.003833973780274391,  -0.072592109441757202,   0.029082214459776878,   0.013656087219715118,   0.020298205316066742,  -0.071136765182018280,  -0.009570364840328693,   0.036910813301801682,   0.047194968909025192,  -0.023383839055895805,  -0.037411909550428391,   0.115051522850990295,  -0.044932529330253601,  -0.052844762802124023,  -0.086463183164596558,  -0.070523388683795929,   0.034459657967090607,   0.027875248342752457,   0.114103473722934723,
		   0.082516878843307495,   0.097433291375637054,   0.109102070331573486,  -0.124711088836193085,   0.074564032256603241,  -0.032760482281446457,  -0.017641045153141022,  -0.088614530861377716,   0.084151148796081543,   0.045337297022342682,   0.132306233048439026,   0.025335423648357391,  -0.088321529328823090,  -0.114813335239887238,   0.062185101211071014,  -0.060672353953123093,  -0.014528783969581127,   0.127250015735626221,  -0.144342750310897827,  -0.087285548448562622,   0.049165047705173492,  -0.022823346778750420,  -0.022251542657613754,   0.123553998768329620,  -0.098407715559005737,  -0.127034440636634827,   0.109987758100032806,   0.034009918570518494,  -0.052848104387521744,  -0.100710436701774597,   0.045657772570848465,  -0.050704959779977798,   0.087048441171646118,  -0.099685221910476685,  -0.147593364119529724,   0.006135788280516863,  -0.070005990564823151,   0.100150711834430695,   0.083315417170524597,  -0.034231610596179962,  -0.100603833794593811,   0.099389880895614624,   0.102300323545932770,  -0.011990698985755444,   0.024278592318296432,  -0.097655124962329865,   0.053228765726089478,   0.043000072240829468,
		   0.083581604063510895,   0.021160950884222984,   0.080684401094913483,   0.046408295631408691,  -0.109089657664299011,   0.058938518166542053,   0.092477835714817047,  -0.027004394680261612,  -0.124979816377162933,  -0.060522720217704773,   0.057629980146884918,  -0.099485434591770172,   0.024947473779320717,  -0.049755793064832687,   0.030550017952919006,  -0.036860819905996323,   0.054317019879817963,  -0.081465736031532288,  -0.074205271899700165,  -0.022987131029367447,   0.060688167810440063,  -0.019830549135804176,   0.147518232464790344,   0.062791228294372559,  -0.072955317795276642,   0.099560029804706573,  -0.019527735188603401,  -0.037845641374588013,   0.039034795016050339,  -0.102538101375102997,   0.057088203728199005,  -0.100764155387878418,  -0.095605529844760895,   0.123972050845623016,   0.047899048775434494,   0.181860134005546570,  -0.043023765087127686,   0.035609204322099686,  -0.085401609539985657,  -0.000766980694606900,  -0.046908218413591385,   0.012853725813329220,   0.069345384836196899,  -0.035086262971162796,   0.151092052459716797,  -0.029293034225702286,   0.039848208427429199,  -0.005126014817506075,
		   0.041223101317882538,  -0.077215358614921570,  -0.033967483788728714,   0.019954251125454903,  -0.005466813221573830,  -0.119103632867336273,  -0.074312917888164520,  -0.059117469936609268,  -0.003475175704807043,  -0.088477283716201782,  -0.101054489612579346,  -0.059465579688549042,  -0.105927601456642151,  -0.010880614630877972,   0.038320451974868774,   0.174250677227973938,  -0.108755126595497131,  -0.120495833456516266,  -0.130315542221069336,  -0.118642717599868774,   0.076239712536334991,  -0.139955461025238037,   0.081897705793380737,   0.036551676690578461,  -0.031393643468618393,  -0.003741451771929860,  -0.120028890669345856,   0.044699814170598984,   0.122387886047363281,  -0.134279325604438782,  -0.097010783851146698,  -0.112734489142894745,   0.035037852823734283,  -0.099860675632953644,   0.000213384002563544,   0.132547587156295776,   0.058906529098749161,   0.154074341058731079,   0.059114001691341400,  -0.143536657094955444,  -0.150340870022773743,   0.134098634123802185,  -0.042188387364149094,  -0.022986164316534996,  -0.132588014006614685,  -0.104546919465065002,  -0.101675570011138916,  -0.091434895992279053,
		   0.062472403049468994,   0.093926310539245605,   0.004977225791662931,   0.097339555621147156,   0.025236133486032486,  -0.141628861427307129,   0.105448722839355469,  -0.092048488557338715,  -0.158298164606094360,   0.130064889788627625,  -0.095239244401454926,  -0.082933373749256134,  -0.065597884356975555,  -0.090997762978076935,  -0.112312458455562592,  -0.053329162299633026,  -0.054647799581289291,  -0.126931071281433105,   0.121434040367603302,  -0.033775743097066879,  -0.050343673676252365,  -0.051503088325262070,  -0.140386238694190979,   0.176748096942901611,  -0.021971784532070160,  -0.013160447590053082,   0.044116772711277008,   0.099450342357158661,  -0.083505704998970032,  -0.166962414979934692,  -0.122307606041431427,   0.126855656504631042,   0.021895142272114754,  -0.185547545552253723,   0.093322210013866425,   0.168690562248229980,  -0.076150350272655487,  -0.057401243597269058,   0.177252531051635742,   0.081217139959335327,  -0.140920445322990417,  -0.094656758010387421,   0.116994641721248627,  -0.038803298026323318,   0.086285702884197235,  -0.057138126343488693,   0.090037450194358826,   0.170250877737998962,
		  -0.105190150439739227,   0.132561102509498596,  -0.068412892520427704,   0.106381386518478394,  -0.102689817547798157,  -0.025300346314907074,   0.008739904500544071,  -0.154356732964515686,  -0.113588482141494751,   0.189294680953025818,   0.079353280365467072,   0.063479013741016388,   0.027586897835135460,   0.137164875864982605,   0.074407391250133514,   0.178299844264984131,   0.020389627665281296,   0.045513112097978592,   0.109914265573024750,  -0.131962925195693970,   0.109006561338901520,   0.089166402816772461,   0.061008121818304062,   0.162024959921836853,  -0.004129787907004356,  -0.054210722446441650,   0.165918022394180298,   0.093199096620082855,   0.054397866129875183,   0.073220610618591309,   0.137195527553558350,  -0.097161166369915009,   0.122781679034233093,  -0.129432693123817444,  -0.082780487835407257,   0.030558953061699867,   0.093370385468006134,   0.142396479845046997,   0.093394517898559570,   0.009184764698147774,  -0.137494757771492004,  -0.008550544269382954,   0.093237154185771942,  -0.024689814075827599,   0.102353774011135101,  -0.039500363171100616,  -0.118840336799621582,   0.157152891159057617,
		   0.010078033432364464,   0.078848764300346375,   0.043945547193288803,  -0.064723387360572815,   0.046883657574653625,   0.164790526032447815,   0.028847841545939445,  -0.096564799547195435,   0.137515053153038025,  -0.036716923117637634,  -0.102762691676616669,   0.012486282736063004,  -0.096360266208648682,  -0.049262329936027527,   0.126657545566558838,  -0.118310414254665375,   0.133407101035118103,   0.177309766411781311,   0.103280745446681976,  -0.072310499846935272,  -0.007695408537983894,  -0.019623197615146637,   0.035716511309146881,  -0.132531777024269104,   0.011630230583250523,  -0.133209854364395142,   0.016467252746224403,  -0.087011113762855530,  -0.109901390969753265,  -0.147567898035049438,  -0.122148834168910980,  -0.018454367294907570,  -0.065978109836578369,  -0.051343526691198349,   0.005467334296554327,   0.035334330052137375,  -0.036731436848640442,   0.112941794097423553,   0.002430999418720603,   0.119121626019477844,   0.012662331573665142,   0.047665070742368698,   0.024409154430031776,  -0.066510803997516632,   0.183768391609191895,  -0.146170556545257568,  -0.071949951350688934,  -0.016807133331894875,
		  -0.026305899024009705,  -0.088773742318153381,   0.164912596344947815,   0.191342830657958984,  -0.042633730918169022,   0.180958345532417297,  -0.118926309049129486,   0.076570831239223480,   0.058658022433519363,  -0.023258395493030548,  -0.023213634267449379,   0.003490535542368889,   0.142969980835914612,   0.082347020506858826,   0.037360996007919312,  -0.123277291655540466,   0.089780822396278381,   0.129936486482620239,   0.106308802962303162,  -0.079318009316921234,   0.095713339745998383,  -0.095802262425422668,  -0.087068304419517517,  -0.039246361702680588,  -0.026682216674089432,  -0.162152767181396484,   0.060259923338890076,  -0.029379656538367271,   0.029798118397593498,  -0.128399521112442017,  -0.035298407077789307,   0.098395869135856628,   0.005219290498644114,   0.163284197449684143,   0.056639011949300766,  -0.049571249634027481,   0.036436349153518677,  -0.092540040612220764,   0.012389689683914185,   0.156160056591033936,   0.016549121588468552,  -0.122984498739242554,   0.088991656899452209,   0.029271773993968964,  -0.026823967695236206,  -0.142487555742263794,  -0.169174015522003174,   0.084976181387901306,
		  -0.166372865438461304,  -0.010543842799961567,   0.094939142465591431,   0.057535465806722641,  -0.082380138337612152,  -0.073231168091297150,  -0.141018196940422058,  -0.120747655630111694,  -0.128972753882408142,   0.034626480191946030,  -0.180021360516548157,  -0.045474365353584290,   0.044280856847763062,   0.155357584357261658,   0.013799446634948254,   0.059844333678483963,   0.078772261738777161,  -0.005743701942265034,  -0.083722457289695740,  -0.158591449260711670,  -0.061245027929544449,   0.027755182236433029,   0.088939487934112549,  -0.018466196954250336,   0.057724040001630783,  -0.059898953884840012,  -0.006188940722495317,   0.066837273538112640,   0.066859498620033264,  -0.113111734390258789,  -0.001839907490648329,   0.127232030034065247,   0.038118854165077209,  -0.113843776285648346,   0.018198251724243164,  -0.023859975859522820,   0.080734916031360626,   0.075402475893497467,  -0.023747088387608528,   0.159571170806884766,  -0.043134454637765884,   0.038153864443302155,   0.194102838635444641,   0.013290262781083584,  -0.101120516657829285,  -0.137595579028129578,   0.017479438334703445,  -0.105014123022556305,
		  -0.045036617666482925,  -0.128195583820343018,   0.040824782103300095,   0.088509865105152130,  -0.114688456058502197,   0.008601055480539799,  -0.078055463731288910,  -0.108858212828636169,   0.086308225989341736,   0.081888921558856964,  -0.091160066425800323,  -0.009791814722120762,   0.201208427548408508,  -0.059570841491222382,   0.076272696256637573,  -0.042143873870372772,  -0.026770994067192078,   0.008241700939834118,   0.074452109634876251,  -0.063893377780914307,   0.011060120537877083,   0.188859328627586365,   0.169861242175102234,   0.066043697297573090,   0.092749111354351044,  -0.046426855027675629,  -0.111660793423652649,   0.168244123458862305,  -0.004453826230019331,   0.053290281444787979,  -0.127940058708190918,  -0.065319456160068512,  -0.012784557417035103,   0.050461973994970322,  -0.082196570932865143,  -0.081310935318470001,   0.009961277246475220,  -0.090398818254470825,   0.019504932686686516,   0.131153196096420288,   0.109243534505367279,  -0.117325559258460999,   0.118145488202571869,  -0.131597757339477539,   0.006788160186260939,  -0.127396166324615479,  -0.135678634047508240,   0.057623971253633499,
		  -0.107852615416049957,   0.092760168015956879,  -0.058215085417032242,   0.114949613809585571,   0.119501627981662750,   0.143803119659423828,   0.014286402612924576,  -0.038841005414724350,   0.067856423556804657,  -0.129376634955406189,  -0.099417932331562042,  -0.057737771421670914,   0.035241764038801193,  -0.126123681664466858,   0.077535167336463928,   0.089179329574108124,   0.013512440025806427,   0.129802688956260681,  -0.096253588795661926,   0.072310477495193481,   0.145099565386772156,   0.053050648421049118,   0.039391990751028061,  -0.090022258460521698,  -0.152690172195434570,  -0.191271007061004639,  -0.104560993611812592,   0.153593212366104126,   0.063446745276451111,   0.068590588867664337,   0.121713422238826752,  -0.151308789849281311,  -0.055480320006608963,   0.067066475749015808,   0.026973040774464607,  -0.012901338748633862,   0.085384219884872437,   0.036917503923177719,  -0.137661501765251160,   0.124814756214618683,  -0.005382271017879248,  -0.035807359963655472,   0.195255562663078308,  -0.149065732955932617,   0.072774879634380341,  -0.119380265474319458,  -0.166021555662155151,   0.006450329907238483,
		  -0.059485614299774170,  -0.012297815643250942,  -0.129102602601051331,  -0.006107541266828775,   0.007585760671645403,  -0.054320525377988815,   0.027783263474702835,   0.144752562046051025,   0.133746892213821411,   0.082512639462947845,  -0.000077691220212728,  -0.088787987828254700,  -0.009829699061810970,  -0.067583933472633362,   0.099462233483791351,   0.069849587976932526,  -0.007925764657557011,  -0.060298167169094086,  -0.128637626767158508,   0.143119066953659058,   0.020558213815093040,   0.100040480494499207,   0.055207382887601852,  -0.096634536981582642,   0.024032292887568474,  -0.092047750949859619,   0.091409921646118164,  -0.089500345289707184,   0.070298664271831512,   0.080443955957889557,  -0.140206709504127502,  -0.116389788687229156,   0.022350987419486046,  -0.107641518115997314,  -0.130954340100288391,  -0.045379541814327240,  -0.086340256035327911,   0.026445895433425903,  -0.121301881968975067,   0.080620303750038147,  -0.072191216051578522,   0.083171159029006958,   0.021442977711558342,  -0.044715806841850281,   0.153500571846961975,  -0.119517974555492401,   0.068203836679458618,  -0.042216423898935318,
		   0.092497795820236206,  -0.066093876957893372,  -0.059179473668336868,  -0.141997501254081726,   0.065144680440425873,  -0.157078713178634644,   0.094672106206417084,  -0.138198018074035645,   0.036942798644304276,  -0.123646616935729980,   0.048771657049655914,  -0.103495284914970398,  -0.004689646419137716,  -0.075689315795898438,   0.031073836609721184,  -0.111808598041534424,  -0.129496648907661438,  -0.174913734197616577,  -0.113828398287296295,  -0.126723527908325195,  -0.046655859798192978,  -0.058841828256845474,  -0.151155740022659302,   0.122473806142807007,   0.070292152464389801,  -0.048066068440675735,  -0.039729651063680649,  -0.151535928249359131,   0.083611451089382172,  -0.069586612284183502,   0.056454610079526901,  -0.069558598101139069,   0.037014462053775787,  -0.090212732553482056,   0.114949032664299011,   0.080877982079982758,  -0.073498308658599854,   0.018141623586416245,  -0.072747424244880676,   0.038898751139640808,  -0.033205494284629822,   0.100025691092014313,  -0.093692503869533539,   0.136461734771728516,  -0.007998095825314522,  -0.043586026877164841,  -0.091019161045551300,   0.039736207574605942,
		   0.070510305464267731,   0.145914763212203979,   0.026430847123265266,   0.159806460142135620,   0.055138137191534042,  -0.068883739411830902,  -0.107179835438728333,  -0.026729630306363106,  -0.024843661114573479,  -0.032342035323381424,   0.030567331239581108,   0.129454091191291809,  -0.079861678183078766,   0.068959385156631470,   0.075188033282756805,   0.082900896668434143,  -0.042657133191823959,   0.152151331305503845,   0.047965969890356064,  -0.037540733814239502,   0.101145662367343903,  -0.085969708859920502,   0.042321469634771347,   0.037394024431705475,   0.052454896271228790,  -0.170847609639167786,   0.028361860662698746,  -0.084629140794277191,   0.190920695662498474,  -0.067152462899684906,   0.054926071316003799,   0.138359576463699341,  -0.092954233288764954,  -0.109653495252132416,   0.062439475208520889,   0.138951182365417480,   0.099066175520420074,  -0.098266437649726868,   0.160601824522018433,   0.021179236471652985,   0.109652340412139893,  -0.154726833105087280,   0.098814003169536591,   0.072886064648628235,  -0.086576312780380249,   0.148279353976249695,  -0.004344379063695669,  -0.025188973173499107,
		   0.025670008733868599,   0.100793540477752686,  -0.106610178947448730,  -0.094971559941768646,   0.094847440719604492,   0.014025923795998096,  -0.023432383313775063,   0.101389013230800629,  -0.079603038728237152,  -0.111052609980106354,  -0.005157291889190674,  -0.019695037975907326,  -0.054845958948135376,  -0.102470189332962036,   0.069110065698623657,  -0.078212380409240723,  -0.133642643690109253,  -0.120925024151802063,  -0.062947548925876617,  -0.084484182298183441,  -0.153557911515235901,   0.119534723460674286,   0.031851634383201599,  -0.108769945800304413,  -0.014809268526732922,   0.013285190798342228,  -0.116486787796020508,  -0.051590945571660995,  -0.056566797196865082,   0.115120232105255127,  -0.151588797569274902,   0.022784609347581863,  -0.105555564165115356,  -0.046459954231977463,   0.090969905257225037,  -0.015195436775684357,  -0.032689616084098816,  -0.143211394548416138,  -0.162102833390235901,  -0.139283418655395508,   0.097148694097995758,   0.076524615287780762,  -0.154776424169540405,  -0.110621787607669830,  -0.149453148245811462,  -0.061849832534790039,  -0.122191734611988068,   0.118347525596618652,
		   0.113959632813930511,  -0.005039382260292768,  -0.089037835597991943,  -0.020861288532614708,  -0.132546573877334595,  -0.091989286243915558,  -0.063518963754177094,  -0.066834121942520142,   0.015043437480926514,   0.140819013118743896,  -0.102910302579402924,  -0.091965943574905396,  -0.001077234046533704,  -0.094411976635456085,   0.084195107221603394,  -0.056289169937372208,  -0.003501441096886992,   0.098935239017009735,   0.028290852904319763,   0.029224712401628494,   0.150717407464981079,   0.081935994327068329,  -0.090515680611133575,  -0.137832447886466980,  -0.044948156923055649,  -0.092659555375576019,  -0.035124383866786957,   0.069373495876789093,   0.148742184042930603,   0.041515827178955078,  -0.084047362208366394,  -0.095683097839355469,   0.138078078627586365,  -0.011834934353828430,  -0.101910263299942017,   0.126701414585113525,  -0.102102927863597870,  -0.096806034445762634,   0.038601785898208618,  -0.003188876202329993,   0.023007167503237724,  -0.025990627706050873,   0.122275426983833313,  -0.133301183581352234,   0.161433681845664978,   0.055446367710828781,   0.004479371942579746,  -0.097920894622802734,
		   0.014867879450321198,  -0.122390054166316986,   0.109198980033397675,  -0.058415614068508148,   0.052300009876489639,  -0.042140256613492966,  -0.031530663371086121,  -0.029547605663537979,   0.140622928738594055,   0.028365142643451691,  -0.051180332899093628,  -0.097276069223880768,  -0.135720178484916687,  -0.037462305277585983,  -0.044170953333377838,  -0.004703877028077841,   0.022597067058086395,  -0.106143802404403687,   0.106744125485420227,   0.045570928603410721,  -0.124833017587661743,  -0.051767997443675995,   0.040374673902988434,  -0.098932966589927673,  -0.048742108047008514,  -0.101370871067047119,   0.044768899679183960,  -0.108665995299816132,   0.028590500354766846,  -0.102683633565902710,   0.083624847233295441,  -0.009912180714309216,  -0.091142147779464722,   0.025746736675500870,  -0.083434388041496277,   0.029115021228790283,   0.129279747605323792,   0.093863010406494141,  -0.014553048647940159,  -0.128344848752021790,   0.015627259388566017,   0.000153727014549077,  -0.093371465802192688,   0.040691461414098740,  -0.050766617059707642,  -0.057630132883787155,   0.136598914861679077,   0.108726546168327332,
		  -0.106526672840118408,  -0.021686006337404251,   0.091789439320564270,  -0.001203757710754871,   0.015790898352861404,   0.032308157533407211,  -0.136402875185012817,  -0.013916173018515110,   0.064375475049018860,  -0.148479431867599487,   0.140926405787467957,  -0.130283534526824951,  -0.085504613816738129,  -0.067784838378429413,  -0.030285228043794632,   0.006092580966651440,  -0.088556520640850067,   0.053507041186094284,  -0.127498358488082886,  -0.067084014415740967,  -0.155893683433532715,   0.014208069071173668,  -0.154621720314025879,   0.031311586499214172,   0.134026855230331421,  -0.114981465041637421,  -0.092112749814987183,   0.049633212387561798,  -0.139403015375137329,   0.040173120796680450,   0.122182108461856842,  -0.117253795266151428,   0.101321503520011902,  -0.023876620456576347,  -0.137769490480422974,  -0.102860748767852783,   0.077303320169448853,  -0.163456723093986511,  -0.039858758449554443,  -0.130815967917442322,  -0.116155408322811127,  -0.045067768543958664,  -0.143070146441459656,   0.076360113918781281,  -0.056146424263715744,   0.087265856564044952,   0.136931285262107849,  -0.014888944104313850,
		   0.132590457797050476,   0.108400329947471619,  -0.095803461968898773,  -0.036918848752975464,  -0.090840213000774384,   0.038957547396421432,   0.022030113264918327,   0.048936627805233002,   0.098359428346157074,   0.115712545812129974,   0.120964810252189636,  -0.159196540713310242,  -0.104191586375236511,  -0.092631608247756958,   0.140672340989112854,  -0.099432319402694702,  -0.010992043651640415,   0.040089450776576996,  -0.066825531423091888,  -0.091246671974658966,   0.111842960119247437,   0.129499480128288269,   0.054871752858161926,  -0.126102373003959656,   0.030512716621160507,  -0.019626997411251068,  -0.123802505433559418,   0.113723225891590118,   0.101660735905170441,   0.149344697594642639,   0.147428721189498901,  -0.059134438633918762,   0.012487921863794327,  -0.093086704611778259,   0.035962864756584167,  -0.021579269319772720,   0.082397274672985077,   0.033095806837081909,   0.114500410854816437,   0.085865594446659088,  -0.038270302116870880,   0.042630214244127274,   0.109217554330825806,  -0.096621066331863403,   0.137658134102821350,  -0.014885521493852139,   0.104577228426933289,  -0.039830952882766724,
		   0.116673946380615234,   0.129115298390388489,   0.145292028784751892,   0.000270446733338758,   0.045109752565622330,  -0.048058640211820602,  -0.030957138165831566,  -0.014060869812965393,  -0.095148704946041107,  -0.035073999315500259,  -0.140767872333526611,   0.128659024834632874,  -0.095434829592704773,   0.074964761734008789,  -0.158707022666931152,   0.146867766976356506,   0.021395206451416016,  -0.056267477571964264,  -0.014116727747023106,  -0.025586644187569618,   0.058898314833641052,   0.087286017835140228,   0.001610490027815104,   0.079942412674427032,  -0.071084439754486084,   0.079210340976715088,   0.129809334874153137,   0.033751562237739563,  -0.095840565860271454,  -0.086912684142589569,   0.110732771456241608,  -0.042285174131393433,   0.058322574943304062,  -0.027439406141638756,   0.000222258706344292,  -0.085392326116561890,  -0.082811959087848663,  -0.058214340358972549,   0.114548735320568085,  -0.000750867708120495,  -0.138264372944831848,  -0.045024484395980835,  -0.032574806362390518,  -0.074581876397132874,  -0.000515218533109874,   0.043635249137878418,  -0.064927600324153900,   0.093004293739795685,
		  -0.045871451497077942,  -0.137496590614318848,  -0.012774482369422913,   0.128726124763488770,   0.050838693976402283,  -0.059749662876129150,  -0.079814396798610687,   0.020649807527661324,  -0.057077299803495407,  -0.062482308596372604,  -0.107785560190677643,   0.013410123065114021,   0.194363102316856384,  -0.126950621604919434,   0.147942543029785156,  -0.070097543299198151,   0.192388832569122314,   0.196792826056480408,  -0.098745316267013550,  -0.086870774626731873,   0.166909024119377136,   0.014774794690310955,   0.023926138877868652,  -0.120654508471488953,   0.101626761257648468,  -0.158558636903762817,   0.067404925823211670,   0.083610862493515015,  -0.010959541425108910,  -0.134592294692993164,  -0.001239604433067143,  -0.032066166400909424,  -0.038270074874162674,  -0.063657172024250031,   0.128671273589134216,  -0.035144452005624771,  -0.022271463647484779,   0.064861454069614410,  -0.042980004101991653,   0.120303586125373840,   0.068851098418235779,  -0.134722217917442322,   0.133566349744796753,  -0.107268072664737701,   0.084760032594203949,   0.047888055443763733,  -0.032407522201538086,  -0.066019520163536072,
		   0.083903774619102478,   0.054793085902929306,  -0.106372296810150146,   0.077187187969684601,  -0.013635170646011829,  -0.059121903032064438,   0.066527597606182098,  -0.018192075192928314,   0.011861535720527172,   0.094326108694076538,   0.073441758751869202,  -0.086919829249382019,   0.119798116385936737,   0.045653972774744034,  -0.089597798883914948,  -0.006280796136707067,   0.182716414332389832,   0.015534870326519012,   0.080918505787849426,  -0.144129589200019836,   0.134787708520889282,   0.096541568636894226,  -0.036487270146608353,  -0.113494679331779480,   0.078296981751918793,   0.070226833224296570,   0.063882730901241302,   0.006136957090348005,  -0.078455477952957153,  -0.012720220722258091,  -0.030197530984878540,   0.003836815478280187,   0.158120229840278625,  -0.057747673243284225,  -0.073220275342464447,  -0.027102297171950340,  -0.011289544403553009,   0.086157478392124176,   0.139921829104423523,  -0.061359833925962448,   0.012168126180768013,  -0.114625081419944763,  -0.075933396816253662,   0.112631738185882568,  -0.059573333710432053,   0.131692335009574890,   0.066332317888736725,  -0.014467576518654823,
		   0.035410404205322266,   0.121438845992088318,   0.066757179796695709,   0.030864911153912544,  -0.060126662254333496,   0.076469689607620239,   0.129765406250953674,   0.110281065106391907,   0.076894588768482208,   0.010380514897406101,  -0.056372802704572678,   0.119468986988067627,  -0.105506598949432373,   0.161128744482994080,  -0.132884338498115540,   0.131989806890487671,   0.041814062744379044,  -0.090684190392494202,   0.102878846228122711,   0.060868375003337860,   0.001250979723408818,   0.009866307489573956,   0.030292632058262825,   0.154416203498840332,   0.024316836148500443,  -0.039043530821800232,   0.084726437926292419,  -0.020080979913473129,  -0.108777284622192383,   0.067240171134471893,  -0.136501327157020569,   0.078355617821216583,   0.131915181875228882,  -0.170345753431320190,   0.015359808690845966,   0.061778880655765533,  -0.011250506155192852,  -0.004709157161414623,  -0.014775454066693783,  -0.056884579360485077,  -0.028921892866492271,  -0.055300705134868622,   0.033920880407094955,   0.004284569993615150,  -0.162725776433944702,   0.084682986140251160,  -0.068660594522953033,  -0.040394391864538193,
		  -0.063085652887821198,   0.022842818871140480,  -0.105553522706031799,   0.101861841976642609,  -0.105041116476058960,   0.047318212687969208,  -0.094290658831596375,  -0.051945324987173080,   0.085154965519905090,  -0.031992565840482712,  -0.002306041307747364,   0.095360100269317627,   0.093953087925910950,  -0.078773185610771179,  -0.028659610077738762,  -0.030579676851630211,  -0.099141180515289307,  -0.138725027441978455,  -0.072687909007072449,   0.015122334472835064,   0.021901184692978859,  -0.032799463719129562,  -0.065797470510005951,   0.177768006920814514,  -0.138961091637611389,  -0.090681686997413635,  -0.054292909801006317,   0.134882569313049316,   0.003870117012411356,  -0.116720028221607208,  -0.137037351727485657,  -0.006020369939506054,   0.040912039577960968,   0.122023403644561768,   0.053596109151840210,   0.039131738245487213,   0.008737677708268166,   0.130591288208961487,  -0.017074944451451302,   0.004736952017992735,  -0.071224749088287354,   0.040882971137762070,  -0.046049751341342926,   0.059297561645507812,  -0.132335931062698364,   0.093432240188121796,  -0.112593971192836761,  -0.049916278570890427,
		  -0.096994549036026001,  -0.036314841359853745,  -0.063349753618240356,   0.086513683199882507,  -0.033885132521390915,   0.188660547137260437,   0.077817574143409729,   0.027294900268316269,   0.132184281945228577,   0.057646300643682480,  -0.002757359063252807,   0.131130069494247437,   0.050460170954465866,  -0.126631528139114380,   0.014375158585608006,   0.025536838918924332,   0.029006320983171463,   0.147541731595993042,   0.022234002128243446,  -0.141524493694305420,   0.038867637515068054,  -0.080550491809844971,  -0.088063195347785950,   0.062092900276184082,   0.074612282216548920,  -0.137233033776283264,   0.100161477923393250,  -0.064298011362552643,  -0.006449468899518251,   0.017832536250352859,  -0.028637314215302467,   0.166580870747566223,  -0.057737905532121658,   0.069764472544193268,  -0.073942109942436218,   0.033861882984638214,   0.135545775294303894,   0.012610123492777348,  -0.070557087659835815,   0.153484299778938293,   0.068453334271907806,  -0.171382471919059753,  -0.025351870805025101,   0.004022095818072557,   0.119526542723178864,  -0.072362072765827179,   0.068671435117721558,   0.018557846546173096,
		  -0.059275973588228226,   0.106923080980777740,  -0.003134105121716857,   0.072457380592823029,   0.045159712433815002,  -0.015158228576183319,  -0.112631067633628845,   0.133551195263862610,   0.088725820183753967,   0.111502453684806824,  -0.105788029730319977,  -0.117615021765232086,  -0.102459013462066650,   0.076796755194664001,  -0.104927562177181244,   0.050550870597362518,  -0.123759783804416656,   0.075814448297023773,  -0.129090994596481323,   0.110716730356216431,  -0.068944007158279419,  -0.119472280144691467,   0.155269220471382141,   0.025531860068440437,   0.036245759576559067,   0.109240338206291199,   0.106785617768764496,   0.032523270696401596,  -0.014185882173478603,   0.046595130115747452,   0.110196955502033234,  -0.042474273592233658,  -0.131212502717971802,   0.035722333937883377,   0.078063860535621643,   0.043714035302400589,  -0.070758454501628876,   0.018278285861015320,  -0.025412198156118393,   0.038106020539999008,   0.067986398935317993,  -0.104125298559665680,  -0.083322681486606598,   0.039033707231283188,  -0.004161659628152847,   0.170622840523719788,   0.126520037651062012,   0.066284164786338806,
		   0.025209177285432816,   0.021498097106814384,   0.061141137033700943,   0.093749567866325378,   0.093772284686565399,  -0.143101170659065247,   0.051170974969863892,   0.032905448228120804,  -0.051955603063106537,   0.170375585556030273,  -0.024563582614064217,  -0.064397506415843964,  -0.094609752297401428,  -0.072276823222637177,   0.062821909785270691,   0.147469937801361084,  -0.058194093406200409,  -0.006339739076793194,   0.007400164846330881,  -0.063245616853237152,   0.172403991222381592,  -0.078834861516952515,  -0.136337876319885254,   0.034254156053066254,  -0.140364617109298706,   0.056835360825061798,   0.008628043346107006,   0.035923816263675690,   0.160743072628974915,   0.044850621372461319,  -0.018198447301983833,   0.142649486660957336,  -0.086056478321552277,   0.046768825501203537,   0.145995408296585083,  -0.048524681478738785,   0.085312061011791229,   0.015235474333167076,   0.137523949146270752,  -0.020060600712895393,  -0.158529639244079590,   0.072438478469848633,   0.032850023359060287,   0.162664502859115601,  -0.006364626344293356,   0.096769988536834717,   0.098615854978561401,   0.084536053240299225,
		   0.090628601610660553,  -0.069508396089076996,  -0.093569897115230560,  -0.119205810129642487,  -0.039234288036823273,   0.043227203190326691,   0.041602157056331635,  -0.053346730768680573,  -0.032099541276693344,  -0.035056430846452713,   0.123217478394508362,  -0.139736995100975037,  -0.084976673126220703,  -0.062107745558023453,  -0.062614470720291138,  -0.002910958835855126,   0.122211843729019165,  -0.007675888948142529,  -0.135741308331489563,  -0.064733698964118958,  -0.079381100833415985,   0.155988246202468872,   0.024475365877151489,   0.025494879111647606,   0.014374759979546070,  -0.020399972796440125,   0.116968773305416107,   0.068923808634281158,  -0.103622190654277802,   0.090212233364582062,   0.068200185894966125,  -0.025073662400245667,   0.047631427645683289,   0.079717934131622314,   0.032405093312263489,  -0.154252335429191589,   0.041744146496057510,  -0.019179914146661758,  -0.069476447999477386,  -0.141147419810295105,  -0.023827731609344482,  -0.008668947964906693,   0.086362794041633606,  -0.151018410921096802,   0.055000878870487213,   0.074102759361267090,   0.124635912477970123,   0.116450116038322449,
	};
	static const double bias04[]=
	{
		  -0.034130595624446869,   0.029054811224341393,  -0.037332169711589813,   0.151589363813400269,   0.003486730391159654,  -0.068020768463611603,  -0.042420197278261185,  -0.037854548543691635,   0.039538435637950897,  -0.048481930047273636,   0.017110921442508698,  -0.005674703512340784,  -0.071218222379684448,  -0.056014560163021088,  -0.034004241228103638,  -0.110808178782463074,   0.145048320293426514,  -0.012399022467434406,   0.018387105315923691,   0.128557085990905762,   0.057746335864067078,  -0.038086026906967163,  -0.069653786718845367,  -0.049930706620216370,   0.023829832673072815,  -0.032151691615581512,   0.030234752222895622,   0.041474666446447372,   0.022298093885183334,   0.108804434537887573,   0.095684424042701721,  -0.111073680222034454,   0.023047225549817085,   0.138166531920433044,  -0.134883478283882141,  -0.099651351571083069,  -0.108709350228309631,  -0.099786058068275452,  -0.068823628127574921,  -0.083544485270977020,   0.032044861465692520,  -0.002699437784031034,   0.101552009582519531,   0.026411810889840126,   0.045338690280914307,  -0.021607872098684311,   0.035388603806495667,   0.029285654425621033,
	};
	static const double weight05[]=
	{
		   0.125760316848754883,   0.057287089526653290,  -0.037013202905654907,  -0.147217214107513428,   0.109005786478519440,   0.126722306013107300,  -0.032967098057270050,   0.103656142950057983,   0.147552564740180969,  -0.081429027020931244,  -0.070205852389335632,   0.064507231116294861,  -0.142636105418205261,   0.023179547861218452,   0.077911585569381714,   0.087157197296619415,   0.003896635957062244,   0.082279272377490997,  -0.023707844316959381,  -0.008571821264922619,   0.028974642977118492,  -0.141181811690330505,  -0.099949195981025696,  -0.057428259402513504,  -0.016780119389295578,  -0.077977880835533142,  -0.131312549114227295,  -0.088852800428867340,   0.152021735906600952,  -0.017162062227725983,  -0.063092492520809174,   0.136064127087593079,   0.006434035953134298,   0.012546149082481861,  -0.038214720785617828,  -0.001300349133089185,  -0.123613081872463226,  -0.099578849971294403,   0.008341981098055840,  -0.115150392055511475,  -0.006442510522902012,   0.150654584169387817,   0.138805493712425232,   0.084267184138298035,   0.025370500981807709,  -0.048859238624572754,  -0.025545949116349220,  -0.031818635761737823,
		  -0.011767120100557804,  -0.011042734608054161,  -0.152901470661163330,   0.059028197079896927,   0.000718465191312134,   0.148165285587310791,  -0.023353848606348038,   0.152457356452941895,  -0.031558677554130554,  -0.108999118208885193,  -0.083295702934265137,  -0.110811188817024231,  -0.001601499621756375,   0.080615207552909851,  -0.130309894680976868,   0.128557115793228149,   0.085119523108005524,   0.016511950641870499,  -0.009116688743233681,   0.077697493135929108,  -0.069430753588676453,   0.043927095830440521,   0.139697656035423279,   0.043525990098714828,   0.000378608121536672,   0.146226018667221069,  -0.002618998987600207,   0.033892933279275894,   0.076588295400142670,   0.167273357510566711,   0.096124671399593353,   0.030816359445452690,  -0.112690895795822144,   0.157034754753112793,  -0.115138895809650421,   0.114653199911117554,  -0.017618484795093536,  -0.141196802258491516,   0.141280442476272583,   0.084109723567962646,  -0.039331167936325073,   0.050409711897373199,  -0.085908845067024231,  -0.181386455893516541,   0.098388656973838806,   0.082935050129890442,   0.049946725368499756,   0.013367812149226665,
		   0.159054338932037354,   0.116161100566387177,   0.117826491594314575,  -0.035864576697349548,  -0.030924653634428978,  -0.060612294822931290,  -0.149438485503196716,   0.034691281616687775,   0.034525923430919647,  -0.040183629840612411,  -0.094226650893688202,   0.040028031915426254,  -0.149315699934959412,  -0.031460307538509369,   0.085604950785636902,  -0.028927486389875412,  -0.043838158249855042,  -0.113250620663166046,   0.165336251258850098,  -0.084235668182373047,  -0.047913718968629837,  -0.072985433042049408,  -0.024715349078178406,   0.052015662193298340,   0.045620232820510864,   0.044920649379491806,  -0.090453073382377625,  -0.056675132364034653,  -0.075610242784023285,  -0.037682741880416870,   0.000771106802858412,  -0.136099681258201599,  -0.044972769916057587,  -0.044180441647768021,   0.086516223847866058,   0.002686975756660104,  -0.134808480739593506,   0.093819335103034973,  -0.108581602573394775,   0.081246979534626007,  -0.127859860658645630,   0.159906968474388123,   0.083793044090270996,   0.149025350809097290,  -0.038457088172435760,   0.028480038046836853,  -0.021358985453844070,   0.103025391697883606,
		  -0.053294252604246140,   0.030950624495744705,   0.033162843436002731,  -0.075693823397159576,   0.082893028855323792,   0.057026050984859467,  -0.018408395349979401,  -0.159536078572273254,   0.012124463915824890,   0.015644552186131477,   0.022814750671386719,   0.076668038964271545,  -0.050678540021181107,   0.017427092418074608,  -0.111445866525173187,   0.067022547125816345,  -0.007997843436896801,  -0.107940956950187683,   0.079921916127204895,  -0.033155988901853561,  -0.136615678668022156,   0.009168016724288464,  -0.127476289868354797,   0.082119241356849670,  -0.095490060746669769,  -0.054202370345592499,   0.090184353291988373,  -0.181567355990409851,  -0.024996308609843254,  -0.188262268900871277,   0.015051767230033875,   0.129204377532005310,  -0.117550030350685120,  -0.022515740245580673,   0.059743858873844147,   0.081841677427291870,   0.123344846069812775,   0.081431105732917786,  -0.103238135576248169,  -0.131896898150444031,  -0.058073405176401138,  -0.007715099956840277,   0.109135419130325317,   0.029417233541607857,   0.091276116669178009,  -0.047808717936277390,  -0.023039985448122025,  -0.107076384127140045,
		   0.179705023765563965,   0.107142597436904907,   0.070926256477832794,  -0.035090383142232895,   0.050963405519723892,   0.060873564332723618,   0.084163859486579895,  -0.105991885066032410,   0.092276461422443390,  -0.087297566235065460,   0.067537397146224976,  -0.085708126425743103,   0.053430084139108658,   0.038399450480937958,   0.066307529807090759,   0.116340972483158112,  -0.071653880178928375,   0.135451704263687134,   0.121991641819477081,  -0.147019311785697937,   0.033363096415996552,  -0.101893298327922821,  -0.021594591438770294,   0.104436069726943970,  -0.006917378399521112,   0.141161352396011353,   0.068416081368923187,   0.092089317739009857,   0.057356555014848709,  -0.167292535305023193,   0.014439661987125874,  -0.160904422402381897,  -0.072386085987091064,   0.156674072146415710,   0.077899634838104248,  -0.051042791455984116,   0.130020886659622192,  -0.096586637198925018,  -0.018593855202198029,  -0.075164563953876495,  -0.109333954751491547,  -0.126118138432502747,  -0.035091176629066467,   0.030615009367465973,  -0.081994615495204926,  -0.032616857439279556,  -0.006896664854139090,  -0.091483324766159058,
		  -0.086111307144165039,  -0.079820677638053894,   0.099121920764446259,   0.064242184162139893,   0.014167348854243755,  -0.146554306149482727,   0.016462935134768486,  -0.057580657303333282,  -0.068327985703945160,   0.017013764008879662,   0.096368923783302307,  -0.131445527076721191,  -0.065159872174263000,  -0.044189535081386566,   0.105709016323089600,  -0.129317224025726318,  -0.087461724877357483,   0.062169000506401062,  -0.043622102588415146,  -0.058918796479701996,  -0.102557346224784851,   0.075904361903667450,   0.081546582281589508,   0.147800043225288391,  -0.114268094301223755,  -0.170745402574539185,  -0.150978386402130127,  -0.139384895563125610,   0.116819210350513458,  -0.090984269976615906,  -0.136886060237884521,  -0.054360426962375641,   0.150512218475341797,   0.006761693861335516,  -0.032655920833349228,   0.028433863073587418,  -0.005386666394770145,  -0.031512100249528885,   0.010808953084051609,  -0.023394098505377769,  -0.123765267431735992,  -0.002951445523649454,   0.078391186892986298,   0.072957366704940796,  -0.047964323312044144,   0.004145558923482895,  -0.108615629374980927,  -0.015290946699678898,
		   0.140158638358116150,  -0.126715049147605896,   0.119885273277759552,   0.047799985855817795,   0.030933408066630363,   0.061071701347827911,   0.066728509962558746,   0.047986563295125961,  -0.064789168536663055,   0.063231863081455231,   0.135060399770736694,   0.095321536064147949,   0.160541266202926636,  -0.028454337269067764,  -0.078014813363552094,   0.001147187314927578,  -0.038757201284170151,  -0.115747250616550446,  -0.101101674139499664,   0.084385879337787628,  -0.005675940308719873,  -0.083552643656730652,   0.143731996417045593,  -0.105688862502574921,  -0.020519465208053589,   0.146625861525535583,  -0.012178278528153896,   0.197279661893844604,  -0.002225597621873021,   0.190262496471405029,   0.015911057591438293,  -0.145113825798034668,   0.060033768415451050,  -0.039842151105403900,  -0.076177924871444702,   0.172163113951683044,  -0.036065582185983658,   0.030333185568451881,  -0.020089108496904373,   0.114231638610363007,   0.022259801626205444,   0.024551492184400558,  -0.084307461977005005,   0.097090989351272583,  -0.069231063127517700,   0.019125120714306831,  -0.051080036908388138,  -0.045473523437976837,
		  -0.051036227494478226,  -0.142077997326850891,   0.000028752059733961,   0.137437403202056885,  -0.008655147626996040,   0.172949165105819702,   0.109385132789611816,   0.059257458895444870,   0.097143530845642090,   0.097031898796558380,   0.058226332068443298,  -0.144350156188011169,  -0.058458134531974792,   0.110753260552883148,   0.104168653488159180,   0.066289469599723816,  -0.011229279451072216,   0.044206749647855759,   0.072242453694343567,   0.082073584198951721,  -0.095894373953342438,  -0.010032482445240021,   0.084254421293735504,  -0.118358179926872253,   0.103866308927536011,   0.010757513344287872,   0.090588584542274475,   0.062295250594615936,   0.149854198098182678,   0.075757697224617004,   0.095280610024929047,   0.107399836182594299,  -0.149755075573921204,  -0.109437040984630585,  -0.071594052016735077,  -0.062475960701704025,  -0.137681826949119568,   0.068430580198764801,  -0.063369803130626678,   0.068871207535266876,   0.125830650329589844,   0.076615750789642334,  -0.009979401715099812,   0.012839121744036674,   0.033222921192646027,   0.074171729385852814,   0.026695169508457184,   0.040418677031993866,
		   0.050623424351215363,   0.027024464681744576,   0.076374061405658722,  -0.104779317975044250,  -0.067616783082485199,  -0.177755385637283325,  -0.068105056881904602,  -0.022948987782001495,  -0.089849650859832764,   0.024861788377165794,  -0.005237287841737270,   0.136749416589736938,  -0.090663544833660126,   0.044155474752187729,   0.023610493168234825,   0.105739414691925049,   0.032987937331199646,  -0.128370061516761780,  -0.133040994405746460,  -0.018420804291963577,  -0.054983753710985184,   0.013782372698187828,  -0.095382995903491974,  -0.078693822026252747,   0.099832631647586823,  -0.008521492592990398,   0.006484949495643377,  -0.186580866575241089,  -0.136243358254432678,   0.020113309845328331,  -0.146012842655181885,  -0.113950014114379883,  -0.037775870412588120,   0.003759699873626232,   0.042438525706529617,  -0.124786064028739929,  -0.078884296119213104,   0.106917858123779297,  -0.118924245238304138,   0.122894540429115295,  -0.180200323462486267,  -0.026625910773873329,  -0.102991700172424316,  -0.061750974506139755,  -0.000034404940379318,  -0.045030094683170319,  -0.058790747076272964,   0.052717823535203934,
		  -0.018409064039587975,   0.003730033524334431,   0.072156026959419250,  -0.015802400186657906,   0.016926223412156105,  -0.007981975562870502,   0.055953454226255417,   0.003952544182538986,   0.123975217342376709,   0.044082920998334885,  -0.016047617420554161,   0.090731650590896606,  -0.134846538305282593,  -0.128743141889572144,  -0.022798592224717140,  -0.016008272767066956,   0.007316206581890583,   0.025832751765847206,   0.103156335651874542,  -0.133431643247604370,  -0.033331312239170074,  -0.137040719389915466,   0.019022919237613678,   0.017094964161515236,  -0.150686055421829224,   0.095918931066989899,  -0.089583158493041992,  -0.175558567047119141,  -0.101868383586406708,  -0.196231111884117126,  -0.173562929034233093,  -0.082950390875339508,  -0.129603743553161621,  -0.086262084543704987,   0.062238231301307678,  -0.003668570425361395,   0.111008420586585999,   0.024703098461031914,  -0.153945818543434143,   0.051826216280460358,  -0.081384591758251190,   0.061594355851411819,  -0.129695236682891846,  -0.104645878076553345,  -0.094656825065612793,   0.036292154341936111,  -0.138896167278289795,   0.145740747451782227,
		   0.030852548778057098,  -0.135064363479614258,   0.050580546259880066,   0.067042976617813110,  -0.061479244381189346,  -0.053837332874536514,  -0.090125925838947296,   0.080729030072689056,  -0.089202508330345154,  -0.009035354480147362,  -0.136633291840553284,  -0.076043143868446350,  -0.113053500652313232,   0.131841033697128296,  -0.043450951576232910,  -0.008918107487261295,  -0.197461843490600586,  -0.078622721135616302,  -0.129093110561370850,   0.000732999411411583,  -0.051871020346879959,  -0.017444645985960960,  -0.073771618306636810,   0.007955367676913738,  -0.002096599899232388,  -0.020397443324327469,   0.105342447757720947,  -0.116219870746135712,   0.045526120811700821,  -0.007355166133493185,  -0.188233375549316406,   0.052462410181760788,  -0.027040062472224236,  -0.040031373500823975,  -0.128367468714714050,  -0.121128983795642853,   0.099104858934879303,  -0.047887429594993591,  -0.073315776884555817,   0.136899888515472412,   0.074488840997219086,   0.027908271178603172,  -0.011935876682400703,   0.066264003515243530,   0.078887552022933960,  -0.080852448940277100,   0.104698807001113892,   0.091398812830448151,
		  -0.147501081228256226,  -0.082837454974651337,   0.000298846251098439,   0.117286890745162964,   0.132723212242126465,   0.185562133789062500,  -0.069883443415164948,   0.013015418313443661,   0.022108407691121101,   0.090225525200366974,   0.141081675887107849,   0.036981679499149323,   0.128392681479454041,   0.039313245564699173,  -0.107747554779052734,   0.013088658452033997,   0.180179581046104431,   0.055889409035444260,  -0.110798768699169159,   0.142810970544815063,  -0.007093240972608328,   0.098680883646011353,   0.029867261648178101,  -0.022457432001829147,  -0.031830389052629471,  -0.088366389274597168,   0.122707106173038483,   0.075110740959644318,  -0.056078385561704636,   0.052966650575399399,   0.184521362185478210,   0.084980472922325134,  -0.019290838390588760,   0.170387014746665955,  -0.067600153386592865,   0.003012976143509150,  -0.099007256329059601,  -0.069718576967716217,  -0.034955359995365143,  -0.196476861834526062,   0.186345532536506653,   0.040341518819332123,   0.040733907371759415,   0.030063044279813766,   0.165199473500251770,  -0.122848570346832275,  -0.105919234454631805,  -0.063829652965068817,
		  -0.034839220345020294,  -0.147128596901893616,  -0.039266858249902725,   0.000948377128224820,  -0.055319562554359436,  -0.018437189981341362,   0.095356687903404236,  -0.093367092311382294,  -0.042269289493560791,  -0.137276411056518555,  -0.028153927996754646,  -0.170537069439888000,   0.159801900386810303,   0.025047190487384796,  -0.173032000660896301,  -0.113434813916683197,   0.119331628084182739,  -0.022112026810646057,  -0.188236922025680542,  -0.015379837714135647,  -0.019853517413139343,  -0.128240063786506653,   0.000081925543781836,  -0.137059748172760010,  -0.113086588680744171,   0.088771812617778778,   0.152997806668281555,   0.096510097384452820,   0.021140407770872116,  -0.026663886383175850,   0.154668480157852173,   0.134269088506698608,  -0.105201102793216705,   0.128918483853340149,  -0.091545090079307556,   0.044274248182773590,  -0.112147480249404907,  -0.071441560983657837,   0.171329945325851440,   0.002156680217012763,   0.139538288116455078,   0.028360124677419662,   0.001426416449248791,  -0.014914260245859623,   0.147268891334533691,   0.085931293666362762,  -0.191676586866378784,  -0.032592494040727615,
		   0.024797316640615463,  -0.038615770637989044,  -0.022177202627062798,   0.016181791201233864,  -0.098260328173637390,   0.079380117356777191,   0.061140380799770355,  -0.035898204892873764,   0.003547305706888437,   0.128203734755516052,  -0.057527396827936172,   0.116001620888710022,  -0.157102599740028381,   0.011041728779673576,  -0.099361419677734375,   0.089021459221839905,   0.006068263668566942,   0.090681344270706177,  -0.016918642446398735,  -0.030457478016614914,  -0.117766730487346649,  -0.056901667267084122,   0.008012713864445686,   0.117388479411602020,   0.087587550282478333,  -0.042024638503789902,   0.122043661773204803,  -0.047060955315828323,  -0.024384835734963417,   0.014654795639216900,  -0.170536324381828308,   0.143213808536529541,  -0.067772440612316132,  -0.110473930835723877,   0.055986642837524414,   0.106329299509525299,  -0.087907500565052032,  -0.026250028982758522,   0.080129668116569519,  -0.019342960789799690,  -0.034858375787734985,  -0.042378839105367661,  -0.025740150362253189,  -0.016874836757779121,   0.026246517896652222,   0.090241238474845886,   0.008445178158581257,  -0.053077809512615204,
		  -0.008196485228836536,   0.104762606322765350,   0.057622767984867096,   0.009077465161681175,   0.038718249648809433,   0.099314421415328979,   0.061126582324504852,  -0.020646376535296440,   0.156143799424171448,   0.015700602903962135,  -0.070480272173881531,  -0.012444704771041870,  -0.143561258912086487,  -0.071781061589717865,   0.061351843178272247,   0.175168350338935852,   0.030965775251388550,   0.056779850274324417,  -0.082375779747962952,  -0.066504232585430145,   0.001152194454334676,   0.013692426495254040,   0.041737817227840424,   0.034015636891126633,   0.074163138866424561,  -0.060905821621417999,   0.023229975253343582,   0.125257804989814758,   0.140498891472816467,  -0.081275612115859985,  -0.199046477675437927,   0.120624750852584839,  -0.054572135210037231,  -0.087325327098369598,  -0.089651763439178467,   0.071907192468643188,  -0.075997270643711090,  -0.088460572063922882,  -0.120847508311271667,  -0.049266900867223740,  -0.182412594556808472,   0.029383135959506035,   0.141400888562202454,  -0.132023200392723083,   0.151763543486595154,  -0.087367519736289978,  -0.079177714884281158,  -0.161659836769104004,
		   0.167663022875785828,   0.146555140614509583,   0.072104312479496002,   0.126662850379943848,  -0.181285485625267029,   0.013775008730590343,   0.131432473659515381,  -0.128823712468147278,   0.186878338456153870,   0.038487125188112259,   0.147447690367698669,   0.064617693424224854,  -0.164578810334205627,   0.068905793130397797,  -0.083415426313877106,   0.176602333784103394,   0.030115382745862007,   0.186118468642234802,  -0.068212650716304779,  -0.028695723041892052,   0.152179971337318420,  -0.121587581932544708,   0.009629280306398869,   0.003397878957912326,   0.057488586753606796,   0.113595478236675262,  -0.165986388921737671,  -0.004767374601215124,  -0.003075294895097613,  -0.167248755693435669,   0.086167886853218079,  -0.019865004345774651,   0.057237189263105392,   0.143848463892936707,  -0.058607134968042374,   0.069938033819198608,  -0.010588489472866058,  -0.049588993191719055,  -0.019171476364135742,   0.161076605319976807,  -0.185722216963768005,   0.060363721102476120,   0.108338981866836548,   0.091192573308944702,   0.048953395336866379,  -0.007463021669536829,   0.028718451038002968,  -0.047020494937896729,
		   0.075647883117198944,   0.174593657255172729,   0.124522335827350616,   0.133163645863533020,  -0.158543065190315247,   0.011133794672787189,  -0.088512286543846130,  -0.085446231067180634,   0.134431600570678711,  -0.116025254130363464,  -0.033927962183952332,   0.073359176516532898,  -0.027250485494732857,  -0.077419951558113098,   0.178007051348686218,  -0.035720877349376678,  -0.108997896313667297,   0.116205029189586639,   0.175346106290817261,  -0.146477758884429932,   0.092967830598354340,   0.025022897869348526,   0.094546668231487274,   0.121871381998062134,   0.157032683491706848,  -0.019838875159621239,   0.060766883194446564,   0.140616595745086670,   0.088973924517631531,  -0.104551210999488831,  -0.158758327364921570,  -0.052907075732946396,   0.098094142973423004,   0.134438142180442810,  -0.123416587710380554,   0.048613268882036209,   0.033563118427991867,  -0.075109302997589111,  -0.125631242990493774,   0.178921997547149658,  -0.169201850891113281,   0.161242201924324036,   0.050105199217796326,   0.054803829640150070,   0.018521595746278763,  -0.046305548399686813,   0.145640805363655090,  -0.120671346783638000,
		  -0.029268816113471985,  -0.076774716377258301,   0.060661494731903076,   0.140582710504531860,  -0.096626169979572296,  -0.042249813675880432,  -0.079874768853187561,   0.017521785572171211,   0.063336715102195740,   0.114520460367202759,  -0.019152639433741570,  -0.025322511792182922,  -0.016071461141109467,  -0.027732027694582939,   0.067740775644779205,   0.150431439280509949,  -0.195687755942344666,   0.089726701378822327,   0.053312826901674271,  -0.031566325575113297,  -0.022428998723626137,   0.010230986401438713,  -0.115978337824344635,   0.059309393167495728,   0.125785231590270996,   0.077855579555034637,   0.097792498767375946,   0.073133036494255066,   0.068751983344554901,  -0.164093792438507080,  -0.141046226024627686,  -0.063467517495155334,  -0.085254579782485962,  -0.096773967146873474,  -0.134500458836555481,   0.104295916855335236,   0.051080215722322464,   0.032099429517984390,  -0.143787726759910583,   0.096866086125373840,   0.005324719473719597,  -0.062327526509761810,  -0.026566522195935249,   0.133585304021835327,  -0.068751312792301178,   0.162116706371307373,   0.041809253394603729,   0.057268213480710983,
		  -0.039881099015474319,   0.058326181024312973,  -0.012415111996233463,  -0.023685617372393608,  -0.081147588789463043,  -0.157947212457656860,  -0.138161823153495789,  -0.140548810362815857,   0.122103624045848846,   0.081510417163372040,  -0.018823025748133659,   0.117471374571323395,   0.109949737787246704,   0.049420028924942017,  -0.030161617323756218,  -0.014501977711915970,   0.044113136827945709,   0.073371842503547668,   0.077196694910526276,  -0.067642197012901306,   0.004706640262156725,  -0.012425979599356651,  -0.148120984435081482,  -0.139974191784858704,   0.009047679603099823,  -0.068094417452812195,  -0.079598061740398407,   0.030776662752032280,  -0.116706036031246185,   0.062647193670272827,   0.091976627707481384,  -0.060500759631395340,  -0.129306092858314514,   0.115793660283088684,   0.101450085639953613,   0.086796127259731293,  -0.044980417937040329,  -0.031396009027957916,  -0.114083096385002136,   0.127208992838859558,   0.011252180673182011,  -0.071401849389076233,   0.110674574971199036,  -0.026034474372863770,   0.066792853176593781,  -0.063218273222446442,   0.079875372350215912,   0.106557734310626984,
		   0.093300029635429382,   0.137672603130340576,  -0.119160629808902740,  -0.001400372013449669,   0.102691695094108582,  -0.000254994607530534,   0.097898609936237335,   0.022466329857707024,  -0.062605336308479309,  -0.047499991953372955,  -0.033712100237607956,   0.028719702735543251,   0.121818087995052338,   0.078475683927536011,  -0.044010140001773834,  -0.048759654164314270,   0.168569609522819519,  -0.009683768264949322,   0.022617451846599579,   0.188333123922348022,  -0.116430789232254028,  -0.069364182651042938,  -0.094999551773071289,  -0.061160203069448471,   0.065513275563716888,  -0.050354443490505219,   0.051989395171403885,   0.043076291680335999,   0.134436175227165222,   0.139587283134460449,  -0.031627252697944641,  -0.014351619407534599,   0.061419658362865448,  -0.032958175987005234,  -0.069268241524696350,   0.003003309480845928,   0.026763625442981720,  -0.050918340682983398,   0.187597781419754028,  -0.150942564010620117,   0.179047673940658569,   0.016298610717058182,   0.067707210779190063,   0.054458230733871460,   0.048453357070684433,   0.075583592057228088,   0.052173152565956116,   0.124360613524913788,
		   0.040738042443990707,  -0.083915792405605316,  -0.058499887585639954,  -0.003046571509912610,   0.038674086332321167,   0.053884573280811310,  -0.122775934636592865,  -0.146699324250221252,   0.068063125014305115,  -0.122895576059818268,   0.120938107371330261,   0.184331566095352173,   0.034311763942241669,  -0.074258044362068176,  -0.071044944226741791,   0.092016726732254028,   0.005602969788014889,   0.021206064149737358,   0.162448301911354065,  -0.138965532183647156,   0.169139996170997620,   0.153301924467086792,  -0.122202783823013306,   0.021159829571843147,   0.065940819680690765,   0.046234440058469772,   0.053494859486818314,  -0.121048428118228912,  -0.043277829885482788,  -0.077095255255699158,  -0.058322221040725708,  -0.045811999589204788,   0.049846161156892776,   0.121447011828422546,   0.098099231719970703,  -0.106174483895301819,   0.009149380959570408,  -0.054848361760377884,  -0.171980232000350952,  -0.026411058381199837,  -0.063853621482849121,  -0.078906357288360596,  -0.093014292418956757,   0.018029056489467621,   0.054189603775739670,   0.061475429683923721,   0.175968736410140991,  -0.102474898099899292,
		  -0.125601202249526978,   0.041156720370054245,  -0.001679612905718386,   0.161721259355545044,   0.064328022301197052,   0.188634738326072693,  -0.028993792831897736,   0.126769155263900757,   0.042873885482549667,   0.068185210227966309,  -0.079691641032695770,  -0.084202088415622711,   0.066483087837696075,  -0.098062224686145782,   0.070133052766323090,  -0.099560759961605072,   0.189087554812431335,  -0.070897854864597321,  -0.061707459390163422,  -0.058560900390148163,   0.025347165763378143,   0.044412210583686829,   0.086882062256336212,  -0.130759567022323608,   0.088762424886226654,   0.016021471470594406,  -0.040564578026533127,   0.099894322454929352,   0.141229778528213501,   0.170862540602684021,   0.049272544682025909,  -0.063817590475082397,  -0.143398597836494446,  -0.028007628396153450,  -0.078125804662704468,   0.165800377726554871,  -0.149801731109619141,  -0.092393748462200165,  -0.056583464145660400,  -0.181469619274139404,   0.175515726208686829,  -0.044638566672801971,   0.010492238216102123,  -0.020491868257522583,  -0.024838706478476524,  -0.094835512340068817,  -0.088542155921459198,  -0.055827386677265167,
		   0.045678243041038513,   0.106584496796131134,   0.059086304157972336,   0.055836349725723267,   0.018881579861044884,  -0.058265760540962219,  -0.017682438716292381,   0.101939931511878967,  -0.054876245558261871,   0.057006984949111938,   0.069810032844543457,   0.001450377865694463,   0.051564089953899384,  -0.031329791992902756,  -0.045216731727123260,   0.048310488462448120,  -0.176824703812599182,   0.024150982499122620,   0.044096868485212326,  -0.125874727964401245,  -0.142623275518417358,  -0.045848570764064789,  -0.149121314287185669,   0.100620366632938385,  -0.026209414005279541,   0.035772647708654404,  -0.103414960205554962,  -0.125314369797706604,  -0.031404018402099609,  -0.075270056724548340,   0.041784223169088364,  -0.062549620866775513,   0.092296518385410309,  -0.066393271088600159,   0.066676378250122070,   0.008178491145372391,   0.152449682354927063,  -0.058458067476749420,  -0.033411256968975067,  -0.099918477237224579,  -0.137074619531631470,  -0.029169352725148201,  -0.151332885026931763,  -0.116569831967353821,  -0.094608329236507416,  -0.047773040831089020,  -0.150825724005699158,   0.072550423443317413,
		  -0.178132504224777222,  -0.149905592203140259,  -0.025020891800522804,  -0.052570994943380356,  -0.020756276324391365,  -0.140626296401023865,  -0.042860556393861771,   0.106042683124542236,   0.105792298913002014,   0.130885824561119080,   0.072039619088172913,   0.063327357172966003,  -0.084474459290504456,  -0.136177822947502136,   0.061509441584348679,  -0.127460718154907227,   0.053839974105358124,   0.030689641833305359,  -0.101876348257064819,  -0.119805485010147095,  -0.079773135483264923,  -0.071658909320831299,   0.083060927689075470,  -0.076266422867774963,  -0.051578599959611893,  -0.149022787809371948,   0.067044481635093689,  -0.137494400143623352,  -0.124949671328067780,   0.015371471643447876,  -0.130700439214706421,  -0.067010737955570221,  -0.057790979743003845,  -0.154663696885108948,  -0.078724093735218048,  -0.115016043186187744,  -0.084454677999019623,   0.066207386553287506,   0.029512492939829826,   0.117424331605434418,  -0.146645486354827881,   0.080992497503757477,   0.004670504946261644,   0.112779572606086731,   0.040695317089557648,   0.028808239847421646,   0.016143769025802612,  -0.002119199372828007,
	};
	static const double bias05[]=
	{
		   0.091837942600250244,   0.018493644893169403,  -0.044381678104400635,  -0.161303818225860596,   0.054912261664867401,   0.012351632118225098,  -0.040440678596496582,  -0.072453923523426056,  -0.033789835870265961,  -0.112131863832473755,  -0.190512239933013916,   0.167947039008140564,  -0.100521005690097809,   0.082151308655738831,  -0.068548835813999176,   0.070039950311183929,   0.039981111884117126,   0.097998611629009247,  -0.108500614762306213,   0.164025068283081055,  -0.075272828340530396,  -0.012734623625874519,  -0.066976718604564667,  -0.163080602884292603,
	};
	static const double weight06[]=
	{
		  -0.060736242681741714,   0.153323277831077576,   0.122839696705341339,  -0.112641021609306335,   0.011321227066218853,  -0.165138378739356995,   0.075010024011135101,   0.099530272185802460,   0.194458737969398499,  -0.010795828886330128,   0.128310844302177429,  -0.055551458150148392,  -0.169627666473388672,   0.150400236248970032,   0.091538339853286743,   0.004510958679020405,   0.113178886473178864,   0.224887087941169739,  -0.039828550070524216,   0.076019227504730225,   0.169811472296714783,  -0.071466416120529175,  -0.089742965996265411,  -0.120014950633049011,
		  -0.001883628196083009,   0.127171814441680908,  -0.022682525217533112,  -0.057474836707115173,  -0.203538089990615845,  -0.049043186008930206,   0.227879583835601807,   0.187988117337226868,  -0.182283997535705566,   0.162554964423179626,   0.128188356757164001,   0.059471074491739273,   0.211086705327033997,  -0.167447954416275024,  -0.193478077650070190,  -0.134302198886871338,   0.112401232123374939,   0.145666167140007019,   0.171184524893760681,   0.044843923300504684,  -0.044107090681791306,  -0.061666283756494522,  -0.175103485584259033,  -0.220827966928482056,
		  -0.009446557611227036,  -0.128526926040649414,  -0.182158753275871277,   0.105544894933700562,  -0.023830182850360870,  -0.220201104879379272,  -0.143511399626731873,  -0.183072417974472046,   0.184226900339126587,  -0.059212040156126022,   0.085451804101467133,   0.223035618662834167,   0.039745677262544632,  -0.003985065966844559,  -0.125120699405670166,  -0.117279559373855591,  -0.074052944779396057,   0.053240731358528137,   0.001879858318716288,   0.198097422719001770,  -0.048692755401134491,  -0.018224447965621948,  -0.135852709412574768,   0.030071953311562538,
		   0.079499892890453339,  -0.143787771463394165,  -0.094907954335212708,  -0.045048952102661133,   0.091880589723587036,  -0.045944310724735260,   0.065075524151325226,  -0.101837977766990662,  -0.117667712271213531,  -0.062454756349325180,  -0.224269449710845947,   0.014090085402131081,  -0.190122991800308228,  -0.065475180745124817,   0.074476987123489380,   0.160296037793159485,   0.229693114757537842,   0.011937637813389301,   0.010383751243352890,   0.036903820931911469,  -0.016133036464452744,  -0.175119608640670776,  -0.122551001608371735,  -0.154621541500091553,
		   0.230703860521316528,  -0.033578835427761078,   0.028684398159384727,  -0.109209641814231873,  -0.156115159392356873,   0.061241209506988525,  -0.031855493783950806,   0.116050422191619873,  -0.054998099803924561,  -0.154515266418457031,  -0.177275031805038452,   0.119638144969940186,  -0.110286936163902283,   0.030766157433390617,  -0.095443420112133026,  -0.016002990305423737,  -0.161034494638442993,  -0.118023060262203217,  -0.166711762547492981,   0.268218368291854858,  -0.008552840910851955,  -0.060334589332342148,  -0.050607461482286453,   0.029333040118217468,
		  -0.170682370662689209,   0.220214799046516418,  -0.081448510289192200,  -0.142466768622398376,  -0.204487711191177368,  -0.103795707225799561,   0.032956287264823914,   0.131082102656364441,  -0.180978700518608093,   0.012058815918862820,  -0.236895814538002014,   0.083818748593330383,  -0.122449010610580444,   0.064453579485416412,  -0.148628070950508118,  -0.173870712518692017,  -0.079489640891551971,   0.141383931040763855,  -0.076363191008567810,  -0.008743142709136009,   0.075082801282405853,   0.139408215880393982,   0.031130045652389526,  -0.197806745767593384,
		  -0.117607317864894867,  -0.085174091160297394,   0.177664160728454590,   0.109899871051311493,   0.226252093911170959,  -0.094461157917976379,   0.157153695821762085,   0.070472680032253265,  -0.079088449478149414,   0.045136738568544388,   0.047609880566596985,  -0.241935610771179199,   0.057852499186992645,  -0.169775187969207764,   0.143618270754814148,   0.018485533073544502,   0.196334168314933777,   0.245058611035346985,   0.233543619513511658,   0.110469155013561249,   0.164045289158821106,  -0.119766242802143097,   0.202312186360359192,   0.133275911211967468,
		   0.109431460499763489,  -0.005984063260257244,   0.190501421689987183,   0.087689474225044250,   0.211055338382720947,  -0.183970719575881958,  -0.155309140682220459,   0.054633472114801407,   0.146974757313728333,  -0.131211385130882263,  -0.118153117597103119,  -0.087778382003307343,   0.075258299708366394,  -0.201875686645507812,   0.198763430118560791,   0.001793141593225300,   0.226198241114616394,   0.084301173686981201,   0.018759205937385559,  -0.051629010587930679,   0.163324236869812012,   0.102497510612010956,  -0.120703116059303284,  -0.189795598387718201,
		   0.131256058812141418,  -0.068303972482681274,  -0.220242008566856384,  -0.081497631967067719,  -0.015013623051345348,  -0.020717598497867584,   0.156696662306785583,   0.228843420743942261,  -0.174195200204849243,   0.121696397662162781,   0.052356924861669540,  -0.042450822889804840,   0.230765089392662048,  -0.033159524202346802,  -0.237592965364456177,  -0.094589039683341980,   0.158853068947792053,  -0.016485735774040222,   0.141790434718132019,  -0.001315505127422512,  -0.166180655360221863,   0.089876122772693634,  -0.006825941614806652,   0.045443415641784668,
		  -0.174149662256240845,   0.170877307653427124,  -0.196628585457801819,  -0.053345352411270142,   0.069475188851356506,  -0.185161814093589783,   0.172141745686531067,   0.020004024729132652,  -0.166785240173339844,  -0.184789821505546570,  -0.235642895102500916,   0.230764314532279968,   0.203990325331687927,  -0.178132236003875732,  -0.101202957332134247,  -0.065075010061264038,   0.031006133183836937,  -0.120820678770542145,  -0.181657597422599792,   0.165842965245246887,  -0.134793609380722046,   0.192087799310684204,  -0.206468626856803894,  -0.129970207810401917,
		  -0.067387424409389496,   0.018011067062616348,   0.107435025274753571,   0.099597699940204620,  -0.011602890677750111,  -0.010809930972754955,   0.097837977111339569,   0.038563929498195648,  -0.053022760897874832,  -0.057825095951557159,  -0.067921906709671021,  -0.060184504836797714,  -0.135116726160049438,   0.143288984894752502,   0.056910578161478043,  -0.028575902804732323,  -0.212423235177993774,   0.031316764652729034,   0.047466490417718887,   0.154441475868225098,  -0.012315839529037476,  -0.111869812011718750,  -0.151728659868240356,   0.169047653675079346,
		  -0.082176581025123596,  -0.121721178293228149,  -0.108682073652744293,  -0.161873340606689453,  -0.233192324638366699,  -0.209179610013961792,   0.243606567382812500,   0.212962016463279724,  -0.119046688079833984,  -0.175966247916221619,   0.085670597851276398,   0.238510906696319580,   0.058666978031396866,   0.049741830676794052,  -0.101705752313137054,  -0.219841927289962769,   0.022478960454463959,   0.051214564591646194,   0.065815150737762451,   0.012016143649816513,  -0.124207139015197754,   0.043554283678531647,   0.158288121223449707,  -0.017975628376007080,
	};
	static const double bias06[]=
	{
		  -0.207499399781227112,  -0.075083129107952118,  -0.115031316876411438,  -0.040588956326246262,   0.124039337038993835,   0.095791935920715332,  -0.089956320822238922,   0.093130916357040405,  -0.146124675869941711,   0.168467402458190918,  -0.032112643122673035,  -0.047444786876440048,
	};
	static const double weight07[]=
	{
		  -0.229891389608383179,   0.278955668210983276,   0.293872445821762085,  -0.251623719930648804,   0.073358006775379181,   0.167853638529777527,  -0.082937791943550110,  -0.079294085502624512,   0.095349296927452087,   0.287303596735000610,   0.124519132077693939,   0.148270815610885620,
	};
	static const double bias07[]=
	{
		  -0.138395205140113831,
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

	double avpred[3]={0};//

	int idx, idx2, pred;
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	memcpy(dst, src, (size_t)res<<2);
	for(int kc=0;kc<3;++kc)//channel loop
	{
		for(int ky=0;ky<ih;++ky)//y loop
		{
			for(int kx=0;kx<iw;++kx)//x loop
			{
				//if(kc==0&&kx==3&&ky==3)
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
							nb[idx+C03_NNB]=(unsigned)(ky+ky2-C03_REACH)<(unsigned)(ih-C03_REACH)&&(unsigned)(kx+kx2-C03_REACH)<(unsigned)(iw-C03_REACH*2)?C03_CVT_PIXEL(errors[idx2]):0;
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
						nb[idx+C03_NNB]=(unsigned)(ky-C03_REACH)<(unsigned)(ih-C03_REACH)&&(unsigned)(kx+kx2-C03_REACH)<(unsigned)(iw-C03_REACH*2)?C03_CVT_PIXEL(errors[idx2]):0;
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

					if(k+1<_countof(nnch))
						c03_activation(vo, co);
					
					C03_DATATYPE *t2;
					SWAPVAR(vi, vo, t2);
				}
				pred=C03_CVT_PRED(*vi);
				//int pred2=(int)(((*vi)>(-1) ? (*vi)<(1) ? (*vi) : (1) : (-1)) * 128.f + 0.5f);

				avpred[kc]+=fabs(*vi*128);//

				idx=(iw*ky+kx)<<2|kc;
				pred^=-fwd;
				pred+=fwd;
				dst[idx]=src[idx]+pred;
			}//x loop
		}//y loop
	}//channel loop
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
	for(int k=0;k<_countof(src_weights);++k)
		_mm_free(coeffs[k]);

	set_window_title("%lf %lf %lf", avpred[0]/res, avpred[1]/res, avpred[2]/res);//
}