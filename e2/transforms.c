#include"e2.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#ifdef __AVX__
#include<immintrin.h>
#else
#include<tmmintrin.h>
#endif
#include<Windows.h>//threads
#include<process.h>
static const char file[]=__FILE__;

void calc_histogram(const unsigned char *buf, ptrdiff_t bytesize, ptrdiff_t stride, int *hist)
{
	memset(hist, 0, 256*sizeof(int));
	for(ptrdiff_t k=0, end=bytesize-(stride-1);k<end;k+=stride)
		++hist[buf[k]];
}

static void print_sbuf(short *buf, int bw, int bh, int kc, int bytestride)
{
	static int call=0;
	printf("call %d ch %d\n", call, kc);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%d\t", buf[bytestride*(bw*ky+kx)+kc]);
		printf("\n");
	}
	printf("\n");
	++call;
}
static void print_fbuf(float *buf, int bw, int bh, int nplaces)
{
	static int call=0;
	printf("call %d\n", call);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%.*f\t", nplaces, buf[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
	++call;
}
void compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward)
{
	ptrdiff_t len=(ptrdiff_t)bytestride*iw*ih;
	int inc=bytestride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-bytestride:0;k>=0&&k<len;k+=inc)
	{
		//{//
		//	ptrdiff_t idx=k>>2, kx=idx%1920, ky=idx/1920;
		//	if(kx==5&&ky==1)
		//		kx=5;
		//}//
		if(memcmp(b1+k, b0+k, symbytes))
		{
			ptrdiff_t idx=k/bytestride, kx=idx%iw, ky=idx/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
			for(int kc=0;kc<symbytes;++kc)
				printf("C%d  0x%02X != 0x%02X    %d != %d\n", kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]);
			//if(backward)
			//	printf("%s error at %d - %d: 0x%02X != 0x%02X\n", name, (int)len-1, (int)(len-1-k), b1[k], b0[k]);
			//else
			//	printf("%s error at %d: 0x%02X != 0x%02X\n", name, (int)k, b1[k], b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}
void compare_bufs_ps(float *b1, float *b0, int iw, int ih, const char *name, int backward)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	int inc=1-(backward<<1);
	for(ptrdiff_t k=backward?res-1:0;k>=0&&k<res;k+=inc)
	{
		if(fabsf(b1[k]-b0[k])>1e-6)
		{
			ptrdiff_t kx=k%iw, ky=k/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0  %f != %f  error=%f\n", name, kx, ky, iw, ih, b1[k], b0[k], b1[k]-b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}


short *g_param_ptr=0;
void apply_transforms_fwd(unsigned char *buf, int bw, int bh)
{
	//ArrayHandle sizes=dwt2d_gensizes(bw, bh, 7, 7, 0);
	//unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));//for DWT

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char
	
	colortransform_ycocb_fwd((char*)buf, bw, bh);
	//colortransform_ycocg_fwd((char*)buf, bw, bh);
	//colortransform_xgz_fwd((char*)buf, bw, bh);
	//colortransform_xyz_fwd((char*)buf, bw, bh);
	
	//pred_opt_apply((char*)buf, bw, bh, 1);
	//pred_w2_apply((char*)buf, bw, bh, pw2_params, 1);
	pred_jxl_apply((char*)buf, bw, bh, g_param_ptr?g_param_ptr:jxlparams_i16, 1);

	//pred_jxl((char*)buf, bw, bh, 0, 4, jxlpred_params   , 1);
	//pred_jxl((char*)buf, bw, bh, 1, 4, jxlpred_params+11, 1);
	//pred_jxl((char*)buf, bw, bh, 2, 4, jxlpred_params+22, 1);

	//pred_grad_fwd((char*)buf, bw, bh, 3, 4);
	//image_differentiate((char*)buf, bw, bh, 3, 4);
	//for(int kc=0;kc<3;++kc)
	//	//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_squeeze_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

	addbuf(buf, bw, bh, 3, 4, 128);
	
	//save_DWT_int8_all("kodim21-cdf53", buf, (DWTSize*)sizes->data, (int)sizes->count);
	//save_DWT_int8("kodim21-squeeze-stage", buf, (DWTSize*)sizes->data, (int)sizes->count, 4);//
	//save_DWT_int8("kodim21-cubic-stage", buf, (DWTSize*)sizes->data, (int)sizes->count, 4);//

	//free(temp);
	//array_free(&sizes);
}
void apply_transforms_inv(unsigned char *buf, int bw, int bh)
{
	//ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	//unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	
	addbuf(buf, bw, bh, 3, 4, 128);
	
	//pred_opt_apply((char*)buf, bw, bh, 0);
	//pred_w2_apply((char*)buf, bw, bh, pw2_params, 0);
	pred_jxl_apply((char*)buf, bw, bh, g_param_ptr?g_param_ptr:jxlparams_i16, 0);

	//pred_jxl((char*)buf, bw, bh, 0, 4, jxlpred_params   , 0);
	//pred_jxl((char*)buf, bw, bh, 1, 4, jxlpred_params+11, 0);
	//pred_jxl((char*)buf, bw, bh, 2, 4, jxlpred_params+22, 0);

	//pred_grad_inv((char*)buf, bw, bh, 3, 4);
	//image_integrate((char*)buf, bw, bh, 3, 4);
	//for(int kc=0;kc<3;++kc)
	//	//dwt2d_haar_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	dwt2d_squeeze_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf53_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf97_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	
	colortransform_ycocb_inv((char*)buf, bw, bh);
	//colortransform_ycocg_inv((char*)buf, bw, bh);
	//colortransform_xgz_inv((char*)buf, bw, bh);
	//colortransform_xyz_inv((char*)buf, bw, bh);

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char

	//free(temp);
	//array_free(&sizes);
}
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
	}
}




void colortransform_ycocg_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=b;
		b+=r>>1;
		g-=b;
		b+=g>>1;

		buf[k  ]=r;//Co
		buf[k|1]=g;//Cg
		buf[k|2]=b;//Y
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

		r-=g;
		g+=r>>1;
		b-=g;
		g+=b>>1;

		//buf[k  ]=g;//Y
		//buf[k|1]=r;//Co
		//buf[k|2]=b;//Cb

		//buf[k  ]=r;//Co
		//buf[k|1]=b;//Cb
		//buf[k|2]=g;//Y

		buf[k  ]=r;//Co	original
		buf[k|1]=g;//Y
		buf[k|2]=b;//Cb
	}
}
void colortransform_ycocb_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		//char r=buf[k|1], g=buf[k], b=buf[k|2];//Y Co Cb
		//char r=buf[k], g=buf[k|2], b=buf[k|1];//Co Cb Y
		char r=buf[k], g=buf[k|1], b=buf[k|2];//Co Y Cb	original
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
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
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub;

				if(topleft>MAXVAR(left, top))//planar prediction (unplane)
					sub=MINVAR(left, top);
				else if(topleft<MINVAR(left, top))
					sub=MAXVAR(left, top);
				else
					sub=left+top-topleft;

				buf[idx]-=sub;
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
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub;

				if(topleft>MAXVAR(left, top))
					sub=MINVAR(left, top);
				else if(topleft<MINVAR(left, top))
					sub=MAXVAR(left, top);
				else
					sub=left+top-topleft;

				buf[idx]+=sub;
			}
		}
	}
}


//	#define ONLY_USE_W2PRED
	#define ONLY_USE_JXLPRED

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
				pred=predictions[8];

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
void pred_opt_apply(char *buf, int iw, int ih, int fwd)
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
	pred_jxl_prealloc(buf, iw, ih, 2, jxlparams_i16+11*2   , fwd, buf2, temp);

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
void pred_w2_opt_v2(const char *buf2, int iw, int ih, short *params, int loud)
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
	//char title0[256];
	//get_window_title(title0, 256);
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

				//set_window_title("Ch%d csize %lf [%d/%d %d/%d]...", kc, csize0, kc*COUNTOF(steps)+ks+1, COUNTOF(steps)*3, idx+1, PW2_NPARAM);//
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
	//set_window_title("%s", title0);//
}

double pred_opt_calcloss(const char *src, int iw, int ih, int kc, short *params, int *temp, char *dst, int *hist)
{
	int res=iw*ih;
#ifdef ONLY_USE_W2PRED
	pred_w2_prealloc(src, iw, ih, kc, params, 1, dst, temp);
#elif defined ONLY_USE_JXLPRED
	pred_jxl_prealloc(src, iw, ih, kc, params, 1, dst, temp);
#else
	switch(kc)
	{
	case 0:
	case 2:
		pred_jxl_prealloc(src, iw, ih, kc, params, 1, dst, temp);
		break;
	case 1:
		pred_w2_prealloc(src, iw, ih, kc, params, 1, dst, temp);
		break;
	}
#endif
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

void pred_opt_printrow(const short *p, int count)
{
	for(int kp=0;kp<count;++kp)
	{
		short val=p[kp];
		printf(" %c0x%04X,", val<0?'-':' ', abs(val));
	}
	printf("\n");
}
void pred_opt_printparam()
{
	const short *params=jxlparams_i16;
	pred_opt_printrow(params, 4);
	pred_opt_printrow(params+4, 7);
	printf("\n");

#ifdef ONLY_USE_JXLPRED
	params=jxlparams_i16+11;
	pred_opt_printrow(params, 4);
	pred_opt_printrow(params+4, 7);
#else
	params=pw2_params+PW2_NPARAM;
	pred_opt_printrow(params, PW2_NPRED);
	pred_opt_printrow(params+PW2_NPRED, PW2_NPARAM-PW2_NPRED);
#endif
	printf("\n");

	params=jxlparams_i16+22;
	pred_opt_printrow(params, 4);
	pred_opt_printrow(params+4, 7);
}

#define O6_MAXTHREADS (PW2_NPARAM<<1)
typedef struct Opt6ContextStruct
{
	const char *src;
	int iw, ih, kc;
	int *temp;
	char *dst;
	int *hist;
	int loud;

	int threadidx;
	
	double loss;
	const short *params;
	int nparam;
	int step;
} Opt6Context;
Opt6Context o6_ctx[O6_MAXTHREADS];
static unsigned __stdcall o6_thread(void *args)
{
	Opt6Context *ctx=(Opt6Context*)args;
	short params[MAXVAR(11, PW2_NPARAM)];
	int idx=ctx->threadidx>>1, sign=ctx->threadidx&1;

	if(!ctx->src)
	{
		LOG_ERROR("Invalid thread");
		return 0;
	}

	memcpy(params, ctx->params, ctx->nparam*sizeof(short));
	if(sign)
		params[idx]-=ctx->step;
	else
		params[idx]+=ctx->step;
	ctx->loss=pred_opt_calcloss(ctx->src, ctx->iw, ctx->ih, ctx->kc, params, ctx->temp, ctx->dst, ctx->hist);
	return 0;
}
void pred_opt_opt_v6(const char *buf2, int iw, int ih, int loud)//multi-threaded
{
	double t_start=time_ms();
	double loss=0;
	int res=iw*ih;
	char *buf3=(char*)malloc((O6_MAXTHREADS+1)*(size_t)res<<2);
	int *temp=(int*)malloc((O6_MAXTHREADS+1)*(size_t)iw*(PW2_NPRED+1)*2*sizeof(int));
	int *hist=(int*)malloc((O6_MAXTHREADS+1)*256*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	size_t temp_idx=0;
	size_t buf3_idx=0, hist_idx=0;
	{
		for(int kt=0;kt<O6_MAXTHREADS;++kt)
		{
			Opt6Context *p=o6_ctx+kt;
			p->src=buf2;
			p->iw=iw;
			p->ih=ih;
			p->kc=0;
			p->temp=temp+temp_idx;
			p->dst=buf3+buf3_idx;
			p->hist=hist+hist_idx;
			p->loud=loud;

			p->threadidx=kt;

			p->loss=0;

			temp_idx+=(size_t)iw*(PW2_NPRED+1)*2;
			buf3_idx+=(size_t)res<<2;
			hist_idx+=256;
		}
	}
	
	int steps[]={256, 128, 64, 32, 16, 8, 4, 2, 1};
	void *threads[O6_MAXTHREADS];
	for(int kc=0;kc<3;++kc)
	{
		short *params=0, nparam=0;
#ifdef ONLY_USE_W2PRED
		switch(kc)
		{
		case 0:params=pw2_params,              nparam=PW2_NPARAM;break;
		case 1:params=pw2_params+PW2_NPARAM,   nparam=PW2_NPARAM;break;
		case 2:params=pw2_params+PW2_NPARAM*2, nparam=PW2_NPARAM;break;
		}
#elif defined ONLY_USE_JXLPRED
		switch(kc)
		{
		case 0:params=jxlparams_i16,    nparam=11;break;
		case 1:params=jxlparams_i16+11, nparam=11;break;
		case 2:params=jxlparams_i16+22, nparam=11;break;
		}
#else
		switch(kc)
		{
		case 0:params=jxlparams_i16,         nparam=11;        break;
		case 1:params=pw2_params+PW2_NPARAM, nparam=PW2_NPARAM;break;
		case 2:params=jxlparams_i16+22,      nparam=11;        break;
		}
#endif
		int nthreads=nparam<<1;
		const int nsteps=COUNTOF(steps), ntrials=nsteps*4;
		for(int ks=0;ks<ntrials;++ks)
		{
			int step=steps[ks%nsteps];
			for(int kt=0;kt<nthreads;++kt)
			{
				Opt6Context *p=o6_ctx+kt;
				p->kc=kc;
				p->params=params;
				p->nparam=nparam;
				p->step=step;
				threads[kt]=(void*)_beginthreadex(0, 0, o6_thread, p, 0, 0);
				if(!threads[kt])
				{
					LOG_ERROR("Allocation error");
					return;
				}
			}
			double csize0=pred_opt_calcloss(buf2, iw, ih, kc, params, temp+temp_idx, buf3+buf3_idx, hist+hist_idx);
			WaitForMultipleObjects(nthreads, threads, TRUE, INFINITE);
			for(int kt=0;kt<nthreads;++kt)
				CloseHandle(threads[kt]);

			double bestresult=csize0;
			int bestthread=nthreads;
			for(int kt=0;kt<nthreads;++kt)
			{
				Opt6Context *p=o6_ctx+kt;
				if(bestresult>p->loss)
					bestresult=p->loss, bestthread=kt;
			}
			if(bestthread<nthreads)
			{
				int idx=bestthread>>1, sign=bestthread&1;
				if(sign)
					params[idx]-=step;
				else
					params[idx]+=step;
			}
			printf("[%d/3 %2d/%2d] C%d csize %14lf CR %9lf\r", kc+1, ks+1, ntrials, kc, bestresult, iw*ih/bestresult);//
			//set_window_title("Ch%d csize %lf [%d/%d]...", kc, csize0, kc*COUNTOF(steps)+ks+1, COUNTOF(steps)*3);//
		}
		if(loud)
			printf("\n");//
	}

	if(loud)
	{
		printf("Pred opt elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}
	free(hist);
	free(temp);
	free(buf3);
}


short jxlparams_i16[33]=//signed fixed 7.8 bit
{
	//CLIC	X
	//0x0779,  0x120C,  0x0C00,  0x0A3D, -0x0051,  0x000E, -0x0196, -0x01E7,  0x0065,  0x0076,  0x00FA,
	//0x0897,  0x1A0F,  0x1697,  0x0A6B, -0x0054, -0x002B, -0x0021, -0x0051,  0x00AD, -0x001C, -0x0028,
	//0x0BB2,  0x186F,  0x1012,  0x0B8F, -0x0053, -0x005F, -0x008B, -0x0085,  0x0006,  0x009D,  0x008C,

	//kodak
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,

	//0x0BD6, 0x10E7, 0x11F9, 0x0BEC,    0x0003,  0x0016, -0x01AC, -0x0127, -0x00C3, -0x0048,  0x00C3,
	//0x0DB8, 0x0E22, 0x181F, 0x0BF3,   -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	//0x0684, 0x0F1F, 0x100D, 0x0C40,   -0x0007, -0x000D, -0x0087, -0x0063, -0x009F, -0x0018,  0x00F2,

	//0x0EEE, 0x103B, 0x0FAC, 0x0FF6, -0x0005,  0x0006, -0x0125, -0x007A, -0x004C,  0x0013,  0x0032,
	//0x0F92, 0x0FF3, 0x0FDF, 0x1006, -0x007A, -0x00C2, -0x000B, -0x006E,  0x007F, -0x004E, -0x00A0,
	//0x0E7E, 0x0FD7, 0x1024, 0x0FF4, -0x0011, -0x0048, -0x0054, -0x0053, -0x0038,  0x0018,  0x007D,

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
void   pred_jxl_prealloc(const char *src, int iw, int ih, int kc, const short *params, int fwd, char *dst, int *temp_w10)
{
	int res=iw*ih, errorbuflen=iw*2, rowlen=iw<<2;
	int *error=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int idx=kc;
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred, curr;
			
			const char *src2=fwd?src:dst;
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
				esumtopleft=etop+eleft;
			int predictions[]=//fixed 23.8 bit
			{
				(cleft+ctopright-ctop)<<8,
				(ctop<<8)-((esumtopleft+etopright)*params[4]>>8),
				(cleft<<8)-((esumtopleft+etopleft)*params[5]>>8),
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
double pred_jxl_calcloss(const char *src, int iw, int ih, int kc, const short *params, int *temp, char *dst, int *hist)
{
	int res=iw*ih;
	pred_jxl_prealloc(src, iw, ih, kc, params, 1, dst, temp);
	//addbuf(dst+kc, iw, ih, 1, 4, 128);//just slows down the optimization
	calc_histogram(dst+kc, (ptrdiff_t)res<<2, 4, hist);

	double entropy=0;
	int freq;
	double p;
	for(int k=0;k<256;++k)
	{
		freq=hist[k];
		if(freq)
		{
			p=(double)freq/res;
			p*=0x10000-255;
			++p;
			p/=0x10000;
			entropy-=p*log2(p);
		}
	}
	double invCR=entropy/8, csize=res*invCR;
	//if(loud)
	//	printf("%4d %14lf\r", it, csize);
	return csize;
}
void   pred_jxl_opt_v2(const char *buf2, int iw, int ih, short *params, int loud)
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
void   pred_jxl_apply(char *buf, int iw, int ih, short *allparams, int fwd)
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