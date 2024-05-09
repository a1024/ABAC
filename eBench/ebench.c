#include"ebench.h"
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
#include"stb_image.h"
#include"lodepng.h"
#define DEBUG_MEMORY_IMPLEMENTATION
#include"intercept_malloc.h"
static const char file[]=__FILE__;

float
	mouse_sensitivity=0.003f,
	key_turn_speed=0.03f;
Camera cam=
{
	10, 10, 10,
	(float)(225*M_PI/180), (float)(324.7356103172454f*M_PI/180),
	1,
	0.04f, (float)(2*M_PI/180),
}, cam0;

ArrayHandle fn=0;
size_t filesize=0;

Image *im0=0, *im1=0;//im0: original image, im1: image with selected transforms
int pred_ma_enabled=1;//modular arithmetic for spatial predictors
int separate_grayscale=1;//separate channels are shown greyscale
unsigned txid_separate_r=0, txid_separate_g=0, txid_separate_b=0, txid_separate_a=0;
unsigned char *im_export=0, *zimage=0;

unsigned txid_jointhist[256]={0};

typedef enum VisModeEnum
{
	//VIS_PLANES,
	//VIS_MESH,

//	VIS_MESH_SEPARATE,
	VIS_IMAGE_TRICOLOR,
	VIS_IMAGE,
	VIS_HISTOGRAM,
	VIS_JOINT_HISTOGRAM,
	VIS_ZIPF,
	//VIS_BAYES,
	//VIS_IMAGE_BLOCK,
	//VIS_DWT_BLOCK,

	VIS_COUNT,
} VisMode;
int mode=VIS_IMAGE;

typedef enum TransformTypeEnum
{
	CT_FWD_SUBGREEN,	CT_INV_SUBGREEN,
	CT_FWD_JPEG2000,	CT_INV_JPEG2000,//	(1997) JPEG2000 RCT
	CT_FWD_YCoCg_R,		CT_INV_YCoCg_R,	//	(2003) AVC, HEVC, VVC
	CT_FWD_YCbCr_R_v1,	CT_INV_YCbCr_R_v1,
	CT_FWD_YCbCr_R_v2,	CT_INV_YCbCr_R_v2,
	CT_FWD_YCbCr_R_v3,	CT_INV_YCbCr_R_v3,
	CT_FWD_YCbCr_R_v4,	CT_INV_YCbCr_R_v4,
	CT_FWD_YCbCr_R_v5,	CT_INV_YCbCr_R_v5,
	CT_FWD_YCbCr_R_v6,	CT_INV_YCbCr_R_v6,
	CT_FWD_YCbCr_R_v7,	CT_INV_YCbCr_R_v7,
	CT_FWD_Pei09,		CT_INV_Pei09,
//	CT_FWD_CrCgCb,		CT_INV_CrCgCb,
//	CT_FWD_YRGB_v1,		CT_INV_YRGB_v1,
//	CT_FWD_YRGB_v2,		CT_INV_YRGB_v2,
//	CT_FWD_CMYK,		CT_INV_CMYK_DUMMY,
//	CT_FWD_MATRIX,		CT_INV_MATRIX,
	CT_FWD_YCbCr,		CT_INV_YCbCr,	//LOSSY	JPEG
	CT_FWD_XYB,		CT_INV_XYB,	//LOSSY	(2021) JPEG XL
	CT_FWD_CUSTOM,		CT_INV_CUSTOM,
	CT_FWD_ADAPTIVE,	CT_INV_ADAPTIVE,

	CST_FWD_SEPARATOR,	CST_INV_SEPARATOR,
	
	ST_FWD_PACKSIGN,	ST_INV_PACKSIGN,
	ST_FWD_T47,		ST_INV_T47,
	ST_FWD_P3,		ST_INV_P3,
	ST_FWD_G2,		ST_INV_G2,
	ST_FWD_MM,		ST_INV_MM,
	ST_FWD_WP,		ST_INV_WP,
	ST_FWD_WPU,		ST_INV_WPU,
	ST_FWD_DEFERRED,	ST_INV_DEFERRED,
	ST_FWD_CUSTOM3,		ST_INV_CUSTOM3,
	ST_FWD_CUSTOM,		ST_INV_CUSTOM,
//	ST_FWD_NBLIC,		ST_INV_NBLIC,
	ST_FWD_CALIC,		ST_INV_CALIC,
	ST_FWD_OLS,		ST_INV_OLS,
	ST_FWD_OLS2,		ST_INV_OLS2,
	ST_FWD_OLS3,		ST_INV_OLS3,
	ST_FWD_OLS4,		ST_INV_OLS4,
	ST_FWD_OLS5,		ST_INV_OLS5,
	ST_FWD_PU,		ST_INV_PU,
	ST_FWD_CG3D,		ST_INV_CG3D,
	ST_FWD_CLAMPGRAD,	ST_INV_CLAMPGRAD,
	ST_FWD_AV2,		ST_INV_AV2,
	ST_FWD_MEDIAN,		ST_INV_MEDIAN,
//	ST_FWD_ECOEFF,		ST_INV_ECOEFF,
//	ST_FWD_AVERAGE,		ST_INV_AVERAGE,
	ST_FWD_MULTISTAGE,	ST_INV_MULTISTAGE,
//	ST_FWD_ZIPPER,		ST_INV_ZIPPER,
//	ST_FWD_DIR,		ST_INV_DIR,
#if 0
	ST_FWD_CTX,		ST_INV_CTX,
	ST_FWD_C20,		ST_INV_C20,
//	ST_FWD_WU97,		ST_INV_WU97,
	ST_FWD_SHUFFLE,		ST_INV_SHUFFLE,
//	ST_FWD_SPLIT,		ST_INV_SPLIT,
#endif
	
//	ST_FWD_LAZY,		ST_INV_LAZY,
	ST_FWD_HAAR,		ST_INV_HAAR,
	ST_FWD_SQUEEZE,		ST_INV_SQUEEZE,
	ST_FWD_LEGALL53,	ST_INV_LEGALL53,
	ST_FWD_CDF97,		ST_INV_CDF97,
//	ST_FWD_EXPDWT,		ST_INV_EXPDWT,
//	ST_FWD_CUSTOM_DWT,	ST_INV_CUSTOM_DWT,

	ST_FWD_DCT4,		ST_INV_DCT4,
	ST_FWD_DCT8,		ST_INV_DCT8,

	T_COUNT,
} TransformType;
int transforms_customenabled=0;
char transforms_mask[T_COUNT]={0};
ArrayHandle transforms=0;//array of chars
float guizoom=1.25f;

double av_rmse=0, g_lr=1e-10;

int profile_idx=0;
double minloss=0, maxloss=0;

int blocksize=16, margin=32;
int blockmx=0, blockmy=0;

float
	pixel_amplitude=10,//4
	mesh_separation=100;//10
int extrainfo=0;

ArrayHandle cpu_vertices=0;
unsigned gpu_vertices=0;

ArrayHandle jointhist=0;
int jointhist_nbits=6;//max
int jhx=0, jhy=0;
double ch_entropy[4]={0};//RGBA/YUVA
//int usage[4]={0};
EContext ec_method=ECTX_HIST;//ECTX_MIN_QN_QW;
int ec_adaptive=0, ec_adaptive_threshold=3200, ec_expbits=5, ec_msb=2, ec_lsb=0;

#define combCRhist_SIZE 128
#define combCRhist_logDX 2
float combCRhist[combCRhist_SIZE][4]={0}, combCRhist_max=1;
int combCRhist_idx=0;

int show_full_image=0;
int space_not_color=0;

//joint hist box contour
float jh_cubesize=64;
int jhc_boxdx=64, jhc_boxdy=64;
int jhc_xbox=0, jhc_ybox=0;//[0, screen_dim-1]
ArrayHandle jhc_mesh=0;
unsigned jhc_gpubuf=0;
float jhc_level=10.5f;
int loud_transforms=1;
//int crop_enable=0, crop_x=0, crop_y=0, crop_dx=128, crop_dy=128;

void calc_csize_ans_separate(Image const *image, size_t *csizes)
{
	if(!csizes)
		return;
	int maxdepth=calc_maxdepth(image, 0), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc((maxlevels+1)*sizeof(int));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	int res=image->iw*image->ih;
	unsigned state=0x10000;
	for(int kc=0;kc<3;++kc)//naive way: one dedicated hist for each channel
	{
		int depth=image->depth[kc], nlevels=1<<depth;
		memset(hist, 0, nlevels*sizeof(int));
		for(int k=0;k<res;++k)//calc histogram
		{
			int sym=image->data[k<<2|kc]+(nlevels>>1);
			if((unsigned)sym>=(unsigned)nlevels)
				LOG_ERROR("Symbol OOB");
			++hist[sym];
		}
		int nusedlevels=0;
		for(int ks=0;ks<nlevels;++ks)
			nusedlevels+=hist[ks]!=0;
		int sum=0;
		for(int ks=0, ks2=0;ks<nlevels;++ks)//quantize & accumulate CDF
		{
			int freq=hist[ks];
			hist[ks]=(int)((long long)sum*(0x10000-nusedlevels)/res)+ks2;
			ks2+=freq!=0;
			sum+=freq;
		}
		hist[nlevels]=0x10000;

		size_t csize=0;
		for(int k=0;k<res;++k)//calc csize
		{
			int sym=image->data[k<<2|kc]+(nlevels>>1);
			int cdf=hist[sym], freq=hist[sym+1]-cdf;
			
			if(state>=(unsigned)(freq<<16))//renorm
			{
				csize+=2;
				state>>=16;
			}
			state=state/freq<<16|(cdf+state%freq);//update
		}
		csizes[kc]=csize;
	}
	csizes[3]=csizes[0]+csizes[1]+csizes[2]+4;
	free(hist);
}

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
static void hybriduint_encode(int val, int *tbn, const int *config)
{
	int exp=config[0], msb=config[1], lsb=config[2];
	int token, bypass, nbits;
	val=(val<<1)^-(val<0);//pack sign
	if(val<(1<<exp))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<exp) + ((lgv-exp)<<(msb+lsb)|(mantissa>>(lgv-msb))<<lsb|(mantissa&((1<<lsb)-1)));
		nbits=lgv-(msb+lsb);
		bypass=val>>lsb&((1LL<<nbits)-1);
	}
	tbn[0]=token;
	tbn[1]=bypass;
	tbn[2]=nbits;
}
#if 0
#define HYBRIDUINT_EXP 4//exponent threshold for bypass
#define HYBRIDUINT_MSB 2//number of bits taken by token from left from sign-packed pixel value
#define HYBRIDUINT_LSB 0//number of bits taken by token from right from sign-packed pixel value
void hybriduint_encode(int val, int *tbn)
{
	int token, bypass, nbits;
	val=(val<<1)^-(val<0);//pack sign
	if(val<(1<<HYBRIDUINT_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<HYBRIDUINT_EXP) + ((lgv-HYBRIDUINT_EXP)<<(HYBRIDUINT_MSB+HYBRIDUINT_LSB)|(mantissa>>(lgv-HYBRIDUINT_MSB))<<HYBRIDUINT_LSB|mantissa&((1<<HYBRIDUINT_LSB)-1));
		nbits=lgv-(HYBRIDUINT_MSB+HYBRIDUINT_LSB);
		bypass=val>>HYBRIDUINT_LSB&((1LL<<nbits)-1);
	}
	tbn[0]=token;
	tbn[1]=bypass;
	tbn[2]=nbits;
}
#endif
typedef union TempHybridStruct
{
	struct
	{
		unsigned char token, nbits;
		unsigned short bypass;
	};
	int data;
} TempHybrid;
static const int calcsize_ans_qlevels[]=
{
	-500, -392, -255, -191, -127, -95, -63, -47, -31,
	-23,  -15,  -11,  -7,   -4,   -3,  -1,  0,   1,
	3,    5,    7,    11,   15,   23,  31,  47,  63,
	95,   127,  191,  255,  392,  500
};
//static const int calcsize_ans_qlevels_u[]=
//{
//	0,  1,  3,  5,   7,   11,  15,  23, 31,
//	47, 63, 95, 127, 191, 255, 392, 500
//};
static const int ans_qlevels_u[]=
{
	//0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16, 18, 20, 22, 24, 28, 32, 36, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 384, 512		//v1	opt 2.414476

	0, 1, 3, 5, 7, 9, 11, 13, 15, 18, 23, 31, 47, 63, 95, 127, 191, 255, 392, 500		//v2	opt 2.411191
};
static int get_ctx(Image const *image, int kc, int kx, int ky)
{
#define CTX_REACH 2
	int energy=0;
	//int nb[2*(CTX_REACH+1)*CTX_REACH];
	for(int ky2=-CTX_REACH;ky2<=0;++ky2)
	{
		for(int kx2=-CTX_REACH;kx2<=CTX_REACH;++kx2)
		{
			if(!ky2&&!kx2)
				break;
			if((unsigned)(ky+ky2)<(unsigned)image->ih&&(unsigned)(kx+kx2)<(unsigned)image->iw)
				energy+=abs(image->data[(image->iw*(ky+ky2)+kx+kx2)<<2|kc])/(kx2*kx2+ky2*ky2);//r2 opt 2.390132  r2 opt 2.410421  r3 opt 2.414187  r4 opt 2.414871
				//energy+=abs(image->data[(image->iw*(ky+ky2)+kx+kx2)<<2|kc])/(abs(kx2)+abs(ky2));//r2 opt  2.394843 r4 opt 2.399578
				//energy+=abs(image->data[(image->iw*(ky+ky2)+kx+kx2)<<2|kc])>>(abs(kx2)+abs(ky2));//r2 opt 2.390961
				//energy+=abs(image->data[(image->iw*(ky+ky2)+kx+kx2)<<2|kc]);//r3 opt 2.368454
		}
	}
#if 0
	int
		idx=image->iw*ky+kx,
		N=ky?image->data[(idx-image->iw)<<2|kc]:0,
		W=kx?image->data[(idx-1)<<2|kc]:0,
		NW=kx&&ky?image->data[(idx-image->iw-1)<<2|kc]:0,
		NE=kx+1<image->iw&&ky?image->data[(idx-image->iw+1)<<2|kc]:0,
		
		NN=ky-2>=0?image->data[(idx-image->iw*2)<<2|kc]:0,
		WW=kx-2>=0?image->data[(idx-2)<<2|kc]:0,
		NNWW=kx-2>=0&&ky-2>=0?image->data[(idx-image->iw*2-2)<<2|kc]:0,
		NNEE=kx+2<image->iw&&ky-2>=0?image->data[(idx-image->iw*2+2)<<2|kc]:0;

	//int g180=abs(W-WW), g135=abs(NW-NNWW), g90=abs(N-NN), g45=abs(NE-NNEE);
	

	int energy=abs(N)+abs(W)+abs(NW)+abs(NE)+((abs(NN)+abs(WW)+abs(NNWW)+abs(NNEE))>>2);//opt 2.404962
	
	//int energy=abs(N)+abs(W)+abs(NW)+abs(NE)+((abs(NN)+abs(WW)+abs(NNWW)+abs(NNEE))>>1);//opt 2.398070

	//int energy=abs(N)+abs(W)+abs(NW)+abs(NE);//opt 2.395828

	//int energy=abs(N)+abs(W)+((abs(NW)+abs(NE))>>1);//opt 2.390142
	
	//int energy=abs(N)+abs(W)+abs(NW)+abs(NE)+abs(NN)+abs(WW)+abs(NNWW)+abs(NNEE);//opt 2.388757

	//int energy=N*N+W*W+abs(NW)+abs(NE);//opt 2.365894

	//int energy=g180;//2.236316
	//if(energy>g135)energy=g135;
	//if(energy>g90)energy=g90;
	//if(energy>g45)energy=g45;

	//int energy=g180+g135+g90+g45;//opt 2.333547

	//int energy=N+W+NW+NE;//opt 2.332700

	//int energy;//opt 2.315842
	//int vmax=N, vmin=N;
	//if(vmax<W)vmax=W;
	//if(vmin>W)vmin=W;
	//if(vmax<NW)vmax=NW;
	//if(vmin>NW)vmin=NW;
	//if(vmax<NE)vmax=NE;
	//if(vmin>NE)vmin=NE;
	//energy=vmax-vmin;

	//int energy=N;//opt 2.262733
	//if(abs(energy)>abs(W))
	//	energy=W;
	//if(abs(energy)>abs(NW))
	//	energy=NW;
	//if(abs(energy)>abs(NE))
	//	energy=NE;
#endif


	int ctx=0;
	for(int k=0;k<_countof(ans_qlevels_u);++k)//TODO: binary search
		ctx+=energy>ans_qlevels_u[k];
	return ctx;
}
void calc_csize_ans_energy(Image const *image, size_t *csizes)
{
	if(!csizes)
		return;
	int maxdepth=calc_maxdepth(image, 0);
	int histsize=1<<maxdepth;//number of possibilities
	int *stats=(int*)malloc((histsize+1LL)*sizeof(int[_countof(calcsize_ans_qlevels)+1]));
	if(!stats)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(stats, 0, (histsize+1LL)*sizeof(int[_countof(calcsize_ans_qlevels)+1]));

	//int res=image->iw*image->ih;
	for(int kc=0;kc<3;++kc)			//step 1: fill histograms
	{
		int depth=image->depth[kc], nlevels=1<<depth;
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				int ctx=get_ctx(image, kc, kx, ky);
				int *hist=stats+(histsize+1)*ctx;
				int curr=image->data[idx<<2|kc]+(nlevels>>1);
				if((unsigned)curr>=(unsigned)histsize)
					LOG_ERROR("Token value OOB");
				++hist[curr];
			}
		}
	}
	//int npx=res*3;
	for(int kh=0;kh<_countof(calcsize_ans_qlevels)+1;++kh)		//step 3: accumulate & quantize CDFs
	{
		int *hist=stats+(histsize+1)*kh;
		int sum2=0;
		int npresent=0;
		for(int ks=0;ks<histsize;++ks)//3.1: count nonzero freq (present) symbols
		{
			sum2+=hist[ks];
			npresent=hist[ks]!=0;
		}
		if(sum2)
		{
			int sum=0;
			for(int ks=0, ks2=0;ks<histsize;++ks)//3.2: accumulate & quantize
			{
				int freq=hist[ks];
				hist[ks]=(int)((long long)sum*(0x10000-npresent)/sum2)+ks2;
				sum+=freq;
				ks2+=freq!=0;
			}
			hist[histsize]=0x10000;
		}
	}
	unsigned state=0x10000;
	int ctxhist[34]={0};
	for(int kc=0;kc<3;++kc)			//step 4: estimate compressed size
	{
		int depth=image->depth[kc], nlevels=1<<depth;
		size_t csize=0;
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				int ctx=get_ctx(image, kc, kx, ky);

				++ctxhist[ctx];

				int *hist=stats+(histsize+1)*ctx;
				int curr=image->data[idx<<2|kc]+(nlevels>>1);
				int cdf=hist[curr], freq=hist[curr+1]-cdf;

				if(!freq)
					LOG_ERROR("ZPS");
			
				if(state>=(unsigned)(freq<<16))//renorm
				{
					csize+=2;
					state>>=16;
				}
				state=state/freq<<16|(cdf+state%freq);//update
			}
		}
		csizes[kc]=csize;
	}
	csizes[3]=csizes[0]+csizes[1]+csizes[2]+4;
	free(stats);
}
static void print_hist_as_hexPDF(const int *hist, int nlevels)
{
	int sum=0, vmax=0;
	for(int k=0;k<nlevels;++k)
	{
		sum+=hist[k];
		if(vmax<hist[k])
			vmax=hist[k];
	}
	if(vmax)
	{
		for(int k=0;k<nlevels;++k)
			console_log("%X", (hist[k]*15+(vmax>>1))/vmax);
	}
	else
	{
		for(int k=0;k<nlevels;++k)//uniform
			console_log("0");
	}
	console_log(" %10d\n", sum);
}
void calc_csize_ans_energy_hybrid(Image const *image, size_t *csizes, const int *config)
{
	if(!csizes)
		return;
	int maxdepth=calc_maxdepth(image, 0), maxlevels=1<<maxdepth;
	int tbn[3];
	hybriduint_encode(-(maxlevels>>1), tbn, config);
	int histsize=tbn[0]+1;//number of possibilities
	int *stats=(int*)malloc((histsize+1LL)*sizeof(int[_countof(calcsize_ans_qlevels)+1]));
	Image *im2=0;
	image_copy_nodata(&im2, image);
	if(!stats||!im2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(stats, 0, (histsize+1LL)*sizeof(int[_countof(calcsize_ans_qlevels)+1]));

	int res=image->iw*image->ih;
	TempHybrid *data=(TempHybrid*)im2->data;
	for(int kc=0;kc<3;++kc)			//step 1: encode hybrid uint
	{
		for(int k=0;k<res;++k)
		{
			int pixel=image->data[k<<2|kc];
			hybriduint_encode(pixel, tbn, config);
			if((unsigned)tbn[0]>=256||(unsigned)tbn[2]>=256||(unsigned)tbn[1]>=0x10000)
				LOG_ERROR("Pixel value OOB");
			TempHybrid *curr=data+(k<<2|kc);
			curr->token=(unsigned char)tbn[0];
			curr->nbits=(unsigned char)tbn[2];
			curr->bypass=(unsigned short)tbn[1];
		}
	}
	for(int kc=0;kc<3;++kc)			//step 2: fill histograms
	{
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				int ctx=get_ctx(image, kc, kx, ky);
				int *hist=stats+(histsize+1)*ctx;
				TempHybrid *curr=data+(idx<<2|kc);
				if((unsigned)curr->token>=(unsigned)histsize)
					LOG_ERROR("Token value OOB");
				++hist[curr->token];
			}
		}
	}
	console_log("\n");//

	//int npx=res*3;
	for(int kh=0;kh<_countof(calcsize_ans_qlevels)+1;++kh)		//step 3: accumulate & quantize CDFs
	{
		int *hist=stats+(histsize+1)*kh;

		print_hist_as_hexPDF(hist, histsize);//

		int sum2=0;
		int npresent=0;
		for(int ks=0;ks<histsize;++ks)//3.1: count nonzero freq (present) symbols
		{
			sum2+=hist[ks];
			npresent+=hist[ks]!=0;
		}
		if(sum2)
		{
			int sum=0;
			for(int ks=0, ks2=0;ks<histsize;++ks)//3.2: accumulate & quantize
			{
				int freq=hist[ks];
				hist[ks]=(int)((long long)sum*(0x10000-npresent)/sum2)+ks2;
				sum+=freq;
				ks2+=freq!=0;
			}
			hist[histsize]=0x10000;
		}
	}
	unsigned state=0x10000;
	for(int kc=0;kc<3;++kc)			//step 4: estimate compressed size
	{
		size_t csize=0;
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				int ctx=get_ctx(image, kc, kx, ky);
				int *hist=stats+(histsize+1)*ctx;
				TempHybrid *curr=data+(idx<<2|kc);
				int cdf, freq;

				if(curr->nbits)
				{
					cdf=(curr->bypass<<16)/(1<<curr->nbits);
					freq=((curr->bypass+1)<<16)/(1<<curr->nbits)-cdf;

					if(state>=(unsigned)(freq<<16))//renorm
					{
						csize+=2;
						state>>=16;
					}
					state=state/freq<<16|(cdf+state%freq);//update
				}

				cdf=hist[curr->token];
				freq=hist[curr->token+1]-cdf;

				if(!freq)
					LOG_ERROR("ZPS");
			
				if(state>=(unsigned)(freq<<16))//renorm
				{
					csize+=2;
					state>>=16;
				}
				state=state/freq<<16|(cdf+state%freq);//update
			}
		}
		csizes[kc]=csize;
	}
	csizes[3]=csizes[0]+csizes[1]+csizes[2]+4;
	free(stats);
	free(im2);
}
static void test_predmask(Image const *image)
{
	//last value in kernel is lg(den)
	int hybrid_uint_config[]=
	{
		4, 2, 0,
	};
	static const int mask_left[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 1, 0,
	};
	static const int mask_top[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0,
	};
	static const int mask_av0[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 1, 1,
	};
	static const int mask_topright[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0,
	};
	static const int mask_topleft[]=
	{
		0, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 0,
	};
	static const int mask_leftleft[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		1, 0, 0,
	};
	static const int mask_av1[]=
	{
		0, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 1, 1,
	};
	static const int mask_av2[]=
	{
		0, 0, 0, 0, 0,
		0, 1, 1, 0, 0,
		0, 0, 1,
	};
	static const int mask_av3[]=
	{
		0, 0, 0, 0, 0,
		0, 0, 1, 1, 0,
		0, 0, 1,
	};
	static const int mask_av4[]=
	{
		0, 0, -2, 0, 0,
		0, 0,  6, 3, 1,
		1, 7,  4,
	};
	Image *residues[13]={0};
	for(int k=0;k<_countof(residues);++k)
		image_copy_nodata(residues+k, image);

	//original image == zero predictor
	pred_linear(image, residues[0], mask_left, mask_left[12], 1, 1);
	pred_linear(image, residues[1], mask_top, mask_top[12], 1, 1);
	pred_linear(image, residues[2], mask_av0, mask_av0[12], 1, 1);
	pred_select(image, residues[3], 1, 1);

	size_t bufsize=image_getbufsize(image);
	memcpy(residues[4], image, bufsize);
	pred_clampgrad(residues[4], 1, 1);

	memcpy(residues[5], image, bufsize);
	pred_jxl_apply(residues[5], 1, 1, jxlparams_i16);
	
	pred_linear(image, residues[6], mask_topright, mask_topright[12], 1, 1);
	pred_linear(image, residues[7], mask_topleft, mask_topleft[12], 1, 1);
	pred_linear(image, residues[8], mask_leftleft, mask_leftleft[12], 1, 1);
	pred_linear(image, residues[9], mask_av1, mask_av1[12], 1, 1);
	pred_linear(image, residues[10], mask_av2, mask_av2[12], 1, 1);
	pred_linear(image, residues[11], mask_av3, mask_av3[12], 1, 1);
	pred_linear(image, residues[12], mask_av4, mask_av4[12], 1, 1);
	
	//reused context histogram
#if 1
	ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
	unsigned char *masks=(unsigned char*)malloc(res*sizeof(int[3]));
	int maxdepth=calc_maxdepth(residues[0], 0), maxlevels=1<<maxdepth;
	int *stats=(int*)malloc(maxlevels*sizeof(int[34]));
	if(!masks||!stats)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(stats, 0, maxlevels*sizeof(int[34]));
	{
		int black=0xFF000000;
		memfill(masks, &black, res*sizeof(int[3]), sizeof(int));
	}
	for(int kp=0;kp<14;++kp)		//1: fill histograms
	{
		Image const *sample=kp?residues[kp-1]:image;
		for(int kc=0;kc<3;++kc)
		{
			int depth=sample->depth[kc], nlevels=1<<depth;
			for(int ky=0;ky<image->ih;++ky)
			{
				for(int kx=0;kx<image->iw;++kx)
				{
					int ctx=get_ctx(sample, kc, kx, ky);
					int *hist=stats+maxlevels*ctx;
					int sym=sample->data[(image->iw*ky+kx)<<2|kc]+(nlevels>>1);
					if(sym<0||sym>=nlevels)
						LOG_ERROR("Symbol overflow");
					++hist[sym];
				}
			}
		}
	}
	int hweight[34]={0};
	for(int ctx=0;ctx<_countof(hweight);++ctx)	//2: get histogram sums
	{
		int *hist=stats+maxlevels*ctx;
		for(int ks=0;ks<maxlevels;++ks)
			hweight[ctx]+=hist[ks];
	}
	double csizes[14][3]={0}, optcsize[3]={0};
	
	static const int predcolors[]=
	{
		0xFFFFFFFF,//	zero		undefined	white
		0xFF0000FF,//	left		180		red
		0xFFFF0000,//	top		90		blue
		0xFF00D000,//	average0	135		green
		0xFF009000,//	select		135		green
		0xFF005000,//	grad		135		green
		0xFF000000,//	weighted	undefined	black
		0xFFFF00FF,//	topright	45		pink
		0xFF00FF00,//	topleft		135		green
		0xFF000080,//	leftleft	180		red
		0xFF008080,//	average1	157.5		mustard
		0xFF808000,//	average2	112.5		marine
		0xFF800080,//	average3	67.5		violet
		0xFF808080,//	average4	undefined	grey
	};
	for(int kc=0;kc<3;++kc)				//3: calc csizes
	{
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx)
			{
				double sizes[14];
				for(int kp=0;kp<14;++kp)
				{
					Image const *sample=kp?residues[kp-1]:image;
					int nlevels=1<<sample->depth[kc];
					int ctx=get_ctx(sample, kc, kx, ky);
					int *hist=stats+maxlevels*ctx;
					int sym=sample->data[(image->iw*ky+kx)<<2|kc]+(nlevels>>1);
					if((unsigned)sym>=(unsigned)nlevels)
						LOG_ERROR("Symbol OOB");
					int freq=hist[sym];
					sizes[kp]=-log2((double)freq/hweight[ctx])/image->src_depth[kc];
				}
				int pbest=0;
				double minsize=sizes[0];
				for(int kp=1;kp<14;++kp)//get min size
				{
					if(minsize>sizes[kp])
						minsize=sizes[kp], pbest=kp;
				}
				((int*)masks)[res*kc+(image->iw*ky+kx)]=predcolors[pbest];
				//for(int kp=0;kp<14;++kp)
				//	masks[(res<<2)*kp+((image->iw*ky+kx)<<2|kc)]=minsize==sizes[kp]?0xFF:0;
				for(int kp=0;kp<14;++kp)
					csizes[kp][kc]+=sizes[kp];
				optcsize[kc]+=minsize;
			}
		}
	}
#endif

	//dedicated naive histogram
#if 0
	ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
	unsigned char *masks=(unsigned char*)malloc(res*sizeof(int[14]));
	int maxdepth=calc_maxdepth(residues[0], 0), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc(maxlevels*sizeof(int[14*3]));
	if(!masks||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, maxlevels*sizeof(int[14*3]));
	{
		int black=0xFF000000;
		memfill(masks, &black, res*sizeof(int[14]), sizeof(int));
	}
	for(int kp=0;kp<14;++kp)//calc histograms
	{
		Image const *sample=kp?residues[kp-1]:image;
		for(int kc=0;kc<3;++kc)
		{
			int depth=sample->depth[kc], nlevels=1<<depth;
			int *hk=hist+maxlevels*(3*kp+kc);
			for(ptrdiff_t k=0;k<res;++k)
			{
				int sym=sample->data[k<<2|kc]+(nlevels>>1);
				if(sym<0||sym>=nlevels)
					LOG_ERROR("Symbol overflow");
				++hk[sym];
			}
		}
	}
	double csizes[14][3]={0}, optcsize[3]={0};
	for(int kc=0;kc<3;++kc)
	{
		for(ptrdiff_t k=0;k<res;++k)//calc csizes
		{
			double bytesizes[14];
			for(int kp=0;kp<14;++kp)
			{
				Image const *sample=kp?residues[kp-1]:image;
				const int *hk=hist+maxlevels*(3*kp+kc);
				int nlevels=1<<sample->depth[kc];
				int sym=sample->data[k<<2|kc]+(nlevels>>1);
				int freq=hk[sym];
				if(!freq)
					LOG_ERROR("ZPS");
				bytesizes[kp]=-log2((double)freq/res)/image->src_depth[kc];
			}
			double minbytesize=bytesizes[0];
			for(int kp=1;kp<14;++kp)
			{
				if(minbytesize>bytesizes[kp])
					minbytesize=bytesizes[kp];
			}
			for(int kp=0;kp<14;++kp)
				masks[(res<<2)*kp+(k<<2|kc)]=minbytesize==bytesizes[kp]?0xFF:0;
			for(int kp=0;kp<14;++kp)
				csizes[kp][kc]+=bytesizes[kp];
			optcsize[kc]+=minbytesize;
		}
	}
#endif

	static const char *predname[]=
	{
		"zero",		//00
		"left",		//01
		"top",		//02
		"av0",		//03
		"sel",		//04
		"grad",		//05
		"wp",		//06
		"topright",	//07
		"topleft",	//08
		"leftleft",	//09
		"av1",		//10
		"av2",		//11
		"av3",		//12
		"av4",		//13
	};
	console_start();
	double usize=image_getBMPsize(image);
	console_log("%-10s %16lf\n\n", "usize", usize);
	for(int kp=0;kp<14;++kp)//print csizes
	{
		Image const *sample=kp?residues[kp-1]:image;
		size_t csizes_ans[4]={0};
		calc_csize_ans_energy_hybrid(sample, csizes_ans, hybrid_uint_config);
		//calc_csize_ans_energy(sample, csizes_ans);
		//calc_csize_ans_separate(sample, csizes_ans);//sanity check

		double csize=csizes[kp][0]+csizes[kp][1]+csizes[kp][2];
		console_log("%-10s TYUV %12.2lf %12.2lf %12.2lf %12.2lf  invCR %lf%%  TYUV %10lld %10lld %10lld %10lld  invCR %lf\n",
			predname[kp], csize, csizes[kp][0], csizes[kp][1], csizes[kp][2], 100.*csize/usize,
			csizes_ans[3], csizes_ans[0], csizes_ans[1], csizes_ans[2], 100.*(double)csizes_ans[3]/usize
		);
	}
	double csize2=optcsize[0]+optcsize[1]+optcsize[2];
	console_log("\n%-10s TRGB %16lf %16lf %16lf %16lf  invCR %lf%%\n", "optsize", csize2, optcsize[0], optcsize[1], optcsize[2], 100.*csize2/usize);

	console_log("\nAbout to save predictor masks\n");
	console_pause();
	ArrayHandle path=dialog_open_folder();
	if(path)
	{
		const char *title;
		int titlelen;
		{//get title without extension from filename
			int titleend=(int)acme_strrchr((char*)fn->data, fn->count, '.');
			int titlestart=titleend;
			for(;titlestart>=0&&fn->data[titlestart]!='/'&&fn->data[titlestart]!='\\';--titlestart);
			++titlestart;

			title=(char*)fn->data+titlestart;
			titlelen=titleend-titlestart;
		}

#if 1
		static const char chnames[]="YUV";
		for(int kc=0;kc<3;++kc)
		{
			snprintf(g_buf, G_BUF_SIZE, "%s%.*s-%d%c.PNG", (char*)path->data, titlelen, title, kc, chnames[kc]);
			unsigned char *mask=masks+(res<<2)*kc;
			lodepng_encode_file(g_buf, mask, image->iw, image->ih, LCT_RGBA, 8);
		}
#endif

#if 0
		for(int kp=0;kp<14;++kp)//save snapshots
		{
			snprintf(g_buf, G_BUF_SIZE, "%s%.*s-%02d-%s.PNG", (char*)path->data, titlelen, title, kp, predname[kp]);
			
#if 0
			Image const *sample=kp?residues[kp-1]:image;		//saves residues
			image_save_uint8(g_buf, sample, 1);
#else
			unsigned char *mk=masks+(res<<2)*kp;			//saves pred masks
			lodepng_encode_file(g_buf, mk, image->iw, image->ih, LCT_RGBA, 8);
#endif
		}
#endif
		array_free(&path);
	}

	for(int k=0;k<13;++k)
		free(residues[k]);
	free(masks);
	free(stats);

	console_log("\nDone.\n");
	console_pause();
	console_end();
}

static void calc_csize_stateful(Image const *image, int *hist_full, double *entropy)
{
	if(ec_method==ECTX_HIST)
	{
		int allocated=0;
		if(!hist_full)
		{
			int maxdepth=calc_maxdepth(im1, 0);
			hist_full=(int*)malloc(sizeof(int)<<maxdepth);
			allocated=1;
		}
		for(int kc=0;kc<4;++kc)
		{
			if(image->depth[kc])
			{
				calc_histogram(image->data, image->iw, image->ih, kc, 0, image->iw, 0, image->ih, image->depth[kc], hist_full, 0);
				entropy[kc]=calc_entropy(hist_full, 1<<image->depth[kc], image->iw*image->ih);
			}
			else
				entropy[kc]=0;
		}
		if(allocated)
			free(hist_full);
		//channel_entropy(image, iw*ih, 3, 4, ch_cr, usage);
	}
	else if(ec_method==ECTX_ABAC)
		calc_csize_abac(image, entropy);
	else
		calc_csize_ec(image, ec_method, ec_adaptive?ec_adaptive_threshold:0, ec_expbits, ec_msb, ec_lsb, entropy);
}

typedef struct ThreadCtxStruct
{
	Image *image;
	double usize, csize[4];
	ptrdiff_t idx;
} ThreadCtx;
static unsigned __stdcall sample_thread(void *param)
{
	ThreadCtx *ctx=(ThreadCtx*)param;
	ctx->usize=image_getBMPsize(ctx->image);
	apply_selected_transforms(ctx->image, 0);
	int maxdepth=calc_maxdepth(ctx->image, 0);
	int nlevels=1<<maxdepth;
	int *hist=(int*)malloc(nlevels*sizeof(int));
	double entropy[4]={0};
	calc_csize_stateful(ctx->image, 0, entropy);
	for(int kc=0;kc<4;++kc)
	{
		int depth=ctx->image->src_depth[kc];
		double invCR=depth?entropy[kc]/depth:0;
		ctx->csize[kc]=invCR*ctx->image->iw*ctx->image->ih*ctx->image->src_depth[kc]/8;
	}
	//for(int kc=0;kc<3;++kc)
	//{
	//	calc_histogram(ctx->image->data, ctx->image->iw, ctx->image->ih, kc, 0, ctx->image->iw, 0, ctx->image->ih, ctx->image->depth[kc], hist, 0);
	//	double entropy=calc_entropy(hist, 1<<ctx->image->depth[kc], ctx->image->iw*ctx->image->ih);
	//	double invCR=entropy/ctx->image->src_depth[kc];
	//	ctx->csize[kc]=invCR*ctx->image->iw*ctx->image->ih*ctx->image->src_depth[kc]/8;
	//}
	free(hist);
	free(ctx->image);
	//_endthreadex(0);
	return 0;
}
static void batch_test(void)
{
	loud_transforms=0;
	ArrayHandle path=dialog_open_folder();
	if(!path)
		return;
	const char *ext[]=
	{
		"PNG",
		"JPG",
		"JPEG",
		"PPM",
		"PGM",
	};
	ArrayHandle filenames=get_filenames((char*)path->data, ext, _countof(ext), 1);
	if(!filenames)
	{
		array_free(&path);
		return;
	}

	DisableProcessWindowsGhosting();
	console_start();
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H:%M:%S");
	console_log("Batch Test  %s  %s\n", g_buf, (char*)path->data);
	array_free(&path);
	console_log("Enter number of threads: ");
	int nthreads=console_scan_int();
	double t=time_sec();
	double total_usize=0, total_csize[4]={0};
	int maxlen=0;
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ArrayHandle q;
	ARRAY_ALLOC(ThreadCtx, q, 0, 0, nthreads, 0);
	for(int k=0;k<(int)filenames->count;++k)
	{
		//multi-threaded
#if 1
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		Image *image=image_load((char*)fn2[0]->data, (int)fn2[0]->count);
		if(!image)
			continue;
		ThreadCtx ctx=
		{
			image,
			0, {0},
			k,
		};
		ARRAY_APPEND(q, &ctx, 1, 1, 0);
		if((int)q->count>=nthreads||k+1>=(int)filenames->count)
		{
			HANDLE *handles=(HANDLE*)malloc(q->count*sizeof(HANDLE));
			if(!handles)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			memset(handles, 0, q->count*sizeof(HANDLE));
			for(int k2=0;k2<(int)q->count;++k2)
			{
				ThreadCtx *ptr=(ThreadCtx*)array_at(&q, k2);
				handles[k2]=(void*)_beginthreadex(0, 0, sample_thread, ptr, 0, 0);
				if(!handles[k2])
				{
					LOG_ERROR("Thread alloc error");
					return;
				}
			}
			WaitForMultipleObjects((int)q->count, handles, TRUE, INFINITE);
			for(int k2=0;k2<(int)q->count;++k2)
			{
				ThreadCtx *ptr=(ThreadCtx*)array_at(&q, k2);
				fn2=(ArrayHandle*)array_at(&filenames, ptr->idx);
				double csize=ptr->csize[0]+ptr->csize[1]+ptr->csize[2]+ptr->csize[3];
				console_log(
					"%5d/%5d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  invCR %8.4lf%%\n",
					(int)(k+1-(int)q->count+k2+1), (int)filenames->count, (char*)fn2[0]->data, (int)(maxlen-fn2[0]->count+1), "",
					ptr->usize, csize, ptr->csize[0], ptr->csize[1], ptr->csize[2],
					100.*csize/ptr->usize
				);
				total_usize+=ptr->usize;
				total_csize[0]+=ptr->csize[0];
				total_csize[1]+=ptr->csize[1];
				total_csize[2]+=ptr->csize[2];
				total_csize[3]+=ptr->csize[3];
			}
			array_clear(&q);
		}
#endif

		//single-threaded
#if 0
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		Image *image=image_load((char*)fn2[0]->data);
		if(!image)
			continue;
		double usize=image_getBMPsize(image), csize[3]={0};
		apply_selected_transforms(image, 0);
		int maxdepth=calc_maxdepth(image, 0);
		int nlevels=1<<maxdepth;
		int *hist=(int*)malloc(nlevels*sizeof(int));
		for(int kc=0;kc<3;++kc)
		{
			calc_histogram(image->data, image->iw, image->ih, kc, 0, image->iw, 0, image->ih, image->depth[kc], hist, 0);
			double entropy=calc_entropy(hist, 1<<image->depth[kc], image->iw*image->ih);
			double invCR=entropy/image->src_depth[kc];
			csize[kc]=invCR*image->iw*image->ih*image->src_depth[kc]/8;
		}
		free(hist);
		console_log(
			"%3d/%3d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf\n",
			k+1, (int)filenames->count, (char*)fn2[0]->data, maxlen-fn2[0]->count+1, "",
			usize, csize[0]+csize[1]+csize[2]+csize[3],
			csize[0], csize[1], csize[2]
		);
		total_usize+=usize;
		total_csize[0]+=csize[0];
		total_csize[1]+=csize[1];
		total_csize[2]+=csize[2];
		total_csize[3]+=csize[3];
		free(image);
#endif
	}
	double ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
	double CR=total_usize/ctotal;
	t=time_sec()-t;
	console_log(
		"Total UTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  invCR %8.4lf%%\n",
		total_usize, ctotal, total_csize[0], total_csize[1], total_csize[2], 100./CR
	);
	timedelta2str(g_buf, G_BUF_SIZE, t);
	console_log("Elapsed %s\n", g_buf);
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d-%H:%M:%S");
	console_log("\nDone.  %s\n", g_buf);
	console_pause();
	console_end();
	loud_transforms=1;
}

static int customtransforms_getflag(unsigned char tid)
{
	return
		tid==CT_FWD_ADAPTIVE||
		tid==CT_INV_ADAPTIVE||
		tid==CT_FWD_CUSTOM||
		tid==CT_INV_CUSTOM||
		tid==ST_FWD_CUSTOM||
		tid==ST_INV_CUSTOM;
	//	tid==ST_FWD_WP||
	//	tid==ST_INV_WP;
	//	tid==ST_FWD_EXPDWT||
	//	tid==ST_INV_EXPDWT||
	//	tid==ST_FWD_GRAD2||
	//	tid==ST_INV_GRAD2||
	//	tid==ST_FWD_CUSTOM_DWT||
	//	tid==ST_INV_CUSTOM_DWT;
}
static void transforms_update(void)
{
	if(transforms)
	{
		memset(transforms_mask, 0, T_COUNT);
		transforms_customenabled=0;
		for(int k=0;k<(int)transforms->count;++k)//update trackers
		{
			unsigned char tid2=transforms->data[k];
			if(tid2<T_COUNT)
				transforms_mask[tid2]|=1;
			transforms_customenabled|=customtransforms_getflag(tid2);
		}
	}
}
static void transforms_removebyid(unsigned tid)
{
	if(transforms)
	{
		for(int k=(int)transforms->count-1;k>=0;--k)
		{
			if(tid==transforms->data[k])
			{
				array_erase(&transforms, k, 1);
				break;//remove only one
			}
		}
		transforms_update();
	}
}
static void transforms_removeall(void)
{
	array_free(&transforms);
	transforms_customenabled=0;
	memset(transforms_mask, 0, T_COUNT);
}
static void transforms_append(unsigned tid)
{
	if(tid<T_COUNT&&tid!=CST_FWD_SEPARATOR&&tid!=CST_INV_SEPARATOR)
	{
		if(!transforms)
		{
			ARRAY_ALLOC(char, transforms, &tid, 1, 0, 0);
			transforms_mask[tid]|=1;
			transforms_customenabled|=customtransforms_getflag((unsigned char)tid);
		}
		else
		{
			int idx=-1;//first idx of a transform of this type
			if(GET_KEY_STATE(KEY_CTRL))//replace all transforms of this type
			{
				for(int k=0;k<(int)transforms->count;)
				{
					if((tid<CST_FWD_SEPARATOR)==(transforms->data[k]<CST_FWD_SEPARATOR))
					{
						if(idx==-1)
							idx=k;
						array_erase(&transforms, k, 1);
						transforms_mask[transforms->data[k]]=0;
					}
					else
						++k;
				}
			}
			if(idx==-1)
				ARRAY_APPEND(transforms, &tid, 1, 1, 0);
			else
				array_insert(&transforms, idx, &tid, 1, 1, 0);
			transforms_update();
		}
	}
}
static void transforms_printname(float x, float y, unsigned tid, int place, long long highlight)
{
	const char *a=0;
	switch(tid)
	{
	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
	case CT_FWD_YCbCr_R_v1:		a="C  Fwd YCbCr-R v1";		break;
	case CT_INV_YCbCr_R_v1:		a="C  Inv YCbCr-R v1";		break;
	case CT_FWD_YCbCr_R_v2:		a="C  Fwd YCbCr-R v2";		break;
	case CT_INV_YCbCr_R_v2:		a="C  Inv YCbCr-R v2";		break;
	case CT_FWD_YCbCr_R_v3:		a="C  Fwd YCbCr-R v3";		break;
	case CT_INV_YCbCr_R_v3:		a="C  Inv YCbCr-R v3";		break;
	case CT_FWD_YCbCr_R_v4:		a="C  Fwd YCbCr-R v4";		break;
	case CT_INV_YCbCr_R_v4:		a="C  Inv YCbCr-R v4";		break;
	case CT_FWD_YCbCr_R_v5:		a="C  Fwd YCbCr-R v5";		break;
	case CT_INV_YCbCr_R_v5:		a="C  Inv YCbCr-R v5";		break;
	case CT_FWD_YCbCr_R_v6:		a="C  Fwd YCbCr-R v6";		break;
	case CT_INV_YCbCr_R_v6:		a="C  Inv YCbCr-R v6";		break;
	case CT_FWD_YCbCr_R_v7:		a="C  Fwd YCbCr-R v7";		break;
	case CT_INV_YCbCr_R_v7:		a="C  Inv YCbCr-R v7";		break;
//	case CT_FWD_CrCgCb:		a="C  Fwd CrCgCb";		break;
//	case CT_INV_CrCgCb:		a="C  Inv CrCgCb";		break;
	case CT_FWD_Pei09:		a="C  Fwd Pei09";		break;
	case CT_INV_Pei09:		a="C  Inv Pei09";		break;
	case CT_FWD_YCoCg_R:		a="C  Fwd YCoCg-R";		break;
	case CT_INV_YCoCg_R:		a="C  Inv YCoCg-R";		break;
	case CT_FWD_JPEG2000:		a="C  Fwd JPEG2000 RCT";	break;
	case CT_INV_JPEG2000:		a="C  Inv JPEG2000 RCT";	break;
	case CT_FWD_SUBGREEN:		a="C  Fwd SubGreen";		break;
	case CT_INV_SUBGREEN:		a="C  Inv SubGreen";		break;
//	case CT_FWD_YRGB_v1:		a="C  Fwd YRGB v1";		break;
//	case CT_INV_YRGB_v1:		a="C  Inv YRGB v1";		break;
//	case CT_FWD_YRGB_v2:		a="C  Fwd YRGB v2";		break;
//	case CT_INV_YRGB_v2:		a="C  Inv YRGB v2";		break;
//	case CT_FWD_CMYK:		a="C  Fwd CMYK";		break;
//	case CT_INV_CMYK_DUMMY:		a="";				break;//no inverse
	case CT_FWD_YCbCr:		a="C  Fwd YCbCr";		break;
	case CT_INV_YCbCr:		a="C  Inv YCbCr";		break;
	case CT_FWD_XYB:		a="C  Fwd XYB";			break;
	case CT_INV_XYB:		a="C  Inv XYB";			break;
//	case CT_FWD_MATRIX:		a="C  Fwd Matrix";		break;
//	case CT_INV_MATRIX:		a="C  Inv Matrix";		break;
//	case CT_FWD_XGZ:		a="C  Fwd XGZ";			break;
//	case CT_INV_XGZ:		a="C  Inv XGZ";			break;
//	case CT_FWD_XYZ:		a="C  Fwd XYZ";			break;
//	case CT_INV_XYZ:		a="C  Inv XYZ";			break;
//	case CT_FWD_EXP:		a="C  Fwd Experimental";	break;
//	case CT_INV_EXP:		a="C  Inv Experimental";	break;
//	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
//	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
	case CT_FWD_CUSTOM:		a="C  Fwd CUSTOM";		break;
	case CT_INV_CUSTOM:		a="C  Inv CUSTOM";		break;
//	case CT_FWD_QUAD:		a="C  Fwd Quad";		break;
//	case CT_INV_QUAD:		a="C  Inv Quad";		break;

	case CST_FWD_SEPARATOR:		a="";				break;
	case CST_INV_SEPARATOR:		a="";				break;

//	case ST_PREPROC_GRAD:		a=" S Preproc Grad";		break;
//	case ST_PREPROC_X:		a=" S Preproc X";		break;
//	case ST_PREPROC_X2:		a=" S Preproc X2";		break;
		
	case ST_FWD_PACKSIGN:		a=" S Fwd PackSign";		break;
	case ST_INV_PACKSIGN:		a=" S Fwd PackSign";		break;
	case ST_FWD_P3:			a=" S Fwd P3";			break;
	case ST_INV_P3:			a=" S Inv P3";			break;
	case ST_FWD_G2:			a=" S Fwd G2";			break;
	case ST_INV_G2:			a=" S Inv G2";			break;
	case ST_FWD_T47:		a=" S Fwd T47";			break;
	case ST_INV_T47:		a=" S Inv T47";			break;
	case ST_FWD_OLS:		a=" S Fwd OLS";			break;
	case ST_INV_OLS:		a=" S Inv OLS";			break;
	case ST_FWD_OLS2:		a=" S Fwd OLS-2";		break;
	case ST_INV_OLS2:		a=" S Inv OLS-2";		break;
	case ST_FWD_OLS3:		a=" S Fwd OLS-3";		break;
	case ST_INV_OLS3:		a=" S Inv OLS-3";		break;
	case ST_FWD_OLS4:		a=" S Fwd OLS-4";		break;
	case ST_INV_OLS4:		a=" S Inv OLS-4";		break;
	case ST_FWD_OLS5:		a=" S Fwd OLS-5";		break;
	case ST_INV_OLS5:		a=" S Inv OLS-5";		break;
	case ST_FWD_PU:			a="CS Fwd PU";			break;
	case ST_INV_PU:			a="CS Inv PU";			break;
	case ST_FWD_CG3D:		a="CS Fwd CG3D";		break;
	case ST_INV_CG3D:		a="CS Inv CG3D";		break;
	case ST_FWD_CLAMPGRAD:		a=" S Fwd ClampGrad";		break;
	case ST_INV_CLAMPGRAD:		a=" S Inv ClampGrad";		break;
	case ST_FWD_AV2:		a=" S Fwd (N+W)>>1";		break;
	case ST_INV_AV2:		a=" S Inv (N+W)>>1";		break;
	case ST_FWD_MEDIAN:		a=" S Fwd Median";		break;
	case ST_INV_MEDIAN:		a=" S Inv Median";		break;
//	case ST_FWD_ECOEFF:		a=" S Fwd E-Coeff";		break;
//	case ST_INV_ECOEFF:		a=" S Inv E-Coeff";		break;
//	case ST_FWD_AVERAGE:		a=" S Fwd Average";		break;
//	case ST_INV_AVERAGE:		a=" S Inv Average";		break;
	case ST_FWD_MULTISTAGE:		a=" S Fwd Multistage";		break;
	case ST_INV_MULTISTAGE:		a=" S Inv Multistage";		break;
//	case ST_FWD_ZIPPER:		a=" S Fwd Zipper";		break;
//	case ST_INV_ZIPPER:		a=" S Inv Zipper";		break;
//	case ST_FWD_DIR:		a=" S Fwd Dir";			break;
//	case ST_INV_DIR:		a=" S Inv Dir";			break;
	case ST_FWD_CUSTOM3:		a="CS Fwd CUSTOM3";		break;
	case ST_INV_CUSTOM3:		a="CS Inv CUSTOM3";		break;
	case ST_FWD_CALIC:		a=" S Fwd CALIC";		break;
	case ST_INV_CALIC:		a=" S Inv CALIC";		break;
//	case ST_FWD_NBLIC:		a=" S Fwd NBLIC";		break;
//	case ST_INV_NBLIC:		a=" S Inv NBLIC";		break;
	case ST_FWD_WP:			a=" S Fwd JXL WP";		break;
	case ST_INV_WP:			a=" S Inv JXL WP";		break;
	case ST_FWD_WPU:		a="CS Fwd WPU";			break;
	case ST_INV_WPU:		a="CS Inv WPU";			break;
	case ST_FWD_DEFERRED:		a=" S Fwd DEFERRED";		break;
	case ST_INV_DEFERRED:		a=" S Inv DEFERRED";		break;
	case ST_FWD_MM:			a=" S Fwd MM";			break;
	case ST_INV_MM:			a=" S Inv MM";			break;
	case ST_FWD_CUSTOM:		a=" S Fwd CUSTOM";		break;
	case ST_INV_CUSTOM:		a=" S Inv CUSTOM";		break;
#if 0
//	case ST_FWD_CUSTOM2:		a=" S Fwd CUSTOM2";		break;
//	case ST_INV_CUSTOM2:		a=" S Inv CUSTOM2";		break;
//	case ST_FWD_CUSTOM4:		a=" S Fwd CUSTOM4";		break;
//	case ST_INV_CUSTOM4:		a=" S Inv CUSTOM4";		break;
//	case ST_FWD_KALMAN:		a=" S Fwd Kalman";		break;
//	case ST_INV_KALMAN:		a=" S Inv Kalman";		break;
//	case ST_FWD_LOGIC:		a=" S Fwd Logic";		break;
//	case ST_INV_LOGIC:		a=" S Inv Logic";		break;
//	case ST_FWD_LEARNED:		a=" S Fwd Learned";		break;
//	case ST_INV_LEARNED:		a=" S Inv Learned";		break;
#ifdef ALLOW_OPENCL
//	case ST_FWD_LEARNED_GPU:	a=" S Fwd Learned GPU";		break;
//	case ST_INV_LEARNED_GPU:	a=" S Inv Learned GPU";		break;
#endif
//	case ST_FWD_CFL:		a=" S Fwd CfL";			break;
//	case ST_INV_CFL:		a=" S Inv CfL";			break;
//	case ST_FWD_JOINT:		a=" S Fwd Joint";		break;
//	case ST_INV_JOINT:		a=" S Inv Joint";		break;
//	case ST_FWD_HYBRID3:		a=" S Fwd Hybrid3";		break;
//	case ST_INV_HYBRID3:		a=" S Inv Hybrid3";		break;

//	case ST_FWD_DIFF2D:		a=" S Fwd 2D derivative";	break;
//	case ST_INV_DIFF2D:		a=" S Inv 2D derivative";	break;
//	case ST_FWD_HPF:		a=" S Fwd HPF";			break;
//	case ST_INV_HPF:		a=" S Inv HPF";			break;
//	case ST_FWD_GRAD2:		a=" S Fwd Grad2";		break;
//	case ST_INV_GRAD2:		a=" S Inv Grad2";		break;
//	case ST_FWD_ADAPTIVE:		a=" S Fwd Adaptive";		break;
//	case ST_INV_ADAPTIVE:		a=" S Inv Adaptive";		break;
//	case ST_FWD_JMJ:		a=" S Fwd JMJ";			break;
//	case ST_INV_JMJ:		a=" S Inv JMJ";			break;
//	case ST_FWD_SORTNB:		a=" S Fwd Sort Nb";		break;
//	case ST_INV_SORTNB:		a=" S Inv Sort Nb";		break;
//	case ST_FWD_MEDIAN:		a=" S Fwd Median";		break;
//	case ST_INV_MEDIAN:		a=" S Inv Median";		break;
//	case ST_FWD_DCT3PRED:		a=" S Fwd DCT3 Predictor";	break;
//	case ST_INV_DCT3PRED:		a=" S Inv DCT3 Predictor";	break;
//	case ST_FWD_PATHPRED:		a=" S Fwd Path Predictor";	break;
//	case ST_INV_PATHPRED:		a=" S Inv Path Predictor";	break;
	case ST_FWD_GRADPRED:		a=" S Fwd Gradient";		break;
	case ST_INV_GRADPRED:		a=" S Inv Gradient";		break;
	case ST_FWD_GRAD2:		a=" S Fwd Grad2";		break;
	case ST_INV_GRAD2:		a=" S Inv Grad2";		break;
	case ST_FWD_CTX:		a=" S Fwd Ctx";			break;
	case ST_INV_CTX:		a=" S Inv Ctx";			break;
//	case ST_FWD_C03:		a=" S Fwd C03";			break;
//	case ST_INV_C03:		a=" S Inv C03";			break;
//	case ST_FWD_C10:		a=" S Fwd C10";			break;
//	case ST_INV_C10:		a=" S Inv C10";			break;
	case ST_FWD_C20:		a=" S Fwd C20";			break;
	case ST_INV_C20:		a=" S Inv C20";			break;
//	case ST_FWD_WU97:		a=" S Fwd Wu 97";		break;
//	case ST_INV_WU97:		a=" S Inv Wu 97";		break;
//	case ST_FWD_BITWISE:		a=" S Fwd Bitwise";		break;
//	case ST_INV_BITWISE:		a=" S Inv Bitwise";		break;
	case ST_FWD_DCT4:		a=" S Fwd DCT4";		break;
	case ST_INV_DCT4:		a=" S Inv DCT4";		break;
//	case ST_FWD_DCT8:		a=" S Fwd DCT8";		break;
//	case ST_INV_DCT8:		a=" S Inv DCT8";		break;
	case ST_FWD_SHUFFLE:		a=" S Fwd Shuffle";		break;
	case ST_INV_SHUFFLE:		a=" S Inv Shuffle";		break;
//	case ST_FWD_SPLIT:		a=" S Fwd Split";		break;
//	case ST_INV_SPLIT:		a=" S Inv Split";		break;
	case ST_FWD_LAZY:		a=" S Fwd Lazy DWT";		break;
	case ST_INV_LAZY:		a=" S Inv Lazy DWT";		break;
#endif
	case ST_FWD_HAAR:		a=" S Fwd Haar";		break;
	case ST_INV_HAAR:		a=" S Inv Haar";		break;
	case ST_FWD_SQUEEZE:		a=" S Fwd Squeeze";		break;
	case ST_INV_SQUEEZE:		a=" S Inv Squeeze";		break;
	case ST_FWD_LEGALL53:		a=" S Fwd LeGall 5/3";		break;
	case ST_INV_LEGALL53:		a=" S Inv LeGall 5/3";		break;
	case ST_FWD_CDF97:		a=" S Fwd CDF 9/7";		break;
	case ST_INV_CDF97:		a=" S Inv CDF 9/7";		break;
//	case ST_FWD_GRAD_DWT:		a=" S Fwd Gradient DWT";	break;
//	case ST_INV_GRAD_DWT:		a=" S Inv Gradient DWT";	break;
//	case ST_FWD_DEC_DWT:		a=" S Fwd Dec. DWT";		break;
//	case ST_INV_DEC_DWT:		a=" S Inv Dec. DWT";		break;
//	case ST_FWD_EXPDWT:		a=" S Fwd Exp DWT";		break;
//	case ST_INV_EXPDWT:		a=" S Inv Exp DWT";		break;
//	case ST_FWD_CUSTOM_DWT:		a=" S Fwd CUSTOM DWT";		break;
//	case ST_INV_CUSTOM_DWT:		a=" S Inv CUSTOM DWT";		break;
	case ST_FWD_DCT4:		a=" S Fwd DCT4";		break;
	case ST_INV_DCT4:		a=" S Inv DCT4";		break;
	case ST_FWD_DCT8:		a=" S Fwd DCT8";		break;
	case ST_INV_DCT8:		a=" S Inv DCT8";		break;
	default:			a="ERROR";			break;
	}
	long long c0=0;
	if(highlight)
		c0=set_text_colors(highlight);
	if(place<0)
		GUIPrint(0, x, y, 1, "%s", a);
	else
		GUIPrint(0, x, y, 1, "%d: %s", place, a);
	if(highlight)
		set_text_colors(c0);
}

static int send_image_separate_subpixels(Image const *image, unsigned *txid_r, unsigned *txid_g, unsigned *txid_b, unsigned *txid_a)
{
	ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
	unsigned char *temp_r=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_g=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_b=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_a=(unsigned char*)malloc(res*sizeof(char[4]));
	if(!temp_r||!temp_g||!temp_b||!temp_a)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	int shift[]=
	{
		MAXVAR(0, image->depth[0]-8),
		MAXVAR(0, image->depth[1]-8),
		MAXVAR(0, image->depth[2]-8),
		MAXVAR(0, image->depth[3]-8),
	};
	for(ptrdiff_t k=0;k<res;++k)
	{
		int
			r=((image->data[k<<2|0]>>shift[0])+128)&0xFF,
			g=((image->data[k<<2|1]>>shift[1])+128)&0xFF,
			b=((image->data[k<<2|2]>>shift[2])+128)&0xFF,
			a=((image->data[k<<2|3]>>shift[3])+128)&0xFF;
		if(separate_grayscale)
		{
			((int*)temp_r)[k]=0xFF000000|r<<16|r<<8|r;
			((int*)temp_g)[k]=0xFF000000|g<<16|g<<8|g;
			((int*)temp_b)[k]=0xFF000000|b<<16|b<<8|b;
		}
		else
		{
			((int*)temp_r)[k]=0xFF000000|r;
			((int*)temp_g)[k]=0xFF000000|g<<8;
			((int*)temp_b)[k]=0xFF000000|b<<16;
		}
		((int*)temp_a)[k]=0xFF000000|a<<16|a<<8|a;
	}
	if(!*txid_r)
	{
		unsigned txids[4];
		glGenTextures(4, txids);
		*txid_r=txids[0];
		*txid_g=txids[1];
		*txid_b=txids[2];
		*txid_a=txids[3];
	}
	send_texture_pot(*txid_r, (int*)temp_r, image->iw, image->ih, 0);
	send_texture_pot(*txid_g, (int*)temp_g, image->iw, image->ih, 0);
	send_texture_pot(*txid_b, (int*)temp_b, image->iw, image->ih, 0);
	send_texture_pot(*txid_a, (int*)temp_a, image->iw, image->ih, 0);
	free(temp_r);
	free(temp_g);
	free(temp_b);
	free(temp_a);
	return 1;
}

static void chart_planes_update(Image const *image, ArrayHandle *cpuv, unsigned *gpuv)
{
	if(image->iw*image->ih>1024*1024)
		return;
	int nv=image->iw*image->ih*3*6, nf=nv*5;//subpixel count * 6 vertices
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};
	for(int ky=0, kv=0;ky<image->ih;++ky)
	{
		for(int kx=0;kx<image->iw;++kx)
		{
			int idx=(image->iw*ky+kx)<<2;
			//int idx=3*(iw*ky+kx);
			for(int kc=0;kc<3;++kc)
			{
				int val=image->data[idx|kc];
				//unsigned char val=im[idx+kc]>>(kc<<3)&0xFF;
				
				//planes
				float
					x1=-(image->iw-1-(kx+(float)kc/3)),
					x2=-(image->iw-1-(kx+(float)(kc+0.9f)/3)),
					y1=-(float)ky,
					y2=-((float)ky+0.9f),
					z=(float)val*pixel_amplitude/(nlevels[kc]-1);
				float
					tx=(3*kx+kc+0.5f)/(3*image->iw),
					ty=(ky+0.5f)/image->ih;
				
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;
				vertices[kv++]=x1, vertices[kv++]=y2, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;
				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;

				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;
				vertices[kv++]=x2, vertices[kv++]=y1, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z, vertices[kv++]=tx, vertices[kv++]=ty;
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
static void chart_mesh_update(Image const *image, ArrayHandle *cpuv, unsigned *gpuv)
{
	if(image->iw*image->ih>1024*1024)
		return;
	int nv=(image->iw-1)*(image->ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 6 vertices * 5 floats
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};
	for(int ky=0, kv=0;ky<image->ih-1;++ky)
	{
		for(int kx=0;kx<image->iw-1;++kx)
		{
			int comp[]=
			{
				image->data[(image->iw* ky   +kx  )<<2  ],
				image->data[(image->iw* ky   +kx  )<<2|1],
				image->data[(image->iw* ky   +kx  )<<2|2],
				image->data[(image->iw* ky   +kx+1)<<2  ],
				image->data[(image->iw*(ky+1)+kx  )<<2  ],
				image->data[(image->iw*(ky+1)+kx  )<<2|1],
				image->data[(image->iw*(ky+1)+kx  )<<2|2],
				image->data[(image->iw*(ky+1)+kx+1)<<2  ],
			};
			for(int kc=0;kc<3;++kc)
			{
				//mesh
				float
					x1=-(image->iw-1-(kx+(float)kc/3)),
					x2=-(image->iw-1-(kx+(float)(kc+0.9f)/3)),
					y1=-((float)ky),
					y2=-((float)ky+0.9f),
					z00=(float)comp[kc  ]*pixel_amplitude/(nlevels[kc]-1),
					z01=(float)comp[kc+1]*pixel_amplitude/(nlevels[kc]-1),
					z10=(float)comp[kc+4]*pixel_amplitude/(nlevels[kc]-1),
					z11=(float)comp[kc+5]*pixel_amplitude/(nlevels[kc]-1);
				float
					tx1=(3*kx+kc+0.5f)/(3*image->iw),
					ty1=(ky+0.5f)/image->ih,
					tx2=(3*kx+kc+0.9f+0.5f)/(3*image->iw),
					ty2=(ky+0.9f+0.5f)/image->ih;
				//float
				//	tx=(3*kx+kc+0.5f)/(3*iw),
				//	ty=(ky+0.5f)/ih;
				
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z00, vertices[kv++]=tx1, vertices[kv++]=ty1;
				vertices[kv++]=x1, vertices[kv++]=y2, vertices[kv++]=z10, vertices[kv++]=tx1, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z11, vertices[kv++]=tx2, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z11, vertices[kv++]=tx2, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y1, vertices[kv++]=z01, vertices[kv++]=tx2, vertices[kv++]=ty1;
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z00, vertices[kv++]=tx1, vertices[kv++]=ty1;
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
static void chart_mesh_sep_update(Image const *image, ArrayHandle *cpuv, unsigned *gpuv)
{
	if(image->iw*image->ih>1024*1024)
		return;
	int nv=(image->iw-1)*(image->ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 2 triangles * 3 vertices * 5 floats
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};
	for(int ky=0, kv=0;ky<image->ih-1;++ky)
	{
		for(int kx=0;kx<image->iw-1;++kx)
		{
			for(int kc=0;kc<3;++kc, kv+=30)
			{
				//mesh
				float
					x1=-(image->iw-1-(float)kx),
					x2=-(image->iw-1-((float)kx+0.9f)),
					y1=-(float)ky,
					y2=-((float)ky+0.9f),
					z00=(float)image->data[(image->iw* ky   +kx  )<<2|kc]*pixel_amplitude/(nlevels[kc]-1)+kc*mesh_separation,
					z01=(float)image->data[(image->iw* ky   +kx+1)<<2|kc]*pixel_amplitude/(nlevels[kc]-1)+kc*mesh_separation,
					z10=(float)image->data[(image->iw*(ky+1)+kx  )<<2|kc]*pixel_amplitude/(nlevels[kc]-1)+kc*mesh_separation,
					z11=(float)image->data[(image->iw*(ky+1)+kx+1)<<2|kc]*pixel_amplitude/(nlevels[kc]-1)+kc*mesh_separation;
				float
					tx1=(kx+0.5f)/image->iw,
					tx2=(kx+0.5f+0.9f)/image->iw,
					ty1=(image->ih*kc+ky+0.5f)/(3*image->ih),
					ty2=(image->ih*kc+ky+0.9f+0.5f)/(3*image->ih);
				
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z00, vertices[kv++]=tx1, vertices[kv++]=ty1;
				vertices[kv++]=x1, vertices[kv++]=y2, vertices[kv++]=z10, vertices[kv++]=tx1, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z11, vertices[kv++]=tx2, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y2, vertices[kv++]=z11, vertices[kv++]=tx2, vertices[kv++]=ty2;
				vertices[kv++]=x2, vertices[kv++]=y1, vertices[kv++]=z01, vertices[kv++]=tx2, vertices[kv++]=ty1;
				vertices[kv++]=x1, vertices[kv++]=y1, vertices[kv++]=z00, vertices[kv++]=tx1, vertices[kv++]=ty1;
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
static int hist[768], histmax[3];
//static int hist2[768], histmax2[3];
static int *hist_full=0, hist_full_size=0;
static float blockCR[3]={0};
static void chart_hist_update(Image const *image, int x1, int x2, int y1, int y2, int *_hist, int *_histmax, float *CR)
{
	x1=MAXVAR(x1, 0);
	x2=MINVAR(x2, image->iw);
	y1=MAXVAR(y1, 0);
	y2=MINVAR(y2, image->ih);
	int xcount=x2-x1, ycount=y2-y1, count=xcount*ycount;
	double entropy;
	if(count)
	{
		for(int kc=0;kc<3;++kc)
		{
			calc_histogram(image->data, image->iw, image->ih, kc, x1, x2, y1, y2, image->depth[kc], hist_full, _hist+(kc<<8));
			for(int k=0;k<256;++k)
			{
				if(_histmax[kc]<_hist[kc<<8|k])
					_histmax[kc]=_hist[kc<<8|k];
			}
			if(CR)
			{
				entropy=calc_entropy(hist_full, 1<<image->depth[kc], count);
				CR[kc]=(float)(image->src_depth[kc]/entropy);
			}
		}
	}
}
static void chart_dwthist_update(Image const *image, int kc, int kband, int x1, int x2, int y1, int y2)
{
	memset(hist+((size_t)kband<<8), 0, 256LL*sizeof(int));
	x1=MAXVAR(x1, 0);
	x2=MINVAR(x2, image->iw);
	y1=MAXVAR(y1, 0);
	y2=MINVAR(y2, image->ih);
	int xcount=x2-x1, ycount=y2-y1, count=xcount*ycount;
	if(count)
	{
		calc_histogram(image->data, image->iw, image->ih, kc, x1, x2, y1, y2, image->depth[kc], hist_full, hist+(kc<<8));
		histmax[kband]=0;
		for(int k=0;k<256;++k)
		{
			if(histmax[kband]<hist[kband<<8|k])
				histmax[kband]=hist[kband<<8|k];
		}
		double entropy=calc_entropy(hist_full, image->depth[kc], count);
		blockCR[kband]=(float)(image->src_depth[kc]/entropy);
	}
}
static void move_box_in_window(int boxcenter, int boxsize, int wx1, int wx2, int *ret_boxstart, int *ret_boxend)
{
	if(wx2-wx1<boxsize)
	{
		*ret_boxstart=wx1;
		*ret_boxend=wx2;
	}
	else
	{
		wx1+=boxsize>>1;
		wx2+=(boxsize>>1)-boxsize;
		boxcenter=CLAMP(wx1, boxcenter, wx2);
		*ret_boxstart=boxcenter-(boxsize>>1);
		*ret_boxend=*ret_boxstart+boxsize;
	}
}
static void jhc_getboxbounds(int xs, int ys, int dx, int dy, int iw, int ih, int *bounds)//bounds: {x1, x2, y1, y2}
{
	int ix=(int)(((long long)xs*iw+(w>>1))/w);
	int iy=(int)(((long long)ys*ih+(h>>1))/h);
	move_box_in_window(ix, dx, 0, iw, bounds+0, bounds+1);
	move_box_in_window(iy, dy, 0, ih, bounds+2, bounds+3);
}
static void jh_calchist(int *jhist, int nbits, Image const *image, int x1, int x2, int y1, int y2)
{
	int half[]=
	{
		1<<(image->depth[0]-1),
		1<<(image->depth[1]-1),
		1<<(image->depth[2]-1),
	};
	x1=CLAMP(0, x1, image->iw);
	x2=CLAMP(0, x2, image->iw);
	y1=CLAMP(0, y1, image->ih);
	y2=CLAMP(0, y2, image->ih);
	switch(space_not_color)
	{
	case 0://show correlation in color
		for(int ky=y1;ky<y2;++ky)
		{
			for(int kx=x1;kx<x2;++kx)
			{
				int
					r=(image->data[(image->iw*ky+kx)<<2|0]+half[0])<<nbits>>image->depth[0],
					g=(image->data[(image->iw*ky+kx)<<2|1]+half[1])<<nbits>>image->depth[1],
					b=(image->data[(image->iw*ky+kx)<<2|2]+half[2])<<nbits>>image->depth[2];
				r=CLAMP(0, r, (1<<nbits)-1);
				g=CLAMP(0, g, (1<<nbits)-1);
				b=CLAMP(0, b, (1<<nbits)-1);
				int idx=(b<<nbits|g)<<nbits|r;

				++jhist[idx];
			}
		}
		break;
	case 1://show correlation in space x (CURR, W, WW)
		for(int ky=y1;ky<y2;++ky)
		{
			for(int kx=x1;kx<x2;++kx)
			{
				unsigned char
					v2=kx-2>=0?(unsigned char)((image->data[(image->iw*ky+kx-2)<<2|1]+half[1])<<nbits>>image->depth[1]):0,//WW
					v1=kx-1>=0?(unsigned char)((image->data[(image->iw*ky+kx-1)<<2|1]+half[1])<<nbits>>image->depth[1]):0,//W
					v0=kx-0>=0?(unsigned char)((image->data[(image->iw*ky+kx-0)<<2|1]+half[1])<<nbits>>image->depth[1]):0;//curr
				v0=CLAMP(0, v0, (1<<nbits)-1);
				v1=CLAMP(0, v1, (1<<nbits)-1);
				v2=CLAMP(0, v2, (1<<nbits)-1);
				int idx=(v2<<nbits|v1)<<nbits|v0;

				++jhist[idx];
			}
		}
		break;
	case 2://show correlation in space x (CURR, N, NN)
		for(int ky=y1;ky<y2;++ky)
		{
			for(int kx=x1;kx<x2;++kx)
			{
				unsigned char
					v2=ky-2>=0?(unsigned char)((image->data[(image->iw*(ky-2)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0,//NN
					v1=ky-1>=0?(unsigned char)((image->data[(image->iw*(ky-1)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0,//N
					v0=ky-0>=0?(unsigned char)((image->data[(image->iw*(ky-0)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0;//curr
				v0=CLAMP(0, v0, (1<<nbits)-1);
				v1=CLAMP(0, v1, (1<<nbits)-1);
				v2=CLAMP(0, v2, (1<<nbits)-1);
				int idx=(v2<<nbits|v1)<<nbits|v0;

				++jhist[idx];
			}
		}
		break;
	case 3://show correlation in space x (CURR, N, W)
		for(int ky=y1;ky<y2;++ky)
		{
			for(int kx=x1;kx<x2;++kx)
			{
				unsigned char
					v2=kx-1>=0?(unsigned char)((image->data[(image->iw* ky   +kx-1)<<2|1]+half[1])<<nbits>>image->depth[1]):0,//W
					v1=ky-1>=0?(unsigned char)((image->data[(image->iw*(ky-1)+kx  )<<2|1]+half[1])<<nbits>>image->depth[1]):0,//N
					v0=        (unsigned char)((image->data[(image->iw* ky   +kx  )<<2|1]+half[1])<<nbits>>image->depth[1])  ;//curr
				v0=CLAMP(0, v0, (1<<nbits)-1);
				v1=CLAMP(0, v1, (1<<nbits)-1);
				v2=CLAMP(0, v2, (1<<nbits)-1);
				int idx=(v2<<nbits|v1)<<nbits|v0;

				++jhist[idx];
			}
		}
		break;
	}
}
static void jhc_marchingcubes(ArrayHandle *edges, const int *data, int gx, int gy, int gz, float level, float cubesize)
{
#define VERTIDX(Z, Y, X) 3*(3*Z+Y)+X
	static const char triangles[][3]=
	{
		//6 triangles per cell, reduced resolution - SKEWED
#if 1
		{VERTIDX(0, 0, 0), VERTIDX(0, 0, 2), VERTIDX(2, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(0, 2, 2), VERTIDX(2, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(0, 2, 0), VERTIDX(2, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(2, 2, 0), VERTIDX(2, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(2, 0, 0), VERTIDX(2, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(2, 0, 2), VERTIDX(2, 2, 2)},

		{VERTIDX(0, 0, 0), VERTIDX(0, 0, 2), VERTIDX(0, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(0, 2, 0), VERTIDX(0, 2, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(0, 2, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(0, 0, 0), VERTIDX(2, 0, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(0, 0, 0), VERTIDX(2, 0, 0), VERTIDX(2, 0, 2)},
		{VERTIDX(0, 0, 0), VERTIDX(0, 0, 2), VERTIDX(2, 0, 2)},
#endif

		//96 triangles per cell, full resolution - BUGGY
#if 0
		{VERTIDX(1, 1, 0), VERTIDX(0, 1, 0), VERTIDX(0, 0, 0)},//x==0
		{VERTIDX(1, 1, 0), VERTIDX(0, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(2, 1, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(2, 1, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 0, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 2, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 0, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 2, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 0)},//y==0
		{VERTIDX(1, 0, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 2), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 2), VERTIDX(2, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 0)},//z==0
		{VERTIDX(0, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 2, 1), VERTIDX(0, 2, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 2, 1), VERTIDX(0, 2, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 2), VERTIDX(0, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 2), VERTIDX(0, 2, 2)},

		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 0, 1)},//x==1
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(0, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(2, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 1, 0)},//y==1
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(0, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(2, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(1, 0, 0)},//z==1
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(1, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(1, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(1, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(1, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(1, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(1, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(1, 2, 2)},
		
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(2, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 0), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 0), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 0), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 0), VERTIDX(2, 2, 2)},
		
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 0), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 0), VERTIDX(2, 2, 2)},
		
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 2, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 0, 1), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 0, 1), VERTIDX(2, 2, 2)},
#endif
	};
#undef  VERTIDX
	if(!*edges)
		ARRAY_ALLOC(float[10], *edges, 0, 0, 0, 0);
	edges[0]->count=0;
	__m256 mlevel=_mm256_set1_ps(level);
	float gains[]=
	{
		-cubesize/gx,
		-cubesize/gy,
		cubesize/gz,
	};
	int psize=gx*gy;
	for(int kz=1;kz<gz;++kz)
	{
		for(int ky=1;ky<gy;++ky)
		{
			for(int kx=1;kx<gx;++kx)
			{
				int idx=gx*(gy*kz+ky)+kx;
				float//	vZYX
					v000=(float)data[idx-psize-gx-1],//8 vertices
					v001=(float)data[idx-psize-gx  ],
					v010=(float)data[idx-psize   -1],
					v011=(float)data[idx-psize     ],
					v100=(float)data[idx      -gx-1],
					v101=(float)data[idx      -gx  ],
					v110=(float)data[idx         -1],
					v111=(float)data[idx           ];

				__m256 vec=_mm256_set_ps(v111, v110, v101, v100, v011, v010, v001, v000);
				int cond=_mm256_movemask_ps(_mm256_cmp_ps(vec, mlevel, _CMP_LT_OQ));
				if(!cond||cond==0xFF)//skip if all vertices are below/above level
					continue;
#if 0
				float
					v0mm=(v000+v001+v010+v011)*0.25f,
					v1mm=(v100+v101+v110+v111)*0.25f,
					vm0m=(v000+v001+v100+v101)*0.25f,
					vm1m=(v010+v011+v110+v111)*0.25f,
					vmm0=(v000+v010+v100+v110)*0.25f,
					vmm1=(v001+v011+v101+v111)*0.25f,
					vmmm=(v0mm+v1mm)*0.5f;
				float comps[][3]=
				{
					{(float)(kx-1), (float)kx-0.5f, (float)kx},
					{(float)(ky-1), (float)ky-0.5f, (float)ky},
					{(float)(kz-1), (float)kz-0.5f, (float)kz},
				};
				float vertices[][4]=
				{
					{comps[0][0], comps[1][0], comps[2][0], v000},
					{comps[0][2], comps[1][0], comps[2][0], v001},
					{comps[0][1], comps[1][1], comps[2][0], v0mm},
					{comps[0][0], comps[1][2], comps[2][0], v010},
					{comps[0][2], comps[1][2], comps[2][0], v011},
					{comps[0][1], comps[1][0], comps[2][1], vm0m},
					{comps[0][0], comps[1][1], comps[2][1], vmm0},
					{comps[0][1], comps[1][1], comps[2][1], vmmm},
					{comps[0][2], comps[1][1], comps[2][1], vmm1},
					{comps[0][1], comps[1][2], comps[2][1], vm1m},
					{comps[0][0], comps[1][0], comps[2][2], v100},
					{comps[0][2], comps[1][0], comps[2][2], v101},
					{comps[0][1], comps[1][1], comps[2][2], v1mm},
					{comps[0][0], comps[1][2], comps[2][2], v110},
					{comps[0][2], comps[1][2], comps[2][2], v111},
				};
#endif
#if 1
				float
					v00m=(v000+v001)*0.5f,//12 edge-bisectors
					v01m=(v010+v011)*0.5f,
					v10m=(v100+v101)*0.5f,
					v11m=(v110+v111)*0.5f,
					v0m0=(v000+v010)*0.5f,
					v0m1=(v001+v011)*0.5f,
					v1m0=(v100+v110)*0.5f,
					v1m1=(v101+v111)*0.5f,
					vm00=(v000+v100)*0.5f,
					vm01=(v001+v101)*0.5f,
					vm10=(v010+v110)*0.5f,
					vm11=(v011+v111)*0.5f,
					
					v0mm=(v00m+v01m)*0.5f,//6 face-centers
					v1mm=(v10m+v11m)*0.5f,
					vmm0=(v0m0+v1m0)*0.5f,
					vmm1=(v0m1+v1m1)*0.5f,
					vm0m=(vm00+vm01)*0.5f,
					vm1m=(vm10+vm11)*0.5f,
					
					vmmm=(v0mm+v1mm)*0.5f;//1 cube center

				float comps[][3]=
				{
					{(float)(kx-1), (float)kx-0.5f, (float)kx},
					{(float)(ky-1), (float)ky-0.5f, (float)ky},
					{(float)(kz-1), (float)kz-0.5f, (float)kz},
				};
				float vertices[][4]=
				{
					{comps[0][0], comps[1][0], comps[2][0], v000},
					{comps[0][1], comps[1][0], comps[2][0], v00m},
					{comps[0][2], comps[1][0], comps[2][0], v001},
					{comps[0][0], comps[1][1], comps[2][0], v0m0},
					{comps[0][1], comps[1][1], comps[2][0], v0mm},
					{comps[0][2], comps[1][1], comps[2][0], v0m1},
					{comps[0][0], comps[1][2], comps[2][0], v010},
					{comps[0][1], comps[1][2], comps[2][0], v01m},
					{comps[0][2], comps[1][2], comps[2][0], v011},
					{comps[0][0], comps[1][0], comps[2][1], vm00},
					{comps[0][1], comps[1][0], comps[2][1], vm0m},
					{comps[0][2], comps[1][0], comps[2][1], vm01},
					{comps[0][0], comps[1][1], comps[2][1], vmm0},
					{comps[0][1], comps[1][1], comps[2][1], vmmm},
					{comps[0][2], comps[1][1], comps[2][1], vmm1},
					{comps[0][0], comps[1][2], comps[2][1], vm10},
					{comps[0][1], comps[1][2], comps[2][1], vm1m},
					{comps[0][2], comps[1][2], comps[2][1], vm11},
					{comps[0][0], comps[1][0], comps[2][2], v100},
					{comps[0][1], comps[1][0], comps[2][2], v10m},
					{comps[0][2], comps[1][0], comps[2][2], v101},
					{comps[0][0], comps[1][1], comps[2][2], v1m0},
					{comps[0][1], comps[1][1], comps[2][2], v1mm},
					{comps[0][2], comps[1][1], comps[2][2], v1m1},
					{comps[0][0], comps[1][2], comps[2][2], v110},
					{comps[0][1], comps[1][2], comps[2][2], v11m},
					{comps[0][2], comps[1][2], comps[2][2], v111},
				};
				for(int kt=0;kt<_countof(triangles);++kt)
				{
					const char (*tr)[3]=triangles+kt;
					float (*v[])[4]=
					{
						vertices+tr[0][0],
						vertices+tr[0][1],
						vertices+tr[0][2],
						0,
					};
					if(v[0][0][3]>v[1][0][3])SWAPVAR(v[0], v[1], v[3]);
					if(v[0][0][3]>v[2][0][3])SWAPVAR(v[0], v[2], v[3]);
					if(v[1][0][3]>v[2][0][3])SWAPVAR(v[1], v[2], v[3]);
					if(v[0][0][3]<level&&level<v[2][0][3])
					{
						float *edge=(float*)ARRAY_APPEND(*edges, 0, 1, 1, 0);
						float t=(level-v[0][0][3])/(v[2][0][3]-v[0][0][3]);
						*edge++=MIX(v[0][0][0], v[2][0][0], t)*gains[0];
						*edge++=MIX(v[0][0][1], v[2][0][1], t)*gains[1];
						*edge++=MIX(v[0][0][2], v[2][0][2], t)*gains[2];
						*edge++=0;
						*edge++=0;
						
						int mid=level>v[1][0][3];
						t=(level-v[mid][0][3])/(v[mid+1][0][3]-v[mid][0][3]);
						*edge++=MIX(v[mid][0][0], v[mid+1][0][0], t)*gains[0];
						*edge++=MIX(v[mid][0][1], v[mid+1][0][1], t)*gains[1];
						*edge++=MIX(v[mid][0][2], v[mid+1][0][2], t)*gains[2];
						*edge++=0;
						*edge++=0;
					}
				}
#endif
			}
		}
	}
}
static void chart_jointhist_update(Image const *image, unsigned *txid)
{
	int nlevels=1<<jointhist_nbits, hsize=nlevels*nlevels*nlevels;

	//jointhistogram(image, iw, ih, jointhist_nbits, &jointhist, space_not_color);
	if(jointhist)
	{
		if((int)jointhist->count<hsize)
			ARRAY_APPEND(jointhist, 0, hsize-jointhist->count, 1, 0);
	}
	else
		ARRAY_ALLOC(int, jointhist, 0, hsize, 0, 0);
	int *jhist=(int*)jointhist->data;
#if 1
	memset(jhist, 0, jointhist->esize*jointhist->count);
	int bounds[4]={0};
	jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, image->iw, image->ih, bounds);
	jh_calchist(jhist, jointhist_nbits, image, bounds[0], bounds[1], bounds[2], bounds[3]);
	jhc_marchingcubes(&jhc_mesh, jhist, 1<<jointhist_nbits, 1<<jointhist_nbits, 1<<jointhist_nbits, jhc_level, jh_cubesize);
	if(!jhc_gpubuf)
	{
		glGenBuffers(1, &jhc_gpubuf);
		if(!jhc_gpubuf)
			LOG_ERROR("Alloc error");
	}
	glBindBuffer(GL_ARRAY_BUFFER, jhc_gpubuf);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, (int)(jhc_mesh->count*jhc_mesh->esize), jhc_mesh->data, GL_STATIC_DRAW);	GL_CHECK(error);
#endif
	memset(jhist, 0, jointhist->esize*jointhist->count);
	jh_calchist(jhist, jointhist_nbits, image, 0, image->iw, 0, image->ih);
	int _histmax=0;
	for(int k=0;k<hsize;++k)
	{
		if(!_histmax||_histmax<jhist[k])
			_histmax=jhist[k];
	}
	if(_histmax)
	{
		for(int k=0;k<hsize;++k)
			jointhist->data[k]=(unsigned char)(jhist[k]*255/_histmax);
	}
	
	if(!*txid)
		glGenTextures(nlevels, txid);
	for(int k=0;k<nlevels;++k)
	{
		const unsigned char *hk=jointhist->data+((size_t)k<<(jointhist_nbits<<1));//k*(2^nbits)^2
		send_texture_pot_grey(txid[k], hk, nlevels, nlevels, 1);
		//printf("%d", hk[0]);
	}
	//send_texture_pot_int16x1(*txid, (unsigned*)jointhist->data, nlevels, nlevels*nlevels, 1);
}

#if 0
int bayes_kc=0;
void bayes_update()
{
	float yoffset=tdy*3, half=blocksize*0.5f;
	float x1=blockmx-half, x2=blockmx+half, y1=blockmy-yoffset-half, y2=blockmy-yoffset+half;
	if(blockmx>=0&&blockmx<im1->iw&&blockmy>=yoffset&&blockmy<yoffset+im1->ih)
	{
		if(im1->iw<blocksize)
			x1=0, x2=(float)im1->iw;
		if(x1<0)
			x1=0, x2=(float)blocksize;
		if(x2>im1->iw)
			x1=(float)(im1->iw-blocksize), x2=(float)im1->iw;
				
		if(im1->ih<blocksize)
			y1=0, y2=(float)im1->ih;
		if(y1<0)
			y1=0, y2=(float)blocksize;
		if(y2>im1->ih)
			y1=(float)(im1->ih-blocksize), y2=(float)im1->ih;

		bayes_estimate(im1, im1->iw, im1->ih, (int)roundf(x1), (int)roundf(x2), (int)roundf(y1), (int)roundf(y2), bayes_kc);//?
	}
}
#endif

void apply_selected_transforms(Image *image, int rct_only)
{
	if(!transforms)
		return;
	for(int k=0;k<(int)transforms->count;++k)
	{
		unsigned char tid=transforms->data[k];
		if(rct_only&&tid>CST_INV_SEPARATOR)
			continue;
		switch(tid)
		{
	//	case CT_FWD_ADAPTIVE:		rct_adaptive((char*)image, iw, ih, 1);			break;
	//	case CT_INV_ADAPTIVE:		rct_adaptive((char*)image, iw, ih, 0);			break;
		case CT_FWD_YCoCg_R:		colortransform_YCoCg_R(image, 1);			break;
		case CT_INV_YCoCg_R:		colortransform_YCoCg_R(image, 0);			break;
		case CT_FWD_YCbCr_R_v1:		colortransform_YCbCr_R_v1(image, 1);			break;
		case CT_INV_YCbCr_R_v1:		colortransform_YCbCr_R_v1(image, 0);			break;
		case CT_FWD_YCbCr_R_v2:		colortransform_YCbCr_R_v2(image, 1);			break;
		case CT_INV_YCbCr_R_v2:		colortransform_YCbCr_R_v2(image, 0);			break;
		case CT_FWD_YCbCr_R_v3:		colortransform_YCbCr_R_v3(image, 1);			break;
		case CT_INV_YCbCr_R_v3:		colortransform_YCbCr_R_v3(image, 0);			break;
		case CT_FWD_YCbCr_R_v4:		colortransform_YCbCr_R_v4(image, 1);			break;
		case CT_INV_YCbCr_R_v4:		colortransform_YCbCr_R_v4(image, 0);			break;
		case CT_FWD_YCbCr_R_v5:		colortransform_YCbCr_R_v5(image, 1);			break;
		case CT_INV_YCbCr_R_v5:		colortransform_YCbCr_R_v5(image, 0);			break;
		case CT_FWD_YCbCr_R_v6:		colortransform_YCbCr_R_v6(image, 1);			break;
		case CT_INV_YCbCr_R_v6:		colortransform_YCbCr_R_v6(image, 0);			break;
		case CT_FWD_YCbCr_R_v7:		colortransform_YCbCr_R_v7(image, 1);			break;
		case CT_INV_YCbCr_R_v7:		colortransform_YCbCr_R_v7(image, 0);			break;
	//	case CT_FWD_CrCgCb:		colortransform_CrCgCb_R(image, 1);			break;
	//	case CT_INV_CrCgCb:		colortransform_CrCgCb_R(image, 0);			break;
		case CT_FWD_Pei09:		colortransform_Pei09(image, 1);				break;
		case CT_INV_Pei09:		colortransform_Pei09(image, 0);				break;
		case CT_FWD_JPEG2000:		colortransform_JPEG2000(image, 1);			break;
		case CT_INV_JPEG2000:		colortransform_JPEG2000(image, 0);			break;
		case CT_FWD_SUBGREEN:		colortransform_subtractgreen(image, 1);			break;
		case CT_INV_SUBGREEN:		colortransform_subtractgreen(image, 0);			break;
	//	case CT_FWD_YRGB_v1:		rct_yrgb_v1(image, 1);					break;
	//	case CT_INV_YRGB_v1:		rct_yrgb_v1(image, 0);					break;
	//	case CT_FWD_YRGB_v2:		rct_yrgb_v2(image, 1);					break;
	//	case CT_INV_YRGB_v2:		rct_yrgb_v2(image, 0);					break;
	//	case CT_FWD_CMYK:		ct_cmyk_fwd(image);					break;
	//	case CT_INV_CMYK_DUMMY:									break;
		case CT_FWD_YCbCr:		colortransform_lossy_YCbCr(image, 1);			break;
		case CT_INV_YCbCr:		colortransform_lossy_YCbCr(image, 0);			break;
		case CT_FWD_XYB:		colortransform_lossy_XYB(image, 1);			break;
		case CT_INV_XYB:		colortransform_lossy_XYB(image, 0);			break;
	//	case CT_FWD_MATRIX:		colortransform_lossy_matrix(image, 1);			break;
	//	case CT_INV_MATRIX:		colortransform_lossy_matrix(image, 0);			break;
	//	case CT_FWD_JPEG2000:		colortransform_jpeg2000_fwd((char*)image, iw, ih);	break;
	//	case CT_INV_JPEG2000:		colortransform_jpeg2000_inv((char*)image, iw, ih);	break;
	//	case CT_FWD_XGZ:		colortransform_xgz_fwd((char*)image, iw, ih);		break;
	//	case CT_INV_XGZ:		colortransform_xgz_inv((char*)image, iw, ih);		break;
	//	case CT_FWD_XYZ:		colortransform_xyz_fwd((char*)image, iw, ih);		break;
	//	case CT_INV_XYZ:		colortransform_xyz_inv((char*)image, iw, ih);		break;
	//	case CT_FWD_EXP:		colortransform_exp_fwd((char*)image, iw, ih);		break;
	//	case CT_INV_EXP:		colortransform_exp_inv((char*)image, iw, ih);		break;
	//	case CT_FWD_ADAPTIVE:		colortransform_adaptive((char*)image, iw, ih, 1);	break;
	//	case CT_INV_ADAPTIVE:		colortransform_adaptive((char*)image, iw, ih, 0);	break;
		case CT_FWD_CUSTOM:		rct_custom(image, 1, rct_custom_params);		break;
		case CT_INV_CUSTOM:		rct_custom(image, 0, rct_custom_params);		break;
	//	case CT_FWD_QUAD:		colortransform_quad((char*)image, iw, ih, 1);		break;
	//	case CT_INV_QUAD:		colortransform_quad((char*)image, iw, ih, 1);		break;

	//	case ST_PREPROC_GRAD:		preproc_grad((char*)image, iw, ih);			break;
	//	case ST_PREPROC_X:		preproc_x((char*)image, iw, ih);			break;
	//	case ST_PREPROC_X2:		preproc_x2((char*)image, iw, ih);			break;

	//	case ST_FWD_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 1);break;
	//	case ST_INV_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 0);break;
	//	case ST_FWD_CFL:		pred_cfl((char*)image, iw, ih, 1);			break;
	//	case ST_INV_CFL:		pred_cfl((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_HYBRID3:		pred_hybrid_fwd((char*)image, iw, ih);			break;
	//	case ST_INV_HYBRID3:		pred_hybrid_inv((char*)image, iw, ih);			break;
				
	//	case ST_FWD_CUSTOM2:		custom2_apply((char*)image, iw, ih, 1, &c2_params);	break;
	//	case ST_INV_CUSTOM2:		custom2_apply((char*)image, iw, ih, 0, &c2_params);	break;
		case ST_FWD_CUSTOM3:		custom3_apply(image, 1, pred_ma_enabled, &c3_params);	break;
		case ST_INV_CUSTOM3:		custom3_apply(image, 0, pred_ma_enabled, &c3_params);	break;
		case ST_FWD_CUSTOM:		pred_custom(image, 1, pred_ma_enabled, custom_params);	break;
		case ST_INV_CUSTOM:		pred_custom(image, 0, pred_ma_enabled, custom_params);	break;
	//	case ST_FWD_CUSTOM4:		custom4_apply((char*)image, iw, ih, 1, &c4_params);	break;
	//	case ST_INV_CUSTOM4:		custom4_apply((char*)image, iw, ih, 0, &c4_params);	break;
	//	case ST_FWD_KALMAN:		kalman_apply((char*)image, iw, ih, 1);			break;
	//	case ST_INV_KALMAN:		kalman_apply((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_CUSTOM2:		pred_custom2_apply((char*)image, iw, ih, 1);		break;
	//	case ST_INV_CUSTOM2:		pred_custom2_apply((char*)image, iw, ih, 0);		break;
	//	case ST_FWD_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 1);break;
	//	case ST_INV_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 0);break;
	//	case ST_FWD_LEARNED:		pred_learned_v4((char*)image, iw, ih, 1);		break;
	//	case ST_INV_LEARNED:		pred_learned_v4((char*)image, iw, ih, 0);		break;
#ifdef ALLOW_OPENCL
	//	case ST_FWD_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 1);			break;
	//	case ST_INV_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 0);			break;
#endif
	//	case ST_FWD_DIFF2D:		pred_diff2d_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_DIFF2D:		pred_diff2d_inv((char*)image, iw, ih, 3, 4);		break;
	//	case ST_FWD_HPF:		pred_hpf_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_HPF:		pred_hpf_inv((char*)image, iw, ih, 3, 4);		break;
	//	case ST_FWD_GRAD2:		pred_grad2_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_GRAD2:		pred_grad2_inv((char*)image, iw, ih, 3, 4);		break;
	//	case ST_FWD_ADAPTIVE:		pred_adaptive((char*)image, iw, ih, 3, 4, 1);		break;
	//	case ST_INV_ADAPTIVE:		pred_adaptive((char*)image, iw, ih, 3, 4, 0);		break;
		case ST_FWD_CALIC:		pred_calic(image, 1, pred_ma_enabled);			break;
		case ST_INV_CALIC:		pred_calic(image, 0, pred_ma_enabled);			break;
	//	case ST_FWD_NBLIC:		pred_nblic((char*)image, iw, ih, 1);			break;
	//	case ST_INV_NBLIC:		pred_nblic((char*)image, iw, ih, 0);			break;
		case ST_FWD_WP:			pred_jxl_apply(image, 1, pred_ma_enabled, jxlparams_i16);break;
		case ST_INV_WP:			pred_jxl_apply(image, 0, pred_ma_enabled, jxlparams_i16);break;
		case ST_FWD_WPU:		pred_WPU(image, 1);					break;
		case ST_INV_WPU:		pred_WPU(image, 0);					break;
		case ST_FWD_DEFERRED:		pred_wp_deferred(image, 1);				break;
		case ST_INV_DEFERRED:		pred_wp_deferred(image, 0);				break;
		case ST_FWD_MM:			pred_w2_apply(image, 1, pred_ma_enabled, pw2_params);	break;
		case ST_INV_MM:			pred_w2_apply(image, 0, pred_ma_enabled, pw2_params);	break;
		case ST_FWD_OLS:		pred_ols(image, 1, pred_ma_enabled);			break;
		case ST_INV_OLS:		pred_ols(image, 0, pred_ma_enabled);			break;
		case ST_FWD_OLS2:		pred_ols2(image, 1, pred_ma_enabled);			break;
		case ST_INV_OLS2:		pred_ols2(image, 0, pred_ma_enabled);			break;
		case ST_FWD_OLS3:		pred_ols3(image, 1, pred_ma_enabled);			break;
		case ST_INV_OLS3:		pred_ols3(image, 0, pred_ma_enabled);			break;
		case ST_FWD_OLS4:		pred_ols4(image, ols4_period, ols4_lr, ols4_mask[0], ols4_mask[1], ols4_mask[2], ols4_mask[3], 1);break;
		case ST_INV_OLS4:		pred_ols4(image, ols4_period, ols4_lr, ols4_mask[0], ols4_mask[1], ols4_mask[2], ols4_mask[3], 0);break;
		case ST_FWD_OLS5:		pred_ols5(image, 1);					break;
		case ST_INV_OLS5:		pred_ols5(image, 0);					break;
		case ST_FWD_PACKSIGN:		packsign(image, 1);					break;
		case ST_INV_PACKSIGN:		packsign(image, 0);					break;
		case ST_FWD_PU:			pred_PU(image, 1);					break;
		case ST_INV_PU:			pred_PU(image, 0);					break;
		case ST_FWD_CG3D:		pred_CG3D(image, 1, pred_ma_enabled);			break;
		case ST_INV_CG3D:		pred_CG3D(image, 0, pred_ma_enabled);			break;
		case ST_FWD_CLAMPGRAD:		pred_clampgrad(image, 1, pred_ma_enabled);		break;
		case ST_INV_CLAMPGRAD:		pred_clampgrad(image, 0, pred_ma_enabled);		break;
		case ST_FWD_AV2:		pred_av2(image, 1);					break;
		case ST_INV_AV2:		pred_av2(image, 0);					break;
		case ST_FWD_MEDIAN:		pred_median(image, 1);					break;
		case ST_INV_MEDIAN:		pred_median(image, 0);					break;
	//	case ST_FWD_ECOEFF:		pred_ecoeff(image, 1, pred_ma_enabled);			break;
	//	case ST_INV_ECOEFF:		pred_ecoeff(image, 0, pred_ma_enabled);			break;
	//	case ST_FWD_AVERAGE:		pred_average(image, 1, pred_ma_enabled);		break;
	//	case ST_INV_AVERAGE:		pred_average(image, 0, pred_ma_enabled);		break;
		case ST_FWD_MULTISTAGE:		pred_multistage(image, 1, pred_ma_enabled);		break;
		case ST_INV_MULTISTAGE:		pred_multistage(image, 0, pred_ma_enabled);		break;
	//	case ST_FWD_ZIPPER:		pred_zipper(&image, 1, pred_ma_enabled);		break;
	//	case ST_INV_ZIPPER:		pred_zipper(&image, 0, pred_ma_enabled);		break;
		case ST_FWD_P3:			pred_separate(image, 1, pred_ma_enabled);		break;
		case ST_INV_P3:			pred_separate(image, 0, pred_ma_enabled);		break;
	//	case ST_FWD_DIR:		pred_dir(image, 1, pred_ma_enabled);			break;
	//	case ST_INV_DIR:		pred_dir(image, 0, pred_ma_enabled);			break;
		case ST_FWD_G2:			pred_grad2(image, 1, pred_ma_enabled);			break;
		case ST_INV_G2:			pred_grad2(image, 0, pred_ma_enabled);			break;
		case ST_FWD_T47:		pred_t47(image, 1, pred_ma_enabled);			break;
		case ST_INV_T47:		pred_t47(image, 0, pred_ma_enabled);			break;
		case ST_FWD_DCT4:		image_dct4_fwd(image);					break;
		case ST_INV_DCT4:		image_dct4_inv(image);					break;
		case ST_FWD_DCT8:		image_dct8_fwd(image);					break;
		case ST_INV_DCT8:		image_dct8_inv(image);					break;
#if 0
	//	case ST_FWD_JMJ:		pred_jmj_apply((char*)image, iw, ih, 1);		break;
	//	case ST_INV_JMJ:		pred_jmj_apply((char*)image, iw, ih, 0);		break;
	//	case ST_FWD_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 1);		break;
	//	case ST_INV_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 0);		break;
	//	case ST_FWD_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 1);		break;
	//	case ST_INV_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 0);		break;
	//	case ST_FWD_MEDIAN:		pred_median_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_MEDIAN:		pred_median_inv((char*)image, iw, ih, 3, 4);		break;
	//	case ST_FWD_DCT3PRED:		pred_dct3_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_DCT3PRED:		pred_dct3_inv((char*)image, iw, ih, 3, 4);		break;
	//	case ST_FWD_PATHPRED:		pred_path_fwd((char*)image, iw, ih, 3, 4);		break;
	//	case ST_INV_PATHPRED:		pred_path_inv((char*)image, iw, ih, 3, 4);		break;
		case ST_FWD_GRADPRED:		pred_grad_fwd((char*)image, iw, ih, 3, 4);		break;
		case ST_INV_GRADPRED:		pred_grad_inv((char*)image, iw, ih, 3, 4);		break;
		case ST_FWD_GRAD2:		pred_grad2((char*)image, iw, ih, 1);			break;
		case ST_INV_GRAD2:		pred_grad2((char*)image, iw, ih, 0);			break;
		case ST_FWD_CTX:		pred_ctx((char*)image, iw, ih, 1);			break;
		case ST_INV_CTX:		pred_ctx((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_C03:		pred_c03((char*)image, iw, ih, 1);			break;
	//	case ST_INV_C03:		pred_c03((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_C10:		pred_c10((char*)image, iw, ih, 1);			break;
	//	case ST_INV_C10:		pred_c10((char*)image, iw, ih, 0);			break;
		case ST_FWD_C20:		pred_c20((char*)image, iw, ih, 1);			break;
		case ST_INV_C20:		pred_c20((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_WU97:		pred_wu97((char*)image, iw, ih, 1);			break;
	//	case ST_INV_WU97:		pred_wu97((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_BITWISE:		pred_bitwise((char*)image, iw, ih, 1);			break;
	//	case ST_INV_BITWISE:		pred_bitwise((char*)image, iw, ih, 0);			break;
		case ST_FWD_SHUFFLE:		shuffle((char*)image, iw, ih, 1);			break;
		case ST_INV_SHUFFLE:		shuffle((char*)image, iw, ih, 0);			break;
	//	case ST_FWD_SPLIT:		image_split_fwd((char*)image, iw, ih);			break;
	//	case ST_INV_SPLIT:		image_split_inv((char*)image, iw, ih);			break;
#endif

	//	case ST_FWD_LAZY:
	//	case ST_INV_LAZY:
		case ST_FWD_HAAR:
		case ST_INV_HAAR:
		case ST_FWD_SQUEEZE:
		case ST_INV_SQUEEZE:
		case ST_FWD_LEGALL53:
		case ST_INV_LEGALL53:
		case ST_FWD_CDF97:
		case ST_INV_CDF97:
	//	case ST_FWD_GRAD_DWT:
	//	case ST_INV_GRAD_DWT:
	//	case ST_FWD_EXPDWT:
	//	case ST_INV_EXPDWT:
	//	case ST_FWD_CUSTOM_DWT:
	//	case ST_INV_CUSTOM_DWT:
			{
				ArrayHandle sizes=dwt2d_gensizes(image->iw, image->ih, 3, 3, 0);
				int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
				for(int kc=0;kc<4;++kc)
				{
					if(!im1->depth[kc])
						continue;
					switch(tid)
					{
				//	case ST_FWD_LAZY:      dwt2d_lazy_fwd   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				//	case ST_INV_LAZY:      dwt2d_lazy_inv   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_FWD_HAAR:      dwt2d_haar_fwd   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_INV_HAAR:      dwt2d_haar_inv   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_FWD_SQUEEZE:   dwt2d_squeeze_fwd(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_INV_SQUEEZE:   dwt2d_squeeze_inv(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_FWD_LEGALL53:  dwt2d_cdf53_fwd  (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_INV_LEGALL53:  dwt2d_cdf53_inv  (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_FWD_CDF97:     dwt2d_cdf97_fwd  (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					case ST_INV_CDF97:     dwt2d_cdf97_inv  (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				//	case ST_FWD_GRAD_DWT:  dwt2d_grad_fwd   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				//	case ST_INV_GRAD_DWT:  dwt2d_grad_inv   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				//	case ST_FWD_EXPDWT:    dwt2d_exp_fwd    (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;//TODO use customdwtparams instead of sharing allcustomparam_st
				//	case ST_INV_EXPDWT:    dwt2d_exp_inv    (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
				//	case ST_FWD_CUSTOM_DWT:dwt2d_custom_fwd (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
				//	case ST_INV_CUSTOM_DWT:dwt2d_custom_inv (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
					}
					int inv=tid&1;
					im1->depth[kc]+=(char)(!inv-inv);
					UPDATE_MAX(im1->depth[kc], im1->src_depth[kc]);
				}
				array_free(&sizes);
				free(temp);
			}
			break;
	//	case ST_FWD_DEC_DWT:   dwt2d_dec_fwd((char*)image, iw, ih);	break;
	//	case ST_INV_DEC_DWT:   dwt2d_dec_inv((char*)image, iw, ih);	break;
		}//switch
		//calc_depthfromdata(image->data, image->iw, image->ih, image->depth, image->src_depth);//X  depth must depend only on src_depth and applied RCTs, so that preds can apply MA
	}//for
}
void update_image(void)//apply selected operations on original image, calculate CRs, and export
{
	if(!im0)
		return;
	image_copy(&im1, im0);
	apply_selected_transforms(im1, 0);

	//do not modify im1 beyond this point

	image_export_uint8(im1, &im_export, 1, 0);

	int maxdepth=calc_maxdepth(im1, 0);
	int maxlevels=1<<maxdepth;
	if(hist_full_size<maxlevels)
	{
		void *p=realloc(hist_full, maxlevels*sizeof(int));
		if(!p)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		hist_full=(int*)p;
		hist_full_size=maxlevels;
	}
	calc_csize_stateful(im1, hist_full, ch_entropy);
#if 0
	if(ec_method==ECTX_HIST)
	{
		for(int kc=0;kc<4;++kc)
		{
			if(im1->depth[kc])
			{
				calc_histogram(im1->data, im1->iw, im1->ih, kc, 0, im1->iw, 0, im1->ih, im1->depth[kc], hist_full, 0);
				ch_entropy[kc]=calc_entropy(hist_full, 1<<im1->depth[kc], im1->iw*im1->ih);
			}
			else
				ch_entropy[kc]=0;
		}
		//channel_entropy(image, iw*ih, 3, 4, ch_cr, usage);
	}
	else if(ec_method==ECTX_ABAC)
		calc_csize_abac(im1, ch_entropy);
	else
		calc_csize_ec(im1, ec_method, ec_adaptive?ec_adaptive_threshold:0, ec_expbits, ec_msb, ec_lsb, ch_entropy);
#endif
	
	combCRhist[combCRhist_idx][0]=1/(float)(im1->src_depth[0]/ch_entropy[0]);
	combCRhist[combCRhist_idx][1]=1/(float)(im1->src_depth[0]/ch_entropy[1]);
	combCRhist[combCRhist_idx][2]=1/(float)(im1->src_depth[0]/ch_entropy[2]);
	combCRhist[combCRhist_idx][3]=1/(float)((im1->src_depth[0]+im1->src_depth[1]+im1->src_depth[2]+im1->src_depth[3])/(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3]));
	for(int k=0;k<4;++k)
	{
		if(combCRhist_max==1||combCRhist_max<combCRhist[combCRhist_idx][k])
			combCRhist_max=combCRhist[combCRhist_idx][k];
	}
	combCRhist_idx=(combCRhist_idx+1)%combCRhist_SIZE;

	//if(im1->iw<1024&&im1->ih<1024)
	{
		if(!send_image_separate_subpixels(im1, &txid_separate_r, &txid_separate_g, &txid_separate_b, &txid_separate_a))
			LOG_ERROR("Failed to send texture to GPU");
	}

	switch(mode)
	{
	//case VIS_PLANES:
	//	chart_planes_update(image, iw, ih, &cpu_vertices, &gpu_vertices);
	//	break;
	//case VIS_MESH:
	//	chart_mesh_update(image, iw, ih, &cpu_vertices, &gpu_vertices);
	//	break;
#if 0
	case VIS_MESH_SEPARATE:
		chart_mesh_sep_update(im1, &cpu_vertices, &gpu_vertices);
		break;
#endif
	//case VIS_BAYES:
	//	bayes_update();
	//	break;
	case VIS_ZIPF:
		{
			ptrdiff_t res=(ptrdiff_t)im1->iw*im1->ih;
			float *fbuf=(float*)malloc(res*sizeof(float[4]));
			void *p=realloc(zimage, res*sizeof(char[4]));
			if(!fbuf||!p)
			{
				LOG_ERROR("Allocation error");
				return;
			}
			zimage=(unsigned char*)p;
			for(ptrdiff_t k=0;k<res;++k)
				zimage[k<<2|3]=0xFF;
			for(int kc=0;kc<3;++kc)
			{
				int nlevels=1<<im1->depth[kc];
				calc_histogram(im1->data, im1->iw, im1->ih, kc, 0, im1->iw, 0, im1->ih, im1->depth[kc], hist_full, 0);
				for(ptrdiff_t k=0;k<res;++k)
				{
					int sym=im1->data[k<<2|kc]+(nlevels>>1);
					sym=CLAMP(0, sym, nlevels-1);
					int freq=hist_full[sym];
					if(freq)
					{
						float bitsize=-log2f((float)freq/res);
						fbuf[k<<2|kc]=bitsize;
						//if(vmax<bitsize)
						//	vmax=bitsize;
					}
					else
						fbuf[k<<2|kc]=0;
				}
			}
			float vmax=0;
			for(ptrdiff_t k=0;k<res;++k)
			{
				if(!vmax||vmax<fbuf[k<<2|0])
					vmax=fbuf[k<<2|0];
				if(vmax<fbuf[k<<2|1])
					vmax=fbuf[k<<2|1];
				if(vmax<fbuf[k<<2|2])
					vmax=fbuf[k<<2|2];
			}
			if(vmax)
			{
				for(ptrdiff_t k=0;k<res;++k)
				{
					for(int kc=0;kc<3;++kc)
					{
						float val=fbuf[k<<2|kc]*255/vmax;
						zimage[k<<2|kc]=(unsigned char)CLAMP(0, val, 255);
					}
				}
			}
			else
			{
				for(ptrdiff_t k=0;k<res;++k)
				{
					zimage[k<<2|0]=0;
					zimage[k<<2|1]=0;
					zimage[k<<2|2]=0;
				}
			}
			free(fbuf);
		}
		break;
	case VIS_HISTOGRAM:
		memset(histmax, 0, sizeof(histmax));
		memset(hist, 0, sizeof(hist));
		chart_hist_update(im1, 0, im1->iw, 0, im1->ih, hist, histmax, 0);
		break;
	case VIS_JOINT_HISTOGRAM:
		chart_jointhist_update(im1, txid_jointhist);
		break;
	}
}

static void draw_AAcuboid_wire(float x1, float x2, float y1, float y2, float z1, float z2, int color)
{
	float cuboid[]=
	{
		x1, y1, z1,
		x2, y1, z1,
		x2, y2, z1,
		x1, y2, z1,
		x1, y1, z2,
		x2, y1, z2,
		x2, y2, z2,
		x1, y2, z2,
	};
	draw_3d_line(&cam, cuboid    , cuboid+1*3, color);
	draw_3d_line(&cam, cuboid+1*3, cuboid+2*3, color);
	draw_3d_line(&cam, cuboid+2*3, cuboid+3*3, color);
	draw_3d_line(&cam, cuboid+3*3, cuboid    , color);
		
	draw_3d_line(&cam, cuboid+(4  )*3, cuboid+(4+1)*3, color);
	draw_3d_line(&cam, cuboid+(4+1)*3, cuboid+(4+2)*3, color);
	draw_3d_line(&cam, cuboid+(4+2)*3, cuboid+(4+3)*3, color);
	draw_3d_line(&cam, cuboid+(4+3)*3, cuboid+(4  )*3, color);
		
	draw_3d_line(&cam, cuboid    , cuboid+(4  )*3, color);
	draw_3d_line(&cam, cuboid+1*3, cuboid+(4+1)*3, color);
	draw_3d_line(&cam, cuboid+2*3, cuboid+(4+2)*3, color);
	draw_3d_line(&cam, cuboid+3*3, cuboid+(4+3)*3, color);
}
//void chart_planes_draw()
//{
//	draw_AAcuboid_wire(0, (float)im0->iw, 0, (float)im0->ih, 0, pixel_amplitude, 0xFF000000);
//
//	draw_3D_triangles(&cam, gpu_vertices, 0, (int)(cpu_vertices->count/5), txid_separate_r);
//	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize/3, (int)(cpu_vertices->count/5), txid_separate_g);
//	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize*2/3, (int)(cpu_vertices->count/5), txid_separate_b);
//}
//void chart_mesh_draw()
//{
//	draw_AAcuboid_wire(0, (float)im0->iw, 0, (float)im0->ih, 0, pixel_amplitude, 0xFF000000);
//
//	draw_3D_triangles(&cam, gpu_vertices, 0, (int)(cpu_vertices->count/5), txid_separate_r);
//	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize/3, (int)(cpu_vertices->count/5), txid_separate_g);
//	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize*2/3, (int)(cpu_vertices->count/5), txid_separate_b);
//}
#if 0
void chart_mesh_sep_draw()
{
	float RMSE1=0, RMSE2=0;
	int RMSE_den=0;
	float reach=5;
	float ix=0, iy=0;
	if(extrainfo)
	{
		float *vertices=(float*)cpu_vertices->data;
		ix=im1->iw-1-cam.x, iy=cam.y;
		float xstart=ix-reach, xend=ix+reach, ystart=iy-reach, yend=iy+reach;
		for(int ky=0, kv=0;ky<im1->ih-1;++ky)
		{
			float error[3]={0};
			for(int kx=0;kx<im1->iw-1;++kx)
			{
				for(int kc=0;kc<3;++kc, kv+=30)
				{
					if(!kc&&kx>=xstart&&kx<xend&&ky>=ystart&&ky<yend)
					{
						float
							x1=vertices[kv  ], x2=vertices[kv+10],
							y1=vertices[kv+1], y2=vertices[kv+ 6],
							topleft=vertices[kv+2], top=vertices[kv+22], left=vertices[kv+7], curr=vertices[kv+12];

						float vmin, vmax, pred;
						if(top<left)
							vmin=top, vmax=left;
						else
							vmin=left, vmax=top;
						if(topleft<vmin)//gradient predictor
							pred=vmax;
						else if(topleft>vmax)
							pred=vmin;
						else
							pred=top+left-topleft;

						//if((int)ix==746&&(int)iy==10&&kx==(int)ix&&ky==(int)ky)//
						//	kx=746;//

						float pred2;//=top+left-topleft;	//=topleft+(top-topleft)+(left-topleft)
						int kx2=kx+1, ky2=ky+1, idx=im1->iw*ky2+kx2;
						char
							ctltl    =kx2-2>=0     &&ky2-2>=0?im1->data[(idx-im1->iw*2-2)<<2|kc]-128:0,
							ctt      =kx2  <im1->iw&&ky2-2>=0?im1->data[(idx-im1->iw*2  )<<2|kc]-128:0,
							ctrtr    =kx2+2<im1->iw&&ky2-2>=0?im1->data[(idx-im1->iw*2+2)<<2|kc]-128:0,

							ctopleft =              im1->data[(idx-im1->iw  -1)<<2|kc]-128  ,
							ctop     =kx2  <im1->iw?im1->data[(idx-im1->iw    )<<2|kc]-128:0,
							ctopright=kx2+1<im1->iw?im1->data[(idx-im1->iw  +1)<<2|kc]-128:0,

							cll      =kx2-2>=0     &&ky2<im1->ih?im1->data[(idx     -2)<<2|kc]-128:0,
							cleft    =               ky2<im1->ih?im1->data[(idx     -1)<<2|kc]-128:0,
							ccurr    =kx2  <im1->iw&&ky2<im1->ih?im1->data[ idx        <<2|kc]-128:0;
						int
							g45tl=ctopleft-ctltl,
							gxtl=ctop-ctopleft,
							gyt=ctop-ctt,
							gxtr=ctop-ctopright,
							g45tr=ctopright-ctrtr,
							gxl=cleft-cll,
							gyl=cleft-ctopleft;

						//gradient x1.25		X
#if 0
						pred2=ctopleft;
						if((gyl<0)!=(gxtl<0))//left & top grad have different signs
							pred2+=cleft+ctopleft-(ctopleft<<1);
						else if(gyl<0)//same signs, changing directions
							pred2+=MINVAR(gyl, gxtl)+MAXVAR(gyl, gxtl)*0.25f;
						else
							pred2+=MAXVAR(gyl, gxtl)+MINVAR(gyl, gxtl)*0.25f;
#endif

						//gradient extended		X
#if 0
						pred2=ctopleft;
						if((gyl<0)!=(gxtl<0))//left & top grad have different signs
							pred2+=cleft+ctopleft-(ctopleft<<1);
						else if((gyl<0)==(gyt<0)&&(gxl<0)==(gxtl<0))//same signs, keeping directions
						{
							if(gyl<0)
								pred2+=MAXVAR(gyl, gyt)+MAXVAR(gxl, gxtl);
							else
								pred2+=MINVAR(gyl, gyt)+MINVAR(gxl, gxtl);
						}
						else if(gyl<0)//same signs, changing directions
							pred2+=MINVAR(gyl, gxtl);
						else
							pred2+=MAXVAR(gyl, gxtl);
#endif

#if 0
						pred2=0;
						if((g45tl<0)!=(g45tr<0))//path45
							pred2+=(ctopleft+g45tl+ctopright+g45tr)*0.5f;
						else if(g45tl<0)
							pred2+=MINVAR(ctopleft, ctopright);
						else
							pred2+=MAXVAR(ctopleft, ctopright);
						
						if((gxl<0)!=(gyt<0))//path
							pred2+=(cleft+gxl+ctop+gyt)*0.5f;
						else if(gxl<0)
							pred2+=MINVAR(ctop, cleft);
						else
							pred2+=MAXVAR(ctop, cleft);

						int temp;
						if(gyl<0)//gamma predictor
						{
							if(gxtl<0)
							{
								if(gxtr<0)	//hole
									pred2+=MINVAR(ctop, cleft)-(gxtl+gxtr)*0.5f;
								else		//bottom-right descends
									temp=ctopleft+ctopright, pred2+=MINVAR(temp, cleft<<1)*0.5f;
									//pred2+=(ctopleft+ctopright+ctop+cleft)*0.25f;
									//pred+=ctopright;
							}
							else
							{
								if(gxtr<0)	//bottom-left descends
									pred2+=(ctopleft+ctopright+ctop*2+cleft*2)*(1.f/6);
									//pred2+=ctop+gyl;
								else		//roof, bottom descends
									pred2+=(ctopleft+ctopright)*0.5f+gyl;
							}
						}
						else
						{
							if(gxtl<0)
							{
								if(gxtr<0)	//valley, top descends
									pred2+=(ctopleft+ctopright)*0.5f+gyl;
								else		//top-right descends
									pred2+=(ctopleft+ctopright+ctop*2+cleft*2)*(1.f/6);
									//pred2+=ctop+gyl;
							}
							else
							{
								if(gxtr<0)	//top-left descends
									temp=ctopleft+ctopright, pred2+=MAXVAR(temp, cleft<<1)*0.5f;
									//pred2+=(ctopleft+ctopright+ctop+cleft)*0.25f;
								else		//peak
									pred2+=MAXVAR(ctop, cleft)+(gxtl+gxtr)*0.5f;
							}
						}
						pred2*=1.f/3;
#endif

						int gx=gxl+gxtl, gy=gyl+gyt, T=44;
						if(gy>T+gx)
							pred2=cleft;
						else if(gy+T<gx)
							pred2=ctop;
						else
							pred2=(float)(ctop+cleft-ctopleft);

						//pred2=pred;//plain gradient

						//pred2=pred2*(255-error[kc])/255;

						//if(gyl<0&&gxtl<0&&gxtr<0)
						//	pred2+=MINVAR(ctop, cleft)-(gxtl+gxtr)*0.5f;
						//else if(gyl>0&&gxtl>0&&gxtr>0)
						//	pred2+=MAXVAR(ctop, cleft)+(gxtl+gxtr)*0.5f;

						error[kc]*=(ccurr-pred2)/128;

						pred2+=128;
						pred2=CLAMP(0, pred2, 255);
						pred2*=pixel_amplitude/255;
						//pred2=(ctopleft+g45tl);//linear combination


						float delta=curr-pred;
						RMSE1+=delta*delta;

						delta=curr-pred2;//derivative
						RMSE2+=delta*delta;

						++RMSE_den;

						float v[]=
						{
							-x2, -y1, top,  //0
							-x1, -y2, left, //1
							-x2, -y2, curr, //2
							-x2, -y2, pred, //3
							-x2, -y2, pred2,//4
						};
						int linecolor=0xFF000000|0xFF<<((kc+1)%3<<3);
						draw_3d_line(&cam, v+3*3, v    , linecolor);
						draw_3d_line(&cam, v+3*3, v+3  , linecolor);
						draw_3d_line(&cam, v+3*3, v+3*2, linecolor);

						linecolor=0xFF000000|0xFF<<((kc+2)%3<<3);
						draw_3d_line(&cam, v+3*4, v    , linecolor);
						draw_3d_line(&cam, v+3*4, v+3  , linecolor);
						draw_3d_line(&cam, v+3*4, v+3*2, linecolor);

						//v[3*2  ]*=-1;
						//v[3*2+1]*=-1;
						cam_world2cam(cam, v+3*2, v, v+3);
						if(v[2]>0)
						{
							cam_cam2screen(cam, v, v+3, w>>1, h>>1);
							if(v[3]>=0&&v[3]<w&&v[4]>=0&&v[4]<h)
								GUIPrint(0, v[3], v[4], 1, "%.0f", curr*(255/pixel_amplitude));
								//GUIPrint(0, v[3], v[4], 1, "%.2f-%.2f-%.2f", curr, pred, pred2);
						}
					}
				}
			}
		}
	}
	draw_AAcuboid_wire(0, (float)im1->iw, 0, (float)im1->ih, 0, pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)im1->iw, 0, (float)im1->ih, mesh_separation, mesh_separation+pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)im1->iw, 0, (float)im1->ih, mesh_separation*2, mesh_separation*2+pixel_amplitude, 0xFF000000);

	draw_3D_triangles(&cam, gpu_vertices, 0, (int)(cpu_vertices->count/5), txid_separate_r);
	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize/3, (int)(cpu_vertices->count/5), txid_separate_g);
	draw_3D_triangles(&cam, gpu_vertices, cpu_vertices->count*cpu_vertices->esize*2/3, (int)(cpu_vertices->count/5), txid_separate_b);
	
	if(extrainfo&&RMSE_den)
	{
		RMSE1=sqrtf(RMSE1/RMSE_den)*255/pixel_amplitude;
		RMSE2=sqrtf(RMSE2/RMSE_den)*255/pixel_amplitude;
		float y=(float)(h>>1);
		GUIPrint(0, 0, y, 1, "imXY %10f %10f", ix, iy); y+=tdy;
		int c=set_text_color(0xFF00FF00);
		GUIPrint(0, 0, y, 1, "RMSE_grad %10lf / 255", RMSE1); y+=tdy;
		set_text_color(0xFFFF0000);
		GUIPrint(0, 0, y, 1, "RMSE_diff %10lf", RMSE2);
		set_text_color(c);

		float x=(float)(w>>2);
		y=(float)(h>>1)+tdy;
		draw_rect(x, x+RMSE1*pixel_amplitude, y, y+tdy, 0xA000FF00); y+=tdy;
		draw_rect(x, x+RMSE2*pixel_amplitude, y, y+tdy, 0xA0FF0000);
	}
}
#endif
static void chart_hist_draw(float x1, float x2, float y1, float y2, int cstart, int cend, int color, unsigned char alpha, int *_hist, int *_histmax)
{
	for(int kc=cstart;kc<cend;++kc)
	{
		if(_histmax[kc])
		{
			float dy=(y2-y1)/3.f, histpx=dy/_histmax[kc];
			int k=1;
			float y=k*histpx*10000;
			for(;y<dy;++k)
			{
				draw_line(x1, y1+(kc+1)*dy-y, x2, y1+(kc+1)*dy-y, color?color:alpha<<24|0xFF<<(kc<<3));//0x40
				y=k*histpx*10000;
			}
			for(int k2=0;k2<256;++k2)
				draw_rect(x1+k2*(x2-x1)/256, x1+(k2+1)*(x2-x1)/256, y1+(kc+1)*dy-_hist[kc<<8|k2]*histpx, y1+(kc+1)*dy, color?color:alpha<<24|0xFF<<(kc<<3));//0x80
		}
	}
}
static void chart_hist_draw2(float x1, float x2, float y1, float y2, int color, int *_hist, int _histmax)
{
	if(_histmax==-1)
	{
		_histmax=0;
		for(int sym=0;sym<256;++sym)
		{
			if(_histmax<_hist[sym])
				_histmax=_hist[sym];
		}
	}
	if(!_histmax)
		return;
	float dy=y2-y1, histpx=dy/_histmax;
	int k=1;
	float y=k*histpx*10000;
	for(;y<dy;++k)
	{
		draw_line(x1, y2-y, x2, y2-y, color);
		y=k*histpx*10000;
	}
	for(int k2=0;k2<256;++k2)
		draw_rect(x1+k2*(x2-x1)/256, x1+(k2+1)*(x2-x1)/256, y2-_hist[k2]*histpx, y2, color);
}
#if 0
static void draw_cloud(int x, int y, int blocksize, float cubesize)
{
	static ArrayHandle vertices=0;
	blocksize>>=1;
	int nlevels[]=
	{
		1<<im1->depth[0],
		1<<im1->depth[1],
		1<<im1->depth[2],
	};
	for(int ky=MAXVAR(y-blocksize, 0);ky<MINVAR(y+blocksize, im1->ih);++ky)
	{
		for(int kx=MAXVAR(x-blocksize, 0);kx<MINVAR(x+blocksize, im1->iw);++kx)
		{
			int idx=(im1->iw*ky+kx)<<2;
			int
				r=im1->data[idx|0]+(nlevels[0]>>1),
				g=im1->data[idx|1]+(nlevels[1]>>1),
				b=im1->data[idx|2]+(nlevels[2]>>1);
			draw_3d_line_enqueue(&vertices, (float)r*cubesize/(nlevels[0]-1), (float)g*cubesize/(nlevels[1]-1), (float)b*cubesize/(nlevels[2]-1));
		}
		int color=0xFF0000FF;
		int texture[]={color, color, color, color};
		draw_3d_flush(vertices, &cam, texture, 2, 2, 0, GL_LINE_STRIP);
	}
	for(int kx=MAXVAR(x-blocksize, 0);kx<MINVAR(x+blocksize, im1->iw);++kx)
	{
		for(int ky=MAXVAR(y-blocksize, 0);ky<MINVAR(y+blocksize, im1->ih);++ky)
		{
			int idx=(im1->iw*ky+kx)<<2;
			int
				r=im1->data[idx|0]+(nlevels[0]>>1),
				g=im1->data[idx|1]+(nlevels[1]>>1),
				b=im1->data[idx|2]+(nlevels[2]>>1);
			draw_3d_line_enqueue(&vertices, (float)r*cubesize/(nlevels[0]-1), (float)g*cubesize/(nlevels[1]-1), (float)b*cubesize/(nlevels[2]-1));
		}
		int color=0xFFFF0000;
		int texture[]={color, color, color, color};
		draw_3d_flush(vertices, &cam, texture, 2, 2, 0, GL_LINE_STRIP);
	}
}
#endif
static void chart_jointhist_draw(void)
{
	draw_AAcuboid_wire(0, jh_cubesize, 0, jh_cubesize, 0, jh_cubesize, 0xFF000000);

	//for(int k=0;k<(int)jhc_mesh->count;++k)//
	//{
	//	float *vertices=(float*)array_at(&jhc_mesh, k);
	//	draw_3d_line(&cam, vertices, vertices+5, 0x80FF00FF);
	//}
	draw_3d_wireframe_gpu(&cam, jhc_gpubuf, (int)(jhc_mesh->count<<1), 0xC0000000, GL_LINES);

#if 0
	draw_cloud(jhx, jhy, blocksize, cubesize);
	GUIPrint(0, (float)(w>>1), tdy*2, 1, "(%d, %d) size %f", jhx, jhy, blocksize);
	const int hstep=2;
	jhx+=hstep;
	if(jhx>=im1->iw)
	{
		jhx=0;
		jhy+=hstep;
		if(jhy>=im1->ih)
			jhy=0;
	}
#endif
#if 0
	const int rowlen=iw;
	//if(iw>rowlen)
	{
		int px=0, py=ih>>1;
		//int px=rand()%(iw-rowlen), py=rand()%ih;
		static ArrayHandle vertices=0;
		for(int kx=px;kx<px+rowlen;++kx)
		{
			int idx=(iw*py+kx)<<2;
			unsigned char r=image[idx|0], g=image[idx|1], b=image[idx|2];
			draw_3d_line_enqueue(&vertices, cubesize*r/255, cubesize*g/255, height*cubesize*b/255);
		}
		int color=0xFF000000;
		int texture[]={color, color, color, color};
		draw_3d_flush(vertices, &cam, texture, 2, 2, 0, GL_LINE_STRIP);
	}
#endif
	draw_contour3d(&cam, 0, jh_cubesize, 0, jh_cubesize, 0, jh_cubesize, txid_jointhist, 1<<jointhist_nbits, 0.8f);
	//draw_contour3d_rect(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[2], 0.8f);

	if(show_full_image)
	{
		int bounds[4]={0};//{x1, x2, y1, y2}
		jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, im1->iw, im1->ih, bounds);
		display_texture_i(0, w, 0, h, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 0.5, 0);
		
		int crosshaircolor=0xFF000000;
		float ratios[]={(float)w/im1->iw, (float)h/im1->ih};
		float sbounds[4];
		for(int k=0;k<4;++k)
			sbounds[k]=bounds[k]*ratios[k>>1];
		float xmid=(sbounds[0]+sbounds[1])*0.5f;
		float ymid=(sbounds[2]+sbounds[3])*0.5f;
		draw_line(sbounds[0], sbounds[2], sbounds[0], sbounds[3], crosshaircolor);//{x1, y1, x2, y2}
		draw_line(sbounds[1], sbounds[2], sbounds[1], sbounds[3], crosshaircolor);
		draw_line(sbounds[0], sbounds[2], sbounds[1], sbounds[2], crosshaircolor);
		draw_line(sbounds[0], sbounds[3], sbounds[1], sbounds[3], crosshaircolor);
		draw_line(xmid, 0, xmid, (float)h, crosshaircolor);
		draw_line(0, ymid, (float)w, ymid, crosshaircolor);
	}
}


//active keys turn on timer
#define ACTIVE_KEY_LIST\
	AK('W') AK('A') AK('S') AK('D') AK('T') AK('G')\
	AK(KEY_LEFT) AK(KEY_RIGHT) AK(KEY_UP) AK(KEY_DOWN)\
	AK(KEY_ENTER) AK(KEY_BKSP)
int active_keys_pressed=0;

//mouse
char drag=0;
int mx0=0, my0=0;

typedef struct AABBStruct
{
	float x1, x2, y1, y2;
} AABB;
AABB buttons[6]={0};//0: CT,  1: ST,  2: list of transforms,  3: jxl params,  4: EC,  5: OLS-4

int io_init(int argc, char **argv)//return false to abort
{
	(void)argc;//FIXME open 1 command argument
	(void)argv;

	set_window_title("Entropy Benchmark");
	glClearColor(1, 1, 1, 1);

	cam_zoomIn(cam, 1);
	cam_turnMouse(cam, 0, 0, mouse_sensitivity);
	memcpy(&cam0, &cam, sizeof(cam));

	set_bk_color(0xA0808080);
	return 1;
}

const int
//	custom_rct_h=2, custom_rct_w=2,
	custom_pred_reach=2;
const int
	gui_custom_rct_w=29, gui_custom_rct_h=5,//characters
	gui_custom_pred_w=103+2, gui_custom_pred_h=4,
	gui_ec_width=27,
	gui_ols4_elementchars=9;
int custom_pred_ch_idx=0;//from {0, 1, 2}
void io_resize(void)
{
	AABB *p=buttons;
	float xstep=tdx*guizoom, ystep=tdy*guizoom;
	p->x1=xstep*2, p->x2=p->x1+xstep*gui_custom_rct_w, p->y1=(float)(h>>1), p->y2=p->y1+ystep*gui_custom_rct_h, ++p;//0: color params - left
	p->x1=(float)(w>>3), p->x2=p->x1+xstep*gui_custom_pred_w, p->y1=(float)((h>>1)+(h>>2))-ystep, p->y2=p->y1+ystep*gui_custom_pred_h, ++p;//1: spatial params - bottom
	p->x1=(float)(w-300), p->x2=(float)w, p->y1=tdy*2, p->y2=p->y1+tdy*T_COUNT/2, ++p;//2: transforms list
	
	//p->x1=(float)(w>>1), p->x2=p->x1+xstep*14, p->y1=(float)((h>>1)+(h>>2))+ystep*4, p->y2=p->y1+ystep, ++p;//3: clamp bounds
	//p->x1=(float)(w>>1), p->x2=p->x1+xstep*21, p->y1=(float)((h>>1)+(h>>2))+ystep*5, p->y2=p->y1+ystep, ++p;//4: learning rate

	p->x1=(float)(w>>2), p->x2=p->x1+tdx*11*6, p->y1=(float)((h>>1)+(h>>2)), p->y2=p->y1+tdy*3, ++p;//3: jxl params

	p->x1=(float)(w-450), p->x2=p->x1+tdx*gui_ec_width, p->y1=tdy, p->y2=p->y1+tdy, ++p;//4: EC method	//H.E.M.L..A.0x0000..XXXX_XXX

	p->x1=(float)(w>>2), p->x2=p->x1+xstep*(OLS4_RMAX<<1|1)*gui_ols4_elementchars, p->y1=(float)(h>>1)+10, p->y2=p->y1+ystep*(1+(OLS4_RMAX+1)*(im1?im1->nch:4)), ++p;//5: OLS-4		DON'T USE p->y2
}
int io_mousemove(void)//return true to redraw
{
#if 0
	if(
		mode==VIS_IMAGE_BLOCK||
		mode==VIS_BAYES||
		//mode==VIS_IMAGE_E24||
		mode==VIS_DWT_BLOCK
	)
	{
		if(drag)
		{
			show_mouse(drag);
			drag=0;
		}
		if(GET_KEY_STATE(KEY_LBUTTON))
		{
			blockmx=mx, blockmy=my;
			//if(mode==VIS_IMAGE_E24)
			//	e24_update();
			return 1;
		}
	}
	else
#endif
	if(drag)
	{
		if(im1&&mode==VIS_JOINT_HISTOGRAM&&show_full_image)
		{
			jhc_xbox=mx;
			jhc_ybox=my;
			if(GET_KEY_STATE(KEY_SHIFT))
			{
				int step;

				step=jhc_boxdx*w/im1->iw;
				jhc_xbox-=jhc_xbox%step-(step>>1);
				step=jhc_boxdy*h/im1->ih;
				jhc_ybox-=jhc_ybox%step-(step>>1);
			}
			chart_jointhist_update(im1, txid_jointhist);
		}
		else
		{
			int X0=w>>1, Y0=h>>1;
			cam_turnMouse(cam, mx-X0, my-Y0, mouse_sensitivity);
			set_mouse(X0, Y0);
		}
		return !timer;
	}
	return 0;
}
static void click_hittest(int _mx, int _my, int *objidx, int *cellx, int *celly, int *cellidx, AABB **p)
{
	*p=buttons;
	*objidx=0;
	for(*objidx=0;*objidx<_countof(buttons);++*objidx, ++*p)
	{
		//when these buttons are inactive they shouldn't block the click
		if(
			(!transforms_customenabled&&(*objidx==0||*objidx==1))||
			(!(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])&&*objidx==3)||
			(!(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])&&*objidx==5)
		)
			continue;
		if(_mx>=p[0]->x1&&_mx<p[0]->x2&&_my>=p[0]->y1&&_my<p[0]->y2)
			break;
	}
	switch(*objidx)
	{
	case 0://color transform params
		*cellx=(int)floorf((_mx-p[0]->x1)*gui_custom_rct_w/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((_my-p[0]->y1)*gui_custom_rct_h/(p[0]->y2-p[0]->y1));
		*cellidx=0;
		//*cellx=mx-p[0]->x1>(p[0]->x2-p[0]->x1)*0.5f;
		//*celly=(int)floorf((my-p[0]->y1)*custom_rct_h/(p[0]->y2-p[0]->y1));
		//*cellidx=custom_rct_w**celly+*cellx;
		break;
	case 1://spatial transform params
		*cellx=(int)floorf((_mx-p[0]->x1)*(custom_pred_reach<<1|1)/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((_my-p[0]->y1)*gui_custom_pred_h/(p[0]->y2-p[0]->y1));
		*cellidx=(custom_pred_reach<<1|1)*(*celly-1)+*cellx;
		break;
	case 2://list of transforms
		*cellx=(int)floorf((_mx-p[0]->x1)*2/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((_my-p[0]->y1)*(T_COUNT/2)/(p[0]->y2-p[0]->y1));
		*cellidx=*celly;
		break;
	//case 3://clamp bounds
	//	*cellx=mx-p[0]->x1>(p[0]->x2-p[0]->x1)*0.5f;
	//	*celly=0;
	//	*cellidx=*cellx;
	//	break;
	//case 4://learning rate
	//	*cellx=0;
	//	*celly=0;
	//	*cellidx=*cellx;
	//	break;
	case 3://jxl params
		*cellx=(int)floorf((_mx-p[0]->x1)*11/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((_my-p[0]->y1)* 3/(p[0]->y2-p[0]->y1));
		*cellidx=11**celly+*cellx;
		break;
	case 4:
		//H.E.M.L..A.0x0000..XXXX_XXX
		*cellx=(int)floorf((_mx-p[0]->x1)*gui_ec_width/(p[0]->x2-p[0]->x1));
		*celly=0;
		*cellidx=*cellx;
		break;
	case 5:
		*cellx=(int)floorf((_mx-p[0]->x1)/(tdx*guizoom));
		*celly=(int)floorf((_my-p[0]->y1)/(tdy*guizoom));
		*cellidx=0;
		break;
	default:
		*objidx=-1;
		*p=0;
		*cellx=-1;
		*celly=-1;
		*cellidx=-1;
		break;
	}
}
int io_mousewheel(int forward)
{
	//if(im1&&(transforms_customenabled||transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP]))//change custom transform params
	if(im1)
	{
		int objidx=0, cellx=0, celly=0, cellidx=0;
		AABB *p=buttons;
		click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &p);
		if(objidx!=-1)
		{
			int sign=(forward>0)-(forward<0);//abs(forward) is 120
			int ch=(int)floorf((mx-p->x1)/(guizoom*tdx));
			switch(objidx)
			{
			case 0://color transform params
				if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
				{
					//P0  rgb
					//r += (-0x0000*g-0x0000*b)>>12
					//g += (-0x0000*r-0x0000*b)>>12
					//b += (-0x0000*r-0x0000*g)>>12
					//g += (-0x0000*r-0x0000*b)>>12
					//012345678901234567890123456789
					if(!celly)
					{
						int x=rct_custom_params[8];
						x+=sign;
						MODVAR(x, x, 6);
						rct_custom_params[8]=(short)x;
					}
					else if((unsigned)(cellx-9)<4||(unsigned)(cellx-18)<4)
					{
						int idx;
						if((unsigned)(cellx-9)<4)
							idx=(celly-1)<<1|0, cellx-=9;
						else
							idx=(celly-1)<<1|1, cellx-=18;
						rct_custom_params[idx]+=(short)(sign<<((3-cellx)<<2));
					}
#if 0
					//0000000000111111111122222222223333
					//0123456789012345678901234567890123
					//r-=g
					//g+=(-0x00.0000*r-0x00.0000*b)>>16
					//b-=g
					//g+=(-0x00.0000*r-0x00.0000*b)>>16
					int
						//col=p?(int)floorf((mx-p->x1)/(guizoom*tdx)):0,
						line=p?(int)floorf((my-p->y1)/(guizoom*tdy)):0;
					if((line&1)==1&&((unsigned)(ch-7)<7||(unsigned)(ch-19)<7))
					{
						int idx=(line&2)|((ch-7)/(19-7));
						int digit=(ch-7)%(19-7);
						//0123456
						//00.0000
						if(digit==2)
							break;
						digit-=digit>2;
						digit=1-digit;
						rct_custom_params[idx]+=sign<<((digit+4)<<2);//hex digit
					}
#endif
					update_image();
				}
#if 0
				//000000000011111111112222222222
				//012345678901234567890123456789
				//r-=(>>nnnN.NNN)g+(  nnnN.NNN)b
				ch2=ch-6;
				MODVAR(ch2, ch2, 14);
				if(ch2>=0&&ch2<8)
				{
					ch2-=4;
					ch2+=ch2<0;//skip point
					double delta=pow(10, -ch2);
					customparam_ct[cellidx]+=sign*delta;
				}
				else
					customparam_ct[cellidx]+=sign*0.05;
#endif
				break;
			case 1://spatial transform params
				if(!transforms_mask[ST_FWD_CUSTOM]&&!transforms_mask[ST_INV_CUSTOM])
					break;
				if(!celly)
				{
					//0123456789012345678901234567
					//Ch 0  Clamp [+W +NW +N +NE]
					if(ch==3)
					{
						custom_pred_ch_idx+=sign;
						MODVAR(custom_pred_ch_idx, custom_pred_ch_idx, 3);
					}
					else if(BETWEEN_EXC(13, ch, 16))//W
						custom_clamp[0]=!custom_clamp[0];
					else if(BETWEEN_EXC(16, ch, 20))//NW
						custom_clamp[1]=!custom_clamp[1];
					else if(BETWEEN_EXC(20, ch, 23))//N
						custom_clamp[2]=!custom_clamp[2];
					else if(BETWEEN_EXC(23, ch, 26))//NE
						custom_clamp[3]=!custom_clamp[3];
				}
				else
				{
					//0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
					//0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
					//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					int col=ch%21, idx=5*(celly-1)+(ch/21);
					if(idx<12&&((unsigned)(col-3)<7||(unsigned)(col-13)<7))
					{
						idx=idx<<1|(col>=10);
						//0123456
						//00.0000
						int digit=col%10-3;
						if(digit==2)
							break;
						digit-=digit>2;
						digit=1-digit;
						custom_params[24*custom_pred_ch_idx+idx]+=sign<<((digit+4)<<2);
					}
				}
				update_image();
#if 0
				//000000000011111111112222222222333333333344444444445555
				//012345678901234567890123456789012345678901234567890123
				//>>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN
				ch2=ch-2;
				MODVAR(ch2, ch2, 11);
				if(ch>=0&&ch<54&&ch2>=3&&ch2<8)
				{
					ch2-=4;
					ch2+=ch2<0;//skip point
					double delta=pow(10, -ch2);
					customparam_st[12*customparam_ch_idx+cellidx]+=sign*delta;
				}
				else
					customparam_st[12*customparam_ch_idx+cellidx]+=sign*0.05;
#endif
				break;
#if 0
			case 3://clamp bounds
			
				//012345678901234567890		ch
				//   1234   1234			ch2
				//[ SNNNN, SNNNN] clamp
				//   3210   3210			power
				ch2=ch-2;
				MODVAR(ch2, ch2, 7);
				if(ch>=3&&ch<14&&ch2>=1&&ch2<5)
				{
					ch2=4-ch2;
					int delta=(int)pow(10, ch2);
					customparam_clamp[cellidx]+=sign*delta;
					if(cellidx)
					{
						if(sign<0)//upper decreased & fell below lower
						{
							if(customparam_clamp[0]>customparam_clamp[1])
								customparam_clamp[0]=customparam_clamp[1];
						}
					}
					else
					{
						if(sign>0)//lower increased & rose above upper
						{
							if(customparam_clamp[1]<customparam_clamp[0])
								customparam_clamp[1]=customparam_clamp[0];
						}
					}
				}
				break;
			case 4://learning rate

				//0123456789012345678901
				//lr -0.NNNNNNNNNNNNNNN
				//     -123456789123456
				if(ch>=6&&ch<21)
				{
					ch=5-ch;
					double delta=pow(10, ch);
					g_lr+=sign*delta;
				}
				break;
#endif
			case 3://jxl params
				if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					int ch2;
					ch=(int)floorf((mx-p->x1)/tdx);
					MODVAR(ch2, ch, 6);
					if(ch2>=2&&ch2<6)
					{
						int delta=sign<<((5-ch2)<<2);
						jxlparams_i16[cellidx]+=(short)delta;
					}
					else if(cellx>=4)
						jxlparams_i16[cellidx]=-jxlparams_i16[cellidx];
					update_image();
				}
				break;
			case 4://EC method
				switch(cellx)
				{
				//0123456789012345678901234567
				//H.E.M.L..A.0x0000..XXXX_XXX
				case 2:
					ec_expbits+=sign;
					MODVAR(ec_expbits, ec_expbits, 10);
					break;
				case 4:
					ec_msb+=sign;
					MODVAR(ec_msb, ec_msb, 10);
					break;
				case 6:
					ec_lsb+=sign;
					MODVAR(ec_lsb, ec_lsb, 10);
					break;
				case 9:
					ec_adaptive=!ec_adaptive;
					//ec_adaptive+=sign;
					//MODVAR(ec_adaptive, ec_adaptive, 3);
					break;
				case 13:
					if(ec_adaptive)
					{
						ec_adaptive_threshold+=sign<<12;
						ec_adaptive_threshold&=0xFFFF;
					}
					break;
				case 14:
					if(ec_adaptive)
					{
						ec_adaptive_threshold+=sign<<8;
						ec_adaptive_threshold&=0xFFFF;
					}
					break;
				case 15:
					if(ec_adaptive)
					{
						ec_adaptive_threshold+=sign<<4;
						ec_adaptive_threshold&=0xFFFF;
					}
					break;
				case 16:
					if(ec_adaptive)
					{
						ec_adaptive_threshold+=sign;
						ec_adaptive_threshold&=0xFFFF;
					}
					break;
				case 19:case 20:case 21:case 22:case 23:case 24:case 25:case 26:
					ec_method+=sign;
					MODVAR(ec_method, ec_method, ECTX_COUNT);
					break;
				}
				if(ec_expbits<ec_msb+ec_lsb)
					ec_msb=ec_lsb=0;
				update_image();
				break;
			case 5://OLS-4
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					if(celly)
					{
						--celly;
						int kx=cellx/gui_ols4_elementchars, kchar=7-cellx%gui_ols4_elementchars, kc=celly/(OLS4_RMAX+1), ky=celly%(OLS4_RMAX+1);
						int idx=(OLS4_RMAX<<1|1)*ky+kx;
						if(idx<OLS4_CTXSIZE+1&&kchar>=0)
						{
							ols4_mask[kc][idx]^=1<<kchar;
							if(idx==OLS4_CTXSIZE)//causality mask
							{
								int cmask=((1<<(kc<<1))-1);
								ols4_mask[kc][idx]&=cmask<<4|cmask;
							}
							if(!GET_KEY_STATE(KEY_CTRL))
								update_image();
						}
					}
					else
					{
						int kx=cellx/gui_ols4_elementchars;
						switch(kx)
						{
						case 0:
							ols4_period=SHIFT_LEFT_SIGNED(ols4_period, sign);
							ols4_period=CLAMP(1, ols4_period, 8192);
							if(!GET_KEY_STATE(KEY_CTRL))
								update_image();
							break;
						case 1:
						case 2:
						case 3:
						case 4:
							{
								//	0	X	2	3	4	...
								//	0	.	-1	-2	-3	...
								int digit=-(cellx%gui_ols4_elementchars);
								if(digit!=-1)
								{
									digit+=digit<=-1;
									double val=ols4_lr[kx-1];
									val+=sign*_10pow(digit);
									ols4_lr[kx-1]=CLAMP(0, val, 1);
									if(!GET_KEY_STATE(KEY_CTRL))
										update_image();
								}
							}
							break;
						}
					}
				}
				break;
			}//switch
		}//if
		else
			goto normal_operation;
	}
	else
	{
	normal_operation:
#if 0
		if(
			mode==VIS_IMAGE_BLOCK||
			mode==VIS_BAYES||
		//	mode==VIS_IMAGE_E24||
			mode==VIS_DWT_BLOCK
		)
		{
			if(GET_KEY_STATE(KEY_CTRL))
			{
				if(forward>0)
				{
					margin<<=1;
					if(margin>1024)
						margin=1024;
				}
				else
				{
					margin>>=1;
					if(margin<1)
						margin=1;
				}
			}
			else
			{
				if(forward>0)
				{
					blocksize<<=1;
					if(blocksize>1024)
						blocksize=1024;
				}
				else
				{
					blocksize>>=1;
					if(blocksize<1)
						blocksize=1;
				}
				//e24_update();
			}
		}
		else
#endif
		if(mode==VIS_JOINT_HISTOGRAM&&show_full_image)
		{
			if(GET_KEY_STATE(KEY_CTRL))
			{
				if(forward>0)
				{
					jhc_boxdx<<=1;
					jhc_boxdy<<=1;
				}
				else
				{
					jhc_boxdx>>=1;
					jhc_boxdy>>=1;
				}
				jhc_boxdx=CLAMP(0, jhc_boxdx, 256);
				jhc_boxdy=CLAMP(0, jhc_boxdy, 256);
			}
			else
			{
				jhc_level+=THREEWAY(forward, 0);
				UPDATE_MAX(jhc_level, 0.5f);
			}
			chart_jointhist_update(im1, txid_jointhist);
		}
		else if(GET_KEY_STATE(KEY_SHIFT))//shift wheel		change cam speed
		{
			if(forward>0)	cam.move_speed*=2;
			else		cam.move_speed*=0.5f;
		}
		else
		{
			if(forward>0)	cam_zoomIn(cam, 1.1f);
			else		cam_zoomOut(cam, 1.1f);
		}
	}
	return !timer;
}
static void count_active_keys(IOKey upkey)
{
	keyboard[upkey]=0;
	active_keys_pressed=0;
#define		AK(KEY)		active_keys_pressed+=keyboard[KEY];
	ACTIVE_KEY_LIST
#undef		AK
	if(!active_keys_pressed)
		timer_stop(TIMER_ID_KEYBOARD);
}
static int parse_nvals_i8(ArrayHandle text, int idx, unsigned char *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		int neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		char *end=(char*)text->data+idx;
		params[k]=(unsigned char)strtol((char*)text->data+idx, &end, 16);
		idx=(int)(end-(char*)text->data);
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(text->data[idx])&&text->data[idx]!='-'&&!isdigit(text->data[idx]);++idx);//skip comma

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static int parse_nvals_f64(ArrayHandle text, int idx, double *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		int neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		//if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
		//	idx+=2;
		char *end=(char*)text->data+idx;
		params[k]=strtod((char*)text->data+idx, &end);
		idx=(int)(end-(char*)text->data);
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(text->data[idx])&&text->data[idx]!='-'&&!isdigit(text->data[idx]);++idx);//skip comma

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static int parse_nvals_i32(ArrayHandle text, int idx, int *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		int neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		char *end=(char*)text->data+idx;
		params[k]=(int)strtol((char*)text->data+idx, &end, 16);
		idx=(int)(end-(char*)text->data);
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(text->data[idx])&&text->data[idx]!='-'&&!isdigit(text->data[idx]);++idx);//skip comma

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static int parse_nvals(ArrayHandle text, int idx, short *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		int neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		char *end=(char*)text->data+idx;
		params[k]=(short)strtol((char*)text->data+idx, &end, 16);
		idx=(int)(end-(char*)text->data);
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(text->data[idx])&&text->data[idx]!='-'&&!isdigit(text->data[idx]);++idx);//skip comma

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static void append_i16_row(ArrayHandle *str, const short *vals, int count)
{
	for(int k=0;k<count;++k)
	{
		short val=vals[k];
		str_append(str, "%c0x%04X,", val<0?'-':' ', abs(val));
	}
	str_append(str, "\n");
}
int io_keydn(IOKey key, char c)
{
	(void)c;

	switch(key)//handle shortcuts involving timed keys before the main switch
	{
	case 'S':
		if(keyboard[KEY_CTRL]&&im1&&im_export)
		{
			char *filename=dialog_save_file(0, 0, "Untitled.PNG");
			if(filename)
			{
				lodepng_encode_file(filename, im_export, im1->iw, im1->ih, LCT_RGBA, 8);
				free(filename);
			}
			return 1;
		}
		break;
	default://to make gcc -Wall happy
		break;
	}
	switch(key)
	{
	case KEY_LBUTTON:
	case KEY_RBUTTON:
		if(im1)
		{
			int objidx=-1, cellx=0, celly=0, cellidx=0;
			AABB *p=buttons;
			click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &p);
			switch(objidx)
			{
			case 0://color transform params
				if(key==KEY_RBUTTON)
				{
					//P0  rgb
					//r += (-0x0000*g-0x0000*b)>>12
					//g += (-0x0000*r-0x0000*b)>>12
					//b += (-0x0000*r-0x0000*g)>>12
					//g += (-0x0000*r-0x0000*b)>>12
					//012345678901234567890123456789
					if(GET_KEY_STATE(KEY_CTRL))
					{
						memset(rct_custom_params, 0, sizeof(rct_custom_params));
						update_image();
						return 1;
					}
					else if(!celly)
					{
						rct_custom_params[8]=0;
						update_image();
						return 1;
					}
					else if((unsigned)(cellx-9)<4||(unsigned)(cellx-18)<4)
					{
						int idx;
						if((unsigned)(cellx-9)<4)
							idx=(celly-1)<<1|0, cellx-=9;
						else
							idx=(celly-1)<<1|1, cellx-=18;
						rct_custom_params[idx]=0;
						update_image();
						return 1;
					}
				}
#if 0
				{
					int
						col=p?(int)floorf((mx-p->x1)/(guizoom*tdx)):0,
						line=p?(int)floorf((my-(p->y1+guizoom*tdy))/(guizoom*tdy)):0;
					if(key==KEY_RBUTTON&&(line==1||line==3)&&((unsigned)(col-7)<7||(unsigned)(col-19)<7))
					{
						//0000000000111111111122222222223333
						//0123456789012345678901234567890123
						//r-=g
						//g+=(-0x00.0000*r-0x00.0000*b)>>16
						//b-=g
						//g+=(-0x00.0000*r-0x00.0000*b)>>16
						int idx=(line>2)<<1|(col>17);
						rct_custom_params[idx]=0;
						update_image();
						return 1;
					}
				}
#endif
				break;
			case 1://spatial transforms params
				{
					int
						col=p?(int)floorf((mx-p->x1)/(guizoom*tdx)):0,
						line=p?(int)floorf((my-(p->y1+guizoom*tdy))/(guizoom*tdy)):0;
					if(key==KEY_RBUTTON&&cellidx<12)
					{
						//0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
						//0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
						//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						int idx=5*line+col/21;
						col%=21;
						idx=idx<<1|(col>=10);
						custom_params[24*custom_pred_ch_idx+idx]=0;
						update_image();
						return 1;
					}
				}
				break;
			case 2://transform list
				{
					int
						//col=p?(int)floorf((mx-p->x1)/tdx):0,
						line=p?(int)floorf((my-p->y1)/tdy):0;
					//if(BETWEEN_EXC(0, cellidx, T_COUNT/2))
					if((unsigned)line<(unsigned)(T_COUNT/2))
					{
						int idx=line<<1|(p?(int)floor((mx-p->x1)*2/(p->x2-p->x1)):0);
						if(key==KEY_LBUTTON)
							transforms_append(idx);
						else
							transforms_removebyid(idx);
						update_image();
						return 1;
					}
				}
				break;
			//case 3://clamp bounds
			//	if(key==KEY_RBUTTON)
			//	{
			//		if(cellidx)
			//			customparam_clamp[1]=127;
			//		else
			//			customparam_clamp[0]=-128;
			//		update_image();
			//		return 1;
			//	}
			//	break;
			//case 4://spatial transforms params
			//	if(key==KEY_RBUTTON)
			//	{
			//		g_lr=1e-10;
			//		return 1;
			//	}
			//	break;
			case 3://jxl params
				if(key==KEY_RBUTTON)
				{
					if(cellx>=4)
						jxlparams_i16[cellidx]=0;
					else
						jxlparams_i16[cellidx]=4096;
					update_image();
					return 1;
				}
				break;
			case 4://EC method
				if(key==KEY_LBUTTON)
					break;
				switch(cellx)
				{
				//012345678901234567
				//H.E.M.L..A.0x0000..XXXX_XXX
				case 0:
				case 1:
				case 2:
				case 3:
				case 4:
				case 5:
				case 6:
					ec_expbits=5;
					ec_msb=2;
					ec_lsb=0;
					break;
				case 9:
					ec_adaptive=0;
					break;
				case 13:
				case 14:
				case 15:
				case 16:
					if(ec_adaptive)
						ec_adaptive_threshold=3200;
					break;
				case 19:case 20:case 21:case 22:case 23:case 24:case 25:case 26:
					ec_method=ECTX_HIST;
					//ec_method=ECTX_MIN_QN_QW;
					break;
				}
				if(ec_expbits<ec_msb+ec_lsb)
					ec_msb=ec_lsb=0;
				update_image();
				return 1;
			case 5://OLS-4
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					if(celly)
					{
						--celly;
						int kx=cellx/gui_ols4_elementchars, kchar=7-cellx%gui_ols4_elementchars, kc=celly/(OLS4_RMAX+1), ky=celly%(OLS4_RMAX+1);
						int idx=(OLS4_RMAX<<1|1)*ky+kx;
						if(idx<OLS4_CTXSIZE+1&&kchar>=0)
						{
							if(key==KEY_LBUTTON)
							{
								ols4_mask[kc][idx]=ols4_cache;
								if(!GET_KEY_STATE(KEY_CTRL))
									update_image();
							}
							else
								ols4_cache=ols4_mask[kc][idx];
							return 1;
						}
					}
					else
					{
						int kx=cellx/gui_ols4_elementchars;
						switch(kx)
						{
						case 0:
							ols4_period=OLS4_DEFAULT_PERIOD;
							if(!GET_KEY_STATE(KEY_CTRL))
								update_image();
							return 1;
						case 1:
						case 2:
						case 3:
						case 4:
							ols4_lr[kx-1]=OLS4_DEFAULT_LR;
							if(!GET_KEY_STATE(KEY_CTRL))
								update_image();
							return 1;
						}
					}
				}
				break;
			default:
				if(key==KEY_LBUTTON)
					goto toggle_drag;
				break;
			}
		}
		else if(key==KEY_LBUTTON)
			goto toggle_drag;
		break;
	case KEY_ESC:
	toggle_drag:
#if 0
		if(
			mode==VIS_IMAGE_BLOCK||
			mode==VIS_BAYES||
		//	mode==VIS_IMAGE_E24||
			mode==VIS_DWT_BLOCK
		)
		{
			if(key==KEY_LBUTTON)
			{
				blockmx=mx, blockmy=my;
				//e24_update();
				return 1;
			}
		}
		else
#endif
		{
			show_mouse(drag);
			drag=!drag;
			if(mode==VIS_JOINT_HISTOGRAM&&show_full_image)
			{
				jhc_xbox=mx;
				jhc_ybox=my;
				chart_jointhist_update(im1, txid_jointhist);
				return 1;
			}
			if(drag)//enter mouse control
			{
				mx0=mx, my0=my;
				set_mouse(w>>1, h>>1);
			}
			else//leave mouse control
				set_mouse(mx0, my0);
		}
		break;
	case KEY_MBUTTON:
		//printf("Click at (%d, %d)\n", mx, my);
		break;

#define AK(KEY) case KEY:
	ACTIVE_KEY_LIST
#undef  AK
		timer_start(10, TIMER_ID_KEYBOARD);
		break;
		
	case KEY_F1:
		messagebox(MBOX_OK, "Controls",
			"Ctrl O:\t\tOpen image\n"
			"Ctrl S:\t\tSave preview as 8-bit PNG\n"
			"Mouse1/Mouse2:\tAdd/remove transforms in the list\n"
			"Ctrl Mouse1:\tReplace all transforms of this type\n"
			"Ctrl R:\t\tDisable all transforms\n"
			//"Ctrl E:\t\tReset custom transform parameters\n"
			"Comma/Period:\tSelect context for size estimation\n"
			"Slash:\t\tToggle adaptive histogram\n"
			"[ ]:\t\t(Custom transforms) Select coefficient page\n"
			"Ctrl 1/2/3:\t(Custom predictor) Populate page from channel N\n"
			"Space:\t\t(Custom transforms) Optimize\n"
			"Ctrl N:\t\tAdd noise to CUSTOM3 params\n"
			//"Shift space:\t(Custom transforms) Optimize blockwise\n"
			//"Ctrl Space\t(Custom transforms) Reset params\n"
			"Ctrl B:\t\tBatch test\n"
			"\n"
			"WASDTG:\tMove cam\n"
			"Arrow keys:\tTurn cam\n"
			"Mouse1/Esc:\tToggle mouse look\n"
			"R:\t\tReset cam\n"
			"\n"
			"H:\t\tReset CR history graph\n"
			"Ctrl C:\t\tCopy data\n"
			"Ctrl V:\t\tPaste data\n"
			"Ctrl B:\t\tBatch test\n"
			"Ctrl SPACE:\tCheck integrity (if restored bit-exact)\n"
		//	"Ctrl P:\t\tTest predictors\n"
			"C:\t\tToggle joint histogram type / fill screen in image view\n"
			"\n"
			"M / Shift M:\tCycles between:\n"
		//	"\t1: 3D View: Levels\n"
		//	"\t2: 3D View: Mesh\n"
		//	"\t3: 3D View: Mesh (separate channels)\n"
			"\t%d: Image tricolor view\n"
			"\t%d: Image view\n"
		//	"\t6: Image block histogram\n"
		//	"\t7: Optimized block compression estimate (E24)\n"
		//	"\t7: DWT block histogram\n"
			"\t%d: Histogram\n"
			"\t%d: Joint histogram\n"
			"\t%d: Zipf view\n",

			1+VIS_IMAGE_TRICOLOR,
			1+VIS_IMAGE,
			1+VIS_HISTOGRAM,
			1+VIS_JOINT_HISTOGRAM,
			1+VIS_ZIPF
		);
		//prof_on=!prof_on;
		return 0;
	case 'R':
		if(GET_KEY_STATE(KEY_CTRL)&&im1)
		{
			transforms_removeall();
			update_image();
		}
		else
			memcpy(&cam, &cam0, sizeof(cam));
		return 1;
	case KEY_COMMA:
	case KEY_PERIOD:
		{
			int fwd=key==KEY_PERIOD;
			fwd-=!fwd;
			ec_method+=fwd;
			MODVAR(ec_method, ec_method, ECTX_COUNT);
			update_image();
		}
		return 1;
	case KEY_SLASH:
		ec_adaptive=!ec_adaptive;
		//ec_adaptive=(ec_adaptive+1)%3;
		update_image();
		return 1;
	//case 'E':
	//	if(im1&&GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)//reset params
	//	{
	//		customtransforms_resetparams();
	//		update_image();
	//		return 1;
	//	}
	//	wireframe=!wireframe;
	//	break;
	case 'O':
		if(GET_KEY_STATE(KEY_CTRL))
		{
			ArrayHandle fn2=dialog_open_file(0, 0, 0);
			if(fn2)
			{
				Image *im2=image_load((char*)fn2->data, (int)fn2->count);
				if(im2)
				{
					if(fn)
						array_free(&fn);
					fn=fn2;
					filesize=get_filesize((char*)fn2->data);

					if(im0)
						free(im0);
					im0=im2;
					update_image();
					set_window_title("%s - eBench", (char*)fn->data);
				}
				else
					array_free(&fn2);
			}
		}
		return 1;
	case 'C':
		if(im1&&GET_KEY_STATE(KEY_CTRL))//copy custom transform value
		{
			if(GET_KEY_STATE(KEY_SHIFT))
			{
				unsigned char *buf=0;
				image_export_uint8(im1, &buf, 1, 1);//swap red & blue for WinAPI
				if(!buf)
				{
					LOG_WARNING("Alloc error");
					return 0;
				}
				int success=copy_bmp_to_clipboard(buf, im1->iw, im1->ih);
				if(!success)
					LOG_WARNING("Failed to copy image to clipboard");
				free(buf);
				return 0;
			}
			ArrayHandle str;
			STR_ALLOC(str, 0);
			if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
			{
				str_append(&str, "%8d %8.6lf %8.6lf %8.6lf %8.6lf\n", ols4_period, ols4_lr[0], ols4_lr[1], ols4_lr[2], ols4_lr[3]);
				for(int kc=0;kc<4;++kc)
				{
					for(int ky=0;ky<OLS4_RMAX+1;++ky)
					{
						for(int kx=0;kx<(OLS4_RMAX<<1|1);++kx)
						{
							str_append(&str, "0x%02X,", ols4_mask[kc][(OLS4_RMAX<<1|1)*ky+kx]);
							if(ky==OLS4_RMAX&&kx==OLS4_RMAX)
								goto next_channel;
						}
						str_append(&str, "\n");
					}
				next_channel:
					str_append(&str, "\n");
				}
			}
			if(mode==VIS_IMAGE||mode==VIS_ZIPF)
			{
				double cr_combined=(im1->src_depth[0]+im1->src_depth[1]+im1->src_depth[2]+im1->src_depth[3])/(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3]);
				str_append(&str, "TRGBA %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf",
					100./cr_combined,
					100.*ch_entropy[0]/im1->src_depth[0],
					100.*ch_entropy[1]/im1->src_depth[1],
					100.*ch_entropy[2]/im1->src_depth[2],
					100.*ch_entropy[3]/im1->src_depth[1]
				);
			}
#if 0
			//else if(mode==VIS_IMAGE_E24)
			//{
			//	for(int kc=0;kc<3;++kc)
			//	{
			//		E24Params const *p=e24_params+kc;
			//		str_append(&str, "%3d  %3d %3d %3d  %3d %3d\n", p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
			//	}
			//}
			else if(transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC])
			{
				for(int kc=0, idx=0;kc<3;++kc)
				{
					for(int ky=0;ky<LOGIC_NF0;++ky)
					{
						for(int kx=0;kx<LOGIC_ROWPARAMS;++kx, ++idx)
						{
							short val=logic_params[idx];
							str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
						}
						str_append(&str, "\n");
					}
#ifdef LOGIC_NF1
					for(int kx=0;kx<LOGIC_NF1;++kx, ++idx)
					{
						short val=logic_params[idx];
						str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
					}
					str_append(&str, "\n");
#else
					for(int kx=0;kx<LOGIC_NF0;++kx, ++idx)
					{
						short val=logic_params[idx];
						str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
					}
					str_append(&str, "\n");
#endif
				}
			}
			else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
			{
				for(int k=0;k<6;++k)
				{
					append_i16_row(&str, c2_params.c0+k*12, 5);
					append_i16_row(&str, c2_params.c0+k*12+5, 5);
					append_i16_row(&str, c2_params.c0+k*12+10, 2);
				}
				str_append(&str, "\n");

				for(int k=0;k<6;++k)
				{
					append_i16_row(&str, c2_params.c1+k*12, 5);
					append_i16_row(&str, c2_params.c1+k*12+5, 5);
					append_i16_row(&str, c2_params.c1+k*12+10, 2);
				}
				append_i16_row(&str, c2_params.c1+6*12, 2);
				str_append(&str, "\n");

				for(int k=0;k<6;++k)
				{
					append_i16_row(&str, c2_params.c2+k*12, 5);
					append_i16_row(&str, c2_params.c2+k*12+5, 5);
					append_i16_row(&str, c2_params.c2+k*12+10, 2);
				}
				append_i16_row(&str, c2_params.c2+6*12, 4);
			}
#endif
			else if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				const int width=C3_REACH<<2|2;//interleaved pixels & errors
				const short *coeffs[]=
				{
					c3_params.c00, c3_params.c01, c3_params.c02,
					c3_params.c10, c3_params.c11, c3_params.c12,
					c3_params.c20, c3_params.c21, c3_params.c22,
				};
				const int tails[]=
				{
					C3_REACH<<1,      C3_REACH<<1,    C3_REACH<<1,
					(C3_REACH<<1)+2,  C3_REACH<<1,    C3_REACH<<1,
					(C3_REACH<<1)+2, (C3_REACH<<1)+2, C3_REACH<<1,
				};
				for(int kc=0;kc<9;++kc)
				{
					for(int k=0;k<C3_REACH;++k)
						append_i16_row(&str, coeffs[kc]+width*k, width);
					append_i16_row(&str, coeffs[kc]+width*C3_REACH, tails[kc]);
					if(kc+1<9)
						str_append(&str, "\n");
				}
			}
			else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
			{
				for(int ky=0;ky<3;++ky)
				{
					for(int kx=0;kx<11;++kx)
					{
						short val=jxlparams_i16[11*ky+kx];
						str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
						if(kx+1==11)
							str_append(&str, "\n");
					}
				}
			}
			else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
			{
				for(int ky=0;ky<3;++ky)
				{
					for(int kx=0;kx<PW2_NPARAM;++kx)
					{
						short val=pw2_params[PW2_NPARAM*ky+kx];
						str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
						if(kx+1==PW2_NPARAM)
							str_append(&str, "\n");
					}
				}
			}
#if 0
			else if(transforms_mask[ST_FWD_JOINT]||transforms_mask[ST_INV_JOINT])
			{
				for(int ky=0;ky<3;++ky)
				{
					for(int kx=0;kx<24;++kx)
					{
						int val=jointpredparams[24*ky+kx];
						str_append(&str, "%c0x%08X,%c", val<0?'-':' ', abs(val), kx+1<24?' ':'\n');
					}
				}
				for(int ky=0;ky<3;++ky)
				{
					for(int kx=0;kx<4;++kx)
					{
						int val=jointpredparams[72+24*ky+kx];
						str_append(&str, "%c0x%08X,%c", val<0?'-':' ', abs(val), kx+1<4?' ':'\n');
					}
				}
			}
			else if(transforms_mask[ST_FWD_HYBRID3]||transforms_mask[ST_INV_HYBRID3])
			{
				for(int ko=0;ko<3;++ko)
				{
					for(int ky=0;ky<3;++ky)
					{
						for(int kx=0;kx<24;++kx)
							printed+=snprintf((char*)str->data+printed, str->count-printed, "%g%c", customparam_hybrid[24*(3*ko+ky)+kx], kx+1<24?'\t':'\n');
					}
					printed+=snprintf((char*)str->data+printed, str->count-printed, "\n");
				}
			}
#endif
			else if(transforms_customenabled)
			{
				int shift=GET_KEY_STATE(KEY_SHIFT);
				for(int k=0;k<RCT_CUSTOM_NPARAMS;++k)
				{
					int val=rct_custom_params[k];
					str_append(&str, " %c0x%04X", val<0?'-':' ', abs(val));
				}
				str_append(&str, "\n");
				//for(int ky=0;ky<custom_rct_h;++ky)
				//{
				//	for(int kx=0;kx<custom_rct_w;++kx)
				//	{
				//		int val=rct_custom_params[custom_rct_w*ky+kx];
				//		str_append(&str, "\t%c0x%06X", val<0?'-':' ', abs(val));
				//	}
				//	str_append(&str, "\n");
				//}
				const int stw=(custom_pred_reach<<1|1)*2;
				for(int kc2=0;kc2<3;++kc2)
				{
					const int np=_countof(custom_params)/3;
					const int *params=custom_params+24*kc2;
					for(int k=0;k<np;++k)
					{
						int x=k%stw;
						//int y=k/stw;
						if(shift)
							str_append(&str, "%g,%c", (double)params[k]/65536, x<stw-1&&k<np-1?'\t':'\n');
						else
						{
							int val=params[k];
							str_append(&str, "%c0x%06X,%c", val<0?'-':' ', abs(val), x<stw-1&&k<np-1?'\t':'\n');
						}
					}
				}
				//for(int k=0;k<_countof(customparam_clamp);++k)
				//	str_append(&str, "%d%c", customparam_clamp[k], k<_countof(customparam_clamp)-1?'\t':'\n');
			}
			copy_to_clipboard((char*)str->data, (int)str->count);
			array_free(&str);
		}
		else if(mode==VIS_IMAGE||mode==VIS_ZIPF||mode==VIS_IMAGE_TRICOLOR||mode==VIS_JOINT_HISTOGRAM)
		{
			if(mode==VIS_JOINT_HISTOGRAM&&GET_KEY_STATE(KEY_SHIFT))
			{
				int shift=GET_KEY_STATE(KEY_SHIFT);
				space_not_color+=1-(shift<<1);
				MODVAR(space_not_color, space_not_color, 4);
				update_image();
			}
			else
				show_full_image=!show_full_image;

		}
		return 1;
	case 'V':
		if(im1&&GET_KEY_STATE(KEY_CTRL))//paste custom transform value
		{
			ArrayHandle text=paste_from_clipboard(0);
			if(text)
			{
				int k, kend, idx;

				idx=0;
#if 0
				if(mode==VIS_IMAGE_E24)
				{
					k=0, kend=sizeof(e24_params)/sizeof(e24_params->alpha);
					for(;k<kend;++k)
					{
						for(;idx<text->count&&isspace(text->data[idx]);++idx);

						int neg=text->data[idx]=='-';
						idx+=neg;//skip sign
						if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
							idx+=2;
						char *end=text->data+idx;
						((short*)e24_params)[k]=(int)strtol(text->data+idx, &end, 10);
						idx=(int)(end-text->data);
						if(neg)
							((short*)e24_params)[k]=-((short*)e24_params)[k];

						for(;idx<text->count&&!isspace(text->data[idx]);++idx);//skip comma

						if(idx>=text->count)
							goto paste_finish;
					}
				}
				else if(transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC])
				{
					parse_nvals(text, idx, logic_params, _countof(logic_params));
				}
				else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
				{
					parse_nvals(text, idx, (short*)&c2_params, sizeof(c2_params)/sizeof(short));
				}
				else
#endif
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					idx=parse_nvals_i32(text, idx, &ols4_period, 1);
					idx=parse_nvals_f64(text, idx, ols4_lr, _countof(ols4_lr));
					idx=parse_nvals_i8(text, idx, (unsigned char*)ols4_mask, sizeof(ols4_mask)/sizeof(char));
				}
				if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
				{
					parse_nvals(text, idx, (short*)&c3_params, sizeof(c3_params)/sizeof(short));
				}
				else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					parse_nvals(text, idx, jxlparams_i16, _countof(jxlparams_i16));
				}
				else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
				{
					parse_nvals(text, idx, pw2_params, _countof(pw2_params));
				}
#if 0
				else if(transforms_mask[ST_FWD_JOINT]||transforms_mask[ST_INV_JOINT])
				{
					k=0, kend=_countof(jointpredparams);
					for(;k<kend;++k)
					{
						for(;idx<text->count&&isspace(text->data[idx]);++idx);

						int neg=text->data[idx]=='-';
						idx+=neg;//skip sign
						if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
							idx+=2;
						char *end=text->data+idx;
						jointpredparams[k]=(int)strtol(text->data+idx, &end, 16);
						idx=(int)(end-text->data);
						if(neg)
							jointpredparams[k]=-jointpredparams[k];

						for(;idx<text->count&&!isspace(text->data[idx]);++idx);//skip comma

						if(idx>=text->count)
							goto paste_finish;
					}
				}
				else if(transforms_mask[ST_FWD_HYBRID3]||transforms_mask[ST_INV_HYBRID3])
				{
					k=0, kend=_countof(customparam_hybrid), idx=0;
					for(;k<kend;++k)
					{
						for(;idx<text->count&&isspace(text->data[idx]);++idx);
						customparam_hybrid[k]=atof(text->data+idx);
						for(;idx<text->count&&!isspace(text->data[idx]);++idx);
						if(idx>=text->count)
							goto paste_finish;
					}
				}
#endif
				else if(transforms_customenabled)
				{
					int shift=GET_KEY_STATE(KEY_SHIFT);
					k=0, kend=_countof(rct_custom_params);
					for(;k<kend;++k)
					{
						char *end=0;
						for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);
						rct_custom_params[k]=(short)strtol((char*)text->data+idx, &end, 16);
						for(;idx<(int)text->count&&!isspace(text->data[idx]);++idx);
						if(idx>=(int)text->count)
							goto paste_finish;
					}
					k=0;
					kend=_countof(custom_params);
					if(shift)
					{
						for(;k<kend;++k)
						{
							for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);
							custom_params[k]=(int)round(atof((char*)text->data+idx)*65536);
							for(;idx<(int)text->count&&!isspace(text->data[idx]);++idx);
							if(idx>=(int)text->count)
								goto paste_finish;
						}
					}
					else
						idx=parse_nvals_i32(text, idx, custom_params, _countof(custom_params));
				}

			paste_finish:
				array_free(&text);
				update_image();
				return 1;
			}
		}
		if(GET_KEY_STATE(KEY_CTRL))//paste bitmap
		{
			Image *im2=paste_bmp_from_clipboard();
			if(im2)
			{
				free(im0);
				im0=im2;
				set_window_title("From clipboard - eBench");
				update_image();
				return 1;
			}
		}
		break;
	case KEY_1:
	case KEY_2:
	case KEY_3:
		if(GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)
		{
			int src=key-KEY_1;
			if(src!=custom_pred_ch_idx)
			{
				memcpy(custom_params+(_countof(custom_params)/3)*custom_pred_ch_idx, custom_params+(_countof(custom_params)/3)*src, sizeof(custom_params)/3);
				update_image();
				return 1;
			}
		}
		break;
	case 'M':
		if(im1)
		{
			int shift=GET_KEY_STATE(KEY_SHIFT);
			MODVAR(mode, mode+1-(shift<<1), VIS_COUNT);
			update_image();
		}
		return 1;
	case 'H':
		//if(GET_KEY_STATE(KEY_CTRL))
		{
			memset(combCRhist, 0, sizeof(combCRhist));
			combCRhist_idx=0;
			combCRhist_max=1;
			return 1;
		}
		break;
	//case 'X':
	//	{
	//		extrainfo=!extrainfo;
	//		return 1;
	//	}
	//	break;
	//case 'Z'://TODO show neighbor pixels around cursor
	//	break;
	//case 'P':
	//	if(GET_KEY_STATE(KEY_CTRL))
	//		test_predmask(im1);
	//	break;
	case 'B'://batch test
		if(GET_KEY_STATE(KEY_CTRL))
		{
			//int nthreads=4;
			//if(GET_KEY_STATE(KEY_1))nthreads=1;
			//else if(GET_KEY_STATE(KEY_2))nthreads=2;
			//else if(GET_KEY_STATE(KEY_3))nthreads=3;
			//else if(GET_KEY_STATE(KEY_4))nthreads=4;
			//else if(GET_KEY_STATE(KEY_5))nthreads=5;
			//else if(GET_KEY_STATE(KEY_6))nthreads=6;
			//else if(GET_KEY_STATE(KEY_7))nthreads=7;
			//else if(GET_KEY_STATE(KEY_8))nthreads=8;
			//else if(GET_KEY_STATE(KEY_9))nthreads=9;
			batch_test();
		}
		break;
	case KEY_LBRACKET:
	case KEY_RBRACKET:
		if(transforms_customenabled)
		{
			custom_pred_ch_idx+=((key==KEY_RBRACKET)<<1)-1;
			MODVAR(custom_pred_ch_idx, custom_pred_ch_idx, 3);
			return 1;
		}
		break;
	case 'N':
		if(GET_KEY_STATE(KEY_CTRL))//Ctrl N	add noise to params
		{
			if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				srand((unsigned)__rdtsc());
				short *p=(short*)&c3_params;
				for(int k=0;k<C3_NPARAMS;++k)
					p[k]+=rand()%3-1;
				update_image();
			}
			return 1;
		}
		break;
	case KEY_SPACE:
	//case 'B':
	//case 'N':
		if(im1)
		{
#if 0
			if(mode==VIS_IMAGE_E24)
			{
				float yoffset=tdy*3, half=blocksize*0.5f;
				float x1=blockmx-half, x2=blockmx+half, y1=blockmy-yoffset-half, y2=blockmy-yoffset+half;
				if(blockmx>=0&&blockmx<iw&&blockmy>=yoffset&&blockmy<yoffset+ih)
				{
					if(iw<blocksize)
						x1=0, x2=(float)iw;
					if(x1<0)
						x1=0, x2=(float)blocksize;
					if(x2>iw)
						x1=(float)(iw-blocksize), x2=(float)iw;
				
					if(ih<blocksize)
						y1=0, y2=(float)ih;
					if(y1<0)
						y1=0, y2=(float)blocksize;
					if(y2>ih)
						y1=(float)(ih-blocksize), y2=(float)ih;

					e24_optimizeall(image, iw, ih, (int)roundf(x1), (int)roundf(x2), (int)roundf(y1), (int)roundf(y2), 1);
				}
			}
			else if(transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC])
			{
				int res=iw*ih;
				char *buf2=(char*)malloc((size_t)res<<2);
				if(!buf2)
				{
					LOG_ERROR("Allocation error");
					return 0;
				}
				memcpy(buf2, im0, (size_t)res<<2);
				addhalf((unsigned char*)buf2, iw, ih, 3, 4);
				colortransform_ycmcb_fwd(buf2, iw, ih);
#if 1
				float info[4]={0};
				logic_opt_checkonthread(info);
				if(info[0]<0)
				{
					int kc=customparam_ch_idx/2;
					logic_opt(buf2, iw, ih, kc, logic_params+LOGIC_PARAMS_PER_CH*kc);
				}
				else
					free(buf2);
#else
				free(buf2);
				srand((unsigned)__rdtsc());
				short *params=logic_params+LOGIC_PARAMS_PER_CH*customparam_ch_idx;
				for(int k=0;k<LOGIC_PARAMS_PER_CH;++k)
					params[k]=rand();
				update_image();
#endif
			}
			else
#endif
			if(GET_KEY_STATE(KEY_CTRL))
			{
				int success=1;
				if(im0->iw!=im1->iw||im0->ih!=im1->ih||im0->nch!=im1->nch)
				{
					messagebox(MBOX_OK, "Dimension Error",
						"Image dimensions changed:\n"
						"\tim0\tCWH %d*%d*%d\n"
						"\tim1\tCWH %d*%d*%d",
						im0->nch, im0->iw, im0->ih,
						im1->nch, im1->iw, im1->ih
					);
					success=0;
				}
				else
				{
					for(ptrdiff_t k=0, res=(ptrdiff_t)im0->iw*im0->ih;k<res;++k)
					{
						if(memcmp(im1->data+(k<<2), im0->data+(k<<2), sizeof(int[4])))
						{
							int kx=(int)(k%im0->iw), ky=(int)(k/im0->iw);
							messagebox(MBOX_OK, "Pixel Error",
								"Difference at XY %d %d:\n"
								"\terror\toriginal\n"
								"C0\t0x%04X\t0x%04X\n"
								"C1\t0x%04X\t0x%04X\n"
								"C2\t0x%04X\t0x%04X\n"
								"C3\t0x%04X\t0x%04X",
								kx, ky,
								(unsigned short)im1->data[k<<2|0], (unsigned short)im0->data[k<<2|0],
								(unsigned short)im1->data[k<<2|1], (unsigned short)im0->data[k<<2|1],
								(unsigned short)im1->data[k<<2|2], (unsigned short)im0->data[k<<2|2],
								(unsigned short)im1->data[k<<2|3], (unsigned short)im0->data[k<<2|3]
							);
							success=0;
							break;
						}
					}
				}
				if(success)
					messagebox(MBOX_OK, "SUCCESS", "The image is bit-exact.");
			}
			else if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
			{
				rct_custom_optimize(im0, rct_custom_params);
				update_image();
			}
			else if(transforms_mask[ST_FWD_CUSTOM]||transforms_mask[ST_INV_CUSTOM])
			{
				Image *im2=0;
				image_copy(&im2, im0);
				if(!im2)
					return 0;
				apply_selected_transforms(im2, 1);
				pred_custom_optimize(im2, custom_params);
				free(im2);
				update_image();
			}
#if 0
			else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
			{
				if(GET_KEY_STATE(KEY_CTRL))
					custom2_opt_batch(&c2_params);
					//memset(&c2_params, 0, sizeof(c2_params));
				else
				{
					int res=iw*ih;
					char *buf2=(char*)malloc((size_t)res<<2);
					if(!buf2)
					{
						LOG_ERROR("Allocation error");
						return 0;
					}
					memcpy(buf2, im0, (size_t)res<<2);
					addhalf((unsigned char*)buf2, iw, ih, 3, 4);
					colortransform_ycmcb_fwd(buf2, iw, ih);//
					//pred_grad_fwd(buf2, iw, ih, 3, 4);//

					if(GET_KEY_STATE(KEY_SHIFT))
						custom2_opt_blocks(buf2, iw, ih, &c2_params);
					else
					{
						int maskbits=10;
							 if(keyboard[KEY_1])maskbits=1;
						else if(keyboard[KEY_2])maskbits=2;
						else if(keyboard[KEY_3])maskbits=3;
						else if(keyboard[KEY_4])maskbits=4;
						else if(keyboard[KEY_5])maskbits=5;
						else if(keyboard[KEY_6])maskbits=6;
						else if(keyboard[KEY_7])maskbits=7;
						else if(keyboard[KEY_8])maskbits=8;
						else if(keyboard[KEY_9])maskbits=9;
						custom2_opt(buf2, iw, ih, &c2_params, 0, maskbits, 1, 0);
					}

					free(buf2);
				}
				update_image();
			}
#endif
			else if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				//int maskbits=10;
				int maskbits=3;
				if(keyboard[KEY_1])maskbits=1;
				else if(keyboard[KEY_2])maskbits=2;
				else if(keyboard[KEY_3])maskbits=3;
				else if(keyboard[KEY_4])maskbits=4;
				else if(keyboard[KEY_5])maskbits=5;
				else if(keyboard[KEY_6])maskbits=6;
				else if(keyboard[KEY_7])maskbits=7;
				else if(keyboard[KEY_8])maskbits=8;
				else if(keyboard[KEY_9])maskbits=9;
				//if(GET_KEY_STATE(KEY_CTRL))
				//	custom3_opt_batch(&c3_params, 0, maskbits, 1, 0);
				//else if(GET_KEY_STATE(KEY_SHIFT))
				//	custom3_opt_batch2(&c3_params, 0, maskbits, 1, 0);
				//else
				{
					Image *im2=0;
					image_copy(&im2, im0);
					if(!im2)
						return 0;
					apply_selected_transforms(im2, 1);
					//if(GET_KEY_STATE(KEY_CTRL))
					//	colortransform_YCbCr_R_v1(im2, 1);
					custom3_opt(im2, &c3_params, 0, maskbits, 1, 0);
					free(im2);
				}
				update_image();
			}
#if 0
			else if(transforms_mask[ST_FWD_CUSTOM4]||transforms_mask[ST_INV_CUSTOM4])
			{
				int maskbits=10;
					 if(keyboard[KEY_1])maskbits=1;
				else if(keyboard[KEY_2])maskbits=2;
				else if(keyboard[KEY_3])maskbits=3;
				else if(keyboard[KEY_4])maskbits=4;
				else if(keyboard[KEY_5])maskbits=5;
				else if(keyboard[KEY_6])maskbits=6;
				else if(keyboard[KEY_7])maskbits=7;
				else if(keyboard[KEY_8])maskbits=8;
				else if(keyboard[KEY_9])maskbits=9;
				//if(GET_KEY_STATE(KEY_CTRL))
				//	custom3_opt_batch(&c3_params, 0, maskbits, 1, 0);
				//else if(GET_KEY_STATE(KEY_SHIFT))
				//	custom3_opt_batch2(&c3_params, 0, maskbits, 1, 0);
				//else
				{
					int res=iw*ih;
					char *buf2=(char*)malloc((size_t)res<<2);
					if(!buf2)
					{
						LOG_ERROR("Allocation error");
						return 0;
					}
					memcpy(buf2, im0, (size_t)res<<2);
					addhalf((unsigned char*)buf2, iw, ih, 3, 4);
					colortransform_ycmcb_fwd(buf2, iw, ih);//

					custom4_opt(buf2, iw, ih, &c4_params, 0, maskbits, 1, 0);

					free(buf2);
				}
				update_image();
			}
#endif
			else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP]||transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
			{
				int jxl=transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP];
				int pw2=transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM];
				if(GET_KEY_STATE(KEY_CTRL))
				{
					if(jxl)
					{
						for(int k=0;k<33;++k)
							jxlparams_i16[k]=(k%11<4)<<12;
					}
					else
						memset(pw2_params, 0, sizeof(pw2_params));
				}
				else
				{
					Image *im2=0;
					image_copy(&im2, im0);
					if(!im2)
						return 0;

					//int step=256;
					//if(keyboard['1'])step=128;
					//else if(keyboard['2'])step=64;
					//else if(keyboard['3'])step=32;
					//else if(keyboard['4'])step=16;
					//else if(keyboard['5'])step= 8;
					//else if(keyboard['6'])step= 4;
					//else if(keyboard['7'])step= 2;
					//else if(keyboard['8'])step= 1;

					if(jxl)
					{
						colortransform_YCbCr_R_v1(im2, 1);
						pred_jxl_opt_v2(im2, jxlparams_i16, 1);
					}
					else if(pw2)
					{
						colortransform_YCbCr_R_v1(im2, 1);
						pred_w2_opt_v2(im2, pw2_params, 1);
					}
					free(im2);
				}
				update_image();
			}
			return 1;
		}
		break;
	default://to make gcc -Wall happy
		break;
	}
	return 0;
}
int io_keyup(IOKey key, char c)
{
	(void)c;

	switch(key)
	{
	//case KEY_LBUTTON:
	//case KEY_MBUTTON:
	//case KEY_RBUTTON:
	//	printf("Declick at (%d, %d)\n", mx, my);
	//	break;
		
	case KEY_SPACE:
#define AK(KEY) case KEY:
	ACTIVE_KEY_LIST
#undef  AK
		count_active_keys(key);
		break;

	default://to make gcc -Wall happy
		break;
	//	printf("%02X %02X=%c up\n", key, c, c);
	//	if(key=='A')
	//		timer_stop();
	//	break;
	}
	return 0;
}
void io_timer(void)
{
	float move_speed=keyboard[KEY_SHIFT]?10*cam.move_speed:cam.move_speed;
	if(keyboard['W'])	cam_moveForward(cam, move_speed);
	if(keyboard['A'])	cam_moveLeft(cam, move_speed);
	if(keyboard['S'])	cam_moveBack(cam, move_speed);
	if(keyboard['D'])	cam_moveRight(cam, move_speed);
	if(keyboard['T'])	cam_moveUp(cam, move_speed);
	if(keyboard['G'])	cam_moveDown(cam, move_speed);
	if(keyboard[KEY_UP])	cam_turnUp(cam, key_turn_speed);
	if(keyboard[KEY_DOWN])	cam_turnDown(cam, key_turn_speed);
	if(keyboard[KEY_LEFT])	cam_turnLeft(cam, key_turn_speed);
	if(keyboard[KEY_RIGHT])	cam_turnRight(cam, key_turn_speed);
	if(keyboard[KEY_ENTER])	cam_zoomIn(cam, 1.1f);
	if(keyboard[KEY_BKSP])	cam_zoomOut(cam, 1.1f);
}
static float print_i16_row(float x, float y, float zoom, const short *row, int count)
{
	int printed=0;
	for(int k=0;k<count;++k)
	{
		short val=row[k];
		printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, " %c%4X", val<0?'-':' ', abs(val));
	}
	return print_line_immediate(0, x, y, zoom, g_buf, printed, -1, 0, 0);
}

static void draw_shape_i(const float *points)
{
	const float *axes=points+3*8;
	draw_3d_line(&cam, points+3*0, points+3*1, 0xFF0000FF);
	draw_3d_line(&cam, points+3*2, points+3*3, 0xFF0000C0);
	draw_3d_line(&cam, points+3*4, points+3*5, 0xFF000080);
	draw_3d_line(&cam, points+3*6, points+3*7, 0xFF000040);
	draw_3d_line(&cam, points+3*0, points+3*2, 0xFF00FF00);
	draw_3d_line(&cam, points+3*1, points+3*3, 0xFF00C000);
	draw_3d_line(&cam, points+3*4, points+3*6, 0xFF008000);
	draw_3d_line(&cam, points+3*5, points+3*7, 0xFF004000);
	draw_3d_line(&cam, points+3*0, points+3*4, 0xFFFF0000);
	draw_3d_line(&cam, points+3*1, points+3*5, 0xFFC00000);
	draw_3d_line(&cam, points+3*2, points+3*6, 0xFF800000);
	draw_3d_line(&cam, points+3*3, points+3*7, 0xFF400000);
	draw_3d_line(&cam, axes, axes+3*1, 0xFF0000FF);
	draw_3d_line(&cam, axes, axes+3*2, 0xFF00FF00);
	draw_3d_line(&cam, axes, axes+3*3, 0xFFFF0000);
	draw_3d_line(&cam, axes, axes+3*4, 0xFF000000);
}
static void draw_YCoCg_R(float x, float y, float z)
{
	float points[]=
	{
		x-1, y-1, z-1,
		x+1, y-1, z-1,
		x-1, y+1, z-1,
		x+1, y+1, z-1,
		x-1, y-1, z+1,
		x+1, y-1, z+1,
		x-1, y+1, z+1,
		x+1, y+1, z+1,

		x, y, z,
		x+1, y, z,
		x, y+1, z,
		x, y, z+1,
		x+2, y+2, z+2,
	};
	draw_shape_i(points);
	for(int k=0;k<_countof(points);k+=3)
	{
		float *p=points+k;
		p[0]-=x;
		p[1]-=y;
		p[2]-=z;

		p[0]-=p[2];
		p[2]+=p[0]*0.5f;
		p[1]-=p[2];
		p[2]+=p[1]*0.5f;

		p[0]+=x;
		p[1]+=y;
		p[2]+=z;
	}
	draw_shape_i(points);
}
static void draw_YCbCr_R(float x, float y, float z)
{
	float points[]=
	{
		x-1, y-1, z-1,
		x+1, y-1, z-1,
		x-1, y+1, z-1,
		x+1, y+1, z-1,
		x-1, y-1, z+1,
		x+1, y-1, z+1,
		x-1, y+1, z+1,
		x+1, y+1, z+1,

		x, y, z,
		x+1, y, z,
		x, y+1, z,
		x, y, z+1,
		x+2, y+2, z+2,
	};
	draw_shape_i(points);
	for(int k=0;k<_countof(points);k+=3)
	{
		float *p=points+k;
		p[0]-=x;
		p[1]-=y;
		p[2]-=z;

		p[0]-=p[1];
		p[1]+=p[0]*0.5f;
		p[2]-=p[1];
		p[1]+=p[2]*0.5f;

		p[0]+=x;
		p[1]+=y;
		p[2]+=z;
	}
	draw_shape_i(points);
}
static void draw_JPEG2000(float x, float y, float z)
{
	float points[]=
	{
		x-1, y-1, z-1,
		x+1, y-1, z-1,
		x-1, y+1, z-1,
		x+1, y+1, z-1,
		x-1, y-1, z+1,
		x+1, y-1, z+1,
		x-1, y+1, z+1,
		x+1, y+1, z+1,

		x, y, z,
		x+1, y, z,
		x, y+1, z,
		x, y, z+1,
		x+2, y+2, z+2,
	};
	draw_shape_i(points);
	for(int k=0;k<_countof(points);k+=3)
	{
		float *p=points+k;
		p[0]-=x;
		p[1]-=y;
		p[2]-=z;

		p[0]-=p[1];
		p[2]-=p[1];
		p[1]+=(p[0]+p[2])*0.25f;

		p[0]+=x;
		p[1]+=y;
		p[2]+=z;
	}
	draw_shape_i(points);
}
static void draw_diffav(float x, float y, float z)
{
	float points[]=
	{
		x-1, y-1, z,
		x+1, y-1, z,
		x+1, y+1, z,
		x-1, y+1, z,
	};
	draw_3d_line(&cam, points+3*0, points+3*1, 0xFF000000);
	draw_3d_line(&cam, points+3*1, points+3*2, 0xFF0000FF);
	draw_3d_line(&cam, points+3*2, points+3*3, 0xFF00FF00);
	draw_3d_line(&cam, points+3*3, points+3*0, 0xFFFF0000);
	for(int k=0;k<_countof(points)-2;k+=3)
	{
		points[k+0]-=points[k+1];
		points[k+1]+=points[k+0]/2;
	}
	draw_3d_line(&cam, points+3*0, points+3*1, 0xFF000000);
	draw_3d_line(&cam, points+3*1, points+3*2, 0xFF0000FF);
	draw_3d_line(&cam, points+3*2, points+3*3, 0xFF00FF00);
	draw_3d_line(&cam, points+3*3, points+3*0, 0xFFFF0000);
}

void io_render(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	if(!h)
		return;

	float axes[]=
	{
		0, 0, 0,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};
	draw_3d_line(&cam, axes, axes+3, 0xFF0000FF);
	draw_3d_line(&cam, axes, axes+6, 0xFF00FF00);
	draw_3d_line(&cam, axes, axes+9, 0xFFFF0000);
#if 0
	if(!im0)
	{
		draw_YCoCg_R(10, 0, 0);
		draw_YCbCr_R(20, 0, 0);
		draw_JPEG2000(30, 0, 0);
		draw_diffav(0, 0, -10);
	}
#endif

	if(im1)
	{
		switch(mode)
		{
		//case VIS_PLANES:		chart_planes_draw();	break;
		//case VIS_MESH:		chart_mesh_draw();	break;
		//case VIS_MESH_SEPARATE:	chart_mesh_sep_draw();	break;
		case VIS_HISTOGRAM:
			{
				float yoffset=tdy*3;
				display_texture_i(0, im1->iw, (int)yoffset, (int)yoffset+im1->ih, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
				chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0, 0x60, hist, histmax);
			}
			break;
		case VIS_JOINT_HISTOGRAM:	chart_jointhist_draw();	break;
		case VIS_IMAGE:
		case VIS_ZIPF:
			{
				//int waitstatus=0;
				//if(ghMutex)
				//	waitstatus=WaitForSingleObject(ghMutex, INFINITE);
				float yoffset=tdy*3;
				if(show_full_image)
					display_texture_i(0, w, 0, h, (int*)(mode==VIS_ZIPF?zimage:im_export), im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
				else
					display_texture_i(0, im1->iw, (int)yoffset, (int)yoffset+im1->ih, (int*)(mode==VIS_ZIPF?zimage:im_export), im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
				//if(waitstatus==WAIT_OBJECT_0)
				//	ReleaseMutex(ghMutex);
			}
			break;
#if 0
		case VIS_BAYES:
			{
				bayes_update();
				if(*bayes_mem)
				{
					static ArrayHandle vertices=0;
					int texture[]=
					{
						0xFF0000FF, 0xFF0000FF,
						0xFF0000FF, 0xFF0000FF,
					};
					for(int kb=7;kb>=0;--kb)
					{
						ArrayHandle mem=bayes_mem[kb];
						int x=7-kb;
						for(int kctx=0, nctx=1<<x;kctx<nctx;++kctx)
						{
							for(int kp=0;kp<256;++kp)
							{
								BayesCounter *node=(BayesCounter*)array_at(&mem, kctx<<8|kp);
								float z=(float)node->n[1]/(node->n[0]+node->n[1]);
								draw_3d_line_enqueue(&vertices, (float)x*4, (float)kp/16, z+kctx);
							}
							draw_3d_flush(vertices, &cam, texture, 2, 2, 0, GL_LINE_STRIP);
						}
					}
				}
				float yoffset=tdy*3, half=blocksize*0.5f;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 0, 1, 0, 1, 0.5f, 0);
				float x1=blockmx-half, x2=blockmx+half, y1=blockmy-yoffset-half, y2=blockmy-yoffset+half;

				//if(range_intersect(x1, x2, 0, iw)&&range_intersect(y1, y2, 0, ih))
				if(blockmx>=0&&blockmx<iw&&blockmy>=yoffset&&blockmy<yoffset+ih)
				{
					if(iw<blocksize)
						x1=0, x2=(float)iw;
					if(x1<0)
						x1=0, x2=(float)blocksize;
					if(x2>iw)
						x1=(float)(iw-blocksize), x2=(float)iw;
				
					if(ih<blocksize)
						y1=0, y2=(float)ih;
					if(y1<0)
						y1=0, y2=(float)blocksize;
					if(y2>ih)
						y1=(float)(ih-blocksize), y2=(float)ih;
					
					int boxcolor=0xFFFFFF00;
					y1+=yoffset;
					y2+=yoffset;
					draw_rect_hollow(x1, x2, y1, y2, boxcolor);
				}
			}
			break;
		case VIS_IMAGE_BLOCK:
			{
				float yoffset=tdy*3, half=blocksize*0.5f;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 0, 1, 0, 1, 1, 0);
				float x1=blockmx-half, x2=blockmx+half, y1=blockmy-yoffset-half, y2=blockmy-yoffset+half;

				//if(range_intersect(x1, x2, 0, iw)&&range_intersect(y1, y2, 0, ih))
				if(blockmx>=0&&blockmx<iw&&blockmy>=yoffset&&blockmy<yoffset+ih)
				{
					if(iw<blocksize)
						x1=0, x2=(float)iw;
					if(x1<0)
						x1=0, x2=(float)blocksize;
					if(x2>iw)
						x1=(float)(iw-blocksize), x2=(float)iw;
				
					if(ih<blocksize)
						y1=0, y2=(float)ih;
					if(y1<0)
						y1=0, y2=(float)blocksize;
					if(y2>ih)
						y1=(float)(ih-blocksize), y2=(float)ih;
					
					memset(histmax, 0, sizeof(histmax));
					memset(histmax2, 0, sizeof(histmax2));
					memset(hist, 0, sizeof(hist));
					memset(hist2, 0, sizeof(hist2));
					chart_hist_update(image, iw, ih, (int)x1, (int)x2, (int)y1, (int)y2, hist, histmax, blockCR);
					chart_hist_update(image, iw, ih, (int)x1-margin, (int)x1, (int)y1, (int)y2, hist2, histmax2, 0);
					chart_hist_update(image, iw, ih, (int)x1-margin, (int)x2+margin, (int)y1-margin, (int)y1, hist2, histmax2, 0);
					chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0x80808080, 0, hist2, histmax2);
					chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0, 0x30, hist, histmax);

					int boxcolor=0xFFFFFF00;
					y1+=yoffset;
					y2+=yoffset;
					draw_rect_hollow(x1, x2, y1, y2, boxcolor);
					draw_line(x1, y2, x1-margin, y2, boxcolor);
					draw_line(x1-margin, y2, x1-margin, y1-margin, boxcolor);
					draw_line(x1-margin, y1-margin, x2+margin, y1-margin, boxcolor);
					draw_line(x2+margin, y1-margin, x2+margin, y1, boxcolor);
					draw_line(x2+margin, y1, x2, y1, boxcolor);
					GUIPrint(0, 0, tdy*2, 1, "TRGB %8f [%8f %8f %8f] block %d margin %d", 3/(1/blockCR[0]+1/blockCR[1]+1/blockCR[2]), blockCR[0], blockCR[1], blockCR[2], blocksize, margin);
				}
			}
			//{
			//	//int x=(int)ceilf(tdx*6.5f);
			//	int x=(int)floorf(tdx*6.5f);
			//	display_texture_i(x, x+iw, (int)(tdy*3), (int)(tdy*3+ih), (int*)image, iw, ih, 1);
			//}
			break;
		case VIS_IMAGE_E24:
			{
				float yoffset=tdy*3, half=blocksize*0.5f;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 1, 0, 1, 0, 1);
				float x1=blockmx-half, x2=blockmx+half, y1=blockmy-yoffset-half, y2=blockmy-yoffset+half;
				if(blockmx>=0&&blockmx<iw&&blockmy>=yoffset&&blockmy<yoffset+ih)
				{
					if(iw<blocksize)
						x1=0, x2=(float)iw;
					if(x1<0)
						x1=0, x2=(float)blocksize;
					if(x2>iw)
						x1=(float)(iw-blocksize), x2=(float)iw;
				
					if(ih<blocksize)
						y1=0, y2=(float)ih;
					if(y1<0)
						y1=0, y2=(float)blocksize;
					if(y2>ih)
						y1=(float)(ih-blocksize), y2=(float)ih;
					
					e24_estimate(image, iw, ih, (int)roundf(x1), (int)roundf(x2), (int)roundf(y1), (int)roundf(y2));

					int boxcolor=0xFFFFFF00;
					y1+=yoffset;
					y2+=yoffset;
					draw_rect_hollow(x1, x2, y1, y2, boxcolor);
					float total_cr=(float)(3/(1/e24_cr[0]+1/e24_cr[1]+1/e24_cr[2])), scale=128;

					draw_line(x1+scale, y1-32, x1+scale, y1, 0xC0000000);
					draw_rect(x1, x1+total_cr*scale, y1-2, y1-8, 0xC0000000);
					draw_rect(x1, x1+(float)e24_cr[0]*scale, y1-10, y1-16, 0xC00000FF);
					draw_rect(x1, x1+(float)e24_cr[1]*scale, y1-18, y1-24, 0xC000FF00);
					draw_rect(x1, x1+(float)e24_cr[2]*scale, y1-26, y1-32, 0xC0FF0000);
#if 0
					draw_line(x1+scale, y1, x1+scale, y2, 0x80000000);
					//draw_rect(x1, x1+scale, y1-8, y1, 0x80FF80FF);
					draw_line(x1, y1-2, x1+total_cr*scale, y1-2, 0xFF000000);
					draw_line(x1, y1-4, x1+(float)e24_cr[0]*scale, y1-4, 0xFF0000FF);
					draw_line(x1, y1-6, x1+(float)e24_cr[1]*scale, y1-6, 0xFF00FF00);
					draw_line(x1, y1-8, x1+(float)e24_cr[2]*scale, y1-8, 0xFFFF0000);
#endif
					GUIPrint(0, 0, tdy*2, 1, "E24 TRGB %8f [%8lf %8lf %8lf]  block %d", total_cr, e24_cr[0], e24_cr[1], e24_cr[2], blocksize);
					for(int k=0;k<3;++k)
					{
						E24Params *p=e24_params+k;
						GUIPrint(0, (float)(w>>1)-100, (float)(h>>2)+(k+1)*tdy, 1, "W %3d  MLTR %3d %3d %3d A 0x%02X I %3d", p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
					}
				}
			}
			break;
		case VIS_DWT_BLOCK:
			{
				int w2=iw, h2=ih;//do not scale
				for(int kc=0;kc<3;++kc)
				{
					int kx=kc&1, ky=kc>>1;
					display_texture(kx*w2, (kx+1)*w2, ky*h2, (ky+1)*h2, image_txid[1], 1, 0, 1, kc/3.f, (kc+1)/3.f);
				}
				float half=blocksize*0.5f;
				for(int kc=0;kc<3;++kc)
				{
					int kx=kc&1, ky=kc>>1;
					int
						x1=kx*w2, x2=(kx+1)*w2,
						y1=ky*h2, y2=(ky+1)*h2;
					if(blockmx>=x1&&blockmx<x2&&blockmy>=y1&&blockmy<y2)
					{
						float
							mx1=blockmx-half, mx2=blockmx+half,
							my1=blockmy-half, my2=blockmy+half;

						if(w2<blocksize)
							mx1=(float)x1, mx2=(float)x2;
						if(mx1<x1)
							mx1=(float)x1, mx2=(float)(x1+blocksize);
						if(mx2>x2)
							mx1=(float)(x2-blocksize), mx2=(float)x2;
				
						if(h2<blocksize)
							my1=(float)y1, my2=(float)y2;
						if(my1<y1)
							my1=(float)y1, my2=(float)(y1+blocksize);
						if(my2>y2)
							my1=(float)(y2-blocksize), my2=(float)y2;
						
						float px=(x1+x2)*0.5f, py=(y1+y2)*0.5f;
						if(mx1<px&&my1<py)
						{
							//if(fabsf(my1/mx1)<fabsf(py/px))
							if(fabsf(my1*px)<fabsf(py*mx1))
								mx1=(float)px, mx2=px+blocksize;
							else
								my1=(float)py, my2=py+blocksize;
						}

						DRAW_LINEI(x1, py, x2, py, 0x80FFFFFF);
						DRAW_LINEI(px, y1, px, y2, 0x80FFFFFF);
						float px2=x1+(px-x1)*0.5f, py2=y1+(py-y1)*0.5f;
						DRAW_LINEI(x1, py2, px, py2, 0x80FFFFFF);
						DRAW_LINEI(px2, y1, px2, py, 0x80FFFFFF);

						int histcolor=0x40FF00FF;
						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 2, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 2, 3, histcolor, 0, hist, histmax);

						mx1=x1+(mx1-x1)*0.5f;
						mx2=x1+(mx2-x1)*0.5f;
						my1=y1+(my1-y1)*0.5f;
						my2=y1+(my2-y1)*0.5f;

						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 1, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 1, 2, histcolor, 0, hist, histmax);
						
						mx1=x1+(mx1-x1)*0.5f;
						mx2=x1+(mx2-x1)*0.5f;
						my1=y1+(my1-y1)*0.5f;
						my2=y1+(my2-y1)*0.5f;

						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 0, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 0, 1, histcolor, 0, hist, histmax);

						GUIPrint(0, 0, tdy*2, 1, "[%8f %8f %8f] %g", blockCR[0], blockCR[1], blockCR[2], blocksize);
						break;
					}
				}
			}
			//display_texture(0,  iw,    0,  ih,    image_txid[1], 1, 0, 1, 0,     1.f/3);
			//display_texture(iw, iw<<1, 0,  ih,    image_txid[1], 1, 0, 1, 1.f/3, 2.f/3);
			//display_texture(0,  iw,    ih, ih<<1, image_txid[1], 1, 0, 1, 2.f/3, 1);
			break;
#endif
		case VIS_IMAGE_TRICOLOR:
			if(show_full_image)
			{
				display_texture(0,     w/3,   0, h, txid_separate_r, 1, 0, 1, 0, 1);
				display_texture(w/3,   w*2/3, 0, h, txid_separate_g, 1, 0, 1, 0, 1);
				display_texture(w*2/3, w,     0, h, txid_separate_b, 1, 0, 1, 0, 1);
			}
			else
			{
				display_texture(0,       im1->iw,    0,       im1->ih,    txid_separate_r, 1, 0, 1, 0, 1);
				display_texture(im1->iw, im1->iw<<1, 0,       im1->ih,    txid_separate_g, 1, 0, 1, 0, 1);
				display_texture(0,       im1->iw,    im1->ih, im1->ih<<1, txid_separate_b, 1, 0, 1, 0, 1);
				if(im1->depth[3])
					display_texture(im1->iw, im1->iw<<1, im1->ih, im1->ih<<1, txid_separate_a, 1, 0, 1, 0, 1);
			}
			break;
		}
	}

	if(transforms_customenabled)
	{
		float ystep=tdy*guizoom, x, y;
		if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
		{
			//custom color transform params
			x=buttons[0].x1;
			y=buttons[0].y1;
			//P0  rgb
			//r += (-0x00*g-0x00*b)>>6
			//g += (-0x00*r-0x00*b)>>6
			//b += (-0x00*r-0x00*g)>>6
			//g += (-0x00*r-0x00*b)>>6
			//0123456789012345678901234
			const char chnames[]="rgb";
			unsigned char per[4]={0};
			rct_custom_unpackpermutation(rct_custom_params[8], per);
			GUIPrint(0, x, y+ystep*0, guizoom, "P%d  %c%c%c", rct_custom_params[8], chnames[per[0]], chnames[per[1]], chnames[per[2]]);
			GUIPrint(0, x, y+ystep*1, guizoom, "%c += (%c0x%04X*%c%c0x%04X*%c)>>%d",
				chnames[per[0]],
				rct_custom_params[0]<0?'-':' ', abs(rct_custom_params[0]), chnames[per[1]],
				rct_custom_params[1]<0?'-':'+', abs(rct_custom_params[1]), chnames[per[2]],
				RCT_CUSTOM_PARAMBITS
			);
			GUIPrint(0, x, y+ystep*2, guizoom, "%c += (%c0x%04X*%c%c0x%04X*%c)>>%d",
				chnames[per[1]],
				rct_custom_params[2]<0?'-':' ', abs(rct_custom_params[2]), chnames[per[0]],
				rct_custom_params[3]<0?'-':'+', abs(rct_custom_params[3]), chnames[per[2]],
				RCT_CUSTOM_PARAMBITS
			);
			GUIPrint(0, x, y+ystep*3, guizoom, "%c += (%c0x%04X*%c%c0x%04X*%c)>>%d",
				chnames[per[2]],
				rct_custom_params[4]<0?'-':' ', abs(rct_custom_params[4]), chnames[per[0]],
				rct_custom_params[5]<0?'-':'+', abs(rct_custom_params[5]), chnames[per[1]],
				RCT_CUSTOM_PARAMBITS
			);
			GUIPrint(0, x, y+ystep*4, guizoom, "%c += (%c0x%04X*%c%c0x%04X*%c)>>%d",
				chnames[per[1]],
				rct_custom_params[6]<0?'-':' ', abs(rct_custom_params[6]), chnames[per[0]],
				rct_custom_params[7]<0?'-':'+', abs(rct_custom_params[7]), chnames[per[2]],
				RCT_CUSTOM_PARAMBITS
			);
#if 0
			//0000000000111111111122222222223333
			//0123456789012345678901234567890123
			//r-=g
			//g+=(-0x00.0000*r-0x00.0000*b)>>16
			//b-=g
			//g+=(-0x00.0000*r-0x00.0000*b)>>16
			GUIPrint(0, x, y+ystep*0, guizoom, "r-=g");//do not change these strings!
			GUIPrint(0, x, y+ystep*1, guizoom, "g+=(%c0x%02X.%04X*r%c0x%02X.%04X*b)>>16", rct_custom_params[0]<0?'-':' ', abs(rct_custom_params[0])>>16, abs(rct_custom_params[0])&0xFFFF, rct_custom_params[1]<0?'-':'+', abs(rct_custom_params[1])>>16, abs(rct_custom_params[1])&0xFFFF);
			GUIPrint(0, x, y+ystep*2, guizoom, "b-=g");
			GUIPrint(0, x, y+ystep*3, guizoom, "g+=(%c0x%02X.%04X*r%c0x%02X.%04X*b)>>16", rct_custom_params[2]<0?'-':' ', abs(rct_custom_params[2])>>16, abs(rct_custom_params[2])&0xFFFF, rct_custom_params[3]<0?'-':'+', abs(rct_custom_params[3])>>16, abs(rct_custom_params[3])&0xFFFF);
#endif

			//012345678901234567890123456789
			//r-=(>>nnnN.NNN)g+(  nnnN.NNN)b
			//GUIPrint(0, x, y        , guizoom, "r-=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 0], sel[ 0], customparam_ct[ 0], sel[ 1], sel[ 1], customparam_ct[ 1]);
			//GUIPrint(0, x, y+ystep  , guizoom, "g-=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 2], sel[ 2], customparam_ct[ 2], sel[ 3], sel[ 3], customparam_ct[ 3]);
			//GUIPrint(0, x, y+ystep*2, guizoom, "b-=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[ 4], sel[ 4], customparam_ct[ 4], sel[ 5], sel[ 5], customparam_ct[ 5]);
			//GUIPrint(0, x, y+ystep*3, guizoom, "r+=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 6], sel[ 6], customparam_ct[ 6], sel[ 7], sel[ 7], customparam_ct[ 7]);
			//GUIPrint(0, x, y+ystep*4, guizoom, "g+=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 8], sel[ 8], customparam_ct[ 8], sel[ 9], sel[ 9], customparam_ct[ 9]);
			//GUIPrint(0, x, y+ystep*5, guizoom, "b+=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[10], sel[10], customparam_ct[10], sel[11], sel[11], customparam_ct[11]);
		}
		if(transforms_mask[ST_FWD_CUSTOM]||transforms_mask[ST_INV_CUSTOM])
		{
			int c0=set_bk_color(0x80FFFFFF);
			//custom spatial transform params
			x=buttons[1].x1;
			y=buttons[1].y1;
			//0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
			//0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			GUIPrint(0, x, y, guizoom, "Ch %d  Clamp [%cW %cNW %cN %cNE]",
				custom_pred_ch_idx,
				custom_clamp[0]?'+':'-',
				custom_clamp[1]?'+':'-',
				custom_clamp[2]?'+':'-',
				custom_clamp[3]?'+':'-'
			);
			//GUIPrint(0, x, y-tdy, 1, "Ch %d", custom_pred_ch_idx);
			int *params=custom_params+24*custom_pred_ch_idx;
			GUIPrint(0, x, y+ystep*1, guizoom, "%c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X ",
				params[0]<0?'-':' ', abs(params[0])>>16, abs(params[0])&0xFFFF,
				params[1]<0?'-':' ', abs(params[1])>>16, abs(params[1])&0xFFFF,
				params[2]<0?'-':' ', abs(params[2])>>16, abs(params[2])&0xFFFF,
				params[3]<0?'-':' ', abs(params[3])>>16, abs(params[3])&0xFFFF,
				params[4]<0?'-':' ', abs(params[4])>>16, abs(params[4])&0xFFFF,
				params[5]<0?'-':' ', abs(params[5])>>16, abs(params[5])&0xFFFF,
				params[6]<0?'-':' ', abs(params[6])>>16, abs(params[6])&0xFFFF,
				params[7]<0?'-':' ', abs(params[7])>>16, abs(params[7])&0xFFFF,
				params[8]<0?'-':' ', abs(params[8])>>16, abs(params[8])&0xFFFF,
				params[9]<0?'-':' ', abs(params[9])>>16, abs(params[9])&0xFFFF
			);
			params+=10;
			GUIPrint(0, x, y+ystep*2, guizoom, "%c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X ",
				params[0]<0?'-':' ', abs(params[0])>>16, abs(params[0])&0xFFFF,
				params[1]<0?'-':' ', abs(params[1])>>16, abs(params[1])&0xFFFF,
				params[2]<0?'-':' ', abs(params[2])>>16, abs(params[2])&0xFFFF,
				params[3]<0?'-':' ', abs(params[3])>>16, abs(params[3])&0xFFFF,
				params[4]<0?'-':' ', abs(params[4])>>16, abs(params[4])&0xFFFF,
				params[5]<0?'-':' ', abs(params[5])>>16, abs(params[5])&0xFFFF,
				params[6]<0?'-':' ', abs(params[6])>>16, abs(params[6])&0xFFFF,
				params[7]<0?'-':' ', abs(params[7])>>16, abs(params[7])&0xFFFF,
				params[8]<0?'-':' ', abs(params[8])>>16, abs(params[8])&0xFFFF,
				params[9]<0?'-':' ', abs(params[9])>>16, abs(params[9])&0xFFFF
			);
			params+=10;
			GUIPrint(0, x, y+ystep*3, guizoom, "%c0x%02X.%04X%c0x%02X.%04X %c0x%02X.%04X%c0x%02X.%04X",
				params[0]<0?'-':' ', abs(params[0])>>16, abs(params[0])&0xFFFF,
				params[1]<0?'-':' ', abs(params[1])>>16, abs(params[1])&0xFFFF,
				params[2]<0?'-':' ', abs(params[2])>>16, abs(params[2])&0xFFFF,
				params[3]<0?'-':' ', abs(params[3])>>16, abs(params[3])&0xFFFF
			);
			set_bk_color(c0);

			//0000000000111111111122222222223333333333444444444455555
			//0123456789012345678901234567890123456789012345678901234
			//>>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN
			//int bk0=set_bk_color(0xA060A060);
			//const double *params=customparam_st+12*customparam_ch_idx;
			//if(customparam_ch_idx&1)
			//	GUIPrint(0, x, y-tdy, 1, "Er %d", customparam_ch_idx/2);
			//else
			//	GUIPrint(0, x, y-tdy, 1, "Ch %d", customparam_ch_idx/2);
			//GUIPrint(0, x, y        , guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[12], sel[12], params[ 0], sel[13], sel[13], params[ 1], sel[14], sel[14], params[ 2], sel[15], sel[15], params[ 3], sel[16], sel[16], params[ 4]);
			//GUIPrint(0, x, y+ystep  , guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[17], sel[17], params[ 5], sel[18], sel[18], params[ 6], sel[19], sel[19], params[ 7], sel[20], sel[20], params[ 8], sel[21], sel[21], params[ 9]);
			//set_bk_color(bk0);
			//GUIPrint(0, x, y+ystep*2, guizoom, "%c%c%8.3lf %c%c%8.3lf",                                  sel[22], sel[22], params[10], sel[23], sel[23], params[11]);//do not change these strings!
			//set_text_colors(prevcolor);
		}
#if 0
		//clamp bounds
		//0123456789012345678901
		//[ SNNNN, SNNNN] clamp
		x=buttons[3].x1;
		y=buttons[3].y1;
		GUIPrint(0, x, y, guizoom, "[ %5d, %5d] clamp", customparam_clamp[0], customparam_clamp[1]);

		//learning rate
		//0123456789012345678901
		//lr -0.NNNNNNNNNNNNNNN
		x=buttons[4].x1;
		y=buttons[4].y1;
		GUIPrint(0, x, y, guizoom, "lr %18.15lf", g_lr);
#endif
	}
	else if(im1&&(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4]))
	{
		int c0=set_bk_color(0x80FFFFFF);
		float ystep=tdy*guizoom, x, y;
		x=buttons[5].x1;
		y=buttons[5].y1;
		if(ols4_cache)
			GUIPrint(0, x, y-ystep, guizoom, "M %d%d%d%d%d%d%d%d",
				ols4_cache>>7&1,
				ols4_cache>>6&1,
				ols4_cache>>5&1,
				ols4_cache>>4&1,
				ols4_cache>>3&1,
				ols4_cache>>2&1,
				ols4_cache>>1&1,
				ols4_cache>>0&1
			);
		GUIPrint(0, x, y, guizoom, "%8d%9.6lf%9.6lf%9.6lf%9.6lf", ols4_period, ols4_lr[0], ols4_lr[1], ols4_lr[2], ols4_lr[3]);
		for(int kc=0;kc<im1->nch;++kc)
		{
			for(int ky=0;ky<=OLS4_RMAX;++ky)
			{
				for(int kx=0;kx<(OLS4_RMAX<<1|1);++kx)
				{
					int val=ols4_mask[kc][(OLS4_RMAX<<1|1)*ky+kx];
					if(ky==OLS4_RMAX&&kx==OLS4_RMAX)
					{
						switch(kc)
						{
						case 0:
							GUIPrint_append(0, 0, 0, 0, 0, "   ?   ?");
							break;
						case 1:
							GUIPrint_append(0, 0, 0, 0, 0, "  ?%d  ?%d",
								val>>4&1,
								val>>0&1
							);
							break;
						case 2:
							GUIPrint_append(0, 0, 0, 0, 0, " ?%d%d ?%d%d",
								val>>5&1, val>>4&1,
								val>>1&1, val>>0&1
							);
							break;
						case 3:
							GUIPrint_append(0, 0, 0, 0, 0, " %d%d%d?%d%d%d",
								val>>6&1, val>>5&1, val>>4&1,
								val>>2&1, val>>1&1, val>>0&1
							);
							break;
						}
						//for(int k=kc-1;k>=0;--k)
						//	GUIPrint_append(0, 0, 0, 0, 0, "%d", val>>k&1);
						break;
					}
					GUIPrint_append(0, 0, 0, 0, 0, "%d%d%d%d%d%d%d%d ",
						val>>7&1,
						val>>6&1,
						val>>5&1,
						val>>4&1,
						val>>3&1,
						val>>2&1,
						val>>1&1,
						val>>0&1
					);
				}
				GUIPrint_append(0, x, y+ystep*((OLS4_RMAX+1)*kc+ky+1), guizoom, 1, 0);
			}
		}
		set_bk_color(c0);
	}
#if 0
	else if(transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC])
	{
		float x, y, info[4]={0};
		int kc=customparam_ch_idx/2;
		const short *params=logic_params+LOGIC_PARAMS_PER_CH*kc;

		x=0;
		y=(float)(h>>1);
		logic_opt_checkonthread(info);
		if(info[0]>=0)
		{
			float cr=1/info[2];
			TimeInfo ti;
			parsetimedelta(info[1]*1000, &ti);
			set_window_title(			"Ch %d: %6.2f%%, %02d-%02d-%f, CR %f, bitidx %d", kc, info[0], ti.hours, ti.mins, ti.secs, cr, ((int*)info)[3]);
			GUIPrint(0, x, y-tdy, 1,	"Ch %d: %6.2f%%, %02d-%02d-%f, CR %f, bitidx %d", kc, info[0], ti.hours, ti.mins, ti.secs, cr, ((int*)info)[3]);
			g_repaint=1;
			//int success=RedrawWindow(ghWnd, 0, 0, RDW_INTERNALPAINT|RDW_ERASENOW);
			//if(!success)
			//	LOG_ERROR("Failed to redraw");
			//swapbuffers();
			//return;
			//InvalidateRect(ghWnd, 0, 0);
			if(info[2]!=0&&ch_cr[kc]<cr)
				ch_cr[kc]=cr;
		}
		else
			GUIPrint(0, x, y-tdy, 1, "Ch %d", kc);
		//LOG_ERROR("Unreachable");

		int waitstatus=-1;
		if(info[0]>=0&&ghMutex)
		{
			waitstatus=WaitForSingleObject(ghMutex, INFINITE);
			//if(waitstatus!=WAIT_OBJECT_0)
			//	LOG_ERROR("WaitForSingleObject error");
		}
		for(int ky=0;ky<LOGIC_NF0;++ky)
		{
			print_i16_row(x, y, 1, params, LOGIC_ROWPARAMS);
			y+=tdy;
			params+=LOGIC_ROWPARAMS;
		}
#ifdef LOGIC_NF1
		print_i16_row(x, y, 1, params, LOGIC_NF1);
#else
		print_i16_row(x, y, 1, params, LOGIC_NF0);
#endif
		if(info[0]>=0&&ghMutex)
		{
			if(waitstatus==WAIT_OBJECT_0)
				ReleaseMutex(ghMutex);
		}
	}
#endif
	else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
	{
		float x, y;

		x=buttons[3].x1;
		y=buttons[3].y1;

		print_i16_row(x, y, 1, jxlparams_i16   , 11);	y+=tdy;
		print_i16_row(x, y, 1, jxlparams_i16+11, 11);	y+=tdy;
		print_i16_row(x, y, 1, jxlparams_i16+22, 11);
	}
	else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
	{
		float x, y;
		
		x=0;
		y=buttons[3].y1;
		print_i16_row(x, y, 1, pw2_params             , PW2_NPARAM);	y+=tdy;
		print_i16_row(x, y, 1, pw2_params+PW2_NPARAM  , PW2_NPARAM);	y+=tdy;
		print_i16_row(x, y, 1, pw2_params+PW2_NPARAM*2, PW2_NPARAM);
	}
#if 0
	else if(transforms_mask[ST_FWD_JOINT]||transforms_mask[ST_INV_JOINT])
	{
		int printed;
		float x, y;

		x=tdx;
		y=(float)((h>>1)+(h>>2));
		for(int kc=0;kc<3;++kc)
		{
			printed=0;
			for(int k=0;k<24;++k)
			{
				int val=jointpredparams[24*kc+k];
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, " %c%08X", val<0?'-':' ', abs(val));
			}
			print_line(0, x, y, 0.75f, g_buf, printed, -1, 0, 0);
			y+=tdy*0.75f;
		}
		for(int kc=0;kc<3;++kc)
		{
			printed=0;
			for(int k=0;k<4;++k)
			{
				int val=jointpredparams[72+4*kc+k];
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, " %c%08X", val<0?'-':' ', abs(val));
			}
			print_line(0, x, y, 0.75f, g_buf, printed, -1, 0, 0);
			y+=tdy*0.75f;
		}
	}
#endif
	const char *mode_str=0;
	switch(mode)
	{
//	case VIS_PLANES:		mode_str="Planes";		break;
//	case VIS_MESH:			mode_str="Combined Mesh";	break;
//	case VIS_MESH_SEPARATE:		mode_str="Separate Mesh";	break;
	case VIS_HISTOGRAM:		mode_str="Histogram";		break;
	case VIS_JOINT_HISTOGRAM:	mode_str="Joint Histogram";	break;
	case VIS_IMAGE:			mode_str="Image View";		break;
//	case VIS_BAYES:			mode_str="Bayes";		break;
	case VIS_ZIPF:			mode_str="Zipf View";		break;
//	case VIS_IMAGE_BLOCK:		mode_str="Image Block";		break;
//	case VIS_IMAGE_E24:		mode_str="Image Exp24";		break;
//	case VIS_DWT_BLOCK:		mode_str="DWT Block";		break;
	case VIS_IMAGE_TRICOLOR:	mode_str="Tricolor";		break;
	}
	if(im1)
	{
		float x=(float)(w-300), y=tdy*2, x2=x;
		for(int k=0;k<T_COUNT;++k)//print available transforms on right
		{
			transforms_printname(x2, y, k, -1, transforms_mask[k]?0xA0FF0000FFFFFFFF:0);
			x2=x+(150&-!(k&1));
			if(k&1)
				y+=tdy;
		}
		x=(float)(w-450);
		y=tdy*2;
		const char *label=ec_method_label(ec_method);
		if(ec_method==ECTX_ABAC)
			GUIPrint(x, x, y-tdy, 1, "H - - -  -         %s", label);
		else if(ec_adaptive)//H.E.M.L..A.0x0000..XXXX_XXX
			GUIPrint(x, x, y-tdy, 1, "H %d %d %d  A 0x%04X  %s", ec_expbits, ec_msb, ec_lsb, ec_adaptive_threshold, label);
		else
			GUIPrint(x, x, y-tdy, 1, "H %d %d %d  Static    %s", ec_expbits, ec_msb, ec_lsb, label);
		if(transforms)
		{
			for(int k=0;k<(int)transforms->count;++k, y+=tdy)//print applied transforms on left
				transforms_printname(x, y, transforms->data[k], k, 0);
		}
		float
			cr_combined=(float)((im1->src_depth[0]+im1->src_depth[1]+im1->src_depth[2]+im1->src_depth[3])/(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3])),
			xstart=20, xend=(float)w-330, ystart=(float)(h-tdy*5),
			scale=xend-xstart;
		
		double usize=image_getBMPsize(im0);
		float crformat=(float)(usize/filesize);
		int RGBspace=1;
		for(int k=0;k<CST_FWD_SEPARATOR;++k)
		{
			if(transforms_mask[k])
			{
				RGBspace=0;
				break;
			}
		}
		if(xstart<xend)
		{
			float maxinvcr=1/cr_combined;
			for(int k=0;k<4;++k)//get max CR
			{
				double invCR=ch_entropy[k]/im1->src_depth[k];
				if(isfinite(invCR)&&maxinvcr<invCR)
					maxinvcr=(float)invCR;
			}
			if(!isfinite(maxinvcr))
				maxinvcr=0;
			if(maxinvcr<1/crformat)
				maxinvcr=1/crformat;
			if(xend-scale*maxinvcr<xstart)
				scale=(xend-xstart)/maxinvcr;
			if(scale<50)
				scale=50;

			xstart=xend-scale*ceilf((maxinvcr+1)*2)*0.5f;
			draw_rect(xstart, xend, ystart, (float)h, 0x80808080);//background
			int ks;
			if(scale>20)
			{
				//double gain=scale/usize;
				//for(double size=100000;;)
				//{
				//	x=(float)(xend-size*gain);
				//	if(x<xstart)
				//		break;
				//	draw_line(x, ystart+1, x, (float)h, 0x70900090);//draw minor scale
				//	size+=100000;
				//}
				ks=1;
				x=(float)(xend-0.1f*scale*ks);
				for(;x>xstart;++ks)
				{
					draw_line(x, ystart+1, x, (float)h, 0x70900090);//draw minor scale
					x=(float)(xend-0.1f*scale*ks);
				}
			}
			ks=1;
			x=(float)(xend-scale*ks);
			for(int ks2=1;x>xstart;++ks2)
			{
				draw_rect(x-1, x+2, ystart, (float)h, 0x70800080);//draw major scale
				x=(float)(xend-scale*ks2);
			}
			float barw=4;
#if 0
			if(
				mode==VIS_IMAGE_BLOCK||
			//	mode==VIS_IMAGE_E24||
				mode==VIS_DWT_BLOCK
			)
			{
				float cr[4];
				//if(mode==VIS_IMAGE_E24)
				//{
				//	cr[0]=(float)e24_cr[0];
				//	cr[1]=(float)e24_cr[1];
				//	cr[2]=(float)e24_cr[2];
				//	cr[3]=3/(1/cr[0]+1/cr[1]+1/cr[2]);
				//}
				//else
				{
					cr[0]=blockCR[0];
					cr[1]=blockCR[1];
					cr[2]=blockCR[2];
					cr[3]=3/(1/blockCR[0]+1/blockCR[1]+1/blockCR[2]);
				}
				draw_rect(xend-scale*cr[3]      , xend, ystart+tdy*0.5f-barw, ystart+tdy*0.5f  , 0xFF404040);
				draw_rect(xend-scale*cr[0]      , xend, ystart+tdy*1.5f-barw, ystart+tdy*1.5f  , 0xFF0000B0);
				draw_rect(xend-scale*cr[1]      , xend, ystart+tdy*2.5f-barw, ystart+tdy*2.5f  , 0xFF00B000);
				draw_rect(xend-scale*cr[2]      , xend, ystart+tdy*3.5f-barw, ystart+tdy*3.5f  , 0xFFB00000);
				draw_rect(xend-scale*cr_combined, xend, ystart+tdy*0.5f, ystart+tdy*0.5f+barw+1, 0xFF000000);
				draw_rect(xend-scale*ch_cr[0]   , xend, ystart+tdy*1.5f, ystart+tdy*1.5f+barw+1, 0xFF0000FF);
				draw_rect(xend-scale*ch_cr[1]   , xend, ystart+tdy*2.5f, ystart+tdy*2.5f+barw+1, 0xFF00FF00);
				draw_rect(xend-scale*ch_cr[2]   , xend, ystart+tdy*3.5f, ystart+tdy*3.5f+barw+1, 0xFFFF0000);
			}
			else
#endif
			draw_rect(xend-scale/cr_combined,				xend, ystart+tdy*0.5f-barw, ystart+tdy*0.5f+barw+1, 0xFF000000);
			draw_rect(xend-scale/(float)(im1->src_depth[0]/ch_entropy[0]),	xend, ystart+tdy*1.5f-barw, ystart+tdy*1.5f+barw+1, RGBspace?0xFF0000FF:0xFF404040);//r or Y
			draw_rect(xend-scale/(float)(im1->src_depth[1]/ch_entropy[1]),	xend, ystart+tdy*2.5f-barw, ystart+tdy*2.5f+barw+1, RGBspace?0xFF00FF00:0xFFC00000);//g or Cb
			draw_rect(xend-scale/(float)(im1->src_depth[2]/ch_entropy[2]),	xend, ystart+tdy*3.5f-barw, ystart+tdy*3.5f+barw+1, RGBspace?0xFFFF0000:0xFF0000C0);//b or Cr
			if(im1->depth[3])
				draw_rect(xend-scale/(float)(im1->src_depth[1]/ch_entropy[3]), xend, ystart+tdy*4.5f-barw, ystart+tdy*4.5f+barw+1, RGBspace?0xFF808080:0xFF00C000);//a or Cg
			//draw_rect(xend-scale*ch_cr[3]   , xend, ystart+tdy*4.5f-barw, ystart+tdy*4.5f+barw+1, 0xFFFF00FF);
			x=xend-scale/crformat;
			draw_line(x, ystart, x-10, ystart-10, 0xFF000000);
			draw_line(x, ystart, x+10, ystart-10, 0xFF000000);
		}
		int prevtxtcolor, prevbkcolor;
		xend+=10;
		prevbkcolor=set_bk_color(0xC0C0C0C0);
		prevtxtcolor=set_text_color(0xFF000000);
		GUIPrint(xend, xend, ystart-tdy*2, 1, "Bitmap Size             %9.0lf", usize);
		//GUIPrint(xend, xend, ystart-tdy*2, 1, "Uncompressed Size       %9.0lf", usize);
		GUIPrint(xend, xend, ystart-tdy  , 1, "Format        %8.4f%% %9lld", 100./crformat, filesize);
		set_bk_color(0xE0FFFFFF);
		prevtxtcolor=set_text_color(0xFF000000);GUIPrint(xend, xend, ystart      , 1, "Combined      %8.4f%% %12.2lf", 100./cr_combined, usize/cr_combined);
		set_bk_color(0xC0C0C0C0);
		set_text_color(RGBspace?0xFF0000FF:0xFF404040);	GUIPrint(xend, xend, ystart+tdy  , 1, "%c     %7d %8.4lf%% %12.2lf", RGBspace?'R':'Y', im1->depth[0], 100.*ch_entropy[0]/im1->src_depth[0], (double)im1->iw*im1->ih*ch_entropy[0]/8);
		set_text_color(RGBspace?0xFF00C000:0xFFC00000);	GUIPrint(xend, xend, ystart+tdy*2, 1, "%c     %7d %8.4lf%% %12.2lf", RGBspace?'G':'U', im1->depth[1], 100.*ch_entropy[1]/im1->src_depth[1], (double)im1->iw*im1->ih*ch_entropy[1]/8);
		set_text_color(RGBspace?0xFFFF0000:0xFF0000C0);	GUIPrint(xend, xend, ystart+tdy*3, 1, "%c     %7d %8.4lf%% %12.2lf", RGBspace?'B':'V', im1->depth[2], 100.*ch_entropy[2]/im1->src_depth[2], (double)im1->iw*im1->ih*ch_entropy[2]/8);
		if(im1->depth[3])
		{
			set_text_color(RGBspace?0xFF404040:0xFF00C000);	GUIPrint(xend, xend, ystart+tdy*4, 1, "%c     %7d %8.4lf%% %12.2lf", RGBspace?'A':'W', im1->depth[3], 100.*ch_entropy[3]/im1->src_depth[1], (double)im1->iw*im1->ih*ch_entropy[3]/8);//src_depth[3] is zero assuming no alpha
		}
		//set_text_color(0xFFFF00FF);	GUIPrint(xend, xend, ystart+tdy*4, 1, "Joint %7d %9f", usage[3], ch_cr[3]);
		set_text_color(prevtxtcolor);
		set_bk_color(prevbkcolor);

		//if(transforms_customenabled)
		//{
		//	//double maxloss=0;
		//	GUIPrint(0, 0, tdy*3, 1, "RMSE %lf", av_rmse);
		//	if(minloss<maxloss)
		//		GUIPrint(0, 200, tdy*3, 1, "[%lf~%lf]", minloss, maxloss);
		//}

		float g2=h/combCRhist_max;
		int idx=combCRhist_idx-1, idx2=combCRhist_idx-2;
		idx+=combCRhist_SIZE&-(idx<0);
		idx2+=combCRhist_SIZE&-(idx2<0);
		xstart=(float)(w-(combCRhist_SIZE<<combCRhist_logDX)-300);
		float cx=xstart+(float)(idx<<combCRhist_logDX), cy=h-combCRhist[idx][3]*g2;
		draw_rect_hollow(cx-10, cx+10, cy-10, cy+10, 0xC0C0C0C0);
		if(combCRhist[idx][3]<combCRhist[idx2][3])//ratio improved
		{
			draw_triangle(cx-9, cy, cx+9, cy, cx, cy+9, 0xC080FF80);
			draw_triangle(cx-5, cy, cx+5, cy, cx, cy+5, 0xC0006000);
		}
		else if(combCRhist[idx][3]>combCRhist[idx2][3])//ratio worsened
		{
			draw_triangle(cx-9, cy, cx+9, cy, cx, cy-9, 0xC08080FF);
			draw_triangle(cx-5, cy, cx+5, cy, cx, cy-5, 0xC0000060);
		}
		else//exact same ratio
		{
			draw_line(cx-5, cy-5, cx+5, cy+5, 0xC0808080);
			draw_line(cx-5, cy+5, cx+5, cy-5, 0xC0808080);
		}
		for(int k=0;k<combCRhist_SIZE-1;++k)
		{
			if(k!=combCRhist_idx-1)
			{
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][0]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][0]*g2), RGBspace?0xFF0000FF:0xFF404040);//r or Y
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][1]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][1]*g2), RGBspace?0xFF00FF00:0xFFC00000);//g or Cb
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][2]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][2]*g2), RGBspace?0xFFFF0000:0xFF0000C0);//b or Cr
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][3]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][3]*g2), 0xC0000000);
			}
		}
		draw_line(xstart, h-g2, xstart+(combCRhist_SIZE<<combCRhist_logDX), h-g2, 0xC0000000);

		//grad2 info
#if 1
#if 0
		if(transforms_mask[ST_PREPROC_GRAD]||transforms_mask[ST_PREPROC_X2])
		{
			extern double lossygrad_rmse[3], lossygrad_psnr[3];
			x=(float)(w>>1);
			y=(float)(h>>1);
			GUIPrint(0, x, y, 1, "RMSE RGB %14lf %14lf %14lf", lossygrad_rmse[0], lossygrad_rmse[1], lossygrad_rmse[2]); y+=tdy;
			GUIPrint(0, x, y, 1, "PSNR RGB %14lf %14lf %14lf", lossygrad_psnr[0], lossygrad_psnr[1], lossygrad_psnr[2]);
		}
		else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
		{
			float width=300.f, zoom=1, ystart=tdy*zoom*3;
			x=(float)(w>>2);
			//x=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16_row(x, y, zoom, c2_params.c0+k*12, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c0+k*12+5, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c0+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}

			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16_row(x, y, zoom, c2_params.c1+k*12, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c1+k*12+5, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c1+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}
			print_i16_row(x, y, zoom, c2_params.c1+6*12, 2);		y+=tdy*zoom;
			
			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16_row(x, y, zoom, c2_params.c2+k*12, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c2+k*12+5, 5);		y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c2+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}
			print_i16_row(x, y, zoom, c2_params.c2+6*12, 4);		y+=tdy*zoom;
		}
		else
#endif
		if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
		{
			const int width=C3_REACH<<2|2;
			const short *coeffs[]=
			{
				c3_params.c00, c3_params.c01, c3_params.c02,
				c3_params.c10, c3_params.c11, c3_params.c12,
				c3_params.c20, c3_params.c21, c3_params.c22,
			};
			const int tails[]=
			{
				 C3_REACH<<1,     C3_REACH<<1,    C3_REACH<<1,
				(C3_REACH<<1)+2,  C3_REACH<<1,    C3_REACH<<1,
				(C3_REACH<<1)+2, (C3_REACH<<1)+2, C3_REACH<<1,
			};
			float zoom=0.75;
			x=(float)(w>>1);
			y=tdy*zoom*4;
			for(int kc=0;kc<9;++kc)
			{
				for(int k=0;k<C3_REACH;++k)
				{
					print_i16_row(x, y, zoom, coeffs[kc]+width*k, width);
					y+=tdy*zoom;
				}
				print_i16_row(x, y, zoom, coeffs[kc]+width*C3_REACH, tails[kc]);
				y+=tdy*zoom*2;
			}
		}
		else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
		{
			x=(float)(w>>1);
			y=0;
			for(int k=0;k<PW2_NPRED;++k)
			{
				GUIPrint(0, x, y, 0.9f, "p%d %10lf", k, pw2_errors[k]);
				y+=tdy*0.9f;
			}
		}
#if 0
		else if(transforms_mask[ST_FWD_GRAD2]||transforms_mask[ST_INV_GRAD2])
		{
			const char *prednames[]=
			{
				"grad    ",
				"grad45  ",
				"path    ",
				"path45  ",
				"gamma   ",
				"select  ",
				"grad2   ",
				"combined",
			};
			x=(float)(w>>1);
			y=(float)(h>>1);
			for(int k=0;k<GRAD2PREDCOUNT;++k)
			{
				GUIPrint(0, x, y, 1, "%d%9lf %7d %s", k, grad2_csize[k], grad2_hits[k], prednames[k]);
				y+=tdy;
			}
			//GUIPrint(0, 0, (float)(h>>1)-tdy, 1, "RMSE grad%10lf:%d grad45%10lf:%d path%10lf:%d path45%10lf:%d gamma%10lf:%d",
			//	grad2_rmse[0], grad2_hits[0],
			//	grad2_rmse[1], grad2_hits[1],
			//	grad2_rmse[2], grad2_hits[2],
			//	grad2_rmse[3], grad2_hits[3],
			//	grad2_rmse[4], grad2_hits[4]);
		}
		else if(transforms_mask[ST_FWD_ADAPTIVE]||transforms_mask[ST_INV_ADAPTIVE])
		{
			const char *prednames[]=
			{
#if 1
				"grad    ",
				"avgall  ",
				"left    ",
				"top     ",
				"topleft ",
				"topright",
				"linx    ",
				"liny    ",
#endif
#if 0
				"hole        ",
				"bottom-right",
				"bottom-left ",
				"roof bottom ",
				"valley top  ",
				"top-right   ",
				"top-left    ",
				"peak        ",
#endif
			};

			//draw histograms
#if 1
			float
				xmargin=(float)w*0.05f, dx=(w-xmargin*2)/4,
				ymargin=(float)h*0.05f, dy=(h-ymargin*2)/2;
			for(int ky=0;ky<2;++ky)
			{
				for(int kx=0;kx<4;++kx)
				{
					float
						x1=xmargin+dx*kx, x2=xmargin+dx*(kx+1)*0.95f,
						y1=ymargin+dy*ky, y2=ymargin+dy*(ky+1)*0.95f;
					int idx=ky<<2|kx;
					chart_hist_draw2(x1, x2, y1, y2, 0x80808080, grad2_hist+(idx<<8), -1);
					draw_rect_hollow(x1, x2, y1, y2, 0x80808080);
					GUIPrint(0, x1, y1, 1, "%d %s", idx, prednames[idx]);
				}
			}
#endif

			x=tdx*guizoom;
			y=(float)((h>>1)+(h>>2));
			//x=(float)(w>>2);
			//y=(float)(h>>1);

			int total=0;
			double csize=0;
			for(int k=0;k<ADAGRADCOUNT;++k)
			{
				total+=adagrad_hits[k];
				csize+=adagrad_csize[k];
			}

			g_printed=0;
			for(int k=0;k<ADAGRADCOUNT;++k)
			{
				float width;
#define PRINTSTRING "%d %s\t%7d %5.2lf%% %10lf %14f CR %5.3lf  E %11.8lf abs %12.9lf", k, prednames[k], adagrad_hits[k], 100.*adagrad_hits[k]/total, adagrad_rmse[k], adagrad_csize[k], adagrad_hits[k]/adagrad_csize[k], adagrad_signederror[k]/(iw*ih*3), adagrad_abserror[k]/(iw*ih*3)
				width=GUIPrint_append(0, x, y, 1, 1, PRINTSTRING);
				GUIPrint_append(0, x, y, 1, 0, "\n");

				//GUIPrint(0, x, y, 1, PRINTSTRING);
#undef  PRINTSTRING

				width=ceilf(width/50)*50;
				if(x+width<w-400)
				{
					float divisor=x+width+(float)(adagrad_csize[k]*(w-400-(x+width))/csize);
					draw_rect(x+width, divisor, y, y+tdy*0.75f, 0x80808080);
					draw_rect(divisor, x+width+(float)(adagrad_hits[k]*(w-400-(x+width))/csize), y, y+tdy*0.75f, 0x80C0C0C0);
				}
				y+=tdy;
			}
			static int copied=0;
			if(!copied)
			{
				copy_to_clipboard(g_buf, g_printed);
				copied=1;
			}
		}
		else if(transforms_mask[ST_FWD_SORTNB]||transforms_mask[ST_INV_SORTNB])
		{
			const char *casenames[]=
			{
				"A",
				"(A+B)/2",
				"B",
				"(B+C)/2",
				"C",
				"(C+D)/2",
				"D",
				"grad",
			};
			int res=iw*ih*3;
			x=(float)(w>>1);
			y=(float)(h>>1);
			g_printed=0;
			for(int k=0;k<SORTNBCASES;++k)
			{
				GUIPrint_append(0, x, y, 1, 1, "%5.2lf%%  E %5.2f %s", 100.*sortnb_cases[k]/res, sortnb_rmse[k], casenames[k]);
				GUIPrint_append(0, x, y, 1, 0, "\n");
				y+=tdy;
			}
			static int copied2=0;
			if(!copied2)
			{
				copy_to_clipboard(g_buf, g_printed);
				copied2=1;
			}
		}
#endif
#endif

		//extern int testhist[3];//
		//GUIPrint(0, 0, tdy*2, 1, "%d %d %d", testhist[0], testhist[1], testhist[2]);//
		//const char *label=ec_method_label(ec_method);
		GUIPrint(0, 0, 0, 1, "WH %dx%d  D0[%d %d %d %d] D[%d %d %d %d]",
			im0->iw, im0->ih,
			im0->src_depth[0], im0->src_depth[1], im0->src_depth[2], im0->src_depth[3],
			im1->depth[0], im1->depth[1], im1->depth[2], im1->depth[3]
		);
	}
#if 0
	if(image)
	{
		float cr_combined=3/(1/ch_cr[0]+1/ch_cr[1]+1/ch_cr[2]), scale=150;
		int xstart=w>>1;
		if(xstart<200)//text width
			xstart=200;
		float vmax=cr_combined;
		for(int k=0;k<4;++k)
		{
			if(isfinite((double)ch_cr[k])&&vmax<ch_cr[k])
				vmax=ch_cr[k];
		}
		if(xstart+scale*vmax>w-20)
			scale=(w-20-xstart)/vmax;
		if(scale<5)
			scale=5;
		float ystart=(float)(h-tdy*5);
		draw_rect((float)xstart, (float)w, ystart, (float)h, 0x80808080);
		for(int ks=1;scale*ks<w;++ks)//draw scale
		{
			float x=(float)(xstart+scale*ks);
			draw_rect(x-1, x+1, ystart, (float)h, 0x80800080);
		}
		float barw=3;
		draw_rect((float)xstart, xstart+scale*ch_cr[0]   , ystart+tdy*0.5f-barw, ystart+tdy*0.5f+barw, 0xFF0000FF);
		draw_rect((float)xstart, xstart+scale*ch_cr[1]   , ystart+tdy*1.5f-barw, ystart+tdy*1.5f+barw, 0xFF00FF00);
		draw_rect((float)xstart, xstart+scale*ch_cr[2]   , ystart+tdy*2.5f-barw, ystart+tdy*2.5f+barw, 0xFFFF0000);
		draw_rect((float)xstart, xstart+scale*cr_combined, ystart+tdy*3.5f-barw, ystart+tdy*3.5f+barw, 0xFF000000);
		draw_rect((float)xstart, xstart+scale*ch_cr[3]   , ystart+tdy*4.5f-barw, ystart+tdy*4.5f+barw, 0xFFFF00FF);
		int prevcolor=set_text_color(0xFF0000FF);
		//int prevbk=set_bk_color(0xA0808080);
		GUIPrint(0, 0, ystart      , 1, "R     %7d %9f", usage[0], ch_cr[0]), set_text_color(0xFF00FF00);
		GUIPrint(0, 0, ystart+tdy  , 1, "G     %7d %9f", usage[1], ch_cr[1]), set_text_color(0xFFFF0000);
		GUIPrint(0, 0, ystart+tdy*2, 1, "B     %7d %9f", usage[2], ch_cr[2]), set_text_color(0xFF000000);
		GUIPrint(0, 0, ystart+tdy*3, 1, "Combined      %9f", cr_combined   ), set_text_color(0xFFFF00FF);
		GUIPrint(0, 0, ystart+tdy*4, 1, "Joint %7d %9f", usage[3], ch_cr[3]), set_text_color(prevcolor);
		//set_bk_color(prevbk);
	}
#endif
#if 0
	if(image)
	{
		float cr_combined=3/(1/ch_cr[0]+1/ch_cr[1]+1/ch_cr[2]);
		int scale=150;
		draw_rect(0, tdx*6.5f, 0, (float)h, 0x80808080);
		for(int ks=0;scale*ks<h;++ks)
		{
			float y=(float)(h-scale*ks);
			draw_rect(0, tdx*6.5f, y-1, y+1, 0xFF800080);
			//for(int kt=0;kt<2;++kt)
			//{
			//	float y=(float)(h-scale*ks+kt);
			//	draw_line(0, y, tdx*6.5f, y, 0xFF800080);
			//}
		}
		draw_rect(tdx  -2, tdx  +2, h-scale*ch_cr[0]   , (float)h, 0xFF0000FF);
		draw_rect(tdx*2-2, tdx*2+2, h-scale*ch_cr[1]   , (float)h, 0xFF00FF00);
		draw_rect(tdx*3-2, tdx*3+2, h-scale*ch_cr[2]   , (float)h, 0xFFFF0000);
		draw_rect(tdx*4-2, tdx*4+2, h-scale*cr_combined, (float)h, 0xFF000000);
		draw_rect(tdx*5-2, tdx*5+2, h-scale*ch_cr[3]   , (float)h, 0xFFFF00FF);
		//for(int k=0;k<4;++k)
		//{
		//	draw_line(tdx  +k, h-scale*ch_cr[0]   , tdx  +k, (float)h, 0xFF0000FF);
		//	draw_line(tdx*2+k, h-scale*ch_cr[1]   , tdx*2+k, (float)h, 0xFF00FF00);
		//	draw_line(tdx*3+k, h-scale*ch_cr[2]   , tdx*3+k, (float)h, 0xFFFF0000);
		//	draw_line(tdx*4+k, h-scale*cr_combined, tdx*4+k, (float)h, 0xFF000000);
		//	draw_line(tdx*5+k, h-scale*ch_cr[3]   , tdx*5+k, (float)h, 0xFFFF00FF);
		//}
		GUIPrint(0, 0, tdy*2, 1, "Usage:CR RGB[%d:%f, %d:%f, %d:%f] %f, joint %d:%f", usage[0], ch_cr[0], usage[1], ch_cr[1], usage[2], ch_cr[2], cr_combined, usage[3], ch_cr[3]);
	}
#endif
	//extern double bestslope;
	//GUIPrint(0, 0, 0, 1, "p(%f, %f, %f) dx %f a(%f, %f) fov %f, bestslope=%lf", cam.x, cam.y, cam.z, cam.move_speed, cam.ax, cam.ay, atan(cam.tanfov)*180/M_PI*2, bestslope);
	//if(im0&&im1)
	//{
	//}
	//GUIPrint(0, 0, 0, 1, "p(%f, %f, %f) dx %f a(%f, %f) fov %f", cam.x, cam.y, cam.z, cam.move_speed, cam.ax, cam.ay, atan(cam.tanfov)*180/M_PI*2);
	
	static double t=0;
	double t2=time_ms();
	GUIPrint(0, 0, tdy, 1, "timer %d, fps%11lf, [%2d/%2d] %s", timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str);
	if(mode==VIS_IMAGE||mode==VIS_ZIPF)
		GUIPrint(0, 0, tdy*2, 1, "%s", show_full_image?"FILL SCREEN":"1:1");
	else if(mode==VIS_JOINT_HISTOGRAM)
	{
		switch(space_not_color)
		{
		case 0:mode_str="COLOR   (R, G, B)";break;
		case 1:mode_str="SPACE X (CURR, W, WW)";break;
		case 2:mode_str="SPACE Y (CURR, N, NN)";break;
		case 3:mode_str="SPACE   (CURR, N, W)";break;
		}
		GUIPrint(0, 0, tdy*2, 1, "[Shift C] %s  [Ctrl Wheel] box %dx%d  [Wheel] contour %f", mode_str, jhc_boxdx, jhc_boxdy, jhc_level);
	}
	t=t2;

	swapbuffers();
}
int io_quit_request(void)//return 1 to exit
{
	//logic_opt_forceclosethread();
	//g_repaint=0;
	//int button_idx=messagebox(MBOX_OKCANCEL, "Are you sure?", "Quit application?");
	//return button_idx==0;

	return 1;
}
void io_cleanup(void)//cleanup
{
	free(im0);
	free(im1);
	array_free(&cpu_vertices);
}