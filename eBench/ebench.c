#include"ebench.h"
#include<stdio.h>//snprintf
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
//#include"stb_image.h"
#include"lodepng.h"
#define DEBUG_MEMORY_IMPLEMENTATION
#include"intercept_malloc.h"
#include"c18.h"
static const char file[]=__FILE__;

	#define ENABLE_L1WEIGHTS
	#define USE_OLS

#ifdef ENABLE_L1WEIGHTS
	#define L1DELAY 50
	#define L1SPEED 20
#if 1
#define L1BIAS0	0
#define L1NPREDS 10
#define L1PREDLIST\
	L1PRED(100000, N)\
	L1PRED(100000, W)\
	L1PRED(100000, 3*(N-NN)+NNN)\
	L1PRED(100000, 3*(W-WW)+WWW)\
	L1PRED(100000, W+NE-N)\
	L1PRED(100000, N+W-NW)\
	L1PRED(100000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	L1PRED(100000, N+NE-NNE)\
	L1PRED(100000, W+NW-NWW)\
	L1PRED(100000, NEEE)
#endif
#if 0
#define L1BIAS0	0
#define L1NPREDS 13
#define L1PREDLIST\
	L1PRED(100000, (N+W-NW))\
	L1PRED(200000, N+W-NW)\
	L1PRED(100000, N)\
	L1PRED(100000, W)\
	L1PRED(100000, 3*(N-NN)+NNN)\
	L1PRED(100000, 3*(W-WW)+WWW)\
	L1PRED(100000, W+NE-N)\
	L1PRED(100000, N+NE-NNE)\
	L1PRED(100000, W+((NEEE+NEEEEE-N-W)>>3))\
	L1PRED( 50000, W+NW-NWW)\
	L1PRED( 50000, N+NW-NNW)\
	L1PRED( 50000, NE+NEE-NNEEE)\
	L1PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)
#endif
#if 0
#define L1BIAS0	0
#define L1NPREDS 12
#define L1PREDLIST\
	L1PRED( 90000, N)\
	L1PRED( 90000, W)\
	L1PRED( 90000, 3*(N-NN)+NNN)\
	L1PRED( 60000, 3*(W-WW)+WWW)\
	L1PRED(160000, N+W-NW)\
	L1PRED( 90000, W+NE-N)\
	L1PRED(100000, N+NE-NNE)\
	L1PRED( 20000, W+NW-NWW)\
	L1PRED( 50000, N+NW-NNW)\
	L1PRED( 70000, NE+NEE-NNEEE)\
	L1PRED( 80000, W+((NEEE+NEEEEE-N-W)>>3))\
	L1PRED( 90000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)
#endif
static const char *l1prednames[]=
{
#define L1PRED(W0, EXPR) #EXPR,
	L1PREDLIST
#undef  L1PRED
};
static int l1weights[L1NPREDS+1]=
{
#define L1PRED(W0, EXPR) W0,
	L1PREDLIST
#undef  L1PRED
	L1BIAS0,
};
#ifdef USE_OLS
static double curr_cov[L1NPREDS*L1NPREDS]={0}, curr_vec[L1NPREDS]={0}, curr_cholesky[L1NPREDS*L1NPREDS]={0}, curr_params[L1NPREDS]={0};
#endif
static int wperrors[L1NPREDS]={0};
static int use_ols=0;//0: L1    1: OLS    2: WP
#endif
float
	mouse_sensitivity=0.003f,
	key_turn_speed=0.03f;
Camera cam=
{
	10, 10, 10,
	(float)(225*M_PI/180), (float)(324.7356103172454f*M_PI/180),
	1,
	0.04f, (float)(2*M_PI/180),
	0, 0, 0, 0,
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
#ifdef ENABLE_L1WEIGHTS
	VIS_L1WEIGHTS,
#endif
	VIS_IMAGE,
	VIS_HISTOGRAM,
	VIS_MODEL,
	VIS_ANALYSIS,
	VIS_JOINT_HISTOGRAM,
	VIS_ZIPF,
	//VIS_BAYES,
	//VIS_IMAGE_BLOCK,
	//VIS_DWT_BLOCK,

	VIS_COUNT,
} VisMode;
int mode=VIS_IMAGE;

typedef enum _ImageViewMode
{
	VIEW_RGB,
	VIEW_C0,
	VIEW_C1,
	VIEW_C2,

	VIEW_COUNT,
} ImageViewMode;
int viewmode=VIEW_RGB, view_ma=0;

typedef enum TransformTypeEnum
{
	CT_FWD_SubG_OPT,	CT_INV_SubG_OPT,
	CT_FWD_SUBGREEN,	CT_INV_SUBGREEN,
	CT_FWD_JPEG2000,	CT_INV_JPEG2000,//	(1997) JPEG2000 RCT
	CT_FWD_JPEG2000_MA,	CT_INV_JPEG2000_MA,
	CT_FWD_NBLI,		CT_INV_NBLI,	//	(2024) NBLI
	CT_FWD_YCoCg_R,		CT_INV_YCoCg_R,	//	(2003) AVC, HEVC, VVC
	CT_FWD_YCbCr_R_v1,	CT_INV_YCbCr_R_v1,
	CT_FWD_YCbCr_R_v2,	CT_INV_YCbCr_R_v2,
	CT_FWD_YCbCr_R_v3,	CT_INV_YCbCr_R_v3,
	CT_FWD_YCbCr_R_v4,	CT_INV_YCbCr_R_v4,
	CT_FWD_YCbCr_R_v5,	CT_INV_YCbCr_R_v5,
	CT_FWD_YCbCr_R_v6,	CT_INV_YCbCr_R_v6,
	CT_FWD_YCbCr_R_v7,	CT_INV_YCbCr_R_v7,
	CT_FWD_Pei09,		CT_INV_Pei09,
//	CT_FWD_J2K2,		CT_INV_J2K2,
//	CT_FWD_CrCgCb,		CT_INV_CrCgCb,
//	CT_FWD_YRGB_v1,		CT_INV_YRGB_v1,
//	CT_FWD_YRGB_v2,		CT_INV_YRGB_v2,
//	CT_FWD_CMYK,		CT_INV_CMYK_DUMMY,
//	CT_FWD_MATRIX,		CT_INV_MATRIX,
	CT_FWD_YCbCr,		CT_INV_YCbCr,	//LOSSY	JPEG
	CT_FWD_XYB,		CT_INV_XYB,	//LOSSY	(2021) JPEG XL
	CT_FWD_CUSTOM,		CT_INV_CUSTOM,
//	CT_FWD_ADAPTIVE,	CT_INV_ADAPTIVE,

	CST_COMPARE,		CST_INV_SEPARATOR,
	
	ST_FWD_MIXN,		ST_INV_MIXN,
	ST_FWD_GRFILT,		ST_INV_GRFILT,
	ST_FWD_L1CRCT,		ST_INV_L1CRCT,
	ST_FWD_OLS7,		ST_INV_OLS7,
	ST_FWD_L1BCRCT,		ST_INV_L1BCRCT,
	ST_FWD_OLS8,		ST_INV_OLS8,
	ST_FWD_OLS9,		ST_INV_OLS9,
	ST_FWD_CGCRCT,		ST_INV_CGCRCT,
	ST_FWD_SUB,		ST_INV_SUB,
	ST_FWD_CLAMPGRAD,	ST_INV_CLAMPGRAD,
	ST_FWD_SELECT,		ST_INV_SELECT,
	ST_FWD_AV2,		ST_INV_AV2,
	ST_FWD_CALIC,		ST_INV_CALIC,
	ST_CONVTEST,		ST_CONVTEST2,
	ST_DIFF,		ST_SSIM,
	ST_FILT_MEDIAN33,	ST_FILT_AV33,
	ST_FILT_DEINT422,	ST_FILT_DEINT420,
	ST_FWD_QUANT,		ST_INV_QUANT,

	ST_FWD_PACKSIGN,	ST_INV_PACKSIGN,
	ST_FWD_PALETTE,		ST_INV_PALETTE,
	ST_FWD_BWTX,		ST_INV_BWTX,
	ST_FWD_BWTY,		ST_INV_BWTY,
	ST_FWD_MTF,		ST_INV_MTF,
	ST_FWD_CLEARTYPE,	ST_INV_CLEARTYPE,
	ST_FWD_AV4,		ST_INV_AV4,
	ST_FWD_SEL4,		ST_INV_SEL4,
	ST_FWD_CGPLUS,		ST_INV_CGPLUS,
	ST_FWD_CG3D,		ST_INV_CG3D,
	ST_FWD_PU,		ST_INV_PU,
	ST_FWD_CG422,		ST_INV_CG422,
	ST_FWD_CG420,		ST_INV_CG420,
	ST_FWD_LEGALLCG,	ST_INV_LEGALLCG,
	ST_FWD_WP,		ST_INV_WP,
	ST_FWD_WGRAD,		ST_INV_WGRAD,
//	ST_FWD_WGRAD2,		ST_INV_WGRAD2,
	ST_FWD_WGRAD3,		ST_INV_WGRAD3,
	ST_FWD_WGRAD4CCRCT,	ST_INV_WGRAD4CCRCT,
	ST_FWD_WGRAD4C,		ST_INV_WGRAD4C,
	ST_FWD_WGRAD4,		ST_INV_WGRAD4,
	ST_FWD_WGRAD5,		ST_INV_WGRAD5,
	ST_FWD_WGRAD6,		ST_INV_WGRAD6,
	ST_FWD_WGRAD7,		ST_INV_WGRAD7,
	ST_FWD_SSE,		ST_INV_SSE,
//	ST_FWD_WMIX,		ST_INV_WMIX,
	ST_FWD_T47,		ST_INV_T47,
	ST_FWD_P3,		ST_INV_P3,
	ST_FWD_G2,		ST_INV_G2,
	ST_FWD_MM,		ST_INV_MM,
//	ST_FWD_WP2,		ST_INV_WP2,
	ST_FWD_WPU,		ST_INV_WPU,
//	ST_FWD_DEFERRED,	ST_INV_DEFERRED,
	ST_FWD_WC,		ST_INV_WC,	//irreversible conv
	ST_FWD_CUSTOM4,		ST_INV_CUSTOM4,	//irreversible conv
	ST_FWD_CUSTOM3,		ST_INV_CUSTOM3,
	ST_FWD_CUSTOM,		ST_INV_CUSTOM,
	ST_FWD_CC,		ST_INV_CC,
	ST_FWD_NBLIC,		ST_INV_NBLIC,
	ST_FWD_OLS,		ST_INV_OLS,
	ST_FWD_OLS2,		ST_INV_OLS2,
	ST_FWD_OLS3,		ST_INV_OLS3,
	ST_FWD_OLS4,		ST_INV_OLS4,
	ST_FWD_OLS5,		ST_INV_OLS5,
	ST_FWD_OLS6,		ST_INV_OLS6,
	ST_FWD_TABLE,		ST_INV_TABLE,
	ST_FWD_LWAV,		ST_INV_LWAV,
	ST_FWD_MIX2,		ST_INV_MIX2,
//	ST_FWD_AV3,		ST_INV_AV3,
//	ST_FWD_ECOEFF,		ST_INV_ECOEFF,
//	ST_FWD_AVERAGE,		ST_INV_AVERAGE,
//	ST_FWD_MULTISTAGE,	ST_INV_MULTISTAGE,
	ST_FWD_ZIPPER,		ST_INV_ZIPPER,
//	ST_FWD_DIR,		ST_INV_DIR,
#if 0
	ST_FWD_CTX,		ST_INV_CTX,
	ST_FWD_C20,		ST_INV_C20,
//	ST_FWD_WU97,		ST_INV_WU97,
	ST_FWD_SHUFFLE,		ST_INV_SHUFFLE,
//	ST_FWD_SPLIT,		ST_INV_SPLIT,
#endif
	
	ST_FWD_LAZY,		ST_INV_LAZY,
	ST_FWD_HAAR,		ST_INV_HAAR,
	ST_FWD_SQUEEZE,		ST_INV_SQUEEZE,
	ST_FWD_LEGALL53,	ST_INV_LEGALL53,
	ST_FWD_CDF97,		ST_INV_CDF97,
//	ST_FWD_EXPDWT,		ST_INV_EXPDWT,
//	ST_FWD_CUSTOM_DWT,	ST_INV_CUSTOM_DWT,

	ST_FWD_DCT4,		ST_INV_DCT4,
	ST_FWD_DCT8,		ST_INV_DCT8,
	ST_FWD_WHT4,		ST_INV_WHT4,
	ST_FWD_WHT8,		ST_INV_WHT8,
	ST_FWD_HAAR8,		ST_INV_HAAR8,
	ST_FWD_FDCT8,		ST_INV_FDCT8,

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

C18Info analysis_info={0};
static long long pred_hist[3][PRED_COUNT]={0};

ArrayHandle jointhist=0;
int jointhist_nbits=6;//max
int jhx=0, jhy=0;
double ch_entropy[4]={0};//RGBA/YUVA
//int usage[4]={0};
EContext ec_method=ECTX_GRCTX;
//EContext ec_method=ECTX_HIST;//ECTX_MIN_QN_QW;
int ec_adaptive=0, ec_adaptive_threshold=3200, ec_expbits=5, ec_msb=2, ec_lsb=0;
int abacvis_low=0, abacvis_range=0;

#define MODELPREC 3
//#define MODELNCTX 256
#define MODELNCTX 17

//#define MODELPREC 3
//#define MODELNCTX (1<<(16+MODELPREC))

//#define MODELCTXBITS 6	//X
//#define MODELNCTX (1<<MODELCTXBITS)
int modelnch=0, modelnctx=0, modeldepth=0, modelhistsize=0, *modelhist=0;
double modelcsizes[4*MODELNCTX]={0};
double modelmeans[4*MODELNCTX]={0}, modelsdevs[4*MODELNCTX]={0};
int *modelmhist=0;
double modelmsizes[4*MODELNCTX]={0};
double modelstatoverhead=0;

#define combCRhist_SIZE 128
#define combCRhist_logDX 2
float combCRhist[combCRhist_SIZE][4]={0}, combCRhist_max=1;
int combCRhist_idx=0;

int imagecentered=0;
double
	imzoom=1,//image pixel size in screen pixels
	wpx=0, wpy=0;//window position (top-left corner) in image coordinates
#define screen2image_x(SX)             (wpx+(SX)/imzoom)
#define screen2image_y(SY)             (wpy+(SY)/imzoom)
#define screen2image_x_int(SX)         (int)floor(screen2image_x(SX))
#define screen2image_y_int(SY)         (int)floor(screen2image_y(SY))
#define screen2image_x_int_rounded(SX) (int)floor(screen2image_x(SX)+0.5)
#define screen2image_y_int_rounded(SY) (int)floor(screen2image_y(SY)+0.5)
#define image2screen_x(IX)             (((IX)-wpx)*imzoom)
#define image2screen_y(IY)             (((IY)-wpy)*imzoom)
#define image2screen_x_int(IX)         (int)floor(image2screen_x(IX))
#define image2screen_y_int(IY)         (int)floor(image2screen_y(IY))
typedef enum ProfilePlotModeEnum
{
	PROFILE_OFF,
	PROFILE_X,
	PROFILE_Y,
} ProfilePlotMode;
ProfilePlotMode profileplotmode=PROFILE_OFF;
int show_full_image=0;
int space_not_color=0;
#define ZOOM_LIMIT_LABEL 48
#define ZOOM_LIMIT_ALPHA 96
int pxlabels_hex=0;

//joint hist box contour
float jh_cubesize=64;
int jhc_boxdx=64, jhc_boxdy=64;
int jhc_xbox=0, jhc_ybox=0;//[0, screen_dim-1]
ArrayHandle jhc_mesh=0;
unsigned jhc_gpubuf=0;
float jhc_level=10.5f;
int loud_transforms=1;
int g_dist=1;
float g_uiscale=1;
//int crop_enable=0, crop_x=0, crop_y=0, crop_dx=128, crop_dy=128;

static void entropy2invcr(const double *entropy, const char *src_depth, int nch, double *invcr)//invcr: {T, R/Y, G/U, ...} nch+1=5 elements
{
	double etotal;
	int dtotal;

	etotal=0;
	dtotal=0;
	for(int k=0;k<nch;++k)
	{
		double e=entropy[k];
		int d=src_depth[k];
		etotal+=e;
		dtotal+=d;
		invcr[k+1]=d?e/d:0;
	}
	invcr[0]=dtotal?etotal/dtotal:0;
}
static double invcr2csizes(const double *invcr, const char *src_depths, int iw, int ih, int nch, double *csizes)//invcr & csizes have nch+1=5 elements
{
	double usize=0;
	int dtotal=0;
	for(int k=0;k<nch;++k)
	{
		int d=src_depths[k];
		csizes[k+1]=iw*ih*d*invcr[k+1]/8;
		dtotal+=d;
		usize+=(double)iw*ih*d/8;
	}
	csizes[0]=(double)iw*ih*dtotal*invcr[0]/8;
	return usize;
}
static void center_image(void)
{
	int wndw2, wndh2;

	if(!im1)
		return;
	wndw2=wndw, wndh2=wndh-17;
	if((double)wndw2/wndh2>=(double)im1->iw/im1->ih)//window AR > image AR: fit height
	{
		if(wndh2>0)
			imzoom=(double)wndh2/im1->ih;
	}
	else//window AR < image AR: fit width
		imzoom=(double)wndw2/im1->iw;
	wpx=(im1->iw-wndw2/imzoom)*0.5;//center image
	wpy=(im1->ih-wndh2/imzoom)*0.5;
	imagecentered=1;
}
static void zoom_at(int xs, int ys, double factor)
{
	const double tolerance=1e-2;
	wpx+=xs/imzoom*(1-1/factor);
	wpy+=ys/imzoom*(1-1/factor);
	imzoom*=factor;
	if(fabs(imzoom-1)<tolerance)
		imzoom=1;

	imagecentered=0;
}
#if 0
static void calc_csize_ans_separate(Image const *image, size_t *csizes)
{
	int maxdepth, maxlevels, *hist, res;
	unsigned state;

	if(!csizes)
		return;
	maxdepth=calc_maxdepth(image, 0), maxlevels=1<<maxdepth;
	hist=(int*)malloc((maxlevels+1)*sizeof(int));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	res=image->iw*image->ih;
	state=0x10000;
	for(int kc=0;kc<3;++kc)//naive way: one dedicated hist for each channel
	{
		int depth, nlevels, nusedlevels, sum;
		size_t csize;

		depth=image->depth[kc], nlevels=1<<depth;
		memset(hist, 0, nlevels*sizeof(int));
		for(int k=0;k<res;++k)//calc histogram
		{
			int sym=image->data[k<<2|kc]+(nlevels>>1);
			if((unsigned)sym>=(unsigned)nlevels)
				LOG_ERROR("Symbol OOB");
			++hist[sym];
		}
		nusedlevels=0;
		for(int ks=0;ks<nlevels;++ks)
			nusedlevels+=hist[ks]!=0;
		sum=0;
		for(int ks=0, ks2=0;ks<nlevels;++ks)//quantize & accumulate CDF
		{
			int freq=hist[ks];
			hist[ks]=(int)((long long)sum*(0x10000-nusedlevels)/res)+ks2;
			ks2+=freq!=0;
			sum+=freq;
		}
		hist[nlevels]=0x10000;

		csize=0;
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
		int lgv=FLOOR_LOG2((unsigned)val);
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

	{
		int ctx=0;
		for(int k=0;k<(int)_countof(ans_qlevels_u);++k)//TODO: binary search
			ctx+=energy>ans_qlevels_u[k];
		return ctx;
	}
}
void calc_csize_ans_energy(Image const *image, size_t *csizes)
{
	int maxdepth, histsize, *stats, ctxhist;
	unsigned state;

	if(!csizes)
		return;
	maxdepth=calc_maxdepth(image, 0);
	histsize=1<<maxdepth;//number of possibilities
	stats=(int*)malloc((histsize+1LL)*sizeof(int[_countof(calcsize_ans_qlevels)+1]));
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
	state=0x10000;
	for(int kc=0;kc<3;++kc)			//step 4: estimate compressed size
	{
		int ctxhist[34]={0};
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
	for(int kh=0;kh<(int)_countof(calcsize_ans_qlevels)+1;++kh)		//step 3: accumulate & quantize CDFs
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
	for(int k=0;k<(int)_countof(residues);++k)
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
	int hweight[34]={0};
	double csizes[14][3]={0}, optcsize[3]={0};
	static const int predcolors[]=
	{
		(int)0xFFFFFFFF,//	zero		undefined	white
		(int)0xFF0000FF,//	left		180		red
		(int)0xFFFF0000,//	top		90		blue
		(int)0xFF00D000,//	average0	135		green
		(int)0xFF009000,//	select		135		green
		(int)0xFF005000,//	grad		135		green
		(int)0xFF000000,//	weighted	undefined	black
		(int)0xFFFF00FF,//	topright	45		pink
		(int)0xFF00FF00,//	topleft		135		green
		(int)0xFF000080,//	leftleft	180		red
		(int)0xFF008080,//	average1	157.5		mustard
		(int)0xFF808000,//	average2	112.5		marine
		(int)0xFF800080,//	average3	67.5		violet
		(int)0xFF808080,//	average4	undefined	grey
	};
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
	for(int ctx=0;ctx<(int)_countof(hweight);++ctx)	//2: get histogram sums
	{
		int *hist=stats+maxlevels*ctx;
		for(int ks=0;ks<maxlevels;++ks)
			hweight[ctx]+=hist[ks];
	}
	
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
#endif
static void calc_csize_stateful(Image const *image, int *hist_full, double *entropy)
{
	if(ec_method==ECTX_GRCTX)
	{
		//long long bitsizes[4*16*2]={0};
		int nch=(image->depth[0]!=0)+(image->depth[1]!=0)+(image->depth[2]!=0)+(image->depth[3]!=0);
		int maxdepth=image->depth[0];
		if(maxdepth<image->depth[1])maxdepth=image->depth[1];
		if(maxdepth<image->depth[2])maxdepth=image->depth[2];
		if(maxdepth<image->depth[3])maxdepth=image->depth[3];
		const int nctx=MODELNCTX;
	//	int nctx=(maxdepth+MODELPREC)<<1;
	//	int nctx=1<<MODELCTXBITS;//X
		int nlevels=1<<maxdepth, half=nlevels>>1, mask=nlevels-1;
		int hsize=(int)sizeof(int)*nch*nctx<<maxdepth;
		int *hists=(int*)malloc(hsize);
		int bufsize=(int)sizeof(short[4*4])*(image->iw+16);//4 padded rows * 4 channels max
		short *pixels=(short*)malloc(bufsize);
		if(!hists||!pixels)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(hists, 0, hsize);
		memset(pixels, 0, bufsize);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+((image->iw+16LL)*((ky-0LL)&3)+8)*4,
				pixels+((image->iw+16LL)*((ky-1LL)&3)+8)*4,
				pixels+((image->iw+16LL)*((ky-2LL)&3)+8)*4,
				pixels+((image->iw+16LL)*((ky-3LL)&3)+8)*4,
			};
			int sW[4]={0};
			int xrun[4]={0};
			int nbypass0[4]={0}, rbypass[4]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=4)
			{
				short
					*NNNE	=rows[3]+1*4,
					*NNEE	=rows[2]+2*4,
					*NW	=rows[1]-1*4,
					*N	=rows[1]+0*4,
					*NE	=rows[1]+1*4,
					*NEE	=rows[1]+2*4,
					*NEEE	=rows[1]+3*4,
					*NEEEE	=rows[1]+4*4,
					*WW	=rows[0]-2*4,
					*W	=rows[0]-1*4,
					*curr	=rows[0]+0*4;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)WW;
				(void)W;
				(void)curr;
				for(int kc=0;kc<4;++kc)
				{
					if(image->depth[kc])
					{
#if 1
						//if(ky==1)//
						//if(!kc&&((ky==0&&kx==2)||(ky==0&&kx==3)||(ky==1&&kx==0)))//
						//	printf("");
					//	int ctx=sW[kc]<<8>>image->depth[kc]&255;
						int ctx=FLOOR_LOG2(sW[kc]*sW[kc]+1);
						//if(!kc)//
						//	printf("");
						if(ctx>MODELNCTX-1)
							ctx=MODELNCTX-1;
						int sym=image->data[idx|kc];
						sym<<=32-image->depth[kc];
						sym>>=32-image->depth[kc];
						++hists[(kc*nctx+ctx)<<maxdepth|((sym+(1<<image->depth[kc]>>1))&((1<<image->depth[kc])-1))];
						sym=sym<<1^sym>>31;
						//if(sym>>image->depth[kc])
						//	LOG_WARNING("sym %d > max %d", sym, (1<<image->depth[kc])-1);
						curr[kc]=sW[kc]=(2*sW[kc]+(sym<<MODELPREC)+MAXVAR(NEE[kc], NEEE[kc]))>>2;
#endif
#if 0
						int delta=image->data[idx<<2|kc];
						int sym=(delta<<1^delta>>31)&mask;
						int ctx=FLOOR_LOG2(sW[kc]*sW[kc]+1);
						++hists[(kc*nctx+ctx)<<maxdepth|sym];
						curr[kc]=sW[kc]=(2*sW[kc]+(sym<<MODELPREC)+MAXVAR(NEE[kc], NEEE[kc]))>>2;
#endif
#if 0
						int delta=image->data[idx<<2|kc];
						delta+=1<<image->depth[kc]>>1;
						delta&=(1<<image->depth[kc])-1;
						++hists[(kc*nctx+sW[kc])<<maxdepth|delta];
						sW[kc]=delta<<MODELCTXBITS>>image->depth[kc];
#endif
					}
				}
				rows[0]+=4;
				rows[1]+=4;
				rows[2]+=4;
				rows[3]+=4;
			}
		}
		free(pixels);
		double ctxsizes[4]={0};
		for(int kc=0;kc<nch;++kc)
		{
			double chsize=0;
			for(int kctx=0;kctx<nctx;++kctx)
			{
				int ctx=kc*nctx+kctx;
				int *curr_hist=hists+((ptrdiff_t)ctx<<maxdepth);
				int sum=0;
				for(int ks=0;ks<nlevels;++ks)
					sum+=curr_hist[ks];
				if(!sum)
					continue;
				//if(sum<nlevels*2)
				//{
				//	chsize+=bitsizes[ctx]/8.;
				//	continue;
				//}

				//exact estimate
#if 0
				double e=0, norm=1./sum;
				for(int ks=0;ks<nlevels;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
						e-=freq*log2(freq*norm);
				}
#endif

				//simulate bin tracking guard
#if 1
				int count=0;
				for(int ks=0;ks<nlevels;++ks)
					count+=curr_hist[ks]!=0;
				double e=0;
				for(int ks=0;ks<nlevels;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
					{
						int prob=(int)((long long)freq*(0x1000LL-count)/sum)+1;
						e-=freq*log2(prob*(1./0x1000));
						//if(!kc&&!kctx)console_log("%3d %7d\n", ks, freq);
					}
				}
				e/=8;
				//console_log("%3d %12.2lf\n", kc*nctx+kctx, e);
#endif
				//simulate unconditional guard		nlevels = 512 with SubGopt
#if 0
				double e=0, norm=1./sum*(0x1000-nlevels)/0x1000;
				for(int ks=0, ks2=0;ks<nlevels;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
						e-=freq*log2(freq*norm+1./0x1000);
				}
				e/=8;
#endif
				if(e<0)
					LOG_ERROR("C%d ctx%d e %lf bytes", kc, kctx, e);
				chsize+=e;
#if 1
				double hsize=0;
				int cdfW=0;
				int sum2=0;
				const int probbits=12;
				int codelen=probbits+1, CDFlevels=1<<probbits;
				int nlevels2=1<<image->depth[kc], half2=nlevels2>>1, mask2=nlevels2-1;
				for(int ks=0, ks2=0;ks<nlevels2;++ks)//calc model stat overhead
				{
					int sym=((ks>>1^-(ks&1))+half2)&mask2;
					int freq=curr_hist[sym];
					int cdf=sum2*((1ULL<<probbits)-count)/sum+ks2;
					ks2+=freq!=0;
					int csym=cdf-cdfW;
					if(ks&&CDFlevels)//CDF[0] is always zero
					{
						//GR
						int nbypass=FLOOR_LOG2(CDFlevels);
						if(ks>1)
							nbypass-=7;
						if(nbypass<0)
							nbypass=0;
						hsize+=(csym>>nbypass)+1+nbypass;

						//variable-base code
						//if(codelen>1)
						//{
						//	int codelen0=codelen-1;
						//	codelen-=(CDFlevels/((1<<codelen0)-1)+1)*codelen0 < (CDFlevels/((1<<codelen)-1)+1)*codelen;
						//}
						//hsize+=(csym/((1<<codelen)-1)+1)*codelen;

						//variable-base code v2
						//if(codelen>1)
						//{
						//	int codelen0=codelen>>1;
						//	codelen-=(CDFlevels/((1<<codelen0)-1)+1)*codelen0 < (CDFlevels/((1<<codelen)-1)+1)*codelen;
						//}
						//hsize+=(csym/((1<<codelen)-1)+1)*codelen;

						//naive pack
						//hsize+=probbits;

						//hsize+=csym+1;//unary code

						//hsize+=log2(csym+1);		//X  no stop bit in binary code
					}
					CDFlevels-=csym;
					cdfW=cdf;
					sum2+=freq;
				}
				chsize+=hsize/8;
				ctxsizes[kc]+=hsize/8;
#endif
			}
			double usize=(double)image->iw*image->ih*image->src_depth[0];
			entropy[kc]=usize?chsize*8*8/usize:0;//entropy = cbitsize*8/ubitsize
		}
		(void)ctxsizes;
		//if(loud_transforms)//
		//	LOG_WARNING("%10.3lf | %10.3lf %10.3lf %10.3lf %10.3lf",
		//		ctxsizes[0]+ctxsizes[1]+ctxsizes[2]+ctxsizes[3],
		//		ctxsizes[0],
		//		ctxsizes[1],
		//		ctxsizes[2],
		//		ctxsizes[3]);
		free(hists);
	}
	else if(ec_method==ECTX_GR)
	{
#if 1
		//long long bitsize_details[8]={0};//4 channels max * {syms, runs}
#endif
		long long bitsizes[4]={0};
		int bufsize=(int)sizeof(short[4*4*4])*(image->iw+16);//4 padded rows * 4 channels max * {pixels, errors, pastyrun, futureyrun}
		short *pixels=(short*)malloc(bufsize);
		ptrdiff_t skipbufsize=(ptrdiff_t)image->iw*image->ih<<2;
		//char *skipbuf=(char*)malloc(skipbufsize);
		if(!pixels)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(pixels, 0, bufsize);
		//memset(skipbuf, 0, skipbufsize);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			short *rows[]=
			{
				pixels+((image->iw+16LL)*((ky-0LL)&3)+8)*4*4,
				pixels+((image->iw+16LL)*((ky-1LL)&3)+8)*4*4,
				pixels+((image->iw+16LL)*((ky-2LL)&3)+8)*4*4,
				pixels+((image->iw+16LL)*((ky-3LL)&3)+8)*4*4,
			};
			int sW[4]={0};
			int xrun[4]={0};
			//int nxruns[4]={0};
			int nbypass0[4]={0}, rbypass[4]={0};
			for(int kx=0;kx<image->iw;++kx, ++idx)
			{
				short
					*NNNE	=rows[3]+1*4*4,
					*NNEE	=rows[2]+2*4*4,
					*NW	=rows[1]-1*4*4,
					*N	=rows[1]+0*4*4,
					*NE	=rows[1]+1*4*4,
					*NEE	=rows[1]+2*4*4,
					*NEEE	=rows[1]+3*4*4,
					*NEEEE	=rows[1]+4*4*4,
					*WW	=rows[0]-2*4*4,
					*W	=rows[0]-1*4*4,
					*curr	=rows[0]+0*4*4;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)WW;
				(void)W;
				(void)curr;
				for(int kc=0;kc<4;++kc)
				{
					if(image->depth[kc])
					{
						//simple xruns
#if 0
						int delta=image->data[idx<<2|kc];
						if(delta||kx>=image->iw-1)
						{
							if(xrun[kc])
							{
								bitsizes[kc]+=(long long)(0>>nbypass0[kc])+nbypass0[kc]+1;

								int nbypass2=FLOOR_LOG2(rbypass[kc]+1);
								bitsizes[kc]+=(long long)((xrun[kc]-1)>>nbypass2)+nbypass2+1;
								rbypass[kc]=(3*rbypass[kc]+xrun[kc])>>2;
								xrun[kc]=0;
							}
							int nbypass=FLOOR_LOG2(sW[kc]+1);
							int sym=delta<<1^delta>>31;
							bitsizes[kc]+=(long long)(sym>>nbypass)+nbypass+1;

							curr[kc+4]=sW[kc]=(2*sW[kc]+sym+NEEE[kc+4])>>2;//
						//	sW[kc]=(3*sW[kc]+sym)>>2;
						}
						else
						{
							if(!xrun[kc])
								nbypass0[kc]=FLOOR_LOG2(sW[kc]+1);
							++xrun[kc];
							curr[kc+4]=sW[kc]=(2*sW[kc]+0+NEEE[kc+4])>>2;//
						//	curr[kc+4]=0;
						}
#endif

						//horizontal and vertical zero-runs
#if 0
						int delta=image->data[idx<<2|kc];
						if(!delta)
						{
							if(!skipbuf[idx<<2|kc])
							{
								int xrun=0, yrun=0;
								for(int kx2=kx-1;kx2>=0&&!image->data[(idx-kx+kx2)<<2|kc];--kx2, ++xrun);//naive  O(n^2)
								for(int ky2=ky-1;ky2>=0&&!image->data[(idx+image->iw*(-ky+ky2))<<2|kc];--ky2, ++yrun);
								if(xrun<yrun)//yrun
								{
									int yscan=1;
									for(int ky2=ky+1;ky2<image->ih&&!image->data[(idx+image->iw*(-ky+ky2))<<2|kc];--ky2, ++yscan);
								}
								else//xrun
								{
								}
							}
						}
						else
						{
						}
#endif

						//horizontal zero-runs only
#if 0
						int delta=image->data[idx<<2|kc];
						if(delta)
						{
							while(xrun[kc])
							{
								int n=FLOOR_LOG2(xrun[kc]);
								int sym2=nxruns[kc]-n;

								int nbypass=FLOOR_LOG2(sW[kc]+1);
								int ncodebits=(sym2>>nbypass)+1LL+nbypass;
								bitsizes[kc]+=ncodebits;//GR-encode run deceleration = FLOOR_LOG2(remainder)
								curr[kc+4]=sW[kc]=(2*sW[kc]+sym2+NEEE[kc+4])>>2;

								//bitsize_details[kc+4]+=ncodebits;//

								xrun[kc]-=1<<n;
								nxruns[kc]=n;
							}
							nxruns[kc]=0;

							int sym=delta<<1^delta>>31;
							int nbypass=FLOOR_LOG2(sW[kc]+1);
							int ncodebits=(sym>>nbypass)+1LL+nbypass;
							bitsizes[kc]+=ncodebits;//GR-encode sym
							curr[kc+4]=sW[kc]=(2*sW[kc]+sym+NEEE[kc+4])>>2;

							//bitsize_details[kc+0]+=ncodebits;//
						}
						else
						{
							++xrun[kc];
							if(xrun[kc]>=(1<<nxruns[kc])||kx>=image->iw-1)
							{
								int nbypass=FLOOR_LOG2(sW[kc]+1);
								int ncodebits=(0>>nbypass)+1LL+nbypass;
								bitsizes[kc]+=ncodebits;//GR-encode run acceleration = 0
								curr[kc+4]=sW[kc]=(2*sW[kc]+0+NEEE[kc+4])>>2;

								//bitsize_details[kc+4]+=ncodebits;//

								nxruns[kc]+=xrun[kc]>1;//double the run
							}
							else
								curr[kc+4]=(2*W[kc+4]+0+NEEE[kc+4])>>2;
						}
#endif

						//X
#if 0
						int nbypass=abs(W[kc+4]);
						nbypass+=nbypass<4;
						nbypass=FLOOR_LOG2(nbypass);
						int run=MAXVAR(xrun, N[kc+8]);

						int val=image->data[idx<<2|kc];
						int sym=val<<1^val>>31;
						xrun+=!sym;

						curr[kc+0]=val;
						curr[kc+4]=(2*W[kc+4]+sym+NEEE[kc+4])>>2;
						curr[kc+8]=!sym?N[kc+8]+1:0;
#endif

						//pure GR estimation
#if 1

						int delta=image->data[idx<<2|kc];
						int sym=delta<<1^delta>>31;

						int nbypass=FLOOR_LOG2((sW[kc]>>6)+1);
						bitsizes[kc]+=(long long)(sym>>nbypass)+nbypass+1;

						curr[kc+4]=sW[kc]=(2*sW[kc]+(sym<<6)+MAXVAR(NEE[kc+4], NEEE[kc+4]))>>2;	//B

						//int vmax=MAXVAR(NEE[kc+4], NEEE[kc+4]);				//C
						//curr[kc+4]=sW[kc]=(2*sW[kc]+sym+MAXVAR(NE[kc+4], vmax))>>2;

					//	curr[kc+4]=sW[kc]=(2*sW[kc]+sym+MAXVAR(NEEE[kc+4], NEEEE[kc+4]))>>2;
					//	curr[kc+4]=sW[kc]=(32*sW[kc]+16*sym+16*NEEE[kc+4])/62;
					//	curr[kc+4]=sW[kc]=(2*sW[kc]+sym+NEEE[kc+4])>>2;				//A
					//	curr[kc+4]=sW[kc]=(2*MAXVAR(WW[kc+4], sW[kc])+sym+NEEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(MAXVAR(WW[kc+4], sW[kc])+sym+NE[kc+4]+MAXVAR(NEE[kc+4], NEEE[kc+4]))>>2;

						//int vmax=NEE[kc+4], vmin=NEEE[kc+4];
						//if(NEE[kc+4]<NEEE[kc+4])
						//	vmin=NEE[kc+4], vmax=NEEE[kc+4];
						//int val2=NE[kc+4];
						//CLAMP2(val2, vmin, vmax);
						//curr[kc+4]=sW[kc]=(2*sW[kc]+sym+val2)>>2;

					//	curr[kc+4]=sW[kc]=(sW[kc]+sym+NEE[kc+4]+NEEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(sW[kc]+sym+NE[kc+4]+NEEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(sW[kc]+sym+NEEE[kc+4])/3;
					//	curr[kc+4]=sW[kc]=(WW[kc+4]+sW[kc]+sym+NEEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(3*sW[kc]+2*sym+NE[kc+4]+2*NEEE[kc+4])>>3;
					//	curr[kc+4]=sW[kc]=(4*sW[kc]+2*sym+NNEE[kc+4]+NEEE[kc+4])>>3;
					//	curr[kc+4]=sW[kc]=(sW[kc]+NNNE[kc+4]+sym+NEEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(2*MINVAR(sW[kc], NE[kc+4])+sym+MAXVAR(NEE[kc+4], NEEE[kc+4]))>>2;
					//	curr[kc+4]=sW[kc]=abs(N[kc+4]-NE[kc+4])+abs(NE[kc+4]-NEE[kc+4])<abs(N[kc+4]-sym)+abs(NE[kc+4]-NEE[kc+4])?(sW[kc]+sym)>>1:(2*NE[kc+4]+sym+NEE[kc+4])>>2;
					//	curr[kc+4]=sW[kc]=(sW[kc]+2*MAXVAR(sW[kc], WW[kc+4])+sym)>>2;
					//	curr[kc+4]=sW[kc]=(3*sW[kc]+sym)>>2;
					//	curr[kc+4]=sW[kc]=(sym+NEEE[kc+4])>>1;

					//	curr[kc+0]=val;
#endif
					}
				}
				rows[0]+=4*4;
				rows[1]+=4*4;
				rows[2]+=4*4;
				rows[3]+=4*4;
			}
		}
		free(pixels);
		//free(skipbuf);
		{
			double usize;

			usize=(double)image->iw*image->ih*image->src_depth[0]; entropy[0]=usize?(bitsizes[0]<<3)/usize:0;//entropy = cbitsize*8/ubitsize
			usize=(double)image->iw*image->ih*image->src_depth[1]; entropy[1]=usize?(bitsizes[1]<<3)/usize:0;
			usize=(double)image->iw*image->ih*image->src_depth[2]; entropy[2]=usize?(bitsizes[2]<<3)/usize:0;
			usize=(double)image->iw*image->ih*image->src_depth[3]; entropy[3]=usize?(bitsizes[3]<<3)/usize:0;
		}
#if 0
		if(loud_transforms)
		{
			char msg[2048]={0};
			int nprinted=snprintf(msg, sizeof(msg)-1,
				"syms, runs\n"
				"T %8lld %8lld\n"
				"Y %8lld %8lld\n"
				"U %8lld %8lld\n"
				"V %8lld %8lld\n"
				"A %8lld %8lld\n",
				bitsize_details[0]+bitsize_details[1]+bitsize_details[2]+bitsize_details[3],
				bitsize_details[4]+bitsize_details[5]+bitsize_details[6]+bitsize_details[7],
				bitsize_details[0], bitsize_details[4],
				bitsize_details[1], bitsize_details[5],
				bitsize_details[2], bitsize_details[6],
				bitsize_details[3], bitsize_details[7]
			);
			copy_to_clipboard(msg, nprinted);
			messagebox(MBOX_OK, "Copied to clipboard", "%s", msg);
		}
#endif
	}
	else if(ec_method==ECTX_STATIC_O0)
	{
		int allocated=0;
		if(!hist_full)
		{
			int maxdepth=calc_maxdepth(image, 0);
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
	else if(ec_method==ECTX_STATIC_O1)
	{
		int nlevels[]=
		{
			1<<image->depth[0],
			1<<image->depth[1],
			1<<image->depth[2],
			1<<image->depth[3],
		};
		int half[]=
		{
			nlevels[0]>>1,
			nlevels[1]>>1,
			nlevels[2]>>1,
			nlevels[3]>>1,
		};
		int mask[]=
		{
			nlevels[0]-1,
			nlevels[1]-1,
			nlevels[2]-1,
			nlevels[3]-1,
		};
		int histsize=(int)sizeof(int)*(
			+nlevels[0]*nlevels[0]
			+nlevels[1]*nlevels[1]
			+nlevels[2]*nlevels[2]
			+nlevels[3]*nlevels[3]
		);
		int *hist=(int*)malloc(histsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(hist, 0, histsize);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int W[4]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=4)
			{
				int *currhist=hist;
				for(int kc=0;kc<image->nch;++kc)
				{
					int idx2=idx+kc;
				//	int N=ky?image->data[idx2-image->iw*4]:0;
					int target=image->data[idx2];
					int ctx=(W[kc]+half[kc])&mask[kc];
				//	int ctx=(N+W[kc]+nlevels[kc])>>1&mask[kc];
					
					++currhist[ctx<<image->depth[kc]|((target+half[kc])&mask[kc])];

					currhist+=nlevels[kc]*nlevels[kc];
					W[kc]=target;
				}
			}
		}
		double csizes[4]={0};
		int *currhist=hist;
		for(int kc=0;kc<4;++kc)
		{
			if(!image->depth[kc])
			{
				currhist+=nlevels[kc]*nlevels[kc];
				continue;
			}
			int totalsum=0;
			for(int ks0=0;ks0<nlevels[kc];++ks0)
			{
			//	int *hist3=currhist+((ptrdiff_t)ks0<<image->depth[kc]);
				int sum=0;
				for(int ks=0;ks<nlevels[kc];++ks)
					sum+=currhist[ks];
				if(!sum)
				{
					currhist+=nlevels[kc];
					continue;
				}
				double e=0, norm=1./sum;
				for(int ks=0;ks<nlevels[kc];++ks)
				{
					int freq=currhist[ks];
					if(freq)
						e-=freq*log2(freq*norm);
				}
				csizes[kc]+=e/8+84;//average GR overhead
				totalsum+=sum;
				currhist+=nlevels[kc];
			}
			if(totalsum!=image->iw*image->ih)
				LOG_ERROR("");
			double usize=(double)image->iw*image->ih*image->src_depth[0];
			entropy[kc]=usize?csizes[kc]*8*8/usize:0;//entropy = cbitsize*8/ubitsize

		//	currhist+=nlevels[kc]*nlevels[kc];
		}
	}
	else if(ec_method==ECTX_DWT)
	{
		int maxdepth=image->depth[0];
		UPDATE_MAX(maxdepth, image->depth[1]);
		UPDATE_MAX(maxdepth, image->depth[2]);
		UPDATE_MAX(maxdepth, image->depth[3]);
		int histsize=(int)sizeof(int)<<maxdepth;
		int *hist=(int*)malloc(histsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		double e8[4][4]={0};
		for(int kc=0;kc<image->nch;++kc)
		{
			int nlevels=1<<image->depth[kc], half=nlevels>>1, mask=nlevels-1;
			for(int kb=0;kb<4;++kb)
			{
				int x1=0, x2=0, y1=0, y2=0;
				switch(kb)
				{
				case 0:x1=image->iw>>1, x2=image->iw, y1=0, y2=image->ih>>1;break;
				case 1:x1=0, x2=image->iw>>1, y1=image->ih>>1, y2=image->ih;break;
				case 2:x1=image->iw>>1, x2=image->iw, y1=image->ih>>1, y2=image->ih;break;
				case 3:x1=0, x2=image->iw>>1, y1=0, y2=image->ih>>1;break;
				}
				memset(hist, 0, histsize);
				for(int ky=y1;ky<y2;++ky)
				{
					for(int kx=x1;kx<x2;++kx)
					{
						int val=image->data[(image->iw*ky+kx)<<2|kc];
						val+=half;
						val&=mask;
						++hist[val];
					}
				}
				double e=0, gain=1./((x2-x1)*(y2-y1));
				for(int ks=0;ks<nlevels;++ks)
				{
					int freq=hist[ks];
					if(freq)
						e-=freq*log2((double)freq*gain);
				}
				e8[kb][kc]=e*gain;
			}
		}
		free(hist);
		entropy[0]=(e8[0][0]+e8[1][0]+e8[2][0]+e8[3][0])/4;
		entropy[1]=(e8[0][1]+e8[1][1]+e8[2][1]+e8[3][1])/4;
		entropy[2]=(e8[0][2]+e8[1][2]+e8[2][2]+e8[3][2])/4;
		entropy[3]=(e8[0][3]+e8[1][3]+e8[2][3]+e8[3][3])/4;
	}
	else if(ec_method==ECTX_INTERLEAVED)
	{
		int maxdepth=image->depth[0];
		UPDATE_MAX(maxdepth, image->depth[1]);
		UPDATE_MAX(maxdepth, image->depth[2]);
		UPDATE_MAX(maxdepth, image->depth[3]);
		int histsize=(int)sizeof(int)<<maxdepth;
		int *hist=(int*)malloc(histsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		double e8[4][4]={0};
		for(int kc=0;kc<image->nch;++kc)
		{
			int nlevels=1<<image->depth[kc], half=nlevels>>1, mask=nlevels-1;
			for(int kb=0;kb<4;++kb)
			{
				int x0=0, y0=0;
				int count=0;
				switch(kb)
				{
				case 0:x0=1; y0=0;break;
				case 1:x0=0; y0=1;break;
				case 2:x0=1; y0=1;break;
				case 3:x0=0; y0=0;break;
				}
				memset(hist, 0, histsize);
				for(int ky=y0;ky<image->ih-1;ky+=2)
				{
					for(int kx=x0;kx<image->iw-1;kx+=2)
					{
						int val=image->data[(image->iw*ky+kx)<<2|kc];
						val+=half;
						val&=mask;
						++hist[val];
						++count;
					}
				}
				double e=0, gain=1./count;
				for(int ks=0;ks<nlevels;++ks)
				{
					int freq=hist[ks];
					if(freq)
						e-=freq*log2((double)freq*gain);
				}
				e8[kb][kc]=e*gain;
			}
		}
		free(hist);
		entropy[0]=(e8[0][0]+e8[1][0]+e8[2][0]+e8[3][0])/4;
		entropy[1]=(e8[0][1]+e8[1][1]+e8[2][1]+e8[3][1])/4;
		entropy[2]=(e8[0][2]+e8[1][2]+e8[2][2]+e8[3][2])/4;
		entropy[3]=(e8[0][3]+e8[1][3]+e8[2][3]+e8[3][3])/4;
	}
	else if(ec_method==ECTX_ABAC0||ec_method==ECTX_ABAC1)
		calc_csize_abac(image, ec_method==ECTX_ABAC1, entropy);
	else if(ec_method==ECTX_YUV422||ec_method==ECTX_YUV420)
	{
		int nlevels[]=
		{
			1<<image->depth[0],
			1<<image->depth[1],
			1<<image->depth[2],
		};
		int hoffsets[]=
		{
			0,
			nlevels[0],
			nlevels[0],
			nlevels[1],
			nlevels[1],
			nlevels[1],
			nlevels[1],
			nlevels[2],
			nlevels[2],
			nlevels[2],
			nlevels[2],
		};
		int hsum[_countof(hoffsets)-1]={0};
		for(int kc=1;kc<_countof(hoffsets);++kc)
			hoffsets[kc]+=hoffsets[kc-1];
		int hsize=hoffsets[_countof(hoffsets)-1]*(int)sizeof(int);
		int *hist=(int*)malloc(hsize);
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(hist, 0, hsize);
		int *hist_y=hist, *hist_u=hist+hoffsets[2], *hist_v=hist+hoffsets[6];
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx)
			{
				for(int kc=0;kc<3;++kc)
				{
					int val=image->data[(image->iw*ky+kx)<<2|kc]+(nlevels[kc]>>1);
					val&=nlevels[kc]-1;
					if(ec_method==ECTX_YUV420)
					{
						switch(kc)
						{
						case 0:
							if(!(kx&7)&&!(ky&7))//DC
								++hist_y[0<<image->depth[0]|val], ++hsum[0];
							else//AC
								++hist_y[1<<image->depth[0]|val], ++hsum[1];
							break;
						case 1:
							if(!(kx>>1&7)&&!(ky>>1&7))//DC
							{
								if(!(kx&1)&&!(ky&1))//base
									++hist_u[0<<image->depth[1]|val], ++hsum[2];
								else//interpolated
									++hist_u[1<<image->depth[1]|val], ++hsum[3];
							}
							else//AC
							{
								if(!(kx&1)&&!(ky&1))//base
									++hist_u[2<<image->depth[1]|val], ++hsum[4];
								else//interpolated
									++hist_u[3<<image->depth[1]|val], ++hsum[5];
							}
							break;
						case 2:
							if(!(kx>>1&7)&&!(ky>>1&7))//DC
							{
								if(!(kx&1)&&!(ky&1))//base
									++hist_v[0<<image->depth[2]|val], ++hsum[6];
								else//interpolated
									++hist_v[1<<image->depth[2]|val], ++hsum[7];
							}
							else//AC
							{
								if(!(kx&1)&&!(ky&1))//base
									++hist_v[2<<image->depth[2]|val], ++hsum[8];
								else//interpolated
									++hist_v[3<<image->depth[2]|val], ++hsum[9];
							}
							break;
						}
					}
					else//YUV422
					{
						switch(kc)
						{
						case 0:
							if(!(kx&7)&&!(ky&7))//DC
								++hist_y[0<<image->depth[0]|val], ++hsum[0];
							else//AC
								++hist_y[1<<image->depth[0]|val], ++hsum[1];
							break;
						case 1:
							if(!(kx>>1&7)&&!(ky&7))//DC
							{
								if(!(kx&1))//base
									++hist_u[0<<image->depth[1]|val], ++hsum[2];
								else//interpolated
									++hist_u[1<<image->depth[1]|val], ++hsum[3];
							}
							else//AC
							{
								if(!(kx&1))//base
									++hist_u[2<<image->depth[1]|val], ++hsum[4];
								else//interpolated
									++hist_u[3<<image->depth[1]|val], ++hsum[5];
							}
							break;
						case 2:
							if(!(kx>>1&7)&&!(ky&7))//DC
							{
								if(!(kx&1))//base
									++hist_v[0<<image->depth[2]|val], ++hsum[6];
								else//interpolated
									++hist_v[1<<image->depth[2]|val], ++hsum[7];
							}
							else//AC
							{
								if(!(kx&1))//base
									++hist_v[2<<image->depth[2]|val], ++hsum[8];
								else//interpolated
									++hist_v[3<<image->depth[2]|val], ++hsum[9];
							}
							break;
						}
					}
				}
			}
		}
		int eidxs[]=
		{
			0, 0,
			1, 1, 1, 1,
			2, 2, 2, 2,
		};
		//double e2[_countof(hoffsets)-1]={0};
		memset(entropy, 0, sizeof(double[4]));
		for(int kc=0;kc<_countof(hoffsets)-1;++kc)
		{
			int eidx=eidxs[kc], nlevels2=hoffsets[kc+1]-hoffsets[kc], sum=hsum[kc];
			for(int ks=0;ks<nlevels2;++ks)
			{
				int freq=hist[hoffsets[kc]+ks];
				if(freq)
					entropy[eidx]-=freq*log2((double)freq/sum);
			}
			//e2[kc]=calc_entropy(hist+hoffsets[kc], hoffsets[kc+1]-hoffsets[kc], hsum[kc]);
		}
		int res=image->iw*image->ih;
		entropy[0]/=res;
		entropy[1]/=res;
		entropy[2]/=res;
		//entropy[0]=(e2[0]+e2[1])*0.5;
		//entropy[1]=(e2[2]+e2[3]+e2[4]+e2[5])*0.25;
		//entropy[2]=(e2[6]+e2[7]+e2[8]+e2[9])*0.25;
		free(hist);
	}
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
	double entropy[4]={0};
	ThreadCtx *ctx=(ThreadCtx*)param;

	ctx->usize=image_getBMPsize(ctx->image);
	apply_selected_transforms(&ctx->image, 0, 1, 1);
	calc_csize_stateful(ctx->image, 0, entropy);
	for(int kc=0;kc<4;++kc)
	{
		int depth=ctx->image->src_depth[kc];
		double invCR=depth?entropy[kc]/depth:0;
		ctx->csize[kc]=invCR*ctx->image->iw*ctx->image->ih*ctx->image->src_depth[kc]/8;
	}
	free(ctx->image);
	return 0;
}
static void batch_test(void)
{
	const char *ext[]=
	{
		"PPM", "PGM", "PNM",
		"PNG",
		"JPG", "JPEG",
		"BMP",
		"TIF", "TIFF",
	};
	ArrayHandle path, filenames, q;
	int nthreads, maxlen;
	double t, total_usize=0, total_csize[4]={0};


	loud_transforms=0;
	path=dialog_open_folder();
	if(!path)
		return;
	filenames=get_filenames((char*)path->data, ext, _countof(ext), 1);
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
	nthreads=console_scan_int();
	t=time_sec();
	total_usize=0;
	maxlen=0;
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ARRAY_ALLOC(ThreadCtx, q, 0, 0, nthreads, 0);
	for(int k=0;k<(int)filenames->count;++k)
	{
		//multi-threaded
#if 1
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		Image *image=image_load((char*)fn2[0]->data, (int)fn2[0]->count);
		ThreadCtx ctx=
		{
			image,
			0, {0},
			k,
		};
		if(!image)
			continue;
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
				double csize=ptr->csize[0]+ptr->csize[1]+ptr->csize[2]+ptr->csize[3];
				fn2=(ArrayHandle*)array_at(&filenames, ptr->idx);
				console_log(
					"%5d/%5d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  BPD %8.4lf\n",
					(int)(k+1-(int)q->count+k2+1), (int)filenames->count, (char*)fn2[0]->data, (int)(maxlen-fn2[0]->count+1), "",
					ptr->usize, csize, ptr->csize[0], ptr->csize[1], ptr->csize[2],
					8.*csize/ptr->usize
				);
				//console_log(
				//	"%5d/%5d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  invCR %8.4lf%%\n",
				//	(int)(k+1-(int)q->count+k2+1), (int)filenames->count, (char*)fn2[0]->data, (int)(maxlen-fn2[0]->count+1), "",
				//	ptr->usize, csize, ptr->csize[0], ptr->csize[1], ptr->csize[2],
				//	100.*csize/ptr->usize
				//);
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
		double usize;
		int maxdepth, maxlevels;
		if(!image)
			continue;
		usize=image_getBMPsize(image), csize[3]={0};
		apply_selected_transforms(&image, 0, 1, 1);
		maxdepth=calc_maxdepth(image, 0);
		maxlevels=1<<maxdepth;
		int *hist=(int*)malloc(maxlevels*sizeof(int));
		for(int kc=0;kc<3;++kc)
		{
			double entropy, invCR;

			calc_histogram(image->data, image->iw, image->ih, kc, 0, image->iw, 0, image->ih, image->depth[kc], hist, 0);
			entropy=calc_entropy(hist, 1<<image->depth[kc], image->iw*image->ih);
			invCR=entropy/image->src_depth[kc];
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
	{
		char str[1024]={0};
		double ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
		double CR=total_usize/ctotal;
		t=time_sec()-t;
		int nprinted=snprintf(str, sizeof(str)-1,
			"%12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  BPD %8.4lf  %12.6lf sec  %12.6lf MB/s",
			total_usize, ctotal, total_csize[0], total_csize[1], total_csize[2], 8./CR, t, total_usize/(t*1024*1024)
		);
		console_log("Total UTYUV %s <- copied\n", str);
		copy_to_clipboard(str, nprinted);
		timedelta2str(g_buf, G_BUF_SIZE, t);
		console_log("Elapsed %s\n", g_buf);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d-%H:%M:%S");
		console_log("\nDone.  %s\n", g_buf);
		console_pause();
		console_end();
		loud_transforms=1;
	}
}

ptrdiff_t calc_psnr(const Image *im0, const Image *im1, double *ret_rmse, double *ret_psnr)
{
	if(im0->iw!=im1->iw||im0->ih!=im1->ih||im0->nch!=im1->nch)
	{
		if(ret_rmse)
		{
			ret_rmse[0]=INFINITE;
			ret_rmse[1]=INFINITE;
			ret_rmse[2]=INFINITE;
			ret_rmse[3]=INFINITE;
		}
		if(ret_psnr)
		{
			ret_psnr[0]=0;
			ret_psnr[1]=0;
			ret_psnr[2]=0;
			ret_psnr[3]=0;
		}
		return -1;
	}
	double rmse[4]={0}, psnr[4]={0};
	ptrdiff_t res=(ptrdiff_t)im0->iw*im0->ih, idx=-1;
	int half[]=
	{
		1<<im0->src_depth[0]>>1,
		1<<im0->src_depth[1]>>1,
		1<<im0->src_depth[2]>>1,
		1<<im0->src_depth[3]>>1,
	};
	for(ptrdiff_t k=0;k<res;++k)
	{
		int
			dr=im1->data[k<<2|0]-im0->data[k<<2|0],
			dg=im1->data[k<<2|1]-im0->data[k<<2|1],
			db=im1->data[k<<2|2]-im0->data[k<<2|2],
			da=im1->data[k<<2|3]-im0->data[k<<2|3];
		if(idx<0&&(dr|dg|db|da))
			idx=k;
		rmse[0]+=dr*dr;
		rmse[1]+=dg*dg;
		rmse[2]+=db*db;
		rmse[3]+=da*da;
	}
	for(int kc=0;kc<4;++kc)
	{
		rmse[kc]=sqrt(rmse[kc]/res);
		psnr[kc]=20*log10(half[kc]/rmse[kc]);
	}
	if(ret_rmse)
		memcpy(ret_rmse, rmse, sizeof(rmse));
	if(ret_psnr)
		memcpy(ret_psnr, psnr, sizeof(psnr));
	return idx;
}
typedef struct LossyCtxStruct
{
	Image *image;
	double usize, csize[4];
	ptrdiff_t idx;
	double rmse[4], psnr[4], sqe[4], ssim[4];
} LossyCtx;
static unsigned __stdcall sample_thread_lossy(void *param)
{
	double entropy[4]={0};
	LossyCtx *ctx=(LossyCtx*)param;
	Image *im0=0;

	ctx->usize=image_getBMPsize(ctx->image);
	image_copy(&im0, ctx->image);//save original
	if(!im0)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}

	apply_selected_transforms(&ctx->image, 0, 1, 0);//fwd

	calc_csize_stateful(ctx->image, 0, entropy);
	for(int kc=0;kc<4;++kc)
	{
		int depth=ctx->image->src_depth[kc];
		double invCR=depth?entropy[kc]/depth:0;
		ctx->csize[kc]=invCR*ctx->image->iw*ctx->image->ih*ctx->image->src_depth[kc]/8;
	}

	apply_selected_transforms(&ctx->image, 0, 0, 1);//inv

	calc_psnr(im0, ctx->image, ctx->rmse, ctx->psnr);
	measure_ssim_avg(im0, ctx->image, ctx->ssim);
	{
		ptrdiff_t res=(ptrdiff_t)ctx->image->iw*ctx->image->ih;
		ctx->sqe[0]=res*ctx->rmse[0]*ctx->rmse[0];
		ctx->sqe[1]=res*ctx->rmse[1]*ctx->rmse[1];
		ctx->sqe[2]=res*ctx->rmse[2]*ctx->rmse[2];
		ctx->sqe[3]=res*ctx->rmse[3]*ctx->rmse[3];
	}
	free(im0);
	free(ctx->image);
	return 0;
}
static void batch_test_lossy(void)
{
	const char *ext[]=
	{
		"PPM", "PGM", "PNM",
		"PNG",
		"JPG", "JPEG",
		"BMP",
		"TIF", "TIFF",
	};
	ArrayHandle path, filenames, q;
	int nthreads, maxlen;
	double t, total_usize=0, total_csize[4]={0}, total_sqe[4]={0};
	double total_ssim[4]={0};


	loud_transforms=0;
	path=dialog_open_folder();
	if(!path)
		return;
	filenames=get_filenames((char*)path->data, ext, _countof(ext), 1);
	if(!filenames)
	{
		array_free(&path);
		return;
	}

	DisableProcessWindowsGhosting();
	console_start();
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H:%M:%S");
	console_log("Lossy Batch Test  %s  %s\n", g_buf, (char*)path->data);
	array_free(&path);
	console_log("Enter number of threads: ");
	nthreads=console_scan_int();
	t=time_sec();
	total_usize=0;
	maxlen=0;
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ARRAY_ALLOC(LossyCtx, q, 0, 0, nthreads, 0);
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		Image *image=image_load((char*)fn2[0]->data, (int)fn2[0]->count);
		LossyCtx ctx=
		{
			image,
			0, {0},
			k,
		};
		if(!image)
			continue;
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
				LossyCtx *ptr=(LossyCtx*)array_at(&q, k2);
				handles[k2]=(void*)_beginthreadex(0, 0, sample_thread_lossy, ptr, 0, 0);
				if(!handles[k2])
				{
					LOG_ERROR("Thread alloc error");
					return;
				}
			}
			WaitForMultipleObjects((int)q->count, handles, TRUE, INFINITE);
			for(int k2=0;k2<(int)q->count;++k2)
			{
				LossyCtx *ptr=(LossyCtx*)array_at(&q, k2);
				double csize=ptr->csize[0]+ptr->csize[1]+ptr->csize[2]+ptr->csize[3];
				fn2=(ArrayHandle*)array_at(&filenames, ptr->idx);
				console_log(
					"%5d/%5d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  BPD %8.4lf  PSNR %10.4lf %10.4lf %10.4lf %10.4lf  SSIM%7.4lf\n",
					(int)(k+1-(int)q->count+k2+1), (int)filenames->count, (char*)fn2[0]->data, (int)(maxlen-fn2[0]->count+1), "",
					ptr->usize, csize, ptr->csize[0], ptr->csize[1], ptr->csize[2],
					8.*csize/ptr->usize,
					ptr->psnr[0], ptr->psnr[1], ptr->psnr[2], ptr->psnr[3],
					ptr->ssim[3]
				);
				//console_log(
				//	"%5d/%5d %s%*sUTYUV %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  invCR %8.4lf%%\n",
				//	(int)(k+1-(int)q->count+k2+1), (int)filenames->count, (char*)fn2[0]->data, (int)(maxlen-fn2[0]->count+1), "",
				//	ptr->usize, csize, ptr->csize[0], ptr->csize[1], ptr->csize[2],
				//	100.*csize/ptr->usize
				//);
				total_usize+=ptr->usize;
				total_csize[0]+=ptr->csize[0];
				total_csize[1]+=ptr->csize[1];
				total_csize[2]+=ptr->csize[2];
				total_csize[3]+=ptr->csize[3];
				total_sqe[0]+=ptr->sqe[0];
				total_sqe[1]+=ptr->sqe[1];
				total_sqe[2]+=ptr->sqe[2];
				total_sqe[3]+=ptr->sqe[3];
				total_ssim[0]+=ptr->ssim[0]*ptr->usize/(1024*1024);
				total_ssim[1]+=ptr->ssim[1]*ptr->usize/(1024*1024);
				total_ssim[2]+=ptr->ssim[2]*ptr->usize/(1024*1024);
				total_ssim[3]+=ptr->ssim[3]*ptr->usize/(1024*1024);
			}
			array_clear(&q);
		}
	}
	{
		char str[1024]={0};
		double ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
		double CR=total_usize/ctotal;
		double rmse=sqrt((total_sqe[0]+total_sqe[1]+total_sqe[2]+total_sqe[3])/total_usize);
		double psnr=20*log10(255/rmse);
		double ssim_norm=1024*1024/(double)total_usize;
		total_ssim[0]*=ssim_norm;
		total_ssim[1]*=ssim_norm;
		total_ssim[2]*=ssim_norm;
		total_ssim[3]*=ssim_norm;
		t=time_sec()-t;
		timedelta2str(g_buf, G_BUF_SIZE, t);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y%m%d_%H%M%S");
		int nprinted=snprintf(str, sizeof(str)-1,
			"%12.2lf %12.2lf %12.2lf %12.2lf %12.2lf  BPD %8.4lf  RMSE %10.4lf PSNR %10.4lf  %10.4lf sec  SSIM%7.4lf  dist%3d  %s",
			total_usize, ctotal, total_csize[0], total_csize[1], total_csize[2], 8./CR, rmse, psnr, t, total_ssim[3], g_dist, g_buf
		);
		console_log("Total UTYUV %s <- copied\n", str);
		copy_to_clipboard(str, nprinted);
		console_log("Elapsed %s\n", g_buf);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d-%H:%M:%S");
		console_log("\nDone.  %s\n", g_buf);
		console_pause();
		console_end();
		loud_transforms=1;
	}
}

static int customtransforms_getflag(unsigned char tid)
{
	return
	//	tid==CT_FWD_ADAPTIVE||
	//	tid==CT_INV_ADAPTIVE||
		tid==CT_FWD_CUSTOM||
		tid==CT_INV_CUSTOM||
		tid==ST_FWD_CUSTOM||
		tid==ST_INV_CUSTOM||
		tid==ST_CONVTEST||
		tid==ST_CONVTEST2||
		tid==ST_FWD_CC||
		tid==ST_INV_CC;
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
	if(tid<T_COUNT&&tid!=CST_INV_SEPARATOR)
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
					if((tid<CST_COMPARE)==(transforms->data[k]<CST_COMPARE))
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
//	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
//	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
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
	case CT_FWD_SubG_OPT:		a="C  Fwd OptSub";		break;
	case CT_INV_SubG_OPT:		a="C  Inv OptSub";		break;
	case CT_FWD_YCoCg_R:		a="C  Fwd YCoCg-R";		break;
	case CT_INV_YCoCg_R:		a="C  Inv YCoCg-R";		break;
	case CT_FWD_SUBGREEN:		a="C  Fwd SubGreen";		break;
	case CT_INV_SUBGREEN:		a="C  Inv SubGreen";		break;
	case CT_FWD_JPEG2000:		a="C  Fwd JPEG2000 RCT";	break;
	case CT_INV_JPEG2000:		a="C  Inv JPEG2000 RCT";	break;
	case CT_FWD_JPEG2000_MA:	a="C  Fwd JPEG2000 MA";		break;
	case CT_INV_JPEG2000_MA:	a="C  Inv JPEG2000 MA";		break;
	case CT_FWD_NBLI:		a="C  Fwd NBLI RCT";		break;
	case CT_INV_NBLI:		a="C  Inv NBLI RCT";		break;
//	case CT_FWD_J2K2:		a="C  Fwd J2K/2";		break;
//	case CT_INV_J2K2:		a="C  Inv J2K/2";		break;
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

	case CST_COMPARE:		a="Compare";			break;
	case CST_INV_SEPARATOR:		a="";				break;

//	case ST_PREPROC_GRAD:		a=" S Preproc Grad";		break;
//	case ST_PREPROC_X:		a=" S Preproc X";		break;
//	case ST_PREPROC_X2:		a=" S Preproc X2";		break;
		
	case ST_CONVTEST:		a=" S ConvTest";		break;
	case ST_CONVTEST2:		a=" S ConvTest2";		break;
	case ST_DIFF:			a="   Difference";		break;
	case ST_SSIM:			a="   SSIM";			break;
	case ST_FWD_PACKSIGN:		a=" S Fwd PackSign";		break;
	case ST_INV_PACKSIGN:		a=" S Inv PackSign";		break;
	case ST_FWD_BWTX:		a=" S Fwd BWT-X";		break;
	case ST_INV_BWTX:		a=" S Inv BWT-X";		break;
	case ST_FWD_BWTY:		a=" S Fwd BWT-Y";		break;
	case ST_INV_BWTY:		a=" S Inv BWT-Y";		break;
	case ST_FWD_MTF:		a=" S Fwd MTF";			break;
	case ST_INV_MTF:		a=" S Inv MTF";			break;
	case ST_FWD_PALETTE:		a=" S Fwd Palette";		break;
	case ST_INV_PALETTE:		a=" S Inv Palette";		break;
	case ST_FILT_MEDIAN33:		a=" S Filt Median33";		break;
	case ST_FILT_AV33:		a=" S Filt Av33";		break;
	case ST_FILT_DEINT422:		a=" S Filt DeInt422";		break;
	case ST_FILT_DEINT420:		a=" S Filt DeInt420";		break;
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
	case ST_FWD_OLS6:		a=" S Fwd OLS-6";		break;
	case ST_INV_OLS6:		a=" S Inv OLS-6";		break;
	case ST_FWD_GRFILT:		a=" S Fwd GR Filt";		break;
	case ST_INV_GRFILT:		a=" S Inv GR Filt";		break;
	case ST_FWD_L1CRCT:		a=" S Fwd L1 cRCT";		break;
	case ST_INV_L1CRCT:		a=" S Inv L1 cRCT";		break;
	case ST_FWD_OLS7:		a=" S Fwd L1";			break;
	case ST_INV_OLS7:		a=" S Inv L1";			break;
	case ST_FWD_L1BCRCT:		a=" S Fwd L1B cRCT";		break;
	case ST_INV_L1BCRCT:		a=" S Inv L1B cRCT";		break;
	case ST_FWD_OLS8:		a=" S Fwd L1B";			break;
	case ST_INV_OLS8:		a=" S Inv L1B";			break;
	case ST_FWD_PU:			a="CS Fwd PU";			break;
	case ST_INV_PU:			a="CS Inv PU";			break;
	case ST_FWD_CG3D:		a="CS Fwd CG3D";		break;
	case ST_INV_CG3D:		a="CS Inv CG3D";		break;
	case ST_FWD_CGCRCT:		a=" S Fwd CG cRCT";		break;
	case ST_INV_CGCRCT:		a=" S Inv CG cRCT";		break;
	case ST_FWD_OLS9:		a=" S Fwd OLS9";		break;
	case ST_INV_OLS9:		a=" S Inv OLS9";		break;
	case ST_FWD_SUB:		a=" S Fwd Sub W";		break;
	case ST_INV_SUB:		a=" S Inv Sub W";		break;
	case ST_FWD_CLAMPGRAD:		a=" S Fwd ClampGrad";		break;
	case ST_INV_CLAMPGRAD:		a=" S Inv ClampGrad";		break;
	case ST_FWD_CLEARTYPE:		a=" S Fwd ClearType";		break;
	case ST_INV_CLEARTYPE:		a=" S Inv ClearType";		break;
	case ST_FWD_QUANT:		a=" S Fwd Quantize";		break;
	case ST_INV_QUANT:		a=" S Inv Quantize";		break;
	case ST_FWD_AV4:		a=" S Fwd AV4";			break;
	case ST_INV_AV4:		a=" S Inv AV4";			break;
	case ST_FWD_SEL4:		a=" S Fwd Sel4";		break;
	case ST_INV_SEL4:		a=" S Inv Sel4";		break;
	case ST_FWD_SELECT:		a=" S Fwd Select";		break;
	case ST_INV_SELECT:		a=" S Inv Select";		break;
	case ST_FWD_CGPLUS:		a=" S Fwd CGplus";		break;
	case ST_INV_CGPLUS:		a=" S Inv CGplus";		break;
	case ST_FWD_CG422:		a=" S Fwd CG422";		break;
	case ST_INV_CG422:		a=" S Inv CG422";		break;
	case ST_FWD_CG420:		a=" S Fwd CG420";		break;
	case ST_INV_CG420:		a=" S Inv CG420";		break;
	case ST_FWD_AV2:		a=" S Fwd (N+W)/2";		break;
	case ST_INV_AV2:		a=" S Inv (N+W)/2";		break;
	case ST_FWD_MIX2:		a=" S Fwd MIX2";		break;
	case ST_INV_MIX2:		a=" S Inv MIX2";		break;
	case ST_FWD_MIXN:		a=" S Fwd MIX N";		break;
	case ST_INV_MIXN:		a=" S Inv MIX N";		break;
//	case ST_FWD_AV3:		a=" S Fwd AV3";			break;
//	case ST_INV_AV3:		a=" S Inv AV3";			break;
	case ST_FWD_WGRAD:		a="CS Fwd WGrad";		break;
	case ST_INV_WGRAD:		a="CS Inv WGrad";		break;
//	case ST_FWD_WGRAD2:		a="CS Fwd WGrad2";		break;
//	case ST_INV_WGRAD2:		a="CS Inv WGrad2";		break;
	case ST_FWD_WGRAD3:		a="CS Fwd WGrad3";		break;
	case ST_INV_WGRAD3:		a="CS Inv WGrad3";		break;
	case ST_FWD_WGRAD4CCRCT:	a=" S Fwd WG4C cRCT";		break;
	case ST_INV_WGRAD4CCRCT:	a=" S Inv WG4C cRCT";		break;
	case ST_FWD_WGRAD4C:		a=" S Fwd WG4C";		break;
	case ST_INV_WGRAD4C:		a=" S Inv WG4C";		break;
	case ST_FWD_WGRAD4:		a=" S Fwd WGrad4";		break;
	case ST_INV_WGRAD4:		a=" S Inv WGrad4";		break;
	case ST_FWD_WGRAD5:		a=" S Fwd WGrad5";		break;
	case ST_INV_WGRAD5:		a=" S Inv WGrad5";		break;
	case ST_FWD_WGRAD6:		a=" S Fwd WGrad6";		break;
	case ST_INV_WGRAD6:		a=" S Inv WGrad6";		break;
	case ST_FWD_WGRAD7:		a=" S Fwd WP7";			break;
	case ST_INV_WGRAD7:		a=" S Inv WP7";			break;
	case ST_FWD_SSE:		a=" S Fwd SSE";			break;
	case ST_INV_SSE:		a=" S Inv SSE";			break;
//	case ST_FWD_WMIX:		a=" S Fwd WMIX";		break;
//	case ST_INV_WMIX:		a=" S Inv WMIX";		break;
	case ST_FWD_TABLE:		a=" S Fwd Table";		break;
	case ST_INV_TABLE:		a=" S Inv Table";		break;
	case ST_FWD_LWAV:		a=" S Fwd LWAV";		break;
	case ST_INV_LWAV:		a=" S Inv LWAV";		break;
//	case ST_FWD_ECOEFF:		a=" S Fwd E-Coeff";		break;
//	case ST_INV_ECOEFF:		a=" S Inv E-Coeff";		break;
//	case ST_FWD_AVERAGE:		a=" S Fwd Average";		break;
//	case ST_INV_AVERAGE:		a=" S Inv Average";		break;
//	case ST_FWD_MULTISTAGE:		a=" S Fwd Multistage";		break;
//	case ST_INV_MULTISTAGE:		a=" S Inv Multistage";		break;
	case ST_FWD_ZIPPER:		a=" S Fwd Zipper";		break;
	case ST_INV_ZIPPER:		a=" S Inv Zipper";		break;
//	case ST_FWD_DIR:		a=" S Fwd Dir";			break;
//	case ST_INV_DIR:		a=" S Inv Dir";			break;
	case ST_FWD_WC:			a=" S Fwd WC";			break;
	case ST_INV_WC:			a=" S Inv WC";			break;
	case ST_FWD_CUSTOM4:		a=" S Fwd CUSTOM4";		break;
	case ST_INV_CUSTOM4:		a=" S Inv CUSTOM4";		break;
	case ST_FWD_CUSTOM3:		a="CS Fwd CUSTOM3";		break;
	case ST_INV_CUSTOM3:		a="CS Inv CUSTOM3";		break;
	case ST_FWD_NBLIC:		a=" S Fwd NBLIC";		break;
	case ST_INV_NBLIC:		a=" S Inv NBLIC";		break;
	case ST_FWD_CALIC:		a=" S Fwd CALIC";		break;
	case ST_INV_CALIC:		a=" S Inv CALIC";		break;
	case ST_FWD_LEGALLCG:		a=" S Fwd LeGallCG";		break;
	case ST_INV_LEGALLCG:		a=" S Inv LeGallCG";		break;
	case ST_FWD_WP:			a=" S Fwd JXL WP";		break;
	case ST_INV_WP:			a=" S Inv JXL WP";		break;
//	case ST_FWD_WP2:		a=" S Fwd WP2";			break;
//	case ST_INV_WP2:		a=" S Inv WP2";			break;
	case ST_FWD_WPU:		a="CS Fwd WPU";			break;
	case ST_INV_WPU:		a="CS Inv WPU";			break;
//	case ST_FWD_DEFERRED:		a=" S Fwd DEFERRED";		break;
//	case ST_INV_DEFERRED:		a=" S Inv DEFERRED";		break;
	case ST_FWD_MM:			a=" S Fwd MM";			break;
	case ST_INV_MM:			a=" S Inv MM";			break;
	case ST_FWD_CUSTOM:		a=" S Fwd CUSTOM";		break;
	case ST_INV_CUSTOM:		a=" S Inv CUSTOM";		break;
	case ST_FWD_CC:			a=" S Fwd CC";			break;
	case ST_INV_CC:			a=" S Inv CC";			break;
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
#endif
	case ST_FWD_LAZY:		a=" S Fwd Lazy DWT";		break;
	case ST_INV_LAZY:		a=" S Inv Lazy DWT";		break;
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
	case ST_FWD_WHT4:		a=" S Fwd WHT4";		break;
	case ST_INV_WHT4:		a=" S Inv WHT4";		break;
	case ST_FWD_WHT8:		a=" S Fwd WHT8";		break;
	case ST_INV_WHT8:		a=" S Inv WHT8";		break;
	case ST_FWD_HAAR8:		a=" S Fwd Haar8";		break;
	case ST_INV_HAAR8:		a=" S Inv Haar8";		break;
	case ST_FWD_FDCT8:		a=" S Fwd FDCT8";		break;
	case ST_INV_FDCT8:		a=" S Inv FDCT8";		break;
	default:			a="ERROR";			break;
	}
	{
		long long c0=0;
		if(highlight)
			c0=set_text_colors(highlight);
		if(place<0)
			GUIPrint(0, x, y, g_uiscale, "%s", a);
		else
			GUIPrint(0, x, y, g_uiscale, "%d: %s", place, a);
		if(highlight)
			set_text_colors(c0);
	}
}

static int send_image_separate_subpixels(Image const *image, unsigned *txid_r, unsigned *txid_g, unsigned *txid_b, unsigned *txid_a)
{
	ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
	int shift[]=
	{
		MAXVAR(0, image->depth[0]-8),
		MAXVAR(0, image->depth[1]-8),
		MAXVAR(0, image->depth[2]-8),
		MAXVAR(0, image->depth[3]-8),
	};
	unsigned char *temp_r=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_g=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_b=(unsigned char*)malloc(res*sizeof(char[4]));
	unsigned char *temp_a=(unsigned char*)malloc(res*sizeof(char[4]));
	if(!temp_r||!temp_g||!temp_b||!temp_a)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
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
	int nv, nf;
	float *vertices;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};

	if(image->iw*image->ih>1024*1024)
		return;
	nv=image->iw*image->ih*3*6, nf=nv*5;//subpixel count * 6 vertices
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	vertices=(float*)cpuv[0]->data;
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
	int nv, nf;
	float *vertices;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};

	if(image->iw*image->ih>1024*1024)
		return;
	nv=(image->iw-1)*(image->ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 6 vertices * 5 floats
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	vertices=(float*)cpuv[0]->data;
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
	int nv, nf;
	float *vertices;
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};

	if(image->iw*image->ih>1024*1024)
		return;
	nv=(image->iw-1)*(image->ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 2 triangles * 3 vertices * 5 floats
	if(!*cpuv||(int)cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	vertices=(float*)cpuv[0]->data;
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
	{
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
}
static void chart_dwthist_update(Image const *image, int kc, int kband, int x1, int x2, int y1, int y2)
{
	memset(hist+((size_t)kband<<8), 0, 256LL*sizeof(int));
	x1=MAXVAR(x1, 0);
	x2=MINVAR(x2, image->iw);
	y1=MAXVAR(y1, 0);
	y2=MINVAR(y2, image->ih);
	{
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
			{
				double entropy=calc_entropy(hist_full, image->depth[kc], count);
				blockCR[kband]=(float)(image->src_depth[kc]/entropy);
			}
		}
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
	int ix=(int)(((long long)xs*iw+(wndw>>1))/wndw);
	int iy=(int)(((long long)ys*ih+(wndh>>1))/wndh);
	move_box_in_window(ix, dx, 0, iw, bounds+0, bounds+1);
	move_box_in_window(iy, dy, 0, ih, bounds+2, bounds+3);
}
static void jh_calchist(int *jhist, int nbits, Image const *image, int x1, int x2, int y1, int y2)
{
	int idx;
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
				idx=(b<<nbits|g)<<nbits|r;

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
				idx=(v2<<nbits|v1)<<nbits|v0;

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
				idx=(v2<<nbits|v1)<<nbits|v0;

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
				idx=(v2<<nbits|v1)<<nbits|v0;

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
	{
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
					{
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
						for(int kt=0;kt<(int)_countof(triangles);++kt)
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
								float *edge, t;
								int mid;

								edge=(float*)ARRAY_APPEND(*edges, 0, 1, 1, 0);
								t=(level-v[0][0][3])/(v[2][0][3]-v[0][0][3]);
								*edge++=MIX(v[0][0][0], v[2][0][0], t)*gains[0];
								*edge++=MIX(v[0][0][1], v[2][0][1], t)*gains[1];
								*edge++=MIX(v[0][0][2], v[2][0][2], t)*gains[2];
								*edge++=0;
								*edge++=0;
						
								mid=level>v[1][0][3];
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
	}
}
static void chart_jointhist_update(Image const *image, unsigned *txid)
{
	int nlevels=1<<jointhist_nbits, hsize=nlevels*nlevels*nlevels;
	int *jhist;
	int bounds[4]={0};

	//jointhistogram(image, iw, ih, jointhist_nbits, &jointhist, space_not_color);
	if(jointhist)
	{
		if((int)jointhist->count<hsize)
			ARRAY_APPEND(jointhist, 0, hsize-jointhist->count, 1, 0);
	}
	else
		ARRAY_ALLOC(int, jointhist, 0, hsize, 0, 0);
	jhist=(int*)jointhist->data;
#if 1
	memset(jhist, 0, jointhist->esize*jointhist->count);
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
	{
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
static void analysis_update(Image const *image)
{
	int bounds[4]={0};//x1, x2, y1, y2

	jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, image->iw, image->ih, bounds);
	c18_analyze(image, bounds[0], bounds[1], bounds[2], bounds[3], &analysis_info);
	
	++pred_hist[0][analysis_info.predidx[0]];
	++pred_hist[1][analysis_info.predidx[1]];
	++pred_hist[2][analysis_info.predidx[2]];
	if(pred_hist[0][analysis_info.predidx[0]]>=512)
	{
		for(int k=0;k<PRED_COUNT;++k)
			pred_hist[0][k]>>=1;
	}
	if(pred_hist[1][analysis_info.predidx[1]]>=512)
	{
		for(int k=0;k<PRED_COUNT;++k)
			pred_hist[1][k]>>=1;
	}
	if(pred_hist[2][analysis_info.predidx[2]]>=512)
	{
		for(int k=0;k<PRED_COUNT;++k)
			pred_hist[2][k]>>=1;
	}
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

void apply_transform(Image **pimage, int tid, int hasRCT)
{
	Image *image=*pimage;
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
	case CT_FWD_SubG_OPT:		colortransform_subg_opt(image, 1);			break;
	case CT_INV_SubG_OPT:		colortransform_subg_opt(image, 0);			break;
	case CT_FWD_SUBGREEN:		colortransform_subtractgreen(image, 1);			break;
	case CT_INV_SUBGREEN:		colortransform_subtractgreen(image, 0);			break;
	case CT_FWD_JPEG2000:		colortransform_JPEG2000(image, 1);			break;
	case CT_INV_JPEG2000:		colortransform_JPEG2000(image, 0);			break;
	case CT_FWD_JPEG2000_MA:	colortransform_JPEG2000_MA(image, 1);			break;
	case CT_INV_JPEG2000_MA:	colortransform_JPEG2000_MA(image, 0);			break;
	case CT_FWD_NBLI:		colortransform_NBLI(image, 1);				break;
	case CT_INV_NBLI:		colortransform_NBLI(image, 0);				break;
//	case CT_FWD_J2K2:		colortransform_J2K2(image, 1);				break;
//	case CT_INV_J2K2:		colortransform_J2K2(image, 0);				break;
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
	case ST_FWD_WC:			pred_WC(image);						break;
	case ST_INV_WC:			pred_WC(image);						break;
	case ST_FWD_CUSTOM4:		pred_lossyconv(image);					break;
	case ST_INV_CUSTOM4:		pred_lossyconv(image);					break;
	case ST_FWD_CUSTOM3:		custom3_apply(image, 1, pred_ma_enabled, &c3_params);	break;
	case ST_INV_CUSTOM3:		custom3_apply(image, 0, pred_ma_enabled, &c3_params);	break;
	case ST_FWD_CUSTOM:		pred_custom(image, 1, pred_ma_enabled, custom_params);	break;
	case ST_INV_CUSTOM:		pred_custom(image, 0, pred_ma_enabled, custom_params);	break;
	case ST_CONVTEST:		pred_convtest(image);					break;
	case ST_CONVTEST2:		pred_convtest(image);					break;
	case ST_DIFF:
		if(im1->iw==im0->iw&&im1->ih==im0->ih)
		{
			for(ptrdiff_t k=0, n=(ptrdiff_t)im0->iw*im0->ih*4;k<n;k+=4)
			{
				im1->data[k+0]-=im0->data[k+0];
				im1->data[k+1]-=im0->data[k+1];
				im1->data[k+2]-=im0->data[k+2];
			}
		}
		break;
	case ST_SSIM:
		measure_ssim_map(im0, im1);
		break;
	case ST_FWD_CC:			conv_custom(image);					break;
	case ST_INV_CC:			conv_custom(image);					break;
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
	case ST_FWD_NBLIC:		pred_nblic(image, 1);					break;
	case ST_INV_NBLIC:		pred_nblic(image, 0);					break;
	case ST_FWD_CALIC:		pred_calic(image, 1, pred_ma_enabled);			break;
	case ST_INV_CALIC:		pred_calic(image, 0, pred_ma_enabled);			break;
	case ST_FWD_LEGALLCG:		pred_LeGallCG(image, 1);				break;
	case ST_INV_LEGALLCG:		pred_LeGallCG(image, 0);				break;
	case ST_FWD_WP:			pred_jxl_apply(image, 1, pred_ma_enabled, jxlparams_i16);break;
	case ST_INV_WP:			pred_jxl_apply(image, 0, pred_ma_enabled, jxlparams_i16);break;
//	case ST_FWD_WP2:		pred_divfreeWP(image, 1);				break;
//	case ST_INV_WP2:		pred_divfreeWP(image, 0);				break;
	case ST_FWD_WPU:		pred_WPU(image, 1);					break;
	case ST_INV_WPU:		pred_WPU(image, 0);					break;
//	case ST_FWD_DEFERRED:		pred_wp_deferred(image, 1);				break;
//	case ST_INV_DEFERRED:		pred_wp_deferred(image, 0);				break;
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
	case ST_FWD_OLS6:		pred_ols6(image, 1);					break;
	case ST_INV_OLS6:		pred_ols6(image, 0);					break;
	case ST_FWD_GRFILT:		pred_grfilt(image, 1);					break;
	case ST_INV_GRFILT:		pred_grfilt(image, 0);					break;
	case ST_FWD_L1CRCT:		pred_l1crct(image, 1);					break;
	case ST_INV_L1CRCT:		pred_l1crct(image, 0);					break;
	case ST_FWD_OLS7:		pred_ols7(image, 1);					break;
	case ST_INV_OLS7:		pred_ols7(image, 0);					break;
	case ST_FWD_L1BCRCT:		pred_ols8_crct(image, 1);				break;
	case ST_INV_L1BCRCT:		pred_ols8_crct(image, 0);				break;
	case ST_FWD_OLS8:		pred_ols8(image, 1);					break;
	case ST_INV_OLS8:		pred_ols8(image, 0);					break;
	case ST_FWD_PACKSIGN:		packsign(image, 1);					break;
	case ST_INV_PACKSIGN:		packsign(image, 0);					break;
	case ST_FWD_BWTX:		prep_BWT_x(pimage, 1);					break;
	case ST_INV_BWTX:		prep_BWT_x(pimage, 0);					break;
	case ST_FWD_BWTY:		prep_BWT_y(pimage, 1);					break;
	case ST_INV_BWTY:		prep_BWT_y(pimage, 0);					break;
	case ST_FWD_MTF:		pred_MTF(image, 1);					break;
	case ST_INV_MTF:		pred_MTF(image, 0);					break;
	case ST_FWD_PALETTE:		pred_palette(image, 1);					break;
	case ST_INV_PALETTE:		pred_palette(image, 0);					break;
	case ST_FILT_MEDIAN33:		filt_median33(image);					break;
	case ST_FILT_AV33:		filt_av33(image);					break;
	case ST_FILT_DEINT422:		filt_deint422(image);					break;
	case ST_FILT_DEINT420:		filt_deint420(image);					break;
	case ST_FWD_PU:			pred_PU(image, 1);					break;
	case ST_INV_PU:			pred_PU(image, 0);					break;
	case ST_FWD_CG3D:		pred_CG3D(image, 1, pred_ma_enabled);			break;
	case ST_INV_CG3D:		pred_CG3D(image, 0, pred_ma_enabled);			break;
	case ST_FWD_CGCRCT:		pred_cg_crct(image, 1, pred_ma_enabled);		break;
	case ST_INV_CGCRCT:		pred_cg_crct(image, 0, pred_ma_enabled);		break;
	case ST_FWD_OLS9:		pred_ols9(image, 1);					break;
	case ST_INV_OLS9:		pred_ols9(image, 0);					break;
	case ST_FWD_SUB:		pred_sub(image, 1);					break;
	case ST_INV_SUB:		pred_sub(image, 0);					break;
	case ST_FWD_CLAMPGRAD:		pred_clampgrad(image, 1, pred_ma_enabled);		break;
	case ST_INV_CLAMPGRAD:		pred_clampgrad(image, 0, pred_ma_enabled);		break;
	case ST_FWD_CLEARTYPE:		pred_cleartype(image, 1);				break;
	case ST_INV_CLEARTYPE:		pred_cleartype(image, 0);				break;
	case ST_FWD_QUANT:		pred_artifact(image, 1);				break;
	case ST_INV_QUANT:		pred_artifact(image, 0);				break;
	case ST_FWD_SEL4:		pred_sel4(image, 1);					break;
	case ST_INV_SEL4:		pred_sel4(image, 0);					break;
	case ST_FWD_SELECT:		pred_select(image, 1);					break;
	case ST_INV_SELECT:		pred_select(image, 0);					break;
	case ST_FWD_CGPLUS:		pred_cgplus(image, 1);					break;
	case ST_INV_CGPLUS:		pred_cgplus(image, 0);					break;
	case ST_FWD_CG422:		pred_CG422(image, 1);					break;
	case ST_INV_CG422:		pred_CG422(image, 0);					break;
	case ST_FWD_CG420:		pred_CG420(image, 1);					break;
	case ST_INV_CG420:		pred_CG420(image, 0);					break;
	case ST_FWD_AV2:		pred_av2(image, 1);					break;
	case ST_INV_AV2:		pred_av2(image, 0);					break;
	case ST_FWD_MIX2:		pred_mix2(image, 1);					break;
	case ST_INV_MIX2:		pred_mix2(image, 0);					break;
	case ST_FWD_MIXN:		pred_mixN(image, 1);					break;
	case ST_INV_MIXN:		pred_mixN(image, 0);					break;
//	case ST_FWD_AV3:		pred_av3(image, 1);					break;
//	case ST_INV_AV3:		pred_av3(image, 0);					break;
	case ST_FWD_WGRAD:		pred_wgrad(image, 1, hasRCT);				break;
	case ST_INV_WGRAD:		pred_wgrad(image, 0, hasRCT);				break;
//	case ST_FWD_WGRAD2:		pred_wgrad2(image, 1);					break;
//	case ST_INV_WGRAD2:		pred_wgrad2(image, 0);					break;
	case ST_FWD_WGRAD3:		pred_wgrad3(image, 1, hasRCT);				break;
	case ST_INV_WGRAD3:		pred_wgrad3(image, 0, hasRCT);				break;
	case ST_FWD_WGRAD4CCRCT:	pred_wgrad4c_crct(image, 1);				break;
	case ST_INV_WGRAD4CCRCT:	pred_wgrad4c_crct(image, 0);				break;
	case ST_FWD_WGRAD4C:		pred_wgrad4c(image, 1);					break;
	case ST_INV_WGRAD4C:		pred_wgrad4c(image, 0);					break;
	case ST_FWD_WGRAD4:		pred_wgrad4(image, 1);					break;
	case ST_INV_WGRAD4:		pred_wgrad4(image, 0);					break;
	case ST_FWD_WGRAD5:		pred_wgrad5(image, 1);					break;
	case ST_INV_WGRAD5:		pred_wgrad5(image, 0);					break;
	case ST_FWD_WGRAD6:		pred_wgrad6(image, 1);					break;
	case ST_INV_WGRAD6:		pred_wgrad6(image, 0);					break;
	case ST_FWD_WGRAD7:		pred_wpred7(image, 1);					break;
	case ST_INV_WGRAD7:		pred_wpred7(image, 0);					break;
	case ST_FWD_SSE:		pred_sse(image, 1);					break;
	case ST_INV_SSE:		pred_sse(image, 0);					break;
//	case ST_FWD_WMIX:		pred_wmix(image, 1);					break;
//	case ST_INV_WMIX:		pred_wmix(image, 0);					break;
	case ST_FWD_TABLE:		pred_table(image, 1);					break;
	case ST_INV_TABLE:		pred_table(image, 0);					break;
	case ST_FWD_LWAV:		pred_lwav(image, 1);					break;
	case ST_INV_LWAV:		pred_lwav(image, 0);					break;
//	case ST_FWD_ECOEFF:		pred_ecoeff(image, 1, pred_ma_enabled);			break;
//	case ST_INV_ECOEFF:		pred_ecoeff(image, 0, pred_ma_enabled);			break;
	case ST_FWD_AV4:		pred_av4(image, 1, pred_ma_enabled);			break;
	case ST_INV_AV4:		pred_av4(image, 0, pred_ma_enabled);			break;
//	case ST_FWD_MULTISTAGE:		pred_multistage(image, 1, pred_ma_enabled);		break;
//	case ST_INV_MULTISTAGE:		pred_multistage(image, 0, pred_ma_enabled);		break;
	case ST_FWD_ZIPPER:		pred_zipper(image, 1, pred_ma_enabled);			break;
	case ST_INV_ZIPPER:		pred_zipper(image, 0, pred_ma_enabled);			break;
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
	case ST_FWD_WHT4:		image_wht4_fwd(image);					break;
	case ST_INV_WHT4:		image_wht4_inv(image);					break;
	case ST_FWD_WHT8:		image_wht8_fwd(image);					break;
	case ST_INV_WHT8:		image_wht8_inv(image);					break;
	case ST_FWD_HAAR8:		image_haar8_fwd(image);					break;
	case ST_INV_HAAR8:		image_haar8_inv(image);					break;
	case ST_FWD_FDCT8:		image_fdct8_fwd(image);					break;
	case ST_INV_FDCT8:		image_fdct8_inv(image);					break;
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

	case ST_FWD_LAZY:
	case ST_INV_LAZY:
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
				if(!image->depth[kc])
					continue;
				switch(tid)
				{
				case ST_FWD_LAZY:      dwt2d_lazy_fwd		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_INV_LAZY:      dwt2d_lazy_inv		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_FWD_HAAR:      dwt2d_haar_fwd		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_INV_HAAR:      dwt2d_haar_inv		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_FWD_SQUEEZE:   dwt2d_squeeze_fwd	(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_INV_SQUEEZE:   dwt2d_squeeze_inv	(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_FWD_LEGALL53:  dwt2d_legall53_fwd	(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_INV_LEGALL53:  dwt2d_legall53_inv	(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_FWD_CDF97:     dwt2d_cdf97_fwd		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
				case ST_INV_CDF97:     dwt2d_cdf97_inv		(image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
			//	case ST_FWD_GRAD_DWT:  dwt2d_grad_fwd   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
			//	case ST_INV_GRAD_DWT:  dwt2d_grad_inv   (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
			//	case ST_FWD_EXPDWT:    dwt2d_exp_fwd    (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;//TODO use customdwtparams instead of sharing allcustomparam_st
			//	case ST_INV_EXPDWT:    dwt2d_exp_inv    (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
			//	case ST_FWD_CUSTOM_DWT:dwt2d_custom_fwd (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
			//	case ST_INV_CUSTOM_DWT:dwt2d_custom_inv (image->data+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
				}
				if(tid!=ST_FWD_LAZY&&tid!=ST_INV_LAZY)
				{
					int inv=tid&1;
					image->depth[kc]+=(char)(!inv-inv);
					UPDATE_MAX(image->depth[kc], image->src_depth[kc]);
					//if(image->depth[kc]>=18)
					//	printf("");
					UPDATE_MIN(image->depth[kc], 18);
				}
			}
			array_free(&sizes);
			free(temp);
#if 0
			if(loud_transforms
			//	&&transforms->count==1
			)
			{
				int maxdepth=image->depth[0];
				UPDATE_MAX(maxdepth, image->depth[1]);
				UPDATE_MAX(maxdepth, image->depth[2]);
				UPDATE_MAX(maxdepth, image->depth[3]);
				int histsize=(int)sizeof(int)<<maxdepth;
				int *hist=(int*)malloc(histsize);
				if(!hist)
				{
					LOG_ERROR("Alloc error");
					return;
				}
				double csizes[3][4]={0}, entropy[3][4]={0};
				for(int kc=0;kc<image->nch;++kc)
				{
					int nlevels=1<<image->depth[kc], half=nlevels>>1, mask=nlevels-1;
					for(int kb=0;kb<3;++kb)
					{
						int x1=0, x2=0, y1=0, y2=0;
						switch(kb)
						{
						case 0:x1=image->iw>>1, x2=image->iw, y1=0, y2=image->ih>>1;break;
						case 1:x1=0, x2=image->iw>>1, y1=image->ih>>1, y2=image->ih;break;
						case 2:x1=image->iw>>1, x2=image->iw, y1=image->ih>>1, y2=image->ih;break;
						}
						memset(hist, 0, histsize);
						for(int ky=y1;ky<y2;++ky)
						{
							for(int kx=x1;kx<x2;++kx)
							{
								int val=image->data[(image->iw*ky+kx)<<2|kc];
								val+=half;
								val&=mask;
								++hist[val];
							}
						}
						double e=0, gain=1./((x2-x1)*(y2-y1));
						for(int ks=0;ks<nlevels;++ks)
						{
							int freq=hist[ks];
							if(freq)
								e-=freq*log2((double)freq*gain);
						}
						csizes[kb][kc]=e/8;
						entropy[kb][kc]=csizes[kb][kc]*gain;
					}
				}
				free(hist);

				const int len=4096;
				char *str=(char*)malloc(len);
				int printed=0;
				const char cnames[]="YUVA";
				printed+=snprintf(str+printed, (size_t)len-printed-1, "T, NE, SW, SE:\n");
				double etotal=0, ctotal=0;
				for(int kc=0;kc<image->nch;++kc)
				{
					for(int kb=0;kb<3;++kb)
					{
						etotal+=entropy[kb][kc];
						ctotal+=csizes[kb][kc];
					}
				}
				etotal/=image->nch*3;
				printed+=snprintf(str+printed, (size_t)len-printed-1, "T %8.4lf%% %8.4lf%% %8.4lf%% %8.4lf%% %12.2lf %12.2lf %12.2lf %12.2lf\n",
					etotal*100.,
					(entropy[0][0]+entropy[0][1]+entropy[0][2]+entropy[0][3])*100./image->nch,
					(entropy[1][0]+entropy[1][1]+entropy[1][2]+entropy[1][3])*100./image->nch,
					(entropy[2][0]+entropy[2][1]+entropy[2][2]+entropy[2][3])*100./image->nch,
					ctotal,
					csizes[0][0]+csizes[0][1]+csizes[0][2]+csizes[0][3],
					csizes[1][0]+csizes[1][1]+csizes[1][2]+csizes[1][3],
					csizes[2][0]+csizes[2][1]+csizes[2][2]+csizes[2][3]
				);
				for(int kc=0;kc<image->nch;++kc)
					printed+=snprintf(str+printed, (size_t)len-printed-1, "%c %8.4lf%% %8.4lf%% %8.4lf%% %8.4lf%% %12.2lf %12.2lf %12.2lf %12.2lf\n",
						cnames[kc],
						(entropy[0][kc]+entropy[1][kc]+entropy[2][kc])*100./3,
						entropy[0][kc]*100., entropy[1][kc]*100., entropy[2][kc]*100.,
						csizes[0][kc]+csizes[1][kc]+csizes[2][kc],
						csizes[0][kc], csizes[1][kc], csizes[2][kc]
					);
				copy_to_clipboard(str, printed);
				messagebox(MBOX_OK, "Copied to clipboard", "%s", str);
				free(str);
			}
#endif
		}
		break;
//	case ST_FWD_DEC_DWT:   dwt2d_dec_fwd((char*)image, iw, ih);	break;
//	case ST_INV_DEC_DWT:   dwt2d_dec_inv((char*)image, iw, ih);	break;
	}//switch
}
void apply_selected_transforms(Image **pimage, int rct_only, int applyfwd, int applyinv)
{
	int hasRCT=0;

	if(!transforms)
		return;
	for(int k=0;k<(int)transforms->count;++k)
	{
		unsigned char tid=transforms->data[k];
		if(tid==CST_COMPARE)
		{
			unsigned char tid2=CST_COMPARE;
			for(k=0;k<(int)transforms->count;++k)//look for two spatial transforms
			{
				tid=transforms->data[k];
				if(tid>CST_COMPARE)
					break;
			}
			for(++k;k<(int)transforms->count;++k)
			{
				tid2=transforms->data[k];
				if(tid2>CST_COMPARE)
					break;
			}
			if(tid>CST_COMPARE&&tid2>CST_COMPARE&&tid!=tid2)
			{
				Image *image2=0;
				for(k=0;k<(int)transforms->count;++k)//apply RCTs
				{
					unsigned char tid3=transforms->data[k];
					hasRCT|=tid<CST_INV_SEPARATOR;
					if(tid3==tid||tid3==tid2)
						break;
					apply_transform(pimage, tid3, hasRCT);
				}
				image_copy(&image2, *pimage);
				apply_transform(pimage, tid, 0);
				apply_transform(&image2, tid2, 0);
				int res=pimage[0]->iw*pimage[0]->ih<<2;
				//int amplitudes[]=
				//{
				//	1<<pimage[0]->depth[0]>>3,
				//	1<<pimage[0]->depth[1]>>3,
				//	1<<pimage[0]->depth[2]>>3,
				//	1<<pimage[0]->depth[3]>>3,
				//};
				for(int kp=0;kp<res;++kp)
					pimage[0]->data[kp]=abs(pimage[0]->data[kp])-abs(image2->data[kp]);
				//{
				//	int val=image->data[kp];
				//	int val2=image2->data[kp];
				//	if(val==val2)
				//		image->data[kp]=0;
				//	else
				//	{
				//		val=abs(val);
				//		val2=abs(val2);
				//		if(val<val2)
				//			image->data[kp]=-amplitudes[kp&3];
				//		else
				//			image->data[kp]=amplitudes[kp&3];
				//	}
				//}
				free(image2);
				for(++k;k<(int)transforms->count;++k)
				{
					unsigned char tid3=transforms->data[k];
					hasRCT|=tid3<CST_INV_SEPARATOR;
					if(rct_only&&tid3>CST_INV_SEPARATOR)
						continue;
					if(tid3!=tid&&tid3!=tid2)
						apply_transform(pimage, tid3, hasRCT);
				}
			}
			return;
		}
	}
	for(int k=0;k<(int)transforms->count;++k)
	{
		unsigned char tid=transforms->data[k];
		hasRCT|=tid<CST_INV_SEPARATOR;
		if(rct_only&&tid>CST_INV_SEPARATOR)
			continue;
		if(tid&1?!applyinv:!applyfwd)
			continue;
		apply_transform(pimage, tid, hasRCT);
		//calc_depthfromdata(image->data, image->iw, image->ih, image->depth, image->src_depth);//X  depth must depend only on src_depth and applied RCTs, so that preds can apply MA
	}//for
}
void update_image(void)//apply selected operations on original image, calculate CRs, and export
{
	if(!im0)
		return;
	image_copy(&im1, im0);
	apply_selected_transforms(&im1, 0, 1, 1);

	//do not modify im1 beyond this point
	
	{
		void *ptr=realloc(im_export, sizeof(char[4])*im1->iw*im1->ih);
		if(!ptr)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		im_export=(unsigned char*)ptr;
	}
	{
		int shift[]=
		{
			MAXVAR(0, im1->depth[0]-8),
			MAXVAR(0, im1->depth[1]-8),
			MAXVAR(0, im1->depth[2]-8),
			MAXVAR(0, im1->depth[3]-8),
		};
		int c0=0, c1=1, c2=2;
		switch(viewmode)
		{
		case VIEW_C0:
			c1=0; c2=0;
			break;
		case VIEW_C1:
			c0=1; c2=1;
			break;
		case VIEW_C2:
			c0=2; c1=2;
			break;
		}
		if(view_ma)
			memset(shift, 0, sizeof(shift));
		if(ec_method==ECTX_ABAC0)
		{
			int maxdepth=im1->depth[c0];
			if(maxdepth<im1->depth[c1])maxdepth=im1->depth[1];
			if(maxdepth<im1->depth[c2])maxdepth=im1->depth[2];
		//	if(maxdepth<im1->depth[3])maxdepth=im1->depth[3];
			if(abacvis_range<2)
			{
				abacvis_low=0;
				abacvis_range=1<<maxdepth;
			}
			int halfs[]=
			{
				1<<im1->depth[c0]>>1,
				1<<im1->depth[c1]>>1,
				1<<im1->depth[c2]>>1,
			//	1<<im1->depth[3]>>1,
			};
			int abacvis_hi=abacvis_low+abacvis_range, abacvis_mid=abacvis_low+(abacvis_range>>1);
			for(ptrdiff_t k=0, res=(ptrdiff_t)im1->iw*im1->ih*4;k<res;k+=4)
			{
				int r=im1->data[k|c0]+halfs[0];
				int g=im1->data[k|c1]+halfs[1];
				int b=im1->data[k|c2]+halfs[2];
				if(r<abacvis_low||r>=abacvis_hi)
					r=128;
				else if(r>=abacvis_mid)
					r=255;
				else
					r=0;

				if(g<abacvis_low||g>=abacvis_hi)
					g=128;
				else if(g>=abacvis_mid)
					g=255;
				else
					g=0;

				if(b<abacvis_low||b>=abacvis_hi)
					b=128;
				else if(b>=abacvis_mid)
					b=255;
				else
					b=0;
				im_export[k|0]=r;
				im_export[k|1]=g;
				im_export[k|2]=b;
				im_export[k|3]=0xFF;
			}
		}
		else
		{
			for(ptrdiff_t k=0, res=(ptrdiff_t)im1->iw*im1->ih*4;k<res;k+=4)
			{
				im_export[k|0]=(unsigned char)((im1->data[k|c0]>>shift[c0])+128);
				im_export[k|1]=(unsigned char)((im1->data[k|c1]>>shift[c1])+128);
				im_export[k|2]=(unsigned char)((im1->data[k|c2]>>shift[c2])+128);
				im_export[k|3]=0xFF;
			}
		}
	}
	//image_export_uint8(im1, &im_export, 1, 0);

	{
		int maxdepth=calc_maxdepth(im1, 0), maxlevels=1<<maxdepth;
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
		//if(combCRhist_max==1||combCRhist_max<combCRhist[combCRhist_idx][k])
		if(fabsf(combCRhist_max-1)<1e-7||combCRhist_max<combCRhist[combCRhist_idx][k])
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
					int sym, freq;

					sym=im1->data[k<<2|kc]+(nlevels>>1);
					sym=CLAMP(0, sym, nlevels-1);
					freq=hist_full[sym];
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
			{
				float vmax=0;
				for(ptrdiff_t k=0;k<res;++k)
				{
					if(fabsf(vmax)<1e-9||vmax<fbuf[k<<2|0])
						vmax=fbuf[k<<2|0];
					if(vmax<fbuf[k<<2|1])
						vmax=fbuf[k<<2|1];
					if(vmax<fbuf[k<<2|2])
						vmax=fbuf[k<<2|2];
				}
				if(fabsf(vmax)>1e-9)
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
			}
			free(fbuf);
		}
		break;
	case VIS_HISTOGRAM:
		memset(histmax, 0, sizeof(histmax));
		memset(hist, 0, sizeof(hist));
		chart_hist_update(im1, 0, im1->iw, 0, im1->ih, hist, histmax, 0);
		break;
	case VIS_MODEL:
		{
			modeldepth=im1->depth[0];
			if(modeldepth<im1->depth[1])modeldepth=im1->depth[1];
			if(modeldepth<im1->depth[2])modeldepth=im1->depth[2];
			if(modeldepth<im1->depth[3])modeldepth=im1->depth[3];
			modelnch=im1->nch;
			modelnctx=MODELNCTX;
		//	modelnctx=(modeldepth+MODELPREC)<<1;
		//	modelnctx=1<<MODELCTXBITS;
			modelhistsize=(int)sizeof(int)*modelnch*modelnctx<<modeldepth;
			void *p=(int*)realloc(modelhist, modelhistsize);
			if(!p)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			modelhist=p;
			memset(modelhist, 0, modelhistsize);
			void *p2=(int*)realloc(modelmhist, modelhistsize);
			if(!p)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			modelmhist=p2;
			memset(modelmhist, 0, modelhistsize);
			int nlevels=1<<modeldepth, half=nlevels>>1, mask=nlevels-1;
			int psize=(int)sizeof(int[4*4])*(im1->iw+16);//4 padded rows * 4 channels max
			int *pixels=(int*)malloc(psize);
			if(!psize)
			{
				LOG_ERROR("Alloc error");
				return;
			}
			memset(pixels, 0, psize);
			for(int ky=0, idx=0;ky<im1->ih;++ky)
			{
				int *rows[]=
				{
					pixels+((im1->iw+16LL)*((ky-0LL)&3)+8)*4,
					pixels+((im1->iw+16LL)*((ky-1LL)&3)+8)*4,
					pixels+((im1->iw+16LL)*((ky-2LL)&3)+8)*4,
					pixels+((im1->iw+16LL)*((ky-3LL)&3)+8)*4,
				};
				int eW[4]={0};
				for(int kx=0;kx<im1->iw;++kx, idx+=4)
				{
#if 1
					for(int kc=0;kc<modelnch;++kc)
					{
						int ctx=FLOOR_LOG2(eW[kc]*eW[kc]+1);
						if(ctx>MODELNCTX-1)
							ctx=MODELNCTX-1;
						int sym=im1->data[idx|kc];
						sym<<=32-im1->depth[kc];
						sym>>=32-im1->depth[kc];
						++modelhist[(modelnctx*kc+ctx)<<modeldepth|((sym+half)&mask)];
						sym=sym<<1^sym>>31;
						rows[0][kc]=eW[kc]=(2*eW[kc]+(sym<<MODELPREC)+MAXVAR(rows[1][kc+2*4], rows[1][kc+3*4]))>>2;
					}
#endif
#if 0
					for(int kc=0;kc<modelnch;++kc)
					{
						int ctx=FLOOR_LOG2(eW[kc]*eW[kc]+1);
						//if(ctx>=1&&ctx<8)
						//	printf("");
						int sym=im1->data[idx|kc];
						sym+=half;
						sym&=mask;
						//if((ptrdiff_t)((modelnctx*kc+ctx)<<modeldepth|sym)>=(ptrdiff_t)(modelhistsize/sizeof(int)))//
						//	LOG_ERROR("");
						++modelhist[(modelnctx*kc+ctx)<<modeldepth|sym];
						sym-=half;
						sym=sym<<1^sym>>31;
						rows[0][kc]=eW[kc]=(2*eW[kc]+(sym<<MODELPREC)+MAXVAR(rows[1][kc+2*4], rows[1][kc+3*4]))>>2;
					}
#endif
#if 0
					for(int kc=0;kc<modelnch;++kc)
					{
						int sym=im1->data[idx|kc];
						sym+=1<<im1->depth[kc]>>1;
						sym&=(1<<im1->depth[kc])-1;
						++modelhist[(modelnctx*kc+eW[kc])<<modeldepth|sym];
						eW[kc]=sym<<MODELCTXBITS>>im1->depth[kc]&((1<<MODELCTXBITS)-1);
					}
#endif
					rows[0]+=4;
					rows[1]+=4;
					rows[2]+=4;
					rows[3]+=4;
				}
			}
			free(pixels);
			modelstatoverhead=0;
			for(int kc=0;kc<modelnch*modelnctx;++kc)
			{
				int *curr_hist=modelhist+((ptrdiff_t)kc<<modeldepth);
				int sum=0;
				for(int ks=0;ks<nlevels;++ks)
					sum+=curr_hist[ks];
				if(!sum)
				{
					modelcsizes[kc]=0;
					modelmsizes[kc]=0;
					continue;
				}
#if 1
				int count=0;
				for(int ks=0;ks<nlevels;++ks)
					count+=curr_hist[ks]!=0;
				double e=0, norm=1./sum;
				for(int ks=0;ks<nlevels;++ks)//simulate 12-bit precision
				{
					int freq=curr_hist[ks];
					if(freq)
					{
						int prob=(int)((long long)freq*(0x1000LL-count)/sum)+1;
						e-=freq*log2(prob*(1./0x1000));
					}
				}
				e/=8;
#endif
#if 0
				double e=0, norm=1./sum;
				for(int ks=0;ks<nlevels;++ks)//raw entropy
				{
					int freq=curr_hist[ks];
					if(freq)
						e-=freq*log2(freq*norm);
				}
				e/=8;
#endif
				//if(e>sum)
				//	LOG_ERROR("Context %d %lf > sum %d", kc, e, sum);
				modelcsizes[kc]=e;


				double mean=0;
				for(int ks=0;ks<nlevels;++ks)//calc mean
				{
					int freq=curr_hist[ks];
					mean+=(double)ks*freq;
				}
				mean*=norm;
				modelmeans[kc]=mean;
				double variance=0;
				for(int ks=0;ks<nlevels;++ks)//calc variance
				{
					int freq=curr_hist[ks];
					double val=ks-mean;
					variance+=(double)freq*val*val;//var = E[(X-mean)^2],    E[X] = (p0*X0+p1*X1+...)/pden
				}
				variance*=norm;
				double sdev=sqrt(variance);
				modelsdevs[kc]=sdev;
				
				int *curr_hist2=modelmhist+((ptrdiff_t)kc<<modeldepth);
				for(int ks=0;ks<nlevels;++ks)//generate distribution
				{
					double fval=(ks-mean)/sdev;
					fval=exp(-fval*fval);
					curr_hist2[ks]=(int)(fval*0x10000);
				}
				int msum=0;
				for(int ks=0;ks<nlevels;++ks)
					msum+=curr_hist2[ks];
				norm=1./msum;
				e=0;
				for(int ks=0;ks<nlevels;++ks)//calc cross-entropy
				{
					int freq=curr_hist[ks];
					int prob=curr_hist2[ks];
					if(prob)
						e-=freq*log2(prob*norm);
				}
				modelmsizes[kc]=e/8;


				//if(sum>nlevels*2)
				{
					double hsize=0;
					int cdfW=0;
					int sum2=0;
					const int probbits=12;
					int codelen=probbits+1, CDFlevels=1<<probbits;
					for(int ks=0, ks2=0;ks<nlevels;++ks)//calc model stat overhead
					{
						int sym=((ks>>1^-(ks&1))+half)&mask;//midpoint->zigzag->edges
						int freq=curr_hist[sym];
						//if(freq>sum)
						//	LOG_ERROR("freq %04X sum %04X", freq, sum);
						int cdf=sum2*((1ULL<<probbits)-count)/sum+ks2;
						ks2+=freq!=0;
						//if((unsigned)cdf>=0x10000)
						//	LOG_ERROR("cdf %04X", cdf);
						//if(ks)
						//{
						//	int csym=cdf-cdfW;
						//	if(csym>0)
						//		hsize+=log2(csym);
						//}
						//int csym=cdf-(2*cdfW-cdfWW);
						int csym=cdf-cdfW;
						//csym=csym<<1^csym>>31;
						//if(csym<0)
						//	LOG_ERROR("0x%04X -> 0x%04X", cdfW, cdf);
						if(ks&&CDFlevels)//CDF[0]=0
						{
							//GR				~8.5 KB
							int nbypass=FLOOR_LOG2(CDFlevels);
							if(ks>1)
								nbypass-=7;
							if(nbypass<0)
								nbypass=0;
							hsize+=(csym>>nbypass)+1+nbypass;

							//variable-base code		~25.9 KB
							//if(codelen>1)
							//{
							//	int codelen0=codelen-1;
							//	codelen-=(CDFlevels/((1<<codelen0)-1)+1)*codelen0 < (CDFlevels/((1<<codelen)-1)+1)*codelen;
							//}
							//hsize+=(csym/((1<<codelen)-1)+1)*codelen;

							//hsize+=csym+1;//unary code	~28.7 KB

							//variable-base code v2		39.3 KB
							//if(codelen>1)
							//{
							//	int codelen0=codelen>>1;
							//	codelen-=(CDFlevels/((1<<codelen0)-1)+1)*codelen0 < (CDFlevels/((1<<codelen)-1)+1)*codelen;
							//}
							//hsize+=(csym/((1<<codelen)-1)+1)*codelen;

							//naive pack			~41.5 KB	~38.2 KB
							//hsize+=probbits;

							//hsize+=log2(csym+1)+1;//	~7.2 KB		//X  no stop bit in binary code
						}
						CDFlevels-=csym;
						//if(hsize!=hsize)
						//	LOG_ERROR("sum %d sum2 %d ks %d nlevels %d  0x%04X -> 0x%04X", sum, sum2, ks, nlevels, cdfW, cdf);
						//	LOG_ERROR("ctx %d sum %d sum2 %d ks %d nlevels %d  0x%04X -> 0x%04X", kc, sum, sum2, ks, nlevels, cdfW, cdf);

						cdfW=cdf;
						sum2+=freq;
					}
					modelstatoverhead+=hsize/8;
				}
			}
		}
		break;
	case VIS_JOINT_HISTOGRAM:
		chart_jointhist_update(im1, txid_jointhist);
		break;
	case VIS_ANALYSIS:
		analysis_update(im1);
		break;
	}
}

static void draw_AAcuboid_wire(float size, int applyRCT)
{
	float
		x1=0, x2=size,
		y1=0, y2=size,
		z1=0, z2=size;
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
	if(applyRCT)
	{
		ALIGN(4) static char buf[148]={0};
		static Image *image=0;
		float half=size*0.5f;
		if(!image)
		{
			image=(Image*)buf;
			image->iw=8;
			image->ih=1;
			memset(image->src_depth, 8, sizeof(image->src_depth));
		}
		for(int k=0;k<8;++k)
		{
			image->data[k<<2|0]=(int)(cuboid[k*3+0]-half);
			image->data[k<<2|1]=(int)(cuboid[k*3+1]-half);
			image->data[k<<2|2]=(int)(cuboid[k*3+2]-half);
		}
		memcpy(image->depth, image->src_depth, sizeof(image->depth));
		apply_selected_transforms(&image, 1, 1, 1);
		for(int k=0;k<8;++k)
		{
			cuboid[k*3+0]=(float)image->data[k<<2|0]+half;
			cuboid[k*3+1]=(float)image->data[k<<2|1]+half;
			cuboid[k*3+2]=(float)image->data[k<<2|2]+half;
		}
		//float half=size*0.5f;
		//unsigned char permutation[4]={0};
		//rct_custom_unpackpermutation(rct_custom_params[8], permutation);
		//for(int k=0;k<_countof(cuboid)-2;k+=3)
		//{
		//	float vrtx[]=
		//	{
		//		cuboid[k+permutation[0]]-half,
		//		cuboid[k+permutation[1]]-half,
		//		cuboid[k+permutation[2]]-half,
		//	};
		//	vrtx[0]+=(rct_custom_params[0]*vrtx[1]+rct_custom_params[1]*vrtx[2])*(1/4096.f);
		//	vrtx[1]+=(rct_custom_params[2]*vrtx[0]+rct_custom_params[3]*vrtx[2])*(1/4096.f);
		//	vrtx[2]+=(rct_custom_params[4]*vrtx[0]+rct_custom_params[5]*vrtx[1])*(1/4096.f);
		//	vrtx[1]+=(rct_custom_params[6]*vrtx[0]+rct_custom_params[7]*vrtx[2])*(1/4096.f);
		//	cuboid[k+0]=vrtx[0]+half;
		//	cuboid[k+1]=vrtx[1]+half;
		//	cuboid[k+2]=vrtx[2]+half;
		//}
	}
	draw_3d_line(&cam, cuboid+(0+0)*3, cuboid+(0+1)*3, 0xFF0000FF);
	draw_3d_line(&cam, cuboid+(0+1)*3, cuboid+(0+2)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(0+2)*3, cuboid+(0+3)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(0+3)*3, cuboid+(0+0)*3, 0xFFFF0000);
		
	draw_3d_line(&cam, cuboid+(4+0)*3, cuboid+(4+1)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(4+1)*3, cuboid+(4+2)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(4+2)*3, cuboid+(4+3)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(4+3)*3, cuboid+(4+0)*3, 0xFF000000);
		
	draw_3d_line(&cam, cuboid+(0+0)*3, cuboid+(4+0)*3, 0xFF00FF00);
	draw_3d_line(&cam, cuboid+(0+1)*3, cuboid+(4+1)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(0+2)*3, cuboid+(4+2)*3, 0xFF000000);
	draw_3d_line(&cam, cuboid+(0+3)*3, cuboid+(4+3)*3, 0xFF000000);
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
			//int k=1;
			//float y=k*histpx*10000;
			//for(;y<dy;++k)
			//{
			//	draw_line(x1, y1+(kc+1)*dy-y, x2, y1+(kc+1)*dy-y, color?color:alpha<<24|0xFF<<(kc<<3));//0x40
			//	y=k*histpx*10000;
			//}
			for(int k2=0;k2<256;++k2)
			{
				draw_rect(x1+k2*(x2-x1)/256, x1+(k2+1)*(x2-x1)/256, y1+(kc+1)*dy-_hist[kc<<8|k2]*histpx, y1+(kc+1)*dy, color?color:alpha<<24|0xFF<<(kc<<3));//0x80
				float x=x1+(k2+0.5f)*(x2-x1)/256, y;
				if(k2<128)
					y=(float)_hist[kc<<8|(k2+0)]/_hist[kc<<8|(k2+1)];
				else if(k2>128)
					y=(float)_hist[kc<<8|(k2+0)]/_hist[kc<<8|(k2-1)];
				else
					y=0;
				draw_line(x, y1+(kc+1)*dy, x, y1+(kc+1)*dy-y*dy*0.25f, color?color:0xFF<<24|0xFF<<(kc<<3));
			}
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
	{
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
	draw_AAcuboid_wire(jh_cubesize, 0);
	draw_AAcuboid_wire(jh_cubesize, 1);

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
		display_texture_i(0, wndw, 0, wndh, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 0.5, 0);
		
		{
			int crosshaircolor=0xFF000000;
			float ratios[]={(float)wndw/im1->iw, (float)wndh/im1->ih};
			float sbounds[4];
			for(int k=0;k<4;++k)
				sbounds[k]=bounds[k]*ratios[k>>1];
			{
				float xmid=(sbounds[0]+sbounds[1])*0.5f;
				float ymid=(sbounds[2]+sbounds[3])*0.5f;
				draw_line(sbounds[0], sbounds[2], sbounds[0], sbounds[3], crosshaircolor);//{x1, y1, x2, y2}
				draw_line(sbounds[1], sbounds[2], sbounds[1], sbounds[3], crosshaircolor);
				draw_line(sbounds[0], sbounds[2], sbounds[1], sbounds[2], crosshaircolor);
				draw_line(sbounds[0], sbounds[3], sbounds[1], sbounds[3], crosshaircolor);
				draw_line(xmid, 0, xmid, (float)wndh, crosshaircolor);
				draw_line(0, ymid, (float)wndw, ymid, crosshaircolor);
			}
		}
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
AABB buttons[7]={0};//0: CT,  1: ST,  2: list of transforms,  3: jxl params,  4: EC,  5: OLS-4

int io_init(int argc, char **argv)//return false to abort
{
	(void)argc;//FIXME open 1 command argument
	(void)argv;

	if(tdy*T_COUNT/2*g_uiscale>wndh*1.2f)
		g_uiscale=wndh*1.2f/(tdy*T_COUNT/2);

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
	AABB *p=buttons;//TODO use enum
	float xstep=tdx*guizoom, ystep=tdy*guizoom;
	p->x1=xstep*2, p->x2=p->x1+xstep*gui_custom_rct_w, p->y1=(float)(wndh>>1), p->y2=p->y1+ystep*gui_custom_rct_h, ++p;//0: color params - left
	p->x1=(float)(wndw>>3), p->x2=p->x1+xstep*gui_custom_pred_w, p->y1=(float)((wndh>>1)+(wndh>>2))-ystep, p->y2=p->y1+ystep*gui_custom_pred_h, ++p;//1: spatial params - bottom
	p->x1=(float)(wndw-300), p->x2=(float)wndw, p->y1=tdy*2, p->y2=p->y1+tdy*T_COUNT/2*g_uiscale, ++p;//2: transforms list
	
	//p->x1=(float)(w>>1), p->x2=p->x1+xstep*14, p->y1=(float)((h>>1)+(h>>2))+ystep*4, p->y2=p->y1+ystep, ++p;//3: clamp bounds
	//p->x1=(float)(w>>1), p->x2=p->x1+xstep*21, p->y1=(float)((h>>1)+(h>>2))+ystep*5, p->y2=p->y1+ystep, ++p;//4: learning rate

	p->x1=(float)(wndw>>2), p->x2=p->x1+tdx*11*6, p->y1=(float)((wndh>>1)+(wndh>>2)), p->y2=p->y1+tdy*3, ++p;//3: jxl params

	p->x1=(float)(wndw-450), p->x2=p->x1+tdx*gui_ec_width, p->y1=tdy, p->y2=p->y1+tdy, ++p;//4: EC method	//H.E.M.L..A.0x0000..XXXX_XXX

	p->x1=(float)(wndw>>2), p->x2=p->x1+xstep*(OLS4_RMAX<<1|1)*gui_ols4_elementchars, p->y1=(float)(wndh>>1)+10, p->y2=p->y1+ystep*(1+(OLS4_RMAX+1)*(im1?im1->nch:4)), ++p;//5: OLS-4		DON'T USE p->y2

	p->x1=(float)(wndw>>2), p->x2=p->x1+xstep*60, p->y1=(float)(wndh>>1)+10, p->y2=p->y1+ystep*8, ++p;//6: CUSTOM4
	
	if(im1&&im1->iw&&im1->ih&&imagecentered)
		center_image();
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

				step=jhc_boxdx*wndw/im1->iw;
				jhc_xbox-=jhc_xbox%step-(step>>1);
				step=jhc_boxdy*wndh/im1->ih;
				jhc_ybox-=jhc_ybox%step-(step>>1);
			}
			chart_jointhist_update(im1, txid_jointhist);
		}
		else
		{
			int X0=wndw>>1, Y0=wndh>>1;
			if(mode==VIS_JOINT_HISTOGRAM)//drag cam
			{
				cam_turnMouse(cam, mx-X0, my-Y0, mouse_sensitivity);
				set_mouse(X0, Y0);
			}
			else if(mode==VIS_ANALYSIS)//drag analysis box
			{
				jhc_xbox=mx;
				jhc_ybox=my;
				analysis_update(im1);
			}
			else//drag image
			{
				wpx-=(mx-X0)/imzoom;
				wpy-=(my-Y0)/imzoom;
				set_mouse(X0, Y0);
			}
		}
		return !timer;
	}
	else if(profileplotmode>PROFILE_OFF)
		return 1;
	return 0;
}
static void click_hittest(int _mx, int _my, int *objidx, int *cellx, int *celly, int *cellidx, AABB **p)
{
	*p=buttons;
	*objidx=0;
	for(*objidx=0;*objidx<(int)_countof(buttons);++*objidx, ++*p)
	{
		//when these buttons are inactive they shouldn't block the click
		if(
			(!transforms_customenabled&&(*objidx==0||*objidx==1))||
			(!(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])&&*objidx==3)||
			(!(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])&&*objidx==5)||
			(!(transforms_mask[ST_FWD_CUSTOM4]||transforms_mask[ST_INV_CUSTOM4]||transforms_mask[ST_FWD_WC]||transforms_mask[ST_INV_WC])&&*objidx==6)
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
	case 6:
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
		if(GET_KEY_STATE(KEY_CTRL))
		{
			g_uiscale+=(2*(forward>0)-1)*0.01f;
			CLAMP2(g_uiscale, 0.01f, 4);
			io_resize();
			return 1;
		}
		int objidx=0, cellx=0, celly=0, cellidx=0;
		AABB *gui_cell=buttons;
		click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &gui_cell);
		if(objidx!=-1)
		{
			int sign=(forward>0)-(forward<0);//abs(forward) is 120
			int ch=(int)floorf((mx-gui_cell->x1)/(guizoom*tdx));
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
						x=(x+sign+6)%6;
						//MODVAR(x, x, 6);
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
						//col=gui_cell?(int)floorf((mx-gui_cell->x1)/(guizoom*tdx)):0,
						line=gui_cell?(int)floorf((my-gui_cell->y1)/(guizoom*tdy)):0;
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
				if(!transforms_mask[ST_FWD_CUSTOM]&&!transforms_mask[ST_INV_CUSTOM]
					&&!transforms_mask[ST_CONVTEST]
					&&!transforms_mask[ST_CONVTEST2]
					&&!transforms_mask[ST_FWD_CC]&&!transforms_mask[ST_INV_CC]
				//	&&!transforms_mask[ST_FWD_OLS6]&&!transforms_mask[ST_INV_OLS6]
				)
					break;
				if(!celly)
				{
					//0123456789012345678901234567
					//Ch 0  Clamp [+W +NW +N +NE]
					if(ch==3)
						custom_pred_ch_idx=(custom_pred_ch_idx+sign+3)%3;
					//{
					//	custom_pred_ch_idx+=sign;
					//	MODVAR(custom_pred_ch_idx, custom_pred_ch_idx, 3);
					//}
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
						{
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
					ch=(int)floorf((mx-gui_cell->x1)/tdx);
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
					{
						int method0=ec_method;
						int *gcc_happy=(int*)&ec_method;
						*gcc_happy+=sign;
						MODVAR(*gcc_happy, *gcc_happy, ECTX_COUNT);

						static int prevmode=0;
						if(ec_method==ECTX_ABAC0)
						{
							prevmode=viewmode;
							viewmode=VIEW_C0;
						}
						else if(method0==ECTX_ABAC0)
							viewmode=prevmode;
					}
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
						{
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
									double val;

									digit+=digit<=-1;
									val=ols4_lr[kx-1];
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
			case 6://CUSTOM4
				if(!celly)
				{
					if(cellx==6)//change layer, stay in same channel
					{
						lossyconv_page+=sign<<2;
						lossyconv_page&=15;
					}
					else if(cellx==13)//change channel, stay in same layer
					{
						lossyconv_page=(lossyconv_page&12)|((lossyconv_page+sign)&3);
					}
					else if(cellx==23)
					{
						unsigned char *cell=lossyconv_causalRCT+(lossyconv_page>>2);
						*cell=!*cell;
						update_image();
					}
				}
				else if(celly==6)
				{
					unsigned char *cell;
					if((unsigned)(cellx-11)<3)
					{
						cell=lossyconv_offset+((size_t)lossyconv_page>>2<<1|0);
						*cell+=sign;
						update_image();
					}
					else if((unsigned)(cellx-15)<3)
					{
						cell=lossyconv_offset+((size_t)lossyconv_page>>2<<1|1);
						*cell+=sign;
						update_image();
					}
				}
				else if(celly==7)
				{
					unsigned char *cell;
					if((unsigned)(cellx-11)<3)
					{
						cell=lossyconv_stride+((size_t)lossyconv_page>>2<<1|0);
						*cell+=sign;
						update_image();
					}
					else if((unsigned)(cellx-15)<3)
					{
						cell=lossyconv_stride+((size_t)lossyconv_page>>2<<1|1);
						*cell+=sign;
						update_image();
					}
				}
				else if((unsigned)(celly-1)<5)
				{
					//    000000000011111111112222222222333333
					//    012345678901234567890123456789012345
					//0  "Layer 1/4  C 1/3  cRCT 1"
					//    0000000000111111111122222222223333333333444444444455555555556
					//    0123456789012345678901234567890123456789012345678901234567890
					//1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
					//4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//6  "XY  Offset XXX YYY"
					//7  "    Stride XXX YYY"
					int c=cellx%12, x=cellx/12;
					short *cell=lossyconv_params+5*5*lossyconv_page+5*(celly-1)+x;
					if(c==10)//toggle clamp mask
					{
						*cell^=1;
						update_image();
					}
					else if((unsigned)(c-1)<9)
					{
						//123456789
						//876.54321
						int bitidx=9-c;
						bitidx+=bitidx<6;
						*cell=((unsigned short)(*cell+(sign<<bitidx))&~1)|(*cell&1);
						*cell<<=16-9;//9-bit params
						*cell>>=16-9;
						update_image();
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
		if(mode==VIS_JOINT_HISTOGRAM)
		{
			if(show_full_image)
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
					jhc_boxdx=CLAMP(1, jhc_boxdx, 256);
					jhc_boxdy=CLAMP(1, jhc_boxdy, 256);
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
		else if(mode==VIS_ANALYSIS)
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
			jhc_boxdx=CLAMP(8, jhc_boxdx, 1024);
			jhc_boxdy=CLAMP(8, jhc_boxdy, 1024);
			analysis_update(im1);
		}
		else
		{
			int mw_fwd=forward>0;
			zoom_at(mx, my, mw_fwd?2:1/2.);//fwd zooms in
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
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		end=(char*)text->data+idx;
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
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		//if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
		//	idx+=2;
		end=(char*)text->data+idx;
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
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		end=(char*)text->data+idx;
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
static int parse_nvals_i16(ArrayHandle text, int idx, short *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(text->data[idx]);++idx);

		neg=text->data[idx]=='-';
		idx+=neg||text->data[idx]=='+';//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		end=(char*)text->data+idx;
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
	case KEY_UP:
	case KEY_DOWN:
		if(im1)
		{
			g_dist+=(key==KEY_UP)*2-1;
			if(g_dist<1)
				g_dist=1;
			if(g_dist>16)
				g_dist=16;
			update_image();
			return 1;
		}
		break;
	case KEY_LEFT:
	case KEY_RIGHT:
		if(mode!=VIS_JOINT_HISTOGRAM)
		{
			if(fn)
			{
				const char *ext[]=
				{
					"PPM", "PGM", "PNM",
					"PNG",
					"JPG", "JPEG",
					"BMP",
					"TIF", "TIFF",
				};
				ArrayHandle path, filenames, filteredfn, *fn2;

				STR_COPY(path, fn->data, fn->count);
				//acme_strrchr((char*)path->data, path->count, '/');//X  what about backslash?
				for(int k=(int)path->count-1;k>=0;--k)
				{
					if(path->data[k]=='/'||path->data[k]=='\\')
					{
						path->data[k+1]=0;
						path->count=(size_t)k+1;
						break;
					}
				}
				filenames=get_filenames((char*)path->data, ext, _countof(ext), 1);
				filteredfn=filter_path((char*)fn->data, (int)fn->count, 0);
				if(filenames&&filenames->count)
				{
					int currentidx=-1;
					for(int k=0;k<(int)filenames->count;++k)
					{
						fn2=(ArrayHandle*)array_at(&filenames, k);
						if(!_stricmp((char*)fn2[0]->data, (char*)filteredfn->data))
						{
							currentidx=k;
							break;
						}
					}
					{
						Image *im2=0;
						int step=key==KEY_RIGHT?1:-1;
						for(int k=currentidx+step;MODVAR(k, k, (int)filenames->count), k!=currentidx;k+=step)
						{
							fn2=(ArrayHandle*)array_at(&filenames, k);
							im2=image_load((char*)fn2[0]->data, (int)fn2[0]->count);
							//if(!load_media(fn2[0]->data, &im2, 0))
							if(im2)
							{
								currentidx=k;
								break;
							}
						}
						if(im2)
						{
							free(im0);
							im0=im2;
							array_free(&fn);
							filesize=get_filesize((char*)fn2[0]->data);
							fn=filter_path((char*)fn2[0]->data, (int)fn2[0]->count, 0);
							update_image();
							set_window_title("%s - eBench", (char*)fn->data);
						}
					}
					array_free(&filenames);
				}
				array_free(&filteredfn);
				array_free(&path);
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
			AABB *guicell=buttons;
			click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &guicell);
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
						col=guicell?(int)floorf((mx-guicell->x1)/(guizoom*tdx)):0,
						line=guicell?(int)floorf((my-(guicell->y1+guizoom*tdy))/(guizoom*tdy)):0;
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
						col=guicell?(int)floorf((mx-guicell->x1)/(guizoom*tdx)):0,
						line=guicell?(int)floorf((my-(guicell->y1+guizoom*tdy))/(guizoom*tdy)):0;
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
						line=celly;
					//	col=guicell?(int)floorf((mx-guicell->x1)/tdx):0,
					//	line=guicell?(int)floorf((my-guicell->y1)/tdy):0;
					//if(BETWEEN_EXC(0, cellidx, T_COUNT/2))
					if((unsigned)line<(unsigned)(T_COUNT/2))
					{
						int idx=line<<1|(guicell?(int)floor((mx-guicell->x1)*2/(guicell->x2-guicell->x1)):0);
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
					ec_method=ECTX_STATIC_O0;
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
						{
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
			case 6://CUSTOM4
				if((unsigned)(celly-1)<5)
				{
					//    000000000011111111112222222222333333
					//    012345678901234567890123456789012345
					//0  "Layer 1/4  C 1/3  cRCT 1"
					//    0000000000111111111122222222223333333333444444444455555555556
					//    0123456789012345678901234567890123456789012345678901234567890
					//1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
					//4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					//6  "XY  Offset XXX YYY"
					//7  "    Stride XXX YYY"
					short *ptr=lossyconv_params+5*5*lossyconv_page+5*(celly-1)+cellx/12;
					if(key==KEY_LBUTTON)
					{
						*ptr=lossyconv_clipboard;
						if(!GET_KEY_STATE(KEY_CTRL))
							update_image();
					}
					else
						lossyconv_clipboard=*ptr;
					return 1;
				}
				else if(celly==6)
				{
					if((unsigned)(cellx-11)<3)
					{
						lossyconv_offset[lossyconv_page>>2<<1|0]=0;
						update_image();
						return 1;
					}
					if((unsigned)(cellx-15)<3)
					{
						lossyconv_offset[lossyconv_page>>2<<1|1]=0;
						update_image();
						return 1;
					}
				}
				else if(celly==7)
				{
					if((unsigned)(cellx-11)<3)
					{
						lossyconv_stride[lossyconv_page>>2<<1|0]=0;
						update_image();
						return 1;
					}
					if((unsigned)(cellx-15)<3)
					{
						lossyconv_stride[lossyconv_page>>2<<1|1]=0;
						update_image();
						return 1;
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
			if(mode==VIS_ANALYSIS)
			{
				jhc_xbox=mx;
				jhc_ybox=my;
				analysis_update(im1);
				return 1;
			}
			if(drag)//enter mouse control
			{
				mx0=mx, my0=my;
				set_mouse(wndw>>1, wndh>>1);
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
			"E:\t\tShow image at 1:1 scale\n"
			"C:\t\tCenter image\n"
			"J / Shift J:\tToggle single-channel view\n"
			"N:\t\tToggle modular arithmetic in image view\n"
			"Ctrl R:\t\tDisable all transforms\n"
			"Ctrl Space:\tMeasure PSNR\n"
		//	"Ctrl E:\t\tReset custom transform parameters\n"
			"Comma/Period:\tSelect context for size estimation\n"
			"Slash:\t\tToggle adaptive histogram\n"
			"[ ]:\t\t(Custom transforms) Select coefficient page\n"
			"Ctrl 1/2/3:\t(Custom predictor) Populate page from channel N\n"
			"Space:\t\t(Custom transforms) Optimize\n"
			"Ctrl N:\t\tAdd noise to CUSTOM3 params\n"
		//	"Shift space:\t(Custom transforms) Optimize blockwise\n"
		//	"Ctrl Space\t(Custom transforms) Reset params\n"
			"Ctrl B:\t\tBatch test\n"
			"Ctrl L:\t\tLossy batch test\n"
			"\n"
		//	"WASDTG:\tMove cam\n"
		//	"Arrow keys:\tTurn cam\n"
		//	"Mouse1/Esc:\tToggle mouse look\n"
		//	"R:\t\tReset cam\n"
		//	"\n"
			"H:\t\tReset invCR history graph\n"
			"Ctrl C:\t\tCopy data\n"
			"Ctrl V:\t\tPaste data\n"
			"Ctrl B:\t\tBatch test\n"
		//	"Ctrl SPACE:\tCheck integrity (if restored bit-exact)\n"
		//	"Ctrl P:\t\tTest predictors\n"
		//	"C:\t\tToggle joint histogram type / fill screen in image view\n"
			"\n"
			"M / Shift M:\tCycles between:\n"
		//	"\t1: 3D View: Levels\n"
		//	"\t2: 3D View: Mesh\n"
		//	"\t3: 3D View: Mesh (separate channels)\n"
			"\t%d: Image tricolor view\n"
#ifdef ENABLE_L1WEIGHTS
			"\t%d: L1 Weights\n"
#endif
			"\t%d: Image view\n"
		//	"\t6: Image block histogram\n"
		//	"\t7: Optimized block compression estimate (E24)\n"
		//	"\t7: DWT block histogram\n"
			"\t%d: Histogram\n"
			"\t%d: Model\n"
			"\t%d: Analysis\n"
			"\t%d: Joint histogram\n"
			"\t%d: Zipf view\n",

			1+VIS_IMAGE_TRICOLOR,
#ifdef ENABLE_L1WEIGHTS
			1+VIS_L1WEIGHTS,
#endif
			1+VIS_IMAGE,
			1+VIS_HISTOGRAM,
			1+VIS_MODEL,
			1+VIS_ANALYSIS,
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
			{
				int *gcc_happy=(int*)&ec_method;
				*gcc_happy+=fwd;
				MODVAR(*gcc_happy, *gcc_happy, ECTX_COUNT);
				update_image();
			}
		}
		return 1;
	case KEY_SLASH:
		ec_adaptive=!ec_adaptive;
		//ec_adaptive=(ec_adaptive+1)%3;
		update_image();
		return 1;
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
	//case 'E':
	//	if(im1&&GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)//reset params
	//	{
	//		customtransforms_resetparams();
	//		update_image();
	//		return 1;
	//	}
	//	wireframe=!wireframe;
	//	break;
	case 'E':
		wpx=0, wpy=0, imzoom=1;
		imagecentered=0;
		return 1;
	case 'C':
		if(im1&&GET_KEY_STATE(KEY_CTRL))
		{
			ArrayHandle str;

			if(GET_KEY_STATE(KEY_SHIFT))
			{
				unsigned char *buf=0;
				image_export_uint8(im1, &buf, 1, 1);//swap red & blue for WinAPI
				if(!buf)
				{
					LOG_WARNING("Alloc error");
					return 0;
				}
				{
					int success=copy_bmp_to_clipboard(buf, im1->iw, im1->ih);
					if(!success)
						LOG_WARNING("Failed to copy image to clipboard");
				}
				free(buf);
				return 0;
			}
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
			}//no 'else' here
			if(mode==VIS_IMAGE||mode==VIS_ZIPF)//copy custom transform value
			{
				double invcr[5]={0}, csizes[5]={0}, usize;
				entropy2invcr(ch_entropy, im0->src_depth, 4, invcr);
				usize=invcr2csizes(invcr, im0->src_depth, im0->iw, im0->ih, 4, csizes);
				str_append(&str,
					"UTYUVA %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf bytes  BPD %8.4lf",
					usize,
					csizes[0],
					csizes[1],
					csizes[2],
					csizes[3],
					csizes[4],
					csizes[0]/usize*8
				);

				//str_append(&str,
				//	"TRGBA %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf%%  TRGBA %12.2lf %12.2lf %12.2lf %12.2lf %12.2lf bytes",
				//	100.*invcr[0],
				//	100.*invcr[1],
				//	100.*invcr[2],
				//	100.*invcr[3],
				//	100.*invcr[4],
				//	csizes[0],
				//	csizes[1],
				//	csizes[2],
				//	csizes[3],
				//	csizes[4]
				//);

				//str_append(&str, "TRGBA %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf",
				//	100.*(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3])/(im1->src_depth[0]+im1->src_depth[1]+im1->src_depth[2]+im1->src_depth[3]),
				//	100.*ch_entropy[0]/im1->src_depth[0],
				//	100.*ch_entropy[1]/im1->src_depth[1],
				//	100.*ch_entropy[2]/im1->src_depth[2],
				//	100.*ch_entropy[3]/im1->src_depth[1]
				//);
			}
#if 0
			//else if(mode==VIS_IMAGE_E24)
			//{
			//	for(int kc=0;kc<3;++kc)
			//	{
			//		E24Params const *ptr=e24_params+kc;
			//		str_append(&str, "%3d  %3d %3d %3d  %3d %3d\n", ptr->gwidth, ptr->mleft, ptr->mtop, ptr->mright, ptr->alpha, ptr->maxinc);
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
				const int stw=(custom_pred_reach<<1|1)*2;
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
				for(int kc2=0;kc2<3;++kc2)
				{
					const int np=(int)(_countof(custom_params)/3);
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
				//for(int k=0;k<(int)_countof(customparam_clamp);++k)
				//	str_append(&str, "%d%c", customparam_clamp[k], k<(int)_countof(customparam_clamp)-1?'\t':'\n');
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
			{
				center_image();
				show_full_image=!show_full_image;
			}

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
					parse_nvals_i16(text, idx, logic_params, (int)_countof(logic_params));
				}
				else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
				{
					parse_nvals_i16(text, idx, (short*)&c2_params, sizeof(c2_params)/sizeof(short));
				}
				else
#endif
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					idx=parse_nvals_i32(text, idx, &ols4_period, 1);
					idx=parse_nvals_f64(text, idx, ols4_lr, (int)_countof(ols4_lr));
					idx=parse_nvals_i8(text, idx, (unsigned char*)ols4_mask, sizeof(ols4_mask)/sizeof(char));
				}
				if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
				{
					parse_nvals_i16(text, idx, (short*)&c3_params, sizeof(c3_params)/sizeof(short));
				}
				else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					parse_nvals_i16(text, idx, jxlparams_i16, (int)_countof(jxlparams_i16));
				}
				else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
				{
					parse_nvals_i16(text, idx, pw2_params, (int)_countof(pw2_params));
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
					k=0, kend=(int)_countof(rct_custom_params);
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
					kend=(int)_countof(custom_params);
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
						idx=parse_nvals_i32(text, idx, custom_params, (int)_countof(custom_params));
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
	case KEY_4:
		if(GET_KEY_STATE(KEY_CTRL))
		{
			if(transforms_customenabled)
			{
				int srcidx=key-KEY_1;
				if(srcidx<3&&srcidx!=custom_pred_ch_idx)
				{
					memcpy(custom_params+(_countof(custom_params)/3)*custom_pred_ch_idx, custom_params+(_countof(custom_params)/3)*srcidx, sizeof(custom_params)/3);
					update_image();
					return 1;
				}
			}
			else if(transforms_mask[ST_FWD_CUSTOM4]||transforms_mask[ST_INV_CUSTOM4]||transforms_mask[ST_FWD_WC]||transforms_mask[ST_INV_WC])
			{
				int srcidx=(lossyconv_page&0xF0)|((key-KEY_1)&15);
				if(srcidx!=lossyconv_page)
				{
					memcpy(lossyconv_params+5*5*lossyconv_page, lossyconv_params+5*5*srcidx, sizeof(short[25]));
					update_image();
					return 1;
				}
			}
		}
		break;
	case 'M':
		if(im1)
		{
			int shift=GET_KEY_STATE(KEY_SHIFT);
			int m0=mode;
			mode+=1-2*shift+VIS_COUNT;
			mode%=VIS_COUNT;
			//int m0=mode;
			//MODVAR(mode, mode+1-(shift<<1), VIS_COUNT);
			update_image();
#ifdef ENABLE_L1WEIGHTS
			if(m0==VIS_L1WEIGHTS)
				timer_stop(11);
			else if(mode==VIS_L1WEIGHTS)
				timer_start(L1DELAY, 11);
#endif
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
	case 'X'://toggle horizontal profile plot
		if(profileplotmode!=PROFILE_X)
			profileplotmode=PROFILE_X;
		else
			profileplotmode=PROFILE_OFF;
		return 1;
	case 'Y'://toggle vertical profile plot
		if(profileplotmode!=PROFILE_Y)
			profileplotmode=PROFILE_Y;
		else
			profileplotmode=PROFILE_OFF;
		return 1;
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
	case 'J':
		if(im0)
		{
			int shift=GET_KEY_STATE(KEY_SHIFT);
			viewmode+=shift?-1:1;
			MODVAR(viewmode, viewmode, VIEW_COUNT);
			update_image();
			return 1;
		}
		break;
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
		else if(mode==VIS_L1WEIGHTS)
		{
			int fwd=!GET_KEY_STATE(KEY_SHIFT);
			fwd=2*fwd-1;
			use_ols=(use_ols+fwd+3)%3;
		}
		break;
	case 'L':
		if(GET_KEY_STATE(KEY_CTRL))
		{
			batch_test_lossy();
		}
		break;
	case KEY_LBRACKET:
	case KEY_RBRACKET:
		if(ec_method==ECTX_ABAC0)
		{
			if(im1)
			{
				if(abacvis_range<=2)
					break;
				//{
				//	abacvis_low=0;
				//	abacvis_range=im1->depth[0];
				//	if(abacvis_range<im1->depth[1])abacvis_range=im1->depth[1];
				//	if(abacvis_range<im1->depth[2])abacvis_range=im1->depth[2];
				//	if(abacvis_range<im1->depth[3])abacvis_range=im1->depth[3];
				//	abacvis_range=1<<abacvis_range;
				//}
				//else
				{
					int floorhalf=abacvis_range>>1;
					int bit=key==KEY_RBRACKET;
					if(bit)
						abacvis_low+=abacvis_range-floorhalf;
					abacvis_range=floorhalf;
				}
				update_image();
				return 1;
			}
		}
		else if(transforms_customenabled)
		{
			custom_pred_ch_idx+=((key==KEY_RBRACKET)<<1)-1;
			MODVAR(custom_pred_ch_idx, custom_pred_ch_idx, 3);
			return 1;
		}
		else if(transforms_mask[ST_FWD_CUSTOM4]||transforms_mask[ST_INV_CUSTOM4]||transforms_mask[ST_FWD_WC]||transforms_mask[ST_INV_WC])
		{
			int sign=((key==KEY_RBRACKET)<<1)-1;
			if(GET_KEY_STATE(KEY_SHIFT))//change layer
			{
				lossyconv_page+=sign<<2;
				lossyconv_page&=15;
			}
			else//change channel
				lossyconv_page=(lossyconv_page&12)|((lossyconv_page+sign)&3);
			return 1;

		}
		break;
	case KEY_QUOTE:
		if(im1&&ec_method==ECTX_ABAC0)
		{
			int maxdepth=im1->depth[0];
			if(maxdepth<im1->depth[1])maxdepth=im1->depth[1];
			if(maxdepth<im1->depth[2])maxdepth=im1->depth[2];
			if(maxdepth<im1->depth[3])maxdepth=im1->depth[3];
			if(abacvis_range<(1<<maxdepth))
			{
				abacvis_low&=~abacvis_range;
				abacvis_range<<=1;
				update_image();
				return 1;
			}
		}
		break;
	case 'N':
		if(GET_KEY_STATE(KEY_CTRL))//Ctrl N	add noise to params
		{
			if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				srand((unsigned)__rdtsc());
				{
					short *ptr=(short*)&c3_params;
					for(int k=0;k<C3_NPARAMS;++k)
						ptr[k]+=rand()%3-1;
				}
				update_image();
			}
			return 1;
		}
		else if(mode==VIS_L1WEIGHTS)
		{
			for(int k=0;k<_countof(l1weights);++k)
				l1weights[k]+=rand()-RAND_MAX/2;
			return 1;
		}
		//{
		//	view_ma=!view_ma;
		//	update_image();
		//	return 1;
		//}
		break;
	case 'Z':
		if(mode==VIS_L1WEIGHTS)
		{
			memset(l1weights, 0, sizeof(l1weights));
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
				//int success=1;
				if(!im1)
					return 0;
				if(mode==VIS_ANALYSIS)
				{
					C18Info *info=(C18Info*)malloc(sizeof(C18Info));
					if(!info)
					{
						LOG_WARNING("Alloc error");
						return 0;
					}

					double esizes[3]={0};
					int nblocks=0;
					for(int ky=0;ky<im1->ih;ky+=jhc_boxdy)
					{
						int ky2=ky+jhc_boxdy;
						if(ky2>im1->ih)
							ky2=im1->ih;
						for(int kx=0;kx<im1->iw;kx+=jhc_boxdx)
						{
							int kx2=kx+jhc_boxdx;
							if(kx2>im1->iw)
								kx2=im1->iw;
							c18_analyze(im1, kx, kx2, ky, ky2, info);
							const unsigned char *group=rct2_combinations[info->bestrct];
							int selected_ch[]=
							{
								group[0]*PRED_COUNT+info->predidx[0],
								group[1]*PRED_COUNT+info->predidx[1],
								group[2]*PRED_COUNT+info->predidx[2],
							};
							esizes[0]+=info->esizes[selected_ch[0]]*(kx2-kx)*(ky2-ky);
							esizes[1]+=info->esizes[selected_ch[1]]*(kx2-kx)*(ky2-ky);
							esizes[2]+=info->esizes[selected_ch[2]]*(kx2-kx)*(ky2-ky);
							++nblocks;
						}
					}
					ptrdiff_t usize=(ptrdiff_t)im1->iw*im1->ih;
					char buf[512]={0};
					int nprinted=snprintf(buf, sizeof(buf)-1,
						"%d*%d block size, %d blocks\n"
						"T %6.2lf%%  %12.2lf bytes\n"
						"Y %6.2lf%%  %12.2lf bytes\n"
						"U %6.2lf%%  %12.2lf bytes\n"
						"V %6.2lf%%  %12.2lf bytes\n",
						jhc_boxdx, jhc_boxdy, nblocks,
						(esizes[0]+esizes[1]+esizes[2])*100./(3*usize), esizes[0]+esizes[1]+esizes[2],
						esizes[0]*100/usize, esizes[0],
						esizes[1]*100/usize, esizes[1],
						esizes[2]*100/usize, esizes[2]
					);
					copy_to_clipboard(buf, nprinted);
					messagebox(MBOX_OK, "Copied to clipboard", "%s", buf);
					free(info);
				}
				else if(im0->iw!=im1->iw||im0->ih!=im1->ih||im0->nch!=im1->nch)
				{
					messagebox(MBOX_OK, "Dimension Error",
						"Image dimensions changed:\n"
						"\tim0\tCWH %d*%d*%d\n"
						"\tim1\tCWH %d*%d*%d",
						im0->nch, im0->iw, im0->ih,
						im1->nch, im1->iw, im1->ih
					);
					//success=0;
				}
				else
				{
					double rmse[4]={0}, psnr[4]={0};
					ptrdiff_t idx=calc_psnr(im0, im1, rmse, psnr);
					int printed=snprintf(g_buf, G_BUF_SIZE-1,
						"RMSE PSNR\n"
						"C0 %12lf %12lf\n"
						"C1 %12lf %12lf\n"
						"C2 %12lf %12lf\n"
						"C3 %12lf %12lf",
						rmse[0], psnr[0],
						rmse[1], psnr[1],
						rmse[2], psnr[2],
						rmse[3], psnr[3]
					);
					if(idx>=0)
					{
						int kx=(int)(idx%im0->iw), ky=(int)(idx/im0->iw);
						printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1,
							"\nError vs Original at XY %d %d:\n"
							"C0\t0x%08X\t0x%08X\n"
							"C1\t0x%08X\t0x%08X\n"
							"C2\t0x%08X\t0x%08X\n"
							"C3\t0x%08X\t0x%08X",
							kx, ky,
							im1->data[idx<<2|0], im0->data[idx<<2|0],
							im1->data[idx<<2|1], im0->data[idx<<2|1],
							im1->data[idx<<2|2], im0->data[idx<<2|2],
							im1->data[idx<<2|3], im0->data[idx<<2|3]
						);
						copy_to_clipboard(g_buf, printed);
					}
					messagebox(MBOX_OK, idx>=0?"Error (Copied to clipboard)":"Bit Exact", "%s", g_buf);
#if 0
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
#endif
				}
				//if(success)
				//	messagebox(MBOX_OK, "SUCCESS", "The image is bit-exact.");
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
				apply_selected_transforms(&im2, 1, 1, 1);
				pred_custom_optimize(im2, custom_params, GET_KEY_STATE(KEY_SHIFT)?2:1);
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
					apply_selected_transforms(&im2, 1, 1, 1);
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
	if(mode==VIS_JOINT_HISTOGRAM)
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
	else if(mode==VIS_ANALYSIS)
	{
		int update=0;
		if(keyboard['W'])	--jhc_ybox, update=1;
		if(keyboard['A'])	--jhc_xbox, update=1;
		if(keyboard['S'])	++jhc_ybox, update=1;
		if(keyboard['D'])	++jhc_xbox, update=1;
		if(im1&&update)
			analysis_update(im1);
	}
	else
	{
		const int delta=10;//screen pixels per frame
		if(keyboard[KEY_ENTER])	zoom_at(wndw>>1, wndh>>1, 1.02);
		if(keyboard[KEY_BKSP])	zoom_at(wndw>>1, wndh>>1, 1/1.02);
		if(keyboard['W'])	wpy-=delta/imzoom;//move window up
		if(keyboard['A'])	wpx-=delta/imzoom;//move window left
		if(keyboard['S'])	wpy+=delta/imzoom;//move window down
		if(keyboard['D'])	wpx+=delta/imzoom;//move window right
	}
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
	for(int k=0;k<(int)_countof(points);k+=3)
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
	for(int k=0;k<(int)_countof(points);k+=3)
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
	for(int k=0;k<(int)_countof(points);k+=3)
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
	for(int k=0;k<(int)_countof(points)-2;k+=3)
	{
		points[k+0]-=points[k+1];
		points[k+1]+=points[k+0]/2;
	}
	draw_3d_line(&cam, points+3*0, points+3*1, 0xFF000000);
	draw_3d_line(&cam, points+3*1, points+3*2, 0xFF0000FF);
	draw_3d_line(&cam, points+3*2, points+3*3, 0xFF00FF00);
	draw_3d_line(&cam, points+3*3, points+3*0, 0xFFFF0000);
}

static ArrayHandle vertices_2d=0;
static void draw_profile_x(int comp, int color)//horizontal cross-section profile		to see the color/spatial correlation
{
	int iy=screen2image_y_int(my);
	if((unsigned)iy<(unsigned)im1->ih)
	{
		int *row=im1->data+(im1->iw*iy<<2|comp);
		int ix;
		float y2;
		float gain=(float)wndh/((1<<im1->depth[comp])-1);
	//	float gain=(float)(wndh>>1)/((1<<im1->depth[comp])-1);
		int offset=1<<im1->depth[comp]>>1;
		for(int kx=0;kx<wndw;++kx)
		{
			ix=screen2image_x_int(kx);
			if((unsigned)ix<(unsigned)im1->iw)
				y2=wndh-tdy-(row[ix<<2]+offset)*gain;
			else
				y2=wndh-tdy;
			draw_curve_enqueue(&vertices_2d, (float)kx, y2);
		}
		draw_2d_flush(vertices_2d, color, GL_LINE_STRIP);
	}
}
static void draw_profile_y(int comp, int color)//vertical cross-section profile
{
	int ix=screen2image_x_int(mx);
	if((unsigned)ix<(unsigned)im1->iw)
	{
		int *col=im1->data+(ix<<2|comp);
		
		float x2;
		int iy;
		float gain=(float)wndw/((1<<im1->depth[comp])-1);
	//	float gain=(float)(wndw>>1)/((1<<im1->depth[comp])-1);
		int offset=1<<im1->depth[comp]>>1;
		int stride=im1->iw<<2;
		for(int ky=0;ky<wndh;++ky)
		{
			iy=screen2image_y_int(ky);
			if((unsigned)iy<(unsigned)im1->ih)
				x2=(col[iy*stride]+offset)*gain;
			else
				x2=0;
			draw_curve_enqueue(&vertices_2d, x2, (float)ky);
		}
		draw_2d_flush(vertices_2d, color, GL_LINE_STRIP);
	}
}

#if 0
#define DSP_REACH 5
//#define DSP_NPARAMS (2*(DSP_REACH+1)*DSP_REACH)
static double dspparams[]=
{
	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0.1,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0.1,	0,	0.2,	0,	0,	0,	0,
	0,	0,	0,	0.1,	0,	0,	0,	0.1,	0.1,	0,	0,
	0,	0,	0.1,	0.1,	0,	0.1,	(0),	(0),	(0),	(0),	(0),

	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0.2,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0.2,	0,	0.1,	0.2,	0,	0,	0,
	//0,	0,	0,	0.2,	0,	0.1,	(0),	(0),	(0),	(0),	(0),
	
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0.05,	0.05,	0.05,	0.05,	0.25,	0,	0,
	//0,	0,	0,	0.3,	0.20,	0.05,	(0),	(0),	(0),	(0),	(0),

	//curr = (2*W+curr+NEEE)>>2:
	//0,	0,	0,	0.3,	0.25,	(0.25),	(0),	(0),	(0),	(0),	(0),
	//0,	0,	0,	0,	0.05,	0.05,	0.05,	0.05,	0.25,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,

	//0,	0,	0,	0,	0.15,	(0.25),	(0),	(0),	(0),	(0),	(0),
	//0,	0,	0,	0,	0.15,	0.15,	0.15,	0.15,	0.25,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,

	//curr = (2*W+curr+NEEE)>>2:
	//0,	0,	0,	0,	0.5,	(0.25),	(0),	(0),	(0),	(0),	(0),
	//0,	0,	0,	0,	0,	0,	0,	0,	0.25,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	//0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
};
static double *dspbuf=0;
static int *exbuf=0;
static const int dspw=1024, dsph=512;
static void dsp_update()
{
	if(!dspbuf)
	{
		dspbuf=(double*)malloc(sizeof(double)*dspw*dsph);
		exbuf=(int*)malloc(sizeof(int)*dspw*dsph);
		if(!dspbuf||!exbuf)
		{
			LOG_ERROR("Alloc error");
			return;
		}
	}
	memset(dspbuf, 0, sizeof(double)*dspw*dsph);
	//dspbuf[dspw>>1]=256;
	for(int ky=0;ky<dsph;++ky)
	{
		for(int kx=0;kx<dspw;++kx)
		{
			int dx=kx-dspw/2, dy=ky-dsph/2, dist=100-abs(dx*dx+dy*dy-100*100);
			if(dist<0)
				dist=0;
			double val=255*dist;
			if((unsigned)(kx-dspw/2+25)<50&&(unsigned)(ky-dsph/2+25)<50)
				val=255;
			if((unsigned)(kx-dspw/2+10)<20&&(unsigned)(ky-dsph/2+10)<20)
				val=0;
			dspbuf[dspw*ky+kx]=val;
			//double val=(unsigned)(kx-dspw/2)<10&&(unsigned)(ky-dsph/4)<10?255:0;
			//double val=kx==dspw/2&&ky==1?255:0;
			//double val=dspbuf[dspw*ky+kx];

			val=0;
			for(int ky2=-DSP_REACH;ky2<=0;++ky2)
			{
				for(int kx2=-DSP_REACH;kx2<=DSP_REACH;++kx2)
				{
					if((unsigned)(ky+ky2)<dsph&&(unsigned)(kx+kx2)<dspw)
						val+=dspbuf[dspw*(ky+ky2)+kx+kx2]*dspparams[(DSP_REACH<<1|1)*(ky2+DSP_REACH)+kx2+DSP_REACH];
					if(!ky2&&!kx2)
						break;
				}
			}
			dspbuf[dspw*ky+kx]=val;
			int ival=(int)val;
			CLAMP2_32(ival, ival, 0, 255);
			exbuf[dspw*ky+kx]=0xFF000000|ival<<16|ival<<8|ival;
		}
	}
	for(int kx=0;kx<dspw;++kx)
		exbuf[dspw*(dsph/2)+kx]^=255;
	for(int ky=0;ky<dsph;++ky)
		exbuf[dspw*ky+dspw/2]^=255;
	//double *ptr=dspparams+DSP_REACH-1;
	//*ptr+=0.001;
	//if(*ptr>0.75)
	//	*ptr=0;
}
static void dsp_render()
{
	int
		sx1=image2screen_x_int(0), sx2=image2screen_x_int(dspw),
		sy1=image2screen_y_int(0), sy2=image2screen_y_int(dsph);
	if(exbuf)
		display_texture_i(sx1, sx2, sy1, sy2, exbuf, dspw, dsph, 0, 1, 0, 1, 1, 0);
	if(drag)
	{
		for(int ky=0;ky<DSP_REACH;++ky)
		{
			for(int kx=0;kx<(DSP_REACH<<1|1);++kx)
			{
				GUIPrint(0, (float)(kx*wndw/(DSP_REACH<<1|1)), (float)(ky*wndh/DSP_REACH), guizoom, "%lf", dspparams[(DSP_REACH<<1|1)*ky+kx]);
			}
		}
	}
}
#endif

static ArrayHandle vertices_text=0;
static void print_pixellabels(int ix1, int ix2, int iy1, int iy2, int component, char label, long long txtcolors, int depth)
{
	const char *format;
	int nlevels=1<<depth, half=nlevels>>1, mask=nlevels-1;
	if(pxlabels_hex)
	{
		if(depth<=4)
			format="%c %01X";
		else if(depth<=8)
			format="%c %02X";
		else if(depth<=12)
			format="%c %03X";
		else
			format="%c%04X";
	}
	else
	{
		format="%c%5d";
		half=0;
		mask=~0;
	}
	txtcolors=set_text_colors(txtcolors);
	int iy=MAXVAR(iy1, 0);
	float fontsize=1, labeloffset=tdy*fontsize*component;
	for(int yend=MINVAR(iy2+2, im1->ih);iy<yend;++iy)
	{
		int ky=image2screen_y_int(iy);
		int ix=MAXVAR(ix1, 0);
		for(int xend=MINVAR(ix2+2, im1->iw);ix<xend;++ix)
		{
			int kx=image2screen_x_int(ix);
			int idx=(im1->iw*iy+ix)<<2;
			GUIPrint_enqueue(&vertices_text, 0, (float)kx, (float)ky+labeloffset, fontsize, format, label, (im1->data[idx+component]+half)&mask);
		}
	}
	print_line_flush(vertices_text, fontsize);
	txtcolors=set_text_colors(txtcolors);
}
void io_render(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	if(!wndh)
		return;

	//{
	//	dsp_update();//
	//	dsp_render();//
	//}
	{
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
	}
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
		int
			sx1=image2screen_x_int(0), sx2=image2screen_x_int(im1->iw),
			sy1=image2screen_y_int(0), sy2=image2screen_y_int(im1->ih);
		switch(mode)
		{
		//case VIS_PLANES:		chart_planes_draw();	break;
		//case VIS_MESH:		chart_mesh_draw();	break;
		//case VIS_MESH_SEPARATE:	chart_mesh_sep_draw();	break;
		//	{
		//		float yoffset=tdy*3;
		//		display_texture_i(0, im1->iw, (int)yoffset, (int)yoffset+im1->ih, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
		//		chart_hist_draw(0, (float)wndw, 0, (float)wndh, 0, 3, 0, 0x60, hist, histmax);
		//	}
		//	break;
		case VIS_JOINT_HISTOGRAM:
			chart_jointhist_draw();
			break;
		case VIS_ANALYSIS:
			{
				const unsigned char *group=rct2_combinations[analysis_info.bestrct];
				int selected_ch[]=
				{
					group[0]*PRED_COUNT+analysis_info.predidx[0],
					group[1]*PRED_COUNT+analysis_info.predidx[1],
					group[2]*PRED_COUNT+analysis_info.predidx[2],
				};
				int bounds[4]={0};//{x1, x2, y1, y2}
				float xmid, ymid;
				jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, im1->iw, im1->ih, bounds);
				
				display_texture_i(0, wndw, 0, wndh, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);

				{
					int crosshaircolor=0xFF000000;
					float ratios[]={(float)wndw/im1->iw, (float)wndh/im1->ih};
					float sbounds[4];
					for(int k=0;k<4;++k)
						sbounds[k]=bounds[k]*ratios[k>>1];
					xmid=(sbounds[0]+sbounds[1])*0.5f;
					ymid=(sbounds[2]+sbounds[3])*0.5f;
					draw_rect_hollow(sbounds[0]+1, sbounds[1]+1, sbounds[2]+1, sbounds[3]+1, 0xFFFFFFFF);
					draw_rect_hollow(sbounds[0], sbounds[1], sbounds[2], sbounds[3], crosshaircolor);
					//draw_line(sbounds[0], sbounds[2], sbounds[0], sbounds[3], crosshaircolor);//{x1, y1, x2, y2}
					//draw_line(sbounds[1], sbounds[2], sbounds[1], sbounds[3], crosshaircolor);
					//draw_line(sbounds[0], sbounds[2], sbounds[1], sbounds[2], crosshaircolor);
					//draw_line(sbounds[0], sbounds[3], sbounds[1], sbounds[3], crosshaircolor);
					draw_line(xmid, 0, xmid, (float)wndh, crosshaircolor);
					draw_line(0, ymid, (float)wndw, ymid, crosshaircolor);
				}
				float y=tdy*2;
				GUIPrint(0, (float)(0.25*wndw), y-tdy*1.5f-2, 1.5f,
					"3*%d*%d: %6.2lf%% + %6.2lf%% + %6.2lf%% = %6.2lf%%  %s %s %s %s  %lf MB/s",
					bounds[1]-bounds[0], bounds[3]-bounds[2],
					analysis_info.esizes[selected_ch[0]]*100,
					analysis_info.esizes[selected_ch[1]]*100,
					analysis_info.esizes[selected_ch[2]]*100,
					(analysis_info.esizes[selected_ch[0]]+analysis_info.esizes[selected_ch[1]]+analysis_info.esizes[selected_ch[2]])*100/3,
					rct2_names[analysis_info.bestrct],
					pred_names[analysis_info.predidx[0]],
					pred_names[analysis_info.predidx[1]],
					pred_names[analysis_info.predidx[2]],
					3*(bounds[1]-bounds[0])*(bounds[3]-bounds[2])/analysis_info.t_analysis/(1024*1024)
				);
				{
					float
						rx1=(float)(0.25*wndw),
						rx2=(float)(0.75*wndw),
						ry1=y,
						ry2=y+RCT2_COUNT*3+OCH2_COUNT*PRED_COUNT*3+OCH2_COUNT*8;
					display_texture_i(
						(int)rx1, (int)rx2, (int)ry1, (int)ry2,
						(int*)im_export,
						im1->iw, im1->ih,
						(float)bounds[0]/im1->iw,
						(float)bounds[1]/im1->iw,
						(float)bounds[2]/im1->ih,
						(float)bounds[3]/im1->ih,
						0.5f, 0
					);
					draw_rect_hollow(rx1-1, rx2+1, ry1+0, ry2+0, 0xC0000000);
					draw_rect_hollow(rx1-2, rx2+0, ry1-1, ry2-1, 0xC0FFFFFF);
				}
				//draw_line((float)(0.25*wndw)+0, y, (float)(0.25*wndw)+0, y+OCH2_COUNT*PRED_COUNT*4+OCH2_COUNT*10, 0xC0000000);
				//draw_line((float)(0.25*wndw)-1, y, (float)(0.25*wndw)-1, y+OCH2_COUNT*PRED_COUNT*4+OCH2_COUNT*10, 0xC0FFFFFF);
				//draw_line((float)(0.75*wndw)+1, y, (float)(0.75*wndw)+1, y+OCH2_COUNT*PRED_COUNT*4+OCH2_COUNT*10, 0xC0000000);
				//draw_line((float)(0.75*wndw)+0, y, (float)(0.75*wndw)+0, y+OCH2_COUNT*PRED_COUNT*4+OCH2_COUNT*10, 0xC0FFFFFF);
				y+=10;
				float centers[6]={0};
				int ncenters=0;
				for(int k=0;k<RCT2_COUNT;++k)//RCTs
				{
					int hit=k==analysis_info.bestrct;
					int color=hit?0xC04040FF:0xC0FFFFFF;
					float x1=(float)(0.25*wndw), x2=(float)(wndw*(0.25+0.5*analysis_info.rctsizes[k]));
					draw_line(x1, y+0, x2, y+0, color);
					draw_line(x1, y+1, x2, y+1, color);
					draw_line(x1, y+2, x2, y+2, 0xC0000000);
					y+=3;
				}
				y+=9;
				for(int k=0;k<OCH2_COUNT*PRED_COUNT;++k)//output channels
				{
					int hit=k==selected_ch[0]||k==selected_ch[1]||k==selected_ch[2];
					int color=hit?0xC04040FF:0xC0FFFFFF;
					float x1=(float)(0.25*wndw), x2=(float)(wndw*(0.25+0.5*analysis_info.esizes[k]));
					draw_line(x1, y+0, x2, y+0, color);
					draw_line(x1, y+1, x2, y+1, color);
					draw_line(x1, y+2, x2, y+2, 0xC0000000);
					if(hit)
					{
						centers[ncenters++]=x2;
						centers[ncenters++]=y+2;
					}
					if(!((k+1)%(PRED_COUNT)))
						GUIPrint(0, (float)(0.25*wndw)-40, y-30, 2, "%s", och_names[k/PRED_COUNT]);
					y+=(k+1)%PRED_COUNT?3:9;
					//if(!((k+1)%(PRED_COUNT*3)))
					//	y+=10;
				}
				//long long hmax=0;
				//for(int k=0;k<PRED_COUNT;++k)
				//{
				//	if(hmax<pred_hist[k])
				//		hmax=pred_hist[k];
				//}
				for(int k=0;k<PRED_COUNT;++k)
				{
					float
						x1=(float)(0.25*wndw),
						yk=y+tdy*1.2f*(k+1);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+1-3, (float)(x1-pred_hist[0][k]), (float)(yk+tdy*(1.2/2))+1-3, 0xC00000FF);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+0-3, (float)(x1-pred_hist[0][k]), (float)(yk+tdy*(1.2/2))+0-3, 0xC00000FF);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+1+0, (float)(x1-pred_hist[1][k]), (float)(yk+tdy*(1.2/2))+1+0, 0xC000FF00);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+0+0, (float)(x1-pred_hist[1][k]), (float)(yk+tdy*(1.2/2))+0+0, 0xC000FF00);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+1+3, (float)(x1-pred_hist[2][k]), (float)(yk+tdy*(1.2/2))+1+3, 0xC0FF0000);
					draw_line(x1, (float)(yk+tdy*(1.2/2))+0+3, (float)(x1-pred_hist[2][k]), (float)(yk+tdy*(1.2/2))+0+3, 0xC0FF0000);
					GUIPrint(0, x1+2, yk, 1.2f, "%4s  %s", pred_names[k], pred_desc[k]);
				}
				draw_ellipse(centers[0]-10, centers[0]+10, centers[1]-10, centers[1]+10, 0x8000FF00);
				draw_ellipse(centers[2]-10, centers[2]+10, centers[3]-10, centers[3]+10, 0x8000FF00);
				draw_ellipse(centers[4]-10, centers[4]+10, centers[5]-10, centers[5]+10, 0x8000FF00);
				GUIPrint(0, xmid, ymid, 1,
					"%s %s %s %s",
					rct2_names[analysis_info.bestrct],
					pred_names[analysis_info.predidx[0]],
					pred_names[analysis_info.predidx[1]],
					pred_names[analysis_info.predidx[2]]
				);
			}
			break;
#ifdef ENABLE_L1WEIGHTS
		case VIS_L1WEIGHTS:
#endif
		case VIS_IMAGE:
		case VIS_HISTOGRAM:
		case VIS_MODEL:
		case VIS_ZIPF:
			{
				display_texture_i(sx1, sx2, sy1, sy2, (int*)(mode==VIS_ZIPF?zimage:im_export), im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
#if 0
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
#endif
				if(imzoom>=ZOOM_LIMIT_LABEL)
				{
					int
						csx1=CLAMP(0, sx1, wndw),
						csx2=CLAMP(0, sx2, wndw),
						csy1=CLAMP(0, sy1, wndh),
						csy2=CLAMP(0, sy2, wndh);
					int
						ix1=screen2image_x_int(csx1),
						ix2=screen2image_x_int(csx2),
						iy1=screen2image_y_int(csy1),
						iy2=screen2image_y_int(csy2);
					long long theme[]=
					{
						0xC00000FF80000000,
						0xC000FF0080000000,
						0xC0FF000080FFFFFF,
						0xC0FFFFFF80000000,
					};
					if(viewmode==VIEW_RGB||viewmode==VIEW_C0)
						print_pixellabels(ix1, ix2, iy1, iy2, 0, 'r', theme[0], im1->depth[0]);
					if(viewmode==VIEW_RGB||viewmode==VIEW_C1)
						print_pixellabels(ix1, ix2, iy1, iy2, 1, 'g', theme[1], im1->depth[1]);
					if(viewmode==VIEW_RGB||viewmode==VIEW_C2)
						print_pixellabels(ix1, ix2, iy1, iy2, 2, 'b', theme[2], im1->depth[2]);
					if(imzoom>=ZOOM_LIMIT_ALPHA&&im1->depth[3])
						print_pixellabels(ix1, ix2, iy1, iy2, 3, 'a', theme[3], im1->depth[3]);
				}
#ifdef ENABLE_L1WEIGHTS
#if 0
				{
					static int chasers[5][2]={0};
					for(int k=0;k<5;++k)
					{
						int updatex=(mx-chasers[k][0]+(1<<(k+1)>>1))>>(k+1);
						int updatey=(my-chasers[k][1]+(1<<(k+1)>>1))>>(k+1);
						updatex+=mx-chasers[k][0];
						updatey+=my-chasers[k][1];
						int speed2=updatex*updatex+updatey*updatey;
						int speedcap=64-k*2;
						if(speed2>speedcap*speedcap)
						{
							float invspeed=speedcap/sqrtf((float)speed2);
							updatex=(int)round(updatex*invspeed);
							updatey=(int)round(updatey*invspeed);
						}
						chasers[k][0]+=updatex;
						chasers[k][1]+=updatey;
						int color=255<<24|rand()<<15|rand();
						GUIPrint(0, (float)(chasers[k][0]+2), (float)(chasers[k][1]+2), 1, "%d", k+1);
						draw_line(
							(float)(chasers[k][0]-10), (float)(chasers[k][1]),
							(float)(chasers[k][0]+10), (float)(chasers[k][1]),
							color
						);
						draw_line(
							(float)(chasers[k][0]), (float)(chasers[k][1]-10),
							(float)(chasers[k][0]), (float)(chasers[k][1]+10),
							color
						);
					}
				}
#endif
				if(mode==VIS_L1WEIGHTS)
				{
	#define L1HOLD 1
	#define L1SCALE 1
	#define L1HISTSIZE 1024
					static int vidx=0;
					static int xpos=0;
					static int whist[L1HISTSIZE][(L1NPREDS+1)+5]={0};//+{pixel, pred, delta, CG, CGdelta}
					ptrdiff_t res=(ptrdiff_t)im1->iw*im1->ih;
					int kx=0, ky=0;
					for(int kpx=0;kpx<L1SPEED;++kpx)
					{
						++vidx;
						if(vidx>=res)
							vidx=0;
						kx=vidx%im1->iw;
						ky=vidx/im1->iw;
						int
							NNNN	=(unsigned)(ky-4)<(unsigned)im1->ih&&(unsigned)(kx+0)<(unsigned)im1->iw?im1->data[((ky-4)*im1->iw+kx+0)*4]:0,
							NNNW	=(unsigned)(ky-3)<(unsigned)im1->ih&&(unsigned)(kx-1)<(unsigned)im1->iw?im1->data[((ky-3)*im1->iw+kx-1)*4]:0,
							NNN	=(unsigned)(ky-3)<(unsigned)im1->ih&&(unsigned)(kx+0)<(unsigned)im1->iw?im1->data[((ky-3)*im1->iw+kx+0)*4]:0,
							NNNE	=(unsigned)(ky-3)<(unsigned)im1->ih&&(unsigned)(kx+1)<(unsigned)im1->iw?im1->data[((ky-3)*im1->iw+kx+1)*4]:0,
							NNW	=(unsigned)(ky-2)<(unsigned)im1->ih&&(unsigned)(kx-1)<(unsigned)im1->iw?im1->data[((ky-2)*im1->iw+kx-1)*4]:0,
							NN	=(unsigned)(ky-2)<(unsigned)im1->ih&&(unsigned)(kx+0)<(unsigned)im1->iw?im1->data[((ky-2)*im1->iw+kx+0)*4]:0,
							NNE	=(unsigned)(ky-2)<(unsigned)im1->ih&&(unsigned)(kx+1)<(unsigned)im1->iw?im1->data[((ky-2)*im1->iw+kx+1)*4]:0,
							NNEEE	=(unsigned)(ky-2)<(unsigned)im1->ih&&(unsigned)(kx+3)<(unsigned)im1->iw?im1->data[((ky-2)*im1->iw+kx+3)*4]:0,
							NWW	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx-2)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx-2)*4]:0,
							NW	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx-1)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx-1)*4]:0,
							N	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+0)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+0)*4]:0,
							NE	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+1)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+1)*4]:0,
							NEE	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+2)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+2)*4]:0,
							NEEE	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+3)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+3)*4]:0,
							NEEEE	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+4)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+4)*4]:0,
							NEEEEE	=(unsigned)(ky-1)<(unsigned)im1->ih&&(unsigned)(kx+5)<(unsigned)im1->iw?im1->data[((ky-1)*im1->iw+kx+5)*4]:0,
							WWWWWW	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-6)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-6)*4]:0,
							WWWWW	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-5)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-5)*4]:0,
							WWWW	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-4)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-4)*4]:0,
							WWW	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-3)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-3)*4]:0,
							WW	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-2)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-2)*4]:0,
							W	=(unsigned)(ky-0)<(unsigned)im1->ih&&(unsigned)(kx-1)<(unsigned)im1->iw?im1->data[((ky-0)*im1->iw+kx-1)*4]:0,
							curr	=im1->data[(ky*im1->iw+kx)*4];
						int preds[]=
						{
	#define L1PRED(W0, EXPR) EXPR,
							L1PREDLIST
	#undef  L1PRED
						};
						int vmax=N, vmin=W;
						if(N<W)vmin=N, vmax=W;
						CLAMP2(preds[0], vmin, vmax);
						int *currw=l1weights;
						int pred=0;
						if(use_ols==2)
						{
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
							static int divlookup[DIVLUTSIZE]={0};
							if(!*divlookup)
							{
								for(int k=0;k<DIVLUTSIZE;++k)
									divlookup[k]=(1<<NUMBITS)/(k+1);
							}
							unsigned wsum=0;
							int ipred=0;
							int sh=0;
							for(int kp=0;kp<L1NPREDS;++kp)
							{
								int e=wperrors[kp];
								int sh=FLOOR_LOG2(e+1);
								currw[kp]=(divlookup[e<<(DENBITS-1)>>sh]<<(DENBITS-1)>>sh)+(1<<DENBITS>>2)/L1NPREDS;
								wsum+=currw[kp];
							}
							sh=FLOOR_LOG2(wsum)-(DENBITS-2);
							wsum=0;
							for(int kp=0;kp<L1NPREDS;++kp)
							{
								int c=(int)(currw[kp]>>sh);
								wsum+=c;
								currw[kp]=c;
							}
							ipred=wsum>>1;
							for(int kp=0;kp<L1NPREDS;++kp)
							{
								ipred+=(int)(currw[kp]*preds[kp]);
								currw[kp]<<=8;
							}
							//if(wsum<0)
							//	goto again;
							pred=ipred*divlookup[wsum-1]>>NUMBITS;
						}
						else if(use_ols==1)
						{
							double fpred=0;
							for(int k=0;k<L1NPREDS;++k)
							{
								fpred+=curr_params[k]*preds[k];
								currw[k]=(int)round(curr_params[k]*0x20000);
							}
							pred=(int)CVTFP64_I64(fpred);
						}
						else
						{
							pred=currw[L1NPREDS];
							for(int k=0;k<L1NPREDS;++k)
								pred+=currw[k]*preds[k];
		#define L1SH 19
	//	#define L1SH 20	//DIV2K
	//	#define L1SH 20	//GDCC
	//	#define L1SH 15	//synth
							pred+=1<<L1SH>>1;
							pred>>=L1SH;
						}
						if(vmin>NE)vmin=NE;
						if(vmax<NE)vmax=NE;
						//if(vmin>NW)vmin=NW;
						//if(vmax<NW)vmax=NW;
						if(vmin>NEEE)vmin=NEEE;
						if(vmax<NEEE)vmax=NEEE;
						CLAMP2(pred, vmin, vmax);

						++xpos;
						if(xpos>=L1HISTSIZE-1||xpos>=wndw-1)
							xpos=0;
						for(int k=0;k<L1NPREDS+1;++k)
							whist[xpos][k]=l1weights[(L1NPREDS+1)-1-k];
						whist[xpos][(L1NPREDS+1)+0]=preds[0]		+(1<<im1->depth[0]>>1)+0*64;
						whist[xpos][(L1NPREDS+1)+1]=curr		+(1<<im1->depth[0]>>1)+1*64;
						whist[xpos][(L1NPREDS+1)+2]=pred		+(1<<im1->depth[0]>>1)+2*64;
						whist[xpos][(L1NPREDS+1)+3]=curr-pred		+(1<<im1->depth[0]>>1);
						whist[xpos][(L1NPREDS+1)+4]=curr-preds[0]	+(1<<im1->depth[0]>>1)+(1<<im1->depth[0]);
						{
							static int ctr=0, use_ols0=0;
							static int hist[256]={0};
							static double e=0;
							++hist[(curr-pred)&255];
							++ctr;
							//if((ctr&0xFFF)==0xFFF||use_ols0!=use_ols)
							//{
								e=0;
								for(int ks=0;ks<256;++ks)
								{
									int freq=hist[ks];
									if(freq)
									{
										double p=(double)freq/ctr;
										e-=p*log2(p);
									}
								}
								e/=8;
								if(use_ols0!=use_ols)
								{
									ctr=0;
									memset(hist, 0, sizeof(hist));
									use_ols0=use_ols;
								}
							//}
							static const char *labels[]=
							{
								"L1",
								"OLS",
								"WP",
							};
							GUIPrint(0, 0, (float)wndh*0.5f, 2, "E %10.6lf%% <- %8d %s", e*100, ctr, labels[use_ols]);
						}
						//for(int k=0, wsum=0;k<L1NPREDS+1;++k)
						//{
						//	int w=l1weights[L1NPREDS-k];
						//	wsum+=w;
						//	whist[xpos][k]=wsum>>11;
						//}
						//memcpy(whist[xpos], l1weights, sizeof(l1weights));

						//update
						if(use_ols==2)
						{
							int e2[L1NPREDS], best=0x7FFFFFFF;
							for(int kp=0;kp<L1NPREDS;++kp)
							{
								int e=abs(curr-preds[kp]);
								if(best>e)
									best=e;
								e2[kp]=e;
							}
							for(int kp=0;kp<L1NPREDS;++kp)
							{
								int e=e2[kp]-best;
								wperrors[kp]+=((e<<6)-wperrors[kp]+(1<<3>>1))>>3;
							}
						}
						else if(use_ols==1)
						{
							double lr=0.00001;
							for(int ky2=0, midx=0;ky2<L1NPREDS;++ky2)
							{
								for(int kx2=0;kx2<L1NPREDS;++kx2, ++midx)
									curr_cov[midx]+=(preds[kx2]*preds[ky2]-curr_cov[midx])*lr;
							}
							double lval=curr*lr, lr_comp=1-lr;
							for(int k=0;k<L1NPREDS;++k)
								curr_vec[k]=lval*preds[k]+lr_comp*curr_vec[k];
							int success=1;
							memcpy(curr_cholesky, curr_cov, sizeof(double[L1NPREDS*L1NPREDS]));
							for(int k=0;k<L1NPREDS*L1NPREDS;k+=L1NPREDS+1)
								curr_cholesky[k]+=0.0075;
							double sum;
							for(int i=0;i<L1NPREDS;++i)
							{
								for(int j=0;j<i;++j)
								{
									sum=curr_cholesky[i*L1NPREDS+j];
									for(int k=0;k<j;++k)
										sum-=curr_cholesky[i*L1NPREDS+k]*curr_cholesky[j*L1NPREDS+k];
									curr_cholesky[i*L1NPREDS+j]=sum/curr_cholesky[j*L1NPREDS+j];
								}
								sum=curr_cholesky[i*L1NPREDS+i];
								for(int k=0;k<i;++k)
									sum-=curr_cholesky[i*L1NPREDS+k]*curr_cholesky[i*L1NPREDS+k];
								if(sum<=1e-8)
								{
									success=0;
									break;
								}
								curr_cholesky[i*L1NPREDS+i]=sqrt(sum);
							}
							if(success)
							{
								for(int i=0;i<L1NPREDS;++i)
								{
									sum=curr_vec[i];
									for(int j=0;j<i;++j)
										sum-=curr_cholesky[i*L1NPREDS+j]*curr_params[j];
									curr_params[i]=sum/curr_cholesky[i*L1NPREDS+i];
								}
								for(int i=L1NPREDS-1;i>=0;--i)
								{
									sum=curr_params[i];
									for(int j=i+1;j<L1NPREDS;++j)
										sum-=curr_cholesky[j*L1NPREDS+i]*curr_params[j];
									curr_params[i]=sum/curr_cholesky[i*L1NPREDS+i];
								}
							}
						}
						else
						{
#if 0
							int e=curr-pred;
							CLAMP2(e, 1-8, 15-8);
							static const int etable[]=
							{
								//   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
								-1, -1, -1, -1, -2, -2, -2, -2,  0, +2, +2, +2, +2, +1, +1, +1,
							//	-4, -4, -4, -5, -6, -8, -9, -6,  0, +6, +9, +8, +6, +5, +4, +4,
							};
							e=etable[e+8];
							//	... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7
							//	... -1 -1 -1 -1 -2 -2 -2 -2  0 +2 +2 +2 +2 +1 +1 +1 ...
							//e=(((e>0)-(e<0))<<1)-((e>4)-(e<-4));
#endif
							int e=(curr>pred)-(curr<pred);
							currw[L1NPREDS]+=e;
							for(int k=0;k<L1NPREDS;++k)
								currw[k]+=e*preds[k];
						}
					}

					//draw
					int dlen=wndw;
					if(dlen>L1HISTSIZE)
						dlen=L1HISTSIZE;
					dlen-=2;
					int vmin2=0, vmax2=0;
					for(int kx2=0;kx2<dlen;++kx2)
					{
						for(int ky2=0;ky2<L1NPREDS+1;++ky2)
						{
							int val=whist[kx2][ky2];
							if(vmin2>val)vmin2=val;
							if(vmax2<val)vmax2=val;
						}
					}
					static int vmin=0, vmax=0;
					vmin+=(vmin2-vmin+(1<<3>>1))>>3;
					vmax+=(vmax2-vmax+(1<<3>>1))>>3;
					if(vmin>=vmax)
					{
						vmin=vmin2;
						vmax=vmax2;
					}
					if(wndh&&vmax>vmin)
					{
						//map {vmin, vmax} -> {wndh*7/8, wndh/8}  ys = C1*yu+C0
						float
							yC1=(float)wndh*6/8/(vmin-vmax),
							yC0=(float)wndh*7/8-yC1*vmin;//(y-wndh*7/8)/(0-vmin) = C1  ->  y = wndh*7/8-C1*vmin
						for(int kx2=0;kx2<dlen;++kx2)
						{
							if(kx2==xpos)
								continue;
							for(int ky2=0;ky2<L1NPREDS+1;++ky2)
								draw_line(
									L1SCALE*(float)kx2, whist[kx2][ky2]*yC1+yC0,
									L1SCALE*(float)(kx2+1), whist[kx2+1][ky2]*yC1+yC0, 0xFF000000
								);
							draw_line(//CG - red
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+0]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+0]), 0xFF0000FF
							);
							draw_line(//curr - pink
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+1]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+1]), 0xFFFF00FF
							);
							draw_line(//pred - green
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+2]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+2]), 0xFF00FF00
							);
							draw_line(//delta - gray
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+3]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+3]), 0xFFC0C0C0
							);
							draw_line(//CGdelta - grayish-red
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+4]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+4]), 0xFF5050B0
							);
						}
						static float labely[L1NPREDS+1]={0};
						for(int ky2=0;ky2<L1NPREDS+1;++ky2)
						{
							float y=whist[xpos][ky2]*yC1+yC0;
							labely[ky2]+=(y-labely[ky2])*(1.f/32);
							draw_line(
								(float)xpos, (float)whist[xpos][ky2]*yC1+yC0,
								L1SCALE*L1HISTSIZE, labely[ky2], 0x80000000
							);
							GUIPrint(0, L1SCALE*L1HISTSIZE, labely[ky2], 0.8f, "%10d %s",
								l1weights[(L1NPREDS+1)-1-ky2],
								ky2?l1prednames[L1NPREDS-1-(ky2-1)]:"bias"
							);
						}
					}
#if 0
					for(int kx2=0;kx2<dlen;++kx2)
					{
						if(kx2==xpos)
							continue;
						for(int ky2=0;ky2<L1NPREDS+1;++ky2)
							draw_line(
								L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][ky2]),
								L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1][ky2]), 0xFF000000
							);
						draw_line(//CG - red
							L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+0]),
							L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+0]), 0xFF0000FF
						);
						draw_line(//curr - pink
							L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+1]),
							L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+1]), 0xFFFF00FF
						);
						draw_line(//pred - green
							L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+2]),
							L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+2]), 0xFF00FF00
						);
						draw_line(//delta - gray
							L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+3]),
							L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+3]), 0xFFC0C0C0
						);
						draw_line(//CGdelta - grayish-red
							L1SCALE*(float)kx2, (float)(wndh*7/8-whist[kx2][(L1NPREDS+1)+4]),
							L1SCALE*(float)(kx2+1), (float)(wndh*7/8-whist[kx2+1*L1HOLD][(L1NPREDS+1)+4]), 0xFF5050B0
						);
					}
					for(int ky2=0;ky2<L1NPREDS+1;++ky2)
						GUIPrint(0, L1SCALE*L1HISTSIZE, (float)(wndh*7/8-whist[xpos][ky2]), 0.8f, "%10d %s",
							l1weights[(L1NPREDS+1)-1-ky2],
							ky2?l1prednames[L1NPREDS-1-(ky2-1)]:"bias"
						);
#endif
					{
						float xs=(float)image2screen_x(kx);
						float ys=(float)image2screen_y(ky);
						draw_rect_hollow(xs-7, xs+7, ys-7, ys+7, 0xFFFF00FF);
					}
				}
				else
#endif
				if(mode==VIS_HISTOGRAM)
					chart_hist_draw(0, (float)wndw, 0, (float)wndh, 0, 3, 0, 0x60, hist, histmax);
				else if(mode==VIS_MODEL)
				{
					static const int modeltheme[]=
					{
						0x600000FF,
						0x6000FF00,
						0x60FF0000,
						0x60C0C0C0,
						0x80C0C0C0,
						0x80FF0000,
						0x800000FF,
						0x80C0C0C0,
					};
					int nhist=modelnch*modelnctx;
					int YUV=0;
					for(int k=0;k<CST_COMPARE;++k)
					{
						if(transforms_mask[k])
						{
							YUV=1;
							break;
						}
					}
					int wndh2=wndh-(int)(tdy*5);
					int xhist=modelnch, yhist=nhist/modelnch;
					//if(wndw<wndh)
					//{
					//	xhist=nhist*wndw/wndh;
					//	CLAMP2(xhist, 1, nhist);
					//	yhist=(nhist+xhist-1)/xhist;
					//}
					//else
					//{
					//	yhist=nhist*wndh/wndw;
					//	CLAMP2(yhist, 1, nhist);
					//	xhist=(nhist+yhist-1)/yhist;
					//}
					double msizes=0;
					double esizes[4]={0};
					for(int kh=0;kh<nhist;++kh)
					{
						int hx=kh/yhist, hy=kh%yhist;
						//int hy=kh/xhist, hx=kh%xhist;
						int xstart=hx*wndw/xhist, xend=(hx+1)*wndw/xhist;
						int ystart=hy*wndh2/yhist, yend=(hy+1)*wndh2/yhist;
						int kc=kh/modelnctx;
						int *curr_hist=modelhist+((ptrdiff_t)kh<<modeldepth);
						int nlevels=1<<modeldepth;
						int fmax=0;
						for(int ks=0;ks<nlevels;++ks)//FIXME
						{
							int freq=curr_hist[ks];
							if(fmax<freq)
								fmax=freq;
						}
						if(fmax)
						{
							for(int ks=0;ks<nlevels;++ks)
							{
								int freq=curr_hist[ks];
								int x1=xstart+(xend-xstart)*ks/nlevels;
								int x2=xstart+(xend-xstart)*(ks+1)/nlevels;
								int y1=yend-(yend-ystart)*freq/fmax;
								int y2=yend;
								draw_rect((float)x1, (float)x2, (float)y1, (float)y2, modeltheme[YUV<<2|kc]);
							}
						}
						int *curr_hist2=modelmhist+((ptrdiff_t)kh<<modeldepth);
						int fmax2=0;
						for(int ks=0;ks<nlevels;++ks)//FIXME
						{
							int freq=curr_hist2[ks];
							if(fmax2<freq)
								fmax2=freq;
						}
						if(fmax2)
						{
							for(int ks=0;ks<nlevels-1;++ks)
							{
								int freq1=curr_hist2[ks+0];
								int freq2=curr_hist2[ks+1];
								int x1=xstart+(xend-xstart)*(2*ks+1)/(2*nlevels);
								int x2=xstart+(xend-xstart)*(2*ks+3)/(2*nlevels);
								int y1=yend-(yend-ystart)*freq1/fmax2;
								int y2=yend-(yend-ystart)*freq2/fmax2;
								draw_line((float)x1, (float)y1, (float)x2, (float)y2, modeltheme[YUV<<2|kc]);
							}
						}
						GUIPrint((float)xstart, (float)xstart, (float)ystart, 1, "%02d %12.2lf %+12.2lf", kh%modelnctx, modelcsizes[kh], modelmsizes[kh]-modelcsizes[kh]);
						GUIPrint((float)xstart, (float)xstart, (float)ystart+tdy, 1, "   %12.2lf", modelmsizes[kh]);
						//GUIPrint((float)xstart, (float)xstart, (float)ystart, 1, "%02d%c", kh%modelnctx, " X"[!fmax]);
						esizes[kc]+=modelcsizes[kh];
						msizes+=modelmsizes[kh];
					}
					double esize=0;
					for(int kc=0;kc<modelnch;++kc)
					{
						float x=wndw*(2*(float)kc+1)/(2*(float)modelnch)-tdx*2*12;
						GUIPrint(x, x, (float)(wndh>>1), 2, "%12.2lf", esizes[kc]);
						esize+=esizes[kc];
					}
					float x=(float)(wndw>>1)-tdx*2*12;
					GUIPrint(x, x, (float)(wndh>>1)+tdy*2, 2, "%12.2lf +%12.2lf =%12.2lf", esize, modelstatoverhead, esize+modelstatoverhead);
					GUIPrint(x, x, (float)(wndh>>1)+tdy*4, 2, "%12.2lf", msizes);
				}
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
				display_texture(0,		wndw/3,   0, wndh, txid_separate_r, 1, 0, 1, 0, 1);
				display_texture(wndw/3,		wndw*2/3, 0, wndh, txid_separate_g, 1, 0, 1, 0, 1);
				display_texture(wndw*2/3,	wndw,     0, wndh, txid_separate_b, 1, 0, 1, 0, 1);
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
		if(profileplotmode>PROFILE_OFF)
		{
			void (*draw_profile)(int comp, int color)=profileplotmode==PROFILE_X?draw_profile_x:draw_profile_y;
			draw_profile(0, 0xFF0000FF);
			draw_profile(1, 0xFF00FF00);
			draw_profile(2, 0xFFFF0000);
			draw_profile(3, 0xFF000000);
		}
	}

	if(transforms_customenabled)
	{
		float ystep=tdy*guizoom, x, y;
		if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
		{
			const char chnames[]="rgb";
			unsigned char per[4]={0};

			//custom color transform params
			x=buttons[0].x1;
			y=buttons[0].y1;
			//P0  rgb
			//r += (-0x00*g-0x00*b)>>6
			//g += (-0x00*r-0x00*b)>>6
			//b += (-0x00*r-0x00*g)>>6
			//g += (-0x00*r-0x00*b)>>6
			//0123456789012345678901234
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
			x=buttons[0].x2;
			y=buttons[0].y1;
			{
				double mfwd[9]={0}, minv[9]={0};
				//double mprod[9]={0};

				rct_custom_getmatrix(mfwd, 1);
				rct_custom_getmatrix(minv, 0);
				//for(int ky=0;ky<3;++ky)
				//{
				//	for(int kx=0;kx<3;++kx)
				//	{
				//		double sum=0;
				//		for(int j=0;j<3;++j)
				//			sum+=minv[3*ky+j]*mfwd[3*j+kx];
				//		mprod[ky*3+kx]=sum;
				//	}
				//}

				GUIPrint(0, x, y, guizoom, "Fwd:");
				y+=ystep;
				for(int k=0;k<9-2;k+=3)
				{
					double *v=mfwd+k;
					GUIPrint(0, x, y, guizoom, "%10.6lf %10.6lf %10.6lf", v[0], v[1], v[2]);
					y+=ystep;
				}
				GUIPrint(0, x, y, guizoom, "Inv:");
				y+=ystep;
				for(int k=0;k<9-2;k+=3)
				{
					double *v=minv+k;
					GUIPrint(0, x, y, guizoom, "%10.6lf %10.6lf %10.6lf", v[0], v[1], v[2]);
					y+=ystep;
				}
				//GUIPrint(0, x, y, guizoom, "Product:");
				//y+=ystep;
				//for(int k=0;k<9-2;k+=3)
				//{
				//	double *v=mprod+k;
				//	GUIPrint(0, x, y, guizoom, "%10.6lf %10.6lf %10.6lf", v[0], v[1], v[2]);
				//	y+=ystep;
				//}
			}
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
		if(
			transforms_mask[ST_FWD_CUSTOM]||transforms_mask[ST_INV_CUSTOM]
			||transforms_mask[ST_CONVTEST]
			||transforms_mask[ST_CONVTEST2]
			||transforms_mask[ST_FWD_CC]||transforms_mask[ST_INV_CC]
		//	||transforms_mask[ST_FWD_OLS6]||transforms_mask[ST_INV_OLS6]
		)
		{
			int c0=set_bk_color(0x80FFFFFF);
			int *params=custom_params+CUSTOM_NNB*2*custom_pred_ch_idx;
			int sum=0, neg=0;
			for(int k=0;k<CUSTOM_NNB*2;k+=2)
				sum+=params[k];
			neg=sum<0;
			sum=abs(sum);
			//custom spatial transform params
			x=buttons[1].x1;
			y=buttons[1].y1;
			//0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
			//0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			//-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			GUIPrint(0, x, y, guizoom, "Ch %d  Clamp [%cW %cNW %cN %cNE] sum%c0x%02X.%04X",
				custom_pred_ch_idx,
				custom_clamp[0]?'+':'-',
				custom_clamp[1]?'+':'-',
				custom_clamp[2]?'+':'-',
				custom_clamp[3]?'+':'-',
				neg?'-':' ',
				sum>>16,
				sum&0xFFFF
			);
			//GUIPrint(0, x, y-tdy, 1, "Ch %d", custom_pred_ch_idx);
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
	else if(transforms_mask[ST_FWD_CUSTOM4]||transforms_mask[ST_INV_CUSTOM4]||transforms_mask[ST_FWD_WC]||transforms_mask[ST_INV_WC])
	{
		int c0=set_bk_color(0x80FFFFFF);
		float x=0, y=0, ystep=tdy*guizoom;//
		const short *layer=lossyconv_params+5*5*lossyconv_page;
		int coeffsum=0;
		for(int k=0;k<25;++k)
			coeffsum+=layer[k]>>1;
		x=buttons[6].x1;
		y=buttons[6].y1;
		//Layer 1/4  Ch 1/3
		GUIPrint_append(0, 0, 0, 0, 0, "Layer %d/4  %c %d/%d  cRCT %d  sum %lf",
			(lossyconv_page>>2)+1,
			"RGBA"[lossyconv_page&3],
			(lossyconv_page&3)+1,
			im0->nch,
			lossyconv_causalRCT[lossyconv_page>>2],
			coeffsum/32.
		);
		if(lossyconv_clipboard)
		{
			short val=lossyconv_clipboard, aval=(short)abs(val);
			GUIPrint_append(0, 0, 0, 0, 0, "  M %c%c%c.%c%c%c%c%c%c",
				val<0?'-':' ',
				'0'+(aval>>7&1),
				'0'+(aval>>6&1),
				'0'+(aval>>5&1),
				'0'+(aval>>4&1),
				'0'+(aval>>3&1),
				'0'+(aval>>2&1),
				'0'+(aval>>1&1),
				val&1?'C':' '
			);
			//char clamp=lossyconv_clipboard&1, val=abs(lossyconv_clipboard>>1);
			//GUIPrint_append(0, 0, 0, 0, 0, "  buffer %c%X.%X%c",
			//	lossyconv_clipboard>>7?'-':' ',
			//	val>>4&15,
			//	val&15,
			//	clamp?'C':' '
			//);
		}
		GUIPrint_append(x, x, y+0*ystep, guizoom, 1, 0);
		
		//    000000000011111111112222222222333333
		//    012345678901234567890123456789012345
		//0  "Layer 1/4  C 1/3  cRCT 1"
		//    0000000000111111111122222222223333333333444444444455555555556
		//    0123456789012345678901234567890123456789012345678901234567890
		//1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		//2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		//3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
		//4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		//5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		//6  "XY  Offset XXX YYY"
		//7  "    Stride XXX YYY"
		for(int ky=0;ky<5;++ky)
		{
			g_printed=0;
			for(int kx=0;kx<5;++kx)
			{
				short val=layer[5*ky+kx], aval=(short)abs(val);
				char lbr=' ', rbr=' ';
				if(kx==2&&ky==2)
					lbr='[', rbr=']';
				GUIPrint_append(0, 0, 0, 0, 0, "%c%c%c%c.%c%c%c%c%c%c%c",
					lbr,
					val<0?'-':' ',
					'0'+(aval>>7&1),
					'0'+(aval>>6&1),
					'0'+(aval>>5&1),
					'0'+(aval>>4&1),
					'0'+(aval>>3&1),
					'0'+(aval>>2&1),
					'0'+(aval>>1&1),
					val&1?'C':' ',
					rbr
				);
#if 0
				const char *p=layer+5*ky+kx;
				char clamp=*p&1, val=abs(*p>>1);
				char lbr, rbr;
				if(kx==2&&ky==2)
					lbr='[', rbr=']';
				else
					lbr=' ', rbr=' ';
				GUIPrint_append(0, 0, 0, 0, 0, "%c%c%X.%X%c%c", lbr, p[0]>>7?'-':' ', val>>4&15, val&15, clamp?'C':' ', rbr);
#endif
			}
			GUIPrint_append(x, x, y+(ky+1)*ystep, guizoom, 1, 0);
		}

		GUIPrint(x, x, y+6*ystep, guizoom, "XY  Offset %3d %3d", lossyconv_offset[lossyconv_page>>2<<1|0], lossyconv_offset[lossyconv_page>>2<<1|1]);
		GUIPrint(x, x, y+7*ystep, guizoom, "    Stride %3d %3d", lossyconv_stride[lossyconv_page>>2<<1|0]+1, lossyconv_stride[lossyconv_page>>2<<1|1]+1);
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
	double usize=0, csize=0;
	if(im1)
	{
		float
			x, y, x2, xstart, xend, ystart, scale,
			crformat, cr_combined;
		int RGBspace;

		x=(float)(wndw-300), y=tdy*2, x2=x;
		for(int k=0;k<T_COUNT;++k)//print available transforms on right
		{
			transforms_printname(x2, y, k, -1, transforms_mask[k]?0xA0FF0000FFFFFFFF:0);
			x2=x+(150&-!(k&1));
			if(k&1)
				y+=tdy*g_uiscale;
		}
		x=(float)(wndw-450);
		y=tdy*2;
		{
			const char *label=ec_method_label(ec_method);
			if(ec_method==ECTX_ABAC0||ec_method==ECTX_ABAC1)
				GUIPrint(x, x, y-tdy, 1, "H - - -  -         %s", label);
			else if(ec_adaptive)//H.E.M.L..A.0x0000..XXXX_XXX
				GUIPrint(x, x, y-tdy, 1, "H %d %d %d  A 0x%04X  %s", ec_expbits, ec_msb, ec_lsb, ec_adaptive_threshold, label);
			else
				GUIPrint(x, x, y-tdy, 1, "H %d %d %d  Static    %s", ec_expbits, ec_msb, ec_lsb, label);
		}
		if(transforms)
		{
			for(int k=0;k<(int)transforms->count;++k, y+=tdy)//print applied transforms on left
				transforms_printname(x, y, transforms->data[k], k, 0);
		}
		cr_combined=(float)((im1->src_depth[0]+im1->src_depth[1]+im1->src_depth[2]+im1->src_depth[3])/(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3]));
		xstart=20, xend=(float)wndw-330, ystart=(float)(wndh-tdy*5);
		scale=xend-xstart;
		
		usize=image_getBMPsize(im0);
		crformat=(float)(usize/filesize);
		RGBspace=1;
		for(int k=0;k<CST_COMPARE;++k)
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
			if(maxinvcr>100)
				maxinvcr=100;
			if(maxinvcr<1/crformat)
				maxinvcr=1/crformat;
			if(xend-scale*maxinvcr<xstart)
				scale=(xend-xstart)/maxinvcr;
			if(scale<50)
				scale=50;

			xstart=xend-scale*ceilf((maxinvcr+1)*2)*0.5f;
			draw_rect(xstart, xend, ystart, (float)wndh, 0x80808080);//background
			{
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
						draw_line(x, ystart+1, x, (float)wndh, 0x70900090);//draw minor scale
						x=(float)(xend-0.1f*scale*ks);
					}
				}
				ks=1;
				x=(float)(xend-scale*ks);
				for(int ks2=1;x>xstart;++ks2)
				{
					draw_rect(x-1, x+2, ystart, (float)wndh, 0x70800080);//draw major scale
					x=(float)(xend-scale*ks2);
				}
			}
			{
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
			}
			x=xend-scale/crformat;
			draw_line(x, ystart, x-10, ystart-10, 0xFF000000);
			draw_line(x, ystart, x+10, ystart-10, 0xFF000000);
		}
		{
			int prevtxtcolor, prevbkcolor;
			xend+=10;
			prevbkcolor=set_bk_color(0xC0C0C0C0);
			prevtxtcolor=set_text_color(0xFF000000);
			GUIPrint(xend, xend, ystart-tdy*2, 1, "Bitmap Size             %9.0lf", usize);
			//GUIPrint(xend, xend, ystart-tdy*2, 1, "Uncompressed Size       %9.0lf", usize);
			GUIPrint(xend, xend, ystart-tdy  , 1, "Format        %8.4f%% %9lld", 100./crformat, filesize);
			set_bk_color(0xE0FFFFFF);
			csize=usize/cr_combined;
			prevtxtcolor=set_text_color(0xFF000000);GUIPrint(xend, xend, ystart      , 1, "Combined      %8.4f%% %12.2lf", 100./cr_combined, csize);
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
		}

		//if(transforms_customenabled)
		//{
		//	//double maxloss=0;
		//	GUIPrint(0, 0, tdy*3, 1, "RMSE %lf", av_rmse);
		//	if(minloss<maxloss)
		//		GUIPrint(0, 200, tdy*3, 1, "[%lf~%lf]", minloss, maxloss);
		//}
		{
			float g2, cx, cy;
			int idx, idx2;

			g2=wndh/combCRhist_max;
			idx=combCRhist_idx-1;
			idx2=combCRhist_idx-2;
			idx+=combCRhist_SIZE&-(idx<0);
			idx2+=combCRhist_SIZE&-(idx2<0);
			xstart=(float)(wndw-(combCRhist_SIZE<<combCRhist_logDX)-300);
			cx=xstart+(float)(idx<<combCRhist_logDX);
			cy=wndh-combCRhist[idx][3]*g2;
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
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][0]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][0]*g2), RGBspace?0xFF0000FF:0xFF404040);//r or Y
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][1]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][1]*g2), RGBspace?0xFF00FF00:0xFFC00000);//g or Cb
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][2]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][2]*g2), RGBspace?0xFFFF0000:0xFF0000C0);//b or Cr
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][3]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][3]*g2), 0xC0000000);
				}
			}
			draw_line(xstart, wndh-g2, xstart+(combCRhist_SIZE<<combCRhist_logDX), wndh-g2, 0xC0000000);
		}
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
				print_i16_row(x, y, zoom, c2_params.c0+k*12, 5);	y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c0+k*12+5, 5);	y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c0+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}

			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16_row(x, y, zoom, c2_params.c1+k*12, 5);	y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c1+k*12+5, 5);	y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c1+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}
			print_i16_row(x, y, zoom, c2_params.c1+6*12, 2);		y+=tdy*zoom;
			
			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16_row(x, y, zoom, c2_params.c2+k*12, 5);	y+=tdy*zoom;
				print_i16_row(x, y, zoom, c2_params.c2+k*12+5, 5);	y+=tdy*zoom;
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

			x=(float)(wndw>>1);
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
			x=(float)(wndw>>1);
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
		GUIPrint(0, 0, 0, 1, "WH %dx%d  D0[%d %d %d %d] D[%d %d %d %d]  RCT%2d  Z %13.6lf",
			im0->iw, im0->ih,
			im0->src_depth[0], im0->src_depth[1], im0->src_depth[2], im0->src_depth[3],
			im1->depth[0], im1->depth[1], im1->depth[2], im1->depth[3],
			im1->rct,
			imzoom
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
	{
		const char *mode_str=0, *view_str=0;
		static double t=0;
		double t2=time_ms();

		switch(mode)
		{
	//	case VIS_PLANES:		mode_str="Planes";		break;
	//	case VIS_MESH:			mode_str="Combined Mesh";	break;
	//	case VIS_MESH_SEPARATE:		mode_str="Separate Mesh";	break;
		case VIS_ANALYSIS:		mode_str="Analysis";		break;
		case VIS_HISTOGRAM:		mode_str="Histogram";		break;
		case VIS_MODEL:			mode_str="Model";		break;
		case VIS_JOINT_HISTOGRAM:	mode_str="Joint Histogram";	break;
#ifdef ENABLE_L1WEIGHTS
		case VIS_L1WEIGHTS:		mode_str="L1 Weights";		break;
#endif
		case VIS_IMAGE:			mode_str="Image View";		break;
	//	case VIS_BAYES:			mode_str="Bayes";		break;
		case VIS_ZIPF:			mode_str="Zipf View";		break;
	//	case VIS_IMAGE_BLOCK:		mode_str="Image Block";		break;
	//	case VIS_IMAGE_E24:		mode_str="Image Exp24";		break;
	//	case VIS_DWT_BLOCK:		mode_str="DWT Block";		break;
		case VIS_IMAGE_TRICOLOR:	mode_str="Tricolor";		break;
		}
		GUIPrint(0, 0, tdy, 1, "timer %d, fps%11lf, [%2d/%2d] %s  dist=%d  BPD %6.4lf"
			, timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str, g_dist, usize?8*csize/usize:0
		);
		if(mode==VIS_IMAGE||mode==VIS_ZIPF)
		{
			switch(viewmode)
			{
			case VIEW_RGB:	view_str="RGB";break;
			case VIEW_C0:	view_str="C0";break;
			case VIEW_C1:	view_str="C1";break;
			case VIEW_C2:	view_str="C2";break;
			}
			GUIPrint(0, 0, tdy*2, 1, "%s %s", show_full_image?"FILL SCREEN":"1:1", view_str);
			if(ec_method==ECTX_ABAC0)
			{
				char bin[]=
				{
					'0'+(abacvis_low>>15&1),
					'0'+(abacvis_low>>14&1),
					'0'+(abacvis_low>>13&1),
					'0'+(abacvis_low>>12&1),
					'0'+(abacvis_low>>11&1),
					'0'+(abacvis_low>>10&1),
					'0'+(abacvis_low>> 9&1),
					'0'+(abacvis_low>> 8&1),
					'0'+(abacvis_low>> 7&1),
					'0'+(abacvis_low>> 6&1),
					'0'+(abacvis_low>> 5&1),
					'0'+(abacvis_low>> 4&1),
					'0'+(abacvis_low>> 3&1),
					'0'+(abacvis_low>> 2&1),
					'0'+(abacvis_low>> 1&1),
					'0'+(abacvis_low>> 0&1),
					0,
				};
				GUIPrint(0, 0, tdy*3, 2, "0b%s/%d", bin, abacvis_range);
			//	GUIPrint(0, 0, tdy*3, 2, "0x%04X/%d", abacvis_low, abacvis_range);
			}
		}
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
	}
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