#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
#include"lodepng.h"
#include"c18.h"
static const char file[]=__FILE__;


	#define ENABLE_L1WEIGHTS
	#define USE_OLS


static const char *rct_names[]=
{
#define RCT(LABEL, OCH0,OCH1,OCH2,P0,P1,P2,U0,V0,V1) #LABEL,
	RCTLIST
#undef  RCT
};

/*c18*/
#if 1
static const char *och_names[OCH2_COUNT]=
{
#define OCH2(LABEL) #LABEL,
	OCH2LIST
#undef  OCH2
};
static const char *rct2_names[RCT2_COUNT]=
{
#define RCT2(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) #LABEL,
	RCT2LIST
#undef  RCT2
};
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL, DESC) #LABEL,
	PREDLIST
#undef  PRED
};
static const char *pred_desc[]=
{
#define PRED(LABEL, DESC) DESC,
	PREDLIST
#undef  PRED
};
#endif

#ifdef ENABLE_L1WEIGHTS
	#define L1DELAY 20
	#define L1SPEED 1
#if 1
#define L1SH 19
#define L1BIAS0	0
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
	L1PRED(100000, NEEE)\

#endif
enum
{
#define L1PRED(W, E) +1
	L1NPREDS=L1PREDLIST
#undef  L1PRED
};
static const char *l1prednames[]=
{
#define L1PRED(W0, EXPR) #EXPR,
	L1PREDLIST
#undef  L1PRED
};
static int weights_l1[L1NPREDS+1]=
{
#define L1PRED(W0, EXPR) (1<<L1SH)/L1NPREDS,
	L1PREDLIST
#undef  L1PRED
	L1BIAS0,
};
static int weights_wp[L1NPREDS+1]={0};
static int weights_ols[L1NPREDS+1]={0};
#ifdef USE_OLS
static double curr_cov[L1NPREDS*L1NPREDS]={0}, curr_vec[L1NPREDS]={0}, curr_cholesky[L1NPREDS*L1NPREDS]={0}, curr_params[L1NPREDS]={0};
#endif
static int wperrors[L1NPREDS]={0};
static int mixmode=0;/*0: L1    2: WP    1: OLS*/
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

ArrayHandle currfn=0;
size_t filesize=0;

Image *im0=0, *im1=0;/*im0: original image, im1: image with selected transforms*/
int pred_ma_enabled=1;/*modular arithmetic for spatial predictors*/
int separate_grayscale=1;/*separate channels are shown greyscale*/
uint32_t txid_separate_r=0, txid_separate_g=0, txid_separate_b=0, txid_separate_a=0;
uint8_t *im_export=0, *zimage=0;

uint32_t txid_jointhist[256]={0};

typedef enum VisModeEnum
{
	VIS_IMAGE_TRICOLOR,
#ifdef ENABLE_L1WEIGHTS
	VIS_L1WEIGHTS,
#endif
	VIS_IMAGE,
	VIS_HISTOGRAM,
	VIS_MODEL,
	VIS_ANALYSIS,
	VIS_JOINT_HISTOGRAM,
	VIS_SHANNON,

	VIS_COUNT
} VisMode;
int mode=VIS_IMAGE;

typedef enum _ImageViewMode
{
	VIEW_RGB,
	VIEW_C0,
	VIEW_C1,
	VIEW_C2,

	VIEW_COUNT
} ImageViewMode;
int viewmode=VIEW_RGB, view_ma=0;

typedef enum _HistMode
{
	HISTMODE_LINEAR,
	HISTMODE_LINEAR_NOSPIKE,
	HISTMODE_LOG2,
	HISTMODE_LOG2_NOSPIKE,

	HISTMODE_COUNT
} HistMode;
int histmode=HISTMODE_LINEAR;
ArrayHandle rectangles=0;

typedef enum TransformTypeEnum
{
	CT_FWD_SubG_OPT,	CT_INV_SubG_OPT,
	CT_FWD_SUBGREEN,	CT_INV_SUBGREEN,
	CT_FWD_JPEG2000,	CT_INV_JPEG2000,/*	(1997) JPEG2000 RCT*/
	CT_FWD_JPEG2000_MA,	CT_INV_JPEG2000_MA,
	CT_FWD_NBLI,		CT_INV_NBLI,	/*	(2024) NBLI*/
	CT_FWD_YCoCg_R,		CT_INV_YCoCg_R,	/*	(2003) AVC, HEVC, VVC*/
	CT_FWD_YCbCr_R_v1,	CT_INV_YCbCr_R_v1,
	CT_FWD_YCbCr_R_v2,	CT_INV_YCbCr_R_v2,
	CT_FWD_YCbCr_R_v3,	CT_INV_YCbCr_R_v3,
	CT_FWD_YCbCr_R_v4,	CT_INV_YCbCr_R_v4,
	CT_FWD_YCbCr_R_v5,	CT_INV_YCbCr_R_v5,
	CT_FWD_YCbCr_R_v6,	CT_INV_YCbCr_R_v6,
	CT_FWD_YCbCr_R_v7,	CT_INV_YCbCr_R_v7,
	CT_FWD_Pei09,		CT_INV_Pei09,
#if 0
	CT_FWD_J2K2,		CT_INV_J2K2,
	CT_FWD_CrCgCb,		CT_INV_CrCgCb,
	CT_FWD_YRGB_v1,		CT_INV_YRGB_v1,
	CT_FWD_YRGB_v2,		CT_INV_YRGB_v2,
	CT_FWD_CMYK,		CT_INV_CMYK_DUMMY,
	CT_FWD_MATRIX,		CT_INV_MATRIX,
#endif
	CT_FWD_YCbCr,		CT_INV_YCbCr,	/*LOSSY	JPEG		*/
	CT_FWD_XYB,		CT_INV_XYB,	/*LOSSY	(2021) JPEG XL	*/
	CT_FWD_CUSTOM,		CT_INV_CUSTOM,
#if 0
	CT_FWD_PCA,		CT_INV_PCA,
#endif
	CT_FWD_HISTMATCH,	CT_INV_HISTMATCH,
#if 0
	CT_FWD_ADAPTIVE,	CT_INV_ADAPTIVE,
#endif

	CST_COMPARE,		CST_INV_SEPARATOR,
	
	ST_FWD_MIXR,		ST_INV_MIXR,
	ST_FWD_MIXN,		ST_INV_MIXN,
#if 0
	ST_FWD_MIXNC,		ST_INV_MIXNC,
	ST_FWD_BESTN,		ST_INV_BESTN,
	ST_FWD_RLS,		ST_INV_RLS,
#endif
	ST_FWD_ADAQUANT,	ST_INV_ADAQUANT,
	ST_FWD_GRAY,		ST_INV_GRAY,
	ST_FWD_GRFILT,		ST_INV_GRFILT,
	ST_FWD_L1CRCT,		ST_INV_L1CRCT,
	ST_FWD_OLS7,		ST_INV_OLS7,
	ST_FWD_L1BCRCT,		ST_INV_L1BCRCT,
	ST_FWD_OLS8,		ST_INV_OLS8,
	ST_FWD_OLS9,		ST_INV_OLS9,
	ST_FWD_PACKSIGN,	ST_INV_PACKSIGN,
	ST_FWD_SIGNBIT,		ST_FWD_MAGNITUDE,
#if 0
	ST_FWD_LPC,		ST_INV_LPC,
#endif
	ST_FWD_CGCRCT,		ST_INV_CGCRCT,
	ST_FWD_CLAMPGRAD,	ST_INV_CLAMPGRAD,
	ST_FWD_SELECT,		ST_INV_SELECT,
	ST_FWD_PAETH,		ST_INV_PAETH,
	ST_FWD_SELECT4,		ST_INV_SELECT4,
	ST_FWD_AV2,		ST_INV_AV2,
	ST_FWD_SUB,		ST_INV_SUB,
	ST_FWD_CALIC,		ST_INV_CALIC,
	ST_CONVTEST,		ST_CONVTEST2,
	ST_DIFF,		ST_SSIM,
	ST_FILT_MEDIAN33,	ST_FILT_AV33,
	ST_FILT_DEINT422,	ST_FILT_DEINT420,
	ST_FWD_EXTRAP,		ST_FWD_L1DITHER,
	ST_FWD_QUANT,		ST_INV_QUANT,

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
#if 0
	ST_FWD_WGRAD2,		ST_INV_WGRAD2,
#endif
	ST_FWD_WGRAD3,		ST_INV_WGRAD3,
	ST_FWD_WGRAD4CCRCT,	ST_INV_WGRAD4CCRCT,
	ST_FWD_WGRAD4C,		ST_INV_WGRAD4C,
	ST_FWD_WGRAD4,		ST_INV_WGRAD4,
	ST_FWD_WGRAD5,		ST_INV_WGRAD5,
	ST_FWD_WGRAD6,		ST_INV_WGRAD6,
	ST_FWD_WGRAD7,		ST_INV_WGRAD7,
	ST_FWD_SSE,		ST_INV_SSE,
#if 0
	ST_FWD_WMIX,		ST_INV_WMIX,
#endif
	ST_FWD_T47,		ST_INV_T47,
	ST_FWD_P3,		ST_INV_P3,
	ST_FWD_G2,		ST_INV_G2,
	ST_FWD_MM,		ST_INV_MM,
#if 0
	ST_FWD_WP2,		ST_INV_WP2,
#endif
	ST_FWD_WPU,		ST_INV_WPU,
#if 0
	ST_FWD_DEFERRED,	ST_INV_DEFERRED,
#endif
	ST_FWD_WC,		ST_INV_WC,	/*irreversible conv*/
	ST_FWD_CUSTOM4,		ST_INV_CUSTOM4,	/*irreversible conv*/
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
#if 0
	ST_FWD_MIX2,		ST_INV_MIX2,
	ST_FWD_AV3,		ST_INV_AV3,
	ST_FWD_ECOEFF,		ST_INV_ECOEFF,
	ST_FWD_AVERAGE,		ST_INV_AVERAGE,
	ST_FWD_MULTISTAGE,	ST_INV_MULTISTAGE,
#endif
	ST_FWD_ZIPPER,		ST_INV_ZIPPER,
#if 0
	ST_FWD_DIR,		ST_INV_DIR,
	ST_FWD_CTX,		ST_INV_CTX,
	ST_FWD_C20,		ST_INV_C20,
	ST_FWD_WU97,		ST_INV_WU97,
	ST_FWD_SHUFFLE,		ST_INV_SHUFFLE,
	ST_FWD_SPLIT,		ST_INV_SPLIT,
#endif
	
	ST_FWD_LAZY,		ST_INV_LAZY,
	ST_FWD_HAAR,		ST_INV_HAAR,
	ST_FWD_SQUEEZE,		ST_INV_SQUEEZE,
	ST_FWD_LEGALL53,	ST_INV_LEGALL53,
	ST_FWD_CDF97,		ST_INV_CDF97,
	ST_FWD_AWAV,		ST_INV_AWAV,
#if 0
	ST_FWD_EXPDWT,		ST_INV_EXPDWT,
	ST_FWD_CUSTOM_DWT,	ST_INV_CUSTOM_DWT,
#endif

	ST_FWD_DCT4,		ST_INV_DCT4,
	ST_FWD_DCT8,		ST_INV_DCT8,
	ST_FWD_WHT4,		ST_INV_WHT4,
	ST_FWD_WHT8,		ST_INV_WHT8,
	ST_FWD_HAAR8,		ST_INV_HAAR8,
	ST_FWD_FDCT8,		ST_INV_FDCT8,

	T_COUNT
} TransformType;
int transforms_customenabled=0;
char transforms_mask[T_COUNT]={0};
ArrayHandle transforms=0;/*array of chars*/
float guizoom=1.25f;

double av_rmse=0, g_lr=1e-10;

int profile_idx=0;
double minloss=0, maxloss=0;

int blocksize=16, margin=32;
int blockmx=0, blockmy=0;

float
	pixel_amplitude=10,/*4*/
	mesh_separation=100;/*10*/
int extrainfo=0;

ArrayHandle cpu_vertices=0;
uint32_t gpu_vertices=0;

C18Info analysis_info={0};
static int64_t pred_hist[3][PRED_COUNT]={0};

ArrayHandle jointhist=0;
int jointhist_nbits=6;/*max*/
int jhx=0, jhy=0;
double ch_entropy[4]={0};/*RGBA/YUVA	[0 ~ 1]*/
/*int usage[4]={0};*/
EContext ec_method=ECTX_GRCTX;
#if 0
EContext ec_method=ECTX_HIST;/*ECTX_MIN_QN_QW;*/
#endif
int ec_adaptive=0, ec_adaptive_threshold=3200, ec_expbits=5, ec_msb=2, ec_lsb=0;
int abacvis_low=0, abacvis_range=0;

enum
{
	MODELNCTX=17,
	MODELPREC=3
};
int modelnch=0, modelnctx=0, modeldepth=0, modelhistsize=0, *modelhist=0;
double modelcsizes[4*MODELNCTX]={0};
double modelmeans[4*MODELNCTX]={0}, modelsdevs[4*MODELNCTX]={0};
int *modelmhist=0;
double modelmsizes[4*MODELNCTX]={0};
double modelstatoverhead=0;
int modelvis=0;

#define combCRhist_SIZE 128
#define combCRhist_logDX 2
float combCRhist[combCRhist_SIZE][4]={0}, combCRhist_max=1;
int combCRhist_idx=0;

int imagecentered=0;
double
	imzoom=1,/*image pixel size in screen pixels*/
	wpx=0, wpy=0;/*window position (top-left corner) in image coordinates*/
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
	PRFILE_COUNT
} ProfilePlotMode;
ProfilePlotMode profileplotmode=PROFILE_OFF;
int show_full_image=0;
int space_not_color=0;
#define ZOOM_LIMIT_LABEL 48
#define ZOOM_LIMIT_ALPHA 96
int pxlabels_hex=0;

/*joint hist box contour*/
float jh_cubesize=64;
int jhc_boxdx=64, jhc_boxdy=64;
int jhc_xbox=0, jhc_ybox=0;/*[0, screen_dim-1]*/
ArrayHandle jhc_mesh=0;
uint32_t jhc_gpubuf=0;
float jhc_level=10.5f;
int loud_transforms=1;
int g_dist=1;
float g_uiscale=1;
int drawCRhistory=0;

static double invcr2csizes(const double *invcr, const int8_t *src_depths, int iw, int ih, int nch, double *csizes)/*invcr & csizes have nch+1=5 elements*/
{
	double usize=0;
	int dtotal=0, k;

	for(k=0;k<nch;++k)
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
	wndw2=wndw;
	wndh2=wndh-(int)(tdy*5);
	if(wndh2<0)
		wndh2=wndh;
	if(wndh2<0)
		return;
	/*wndh2=wndh-17;*/
	if((double)wndw2/wndh2>=(double)im1->iw/im1->ih)/*window AR > image AR: fit height*/
	{
		if(wndh2>0)
			imzoom=(double)wndh2/im1->ih;
	}
	else/*window AR < image AR: fit width*/
		imzoom=(double)wndw2/im1->iw;
	wpx=(im1->iw-wndw2/imzoom)*0.5;/*center image*/
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

void calc_csize_stateful(Image const *image, int *hist_full, double *entropy)
{
	if(ec_method==ECTX_C2)
	{
		enum
		{
			C2_NCTX=24,

			XPAD=8,
			NROWS=4,
			NCH=4,
			NVAL=1
		};
		uint8_t sh[4];
		int halfs[4], masks[4];
		int psize;
		int32_t *pixels;
		int nch, maxdepth, nlevels;
		int hsize;
		uint32_t *hists;
		int kx, ky, kc, idx, kctx;
		double ctxsizes[4]={0};

		sh[0]=32-image->depth[0];
		sh[1]=32-image->depth[1];
		sh[2]=32-image->depth[2];
		sh[3]=32-image->depth[3];
		halfs[0]=1<<image->depth[0]>>1;
		halfs[1]=1<<image->depth[1]>>1;
		halfs[2]=1<<image->depth[2]>>1;
		halfs[3]=1<<image->depth[3]>>1;
		masks[0]=(1<<image->depth[0])-1;
		masks[1]=(1<<image->depth[1])-1;
		masks[2]=(1<<image->depth[2])-1;
		masks[3]=(1<<image->depth[3])-1;
		psize=(image->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
		pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
		nch=(image->depth[0]!=0)+(image->depth[1]!=0)+(image->depth[2]!=0)+(image->depth[3]!=0);
		maxdepth=image->depth[0];
		if(maxdepth<image->depth[1])maxdepth=image->depth[1];
		if(maxdepth<image->depth[2])maxdepth=image->depth[2];
		if(maxdepth<image->depth[3])maxdepth=image->depth[3];
		nlevels=1<<maxdepth;
		hsize=(int)sizeof(uint32_t[C2_NCTX])*nch<<maxdepth;
		hists=(uint32_t*)malloc(hsize);
		if(!hists||!pixels)
		{
			CRASH("Alloc error");
			return;
		}
		memset(hists, 0, hsize);
		memset(pixels, 0, psize);
		for(ky=0, idx=0;ky<image->ih;++ky)
		{
			ALIGN(32) int32_t *rows[4]={0};
			int32_t xnoise[4]={0};

			rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
			rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
			rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
			rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
			for(kx=0;kx<image->iw;++kx)
			{
				for(kc=0;kc<4;++kc, ++idx)
				{
					int ctx, sym;
					int32_t eN;

					rows[0]+=NROWS*NVAL;
					rows[1]+=NROWS*NVAL;
					rows[2]+=NROWS*NVAL;
					rows[3]+=NROWS*NVAL;
					if(!image->depth[kc])
						continue;
					eN	=rows[1][0+0*NCH*NROWS*NVAL];

					ctx=eN<xnoise[kc]?eN:xnoise[kc];
					ctx=31^LZCNT32(ctx*ctx+1);
					if(ctx>C2_NCTX/2-1)
						ctx=C2_NCTX/2-1;
					ctx+=C2_NCTX/2*(eN<xnoise[kc]);
					sym=((const int32_t*)HDR_DATC(image))[idx];
					sym<<=sh[kc];
					sym>>=sh[kc];
					++hists[(kc*C2_NCTX+ctx)<<maxdepth|((sym+halfs[kc])&masks[kc])];
					sym=sym<<1^sym>>31;
					xnoise[kc]+=((sym<<MODELPREC)-xnoise[kc])>>3;
					rows[0][0]=eN+(((sym<<MODELPREC)-eN)>>3);
				}
			}
		}
		_mm_free(pixels);
		for(kc=0;kc<nch;++kc)
		{
			double chsize=0;
			for(kctx=0;kctx<C2_NCTX;++kctx)
			{
				double e;
				uint32_t *curr_hist;
				int ctx, ks, ks2, sum;

				ctx=kc*C2_NCTX+kctx;
				curr_hist=hists+((intptr_t)ctx<<maxdepth);
				sum=0;
				for(ks=0;ks<nlevels;++ks)
					sum+=curr_hist[ks];
				if(!sum)
					continue;
				e=0;
				if(maxdepth>9)/*exact estimate*/
				{
					double norm=1./sum;
					for(ks=0;ks<nlevels;++ks)
					{
						int freq=curr_hist[ks];
						if(freq)
							e-=freq*log(freq*norm);
					}
					chsize+=e*(1./(8*M_LN2));
				}
				else
				{
					/*simulate bin tracking guard*/
					const int probbits=12;
					int count=0;
					double hsize2=0;
					int cdfW=0;
					int sum2=0;
					int CDFlevels, nlevels2, half2, mask2;

					for(ks=0;ks<nlevels;++ks)
						count+=curr_hist[ks]!=0;
					for(ks=0;ks<nlevels;++ks)
					{
						int freq=curr_hist[ks];
						if(freq)
						{
							int prob=(int)((int64_t)freq*((int64_t)0x1000-count)/sum)+1;
							e-=freq*log(prob*(1./0x1000));
						}
					}
					e*=1./(8*M_LN2);
					if(e<0||acme_isnan(e))
						CRASH("");
					chsize+=e;
					CDFlevels=1<<probbits;
					nlevels2=1<<image->depth[kc], half2=nlevels2>>1, mask2=nlevels2-1;
					for(ks=0, ks2=0;ks<nlevels2;++ks)/*calc model stat overhead*/
					{
						int sym, freq, cdf, csym;

						sym=((ks>>1^-(ks&1))+half2)&mask2;
						freq=curr_hist[sym];
						cdf=sum2*(((uint64_t)1<<probbits)-count)/sum+ks2;
						ks2+=freq!=0;
						csym=cdf-cdfW;
						if(ks&&CDFlevels)/*CDF[0] is always zero*/
						{
							/*GR*/
							int nbypass=31^LZCNT32(CDFlevels);
							if(ks>1)
								nbypass-=7;
							if(nbypass<0)
								nbypass=0;
							hsize2+=(csym>>nbypass)+1+nbypass;
						}
						CDFlevels-=csym;
						cdfW=cdf;
						sum2+=freq;
					}
					chsize+=hsize2/8;
					ctxsizes[kc]+=hsize2/8;
				}
			}
			entropy[kc]=image->src_depth[kc]?chsize*8/((double)image->iw*image->ih*image->src_depth[kc]):0;
		}
		(void)ctxsizes;
		free(hists);
	}
	else if(ec_method==ECTX_GRCTX)
	{
		enum
		{
			XPAD=8,
			NROWS=4,
			NCH=4,
			NVAL=1
		};
#if 0
		int64_t bitsizes[4*16*2]={0};
#endif
		const int nctx=MODELNCTX;
		int nch, maxdepth, nlevels;
		int hsize, *hists;
		int psize;
		int32_t *pixels;
		uint8_t sh[4]={0};
		int halfs[4]={0}, masks[4]={0};
		int kx, ky, kc, idx;
		double ctxsizes[4]={0};

		nch=(image->depth[0]!=0)+(image->depth[1]!=0)+(image->depth[2]!=0)+(image->depth[3]!=0);
		maxdepth=image->depth[0];
		if(maxdepth<image->depth[1])maxdepth=image->depth[1];
		if(maxdepth<image->depth[2])maxdepth=image->depth[2];
		if(maxdepth<image->depth[3])maxdepth=image->depth[3];
		nlevels=1<<maxdepth;
		hsize=(int)sizeof(int)*nch*nctx<<maxdepth;
		hists=(int*)malloc(hsize);
		psize=(image->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
		pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
		if(!hists||!pixels)
		{
			CRASH("Alloc error");
			return;
		}
		memset(hists, 0, hsize);
		memset(pixels, 0, psize);
		sh[0]=32-image->depth[0];
		sh[1]=32-image->depth[1];
		sh[2]=32-image->depth[2];
		sh[3]=32-image->depth[3];
		halfs[0]=1<<image->depth[0]>>1;
		halfs[1]=1<<image->depth[1]>>1;
		halfs[2]=1<<image->depth[2]>>1;
		halfs[3]=1<<image->depth[3]>>1;
		masks[0]=(1<<image->depth[0])-1;
		masks[1]=(1<<image->depth[1])-1;
		masks[2]=(1<<image->depth[2])-1;
		masks[3]=(1<<image->depth[3])-1;
		for(ky=0, idx=0;ky<image->ih;++ky)
		{
			int32_t *rows[4];

			rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
			rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
			rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
			rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
			for(kx=0;kx<image->iw;++kx)
			{
				for(kc=0;kc<4;++kc, ++idx)
				{
					rows[0]+=NROWS*NVAL;
					rows[1]+=NROWS*NVAL;
					rows[2]+=NROWS*NVAL;
					rows[3]+=NROWS*NVAL;
					if(image->depth[kc])
					{
						int32_t
							eNEE	=rows[1][0+2*NCH*NROWS*NVAL],
							eNEEE	=rows[1][0+3*NCH*NROWS*NVAL],
							eW	=rows[0][0-1*NCH*NROWS*NVAL];
						int ctx, sym;

						ctx=31^LZCNT32(eW*eW+1);
						if(ctx>MODELNCTX-1)
							ctx=MODELNCTX-1;
						sym=((const int32_t*)HDR_DATC(image))[idx];
						sym<<=sh[kc];
						sym>>=sh[kc];
						++hists[(kc*nctx+ctx)<<maxdepth|((sym+halfs[kc])&masks[kc])];
						sym=sym<<1^sym>>31;
						rows[0][0]=(2*eW+(sym<<MODELPREC)+MAXVAR(eNEE, eNEEE))>>2;
					}
				}
			}
		}
		_mm_free(pixels);
		for(kc=0;kc<nch;++kc)
		{
			double chsize=0;
			int kctx=0, ks, ks2;

			for(kctx=0;kctx<nctx;++kctx)
			{
				int ctx=kc*nctx+kctx;
				int *curr_hist=hists+((ptrdiff_t)ctx<<maxdepth);
				int sum=0;
				double e=0;

				for(ks=0;ks<nlevels;++ks)
					sum+=curr_hist[ks];
				if(!sum)
					continue;
				if(maxdepth>9)/*exact estimate*/
				{
					double norm=1./sum;
					for(ks=0;ks<nlevels;++ks)
					{
						int freq=curr_hist[ks];
						if(freq)
							e-=freq*log(freq*norm);
					}
					chsize+=e*(1./(8*M_LN2));
				}
				else
				{
					/*simulate bin tracking guard*/
					const int probbits=12;
					int CDFlevels=1<<probbits;
					int nlevels2=1<<image->depth[kc], half2=nlevels2>>1, mask2=nlevels2-1;
					double hsize2=0;
					int count=0;
					int cdfW=0;
					int sum2=0;

					for(ks=0;ks<nlevels;++ks)
						count+=curr_hist[ks]!=0;
					for(ks=0;ks<nlevels;++ks)
					{
						int freq=curr_hist[ks];
						if(freq)
						{
							int prob=(int)((int64_t)freq*((uint64_t)0x1000-count)/sum)+1;
							e-=freq*log(prob*(1./0x1000));
						}
					}
					if(e<0||acme_isnan(e))
						CRASH("");
					e*=1./(8*M_LN2);
					chsize+=e;
					for(ks=0, ks2=0;ks<nlevels2;++ks)/*calc model stat overhead*/
					{
						int sym, freq, cdf, csym;

						sym=((ks>>1^-(ks&1))+half2)&mask2;
						freq=curr_hist[sym];
						cdf=sum2*(((uint64_t)1<<probbits)-count)/sum+ks2;
						ks2+=freq!=0;
						csym=cdf-cdfW;
						if(ks&&CDFlevels)/*CDF[0] is always zero*/
						{
							/*GR*/
							int nbypass=31^LZCNT32(CDFlevels);
							if(ks>1)
								nbypass-=7;
							if(nbypass<0)
								nbypass=0;
							hsize2+=(csym>>nbypass)+1+nbypass;
						}
						CDFlevels-=csym;
						cdfW=cdf;
						sum2+=freq;
					}
					chsize+=hsize2/8;
					ctxsizes[kc]+=hsize2/8;
				}
			}
			entropy[kc]=image->src_depth[kc]?chsize*8/((double)image->iw*image->ih*image->src_depth[kc]):0;
		}
		(void)ctxsizes;
		free(hists);
	}
	else if(ec_method==ECTX_GR)
	{
		int64_t bitsizes[4]={0};
		int bufsize;
		int16_t *pixels;
		int kx, ky, idx, kc;

		bufsize=(int)sizeof(int16_t[4*4*4])*(image->iw+16);/*4 padded rows * 4 channels max * {pixels, errors, pastyrun, futureyrun}*/
		pixels=(int16_t*)malloc(bufsize);
		if(!pixels)
		{
			CRASH("Alloc error");
			return;
		}
		memset(pixels, 0, bufsize);
		for(ky=0, idx=0;ky<image->ih;++ky)
		{
			int16_t *rows[4];
			int sW[4]={0};

			rows[0]=pixels+((image->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*4;
			rows[1]=pixels+((image->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*4;
			rows[2]=pixels+((image->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*4;
			rows[3]=pixels+((image->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*4;
			for(kx=0;kx<image->iw;++kx, ++idx)
			{
				int16_t
					*NW	=rows[1]-1*4*4,
					*N	=rows[1]+0*4*4,
					*NE	=rows[1]+1*4*4,
					*NEE	=rows[1]+2*4*4,
					*NEEE	=rows[1]+3*4*4,
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
				for(kc=0;kc<4;++kc)
				{
					if(image->depth[kc])
					{
						int delta=((const int32_t*)HDR_DATC(image))[idx<<2|kc];
						int sym=delta<<1^delta>>31;
						int nbypass=31^LZCNT32((sW[kc]>>6)+1);
						bitsizes[kc]+=(int64_t)(sym>>nbypass)+nbypass+1;
						curr[kc+4]=sW[kc]=(2*sW[kc]+(sym<<6)+NEEE[kc+4])>>2;
					}
				}
				rows[0]+=4*4;
				rows[1]+=4*4;
				rows[2]+=4*4;
				rows[3]+=4*4;
			}
		}
		free(pixels);
		entropy[0]=image->src_depth[0]?bitsizes[0]/((double)image->iw*image->ih*image->src_depth[0]):0;
		entropy[1]=image->src_depth[1]?bitsizes[1]/((double)image->iw*image->ih*image->src_depth[1]):0;
		entropy[2]=image->src_depth[2]?bitsizes[2]/((double)image->iw*image->ih*image->src_depth[2]):0;
		entropy[3]=image->src_depth[3]?bitsizes[3]/((double)image->iw*image->ih*image->src_depth[3]):0;
	}
	else if(ec_method==ECTX_STATIC_O0)
	{
		int allocated, kc;

		allocated=0;
		if(!hist_full)
		{
			int maxdepth=calc_maxdepth(image, 0);
			hist_full=(int*)malloc(sizeof(int)<<maxdepth);
			allocated=1;
		}
		for(kc=0;kc<4;++kc)
		{
			if(image->depth[kc])
			{
				calc_histogram(((const int32_t*)HDR_DATC(image)), image->iw, image->ih, kc, 0, image->iw, 0, image->ih, image->depth[kc], hist_full, 0);
				entropy[kc]=calc_entropy(hist_full, 1<<image->depth[kc], image->iw*image->ih);
				if(image->src_depth[kc])
					entropy[kc]/=image->src_depth[kc];
			}
			else
				entropy[kc]=0;
		}
		if(allocated)
			free(hist_full);
	}
	else if(ec_method==ECTX_STATIC_O1)
	{
		double csizes[4]={0};
		int nlevels[4], half[4], mask[4], histsize, *hist, idx, ky, kx, kc, *currhist;

		nlevels[0]=1<<image->depth[0];
		nlevels[1]=1<<image->depth[1];
		nlevels[2]=1<<image->depth[2];
		nlevels[3]=1<<image->depth[3];
		half[0]=nlevels[0]>>1;
		half[1]=nlevels[1]>>1;
		half[2]=nlevels[2]>>1;
		half[3]=nlevels[3]>>1;
		mask[0]=nlevels[0]-1;
		mask[1]=nlevels[1]-1;
		mask[2]=nlevels[2]-1;
		mask[3]=nlevels[3]-1;
		histsize=(int)sizeof(int)*(
			+nlevels[0]*nlevels[0]
			+nlevels[1]*nlevels[1]
			+nlevels[2]*nlevels[2]
			+nlevels[3]*nlevels[3]
		);
		hist=(int*)malloc(histsize);
		if(!hist)
		{
			CRASH("Alloc error");
			return;
		}
		memset(hist, 0, histsize);
		for(ky=0, idx=0;ky<image->ih;++ky)
		{
			int W[4]={0};
			for(kx=0;kx<image->iw;++kx, idx+=4)
			{
				currhist=hist;
				for(kc=0;kc<image->nch;++kc)
				{
					int idx2=idx+kc;
					int target=((const int32_t*)HDR_DATC(image))[idx2];
					int ctx=(W[kc]+half[kc])&mask[kc];

					++currhist[ctx<<image->depth[kc]|((target+half[kc])&mask[kc])];
					currhist+=nlevels[kc]*nlevels[kc];
					W[kc]=target;
				}
			}
		}
		currhist=hist;
		for(kc=0;kc<4;++kc)
		{
			int totalsum, ks0, ks;
			double usize;

			if(!image->depth[kc])
			{
				currhist+=nlevels[kc]*nlevels[kc];
				continue;
			}
			totalsum=0;
			for(ks0=0;ks0<nlevels[kc];++ks0)
			{
				double e, norm;
				int sum;

				sum=0;
				for(ks=0;ks<nlevels[kc];++ks)
					sum+=currhist[ks];
				if(!sum)
				{
					currhist+=nlevels[kc];
					continue;
				}
				e=0, norm=1./sum;
				for(ks=0;ks<nlevels[kc];++ks)
				{
					int freq=currhist[ks];
					if(freq)
						e-=freq*log(freq*norm);
				}
				csizes[kc]+=e*(1./(8*M_LN2))+84;/*average GR overhead*/
				totalsum+=sum;
				currhist+=nlevels[kc];
			}
			if(totalsum!=image->iw*image->ih)
				CRASH("");
			usize=(double)image->iw*image->ih*image->src_depth[0];
			entropy[kc]=usize>1e-9?csizes[kc]*8/usize:0;/*entropy = cbitsize/ubitsize*/
		}
	}
	else if(ec_method==ECTX_DWT)
	{
		double e8[4][4]={0};
		int maxdepth, histsize, *hist, kc, kb, kx, ky, ks;

		maxdepth=image->depth[0];
		UPDATE_MAX(maxdepth, image->depth[1]);
		UPDATE_MAX(maxdepth, image->depth[2]);
		UPDATE_MAX(maxdepth, image->depth[3]);
		histsize=(int)sizeof(int)<<maxdepth;
		hist=(int*)malloc(histsize);
		if(!hist)
		{
			CRASH("Alloc error");
			return;
		}
		for(kc=0;kc<image->nch;++kc)
		{
			int nlevels=1<<image->depth[kc], half=nlevels>>1, mask=nlevels-1;
			for(kb=0;kb<4;++kb)
			{
				double e, gain;
				int x1=0, x2=0, y1=0, y2=0;
				switch(kb)
				{
				case 0:x1=image->iw>>1, x2=image->iw, y1=0, y2=image->ih>>1;break;
				case 1:x1=0, x2=image->iw>>1, y1=image->ih>>1, y2=image->ih;break;
				case 2:x1=image->iw>>1, x2=image->iw, y1=image->ih>>1, y2=image->ih;break;
				case 3:x1=0, x2=image->iw>>1, y1=0, y2=image->ih>>1;break;
				}
				memset(hist, 0, histsize);
				for(ky=y1;ky<y2;++ky)
				{
					for(kx=x1;kx<x2;++kx)
					{
						int val=((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|kc];
						val+=half;
						val&=mask;
						++hist[val];
					}
				}
				e=0, gain=1./((x2-x1)*(y2-y1));
				for(ks=0;ks<nlevels;++ks)
				{
					int freq=hist[ks];
					if(freq)
						e-=freq*log((double)freq*gain);
				}
				e8[kb][kc]=e*gain*(1./M_LN2);
			}
		}
		free(hist);
		entropy[0]=image->src_depth[0]?(e8[0][0]+e8[1][0]+e8[2][0]+e8[3][0])/(4*image->src_depth[0]):0;
		entropy[1]=image->src_depth[1]?(e8[0][1]+e8[1][1]+e8[2][1]+e8[3][1])/(4*image->src_depth[1]):0;
		entropy[2]=image->src_depth[2]?(e8[0][2]+e8[1][2]+e8[2][2]+e8[3][2])/(4*image->src_depth[2]):0;
		entropy[3]=image->src_depth[3]?(e8[0][3]+e8[1][3]+e8[2][3]+e8[3][3])/(4*image->src_depth[3]):0;
	}
	else if(ec_method==ECTX_INTERLEAVED)
	{
		double e8[4][4]={0};
		int maxdepth, histsize, *hist, kc, kb, ky, kx, ks;

		maxdepth=image->depth[0];
		UPDATE_MAX(maxdepth, image->depth[1]);
		UPDATE_MAX(maxdepth, image->depth[2]);
		UPDATE_MAX(maxdepth, image->depth[3]);
		histsize=(int)sizeof(int)<<maxdepth;
		hist=(int*)malloc(histsize);
		if(!hist)
		{
			CRASH("Alloc error");
			return;
		}
		for(kc=0;kc<image->nch;++kc)
		{
			int nlevels=1<<image->depth[kc], half=nlevels>>1, mask=nlevels-1;
			for(kb=0;kb<4;++kb)
			{
				double e, gain;
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
				for(ky=y0;ky<image->ih-1;ky+=2)
				{
					for(kx=x0;kx<image->iw-1;kx+=2)
					{
						int val=((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|kc];
						val+=half;
						val&=mask;
						++hist[val];
						++count;
					}
				}
				e=0, gain=1./count;
				for(ks=0;ks<nlevels;++ks)
				{
					int freq=hist[ks];
					if(freq)
						e-=freq*log((double)freq*gain);
				}
				e8[kb][kc]=e*gain*(1./M_LN2);
			}
		}
		free(hist);
		entropy[0]=image->src_depth[0]?(e8[0][0]+e8[1][0]+e8[2][0]+e8[3][0])/(4*image->src_depth[0]):0;
		entropy[1]=image->src_depth[1]?(e8[0][1]+e8[1][1]+e8[2][1]+e8[3][1])/(4*image->src_depth[1]):0;
		entropy[2]=image->src_depth[2]?(e8[0][2]+e8[1][2]+e8[2][2]+e8[3][2])/(4*image->src_depth[2]):0;
		entropy[3]=image->src_depth[3]?(e8[0][3]+e8[1][3]+e8[2][3]+e8[3][3])/(4*image->src_depth[3]):0;
	}
	else if(ec_method==ECTX_ABAC0||ec_method==ECTX_ABAC1)
		calc_csize_abac(image, ec_method==ECTX_ABAC1, entropy);
	else if(ec_method==ECTX_YUV422||ec_method==ECTX_YUV420)
	{
		int nlevels[4]={0}, hoffsets[11]={0}, hsum[_countof(hoffsets)-1]={0}, ky, kx, kc, hsize, *hist, *hist_y, *hist_u, *hist_v;

		nlevels[0]=1<<image->depth[0];
		nlevels[1]=1<<image->depth[1];
		nlevels[2]=1<<image->depth[2];
		hoffsets[0x0]=0;
		hoffsets[0x1]=nlevels[0];
		hoffsets[0x2]=nlevels[0];
		hoffsets[0x3]=nlevels[1];
		hoffsets[0x4]=nlevels[1];
		hoffsets[0x5]=nlevels[1];
		hoffsets[0x6]=nlevels[1];
		hoffsets[0x7]=nlevels[2];
		hoffsets[0x8]=nlevels[2];
		hoffsets[0x9]=nlevels[2];
		hoffsets[0xA]=nlevels[2];
		for(kc=1;kc<(int)_countof(hoffsets);++kc)
			hoffsets[kc]+=hoffsets[kc-1];
		hsize=hoffsets[_countof(hoffsets)-1]*(int)sizeof(int);
		hist=(int*)malloc(hsize);
		if(!hist)
		{
			CRASH("Alloc error");
			return;
		}
		memset(hist, 0, hsize);
		hist_y=hist;
		hist_u=hist+hoffsets[2];
		hist_v=hist+hoffsets[6];
		for(ky=0;ky<image->ih;++ky)
		{
			for(kx=0;kx<image->iw;++kx)
			{
				for(kc=0;kc<3;++kc)
				{
					int val=((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|kc]+(nlevels[kc]>>1);

					val&=nlevels[kc]-1;
					if(ec_method==ECTX_YUV420)
					{
						switch(kc)
						{
						case 0:
							if(!(kx&7)&&!(ky&7))/*DC*/
								++hist_y[0<<image->depth[0]|val], ++hsum[0];
							else/*AC*/
								++hist_y[1<<image->depth[0]|val], ++hsum[1];
							break;
						case 1:
							if(!(kx>>1&7)&&!(ky>>1&7))/*DC*/
							{
								if(!(kx&1)&&!(ky&1))/*base*/
									++hist_u[0<<image->depth[1]|val], ++hsum[2];
								else/*interpolated*/
									++hist_u[1<<image->depth[1]|val], ++hsum[3];
							}
							else/*AC*/
							{
								if(!(kx&1)&&!(ky&1))/*base*/
									++hist_u[2<<image->depth[1]|val], ++hsum[4];
								else/*interpolated*/
									++hist_u[3<<image->depth[1]|val], ++hsum[5];
							}
							break;
						case 2:
							if(!(kx>>1&7)&&!(ky>>1&7))/*DC*/
							{
								if(!(kx&1)&&!(ky&1))/*base*/
									++hist_v[0<<image->depth[2]|val], ++hsum[6];
								else/*interpolated*/
									++hist_v[1<<image->depth[2]|val], ++hsum[7];
							}
							else/*AC*/
							{
								if(!(kx&1)&&!(ky&1))/*base*/
									++hist_v[2<<image->depth[2]|val], ++hsum[8];
								else/*interpolated*/
									++hist_v[3<<image->depth[2]|val], ++hsum[9];
							}
							break;
						}
					}
					else/*YUV422*/
					{
						switch(kc)
						{
						case 0:
							if(!(kx&7)&&!(ky&7))/*DC*/
								++hist_y[0<<image->depth[0]|val], ++hsum[0];
							else/*AC*/
								++hist_y[1<<image->depth[0]|val], ++hsum[1];
							break;
						case 1:
							if(!(kx>>1&7)&&!(ky&7))/*DC*/
							{
								if(!(kx&1))/*base*/
									++hist_u[0<<image->depth[1]|val], ++hsum[2];
								else/*interpolated*/
									++hist_u[1<<image->depth[1]|val], ++hsum[3];
							}
							else/*AC*/
							{
								if(!(kx&1))/*base*/
									++hist_u[2<<image->depth[1]|val], ++hsum[4];
								else/*interpolated*/
									++hist_u[3<<image->depth[1]|val], ++hsum[5];
							}
							break;
						case 2:
							if(!(kx>>1&7)&&!(ky&7))/*DC*/
							{
								if(!(kx&1))/*base*/
									++hist_v[0<<image->depth[2]|val], ++hsum[6];
								else/*interpolated*/
									++hist_v[1<<image->depth[2]|val], ++hsum[7];
							}
							else/*AC*/
							{
								if(!(kx&1))/*base*/
									++hist_v[2<<image->depth[2]|val], ++hsum[8];
								else/*interpolated*/
									++hist_v[3<<image->depth[2]|val], ++hsum[9];
							}
							break;
						}
					}
				}
			}
		}
		{
			int eidxs[]=
			{
				0, 0,
				1, 1, 1, 1,
				2, 2, 2, 2,
			};
			int res, ks;

			memset(entropy, 0, sizeof(double[4]));
			for(kc=0;kc<(int)_countof(hoffsets)-1;++kc)
			{
				int eidx=eidxs[kc], nlevels2=hoffsets[kc+1]-hoffsets[kc], sum=hsum[kc];
				for(ks=0;ks<nlevels2;++ks)
				{
					int freq=hist[hoffsets[kc]+ks];
					if(freq)
						entropy[eidx]-=freq*log((double)freq/sum);
				}
			}
			entropy[0]*=1./M_LN2;
			entropy[1]*=1./M_LN2;
			entropy[2]*=1./M_LN2;
			res=image->iw*image->ih;
			entropy[0]=image->src_depth[0]?entropy[0]/(res*image->src_depth[0]):0;
			entropy[1]=image->src_depth[0]?entropy[1]/(res*image->src_depth[1]):0;
			entropy[2]=image->src_depth[0]?entropy[2]/(res*image->src_depth[2]):0;
			free(hist);
		}
	}
	else
		calc_csize_ec(image, ec_method, ec_adaptive?ec_adaptive_threshold:0, ec_expbits, ec_msb, ec_lsb, entropy);
}

static ptrdiff_t calc_psnr(const Image *ima, const Image *imb, double *ret_rmse, double *ret_psnr)
{
	double rmse[5]={0}, psnr[5]={0};
	intptr_t res, idx, k, kc, usize;
	double vmax[4]={0};

	if(ima->iw!=imb->iw||ima->ih!=imb->ih||ima->nch!=imb->nch)
	{
		if(ret_rmse)
		{
			ret_rmse[0]=INFINITY;
			ret_rmse[1]=INFINITY;
			ret_rmse[2]=INFINITY;
			ret_rmse[3]=INFINITY;
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
	res=(intptr_t)ima->iw*ima->ih;
	idx=-1;
	vmax[0]=1./((1<<ima->src_depth[0])-1);
	vmax[1]=1./((1<<ima->src_depth[1])-1);
	vmax[2]=1./((1<<ima->src_depth[2])-1);
	vmax[3]=1./((1<<ima->src_depth[3])-1);
	for(k=0;k<res;++k)
	{
		int
			dr=((const int32_t*)HDR_DATC(imb))[k<<2|0]-((const int32_t*)HDR_DATC(ima))[k<<2|0],
			dg=((const int32_t*)HDR_DATC(imb))[k<<2|1]-((const int32_t*)HDR_DATC(ima))[k<<2|1],
			db=((const int32_t*)HDR_DATC(imb))[k<<2|2]-((const int32_t*)HDR_DATC(ima))[k<<2|2],
			da=((const int32_t*)HDR_DATC(imb))[k<<2|3]-((const int32_t*)HDR_DATC(ima))[k<<2|3];
		if(idx<0&&(dr|dg|db|da))
			idx=k;
		rmse[0]+=dr*dr;
		rmse[1]+=dg*dg;
		rmse[2]+=db*db;
		rmse[3]+=da*da;
	}
	usize=0;
	for(kc=0;kc<4;++kc)
	{
		if(ima->src_depth[kc])
			usize+=res;
	}
	rmse[4]=sqrt((rmse[0]+rmse[1]+rmse[2]+rmse[3])/usize);
	psnr[4]=-20*log10(rmse[4]*vmax[0]);
	for(kc=0;kc<4;++kc)
	{
		rmse[kc]=sqrt(rmse[kc]/res);
		psnr[kc]=-20*log10(rmse[kc]*vmax[kc]);
		if(ima->src_depth[kc])
			usize+=res;
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
	double rmse[5], psnr[5], sqe[4], ssim[4];
} LossyCtx;
static uint32_t __stdcall sample_thread_lossy(void *param)
{
	double entropy[4]={0};
	LossyCtx *ctx;
	Image *im2=0;
	int kc;

	ctx=(LossyCtx*)param;
	ctx->usize=image_getBMPsize(ctx->image);
	image_copy(&im2, ctx->image);/*save original*/
	if(!im2)
	{
		CRASH("Alloc error");
		return 0;
	}

	apply_selected_transforms(&ctx->image, 0, 1, 0);/*fwd*/

	calc_csize_stateful(ctx->image, 0, entropy);
	for(kc=0;kc<4;++kc)
	{
		int depth=ctx->image->src_depth[kc];
		double invCR=depth?entropy[kc]:0;
		ctx->csize[kc]=invCR*ctx->image->iw*ctx->image->ih*ctx->image->src_depth[kc]/8;
	}

	apply_selected_transforms(&ctx->image, 0, 0, 1);/*inv*/

	calc_psnr(im2, ctx->image, ctx->rmse, ctx->psnr);
	measure_ssim_avg(im2, ctx->image, ctx->ssim);
	{
		intptr_t res=(intptr_t)ctx->image->iw*ctx->image->ih;
		ctx->sqe[0]=res*ctx->rmse[0]*ctx->rmse[0];
		ctx->sqe[1]=res*ctx->rmse[1]*ctx->rmse[1];
		ctx->sqe[2]=res*ctx->rmse[2]*ctx->rmse[2];
		ctx->sqe[3]=res*ctx->rmse[3]*ctx->rmse[3];
	}
	free(im2);
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
	int k;

	loud_transforms=0;
	path=dialog_open_folder();
	if(!path)
		return;
	filenames=get_filenames((char*)HDR_DATA(path), ext, _countof(ext), 1);
	if(!filenames)
	{
		array_free(&path);
		return;
	}

	DisableProcessWindowsGhosting();
	console_start();
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H:%M:%S");
	console_log("Lossy Batch Test  %s  %s\n", g_buf, (char*)HDR_DATA(path));
	array_free(&path);
	console_log("Enter number of threads: ");
	nthreads=console_scan_int();
	t=time_sec();
	total_usize=0;
	maxlen=0;
	for(k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ARRAY_ALLOC(LossyCtx, q, 0, 0, nthreads, 0);
	for(k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		Image *image=image_load((char*)HDR_DATA(fn2[0]), (int)fn2[0]->count);
		LossyCtx ctx={0};
		int k2;

		ctx.image=image;
		ctx.idx=k;
		if(!image)
			continue;
		ARRAY_APPEND(q, &ctx, 1, 1, 0);
		if((int)q->count>=nthreads||k+1>=(int)filenames->count)
		{
			HANDLE *handles=(HANDLE*)malloc(q->count*sizeof(HANDLE));
			if(!handles)
			{
				CRASH("Alloc error");
				return;
			}
			memset(handles, 0, q->count*sizeof(HANDLE));
			for(k2=0;k2<(int)q->count;++k2)
			{
				LossyCtx *ptr=(LossyCtx*)array_at(&q, k2);
				handles[k2]=(void*)_beginthreadex(0, 0, sample_thread_lossy, ptr, 0, 0);
				if(!handles[k2])
				{
					CRASH("Alloc error");
					return;
				}
			}
			WaitForMultipleObjects((int)q->count, handles, TRUE, INFINITE);
			for(k2=0;k2<(int)q->count;++k2)
			{
				LossyCtx *ptr=(LossyCtx*)array_at(&q, k2);
				double csize=ptr->csize[0]+ptr->csize[1]+ptr->csize[2]+ptr->csize[3];
				fn2=(ArrayHandle*)array_at(&filenames, ptr->idx);
				console_log(
					"%5d/%5d %s%*sUTYUV %13.2f %13.2f %13.2f %13.2f %13.2f  BPD %8.4f  PSNR %10.4f %10.4f %10.4f %10.4f  SSIM%7.4f\n"
					, (int)(k+1-(int)q->count+k2+1)
					, (int)filenames->count
					, (char*)HDR_DATA(fn2[0])
					, (int)(maxlen-fn2[0]->count+1), ""
					, ptr->usize
					, csize
					, ptr->csize[0]
					, ptr->csize[1]
					, ptr->csize[2]
					, 8.*csize/ptr->usize
					, ptr->psnr[0]
					, ptr->psnr[1]
					, ptr->psnr[2]
					, ptr->psnr[3]
					, ptr->ssim[3]
				);
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
		double ctotal, CR, rmse, psnr, ssim_norm;
		int nprinted;

		ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
		CR=total_usize/ctotal;
		rmse=sqrt((total_sqe[0]+total_sqe[1]+total_sqe[2]+total_sqe[3])/total_usize);
		psnr=20*log10(255/rmse);
		ssim_norm=1024*1024/(double)total_usize;
		total_ssim[0]*=ssim_norm;
		total_ssim[1]*=ssim_norm;
		total_ssim[2]*=ssim_norm;
		total_ssim[3]*=ssim_norm;
		t=time_sec()-t;
		timedelta2str(g_buf, G_BUF_SIZE, t);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y%m%d_%H%M%S");
		nprinted=snprintf(str, sizeof(str)-1,
			"%13.2f %13.2f %13.2f %13.2f %13.2f  BPD %8.4f  RMSE %10.4f PSNR %10.4f  %10.4f sec  SSIM%7.4f  dist%3d  %s"
			, total_usize
			, ctotal
			, total_csize[0]
			, total_csize[1]
			, total_csize[2]
			, 8./CR
			, rmse
			, psnr
			, t
			, total_ssim[3]
			, g_dist
			, g_buf
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

static int customtransforms_getflag(uint8_t tid)
{
	return
		tid==CT_FWD_CUSTOM||
		tid==CT_INV_CUSTOM||
		tid==ST_FWD_CUSTOM||
		tid==ST_INV_CUSTOM||
		tid==ST_CONVTEST||
		tid==ST_CONVTEST2||
		tid==ST_FWD_CC||
		tid==ST_INV_CC
	;
}
static void transforms_update(void)
{
	if(transforms)
	{
		int k;

		memset(transforms_mask, 0, T_COUNT);
		transforms_customenabled=0;
		for(k=0;k<(int)transforms->count;++k)/*update trackers*/
		{
			uint8_t tid2=HDR_DATA(transforms)[k];
			if(tid2<T_COUNT)
				transforms_mask[tid2]|=1;
			transforms_customenabled|=customtransforms_getflag(tid2);
		}
	}
}
static void transforms_removebyid(uint32_t tid)
{
	if(transforms)
	{
		int k;

		for(k=(int)transforms->count-1;k>=0;--k)
		{
			if(tid==HDR_DATA(transforms)[k])
			{
				array_erase(&transforms, k, 1);
				break;/*remove only one*/
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
static void transforms_append(uint32_t tid)
{
	if(tid<T_COUNT&&tid!=CST_INV_SEPARATOR)
	{
		if(!transforms)
		{
			ARRAY_ALLOC(char, transforms, &tid, 1, 0, 0);
			transforms_mask[tid]|=1;
			transforms_customenabled|=customtransforms_getflag((uint8_t)tid);
		}
		else
		{
			int idx, k;

			idx=-1;/*first idx of a transform of this type*/
			if(GET_KEY_STATE(KEY_CTRL))/*replace all transforms of this type*/
			{
				for(k=0;k<(int)transforms->count;)
				{
					if((tid<CST_COMPARE)==(HDR_DATA(transforms)[k]<CST_COMPARE))
					{
						if(idx==-1)
							idx=k;
						array_erase(&transforms, k, 1);
						transforms_mask[HDR_DATA(transforms)[k]]=0;
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
static void transforms_printname(float x, float y, uint32_t tid, int place, int64_t highlight)
{
	const char *a=0;
	switch(tid)
	{
#if 0
	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
#endif
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
#if 0
	case CT_FWD_CrCgCb:		a="C  Fwd CrCgCb";		break;
	case CT_INV_CrCgCb:		a="C  Inv CrCgCb";		break;
#endif
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
#if 0
	case CT_FWD_J2K2:		a="C  Fwd J2K/2";		break;
	case CT_INV_J2K2:		a="C  Inv J2K/2";		break;
	case CT_FWD_YRGB_v1:		a="C  Fwd YRGB v1";		break;
	case CT_INV_YRGB_v1:		a="C  Inv YRGB v1";		break;
	case CT_FWD_YRGB_v2:		a="C  Fwd YRGB v2";		break;
	case CT_INV_YRGB_v2:		a="C  Inv YRGB v2";		break;
	case CT_FWD_CMYK:		a="C  Fwd CMYK";		break;
	case CT_INV_CMYK_DUMMY:		a="";				break;/*no inverse*/
#endif
	case CT_FWD_YCbCr:		a="C  Fwd YCbCr";		break;
	case CT_INV_YCbCr:		a="C  Inv YCbCr";		break;
	case CT_FWD_XYB:		a="C  Fwd XYB";			break;
	case CT_INV_XYB:		a="C  Inv XYB";			break;
#if 0
	case CT_FWD_MATRIX:		a="C  Fwd Matrix";		break;
	case CT_INV_MATRIX:		a="C  Inv Matrix";		break;
	case CT_FWD_XGZ:		a="C  Fwd XGZ";			break;
	case CT_INV_XGZ:		a="C  Inv XGZ";			break;
	case CT_FWD_XYZ:		a="C  Fwd XYZ";			break;
	case CT_INV_XYZ:		a="C  Inv XYZ";			break;
	case CT_FWD_EXP:		a="C  Fwd Experimental";	break;
	case CT_INV_EXP:		a="C  Inv Experimental";	break;
	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
#endif
	case CT_FWD_CUSTOM:		a="C  Fwd CUSTOM";		break;
	case CT_INV_CUSTOM:		a="C  Inv CUSTOM";		break;
#if 0
	case CT_FWD_PCA:		a="C  Fwd PCA";			break;
	case CT_INV_PCA:		a="C  Inv PCA";			break;
#endif
	case CT_FWD_HISTMATCH:		a="C  Fwd HistMatch";		break;
	case CT_INV_HISTMATCH:		a="C  Inv HistMatch";		break;
#if 0
	case CT_FWD_QUAD:		a="C  Fwd Quad";		break;
	case CT_INV_QUAD:		a="C  Inv Quad";		break;
#endif
	case CST_COMPARE:		a="Compare";			break;
	case CST_INV_SEPARATOR:		a="";				break;
#if 0
	case ST_PREPROC_GRAD:		a=" S Preproc Grad";		break;
	case ST_PREPROC_X:		a=" S Preproc X";		break;
	case ST_PREPROC_X2:		a=" S Preproc X2";		break;
#endif
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
#if 0
	case ST_FWD_LPC:		a=" S Fwd LPC";			break;
	case ST_INV_LPC:		a=" S Inv LPC";			break;
#endif
	case ST_FWD_SIGNBIT:		a=" S Fwd Sign Bit";		break;
	case ST_FWD_MAGNITUDE:		a=" S Fwd Magnitude";		break;
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
	case ST_FWD_EXTRAP:		a=" S Fwd Extrapolate";		break;
	case ST_FWD_L1DITHER:		a=" S Fwd L1 Dither";		break;
	case ST_FWD_QUANT:		a=" S Fwd Quantize";		break;
	case ST_INV_QUANT:		a=" S Inv Quantize";		break;
	case ST_FWD_ADAQUANT:		a=" S Fwd AdaQuant";		break;
	case ST_INV_ADAQUANT:		a=" S Inv AdaQuant";		break;
	case ST_FWD_AV4:		a=" S Fwd AV4";			break;
	case ST_INV_AV4:		a=" S Inv AV4";			break;
	case ST_FWD_SEL4:		a=" S Fwd Sel4";		break;
	case ST_INV_SEL4:		a=" S Inv Sel4";		break;
	case ST_FWD_SELECT:		a=" S Fwd Select";		break;
	case ST_INV_SELECT:		a=" S Inv Select";		break;
	case ST_FWD_PAETH:		a=" S Fwd Paeth";		break;
	case ST_INV_PAETH:		a=" S Inv Paeth";		break;
	case ST_FWD_SELECT4:		a=" S Fwd SelGrad cRCT";	break;
	case ST_INV_SELECT4:		a=" S Inv SelGrad cRCT";	break;
	case ST_FWD_CGPLUS:		a=" S Fwd CGplus";		break;
	case ST_INV_CGPLUS:		a=" S Inv CGplus";		break;
	case ST_FWD_CG422:		a=" S Fwd CG422";		break;
	case ST_INV_CG422:		a=" S Inv CG422";		break;
	case ST_FWD_CG420:		a=" S Fwd CG420";		break;
	case ST_INV_CG420:		a=" S Inv CG420";		break;
	case ST_FWD_AV2:		a=" S Fwd (N+W)/2";		break;
	case ST_INV_AV2:		a=" S Inv (N+W)/2";		break;
#if 0
	case ST_FWD_MIX2:		a=" S Fwd MIX2";		break;
	case ST_INV_MIX2:		a=" S Inv MIX2";		break;
#endif
	case ST_FWD_GRAY:		a=" S Fwd Gray";		break;
	case ST_INV_GRAY:		a=" S Inv Gray";		break;
	case ST_FWD_MIXR:		a=" S Fwd Mix RCT";		break;
	case ST_INV_MIXR:		a=" S Inv Mix RCT";		break;
	case ST_FWD_MIXN:		a=" S Fwd MIX N";		break;
	case ST_INV_MIXN:		a=" S Inv MIX N";		break;
#if 0
	case ST_FWD_MIXNC:		a=" S Fwd MIX NC";		break;
	case ST_INV_MIXNC:		a=" S Inv MIX NC";		break;
	case ST_FWD_BESTN:		a=" S Fwd Best N";		break;
	case ST_INV_BESTN:		a=" S Inv Best N";		break;
	case ST_FWD_RLS:		a=" S Fwd RLS";			break;
	case ST_INV_RLS:		a=" S Inv RLS";			break;
	case ST_FWD_AV3:		a=" S Fwd AV3";			break;
	case ST_INV_AV3:		a=" S Inv AV3";			break;
#endif
	case ST_FWD_WGRAD:		a="CS Fwd WGrad";		break;
	case ST_INV_WGRAD:		a="CS Inv WGrad";		break;
#if 0
	case ST_FWD_WGRAD2:		a="CS Fwd WGrad2";		break;
	case ST_INV_WGRAD2:		a="CS Inv WGrad2";		break;
#endif
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
#if 0
	case ST_FWD_WMIX:		a=" S Fwd WMIX";		break;
	case ST_INV_WMIX:		a=" S Inv WMIX";		break;
#endif
	case ST_FWD_TABLE:		a=" S Fwd Table";		break;
	case ST_INV_TABLE:		a=" S Inv Table";		break;
	case ST_FWD_AWAV:		a=" S Fwd AWAV";		break;
	case ST_INV_AWAV:		a=" S Inv AWAV";		break;
	case ST_FWD_LWAV:		a=" S Fwd LWAV";		break;
	case ST_INV_LWAV:		a=" S Inv LWAV";		break;
#if 0
	case ST_FWD_ECOEFF:		a=" S Fwd E-Coeff";		break;
	case ST_INV_ECOEFF:		a=" S Inv E-Coeff";		break;
	case ST_FWD_AVERAGE:		a=" S Fwd Average";		break;
	case ST_INV_AVERAGE:		a=" S Inv Average";		break;
	case ST_FWD_MULTISTAGE:		a=" S Fwd Multistage";		break;
	case ST_INV_MULTISTAGE:		a=" S Inv Multistage";		break;
#endif
	case ST_FWD_ZIPPER:		a=" S Fwd Zipper";		break;
	case ST_INV_ZIPPER:		a=" S Inv Zipper";		break;
#if 0
	case ST_FWD_DIR:		a=" S Fwd Dir";			break;
	case ST_INV_DIR:		a=" S Inv Dir";			break;
#endif
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
#if 0
	case ST_FWD_WP2:		a=" S Fwd WP2";			break;
	case ST_INV_WP2:		a=" S Inv WP2";			break;
#endif
	case ST_FWD_WPU:		a="CS Fwd WPU";			break;
	case ST_INV_WPU:		a="CS Inv WPU";			break;
#if 0
	case ST_FWD_DEFERRED:		a=" S Fwd DEFERRED";		break;
	case ST_INV_DEFERRED:		a=" S Inv DEFERRED";		break;
#endif
	case ST_FWD_MM:			a=" S Fwd MM";			break;
	case ST_INV_MM:			a=" S Inv MM";			break;
	case ST_FWD_CUSTOM:		a=" S Fwd CUSTOM";		break;
	case ST_INV_CUSTOM:		a=" S Inv CUSTOM";		break;
	case ST_FWD_CC:			a=" S Fwd CC";			break;
	case ST_INV_CC:			a=" S Inv CC";			break;
#if 0
	case ST_FWD_CUSTOM2:		a=" S Fwd CUSTOM2";		break;
	case ST_INV_CUSTOM2:		a=" S Inv CUSTOM2";		break;
	case ST_FWD_CUSTOM4:		a=" S Fwd CUSTOM4";		break;
	case ST_INV_CUSTOM4:		a=" S Inv CUSTOM4";		break;
	case ST_FWD_KALMAN:		a=" S Fwd Kalman";		break;
	case ST_INV_KALMAN:		a=" S Inv Kalman";		break;
	case ST_FWD_LOGIC:		a=" S Fwd Logic";		break;
	case ST_INV_LOGIC:		a=" S Inv Logic";		break;
	case ST_FWD_LEARNED:		a=" S Fwd Learned";		break;
	case ST_INV_LEARNED:		a=" S Inv Learned";		break;
#ifdef ALLOW_OPENCL
	case ST_FWD_LEARNED_GPU:	a=" S Fwd Learned GPU";		break;
	case ST_INV_LEARNED_GPU:	a=" S Inv Learned GPU";		break;
#endif
	case ST_FWD_CFL:		a=" S Fwd CfL";			break;
	case ST_INV_CFL:		a=" S Inv CfL";			break;
	case ST_FWD_JOINT:		a=" S Fwd Joint";		break;
	case ST_INV_JOINT:		a=" S Inv Joint";		break;
	case ST_FWD_HYBRID3:		a=" S Fwd Hybrid3";		break;
	case ST_INV_HYBRID3:		a=" S Inv Hybrid3";		break;

	case ST_FWD_DIFF2D:		a=" S Fwd 2D derivative";	break;
	case ST_INV_DIFF2D:		a=" S Inv 2D derivative";	break;
	case ST_FWD_HPF:		a=" S Fwd HPF";			break;
	case ST_INV_HPF:		a=" S Inv HPF";			break;
	case ST_FWD_GRAD2:		a=" S Fwd Grad2";		break;
	case ST_INV_GRAD2:		a=" S Inv Grad2";		break;
	case ST_FWD_ADAPTIVE:		a=" S Fwd Adaptive";		break;
	case ST_INV_ADAPTIVE:		a=" S Inv Adaptive";		break;
	case ST_FWD_JMJ:		a=" S Fwd JMJ";			break;
	case ST_INV_JMJ:		a=" S Inv JMJ";			break;
	case ST_FWD_SORTNB:		a=" S Fwd Sort Nb";		break;
	case ST_INV_SORTNB:		a=" S Inv Sort Nb";		break;
	case ST_FWD_MEDIAN:		a=" S Fwd Median";		break;
	case ST_INV_MEDIAN:		a=" S Inv Median";		break;
	case ST_FWD_DCT3PRED:		a=" S Fwd DCT3 Predictor";	break;
	case ST_INV_DCT3PRED:		a=" S Inv DCT3 Predictor";	break;
	case ST_FWD_PATHPRED:		a=" S Fwd Path Predictor";	break;
	case ST_INV_PATHPRED:		a=" S Inv Path Predictor";	break;
	case ST_FWD_GRADPRED:		a=" S Fwd Gradient";		break;
	case ST_INV_GRADPRED:		a=" S Inv Gradient";		break;
	case ST_FWD_GRAD2:		a=" S Fwd Grad2";		break;
	case ST_INV_GRAD2:		a=" S Inv Grad2";		break;
	case ST_FWD_CTX:		a=" S Fwd Ctx";			break;
	case ST_INV_CTX:		a=" S Inv Ctx";			break;
	case ST_FWD_C03:		a=" S Fwd C03";			break;
	case ST_INV_C03:		a=" S Inv C03";			break;
	case ST_FWD_C10:		a=" S Fwd C10";			break;
	case ST_INV_C10:		a=" S Inv C10";			break;
	case ST_FWD_C20:		a=" S Fwd C20";			break;
	case ST_INV_C20:		a=" S Inv C20";			break;
	case ST_FWD_WU97:		a=" S Fwd Wu 97";		break;
	case ST_INV_WU97:		a=" S Inv Wu 97";		break;
	case ST_FWD_BITWISE:		a=" S Fwd Bitwise";		break;
	case ST_INV_BITWISE:		a=" S Inv Bitwise";		break;
	case ST_FWD_DCT4:		a=" S Fwd DCT4";		break;
	case ST_INV_DCT4:		a=" S Inv DCT4";		break;
	case ST_FWD_DCT8:		a=" S Fwd DCT8";		break;
	case ST_INV_DCT8:		a=" S Inv DCT8";		break;
	case ST_FWD_SHUFFLE:		a=" S Fwd Shuffle";		break;
	case ST_INV_SHUFFLE:		a=" S Inv Shuffle";		break;
	case ST_FWD_SPLIT:		a=" S Fwd Split";		break;
	case ST_INV_SPLIT:		a=" S Inv Split";		break;
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
#if 0
	case ST_FWD_GRAD_DWT:		a=" S Fwd Gradient DWT";	break;
	case ST_INV_GRAD_DWT:		a=" S Inv Gradient DWT";	break;
	case ST_FWD_DEC_DWT:		a=" S Fwd Dec. DWT";		break;
	case ST_INV_DEC_DWT:		a=" S Inv Dec. DWT";		break;
	case ST_FWD_EXPDWT:		a=" S Fwd Exp DWT";		break;
	case ST_INV_EXPDWT:		a=" S Inv Exp DWT";		break;
	case ST_FWD_CUSTOM_DWT:		a=" S Fwd CUSTOM DWT";		break;
	case ST_INV_CUSTOM_DWT:		a=" S Inv CUSTOM DWT";		break;
#endif
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
		int64_t c0=0;
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

static int send_image_separate_subpixels(Image const *image, uint32_t *txid_r, uint32_t *txid_g, uint32_t *txid_b, uint32_t *txid_a)
{
	intptr_t res, k;
	int shift[4]={0};
	uint8_t *temp_r, *temp_g, *temp_b, *temp_a;

	res=(intptr_t)image->iw*image->ih;
	shift[0]=MAXVAR(0, image->depth[0]-8);
	shift[1]=MAXVAR(0, image->depth[1]-8);
	shift[2]=MAXVAR(0, image->depth[2]-8);
	shift[3]=MAXVAR(0, image->depth[3]-8);
	temp_r=(uint8_t*)malloc(res*sizeof(uint8_t[4]));
	temp_g=(uint8_t*)malloc(res*sizeof(uint8_t[4]));
	temp_b=(uint8_t*)malloc(res*sizeof(uint8_t[4]));
	temp_a=(uint8_t*)malloc(res*sizeof(uint8_t[4]));
	if(!temp_r||!temp_g||!temp_b||!temp_a)
	{
		CRASH("Alloc error");
		return 0;
	}
	for(k=0;k<res;++k)
	{
		int
			r=(uint8_t)((((const int32_t*)HDR_DATC(image))[k<<2|0]>>shift[0]))^128,
			g=(uint8_t)((((const int32_t*)HDR_DATC(image))[k<<2|1]>>shift[1]))^128,
			b=(uint8_t)((((const int32_t*)HDR_DATC(image))[k<<2|2]>>shift[2]))^128,
			a=(uint8_t)((((const int32_t*)HDR_DATC(image))[k<<2|3]>>shift[3]))^128;
		if(separate_grayscale)
		{
			((int*)temp_r)[k]=0xFF000000|r<<16|r<<8|r;
			((int*)temp_g)[k]=0xFF000000|g<<16|g<<8|g;
			((int*)temp_b)[k]=0xFF000000|b<<16|b<<8|b;
		}
		else
		{
			((int*)temp_r)[k]=0xFF000000|r<< 0;
			((int*)temp_g)[k]=0xFF000000|g<< 8;
			((int*)temp_b)[k]=0xFF000000|b<<16;
		}
		((int*)temp_a)[k]=0xFF000000|a<<16|a<<8|a;
	}
	if(!*txid_r)
	{
		uint32_t txids[4];

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
static int hist[768], histmax[3];
static double hist_mean[3]={0}, hist_sdev[3]={0}, hist_mad[3]={0}, hist_gausslike[3]={0}, hist_laplacelike[3]={0};/*maintained by chart_hist_update*/
static int *hist_full=0, hist_full_size=0;
static void chart_hist_update(Image const *image, int x1, int x2, int y1, int y2, int *_hist, int *_histmax, float *CR)
{
	double entropy;
	int xcount, ycount, count;

	x1=MAXVAR(x1, 0);
	x2=MINVAR(x2, image->iw);
	y1=MAXVAR(y1, 0);
	y2=MINVAR(y2, image->ih);
	xcount=x2-x1;
	ycount=y2-y1;
	count=xcount*ycount;
	if(count)
	{
		double norm;
		int kc;

		norm=1./count;
		for(kc=0;kc<3;++kc)
		{
			double invsdev, invb, gnorm, lnorm;
			int *currhist, ks;

			currhist=_hist+((ptrdiff_t)kc<<8);
			calc_histogram((const int32_t*)HDR_DATC(image), image->iw, image->ih, kc, x1, x2, y1, y2, image->depth[kc], hist_full, currhist);
			hist_mean[kc]=0;
			for(ks=0;ks<256;++ks)
				hist_mean[kc]+=(double)ks*currhist[ks];
			hist_mean[kc]*=norm;
			hist_sdev[kc]=0;
			hist_mad[kc]=0;
			for(ks=0;ks<256;++ks)
			{
				int freq=currhist[ks];
				double val=ks-hist_mean[kc];
				hist_sdev[kc]+=val*val*freq;
				hist_mad[kc]+=fabs(val)*freq;
			}
			hist_sdev[kc]=sqrt(hist_sdev[kc]/(count-1));
			hist_mad[kc]/=count-1;
			hist_gausslike[kc]=0;
			hist_laplacelike[kc]=0;
			invsdev=1/hist_sdev[kc];
			invb=1/hist_mad[kc];
			gnorm=1/(sqrt(2*M_PI)*hist_sdev[kc]);
			lnorm=1/(2*hist_mad[kc]);
			for(ks=0;ks<256;++ks)
			{
				int freq=currhist[ks];
				if(freq)
				{
					double pgauss, plaplace;

					pgauss=(ks-hist_mean[kc])*invsdev;
					pgauss=gnorm*exp(-pgauss*pgauss);
					plaplace=(ks-hist_mean[kc])*invb;
					plaplace=lnorm*exp(-fabs(plaplace));

					hist_gausslike[kc]+=freq*(pgauss>1e-9?-log(pgauss):256);/*finite codelength*/
					hist_laplacelike[kc]+=freq*(plaplace>1e-9?-log(plaplace):256);
				}
			}
			hist_gausslike[kc]*=1./(8*M_LN2);
			hist_laplacelike[kc]*=1./(8*M_LN2);
			for(ks=0;ks<256;++ks)
			{
				if(_histmax[kc]<currhist[ks])
					_histmax[kc]=currhist[ks];
			}
			if(CR)
			{
				entropy=calc_entropy(hist_full, 1<<image->depth[kc], count);
				CR[kc]=(float)(1/entropy);
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
static void jhc_getboxbounds(int xs, int ys, int dx, int dy, int iw, int ih, int *bounds)/*bounds: {x1, x2, y1, y2}*/
{
	int ix=(int)(((int64_t)xs*iw+(wndw>>1))/wndw);
	int iy=(int)(((int64_t)ys*ih+(wndh>>1))/wndh);
	move_box_in_window(ix, dx, 0, iw, bounds+0, bounds+1);
	move_box_in_window(iy, dy, 0, ih, bounds+2, bounds+3);
}
static void jh_calchist(int *jhist, int nbits, Image const *image, int x1, int x2, int y1, int y2)
{
	int idx, half[3], kx, ky;

	half[0]=1<<(image->depth[0]-1);
	half[1]=1<<(image->depth[1]-1);
	half[2]=1<<(image->depth[2]-1);
	x1=CLAMP(0, x1, image->iw);
	x2=CLAMP(0, x2, image->iw);
	y1=CLAMP(0, y1, image->ih);
	y2=CLAMP(0, y2, image->ih);
	switch(space_not_color)
	{
	case 0:/*show correlation in color*/
		for(ky=y1;ky<y2;++ky)
		{
			for(kx=x1;kx<x2;++kx)
			{
				int
					r=(((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|0]+half[0])<<nbits>>image->depth[0],
					g=(((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|1]+half[1])<<nbits>>image->depth[1],
					b=(((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx)<<2|2]+half[2])<<nbits>>image->depth[2];
				r=CLAMP(0, r, (1<<nbits)-1);
				g=CLAMP(0, g, (1<<nbits)-1);
				b=CLAMP(0, b, (1<<nbits)-1);
				idx=(b<<nbits|g)<<nbits|r;

				++jhist[idx];
			}
		}
		break;
	case 1:/*show correlation in space x (CURR, W, WW)*/
		for(ky=y1;ky<y2;++ky)
		{
			for(kx=x1;kx<x2;++kx)
			{
				uint8_t
					v2=kx-2>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx-2)<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*WW	*/
					v1=kx-1>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx-1)<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*W	*/
					v0=kx-0>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*ky+kx-0)<<2|1]+half[1])<<nbits>>image->depth[1]):0;/*curr	*/
				v0=CLAMP(0, v0, (1<<nbits)-1);
				v1=CLAMP(0, v1, (1<<nbits)-1);
				v2=CLAMP(0, v2, (1<<nbits)-1);
				idx=(v2<<nbits|v1)<<nbits|v0;

				++jhist[idx];
			}
		}
		break;
	case 2:/*show correlation in space x (CURR, N, NN)*/
		for(ky=y1;ky<y2;++ky)
		{
			for(kx=x1;kx<x2;++kx)
			{
				uint8_t
					v2=ky-2>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*(ky-2)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*NN	*/
					v1=ky-1>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*(ky-1)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*N	*/
					v0=ky-0>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*(ky-0)+kx)<<2|1]+half[1])<<nbits>>image->depth[1]):0;/*curr	*/
				v0=CLAMP(0, v0, (1<<nbits)-1);
				v1=CLAMP(0, v1, (1<<nbits)-1);
				v2=CLAMP(0, v2, (1<<nbits)-1);
				idx=(v2<<nbits|v1)<<nbits|v0;

				++jhist[idx];
			}
		}
		break;
	case 3:/*show correlation in space x (CURR, N, W)*/
		for(ky=y1;ky<y2;++ky)
		{
			for(kx=x1;kx<x2;++kx)
			{
				uint8_t
					v2=kx-1>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw* ky   +kx-1)<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*W		*/
					v1=ky-1>=0?(uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw*(ky-1)+kx  )<<2|1]+half[1])<<nbits>>image->depth[1]):0,/*N		*/
					v0=        (uint8_t)((((const int32_t*)HDR_DATC(image))[(image->iw* ky   +kx  )<<2|1]+half[1])<<nbits>>image->depth[1])  ;/*curr	*/
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
	static const uint8_t triangles[][3]=
	{
		/*6 triangles per cell, reduced resolution - SKEWED*/
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

		/*96 triangles per cell, full resolution - BUGGY*/
#if 0
		{VERTIDX(1, 1, 0), VERTIDX(0, 1, 0), VERTIDX(0, 0, 0)},/*x==0*/
		{VERTIDX(1, 1, 0), VERTIDX(0, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(2, 1, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(2, 1, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 0, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 2, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 0, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 1, 0), VERTIDX(1, 2, 0), VERTIDX(2, 2, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 0)},/*y==0*/
		{VERTIDX(1, 0, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(2, 0, 1), VERTIDX(2, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 0), VERTIDX(2, 0, 0)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 2), VERTIDX(0, 0, 2)},
		{VERTIDX(1, 0, 1), VERTIDX(1, 0, 2), VERTIDX(2, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 0)},/*z==0*/
		{VERTIDX(0, 1, 1), VERTIDX(0, 0, 1), VERTIDX(0, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 2, 1), VERTIDX(0, 2, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 2, 1), VERTIDX(0, 2, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 0, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 0), VERTIDX(0, 2, 0)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 2), VERTIDX(0, 0, 2)},
		{VERTIDX(0, 1, 1), VERTIDX(0, 1, 2), VERTIDX(0, 2, 2)},

		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 0, 1)},/*x==1*/
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(0, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(0, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(2, 0, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 2, 1), VERTIDX(2, 2, 1)},
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 1, 0)},/*y==1*/
		{VERTIDX(1, 1, 1), VERTIDX(0, 1, 1), VERTIDX(0, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(2, 1, 1), VERTIDX(2, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(0, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 0), VERTIDX(2, 1, 0)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(0, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 1, 2), VERTIDX(2, 1, 2)},
		{VERTIDX(1, 1, 1), VERTIDX(1, 0, 1), VERTIDX(1, 0, 0)},/*z==1*/
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
		int kx, ky, kz, kt;
		__m256 mlevel=_mm256_set1_ps(level);
		float gains[3]={0};
		int psize;

		gains[0]=-cubesize/gx;
		gains[0]=-cubesize/gy;
		gains[0]=+cubesize/gz;
		psize=gx*gy;
		for(kz=1;kz<gz;++kz)
		{
			for(ky=1;ky<gy;++ky)
			{
				for(kx=1;kx<gx;++kx)
				{
					int idx=gx*(gy*kz+ky)+kx;
					float/*	vZYX*/
						v000=(float)data[idx-psize-gx-1],/*8 vertices*/
						v001=(float)data[idx-psize-gx  ],
						v010=(float)data[idx-psize   -1],
						v011=(float)data[idx-psize     ],
						v100=(float)data[idx      -gx-1],
						v101=(float)data[idx      -gx  ],
						v110=(float)data[idx         -1],
						v111=(float)data[idx           ];

					__m256 vec=_mm256_set_ps(v111, v110, v101, v100, v011, v010, v001, v000);
					int cond=_mm256_movemask_ps(_mm256_cmp_ps(vec, mlevel, _CMP_LT_OQ));
					if(!cond||cond==0xFF)/*skip if all vertices are below/above level*/
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
							v00m=(v000+v001)*0.5f,/*12 edge-bisectors*/
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
					
							v0mm=(v00m+v01m)*0.5f,/*6 face-centers*/
							v1mm=(v10m+v11m)*0.5f,
							vmm0=(v0m0+v1m0)*0.5f,
							vmm1=(v0m1+v1m1)*0.5f,
							vm0m=(vm00+vm01)*0.5f,
							vm1m=(vm10+vm11)*0.5f,
					
							vmmm=(v0mm+v1mm)*0.5f;/*1 cube center*/
						float comps[3][3]={0}, vertices[27][4];
						int j;

						comps[0][0]=(float)(kx-1); comps[0][1]=(float)kx-0.5f; comps[0][2]=(float)kx;
						comps[1][0]=(float)(ky-1); comps[1][1]=(float)ky-0.5f; comps[1][2]=(float)ky;
						comps[2][0]=(float)(kz-1); comps[2][1]=(float)kz-0.5f; comps[2][2]=(float)kz;

						j=0;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][0]; vertices[j][3]=v000; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][0]; vertices[j][3]=v00m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][0]; vertices[j][3]=v001; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][0]; vertices[j][3]=v0m0; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][0]; vertices[j][3]=v0mm; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][0]; vertices[j][3]=v0m1; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][0]; vertices[j][3]=v010; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][0]; vertices[j][3]=v01m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][0]; vertices[j][3]=v011; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm00; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm0m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm01; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][1]; vertices[j][3]=vmm0; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][1]; vertices[j][3]=vmmm; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][1]; vertices[j][3]=vmm1; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm10; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm1m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][1]; vertices[j][3]=vm11; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][2]; vertices[j][3]=v100; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][2]; vertices[j][3]=v10m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][0]; vertices[j][2]=comps[2][2]; vertices[j][3]=v101; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][2]; vertices[j][3]=v1m0; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][2]; vertices[j][3]=v1mm; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][1]; vertices[j][2]=comps[2][2]; vertices[j][3]=v1m1; ++j;
						vertices[j][0]=comps[0][0]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][2]; vertices[j][3]=v110; ++j;
						vertices[j][0]=comps[0][1]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][2]; vertices[j][3]=v11m; ++j;
						vertices[j][0]=comps[0][2]; vertices[j][1]=comps[1][2]; vertices[j][2]=comps[2][2]; vertices[j][3]=v111; ++j;
						for(kt=0;kt<(int)_countof(triangles);++kt)
						{
							const uint8_t (*tr)[3];
							float (*v[4])[4], *edge, t;
							int mid;

							tr=triangles+kt;
							v[0]=vertices+tr[0][0];
							v[1]=vertices+tr[0][1];
							v[2]=vertices+tr[0][2];
							v[3]=0;
							if(v[0][0][3]>v[1][0][3])SWAPVAR(v[0], v[1], v[3]);
							if(v[0][0][3]>v[2][0][3])SWAPVAR(v[0], v[2], v[3]);
							if(v[1][0][3]>v[2][0][3])SWAPVAR(v[1], v[2], v[3]);
							if(v[0][0][3]<level&&level<v[2][0][3])
							{
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
static void chart_jointhist_update(Image const *image, uint32_t *txid)
{
	int nlevels=1<<jointhist_nbits, hsize=nlevels*nlevels*nlevels;
	int *jhist;
	int bounds[4]={0};
	int _histmax, k;

	if(jointhist)
	{
		if((int)jointhist->count<hsize)
			ARRAY_APPEND(jointhist, 0, hsize-jointhist->count, 1, 0);
	}
	else
		ARRAY_ALLOC(int, jointhist, 0, hsize, 0, 0);
	jhist=(int*)HDR_DATA(jointhist);
	memset(jhist, 0, jointhist->esize*jointhist->count);
	jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, image->iw, image->ih, bounds);
	jh_calchist(jhist, jointhist_nbits, image, bounds[0], bounds[1], bounds[2], bounds[3]);
	jhc_marchingcubes(&jhc_mesh, jhist, 1<<jointhist_nbits, 1<<jointhist_nbits, 1<<jointhist_nbits, jhc_level, jh_cubesize);
	if(!jhc_gpubuf)
	{
		glGenBuffers(1, &jhc_gpubuf);
		if(!jhc_gpubuf)
			CRASH("Alloc error");
	}
	glBindBuffer(GL_ARRAY_BUFFER, jhc_gpubuf);	CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, (int)(jhc_mesh->count*jhc_mesh->esize), HDR_DATA(jhc_mesh), GL_STATIC_DRAW);	CHECK_GL(error);
	memset(jhist, 0, jointhist->esize*jointhist->count);
	jh_calchist(jhist, jointhist_nbits, image, 0, image->iw, 0, image->ih);
	_histmax=0;
	for(k=0;k<hsize;++k)
	{
		if(!_histmax||_histmax<jhist[k])
			_histmax=jhist[k];
	}
	if(_histmax)
	{
		for(k=0;k<hsize;++k)
			((int*)HDR_DATA(jointhist))[k]=(uint8_t)(jhist[k]*255/_histmax);
	}
	if(!*txid)
		glGenTextures(nlevels, txid);
	for(k=0;k<nlevels;++k)
	{
		const uint8_t *hk=HDR_DATC(jointhist)+((size_t)k<<(jointhist_nbits<<1));/*k*(2^nbits)^2*/
		send_texture_pot_grey(txid[k], hk, nlevels, nlevels, 1);
	}
}
static void analysis_update(Image const *image)
{
	int bounds[4]={0};/*x1, x2, y1, y2*/
	int k;

	jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, image->iw, image->ih, bounds);
	c18_analyze(image, bounds[0], bounds[1], bounds[2], bounds[3], &analysis_info);
	
	++pred_hist[0][analysis_info.predidx[0]];
	++pred_hist[1][analysis_info.predidx[1]];
	++pred_hist[2][analysis_info.predidx[2]];
	if(pred_hist[0][analysis_info.predidx[0]]>=512)
	{
		for(k=0;k<PRED_COUNT;++k)
			pred_hist[0][k]>>=1;
	}
	if(pred_hist[1][analysis_info.predidx[1]]>=512)
	{
		for(k=0;k<PRED_COUNT;++k)
			pred_hist[1][k]>>=1;
	}
	if(pred_hist[2][analysis_info.predidx[2]]>=512)
	{
		for(k=0;k<PRED_COUNT;++k)
			pred_hist[2][k]>>=1;
	}
}

void apply_transform(Image **pimage, int tid, int hasRCT)
{
	Image *image=*pimage;

	switch(tid)
	{
#if 0
	case CT_FWD_ADAPTIVE:		rct_adaptive((char*)image, iw, ih, 1);			break;
	case CT_INV_ADAPTIVE:		rct_adaptive((char*)image, iw, ih, 0);			break;
#endif
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
#if 0
	case CT_FWD_CrCgCb:		colortransform_CrCgCb_R(image, 1);			break;
	case CT_INV_CrCgCb:		colortransform_CrCgCb_R(image, 0);			break;
#endif
	case CT_FWD_Pei09:		colortransform_Pei09(image, 1);				break;
	case CT_INV_Pei09:		colortransform_Pei09(image, 0);				break;
	case CT_FWD_SubG_OPT:		colortransform_optsub(image, 1);			break;
	case CT_INV_SubG_OPT:		colortransform_optsub(image, 0);			break;
	case CT_FWD_SUBGREEN:		colortransform_subtractgreen(image, 1);			break;
	case CT_INV_SUBGREEN:		colortransform_subtractgreen(image, 0);			break;
	case CT_FWD_JPEG2000:		colortransform_JPEG2000(image, 1);			break;
	case CT_INV_JPEG2000:		colortransform_JPEG2000(image, 0);			break;
	case CT_FWD_JPEG2000_MA:	colortransform_JPEG2000_MA(image, 1);			break;
	case CT_INV_JPEG2000_MA:	colortransform_JPEG2000_MA(image, 0);			break;
	case CT_FWD_NBLI:		colortransform_NBLI(image, 1);				break;
	case CT_INV_NBLI:		colortransform_NBLI(image, 0);				break;
#if 0
	case CT_FWD_J2K2:		colortransform_J2K2(image, 1);				break;
	case CT_INV_J2K2:		colortransform_J2K2(image, 0);				break;
	case CT_FWD_YRGB_v1:		rct_yrgb_v1(image, 1);					break;
	case CT_INV_YRGB_v1:		rct_yrgb_v1(image, 0);					break;
	case CT_FWD_YRGB_v2:		rct_yrgb_v2(image, 1);					break;
	case CT_INV_YRGB_v2:		rct_yrgb_v2(image, 0);					break;
	case CT_FWD_CMYK:		ct_cmyk_fwd(image);					break;
	case CT_INV_CMYK_DUMMY:									break;
#endif
	case CT_FWD_YCbCr:		colortransform_lossy_YCbCr(image, 1);			break;
	case CT_INV_YCbCr:		colortransform_lossy_YCbCr(image, 0);			break;
	case CT_FWD_XYB:		colortransform_lossy_XYB(image, 1);			break;
	case CT_INV_XYB:		colortransform_lossy_XYB(image, 0);			break;
#if 0
	case CT_FWD_MATRIX:		colortransform_lossy_matrix(image, 1);			break;
	case CT_INV_MATRIX:		colortransform_lossy_matrix(image, 0);			break;
	case CT_FWD_JPEG2000:		colortransform_jpeg2000_fwd((char*)image, iw, ih);	break;
	case CT_INV_JPEG2000:		colortransform_jpeg2000_inv((char*)image, iw, ih);	break;
	case CT_FWD_XGZ:		colortransform_xgz_fwd((char*)image, iw, ih);		break;
	case CT_INV_XGZ:		colortransform_xgz_inv((char*)image, iw, ih);		break;
	case CT_FWD_XYZ:		colortransform_xyz_fwd((char*)image, iw, ih);		break;
	case CT_INV_XYZ:		colortransform_xyz_inv((char*)image, iw, ih);		break;
	case CT_FWD_EXP:		colortransform_exp_fwd((char*)image, iw, ih);		break;
	case CT_INV_EXP:		colortransform_exp_inv((char*)image, iw, ih);		break;
	case CT_FWD_ADAPTIVE:		colortransform_adaptive((char*)image, iw, ih, 1);	break;
	case CT_INV_ADAPTIVE:		colortransform_adaptive((char*)image, iw, ih, 0);	break;
#endif
	case CT_FWD_CUSTOM:		rct_custom(image, 1, rct_custom_params);		break;
	case CT_INV_CUSTOM:		rct_custom(image, 0, rct_custom_params);		break;
	case CT_FWD_HISTMATCH:		ct_histmatch(image, 1);					break;
	case CT_INV_HISTMATCH:		ct_histmatch(image, 0);					break;
#if 0
	case CT_FWD_QUAD:		colortransform_quad((char*)image, iw, ih, 1);		break;
	case CT_INV_QUAD:		colortransform_quad((char*)image, iw, ih, 1);		break;

	case ST_PREPROC_GRAD:		preproc_grad((char*)image, iw, ih);			break;
	case ST_PREPROC_X:		preproc_x((char*)image, iw, ih);			break;
	case ST_PREPROC_X2:		preproc_x2((char*)image, iw, ih);			break;

	case ST_FWD_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 1);break;
	case ST_INV_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 0);break;
	case ST_FWD_CFL:		pred_cfl((char*)image, iw, ih, 1);			break;
	case ST_INV_CFL:		pred_cfl((char*)image, iw, ih, 0);			break;
	case ST_FWD_HYBRID3:		pred_hybrid_fwd((char*)image, iw, ih);			break;
	case ST_INV_HYBRID3:		pred_hybrid_inv((char*)image, iw, ih);			break;
				
	case ST_FWD_CUSTOM2:		custom2_apply((char*)image, iw, ih, 1, &c2_params);	break;
	case ST_INV_CUSTOM2:		custom2_apply((char*)image, iw, ih, 0, &c2_params);	break;
#endif
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
			intptr_t k, n;

			for(k=0, n=(ptrdiff_t)im0->iw*im0->ih*4;k<n;k+=4)
			{
				((int32_t*)HDR_DATA(im1))[k+0]-=((int32_t*)HDR_DATA(im0))[k+0];
				((int32_t*)HDR_DATA(im1))[k+1]-=((int32_t*)HDR_DATA(im0))[k+1];
				((int32_t*)HDR_DATA(im1))[k+2]-=((int32_t*)HDR_DATA(im0))[k+2];
			}
		}
		break;
	case ST_SSIM:
		measure_ssim_map(im0, im1);
		break;
	case ST_FWD_CC:			conv_custom(image);					break;
	case ST_INV_CC:			conv_custom(image);					break;
#if 0
	case ST_FWD_CUSTOM4:		custom4_apply((char*)image, iw, ih, 1, &c4_params);	break;
	case ST_INV_CUSTOM4:		custom4_apply((char*)image, iw, ih, 0, &c4_params);	break;
	case ST_FWD_KALMAN:		kalman_apply((char*)image, iw, ih, 1);			break;
	case ST_INV_KALMAN:		kalman_apply((char*)image, iw, ih, 0);			break;
	case ST_FWD_CUSTOM2:		pred_custom2_apply((char*)image, iw, ih, 1);		break;
	case ST_INV_CUSTOM2:		pred_custom2_apply((char*)image, iw, ih, 0);		break;
	case ST_FWD_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 1);break;
	case ST_INV_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 0);break;
	case ST_FWD_LEARNED:		pred_learned_v4((char*)image, iw, ih, 1);		break;
	case ST_INV_LEARNED:		pred_learned_v4((char*)image, iw, ih, 0);		break;
#ifdef ALLOW_OPENCL
	case ST_FWD_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 1);			break;
	case ST_INV_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 0);			break;
#endif
	case ST_FWD_DIFF2D:		pred_diff2d_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_DIFF2D:		pred_diff2d_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_HPF:		pred_hpf_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_HPF:		pred_hpf_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_GRAD2:		pred_grad2_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_GRAD2:		pred_grad2_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_ADAPTIVE:		pred_adaptive((char*)image, iw, ih, 3, 4, 1);		break;
	case ST_INV_ADAPTIVE:		pred_adaptive((char*)image, iw, ih, 3, 4, 0);		break;
#endif
	case ST_FWD_NBLIC:		pred_nblic(image, 1);					break;
	case ST_INV_NBLIC:		pred_nblic(image, 0);					break;
	case ST_FWD_CALIC:		pred_calic(image, 1, pred_ma_enabled);			break;
	case ST_INV_CALIC:		pred_calic(image, 0, pred_ma_enabled);			break;
	case ST_FWD_LEGALLCG:		pred_LeGallCG(image, 1);				break;
	case ST_INV_LEGALLCG:		pred_LeGallCG(image, 0);				break;
	case ST_FWD_WP:			pred_jxl_apply(image, 1, pred_ma_enabled, jxlparams_i16);break;
	case ST_INV_WP:			pred_jxl_apply(image, 0, pred_ma_enabled, jxlparams_i16);break;
#if 0
	case ST_FWD_WP2:		pred_divfreeWP(image, 1);				break;
	case ST_INV_WP2:		pred_divfreeWP(image, 0);				break;
#endif
	case ST_FWD_WPU:		pred_WPU(image, 1);					break;
	case ST_INV_WPU:		pred_WPU(image, 0);					break;
#if 0
	case ST_FWD_DEFERRED:		pred_wp_deferred(image, 1);				break;
	case ST_INV_DEFERRED:		pred_wp_deferred(image, 0);				break;
#endif
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
#if 0
	case ST_FWD_LPC:		pred_lpc(image, 1);					break;
	case ST_INV_LPC:		pred_lpc(image, 0);					break;
#endif
	case ST_FWD_SIGNBIT:		pred_signbit(image, 1);					break;
	case ST_FWD_MAGNITUDE:		pred_signbit(image, 0);					break;
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
	case ST_FWD_EXTRAP:		pred_extrap(image);					break;
	case ST_FWD_L1DITHER:		pred_l1dither(image);					break;
	case ST_FWD_QUANT:		pred_quantize(image, 1);				break;
	case ST_INV_QUANT:		pred_quantize(image, 0);				break;
	case ST_FWD_ADAQUANT:		pred_adaquant(image, 1);				break;
	case ST_INV_ADAQUANT:		pred_adaquant(image, 0);				break;
	case ST_FWD_SEL4:		pred_sel4(image, 1);					break;
	case ST_INV_SEL4:		pred_sel4(image, 0);					break;
	case ST_FWD_SELECT:		pred_select(image, 1);					break;
	case ST_INV_SELECT:		pred_select(image, 0);					break;
	case ST_FWD_PAETH:		pred_paeth(image, 1);					break;
	case ST_INV_PAETH:		pred_paeth(image, 0);					break;
	case ST_FWD_SELECT4:		pred_select4(image, 1);					break;
	case ST_INV_SELECT4:		pred_select4(image, 0);					break;
	case ST_FWD_CGPLUS:		pred_cgplus(image, 1);					break;
	case ST_INV_CGPLUS:		pred_cgplus(image, 0);					break;
	case ST_FWD_CG422:		pred_CG422(image, 1);					break;
	case ST_INV_CG422:		pred_CG422(image, 0);					break;
	case ST_FWD_CG420:		pred_CG420(image, 1);					break;
	case ST_INV_CG420:		pred_CG420(image, 0);					break;
	case ST_FWD_AV2:		pred_av2(image, 1);					break;
	case ST_INV_AV2:		pred_av2(image, 0);					break;
#if 0
	case ST_FWD_MIX2:		pred_mix2(image, 1);					break;
	case ST_INV_MIX2:		pred_mix2(image, 0);					break;
#endif
	case ST_FWD_GRAY:		pred_gray(image, 1);					break;
	case ST_INV_GRAY:		pred_gray(image, 0);					break;
	case ST_FWD_MIXR:		pred_mixR(image, 1);					break;
	case ST_INV_MIXR:		pred_mixR(image, 0);					break;
	case ST_FWD_MIXN:		pred_mixN(image, 1);					break;
	case ST_INV_MIXN:		pred_mixN(image, 0);					break;
#if 0
	case ST_FWD_MIXNC:		pred_mixNC(image, 1);					break;
	case ST_INV_MIXNC:		pred_mixNC(image, 0);					break;
	case ST_FWD_BESTN:		pred_bestN(image, 1);					break;
	case ST_INV_BESTN:		pred_bestN(image, 0);					break;
	case ST_FWD_RLS:		pred_rls(image, 1);					break;
	case ST_INV_RLS:		pred_rls(image, 0);					break;
	case ST_FWD_AV3:		pred_av3(image, 1);					break;
	case ST_INV_AV3:		pred_av3(image, 0);					break;
#endif
	case ST_FWD_WGRAD:		pred_wgrad(image, 1, hasRCT);				break;
	case ST_INV_WGRAD:		pred_wgrad(image, 0, hasRCT);				break;
#if 0
	case ST_FWD_WGRAD2:		pred_wgrad2(image, 1);					break;
	case ST_INV_WGRAD2:		pred_wgrad2(image, 0);					break;
#endif
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
#if 0
	case ST_FWD_WMIX:		pred_wmix(image, 1);					break;
	case ST_INV_WMIX:		pred_wmix(image, 0);					break;
#endif
	case ST_FWD_TABLE:		pred_table(image, 1);					break;
	case ST_INV_TABLE:		pred_table(image, 0);					break;
	case ST_FWD_AWAV:		pred_awav(image, 1);					break;
	case ST_INV_AWAV:		pred_awav(image, 0);					break;
	case ST_FWD_LWAV:		pred_lwav(image, 1);					break;
	case ST_INV_LWAV:		pred_lwav(image, 0);					break;
#if 0
	case ST_FWD_ECOEFF:		pred_ecoeff(image, 1, pred_ma_enabled);			break;
	case ST_INV_ECOEFF:		pred_ecoeff(image, 0, pred_ma_enabled);			break;
#endif
	case ST_FWD_AV4:		pred_av4(image, 1, pred_ma_enabled);			break;
	case ST_INV_AV4:		pred_av4(image, 0, pred_ma_enabled);			break;
#if 0
	case ST_FWD_MULTISTAGE:		pred_multistage(image, 1, pred_ma_enabled);		break;
	case ST_INV_MULTISTAGE:		pred_multistage(image, 0, pred_ma_enabled);		break;
#endif
	case ST_FWD_ZIPPER:		pred_zipper(image, 1, pred_ma_enabled);			break;
	case ST_INV_ZIPPER:		pred_zipper(image, 0, pred_ma_enabled);			break;
	case ST_FWD_P3:			pred_separate(image, 1, pred_ma_enabled);		break;
	case ST_INV_P3:			pred_separate(image, 0, pred_ma_enabled);		break;
#if 0
	case ST_FWD_DIR:		pred_dir(image, 1, pred_ma_enabled);			break;
	case ST_INV_DIR:		pred_dir(image, 0, pred_ma_enabled);			break;
#endif
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
	case ST_FWD_JMJ:		pred_jmj_apply((char*)image, iw, ih, 1);		break;
	case ST_INV_JMJ:		pred_jmj_apply((char*)image, iw, ih, 0);		break;
	case ST_FWD_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 1);		break;
	case ST_INV_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 0);		break;
	case ST_FWD_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 1);		break;
	case ST_INV_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 0);		break;
	case ST_FWD_MEDIAN:		pred_median_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_MEDIAN:		pred_median_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_DCT3PRED:		pred_dct3_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_DCT3PRED:		pred_dct3_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_PATHPRED:		pred_path_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_PATHPRED:		pred_path_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_GRADPRED:		pred_grad_fwd((char*)image, iw, ih, 3, 4);		break;
	case ST_INV_GRADPRED:		pred_grad_inv((char*)image, iw, ih, 3, 4);		break;
	case ST_FWD_GRAD2:		pred_grad2((char*)image, iw, ih, 1);			break;
	case ST_INV_GRAD2:		pred_grad2((char*)image, iw, ih, 0);			break;
	case ST_FWD_CTX:		pred_ctx((char*)image, iw, ih, 1);			break;
	case ST_INV_CTX:		pred_ctx((char*)image, iw, ih, 0);			break;
	case ST_FWD_C03:		pred_c03((char*)image, iw, ih, 1);			break;
	case ST_INV_C03:		pred_c03((char*)image, iw, ih, 0);			break;
	case ST_FWD_C10:		pred_c10((char*)image, iw, ih, 1);			break;
	case ST_INV_C10:		pred_c10((char*)image, iw, ih, 0);			break;
	case ST_FWD_C20:		pred_c20((char*)image, iw, ih, 1);			break;
	case ST_INV_C20:		pred_c20((char*)image, iw, ih, 0);			break;
	case ST_FWD_WU97:		pred_wu97((char*)image, iw, ih, 1);			break;
	case ST_INV_WU97:		pred_wu97((char*)image, iw, ih, 0);			break;
	case ST_FWD_BITWISE:		pred_bitwise((char*)image, iw, ih, 1);			break;
	case ST_INV_BITWISE:		pred_bitwise((char*)image, iw, ih, 0);			break;
	case ST_FWD_SHUFFLE:		shuffle((char*)image, iw, ih, 1);			break;
	case ST_INV_SHUFFLE:		shuffle((char*)image, iw, ih, 0);			break;
	case ST_FWD_SPLIT:		image_split_fwd((char*)image, iw, ih);			break;
	case ST_INV_SPLIT:		image_split_inv((char*)image, iw, ih);			break;
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
#if 0
	case ST_FWD_GRAD_DWT:
	case ST_INV_GRAD_DWT:
	case ST_FWD_EXPDWT:
	case ST_INV_EXPDWT:
	case ST_FWD_CUSTOM_DWT:
	case ST_INV_CUSTOM_DWT:
#endif
		{
			ArrayHandle sizes;
			int *temp, kc, vmin, vmax;

			sizes=dwt2d_gensizes(image->iw, image->ih, 3, 3, 3);
			temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
			for(kc=0;kc<4;++kc)
			{
				vmin=-(1<<(image->depth[kc]-1)>>1);
				vmax=(1<<(image->depth[kc]-1)>>1)-1;
				if(!image->depth[kc])
					continue;
				switch(tid)
				{
				case ST_FWD_LAZY:      dwt2d_lazy_fwd		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_LAZY:      dwt2d_lazy_inv		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, vmin, vmax);break;
				case ST_FWD_HAAR:      dwt2d_haar_fwd		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_HAAR:      dwt2d_haar_inv		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, vmin, vmax);break;
				case ST_FWD_SQUEEZE:   dwt2d_squeeze_fwd	(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_SQUEEZE:   dwt2d_squeeze_inv	(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, vmin, vmax);break;
				case ST_FWD_LEGALL53:  dwt2d_legall53_fwd	(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_LEGALL53:  dwt2d_legall53_inv	(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, vmin, vmax);break;
				case ST_FWD_CDF97:     dwt2d_cdf97_fwd		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_CDF97:     dwt2d_cdf97_inv		(((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, vmin, vmax);break;
#if 0
				case ST_FWD_GRAD_DWT:  dwt2d_grad_fwd   (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_INV_GRAD_DWT:  dwt2d_grad_inv   (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp);break;
				case ST_FWD_EXPDWT:    dwt2d_exp_fwd    (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, customparam_st);break;/*TODO use customdwtparams instead of sharing allcustomparam_st*/
				case ST_INV_EXPDWT:    dwt2d_exp_inv    (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, customparam_st);break;
				case ST_FWD_CUSTOM_DWT:dwt2d_custom_fwd (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, customparam_st);break;
				case ST_INV_CUSTOM_DWT:dwt2d_custom_inv (((int32_t*)HDR_DATA(image))+kc, (DWTSize*)HDR_DATA(sizes), 0, (int)sizes->count, 4, temp, customparam_st);break;
#endif
				}
				if(tid!=ST_FWD_LAZY&&tid!=ST_INV_LAZY)
				{
					int inv=tid&1;
					image->depth[kc]+=(char)(!inv-inv);
					UPDATE_MAX(image->depth[kc], image->src_depth[kc]);
					UPDATE_MIN(image->depth[kc], 18);
				}
			}
			array_free(&sizes);
			free(temp);
		}
		break;
#if 0
	case ST_FWD_DEC_DWT:   dwt2d_dec_fwd((char*)image, iw, ih);	break;
	case ST_INV_DEC_DWT:   dwt2d_dec_inv((char*)image, iw, ih);	break;
#endif
	}/*switch*/
}
static char wnd_title[4096];
void apply_selected_transforms(Image **pimage, int rct_only, int applyfwd, int applyinv)
{
	volatile double t=0;
	int hasRCT=0;
	int k;

	if(loud_transforms)
		t=time_sec();
	if(!transforms)
		return;
	for(k=0;k<(int)transforms->count;++k)
	{
		uint8_t tid=HDR_DATA(transforms)[k];
		if(tid==CST_COMPARE)
		{
			uint8_t tid2=CST_COMPARE;
			for(k=0;k<(int)transforms->count;++k)/*look for two spatial transforms*/
			{
				tid=HDR_DATA(transforms)[k];
				if(tid>CST_COMPARE)
					break;
			}
			for(++k;k<(int)transforms->count;++k)
			{
				tid2=HDR_DATA(transforms)[k];
				if(tid2>CST_COMPARE)
					break;
			}
			if(tid>CST_COMPARE&&tid2>CST_COMPARE&&tid!=tid2)
			{
				Image *image2=0;
				int res, kp;

				for(k=0;k<(int)transforms->count;++k)/*apply RCTs*/
				{
					uint8_t tid3=HDR_DATA(transforms)[k];
					hasRCT|=tid<CST_INV_SEPARATOR;
					if(tid3==tid||tid3==tid2)
						break;
					apply_transform(pimage, tid3, hasRCT);
				}
				image_copy(&image2, *pimage);
				apply_transform(pimage, tid, 0);
				apply_transform(&image2, tid2, 0);
				res=pimage[0]->iw*pimage[0]->ih<<2;
				for(kp=0;kp<res;++kp)
					((int32_t*)HDR_DATA(pimage[0]))[kp]=abs(((int32_t*)HDR_DATA(pimage[0]))[kp])-abs(((int32_t*)HDR_DATA(image2))[kp]);
				free(image2);
				for(++k;k<(int)transforms->count;++k)
				{
					uint8_t tid3=HDR_DATA(transforms)[k];
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
	for(k=0;k<(int)transforms->count;++k)
	{
		uint8_t tid=HDR_DATA(transforms)[k];
		hasRCT|=tid<CST_INV_SEPARATOR;
		if(rct_only&&tid>CST_INV_SEPARATOR)
			continue;
		if(tid&1?!applyinv:!applyfwd)
			continue;
		apply_transform(pimage, tid, hasRCT);
	}
	if(loud_transforms)
	{
		intptr_t size;
		int nch;

		t=time_sec()-t;
		nch=(im1->depth[0]!=0)+(im1->depth[1]!=0)+(im1->depth[2]!=0)+(im1->depth[3]!=0);
		size=(intptr_t)nch*im1->iw*im1->ih;
		get_window_title(wnd_title, sizeof(wnd_title)-1);

		for(k=0;wnd_title[k]&&wnd_title[k]!='|';++k);

		if(transforms->count)
			set_window_title(
				"%.*s | %12.6f sec  %12.6f MB/s  %12.6f ms/MB"
				, k-(wnd_title[k]=='|')
				, wnd_title
				, t
				, size/(t*1024*1024)
				, t*1024*1024*1000/size
			);
		else
			set_window_title(
				"%.*s"
				, k-(wnd_title[k]=='|')
				, wnd_title
			);
	}
}
static void export_image(Image *im2, uint8_t **pim_export)
{
	intptr_t k, res;
	uint8_t *im_export2;
	void *ptr;
	int shift[4]={0}, c0, c1, c2;

	im_export2=*pim_export;
	ptr=realloc(im_export2, sizeof(char[4])*im2->iw*im2->ih);
	if(!ptr)
	{
		CRASH("Alloc error");
		return;
	}
	im_export2=(uint8_t*)ptr;
	shift[0]=MAXVAR(0, im2->depth[0]-8);
	shift[1]=MAXVAR(0, im2->depth[1]-8);
	shift[2]=MAXVAR(0, im2->depth[2]-8);
	shift[3]=MAXVAR(0, im2->depth[3]-8);
	c0=0;
	c1=1;
	c2=2;
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
		int maxdepth, halfs[4]={0}, abacvis_hi, abacvis_mid;

		maxdepth=im2->depth[c0];
		if(maxdepth<im2->depth[c1])maxdepth=im2->depth[1];
		if(maxdepth<im2->depth[c2])maxdepth=im2->depth[2];
	/*	if(maxdepth<im2->depth[3])maxdepth=im2->depth[3];*/
		if(abacvis_range<2)
		{
			abacvis_low=0;
			abacvis_range=1<<maxdepth;
		}
		halfs[0]=1<<im2->depth[c0]>>1;
		halfs[1]=1<<im2->depth[c1]>>1;
		halfs[2]=1<<im2->depth[c2]>>1;
	/*	halfs[3]=1<<im2->depth[3]>>1;*/
		abacvis_hi=abacvis_low+abacvis_range;
		abacvis_mid=abacvis_low+(abacvis_range>>1);
		for(k=0, res=(intptr_t)im2->iw*im2->ih*4;k<res;k+=4)
		{
			int r=((const int32_t*)HDR_DATC(im2))[k|c0]+halfs[0];
			int g=((const int32_t*)HDR_DATC(im2))[k|c1]+halfs[1];
			int b=((const int32_t*)HDR_DATC(im2))[k|c2]+halfs[2];
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
			im_export2[k|0]=r;
			im_export2[k|1]=g;
			im_export2[k|2]=b;
			im_export2[k|3]=0xFF;
		}
	}
	else
	{
		for(k=0, res=(intptr_t)im2->iw*im2->ih*4;k<res;k+=4)
		{
			im_export2[k|0]=(uint8_t)((((const int32_t*)HDR_DATC(im2))[k|c0]>>shift[c0]))^128;
			im_export2[k|1]=(uint8_t)((((const int32_t*)HDR_DATC(im2))[k|c1]>>shift[c1]))^128;
			im_export2[k|2]=(uint8_t)((((const int32_t*)HDR_DATC(im2))[k|c2]>>shift[c2]))^128;
			im_export2[k|3]=0xFF;
		}
	}
	*pim_export=im_export2;
}
void update_image(void)/*apply selected operations on original image, calculate CRs, and export*/
{
	void *p;
	int maxdepth, maxlevels, k;

	if(!im0)
		return;
	image_copy(&im1, im0);
	apply_selected_transforms(&im1, 0, 1, 1);

	/*do not modify im1 beyond this point*/
	
	export_image(im1, &im_export);
	export_image(im0, &zimage);
	/*image_export_uint8(im1, &im_export, 1, 0);*/

	maxdepth=calc_maxdepth(im1, 0);
	maxlevels=1<<maxdepth;
	if(hist_full_size<maxlevels)
	{
		p=realloc(hist_full, maxlevels*sizeof(int));
		if(!p)
		{
			CRASH("Alloc error");
			return;
		}
		hist_full=(int*)p;
		hist_full_size=maxlevels;
	}
	calc_csize_stateful(im1, hist_full, ch_entropy);
	combCRhist[combCRhist_idx][0]=(float)ch_entropy[0];
	combCRhist[combCRhist_idx][1]=(float)ch_entropy[1];
	combCRhist[combCRhist_idx][2]=(float)ch_entropy[2];
	combCRhist[combCRhist_idx][3]=(float)((ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3])
		/((im1->src_depth[0]>0)+(im1->src_depth[1]>0)+(im1->src_depth[2]>0)+(im1->src_depth[3]>0)));
	for(k=0;k<4;++k)
	{
		/*if(combCRhist_max==1||combCRhist_max<combCRhist[combCRhist_idx][k])*/
		if(fabsf(combCRhist_max-1)<1e-7||combCRhist_max<combCRhist[combCRhist_idx][k])
			combCRhist_max=combCRhist[combCRhist_idx][k];
	}
	combCRhist_idx=(combCRhist_idx+1)%combCRhist_SIZE;
	if(!send_image_separate_subpixels(im1, &txid_separate_r, &txid_separate_g, &txid_separate_b, &txid_separate_a))
		CRASH("Failed to send texture to GPU");
	switch(mode)
	{
	case VIS_SHANNON:
		{
			intptr_t res;
			float *fbuf, vmax;
			double norm;
			int kc, nlevels;

			res=(intptr_t)im1->iw*im1->ih;
			norm=1./(double)res;
			fbuf=(float*)malloc(res*sizeof(float[4]));
			p=realloc(zimage, res*sizeof(char[4]));
			if(!fbuf||!p)
			{
				CRASH("Alloc error");
				return;
			}
			zimage=(uint8_t*)p;
			for(k=0;k<res;++k)
				zimage[k<<2|3]=0xFF;
			for(kc=0;kc<3;++kc)
			{
				nlevels=1<<im1->depth[kc];
				calc_histogram((const int32_t*)HDR_DATC(im1), im1->iw, im1->ih, kc, 0, im1->iw, 0, im1->ih, im1->depth[kc], hist_full, 0);
				for(k=0;k<res;++k)
				{
					int sym, freq;

					sym=((const int32_t*)HDR_DATC(im1))[k<<2|kc]+(nlevels>>1);
					sym=CLAMP(0, sym, nlevels-1);
					freq=hist_full[sym];
					if(freq)
						fbuf[k<<2|kc]=-(float)(log(freq*norm)*(1./M_LN2));
					else
						fbuf[k<<2|kc]=0;
				}
			}
			vmax=0;
			for(k=0;k<res;++k)
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
				for(k=0;k<res;++k)
				{
					for(kc=0;kc<3;++kc)
					{
						float val=fbuf[k<<2|kc]*255/vmax;
						zimage[k<<2|kc]=(uint8_t)CLAMP(0, val, 255);
					}
				}
			}
			else
			{
				for(k=0;k<res;++k)
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
	case VIS_MODEL:
		{
#if 0
	#define PRINT_RESIDUALS
#endif
			void *p2;
			int nlevels, half, mask, kx, ky, kc, idx;
			int psize;
			int *pixels;
#ifdef PRINT_RESIDUALS
			int count=0;
#endif

			modeldepth=im1->depth[0];
			if(modeldepth<im1->depth[1])modeldepth=im1->depth[1];
			if(modeldepth<im1->depth[2])modeldepth=im1->depth[2];
			if(modeldepth<im1->depth[3])modeldepth=im1->depth[3];
			modelnch=im1->nch;
			modelnctx=MODELNCTX;
			modelhistsize=(int)sizeof(int)*modelnch*modelnctx<<modeldepth;
			p=realloc(modelhist, modelhistsize);
			if(!p)
			{
				CRASH("Alloc error");
				return;
			}
			modelhist=(int*)p;
			memset(modelhist, 0, modelhistsize);
			p2=realloc(modelmhist, modelhistsize);
			if(!p2)
			{
				CRASH("Alloc error");
				return;
			}
			modelmhist=(int*)p2;
			memset(modelmhist, 0, modelhistsize);
			nlevels=1<<modeldepth;
			half=nlevels>>1;
			mask=nlevels-1;
			psize=(int)sizeof(int[4*4])*(im1->iw+16);/*4 padded rows * 4 channels max*/
			pixels=(int*)malloc(psize);
			if(!psize)
			{
				CRASH("Alloc error");
				return;
			}
			memset(pixels, 0, psize);
#ifdef PRINT_RESIDUALS
			console_start();
#endif
			for(ky=0, idx=0;ky<im1->ih;++ky)
			{
				int *rows[4]={0};
				int eW[4]={0};

				rows[0]=pixels+((im1->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
				rows[1]=pixels+((im1->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
				rows[2]=pixels+((im1->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
				rows[3]=pixels+((im1->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
				for(kx=0;kx<im1->iw;++kx, idx+=4)
				{
					for(kc=0;kc<modelnch;++kc)
					{
						int ctx, sym;

						ctx=31^LZCNT32(eW[kc]*eW[kc]+1);
						if(ctx>MODELNCTX-1)
							ctx=MODELNCTX-1;
						sym=((const int32_t*)HDR_DATC(im1))[idx|kc];
						sym<<=32-im1->depth[kc];
						sym>>=32-im1->depth[kc];
#ifdef PRINT_RESIDUALS
						if(!kc&&ctx==6&&count++<4096)
							console_log(" %d", sym);
#endif
						++modelhist[(modelnctx*kc+ctx)<<modeldepth|((sym+half)&mask)];
						sym=sym<<1^sym>>31;
						rows[0][kc]=eW[kc]=(2*eW[kc]+(sym<<MODELPREC)+MAXVAR(rows[1][kc+2*4], rows[1][kc+3*4]))>>2;
					}
					rows[0]+=4;
					rows[1]+=4;
					rows[2]+=4;
					rows[3]+=4;
				}
			}
			free(pixels);
			modelstatoverhead=0;
			for(kc=0;kc<modelnch*modelnctx;++kc)
			{
				enum
				{
					PROBBITS=12
				};
				double e, norm, mean, variance, sdev, hsize;
				int *curr_hist, sum, ks, ks2, count, *curr_hist2, msum, cdfW, sum2, CDFlevels;

				curr_hist=modelhist+((ptrdiff_t)kc<<modeldepth);
				sum=0;
				for(ks=0;ks<nlevels;++ks)
					sum+=curr_hist[ks];
				if(!sum)
				{
					modelcsizes[kc]=0;
					modelmsizes[kc]=0;
					continue;
				}
				count=0;
				for(ks=0;ks<nlevels;++ks)
					count+=curr_hist[ks]!=0;
				e=0;
				norm=1./sum;
				for(ks=0;ks<nlevels;++ks)/*simulate 12-bit precision*/
				{
					int freq=curr_hist[ks];
					if(freq)
					{
						int prob=(int)(freq*((int64_t)0x1000-count)/sum)+1;
						e-=freq*log(prob*(1./0x1000));
					}
				}
				e*=1./(8*M_LN2);
				modelcsizes[kc]=e;


				mean=0;
				for(ks=0;ks<nlevels;++ks)/*calc mean*/
				{
					int freq=curr_hist[ks];
					mean+=(double)ks*freq;
				}
				mean*=norm;
				modelmeans[kc]=mean;
				variance=0;
				for(ks=0;ks<nlevels;++ks)/*calc variance*/
				{
					int freq=curr_hist[ks];
					double val=ks-mean;
					variance+=(double)freq*val*val;/*var = E[(X-mean)^2],    E[X] = (p0*X0+p1*X1+...)/pden*/
				}
				variance*=norm;
				sdev=sqrt(variance);
				modelsdevs[kc]=sdev;
				
				curr_hist2=modelmhist+((ptrdiff_t)kc<<modeldepth);
				for(ks=0;ks<nlevels;++ks)/*generate distribution*/
				{
					double fval=(ks-mean)/sdev;
					fval=exp(-fval*fval);
					curr_hist2[ks]=(int)(fval*0x10000);
				}
				msum=0;
				for(ks=0;ks<nlevels;++ks)
					msum+=curr_hist2[ks];
				norm=1./msum;
				e=0;
				for(ks=0;ks<nlevels;++ks)/*calc cross-entropy*/
				{
					int freq=curr_hist[ks];
					int prob=curr_hist2[ks];
					if(prob)
						e-=freq*log(prob*norm);
				}
				modelmsizes[kc]=e*(1./(8*M_LN2));


				hsize=0;
				cdfW=0;
				sum2=0;
				CDFlevels=1<<PROBBITS;
				for(ks=0, ks2=0;ks<nlevels;++ks)/*calc model stat overhead*/
				{
					int sym, freq, cdf, csym;

					sym=((ks>>1^-(ks&1))+half)&mask;/*midpoint->zigzag->edges*/
					freq=curr_hist[sym];
					cdf=sum2*(((uint64_t)1<<PROBBITS)-count)/sum+ks2;
					ks2+=freq!=0;
					csym=cdf-cdfW;
					if(ks&&CDFlevels)/*CDF[0]=0*/
					{
						/*GR*/				/*~8.5 KB*/
						int nbypass=31^LZCNT32(CDFlevels);
						if(ks>1)
							nbypass-=7;
						if(nbypass<0)
							nbypass=0;
						hsize+=(csym>>nbypass)+1+nbypass;
					}
					CDFlevels-=csym;
					cdfW=cdf;
					sum2+=freq;
				}
				modelstatoverhead+=hsize/8;
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
	float x1, x2, y1, y2, z1, z2, cuboid[24];
	int j, k;

	x1=0; x2=size;
	y1=0; y2=size;
	z1=0; z2=size;
	j=0;
	cuboid[j++]=x1; cuboid[j++]=y1; cuboid[j++]=z1;
	cuboid[j++]=x2; cuboid[j++]=y1; cuboid[j++]=z1;
	cuboid[j++]=x2; cuboid[j++]=y2; cuboid[j++]=z1;
	cuboid[j++]=x1; cuboid[j++]=y2; cuboid[j++]=z1;
	cuboid[j++]=x1; cuboid[j++]=y1; cuboid[j++]=z2;
	cuboid[j++]=x2; cuboid[j++]=y1; cuboid[j++]=z2;
	cuboid[j++]=x2; cuboid[j++]=y2; cuboid[j++]=z2;
	cuboid[j++]=x1; cuboid[j++]=y2; cuboid[j++]=z2;
	if(applyRCT)
	{
		ALIGN(4) static char buf[148]={0};
		static Image *image=0;
		float half;

		half=size*0.5f;
		if(!image)
		{
			image=(Image*)buf;
			image->iw=8;
			image->ih=1;
			memset(image->src_depth, 8, sizeof(image->src_depth));
		}
		for(k=0;k<8;++k)
		{
			((int32_t*)HDR_DATA(image))[k<<2|0]=(int)(cuboid[k*3+0]-half);
			((int32_t*)HDR_DATA(image))[k<<2|1]=(int)(cuboid[k*3+1]-half);
			((int32_t*)HDR_DATA(image))[k<<2|2]=(int)(cuboid[k*3+2]-half);
		}
		memcpy(image->depth, image->src_depth, sizeof(image->depth));
		apply_selected_transforms(&image, 1, 1, 1);
		for(k=0;k<8;++k)
		{
			cuboid[k*3+0]=(float)((int32_t*)HDR_DATA(image))[k<<2|0]+half;
			cuboid[k*3+1]=(float)((int32_t*)HDR_DATA(image))[k<<2|1]+half;
			cuboid[k*3+2]=(float)((int32_t*)HDR_DATA(image))[k<<2|2]+half;
		}
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
static void chart_hist_draw(float x1, float x2, float y1, float y2, int cstart, int cend, int color, uint8_t alpha, int *_hist, int *_histmax)
{
	int kc, k2, ks;

	switch(histmode)
	{
	case HISTMODE_LINEAR:
		for(kc=cstart;kc<cend;++kc)
		{
			if(_histmax[kc])
			{
				float dy=(y2-y1)/3.f, histpx=dy/_histmax[kc];
				for(k2=0;k2<256;++k2)
					draw_rect(
						x1+k2*(x2-x1)/256,
						x1+(k2+1)*(x2-x1)/256,
						y1+(kc+1)*dy-_hist[kc<<8|k2]*histpx,
						y1+(kc+1)*dy,
						color?color:alpha<<24|0xFF<<(kc<<3)/*0x80*/
					);
				for(k2=0;k2<256;++k2)
				{
					if(_hist[kc<<8|k2]>(k2>0?_hist[kc<<8|(k2-1)]:0)&&_hist[kc<<8|k2]>(k2<255?_hist[kc<<8|(k2+1)]:0))
						GUIPrint(0, x1+k2*(x2-x1)/256, y1+(kc+1)*dy-_hist[kc<<8|k2]*histpx, 0.8f, "%f", k2*(1./256));
				}
			}
		}
		break;
	case HISTMODE_LINEAR_NOSPIKE:
		for(kc=cstart;kc<cend;++kc)
		{
			int histmax3=0, histmax2=0;
			for(ks=0;ks<256;++ks)
			{
				int freq=_hist[kc<<8|ks];
				if(histmax3<freq)
				{
					histmax2=histmax3;
					histmax3=freq;
				}
				else if(histmax2<freq)
					histmax2=freq;
			}
			if(histmax2)/*2nd max*/
			{
				float dy=(y2-y1)/3.f, histpx=dy/histmax2;
				for(k2=0;k2<256;++k2)
				{
					draw_rect_enqueue(
						&rectangles,
						x1+k2*(x2-x1)/256,
						x1+(k2+1)*(x2-x1)/256,
						y1+(kc+1)*dy-_hist[kc<<8|k2]*histpx,
						y1+(kc+1)*dy
					);
				}
				draw_2d_flush(rectangles, color?color:alpha<<24|0xFF<<(kc<<3), GL_TRIANGLES);
			}
		}
		break;
	case HISTMODE_LOG2:
	case HISTMODE_LOG2_NOSPIKE:
		for(kc=cstart;kc<cend;++kc)
		{
			float hist2[256];
			int hsum=0;
			double invsum=0;
			float histmax3=0, histmax2=0;

			for(ks=0;ks<256;++ks)
				hsum+=_hist[kc<<8|ks];
			invsum=1./hsum;
			for(ks=0;ks<256;++ks)
			{
				int freq=_hist[kc<<8|ks];
				hist2[ks]=0;
				if(freq)
				{
					float f2=hist2[ks]=-(float)(freq*log(freq*invsum)*(1./M_LN2));
					if(histmax3<f2)
					{
						histmax2=histmax3;
						histmax3=f2;
					}
					else if(histmax2<f2)
						histmax2=f2;
				}
			}
			if(histmode==HISTMODE_LOG2_NOSPIKE)
				histmax3=histmax2;
			if(histmax3>0)
			{
				float dy=(y2-y1)/3.f, histpx=dy/histmax3;
				for(k2=0;k2<256;++k2)
				{
					draw_rect_enqueue(
						&rectangles,
						x1+k2*(x2-x1)/256,
						x1+(k2+1)*(x2-x1)/256,
						y1+(kc+1)*dy-hist2[k2]*histpx,
						y1+(kc+1)*dy
					);
				}
				draw_2d_flush(rectangles, color?color:alpha<<24|0xFF<<(kc<<3), GL_TRIANGLES);
			}
		}
		break;
	}
}
static int world2screen(const Camera *pcam, float *w, float *s)
{
	float cp[3]={0}, tmp[3]={0};

	/*
	(tmp)[0] = (w)[0] - (&(*cam).x)[0];
	(tmp)[1] = (w)[1] - (&(*cam).x)[1];
	(tmp)[2] = (w)[2] - (&(*cam).x)[2];
	(cp)[2] = (tmp)[0] * (*cam).cax + (tmp)[1] * (*cam).sax;
	(cp)[0] = (tmp)[0] * (*cam).sax - (tmp)[1] * (*cam).cax;
	(cp)[1] = (cp)[2] * (*cam).say - (tmp)[2] * (*cam).cay;
	(cp)[2] = (cp)[2] * (*cam).cay + (tmp)[2] * (*cam).say;
	*/
	cam_world2cam(*pcam, w, cp, tmp);
	if(cp[2]<=0)
		return 0;
	cam_cam2screen(*pcam, cp, s, wndw*0.5f, wndh*0.5f);
	return 1;
}
typedef struct _NBInfo
{
	int dx, dy;
	const char *label;
} NBInfo;
static void chart_jointhist_draw(void)
{
	draw_AAcuboid_wire(jh_cubesize, 0);
	draw_AAcuboid_wire(jh_cubesize, 1);
	if(wndh)
	{
		static const NBInfo nbs[]=
		{
			{-2, -2, "NNWW"},	/* 0*/
			{-1, -2, "NNW"},	/* 1*/
			{ 0, -2, "NN"},		/* 2*/
			{ 1, -2, "NNE"},	/* 3*/
			{ 2, -2, "NNEE"},	/* 4*/
			{-2, -1, "NWW"},	/* 5*/
			{-1, -1, "NW"},		/* 6*/
			{ 0, -1, "N"},		/* 7*/
			{ 1, -1, "NE"},		/* 8*/
			{ 2, -1, "NEE"},	/* 9*/
			{-2,  0, "WW"},		/*10*/
			{-1,  0, "W"},		/*11*/
			{ 0,  0, "?"},		/*12*/
			{ 0, -3, "NNN"},	/*13*/
			{-3,  0, "WWW"},	/*14*/
		};
		static const int lineindices[][2]=
		{
			/*horizontal lines*/
			{0+0*5, 1+0*5},
			{1+0*5, 2+0*5},
			{2+0*5, 3+0*5},
			{3+0*5, 4+0*5},

			{0+1*5, 1+1*5},
			{1+1*5, 2+1*5},
			{2+1*5, 3+1*5},
			{3+1*5, 4+1*5},
			
			{0+2*5, 1+2*5},
			{1+2*5, 2+2*5},

			/*vertical lines*/
			{0+0*5, 0+1*5},
			{1+0*5, 1+1*5},
			{2+0*5, 2+1*5},
			{3+0*5, 3+1*5},
			{4+0*5, 4+1*5},
			
			{0+1*5, 0+2*5},
			{1+1*5, 1+2*5},
			{2+1*5, 2+2*5},

			{ 2, 13},
			{10, 14},
		};
		int x, y, k;
		float half[3], gain[3], s[_countof(nbs)+1][6]={0}, *pred, vmin, vmax, *s1, *s2, *s3;

		x=mx*im1->iw/wndw;
		y=my*im1->ih/wndh;
		half[0]=(float)(1<<im1->depth[0]>>1);
		half[1]=(float)(1<<im1->depth[1]>>1);
		half[2]=(float)(1<<im1->depth[2]>>1);
		gain[0]=jh_cubesize/(2*half[0]);
		gain[1]=jh_cubesize/(2*half[1]);
		gain[2]=jh_cubesize/(2*half[2]);
		for(k=0;k<(int)_countof(nbs);++k)
		{
			const NBInfo *nb=nbs+k;
			if((uint32_t)(x+nb->dx)<(uint32_t)im1->iw&&(uint32_t)(y+nb->dy)<(uint32_t)im1->ih)
			{
				s[k][0]=gain[0]*((float)((const int32_t*)HDR_DATC(im1))[4*(im1->iw*(y+nb->dy)+x+nb->dx)+0]+half[0]);
				s[k][1]=gain[1]*((float)((const int32_t*)HDR_DATC(im1))[4*(im1->iw*(y+nb->dy)+x+nb->dx)+1]+half[1]);
				s[k][2]=gain[2]*((float)((const int32_t*)HDR_DATC(im1))[4*(im1->iw*(y+nb->dy)+x+nb->dx)+2]+half[2]);
				if(world2screen(&cam, s[k]+0, s[k]+3))
					s[k][5]=1;
			}
		}
#define PRED_CG(P, N, W, NW) P=N+W-NW; vmin=N; vmax=W; if(vmax<vmin)vmax=N, vmin=W; CLAMP2(P, vmin, vmax);
		pred=s[_countof(nbs)];
		PRED_CG(pred[0], s[7][0], s[11][0], s[6][0]);
		PRED_CG(pred[1], s[7][1], s[11][1], s[6][1]);
		PRED_CG(pred[2], s[7][2], s[11][2], s[6][2]);
		if(world2screen(&cam, pred+0, pred+3))
			pred[5]=1;
		for(k=0;k<(int)_countof(lineindices);++k)
		{
			const int *idx=lineindices[k];

			s1=s[idx[0]];
			s2=s[idx[1]];
			if(s1[5]>0.5f&&s2[5]>0.5f)
				draw_line(s1[3], s1[4], s2[3], s2[4], 0xC0000000);
		}
		{
			s1=s[7];
			s2=s[_countof(nbs)];
			s3=s[11];
			if(s1[5]>0.5f&&s2[5]>0.5f)
				draw_line(s1[3], s1[4], s2[3], s2[4], 0xC0F000F0);
			if(s3[5]>0.5f&&s2[5]>0.5f)
				draw_line(s3[3], s3[4], s2[3], s2[4], 0xC0F000F0);
		}
		for(k=0;k<(int)_countof(nbs);++k)
		{
			if(s[k][5]>0.5f)
				GUIPrint(0, s[k][3], s[k][4], 1, "%s", nbs[k].label);
		}
		pred=s[_countof(nbs)];
		if(pred[5]>0.5f)
			GUIPrint(0, pred[3], pred[4], 1, "cg");
	}
	draw_contour3d(&cam, 0, jh_cubesize, 0, jh_cubesize, 0, jh_cubesize, txid_jointhist, 1<<jointhist_nbits, 0.8f);

	if(show_full_image)
	{
		const int crosshaircolor=0xFF000000;
		int bounds[4]={0};/*{x1, x2, y1, y2}*/
		int k;
		float sbounds[4], ratios[2], xmid, ymid;

		jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, im1->iw, im1->ih, bounds);
		display_texture_i(0, wndw, 0, wndh, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 0.5, 0);
		ratios[0]=(float)wndw/im1->iw;
		ratios[1]=(float)wndh/im1->ih;
		for(k=0;k<4;++k)
			sbounds[k]=bounds[k]*ratios[k>>1];
		xmid=(sbounds[0]+sbounds[1])*0.5f;
		ymid=(sbounds[2]+sbounds[3])*0.5f;
		draw_line(sbounds[0], sbounds[2], sbounds[0], sbounds[3], crosshaircolor);/*{x1, y1, x2, y2}*/
		draw_line(sbounds[1], sbounds[2], sbounds[1], sbounds[3], crosshaircolor);
		draw_line(sbounds[0], sbounds[2], sbounds[1], sbounds[2], crosshaircolor);
		draw_line(sbounds[0], sbounds[3], sbounds[1], sbounds[3], crosshaircolor);
		draw_line(xmid, 0, xmid, (float)wndh, crosshaircolor);
		draw_line(0, ymid, (float)wndw, ymid, crosshaircolor);
	}
}


/*active keys turn on timer*/
#define ACTIVE_KEY_LIST\
	AK('W') AK('A') AK('S') AK('D') AK('T') AK('G')\
	AK(KEY_LEFT) AK(KEY_RIGHT) AK(KEY_UP) AK(KEY_DOWN)\
	AK(KEY_ENTER) AK(KEY_BKSP)
int active_keys_pressed=0;

/*mouse*/
char drag=0;
int mx0=0, my0=0;

typedef struct AABBStruct
{
	float x1, x2, y1, y2;
} AABB;
AABB buttons[7]={0};/*0: CT,  1: ST,  2: list of transforms,  3: jxl params,  4: EC,  5: OLS-4*/

int io_init(int argc, char **argv)/*return false to abort*/
{
	(void)argc;/*FIXME open 1 command argument*/
	(void)argv;

	if(tdy*T_COUNT/2*g_uiscale>wndh*1.2f)
		g_uiscale=wndh*1.2f/(tdy*T_COUNT/2);

	set_window_title("Entropy Benchmark");
	glClearColor(1, 1, 1, 1);

	cam_zoomIn(cam, 1);
	cam_turnMouse(cam, 0, 0, mouse_sensitivity);
	memcpy(&cam0, &cam, sizeof(cam));

	set_bk_color(0xA0808080);
#if 1
	{
		const char *str=

			"C:/dataset-GDCC2020-ppm/photo-72.ppm"
		/*	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"*/

		;
		int len=(int)strlen(str);
		Image *im2=image_load(str, len);
		if(im2)
		{
			STR_COPY(currfn, str, len);
			filesize=get_filesize((char*)HDR_DATA(currfn));

			if(im0)
				free(im0);
			im0=im2;
			set_window_title("%s - eBench", (char*)HDR_DATA(currfn));
			update_image();
			if(imagecentered)
				center_image();
		}
	}
#endif
	return 1;
}

const int
	custom_pred_reach=2;
const int
	gui_custom_rct_w=29, gui_custom_rct_h=5,/*characters*/
	gui_custom_pred_w=103+2, gui_custom_pred_h=4,
	gui_ec_width=27,
	gui_ols4_elementchars=9;
int custom_pred_ch_idx=0;/*from {0, 1, 2}*/
void io_resize(void)
{
	AABB *p=buttons;/*FIXME enum*/
	float xstep=tdx*guizoom, ystep=tdy*guizoom;
	p->x1=xstep*2, p->x2=p->x1+xstep*gui_custom_rct_w, p->y1=(float)(wndh>>1), p->y2=p->y1+ystep*gui_custom_rct_h, ++p;/*0: color params - left*/
	p->x1=(float)(wndw>>3), p->x2=p->x1+xstep*gui_custom_pred_w, p->y1=(float)((wndh>>1)+(wndh>>2))-ystep, p->y2=p->y1+ystep*gui_custom_pred_h, ++p;/*1: spatial params - bottom*/
	p->x1=(float)(wndw-300), p->x2=(float)wndw, p->y1=tdy*2, p->y2=p->y1+tdy*T_COUNT/2*g_uiscale, ++p;/*2: transforms list*/
#if 0
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*14, p->y1=(float)((h>>1)+(h>>2))+ystep*4, p->y2=p->y1+ystep, ++p;/*3: clamp bounds*/
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*21, p->y1=(float)((h>>1)+(h>>2))+ystep*5, p->y2=p->y1+ystep, ++p;/*4: learning rate*/
#endif

	p->x1=(float)(wndw>>2), p->x2=p->x1+tdx*11*6, p->y1=(float)((wndh>>1)+(wndh>>2)), p->y2=p->y1+tdy*3, ++p;/*3: jxl params*/

	p->x1=(float)(wndw-450), p->x2=p->x1+tdx*gui_ec_width, p->y1=tdy, p->y2=p->y1+tdy, ++p;/*4: EC method	H.E.M.L..A.0x0000..XXXX_XXX*/

	p->x1=(float)(wndw>>2), p->x2=p->x1+xstep*(OLS4_RMAX<<1|1)*gui_ols4_elementchars, p->y1=(float)(wndh>>1)+10, p->y2=p->y1+ystep*(1+(OLS4_RMAX+1)*(im1?im1->nch:4)), ++p;/*5: OLS-4	DON'T USE p->y2*/

	p->x1=(float)(wndw>>2), p->x2=p->x1+xstep*60, p->y1=(float)(wndh>>1)+10, p->y2=p->y1+ystep*8, ++p;/*6: CUSTOM4*/
	
	if(im1&&im1->iw&&im1->ih&&imagecentered)
		center_image();
}
int io_mousemove(void)/*return true to redraw*/
{
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
			if(mode==VIS_JOINT_HISTOGRAM)/*drag cam*/
			{
				cam_turnMouse(cam, mx-X0, my-Y0, mouse_sensitivity);
				set_mouse(X0, Y0);
			}
			else if(mode==VIS_ANALYSIS)/*drag analysis box*/
			{
				jhc_xbox=mx;
				jhc_ybox=my;
				analysis_update(im1);
			}
			else/*drag image*/
			{
				wpx-=(mx-X0)/imzoom;
				wpy-=(my-Y0)/imzoom;
				set_mouse(X0, Y0);
			}
		}
		return !timer;
	}
	if(profileplotmode>PROFILE_OFF)
		return 1;
	if(mode==VIS_IMAGE&&GET_KEY_STATE('Q'))
		return 1;
	if(mode==VIS_JOINT_HISTOGRAM)
		return !timer;
	return 0;
}
static void click_hittest(int _mx, int _my, int *objidx, int *cellx, int *celly, int *cellidx, AABB **p)
{
	*p=buttons;
	*objidx=0;
	for(*objidx=0;*objidx<(int)_countof(buttons);++*objidx, ++*p)
	{
		/*when these buttons are inactive they shouldn't block the click*/
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
	case 0:/*color transform params*/
		*cellx=(int)floor((_mx-p[0]->x1)*gui_custom_rct_w/(p[0]->x2-p[0]->x1));
		*celly=(int)floor((_my-p[0]->y1)*gui_custom_rct_h/(p[0]->y2-p[0]->y1));
		*cellidx=0;
		break;
	case 1:/*spatial transform params*/
		*cellx=(int)floor((_mx-p[0]->x1)*(custom_pred_reach<<1|1)/(p[0]->x2-p[0]->x1));
		*celly=(int)floor((_my-p[0]->y1)*gui_custom_pred_h/(p[0]->y2-p[0]->y1));
		*cellidx=(custom_pred_reach<<1|1)*(*celly-1)+*cellx;
		break;
	case 2:/*list of transforms*/
		*cellx=(int)floor((_mx-p[0]->x1)*2/(p[0]->x2-p[0]->x1));
		*celly=(int)floor((_my-p[0]->y1)*(T_COUNT/2)/(p[0]->y2-p[0]->y1));
		*cellidx=*celly;
		break;
	case 3:/*jxl params*/
		*cellx=(int)floor((_mx-p[0]->x1)*11/(p[0]->x2-p[0]->x1));
		*celly=(int)floor((_my-p[0]->y1)* 3/(p[0]->y2-p[0]->y1));
		*cellidx=11**celly+*cellx;
		break;
	case 4:
		/*H.E.M.L..A.0x0000..XXXX_XXX*/
		*cellx=(int)floor((_mx-p[0]->x1)*gui_ec_width/(p[0]->x2-p[0]->x1));
		*celly=0;
		*cellidx=*cellx;
		break;
	case 5:
		*cellx=(int)floor((_mx-p[0]->x1)/(tdx*guizoom));
		*celly=(int)floor((_my-p[0]->y1)/(tdy*guizoom));
		*cellidx=0;
		break;
	case 6:
		*cellx=(int)floor((_mx-p[0]->x1)/(tdx*guizoom));
		*celly=(int)floor((_my-p[0]->y1)/(tdy*guizoom));
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
	if(im1)
	{
		int objidx=0, cellx=0, celly=0, cellidx=0;
		AABB *gui_cell;

		if(GET_KEY_STATE(KEY_CTRL))
		{
			g_uiscale+=(2*(forward>0)-1)*0.01f;
			CLAMP2(g_uiscale, 0.01f, 4);
			io_resize();
			return 1;
		}
		gui_cell=buttons;
		click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &gui_cell);
		if(objidx!=-1)
		{
			int sign=(forward>0)-(forward<0);/*abs(forward) is 120*/
			int ch=(int)floor((mx-gui_cell->x1)/(guizoom*tdx));
			switch(objidx)
			{
			case 0:/*color transform params*/
				if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
				{
					/*
					P0  rgb
					r += (-0x0000*g-0x0000*b)>>12
					g += (-0x0000*r-0x0000*b)>>12
					b += (-0x0000*r-0x0000*g)>>12
					g += (-0x0000*r-0x0000*b)>>12
					012345678901234567890123456789
					*/
					if(!celly)
					{
						int x=rct_custom_params[8];
						x=(x+sign+6)%6;
						rct_custom_params[8]=(int16_t)x;
					}
					else if((uint32_t)(cellx-9)<4||(uint32_t)(cellx-18)<4)
					{
						int idx;
						if((uint32_t)(cellx-9)<4)
							idx=(celly-1)<<1|0, cellx-=9;
						else
							idx=(celly-1)<<1|1, cellx-=18;
						rct_custom_params[idx]+=(int16_t)(sign<<((3-cellx)<<2));
					}
					update_image();
				}
				break;
			case 1:/*spatial transform params*/
				if(!transforms_mask[ST_FWD_CUSTOM]&&!transforms_mask[ST_INV_CUSTOM]
					&&!transforms_mask[ST_CONVTEST]
					&&!transforms_mask[ST_CONVTEST2]
					&&!transforms_mask[ST_FWD_CC]&&!transforms_mask[ST_INV_CC]
				)
					break;
				if(!celly)
				{
					/*
					0123456789012345678901234567
					Ch 0  Clamp [+W +NW +N +NE]
					*/
					if(ch==3)
						custom_pred_ch_idx=(custom_pred_ch_idx+sign+3)%3;
					else if(BETWEEN_EXC(13, ch, 16))/*W*/
						custom_clamp[0]=!custom_clamp[0];
					else if(BETWEEN_EXC(16, ch, 20))/*NW*/
						custom_clamp[1]=!custom_clamp[1];
					else if(BETWEEN_EXC(20, ch, 23))/*N*/
						custom_clamp[2]=!custom_clamp[2];
					else if(BETWEEN_EXC(23, ch, 26))/*NE*/
						custom_clamp[3]=!custom_clamp[3];
				}
				else
				{
					/*
					0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
					0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
					-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
					*/
					int col=ch%21, idx=5*(celly-1)+(ch/21);
					if(idx<12&&((uint32_t)(col-3)<7||(uint32_t)(col-13)<7))
					{
						idx=idx<<1|(col>=10);
						{
							/*
							0123456
							00.0000
							*/
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
				break;
			case 3:/*jxl params*/
				if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					int ch2;
					ch=(int)floor((mx-gui_cell->x1)/tdx);
					MODVAR(ch2, ch, 6);
					if(ch2>=2&&ch2<6)
					{
						int delta=sign<<((5-ch2)<<2);
						jxlparams_i16[cellidx]+=(int16_t)delta;
					}
					else if(cellx>=4)
						jxlparams_i16[cellidx]=-jxlparams_i16[cellidx];
					update_image();
				}
				break;
			case 4:/*EC method*/
				switch(cellx)
				{
				/*
				0123456789012345678901234567
				H.E.M.L..A.0x0000..XXXX_XXX
				*/
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
					ec_adaptive^=1;
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
						static int prevmode=0;
						int method1, method2;

						method1=ec_method;
						method2=ec_method;
						method2=(method2+sign+ECTX_COUNT)%ECTX_COUNT;
						ec_method=(EContext)method2;
						if(ec_method==ECTX_ABAC0)
						{
							prevmode=viewmode;
							viewmode=VIEW_C0;
						}
						else if(method1==ECTX_ABAC0)
							viewmode=prevmode;
					}
					break;
				}
				if(ec_expbits<ec_msb+ec_lsb)
					ec_msb=ec_lsb=0;
				update_image();
				break;
			case 5:/*OLS-4*/
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
								if(idx==OLS4_CTXSIZE)/*causality mask*/
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
								/*
								0	X	2	3	4	...
								0	.	-1	-2	-3	...
								*/
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
			case 6:/*CUSTOM4*/
				if(!celly)
				{
					if(cellx==6)/*change layer, stay in same channel*/
					{
						lossyconv_page+=sign<<2;
						lossyconv_page&=15;
					}
					else if(cellx==13)/*change channel, stay in same layer*/
					{
						lossyconv_page=(lossyconv_page&12)|((lossyconv_page+sign)&3);
					}
					else if(cellx==23)
					{
						uint8_t *cell=lossyconv_causalRCT+(lossyconv_page>>2);
						*cell=!*cell;
						update_image();
					}
				}
				else if(celly==6)
				{
					uint8_t *cell;
					if((uint32_t)(cellx-11)<3)
					{
						cell=lossyconv_offset+((size_t)lossyconv_page>>2<<1|0);
						*cell+=sign;
						update_image();
					}
					else if((uint32_t)(cellx-15)<3)
					{
						cell=lossyconv_offset+((size_t)lossyconv_page>>2<<1|1);
						*cell+=sign;
						update_image();
					}
				}
				else if(celly==7)
				{
					uint8_t *cell;
					if((uint32_t)(cellx-11)<3)
					{
						cell=lossyconv_stride+((size_t)lossyconv_page>>2<<1|0);
						*cell+=sign;
						update_image();
					}
					else if((uint32_t)(cellx-15)<3)
					{
						cell=lossyconv_stride+((size_t)lossyconv_page>>2<<1|1);
						*cell+=sign;
						update_image();
					}
				}
				else if((uint32_t)(celly-1)<5)
				{
					/*
					|    000000000011111111112222222222333333
					|    012345678901234567890123456789012345
					|0  "Layer 1/4  C 1/3  cRCT 1"
					|    0000000000111111111122222222223333333333444444444455555555556
					|    0123456789012345678901234567890123456789012345678901234567890
					|1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
					|4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|6  "XY  Offset XXX YYY"
					|7  "    Stride XXX YYY"
					*/
					int c=cellx%12, x=cellx/12;
					int16_t *cell=lossyconv_params+5*5*lossyconv_page+5*(celly-1)+x;
					if(c==10)/*toggle clamp mask*/
					{
						*cell^=1;
						update_image();
					}
					else if((uint32_t)(c-1)<9)
					{
						/*
						123456789
						876.54321
						*/
						int bitidx=9-c;
						bitidx+=bitidx<6;
						*cell=((uint16_t)(*cell+(sign<<bitidx))&~1)|(*cell&1);
						*cell<<=16-9;/*9-bit params*/
						*cell>>=16-9;
						update_image();
					}
				}
				break;
			}/*switch*/
		}/*if*/
		else
			goto normal_operation;
	}
	else
	{
	normal_operation:
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
			else if(GET_KEY_STATE(KEY_SHIFT))/*shift wheel		change cam speed*/
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
			zoom_at(mx, my, mw_fwd?2:1/2.);/*fwd zooms in*/
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
static int parse_nvals_i8(ArrayHandle text, int idx, uint8_t *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);

		neg=HDR_DATC(text)[idx]=='-';
		idx+=neg||HDR_DATC(text)[idx]=='+';/*skip sign*/
		if(HDR_DATC(text)[idx]=='0'&&(HDR_DATC(text)[idx]&0xDF)=='X')/*skip hex prefix*/
			idx+=2;
		end=(char*)HDR_DATA(text)+idx;
		params[k]=(uint8_t)strtol((char*)HDR_DATA(text)+idx, &end, 16);
		idx=(int)(end-(char*)HDR_DATA(text));
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx])&&HDR_DATC(text)[idx]!='-'&&!isdigit(HDR_DATC(text)[idx]);++idx);/*skip comma*/

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

		for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);

		neg=HDR_DATC(text)[idx]=='-';
		idx+=neg||HDR_DATC(text)[idx]=='+';/*skip sign*/
		end=(char*)HDR_DATA(text)+idx;
		params[k]=strtod((char*)HDR_DATA(text)+idx, &end);
		idx=(int)(end-(char*)HDR_DATA(text));
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx])&&HDR_DATC(text)[idx]!='-'&&!isdigit(HDR_DATC(text)[idx]);++idx);/*skip comma*/

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

		for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);

		neg=HDR_DATC(text)[idx]=='-';
		idx+=neg||HDR_DATC(text)[idx]=='+';/*skip sign*/
		if(HDR_DATC(text)[idx]=='0'&&(HDR_DATC(text)[idx]&0xDF)=='X')/*skip hex prefix*/
			idx+=2;
		end=(char*)HDR_DATA(text)+idx;
		params[k]=(int)strtol((char*)HDR_DATA(text)+idx, &end, 16);
		idx=(int)(end-(char*)HDR_DATA(text));
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx])&&HDR_DATC(text)[idx]!='-'&&!isdigit(HDR_DATC(text)[idx]);++idx);/*skip comma*/

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static int parse_nvals_i16(ArrayHandle text, int idx, int16_t *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		int neg;
		char *end;

		for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);

		neg=HDR_DATC(text)[idx]=='-';
		idx+=neg||HDR_DATC(text)[idx]=='+';/*skip sign*/
		if(HDR_DATC(text)[idx]=='0'&&(HDR_DATC(text)[idx]&0xDF)=='X')/*skip hex prefix*/
			idx+=2;
		end=(char*)HDR_DATA(text)+idx;
		params[k]=(int16_t)strtol((char*)HDR_DATA(text)+idx, &end, 16);
		idx=(int)(end-(char*)HDR_DATA(text));
		if(neg)
			params[k]=-params[k];

		for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx])&&HDR_DATC(text)[idx]!='-'&&!isdigit(HDR_DATC(text)[idx]);++idx);/*skip comma*/

		++k;

		if(idx>=(int)text->count&&k<count)
			return idx;
	}
	return idx;
}
static void append_i16_row(ArrayHandle *str, const int16_t *vals, int count)
{
	int k;

	for(k=0;k<count;++k)
	{
		int16_t val=vals[k];
		str_append(str, "%c0x%04X,", val<0?'-':' ', abs(val));
	}
	str_append(str, "\n");
}
int io_keydn(IOKey key, char c)
{
	(void)c;

	switch(key)/*handle shortcuts involving timed keys before the main switch*/
	{
	case 'S':
		if(keyboard[KEY_CTRL]&&im1&&im_export)
		{
			Filter ext[]=
			{
				{"Portable Network Graphics (*.PNG) File", "*.PNG"},
				{"Portable Pixmap (*.PPM) File", "*.PPM"},
			};
			char *filename=dialog_save_file(ext, _countof(ext), "Untitled.PNG");
			if(filename)
			{
				int len=(int)strlen(filename);

				if(!_stricmp(filename+len-4, ".PNG"))
					lodepng_encode_file(filename, im_export, im1->iw, im1->ih, LCT_RGBA, 8);
				else if(!_stricmp(filename+len-4, ".PPM"))
				{
					size_t k, res;
					FILE *fdst;

					fdst=fopen(filename, "wb");
					if(!fdst)
					{
						crash(__FILE__, __LINE__, 0, "Cannot open \"%s\" for writing", currfn);
						return 1;
					}
					fprintf(fdst, "P6\n%d %d\n255\n", im1->iw, im1->ih);
					for(k=0, res=(size_t)4*im1->iw*im1->ih;k<res;k+=4)
						fwrite(im_export+k, 1, 3, fdst);
					fclose(fdst);
				}
				free(filename);
			}
			return 1;
		}
		break;
	case KEY_UP:
	case KEY_DOWN:
		if(mode!=VIS_JOINT_HISTOGRAM)
		{
			if(im1)
			{
				g_dist+=(key==KEY_UP)*2-1;
				if(g_dist<1)
					g_dist=1;
				if(g_dist>31)
					g_dist=31;
				update_image();
				return 1;
			}
		}
		break;
	case KEY_LEFT:
	case KEY_RIGHT:
		if(mode!=VIS_JOINT_HISTOGRAM)
		{
			if(currfn)
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
				int k, currentidx;

				STR_COPY(path, HDR_DATA(currfn), currfn->count);
				for(k=(int)path->count-1;k>=0;--k)
				{
					if(HDR_DATA(path)[k]=='/'||HDR_DATA(path)[k]=='\\')
					{
						HDR_DATA(path)[k+1]=0;
						path->count=(size_t)k+1;
						break;
					}
				}
				filenames=get_filenames((char*)HDR_DATA(path), ext, _countof(ext), 1);
				filteredfn=filter_path((char*)HDR_DATA(currfn), (int)currfn->count, 0);
				if(filenames&&filenames->count)
				{
					currentidx=-1;
					for(k=0;k<(int)filenames->count;++k)
					{
						fn2=(ArrayHandle*)array_at(&filenames, k);
						if(!_stricmp((char*)HDR_DATA(fn2[0]), (char*)HDR_DATA(filteredfn)))
						{
							currentidx=k;
							break;
						}
					}
					{
						Image *im2=0;
						int step=key==KEY_RIGHT?1:-1;
						for(k=currentidx+step;(k=(k+(int)filenames->count)%(int)filenames->count)!=currentidx;k+=step)
						{
							fn2=(ArrayHandle*)array_at(&filenames, k);
							im2=image_load((char*)HDR_DATA(fn2[0]), (int)fn2[0]->count);
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
							array_free(&currfn);
							filesize=get_filesize((char*)HDR_DATA(fn2[0]));
							currfn=filter_path((char*)HDR_DATA(fn2[0]), (int)fn2[0]->count, 0);
							set_window_title("%s - eBench", (char*)HDR_DATA(currfn));
							update_image();
							if(imagecentered)
								center_image();
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
	case 'G':
		if(GET_KEY_STATE(KEY_CTRL))
		{
			ArrayHandle text=paste_from_clipboard(0);
			char *p=(char*)HDR_DATA(text);

			int kx=strtol(p, &p, 10);
			int ky=strtol(p, &p, 10);
			imzoom=64;
			wpx=kx-wndw/(2*imzoom);
			wpy=ky-wndh/(2*imzoom);
			imagecentered=0;
			return 1;
		}
		break;
	default:break;
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
			case 0:/*color transform params*/
				if(key==KEY_RBUTTON)
				{
					/*
					|P0  rgb
					|r += (-0x0000*g-0x0000*b)>>12
					|g += (-0x0000*r-0x0000*b)>>12
					|b += (-0x0000*r-0x0000*g)>>12
					|g += (-0x0000*r-0x0000*b)>>12
					|012345678901234567890123456789
					*/
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
					else if((uint32_t)(cellx-9)<4||(uint32_t)(cellx-18)<4)
					{
						int idx;
						if((uint32_t)(cellx-9)<4)
							idx=(celly-1)<<1|0, cellx-=9;
						else
							idx=(celly-1)<<1|1, cellx-=18;
						rct_custom_params[idx]=0;
						update_image();
						return 1;
					}
				}
				break;
			case 1:/*spatial transforms params*/
				{
					int
						col=guicell?(int)floor((mx-guicell->x1)/(guizoom*tdx)):0,
						line=guicell?(int)floor((my-(guicell->y1+guizoom*tdy))/(guizoom*tdy)):0;
					if(key==KEY_RBUTTON&&cellidx<12)
					{
						/*
						|0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
						|0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
						|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
						*/
						int idx=5*line+col/21;
						col%=21;
						idx=idx<<1|(col>=10);
						custom_params[24*custom_pred_ch_idx+idx]=0;
						update_image();
						return 1;
					}
				}
				break;
			case 2:/*transform list*/
				{
					int
						line=celly;
					if((uint32_t)line<(uint32_t)(T_COUNT/2))
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
			case 3:/*jxl params*/
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
			case 4:/*EC method*/
				if(key==KEY_LBUTTON)
					break;
				switch(cellx)
				{
				/*
				012345678901234567
				H.E.M.L..A.0x0000..XXXX_XXX
				*/
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
					break;
				}
				if(ec_expbits<ec_msb+ec_lsb)
					ec_msb=ec_lsb=0;
				update_image();
				return 1;
			case 5:/*OLS-4*/
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
			case 6:/*CUSTOM4*/
				if((uint32_t)(celly-1)<5)
				{
					/*
					|    000000000011111111112222222222333333
					|    012345678901234567890123456789012345
					|0  "Layer 1/4  C 1/3  cRCT 1"
					|    0000000000111111111122222222223333333333444444444455555555556
					|    0123456789012345678901234567890123456789012345678901234567890
					|1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
					|4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
					|6  "XY  Offset XXX YYY"
					|7  "    Stride XXX YYY"
					*/
					int16_t *ptr=lossyconv_params+5*5*lossyconv_page+5*(celly-1)+cellx/12;
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
					if((uint32_t)(cellx-11)<3)
					{
						lossyconv_offset[lossyconv_page>>2<<1|0]=0;
						update_image();
						return 1;
					}
					if((uint32_t)(cellx-15)<3)
					{
						lossyconv_offset[lossyconv_page>>2<<1|1]=0;
						update_image();
						return 1;
					}
				}
				else if(celly==7)
				{
					if((uint32_t)(cellx-11)<3)
					{
						lossyconv_stride[lossyconv_page>>2<<1|0]=0;
						update_image();
						return 1;
					}
					if((uint32_t)(cellx-15)<3)
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
			if(drag)/*enter mouse control*/
			{
				mx0=mx, my0=my;
				set_mouse(wndw>>1, wndh>>1);
			}
			else/*leave mouse control*/
				set_mouse(mx0, my0);
		}
		break;
	case KEY_MBUTTON:
#if 0
		messageboxa("Click at (%d, %d)\n", mx, my);
#endif
		break;

#define AK(KEY) case KEY:
	ACTIVE_KEY_LIST
#undef  AK
		timer_start(10, TIMER_ID_KEYBOARD);
		break;
		
	case KEY_F1:
		{
			static char msgbuf[4096]={0};
			static int printed=0;

			assert(VIS_IMAGE_TRICOLOR	==0);
#ifdef ENABLE_L1WEIGHTS
			assert(VIS_L1WEIGHTS		==1);
#endif
			assert(VIS_IMAGE		==2);
			assert(VIS_HISTOGRAM		==3);
			assert(VIS_MODEL		==4);
			assert(VIS_ANALYSIS		==5);
			assert(VIS_JOINT_HISTOGRAM	==6);
			assert(VIS_SHANNON		==7);
			if(!printed)/*workaround for 509 limit for string literals*/
			{
#define STR2BUF(STR)\
	do\
	{\
		const char mystr[]=STR;\
		const int mylen=(int)sizeof(mystr)-1;\
		\
		memcpy(msgbuf+printed, mystr, mylen);\
		printed+=mylen;\
		msgbuf[printed]='\0';\
	}while(0)

				STR2BUF(
					"Ctrl O:\t\tOpen image\n"
					"Ctrl S:\t\tSave preview as 8-bit PNG\n"
					"Mouse1/Mouse2:\tAdd/remove transforms in the list\n"
					"Ctrl Mouse1:\tReplace all transforms of this type\n"
					"E:\t\tShow image at 1:1 scale\n"
					"C:\t\tCenter image\n"
					"Ctrl C:\t\tToggle histogram type\n"
					"Ctrl G:\t\tGoto pixel from clipboard\n"
					"Ctrl H:\t\tToggle hex pixel labels\n"
					"J / Shift J:\tToggle single-channel view\n"
					"N:\t\tToggle modular arithmetic in image view\n"
					"Ctrl R:\t\tDisable all transforms\n"
					"Ctrl Space:\tMeasure PSNR\n"
				);
				STR2BUF(
					"Comma/Period:\tSelect context for size estimation\n"
					"Slash:\t\tToggle adaptive histogram\n"
					"[ ]:\t\t(Custom transforms) Select coefficient page\n"
					"Ctrl 1/2/3:\t(Custom predictor) Populate page from channel N\n"
					"Space:\t\t(Custom transforms) Optimize\n"
					"B:\t\tWeights view: Toggle mixing type\n"
					"P:\t\tWeights view: Hold to pause graph\n"
					"L:\t\tWeights view: Set cursor to mouse pointer\n"
					"Ctrl B:\t\tBatch test\n"
					"Ctrl L:\t\tLossy batch test\n"
					"Hold Q:\t\tLossy comparison\n"
					"\n"
					"H:\t\tReset invCR history graph\n"
					"Ctrl C:\t\tCopy data\n"
					"Ctrl V:\t\tPaste data\n"
				);
				STR2BUF(
					"\n"
					"M / Shift M:\tCycles between:\n"
					"\t1: Image tricolor view\n"
					"\t2: Weights view\n"
					"\t3: Image view\n"
					"\t4: Histogram\n"
					"\t5: Model\n"
					"\t6: Analysis\n"
					"\t7: Joint histogram\n"
					"\t8: Shannon view\n"
				);
#undef  STR2BUF
			}
			messagebox(MBOX_OK, "Controls", "%s", msgbuf);
		}
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
		ec_adaptive^=1;
		update_image();
		return 1;
	case 'O':
		if(GET_KEY_STATE(KEY_CTRL))
		{
			ArrayHandle fn2=dialog_open_file(0, 0, 0);
			if(fn2)
			{
				Image *im2=image_load((char*)HDR_DATA(fn2), (int)fn2->count);
				if(im2)
				{
					if(currfn)
						array_free(&currfn);
					currfn=fn2;
					filesize=get_filesize((char*)HDR_DATA(fn2));

					if(im0)
						free(im0);
					im0=im2;
					set_window_title("%s - eBench", (char*)HDR_DATA(currfn));
					update_image();
					if(imagecentered)
						center_image();
				}
				else
					array_free(&fn2);
			}
		}
		return 1;
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
				uint8_t *buf=0;
				image_export_uint8(im1, &buf, 1, 1);/*swap red & blue for WinAPI*/
				if(!buf)
				{
					WARNING("Alloc error");
					return 0;
				}
				{
					int success=copy_bmp_to_clipboard(buf, im1->iw, im1->ih);
					if(!success)
						WARNING("Failed to copy image to clipboard");
				}
				free(buf);
				return 0;
			}
			STR_ALLOC(str, 0);
			if(mode==VIS_MODEL)
			{
				int nhist, nlevels, kh, ks;
				
				nhist=modelnch*modelnctx;
				nlevels=1<<modeldepth;
				str_append(&str, "%-9s", "");
				for(kh=0;kh<nhist;++kh)
					str_append(&str, "   C%d CTX%02d", kh/modelnctx, kh%modelnctx);
				str_append(&str, "\n\n");
				for(ks=0;ks<nlevels;++ks)
				{
					str_append(&str, "[%4d]:  ", ks);
				/*	str_append(&str, "C%d CTX%d:\n", kh/modelnctx, kh%modelnctx);*/
					for(kh=0;kh<nhist;++kh)
						str_append(&str, " %10d", modelhist[(kh<<modeldepth)+ks]);
					str_append(&str, "\n");
				}
			}
			else
			{
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					int kx, ky, kc;

					str_append(&str, "%8d %8.6f %8.6f %8.6f %8.6f\n", ols4_period, ols4_lr[0], ols4_lr[1], ols4_lr[2], ols4_lr[3]);
					for(kc=0;kc<4;++kc)
					{
						for(ky=0;ky<OLS4_RMAX+1;++ky)
						{
							for(kx=0;kx<(OLS4_RMAX<<1|1);++kx)
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
				}/*no 'else' here*/
				if(mode==VIS_IMAGE||mode==VIS_SHANNON)/*copy custom transform value*/
				{
					double invcr[5]={0}, csizes[5]={0}, usize=0;

					invcr[0]=(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3])
						/((im1->src_depth[0]!=0)+(im1->src_depth[1]!=0)+(im1->src_depth[2]!=0)+(im1->src_depth[3]!=0));
					invcr[1]=ch_entropy[0];
					invcr[2]=ch_entropy[1];
					invcr[3]=ch_entropy[2];
					invcr[4]=ch_entropy[3];
					usize=invcr2csizes(invcr, im0->src_depth, im0->iw, im0->ih, 4, csizes);
					str_append(&str,
						"UTYUVA %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f bytes  BPD %8.4f",
						usize,
						csizes[0],
						csizes[1],
						csizes[2],
						csizes[3],
						csizes[4],
						csizes[0]/usize*8
					);
				}
				else if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
				{
					const int width=C3_REACH<<2|2;/*interleaved pixels & errors*/
					const int16_t *coeffs[]=
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
					int kc, k;

					for(kc=0;kc<9;++kc)
					{
						for(k=0;k<C3_REACH;++k)
							append_i16_row(&str, coeffs[kc]+width*k, width);
						append_i16_row(&str, coeffs[kc]+width*C3_REACH, tails[kc]);
						if(kc+1<9)
							str_append(&str, "\n");
					}
				}
				else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					int kx, ky;

					for(ky=0;ky<3;++ky)
					{
						for(kx=0;kx<11;++kx)
						{
							int16_t val=jxlparams_i16[11*ky+kx];
							str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
							if(kx+1==11)
								str_append(&str, "\n");
						}
					}
				}
				else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
				{
					int kx, ky;

					for(ky=0;ky<3;++ky)
					{
						for(kx=0;kx<PW2_NPARAM;++kx)
						{
							int16_t val=pw2_params[PW2_NPARAM*ky+kx];
							str_append(&str, "%c0x%04X,", val<0?'-':' ', abs(val));
							if(kx+1==PW2_NPARAM)
								str_append(&str, "\n");
						}
					}
				}
				else if(transforms_customenabled)
				{
					const int stw=(custom_pred_reach<<1|1)*2;
					int shift=GET_KEY_STATE(KEY_SHIFT);
					int k, kc2;

					for(k=0;k<RCT_CUSTOM_NPARAMS;++k)
					{
						int val=rct_custom_params[k];
						str_append(&str, " %c0x%04X", val<0?'-':' ', abs(val));
					}
					str_append(&str, "\n");
					for(kc2=0;kc2<3;++kc2)
					{
						const int np=(int)(_countof(custom_params)/3);
						const int *params=custom_params+24*kc2;
						for(k=0;k<np;++k)
						{
							int x=k%stw;
							if(shift)
								str_append(&str, "%g,%c", (double)params[k]/65536, x<stw-1&&k<np-1?'\t':'\n');
							else
							{
								int val=params[k];
								str_append(&str, "%c0x%06X,%c", val<0?'-':' ', abs(val), x<stw-1&&k<np-1?'\t':'\n');
							}
						}
					}
				}
			}
			copy_to_clipboard((char*)HDR_DATA(str), (int)str->count);
			array_free(&str);
		}
		else/*no Ctrl*/
		{
			if(mode==VIS_HISTOGRAM)
			{
				int shift=GET_KEY_STATE(KEY_SHIFT);
				histmode+=1-(shift<<1);
				MODVAR(histmode, histmode, HISTMODE_COUNT);
			}
			else if(mode==VIS_MODEL)
			{
				modelvis^=1;
			}
			else
			{
				center_image();
				show_full_image=!show_full_image;
			}
		}
		return 1;
	case 'V':
		if(im1&&GET_KEY_STATE(KEY_CTRL))/*paste custom transform value*/
		{
			ArrayHandle text=paste_from_clipboard(0);
			if(text)
			{
				int k, kend, idx;

				idx=0;
				if(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4])
				{
					idx=parse_nvals_i32(text, idx, &ols4_period, 1);
					idx=parse_nvals_f64(text, idx, ols4_lr, (int)_countof(ols4_lr));
					idx=parse_nvals_i8(text, idx, (uint8_t*)ols4_mask, sizeof(ols4_mask)/sizeof(char));
				}
				if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
				{
					parse_nvals_i16(text, idx, (int16_t*)&c3_params, sizeof(c3_params)/sizeof(int16_t));
				}
				else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP])
				{
					parse_nvals_i16(text, idx, jxlparams_i16, (int)_countof(jxlparams_i16));
				}
				else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
				{
					parse_nvals_i16(text, idx, pw2_params, (int)_countof(pw2_params));
				}
				else if(transforms_customenabled)
				{
					int shift=GET_KEY_STATE(KEY_SHIFT);

					k=0, kend=(int)_countof(rct_custom_params);
					for(;k<kend;++k)
					{
						char *end=0;
						for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);
						rct_custom_params[k]=(int16_t)strtol((char*)HDR_DATA(text)+idx, &end, 16);
						for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx]);++idx);
						if(idx>=(int)text->count)
							goto paste_finish;
					}
					k=0;
					kend=(int)_countof(custom_params);
					if(shift)
					{
						for(;k<kend;++k)
						{
							for(;idx<(int)text->count&&isspace(HDR_DATC(text)[idx]);++idx);
							custom_params[k]=(int)round(atof((char*)HDR_DATA(text)+idx)*65536);
							for(;idx<(int)text->count&&!isspace(HDR_DATC(text)[idx]);++idx);
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
		if(GET_KEY_STATE(KEY_CTRL))/*paste bitmap*/
		{
			Image *im2=paste_bmp_from_clipboard();
			if(im2)
			{
				free(im0);
				im0=im2;
				filesize=((size_t)im2->src_depth[0]+im2->src_depth[1]+im2->src_depth[2]+im2->src_depth[3])*im2->iw*im2->ih>>3;
				set_window_title("From clipboard - eBench");
				update_image();
				if(imagecentered)
					center_image();
				return 1;
			}
		}
		else if(mode==VIS_JOINT_HISTOGRAM)
		{
			int shift=GET_KEY_STATE(KEY_SHIFT);
			space_not_color+=1-(shift<<1);
			MODVAR(space_not_color, space_not_color, 4);
			update_image();
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
					memcpy(lossyconv_params+5*5*lossyconv_page, lossyconv_params+5*5*srcidx, sizeof(int16_t[25]));
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
		if(GET_KEY_STATE(KEY_CTRL))
			pxlabels_hex=!pxlabels_hex;
		else
		{
			drawCRhistory^=1;
			//memset(combCRhist, 0, sizeof(combCRhist));
			//combCRhist_idx=0;
			//combCRhist_max=1;
		}
		return 1;
	case 'Q':
		return 1;
	case 'X':/*toggle horizontal profile plot*/
		if(profileplotmode!=PROFILE_X)
			profileplotmode=PROFILE_X;
		else
			profileplotmode=PROFILE_OFF;
		return 1;
	case 'Y':/*toggle vertical profile plot*/
		if(profileplotmode!=PROFILE_Y)
			profileplotmode=PROFILE_Y;
		else
			profileplotmode=PROFILE_OFF;
		return 1;
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
	case 'B':/*batch test*/
		if(GET_KEY_STATE(KEY_CTRL))
		{
			batch_test2(!GET_KEY_STATE(KEY_SHIFT));
		}
		else if(mode==VIS_L1WEIGHTS)
		{
			int fwd=!GET_KEY_STATE(KEY_SHIFT);
			fwd=2*fwd-1;
			mixmode=(mixmode+fwd+3)%3;
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
				int floorhalf, bit;

				if(abacvis_range<=2)
					break;
				floorhalf=abacvis_range>>1;
				bit=key==KEY_RBRACKET;
				if(bit)
					abacvis_low+=abacvis_range-floorhalf;
				abacvis_range=floorhalf;
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
			if(GET_KEY_STATE(KEY_SHIFT))/*change layer*/
			{
				lossyconv_page+=sign<<2;
				lossyconv_page&=15;
			}
			else/*change channel*/
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
		if(GET_KEY_STATE(KEY_CTRL))/*Ctrl N	add noise to params*/
		{
			if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				srand((uint32_t)__rdtsc());
				{
					int16_t *ptr;
					int k;

					ptr=(int16_t*)&c3_params;
					for(k=0;k<C3_NPARAMS;++k)
						ptr[k]+=rand()%3-1;
				}
				update_image();
			}
			return 1;
		}
		else
		{
			view_ma=!view_ma;
			update_image();
			return 1;
		}
		break;
	case KEY_SPACE:
		if(im1)
		{
			if(GET_KEY_STATE(KEY_CTRL))
			{
				if(!im1)
					return 0;
				if(mode==VIS_ANALYSIS)
				{
					char buf[512]={0};
					double esizes[3]={0};
					C18Info *info;
					intptr_t usize;
					int nblocks, kx, ky, nprinted;

					info=(C18Info*)malloc(sizeof(C18Info));
					if(!info)
					{
						WARNING("Alloc error");
						return 0;
					}
					nblocks=0;
					for(ky=0;ky<im1->ih;ky+=jhc_boxdy)
					{
						int ky2;

						ky2=ky+jhc_boxdy;
						if(ky2>im1->ih)
							ky2=im1->ih;
						for(kx=0;kx<im1->iw;kx+=jhc_boxdx)
						{
							const uint8_t *group;
							int kx2, selected_ch[3];

							kx2=kx+jhc_boxdx;
							if(kx2>im1->iw)
								kx2=im1->iw;
							c18_analyze(im1, kx, kx2, ky, ky2, info);
							group=rct2_combinations[info->bestrct];
							selected_ch[0]=group[0]*PRED_COUNT+info->predidx[0];
							selected_ch[1]=group[1]*PRED_COUNT+info->predidx[1];
							selected_ch[2]=group[2]*PRED_COUNT+info->predidx[2];
							esizes[0]+=info->esizes[selected_ch[0]]*(kx2-kx)*(ky2-ky);
							esizes[1]+=info->esizes[selected_ch[1]]*(kx2-kx)*(ky2-ky);
							esizes[2]+=info->esizes[selected_ch[2]]*(kx2-kx)*(ky2-ky);
							++nblocks;
						}
					}
					usize=(intptr_t)im1->iw*im1->ih;
					nprinted=acme_snprintf(buf, sizeof(buf)-1,
						"%d*%d block size, %d blocks\n"
						"T %6.2f%%  %12.2f bytes\n"
						"Y %6.2f%%  %12.2f bytes\n"
						"U %6.2f%%  %12.2f bytes\n"
						"V %6.2f%%  %12.2f bytes\n",
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
				}
				else
				{
					double rmse[5]={0}, psnr[5]={0};
					intptr_t idx;
					int printed, kx, ky;

					idx=calc_psnr(im0, im1, rmse, psnr);
					printed=acme_snprintf(g_buf, (size_t)G_BUF_SIZE-1,
						"RMSE PSNR\n"
						"T  %12.6f %12.6f\n"
						"C0 %12.6f %12.6f\n"
						"C1 %12.6f %12.6f\n"
						"C2 %12.6f %12.6f\n"
						"C3 %12.6f %12.6f"
						, rmse[4], psnr[4]
						, rmse[0], psnr[0]
						, rmse[1], psnr[1]
						, rmse[2], psnr[2]
						, rmse[3], psnr[3]
					);
					if(idx>=0)
					{
						kx=(int)(idx%im0->iw);
						ky=(int)(idx/im0->iw);
						printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed-1,
							"\nError vs Original at XY %d %d:\n"
							"C0\t0x%08X\t0x%08X\n"
							"C1\t0x%08X\t0x%08X\n"
							"C2\t0x%08X\t0x%08X\n"
							"C3\t0x%08X\t0x%08X",
							kx, ky,
							((int32_t*)HDR_DATA(im1))[idx<<2|0], ((int32_t*)HDR_DATA(im0))[idx<<2|0],
							((int32_t*)HDR_DATA(im1))[idx<<2|1], ((int32_t*)HDR_DATA(im0))[idx<<2|1],
							((int32_t*)HDR_DATA(im1))[idx<<2|2], ((int32_t*)HDR_DATA(im0))[idx<<2|2],
							((int32_t*)HDR_DATA(im1))[idx<<2|3], ((int32_t*)HDR_DATA(im0))[idx<<2|3]
						);
						copy_to_clipboard(g_buf, printed);
					}
					messagebox(MBOX_OK, idx>=0?"Error (Copied to clipboard)":"Bit Exact", "%s", g_buf);
				}
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
			else if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
			{
				int maskbits=3;
				Image *im2=0;

				if(keyboard[KEY_1])maskbits=1;
				else if(keyboard[KEY_2])maskbits=2;
				else if(keyboard[KEY_3])maskbits=3;
				else if(keyboard[KEY_4])maskbits=4;
				else if(keyboard[KEY_5])maskbits=5;
				else if(keyboard[KEY_6])maskbits=6;
				else if(keyboard[KEY_7])maskbits=7;
				else if(keyboard[KEY_8])maskbits=8;
				else if(keyboard[KEY_9])maskbits=9;
				image_copy(&im2, im0);
				if(!im2)
					return 0;
				apply_selected_transforms(&im2, 1, 1, 1);
				custom3_opt(im2, &c3_params, 0, maskbits, 1, 0);
				free(im2);
				update_image();
			}
			else if(transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP]||transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
			{
				int jxl=transforms_mask[ST_FWD_WP]||transforms_mask[ST_INV_WP];
				int pw2=transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM];
				if(GET_KEY_STATE(KEY_CTRL))
				{
					if(jxl)
					{
						int k;

						for(k=0;k<33;++k)
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
	default:break;
	}
	return 0;
}
int io_keyup(IOKey key, char c)
{
	(void)c;

	switch(key)
	{
	case KEY_SPACE:
#define AK(KEY) case KEY:
	ACTIVE_KEY_LIST
#undef  AK
		count_active_keys(key);
		break;

	case 'Q':
		return mode==VIS_IMAGE;
	default:break;
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
		const int delta=10;/*screen pixels per frame*/
		if(keyboard[KEY_ENTER])	zoom_at(wndw>>1, wndh>>1, 1.02);
		if(keyboard[KEY_BKSP])	zoom_at(wndw>>1, wndh>>1, 1/1.02);
		if(keyboard['W'])	wpy-=delta/imzoom;/*move window up	*/
		if(keyboard['A'])	wpx-=delta/imzoom;/*move window left	*/
		if(keyboard['S'])	wpy+=delta/imzoom;/*move window down	*/
		if(keyboard['D'])	wpx+=delta/imzoom;/*move window right	*/
	}
}
static float print_i16_row(float x, float y, float zoom, const int16_t *row, int count)
{
	int printed, k;

	for(printed=k=0;k<count;++k)
	{
		int16_t val=row[k];
		printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, " %c%4X", val<0?'-':' ', abs(val));
	}
	return print_line_immediate(0, x, y, zoom, g_buf, printed, -1, 0, 0);
}

static ArrayHandle vertices_2d=0;
static void draw_profile_x(int comp, int color)/*horizontal cross-section profile		to see the color/spatial correlation*/
{
	int iy, *row, ix, offset, kx;
	float y2, gain;

	if((float)wndh-tdy-1<=0)
		return;
	iy=screen2image_y_int(my);
	if((uint32_t)iy<(uint32_t)im1->ih)
	{
		row=((int32_t*)HDR_DATA(im1))+(im1->iw*iy<<2|comp);
		gain=((float)wndh-tdy-1)/((1<<im1->depth[comp])-1);
		offset=1<<im1->depth[comp]>>1;
		for(kx=0;kx<wndw;++kx)
		{
			ix=screen2image_x_int(kx);
			if((uint32_t)ix<(uint32_t)im1->iw)
				y2=wndh-tdy-(row[ix<<2]+offset)*gain;
			else
				y2=wndh-tdy;
			draw_curve_enqueue(&vertices_2d, (float)kx, y2);
		}
		draw_2d_flush(vertices_2d, color, GL_LINE_STRIP);
	}
}
static void draw_profile_y(int comp, int color)/*vertical cross-section profile*/
{
	int ix, *col, iy, offset, stride, ky;
	float x2, gain;

	ix=screen2image_x_int(mx);
	if((uint32_t)ix<(uint32_t)im1->iw)
	{
		col=((int32_t*)HDR_DATA(im1))+(ix<<2|comp);
		
		gain=((float)wndw-1)/((1<<im1->depth[comp])-1);
		offset=1<<im1->depth[comp]>>1;
		stride=im1->iw<<2;
		for(ky=0;ky<wndh;++ky)
		{
			iy=screen2image_y_int(ky);
			if((uint32_t)iy<(uint32_t)im1->ih)
				x2=(col[iy*stride]+offset)*gain;
			else
				x2=0;
			draw_curve_enqueue(&vertices_2d, x2, (float)ky);
		}
		draw_2d_flush(vertices_2d, color, GL_LINE_STRIP);
	}
}
static void draw_profile_xy(int comp1, int comp2, int color)/*horizontal cross-section profile		to see the color/spatial correlation*/
{
	int iy, *row1, *row2, ix, k;

	iy=screen2image_y_int(my);
	if((uint32_t)iy<(uint32_t)im1->ih)
	{
		row1=((int32_t*)HDR_DATA(im1))+(im1->iw*iy<<2|comp1);
		row2=((int32_t*)HDR_DATA(im1))+(im1->iw*iy<<2|comp2);
		ix=screen2image_x_int(mx);
		for(k=ix-10;k<=ix+10;++k)
		{
			if((uint32_t)k<(uint32_t)im1->iw)
				draw_curve_enqueue(&vertices_2d, wndw*0.5f+2*(float)row1[k<<2], wndh*0.5f+2*(float)row2[k<<2]);
		}
		draw_2d_flush(vertices_2d, color, GL_LINE_STRIP);
	}
}

static ArrayHandle vertices_text=0;
static void print_pixellabels(int ix1, int ix2, int iy1, int iy2, int component, char label, int64_t txtcolors, int depth)
{
	const char *format;
	int nlevels, half, mask, ix, kx, xend, iy, ky, yend, idx;
	float fontsize, labeloffset;

	nlevels=1<<depth;
	half=nlevels>>1;
	mask=nlevels-1;
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
	iy=MAXVAR(iy1, 0);
	fontsize=1;
	labeloffset=tdy*fontsize*component;
	for(yend=MINVAR(iy2+2, im1->ih);iy<yend;++iy)
	{
		ky=image2screen_y_int(iy);
		ix=MAXVAR(ix1, 0);
		for(xend=MINVAR(ix2+2, im1->iw);ix<xend;++ix)
		{
			kx=image2screen_x_int(ix);
			idx=(im1->iw*iy+ix)<<2;
			GUIPrint_enqueue(&vertices_text, 0, (float)kx, (float)ky+labeloffset, fontsize, format, label, (((int32_t*)HDR_DATA(im1))[idx+component]+half)&mask);
		}
	}
	print_line_flush(vertices_text, fontsize);
	txtcolors=set_text_colors(txtcolors);
}
void io_render(void)
{
	double usize, csize;

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	if(!wndh)
		return;
	{
		const float axes[]=
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
	if(im1)
	{
		int
			sx1=image2screen_x_int(0), sx2=image2screen_x_int(im1->iw),
			sy1=image2screen_y_int(0), sy2=image2screen_y_int(im1->ih);
		switch(mode)
		{
		case VIS_JOINT_HISTOGRAM:
			chart_jointhist_draw();
			break;
		case VIS_ANALYSIS:
			{
				const uint8_t *group;
				int selected_ch[3];
				int bounds[4]={0};/*{x1, x2, y1, y2}*/
				float xmid, ymid, y, centers[6]={0};
				int ncenters, k;

				group=rct2_combinations[analysis_info.bestrct];
				selected_ch[0]=group[0]*PRED_COUNT+analysis_info.predidx[0];
				selected_ch[1]=group[1]*PRED_COUNT+analysis_info.predidx[1];
				selected_ch[2]=group[2]*PRED_COUNT+analysis_info.predidx[2];
				jhc_getboxbounds(jhc_xbox, jhc_ybox, jhc_boxdx, jhc_boxdy, im1->iw, im1->ih, bounds);
				
				display_texture_i(0, wndw, 0, wndh, (int*)im_export, im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);

				{
					const int crosshaircolor=0xFF000000;
					float ratios[2], sbounds[4];

					ratios[0]=(float)wndw/im1->iw;
					ratios[1]=(float)wndh/im1->ih;
					for(k=0;k<4;++k)
						sbounds[k]=bounds[k]*ratios[k>>1];
					xmid=(sbounds[0]+sbounds[1])*0.5f;
					ymid=(sbounds[2]+sbounds[3])*0.5f;
					draw_rect_hollow(sbounds[0]+1, sbounds[1]+1, sbounds[2]+1, sbounds[3]+1, 0xFFFFFFFF);
					draw_rect_hollow(sbounds[0], sbounds[1], sbounds[2], sbounds[3], crosshaircolor);
					draw_line(xmid, 0, xmid, (float)wndh, crosshaircolor);
					draw_line(0, ymid, (float)wndw, ymid, crosshaircolor);
				}
				y=tdy*2;
				GUIPrint(0, (float)(0.25*wndw), y-tdy*1.5f-2, 1.5f,
					"3*%d*%d: %6.2f%% + %6.2f%% + %6.2f%% = %6.2f%%  %s %s %s %s  %f MB/s",
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
				y+=10;
				ncenters=0;
				for(k=0;k<RCT2_COUNT;++k)/*RCTs*/
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
				for(k=0;k<OCH2_COUNT*PRED_COUNT;++k)/*output channels*/
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
				}
				for(k=0;k<PRED_COUNT;++k)
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
		case VIS_SHANNON:
			{
				if(mode==VIS_IMAGE&&GET_KEY_STATE('Q'))
				{
					int mx2, xstart, xend;
					double divx;

					mx2=mx;
					xstart=image2screen_x_int(0);
					xend=image2screen_x_int(im0->iw);
					CLAMP2(mx2, xstart, xend);
					divx=screen2image_x(mx2)/im0->iw;
					CLAMP2(divx, 0, 1);
					if(divx>0)
						display_texture_i(sx1, mx2, sy1, sy2, (int*)zimage, im1->iw, im1->ih, 0, (float)divx, 0, 1, 1, 0);
					if(divx<1)
						display_texture_i(mx2, sx2, sy1, sy2, (int*)im_export, im1->iw, im1->ih, (float)divx, 1, 0, 1, 1, 0);
				}
				else
					display_texture_i(sx1, sx2, sy1, sy2, (int*)(mode==VIS_SHANNON?zimage:im_export), im1->iw, im1->ih, 0, 1, 0, 1, 1, 0);
				if(imzoom>=ZOOM_LIMIT_LABEL)
				{
					const int64_t theme[]=
					{
						(uint64_t)0xC00000FF<<32|0x80000000,
						(uint64_t)0xC000FF00<<32|0x80000000,
						(uint64_t)0xC0FF0000<<32|0x80FFFFFF,
						(uint64_t)0xC0FFFFFF<<32|0x80000000,
					};
					int csx1, csx2, csy1, csy2, ix1, ix2, iy1, iy2;

					csx1=CLAMP(0, sx1, wndw);
					csx2=CLAMP(0, sx2, wndw);
					csy1=CLAMP(0, sy1, wndh);
					csy2=CLAMP(0, sy2, wndh);
					ix1=screen2image_x_int(csx1);
					ix2=screen2image_x_int(csx2);
					iy1=screen2image_y_int(csy1);
					iy2=screen2image_y_int(csy2);
					if(viewmode==VIEW_RGB||viewmode==VIEW_C0)
						print_pixellabels(ix1, ix2, iy1, iy2, 0, 'r', theme[0], im1->depth[0]);
					if(viewmode==VIEW_RGB||viewmode==VIEW_C1)
						print_pixellabels(ix1, ix2, iy1, iy2, 1, 'g', theme[1], im1->depth[1]);
					if(viewmode==VIEW_RGB||viewmode==VIEW_C2)
						print_pixellabels(ix1, ix2, iy1, iy2, 2, 'b', theme[2], im1->depth[2]);
					if(imzoom>=ZOOM_LIMIT_ALPHA&&im1->depth[3])
						print_pixellabels(ix1, ix2, iy1, iy2, 3, 'a', theme[3], im1->depth[3]);
					{
						double kx=screen2image_x(wndw*0.5f);
						double ky=screen2image_y(wndh*0.5f);
						GUIPrint(0, (float)wndw/2, (float)wndh/2, 1, "(%.2f, %.2f)", kx, ky);
					}
				}
#ifdef ENABLE_L1WEIGHTS
				if(mode==VIS_L1WEIGHTS)
				{
					enum
					{
						SIGNAL_WEIGHT_START,
						SIGNAL_BIAS=SIGNAL_WEIGHT_START+L1NPREDS,
						SIGNAL_ESTIM_START,
						SIGNAL_CURR=SIGNAL_ESTIM_START+L1NPREDS,
						SIGNAL_PRED_L1,
						SIGNAL_PRED_WP,
						SIGNAL_PRED_OLS,
						SIGNAL_COUNT,

						L1HOLD=1,
						L1XSCALE=16,
						L1HISTSIZE=128,
						L1SMOOTH=16,
						L1YSCALE=4	/*256*4 = 1024 pixels*/
					};
					static int posctr=0;/*image position*/
					static int gpos=0;/*graph position*/
					static int whist[L1HISTSIZE][SIGNAL_COUNT]={0};/*{weights, bias, estimators, pred, curr, delta}*/
					static int ctr=0;
					static int localhist[3][256]={0};
					static int64_t mad[3]={0}, rmse[3]={0};
					static int prevmx=0, prevmy=0;
					static int preds[L1NPREDS]={0};
					int posctr0, kx, ky, speed, smooth, posidx0, posidx, *currw, kpos;

					posctr0=0;
					kx=0;
					ky=0;
					speed=GET_KEY_STATE(KEY_SHIFT)?20*L1SMOOTH*L1SPEED:L1SPEED;
					if(GET_KEY_STATE('L'))
					{
						kx=screen2image_x_int(mx);
						ky=screen2image_y_int(my);
						CLAMP2(kx, 0, im1->iw-1);
						CLAMP2(ky, 0, im1->ih-1);
						posctr0=posctr=L1SMOOTH*(im1->iw*ky+kx);
						if(prevmx!=mx||prevmy!=my)
						{
							speed=L1SMOOTH;
							prevmx=mx;
							prevmy=my;
						}
						else
							speed=0;
					}
					posctr0=posctr;
					if(!GET_KEY_STATE('P'))
						posctr+=speed;
					smooth=(posctr+L1SMOOTH)%L1SMOOTH;
					posidx0=posctr0/L1SMOOTH;
					posidx=posctr/L1SMOOTH;
					if(!GET_KEY_STATE('L'))/*pan to cursor (wpx, wpy)*/
					{
						int x1=screen2image_x_int(wndw*0/8);
						int x2=screen2image_x_int(wndw*7/8);
						int y1=screen2image_y_int(wndh*0/8);
						int y2=screen2image_y_int(wndh*7/8);
						int ix=posidx%im1->iw;
						int iy=posidx/im1->iw;
						if(ix<x1)wpx+=ix-x1;
						if(ix>x2)wpx+=ix-x2;
						if(iy<y1)wpy+=iy-y1;
						if(iy>y2)wpy+=iy-y2;
					}
					currw=weights_l1;
					if(mixmode==1)
						currw=weights_wp;
					else if(mixmode==2)
						currw=weights_ols;
					for(kpos=posidx0;kpos<posidx;++kpos)
					{
						int NNN, NN, NNE, NNEEE, NWW, NW, N, NE, NEEE, NEEEEE, WWWWW, WWW, WW, W, curr, j2, k;
						int pred_l1, pred_l10, pred_wp, pred_ols;
						int vmax, vmin;

						kx=kpos%im1->iw;
						ky=kpos/im1->iw;
						if(GET_KEY_STATE('L'))
						{
							kx=screen2image_x_int(mx);
							ky=screen2image_y_int(my);
						}
						CLAMP2(kx, 0, im1->iw-1);
						CLAMP2(ky, 0, im1->ih-1);
						NNN	=(uint32_t)(ky-3)<(uint32_t)im1->ih&&(uint32_t)(kx+0)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-3)*im1->iw+kx+0)*4]:0;
						NN	=(uint32_t)(ky-2)<(uint32_t)im1->ih&&(uint32_t)(kx+0)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-2)*im1->iw+kx+0)*4]:0;
						NNE	=(uint32_t)(ky-2)<(uint32_t)im1->ih&&(uint32_t)(kx+1)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-2)*im1->iw+kx+1)*4]:0;
						NNEEE	=(uint32_t)(ky-2)<(uint32_t)im1->ih&&(uint32_t)(kx+3)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-2)*im1->iw+kx+3)*4]:0;
						NWW	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx-2)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx-2)*4]:0;
						NW	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx-1)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx-1)*4]:0;
						N	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx+0)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx+0)*4]:0;
						NE	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx+1)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx+1)*4]:0;
						NEEE	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx+3)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx+3)*4]:0;
						NEEEEE	=(uint32_t)(ky-1)<(uint32_t)im1->ih&&(uint32_t)(kx+5)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-1)*im1->iw+kx+5)*4]:0;
						WWWWW	=(uint32_t)(ky-0)<(uint32_t)im1->ih&&(uint32_t)(kx-5)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-0)*im1->iw+kx-5)*4]:0;
						WWW	=(uint32_t)(ky-0)<(uint32_t)im1->ih&&(uint32_t)(kx-3)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-0)*im1->iw+kx-3)*4]:0;
						WW	=(uint32_t)(ky-0)<(uint32_t)im1->ih&&(uint32_t)(kx-2)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-0)*im1->iw+kx-2)*4]:0;
						W	=(uint32_t)(ky-0)<(uint32_t)im1->ih&&(uint32_t)(kx-1)<(uint32_t)im1->iw?((int32_t*)HDR_DATA(im1))[((ky-0)*im1->iw+kx-1)*4]:0;
						curr	=((int32_t*)HDR_DATA(im1))[(ky*im1->iw+kx)*4];
						(void)NNN;
						(void)NN;
						(void)NNE;
						(void)NNEEE;
						(void)NWW;
						(void)NW;
						(void)N;
						(void)NE;
						(void)NEEE;
						(void)NEEEEE;
						(void)WWWWW;
						(void)WWW;
						(void)WW;
						(void)W;
#define L1PRED(W0, EXPR) preds[j2++]=EXPR;
						j2=0;
						L1PREDLIST;
#undef  L1PRED
						vmax=N; vmin=W;
						if(N<W)vmin=N, vmax=W;
						CLAMP2(preds[5], vmin, vmax);/*CG = median(N, W, N+W-NW)*/
						{
							enum
							{
								NUMBITS=15,
								DENBITS=7,
								DIVLUTSIZE=(1<<DENBITS)
							};
							static int divlookup[DIVLUTSIZE]={0};
							uint32_t wsum;
							int ipred, shift, kp;

							if(!*divlookup)
							{
								for(k=0;k<DIVLUTSIZE;++k)
									divlookup[k]=(1<<NUMBITS)/(k+1);
							}
							wsum=0;
							ipred=0;
							shift=0;
							for(kp=0;kp<L1NPREDS;++kp)
							{
								int e=wperrors[kp];
								shift=31^LZCNT32(e+1);
								weights_wp[kp]=(divlookup[e<<(DENBITS-1)>>shift]<<(DENBITS-1)>>shift)+(1<<DENBITS>>2)/L1NPREDS;
								wsum+=weights_wp[kp];
							}
							shift=(31-(DENBITS-2))-LZCNT32(wsum);
							wsum=0;
							for(kp=0;kp<L1NPREDS;++kp)
							{
								int c=(int)(weights_wp[kp]>>shift);
								wsum+=c;
								weights_wp[kp]=c;
							}
							ipred=wsum>>1;
							for(kp=0;kp<L1NPREDS;++kp)
							{
								ipred+=(int)(weights_wp[kp]*preds[kp]);
								weights_wp[kp]<<=8;
							}
							pred_wp=ipred*divlookup[wsum-1]>>NUMBITS;
						}
						{
							double fpred;

							fpred=0;
							for(k=0;k<L1NPREDS;++k)
							{
								fpred+=curr_params[k]*preds[k];
								weights_ols[k]=(int)round(curr_params[k]*0x20000);
							}
							pred_ols=(int)ROUND64(fpred);
						}
						{
							pred_l1=weights_l1[L1NPREDS];
							for(k=0;k<L1NPREDS;++k)
								pred_l1+=weights_l1[k]*preds[k];
							pred_l1+=1<<L1SH>>1;
							pred_l1>>=L1SH;
						}
						if(vmin>NE)vmin=NE;
						if(vmax<NE)vmax=NE;
						if(vmin>NEEE)vmin=NEEE;
						if(vmax<NEEE)vmax=NEEE;
						pred_l10=pred_l1;
						CLAMP2(pred_l1, vmin, vmax);
						CLAMP2(pred_wp, vmin, vmax);
						CLAMP2(pred_ols, vmin, vmax);

						gpos=(gpos+1)%L1HISTSIZE;
						for(k=0;k<L1NPREDS+1;++k)
							whist[gpos][SIGNAL_WEIGHT_START+k]=currw[k];
						for(k=0;k<L1NPREDS;++k)
							whist[gpos][SIGNAL_ESTIM_START+k]=preds[k];
						whist[gpos][SIGNAL_PRED_L1	]=pred_l1;
						whist[gpos][SIGNAL_PRED_WP	]=pred_wp;
						whist[gpos][SIGNAL_PRED_OLS	]=pred_ols;
						whist[gpos][SIGNAL_CURR		]=curr;
						{
							int delta;

							delta=curr-pred_l1;
							mad[0]+=abs(delta); rmse[0]+=delta*delta;
							delta=curr-pred_wp;
							mad[1]+=abs(delta); rmse[1]+=delta*delta;
							delta=curr-pred_ols;
							mad[2]+=abs(delta); rmse[2]+=delta*delta;
						}
						++localhist[0][(curr-pred_l1+128)&255];
						++localhist[1][(curr-pred_wp+128)&255];
						++localhist[2][(curr-pred_ols+128)&255];
						++ctr;

						/*update*/
						{
							int e2[L1NPREDS], best=0x7FFFFFFF;
							int kp;

							for(kp=0;kp<L1NPREDS;++kp)
							{
								int e=abs(curr-preds[kp]);
								if(best>e)
									best=e;
								e2[kp]=e;
							}
							for(kp=0;kp<L1NPREDS;++kp)
							{
								int e=e2[kp]-best;
								wperrors[kp]+=((e<<6)-wperrors[kp]+(1<<3>>1))>>3;
							}
						}
						{
							const double lr=0.00001;
							double lval, lr_comp, sum;
							int ky2, kx2, midx, i, j, success;

							for(ky2=0, midx=0;ky2<L1NPREDS;++ky2)
							{
								for(kx2=0;kx2<L1NPREDS;++kx2, ++midx)
									curr_cov[midx]+=(preds[kx2]*preds[ky2]-curr_cov[midx])*lr;
							}
							lval=curr*lr;
							lr_comp=1-lr;
							for(k=0;k<L1NPREDS;++k)
								curr_vec[k]=lval*preds[k]+lr_comp*curr_vec[k];
							success=1;
							memcpy(curr_cholesky, curr_cov, sizeof(double[L1NPREDS*L1NPREDS]));
							for(k=0;k<L1NPREDS*L1NPREDS;k+=L1NPREDS+1)
								curr_cholesky[k]+=0.0075;
							for(i=0;i<L1NPREDS;++i)
							{
								for(j=0;j<i;++j)
								{
									sum=curr_cholesky[i*L1NPREDS+j];
									for(k=0;k<j;++k)
										sum-=curr_cholesky[i*L1NPREDS+k]*curr_cholesky[j*L1NPREDS+k];
									curr_cholesky[i*L1NPREDS+j]=sum/curr_cholesky[j*L1NPREDS+j];
								}
								sum=curr_cholesky[i*L1NPREDS+i];
								for(k=0;k<i;++k)
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
								for(i=0;i<L1NPREDS;++i)
								{
									sum=curr_vec[i];
									for(j=0;j<i;++j)
										sum-=curr_cholesky[i*L1NPREDS+j]*curr_params[j];
									curr_params[i]=sum/curr_cholesky[i*L1NPREDS+i];
								}
								for(i=L1NPREDS-1;i>=0;--i)
								{
									sum=curr_params[i];
									for(j=i+1;j<L1NPREDS;++j)
										sum-=curr_cholesky[j*L1NPREDS+i]*curr_params[j];
									curr_params[i]=sum/curr_cholesky[i*L1NPREDS+i];
								}
							}
						}
						{
							int e=(curr>pred_l10)-(curr<pred_l10);
							weights_l1[L1NPREDS]+=e;
							for(k=0;k<L1NPREDS;++k)
								weights_l1[k]+=e*preds[k];
						}
					}

					/*draw*/
					if(wndh)
					{
						static ArrayHandle vertices=0;
						const float charty[]={200, 400, 600};
						const int colors_estim[]=
						{
							0x80FF4080,
							0x8080FF40,
							0x804080FF,
							0x80FFC080,
							0x8080FFC0,
							0x80C080FF,
						};
						const int colors_pred[]=
						{
							0xFF0000FF,
							0xFF00FF00,
							0xFFFF0000,
						};
						const int colors_delta[]=
						{
							0xFF8080FF,
							0xFF80FF80,
							0xFFFF8080,
						};
						const int ncolors=_countof(colors_estim);
						int dlen, color, gheight_signal, gheight_delta, k;

						dlen=wndw;
						if(dlen>L1HISTSIZE)
							dlen=L1HISTSIZE;
						color=0;
						/*draw mix L1*/
						{
							int sum=1<<L1SH;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_l1[k]/sum;
								draw_rect(xs0, xs, charty[0]+50, charty[0]+50-preds[k], colors_estim[k%ncolors]);
							}
						}

						/*draw mix WP*/
						{
							int sum=0;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
								sum+=weights_wp[k];
							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_wp[k]/sum;
								draw_rect(xs0, xs, charty[1]+50, charty[1]+50-preds[k], colors_estim[k%ncolors]);
							}
						}

						/*draw mix OLS*/
						{
							int sum=0x20000;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_ols[k]/sum;
								draw_rect(xs0, xs, charty[2]+50, charty[2]+50-preds[k], colors_estim[k%ncolors]);
							}
						}
						
						/*label mix L1*/
						{
							int sum=1<<L1SH;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_l1[k]/sum;
								GUIPrint(0, xs0, charty[0]+k*tdy, 1, "%s", l1prednames[k]);
							}
						}

						/*label mix WP*/
						{
							int sum=0;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
								sum+=weights_wp[k];
							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_wp[k]/sum;
								GUIPrint(0, xs0, charty[1]+k*tdy, 1, "%s", l1prednames[k]);
							}
						}

						/*label mix OLS*/
						{
							int sum=0x20000;
							float xs=0;

							for(k=0;k<L1NPREDS;++k)
							{
								float xs0=xs;
								xs+=1000.f*weights_ols[k]/sum;
								GUIPrint(0, xs0, charty[2]+k*tdy, 1, "%s", l1prednames[k]);
							}
						}
						{
							int xs, kp;

							gheight_signal=-128;
							gheight_delta=128;

							/*curr - black bars*/
							color=0x80C0C0C0;
							for(xs=1;xs<dlen;++xs)
							{
								int kx2, c0;
								float ys;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								ys=(float)(wndh*7/8-512-gheight_signal-whist[kx2][SIGNAL_CURR]*L1YSCALE);
								draw_rect_enqueue(&vertices
									, L1XSCALE*((float)xs-0.4f-(float)smooth/L1SMOOTH)
									, L1XSCALE*((float)xs+0.4f-(float)smooth/L1SMOOTH)
									, ys
									, (float)(wndh*7/8-512-gheight_signal-0*L1YSCALE)
								);
								if(xs==dlen-1)
								{
									c0=set_text_color(color);
									GUIPrint(0, L1XSCALE*(float)xs+10, ys, 1, "%4d", whist[kx2][SIGNAL_CURR]);
									set_text_color(c0);
								}
							}
							draw_2d_flush(vertices, color, GL_TRIANGLES);

							/*pred L1 - red*/
							color=colors_pred[0];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_signal-whist[kx2][SIGNAL_PRED_L1]*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);
							
							/*pred WP - green*/
							color=colors_pred[1];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_signal-whist[kx2][SIGNAL_PRED_WP]*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);
							
							/*pred OLS - blue*/
							color=colors_pred[2];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_signal-whist[kx2][SIGNAL_PRED_OLS]*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);

							/*estim - orange*/
							for(kp=1;kp<MINVAR(10, L1NPREDS);++kp)
							{
								if(GET_KEY_STATE('0'+kp)&&kp-1<L1NPREDS)
								{
									int c0;

									color=colors_estim[(kp-1)%ncolors];
									for(xs=1;xs<dlen;++xs)/*estim*/
									{
										int kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
										draw_curve_enqueue(&vertices/*estim[kp-1]*/
											, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_signal-whist[kx2][SIGNAL_ESTIM_START+kp-1]*L1YSCALE)
										);
									}
									draw_2d_flush(vertices, color, GL_LINE_STRIP);
									for(xs=1;xs<dlen;++xs)/*estim delta*/
									{
										int kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
										draw_curve_enqueue(&vertices/*estim[kp-1]*/
											, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_delta-(whist[kx2][SIGNAL_CURR]-whist[kx2][SIGNAL_ESTIM_START+kp-1])*L1YSCALE)
										);
									}
									draw_2d_flush(vertices, color, GL_LINE_STRIP);
									c0=set_text_color(color);
									GUIPrint(0, wndw/2.f, wndh/2.f, 4, "%s", kp-1<L1NPREDS+1-1?l1prednames[kp-1]:"bias");
									set_text_color(c0);
								}
							}


							/*delta L1 - gray red*/
							color=colors_delta[0];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_delta-(whist[kx2][SIGNAL_CURR]-whist[kx2][SIGNAL_PRED_L1])*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);

							/*delta WP - gray green*/
							color=colors_delta[1];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_delta-(whist[kx2][SIGNAL_CURR]-whist[kx2][SIGNAL_PRED_WP])*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);

							/*delta OLS - gray blue*/
							color=colors_delta[2];
							for(xs=1;xs<dlen;++xs)
							{
								int kx2;

								kx2=(xs+gpos+L1HISTSIZE)%L1HISTSIZE;
								draw_curve_enqueue(&vertices
									, L1XSCALE*((float)xs-(float)smooth/L1SMOOTH), (float)(wndh*7/8-512-gheight_delta-(whist[kx2][SIGNAL_CURR]-whist[kx2][SIGNAL_PRED_OLS])*L1YSCALE)
								);
							}
							draw_2d_flush(vertices, color, GL_LINE_STRIP);
						}
						{
							const char *labels[]=
							{
								"L1",
								"WP",
								"OLS",
							};
							double e[3]={0};
							int km, ks, c0;

							for(km=0;km<3;++km)
							{
								for(ks=0;ks<256;++ks)
								{
									int freq=localhist[km][ks];
									if(freq)
									{
										double p=(double)freq/ctr;
										e[km]-=p*log(p);
									}
								}
								e[km]*=1./(8*M_LN2);
								c0=set_text_color(colors_pred[km]);
								GUIPrint(0, 0, charty[km]-2*tdy, 2, "%-3s MAD %7.2f RMSE %7.2f H%8d E%10.6f%%"
									, labels[km]
									, (double)mad[km]/ctr
									, sqrt((double)rmse[km]/ctr)
									, ctr
									, e[km]*100
								);
								set_text_color(c0);
							}
							GUIPrint(0, wndw*0.5f, wndh*0.5f, 2, "X%5d Y%5d"
								, posidx%im1->iw
								, posidx/im1->iw
							);
							{
								float xs, ys;

								kx=posidx%im1->iw;
								ky=posidx/im1->iw;
								if(GET_KEY_STATE('L'))
								{
									kx=screen2image_x_int(mx);
									ky=screen2image_y_int(my);
								}
								CLAMP2(kx, 0, im1->iw-1);
								CLAMP2(ky, 0, im1->ih-1);
								xs=(float)image2screen_x(kx+0.5f);
								ys=(float)image2screen_y(ky+0.5f);
								draw_rect_hollow(xs-10, xs+10, ys-10, ys+10, 0xFFFF00FF);
								draw_line(wndw*0.5f, wndh*0.5f, xs, ys, 0xFFFF00FF);
							}
						}
					}
				}
				else
#endif
				if(mode==VIS_HISTOGRAM)
				{
					if(histmode==HISTMODE_LINEAR_NOSPIKE||histmode==HISTMODE_LOG2_NOSPIKE)
					{
						if(wndh*0.75f>5*tdy)
							chart_hist_draw(0, (float)wndw, (float)wndh*0.25f, (float)wndh-5*tdy
								, 0, 3, 0, 0x60, hist, histmax
							);
					}
					else
					{
						if(wndh>5*tdy)
							chart_hist_draw(0, (float)wndw, 0, (float)wndh-5*tdy
								, 0, 3, 0, 0x60, hist, histmax
							);
					}
					GUIPrint(0, 0, 5*tdy, 1, "Mean, RMSE, MAD, RMSE/MAD, Gsize, Lsize");
					GUIPrint(0, 0, 6*tdy, 1, "Y %10.4f %10.4f %10.4f %10.4f %12.2f %12.2f", hist_mean[0]-128, hist_sdev[0], hist_mad[0], hist_sdev[0]/hist_mad[0], hist_gausslike[0], hist_laplacelike[0]);
					GUIPrint(0, 0, 7*tdy, 1, "U %10.4f %10.4f %10.4f %10.4f %12.2f %12.2f", hist_mean[1]-128, hist_sdev[1], hist_mad[1], hist_sdev[1]/hist_mad[1], hist_gausslike[1], hist_laplacelike[1]);
					GUIPrint(0, 0, 8*tdy, 1, "V %10.4f %10.4f %10.4f %10.4f %12.2f %12.2f", hist_mean[2]-128, hist_sdev[2], hist_mad[2], hist_sdev[2]/hist_mad[2], hist_gausslike[2], hist_laplacelike[2]);
				}
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
					double msizes, esizes[4]={0}, esize;
					int nhist, YUV, k, wndh2, xhist, yhist, smax, amax, nlevels, kh, kc, ks;
					float x4;

					nhist=modelnch*modelnctx;
					YUV=0;
					for(k=0;k<CST_COMPARE;++k)
					{
						if(transforms_mask[k])
						{
							YUV=1;
							break;
						}
					}
					wndh2=wndh-(int)(tdy*5);
					xhist=modelnch;
					yhist=nhist/modelnch;
					msizes=0;
					smax=0;
					amax=0;
					nlevels=1<<modeldepth;
					for(kh=0;kh<nhist;++kh)
					{
						int *curr_hist, hsum;

						curr_hist=modelhist+((ptrdiff_t)kh<<modeldepth);
						hsum=0;
						for(ks=0;ks<nlevels;++ks)
						{
							int f=curr_hist[ks];
							if(amax<f)amax=f;
							hsum+=f;
						}
						if(!kh||smax<hsum)
							smax=hsum;
					}
					for(kh=0;kh<nhist;++kh)
					{
						int hx, hy, xstart, xend, ystart, yend, *curr_hist, fmax, hsum/*, *curr_hist2, fmax2*/;

						hx=kh/yhist;
						hy=kh%yhist;
						xstart=hx*wndw/xhist;
						ystart=hy*wndh2/yhist;
						xend=(hx+1)*wndw/xhist;
						yend=(hy+1)*wndh2/yhist;
						kc=kh/modelnctx;
						curr_hist=modelhist+((ptrdiff_t)kh<<modeldepth);
						fmax=0;
						hsum=0;
						for(ks=0;ks<nlevels;++ks)/*FIXME*/
						{
							int f=curr_hist[ks];
							if(fmax<f)fmax=f;
							hsum+=f;
						}
						if(fmax)
						{
							if(modelvis)
							{
								draw_rect((float)xstart
									, (float)xstart+(((float)xend-xstart)*hsum*0.9f/smax)
									, (float)ystart+(yend-ystart)*0.25f
									, (float)ystart+(yend-ystart)*0.75f
									, 0x80000000|128<<kc*8
								);
								for(ks=0;ks<nlevels;++ks)
								{
									int f;

									f=curr_hist[ks];
									if(f)
										draw_rect(
											(float)xstart+(xend-xstart)*ks/nlevels,
											(float)xstart+(xend-xstart)*(ks+1)/nlevels,
											(float)yend-ceilf((float)(yend-ystart)*f/fmax),
											(float)yend,
											modeltheme[YUV<<2|kc]
										);
								}
							}
							else
							{
								for(ks=0;ks<nlevels;++ks)
								{
									int f;

									f=curr_hist[ks];
									if(f)
									{
										int c=modeltheme[YUV<<2|kc];
#if 0
										if((ks^kh)&1)
											c^=0x001F1F1F;
#endif
#if 0
										{
											uint8_t *p=(uint8_t*)&c;
											int t;
											
											t=p[0]+(rand()>>(15-5))-16;
											CLAMP2(t, 0, 255);
											p[0]=t;
											t=p[1]+(rand()>>(15-5))-16;
											CLAMP2(t, 0, 255);
											p[1]=t;
											t=p[2]+(rand()>>(15-5))-16;
											CLAMP2(t, 0, 255);
											p[2]=t;
										}
#endif
										draw_rect(
											(float)xstart+(xend-xstart)*ks/nlevels,
											(float)xstart+(xend-xstart)*(ks+1)/nlevels,
											(float)yend-ceilf((float)(yend-ystart)*f/amax),
											(float)yend,
											c
										);
									}
								}
							}
						}
#if 0
						curr_hist2=modelmhist+((ptrdiff_t)kh<<modeldepth);
						fmax2=0;
						for(ks=0;ks<nlevels;++ks)/*FIXME*/
						{
							int f=curr_hist2[ks];
							if(fmax2<f)fmax2=f;
						}
						if(fmax2)
						{
							for(ks=0;ks<nlevels-1;++ks)
							{
								int freq1=curr_hist2[ks+0];
								int freq2=curr_hist2[ks+1];
								draw_line((float)xstart+(xend-xstart)*(2*ks+1)/(2*nlevels)
									, (float)yend-(yend-ystart)*freq1/fmax2
									, (float)xstart+(xend-xstart)*(2*ks+3)/(2*nlevels)
									, (float)yend-(yend-ystart)*freq2/fmax2
									, modeltheme[YUV<<2|kc]
								);
							}
						}
#endif
						GUIPrint((float)xstart, (float)xstart, (float)ystart, 1, "%02d %12.2f"
							, kh%modelnctx
							, modelcsizes[kh]
						);
#if 0
						GUIPrint((float)xstart, (float)xstart, (float)ystart, 1, "%02d %12.2f %+12.2f"
							, kh%modelnctx
							, modelcsizes[kh]
							, modelmsizes[kh]-modelcsizes[kh]
						);
						GUIPrint((float)xstart, (float)xstart, (float)ystart+tdy, 1, "   %12.2f"
							, modelmsizes[kh]
						);
#endif
						esizes[kc]+=modelcsizes[kh];
						msizes+=modelmsizes[kh];
					}
					esize=0;
					for(kc=0;kc<modelnch;++kc)
					{
						x4=wndw*(2*(float)kc+1)/(2*(float)modelnch)-tdx*2*12;
						GUIPrint(x4, x4, (float)(wndh>>1), 2, "%12.2f", esizes[kc]);
						esize+=esizes[kc];
					}
					x4=(float)(wndw>>1)-tdx*2*12;
					GUIPrint(x4, x4, (float)(wndh>>1)+tdy*2, 2, "%12.2f +%12.2f =%12.2f"
						, esize
						, modelstatoverhead
						, esize+modelstatoverhead
					);
					GUIPrint(x4, x4, (float)(wndh>>1)+tdy*4, 2, "%12.2f"
						, msizes
					);
				}
			}
			break;
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
#if 0
			draw_profile_xy(0, 1, 0xFF00FFFF);
			draw_profile_xy(1, 2, 0xFFFFFF00);
			draw_profile_xy(2, 0, 0xFFFF00FF);
#endif
		}
	}

	if(transforms_customenabled)
	{
		float ystep, x, y;
		int k;

		ystep=tdy*guizoom;
		if(transforms_mask[CT_FWD_CUSTOM]||transforms_mask[CT_INV_CUSTOM])
		{
			const char chnames[]="rgb";
			uint8_t per[4]={0};
			double mfwd[9]={0}, minv[9]={0};

			/*custom color transform params*/
			x=buttons[0].x1;
			y=buttons[0].y1;
			/*
			P0  rgb
			r += (-0x00*g-0x00*b)>>6
			g += (-0x00*r-0x00*b)>>6
			b += (-0x00*r-0x00*g)>>6
			g += (-0x00*r-0x00*b)>>6
			0123456789012345678901234
			*/
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
			rct_custom_getmatrix(mfwd, 1);
			rct_custom_getmatrix(minv, 0);
			GUIPrint(0, x, y, guizoom, "Fwd:");
			y+=ystep;
			for(k=0;k<9-2;k+=3)
			{
				double *v=mfwd+k;
				GUIPrint(0, x, y, guizoom, "%10.6f %10.6f %10.6f", v[0], v[1], v[2]);
				y+=ystep;
			}
			GUIPrint(0, x, y, guizoom, "Inv:");
			y+=ystep;
			for(k=0;k<9-2;k+=3)
			{
				double *v=minv+k;
				GUIPrint(0, x, y, guizoom, "%10.6f %10.6f %10.6f", v[0], v[1], v[2]);
				y+=ystep;
			}
		}
		if(
			transforms_mask[ST_FWD_CUSTOM]||transforms_mask[ST_INV_CUSTOM]
			||transforms_mask[ST_CONVTEST]
			||transforms_mask[ST_CONVTEST2]
			||transforms_mask[ST_FWD_CC]||transforms_mask[ST_INV_CC]
		)
		{
			int c0, *params, sum, neg;

			c0=set_bk_color(0x80FFFFFF);
			params=custom_params+CUSTOM_NNB*2*custom_pred_ch_idx;
			sum=0;
			neg=0;
			for(k=0;k<CUSTOM_NNB*2;k+=2)
				sum+=params[k];
			neg=sum<0;
			sum=abs(sum);
			/*custom spatial transform params*/
			x=buttons[1].x1;
			y=buttons[1].y1;
			/*
			|0000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000
			|0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
			|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			|-0x00.0000-0x00.0000 -0x00.0000-0x00.0000
			*/
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
		}
	}
	else if(im1&&(transforms_mask[ST_FWD_OLS4]||transforms_mask[ST_INV_OLS4]))
	{
		int c0, kc, ky, kx;
		float ystep, x, y;

		c0=set_bk_color(0x80FFFFFF);
		ystep=tdy*guizoom;
		x=buttons[5].x1;
		y=buttons[5].y1;
		if(ols4_cache)
			GUIPrint(0, x, y-ystep, guizoom, "M %d%d%d%d%d%d%d%d"
				, ols4_cache>>7&1
				, ols4_cache>>6&1
				, ols4_cache>>5&1
				, ols4_cache>>4&1
				, ols4_cache>>3&1
				, ols4_cache>>2&1
				, ols4_cache>>1&1
				, ols4_cache>>0&1
			);
		GUIPrint(0, x, y, guizoom, "%8d%9.6f%9.6f%9.6f%9.6f"
			, ols4_period
			, ols4_lr[0]
			, ols4_lr[1]
			, ols4_lr[2]
			, ols4_lr[3]
		);
		for(kc=0;kc<im1->nch;++kc)
		{
			for(ky=0;ky<=OLS4_RMAX;++ky)
			{
				for(kx=0;kx<(OLS4_RMAX<<1|1);++kx)
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
		const int16_t *layer;
		int c0, coeffsum, k, ky, kx;
		float x, y, ystep;

		c0=set_bk_color(0x80FFFFFF);
		x=0;
		y=0;
		ystep=tdy*guizoom;
		layer=lossyconv_params+5*5*lossyconv_page;
		coeffsum=0;
		for(k=0;k<25;++k)
			coeffsum+=layer[k]>>1;
		x=buttons[6].x1;
		y=buttons[6].y1;
		/*Layer 1/4  Ch 1/3*/
		GUIPrint_append(0, 0, 0, 0, 0, "Layer %d/4  %c %d/%d  cRCT %d  sum %f"
			, (lossyconv_page>>2)+1
			, "RGBA"[lossyconv_page&3]
			, (lossyconv_page&3)+1
			, im0->nch
			, lossyconv_causalRCT[lossyconv_page>>2]
			, coeffsum/32.
		);
		if(lossyconv_clipboard)
		{
			int16_t val=lossyconv_clipboard, aval=(int16_t)abs(val);
			GUIPrint_append(0, 0, 0, 0, 0, "  M %c%c%c.%c%c%c%c%c%c"
				, val<0?'-':' '
				, '0'+(aval>>7&1)
				, '0'+(aval>>6&1)
				, '0'+(aval>>5&1)
				, '0'+(aval>>4&1)
				, '0'+(aval>>3&1)
				, '0'+(aval>>2&1)
				, '0'+(aval>>1&1)
				, val&1?'C':' '
			);
		}
		GUIPrint_append(x, x, y+0*ystep, guizoom, 1, 0);
		/*
		|    000000000011111111112222222222333333
		|    012345678901234567890123456789012345
		|0  "Layer 1/4  C 1/3  cRCT 1"
		|    0000000000111111111122222222223333333333444444444455555555556
		|    0123456789012345678901234567890123456789012345678901234567890
		|1  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		|2  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		|3  " +II.FFFFFC  +II.FFFFFC [+II.FFFFFC] +II.FFFFFC  +II.FFFFFC "
		|4  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		|5  " +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC  +II.FFFFFC "
		|6  "XY  Offset XXX YYY"
		|7  "    Stride XXX YYY"
		*/
		for(ky=0;ky<5;++ky)
		{
			g_printed=0;
			for(kx=0;kx<5;++kx)
			{
				int16_t val=layer[5*ky+kx], aval=(int16_t)abs(val);
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
			}
			GUIPrint_append(x, x, y+(ky+1)*ystep, guizoom, 1, 0);
		}

		GUIPrint(x, x, y+6*ystep, guizoom, "XY  Offset %3d %3d", lossyconv_offset[lossyconv_page>>2<<1|0], lossyconv_offset[lossyconv_page>>2<<1|1]);
		GUIPrint(x, x, y+7*ystep, guizoom, "    Stride %3d %3d", lossyconv_stride[lossyconv_page>>2<<1|0]+1, lossyconv_stride[lossyconv_page>>2<<1|1]+1);
		set_bk_color(c0);
	}
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
	usize=0;
	csize=0;
	if(im1)
	{
		float
			x, y, x2, xstart, xend, ystart, scale,
			crformat, cr_combined;
		int RGBspace, k;

		x=(float)(wndw-300), y=tdy*2, x2=x;
		for(k=0;k<T_COUNT;++k)/*print available transforms on right*/
		{
			transforms_printname(x2, y, k, -1, transforms_mask[k]?(uint64_t)0xA0FF0000<<32|0xFFFFFFFF:0);
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
			else if(ec_adaptive)/*H.E.M.L..A.0x0000..XXXX_XXX*/
				GUIPrint(x, x, y-tdy, 1, "H %d %d %d  A 0x%04X  %s"
					, ec_expbits
					, ec_msb
					, ec_lsb
					, ec_adaptive_threshold
					, label
				);
			else
				GUIPrint(x, x, y-tdy, 1, "H %d %d %d  Static    %s"
					, ec_expbits
					, ec_msb
					, ec_lsb
					, label
				);
		}
		if(transforms)
		{
			for(k=0;k<(int)transforms->count;++k, y+=tdy)/*print applied transforms on left*/
				transforms_printname(x, y, HDR_DATA(transforms)[k], k, 0);
		}
		cr_combined=(float)(((im1->src_depth[0]>0)+(im1->src_depth[1]>0)+(im1->src_depth[2]>0)+(im1->src_depth[3]>0))/(ch_entropy[0]+ch_entropy[1]+ch_entropy[2]+ch_entropy[3]));
		xstart=20, xend=(float)wndw-330, ystart=(float)(wndh-tdy*5);
		scale=xend-xstart;
		
		usize=image_getBMPsize(im0);
		crformat=(float)(usize/filesize);
		RGBspace=1;
		for(k=0;k<CST_COMPARE;++k)
		{
			if(transforms_mask[k])
			{
				RGBspace=0;
				break;
			}
		}
		if(xstart<xend)
		{
			float maxinvcr, barw;
			int ks, ks2;

			maxinvcr=1/cr_combined;
			for(k=0;k<4;++k)/*get max CR*/
			{
				double invCR=ch_entropy[k];
				if(!acme_isinf(invCR)&&maxinvcr<invCR)
					maxinvcr=(float)invCR;
			}
			if(acme_isinf(maxinvcr))
				maxinvcr=0;
			if(maxinvcr>100)
				maxinvcr=100;
			if(maxinvcr<1/crformat)
				maxinvcr=1/crformat;
			if(xend-scale*maxinvcr<xstart)
				scale=(xend-xstart)/maxinvcr;
			if(scale<50)
				scale=50;
			xstart=xend-scale*(float)ceil(((double)maxinvcr+1)*2)*0.5f;
			draw_rect(xstart, xend, ystart, (float)wndh, 0x80808080);/*background*/
			if(scale>20)
			{
				ks=1;
				x=(float)(xend-0.1f*scale*ks);
				for(;x>xstart;++ks)
				{
					draw_line(x, ystart+1, x, (float)wndh, 0x70900090);/*draw minor scale*/
					x=(float)(xend-0.1f*scale*ks);
				}
			}
			ks=1;
			x=(float)(xend-scale*ks);
			for(ks2=1;x>xstart;++ks2)
			{
				draw_rect(x-1, x+2, ystart, (float)wndh, 0x70800080);/*draw major scale*/
				x=(float)(xend-scale*ks2);
			}
			barw=4;
			draw_rect(xend-scale/cr_combined,		xend, ystart+tdy*0.5f-barw, ystart+tdy*0.5f+barw+1, 0xFF000000);
			draw_rect(xend-scale*(float)ch_entropy[0],	xend, ystart+tdy*1.5f-barw, ystart+tdy*1.5f+barw+1, RGBspace?0xFF0000FF:0xFF404040);/*r or Y	*/
			draw_rect(xend-scale*(float)ch_entropy[1],	xend, ystart+tdy*2.5f-barw, ystart+tdy*2.5f+barw+1, RGBspace?0xFF00FF00:0xFFC00000);/*g or Cb	*/
			draw_rect(xend-scale*(float)ch_entropy[2],	xend, ystart+tdy*3.5f-barw, ystart+tdy*3.5f+barw+1, RGBspace?0xFFFF0000:0xFF0000C0);/*b or Cr	*/
			if(im1->depth[3])
				draw_rect(xend-scale*(float)ch_entropy[3], xend, ystart+tdy*4.5f-barw, ystart+tdy*4.5f+barw+1, RGBspace?0xFF808080:0xFF00C000);/*a or Cg	*/
			x=xend-scale/crformat;
			draw_line(x, ystart, x-10, ystart-10, 0xFF000000);
			draw_line(x, ystart, x+10, ystart-10, 0xFF000000);
		}
		{
			int prevtxtcolor, prevbkcolor;

			xend+=10;
			prevbkcolor=set_bk_color(0xC0C0C0C0);
			prevtxtcolor=set_text_color(0xFF000000);
			GUIPrint(xend, xend, ystart-tdy*2, 1, "Bitmap Size             %9.0f", usize);
			GUIPrint(xend, xend, ystart-tdy  , 1, "Format        %8.4f%% %9lld", 100./crformat, filesize);
			set_bk_color(0xE0FFFFFF);
			csize=usize/cr_combined;
			prevtxtcolor=set_text_color(0xFF000000);GUIPrint(xend, xend, ystart      , 1, "Combined      %8.4f%% %12.2f", 100./cr_combined, csize);
			set_bk_color(0xC0C0C0C0);
			set_text_color(RGBspace?0xFF0000FF:0xFF404040);
			GUIPrint(xend, xend, ystart+tdy  , 1, "%c     %7d %8.4f%% %12.2f"
				, RGBspace?'R':'Y'
				, im1->depth[0]
				, 100.*ch_entropy[0]
				, (double)im1->iw*im1->ih*im1->src_depth[0]*ch_entropy[0]/8
			);
			set_text_color(RGBspace?0xFF00C000:0xFFC00000);
			GUIPrint(xend, xend, ystart+tdy*2, 1, "%c     %7d %8.4f%% %12.2f"
				, RGBspace?'G':'U'
				, im1->depth[1]
				, 100.*ch_entropy[1]
				, (double)im1->iw*im1->ih*im1->src_depth[1]*ch_entropy[1]/8
			);
			set_text_color(RGBspace?0xFFFF0000:0xFF0000C0);
			GUIPrint(xend, xend, ystart+tdy*3, 1, "%c     %7d %8.4f%% %12.2f"
				, RGBspace?'B':'V'
				, im1->depth[2]
				, 100.*ch_entropy[2]
				, (double)im1->iw*im1->ih*im1->src_depth[2]*ch_entropy[2]/8
			);
			if(im1->depth[3])/*src_depth[3] is zero assuming no alpha*/
			{
				set_text_color(RGBspace?0xFF404040:0xFF00C000);
				GUIPrint(xend, xend, ystart+tdy*4, 1, "%c     %7d %8.4f%% %12.2f"
					, RGBspace?'A':'W'
					, im1->depth[3]
					, 100.*ch_entropy[3]
					, (double)im1->iw*im1->ih*im1->src_depth[3]*ch_entropy[3]/8
				);
			}
			set_text_color(prevtxtcolor);
			set_bk_color(prevbkcolor);
		}
		if(drawCRhistory)
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
			if(combCRhist[idx][3]<combCRhist[idx2][3])/*ratio improved*/
			{
				draw_triangle(cx-9, cy, cx+9, cy, cx, cy+9, 0xC080FF80);
				draw_triangle(cx-5, cy, cx+5, cy, cx, cy+5, 0xC0006000);
			}
			else if(combCRhist[idx][3]>combCRhist[idx2][3])/*ratio worsened*/
			{
				draw_triangle(cx-9, cy, cx+9, cy, cx, cy-9, 0xC08080FF);
				draw_triangle(cx-5, cy, cx+5, cy, cx, cy-5, 0xC0000060);
			}
			else/*exact same ratio*/
			{
				draw_line(cx-5, cy-5, cx+5, cy+5, 0xC0808080);
				draw_line(cx-5, cy+5, cx+5, cy-5, 0xC0808080);
			}
			for(k=0;k<combCRhist_SIZE-1;++k)
			{
				if(k!=combCRhist_idx-1)
				{
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][0]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][0]*g2), RGBspace?0xFF0000FF:0xFF404040);/*r or Y	*/
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][1]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][1]*g2), RGBspace?0xFF00FF00:0xFFC00000);/*g or Cb	*/
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][2]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][2]*g2), RGBspace?0xFFFF0000:0xFF0000C0);/*b or Cr	*/
					draw_line(xstart+(float)(k<<combCRhist_logDX), wndh-(float)(combCRhist[k][3]*g2), xstart+(float)((k+1)<<combCRhist_logDX), wndh-(float)(combCRhist[k+1][3]*g2), 0xC0000000);
				}
			}
			draw_line(xstart, wndh-g2, xstart+(combCRhist_SIZE<<combCRhist_logDX), wndh-g2, 0xC0000000);
		}
		/*grad2 info*/
#if 1
		if(transforms_mask[ST_FWD_CUSTOM3]||transforms_mask[ST_INV_CUSTOM3])
		{
			const int width=C3_REACH<<2|2;
			const int16_t *coeffs[]=
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
			int kc;

			x=(float)(wndw>>1);
			y=tdy*zoom*4;
			for(kc=0;kc<9;++kc)
			{
				for(k=0;k<C3_REACH;++k)
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
			for(k=0;k<PW2_NPRED;++k)
			{
				GUIPrint(0, x, y, 0.9f, "p%d %10f", k, pw2_errors[k]);
				y+=tdy*0.9f;
			}
		}
#endif
		if(im1->rct&0x8000)
		{
			int uc0=0, vc0=0, vc1=0;
			const uint8_t *perm=0;

			perm=crct2_unpack(im1->rct, &uc0, &vc0, &vc1);
			GUIPrint(0, 0, 0, 1, "WH %dx%d  D0[%d %d %d %d] D[%d %d %d %d]  RCT%c%c%c_%02d_%02d%02d/%02d  Z %13.6f"
				, im0->iw
				, im0->ih
				, im0->src_depth[0]
				, im0->src_depth[1]
				, im0->src_depth[2]
				, im0->src_depth[3]
				, im1->depth[0]
				, im1->depth[1]
				, im1->depth[2]
				, im1->depth[3]
				, '0'+perm[0]
				, '0'+perm[1]
				, '0'+perm[2]
				, uc0
				, vc0
				, vc1
				, 1<<RCTBITS
				, imzoom
			);
		}
		else
		{
			GUIPrint(0, 0, 0, 1, "WH %dx%d  D0[%d %d %d %d] D[%d %d %d %d]  RCT%2d %s  Z %13.6f"
				, im0->iw
				, im0->ih
				, im0->src_depth[0]
				, im0->src_depth[1]
				, im0->src_depth[2]
				, im0->src_depth[3]
				, im1->depth[0]
				, im1->depth[1]
				, im1->depth[2]
				, im1->depth[3]
				, im1->rct
				, rct_names[im1->rct]
				, imzoom
			);
		}
	}
	{
		const char *mode_str=0, *view_str=0;
		static double t=0;
		double t2=time_ms();

		switch(mode)
		{
		case VIS_ANALYSIS:		mode_str="Analysis";		break;
		case VIS_HISTOGRAM:		mode_str="Histogram";		break;
		case VIS_MODEL:			mode_str="Model";		break;
		case VIS_JOINT_HISTOGRAM:	mode_str="Joint Histogram";	break;
#ifdef ENABLE_L1WEIGHTS
		case VIS_L1WEIGHTS:		mode_str="L1 Weights";		break;
#endif
		case VIS_IMAGE:			mode_str="Image View";		break;
		case VIS_SHANNON:		mode_str="Zipf View";		break;
		case VIS_IMAGE_TRICOLOR:	mode_str="Tricolor";		break;
		}
		GUIPrint(0, 0, tdy, 1, "timer 0x%04X, fps%11f, [%2d/%2d] %s  dist=%d  BPD %6.4f"
			, timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str, g_dist, usize>1e-9?8*csize/usize:0
		);
		if(mode==VIS_IMAGE||mode==VIS_SHANNON)
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
				char bin[17]={0};
				int j;

				j=0;
				bin[j++]='0'+(abacvis_low>>15&1);
				bin[j++]='0'+(abacvis_low>>14&1);
				bin[j++]='0'+(abacvis_low>>13&1);
				bin[j++]='0'+(abacvis_low>>12&1);
				bin[j++]='0'+(abacvis_low>>11&1);
				bin[j++]='0'+(abacvis_low>>10&1);
				bin[j++]='0'+(abacvis_low>> 9&1);
				bin[j++]='0'+(abacvis_low>> 8&1);
				bin[j++]='0'+(abacvis_low>> 7&1);
				bin[j++]='0'+(abacvis_low>> 6&1);
				bin[j++]='0'+(abacvis_low>> 5&1);
				bin[j++]='0'+(abacvis_low>> 4&1);
				bin[j++]='0'+(abacvis_low>> 3&1);
				bin[j++]='0'+(abacvis_low>> 2&1);
				bin[j++]='0'+(abacvis_low>> 1&1);
				bin[j++]='0'+(abacvis_low>> 0&1);
				bin[j++]='\0';
				GUIPrint(0, 0, tdy*3, 2, "0b%s/%d", bin, abacvis_range);
			}
		}
		else if(mode==VIS_HISTOGRAM)
		{
			switch(histmode)
			{
			case HISTMODE_LINEAR:		mode_str="Linear";break;
			case HISTMODE_LINEAR_NOSPIKE:	mode_str="Linear (no spike)";break;
			case HISTMODE_LOG2:		mode_str="Shannon";break;
			case HISTMODE_LOG2_NOSPIKE:	mode_str="Shannon (no spike)";break;
			}
			GUIPrint(0, 0, tdy*2, 1, "%d/%d %s", histmode+1, HISTMODE_COUNT, mode_str);
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
int io_quit_request(void)/*return 1 to exit*/
{
	return 1;
}
void io_cleanup(void)/*cleanup*/
{
	free(im0);
	free(im1);
	array_free(&cpu_vertices);
}
