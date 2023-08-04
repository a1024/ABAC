#include"pxview3d.h"
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
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
int iw=0, ih=0, nch0=0;
unsigned char *im0, *image=0;//stride=4
int *im2=0, *im3=0;
unsigned image_txid[3]={0};//0: interleaved channels, 1: separate channels, 2: histogram as texture

typedef enum VisModeEnum
{
	VIS_PLANES,
	VIS_MESH,
	VIS_MESH_SEPARATE,
	VIS_IMAGE_TRICOLOR,
	VIS_IMAGE,
	VIS_IMAGE_BLOCK,
//	VIS_IMAGE_E24,//experiment 24
	VIS_DWT_BLOCK,
	VIS_HISTOGRAM,
	VIS_JOINT_HISTOGRAM,

	VIS_COUNT,
} VisMode;
int mode=VIS_IMAGE;

#if 0
typedef enum ColorTransformTypeEnum
{
	CT_NONE,
	CT_YCoCg,
	CT_YCoCgT,
	CT_XGZ,
	CT_XYZ,
	CT_EXP,
	CT_LEARNED,
	CT_CUSTOM,

	CT_COUNT,
} ColorTransformType;
int color_transform=CT_NONE;

typedef enum SpatialTransformTypeEnum
{
	ST_NONE,
	ST_DIFF2D,
	ST_UNPLANE,
	ST_LAZY,
	ST_HAAR,
	ST_SQUEEZE,
	ST_CDF53,
	ST_CDF97,
	ST_CUSTOM,

	ST_COUNT,
} SpatialTransformType;
int spatialtransform=ST_NONE;

typedef struct TransformTypeStruct
{
	short is_spatial, id;
} TransformType;
#endif
typedef enum TransformTypeEnum
{
	T_NONE,

	CT_FWD_YCoCg,		CT_INV_YCoCg,//HEVC, VVC
	CT_FWD_YCoCb,		CT_INV_YCoCb,
	CT_FWD_XGZ,			CT_INV_XGZ,
	CT_FWD_XYZ,			CT_INV_XYZ,
//	CT_FWD_EXP,			CT_INV_EXP,
	CT_FWD_ADAPTIVE,	CT_INV_ADAPTIVE,
	CT_FWD_CUSTOM,		CT_INV_CUSTOM,

	CST_SEPARATOR,
	
	ST_FWD_CUSTOM2,		ST_INV_CUSTOM2,
//	ST_FWD_LOGIC,		ST_INV_LOGIC,
	ST_FWD_LEARNED,		ST_INV_LEARNED,
#ifdef ALLOW_OPENCL
	ST_FWD_LEARNED_GPU,	ST_INV_LEARNED_GPU,
#endif
//	ST_FWD_CFL,			ST_INV_CFL,
//	ST_FWD_JOINT,		ST_INV_JOINT,
//	ST_FWD_HYBRID3,		ST_INV_HYBRID3,
	
//	ST_FWD_DIFF2D,		ST_INV_DIFF2D,
//	ST_FWD_GRAD2,		ST_INV_GRAD2,
//	ST_FWD_ADAPTIVE,	ST_INV_ADAPTIVE,
	ST_FWD_JXL,			ST_INV_JXL,
	ST_FWD_MM,			ST_INV_MM,
	ST_FWD_JMJ,			ST_INV_JMJ,
//	ST_FWD_MEDIAN,		ST_INV_MEDIAN,
//	ST_FWD_DCT3PRED,	ST_INV_DCT3PRED,
//	ST_FWD_PATHPRED,	ST_INV_PATHPRED,
//	ST_FWD_HPF,			ST_INV_HPF,
	ST_FWD_CUSTOM,		ST_INV_CUSTOM,
	ST_FWD_GRADPRED,	ST_INV_GRADPRED,
	ST_FWD_SORTNB,		ST_INV_SORTNB,
//	ST_FWD_BITWISE,		ST_INV_BITWISE,
	ST_FWD_DCT4,		ST_INV_DCT4,
//	ST_FWD_DCT8,		ST_INV_DCT8,
	ST_FWD_SHUFFLE,		ST_INV_SHUFFLE,
//	ST_FWD_SPLIT,		ST_INV_SPLIT,
	
	ST_FWD_LAZY,		ST_INV_LAZY,
	ST_FWD_HAAR,		ST_INV_HAAR,
	ST_FWD_SQUEEZE,		ST_INV_SQUEEZE,
	ST_FWD_LEGALL53,	ST_INV_LEGALL53,
//	ST_FWD_CDF97,		ST_INV_CDF97,
//	ST_FWD_GRAD_DWT,	ST_INV_GRAD_DWT,
//	ST_FWD_DEC_DWT,		ST_INV_DEC_DWT,
	ST_FWD_EXPDWT,		ST_INV_EXPDWT,
	ST_FWD_CUSTOM_DWT,	ST_INV_CUSTOM_DWT,

	T_COUNT,
} TransformType;
int transforms_customenabled=0;
char transforms_mask[T_COUNT]={0};
ArrayHandle transforms=0;//array of chars
float guizoom=1.25f;

double av_rmse=0, g_lr=1e-10;

int profile_idx=0;
ArrayHandle losshist=0;
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
float ch_cr[4]={0};
int usage[4]={0};

#define combCRhist_SIZE 128
#define combCRhist_logDX 2
float combCRhist[combCRhist_SIZE][4]={0}, combCRhist_max=1;
int combCRhist_idx=0;

int show_full_image=0;
int space_not_color=0;

void transforms_update()
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
			transforms_customenabled|=
				tid2==CT_FWD_CUSTOM||
				tid2==CT_INV_CUSTOM||
				tid2==ST_FWD_CUSTOM||
				tid2==ST_INV_CUSTOM||
				tid2==ST_FWD_EXPDWT||
				tid2==ST_INV_EXPDWT||
				tid2==ST_FWD_CUSTOM_DWT||
				tid2==ST_INV_CUSTOM_DWT||
			//	tid2==ST_FWD_HYBRID3||
			//	tid2==ST_INV_HYBRID3||
			//	tid2==ST_FWD_ADAPTIVE||
			//	tid2==ST_INV_ADAPTIVE||
			//	tid2==ST_FWD_JXL||
			//	tid2==ST_INV_JXL||
				tid2==ST_FWD_SORTNB||
				tid2==ST_INV_SORTNB;
		}
	}
}
void transforms_removebyid(unsigned tid)
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
void transforms_removeall()
{
	array_free(&transforms);
	transforms_customenabled=0;
	memset(transforms_mask, 0, T_COUNT);
}
void transforms_append(unsigned tid)
{
	if(tid<T_COUNT&&tid!=CST_SEPARATOR)
	{
		if(!transforms)
		{
			ARRAY_ALLOC(char, transforms, &tid, 1, 0, 0);
			transforms_mask[tid]|=1;
			transforms_customenabled|=
				tid==CT_FWD_CUSTOM||
				tid==CT_INV_CUSTOM||
				tid==ST_FWD_CUSTOM||
				tid==ST_INV_CUSTOM||
				tid==ST_FWD_EXPDWT||
				tid==ST_INV_EXPDWT||
				tid==ST_FWD_CUSTOM_DWT||
				tid==ST_INV_CUSTOM_DWT||
			//	tid==ST_FWD_HYBRID3||
			//	tid==ST_INV_HYBRID3||
			//	tid==ST_FWD_ADAPTIVE||
			//	tid==ST_INV_ADAPTIVE||
			//	tid==ST_FWD_JXL||
			//	tid==ST_INV_JXL||
				tid==ST_FWD_SORTNB||
				tid==ST_INV_SORTNB;
		}
		else
		{
			int idx=-1;//first idx of a transform of this type
			if(GET_KEY_STATE(KEY_CTRL))//replace all transforms of this type
			{
				for(int k=0;k<transforms->count;)
				{
					if((tid<CST_SEPARATOR)==(transforms->data[k]<CST_SEPARATOR))
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
void transforms_printname(float x, float y, unsigned tid, int place, long long highlight)
{
	const char *a=0;
	switch(tid)
	{
	case T_NONE:				a="NONE";					break;
	case CT_FWD_YCoCg:			a="C  Fwd YCoCg-R";			break;
	case CT_INV_YCoCg:			a="C  Inv YCoCg-R";			break;
	case CT_FWD_YCoCb:			a="C  Fwd YCoCb-R";			break;
	case CT_INV_YCoCb:			a="C  Inv YCoCb-R";			break;
	case CT_FWD_XGZ:			a="C  Fwd XGZ";				break;
	case CT_INV_XGZ:			a="C  Inv XGZ";				break;
	case CT_FWD_XYZ:			a="C  Fwd XYZ";				break;
	case CT_INV_XYZ:			a="C  Inv XYZ";				break;
//	case CT_FWD_EXP:			a="C  Fwd Experimental";	break;
//	case CT_INV_EXP:			a="C  Inv Experimental";	break;
	case CT_FWD_ADAPTIVE:		a="C  Fwd Adaptive";		break;
	case CT_INV_ADAPTIVE:		a="C  Inv Adaptive";		break;
	case CT_FWD_CUSTOM:			a="C  Fwd CUSTOM";			break;
	case CT_INV_CUSTOM:			a="C  Inv CUSTOM";			break;

	case CST_SEPARATOR:			a="";						break;

	case ST_FWD_CUSTOM2:		a=" S Fwd CUSTOM2";			break;
	case ST_INV_CUSTOM2:		a=" S Inv CUSTOM2";			break;
//	case ST_FWD_LOGIC:			a=" S Fwd Logic";			break;
//	case ST_INV_LOGIC:			a=" S Inv Logic";			break;
	case ST_FWD_LEARNED:		a=" S Fwd Learned";			break;
	case ST_INV_LEARNED:		a=" S Inv Learned";			break;
#ifdef ALLOW_OPENCL
	case ST_FWD_LEARNED_GPU:	a=" S Fwd Learned GPU";		break;
	case ST_INV_LEARNED_GPU:	a=" S Inv Learned GPU";		break;
#endif
//	case ST_FWD_CFL:			a=" S Fwd CfL";				break;
//	case ST_INV_CFL:			a=" S Inv CfL";				break;
//	case ST_FWD_JOINT:			a=" S Fwd Joint";			break;
//	case ST_INV_JOINT:			a=" S Inv Joint";			break;
//	case ST_FWD_HYBRID3:		a=" S Fwd Hybrid3";			break;
//	case ST_INV_HYBRID3:		a=" S Inv Hybrid3";			break;

//	case ST_FWD_DIFF2D:			a=" S Fwd 2D derivative";	break;
//	case ST_INV_DIFF2D:			a=" S Inv 2D derivative";	break;
//	case ST_FWD_HPF:			a=" S Fwd HPF";				break;
//	case ST_INV_HPF:			a=" S Inv HPF";				break;
//	case ST_FWD_GRAD2:			a=" S Fwd Grad2";			break;
//	case ST_INV_GRAD2:			a=" S Inv Grad2";			break;
//	case ST_FWD_ADAPTIVE:		a=" S Fwd Adaptive";		break;
//	case ST_INV_ADAPTIVE:		a=" S Inv Adaptive";		break;
	case ST_FWD_JXL:			a=" S Fwd JXL";				break;
	case ST_INV_JXL:			a=" S Inv JXL";				break;
	case ST_FWD_MM:				a=" S Fwd MM";				break;
	case ST_INV_MM:				a=" S Inv MM";				break;
	case ST_FWD_JMJ:			a=" S Fwd JMJ";				break;
	case ST_INV_JMJ:			a=" S Inv JMJ";				break;
	case ST_FWD_SORTNB:			a=" S Fwd Sort Nb";			break;
	case ST_INV_SORTNB:			a=" S Inv Sort Nb";			break;
//	case ST_FWD_MEDIAN:			a=" S Fwd Median";			break;
//	case ST_INV_MEDIAN:			a=" S Inv Median";			break;
//	case ST_FWD_DCT3PRED:		a=" S Fwd DCT3 Predictor";	break;
//	case ST_INV_DCT3PRED:		a=" S Inv DCT3 Predictor";	break;
//	case ST_FWD_PATHPRED:		a=" S Fwd Path Predictor";	break;
//	case ST_INV_PATHPRED:		a=" S Inv Path Predictor";	break;
	case ST_FWD_GRADPRED:		a=" S Fwd Gradient";		break;
	case ST_INV_GRADPRED:		a=" S Inv Gradient";		break;
//	case ST_FWD_BITWISE:		a=" S Fwd Bitwise";			break;
//	case ST_INV_BITWISE:		a=" S Inv Bitwise";			break;
	case ST_FWD_DCT4:			a=" S Fwd DCT4";			break;
	case ST_INV_DCT4:			a=" S Inv DCT4";			break;
//	case ST_FWD_DCT8:			a=" S Fwd DCT8";			break;
//	case ST_INV_DCT8:			a=" S Inv DCT8";			break;
	case ST_FWD_SHUFFLE:		a=" S Fwd Shuffle";			break;
	case ST_INV_SHUFFLE:		a=" S Inv Shuffle";			break;
//	case ST_FWD_SPLIT:			a=" S Fwd Split";			break;
//	case ST_INV_SPLIT:			a=" S Inv Split";			break;
	case ST_FWD_CUSTOM:			a=" S Fwd CUSTOM";			break;
	case ST_INV_CUSTOM:			a=" S Inv CUSTOM";			break;
	case ST_FWD_LAZY:			a=" S Fwd Lazy DWT";		break;
	case ST_INV_LAZY:			a=" S Inv Lazy DWT";		break;
	case ST_FWD_HAAR:			a=" S Fwd Haar";			break;
	case ST_INV_HAAR:			a=" S Inv Haar";			break;
	case ST_FWD_SQUEEZE:		a=" S Fwd Squeeze";			break;
	case ST_INV_SQUEEZE:		a=" S Inv Squeeze";			break;
	case ST_FWD_LEGALL53:		a=" S Fwd LeGall 5/3";		break;
	case ST_INV_LEGALL53:		a=" S Inv LeGall 5/3";		break;
//	case ST_FWD_CDF97:			a=" S Fwd CDF 9/7";			break;
//	case ST_INV_CDF97:			a=" S Inv CDF 9/7";			break;
//	case ST_FWD_GRAD_DWT:		a=" S Fwd Gradient DWT";	break;
//	case ST_INV_GRAD_DWT:		a=" S Inv Gradient DWT";	break;
//	case ST_FWD_DEC_DWT:		a=" S Fwd Dec. DWT";		break;
//	case ST_INV_DEC_DWT:		a=" S Inv Dec. DWT";		break;
	case ST_FWD_EXPDWT:			a=" S Fwd Exp DWT";			break;
	case ST_INV_EXPDWT:			a=" S Inv Exp DWT";			break;
	case ST_FWD_CUSTOM_DWT:		a=" S Fwd CUSTOM DWT";		break;
	case ST_INV_CUSTOM_DWT:		a=" S Inv CUSTOM DWT";		break;
	default:					a="ERROR";					break;
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

int send_image_separate_subpixels(unsigned char *image, int iw, int ih, unsigned *txid)
{
	int res=iw*ih;
	if(im2)
		free(im2);
	im2=(int*)malloc((size_t)res*3*sizeof(int));
	if(!im2)
		return 0;
	for(int k=0;k<res;++k)
	{
		unsigned char r=image[k<<2], g=image[k<<2|1], b=image[k<<2|2];
		int idx=3*k;
		im2[idx  ]=0xC0000000|r;
		im2[idx+1]=0xC0000000|g<<8;
		im2[idx+2]=0xC0000000|b<<16;
	}

	//char buf[MAX_PATH];
	//GetCurrentDirectoryA(MAX_PATH, buf);
	//lodepng_encode_file("out.PNG", im2, iw*3, ih, LCT_RGBA, 8);//

	
	if(im3)
		free(im3);
	im3=(int*)malloc((size_t)res*3*sizeof(int));
	if(!im3)
		return 0;
	for(int k=0;k<res;++k)
	{
		im3[k]=0xC0000000|image[k<<2  ];
		im3[res+k]=0xC0000000|image[k<<2|1]<<8;
		im3[(res<<1)+k]=0xC0000000|image[k<<2|2]<<16;
	}

	if(!*txid)
		glGenTextures(2, txid);
	send_texture_pot(txid[0], im2, iw*3, ih, 1);
	send_texture_pot(txid[1], im3, iw, ih*3, 1);

	//free(im2);
	return 1;
}

void chart_planes_update(unsigned char *im, int iw, int ih, ArrayHandle *cpuv, unsigned *gpuv)
{
	int nv=iw*ih*3*6, nf=nv*5;//subpixel count * 6 vertices
	if(!*cpuv||cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;

	for(int ky=0, kv=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			//int idx=3*(iw*ky+kx);
			for(int kc=0;kc<3;++kc, kv+=30)
			{
				unsigned char val=im[idx|kc];
				//unsigned char val=im[idx+kc]>>(kc<<3)&0xFF;
				
				//planes
#if 1
				float
					x1=-(iw-1-(kx+(float)kc/3)),
					x2=-(iw-1-(kx+(float)(kc+0.9f)/3)),
					y1=-(float)ky,
					y2=-((float)ky+0.9f),
					z=(float)val*pixel_amplitude/255;
				float
					tx=(3*kx+kc+0.5f)/(3*iw),
					ty=(ky+0.5f)/ih;
				
				vertices[kv   ]=x1, vertices[kv+ 1]=y1, vertices[kv+ 2]=z, vertices[kv+ 3]=tx, vertices[kv+ 4]=ty;
				vertices[kv+ 5]=x1, vertices[kv+ 6]=y2, vertices[kv+ 7]=z, vertices[kv+ 8]=tx, vertices[kv+ 9]=ty;
				vertices[kv+10]=x2, vertices[kv+11]=y2, vertices[kv+12]=z, vertices[kv+13]=tx, vertices[kv+14]=ty;
				
				vertices[kv+15]=x2, vertices[kv+16]=y2, vertices[kv+17]=z, vertices[kv+18]=tx, vertices[kv+19]=ty;
				vertices[kv+20]=x2, vertices[kv+21]=y1, vertices[kv+22]=z, vertices[kv+23]=tx, vertices[kv+24]=ty;
				vertices[kv+25]=x1, vertices[kv+26]=y1, vertices[kv+27]=z, vertices[kv+28]=tx, vertices[kv+29]=ty;
#endif

				//rectangles
#if 0
				float
					x1=kx+(float)kc/3,
					x2=kx+(float)(kc+0.9f)/3,
					y1=(float)ky,
					y2=(float)ky+0.9f,
					z1=0,
					z2=(float)val*pixel_amplitude/255;
				float
					tx=(3*kx+kc+0.5f)/(3*iw-1),
					ty=(ky+0.5f)/(ih-1);
				
				vertices[kv   ]=x1, vertices[kv+ 1]=y1, vertices[kv+ 2]=z1, vertices[kv+ 3]=tx, vertices[kv+ 4]=ty;
				vertices[kv+ 5]=x2, vertices[kv+ 6]=y2, vertices[kv+ 7]=z1, vertices[kv+ 8]=tx, vertices[kv+ 9]=ty;
				vertices[kv+10]=x2, vertices[kv+11]=y2, vertices[kv+12]=z2, vertices[kv+13]=tx, vertices[kv+14]=ty;
				
				vertices[kv+15]=x2, vertices[kv+16]=y2, vertices[kv+17]=z2, vertices[kv+18]=tx, vertices[kv+19]=ty;
				vertices[kv+20]=x1, vertices[kv+21]=y1, vertices[kv+22]=z2, vertices[kv+23]=tx, vertices[kv+24]=ty;
				vertices[kv+25]=x1, vertices[kv+26]=y1, vertices[kv+27]=z1, vertices[kv+28]=tx, vertices[kv+29]=ty;
#endif
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
void chart_mesh_update(unsigned char *im, int iw, int ih, ArrayHandle *cpuv, unsigned *gpuv)
{
	int nv=(iw-1)*(ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 6 vertices * 5 floats
	if(!*cpuv||cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;

	for(int ky=0, kv=0;ky<ih-1;++ky)
	{
		for(int kx=0;kx<iw-1;++kx)
		{
			unsigned char comp[]=
			{
				im[(iw* ky   +kx  )<<2  ],
				im[(iw* ky   +kx  )<<2|1],
				im[(iw* ky   +kx  )<<2|2],
				im[(iw* ky   +kx+1)<<2  ],
				im[(iw*(ky+1)+kx  )<<2  ],
				im[(iw*(ky+1)+kx  )<<2|1],
				im[(iw*(ky+1)+kx  )<<2|2],
				im[(iw*(ky+1)+kx+1)<<2  ],
			};
			for(int kc=0;kc<3;++kc, kv+=30)
			{
				//mesh
				float
					x1=-(iw-1-(kx+(float)kc/3)),
					x2=-(iw-1-(kx+(float)(kc+0.9f)/3)),
					y1=-((float)ky),
					y2=-((float)ky+0.9f),
					z00=(float)comp[kc  ]*pixel_amplitude/255,
					z01=(float)comp[kc+1]*pixel_amplitude/255,
					z10=(float)comp[kc+4]*pixel_amplitude/255,
					z11=(float)comp[kc+5]*pixel_amplitude/255;
				float
					tx1=(3*kx+kc+0.5f)/(3*iw),
					ty1=(ky+0.5f)/ih,
					tx2=(3*kx+kc+0.9f+0.5f)/(3*iw),
					ty2=(ky+0.9f+0.5f)/ih;
				//float
				//	tx=(3*kx+kc+0.5f)/(3*iw),
				//	ty=(ky+0.5f)/ih;
				
				vertices[kv   ]=x1, vertices[kv+ 1]=y1, vertices[kv+ 2]=z00, vertices[kv+ 3]=tx1, vertices[kv+ 4]=ty1;
				vertices[kv+ 5]=x1, vertices[kv+ 6]=y2, vertices[kv+ 7]=z10, vertices[kv+ 8]=tx1, vertices[kv+ 9]=ty2;
				vertices[kv+10]=x2, vertices[kv+11]=y2, vertices[kv+12]=z11, vertices[kv+13]=tx2, vertices[kv+14]=ty2;
				
				vertices[kv+15]=x2, vertices[kv+16]=y2, vertices[kv+17]=z11, vertices[kv+18]=tx2, vertices[kv+19]=ty2;
				vertices[kv+20]=x2, vertices[kv+21]=y1, vertices[kv+22]=z01, vertices[kv+23]=tx2, vertices[kv+24]=ty1;
				vertices[kv+25]=x1, vertices[kv+26]=y1, vertices[kv+27]=z00, vertices[kv+28]=tx1, vertices[kv+29]=ty1;
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
void chart_mesh_sep_update(unsigned char *im, int iw, int ih, ArrayHandle *cpuv, unsigned *gpuv)
{
	int nv=(iw-1)*(ih-1)*3*6, nf=nv*5;//pixel count * 3 colors * 2 triangles * 3 vertices * 5 floats
	if(!*cpuv||cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;

	for(int ky=0, kv=0;ky<ih-1;++ky)
	{
		for(int kx=0;kx<iw-1;++kx)
		{
			for(int kc=0;kc<3;++kc, kv+=30)
			{
				//mesh
				float
					x1=-(iw-1-(float)kx),
					x2=-(iw-1-((float)kx+0.9f)),
					y1=-(float)ky,
					y2=-((float)ky+0.9f),
					z00=(float)im[(iw* ky   +kx  )<<2|kc]*pixel_amplitude/255+kc*mesh_separation,
					z01=(float)im[(iw* ky   +kx+1)<<2|kc]*pixel_amplitude/255+kc*mesh_separation,
					z10=(float)im[(iw*(ky+1)+kx  )<<2|kc]*pixel_amplitude/255+kc*mesh_separation,
					z11=(float)im[(iw*(ky+1)+kx+1)<<2|kc]*pixel_amplitude/255+kc*mesh_separation;
				float
					tx1=(kx+0.5f)/iw,
					tx2=(kx+0.5f+0.9f)/iw,
					ty1=(ih*kc+ky+0.5f)/(3*ih),
					ty2=(ih*kc+ky+0.9f+0.5f)/(3*ih);
				
				vertices[kv   ]=x1, vertices[kv+ 1]=y1, vertices[kv+ 2]=z00, vertices[kv+ 3]=tx1, vertices[kv+ 4]=ty1;
				vertices[kv+ 5]=x1, vertices[kv+ 6]=y2, vertices[kv+ 7]=z10, vertices[kv+ 8]=tx1, vertices[kv+ 9]=ty2;
				vertices[kv+10]=x2, vertices[kv+11]=y2, vertices[kv+12]=z11, vertices[kv+13]=tx2, vertices[kv+14]=ty2;
				
				vertices[kv+15]=x2, vertices[kv+16]=y2, vertices[kv+17]=z11, vertices[kv+18]=tx2, vertices[kv+19]=ty2;
				vertices[kv+20]=x2, vertices[kv+21]=y1, vertices[kv+22]=z01, vertices[kv+23]=tx2, vertices[kv+24]=ty1;
				vertices[kv+25]=x1, vertices[kv+26]=y1, vertices[kv+27]=z00, vertices[kv+28]=tx1, vertices[kv+29]=ty1;
			}
		}
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
static int hist[768], histmax[3], hist2[768], histmax2[3];
static float blockCR[3]={0};
double calc_entropy(int *hist, int sum)
{
	double entropy=0;
	if(sum==-1)
	{
		sum=0;
		for(int sym=0;sym<256;++sym)
			sum+=hist[sym];
	}
	if(!sum)
		return 0;
	for(int k=0;k<256;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/sum;
			p*=0xFF00;
			++p;
			p/=0x10000;
			entropy-=p*log2(p);
		}
	}
	return entropy;
}
void chart_hist_update(unsigned char *im, int iw, int ih, int x1, int x2, int y1, int y2, int *hist, int *histmax, float *CR)
{
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	int xcount=x2-x1, ycount=y2-y1, count=xcount*ycount;
	double entropy;
	if(count)
	{
		for(int kc=0;kc<3;++kc)
		{
			//histmax[kc]=0;
			for(int ky=y1;ky<y2;++ky)
			{
				for(int kx=x1;kx<x2;++kx)
				{
					unsigned char sym=im[(iw*ky+kx)<<2|kc];
					++hist[kc<<8|sym];
				}
			}
			//for(int k=0;k<res;++k)
			//{
			//	unsigned char sym=im[k<<2|kc];
			//	++hist[kc<<8|sym];
			//}
			for(int k=0;k<256;++k)
			{
				if(histmax[kc]<hist[kc<<8|k])
					histmax[kc]=hist[kc<<8|k];
			}
			if(CR)
			{
				entropy=calc_entropy(hist+(kc<<8), count);
				CR[kc]=(float)(8/entropy);
			}
		}
	}
}
void chart_dwthist_update(unsigned char *im, int iw, int ih, int kc, int kband, int x1, int x2, int y1, int y2)
{
	memset(hist+((size_t)kband<<8), 0, 256LL*sizeof(int));
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	int xcount=x2-x1, ycount=y2-y1, count=xcount*ycount;
	if(count)
	{
		for(int ky=y1;ky<y2;++ky)
		{
			for(int kx=x1;kx<x2;++kx)
			{
				unsigned char sym=im[(iw*ky+kx)<<2|kc];
				++hist[kband<<8|sym];
			}
		}
		histmax[kband]=0;
		for(int k=0;k<256;++k)
		{
			if(histmax[kband]<hist[kband<<8|k])
				histmax[kband]=hist[kband<<8|k];
		}
		double entropy=0;
		for(int k=0;k<256;++k)
		{
			int freq=hist[kband<<8|k];
			if(freq)
			{
				double p=(double)freq/count;
				p*=0xFF00;
				++p;
				p/=0x10000;
				entropy-=p*log2(p);
			}
		}
		blockCR[kband]=(float)(8/entropy);
	}
}
void chart_jointhist_update(unsigned char *im, int iw, int ih, ArrayHandle *cpuv, unsigned *gpuv, unsigned *txid)
{
	int nbits=6;

	int nlevels=1<<nbits, th=nlevels*nlevels;
	jointhistogram(image, iw, ih, nbits, &jointhist, space_not_color);
	
	if(!*txid)
		glGenTextures(1, txid);
	send_texture_pot_grey(*txid, jointhist->data, nlevels, th, 1);
	//send_texture_pot_int16x1(*txid, (unsigned*)jointhist->data, nlevels, nlevels*nlevels, 1);

	int nv=nlevels*2*3, nf=nv*5;//nlevels * 2 triangles * 3 vertices * 5 floats
	if(!*cpuv||cpuv[0]->count!=nf)
	{
		if(*cpuv)
			array_free(cpuv);
		ARRAY_ALLOC(float, *cpuv, 0, nf, 0, 0);
	}
	float *vertices=(float*)cpuv[0]->data;

	float
		x1=0, x2=-(float)nlevels,
		y1=0, y2=-(float)nlevels;
	for(int kz=0, kv=0;kz<nlevels;++kz, kv+=30)
	{
		float z=(float)kz;
		float
			tx1=0,
			tx2=(float)(nlevels-1)/nlevels,
			ty1=(float)(nlevels*kz+0)/th,
			ty2=(float)(nlevels*kz+nlevels-1)/th;
		vertices[kv   ]=x1, vertices[kv+ 1]=y1, vertices[kv+ 2]=z, vertices[kv+ 3]=tx1, vertices[kv+ 4]=ty1;
		vertices[kv+ 5]=x1, vertices[kv+ 6]=y2, vertices[kv+ 7]=z, vertices[kv+ 8]=tx1, vertices[kv+ 9]=ty2;
		vertices[kv+10]=x2, vertices[kv+11]=y2, vertices[kv+12]=z, vertices[kv+13]=tx2, vertices[kv+14]=ty2;
				
		vertices[kv+15]=x2, vertices[kv+16]=y2, vertices[kv+17]=z, vertices[kv+18]=tx2, vertices[kv+19]=ty2;
		vertices[kv+20]=x2, vertices[kv+21]=y1, vertices[kv+22]=z, vertices[kv+23]=tx2, vertices[kv+24]=ty1;
		vertices[kv+25]=x1, vertices[kv+26]=y1, vertices[kv+27]=z, vertices[kv+28]=tx1, vertices[kv+29]=ty1;
	}
	if(!*gpuv)
		glGenBuffers(1, gpuv);
	glBindBuffer(GL_ARRAY_BUFFER, *gpuv);	GL_CHECK(error);
	glBufferData(GL_ARRAY_BUFFER, nf*sizeof(float), cpuv[0]->data, GL_STATIC_DRAW);	GL_CHECK(error);
}
void update_image()
{
	if(image)
	{
		void *p2=realloc(image, (size_t)iw*ih<<2);
		if(!p2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		image=(unsigned char*)p2;
	}
	else
	{
		image=(unsigned char*)malloc((size_t)iw*ih<<2);
		if(!image)
		{
			LOG_ERROR("Allocation error");
			return;
		}
	}
	memcpy(image, im0, (size_t)iw*ih<<2);
	if(transforms)
	{
		addhalf(image, iw, ih, 3, 4);
		for(int k=0;k<(int)transforms->count;++k)
		{
			unsigned char tid=transforms->data[k];
			switch(tid)
			{
			case CT_FWD_YCoCg:		colortransform_ycocg_fwd((char*)image, iw, ih);		break;
			case CT_INV_YCoCg:		colortransform_ycocg_inv((char*)image, iw, ih);		break;
			case CT_FWD_YCoCb:		colortransform_ycocb_fwd((char*)image, iw, ih);		break;
			case CT_INV_YCoCb:		colortransform_ycocb_inv((char*)image, iw, ih);		break;
			case CT_FWD_XGZ:		colortransform_xgz_fwd((char*)image, iw, ih);		break;
			case CT_INV_XGZ:		colortransform_xgz_inv((char*)image, iw, ih);		break;
			case CT_FWD_XYZ:		colortransform_xyz_fwd((char*)image, iw, ih);		break;
			case CT_INV_XYZ:		colortransform_xyz_inv((char*)image, iw, ih);		break;
		//	case CT_FWD_EXP:		colortransform_exp_fwd((char*)image, iw, ih);		break;
		//	case CT_INV_EXP:		colortransform_exp_inv((char*)image, iw, ih);		break;
			case CT_FWD_ADAPTIVE:	colortransform_adaptive((char*)image, iw, ih, 1);	break;
			case CT_INV_ADAPTIVE:	colortransform_adaptive((char*)image, iw, ih, 0);	break;
			case CT_FWD_CUSTOM:		colortransform_custom_fwd((char*)image, iw, ih);	break;
			case CT_INV_CUSTOM:		colortransform_custom_inv((char*)image, iw, ih);	break;

		//	case ST_FWD_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 1);break;
		//	case ST_INV_JOINT:		pred_joint_apply((char*)image, iw, ih, jointpredparams, 0);break;
		//	case ST_FWD_CFL:		pred_cfl((char*)image, iw, ih, 1);					break;
		//	case ST_INV_CFL:		pred_cfl((char*)image, iw, ih, 0);					break;
		//	case ST_FWD_HYBRID3:	pred_hybrid_fwd((char*)image, iw, ih);				break;
		//	case ST_INV_HYBRID3:	pred_hybrid_inv((char*)image, iw, ih);				break;
				
			case ST_FWD_CUSTOM2:	custom2_apply((char*)image, iw, ih, 1, &c2_params);	break;
			case ST_INV_CUSTOM2:	custom2_apply((char*)image, iw, ih, 0, &c2_params);	break;
		//	case ST_FWD_CUSTOM2:	pred_custom2_apply((char*)image, iw, ih, 1);		break;
		//	case ST_INV_CUSTOM2:	pred_custom2_apply((char*)image, iw, ih, 0);		break;
		//	case ST_FWD_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 1);break;
		//	case ST_INV_LOGIC:		pred_logic_apply((char*)image, iw, ih, logic_params, 0);break;
			case ST_FWD_LEARNED:	pred_learned_v4((char*)image, iw, ih, 1);			break;
			case ST_INV_LEARNED:	pred_learned_v4((char*)image, iw, ih, 0);			break;
#ifdef ALLOW_OPENCL
			case ST_FWD_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 1);			break;
			case ST_INV_LEARNED_GPU:pred_learned_gpu((char*)image, iw, ih, 0);			break;
#endif
		//	case ST_FWD_DIFF2D:		pred_diff2d_fwd((char*)image, iw, ih, 3, 4);		break;
		//	case ST_INV_DIFF2D:		pred_diff2d_inv((char*)image, iw, ih, 3, 4);		break;
		//	case ST_FWD_HPF:		pred_hpf_fwd((char*)image, iw, ih, 3, 4);			break;
		//	case ST_INV_HPF:		pred_hpf_inv((char*)image, iw, ih, 3, 4);			break;
		//	case ST_FWD_GRAD2:		pred_grad2_fwd((char*)image, iw, ih, 3, 4);			break;
		//	case ST_INV_GRAD2:		pred_grad2_inv((char*)image, iw, ih, 3, 4);			break;
		//	case ST_FWD_ADAPTIVE:	pred_adaptive((char*)image, iw, ih, 3, 4, 1);		break;
		//	case ST_INV_ADAPTIVE:	pred_adaptive((char*)image, iw, ih, 3, 4, 0);		break;
			case ST_FWD_JXL:		pred_jxl_apply((char*)image, iw, ih, jxlparams_i16, 1);break;
			case ST_INV_JXL:		pred_jxl_apply((char*)image, iw, ih, jxlparams_i16, 0);break;
			case ST_FWD_MM:			pred_w2_apply((char*)image, iw, ih, pw2_params, 1);	break;
			case ST_INV_MM:			pred_w2_apply((char*)image, iw, ih, pw2_params, 0);	break;
			case ST_FWD_JMJ:		pred_jmj_apply((char*)image, iw, ih, 1);			break;
			case ST_INV_JMJ:		pred_jmj_apply((char*)image, iw, ih, 0);			break;
		//	case ST_FWD_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 1);			break;
		//	case ST_INV_JXL:		pred_jxl((char*)image, iw, ih, 3, 4, 0);			break;
			case ST_FWD_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 1);			break;
			case ST_INV_SORTNB:		pred_sortnb((char*)image, iw, ih, 3, 4, 0);			break;
		//	case ST_FWD_MEDIAN:		pred_median_fwd((char*)image, iw, ih, 3, 4);		break;
		//	case ST_INV_MEDIAN:		pred_median_inv((char*)image, iw, ih, 3, 4);		break;
		//	case ST_FWD_DCT3PRED:	pred_dct3_fwd((char*)image, iw, ih, 3, 4);			break;
		//	case ST_INV_DCT3PRED:	pred_dct3_inv((char*)image, iw, ih, 3, 4);			break;
		//	case ST_FWD_PATHPRED:	pred_path_fwd((char*)image, iw, ih, 3, 4);			break;
		//	case ST_INV_PATHPRED:	pred_path_inv((char*)image, iw, ih, 3, 4);			break;
			case ST_FWD_GRADPRED:	pred_grad_fwd((char*)image, iw, ih, 3, 4);			break;
			case ST_INV_GRADPRED:	pred_grad_inv((char*)image, iw, ih, 3, 4);			break;
		//	case ST_FWD_BITWISE:	pred_bitwise((char*)image, iw, ih, 1);				break;
		//	case ST_INV_BITWISE:	pred_bitwise((char*)image, iw, ih, 0);				break;
			case ST_FWD_CUSTOM:		pred_custom_apply((char*)image, iw, ih, 1, customparam_st);break;
			case ST_INV_CUSTOM:		pred_custom_apply((char*)image, iw, ih, 0, customparam_st);break;
			case ST_FWD_DCT4:		image_dct4_fwd((char*)image, iw, ih);				break;
			case ST_INV_DCT4:		image_dct4_inv((char*)image, iw, ih);				break;
		//	case ST_FWD_DCT8:		image_dct8_fwd((char*)image, iw, ih);				break;
		//	case ST_INV_DCT8:		image_dct8_inv((char*)image, iw, ih);				break;
			case ST_FWD_SHUFFLE:	shuffle((char*)image, iw, ih, 1);					break;
			case ST_INV_SHUFFLE:	shuffle((char*)image, iw, ih, 0);					break;
		//	case ST_FWD_SPLIT:		image_split_fwd((char*)image, iw, ih);				break;
		//	case ST_INV_SPLIT:		image_split_inv((char*)image, iw, ih);				break;

			case ST_FWD_LAZY:
			case ST_INV_LAZY:
			case ST_FWD_HAAR:
			case ST_INV_HAAR:
			case ST_FWD_SQUEEZE:
			case ST_INV_SQUEEZE:
			case ST_FWD_LEGALL53:
			case ST_INV_LEGALL53:
		//	case ST_FWD_CDF97:
		//	case ST_INV_CDF97:
		//	case ST_FWD_GRAD_DWT:
		//	case ST_INV_GRAD_DWT:
			case ST_FWD_EXPDWT:
			case ST_INV_EXPDWT:
			case ST_FWD_CUSTOM_DWT:
			case ST_INV_CUSTOM_DWT:
				{
					ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
					char *temp=(char*)malloc(MAXVAR(iw, ih));
					for(int kc=0;kc<3;++kc)
					{
						switch(tid)
						{
						case ST_FWD_LAZY:      dwt2d_lazy_fwd   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_LAZY:      dwt2d_lazy_inv   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_HAAR:      dwt2d_haar_fwd   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_HAAR:      dwt2d_haar_inv   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_SQUEEZE:   dwt2d_squeeze_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_SQUEEZE:   dwt2d_squeeze_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_LEGALL53:  dwt2d_cdf53_fwd  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_LEGALL53:  dwt2d_cdf53_inv  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					//	case ST_FWD_CDF97:     dwt2d_cdf97_fwd  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					//	case ST_INV_CDF97:     dwt2d_cdf97_inv  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					//	case ST_FWD_GRAD_DWT:  dwt2d_grad_fwd   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
					//	case ST_INV_GRAD_DWT:  dwt2d_grad_inv   ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_EXPDWT:    dwt2d_exp_fwd    ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;//TODO use customdwtparams instead of sharing allcustomparam_st
						case ST_INV_EXPDWT:    dwt2d_exp_inv    ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
						case ST_FWD_CUSTOM_DWT:dwt2d_custom_fwd ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
						case ST_INV_CUSTOM_DWT:dwt2d_custom_inv ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp, customparam_st);break;
						}
					}
					array_free(&sizes);
					free(temp);
				}
				break;
		//	case ST_FWD_DEC_DWT:   dwt2d_dec_fwd((char*)image, iw, ih);	break;
		//	case ST_INV_DEC_DWT:   dwt2d_dec_inv((char*)image, iw, ih);	break;
			}
		}//for
		addhalf(image, iw, ih, 3, 4);
	}//if transforms

	channel_entropy(image, iw*ih, 3, 4, ch_cr, usage);
	
	combCRhist[combCRhist_idx][0]=ch_cr[0];
	combCRhist[combCRhist_idx][1]=ch_cr[1];
	combCRhist[combCRhist_idx][2]=ch_cr[2];
	combCRhist[combCRhist_idx][3]=3/(1/ch_cr[0]+1/ch_cr[1]+1/ch_cr[2]);
	for(int k=0;k<4;++k)
	{
		if(combCRhist_max==1||combCRhist_max<combCRhist[combCRhist_idx][k])
			combCRhist_max=combCRhist[combCRhist_idx][k];
	}
	combCRhist_idx=(combCRhist_idx+1)%combCRhist_SIZE;

	if(!send_image_separate_subpixels(image, iw, ih, image_txid))
		LOG_ERROR("Failed to send texture to GPU");

	switch(mode)
	{
	case VIS_PLANES:
		chart_planes_update(image, iw, ih, &cpu_vertices, &gpu_vertices);
		break;
	case VIS_MESH:
		chart_mesh_update(image, iw, ih, &cpu_vertices, &gpu_vertices);
		break;
	case VIS_MESH_SEPARATE:
		chart_mesh_sep_update(image, iw, ih, &cpu_vertices, &gpu_vertices);
		break;
	case VIS_HISTOGRAM:
		memset(histmax, 0, sizeof(histmax));
		memset(hist, 0, sizeof(hist));
		chart_hist_update(image, iw, ih, 0, iw, 0, ih, hist, histmax, 0);
		break;
	case VIS_JOINT_HISTOGRAM:
		chart_jointhist_update(image, iw, ih, &cpu_vertices, &gpu_vertices, image_txid+2);
		break;
	}
}

void draw_AAcuboid_wire(float x1, float x2, float y1, float y2, float z1, float z2, int color)
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
void chart_planes_draw()
{
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, 0, pixel_amplitude, 0xFF000000);

	draw_3D_triangles(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[0]);
}
void chart_mesh_draw()
{
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, 0, pixel_amplitude, 0xFF000000);

	draw_3D_triangles(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[0]);
}
void chart_mesh_sep_draw()
{
	float RMSE1=0, RMSE2=0;
	int RMSE_den=0;
	float reach=5;
	float ix=0, iy=0;
	if(extrainfo)
	{
		float *vertices=(float*)cpu_vertices->data;
		ix=iw-1-cam.x, iy=cam.y;
		float xstart=ix-reach, xend=ix+reach, ystart=iy-reach, yend=iy+reach;
		for(int ky=0, kv=0;ky<ih-1;++ky)
		{
			float error[3]={0};
			for(int kx=0;kx<iw-1;++kx)
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
						int kx2=kx+1, ky2=ky+1, idx=iw*ky2+kx2;
						char
							ctltl    =kx2-2>=0&&ky2-2>=0?image[(idx-iw*2-2)<<2|kc]-128:0,
							ctt      =kx2  <iw&&ky2-2>=0?image[(idx-iw*2  )<<2|kc]-128:0,
							ctrtr    =kx2+2<iw&&ky2-2>=0?image[(idx-iw*2+2)<<2|kc]-128:0,

							ctopleft =                   image[(idx-iw  -1)<<2|kc]-128  ,
							ctop     =kx2  <iw          ?image[(idx-iw    )<<2|kc]-128:0,
							ctopright=kx2+1<iw          ?image[(idx-iw  +1)<<2|kc]-128:0,

							cll      =kx2-2>=0&&ky2  <ih?image[(idx     -2)<<2|kc]-128:0,
							cleft    =          ky2  <ih?image[(idx     -1)<<2|kc]-128:0,
							ccurr    =kx2  <iw&&ky2  <ih?image[ idx        <<2|kc]-128:0;
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
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, 0, pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, mesh_separation, mesh_separation+pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, mesh_separation*2, mesh_separation*2+pixel_amplitude, 0xFF000000);

	draw_3D_triangles(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[1]);
	
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
void chart_hist_draw(float x1, float x2, float y1, float y2, int cstart, int cend, int color, unsigned char alpha, int *hist, int *histmax)
{
	for(int kc=cstart;kc<cend;++kc)
	{
		if(histmax[kc])
		{
			float dy=(y2-y1)/3.f, histpx=dy/histmax[kc];
			int k=1;
			float y=k*histpx*10000;
			for(;y<dy;++k)
			{
				draw_line(x1, y1+(kc+1)*dy-y, x2, y1+(kc+1)*dy-y, color?color:alpha<<24|0xFF<<(kc<<3));//0x40
				y=k*histpx*10000;
			}
			for(int k=0;k<256;++k)
				draw_rect(x1+k*(x2-x1)/256, x1+(k+1)*(x2-x1)/256, y1+(kc+1)*dy-hist[kc<<8|k]*histpx, y1+(kc+1)*dy, color?color:alpha<<24|0xFF<<(kc<<3));//0x80
		}
	}
}
void chart_hist_draw2(float x1, float x2, float y1, float y2, int color, int *hist, int histmax)
{
	if(histmax==-1)
	{
		histmax=0;
		for(int sym=0;sym<256;++sym)
		{
			if(histmax<hist[sym])
				histmax=hist[sym];
		}
	}
	if(!histmax)
		return;
	float dy=y2-y1, histpx=dy/histmax;
	int k=1;
	float y=k*histpx*10000;
	for(;y<dy;++k)
	{
		draw_line(x1, y2-y, x2, y2-y, color);
		y=k*histpx*10000;
	}
	for(int k=0;k<256;++k)
		draw_rect(x1+k*(x2-x1)/256, x1+(k+1)*(x2-x1)/256, y2-hist[k]*histpx, y2, color);
}
void chart_jointhist_draw()
{
	draw_AAcuboid_wire(0, 64, 0, 64, 0, 64, 0xFF000000);
	draw_contour3d_rect(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[2], 0.8f);
}
void e24_update()
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

		e24_estimate(image, iw, ih, (int)roundf(x1), (int)roundf(x2), (int)roundf(y1), (int)roundf(y2));
	}
}

#if 0
void applycolortransform(int tidx)
{
	addhalf(image, iw, ih, 3, 4);
	switch(tidx)
	{
	case CT_NONE:break;
	case CT_YCoCg:
		colortransform_ycocg_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_YCoCgT:
		colortransform_ycocb_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_XGZ:
		colortransform_xgz_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_XYZ:
		colortransform_xyz_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_EXP:
		colortransform_exp_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_LEARNED:
		colortransform_learned_fwd((unsigned char*)image, iw, ih);
		break;
	case CT_CUSTOM:
		colortransform_custom_fwd((unsigned char*)image, iw, ih);
		break;
	}
	addhalf(image, iw, ih, 3, 4);
}
void undocolortransform(int tidx)
{
	addhalf(image, iw, ih, 3, 4);
	switch(tidx)
	{
	case CT_NONE:break;
	case CT_YCoCg:
		colortransform_ycocg_inv((unsigned char*)image, iw, ih);
		break;
	case CT_YCoCgT:
		colortransform_ycocb_inv((unsigned char*)image, iw, ih);
		break;
	case CT_XGZ:
		colortransform_xgz_inv((unsigned char*)image, iw, ih);
		break;
	case CT_XYZ:
		colortransform_xyz_inv((unsigned char*)image, iw, ih);
		break;
	case CT_EXP:
		colortransform_exp_inv((unsigned char*)image, iw, ih);
		break;
	case CT_LEARNED:
		colortransform_learned_inv((unsigned char*)image, iw, ih);
		break;
	case CT_CUSTOM:
		colortransform_custom_inv((unsigned char*)image, iw, ih);
		break;
	}
	addhalf(image, iw, ih, 3, 4);
}
void applyspatialtransform(int tidx)
{
	addhalf(image, iw, ih, 3, 4);
	switch(tidx)
	{
	case ST_NONE:
		break;
	case ST_DIFF2D:
		image_differentiate(image, iw, ih, 3, 4);
		break;
	case ST_UNPLANE:
		image_unplane(image, iw, ih, 3, 4);
		break;
	case ST_LAZY:
	case ST_HAAR:
	case ST_SQUEEZE:
	case ST_CDF53:
	case ST_CDF97:
		{
			ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
			char *temp=(char*)malloc(MAXVAR(iw, ih));
			for(int kc=0;kc<3;++kc)
			{
				switch(tidx)
				{
				case ST_LAZY:
					dwt2d_lazy_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_HAAR:
					dwt2d_haar_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_SQUEEZE:
					dwt2d_squeeze_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_CDF53:
					dwt2d_cdf53_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_CDF97:
					dwt2d_cdf97_fwd((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				}
			}
			array_free(&sizes);
			free(temp);
		}
		break;
	case ST_CUSTOM:
		image_customst_fwd(image, iw, ih, 3, 4);
		break;
	}
	addhalf(image, iw, ih, 3, 4);
}
void undospatialtransform(int tidx)
{
	addhalf(image, iw, ih, 3, 4);
	switch(tidx)
	{
	case ST_NONE:
		break;
	case ST_DIFF2D:
		image_integrate(image, iw, ih, 3, 4);
		break;
	case ST_UNPLANE:
		image_replane(image, iw, ih, 3, 4);
		break;
	case ST_LAZY:
	case ST_HAAR:
	case ST_SQUEEZE:
	case ST_CDF53:
	case ST_CDF97:
		{
			ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
			char *temp=(char*)malloc(MAXVAR(iw, ih));
			for(int kc=0;kc<3;++kc)
			{
				switch(tidx)
				{
				case ST_LAZY:
					dwt2d_lazy_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_HAAR:
					dwt2d_haar_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_SQUEEZE:
					dwt2d_squeeze_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_CDF53:
					dwt2d_cdf53_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				case ST_CDF97:
					dwt2d_cdf97_inv((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
					break;
				}
			}
			array_free(&sizes);
			free(temp);
		}
		break;
	case ST_CUSTOM:
		image_customst_inv(image, iw, ih, 3, 4);
		break;
	}
	addhalf(image, iw, ih, 3, 4);
}
#endif


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
AABB buttons[6]={0};//0: CT,  1: ST,  2: list of transforms,  3: clamp bounds,  4: learning rate

int io_init(int argc, char **argv)//return false to abort
{
	set_window_title("pxView3D");
	glClearColor(1, 1, 1, 1);

	cam_zoomIn(cam, 1);
	cam_turnMouse(cam, 0, 0, mouse_sensitivity);
	memcpy(&cam0, &cam, sizeof(cam));

	set_bk_color(0xA0808080);
	return 1;
}
void io_resize()
{
	AABB *p=buttons;
	float xstep=tdx*guizoom, ystep=tdy*guizoom;
	p->x1=xstep*2, p->x2=p->x1+xstep*30, p->y1=(float)(h>>1), p->y2=p->y1+customparam_ct_h*ystep, ++p;//0: color params - left
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*54, p->y1=(float)((h>>1)+(h>>2)), p->y2=p->y1+ystep*3, ++p;//1: spatial params
	p->x1=(float)(w-200), p->x2=(float)w, p->y1=tdy*2, p->y2=p->y1+tdy*T_COUNT, ++p;//2: transforms list
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*14, p->y1=(float)((h>>1)+(h>>2))+ystep*4, p->y2=p->y1+ystep, ++p;//3: clamp bounds
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*21, p->y1=(float)((h>>1)+(h>>2))+ystep*5, p->y2=p->y1+ystep, ++p;//4: learning rate
	p->x1=(float)(w>>2), p->x2=p->x1+tdx*11*6, p->y1=(float)((h>>1)+(h>>2)), p->y2=p->y1+tdy*3;//5: jxl params
}
int io_mousemove()//return true to redraw
{
	if(mode==VIS_IMAGE_BLOCK
		//||mode==VIS_IMAGE_E24
		||mode==VIS_DWT_BLOCK)
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
	else if(drag)
	{
		int X0=w>>1, Y0=h>>1;
		cam_turnMouse(cam, mx-X0, my-Y0, mouse_sensitivity);
		set_mouse(X0, Y0);
		return !timer;
	}
	return 0;
}
void click_hittest(int mx, int my, int *objidx, int *cellx, int *celly, int *cellidx, AABB **p)
{
	*p=buttons;
	*objidx=0;
	for(*objidx=0;*objidx<COUNTOF(buttons);++*objidx, ++*p)
	{
		if(!transforms_customenabled&&(*p-buttons==0||*p-buttons==1||*p-buttons==3||*p-buttons==4))//when these buttons are inactive they shouldn't block the click
			continue;
		if(mx>=p[0]->x1&&mx<p[0]->x2&&my>=p[0]->y1&&my<p[0]->y2)
			break;
	}
	switch(*objidx)
	{
	case 0://color transform
		*cellx=mx-p[0]->x1>(p[0]->x2-p[0]->x1)*0.5f;
		*celly=(int)floorf((my-p[0]->y1)*customparam_ct_h/(p[0]->y2-p[0]->y1));
		*cellidx=customparam_ct_w**celly+*cellx;
		break;
	case 1://spatial transform
		*cellx=(int)floorf((mx-p[0]->x1)*(customparam_st_reach<<1|1)/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((my-p[0]->y1)*(customparam_st_reach+1)/(p[0]->y2-p[0]->y1));
		*cellidx=(customparam_st_reach<<1|1)**celly+*cellx;
		break;
	case 2://list of transforms
		*cellx=0;
		*celly=(int)floorf((my-p[0]->y1)*T_COUNT/(p[0]->y2-p[0]->y1))+1;//zero idx is T_NONE
		*cellidx=*celly;
		break;
	case 3://clamp bounds
		*cellx=mx-p[0]->x1>(p[0]->x2-p[0]->x1)*0.5f;
		*celly=0;
		*cellidx=*cellx;
		break;
	case 4://learning rate
		*cellx=0;
		*celly=0;
		*cellidx=*cellx;
		break;
	case 5://jxl params
		*cellx=(int)floorf((mx-p[0]->x1)*11/(p[0]->x2-p[0]->x1));
		*celly=(int)floorf((my-p[0]->y1)* 3/(p[0]->y2-p[0]->y1));
		*cellidx=11**celly+*cellx;
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
	if(image&&(transforms_customenabled||transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL]))//change custom transform params
	{
		int objidx=0, cellx=0, celly=0, cellidx=0;
		AABB *p=buttons;
		click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &p);
		if(objidx!=-1)
		{
			int sign=(forward>0)-(forward<0);//abs(forward) is 120
			int ch=(int)floorf((mx-p->x1)/(guizoom*tdx)), ch2;
			switch(objidx)
			{
			case 0://color transform

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
				break;
			case 1://spatial transform

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
				break;
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
			case 5://jxl params
				ch=(int)floorf((mx-p->x1)/tdx);
				MODVAR(ch2, ch, 6);
				if(ch2>=2&&ch2<6)
				{
					int delta=sign<<((5-ch2)<<2);
					jxlparams_i16[cellidx]+=delta;
				}
				else if(cellx>=4)
					jxlparams_i16[cellidx]=-jxlparams_i16[cellidx];
				break;
			}
			update_image();
		}
		else
			goto normal_operation;
	}
	else
	{
	normal_operation:
		if(mode==VIS_IMAGE_BLOCK
			//||mode==VIS_IMAGE_E24
			||mode==VIS_DWT_BLOCK)
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
		else if(GET_KEY_STATE(KEY_SHIFT))//shift wheel		change cam speed
		{
				 if(forward>0)	cam.move_speed*=2;
			else				cam.move_speed*=0.5f;
		}
		else
		{
				 if(forward>0)	cam_zoomIn(cam, 1.1f);
			else				cam_zoomOut(cam, 1.1f);
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
int parse_nvals(ArrayHandle text, int idx, short *params, int count)
{
	int k;

	k=0;
	while(k<count)
	{
		for(;idx<text->count&&isspace(text->data[idx]);++idx);

		int neg=text->data[idx]=='-';
		idx+=neg;//skip sign
		if(text->data[idx]=='0'&&(text->data[idx]&0xDF)=='X')//skip hex prefix
			idx+=2;
		char *end=text->data+idx;
		params[k]=(int)strtol(text->data+idx, &end, 16);
		idx=(int)(end-text->data);
		if(neg)
			params[k]=-params[k];

		for(;idx<text->count&&!isspace(text->data[idx])&&text->data[idx]!='-'&&!isdigit(text->data[idx]);++idx);//skip comma

		++k;

		if(idx>=text->count&&k<count)
			return idx;
	}
	return idx;
}
static void append_i16_row(ArrayHandle *str, short *vals, int count)
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
	//switch(key)
	//{
	//case 'S':
	//	if(keyboard[KEY_CTRL])
	//	{
	//	//	savedfile=dialog_save_file(file_filters, SIZEOF(file_filters), initialname);
	//		return 1;
	//	}
	//	break;
	//}
#if 0
	if(mode==VIS_IMAGE_E24)
	{
		switch(key)
		{
		case KEY_LEFT:
			blockmx-=blocksize;
			//e24_update();
			return 1;
		case KEY_RIGHT:
			blockmx+=blocksize;
			//e24_update();
			return 1;
		case KEY_UP:
			blockmy-=blocksize;
			//e24_update();
			return 1;
		case KEY_DOWN:
			blockmy+=blocksize;
			//e24_update();
			return 1;
		}
	}
	else
#endif
	if(transforms_customenabled)
	{
		switch(key)
		{
#if 0
		case KEY_LBUTTON:
			{
				int objidx=0, cellx=0, celly=0, cellidx=0;
				AABB *p;
				click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &p);
				if(objidx==1)
				{
					if(!losshist)
						ARRAY_ALLOC(double, losshist, 0, 1000, 0, 0);
					double *p2=(double*)losshist->data;
					profile_idx=cellidx;
					double temp=customparam_st[12*customparam_ch_idx+cellidx];
					for(int k=0;k<1000;++k)
					{
						customparam_st[12*customparam_ch_idx+cellidx]=customparam_ct[8]+(customparam_ct[9]-customparam_ct[8])*k/1000;
						p2[k]=opt_causal_reach2(image, iw, ih, 1, customparam_st+12*customparam_ch_idx, customparam_ct+11, g_lr, 1);
					}
					customparam_st[12*customparam_ch_idx+cellidx]=temp;
					//av_rmse=opt_causal_reach2(image, iw, ih, 1, customparam_st, customparam_ct+11, 1e-10, 1);
						
					for(int k=0;k<1000;++k)
					{
						if(!k||minloss>p2[k])
							minloss=p2[k];
						if(!k||maxloss<p2[k])
							maxloss=p2[k];
					}
					//update_image();
					return 1;
				}
			}
			break;
#endif
		case KEY_UP:
		case KEY_DOWN:
		case KEY_LEFT:
		case KEY_RIGHT:
			{
				const int idx_limit=COUNTOF(customparam_ct)+COUNTOF(customparam_st)/6;
				const int st_w=customparam_st_reach<<1|1, st_h=customparam_st_reach+1, st_yoffset=customparam_ct_h-st_h,
					w_total=customparam_ct_w+st_w, h_total=customparam_ct_h;
				int x, y;
				int dx=0, dy=0;
				switch(key)
				{
				case KEY_UP:	dy=-1;	break;
				case KEY_DOWN:	dy= 1;	break;
				case KEY_LEFT:	dx=-1;	break;
				case KEY_RIGHT:	dx= 1;	break;
				}
				MODVAR(customparam_sel, customparam_sel, idx_limit);
				if(customparam_sel<COUNTOF(customparam_ct))
					x=customparam_sel%customparam_ct_w, y=customparam_sel/customparam_ct_w;
				else
					x=customparam_ct_w+(customparam_sel-COUNTOF(customparam_ct))%st_w, y=st_yoffset+(customparam_sel-COUNTOF(customparam_ct))/st_w;
				x+=dx;
				y+=dy;

				if(dx)
				{
					if(dx>0&&y==2&&x==2)
						++y;
					else if(dx>0&&x==4&&y==5)
						--y;
					else if(y<st_yoffset)
						MODVAR(x, x, customparam_ct_w);
					else if(y<h_total-1)
						MODVAR(x, x, w_total);
					else
						MODVAR(x, x, customparam_ct_w+customparam_st_reach);
				}
				else if(dy)
				{
					if(dy>0&&x==4&&y==5)
						--x;
					else if(dy<0&&x==2&&y==2)
						--x;
					else if(x<customparam_ct_w)
						MODVAR(y, y, h_total);
					else if(x<customparam_ct_w+customparam_st_reach)
					{
						y-=st_yoffset;
						MODVAR(y, y, h_total-st_yoffset);
						y+=st_yoffset;
					}
					else
					{
						y-=st_yoffset;
						MODVAR(y, y, h_total-st_yoffset-1);
						y+=st_yoffset;
					}
				}
#if 0
				if(!(BETWEEN_EXC(0, x, customparam_ct_w)&&BETWEEN_EXC(0, y, customparam_ct_h)||BETWEEN_EXC(customparam_ct_w, x, w_total)&&BETWEEN_EXC(st_yoffset, y, h_total)))//bring sel inside
				{
					MODVAR(x, x, w_total);//bring sel inside AABB
					MODVAR(y, y, h_total);


					//if(y<st_yoffset)
					//	MODVAR(x, x, customparam_ct_w);
					//if(x>=customparam_ct_w+customparam_st_reach)
					//{
					//	y-=st_yoffset;
					//	MODVAR(y, y, customparam_st_reach);
					//	y+=st_yoffset;
					//}

					//if(x>=customparam_ct_w)//clamp to region edges
					//{
					//	if(y<st_yoffset)
					//		y=st_yoffset;
					//	if(x>=customparam_ct_w+customparam_st_reach)
					//	{
					//		if(y>st_yoffset+customparam_st_reach-1)
					//			y=st_yoffset+customparam_st_reach-1;
					//	}
					//}
				}
#endif
				if(x<customparam_ct_w)
					customparam_sel=customparam_ct_w*y+x;//color transform
				else
					customparam_sel=COUNTOF(customparam_ct)+st_w*(y-st_yoffset)+x-customparam_ct_w;//spatial filter
				MODVAR(customparam_sel, customparam_sel, idx_limit);
				//if((size_t)customparam_sel>=(size_t)idx_limit)//assertion
				//	LOG_ERROR("Index error, idx %d", customparam_sel);
			}
			return 1;
		}
	}
	switch(key)
	{
	case KEY_LBUTTON:
	case KEY_RBUTTON:
		if(image)
		{
			int objidx=-1, cellx=0, celly=0, cellidx=0;
			AABB *p=buttons;
			click_hittest(mx, my, &objidx, &cellx, &celly, &cellidx, &p);
			switch(objidx)
			{
			case 0://color transform params
				if(key==KEY_RBUTTON)
				{
					customparam_ct[cellidx]=0;
					update_image();
					return 1;
				}
				break;
			case 1://spatial transforms params
				if(key==KEY_RBUTTON)
				{
					customparam_st[12*customparam_ch_idx+cellidx]=0;
					update_image();
					return 1;
				}
				break;
			case 2://transform list
				{
					if(BETWEEN_EXC(1, cellidx, T_COUNT))
					{
						if(key==KEY_LBUTTON)
						{
							//if(GET_KEY_STATE(KEY_CTRL))
							//	transforms_removeall();
							transforms_append(cellidx);
						}
						else
							transforms_removebyid(cellidx);
						update_image();
						return 1;
					}
				}
				break;
			case 3://clamp bounds
				if(key==KEY_RBUTTON)
				{
					if(cellidx)
						customparam_clamp[1]=127;
					else
						customparam_clamp[0]=-128;
					update_image();
					return 1;
				}
				break;
			case 4://spatial transforms params
				if(key==KEY_RBUTTON)
				{
					g_lr=1e-10;
					return 1;
				}
				break;
			case 5://jxl params
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
			default:
				if(key==KEY_LBUTTON)
					goto toggle_drag;
				break;
			}
		}
		else
			goto toggle_drag;
		break;
	case KEY_ESC:
toggle_drag:
		if(mode==VIS_IMAGE_BLOCK
			//||mode==VIS_IMAGE_E24
			||mode==VIS_DWT_BLOCK)
		{
			if(key==KEY_LBUTTON)
			{
				blockmx=mx, blockmy=my;
				//e24_update();
				return 1;
			}
		}
		else
		{
			show_mouse(drag);
			drag=!drag;
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

#define		AK(KEY)		case KEY:
	ACTIVE_KEY_LIST
#undef		AK
		timer_start(10, TIMER_ID_KEYBOARD);
		break;
		
	case KEY_F1:
		messagebox(MBOX_OK, "Controls",
			"WASDTG:\tMove cam\n"
			"Arrow keys:\tTurn cam\n"
			"Mouse1:\t\tToggle mouse look\n"
			"Ctrl O:\t\tOpen image\n"
			"\n"
			"R:\t\tReset cam\n"
			"Ctrl R:\t\tDisable all transforms\n"
			"Ctrl E:\t\tReset custom transform parameters\n"
			"[]:\t\tSwitch channel (custom transform)\n"
			"H:\t\tReset CR history graph\n"
			"Ctrl C:\t\tCopy data\n"
			"Ctrl V:\t\tPaste data\n"
			"C:\t\tToggle joint histogram type / fill screen in image view\n"
			"\n"
			"Mouse1/Mouse2:\tAdd/remove transforms to the list\n"
			"Ctrl Mouse1:\tReplace all transforms of this type\n"
			"\n"
			"M / Shift M:\tCycles between:\n"
			"\t1: 3D View: Levels\n"
			"\t2: 3D View: Mesh\n"
			"\t3: 3D View: Mesh (separate channels)\n"
			"\t4: Image tricolor view\n"
			"\t5: Image view\n"
			"\t6: Image block histogram\n"
		//	"\t7: Optimized block compression estimate (E24)\n"
			"\t7: DWT block histogram\n"
			"\t8: Histogram\n"
			"\t9: Joint histogram\n"
		);
		//prof_on=!prof_on;
		return 0;
	case 'R':
		if(GET_KEY_STATE(KEY_CTRL)&&image)
		{
			transforms_removeall();
			update_image();
		}
		else
			memcpy(&cam, &cam0, sizeof(cam));
		return 1;
	case 'E':
		if(image&&GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)//reset params
		{
			customtransforms_resetparams();
			update_image();
			return 1;
		}
	//	wireframe=!wireframe;
		break;
	case 'O':
		if(GET_KEY_STATE(KEY_CTRL))
		{
			ArrayHandle fn2=dialog_open_file(0, 0, 0);
			if(fn2)
			{
				if(fn)
					array_free(&fn);
				fn=fn2;
				if(im0)
					free(im0);
				im0=stbi_load((char*)fn->data, &iw, &ih, &nch0, 4);
				//array_free(&fn);
				if(im0)
				{
					filesize=get_filesize((char*)fn->data);
					set_window_title("%s - pxView3D", (char*)fn->data);
					//color_transform=CT_NONE;
					//spatialtransform=CT_NONE;
					//array_free(&transforms);
					update_image();
				}
			}
		}
		return 1;
	case 'C':
		if(image&&GET_KEY_STATE(KEY_CTRL))//copy custom transform value
		{
			ArrayHandle str;
			STR_ALLOC(str, 0);
			if(mode==VIS_IMAGE)
			{
				float cr_combined=3/(1/ch_cr[0]+1/ch_cr[1]+1/ch_cr[2]);
				str_append(&str, "T %f\tR %f\tG %f\tB %f\tJ %f", cr_combined, ch_cr[0], ch_cr[1], ch_cr[2], ch_cr[3]);
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
#endif
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
			else if(transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL])
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
				for(int ky=0;ky<customparam_ct_h;++ky)
				{
					for(int kx=0;kx<customparam_ct_w;++kx)
						str_append(&str, "\t%g", customparam_ct[customparam_ct_w*ky+kx]);
					str_append(&str, "\n");
				}
				const int stw=customparam_st_reach<<1|1;
				for(int kc2=0;kc2<6;++kc2)
				{
					const int np=_countof(customparam_st)/6;
					const double *params=customparam_st+12*kc2;
					for(int k=0;k<np;++k)
					{
						int x=k%stw, y=k/stw;
						if(shift)
						{
							int val=(int)(params[k]*0x1000);//fixed 3.12 bit
							str_append(&str, "%c0x%04X,%c", val<0?'-':' ', abs(val), x<stw-1&&k<np-1?'\t':'\n');
						}
						else
							str_append(&str, "%g%c", params[k], x<stw-1&&k<np-1?'\t':'\n');
					}
				}
				for(int k=0;k<COUNTOF(customparam_clamp);++k)
					str_append(&str, "%d%c", customparam_clamp[k], k<COUNTOF(customparam_clamp)-1?'\t':'\n');
			}
			copy_to_clipboard((char*)str->data, (int)str->count);
			array_free(&str);
		}
		else
		{
			if(mode==VIS_IMAGE)
				show_full_image=!show_full_image;
			else if(mode==VIS_JOINT_HISTOGRAM)
			{
				int shift=GET_KEY_STATE(KEY_SHIFT);
				space_not_color+=1-(shift<<1);
				MODVAR(space_not_color, space_not_color, 4);
				update_image();
			}
		}
		return 1;
	case 'V':
		if(image&&GET_KEY_STATE(KEY_CTRL))//paste custom transform value
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
					parse_nvals(text, idx, logic_params, COUNTOF(logic_params));
				}
				else
#endif
				if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
				{
					parse_nvals(text, idx, (short*)&c2_params, _countof(c2_params.c0)+_countof(c2_params.c1)+_countof(c2_params.c2));
				}
				else if(transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL])
				{
					parse_nvals(text, idx, jxlparams_i16, COUNTOF(jxlparams_i16));
				}
				else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
				{
					parse_nvals(text, idx, pw2_params, COUNTOF(pw2_params));
				}
#if 0
				else if(transforms_mask[ST_FWD_JOINT]||transforms_mask[ST_INV_JOINT])
				{
					k=0, kend=COUNTOF(jointpredparams);
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
					k=0, kend=COUNTOF(customparam_hybrid), idx=0;
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
					k=0, kend=COUNTOF(customparam_ct);
					for(;k<kend;++k)
					{
						for(;idx<text->count&&isspace(text->data[idx]);++idx);
						customparam_ct[k]=atof(text->data+idx);
						for(;idx<text->count&&!isspace(text->data[idx]);++idx);
						if(idx>=text->count)
							goto paste_finish;
					}
					k=0;
					kend=COUNTOF(customparam_st);
					if(shift)
					{
						short params[COUNTOF(customparam_st)];
						idx=parse_nvals(text, idx, params, COUNTOF(params));
						for(int k2=0;k2<COUNTOF(customparam_st);++k2)
							customparam_st[k2]=(double)params[k2]/0x1000;
					}
					else
					{
						for(;k<kend;++k)
						{
							for(;idx<text->count&&isspace(text->data[idx]);++idx);
							customparam_st[k]=atof(text->data+idx);
							for(;idx<text->count&&!isspace(text->data[idx]);++idx);
							if(idx>=text->count)
								goto paste_finish;
						}
					}
					k=0;
					kend=COUNTOF(customparam_clamp);
					for(;k<kend;++k)
					{
						for(;idx<text->count&&isspace(text->data[idx]);++idx);
						customparam_clamp[k]=atoi(text->data+idx);
						for(;idx<text->count&&!isspace(text->data[idx]);++idx);
						if(idx>=text->count)
							goto paste_finish;
					}
				}

			paste_finish:
				array_free(&text);
				update_image();
				return 1;
			}
		}
		break;
	case 'M':
		if(image)
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
	case 'X':
		{
			extrainfo=!extrainfo;
			return 1;
		}
		break;
	case 'Z'://TODO show neighbor pixels around cursor
		break;
	case KEY_LBRACKET:
	case KEY_RBRACKET:
		if(transforms_customenabled
			//||transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC]
			)
		{
			//int logic=transforms_mask[ST_FWD_LOGIC]||transforms_mask[ST_INV_LOGIC];
			//const int nval=_countof(customparam_st);
			//memcpy(allcustomparam_st+nval*customparam_ch_idx, customparam_st, sizeof(customparam_st));//save
			customparam_ch_idx+=((key==KEY_RBRACKET)<<1)-1;
			//customparam_ch_idx+=(((key==KEY_RBRACKET)<<1)-1)<<logic;
			MODVAR(customparam_ch_idx, customparam_ch_idx, 6);
			//customparam_ch_idx&=~logic;
			//memcpy(customparam_st, allcustomparam_st+nval*customparam_ch_idx, sizeof(customparam_st));//load
			return 1;
		}
		break;
	case KEY_SPACE:
	//case 'B':
	//case 'N':
		if(image)
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
				colortransform_ycocb_fwd(buf2, iw, ih);
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
			if(transforms_customenabled)
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
				colortransform_ycocb_fwd(buf2, iw, ih);
				opt_cr2_v2(buf2, iw, ih, customparam_ch_idx/2);
				free(buf2);
				update_image();
			}
			//	timer_start(50);
			else if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
			{
				if(GET_KEY_STATE(KEY_CTRL))
					memset(&c2_params, 0, sizeof(c2_params));
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
					//colortransform_ycocb_fwd(buf2, iw, ih);
					custom2_opt(buf2, iw, ih, &c2_params);
					free(buf2);
				}
				update_image();
			}
			else if(transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL]||transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM]
			//	||transforms_mask[ST_FWD_JOINT]||transforms_mask[ST_INV_JOINT]
			)
			{
				int jxl=transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL];
				int pw2=transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM];
				if(GET_KEY_STATE(KEY_CTRL))
				{
					if(jxl)
					{
						for(int k=0;k<33;++k)
							jxlparams_i16[k]=(k%11<4)<<12;
					}
					else
					{
						memset(jointpredparams, 0, sizeof(jointpredparams));
						//jointpredparams[0]=0x1000;
						//jointpredparams[25]=0x1000;
						//jointpredparams[50]=0x1000;
						//for(int k=0;k<12;++k)
						//	jointpredparams[72+k]=0x1000;
					}
				}
				else
				{
					//static int press=0;
					int res=iw*ih;
					char *buf2=(char*)malloc((size_t)res<<2);
					char *buf3=(char*)malloc((size_t)res<<2);
					if(!buf2||!buf3)
					{
						LOG_ERROR("Allocation error");
						return 0;
					}
					memcpy(buf2, im0, (size_t)res<<2);
					addhalf((unsigned char*)buf2, iw, ih, 3, 4);

					int step=256;
						 if(keyboard['1'])step=128;
					else if(keyboard['2'])step= 64;
					else if(keyboard['3'])step= 32;
					else if(keyboard['4'])step= 16;
					else if(keyboard['5'])step=  8;
					else if(keyboard['6'])step=  4;
					else if(keyboard['7'])step=  2;
					else if(keyboard['8'])step=  1;
					;
					//int step=key=='N'?1:(key=='B'?8:64);

					//int step=0x40/(press+1);
					//step+=!step;

					//int idx=press%33;
					if(jxl)
					{
						colortransform_ycocb_fwd(buf2, iw, ih);
						pred_jxl_opt_v2(buf2, iw, ih, jxlparams_i16, 0);
						//for(int idx=0;idx<33;++idx)
						//{
						//	int kc=idx/11;
						//	pred_jxl_optimize(buf2, iw, ih, kc, jxlparams_i16+11*kc, step, idx%11, buf3, 0);
						//}
					}
					else if(pw2)
					{
						colortransform_ycocb_fwd(buf2, iw, ih);
						pred_w2_opt_v2(buf2, iw, ih, pw2_params, 0);
					}
					else
					{
						memset(buf3, 0xFF, (size_t)iw*ih<<2);
						pred_joint_optimize(buf2, iw, ih, jointpredparams, step, buf3, 1);
					}

					//pred_jxl_optimize(buf2, iw, ih, 0, jxlparams_i16   , step, 1, buf3, 1);
					//pred_jxl_optimize(buf2, iw, ih, 1, jxlparams_i16+11, step, 1, buf3, 1);
					//pred_jxl_optimize(buf2, iw, ih, 2, jxlparams_i16+22, step, 1, buf3, 1);

					free(buf2);
					free(buf3);
					//++press;
				}
				update_image();
			}
			//double l0=av_rmse, t0=time_ms(), tnow;
			//int it=0;
			//do
			//{
				//av_rmse=opt_causal_reach2(image, iw, ih, 1, customparam_st, customparam_ct+11, g_lr, 0);
				//++it;
				//tnow=time_ms();
			//}
			//while(fabs(av_rmse-l0)>0.01&&it<10000&&tnow-t0<2000);

			//update_image();
			return 1;
		}
		break;
	//default:
	//	printf("%02X %02X=%c down\n", key, c, c);
	//	if(key=='A')
	//		timer_start();
	//	break;
	}
	return 0;
}
int io_keyup(IOKey key, char c)
{
	switch(key)
	{
	//case KEY_LBUTTON:
	//case KEY_MBUTTON:
	//case KEY_RBUTTON:
	//	printf("Declick at (%d, %d)\n", mx, my);
	//	break;
		
	case KEY_SPACE:
#define		AK(KEY)		case KEY:
	ACTIVE_KEY_LIST
#undef		AK
		count_active_keys(key);
		break;

	//default:
	//	printf("%02X %02X=%c up\n", key, c, c);
	//	if(key=='A')
	//		timer_stop();
	//	break;
	}
	return 0;
}
void io_timer()
{
	float move_speed=keyboard[KEY_SHIFT]?10*cam.move_speed:cam.move_speed;
	if(keyboard['W'])		cam_moveForward(cam, move_speed);
	if(keyboard['A'])		cam_moveLeft(cam, move_speed);
	if(keyboard['S'])		cam_moveBack(cam, move_speed);
	if(keyboard['D'])		cam_moveRight(cam, move_speed);
	if(keyboard['T'])		cam_moveUp(cam, move_speed);
	if(keyboard['G'])		cam_moveDown(cam, move_speed);
	if(keyboard[KEY_UP])	cam_turnUp(cam, key_turn_speed);
	if(keyboard[KEY_DOWN])	cam_turnDown(cam, key_turn_speed);
	if(keyboard[KEY_LEFT])	cam_turnLeft(cam, key_turn_speed);
	if(keyboard[KEY_RIGHT])	cam_turnRight(cam, key_turn_speed);

#if 0
	if(keyboard[KEY_SPACE])
	{
		//if(transforms_mask[ST_FWD_HYBRID3]||transforms_mask[ST_INV_HYBRID3])
		//{
		//	opt_causal_hybrid_r3(image, iw, ih, g_lr);
		//	update_image();
		//}
		//else
		if(transforms_mask[ST_FWD_CUSTOM]||transforms_mask[ST_INV_CUSTOM])
		{
			av_rmse=opt_causal_reach2(image, iw, ih, 1, customparam_st, customparam_ct+11, g_lr, 0);
			update_image();
		}
	}
#endif
	if(transforms_customenabled)
	{
		int update=keyboard[KEY_ENTER]-keyboard[KEY_BKSP];
		if(update)
		{
			double speed;
			if(GET_KEY_STATE(KEY_SHIFT))//fast
				speed=0.1;
			else if(GET_KEY_STATE(KEY_CTRL))//slow
				speed=0.001;
			else//normal speed
				speed=0.01;
			if(customparam_sel<COUNTOF(customparam_ct))
			{
				//undocolortransform(color_transform);
				customparam_ct[customparam_sel]+=speed*update;
				//applycolortransform(color_transform);
			}
			else
			{
				//undospatialtransform(spatialtransform);
				customparam_st[12*customparam_ch_idx+customparam_sel-COUNTOF(customparam_ct)]+=speed*update;
				//applyspatialtransform(spatialtransform);
			}
			update_image();
		}
	}
	else
	{
		if(keyboard[KEY_ENTER])	cam_zoomIn(cam, 1.1f);
		if(keyboard[KEY_BKSP])	cam_zoomOut(cam, 1.1f);
	}
}
void print_i16s(float x, float y, float zoom, const short *row, int count)
{
	int printed=0;
	for(int k=0;k<count;++k)
	{
		short val=row[k];
		printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, " %c%4X", val<0?'-':' ', abs(val));
	}
	print_line(0, x, y, zoom, g_buf, printed, -1, 0, 0);
}
void io_render()
{
	//prof_add("entry");
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
	//prof_add("model");

	if(image)
	{
		switch(mode)
		{
		case VIS_PLANES:			chart_planes_draw();	break;
		case VIS_MESH:				chart_mesh_draw();		break;
		case VIS_MESH_SEPARATE:		chart_mesh_sep_draw();	break;
		case VIS_HISTOGRAM:
			{
				float yoffset=tdy*3;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 1, 0, 1, 0, 1);
				chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0, 0x60, hist, histmax);
			}
			break;
		case VIS_JOINT_HISTOGRAM:	chart_jointhist_draw();	break;
		case VIS_IMAGE:
			{
				int waitstatus=0;
				if(ghMutex)
					waitstatus=WaitForSingleObject(ghMutex, INFINITE);
				float yoffset=tdy*3;
				if(show_full_image)
					display_texture_i(0, w, 0, h, (int*)image, iw, ih, 1, 0, 1, 0, 1);
				else
					display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 1, 0, 1, 0, 1);
				if(waitstatus==WAIT_OBJECT_0)
					ReleaseMutex(ghMutex);
			}
			break;
		case VIS_IMAGE_BLOCK:
			{
				float yoffset=tdy*3, half=blocksize*0.5f;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 1, 0, 1, 0, 1);
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
#if 0
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
#endif
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
		case VIS_IMAGE_TRICOLOR:
			display_texture(0,  iw,    0,  ih,    image_txid[1], 1, 0, 1, 0,     1.f/3);
			display_texture(iw, iw<<1, 0,  ih,    image_txid[1], 1, 0, 1, 2.f/3, 1);
			display_texture(0,  iw,    ih, ih<<1, image_txid[1], 1, 0, 1, 1.f/3, 2.f/3);
			break;
		}
	}
	
	//if(image_txid)
	//	display_texture(0, iw, 0, ih, image_txid, 0.4f);//
	//if(image)
	//	display_texture_i(0, iw, 0, ih, (int*)image, iw, ih, 0.4f);
	//if(im2)
	//	display_texture_i(0, iw*3, h>>1, (h>>1)+ih, im2, iw*2, ih, 0.4f);

	if(transforms_customenabled)
	//if(color_transform==CT_CUSTOM||spatialtransform==ST_CUSTOM)
	{
		char sel[_countof(customparam_ct)+_countof(customparam_st)/6]={0};
		for(int k=0;k<_countof(sel);++k)
			sel[k]=' ';
		sel[customparam_sel]='>';
		float xstep=tdx*guizoom, ystep=tdy*guizoom, x, y;
		//long long prevcolor=set_text_colors(0xFF000000);
		//float x=xstep*2, y=(float)(h>>1);

		//if(!(transforms_mask[ST_FWD_HYBRID3]||transforms_mask[ST_INV_HYBRID3]))
		{
			//color transforms
			//012345678901234567890123456789
			//r-=(>>nnnN.NNN)g+(  nnnN.NNN)b
			x=buttons[0].x1;
			y=buttons[0].y1;
			GUIPrint(0, x, y        , guizoom, "r-=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 0], sel[ 0], customparam_ct[ 0], sel[ 1], sel[ 1], customparam_ct[ 1]);//do not change these strings!
			GUIPrint(0, x, y+ystep  , guizoom, "g-=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 2], sel[ 2], customparam_ct[ 2], sel[ 3], sel[ 3], customparam_ct[ 3]);
			GUIPrint(0, x, y+ystep*2, guizoom, "b-=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[ 4], sel[ 4], customparam_ct[ 4], sel[ 5], sel[ 5], customparam_ct[ 5]);
			GUIPrint(0, x, y+ystep*3, guizoom, "r+=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 6], sel[ 6], customparam_ct[ 6], sel[ 7], sel[ 7], customparam_ct[ 7]);
			GUIPrint(0, x, y+ystep*4, guizoom, "g+=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 8], sel[ 8], customparam_ct[ 8], sel[ 9], sel[ 9], customparam_ct[ 9]);
			GUIPrint(0, x, y+ystep*5, guizoom, "b+=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[10], sel[10], customparam_ct[10], sel[11], sel[11], customparam_ct[11]);

			//spatial transforms
			//0000000000111111111122222222223333333333444444444455555
			//0123456789012345678901234567890123456789012345678901234
			//>>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN
			x=buttons[1].x1;
			y=buttons[1].y1;
			int bk0=set_bk_color(0xA060A060);
			const double *params=customparam_st+12*customparam_ch_idx;
			if(customparam_ch_idx&1)
				GUIPrint(0, x, y-tdy, 1, "Er %d", customparam_ch_idx/2);
			else
				GUIPrint(0, x, y-tdy, 1, "Ch %d", customparam_ch_idx/2);
			GUIPrint(0, x, y        , guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[12], sel[12], params[ 0], sel[13], sel[13], params[ 1], sel[14], sel[14], params[ 2], sel[15], sel[15], params[ 3], sel[16], sel[16], params[ 4]);
			GUIPrint(0, x, y+ystep  , guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[17], sel[17], params[ 5], sel[18], sel[18], params[ 6], sel[19], sel[19], params[ 7], sel[20], sel[20], params[ 8], sel[21], sel[21], params[ 9]);
			set_bk_color(bk0);
			GUIPrint(0, x, y+ystep*2, guizoom, "%c%c%8.3lf %c%c%8.3lf",                                  sel[22], sel[22], params[10], sel[23], sel[23], params[11]);//do not change these strings!
			//set_text_colors(prevcolor);
		}
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
			print_i16s(x, y, 1, params, LOGIC_ROWPARAMS);
			y+=tdy;
			params+=LOGIC_ROWPARAMS;
		}
#ifdef LOGIC_NF1
		print_i16s(x, y, 1, params, LOGIC_NF1);
#else
		print_i16s(x, y, 1, params, LOGIC_NF0);
#endif
		if(info[0]>=0&&ghMutex)
		{
			if(waitstatus==WAIT_OBJECT_0)
				ReleaseMutex(ghMutex);
		}
	}
#endif
	else if(transforms_mask[ST_FWD_JXL]||transforms_mask[ST_INV_JXL])
	{
		float x, y;

		x=buttons[5].x1;
		y=buttons[5].y1;

		print_i16s(x, y, 1, jxlparams_i16   , 11);	y+=tdy;
		print_i16s(x, y, 1, jxlparams_i16+11, 11);	y+=tdy;
		print_i16s(x, y, 1, jxlparams_i16+22, 11);
	}
	else if(transforms_mask[ST_FWD_MM]||transforms_mask[ST_INV_MM])
	{
		float x, y;
		
		x=0;
		y=buttons[5].y1;
		print_i16s(x, y, 1, pw2_params             , PW2_NPARAM);	y+=tdy;
		print_i16s(x, y, 1, pw2_params+PW2_NPARAM  , PW2_NPARAM);	y+=tdy;
		print_i16s(x, y, 1, pw2_params+PW2_NPARAM*2, PW2_NPARAM);
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
	case VIS_PLANES:			mode_str="Planes";			break;
	case VIS_MESH:				mode_str="Combined Mesh";	break;
	case VIS_MESH_SEPARATE:		mode_str="Separate Mesh";	break;
	case VIS_HISTOGRAM:			mode_str="Histogram";		break;
	case VIS_JOINT_HISTOGRAM:	mode_str="Joint Histogram";	break;
	case VIS_IMAGE:				mode_str="Image View";		break;
	case VIS_IMAGE_BLOCK:		mode_str="Image Block";		break;
//	case VIS_IMAGE_E24:			mode_str="Image Exp24";		break;
	case VIS_DWT_BLOCK:			mode_str="DWT Block";		break;
	case VIS_IMAGE_TRICOLOR:	mode_str="Tricolor";		break;
	}
#if 0
	switch(color_transform)
	{
	case CT_NONE:				color_str="NONE";			break;
	case CT_YCoCg:				color_str="YCoCg";			break;
	case CT_YCoCgT:				color_str="YCoCgT";			break;
	case CT_XGZ:				color_str="XGZ";			break;
	case CT_XYZ:				color_str="XYZ";			break;
	case CT_EXP:				color_str="Experimental";	break;
	case CT_LEARNED:			color_str="Learned";		break;
	case CT_CUSTOM:				color_str="Custom";			break;
	}
	switch(spatialtransform)
	{
	case ST_NONE:				space_str="NONE";			break;
	case ST_DIFF2D:				space_str="2D derivative";	break;
	case ST_UNPLANE:			space_str="Unplane";		break;
	case ST_LAZY:				space_str="Lazy DWT";		break;
	case ST_HAAR:				space_str="Haar";	break;
	case ST_SQUEEZE:			space_str="Squeeze";break;
	case ST_CDF53:				space_str="CDF 5/3";break;
	case ST_CDF97:				space_str="CDF 9/7";break;
	case ST_CUSTOM:				space_str="Custom"; break;
	}
#endif
	if(image)
	{
		float x=(float)(w-200), y=tdy*2;
		for(int k=T_NONE+1;k<T_COUNT;++k, y+=tdy)//print available transforms on right
			transforms_printname(x, y, k, -1, transforms_mask[k]?0xA0FF0000FFFFFFFF:0);
		if(transforms)
		{
			x=(float)(w-400);
			y=tdy*2;
			for(int k=0;k<(int)transforms->count;++k, y+=tdy)//print applied transforms on left
				transforms_printname(x, y, transforms->data[k], k, 0);
		}
		float
			cr_combined=3/(1/ch_cr[0]+1/ch_cr[1]+1/ch_cr[2]),
			scale=400,
			xstart=20, xend=(float)w-210, ystart=(float)(h-tdy*5);
		
		float crformat=(float)iw*ih*3/filesize;
		if(xstart<xend)
		{
			float crmax=cr_combined;
			for(int k=0;k<4;++k)//get max CR
			{
				if(isfinite((double)ch_cr[k])&&crmax<ch_cr[k])
					crmax=ch_cr[k];
			}
			if(crmax<crformat)
				crmax=crformat;
			if(xend-scale*crmax<xstart)
				scale=(xend-xstart)/crmax;
			if(scale<5)
				scale=5;

			xstart=xend-scale*ceilf((crmax+1)*2)*0.5f;
			draw_rect(xstart, xend, ystart, (float)h, 0x80808080);//background
			int ks;
			float x;
			if(scale>20)
			{
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
			for(int ks=1;x>xstart;++ks)
			{
				draw_rect(x-1, x+2, ystart, (float)h, 0x70800080);//draw major scale
				x=(float)(xend-scale*ks);
			}
			float barw=4;
			if(mode==VIS_IMAGE_BLOCK
			//	||mode==VIS_IMAGE_E24
				||mode==VIS_DWT_BLOCK)
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
			{
				draw_rect(xend-scale*cr_combined, xend, ystart+tdy*0.5f-barw, ystart+tdy*0.5f+barw+1, 0xFF000000);
				draw_rect(xend-scale*ch_cr[0]   , xend, ystart+tdy*1.5f-barw, ystart+tdy*1.5f+barw+1, 0xFF0000FF);
				draw_rect(xend-scale*ch_cr[1]   , xend, ystart+tdy*2.5f-barw, ystart+tdy*2.5f+barw+1, 0xFF00FF00);
				draw_rect(xend-scale*ch_cr[2]   , xend, ystart+tdy*3.5f-barw, ystart+tdy*3.5f+barw+1, 0xFFFF0000);
			}
			draw_rect(xend-scale*ch_cr[3]   , xend, ystart+tdy*4.5f-barw, ystart+tdy*4.5f+barw+1, 0xFFFF00FF);
			x=xend-crformat*scale;
			draw_line(x, ystart, x-10, ystart-10, 0xFF000000);
			draw_line(x, ystart, x+10, ystart-10, 0xFF000000);
		}
		int prevtxtcolor, prevbkcolor;
		xend=(float)w-200;
		prevbkcolor=set_bk_color(0xC0C0C0C0);
		prevtxtcolor=set_text_color(0xFF000000);GUIPrint(xend, xend, ystart-tdy  , 1, "Format        %9f", crformat);
		set_bk_color(0xE0FFFFFF);
		prevtxtcolor=set_text_color(0xFF000000);GUIPrint(xend, xend, ystart      , 1, "Combined      %9f", cr_combined);
		set_bk_color(0xC0C0C0C0);
		set_text_color(0xFF0000FF);				GUIPrint(xend, xend, ystart+tdy  , 1, "R     %7d %9f", usage[0], ch_cr[0]);
		set_text_color(0xFF00C000);				GUIPrint(xend, xend, ystart+tdy*2, 1, "G     %7d %9f", usage[1], ch_cr[1]);
		set_text_color(0xFFFF0000);				GUIPrint(xend, xend, ystart+tdy*3, 1, "B     %7d %9f", usage[2], ch_cr[2]);
		set_text_color(0xFFFF00FF);				GUIPrint(xend, xend, ystart+tdy*4, 1, "Joint %7d %9f", usage[3], ch_cr[3]);
		set_text_color(prevtxtcolor);
		set_bk_color(prevbkcolor);

		if(transforms_customenabled)
		{
			//double maxloss=0;
			if(losshist&&minloss<maxloss)
			{
				double *data=(double*)losshist->data;
				//for(int k=0;k<(int)losshist->count-1;++k)
				//{
				//	if(maxloss<data[k])
				//		maxloss=data[k];
				//}
				double gain=h/(maxloss-minloss);
				for(int k=0;k<(int)losshist->count-1;++k)
					draw_line((float)k, h-(float)((data[k]-minloss)*gain), (float)(k+1), h-(float)((data[k+1]-minloss)*gain), 0xFFFF00FF);
					//draw_line((float)k*w/losshist->count, (float)((data[k]-minloss)*gain), (float)(k+1)*w/losshist->count, (float)((data[k+1]-minloss)*gain), 0xFFFF00FF);
				double xmark=(customparam_st[12*customparam_ch_idx+profile_idx]-customparam_ct[8])*1000/(customparam_ct[9]-customparam_ct[8]);
				//double xmark=(customparam_st[profile_idx]-customparam_ct[8])*w/(customparam_ct[9]-customparam_ct[8]);
				draw_line((float)xmark, 0, (float)xmark, (float)h, 0xFFFFFF00);
			}
			GUIPrint(0, 0, tdy*3, 1, "RMSE %lf", av_rmse);
			if(minloss<maxloss)
				GUIPrint(0, 200, tdy*3, 1, "[%lf~%lf]", minloss, maxloss);
		}

		float g2=h/combCRhist_max;
		int idx=combCRhist_idx-1, idx2=combCRhist_idx-2;
		idx+=combCRhist_SIZE&-(idx<0);
		idx2+=combCRhist_SIZE&-(idx2<0);
		xstart=(float)(w-(combCRhist_SIZE<<combCRhist_logDX)-200);
		float cx=xstart+(float)(idx<<combCRhist_logDX), cy=h-combCRhist[idx][3]*g2;
		draw_rect_hollow(cx-10, cx+10, cy-10, cy+10, 0xC0C0C0C0);
		if(combCRhist[idx][3]<combCRhist[idx2][3])
			draw_rect(cx-5, cx+5, cy, cy+5, 0xC0C0C0C0);
		else
			draw_rect(cx-5, cx+5, cy-5, cy, 0xC0C0C0C0);
		//draw_ellipse(cx-10, cx+10, cy-5, cy+5, 0xC0C0C0C0);
		for(int k=0;k<combCRhist_SIZE-1;++k)
		{
			if(k!=combCRhist_idx-1)
			{
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][0]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][0]*g2), 0xC00000FF);
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][1]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][1]*g2), 0xC000FF00);
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][2]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][2]*g2), 0xC0FF0000);
				draw_line(xstart+(float)(k<<combCRhist_logDX), h-(float)(combCRhist[k][3]*g2), xstart+(float)((k+1)<<combCRhist_logDX), h-(float)(combCRhist[k+1][3]*g2), 0xC0000000);
			}
		}
		draw_line(xstart, h-g2, xstart+(combCRhist_SIZE<<combCRhist_logDX), h-g2, 0xC0000000);

		//grad2 info
#if 1
		if(transforms_mask[ST_FWD_CUSTOM2]||transforms_mask[ST_INV_CUSTOM2])
		{
			float width=300.f, zoom=1, ystart=tdy*zoom*3;
			x=(float)(w>>2);
			//x=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16s(x, y, zoom, c2_params.c0+k*12, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c0+k*12+5, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c0+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}

			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16s(x, y, zoom, c2_params.c1+k*12, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c1+k*12+5, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c1+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}
			print_i16s(x, y, zoom, c2_params.c1+6*12, 2);		y+=tdy*zoom;
			
			x+=width;
			y=ystart;
			for(int k=0;k<6;++k)
			{
				print_i16s(x, y, zoom, c2_params.c2+k*12, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c2+k*12+5, 5);		y+=tdy*zoom;
				print_i16s(x, y, zoom, c2_params.c2+k*12+10, 2);	y+=tdy*zoom;
				y+=tdy*zoom;
			}
			print_i16s(x, y, zoom, c2_params.c2+6*12, 4);		y+=tdy*zoom;
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
#endif
#if 0
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
#endif
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

		//extern int testhist[3];//
		//GUIPrint(0, 0, tdy*2, 1, "%d %d %d", testhist[0], testhist[1], testhist[2]);//
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
	GUIPrint(0, 0, 0, 1, "p(%f, %f, %f) dx %f a(%f, %f) fov %f", cam.x, cam.y, cam.z, cam.move_speed, cam.ax, cam.ay, atan(cam.tanfov)*180/M_PI*2);
	
	static double t=0;
	double t2=time_ms();
	GUIPrint(0, 0, tdy, 1, "timer %d, fps%11lf, [%2d/%2d] %s", timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str);
	if(mode==VIS_IMAGE)
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
		GUIPrint(0, 0, tdy*2, 1, "%s", mode_str);
	}
	//GUIPrint(0, 0, tdy, 1, "timer %d, fps %10lf, [%d/%d] %s,\tCT: [%d/%d] %s,\tST: [%d/%d] %s", timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str, color_transform+1, CT_COUNT, color_str, spatialtransform+1, ST_COUNT, space_str);
	//if(joint_CR)
	//	GUIPrint(0, 0, tdy*3, 1, "Joint CR %f", joint_CR);
	t=t2;

	//prof_add("finish");
	swapbuffers();
}
int io_quit_request()//return 1 to exit
{
	logic_opt_forceclosethread();
	g_repaint=0;
	//int button_idx=messagebox(MBOX_OKCANCEL, "Are you sure?", "Quit application?");
	//return button_idx==0;

	return 1;
}
void io_cleanup()//cleanup
{
	free(image);
	free(im2);
	array_free(&cpu_vertices);
}