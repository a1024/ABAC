#include"pxview3d.h"
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include"stb_image.h"
#include"lodepng.h"
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

	CT_FWD_YCoCg,		CT_INV_YCoCg,
	CT_FWD_YCoCgT,		CT_INV_YCoCgT,
	CT_FWD_XGZ,			CT_INV_XGZ,
	CT_FWD_XYZ,			CT_INV_XYZ,
	CT_FWD_EXP,			CT_INV_EXP,
	CT_FWD_LEARNED,		CT_INV_LEARNED,
	CT_FWD_CUSTOM,		CT_INV_CUSTOM,

	ST_FWD_DIFF2D,		ST_INV_DIFF2D,
	ST_FWD_UNPLANE,		ST_INV_UNPLANE,
	ST_FWD_CUSTOM,		ST_INV_CUSTOM,
	ST_FWD_DCT4,		ST_INV_DCT4,
	ST_FWD_LAZY,		ST_INV_LAZY,
	ST_FWD_HAAR,		ST_INV_HAAR,
	ST_FWD_SQUEEZE,		ST_INV_SQUEEZE,
	ST_FWD_CDF53,		ST_INV_CDF53,
	ST_FWD_CDF97,		ST_INV_CDF97,
	ST_FWD_DEC_DWT,		ST_INV_DEC_DWT,
	ST_FWD_EXPDWT,		ST_INV_EXPDWT,
	ST_FWD_CUSTOM_DWT,	ST_INV_CUSTOM_DWT,

	T_COUNT,
} TransformType;
int transforms_customenabled=0;
char transforms_mask[T_COUNT]={0};
ArrayHandle transforms=0;//array of chars
float guizoom=1.25f;

float blocksize=64;
int blockmx=0, blockmy=0;

float
	pixel_amplitude=10,//4
	mesh_separation=100;//10

ArrayHandle cpu_vertices=0;
unsigned gpu_vertices=0;

ArrayHandle jointhist=0;
float ch_cr[4]={0};
int usage[4]={0};

//int range_intersect(float a1, float a2, float b1, float b2)
//{
//	return a2>b1&&a1<b2;
//}

void transforms_append(unsigned tid)
{
	if(tid<T_COUNT)
	{
		if(!transforms)
			ARRAY_ALLOC(char, transforms, &tid, 1, 0, 0);
		else
			ARRAY_APPEND(transforms, &tid, 1, 1, 0);
		transforms_mask[tid]|=1;
		transforms_customenabled|=
			tid==CT_FWD_CUSTOM||
			tid==CT_INV_CUSTOM||
			tid==ST_FWD_CUSTOM||
			tid==ST_INV_CUSTOM||
			tid==ST_FWD_EXPDWT||
			tid==ST_INV_EXPDWT||
			tid==ST_FWD_CUSTOM_DWT||
			tid==ST_INV_CUSTOM_DWT;
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
				tid2==ST_INV_CUSTOM_DWT;
		}
	}
}
void transforms_removeall()
{
	array_free(&transforms);
	transforms_customenabled=0;
	memset(transforms_mask, 0, T_COUNT);
}
void printtransform(float x, float y, unsigned tid, int place, long long highlight)
{
	const char *a=0;
	switch(tid)
	{
	case T_NONE:				a="NONE";					break;
	case CT_FWD_YCoCg:			a="C  Fwd YCoCg";			break;
	case CT_INV_YCoCg:			a="C  Inv YCoCg";			break;
	case CT_FWD_YCoCgT:			a="C  Fwd YCoCgT";			break;
	case CT_INV_YCoCgT:			a="C  Inv YCoCgT";			break;
	case CT_FWD_XGZ:			a="C  Fwd XGZ";				break;
	case CT_INV_XGZ:			a="C  Inv XGZ";				break;
	case CT_FWD_XYZ:			a="C  Fwd XYZ";				break;
	case CT_INV_XYZ:			a="C  Inv XYZ";				break;
	case CT_FWD_EXP:			a="C  Fwd Experimental";	break;
	case CT_INV_EXP:			a="C  Inv Experimental";	break;
	case CT_FWD_LEARNED:		a="C  Fwd Learned";			break;
	case CT_INV_LEARNED:		a="C  Inv Learned";			break;
	case CT_FWD_CUSTOM:			a="C  Fwd CUSTOM";			break;
	case CT_INV_CUSTOM:			a="C  Inv CUSTOM";			break;
	case ST_FWD_DIFF2D:			a=" S Fwd 2D derivative";	break;
	case ST_INV_DIFF2D:			a=" S Inv 2D derivative";	break;
	case ST_FWD_UNPLANE:		a=" S Fwd Unplane";			break;
	case ST_INV_UNPLANE:		a=" S Inv Unplane";			break;
	case ST_FWD_DCT4:			a=" S Fwd DCT4";			break;
	case ST_INV_DCT4:			a=" S Inv DCT4";			break;
	case ST_FWD_CUSTOM:			a=" S Fwd CUSTOM";			break;
	case ST_INV_CUSTOM:			a=" S Inv CUSTOM";			break;
	case ST_FWD_LAZY:			a=" S Fwd Lazy DWT";		break;
	case ST_INV_LAZY:			a=" S Inv Lazy DWT";		break;
	case ST_FWD_HAAR:			a=" S Fwd Haar";			break;
	case ST_INV_HAAR:			a=" S Inv Haar";			break;
	case ST_FWD_SQUEEZE:		a=" S Fwd Squeeze";			break;
	case ST_INV_SQUEEZE:		a=" S Inv Squeeze";			break;
	case ST_FWD_CDF53:			a=" S Fwd CDF 5/3";			break;
	case ST_INV_CDF53:			a=" S Inv CDF 5/3";			break;
	case ST_FWD_CDF97:			a=" S Fwd CDF 9/7";			break;
	case ST_INV_CDF97:			a=" S Inv CDF 9/7";			break;
	case ST_FWD_DEC_DWT:		a=" S Fwd Dec. DWT";		break;
	case ST_INV_DEC_DWT:		a=" S Inv Dec. DWT";		break;
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
static int hist[768], histmax;
static float blockCR[3]={0};
void chart_hist_update(unsigned char *im, int iw, int ih, int x1, int x2, int y1, int y2)
{
	memset(hist, 0, 768LL*sizeof(int));
	histmax=0;
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
		for(int kc=0;kc<3;++kc)
		{
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
				if(histmax<hist[kc<<8|k])
					histmax=hist[kc<<8|k];
			}
			double entropy=0;
			for(int k=0;k<256;++k)
			{
				int freq=hist[kc<<8|k];
				if(freq)
				{
					double p=(double)freq/count;
					p*=0xFF00;
					++p;
					p/=0x10000;
					entropy-=p*log2(p);
				}
			}
			blockCR[kc]=(float)(8/entropy);
		}
	}
}
void chart_dwthist_update(unsigned char *im, int iw, int ih, int kc, int kband, int x1, int x2, int y1, int y2)
{
	memset(hist+((size_t)kband<<8), 0, 256LL*sizeof(int));
	//histmax=0;
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
		for(int k=0;k<256;++k)
		{
			if(histmax<hist[kband<<8|k])
				histmax=hist[kband<<8|k];
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
	jointhistogram(image, iw*ih, nbits, &jointhist);
	
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
			case CT_FWD_YCoCgT:		colortransform_ycocgt_fwd((char*)image, iw, ih);	break;
			case CT_INV_YCoCgT:		colortransform_ycocgt_inv((char*)image, iw, ih);	break;
			case CT_FWD_XGZ:		colortransform_xgz_fwd((char*)image, iw, ih);		break;
			case CT_INV_XGZ:		colortransform_xgz_inv((char*)image, iw, ih);		break;
			case CT_FWD_XYZ:		colortransform_xyz_fwd((char*)image, iw, ih);		break;
			case CT_INV_XYZ:		colortransform_xyz_inv((char*)image, iw, ih);		break;
			case CT_FWD_EXP:		colortransform_exp_fwd((char*)image, iw, ih);		break;
			case CT_INV_EXP:		colortransform_exp_inv((char*)image, iw, ih);		break;
			case CT_FWD_LEARNED:	colortransform_learned_fwd((char*)image, iw, ih);	break;
			case CT_INV_LEARNED:	colortransform_learned_inv((char*)image, iw, ih);	break;
			case CT_FWD_CUSTOM:		colortransform_custom_fwd((char*)image, iw, ih);	break;
			case CT_INV_CUSTOM:		colortransform_custom_inv((char*)image, iw, ih);	break;

			case ST_FWD_DIFF2D:		image_differentiate((char*)image, iw, ih, 3, 4);	break;
			case ST_INV_DIFF2D:		image_integrate((char*)image, iw, ih, 3, 4);		break;
			case ST_FWD_UNPLANE:	image_unplane((char*)image, iw, ih, 3, 4);			break;
			case ST_INV_UNPLANE:	image_replane((char*)image, iw, ih, 3, 4);			break;
			case ST_FWD_CUSTOM:		image_customst_fwd((char*)image, iw, ih, 3, 4);		break;
			case ST_INV_CUSTOM:		image_customst_inv((char*)image, iw, ih, 3, 4);		break;
			case ST_FWD_DCT4:		image_dct4_fwd((char*)image, iw, ih);				break;
			case ST_INV_DCT4:		image_dct4_inv((char*)image, iw, ih);				break;
			case ST_FWD_LAZY:
			case ST_INV_LAZY:
			case ST_FWD_HAAR:
			case ST_INV_HAAR:
			case ST_FWD_SQUEEZE:
			case ST_INV_SQUEEZE:
			case ST_FWD_CDF53:
			case ST_INV_CDF53:
			case ST_FWD_CDF97:
			case ST_INV_CDF97:
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
						case ST_FWD_CDF53:     dwt2d_cdf53_fwd  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_CDF53:     dwt2d_cdf53_inv  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_CDF97:     dwt2d_cdf97_fwd  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_CDF97:     dwt2d_cdf97_inv  ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_EXPDWT:    dwt2d_exp_fwd    ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_EXPDWT:    dwt2d_exp_inv    ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_FWD_CUSTOM_DWT:dwt2d_custom_fwd ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						case ST_INV_CUSTOM_DWT:dwt2d_custom_inv ((char*)image+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);break;
						}
					}
					array_free(&sizes);
					free(temp);
				}
				break;
			case ST_FWD_DEC_DWT:   dwt2d_dec_fwd((char*)image, iw, ih);	break;
			case ST_INV_DEC_DWT:   dwt2d_dec_inv((char*)image, iw, ih);	break;
			}
		}//for
		addhalf(image, iw, ih, 3, 4);
	}//if transforms

	channel_entropy(image, iw*ih, 3, 4, ch_cr, usage);

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
		chart_hist_update(image, iw, ih, 0, iw, 0, ih);
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
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, 0, pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, mesh_separation, mesh_separation+pixel_amplitude, 0xFF000000);
	draw_AAcuboid_wire(0, (float)iw, 0, (float)ih, mesh_separation*2, mesh_separation*2+pixel_amplitude, 0xFF000000);

	draw_3D_triangles(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[1]);
}
void chart_hist_draw(float x1, float x2, float y1, float y2, int cstart, int cend, int color, unsigned char alpha)
{
	if(histmax)
	{
		float dy=(y2-y1)/3.f, histpx=dy/histmax;
		for(int kc=cstart;kc<cend;++kc)
		{
			int k=1;
			float y=k*histpx*10000;
			for(;y<dy;++k)
			{
				draw_line(x1, y1+(kc+1)*dy-y, x2, y1+(kc+1)*dy-y, color?color:alpha<<24|0xFF<<(kc<<3));//0x40
				y=k*histpx*10000;
			}
		}
		for(int kc=cstart;kc<cend;++kc)
		{
			for(int k=0;k<256;++k)
				draw_rect(x1+k*(x2-x1)/256, x1+(k+1)*(x2-x1)/256, y1+(kc+1)*dy-hist[kc<<8|k]*histpx, y1+(kc+1)*dy, color?color:alpha<<24|0xFF<<(kc<<3));//0x80
		}
	}
}
void chart_jointhist_draw()
{
	draw_AAcuboid_wire(0, 64, 0, 64, 0, 64, 0xFF000000);
	draw_contour3d_rect(&cam, gpu_vertices, (int)(cpu_vertices->count/5), image_txid[2], 0.8f);
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
		colortransform_ycocgt_fwd((unsigned char*)image, iw, ih);
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
		colortransform_ycocgt_inv((unsigned char*)image, iw, ih);
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
AABB buttons[3]={0};//CT left, CT right, ST grid

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
	p->x1=xstep*2, p->x2=p->x1+xstep*30, p->y1=(float)(h>>1), p->y2=p->y1+customparam_ct_h*ystep, ++p;//color params - left
	p->x1=(float)(w>>1), p->x2=p->x1+xstep*54, p->y1=(float)((h>>1)+(h>>2)), p->y2=p->y1+ystep*3, ++p;//spatial params
	p->x1=(float)(w-200), p->x2=(float)w, p->y1=tdy*2, p->y2=p->y1+tdy*T_COUNT, ++p;//transforms list
}
int io_mousemove()//return true to redraw
{
	if(mode==VIS_IMAGE_BLOCK||mode==VIS_DWT_BLOCK)
	{
		if(drag)
		{
			show_mouse(drag);
			drag=0;
		}
		if(GET_KEY_STATE(KEY_LBUTTON))
		{
			blockmx=mx, blockmy=my;
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
int io_mousewheel(int forward)
{
	if(image&&transforms_customenabled)//change custom transform params
	{
		int idx;
		AABB *p=buttons;
		for(idx=0;idx<COUNTOF(buttons)&&!(mx>=p->x1&&mx<p->x2&&my>=p->y1&&my<p->y2);++idx, ++p);
		if(idx<2)
		{
			int x, y;
			switch(idx)
			{
			case 0://color transform
				x=idx<<1|(mx-p->x1>(p->x2-p->x1)*0.5f);
				y=(int)floorf((my-p->y1)*customparam_ct_h/(p->y2-p->y1));
				idx=customparam_ct_w*y+x;
				break;
			case 1://spatial transform
				x=(int)floorf((mx-p->x1)*(customparam_st_reach<<1|1)/(p->x2-p->x1));
				y=(int)floorf((my-p->y1)*(customparam_st_reach+1)/(p->y2-p->y1));
				idx=COUNTOF(customparam_ct)+(customparam_st_reach<<1|1)*y+x;
				break;
			}
			if(idx<COUNTOF(customparam_ct)+COUNTOF(customparam_st))
			{
				int sign=(forward>0)-(forward<0);//abs(forward) is 120
				int ch=(int)floorf((mx-p->x1)/(guizoom*tdx)), ch2;
				//if(GET_KEY_STATE(KEY_SHIFT))//fast
				//	speed=0.1;
				//else if(GET_KEY_STATE(KEY_CTRL))//slow
				//	speed=0.001;
				//else//normal speed
				//	speed=0.01;
				if(idx<COUNTOF(customparam_ct))//color transform
				{
					//000000000011111111112222222222
					//012345678901234567890123456789
					//r-=(>>nnnN.NNN)g+(  nnnN.NNN)b
					ch2=ch-6;
					MODVAR(ch2, ch2, 14);
					if(ch2>=0&&ch2<8)
					{
						ch2-=4;
						ch2+=ch2<0;//skip point
						double speed=pow(10, -ch2);
						customparam_ct[idx]+=sign*speed;
					}
					else
						customparam_ct[idx]+=sign*0.05;
					//ch2=(ch-6)%14;
					//if(ch>=6&&ch<28&&ch2>=3&&ch2<8)
					//{
					//	ch2-=3;
					//	ch2+=ch2>0;//skip point
					//	double speed=pow(10, -ch2);
					//	customparam_ct[idx]+=sign*speed;
					//}
					//else
					//	customparam_ct[idx]+=sign*0.05;

					//undocolortransform(color_transform);
					//switch(ch)
					//{
					//case  9:case 23:customparam_ct[idx]+=sign;break;
					//case 10:case 24:
					//case 11:case 25:customparam_ct[idx]+=sign*0.1;break;
					//case 12:case 26:customparam_ct[idx]+=sign*0.01;break;
					//case 13:case 27:customparam_ct[idx]+=sign*0.001;break;
					//default:		customparam_ct[idx]+=sign*0.01;break;
					//}
					//applycolortransform(color_transform);
				}
				else//spatial transform
				{
					//000000000011111111112222222222333333333344444444445555
					//012345678901234567890123456789012345678901234567890123
					//>>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN
					ch2=ch-2;
					MODVAR(ch2, ch2, 11);
					if(ch>=0&&ch<54&&ch2>=3&&ch2<8)
					{
						ch2-=4;
						ch2+=ch2<0;//skip point
						double speed=pow(10, -ch2);
						customparam_st[idx-COUNTOF(customparam_ct)]+=sign*speed;
					}
					else
						customparam_st[idx-COUNTOF(customparam_ct)]+=sign*0.05;

					//undospatialtransform(spatialtransform);
					//customparam_st[idx-COUNTOF(customparam_ct)]+=sign*speed;
					//applyspatialtransform(spatialtransform);
				}
				update_image();
			}
			else
				goto normal_operation;
		}
		else
			goto normal_operation;
	}
	else
	{
	normal_operation:
		if(mode==VIS_IMAGE_BLOCK||mode==VIS_DWT_BLOCK)
		{
			if(forward>0)
			{
				blocksize*=2;
				if(blocksize>1024)
					blocksize=1024;
			}
			else
			{
				blocksize*=0.5f;
				if(blocksize<1)
					blocksize=1;
			}
		}
		else if(keyboard[KEY_SHIFT])//shift wheel		change cam speed
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
		timer_stop();
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
	if(transforms_customenabled)
	{
		switch(key)
		{
		case KEY_UP:
		case KEY_DOWN:
		case KEY_LEFT:
		case KEY_RIGHT:
			{
				const int idx_limit=COUNTOF(customparam_ct)+COUNTOF(customparam_st);
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
			int idx;
			AABB *p=buttons;
			for(idx=0;idx<COUNTOF(buttons)&&!(mx>=p->x1&&mx<p->x2&&my>=p->y1&&my<p->y2);++idx, ++p);
			switch(idx)
			{
			case 0://color transform params
			case 1://spatial transforms params
				if(transforms_customenabled&&key==KEY_RBUTTON)
				{
					int x, y;
					if(idx)//spatial transform
					{
						x=(int)floorf((mx-p->x1)*(customparam_st_reach<<1|1)/(p->x2-p->x1));
						y=(int)floorf((my-p->y1)*(customparam_st_reach+1)/(p->y2-p->y1));
						idx=COUNTOF(customparam_ct)+((size_t)customparam_st_reach<<1|1)*y+x;
						customparam_st[idx-COUNTOF(customparam_ct)]=0;
					}
					else//color transform
					{
						x=idx<<1|(mx-p->x1>(p->x2-p->x1)*0.5f);
						y=(int)floorf((my-p->y1)*customparam_ct_h/(p->y2-p->y1));
						idx=customparam_ct_w*y+x;
						customparam_ct[idx]=0;
					}
					update_image();
					return 1;
				}
				break;
			case 2://transform list
				{
					int idx=(int)floorf((my-p->y1)*T_COUNT/(p->y2-p->y1))+1;//zero idx is T_NONE
					if(BETWEEN_EXC(1, idx, T_COUNT))
					{
						if(key==KEY_LBUTTON)
						{
							if(GET_KEY_STATE(KEY_CTRL))
								transforms_removeall();
							transforms_append(idx);
						}
						else
							transforms_removebyid(idx);
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
#if 0
		{
			AABB *p=buttons+2;
			if(image&&mx>=p->x1&&mx<p->x2&&my>=p->y1&&my<p->y2)
			{
				int idx=(int)(my-p->y1)*T_COUNT/(int)(p->y2-p->y1)+1;//zero idx is T_NONE
				if(BETWEEN_EXC(1, idx, T_COUNT))
				{
					if(key==KEY_LBUTTON)
						transforms_append(idx);
					else
						transforms_removebyid(idx);
					update_image();
					return 1;
				}
			}
			else if(key==KEY_LBUTTON)
				goto esc;
		}
#endif
		break;
	case KEY_ESC:
toggle_drag:
		if(mode==VIS_IMAGE_BLOCK||mode==VIS_DWT_BLOCK)
		{
			if(key==KEY_LBUTTON)
			{
				blockmx=mx, blockmy=my;
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
		timer_start();
		break;
		
	case KEY_F1:
		messagebox(MBOX_OK, "Controls",
			"WASDTG:\tMove cam\n"
			"Arrow keys:\tTurn cam\n"
			"Mouse1:\tToggle mouse look\n"
			"Ctrl O:\t\tOpen image\n"
			"\n"
			"R:\t\tReset cam\n"
			"Ctrl R:\t\tDisable all transforms\n"
			"Ctrl E:\t\tReset custom transform parameters\n"
			"\n"
			"Mouse1/Mouse2:\tAdd/remove transforms to the list\n"
			"Ctrl Mouse1:\tReplace all transform with this one\n"
			"\n"
			//"Ctrl 1:\t\tApply YCoCg color transform\n"
			//"Ctrl 2:\t\tApply YCoCgT color transform\n"
			//"Ctrl 3:\t\tApply XGZ color transform\n"
			//"Ctrl 4:\t\tApply XYZ color transform\n"
			//"Ctrl 5:\t\tApply experimental transform\n"
			//"Ctrl 6:\t\tApply learned transform\n"
			//"Ctrl 7:\t\tApply custom color transform\n"
			//"Ctrl 0:\t\tUndo color transform\n"
			//"\n"
			//"1:\t\tApply 2D derivative\n"
			//"2:\t\tApply Unplane\n"
			//"3:\t\tApply Lazy DWT\n"
			//"4:\t\tApply Haar transform\n"
			//"5:\t\tApply Squeeze transform from JPEG XL\n"
			//"6:\t\tApply CDF 5/3 transform\n"
			//"7:\t\tApply CDF 9/7 transform\n"
			//"8:\t\tApply custom spatial transform\n"
			//"0:\t\tUndo spatial transform\n"
			//"\n"
			"M / Shift M:\tSelect among:\n"
			"\t1: Levels\n"
			"\t2: Mesh\n"
			"\t3: Separate Mesh\n"
			"\t4: Histogram\n"
			"\t5: Joint Histogram\n"
			"\t6: Image\n"
			"\t7: Separate components\n"
			//"[Levels / Mesh / SeparateMesh / Histogram / Image]\n"
		);
		//prof_on=!prof_on;
		return 0;
	case 'R':
		if(GET_KEY_STATE(KEY_CTRL)&&image)
		{
			transforms_removeall();
			update_image();
#if 0
			if(fn)
			{
				if(im0)
					free(im0);
				im0=stbi_load((char*)fn->data, &iw, &ih, &nch0, 4);
				if(im0)
				{
					array_free(&transforms);
					memset(transforms_mask, 0, T_COUNT);
					//color_transform=CT_NONE;
					//spatialtransform=CT_NONE;
					update_image();
				}
			}
#endif
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
		//if(keyboard[KEY_CTRL])
		//{
		//	if(openfiles)
		//	{
		//		int temp1[2]={};
		//		FREE_ARRAY_OF_POINTERS(openfiles, temp1);
		//	}
		//	if(keyboard[KEY_SHIFT])
		//		openfiles=dialog_open_folder(1);
		//	else
		//		openfiles=dialog_open_file(file_filters, SIZEOF(file_filters), 1);
		//}
		return 1;
#if 0
	case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':case '8':case '9':
		if(image)
		{
			int tidx=key-'0';
			if(GET_KEY_STATE(KEY_CTRL))//color transform
			{
				if(tidx<CT_COUNT&&tidx!=color_transform)
				{
					transforms_append(0, tidx);
					//undocolortransform(color_transform);
					//color_transform=tidx;
					//applycolortransform(color_transform);
					update_image();
					return 1;
				}
			}
			else//spatial transform
			{
				if(tidx<ST_COUNT&&tidx!=spatialtransform)
				{
					transforms_append(0, tidx);
					//undospatialtransform(spatialtransform);
					//spatialtransform=tidx;
					//applyspatialtransform(spatialtransform);
					update_image();
					return 1;
				}
			}
		}
		break;
#endif
	//case 'E':
	//	if(image)
	//	{
	//		differentiate_image(image, iw, ih, 3, 4);
	//		update_image();
	//	}
	//	return 1;
	case 'C':
		if(image&&GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)//copy custom transform value
		{
			int printed=0;
			for(int ky=0;ky<customparam_ct_h;++ky)
			{
				for(int kx=0;kx<customparam_ct_w;++kx)
					printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "\t%g", customparam_ct[customparam_ct_w*ky+kx]);
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "\n");
			}
			int stw=customparam_st_reach<<1|1;
			for(int k=0;k<COUNTOF(customparam_st);++k)
			{
				int x=k%stw, y=k/stw;
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "%g%c", customparam_st[k], x<stw-1&&k<COUNTOF(customparam_st)-1?'\t':'\n');
			}

			//for(int k=0;k<COUNTOF(customparam_ct);++k)
			//	printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "%15.13lf%c", customparam_ct[k], k<COUNTOF(customparam_ct)-1?'\t':'\n');
			//for(int k=0;k<COUNTOF(customparam_st);++k)
			//	printed+=snprintf(g_buf+printed, G_BUF_SIZE-printed, "%15.13lf%c", customparam_st[k], k<COUNTOF(customparam_st)-1?'\t':'\n');

			//int printed=snprintf(g_buf, G_BUF_SIZE, "%15.13lf", customparam_sel<COUNTOF(customparam_ct)?customparam_ct[customparam_sel]:customparam_st[customparam_sel-COUNTOF(customparam_ct)]);
			copy_to_clipboard(g_buf, printed);

			//color_transform(image, iw, ih);
			//update_image();
		}
		return 1;
	case 'V':
		if(image&&GET_KEY_STATE(KEY_CTRL)&&transforms_customenabled)//paste custom transform value
		{
			ArrayHandle text=paste_from_clipboard(0);
			if(text)
			{
				int k=0, idx=0;
				for(;k<COUNTOF(customparam_ct);++k)
				{
					for(;idx<text->count&&isspace(text->data[idx]);++idx);
					customparam_ct[k]=atof(text->data+idx);
					for(;idx<text->count&&!isspace(text->data[idx]);++idx);
				}
				for(;k<COUNTOF(customparam_ct)+COUNTOF(customparam_st);++k)
				{
					for(;idx<text->count&&isspace(text->data[idx]);++idx);
					customparam_st[k-COUNTOF(customparam_ct)]=atof(text->data+idx);
					for(;idx<text->count&&!isspace(text->data[idx]);++idx);
				}

				//double val=atof(text->data);
				//if(customparam_sel<COUNTOF(customparam_ct))
				//	customparam_ct[customparam_sel]=val;
				//else
				//	customparam_st[customparam_sel-COUNTOF(customparam_ct)]=val;

				array_free(&text);
				update_image();
				return 1;
			}
		}
		break;
	case 'M':
		if(image)
		{
			MODVAR(mode, mode+1-(GET_KEY_STATE(KEY_SHIFT)<<1), VIS_COUNT);
			update_image();
		}
		return 1;
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
	if(transforms_customenabled)
	//if(color_transform==CT_CUSTOM||spatialtransform==ST_CUSTOM)
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
				customparam_st[customparam_sel-COUNTOF(customparam_ct)]+=speed*update;
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
				chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0, 0x60);
			}
			break;
		case VIS_JOINT_HISTOGRAM:	chart_jointhist_draw();	break;
		case VIS_IMAGE:
			{
				float yoffset=tdy*3;
				display_texture_i(0, iw, (int)yoffset, (int)yoffset+ih, (int*)image, iw, ih, 1, 0, 1, 0, 1);
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
						x1=0, x2=blocksize;
					if(x2>iw)
						x1=iw-blocksize, x2=(float)iw;
				
					if(ih<blocksize)
						y1=0, y2=(float)ih;
					if(y1<0)
						y1=0, y2=blocksize;
					if(y2>ih)
						y1=ih-blocksize, y2=(float)ih;

					draw_rect_hollow(x1, x2, yoffset+y1, yoffset+y2, 0xFFFF00FF);
					chart_hist_update(image, iw, ih, (int)x1, (int)x2, (int)y1, (int)y2);
					chart_hist_draw(0, (float)w, 0, (float)h, 0, 3, 0, 0x30);
					GUIPrint(0, 0, tdy*2, 1, "RGB[%8f %8f %8f] %8f %g", blockCR[0], blockCR[1], blockCR[2], 3/(1/blockCR[0]+1/blockCR[1]+1/blockCR[2]), blocksize);
				}
			}
			//{
			//	//int x=(int)ceilf(tdx*6.5f);
			//	int x=(int)floorf(tdx*6.5f);
			//	display_texture_i(x, x+iw, (int)(tdy*3), (int)(tdy*3+ih), (int*)image, iw, ih, 1);
			//}
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
							mx1=(float)x1, mx2=x1+blocksize;
						if(mx2>x2)
							mx1=x2-blocksize, mx2=(float)x2;
				
						if(h2<blocksize)
							my1=(float)y1, my2=(float)y2;
						if(my1<y1)
							my1=(float)y1, my2=y1+blocksize;
						if(my2>y2)
							my1=y2-blocksize, my2=(float)y2;
						
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
						histmax=0;
						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 2, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 2, 3, histcolor, 0);

						mx1=x1+(mx1-x1)*0.5f;
						mx2=x1+(mx2-x1)*0.5f;
						my1=y1+(my1-y1)*0.5f;
						my2=y1+(my2-y1)*0.5f;

						histmax=0;
						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 1, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 1, 2, histcolor, 0);
						
						mx1=x1+(mx1-x1)*0.5f;
						mx2=x1+(mx2-x1)*0.5f;
						my1=y1+(my1-y1)*0.5f;
						my2=y1+(my2-y1)*0.5f;

						histmax=0;
						draw_rect_hollow(mx1, mx2, my1, my2, 0xFFFF00FF);
						chart_dwthist_update(image, iw, ih, kc, 0, (int)mx1-x1, (int)mx2-x1, (int)my1-y1, (int)my2-y1);
						chart_hist_draw(0, (float)w, 0, (float)h, 0, 1, histcolor, 0);

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
			display_texture(iw, iw<<1, 0,  ih,    image_txid[1], 1, 0, 1, 1.f/3, 2.f/3);
			display_texture(0,  iw,    ih, ih<<1, image_txid[1], 1, 0, 1, 2.f/3, 1);
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
		char sel[COUNTOF(customparam_ct)+COUNTOF(customparam_st)]={0};
		for(int k=0;k<COUNTOF(sel);++k)
			sel[k]=' ';
		sel[customparam_sel]='>';
		float xstep=tdx*guizoom, ystep=tdy*guizoom;
		//long long prevcolor=set_text_colors(0xFF000000);
		float x=xstep*2, y=(float)(h>>1);
		//012345678901234567890123456789
		//r-=(>>nnnN.NNN)g+(  nnnN.NNN)b
		GUIPrint(0, x, y        , guizoom, "r-=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 0], sel[ 0], customparam_ct[ 0], sel[ 1], sel[ 1], customparam_ct[ 1]);//do not change these strings!
		GUIPrint(0, x, y+ystep  , guizoom, "g-=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 2], sel[ 2], customparam_ct[ 2], sel[ 3], sel[ 3], customparam_ct[ 3]);
		GUIPrint(0, x, y+ystep*2, guizoom, "b-=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[ 4], sel[ 4], customparam_ct[ 4], sel[ 5], sel[ 5], customparam_ct[ 5]);
		GUIPrint(0, x, y+ystep*3, guizoom, "r+=(%c%c%8.3lf)g+(%c%c%8.3lf)b", sel[ 6], sel[ 6], customparam_ct[ 6], sel[ 7], sel[ 7], customparam_ct[ 7]);
		GUIPrint(0, x, y+ystep*4, guizoom, "g+=(%c%c%8.3lf)r+(%c%c%8.3lf)b", sel[ 8], sel[ 8], customparam_ct[ 8], sel[ 9], sel[ 9], customparam_ct[ 9]);
		GUIPrint(0, x, y+ystep*5, guizoom, "b+=(%c%c%8.3lf)r+(%c%c%8.3lf)g", sel[10], sel[10], customparam_ct[10], sel[11], sel[11], customparam_ct[11]);

		x=(float)(w>>1);
		y=(float)((h>>1)+(h>>2));
		//000000000011111111112222222222333333333344444444445555
		//012345678901234567890123456789012345678901234567890123
		//>>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN >>nnnN.NNN
		int bk0=set_bk_color(0xA060A060);
		GUIPrint(0, (float)x, (float)(y        ), guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[12], sel[12], customparam_st[ 0], sel[13], sel[13], customparam_st[ 1], sel[14], sel[14], customparam_st[ 2], sel[15], sel[15], customparam_st[ 3], sel[16], sel[16], customparam_st[ 4]);
		GUIPrint(0, (float)x, (float)(y+ystep  ), guizoom, "%c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf %c%c%8.3lf", sel[17], sel[17], customparam_st[ 5], sel[18], sel[18], customparam_st[ 6], sel[19], sel[19], customparam_st[ 7], sel[20], sel[20], customparam_st[ 8], sel[21], sel[21], customparam_st[ 9]);
		set_bk_color(bk0);
		GUIPrint(0, (float)x, (float)(y+ystep*2), guizoom, "%c%c%8.3lf %c%c%8.3lf",                                  sel[22], sel[22], customparam_st[10], sel[23], sel[23], customparam_st[11]);//do not change these strings!
		//set_text_colors(prevcolor);
	}
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
			printtransform(x, y, k, -1, transforms_mask[k]?0xA0FF0000FFFFFFFF:0);
		if(transforms)
		{
			x=(float)(w-400);
			y=tdy*2;
			for(int k=0;k<(int)transforms->count;++k, y+=tdy)//print applied transforms on left
				printtransform(x, y, transforms->data[k], k, 0);
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
			if(mode==VIS_IMAGE_BLOCK||mode==VIS_DWT_BLOCK)
			{
				float blockcombined=3/(1/blockCR[0]+1/blockCR[1]+1/blockCR[2]);
				draw_rect(xend-scale*blockcombined, xend, ystart+tdy*0.5f-barw, ystart+tdy*0.5f  , 0xFF404040);
				draw_rect(xend-scale*blockCR[0]   , xend, ystart+tdy*1.5f-barw, ystart+tdy*1.5f  , 0xFF0000B0);
				draw_rect(xend-scale*blockCR[1]   , xend, ystart+tdy*2.5f-barw, ystart+tdy*2.5f  , 0xFF00B000);
				draw_rect(xend-scale*blockCR[2]   , xend, ystart+tdy*3.5f-barw, ystart+tdy*3.5f  , 0xFFB00000);
				draw_rect(xend-scale*cr_combined  , xend, ystart+tdy*0.5f, ystart+tdy*0.5f+barw+1, 0xFF000000);
				draw_rect(xend-scale*ch_cr[0]     , xend, ystart+tdy*1.5f, ystart+tdy*1.5f+barw+1, 0xFF0000FF);
				draw_rect(xend-scale*ch_cr[1]     , xend, ystart+tdy*2.5f, ystart+tdy*2.5f+barw+1, 0xFF00FF00);
				draw_rect(xend-scale*ch_cr[2]     , xend, ystart+tdy*3.5f, ystart+tdy*3.5f+barw+1, 0xFFFF0000);
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
		set_text_color(0xFF00FF00);				GUIPrint(xend, xend, ystart+tdy*2, 1, "G     %7d %9f", usage[1], ch_cr[1]);
		set_text_color(0xFFFF0000);				GUIPrint(xend, xend, ystart+tdy*3, 1, "B     %7d %9f", usage[2], ch_cr[2]);
		set_text_color(0xFFFF00FF);				GUIPrint(xend, xend, ystart+tdy*4, 1, "Joint %7d %9f", usage[3], ch_cr[3]);
		set_text_color(prevtxtcolor);
		set_bk_color(prevbkcolor);
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
	GUIPrint(0, 0, 0, 1, "p(%f, %f, %f) dx %f a(%f, %f) fov %f", cam.x, cam.y, cam.z, cam.move_speed, cam.ax, cam.ay, atan(cam.tanfov)*180/M_PI*2);
	
	static double t=0;
	double t2=time_ms();
	GUIPrint(0, 0, tdy, 1, "timer %d, fps %10lf, [%d/%d] %s", timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str);
	//GUIPrint(0, 0, tdy, 1, "timer %d, fps %10lf, [%d/%d] %s,\tCT: [%d/%d] %s,\tST: [%d/%d] %s", timer, 1000./(t2-t), mode+1, VIS_COUNT, mode_str, color_transform+1, CT_COUNT, color_str, spatialtransform+1, ST_COUNT, space_str);
	//if(joint_CR)
	//	GUIPrint(0, 0, tdy*3, 1, "Joint CR %f", joint_CR);
	t=t2;

	//prof_add("finish");
}
int io_quit_request()//return 1 to exit
{
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