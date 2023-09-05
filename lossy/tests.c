#include"lossy.h"
#define AC_IMPLEMENTATION
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"lodepng.h"
static const char file[]=__FILE__;


void save_ps1(const char *filename, float *fbuf, int iw, int ih, int loud)
{
	int res=iw*ih;
	float vmin=0, vmax=0;
	for(int k=0;k<res;++k)
	{
		if(vmin>fbuf[k])
			vmin=fbuf[k];
		if(vmax<fbuf[k])
			vmax=fbuf[k];
	}
	if(loud)
		printf("[%f ~ %f]\n", vmin, vmax);
	unsigned char *b3=(unsigned char*)malloc(res);
	if(!b3)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	float gain=vmax-vmin;
	if(gain)
	{
		gain=255/gain;
		for(int k=0;k<res;++k)
			b3[k]=(unsigned char)((fbuf[k]-vmin)*gain);
	}
	else
		memset(b3, 0, res);
	lodepng_encode_file(filename, b3, iw, ih, LCT_GREY, 8);
	free(b3);
}

//T44: pseudo-JPEG float (encoder only)
#if 0
#define T44_NBITS 9
static void t44_enc(ABACEncContext *ac, const float *buf, int iw, int ih, unsigned short *hist, float *csizes)
{
	for(int k=0, n=64*(1<<T44_NBITS);k<n;++k)
		hist[k]=0x8000;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			float N, W, NW, curr;
			N=ky-8>=0?buf[idx-iw*8]:0;
			W=kx-8>=0?buf[idx-8]:0;
			NW=kx-8>=0&&ky-8>=0?buf[idx-iw*8-8]:0;
			curr=buf[idx];

			float vmin, vmax;
			if(N<W)
				vmin=N, vmax=W;
			else
				vmin=W, vmax=N;
			float pred=N+W-NW;
			pred=CLAMP(vmin, pred, vmax);
			//int sym=(int)roundf(curr-pred)+(1<<(T44_NBITS-1));
			int sym=(int)roundf(curr)+(1<<(T44_NBITS-1));
			if(sym<0||sym>=(1<<T44_NBITS))
				LOG_ERROR("Range error");
			unsigned short *hk=hist+(1<<T44_NBITS)*((ky&7)<<3|kx&7);
			for(int kb=T44_NBITS-1, hidx=0;kb>=0;--kb)//16 bit
			{
				int bpos=T44_NBITS-1-kb;
				int p0=hk[hidx];
				
				int bit=sym>>kb&1;
				abac_enc(ac, p0, bit);
				

				int prob=bit?0x10000-p0:p0;//
				float bitsize=-log2f((float)prob*(1.f/0x10000));
				csizes[kb]+=bitsize;//


				p0+=((!bit<<16)-p0)>>5;
				p0=CLAMP(1, p0, 0xFFFF);
				hk[hidx]=p0;
				hidx+=1<<bpos;
			}
		}
	}
	abac_enc_flush(ac);
}
int t44_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T44: pseudo-JPEG   Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int res=iw*ih;
	float *buf=(float*)malloc((size_t)res*3*sizeof(float));
	unsigned short *hist=(unsigned short*)malloc(64*(1<<T44_NBITS)*sizeof(short));
	if(!buf||!hist)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	float *luma=buf, *co=buf+res, *cb=buf+res*2;
	cvt_u8_ps(co  , src  , res, 4, 1);//r -> Co
	cvt_u8_ps(luma, src+1, res, 4, 1);//g -> Y
	cvt_u8_ps(cb  , src+2, res, 4, 1);//b -> Cb
	colortransform_ycocb_ps_fwd(co, luma, cb, iw, ih, 1);
	dnsample(co, iw, ih, co);//YUV 4:2:0
	dnsample(cb, iw, ih, cb);
	int w2=iw>>1;
	int h2=ih>>2;
	DCT2_8x8_ps_buf(luma, iw, ih);
	DCT2_8x8_ps_buf(co, w2, h2);
	DCT2_8x8_ps_buf(cb, w2, h2);

	float qmatrix[]=
	{
		16, 11, 10, 16, 24, 40, 51, 61,//https://en.wikipedia.org/w/index.php?title=Quantization_(image_processing)&useskin=monobook#Quantization_matrices
		12, 12, 14, 19, 26, 58, 60, 55,
		14, 13, 16, 24, 40, 57, 69, 56,
		14, 17, 22, 29, 51, 87, 80, 62,
		18, 22, 37, 56, 68, 109, 103, 77,
		24, 35, 55, 64, 81, 104, 113, 92,
		49, 64, 78, 87, 103, 121, 120, 101,
		72, 92, 95, 98, 112, 100, 103, 99,
	};
	quantize(luma, iw, ih, qmatrix);
	quantize(co, w2, h2, qmatrix);
	quantize(cb, w2, h2, qmatrix);

	DList list;
	dlist_init(&list, 1, 1024, 0);

	ABACEncContext ac;
	abac_enc_init(&ac, &list);

	float csizes[T44_NBITS*3]={0};
	int bm[3];
	dlist_push_back(&list, 0, 12);

	t44_enc(&ac, luma, iw, ih, hist, csizes+T44_NBITS);
	bm[0]=(int)list.nobj;

	t44_enc(&ac, co, w2, h2, hist, csizes);
	bm[1]=(int)list.nobj;

	t44_enc(&ac, cb, w2, h2, hist, csizes+T44_NBITS*2);
	bm[2]=(int)list.nobj;

	size_t startidx=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+startidx, bm, 12);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Enc ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
		float chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=T44_NBITS-1;kb>=0;--kb)
		{
			printf("B%2d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc*T44_NBITS+kb;
				float size=csizes[idx];
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);
	}

	dlist_clear(&list);
	free(hist);
	free(buf);
	return 0;
}
#endif


//T45: pseudo-JPEG (fixed precision)

	#define T45_USE_BINARYCODER	//this must be enabled
	#define T45_DISABLE_DC_PRED
//	#define T45_DISABLE_DCT

#define T45_NBITS 16
float t45_q=3;
static void t45_enc(DList *list, const short *buf, int iw, int ih, unsigned short *hist, float *csizes)
{
	ABACEncContext ac;
	abac_enc_init(&ac, list);
	int idx;
#ifndef T45_DISABLE_DC_PRED
	int pred;
#endif
	for(int k=0, n=64*(1<<T45_NBITS);k<n;++k)
		hist[k]=0x8000;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			idx=iw*ky+kx;
			short curr=buf[idx];
#ifndef T45_DISABLE_DC_PRED
			int is_ac=(ky&7)||(kx&7);
			short N, W, NW;
			if(is_ac)//AC
				pred=0;
			else//DC
			{
				N=ky-8>=0?buf[idx-iw*8]:0;
				W=kx-8>=0?buf[idx-8]:0;
				NW=kx-8>=0&&ky-8>=0?buf[idx-iw*8-8]:0;

				short vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
			}
			int sym=(curr-pred+(1<<(T45_NBITS-1)))&((1<<T45_NBITS)-1);
#else
			int sym=(curr+(1<<(T45_NBITS-1)))&((1<<T45_NBITS)-1);
#endif
			//int sym=curr+(1<<(T45_NBITS-1));
			//if(sym<0||sym>=(1<<T45_NBITS))
			//	LOG_ERROR("Range error");
			unsigned short *hk=hist+(1<<T45_NBITS)*((ky&7)<<3|kx&7);
			//unsigned short *hk=hist+((1<<T45_NBITS)&-is_ac);
			for(int kb=T45_NBITS-1, hidx=0;kb>=0;--kb)
			{
				//int bpos=T45_NBITS-1-kb;
				int p0=hk[hidx];
				
				int bit=sym>>kb&1;
				abac_enc(&ac, p0, bit);
				

				int prob=bit?0x10000-p0:p0;//
				float bitsize=-log2f((float)prob*(1.f/0x10000));
				csizes[kb]+=bitsize;//


				p0+=((!bit<<16)-p0)>>4;	//p0 = p0*(1 - 2^-P) + (!bit<<16)*2^-P,		P=15
				p0=CLAMP(1, p0, 0xFFFF);
				hk[hidx]=p0;
				hidx=(hidx<<1|1)+bit;
				//printf("");
			}
		}
	}
	//for(int k=0;k<(1<<T45_NBITS);++k)
	//	printf("%3d 0x%04X\n", k, hist[k]);
	abac_enc_flush(&ac);
}
static void t45_enc2(DList *list, const short *buf, int iw, int ih, int *CDF, float *csizes)
{
	AC32 ac;
	ac32_enc_init(&ac, list);
	int idx;
#ifndef T45_DISABLE_DC_PRED
	int pred;
#endif
	for(int k=0, n=64*(1<<T45_NBITS|1);k<n;++k)
		CDF[k]=k%(1<<T45_NBITS|1);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==384&&ky==0)//
			//	printf("");

			idx=iw*ky+kx;
			int is_ac=(ky&7)||(kx&7);
			short curr=buf[idx];
#ifndef T45_DISABLE_DC_PRED
			short N, W, NW;
			if(is_ac)//AC
				pred=0;
			else//DC
			{
				N=ky-8>=0?buf[idx-iw*8]:0;
				W=kx-8>=0?buf[idx-8]:0;
				NW=kx-8>=0&&ky-8>=0?buf[idx-iw*8-8]:0;

				short vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
			}
			int sym=(curr-pred+(1<<(T45_NBITS-1)))&((1<<T45_NBITS)-1);
#else
			int sym=(curr+(1<<(T45_NBITS-1)))&((1<<T45_NBITS)-1);
#endif
			int *hk=CDF+(1<<T45_NBITS|1)*((ky&7)<<3|kx&7);

			ac32_enc(&ac, hk, 256, sym);
			
			float prob=(float)(hk[sym+1]-hk[sym])/hk[1<<T45_NBITS];//
			float bitsize=-log2f(prob);
			csizes[0]+=bitsize;//

			for(int k=sym+1;k<=(1<<T45_NBITS);++k)
				++hk[k];

			//int sym=curr+(1<<(T45_NBITS-1));
			//if(sym<0||sym>=(1<<T45_NBITS))
			//	LOG_ERROR("Range error");
			//for(int kb=T45_NBITS-1, hidx=0;kb>=0;--kb)//16 bit
			//{
			//	int bpos=T45_NBITS-1-kb;
			//	int p0=hk[hidx];
			//	
			//	int bit=sym>>kb&1;
			//	ac_enc(&ac, p0, bit);
			//	
			//
			//	int prob=bit?0x10000-p0:p0;//
			//	float bitsize=-log2f((float)prob*(1.f/0x10000));
			//	csizes[kb]+=bitsize;//
			//
			//
			//	p0+=((!bit<<16)-p0)>>4;	//p0 = p0*(1 - 2^-P) + (!bit<<16)*2^-P,		P=15
			//	p0=CLAMP(1, p0, 0xFFFF);
			//	hk[hidx]=p0;
			//	hidx=(hidx<<1|1)+bit;
			//}
		}
	}
	ac32_enc_flush(&ac);
}
const unsigned char qmatrices[2][64]=
{
	//default luma/chroma matrices (JPEG spec section K.1)
	{
		16,  11,  10,  16,  24,  40,  51,  61,
		12,  12,  14,  19,  26,  58,  60,  55,
		14,  13,  16,  24,  40,  57,  69,  56,
		14,  17,  22,  29,  51,  87,  80,  62,
		18,  22,  37,  56,  68, 109, 103,  77,
		24,  35,  55,  64,  81, 104, 113,  92,
		49,  64,  78,  87, 103, 121, 120, 101,
		72,  92,  95,  98, 112, 100, 103,  99,
	},
	{
		17,  18,  24,  47,  99,  99,  99,  99,
		18,  21,  26,  66,  99,  99,  99,  99,
		24,  26,  56,  99,  99,  99,  99,  99,
		47,  66,  99,  99,  99,  99,  99,  99,
		99,  99,  99,  99,  99,  99,  99,  99,
		99,  99,  99,  99,  99,  99,  99,  99,
		99,  99,  99,  99,  99,  99,  99,  99,
		99,  99,  99,  99,  99,  99,  99,  99,
	},
};
int t45_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T45: pseudo-JPEG   Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=((res+res2)<<1);
	short *buf=(short*)malloc(buflen*sizeof(short));
#ifdef T45_USE_BINARYCODER
	unsigned short *hist=(unsigned short*)malloc(64*(1<<T45_NBITS)*sizeof(short));
#else
	int *hist=(int*)malloc(64*(1<<T45_NBITS|1)*sizeof(int));
#endif
	if(!buf||!hist)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	short *luma=buf, *co=buf+res, *cb=co+res2, *coeffs=cb+res2;
	ycocb_fwd_subsample_separate(src, iw, ih, luma, co, cb);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
#ifdef T45_USE_BINARYCODER
	float csizes[T45_NBITS*3]={0};
#else
	float csizes[3]={0};
#endif
	int bm[3];
	dlist_push_back(&list, 0, 12);
	
	//short qmatrix[64]=
	//{
	//	16, 11, 10, 16, 24, 40, 51, 61,//https://en.wikipedia.org/w/index.php?title=Quantization_(image_processing)&useskin=monobook#Quantization_matrices
	//	12, 12, 14, 19, 26, 58, 60, 55,
	//	14, 13, 16, 24, 40, 57, 69, 56,
	//	14, 17, 22, 29, 51, 87, 80, 62,
	//	18, 22, 37, 56, 68, 109, 103, 77,
	//	24, 35, 55, 64, 81, 104, 113, 92,
	//	49, 64, 78, 87, 103, 121, 120, 101,
	//	72, 92, 95, 98, 112, 100, 103, 99,
	//};
	//memfill(qmatrix+1, qmatrix, 63*sizeof(short), sizeof(short));
//#ifdef T45_DISABLE_QUANTIZATION
//	float q=0;
//#else
//	float q=4;
//#endif

#ifdef T45_DISABLE_DCT
	for(int k=0;k<res;++k)
		coeffs[k]=luma[k];
#else
	DCT2_8x8_i16_buf(luma, iw, ih, coeffs, t45_q, qmatrices[0]);
#endif
#ifdef T45_USE_BINARYCODER
	t45_enc(&list, coeffs, iw, ih, hist, csizes+T45_NBITS);
#else
	t45_enc2(&list, coeffs, iw, ih, hist, csizes+1);
#endif
	bm[0]=(int)list.nobj;
	
#ifdef T45_DISABLE_DCT
	for(int k=0;k<res2;++k)
		coeffs[k]=co[k];
#else
	DCT2_8x8_i16_buf(co, w2, h2, coeffs, t45_q, qmatrices[1]);
#endif
#ifdef T45_USE_BINARYCODER
	t45_enc(&list, coeffs, w2, h2, hist, csizes);
#else
	t45_enc2(&list, coeffs, w2, h2, hist, csizes);
#endif
	bm[1]=(int)list.nobj;
	
#ifdef T45_DISABLE_DCT
	for(int k=0;k<res2;++k)
		coeffs[k]=cb[k];
#else
	DCT2_8x8_i16_buf(cb, w2, h2, coeffs, t45_q, qmatrices[1]);
#endif
#ifdef T45_USE_BINARYCODER
	t45_enc(&list, coeffs, w2, h2, hist, csizes+T45_NBITS*2);
#else
	t45_enc2(&list, coeffs, w2, h2, hist, csizes+2);
#endif
	bm[2]=(int)list.nobj;

	size_t startidx=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+startidx, bm, 12);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Enc ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
#ifdef T45_USE_BINARYCODER
		float chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=T45_NBITS-1;kb>=0;--kb)
		{
			printf("B%2d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc*T45_NBITS+kb;
				float size=csizes[idx];
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);
#else
		float csize=0;
		for(int kc=0;kc<3;++kc)
		{
			printf("C%d  %15.6lf  %15.6lf\n", kc, csizes[kc]/8, iw*ih*8/csizes[kc]);
			csize+=csizes[kc];
		}
		printf("Total %lf, actual %lld  %lf\n\n", csize/8, list.nobj, (double)iw*ih*3/list.nobj);
#endif
	}

	dlist_clear(&list);
	free(hist);
	free(buf);
	return 0;
}
static void t45_dec(const unsigned char *data, size_t srclen, int iw, int ih, short *buf, unsigned short *hist)
{
	ABACDecContext ac;
	abac_dec_init(&ac, data, data+srclen);
	int idx;
#ifndef T45_DISABLE_DC_PRED
	int pred;
#endif
	for(int k=0, n=64*(1<<T45_NBITS);k<n;++k)
		hist[k]=0x8000;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==384&&ky==0)//
			//	printf("");

			idx=iw*ky+kx;
#ifndef T45_DISABLE_DC_PRED
			int is_ac=(ky&7)||(kx&7);
			short N, W, NW;
			if(is_ac)//AC
				pred=0;
			else//DC
			{
				N=ky-8>=0?buf[idx-iw*8]:0;
				W=kx-8>=0?buf[idx-8]:0;
				NW=kx-8>=0&&ky-8>=0?buf[idx-iw*8-8]:0;

				short vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
			}
#endif
			unsigned short *hk=hist+(1<<T45_NBITS)*((ky&7)<<3|kx&7);
			
			int sym=0;
			for(int kb=T45_NBITS-1, hidx=0;kb>=0;--kb)
			{
				//int bpos=T45_NBITS-1-kb;
				int p0=hk[hidx];

				int bit=abac_dec(&ac, p0);
				sym|=bit<<kb;
				
				p0+=((!bit<<16)-p0)>>4;	//p0 = p0*(1 - 2^-P) + (!bit<<16)*2^-P,		P=15
				p0=CLAMP(1, p0, 0xFFFF);
				hk[hidx]=p0;
				hidx=(hidx<<1|1)+bit;
			}
#ifdef T45_DISABLE_DC_PRED
			int curr=sym-(1<<(T45_NBITS-1));
#else
			int curr=sym+pred-(1<<(T45_NBITS-1));
#endif
			buf[idx]=curr;
		}
	}
}
static void t45_dec2(const unsigned char *data, size_t srclen, int iw, int ih, short *buf, int *CDF)
{
	AC32 ac;
	ac32_dec_init(&ac, data, data+srclen);
	int idx;
#ifndef T45_DISABLE_DC_PRED
	int pred;
#endif
	for(int k=0, n=64*(1<<T45_NBITS|1);k<n;++k)
		CDF[k]=k%(1<<T45_NBITS|1);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(kx==384&&ky==0)//
			//	printf("");

			idx=iw*ky+kx;
			int is_ac=(ky&7)||(kx&7);
#ifndef T45_DISABLE_DC_PRED
			short N, W, NW;
			if(is_ac)//AC
				pred=0;
			else//DC
			{
				N=ky-8>=0?buf[idx-iw*8]:0;
				W=kx-8>=0?buf[idx-8]:0;
				NW=kx-8>=0&&ky-8>=0?buf[idx-iw*8-8]:0;

				short vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
			}
#endif
			int *hk=CDF+(1<<T45_NBITS|1)*((ky&7)<<3|kx&7);

			//if(kx==8&&ky==0)//
			//	printf("");

			int sym=ac32_dec(&ac, hk, 256);
#ifdef T45_DISABLE_DC_PRED
			int curr=sym-(1<<(T45_NBITS-1));
#else
			int curr=sym+pred-(1<<(T45_NBITS-1));
#endif

			buf[idx]=curr;

			for(int k=sym+1;k<=(1<<T45_NBITS);++k)
				++hk[k];
		}
	}
}
int t45_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud)
{
	double t_start=time_ms();
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=((res+res2)<<1);
	short *buf=(short*)malloc(buflen*sizeof(short));
#ifdef T45_USE_BINARYCODER
	unsigned short *hist=(unsigned short*)malloc(64*(1<<T45_NBITS)*sizeof(short));
#else
	int *hist=(int*)malloc(64*(1<<T45_NBITS|1)*sizeof(int));
#endif
	if(!buf||!hist)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}

	int bm[3];
	if(srclen<12)
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	memcpy(bm, data, 12);
	if(srclen<bm[2])
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	short *luma=buf, *co=buf+res, *cb=co+res2, *coeffs=cb+res2;
//#ifdef T45_DISABLE_QUANTIZATION
//	float q=0;
//#else
//	float q=4;
//#endif

#ifdef T45_USE_BINARYCODER
	t45_dec(data+12, bm[0]-12, iw, ih, coeffs, hist);
#else
	t45_dec2(data+12, bm[0]-12, iw, ih, coeffs, hist);
#endif
#ifdef T45_DISABLE_DCT
	for(int k=0;k<res;++k)
		((char*)luma)[k]=(char)coeffs[k];
#else
	DCT3_8x8_i16_buf(coeffs, iw, ih, luma, t45_q, qmatrices[0]);
#endif
	
#ifdef T45_USE_BINARYCODER
	t45_dec(data+bm[0], bm[1]-bm[0], w2, h2, coeffs, hist);
#else
	t45_dec2(data+bm[0], bm[1]-bm[0], w2, h2, coeffs, hist);
#endif
#ifdef T45_DISABLE_DCT
	for(int k=0;k<res2;++k)
		((char*)co)[k]=(char)coeffs[k];
#else
	DCT3_8x8_i16_buf(coeffs, w2, h2, co, t45_q, qmatrices[1]);
#endif
	
#ifdef T45_USE_BINARYCODER
	t45_dec(data+bm[1], bm[2]-bm[1], w2, h2, coeffs, hist);
#else
	t45_dec2(data+bm[1], bm[2]-bm[1], w2, h2, coeffs, hist);
#endif
#ifdef T45_DISABLE_DCT
	for(int k=0;k<res2;++k)
		((char*)cb)[k]=(char)coeffs[k];
#else
	DCT3_8x8_i16_buf(coeffs, w2, h2, cb, t45_q, qmatrices[1]);
#endif
	
	ycocb_inv_upsample_i8(luma, co, cb, iw, ih, dst);
	
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}

	free(hist);
	free(buf);
	return 0;
}

#if 0
#include<conio.h>
void print_matrix_i16(short *buf)
{
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			printf(" %d,", buf[ky<<3|kx]);
		printf("\n");
	}
	printf("\n");
}
void print_matrix_i8(char *buf)
{
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			printf(" %d,", buf[ky<<3|kx]);
		printf("\n");
	}
	printf("\n");
}
void DCT_test()
{
	short buf[64], b2[64];
	char b3[64];

	for(;;)
	{
		srand((unsigned)__rdtsc());
		for(int k=0;k<64;++k)
			buf[k]=(rand()&0xFF)-0x80;

		DCT2_8x8_i16_buf(buf, 8, 8, b2, 1, qmatrices[1]);
		DCT3_8x8_i16_buf(b2, 8, 8, b3, 1, qmatrices[1]);

		print_matrix_i16(buf);
		print_matrix_i16(b2);
		print_matrix_i8(b3);

		_getch();
	}

	pause();
	exit(0);
}
#endif

#if 0
typedef struct opj_mqc_state
{
	unsigned qeval;//the probability of the Least Probable Symbol (0.75->0x8000, 1.5->0xffff)
	unsigned mps;//the Most Probable Symbol (0 or 1)
	unsigned idx_mps;//next state if the next encoded symbol is the MPS
	unsigned idx_lps;//next state if the next encoded symbol is the LPS
} opj_mqc_state_t;
static const opj_mqc_state_t mqc_states[47 * 2] =
{
    {0x5601, 0,  2,  3},// 0
    {0x5601, 1,  3,  2},// 1
    {0x3401, 0,  4, 12},// 2
    {0x3401, 1,  5, 13},// 3
    {0x1801, 0,  6, 18},// 4
    {0x1801, 1,  7, 19},// 5
    {0x0ac1, 0,  8, 24},// 6
    {0x0ac1, 1,  9, 25},// 7
    {0x0521, 0, 10, 58},// 8
    {0x0521, 1, 11, 59},// 9
    {0x0221, 0, 76, 66},//10
    {0x0221, 1, 77, 67},//11
    {0x5601, 0, 14, 13},//12
    {0x5601, 1, 15, 12},//13
    {0x5401, 0, 16, 28},//14
    {0x5401, 1, 17, 29},//15
    {0x4801, 0, 18, 28},//16
    {0x4801, 1, 19, 29},//17
    {0x3801, 0, 20, 28},//18
    {0x3801, 1, 21, 29},//19
    {0x3001, 0, 22, 34},//20
    {0x3001, 1, 23, 35},//21
    {0x2401, 0, 24, 36},//22
    {0x2401, 1, 25, 37},//23
    {0x1c01, 0, 26, 40},//24
    {0x1c01, 1, 27, 41},//25
    {0x1601, 0, 58, 42},//26
    {0x1601, 1, 59, 43},//27
    {0x5601, 0, 30, 29},//28
    {0x5601, 1, 31, 28},//29
    {0x5401, 0, 32, 28},//30
    {0x5401, 1, 33, 29},//31
    {0x5101, 0, 34, 30},//32
    {0x5101, 1, 35, 31},//33
    {0x4801, 0, 36, 32},//34
    {0x4801, 1, 37, 33},//35
    {0x3801, 0, 38, 34},//36
    {0x3801, 1, 39, 35},//37
    {0x3401, 0, 40, 36},//38
    {0x3401, 1, 41, 37},//39
    {0x3001, 0, 42, 38},//40
    {0x3001, 1, 43, 39},//41
    {0x2801, 0, 44, 38},//42
    {0x2801, 1, 45, 39},//43
    {0x2401, 0, 46, 40},//44
    {0x2401, 1, 47, 41},//45
    {0x2201, 0, 48, 42},//46
    {0x2201, 1, 49, 43},//47
    {0x1c01, 0, 50, 44},//48
    {0x1c01, 1, 51, 45},//49
    {0x1801, 0, 52, 46},//50
    {0x1801, 1, 53, 47},//51
    {0x1601, 0, 54, 48},//52
    {0x1601, 1, 55, 49},//53
    {0x1401, 0, 56, 50},//54
    {0x1401, 1, 57, 51},//55
    {0x1201, 0, 58, 52},//56
    {0x1201, 1, 59, 53},//57
    {0x1101, 0, 60, 54},//58
    {0x1101, 1, 61, 55},//59
    {0x0ac1, 0, 62, 56},//60
    {0x0ac1, 1, 63, 57},//61
    {0x09c1, 0, 64, 58},//62
    {0x09c1, 1, 65, 59},//63
    {0x08a1, 0, 66, 60},//64
    {0x08a1, 1, 67, 61},//65
    {0x0521, 0, 68, 62},//66
    {0x0521, 1, 69, 63},//67
    {0x0441, 0, 70, 64},//68
    {0x0441, 1, 71, 65},//69
    {0x02a1, 0, 72, 66},//70
    {0x02a1, 1, 73, 67},//71
    {0x0221, 0, 74, 68},//72
    {0x0221, 1, 75, 69},//73
    {0x0141, 0, 76, 70},//74
    {0x0141, 1, 77, 71},//75
    {0x0111, 0, 78, 72},//76
    {0x0111, 1, 79, 73},//77
    {0x0085, 0, 80, 74},//78
    {0x0085, 1, 81, 75},//79
    {0x0049, 0, 82, 76},//80
    {0x0049, 1, 83, 77},//81
    {0x0025, 0, 84, 78},//82
    {0x0025, 1, 85, 79},//83
    {0x0015, 0, 86, 80},//84
    {0x0015, 1, 87, 81},//85
    {0x0009, 0, 88, 82},//86
    {0x0009, 1, 89, 83},//87
    {0x0005, 0, 90, 84},//88
    {0x0005, 1, 91, 85},//89
    {0x0001, 0, 90, 86},//90
    {0x0001, 1, 91, 87},//91
    {0x5601, 0, 92, 92},//92
    {0x5601, 1, 93, 93},//93
};
typedef struct MQInfoStruct
{
	int idx, nbits;
	unsigned long long seq;//MSB oldest, LSB newest
} MQInfo;
char mq_visited[_countof(mqc_states)];
void test_mq()
{
	SList queue;
	slist_init(&queue, sizeof(MQInfo), 0);
	MQInfo *info=(MQInfo*)QUEUE_ENQUEUE(&queue, 0), *i2;
	info->idx=0;
	info->nbits=0;
	info->seq=0;
	memset(mq_visited, 0, _countof(mqc_states));
	while(queue.count)
	{
		info=(MQInfo*)QUEUE_FRONT(&queue);
		mq_visited[info->idx]=1;

		//print node
		opj_mqc_state_t *s=mqc_states+info->idx;
		int p0=(int)((long long)s->qeval*0x10000/0xAAAA);
		int idx0, idx1;
		if(s->mps)
			idx0=s->idx_lps, idx1=s->idx_mps;
		else
		{
			idx0=s->idx_mps, idx1=s->idx_lps;
			p0=0x10000-p0;
		}
		printf("%3d  %3d  %3d  p0 0x%04X  %6.2lf%%  ", info->idx, idx0, idx1, p0, 100.*p0/0x10000);
		if(info->nbits)
		{
			for(int k=0;k<info->nbits;++k)
				printf("%d", (int)(info->seq>>(info->nbits-1-k)&1));
		}
		else
			printf("-");
		printf("\n");
		mq_visited[info->idx]=1;

		//enqueue unvisited nodes
		if(info->nbits<8)
		//if(!mq_visited[idx0])
		{
			i2=(MQInfo*)QUEUE_ENQUEUE(&queue, 0);
			i2->idx=idx0;
			i2->nbits=info->nbits+1;
			i2->seq=info->seq<<1|0;
		//}
		////if(!mq_visited[idx1])
		//{
			i2=(MQInfo*)QUEUE_ENQUEUE(&queue, 0);
			i2->idx=idx1;
			i2->nbits=info->nbits+1;
			i2->seq=info->seq<<1|1;
		}
		QUEUE_DEQUEUE(&queue);
	}
	printf("Done.\n");
	pause();
	exit(0);
}
#endif


//T46	pseudo-J2K (fixed precision)
#if 0
#define T46_NBITS 16
#define T46_ENABLE_DWT
#define T46_DWT_LEVELS 4
#define T46_QLUMA   0x10001
#define T46_QCHROMA 0x10001
//int t46_hist[256]={0};
short t46_weights[]=
{
	0x01D0, 0x0821, 0x0D67, 0x0821, 0x01D0,
	0x0821, 0x2472, 0x3C16, 0x2470, 0x0821,
	0x0D67, 0x3C16,
};
static void t46_enc(DList *list, const short *buf, int iw, int ih, float *csizes)
{
	//memset(t46_hist, 0, 256*sizeof(int));
	int idx, p0;
	ABACEncContext ac;
	abac_enc_init(&ac, list);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?buf[iw*(ky+(Y))+kx+(X)]:0
			short nb[12]=
			{
				LOAD(-2, -2), LOAD(-1, -2), LOAD(0, -2), LOAD(1, -2), LOAD(2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD(0, -1), LOAD(1, -1), LOAD(2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
			};
#undef  LOAD
			idx=iw*ky+kx;
			int sym=buf[idx];
			int vmin=-0x8000, vmax=0x8000;
			for(int kb=15;kb>=0;--kb)
			{
				int threshold=(vmin+vmax)>>1;//floor half
				p0=0;
				for(int k=0;k<12;++k)
					p0+=(int)(t46_weights[k]/(1+exp(-(nb[k]-threshold)*(11./0x10000))));
					//p0+=(t46_weights[k]>>1)+(int)((long long)(nb[k]-threshold)*(t46_weights[k]>>1)/0x10000);
				//{
				//	p0+=t46_weights[k]&-(nb[k]<threshold);
				//	p0+=(t46_weights[k]>>1)&-(nb[k]==threshold);
				//}
				
				//if(kx==1&&ky==0)
				//	printf("");
				//if(kx==(iw>>1)&&ky==(ih>>1))
				//	printf("");

				p0=CLAMP(1, p0, 0xFFFF);

				int bit=sym>>kb&1;
				abac_enc(&ac, p0, bit);
				
				//++t46_hist[p0>>8];
				int prob=bit?0x10000-p0:p0;//
				float bitsize=-log2f((float)prob*(1.f/0x10000));
				csizes[kb]+=bitsize;//

				if(bit)
					vmin+=1<<kb;
				else
					vmax-=1<<kb;
			}
		}
	}
	abac_enc_flush(&ac);
	//for(int k=0;k<256;++k)
	//	printf("%3d  %d\n", k, t46_hist[k]);
	//printf("\n");
}
int t46_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T46: pseudo-J2K   Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=res+(res2<<1);
	short *buf=(short*)malloc(buflen*sizeof(short));
	short *temp=(short*)malloc(MAXVAR(iw, ih)*sizeof(short));
	if(!buf||!temp)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	short *luma=buf, *co=buf+res, *cb=co+res2;
	ycocb_fwd_subsample_separate(src, iw, ih, luma, co, cb);

	//int vmin, vmax;
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
#ifdef T46_ENABLE_DWT
	ArrayHandle sizes;
	DWTSize *psizes;
	int nsizes;
	sizes=dwt2d_gensizes(iw, ih, 0, 0, T46_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	dwt2d_cdf97_fwd(luma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);
	quantize_dwt_fwd(luma, psizes[nsizes-1].w, psizes[nsizes-1].h, iw, ih, T46_QLUMA);
	quantize_dwt_inv(luma, psizes[nsizes-1].w, psizes[nsizes-1].h, iw, ih, T46_QLUMA);//
	dwt2d_cdf97_inv(luma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);//
	array_free(&sizes);

	sizes=dwt2d_gensizes(w2, h2, 0, 0, T46_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	dwt2d_cdf97_fwd(co, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);
	dwt2d_cdf97_fwd(cb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);
	quantize_dwt_fwd(co, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, T46_QCHROMA);
	quantize_dwt_fwd(cb, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, T46_QCHROMA);
	quantize_dwt_inv(co, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, T46_QCHROMA);//
	quantize_dwt_inv(cb, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, T46_QCHROMA);//
	dwt2d_cdf97_inv(co, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);//
	dwt2d_cdf97_inv(cb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 0x10000);//
	array_free(&sizes);
#endif
	
#if 1
	unsigned char *dst=(unsigned char*)malloc((size_t)res<<2);
	ycocb_inv_upsample_i8(luma, co, cb, iw, ih, dst);
	lodepng_encode_file("dump2.PNG", dst, iw, ih, LCT_RGBA, 8);
	pause();
	exit(0);
#endif
	
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	//for(int k=0;k<buflen;++k)//simple quantization
	//	buf[k]/=T46_Q;//rounded towards zero
	
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	
	float csizes[T46_NBITS*3]={0};
	int bm[3];
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);

	t46_enc(&list, luma, iw, ih, csizes+T46_NBITS);
	bm[0]=(int)list.nobj;

	t46_enc(&list, co, w2, h2, csizes);
	bm[1]=(int)list.nobj;

	t46_enc(&list, cb, w2, h2, csizes+T46_NBITS*2);
	bm[2]=(int)list.nobj;
	
	size_t startidx=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+startidx, bm, 12);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Enc ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
		float chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=T45_NBITS-1;kb>=0;--kb)
		{
			printf("B%2d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc*T45_NBITS+kb;
				float size=csizes[idx];
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);
	}
	dlist_clear(&list);
	free(temp);
	free(buf);
	return 0;
}
static void t46_dec(const unsigned char *data, size_t srclen, int iw, int ih, short *buf)
{
	int idx, p0;
	ABACDecContext ac;
	abac_dec_init(&ac, data, data+srclen);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?buf[iw*(ky+(Y))+kx+(X)]:0
			short nb[12]=
			{
				LOAD(-2, -2), LOAD(-1, -2), LOAD(0, -2), LOAD(1, -2), LOAD(2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD(0, -1), LOAD(1, -1), LOAD(2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
			};
#undef  LOAD
			idx=iw*ky+kx;
			short sym=0;
			int vmin=-0x8000, vmax=0x8000;
			for(int kb=15;kb>=0;--kb)
			{
				int threshold=(vmin+vmax)>>1;//floor half
				p0=0;
				for(int k=0;k<12;++k)
					p0+=(int)(t46_weights[k]/(1+exp(-(nb[k]-threshold)*(11./0x10000))));
					//p0+=(t46_weights[k]>>1)+(int)((long long)(nb[k]-threshold)*(t46_weights[k]>>1)/0x10000);
				//{
				//	p0+=t46_weights[k]&-(nb[k]<threshold);
				//	p0+=(t46_weights[k]>>1)&-(nb[k]==threshold);
				//}

				//if(kx==(iw>>1)&&ky==(ih>>1))
				//if(kx==1&&ky==0)
				//	printf("");

				p0=CLAMP(1, p0, 0xFFFF);

				int bit=abac_dec(&ac, p0);
				sym|=bit<<kb;

				if(bit)
					vmin+=1<<kb;
				else
					vmax-=1<<kb;
			}
			buf[idx]=sym;
		}
	}
}
int t46_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud)
{
	double t_start=time_ms();
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=res+(res2<<1);
	short *buf=(short*)malloc(buflen*sizeof(short));
	short *temp=(short*)malloc(MAXVAR(iw, ih)*sizeof(short));
	if(!buf||!temp)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	int bm[3];
	if(srclen<12)
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	memcpy(bm, data, 12);
	if(srclen<bm[2])
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	short *luma=buf, *co=buf+res, *cb=co+res2;
	t46_dec(data+12, bm[0]-12, iw, ih, luma);
	t46_dec(data+bm[0], bm[1]-bm[0], w2, h2, co);
	t46_dec(data+bm[1], bm[2]-bm[1], w2, h2, cb);

	//int vmin, vmax;
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif

	//for(int k=0;k<buflen;++k)//simple quantization
	//	buf[k]*=T46_Q;//rounded towards zero
	
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
#ifdef T46_ENABLE_DWT
	ArrayHandle sizes;
	sizes=dwt2d_gensizes(iw, ih, 0, 0, T46_DWT_LEVELS);
	dwt2d_cdf97_inv(luma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, T46_QLUMA);
	array_free(&sizes);

	sizes=dwt2d_gensizes(w2, h2, 0, 0, T46_DWT_LEVELS);
	dwt2d_cdf97_inv(co, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, T46_QCHROMA);
	dwt2d_cdf97_inv(cb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, T46_QCHROMA);
	array_free(&sizes);
#endif

#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	ycocb_inv_upsample_i8(luma, co, cb, iw, ih, dst);
	
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}

	free(temp);
	free(buf);
	return 0;
}
#endif


//T47 pseudo-J2K float
#define T47_NBITS 16
#define T47_BAYES_ESTIMATOR
//#define T47_DISABLE_DWT
#if 1
#define DWT_FUNC_FWD cdf97ps2d_fwd
#define DWT_FUNC_INV cdf97ps2d_inv
#else
#define DWT_FUNC_FWD lg53ps2d_fwd
#define DWT_FUNC_INV lg53ps2d_inv
#endif
#define T47_DWT_LEVELS 4
float t47_q=32, t47_z=1;
//#define T47_LUMA_Q   32
//#define T47_LUMA_Z   1
//#define T47_CHROMA_Q 32
//#define T47_CHROMA_Z 1
#ifndef T47_BAYES_ESTIMATOR
short t47_weights[]=
{
	0x01D0, 0x0821, 0x0D67, 0x0821, 0x01D0,
	0x0821, 0x2472, 0x3C16, 0x2470, 0x0821,
	0x0D67, 0x3C16,
};
#endif
static void t47_enc(DList *list, const short *buf, int iw, int ih, unsigned short *hist, float *csizes)
{
	int idx, p0;
	ABACEncContext ac;
	abac_enc_init(&ac, list);
#ifdef T47_BAYES_ESTIMATOR
	for(int k=0, n=1<<T47_NBITS;k<n;++k)
		hist[k]=0x8000;
#endif
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
#ifdef T47_BAYES_ESTIMATOR
			int hidx=0;
#else
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?buf[iw*(ky+(Y))+kx+(X)]:0
			short nb[12]=
			{
				LOAD(-2, -2), LOAD(-1, -2), LOAD(0, -2), LOAD(1, -2), LOAD(2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD(0, -1), LOAD(1, -1), LOAD(2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
			};
#undef  LOAD
			int vmin=-0x8000, vmax=0x8000;
#endif
			idx=iw*ky+kx;
			int sym=buf[idx];
			for(int kb=15;kb>=0;--kb)
			{
#ifdef T47_BAYES_ESTIMATOR
				p0=hist[hidx];
#else
				int threshold=(vmin+vmax)>>1;//floor half
				p0=0;
				for(int k=0;k<12;++k)
					p0+=(int)(t47_weights[k]/(1+exp(-(nb[k]-threshold)*(11./0x10000))));

				//if(kx==1&&ky==0)
				//	printf("");
				//if(kx==(iw>>1)&&ky==(ih>>1))
				//	printf("");

				p0=CLAMP(1, p0, 0xFFFF);
#endif

				int bit=sym>>kb&1;
				abac_enc(&ac, p0, bit);
				
				int prob=bit?0x10000-p0:p0;//
				float bitsize=-log2f((float)prob*(1.f/0x10000));
				csizes[kb]+=bitsize;//
				
#ifdef T47_BAYES_ESTIMATOR
				p0+=((!bit<<16)-p0)>>4;	//p0 = p0*(1 - 2^-P) + (!bit<<16)*2^-P,		P=4
				p0=CLAMP(1, p0, 0xFFFF);
				hist[hidx]=p0;
				hidx=(hidx<<1|1)+bit;
#else
				if(bit)
					vmin+=1<<kb;
				else
					vmax-=1<<kb;
#endif
			}
		}
	}
	abac_enc_flush(&ac);
}
int t47_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud)
{
	double t_start=time_ms();
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T47: pseudo-J2K float   Enc  %s  WH %dx%d\n", g_buf, iw, ih);
	}
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=res+(res2<<1);
	float *fbuf=(float*)malloc(buflen*sizeof(float));
	float *temp=(float*)malloc(MAXVAR(iw, ih)*sizeof(float));
	short *buf=(short*)malloc(buflen*sizeof(short));
	unsigned short *hist=(unsigned short*)malloc((1<<T47_NBITS)*sizeof(short));
	if(!fbuf||!buf||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	float *fluma=fbuf, *fco=fluma+res, *fcb=fco+res2;
	short *luma=buf, *co=luma+res, *cb=co+res2;
	ycocb_fwd_subsample_ps(src, iw, ih, fluma, fco, fcb);

	//float vmin, vmax;
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
#ifndef T47_DISABLE_DWT
	//save_ps1("f_luma.PNG", fluma, iw, ih, 1);//

	ArrayHandle sizes;
	DWTSize *psizes;
	int nsizes;
	sizes=dwt2d_gensizes(iw, ih, 0, 0, T47_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	DWT_FUNC_FWD(fluma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	//save_ps1("f_luma_DWT.PNG", fluma, iw, ih, 1);//
	quantize_dwt_ps_i16(fluma, luma, psizes[nsizes-1].w, psizes[nsizes-1].h, iw, ih, t47_q, t47_z);
//	dequantize_dwt_i16_ps(luma, fluma, psizes[nsizes-1].w, psizes[nsizes-1].h, iw, ih, t47_q, t47_z);//
	//save_ps1("f_luma_DWT_Q.PNG", fluma, iw, ih, 1);//
//	DWT_FUNC_INV(fluma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);//
	//save_ps1("f_luma_DWT_Q_IDWT.PNG", fluma, iw, ih, 1);//
	//save_ps1("f_luma_DWT_IDWT.PNG", fluma, iw, ih, 1);//
	array_free(&sizes);

	sizes=dwt2d_gensizes(w2, h2, 0, 0, T47_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	DWT_FUNC_FWD(fco, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	DWT_FUNC_FWD(fcb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	quantize_dwt_ps_i16(fco, co, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);
	quantize_dwt_ps_i16(fcb, cb, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);
//	dequantize_dwt_i16_ps(co, fco, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);//
//	dequantize_dwt_i16_ps(cb, fcb, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);//
//	DWT_FUNC_INV(fco, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);//
//	DWT_FUNC_INV(fcb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);//
	array_free(&sizes);
#endif
	
#if 0
	unsigned char *dst=(unsigned char*)malloc((size_t)res<<2);
	ycocb_inv_upsample_ps(fluma, fco, fcb, iw, ih, dst);
	lodepng_encode_file("dump2.PNG", dst, iw, ih, LCT_RGBA, 8);
	pause();
	exit(0);
#endif
	
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	//for(int k=0;k<buflen;++k)//simple quantization
	//	buf[k]/=T47_Q;//rounded towards zero
	
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	
	float csizes[T47_NBITS*3]={0};
	int bm[3];
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);

	t47_enc(&list, luma, iw, ih, hist, csizes+T47_NBITS);
	bm[0]=(int)list.nobj;

	t47_enc(&list, co, w2, h2, hist, csizes);
	bm[1]=(int)list.nobj;

	t47_enc(&list, cb, w2, h2, hist, csizes+T47_NBITS*2);
	bm[2]=(int)list.nobj;
	
	size_t startidx=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+startidx, bm, 12);
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Enc ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		
		float chsizes[4]={0};
		printf("\tC0\t\tC1\t\tC2\n\n");
		for(int kb=T45_NBITS-1;kb>=0;--kb)
		{
			printf("B%2d  ", kb);
			for(int kc=0;kc<3;++kc)
			{
				int idx=kc*T45_NBITS+kb;
				float size=csizes[idx];
				printf(" %15.6f", iw*ih/size);
				chsizes[kc]+=size;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[3]=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total%15.6f %15.6f %15.6f %15.6f\n", iw*ih*8/chsizes[0], iw*ih*8/chsizes[1], iw*ih*8/chsizes[2], iw*ih*24/chsizes[3]);
		printf("Total size\t%8d\t\t\t     %15.6f\n", (int)list.nobj, iw*ih*3./list.nobj);
	}
	dlist_clear(&list);
	free(hist);
	free(buf);
	free(temp);
	free(fbuf);
	return 0;
}
static void t47_dec(const unsigned char *data, size_t srclen, int iw, int ih, unsigned short *hist, short *buf)
{
	int idx, p0;
	ABACDecContext ac;
	abac_dec_init(&ac, data, data+srclen);
#ifdef T47_BAYES_ESTIMATOR
	for(int k=0, n=1<<T47_NBITS;k<n;++k)
		hist[k]=0x8000;
#endif
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
#ifdef T47_BAYES_ESTIMATOR
			int hidx=0;
#else
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?buf[iw*(ky+(Y))+kx+(X)]:0
			short nb[12]=
			{
				LOAD(-2, -2), LOAD(-1, -2), LOAD(0, -2), LOAD(1, -2), LOAD(2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD(0, -1), LOAD(1, -1), LOAD(2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
			};
#undef  LOAD
			int vmin=-0x8000, vmax=0x8000;
#endif
			idx=iw*ky+kx;
			short sym=0;
			for(int kb=15;kb>=0;--kb)
			{
#ifdef T47_BAYES_ESTIMATOR
				p0=hist[hidx];
#else
				int threshold=(vmin+vmax)>>1;//floor half
				p0=0;
				for(int k=0;k<12;++k)
					p0+=(int)(t47_weights[k]/(1+exp(-(nb[k]-threshold)*(11./0x10000))));

				//if(kx==(iw>>1)&&ky==(ih>>1))
				//if(kx==1&&ky==0)
				//	printf("");

				p0=CLAMP(1, p0, 0xFFFF);
#endif

				int bit=abac_dec(&ac, p0);
				sym|=bit<<kb;
				
#ifdef T47_BAYES_ESTIMATOR
				p0+=((!bit<<16)-p0)>>4;	//p0 = p0*(1 - 2^-P) + (!bit<<16)*2^-P,		P=4
				p0=CLAMP(1, p0, 0xFFFF);
				hist[hidx]=p0;
				hidx=(hidx<<1|1)+bit;
#else
				if(bit)
					vmin+=1<<kb;
				else
					vmax-=1<<kb;
#endif
			}
			buf[idx]=sym;
		}
	}
}
int t47_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *dst, int loud)
{
	double t_start=time_ms();
	int w2=(iw+1)>>1;
	int h2=(ih+1)>>1;
	int res=iw*ih, res2=w2*h2;
	int buflen=res+(res2<<1);
	short *buf=(short*)malloc(buflen*sizeof(short));
	float *fbuf=(float*)malloc(buflen*sizeof(float));
	float *temp=(float*)malloc(MAXVAR(iw, ih)*sizeof(float));
	unsigned short *hist=(unsigned short*)malloc((1<<T47_NBITS)*sizeof(short));
	if(!buf||!fbuf||!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return 1;
	}
	int bm[3];
	if(srclen<12)
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	memcpy(bm, data, 12);
	if(srclen<bm[2])
	{
		LOG_ERROR("Invalid data");
		return 1;
	}
	float *fluma=fbuf, *fco=fluma+res, *fcb=fco+res2;
	short *luma=buf, *co=luma+res, *cb=co+res2;
	t47_dec(data+12, bm[0]-12, iw, ih, hist, luma);
	t47_dec(data+bm[0], bm[1]-bm[0], w2, h2, hist, co);
	t47_dec(data+bm[1], bm[2]-bm[1], w2, h2, hist, cb);

	//float vmin, vmax;
#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
#ifndef T47_DISABLE_DWT
	ArrayHandle sizes;
	DWTSize *psizes;
	int nsizes;
	sizes=dwt2d_gensizes(iw, ih, 0, 0, T47_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	dequantize_dwt_i16_ps(luma, fluma, psizes[nsizes-1].w, psizes[nsizes-1].h, iw, ih, t47_q, t47_z);//
	DWT_FUNC_INV(fluma, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	array_free(&sizes);

	sizes=dwt2d_gensizes(w2, h2, 0, 0, T47_DWT_LEVELS);
	psizes=(DWTSize*)sizes->data;
	nsizes=(int)sizes->count;
	dequantize_dwt_i16_ps(co, fco, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);//
	dequantize_dwt_i16_ps(cb, fcb, psizes[nsizes-1].w, psizes[nsizes-1].h, w2, h2, t47_q, t47_z);//
	DWT_FUNC_INV(fco, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	DWT_FUNC_INV(fcb, (DWTSize*)sizes->data, 0, (int)sizes->count, 1, temp, 1);
	array_free(&sizes);
#endif

#if 0
	vmin=0, vmax=0;
	for(int k=0;k<buflen;++k)
	{
		//if((unsigned short)buf[k]==0xCDCD)
		//	LOG_ERROR("");
		if(vmin>buf[k])
			vmin=buf[k];
		if(vmax<buf[k])
			vmax=buf[k];
	}
#endif
	ycocb_inv_upsample_ps(fluma, fco, fcb, iw, ih, dst);
	
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
	}

	free(hist);
	free(temp);
	free(buf);
	free(fbuf);
	return 0;
}