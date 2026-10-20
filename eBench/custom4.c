#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#if 0
#	define _USE_MATH_DEFINES
#endif
#include<math.h>
#include<immintrin.h>

short lossyconv_clipboard=0;
int lossyconv_page=0;/*[0~15]: 4 channels max * 4 layers	layer<<2|ch*/
short lossyconv_params[5*5*4*4]={0};/*(r2 = 5*5) * 4 channels max * 4 layers			9-bit precision: sign_bit.2.5.clamp_bit*/
uint8_t lossyconv_stride[2*4]={0}, lossyconv_offset[2*4]={0};/*2 dimensions * 4 layers*/
uint8_t lossyconv_causalRCT[4]={0};
#if 0
#define LOAD(BUF, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
#endif
void pred_lossyconv(Image *src)
{
	Image *dst;
	int *p1, *p2, kb, k, kc, ky, kx, ky2, kx2, idx2;

	dst=0;
	image_copy(&dst, src);
#if 0
	if(!dst)
	{
		CRASH("Alloc error");
		return;
	}
#endif
	p1=(int32_t*)HDR_DATA(src);
	p2=(int32_t*)HDR_DATA(dst);
	for(kb=0;kb<4;++kb)/*no need for memcpy due to even number of stages:  src->dst->src*/
	{
		int
			xstart=lossyconv_offset[kb<<1|0], xstride=lossyconv_stride[kb<<1|0]+1,
			ystart=lossyconv_offset[kb<<1|1], ystride=lossyconv_stride[kb<<1|1]+1;
		int skip_layer=1;

		for(k=0;k<5*5*4;++k)
		{
			if(lossyconv_params[5*5*4*kb+k])
				skip_layer=0;
		}
		if(skip_layer)
			memcpy(p2, p1, sizeof(int[4])*src->iw*src->ih);
		else
		{
			for(kc=0;kc<4;++kc)
			{
				int depth, half;
				short *curr_filt;

				depth=src->depth[kc];
				if(!depth)
					continue;
				half=1<<depth>>1;
				curr_filt=lossyconv_params+5*5*(4*kb+kc);
				for(ky=ystart;ky<src->ih-(ystride-1);ky+=ystride)
				{
					for(kx=xstart;kx<src->iw-(xstride-1);kx+=xstride)
					{
						int idx=src->iw*ky+kx;
						int offset=0, pred=0, vmin=-half, vmax=half-1, uninit=1;
						for(ky2=-2, idx2=0;ky2<=2;++ky2)
						{
							for(kx2=-2;kx2<=2;++kx2, ++idx2)
							{
								if((uint32_t)(ky+ky2)<(uint32_t)src->ih&&(uint32_t)(kx+kx2)<(uint32_t)src->iw)
								{
									short param=curr_filt[idx2];
									int idx3=src->iw*(ky+ky2)+kx+kx2;
									int nb=p1[idx3<<2|kc];
									if(kc!=1&&lossyconv_causalRCT[kb])
										nb-=p1[idx3<<2|1];
									
									pred+=(param>>1)*nb;
									if(param&1)
									{
										if(uninit)
										{
											uninit=0;
											vmin=nb, vmax=nb;
										}
										else
										{
											UPDATE_MIN(vmin, nb);
											UPDATE_MAX(vmax, nb);
										}
									}
								}
							}
						}
						if(kc!=1&&lossyconv_causalRCT[kb])
							offset+=p1[idx<<2|1];
						pred+=16;
						pred>>=5;
						pred=CLAMP(vmin, pred, vmax);
						pred+=offset;
						pred=CLAMP(-half, pred, half);
#if 0
						int nb[]=
						{
							LOAD(p1, -2, -2),
							LOAD(p1, -1, -2),
							LOAD(p1,  0, -2),
							LOAD(p1,  1, -2),
							LOAD(p1,  2, -2),
							LOAD(p1, -2, -1),
							LOAD(p1, -1, -1),
							LOAD(p1,  0, -1),
							LOAD(p1,  1, -1),
							LOAD(p1,  2, -1),
							LOAD(p1, -2,  0),
							LOAD(p1, -1,  0),
							LOAD(p1,  0,  0),
							LOAD(p1,  1,  0),
							LOAD(p1,  2,  0),
							LOAD(p1, -2,  1),
							LOAD(p1, -1,  1),
							LOAD(p1,  0,  1),
							LOAD(p1,  1,  1),
							LOAD(p1,  2,  1),
							LOAD(p1, -2,  2),
							LOAD(p1, -1,  2),
							LOAD(p1,  0,  2),
							LOAD(p1,  1,  2),
							LOAD(p1,  2,  2),
						};
#if 0
						if(ky==10&&kx==10)
							printf("");
#endif
						int pred=0;
						int vmin=-half, vmax=half-1, uninit=1;
						for(k=0;k<5*5;++k)
						{
							pred+=(curr_filt[k]>>1)*nb[k];
							if(curr_filt[k]&1)
							{
								if(uninit)
								{
									uninit=0;
									vmin=nb[k], vmax=nb[k];
								}
								else
								{
									UPDATE_MIN(vmin, nb[k]);
									UPDATE_MAX(vmax, nb[k]);
								}
							}
						}
						pred+=8;
						pred>>=4;
						pred=CLAMP(vmin, pred, vmax);
#endif
						{
							int curr=p1[idx<<2|kc];
							if(keyboard[KEY_ALT])
								curr=pred;
							else
							{
								curr-=pred;
								curr<<=32-depth;/*signed MA*/
								curr>>=32-depth;
							}
							p2[idx<<2|kc]=curr;
						}
					}
				}
			}
		}
		{
			int *temp;
			SWAPVAR(p1, p2, temp);
		}
	}
	free(dst);
}
void pred_WC(Image *src)
{
	int wpmask;
	Image *dst;
	int *buf;
	ptrdiff_t bufsize;
	int k1, k2, kc, ky, kx, idx, ky2, idx2, kx2, kp;

	wpmask=0;
	for(k1=0;k1<16;++k1)
	{
		short *curr_filt=lossyconv_params+5*5*k1;
		for(k2=0;k2<25;++k2)
		{
			if(curr_filt[k2])
			{
				wpmask|=1<<k1;
				break;
			}
		}
	}
	if(!wpmask)
		return;

	bufsize=(src->iw+(intptr_t)16)*sizeof(int[4*16]);/*4 padded rows * 16 pred errors*/
	buf=(int*)malloc(bufsize);
	if(!buf)
	{
		CRASH("Alloc error");
		return;
	}
	dst=0;
	image_copy(&dst, src);
	for(kc=0;kc<4;++kc)
	{
		int depth, half;
		int preds[16]={0}, nb[25]={0};
#if 0
		int weights[16]={0};
#endif

		depth=src->depth[kc];
		if(!depth)
			continue;
		half=1<<depth>>1;
		memset(buf, 0, bufsize);
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			int32_t *rows[4]={0};

			rows[0]=buf+(((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)<<4);
			rows[1]=buf+(((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)<<4);
			rows[2]=buf+(((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)<<4);
			rows[3]=buf+(((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)<<4);
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int64_t lpred, wsum;
				int weight;
				int offset=0, pred=0, vmin=-half, vmax=half-1, uninit=1;
				int
					*eNW	=rows[1]-1*16,
					*eN	=rows[1]+0*16,
					*eNE	=rows[1]+1*16,
					*eNEE	=rows[1]+2*16,
					*eWW	=rows[0]-2*16,
					*eW	=rows[0]-1*16,
					*ecurr	=rows[0]+0*16;

				for(ky2=-2, idx2=0;ky2<=2;++ky2)
				{
					for(kx2=-2;kx2<=2;++kx2, ++idx2)
					{
						if((uint32_t)(ky+ky2)<(uint32_t)src->ih&&(uint32_t)(kx+kx2)<(uint32_t)src->iw)
						{
							int idx3=src->iw*(ky+ky2)+kx+kx2;
							nb[idx2]=((int32_t*)HDR_DATA(dst))[idx3<<2|kc];
							if(kc!=1&&lossyconv_causalRCT[0])
								nb[idx2]-=((int32_t*)HDR_DATA(dst))[idx3<<2|1];
						}
					}
				}
				if(kc!=1&&lossyconv_causalRCT[0])
					offset+=((int32_t*)HDR_DATA(dst))[idx<<2|1];

				lpred=0;
				wsum=0;
				for(kp=0;kp<16;++kp)
				{
					if(wpmask>>kp&1)
					{
						preds[kp]=0;
						uninit=1;
						for(idx2=0;idx2<25;++idx2)
						{
							short param=lossyconv_params[25*kp+idx2];
							preds[kp]+=(param>>1)*nb[idx2];
							if(param&1)
							{
								if(uninit)
								{
									uninit=0;
									vmin=nb[idx2], vmax=nb[idx2];
								}
								else
								{
									UPDATE_MIN(vmin, nb[idx2]);
									UPDATE_MAX(vmax, nb[idx2]);
								}
							}
						}
						if(!uninit)
						{
							MEDIAN3_32(preds[kp], vmin<<5, vmax<<5, preds[kp]);
						}
						weight=eNW[kp]+eN[kp]+eNE[kp]+eNEE[kp]+eWW[kp]+eW[kp];
						weight=0x1000000/(weight+1);
						lpred+=(int64_t)weight*preds[kp];
						wsum+=weight;
					}
				}
				lpred/=wsum+1;
				lpred+=16;
				lpred>>=5;
				pred=(int)lpred;
#if 0
				CLAMP3_32(pred, (int)lpred, nb[10], nb[7], nb[8]);/*clamp(N, W, NE)*/
#endif

				pred+=offset;
				CLAMP2_32(pred, pred, -half, half);
				{
					int curr=((int32_t*)HDR_DATA(dst))[idx<<2|kc];
					curr-=pred;
					curr<<=32-depth;/*signed-MA*/
					curr>>=32-depth;
					((int32_t*)HDR_DATA(src))[idx<<2|kc]=curr;
				}
				for(kp=0;kp<16;++kp)
				{
					if(wpmask>>kp&1)
						ecurr[kc]=abs(nb[12]-preds[kp]);
				}
				rows[0]+=16;
				rows[1]+=16;
				rows[2]+=16;
				rows[3]+=16;
			}
		}
	}
	free(dst);
	free(buf);
}
void filt_median33(Image *src)
{
	Image *dst;
	int kc, ky, idx, kx;

	dst=0;
	image_copy(&dst, src);
	for(kc=0;kc<4;++kc)
	{
		int depth;

		depth=src->depth[kc];
		if(!depth)
			continue;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[9]={0}, j, medians[3], pred;

#define LOAD(X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?((int32_t*)HDR_DATA(dst))[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				j=0;
				nb[j++]=LOAD(-1, -1);
				nb[j++]=LOAD(-1,  0);
				nb[j++]=LOAD(-1,  1);
				nb[j++]=LOAD( 0, -1);
				nb[j++]=LOAD( 0,  0);
				nb[j++]=LOAD( 0,  1);
				nb[j++]=LOAD( 1, -1);
				nb[j++]=LOAD( 1,  0);
				nb[j++]=LOAD( 1,  1);
#undef  LOAD

				MEDIAN3_32(medians[0], nb[0], nb[1], nb[2]);
				MEDIAN3_32(medians[1], nb[3], nb[4], nb[5]);
				MEDIAN3_32(medians[2], nb[6], nb[7], nb[8]);
				MEDIAN3_32(pred, medians[0], medians[1], medians[2]);

				((int32_t*)HDR_DATA(src))[idx<<2|kc]=pred;
			}
		}
	}
	free(dst);
}
void filt_av33(Image *src)
{
	Image *dst;
	int kc, ky, kx, idx, k;

	dst=0;
	image_copy(&dst, src);
	for(kc=0;kc<4;++kc)
	{
		int depth;

		depth=src->depth[kc];
		if(!depth)
			continue;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int pred, nb[9], j;

#define LOAD(X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?((int32_t*)HDR_DATA(dst))[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				j=0;
				nb[j++]=LOAD(-1, -1);
				nb[j++]=LOAD(-1,  0);
				nb[j++]=LOAD(-1,  1);
				nb[j++]=LOAD( 0, -1);
				nb[j++]=LOAD( 0,  0);
				nb[j++]=LOAD( 0,  1);
				nb[j++]=LOAD( 1, -1);
				nb[j++]=LOAD( 1,  0);
				nb[j++]=LOAD( 1,  1);
#undef  LOAD
#if 0
				if(kx==10&&ky==10)
					printf("");
#endif
				pred=0;
				for(k=0;k<(int)_countof(nb);++k)
					pred+=nb[k];
				pred/=(int)_countof(nb);

				((int32_t*)HDR_DATA(src))[idx<<2|kc]=pred;
			}
		}
	}
#if 0
	for(kc=0;kc<4;++kc)
	{
		int depth;

		depth=src->depth[kc];
		if(!depth)
			continue;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[9]={0}, j, pred;

#define LOAD(X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?((int32_t*)HDR_DATA(dst))[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				j=0;
				nb[j++]=  -LOAD(-1, -1),
				nb[j++]=  -LOAD(-1,  0),
				nb[j++]=  -LOAD(-1,  1),
				nb[j++]=  -LOAD( 0, -1),
				nb[j++]=17*LOAD( 0,  0),
				nb[j++]=  -LOAD( 0,  1),
				nb[j++]=  -LOAD( 1, -1),
				nb[j++]=  -LOAD( 1,  0),
				nb[j++]=  -LOAD( 1,  1),
#undef  LOAD
#if 0
				if(kx==10&&ky==10)
					printf("");
#endif
				pred=0;
				for(k=0;k<(int)_countof(nb);++k)
					pred+=nb[k];
				pred/=(int)_countof(nb);

				((int32_t*)HDR_DATA(src))[idx<<2|kc]=pred;
			}
		}
	}
#endif
	free(dst);
}