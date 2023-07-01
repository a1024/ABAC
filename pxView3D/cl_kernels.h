#ifdef _MSC_VER
#include<stdio.h>
#include<math.h>
#define __kernel
#define __global
#define __constant
#define get_global_id(X) X
#define ABS fabsf
#else
#define ABS fabs
#endif

#define CLAMP(LO, X, HI)    ((X)>(LO)?(X)<(HI)?(X):(HI):(LO))
#define COUNTOF(A) (sizeof(A)/sizeof(*(A)))
#define SGN(X) ((X)<0?-1:((X)>0?1:0))


	#define RESNET2
//	#define RECYCLE_PARAMS//writers juggle params between neighbor blocks	X  use V1_XBLOCKS==1 instead


#define LR 0.001f


//SYNC THIS WITH TRANSFORMS_GPU.C:

//	#define ZIGZAG_TRAVERSAL//X

#define V1_XBLOCKS 1//leave XBLOCKS at 1 for optimal performance
#define V1_YBLOCKS 8
#define V1_REACH 7
#define V1_NITER 7
#define V1_XYWRITERS (V1_XBLOCKS*V1_YBLOCKS)
#define V1_TOTALWRITERS (3*V1_XBLOCKS*V1_YBLOCKS)
#define V1_KERNELSIZE (V1_REACH*(V1_REACH+1)*2)
#define V1_NLANES (V1_TOTALWRITERS*V1_KERNELSIZE)

#define NF0_MAIN	55		//55
#define NF0_ERRORS	41		//41
#define NF0_PX		0		//41
#define NF0 (NF0_MAIN+NF0_ERRORS+NF0_PX)
#define NF1 16
#define NF2 16
#define NF3 96	//96
typedef struct TempsStruct
{
	float
		nb[NF0],
		net1[NF1], x1[NF1],
		net2[NF2], x2[NF2],
		net3[NF3], x3[NF3],
		pred, expr, loss;
} Temps;
typedef struct BwdTempsStruct
{
	float dL_dp, dL_dx3[NF3], dL_dx2[NF2], dL_dx1[NF1];
} BwdTemps;
typedef struct ParamsStruct
{
	float
		weight1[NF1*NF0], bias1[NF1],
		weight2[NF2*NF1], bias2[NF2],
		weight3[NF3*NF2], bias3[NF3],
		weight4[1*NF3];
} Params;
#define NPARAMS (sizeof(Params)/sizeof(float))

__constant int xoffsets_pos[]=
{
	                   0,
	                   0,
	               -3, 0, 3,
	   -5,-4,-3,-2,-1, 0, 1, 2, 3, 4,
	         -3,-2,-1, 0, 1, 2, 3, 4,
	      -4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
	-6,-5,-4,-3,-2,-1,
};
__constant int xoffsets_neg[]=
{
	                   0,
	                   0,
	                3, 0,-3,
	    5, 4, 3, 2, 1, 0,-1,-2,-3,-4,
	          3, 2, 1, 0,-1,-2,-3,-4,
	       4, 3, 2, 1, 0,-1,-2,-3,-4,-5,-6,-7,
	 6, 5, 4, 3, 2, 1,
};
static float clamp4(float p, float a, float b, float c, float d)
{
	float vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmin>c)vmin=c;
	if(vmin>d)vmin=d;
	if(vmax<b)vmax=b;
	if(vmax<c)vmax=c;
	if(vmax<d)vmax=d;
	p=CLAMP(vmin, p, vmax);
	return p;
}
static float clip(float x)
{
	x=CLAMP(-1, x, 1);
	return x;
}
static void get_nb2(__global const float *buf, __global const float *errors, int iw, int ih, int kc, int kx, int ky, int xstart, int xend, int ystart, int yend, int flip_east_west, __global float *ctx)
{
	__constant int *offsets=flip_east_west?xoffsets_pos:xoffsets_neg;
	unsigned dx=xend-xstart, dy=yend-ystart;
	
	//if(xstart>=xend||ystart>=yend)
	//if(kx-xstart==32&&ky-ystart==32)
	//{
	//	printf("get_nb2() CXY %d %3d %3d  X [%3d, %3d], Y [%3d, %3d]\n", kc, kx, ky, xstart, xend, ystart, yend);
	//}

#define LOAD(BUF, C, X, Y) ((X)>=xstart&&(X)<xend&&(Y)>=ystart&&(Y)<yend?BUF[iw*(ih*C+Y)+X]:0)
//#define LOAD(BUF, C, X, Y) ((unsigned)(X)<iw&&(unsigned)(Y)<ih?BUF[iw*(ih*C+Y)+X]:0)
//#define LOAD(BUF, C, X, Y) ((unsigned)(X-xstart)<dx&&(unsigned)(Y-ystart)<dy?BUF[iw*(ih*C+Y)+X]:0)
	float
		NNNNNN  =LOAD(buf, kc, kx+offsets[ 0], ky-6),
		NNNNN   =LOAD(buf, kc, kx+offsets[ 1], ky-5),
		NNNNWWW =LOAD(buf, kc, kx+offsets[ 2], ky-4),
		NNNN    =LOAD(buf, kc, kx+offsets[ 3], ky-4),
		NNNNEEE =LOAD(buf, kc, kx+offsets[ 4], ky-4),
		NNNWWWWW=LOAD(buf, kc, kx+offsets[ 5], ky-3),
		NNNWWWW =LOAD(buf, kc, kx+offsets[ 6], ky-3),
		NNNWWW  =LOAD(buf, kc, kx+offsets[ 7], ky-3),
		NNNWW   =LOAD(buf, kc, kx+offsets[ 8], ky-3),
		NNNW    =LOAD(buf, kc, kx+offsets[ 9], ky-3),
		NNN     =LOAD(buf, kc, kx+offsets[10], ky-3),
		NNNE    =LOAD(buf, kc, kx+offsets[11], ky-3),
		NNNEE   =LOAD(buf, kc, kx+offsets[12], ky-3),
		NNNEEE  =LOAD(buf, kc, kx+offsets[13], ky-3),
		NNNEEEE =LOAD(buf, kc, kx+offsets[14], ky-3),
		NNWWW   =LOAD(buf, kc, kx+offsets[15], ky-2),
		NNWW    =LOAD(buf, kc, kx+offsets[16], ky-2),
		NNW     =LOAD(buf, kc, kx+offsets[17], ky-2),
		NN      =LOAD(buf, kc, kx+offsets[18], ky-2),
		NNE     =LOAD(buf, kc, kx+offsets[19], ky-2),
		NNEE    =LOAD(buf, kc, kx+offsets[20], ky-2),
		NNEEE   =LOAD(buf, kc, kx+offsets[21], ky-2),
		NNEEEE  =LOAD(buf, kc, kx+offsets[22], ky-2),
		NWWWW   =LOAD(buf, kc, kx+offsets[23], ky-1),
		NWWW    =LOAD(buf, kc, kx+offsets[24], ky-1),
		NWW     =LOAD(buf, kc, kx+offsets[25], ky-1),
		NW      =LOAD(buf, kc, kx+offsets[26], ky-1),
		N       =LOAD(buf, kc, kx+offsets[27], ky-1),
		NE      =LOAD(buf, kc, kx+offsets[28], ky-1),
		NEE     =LOAD(buf, kc, kx+offsets[29], ky-1),
		NEEE    =LOAD(buf, kc, kx+offsets[30], ky-1),
		NEEEE   =LOAD(buf, kc, kx+offsets[31], ky-1),
		NEEEEE  =LOAD(buf, kc, kx+offsets[32], ky-1),
		NEEEEEE =LOAD(buf, kc, kx+offsets[33], ky-1),
		NEEEEEEE=LOAD(buf, kc, kx+offsets[34], ky-1),
		WWWWWW  =LOAD(buf, kc, kx+offsets[35], ky  ),
		WWWWW   =LOAD(buf, kc, kx+offsets[36], ky  ),
		WWWW    =LOAD(buf, kc, kx+offsets[37], ky  ),
		WWW     =LOAD(buf, kc, kx+offsets[38], ky  ),
		WW      =LOAD(buf, kc, kx+offsets[39], ky  ),
		W       =LOAD(buf, kc, kx+offsets[40], ky  );
#if 0
	float
		NNNNNN  =LOAD(buf, kc, kx  , ky-6),
		NNNNN   =LOAD(buf, kc, kx  , ky-5),

		NNNNWWW =LOAD(buf, kc, kx-3, ky-4),
		NNNN    =LOAD(buf, kc, kx  , ky-4),
		NNNNEEE =LOAD(buf, kc, kx+3, ky-4),

		NNNWWWWW=LOAD(buf, kc, kx-5, ky-3),
		NNNWWWW =LOAD(buf, kc, kx-4, ky-3),
		NNNWWW  =LOAD(buf, kc, kx-3, ky-3),
		NNNWW   =LOAD(buf, kc, kx-2, ky-3),
		NNNW    =LOAD(buf, kc, kx-1, ky-3),
		NNN     =LOAD(buf, kc, kx  , ky-3),
		NNNE    =LOAD(buf, kc, kx+1, ky-3),
		NNNEE   =LOAD(buf, kc, kx+2, ky-3),
		NNNEEE  =LOAD(buf, kc, kx+3, ky-3),
		NNNEEEE =LOAD(buf, kc, kx+4, ky-3),

		NNWWW   =LOAD(buf, kc, kx-3, ky-2),
		NNWW    =LOAD(buf, kc, kx-2, ky-2),
		NNW     =LOAD(buf, kc, kx-1, ky-2),
		NN      =LOAD(buf, kc, kx  , ky-2),
		NNE     =LOAD(buf, kc, kx+1, ky-2),
		NNEE    =LOAD(buf, kc, kx+2, ky-2),
		NNEEE   =LOAD(buf, kc, kx+3, ky-2),
		NNEEEE  =LOAD(buf, kc, kx+4, ky-2),
		
		NWWWW   =LOAD(buf, kc, kx-4, ky-1),
		NWWW    =LOAD(buf, kc, kx-3, ky-1),
		NWW     =LOAD(buf, kc, kx-2, ky-1),
		NW      =LOAD(buf, kc, kx-1, ky-1),
		N       =LOAD(buf, kc, kx  , ky-1),
		NE      =LOAD(buf, kc, kx+1, ky-1),
		NEE     =LOAD(buf, kc, kx+2, ky-1),
		NEEE    =LOAD(buf, kc, kx+3, ky-1),
		NEEEE   =LOAD(buf, kc, kx+4, ky-1),
		NEEEEE  =LOAD(buf, kc, kx+5, ky-1),
		NEEEEEE =LOAD(buf, kc, kx+6, ky-1),
		NEEEEEEE=LOAD(buf, kc, kx+7, ky-1),

		WWWWWW  =LOAD(buf, kc, kx-6, ky),
		WWWWW   =LOAD(buf, kc, kx-5, ky),
		WWWW    =LOAD(buf, kc, kx-4, ky),
		WWW     =LOAD(buf, kc, kx-3, ky),
		WW      =LOAD(buf, kc, kx-2, ky),
		W       =LOAD(buf, kc, kx-1, ky);
#endif
	
#if NF0_MAIN==55
	*ctx++ = clamp4(W + N - NW, W, NW, N, NE);//0
	*ctx++ = clip(W + N - NW);//1
	*ctx++ = clamp4(W + NE - N, W, NW, N, NE);//2
	*ctx++ = clip(W + NE - N);//3
	*ctx++ = clamp4(N + NW - NNW, W, NW, N, NE);//4
	*ctx++ = clip(N + NW - NNW);//5
	*ctx++ = clamp4(N + NE - NNE, W, N, NE, NEE);//6
	*ctx++ = clip(N + NE - NNE);//7
	*ctx++ = (W + NEE) / 2;//8
	*ctx++ = clip(N * 3 - NN * 3 + NNN);//9
	*ctx++ = clip(W * 3 - WW * 3 + WWW);//10
	*ctx++ = (W + clip(NE * 3 - NNE * 3 + NNNE)) / 2;//11
	*ctx++ = (W + clip(NEE * 3 - NNEEE * 3 + NNNEEEE)) / 2;//12
	*ctx++ = clip(NN + NNNN - NNNNNN);//13
	*ctx++ = clip(WW + WWWW - WWWWWW);//14
	*ctx++ = clip((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 2 - NWW, W, NW, N, NN)) / 6);//15
	*ctx++ = clip((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);//16
	*ctx++ = clip(NN + NW - NNNW);//17
	*ctx++ = clip(NN + NE - NNNE);//18
	*ctx++ = clip((W * 2 + NW) - (WW + 2 * NWW) + NWWW);//19
	*ctx++ = clip(((NW + NWW) / 2) * 3 - NNWWW * 3 + (NNNWWWW + NNNWWWWW) / 2);//20
	*ctx++ = clip(NEE + NE - NNEEE);//21
	*ctx++ = clip(NWW + WW - NWWWW);//22
	*ctx++ = clip(((W + NW) * 3 - NWW * 6 + NWWW + NNWWW) / 2);//23
	*ctx++ = clip((NE * 2 + NNE) - (NNEE + NNNEE * 2) + NNNNEEE);//24
	*ctx++ = NNNNNN;//25
	*ctx++ = (NEEEE + NEEEEEE) / 2;//26
	*ctx++ = (WWWW + WWWWWW) / 2;//27
	*ctx++ = (W + N + NEEEEE + NEEEEEEE) / 4;//28
	*ctx++ = clip(NEEE + W - NEE);//29
	*ctx++ = clip(4 * NNN - 3 * NNNN);//30
	*ctx++ = clip(N + NN - NNN);//31
	*ctx++ = clip(W + WW - WWW);//32
	*ctx++ = clip(W + NEE - NE);//33
	*ctx++ = clip(WW + NEE - N);//34
	*ctx++ = (clip(W * 2 - NW) + clip(W * 2 - NWW) + N + NE) / 4;//35
	*ctx++ = clamp4(N * 2 - NN, W, N, NE, NEE);//36
	*ctx++ = (N + NNN) / 2;//37
	*ctx++ = clip(NN + W - NNW);//38
	*ctx++ = clip(NWW + N - NNWW);//39
	*ctx++ = clip((4 * WWW - 15 * WW + 20 * W + clip(NEE * 2 - NNEE)) / 10);//40
	*ctx++ = clip((NNNEEE - 4 * NNEE + 6 * NE + clip(W * 3 - NW * 3 + NNW)) / 4);//41
	*ctx++ = clip((N * 2 + NE) - (NN + 2 * NNE) + NNNE);//42
	*ctx++ = clip((NW * 2 + NNW) - (NNWW + NNNWW * 2) + NNNNWWW);//43
	*ctx++ = clip(NNWW + W - NNWWW);//44
	*ctx++ = clip((-NNNN + 5 * NNN - 10 * NN + 10 * N + clip(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW)) / 5);//45
	*ctx++ = clip(NEE + clip(NEEE * 2 - NNEEEE) - NEEEE);//46
	*ctx++ = clip(NW + W - NWW);//47
	*ctx++ = clip((N * 2 + NW) - (NN + 2 * NNW) + NNNW);//48
	*ctx++ = clip(NN + clip(NEE * 2 - NNEEE) - NNE);//49
	*ctx++ = clip((-WWWW + 5 * WWW - 10 * WW + 10 * W + clip(NE * 2 - NNE)) / 5);//50
	*ctx++ = clip((-WWWWW + 4 * WWWW - 5 * WWW + 5 * W + clip(NE * 2 - NNE)) / 4);//51
	*ctx++ = clip((WWW - 4 * WW + 6 * W + clip(NE * 3 - NNE * 3 + NNNE)) / 4);//52
	*ctx++ = clip((-NNEE + 3 * NE + clip(W * 4 - NW * 6 + NNW * 4 - NNNW)) / 3);//53
	*ctx++ = ((W + N) * 3 - NW * 2) / 4;//54
#endif
#if NF0_ERRORS==41
	*ctx++ = LOAD(errors, kc, kx+offsets[ 0], ky-6);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 1], ky-5);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 2], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 3], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 4], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 5], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 6], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 7], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 8], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 9], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[10], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[11], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[12], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[13], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[14], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[15], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[16], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[17], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[18], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[19], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[20], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[21], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[22], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[23], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[24], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[25], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[26], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[27], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[28], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[29], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[30], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[31], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[32], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[33], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[34], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[35], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[36], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[37], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[38], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[39], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[40], ky  );
#endif
#if NF0_PX==41
	*ctx++ = LOAD(buf, kc, kx+offsets[ 0], ky-6);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 1], ky-5);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 2], ky-4);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 3], ky-4);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 4], ky-4);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 5], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 6], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 7], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 8], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[ 9], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[10], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[11], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[12], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[13], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[14], ky-3);
	*ctx++ = LOAD(buf, kc, kx+offsets[15], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[16], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[17], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[18], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[19], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[20], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[21], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[22], ky-2);
	*ctx++ = LOAD(buf, kc, kx+offsets[23], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[24], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[25], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[26], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[27], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[28], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[29], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[30], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[31], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[32], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[33], ky-1);
	*ctx++ = LOAD(buf, kc, kx+offsets[34], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[35], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[36], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[37], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[38], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[39], ky  );
	*ctx++ = LOAD(buf, kc, kx+offsets[40], ky  );
#endif
}


static void add_vec(__global float *dst, __global const float *a, __global const float *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b[k];
}
static void act(__global float *dst, __global const float *src, int count)
{
	for(int k=0;k<count;++k)
	{
		float val=src[k];
		float negpart=val*0.01f;
		val=val>negpart?val:negpart;
		val=CLAMP(-10, val, 10);
		dst[k]=val;
	}
}
static void act_dash(__global float *dst, __global const float *src, int count)
{
	for(int k=0;k<count;++k)
	{
		float val=src[k];
		if(val<-10||val>10)
			val=0;
		else
			val=val<0?0.01f:1;
		dst[k]=val;
	}
}
static void mul_vec_scalar(__global float *dst, __global float *vec, float scalar, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=vec[k]*scalar;
}
static void mul_vec_ew(__global float *dst, __global const float *v1, __global const float *v2, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=v1[k]*v2[k];
}
static void mul_vec_outer(__global float *dst, __global const float *left, __global const float *right, int lh, int rw)
{
	for(int ky=0;ky<lh;++ky)
	{
		for(int kx=0;kx<rw;++kx)
			dst[rw*ky+kx]=left[ky]*right[kx];
	}
}
static void mul_mat(__global float *dst, __global const float *m1, __global const float *m2, int h1, int w1h2, int w2)
{
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			float sum=0;
			for(int k=0;k<w1h2;++k)
				sum+=m1[w1h2*ky+k]*m2[w2*k+kx];
			dst[w2*ky+kx]=sum;
		}
	}
}
static void linear(__global float *dst, __global const float *mat, __global const float *vec, __global const float *bias, int win, int hout)//fixed 15.16 bit
{
	for(int ko=0;ko<hout;++ko)
	{
		float temp=bias?bias[ko]:0;
		for(int ki=0;ki<win;++ki)
			temp+=mat[win*ko+ki]*vec[ki];
		dst[ko]=temp;
	}
}


//fwd
//[NF1]	x1 = act(w1*nb + b1)		//+ nb
//[NF2]	x2 = act(w2*x1 + b2)		//+ x1
//[NF3]	x3 = act(w3*x2 + b3) + nb	//+ x2
//[1]	pred = w4*x3
//[1]	L = abs(pred-val)
//
//bwd
//[1]		dL_dp = sgn(p-x)
//[NF3]		g_w4 = dL_dp * x3
//
//[NF3]		dL_dx3 = dL_dp * w4
//[NF3]		g_b3 = dL_dx3 .* act'(n3)
//[NF3*NF2]	g_w3 = g_b3 o x2
//
//[NF2]		dL_dx2 = dL_dx3 * (act'(n3) * w3 [+ I])		= g_b3 * w3 + dL_dx3
//[NF2]		g_b2 = dL_dx2 * act'(n2)
//[NF2*NF1]	g_w2 = g_b2 o x1
//
//[NF1]		dL_dx1 = dL_dx2 * (act'(n2) * w2 [+ I])		= g_b2 * w2 + dL_dx2
//[NF1]		g_b1 = dL_dx1 * act'(n1)
//[NF1*NF0]	g_w1 = g_b1 o nb
static void eval_fwd_v1(__global float *pixels, __global float *errors, int iw, int ih, int kc, int kx, int ky, int xstart, int xend, int ystart, int yend, int flip_east_west, __global Params *p, __global Temps *t)
{
	__constant int *offsets=flip_east_west?xoffsets_pos:xoffsets_neg;
	unsigned dx=xend-xstart, dy=yend-ystart;
	__global float *ctx=t->nb;
	
	//if(xstart>=xend||ystart>=yend)
	//if(kx-xstart==32&&ky-ystart==32)
	//{
	//	printf("get_nb2() CXY %d %3d %3d  X [%3d, %3d], Y [%3d, %3d]\n", kc, kx, ky, xstart, xend, ystart, yend);
	//}

#define LOAD(BUF, C, X, Y) ((X)>=xstart&&(X)<xend&&(Y)>=ystart&&(Y)<yend?BUF[iw*(ih*C+Y)+X]:0)
//#define LOAD(BUF, C, X, Y) ((unsigned)(X)<iw&&(unsigned)(Y)<ih?BUF[iw*(ih*C+Y)+X]:0)
//#define LOAD(BUF, C, X, Y) ((unsigned)(X-xstart)<dx&&(unsigned)(Y-ystart)<dy?BUF[iw*(ih*C+Y)+X]:0)
#if 1
	float
		NNNNNN  =LOAD(pixels, kc, kx+offsets[ 0], ky-6),
		NNNNN   =LOAD(pixels, kc, kx+offsets[ 1], ky-5),
		NNNNWWW =LOAD(pixels, kc, kx+offsets[ 2], ky-4),
		NNNN    =LOAD(pixels, kc, kx+offsets[ 3], ky-4),
		NNNNEEE =LOAD(pixels, kc, kx+offsets[ 4], ky-4),
		NNNWWWWW=LOAD(pixels, kc, kx+offsets[ 5], ky-3),
		NNNWWWW =LOAD(pixels, kc, kx+offsets[ 6], ky-3),
		NNNWWW  =LOAD(pixels, kc, kx+offsets[ 7], ky-3),
		NNNWW   =LOAD(pixels, kc, kx+offsets[ 8], ky-3),
		NNNW    =LOAD(pixels, kc, kx+offsets[ 9], ky-3),
		NNN     =LOAD(pixels, kc, kx+offsets[10], ky-3),
		NNNE    =LOAD(pixels, kc, kx+offsets[11], ky-3),
		NNNEE   =LOAD(pixels, kc, kx+offsets[12], ky-3),
		NNNEEE  =LOAD(pixels, kc, kx+offsets[13], ky-3),
		NNNEEEE =LOAD(pixels, kc, kx+offsets[14], ky-3),
		NNWWW   =LOAD(pixels, kc, kx+offsets[15], ky-2),
		NNWW    =LOAD(pixels, kc, kx+offsets[16], ky-2),
		NNW     =LOAD(pixels, kc, kx+offsets[17], ky-2),
		NN      =LOAD(pixels, kc, kx+offsets[18], ky-2),
		NNE     =LOAD(pixels, kc, kx+offsets[19], ky-2),
		NNEE    =LOAD(pixels, kc, kx+offsets[20], ky-2),
		NNEEE   =LOAD(pixels, kc, kx+offsets[21], ky-2),
		NNEEEE  =LOAD(pixels, kc, kx+offsets[22], ky-2),
		NWWWW   =LOAD(pixels, kc, kx+offsets[23], ky-1),
		NWWW    =LOAD(pixels, kc, kx+offsets[24], ky-1),
		NWW     =LOAD(pixels, kc, kx+offsets[25], ky-1),
		NW      =LOAD(pixels, kc, kx+offsets[26], ky-1),
		N       =LOAD(pixels, kc, kx+offsets[27], ky-1),
		NE      =LOAD(pixels, kc, kx+offsets[28], ky-1),
		NEE     =LOAD(pixels, kc, kx+offsets[29], ky-1),
		NEEE    =LOAD(pixels, kc, kx+offsets[30], ky-1),
		NEEEE   =LOAD(pixels, kc, kx+offsets[31], ky-1),
		NEEEEE  =LOAD(pixels, kc, kx+offsets[32], ky-1),
		NEEEEEE =LOAD(pixels, kc, kx+offsets[33], ky-1),
		NEEEEEEE=LOAD(pixels, kc, kx+offsets[34], ky-1),
		WWWWWW  =LOAD(pixels, kc, kx+offsets[35], ky  ),
		WWWWW   =LOAD(pixels, kc, kx+offsets[36], ky  ),
		WWWW    =LOAD(pixels, kc, kx+offsets[37], ky  ),
		WWW     =LOAD(pixels, kc, kx+offsets[38], ky  ),
		WW      =LOAD(pixels, kc, kx+offsets[39], ky  ),
		W       =LOAD(pixels, kc, kx+offsets[40], ky  );
#endif
#if 0
	float
		NNNNNN  =LOAD(pixels, kc, kx  , ky-6),
		NNNNN   =LOAD(pixels, kc, kx  , ky-5),

		NNNNWWW =LOAD(pixels, kc, kx-3, ky-4),
		NNNN    =LOAD(pixels, kc, kx  , ky-4),
		NNNNEEE =LOAD(pixels, kc, kx+3, ky-4),

		NNNWWWWW=LOAD(pixels, kc, kx-5, ky-3),
		NNNWWWW =LOAD(pixels, kc, kx-4, ky-3),
		NNNWWW  =LOAD(pixels, kc, kx-3, ky-3),
		NNNWW   =LOAD(pixels, kc, kx-2, ky-3),
		NNNW    =LOAD(pixels, kc, kx-1, ky-3),
		NNN     =LOAD(pixels, kc, kx  , ky-3),
		NNNE    =LOAD(pixels, kc, kx+1, ky-3),
		NNNEE   =LOAD(pixels, kc, kx+2, ky-3),
		NNNEEE  =LOAD(pixels, kc, kx+3, ky-3),
		NNNEEEE =LOAD(pixels, kc, kx+4, ky-3),

		NNWWW   =LOAD(pixels, kc, kx-3, ky-2),
		NNWW    =LOAD(pixels, kc, kx-2, ky-2),
		NNW     =LOAD(pixels, kc, kx-1, ky-2),
		NN      =LOAD(pixels, kc, kx  , ky-2),
		NNE     =LOAD(pixels, kc, kx+1, ky-2),
		NNEE    =LOAD(pixels, kc, kx+2, ky-2),
		NNEEE   =LOAD(pixels, kc, kx+3, ky-2),
		NNEEEE  =LOAD(pixels, kc, kx+4, ky-2),
		
		NWWWW   =LOAD(pixels, kc, kx-4, ky-1),
		NWWW    =LOAD(pixels, kc, kx-3, ky-1),
		NWW     =LOAD(pixels, kc, kx-2, ky-1),
		NW      =LOAD(pixels, kc, kx-1, ky-1),
		N       =LOAD(pixels, kc, kx  , ky-1),
		NE      =LOAD(pixels, kc, kx+1, ky-1),
		NEE     =LOAD(pixels, kc, kx+2, ky-1),
		NEEE    =LOAD(pixels, kc, kx+3, ky-1),
		NEEEE   =LOAD(pixels, kc, kx+4, ky-1),
		NEEEEE  =LOAD(pixels, kc, kx+5, ky-1),
		NEEEEEE =LOAD(pixels, kc, kx+6, ky-1),
		NEEEEEEE=LOAD(pixels, kc, kx+7, ky-1),

		WWWWWW  =LOAD(pixels, kc, kx-6, ky),
		WWWWW   =LOAD(pixels, kc, kx-5, ky),
		WWWW    =LOAD(pixels, kc, kx-4, ky),
		WWW     =LOAD(pixels, kc, kx-3, ky),
		WW      =LOAD(pixels, kc, kx-2, ky),
		W       =LOAD(pixels, kc, kx-1, ky);
#endif
	
#if NF0_MAIN==55
	*ctx++ = clamp4(W + N - NW, W, NW, N, NE);//0
	*ctx++ = clip(W + N - NW);//1
	*ctx++ = clamp4(W + NE - N, W, NW, N, NE);//2
	*ctx++ = clip(W + NE - N);//3
	*ctx++ = clamp4(N + NW - NNW, W, NW, N, NE);//4
	*ctx++ = clip(N + NW - NNW);//5
	*ctx++ = clamp4(N + NE - NNE, W, N, NE, NEE);//6
	*ctx++ = clip(N + NE - NNE);//7
	*ctx++ = (W + NEE) / 2;//8
	*ctx++ = clip(N * 3 - NN * 3 + NNN);//9
	*ctx++ = clip(W * 3 - WW * 3 + WWW);//10
	*ctx++ = (W + clip(NE * 3 - NNE * 3 + NNNE)) / 2;//11
	*ctx++ = (W + clip(NEE * 3 - NNEEE * 3 + NNNEEEE)) / 2;//12
	*ctx++ = clip(NN + NNNN - NNNNNN);//13
	*ctx++ = clip(WW + WWWW - WWWWWW);//14
	*ctx++ = clip((NNNNN - 6 * NNNN + 15 * NNN - 20 * NN + 15 * N + clamp4(W * 2 - NWW, W, NW, N, NN)) / 6);//15
	*ctx++ = clip((-3 * WW + 8 * W + clamp4(NEE * 3 - NNEE * 3 + NNNEE, NE, NEE, NEEE, NEEEE)) / 6);//16
	*ctx++ = clip(NN + NW - NNNW);//17
	*ctx++ = clip(NN + NE - NNNE);//18
	*ctx++ = clip((W * 2 + NW) - (WW + 2 * NWW) + NWWW);//19
	*ctx++ = clip(((NW + NWW) / 2) * 3 - NNWWW * 3 + (NNNWWWW + NNNWWWWW) / 2);//20
	*ctx++ = clip(NEE + NE - NNEEE);//21
	*ctx++ = clip(NWW + WW - NWWWW);//22
	*ctx++ = clip(((W + NW) * 3 - NWW * 6 + NWWW + NNWWW) / 2);//23
	*ctx++ = clip((NE * 2 + NNE) - (NNEE + NNNEE * 2) + NNNNEEE);//24
	*ctx++ = NNNNNN;//25
	*ctx++ = (NEEEE + NEEEEEE) / 2;//26
	*ctx++ = (WWWW + WWWWWW) / 2;//27
	*ctx++ = (W + N + NEEEEE + NEEEEEEE) / 4;//28
	*ctx++ = clip(NEEE + W - NEE);//29
	*ctx++ = clip(4 * NNN - 3 * NNNN);//30
	*ctx++ = clip(N + NN - NNN);//31
	*ctx++ = clip(W + WW - WWW);//32
	*ctx++ = clip(W + NEE - NE);//33
	*ctx++ = clip(WW + NEE - N);//34
	*ctx++ = (clip(W * 2 - NW) + clip(W * 2 - NWW) + N + NE) / 4;//35
	*ctx++ = clamp4(N * 2 - NN, W, N, NE, NEE);//36
	*ctx++ = (N + NNN) / 2;//37
	*ctx++ = clip(NN + W - NNW);//38
	*ctx++ = clip(NWW + N - NNWW);//39
	*ctx++ = clip((4 * WWW - 15 * WW + 20 * W + clip(NEE * 2 - NNEE)) / 10);//40
	*ctx++ = clip((NNNEEE - 4 * NNEE + 6 * NE + clip(W * 3 - NW * 3 + NNW)) / 4);//41
	*ctx++ = clip((N * 2 + NE) - (NN + 2 * NNE) + NNNE);//42
	*ctx++ = clip((NW * 2 + NNW) - (NNWW + NNNWW * 2) + NNNNWWW);//43
	*ctx++ = clip(NNWW + W - NNWWW);//44
	*ctx++ = clip((-NNNN + 5 * NNN - 10 * NN + 10 * N + clip(W * 4 - NWW * 6 + NNWWW * 4 - NNNWWWW)) / 5);//45
	*ctx++ = clip(NEE + clip(NEEE * 2 - NNEEEE) - NEEEE);//46
	*ctx++ = clip(NW + W - NWW);//47
	*ctx++ = clip((N * 2 + NW) - (NN + 2 * NNW) + NNNW);//48
	*ctx++ = clip(NN + clip(NEE * 2 - NNEEE) - NNE);//49
	*ctx++ = clip((-WWWW + 5 * WWW - 10 * WW + 10 * W + clip(NE * 2 - NNE)) / 5);//50
	*ctx++ = clip((-WWWWW + 4 * WWWW - 5 * WWW + 5 * W + clip(NE * 2 - NNE)) / 4);//51
	*ctx++ = clip((WWW - 4 * WW + 6 * W + clip(NE * 3 - NNE * 3 + NNNE)) / 4);//52
	*ctx++ = clip((-NNEE + 3 * NE + clip(W * 4 - NW * 6 + NNW * 4 - NNNW)) / 3);//53
	*ctx++ = ((W + N) * 3 - NW * 2) / 4;//54
#endif
#if NF0_ERRORS==41
	*ctx++ = LOAD(errors, kc, kx+offsets[ 0], ky-6);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 1], ky-5);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 2], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 3], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 4], ky-4);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 5], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 6], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 7], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 8], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[ 9], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[10], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[11], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[12], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[13], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[14], ky-3);
	*ctx++ = LOAD(errors, kc, kx+offsets[15], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[16], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[17], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[18], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[19], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[20], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[21], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[22], ky-2);
	*ctx++ = LOAD(errors, kc, kx+offsets[23], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[24], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[25], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[26], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[27], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[28], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[29], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[30], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[31], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[32], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[33], ky-1);
	*ctx++ = LOAD(errors, kc, kx+offsets[34], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[35], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[36], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[37], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[38], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[39], ky  );
	*ctx++ = LOAD(errors, kc, kx+offsets[40], ky  );
#endif
#if NF0_PX==41
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 0], ky-6);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 1], ky-5);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 2], ky-4);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 3], ky-4);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 4], ky-4);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 5], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 6], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 7], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 8], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[ 9], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[10], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[11], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[12], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[13], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[14], ky-3);
	*ctx++ = LOAD(pixels, kc, kx+offsets[15], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[16], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[17], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[18], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[19], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[20], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[21], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[22], ky-2);
	*ctx++ = LOAD(pixels, kc, kx+offsets[23], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[24], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[25], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[26], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[27], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[28], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[29], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[30], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[31], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[32], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[33], ky-1);
	*ctx++ = LOAD(pixels, kc, kx+offsets[34], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[35], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[36], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[37], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[38], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[39], ky  );
	*ctx++ = LOAD(pixels, kc, kx+offsets[40], ky  );
#endif

	//if(xstart>=xend||ystart>=yend)
	//if(kx-xstart==32&&ky-ystart==32)
	//{
	//	printf("eval_fwd_v1() CXY %d %3d %3d  X [%3d, %3d], Y [%3d, %3d]\n", kc, kx, ky, xstart, xend, ystart, yend);
	//}
	//get_nb2(pixels, errors, iw, ih, kc, kx, ky, xstart, xend, ystart, xend, flip_east_west, t->nb);
	
	//t->pred=t->nb[0];//
	//return;//

	linear(t->net1, p->weight1, t->nb, p->bias1, NF0, NF1);	//n1 = w1*nb + b1
	act(t->x1, t->net1, NF1);								//x1 = act(n1)
#ifdef RESNET
	add_vec(t->x1, t->x1, t->nb, NF0);						//x1r = x1 + nb		= act(w1*nb + b1) + nb												dx1r_dw1 = act'(n1)*nb		dx1r_db1 = act'(n1)
#endif

	linear(t->net2, p->weight2, t->x1, p->bias2, NF1, NF2);//n2 = w2*x1r + b2
	act(t->x2, t->net2, NF2);								//x2 = act(n2)
#ifdef RESNET
	add_vec(t->x2, t->x2, t->x1, NF0);						//x2r = x2 + x1r	= act(w2*x1r + b2) + x1r		dx2r_dx1r = act'(n2)*w2 + I			dx2r_dw2 = act'(n2)*x1r		dx2r_db2 = act'(n2)
#endif

	linear(t->net3, p->weight3, t->x2, p->bias3, NF2, NF3);	//n3 = w3*x2r + b3
	act(t->x3, t->net3, NF3);								//x3 = act(n3)
#ifdef RESNET2
	add_vec(t->x3, t->x3, t->nb, NF0);
#endif
#ifdef RESNET
	add_vec(t->x3, t->x3, t->x2, NF0);						//x3r = x3 + x2r	= act(w3*x2r + b3) + x2r		dx3r_dx2r = act'(n3)*w3 + I			dx3r_dw3 = act'(n3)*x2r		dx3r_db3 = act'(n3)		<- read right to left
#endif

	linear(&t->pred, p->weight4, t->x3, 0, NF3, 1);			//pred = w4*x3r				dL_dp = sgn(p-x)		dp_dx3r = w4						dp_dw4 = x3r
}
static void eval_bwd_v1(__global float *pixels, __global float *errors, int iw, int ih, int kc, int kx, int ky, int fwd, __global Params *p, __global Temps *t, __global BwdTemps *b, __global Params *g)
{
	float curr=pixels[iw*(ih*kc+ky)+kx];
	t->expr=t->pred-curr;
	t->loss=ABS(t->expr);
	
	//bwd
	b->dL_dp=SGN(t->expr);//dL_dp = sgn(p-x)
	mul_vec_scalar(g->weight4, t->x3, b->dL_dp, NF3);//dL_dw4 = dL_dp * dp_dw4		= sgn(p-x) * x3T


	mul_vec_scalar(b->dL_dx3, p->weight4, b->dL_dp, NF3);
	act_dash(g->bias3, t->net3, NF3);
	mul_vec_ew(g->bias3, g->bias3, b->dL_dx3, NF3);
	mul_vec_outer(g->weight3, g->bias3, t->x2, NF3, NF2);


	mul_mat(b->dL_dx2, g->bias3, p->weight3, 1, NF3, NF2);
#ifdef RESNET
	add_vec(b->dL_dx2, b->dL_dx2, b->dL_dx3, NF2);
#endif
	act_dash(g->bias2, t->net2, NF2);
	mul_vec_ew(g->bias2, g->bias2, b->dL_dx2, NF2);
	mul_vec_outer(g->weight2, g->bias2, t->x1, NF2, NF1);


	mul_mat(b->dL_dx1, g->bias2, p->weight2, 1, NF2, NF1);
#ifdef RESNET
	add_vec(b->dL_dx1, b->dL_dx1, b->dL_dx2, NF1);
#endif
	act_dash(g->bias1, t->net1, NF1);
	mul_vec_ew(g->bias1, g->bias1, b->dL_dx1, NF1);
	mul_vec_outer(g->weight1, g->bias1, t->nb, NF1, NF0);
}

//worksize: V1_NLANES				indices[3]={iw, ih, kx, ky, fwd}, params[V1_TOTALWRITERS], temps[V1_NLANES], bwdtemps[V1_NLANES], gradients[V1_NLANES]
__kernel void train(__global float *pixels, __global float *errors, __constant int *indices, __global Params *params, __global Temps *temps, __global BwdTemps *bwdtemps, __global Params *gradients)
{
	int threadidx=get_global_id(0);//threadidx = V1_KERNELSIZE*(V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xblockidx) + V1_KERNELWIDTH*unsigned_ry + unsigned_rx			note: V1_KERNELSIZE is not divisible by V1_KERNELWIDTH
	int kernelidx=threadidx%V1_KERNELSIZE,
		writeridx=threadidx/V1_KERNELSIZE,//writeridx = V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xblockidx

		rx=kernelidx%(V1_REACH<<1|1),
		ry=kernelidx/(V1_REACH<<1|1),

		xblockidx=writeridx%V1_XBLOCKS,
		rowidx   =writeridx/V1_XBLOCKS,
		yblockidx=rowidx%V1_YBLOCKS,
		kc       =rowidx/V1_YBLOCKS;

#if 0
	int threadidx=get_global_id(0);//threadidx = V1_KERNELSIZE*writeridx + kernelidx				= V1_KERNELSIZE*(V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xblockidx) + V1_KERNELWIDTH*unsigned_ry + unsigned_rx
	int kernelidx=threadidx%V1_KERNELSIZE,//kernelidx = V1_KERNELWIDTH*unsigned_ry + unsigned_rx
		ry=kernelidx/(V1_REACH<<1|1),
		rx=kernelidx%(V1_REACH<<1|1),
		writeridx=threadidx/V1_KERNELSIZE,//writeridx = V1_XYWRITERS*kc + V1_XBLOCKS*yblockidx + xblockidx
		kc=writeridx/V1_XYWRITERS,
		spaceidx=writeridx%V1_XYWRITERS,
		yblockidx=spaceidx/V1_XBLOCKS,
		xblockidx=spaceidx%V1_XBLOCKS;
#endif

	int kx0=indices[0], ky0=indices[1], iw=indices[2], ih=indices[3], fwd=indices[4], flip_east_west=indices[5];
	int blockw=(iw+V1_XBLOCKS-1)/V1_XBLOCKS;
	int blockh=(ih+V1_YBLOCKS-1)/V1_YBLOCKS;
	int kx, ky, xstart, xend, ystart, yend;
#ifdef RECYCLE_PARAMS
	int xparam, yparam, recycled_writeridx;
#endif
	__global Params *p, *g;
	__global Temps *t=temps+threadidx;
	__global BwdTemps *b=bwdtemps+threadidx;

	xstart=blockw*xblockidx;
	ystart=blockh*yblockidx;
	xend=xstart+blockw;
	yend=ystart+blockh;
	if(xend>iw)xend=iw;
	if(yend>ih)yend=ih;
	rx-=V1_REACH;
	ry-=V1_REACH;
	kx0+=xstart;
	ky0+=ystart;

#ifdef RECYCLE_PARAMS
	xparam=(xblockidx-ky0)%V1_XBLOCKS, xparam+=V1_XBLOCKS&-(xparam<0);
	recycled_writeridx=V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xparam;
	p=params+recycled_writeridx;
	g=gradients+V1_KERNELSIZE*recycled_writeridx+kernelidx;
#else
	p=params+writeridx;
	g=gradients+threadidx;
#endif
	
	if((unsigned)kx0<(unsigned)iw&&(unsigned)ky0<(unsigned)ih)
	{
		kx=kx0+rx;
		ky=ky0+ry;
		eval_fwd_v1(pixels, errors, iw, ih, kc, kx, ky, xstart, xend, ystart, yend, flip_east_west, p, t);

		eval_bwd_v1(pixels, errors, iw, ih, kc, kx, ky, fwd, p, t, b, g);
	}
}

static float calc_delta(float grad, int sqdist)
{
	//float lr=LR/sqrt((float)sqdist);
	float lr=LR/(float)sqdist;
	//float lr=LR;
	return lr*grad;
}

//worksize: V1_TOTALWRITERS*NPARAMS			params[V1_TOTALWRITERS], gradient[V1_NLANES]
__kernel void update(__global Params *params, __global Params *gradients)
{
	int threadidx=get_global_id(0);//threadidx = NPARAMS*writeridx + paramidx
	int paramidx =threadidx%NPARAMS,
		writeridx=threadidx/NPARAMS;
	//__global float *param=(__global float*)(params+writeridx)+paramidx;
	__global float *param=(__global float*)params+threadidx;
	int kx, ky;

	gradients+=V1_KERNELSIZE*writeridx;
	for(ky=-V1_REACH;ky<0;++ky)
	{
		for(kx=-V1_REACH;kx<=V1_REACH;++kx, ++gradients)
		{
			float gk=((__global float*)gradients)[paramidx];
			*param-=calc_delta(gk, kx*kx+ky*ky);
		}
	}
	for(kx=-V1_REACH;kx<0;++kx, ++gradients)
	{
		float gk=((__global float*)gradients)[paramidx];
		*param-=calc_delta(gk, kx*kx);
	}
}

//worksize: V1_TOTALWRITERS		params[V1_TOTALWRITERS], temps[V1_NLANES] (using only first V1_TOTALWRITERS temps)
__kernel void predict(__global float *pixels, __global float *errors, __constant int *indices, __global Params *params, __global Temps *temps)
{
	int threadidx=get_global_id(0);//threadidx = V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xblockidx
	int xblockidx=threadidx%V1_XBLOCKS,
		rowidx   =threadidx/V1_XBLOCKS,
		yblockidx=rowidx%V1_YBLOCKS,
		kc       =rowidx/V1_YBLOCKS;

	//int threadidx=get_global_id(0);//threadidx = V1_XYWRITERS*kc + spaceidx
	//int kc=threadidx/V1_XYWRITERS,
	//	spaceidx=threadidx%V1_XYWRITERS,//spaceidx = V1_XBLOCKS*yblockidx + xblockidx
	//	yblockidx=spaceidx/V1_XBLOCKS,
	//	xblockidx=spaceidx%V1_XBLOCKS;

	int kx=indices[0], ky=indices[1], iw=indices[2], ih=indices[3], fwd=indices[4], flip_east_west=indices[5];
	int blockw=(iw+V1_XBLOCKS-1)/V1_XBLOCKS;
	int blockh=(ih+V1_YBLOCKS-1)/V1_YBLOCKS;
	int xstart, xend, ystart, yend;
#ifdef RECYCLE_PARAMS
	int xparamidx, recycled_writeridx;
#endif
	__global Params *p;
	__global Temps *t=temps+threadidx;//note: temps is size V1_NLANES > V1_TOTALWRITERS
	int idx=0;

	xstart=blockw*xblockidx;
	ystart=blockh*yblockidx;
	xend=xstart+blockw;
	yend=ystart+blockh;
	if(xend>iw)xend=iw;
	if(yend>ih)yend=ih;

#if 0
	int threadidx=get_global_id(0);//threadidx = V1_XBLOCKS*kc + sblockidx
	int kc=threadidx/V1_XBLOCKS, sblockidx=threadidx%V1_XBLOCKS;
	int kx=indices[0], ky=indices[1], iw=indices[2], ih=indices[3], fwd=indices[4], flip_east_west=indices[5];
	int blockw=(iw+V1_XBLOCKS-1)/V1_XBLOCKS;
	int xstart, xend;
	__global Params *p, *g;
	__global Temps *t;
	__global BwdTemps *b;

	kx+=blockw*sblockidx;
	xstart=blockw*sblockidx;
	xend=xstart+blockw;
	if(xend>iw)
		xend=iw;

	p=params+threadidx;
	t=temps+threadidx;
#endif
	
#ifdef RECYCLE_PARAMS
	xparamidx=xblockidx-ky;
	xparamidx%=V1_XBLOCKS;
	xparamidx+=V1_XBLOCKS&-(xparamidx<0);
	recycled_writeridx=V1_XBLOCKS*(V1_YBLOCKS*kc + yblockidx) + xparamidx;
	p=params+recycled_writeridx;
#else
	p=params+threadidx;
#endif

	//if(xstart>=xend||ystart>=yend)
	//{
	//	if(kx==32&&ky==32)
	//		printf("T%d CXY %d %d %d  [%d, %d], [%d, %d]\n", threadidx, kc, kx, ky, xstart, xend, ystart, yend);
	//}
	if((unsigned)kx<(unsigned)(xend-xstart)&&(unsigned)ky<(unsigned)(yend-ystart))
	{
		kx+=xstart;
		ky+=ystart;
		idx=iw*(ih*kc+ky)+kx;

		eval_fwd_v1(pixels, errors, iw, ih, kc, kx, ky, xstart, xend, ystart, yend, flip_east_west, p, t);

		if(fwd)
			errors[idx]=pixels[idx]-t->pred;
			//errors[idx]=t->pred;
		else
			pixels[idx]=errors[idx]+t->pred;

		//if(!kc&&xblockidx==0&&yblockidx==1&&!kx&&!ky)//
		//if(!kc&&xblockidx==0&&yblockidx==1&&(unsigned)(kx-xstart-100)<100&&(unsigned)(ky-ystart-100)<100)//
		//	printf("T%d: Condition %d CXY %d %d %d, X[%d, %d], Y[%d, %d]  %f - %f = %f\n", threadidx, (unsigned)kx<(unsigned)(xend-xstart)&&(unsigned)ky<(unsigned)(yend-ystart), kc, kx, ky, xstart, xend, ystart, yend, pixels[idx], t->pred, errors[idx]);
	}
	//kx-=xstart;
	//ky-=ystart;
	//if(kx==64&&ky==64)
	//	printf("T%4d CXY %4d %4d %4d (+64)  %10f - %10f = %10f\n", threadidx, kc, xstart, ystart, pixels[idx], t->pred, errors[idx]);
}