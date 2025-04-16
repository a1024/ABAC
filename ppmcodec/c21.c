#include"codec.h"
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;

//	#define LOUD
//	#define ESTIMATE_SIZE
//	#define ENABLE_GUIDE

//	#define USE_RANGECODER

#define NCTX 32
#define NSTREAMS 24

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef ESTIMATE_SIZE
static int g_hist[3][256]={0};
#endif
int c21_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
#endif
#ifdef ESTIMATE_SIZE
	double codecsizes[3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int iw=0, ih=0, fwd=0;
	ptrdiff_t srcsize, res;
	unsigned char *srcbuf=0, *srcptr=0;
//	unsigned char *srcend=0;
	unsigned char *decbuf=0;
	unsigned char *image=0;
	int headersize=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<0)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		if(srcsize<=2)
		{
			LOG_ERROR("File is empty");
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		unsigned char header[128];
		fread(header, 1, sizeof(header)-1, fsrc);
		srcptr=header;
		fwd=!memcmp(srcptr, "P6\n", 3);
		if(!fwd&&memcmp(srcptr, "CH", 2))
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		if(fwd)
		{
			srcptr+=3;
			while((unsigned)(*srcptr-'0')<10)
				iw=10*iw+*srcptr++ - '0';
			if(iw<1||*srcptr++ != ' ')
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			while((unsigned)(*srcptr-'0')<10)
				ih=10*ih+*srcptr++ - '0';
			if(memcmp(srcptr, "\n255\n", 5))
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			srcptr+=5;
			if(iw<1||ih<1)
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			res=(ptrdiff_t)iw*ih;

			headersize=(int)(srcptr-header);
			fseek(fsrc, 0, SEEK_SET);
			srcbuf=(unsigned char*)malloc(headersize+sizeof(char[3])*iw*(ih+3LL));
			if(!srcbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			image=srcbuf+headersize+sizeof(char[3])*iw*3;
			memset(srcbuf, 0, image-srcbuf);
			fread(image-headersize, 1, headersize+3LL*res, fsrc);
			srcptr=image;
		//	srcend=image+3LL*res;
		}
		else
		{
			srcptr+=2;
			memcpy(&iw, srcptr, 4); srcptr+=4;
			memcpy(&ih, srcptr, 4); srcptr+=4;
			if(iw<1||ih<1)
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			res=(ptrdiff_t)iw*ih;

			fseek(fsrc, 0, SEEK_SET);
			srcbuf=(unsigned char*)malloc(srcsize+16);
			if(!srcbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			fread(srcbuf, 1, srcsize, fsrc);
			srcbuf[srcsize]=0;
			srcptr=srcbuf+10;//CH header size = 10
		//	srcend=srcbuf+srcsize;
			
			headersize=snprintf((char*)header, sizeof(header)-1, "P6\n%d %d\n255\n", iw, ih);
			decbuf=(unsigned char*)malloc(headersize+sizeof(char[3])*iw*(ih+3LL));
			if(!decbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			memset(decbuf, 0, headersize+iw*9LL);
			image=decbuf+headersize+iw*9LL;
			memcpy(image-headersize, header, headersize);
		}
		fclose(fsrc);
	}
	unsigned char *NNNptr	=image-3*3*iw;
	unsigned char *NNptr	=image-2*3*iw;
	unsigned char *Nptr	=image-1*3*iw;
	unsigned char *currptr	=image+0*3*iw;

	int energysize=(int)sizeof(char[3])*(iw+8);
	unsigned char *energy=(unsigned char*)malloc(energysize);
	int hsize=(int)sizeof(int[3*NCTX*257]);
	unsigned short *hist=(unsigned short*)malloc(hsize);
	unsigned short *CDF=(unsigned short*)malloc(hsize);
	if(!energy||!hist||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(energy, 0, energysize);
	memset(hist, 0, hsize);
	for(int k=0;k<3*NCTX;++k)
	{
		for(int ks=0;ks<257;++ks)
			CDF[257*k+ks]=ks<<(12-8);
	}
	//unsigned short hist[3][257]={0}, CDF[3][257]={0};
	//for(int ks=0;ks<257;++ks)
	//{
	//	int c=ks<<(12-8);
	//	CDF[0][ks]=c;
	//	CDF[1][ks]=c;
	//	CDF[2][ks]=c;
	//}
#ifdef ESTIMATE_SIZE
	memset(g_hist, 0, sizeof(g_hist));
#endif

	ALIGN(32) unsigned low[NSTREAMS]={0}, range[NSTREAMS]={0}, code[NSTREAMS]={0}, code2[NSTREAMS]={0};
	memset(range, -1, sizeof(range));
	//for(int k=0;k<NSTREAMS;++k)
	//	range[k]=0xFFFFFFFF;
	//__m256i mlow0=_mm256_load_si256((__m256i*)low+0);
	//__m256i mlow1=_mm256_load_si256((__m256i*)low+1);
	//__m256i mlow2=_mm256_load_si256((__m256i*)low+2);
	//__m256i mrange0=_mm256_load_si256((__m256i*)range+0);
	//__m256i mrange1=_mm256_load_si256((__m256i*)range+1);
	//__m256i mrange2=_mm256_load_si256((__m256i*)range+2);

	ALIGN(32) int weights[3*8]={0}, errors[3*8]={0}, circlebuf[3*8]={0};
	FILLMEM(weights, 0x10000/8, sizeof(weights), sizeof(int));

	ptrdiff_t maxemitts=res/12;
	unsigned short *dstbufs[NSTREAMS]={0}, *streamptrs[NSTREAMS]={0};
	int streamsizes[NSTREAMS]={0};
	if(fwd)
	{
		guide_save(image, iw, ih);
		for(int k=0;k<NSTREAMS;++k)
		{
			dstbufs[k]=(unsigned short*)malloc(sizeof(short)*maxemitts);
			if(!dstbufs[k])
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			streamptrs[k]=dstbufs[k];
		}
	}
	else
	{
		memcpy(streamsizes, srcptr, sizeof(streamsizes)); srcptr+=sizeof(streamsizes);
		for(int k=0;k<NSTREAMS;++k)
		{
			streamptrs[k]=(unsigned short*)srcptr;
			srcptr+=streamsizes[k]*sizeof(short);
		}
		for(int k=0;k<NSTREAMS;++k)
		{
			code[k]=*streamptrs[k]++;
			code[k]=code[k]<<16|*streamptrs[k]++;
		}
	}
	int yctx=0, uctx=0, vctx=0;
	unsigned short *CDFs[]=
	{
		CDF+0*NCTX*257,
		CDF+1*NCTX*257,
		CDF+2*NCTX*257,
	};
	unsigned short *hists[]=
	{
		hist+0*NCTX*257,
		hist+1*NCTX*257,
		hist+2*NCTX*257,
	};
	for(int kp=0, kx=0;kp<res;++kp)
	{
		int *curr_weights=weights, *curr_errors=errors;
		int update=fwd?((kp&7)==7||kp==res-1):!(kp&7);
		
		//if(kp==127048)//
		//if(kp==17544)//
		//if(kp==16801)//
		//if(kp==14602)//
		//if(kp==13920)//
		//if(kp==15768)//
		//if(kp==7776)//
		//if(kp==2027)//
		//	printf("");

		int
			yN	=Nptr	[1+0*3],
			uN	=Nptr	[2+0*3]-yN,
			vN	=Nptr	[0+0*3]-yN,
			yNE	=Nptr	[1+1*3],
			uNE	=Nptr	[2+1*3]-yNE,
			vNE	=Nptr	[0+1*3]-yNE,
			yW	=currptr[1-1*3],
			uW	=currptr[2-1*3]-yW,
			vW	=currptr[0-1*3]-yW;
		int ymax=yN, ymin=yW;
		int umax=uN, umin=uW;
		int vmax=vN, vmin=vW;
		if(yN<yW)ymin=yN, ymax=yW;
		if(uN<uW)umin=uN, umax=uW;
		if(vN<vW)vmin=vN, vmax=vW;
		if(ymin>yNE)ymin=yNE;
		if(ymax<yNE)ymax=yNE;
		if(umin>uNE)umin=uNE;
		if(umax<uNE)umax=uNE;
		if(vmin>vNE)vmin=vNE;
		if(vmax<vNE)vmax=vNE;
		int preds[]=
		{
			Nptr[1+0*3],//N
			Nptr[2+0*3],
			Nptr[0+0*3],
			currptr[1-1*3],//W
			currptr[2-1*3],
			currptr[0-1*3],
			3*(Nptr[1+0*3]-NNptr[1+0*3])+NNNptr[1+0*3],//3*(N-NN)+NNN
			3*(Nptr[2+0*3]-NNptr[2+0*3])+NNNptr[2+0*3],
			3*(Nptr[0+0*3]-NNptr[0+0*3])+NNNptr[0+0*3],
			3*(currptr[1-1*3]-currptr[1-2*3])+currptr[1-3*3],//3*(W-WW)+WWW
			3*(currptr[2-1*3]-currptr[2-2*3])+currptr[2-3*3],
			3*(currptr[0-1*3]-currptr[0-2*3])+currptr[0-3*3],
			currptr[1-1*3]+Nptr[1+1*3]-Nptr[1+0*3],//W+NE-N
			currptr[2-1*3]+Nptr[2+1*3]-Nptr[2+0*3],
			currptr[0-1*3]+Nptr[0+1*3]-Nptr[0+0*3],
			(currptr[1-4*3]+currptr[1-3*3]+NNNptr[1+0*3]+Nptr[1+2*3]+Nptr[1+3*3]+Nptr[1+4*3]-2*Nptr[1-1*3]+2)>>2,//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW+2)>>2
			(currptr[2-4*3]+currptr[2-3*3]+NNNptr[2+0*3]+Nptr[2+2*3]+Nptr[2+3*3]+Nptr[2+4*3]-2*Nptr[2-1*3]+2)>>2,
			(currptr[0-4*3]+currptr[0-3*3]+NNNptr[0+0*3]+Nptr[0+2*3]+Nptr[0+3*3]+Nptr[0+4*3]-2*Nptr[0-1*3]+2)>>2,
			Nptr[1+0*3]+currptr[1-1*3]-Nptr[1-1*3],//N+W-NW
			Nptr[2+0*3]+currptr[2-1*3]-Nptr[2-1*3],
			Nptr[0+0*3]+currptr[0-1*3]-Nptr[0-1*3],
			Nptr[1+0*3]+Nptr[1+1*3]-NNptr[1+1*3],//N+NE-NNE
			Nptr[2+0*3]+Nptr[2+1*3]-NNptr[2+1*3],
			Nptr[0+0*3]+Nptr[0+1*3]-NNptr[0+1*3],
		};
		preds[1+0*3]-=preds[0+0*3];
		preds[2+0*3]-=preds[0+0*3];
		preds[1+1*3]-=preds[0+1*3];
		preds[2+1*3]-=preds[0+1*3];
		preds[1+2*3]-=preds[0+2*3];
		preds[2+2*3]-=preds[0+2*3];
		preds[1+3*3]-=preds[0+3*3];
		preds[2+3*3]-=preds[0+3*3];
		preds[1+4*3]-=preds[0+4*3];
		preds[2+4*3]-=preds[0+4*3];
		preds[1+5*3]-=preds[0+5*3];
		preds[2+5*3]-=preds[0+5*3];
		preds[1+6*3]-=preds[0+6*3];
		preds[2+6*3]-=preds[0+6*3];
		preds[1+7*3]-=preds[0+7*3];
		preds[2+7*3]-=preds[0+7*3];
		int ypred=(
			+curr_weights[0+0*3]*preds[0+0*3]
			+curr_weights[0+1*3]*preds[0+1*3]
			+curr_weights[0+2*3]*preds[0+2*3]
			+curr_weights[0+3*3]*preds[0+3*3]
			+curr_weights[0+4*3]*preds[0+4*3]
			+curr_weights[0+5*3]*preds[0+5*3]
			+curr_weights[0+6*3]*preds[0+6*3]
			+curr_weights[0+7*3]*preds[0+7*3]
			+0x8000
		)>>16;
		int upred=(
			+curr_weights[1+0*3]*preds[1+0*3]
			+curr_weights[1+1*3]*preds[1+1*3]
			+curr_weights[1+2*3]*preds[1+2*3]
			+curr_weights[1+3*3]*preds[1+3*3]
			+curr_weights[1+4*3]*preds[1+4*3]
			+curr_weights[1+5*3]*preds[1+5*3]
			+curr_weights[1+6*3]*preds[1+6*3]
			+curr_weights[1+7*3]*preds[1+7*3]
			+0x8000
		)>>16;
		int vpred=(
			+curr_weights[2+0*3]*preds[2+0*3]
			+curr_weights[2+1*3]*preds[2+1*3]
			+curr_weights[2+2*3]*preds[2+2*3]
			+curr_weights[2+3*3]*preds[2+3*3]
			+curr_weights[2+4*3]*preds[2+4*3]
			+curr_weights[2+5*3]*preds[2+5*3]
			+curr_weights[2+6*3]*preds[2+6*3]
			+curr_weights[2+7*3]*preds[2+7*3]
			+0x8000
		)>>16;
		CLAMP2(ypred, ymin, ymax);
		CLAMP2(upred, umin, umax);
		CLAMP2(vpred, vmin, vmax);
		//yctx+=Nptr[0+8*3];
		//uctx+=Nptr[1+8*3]-Nptr[0+8*3];
		//vctx+=Nptr[2+8*3]-Nptr[0+8*3];
		//pred+=offset;
		//CLAMP2(pred, 0, 255);
		
		if(!(kp&7))
		{
			yctx=energy[kx+0+0*3]+energy[kx+0+1*3]+energy[kx+0+2*3]+energy[kx+0+3*3]+energy[kx+0+4*3]+energy[kx+0+5*3]+energy[kx+0+6*3]+energy[kx+0+7*3];
			uctx=energy[kx+1+0*3]+energy[kx+1+1*3]+energy[kx+1+2*3]+energy[kx+1+3*3]+energy[kx+1+4*3]+energy[kx+1+5*3]+energy[kx+1+6*3]+energy[kx+1+7*3];
			vctx=energy[kx+2+0*3]+energy[kx+2+1*3]+energy[kx+2+2*3]+energy[kx+2+3*3]+energy[kx+2+4*3]+energy[kx+2+5*3]+energy[kx+2+6*3]+energy[kx+2+7*3];
			//yctx*=yctx;
			//uctx*=uctx;
			//vctx*=vctx;
			yctx=FLOOR_LOG2_P1(yctx);
			uctx=FLOOR_LOG2_P1(uctx);
			vctx=FLOOR_LOG2_P1(vctx);
			yctx&=255;
			uctx&=255;
			vctx&=255;
			CDFs[0]=CDF+257*yctx+0*NCTX*257;
			CDFs[1]=CDF+257*uctx+1*NCTX*257;
			CDFs[2]=CDF+257*vctx+2*NCTX*257;
			hists[0]=hist+257*yctx+0*NCTX*257;
			hists[1]=hist+257*uctx+1*NCTX*257;
			hists[2]=hist+257*vctx+2*NCTX*257;
		}
		if(update)
		{
			if(fwd)
			{
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					if(range[k]<0x1000)//enc renorm
					{
						*streamptrs[k]++=low[k]>>16;
						low[k]<<=16;
						range[k]=range[k]<<16|0xFFFF;
						unsigned rmax=~low[k];
						if(range[k]>rmax)
							range[k]=rmax;
					}
				}
			}
			else
			{
				//if(kp==5760)//
				//if(kp==5752)//
				//	printf("");//
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					if(range[k]<0x1000)//dec renorm
					{
						code[k]=code[k]<<16|*streamptrs[k]++;
						low[k]<<=16;
						range[k]=range[k]<<16|0xFFFF;
						unsigned rmax=~low[k];
						if(range[k]>rmax)
							range[k]=rmax;
					}
				}
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
				for(int k=0;k<NSTREAMS;++k)
#ifdef USE_RANGECODER
				//	code2[k]=code[k]-low[k];
					code2[k]=(code[k]-low[k])/(range[k]>>12);
#else
					code2[k]=(unsigned)(((unsigned long long)(code[k]-low[k])<<12|0xFFF)/range[k]);
#endif
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					int kc=k>>3;
					int sym;

					sym=0;
					for(;;)//dec search
					{
					//	if((range[k]>>12)*CDFs[kc][sym+2]>code2[k])
						if(CDFs[kc][sym+2]>code2[k])
						{
						//	sym+=(range[k]>>12)*CDFs[kc][sym+1]<=code2[k];
							sym+=CDFs[kc][sym+1]<=code2[k];
							break;
						}
						sym+=2;
#ifdef _DEBUG
						if(sym>255)
							LOG_ERROR("Decode error at %d", kp);
#endif
					}
					circlebuf[k]=sym;
				}
			}
		}
		
		int curr[3], uerror[3], serror[3];
		//if(kp==5762)//
		//	printf("");
		if(fwd)
		{
			curr[0]=currptr[1];
			curr[1]=currptr[2];
			curr[2]=currptr[0];

			upred+=curr[0];
			CLAMP2(upred, 0, 255);
			vpred+=curr[0];
			CLAMP2(vpred, 0, 255);

			serror[0]=(char)(curr[0]-ypred);
			serror[1]=(char)(curr[1]-upred);
			serror[2]=(char)(curr[2]-vpred);
			uerror[0]=serror[0]<<1^serror[0]>>31;//pack sign
			uerror[1]=serror[1]<<1^serror[1]>>31;
			uerror[2]=serror[2]<<1^serror[2]>>31;
			circlebuf[0<<3|(kp&7)]=uerror[0];
			circlebuf[1<<3|(kp&7)]=uerror[1];
			circlebuf[2<<3|(kp&7)]=uerror[2];
		}
		else
		{
			uerror[0]=circlebuf[0<<3|(kp&7)];
			serror[0]=uerror[0]>>1^-(uerror[0]&1);//unpack sign
			currptr[1]=curr[0]=(unsigned char)(serror[0]+ypred);

			upred+=curr[0];
			CLAMP2(upred, 0, 255);
			vpred+=curr[0];
			CLAMP2(vpred, 0, 255);

			uerror[1]=circlebuf[1<<3|(kp&7)];
			uerror[2]=circlebuf[2<<3|(kp&7)];
			serror[1]=uerror[1]>>1^-(uerror[1]&1);
			serror[2]=uerror[2]>>1^-(uerror[2]&1);
			currptr[2]=curr[1]=(unsigned char)(serror[1]+upred);
			currptr[0]=curr[2]=(unsigned char)(serror[2]+vpred);
		}
		energy[kx+0]=abs(serror[0]);
		energy[kx+1]=abs(serror[1]);
		energy[kx+2]=abs(serror[2]);
		//energy[kx+0]=uerror[0];
		//energy[kx+1]=uerror[1];
		//energy[kx+2]=uerror[2];
		++hists[0][uerror[0]];
		++hists[1][uerror[1]];
		++hists[2][uerror[2]];
		++hists[0][256];
		++hists[1][256];
		++hists[2][256];
#ifdef ESTIMATE_SIZE
		++g_hist[0][uerror[0]];
		++g_hist[1][uerror[1]];
		++g_hist[2][uerror[2]];
#endif

		curr_errors[ 0]+=abs(preds[ 0]-curr[0]);
		curr_errors[ 1]+=abs(preds[ 1]-curr[1]);
		curr_errors[ 2]+=abs(preds[ 2]-curr[2]);
		curr_errors[ 3]+=abs(preds[ 3]-curr[0]);
		curr_errors[ 4]+=abs(preds[ 4]-curr[1]);
		curr_errors[ 5]+=abs(preds[ 5]-curr[2]);
		curr_errors[ 6]+=abs(preds[ 6]-curr[0]);
		curr_errors[ 7]+=abs(preds[ 7]-curr[1]);
		curr_errors[ 8]+=abs(preds[ 8]-curr[2]);
		curr_errors[ 9]+=abs(preds[ 9]-curr[0]);
		curr_errors[10]+=abs(preds[10]-curr[1]);
		curr_errors[11]+=abs(preds[11]-curr[2]);
		curr_errors[12]+=abs(preds[12]-curr[0]);
		curr_errors[13]+=abs(preds[13]-curr[1]);
		curr_errors[14]+=abs(preds[14]-curr[2]);
		curr_errors[15]+=abs(preds[15]-curr[0]);
		curr_errors[16]+=abs(preds[16]-curr[1]);
		curr_errors[17]+=abs(preds[17]-curr[2]);
		curr_errors[18]+=abs(preds[18]-curr[0]);
		curr_errors[19]+=abs(preds[19]-curr[1]);
		curr_errors[20]+=abs(preds[20]-curr[2]);
		curr_errors[21]+=abs(preds[21]-curr[0]);
		curr_errors[22]+=abs(preds[22]-curr[1]);
		curr_errors[23]+=abs(preds[23]-curr[2]);
		if(!(kp&127))
		{
			//FIXME use _mm_div_ps
			unsigned ww[]=
			{
				0x8000/(curr_errors[ 0]+1),//y		//FIXME try curr_weights[0]/(curr_errors[0]+1)
				0x8000/(curr_errors[ 1]+1),//u
				0x8000/(curr_errors[ 2]+1),//v
				0x8000/(curr_errors[ 3]+1),
				0x8000/(curr_errors[ 4]+1),
				0x8000/(curr_errors[ 5]+1),
				0x8000/(curr_errors[ 6]+1),
				0x8000/(curr_errors[ 7]+1),
				0x8000/(curr_errors[ 8]+1),
				0x8000/(curr_errors[ 9]+1),
				0x8000/(curr_errors[10]+1),
				0x8000/(curr_errors[11]+1),
				0x8000/(curr_errors[12]+1),
				0x8000/(curr_errors[13]+1),
				0x8000/(curr_errors[14]+1),
				0x8000/(curr_errors[15]+1),
				0x8000/(curr_errors[16]+1),
				0x8000/(curr_errors[17]+1),
				0x8000/(curr_errors[18]+1),
				0x8000/(curr_errors[19]+1),
				0x8000/(curr_errors[20]+1),
				0x8000/(curr_errors[21]+1),
				0x8000/(curr_errors[22]+1),
				0x8000/(curr_errors[23]+1),
			};
			unsigned ysum=ww[0+0*3]+ww[0+1*3]+ww[0+2*3]+ww[0+3*3]+ww[0+4*3]+ww[0+5*3]+ww[0+6*3]+ww[0+7*3]+1;
			unsigned usum=ww[1+0*3]+ww[1+1*3]+ww[1+2*3]+ww[1+3*3]+ww[1+4*3]+ww[1+5*3]+ww[1+6*3]+ww[1+7*3]+1;
			unsigned vsum=ww[2+0*3]+ww[2+1*3]+ww[2+2*3]+ww[2+3*3]+ww[2+4*3]+ww[2+5*3]+ww[2+6*3]+ww[2+7*3]+1;
			curr_weights[ 0]=(ww[ 0]<<16)/ysum;
			curr_weights[ 1]=(ww[ 1]<<16)/usum;
			curr_weights[ 2]=(ww[ 2]<<16)/vsum;
			curr_weights[ 3]=(ww[ 3]<<16)/ysum;
			curr_weights[ 4]=(ww[ 4]<<16)/usum;
			curr_weights[ 5]=(ww[ 5]<<16)/vsum;
			curr_weights[ 6]=(ww[ 6]<<16)/ysum;
			curr_weights[ 7]=(ww[ 7]<<16)/usum;
			curr_weights[ 8]=(ww[ 8]<<16)/vsum;
			curr_weights[ 9]=(ww[ 9]<<16)/ysum;
			curr_weights[10]=(ww[10]<<16)/usum;
			curr_weights[11]=(ww[11]<<16)/vsum;
			curr_weights[12]=(ww[12]<<16)/ysum;
			curr_weights[13]=(ww[13]<<16)/usum;
			curr_weights[14]=(ww[14]<<16)/vsum;
			curr_weights[15]=(ww[15]<<16)/ysum;
			curr_weights[16]=(ww[16]<<16)/usum;
			curr_weights[17]=(ww[17]<<16)/vsum;
			curr_weights[18]=(ww[18]<<16)/ysum;
			curr_weights[19]=(ww[19]<<16)/usum;
			curr_weights[20]=(ww[20]<<16)/vsum;
			curr_weights[21]=0x10000-(
				+curr_weights[0+0*3]
				+curr_weights[0+1*3]
				+curr_weights[0+2*3]
				+curr_weights[0+3*3]
				+curr_weights[0+4*3]
				+curr_weights[0+5*3]
				+curr_weights[0+6*3]
			);
			curr_weights[22]=0x10000-(
				+curr_weights[1+0*3]
				+curr_weights[1+1*3]
				+curr_weights[1+2*3]
				+curr_weights[1+3*3]
				+curr_weights[1+4*3]
				+curr_weights[1+5*3]
				+curr_weights[1+6*3]
			);
			curr_weights[23]=0x10000-(
				+curr_weights[2+0*3]
				+curr_weights[2+1*3]
				+curr_weights[2+2*3]
				+curr_weights[2+3*3]
				+curr_weights[2+4*3]
				+curr_weights[2+5*3]
				+curr_weights[2+6*3]
			);
#ifdef _DEBUG
			if(
				curr_weights[ 0]<0||
				curr_weights[ 1]<0||
				curr_weights[ 2]<0||
				curr_weights[ 3]<0||
				curr_weights[ 4]<0||
				curr_weights[ 5]<0||
				curr_weights[ 6]<0||
				curr_weights[ 7]<0||
				curr_weights[ 8]<0||
				curr_weights[ 9]<0||
				curr_weights[10]<0||
				curr_weights[11]<0||
				curr_weights[12]<0||
				curr_weights[13]<0||
				curr_weights[14]<0||
				curr_weights[15]<0||
				curr_weights[16]<0||
				curr_weights[17]<0||
				curr_weights[18]<0||
				curr_weights[19]<0||
				curr_weights[20]<0||
				curr_weights[21]<0||
				curr_weights[22]<0||
				curr_weights[23]<0||
				curr_weights[0+0*3]+curr_weights[0+1*3]+curr_weights[0+2*3]+curr_weights[0+3*3]+curr_weights[0+4*3]+curr_weights[0+5*3]+curr_weights[0+6*3]+curr_weights[0+7*3]!=0x10000||
				curr_weights[1+0*3]+curr_weights[1+1*3]+curr_weights[1+2*3]+curr_weights[1+3*3]+curr_weights[1+4*3]+curr_weights[1+5*3]+curr_weights[1+6*3]+curr_weights[1+7*3]!=0x10000||
				curr_weights[2+0*3]+curr_weights[2+1*3]+curr_weights[2+2*3]+curr_weights[2+3*3]+curr_weights[2+4*3]+curr_weights[2+5*3]+curr_weights[2+6*3]+curr_weights[2+7*3]!=0x10000
			)
				LOG_ERROR("");
#endif
			memset(curr_errors, 0, sizeof(int[3*8]));//FIXME try halving errors
		}
		if(!fwd)
			guide_check(image, kp%iw, kp/iw);//

		if(update)
		{
			//if(kp==5759)//
			//	printf("");

			int cdf[3*8]=
			{
				CDFs[0][circlebuf[0+0*8]],
				CDFs[0][circlebuf[1+0*8]],
				CDFs[0][circlebuf[2+0*8]],
				CDFs[0][circlebuf[3+0*8]],
				CDFs[0][circlebuf[4+0*8]],
				CDFs[0][circlebuf[5+0*8]],
				CDFs[0][circlebuf[6+0*8]],
				CDFs[0][circlebuf[7+0*8]],
				CDFs[1][circlebuf[0+1*8]],
				CDFs[1][circlebuf[1+1*8]],
				CDFs[1][circlebuf[2+1*8]],
				CDFs[1][circlebuf[3+1*8]],
				CDFs[1][circlebuf[4+1*8]],
				CDFs[1][circlebuf[5+1*8]],
				CDFs[1][circlebuf[6+1*8]],
				CDFs[1][circlebuf[7+1*8]],
				CDFs[2][circlebuf[0+2*8]],
				CDFs[2][circlebuf[1+2*8]],
				CDFs[2][circlebuf[2+2*8]],
				CDFs[2][circlebuf[3+2*8]],
				CDFs[2][circlebuf[4+2*8]],
				CDFs[2][circlebuf[5+2*8]],
				CDFs[2][circlebuf[6+2*8]],
				CDFs[2][circlebuf[7+2*8]],
			};
			int freq[3*8]=
			{
				CDFs[0][circlebuf[0+0*8]+1]-cdf[0+0*8],
				CDFs[0][circlebuf[1+0*8]+1]-cdf[1+0*8],
				CDFs[0][circlebuf[2+0*8]+1]-cdf[2+0*8],
				CDFs[0][circlebuf[3+0*8]+1]-cdf[3+0*8],
				CDFs[0][circlebuf[4+0*8]+1]-cdf[4+0*8],
				CDFs[0][circlebuf[5+0*8]+1]-cdf[5+0*8],
				CDFs[0][circlebuf[6+0*8]+1]-cdf[6+0*8],
				CDFs[0][circlebuf[7+0*8]+1]-cdf[7+0*8],
				CDFs[1][circlebuf[0+1*8]+1]-cdf[0+1*8],
				CDFs[1][circlebuf[1+1*8]+1]-cdf[1+1*8],
				CDFs[1][circlebuf[2+1*8]+1]-cdf[2+1*8],
				CDFs[1][circlebuf[3+1*8]+1]-cdf[3+1*8],
				CDFs[1][circlebuf[4+1*8]+1]-cdf[4+1*8],
				CDFs[1][circlebuf[5+1*8]+1]-cdf[5+1*8],
				CDFs[1][circlebuf[6+1*8]+1]-cdf[6+1*8],
				CDFs[1][circlebuf[7+1*8]+1]-cdf[7+1*8],
				CDFs[2][circlebuf[0+2*8]+1]-cdf[0+2*8],
				CDFs[2][circlebuf[1+2*8]+1]-cdf[1+2*8],
				CDFs[2][circlebuf[2+2*8]+1]-cdf[2+2*8],
				CDFs[2][circlebuf[3+2*8]+1]-cdf[3+2*8],
				CDFs[2][circlebuf[4+2*8]+1]-cdf[4+2*8],
				CDFs[2][circlebuf[5+2*8]+1]-cdf[5+2*8],
				CDFs[2][circlebuf[6+2*8]+1]-cdf[6+2*8],
				CDFs[2][circlebuf[7+2*8]+1]-cdf[7+2*8],
			};
#ifdef USE_RANGECODER
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
			for(int k=0;k<NSTREAMS;++k)
				low[k]+=(range[k]>>12)*cdf[k];
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
			for(int k=0;k<NSTREAMS;++k)
				range[k]=(range[k]>>12)*freq[k]-1;
#else
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
			for(int k=0;k<NSTREAMS;++k)
				low[k]+=(unsigned)((unsigned long long)range[k]*cdf[k]>>12);
#ifdef __GNUC__
#pragma GCC unroll 24
#endif
			for(int k=0;k<NSTREAMS;++k)
				range[k]=(unsigned)((unsigned long long)range[k]*freq[k]>>12)-1;
#endif
#ifdef ESTIMATE_SIZE
			for(int k=0;k<24;++k)
				codecsizes[k>>3]-=log2(freq[k]*(1./0x1000));
#endif
		}
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
		for(int kc=0;kc<3;++kc)
		{
			unsigned short *curr_hist=hists[kc];
			if(curr_hist[256]==4096-256)//snapshot-CDF
			{
				//if(kp==5765&&kc==1)//
				//	printf("");
				unsigned short *curr_CDF=CDFs[kc];
				//int hsum=0;
				for(int ks=0, sum=0;ks<256;++ks)
				{
					curr_CDF[ks]=sum+ks;
					sum+=curr_hist[ks];
					curr_hist[ks]=0;
					//hsum+=curr_hist[ks]>>=1;	//X  enc shoots early, dec shoots late
				}
				curr_hist[256]=0;
				//curr_hist[256]=hsum;
			}
		}
		kx+=3;
		if(kx>=iw*3)
			kx=0;
		NNptr	+=3;
		Nptr	+=3;
		currptr	+=3;
	}
	free(energy);
	free(hist);
	free(CDF);
#ifdef LOUD
	ptrdiff_t usize=res*3;
#endif
	ptrdiff_t csize_actual=0;
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		for(int k=0;k<NSTREAMS;++k)//flush
		{
			*streamptrs[k]++=low[k]>>16;
			*streamptrs[k]++=low[k];
		}
		csize_actual+=fwrite("CH", 1, 2, fdst);
		csize_actual+=fwrite(&iw, 1, 4, fdst);
		csize_actual+=fwrite(&ih, 1, 4, fdst);
		for(int k=0;k<NSTREAMS;++k)
			streamsizes[k]=(int)(streamptrs[k]-dstbufs[k]);
		csize_actual+=fwrite(streamsizes, 1, sizeof(int[NSTREAMS]), fdst);
		for(int k=0;k<NSTREAMS;++k)
			csize_actual+=fwrite(dstbufs[k], 1, streamsizes[k]*sizeof(short), fdst);

		for(int k=0;k<NSTREAMS;++k)
			free(dstbufs[k]);
	}
	else
	{
		fwrite(image-headersize, 1, headersize+3*res, fdst);
		free(decbuf);
	}
	fclose(fdst);
	free(srcbuf);
	
#ifdef ESTIMATE_SIZE
	codecsizes[0]/=8;
	codecsizes[1]/=8;
	codecsizes[2]/=8;
	if(fwd)
	{
		double o0sizes[3]={0};
		for(int kc=0;kc<3;++kc)
		{
			int *curr_hist=g_hist[kc];
			double e=0;
			int hsum=0;
			for(int ks=0;ks<256;++ks)
				hsum+=curr_hist[ks];
			double gain=1./hsum;
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					e-=freq*log2((double)freq*gain);
			}
			o0sizes[kc]=e/8;
		}
		printf("T %12.2lf -> %12.2lf\n",
			o0sizes[0]+o0sizes[1]+o0sizes[2],
			codecsizes[0]+codecsizes[1]+codecsizes[2]
		);
		printf("Y %12.2lf -> %12.2lf\n", o0sizes[0], codecsizes[0]);
		printf("U %12.2lf -> %12.2lf\n", o0sizes[1], codecsizes[1]);
		printf("V %12.2lf -> %12.2lf\n", o0sizes[2], codecsizes[2]);
	}
#endif
#ifdef LOUD
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	if(fwd)
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	printf("%c%12lf sec %12lf MB/s  %12lld cycles %12lf C/B\n",
		'D'+fwd,
		elapsed,
		usize/(elapsed*1024*1024),
		cycles,
		(double)cycles/usize
	);
#endif
	return 0;
}