#include"ebench.h"
#include<stdint.h>
#include<stdio.h>//snprintf
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;

typedef struct SampleCtxStruct
{
	ArrayHandle fn;
	double usize, csize[4], bpd0;
	ptrdiff_t idx, bpdrank;
} SampleCtx;
static unsigned __stdcall analysis(void *param)
{
	enum
	{
		NCTX=18,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,
	};
	int nch=0, depth=0;
	int hsize=0;
	int32_t *hists=0;
	int psize=0;
	int16_t *pixels=0;
	SampleCtx *tctx=(SampleCtx*)param;
	Image *src=image_load((char*)tctx->fn->data, (int)tctx->fn->count);
	int bestrct=0, yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;

	if(!src)
		return 0;
	int amin[]=
	{
		-(1<<src->depth[0]>>1),
		-(1<<src->depth[1]>>1),
		-(1<<src->depth[2]>>1),
		-(1<<src->depth[3]>>1),
	};
	int amax[]=
	{
		(1<<src->depth[0]>>1)-1,
		(1<<src->depth[1]>>1)-1,
		(1<<src->depth[2]>>1)-1,
		(1<<src->depth[3]>>1)-1,
	};

	depth=src->depth[0];
	if(depth<src->depth[1])depth=src->depth[1];
	if(depth<src->depth[2])depth=src->depth[2];
	if(depth<src->depth[3])depth=src->depth[3];
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	hsize=(int)sizeof(int32_t[NCTX])*nch<<depth;
	hists=(int32_t*)malloc(hsize);
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!hists||!pixels)
	{
		LOG_ERROR("ALloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	memset(hists, 0, hsize);
	bestrct=crct_analysis(src);
	yidx=rct_combinations[bestrct][II_PERM_Y];
	uidx=rct_combinations[bestrct][II_PERM_U];
	vidx=rct_combinations[bestrct][II_PERM_V];
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0;
			int yuv[]=
			{
				src->data[idx+yidx],
				src->data[idx+uidx],
				src->data[idx+vidx],
				src->data[idx+3],
			};
			for(int kc=0;kc<4;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				
				int16_t
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					eNEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int ctx=FLOOR_LOG2(eW*eW+1);
				int pred=N+W-NW, vmax=N, vmin=W;
				int error, sym;

				if(ctx>NCTX-1)
					ctx=NCTX-1;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, amin[kc], amax[kc]);

				//fwd
				error=yuv[kc]-pred;
				error<<=32-src->depth[kc];
				error>>=32-src->depth[kc];
				sym=error<<1^error>>31;
				++hists[(NCTX*kc+ctx)<<depth|sym];

				rows[0][0]=yuv[kc]-offset;
				rows[0][1]=(2*eW+8*sym+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				if(!kc)
					offset=uc0*yuv[0]>>2;
				else if(kc==1)
					offset=(vc0*yuv[0]+vc1*yuv[1])>>2;
				else
					offset=0;
			}
		}
	}
	free(pixels);

	int nlevels=1<<depth;
	double csize=0;
	for(int kc=0;kc<NCTX*nch;++kc)
	{
		int32_t *currhist=hists+((ptrdiff_t)kc<<depth), sum=0;
		for(int ks=0;ks<nlevels;++ks)
			sum+=currhist[ks];
		if(!sum)
			continue;
		double invsum=1./sum;
		for(int ks=0;ks<nlevels;++ks)
		{
			int freq=currhist[ks];
			if(freq)
				csize-=freq*log2(freq*invsum);
		}
	}
	free(hists);

	//BPD = (csize / usize) * depth
	tctx->usize=(double)src->iw*src->ih*nch*depth;
	tctx->bpd0=csize/((double)src->iw*src->ih*nch);
	free(src);
	return 0;
}
static unsigned __stdcall sample_thread(void *param)
{
	double entropy[4]={0};
	SampleCtx *tctx=(SampleCtx*)param;
	Image *src=image_load((char*)tctx->fn->data, (int)tctx->fn->count);
	if(!src)
		return 0;

	tctx->usize=image_getBMPsize(src);
	apply_selected_transforms(&src, 0, 1, 1);
	calc_csize_stateful(src, 0, entropy);
	for(int kc=0;kc<4;++kc)
	{
		int depth=src->src_depth[kc];
		double invCR=depth?entropy[kc]/depth:0;
		tctx->csize[kc]=invCR*src->iw*src->ih*src->src_depth[kc]/8;
	}
	free(src);
	return 0;
}
//static int cmp_bpd(const void *left, const void *right)
//{
//	const SampleCtx *pL=(const SampleCtx*)left;
//	const SampleCtx *pR=(const SampleCtx*)left;
//	return pL->bpd0<pR->bpd0;
//}
void batch_test2(void)
{
	enum
	{
		NTIERS=5,
	};
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
	double t=0, total_usize=0, total_csize[4]={0};


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
	if(nthreads<=0)
		nthreads=query_cpu_cores();
	total_usize=0;
	maxlen=0;
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ARRAY_ALLOC(SampleCtx, q, 0, 0, filenames->count, 0);
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		SampleCtx *tctx=(SampleCtx*)ARRAY_APPEND(q, 0, 1, 1, 0);
		tctx->fn=*fn2;
		tctx->idx=k;
		tctx->bpdrank=-1;
	}
	HANDLE *handles=(HANDLE*)malloc(nthreads*sizeof(HANDLE));
	if(!handles)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(handles, 0, nthreads*sizeof(HANDLE));
	for(int k=0;k<(int)filenames->count;k+=nthreads)
	{
		int nt2=k+nthreads<(int)filenames->count?nthreads:(int)filenames->count-k;
		console_log(
			"Analysis %4d/%4d\r"
			, k+nt2
			, (int)filenames->count
		);
		for(int k2=0;k2<(int)nt2;++k2)
		{
			SampleCtx *tctx=(SampleCtx*)array_at(&q, (ptrdiff_t)k+k2);
			handles[k2]=(void*)_beginthreadex(0, 0, analysis, tctx, 0, 0);
			if(!handles[k2])
			{
				LOG_ERROR("Alloc error");
				return;
			}
		}
		WaitForMultipleObjects(nt2, handles, TRUE, INFINITE);
		for(int k2=0;k2<(int)nt2;++k2)
			CloseHandle(handles[k2]);
	}
	console_log("\n");

	//rank images from most to least compressible
	int *indices=(int*)malloc(sizeof(int)*q->count);
	if(!indices)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int k=0;k<(int)q->count;++k)
		indices[k]=k;
	for(int k=0;k<(int)q->count-1;++k)
	{
		int best=k;
		for(int k2=k+1;k2<(int)q->count;++k2)
		{
			SampleCtx *rb=(SampleCtx*)array_at(&q, indices[best]);
			SampleCtx *r2=(SampleCtx*)array_at(&q, indices[k2]);
			if(rb->bpd0>r2->bpd0||(rb->bpd0==r2->bpd0&&indices[k2]<indices[best]))
				best=k2;
		}
		if(best!=k)
		{
			int tmp=indices[best];
			indices[best]=indices[k];
			indices[k]=tmp;
		}
	}
	for(int kr=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		tctx->bpdrank=kr;
	}

	//divide the ranked test into N tiers per usize
	double tsum=0, csum=0;
	for(int kr=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		tsum+=tctx->usize;
	}
	int bpdgroup=0;
	int checkpoints[NTIERS]={0};
	for(int kr=0, kr2=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		int bpdgroup2=(int)(NTIERS*csum/tsum);
		if(bpdgroup2!=bpdgroup)
			checkpoints[kr2++]=kr;
		bpdgroup=bpdgroup2;
		csum+=tctx->usize;
	}
	checkpoints[NTIERS-1]=(int)filenames->count;

	//the batch test
	double tier_usizes[NTIERS]={0};
	double tier_csizes[NTIERS]={0};
	t=time_sec();
	for(int kr=0, kr2=0;kr<(int)q->count;kr+=nthreads)
	{
		int nt2=kr+nthreads<(int)filenames->count?nthreads:(int)filenames->count-kr;
		if(nt2==1)
		{
			SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
			sample_thread(tctx);
		}
		else
		{
			for(int k2=0;k2<(int)nt2;++k2)
			{
				SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr+k2]);
				handles[k2]=(void*)_beginthreadex(0, 0, sample_thread, tctx, 0, 0);
				if(!handles[k2])
				{
					LOG_ERROR("Alloc error");
					return;
				}
			}
			WaitForMultipleObjects(nt2, handles, TRUE, INFINITE);
			for(int k2=0;k2<nt2;++k2)
				CloseHandle(handles[k2]);
		}
		for(int k2=0;k2<nt2;++k2)
		{
			int idx=indices[kr+k2];
			SampleCtx *tctx=(SampleCtx*)array_at(&q, idx);
			double csize=tctx->csize[0]+tctx->csize[1]+tctx->csize[2]+tctx->csize[3];
			ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, tctx->idx);
			console_log(
				"%5d/%5d %s%*sUTYUV %13.2lf %13.2lf %13.2lf %13.2lf %13.2lf  BPD %8.4lf->%8.4lf\n"
				, (int)(idx+1)
				, (int)filenames->count
				, (char*)fn2[0]->data
				, (int)(maxlen-fn2[0]->count+1), ""
				, tctx->usize
				, csize
				, tctx->csize[0]
				, tctx->csize[1]
				, tctx->csize[2]
				, tctx->bpd0
				, 8.*csize/tctx->usize
			);
			tier_usizes[kr2]+=tctx->usize;
			tier_csizes[kr2]+=tctx->csize[0]+tctx->csize[1]+tctx->csize[2]+tctx->csize[3];
			total_usize+=tctx->usize;
			total_csize[0]+=tctx->csize[0];
			total_csize[1]+=tctx->csize[1];
			total_csize[2]+=tctx->csize[2];
			total_csize[3]+=tctx->csize[3];
			if(kr+k2>=checkpoints[kr2]-1)
			{
				++kr2;
				if(kr2<NTIERS)
					console_log("Checkpoint %d\n", kr2);
			}
		}
	}
	t=time_sec()-t;
	{
		char str[1024]={0};
		double ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
		double CR=total_usize/ctotal;
		int nprinted=snprintf(str, sizeof(str)-1
			, "%13.2lf %13.2lf %13.2lf %13.2lf %13.2lf  BPD %8.4lf {%8.4lf %8.4lf %8.4lf %8.4lf %8.4lf}  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB"
			, total_usize
			, ctotal
			, total_csize[0]
			, total_csize[1]
			, total_csize[2]
			, 8./CR
			, tier_csizes[0]/tier_usizes[0]*8
			, tier_csizes[1]/tier_usizes[1]*8
			, tier_csizes[2]/tier_usizes[2]*8
			, tier_csizes[3]/tier_usizes[3]*8
			, tier_csizes[4]/tier_usizes[4]*8
			, t
			, total_usize/(t*1024*1024)
			, t*1024*1024*1000/total_usize
		);
		copy_to_clipboard(str, nprinted);
		console_log("Total UTYUV %s <- copied\n", str);
		timedelta2str(g_buf, G_BUF_SIZE, t);
		console_log("Elapsed %s\n", g_buf);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d-%H:%M:%S");
		console_log("\nDone.  %s\n", g_buf);
		console_pause();
		console_end();
		loud_transforms=1;
	}
	free(handles);
	free(indices);
	array_free(&q);
	array_free(&filenames);
}
