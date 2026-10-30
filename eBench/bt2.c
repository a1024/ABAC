#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
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
		NVAL=2
	};
	int amin[4]={0}, amax[4]={0};
	int nch=0, depth=0;
	int hsize=0;
	int32_t *hists=0;
	int psize=0;
	int16_t *pixels=0;
	SampleCtx *tctx=0;
	Image *src=0;
	int bestrct=0, yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
	int idx=0, ky=0, kx=0, kc=0;
	int32_t *data=0;
	int yuv[4]={0};

	tctx=(SampleCtx*)param;
	src=image_load((char*)HDR_DATA(tctx->fn), (int)tctx->fn->count);
	if(!src)
		return 0;
	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
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
		CRASH("ALloc error");
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
	data=(int32_t*)HDR_DATA(src);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4]={0};

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0;

			yuv[0]=data[idx+yidx];
			yuv[1]=data[idx+uidx];
			yuv[2]=data[idx+vidx];
			yuv[3]=data[idx+3];
			for(kc=0;kc<4;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int16_t
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						eNEE	=rows[1][1+3*NCH*NROWS*NVAL],
						eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
						eW	=rows[0][1-1*NCH*NROWS*NVAL];
					int ctx=31^LZCNT32(eW*eW+1);
					int pred=N+W-NW, vmax=N, vmin=W;
					int error2, sym;

					if(ctx>NCTX-1)
						ctx=NCTX-1;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(pred, vmin, vmax);
					pred+=offset;
					CLAMP2(pred, amin[kc], amax[kc]);

					/*fwd*/
					error2=yuv[kc]-pred;
					error2<<=32-src->depth[kc];
					error2>>=32-src->depth[kc];
					sym=error2<<1^error2>>31;
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
	}
	free(pixels);
	{
		double csize=0;
		int nlevels=1<<depth;
		int ks=0;

		for(kc=0;kc<NCTX*nch;++kc)
		{
			int32_t *currhist=0, sum=0;
			double invsum=0;

			currhist=hists+((ptrdiff_t)kc<<depth);
			for(ks=0;ks<nlevels;++ks)
				sum+=currhist[ks];
			if(!sum)
				continue;
			invsum=1./sum;
			for(ks=0;ks<nlevels;++ks)
			{
				int freq=currhist[ks];
				if(freq)
					csize-=freq*log(freq*invsum);
			}
		}
		free(hists);
		csize*=1./M_LN2;

		/*BPD = (csize / usize) * depth*/
		tctx->usize=(double)src->iw*src->ih*nch*depth;
		tctx->bpd0=csize/((double)src->iw*src->ih*nch);
		free(src);
	}
	return 0;
}
static unsigned __stdcall sample_thread(void *param)
{
	double entropy[4]={0};
	SampleCtx *tctx=0;
	Image *src=0;
	int kc=0;

	tctx=(SampleCtx*)param;
	src=image_load((char*)HDR_DATA(tctx->fn), (int)tctx->fn->count);
	if(!src)
		return 0;
	tctx->usize=image_getBMPsize(src);
	apply_selected_transforms(&src, 0, 1, 1);
	calc_csize_stateful(src, 0, entropy);
	for(kc=0;kc<4;++kc)
	{
		int depth=src->src_depth[kc];
		double invCR=depth?entropy[kc]:0;
		tctx->csize[kc]=invCR*src->iw*src->ih*src->src_depth[kc]/8;
	}
	free(src);
	return 0;
}
#if 0
static int cmp_bpd(const void *left, const void *right)
{
	const SampleCtx *pL=(const SampleCtx*)left;
	const SampleCtx *pR=(const SampleCtx*)left;
	return pL->bpd0<pR->bpd0;
}
#endif
void batch_test2(int custom)
{
	enum
	{
		NTIERS=5
	};
	const char cfgfn[]="ebench.txt";
	const char *ext[]=
	{
		"PPM", "PGM", "PNM",
		"PNG",
		"JPG", "JPEG",
		"BMP",
		"TIF", "TIFF",
	};
	ArrayHandle path=0, filenames=0, q=0, paths=0;
	int nthreads=0, maxlen=0;
	double t=0, total_usize=0, total_csize[4]={0};
	int user=0;
	HANDLE *handles=0;
	int *indices=0;
	int k=0, kr=0, kr2=0;
	double tsum=0, csum=0;
	int bpdgroup=0;
	int checkpoints[NTIERS]={0};
	double tier_usizes[NTIERS]={0};
	double tier_csizes[NTIERS]={0};

	DisableProcessWindowsGhosting();
	console_start();
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H:%M:%S");
	console_log("Batch Test  %s\n", g_buf);
	if(!custom)
	{
		FILE *f=fopen(cfgfn, "r");
		if(f)
		{
			for(;;)
			{
				char buf[1024], *ptr;
				ArrayHandle fns;

				memset(buf, '\n', sizeof(buf));
				if(!fgets(buf, sizeof(buf)-1, f))
					break;
				ptr=buf;
				while(*ptr!='\n')++ptr;
				*ptr=0;
				fns=get_filenames(buf, ext, _countof(ext), 1);
				if(fns)
				{
					array_append(&filenames, HDR_DATA(fns), sizeof(ArrayHandle*), fns->count, 1, 0, free_str);
					free(fns);
				}
			}
			fclose(f);
			if(!filenames||!filenames->count)
			{
				console_log("No data found in \"%s\", please add data paths there or remove the file\n", cfgfn);
				console_pause();
				console_end();
				return;
			}
			goto start;
		}
	}
	for(filenames=0;;)
	{
		ArrayHandle fns;

		path=dialog_open_folder();
		if(!path)
			break;
		fns=get_filenames((char*)HDR_DATA(path), ext, _countof(ext), 1);
		if(!fns)
		{
			array_free(&path);
			break;
		}
		console_log("\"%s\"\n", (char*)HDR_DATA(path));
		array_append(&paths, &path, sizeof(ArrayHandle*), 1, 1, 0, free_str);
#if 0
		array_free(&path);
#endif
		array_append(&filenames, HDR_DATA(fns), sizeof(ArrayHandle*), fns->count, 1, 0, free_str);
		free(fns);/*shallow free*/
	}
	if(!filenames||!filenames->count)
	{
		console_end();
		return;
	}
	user=1;
start:
	loud_transforms=0;
	console_log("Enter number of threads: ");
	nthreads=console_scan_int();
	if(nthreads<=0)
		nthreads=query_cpu_cores();

	total_usize=0;
	maxlen=0;
	for(k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		if(maxlen<(int)fn2[0]->count)
			maxlen=(int)fn2[0]->count;
	}
	ARRAY_ALLOC(SampleCtx, q, 0, 0, filenames->count, 0);
	for(k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, k);
		SampleCtx *tctx=(SampleCtx*)ARRAY_APPEND(q, 0, 1, 1, 0);
		tctx->fn=*fn2;
		tctx->idx=k;
		tctx->bpdrank=-1;
	}
	handles=(HANDLE*)malloc(nthreads*sizeof(HANDLE));
	if(!handles)
	{
		CRASH("Alloc error");
		return;
	}
	memset(handles, 0, nthreads*sizeof(HANDLE));
	for(k=0;k<(int)filenames->count;k+=nthreads)
	{
		int nt2, k2;

		nt2=k+nthreads<(int)filenames->count?nthreads:(int)filenames->count-k;
		console_log(
			"Analysis %4d/%4d\r"
			, k+nt2
			, (int)filenames->count
		);
		for(k2=0;k2<(int)nt2;++k2)
		{
			SampleCtx *tctx=(SampleCtx*)array_at(&q, (ptrdiff_t)k+k2);
			handles[k2]=(void*)_beginthreadex(0, 0, analysis, tctx, 0, 0);
			if(!handles[k2])
			{
				CRASH("Alloc error");
				return;
			}
		}
		WaitForMultipleObjects(nt2, handles, TRUE, INFINITE);
		for(k2=0;k2<(int)nt2;++k2)
			CloseHandle(handles[k2]);
	}
	console_log("\n");

	/*rank images from most to least compressible*/
	indices=(int*)malloc(sizeof(int)*q->count);
	if(!indices)
	{
		CRASH("Alloc error");
		return;
	}
	for(k=0;k<(int)q->count;++k)
		indices[k]=k;
	for(k=0;k<(int)q->count-1;++k)
	{
		int best, k2;

		best=k;
		for(k2=k+1;k2<(int)q->count;++k2)
		{
			SampleCtx *rb=(SampleCtx*)array_at(&q, indices[best]);
			SampleCtx *r2=(SampleCtx*)array_at(&q, indices[k2]);
			if(rb->bpd0>r2->bpd0||(fabs(rb->bpd0-r2->bpd0)<1e-9&&indices[k2]<indices[best]))
				best=k2;
		}
		if(best!=k)
		{
			int tmp=indices[best];
			indices[best]=indices[k];
			indices[k]=tmp;
		}
	}
	for(kr=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		tctx->bpdrank=kr;
	}

	/*divide the ranked test into N tiers per usize*/
	for(kr=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		tsum+=tctx->usize;
	}
	for(kr=0, kr2=0;kr<(int)q->count;++kr)
	{
		SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
		int bpdgroup2=(int)(NTIERS*csum/tsum);
		if(bpdgroup2!=bpdgroup)
			checkpoints[kr2++]=kr;
		bpdgroup=bpdgroup2;
		csum+=tctx->usize;
	}
	checkpoints[NTIERS-1]=(int)filenames->count;

	/*the batch test*/
	g_automated=1;
	t=time_sec();
	for(kr2=kr=0;kr<(int)q->count;kr+=nthreads)
	{
		int nt2, k2;

		nt2=kr+nthreads<(int)filenames->count?nthreads:(int)filenames->count-kr;
		if(nt2==1)
		{
			SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr]);
			sample_thread(tctx);
		}
		else
		{
			for(k2=0;k2<(int)nt2;++k2)
			{
				SampleCtx *tctx=(SampleCtx*)array_at(&q, indices[kr+k2]);
				handles[k2]=(void*)_beginthreadex(0, 0, sample_thread, tctx, 0, 0);
				if(!handles[k2])
				{
					CRASH("Alloc error");
					return;
				}
			}
			WaitForMultipleObjects(nt2, handles, TRUE, INFINITE);
			for(k2=0;k2<nt2;++k2)
				CloseHandle(handles[k2]);
		}
		for(k2=0;k2<nt2;++k2)
		{
			int idx=indices[kr+k2];
			SampleCtx *tctx=(SampleCtx*)array_at(&q, idx);
			double csize=tctx->csize[0]+tctx->csize[1]+tctx->csize[2]+tctx->csize[3];
			ArrayHandle *fn2=(ArrayHandle*)array_at(&filenames, tctx->idx);
			console_log(
				"%5d %5d/%5d %s%*sUTYUV %13.2f %13.2f %13.2f %13.2f %13.2f  BPD %8.4f->%8.4f\n"
				, (int)(idx+1)
				, kr+k2
				, (int)filenames->count
				, (char*)HDR_DATA(fn2[0])
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
					console_log("Checkpoint %d/%d\n", kr2, NTIERS);
			}
		}
	}
	g_automated=0;
	t=time_sec()-t;
	{
		char str[1024]={0};
		double ctotal=total_csize[0]+total_csize[1]+total_csize[2]+total_csize[3];
		double CR=total_usize/ctotal;
		int nprinted=snprintf(str, sizeof(str)-1
			, "%13.2f %13.2f %13.2f %13.2f %13.2f  BPD %8.6f {%8.6f %8.6f %8.6f %8.6f %8.6f}  %12.6f sec  %12.6f MB/s  %12.6f ms/MB"
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
		loud_transforms=1;
	}
	if(!custom&&user)
	{
		FILE *f=fopen(cfgfn, "w");
		if(!f)
		{
			console_log("Cannot open \"%s\" for writing\n", cfgfn);
			console_pause();
			console_end();
			return;
		}
		for(k=0;k<(int)paths->count;++k)
		{
			ArrayHandle *fn2=(ArrayHandle*)array_at(&paths, k);
			fwrite((char*)HDR_DATA(fn2[0]), 1, fn2[0]->count, f);
			fputc('\n', f);
		}
		fclose(f);
	}
	array_free(&paths);
	console_end();
	free(handles);
	free(indices);
	array_free(&q);
	array_free(&filenames);
}
