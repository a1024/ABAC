#define AC_IMPLEMENTATION
#include"ac.h"
#include<stdio.h>
#include<string.h>
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
static const char file[]=__FILE__;


//	#define AC_PRINT_HITCOUNT
//	#define AC_PRINT_PROB
//	#define AC_DISABLE_BYPASS


static const int tag_ac00='A'|'C'<<8|'0'<<16|'0'<<24;

#ifdef DEBUG_ABAC2
static int examined_plane=9, examined_start=0, examined_end=100;
#endif
long long ac0_encode(const void *src, ptrdiff_t nbytes, int bitoffset, int bitdepth, int bytestride, ArrayHandle *out, int loud)
{
	long long cycles=__rdtsc();

	if(!src||!nbytes||!bitdepth||!bytestride||!out)
		LOG_ERROR("binac0_encode(src=%p, nbytes=%d, bitdepth=%d, stride=%d, out=%p)", src, nbytes, bitdepth, bytestride, out);

	const unsigned char *buf=(const unsigned char*)src;
	unsigned char *prob=(unsigned char*)malloc(bitdepth);
	size_t nsym=nbytes/bytestride;
	for(int kp=0;kp<bitdepth;++kp)//analyze bitplanes
	{
		int startbitidx=bitoffset+kp;
		const unsigned char
			*srcptr=buf+(startbitidx>>3),
			*srcend=buf+nbytes;
		size_t freq=nsym;
		//ptrdiff_t step=bytestride<<3;//needed when bytestride < 1 byte
		for(ptrdiff_t sh=startbitidx&7;srcptr<srcend;srcptr+=bytestride)
		{
			unsigned char bit=*srcptr>>sh&1;
			freq-=bit;//freq=P(0)

			//sh=(sh+step)&7;//needed when bytestride < 1 byte
		}
		int qfreq=(int)(((freq<<8)+(nsym>>1))/nsym);//round(p0/size)
		qfreq+=(qfreq<0x80)-(qfreq>0x80);//0x80 means bypass

		//qfreq-=(freq<nsym)&(qfreq>255);
		//--qfreq;//add 1 when assigning p0
		//if(qfreq<1)//0x80 qfreq means bypass, which is triggered when ratio < 1
		//	qfreq=1;

		prob[kp]=qfreq;//adaptive prediction must be opposite to encoding direction
#ifdef AC_PRINT_PROB
		if(loud)
			printf("bit %d: P(0) = %7d/%7d = 0x%02X\n", kp, (int)freq, (int)nsym, qfreq);
#endif
	}

	size_t dststartidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return 0;
		ARRAY_APPEND(*out, 0, 0, 1, nbytes);
		dststartidx=out[0]->count;
	}
	else
	{
		ARRAY_ALLOC(char, *out, 0, 0, nbytes, 0);
		dststartidx=0;
	}
	ARRAY_APPEND(*out, 0, 12+bitdepth, 1, 0);//{tag, total written bytes, prob array}
	memcpy(out[0]->data+dststartidx, &tag_ac00, 4);

#ifdef AC_PRINT_HITCOUNT
	size_t hitcount=0;
#endif
	DList list={0};
	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
#ifdef AC_PRINT_HITCOUNT
		size_t hitcount_p=0;
#endif
		unsigned char p0=prob[kp];
		int startidx=(bitoffset+kp)>>3, sh=(bitoffset+kp)&7;
#ifndef AC_DISABLE_BYPASS
		if(p0!=0x80)
#endif
		{
			unsigned r_start=0, r_end=0xFFFFFFFF;
			dlist_init(&list, 1, 1024, 0);
			for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
			{
				unsigned middle=r_start+(unsigned)((unsigned long long)(r_end-r_start)*p0>>8);

				int bit=buf[startidx+ks]>>sh&1;
				if(bit)
					r_start=middle+1;
				else
					r_end=middle-1;
#ifdef AC_PRINT_HITCOUNT
				int correct=bit^(p0>=0x80);
				hitcount_p+=correct;
#endif

				while((r_start^r_end)<0x1000000)//Matt Mahoney
				{
					dlist_push_back1(&list, (char*)&r_start+3);//shift out most-significant byte
					r_start<<=8;
					r_end=r_end<<8|0xFF;
				}

				if(r_start+3<r_start||r_start+3>r_end)
				{
					dlist_push_back1(&list, (char*)&r_start+3);//big endian
					dlist_push_back1(&list, (char*)&r_start+2);
					dlist_push_back1(&list, (char*)&r_start+1);
					dlist_push_back1(&list, (char*)&r_start);

					r_start=0, r_end=0xFFFFFFFF;//because 1=0.9999...
				}
			}
			dlist_push_back1(&list, (char*)&r_start+3);//big endian
			dlist_push_back1(&list, (char*)&r_start+2);
			dlist_push_back1(&list, (char*)&r_start+1);
			dlist_push_back1(&list, (char*)&r_start);
		}
		size_t csize;
#ifndef AC_DISABLE_BYPASS
		if(p0==0x80||list.nobj>nsym)//ratio > 1, bypass
		{
			prob[kp]=0x80;
			size_t bypass_idx=out[0]->count, bypass_len=(nsym+7)>>3;
			ARRAY_APPEND(*out, 0, bypass_len, 1, 0);
			unsigned char frag;
			ptrdiff_t ks=0, end=nbytes-bytestride*7, kb=0;
			for(;ks<end;++kb)
			{
				frag=0;
				frag|= buf[startidx+ks]>>sh&1,     ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<1, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<2, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<3, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<4, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<5, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<6, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<7, ks+=bytestride;
				out[0]->data[bypass_idx+kb]=frag;
			}
			if(ks<nbytes)
			{
				frag=0;
				int sh2=0;
				do
				{
					frag|=(buf[startidx+ks]>>sh&1)<<sh2;
					ks+=bytestride;
					++sh2;
				}
				while(ks<nbytes);
				out[0]->data[bypass_idx+kb]=frag;
			}
			csize=nsym;
		}
		else
#endif
		{
			csize=list.nobj;
			dlist_appendtoarray(&list, out);
		}
		if(p0!=0x80)
			dlist_clear(&list);

#ifdef AC_PRINT_HITCOUNT
		if(loud)
			printf("bit %d: r =%6d /%6d = %lf, hit=%6d=%lf%%\n", kp, (int)nsym, (int)(csize<<3), (double)nsym/(csize<<3), (int)hitcount_p, 100.*hitcount_p/nsym);
		hitcount+=hitcount_p;
#else
		if(loud)
			printf("bit %d: r =%6d /%6d = %lf\n", kp, (int)nsym, (int)csize, (double)nsym/csize);
#endif
	}
	size_t byteswritten=out[0]->count-dststartidx;
	memcpy(out[0]->data+dststartidx+4, &byteswritten, 8);
	memcpy(out[0]->data+dststartidx+12, prob, bitdepth);
	cycles=__rdtsc()-cycles;

	if(loud)
	{
		size_t original_bitsize=nsym*bitdepth, compressed_bitsize=byteswritten<<3;
		printf("AC encode:  %lld cycles, %lf c/byte\n", cycles, (double)(cycles<<3)/original_bitsize);
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", (int)original_bitsize>>3, (int)compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize, (double)compressed_bitsize/nsym);
#ifdef AC_PRINT_HITCOUNT
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitcount, original_bitsize, 100.*hitcount/original_bitsize);
#endif
		
		printf("Preview:\n");
		int kprint=200;
		if(byteswritten<kprint)
			kprint=(int)byteswritten;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out[0]->data[dststartidx+k]&0xFF);
		printf("\n");
	}
	return byteswritten;
}
const void* ac0_decode(const void *srcstart, const void *srcend, void *dst, size_t nbytes, int bitoffset, int bitdepth, int bytestride, int loud
#ifdef ENABLE_GUIDE
	, unsigned char *guide
#endif
)//set the dst buffer to zero
{
	if(!srcstart||!srcend||!dst||!nbytes||!bitdepth||!bytestride)
		LOG_ERROR("abac4_decode(in_start=%p, in_end=%p, dst=%p, imsize=%d, bitdepth=%d, stride=%d)", srcstart, srcend, dst, nbytes, bitdepth, bytestride);

	long long cycles=__rdtsc();

	const unsigned char
		*data=(const unsigned char*)srcstart,
		*srcptr=data+12+bitdepth,
		*dataend=(const unsigned char*)srcend;
	size_t datalen=dataend-data;
	unsigned char *buf=(unsigned char*)dst;

	if(datalen<12)
		LOG_ERROR("AC: datalen = %p", datalen);
	int tag;
	memcpy(&tag, data, 4);
	if(tag!=tag_ac00)
		LOG_ERROR("AC: invalid tag 0x%08X, expected 0x%08X", tag, tag_ac00);

	size_t csize;
	memcpy(&csize, data+4, 8);
	if(datalen<csize)
		LOG_ERROR("AC: Unexpected EOF, datalen = %p", datalen);
	dataend=data+csize;

	unsigned char *prob=(unsigned char*)malloc(bitdepth);
	memcpy(prob, data+12, bitdepth);

	size_t nsym=nbytes/bytestride;
	for(int kp=bitdepth-1;kp>=0;--kp)
	{
		int startidx=(bitoffset+kp)>>3, sh=(bitoffset+kp)&7;
		unsigned char p0=prob[kp];
#ifndef AC_DISABLE_BYPASS
		if(p0==0x80)//bypass
		{
			ptrdiff_t ks=0, end=nbytes-bytestride*7;
			for(;ks<end;++srcptr)
			{
				unsigned char frag=*srcptr;
				buf[startidx+ks]|=(frag>>0&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>1&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>2&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>3&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>4&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>5&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>6&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>7&1)<<sh, ks+=bytestride;
			}
			if(ks<(ptrdiff_t)nbytes)
			{
				unsigned char frag=*srcptr;
				int sh2=0;
				do
				{
					buf[startidx+ks]|=(frag>>sh2&1)<<sh;
					ks+=bytestride;
					++sh2;
				}while(ks<(ptrdiff_t)nbytes);
				++srcptr;
			}
#ifdef ENABLE_GUIDE
			for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
			{
				unsigned char
					b1=guide[startidx+ks]>>sh&1,
					b2=buf[startidx+ks]>>sh&1;
				if(b2!=b1)
					LOG_ERROR("AC bypass error kp%2d ks%7d dec %d != %d", kp, ks, b2, b1);
			}
#endif
		}
		else
#endif
		{
			unsigned r_start=0, r_end=0xFFFFFFFF, code;

#ifdef ENABLE_GUIDE
			unsigned r_start2=0, r_end2=0xFFFFFFFF;
			size_t enc_renorms=0, dec_renorms=0;
#endif
			srcptr+=4;
			if(srcptr>dataend)
				LOG_ERROR("AC OOB [start] srcptr %p dataend %p", srcptr, dataend);
			code=srcptr[-4]<<24|srcptr[-3]<<16|srcptr[-2]<<8|srcptr[-1];//big endian

			//for(ptrdiff_t ks=0;;)
			for(ptrdiff_t ks=0;ks<(ptrdiff_t)nbytes;ks+=bytestride)
			{
				unsigned middle=r_start+(unsigned)((unsigned long long)(r_end-r_start)*p0>>8);

				unsigned char bit=code>middle;
#ifdef ENABLE_GUIDE
				unsigned mid2=r_start2+(unsigned)((unsigned long long)(r_end2-r_start2)*p0>>8);
				unsigned char b0=guide[startidx+ks]>>sh&1;
				if(bit!=b0||enc_renorms!=dec_renorms||r_start!=r_start2||r_end!=r_end2)
					LOG_ERROR("kp%2d ks%7d dec %d != %d", kp, ks, bit, b0);
				if(b0)
					r_start2=mid2+1;
				else
					r_end2=mid2-1;
				while((r_start2^r_end2)<0x1000000)
				{
					r_start2<<=8;
					r_end2=r_end2<<8|0xFF;
					++enc_renorms;
				}
				if(r_start2+3<r_start2||r_start2+3>r_end2)
				{
					r_start2=0, r_end2=0xFFFFFFFF;//because 1=0.9999...
					enc_renorms+=4;
				}
#endif

				buf[startidx+ks]|=bit<<sh;

				//ks+=bytestride;
				//if(ks>=nbytes)
				//	break;

				if(bit)
					r_start=middle+1;
				else
					r_end=middle-1;

				//if(srcptr<dataend)
				//{
					while((r_start^r_end)<0x1000000)
					{
						++srcptr;
						if(srcptr>dataend)
							LOG_ERROR("AC OOB kp%2d ks%7d srcptr %p dataend %p", kp, ks, srcptr, dataend);

						//++srcptr;
						code=code<<8|srcptr[-1];//shift out most-significant byte
						r_start<<=8;
						r_end=r_end<<8|0xFF;
#ifdef ENABLE_GUIDE
						++dec_renorms;
#endif
					}
				//}
				//else
				//	srcptr=dataend;

				if(r_start+3<r_start||r_start+3>r_end)
				{
					srcptr+=4;
					if(srcptr>dataend)
						LOG_ERROR("AC OOB kp%2d ks%7d srcptr %p dataend %p", kp, ks, srcptr, dataend);

					code=srcptr[-4]<<24|srcptr[-3]<<16|srcptr[-2]<<8|srcptr[-1];//big endian

					r_start=0, r_end=0xFFFFFFFF;//because 1=0.9999...
#ifdef ENABLE_GUIDE
					dec_renorms+=4;
#endif
				}
			}
		}
	}
	free(prob);
	cycles=__rdtsc()-cycles;

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", cycles, (double)(cycles<<3)/(nsym*bitdepth));
	return srcptr;
}




	#define T26_OPT_PRED

//test26: T16 with arithmetic range coder
void t25_calchist(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, unsigned *hist);
void t26_normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	if(nsymbols)
	{
		int sum=0;
		for(int sym=0;sym<nlevels;++sym)
		{
			int qfreq=((long long)srchist[sym]*0xFFFF)/nsymbols;
			CDF[sym]=sum;
			sum+=qfreq;
		}
	}
	else//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
	}
}
void t25_addhist(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0a, int x0b, int y0, int maxinc, unsigned *CDF2);
int t26_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int iw, int ih, int kc, int x1, int x2, int y, T26Params const *p, unsigned *CDF2, int rec)
{
	int overflow=0;
	int sum, cdf1, f1, f2, freq;
	memset(CDF2, 0, 257*sizeof(unsigned));
	if(p->mtop)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x2+p->mright, y-p->mtop, y, x1, x2, y, p->maxinc, CDF2);
	if(p->mleft)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x1, y, y+1, x1, x2, y, p->maxinc, CDF2);

	if(CDF2[256])
	{
		sum=0;
		for(int sym=0;sym<256;++sym)
		{
			if(sym==0x80)//
				sym=0x80;

			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize
			
			if(f2<f1)
				freq=f2+(int)(((long long)f1-f2)*(0xFF-p->alpha)/0xFF);//blend
			else
				freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

			int f3=freq;//

			freq=(int)((long long)freq*0xFF00>>16)+1;//guard
			//freq=CLAMP(0, freq, 0xFF01);

			if(rec)//
				printf("%3d 0x%04X 0x%04X 0x%04X 0x%04X\n", sym, f1, f2, f3, freq);

			if(freq<0||freq>0xFF01)
			{
				printf("Impossible freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
				return 0;
			}
				//LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000&&sym<255)
			{
				if(!rec)//
				{
					t26_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, y, p, CDF2, 1);//
					//for(int k=0;k<=sym;++k)
					//	printf("%3d 0x%04X 0x%04X\n", k, CDF0[k], CDF2[k]);
					printf("ANS CDF sym 0x%02X sum 0x%04X  freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", sym, sum, freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
					//printf("ANS CDF sym 0x%02X sum 0x%04X freq 0x%04X\n", sym, sum, freq);
				}
				return 0;
			}
		}
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
	}
	CDF2[256]=0x10000;
	return 1;
}
int t26_encode(const unsigned char *src, int iw, int ih, T26Params const *params, int use_ans, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	double t_start=time_ms();
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifndef T26_OPT_PRED
	apply_transforms_fwd(buf2, iw, ih);
#else
	addbuf(buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd((char*)buf2, iw, ih);
	pred_opt_opt_v6((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v5((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v4((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v3((char*)buf2, iw, ih, loud);
	//pred_opt_opt_v2((char*)buf2, iw, ih);
#endif
	short predparams[22+PW2_NPARAM];
	const short
		predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM},
		predidx[]={0, 11, 11+PW2_NPARAM};
	memcpy(predparams+predidx[0], jxlparams_i16,         predlen[0]*sizeof(short));
	memcpy(predparams+predidx[1], pw2_params+PW2_NPARAM, predlen[1]*sizeof(short));
	memcpy(predparams+predidx[2], jxlparams_i16+22,      predlen[2]*sizeof(short));
#ifdef T26_OPT_PRED
	//pred_jxl_opt_v2((char*)buf2, iw, ih, jxlparams_i16, loud);
#if 1
	if(loud)
		pred_opt_printparam();
#endif
	pred_opt_apply((char*)buf2, iw, ih, 1);
	//pred_jxl_apply((char*)buf2, iw, ih, jxlparams_i16, 1);
	addbuf(buf2, iw, ih, 3, 4, 128);
#endif
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, CDF2);
		t26_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int bookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, predparams, (predlen[3])*sizeof(short));
	//dlist_push_back(&list, jxlparams_i16, 33*sizeof(short));
	dlist_push_back(&list, CDF0, 768*sizeof(short));
	
	for(int kc=0;kc<3;++kc)//for each channel
	{
		T26Params const *p=params+kc;
		int gxcount=(iw+p->gwidth-1)/p->gwidth;

		if(use_ans)
		{
			ANSEncContext ctx;
			ans_enc_init(&ctx, CDF2, &list);

			for(int ky=ih-1;ky>=0;--ky)//for each row
			{
				for(int bx=gxcount-1;bx>=0;--bx)//for each group
				{
					//if(kc==0&&bx==1&&ky==0)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					t26_prepblock(buf2, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);

					//if(kc==0&&bx==0&&by==0)
					//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

					//encode group
					for(int kx=x2-1;kx>=x1;--kx)//for each pixel
					{
						if(kc==0&&ky==0&&bx==0&&kx==0)//
							printf("");

						ans_enc(&ctx, buf2[(iw*ky+kx)<<2|kc], kc);
					}
				}
			}
			ans_enc_flush(&ctx);
		}
		else
		{
			ACEncContext ctx;
			ac_enc_init(&ctx, CDF2, &list);

			//unsigned state_lo=0, state_hi=0xFFFFFFFF;
			//unsigned cache=0;
			//int nbits=0;
			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					//if(kc==0&&bx==1&&ky==0)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					t26_prepblock(buf2, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);

					//if(kc==0&&bx==0&&by==0)
					//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

					//encode group
					for(int kx=x1;kx<x2;++kx)//for each pixel
					{
						//if(kc==2&&ky2==0&&kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);
					
						//if(acval&&acval->count==2525)//
						//	printf("");
						//ac_enc_renorm(&state_lo, &state_hi, &cache, &nbits, &list);
	#if 0
						unsigned rlo, rhi, mingap;
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							//if(mingap>0&&state_lo+0x10000>state_hi)
							//	printf("");
							if(state_lo<state_hi&&mingap>0)
								break;
							if(nbits>=32)
							{
								dlist_push_back(&list, &cache, 4);
								cache=0;
								nbits=0;
							}
							cache|=(state_lo&0x80000000)>>nbits;//cache is written MSB -> LSB
							++nbits;
							//cache|=(state_lo>>31)<<32-nbits;

							state_lo<<=1;//shift out MSB
							state_hi<<=1;
							state_hi|=1;

	#if 0
							cache<<=1;
							cache|=state_lo>>31;
							state_lo<<=1;//shift out MSB
							state_hi<<=1;
							state_hi|=1;
							++nbits;
	#endif
						}
	#endif
					
						ac_enc(&ctx, buf2[(iw*ky+kx)<<2|kc]);
					}
				}
			}
			ac_enc_flush(&ctx);
#if 0
			int k2=0;
			do//flush
			{
				while(nbits<32)
				{
					cache|=(state_lo&0x80000000)>>nbits;//cache is written MSB -> LSB
					++nbits;
					++k2;

					state_lo<<=1;//shift out MSB
					state_hi<<=1;
					state_hi|=1;
				}
				dlist_push_back(&list, &cache, 4);
				cache=0;
				nbits=0;
			}while(k2<32);
#endif
		}
		bookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, bookmarks, 12);
	
	int overhead=12+(int)(768*sizeof(short));
	int ch[]=
	{
		bookmarks[0]-overhead,
		bookmarks[1]-bookmarks[0],
		bookmarks[2]-bookmarks[1],
	};
	//if(csizes)
	//{
	//	csizes[0]=ch[0];
	//	csizes[1]=ch[1];
	//	csizes[2]=ch[2];
	//}
	if(loud)
	{
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_ms()-t_start);
		printf("\n");
		printf("Total    %7d  %lf\n", bookmarks[2], 3.*res/bookmarks[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf\n", ch[0], (double)res/ch[0]);
		printf("Green    %7d  %lf\n", ch[1], (double)res/ch[1]);
		printf("Blue     %7d  %lf\n", ch[2], (double)res/ch[2]);
	}

	dlist_clear(&list);
	free(CDF2);
	free(CDF0);
	free(buf2);
	return 1;
}
int t26_decode(const unsigned char *data, size_t srclen, int iw, int ih, T26Params const *params, int use_ans, unsigned char *buf, int loud)
{
	const short predidx[]={0, 11, 11+PW2_NPARAM}, predlen[]={11, PW2_NPARAM, 11, 22+PW2_NPARAM};
	const int cdflen=768LL*sizeof(short), overhead=12LL+predlen[3]*sizeof(short)+cdflen;
	int res=iw*ih;
	
	const unsigned char *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}

	unsigned bookmarks[3];
	memcpy(bookmarks, data, 12);
	if(bookmarks[2]<(unsigned)overhead||bookmarks[2]>srclen)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short predparams[22+PW2_NPARAM];
	memcpy(predparams, data+12, predlen[3]*sizeof(short));
	memcpy(CDF0, data+12+predlen[3]*sizeof(short), cdflen);

#ifdef AC_VALIDATE
	acval_idx=0;
#endif
	for(int kc=0;kc<3;++kc)//for each channel
	{
		T26Params const *p=params+kc;
		int gxcount=(iw+p->gwidth-1)/p->gwidth;

		if(use_ans)
		{
			ANSDecContext ctx;
			ans_dec_init(&ctx, CDF2, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			//unsigned state_lo=0, state_hi=0xFFFFFFFF, code, cache;
			//int nbits=32;
			//srcptr=kc?data+bookmarks[kc-1]:data+overhead;
			//srcend=data+bookmarks[kc];

			//if(ctx.srcend-ctx.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ctx.code, ctx.srcptr, 4);
			//ctx.srcptr+=4;
			//
			//if(ctx.srcend-ctx.srcptr<4)
			//	LOG_ERROR("buffer overflow");
			//memcpy(&ctx.cache, ctx.srcptr, 4);
			//ctx.srcptr+=4;

			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					//if(kc==0&&ky==0&&bx==1)//
					//	kc=0;
					//if(kc==1)//
					//	kc=1;

					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					int success=t26_prepblock(buf, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);
					if(!success)
						LOG_ERROR("t26_prepblock error");
					for(int kx=x1;kx<x2;++kx)//for each pixel
					{
						buf[(iw*ky+kx)<<2|kc]=ans_dec(&ctx, kc);
#if 0
						unsigned rlo, rhi, mingap;

						//if(acval_idx==2525)//
						//	printf("");
					
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							if(state_lo<state_hi&&mingap>0)
								break;
							if(!nbits)
							{
								if(srcend-srcptr<4)
								{
#ifdef AC_VALIDATE
									printf("buffer overflow\n");
									acval_dump();
#endif
									LOG_ERROR("buffer overflow");
								}
								memcpy(&cache, srcptr, 4);
								srcptr+=4;

								nbits=32;
							}
							--nbits;
							code<<=1;//shift out MSB		cache is read MSB -> LSB
							code|=(unsigned)(cache>>nbits&1);

							state_lo<<=1;
							state_hi<<=1;
							state_hi|=1;
						}
#if 0
						for(;;)//renorm
						{
							mingap=(state_hi-state_lo)>>16;
							if(mingap>0)
								break;
							if(!nbits)
							{
								nbits+=32;
								cache<<=32;

								if(srcptr+4>=srcend)
									LOG_ERROR("buffer overflow");
								memcpy(&cache, srcptr, 4);
								srcptr+=4;
							}
							--nbits;
							code<<=1;//shift out MSB
							code|=(unsigned)(cache>>nbits&1);

							state_lo<<=1;
							state_hi<<=1;
							state_hi|=1;
						}
#endif
						int sym=0;
						int L=0, R=256, found=0;
						unsigned code2;
						while(L<=R)
						{
							sym=(L+R)>>1;
							code2=state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym]>>16);
							if(code2<code)
								L=sym+1;
							else if(code2>code)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(found)
							for(;sym<256-1&&state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym+1]>>16)==code;++sym);
						else
							sym=L+(L<256&&state_lo+((unsigned long long)(state_hi-state_lo)*CDF2[sym+1]>>16)<code)-(L!=0);
#if 0
						//unsigned c=(unsigned)((((long long)(code-state_lo)<<16)+((state_hi-state_lo)>>1))/(state_hi-state_lo));
						unsigned c=(unsigned)(((long long)(code-state_lo)<<16)/(state_hi-state_lo));
						int sym=0;
					
						//if(kc==0&&ky==2&&kx==512)//
						//	printf("");

						int L=0, R=256, found=0;
						while(L<=R)
						{
							sym=(L+R)>>1;
							if(CDF2[sym]<c)
								L=sym+1;
							else if(CDF2[sym]>c)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(found)
							for(;sym<256-1&&CDF2[sym+1]==c;++sym);
						else
							sym=L+(L<256&&CDF2[L]<c)-(L!=0);
#endif
						buf[(iw*ky+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf_start=CDF2[sym], cdf_end=CDF2[sym+1];
					
						rlo=state_lo+((unsigned long long)(state_hi-state_lo)*cdf_start>>16);
						rhi=state_lo+((unsigned long long)(state_hi-state_lo)*cdf_end  >>16);
						acval_dec(sym, cdf_start, cdf_end, state_lo, state_hi, rlo, rhi, cache, nbits, code);//
						state_lo=rlo;
						state_hi=rhi-1;//OBLIGATORY range leak guard
#endif
						//debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
						//state=freq*(state>>16)+c-cdf;//update
						//if(state<0x10000)//renorm
						//{
						//	state<<=16;
						//	if(srcptr-2>=srcstart)
						//	{
						//		srcptr-=2;
						//		memcpy(&state, srcptr, 2);
						//	}
						//}
					}
				}
			}
		}
		else
		{
			ACDecContext ctx;
			ac_dec_init(&ctx, CDF2, kc?data+bookmarks[kc-1]:data+overhead, data+bookmarks[kc]);

			for(int ky=0;ky<ih;++ky)//for each row
			{
				for(int bx=0;bx<gxcount;++bx)//for each group
				{
					int x1=bx*p->gwidth, x2=MINVAR(x1+p->gwidth, iw);
					int success=t26_prepblock(buf, CDF0+((size_t)kc<<8), iw, ih, kc, x1, x2, ky, p, CDF2, 0);
					if(!success)
						LOG_ERROR("t26_prepblock error");
					for(int kx=x1;kx<x2;++kx)//for each pixel
						buf[(iw*ky+kx)<<2|kc]=ac_dec(&ctx);
				}
			}
		}
	}
	free(CDF0);
	free(CDF2);

#ifndef T25_OPT_PRED
	apply_transforms_inv(buf, iw, ih);
#else
	addbuf(buf, iw, ih, 3, 4, 128);
	memcpy(jxlparams_i16,         predparams+predidx[0], predlen[0]*sizeof(short));
	memcpy(pw2_params+PW2_NPARAM, predparams+predidx[1], predlen[1]*sizeof(short));
	memcpy(jxlparams_i16+22,      predparams+predidx[2], predlen[2]*sizeof(short));
	pred_opt_apply((char*)buf, iw, ih, 0);
	//pred_jxl_apply((char*)buf, iw, ih, jxlparams, 0);
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);
#endif

	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
#ifdef AC_VALIDATE
	array_free(&acval);
#endif
	return 1;
}