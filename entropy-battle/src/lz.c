#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<math.h>
#include"lodepng.h"//for testing
static const char file[]=__FILE__;

static void print_hist(int *hist, int hstride, int nlevels)
{
	int vmax=0;
	for(int k=0;k<nlevels;++k)
	{
		if(vmax<hist[hstride*k])
			vmax=hist[hstride*k];
	}
	if(!vmax)
		return;
	for(int k=0;k<nlevels;++k)
	{
		printf("%3d %6d ", k, hist[hstride*k]);
		for(int k2=0, end=hist[hstride*k]*64/vmax;k2<end;++k2)
			printf("*");
		printf("\n");
	}
}
typedef struct PaletteStruct
{
	int sym, freq;
} Palette;
static int palette_cmp(const void *left, const void *right)
{
	Palette const *a, *b;

	a=(Palette const*)left;
	b=(Palette const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static void palette_add(ArrayHandle *palette, int val)
{
	size_t idx=0;
	Palette temp={val, 1};
	int found=binary_search(palette[0]->data, palette[0]->count, palette[0]->esize, palette_cmp, &temp, &idx);
	if(!found)
		array_insert(palette, idx, &temp, 1, 1, 0);
	else
	{
		Palette *t2=(Palette*)array_at(palette, idx);
		++t2->freq;
	}
}
static int palette_replace(ArrayHandle *palette, int val)
{
	size_t idx=0;
	Palette temp={val, 1};
	int found=binary_search(palette[0]->data, palette[0]->count, palette[0]->esize, palette_cmp, &temp, &idx);
	if(!found)
	{
		LOG_ERROR("val = %d not found in palette", val);
		return -1;
	}
	return (int)idx;
	//Palette *t2=(Palette*)array_at(palette, idx);
	//return t2->sym;
}

typedef struct LZInfoStruct
{
	int idx,//-1 means bypass
		len;
} LZInfo;
static ArrayHandle table[256];
//static int table[256];

typedef struct ContextInfoStruct
{
	unsigned sym, type;
} ContextInfo;

typedef struct FreqInfoStruct
{
	unsigned
		freq,	//original freq
		sym,    //symbol
		qfreq;  //quantized freq
} FreqInfo;
static int freqinfo_byfreq(const void *left, const void *right)
{
	FreqInfo const *a, *b;

	a=(FreqInfo const*)left;
	b=(FreqInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
static int freqinfo_bysym(const void *left, const void *right)
{
	FreqInfo const *a, *b;

	a=(FreqInfo const*)left;
	b=(FreqInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static void quantize_hist(const int *hist, int hstride, int nlevels, unsigned *CDF, int ostride)
{
	int count=0;
	size_t len=(size_t)nlevels*hstride;
	for(int k=0;k<len;k+=hstride)
		count+=hist[k];

	FreqInfo *buf=(FreqInfo*)malloc(nlevels*sizeof(FreqInfo));
	if(!buf)
	{
		LOG_ERROR("malloc error");
		return;
	}
#if 0
	double dsum=0, x=0;
#endif
	long long error=-0x100000000;
	for(int k=0;k<nlevels;++k)
	{
		int f=hist[hstride*k];
#if 0
		x=(double)f/count;
		dsum+=x;
#endif
		unsigned qf=(unsigned)(((unsigned long long)f<<32)/count);
		buf[k].freq=f;
		buf[k].sym=k;
		buf[k].qfreq=qf;
		error+=qf;
	}
	isort(buf, nlevels, sizeof(FreqInfo), freqinfo_byfreq);
	
	if(error>0)
	{
		while(error)
		{
			for(int k=0;k<nlevels&&error;++k)
			{
				int dec=buf[k].qfreq>1;
				buf[k].qfreq-=dec, error-=dec;
			}
		}
	}
	else
	{
		while(error)
		{
			for(int k=nlevels-1;k>=0&&error;--k)
			{
				int inc=buf[k].qfreq<0xFFFFFFFF;
				buf[k].qfreq+=inc, error+=inc;
			}
		}
	}
	isort(buf, nlevels, sizeof(FreqInfo), freqinfo_bysym);

	unsigned long long sum=0;
	for(int k=0;k<nlevels;++k)
	{
		unsigned f=buf[k].qfreq;
		CDF[ostride*k]=(unsigned)sum;
		sum+=f;
	}
}
long long lz_encode(const void *src, int bw, int bh, int bytestride, ArrayHandle *data, size_t *ret_overhead)//FIXME encoding channels separately
{
#if 1
	long long cycles=__rdtsc();
	
	size_t count=(size_t)bw*bh;
	unsigned char *buf=(unsigned char*)malloc(count);
	for(ptrdiff_t k=0;k<(ptrdiff_t)count;++k)
		buf[k]=((unsigned char*)src)[bytestride*k];
	for(int ky=0;ky<bh;++ky)//differentiate horizontally
	{
		for(int kx=bw-1;kx>0;--kx)
			buf[bw*ky+kx]-=buf[bw*ky+kx-1];
	}
	for(int ky=bh-1;ky>0;--ky)//differentiate vertically
	{
		for(int kx=0;kx<bw;++kx)
			buf[bw*ky+kx]-=buf[bw*(ky-1)+kx];
	}
	//const unsigned char *buf=(const unsigned char*)src;

	DList list;
	dlist_init(&list, sizeof(LZInfo), 1024, 0);

	for(int k=0;k<256;++k)
		ARRAY_ALLOC(ptrdiff_t, table[k], 0, 0, 0, 0);
	
	LZInfo temp={-1, 0};
	for(ptrdiff_t ks=0;ks<(ptrdiff_t)count;)
	{
		unsigned char sym=buf[ks];
		ArrayHandle *key=table+sym;
		ptrdiff_t seq_idx=-1, longest_match=0;
		for(ptrdiff_t k2=0;k2<(ptrdiff_t)key[0]->count;++k2)
		{
			ptrdiff_t idx=*(ptrdiff_t*)array_at(key, k2);
			ptrdiff_t km=0;
			for(ptrdiff_t end=count-ks;km<end&&buf[idx+km]==buf[ks+km];++km);
			if(longest_match<km)
			{
				longest_match=km;
				//seq_idx=k2;//sequence number
				seq_idx=idx;
			}
		}
		ARRAY_APPEND(*key, &ks, 1, 1, 0);
		if(seq_idx==-1)//bypass
		{
			if(temp.idx!=-1)
			{
				temp.idx=-1;
				temp.len=1;
			}
			else
				++temp.len;
			//int bypass_flag=-1;
			//dlist_push_back(&list, &bypass_flag, 1);//bypass: {0xFF, byte}
			//dlist_push_back(&list, &sym, 1);
			++ks;
		}
		else//repeat
		{
			if(temp.idx==-1)
				dlist_push_back1(&list, &temp);
			//temp.idx=(int)seq_idx;//sequence number
			temp.idx=(int)(ks-seq_idx);//backtrack
			temp.len=(int)longest_match;
			dlist_push_back1(&list, &temp);
			//dlist_push_back(&list, &seq_idx, 4);
			//dlist_push_back(&list, &longest_match, 4);
			ks+=longest_match;
		}
	}
	for(int k=0;k<256;++k)
		array_free(table+k);

	ArrayHandle t2=0;
	dlist_appendtoarray(&list, &t2);
	dlist_clear(&list);
	LZInfo *tags=(LZInfo*)t2->data;
	int ntags=(int)t2->count;

	ArrayHandle palette=0;
	ARRAY_ALLOC(Palette, palette, 0, 0, 0, 0);
	for(int k=0;k<ntags;++k)
	{
		LZInfo *tag=tags+k;
		palette_add(&palette, tag->idx);
		palette_add(&palette, tag->len);
	}

	//for(int k=0;k<(int)palette->count;++k)//
	//{
	//	Palette *pk=(Palette*)array_at(&palette, k);
	//	printf("%3d %5d %6d\n", k, pk->sym, pk->freq);
	//}

	//quantize histogram
	unsigned *palette_CDF=(unsigned*)malloc(palette->count*sizeof(unsigned));
	if(!palette_CDF)
	{
		LOG_ERROR("malloc error");
		return 0;
	}
	quantize_hist((int*)palette->data+1, 2, (int)palette->count, palette_CDF, 1);

	//push symbols
	dlist_init(&list, sizeof(ContextInfo), 1024, 0);
	ContextInfo ctx={0};
	for(ptrdiff_t kt=0, ks=0;kt<(ptrdiff_t)ntags;++kt)
	{
		ctx.sym=palette_replace(&palette, tags[kt].idx);
		ctx.type=1;
		dlist_push_back1(&list, &ctx);

		ctx.sym=palette_replace(&palette, tags[kt].len);
		ctx.type=1;
		dlist_push_back1(&list, &ctx);

		if(tags[kt].idx==-1)//bypass
		{
			ptrdiff_t end=ks+tags[kt].len;
			for(;ks<end;++ks)
			{
				ctx.sym=buf[ks];
				ctx.type=0;
				dlist_push_back1(&list, &ctx);
			}
		}
		else
			ks+=tags[kt].len;
	}
	array_free(&t2);
	tags=0, ntags=0;
	dlist_appendtoarray(&list, &t2);
	dlist_clear(&list);

	dlist_init(&list, 1, 1024, 0);
	//write header
	dlist_push_back(&list, &palette->count, 4);		//nlevels
	for(int k=0;k<(int)palette->count;++k)			//palette
	{
		Palette *p=(Palette*)array_at(&palette, k);
		dlist_push_back(&list, &p->sym, 4);
	}
	for(int k=0;k<(int)palette->count;++k)			//CDF
		dlist_push_back(&list, palette_CDF+k, 4);

	if(ret_overhead)
		*ret_overhead=list.nobj;
	//printf("LZO+rANS Overhead %lld\n", list.nobj);

	//encode symbols backwards
	unsigned long long state=0x100000000;
	for(int ks=(int)t2->count-1;ks>=0;--ks)
	{
		ContextInfo *sym=(ContextInfo*)array_at(&t2, ks);
		unsigned CDF, freq;
		if(sym->type)
		{
			CDF=palette_CDF[sym->sym];
			freq=(unsigned)((sym->sym+1<palette->count?(unsigned long long)palette_CDF[sym->sym+1]:0x100000000ull)-CDF);
		}
		else//bypass is uniform over 256 levels
		{
			CDF=sym->sym<<24;
			freq=1<<24;
		}
		if(state>=(unsigned long long)freq<<32)
		{
			dlist_push_back(&list, &state, 4);
			state>>=32;
		}
		state=state/freq<<32|(CDF+state%freq);
	}
	dlist_push_back(&list, &state, 8);

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);

	array_free(&palette);
	free(palette_CDF);
	free(buf);

	cycles=__rdtsc()-cycles;
	return cycles;
#endif
	
#if 0
	int vmax=0;
	for(int k=0;k<ntags;++k)
	{
		LZInfo *tag=tags+k;
		if(tag->idx!=-1&&vmax<tag->idx)
			vmax=tag->idx;
		if(vmax<tag->len)
			vmax=tag->len;
	}
	++vmax;
	int *hist=(int*)malloc(vmax*sizeof(int));
	if(!hist)
		return 0;
	memset(hist, 0, vmax*sizeof(int));
	for(int k=0;k<ntags;++k)
	{
		LZInfo *tag=tags+k;
		if(tag->idx!=-1)
			++hist[tag->idx];
		++hist[tag->len];
	}

	print_hist(hist, 1, vmax);//

	free(hist);
#endif

#if 0
	const unsigned char *buf=(const unsigned char*)src;

	long long cycles=__rdtsc();
	size_t csize=0;

	//int tablesize=symbytes*256;
	//Offset2d *table=(Offset2d*)malloc((size_t)tablesize*sizeof(Offset2d));
	//if(!table)
	//	return 0;
	//memset(table, -1, (size_t)tablesize*sizeof(Offset2d));

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int state[]={-1, 0};//{srcidx (-1 means bypass), copysize}
	memset(table, -1, (size_t)256*sizeof(int));
	for(int ks=0;ks<count;)
	{
		unsigned char sym=buf[bytestride*ks];
		int *entry=table+sym;

		if(*entry==-1)//first time seen, emit bypass
		{
			*entry=ks;
			++state[1];
			++csize;
			++ks;
		}
		else
		{
			dlist_push_back(&list, state, 2*sizeof(int));
			unsigned char s2=0;
			state[0]=*entry;
			state[1]=0;
			*entry=ks;//ambiguous part		FIXME should choose longest match
			do
			{
				s2=buf[bytestride*(state[0]+state[1])];
				sym=buf[bytestride*ks];
				++state[1];
				++ks;
			}
			while(s2==sym&&ks<count);
			dlist_push_back(&list, state, 2*sizeof(int));
			state[0]=-1;
			state[1]=0;
#if 0
		if(*entry==-1)//first time seen, emit bypass
		{
			*entry=ks;
			if(state[0]==-1)
			{
				++state[1];
				++csize;
			}
			else
			{
				dlist_push_back(&list, state, 2*sizeof(int));
				state[0]=-1;
				state[1]=1;
			}
			++ks;
		}
		else
		{
			if(state[0]==-1)
				dlist_push_back(&list, state, 2*sizeof(int));
			unsigned char s2=0;
			state[0]=*entry;
			state[1]=0;
			*entry=ks;//ambiguous part
			do
			{
				s2=buf[bytestride*(state[0]+state[1])];
				sym=buf[bytestride*ks];
				++state[1];
				++ks;
			}
			while(s2==sym&&ks<count);
			dlist_push_back(&list, state, 2*sizeof(int));
			state[0]=-1;
			state[1]=0;
			//unsigned char s2=buf[state[0]+state[1]];
			//if(sym==s2)//continue copying
			//	++state[1];
			//else
			//{
			//	dlist_push_back(&list, state, 2*sizeof(int));
			//
			//}
#endif
		}
	}
	//free(table);
	dlist_appendtoarray(&list, data);

	cycles=__rdtsc()-cycles;

	if(ret_csize)
		*ret_csize=csize;

	int *tags=(int*)data[0]->data;
	int ntags=(int)(data[0]->count/(2*sizeof(int)));
	int max_bypass=0, max_repeat=0;
	for(int k=0;k<ntags;++k)
	{
		int *tag=tags+((size_t)k<<1);
		if(tag[0]==-1)//bypass
		{
			if(max_bypass<tag[1])
				max_bypass=tag[1];
		}
		else
		{
			if(max_repeat<tag[1])
				max_repeat=tag[1];
		}
	}
	++max_bypass;
	++max_repeat;
	int *hist_bypass=(int*)malloc(max_bypass*sizeof(int)),
		*hist_repeat=(int*)malloc(max_repeat*sizeof(int));
	if(!hist_bypass||!hist_repeat)
		return 0;
	memset(hist_bypass, 0, max_bypass*sizeof(int));
	memset(hist_repeat, 0, max_repeat*sizeof(int));
	for(int k=0;k<ntags;++k)
	{
		int *tag=tags+((size_t)k<<1);
		if(tag[0]==-1)//bypass
			++hist_bypass[tag[1]];
		else
			++hist_repeat[tag[1]];
	}
	printf("bypass:\n");
	print_hist(hist_bypass, 1, max_bypass);
	printf("repeat:\n");
	print_hist(hist_repeat, 1, max_repeat);
	free(hist_bypass);
	free(hist_repeat);
	return cycles;
#endif
}

//LZ2 UNDECODABLE
#if 0
int lz2_limit=1;
static int history[65536];
void lz2_encode(unsigned char *buf, int len, ArrayHandle *coeff, ArrayHandle *bypass)
{
	DList list;
	dlist_init(&list, 1, 1024, 0);
	memset(history, -1, 65536*sizeof(int));
	ARRAY_ALLOC(char, *bypass, 0, 0, 0, 0);
	LZInfo temp={-1, 0};
	for(int ks=0;ks<len;)
	{
		unsigned char sym=buf[ks];
		int *key=history+((size_t)sym<<8);
		unsigned char seq_idx=0xFF;
		int longest_match=0;
		int k2=0;
		for(;k2<255&&key[k2]!=-1;++k2)
		{
			int idx=key[k2];
			int len2=0;
			for(int end=len-ks;len2<end&&len2<255&&buf[idx+len2]==buf[ks+len2];++len2);
			if(longest_match<len2)
			{
				longest_match=len2;
				seq_idx=k2;//sequence number
			}
		}

		//if(longest_match<=lz2_limit)
		//	seq_idx=0xFF;

		if(seq_idx!=0xFF&&seq_idx)//place idx at seq_idx to front
		{
			int temp=key[seq_idx];
			key[seq_idx]=key[0];
			key[0]=temp;
		}
		key[k2-(k2==255)]=ks;//append current position to key

		////prepend ks, shift-out key
		//int ins=ks;
		//for(int k2=0;k2<255&&key[k2]!=-1;++k2)
		//{
		//	int temp=key[k2];
		//	key[k2]=temp;
		//	ins=temp;
		//}

		if(seq_idx==0xFF)//emit bypass
		{
			if(temp.idx!=-1)
			{
				temp.idx=-1;
				temp.len=1;
			}
			else
			{
				++temp.len;
				if(temp.len==0xFF)//end bypass
				{
					dlist_push_back1(&list, &temp.idx);
					dlist_push_back1(&list, &temp.len);
					ARRAY_APPEND(*bypass, buf+ks-temp.len, temp.len, 1, 0);
					temp.len=0;
				}
			}
			++ks;
		}
		else//emit repeat
		{
			if(temp.idx==-1)//end bypass
			{
				dlist_push_back1(&list, &temp.idx);
				dlist_push_back1(&list, &temp.len);
				ARRAY_APPEND(*bypass, buf+ks-temp.len, temp.len, 1, 0);
			}
			temp.idx=0;
			dlist_push_back1(&list, &seq_idx);
			dlist_push_back1(&list, &longest_match);
			ks+=longest_match;
		}
	}
	dlist_appendtoarray(&list, coeff);
	dlist_clear(&list);
}
#endif
#if 0
void lz2_decode(const unsigned char *coeff, const unsigned char *bypass, int coefflen, int bypasslen, int dstlen, unsigned char *dst)
{
	memset(history, -1, 65536*sizeof(int));
	for(int kc=0, kb=0, kd=0;kd<dstlen;)
	{
		unsigned char seq_idx=coeff[kc], len, sym;
		++kc;

		len=coeff[kc];
		++kc;

		if(seq_idx==0xFF)
		{
			for(int k2=0;k2<len;++k2)
			{
				sym=bypass[kb];
				++kb;

				int *key=history+((size_t)sym<<8);
				key[0]=kd;

				dst[kd]=sym;
				++kd;
			}
		}
		else
		{
			for(int k2=0;k2<len;++k2)
			{
				sym=
			}
		}
	}
}
#endif
typedef struct HistInfoStruct
{
	int	sym,  //symbol
		freq, //original freq
		qfreq;//quantized freq
} HistInfo;
static int histinfo_byfreq(const void *left, const void *right)
{
	HistInfo const *a, *b;

	a=(HistInfo const*)left;
	b=(HistInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
static int histinfo_bysym(const void *left, const void *right)
{
	HistInfo const *a, *b;

	a=(HistInfo const*)left;
	b=(HistInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static void rans_calc_CDF(const unsigned char *buffer, int len, int stride, unsigned short *CDF)
{
	const int nlevels=256;
	if(!len||!stride)
	{
		LOG_ERROR("rans_calc_CDF(): Invalid args");
		return;
	}
	HistInfo hist[256];
	for(int k=0;k<nlevels;++k)
	{
		hist[k].sym=k;
		hist[k].freq=0;
	}
	for(int k=0;k<len;k+=stride)//this loop takes 73% of encode time
		++hist[buffer[k]].freq;

#if 0
	static int call=0;
	if(call==3)
		print_hist(&hist->freq, sizeof(HistInfo)/sizeof(int), 256);//
	++call;
#endif

	int count=len/stride;
	for(int k=0;k<nlevels;++k)
		hist[k].qfreq=((long long)hist[k].freq<<16)/count;
	
	if(count!=0x10000)
	{
		const int prob_max=0xFFFF;

		isort(hist, nlevels, sizeof(HistInfo), histinfo_byfreq);
		int idx=0;
		for(;idx<nlevels&&!hist[idx].freq;++idx);
		for(;idx<nlevels&&!hist[idx].qfreq;++idx)
			++hist[idx].qfreq;
		//for(idx=nlevels-1;idx>=0&&hist[idx].qfreq>=prob_max;--idx);
		//for(++idx;idx<nlevels;++idx)
		//	hist[idx].qfreq=prob_max;

		int error=-0x10000;//too much -> +ve error & vice versa
		for(int k=0;k<nlevels;++k)
			error+=hist[k].qfreq;
		if(error>0)
		{
			while(error)
			{
				for(int k=0;k<nlevels&&error;++k)
				{
					int dec=hist[k].qfreq>1;
					hist[k].qfreq-=dec, error-=dec;
				}
			}
		}
		else
		{
			while(error)
			{
				for(int k=nlevels-1;k>=0&&error;--k)
				{
					int inc=hist[k].qfreq<prob_max;
					hist[k].qfreq+=inc, error+=inc;
				}
			}
		}
		isort(hist, nlevels, sizeof(HistInfo), histinfo_bysym);
	}
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		CDF[k]=sum;
		sum+=hist[k].qfreq;
	}
}
void rans_bytes_encode(unsigned char *buf, int len, int symbytes, unsigned short *CDF, int nlevels, DList *list)
{
	unsigned state=0x10000;
	for(int k=len-1;k>=0;--k)
	{
		unsigned char sym=buf[k];
		unsigned short c, freq;
		if(CDF)
		{
			int idx=(k%symbytes<<8)+sym;
			c=CDF[idx], freq=(sym<255?CDF[idx+1]:0x10000)-c;
		}
		else
			c=(sym<<16)/nlevels, freq=0x10000/nlevels;

		if(state>=(unsigned)(freq<<16))//renorm
		{
			dlist_push_back(list, &state, 2);
			state>>=16;
		}

		state=state/freq<<16|(c+state%freq);//update
	}
	dlist_push_back(list, &state, 4);
}
int rans_bytes_decode(const unsigned char *data, int srclen, int count, int symbytes, int nlevels, unsigned short *CDF, unsigned char *CDF2sym, unsigned char *dst)
{
	const unsigned char *srcptr=data+srclen;
	unsigned state;
	srcptr-=4;
	if(srcptr<data)
	{
		LOG_ERROR("rans_decode_bytes: idx %lld", srcptr-data);
		return 0;
	}
	memcpy(&state, srcptr, 4);
	for(int k=0;k<count;++k)
	{
		unsigned short c=(unsigned short)state;
		unsigned char sym;
		unsigned short cdf, freq;
		if(CDF)				//fetch
		{
			int kc=k%symbytes;
			sym=CDF2sym[kc<<16|c];
			cdf=CDF[kc<<8|sym], freq=(sym<255?CDF[kc<<8|(sym+1)]:0x10000)-cdf;
		}
		else
		{
			sym=c*nlevels>>16;
			cdf=(sym<<16)/nlevels, freq=0x10000/nlevels;
		}
		dst[k]=sym;

		state=freq*(state>>16)+c-cdf;//update

		if(state<0x10000)//renorm
		{
			srcptr-=2;
			if(srcptr<data)
			{
				LOG_ERROR("rans_decode_bytes: idx %lld", srcptr-data);
				return 0;
			}
			state<<=16;
			memcpy(&state, srcptr, 2);
		}
	}
	return 1;
}
double calc_sdev(unsigned char *buf, int len, int stride)
{
	int sum=0;
	for(int k=0;k<len;k+=stride)
		sum+=buf[k];
	double mean=(double)sum/len, s2=0;
	for(int k=0;k<len;k+=stride)
	{
		double x=buf[k]-mean;
		s2+=x*x;
	}
	s2=sqrt(s2/len);
	return s2;
}
static const int tag_lz02='L'|'Z'<<8|'0'<<16|'2'<<24;
#if 0
long long test1_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data)
{
	long long cycles=__rdtsc();
	
	unsigned short CDF[256];
	DList list;
	dlist_init(&list, 1, 1024, 0);
	const unsigned char *buf=(const unsigned char*)src;
	unsigned char *testrow=(unsigned char*)malloc((size_t)bw*5);
	unsigned char *b_type=(unsigned char*)malloc(bh);
	size_t res=(size_t)bw*bh;
	unsigned char *b2=(unsigned char*)malloc(res);
	for(int kc=symbytes-1;kc>=0;--kc)
	{
		for(int ky=bh-1;ky>=0;--ky)
		{
			int type=0;
			double sdev, temp;
			for(int kx=0;kx<bw;++kx)//bypass
				testrow[kx]=buf[bytestride*(bw*ky+kx)+kc];
			sdev=calc_sdev(testrow, bw, 1);

			for(int kx=0;kx<bw;++kx)//h-diff
				testrow[bw+kx]=buf[bytestride*(bw*ky+kx)+kc]-(kx?buf[bytestride*(bw*ky+kx-1)+kc]:0);
			temp=calc_sdev(testrow+bw, bw, 1);
			if(sdev>temp)
				sdev=temp, type=1;

			for(int kx=0;kx<bw;++kx)//v-diff
				testrow[bw*2+kx]=buf[bytestride*(bw*ky+kx)+kc]-(ky?buf[bytestride*(bw*(ky-1)+kx)+kc]:0);
			temp=calc_sdev(testrow+bw*2, bw, 1);
			if(sdev>temp)
				sdev=temp, type=2;

			for(int kx=0;kx<bw;++kx)//av-diff
				testrow[bw*3+kx]=buf[bytestride*(bw*ky+kx)+kc]-(((kx?buf[bytestride*(bw*ky+kx-1)+kc]:0)+(ky?buf[bytestride*(bw*(ky-1)+kx)+kc]:0))>>1);
			temp=calc_sdev(testrow+bw*3, bw, 1);
			if(sdev>temp)
				sdev=temp, type=3;

			for(int kx=0;kx<bw;++kx)//Paeth predictor
			{
				unsigned char
					A=kx?buf[bytestride*(bw*ky+kx-1)+kc]:0,
					B=ky?buf[bytestride*(bw*(ky-1)+kx)+kc]:0,
					C=kx&&ky?buf[bytestride*(bw*(ky-1)+kx-1)+kc]:0,
					p=A+B-C,
					dist, d2, sub;

				dist=abs(p-A);
				sub=A;

				d2=abs(p-B);
				if(dist>d2)
					dist=d2, sub=B;

				d2=abs(p-C);
				if(dist>d2)
					dist=d2, sub=C;

				testrow[bw*4+kx]=buf[bytestride*(bw*ky+kx)+kc]-sub;
			}
			temp=calc_sdev(testrow+bw*4, bw, 1);
			if(sdev>temp)
				sdev=temp, type=4;

			b_type[ky]=type;
			memcpy(b2+bw*ky, testrow+bw*type, bw);
		}
		ArrayHandle coeff=0, bypass=0;
		lz2_encode(b2, (int)res, &coeff, &bypass);

		rans_calc_CDF(coeff->data, (int)coeff->count, 1, CDF);
		dlist_push_back(&list, &tag_lz02, 4);
		dlist_push_back(&list, &coeff->count, 4);
		dlist_push_back(&list, CDF, 256*sizeof(short));
		rans_bytes_encode(coeff->data, (int)coeff->count, 1, CDF, 256, &list);
		dlist_push_back(&list, bypass->data, bypass->count);
		dlist_push_back(&list, b_type, bh);
		//rans_bytes_encode(bypass->data, (int)bypass->count, 0, 256, &list);
		//rans_bytes_encode(b_type, bh, 0, 5, &list);

		array_free(&coeff);
		array_free(&bypass);
	}
	dlist_appendtoarray(&list, data);
	dlist_clear(&list);

	cycles=__rdtsc()-cycles;
	return cycles;
}
#endif

static const int tag_lz03='L'|'Z'<<8|'0'<<16|'3'<<24;
static unsigned short g_CDF[256*4];
static unsigned char g_CDF2sym[65536*4];
#if 0
long long test2_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data)
{
	if(symbytes>4)
		return 0;
	long long cycles=__rdtsc();
	
	DList list;
	dlist_init(&list, 1, 1024, 0);
	const unsigned char *buf=(const unsigned char*)src;
	int rowlen=symbytes*bw, srcrowlen=bytestride*bw;
	unsigned char *testrow=(unsigned char*)malloc((size_t)rowlen*5);
	unsigned char *b_type=(unsigned char*)malloc(bh);
	size_t res=(size_t)rowlen*bh;
	unsigned char *b2=(unsigned char*)malloc(res);
	for(int ky=bh-1;ky>=0;--ky)
	{
		int type;
		double sdev, temp;

		for(int kx=0;kx<bw;++kx)//bypass
		{
			int srcidx=bytestride*(bw*ky+kx), dstidx=symbytes*kx;
			for(int kc=0;kc<symbytes;++kc)
				testrow[dstidx+kc]=buf[srcidx+kc];
		}
		sdev=0;
		for(int kc=0;kc<symbytes;++kc)
			sdev+=calc_sdev(testrow+kc, rowlen, symbytes);
		type=0;
		//omit division of sdev by symbytes

		for(int kx=0;kx<bw;++kx)//h-diff
		{
			int srcidx=bytestride*(bw*ky+kx), dstidx=rowlen+symbytes*kx;
			for(int kc=0;kc<symbytes;++kc)
				testrow[dstidx+kc]=buf[srcidx+kc]-(kx?buf[srcidx-bytestride+kc]:0);
		}
		temp=0;
		for(int kc=0;kc<symbytes;++kc)
			temp+=calc_sdev(testrow+rowlen+kc, rowlen, symbytes);
		if(sdev>temp)
			sdev=temp, type=1;

		for(int kx=0;kx<bw;++kx)//v-diff
		{
			int srcidx=bytestride*(bw*ky+kx), dstidx=rowlen*2+symbytes*kx;
			for(int kc=0;kc<symbytes;++kc)
				testrow[dstidx+kc]=buf[srcidx]-(ky?buf[srcidx-srcrowlen+kc]:0);
		}
		temp=0;
		for(int kc=0;kc<symbytes;++kc)
			temp+=calc_sdev(testrow+rowlen*2+kc, rowlen, symbytes);
		if(sdev>temp)
			sdev=temp, type=2;

		for(int kx=0;kx<bw;++kx)//av-diff
		{
			int srcidx=bytestride*(bw*ky+kx), dstidx=rowlen*3+symbytes*kx;
			for(int kc=0;kc<symbytes;++kc)
				testrow[dstidx+kc]=buf[srcidx+kc]-(((kx?buf[srcidx-bytestride+kc]:0)+(ky?buf[srcidx-srcrowlen+kc]:0))>>1);
		}
		temp=0;
		for(int kc=0;kc<symbytes;++kc)
			temp+=calc_sdev(testrow+rowlen*3+kc, rowlen, symbytes);
		if(sdev>temp)
			sdev=temp, type=3;

		for(int kx=0;kx<bw;++kx)//Paeth predictor
		{
			int srcidx=bytestride*(bw*ky+kx), dstidx=rowlen*4+symbytes*kx;
			for(int kc=0;kc<symbytes;++kc)
			{
				unsigned char
					A=kx?buf[srcidx-bytestride+kc]:0,
					B=ky?buf[srcidx-srcrowlen+kc]:0,
					C=kx&&ky?buf[srcidx-srcrowlen-bytestride+kc]:0,
					p=A+B-C,
					dist, d2, sub;

				dist=abs(p-A);
				sub=A;

				d2=abs(p-B);
				if(dist>d2)
					dist=d2, sub=B;

				d2=abs(p-C);
				if(dist>d2)
					dist=d2, sub=C;

				testrow[dstidx+kc]=buf[srcidx+kc]-sub;
			}
		}
		temp=0;
		for(int kc=0;kc<symbytes;++kc)
			temp+=calc_sdev(testrow+rowlen*4+kc, rowlen, symbytes);
		if(sdev>temp)
			sdev=temp, type=4;

		//printf("%3d %d\n", ky, type);//natural: 0		synthetic: variable

		b_type[ky]=type;
		memcpy(b2+rowlen*ky, testrow+rowlen*type, rowlen);
	}
	ArrayHandle coeff=0, bypass=0;
	lz2_encode(b2, (int)res, &coeff, &bypass);

	printf("coeff %lld bypass %lld\n", coeff->count, bypass->count);

	//lodepng_encode_file("out.PNG", b2, bw, bh, bytestride==4?LCT_RGBA:LCT_RGB, 8);//DEBUG SAVE

	for(int kc=0;kc<symbytes;++kc)
		rans_calc_CDF(coeff->data+kc, (int)coeff->count, symbytes, g_CDF+((size_t)kc<<8));
	dlist_push_back(&list, &tag_lz03, 4);
	dlist_push_back(&list, &coeff->count, 4);
	dlist_push_back(&list, &bypass->count, 4);
	dlist_push_back(&list, g_CDF, ((size_t)symbytes<<8)*sizeof(short));
	rans_bytes_encode(coeff->data, (int)coeff->count, symbytes, g_CDF, 256, &list);
	dlist_push_back(&list, bypass->data, bypass->count);
	//rans_bytes_encode(bypass->data, (int)bypass->count, 0, 256, &list);
	//dlist_push_back(&list, b_type, bh);
	rans_bytes_encode(b_type, bh, 1, 0, 5, &list);

	array_free(&coeff);
	array_free(&bypass);

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);

	cycles=__rdtsc()-cycles;
	return cycles;
}
#endif
#if 0
long long test2_decode(const void *src, size_t srclen, int bw, int bh, int symbytes, int bytestride, void *dst)
{
	long long cycles=__rdtsc();
	const unsigned char
		*data=(const unsigned char*)src,
		*srcptr=data,
		*srcend=data+srclen;
	ptrdiff_t coeffcount=0, bypasscount=0;

	if(srcptr+4>srcend||memcmp(srcptr, &tag_lz03, 4))
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}
	srcptr+=4;

	if(srcptr+4>srcend)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}
	memcpy(&coeffcount, srcptr, 4);
	srcptr+=4;

	if(srcptr+4>srcend)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}
	memcpy(&bypasscount, srcptr, 4);
	srcptr+=4;

	if(srcptr+((size_t)symbytes<<8)*sizeof(short)>srcend)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}
	memcpy(g_CDF, srcptr, ((size_t)symbytes<<8)*sizeof(short));
	srcptr+=((size_t)symbytes<<8)*sizeof(short);

	if(srcptr+coeffcount>srcend)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}

	for(int kc=0;kc<symbytes;++kc)
	{
		for(int sym=0;sym<256;++sym)
		{
			for(int k2=g_CDF[sym], end=(sym<255?g_CDF[sym+1]:0x10000);k2<end;++k2)
				g_CDF2sym[kc<<16|k2]=sym;
		}
	}
	unsigned char *coeff=(unsigned char*)malloc(coeffcount);
	if(!coeff)
	{
		LOG_ERROR("malloc fail");
		return 0;
	}
	rans_bytes_decode(srcptr, (int)(srcend-srcptr), (int)coeffcount, symbytes, 256, g_CDF, g_CDF2sym, coeff);//FIXME need to store rANS buffer length

	cycles=__rdtsc()-cycles;
	return cycles;
}
#endif


//LZ2D
typedef struct LZ2DSrcInfoStruct
{
	short srcx, srcy;
} LZ2DSrcInfo;
typedef unsigned char DeltaType;
//typedef unsigned short DeltaType;
typedef struct LZ2DInfoStruct//10 bytes
{
	unsigned short srcx, srcy;
	DeltaType w, h;
	unsigned short dstx, dsty;
} LZ2DInfo;
typedef struct LZ2DRLEInfoStruct//8 bytes
{
	unsigned short x, y, w, h;
} LZ2DRLEInfo;
LZ2DSrcInfo g_hist[256*16];
#if 0
int strided_valcmp(const unsigned char *dst, const unsigned char *src, int symbytes, int bytestride, int pxstride, int count, int pad1, int pad2, const unsigned char *dstmask)
{
	int fullstride=bytestride*pxstride;
	const unsigned char *end=dst+(size_t)fullstride*count;
	if(pad1&&!dstmask[-pxstride]&&!memcmp(dst-fullstride, src, symbytes))//check if entering a sea
		return 0;
	for(;dst<end;dst+=fullstride, dstmask+=pxstride)
	{
		if(*dstmask||memcmp(dst, src, symbytes))
			return 0;
	}
	if(pad2&&!*dstmask&&!memcmp(dst, src, symbytes))//check if entering a sea
		return 0;
	return 1;
}
#endif
int strided_memcmp(const unsigned char *dst, const unsigned char *src, int symbytes, int bytestride, int pxstride, int count, const unsigned char *dstmask, const unsigned char *srcmask)
{
	int fullstride=bytestride*pxstride;
	const unsigned char *end=src+(size_t)fullstride*count;
	for(;src<end;src+=fullstride, dst+=fullstride, srcmask+=pxstride, dstmask+=pxstride)
	{
		if(*srcmask||*dstmask||memcmp(dst, src, symbytes))
			return 0;
		//for(int k=0;k<symbytes;++k)
		//{
		//	if(src[k]!=dst[k])
		//		return 0;
		//}
	}
	return 1;
}
#if 0
size_t lz2d_encode(const unsigned char *buf, int bw, int bh, int symbytes, int bytestride, ArrayHandle *mask, ArrayHandle *coeff)
{
	size_t res=(size_t)bw*bh;
	ARRAY_ALLOC(char, *mask, 0, res, 0, 0);
	if(!*mask)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	memset(mask[0]->data, 0, res);
	memset(g_hist, -1, sizeof(g_hist));

	ARRAY_ALLOC(LZ2DInfo, *coeff, 0, 0, 0, 0);
	if(!*coeff)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}

	int maxw=0, maxh=0, maxaw=0, maxah=0;//

	size_t savedbytes=0;
	for(int ky=0;ky<bh;++ky)
	{
#if 0
		if(!(ky&15))//
			printf("\r%d / %d = %lf", ky+1, bh, 100.*(ky+1)/bh);//
#endif
		for(int kx=0;kx<bw;)
		{
			LZ2DSrcInfo *p;
			int idx=bw*ky+kx;
			if(mask[0]->data[idx])
			{
				++kx;
				continue;
			}
			unsigned char sym=buf[bytestride*idx], bestmatch=0xFF;
			int bestw=0, besth=0;
			int kh=0;
			for(;kh<16;++kh)
			{
				p=g_hist+(sym<<4|kh);
				if(p->srcx==-1)
					break;
				int blockw=1, blockh=1, updatedw=1, updatedh=1;
				do
				{
					int srcidx, dstidx;
					if(updatedh&&blockh<(1<<(sizeof(DeltaType)<<3))&&blockh<blockw*4&&ky+blockh<bh&&p->srcy+blockh<bh)
					{
						srcidx=bw*(p->srcy+blockh)+p->srcx, dstidx=bw*(ky+blockh)+kx;
						if(strided_memcmp(buf+bytestride*dstidx, buf+bytestride*srcidx, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx, mask[0]->data+srcidx))
							++blockh;
						else
							updatedh=0;
					}
					else
						updatedh=0;
					if(updatedw&&blockw<(1<<(sizeof(DeltaType)<<3))&&blockw<blockh*4&&kx+blockw<bw&&p->srcx+blockw<bw)
					{
						srcidx=bw*p->srcy+p->srcx+blockw, dstidx=bw*ky+kx+blockw;
						if(strided_memcmp(buf+bytestride*dstidx, buf+bytestride*srcidx, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx, mask[0]->data+srcidx))
							++blockw;
						else
							updatedw=0;
					}
					else
						updatedw=0;
				}while(updatedh||updatedw);
				if(bestw*besth<blockw*blockh)
					bestw=blockw, besth=blockh, bestmatch=kh;
			}

			if(maxw<bestw)
				maxw=bestw;
			if(maxh<besth)
				maxh=besth;
			if(maxaw*maxah<bestw*besth)
				maxaw=bestw, maxah=besth;

			int redundantbytes=symbytes*bestw*besth, success=bestmatch!=0xFF&&redundantbytes>sizeof(LZ2DInfo);
			if(success)//ignore if repeated bytes <= emission bytes
			{
				LZ2DInfo *emit=(LZ2DInfo*)ARRAY_APPEND(*coeff, 0, 1, 1, 0);//emit repeated block
				p=g_hist+(sym<<4|bestmatch);
				emit->srcx=p->srcx;
				emit->srcy=p->srcy;
				emit->w=bestw-1;
				emit->h=besth-1;
				emit->dstx=kx;
				emit->dsty=ky;

				unsigned char val=1+rand()%255;
				for(int ky2=0;ky2<besth;++ky2)//mask block
					//memset(mask[0]->data+bw*(ky+ky2)+kx, 0xFF, bestw);
				{
					for(int kx2=0;kx2<bestw;++kx2)
					{
						unsigned char *p3=mask[0]->data+bw*(ky+ky2)+kx+kx2;
						*p3=val;
						//*p3+=*p3<255;
					}
				}

				savedbytes+=redundantbytes;
			}

			//update g_hist
			if(kh>15)
				kh=15;
			p=g_hist+((size_t)sym<<4|kh);
			p->srcx=kx;
			p->srcy=ky;
			if(kh)
			{
				LZ2DSrcInfo temp, *p2=g_hist+((size_t)sym<<4);
				SWAPVAR(*p, *p2, temp);
			}

			if(success)
				kx+=bestw;
			else
				++kx;
		}
	}
#if 0
	printf("\n");//
#endif
	printf("maxw %d maxh %d maxA %d*%d\n", maxw, maxh, maxaw, maxah);
#if 0
	if(savedbytes*256<res*symbytes)//if saved bytes < (1/256) image bytes
	{
		array_free(mask);
		array_free(coeff);
		return 0;
	}
#endif
	return savedbytes;
}
size_t lz2d2_encode(const unsigned char *buf, int bw, int bh, int symbytes, int bytestride, ArrayHandle *mask, ArrayHandle *rle, ArrayHandle *lz)
{
	size_t res=(size_t)bw*bh;
	ARRAY_ALLOC(char, *mask, 0, res, 0, 0);
	if(!*mask)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	memset(mask[0]->data, 0, res);
	
	memset(g_hist, -1, sizeof(g_hist));

	ARRAY_ALLOC(LZ2DRLEInfo, *rle, 0, 0, 0, 0);
	if(!*rle)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}

	ARRAY_ALLOC(LZ2DInfo, *lz, 0, 0, 0, 0);
	if(!*lz)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	
	int maxw=0, maxh=0, maxaw=0, maxah=0;//

	size_t savedbytes=0;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;)
		{
			//if(kx==1&&ky==54)//
			//	kx=1;//

			int idx=bw*ky+kx;
			if(mask[0]->data[idx])
			{
				++kx;
				continue;
			}
			int blockw=1, blockh=1, updatew=1, updateh=1;
			do
			{
				int dstidx;
				if(updatew)
				{
					if(blockw<(1<<(sizeof(DeltaType)<<3))&&kx+blockw<bw)//try to add next column
					{
						dstidx=bw*ky+kx+blockw;
						if(strided_valcmp(buf+bytestride*dstidx, buf+bytestride*idx, symbytes, bytestride, bw, blockh, 0, !updateh&&ky+blockh<bh, mask[0]->data+dstidx))
							++blockw;
						else
							updatew=0;
					}
					else
						updatew=0;
				}
				if(updateh)
				{
					if(blockh<(1<<(sizeof(DeltaType)<<3))&&ky+blockh<bh)//try to add next row
					{
						dstidx=bw*(ky+blockh)+kx;
						if(strided_valcmp(buf+bytestride*dstidx, buf+bytestride*idx, symbytes, bytestride, 1, blockw, kx>0, !updatew&&kx+blockw<bw, mask[0]->data+dstidx))
							++blockh;
						else
							updateh=0;
					}
					else
						updateh=0;
				}
			}while(updateh||updatew);
			
			if(maxw<blockw)
				maxw=blockw;
			if(maxh<blockh)
				maxh=blockh;
			if(maxaw*maxah<blockw*blockh)
				maxaw=blockw, maxah=blockh;

			size_t redundant=((size_t)blockw*blockh-1)*symbytes;
			if(redundant>sizeof(LZ2DRLEInfo))
			{
				LZ2DRLEInfo *emit=(LZ2DRLEInfo*)ARRAY_APPEND(*rle, 0, 1, 1, 0);//emit RLE block
				emit->x=kx;
				emit->y=ky;
				emit->w=blockw;
				emit->h=blockh;
				
				//mask out redundant pixels
				unsigned char val=64+rand()%31;
				if(blockw>1)
					memset(mask[0]->data+bw*ky+kx+1, val, blockw-1);//leave first pixel unmasked
				for(int ky2=1;ky2<blockh;++ky2)
					memset(mask[0]->data+bw*(ky+ky2)+kx, val, blockw);

				savedbytes+=redundant;
				kx+=blockw;
			}
			else
				++kx;
		}
	}
	printf("RLE maxw %d maxh %d maxA %d*%d\n", maxw, maxh, maxaw, maxah);

	maxw=0, maxh=0, maxaw=0, maxah=0;

	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;)
		{
			LZ2DSrcInfo *p;
			int idx=bw*ky+kx;
			if(mask[0]->data[idx])
			{
				++kx;
				continue;
			}
			unsigned char sym=buf[bytestride*idx], bestmatch=0xFF;
			int bestw=0, besth=0;
			int kh=0;
			for(;kh<16;++kh)
			{
				p=g_hist+(sym<<4|kh);
				if(p->srcx==-1)
					break;
				int blockw=1, blockh=1, updatedw=1, updatedh=1;
				do
				{
					int srcidx, dstidx;
					if(updatedh&&blockh<(1<<(sizeof(DeltaType)<<3))&&blockh<blockw*4&&ky+blockh<bh&&p->srcy+blockh<bh)
					{
						srcidx=bw*(p->srcy+blockh)+p->srcx, dstidx=bw*(ky+blockh)+kx;
						if(strided_memcmp(buf+bytestride*dstidx, buf+bytestride*srcidx, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx, mask[0]->data+srcidx))
							++blockh;
						else
							updatedh=0;
					}
					else
						updatedh=0;
					if(updatedw&&blockw<(1<<(sizeof(DeltaType)<<3))&&blockw<blockh*4&&kx+blockw<bw&&p->srcx+blockw<bw)
					{
						srcidx=bw*p->srcy+p->srcx+blockw, dstidx=bw*ky+kx+blockw;
						if(strided_memcmp(buf+bytestride*dstidx, buf+bytestride*srcidx, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx, mask[0]->data+srcidx))
							++blockw;
						else
							updatedw=0;
					}
					else
						updatedw=0;
				}while(updatedh||updatedw);
				if(bestw*besth<blockw*blockh)
					bestw=blockw, besth=blockh, bestmatch=kh;
			}

			if(maxw<bestw)
				maxw=bestw;
			if(maxh<besth)
				maxh=besth;
			if(maxaw*maxah<bestw*besth)
				maxaw=bestw, maxah=besth;

			int redundantbytes=symbytes*bestw*besth, success=bestmatch!=0xFF&&redundantbytes>sizeof(LZ2DInfo);
			if(success)//ignore if repeated bytes <= emission bytes
			{
				LZ2DInfo *emit=(LZ2DInfo*)ARRAY_APPEND(*lz, 0, 1, 1, 0);//emit repeated block
				p=g_hist+(sym<<4|bestmatch);
				emit->srcx=p->srcx;
				emit->srcy=p->srcy;
				emit->w=bestw-1;
				emit->h=besth-1;
				emit->dstx=kx;
				emit->dsty=ky;

				unsigned char val=224+rand()%31;
				for(int ky2=0;ky2<besth;++ky2)//mask block
					//memset(mask[0]->data+bw*(ky+ky2)+kx, 0xFF, bestw);
				{
					for(int kx2=0;kx2<bestw;++kx2)
					{
						unsigned char *p3=mask[0]->data+bw*(ky+ky2)+kx+kx2;
						*p3=val;
						//*p3+=*p3<255;
					}
				}

				savedbytes+=redundantbytes;
			}

			//update g_hist
			if(kh>15)
				kh=15;
			p=g_hist+((size_t)sym<<4|kh);
			p->srcx=kx;
			p->srcy=ky;
			if(kh)
			{
				LZ2DSrcInfo temp, *p2=g_hist+((size_t)sym<<4);
				SWAPVAR(*p, *p2, temp);
			}

			if(success)
				kx+=bestw;
			else
				++kx;
		}
	}
	printf("LZ  maxw %d maxh %d maxA %d*%d\n", maxw, maxh, maxaw, maxah);
	return savedbytes;
}
#endif

typedef struct CNCPointStruct
{
	short x, y;
} CNCPoint;
size_t lz2d3_encode(const unsigned char *buf, int bw, int bh, int symbytes, int bytestride, ArrayHandle *mask, ArrayHandle *rle, ArrayHandle *lz)
{
	ARRAY_ALLOC(LZ2DRLEInfo, *rle, 0, 0, 0, 0);
	if(!*rle)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	
	ARRAY_ALLOC(LZ2DInfo, *lz, 0, 0, 0, 0);
	if(!*lz)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	
	size_t res=(size_t)bw*bh;
	ARRAY_ALLOC(char, *mask, 0, res, 0, 0);
	if(!*mask)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	memset(mask[0]->data, 1, res);//lgstep+1	initially all RLE blocks are 1x1

	CNCPoint *temp=(CNCPoint*)malloc(res*sizeof(CNCPoint));
	if(!temp)
	{
		LOG_ERROR("lz2d_encode(): Allocation failed");
		return 0;
	}
	temp->x=1;
	temp->y=1;
	memfill(temp+1, temp, (res-1)*sizeof(CNCPoint), sizeof(CNCPoint));

	size_t savedbytes=0;
#if 1
	for(int lgstep=1;lgstep<=8;++lgstep)
	{
		int step=1<<lgstep;
		int halfx=step>>1, halfy=bw*halfx;
		for(int ky=0, yend=bh-halfx;ky<yend;ky+=step)
		{
			for(int kx=0, xend=bw-halfx;kx<xend;kx+=step)
			{
				int idx=bw*ky+kx;
				if(temp[idx].x==halfx&&temp[idx+halfx].x==halfx&&temp[idx+halfy].x==halfx&&temp[idx+halfy+halfx].x==halfx)
				{
					int r1=memcmp(buf+bytestride*idx, buf+bytestride*(idx+halfx), symbytes),
						r2=memcmp(buf+bytestride*idx, buf+bytestride*(idx+halfy), symbytes),
						r3=memcmp(buf+bytestride*idx, buf+bytestride*(idx+halfy+halfx), symbytes);
					if(!r1&&!r2&&!r3)
					{
						temp[idx].x=step;
						temp[idx].y=step;

						for(int ky2=0;ky2<temp[idx].y;++ky2)
							memset(mask[0]->data+bw*(ky+ky2)+kx, lgstep+1, temp[idx].x);
					}
				}
			}
		}
	}
#if 1
	for(int ky=0;ky<bh;++ky)//merge adjacent blocks with same height
	{
		for(int kx=0;kx<bw;)
		{
			int idx=bw*ky+kx;
			if(mask[0]->data[idx]==0x10)
			{
				++kx;
				continue;
			}
			CNCPoint *p=temp+idx, *p2=temp+idx+p->x;
			//while there is room to the right && p2 is responsible for its block && same height && same content && no overflow
			while(kx+p->x<bw&&mask[0]->data[idx+p->x]&15&&p->y==p2->y&&!memcmp(buf+bytestride*idx, buf+bytestride*(idx+p->x), symbytes)&&p->x+p2->x<0x8000)
			{
				p->x+=p2->x;
				p2=temp+idx+p->x;
			}

			//mask out merged block
			for(int ky2=0;ky2<temp[idx].y;++ky2)
				memset(mask[0]->data+bw*(ky+ky2)+kx, 0x10, temp[idx].x);
			mask[0]->data[idx]|=floor_log2(p->y)+1;//0~7 + 1

			kx+=p->x;
		}
	}
#endif
#if 1
	for(int ky=0;ky<bh;++ky)//merge stacked blocks with same width
	{
		for(int kx=0;kx<bw;)
		{
			int idx=bw*ky+kx;
			if(mask[0]->data[idx]==0x20)
			{
				++kx;
				continue;
			}
			CNCPoint *p=temp+idx, *p2=temp+idx+bw*p->y;
			//while there is room downwards && p2 is responsible for its block && same height && same content && no overflow
			while(ky+p->y<bh&&mask[0]->data[idx+bw*p->y]&15&&p->x==p2->x&&!memcmp(buf+bytestride*idx, buf+bytestride*(idx+bw*p->y), symbytes)&&p->y+p2->y<0x8000)
			{
				p->y+=p2->y;
				p2=temp+idx+bw*p->y;
			}

			//mask out merged block
			for(int ky2=0;ky2<temp[idx].y;++ky2)
				memset(mask[0]->data+bw*(ky+ky2)+kx, 0x20, temp[idx].x);
			mask[0]->data[idx]|=floor_log2(p->y)+1;//0~14 + 1

			kx+=p->x;
		}
	}
#endif
#if 1
	for(int ky=0;ky<bh;++ky)//emit RLE codes
	{
		for(int kx=0;kx<bw;)
		{
			//if(kx==1910&&ky==1074)//
			//	kx=1910;//
			//if(kx==1908&&ky==1076)//
			//	kx=1908;//

			int idx=bw*ky+kx;
			if(!mask[0]->data[idx]||mask[0]->data[idx]>=0xC0)
			{
				++kx;
				continue;
			}
			size_t area=symbytes*((size_t)temp[idx].x*temp[idx].y-1);
		//	if(area>minlzbytes)
			if(area>sizeof(LZ2DRLEInfo))
			{
				//mask out block
				unsigned char val=0xC0+rand()%63;//192~255, for debugging
				mask[0]->data[idx]=0;//unmask top-left pixel
				if(temp[idx].x>1)
					memset(mask[0]->data+bw*ky+kx+1, val, temp[idx].x-1);
				for(int ky2=1;ky2<temp[idx].y;++ky2)
					memset(mask[0]->data+bw*(ky+ky2)+kx, val, temp[idx].x);

				LZ2DRLEInfo *emit=(LZ2DRLEInfo*)ARRAY_APPEND(*rle, 0, 1, 1, 0);
				emit->x=kx;
				emit->y=ky;
				emit->w=temp[idx].x;
				emit->h=temp[idx].y;

				savedbytes+=area;
			}
			else
			{
				for(int ky2=0;ky2<temp[idx].y;++ky2)
					memset(mask[0]->data+bw*(ky+ky2)+kx, 0, temp[idx].x);
			}
			kx+=temp[idx].x;
		}
	}
#endif
#else//no RLE
	memset(mask[0]->data, 0, res);
#endif

	//LZ2D
#if 1
	memset(g_hist, -1, sizeof(g_hist));
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;)
		{
			//if(kx==1704&&ky==1073)//
			//	kx=1704;//

			LZ2DSrcInfo *p;
			int idx=bw*ky+kx;
			if(mask[0]->data[idx])
			{
				++kx;
				continue;
			}
			unsigned char sym=buf[bytestride*idx], bestmatch=0xFF;
			int bestw=0, besth=0;
			int kh=0;
			for(;kh<16;++kh)
			{
				//if(kx==1182&&ky==1078&&kh==15)//
				//	kx=1182;//

				p=g_hist+(sym<<4|kh);
				if(p->srcx==-1)//first time
					break;

				ptrdiff_t srcidx=(ptrdiff_t)bw*p->srcy+p->srcx, dstidx=(ptrdiff_t)bw*ky+kx;
				if(mask[0]->data[srcidx]&0x7F||symbytes>1&&memcmp(buf+bytestride*srcidx+1, buf+bytestride*dstidx+1, (size_t)symbytes-1))//check if first src pixel is clean, first pixel must match
					continue;
				int blockw=1, blockh=1, updatew=1, updateh=1;
				do
				{
					if(updateh)
					{
						if(blockh<(1<<(sizeof(DeltaType)<<3))&&ky+blockh<bh&&p->srcy+blockh<bh)
						{
							ptrdiff_t offset=(ptrdiff_t)bw*blockh, srcidx2=srcidx+offset, dstidx2=dstidx+offset;
							if(strided_memcmp(buf+bytestride*dstidx2, buf+bytestride*srcidx2, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx2, mask[0]->data+srcidx2))
								++blockh;
							else
								updateh=0;
						}
						else
							updateh=0;
					}
					if(updatew)
					{
						if(blockw<(1<<(sizeof(DeltaType)<<3))&&kx+blockw<bw&&p->srcx+blockw<bw)
						{
							ptrdiff_t srcidx2=srcidx+(ptrdiff_t)blockw, dstidx2=dstidx+(ptrdiff_t)blockw;
							if(strided_memcmp(buf+bytestride*dstidx2, buf+bytestride*srcidx2, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx2, mask[0]->data+srcidx2))
								++blockw;
							else
								updatew=0;
						}
						else
							updatew=0;
					}
				}while(updateh||updatew);
				if(bestw*besth<blockw*blockh)
					bestw=blockw, besth=blockh, bestmatch=kh;
			}

			int redundantbytes=symbytes*bestw*besth,
			//	success=bestmatch!=0xFF&&redundantbytes>minlzbytes;
				success=bestmatch!=0xFF&&redundantbytes>sizeof(LZ2DInfo);
			if(success)//ignore if repeated bytes <= emission bytes
			{
				LZ2DInfo *emit=(LZ2DInfo*)ARRAY_APPEND(*lz, 0, 1, 1, 0);//emit repeated block
				p=g_hist+(sym<<4|bestmatch);
				emit->srcx=p->srcx;
				emit->srcy=p->srcy;
				emit->w=bestw-1;//1~256 - 1
				emit->h=besth-1;
				emit->dstx=kx;
				emit->dsty=ky;
				
#if 0
				{//
					for(int ky2=0;ky2<besth;++ky2)//DEBUG check if src clean
					{
						for(int kx2=0;kx2<bestw;++kx2)
						{
							unsigned char src0=mask[0]->data[bw*(p->srcy+ky2)+p->srcx+kx2];
							if(src0&0x7F)//LZ src'es can overlap
							{
								LOG_ERROR("Debug check: LZ src is not clean");
								return 0;
							}
						}
					}
				}//
#endif

				//false-mask src block, avoiding LZ chains (but with worse ratio)
				for(int ky2=0;ky2<besth;++ky2)
					memset(mask[0]->data+bw*(p->srcy+ky2)+p->srcx, 128, bestw);

				//mask dst block (can overlap with src)
				unsigned char val=32+rand()%31;//32~63, for debugging
				for(int ky2=0;ky2<besth;++ky2)
					memset(mask[0]->data+bw*(ky+ky2)+kx, val, bestw);

				savedbytes+=redundantbytes;
			}

			//update g_hist
			if(kh>15)
				kh=15;
			p=g_hist+((size_t)sym<<4|kh);
			p->srcx=kx;
			p->srcy=ky;
			if(kh)
			{
				LZ2DSrcInfo temp, *p2=g_hist+((size_t)sym<<4);
				SWAPVAR(*p, *p2, temp);
			}

			if(success)
				kx+=bestw;
			else
				++kx;
		}
	}
#endif
	return savedbytes;
}

//unsigned char *buf0=0;
static void image_diff_masked(unsigned char *buf, int iw, int ih, int nch, int bytestride, const unsigned char *mask)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx2=iw*ih-1,//differeniation is last -> first
			idx=idx2*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, --idx2, idx-=bytestride)
			{
				if(!(mask[idx2]&0x7F))//128 is a false mask for LZ src to avoid LZ chaining
				{
					unsigned char
						left=kx&&!(mask[idx2-1]&0x7F)?buf[idx-bytestride]:0,
						top=ky&&!(mask[idx2-iw]&0x7F)?buf[idx-rowlen]:0,
						topleft=kx&&ky&&!(mask[idx2-iw-1]&0x7F)?buf[idx-rowlen-bytestride]:0,
						sub=left+top-topleft;
					//if(kx||ky)
					if(sub)
						sub-=128;

					//if(kx==6&&ky==3)//DEBUG
					////if(kx==1919&&ky==1079)//
					//{
					//	printf("XYC %5d %5d %d:\n", kx, ky, kc);
					//	printf("  %3d %3d\n", topleft, top);
					//	printf("  %3d %3d -> %3d\n", left, buf[idx], (unsigned)(unsigned char)(buf[idx]-sub));
					//}

					buf[idx]-=sub;
				}
				//else//DEBUG
				//	buf[idx]=0;//
			}
		}
	}
}
static void image_int_masked(unsigned char *buf, int iw, int ih, int nch, int bytestride, const unsigned char *mask)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx2=0,//integration is first ->last
			idx=idx2*bytestride+kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx2, idx+=bytestride)
			{
				if(!mask[idx2])
				{
					unsigned char
						left=kx&&!mask[idx2-1]?buf[idx-bytestride]:0,
						top=ky&&!mask[idx2-iw]?buf[idx-rowlen]:0,
						topleft=kx&&ky&&!mask[idx2-iw-1]?buf[idx-rowlen-bytestride]:0,
						sub=left+top-topleft;
					//if(kx||ky)
					if(sub)
						sub-=128;

					//if(kx==6&&ky==3)//DEBUG
					////if(kx==1919&&ky==1079)//
					//{
					//	printf("XYC %5d,%5d,%2d:\n", kx, ky, kc);
					//	printf("  %3d %3d\n", topleft, top);
					//	printf("  %3d %3d -> %3d\n", left, buf[idx], (unsigned)(unsigned char)(buf[idx]+sub));
					//}

					buf[idx]+=sub;//addition instead of subtraction

					//if(buf0&&buf[idx]!=buf0[idx])//DEBUG
					//{
					//	printf("XYC %5d,%5d,%2d:  dec != buf0  0x%02X != 0x%02X  %d != %d\n", kx, ky, kc, buf[idx], buf0[idx], buf[idx], buf0[idx]);
					//	LOG_ERROR("Integration error");
					//}
				}
			}
		}
	}
}

static int rans_calc_CDF_masked(const unsigned char *buffer, int len, int stride, const unsigned char *mask, unsigned short *CDF)
{
	const int nlevels=256;
	if(!len||!stride)
	{
		LOG_ERROR("rans_calc_CDF_masked(): Invalid args");
		return 0;
	}
	HistInfo hist[256];
	for(int k=0;k<nlevels;++k)
	{
		hist[k].sym=k;
		hist[k].freq=0;
	}
	int count=0;
	for(int k=0, km=0;k<len;k+=stride, ++km)
	{
		int inc=!(mask[km]&0x7F);//128 is a false mask for LZ src to avoid LZ chaining
		hist[buffer[k]].freq+=inc;
		count+=inc;
	}
	if(!count)
		return 0;

#if 0
	static int call=0;
	if(call==3)
		print_hist(&hist->freq, sizeof(HistInfo)/sizeof(int), 256);//
	++call;
#endif

	//int count=len/stride;
	for(int k=0;k<nlevels;++k)
		hist[k].qfreq=((long long)hist[k].freq<<16)/count;
	
	if(count!=0x10000)
	{
		const int prob_max=0xFFFF;

		isort(hist, nlevels, sizeof(HistInfo), histinfo_byfreq);
		int idx=0;
		for(;idx<nlevels&&!hist[idx].freq;++idx);//skip absent symbols
		for(;idx<nlevels&&!hist[idx].qfreq;++idx)//restore rare symbols
			++hist[idx].qfreq;
		//for(idx=nlevels-1;idx>=0&&hist[idx].qfreq>=prob_max;--idx);
		//for(++idx;idx<nlevels;++idx)
		//	hist[idx].qfreq=prob_max;

		int error=-0x10000;//too much -> +ve error & vice versa
		for(int k=0;k<nlevels;++k)
			error+=hist[k].qfreq;
		if(error>0)
		{
			while(error)
			{
				for(int k=0;k<nlevels&&error;++k)
				{
					int dec=hist[k].qfreq>1;
					hist[k].qfreq-=dec, error-=dec;
				}
			}
		}
		else
		{
			while(error)
			{
				for(int k=nlevels-1;k>=0&&error;--k)
				{
					int inc=hist[k].qfreq<prob_max;
					hist[k].qfreq+=inc, error+=inc;
				}
			}
		}
		isort(hist, nlevels, sizeof(HistInfo), histinfo_bysym);
	}
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		//freq[k]=hist[k].qfreq;
		CDF[k]=sum;
		sum+=hist[k].qfreq;
		//if(sum>0xFFFF)
		//	sum=0xFFFF;
	}
	return count;
}


//	#define SAVE_DIFF2
	#define SAVE_MASK_ENC		//just use this for debugging
//	#define SAVE_DIFF1
//	#define SAVE_MASK_DEC
//	#define SAVE_OUT

//	#define PRINT_TEST3_rANS


//unsigned short g_freq[256*4];
static const int tag_lz04='L'|'Z'<<8|'0'<<16|'4'<<24;
size_t test3_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data, int diff)
{
	size_t start=0;
	if(*data)
	{
		if(data[0]->esize!=1)
			return 0;
		start=data[0]->count;
	}
	else
	{
		ARRAY_ALLOC(char, *data, 0, 0, 0, 0);
		start=0;
	}

	const unsigned char *buf=(const unsigned char*)src;
	ptrdiff_t srcbytes=(ptrdiff_t)bytestride*bw*bh;
	unsigned char *b2=(unsigned char*)malloc(srcbytes);
	if(!b2)
	{
		LOG_ERROR("test3_encode(): Allocation error");
		return 0;
	}
	memcpy(b2, buf, srcbytes);

	if(diff==2)
	{
		differentiate_image(b2, bw, bh, symbytes, bytestride);
#ifdef SAVE_DIFF2
		lodepng_encode_file("diff.PNG", b2, bw, bh, LCT_RGBA, 8);//DEBUG SAVE
#endif
	}

	ArrayHandle mask=0, rle=0, lz=0;
	size_t savedbytes=lz2d3_encode(b2, bw, bh, symbytes, bytestride, &mask, &rle, &lz);
	if(!mask||!rle||!lz)
		return 0;
	
#ifdef SAVE_MASK_ENC
	lodepng_encode_file("mask_enc.PNG", mask->data, bw, bh, LCT_GREY, 8);//DEBUG SAVE
#endif

	if(diff==1)
	{
		image_diff_masked(b2, bw, bh, symbytes, bytestride, mask->data);
#ifdef SAVE_DIFF1
		lodepng_encode_file("diff.PNG", b2, bw, bh, LCT_RGBA, 8);//DEBUG SAVE
#endif
	}

	int rANS_pixelcount=0;
#ifdef PRINT_TEST3_rANS
	{//
		for(int km=0;km<mask->count;++km)
			rANS_pixelcount+=!(mask->data[km]&0x7F);
		printf("Unmasked count %d\n", rANS_pixelcount);
		rANS_pixelcount=0;
	}//
#endif
	for(int kc=0;kc<symbytes;++kc)
	{
		int count=rans_calc_CDF_masked(b2+kc, (int)srcbytes, bytestride, mask->data, g_CDF+((size_t)kc<<8));
		//if(!rans_calc_CDF_masked(b2+kc, srcbytes, bytestride, mask->data, g_freq+((size_t)kc<<8), g_CDF+((size_t)kc<<8)))
		if(!count)
			return 0;
		if(rANS_pixelcount&&rANS_pixelcount!=count)
		{
			LOG_ERROR("rANS mask channel mismatch");
			return 0;
		}
		rANS_pixelcount=count;
	}

	//printf("rANS pixel count %d\n", rANS_pixelcount);//

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, (size_t)4*4);//tag, rle_count, lz_count, rANS_count
	dlist_push_back(&list, rle->data, rle->esize*rle->count);
	dlist_push_back(&list, lz->data, lz->esize*lz->count);
	dlist_push_back(&list, g_CDF, ((size_t)symbytes<<8)*sizeof(short));
	size_t rANS_bytes=list.nobj;

	//unsigned s0=0;//

	unsigned state=0x10000;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
		{
			int idx=bw*ky+kx;
			if(!(mask->data[idx]&0x7F))//128 is a false mask for LZ src to avoid LZ chaining
			{
				--rANS_pixelcount;
				//s0=state;//

				idx*=bytestride;
				for(int kc=0;kc<symbytes;++kc)
				{
					unsigned char sym=b2[idx+kc];
					//unsigned short cdf=g_CDF[kc<<8|sym], freq=g_freq[kc<<8|sym];
					unsigned short cdf=g_CDF[kc<<8|sym], freq;
					{
						int next;
						if(sym==255)
							next=0x10000;
						else
						{
							next=g_CDF[kc<<8|(sym+1)];
							if(next<cdf)//CDF will overflow with trailing zero frequency symbols  (uint16)0x10000==0
								next=0x10000;
						}
						freq=(unsigned short)(next-cdf);
#ifdef _DEBUG
						if(!freq)
							LOG_ERROR("rANS zero freq sym 0x%02X at XYC %5d,%5d,%2d", sym, kx, ky, kc);
#endif
					}

					if(state>=(unsigned)(freq<<16))//renorm
					{
						dlist_push_back(&list, &state, 2);
						state>>=16;
					}

					state=state/freq<<16|(cdf+state%freq);//update
					
#ifdef PRINT_TEST3_rANS
					if(rANS_pixelcount<2)//
						printf("enc [%d] XYC %5d %5d %d%c 0x%02X %3d  cdf 0x%04X f 0x%04X  x2 0x%08X\n", rANS_pixelcount, kx, ky, kc, "RGB"[kc], sym, sym, cdf, freq, state);//
#endif
				}
			}
		}
	}
	if(rANS_pixelcount)
	{
		LOG_ERROR("rANS pixel count mismatch");
		return 0;
	}
	//printf("Last states: end-3 0x%08X, end 0x%08X\n", s0, state);//
	dlist_push_back(&list, &state, 4);
	rANS_bytes=list.nobj-rANS_bytes;

	dlist_appendtoarray(&list, data);

	memcpy(data[0]->data+start, &tag_lz04, 4);
	start+=4;

	memcpy(data[0]->data+start, &rle->count, 4);
	start+=4;

	memcpy(data[0]->data+start, &lz->count, 4);
	start+=4;

	memcpy(data[0]->data+start, &rANS_bytes, 4);
	start+=4;

	dlist_clear(&list);
	array_free(&mask);
	array_free(&rle);
	array_free(&lz);
	free(b2);
	return savedbytes;
}
int    test3_decode(const void *src, ptrdiff_t srclen, int iw, int ih, int symbytes, int bytestride, void *dst, int diff)
{
	const unsigned char
		*data=(const unsigned char*)src,
		*srcptr=data,
		*srcend=data+srclen;
	if(srcptr+(ptrdiff_t)4*4>=srcend||memcmp(data, &tag_lz04, 4))
		return 0;
	srcptr+=4;

	ptrdiff_t rle_count=0, lz_count=0, rANS_bytes=0;
	memcpy(&rle_count, srcptr, 4);
	srcptr+=4;

	memcpy(&lz_count, srcptr, 4);
	srcptr+=4;

	memcpy(&rANS_bytes, srcptr, 4);
	srcptr+=4;

	ptrdiff_t res=(ptrdiff_t)iw*ih;
	unsigned char *mask=(unsigned char*)malloc(res);
	if(!mask)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(mask, 0, res);
	
	//1. reconstruct pixel mask
	if(srcptr+rle_count*sizeof(LZ2DRLEInfo)>srcend)
	{
		LOG_ERROR("Unexpected EOF at RLE");
		return 0;
	}
	for(int k=0;k<rle_count;++k)
	{
		LZ2DRLEInfo temp;
		memcpy(&temp, srcptr, sizeof(LZ2DRLEInfo));
		srcptr+=sizeof(LZ2DRLEInfo);
		
		unsigned char val=0xC0+rand()%63;//for debugging
		memset(mask+iw*temp.y+temp.x+1, val, temp.w-1);//leave top-left pixel unmasked in RLE blocks
		for(int ky=1;ky<temp.h;++ky)
			memset(mask+iw*(temp.y+ky)+temp.x, val, temp.w);
	}

	if(srcptr+lz_count*sizeof(LZ2DInfo)>srcend)
	{
		LOG_ERROR("Unexpected EOF at LZ");
		return 0;
	}
	for(int k=0;k<lz_count;++k)
	{
		LZ2DInfo t2;
		memcpy(&t2, srcptr, sizeof(LZ2DInfo));
		srcptr+=sizeof(LZ2DInfo);

		int blockw=t2.w+1, blockh=t2.h+1;//0~255 + 1
		
		unsigned char val=32+rand()%31;//for debugging
		for(int ky=0;ky<blockh;++ky)
			memset(mask+iw*(t2.dsty+ky)+t2.dstx, val, blockw);
	}

#ifdef SAVE_MASK_DEC
	lodepng_encode_file("mask_dec.PNG", mask, iw, ih, LCT_GREY, 8);//DEBUG SAVE
#endif

	ptrdiff_t CDF_bytes=((size_t)symbytes<<8)*sizeof(short);//read CDF
	if(srcptr+CDF_bytes>srcend)
	{
		LOG_ERROR("Unexpected EOF at CDF");
		return 0;
	}
	memcpy(g_CDF, srcptr, CDF_bytes);
	srcptr+=CDF_bytes;

	for(int kc=0;kc<symbytes;++kc)//fill g_CDF2sym[]
	{
		int k2=0, end=0;
		for(int sym=0;sym<256;++sym, k2=end)
		{
			if(sym==255)
				end=0x10000;
			else
			{
				end=g_CDF[kc<<8|(sym+1)];
				if(end<k2)
					end=0x10000;
			}

			for(;k2<end;++k2)
				g_CDF2sym[kc<<16|k2]=sym;

			if(end==0x10000)//now redundant because	'k2' is recycled 'end'		//to avoid overwriting whole CDF2sym with 0xFF with trailing zero-freq symbols
				break;
		}
	}

#ifdef PRINT_TEST3_rANS
	int rANS_pixelcount=0;//
	{//
		for(ptrdiff_t km=0;km<res;++km)
			rANS_pixelcount+=!mask[km];
		printf("Unmasked count %d\n", rANS_pixelcount);
		rANS_pixelcount=0;
	}//
#endif
	
	//2. decode unmasked pixels
	unsigned char *buf=(unsigned char*)dst;
	memset(buf, 0, (size_t)bytestride*iw*ih);
	const unsigned char *srcstart=srcptr;
	srcptr+=rANS_bytes;
	if(srcptr>srcend)
	{
		LOG_ERROR("Unexpected EOF at rANS");
		return 0;
	}
	if(srcptr!=srcend)
	{
		LOG_ERROR("Dangling data");
		return 0;
	}
	unsigned state;
	srcptr-=4;
	if(srcptr<srcstart)
	{
		LOG_ERROR("rANS OOB");
		return 0;
	}
	memcpy(&state, srcptr, 4);
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx)
		{
			int idx=iw*ky+kx;
			if(!mask[idx])
			{
#ifdef PRINT_TEST3_rANS
				++rANS_pixelcount;//
#endif

				int idx2=idx*bytestride;
				for(int kc=symbytes-1;kc>=0;--kc)
				{
					unsigned short lo=(unsigned short)state;
					unsigned char sym=g_CDF2sym[kc<<16|lo];
					unsigned short cdf=g_CDF[kc<<8|sym], freq;
					{
						int next;
						if(sym==255)
							next=0x10000;
						else
						{
							next=g_CDF[kc<<8|(sym+1)];
							if(next<cdf)
								next=0x10000;
						}
						freq=(unsigned short)(next-cdf);
#ifdef _DEBUG
						if(!freq)
							LOG_ERROR("rANS zero freq sym 0x%02X at XYC %5d,%5d,%2d", sym, kx, ky, kc);
#endif
					}
					buf[idx2+kc]=sym;
					
#ifdef PRINT_TEST3_rANS
					if(rANS_pixelcount<3)//
						printf("dec [%d] XYC %5d %5d %d%c 0x%02X %3d  cdf 0x%04X f 0x%04X  x1 0x%08X\n", rANS_pixelcount-1, kx, ky, kc, "RGB"[kc], sym, sym, cdf, freq, state);//
#endif

					state=freq*(state>>16)+lo-cdf;//update

					if(state<0x10000)//renorm
					{
						srcptr-=2;
						if(srcptr<srcstart)
						{
							LOG_ERROR("rANS OOB");
							return 0;
						}
						state<<=16;
						memcpy(&state, srcptr, 2);
					}
				}
			}
		}
	}

	//3. masked integration (optional)
	if(diff==1)
		image_int_masked(buf, iw, ih, symbytes, bytestride, mask);
	
	//4. decode RLE & LZ simultaneously			decoded slower than encoded
#if 1
	const unsigned char
		*rleptr=data+(ptrdiff_t)4*4,
		*rleend=rleptr+rle_count*sizeof(LZ2DRLEInfo),
		*lzptr=rleend,
		*lzend=lzptr+lz_count*sizeof(LZ2DInfo);
	if(rle_count||lz_count)
	{
		LZ2DRLEInfo rle;
		LZ2DInfo lz;

		if(rleptr+sizeof(LZ2DRLEInfo)<=rleend)
		{
			memcpy(&rle, rleptr, sizeof(LZ2DRLEInfo));
			rleptr+=sizeof(LZ2DRLEInfo);
		}
		else
			memset(&rle, -1, sizeof(LZ2DRLEInfo));

		if(lzptr+sizeof(LZ2DInfo)<=lzend)
		{
			memcpy(&lz, lzptr, sizeof(LZ2DInfo));
			lzptr+=sizeof(LZ2DInfo);
		}
		else
			memset(&lz, -1, sizeof(LZ2DInfo));

		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				if(kx==0&&ky==511)//
					kx=0;//

				if(lz.dstx!=0xFFFF&&ky==lz.dsty&&kx==lz.dstx)//4.1. first check LZ
				{
					//if(kx==1895&&ky==1072)//
					//	kx=1895;//

					int blockw=lz.w+1, blockh=lz.h+1;//0~255 + 1

					int delta=bytestride*(iw-blockw);
					unsigned char
						*pd=buf+bytestride*(iw*lz.dsty+lz.dstx),
						*ps=buf+bytestride*(iw*lz.srcy+lz.srcx);
					for(int ky=0;ky<blockh;++ky)
					{
						for(int kx=0;kx<blockw;++kx)//copy pixels one by one because buffers may overlap
						{
							memcpy(pd, ps, bytestride);
							pd+=bytestride;
							ps+=bytestride;
						}
						pd+=delta;
						ps+=delta;
					}
					
					if(lzptr+sizeof(LZ2DInfo)<=lzend)
					{
						memcpy(&lz, lzptr, sizeof(LZ2DInfo));
						lzptr+=sizeof(LZ2DInfo);
					}
					else//LZ is over
					{
						memset(&lz, -1, sizeof(LZ2DInfo));
						if(rle.x==0xFFFF)
							break;
					}
				}
				if(rle.x!=0xFFFF&&ky==rle.y&&kx==rle.x)//4.2. then check RLE
				{
					unsigned char *start=buf+bytestride*(iw*rle.y+rle.x);
					if(rle.w>1)
						memfill(start+bytestride, start, bytestride*((size_t)rle.w-1), bytestride);
					for(int ky=1;ky<rle.h;++ky)
						memcpy(buf+bytestride*(iw*(rle.y+ky)+rle.x), start, bytestride*rle.w);
					
					if(rleptr+sizeof(LZ2DRLEInfo)<=rleend)
					{
						memcpy(&rle, rleptr, sizeof(LZ2DRLEInfo));
						rleptr+=sizeof(LZ2DRLEInfo);
					}
					else//RLE is over
					{
						memset(&rle, -1, sizeof(LZ2DRLEInfo));
						if(lz.dstx==0xFFFF)
							break;
					}
				}
			}
			if(rle.x==0xFFFF&&lz.dstx==0xFFFF)//both RLE and LZ are over
				break;
		}
	}
#endif

	//decode RLE & LZ separately				decoded faster than encoded
#if 0
	//3. copy LZ blocks
#if 1
	srcptr=data+(ptrdiff_t)4*4+rle_count*sizeof(LZ2DRLEInfo);
	for(int k=0;k<lz_count;++k)
	{
		LZ2DInfo t2;
		memcpy(&t2, srcptr, sizeof(LZ2DInfo));
		srcptr+=sizeof(LZ2DInfo);

		int blockw=t2.w+1, blockh=t2.h+1;//0~255 + 1

		int delta=bytestride*(iw-blockw);
		unsigned char
			*pd=buf+bytestride*(iw*t2.dsty+t2.dstx),
			*ps=buf+bytestride*(iw*t2.srcy+t2.srcx);
		for(int ky=0;ky<blockh;++ky)
		{
			for(int kx=0;kx<blockw;++kx)//copy pixels one by one because buffers may overlap
			{
				memcpy(pd, ps, bytestride);
				pd+=bytestride;
				ps+=bytestride;
			}
			pd+=delta;
			ps+=delta;
		}
		//	memcpy(buf+bytestride*(iw*(t2.dsty+ky)+t2.dstx), buf+bytestride*(iw*(t2.srcy+ky)+t2.srcx), t2.w);//X  can overlap
	}
#endif

	//4. fill RLE blocks
#if 1
	srcptr=data+(ptrdiff_t)4*4;
	for(int k=0;k<rle_count;++k)
	{
		LZ2DRLEInfo temp;
		memcpy(&temp, srcptr, sizeof(LZ2DRLEInfo));
		srcptr+=sizeof(LZ2DRLEInfo);

		//if(temp.x==128&&temp.y==128)//
		//	temp.x=128;//

		unsigned char *start=buf+bytestride*(iw*temp.y+temp.x);
		if(temp.w>1)
			memfill(start+bytestride, start, bytestride*((size_t)temp.w-1), bytestride);
		for(int ky=1;ky<temp.h;++ky)
			memcpy(buf+bytestride*(iw*(temp.y+ky)+temp.x), start, bytestride*temp.w);
	}
#endif
#endif

	if(symbytes+1==bytestride)//set alpha
	{
		size_t dstsize=bytestride*res;
		for(ptrdiff_t k=symbytes;k<(ptrdiff_t)dstsize;k+=bytestride)
			buf[k]=0xFF;
	}

	if(diff==2)
	{
#ifdef SAVE_DIFF2
		lodepng_encode_file("diff2.PNG", buf, iw, ih, LCT_RGBA, 8);//DEBUG SAVE
#endif
		integrate_image(buf, iw, ih, symbytes, bytestride);
	}
#ifdef SAVE_OUT
	lodepng_encode_file("out.PNG", buf, iw, ih, LCT_RGBA, 8);//DEBUG SAVE
#endif

	return 1;
}
void   test3_printsummary(ArrayHandle cdata, size_t savedbytes, size_t usize, int symbytes)
{
	int *p=(int*)cdata->data;
	const char *names[6]={"header", "RLE", "LZ", "CDF", "rANS", "total"};
	size_t components[6]=
	{
		(size_t)4*4,//header {tag, count_RLE, count_LZ, ransbytes}
		(size_t)p[1]*8,//RLE
		(size_t)p[2]*10,//LZ
		(size_t)symbytes<<9,//CDF		256*sizeof(short) per channel			can be replaced with sizeof(int) per channel using differentiation & fitting
		p[3],//rANS
		0,//total
	};
	const int ncomp=COUNTOF(components);
	for(int k=0;k<ncomp-1;++k)
		components[ncomp-1]+=components[k];

	//printf("test3 summary:\n");
	printf("masked %lld / %lld bytes %lf%%\n", savedbytes, usize, 100.*savedbytes/usize);
	for(int k=0;k<COUNTOF(components);++k)		//Windows 11 / MSVC 2022 console text gets messed up
	{
		int printed=(int)strlen(names[k]);
		//printf("%s", names[k]);
		fputs(names[k], stdout);
		if(printed<8)
			printf("%*s", 8-printed, "");

		printed=snprintf(g_buf, G_BUF_SIZE, "%lld", components[k]);
		if(printed<8)
			printf("%*s", 8-printed, "");
		printf("%lld bytes %6.2lf%% ", components[k], 100.*components[k]/components[ncomp-1]);

		const int graphwidth=75-(8+8+8);
		for(int k2=0, end=(int)((graphwidth*components[k]+(components[ncomp-1]>>1))/components[ncomp-1]);k2<end;++k2)
			fputc('*', stdout);

		printf("\n");
	}

	//printf("\ntotal %lld bytes  ", total);
	//if(total==cdata->count)
	//	printf("CORRECT\n");
	//else
	//	printf("WRONG: cdata is %lld bytes (%d bytes %s)\n", cdata->count, abs((int)(total-cdata->count)), cdata->count>total?"more":"less");
}