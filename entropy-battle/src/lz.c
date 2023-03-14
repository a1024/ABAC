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

//LZ2
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

static const int tag_lz03='L'|'Z'<<8|'0'<<16|'3'<<24;
static unsigned short g_CDF[256*4];
static unsigned char g_CDF2sym[65536*4];
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

	//lodepng_encode_file("out.PNG", b2, bw, bh, bytestride==4?LCT_RGBA:LCT_RGB, 8);//

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
		unsigned sum=0;
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
typedef struct LZ2DInfoStruct
{
	short srcx, srcy;
	DeltaType w, h;
	short dstx, dsty;
} LZ2DInfo;
typedef struct LZ2DRLEInfoStruct
{
	short x, y, w, h;
} LZ2DRLEInfo;
LZ2DSrcInfo g_hist[256*16];
int strided_valcmp(const unsigned char *p1, const unsigned char *p2, int symbytes, int bytestride, int pxstride, int count, int pad1, int pad2, const unsigned char *dstmask)
{
	int fullstride=bytestride*pxstride;
	const unsigned char *end=p1+(size_t)fullstride*count;
	if(pad1&&!dstmask[-pxstride]&&!memcmp(p1-fullstride, p2, symbytes))//check if entering a sea
		return 0;
	for(;p1<end;p1+=fullstride, dstmask+=pxstride)
	{
		if(*dstmask||memcmp(p1, p2, symbytes))
			return 0;
	}
	if(pad2&&!*dstmask&&!memcmp(p1, p2, symbytes))//check if entering a sea
		return 0;
	return 1;
}
int strided_memcmp(const unsigned char *p1, const unsigned char *p2, int symbytes, int bytestride, int pxstride, int count, const unsigned char *dstmask)
{
	int fullstride=bytestride*pxstride;
	const unsigned char *end=p1+(size_t)fullstride*count;
	for(;p1<end;p1+=fullstride, p2+=fullstride, dstmask+=pxstride)
	{
		if(*dstmask)
			return 0;
		for(int k=0;k<symbytes;++k)
		{
			if(p1[k]!=p2[k])
				return 0;
		}
	}
	return 1;
}
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
						if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx))
							++blockh;
						else
							updatedh=0;
					}
					else
						updatedh=0;
					if(updatedw&&blockw<(1<<(sizeof(DeltaType)<<3))&&blockw<blockh*4&&kx+blockw<bw&&p->srcx+blockw<bw)
					{
						srcidx=bw*p->srcy+p->srcx+blockw, dstidx=bw*ky+kx+blockw;
						if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx))
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
						if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx))
							++blockh;
						else
							updatedh=0;
					}
					else
						updatedh=0;
					if(updatedw&&blockw<(1<<(sizeof(DeltaType)<<3))&&blockw<blockh*4&&kx+blockw<bw&&p->srcx+blockw<bw)
					{
						srcidx=bw*p->srcy+p->srcx+blockw, dstidx=bw*ky+kx+blockw;
						if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx))
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

typedef struct CNCPointStruct
{
	short x, y;
} CNCPoint;
size_t lz2d3_encode(const unsigned char *buf, int bw, int bh, int symbytes, int bytestride, ArrayHandle *mask, ArrayHandle *rle, ArrayHandle *lz, int minlzbytes)
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
	for(int lgstep=1;lgstep<=8;++lgstep)
	{
		int step=1<<lgstep;
		int halfx=step>>1, halfy=bw*halfx;
		for(int ky=0;ky<bh;ky+=step)
		{
			for(int kx=0;kx<bw;kx+=step)
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
				unsigned char val=0xC0+rand()%63;
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

	//LZ2D
#if 1
	memset(g_hist, -1, sizeof(g_hist));
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
				int blockw=1, blockh=1, updatew=1, updateh=1;
				do
				{
					int srcidx, dstidx;
					if(updateh)
					{
						if(blockh<(1<<(sizeof(DeltaType)<<3))&&ky+blockh<bh&&p->srcy+blockh<bh)
						{
							srcidx=bw*(p->srcy+blockh)+p->srcx, dstidx=bw*(ky+blockh)+kx;
							if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, 1, blockw, mask[0]->data+dstidx))
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
							srcidx=bw*p->srcy+p->srcx+blockw, dstidx=bw*ky+kx+blockw;
							if(strided_memcmp(buf+bytestride*srcidx, buf+bytestride*dstidx, symbytes, bytestride, bw, blockh, mask[0]->data+dstidx))
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
				emit->w=bestw-1;
				emit->h=besth-1;
				emit->dstx=kx;
				emit->dsty=ky;

				unsigned char val=32+rand()%31;
				for(int ky2=0;ky2<besth;++ky2)//mask block
					memset(mask[0]->data+bw*(ky+ky2)+kx, val, bestw);
				//{
				//	for(int kx2=0;kx2<bestw;++kx2)
				//	{
				//		unsigned char *p3=mask[0]->data+bw*(ky+ky2)+kx+kx2;
				//		*p3=val;
				//		//*p3+=*p3<255;
				//	}
				//}

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

static void image_diff_masked(unsigned char *buf, int iw, int ih, int nch, int bytestride, const unsigned char *mask)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx2=iw*ih-1,
			idx=idx2*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, --idx2, idx-=bytestride)
			{
				if(!mask[idx2])
				{
					unsigned char
						left=kx&&!mask[idx2-1]?buf[idx-bytestride]:0,
						top=ky&&!mask[idx2-iw]?buf[idx-rowlen]:0,
						topleft=kx&&ky&&!mask[idx2-iw-1]?buf[idx-rowlen-bytestride]:0,
						sub=left+top-topleft;
					if(kx||ky)
						sub-=128;
					buf[idx]-=sub;
				}
			}
		}
	}
}
static void image_int_masked(unsigned char *buf, int iw, int ih, int nch, int bytestride, const unsigned char *mask)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx2=iw*ih-1,
			idx=idx2*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, --idx2, idx-=bytestride)
			{
				if(!mask[idx2])
				{
					unsigned char
						left=kx&&!mask[idx2-1]?buf[idx-bytestride]:0,
						top=ky&&!mask[idx2-iw]?buf[idx-rowlen]:0,
						topleft=kx&&ky&&!mask[idx2-iw-1]?buf[idx-rowlen-bytestride]:0,
						sub=left+top-topleft;
					if(kx||ky)
						sub-=128;
					buf[idx]+=sub;
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
	for(int k=0, km=0;k<len;k+=stride, ++km)//this loop takes 73% of encode time
	{
		int inc=!mask[km];
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
	return 1;
}

static const int tag_lz04='L'|'Z'<<8|'0'<<16|'4'<<24;
unsigned short g_freq[256*4];
int test3_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data, int minlzbytes)
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
	ArrayHandle mask=0, rle=0, lz=0;
	lz2d3_encode(buf, bw, bh, symbytes, bytestride, &mask, &rle, &lz, minlzbytes);
	if(!mask||!rle||!lz)
		return 0;

	unsigned char *b2=(unsigned char*)malloc((size_t)bw*bh*bytestride);
	if(!b2)
	{
		LOG_ERROR("test3_encode(): Allocation error");
		return 0;
	}
	memcpy(b2, buf, (size_t)bw*bh*bytestride);
	image_diff_masked(b2, bw, bh, symbytes, bytestride, mask->data);

	for(int kc=0;kc<symbytes;++kc)
	{
		if(!rans_calc_CDF_masked(b2+kc, bytestride*bw*bh, bytestride, mask->data, g_CDF+((size_t)kc<<8)))
		//if(!rans_calc_CDF_masked(b2+kc, bytestride*bw*bh, bytestride, mask->data, g_freq+((size_t)kc<<8), g_CDF+((size_t)kc<<8)))
			return 0;
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, (size_t)4*4);//tag, rle_count, lz_count, rANS_count
	dlist_push_back(&list, rle->data, rle->esize*rle->count);
	dlist_push_back(&list, lz->data, lz->esize*lz->count);
	dlist_push_back(&list, g_CDF, ((size_t)symbytes<<8)*sizeof(short));
	size_t rANS_bytes=list.nobj;

	unsigned state=0x10000;
	int idx=0;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bh;++kx, idx+=bytestride)
		{
			if(!mask->data[bw*ky+kx])
			{
				for(int kc=0;kc<symbytes;++kc)
				{
					unsigned char sym=b2[idx+kc];
					//unsigned short c=g_CDF[kc<<8|sym], freq=g_freq[kc<<8|sym];
					unsigned short c=g_CDF[kc<<8|sym], freq=(sym==255||c>g_CDF[kc<<8|(sym+1)]?0x10000:g_CDF[kc<<8|(sym+1)])-c;//CDF can overflow with zero frequency symbols at the end

					if(state>=(unsigned)(freq<<16))//renorm
					{
						dlist_push_back(&list, &state, 2);
						state>>=16;
					}

					state=state/freq<<16|(c+state%freq);//update
				}
			}
		}
	}
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
	return 1;
}