#include"battle.h"
#include<stdio.h>//for debugging
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<math.h>
#include"lodepng.h"//for testing
static const char file[]=__FILE__;

static void print_hist(int *hist, int nlevels)
{
	int vmax=0;
	for(int k=0;k<nlevels;++k)
	{
		if(vmax<hist[k])
			vmax=hist[k];
	}
	for(int k=0;k<nlevels;++k)
	{
		printf("%3d %6d ", k, hist[k]);
		for(int k2=0, end=hist[k]*64/vmax;k2<end;++k2)
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

	print_hist(hist, vmax);//

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
	print_hist(hist_bypass, max_bypass);
	printf("repeat:\n");
	print_hist(hist_repeat, max_repeat);
	free(hist_bypass);
	free(hist_repeat);
	return cycles;
#endif
}

//LZ2
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
					ARRAY_APPEND(*bypass, &sym, 1, 1, 0);
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
				ARRAY_APPEND(*bypass, &sym, 1, 1, 0);
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
			c=sym<<16/nlevels, freq=0x10000/nlevels;

		if(state>=(unsigned)(freq<<16))//renorm
		{
			dlist_push_back(list, &state, 2);
			state>>=16;
		}

		state=state/freq<<16|(c+state%freq);//update
	}
	dlist_push_back(list, &state, 4);
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
long long test2_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data)
{
	if(symbytes>4)
		return 0;
	long long cycles=__rdtsc();
	
	unsigned short CDF[256*4];
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

		b_type[ky]=type;
		memcpy(b2+rowlen*ky, testrow+rowlen*type, rowlen);
	}
	ArrayHandle coeff=0, bypass=0;
	lz2_encode(b2, (int)res, &coeff, &bypass);

	//lodepng_encode_file("out.PNG", b2, bw, bh, bytestride==4?LCT_RGBA:LCT_RGB, 8);//

	for(int kc=0;kc<symbytes;++kc)
		rans_calc_CDF(coeff->data+kc, (int)coeff->count, symbytes, CDF+((size_t)kc<<8));
	dlist_push_back(&list, &tag_lz03, 4);
	dlist_push_back(&list, &coeff->count, 4);
	dlist_push_back(&list, CDF, ((size_t)symbytes<<8)*sizeof(short));
	rans_bytes_encode(coeff->data, (int)coeff->count, symbytes, CDF, 256, &list);
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
long long test3_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *data)
{
	if(symbytes>4)
		return 0;
	long long cycles=__rdtsc();
	
	unsigned short CDF[256*4];
	DList list;
	dlist_init(&list, 1, 1024, 0);
	//const unsigned char *buf=(const unsigned char*)src;
	unsigned short *buf;
	haar_2d_fwd((const unsigned char*)src, bw, bh, symbytes, bytestride, 0, &buf);
	int rowlen=symbytes*bw, srcrowlen=bytestride*bw;
	unsigned char *testrow=(unsigned char*)malloc((size_t)rowlen*5);
	unsigned char *b_type=(unsigned char*)malloc(bh);
	size_t res=(size_t)rowlen*bh;
	//unsigned char *b2=(unsigned char*)malloc(res);
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

		b_type[ky]=type;
		memcpy(b2+rowlen*ky, testrow+rowlen*type, rowlen);
	}
	ArrayHandle coeff=0, bypass=0;
	lz2_encode(b2, (int)res, &coeff, &bypass);

	//lodepng_encode_file("out.PNG", b2, bw, bh, symbytes==4?LCT_RGBA:LCT_RGB, 8);//

	for(int kc=0;kc<symbytes;++kc)
		rans_calc_CDF(coeff->data+kc, (int)coeff->count, symbytes, CDF+((size_t)kc<<8));
	dlist_push_back(&list, &tag_lz03, 4);
	dlist_push_back(&list, &coeff->count, 4);
	dlist_push_back(&list, CDF, ((size_t)symbytes<<8)*sizeof(short));
	rans_bytes_encode(coeff->data, (int)coeff->count, symbytes, CDF, 256, &list);
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