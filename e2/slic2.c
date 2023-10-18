//SLIC2: A simple lossless image codec v2
//By:  Ayman Wagih Mohsen
//To compile:
//  gcc -O3 -DMAIN slic2.c lodepng.c -o slic2


//START OF HEADER
#ifndef _INC_SLIC_H
#define _INC_SLIC_H
#ifdef __cplusplus
extern "C"
{
#endif
	

//nch:   Must be from 1 to 4. The channels must be interleaved and packed.
//depth: Must be from [1~16]. If depth<=8, data must be in bytes, otherwise data must be in little-endian uint16's (shorts).
//src:   Must be unsigned integers shifted leftmost. For example:
//	A 5-bit subpixel must be stored like this: 0bXXXX_X000
//	A 14-bit subpixel must be stored like this: 0bXXXX_XXXX_XXXX_XX00
//ret_dummy_alpha:  Tells if image has redundant alpha channel
//Don't forget to free returned buffers
unsigned char* slic2_encode(int iw, int ih, int nch, int depth, const void *src, long long *ret_size);
void*          slic2_decode(const unsigned char *data, long long len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha);

//I/O wrappers, return FALSE on error
int   slic2_save(const char *filename, int iw, int ih, int nch, int depth, const void *src);
void* slic2_load(const char *filename, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha);
	

#ifdef __cplusplus
}
#endif
#endif//_INC_SLIC_H
//END OF HEADER


//START OF IMPLEMENTATION
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdio.h>   //FILE/fopen/fread/fwrite/fclose
#include<stdlib.h>  //malloc/free		(optional) exit() in slic2_error()
#include<string.h>  //memset/memcpy
#include<stdarg.h>  //just for for slic2_error(): va_list
#include<sys/stat.h>//stat: (standard) to get file size fast


//	#define DISABLE_PALETTE
//	#define DISABLE_RCT
//	#define DISABLE_PREDICTOR


static void slic2_memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, copied);
	while(copied<<1<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
static void slic2_rct(short *buf, int iw, int ih, int depth, int fwd)//reversible color transform: YCmCb-R
{
	size_t res=(size_t)iw*ih;
	int mask=(short)(0xFFFF0000>>depth);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			short
				r=buf[idx],
				g=buf[idx+res],
				b=buf[idx+(res<<1)];

			if(fwd)
			{
				r-=g;
				g+=r>>1&mask;
				b-=g;
				g+=b>>1&mask;
			}
			else
			{
				g-=b>>1&mask;
				b+=g;
				g-=r>>1&mask;
				r+=g;
			}

			buf[idx]=r;
			buf[idx+res]=g;
			buf[idx+(res<<1)]=b;
		}
	}
}
static void slic_predict(const short *src, short *dst, int iw, int ih, int fwd)//one channel pixels are packed
{
	const short *pixels=fwd?src:dst;
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ++idx)
		{
			short
				N =ky?pixels[idx-iw]:0,			//SLIC_LOAD(pixels,  0, -1),
				W =kx?pixels[idx-1]:0,			//SLIC_LOAD(pixels, -1,  0),
				NW=kx&&ky?pixels[idx-iw-1]:0;	//SLIC_LOAD(pixels, -1, -1);
			short vmin, vmax, pred;
			if(N<W)
				vmin=N, vmax=W;
			else
				vmin=W, vmax=N;

			if(NW<vmin)
				pred=vmax;
			else if(NW>vmax)
				pred=vmin;
			else
				pred=(short)(N+W-NW);//shouldn't overflow, this is clamped gradient

			pred^=-fwd;//negate pred if fwd
			pred+=fwd;
			dst[idx]=src[idx]+pred;
		}
	}
}

static int slic2_error(int line, const char *msg, ...)
{
	int k;
	va_list args;

	printf("SLIC2(%d): ", line);
	if(msg)
	{
		va_start(args, msg);
		vprintf(msg, args);
		va_end(args);
	}
	else
		printf("Error");
	printf("\n");

	printf("Enter 0 to continue ...");
	while(!scanf(" %d", &k));
	exit(0);//for debug reasons
	return k;
}
#define ERROR(MSG, ...) slic2_error(__LINE__, MSG, ##__VA_ARGS__)


//simple double-linked list
#if 1
typedef struct C2DNodeStruct
{
	struct C2DNodeStruct *prev, *next;
	unsigned char data[];
} C2DNodeHeader, *C2DNodeHandle;
typedef struct C2DListStruct
{
	C2DNodeHandle i, f;
	size_t          //objsize == 1 byte
		nodebytes,  //objpernode
		nnodes,		//node count
		nbytes;		//nobj: total byte count in list
} C2DList, *C2DListHandle;
static void slic2_dlist_init(C2DListHandle list, size_t nodebytes)
{
	list->i=list->f=0;
	list->nodebytes=nodebytes;
	list->nnodes=list->nbytes=0;//empty
}
static void slic2_dlist_clear(C2DListHandle list)
{
	C2DNodeHandle it;

	it=list->i;
	if(it)
	{
		while(it->next)
		{
			it=it->next;
			free(it->prev);
		}
		free(it);
		list->i=list->f=0;
		list->nbytes=list->nnodes=0;
	}
}
static unsigned char* slic2_dlist_appendtoarray(C2DListHandle list, size_t *ret_size)
{
	C2DNodeHandle it;
	unsigned char *dst=(unsigned char*)malloc(list->nnodes*list->nodebytes);
	if(!dst)
	{
		ERROR("Allocation error");
		return 0;
	}
	*ret_size=list->nbytes;

	it=list->i;
	for(size_t offset=0;it;)
	{
		memcpy(dst+offset, it->data, list->nodebytes);
		offset+=list->nodebytes;
		it=it->next;
	}
	return dst;
}
static void slic2_dlist_append_node(C2DListHandle list)
{
	C2DNodeHandle temp;

	temp=(C2DNodeHandle)malloc(sizeof(C2DNodeHeader)+list->nodebytes);
	if(!temp)
	{
		ERROR("Allocation error");
		return;
	}
	temp->next=0;
	if(list->nnodes)
	{
		temp->prev=list->f;
		list->f->next=temp;
	}
	else
	{
		temp->prev=0;
		list->i=temp;
	}
	list->f=temp;
	++list->nnodes;
}
static void slic2_dlist_fill_node(C2DListHandle list, size_t copysize, char **src, void *dst)
{
	list->nbytes+=copysize;

	if(*src)
		memcpy(dst, *src, copysize), *src+=copysize;
	else
		memset(dst, 0, copysize);
}
void* slic2_dlist_push_back(C2DListHandle list, const void *data, size_t count)
{
	size_t obj_idx, copysize;
	char *buffer;
	void *ret;
	
	buffer=(char*)data;
	ret=0;
	obj_idx=list->nbytes%list->nodebytes;
	if(obj_idx)
	{
		copysize=list->nodebytes<obj_idx+count?list->nodebytes-obj_idx:count;
		count-=copysize;
		ret=list->f->data+obj_idx;

		slic2_dlist_fill_node(list, copysize, &buffer, ret);
	}
	while(count)
	{
		slic2_dlist_append_node(list);
		
		copysize=list->nodebytes<count?list->nodebytes:count;
		count-=copysize;

		if(!ret)
			ret=list->f->data;
		slic2_dlist_fill_node(list, copysize, &buffer, list->f->data);
	}
	return ret;
}
#endif

//ANS validation (should be removed)

//	#define ENABLE_EC_VALIDATION

#ifdef ENABLE_EC_VALIDATION
//single-linked list  for EC validation
#if 1
typedef struct C2SNodeStruct
{
	struct C2SNodeStruct *prev;
	unsigned char data[];//4-byte aligned on 32-bit, not suitable for double on 32-bit
} C2SNode, *C2SNodeHandle;
typedef struct C2SListStruct
{
	//[front] -> ... -> [back] -> nullptr
	size_t esize, count;
	void (*destructor)(void*);
	C2SNodeHandle
		front,	//can remove from or append to front
		back;	//prev always nullptr, can only append to back
} C2SList, *C2SListHandle;
static void c2_slist_init(C2SListHandle list, size_t esize, void (*destructor)(void*))
{
	list->esize=esize;
	list->count=0;
	list->destructor=destructor;
	list->front=list->back=0;
}
static void c2_slist_clear(C2SListHandle list)
{
	C2SNodeHandle temp;

	while(list->front)
	{
		temp=list->front;//copy pointer
		list->front=temp->prev;//advance

		if(list->destructor)//destroy & free
			list->destructor(temp->data);
		free(temp);
	}
	list->back=0;
	list->count=0;
}
static C2SNodeHandle c2_slist_alloc_node(C2SListHandle list, C2SNodeHandle prev, const void *data)
{
	C2SNodeHandle temp;

	//allocate new node
	temp=(C2SNodeHandle)malloc(sizeof(C2SNode)+list->esize);
	if(!temp)
	{
		ERROR("Allocation error");
		return 0;
	}
	temp->prev=prev;
	if(data)
		memcpy(temp->data, data, list->esize);
	else
		memset(temp->data, 0, list->esize);
	return temp;
}
static void* c2_slist_push_front(C2SListHandle list, const void *data)
{
	C2SNodeHandle temp;

	//allocate new node
	temp=c2_slist_alloc_node(list, list->front, data);

	if(list->count)//if front is not nullptr
		list->front=temp;//assign the new front
	else//list was empty
		list->front=list->back=temp;//initialize front and back with same node
	
	++list->count;
	return temp->data;
}
static void* c2_slist_push_back(C2SListHandle list, const void *data)
{
	C2SNodeHandle temp;

	//allocate new node
	temp=c2_slist_alloc_node(list, 0, data);//back->prev is always nullptr

	if(list->count)//if back is not nullptr
	{
		list->back->prev=temp;//make old back point to new note (new back)
		list->back=temp;//assign the new back
	}
	else//list was empty
		list->front=list->back=temp;//initialize front and back with same node
	
	++list->count;
	return temp->data;
}
static void* c2_slist_front(C2SListHandle list)
{
	if(!list->front)
		return 0;
	return list->front->data;
}
static void* c2_slist_back(C2SListHandle list)
{
	if(!list->back)
		return 0;
	return list->back->data;
}
static void c2_slist_pop_front(C2SListHandle list)
{
	C2SNodeHandle temp;

	if(!list->count)
		return;
	temp=list->front;//copy front pointer
	list->front=temp->prev;//advance

	if(list->destructor)//destroy & free
		list->destructor(temp->data);
	free(temp);
	
	--list->count;
}
static void c2_slist_print(C2SListHandle list, void (*printer)(const void*))
{
	C2SNodeHandle node=list->front;
	while(node)
	{
		printer(node->data);
		node=node->prev;
	}
}

//list-based stack
#define C2_STACK_PUSH(LIST, DATA) c2_slist_push_front(LIST, DATA)
#define C2_STACK_TOP(LIST)        c2_slist_front(LIST)
#define C2_STACK_POP(LIST)        c2_slist_pop_front(LIST)

//list-based queue
#define C2_QUEUE_ENQUEUE(LIST, DATA) c2_slist_push_back(LIST, DATA)
#define C2_QUEUE_FRONT(LIST)         c2_slist_front(LIST)
#define C2_QUEUE_DEQUEUE(LIST)       c2_slist_pop_front(LIST)
#endif

typedef struct DebugANSInfoStruct
{
	unsigned s0, state, cdf, freq, id, kx, ky;
	unsigned char kq, kc;
	unsigned short sym;
} DebugANSInfo;
//extern C2SList states;
static void debug_enc_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym);
static void debug_dec_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym);

static C2SList states={0};
static int debug_channel=0;
static void debug_ans_print(DebugANSInfo *info, int dec)
{
	printf("%6d state 0x%08X%s0x%08X cdf 0x%04X freq 0x%04X sym 0x%02X\n", info->id, info->s0, dec?"<-":"->", info->state, info->cdf, info->freq, info->sym);
}
static void debug_enc_dump(DebugANSInfo *i2)
{
	debug_ans_print(i2, 1);
	printf("\n");
	for(int k=0;k<20&&states.count;++k)
	{
		DebugANSInfo *i0=(DebugANSInfo*)C2_STACK_TOP(&states);
		debug_ans_print(i0, 0);
		C2_STACK_POP(&states);
	}
}
static void debug_enc_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		unsigned s0=state;

		if(!states.count)
			c2_slist_init(&states, sizeof(DebugANSInfo), 0);

		state=state/freq<<16|(cdf+state%freq);//enc update

		DebugANSInfo info={s0, state, cdf, freq, (unsigned)states.count, kx, ky, kq, kc, sym};
		C2_STACK_PUSH(&states, &info);
	}
}
static void debug_dec_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		if(!states.count)
			ERROR("Nothing to decode");
		DebugANSInfo *i0=(DebugANSInfo*)C2_STACK_TOP(&states), info;
		memcpy(&info, i0, sizeof(info));

		unsigned s0=freq*(state>>16)+(unsigned short)state-cdf;//dec update

		if(info.s0!=s0||info.state!=state||info.cdf!=cdf||info.freq!=freq||kx!=info.kx||ky!=info.ky||kq!=info.kq||kc!=info.kc||info.sym!=sym)
		{
			DebugANSInfo i2={s0, state, cdf, freq, (unsigned)states.count-1, kx, ky, kq, kc, sym};
			debug_enc_dump(&i2);
			ERROR("Decode error  (%d decodes remaining)", info.id);
		}

		C2_STACK_POP(&states);
	}
}
#else
#define debug_enc_update(...)
#define debug_dec_update(...)
#endif

//asymmetric numeral systems coder
typedef struct ANSEncoderStruct
{
	unsigned state;
	C2DList *list;
} ANSEncoder;
static void slic2_ans_enc_init(ANSEncoder *ec, C2DList *list)
{
	ec->state=0x10000;
	ec->list=list;
}
static void slic2_ans_enc(ANSEncoder *ec, int sym, const unsigned *CDF, int nlevels)//CDF length is nlevels+1
{
	int cdf, freq;
	if(CDF)
		cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	else//bypass
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	if(!freq)
		ERROR("ZPS");
	if(ec->state>=(unsigned)(freq<<16))//renorm
	{
		slic2_dlist_push_back(ec->list, &ec->state, 2);
		ec->state>>=16;
	}
	debug_enc_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=ec->state/freq<<16|(cdf+ec->state%freq);//update
}
static void slic2_ans_enc_flush(ANSEncoder *ec)
{
	slic2_dlist_push_back(ec->list, &ec->state, 4);
}
typedef struct ANSDecoderStruct
{
	unsigned state;
	const unsigned char *srcptr, *srcstart;
} ANSDecoder;
static void slic2_ans_dec_init(ANSDecoder *ec, const unsigned char *start, const unsigned char *end)
{
	ec->srcptr=end;
	ec->srcstart=start;
	
	ec->srcptr-=4;
	if(ec->srcptr<ec->srcstart)
		ERROR("ANS buffer overflow");
	memcpy(&ec->state, ec->srcptr, 4);
}
static int slic2_ans_dec(ANSDecoder *ec, const unsigned *CDF, int nlevels)//CDF length is nlevels+1
{
	unsigned c=(unsigned short)ec->state;
	int sym=0;
	unsigned cdf, freq;
	if(CDF)
	{
		int L=0, R=nlevels, found=0;
		while(L<=R)
		{
			sym=(L+R)>>1;
			if(CDF[sym]<c)
				L=sym+1;
			else if(CDF[sym]>c)
				R=sym-1;
			else
			{
				found=1;
				break;
			}
		}
		if(found)
			for(;sym<256-1&&CDF[sym+1]==c;++sym);
		else
			sym=L+(L<256&&CDF[L]<c)-(L!=0);

		cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	}
	else//bypass
	{
		sym=c*nlevels>>16;
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	}
	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	if(!freq)
		ERROR("ZPS");
	ec->state=freq*(ec->state>>16)+c-cdf;//update
	if(ec->state<0x10000)//renorm
	{
		ec->state<<=16;
		if(ec->srcptr-2>=ec->srcstart)
		{
			ec->srcptr-=2;
			memcpy(&ec->state, ec->srcptr, 2);
		}
	}
	return sym;
}
#define SLIC2_ENCODER   ANSEncoder
#define SLIC2_ENC_INIT  slic2_ans_enc_init
#define SLIC2_ENC       slic2_ans_enc
#define SLIC2_ENC_FLUSH slic2_ans_enc_flush
#define SLIC2_DECODER   ANSDecoder
#define SLIC2_DEC_INIT  slic2_ans_dec_init
#define SLIC2_DEC       slic2_ans_dec


#define SLIC2_HISTLEN 40//TODO this should depend on given bitdepth
typedef struct SLI2HeaderStruct
{
	char tag[4];//"SLI2"
	int iw, ih;

	unsigned char nch, depth;//nch: 1, 2, 3 or 4, depth: [1~16]
	short alpha;

	unsigned short palettesizes[4];//per channel, nonzero means enabled (depth must be > 8)
	unsigned short hist[4][SLIC2_HISTLEN];
} SLI2Header;
static SLI2Header slic2_header;
#define FLOOR_LOG2_16BIT(LOGN, N, SH, TEMP)\
	TEMP=N,\
	SH=(TEMP>=1<<8)<<3,	LOGN =SH, TEMP>>=SH,\
	SH=(TEMP>=1<<4)<<2,	LOGN+=SH, TEMP>>=SH,\
	SH=(TEMP>=1<<2)<<1,	LOGN+=SH, TEMP>>=SH,\
	SH= TEMP>=1<<1,		LOGN+=SH;
unsigned char* slic2_encode(int iw, int ih, int nch, int depth, const void *src, long long *ret_size)
{
	if(iw<1||ih<1 || nch<1||nch>4 || depth<1||depth>16 || !src)
		return 0;
	size_t res=(size_t)iw*ih;
	short *buf=(short*)malloc(res*(nch+2LL)*sizeof(short));//the image (separate channels), followed by current residues, then the residues are replaced with {bypass}, followed by {nbits, token} channel
	unsigned short *palettes[4]={0};
	if(!buf)
	{
		ERROR("Allocation error");
		return 0;
	}
	memcpy(slic2_header.tag, "SLI2", 4);
	slic2_header.iw=iw;
	slic2_header.ih=ih;
	slic2_header.nch=nch;
	slic2_header.depth=depth;
	unsigned char truedepth[]={depth, depth, depth, depth};
	if(depth<=8)
	{
		const unsigned char *src0=(const unsigned char*)src;
		memset(slic2_header.palettesizes, 0, sizeof(slic2_header.palettesizes));
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				buf[res*kc+k]=(src0[nch*k+kc]<<8)-0x8000;
		}
	}
	else
	{
#ifndef DISABLE_PALETTE
		int histlen=1<<depth;
		int *hist=(int*)malloc(histlen*sizeof(int));
		if(!hist)
		{
			ERROR("Allocation error");
			return 0;
		}
#endif
		const unsigned short *src0=(const unsigned short*)src;
		for(int kc=0;kc<nch;++kc)
		{
#ifndef DISABLE_PALETTE
			memset(hist, 0, histlen*sizeof(int));
			for(int k=0;k<res;++k)
			{
				unsigned short val=src0[nch*k+kc];
				++hist[val];
			}
			int palettesize=0;
			for(int k=0;k<histlen;++k)
				palettesize+=hist[k]!=0;
			if(palettesize<=(1<<(depth-8)))
			{
				{
					int sh, temp;
					FLOOR_LOG2_16BIT(truedepth[kc], palettesize-1, sh, temp);//truedepth = number of bits in maximum unsigned symbol
					++truedepth[kc];
				}
				slic2_header.palettesizes[kc]=palettesize;
				palettes[kc]=(unsigned short*)malloc(palettesize*sizeof(short));
				for(int k=0, idx=0;k<histlen;++k)
				{
					if(hist[k])
					{
						palettes[kc][idx]=k;
						++idx;
					}
				}
				for(int k=0;k<res;++k)
				{
					unsigned short val=src0[nch*k+kc];
					int idx=0;
					int L=0, R=palettesize;
					while(L<=R)
					{
						idx=(L+R)>>1;
						if(palettes[kc][idx]<val)
							L=idx+1;
						else if(palettes[kc][idx]>val)
							R=idx-1;
						else
							break;
					}
					buf[res*kc+k]=(short)((idx-((palettesize+1)>>1))<<(16-truedepth[kc]));
				}
			}
			else
#endif
			{
				for(int k=0;k<res;++k)
					buf[res*kc+k]=src0[nch*k+kc]-0x8000;
			}
		}
#ifndef DISABLE_PALETTE
		free(hist);
#endif
	}
	short *pixels, *tokens;
	int has_alpha=0;
	if(nch==2||nch==4)//if there is alpha: check if alpha is redundant to store it in the header
	{
		pixels=buf+res*(nch-1LL);//last channel
		short alpha=pixels[0];
		for(int k=1;k<res;++k)
		{
			if(pixels[k]!=alpha)
			{
				has_alpha=1;
				break;
			}
		}
		if(!has_alpha)
			slic2_header.alpha=alpha;
	}
#ifndef DISABLE_RCT
	if(nch==3||nch==4)
	{
		int mindepth=truedepth[0];//don't let RCT increase true depth in case of using palette
		if(mindepth>truedepth[1])mindepth=truedepth[1];
		if(mindepth>truedepth[2])mindepth=truedepth[2];
		slic2_rct(buf, iw, ih, mindepth, 1);
	}
#endif
	
	unsigned CDF2[SLIC2_HISTLEN+1];
	C2DList list;
	slic2_dlist_init(&list, 1024);
	slic2_dlist_push_back(&list, &slic2_header, sizeof(slic2_header));
	for(int kc=0;kc<4;++kc)//insert palettes, if any
	{
		if(palettes[kc])
			slic2_dlist_push_back(&list, palettes[kc], slic2_header.palettesizes[kc]*sizeof(short));
	}
	pixels=buf+res*nch;
	tokens=pixels+res;
	SLIC2_ENCODER ec;
	SLIC2_ENC_INIT(&ec, &list);
	for(int kc=nch-1;kc>=0;--kc)//for each channel
	{
		if((nch==2||nch==4)&&!has_alpha&&kc==nch-1)
			continue;
		const unsigned short half=1<<(truedepth[kc]-1);
#ifdef DISABLE_PREDICTOR
		pixels=buf+res*kc;
#else
		slic_predict(buf+res*kc, pixels, iw, ih, 1);
#endif
		memset(CDF2, 0, SLIC2_HISTLEN*sizeof(int));
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx)
			{
				//if(kc==1&&idx==iw*2+506)
				//if(kc==0&&idx==2908)//
				//	printf("");

				short pixel=pixels[idx]>>(16-truedepth[kc]);

				unsigned short sym=(pixel<<1)^-((unsigned short)pixel>=half);//L-permutation
				sym&=(1<<truedepth[kc])-1;
				
				int token, bypass;//the scheme for high bitdepth from JPEG XL
				int nbits;
				if(sym)
				{
					int sh, temp;
					FLOOR_LOG2_16BIT(nbits, sym, sh, temp);
					++nbits;
				}
				else
					nbits=0;
				//int nbits=sym?floor_log2(sym)+1:0;
				if(nbits<=4)
					token=sym, bypass=0, nbits=0;
				else
				{
					nbits-=2;
					token=16+((nbits-3)<<1)+(sym>>nbits)-2;
					bypass=sym&((1<<nbits)-1);
				}
				pixels[idx]=bypass;
				tokens[idx]=nbits<<8|token;
				//if(token>=SLIC2_HISTLEN)//
				//{
				//	ERROR("Token error");
				//	return 0;
				//}
				++CDF2[token];
			}//x-loop
		}//y-loop

		for(int kt=0, sum=0, overflow=0;kt<SLIC2_HISTLEN;++kt)
		{
			unsigned freq=CDF2[kt];
			CDF2[kt]=(unsigned)((unsigned long long)sum*(0x10000LL-SLIC2_HISTLEN)/res)+kt;
			sum+=freq;

			slic2_header.hist[kc][kt]=(unsigned short)CDF2[kt];
		}
		CDF2[SLIC2_HISTLEN]=0x10000;

		for(int k=(int)res-1;k>=0;--k)
		{
			int token=tokens[k]&0xFF, nbits=tokens[k]>>8, bypass=pixels[k];
			if(nbits)
				SLIC2_ENC(&ec, bypass, 0, 1<<nbits);
			SLIC2_ENC(&ec, token, CDF2, SLIC2_HISTLEN);
		}
	}//channel-loop
	SLIC2_ENC_FLUSH(&ec);
	unsigned char *ret=slic2_dlist_appendtoarray(&list, ret_size);
	slic2_dlist_clear(&list);
	free(buf);
	for(int kc=0;kc<4;++kc)
	{
		if(palettes[kc])
			free(palettes[kc]);
	}
	SLI2Header *header2=(SLI2Header*)ret;
	memcpy(header2->hist, slic2_header.hist, sizeof(header2->hist));
	return ret;
}
void* slic2_decode(const unsigned char *src, long long len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha)
{
	if(!src||len<sizeof(SLI2Header)||!ret_iw||!ret_ih||!ret_nch||!ret_depth)
	{
		ERROR("Invalid file/args");
		return 0;
	}
	memcpy(&slic2_header, src, sizeof(slic2_header));
	int iw=slic2_header.iw, ih=slic2_header.ih, nch=slic2_header.nch, depth=slic2_header.depth, res=iw*ih;
	if(memcmp(slic2_header.tag, "SLI2", 4)||(unsigned)nch-1>4-1||(unsigned)depth-1>16-1)
	{
		ERROR("Invalid file");
		return 0;
	}
	unsigned char truedepth[4];
	const unsigned char *srcstart=src+sizeof(slic2_header);
	unsigned short *palettes[4]={0};
	for(int kc=0;kc<4;++kc)
	{
		if(slic2_header.palettesizes[kc]>0)
		{
			int bytesize=slic2_header.palettesizes[kc]*sizeof(short);
			palettes[kc]=(unsigned short*)malloc(bytesize);
			if(!palettes[kc])
			{
				ERROR("Allocation error");
				return 0;
			}
			memcpy(palettes[kc], srcstart, bytesize);
			srcstart+=bytesize;

			{
				int sh, temp;
				FLOOR_LOG2_16BIT(truedepth[kc], slic2_header.palettesizes[kc]-1, sh, temp);
				++truedepth[kc];
			}
		}
		else
			truedepth[kc]=slic2_header.depth;
	}
	short *dst=(short*)malloc(res*(nch+1LL)*sizeof(short));
	if(!dst)
	{
		ERROR("Allocation error");
		return 0;
	}
	short *buf=dst+res*nch;
	unsigned CDF2[SLIC2_HISTLEN+1]={0};
	SLIC2_DECODER ec;
	SLIC2_DEC_INIT(&ec, srcstart, src+len);
	for(int kc=0;kc<nch;++kc)
	{
		if((nch==2||nch==4)&&kc==nch-1&&ec.srcptr==ec.srcstart)
		{
			slic2_memfill(dst+res*kc, &slic2_header.alpha, res*sizeof(short), sizeof(short));
			continue;
		}
#ifdef DISABLE_PREDICTOR
		buf=dst+res*kc;
#endif

		for(int kt=0, overflow=0;kt<SLIC2_HISTLEN;++kt)
		{
			if(overflow)
				CDF2[kt]=0x10000-SLIC2_HISTLEN+kt;
			else
			{
				unsigned cdf=slic2_header.hist[kc][kt];
				CDF2[kt]=cdf;
				if(kt<SLIC2_HISTLEN-1)
					overflow|=cdf>slic2_header.hist[kc][kt+1];
			}
		}
		CDF2[SLIC2_HISTLEN]=0x10000;

		for(int k=0;k<res;++k)
		{
			//if(kc==1&&k==iw*2+506)
			//if(kc==0&&k==2908)//
			//if(kc==0&&k==2633)//
			//	printf("");

			int token=SLIC2_DEC(&ec, CDF2, SLIC2_HISTLEN);
			unsigned short sym;
			if(token<16)
				sym=token;
			else
			{
				int nbits=((token-16)>>1)+3;
				int bypass=SLIC2_DEC(&ec, 0, 1<<nbits);
				sym=1<<(nbits+1)|(token&1)<<nbits|bypass;
			}

			short pixel=(sym>>1)^-(sym&1);//inverse L-permutation
			buf[k]=pixel<<(16-truedepth[kc]);
		}
#ifndef DISABLE_PREDICTOR
		slic_predict(buf, dst+res*kc, iw, ih, 0);
#endif
	}//channel-loop
	
#ifndef DISABLE_RCT
	if(nch==3||nch==4)
	{
		int mindepth=truedepth[0];//don't let RCT increase true depth in case of using palette
		if(mindepth>truedepth[1])mindepth=truedepth[1];
		if(mindepth>truedepth[2])mindepth=truedepth[2];
		slic2_rct(dst, iw, ih, mindepth, 0);
	}
#endif

	for(int kc=0;kc<4;++kc)
	{
		if(palettes[kc])
		{
			short *pixels=dst+kc*res;
			for(int k=0;k<res;++k)
			{
				short val=pixels[k];
				val>>=16-truedepth[kc];
				val+=(slic2_header.palettesizes[kc]+1)>>1;
				if((unsigned)val>=slic2_header.palettesizes[kc])
				{
					ERROR("Palette error");
					return 0;
				}
				pixels[k]=(short)(palettes[kc][val]-0x8000);
			}
			free(palettes[kc]);
		}
	}

	int pxsize=depth<=8?sizeof(char):sizeof(short);
	int add_alpha=nch==3&&force_alpha;
	void *ret=malloc((size_t)res*(nch+add_alpha)*pxsize);
	if(!ret)
	{
		ERROR("Allocation error");
		return 0;
	}

	//pack pixels, interleave channels, and add half
	if(depth<=8)
	{
		unsigned char *dst2=(unsigned char*)ret;
		if(add_alpha)
		{
			int black=0xFF000000;
			slic2_memfill(dst2, &black, (size_t)res*4*sizeof(char), sizeof(black));
		}
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=(dst[res*kc+k]>>8)+0x80;
		}
	}
	else
	{
		unsigned short *dst2=(unsigned short*)ret;
		if(add_alpha)
		{
			long long black=0xFFFF000000000000;
			slic2_memfill(dst2, &black, (size_t)res*4*sizeof(short), sizeof(black));
		}
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=dst[res*kc+k]+0x8000;
		}
	}
	free(dst);
	if(ret_iw           )*ret_iw=iw;
	if(ret_ih           )*ret_ih=ih;
	if(ret_nch          )*ret_nch=add_alpha?4:nch;
	if(ret_depth        )*ret_depth=depth;
	if(ret_dummy_alpha  )*ret_dummy_alpha=add_alpha;
	return ret;
}

static ptrdiff_t get_filesize(const char *filename)//-1 not found,  0: not a file,  ...: regular file size
{
	struct stat info={0};
	int error=stat(filename, &info);
	if(error)
		return -1;
	if((info.st_mode&S_IFMT)==S_IFREG)
		return info.st_size;
	return 0;
}
int slic2_save(const char *filename, int iw, int ih, int nch, int depth, const void *src)
{
	long long ret_size=0;
	unsigned char *ret=slic2_encode(iw, ih, 4, depth, src, &ret_size);
	if(!ret)
		return 0;
	FILE *f=fopen(filename, "wb");
	size_t byteswritten=fwrite(ret, 1, ret_size, f);

	if(byteswritten!=ret_size)//OPTIONAL CHECK
		printf("Error saving \'%s\'\n", filename);
	
	fclose(f);
	free(ret);
	return 1;
}
void* slic2_load(const char *filename, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha)
{
	ptrdiff_t size=get_filesize(filename);
	FILE *f=fopen(filename, "rb");
	if(!f)
		return 0;
	unsigned char *cdata=(unsigned char*)malloc(size);
	if(!cdata)
	{
		ERROR("Allocation error");
		return 0;
	}
	size_t bytesread=fread(cdata, 1, size, f);

	if(bytesread!=size)//OPTIONAL CHECK
		printf("Error reading \'%s\'\n", filename);

	void *udata=slic2_decode(cdata, bytesread, ret_iw, ret_ih, ret_nch, ret_depth, ret_dummy_alpha, force_alpha);
	free(cdata);
	return udata;
}
//END OF IMPLEMENTATION


//TEST CODE		gcc -O3 -DMAIN slic2.c lodepng.c -o slic2
#ifdef MAIN
#include<time.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
#include"lodepng.h"
void print_usage(const char *programname)
{
	printf(
		"Usage:  \"%s\" operation input [output]\n"
		"Operation:\n"
		"  e  Encode image to .SLI\n"
		"  d  Decode .SLI to .PNG\n"
		"  t  Test an image wihtout saving (no output)",
		programname
	);
}
void* load_image(const char *filename, int *iw, int *ih, int *nch0, int *depth, clock_t *elapsed)//returned number of channels is always 4 for simplicity, nch0 is just for CR calculation
{
	*elapsed=clock();
	unsigned short *udata=stbi_load_16(filename, iw, ih, nch0, 4);
	*elapsed=clock()-*elapsed;
	if(!udata)
		return 0;
	int actually8bit=1;
	ptrdiff_t npx=(ptrdiff_t)*iw**ih*4;
	for(ptrdiff_t k=0;k<npx;++k)//check for 8-bit image: stbi_load_16() duplicates bytes in this case
	{
		unsigned short val=udata[k];
		if((val>>8)!=(val&0xFF))
		{
			actually8bit=0;
			break;
		}
	}
	if(!actually8bit)
	{
		*depth=16;
		return udata;
	}
	*depth=8;
	unsigned char *udata8=(unsigned char*)malloc(npx);
	if(!udata8)
	{
		ERROR("Allocaiton error");
		return 0;
	}
	for(ptrdiff_t k=0;k<npx;++k)
		udata8[k]=(unsigned char)udata[k];
	free(udata);
	return udata8;
}
int main(int argc, char **argv)
{
#ifndef _DEBUG
	if(argc!=3&&argc!=4)
	{
		print_usage(argv[0]);
		return 1;
	}
	const char *op=argv[1], *fn1=argv[2], *fn2=argv[3];
#else
	const char *op="T", *fn1="C:/Projects/datasets/dataset-kodak/kodim01.png", *fn2=0;
#endif
	if(strlen(op)!=1)
	{
		printf("INVALID operation \'%s\'\n", op);
		print_usage(argv[0]);
		return 1;
	}
	int iw=0, ih=0, nch0=0, nch=4, depth=16;
	clock_t t;
	ptrdiff_t size1=get_filesize(fn1), size2;
	int dummy_alpha=0;
	switch(*op&0xDF)
	{
	case 'E'://encode SLI
		{
			printf("Encode SLI\n");
			printf("Opening \'%s\'...\n", fn1);
			void *udata;

			udata=load_image(fn1, &iw, &ih, &nch0, &depth, &t);
			if(!udata)
			{
				printf("CANNOT open \'%s\'\n", fn1);
				return 1;
			}
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			if(nch0<4)
				dummy_alpha=1;
			printf("C*W*H*D  %d*%d*%d*%d\n", nch0, iw, ih, depth);
		
			printf("Saving \'%s\'...\n", fn2);
			t=clock();
			long long ret_size=0;
			if(!slic2_save(fn2, iw, ih, 4, depth, udata))
				printf("ERROR saving \'%s\'\n", fn2);
			t=clock()-t;
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			free(udata);
			size2=get_filesize(fn2);
		}
		break;
	case 'D'://decode SLI
		{
			printf("Decode SLI\n");
			printf("Opening \'%s\'...\n", fn1);
			t=clock();
			void *udata=slic2_load(fn1, &iw, &ih, &nch, &depth, &dummy_alpha, 1);
			t=clock()-t;
			if(!udata)
			{
				printf("CANNOT open \'%s\'\n", fn1);
				return 1;
			}
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			printf("C*W*H*D  %d*%d*%d*%d\n", nch-dummy_alpha, iw, ih, depth);

			if(depth>8)//convert LE -> BE
			{
				printf("Converting LE to BE...\n");
				t=clock();
				ptrdiff_t npx=(ptrdiff_t)iw*ih*nch;
				unsigned short *ptr=(unsigned short*)udata;
				for(int k=0;k<npx;++k)
					ptr[k]=ptr[k]<<8|ptr[k]>>8;
				t=clock()-t;
				printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			}

			printf("Saving \'%s\'...\n", fn2);
			t=clock();
			LodePNGColorType type=LCT_RGBA;
			switch(nch)
			{
			case 1:type=LCT_GREY;break;
			case 2:type=LCT_GREY_ALPHA;break;
			case 3:type=LCT_RGB;break;
			case 4:type=LCT_RGBA;break;
			default:
				printf("Unsupported number of channels\n");
				break;
			}
			if(depth>8)
				depth=16;
			else
				depth=8;
			int error=lodepng_encode_file(fn2, (unsigned char*)udata, iw, ih, type, depth);
			t=clock()-t;
			if(error)
				printf("ERROR saving \'%s\'\n", fn2);
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			size2=get_filesize(fn2);
		}
		break;
	case 'T'://test SLI
		{
			printf("Test SLI\n");
			printf("Opening \'%s\'...\n", fn1);
			void *udata;
			
			udata=load_image(fn1, &iw, &ih, &nch0, &depth, &t);
			if(!udata)
			{
				printf("CANNOT open \'%s\'\n", fn1);
				return 1;
			}
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			printf("C*W*H*D  %d*%d*%d*%d\n", nch0, iw, ih, depth);
			
			printf("Test encode...\n");
			t=clock();
			long long ret_size=0;
			unsigned char *cdata=slic2_encode(iw, ih, 4, depth, udata, &ret_size);
			t=clock()-t;
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);
			
			printf("Test decode...\n");
			t=clock();
			int iw2=0, ih2=0, nch2=0, depth2=0, dummy_alpha2=0;
			unsigned char *udata2=slic2_decode(cdata, ret_size, &iw2, &ih2, &nch2, &depth2, &dummy_alpha2, 1);
			t=clock()-t;
			printf("Took %lf sec\n", (double)t/CLOCKS_PER_SEC);

			printf("Checking...\n");
			if(iw2!=iw||ih2!=ih||nch2!=4||depth2!=depth)
			{
				printf(
					"ERROR wrong image parameters:\n"
					"C*W*H*D:  %d*%d*%d*%d  !=  %d*%d*%d*%d\n",
					4, iw, ih, depth,
					nch2, iw2, ih2, depth2
				);
				return 1;
			}
			int success=1;
			for(ptrdiff_t k=0, npx=(ptrdiff_t)nch*iw*ih*depth>>3;k<npx;++k)
			{
				unsigned char
					v1=((unsigned char*)udata)[k],
					v2=udata2[k];
				if(v2!=v1)
				{
					int kc, kx, ky;
					if(depth>8)
						k>>=1;
					kc=k%nch;
					k/=nch;
					kx=k%iw;
					ky=(int)(k/iw);
					printf("ERROR  idx %lld  CXY %d %d %d\n", k, kc, kx, ky);
					printf("  original 0x%04X = %d\n", v1, v1);
					printf("  !=       0x%04X = %d\n", v2, v2);
					success=0;
					break;
				}
			}
			if(success)
				printf("SUCCESS\n");
			free(cdata);
			free(udata2);
			free(udata);
			size2=ret_size;
		}
		break;
	default:
		printf("Invalid operation \'%s\'\n", op);
		print_usage(argv[0]);
		return 1;
	}
	if(size1>0&&size2>0)
	{
		ptrdiff_t usize=(ptrdiff_t)iw*ih*nch0*depth>>3;
		printf(
			"Uncompressed/Compressed ratio (16-bit):\n"
			"  Uncompressed size C*W*H*D %d*%d*%d*%d = %lld bytes\n"
			"  Compressed size   %lld -> %lld bytes\n"
			"  Compression ratio %lf -> %lf\n",
			nch0, iw, ih, depth, usize, size1, size2, (double)usize/size1, (double)usize/size2
		);
		if(dummy_alpha)
			printf("  Dummy alpha\n");
	}
	else
		printf("Operation failed\n");
	return 0;
}
#endif
#undef  ERROR
//END OF TEST CODE
