//double-linked LIST of identical size buffers
#ifndef INC_BLIST_H
#define INC_BLIST_H
#ifdef __GNUC__
#define BLIST_INLINE __attribute__((always_inline)) inline
#else
#define BLIST_INLINE __forceinline
#endif

//power of two
//	#define BLIST_PAYLOADSIZE 65536
	#define BLIST_PAYLOADSIZE 16384
//	#define BLIST_PAYLOADSIZE 4096
//	#define BLIST_PAYLOADSIZE 1024
//	#define BLIST_PAYLOADSIZE 256

typedef struct BNodeStruct
{
	struct BNodeStruct *next, *prev;
	unsigned char data[BLIST_PAYLOADSIZE];
} BNode;
typedef struct _BList
{
	BNode *i, *f;
	size_t
		nnodes,
		nbytes;
} BList;
//void blist_init(BList *list);
//size_t blist_appendtoarray(BList *list, ArrayHandle *dst);
//void blist_clear(BList *list);
//void* blist_push_back1(BList *list, const void *src);
//void* blist_push_back(BList *list, const void *data, size_t size);

BLIST_INLINE static void blist_init(BList *list)
{
	list->i=list->f=0;
	list->nnodes=list->nbytes=0;//empty
}
BLIST_INLINE static size_t blist_appendtoarray(BList *list, ArrayHandle *dst)
{
	BNode *it;
	size_t start;

	if(!*dst)
	{
		start=0;
		*dst=array_construct(0, 1, 0, 0, list->nnodes*BLIST_PAYLOADSIZE, 0);
	}
	else
	{
		if(dst[0]->esize!=1)
		{
			LOG_ERROR("blist_appendtoarray(): dst->esize=%d != 1", dst[0]->esize);
			return 0;
		}
		start=dst[0]->count;
		ARRAY_APPEND(*dst, 0, 0, 0, list->nnodes*BLIST_PAYLOADSIZE);
	}
	it=list->i;
	for(size_t offset=dst[0]->count;it;)
	{
		memcpy(dst[0]->data+offset, it->data, BLIST_PAYLOADSIZE);
		offset+=BLIST_PAYLOADSIZE;
		it=it->next;
	}
	dst[0]->count+=list->nbytes;
	return start;
}
BLIST_INLINE static void blist_appendtofile(BList *list, FILE *f)
{
	ptrdiff_t remaining=list->nbytes;
	BNode *it;

	it=list->i;
	for(;it->next;)
	{
		fwrite(it->data, 1, BLIST_PAYLOADSIZE, f);
		remaining-=BLIST_PAYLOADSIZE;
		it=it->next;
	}
	fwrite(it->data, 1, remaining, f);
}
BLIST_INLINE static void blist_clear(BList *list)
{
	BNode *it;

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
BLIST_INLINE static void blist_append_node(BList *list)
{
	BNode *temp=(BNode*)malloc(sizeof(BNode));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
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
BLIST_INLINE static void* blist_push_back1(BList *list, const void *src)
{
	unsigned char *p;
	size_t obj_idx=list->nbytes%BLIST_PAYLOADSIZE;//index of next object
	if(!obj_idx)//need a new node
		blist_append_node(list);
	p=list->f->data+obj_idx;
	*p=src?*(unsigned char*)src:0;
	++list->nbytes;
	return p;
}
BLIST_INLINE static void blist_fill_node(BList *list, size_t copysize, const unsigned char **src, void *dst)
{
	if(*src)
	{
		memcpy(dst, *src, copysize);
		*src+=copysize;
	}
	else
		memset(dst, 0, copysize);
	list->nbytes+=copysize;
}
BLIST_INLINE static void* blist_push_back(BList *list, const void *data, size_t size)
{
	size_t obj_idx, copysize;
	const unsigned char *buffer;
	void *ret;
	
	buffer=(const unsigned char*)data;
	ret=0;
	obj_idx=list->nbytes%BLIST_PAYLOADSIZE;
	if(obj_idx)
	{
		copysize=BLIST_PAYLOADSIZE<obj_idx+size?BLIST_PAYLOADSIZE-obj_idx:size;
		size-=copysize;
		ret=list->f->data+obj_idx;

		blist_fill_node(list, copysize, &buffer, ret);
	}
	while(size)
	{
		blist_append_node(list);
		
		copysize=BLIST_PAYLOADSIZE<size?BLIST_PAYLOADSIZE:size;
		size-=copysize;

		if(!ret)
			ret=list->f->data;
		blist_fill_node(list, copysize, &buffer, list->f->data);
	}
	return ret;
}
#endif
