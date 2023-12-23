#pragma once
#ifdef DEBUG_MEMORY
#include<stdlib.h>

void* my_malloc(size_t size);
void* my_realloc(void *p, size_t size);
void my_free(void *p);

#ifdef DEBUG_MEMORY_IMPLEMENTATION
#define NPOINTERS 1024
typedef struct PointerInfoStruct
{
	void *p0, *p;
	size_t size;
	char msg[8];
} PointerInfo;
PointerInfo pointers[NPOINTERS]={0};
int pointer_idx=0;
void* my_malloc(size_t size)
{
	void *p=malloc(size);
	PointerInfo *ptr=pointers+pointer_idx;
	ptr->p0=0;
	ptr->p=p;
	ptr->size=size;
	memset(ptr->msg, 0, sizeof(ptr->msg));
	memcpy(ptr->msg, "malloc", 6);
	pointer_idx=(pointer_idx+1)%NPOINTERS;
	return p;
}
void* my_realloc(void *p, size_t size)
{
	void *p2=realloc(p, size);
	PointerInfo *ptr=pointers+pointer_idx;
	ptr->p0=p;
	ptr->p=p2;
	ptr->size=size;
	memset(ptr->msg, 0, sizeof(ptr->msg));
	memcpy(ptr->msg, "realloc", 7);
	pointer_idx=(pointer_idx+1)%NPOINTERS;
	return p2;
}
void my_free(void *p)
{
	PointerInfo *ptr=pointers+pointer_idx;
	ptr->p0=p;
	ptr->p=0;
	ptr->size=0;
	memset(ptr->msg, 0, sizeof(ptr->msg));
	memcpy(ptr->msg, "free", 4);
	pointer_idx=(pointer_idx+1)%NPOINTERS;
	if(p==im0)
		LOG_ERROR2("");
	free(p);
}
#endif

#define free my_free
#define malloc my_malloc
#define realloc my_realloc
#endif