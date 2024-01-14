//util.h - Utilities declarations
//Copyright (C) 2023  Ayman Wagih Mohsen
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once
#ifndef INC_UTIL_H
#define INC_UTIL_H
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stddef.h>//for size_t
#ifdef __cplusplus
extern "C"
{
#endif

//utility
#define COUNTOF(ARR)        (sizeof(ARR)/sizeof(*(ARR)))		//stdlib defines _countof
#define BETWEEN_INC(LO, X, HI) ((unsigned)((X)-LO)<(unsigned)(HI+1-LO))
#define BETWEEN_EXC(LO, X, HI) ((unsigned)((X)-LO)<(unsigned)(HI-LO))
#define SWAPVAR(A, B, TEMP)	TEMP=A, A=B, B=TEMP
#define SWAPMEM(A, B, TEMP, SIZE) memcpy(TEMP, A, SIZE), memcpy(A, B, SIZE), memcpy(B, TEMP, SIZE)
#define MINVAR(A, B)        ((A)<(B)?(A):(B))
#define MAXVAR(A, B)        ((A)>(B)?(A):(B))
#define CLAMP(LO, X, HI)    ((X)>(LO)?(X)<(HI)?(X):(HI):(LO))
//#define CLAMP(LO, X, HI)    ((X)<(LO)?(LO):((X)>(HI)?(HI):(X)))
#define MOVEOBJ(SRC, DST, SIZE) memcpy(DST, SRC, SIZE), memset(SRC, 0, SIZE)
#define MODVAR(DST, SRC, N) DST=(SRC)%(N), DST+=(N)&-(DST<0)

#ifdef _MSC_VER
#define	ALIGN(N)	__declspec(align(N))
#else
#define	ALIGN(N)	__attribute__((aligned(N)))
#endif

#ifndef _MSC_VER
#define sprintf_s	snprintf
#endif
#define G_BUF_SIZE	4096
extern char g_buf[G_BUF_SIZE];

void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes);
void memswap_slow(void *p1, void *p2, size_t size);
void memswap(void *p1, void *p2, size_t size, void *temp);
void memreverse(void *p, size_t count, size_t esize);//calls memswap
void memrotate(void *p, size_t byteoffset, size_t bytesize, void *temp);//temp buffer is min(byteoffset, bytesize-byteoffset)
int binary_search(const void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*), const void *val, size_t *idx);//returns true if found, otherwise the idx is where val should be inserted, standard bsearch doesn't do this
void isort(void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*));//binary insertion sort

typedef enum GetOptRetEnum
{
	OPT_ENDOFARGS=-3,
	OPT_INVALIDARG,
	OPT_NOMATCH,
} GetOptRet;
int acme_getopt(int argc, char **argv, int *start, const char **keywords, int kw_count);//keywords[i]: shortform char, followed by longform null-terminated string, returns 

int floor_log2(unsigned long long n);
int ceil_log2(unsigned long long n);
int floor_log10(double x);
double power(double x, int y);
double _10pow(int n);
int minimum(int a, int b);
int maximum(int a, int b);
int acme_isdigit(char c, char base);

double time_ms();

typedef struct TimeInfoStruct
{
	int days, hours, mins;
	float secs;
} TimeInfo;
void parsetimedelta(double delta, TimeInfo *ti);
int timedelta2str(char *buf, size_t len, double ms);
int acme_strftime(char *buf, size_t len, const char *format);//prints current time to string

int print_bin8(int x);
int print_bin32(unsigned x);
int print_binn(unsigned long long x, int nbits);

//error handling
int log_error(const char *file, int line, int quit, const char *format, ...);//doesn't stop execution
#define LOG_ERROR(format, ...)   log_error(file, __LINE__, 1, format, ##__VA_ARGS__)
#define LOG_ERROR2(format, ...)  log_error(__FILE__, __LINE__, 1, format, ##__VA_ARGS__)
#define LOG_WARNING(format, ...) log_error(file, __LINE__, 0, format, ##__VA_ARGS__)
int valid(const void *p);
int pause();
#ifdef _MSC_VER
int pause1();
#endif
int pause_abort(const char *file, int lineno, const char *extraInfo);
#define PANIC() pause_abort(file, __LINE__, 0)
#define ASSERT(SUCCESS)   ((SUCCESS)!=0||pause_abort(file, __LINE__, #SUCCESS))
#define ASSERT_P(POINTER) (valid(POINTER)||pause_abort(file, __LINE__, #POINTER " == 0"))


//ARRAY
#if 1
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4200)//no default-constructor for struct with zero-length array
#endif
typedef struct ArrayHeaderStruct//32 bytes on 64 bit system, or 16 bytes on 32 bit system
{
	size_t count,
		esize, cap;//in bytes
	void (*destructor)(void*);
	unsigned char data[];
} ArrayHeader, *ArrayHandle;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
ArrayHandle array_construct(const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*));
ArrayHandle array_copy(ArrayHandle *arr);//shallow
void  array_clear(ArrayHandle *arr);//keeps allocation
void  array_free(ArrayHandle *arr);
void  array_fit(ArrayHandle *arr, size_t pad);

void* array_insert(ArrayHandle *arr, size_t idx, const void *data, size_t count, size_t rep, size_t pad);//cannot be nullptr
void* array_erase(ArrayHandle *arr, size_t idx, size_t count);//does not reallocate
void* array_replace(ArrayHandle *arr, size_t idx, size_t rem_count, const void *data, size_t ins_count, size_t rep, size_t pad);

void* array_at(ArrayHandle *arr, size_t idx);
void* array_back(ArrayHandle *arr);

int str_append(ArrayHandle *str, const char *format, ...);//requires C99, prints twice

#define ARRAY_ALLOC(ELEM_TYPE, ARR, DATA, COUNT, PAD, DESTRUCTOR) ARR=array_construct(DATA, sizeof(ELEM_TYPE), COUNT, 1, PAD, DESTRUCTOR)
#define ARRAY_APPEND(ARR, DATA, COUNT, REP, PAD)  array_insert(&(ARR), (ARR)->count, DATA, COUNT, REP, PAD)
#define ARRAY_APPEND_OFFSET(ARR, DATA, COUNT, REP, PAD) (((char*)array_insert(&(ARR), (ARR)->count, DATA, COUNT, REP, PAD)-(ARR)->data)/(ARR)->esize)
//#define ARRAY_DATA(ARR) (ARR)->data
//#define ARRAY_I(ARR, IDX) *(int*)array_at(&ARR, IDX)
//#define ARRAY_U(ARR, IDX) *(unsigned*)array_at(&ARR, IDX)
//#define ARRAY_F(ARR, IDX) *(double*)array_at(&ARR, IDX)


//null terminated array
#define ESTR_ALLOC(TYPE, STR, DATA, LEN) STR=array_construct(DATA, sizeof(TYPE), LEN, 1, 1, 0)
#define STR_APPEND(STR, SRC, LEN, REP)   array_insert(&(STR), (STR)->count, SRC, LEN, REP, 1)
#define STR_POPBACK(STR, COUNT)          memset(array_erase(&(STR), (STR)->count-(COUNT), COUNT), 0, (STR)->esize)
#define STR_FIT(STR) array_fit(&STR, 1)
#define ESTR_AT(TYPE, STR, IDX) *(TYPE*)array_at(&(STR), IDX)

#define STR_ALLOC(STR, LEN)      ESTR_ALLOC(char, STR, 0, LEN)
#define STR_COPY(STR, DATA, LEN) ESTR_ALLOC(char, STR, DATA, LEN)
#define STR_AT(STR, IDX)         ESTR_AT(char, STR, IDX)

#define WSTR_ALLOC(STR, LEN)      ESTR_ALLOC(wchar_t, STR, 0, LEN)
#define WSTR_COPY(STR, DATA, LEN) ESTR_ALLOC(wchar_t, STR, DATA, LEN)
#define WSTR_AT(STR, IDX)         ESTR_AT(wchar_t, STR, IDX)
#endif


//double-linked LIST of identical size arrays,		append-only, no mid-insertion
#if 1
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4200)//no default-constructor for struct with zero-length array
#endif
typedef struct DNodeStruct
{
	struct DNodeStruct *prev, *next;
	unsigned char data[];
} DNodeHeader, *DNodeHandle;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
typedef struct DListStruct
{
	DNodeHandle i, f;
	size_t
		objsize,	//size of one contained object
		objpernode,	//object count per node,		recommended value 128
		nnodes,		//node count
		nobj;		//total object count
	void (*destructor)(void*);
} DList, *DListHandle;
void dlist_init(DListHandle list, size_t objsize, size_t objpernode, void (*destructor)(void*));
void dlist_copy(DListHandle dst, DListHandle src);
void dlist_clear(DListHandle list);
void dlist_appendtoarray(DListHandle list, ArrayHandle *dst);
//void dlist_appendtoarrayandclear(DListHandle list, ArrayHandle *dst);

void* dlist_push_back1(DListHandle list, const void *obj);//shallow copy of obj
void* dlist_push_back(DListHandle list, const void *data, size_t count);
void* dlist_back(DListHandle list);//returns address of last object
void  dlist_pop_back1(DListHandle list);
void  dlist_pop_back(DListHandle list, size_t count);

//iterator: seamlessly iterate through contained objects
typedef struct DListIteratorStruct
{
	DListHandle list;
	DNodeHandle node;
	size_t obj_idx;
} DListIterator, *DListItHandle;
void  dlist_first(DListHandle list, DListItHandle it);
void  dlist_last(DListHandle list, DListItHandle it);
void* dlist_it_deref(DListItHandle it);
int   dlist_it_inc(DListItHandle it);
int   dlist_it_dec(DListItHandle it);
#endif


//ordered MAP (implemented as a (self-balancing) red-black tree)
#if 1
typedef struct RBNodeStruct
{
	struct RBNodeStruct *parent, *left, *right;
	size_t is_red;
	unsigned char data[];//key then value
} RBNodeHeader, *RBNodeHandle;
typedef enum CmpResEnum
{
	RESULT_LESS=-1,
	RESULT_EQUAL,
	RESULT_GREATER,
} CmpRes;
typedef CmpRes (*MapCmpFn)(const void *key, const void *candidate);//the search key is always on left
typedef struct MapStruct
{
	size_t
		esize,	//object size in bytes
		nnodes;	//object count
	RBNodeHandle root;
	MapCmpFn comparator;
	void (*destructor)(void*);//key and value are packed consequtively
} Map, *MapHandle;
typedef Map const *MapConstHandle;
void          map_init(MapHandle map, size_t esize, MapCmpFn comparator, void (*destructor)(void*));
RBNodeHandle* map_find(MapHandle map, const void *key);
RBNodeHandle* map_insert(MapHandle map, const void *data, int *found);//the map doesn't know where the object->key member(s) is/are, initialize entire object yourself, including the passed key
int           map_erase(MapHandle map, const void *data, RBNodeHandle node);//either pass data object or node

void map_clear_r(MapHandle map, RBNodeHandle node);
void map_debugprint_r(RBNodeHandle *node, int depth, void (*printer)(RBNodeHandle *node, int depth));

#define MAP_INIT(MAP, ETYPE, CMP, DESTRUCTOR) map_init(MAP, sizeof(ETYPE), CMP, DESTRUCTOR)
#define MAP_ERASE_DATA(MAP, DATA)             map_erase(MAP, DATA, 0)
#define MAP_ERASE_NODE(MAP, NODE)             map_erase(MAP, 0, NODE)
#define MAP_CLEAR(MAP)                        map_clear_r(MAP, (MAP)->root), (MAP)->root=0, (MAP)->nnodes=0
#define MAP_DEBUGPRINT(MAP, PRINTER)          map_debugprint_r(&(MAP)->root, 0, PRINTER)
#endif


//single-linked list, queue and stack
#if 1
typedef struct SNodeStruct
{
	struct SNodeStruct *prev;
	unsigned char data[];//4-byte aligned on 32-bit, not suitable for double on 32-bit
} SNode, *SNodeHandle;
typedef struct SListStruct
{
	//[front] -> ... -> [back] -> nullptr
	size_t esize, count;
	void (*destructor)(void*);
	SNodeHandle
		front,	//can remove from or append to front
		back;	//prev always nullptr, can only append to back
} SList, *SListHandle;

//API
void slist_init(SListHandle list, size_t esize, void (*destructor)(void*));
void slist_clear(SListHandle list);
void* slist_push_front(SListHandle list, const void *data);
void* slist_push_back(SListHandle list, const void *data);
void* slist_front(SListHandle list);
void* slist_back(SListHandle list);
void slist_pop_front(SListHandle list);
void slist_print(SListHandle list, void (*printer)(const void*));

//list-based stack
#define STACK_PUSH(LIST, DATA) slist_push_front(LIST, DATA)
#define STACK_TOP(LIST)        slist_front(LIST)
#define STACK_POP(LIST)        slist_pop_front(LIST)

//list-based queue
#define QUEUE_ENQUEUE(LIST, DATA) slist_push_back(LIST, DATA)
#define QUEUE_FRONT(LIST)         slist_front(LIST)
#define QUEUE_DEQUEUE(LIST)       slist_pop_front(LIST)
#endif


//bit-string
#if 1
typedef struct BitstringStruct
{
	size_t bitCount, byteCap;
	unsigned char data[];//8bits/element
	//unsigned data[];//32bits/element
} BitstringHeader, *BitstringHandle;
//-> [bitCount], [byteCap], data...
//arr->data		((char*)arr+sizeof(Header))

//bitstring_construct: allocates a new bit-string
//src: If 0, initialize memory with 0
//bitCount:  The number of bits in resulting array
//bytePad:  Number of zero-pad bytes
BitstringHandle bitstring_construct(const void *src, size_t bitCount, size_t bitOffset, size_t bytePad);

//bitstring_free: Frees allocated buffer and sets pointer to zero
//str: The bit-string to be freed
void bitstring_free(BitstringHandle *str);

//bitstring_append: Appends bits to bit-string
//str: The destination bit-string
//src: A pointer to the source data (bytes)
//bitCount:  The number of bits to read from src
//  Actually, this will read (bitCount+7)/8 bytes
//bitOffset: The number of bits to skip at the start of src. Must be between 0 and 7.
//For example, if s1 & s2 are bit-strings, the following should work properly:
//  bitstring_append(&s2, s1->data, s1->bitCount, 0);
void bitstring_append(BitstringHandle *str, const void *src, size_t bitCount, size_t bitOffset);

//bitstring_get: Get bit from bit-string
//str: The source bit-string
//bitIdx:  The index of bit to get
int bitstring_get(BitstringHandle *str, size_t bitIdx);

//bitstring_set: Set a certain bit in bit-string. Does not resize the buffer and calls some error function on out-of-bounds access.
//str: The destination bit-string
//bitIdx:  The index of bit to set
//bit: The bit value to set
void bitstring_set(BitstringHandle *str, size_t bitIdx, int bit);

//bitstring_print: Prints the bits in the bit-string
void bitstring_print(BitstringHandle str);
#endif


//Max-heap-based priority queue
#if 1
typedef struct PQueueStruct
{
	size_t count, //number of elements
	esize, //total element size in bytes
	byteCap;  //allocated buffer size in bytes
	int (*less)(const void*, const void*);//comparator function pointer
	void (*destructor)(void*);//element destructor, can be 0
	unsigned char data[];//each contained object begins with the key
} PQueueHeader, *PQueueHandle;

//pqueue_construct: allocates a new bit-string
//src:		If 0, initialize elements with 0
//esize:	Element size
//pad:		Initial array size
//less:		The comparator, returns 1 if the left argument is less
//destructor: Optional element destructor
PQueueHandle pqueue_construct(
	size_t esize,
	size_t pad,
	int (*less)(const void*, const void*),
	void (*destructor)(void*)
);
#define PQUEUE_ALLOC(TYPE, Q, PAD, LESS, DESTRUCTOR)	Q=pqueue_construct(sizeof(TYPE), PAD, LESS, DESTRUCTOR)

//pqueue_free: Frees allocated buffer and sets pointer to zero
void pqueue_free(PQueueHandle *pq);

//pqueue_enqueue:  Inserts one element to the queue
//pq: The destination priority queue
//src: Pointer to the source data
void pqueue_enqueue(PQueueHandle *pq, const void *src);

//pqueue_front: Returns address of the maximum element in priority queue
void* pqueue_front(PQueueHandle *pq);

//pqueue_dequeue: Removes heap root from priority queue
//Return type is void
void pqueue_dequeue(PQueueHandle *pq);

//pqueue_print:	Prints the contents of the heap array given an element printing function
void pqueue_print(PQueueHandle *pq, void (*printer)(const void*));

//pqueue_print_heap: Prints contents as heap
void pqueue_print_heap(PQueueHandle *pq, void (*printer)(const void*));
#endif


ptrdiff_t get_filesize(const char *filename);//-1 not found,  0: folder (probably),  ...: regular file size

int acme_stricmp(const char *a, const char *b);//case insensitive strcmp
ptrdiff_t acme_strrchr(const char *str, ptrdiff_t len, char c);//find last occurrence, with known length for backward search
ArrayHandle filter_path(const char *path);//replaces back slashes with slashes, and adds trailing slash if missing, as ArrayHandle
ArrayHandle get_filenames(const char *path, const char **extensions, int extCount, int fullyqualified);//returns array of strings

ArrayHandle load_file(const char *filename, int bin, int pad, int erroronfail);
int save_file_bin(const char *filename, const unsigned char *src, size_t srcSize);

	
#ifdef __cplusplus
}
#endif
#endif//INC_UTIL_H
