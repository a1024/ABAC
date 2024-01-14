//awm_util.h - Utilities declarations
//Copyright (C) 2022  Ayman Wagih Mohsen
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
#ifndef AWM_UTIL_H
#define AWM_UTIL_H
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stddef.h>//for size_t
#ifdef __cplusplus
extern "C"
{
#endif

//utility
#define			COUNTOF(ARR)		(sizeof(ARR)/sizeof(*(ARR)))		//stdlib defines _countof
#define			BETWEEN(LO, X, HI)	((unsigned)((X)-LO)<(unsigned)(HI+1-LO))
#ifndef _MSC_VER
#define			sprintf_s	snprintf
#endif
#define			G_BUF_SIZE	4096
extern char		g_buf[G_BUF_SIZE];

void			memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes);
void			memswap_slow(void *p1, void *p2, size_t size);
void 			memswap(void *p1, void *p2, size_t size, void *temp);
void			memreverse(void *p, size_t count, size_t esize);//calls memswap
void 			memrotate(void *p, size_t byteoffset, size_t bytesize, void *temp);//temp buffer is min(byteoffset, bytesize-byteoffset)
int 			binary_search(const void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*), const void *val, size_t *idx);//returns true if found, otherwise the idx is where val should be inserted, standard bsearch doesn't do this
void 			isort(void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*));//binary insertion sort

typedef enum GetOptRetEnum
{
	OPT_ENDOFARGS=-3,
	OPT_INVALIDARG,
	OPT_NOMATCH,
} GetOptRet;
int				acme_getopt(int argc, char **argv, int *start, const char **keywords, int kw_count);//keywords[i]: shortform char, followed by longform null-terminated string, returns 

int				floor_log2(unsigned long long n);
int				ceil_log2(unsigned long long n);
int				floor_log10(double x);
double			power(double x, int y);
double			_10pow(int n);
int				minimum(int a, int b);
int				maximum(int a, int b);
int				acme_isdigit(char c, char base);
double			time_ms();

//error handling
int				log_error(const char *file, int line, const char *format, ...);//doesn't stop execution
#define			LOG_ERROR(format, ...)	log_error(file, __LINE__, format, ##__VA_ARGS__)
int				valid(const void *p);
void			pause();
int				pause_abort(const char *file, int lineno, const char *extraInfo);
#define			PANIC()					pause_abort(file, __LINE__, 0)
#define			ASSERT(SUCCESS)			((SUCCESS)!=0||pause_abort(file, __LINE__, #SUCCESS))
#define			ASSERT_P(POINTER)		(valid(POINTER)||pause_abort(file, __LINE__, #POINTER " == 0"))


//ARRAY
#if 1
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4200)//no default-constructor for struct with zero-length array
#endif
typedef struct ArrayHeaderStruct
{
	size_t count, esize, cap;//cap is in bytes
	void (*destructor)(void*);
	unsigned char data[];
} ArrayHeader, *ArrayHandle;
//typedef const ArrayHeader *ArrayConstHandle;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
ArrayHandle		array_construct(const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*));
ArrayHandle		array_copy(ArrayHandle *arr);//shallow
void			array_clear(ArrayHandle *arr);//keeps allocation
void			array_free(ArrayHandle *arr);
void			array_fit(ArrayHandle *arr, size_t pad);

void*			array_insert(ArrayHandle *arr, size_t idx, const void *data, size_t count, size_t rep, size_t pad);//cannot be nullptr
void*			array_erase(ArrayHandle *arr, size_t idx, size_t count);
void*			array_replace(ArrayHandle *arr, size_t idx, size_t rem_count, const void *data, size_t ins_count, size_t rep, size_t pad);

void*			array_at(ArrayHandle *arr, size_t idx);
void*			array_back(ArrayHandle *arr);

#define			ARRAY_ALLOC(ELEM_TYPE, ARR, DATA, COUNT, PAD, DESTRUCTOR)	ARR=array_construct(DATA, sizeof(ELEM_TYPE), COUNT, 1, PAD, DESTRUCTOR)
#define			ARRAY_APPEND(ARR, DATA, COUNT, REP, PAD)					array_insert(&(ARR), (ARR)->count, DATA, COUNT, REP, PAD)
#define			ARRAY_DATA(ARR)			(ARR)->data
#define			ARRAY_I(ARR, IDX)		*(int*)array_at(&ARR, IDX)
#define			ARRAY_U(ARR, IDX)		*(unsigned*)array_at(&ARR, IDX)
#define			ARRAY_F(ARR, IDX)		*(double*)array_at(&ARR, IDX)


//null terminated array
#define			ESTR_ALLOC(TYPE, STR, DATA, LEN)	STR=array_construct(DATA, sizeof(TYPE), LEN, 1, 1, 0)
#define			STR_APPEND(STR, SRC, LEN, REP)		array_insert(&(STR), (STR)->count, SRC, LEN, REP, 1)
#define			STR_FIT(STR)						array_fit(&STR, 1)
#define			ESTR_AT(TYPE, STR, IDX)				*(TYPE*)array_at(&(STR), IDX)

#define			STR_ALLOC(STR, LEN)				ESTR_ALLOC(char, STR, 0, LEN)
#define			STR_COPY(STR, DATA, LEN)		ESTR_ALLOC(char, STR, DATA, LEN)
#define			STR_AT(STR, IDX)				ESTR_AT(char, STR, IDX)

#define			WSTR_ALLOC(STR, LEN)			ESTR_ALLOC(wchar_t, STR, 0, LEN)
#define			WSTR_COPY(STR, DATA, LEN)		ESTR_ALLOC(wchar_t, STR, DATA, LEN)
#define			WSTR_AT(STR, IDX)				ESTR_AT(wchar_t, STR, IDX)
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
void			dlist_init(DListHandle list, size_t objsize, size_t objpernode, void (*destructor)(void*));
void			dlist_copy(DListHandle dst, DListHandle src);
void			dlist_clear(DListHandle list);
void			dlist_appendtoarray(DListHandle list, ArrayHandle *dst);
void			dlist_appendtoarrayandclear(DListHandle list, ArrayHandle *dst);

void*			dlist_push_back1(DListHandle list, const void *obj);//shallow copy of obj
void*			dlist_push_back(DListHandle list, const void *data, size_t count);
void*			dlist_back(DListHandle list);//returns address of last object
void			dlist_pop_back(DListHandle list);

//iterator: seamlessly iterate through contained objects
typedef struct DListIteratorStruct
{
	DListHandle list;
	DNodeHandle node;
	size_t obj_idx;
} DListIterator, *DListItHandle;
void			dlist_first(DListHandle list, DListItHandle it);
void			dlist_last(DListHandle list, DListItHandle it);
void*			dlist_it_deref(DListItHandle it);
int				dlist_it_inc(DListItHandle it);
int				dlist_it_dec(DListItHandle it);
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
void			map_init(MapHandle map, size_t esize, MapCmpFn comparator, void (*destructor)(void*));
RBNodeHandle*	map_find(MapHandle map, const void *key);
RBNodeHandle*	map_insert(MapHandle map, const void *data, int *found);//the map doesn't know where the object->key member(s) is/are, initialize entire object yourself, including the passed key
int				map_erase(MapHandle map, const void *data, RBNodeHandle node);//either pass data object or node

void			map_clear_r(MapHandle map, RBNodeHandle node);
void			map_debugprint_r(RBNodeHandle *node, int depth, void (*printer)(RBNodeHandle *node, int depth));

#define			MAP_INIT(MAP, ETYPE, CMP, DESTRUCTOR)	map_init(MAP, sizeof(ETYPE), CMP, DESTRUCTOR)
#define			MAP_ERASE_DATA(MAP, DATA)				map_erase(MAP, DATA, 0)
#define			MAP_ERASE_NODE(MAP, NODE)				map_erase(MAP, 0, NODE)
#define			MAP_CLEAR(MAP)							map_clear_r(MAP, (MAP)->root), (MAP)->root=0, (MAP)->nnodes=0
#define			MAP_DEBUGPRINT(MAP, PRINTER)			map_debugprint_r(&(MAP)->root, 0, PRINTER)
#endif

	
#ifdef __cplusplus
}
#endif
#endif//AWM_UTIL_H