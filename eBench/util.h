#pragma once
#ifndef INC_UTIL_H
#define INC_UTIL_H
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stddef.h>
#include<stdint.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __cplusplus
extern "C"
{
#endif

/*utility*/
#ifdef _MSC_VER
#	define	ALIGN(N) __declspec(align(N))
#	define INLINE static
#	define FORCE_INLINE __forceinline static
#	if 0
#		define sprintf_s snprintf
#	endif
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define INLINE static inline
#	define FORCE_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#	if 0
/*moved to source	because this interferes with later includes*/
#		define _stricmp strcasecmp
#	endif
#endif

#ifndef INFINITY/*C89*/
#include<stdio.h>
#include<stdarg.h>
#include<float.h>
#define INFINITY (float)DBL_MAX
#define round(x) floor((x)+0.5)
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif
static int acme_snprintf(char *buf, size_t len, const char *format, ...)
{
	va_list args;
	int ret;

	va_start(args, format);
	ret=vsprintf(buf, format, args);
	va_end(args);
	(void)len;
	return ret;
}
static int acme_isinf(double x)
{
	return fabs(x)<=DBL_MAX;
}
static int acme_isinff(float x)
{
	return fabsf(x)<=FLT_MAX;
}
#ifdef _MSC_VER
#	define acme_isnan _isnan
#	define inline __inline
#elif defined __GNUC__
#	define acme_isnan __builtin_isnan
#	define inline __inline__
#endif
#if 0
#undef  FORCE_INLINE
#define FORCE_INLINE static
#endif
#else
#define acme_snprintf snprintf
#	define acme_isnan isnan
#	define acme_isinf isinf
#endif

#define BETWEEN_INC(LO, X, HI) ((unsigned)((X)-LO)<(unsigned)(HI+1-LO))
#define BETWEEN_EXC(LO, X, HI) ((unsigned)((X)-LO)<(unsigned)(HI-LO))
#define SWAPVAR(A, B, TEMP) TEMP=A, A=B, B=TEMP
#define SWAPMEM(A, B, TEMP) memcpy(TEMP, A, sizeof(*(TEMP))), memcpy(A, B, sizeof(*(TEMP))), memcpy(B, TEMP, sizeof(*(TEMP)))
#define ROTATE3(A, B, C, TEMP) TEMP=A, A=B, B=C, C=TEMP
#define MINVAR(A, B) ((A)<(B)?(A):(B))/*_mm_min/max are faster*/
#define MAXVAR(A, B) ((A)>(B)?(A):(B))
#define CLAMP2(X, L, H) (X=X<(L)?L:X, X=X>(H)?H:X)
#define CLAMP(L, X, H) ((X)>(L)?(X)<(H)?(X):(H):(L))/*use MEDIAN3_32*/
#ifdef _MSC_VER
#	define LZCNT32 _lzcnt_u32
#	define LZCNT64 _lzcnt_u64
#	define TZCNT32 _tzcnt_u32
#	define TZCNT64 _tzcnt_u64
#else
/*zero is UB!*/
#	define LZCNT32 __builtin_clz
#	define LZCNT64 __builtin_clzll
#	define TZCNT32 __builtin_ctz
#	define TZCNT64 __builtin_ctzll
#endif

/*<smmintrin.h>	SSE4.1*/
#define MEDIAN3_32(DST, A, B, C)\
	do\
	{\
		ALIGN(16) int _arr[4];\
		__m128i _va=_mm_set_epi32(0, 0, 0, A);\
		__m128i _vb=_mm_set_epi32(0, 0, 0, B);\
		__m128i _vc=_mm_set_epi32(0, 0, 0, C);\
		__m128i _vmin=_mm_min_epi32(_va, _vb);\
		__m128i _vmax=_mm_max_epi32(_va, _vb);\
		_vc=_mm_max_epi32(_vc, _vmin);\
		_vc=_mm_min_epi32(_vc, _vmax);\
		_mm_store_si128((__m128i*)_arr, _vc);\
		DST=_arr[0];\
	}while(0)

/*include<emmintrin.h>	SSE2*/
#define MEDIAN3_16(DST, A, B, C)\
	do\
	{\
		ALIGN(16) short arr[8];\
		__m128i va=_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, A);\
		__m128i vb=_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, B);\
		__m128i vc=_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, C);\
		__m128i vmin=_mm_min_epi16(va, vb);\
		__m128i vmax=_mm_max_epi16(va, vb);\
		vc=_mm_max_epi16(vc, vmin);\
		vc=_mm_min_epi16(vc, vmax);\
		_mm_store_si128((__m128i*)arr, vc);\
		DST=arr[0];\
	}while(0)
#define MEDIAN3_fp64(DST, A, B, C)\
	do\
	{\
		ALIGN(16) double arr[2];\
		__m128d va=_mm_set_pd(0, A);\
		__m128d vb=_mm_set_pd(0, B);\
		__m128d vc=_mm_set_pd(0, C);\
		__m128d vmin=_mm_min_pd(va, vb);\
		__m128d vmax=_mm_max_pd(va, vb);\
		vc=_mm_max_pd(vc, vmin);\
		vc=_mm_min_pd(vc, vmax);\
		_mm_store_pd(arr, vc);\
		DST=arr[0];\
	}while(0)
#define MEDIAN3_CMOV(DST, A, B, C)\
	do\
	{\
		int _vmax=A, _vmin=B;\
		if((A)<(B))_vmin=A, _vmax=B;\
		DST=C;\
		CLAMP2(DST, _vmin, _vmax);\
	}while(0)
#define MEDIAN3(A, B, C) (B<A?B<C?C<A?C:A:B:A<C?C<B?C:B:A)/*SLOW*/
	
/*include<smmintrin.h>	SSE4.1*/
#define CLAMP2_32(DST, X, LO, HI)\
	do\
	{\
		ALIGN(16) int _dst[4];\
		__m128i _mx=_mm_set_epi32(0, 0, 0, X);\
		__m128i _vmin=_mm_set_epi32(0, 0, 0, LO);\
		__m128i _vmax=_mm_set_epi32(0, 0, 0, HI);\
		_mx=_mm_max_epi32(_mx, _vmin);\
		_mx=_mm_min_epi32(_mx, _vmax);\
		_mm_store_si128((__m128i*)_dst, _mx);\
		DST=_dst[0];\
	}while(0)
#define CLAMP2_fp64(DST, X, LO, HI)\
	do\
	{\
		ALIGN(16) double _dst[2];\
		__m128d _mx=_mm_set_pd(0, X);\
		__m128d _vmin=_mm_set_pd(0, LO);\
		__m128d _vmax=_mm_set_pd(0, HI);\
		_mx=_mm_max_pd(_mx, _vmin);\
		_mx=_mm_min_pd(_mx, _vmax);\
		_mm_store_pd(_dst, _mx);\
		DST=_dst[0];\
	}while(0)
#define CLAMP3_32(DST, X, A, B, C)\
	do\
	{\
		ALIGN(16) int _dst[4];\
		__m128i _mx=_mm_set_epi32(0, 0, 0, X);\
		__m128i ma=_mm_set_epi32(0, 0, 0, A);\
		__m128i mb=_mm_set_epi32(0, 0, 0, B);\
		__m128i mc=_mm_set_epi32(0, 0, 0, C);\
		__m128i vmin=_mm_min_epi32(ma, mb);\
		__m128i vmax=_mm_max_epi32(ma, mb);\
		vmin=_mm_min_epi32(vmin, mc);\
		vmax=_mm_max_epi32(vmax, mc);\
		_mx=_mm_max_epi32(_mx, vmin);\
		_mx=_mm_min_epi32(_mx, vmax);\
		_mm_store_si128((__m128i*)_dst, _mx);\
		DST=_dst[0];\
	}while(0)

#define RCPSS(X) _mm_cvtss_f32(_mm_rcp_ss(_mm_set_ss(X)))/*RCPSS is inaccurate across platforms*/
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))

#define IDIV23(DST, NUM, DEN)\
	do\
	{\
		__m128 mnum=_mm_set_ss((float)NUM);\
		__m128 mden=_mm_set_ss((float)DEN);\
		mnum=_mm_div_ss(mnum, mden);\
		__m128i result=_mm_cvtps_epi32(mnum);\
		DST=_mm_extract_epi32(result, 0);\
	}while(0)

#define IDIV52(DST, NUM, DEN)\
	do\
	{\
		__m128d mnum=_mm_set_sd((double)NUM);\
		__m128d mden=_mm_set_sd((double)DEN);\
		mnum=_mm_div_sd(mnum, mden);\
		__m128i result=_mm_cvtpd_epi32(mnum);\
		DST=_mm_extract_epi32(result, 0);\
	}while(0)

#define MOVEOBJ(SRC, DST, SIZE) memcpy(DST, SRC, SIZE), memset(SRC, 0, SIZE)
#define MODVAR(DST, SRC, N) (DST=(SRC)%(N), DST+=(N)&-(DST<0))
#define SHIFT_LEFT_SIGNED(X, SH) ((SH)<0?(X)>>-(SH):(X)<<(SH))
#define SHIFT_RIGHT_SIGNED(X, SH) ((SH)<0?(X)<<-(SH):(X)>>(SH))
#define UPDATE_MIN(M, X) if(M>X)M=X
#define UPDATE_MAX(M, X) if(M<X)M=X
#define THREEWAY(L, R) (((L)>(R))-((L)<(R)))
#define MIX(V0, V1, X) ((V0)+((V1)-(V0))*(X))
#define LSB_IDX_64(X)	(int)_tzcnt_u64(X)
#define LSB_IDX_32(X)	(int)_tzcnt_u32(X)
#define LSB_IDX_16(X)	(int)_tzcnt_u16(X)

enum
{
	G_BUF_SIZE=4096
};
extern char g_buf[G_BUF_SIZE];

void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes);
#define FILLMEM(PTR, ELEM, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(ELEM);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
void memswap_slow(void *p1, void *p2, size_t size);
void memswap(void *p1, void *p2, size_t size, void *temp);
void memreverse(void *p, size_t count, size_t esize);/*calls memswap*/
void reverse16(void *start, void *end);
void memrotate(void *p, size_t byteoffset, size_t bytesize, void *temp);/*temp buffer is min(byteoffset, bytesize-byteoffset)*/

/*returns true if found, otherwise the idx is where val should be inserted, standard bsearch doesn't do this*/
int binary_search(const void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*), const void *val, size_t *idx);

/*binary insertion sort*/
void isort(void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*));

typedef enum GetOptRetEnum
{
	OPT_ENDOFARGS=-3,
	OPT_INVALIDARG,
	OPT_NOMATCH
} GetOptRet;
/*keywords[i]: shortform char, followed by longform null-terminated string, returns FIXME*/
int acme_getopt(int argc, char **argv, int *start, const char **keywords, int kw_count);

int hammingweight16(unsigned short x);
int hammingweight32(unsigned x);
int hammingweight64(uint64_t x);
#if 0
int floor_log2_p1(uint64_t n);
int floor_log2(uint64_t n);		/*use (31^_lzcnt_u64(n)) instead*/
int floor_log2_32(unsigned n);			/*use (31^_lzcnt_u32(n)) instead*/
int ceil_log2(uint64_t n);
int ceil_log2_32(unsigned n);
int get_lsb_index(uint64_t n);	/*returns lsb position + 1,  returns register bit count if n is zero*/
int get_lsb_index32(unsigned n);
int get_lsb_index16(unsigned short n);
#endif
int floor_log10(double x);
unsigned floor_sqrt(uint64_t x);
unsigned exp2_fix24_neg(unsigned x);
unsigned exp2_neg_fix24_avx2(unsigned x);
uint64_t exp2_fix24(int x);
int log2_fix24(uint64_t x);
#define POW_FIX24(BASE, EXP) exp2_fix24((int)((long long)(EXP)*log2_fix24(BASE)>>24))
double power(double x, int y);
double _10pow(int n);
int acme_isdigit(char c, char base);

double time_ms(void);
double time_sec(void);

typedef struct TimeInfoStruct
{
	int days, hours, mins;
	float secs;
} TimeInfo;
void parsetimedelta(double secs, TimeInfo *ti);
int timedelta2str(char *buf, size_t len, double secs);

/*prints current time to string*/
int acme_strftime(char *buf, size_t len, const char *format);

int print_bin8(int x);
int print_bin32(unsigned x);
int print_binn(uint64_t x, int nbits);

double convert_size(double bytesize, int *log1024);
int print_size(double bytesize, int ndigits, int pdigits, char *str, int len);

/*error handling*/
int crash(const char *file, int line, int quit, const char *format, ...);
#define CRASH(LABEL) crash(__FILE__, __LINE__, 1, "%s", LABEL)

/*doesn't stop execution*/
#define WARNING(LABEL) crash(file, __LINE__, 0, "%s", LABEL)
#define ASSERT_MSG(COND) ((COND)!=0||crash(file, __LINE__, 1, "%s", #COND))

/*int valid(const void *p);*/
int pause(void);
#if defined _MSC_VER && 0
int pause1(void);
#endif
int pause_abort(const char *file, int lineno, const char *extraInfo);
#define PANIC() pause_abort(file, __LINE__, 0)
#define ASSERT(SUCCESS) ((SUCCESS)!=0||pause_abort(file, __LINE__, #SUCCESS))
#if 0
#define ASSERT_P(POINTER) (void)(valid(POINTER)||pause_abort(file, __LINE__, #POINTER " == 0"))
#endif


#define HDR_DATA(ARROBJ) ((uint8_t*)((ARROBJ)+1))
#define HDR_DATC(ARROBJ) ((const uint8_t*)((ARROBJ)+1))
#if 0
#define HDR_DATA(ARROBJ) ((uint8_t*)(ARROBJ)+sizeof(*(ARROBJ)))
#define HDR_DATC(ARROBJ) ((const uint8_t*)(ARROBJ)+sizeof(*(ARROBJ)))
#endif


/*ARRAY*/
#if 1
/*
32 bytes on 64 bit system, or 16 bytes on 32 bit system
esize, cap: in bytes
*/
typedef struct ArrayHeaderStruct
{
	size_t count,
		esize, cap;
	void (*destructor)(void*);
	/*unsigned char data[];*/
} ArrayHeader, *ArrayHandle;
ArrayHandle array_construct(const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*));
size_t array_append(ArrayHandle *dst, const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*));/*arr can be 0, returns original array size*/
ArrayHandle array_copy(ArrayHandle *arr);/*shallow*/
void  array_clear(ArrayHandle *arr);/*keeps allocation*/
void  array_free(ArrayHandle *arr);
void  array_fit(ArrayHandle *arr, size_t pad);

void* array_insert(ArrayHandle *arr, size_t idx, const void *data, size_t count, size_t rep, size_t pad);/*cannot be nullptr*/
void* array_erase(ArrayHandle *arr, size_t idx, size_t count);/*does not reallocate*/
void* array_replace(ArrayHandle *arr, size_t idx, size_t rem_count, const void *data, size_t ins_count, size_t rep, size_t pad);

#define ARRAY_AT(ETYPE, ARR, IDX) (ETYPE*)((ARR)&&IDX<(ARR)->count?HDR_DATA(ARR)+(IDX)*(ARR)->esize:LOG_ERROR("OOB"))
void* array_at(ArrayHandle *arr, size_t idx);
void* array_back(ArrayHandle *arr);

int str_append(ArrayHandle *str, const char *format, ...);/*requires C99, calls vsnprintf twice*/

#define ARRAY_ALLOC(ELEM_TYPE, ARR, DATA, COUNT, PAD, DESTRUCTOR) ARR=array_construct(DATA, sizeof(ELEM_TYPE), COUNT, 1, PAD, DESTRUCTOR)
#define ARRAY_APPEND(ARR, DATA, COUNT, REP, PAD) array_insert(&(ARR), (ARR)->count, DATA, COUNT, REP, PAD)
#define ARRAY_APPEND_OFFSET(ARR, DATA, COUNT, REP, PAD) (((char*)array_insert(&(ARR), (ARR)->count, DATA, COUNT, REP, PAD)-HDR_DATA(ARR))/(ARR)->esize)
#define ARRAY_DATA(ARR) (void*)((uint8_t*)(ARR)+sizeof(*(ARR)))
#if 0
#define ARRAY_U32(ARR, IDX) *(uint32_t*)array_at(&ARR, IDX)
#define ARRAY_U(ARR, IDX) *(unsigned*)array_at(&ARR, IDX)
#define ARRAY_F(ARR, IDX) *(double*)array_at(&ARR, IDX)
#endif


/*null terminated array*/
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


/*double-linked LIST of identical size arrays,		append-only, no mid-insertion*/
#if 1
typedef struct DNodeStruct
{
	struct DNodeStruct *prev, *next;
	/*unsigned char data[];*/
} DNodeHeader, *DNodeHandle;
/*
objsize	size of one contained object
objpernode	object count per node,		recommended value 128
nnodes	node count
nobj	total object count
*/
typedef struct DListStruct
{
	DNodeHandle i, f;
	size_t
		objsize,
		objpernode,
		nnodes,
		nobj;
	void (*destructor)(void*);
	/*unsigned char data[];*/
} DList, *DListHandle;
void dlist_init(DListHandle list, size_t objsize, size_t objpernode, void (*destructor)(void*));
void dlist_copy(DListHandle dst, DListHandle src);
void dlist_clear(DListHandle list);
size_t dlist_appendtoarray(DListHandle list, ArrayHandle *dst);
#if 0
void dlist_appendtoarrayandclear(DListHandle list, ArrayHandle *dst);
#endif

/*shallow copy of obj*/
void* dlist_push_back1(DListHandle list, const void *obj);
void* dlist_push_back(DListHandle list, const void *data, size_t count);
/*returns address of last object*/
void* dlist_back(DListHandle list);
void  dlist_pop_back1(DListHandle list);
void  dlist_pop_back(DListHandle list, size_t count);

/*iterator: seamlessly iterate through contained objects*/
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


/*ordered MAP/SET (implemented as a (self-balancing) red-black tree)*/
#if 1
typedef struct RBNodeStruct
{
	struct RBNodeStruct *parent, *left, *right;
	size_t is_red;
	/*unsigned char data[];*/
} RBNodeHeader, *RBNodeHandle;
typedef enum CmpResEnum
{
	RESULT_LESS=-1,
	RESULT_EQUAL,
	RESULT_GREATER
} CmpRes;
typedef CmpRes (*MapCmpFn)(const void *key, const void *candidate);
/*
esize: object size in bytes
nnodes: object count
*/
typedef struct MapStruct
{
	size_t
		esize,
		nnodes;
	RBNodeHandle root;
	MapCmpFn comparator;
	void (*destructor)(void*);
	/*unsigned char data[];*/
} Map, *MapHandle;
typedef Map const *MapConstHandle;
void          map_init(MapHandle map, size_t esize, MapCmpFn comparator, void (*destructor)(void*));
RBNodeHandle* map_find(MapHandle map, const void *key);
/*the map doesn't know where the object->key member(s) is/are, initialize entire object yourself, including the passed key*/
RBNodeHandle* map_insert(MapHandle map, const void *data, int *found);
/*either pass data object or node*/
int map_erase(MapHandle map, const void *data, RBNodeHandle node);

void map_clear_r(MapHandle map, RBNodeHandle node);
void map_debugprint_r(RBNodeHandle *node, int depth, void (*printer)(RBNodeHandle *node, int depth));

#define MAP_INIT(MAP, ETYPE, CMP, DESTRUCTOR) map_init(MAP, sizeof(ETYPE), CMP, DESTRUCTOR)
#define MAP_ERASE_DATA(MAP, DATA)             map_erase(MAP, DATA, 0)
#define MAP_ERASE_NODE(MAP, NODE)             map_erase(MAP, 0, NODE)
#define MAP_CLEAR(MAP)                        map_clear_r(MAP, (MAP)->root), (MAP)->root=0, (MAP)->nnodes=0
#define MAP_DEBUGPRINT(MAP, PRINTER)          map_debugprint_r(&(MAP)->root, 0, PRINTER)
#endif


/*single-linked list, queue and stack*/
#if 1
/*4-byte aligned on 32-bit, not suitable for double on 32-bit*/
typedef struct SNodeStruct
{
	struct SNodeStruct *prev;
	/*unsigned char data[];*/
} SNode, *SNodeHandle;
/*
[front] -> ... -> [back] -> nullptr
front: can remove from or append to front
back: prev always nullptr, can only append to back
*/
typedef struct SListStruct
{
	size_t esize, count;
	void (*destructor)(void*);
	SNodeHandle
		front,
		back;
} SList, *SListHandle;

/*single-linked list API*/
void slist_init(SListHandle list, size_t esize, void (*destructor)(void*));
void slist_clear(SListHandle list);
void* slist_push_front(SListHandle list, const void *data);
void* slist_push_back(SListHandle list, const void *data);
void* slist_front(SListHandle list);
void* slist_back(SListHandle list);
void slist_pop_front(SListHandle list);
void slist_print(SListHandle list, void (*printer)(const void*));

/*list-based stack*/
#define STACK_PUSH(LIST, DATA) slist_push_front(LIST, DATA)
#define STACK_TOP(LIST)        slist_front(LIST)
#define STACK_POP(LIST)        slist_pop_front(LIST)

/*list-based queue*/
#define QUEUE_ENQUEUE(LIST, DATA) slist_push_back(LIST, DATA)
#define QUEUE_FRONT(LIST)         slist_front(LIST)
#define QUEUE_DEQUEUE(LIST)       slist_pop_front(LIST)
#endif


/*bit-string*/
#if 1
/*8 bits/element*/
typedef struct BitstringStruct
{
	size_t bitCount, byteCap;
	/*unsigned char data[];*/
} BitstringHeader, *BitstringHandle;
/*
-> [bitCount], [byteCap], data...
arr->data		((char*)arr+sizeof(Header))

bitstring_construct: allocates a new bit-string
src: If 0, initialize memory with 0
bitCount:  The number of bits in resulting array
bytePad:  Number of zero-pad bytes
*/
BitstringHandle bitstring_construct(const void *src, size_t bitCount, size_t bitOffset, size_t bytePad);

/*
bitstring_free: Frees allocated buffer and sets pointer to zero
str: The bit-string to be freed
*/
void bitstring_free(BitstringHandle *str);

/*
bitstring_append: Appends bits to bit-string
str: The destination bit-string
src: A pointer to the source data (bytes)
bitCount:  The number of bits to read from src
  Actually, this will read (bitCount+7)/8 bytes
bitOffset: The number of bits to skip at the start of src. Must be between 0 and 7.
For example, if s1 & s2 are bit-strings, the following should work properly:
  bitstring_append(&s2, s1->data, s1->bitCount, 0);
*/
void bitstring_append(BitstringHandle *str, const void *src, size_t bitCount, size_t bitOffset);

/*
bitstring_get: Get bit from bit-string
str: The source bit-string
bitIdx:  The index of bit to get
*/
int bitstring_get(BitstringHandle *str, size_t bitIdx);

/*
bitstring_set: Set a certain bit in bit-string. Does not resize the buffer and calls some error function on out-of-bounds access.
str: The destination bit-string
bitIdx:  The index of bit to set
bit: The bit value to set
*/
void bitstring_set(BitstringHandle *str, size_t bitIdx, int bit);

void bitstring_print(BitstringHandle str);/*bitstring_print: Prints the bits in the bit-string*/
#endif


/*Priority Queue (Max-heap-based)*/
#if 1
/*
count: number of elements
esize: total element size in bytes
buteCap: allocated buffer size in bytes
less: comparator function pointer
destructor: element destructor, can be 0
*/
typedef struct PQueueStruct
{
	size_t count,
		esize,
		byteCap;
	int (*less)(const void*, const void*);
	void (*destructor)(void*);
} PQueueHeader, *PQueueHandle;

/*pqueue_construct: allocates a new bit-string
src:	If 0, initialize elements with 0
esize:	Element size
pad:	Initial array size
less:	The comparator, returns 1 if the left argument is less
destructor: Optional element destructor
*/
PQueueHandle pqueue_construct(
	size_t esize,
	size_t pad,
	int (*less)(const void*, const void*),
	void (*destructor)(void*)
);
#define PQUEUE_ALLOC(TYPE, Q, PAD, LESS, DESTRUCTOR)	Q=pqueue_construct(sizeof(TYPE), PAD, LESS, DESTRUCTOR)

/*pqueue_free: Frees allocated buffer and sets pointer to zero*/
void pqueue_free(PQueueHandle *pq);

/*
pqueue_enqueue:  Inserts one element to the queue
pq: The destination priority queue
src: Pointer to the source data
*/
void pqueue_enqueue(PQueueHandle *pq, const void *src);

/*pqueue_front: Returns address of the maximum element in priority queue*/
void* pqueue_front(PQueueHandle *pq);

/*
pqueue_dequeue: Removes heap root from priority queue
Return type is void
*/
void pqueue_dequeue(PQueueHandle *pq);

/*pqueue_print:	Prints the contents of the heap array given an element printing function*/
void pqueue_print(PQueueHandle *pq, void (*printer)(const void*));

/*pqueue_print_heap: Prints contents as heap*/
void pqueue_print_heap(PQueueHandle *pq, void (*printer)(const void*));
#endif


/*
returns:
-1	not found
0	folder (probably)
...	regular file size
*/
ptrdiff_t get_filesize(const char *filename);

/*case insensitive strcmp*/
int acme_stricmp(const char *a, const char *b);

/*find last occurrence, with known length for backward search*/
ptrdiff_t acme_strrchr(const char *str, ptrdiff_t len, char c);

/*replaces back slashes with slashes, and adds trailing slash if missing, as ArrayHandle*/
ArrayHandle filter_path(const char *path, int len, int appendslash);

/*pass -1 for len if unknown*/
void get_filetitle(const char *fn, int len, int *idx_start, int *idx_end);
void free_str(void *p);

/*returns array of strings, extensions without period '.'*/
ArrayHandle get_filenames(const char *path, const char **extensions, int extCount, int fullyqualified);

ArrayHandle load_file(const char *filename, int bin, int pad, int erroronfail);
int save_file(const char *filename, const unsigned char *src, size_t srcSize, int is_bin);

ArrayHandle searchfor_file(const char *searchpath, const char *filetitle);

/*
returns:
0  old CPU
1  AVX2
3  AVX512
*/
int get_cpu_features(void);
int query_cpu_cores(void);


#if defined _MSC_VER && 0
#pragma warning(pop)
#endif
#ifdef __cplusplus
}
#endif
#endif/*INC_UTIL_H*/
