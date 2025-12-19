#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<ctype.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include<immintrin.h>
#include<sys/stat.h>
#if defined _WIN32 || defined WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#include<Psapi.h>
#include<tlhelp32.h>
#endif


//GDCC score
#define CALCSCORE(CSIZE, ENC, DEC) ((CSIZE)*(1./(1024*1024))+(ENC)+(DEC)*2)


static char g_buf[8192]={0}, g_buf2[8192]={0};

//runtime
#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define AWM_INLINE __forceinline static
#	if _MSC_VER<1900
#		define snprintf sprintf_s
#	endif
#else
#	define ALIGN(N) __attribute__((aligned(N)))
#	define AWM_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#define SWAPVAR(A, B, TEMP) TEMP=A, A=B, B=TEMP
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
static double time_sec(void)//time since first call
{
#if defined _WIN32 || defined WIN32
	static int64_t t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);//computer up time
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}

//array
#if 1
#define ARRAY_DECL(TYPENAME_ARRAY, TYPENAME_ELEMENT)\
	typedef struct _##TYPENAME_ARRAY\
	{\
		ptrdiff_t count, cap, esize;\
		void (*destructor)(void*);\
		TYPENAME_ELEMENT data[];\
	} *TYPENAME_ARRAY

//C-only (why would anyone use this in C++?)		SET COUNT MANUALLY
#define ARRAY_ALLOC(ARRAYNAME, PAD, DESTRUCTOR, ERR_RET, ZEROMEM)\
	do\
	{\
		ptrdiff_t _cap=(PAD)*sizeof(ARRAYNAME->data[0]);\
		ARRAYNAME=malloc(sizeof(*ARRAYNAME)+_cap);\
		if(!ARRAYNAME)\
		{\
			CRASH("Alloc error");\
			return ERR_RET;\
		}\
		ARRAYNAME->count=0;\
		ARRAYNAME->cap=_cap;\
		ARRAYNAME->esize=sizeof(ARRAYNAME->data[0]);\
		ARRAYNAME->destructor=DESTRUCTOR;\
		if(ZEROMEM)\
			memset(ARRAYNAME->data, 0, _cap);\
	}while(0)

//flat copy, set the destructor manually if required
#define ARRAY_COPY(ARRAYNAME, SRC, LEN, PAD, ERR_RET)\
	do\
	{\
		ptrdiff_t _count=(LEN), _pad=(PAD), _cap=(_count+_pad)*sizeof(ARRAYNAME->data[0]);\
		ARRAYNAME=malloc(sizeof(*ARRAYNAME)+_cap);\
		if(!ARRAYNAME)\
		{\
			CRASH("Alloc error");\
			return ERR_RET;\
		}\
		ARRAYNAME->count=_count;\
		ARRAYNAME->cap=_cap;\
		ARRAYNAME->esize=sizeof(ARRAYNAME->data[0]);\
		ARRAYNAME->destructor=0;\
		memcpy(ARRAYNAME->data, SRC, _count*sizeof(ARRAYNAME->data[0]));\
		memset(ARRAYNAME->data+_count, 0, _pad*sizeof(ARRAYNAME->data[0]));\
	}while(0)

#define ARRAY_FREE(ARRAYNAME)\
	do\
	{\
		if((ARRAYNAME)->destructor)\
		{\
			ptrdiff_t k;\
			\
			for(k=0;k<(ARRAYNAME)->count;++k)\
				(ARRAYNAME)->destructor(ARRAYNAME->data+k);\
		}\
		free(ARRAYNAME);\
		ARRAYNAME=0;\
	}while(0)

//does not affect array->count, only affects array->cap
#define ARRAY_REALLOC(ARRAYNAME, PAD, ERR_RET, ZEROMEM)\
	do\
	{\
		ptrdiff_t _newcap=(PAD)*sizeof(ARRAYNAME->data[0]);\
		if(_newcap>ARRAYNAME->cap)\
		{\
			void *_p=realloc(ARRAYNAME, sizeof(*ARRAYNAME)+_newcap);\
			if(!_p)\
			{\
				CRASH("Alloc error");\
				return ERR_RET;\
			}\
			ARRAYNAME=_p;\
			if(ZEROMEM)\
				memset(ARRAYNAME->data+ARRAYNAME->count, 0, _newcap-ARRAYNAME->cap);\
			ARRAYNAME->cap=_newcap;\
		}\
	}while(0)

#define ARRAY_APPEND1(ARRAYNAME, ERR_RET)\
	ARRAY_REALLOC(ARRAYNAME, (ARRAYNAME->count+64)&~31, ERR_RET, 1)

#define ARRAY_APPEND(ARRAYNAME, DELTA, ERR_RET, ZEROMEM)\
	ARRAY_REALLOC(ARRAYNAME, ARRAYNAME->count+(DELTA), ERR_RET, ZEROMEM)

#define ARRAY_BACK(ARRAYNAME, IDX)\
	ARRAYNAME->data[ARRAYNAME->count-(IDX)]

#define ARRAY_PUSHBACKN(ARRAYNAME, DATA, COUNT)\
	do\
	{\
		ptrdiff_t _delta=COUNT;\
		memcpy(&ARRAY_BACK(ARRAYNAME, 0), DATA, _delta*sizeof(ARRAYNAME->data[0]));\
		ARRAYNAME->count+=_delta;\
	}while(0)

#define ARRAY_POPBACKN(ARRAYNAME, COUNT)\
	do\
	{\
		ptrdiff_t _delta=COUNT;\
		ARRAYNAME->count-=_delta;\
		memset(&ARRAY_BACK(ARRAYNAME), 0, _delta*sizeof(ARRAYNAME->data[0]));\
	}while(0);

#define ARRAY_PUSHBACK(ARRAYNAME, DATA)\
	(ARRAYNAME->data[ARRAYNAME->count++]=DATA)

#define ARRAY_POPBACK(ARRAYNAME)\
	(ARRAYNAME->data[--ARRAYNAME->count]=0)
#endif

ARRAY_DECL(String, char);

static String filter_path(const char *path, int len, int trailingslash, int pad)//replaces back slashes with slashes, control trailing slash
{
	String path2=0;
	ptrdiff_t k;

	if(len<0)
		return 0;
	if(!len)
		len=(int)strlen(path);
	ARRAY_COPY(path2, path, len, pad+1, 0);
	for(k=0;k<path2->count;++k)//replace back slashes
	{
		if(path2->data[k]=='\\')
			path2->data[k]='/';
	}
	k=ARRAY_BACK(path2, 1);
	if(trailingslash)
	{
		if(k!='/')//ensure trailing slash
			ARRAY_PUSHBACK(path2, '/');
	}
	else if(k=='/')
		ARRAY_POPBACK(path2);
	return path2;
}
static int acme_stricmp(const char *a, const char *b)//case insensitive strcmp
{
	if(!a||!b)
		return !a&&!b;
	while(*a&&tolower(*a)==tolower(*b))
		++a, ++b;
	return (*a>*b)-(*a<*b);
}
static const char* get_extension(const char *filename, ptrdiff_t len)//excludes the dot
{
	const char *ptr=filename+len;

	while(--ptr>=filename)
	{
		if(*ptr=='.')
			return ptr+1;
	}
	return 0;
#if 0
	const char *dot=strrchr(filename, '.');//https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
	if(!dot||dot==filename)
		return "";
	return dot+1;
#endif
}
static void get_filetitle(const char *fn, int len, int *idx_start, int *idx_end)//pass -1 for len if unknown
{
	int kpoint, kslash;
	if(len<0)
		len=(int)strlen(fn);
	for(kpoint=(int)len-1;kpoint>=0&&fn[kpoint]!='.';--kpoint);
	if(kpoint<0)
		kpoint=len;
	for(kslash=kpoint-1;kslash>=0&&fn[kslash]!='/'&&fn[kslash]!='\\';--kslash);
	++kslash;
	if(idx_start)
		*idx_start=kslash;
	if(idx_end)
		*idx_end=kpoint;
}
#if 0
ARRAY_DECL(Filenames, String);
static Filenames get_filenames(const char *path, const char **extensions, int extCount, int fullyqualified)
{
#if defined _MSC_VER || defined _WIN32
	String searchpath=0;
	Filenames filenames=0;
	WIN32_FIND_DATAA data={0};
	void *hSearch=0;
	int success=0;
	const char *extension=0;
	ptrdiff_t len=0;
	
	searchpath=filter_path(path, 0, 1, 10);
	ARRAY_PUSHBACK(searchpath, '*');

	hSearch=FindFirstFileA(searchpath->data, &data);//skip .
	if(hSearch==INVALID_HANDLE_VALUE)
	{
		free(searchpath);
		return 0;
	}
	success=FindNextFileA(hSearch, &data);//skip ..
	ARRAY_POPBACK(searchpath);
	ARRAY_ALLOC(filenames, 20, free, 0, 1);
	for(;;)
	{
		success=FindNextFileA(hSearch, &data);
		if(!success)
			break;
		len=strlen(data.cFileName);
		extension=get_extension(data.cFileName, len);
		if(!(data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY))
		{
			for(int k=0;k<extCount;++k)
			{
				if(!acme_stricmp(extension, extensions[k]))
				{
					String fn;

					ARRAY_ALLOC(fn, (fullyqualified?searchpath->count:0)+len+5, 0, 0, 1);
					if(fullyqualified)
						ARRAY_PUSHBACKN(fn, searchpath->data, searchpath->count);
					ARRAY_PUSHBACKN(fn, data.cFileName, len);

					ARRAY_APPEND1(filenames, 0);
					ARRAY_PUSHBACK(filenames, fn);
					break;
				}
			}
		}
	}
	success=FindClose(hSearch);
	free(searchpath);
	return filenames;
#elif defined __linux__
	String searchpath=0;
	Filenames filenames=0;
	struct dirent *dir=0;
	DIR *d=0;
	
	searchpath=filter_path(path, 0, 1, 10);
	ARRAY_PUSHBACK(searchpath, '*');

	d=opendir(searchpath->data);
	if(!d)
		return 0;
	ARRAY_ALLOC(filenames, 20, free, 0, 1);
	while((dir=readdir(d)))
	{
		if(dir->d_type==DT_REG)//regular file
		{
			const char *name=dir->d_name;
			ptrdiff_t len=strlen(name);
			const char *extension=get_extension(name, len);
			for(int k=0;k<extCount;++k)
			{
				if(!acme_stricmp(extension, extensions[k]))
				{
					String fn;
					
					ARRAY_ALLOC(fn, (fullyqualified?len0:0)+len+5, 0, 0, 1);
					if(fullyqualified)
						ARRAY_PUSHBACKN(fn, searchpath->data, searchpath->count);
					ARRAY_PUSHBACKN(fn, name, len);

					ARRAY_APPEND1(filenames, 1);
					ARRAY_PUSHBACK(filenames, fn);
					break;
				}
			}
		}
	}
	closedir(d);
	free(searchpath);
	qsort(filenames->data, filenames->count, sizeof(filenames->data[0]), (int(*)(const void*, const void*))strcmp);
	return filenames;
#endif
}
#endif
static int acme_getline(char *buf, int len, FILE *f)
{
	int k;

	memset(buf, '\n', len);
	fgets(buf, len, f);
	--len;
	for(k=0;k<len&&buf[k]!='\n';++k);
	buf[k]=0;
	return k;
}
static int acme_strnimatch(const char *s1, ptrdiff_t len1, const char *s2, ptrdiff_t len2)//return 1: match		FIXME return ASCII order & index of first difference
{
	if(len1!=len2)
		return 0;
	const char *end1=s1+len1;
	while(s1<end1&&tolower(*s1)==tolower(*s2))++s1, ++s2;//check then increment  (blind increment misses the last character at s1==end1)
	return s1==end1;
}
#if 0
ARRAY_DECL(Buffer, uint8_t);
static Buffer file_load(const char *fn)
{
	Buffer buf=0;
	struct stat info={0};
	FILE *f=0;

	stat(fn, &info);
	f=fopen(fn, "rb");
	if(!f)
		return 0;
	ARRAY_ALLOC(buf, info.st_size+16, 0, 0, 1);
	buf->count=fread(buf, 1, info.st_size, f);
	fclose(f);
	return buf;
}
#endif

typedef struct _Image
{
	int iw, ih;
	uint8_t data[];
} *Image;
static Image ppm_load(const char *fn)
{
	int iw=0, ih=0, c=0;
	FILE *fsrc=0;
	ptrdiff_t size=0, nread=0;
	Image image=0;
	
	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	c|=fgetc(fsrc)<<8;
	if(c!=('P'|'6'<<8))
	{
		CRASH("Unsupported PPM file \"%s\"", fn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		CRASH("Unsupported PPM file");
		return 0;
	}
	c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	iw=0;
	while((uint32_t)(c-'0')<10)
	{
		iw=10*iw+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	ih=0;
	while((uint32_t)(c-'0')<10)
	{
		ih=10*ih+c-'0';
		c=fgetc(fsrc);
	}
	while(c<=' ')
		c=fgetc(fsrc);
	while(c=='#')
	{
		c=fgetc(fsrc);
		while(c!='\n')
			c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	c=c<<8|fgetc(fsrc);
	if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
	{
		CRASH("Unsupported PPM file");
		return 0;
	}
	size=(ptrdiff_t)3*iw*ih;
	image=(Image)malloc(sizeof(*image)+size+16);
	if(!image)
	{
		CRASH("Alloc error");
		return 0;
	}
	image->iw=iw;
	image->ih=ih;
	nread=fread(image->data, 1, size, fsrc);
	fclose(fsrc);
	if(nread!=size)
		printf("Truncated file  expected %td  got %td", size, nread);
	return image;
}
static int verify_files(const char *fn1, const char *fn2, double *ret_rmse)//0  identical    1  different
{
	Image im1=0, im2=0;
	int iw=0, ih=0;
	ptrdiff_t size=0, k, diffidx;
	const uint8_t *ptr1=0, *ptr2=0;
	int ret=0;

	im1=ppm_load(fn1);
	im2=ppm_load(fn2);
	if(im1->iw!=im2->iw||im1->ih!=im2->ih)
		printf("WH %5d*%5d != %5d*%5d"
			, im1->iw, im1->ih
			, im2->iw, im2->ih
		);
	iw=im1->iw<im2->iw?im1->iw:im2->iw;
	ih=im1->ih<im2->ih?im1->ih:im2->ih;
	size=(ptrdiff_t)3*iw*ih;
	ptr1=im1->data;
	ptr2=im2->data;
	{
		__m128i acc=_mm_setzero_si128();
		diffidx=-1;
		for(k=0;k<size-63;k+=64)
		{
			__m128i a0=_mm_loadu_si128((__m128i*)ptr1+0);
			__m128i a1=_mm_loadu_si128((__m128i*)ptr1+1);
			__m128i a2=_mm_loadu_si128((__m128i*)ptr1+2);
			__m128i a3=_mm_loadu_si128((__m128i*)ptr1+3);
			a0=_mm_xor_si128(a0, _mm_loadu_si128((__m128i*)ptr2+0));
			a1=_mm_xor_si128(a1, _mm_loadu_si128((__m128i*)ptr2+1));
			a2=_mm_xor_si128(a2, _mm_loadu_si128((__m128i*)ptr2+2));
			a3=_mm_xor_si128(a3, _mm_loadu_si128((__m128i*)ptr2+3));
			acc=_mm_or_si128(acc, a0);
			acc=_mm_or_si128(acc, a1);
			acc=_mm_or_si128(acc, a2);
			acc=_mm_or_si128(acc, a3);
			{
				int mask=_mm_movemask_epi8(acc);
				if(mask)
				{
					int k2;
					for(k2=0;k2<64;++k2)
					{
						if(ptr1[k2]!=ptr2[k2])
						{
							diffidx=k+k2;
							break;
						}
					}
					break;
				}
			}
			ptr1+=64;
			ptr2+=64;
		}
	}
	if(diffidx==-1)
	{
		memset(ret_rmse, 0, sizeof(double[3]));
		ret=0;
	}
	else
	{
		static const double jpeg_matrix[]=
		{
			+0.299		,//0
			+0.587		,//1
			+0.114		,//2

			-0.168736	,//3
			-0.331264	,//4
			+0.5		,//5

			+0.5		,//6
			-0.418688	,//7
			-0.081312	,//8
		};
		double e[3]={0};
		ptrdiff_t res=(ptrdiff_t)iw*ih;
		int kx, ky;

		ptr1=im1->data;
		ptr2=im2->data;
		for(ky=0;ky<ih;++ky)
		{
			for(kx=0;kx<iw;++kx)
			{
				int dr=*ptr1++-*ptr2++;
				int dg=*ptr1++-*ptr2++;
				int db=*ptr1++-*ptr2++;
				double yuv[]=
				{
					jpeg_matrix[0]*dr+jpeg_matrix[1]*dg+jpeg_matrix[2]*db,
					jpeg_matrix[3]*dr+jpeg_matrix[4]*dg+jpeg_matrix[5]*db,
					jpeg_matrix[6]*dr+jpeg_matrix[7]*dg+jpeg_matrix[8]*db,
				};
				e[0]+=yuv[0]*yuv[0];//int8 * int8 = int16	static analysis is retarded
				e[1]+=yuv[1]*yuv[1];
				e[2]+=yuv[2]*yuv[2];
			}
		}
		ret_rmse[0]=sqrt((double)e[0]/res);
		ret_rmse[1]=sqrt((double)e[1]/res);
		ret_rmse[2]=sqrt((double)e[2]/res);
		ret=1;
	}
	free(im1);
	free(im2);
	return ret;
}

typedef union _DateTime
{
	struct
	{
		uint8_t cs, second, minute, hour, day, month;
		uint16_t year;
	};
	uint64_t timestamp;
} DateTime;
typedef struct _UInfo//information about an image
{
	int64_t usize;
	String filename;
} UInfo;
ARRAY_DECL(Dataset, UInfo);
static void free_uinfo(void *p)
{
	UInfo *info=(UInfo*)p;

	free(info->filename);
}
typedef struct _ImCodecResult//information about an image compressed with a codec
{
	int64_t csize;
	double etime, dtime;
	double emem, dmem;
	double rmse[3];
} ImCodecResult;
ARRAY_DECL(TableRow, ImCodecResult);
typedef struct _TestInfo//each test has a codecname, date, and many compression results
{
	String codecname;
	DateTime datetime;
	ImCodecResult total;
	TableRow cells;
} TestInfo;
ARRAY_DECL(Tests, TestInfo);
static void free_testinfo(void *p)
{
	TestInfo *info=(TestInfo*)p;
	free(info->codecname);
	free(info->cells);
}
typedef struct _Range
{
	int issrc, start, end;
} Range;
ARRAY_DECL(Ranges, Range);
typedef struct _CommandFormat
{
	String format;
	Ranges bounds;
} CommandFormat;
ARRAY_DECL(IntArray, int);
static String txt_load(const char *fn)
{
	String txt=0;
	struct stat info={0};
	FILE *f=0;

	stat(fn, &info);
	f=fopen(fn, "r");
	if(!f)
		return 0;
	ARRAY_ALLOC(txt, info.st_size+16, 0, 0, 1);
	txt->count=fread(txt->data, 1, info.st_size, f);
	fclose(f);
	return txt;
}
static int txt_save(const char *fn, const char *txt, int len)
{
	FILE *f=fopen(fn, "w");
	if(!f)
	{
		printf("Cannot open \"%s\" for writing\n", fn);
		return 1;
	}
	fwrite(txt, 1, len, f);
	fclose(f);
	return 0;
}
static int getlineno(const char *start, const char *ptr)
{
	int lineno=1;
	while(ptr<start)
		lineno+=*ptr++=='\n';
	return lineno;
}
static void skipspace(const char **ptr)
{
	const char *ptr2=*ptr;
	while(*ptr2&&isspace(*ptr2))++ptr2;
	*ptr=ptr2;
}
static void skip2space(const char **ptr)
{
	const char *ptr2=*ptr;
	while(*ptr2&&!isspace(*ptr2))++ptr2;
	*ptr=ptr2;
}
static int peeklabel(const char *ptr, const char *label)
{
	const char *lptr=label;
	while(*lptr&&*ptr++==*lptr++);
	return !*lptr;
}
static void skiplabel(const char *start, const char **ptr, const char *label)
{
	if(!*ptr)
	{
		CRASH("");
		return;
	}
	const char *ptr2=*ptr, *lptr=label;
	while(*lptr&&*ptr2++==*lptr++);
	if(*lptr)
	{
		int lineno=getlineno(start, *ptr);
		int len=(int)strlen(*ptr);
		if(len>(int)(lptr-label))
			len=(int)(lptr-label);
		printf("Parser line %d expected \"%s\" got \"%.*s\"", lineno, label, len, *ptr);
		CRASH("");
	}
	*ptr=ptr2;
}
static String parse_str(const char *start, const char **ptr, char delim, int allow_eof)//delim must be one of the whitespace characters
{
	String str=0;
	const char *ptr2=*ptr;
	for(;*ptr2&&*ptr2!=delim;++ptr2);
	if(allow_eof&&!*ptr2)
		goto skipcheck;
	if(*ptr2!=delim)
	{
		int lineno=getlineno(start, *ptr);
		int len=(int)strlen(*ptr);
		if(len>10)
			len=10;
		printf("Parser line %d expected \"%c\"=0x%02X got \"%.*s\"", lineno, delim, delim, len, *ptr);
		CRASH("");
	}
skipcheck:
	{
		const char *ptr3=ptr2;
		while(ptr3>*ptr&&isspace(*(ptr3-1)))--ptr3;
		ARRAY_COPY(str, *ptr, ptr3-*ptr, 1, 0);
	}
	*ptr=ptr2;
	return str;
}
static int64_t parse_uint(const char **ptr)
{
	const char *ptr2=*ptr;
	int64_t val=0;
	skipspace(&ptr2);
	while((uint32_t)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	*ptr=ptr2;
	return val;
}
static double parse_float(const char **ptr)
{
	const char *ptr2=*ptr;
	double val=0;
	skipspace(&ptr2);
	while((uint32_t)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	if(*ptr2=='.')
	{
		++ptr2;
		double p=0.1;
		while((uint32_t)(*ptr2-'0')<10)
		{
			val+=(*ptr2++-'0')*p;
			p*=0.1;
		}
	}
	*ptr=ptr2;
	return val;
}
static void parse_cell(const char **ptr, ImCodecResult *cell)
{
	//{int csize; float etime, dtime, emem, dmem, rmse_y, rmse_u, rmse_v;}		at .4 digits of precision
	cell->csize=parse_uint(ptr);
	cell->etime=parse_float(ptr);
	cell->dtime=parse_float(ptr);
	cell->emem=parse_float(ptr);
	cell->dmem=parse_float(ptr);
	cell->rmse[0]=parse_float(ptr);
	cell->rmse[1]=parse_float(ptr);
	cell->rmse[2]=parse_float(ptr);
}
static int strimatch(const char *text, const char *label)//return 1: match
{
	while(*label&&tolower(*text)==tolower(*label))++text, ++label;
	return !*label;
}
static Ranges parse_cmdformatext(const char *text, int len, const char *srcext, int requiredst, int *srcdstbounds)
{
	Ranges bounds=0;
	const char *ptr=text, *end=text+len-1;
	memset(srcdstbounds, -1, sizeof(int[4]));
	ARRAY_ALLOC(bounds, 2, 0, 0, 1);
	while(ptr<end)
	{
		short c=*(short*)ptr;
		if(c==('*'|'.'<<8))
		{
			Range *range=0;
			int idx;
			
			ARRAY_APPEND1(bounds, 0);
			++bounds->count;
			range=&ARRAY_BACK(bounds, 1);

			idx=(int)(ptr-text);
			range->issrc=strimatch(text+idx+2, srcext);//add 2 to skip "*."
			range->start=idx;
			skip2space(&ptr);
			range->end=(int)(ptr-text);
			int *ridx=range->issrc?srcdstbounds:srcdstbounds+2;
			if(*ridx==-1)
			{
				ridx[0]=range->start;
				ridx[1]=range->end;
			}
			else if(!acme_strnimatch(
				text+range->start, (ptrdiff_t)range->end-range->start,
				text+ridx[0], (ptrdiff_t)ridx[1]-ridx[0]
			))
			{
				printf(
					"\n"
					"Invalid template:  %s\n"
					"Expected\n"
					"  [%d] \"%.*s\"  to match\n"
					"  [%d] \"%.*s\"  (case-insensitive, no spaces)\n",
					text,
					ridx[0], ridx[1]-ridx[0], text+ridx[0],
					range->start, range->end-range->start, text+range->start
				);
				CRASH("");
			}
			continue;
		}
		++ptr;
	}
	if(!bounds)
	{
		printf(
			"\n"
			"Invalid template:  %s\n"
			"Missing source%s \"*.extension%s\"\n",
			text,
			requiredst?" and codec":"",
			requiredst?"s":""
		);
		CRASH("");
	}
	if(requiredst&&srcdstbounds[2]==-1)
	{
		printf(
			"\n"
			"Invalid template:  %s\n"
			"Missing codec \"*.extension\"\n",
			text
		);
		CRASH("");
	}
	return bounds;
}
static void substitute_cmdplaceholders(char **pptr, const char *end, const CommandFormat *cmd, const char *t1fn, const char *t2fn)
{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wrestrict"
#endif
	char *ptr=*pptr;
	int lastidx=0;
	for(int k=0;k<(int)cmd->bounds->count;++k)
	{
		const Range *range=cmd->bounds->data+k;
		ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" ",
			range->start-lastidx, (char*)cmd->format->data+lastidx,
			range->issrc?t1fn:t2fn
		);
		lastidx=range->end;
	}
	ptr+=snprintf(ptr, end-ptr, "%s",
		(char*)cmd->format->data+lastidx
	)+1;
	*pptr=ptr;
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
}
static Dataset get_uinfo(const char *path, const char *ext)//extension without '.'
{
	String searchpath=0;
	void *hSearch=0;
	WIN32_FIND_DATAA data={0};
	Dataset dataset=0;

	//prepare searchpath
	searchpath=filter_path(path, 0, 1, 10);
	ARRAY_PUSHBACK(searchpath, '*');

	hSearch=FindFirstFileA((char*)searchpath->data, &data);//skip .
	if(hSearch==INVALID_HANDLE_VALUE)
		return 0;
	FindNextFileA(hSearch, &data);//skip ..
	ARRAY_POPBACK(searchpath);
	ARRAY_ALLOC(dataset, 20, free_uinfo, 0, 1);
	while(FindNextFileA(hSearch, &data))
	{
		ptrdiff_t len=strlen(data.cFileName);
		const char *curr_ext=get_extension(data.cFileName, len);
		if(!(data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)&&!acme_stricmp(curr_ext, ext))
		{
			UInfo *info=0;
			String fn=0;

			ARRAY_APPEND1(dataset, 0);
			++dataset->count;
			info=&ARRAY_BACK(dataset, 1);
			info->usize=(int64_t)data.nFileSizeHigh<<32|data.nFileSizeLow;
			ARRAY_COPY(fn, searchpath->data, searchpath->count, len+1, 0);
			memcpy(fn->data+fn->count, data.cFileName, len);
			fn->count+=len;
			info->filename=fn;
		}
	}
	int success=FindClose(hSearch);
	if(!success)
	{
		CRASH("FindClose %d", GetLastError());
		return 0;
	}
	free(searchpath);
	return dataset;
}
static void write_cell(FILE *fdst, ImCodecResult *cell)//tabs for excel
{
	fprintf(fdst, "\t%10lld\t%12.6lf\t%12.6lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf"
		, cell->csize
		, cell->etime
		, cell->dtime
		, cell->emem
		, cell->dmem
		, cell->rmse[0]
		, cell->rmse[1]
		, cell->rmse[2]
	);
}
static int qualify_command(const char *cmd, int len, char *ret)
{
	const char *ptr=0;
	char fn1[MAX_PATH+1]={0}, fn2[MAX_PATH+1]={0};
	int len1=0, len2=0;
	struct stat info={0};
	int remlen=0;

	ptr=cmd;
	skip2space(&ptr);//assume no spaces

	len1=(int)(ptr-cmd);
	memcpy(fn1, cmd, len1);
	if(!strimatch(fn1+len1-4, ".exe"))
	{
		memcpy(fn1+len1, ".exe", 4);
		len1+=4;
	}
	len2=SearchPathA(0, fn1, 0, MAX_PATH, fn2+1, 0);
	if(!len2)
	{
		CRASH("GetFullPathNameA %d", GetLastError());
		return -1;
	}
	stat(fn2+1, &info);
#ifdef _MSC_VER
	if(!(info.st_mode&S_IFREG))
#else
	if(!S_ISREG(info.st_mode))
#endif
	{
		printf("Program is inaccessible:\n");
		printf("  \"%s\"\n", fn1);
		printf("  \"%s\"\n", fn2+1);
		printf("  size %lld bytes", (int64_t)info.st_size);
		CRASH("");
		return -1;
	}
	len2+=2;
	fn2[0]='\"';
	fn2[len2-1]='\"';
	remlen=len-(int)(ptr-cmd);
	if(ret)
	{
		memcpy(ret, fn2, len2);
		memcpy(ret+len2, ptr, remlen);
		ret[len2+remlen]=0;
	}
	return len2+remlen;
}
static void exec_process2(char *cmd, const char *currdir, int loud, double *elapsed, int64_t *maxmem)
{
	int success;
	STARTUPINFOA si={0};
	PROCESS_INFORMATION pi={0};
	si.cb=sizeof(si);
	if(!loud)
	{
		si.dwFlags=STARTF_USESTDHANDLES;
		si.hStdOutput=CreateFileA("NUL", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdOutput==INVALID_HANDLE_VALUE)
		{
			CRASH("CreateFileA %d", GetLastError());
			return;
		}
		si.hStdError=CreateFileA("NUL", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdError==INVALID_HANDLE_VALUE)
		{
			CRASH("CreateFileA %d", GetLastError());
			return;
		}
		si.hStdInput=CreateFileA("NUL", GENERIC_READ, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdInput==INVALID_HANDLE_VALUE)
		{
			CRASH("CreateFileA %d", GetLastError());
			return;
		}
	}
	success=CreateProcessA(0, cmd, 0, 0, 0, CREATE_SUSPENDED, 0, currdir, &si, &pi);
	if(!success)
	{
		CRASH("CreateProcessA %d", GetLastError());
		return;
	}
	ptrdiff_t memusage=0;
	double t=time_sec();
	int suspendcount=ResumeThread(pi.hThread);
	if(suspendcount==(DWORD)-1)
	{
		CRASH("CreateProcessA %d", GetLastError());
		return;
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
	if(elapsed)*elapsed=time_sec()-t;
	FILETIME tstart={0}, tfinish={0}, tkernel={0}, tuser={0};
	//timer7
	GetProcessTimes(pi.hProcess, &tstart, &tfinish, &tkernel, &tuser);
	PROCESS_MEMORY_COUNTERS counters={sizeof(PROCESS_MEMORY_COUNTERS)};
	GetProcessMemoryInfo(pi.hProcess, &counters, sizeof(counters));
	memusage=counters.PeakWorkingSetSize;
	if(maxmem)*maxmem=memusage;
	success=CloseHandle(pi.hThread);
	if(!success)
	{
		CRASH("CloseHandle %d", GetLastError());
		return;
	}
	success=CloseHandle(pi.hProcess);
	if(!success)
	{
		CRASH("CloseHandle %d", GetLastError());
		return;
	}
	if(!loud)
	{
		success=CloseHandle(si.hStdOutput);
		if(!success)
		{
			CRASH("CloseHandle %d", GetLastError());
			return;
		}
		success=CloseHandle(si.hStdError);
		if(!success)
		{
			CRASH("CloseHandle %d", GetLastError());
			return;
		}
		success=CloseHandle(si.hStdInput);
		if(!success)
		{
			CRASH("CloseHandle %d", GetLastError());
			return;
		}
	}
}
static void print_cmpresult(double *rmse)
{
	printf("PSNR %10.4lf"
		, 20*log10(255*8/(6*rmse[0]+rmse[1]+rmse[2]))
	);
}
static void print_summary(IntArray besttestidxs, Tests testinfo, ptrdiff_t usize, int special)
{
	double ebest=0, dbest=0;
	for(int k2=0;k2<(int)besttestidxs->count;++k2)
	{
		int idx=besttestidxs->data[k2];
		TestInfo *test=testinfo->data+idx;
		int epareto=0, dpareto=0;

		if(!k2||ebest>test->total.etime)
			ebest=test->total.etime, epareto=1;
		if(!k2||dbest>test->total.dtime)
			dbest=test->total.dtime, dpareto=1;
		if(k2&&k2==special)
			printf("\n");
		printf("%10lld B  %12.6lf%c %12.6lf%c sec  %10.4lf %10.4lf MB/s %8.2lf %8.2lf MB  "
			, test->total.csize
			, test->total.etime, epareto?'*':' '
			, test->total.dtime, dpareto?'*':' '
			, usize/(test->total.etime*1024*1024)
			, usize/(test->total.dtime*1024*1024)
			, test->total.emem
			, test->total.dmem
		);
		print_cmpresult(test->total.rmse);
		printf("  %04d%02d%02d_%02d%02d%02d"
			, test->datetime.year
			, test->datetime.month
			, test->datetime.day
			, test->datetime.hour
			, test->datetime.minute
			, test->datetime.second
		);
		if((uint32_t)special<(uint32_t)besttestidxs->count)
		{
			if(k2==special)
				printf("  %X <- *", (k2+1)&15);
			else if(k2<special)
				printf("  %X     ", (k2+1)&15);
			else
				printf("  %X <- %X", (k2+1)&15, k2&15);
		}
		else
			printf("  %X", (k2+1)&15);
		printf(" %s\n", (char*)test->codecname->data);
		if(k2==special)
			printf("\n");
	}
}
static void ascii_deletefile(const char *fn)
{
	struct stat info={0};
	int error=stat(fn, &info);
	if(error)
	{
		printf("Cannot stat delete target \"%s\"\n", fn);
		CRASH("");
	}
	if(!error)
	{
		int success=DeleteFileA(fn);
		if(!success)
			CRASH("DeleteFileA", GetLastError());
	}
}
#define COLORPRINTF_BK_DEFAULT 0x0C0C0C
#define COLORPRINTF_TXT_DEFAULT 0xF2F2F2
static int colorprintf(int textcolor, int bkcolor, const char *format, ...)//0x00BBGGRR
{
	int printed=0, msg=0;
	va_list args;

	printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "\33[48;2;%d;%d;%d;38;2;%d;%d;%dm",
		bkcolor&255, bkcolor>>8&255, bkcolor>>16&255,
		textcolor&255, textcolor>>8&255, textcolor>>16&255
	);
	va_start(args, format);
	msg=vsnprintf(g_buf+printed, sizeof(g_buf)-1-printed, format, args);
	printed+=msg;
	va_end(args);
	printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "\33[0m");

	printf("%s", g_buf);

	return msg;
}
static int acme_strftime(char *buf, size_t len, const char *format)
{
#ifdef _MSC_VER
	time_t tstamp;
	struct tm tformat;

	tstamp=time(0);
	localtime_s(&tformat, &tstamp);
	return (int)strftime(buf, len, format, &tformat);
#else
	time_t tstamp=time(0);
	struct tm *tformat=localtime(&tstamp);
	return (int)strftime(buf, len, format, tformat);
#endif
}
static int print_currtimestamp(const char *format)
{
	acme_strftime(g_buf, sizeof(g_buf)-1, format);
	return printf("%s", g_buf);
}

#if 0
static IntArray get_paretofront(IntArray besttestidxs, Tests testinfo, int dec)
{
	IntArray sorted=0, stack=0;
	int lowest=0;
	TestInfo *info=testinfo->data, *p0=0;

	//Graham's scan		https://en.wikipedia.org/wiki/Graham_scan
	ARRAY_ALLOC(sorted, besttestidxs->count, 0, 0, 1);
	ARRAY_ALLOC(stack, besttestidxs->count, 0, 0, 1);

	//get lowest (and leftmost) point P0
	for(int ki=1;ki<(int)besttestidxs->count;++ki)
	{
		int idx=besttestidxs->data[ki];
		double elowest	=(&info[lowest	].total.etime)[dec];
		double e2	=(&info[idx	].total.etime)[dec];
		if(elowest>e2||(elowest==e2&&info[lowest].total.csize<info[idx].total.csize))
			lowest=idx;
	}
	p0=info+lowest;

	//sort by decreasing polar angle with P0
	for(int ki=0;ki<(int)besttestidxs->count;++ki)
	{
		int idx=0, ki2=0;
		TestInfo *p1=0;
		double dx=0, dy=0, angle1=0, dist1=0;
		int insert=1;

		idx=besttestidxs->data[ki];
		p1=info+idx;
		if(p1->total.csize>p0->total.csize)//reject anything to the right of p0
			continue;
		dx=(double)(p1->total.csize-p0->total.csize)/(1024.*1024);
		dy=(&p1->total.etime)[dec]-(&p0->total.etime)[dec];
		angle1=idx==lowest?0:atan2(dy, dx)*(180/M_PI);
		dist1=dx*dx+dy*dy;

		for(ki2=0;ki2<sorted->count;++ki2)
		{
			int idx2=0;
			TestInfo *p2=0;
			double angle2=0, dist2=0;

			idx2=sorted->data[ki2];
			p2=info+idx2;
			dx=(double)(p2->total.csize-p0->total.csize)/(1024.*1024);
			dy=(&p2->total.etime)[dec]-(&p0->total.etime)[dec];
			angle2=idx2==lowest?0:atan2(dy, dx)*(180/M_PI);
			dist2=dx*dx+dy*dy;

			/*
			gralic		close to 90 deg
			lea
			flic
			c32
			qlic2
			halic072fast
			fox		close to 180 deg
			copy	= p0	undefined
			*/
			if(angle2<angle1)//smaller angle first
				break;
			
			if(angle2==angle1)//select furthest on tie
			{
				insert=0;
				if(dist2<dist1)//replace idx2 with idx1 if dist1 further
					sorted->data[ki2]=idx;
				break;
			}
		}
		if(insert)
		{
			if(sorted->count>ki2)
				memmove(sorted->data+ki2+1, sorted->data+ki2, (sorted->count-ki2)*sizeof(sorted->data[0]));
			++sorted->count;
			sorted->data[ki2]=idx;
		}
	}

	for(int ki=0;ki<(int)sorted->count;++ki)
	{
		int curr=sorted->data[ki];
		TestInfo *p0=info+curr;
		double x0=(double)p0->total.csize;
		double y0=(&p0->total.etime)[dec];
		while(stack->count>1)
		{
			int top1=stack->data[stack->count-1];
			int top2=stack->data[stack->count-2];
			
			TestInfo *p1=info+top1;
			double x1=(double)p1->total.csize;
			double y1=(&p1->total.etime)[dec];
			
			TestInfo *p2=info+top2;
			double x2=(double)p2->total.csize;
			double y2=(&p2->total.etime)[dec];
			/*
			ccw(top2, top1, curr) := det <= 0
			x2	y2	1	top2
			x1	y1	1	top1
			x0	y0	1	curr

			x0*y2+y0*x1+y1*x2-(x0*y0+y0*x2+x1*y2)
			*/
			int ccw=x2*y1+x1*y0+x0*y2 <= x2*y0+x1*y2+x0*y1;
			if(!ccw)
				break;
			--stack->count;
		}
		stack->data[stack->count++]=curr;
	}
	//free(idxs);
	free(sorted);
	return stack;
}
#endif
#define SETPIXEL(X, Y, V) do{if((uint32_t)(X)<(uint32_t)cw&&(uint32_t)(Y)<(uint32_t)ch)buf[cw*(Y)+(X)]=(V);}while(0)
static void draw_curve(IntArray idxs, Tests testinfo, char *buf, int cw, int ch, double px, double xmin, double py, double ymin, int xmargin, int ymargin, int dec, char c)
{
	for(int k=1;k<(int)idxs->count;++k)
	{
		int idx1=idxs->data[k-1], idx2=idxs->data[k];
		TestInfo *p1=testinfo->data+idx1;
		TestInfo *p2=testinfo->data+idx2;
		int32_t x1=(int32_t)(((double)p1->total.csize-xmin)/px)+xmargin;
		int32_t x2=(int32_t)(((double)p2->total.csize-xmin)/px)+xmargin;
		int32_t y1=(int32_t)(((double)(&p1->total.etime)[dec]-ymin)/py)+ymargin;
		int32_t y2=(int32_t)(((double)(&p2->total.etime)[dec]-ymin)/py)+ymargin;

		//printf("line XY1 %3d %3d,  XY2 %3d %3d\n", x1, y1, x2, y2);

		{//https://every-algorithm.github.io/2024/10/11/bresenhams_line_algorithm.html
			int32_t dx=x2-x1, dy=y2-y1;
			int sx=(dx>0)-(dx<0), sy=(dy>0)-(dy<0);
			int x=x1, y=y1, error;

			SETPIXEL(x, y, c);
			dx=abs(dx);
			dy=abs(dy);
			if(dx>dy)
			{
				error=dx-dy;
				while(x!=x2)
				{
					int e2=2*error;
					if(e2>-dy)
						error-=dy, x+=sx;
					if(e2<dx)
						error+=dx, y+=sy;
					SETPIXEL(x, y, c);
				}
			}
			else
			{
				error=dy-dx;
				while(y!=y2)
				{
					int e2=2*error;
					if(e2>-dx)
						error-=dx, y+=sy;
					if(e2<dy)
						error+=dy, x+=sx;
					SETPIXEL(x, y, c);
				}
			}
		}
	}
}
static int get_paretorank(IntArray pareto, Tests testinfo, int specialidx, int dec)//rank is in fixed precision
{
	TestInfo *curr=testinfo->data+specialidx;
	for(int k=(int)pareto->count-1;k>=0;--k)
	{
		int idx=pareto->data[k];
		TestInfo *p=testinfo->data+idx;
		if(idx==specialidx)
			return 2*k;
		if(curr->total.csize<p->total.csize||(curr->total.csize==p->total.csize&&(&curr->total.etime)[dec]<(&p->total.etime)[dec]))
			return 2*k+1;
	}
	return 0;
}
static void update_bounds(IntArray idxs, Tests testinfo, int rank, int dec, double *bounds)
{
	static const int reach=1;
	int start=rank-2*reach;
	int end=rank+2*reach;
	if(start<0)
		start=0;
	if(end>2*(int)idxs->count)
		end=2*(int)idxs->count;
	for(int ki=start;ki<=end;ki+=2)
	{
		int ka=(ki+1)>>1;
		int idx=idxs->data[ka<=(int)idxs->count-1?ka:idxs->count-1];
		TestInfo *info=testinfo->data+idx;
		double t=(&info->total.etime)[dec];
		if(bounds[0]>(double)info->total.csize)bounds[0]=(double)info->total.csize;
		if(bounds[1]<(double)info->total.csize)bounds[1]=(double)info->total.csize;
		if(bounds[2]>t)bounds[2]=t;
		if(bounds[3]<t)bounds[3]=t;
	}
}
static IntArray get_paretofront2(IntArray besttestidxs, Tests testinfo, int dec)
{
	IntArray sorted=0, pareto=0;
	int klowest=0, lowest=0;
	TestInfo *info=testinfo->data;
	
	ARRAY_ALLOC(sorted, besttestidxs->count, 0, 0, 1);
	ARRAY_ALLOC(pareto, besttestidxs->count, 0, 0, 1);
	
	/*
	staircase pareto front:
	
	besttestidxs is already sorted in descending size
	p0 = find fastest algorithm regardless of size (eg: copy)
	for each point
		p1 = immediate next point to left of p0
		for p2 from remaining points:
			if p2 better than p1 in both metrics
				p1 = p2
		p0 = add p1 to pareto front
	*/

	//p0 = find fastest algorithm
	klowest=0, lowest=0;
	for(int ki=1;ki<(int)besttestidxs->count;++ki)
	{
		int idx=besttestidxs->data[ki];
		double elowest	=(&info[lowest	].total.etime)[dec];
		double e2	=(&info[idx	].total.etime)[dec];
		if(elowest>e2||(elowest==e2&&info[lowest].total.csize<info[idx].total.csize))
			klowest=ki, lowest=idx;
	}
	pareto->data[pareto->count++]=lowest;

	for(int ki=klowest-1;ki>=0;)
	{
		int idx1=besttestidxs->data[ki];//get immediate next
		TestInfo *p1=info+idx1;
		int kbest=ki;
		for(int k2=ki-1;k2>=0;--k2)//find a better algorithm (in both metrics)
		{
			int idx2=besttestidxs->data[k2];
			TestInfo *p2=info+idx2;
			if(p2->total.csize<=p1->total.csize&&(&p2->total.etime)[dec]<(&p1->total.etime)[dec])
				kbest=k2, idx1=idx2, p1=p2;
		}
		pareto->data[pareto->count++]=idx1;//add to pareto front
		ki=kbest-1;
	}

	//[debug] print pareto front
#if defined _MSC_VER || 0
	printf("%s Pareto front:\n", dec?"Decode":"Encode");
	for(int ki=pareto->count-1;ki>=0;--ki)
	{
		int idx=pareto->data[ki];
		TestInfo *p=info+idx;
		printf("%2d %10lld %10.4lf %s\n"
			, ki
			, p->total.csize
			, (&p->total.etime)[dec]
			, p->codecname->data
		);
	}
	printf("\n");
#endif
	return pareto;
}
static double get_labelstep(double xmin, double xmax, double chartsize, int *ret_ndecimals)
{
	double step=(xmax-xmin)/chartsize;
	double prec=floor(log10(step));
	double order=pow(10, prec);
	double msd=step/order;
	if(msd<1.5)
		step=1*order;
	else if(msd<3)
		step=2*order;
	else if(msd<4.5)
		step=2.5*order;
	else if(msd<7)
		step=5*order;
	else
		step=10*order;
	if(ret_ndecimals)
	{
		prec=-prec;
		CLAMP2(prec, 0, 8);
		*ret_ndecimals=(int)prec;
	}
	return step;
}
static void print_pareto(IntArray besttestidxs, Tests testinfo, ptrdiff_t usize, int special)
{
	static const int cw=120, ch=60, leftmargin=2, rightmargin=12, ymargin=1;

	const int bsize=cw*ch;
	char *buf=0;
	double bounds[4]={0};//cmin, cmax, tmin, tmax
	double cmin=0, cmax=0, tmin=0, tmax=0;
	int erank=0, drank=0;

	IntArray penc=get_paretofront2(besttestidxs, testinfo, 0);
	IntArray pdec=get_paretofront2(besttestidxs, testinfo, 1);
	
	{
		int idx=besttestidxs->data[special];
		TestInfo *p=testinfo->data+idx;
		bounds[0]=(double)p->total.csize;
		bounds[1]=(double)p->total.csize;
		bounds[2]=p->total.etime;
		bounds[3]=p->total.etime;
		if(bounds[2]>p->total.dtime)bounds[2]=p->total.dtime;
		if(bounds[3]<p->total.dtime)bounds[3]=p->total.dtime;
		erank=get_paretorank(penc, testinfo, idx, 0);
		drank=get_paretorank(pdec, testinfo, idx, 1);
#if defined _MSC_VER || 0
		printf("Rank E%d.%d D%d.%d\n"
			, erank>>1, 5&-(erank&1)
			, drank>>1, 5&-(drank&1)
		);
#endif
		update_bounds(penc, testinfo, erank, 0, bounds);
		update_bounds(pdec, testinfo, drank, 1, bounds);
		cmin=bounds[0];
		cmax=bounds[1];
		tmin=bounds[2];
		tmax=bounds[3];
	}
	if(tmin<tmax&&cmin<cmax)
	{
		double px, py;//character width & height
		int xdecimals=0, ydecimals=0;
		double xstep=get_labelstep(cmin/(1024.*1024), cmax/(1024.*1024), cw/50., &xdecimals);
		double ystep=get_labelstep(tmin, tmax, ch/10., &ydecimals);

		buf=(char*)malloc(bsize);
		if(!buf)
		{
			CRASH("Alloc error");
			return;
		}
		memset(buf, ' ', bsize);
		px=(double)(cmax-cmin)/(cw-(leftmargin+rightmargin));
		py=(tmax-tmin)/(ch-ymargin*2);
		
		//draw grid
#if 1
		{
			double t0=floor(((0-ymargin)*py+tmin)/ystep)*ystep;
			for(int ky=0;ky<ch;++ky)
			{
				double t=floor(((ky+1-ymargin)*py+tmin)/ystep)*ystep;
				int ygrid=fabs(t-t0)>1e-6;
				double x0=floor(((0-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
				for(int kx=0;kx<cw;++kx)
				{
					double x=floor(((kx+1-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
					int xgrid=fabs(x-x0)>1e-6;
					if(xgrid&&ygrid)
						buf[cw*ky+kx]='+';
					else if(xgrid)
						buf[cw*ky+kx]='|';
					else if(ygrid)
						buf[cw*ky+kx]='-';
					//else
					//	buf[cw*ky+kx]=' ';
					x0=x;
				}
				t0=t;
			}
		}
#endif
		draw_curve(penc, testinfo, buf, cw, ch, px, cmin, py, tmin, leftmargin, ymargin, 0, '.');
		draw_curve(pdec, testinfo, buf, cw, ch, px, cmin, py, tmin, leftmargin, ymargin, 1, '*');

		for(int ki=0;ki<(int)besttestidxs->count;++ki)
		{
			int idx=besttestidxs->data[ki];
			TestInfo *info=testinfo->data+idx;
			int x=(int)((info->total.csize-cmin)/px)+leftmargin;
			int ye=(int)((info->total.etime-tmin)/py)+ymargin;
			int yd=(int)((info->total.dtime-tmin)/py)+ymargin;

			if(ye==yd)
				SETPIXEL(x, ye, 'X');
			else
			{
				SETPIXEL(x, ye, 'e');
				SETPIXEL(x, yd, 'd');
			}
			if(ki==special&&x+1<cw)
			{
				SETPIXEL(x+1, ye, '[');
				SETPIXEL(x+1, yd, '[');
			}
			for(int kx=0;kx<cw-(x+2)&&kx<info->codecname->count;++kx)
			{
				SETPIXEL(x+2+kx, ye, info->codecname->data[kx]);
				SETPIXEL(x+2+kx, yd, info->codecname->data[kx]);
			}
			if(ki==special&&(ptrdiff_t)x+2+info->codecname->count<cw)
			{
				SETPIXEL(x+2+info->codecname->count, ye, ']');
				SETPIXEL(x+2+info->codecname->count, yd, ']');
			}
		}

		printf("+");
		for(int k=0;k<cw;++k)
			printf("-");
		printf("+\n");

		{
			double t0=floor(((ch-0-ymargin)*py+tmin)/ystep)*ystep;
			int unitprinted=0;
			for(int ky=0;ky<ch;++ky)
			{
				double t=floor(((ch-(ky+1)-ymargin)*py+tmin)/ystep)*ystep;
				printf("|%.*s|", cw, buf+cw*(ch-1-ky));
				if(fabs(t-t0)>1e-6)
				{
					printf("%10.*lf", ydecimals, t0);
					if(!unitprinted)
					{
						printf(" sec");
						unitprinted=1;
					}
				}
				printf("\n");
				t0=t;
			}
		}

		printf("+");
		for(int k=0;k<cw;++k)
			printf("-");
		printf("+\n");

		{
			double x0=floor(((0-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
			printf(" ");//chart xborder
			for(int k=0;k<cw;)
			{
				double x=floor(((k+1-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
				int cond=fabs(x-x0)>1e-6;
				k+=printf("%c", cond?'^':' ');
				x0=x;
			}
			printf("\n");
		}
		{
			double x0=floor(((0+10-1-xdecimals-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
			printf(" ");//chart xborder
			for(int k=0;k<cw;)
			{
				double x=floor(((k+1+10-1-xdecimals-leftmargin)*px+cmin)/(1024.*1024*xstep))*xstep;
				int cond=fabs(x-x0)>1e-6;
				if(cond)
					k+=printf("%10.*lf", xdecimals, x0);
				else
					k+=printf(" ");
				x0=x;
			}
			printf(" MB\n");
		}

		free(buf);
	}
	free(penc);
	free(pdec);
}
int main(int argc, char **argv)
{
	static const char placeholdertag[]="*.";
	static const int placeholderlen=sizeof(placeholdertag)-1;
	const char *datasetname=0, *codecname=0;
	char programpath[MAX_PATH+1]={0};
	String currdir=0, tmpfn1=0, tmpfn2=0, srcpath=0, ext=0;
	Dataset uinfo=0;
	Tests testinfo=0;
	CommandFormat enccmd={0}, deccmd={0};
	int srcdstbounds[4]={0};
	char *srctitle=0, *dsttitle=0;
	IntArray besttestidxs=0;
	int titlecolwidth=0;

#ifdef __GNUC__
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:    \"%s\"  DATASET  CODEC\n"
			"You will be prompted to define DATASET and CODEC.\n"
			, argv[0]
		);
		return 0;
	}
	datasetname=argv[1];
	codecname=argv[2];
#else
	datasetname="lpcb";
	codecname="qlic2";
#endif

	//1. get program path
	{
		int k, len;

		len=GetModuleFileNameA(0, programpath, sizeof(programpath)-1);
		if(!len||len==sizeof(programpath)-1)
		{
			CRASH("GetModuleFileNameA", GetLastError());
			return 0;
		}
		k=len-1;
		for(;k>=0;--k)
		{
			if(programpath[k]=='/'||programpath[k]=='\\')
			{
				++k;
				break;
			}
		}
		programpath[k]=0;
	//	printf("%s\n", programpath);//
	}

	//2. get temp filenames
	{
		int len, val;

#ifdef __GNUC__
		len=GetTempPathA(sizeof(g_buf)-1, g_buf);
#elif defined _MSC_VER
		len=GetTempPath2A(sizeof(g_buf)-1, g_buf);
#endif
		if(!len)
		{
			CRASH("GetTempPath2A %d", GetLastError());
			return 0;
		}
		ARRAY_COPY(currdir, g_buf, len, 1, 1);
		ARRAY_ALLOC(tmpfn1, MAX_PATH+1, 0, 1, 1);
		val=GetTempFileNameA(g_buf, "t1_", 0, tmpfn1->data);
		if(!val)
		{
			CRASH("GetTempFileNameA %d", GetLastError());
			return 0;
		}
		tmpfn1->count=strlen(tmpfn1->data);
		ARRAY_ALLOC(tmpfn2, MAX_PATH+1, 0, 1, 1);
		val=GetTempFileNameA(g_buf, "t2_", 0, tmpfn2->data);
		if(!val)
		{
			CRASH("GetTempFileNameA %d", GetLastError());
			return 0;
		}
		tmpfn2->count=strlen(tmpfn2->data);
		ascii_deletefile(tmpfn1->data);//these were generated by GetTempFileNameA()
		ascii_deletefile(tmpfn2->data);
	}

	//3. get dataset & previous tests
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
	{
		String text=txt_load(g_buf);//null terminated
		if(text)
		{
			const char *start, *ptr, *end;
/*
zzzdata_DATASET.TXT
path	extension
"files:"
{usize  filetitle}[NFILES]
"tests:"
{codec  timestamp  {ctotal  etotal  dtotal  emax  dmax  RMSE[3]}  {csize  etime  dtime  emem  dmem  RMSE[3]}[NFILES]}

zzzcode_CODEC.TXT
enc_command_template
dec_command_template
*/
			start=(char*)text->data;
			ptr=start;
			end=start+text->count;
			srcpath=parse_str(start, &ptr, '\t', 0);
			skipspace(&ptr);
			ext=parse_str(start, &ptr, '\n', 0);
			skipspace(&ptr);
			
			ARRAY_ALLOC(uinfo, 0, free_uinfo, 1, 1);
			skiplabel(start, &ptr, "files:");
			while(!peeklabel(ptr, "tests:"))
			{
				UInfo *info;

				ARRAY_APPEND1(uinfo, 1);
				info=uinfo->data+uinfo->count++;
				//UInfo *info=(UInfo*)ARRAY_APPEND(uinfo, 0, 1, 1, 0);
				info->usize=parse_uint(&ptr);
				skipspace(&ptr);
				info->filename=parse_str(start, &ptr, '\n', 0);
				skipspace(&ptr);
			}
			if(!uinfo->count)
			{
				CRASH("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			//{
			//	time_t current=time(0);
			//	struct tm *local=localtime(&current);
			//	dst=local->tm_isdst;
			//}
			ARRAY_ALLOC(testinfo, 0, free_testinfo, 1, 1);
			skiplabel(start, &ptr, "tests:");
			skipspace(&ptr);
			while(ptr<end)
			{
				TestInfo *info=0;

				ARRAY_APPEND1(testinfo, 1);
				info=testinfo->data+testinfo->count++;
				info->codecname=parse_str(start, &ptr, '\t', 0);

				//YYYYmmdd_HHMMSS
				skipspace(&ptr);
				info->datetime.year=10*(10*(10*(ptr[0]-'0')+ptr[1]-'0')+ptr[2]-'0')+ptr[3]-'0';
				ptr+=4;
				info->datetime.month=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.day=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=3;//skip '_'
				info->datetime.hour=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.minute=10*(ptr[0]-'0')+ptr[1]-'0';
				ptr+=2;
				info->datetime.second=10*(ptr[0]-'0')+ptr[1]-'0';
				info->datetime.cs=0;
				ptr+=2;

				parse_cell(&ptr, &info->total);
				ARRAY_ALLOC(info->cells, uinfo->count, 0, 1, 1);
				info->cells->count=uinfo->count;
				for(int k=0;k<(int)uinfo->count;++k)
				{
					if(ptr>=end)
						CRASH("Unexpected EOF");
					ImCodecResult *cell=info->cells->data+k;
					parse_cell(&ptr, cell);
				}

				skipspace(&ptr);
			}
			free(text);
		}
		else
		{
			int len=0;

			for(;;)
			{
				printf("Define %s path:  ", datasetname);
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				srcpath=filter_path(g_buf, len, 1, 16);
				{
					struct stat info={0};
					stat(srcpath->data, &info);
					if(info.st_mode&S_IFDIR)
						break;
				}
				free(srcpath);
			}
			for(;;)
			{
				printf("Extension:  ");
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				if(!len)
					continue;
				{
					int valid=1;
					for(int k=0;k<len;++k)
					{
						if(g_buf[k]<=' ')
						{
							valid=0;
							break;
						}
					}
					if(valid)
						break;
				}
			}
			ARRAY_COPY(ext, g_buf, len, 1, 1);

			uinfo=get_uinfo((char*)srcpath->data, (char*)ext->data);
			if(!uinfo||!uinfo->count)
			{
				CRASH("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			ARRAY_ALLOC(testinfo, 0, free_testinfo, 1, 1);
		}
	}
	titlecolwidth=(int)strlen(codecname);
	for(int k=0;k<(int)uinfo->count;++k)//get filetitle column width
	{
		UInfo *info=uinfo->data+k;
		int start=0, end=0;
		get_filetitle(info->filename->data, (int)info->filename->count, &start, &end);
		int width=end-start;
		if(titlecolwidth<width)
			titlecolwidth=width;
	}
	ARRAY_ALLOC(besttestidxs, testinfo->count, 0, 1, 1);
	if(testinfo->count)
		besttestidxs->data[besttestidxs->count++]=0;
	for(int k=1;k<(int)testinfo->count;++k)//select best prev tests
	{
		TestInfo *test1=testinfo->data+k;
		double score1=CALCSCORE(test1->total.csize, test1->total.etime, test1->total.dtime);
		int encountered=0, bestsofar=1;
		for(int k2=0;k2<k;++k2)//compare with previous tests
		{
			TestInfo *test2=testinfo->data+k2;
			if(acme_strnimatch((char*)test2->codecname->data, test2->codecname->count, (char*)test1->codecname->data, test1->codecname->count))
			{
				double score2=CALCSCORE(test2->total.csize, test2->total.etime, test2->total.dtime);
				encountered=1;
				if(score2<score1)//new test loses
				{
					bestsofar=0;
					break;
				}
				for(int k3=0;k3<(int)besttestidxs->count;++k3)//new test surpasses previous test
				{
					int idx=besttestidxs->data[k3];
					if(idx==k2)
					{
						if(besttestidxs->count)
						{
							memmove(besttestidxs->data+k3,
								besttestidxs->data+k3+1,
								(besttestidxs->count-(k3+1LL))*sizeof(besttestidxs->data[0])
							);
							--besttestidxs->count;
						}
						//array_erase(&besttestidxs, k3, 1);
						break;
					}
				}
			}
		}
		if(!encountered||bestsofar)
			besttestidxs->data[besttestidxs->count++]=k;
	}
	for(int k=0;k<(int)besttestidxs->count-1;++k)//rank besttestidxs by csize (insertion sort for simplicity)
	{
		int *idx1=besttestidxs->data+k;
		TestInfo *test1=testinfo->data+*idx1;
		for(int k2=k+1;k2<(int)besttestidxs->count;++k2)
		{
			int *idx2=besttestidxs->data+k2;
			TestInfo *test2=testinfo->data+*idx2;
			if(test1->total.csize>test2->total.csize)//test2 is tighter
			{
				int temp;
				SWAPVAR(*idx1, *idx2, temp);
				test1=test2;
			}
		}
	}

	//4. get command templates
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzcode_%s.txt", programpath, codecname);
	{
		String text=txt_load(g_buf);//null terminated
		String enc0=0, dec0=0;
		int is_new=!text;
		if(text)
		{
			const char *start=(char*)text->data, *ptr=start;
			enccmd.format=parse_str(start, &ptr, '\n', 0);
			skipspace(&ptr);
			deccmd.format=parse_str(start, &ptr, '\n', 1);
		}
		else//new
		{
			int len=0;
			printf("\n");
			printf("  Use \"*.extension\" as filename placeholders.\n");
			printf("Define %s encode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
			ARRAY_COPY(enccmd.format, g_buf, len, 1, 1);
			printf("Define %s decode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
			ARRAY_COPY(deccmd.format, g_buf, len, 1, 1);
			ARRAY_COPY(enc0, enccmd.format->data, enccmd.format->count, 1, 1);
			ARRAY_COPY(dec0, deccmd.format->data, deccmd.format->count, 1, 1);
		}
		{
			int len=0;

			len=qualify_command(enccmd.format->data, (int)enccmd.format->count, g_buf);
			free(enccmd.format);
			ARRAY_COPY(enccmd.format, g_buf, len, 1, 1);

			len=qualify_command(deccmd.format->data, (int)deccmd.format->count, g_buf);
			free(deccmd.format);
			ARRAY_COPY(deccmd.format, g_buf, len, 1, 1);
		}
		int srcdstbounds2[4]={0};
		enccmd.bounds=parse_cmdformatext((char*)enccmd.format->data, (int)enccmd.format->count, (char*)ext->data, 1, srcdstbounds);
		deccmd.bounds=parse_cmdformatext((char*)deccmd.format->data, (int)deccmd.format->count, (char*)ext->data, 0, srcdstbounds2);
		int valid=acme_strnimatch(
			(char*)enccmd.format->data+srcdstbounds[2], (ptrdiff_t)srcdstbounds[3]-srcdstbounds[2],
			(char*)deccmd.format->data+srcdstbounds2[2], (ptrdiff_t)srcdstbounds2[3]-srcdstbounds2[2]
		);
		if(srcdstbounds2[0]!=-1)
		{
			valid+=acme_strnimatch(
				(char*)enccmd.format->data+srcdstbounds[0], (ptrdiff_t)srcdstbounds[1]-srcdstbounds[0],
				(char*)deccmd.format->data+srcdstbounds2[0], (ptrdiff_t)srcdstbounds2[1]-srcdstbounds2[0]
			)*2;
		}
		if(!(valid&1)||(srcdstbounds2[0]!=-1&&!(valid&2)))//dec dst can be omitted, eg 7zip
		{
			printf("Template \"*.extension\" mismatch\n");
			printf("Enc:  %s\n", (char*)enccmd.format->data);
			printf("Dec:  %s\n", (char*)deccmd.format->data);
			printf("\n");
			printf("Enc template requires source and destination placeholders\n");
			printf("Dec template requires source placeholder at least\n");
			printf("Source      extensions must match.\n");
			printf("Destination extensions must match.\n");
			if(!is_new)
				printf("\nCheck \"%s\"\n", g_buf);
			CRASH("");
		}
		if(is_new)
		{
			int printed=snprintf(g_buf, sizeof(g_buf)-1, "%s\n%s\n", enc0->data, dec0->data);
			char fn[128]={0};
			snprintf(fn, sizeof(fn)-1, "%szzzcode_%s.txt", programpath, codecname);
			txt_save(fn, g_buf, printed);
			free(enc0);
			free(dec0);
		}
		else
			free(text);
	}
	char *t1fn=0, *t2fn, *encline=0, *decline=0;
	//int t1len=0, t2len=0, enclen=0, declen=0;
	{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wrestrict"
#endif
		char *ptr=g_buf2, *end=g_buf2+sizeof(g_buf2)-1;
		t1fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn1->data, srcdstbounds[1]-(srcdstbounds[0]+placeholderlen), (char*)enccmd.format->data+srcdstbounds[0]+placeholderlen
		)+1;//skip null terminator
		srctitle=ptr;
		while(srctitle>t1fn&&srctitle[-1]!='/'&&srctitle[-1]!='\\')--srctitle;

		t2fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn2->data, srcdstbounds[3]-(srcdstbounds[2]+placeholderlen), (char*)enccmd.format->data+srcdstbounds[2]+placeholderlen
		)+1;
		dsttitle=ptr;
		while(dsttitle>t2fn&&dsttitle[-1]!='/'&&dsttitle[-1]!='\\')--dsttitle;
		
		encline=ptr;
		substitute_cmdplaceholders(&ptr, end, &enccmd, srctitle, dsttitle);
		decline=ptr;
		substitute_cmdplaceholders(&ptr, end, &deccmd, srctitle, dsttitle);
#if 0
		char *ptr=g_buf2, *end=g_buf2+sizeof(g_buf2)-1;
		t1fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn1->data, enccmd.srcbounds[1]-(enccmd.srcbounds[0]+3), (char*)enccmd.format->data+enccmd.srcbounds[0]+3
		)+1;//skip null terminator
		t1len=(int)(ptr-t1fn-1);

		t2fn=ptr;
		ptr+=snprintf(ptr, end-ptr, "%s.%.*s",
			(char*)tmpfn2->data, enccmd.dstbounds[1]-(enccmd.dstbounds[0]+3), (char*)enccmd.format->data+enccmd.dstbounds[0]+3
		)+1;
		t2len=(int)(ptr-t2fn-1);

		encline=ptr;
		if(enccmd.srcbounds[0]<enccmd.dstbounds[0])
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.srcbounds[0], (char*)enccmd.format->data,
				t1fn,
				enccmd.dstbounds[0]-enccmd.srcbounds[1], (char*)enccmd.format->data+enccmd.srcbounds[1],
				t2fn,
				(char*)enccmd.format->data+enccmd.dstbounds[1]
			)+1;
		}
		else
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.dstbounds[0], (char*)enccmd.format->data,
				t2fn,
				enccmd.srcbounds[0]-enccmd.dstbounds[1], (char*)enccmd.format->data+enccmd.dstbounds[1],
				t1fn,
				(char*)enccmd.format->data+enccmd.srcbounds[1]
			)+1;
		}
		enclen=(int)(ptr-encline-1);
		
		decline=ptr;
		if(deccmd.dstbounds[0]==-1)
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %s",
				deccmd.srcbounds[0], (char*)deccmd.format->data,
				t2fn,
				(char*)deccmd.format->data+deccmd.srcbounds[1]
			)+1;
		}
		else if(deccmd.srcbounds[0]<deccmd.dstbounds[0])
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.srcbounds[0], (char*)deccmd.format->data,
				t2fn,
				deccmd.dstbounds[0]-deccmd.srcbounds[1], (char*)deccmd.format->data+deccmd.srcbounds[1],
				t1fn,
				(char*)deccmd.format->data+deccmd.dstbounds[1]
			)+1;
		}
		else
		{
			ptr+=snprintf(ptr, end-ptr, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.dstbounds[0], (char*)deccmd.format->data,
				t1fn,
				deccmd.srcbounds[0]-deccmd.dstbounds[1], (char*)deccmd.format->data+deccmd.dstbounds[1],
				t2fn,
				(char*)deccmd.format->data+deccmd.srcbounds[1]
			)+1;
		}
		declen=(int)(ptr-decline-1);
		(void)t1len;
		(void)t2len;
		(void)enclen;
		(void)declen;
#endif
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
		printf("Temp filenames:\n");//
		printf("  \"%s\"\n", t1fn);
		printf("  \"%s\"\n", t2fn);
		printf("Working directory:\n");
		printf("  \"%s\"\n", (char*)currdir->data);
		printf("Commands:\n");
		printf("  %s\n", encline);//
		printf("  %s\n", decline);//
		printf("\n");
		if(ptr>end)
		{
			printf("\n\nsnprintf OOB  ptr %016zX > %016zX\n", (size_t)ptr, (size_t)end);
			CRASH("");
		}
	}

	//5. test		DON'T MODIFY g_buf2 BELOW THIS POINT
	ptrdiff_t usize=0;
	TestInfo *currtest=0;
	for(int k=0;k<(int)uinfo->count;++k)
	{
		UInfo *info=uinfo->data+k;
		usize+=info->usize;
	}
	ARRAY_APPEND1(testinfo, 1);
	++testinfo->count;
	currtest=testinfo->data+testinfo->count-1;
	ARRAY_COPY(currtest->codecname, codecname, strlen(codecname), 1, 1);
	ARRAY_ALLOC(currtest->cells, uinfo->count, 0, 1, 1);
	currtest->cells->count=uinfo->count;
	print_summary(besttestidxs, testinfo, usize, -1);
	printf("\n");
	print_currtimestamp("%Y-%m-%d_%H%M%S");
	printf("  ");
	colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xC060FF, "%s", datasetname);
	printf(" ");
	colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xFFC060, "%s", codecname);
	printf("  %d images\n", (int)uinfo->count);
//	printf("  %s %d  %s\n", datasetname, (int)uinfo->count, codecname);
	{
		time_t t=time(0);
		struct tm *date=localtime(&t);
		currtest->datetime.year		=date->tm_year+1900;
		currtest->datetime.month	=date->tm_mon+1;
		currtest->datetime.day		=date->tm_mday;
		currtest->datetime.hour		=date->tm_hour;
		currtest->datetime.minute	=date->tm_min;
		currtest->datetime.second	=date->tm_sec;
	}
	for(int k=0;k<(int)uinfo->count;++k)
	{
		struct stat finf={0};
		UInfo *info=uinfo->data+k;
		ImCodecResult *currcell=currtest->cells->data+k;
		int titlestart=0, titleend=0;

		get_filetitle(info->filename->data, (int)info->filename->count, &titlestart, &titleend);

		//Print 1:  idx filetitle usize			//print filetitle first in case of CRASH
		printf("%5d %-*.*s %10lld",
			k+1,
			titlecolwidth, titleend-titlestart, info->filename->data+titlestart,
			info->usize
		);

		int success=CopyFileA(info->filename->data, t1fn, 1);
		if(!success)
		{
			printf("Source:       %s\n", (char*)info->filename->data);
			printf("Destination:  %s\n", t1fn);
			CRASH("CopyFileA, %d", GetLastError());
			return 0;
		}

		{
			int64_t emem=0;
			exec_process2(encline, (char*)currdir->data, 0, &currcell->etime, &emem);
			currcell->emem=emem/(1024.*1024);
		}
		stat(t2fn, &finf);
		currcell->csize=finf.st_size;

		//Print 2:  -> csize etime
		printf(" -> %10lld B  %12.6lf", currcell->csize, currcell->etime);
		{
			if(!(finf.st_mode&S_IFREG)||!finf.st_size)
			{
				printf("\n");
				printf("Temp. file #1 is not found:\n");
				printf("  \"%s\"\n", t2fn);
				printf("The encode command was:\n");
				printf("  %s\n", encline);
				CRASH("");
			}
		}
		
		ascii_deletefile(t1fn);
		{
			int64_t dmem=0;
			exec_process2(decline, (char*)currdir->data, 0, &currcell->dtime, &dmem);
			currcell->dmem=dmem/(1024.*1024);
		}
		
		//Print 3:  dtime espeed dspeed emem dmem
		printf(
			" %12.6lf sec  %10.4lf %10.4lf MB/s %8.2lf %8.2lf MB  "
			, currcell->dtime
			, info->usize/(currcell->etime*1024*1024)
			, info->usize/(currcell->dtime*1024*1024)
			, currcell->emem
			, currcell->dmem
		);

		verify_files((char*)info->filename->data, t1fn, currcell->rmse);
		print_cmpresult(currcell->rmse);
		
		ascii_deletefile(t1fn);
		ascii_deletefile(t2fn);
		printf("\n");

		currtest->total.csize+=currcell->csize;
		currtest->total.etime+=currcell->etime;
		if(currtest->total.emem<currcell->emem)
			currtest->total.emem=currcell->emem;
		currtest->total.dtime+=currcell->dtime;
		if(currtest->total.dmem<currcell->dmem)
			currtest->total.dmem=currcell->dmem;
		currtest->total.rmse[0]+=currcell->rmse[0]*info->usize/3;
		currtest->total.rmse[1]+=currcell->rmse[1]*info->usize/3;
		currtest->total.rmse[2]+=currcell->rmse[2]*info->usize/3;
	}
	printf("\n");
	//print summary
	{
		currtest->total.rmse[0]/=usize/3.;
		currtest->total.rmse[1]/=usize/3.;
		currtest->total.rmse[2]/=usize/3.;
		printf(
			"%5d %*s %10lld -> %10lld B  %12.6lf %12.6lf sec  %10.4lf %10.4lf MB/s %8.2lf %8.2lf MB  "
			, (int)uinfo->count
			, titlecolwidth, ""
			, usize
			, currtest->total.csize
			, currtest->total.etime
			, currtest->total.dtime
			, usize/(currtest->total.etime*1024*1024)
			, usize/(currtest->total.dtime*1024*1024)
			, currtest->total.emem
			, currtest->total.dmem
		);
		print_cmpresult(currtest->total.rmse);
		printf("  BPD %10.6lf", (double)currtest->total.csize/usize*8);
		printf("\n");

		print_currtimestamp("%Y-%m-%d_%H%M%S");
		printf("  ");
		colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xC060FF, "%s", datasetname);
		printf(" ");
		colorprintf(COLORPRINTF_TXT_DEFAULT^0xFFFFFF, 0xFFC060, "%s", codecname);
		printf("  %lf sec\n", currtest->total.etime+currtest->total.dtime);
		
		printf("\n");
		int rank=0, *idx=0;
		for(;rank<(int)besttestidxs->count;++rank)
		{
			idx=besttestidxs->data+rank;
			TestInfo *test=testinfo->data+*idx;
			int won=0;
			if(currtest->total.csize==test->total.csize)
				won=currtest->total.etime+currtest->total.dtime*2<test->total.etime+test->total.dtime*2;
			else
				won=currtest->total.csize<test->total.csize;
			if(won)
				break;
		}
		ARRAY_APPEND1(besttestidxs, 1);
		memmove(
			besttestidxs->data+rank+1,
			besttestidxs->data+rank,
			(besttestidxs->count-rank)*sizeof(besttestidxs->data[0])
		);
		++besttestidxs->count;
		besttestidxs->data[rank]=(int)testinfo->count-1;
		print_summary(besttestidxs, testinfo, usize, rank);
		print_pareto(besttestidxs, testinfo, usize, rank);
	}

	//6. save		g_buf2 can be modified now
	{
		FILE *fdst=0;

		snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
		fdst=fopen(g_buf, "w");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", g_buf);
			return 0;
		}
		fprintf(fdst, "%s\t%s\nfiles:\n", srcpath->data, ext->data);
		for(int k=0;k<(int)uinfo->count;++k)
		{
			UInfo *info=uinfo->data+k;
			fprintf(fdst, "%10lld\t%s\n", info->usize, info->filename->data);
		}
		fprintf(fdst, "tests:\n");
		for(int k=0;k<(int)testinfo->count;++k)
		{
			TestInfo *info=testinfo->data+k;
			fprintf(fdst, "%-20s\t%04d%02d%02d_%02d%02d%02d"
				, (char*)info->codecname->data
				, info->datetime.year
				, info->datetime.month
				, info->datetime.day
				, info->datetime.hour
				, info->datetime.minute
				, info->datetime.second
			);
			write_cell(fdst, &info->total);
			for(int k2=0;k2<(int)info->cells->count;++k2)
				write_cell(fdst, info->cells->data+k2);
			fprintf(fdst, "\n");
		}
		fclose(fdst);
	}

	//7. free
	free(enccmd.format);
	free(deccmd.format);
	free(srcpath);
	free(ext);
	ARRAY_FREE(uinfo);
	ARRAY_FREE(testinfo);
	free(tmpfn1);
	free(tmpfn2);
	free(besttestidxs);
	free(currdir);
	return 0;
}
