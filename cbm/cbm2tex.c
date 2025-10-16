#include"util.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<time.h>
#include<Windows.h>
#include<Psapi.h>
#include<tlhelp32.h>
#include<immintrin.h>
static const char file[]=__FILE__;


typedef union _DateTime
{
	struct
	{
		uint8_t cs, second, minute, hour, day, month;
		uint16_t year;
	};
	uint64_t timestamp;
} DateTime;

static int acme_getline(char *buf, int len, FILE *f)
{
	memset(buf, '\n', len);
	fgets(buf, len, f);
	int k=0;
	for(;k<len-1&&buf[k]!='\n';++k);
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

//GDCC score
#define CALCSCORE(CSIZE, ENC, DEC) ((CSIZE)*(1./(1024*1024))+(ENC)+(DEC)*2)

typedef struct _UInfo
{
	long long usize;
	ArrayHandle filename;
} UInfo;
static void free_uinfo(void *p)
{
	UInfo *info=(UInfo*)p;
	array_free(&info->filename);
}
typedef struct _CellInfo
{
	long long csize;
	double etime, dtime;
	long long emem, dmem;
} CellInfo;
typedef struct _TestInfo
{
	ArrayHandle codecname;
	DateTime datetime;
	CellInfo total;
	ArrayHandle cells;
} TestInfo;
static void free_testinfo(void *p)
{
	TestInfo *info=(TestInfo*)p;
	array_free(&info->codecname);
	array_free(&info->cells);
}
typedef struct _Range
{
	int issrc, start, end;
} Range;
typedef struct _CommandFormat
{
	ArrayHandle format;
	ArrayHandle bounds;//Range
} CommandFormat;
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
		LOG_ERROR("");
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
		LOG_ERROR("");
	}
	*ptr=ptr2;
}
static ArrayHandle parse_str(const char *start, const char **ptr, char delim, int allow_eof)//delim must be one of the whitespace characters
{
	ArrayHandle str=0;
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
		LOG_ERROR("");
	}
skipcheck:
	{
		const char *ptr3=ptr2;
		while(ptr3>*ptr&&isspace(*(ptr3-1)))--ptr3;
		STR_COPY(str, *ptr, ptr3-*ptr);
	}
	//skipspace(&ptr2);
	//while(*ptr2&&isspace(*ptr2++));
	//ptr2-=*ptr2!=0;
	*ptr=ptr2;
	return str;
}
static long long parse_uint(const char **ptr)
{
	const char *ptr2=*ptr;
	long long val=0;
	skipspace(&ptr2);
	while((unsigned)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	*ptr=ptr2;
	return val;
}
static double parse_float(const char **ptr)
{
	const char *ptr2=*ptr;
	double val=0;
	skipspace(&ptr2);
	while((unsigned)(*ptr2-'0')<10)
		val=10*val+*ptr2++-'0';
	if(*ptr2=='.')
	{
		++ptr2;
		double p=0.1;
		while((unsigned)(*ptr2-'0')<10)
		{
			val+=(*ptr2++-'0')*p;
			p*=0.1;
		}
	}
	*ptr=ptr2;
	return val;
}
static void parse_cell(const char **ptr, CellInfo *cell)
{
	//{csize  etime  dtime  emem  dmem}
	cell->csize=parse_uint(ptr);
	cell->etime=parse_float(ptr);
	cell->dtime=parse_float(ptr);
	cell->emem=parse_uint(ptr);
	cell->dmem=parse_uint(ptr);
}
static const char* get_extension(const char *filename, ptrdiff_t len)//excludes the dot
{
	ptrdiff_t idx;

	idx=acme_strrchr(filename, len, '.');
	if(idx==-1)
		return 0;
	return filename+idx+1;
}
static ArrayHandle get_uinfo(const char *path, const char *ext)//extension without '.'
{
	//prepare searchpath
	ArrayHandle searchpath=filter_path(path, -1);
	STR_APPEND(searchpath, "*", 1, 1);

	WIN32_FIND_DATAA data={0};
	void *hSearch=FindFirstFileA((char*)searchpath->data, &data);//skip .
	if(hSearch==INVALID_HANDLE_VALUE)
		return 0;
	FindNextFileA(hSearch, &data);//skip ..

	ArrayHandle uinfo;
	ARRAY_ALLOC(UInfo, uinfo, 0, 0, 0, free_uinfo);
	while(FindNextFileA(hSearch, &data))
	{
		ptrdiff_t len=strlen(data.cFileName);
		const char *curr_ext=get_extension(data.cFileName, len);
		if(!(data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)&&!acme_stricmp(curr_ext, ext))
		{
			UInfo *info=(UInfo*)ARRAY_APPEND(uinfo, 0, 1, 1, 0);
			info->usize=(long long)data.nFileSizeHigh<<32|data.nFileSizeLow;
			STR_COPY(info->filename, searchpath->data, searchpath->count-1);
			STR_APPEND(info->filename, data.cFileName, len, 1);
		}
	}
	int success=FindClose(hSearch);
	if(!success)
	{
		SYSTEMERROR("FindClose");
		return 0;
	}
	array_free(&searchpath);
	return uinfo;
}
int main(int argc, char **argv)
{
	const char *datasetname=0;
//#ifndef _DEBUG
#ifdef __GNUC__
	if(argc!=2)
	{
		printf("Usage:    \"%s\"  DATASET\n", argv[0]);
		return 0;
	}
	datasetname=argv[1];
#else
	datasetname="div2k";
#endif
	char programpath[MAX_PATH+1]={0};
	ArrayHandle currdir=0;
	ArrayHandle srcpath=0, ext=0, uinfo=0, testinfo=0;

	ArrayHandle besttestidxs=0;

	//1. get program path
	{
		int len=GetModuleFileNameA(0, programpath, sizeof(programpath)-1);
		if(!len||len==sizeof(programpath)-1)
		{
			SYSTEMERROR("GetModuleFileNameA");
			return 0;
		}
		int k=len-1;
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

	//3. get dataset & previous tests
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
	{
		ArrayHandle text=load_file(g_buf, 0, 16, 0);//null terminated
		if(text)
		{
/*
dataset_DATASET.TXT
path	extension
"files:"
{usize  filetitle}[NFILES]
"tests:"
{codec  timestamp  {ctotal  etotal  dtotal  emax  dmax}  {csize  etime  dtime  emem  dmem}[NFILES]}

codec_CODEC.TXT
enc command template
dec command template
*/
			//int dst=0;
			const char *start=(char*)text->data, *ptr=start, *end=start+text->count;
			srcpath=parse_str(start, &ptr, '\t', 0);
			skipspace(&ptr);
			ext=parse_str(start, &ptr, '\n', 0);
			skipspace(&ptr);

			ARRAY_ALLOC(UInfo, uinfo, 0, 0, 0, free_uinfo);
			skiplabel(start, &ptr, "files:");
			while(!peeklabel(ptr, "tests:"))
			{
				UInfo *info=(UInfo*)ARRAY_APPEND(uinfo, 0, 1, 1, 0);
				info->usize=parse_uint(&ptr);
				skipspace(&ptr);
				info->filename=parse_str(start, &ptr, '\n', 0);
				skipspace(&ptr);
			}
			if(!uinfo->count)
			{
				LOG_ERROR("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			//{
			//	time_t current=time(0);
			//	struct tm *local=localtime(&current);
			//	dst=local->tm_isdst;
			//}
			ARRAY_ALLOC(TestInfo, testinfo, 0, 0, 0, free_testinfo);
			skiplabel(start, &ptr, "tests:");
			skipspace(&ptr);
			while(ptr<end)
			{
				TestInfo *info=(TestInfo*)ARRAY_APPEND(testinfo, 0, 1, 1, 0);
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
				ARRAY_ALLOC(CellInfo, info->cells, 0, 0, 0, 0);
				for(int k=0;k<(int)uinfo->count;++k)
				{
					if(ptr>=end)
						LOG_ERROR("Unexpected EOF");
					CellInfo *cell=(CellInfo*)ARRAY_APPEND(info->cells, 0, 1, 1, 0);
					parse_cell(&ptr, cell);
				}

				skipspace(&ptr);
			}
			array_free(&text);
		}
		else
		{
			ptrdiff_t size=0;
			int len=0;
			for(;;)
			{
				printf("Define %s path:  ", datasetname);
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				srcpath=filter_path(g_buf, len);
				size=get_filesize((char*)srcpath->data);
				if(size==FSIZE_FOLDER)
					break;
				array_free(&srcpath);
			}
			for(;;)
			{
				printf("Extension:  ");
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				if(!len)
					continue;
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
			STR_COPY(ext, g_buf, len);

			uinfo=get_uinfo((char*)srcpath->data, (char*)ext->data);
			if(!uinfo||!uinfo->count)
			{
				LOG_ERROR("No %s files in \"%s\"", ext->data, srcpath->data);
				return 0;
			}
			ARRAY_ALLOC(TestInfo, testinfo, 0, 0, 0, free_testinfo);
		}
	}
	int titlecolwidth=0;
	for(int k=0;k<(int)uinfo->count;++k)//get filetitle column width
	{
		UInfo *info=(UInfo*)array_at(&uinfo, k);
		int start=0, end=0;
		get_filetitle((char*)info->filename->data, (int)info->filename->count, &start, &end);
		int width=end-start;
		if(titlecolwidth<width)
			titlecolwidth=width;
	}
	ARRAY_ALLOC(int, besttestidxs, 0, 0, testinfo->count, 0);
	{
		int *idx=(int*)ARRAY_APPEND(besttestidxs, 0, 1, 1, 0);
		*idx=0;
	}
	for(int k=1;k<(int)testinfo->count;++k)//select best prev tests
	{
		TestInfo *test1=(TestInfo*)array_at(&testinfo, k);
		double score1=CALCSCORE(test1->total.csize, test1->total.etime, test1->total.dtime);
		int encountered=0, bestsofar=1;
		for(int k2=0;k2<k;++k2)//compare with previous tests
		{
			TestInfo *test2=(TestInfo*)array_at(&testinfo, k2);
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
					int *idx=(int*)array_at(&besttestidxs, k3);
					if(*idx==k2)
					{
						array_erase(&besttestidxs, k3, 1);
						break;
					}
				}
			}
		}
		if(!encountered||bestsofar)
		{
			int *idx=(int*)ARRAY_APPEND(besttestidxs, 0, 1, 1, 0);
			*idx=k;
		}
	}
	for(int k=0;k<(int)besttestidxs->count-1;++k)//rank besttestidxs by csize (insertion sort for simplicity)
	{
		int *idx1=(int*)array_at(&besttestidxs, k);
		TestInfo *test1=(TestInfo*)array_at(&testinfo, *idx1);
		for(int k2=k+1;k2<(int)besttestidxs->count;++k2)
		{
			int *idx2=(int*)array_at(&besttestidxs, k2);
			TestInfo *test2=(TestInfo*)array_at(&testinfo, *idx2);
			if(test1->total.csize>test2->total.csize)//test2 is tighter
			{
				int temp;
				SWAPVAR(*idx1, *idx2, temp);
				test1=test2;
			}
		}
	}

	//6. save
	{
		snprintf(g_buf, sizeof(g_buf)-1, "%sfig_%s.tex", programpath, datasetname);
		FILE *fdst=fopen(g_buf, "w");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", g_buf);
			return 0;
		}
		fprintf(fdst, "\tx\ty\tlabel\n");
		for(int k2=0;k2<(int)besttestidxs->count;++k2)//enc
		{
			int *idx=(int*)array_at(&besttestidxs, k2);
			TestInfo *test=(TestInfo*)array_at(&testinfo, *idx);
			fprintf(fdst, "\t%8.3lf\t%11.6lf\t%s\n"
				, (double)test->total.csize/(1024*1024)
				, test->total.etime
				, (char*)test->codecname->data
			);
		}
		fprintf(fdst, "\n\n");

		fprintf(fdst, "\tx\ty\tlabel\n");
		for(int k2=0;k2<(int)besttestidxs->count;++k2)//dec
		{
			int *idx=(int*)array_at(&besttestidxs, k2);
			TestInfo *test=(TestInfo*)array_at(&testinfo, *idx);
			fprintf(fdst, "\t%8.3lf\t%11.6lf\t%s\n"
				, (double)test->total.csize/(1024*1024)
				, test->total.dtime
				, (char*)test->codecname->data
			);
		}
		fprintf(fdst, "\n\n");
		fclose(fdst);
	}

	//7. free
	array_free(&srcpath);
	array_free(&ext);
	array_free(&uinfo);
	array_free(&testinfo);
	array_free(&besttestidxs);
	array_free(&currdir);
	return 0;
}