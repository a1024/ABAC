#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<time.h>
#include<Windows.h>
#include<Psapi.h>
static const char file[]=__FILE__;

static int acme_getline(char *buf, int len, FILE *f)
{
	memset(buf, '\n', len);
	fgets(buf, len, f);
	int k=0;
	for(;k<len&&buf[k]!='\n';++k);
	buf[k]=0;
	return k;
}
static int acme_strnimatch(const char *s1, ptrdiff_t len1, const char *s2, ptrdiff_t len2)//1: match	FIXME return ASCII order & index of first difference
{
	if(len1!=len2)
		return 0;
	const char *end1=s1+len1;
	while(s1<end1&&tolower(*s1++)==tolower(*s2++));
	return s1==end1;
}

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
	time_t timestamp;
	CellInfo total;
	ArrayHandle cells;
} TestInfo;
static void free_testinfo(void *p)
{
	TestInfo *info=(TestInfo*)p;
	array_free(&info->codecname);
	array_free(&info->cells);
}
typedef struct _CommandFormat
{
	ArrayHandle format;
	int srcbounds[2], dstbounds[2];
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
	while(*ptr2&&isspace(*ptr2++));
	ptr2-=*ptr2!=0;
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
	STR_COPY(str, *ptr, ptr2-*ptr);
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
	//while(*ptr2&&isspace(*ptr2++));
	//ptr2-=*ptr2!=0;
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
	//while(*ptr2&&isspace(*ptr2++));
	//ptr2-=*ptr2!=0;
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
static int strimatch(const char *text, const char *label)//1: match
{
	while(*label&&tolower(*text++)==tolower(*label++));
	return !*label;
}
static int parse_ext(const char *text, int len, int *bounds, const char *targetext, const int *excludebounds)
{
	const int ndashes=3;
	int found=0;
	int count=0;
	const char *ptr=text+len;
	while(--ptr>=text)//search backward, in case the path itself contains dashes
	{
		char c=*ptr;
		if(c=='-')
		{
			++count;
			if(count==ndashes)
			{
				int idx=(int)(ptr-text);
				if(excludebounds&&idx==excludebounds[0])
					count=0;
				else if(!targetext||strimatch(ptr+ndashes, targetext))
				{
					found=1;
					bounds[0]=idx;
					idx+=ndashes;
					while(text[idx]&&!isspace(text[idx++]));
					idx-=text[idx]!=0;
					bounds[1]=idx;
					break;
				}
			}
		}
		else
			count=0;
	}
	return found;
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
		LOG_ERROR("FindClose GetLastError %d", (int)GetLastError());
		return 0;
	}
	array_free(&searchpath);
	return uinfo;
}
static void write_cell(FILE *fdst, CellInfo *cell)//tabs for excel
{
	fprintf(fdst, "\t%10lld\t%12.6lf\t%12.6lf\t%10lld\t%10lld", cell->csize, cell->etime, cell->dtime, cell->emem, cell->dmem);
}
#if 0
static void verify_command(const char *cmd)
{
	if(!cmd)
	{
		LOG_ERROR("");
		return;
	}
	const char *ptr=cmd;
	while(*ptr&&!isspace(*ptr++));//assume no spaces
	ptr-=*ptr!=0;
	char fn1[MAX_PATH+1]={0}, fn2[MAX_PATH+1]={0};
	memcpy(fn1, cmd, ptr-cmd);
	int len=GetFullPathNameA(fn1, sizeof(fn2)-1, fn2, 0);
	if(!len)
	{
		LOG_ERROR("GetFullPathNameA GetLastError %d", (int)GetLastError());
		return;
	}
	ptrdiff_t size=get_filesize(fn2);
	if(size<1)
	{
		printf("Program is inaccessible:\n");
		printf("  \"%s\"\n", fn1);
		printf("  \"%s\"\n", fn2);
		printf("  size %td bytes", size);
		LOG_ERROR("");
	}
}
#endif
static void qualify_command(ArrayHandle *cmd)
{
	const char *ptr=(char*)cmd[0]->data;
	while(*ptr&&!isspace(*ptr++));//assume no spaces
	ptr-=*ptr!=0;
	char fn1[MAX_PATH+1]={0}, fn2[MAX_PATH+2]={0};
	int len1=(int)(ptr-(char*)cmd[0]->data);
	memcpy(fn1, cmd[0]->data, len1);
	if(!strimatch(fn1+len1-4, ".exe"))
	{
		memcpy(fn1+len1, ".exe", 4);
		len1+=4;
	}
	int len2=SearchPathA(0, fn1, 0, sizeof(fn2)-2, fn2+1, 0);
	//int len=GetFullPathNameA(fn1, sizeof(fn2)-2, fn2+1, 0);
	if(!len2)
	{
		LOG_ERROR("GetFullPathNameA GetLastError %d", (int)GetLastError());
		return;
	}
	ptrdiff_t size=get_filesize(fn2+1);
	if(size<1)
	{
		printf("Program is inaccessible:\n");
		printf("  \"%s\"\n", fn1);
		printf("  \"%s\"\n", fn2+1);
		printf("  size %td bytes", size);
		LOG_ERROR("");
	}
	len2+=2;
	fn2[0]='\"';
	fn2[len2-1]='\"';
//	char *cmd2=(char*)
	array_replace(cmd, 0, ptr-(char*)cmd[0]->data, fn2, len2, 1, 1);
//	printf("  %s\n", cmd2);//
}
static void exec_process(char *cmd, int loud, double *elapsed, long long *maxmem)
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
			LOG_ERROR("CreateFileA GetLastError %d", (int)GetLastError());
			return;
		}
		si.hStdError=CreateFileA("NUL", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdError==INVALID_HANDLE_VALUE)
		{
			LOG_ERROR("CreateFileA GetLastError %d", (int)GetLastError());
			return;
		}
		si.hStdInput=CreateFileA("NUL", GENERIC_READ, 0, NULL, CREATE_ALWAYS, 0, NULL);
		if(si.hStdInput==INVALID_HANDLE_VALUE)
		{
			LOG_ERROR("CreateFileA GetLastError %d", (int)GetLastError());
			return;
		}
	}
	success=CreateProcessA(0, cmd, 0, 0, 0, CREATE_SUSPENDED, 0, 0, &si, &pi);
	if(!success)
	{
		LOG_ERROR("CreateProcessA GetLastError %d", (int)GetLastError());
		return;
	}
	ptrdiff_t memusage=0;
	double t=time_sec();
	int suspendcount=ResumeThread(pi.hThread);
	if(suspendcount==(DWORD)-1)
	{
		LOG_ERROR("ResumeThread GetLastError %d", (int)GetLastError());
		return;
	}
	while(WaitForSingleObject(pi.hProcess, 10)==WAIT_TIMEOUT)
	{
		PROCESS_MEMORY_COUNTERS pmc={0};
		pmc.cb=sizeof(pmc);
		success=GetProcessMemoryInfo(pi.hProcess, &pmc, sizeof(pmc));
		if(!success)
		{
			LOG_ERROR("GetProcessMemoryInfo GetLastError %d", (int)GetLastError());
			return;
		}
		if(memusage<(ptrdiff_t)pmc.WorkingSetSize)
			memusage=pmc.WorkingSetSize;
	}
	if(elapsed)*elapsed=time_sec()-t;
	if(maxmem)*maxmem=memusage;
	success=CloseHandle(pi.hProcess);
	if(!success)
	{
		LOG_ERROR("CloseHandle GetLastError %d", (int)GetLastError());
		return;
	}
	success=CloseHandle(pi.hThread);
	if(!success)
	{
		LOG_ERROR("CloseHandle GetLastError %d", (int)GetLastError());
		return;
	}
	if(!loud)
	{
		success=CloseHandle(si.hStdOutput);
		if(!success)
		{
			LOG_ERROR("CloseHandle GetLastError %d", (int)GetLastError());
			return;
		}
		success=CloseHandle(si.hStdError);
		if(!success)
		{
			LOG_ERROR("CloseHandle GetLastError %d", (int)GetLastError());
			return;
		}
		success=CloseHandle(si.hStdInput);
		if(!success)
		{
			LOG_ERROR("CloseHandle GetLastError %d", (int)GetLastError());
			return;
		}
	}
}
static void print_usage(const char *argv0)
{
	printf(
		"Usage:    %s  DATASET  CODEC\n"
		"Example:  %s  div2k    libjxl\n"
		"You will be prompted to define the DATASET and CODEC macros.\n"
		, argv0, argv0
	);
}
int main(int argc, char **argv)
{
	const char *datasetname=0, *codecname=0;
//#ifndef _DEBUG
#ifdef __GNUC__
	if(argc!=3)
	{
		print_usage(argv[0]);
		return 0;
	}
	datasetname=argv[1];
	codecname=argv[2];
#else
	datasetname="div2k";
	codecname="halic072fast";
#endif
	char programpath[MAX_PATH+1]={0};
	ArrayHandle tmpfn1=0, tmpfn2=0;
	ArrayHandle srcpath=0, ext=0, uinfo=0, testinfo=0;
	CommandFormat enccmd={0}, deccmd={0};

	//1. get program path
	{
		int len=GetModuleFileNameA(0, programpath, sizeof(programpath)-1);
		if(!len||len==sizeof(programpath)-1)
		{
			LOG_ERROR("GetModuleFileNameA GetLastError %s", (int)GetLastError());
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
			LOG_ERROR("GetTempPath2A GetLastError %d", (int)GetLastError());
			return 0;
		}
		STR_ALLOC(tmpfn1, MAX_PATH+1);
		val=GetTempFileNameA(g_buf, "t1_", 0, (char*)tmpfn1->data);
		if(!val)
		{
			LOG_ERROR("GetTempFileNameA GetLastError %d", (int)GetLastError());
			return 0;
		}
		tmpfn1->count=strlen((char*)tmpfn1->data);
		STR_ALLOC(tmpfn2, MAX_PATH+1);
		val=GetTempFileNameA(g_buf, "t2_", 0, (char*)tmpfn2->data);
		if(!val)
		{
			LOG_ERROR("GetTempFileNameA GetLastError %d", (int)GetLastError());
			return 0;
		}
		tmpfn2->count=strlen((char*)tmpfn2->data);
		printf("Temp filenames:\n");
		printf("  %s\n", (char*)tmpfn1->data);
		printf("  %s\n", (char*)tmpfn2->data);
		printf("\n");
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

			ARRAY_ALLOC(TestInfo, testinfo, 0, 0, 0, free_testinfo);
			skiplabel(start, &ptr, "tests:");
			skipspace(&ptr);
			while(ptr<end)
			{
				TestInfo *info=(TestInfo*)ARRAY_APPEND(testinfo, 0, 1, 1, 0);
				info->codecname=parse_str(start, &ptr, '\t', 0);

				//YYYYmmdd_HHMMSS
				int d=(int)parse_uint(&ptr);
				ptr+=*ptr=='_';
				int t=(int)parse_uint(&ptr);
#ifdef _DEBUG
				int d0=d, t0=t;
				(void)d0;
				(void)t0;
#endif
				struct tm date={0};
				date.tm_sec=t%100;
				t/=100;
				date.tm_min=t%100;
				t/=100;
				date.tm_hour=t;
				date.tm_mday=d%100;
				d/=100;
				date.tm_mon=d%100-1;
				d/=100;
				date.tm_year=d-1900;
				if((unsigned)date.tm_sec>=60
					||(unsigned)date.tm_min>=60
					||(unsigned)date.tm_hour>=24
					||(unsigned)(date.tm_mday-1)>=(unsigned)(31-1)
					||(unsigned)date.tm_mon>=12
					||date.tm_year<2025-1900
				)
					LOG_ERROR("Invalid timestamp %04d-%02d-%02d_%02d-%02d-%02d", date.tm_year+1900, date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
				info->timestamp=mktime(&date);

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
			do
			{
				printf("Define %s path:  ", datasetname);
				len=acme_getline(g_buf, sizeof(g_buf), stdin);
				srcpath=filter_path(g_buf, len);
				size=get_filesize((char*)srcpath->data);
			}while(size);
			printf("Extension:  ");
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
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

	//4. get command templates
	snprintf(g_buf, sizeof(g_buf)-1, "%szzzcode_%s.txt", programpath, codecname);
	{
		ArrayHandle text=load_file(g_buf, 0, 16, 0);//null terminated
		ArrayHandle enc0=0, dec0=0;
		int is_new=!text;
		if(text)
		{
			const char *start=(char*)text->data, *ptr=start;
			enccmd.format=parse_str(start, &ptr, '\n', 0);
		//	verify_command((char*)enccmd.format->data);
			skipspace(&ptr);
			deccmd.format=parse_str(start, &ptr, '\n', 1);
		//	verify_command((char*)deccmd.format->data);
			array_free(&text);
		}
		else
		{
			int len=0;
			printf("\n");
			printf("  Use \"---extension\" as filename placeholders.\n");
			printf("Define %s encode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
		//	verify_command(g_buf);
			STR_COPY(enccmd.format, g_buf, len);
			printf("Define %s decode:  ", codecname);
			len=acme_getline(g_buf, sizeof(g_buf), stdin);
		//	verify_command(g_buf);
			STR_COPY(deccmd.format, g_buf, len);
			STR_COPY(enc0, enccmd.format->data, enccmd.format->count);
			STR_COPY(dec0, deccmd.format->data, deccmd.format->count);
		}
		qualify_command(&enccmd.format);
		qualify_command(&deccmd.format);
		enccmd.srcbounds[0]=-1;
		enccmd.srcbounds[1]=-1;
		deccmd.srcbounds[0]=-1;
		deccmd.srcbounds[1]=-1;
		parse_ext((char*)enccmd.format->data, (int)enccmd.format->count, enccmd.srcbounds, (char*)ext->data, 0);
		parse_ext((char*)enccmd.format->data, (int)enccmd.format->count, enccmd.dstbounds, 0, enccmd.srcbounds);
		parse_ext((char*)deccmd.format->data, (int)deccmd.format->count, deccmd.dstbounds, (char*)ext->data, 0);
		parse_ext((char*)deccmd.format->data, (int)deccmd.format->count, deccmd.srcbounds, 0, deccmd.dstbounds);
		if(
			enccmd.srcbounds[0]==-1
		||	enccmd.srcbounds[1]==-1
		||	deccmd.srcbounds[0]==-1
		||	deccmd.srcbounds[1]==-1
		||	!acme_strnimatch((char*)enccmd.format->data+enccmd.srcbounds[0]+3, enccmd.srcbounds[1]-((ptrdiff_t)enccmd.srcbounds[0]+3), (char*)ext->data, ext->count)
		||	!acme_strnimatch((char*)deccmd.format->data+deccmd.dstbounds[0]+3, deccmd.dstbounds[1]-((ptrdiff_t)deccmd.dstbounds[0]+3), (char*)ext->data, ext->count)
		||	!acme_strnimatch(
				(char*)enccmd.format->data+enccmd.dstbounds[0], enccmd.dstbounds[1]-(ptrdiff_t)enccmd.dstbounds[0],
				(char*)deccmd.format->data+deccmd.srcbounds[0], deccmd.srcbounds[1]-(ptrdiff_t)deccmd.srcbounds[0]
			)
		)
		{
			printf("\n");
			printf("Invalid templates\n");
			printf("Each template should contain 2 placeholders of the form ---extension\n");
			printf("Example:\n");
			printf("  program ---ppm ---ext\n");
			printf("  program ---ext ---ppm\n");
			printf("\n");
			printf("Enc template:  %s\n", (char*)enccmd.format->data);
			printf("Dec template:  %s\n", (char*)deccmd.format->data);
			printf("\n");
			printf("Enc src \"%.*s\" [%d ~ %d]\n", enccmd.srcbounds[1]-enccmd.srcbounds[0], (char*)enccmd.format->data+enccmd.srcbounds[0], enccmd.srcbounds[0], enccmd.srcbounds[1]);
			printf("Enc dst \"%.*s\" [%d ~ %d]\n", enccmd.dstbounds[1]-enccmd.dstbounds[0], (char*)enccmd.format->data+enccmd.dstbounds[0], enccmd.dstbounds[0], enccmd.dstbounds[1]);
			printf("Dec src \"%.*s\" [%d ~ %d]\n", deccmd.srcbounds[1]-deccmd.srcbounds[0], (char*)deccmd.format->data+deccmd.srcbounds[0], deccmd.srcbounds[0], deccmd.srcbounds[1]);
			printf("Dec dst \"%.*s\" [%d ~ %d]\n", deccmd.dstbounds[1]-deccmd.dstbounds[0], (char*)deccmd.format->data+deccmd.dstbounds[0], deccmd.dstbounds[0], deccmd.dstbounds[1]);
			LOG_ERROR("");
		}
		if(is_new)
		{
			int printed=snprintf(g_buf, sizeof(g_buf)-1, "%s\n%s\n", enc0->data, dec0->data);
			char fn[128]={0};
			snprintf(fn, sizeof(fn)-1, "%szzzcode_%s.txt", programpath, codecname);
			save_file(fn, g_buf, printed, 0);
			array_free(&enc0);
			array_free(&dec0);
		}
	}

	//5. test
	ptrdiff_t usize=0;
	TestInfo *currtest=(TestInfo*)ARRAY_APPEND(testinfo, 0, 1, 1, 0);
	STR_COPY(currtest->codecname, codecname, strlen(codecname));
	ARRAY_ALLOC(CellInfo, currtest->cells, 0, uinfo->count, 0, 0);
	print_timestamp("%Y-%m-%d_%H%M%S\n");
	currtest->timestamp=time(0);
	for(int k=0;k<(int)uinfo->count;++k)
	{
		UInfo *info=(UInfo*)array_at(&uinfo, k);
		CellInfo *currcell=(CellInfo*)array_at(&currtest->cells, k);
		int printed=0;
		if(enccmd.srcbounds[0]<enccmd.dstbounds[0])
		{
			printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.srcbounds[0], (char*)enccmd.format->data,
				(char*)info->filename->data,
				enccmd.dstbounds[0]-enccmd.srcbounds[1], (char*)enccmd.format->data+enccmd.srcbounds[1],
				(char*)tmpfn1->data,
				(char*)enccmd.format->data+enccmd.dstbounds[1]
			);
		}
		else
		{
			printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "%.*s \"%s\" %.*s \"%s\" %s",
				enccmd.dstbounds[0], (char*)enccmd.format->data,
				(char*)tmpfn1->data,
				enccmd.srcbounds[0]-enccmd.dstbounds[1], (char*)enccmd.format->data+enccmd.dstbounds[1],
				(char*)info->filename->data,
				(char*)enccmd.format->data+enccmd.srcbounds[1]
			);
		}
		++printed;
		int decoffset=printed;
		if(deccmd.srcbounds[0]<deccmd.dstbounds[0])
		{
			printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.srcbounds[0], (char*)deccmd.format->data,
				(char*)tmpfn1->data,
				deccmd.dstbounds[0]-deccmd.srcbounds[1], (char*)deccmd.format->data+deccmd.srcbounds[1],
				(char*)tmpfn2->data,
				(char*)deccmd.format->data+deccmd.dstbounds[1]
			);
		}
		else
		{
			printed+=snprintf(g_buf+printed, sizeof(g_buf)-1-printed, "%.*s \"%s\" %.*s \"%s\" %s",
				deccmd.dstbounds[0], (char*)deccmd.format->data,
				(char*)tmpfn2->data,
				deccmd.srcbounds[0]-deccmd.dstbounds[1], (char*)deccmd.format->data+deccmd.dstbounds[1],
				(char*)tmpfn1->data,
				(char*)deccmd.format->data+deccmd.srcbounds[1]
			);
		}
//#ifdef _DEBUG
//		printf("  %s\n", g_buf+0);//
//		printf("  %s\n", g_buf+decoffset);//
//#endif

		printf("%5d %10lld", k+1, info->usize);
	//	printf("%s\n", g_buf);//
		exec_process(g_buf+0, 0, &currcell->etime, &currcell->emem);
		
		currcell->csize=get_filesize((char*)tmpfn1->data);
		printf(" -> %10lld B  %12.6lf sec %12.6lf MB/s %10lld B", currcell->csize, currcell->etime, info->usize/(currcell->etime*1024*1024), currcell->emem);
		
		exec_process(g_buf+decoffset, 0, &currcell->dtime, &currcell->dmem);
		int kslash=(int)info->filename->count-1;
		while(kslash>=0&&info->filename->data[kslash]!='/'&&info->filename->data[kslash]!='\\')--kslash;
		kslash+=kslash!=0;
		int kdot=(int)info->filename->count-1;
		while(kdot>=0&&info->filename->data[kdot]!='.')--kdot;
		if(kdot<=kslash)
			kdot=(int)info->filename->count;
		printf(" -> %12.6lf sec %12.6lf MB/s %10lld B  %.*s\n", currcell->dtime, info->usize/(currcell->dtime*1024*1024), currcell->dmem, kdot-kslash, (char*)info->filename->data+kslash);//
	//	printf("%s\n", g_buf);//
	//	printf("\n");//

		usize+=info->usize;
		currtest->total.csize+=currcell->csize;
		currtest->total.etime+=currcell->etime;
		if(currtest->total.emem<currcell->emem)
			currtest->total.emem=currcell->emem;
		currtest->total.dtime+=currcell->dtime;
		if(currtest->total.dmem<currcell->dmem)
			currtest->total.dmem=currcell->dmem;
	}
	printf("\n");
	printf("%5d %10lld -> %10lld B  %12.6lf sec %12.6lf MB/s %10lld B -> %12.6lf sec %12.6lf MB/s %10lld B\n",
		(int)uinfo->count,
		usize, currtest->total.csize,
		currtest->total.etime, usize/(currtest->total.etime*1024*1024), currtest->total.emem,
		currtest->total.dtime, usize/(currtest->total.dtime*1024*1024), currtest->total.dmem
	);
#if 1
	{
		ptrdiff_t size;

		size=get_filesize((char*)tmpfn1->data);
		if(size>0)
		{
			int success=DeleteFileA((char*)tmpfn1->data);
			if(!success)
				printf("DeleteFileA GetLastError %d\n", (int)GetLastError());
		}

		size=get_filesize((char*)tmpfn2->data);
		if(size>0)
		{
			int success=DeleteFileA((char*)tmpfn2->data);
			if(!success)
				printf("DeleteFileA GetLastError %d\n", (int)GetLastError());
		}
	}
#endif
	print_timestamp("%Y-%m-%d_%H%M%S\n");

	//6. save
#if 1
	{
		snprintf(g_buf, sizeof(g_buf)-1, "%szzzdata_%s.txt", programpath, datasetname);
		FILE *fdst=fopen(g_buf, "w");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", g_buf);
			return 0;
		}
		fprintf(fdst, "%s\t%s\nfiles:\n", srcpath->data, ext->data);
		for(int k=0;k<(int)uinfo->count;++k)
		{
			UInfo *info=(UInfo*)array_at(&uinfo, k);
			fprintf(fdst, "%10lld\t%s\n", info->usize, info->filename->data);
		}
		fprintf(fdst, "tests:\n");
		for(int k=0;k<(int)testinfo->count;++k)
		{
			TestInfo *info=(TestInfo*)array_at(&testinfo, k);
			struct tm *t=localtime(&info->timestamp);
			fprintf(fdst, "%-20s\t%04d%02d%02d_%02d%02d%02d",
				(char*)info->codecname->data,
				t->tm_year+1900,
				t->tm_mon+1,
				t->tm_mday,
				t->tm_hour,
				t->tm_min,
				t->tm_sec
			);
			write_cell(fdst, &info->total);
			for(int k2=0;k2<(int)info->cells->count;++k2)
			{
				CellInfo *cell=(CellInfo*)array_at(&info->cells, k2);
				write_cell(fdst, cell);
			}
			fprintf(fdst, "\n");
		}
		fclose(fdst);
	}
#endif

	//7. free
	array_free(&enccmd.format);
	array_free(&deccmd.format);
	array_free(&srcpath);
	array_free(&ext);
	array_free(&uinfo);
	array_free(&testinfo);
	array_free(&tmpfn1);
	array_free(&tmpfn2);
	return 0;
}