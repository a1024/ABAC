#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<errno.h>
#if defined _WIN32 || defined WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
static double time_sec(void)
{
#if defined _WIN32 || defined WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);/*<time.h>*/
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}

typedef struct _Command
{
	const char *name, *fmt;
	int dstfirst;

	intptr_t csize;
	double ebest, eav;
} Command;
static int cmp_commands(const void *left, const void *right)
{
	const Command *p1, *p2;

	p1=(const Command*)left;
	p2=(const Command*)right;
	if(p1->csize==p2->csize)
		return (p1->ebest>p2->ebest)-(p1->ebest<p2->ebest);
	return (p1->csize>p2->csize)-(p1->csize<p2->csize);
}
static Command commands[]=
{/*	name			fmt								dstfirst	csize	ebest	eav*/
	{"EMMA",		"\"emma_c \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"Gralic 1.11 demo",	"\"gralic111d c \"%s\" \"%s\"\"",				1,		0,	0,	0},
	{"X14",			"\"x14 C \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"LEA 6.0 beta",	"\"clea6 \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"LEA 5.0",		"\"clea5 \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"FLIC 2.1 demo",	"\"flic c \"%s\" \"%s\"\"",					1,		0,	0,	0},
	{"C32",			"\"c32 \"%s\" \"%s\" 2\"",					0,		0,	0,	0},
	{"fNBLI 0.4",		"\"fnbli -f \"%s\" -o \"%s\"\"",				0,		0,	0,	0},
	{"Qlic 2",		"\"qlic2 c \"%s\" \"%s\"\"",					1,		0,	0,	0},
	{"L1C2",		"\"l1c-avx2 \"%s\" \"%s\" 2\"",					0,		0,	0,	0},
	{"L1C1",		"\"l1c-avx2 \"%s\" \"%s\" 1\"",					0,		0,	0,	0},
	{"L1C0",		"\"l1c-avx2 \"%s\" \"%s\" 0\"",					0,		0,	0,	0},
	{"Kvick",		"\"kvick c i \"%s\" \"%s\"\"",					1,		0,	0,	0},
	{"Qlic",		"\"qlic c \"%s\" \"%s\"\"",					1,		0,	0,	0},
	{"HALIC072",		"\"HALIC_ENCODE_V.0.7.2_ST_AVX2.exe \"%s\" \"%s\"\"",		0,		0,	0,	0},
	{"Zstd-7",		"\"zstd -f --single-thread -7 \"%s\" -o \"%s\"\"",		0,		0,	0,	0},
	{"Zstd-5",		"\"zstd -f --single-thread -5 \"%s\" -o \"%s\"\"",		0,		0,	0,	0},
	{"Zstd-1",		"\"zstd -f --single-thread -1 \"%s\" -o \"%s\"\"",		0,		0,	0,	0},
	{"C61",			"\"c61 \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"Qic",			"\"qic c \"%s\" \"%s\"\"",					1,		0,	0,	0},
	{"HALIC072FAST",	"\"HALIC_ENCODE_V.0.7.2_ST_FAST_AVX2.exe \"%s\" \"%s\"\"",	0,		0,	0,	0},
	{"ImageZero 0.1",	"\"imagezero c \"%s\" \"%s\"\"",				0,		0,	0,	0},
	{"FOX",			"\"cfox \"%s\" \"%s\"\"",					0,		0,	0,	0},
	{"COPY",		"\"cmd /c copy /b \"%s\" \"%s\"\"",				0,		0,	0,	0},
};
enum
{
	COUNT=30,
	BUFSIZE=8192
};
static char buf[BUFSIZE]={0};
static void print_cmd(Command *cmd, const char *fn, const char *dstfn)
{
	if(cmd->dstfirst)
		snprintf(buf, BUFSIZE, cmd->fmt, dstfn, fn);
	else
		snprintf(buf, BUFSIZE, cmd->fmt, fn, dstfn);
}
int main(int argc, char **argv)
{
	const char *fn, *dstfn;
	struct stat info={0};
	int e, k, k2;
	intptr_t usize;
	double bestenc=0;

	if(argc!=2)
	{
		printf("Usage:  \"%s\"  filename\n", argv[0]);
		return 1;
	}
	fn=argv[1];
	dstfn="20261020Tu_125612_456am.dat";
	e=stat(fn, &info);
	if(e)
	{
		printf("Cannot stat src \"%s\"\n", fn);
		return 1;
	}
	usize=info.st_size;
	for(k=0;k<(int)_countof(commands);++k)
	{
		Command *cmd=commands+k;
		double ttotal=0;
		int failed=0;

		print_cmd(cmd, fn, dstfn);
		printf("%s\n", buf);
		for(k2=0;k2<COUNT;++k2)
		{
			double t;

			t=time_sec();
			e=system(buf);
			t=time_sec()-t;
			if(e)
			{
				printf("Command returned %d\n", errno);
				failed=1;
				break;
			}
			if(!k2||cmd->ebest>t)
				cmd->ebest=t;
			cmd->eav+=t;
			ttotal+=t;
			if(ttotal>90)
				break;
		}
		if(failed)
			cmd->csize=0;
		else
		{
			e=stat(dstfn, &info);
			if(e)
			{
				printf("Cannot stat dst \"%s\"\n", dstfn);
				continue;
			}
			cmd->csize=info.st_size;
			cmd->eav/=k2;
		}
		printf("\n");
#if 0
		Sleep(100);
#endif
	}
	qsort(commands, _countof(commands), sizeof(Command), cmp_commands);
	printf("\"%s\"\n", fn);
	for(k=0;k<(int)_countof(commands);++k)
	{
		Command *cmd=commands+k;

		if(cmd->csize)
		{
			int pareto=0;

			if(!bestenc||bestenc>cmd->ebest)
			{
				bestenc=cmd->ebest;
				pareto=1;
			}
			printf("%12zd bytes  av %12.6f  best %12.6f sec  %10.4f MB/s  %10.4f ms/MB  %c  %s\n"
				, cmd->csize
				, cmd->eav
				, cmd->ebest
				, usize/(cmd->ebest*(1<<20))
				, cmd->ebest*(1000<<20)/usize
				, pareto?'*':' '
				, cmd->name
			);
		}
	}
	printf("%12zd bytes\n", usize);

	snprintf(buf, BUFSIZE, "\"del \"%s\"\"", dstfn);
	printf("%s\n", buf);
	e=system(buf);
	if(e)
		printf("Command returned %d\n", errno);
#if 0
	system("pause");
#endif
	return 0;
}
