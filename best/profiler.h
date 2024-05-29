/*
profiler.h

1. Choose timer function:

#define PROFILER 0	//for __rdtsc()
#define PROFILER 1	//for time_sec()

2. Define your checkpoints like this:

#define CHECKPOINTLIST\
	CHECKPOINT(ALLOCS)\
	CHECKPOINT(INIT)\
	CHECKPOINT(GETCTX)\
	CHECKPOINT(QUANTIZE)\
	CHECKPOINT(ENTROPYCODER)\
	CHECKPOINT(RESCALE)\
	CHECKPOINT(UPDATE_CDF)\
	CHECKPOINT(TOARRAY)

3. Place checkpoints in your code like this:
	PROF_START();
	...
	PROF(ALLOCS);
	...
	PROF(INIT);
	...
	for(;;)
	{
		PROF(GETCTX);
		...
		PROF(QUANTIZE);
	}
	prof_print();
*/

#ifdef PROFILER
typedef enum ProfilerLabelEnum
{
#define CHECKPOINT(X) PROF_##X,
	CHECKPOINTLIST
#undef  CHECKPOINT
	PROF_COUNT,
} ProfilerLabel;
static const char *prof_labels[]=
{
#define CHECKPOINT(X) #X,
	CHECKPOINTLIST
#undef  CHECKPOINT
};
#if PROFILER==0
#define PROFPRINT "lld"
#define PROFTYPE long long
#define TIMEFUNC __rdtsc
#elif PROFILER==1
#define PROFPRINT "lf"
#define PROFTYPE double
#define TIMEFUNC time_sec
#endif
static volatile PROFTYPE prof_timestamp=0, prof_cycles[PROF_COUNT]={0};
#define PROF_START() memset((void*)prof_cycles, 0, sizeof(prof_cycles)), prof_timestamp=TIMEFUNC()
#define PROF(X) prof_cycles[PROF_##X]+=TIMEFUNC()-prof_timestamp, prof_timestamp=TIMEFUNC()
static void prof_print()
{
	PROFTYPE sum=0;
	int maxlen=0;
	for(int k=0;k<_countof(prof_labels);++k)
	{
		int len=(int)strlen(prof_labels[k]);
		UPDATE_MAX(maxlen, len);
		sum+=prof_cycles[k];
	}
	printf("Profiler:\n");
	for(int k=0;k<_countof(prof_labels);++k)
	{
		double percent=100.*prof_cycles[k]/sum;
		printf("%-*s %16" PROFPRINT " %6.2lf%%  ", maxlen, prof_labels[k], prof_cycles[k], percent);
		for(int k2=0, npoints=(int)percent;k2<npoints;++k2)
			printf("*");
		printf("\n");
	}
}
#else
#define PROF_START()
#define PROF(...)
#define prof_print()
#endif
