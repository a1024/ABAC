#ifdef _MSC_VER
#	define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<wchar.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#include<Shellapi.h>
#include<sys/stat.h>

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

/*
 * command:    Full command line, including executable and arguments.
 * result:     Output buffer, in wchar_t characters.
 * result_cap: Number of wchar_t elements in result.
 *
 * Returns the exit code
 *
 * On failure, result may contain partial output.
 */
static int RunCommandW(const wchar_t *command, wchar_t *result, size_t result_cap, double *ret_elapsed, int *ret_didrun)
{
	SECURITY_ATTRIBUTES sa={0};
	STARTUPINFOW si={0};
	PROCESS_INFORMATION pi={0};
	HANDLE read_pipe=NULL, write_pipe=NULL;
	wchar_t *cmdline=NULL;
	char *bytes=NULL;
	size_t bytes_cap=0, bytes_len=0;
	DWORD exit_code=1;
	int didrun=0;
	double t=0;

	if(!command || !result || result_cap == 0)
		return 0;

	result[0]=L'\0';

	sa.nLength=sizeof(sa);
	sa.lpSecurityDescriptor=NULL;
	sa.bInheritHandle=TRUE;

	if(!CreatePipe(&read_pipe, &write_pipe, &sa, 0))
		goto cleanup;

	/* The parent must not inherit the read end. */
	if(!SetHandleInformation(read_pipe, HANDLE_FLAG_INHERIT, 0))
		goto cleanup;

	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));

	si.cb=sizeof(si);
	si.dwFlags=STARTF_USESTDHANDLES;
	si.hStdOutput=write_pipe;
	si.hStdError=write_pipe;
	si.hStdInput=GetStdHandle(STD_INPUT_HANDLE);

	cmdline=_wcsdup(command);
	if(!cmdline)
		goto cleanup;

	t=time_sec();
	if(!CreateProcessW(NULL, cmdline, NULL, NULL, TRUE /*Inherit pipe handles*/, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
		goto cleanup;

	/* Parent must close its write end or reading may never reach EOF. */
	CloseHandle(write_pipe);
	write_pipe=NULL;

	/*
	* Drain output while the child runs. Do not wait first:
	* a full pipe can block the child and deadlock the parent.
	*/
	for(;;)
	{
		char chunk[4096];
		DWORD nread=0;

		if(!ReadFile(read_pipe, chunk, sizeof(chunk), &nread, NULL))
		{
			if(GetLastError()==ERROR_BROKEN_PIPE)
				break;
			goto cleanup;
		}
		if(!nread)
			break;
		if(bytes_len+nread>bytes_cap)
		{
			size_t new_cap=bytes_cap ? bytes_cap * 2 : 8192;
			while (new_cap < bytes_len+nread)
				new_cap *= 2;

			char *tmp=(char *)realloc(bytes, new_cap);
			if(!tmp)
				goto cleanup;

			bytes=tmp;
			bytes_cap=new_cap;
		}
		memcpy(bytes+bytes_len, chunk, nread);
		bytes_len+=nread;
	}

	WaitForSingleObject(pi.hProcess, INFINITE);
	t=time_sec()-t;

	if(!GetExitCodeProcess(pi.hProcess, &exit_code))
		goto cleanup;

	/* Convert captured UTF-8 output to wide characters. */
	if(bytes_len>0)
	{
		int needed=MultiByteToWideChar(CP_UTF8, 0, bytes, (int)bytes_len, NULL, 0);

		if(needed<=0 || (size_t)needed>=result_cap)
			goto cleanup;

		if(!MultiByteToWideChar(CP_UTF8, 0, bytes, (int)bytes_len, result, needed))
			goto cleanup;

		result[needed]=L'\0';
	}
	didrun=1;
cleanup:
	if(write_pipe) CloseHandle(write_pipe);
	if(read_pipe)  CloseHandle(read_pipe);

	if(pi.hThread)  CloseHandle(pi.hThread);
	if(pi.hProcess) CloseHandle(pi.hProcess);

	free(cmdline);
	free(bytes);
	
	if(ret_didrun)*ret_didrun=didrun;
	if(ret_elapsed)*ret_elapsed=didrun?t:0;

	return (int)exit_code;
}
static int GetTempFn(wchar_t *filename, DWORD capacity)
{
	wchar_t temp_path[MAX_PATH*2];

	DWORD len=GetTempPathW(_countof(temp_path), temp_path);
	if(!len||len>=_countof(temp_path))
		return 0;
	if(!GetTempFileNameW(temp_path, L"rdc", 0, filename))
		return 0;
	(void)capacity;
	return 1;
}
typedef struct _Result
{
	size_t csize;
	double e, d, lossyscores[2];
	char preset[88];
} Result;
int main(int argc, char **argv)
{
	enum
	{
		CODEC_JPG,
		CODEC_J2K,
		CODEC_JXL,
		CODEC_C47,
		CODEC_C60,
		CODEC_L1C,
		CODEC_COUNT,
	};
	enum
	{
		BUFSIZE=1<<16
	};
	const int jpg_scale[]=/*lower is better*/
	{
		2, 3, 4, 5, 7, 9, 11, 15, 19, 23, 27,
	};
	const int j2k_ratio[]=/*lower is better*/
	{
		2, 4, 8, 16, 24, 32, 48, 64, 96, 128,
	};
	const int jxl_quality[]=/*higher is better*/
	{
		100,
		95, 90, 85, 80, 75, 70, 65, 60, 50, 40, 30, 20, 10,
	//	0,/*lossless*/
	//	0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 10
	};
	const int c47_quality[]=/*?*/
	{
		75,
	};
	const int c60_dist[]=/*lower is better*/
	{
		1, 2, 3, 4, 5, 6, 7,
	};
	const int l1c_dist[]=/*lower is better*/
	{
		1, 3, 4, 5, 7, 9, 11, 15,
	};
	const int counts[]=
	{
		_countof(jpg_scale),
		_countof(j2k_ratio),
		_countof(jxl_quality),
		_countof(c47_quality),
		_countof(c60_dist),
		_countof(l1c_dist),
	};
	static wchar_t buf[BUFSIZE]={0}, fn1[MAX_PATH*2]={0}, fn2[MAX_PATH*2]={0};
	static Result results[128]={0};
	struct _stat64 info={0};
	const wchar_t *srcfn=0;
	wchar_t *cmdline, **wargv, *ext;
	size_t srcsize;
	int wargc, error, usessimulacra, len, kc, kp, nresults, success, totalpresets;
	
	cmdline=GetCommandLineW();
	wargv=CommandLineToArgvW(cmdline, &wargc);
#if defined _MSC_VER
	srcfn=
		L"C:/dataset-GDCC2020-ppm/photo-72.ppm"
	;
#elif defined __GNUC__
	(void)argc;
	(void)argv;
	if(wargc!=2)
	{
		wprintf(L"Usage:  \"%ls\"  image.ppm\n", wargv[0]);
		return 1;
	}
	srcfn=wargv[1];
#endif
	error=_wstat64(srcfn, &info);
	if(error)
	{
		wprintf(L"Cannot stat \"%ls\"\n", srcfn);
		return 1;
	}
	srcsize=info.st_size;
	if(!GetTempFn(fn1, _countof(fn1))||!GetTempFn(fn2, _countof(fn2)))
	{
		CRASH("Error creating temp file(s)");
		return 1;
	}
	ext=wcsrchr(srcfn, L'.');
	len=(int)wcslen(srcfn);
	if(!ext||ext+1>=srcfn+len)
	{
		wprintf(L"No extension in \"%ls\"\n", srcfn);
		return 1;
	}
	if(_wcsicmp(ext, L".PPM")&&_wcsicmp(ext, L".PGM")&&_wcsicmp(ext, L".PNM"))
	{
		wprintf(L"Uncompressed image required for PSNR.\n\t%ls\n", srcfn);
		return 1;
	}
	wcsncat(fn2, ext, _countof(fn2)-wcslen(fn2));
	wprintf(L"%12zu  \"%ls\"\n", srcsize, srcfn);

	RunCommandW(L"ssimulacra2", buf, BUFSIZE, 0, &usessimulacra);
	if(!usessimulacra)
		printf("ssimulacra2 not found. Defaulting to PSNR.\n");
	nresults=0;
	totalpresets=_countof(counts);
	for(kp=0;kp<(int)_countof(counts);++kp)
		totalpresets+=counts[kp];
	for(kc=0;kc<CODEC_COUNT;++kc)
	{
		Result *res;
		wchar_t *p, *ext1;
		int q, printed, i;
		
		ext1=0;
		for(kp=counts[kc]-1;kp>=0;--kp)
		{
			if(nresults>=(int)_countof(results))
			{
				printf("OUT OF MEMORY\n");
				break;
			}
			res=results+nresults++;
			switch(kc)
			{
			case CODEC_JPG:
				ext1=L"jpg";
				q=jpg_scale[kp];
				snprintf(res->preset, _countof(res->preset)-1, "JPG_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"ffmpeg -hide_banner -loglevel error -i \"%ls\" -frames:v 1 -c:v mjpeg -q:v %d -y \"%ls.%ls\"", srcfn, q, fn1, ext1);
				break;
			case CODEC_J2K:
				ext1=L"j2k";
				q=j2k_ratio[kp];
				snprintf(res->preset, _countof(res->preset)-1, "J2K_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"opj_compress -quiet -i \"%ls\" -o \"%ls.%ls\" -I -r %d", srcfn, fn1, ext1, q);
				break;
			case CODEC_JXL:
				ext1=L"jxl";
				q=jxl_quality[kp];
				snprintf(res->preset, _countof(res->preset)-1, "JXL_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"cjxl --quiet -q %d \"%ls\" \"%ls.%ls\"", q, srcfn, fn1, ext1);
				break;
			case CODEC_C47:
				ext1=L"c47";
				q=c47_quality[kp];
				snprintf(res->preset, _countof(res->preset)-1, "C47_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"c47 \"%ls\" \"%ls.%ls\" %d", srcfn, fn1, ext1, q);
				break;
			case CODEC_C60:
				ext1=L"c60";
				q=c60_dist[kp];
				snprintf(res->preset, _countof(res->preset)-1, "C60_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"c60 \"%ls\" \"%ls.%ls\" %d", srcfn, fn1, ext1, q);
				break;
			case CODEC_L1C:
				ext1=L"l1c";
				q=l1c_dist[kp];
				snprintf(res->preset, _countof(res->preset)-1, "L1C_%03d", q);
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"l1c-avx2 \"%ls\" \"%ls.%ls\" 2 %d", srcfn, fn1, ext1, q);
				break;
			}
			printf("[%3d/%3d] %s\r", nresults, totalpresets, res->preset);
			if(printed>=BUFSIZE-8192)
			{
				CRASH("Buffer OOB");
				return 1;
			}
			buf[printed++]=L'\0';
			RunCommandW(buf, buf+printed, (size_t)BUFSIZE-1-printed, &res->e, &success);
			if(!success)
			{
				--nresults;
				wprintf(L"\n%ls\n\tCommand failed\n", buf);
				break;
			}

			printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"%ls.%ls", fn1, ext1);
			error=_wstat64(buf, &info);
			if(error)
			{
				wprintf(L"Cannot stat \"%ls\"\n", buf);
				wprintf(L"Enc:\n%ls\n", buf+printed);
				CRASH("Error");
				return 1;
			}
			res->csize=info.st_size;
			
			switch(kc)
			{
			case CODEC_JPG:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"ffmpeg -hide_banner -loglevel error -i \"%ls.%ls\" -frames:v 1 -y \"%ls\"", fn1, ext1, fn2);
				break;
			case CODEC_J2K:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"opj_decompress -quiet -i \"%ls.%ls\" -o \"%ls\"", fn1, ext1, fn2);
				break;
			case CODEC_JXL:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"djxl --quiet \"%ls.%ls\" \"%ls\"", fn1, ext1, fn2);
				break;
			case CODEC_C47:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"c47 \"%ls.%ls\" \"%ls\"", fn1, ext1, fn2);
				break;
			case CODEC_C60:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"c60 \"%ls.%ls\" \"%ls\"", fn1, ext1, fn2);
				break;
			case CODEC_L1C:
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"l1c-avx2 \"%ls.%ls\" \"%ls\"", fn1, ext1, fn2);
				break;
			}
			if(printed>=BUFSIZE-1)
			{
				CRASH("Buffer OOB");
				return 1;
			}
			buf[printed++]=L'\0';
			RunCommandW(buf, buf+printed, (size_t)BUFSIZE-1-printed, &res->d, 0);

			error=_wstat64(fn2, &info);
			if(error)
			{
				wprintf(L"Cannot stat \"%ls\"\n", fn2);
				wprintf(L"Dec:\n%ls\n", buf+printed);
				CRASH("Error");
				return 1;
			}
			for(i=0;i<=usessimulacra;++i)
			{
				printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"%ls \"%ls\" \"%ls\""
					, i?L"ssimulacra2":L"psnr"
					, srcfn
					, fn2
				);
				if(printed>=BUFSIZE-1)
				{
					CRASH("Buffer OOB");
					return 1;
				}
				buf[printed++]=L'\0';
				RunCommandW(buf, buf+printed, (size_t)BUFSIZE-1-printed, 0, 0);
				p=0;
				res->lossyscores[i]=wcstod(buf+printed, &p);
			}
		}

		printed=_snwprintf(buf, (size_t)BUFSIZE-1, L"%ls.%ls", fn1, ext1);
		buf[printed]=L'\0';
		success=DeleteFileW(buf);
		if(!success)
		{
			CRASH("DeleteFileW %d", GetLastError());
			return 1;
		}
		++nresults;
	}
	success=DeleteFileW(fn2);
	if(!success)
	{
		CRASH("DeleteFileW %d", GetLastError());
		return 1;
	}

	printf("CSIZE,  BPS,  PSNR%s,  ENC  DEC (sec),  CODEC\n", usessimulacra?",  SSIMULACRA2":"");
	for(kp=0;kp<nresults;++kp)
	{
		Result *res;

		res=results+kp;
		if(!res->csize)
		{
			printf("\n");
			continue;
		}
		printf("%12zu  %10.6f  %10.6f"
			, res->csize
			, res->csize*8./srcsize
			, res->lossyscores[0]
		);
		if(usessimulacra)
			printf("  %10.6f"
				, res->lossyscores[1]
			);
		printf("  %12.6f  %12.6f  %s\n"
			, res->e
			, res->d
			, res->preset
		);
	}
	return 0;
}
