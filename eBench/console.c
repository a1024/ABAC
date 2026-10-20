#include<stdint.h>
#include<io.h>
#include<stdio.h>
#include<fcntl.h>
#ifdef __cplusplus
#include<iostream>
#endif
#include<Windows.h>
#include<stdarg.h>


void   console_start(void);
void   console_end(void);
void   console_buffer_size(short x, short y);
int    console_log(const char *format, ...);
void   console_pause(void);
int    console_scan(char *buf, int len);
int    console_scan_int(void);
double console_scan_float(void);


int consoleactive=0;
#define G_BUF_SIZE 1024
static char g_buf2[G_BUF_SIZE]={0};

int console_log(const char *format, ...)
{
	unsigned long printed;

	printed=0;
	if(format)
	{
		va_list args;

		va_start(args, format);
		printed=vsnprintf(g_buf2, G_BUF_SIZE, format, args);
		va_end(args);
		WriteConsoleA(GetStdHandle(STD_OUTPUT_HANDLE), g_buf2, printed, &printed, 0);
	}
	return printed;
}
/*
https://stackoverflow.com/questions/191842/how-do-i-get-console-output-in-c-with-a-windows-program
*/
void console_start(void)
{
	if(!consoleactive)
	{
		int hConHandle;
		size_t lStdHandle;
		FILE *fp;
#if 0
		CONSOLE_SCREEN_BUFFER_INFO coninfo;
#endif

		consoleactive=1;

		/*allocate a console for current process*/
		AllocConsole();

		/*set the screen buffer to be big enough to let us scroll text*/
#if 0
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);
		coninfo.dwSize.Y=1000;
		SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);
#endif

		/*redirect unbuffered STDOUT to the console*/
		lStdHandle=(size_t)GetStdHandle(STD_OUTPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stdout=*fp;
		setvbuf(stdout, 0, _IONBF, 0);

		/*redirect unbuffered STDIN to the console*/
		lStdHandle=(size_t)GetStdHandle(STD_INPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "r");
		*stdin=*fp;
		setvbuf(stdin, 0, _IONBF, 0);

		/*redirect unbuffered STDERR to the console*/
		lStdHandle=(size_t)GetStdHandle(STD_ERROR_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stderr=*fp;
		setvbuf(stderr, 0, _IONBF, 0);

#ifdef __cplusplus
		/*
		make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
		point to console as well
		*/
		std::ios::sync_with_stdio();
#endif

		console_log("\n\tWARNING: CLOSING THIS WINDOW WILL CLOSE THE PROGRAM\n\n");
	}
}
void console_end(void)
{
	if(consoleactive)
	{
		FreeConsole();
		consoleactive=0;
	}
}
void console_buffer_size(short x, short y)
{
#if 0
	COORD coord={x, y};
	SMALL_RECT rect={0, 0, (short)(x-1), (short)(y-1)};
#endif
	COORD coord={0};
	SMALL_RECT rect={0};
	HANDLE Handle=0;

	coord.X=x;
	coord.Y=y;
	rect.Right=(int16_t)(x-1);
	rect.Bottom=(int16_t)(y-1);
	Handle=GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleScreenBufferSize(Handle, coord);
	SetConsoleWindowInfo(Handle, TRUE, &rect);
}
#if 0
void console_start(int cx, int cy)
{
	int was_active=consoleactive;
	console_start();
	console_buffer_size(cx, cy);
}
void console_start_good()
{
	int was_active=consoleactive;
	console_start();
	if(!was_active)
		console_buffer_size(100, 9000);
#if 0
		console_buffer_size(160, 9000);
#endif
}
#endif

int console_scan(char *buf, int len)
{
	DWORD ret_len;
	int success;

	ret_len=0;
	success=ReadConsoleA(GetStdHandle(STD_INPUT_HANDLE), buf, len, &ret_len, 0);
	if(!success)
	{
		int error=GetLastError();
		console_log("I/O ERROR %d\n", error);
	}
	return ret_len;
}
int console_scan_int(void)
{
	console_scan(g_buf2, G_BUF_SIZE);
	return atoi(g_buf2);
}
double console_scan_float(void)
{
	console_scan(g_buf2, G_BUF_SIZE);
	return atof(g_buf2);
}
void console_pause(void)
{
	console_log("Enter 0 to continue: ");
	console_scan(g_buf2, G_BUF_SIZE);
#if 0
	scanf_s("%d", &k);
	system("pause");
#endif
}
