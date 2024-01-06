//h_console.cpp - allocates and controls the console window
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

#include<io.h>
#include<stdio.h>
#include<fcntl.h>
#ifdef __cplusplus
#include<iostream>
#endif
#include<Windows.h>
#include<stdarg.h>
int consoleactive=0;
#define G_BUF_SIZE	1024
static char g_buf2[G_BUF_SIZE]={0};
int console_log(const char *format, ...)
{
	va_list args;
	int printed=0;
	if(format)
	{
		va_start(args, format);
		printed=vsnprintf(g_buf2, G_BUF_SIZE, format, args);
		va_end(args);
		WriteConsoleA(GetStdHandle(STD_OUTPUT_HANDLE), g_buf2, printed, &printed, 0);
	}
	return printed;
}
void console_start()//https://stackoverflow.com/questions/191842/how-do-i-get-console-output-in-c-with-a-windows-program
{
	if(!consoleactive)
	{
		consoleactive=1;
		int hConHandle;
		size_t lStdHandle;
		//CONSOLE_SCREEN_BUFFER_INFO coninfo;
		FILE *fp;

		// allocate a console for this app
		AllocConsole();

		// set the screen buffer to be big enough to let us scroll text
		//GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);
		//coninfo.dwSize.Y=1000;
		//SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);

		// redirect unbuffered STDOUT to the console
		lStdHandle=(size_t)GetStdHandle(STD_OUTPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stdout=*fp;
		setvbuf(stdout, 0, _IONBF, 0);

		// redirect unbuffered STDIN to the console
		lStdHandle=(size_t)GetStdHandle(STD_INPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "r");
		*stdin=*fp;
		setvbuf(stdin, 0, _IONBF, 0);

		// redirect unbuffered STDERR to the console
		lStdHandle=(size_t)GetStdHandle(STD_ERROR_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stderr=*fp;
		setvbuf(stderr, 0, _IONBF, 0);

#ifdef __cplusplus
		// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
		// point to console as well
		std::ios::sync_with_stdio();
#endif

		console_log("\n\tWARNING: CLOSING THIS WINDOW WILL CLOSE THE PROGRAM\n\n");
	}
}
void console_end()
{
	if(consoleactive)
	{
		FreeConsole();
		consoleactive=0;
	}
}
void console_buffer_size(short x, short y)
{
	COORD coord={x, y};
	SMALL_RECT Rect={0, 0, x-1, y-1};

    HANDLE Handle=GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleScreenBufferSize(Handle, coord);
    SetConsoleWindowInfo(Handle, TRUE, &Rect);
}
//void			console_start(int cx, int cy)
//{
//	int was_active=consoleactive;
//	console_start();
//	console_buffer_size(cx, cy);
//}
//void console_start_good()
//{
//	int was_active=consoleactive;
//	console_start();
//	if(!was_active)
//		console_buffer_size(100, 9000);
//		//console_buffer_size(160, 9000);
//}

int console_scan(char *buf, int len)
{
	unsigned long ret_len=0;
	int success=ReadConsoleA(GetStdHandle(STD_INPUT_HANDLE), buf, len, &ret_len, 0);
	if(!success)
	{
		int error=GetLastError();
		console_log("I/O ERROR %d\n", error);
	}
	return ret_len;
}
int console_scan_int()
{
	int read=console_scan(g_buf2, G_BUF_SIZE);
	return atoi(g_buf2);
}
double console_scan_float()
{
	int read=console_scan(g_buf2, G_BUF_SIZE);
	return atof(g_buf2);
}
void console_pause()
{
	console_log("Enter 0 to continue: ");
	int k=console_scan(g_buf2, G_BUF_SIZE);
	//scanf_s("%d", &k);
	//system("pause");
}