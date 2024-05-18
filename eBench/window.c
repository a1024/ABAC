#include<stdio.h>
#include<stdlib.h>
#include<Windowsx.h>//for GET_X_LPARAM

#define STRICT_TYPED_ITEMIDS
#include<shlobj.h>
#include<objbase.h>	// For COM headers
#include<shobjidl.h>	// for IFileDialogEvents and IFileDialogControlEvents
#include<shlwapi.h>
#include<knownfolders.h>// for KnownFolder APIs/datatypes/function headers
#include<propvarutil.h>	// for PROPVAR-related functions
#include<propkey.h>	// for the Property key APIs/datatypes
#include<propidl.h>	// for the Property System APIs
#include<strsafe.h>	// for StringCchPrintfW
#include<shtypes.h>	// for COMDLG_FILTERSPEC
#include<combaseapi.h>	// for IID_PPV_ARGS

#include"ebench.h"
static const char file[]=__FILE__;

int wndw=0, wndh=0,
	mx=0, my=0, mouse_bypass=0;
HWND ghWnd=0;
HDC ghDC=0;
HGLRC hRC=0;
char keyboard[256]={0}, timer=0;
int g_repaint=0;
RECT R={0};
wchar_t g_wbuf[G_BUF_SIZE]={0};
ArrayHandle exedir=0;

#define UTF8TOWCHAR(U8, LEN, RET_U16, BUF_SIZE, RET_LEN)    RET_LEN=MultiByteToWideChar(CP_UTF8, 0, U8, LEN, RET_U16, BUF_SIZE)
#define WCHARTOUTF8(WSTR, LEN, RET_UTF8, BUF_SIZE, RET_LEN) RET_LEN=WideCharToMultiByte(CP_UTF8, 0, WSTR, LEN, RET_UTF8, BUF_SIZE, 0, 0)

const char* wm2str(int message)
{
	const char *a="???";
	switch(message)
	{//message case
#define		MC(x)	case x:a=#x;break;
	MC(WM_NULL)
	MC(WM_CREATE)
	MC(WM_DESTROY)
	MC(WM_MOVE)
	MC(WM_SIZE)

	MC(WM_ACTIVATE)
	MC(WM_SETFOCUS)
	MC(WM_KILLFOCUS)
	MC(WM_ENABLE)
	MC(WM_SETREDRAW)
	MC(WM_SETTEXT)
	MC(WM_GETTEXT)
	MC(WM_GETTEXTLENGTH)
	MC(WM_PAINT)
	MC(WM_CLOSE)

	MC(WM_QUERYENDSESSION)
	MC(WM_QUERYOPEN)
	MC(WM_ENDSESSION)

	MC(WM_QUIT)
	MC(WM_ERASEBKGND)
	MC(WM_SYSCOLORCHANGE)
	MC(WM_SHOWWINDOW)
	MC(WM_WININICHANGE)
//	MC(WM_SETTINGCHANGE)//==WM_WININICHANGE

	MC(WM_DEVMODECHANGE)
	MC(WM_ACTIVATEAPP)
	MC(WM_FONTCHANGE)
	MC(WM_TIMECHANGE)
	MC(WM_CANCELMODE)
	MC(WM_SETCURSOR)
	MC(WM_MOUSEACTIVATE)
	MC(WM_CHILDACTIVATE)
	MC(WM_QUEUESYNC)

	MC(WM_GETMINMAXINFO)

	MC(WM_PAINTICON)
	MC(WM_ICONERASEBKGND)
	MC(WM_NEXTDLGCTL)
	MC(WM_SPOOLERSTATUS)
	MC(WM_DRAWITEM)
	MC(WM_MEASUREITEM)
	MC(WM_DELETEITEM)
	MC(WM_VKEYTOITEM)
	MC(WM_CHARTOITEM)
	MC(WM_SETFONT)
	MC(WM_GETFONT)
	MC(WM_SETHOTKEY)
	MC(WM_GETHOTKEY)
	MC(WM_QUERYDRAGICON)
	MC(WM_COMPAREITEM)

	MC(WM_GETOBJECT)

	MC(WM_COMPACTING)
	MC(WM_COMMNOTIFY)
	MC(WM_WINDOWPOSCHANGING)
	MC(WM_WINDOWPOSCHANGED)

	MC(WM_POWER)

	MC(WM_COPYDATA)
	MC(WM_CANCELJOURNAL)

	MC(WM_NOTIFY)
	MC(WM_INPUTLANGCHANGEREQUEST)
	MC(WM_INPUTLANGCHANGE)
	MC(WM_TCARD)
	MC(WM_HELP)
	MC(WM_USERCHANGED)
	MC(WM_NOTIFYFORMAT)

	MC(WM_CONTEXTMENU)
	MC(WM_STYLECHANGING)
	MC(WM_STYLECHANGED)
	MC(WM_DISPLAYCHANGE)
	MC(WM_GETICON)
	MC(WM_SETICON)

	MC(WM_NCCREATE)
	MC(WM_NCDESTROY)
	MC(WM_NCCALCSIZE)
	MC(WM_NCHITTEST)
	MC(WM_NCPAINT)
	MC(WM_NCACTIVATE)
	MC(WM_GETDLGCODE)

	MC(WM_SYNCPAINT)

	MC(WM_NCMOUSEMOVE)
	MC(WM_NCLBUTTONDOWN)
	MC(WM_NCLBUTTONUP)
	MC(WM_NCLBUTTONDBLCLK)
	MC(WM_NCRBUTTONDOWN)
	MC(WM_NCRBUTTONUP)
	MC(WM_NCRBUTTONDBLCLK)
	MC(WM_NCMBUTTONDOWN)
	MC(WM_NCMBUTTONUP)
	MC(WM_NCMBUTTONDBLCLK)

	MC(WM_NCXBUTTONDOWN  )
	MC(WM_NCXBUTTONUP    )
	MC(WM_NCXBUTTONDBLCLK)

	MC(WM_INPUT_DEVICE_CHANGE)

	MC(WM_INPUT)

//	MC(WM_KEYFIRST   )//==WM_KEYDOWN
	MC(WM_KEYDOWN    )
	MC(WM_KEYUP      )
	MC(WM_CHAR       )
	MC(WM_DEADCHAR   )
	MC(WM_SYSKEYDOWN )
	MC(WM_SYSKEYUP   )
	MC(WM_SYSCHAR    )
	MC(WM_SYSDEADCHAR)

	MC(WM_UNICHAR)
//	MC(WM_KEYLAST)		//==WM_UNICHAR
	MC(UNICODE_NOCHAR)	//0xFFFF

	MC(WM_IME_STARTCOMPOSITION)
	MC(WM_IME_ENDCOMPOSITION)
	MC(WM_IME_COMPOSITION)
//	MC(WM_IME_KEYLAST)	//==WM_IME_KEYLAST

	MC(WM_INITDIALOG   )
	MC(WM_COMMAND      )
	MC(WM_SYSCOMMAND   )
	MC(WM_TIMER        )
	MC(WM_HSCROLL      )
	MC(WM_VSCROLL      )
	MC(WM_INITMENU     )
	MC(WM_INITMENUPOPUP)

	MC(WM_GESTURE      )
	MC(WM_GESTURENOTIFY)

	MC(WM_MENUSELECT)
	MC(WM_MENUCHAR  )
	MC(WM_ENTERIDLE )

	MC(WM_MENURBUTTONUP  )
	MC(WM_MENUDRAG       )
	MC(WM_MENUGETOBJECT  )
	MC(WM_UNINITMENUPOPUP)
	MC(WM_MENUCOMMAND    )

	MC(WM_CHANGEUISTATE)
	MC(WM_UPDATEUISTATE)
	MC(WM_QUERYUISTATE )

	MC(WM_CTLCOLORMSGBOX   )
	MC(WM_CTLCOLOREDIT     )
	MC(WM_CTLCOLORLISTBOX  )
	MC(WM_CTLCOLORBTN      )
	MC(WM_CTLCOLORDLG      )
	MC(WM_CTLCOLORSCROLLBAR)
	MC(WM_CTLCOLORSTATIC   )
	MC(MN_GETHMENU         )

//	MC(WM_MOUSEFIRST   )
	MC(WM_MOUSEMOVE    )
	MC(WM_LBUTTONDOWN  )
	MC(WM_LBUTTONUP    )
	MC(WM_LBUTTONDBLCLK)
	MC(WM_RBUTTONDOWN  )
	MC(WM_RBUTTONUP    )
	MC(WM_RBUTTONDBLCLK)
	MC(WM_MBUTTONDOWN  )
	MC(WM_MBUTTONUP    )
	MC(WM_MBUTTONDBLCLK)

	MC(WM_MOUSEWHEEL)

	MC(WM_XBUTTONDOWN  )
	MC(WM_XBUTTONUP    )
	MC(WM_XBUTTONDBLCLK)

//	MC(WM_MOUSELAST)	//==WM_MOUSEWHEEL

	MC(WM_PARENTNOTIFY )
	MC(WM_ENTERMENULOOP)
	MC(WM_EXITMENULOOP )

	MC(WM_NEXTMENU      )
	MC(WM_SIZING        )
	MC(WM_CAPTURECHANGED)
	MC(WM_MOVING        )

	MC(WM_POWERBROADCAST)

	MC(WM_DEVICECHANGE)

	MC(WM_MDICREATE     )
	MC(WM_MDIDESTROY    )
	MC(WM_MDIACTIVATE   )
	MC(WM_MDIRESTORE    )
	MC(WM_MDINEXT       )
	MC(WM_MDIMAXIMIZE   )
	MC(WM_MDITILE       )
	MC(WM_MDICASCADE    )
	MC(WM_MDIICONARRANGE)
	MC(WM_MDIGETACTIVE  )

	MC(WM_MDISETMENU    )
	MC(WM_ENTERSIZEMOVE )
	MC(WM_EXITSIZEMOVE  )
	MC(WM_DROPFILES     )
	MC(WM_MDIREFRESHMENU)

//	MC(WM_POINTERDEVICECHANGE    )
//	MC(WM_POINTERDEVICEINRANGE   )
//	MC(WM_POINTERDEVICEOUTOFRANGE)

	MC(WM_TOUCH)

//	MC(WM_NCPOINTERUPDATE      )
//	MC(WM_NCPOINTERDOWN        )
//	MC(WM_NCPOINTERUP          )
//	MC(WM_POINTERUPDATE        )
//	MC(WM_POINTERDOWN          )
//	MC(WM_POINTERUP            )
//	MC(WM_POINTERENTER         )
//	MC(WM_POINTERLEAVE         )
//	MC(WM_POINTERACTIVATE      )
//	MC(WM_POINTERCAPTURECHANGED)
//	MC(WM_TOUCHHITTESTING      )
//	MC(WM_POINTERWHEEL         )
//	MC(WM_POINTERHWHEEL        )
//	MC(DM_POINTERHITTEST       )

	MC(WM_IME_SETCONTEXT     )
	MC(WM_IME_NOTIFY         )
	MC(WM_IME_CONTROL        )
	MC(WM_IME_COMPOSITIONFULL)
	MC(WM_IME_SELECT         )
	MC(WM_IME_CHAR           )

	MC(WM_IME_REQUEST)

	MC(WM_IME_KEYDOWN)
	MC(WM_IME_KEYUP  )

	MC(WM_MOUSEHOVER)
	MC(WM_MOUSELEAVE)

	MC(WM_NCMOUSEHOVER)
	MC(WM_NCMOUSELEAVE)

	MC(WM_WTSSESSION_CHANGE)

	MC(WM_TABLET_FIRST)
	MC(WM_TABLET_LAST )

//	MC(WM_DPICHANGED)

	MC(WM_CUT              )
	MC(WM_COPY             )
	MC(WM_PASTE            )
	MC(WM_CLEAR            )
	MC(WM_UNDO             )
	MC(WM_RENDERFORMAT     )
	MC(WM_RENDERALLFORMATS )
	MC(WM_DESTROYCLIPBOARD )
	MC(WM_DRAWCLIPBOARD    )
	MC(WM_PAINTCLIPBOARD   )
	MC(WM_VSCROLLCLIPBOARD )
	MC(WM_SIZECLIPBOARD    )
	MC(WM_ASKCBFORMATNAME  )
	MC(WM_CHANGECBCHAIN    )
	MC(WM_HSCROLLCLIPBOARD )
	MC(WM_QUERYNEWPALETTE  )
	MC(WM_PALETTEISCHANGING)
	MC(WM_PALETTECHANGED   )
	MC(WM_HOTKEY           )

	MC(WM_PRINT      )
	MC(WM_PRINTCLIENT)

	MC(WM_APPCOMMAND)

	MC(WM_THEMECHANGED)

	MC(WM_CLIPBOARDUPDATE)

	MC(WM_DWMCOMPOSITIONCHANGED      )
	MC(WM_DWMNCRENDERINGCHANGED      )
	MC(WM_DWMCOLORIZATIONCOLORCHANGED)
	MC(WM_DWMWINDOWMAXIMIZEDCHANGE   )

	MC(WM_DWMSENDICONICTHUMBNAIL        )
	MC(WM_DWMSENDICONICLIVEPREVIEWBITMAP)

	MC(WM_GETTITLEBARINFOEX)

	MC(WM_HANDHELDFIRST)
	MC(WM_HANDHELDLAST )

	MC(WM_AFXFIRST)
	MC(WM_AFXLAST )

	MC(WM_PENWINFIRST)
	MC(WM_PENWINLAST )

	MC(WM_APP)

	MC(WM_USER)
#undef		MC
	}
	return a;
}
int sys_check(const char *fn, int line, const char *info)
{
	int e2=GetLastError();
	if(e2)
	{
		char *messageBuffer=0;
		/*size_t size=*/FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER|FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS, NULL, e2, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);
		log_error(fn, line, 1, "%s%sGetLastError() returned %d: %s", info?info:"", info?"\n":"", e2, messageBuffer);
		LocalFree(messageBuffer);
	}
	return 0;
}

static int format_utf8_message(const char *title, const char *format, char *args)//returns idx of title in g_wbuf
{
	int len=vsnprintf(g_buf, G_BUF_SIZE, format, args);
	len=MultiByteToWideChar(CP_UTF8, 0, g_buf, len, g_wbuf, G_BUF_SIZE);	SYS_ASSERT(len);
	g_wbuf[len]='\0';
	++len;
	int len2=MultiByteToWideChar(CP_UTF8, 0, title, (int)strlen(title), g_wbuf+len, G_BUF_SIZE-len);	SYS_ASSERT(len2);
	g_wbuf[len+len2]='\0';
	return len;
}
int messagebox(MessageBoxType type, const char *title, const char *format, ...)//returns index of pressed button
{
	int len=format_utf8_message(title, format, (char*)(&format+1));
	int wintypes[]={MB_OK, MB_OKCANCEL, MB_YESNOCANCEL};
	int result=MessageBoxW(ghWnd, g_wbuf, g_wbuf+len, wintypes[type]);
	switch(type)
	{
	case MBOX_OK:result=0;
	case MBOX_OKCANCEL:
		switch(result)
		{
		case IDOK:
			result=0;
			break;
		case IDCANCEL:
		default:
			result=1;
			break;
		}
		break;
	case MBOX_YESNOCANCEL:
		switch(result)
		{
		case IDYES:		result=0;	break;
		case IDNO:		result=1;	break;
		default:
		case IDCANCEL:	result=2;	break;
		}
		break;
	}
	return result;
}

//static void free_pchar(void *data)
//{
//	char **str=(char**)data;
//	free(*str);
//	*str=0;
//}
ArrayHandle dialog_open_folder(void)
{
	ArrayHandle arr=0;
	IFileOpenDialog *pFileOpenDialog=0;
	HRESULT hr=OleInitialize(0);
	if(hr!=S_OK)
	{
		OleUninitialize();
		return 0;
	}
	IID fileOpenDialogIID={0xD57C7288, 0xD4AD, 0x4768, {0xBE, 0x02, 0x9D, 0x96, 0x95, 0x32, 0xD9, 0x60}};//IFileOpenDialog
	hr=CoCreateInstance((const IID*)&CLSID_FileOpenDialog, 0, CLSCTX_INPROC_SERVER, &fileOpenDialogIID, (void*)&pFileOpenDialog);
	if(SUCCEEDED(hr))
	{
		int success=0, len=0;
#ifdef __cplusplus
#define	CALL_METHOD(OBJ, METHOD, ...)	(OBJ)->METHOD(__VA_ARGS__)
#else
#define	CALL_METHOD(OBJ, METHOD, ...)	(OBJ)->lpVtbl->METHOD(OBJ, ##__VA_ARGS__)
#endif
		hr=CALL_METHOD(pFileOpenDialog, SetOptions, FOS_PICKFOLDERS|FOS_FORCEFILESYSTEM);
		hr=CALL_METHOD(pFileOpenDialog, Show, 0);
		success=SUCCEEDED(hr);
		if(success)
		{
			IShellItem *pShellItem=0;
			wchar_t *fullpath=0;
			hr=CALL_METHOD(pFileOpenDialog, GetResult, &pShellItem);
			CALL_METHOD(pShellItem, GetDisplayName, SIGDN_FILESYSPATH, &fullpath);
			WCHARTOUTF8(fullpath, (int)wcslen(fullpath), g_buf, G_BUF_SIZE, len);

			g_buf[len]=0;
			arr=filter_path(g_buf, len, 1);
			//STR_COPY(arr, g_buf, len);

			CoTaskMemFree(fullpath);
		}
		CALL_METHOD(pFileOpenDialog, Release);
#undef	CALL_METHOD
	}
	OleUninitialize();
	return arr;
}
static ArrayHandle	prep_filters(Filter *filters, int nfilters)
{
	ArrayHandle winfilts=0;
	WSTR_ALLOC(winfilts, 0);
	for(int k=0;k<nfilters;++k)
	{
		int len=0;

		UTF8TOWCHAR(filters[k].comment, (int)strlen(filters[k].comment)+1, g_wbuf, G_BUF_SIZE, len);
		if(!len)
			break;
		STR_APPEND(winfilts, g_wbuf, len, 1);

		UTF8TOWCHAR(filters[k].ext, (int)strlen(filters[k].ext)+1, g_wbuf, G_BUF_SIZE, len);
		if(!len)
			break;
		STR_APPEND(winfilts, g_wbuf, len, 1);
	}
	STR_APPEND(winfilts, 0, 2, 1);
	return winfilts;
}
ArrayHandle dialog_open_file(Filter *filters, int nfilters, int multiple)//TODO: multiple
{
	ArrayHandle winfilts=prep_filters(filters, nfilters), result=0;

	g_wbuf[0]=0;
	OPENFILENAMEW ofn=
	{
		sizeof(OPENFILENAMEW), ghWnd, 0,
		&WSTR_AT(winfilts, 0), 0, 0, 1,
		g_wbuf, G_BUF_SIZE,
		0, 0,//initial filename
		0,
		0,//dialog title
		OFN_CREATEPROMPT|OFN_PATHMUSTEXIST|OFN_NOCHANGEDIR,//flags
		0,//file offset
		0,//extension offset
		L"txt",//default extension
		0, 0,//data & hook
		0,//template name
		0,//reserved
		0,//reserved
		0,//flags ex
	};
	//if(multiple)//CRASHES
	//	ofn.Flags|=OFN_ALLOWMULTISELECT;
	int success=GetOpenFileNameW(&ofn);
	array_free(&winfilts);
	if(!success)
		return 0;

	int wlen=0;
	if(multiple)
	{
		for(;wlen<G_BUF_SIZE;++wlen)
			if(ofn.lpstrFile[wlen]=='\0')
				break;
	}
	else
		wlen=(int)wcslen(ofn.lpstrFile);

	int len=0;
	WCHARTOUTF8(ofn.lpstrFile, wlen, g_buf, G_BUF_SIZE, len);
	STR_COPY(result, g_buf, len);

	return result;
}
//const wchar_t		initialname[]=L"Untitled.txt";
char* dialog_save_file(Filter *filters, int nfilters, const char *initialname)
{
	ArrayHandle winfilts=prep_filters(filters, nfilters);
	
	int len0=(int)strlen(initialname), ext_offset=0;
	int len=0;
	wchar_t def_ext[16]={0};//default extension
	UTF8TOWCHAR(initialname, len0+1, g_wbuf, G_BUF_SIZE, len);
	for(ext_offset=len0-1;ext_offset>=0&&initialname[ext_offset]!='.';--ext_offset);
	memcpy(def_ext, g_wbuf+ext_offset, (len0+1-ext_offset)*sizeof(wchar_t));
	//memcpy(g_wbuf, initialname, sizeof(initialname));

	OPENFILENAMEW ofn=
	{
		sizeof(OPENFILENAMEW), ghWnd, 0,
		
		&WSTR_AT(winfilts, 0),	//<- filter
		
		0, 0,//custom filter & count
		1,								//<- initial filter index
		g_wbuf, G_BUF_SIZE,				//<- output filename
		0, 0,//initial filename
		0,
		0,//dialog title
		OFN_NOTESTFILECREATE|OFN_PATHMUSTEXIST|OFN_EXTENSIONDIFFERENT|OFN_OVERWRITEPROMPT,
		0, (unsigned short)ext_offset,				//<- file offset & extension offset
		def_ext,						//<- default extension (if user didn't type one)
		0, 0,//data & hook
		0,//template name
		0, 0,//reserved
		0,//flags ex
	};
	int success=GetSaveFileNameW(&ofn);
	array_free(&winfilts);
	if(!success)
		return 0;

	int retlen=(int)wcslen(ofn.lpstrFile)+1;
	char *ret=(char*)malloc(retlen+16);
	WCHARTOUTF8(ofn.lpstrFile, retlen, ret, retlen+16, len);
	if(!len)
		return 0;
	return ret;

	//int len=WideCharToMultiByte(CP_UTF8, 0, ofn.lpstrFile, wcslen(ofn.lpstrFile), g_buf, G_BUF_SIZE, 0, 0);	SYS_ASSERT(len);
	//if(!len)
	//	return 0;
	//g_buf[len]='\0';
	//return g_buf;
}

void get_window_title(char *buf, int len)
{
	GetWindowTextA(ghWnd, buf, len);
}
void set_window_title(const char *format, ...)
{
	va_list args;
	va_start(args, format);
	vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);
	int success=SetWindowTextA(ghWnd, g_buf);
	if(!success)
		LOG_ERROR("Error setting window title");
}

int copy_to_clipboard(const char *a, int size)
{
	char *clipboard=(char*)LocalAlloc(LMEM_FIXED, (size+1)*sizeof(char));
	if(!clipboard)
		return 0;
	//	LOG_ERROR("copy_to_clipboard(): LocalAlloc() error");

	memcpy(clipboard, a, (size+1)*sizeof(char));
	clipboard[size]='\0';

	OpenClipboard(ghWnd);
	EmptyClipboard();
	SetClipboardData(CF_OEMTEXT, (void*)clipboard);
	CloseClipboard();
	return 1;
}
ArrayHandle paste_from_clipboard(int loud)
{
	OpenClipboard(ghWnd);
	char *a=(char*)GetClipboardData(CF_OEMTEXT);
	if(!a)
	{
		CloseClipboard();
		if(loud)
			messagebox(MBOX_OK, "Error", "Failed to paste from clipboard");
		return 0;
	}
	int len0=(int)strlen(a);

	ArrayHandle ret;
	STR_COPY(ret, a, len0);

	CloseClipboard();

	return ret;
}

int copy_bmp_to_clipboard(const unsigned char *rgba, int iw, int ih)
{
	int res=iw*ih;
	ptrdiff_t size=sizeof(BITMAPINFO)+((size_t)res<<2);
	char *clipboard=(char*)LocalAlloc(LMEM_FIXED, size);
	if(!clipboard)
		return 0;
	BITMAPINFO bmi={{sizeof(BITMAPINFOHEADER), iw, -ih, 1, 32, BI_RGB, (DWORD)(res<<2), 0, 0, 0, 0}};
	memcpy(clipboard, &bmi, sizeof(BITMAPINFOHEADER));
	memcpy(clipboard+sizeof(BITMAPINFOHEADER), rgba, (size_t)res<<2);
	int success=OpenClipboard(ghWnd);
	if(!success)
		return 0;
	EmptyClipboard();
	SetClipboardData(CF_DIB, clipboard);
	CloseClipboard();
	return 1;
}
Image* paste_bmp_from_clipboard(void)
{
	unsigned clipboardformats[]={CF_DIB, CF_DIBV5, CF_BITMAP, CF_TIFF, CF_METAFILEPICT};
	int format=GetPriorityClipboardFormat(clipboardformats, sizeof clipboardformats>>2);
	if(format<=0)
		return 0;
	int success=OpenClipboard(ghWnd);
	if(!success)
		return 0;
	Image *image=0;
	void *hData=0;
	ptrdiff_t size=0;
	switch(format)
	{
	case CF_BITMAP://20231016
		{
			HBITMAP hBm2=(HBITMAP)GetClipboardData(CF_BITMAP);		SYS_ASSERT(hBm2);
			SIZE s;
			if(!hBm2)
				return 0;
			GetBitmapDimensionEx(hBm2, &s);
			if(!s.cx||!s.cy)
				break;

			image=image_from_uint8(0, s.cx, s.cy, 4, 8, 8, 8, 8);
			BITMAPINFO bmh={{sizeof(BITMAPINFOHEADER), s.cx, -s.cy, 1, 32, BI_RGB, 0}};
			HDC hDC2=CreateCompatibleDC(0);		SYS_ASSERT(hDC2);
			if(!hDC2)
				return 0;
			hBm2=(HBITMAP)SelectObject(hDC2, hBm2);
			int copied=GetDIBits(hDC2, hBm2, 0, s.cy, image->data, &bmh, DIB_RGB_COLORS);
			(void)copied;
			hBm2=(HBITMAP)SelectObject(hDC2, hBm2);
			DeleteDC(hDC2);
			for(ptrdiff_t k=((ptrdiff_t)image->iw*image->ih<<2)-1;k>=0;--k)//untested
				image->data[k]=((unsigned char*)image->data)[k]-128;
		}
		break;
	case CF_DIBV5://20210402
		{
			//hData=GetClipboardData(CF_DIB);		SYS_ASSERT(hData);
			//size=GlobalSize(hData);
			//BITMAPV5HEADER *bm5=(BITMAPV5HEADER*)GlobalLock(hData);
			//GlobalUnlock(hData);
			LOG_WARNING("Unsupported format CF_DIBV5");
		}
		break;
	case CF_DIB:
		hData=GetClipboardData(CF_DIB);						SYS_ASSERT(hData);
		if(hData)
		{
			size=GlobalSize(hData);
			(void)size;
			BITMAPINFO *bmi=(BITMAPINFO*)GlobalLock(hData);		SYS_ASSERT(bmi);
			if(bmi->bmiHeader.biCompression==BI_RGB||bmi->bmiHeader.biCompression==BI_BITFIELDS)//uncompressed
			{
				int iw=bmi->bmiHeader.biWidth, ih=abs(bmi->bmiHeader.biHeight);
				image=image_from_uint8(0, iw, ih, 4, 8, 8, 8, 8);
				byte *data=(byte*)bmi->bmiColors;
				int idx_b=0, idx_r=2;
				if(bmi->bmiHeader.biCompression==BI_BITFIELDS)//first 3 DWORDs are red, green, blue bitfields
				{
					if(((int*)data)[0]==0x000000FF)//swap: 0xAABBGGRR -> 0xAARRGGBB
						idx_b=2, idx_r=0, data+=3*(bmi->bmiHeader.biBitCount>>3);
					else if(((int*)data)[0]==0x00FF0000)//as it is: 0xAARRGGBB
						idx_b=0, idx_r=2, data+=3*(bmi->bmiHeader.biBitCount>>3);
					else
						LOG_WARNING("Invalid bitmap bitfield");
				}
				int data_size;
				switch(bmi->bmiHeader.biBitCount)
				{
				case 24:
					{
						int row_bytes=iw*3;
						row_bytes+=(4-(row_bytes&3))&3;
						data_size=row_bytes*ih;
						if(bmi->bmiHeader.biSizeImage)//skip palette bits
							data+=bmi->bmiHeader.biSizeImage-data_size;
						for(int ky=0;ky<ih&&row_bytes*(ky+1)<=data_size;++ky)
						{
							int *dstrow=image->data+((size_t)iw*ky<<2);
							byte *srcrow=data+row_bytes*(bmi->bmiHeader.biHeight<0?ky:bmi->bmiHeader.biHeight-1-ky);
							for(int kx=0;kx<iw;++kx)
							{
								byte *src=(byte*)(srcrow+3*kx);
								dstrow[kx<<2|0]=src[idx_r	]-128;
								dstrow[kx<<2|1]=src[1		]-128;
								dstrow[kx<<2|2]=src[idx_b	]-128;
								dstrow[kx<<2|3]=0;
							}
						}
						//for(int ky=0;ky<ih&&row_bytes*(ky+1)<=data_size;++ky)
						//{
						//	int *dstrow=(int*)image->data+iw*ky;
						//	byte *srcrow=data+row_bytes*(bmi->bmiHeader.biHeight<0?ky:bmi->bmiHeader.biHeight-1-ky);
						//	for(int kx=0;kx<iw;++kx)
						//	{
						//		byte *dst=(byte*)(dstrow+kx), *src=(byte*)(srcrow+3*kx);
						//		dst[0]=src[idx_r];
						//		dst[1]=src[1];
						//		dst[2]=src[idx_b];
						//		dst[3]=0xFF;
						//	}
						//}
						image->nch=3;
					}
					break;
				case 32://screenshots
					data_size=iw*ih<<2;
					if(bmi->bmiHeader.biSizeImage)//skip palette bits
						data+=bmi->bmiHeader.biSizeImage-data_size;
					//for(int ky=0;ky<ih&&(iw*(ky+1)<<2)<=data_size;++ky)
					for(int ky=0;ky<ih;++ky)
					{
						int *dst=(int*)image->data+((size_t)iw*ky<<2),
							*src=(int*)data+iw*(bmi->bmiHeader.biHeight<0?ky:bmi->bmiHeader.biHeight-1-ky);
						for(int kx=0;kx<iw;++kx)
						{
							unsigned char *p=(unsigned char*)(src+kx);
							dst[kx<<2|0]=p[2]-128;
							dst[kx<<2|1]=p[1]-128;
							dst[kx<<2|2]=p[0]-128;
							dst[kx<<2|3]=p[3]-128;
						}
					}
					//if(idx_r==0)
					//{
					//	for(int ky=0;ky<ih&&(iw*(ky+1)<<2)<=data_size;++ky)
					//		memcpy((int*)image->data+iw*ky, (int*)data+iw*(bmi->bmiHeader.biHeight<0?ky:bmi->bmiHeader.biHeight-1-ky), (size_t)iw<<2);
					//}
					//else//swap rb
					//{
					//	for(int ky=0;ky<ih&&(iw*(ky+1)<<2)<=data_size;++ky)
					//	{
					//		int *dst=(int*)image->data+iw*ky,
					//			*src=(int*)data+iw*(bmi->bmiHeader.biHeight<0?ky:bmi->bmiHeader.biHeight-1-ky);
					//		for(int kx=0;kx<iw;++kx)
					//		{
					//			unsigned char *p1=(unsigned char*)(src+kx), *p2=(unsigned char*)(dst+kx);
					//			p2[0]=p1[2];
					//			p2[1]=p1[1];
					//			p2[2]=p1[0];
					//			p2[3]=p1[3];
					//		}
					//	}
					//}
					image->nch=4;
					break;
				}
			}
			else
				LOG_WARNING("Unsupported compression: bmiHeader.biCompression == %d", bmi->bmiHeader.biCompression);
			GlobalUnlock(hData);
		}
		break;
	case CF_TEXT:
	case CF_SYLK:
	case CF_DIF:
	case CF_OEMTEXT:
	case CF_PENDATA:
	case CF_RIFF:
	case CF_WAVE:
	case CF_UNICODETEXT:
	case CF_ENHMETAFILE:
	case CF_HDROP:
	case CF_LOCALE:
		break;
	case CF_METAFILEPICT:
		LOG_WARNING("Unsupported format %d CF_METAFILEPICT", format);
		break;
	case CF_TIFF:
		LOG_WARNING("Unsupported format %d CF_TIFF", format);
		break;
	case CF_PALETTE:
		LOG_WARNING("Unsupported format %d CF_PALETTE", format);
		break;
	}
	CloseClipboard();
	return image;
}

static void update_main_key_states(void)
{
	keyboard[VK_CONTROL]=GET_KEY_STATE(VK_LCONTROL)<<1|GET_KEY_STATE(VK_RCONTROL);
	keyboard[VK_SHIFT]=GET_KEY_STATE(VK_LSHIFT)<<1|GET_KEY_STATE(VK_RSHIFT);
	keyboard[VK_MENU]=GET_KEY_STATE(VK_LMENU)<<1|GET_KEY_STATE(VK_RMENU);
}

void timer_start(int ms, int id)
{
	if(!timer)
		SetTimer(ghWnd, id, ms, 0), timer=1;
}
void timer_stop(int id)
{
	if(timer&&!g_repaint)
		KillTimer(ghWnd, id), timer=0;
}

void set_mouse(int x, int y)
{
	POINT p={x, y};
	ClientToScreen(ghWnd, &p);
	SetCursorPos(p.x, p.y);
	mouse_bypass=1;
}
void get_mouse(int *px, int *py)
{
	POINT p;
	GetCursorPos(&p);
	ScreenToClient(ghWnd, &p);
	*px=p.x, *py=p.y;
}
void show_mouse(int show)
{
	ShowCursor(show);
}

void swapbuffers(void)
{
	SwapBuffers(ghDC);
}
LRESULT __stdcall WndProc(HWND hWnd, unsigned message, WPARAM wParam, LPARAM lParam)
{
	switch(message)
	{
	case WM_SIZE:
		GetClientRect(ghWnd, &R);
		wndw=R.right-R.left, wndh=R.bottom-R.top;
		set_region_immediate(0, wndw, 0, wndh);
		//InvalidateRect(hWnd, 0, 0);
		io_resize();
		break;
	case WM_TIMER:
		g_repaint=0;
		io_timer();
		io_render();
		//prof_print();
		//SwapBuffers(ghDC);
		break;
	case WM_PAINT:
		g_repaint=0;
		io_render();
		//prof_print();
		//SwapBuffers(ghDC);
		break;
	case WM_ACTIVATE:
		update_main_key_states();
		break;

	case WM_MOUSEMOVE:
		mx=GET_X_LPARAM(lParam), my=GET_Y_LPARAM(lParam);
		if(mouse_bypass)
			mouse_bypass=0;
		else
		{
			if(io_mousemove())
				InvalidateRect(hWnd, 0, 0);
		}
		break;
	case WM_MOUSEWHEEL:
		if(io_mousewheel(GET_Y_LPARAM(wParam)))
			InvalidateRect(hWnd, 0, 0);
		break;

	case WM_LBUTTONDOWN:
	case WM_LBUTTONDBLCLK:
		if(io_keydn(KEY_LBUTTON, message==WM_LBUTTONDBLCLK))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_LBUTTON]=1;
		break;
	case WM_LBUTTONUP:
		if(io_keyup(KEY_LBUTTON, 0))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_LBUTTON]=0;
		break;
		
	case WM_MBUTTONDOWN:
	case WM_MBUTTONDBLCLK:
		if(io_keydn(KEY_MBUTTON, message==WM_LBUTTONDBLCLK))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_MBUTTON]=1;
		break;
	case WM_MBUTTONUP:
		if(io_keyup(KEY_MBUTTON, 0))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_MBUTTON]=0;
		break;
		
	case WM_RBUTTONDOWN:
	case WM_RBUTTONDBLCLK:
		if(io_keydn(KEY_RBUTTON, message==WM_LBUTTONDBLCLK))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_RBUTTON]=1;
		break;
	case WM_RBUTTONUP:
		if(io_keyup(KEY_RBUTTON, 0))
			InvalidateRect(hWnd, 0, 0);
		keyboard[VK_RBUTTON]=0;
		break;
		
	case WM_KEYDOWN:
	case WM_SYSKEYDOWN:
		if(wParam==KEY_F4&&keyboard[KEY_ALT])
		{
			if(!io_quit_request())
				return 0;
			PostQuitMessage(0);
		}
		else if(io_keydn((IOKey)wParam, 0))
			InvalidateRect(hWnd, 0, 0);
		keyboard[wParam]=GET_KEY_STATE((int)wParam);
		break;
	case WM_KEYUP:
	case WM_SYSKEYUP:
		if(io_keyup((IOKey)wParam, 0))
			InvalidateRect(hWnd, 0, 0);
		keyboard[wParam]=0;
		update_main_key_states();
		break;

	case WM_CLOSE:
		if(!io_quit_request())
			return 0;
		PostQuitMessage(0);
		break;
	}
	return DefWindowProcA(hWnd, message, wParam, lParam);
}
int __stdcall WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrev, _In_ char *pCmdLine, _In_ int nShowCmd)
{
	WNDCLASSEXA wndClassEx=
	{
		sizeof(WNDCLASSEXA), CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS,
		WndProc, 0, 0, hInstance,
		LoadIconA(0, (char*)0x00007F00),
		LoadCursorA(0, (char*)0x00007F00),

		0,
	//	(HBRUSH)(COLOR_WINDOW+1),

		0, "New format", 0
	};
	MSG msg;

	(void)hPrev;
	(void)pCmdLine;

	int len=GetModuleFileNameW(0, g_wbuf, G_BUF_SIZE);
	int len2=WideCharToMultiByte(CP_UTF8, 0, g_wbuf, len, g_buf, G_BUF_SIZE, 0, 0);
	STR_COPY(exedir, g_buf, len2);
	for(int k=(int)exedir->count-1;k>=0;--k)
	{
		char c=exedir->data[k];
		if(c=='/'||c=='\\')
		{
			STR_POPBACK(exedir, exedir->count-(k+1));
			break;
		}
	}

	int success=RegisterClassExA(&wndClassEx);	SYS_ASSERT(success);
	ghWnd=CreateWindowExA(WS_EX_ACCEPTFILES, wndClassEx.lpszClassName, "", WS_CAPTION|WS_SYSMENU|WS_THICKFRAME|WS_MINIMIZEBOX|WS_MAXIMIZEBOX|WS_CLIPCHILDREN, CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 0, 0, hInstance, 0);	SYS_ASSERT(ghWnd);//2023-04-11
	if(!ghWnd)
		return 0;

	GetClientRect(ghWnd, &R);
	wndw=R.right-R.left, wndh=R.bottom-R.top;
	ghDC=GetDC(ghWnd);
	gl_init();
	glClearColor(1, 1, 1, 1);
	if(!io_init(__argc-1, __argv+1))
		return EXIT_FAILURE;

	ShowWindow(ghWnd, nShowCmd);

	int ret;
	for(;;)
	{
		ret=GetMessageA(&msg, ghWnd, 0, 0);
		if(!ret||ret==-1)//GetMessageA returns 0 at WM_QUIT, and -1 when ghWnd was closed
			break;
		TranslateMessage(&msg);
		DispatchMessageA(&msg);
	}

		//finish
		io_cleanup();
		ReleaseDC(ghWnd, ghDC);

	return (int)msg.wParam;
}