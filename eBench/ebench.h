#pragma once
#ifndef INC_EBENCH_H
#define INC_EBENCH_H
#include"util.h"
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#include<GL/gl.h>
#ifdef __cplusplus
extern "C"
{
#endif
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4200)
#endif


//	#define ALLOW_OPENCL
//	#define DEBUG_MEMORY
//	#define ENABLE_CONSOLE_MAIN_TEST


#define TIMER_ID_KEYBOARD 1
#define TIMER_ID_MONITOR 2
#define TIMER_MONITOR_MS 100

extern int wndw, wndh, mx, my, mouse_bypass;
extern HWND ghWnd;
extern HDC ghDC;
extern HGLRC hRC;
extern char keyboard[256], timer;
extern int g_repaint;//set to 1 to paint again
extern RECT R;
extern wchar_t g_wbuf[G_BUF_SIZE];
extern ArrayHandle exedir;

extern float
	SN_x0, SN_x1, SN_y0, SN_y1,//screen-NDC conversion
	NS_x0, NS_x1, NS_y0, NS_y1;
extern float tdx, tdy;

typedef enum IOKeyEnum
{
#if defined _MSC_VER || defined _WINDOWS || defined _WIN32
#	define	IOKEY(LinVAL, VAL, LABEL)	KEY_##LABEL=VAL,
#elif defined __linux__
#	define	IOKEY(VAL, WinVAL, LABEL)	KEY_##LABEL=VAL,
#endif
	
//IO value, win32 value, label
IOKEY(0x00, 0x00, UNKNOWN)


//direct map keys
IOKEY(0x01, 0x01, LBUTTON)	//inserted
IOKEY(0x02, 0x04, MBUTTON)	//inserted
IOKEY(0x03, 0x02, RBUTTON)	//inserted

IOKEY(0x08, 0x08, BKSP)
IOKEY(0x09, 0x09, TAB)

IOKEY(0x0D, 0x0D, ENTER)

IOKEY(0x1B, 0x1B, ESC)

IOKEY(0x20, 0x20, SPACE)
IOKEY(0x30, 0x30, 0)
IOKEY(0x31, 0x31, 1)
IOKEY(0x32, 0x32, 2)
IOKEY(0x33, 0x33, 3)
IOKEY(0x34, 0x34, 4)
IOKEY(0x35, 0x35, 5)
IOKEY(0x36, 0x36, 6)
IOKEY(0x37, 0x37, 7)
IOKEY(0x38, 0x38, 8)
IOKEY(0x39, 0x39, 9)
IOKEY(0x41, 0x41, A)
IOKEY(0x42, 0x42, B)
IOKEY(0x43, 0x43, C)
IOKEY(0x44, 0x44, D)
IOKEY(0x45, 0x45, E)
IOKEY(0x46, 0x46, F)
IOKEY(0x47, 0x47, G)
IOKEY(0x48, 0x48, H)
IOKEY(0x49, 0x49, I)
IOKEY(0x4A, 0x4A, J)
IOKEY(0x4B, 0x4B, K)
IOKEY(0x4C, 0x4C, L)
IOKEY(0x4D, 0x4D, M)
IOKEY(0x4E, 0x4E, N)
IOKEY(0x4F, 0x4F, O)
IOKEY(0x50, 0x50, P)
IOKEY(0x51, 0x51, Q)
IOKEY(0x52, 0x52, R)
IOKEY(0x53, 0x53, S)
IOKEY(0x54, 0x54, T)
IOKEY(0x55, 0x55, U)
IOKEY(0x56, 0x56, V)
IOKEY(0x57, 0x57, W)
IOKEY(0x58, 0x58, X)
IOKEY(0x59, 0x59, Y)
IOKEY(0x5A, 0x5A, Z)


//other keys
IOKEY(0x13, 0x13, PAUSE)
IOKEY(0x14, 0x91, SCROLLLOCK)
IOKEY(0x15, 0x2C, PRINTSCR)

IOKEY(0x27, 0xDE, QUOTE)	//inserted '\'' with '\"' 0x22
IOKEY(0x2B, 0xBB, PLUS)		//inserted '+' with '=' 0x3D
IOKEY(0x2C, 0xBC, COMMA)	//inserted ',' with '<' 0x3C
IOKEY(0x2D, 0xBD, MINUS)	//inserted '-' with '_' 0x5F
IOKEY(0x2E, 0xBE, PERIOD)	//inserted '.' with '>' 0x3E
IOKEY(0x2F, 0xBF, SLASH)	//inserted '/' with '?' 0x3F

IOKEY(0x3B, 0xBA, SEMICOLON)	//inserted ';' with ':' 0x3A

IOKEY(0x5B, 0xDB, LBRACKET)	//inserted '[' with '{' 0x7B
IOKEY(0x5C, 0xDC, BACKSLASH)	//inserted '\\' with '|' 0x7C
IOKEY(0x5D, 0xDD, RBRACKET)	//inserted ']' with '}' 0x7D
IOKEY(0x60, 0xC0, GRAVEACCENT)	//inserted '`' with '~' 0x7E

IOKEY(0x7F, 0x2E, DEL)

IOKEY(0x84, 0x10, SHIFT)	//inserted
IOKEY(0x85, 0x11, CTRL)		//inserted
IOKEY(0x86, 0x12, ALT)		//inserted
IOKEY(0x87, 0x00, START)	//inserted

IOKEY(0xA0, 0x60, NP_0)
IOKEY(0xA1, 0x61, NP_1)
IOKEY(0xA2, 0x62, NP_2)
IOKEY(0xA3, 0x63, NP_3)
IOKEY(0xA4, 0x64, NP_4)
IOKEY(0xA5, 0x65, NP_5)
IOKEY(0xA6, 0x66, NP_6)
IOKEY(0xA7, 0x67, NP_7)
IOKEY(0xA8, 0x68, NP_8)
IOKEY(0xA9, 0x69, NP_9)
IOKEY(0xAA, 0x6A, NP_MUL)
IOKEY(0xAB, 0x6B, NP_PLUS)
IOKEY(0xAC, 0x6D, NP_MINUS)
IOKEY(0xAD, 0x6E, NP_PERIOD)
IOKEY(0xAE, 0x6F, NP_DIV)

IOKEY(0xBE, 0x70, F1)
IOKEY(0xBF, 0x71, F2)
IOKEY(0xC0, 0x72, F3)
IOKEY(0xC1, 0x73, F4)
IOKEY(0xC2, 0x74, F5)
IOKEY(0xC3, 0x75, F6)
IOKEY(0xC4, 0x76, F7)
IOKEY(0xC5, 0x77, F8)
IOKEY(0xC6, 0x78, F9)
IOKEY(0xC7, 0x79, F10)
IOKEY(0xC8, 0x7A, F11)
IOKEY(0xC9, 0x7B, F12)

IOKEY(0xD0, 0x24, HOME)
IOKEY(0xD1, 0x25, LEFT)
IOKEY(0xD2, 0x26, UP)
IOKEY(0xD3, 0x27, RIGHT)
IOKEY(0xD4, 0x28, DOWN)
IOKEY(0xD5, 0x21, PGUP)
IOKEY(0xD6, 0x22, PGDN)
IOKEY(0xD7, 0x23, END)

IOKEY(0xE5, 0x14, CAPSLOCK)

IOKEY(0xE3, 0x2D, INSERT)
IOKEY(0xEB, 0x03, BREAK)

IOKEY(0xFF, 0x90, NUMLOCK)

#undef	IOKEY
} IOKey;


//callbacks - implement these in your application:
int io_init(int argc, char **argv);//return false to abort
void io_resize(void);
int io_mousemove(void);//return true to redraw
int io_mousewheel(int forward);
int io_keydn(IOKey key, char c);
int io_keyup(IOKey key, char c);
void io_timer(void);
void io_render(void);
int io_quit_request(void);//return 1 to exit
void io_cleanup(void);//cleanup


void   console_start(void);
void   console_end(void);
void   console_buffer_size(short x, short y);
int    console_log(const char *format, ...);
void   console_pause(void);
int    console_scan(char *buf, int len);
int    console_scan_int(void);
double console_scan_float(void);

int sys_check(const char *file, int line, const char *info);
#define SYS_ASSERT(SUCCESS)	((void)((SUCCESS)!=0||sys_check(file, __LINE__, 0)))

typedef enum MessageBoxTypeEnum
{
	MBOX_OK,
	MBOX_OKCANCEL,
	MBOX_YESNOCANCEL,
} MessageBoxType;
int messagebox(MessageBoxType type, const char *title, const char *format, ...);//returns index of pressed button

typedef struct FilterStruct
{
	const char *comment, *ext;
} Filter;
ArrayHandle dialog_open_folder(void);
ArrayHandle dialog_open_file(Filter *filters, int nfilters, int multiple);
char* dialog_save_file(Filter *filters, int nfilters, const char *initialname);

void get_window_title(char *buf, int len);
void set_window_title(const char *format, ...);

int copy_to_clipboard(const char *a, int size);
ArrayHandle paste_from_clipboard(int loud);

#define GET_KEY_STATE(KEY)	(keyboard[KEY]=(GetAsyncKeyState(KEY)>>15)!=0)

void timer_start(int ms, int id);
void timer_stop(int id);

void set_mouse(int x, int y);
void get_mouse(int *px, int *py);
void show_mouse(int show);

void swapbuffers(void);


//OpenGL standard macros & types
#if 1
#define GL_FUNC_ADD		0x8006//GL/glew.h
#define GL_MIN			0x8007
#define GL_MAX			0x8008
#define GL_MAJOR_VERSION	0x821B
#define GL_MINOR_VERSION	0x821C
#define GL_TEXTURE0		0x84C0
#define GL_TEXTURE1		0x84C1
#define GL_TEXTURE2		0x84C2
#define GL_TEXTURE3		0x84C3
#define GL_TEXTURE4		0x84C4
#define GL_TEXTURE5		0x84C5
#define GL_TEXTURE6		0x84C6
#define GL_TEXTURE7		0x84C7
#define GL_TEXTURE8		0x84C8
#define GL_TEXTURE9		0x84C9
#define GL_TEXTURE10		0x84CA
#define GL_TEXTURE11		0x84CB
#define GL_TEXTURE12		0x84CC
#define GL_TEXTURE13		0x84CD
#define GL_TEXTURE14		0x84CE
#define GL_TEXTURE15		0x84CF
#define GL_TEXTURE_RECTANGLE	0x84F5
#define GL_PROGRAM_POINT_SIZE	0x8642
#define GL_BUFFER_SIZE		0x8764
#define GL_ARRAY_BUFFER		0x8892
#define GL_ELEMENT_ARRAY_BUFFER	0x8893
#define GL_STATIC_DRAW		0x88E4
#define GL_FRAGMENT_SHADER	0x8B30
#define GL_VERTEX_SHADER	0x8B31
#define GL_COMPILE_STATUS	0x8B81
#define GL_LINK_STATUS		0x8B82
#define GL_INFO_LOG_LENGTH	0x8B84
#define GL_DEBUG_OUTPUT		0x92E0//OpenGL 4.3+

#define GLFUNCLIST\
	GLFUNC(glBlendEquation)\
	GLFUNC(glBindVertexArray)\
	GLFUNC(glGenBuffers)\
	GLFUNC(glBindBuffer)\
	GLFUNC(glBufferData)\
	GLFUNC(glBufferSubData)\
	GLFUNC(glEnableVertexAttribArray)\
	GLFUNC(glVertexAttribPointer)\
	GLFUNC(glDisableVertexAttribArray)\
	GLFUNC(glCreateShader)\
	GLFUNC(glShaderSource)\
	GLFUNC(glCompileShader)\
	GLFUNC(glGetShaderiv)\
	GLFUNC(glGetShaderInfoLog)\
	GLFUNC(glCreateProgram)\
	GLFUNC(glAttachShader)\
	GLFUNC(glLinkProgram)\
	GLFUNC(glGetProgramiv)\
	GLFUNC(glGetProgramInfoLog)\
	GLFUNC(glDetachShader)\
	GLFUNC(glDeleteShader)\
	GLFUNC(glUseProgram)\
	GLFUNC(glGetAttribLocation)\
	GLFUNC(glDeleteProgram)\
	GLFUNC(glDeleteBuffers)\
	GLFUNC(glGetUniformLocation)\
	GLFUNC(glUniformMatrix3fv)\
	GLFUNC(glUniformMatrix4fv)\
	GLFUNC(glGetBufferParameteriv)\
	GLFUNC(glActiveTexture)\
	GLFUNC(glUniform1i)\
	GLFUNC(glUniform2i)\
	GLFUNC(glUniform1f)\
	GLFUNC(glUniform2f)\
	GLFUNC(glUniform3f)\
	GLFUNC(glUniform3fv)\
	GLFUNC(glUniform4f)\
	GLFUNC(glUniform4fv)
//	GLFUNC(glGenVertexArrays)
//	GLFUNC(glDeleteVertexArrays)
typedef void     (__stdcall *t_glBlendEquation)(unsigned mode);
typedef void     (__stdcall *t_glBindVertexArray)(unsigned arr);
typedef void     (__stdcall *t_glGenBuffers)(int n, unsigned *buffers);
typedef void     (__stdcall *t_glBindBuffer)(unsigned target, unsigned buffer);
typedef void     (__stdcall *t_glBufferData)(unsigned target, int size, const void *data, unsigned usage);
typedef void     (__stdcall *t_glBufferSubData)(unsigned target, int offset, int size, const void *data);
typedef void     (__stdcall *t_glEnableVertexAttribArray)(unsigned index);
typedef void     (__stdcall *t_glVertexAttribPointer)(unsigned index, int size, unsigned type, unsigned char normalized, int stride, const void *pointer);
typedef void     (__stdcall *t_glDisableVertexAttribArray)(unsigned index);
typedef unsigned (__stdcall *t_glCreateShader)(unsigned shaderType);
typedef void     (__stdcall *t_glShaderSource)(unsigned shader, int count, const char **string, const int *length);
typedef void     (__stdcall *t_glCompileShader)(unsigned shader);
typedef void     (__stdcall *t_glGetShaderiv)(unsigned shader, unsigned pname, int *params);
typedef void     (__stdcall *t_glGetShaderInfoLog)(unsigned shader, int maxLength, int *length, char *infoLog);
typedef unsigned (__stdcall *t_glCreateProgram)(void);
typedef void     (__stdcall *t_glAttachShader)(unsigned program, unsigned shader);
typedef void     (__stdcall *t_glLinkProgram)(unsigned program);
typedef void     (__stdcall *t_glGetProgramiv)(unsigned program, unsigned pname, int *params);
typedef void     (__stdcall *t_glGetProgramInfoLog)(unsigned program, int maxLength, int *length, char *infoLog);
typedef void     (__stdcall *t_glDetachShader)(unsigned program, unsigned shader);
typedef void     (__stdcall *t_glDeleteShader)(unsigned shader);
typedef void     (__stdcall *t_glUseProgram)(unsigned program);
typedef int      (__stdcall *t_glGetAttribLocation)(unsigned program, const char *name);
typedef void     (__stdcall *t_glDeleteProgram)(unsigned program);
typedef void     (__stdcall *t_glDeleteBuffers)(int n, const unsigned *buffers);
typedef int      (__stdcall *t_glGetUniformLocation)(unsigned program, const char *name);
typedef void     (__stdcall *t_glUniformMatrix3fv)(int location, int count, unsigned char transpose, const float *value);
typedef void     (__stdcall *t_glUniformMatrix4fv)(int location, int count, unsigned char transpose, const float *value);
typedef void     (__stdcall *t_glGetBufferParameteriv)(unsigned target, unsigned value, int *data);
typedef void     (__stdcall *t_glActiveTexture)(unsigned texture);
typedef void     (__stdcall *t_glUniform1i)(int location, int v0);
typedef void     (__stdcall *t_glUniform2i)(int location, int v0, int v1);
typedef void     (__stdcall *t_glUniform1f)(int location, float v0);
typedef void     (__stdcall *t_glUniform2f)(int location, float v0, float v1);
typedef void     (__stdcall *t_glUniform3f)(int location, float v0, float v1, float v2);
typedef void     (__stdcall *t_glUniform3fv)(int location, int count, const float *value);
typedef void     (__stdcall *t_glUniform4f)(int location, float v0, float v1, float v2, float v3);
typedef void     (__stdcall *t_glUniform4fv)(int location, int count, float *value);
//typedef void   (__stdcall *t_glGenVertexArrays)(int n, unsigned *arrays);//OpenGL 3.0
//typedef void   (__stdcall *t_glDeleteVertexArrays)(int n, unsigned *arrays);//OpenGL 3.0
#endif

#define GLFUNC(X) extern t_##X X;
GLFUNCLIST
#undef  GLFUNC

//Macros for 3D		DO NOT NEST THE MACROS LIKE f(g(x), h(y))
#if 1
#define		vec2_copy(DST, SRC)		(DST)[0]=(SRC)[0], (DST)[1]=(SRC)[1]
#define		vec2_add(DST, A, B)		(DST)[0]=(A)[0]+(B)[0], (DST)[1]=(A)[1]+(B)[1]
#define		vec2_sub(DST, A, B)		(DST)[0]=(A)[0]-(B)[0], (DST)[1]=(A)[1]-(B)[1]
#define		vec2_add1(DST, V, S)	(DST)[0]=(V)[0]+(S), (DST)[1]=(V)[1]+(S)
#define		vec2_sub1(DST, V, S)	(DST)[0]=(V)[0]-(S), (DST)[1]=(V)[1]-(S)
#define		vec2_mul1(DST, V, S)	(DST)[0]=(V)[0]*(S), (DST)[1]=(V)[1]*(S)
#define		vec2_div1(DST, V, S)	(DST)[0]=(V)[0]/(S), (DST)[1]=(V)[1]/(S)
#define		vec2_dot(A, B)			((A)[0]*(B)[0]+(A)[1]*(B)[1])
#define		vec2_cross(DST, A, B)	((A)[0]*(B)[1]-(A)[1]*(B)[0])
#define		vec2_abs(A)				sqrtf(vec2_dot(A, A))
#define		vec2_abs2(A)			vec2_dot(A, A)
#define		vec2_arg(A)				atan((A)[1]/(A)[0])
#define		vec2_arg2(A)			atan2((A)[1], (A)[0])
#define		vec2_eq(A, B)			((A)[0]==(B)[0]&&(A)[1]==(B)[1])
#define		vec2_ne(A, B)			((A)[0]!=(B)[0]||(A)[1]!=(B)[1])
#define		vec2_neg(DST, A)		(DST)[0]=-(A)[0], (DST)[1]=-(A)[1]

#define		mat2_mul_vec2(DST, M2, V2)	(DST)[0]=(M2)[0]*(V2)[0]+(M2)[1]*(V2)[1], (DST)[1]=(M2)[2]*(V2)[0]+(M2)[3]*(V2)[1]

#define		vec3_copy(DST, SRC)		(DST)[0]=(SRC)[0], (DST)[1]=(SRC)[1], (DST)[2]=(SRC)[2]
#define		vec3_set1(V3, GAIN)		(V3)[0]=(V3)[1]=(V3)[2]=GAIN
#define		vec3_setp(V3, POINTER)	(V3)[0]=(POINTER)[0], (V3)[1]=(POINTER)[1], (V3)[2]=(POINTER)[2]
#define		vec3_seti(V3, X, Y, Z)	(V3)[0]=X, (V3)[1]=Y, (V3)[2]=Z
#define		vec3_add(DST, A, B)		(DST)[0]=(A)[0]+(B)[0], (DST)[1]=(A)[1]+(B)[1], (DST)[2]=(A)[2]+(B)[2]
#define		vec3_sub(DST, A, B)		(DST)[0]=(A)[0]-(B)[0], (DST)[1]=(A)[1]-(B)[1], (DST)[2]=(A)[2]-(B)[2]
#define		vec3_add1(DST, V, S)	(DST)[0]=(V)[0]+(S), (DST)[1]=(V)[1]+(S), (DST)[2]=(V)[2]+(S)
#define		vec3_sub1(DST, V, S)	(DST)[0]=(V)[0]-(S), (DST)[1]=(V)[1]-(S), (DST)[2]=(V)[2]-(S)
#define		vec3_mul1(DST, V, S)	(DST)[0]=(V)[0]*(S), (DST)[1]=(V)[1]*(S), (DST)[2]=(V)[2]*(S)
#define		vec3_div1(DST, V, S)	(DST)[0]=(V)[0]/(S), (DST)[1]=(V)[1]/(S), (DST)[2]=(V)[2]/(S)
#define		vec3_dot(A, B)			((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2])
#define		vec3_cross(DST, A, B)\
	(DST)[0]=(A)[1]*(B)[2]-(A)[2]*(B)[1],\
	(DST)[1]=(A)[2]*(B)[0]-(A)[0]*(B)[2],\
	(DST)[2]=(A)[0]*(B)[1]-(A)[1]*(B)[0]
#define		vec3_triple_product(DST, A, B, C, TEMP_F1, TEMP_F2)\
	TEMP_F1=vec3_dot(A, C), TEMP_F2=vec3_dot(B, C), (DST)[0]=TEMP_F1*(B)[0]-TEMP_F2*(C)[0], (DST)[1]=TEMP_F1*(B)[1]-TEMP_F2*(C)[1], (DST)[2]=TEMP_F1*(B)[2]-TEMP_F2*(C)[2]
#define		vec3_abs(A)				sqrtf(vec3_dot(A, A))
#define		vec3_abs2(A)			vec3_dot(A, A)
#define		vec3_theta(A)			atan((A)[2]/sqrtf((A)[0]*(A)[0]+(A)[1]*(A)[1]))
#define		vec3_phi(A)				atan((A)[1]/(A)[0])
#define		vec3_phi2(A)			atan2((A)[1], (A)[0])
#define		vec3_isnan(A)			((A)[0]!=(A)[0]||(A)[1]!=(A)[1]||(A)[2]!=(A)[2])
#define		vec3_isnan_or_inf(A)	(vec3_isnan(A)||fabsf((A)[0])==infinity||fabsf((A)[1])==infinity||fabsf((A)[2])==infinity)
#define		vec3_eq(A, B)			((A)[0]==(B)[0]&&(A)[1]==(B)[1]&&(A)[2]==(B)[2])
#define		vec3_ne(A, B)			((A)[0]!=(B)[0]||(A)[1]!=(B)[1]||(A)[2]!=(B)[2])
#define		vec3_neg(DST, A)		(DST)[0]=-(A)[0], (DST)[1]=-(A)[1], (DST)[2]=-(A)[2]
#define		vec3_normalize(DST, A, TEMP_F)		TEMP_F=1/vec3_abs(A), vec3_div1(DST, A, TEMP_F)
#define		vec3_mix(DST, A, B, X)\
	(DST)[0]=(A)[0]+((B)[0]-(A)[0])*(X),\
	(DST)[1]=(A)[1]+((B)[1]-(A)[1])*(X),\
	(DST)[2]=(A)[2]+((B)[2]-(A)[2])*(X)

//column-major
#define		mat3_diag(MAT3, GAIN)	memset(MAT3, 0, 9*sizeof(float)), (MAT)[0]=(MAT)[1]=(MAT)[2]=GAIN
//#define	mat3_diag(MAT3, GAIN)	(MAT3)[0]=GAIN, (MAT3)[1]=0, (MAT3)[2]=0, (MAT3)[3]=0, (MAT3)[4]=GAIN, (MAT3)[5]=0, (MAT3)[6]=0, (MAT3)[7]=0, (MAT3)[8]=GAIN

#define		vec4_copy(DST, SRC)		(DST)[0]=(SRC)[0], (DST)[1]=(SRC)[1], (DST)[2]=(SRC)[2], (DST)[3]=(SRC)[3]
#define		vec4_dot(DST, A, B, TEMP_V1, TEMP_V2)		TEMP_V1=_mm_loadu_ps(A), TEMP_V1=_mm_mul_ps(TEMP_V1, _mm_loadu_ps(B)), TEMP_V1=_mm_hadd_ps(TEMP_V1, TEMP_V1), TEMP_V1=_mm_hadd_ps(TEMP_V1, TEMP_V1), _mm_store_ss(DST, TEMP_V1)
#define		vec4_add(DST, A, B)		(DST)[0]=(A)[0]+(B)[0], (DST)[1]=(A)[1]+(B)[1], (DST)[2]=(A)[2]+(B)[2], (DST)[3]=(A)[3]+(B)[3]
#define		vec4_sub(DST, A, B)		(DST)[0]=(A)[0]-(B)[0], (DST)[1]=(A)[1]-(B)[1], (DST)[2]=(A)[2]-(B)[2], (DST)[3]=(A)[3]-(B)[3]
#define		vec4_mul1(DST, V, S)	(DST)[0]=(V)[0]*(S), (DST)[1]=(V)[1]*(S), (DST)[2]=(V)[2]*(S), (DST)[3]=(V)[3]*(S)
//#define		vec4_add(DST, A, B)		{_mm_storeu_ps(DST, _mm_add_ps(_mm_loadu_ps(A), _mm_loadu_ps(B)));}
//#define		vec4_sub(DST, A, B)		{_mm_storeu_ps(DST, _mm_sub_ps(_mm_loadu_ps(A), _mm_loadu_ps(B)));}
//#define		vec4_mul1(DST, A, S)	{_mm_storeu_ps(DST, _mm_sub_ps(_mm_loadu_ps(A), _mm_set1_ps(S)));}

//column-major
#define		mat4_copy(DST, SRC)		memcpy(DST, SRC, 16*sizeof(float));
#define		mat4_identity(M4, GAIN)	memset(M4, 0, 16*sizeof(float)), (M4)[0]=(M4)[5]=(M4)[10]=(M4)[15]=GAIN
#define		mat4_data(M4, X, Y)		(M4)[(X)<<2|(Y)]
#define		mat4_mat3(DST, M4)\
	(DST)[0]=(M4)[0], (DST)[1]=(M4)[1], (DST)[2]=(M4)[2],\
	(DST)[3]=(M4)[4], (DST)[4]=(M4)[5], (DST)[5]=(M4)[6],\
	(DST)[6]=(M4)[8], (DST)[7]=(M4)[9], (DST)[8]=(M4)[10],
#define		mat4_transpose(DST, M4, TEMP_8V)\
	(TEMP_8V)[0]=_mm_loadu_ps(M4),\
	(TEMP_8V)[1]=_mm_loadu_ps((M4)+4),\
	(TEMP_8V)[2]=_mm_loadu_ps((M4)+8),\
	(TEMP_8V)[3]=_mm_loadu_ps((M4)+12),\
	(TEMP_8V)[4]=_mm_unpacklo_ps((TEMP_8V)[0], (TEMP_8V)[1]),\
	(TEMP_8V)[5]=_mm_unpacklo_ps((TEMP_8V)[2], (TEMP_8V)[3]),\
	(TEMP_8V)[6]=_mm_unpackhi_ps((TEMP_8V)[0], (TEMP_8V)[1]),\
	(TEMP_8V)[7]=_mm_unpackhi_ps((TEMP_8V)[2], (TEMP_8V)[3]),\
	_mm_storeu_ps(DST, _mm_movelh_ps((TEMP_8V)[4], (TEMP_8V)[5])),\
	_mm_storeu_ps((DST)+4, _mm_movehl_ps((TEMP_8V)[5], (TEMP_8V)[4])),\
	_mm_storeu_ps((DST)+8, _mm_movelh_ps((TEMP_8V)[6], (TEMP_8V)[7])),\
	_mm_storeu_ps((DST)+12, _mm_movehl_ps((TEMP_8V)[7], (TEMP_8V)[6]))
#define		mat4_mul_vec4(DST, M4, V4, TEMP_V)\
	TEMP_V=_mm_mul_ps(_mm_loadu_ps(M4), _mm_set1_ps((V4)[0])),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_mul_ps(_mm_loadu_ps(M4+4), _mm_set1_ps((V4)[1]))),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_mul_ps(_mm_loadu_ps(M4+8), _mm_set1_ps((V4)[2]))),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_mul_ps(_mm_loadu_ps(M4+12), _mm_set1_ps((V4)[3]))),\
	_mm_storeu_ps(DST, TEMP_V)
#define		mat4_mul_mat4(DST_NEW, M4A, M4B, TEMP_V)\
	mat4_mul_vec4(DST_NEW,		M4A, M4B,		TEMP_V),\
	mat4_mul_vec4(DST_NEW+4,	M4A, M4B+4,		TEMP_V),\
	mat4_mul_vec4(DST_NEW+8,	M4A, M4B+8,		TEMP_V),\
	mat4_mul_vec4(DST_NEW+12,	M4A, M4B+12,	TEMP_V)
#define		mat4_translate(M4, V3, TEMP_V)\
	TEMP_V=_mm_mul_ps(_mm_loadu_ps(M4), _mm_set1_ps((V3)[0])),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_mul_ps(_mm_loadu_ps(M4+4), _mm_set1_ps((V3)[1]))),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_mul_ps(_mm_loadu_ps(M4+8), _mm_set1_ps((V3)[2]))),\
	TEMP_V=_mm_add_ps(TEMP_V, _mm_loadu_ps(M4+12)),\
	_mm_storeu_ps(M4+12, TEMP_V)
#define		mat4_rotate(DST_NEW, M4, ANGLE, DIR, TEMP_VEC2, TEMP_VEC3A, TEMP_VEC3B)\
	vec3_normalize(TEMP_VEC3A, DIR, (TEMP_VEC2)[0]),\
	(TEMP_VEC2)[0]=cosf(ANGLE),\
	(TEMP_VEC2)[1]=1-(TEMP_VEC2)[0],\
	vec3_mul1(TEMP_VEC3B, TEMP_VEC3A, (TEMP_VEC2)[1]),\
	(TEMP_VEC2)[1]=sinf(ANGLE),\
	(DST_NEW)[0]=(TEMP_VEC3B)[0]*(TEMP_VEC3A)[0]+(TEMP_VEC2)[0],\
	(DST_NEW)[1]=(TEMP_VEC3B)[0]*(TEMP_VEC3A)[1]+(TEMP_VEC2)[1]*(TEMP_VEC3A)[2],\
	(DST_NEW)[2]=(TEMP_VEC3B)[0]*(TEMP_VEC3A)[2]-(TEMP_VEC2)[1]*(TEMP_VEC3A)[1],\
	(DST_NEW)[3]=0,\
	\
	(DST_NEW)[4]=(TEMP_VEC3B)[1]*(TEMP_VEC3A)[0]-(TEMP_VEC2)[1]*(TEMP_VEC3A)[2],\
	(DST_NEW)[5]=(TEMP_VEC3B)[1]*(TEMP_VEC3A)[1]+(TEMP_VEC2)[0],\
	(DST_NEW)[6]=(TEMP_VEC3B)[1]*(TEMP_VEC3A)[2]+(TEMP_VEC2)[1]*(TEMP_VEC3A)[0],\
	(DST_NEW)[7]=0,\
	\
	(DST_NEW)[8]=(TEMP_VEC3B)[2]*(TEMP_VEC3A)[0]+(TEMP_VEC2)[1]*(TEMP_VEC3A)[1],\
	(DST_NEW)[9]=(TEMP_VEC3B)[2]*(TEMP_VEC3A)[1]-(TEMP_VEC2)[1]*(TEMP_VEC3A)[0],\
	(DST_NEW)[10]=(TEMP_VEC3B)[2]*(TEMP_VEC3A)[2]+(TEMP_VEC2)[0],\
	(DST_NEW)[11]=0,\
	\
	(DST_NEW)[12]=0,\
	(DST_NEW)[13]=0,\
	(DST_NEW)[14]=0,\
	(DST_NEW)[15]=1
#define		mat4_scale(M4, AMMOUNT, TEMP_V0)\
		_mm_storeu_ps(M4, _mm_mul_ps(_mm_loadu_ps(M4), _mm_set1_ps((AMMOUNT)[0]))),\
		_mm_storeu_ps(M4+4, _mm_mul_ps(_mm_loadu_ps(M4+4), _mm_set1_ps((AMMOUNT)[1]))),\
		_mm_storeu_ps(M4+8, _mm_mul_ps(_mm_loadu_ps(M4+8), _mm_set1_ps((AMMOUNT)[2]))),\

typedef struct CameraStruct
{
	float
		x, y, z,//position
		ax, ay,//yaw/phi, pitch/theta
		tanfov,
		move_speed, turn_speed,
		cax, sax, cay, say;
} Camera;
#define cam_copy(DSTCAM, SRCCAM) memcpy(&(DSTCAM), &(SRCCAM), sizeof(DSTCAM))
#define cam_moveForward(CAM, SPEED) (CAM).x+=(SPEED)*(CAM).cax*(CAM).cay, (CAM).y+=(SPEED)*(CAM).sax*(CAM).cay, (CAM).z+=(SPEED)*(CAM).say
#define cam_moveBack(CAM, SPEED)    (CAM).x-=(SPEED)*(CAM).cax*(CAM).cay, (CAM).y-=(SPEED)*(CAM).sax*(CAM).cay, (CAM).z-=(SPEED)*(CAM).say
#define cam_moveLeft(CAM, SPEED)    (CAM).x-=(SPEED)*(CAM).sax, (CAM).y+=(SPEED)*(CAM).cax
#define cam_moveRight(CAM, SPEED)   (CAM).x+=(SPEED)*(CAM).sax, (CAM).y-=(SPEED)*(CAM).cax
#define cam_moveUp(CAM, SPEED)      (CAM).z+=SPEED
#define cam_moveDown(CAM, SPEED)    (CAM).z-=SPEED
#define cam_update_ax(CAM)          (CAM).ax=fmodf((CAM).ax, (float)(2*M_PI)), (CAM).cax=cosf((CAM).ax), (CAM).sax=sinf((CAM).ax)
#define cam_update_ay(CAM)          (CAM).ay=fmodf((CAM).ay, (float)(2*M_PI)), (CAM).cay=cosf((CAM).ay), (CAM).say=sinf((CAM).ay)
#define cam_turnUp(CAM, SPEED)      (CAM).ay+=(SPEED)*(CAM).turn_speed, cam_update_ay(CAM)
#define cam_turnDown(CAM, SPEED)    (CAM).ay-=(SPEED)*(CAM).turn_speed, cam_update_ay(CAM)
#define cam_turnLeft(CAM, SPEED)    (CAM).ax+=(SPEED)*(CAM).turn_speed, cam_update_ax(CAM)
#define cam_turnRight(CAM, SPEED)   (CAM).ax-=(SPEED)*(CAM).turn_speed, cam_update_ax(CAM)
#define cam_turnMouse(CAM, DX, DY, SENSITIVITY)\
	(CAM).ax-=(SENSITIVITY)*(CAM).turn_speed*(DX), cam_update_ax(CAM),\
	(CAM).ay-=(SENSITIVITY)*(CAM).turn_speed*(DY), cam_update_ay(CAM)
#define		cam_zoomIn(CAM, RATIO)  (CAM).tanfov/=RATIO, (CAM).turn_speed=(CAM).tanfov>1?1:(CAM).tanfov
#define		cam_zoomOut(CAM, RATIO) (CAM).tanfov*=RATIO, (CAM).turn_speed=(CAM).tanfov>1?1:(CAM).tanfov
#define		cam_accelerate(GAIN)    (CAM).move_speed*=GAIN

#define		cam_relworld2cam(CAM, DISP, DST_CP)\
	(DST_CP)[2]=(DISP)[0]*(CAM).cax+(DISP)[1]*(CAM).sax,\
	(DST_CP)[0]=(DISP)[0]*(CAM).sax-(DISP)[1]*(CAM).cax,\
	(DST_CP)[1]=(DST_CP)[2]*(CAM).say-(DISP)[2]*(CAM).cay,\
	(DST_CP)[2]=(DST_CP)[2]*(CAM).cay+(DISP)[2]*(CAM).say
#define		cam_world2cam(CAM, P, DST_CP, TEMP_3F)\
	vec3_sub(TEMP_3F, P, &(CAM).x),\
	cam_relworld2cam(CAM, TEMP_3F, DST_CP)
#define		cam_cam2screen(CAM, CP, DST_S, X0, Y0)\
	(DST_S)[1]=(X0)/((CP)[2]*(CAM).tanfov),\
	(DST_S)[0]=(X0)+(CP)[0]*(DST_S)[1],\
	(DST_S)[1]=(Y0)+(CP)[1]*(DST_S)[1]
#endif

#define screen2NDC_x(Xs)      (SN_x1*(Xs)+SN_x0)
#define screen2NDC_y(Ys)      (SN_y1*(Ys)+SN_y0)
#define screen2NDC_x_bias(Xs) (SN_x1*(Xs+0.5f)+SN_x0)
#define screen2NDC_y_bias(Ys) (SN_y1*(Ys+0.5f)+SN_y0)
#define NDC2screen_x(X)       (NS_x1*(X)+NS_x0)
#define NDC2screen_y(Y)       (NS_y1*(Y)+NS_y0)
#define NDC2screen_x_bias(X)  (NS_x1*(X)+NS_x0-0.5f)
#define NDC2screen_y_bias(Y)  (NS_y1*(Y)+NS_y0-0.5f)


//The Graphics API
void gl_init(void);

extern int error;
extern const char *gl_error_msg;
const char* glerr2str(int error);
#define GL_CHECK(E) (void)((E=glGetError())==0||log_error(file, __LINE__, 1, gl_error_msg, E, glerr2str(E)))

void set_region_immediate(int x1, int x2, int y1, int y2);//calls glViewport

void send_texture_pot(unsigned gl_texture, const int *rgba, int txw, int txh, int linear);
void send_texture_pot_grey(unsigned gl_texture, const unsigned char *bmp, int txw, int txh, int linear);
void send_texture_pot_int16x1(unsigned gl_texture, const unsigned *texture, int txw, int txh, int linear);
void select_texture(unsigned tx_id, int u_location);

void draw_line		(float x1, float y1, float x2, float y2, int color);
void draw_rect		(float x1, float x2, float y1, float y2, int color);
void draw_rect_hollow	(float x1, float x2, float y1, float y2, int color);
void draw_ellipse	(float x1, float x2, float y1, float y2, int color);
void draw_triangle	(float x1, float y1, float x2, float y2, float x3, float y3, int color);
#define        DRAW_LINEI(X1, Y1, X2, Y2, COLOR)        draw_line((float)(X1), (float)Y1, (float)X2, (float)Y2, COLOR)
#define        DRAW_RECTI(X1, X2, Y1, Y2, COLOR)        draw_rect((float)(X1), (float)X2, (float)Y1, (float)Y2, COLOR)
#define DRAW_RECT_HOLLOWI(X1, X2, Y1, Y2, COLOR) draw_rect_hollow((float)(X1), (float)X2, (float)Y1, (float)Y2, COLOR)
#define     DRAW_ELLIPSEI(X1, X2, Y1, Y2, COLOR)     draw_ellipse((float)(X1), (float)X2, (float)Y1, (float)Y2, COLOR)

void draw_rect_enqueue(ArrayHandle *vertices, float x1, float x2, float y1, float y2);
void draw_curve_enqueue(ArrayHandle *vertices, float x, float y);
void draw_2d_flush(ArrayHandle vertices, int color, unsigned primitive);

int toggle_sdftext(void);
int set_text_color(int color_txt);
int set_bk_color(int color_bk);
long long set_text_colors(long long colors);//0xBKBKBKBK_TXTXTXTX
float print_line_immediate(float tab_origin, float x, float y, float zoom, const char *msg, int msg_length, int req_cols, int *ret_idx, int *ret_cols);//returns text width
float GUIPrint(float tab_origin, float x, float y, float zoom, const char *format, ...);//returns text width
extern int g_printed;
float GUIPrint_append(float tab_origin, float x, float y, float zoom, int show, const char *format, ...);

float GUIPrint_enqueue(ArrayHandle *vertices, float tab_origin, float x, float y, float zoom, const char *format, ...);
float print_line_enqueue(ArrayHandle *vertices, float tab_origin, float x, float y, float zoom, const char *msg, int msg_length, int req_cols, int *ret_idx, int *ret_cols);
void print_line_flush(ArrayHandle vertices, float zoom);

void display_texture(int x1, int x2, int y1, int y2, unsigned txid, float alpha, float tx1, float tx2, float ty1, float ty2);
void display_texture_i(int x1, int x2, int y1, int y2, int *rgb, int txw, int txh, float tx1, float tx2, float ty1, float ty2, float alpha, int linear);

//3D
void mat4_lookAt(float *dst, const float *cam, const float *center, const float *up);
void mat4_FPSView(float *dst, const float *campos, float yaw, float pitch);
void mat4_perspective(float *dst, float tanfov, float w_by_h, float znear, float zfar);
void mat4_normalmat3(float *dst, float *m4);//inverse transpose of top left 3x3 submatrix

void draw_3D_triangles(Camera const *cam, unsigned vbuf, size_t offset, int nvertices, unsigned txid);

typedef struct GPUModelStruct
{
	unsigned VBO, EBO, txid;
	int n_elements, stride, vertices_start, normals_start, txcoord_start;
} GPUModel;
void gpubuf_send_VNT(GPUModel *dst, const float *VVVNNNTT, int n_floats, const int *indices, int n_ints);
void draw_L3D(Camera const *cam, GPUModel const *model, const float *modelpos, const float *lightpos, int lightcolor);

void draw_3d_line(Camera const *cam, const float *w1, const float *w2, int color);//world coordinates

void draw_3d_line_enqueue(ArrayHandle *vertices, float x, float y, float z);
void draw_3d_flush(ArrayHandle vertices, Camera const *cam, int *texture, int tw, int th, int linear, unsigned primitive);
void draw_3d_wireframe_gpu(Camera const *cam, unsigned gpubuf, int nvertices, int color, unsigned primitive);

void draw_contour3d_rect(Camera const *cam, unsigned vbuf, int nvertices, unsigned txid, float alpha);
void draw_contour3d(Camera const *cam, float x1, float x2, float y1, float y2, float z1, float z2, unsigned *textures, int count, float alpha);

void depth_test(int enable);

unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
unsigned short* stbi_load_16(char const *filename, int *x, int *y, int *comp, int req_comp);




//Entropy Benchmark (eBench)
typedef struct ImageStruct
{
	int iw, ih;
	short
		rct,
		nch;//{greyscale, greyscale+alpha, RGB, RGB+alpha}	alpha can be ignored for now
	char depth[4];
	char src_depth[4];//for entropy calculations
	int data[];//stride always sizeof(int[4])
} Image;
Image* image_load(const char *fn, int fnlen);
int image_save_uint8(const char *fn, Image const *image, int override_alpha);
Image* image_from_uint8(const unsigned char *src, int iw, int ih, int nch, char rdepth, char gdepth, char bdepth, char adepth);
Image* image_from_uint16(const unsigned short *src, int iw, int ih, int nch, char *src_depths, char *dst_depths);
void image_export_uint8(Image const *image, unsigned char **dst, int override_alpha, int swap_rb);
double image_getBMPsize(Image const *image);
size_t image_getbufsize(Image const *image);
void image_copy_nodata(Image **dst, Image const *src);
void image_copy(Image **dst, Image const *src);
int copy_bmp_to_clipboard(const unsigned char *rgba, int iw, int ih);
Image* paste_bmp_from_clipboard(void);//window.c


//extern Image *im0, *im1;
extern ArrayHandle fn;
extern double ch_entropy[4];
extern int loud_transforms;
extern int g_dist;
void update_image(void);


//transforms
void apply_transform(Image **image, int tid, int hasRCT);
void apply_selected_transforms(Image **pimage, int rct_only, int applyfwd, int applyinv);


//aux func
void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8);
double calc_entropy(const int *hist, int nlevels, int sum);//pass -1 if sum is unknown
int calc_maxdepth(Image const *image, int *inflation);
void calc_depthfromdata(int *image, int iw, int ih, char *depths, const char *src_depths);
void calc_csize_ans_separate(Image const *image, size_t *csizes);
void calc_csize_ans_energy(Image const *image, size_t *csizes);
void calc_csize_ans_energy_hybrid(Image const *image, size_t *csizes, const int *config);


//color transforms
void colortransform_YCoCg_R(Image *image, int fwd);
void colortransform_YCbCr_R_v1(Image *image, int fwd);
void colortransform_YCbCr_R_v2(Image *image, int fwd);
void colortransform_YCbCr_R_v3(Image *image, int fwd);
void colortransform_YCbCr_R_v4(Image *image, int fwd);
void colortransform_YCbCr_R_v5(Image *image, int fwd);
void colortransform_YCbCr_R_v6(Image *image, int fwd);
void colortransform_YCbCr_R_v7(Image *image, int fwd);
void colortransform_CrCgCb_R(Image *image, int fwd);
void colortransform_Pei09(Image *image, int fwd);
void colortransform_subg_opt(Image *image, int fwd);
void colortransform_JPEG2000(Image *image, int fwd);
void colortransform_JPEG2000_MA(Image *image, int fwd);
void colortransform_NBLI(Image *image, int fwd);
void colortransform_subtractgreen(Image *image, int fwd);
void colortransform_J2K2(Image *image, int fwd);
void rct_yrgb_v1(Image *image, int fwd);
void rct_yrgb_v2(Image *image, int fwd);
void ct_cmyk_fwd(Image *image);

void colortransform_lossy_YCbCr(Image *image, int fwd);
void colortransform_lossy_XYB(Image *image, int fwd);
void colortransform_lossy_matrix(Image *image, int fwd);
void colortransform_lossy_pca(Image *image, int fwd);

#define RCT_CUSTOM_NPARAMS 9
#define RCT_CUSTOM_PARAMBITS 12
extern short rct_custom_params[RCT_CUSTOM_NPARAMS];
void rct_custom_unpackpermutation(short p, unsigned char *permutation);
void rct_custom(Image *image, int fwd, const short *params);
void rct_custom_optimize(Image const *image, short *params);
void rct_custom_getmatrix(double *matrix, int fwd);

void rct_adaptive(Image *src, int fwd);


#define ENABLE_EXTENDED_RCT
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_COEFF_U_SUB_Y,
	II_COEFF_V_SUB_Y,
	II_COEFF_V_SUB_U,

	II_COUNT,
} RCTInfoIdx;
#ifdef ENABLE_EXTENDED_RCT
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
} OCHIndex;
#endif
#ifndef ENABLE_EXTENDED_RCT
typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_COUNT,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,
} OCHIndex;
#endif
#ifdef ENABLE_EXTENDED_RCT
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_400_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_400_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_400_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_400_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_040_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_040_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_040_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_004_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_004_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_400_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_400_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_400_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_400_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_040_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_040_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_040_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_004_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_004_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_400_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_400_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_400_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_400_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_040_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_040_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_040_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_004_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_004_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
#endif
#ifndef ENABLE_EXTENDED_RCT
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
int crct_analysis(Image *src);




//chroma from luma (CfL)

//void pred_cfl(char *buf, int iw, int ih, int fwd);


//pre-processors
void prep_BWT_x(Image **psrc, int fwd);
void prep_BWT_y(Image **psrc, int fwd);


//spatial transforms
void filt_deint422(Image *src);
void filt_deint420(Image *src);
void packsign(Image *src, int fwd);
void pred_cg_crct(Image *image, int fwd, int enable_ma);
void pred_clampgrad(Image *image, int fwd, int enable_ma);
void pred_cleartype(Image *src, int fwd);
void pred_sse(Image *src, int fwd);
void pred_artifact(Image *src, int fwd);
void pred_LeGallCG(Image *src, int fwd);
void pred_cgplus(Image *src, int fwd);
void pred_CG420(Image *src, int fwd);
void pred_CG422(Image *src, int fwd);
void pred_sel4(Image *src, int fwd);
void pred_select(Image *src, int fwd);
void pred_sub(Image *src, int fwd);
void pred_av2(Image *src, int fwd);
void pred_mix2(Image *src, int fwd);
void pred_mixN(Image *src, int fwd);
void pred_gray(Image *src, int fwd);
void pred_MTF(Image *src, int fwd);
void pred_palette(Image *src, int fwd);

void pred_CG3D(Image *src, int fwd, int enable_ma);
void pred_PU(Image *src, int fwd);
void pred_divfreeWP(Image *src, int fwd);
void pred_nblic(Image *src, int fwd);
void pred_wgrad(Image *src, int fwd, int hasRCT);
void pred_wgrad2(Image *src, int fwd);
void pred_wgrad3(Image *src, int fwd, int hasRCT);
void pred_wgrad4(Image *src, int fwd);
void pred_wgrad4c(Image *src, int fwd);
void pred_wgrad4c_crct(Image *src, int fwd);
void pred_wgrad5(Image *src, int fwd);
void pred_wgrad6(Image *src, int fwd);
void pred_wpred7(Image *src, int fwd);
void pred_wmix(Image *src, int fwd);
void pred_table(Image *src, int fwd);
void pred_convtest(Image *src);
void pred_lwav(Image *src, int fwd);
void pred_WPU(Image *src, int fwd);

void pred_av4(Image *src, int fwd, int enable_ma);
void pred_wp_deferred(Image *src, int fwd);
void pred_ecoeff(Image *src, int fwd, int enable_ma);
void pred_zipper(Image *src, int fwd, int enable_ma);
void pred_multistage(Image *src, int fwd, int enable_ma);
void pred_separate(Image *src, int fwd, int enable_ma);

void pred_dir(Image *src, int fwd, int enable_ma);

void pred_ols(Image *src, int fwd, int enable_ma);
void pred_ols2(Image *src, int fwd, int enable_ma);
void pred_ols3(Image *src, int fwd, int enable_ma);
#define OLS4_DEFAULT_PERIOD 128
#define OLS4_DEFAULT_LR 0.0018
//#define OLS4_DEFAULT_PERIOD 32
//#define OLS4_DEFAULT_LR 0.0003
#define OLS4_RMAX 4
#define OLS4_CTXSIZE (2*(OLS4_RMAX+1)*OLS4_RMAX)
extern int ols4_period;
extern double ols4_lr[4];
extern unsigned char ols4_cache;
extern unsigned char ols4_mask[4][OLS4_CTXSIZE+1];
void pred_ols4(Image *src, int period, double *lrs, unsigned char *mask0, unsigned char *mask1, unsigned char *mask2, unsigned char *mask3, int fwd);
void pred_ols5(Image *src, int fwd);
void pred_ols6(Image *src, int fwd);
void pred_ols7(Image *src, int fwd);		//L1
void pred_l1crct(Image *src, int fwd);		//
void pred_grfilt(Image *src, int fwd);		//
void pred_ols8(Image *src, int fwd);		//
void pred_ols8_crct(Image *src, int fwd);	//
void pred_ols9(Image *src, int fwd);

//void pred_select(Image const *src, Image *dst, int fwd, int enable_ma);

#define LIN_REACH 2
void pred_linear(Image const *src, Image *dst, const int *coeffs, int lgden, int fwd, int enable_ma);

void pred_t47(Image *src, int fwd, int enable_ma);

#define CUSTOM_REACH 2
#define CUSTOM_NNB (2*(CUSTOM_REACH+1)*CUSTOM_REACH)
#define CUSTOM_NPARAMS (2*(CUSTOM_REACH+1)*CUSTOM_REACH*6)
extern int custom_params[CUSTOM_NPARAMS];
extern int custom_clamp[4];//W NW N NE
extern int custom_pred_ch_idx;
void pred_custom(Image *src, int fwd, int enable_ma, const int *params);
void pred_custom_optimize(Image const *image, int *params, int reach);
void conv_custom(Image *src);

#define PW2_NPRED 20	//63
#define PW2_NPARAM (PW2_NPRED+11)
extern double pw2_errors[PW2_NPRED];
extern short pw2_params[PW2_NPARAM*3];
void pred_w2_opt_v2(Image *src, short *params, int loud);
void pred_w2_apply(Image *src, int fwd, int enable_ma, short *params);

extern short jxlparams_i16[33];
void pred_jxl_opt_v2(Image *src, short *params, int loud);
void pred_jxl_apply(Image *src, int fwd, int enable_ma, short *params);
void pred_jmj_apply(Image *src, int fwd, int enable_ma);

//CUSTOM3
#define C3_REACH 3
#define C3_NNB (C3_REACH*(C3_REACH+1)*4)		//for pixels & errors
#define C3_NPARAMS (C3_NNB*9+6)
typedef struct Custom3ParamsStruct
{
	short c00[C3_NNB  ], c01[C3_NNB  ], c02[C3_NNB];//fixed 1.14 bit
	short c10[C3_NNB+2], c11[C3_NNB  ], c12[C3_NNB];
	short c20[C3_NNB+2], c21[C3_NNB+2], c22[C3_NNB];
} Custom3Params;
extern Custom3Params c3_params;
void custom3_apply(Image *src, int fwd, int enable_ma, Custom3Params const *params);
void custom3_opt(Image const *src, Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss);

//CUSTOM4
extern short lossyconv_clipboard;
extern int lossyconv_page;
extern short lossyconv_params[5*5*4*4];//(r2 = 5*5) * 4 channels max * 4 banks
extern unsigned char lossyconv_stride[2*4], lossyconv_offset[2*4];//2 dimensions * 4 stages
extern unsigned char lossyconv_causalRCT[4];
void pred_lossyconv(Image *src);//irreversible conv
void pred_WC(Image *src);
void filt_median33(Image *src);
void filt_av33(Image *src);

void pred_calic(Image *src, int fwd, int enable_ma);
void pred_grad2(Image *src, int fwd, int enable_ma);

extern unsigned long long xoroshiro128_state[2];
#define XOROSHIRO128_RESET() xoroshiro128_state[0]=0xDF900294D8F554A5, xoroshiro128_state[1]=0x170865DF4B3201FC
unsigned long long xoroshiro128_next(void);



//DWTs			temp is size max(w, h)
typedef struct DWTSizeStruct
{
	unsigned short w, h;
} DWTSize;
ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override);//calculate dimensions of each DWT stage in descending order

void dwt2d_lazy_fwd	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_lazy_inv	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_haar_fwd	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_haar_inv	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_squeeze_fwd	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_squeeze_inv	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_legall53_fwd	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_legall53_inv	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_cdf97_fwd	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_cdf97_inv	(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
#if 0
void dwt2d_grad_fwd   (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_grad_inv   (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp);
void dwt2d_exp_fwd    (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, const double *params);
void dwt2d_exp_inv    (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, const double *params);
void dwt2d_custom_fwd (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, const double *params);
void dwt2d_custom_inv (int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, const double *params);

void dwt2d_dec_fwd(char *buffer, int iw, int ih);
void dwt2d_dec_inv(char *buffer, int iw, int ih);

void pred_wu97(char *buf, int iw, int ih, int fwd);
#endif


//DCTs

void image_dct4_fwd(Image *image);
void image_dct4_inv(Image *image);

void image_dct8_fwd(Image *image);
void image_dct8_inv(Image *image);

//Walsh–Hadamard transform
void image_wht4_fwd(Image *image);
void image_wht4_inv(Image *image);
void image_wht8_fwd(Image *image);
void image_wht8_inv(Image *image);
void image_haar8_fwd(Image *image);
void image_haar8_inv(Image *image);
void image_fdct8_fwd(Image *image);
void image_fdct8_inv(Image *image);

#if 0
void image_split_fwd(char *image, int iw, int ih);
void image_split_inv(char *image, int iw, int ih);

void pred_dct3_fwd(char *buf, int iw, int ih, int nch, int bytestride);
void pred_dct3_inv(char *buf, int iw, int ih, int nch, int bytestride);

double calc_entropy(int *hist, int sum);//pass -1 if sum is unknown
void channel_entropy(unsigned char *buf, int resolution, int nch, int bytestride, float *cr, int *usage);
void jointhistogram(unsigned char *buf, int iw, int ih, int nbits, ArrayHandle *hist, int space_not_color);


//experiment 24: test 16 optimization per block
#if 0
typedef struct E24ParamsStruct
{
	short
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} E24Params;
extern E24Params e24_params[3];
extern double e24_cr[3];
int e24_optimizeall(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int loud);
void e24_estimate(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2);
#endif

//Bayes memory experiment
typedef struct BayesCounterStruct
{
	int n[2];
} BayesCounter;
extern ArrayHandle bayes_mem[8];
void bayes_estimate(unsigned char *image, int iw, int ih, int x1, int x2, int y1, int y2, int kc);


//transforms pt2

void pred_learned(char *buf, int iw, int ih, int fwd);
void pred_learned_v2(char *buf, int iw, int ih, int fwd);
void pred_learned_v3(char *buf, int iw, int ih, int fwd);
void pred_learned_v4(char *buf, int iw, int ih, int fwd);//kodim13 2.007026


//CPU transforms
void pred_learned_cpu(char *buf, int iw, int ih, int fwd);


//GPU trainsforms
#ifdef ALLOW_OPENCL
void pred_learned_gpu(char *buf, int iw, int ih, int fwd);
#endif


void pred_c03(char *src, int iw, int ih, int fwd);
void pred_c10(char *src, int iw, int ih, int fwd);
void pred_c20(char *src, int iw, int ih, int fwd);
#endif


typedef enum EContextEnum
{
	ECTX_QNW,
	ECTX_MIN_QN_QW,	//#2  prev best
	ECTX_MAX_QN_QW,
	ECTX_MIN_N_W_NW_NE,
	ECTX_ARGMIN_N_W_NW_NE,
	ECTX_ABAC0,
	ECTX_ABAC1,
	ECTX_GRCTX,	//#1  the best by far
	ECTX_GR,
	ECTX_STATIC_O0,
	ECTX_STATIC_O1,
	ECTX_INTERLEAVED,
	ECTX_DWT,
	ECTX_YUV422,
	ECTX_YUV420,
	ECTX_ZERO,

	ECTX_COUNT,
} EContext;
const char* ec_method_label(EContext ec_method);
void calc_csize_ec(Image const *src, EContext method, int adaptive, int expbits, int msb, int lsb, double *entropy);
void calc_csize_abac(Image const *src, int order, double *entropy);

void measure_ssim_avg(const Image *im0, const Image *im1, double *ret_ssim);
void measure_ssim_map(const Image *im0, Image *im1);


#ifdef _MSC_VER
#pragma warning(pop)
#endif
#ifdef __cplusplus
}
#endif
#endif//INC_EBENCH_H
