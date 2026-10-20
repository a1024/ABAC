#include"ebench.h"
#include<stdio.h>
#include<tmmintrin.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
#ifdef _MSC_VER
#pragma comment(lib, "OpenGL32.lib")
#endif

#define GLFUNC(X) t_##X X=0;
GLFUNCLIST
#undef  GLFUNC

int error=0;
const char *GLversion=0;
int rx0=0, ry0=0, rdx=0, rdy=0;/*current OpenGL region*/
float
	SN_x0=0, SN_x1=0, SN_y0=0, SN_y1=0,/*screen-NDC conversion*/
	NS_x0=0, NS_x1=0, NS_y0=0, NS_y1=0;
unsigned vertex_buffer=0;

unsigned current_program=0;
typedef struct ShaderVarStruct
{
	int *pvar;/*initialize to -1*/
	const char *name;
} ShaderVar;
typedef struct ShaderProgramStruct
{
	const char *name,/*program name for error reporting*/
		*vsrc, *fsrc;
	ShaderVar *attributes, *uniforms;
	int n_attr, n_unif;
	unsigned program;/*initialize to 0*/
} ShaderProgram;

/*text globals*/
char sdf_available=0, sdf_active=0;
short tab_count=8;
float tdx=0, tdy=0;/*non-tab character dimensions*/
float sdf_dx=0, sdf_dy=0, sdf_txh=0;
float font_zoom=1, font_zoom_min=1, font_zoom_max=32, sdf_dzoom=1.02f, sdf_slope=0.062023f;
typedef struct _QuadCoords
{
	float x1, x2, y1, y2;
} QuadCoords;
QuadCoords font_coords[128-32]={0}, sdf_glyph_coords[128-32]={0};
unsigned   font_txid=0, sdf_atlas_txid=0;
typedef struct SDFTextureHeaderStruct
{
	double slope;
	char
		grid_start_x, grid_start_y,
		cell_size_x, cell_size_y,
		csize_x, csize_y,
		reserved[2];
} SDFTextureHeader;
uint64_t colors_text=(uint64_t)0xFFABABAB<<32|0xFF000000;/*0xBKBKBKBK_TXTXTXTX*/

/*shaders*/
#define SHADER_LIST SHADER(2D) SHADER(texture) SHADER(text) SHADER(sdftext) SHADER(3D) SHADER(L3D) SHADER(contour3D)
#if 1
/*shader_2D*/
#define ATTR_2D ATTR(2D, coords)
#define UNIF_2D UNIF(2D, color)
const char
	src_vert_2D[]=
		"#version 120\n"
		"attribute vec2 a_coords;\n"		/*attributes: a_coords*/
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(a_coords, 0., 1.);\n"
		"}",
	src_frag_2D[]=
		"#version 120\n"
		"uniform vec4 u_color;\n"		/*uniforms: u_color*/
		"void main()\n"
		"{\n"
		"    gl_FragColor=u_color;\n"
		"}";

/*shader_texture*/
#define ATTR_texture ATTR(texture, coords)
#define UNIF_texture UNIF(texture, texture) UNIF(texture, alpha)
const char
	src_vert_texture[]=
		"#version 120\n"
		"attribute vec4 a_coords;"		/*attributes: a_coords*/
		"varying vec2 v_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(a_coords.xy, 0., 1.);\n"
		"    v_texcoord=a_coords.zw;\n"
		"}",
	src_frag_texture[]=
		"#version 120\n"
		"varying vec2 v_texcoord;\n"
		"uniform sampler2D u_texture;\n"	/*uniforms: u_texture, u_alpha*/
		"uniform float u_alpha;\n"
		"void main()\n"
		"{\n"
		"    gl_FragColor=texture2D(u_texture, v_texcoord);\n"
		"    gl_FragColor.a*=u_alpha;\n"
		"}";

/*shader_text*/
#define ATTR_text ATTR(text, coords)
#define UNIF_text UNIF(text, atlas) UNIF(text, txtColor) UNIF(text, bkColor)
const char
	src_vert_text[]=
		"#version 120\n"
		"attribute vec4 a_coords;"		/*attributes: a_coords*/
		"varying vec2 v_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(a_coords.xy, 0., 1.);\n"
		"    v_texcoord=a_coords.zw;\n"
		"}",
	src_frag_text[]=
		"#version 120\n"
		"varying vec2 v_texcoord;\n"
		"uniform sampler2D u_atlas;\n"		/*uniforms: u_atlas, u_txtColor, u_bkColor*/
		"uniform vec4 u_txtColor, u_bkColor;\n"
		"void main()\n"
		"{\n"
		"    vec4 region=texture2D(u_atlas, v_texcoord);\n"
		"    gl_FragColor=mix(u_txtColor, u_bkColor, region.r);\n"/*u_txtColor*(1-region.r) + u_bkColor*region.r*/
		"}";

/*shader_sdftext*/
#define ATTR_sdftext ATTR(sdftext, coords)
#define UNIF_sdftext UNIF(sdftext, atlas) UNIF(sdftext, txtColor) UNIF(sdftext, bkColor) UNIF(sdftext, zoom)
const char
	src_vert_sdftext[]=
		"#version 120\n"
		"attribute vec4 a_coords;"		/*attributes: a_coords*/
		"varying vec2 v_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(a_coords.xy, 0., 1.);\n"
		"    v_texcoord=a_coords.zw;\n"
		"}",
	src_frag_sdftext[]=
		"#version 120\n"
		"varying vec2 v_texcoord;\n"
		"uniform sampler2D u_atlas;\n"		/*uniforms: u_atlas, u_txtColor, u_bkColor, u_zoom*/
		"uniform vec4 u_txtColor, u_bkColor;\n"
		"uniform float u_zoom;\n"
		"void main()\n"
		"{\n"
		"    vec4 region=texture2D(u_atlas, v_texcoord);\n"

		"    float temp=clamp(u_zoom*(0.5f+0.45f/u_zoom-region.r), 0, 1);\n"
	/*	"    float temp=clamp(u_zoom*(0.5f+0.001f*u_zoom-region.r), 0, 1);\n"*/
		"    gl_FragColor=mix(u_txtColor, u_bkColor, temp);\n"

	/*	"    gl_FragColor=region.r>=0.5f?u_txtColor:u_bkColor;\n"*/ /*no anti-aliasing*/
		"}";

/*shader_3D*/
#define ATTR_3D ATTR(3D, vertex) ATTR(3D, texcoord)
#define UNIF_3D UNIF(3D, matrix) UNIF(3D, texture)
const char
	src_vert_3D[]=
		"#version 120\n"
		"uniform mat4 u_matrix;\n"		/*attributes: a_vertex, a_texcoord*/
		"attribute vec3 a_vertex;\n"
		"attribute vec2 a_texcoord;\n"
		"varying vec2 v_texcoord;\n"
		"varying vec4 v_glposition;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=u_matrix*vec4(a_vertex, 1.);\n"
		"    v_glposition=gl_Position;\n"
		"    gl_Position.z=0.;\n"
		"    v_texcoord=a_texcoord;\n"
		"}",
	src_frag_3D[]=
		"#version 120\n"
		"varying vec2 v_texcoord;\n"		/*uniforms: u_matrix, u_texture*/
		"varying vec4 v_glposition;\n"
		"uniform sampler2D u_texture;\n"/*use 1x1 texture for solid color*/
		"void main()\n"
		"{\n"
		"	 gl_FragColor=texture2D(u_texture, v_texcoord);\n"/*alpha is in the texture*/
		"    gl_FragDepth=(-(1000.+0.1)*(-v_glposition.w)-2.*1000.*0.1)/((1000.-0.1)*v_glposition.w);\n"/*USE mProj.znear=0.1f, mProj.zfar=1000.f*/
		"}";

/*shader_L3D*/
#define ATTR_L3D ATTR(L3D, vertex) ATTR(L3D, normal) ATTR(L3D, texcoord)
#define UNIF_L3D UNIF(L3D, matVP_Model) UNIF(L3D, matNormal) UNIF(L3D, texture) UNIF(L3D, sceneInfo)
const char
	src_vert_L3D[]=
		"#version 120\n"
		"uniform mat4 u_matVP_Model[2];\n"	/*attributes: a_vertex, a_normal, a_texcoord*/
		"uniform mat3 u_matNormal;\n"
		"attribute vec3 a_vertex;\n"
		"attribute vec3 a_normal;\n"
		"attribute vec2 a_texcoord;\n"
		"varying vec3 v_fragpos;\n"
		"varying vec3 v_normal;\n"
		"varying vec2 v_texcoord;\n"
		"varying vec4 v_glposition;\n"
		"void main()\n"
		"{\n"
		"    vec4 fullpos=vec4(a_vertex, 1.);\n"
		"    gl_Position=u_matVP_Model[0]*fullpos;\n"
		"    v_glposition=gl_Position;\n"
		"    gl_Position.z=0.;\n"
		"    v_fragpos=vec3(u_matVP_Model[1]*fullpos);\n"
		"    v_normal=u_matNormal*a_normal;\n"
		"    v_texcoord=a_texcoord;\n"
		"}",
	src_frag_L3D[]=
		"#version 120\n"
		"varying vec3 v_fragpos;\n"		/*uniforms: u_matVP_Model, u_matNormal,*/
		"varying vec3 v_normal;\n"		/*	u_texture, u_sceneInfo*/
		"varying vec4 v_glposition;\n"
		"varying vec2 v_texcoord;\n"
		"uniform sampler2D u_texture;\n"/*use 1x1 texture for solid color*/
		"uniform vec3 u_sceneInfo[3];\n"
		"void main()\n"
		"{\n"
		"    vec3 lightPos=u_sceneInfo[0], lightColor=u_sceneInfo[1], viewPos=u_sceneInfo[1];\n"
		"	 vec4 objectColor=texture2D(u_texture, v_texcoord);\n"/*alpha is in the texture*/

		"    vec3 normal=normalize(v_normal);\n"
		"    vec3 lightdir=normalize(lightPos-v_fragpos);\n"
				
		"    float specularstrength=0.5;\n"
		"    vec3 viewdir=normalize(viewPos-v_fragpos), reflectdir=reflect(-lightdir, normal);\n"
		"    vec3 specular=specularstrength*lightColor*pow(max(dot(viewdir, reflectdir), 0.), 32);\n"

		"    vec3 diffuse=max(dot(normal, lightdir), 0.)*lightColor;\n"
		"    gl_FragColor=vec4((0.1*lightColor+diffuse+specular)*objectColor.rgb, objectColor.a);\n"

		"    gl_FragDepth=(-(1000.+0.1)*(-v_glposition.w)-2.*1000.*0.1)/((1000.-0.1)*v_glposition.w);\n"/*USE mProj.znear=0.1f, mProj.zfar=1000.f*/
		"}";

/*shader_contour3D*/
#define ATTR_contour3D ATTR(contour3D, vertex) ATTR(contour3D, texcoord)
#define UNIF_contour3D UNIF(contour3D, matrix) UNIF(contour3D, alpha) UNIF(contour3D, texture)
const char
	src_vert_contour3D[]=
		"#version 120\n"
		"uniform mat4 u_matrix;\n"			/*attributes: a_vertex, a_texcoord*/
		"attribute vec3 a_vertex;\n"
		"attribute vec2 a_texcoord;\n"
		"varying vec2 v_texcoord;\n"
		"varying vec4 v_glposition;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=u_matrix*vec4(a_vertex, 1.);\n"
		"    v_glposition=gl_Position;\n"
		"    gl_Position.z=0.;\n"
		"    v_texcoord=a_texcoord;\n"
		"}",
	src_frag_contour3D[]=
		"#version 120\n"
		"varying vec2 v_texcoord;\n"		/*uniforms: u_matrix, u_alpha, u_texture*/
		"varying vec4 v_glposition;\n"
		"uniform float u_alpha;\n"
		"uniform sampler2D u_texture;\n"/*use 1x1 texture for solid color*/
		"void main()\n"
		"{\n"
		"    float val=texture2D(u_texture, v_texcoord).r;\n"
		"    gl_FragColor=vec4(0.5+0.5*cos(20*val-4.1887902), 0.5+0.5*cos(20*val-2.0943951), 0.5+0.5*cos(20*val), u_alpha*sqrt(abs(val)));\n"
	/*	"    gl_FragColor=vec4(val, val, u_alpha, 1.);\n"*/
		"    gl_FragDepth=(-(1000.+0.1)*(-v_glposition.w)-2.*1000.*0.1)/((1000.-0.1)*v_glposition.w);\n"/*USE mProj.znear=0.1f, mProj.zfar=1000.f*/
		"}";
#endif

/*shader declarations*/
#define ATTR(NAME, LABEL) int a_##NAME##_##LABEL=-1;
#define SHADER(NAME) ATTR_##NAME
SHADER_LIST
#undef SHADER
#undef ATTR

#define UNIF(NAME, LABEL) int u_##NAME##_##LABEL=-1;
#define SHADER(NAME) UNIF_##NAME
SHADER_LIST
#undef SHADER
#undef UNIF

#define ATTR(NAME, LABEL) {&a_##NAME##_##LABEL, "a_" #LABEL},
#define SHADER(NAME) ShaderVar attr_##NAME[]={ATTR_##NAME};
SHADER_LIST
#undef SHADER
#undef ATTR

#define UNIF(NAME, LABEL) {&u_##NAME##_##LABEL, "u_" #LABEL},
#define SHADER(NAME) ShaderVar unif_##NAME[]={UNIF_##NAME};
SHADER_LIST
#undef SHADER
#undef UNIF

#define SHADER(NAME) ShaderProgram shader_##NAME={"shader_" #NAME, src_vert_##NAME, src_frag_##NAME, attr_##NAME, unif_##NAME, _countof(attr_##NAME), _countof(unif_##NAME), 0};
SHADER_LIST
#undef SHADER

const char* glerr2str(int _error)
{
#define 			EC(x)	case x:a=(const char*)#x;break
	const char *a=0;
	switch(_error)
	{
	case 0:a="SUCCESS";break;
	EC(GL_INVALID_ENUM);
	EC(GL_INVALID_VALUE);
	EC(GL_INVALID_OPERATION);
	case 0x0503:a="GL_STACK_OVERFLOW";break;
	case 0x0504:a="GL_STACK_UNDERFLOW";break;
	EC(GL_OUT_OF_MEMORY);
	case 0x0506:a="GL_INVALID_FRAMEBUFFER_OPERATION";break;
	case 0x0507:a="GL_CONTEXT_LOST";break;
	case 0x8031:a="GL_TABLE_TOO_LARGE";break;
	default:a="???";break;
	}
	return a;
#undef				EC
}
const char *gl_error_msg="GL %d: %s";
#if 0
#define CHECK_GL(E) (void)((E=glGetError())==0||crash(__FILE__, __LINE__, 1, gl_error_msg, E, glerr2str(E)))
#endif

/*shader API*/
static unsigned CompileShader(const char *src, unsigned type, const char *programname)
{
	int success=0;
	unsigned shaderID=glCreateShader(type);
	glShaderSource(shaderID, 1, &src, 0);
	glCompileShader(shaderID);
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		char *errorMessage;

		glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		errorMessage=(char*)malloc(infoLogLength+1);
		if(errorMessage)
		{
			glGetShaderInfoLog(shaderID, infoLogLength, 0, errorMessage);
			copy_to_clipboard(errorMessage, infoLogLength);
			free(errorMessage);
		}
		else
		{
			glGetShaderInfoLog(shaderID, G_BUF_SIZE, 0, g_buf);
			copy_to_clipboard(g_buf, MINVAR(G_BUF_SIZE, infoLogLength));
		}
		if(programname)
			crash(__FILE__, __LINE__, 1, "%s shader compilation failed. Output copied to clipboard.", programname);
		else
			CRASH("Shader compilation failed. Output copied to clipboard.");
		return 0;
	}
	return shaderID;
}
static unsigned make_gl_program_impl(const char *vertSrc, const char *fragSrc, const char *programname)
{
	int success=0;
	unsigned
		vertShaderID=CompileShader(vertSrc, GL_VERTEX_SHADER, programname),
		fragShaderID=CompileShader(fragSrc, GL_FRAGMENT_SHADER, programname);

#if 0
	if(!vertShaderID||!fragShaderID)
		return 0;
#endif
	unsigned ProgramID=glCreateProgram();
	glAttachShader(ProgramID, vertShaderID);
	glAttachShader(ProgramID, fragShaderID);
	glLinkProgram(ProgramID);

	glGetProgramiv(ProgramID, GL_LINK_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		char *errorMessage;

		glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &infoLogLength);
		errorMessage=(char*)malloc(infoLogLength+(size_t)1);
		glGetProgramInfoLog(ProgramID, infoLogLength, 0, errorMessage);
		copy_to_clipboard(errorMessage, infoLogLength);
		if(programname)
			crash(__FILE__, __LINE__, 1, "%s shader link failed. Output copied to cipboard.", programname);
		else
			CRASH("Shader link failed. Output copied to cipboard.");
		free(errorMessage);
		return 0;
	}
	glDetachShader(ProgramID, vertShaderID);
	glDetachShader(ProgramID, fragShaderID);
	glDeleteShader(vertShaderID);
	glDeleteShader(fragShaderID);

	CHECK_GL(error);
	return ProgramID;
}
static int make_gl_program(ShaderProgram *p)
{
	int ka, ku;

	p->program=make_gl_program_impl(p->vsrc, p->fsrc, p->name);
	if(!p->program)
		return 0;
	for(ka=0;ka<p->n_attr;++ka)
	{
		ShaderVar *attr=p->attributes+ka;
		if((*attr->pvar=glGetAttribLocation(p->program, attr->name))==-1)
		{
			crash(__FILE__, __LINE__, 1, "%s: attribute %s == -1", p->name, attr->name);
			return 0;
		}
	}
	for(ku=0;ku<p->n_unif;++ku)
	{
		ShaderVar *unif=p->uniforms+ku;
		if((*unif->pvar=glGetUniformLocation(p->program, unif->name))==-1)
		{
			crash(__FILE__, __LINE__, 1, "%s: uniform %s == -1", p->name, unif->name);
			return 0;
		}
	}
	return 1;
}

void set_region_immediate(int x1, int x2, int y1, int y2)
{
	rx0=x1, ry0=y1, rdx=x2-x1, rdy=y2-y1;
	glViewport(rx0, wndh-y2, rdx, rdy);

	/*frac bias*/
#if 1
	SN_x1=2.f/rdx, SN_x0=(float)(-2*rx0-(rdx-1))/rdx;
	SN_y1=-2.f/rdy, SN_y0=(float)(rdy-1+2*ry0)/rdy;
#endif

	/*frac bias*/
#if 0
	SN_x1=2.f/rdx, SN_x0=-rx0*SN_x1-(float)(rdx-1)/rdx;
	SN_y1=-2.f/rdy, SN_y0=(float)(rdy-1)/rdy-ry0*SN_y1;
#endif
	/*X vernier lines*/
#if 0
	SN_x1=2.f/(rdx-(rdx!=0)), SN_x0=-rx0*SN_x1-1;
	SN_y1=-2.f/(rdy-(rdy!=0)), SN_y0=1-ry0*SN_y1;
#endif

	/*X old equation*/
#if 0
	SN_x1=2.f/rdx, SN_x0=-rx0*SN_x1-1;
	SN_y1=-2.f/rdy, SN_y0=1-ry0*SN_y1;
#endif

	NS_x1=rdx/2.f, NS_x0=rx0+NS_x1;
	NS_y1=-rdy/2.f, NS_y0=ry0-NS_x1;
}
static void setGLProgram(unsigned program)
{
	if(current_program!=program)
	{
		glUseProgram(current_program=program);
		CHECK_GL(error);
	}
}
static const float inv255=1.f/255;
static void send_color(unsigned location, int color)
{
	__m128 m_255, c;
	ALIGN(16) float comp[4];
	uint8_t *p;

	p=(uint8_t*)&color;
	m_255=_mm_set1_ps(inv255);

	c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));
	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	_mm_store_ps(comp, c);
	glUniform4fv(location, 1, comp);
#if 0
	glUniform4f(location, p[0]*inv255, p[1]*inv255, p[2]*inv255, p[3]*inv255);
#endif
	CHECK_GL(error);
}
#if 0
static void send_color_rgb(unsigned location, int color)
{
	__m128 m_255=_mm_set1_ps(inv255);

	uint8_t *p=(uint8_t*)&color;
	__m128 c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));

	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	{
		ALIGN(16) float comp[4];
		_mm_store_ps(comp, c);
		glUniform3fv(location, 1, comp);
	}
#if 0
	glUniform3f(location, p[0]*inv255, p[1]*inv255, p[2]*inv255);
#endif
	CHECK_GL(error);
}
#endif
static void set_texture_params(int linear)
{
	if(linear)
	{
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	CHECK_GL(error);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	CHECK_GL(error);
	}
	else
	{
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);	CHECK_GL(error);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);	CHECK_GL(error);
	}
}
void send_texture_pot(unsigned gl_texture, const int *rgba, int txw, int txh, int linear)
{
	glBindTexture(GL_TEXTURE_2D, gl_texture);	CHECK_GL(error);
	set_texture_params(linear);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txw, txh, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba);	CHECK_GL(error);
}
void send_texture_pot_int16x1(unsigned gl_texture, const unsigned *texture, int txw, int txh, int linear)
{
	glBindTexture(GL_TEXTURE_2D, gl_texture);	CHECK_GL(error);
	set_texture_params(linear);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, txw, txh, 0, GL_RED, GL_UNSIGNED_INT, texture);	CHECK_GL(error);
}
void send_texture_pot_grey(unsigned gl_texture, const uint8_t *bmp, int txw, int txh, int linear)
{
	glBindTexture(GL_TEXTURE_2D, gl_texture);	CHECK_GL(error);
	set_texture_params(linear);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, txw, txh, 0, GL_RED, GL_UNSIGNED_BYTE, bmp);	CHECK_GL(error);
}
void select_texture(unsigned tx_id, int u_location)
{
	glActiveTexture(GL_TEXTURE0);		CHECK_GL(error);
	glBindTexture(GL_TEXTURE_2D, tx_id);	CHECK_GL(error);/*select texture*/
	glUniform1i(u_location, 0);		CHECK_GL(error);
}

int gl_init(void)
{
	PIXELFORMATDESCRIPTOR pfd=
	{
		sizeof(PIXELFORMATDESCRIPTOR), 1,
		PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER,
		PFD_TYPE_RGBA, 32,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		16,/*depth bits*/
		0, 0, PFD_MAIN_PLANE, 0, 0, 0, 0
	};
	int PixelFormat;

	PixelFormat=ChoosePixelFormat(ghDC, &pfd);
	SetPixelFormat(ghDC, PixelFormat, &pfd);
	hRC=wglCreateContext(ghDC);
	wglMakeCurrent(ghDC, hRC);

	GLversion=(const char*)glGetString(GL_VERSION);

	{
		static const char msg2[]="%s is NULL";
#define GLFUNC(X) (void)((X=(t_##X)(size_t)wglGetProcAddress(#X))!=0||crash(__FILE__, __LINE__, 1, msg2, X));
		GLFUNCLIST
#undef  GLFUNC
	}
		
	glGenBuffers(1, &vertex_buffer);
	glEnable(GL_BLEND);					CHECK_GL(error);/*vast majority of applications need alpha blend*/
	glBlendEquation(GL_FUNC_ADD);				CHECK_GL(error);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	CHECK_GL(error);
#if 0
	glDisable(GL_BLEND);		CHECK_GL();

	glDisable(GL_LINE_SMOOTH);	CHECK_GL();/*does nothing*/
#endif
	glEnable(GL_DEPTH_TEST);

#define	SHADER(NAME)	make_gl_program(&shader_##NAME);
	SHADER_LIST
#undef	SHADER
	
	/*load fonts*/
#if 1
	{
		int iw, ih, bytespp, *rgb, k, size, c;
		uint8_t *bmp;

		iw=0, ih=0, bytespp=0;
		snprintf(g_buf, G_BUF_SIZE, "%sfont.PNG", (char*)HDR_DATA(exedir));
		rgb=(int*)stbi_load(g_buf, &iw, &ih, &bytespp, 4);
		if(!rgb)
		{
			CRASH("Font texture not found.\nPlace a \'font.PNG\' file with the program.\n");
			return 1;
		}
		tdx=(float)(rgb[0]&0xFF);
		tdy=(float)(rgb[1]&0xFF);
		if(fabsf(tdx)<1e-9||fabsf(tdy)<1e-9)
		{
			crash(__FILE__, __LINE__, 1, "Invalid font texture character dimensions: dx=%d, dy=%d", tdx, tdy);
			return 1;
		}
		for(k=0, size=iw*ih;k<size;++k)
		{
			if(rgb[k]&0x00FFFFFF)
				rgb[k]=0xFFFFFFFF;
		}
		for(c=32;c<127;++c)
		{
			QuadCoords *rect=font_coords+c-32;
			int px=(iw>>3)*(c&7), py=(ih>>4)*(c>>3);

			rect->x1=(float)px/iw, rect->x2=(float)(px+tdx)/iw;
			rect->y1=(float)py/ih, rect->y2=(float)(py+tdy)/ih;
		}
		glGenTextures(1, &font_txid);
		send_texture_pot(font_txid, rgb, iw, ih, 0);
		stbi_image_free(rgb);
	
		snprintf(g_buf, G_BUF_SIZE, "%sfont_sdf.PNG", (char*)HDR_DATA(exedir));
		bmp=(uint8_t*)stbi_load(g_buf, &iw, &ih, &bytespp, 1);
		if(bmp)
		{
			SDFTextureHeader header;

			sdf_available=1;
			memcpy(&header, bmp, sizeof(header));
			sdf_dx=header.csize_x;
			sdf_dy=header.csize_y;

			sdf_slope=(float)header.slope;
			for(c=32;c<127;++c)
			{
				QuadCoords *rect=sdf_glyph_coords+c-32;
				int px=header.grid_start_x+header.cell_size_x*(c&7);
				int py=header.grid_start_y+header.cell_size_y*((c>>3)-4);

				rect->x1=(float)px/iw;
				rect->x2=(float)(px+sdf_dx)/iw;
				rect->y1=(float)py/ih;
				rect->y2=(float)(py+sdf_dy)/ih;
			}
			sdf_txh=sdf_dy;
			sdf_dx*=16.f/sdf_dy;
			sdf_dy=16;

			tdx*=0.85f;
			tdy*=0.85f;
			sdf_dx*=0.85f;
			sdf_dy*=0.85f;

			glGenTextures(1, &sdf_atlas_txid);
			send_texture_pot_grey(sdf_atlas_txid, bmp, iw, ih, 1);
			stbi_image_free(bmp);
			sdf_active=1;
			set_text_colors(colors_text);
			sdf_active=0;
		}
		set_text_colors(colors_text);
		toggle_sdftext();
	}
#endif
	return 0;
}
void depth_test(int enable)
{
	if(enable)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);
}

/*immediate drawing functions*/
static float g_fbuf[16]={0};
void draw_line(float x1, float y1, float x2, float y2, int color)
{
	g_fbuf[0]=screen2NDC_x(x1), g_fbuf[1]=screen2NDC_y(y1);
	g_fbuf[2]=screen2NDC_x(x2), g_fbuf[3]=screen2NDC_y(y2);
	setGLProgram(shader_2D.program);	CHECK_GL(error);
	send_color(u_2D_color, color);		CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);	CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, 4*sizeof(float), g_fbuf, GL_STATIC_DRAW);	CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);	CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	CHECK_GL(error);
#endif

#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(GL_LINES, 0, 2);		CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}
void draw_rect(float x1, float x2, float y1, float y2, int color)
{
	float X1, X2, Y1, Y2;

	X1=screen2NDC_x(x1); Y1=screen2NDC_y(y1);
	X2=screen2NDC_x(x2); Y2=screen2NDC_y(y2);
#if 0
	X1=screen2NDC_x_bias(x1); Y1=screen2NDC_y_bias(y1);
	X2=screen2NDC_x_bias(x2); Y2=screen2NDC_y_bias(y2);
#endif
	g_fbuf[0]=X1, g_fbuf[1]=Y1;
	g_fbuf[2]=X2, g_fbuf[3]=Y1;
	g_fbuf[4]=X2, g_fbuf[5]=Y2;
	g_fbuf[6]=X1, g_fbuf[7]=Y2;
	g_fbuf[8]=X1, g_fbuf[9]=Y1;
	setGLProgram(shader_2D.program);	CHECK_GL(error);
	send_color(u_2D_color, color);		CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);	CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, 10*sizeof(float), g_fbuf, GL_STATIC_DRAW);CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);	CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	CHECK_GL(error);
#endif
#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(GL_TRIANGLE_FAN, 0, 5);	CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}
void draw_rect_hollow(float x1, float x2, float y1, float y2, int color)
{
	float
		X1=screen2NDC_x(x1), Y1=screen2NDC_y(y1),
		X2=screen2NDC_x(x2), Y2=screen2NDC_y(y2);
	g_fbuf[0]=X1, g_fbuf[1]=Y1;
	g_fbuf[2]=X1, g_fbuf[3]=Y2;
	g_fbuf[4]=X2, g_fbuf[5]=Y2;
	g_fbuf[6]=X2, g_fbuf[7]=Y1;
	setGLProgram(shader_2D.program);	CHECK_GL(error);
	send_color(u_2D_color, color);		CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);	CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, 8*sizeof(float), g_fbuf, GL_STATIC_DRAW);	CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);	CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	CHECK_GL(error);
#endif
#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(GL_LINE_LOOP, 0, 4);	CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}
void draw_triangle(float x1, float y1, float x2, float y2, float x3, float y3, int color)
{
	float
		X1=screen2NDC_x(x1), Y1=screen2NDC_y(y1),
		X2=screen2NDC_x(x2), Y2=screen2NDC_y(y2),
		X3=screen2NDC_x(x3), Y3=screen2NDC_y(y3);
	g_fbuf[0]=X1, g_fbuf[1]=Y1;
	g_fbuf[2]=X2, g_fbuf[3]=Y2;
	g_fbuf[4]=X3, g_fbuf[5]=Y3;
	setGLProgram(shader_2D.program);	CHECK_GL(error);
	send_color(u_2D_color, color);		CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);	CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float[6]), g_fbuf, GL_STATIC_DRAW);CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);	CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	CHECK_GL(error);
#endif
#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(GL_TRIANGLES, 0, 3);	CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}

void draw_rect_enqueue(ArrayHandle *vertices, float x1, float x2, float y1, float y2)/*primitive GL_TRIANGLES*/
{
	float X1, X2, Y1, Y2, *vptr;

	X1=screen2NDC_x(x1); Y1=screen2NDC_y(y1);
	X2=screen2NDC_x(x2); Y2=screen2NDC_y(y2);
#if 0
	X1=screen2NDC_x_bias(x1); Y1=screen2NDC_y_bias(y1),
	X2=screen2NDC_x_bias(x2); Y2=screen2NDC_y_bias(y2);
#endif
	if(!*vertices)
	{
		ARRAY_ALLOC(float[2], *vertices, 0, 6, 0, 0);
		vptr=(float*)HDR_DATA(vertices[0]);
	}
	else
		vptr=(float*)ARRAY_APPEND(*vertices, 0, 6, 1, 0);
	*vptr++=X1; *vptr++=Y1;
	*vptr++=X1; *vptr++=Y2;
	*vptr++=X2; *vptr++=Y2;
	*vptr++=X2; *vptr++=Y2;
	*vptr++=X2; *vptr++=Y1;
	*vptr++=X1; *vptr++=Y1;
}
void draw_curve_enqueue(ArrayHandle *vertices, float x, float y)/*primitive GL_LINE_STRIP*/
{
	float *vptr;
	if(!*vertices)
	{
		ARRAY_ALLOC(float[2], *vertices, 0, 1, 0, 0);
		vptr=(float*)HDR_DATA(vertices[0]);
	}
	else
		vptr=(float*)ARRAY_APPEND(*vertices, 0, 1, 1, 0);
	*vptr++=screen2NDC_x(x); *vptr++=screen2NDC_y(y);
}
void draw_2d_flush(ArrayHandle vertices, int color, unsigned primitive)
{
	if(!vertices||!vertices->count)
		return;

	setGLProgram(shader_2D.program);		CHECK_GL(error);
	send_color(u_2D_color, color);			CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);		CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);	CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, (int)(vertices->count*vertices->esize), (float*)HDR_DATA(vertices), GL_STATIC_DRAW);	CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);	CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	CHECK_GL(error);
#endif
#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(primitive, 0, (int)vertices->count);	CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);		CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
	vertices->count=0;
}

static ArrayHandle vrtx=0;
static void vrtx_resize(int vcount, int floatspervertex)
{
	int nfloats=vcount*floatspervertex;
	if(!vrtx)
		ARRAY_ALLOC(float, vrtx, 0, nfloats, 0, 0);
	else if(nfloats>(int)vrtx->count)
		ARRAY_APPEND(vrtx, 0, nfloats-vrtx->count, 1, 0);
}
void draw_ellipse(float x1, float x2, float y1, float y2, int color)
{
	float ya, yb, x0, y0, rx, ry, ygain, *vptr;
	int line_count, nlines, kl;

	if(y1<y2)
		ya=y1, yb=y2;
	else
		ya=y2, yb=y1;
	line_count=(int)ceil(yb)-(int)floor(ya);
	vrtx_resize(line_count*2, 2);
	x0=(x1+x2)*0.5f, y0=(y1+y2)*0.5f, rx=fabsf(x2-x0), ry=yb-y0;
	nlines=0;
	ygain=2.f/line_count;
	vptr=(float*)HDR_DATA(vrtx);
	for(kl=0;kl<line_count;++kl)/*pixel-perfect ellipse (no anti-aliasing) drawn as horizontal lines*/
	{
		/*ellipse equation: sq[(x-x0)/rx] + sq[(y-y0)/ry] = 1	->	x0 +- rx*sqrt(1 - sq[(y-y0)/ry])*/
		float
			y=kl*ygain-1,
			x=1-y*y;
		if(x<0)
			continue;
		x=rx*(float)sqrt((double)x);
		y=y*ry+y0;

		y=screen2NDC_y(y);
		*vptr++=screen2NDC_x(x0-x  ), *vptr++=y;
		*vptr++=screen2NDC_x(x0+x+1), *vptr++=y;
		++nlines;
	}
	setGLProgram(shader_2D.program);	CHECK_GL(error);
	send_color(u_2D_color, color);		CHECK_GL(error);

	glEnableVertexAttribArray(a_2D_coords);	CHECK_GL(error);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);			CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, nlines*sizeof(float[4]), (float*)HDR_DATA(vrtx), GL_STATIC_DRAW);	CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, 0);CHECK_GL(error);
#if 0
	glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
	glVertexAttribPointer(a_2D_coords, 2, GL_FLOAT, GL_FALSE, 0, vrtx);	CHECK_GL(error);
#endif
#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glDrawArrays(GL_LINES, 0, nlines*2);	CHECK_GL(error);
	glDisableVertexAttribArray(a_2D_coords);CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}

/*text API*/
int toggle_sdftext(void)
{
	if(sdf_available)
	{
		sdf_active=!sdf_active;
		{
			float temp;
			SWAPVAR(tdx, sdf_dx, temp);
			SWAPVAR(tdy, sdf_dy, temp);
		}
		if(sdf_active)
			font_zoom_min=0.25, font_zoom_max=64;
		else
		{
			if(font_zoom<1)
				font_zoom=1;
			else
				font_zoom=(float)(1<<(31^LZCNT32((int)font_zoom)));
			font_zoom_min=1, font_zoom_max=32;
		}
	}
	return sdf_available;
}
int set_text_color(int color_txt)
{
	if(sdf_active)
	{
		setGLProgram(shader_sdftext.program);
		send_color(u_sdftext_txtColor, color_txt);
	}
	else
	{
		setGLProgram(shader_text.program);
		send_color(u_text_txtColor, color_txt);
	}
	{
		int *comp=(int*)&colors_text;
		int prev=comp[0];
		comp[0]=color_txt;
		return prev;
	}
}
int set_bk_color(int color_bk)
{
	if(sdf_active)
	{
		setGLProgram(shader_sdftext.program);
		send_color(u_sdftext_bkColor, color_bk);
	}
	else
	{
		setGLProgram(shader_text.program);
		send_color(u_text_bkColor, color_bk);
	}
	{
		int *comp=(int*)&colors_text;
		int prev=comp[1];
		comp[1]=color_bk;
		return prev;
	}
}
uint64_t set_text_colors(uint64_t colors)/*0xBKBKBKBK_TXTXTXTX*/
{
	int *comp=(int*)&colors_text;
	{
		uint64_t temp;
		SWAPVAR(colors_text, colors, temp);
	}
	if(sdf_active)
	{
		setGLProgram(shader_sdftext.program);
		send_color(u_sdftext_txtColor, comp[0]);
		send_color(u_sdftext_bkColor, comp[1]);
	}
	else
	{
		setGLProgram(shader_text.program);
		send_color(u_text_txtColor, comp[0]);
		send_color(u_text_bkColor, comp[1]);
	}
	return colors;
}
float print_line_enqueue(ArrayHandle *vertices, float tab_origin, float x, float y, float zoom, const char *msg, int msg_length, int req_cols, int *ret_idx, int *ret_cols)
{
	QuadCoords *txc, *atlas;
	float
		width, height,
		rect[4]={0},/*{x1, y1, x2, y2}*/
		CX1, CX0;
	int tab_origin_cols, cursor_cols, advance_cols, k;

	if(msg_length<1)
		return 0;
	txc=0;
	atlas=sdf_active?sdf_glyph_coords:font_coords;
#if 0
	if(sdf_active)
		zoom*=16.f/sdf_txh;
#endif
	width=tdx*zoom;
	height=tdy*zoom;
	tab_origin_cols=(int)(tab_origin/width);
#if 0
	printable_count=0;
#endif
	cursor_cols=ret_cols?*ret_cols:0;
	if(y+height<ry0||y>=ry0+rdy)/*off-screen optimization*/
		return 0;/*no need to do anything, this line is outside the screen*/
#if 0
		return idx2col(msg, msg_length, (int)(tab_origin/width))*width;
#endif
	CX1=2.f/rdx;
	CX0=CX1*(x-rx0)-1;/*delta optimization*/
	rect[1]=1-(y       -ry0)*2.f/rdy;/*y1*/
	rect[3]=1-(y+height-ry0)*2.f/rdy;/*y2*/

	if(!*vertices)
		ARRAY_ALLOC(float[4], *vertices, 0, 0, (size_t)msg_length<<2, 0);/*vx, vy, txx, txy		x4 vertices/char*/
	else
		ARRAY_APPEND(*vertices, 0, 0, 1, (size_t)msg_length<<2);
	/*vrtx_resize(msg_length<<2, 4);*/

	k=ret_idx?*ret_idx:0;
	if(req_cols<0||cursor_cols<req_cols)
	{
		CX1*=width;
		for(;k<msg_length;++k)
		{
			char c=msg[k];
			if(c>=' '&&c<='~')
				advance_cols=1;
			else if(c=='\t')
			{
				MODVAR(advance_cols, cursor_cols-tab_origin_cols, tab_count);
				advance_cols=tab_count-advance_cols;
				c=' ';
				/*advance_cols=tab_count-mod(cursor_cols-tab_origin_cols, tab_count), c=' ';*/
			}
			else
				advance_cols=0;
			if(advance_cols)
			{
				if(x+(cursor_cols+advance_cols)*width>=rx0&&x+cursor_cols*width<rx0+rdx)/*off-screen optimization*/
				{
					rect[0]=CX1*cursor_cols+CX0;/*xn1*/
					cursor_cols+=advance_cols;
					rect[2]=CX1*cursor_cols+CX0;/*xn2*/
#if 0
					rect[0]=(x+msg_width-rx0)*2.f/rdx-1;/*xn1*/
					rect[1]=1-(y-ry0)*2.f/rdy;/*yn1*/
					rect[2]=(x+msg_width+width-rx0)*2.f/rdx-1;/*xn2*/
					rect[3]=1-(y+height-ry0)*2.f/rdy;/*yn2*/
#endif
#if 0
					toNDC_nobias(float(x+msg_width		), float(y		), rect[0], rect[1]);
					toNDC_nobias(float(x+msg_width+width	), float(y+height	), rect[2], rect[3]);/*y2<y1*/
#endif

					txc=atlas+c-32;
					{
						float *vptr=(float*)ARRAY_APPEND(*vertices, 0, 4, 1, 0);
						*vptr++=rect[0], *vptr++=rect[1], *vptr++=txc->x1, *vptr++=txc->y1;/*top left		*/
						*vptr++=rect[0], *vptr++=rect[3], *vptr++=txc->x1, *vptr++=txc->y2;/*bottom left	*/
						*vptr++=rect[2], *vptr++=rect[3], *vptr++=txc->x2, *vptr++=txc->y2;/*bottom right	*/
						*vptr++=rect[2], *vptr++=rect[1], *vptr++=txc->x2, *vptr++=txc->y1;/*top right		*/
					}
#if 0
					++printable_count;
#endif
				}
				else
					cursor_cols+=advance_cols;
				if(req_cols>=0&&cursor_cols>=req_cols)
				{
					++k;
					break;
				}
			}
		}
	}
	if(ret_idx)
		*ret_idx=k;
	if(ret_cols)
		*ret_cols=cursor_cols;
	return cursor_cols*width;
}
void print_line_flush(ArrayHandle vertices, float zoom)
{
	if(vertices&&vertices->count)
	{
#ifndef NO_3D
		glDisable(GL_DEPTH_TEST);
#endif
		if(sdf_active)
		{
			setGLProgram(shader_sdftext.program);
			glUniform1f(u_sdftext_zoom, zoom*16.f/(sdf_txh*sdf_slope));
			select_texture(sdf_atlas_txid, u_sdftext_atlas);

			glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
			glBufferData(GL_ARRAY_BUFFER, (int)(vertices->count*vertices->esize), (float*)HDR_DATA(vertices), GL_STATIC_DRAW);CHECK_GL(error);
			glVertexAttribPointer(a_sdftext_coords, 4, GL_FLOAT, GL_TRUE, 0, 0);	CHECK_GL(error);
#if 0
			glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
			glVertexAttribPointer(a_sdftext_coords, 4, GL_FLOAT, GL_TRUE, 0, vptr);	CHECK_GL(error);
#endif
			glEnableVertexAttribArray(a_sdftext_coords);	CHECK_GL(error);
			glDrawArrays(GL_QUADS, 0, (int)vertices->count);CHECK_GL(error);
			glDisableVertexAttribArray(a_sdftext_coords);	CHECK_GL(error);
		}
		else
		{
			setGLProgram(shader_text.program);
			select_texture(font_txid, u_text_atlas);

			glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
			glBufferData(GL_ARRAY_BUFFER, (int)(vertices->count*vertices->esize), (float*)HDR_DATA(vertices), GL_STATIC_DRAW);CHECK_GL(error);/*set vertices & texcoords*/
			glVertexAttribPointer(a_text_coords, 4, GL_FLOAT, GL_TRUE, 0, 0);	CHECK_GL(error);
#if 0
			glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
			glVertexAttribPointer(a_text_coords, 4, GL_FLOAT, GL_TRUE, 0, vptr);	CHECK_GL(error);
#endif
			glEnableVertexAttribArray(a_text_coords);		CHECK_GL(error);
			glDrawArrays(GL_QUADS, 0, (int)vertices->count);	CHECK_GL(error);/*draw the quads: 4 vertices per character quad*/
			glDisableVertexAttribArray(a_text_coords);		CHECK_GL(error);
		}
#ifndef NO_3D
		glEnable(GL_DEPTH_TEST);
#endif
		vertices->count=0;
	}
}
float GUIPrint_enqueue(ArrayHandle *vertices, float tab_origin, float x, float y, float zoom, const char *format, ...)
{
	int len;
	va_list args;

	va_start(args, format);
	len=vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);
	return print_line_enqueue(vertices, tab_origin, x, y, zoom, g_buf, len, -1, 0, 0);
}

float print_line_immediate(float tab_origin, float x, float y, float zoom, const char *msg, int msg_length, int req_cols, int *ret_idx, int *ret_cols)
{
	QuadCoords *txc, *atlas;
	float
		width, height,
		rect[4]={0},/*{x1, y1, x2, y2}*/
		CX1, CX0,
		*vptr;
	int tab_origin_cols, printable_count, cursor_cols, advance_cols, k;

	if(msg_length<1)
		return 0;
	txc=0;
	atlas=sdf_active?sdf_glyph_coords:font_coords;
#if 0
	if(sdf_active)
		zoom*=16.f/sdf_txh;
#endif
	width=tdx*zoom;
	height=tdy*zoom;
	tab_origin_cols=(int)(tab_origin/width);
	printable_count=0;
	cursor_cols=ret_cols?*ret_cols:0;
	if(y+height<ry0||y>=ry0+rdy)/*off-screen optimization*/
		return 0;/*no need to do anything, this line is outside the screen*/
#if 0
		return idx2col(msg, msg_length, (int)(tab_origin/width))*width;
#endif
	CX1=2.f/rdx, CX0=CX1*(x-rx0)-1;/*delta optimization*/
	rect[1]=1-(y       -ry0)*2.f/rdy;/*y1*/
	rect[3]=1-(y+height-ry0)*2.f/rdy;/*y2*/
	vrtx_resize(msg_length<<2, 4);/*vx, vy, txx, txy		x4 vertices/char*/
	k=ret_idx?*ret_idx:0;
	vptr=(float*)HDR_DATA(vrtx);
	if(req_cols<0||cursor_cols<req_cols)
	{
		CX1*=width;
		for(;k<msg_length;++k)
		{
			char c=msg[k];
			if(c>=' '&&c<='~')
				advance_cols=1;
			else if(c=='\t')
			{
				MODVAR(advance_cols, cursor_cols-tab_origin_cols, tab_count);
				advance_cols=tab_count-advance_cols;
				c=' ';
#if 0
				advance_cols=tab_count-mod(cursor_cols-tab_origin_cols, tab_count), c=' ';
#endif
			}
			else
				advance_cols=0;
			if(advance_cols)
			{
				if(x+(cursor_cols+advance_cols)*width>=rx0&&x+cursor_cols*width<rx0+rdx)/*off-screen optimization*/
				{
					rect[0]=CX1*cursor_cols+CX0;/*xn1*/
					cursor_cols+=advance_cols;
					rect[2]=CX1*cursor_cols+CX0;/*xn2*/
#if 0
					rect[0]=(x+msg_width-rx0)*2.f/rdx-1;/*xn1*/
					rect[1]=1-(y-ry0)*2.f/rdy;/*yn1*/
					rect[2]=(x+msg_width+width-rx0)*2.f/rdx-1;/*xn2*/
					rect[3]=1-(y+height-ry0)*2.f/rdy;/*yn2*/

					toNDC_nobias(float(x+msg_width		), float(y			), rect[0], rect[1]);
					toNDC_nobias(float(x+msg_width+width	), float(y+height	), rect[2], rect[3]);/*y2<y1*/
#endif
					txc=atlas+c-32;
					*vptr++=rect[0], *vptr++=rect[1],	*vptr++=txc->x1, *vptr++=txc->y1;/*top left	*/
					*vptr++=rect[0], *vptr++=rect[3],	*vptr++=txc->x1, *vptr++=txc->y2;/*bottom left	*/
					*vptr++=rect[2], *vptr++=rect[3],	*vptr++=txc->x2, *vptr++=txc->y2;/*bottom right	*/
					*vptr++=rect[2], *vptr++=rect[1],	*vptr++=txc->x2, *vptr++=txc->y1;/*top right	*/

					++printable_count;
				}
				else
					cursor_cols+=advance_cols;
				if(req_cols>=0&&cursor_cols>=req_cols)
				{
					++k;
					break;
				}
			}
		}
		if(printable_count)
		{
#ifndef NO_3D
			glDisable(GL_DEPTH_TEST);
#endif
			if(sdf_active)
			{
				setGLProgram(shader_sdftext.program);
				glUniform1f(u_sdftext_zoom, zoom*16.f/(sdf_txh*sdf_slope));
				select_texture(sdf_atlas_txid, u_sdftext_atlas);

				glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
				glBufferData(GL_ARRAY_BUFFER, printable_count<<6, (float*)HDR_DATA(vrtx), GL_STATIC_DRAW);CHECK_GL(error);
				glVertexAttribPointer(a_sdftext_coords, 4, GL_FLOAT, GL_TRUE, 0, 0);	CHECK_GL(error);
#if 0
				glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
				glVertexAttribPointer(a_sdftext_coords, 4, GL_FLOAT, GL_TRUE, 0, vptr);	CHECK_GL(error);
#endif
				glEnableVertexAttribArray(a_sdftext_coords);	CHECK_GL(error);
				glDrawArrays(GL_QUADS, 0, printable_count<<2);	CHECK_GL(error);
				glDisableVertexAttribArray(a_sdftext_coords);	CHECK_GL(error);
			}
			else
			{
				setGLProgram(shader_text.program);
				select_texture(font_txid, u_text_atlas);

				glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);				CHECK_GL(error);
				glBufferData(GL_ARRAY_BUFFER, printable_count<<6, (float*)HDR_DATA(vrtx), GL_STATIC_DRAW);CHECK_GL(error);/*set vertices & texcoords*/
				glVertexAttribPointer(a_text_coords, 4, GL_FLOAT, GL_TRUE, 0, 0);	CHECK_GL(error);
#if 0
				glBindBuffer(GL_ARRAY_BUFFER, 0);					CHECK_GL(error);
				glVertexAttribPointer(a_text_coords, 4, GL_FLOAT, GL_TRUE, 0, vptr);	CHECK_GL(error);
#endif
				glEnableVertexAttribArray(a_text_coords);	CHECK_GL(error);
				glDrawArrays(GL_QUADS, 0, printable_count<<2);	CHECK_GL(error);/*draw the quads: 4 vertices per character quad*/
				glDisableVertexAttribArray(a_text_coords);	CHECK_GL(error);
			}
#ifndef NO_3D
			glEnable(GL_DEPTH_TEST);
#endif
		}
	}
	if(ret_idx)
		*ret_idx=k;
	if(ret_cols)
		*ret_cols=cursor_cols;
	return cursor_cols*width;
}
float GUIPrint(float tab_origin, float x, float y, float zoom, const char *format, ...)
{
	int len;
	va_list args;

	va_start(args, format);
	len=vsnprintf(g_buf, G_BUF_SIZE, format, args);
	va_end(args);
	return print_line_immediate(tab_origin, x, y, zoom, g_buf, len, -1, 0, 0);
}
int g_printed=0;
float GUIPrint_append(float tab_origin, float x, float y, float zoom, int show, const char *format, ...)
{
	float width=0;

	if(format)
	{
		va_list args;
		va_start(args, format);
		g_printed+=vsnprintf(g_buf+g_printed, G_BUF_SIZE-g_printed, format, args);
		va_end(args);
	}
	if(show)
	{
		width=print_line_immediate(tab_origin, x, y, zoom, g_buf, g_printed, -1, 0, 0);
		g_printed=0;
	}
	return width;
}

void display_texture(int x1, int x2, int y1, int y2, unsigned txid, float alpha, float tx1, float tx2, float ty1, float ty2)
{
	float _2_w, _2_h, rect[4], _vrtx[16];
	int j;

	_2_w=2.f/wndw, _2_h=2.f/wndh;
	j=0;
	rect[j++]=x1*_2_w-1;
	rect[j++]=1-y1*_2_h,
	rect[j++]=x2*_2_w-1;
	rect[j++]=1-y2*_2_h;/*y2<y1*/
	j=0;
	_vrtx[j++]=rect[0]; _vrtx[j++]=rect[1];  _vrtx[j++]=tx1; _vrtx[j++]=ty1;/*top left	*/
	_vrtx[j++]=rect[0]; _vrtx[j++]=rect[3];  _vrtx[j++]=tx1; _vrtx[j++]=ty2;/*bottom left	*/
	_vrtx[j++]=rect[2]; _vrtx[j++]=rect[3];  _vrtx[j++]=tx2; _vrtx[j++]=ty2;/*bottom right	*/
	_vrtx[j++]=rect[2]; _vrtx[j++]=rect[1];  _vrtx[j++]=tx2; _vrtx[j++]=ty1;/*top right	*/
	setGLProgram(shader_texture.program);	CHECK_GL(error);
	glUniform1f(u_texture_alpha, alpha);	CHECK_GL(error);

	select_texture(txid, u_texture_texture);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);	CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float[16]), _vrtx, GL_STATIC_DRAW);CHECK_GL(error);/*send vertices & texcoords*/

	glVertexAttribPointer(a_texture_coords, 4, GL_FLOAT, GL_FALSE, sizeof(float[4]), (void*)0);CHECK_GL(error);/*select vertices & texcoord*/

#ifndef NO_3D
	glDisable(GL_DEPTH_TEST);
#endif
	glEnableVertexAttribArray(a_texture_coords);	CHECK_GL(error);
	glDrawArrays(GL_QUADS, 0, 4);			CHECK_GL(error);/*draw the quad*/
	glDisableVertexAttribArray(a_texture_coords);	CHECK_GL(error);
#ifndef NO_3D
	glEnable(GL_DEPTH_TEST);
#endif
}
void display_texture_i(int x1, int x2, int y1, int y2, int *rgb, int txw, int txh, float tx1, float tx2, float ty1, float ty2, float alpha, int linear)
{
	static unsigned tx_id=0;
	if(rgb)
	{
#if 0
#	define NPOT_ATIX2300_FIX
#endif
		int *rgb2, w2, h2;
#ifdef NPOT_ATIX2300_FIX
		int expand;
		int logw=floor_log2(txw), logh=floor_log2(txh);

		if(expand=(!GLversion||GLversion[0]<'3')&&(txw>1<<logw||txh>1<<logh))
		{
			int size, ky;

			w2=txw>1<<logw?1<<(logw+1):txw;
			h2=txh>1<<logh?1<<(logh+1):txh;
			int size=w2*h2;
			rgb2=(int*)malloc((size_t)size<<2);
			memset(rgb2, 0, (size_t)size<<2);
			for(ky=0;ky<txh;++ky)
				memcpy(rgb2+w2*ky, rgb+txw*ky, (size_t)txw<<2);
			float nw=(float)txw/w2, nh=(float)txh/h2;
			vrtx[ 2]=0,		vrtx[ 3]=0;
			vrtx[ 6]=0,		vrtx[ 7]=nh;
			vrtx[10]=nw,	vrtx[11]=nh;
			vrtx[14]=nw,	vrtx[15]=0;
		}
		else
#endif
			rgb2=rgb, w2=txw, h2=txh;

		if(!tx_id)/*generate texture id once*/
		{
			glGenTextures(1, &tx_id);
			CHECK_GL(error);
		}
		send_texture_pot(tx_id, rgb2, w2, h2, linear);
#if 0
		glBindTexture(GL_TEXTURE_2D, tx_id);					CHECK_GL(error);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);	CHECK_GL(error);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);	CHECK_GL(error);
		glTexImage2D(GL_TEXTURE_2D,			0, GL_RGBA, w2, h2, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb2);	CHECK_GL(error);/*send bitmap to GPU*/
#endif

		display_texture(x1, x2, y1, y2, tx_id, alpha, tx1, tx2, ty1, ty2);

#ifdef NPOT_ATIX2300_FIX
		if(expand)
			free(rgb2);
#endif
	}
}

/*3D*/
void mat4_lookAt(float *dst, const float *cam, const float *center, const float *up)
{
	float f[3], t1, u[3], s[3];

	vec3_sub(f, center, cam);
	vec3_normalize(f, f, t1);
	vec3_normalize(u, up, t1);
	vec3_cross(s, f, u);
	vec3_normalize(s, s, t1);
	vec3_cross(u, s, f);
	dst[0]=s[0], dst[4]=s[1], dst[8]=s[2], dst[12]=-vec3_dot(s, cam);
	dst[1]=u[0], dst[5]=u[1], dst[9]=u[2], dst[13]=-vec3_dot(u, cam);
	dst[2]=f[0], dst[6]=f[1], dst[10]=f[2], dst[14]=-vec3_dot(f, cam);
	dst[3]=0, dst[7]=0, dst[11]=0, dst[15]=1;
}
void mat4_FPSView(float *dst, const float *campos, float yaw, float pitch)
{
	float _cam[3], cos_p, sin_p, cos_y, sin_y, xaxis[3], yaxis[3], zaxis[3];
	int j;

#if 0
	cos_p=(float)cos(pitch);
	sin_p=(float)sin(pitch);
	cos_y=(float)cos(yaw);
	sin_y=(float)sin(yaw);
	j=0; _cam[j++]=campos[0];	_cam[j++]campos[1];	_cam[j++]=campos[2];
	j=0; xaxis[j++]=cos_y;		xaxis[j++]=0;		xaxis[j++]=-sin_y;
	j=0; yaxis[j++]=sin_y*sin_p;	yaxis[j++]=cos_p;	yaxis[j++]=cos_y*sin_p;
	j=0; zaxis[j++]=sin_y*cos_p;	zaxis[j++]=-sin_p;	zaxis[j++]=cos_p*cos_y;

	j=0;
	dst[j++]=xaxis[0], dst[j++]=yaxis[0], dst[j++]=zaxis[0], dst[j++]=0;
	dst[j++]=xaxis[1], dst[j++]=yaxis[1], dst[j++]=zaxis[1], dst[j++]=0;
	dst[j++]=xaxis[2], dst[j++]=yaxis[2], dst[j++]=zaxis[2], dst[j++]=0;
	dst[j++]=-vec3_dot(xaxis, _cam), dst[j++]=-vec3_dot(yaxis, _cam), dst[j++]=-vec3_dot(zaxis, _cam), dst[j++]=1;
#endif

#if 1
	cos_p=(float)cos(pitch);
	sin_p=(float)sin(pitch);
	cos_y=(float)cos(yaw);
	sin_y=(float)sin(yaw);
	j=0; _cam[j++]=-campos[1];	_cam[j++]=+campos[2];	_cam[j++]=-campos[0];
	j=0; xaxis[j++]=cos_y;		xaxis[j++]=0;		xaxis[j++]=-sin_y;
	j=0; yaxis[j++]=sin_y*sin_p;	yaxis[j++]=cos_p;	yaxis[j++]=cos_y*sin_p;
	j=0; zaxis[j++]=sin_y*cos_p;	zaxis[j++]=-sin_p;	zaxis[j++]=cos_p*cos_y;

	j=0;
	dst[j++]=xaxis[2], dst[j++]=yaxis[2], dst[j++]=zaxis[2], dst[j++]=0;
	dst[j++]=xaxis[0], dst[j++]=yaxis[0], dst[j++]=zaxis[0], dst[j++]=0;
	dst[j++]=xaxis[1], dst[j++]=yaxis[1], dst[j++]=zaxis[1], dst[j++]=0;
	dst[j++]=-vec3_dot(xaxis, _cam), dst[j++]=-vec3_dot(yaxis, _cam), dst[j++]=-vec3_dot(zaxis, _cam), dst[j++]=1;
#endif
}
void mat4_perspective(float *dst, float tanfov, float w_by_h, float znear, float zfar)
{
	int j;

	j=0;
	dst[j++]=1/tanfov,	dst[j++]=0,		dst[j++]=0,				dst[j++]=0;
	dst[j++]=0,		dst[j++]=w_by_h/tanfov,	dst[j++]=0,				dst[j++]=0;
	dst[j++]=0,		dst[j++]=0,		dst[j++]=-(zfar+znear)/(zfar-znear),	dst[j++]=-1;
	dst[j++]=0,		dst[j++]=0,		dst[j++]=-2*zfar*znear/(zfar-znear),	dst[j++]=0;
}
static void GetTransformInverseNoScale(float *dst, const float *src)/*Requires this matrix to be transform matrix, NoScale version requires this matrix be of scale 1*/
{
	/*
	https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html
	*/
#define MakeShuffleMask(x,y,z,w)           (x | (y<<2) | (z<<4) | (w<<6))

/* vec(0, 1, 2, 3) -> (vec[x], vec[y], vec[z], vec[w])*/
#define VecSwizzleMask(vec, mask)          _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(vec), mask))
#define VecSwizzle(vec, x, y, z, w)        VecSwizzleMask(vec, MakeShuffleMask(x,y,z,w))
#define VecSwizzle1(vec, x)                VecSwizzleMask(vec, MakeShuffleMask(x,x,x,x))
/*special swizzle*/
#define VecSwizzle_0022(vec)               _mm_moveldup_ps(vec)
#define VecSwizzle_1133(vec)               _mm_movehdup_ps(vec)

/*return (vec1[x], vec1[y], vec2[z], vec2[w])*/
#define VecShuffle(vec1, vec2, x,y,z,w)    _mm_shuffle_ps(vec1, vec2, MakeShuffleMask(x,y,z,w))
/*special shuffle*/
#define VecShuffle_0101(vec1, vec2)        _mm_movelh_ps(vec1, vec2)
#define VecShuffle_2323(vec1, vec2)        _mm_movehl_ps(vec2, vec1)
	__m128 inM[4], ret[4];
	
	inM[0]=_mm_loadu_ps(src+ 0);
	inM[1]=_mm_loadu_ps(src+ 4);
	inM[2]=_mm_loadu_ps(src+ 8);
	inM[3]=_mm_loadu_ps(src+12);

	/*transpose 3x3, we know m03 = m13 = m23 = 0*/
	{
		__m128 t0=VecShuffle_0101(inM[0], inM[1]); /* 00, 01, 10, 11*/
		__m128 t1=VecShuffle_2323(inM[0], inM[1]); /* 02, 03, 12, 13*/
		ret[0]=VecShuffle(t0, inM[2], 0,2,0,3);    /* 00, 10, 20, 23(=0)*/
		ret[1]=VecShuffle(t0, inM[2], 1,3,1,3);    /* 01, 11, 21, 23(=0)*/
		ret[2]=VecShuffle(t1, inM[2], 0,2,2,3);    /* 02, 12, 22, 23(=0)*/

		/*last line*/
		ret[3]=                   _mm_mul_ps(ret[0], VecSwizzle1(inM[3], 0));
		ret[3]=_mm_add_ps(ret[3], _mm_mul_ps(ret[1], VecSwizzle1(inM[3], 1)));
		ret[3]=_mm_add_ps(ret[3], _mm_mul_ps(ret[2], VecSwizzle1(inM[3], 2)));
		ret[3]=_mm_sub_ps(_mm_setr_ps(0.f, 0.f, 0.f, 1.f), ret[3]);
	}

	_mm_storeu_ps(dst+ 0, ret[0]);
	_mm_storeu_ps(dst+ 4, ret[1]);
	_mm_storeu_ps(dst+ 8, ret[2]);
	_mm_storeu_ps(dst+12, ret[3]);
#undef MakeShuffleMask
#undef VecSwizzleMask
#undef VecSwizzle
#undef VecSwizzle1
#undef VecSwizzle_0022
#undef VecSwizzle_1133
#undef VecShuffle
#undef VecShuffle_0101
#undef VecShuffle_2323
}
void mat4_normalmat3(float *dst, float *m4)/*inverse transpose of top left 3x3 submatrix*/
{
	float temp[16];

	GetTransformInverseNoScale(temp, m4);
	vec3_copy(dst, temp);
	vec3_copy(dst+3, temp+4);
	vec3_copy(dst+6, temp+8);
}

void draw_3D_triangles(Camera const *cam, unsigned vbuf, size_t offset, int nvertices, unsigned txid)
{
	float mView[16], mProj[16], mVP[16];
	__m128 temp1;

	setGLProgram(shader_3D.program);

	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);
	
	glUniformMatrix4fv(u_3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	select_texture(txid, u_3D_texture);
	
	glBindBuffer(GL_ARRAY_BUFFER, vbuf);		CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_vertex);		CHECK_GL(error);
	glVertexAttribPointer    (a_3D_vertex, 3, GL_FLOAT, GL_FALSE, sizeof(float[5]), (void*)offset);	CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_texcoord);	CHECK_GL(error);
	glVertexAttribPointer    (a_3D_texcoord, 2, GL_FLOAT, GL_FALSE, sizeof(float[5]), (void*)(offset+sizeof(float[3])));	CHECK_GL(error);

	glDrawArrays(GL_TRIANGLES, 0, nvertices);	CHECK_GL(error);
	glDisableVertexAttribArray(a_3D_vertex);	CHECK_GL(error);
	glDisableVertexAttribArray(a_3D_texcoord);	CHECK_GL(error);
}

void gpubuf_send_VNT(GPUModel *dst, const float *VVVNNNTT, int n_floats, const int *indices, int n_ints)
{
	dst->n_elements=n_ints, dst->stride=8*sizeof(float);
	dst->vertices_start=0;
	dst->normals_start=3*sizeof(float);
	dst->txcoord_start=6*sizeof(float);
	glGenBuffers(2, &dst->VBO);
	glBindBuffer(GL_ARRAY_BUFFER, dst->VBO);						CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, n_floats*sizeof(float), VVVNNNTT, GL_STATIC_DRAW);	CHECK_GL(error);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, dst->EBO);					CHECK_GL(error);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_ints*sizeof(int), indices, GL_STATIC_DRAW);	CHECK_GL(error);
}
void draw_L3D(Camera const *cam, GPUModel const *model, const float *modelpos, const float *lightpos, int lightcolor)
{
	float
		mView[16], mProj[16], mVP[16], mMVP_Model[32], mNormal[9], sceneInfo[9],
		*mMVP, *mModel;
	__m128 temp1;
	uint8_t *ptr;

	setGLProgram(shader_L3D.program);

	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);

	mMVP=mMVP_Model;
	mModel=mMVP_Model+16;
	mat4_identity(mModel, 1);
	mat4_translate(mModel, modelpos, temp1);
	mat4_normalmat3(mNormal, mModel);
	mat4_mul_mat4(mMVP, mVP, mModel, temp1);

	ptr=(uint8_t*)&lightcolor;
	memcpy(sceneInfo, lightpos, 3*sizeof(float));
	sceneInfo[3]=ptr[0]*inv255;
	sceneInfo[4]=ptr[1]*inv255;
	sceneInfo[5]=ptr[2]*inv255;
	memcpy(sceneInfo+6, &cam->x, 3*sizeof(float));
	
	glUniformMatrix4fv(u_L3D_matVP_Model, 2, GL_FALSE, mMVP_Model);	CHECK_GL(error);
	glUniformMatrix3fv(u_L3D_matNormal, 1, GL_FALSE, mNormal);	CHECK_GL(error);
	glUniform3fv(u_L3D_sceneInfo, 3, sceneInfo);			CHECK_GL(error);
#if 0
	glBindVertexArray(s_VAO);					CHECK_GL(error);
#endif

	select_texture(model->txid, u_L3D_texture);
	
	glBindBuffer(GL_ARRAY_BUFFER, model->VBO);		CHECK_GL(error);
	glEnableVertexAttribArray(a_L3D_vertex);		CHECK_GL(error);
	glVertexAttribPointer(a_L3D_vertex, 3, GL_FLOAT, GL_FALSE, model->stride, (void*)(size_t)model->vertices_start);	CHECK_GL(error);
	glEnableVertexAttribArray(a_L3D_normal);		CHECK_GL(error);
	glVertexAttribPointer(a_L3D_normal, 3, GL_FLOAT, GL_FALSE, model->stride, (void*)(size_t)model->normals_start);		CHECK_GL(error);
	glEnableVertexAttribArray(a_L3D_texcoord);		CHECK_GL(error);
	glVertexAttribPointer(a_L3D_texcoord, 2, GL_FLOAT, GL_FALSE, model->stride, (void*)(size_t)model->txcoord_start);	CHECK_GL(error);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model->EBO);	CHECK_GL(error);

	glDrawElements(GL_TRIANGLES, model->n_elements, GL_UNSIGNED_INT, 0);	CHECK_GL(error);
	glDisableVertexAttribArray(a_L3D_vertex);		CHECK_GL(error);
	glDisableVertexAttribArray(a_L3D_normal);		CHECK_GL(error);
	glDisableVertexAttribArray(a_L3D_texcoord);		CHECK_GL(error);
}

void draw_3d_line(Camera const *cam, const float *w1, const float *w2, int color)
{
	static unsigned txid=0;
	__m128 temp1;
	int bitmap[4], j;
	float
		mView[16], mProj[16], mVP[16],
		*vptr;

	j=0;
	bitmap[j++]=color;
	bitmap[j++]=color;
	bitmap[j++]=color;
	bitmap[j++]=color;

	/*prepare texture*/
	if(!txid)
		glGenTextures(1, &txid);
	send_texture_pot(txid, bitmap, 2, 2, 0);
	
	/*prepare coords*/
	vrtx_resize(2, 5);
	vptr=(float*)HDR_DATA(vrtx);
	vptr[0]=-w1[0];
	vptr[1]=-w1[1];
	vptr[2]=w1[2];
#if 0
	memcpy(vptr, w1, sizeof(float[3]));
#endif
	vptr[3]=0;
	vptr[4]=0;
	
	vptr[5]=-w2[0];
	vptr[6]=-w2[1];
	vptr[7]=w2[2];
#if 0
	memcpy(vptr+5, w2, sizeof(float[3]));
#endif
	vptr[8]=0;
	vptr[9]=0;
	
	/*prepare matrix*/
	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);

	setGLProgram(shader_3D.program);

	select_texture(txid, u_3D_texture);
	
	glUniformMatrix4fv(u_3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);		CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, 2*5*sizeof(float), vptr, GL_STATIC_DRAW);					CHECK_GL(error);/*send vertices & texcoords*/
	glVertexAttribPointer(a_3D_vertex, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)0);			CHECK_GL(error);/*select vertices & texcoords*/
	glVertexAttribPointer(a_3D_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(3*sizeof(float)));	CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_vertex);			CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);

	glDrawArrays(GL_LINES, 0, 2);				CHECK_GL(error);

	glDisableVertexAttribArray(a_3D_vertex);		CHECK_GL(error);
	glDisableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);
}

void draw_3d_line_enqueue(ArrayHandle *vertices, float x, float y, float z)
{
	float *vptr;
	if(!*vertices)
	{
		ARRAY_ALLOC(float[5], *vertices, 0, 1, 0, 0);
		vptr=(float*)HDR_DATA(vertices[0]);
	}
	else
		vptr=(float*)ARRAY_APPEND(*vertices, 0, 1, 1, 0);

	vptr[0]=-x;
	vptr[1]=-y;
	vptr[2]=z;
	vptr[3]=0;
	vptr[4]=0;
}
void draw_3d_flush(ArrayHandle vertices, Camera const *cam, int *texture, int tw, int th, int linear, unsigned primitive)
{
	static unsigned txid=0;

	float mView[16], mProj[16], mVP[16];
	__m128 temp1;

	/*prepare texture*/
	if(!txid)
		glGenTextures(1, &txid);
	send_texture_pot(txid, texture, tw, th, linear);
	
	/*prepare matrix*/
	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);

	setGLProgram(shader_3D.program);

	select_texture(txid, u_3D_texture);
	
	glUniformMatrix4fv(u_3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);		CHECK_GL(error);
	glBufferData(GL_ARRAY_BUFFER, (int)(vertices->count*vertices->esize), (float*)HDR_DATA(vertices), GL_STATIC_DRAW);	CHECK_GL(error);/*send vertices & texcoords*/
	glVertexAttribPointer(a_3D_vertex, 3, GL_FLOAT, GL_FALSE, (int)vertices->esize, (void*)0);				CHECK_GL(error);/*select vertices & texcoords*/
	glVertexAttribPointer(a_3D_texcoord, 2, GL_FLOAT, GL_FALSE, (int)vertices->esize, (void*)(3*sizeof(float)));		CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_vertex);			CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);

	glDrawArrays(primitive, 0, (int)vertices->count);	CHECK_GL(error);

	glDisableVertexAttribArray(a_3D_vertex);		CHECK_GL(error);
	glDisableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);

	vertices->count=0;
}
void draw_3d_wireframe_gpu(Camera const *cam, unsigned gpubuf, int nvertices, int color, unsigned primitive)
{
	static unsigned txid=0;
	int texture[4], j;
	float mView[16], mProj[16], mVP[16];
	__m128 temp1;

	j=0;
	texture[j++]=color;
	texture[j++]=color;
	texture[j++]=color;
	texture[j++]=color;
	if(!txid)
		glGenTextures(1, &txid);
	send_texture_pot(txid, texture, 2, 2, 0);

	/*prepare matrix*/
	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);

	setGLProgram(shader_3D.program);

	select_texture(txid, u_3D_texture);
	
	glUniformMatrix4fv(u_3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	
	glBindBuffer(GL_ARRAY_BUFFER, gpubuf);			CHECK_GL(error);
	glVertexAttribPointer(a_3D_vertex, 3, GL_FLOAT, GL_FALSE, (int)sizeof(float[5]), (void*)0);			CHECK_GL(error);/*select vertices & texcoords*/
	glVertexAttribPointer(a_3D_texcoord, 2, GL_FLOAT, GL_FALSE, (int)sizeof(float[5]), (void*)(3*sizeof(float)));	CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_vertex);			CHECK_GL(error);
	glEnableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);

	glDrawArrays(primitive, 0, nvertices);			CHECK_GL(error);

	glDisableVertexAttribArray(a_3D_vertex);		CHECK_GL(error);
	glDisableVertexAttribArray(a_3D_texcoord);		CHECK_GL(error);
}

void draw_contour3d_rect(Camera const *cam, unsigned vbuf, int nvertices, unsigned txid, float alpha)
{
	float mView[16], mProj[16], mVP[16];
	__m128 temp1;

	setGLProgram(shader_contour3D.program);

	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);
	
	glUniform1f(u_contour3D_alpha, alpha);

	glUniformMatrix4fv(u_contour3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	select_texture(txid, u_contour3D_texture);
	
	glBindBuffer(GL_ARRAY_BUFFER, vbuf);				CHECK_GL(error);
	glEnableVertexAttribArray(a_contour3D_vertex);			CHECK_GL(error);
	glVertexAttribPointer    (a_contour3D_vertex, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)0);			CHECK_GL(error);
	glEnableVertexAttribArray(a_contour3D_texcoord);		CHECK_GL(error);
	glVertexAttribPointer    (a_contour3D_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(3*sizeof(float)));	CHECK_GL(error);

	glDrawArrays(GL_TRIANGLES, 0, nvertices);			CHECK_GL(error);
	glDisableVertexAttribArray(a_contour3D_vertex);			CHECK_GL(error);
	glDisableVertexAttribArray(a_contour3D_texcoord);		CHECK_GL(error);
}

typedef struct ConVertStruct
{
	float x, y, z, tx, ty;
} ConVert;
void draw_contour3d(Camera const *cam, float x1, float x2, float y1, float y2, float z1, float z2, unsigned *textures, int count, float alpha)
{
	static unsigned buffer=0;
	float mView[16], mProj[16], mVP[16];
	__m128 temp1;
	ConVert v[6]={0}, *pv;
	int kz, k2;

	
	pv=v;
	pv->x=-x1; pv->y=-y1, pv->z=0, pv->tx=0, pv->ty=0; ++pv;
	pv->x=-x2; pv->y=-y1, pv->z=0, pv->tx=1, pv->ty=0; ++pv;
	pv->x=-x2; pv->y=-y2, pv->z=0, pv->tx=1, pv->ty=1; ++pv;
	pv->x=-x2; pv->y=-y2, pv->z=0, pv->tx=1, pv->ty=1; ++pv;
	pv->x=-x1; pv->y=-y2, pv->z=0, pv->tx=0, pv->ty=1; ++pv;
	pv->x=-x1; pv->y=-y1, pv->z=0, pv->tx=0, pv->ty=0; ++pv;
	if(!buffer)
		glGenBuffers(1, &buffer);
	setGLProgram(shader_contour3D.program);
	mat4_FPSView(mView, &cam->x, cam->ax, cam->ay);
	mat4_perspective(mProj, cam->tanfov, (float)(wndw)/wndh, 0.1f, 1000.f);
	mat4_mul_mat4(mVP, mProj, mView, temp1);
	glUniformMatrix4fv(u_contour3D_matrix, 1, GL_FALSE, mVP);	CHECK_GL(error);
	
	glUniform1f(u_contour3D_alpha, alpha);

	glBindBuffer(GL_ARRAY_BUFFER, buffer);			CHECK_GL(error);
	glEnableVertexAttribArray(a_contour3D_vertex);		CHECK_GL(error);
	glVertexAttribPointer    (a_contour3D_vertex, 3, GL_FLOAT, GL_FALSE, sizeof(float[5]), (void*)0);	CHECK_GL(error);/*use offsetof*/
	glEnableVertexAttribArray(a_contour3D_texcoord);	CHECK_GL(error);
	glVertexAttribPointer    (a_contour3D_texcoord, 2, GL_FLOAT, GL_FALSE, sizeof(float[5]), (void*)(3*sizeof(float)));	CHECK_GL(error);
	for(kz=0;kz<count;++kz)
	{
		float z=z1+(z2-z1)*kz/count;
		for(k2=0;k2<6;++k2)
			v[k2].z=z;

		select_texture(textures[kz], u_contour3D_texture);

		glBufferData(GL_ARRAY_BUFFER, sizeof(v), v, GL_STATIC_DRAW);	CHECK_GL(error);
		glDrawArrays(GL_TRIANGLES, 0, 6);				CHECK_GL(error);
	}
	glDisableVertexAttribArray(a_contour3D_vertex);		CHECK_GL(error);
	glDisableVertexAttribArray(a_contour3D_texcoord);	CHECK_GL(error);
}
