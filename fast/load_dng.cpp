#include"fast.h"
#include<stdint.h>
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wpointer-arith"
#pragma GCC diagnostic ignored "-Winit-self"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wformat=2"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#ifdef __cplusplus
//C++ warnings
#pragma GCC diagnostic ignored "-Wshift-negative-value"
//#pragma GCC diagnostic ignored "-Wcalloc-transposed-args"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#else
//C warnings:
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#pragma GCC diagnostic ignored "-Wdeclaration-after-statement"
#pragma GCC diagnostic ignored "-Wc++-compat"
#endif
#endif
#define TINY_DNG_LOADER_IMPLEMENTATION
#include"tiny_dng_loader.h"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
static const char file[]=__FILE__;

extern "C" int load_dng(const char *fn, Image *dst)
{
	std::vector<tinydng::FieldInfo> custom_fields;
	std::vector<tinydng::DNGImage> images;
	std::string warnings, errors;
	bool success=tinydng::LoadDNG(fn, custom_fields, &images, &warnings, &errors);
	if(!success||!images.size())
		return 0;
	int klargest=0;
	size_t maxres=0;
	for(int ki=0;ki<(int)images.size();++ki)
	{
		auto &im=images[ki];
		size_t res=(size_t)im.width*im.height;
		if(!maxres||maxres<res)
			maxres=res, klargest=ki;
	}
	auto &src=images[klargest];
	memset(dst, 0, sizeof(*dst));
	dst->iw=src.width;
	dst->ih=src.height;
	dst->nch=src.samples_per_pixel==4?-4:1;
	dst->depth=src.bits_per_sample;
	size_t bufsize=sizeof(short)*dst->iw*dst->ih;
	dst->data=(short*)malloc(bufsize);
	if(!dst->data)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(dst->data, 0, bufsize);
	const unsigned char *data=src.data.data();
	int nlevels=1<<src.bits_per_sample, half=nlevels>>1;
	for(int ky=0, idx=0, kb2=0;ky<dst->ih;++ky)
	{
		for(int kx=0;kx<dst->iw;++kx, ++idx)
		{
			unsigned val=0;
			for(int kb=0;kb<dst->depth;++kb, ++kb2)
			{
				int bit=data[kb2>>3]>>(kb2&7)&1;
				val|=bit<<kb;
			}
			dst->data[idx]=(short)(val-half);
		}
	}
	return 1;
}