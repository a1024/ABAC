#include"fast.h"
#include<stdint.h>
#define TINY_DNG_LOADER_IMPLEMENTATION
#include"tiny_dng_loader.h"
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