#include"codec.h"
#include<stdlib.h>
#include<math.h>//abs
static const char file[]=__FILE__;


//	#define ENABLE_FILEGUARD


#include"entropy.h"
int c08_codec(const char *srcfn, const char *dstfn)
{
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec C08 requires both source and destination filenames");
		return 1;
	}
#ifdef ENABLE_FILEGUARD
	ptrdiff_t dstsize=get_filesize(dstfn);
	if(dstsize>=0)
	{
		LOG_ERROR("Destination file already exists");
		return 1;
	}
#endif
	FILE *fsrc=fopen(srcfn, "rb");
	int tag=0;
	size_t nread=fread(&tag, 1, 2, fsrc), nwritten=0;
	if(nread!=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
	}
	else
	{
		nread=fread(&iw, 1, 4, fsrc);
		nread+=fread(&ih, 1, 4, fsrc);
		if(nread!=8)
		{
			LOG_ERROR("Unsupported archive");
			return 1;
		}
	}
	FILE *fdst=fopen(dstfn, "wb");
	AC5 ec;
	nwritten=0;
	if(fwd)
	{
		nwritten+=fwrite("CH", 1, 2, fdst);
		nwritten+=fwrite(&iw, 1, 4, fdst);
		nwritten+=fwrite(&ih, 1, 4, fdst);
		ac5_enc_init(&ec, fdst);
	}
	else
	{
		nwritten+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		ac5_dec_init(&ec, fsrc);
	}
	unsigned short stats0[3][256]={0};
	FILLMEM((unsigned short*)stats0, 0x8000, sizeof(stats0), sizeof(short));
	int psize=(iw+16LL)*sizeof(short[4*4]);//4 padded rows * 4 channels max
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		unsigned char rgb[3]={0};
		int yuv[3]={0};
		for(int kx=0;kx<iw;++kx)
		{
			if(fwd)
			{
				fread(rgb, 1, 3, fsrc);
				yuv[0]=rgb[1]-128;
				yuv[1]=rgb[2]-128;
				yuv[2]=rgb[0]-128;
			}
			for(int kc=0;kc<3;++kc)
			{
				int offset=kc?yuv[0]:0;
				short
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					W	=rows[0][kc-1*4],
					*curr	=rows[0]+kc;
				int pred;
				MEDIAN3_32(pred, N, W, N+W-NW);
				pred+=offset;
				CLAMP2_32(pred, pred, -128, 127);
				int error=0;
				if(fwd)
				{
					*curr=yuv[kc];
					error=yuv[kc]-pred;
					error<<=24;
					error>>=24;
				}
				unsigned short *curr_stats=stats0[kc];
				int bit=0;
				for(int kb=7, tidx=1;kb>=0;--kb)
				{
					int p1=curr_stats[tidx];
					if(fwd)
					{
						bit=error>>kb&1;
						
						while((ec.lo^ec.hi)<0x1000000)
						{
							fputc(ec.lo>>24, ec.f);
							ec.lo<<=8;
							ec.hi=ec.hi<<8|255;
						}
						{
							unsigned mid=ec.lo+(unsigned)((unsigned long long)(ec.hi-ec.lo)*p1>>16);
						//	unsigned mid=ec.lo+((ec.hi-ec.lo)>>AC5_PROB_BITS)*(p1>>(16-AC5_PROB_BITS));
							ec.lo=bit?ec.lo:mid+1;
							ec.hi=bit?mid:ec.hi;
						}
						//ac5_enc_bin(&ec, p1>>(16-AC5_PROB_BITS), bit);
					}
					else
					{
						while((ec.lo^ec.hi)<0x1000000)
						{
							ec.code=ec.code<<8|(fgetc(ec.f)&255);
							ec.lo<<=8;
							ec.hi=ec.hi<<8|255;
						}
						{
							unsigned mid=ec.lo+(unsigned)((unsigned long long)(ec.hi-ec.lo)*p1>>16);
						//	unsigned mid=ec.lo+((ec.hi-ec.lo)>>AC5_PROB_BITS)*(p1>>(16-AC5_PROB_BITS));
							bit=ec.code<=mid;
							ec.lo=bit?ec.lo:mid+1;
							ec.hi=bit?mid:ec.hi;
						}
						//bit=ac5_dec_bin(&ec, p1>>(16-AC5_PROB_BITS));
						error|=bit<<kb;
					}
					p1+=((bit<<16)-p1+(1<<7>>1))>>7;
					curr_stats[tidx]=p1;
					tidx+=tidx+bit;
				}
				if(!fwd)
				{
					error+=pred;
					error<<=24;
					error>>=24;
					*curr=yuv[kc]=error;
				}
				*curr-=offset;
			}
			if(!fwd)
			{
				rgb[0]=yuv[2]+128;
				rgb[1]=yuv[0]+128;
				rgb[2]=yuv[1]+128;
				fwrite(rgb, 1, 3, fdst);
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	if(fwd)
		ac5_enc_flush(&ec);
	(void)nwritten;
	return 0;
}