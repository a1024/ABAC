#include"lz77.h"
#include<string.h>
#include<stdarg.h>
#ifdef __linux__
#define __rdtsc __builtin_ia32_rdtsc
#define scanf_s scanf
#endif

static void		print_hex(const char *buffer, int size)
{
	for(int k=0;k<size;++k)
		printf("%02X-", buffer[k]&0xFF);
	printf("\n");
}
static void		breakpoint()
{
	printf("BREAKPOINT\n");
	int x=0;
	scanf_s("%d", &x);
}
static void		error(const char *data, int size, int idx, const char *format, ...)
{
	va_list args;
	va_start(args, format);
	vprintf(format, args);
	va_end(args);
	int start=idx-100;
	if(start<0)
		start=0;
	int end=start+200;
	if(end>size)
		end=size;
	printf("Before the error index:\n");
	print_hex(data+start, idx-start);
	printf("After the error index:\n");
	print_hex(data+idx, end-idx);
	breakpoint();
}

int data_count=0, int_count=0;

const int		copyflag=0x80;
//const int		smax=0x80,//0xFE, 0x7F
//				copyflag=0xFF;
inline void		int2_push(std::string &data, int x)
{
	auto p=(char*)&x;
	data.insert(data.end(), p, p+4);
	int_count+=4;
}
inline void		int2_add(unsigned char *data, int k, int value)
{
	int v0=0;
	memcpy(&v0, data+k, sizeof(int));
	v0+=value;
	memcpy(data+k, &v0, sizeof(int));
}
inline int		int2_get(const unsigned char *data, int &ks)
{
	int value=0;
	memcpy(&value, data+ks, sizeof(int));
	ks+=4;
	return value;
}
void			int1_push(std::string &data, int x)
{
	if(x)
	{
		while(x)
			data.push_back(x&0x7F), x>>=7, ++int_count;
		data.back()|=0x80;//set bit 7, int last byte
	}
	else
		data.push_back(0), data.push_back((char)0x80);

	//x=x&0x7F|(x>>7)<<8;//bit 7 should be clear	//little-endian
	//while(x)
	//	data.push_back(x&0xFF), x>>=8, ++int_count;
	//data.push_back(0), ++int_count;

	//auto p=(char*)&x;
	//char p2[]={p[3], p[2], p[1], p[0]};//big-endian
	//data.insert(data.end(), p2, p2+4);
	//int_count+=4;

	//while(x)
	//	data.push_back((char)x%smax), x/=smax, ++int_count;
	//data.push_back(0), ++int_count;

	//while(x>=smax)
	//	data.push_back((char)smax), x-=smax, ++int_count;
	//data.push_back(x), ++int_count;
}
int				int1_get(const unsigned char *data, int &ks)
{
	int value=0;

	int i=0;
	for(;!(data[ks]&0x80);++ks, i+=7)
		value|=data[ks]<<i;
	value|=(data[ks]&0x7F)<<i, ++ks;

	//auto p=(char*)&value;
	//p[3]=data[ks];//big-endian
	//p[2]=data[ks+1];
	//p[1]=data[ks+2];
	//p[0]=data[ks+3];
	//ks+=4;

	//for(;data[ks];++ks)
	//	value*=smax, value+=data[ks];
	//++ks;

	//while((data[ks]&0xFF)==smax)
	//	value+=smax, ++ks;
	//value+=data[ks]&0xFF, ++ks;

	return value;
}
void			lz77_encode(const void *src, int imsize, std::string &data, bool loud)
{
	const int min_backtrack_size=8;
	auto buffer=(const unsigned char*)src;
	auto t1=__rdtsc();
	data_count=0, int_count=0;
	int *table=new int[256];
	memset(table, 0xFF, 256*sizeof(int));
	data.reserve(imsize);
	for(int k=0, append=0, block_start=0;k<imsize;)
	{
		auto &entry=table[buffer[k]];
		if(entry==-1)//copy as it is
		{
			entry=k;
			int start=k;
			++k;
			for(;k<imsize&&table[buffer[k]]==-1;++k)
				table[buffer[k]]=k;
			if(append)
				int2_add((unsigned char*)data.c_str(), block_start, k-start);
			else
			{
				append=1;
				data.push_back((char)copyflag), ++int_count;
				block_start=data.size();
				int2_push(data, k-start);
			}
			data.insert(data.end(), buffer+start, buffer+k);
			data_count+=k-start;
		}
		else//backtrack
		{
			int start=entry, backtrack=k-entry, offset=1;
			++k;
			for(;k<imsize&&buffer[k]==buffer[start+offset];++k, ++offset)
				table[buffer[k]]=k;//replace entries for smallest backtrack

			if(offset>=min_backtrack_size)
			{
				int1_push(data, backtrack);
				int1_push(data, offset);
				append=0;
			}
			else//copy as it is
			{
				if(!append)
				{
					append=1;
					data.push_back((char)copyflag), ++int_count;
					block_start=data.size();
					int2_push(data, offset);
				}
				else
					int2_add((unsigned char*)data.c_str(), block_start, offset);
				data.insert(data.end(), buffer+k-offset, buffer+k);
				data_count+=offset;
			}
		}
		//if(loud&&!(k&0xFF))
		//	printf("\r%6d / %6d = %lf%%, csize: %d, ratio: %lf\t\t", k, imsize, 100.*k/imsize, (int)data.size(), (double)k/data.size());
	}
	auto t2=__rdtsc();
	if(loud)
	{
		printf("\n");
		printf("LZ77 Encode: %lld\n", t2-t1);
		printf("Size: %d -> %d, ratio: %lf\n", imsize, (int)data.size(), (double)imsize/data.size());
		printf("Data usage:\n");
		printf("\tData: %5d\t%lf%%\n", data_count, 100.*data_count/data.size());
		printf("\tInt:  %5d\t%lf%%\n", int_count, 100.*int_count/data.size());
		int preview=data.size();
#if 1
		if(preview>=200)
			preview=200;
		printf("Preview:\n");
		print_hex(data.data(), preview);
#endif
#if 0
		printf("Explanation:\n");
		for(int ks=0, kd=0;ks<imsize;)
		{
			int ks0=ks;
			if((unsigned char)data[ks]==copyflag)
			{
				++ks;
				int size=int2_get((unsigned char*)data.c_str(), ks);
				print_hex(data.c_str()+ks0, (ks-ks0)+size);
				if(size<0||size>imsize)
					error(data.data(), data.size(), ks0, "ERROR: copy size = %d, original size = %d\n", size, imsize);
				//{
				//	printf("ERROR: copy size = %d, original size = %d\n", size, imsize);
				//	breakpoint();
				//}
				printf("ks=%d, kd=%d: copy %d:\n%.*s\n\n", ks, kd, size, size, data.data()+ks);
				ks+=size;
				//for(int count=0;count<size&&ks<preview;++count, ++ks)
				//	printf("%c", data[ks]);
				//printf("\n");
				kd+=size;
			}
			else
			{
				int backtrack=int1_get((unsigned char*)data.c_str(), ks), offset=int1_get((unsigned char*)data.c_str(), ks);
				print_hex(data.c_str()+ks0, ks-ks0);
				if(backtrack<0||backtrack>kd)//
					error(data.data(), data.size(), ks0, "ERROR: backtrack = %d, kdst = %d\n", backtrack, kd);
				//{
				//	printf("ERROR: backtrack = %d, kdst = %d\n", backtrack, kd);
				//	breakpoint();
				//}
				if(offset<0||offset>imsize)//
					error(data.data(), data.size(), ks0, "ERROR: offset = %d, original size = %d\n", offset, imsize);
				//{
				//	printf("ERROR: offset = %d, original size = %d\n", offset, imsize);
				//	breakpoint();
				//}
				printf("ks=%d, kd=%d: back: %d, size: %d\n%.*s\n\n", ks, kd, backtrack, offset, offset, buffer+kd-backtrack);
				kd+=offset;
			}
		}
#endif
	}
}
void			lz77_decode(const void *src, int csize, int imsize, char *buffer, bool loud)
{
	auto data=(const unsigned char*)src;
	auto t1=__rdtsc();
	for(int ks=0, kd=0;ks<csize&&kd<imsize;)
	{
		int ks0=ks;
		if(data[ks]==copyflag)//new block
		{
			++ks;
			int blocksize=int2_get(data, ks);
			if(blocksize<0||blocksize>imsize)//
				error((char*)data, csize, ks0, "ERROR: blocksize = %d, original size = %d\n", blocksize, imsize);
			//{
			//	printf("ERROR: blocksize = %d, original size = %d\n", blocksize, imsize);
			//	breakpoint();
			//}
			memcpy(buffer+kd, data+ks, blocksize);
			ks+=blocksize, kd+=blocksize;
		}
		else//backtrack
		{
			int backtrack=int1_get(data, ks);
			int offset=int1_get(data, ks);
			if(backtrack<0||backtrack>kd)//
				error((char*)data, csize, ks0, "ERROR: backtrack = %d, kdst = %d\n", backtrack, kd);
			//{
			//	printf("ERROR: backtrack = %d, kdst = %d\n", backtrack, kd);
			//	breakpoint();
			//}
			if(offset<0||offset>imsize)//
				error((char*)data, csize, ks0, "ERROR: offset = %d, original size = %d\n", offset, imsize);
			//{
			//	printf("ERROR: offset = %d, original size = %d\n", offset, imsize);
			//	breakpoint();
			//}
			for(int k2=0;k2<offset;++k2, ++kd)
				buffer[kd]=buffer[kd-backtrack];
			//memmove(buffer+kd, buffer+kd-backtrack, offset);
			//kd+=offset;
		}
		//if(loud&&!(ks&0xFF))
		//	printf("\r%6d / %6d = %lf%%\t\t", ks, csize, 100.*ks/csize);
	}
	printf("\n");
	auto t2=__rdtsc();
	if(loud)
	{
		printf("LZ77 Decode: %lld\n", t2-t1);
	}
}
#if 0
const int smax=0xFFFE;
inline void		push_int(std::vector<short> &data, int x)
{
	while(x>=smax)
		data.push_back(smax), x-=smax;
	data.push_back(x);
}
void			lz77_encode(const short *buffer, int imsize, int depth, std::vector<short> &data, bool loud)
{
	auto t1=__rdtsc();
	int nlevels=1<<depth;
	int *table=new int[nlevels];
	memset(table, 0xFF, nlevels*sizeof(int));
	for(int k=0;k<imsize;)
	{
		auto &entry=table[buffer[k]];
		if(entry==-1)
		{
			entry=k;
			int start=k;
			++k;
			for(;k<imsize&&table[buffer[k]]==-1;++k);
			int count=k-start;
			data.push_back(0xFFFF);
			push_int(data, count);
			data.insert(data.end(), buffer+start, buffer+k);
		}
		else
		{
			int start=k, offset=1;
			++k;
			for(;k<imsize&&buffer[k]==buffer[entry+offset];++k, ++offset);
			int backtrack=k-start;

			push_int(data, backtrack);
			push_int(data, offset);
		}
	}
	auto t2=__rdtsc();
	if(loud)
	{
		printf("LZ77 Encode: %lld\n", t2-t1);
	}
}
int				get_int(const short *data, int ks)
{
	int value=0;
	while(*data==smax)
		value+=smax, ++data, ++ks;
	value+=*data, ++ks;
	return value;
}
void			lz77_decode(const short *data, int csize, int imsize, int depth, short *buffer, bool loud)
{
	for(int ks=0, kd=0;ks<csize;++ks)
	{
		if(data[ks]==0xFFFF)//new block
		{
			++ks;
			int blocksize=get_int(data, ks);
			memcpy(buffer+kd, data+ks, blocksize*sizeof(short));
			ks+=blocksize, kd+=blocksize;
		}
		else
		{
			int backtrack=get_int(data, ks);
			int offset=get_int(data, ks);
			memmove(buffer+kd, buffer+kd-backtrack, offset*sizeof(short));
			kd+=offset;
		}
	}
}
#endif