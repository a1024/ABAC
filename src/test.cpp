#include"ac.h"
#include"huffman.h"
#include"lz77.h"

#include"lodepng.h"
#ifndef __linux__
#include<Windows.h>
#endif
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<string>
#include<fstream>

#ifdef __linux__
#define	__rdtsc	__builtin_ia32_rdtsc
#else
#include<intrin.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"

#define	SIZEOF(STATIC_ARRAY)	(sizeof(STATIC_ARRAY)/sizeof(*STATIC_ARRAY))

#ifdef __linux__
typedef unsigned char byte;
#define		set_console_buffer_size(...)
#else
int			set_console_buffer_size(short w, short h)
{
	COORD coords={w, h};
	HANDLE handle=GetStdHandle(STD_OUTPUT_HANDLE);
	int success=SetConsoleScreenBufferSize(handle, coords);
	if(!success)
		printf("Failed to resize console buffer: %d\n\n", GetLastError());
	return success;
}
#endif
bool			open_text(const char *filename, std::string &data)
{
	std::ifstream input(filename);
	if(!input.is_open())
		return false;
	data.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
	input.close();
	return true;
}

void			print_hex(short *buffer, int imsize, int depth)
{
	for(int k=0;k<imsize;++k)
		printf("%0*X-", (depth+3)>>2, buffer[k]);
	printf("\n");
}
void			print_rgba8(unsigned *data, int npixels)
{
	printf("AABBGGRR:\n");
	for(int k=0;k<npixels;++k)
		printf("%08X-", data[k]);
	printf("\n");
}
void			print_bin(int x, int nbits)
{
	//printf("0b");
	for(int k=nbits-1;k>=0;--k)
	{
		printf("%d", x>>k&1);
		if(k&&!(k&3))
			printf("_");
	}
}
void			print_bitplane(const short *buffer, int imsize, int bitplane)
{
	for(int k=0;k<imsize;++k)
		printf("%d", buffer[k]>>bitplane&1);
	printf("\n");
}

#if 0
const int	depth=8,
			nlevels=1<<depth;
int			histogram[nlevels]={};
#endif
int			main(int argc, char **argv)
{
	set_console_buffer_size(120, 4000);
	if(argc<2)
	{
		printf("Pass filename as command argument\n");
		return 1;
	}

	//LZ77
#if 0
	int iw=0, ih=0, nch=0;
	//byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	byte *original_image=stbi_load("20211129 1 confidence.PNG", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Couldn't open image\n");
		return 1;
	}
	auto buffer=(char*)original_image;
	int imsize=iw*ih;//*/

/*	std::string str;
//	open_text("D:/C/Compression Sandbox/Compression Sandbox/ac.cpp", str);
	open_text("D:/C/Compression Sandbox/g2.cpp", str);
	auto buffer=str.c_str();
	int imsize=str.size();//*/

/*	const char buffer[]="01234555555555555555555555555555555555555556789";
	//const char buffer[]="Hello World! Sample Text. What is going on??!";
	//const char buffer[]="0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_0101_";
	printf("Original:\n%s\n", buffer);
	int imsize=sizeof(buffer);//*/

	//const int depth=8;
	char *b2=new char[imsize];
	memset(b2, 0, imsize);

	std::string data;
	lz77_encode(buffer, imsize, data, true);
	lz77_decode(data.data(), data.size(), imsize, b2, true);
	bool error=false;
	for(int k=0;k<imsize;++k)
	{
		if(b2[k]!=buffer[k])
		{
			error=true;
			printf("Error %d: %d -> %d, %c -> %c\n", k, b2[k]&0xFF, buffer[k]&0xFF, b2[k], buffer[k]);
			int start=k-100;
			if(start<0)
				start=0;
			int end=start+200;
			if(end>imsize)
				end=imsize;
			printf("Before:\n%.*s\nAfter:\n%.*s\n", end-start, buffer+start, end-start, b2+start);
			break;
		}
	}
	if(!error)
	{
		int preview=imsize;
		if(preview>200)
			preview=200;
		printf("Decoded:\n%.*s\n", preview, b2);
	}
	
	//STBI_FREE(original_image);
	delete[] b2;
#endif

	//Huffman - text
#if 0
	std::string text;
	open_text("abac.cpp", text);
	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];
	
	vector_bool bits;
	huff_encode(buffer, imsize, 8, bits, true);
/*	auto t1=__rdtsc();
	calculate_histogram(buffer, imsize, histogram, nlevels);
	auto t_hist=__rdtsc();
	build_tree(histogram, nlevels);
	auto t_tree=__rdtsc();
	std::vector<vector_bool> alphabet;
	make_alphabet(alphabet);
	auto t_alpha=__rdtsc();
	int size_estimated=huff_rough_size_estimate(alphabet, imsize);
	vector_bool bits;
	bits.data.reserve(size_estimated);
	for(int k=0;k<imsize;++k)
		bits.push_back(alphabet[buffer[k]]);
	bits.clear_tail();
	auto t2=__rdtsc();

	int compressed_bytesize=bits.data.size()*sizeof(int);
	printf("Huffman:  %lld cycles\n", t2-t1);
	printf("Size: %d -> %d bytes,  ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
	printf("Histogram:\t%lld\n", t_hist-t1);
	printf("Build tree:\t%lld\n", t_tree-t_hist);
	printf("Alphabet:\t%lld\n", t_alpha-t_tree);
	printf("Pack:\t%lld\n", t2-t_alpha);//*/
#endif

	//Huffman - image
#if 0
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Failed to open image\n");
		return 1;
	}
	auto image=(int*)original_image;
	//unsigned iw=0, ih=0;
	//std::vector<unsigned char> vimage;
	//int error=lodepng::decode(vimage, iw, ih, "example.png");
	//if(error)
	//{
	//	printf("%s\n", lodepng_error_text(error));
	//	return 1;
	//}
	//auto image=(int*)vimage.data();

	int imsize=iw*ih;
//	int imsize=800;
	short *buffer=new short[imsize];
	const int depth=8;
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;
	short *b2=new short[imsize];
	STBI_FREE(original_image);


	int histogram[256]={};
	vector_bool data;
	huff_encode(buffer, imsize, 8, histogram, data, true);

	huff_decode(data.data.data(), data.bitSize, imsize, 8, histogram, b2, true);


	delete[] buffer;
	delete[] b2;
#endif

	//AC - image
#if 1
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load(argv[1], &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("20211129 1 confidence.PNG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("2005-12-29 Empire state building 29122005.JPG", &iw, &ih, &nch, 4);
	//byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(original_image)
		printf("Opened \'%s\'\n", argv[1]);
	else
	{
		printf("Failed to open \'%s\'\n", argv[1]);
		return 1;
	}
	auto image=(int*)original_image;
	const int depth=8;
	int imsize=iw*ih;//*/

/*	const char text[]="Sample text";
	const int depth=8;
	int imsize=sizeof(text);
	const char *image=text;//*/

	short *buffer=new short[imsize];
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;//*/
	//stbi_image_free(original_image);

/*	printf("Opening \'%s\'\n", argv[1]);
	std::string text;
	open_text(argv[1], text);
	//open_text("D:/C/ABAC Linux/Makefile", text);
	const int depth=8;
	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];//*/


/*	int imsize=1024;
	const int depth=1;
	short *buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=rand()&1;
		//buffer[k]=!(k&1);
		//buffer[k]=0;//*/


/*	unsigned iw=0, ih=0;
	std::vector<unsigned char> vimage;
	int error=lodepng::decode(vimage, iw, ih, "example.png");
	if(error)
	{
		printf("%s\n", lodepng_error_text(error));
		return 1;
	}
	auto image=(int*)vimage.data();
	int imsize=iw*ih;
//	int imsize=800;
	short *buffer=new short[imsize];
	const int depth=8;
	for(int k=0;k<imsize;++k)//extract red channel
		buffer[k]=image[k]&0xFF;//*/
	//memset(buffer, 0, imsize*sizeof(short));
	//for(int k=0;k<100000;++k)
	//	buffer[rand()%imsize]=rand();
	auto b2=new short[imsize];

	//print_rgba8((unsigned*)image.data(), imsize);
	//print_hex(buffer, imsize, depth);
	
#if 0
	vector_bool bits;
	int histogram[256];
	huff_encode(buffer, imsize, 8, histogram, bits, true);
	huff_decode(bits.data.data(), bits.bitSize, imsize, 8, histogram, b2, true);
#endif
#if 1
	std::string data;
	int sizes[depth]={};
	int conf[depth]={};
//	int probs[depth]={};

	//abac_estimate(buffer, imsize, depth, 2, true);

	abac2_encode(buffer, imsize, depth, data, sizes, conf, true);
	abac2_decode(data.data(), sizes, conf, b2, imsize, depth, true);

	//abac_encode(buffer, imsize, depth, data, sizes, true);
	//abac_decode(data.data(), sizes, b2, imsize, depth, true);

	//abac_encode_sse2(buffer, imsize, depth, data, sizes, true);
	//abac_decode_sse2(data.data(), sizes, b2, imsize, depth, true);

	//abac_encode_avx2(buffer, imsize, depth, data, sizes, true);
	//abac_decode_avx2(data.data(), sizes, b2, imsize, depth, true);

	//ac_encode(buffer, imsize, depth, data, sizes, probs, true);
	//ac_decode(data.data(), sizes, probs, b2, imsize, depth, true);
	
	//ac_debug(buffer, imsize, depth, data, sizes, probs, b2, true);
#endif


	//print_hex(b2, imsize, depth);

	int nerrors=0, kp=0, kb=0;
	int depthmask=(1<<depth)-1;
	//printf("Depthmask: %04X\n", depthmask);
	for(int k=0;k<imsize;++k)
	{
		if((b2[k]^buffer[k])&depthmask)
		{
			if(!nerrors)
			{
				kb=k;
				for(kp=0;kp<depth&&!((b2[k]^buffer[k])>>kp&1);++kp);
			}
			++nerrors;
			printf("%d: %04X != %04X\n", k, (int)buffer[k], (int)b2[k]);
			//if(nerrors>=100)
				break;
		}
	}
	//printf("original, recovered, xor:\n");
	//for(int k=0;k<imsize;++k)
	//{
	//	printf("%3d ", k);
	//	print_bin(buffer[k], depth);
	//	printf(", ");
	//	print_bin(b2[k], depth);
	//	printf(", ");
	//	print_bin(buffer[k]^b2[k], depth);
	//	printf("\n");
	//}

	//for(int k=0;k<depth;++k)//print all bitplanes
	//{
	//	printf("Plane %d:\n", k);
	//	print_bitplane(buffer, imsize, k);//
	//	print_bitplane(b2, imsize, k);//
	//}
	//printf("\n");

	if(nerrors)
	{
		const int width=64;
		printf("Error in bit plane %d, pixel %d:\n", kp, kb);
		int start=kb-(width>>1);
		if(start<0)
			start=0;
		int end=start+width;
		if(end>imsize)
			end=imsize;
		for(int k=start;k<end;++k)
			printf("%d", k%10);
		printf("\n");
		print_bitplane(buffer+start, end-start, kp);//
		print_bitplane(b2+start, end-start, kp);//
	}

	//print_bitplane(buffer, imsize, 0);//
	//print_bitplane(b2, imsize, 0);//

	//print_bitplane(buffer+2800, 100, 7);//
	//print_bitplane(b2+2800, 100, 7);//
	//print_bitplane(buffer+700, 100, 0);//
	//print_bitplane(b2+700, 100, 0);//
	//print_bitplane(buffer+215000, 100, 0);//
	//print_bitplane(b2+215000, 100, 0);//
	
	delete[] b2;
	delete[] buffer;
#endif

	//AC text
#if 0
	const char str[]="111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100000111111";
//	const char str[]="000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="111100000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="11111111111111110100000001";
//	const char str[]="00111111100000011111001000000100010001100111001100001010001000100100101101000001011111011101111111110111111111101111101110011110101111111111111111010000000110011000000000000101100100100000000011100000";
//	const char str[]="011110100000100011101110000011010010110100110100101011000110111101111011000111001110110111001100111";
//	const char str[]="1011110011010110000010110001111000111010111101001010100100011101";
//	const char str[]="0000000100000001000000010000000100000001000000010000000100000001";
	const int imsize=sizeof(str);
	short buffer[imsize]={};
	for(int k=0;k<imsize;++k)
		buffer[k]=str[k]-'0';
	//	buffer[k]=rand()&1;

	//int dmask=0;
	//ac_test_bitplane_differentiation(buffer, imsize, 1, dmask);
	//if(dmask)
	//	ac_differentiate_bitplanes(buffer, imsize, 1, dmask);
	
	const int depth=1;
	std::string data;
	int sizes[depth], probs[depth];
	short b2[imsize]={};

	//for(int k=0;k<imsize-1;++k)
	//{
	//	printf("k=%d\n", k);
	//	memset(b2, 0, sizeof(b2));
	//	ac_debug(buffer+k, imsize-k, 1, data, sizes, probs, b2, true);
	//}
	ac_debug(buffer, imsize, 1, data, sizes, probs, b2, true);
//	ac_encode(buffer, imsize, 1, data, sizes, true);
//	ac_decode(data.c_str(), sizes.data(), b2, imsize, 1, true);

	//if(dmask)
	//{
	//	ac_integrate_bitplanes(buffer, imsize, 8, dmask);
	//	ac_integrate_bitplanes(b2, imsize, 8, dmask);
	//}

	print_bitplane(buffer, imsize, 0);
	print_bitplane(b2, imsize, 0);//*/
#endif

	//AC
#if 0
	const int imsize=100;
	short buffer[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=1;
	buffer[62]&=0;
	std::string data;
	int sizes[imsize], probs[imsize];
	short b2[imsize]={};
	ac_debug(buffer, imsize, 1, data, sizes, probs, b2, true);
	print_bitplane(buffer, imsize, 0);
	print_bitplane(b2, imsize, 0);
#endif

	//AC - text
#if 0
	//std::string text=R"(void			ac_print_summary())";//*/
	
	std::string text;
	for(int k=0;k<100;++k)//100
		text.push_back(0xFF);
	text[62]&=0;
//	text[62]&=0xDF;//62
	
	//std::string text;
	//if(!open_text("Makefile", text))
	//{
	//	printf("Failed to open file\n");
	//	return 1;
	//}

	//text.resize(100);
	//for(int k=0;k<text.size();++k)
	//	text[k]=rand()&1;

	//memset(&text[0], 'A', text.size());

	int imsize=text.size();
	auto buffer=new short[imsize];
	for(int k=0;k<imsize;++k)
		buffer[k]=text[k];
	
//	const int depth=8;
	const int depth=1;

	//printf("Bitplane 2:\n");
	//print_bitplane(buffer, imsize, 2);
	int dmask=0;
	ac_test_bitplane_differentiation(buffer, imsize, depth, dmask);
	dmask=0;//
	if(dmask)
	{
		printf("dmask: ");
		print_bin(dmask, depth);
		printf("\n");
	//	ac_differentiate_bitplanes(buffer, imsize, depth, dmask);
	}//*/
/*	int bitplane=0;
	int p1[depth]={};
	for(int k=0;k<imsize;++k)
		for(int k2=0;k2<depth;++k2)
			p1[k2]+=buffer[k]>>k2&1;
	for(int k=0, latch=0;k<imsize;++k)//differentiator
	{
		if(k<imsize-1)
			latch=buffer[k]^buffer[k+1];
		else
			latch=buffer[k];
		buffer[k]=latch;
	}
	int p2[depth]={};
	for(int k=0;k<imsize;++k)
		for(int k2=0;k2<depth;++k2)
			p2[k2]+=buffer[k]>>k2&1;
		printf("One\'s count: (size=%d)\n", imsize);
	for(int k2=0;k2<depth;++k2)
		printf("%d: %d -> %d (50 + %lf%% -> 50 + %lf%%)\n", k2, p1[k2], p2[k2], 50-100.*p1[k2]/imsize, 50-100.*p2[k2]/imsize);
	print_bitplane(buffer, imsize, bitplane);//*/

	std::string data;
	int sizes[depth]={}, prob[depth]={};
	ac_encode(buffer, imsize, depth, data, sizes, prob, true);

	auto b2=new short[imsize];
	ac_decode(data.c_str(), sizes, prob, b2, imsize, depth, true);
	if(dmask)
	{
	//	ac_integrate_bitplanes(buffer, imsize, depth, dmask);
	//	ac_integrate_bitplanes(b2, imsize, depth, dmask);
	}
	for(int k=0, error_count=0;k<imsize;++k)
	{
		int error=(buffer[k]^b2[k])&((1<<depth)-1);
		if(error)
		{
			printf("Error %d: ", k);
			print_bin(buffer[k], depth);
			printf("!=");
			print_bin(b2[k], depth);
			printf(", xor=");
			print_bin(error, depth);
			printf(",%c!=%c\n", (char)buffer[k], (char)b2[k]);
			//printf("Error %d: %04X!=%04X,%3d!=%3d,%c!=%c\n", k, (int)buffer[k], (int)b2[k], (int)buffer[k], (int)b2[k], (char)buffer[k], (char)b2[k]);
			++error_count;
			if(error_count>=100)
				break;
		}
	}
	int printsize=imsize;
	if(printsize>100)
		printsize=100;
	printf("Original: ");
	print_bitplane(buffer, printsize, 0);//5
	printf("Result:   ");
	print_bitplane(b2, printsize, 0);//5

	//std::vector<Symbol> data;
	//ac_encode(buffer, imsize, depth, data, true);
	//
	//ac_decode(data.data(), buffer, imsize, depth, true);
	std::string str(printsize, '\0');
	for(int k=0;k<(int)str.size();++k)
		str[k]=b2[k];
	printf("Result:\n%.*s\n", printsize, str.c_str());
/*	for(int k=0;k<head;++k)
	{
		printf("text:\t");
		print_bin(text[k]&((1<<depth)-1), depth);
		printf(" %d %c\n", text[k], text[k]);

		printf("dec:\t");
		print_bin(str[k]&((1<<depth)-1), depth);
		printf(" %d %c\n\n", str[k], str[k]);
	}
	//	printf("%02X: %c\n", str[k]&255, str[k]);//*/

	delete[] buffer;
#endif

	//AC v1 - raw PNG
#if 0
	int iw=0, ih=0, nch=0;
	byte *original_image=stbi_load("example.png", &iw, &ih, &nch, 4);
	if(!original_image)
	{
		printf("Failed to open image\n");
		return 1;
	}
	int imsize=iw*ih;
	short *red=new short[imsize];
	for(int k=0;k<imsize;++k)//extract red channel
		red[k]=original_image[k]&0xFF;
	stbi_image_free(original_image);
	
	std::vector<Symbol> data;
	
	//auto t1=__rdtsc();
	ac_encode(red, imsize, 8, data, true);
	//auto t2=__rdtsc();

/*	int compressed_bytesize=(int)data.size()*sizeof(Symbol);
	printf("u%d/u%d:\n", SYMBOL_BITS, SYMBOL_BITS*2);
	printf("Size: %d -> %d, ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
	//printf("Size: %dx%d=%d -> %d bytes\n", iw, ih, imsize, compressed_bytesize);
	//printf("Compresison ratio: %lf\n", (double)imsize/compressed_bytesize);
	printf("Elapsed: %lld cycles\n", t2-t1);//*/
	delete[] red;
#endif

#ifdef _MSC_VER
	int i=0;
	scanf_s("%d", &i);
#endif
	return 0;
}