//huffman.cpp - Huffman Coder implementation
//Copyright (C) 2021  Ayman Wagih Mohsen
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

#include"huffman.h"
#include<stdio.h>
#include<stdarg.h>
#include<string.h>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>

#ifdef __linux__
#define	__rdtsc	__builtin_ia32_rdtsc
#else
#include<intrin.h>
#endif

//Huffman coder
void 			print(const char *format, ...)
{
	va_list args;
	va_start(args,format);
	vprintf(format, args);
	va_end(args);
}
void 			print_flush()
{
	printf("\n");
}
void			print_bin(const byte *data, int bytesize)
{
	for(int kb=0;kb<(bytesize<<3);++kb)
	{
		print("%d", data[kb>>3]>>(kb&7)&1);
		if((kb&7)==7)
			print("-");
	}
}
struct			AlphabetComparator
{
	vector_bool const *alphabet;
	explicit AlphabetComparator(vector_bool const *alphabet):alphabet(alphabet){}
	bool operator()(int idx1, int idx2)const
	{
		auto &s1=alphabet[idx1], &s2=alphabet[idx2];
		for(int kb=0;kb<s1.bitSize&&kb<s2.bitSize;++kb)
		{
			int bit1=s1.get(kb), bit2=s2.get(kb);
			if(bit1!=bit2)
				return bit1<bit2;
		}
		return s1.bitSize<s2.bitSize;//shortest symbol first
		//return s1.bitSize>s2.bitSize;//longest symbol first
	}
};
static void		sort_alphabet(vector_bool const *alphabet, int nLevels, std::vector<int> &idx)
{
	idx.resize(nLevels);
	for(int k=0;k<nLevels;++k)
		idx[k]=k;
	std::sort(idx.begin(), idx.end(), AlphabetComparator(alphabet));
}
static void		print_alphabet(vector_bool const *alphabet, const int *histogram, int nlevels, int symbols_to_compress, const int *sort_idx)
{
	print("symbol");
	if(histogram)
		print(", freq, %%");
	print_flush();
	for(int k=0;k<nlevels;++k)//print alphabet
	{
		int symbol=sort_idx?sort_idx[k]:k;
		if(histogram&&!histogram[symbol])
			continue;
		print("%4d ", symbol);
		if(histogram)
		{
			print("%6d ", histogram[symbol]);
			print("%2d ", histogram[symbol]*100/symbols_to_compress);
		}
		auto &code=alphabet[symbol];
		for(int k2=0, k2End=code.bitSize;k2<k2End;++k2)
			print("%c", char('0'+code.get(k2)));
		print_flush();
	}
}
void			print_histogram(int *histogram, int nlevels, int scanned_count, int *sort_idx)
{
	int histmax=0;
	for(int k=0;k<nlevels;++k)
		if(histmax<histogram[k])
			histmax=histogram[k];
	const int consolechars=79-15-5*(sort_idx!=0);
	if(!histmax)
		return;

	if(sort_idx)
		print("idx, ");
	print("symbol, freq, %%");
	print_flush();
	for(int k=0;k<nlevels;++k)//print histogram
	{
		int symbol=sort_idx?sort_idx[k]:k;
		if(!histogram[symbol])
			continue;
		if(sort_idx)
			print("%4d ", k);
		print("%4d %6d %2d ", symbol, histogram[symbol], histogram[symbol]*100/scanned_count);
		for(int kr=0, count=histogram[symbol]*consolechars/histmax;kr<count;++kr)
			print("*");
		print_flush();
	}
}


struct			Node
{
	int branch[2];
	unsigned short value;
	int freq;
};
static std::vector<Node> tree;//root is at the end of array
static int		nLevels;
static int		make_node(int symbol, int freq, int left, int right)//https://gist.github.com/pwxcoo/72d7d3c5c3698371c21e486722f9b34b
{
	int idx=tree.size();
	tree.push_back(Node());
	auto &n=*tree.rbegin();
	n.value=symbol, n.freq=freq;
	n.branch[0]=left, n.branch[1]=right;
	return idx;
}
struct			compare_nodes
{
	bool operator()(int idx1, int idx2)
	{
		if(tree[idx1].freq==tree[idx2].freq)
			return idx1<idx2;
		return tree[idx1].freq>tree[idx2].freq;
	}
};

static void		print_tree()
{
	for(int k=tree.size()-1;k>=0;--k)
	{
		auto &node=tree[k];
		if(!node.freq)
			continue;
		print("[%d] 0:%d,1:%d, freq=%d, val=%d", k, node.branch[0], node.branch[1], node.freq, node.value);
		print_flush();
	}
}
static void		build_tree(const int *histogram, int nLevels)
{
	::nLevels=nLevels;
	tree.clear();
	tree.reserve(nLevels);
	std::priority_queue<int, std::vector<int>, compare_nodes> pq((compare_nodes()));
	for(int k=0;k<nLevels;++k)
		pq.push(make_node(k, histogram[k], -1, -1));
	while(pq.size()>1)//build Huffman tree
	{
		int left=pq.top();	pq.pop();
		int right=pq.top();	pq.pop();
		pq.push(make_node(0, tree[left].freq+tree[right].freq, left, right));
	}
}
static void		make_alphabet(std::vector<vector_bool> &alphabet)
{
	alphabet.resize(nLevels);
	typedef std::pair<int, vector_bool> TraverseInfo;
	std::stack<TraverseInfo> s;
	s.push(TraverseInfo(tree.size()-1, vector_bool()));
	vector_bool left, right;
	while(s.size())//depth-first
	{
		auto &info=s.top();
		int idx=info.first;
		if(idx==-1)
		{
			s.pop();
			continue;
		}
		auto &r2=tree[idx];
		if(r2.branch[0]==-1&&r2.branch[1]==-1)
		{
			alphabet[r2.value]=std::move(info.second);
			s.pop();
			continue;
		}
		left=std::move(info.second);
		right=left;
		s.pop();
		if(r2.branch[1]!=-1)
		{
			right.push_back(true);
			s.push(TraverseInfo(r2.branch[1], std::move(right)));
		}
		if(r2.branch[0]!=-1)
		{
			left.push_back(false);
			s.push(TraverseInfo(r2.branch[0], std::move(left)));
		}
	}
}
void			calculate_histogram(const short *image, int size, int *histogram, int nLevels)
{
	memset(histogram, 0, nLevels*sizeof(int));
	for(int k=0;k<size;++k)
		++histogram[image[k]];
}
typedef std::vector<vector_bool> Alphabet;
int				huff_rough_size_estimate(Alphabet const &alphabet, int imsize)
{
	int length=0;
	for(int k=0;k<(int)alphabet.size();++k)
		length+=alphabet[k].bitSize;
	length/=nLevels;
	return length*imsize>>5;
}

void			huff_encode(const short *buffer, int imsize, int depth, int *histogram, vector_bool &bits, bool loud)
{
	auto t1=__rdtsc();
	int nlevels=1<<depth;
	calculate_histogram(buffer, imsize, histogram, nlevels);
	auto t_hist=__rdtsc();

	build_tree(histogram, nlevels);
	auto t_tree=__rdtsc();

	std::vector<vector_bool> alphabet;
	make_alphabet(alphabet);
	auto t_alpha=__rdtsc();

	int size_estimated=huff_rough_size_estimate(alphabet, imsize);
	bits.data.clear(), bits.bitSize=0;
	bits.data.reserve(size_estimated);
	for(int k=0;k<imsize;++k)
		bits.push_back(alphabet[buffer[k]]);
	bits.clear_tail();
	auto t2=__rdtsc();

	if(loud)
	{
		int compressed_bytesize=bits.data.size()*sizeof(int);
		printf("Huffman encode:  %lld cycles\n", t2-t1);
		printf("Size: %d -> %d bytes,  ratio: %lf\n", imsize, compressed_bytesize, (double)imsize/compressed_bytesize);
		printf("Histogram:\t%lld\n", t_hist-t1);
		printf("Build tree:\t%lld\n", t_tree-t_hist);
		printf("Alphabet:\t%lld\n", t_alpha-t_tree);
		printf("Pack:\t%lld\n", t2-t_alpha);
	}
}
void			huff_decode(const int *src, long long bitsize, int imsize, int depth, const int *histogram, short *dst, bool loud)
{
	auto t1=__rdtsc();

	int nlevels=1<<depth;
	build_tree(histogram, nlevels);

	int bit_idx=0, kd=0;//naive 1bit decode
	for(;kd<imsize&&bit_idx<bitsize;++kd)
	{
		int prev=tree.size()-1, node=prev;
		while(bit_idx<bitsize&&(tree[node].branch[0]!=-1||tree[node].branch[1]!=-1))
		{
			int ex_idx=bit_idx>>5, in_idx=bit_idx&31;
			int bit=src[ex_idx]>>in_idx&1;
			prev=node;
			node=tree[node].branch[bit];
			++bit_idx;
		}
		dst[kd]=tree[node].value;
	}
	auto t2=__rdtsc();
	if(loud)
	{
		printf("Huffman decode:  %lld cycles\n", t2-t1);
	}
	if(bit_idx!=bitsize)
	{
		print("Decompression error:\n"), print_flush();
		print("\tbit_idx = %d", bit_idx), print_flush();
		print("\tbitsize = %lld", bitsize), print_flush();
		print_flush();
	}
}