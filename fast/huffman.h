#pragma once
#ifndef INC_HUFFMAN_H
#define INC_HUFFMAN_H
#include"util.h"
#include<string.h>

#define HUFF_DEC_BITS 4
typedef struct HuffmanCoderStruct
{
	unsigned long long cache;
	unsigned nbits;//enc: number of free bits in cache [0 ~ 64], dec: number of unread bits in cache [0 ~ 32]
	int is_enc;
	union
	{
		struct
		{
			const unsigned char *srcptr, *srcend, *srcstart;
		};
#ifdef EC_USE_STATIC_BUFFER
		struct
		{
			unsigned char *dstptr, *dstend, *dststart;
		};
#else
		DList *list;
#endif
	};
} HuffmanCoder;
typedef struct HuffmanNodeStruct
{
	struct HuffmanNodeStruct *ch[1<<HUFF_DEC_BITS];
	int sym, nbits;
} HuffmanNode;
typedef struct HuffInternalNodeStruct
{
	struct HuffInternalNodeStruct *left, *right;//0 and 1
	size_t symbol, freq;
	BitstringHandle code;
} HuffInternalNode, *HuffInternalNodeHandle;
HuffInternalNodeHandle huff_tree_gen(unsigned *hist, int nlevels, int exclude_zps);
void huff_tree_free(HuffInternalNodeHandle *root);
void huff_alphabet_gen(HuffInternalNodeHandle root, BitstringHandle *alphabet);//alphabet has nlevels BitstringHandles
void huff_alphabet_free(BitstringHandle *alphabet, int nlevels);
void huff_alphabet_print(unsigned *hist, BitstringHandle const *alphabet, int nlevels);

INLINE void huff_enc_init(HuffmanCoder *ec, DList *list)
{
	ec->cache=0;
	ec->nbits=sizeof(ec->cache)<<3;
	ec->is_enc=1;
	ec->list=list;
}
INLINE int huff_enc(HuffmanCoder *ec, int sym, BitstringHandle const *alphabet)
{
	//BitstringHandle code=alphabet[sym];
	//int bitidx=0;
	//if(code->bitCount>ec->nbits)
	//{
	//	//TODO
	//}
	return 0;
}

HuffmanNode* huff_dectree_gen(HuffInternalNodeHandle root);
void huff_dectree_free(HuffmanNode **root);
void huff_dec_init(HuffmanCoder *ec, const unsigned char *start, const unsigned char *end);
INLINE int huff_dec(HuffmanCoder *ec, HuffmanNode *root)
{
	HuffmanNode *node;
	for(;;)
	{
		if(ec->nbits<HUFF_DEC_BITS)
		{
			if(ec->srcptr+4>ec->srcend)
			{
				LOG_ERROR2("Huffman buffer overflow");
				return 0;
			}
			ec->nbits+=32;
			ec->cache<<=32;
			memcpy(&ec->cache, ec->srcptr, 4);
			ec->srcptr+=4;
		}
		node=root->ch[ec->cache>>(ec->nbits-HUFF_DEC_BITS)];
		if(!node)
			break;
		ec->nbits-=root->nbits;
		ec->cache&=(1LL<<ec->nbits)-1;
		root=node;
	}
	return root->sym;
}

void huff_test(Image const *image);


#endif