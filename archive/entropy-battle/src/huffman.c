#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
static const char file[]=__FILE__;

//	#define HUFF_LATE_FREE_TREE
//	#define HUFF_NODEHIST
//	#define HUFF_ALPHAGUARD
//	#define HUFF_PROGRESS
//	#define DEBUG_HUFF_ENC
//	#define DEBUG_HUFF_DEC

const int tag_huff='H'|'U'<<8|'F'<<16|'F'<<24;
typedef struct HuffFileHeaderStruct
{
    int tag, zeros;
    size_t
		u_size,//uncompressed (original) size
		c_size;//compressed size, excluding the header
    size_t hist[256];
    unsigned char data[];
} HuffFileHeader, *HuffFileHandle;

typedef struct HuffNodeStruct
{
	struct HuffNodeStruct *left, *right;//0 and 1
	size_t symbol, freq;
	BitstringHandle code;
} HuffNode, *HuffNodeHandle;

int huffnode_less(const void *left, const void *right)
{
	HuffNode *L=*(HuffNode**)left, *R=*(HuffNode**)right;
	return L->freq>R->freq;//for min heap
}
void huffnode_print(const void *data)
{
	HuffNodeHandle p=*(HuffNodeHandle*)data;

	if(p->symbol!=-1)
		printf("%d:%d", (int)p->freq, (int)p->symbol);
	else
		printf("%d", (int)p->freq);
	if(p->code)
	{
		printf(" code=");
		bitstring_print(p->code);
	}
	printf("\n");
	
	//printf("%d:%d ", p->symbol, p->freq);
	//printf("%c:%d ", (char)p->symbol, p->freq);
	//printf("0x%02X:%d ", p->symbol, p->freq);
	//printf("sym=0x%02X freq=0x%08X left=%p right=%p\n", p->symbol, p->freq, p->left, p->right);
}
void print_tree(HuffNodeHandle root, int depth)
{
	if(root)
	{
		print_tree(root->left, depth+1);

		printf("%3d %*s", depth, depth, "");
		huffnode_print(&root);
		//printf("%3d %*s%d: %d\n", depth, depth, "", root->symbol, root->freq);

		print_tree(root->right, depth+1);
	}
}
void free_tree(HuffNodeHandle *root)
{
	if(*root)
	{
		free_tree(&root[0]->left);
		free_tree(&root[0]->right);
		bitstring_free(&root[0]->code);
		free(*root);
		*root=0;
	}
}

void gen_histogram(const unsigned char *data, size_t nBytes, size_t *hist)
{
	memset(hist, 0, 256*sizeof(size_t));
	for(size_t i=0;i<nBytes;++i)
		++hist[data[i]];
}
HuffNodeHandle gen_tree(size_t *hist)
{
	PQueueHandle q;
	HuffNodeHandle node, L, R;

	PQUEUE_ALLOC(HuffNodeHandle, q, 0, huffnode_less, 0);
	for(int k=0;k<256;++k)
	{
		if(hist[k])
		{
			node=(HuffNodeHandle)malloc(sizeof(HuffNode));
			if(!node)
				LOG_ERROR("Huffman gen_tree: malloc failed");
			node->left=node->right=0;
			node->symbol=k;
			node->freq=hist[k];
			node->code=0;
			pqueue_enqueue(&q, &node);
		}
	}

		//pqueue_print_heap(&q, huffnode_print);//

	for(int it=0;q->count>1;++it)
	{
		//printf("%d\n", it);
		//printf("q->count = %d\n", q->count);//

		L=*(HuffNodeHandle*)q->data;
		pqueue_dequeue(&q);

		//pqueue_print_heap(&q, huffnode_print);//
		//printf("q->count = %d\n", q->count);//
		
		R=*(HuffNodeHandle*)q->data;
		pqueue_dequeue(&q);

		//pqueue_print_heap(&q, huffnode_print);//
		//printf("q->count = %d\n", q->count);//

		node=(HuffNodeHandle)malloc(sizeof(HuffNode));
		if(!node)
			LOG_ERROR("Huffman gen_tree: malloc failed");
		node->symbol=-1;//don't care for symbol because not leaf
		node->freq=L->freq+R->freq;
		node->left=L;
		node->right=R;
		node->code=0;

		//printf("%d+%d = %d\n", L->freq, R->freq, n2->freq);//
		
		pqueue_enqueue(&q, &node);

		//pqueue_print_heap(&q, huffnode_print);//
		//printf("q->count = %d\n\n", q->count);//
	}
	node=*(HuffNodeHandle*)q->data;
	pqueue_free(&q);
	return node;
}
void gen_alphabet(HuffNodeHandle root, BitstringHandle *alphabet)
{
	SList s;//stack
	HuffNodeHandle node;
	char bit;

	slist_init(&s, sizeof(HuffNodeHandle), 0);
	STACK_PUSH(&s, &root);
	while(s.count)
	{
		//slist_print(&s, huffnode_print);//
		//printf("\n");

		node=*(HuffNodeHandle*)STACK_TOP(&s);
		STACK_POP(&s);

		if(!node->left&&!node->right)//leaf node: move leaf code to alphabet
		{
			alphabet[node->symbol]=node->code;
#ifndef HUFF_LATE_FREE_TREE
			node->code=0;
#endif
			continue;
		}
		if(node->left)//append false-bit
		{
			bit=0;
			if(node->code)
				node->left->code=bitstring_construct(node->code->data, node->code->bitCount, 0, 1);
			else
				node->left->code=bitstring_construct(0, 0, 0, 1);
			bitstring_append(&node->left->code, &bit, 1, 0);
			STACK_PUSH(&s, &node->left);
		}
		if(node->right)//append true-bit
		{
			bit=1;
			if(node->code)
				node->right->code=bitstring_construct(node->code->data, node->code->bitCount, 0, 1);
			else
				node->right->code=bitstring_construct(0, 0, 0, 1);
			bitstring_append(&node->right->code, &bit, 1, 0);
			STACK_PUSH(&s, &node->right);
		}
	}
	slist_clear(&s);
}

void print_hist(size_t *hist)
{
	for(int i=0;i<256;++i)
	{
		if(hist[i])
			printf("0x%02X \'%c\':  %lld\n", i, i, (long long)hist[i]);
	}
	printf("\n");
}
void print_alphabet(BitstringHandle *alphabet, size_t *hist, size_t nsymbols, int print_chars)
{
	for(size_t sym=0;sym<nsymbols;++sym)
	{
		if(alphabet[sym])
		{
			if(print_chars)
				printf("\'%c\' %8d ", (char)sym, (int)hist[sym]);
			else
				printf("%3d %8d ", (int)sym, (int)hist[sym]);
			bitstring_print(alphabet[sym]);
			printf("\n");
		}
	}
}
#ifdef HUFF_ALPHAGUARD
size_t maxlen=0;//
#endif

HuffFileHeader	header={0};
BitstringHandle alphabet[256]={0};
int huff_compress(const void *src, size_t srcSize, ArrayHandle *dst)
{
	const unsigned char *srcbuf=(const unsigned char*)src;
	HuffNodeHandle root;
	BitstringHandle d2;
	int success;

	success=1;
#ifdef HUFF_PROGRESS
	printf("Huffman Enc: gen histogram...\n");
#endif
	gen_histogram(srcbuf, srcSize, header.hist);
#ifdef HUFF_PROGRESS
	printf("Huffman Enc: gen tree...\n");
#endif
	root=gen_tree(header.hist);
#ifdef HUFF_PROGRESS
	printf("Huffman Enc: gen alphabet...\n");
#endif
	gen_alphabet(root, alphabet);
	
	//print_tree(root, 0);//
	//print_alphabet(alphabet, header.hist, 256, 0);//
#ifdef HUFF_ALPHAGUARD
	maxlen=0;
	for(int k=0;k<256;++k)
	{
		if(alphabet[k]&&maxlen<alphabet[k]->bitCount)
			maxlen=alphabet[k]->bitCount;
	}
	printf("Huffman: Max alphabet len: %lld bits\n", maxlen);
#endif
#ifdef HUFF_PROGRESS
	printf("Huffman Enc: free tree...\n");
#endif
#ifndef HUFF_LATE_FREE_TREE
	free_tree(&root);
#endif

	d2=bitstring_construct(0, 0, 0, srcSize<<1);
	for(size_t i=0;i<srcSize;++i)
	{
#ifdef HUFF_PROGRESS
		if(srcSize<=10000||!(i%(srcSize/10000)))
			printf("\rEnc: %7d/%7d = %.2f%%...", (int)i+1, (int)srcSize, 100.*(i+1)/srcSize);
#endif
		int sym=srcbuf[i];
		BitstringHandle code=alphabet[sym];

#ifdef DEBUG_HUFF_ENC
		printf("[%d] 0x%02X: ", (int)i, sym);
		bitstring_print(code);
		printf("\n");//
#endif
		
		if(d2->bitCount&7)
		{
			for(int k=0;k<code->bitCount;k+=8)
			{
				unsigned char bits=code->data[k>>3];
				size_t byteidx=(d2->bitCount+k)>>3, bitidx=(d2->bitCount+k)&7;
				d2->data[byteidx  ]|=bits<<bitidx;
				d2->data[byteidx+1]|=bits>>(8-bitidx);
			}
		}
		else
		{
			for(int k=0;k<code->bitCount;k+=8)
				d2->data[(d2->bitCount+k)>>3]=code->data[k>>3];
		}
		d2->bitCount+=code->bitCount;
		//bitstring_append(&d2, code->data, code->bitCount, 0);

#ifdef DEBUG_HUFF_ENC
		bitstring_print(d2), printf("\n");//
#endif
	}
#ifdef HUFF_PROGRESS
	printf("\n");
#endif

#ifdef DEBUG_HUFF_ENC
	printf("\nFinal compressed bitstring:\n");
	bitstring_print(d2), printf("\n");//
#endif

	//write to dst
	header.tag=tag_huff;
	header.zeros=0;
	header.u_size=srcSize;
	header.c_size=(d2->bitCount+7)/8;
	if(!*dst)
		//ARRAY_ALLOC(unsigned char, *dst, 0, sizeof(HuffFileHeader)+header.c_size, 0, 0);
		ARRAY_ALLOC(unsigned char, *dst, 0, 0, 0, 0);
	if(dst[0]->esize!=1)
	{
		success=0;
		goto finish;
	}
	ARRAY_APPEND(*dst, &header, sizeof(header), 1, 0);
	ARRAY_APPEND(*dst, d2->data, header.c_size, 1, 0);
finish:
#ifdef HUFF_LATE_FREE_TREE
	free_tree(&root);
	memset(alphabet, 0, sizeof(alphabet));
#else
	for(int i=0;i<256;++i)
		bitstring_free(alphabet+i);
#endif
	bitstring_free(&d2);
	return success;
}
int huff_decompress(const unsigned char *src, size_t srcSize, ArrayHandle *dst)
{
	HuffNodeHandle root;
	size_t bitIdx, nSym;
	int success;

	success=1;
	if(srcSize<sizeof(header))//compressed file must contain the header with histogram
		return 0;
	
	memcpy(&header, src, sizeof(header));
	if(header.tag!=tag_huff)
		LOG_ERROR("Invalid Huffman tag 0x%08X != 0x%08X", header.tag, tag_huff);
#ifdef HUFF_PROGRESS
	printf("Huffman Dec: gen tree...\n");
#endif
	root=gen_tree(header.hist);

#if 0
	gen_alphabet(root, alphabet);//
	print_alphabet(alphabet, header.hist, 256, 0);//
#endif
	size_t dststart;
	if(*dst)
	{
		if(dst[0]->esize!=1)
		{
			success=0;
			return 0;
		}
		dststart=dst[0]->count;
		ARRAY_APPEND(*dst, 0, header.u_size, 1, 0);
	}
	else
	{
		dststart=0;
		ARRAY_ALLOC(unsigned char, *dst, 0, header.u_size, 0, 0);
	}

	src+=sizeof(header);
	srcSize-=sizeof(header);
	for(bitIdx=0, nSym=0;bitIdx<(srcSize<<3)&&nSym<header.u_size;++nSym)//naive one-bit reading
	{
#ifdef HUFF_PROGRESS
		if(header.u_size<=10000||!(nSym%(header.u_size/10000)))
			printf("\rDec: %7d/%7d = %.2f%%...", (int)nSym+1, (int)header.u_size, 100.*(nSym+1)/header.u_size);
#endif
		HuffNodeHandle node=root;
#ifdef HUFF_NODEHIST
#define HISTLEN 20
		HuffNodeHandle nhist[HISTLEN]={0};
#endif
#ifdef HUFF_ALPHAGUARD
		size_t len=0;
#endif
		//while(node->left||node->right)
		while(node->symbol==-1)//while not leaf
		{
#ifdef HUFF_ALPHAGUARD
			if(len>maxlen)
			{
				printf("Huffman: Decode error: bitIdx %lld", bitIdx);
				pause();
				exit(0);
			}
#endif
			int bit=src[bitIdx>>3]>>(bitIdx&7)&1;

#ifdef DEBUG_HUFF_DEC
			printf("%d", bit);//
#endif
#ifdef HUFF_NODEHIST
			for(int k=0;k<HISTLEN-1;++k)
				nhist[k]=nhist[k+1];
			nhist[HISTLEN-1]=node;
#endif
			node=(&node->left)[bit];
			++bitIdx;
#ifdef HUFF_ALPHAGUARD
			++len;
#endif
		}

#ifdef DEBUG_HUFF_DEC
		printf("\n0x%02X\n", (int)node->symbol);//
#endif

		dst[0]->data[dststart+nSym]=(unsigned char)node->symbol;
		//ARRAY_APPEND(*dst, &node->symbol, 1, 1, 0);
	}
#ifdef HUFF_PROGRESS
	printf("\n");
#endif
	free_tree(&root);
	return 1;
}