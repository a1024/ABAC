#include"fast.h"
#include"huffman.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
static const char file[]=__FILE__;

#if 0
static int huffnode_less(const void *left, const void *right)
{
	HuffInternalNodeHandle L=*(HuffInternalNodeHandle*)left, R=*(HuffInternalNodeHandle*)right;
	if(L->freq==R->freq)
		return L->symbol<R->symbol;
	return L->freq>R->freq;//for min heap
}
HuffInternalNodeHandle huff_tree_gen(unsigned *hist, int nlevels, int exclude_zps)
{
	PQueueHandle q;
	HuffInternalNodeHandle node, L, R;
	PQUEUE_ALLOC(HuffInternalNodeHandle, q, 0, huffnode_less, 0);
	for(int ks=0;ks<nlevels;++ks)
	{
		unsigned freq=hist[ks];
		if(!exclude_zps||freq)
		{
			node=(HuffInternalNodeHandle)malloc(sizeof(HuffInternalNode));
			if(!node)
			{
				LOG_ERROR2("Alloc error");
				return 0;
			}
			node->left=0;
			node->right=0;
			node->symbol=ks;
			node->freq=freq;
			node->code=0;
			pqueue_enqueue(&q, &node);
		}
	}
	while(q->count>1)
	{
		L=*(HuffInternalNodeHandle*)q->data;
		pqueue_dequeue(&q);
		R=*(HuffInternalNodeHandle*)q->data;
		pqueue_dequeue(&q);
		node=(HuffInternalNodeHandle)malloc(sizeof(HuffInternalNode));
		if(!node)
		{
			LOG_ERROR2("Alloc error");
			return 0;
		}
		node->left=L;
		node->right=R;
		node->symbol=MINVAR(L->symbol, R->symbol);
		node->freq=L->freq+R->freq;
		node->code=0;
		pqueue_enqueue(&q, &node);
	}
	node=*(HuffInternalNodeHandle*)q->data;
	pqueue_free(&q);
	return node;
}
void huff_tree_free(HuffInternalNodeHandle *root)
{
	if(*root)
	{
		huff_tree_free(&root[0]->left);
		huff_tree_free(&root[0]->right);
		bitstring_free(&root[0]->code);
		free(*root);
		*root=0;
	}
}
void huff_alphabet_gen(HuffInternalNodeHandle root, BitstringHandle *alphabet)
{
	SList s;//stack
	HuffInternalNodeHandle node;
	char bit;
	
	slist_init(&s, sizeof(HuffInternalNodeHandle), 0);
	STACK_PUSH(&s, &root);
	while(s.count)
	{
		//slist_print(&s, huffnode_print);//
		//printf("\n");

		node=*(HuffInternalNodeHandle*)STACK_TOP(&s);
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
void huff_alphabet_free(BitstringHandle *alphabet, int nlevels)
{
	for(int ks=0;ks<nlevels;++ks)
		bitstring_free(alphabet+ks);
}
void huff_alphabet_print(unsigned *hist, BitstringHandle const *alphabet, int nlevels)
{
	for(int ks=0;ks<256;++ks)
	{
		if(hist[ks])
		{
			printf("%4d  %7d  ", ks, hist[ks]);
			bitstring_print(alphabet[ks]);
			printf("\n");
		}
	}
	printf("\n");
}

HuffmanNode* huff_dectree_gen(HuffInternalNodeHandle root)
{
	return 0;//TODO
}
void huff_dectree_free(HuffmanNode **root)
{
	if(*root)
	{
		for(int k=0;k<(1<<HUFF_DEC_BITS);++k)
			huff_dectree_free(root[0]->ch+k);
		free(*root);
		*root=0;
	}
}
void huff_dec_init(HuffmanCoder *ec, const unsigned char *start, const unsigned char *end)
{
	ec->cache=0;
	ec->nbits=0;
	ec->is_enc=0;
	ec->srcptr=start;
	ec->srcend=end;
	ec->srcstart=start;
}

void huff_test(Image const *image)
{
	int nlevels=1<<image->depth, half=nlevels>>1;
	int clevels=floor_log2_32(nlevels)+1;
	int *hist=(int*)malloc(sizeof(int)*image->nch*clevels);
	BitstringHandle *alphabet=(BitstringHandle*)malloc(sizeof(BitstringHandle)*image->nch*clevels);
	if(!hist||!alphabet)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, sizeof(int)*image->nch*clevels);
	ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
	int W[4]={0};
	for(ptrdiff_t k=0;k<res;++k)
	{
		short comp[4]={0};
		memcpy(comp, image->data+image->nch*k, sizeof(short)*image->nch);
		for(int kc=1;kc<image->nch;++kc)
			comp[kc]-=comp[0];
		for(int kc=0;kc<image->nch;++kc)
		{
			//int *curr_hist=hist+((size_t)kc<<image->depth);
			int val=comp[kc]-W[kc];
			val+=half;
			val&=nlevels-1;
			val-=half;
			val=val<<1^-(val<0);
			val=floor_log2_32(val)+1;
			//if(val>=9)
			//	LOG_ERROR("");
			++hist[clevels*kc+val];
			//printf("%7d 0x%04X", (int)(image->nch*k+kc), val);
			//if((kc<<image->depth|val)>=(image->nch<<image->depth))
			//	LOG_ERROR("");
			//++hist[kc<<image->depth|val];
			//printf(" %7d\n", hist[kc<<image->depth|val]);
		}
		memcpy(W, comp, sizeof(short)*image->nch);
	}
	memset(alphabet, 0, sizeof(BitstringHandle)*image->nch*clevels);
	for(int kc=0;kc<image->nch;++kc)
	{
		int *curr_hist=hist+(size_t)clevels*kc;
		BitstringHandle *curr_alphabet=alphabet+(size_t)clevels*kc;

		double ttree=time_sec();
		HuffInternalNodeHandle tree=huff_tree_gen((unsigned*)curr_hist, clevels, 1);
		ttree=time_sec()-ttree;

		double talphabet=time_sec();
		huff_alphabet_gen(tree, curr_alphabet);
		talphabet=time_sec()-talphabet;

		double tfree=time_sec();
		huff_tree_free(&tree);
		tfree=time_sec()-tfree;


		printf("C%d  TAF %lf %lf %lf sec\n", kc, ttree, talphabet, tfree);
		huff_alphabet_print((unsigned*)curr_hist, curr_alphabet, clevels);
		printf("\n");
		huff_alphabet_free(curr_alphabet, clevels);
	}
	free(hist);
	free(alphabet);
}
#endif