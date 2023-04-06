#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>
#include<ctype.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include"lodepng.h"//for debugging
static const char file[]=__FILE__;

typedef enum TensorTypeEnum
{
	TENSOR_WEIGHT,
	TENSOR_BIAS,
} TensorType;
typedef struct TensorStruct
{
	TensorType type;
	int shape[4];
	ArrayHandle name, data;
} Tensor;
static void free_tensor(void *p)
{
	Tensor *t=(Tensor*)p;
	array_free(&t->name);
	array_free(&t->data);
}

void print_param_shapes(ArrayHandle params)
{
	for(int k=0;k<(int)params->count;++k)
	{
		Tensor *t=(Tensor*)array_at(&params, k);
		printf("%s\t[%d, %d, %d, %d]\n", t->name->data, t->shape[0], t->shape[1], t->shape[2], t->shape[3]);
	}
}

static void parse_error(ArrayHandle text, ptrdiff_t idx, const char *msg)
{
	int line=0;
	ptrdiff_t linestart=0;
	for(ptrdiff_t k=0;k<idx;++k)
	{
		if(text->data[k]=='\n')
		{
			++line;
			linestart=k;
		}
	}
	LOG_ERROR("(%d : %d) %s", line+1, (int)linestart, msg);
}
static int skip_ws(ArrayHandle text, ptrdiff_t *idx)
{
	for(;*idx<(ptrdiff_t)text->count&&isspace(text->data[*idx]);++*idx);
	return *idx>=(ptrdiff_t)text->count;
}
static int skip_nonws(ArrayHandle text, ptrdiff_t *idx)
{
	for(;*idx<(ptrdiff_t)text->count&&!isspace(text->data[*idx]);++*idx);
	return *idx>=(ptrdiff_t)text->count;
}
static long long acme_atoll(ArrayHandle text, ptrdiff_t *idx, int is_signed, int *ndigits)//string must be null-terminated
{
	long long val=0;
	int neg=0;
	unsigned char c;
	if(is_signed)
	{
		c=text->data[*idx];//no need to check for *idx>=text->count
		if(c=='+')
			++*idx;
		else if(c=='-')
		{
			++*idx;
			neg=1;
		}
	}
	ptrdiff_t start=*idx;
	for(;;++*idx)
	{
		c=text->data[*idx]-'0';
		if(c>=10)
			break;
		val*=10;
		val+=c;
	}
	if(start>=*idx)
		parse_error(text, *idx, "Expected a number");
	if(neg)
		val=-val;
	if(ndigits)
		*ndigits=(int)(*idx-start);
	return val;
}
static double acme_atof(ArrayHandle text, ptrdiff_t *idx)
{
	double val=(double)acme_atoll(text, idx, 1, 0);
	int tail_len=0;
	if(text->data[*idx]=='.')
	{
		++*idx;
		long long tail=acme_atoll(text, idx, 0, &tail_len);
		val+=(double)tail*pow(10., (double)-tail_len);
	}
	if(text->data[*idx]=='e')
	{
		++*idx;
		int exponent=(int)acme_atoll(text, idx, 1, 0);
		val*=pow(10., (double)exponent);
	}
	return val;
}
ArrayHandle codec_loadweights(const char *filename)
{
	ArrayHandle text=load_file(filename, 0, 16);
	if(!text)
		return 0;
	const char *str=text->data;
	ptrdiff_t len=text->count, idx=0;
	ArrayHandle params;
	ARRAY_ALLOC(Tensor, params, 0, 0, 0, free_tensor);
	//DList list;
	//dlist_init(&list, sizeof(float), 256, 0);
	for(;;)
	{
		if(skip_ws(text, &idx))
			break;
		ptrdiff_t start=idx;
		if(skip_nonws(text, &idx))
		{
			parse_error(text, idx, "Expected a declaration");
			break;
		}
		if(start>=idx)
			parse_error(text, idx, "Expected a name");
		Tensor *t=(Tensor*)ARRAY_APPEND(params, 0, 1, 1, 0);

		//if(params->count==37)//
		//	params->count=37;//

		STR_COPY(t->name, text->data+start, idx-start);
		if(skip_ws(text, &idx))
		{
			parse_error(text, idx, "Expected dimensions");
			break;
		}

		//if(params->count==3)//
		//	params->count=3;//
		//if(!strcmp((char*)t->name->data, "predx.conv01.bias"))//
		//	t->name->data[0]='p';

		for(int kd=0;kd<4;)
		{
			t->shape[kd]=(int)acme_atoll(text, &idx, 0, 0);
			//t->shape[kd]=atoi(text->data+idx);//how many digits?
			if(skip_ws(text, &idx))
			{
				parse_error(text, idx, "Expected colon \':\'");
				break;
			}
			++kd;
			if(text->data[idx]==':')
			{
				++idx;
				if(kd<4)
				{
					for(;kd<4;++kd)//fill the rest of shape with ones
						t->shape[kd]=1;

					//for(int ksrc=kd, kdst=3;ksrc>=0;--ksrc, --kdst)//shift sizes
					//	t->shape[kdst]=t->shape[ksrc];
					//for(int k=0;k<3-kd;++k)//fill other sizes with one
					//	t->shape[k]=1;
				}
				break;
			}
		}
		
		//if(params->count==2)//
		//	params->count=2;//

		int nelem=t->shape[0]*t->shape[1]*t->shape[2]*t->shape[3];
		ARRAY_ALLOC(float, t->data, 0, nelem, 0, 0);
		int kv;
		for(kv=0;kv<nelem;++kv)
		{
			if(skip_ws(text, &idx))
			{
				parse_error(text, idx, "Missing values");
				break;
			}
			float *val=(float*)array_at(&t->data, kv);
			*val=(float)acme_atof(text, &idx);
		}
		if(kv!=nelem)
			parse_error(text, idx, "Wrong number of arguments");
	}
	array_free(&text);
	return params;
}

#if 0
static int match_name(const char *name, const char *obj, const char *member, const char *type)
{
	int p1, p2, match=1, k;
	const char *temp;

	for(p1=0;name[p1]&&name[p1]!='.';++p1);
	for(p2=p1+1;name[p2]&&name[p2]!='.';++p2);
	if(obj)
	{
		for(k=0;obj[k]&&name[k]!='.'&&obj[k]==name[k];++k);
		match&=obj[k]||name[k]!='.';
	}
	if(member)
	{
		temp=name+p1+1;
		for(k=0;member[k]&&temp[k]!='.'&&member[k]==temp[k];++k);
		match&=member[k]||temp[k]!='.';
	}
	if(type)
	{
		temp=name+p2+1;
		for(k=0;type[k]&&temp[k]!='.'&&type[k]==temp[k];++k);
		match&=type[k]||temp[k]!='.';
	}
	return match;
}
static int find_tensor(Tensor *ptr, int start, int ntensors, const char *obj, const char *member, const char *type)
{
	int ret;
	for(ret=start;ret<ntensors;++ret)
	{
		Tensor *t=ptr+ret;
		if(match_name(t->name->data, obj, member, type))
			break;
	}
	return ret;
}
#endif
static int find_first_tensor(Tensor *ptr, int ntensors, int start, const char *objname, int namelen)//linear search
{
	int ret;
	for(ret=start;ret<ntensors;++ret)
	{
		Tensor *t=ptr+ret;
		if(!strncmp(t->name->data, objname, namelen))
			break;
	}
	return ret;
}
static int find_tensor2(Tensor *ptr, int ntensors, const char *format, ...)
{
	char buf[128];
	int ret;
	va_list args;

	va_start(args, format);
	vsnprintf(buf, 128, format, args);
	va_end(args);

	for(ret=0;ret<ntensors;++ret)
	{
		Tensor *t=ptr+ret;
		if(t->name&&!strcmp(t->name->data, buf))
			break;
	}
	if(ret==ntensors)
		LOG_ERROR("Couln't find %s", buf);
	return ret;
}
typedef struct Conv2dStageStruct
{
	int co, ci, kh, kw;
	ArrayHandle weight, bias;
} Conv2dStage;
typedef struct C33PredStruct
{
	ArrayHandle//all arrays are of Conv2dStage
		preproc_x, mainpart_x,
		preproc_y, mainpart_y,
		preproc_xy, mainpart_xy;
} C33Pred;
void free_conv2d(void *p)
{
	Conv2dStage *c=(Conv2dStage*)p;
	array_free(&c->weight);
	array_free(&c->bias);
}
static void c33_init_pred(Tensor *ptr, int ntensors, const char *name, ArrayHandle *preproc, ArrayHandle *mainpart)
{
	int npreprocstages=2, nmainstages=8;
	ARRAY_ALLOC(Conv2dStage, *preproc, 0, npreprocstages, 0, free_conv2d);
	for(int k=0;k<npreprocstages;++k)
	{
		int idx_conv_weight=find_tensor2(ptr, ntensors, "%s.conv%02d.weight", name, k+1),
			idx_conv_bias  =find_tensor2(ptr, ntensors, "%s.conv%02d.bias", name, k+1);
		//if(idx_conv_weight==ntensors)
		//	LOG_ERROR("%s.conv%02d.weight", name, k);
		//if(idx_conv_bias==ntensors)
		//	LOG_ERROR("%s.conv%02d.bias", name, k);
		Tensor
			*conv_weight=ptr+idx_conv_weight,
			*conv_bias=ptr+idx_conv_bias;
		Conv2dStage *conv2d=(Conv2dStage*)array_at(preproc, k);
		memcpy(&conv2d->co, conv_weight->shape, sizeof(int[4]));

		conv2d->weight=conv_weight->data;
		conv_weight->data=0;
		
		conv2d->bias=conv_bias->data;
		conv_bias->data=0;
	}
	ARRAY_ALLOC(Conv2dStage, *mainpart, 0, nmainstages, 0, free_conv2d);
	for(int k=0;k<nmainstages;++k)
	{
		//int idx_ct_weight=find_tensor2(ptr, ntensors, "%s.ct%02d.weight", name, k),
		//	idx_cl_weight=find_tensor2(ptr, ntensors, "%s.cl%02d.weight", name, k),
		//	idx_bn_weight=find_tensor2(ptr, ntensors, "%s.b%02d.weight", name, k),
		//	idx_ct_bias=find_tensor2(ptr, ntensors, "%s.ct%02d.bias", name, k),
		//	idx_cl_bias=find_tensor2(ptr, ntensors, "%s.cl%02d.bias", name, k),
		//	idx_bn_bias=find_tensor2(ptr, ntensors, "%s.b%02d.bias", name, k);
		Tensor
			*ct_weight=ptr+find_tensor2(ptr, ntensors, "%s.ct%02d.weight", name, k+1),
			*cl_weight=ptr+find_tensor2(ptr, ntensors, "%s.cl%02d.weight", name, k+1),
			*bn_weight=ptr+find_tensor2(ptr, ntensors, "%s.b%02d.weight", name, k+1),
			*ct_bias  =ptr+find_tensor2(ptr, ntensors, "%s.ct%02d.bias", name, k+1),
			*cl_bias  =ptr+find_tensor2(ptr, ntensors, "%s.cl%02d.bias", name, k+1),
			*bn_bias  =ptr+find_tensor2(ptr, ntensors, "%s.b%02d.bias", name, k+1);
		Conv2dStage *conv2d=(Conv2dStage*)array_at(mainpart, k);

		if(ct_weight->shape[0]!=ct_bias->shape[0]||
			cl_weight->shape[0]!=cl_bias->shape[0]||
			ct_weight->shape[0]!=cl_weight->shape[0]||
			ct_weight->shape[1]!=cl_weight->shape[1]||
			ct_bias  ->shape[1]!=1||ct_bias  ->shape[2]!=1||ct_bias  ->shape[3]!=1||
			cl_bias  ->shape[1]!=1||cl_bias  ->shape[2]!=1||cl_bias  ->shape[3]!=1||
			bn_weight->shape[1]!=1||bn_weight->shape[2]!=1||bn_weight->shape[3]!=1||
			bn_bias  ->shape[1]!=1||bn_bias  ->shape[2]!=1||bn_bias  ->shape[3]!=1||
			bn_weight->shape[0]!=ct_weight->shape[0]||bn_bias->shape[0]!=ct_weight->shape[0]||
			cl_weight->shape[2]!=1)
			LOG_ERROR("Wrong filter shapes CT.W[%d, %d, %d, %d], CT.B[%d, %d, %d, %d], CL.W[%d, %d, %d, %d], CL.B[%d, %d, %d, %d]",
				ct_weight->shape[0], ct_weight->shape[1], ct_weight->shape[2], ct_weight->shape[3],
				ct_bias  ->shape[0], ct_bias  ->shape[1], ct_bias  ->shape[2], ct_bias  ->shape[3],
				cl_weight->shape[0], cl_weight->shape[1], cl_weight->shape[2], cl_weight->shape[3],
				cl_bias  ->shape[0], cl_bias  ->shape[1], cl_bias  ->shape[2], cl_bias  ->shape[3]);

		conv2d->co=ct_weight->shape[0];
		conv2d->ci=ct_weight->shape[1];
		conv2d->kh=ct_weight->shape[2]+cl_weight->shape[2];
		conv2d->kw=ct_weight->shape[3];
		
		//(ct0X.weight | cl0X.weight)*b0X.weight  /  (ct0X.bias + cl0X.bias)*b0X.weight + b0X.bias
		int nweights=conv2d->co*conv2d->ci*conv2d->kh*conv2d->kw;
		ARRAY_ALLOC(int, conv2d->weight, 0, nweights, 0, 0);
		ARRAY_ALLOC(int, conv2d->bias, 0, conv2d->co, 0, 0);
		for(int ko=0, srcidx=0, srcidx2=0, dstidx=0;ko<conv2d->co;++ko)//for each filter
		{
			float
				ct_bk=*(float*)array_at(&ct_bias->data, ko),
				cl_bk=*(float*)array_at(&cl_bias->data, ko),
				b_wk=*(float*)array_at(&bn_weight->data, ko),
				b_bk=*(float*)array_at(&bn_bias->data, ko);

			int *bias=(int*)array_at(&conv2d->bias, ko);
			float val=(ct_bk+cl_bk)*b_wk+b_bk;
			*bias=(int)roundf(val*0x10000);

			for(int ki=0;ki<conv2d->ci;++ki)//for each kernel
			{
				for(int ky=0;ky<ct_weight->shape[2];++ky)
				{
					for(int kx=0;kx<conv2d->kw;++kx, ++dstidx, ++srcidx)
					{
						float *src=(float*)array_at(&ct_weight->data, srcidx);
						int *dst=(int*)array_at(&conv2d->weight, dstidx);
						//float *src=(float*)array_at(&ct_weight->data, conv2d->kw*(ct_weight->shape[2]*((size_t)conv2d->ci*ko+ki)+ky)+kx);
						//float *dst=(float*)array_at(&conv2d->weight, conv2d->kw*(conv2d->kh*((size_t)conv2d->ci*ko+ki)+ky)+kx);
						float val=*src*b_wk;
						*dst=(int)roundf(val*0x10000);
					}
				}
				for(int ky=0;ky<cl_weight->shape[2];++ky)//left filter height is always 1
				{
					int kx=0;
					for(;kx<cl_weight->shape[3];++kx, ++dstidx, ++srcidx2)
					{
						float *src=(float*)array_at(&cl_weight->data, srcidx2);
						int *dst=(int*)array_at(&conv2d->weight, dstidx);
						float val=*src*b_wk;
						*dst=(int)roundf(val*0x10000);
					}
					for(;kx<conv2d->kw;++kx, ++dstidx)
					{
						int *dst=(int*)array_at(&conv2d->weight, dstidx);
						*dst=0;
					}
				}
			}
		}
	}
}
static void find_minmax(ArrayHandle params, float *mm_abs)
{
	Tensor *ptr=(Tensor*)params->data;
	int ntensors=(int)params->count;
	mm_abs[0]=0;
	mm_abs[1]=0;
	for(int k=0;k<ntensors;++k)
	{
		Tensor *t=ptr+k;
		int nelem=t->shape[0]*t->shape[1]*t->shape[2]*t->shape[3];
		for(int k2=0;k2<nelem;++k2)
		{
			float val=*(float*)array_at(&t->data, k2);
			val=fabsf(val);
			if(mm_abs[0]>val)
				mm_abs[0]=val;
			if(mm_abs[1]<val)
				mm_abs[1]=val;
		}
	}
}
static C33Pred ctx={0};
void codec33_init(const char *filename)
{
	if(!ctx.mainpart_x)
	{
		ArrayHandle params=codec_loadweights(filename);

		//float minmax[2];
		//find_minmax(params, minmax);
		//printf("min %f max %f\n", minmax[0], minmax[1]);

		Tensor *ptr=(Tensor*)params->data;
		int ntensors=(int)params->count;
		int bookmarks[4];

		bookmarks[0]=0;
		bookmarks[1]=find_first_tensor(ptr, ntensors, bookmarks[0], "predy.", 6);
		bookmarks[2]=find_first_tensor(ptr, ntensors, bookmarks[1], "predxy.", 7);
		bookmarks[3]=ntensors;

		c33_init_pred(ptr+bookmarks[0], bookmarks[1]-bookmarks[0], "predx" , &ctx.preproc_x , &ctx.mainpart_x);
		c33_init_pred(ptr+bookmarks[1], bookmarks[2]-bookmarks[1], "predy" , &ctx.preproc_y , &ctx.mainpart_y);
		c33_init_pred(ptr+bookmarks[2], bookmarks[3]-bookmarks[2], "predxy", &ctx.preproc_xy, &ctx.mainpart_xy);
		array_free(&params);
	}
}
static void c33_enc_pred(const short *buf, int bw0, int bw, int bh, int *aux, int aw, int ah, ArrayHandle preproc, ArrayHandle mainpart, int *mean, int *conf)
{
	Conv2dStage *pconv=(Conv2dStage*)preproc->data;
	int npconv=(int*)preproc->count;
	for(int k=0;k<(int)preproc->count;++k)
	{
		Conv2dStage *pconv=array_at(&preproc, k);
	}
}
size_t codec33_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	if(!ctx.mainpart_x)
	{
		LOG_ERROR("Codec33 is not initialized");
		return 0;
	}
	size_t res=(size_t)bw*bh, srclen=res<<2;
	int maxdim=MAXVAR(bw, bh);
	short
		*buf=(short*)malloc(res*3*sizeof(short)),
		*temprow=(short*)malloc(maxdim*sizeof(short));
	int *buf_mean=(int*)malloc(res*3*sizeof(int)),
		*buf_conf=(int*)malloc(res*3*sizeof(int));
	if(!buf||!temprow||!buf_mean||!buf_conf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	for(int ks=0, kd=0;ks<srclen-3;ks+=4, kd+=3)//color transform
	{
		buf[kd  ]=src[ks  ]-src[ks|1];//XGZ
		buf[kd+1]=src[ks|1]-128;
		buf[kd+2]=src[ks|2]-src[ks|1];

		//buf[kd  ]=src[ks  ]-128;//identity
		//buf[kd+1]=src[ks|1]-128;
		//buf[kd+2]=src[ks|2]-128;
	}

	int dstop=3;
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, dstop, dstop, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;

	int aw=(bw>>1)+2, ah=(bh>>1)+2;
	size_t auxlen=(size_t)aw*ah*24*2*sizeof(int);
	int *aux=(int*)malloc(auxlen);
	if(!aux)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(aux, 0, auxlen);
	int alen=aw*ah*3;
	
	for(int ks=0;ks<nsizes-1;++ks)
	{
		int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
		for(int kc=0;kc<3;++kc)
			squeeze_2d_fwd(buf+kc, psizes, ks, ks+2, 3, 0, temprow);//squeeze transform from JPEG XL

		//diffx
		memset(aux, 0, alen);
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				for(int kc=0;kc<3;++kc)
				{
					int idx=3*(bw*ky+kx)+kc;
					short
						curr=buf[idx],
						left=kx?buf[idx-3]:0;
					aux[aw*(ky+1)+kx+1]=curr-left;
				}
			}
		}
		c33_enc_pred(buf+3*px, bw, px, py, aux, aw, ah, ctx.preproc_x, ctx.mainpart_x, buf_mean, buf_conf);

		//diffy
		memset(aux, 0, alen);
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				for(int kc=0;kc<3;++kc)
				{
					int idx=3*(bw*ky+kx)+kc;
					short
						curr=buf[idx],
						top=ky?buf[idx-3*bw]:0;
					aux[aw*(ky+1)+kx+1]=curr-top;
				}
			}
		}
		c33_enc_pred(buf+3*bw*py, bw, px, py, aux, aw, ah, ctx.preproc_y, ctx.mainpart_y, buf_mean, buf_conf);
			
		//diffxy
		memset(aux, 0, alen);
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				for(int kc=0;kc<3;++kc)
				{
					int idx=3*(bw*ky+kx)+kc;
					short
						curr=buf[idx],
						left=kx?buf[idx-3]:0,
						top=ky?buf[idx-3*bw]:0,
						topleft=kx&&ky?buf[idx-3*bw-3]:0;
					aux[aw*(ky+1)+kx+1]=curr-left-top+topleft;
				}
			}
		}
		c33_enc_pred(buf+3*(bw*py+px), bw, px, py, aux, aw, ah, ctx.preproc_xy, ctx.mainpart_xy, buf_mean, buf_conf);
	}
	array_free(&sizes);
	free(temprow);
	free(aux);


	DList list;
	dlist_init(&list, 1, 1024, 0);

	//encode


	dlist_appendtoarray(&list, data);
	dlist_clear(&list);

	free(buf_mean);
	free(buf_conf);
	free(buf);
	return 1;
}