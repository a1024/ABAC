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

#define FIXED_POINT		20
#define FIXED_ONE		(1<<FIXED_POINT)

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

		ArrayHandle src, dst;

		src=conv_weight->data;
		ARRAY_ALLOC(int, conv2d->weight, 0, src->count, 0, 0);
		dst=conv2d->weight;
		for(int kv=0;kv<(int)src->count;++kv)
			((int*)dst->data)[kv]=(int)(((float*)src->data)[kv]*FIXED_ONE);

		src=conv_bias->data;
		ARRAY_ALLOC(int, conv2d->bias, 0, src->count, 0, 0);
		dst=conv2d->bias;
		for(int kv=0;kv<(int)src->count;++kv)
			((int*)dst->data)[kv]=(int)(((float*)src->data)[kv]*FIXED_ONE);

		//conv2d->weight=conv_weight->data;
		//conv_weight->data=0;
		//
		//conv2d->bias=conv_bias->data;
		//conv_bias->data=0;
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
			*bias=(int)roundf(val*FIXED_ONE);

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
						*dst=(int)roundf(val*FIXED_ONE);
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
						*dst=(int)roundf(val*FIXED_ONE);
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
static void conv2d(const int *src, const int *weight, const int *bias, int co, int ci, int kw, int kh, int bw, int bh, int add, int *dst)
{
	printf("conv2d [%d, %d, %d, %d] %d*%d\n", co, ci, kw, kh, bw, bh);//

	int padx=kw>>1, pady=kh>>1;
	for(int ko=0;ko<co;++ko)//for each output channel
	{
		for(int ky=pady;ky<bh-pady;++ky)//for each output row
		{
			for(int kx=padx;kx<bw-padx;++kx)//for each output pixel
			{
				int sum=bias[ko];
				for(int ki=0;ki<ci;++ki)//for each input channel
				{
					const int *kernel=weight+kw*kh*(ci*ko+ki);
					const int *channel=src+bw*bh*ki;
					for(int ky2=0;ky2<kh;++ky2)//for each kernel row
					{
						for(int kx2=0;kx2<kw;++kx2)//for each kernel coeff
							sum+=(int)((long long)kernel[kw*ky2+kx2]*channel[bw*(ky+ky2-pady)+kx+kx2-padx]>>FIXED_POINT);
					}
				}
				if(add)
					sum+=dst[bw*(bh*ko+ky)+kx];
				if(sum<0)//LeakyReLU
					sum=(long long)sum*0x28F6>>FIXED_POINT;
					//sum=(((long long)sum+50)<<FIXED_POINT)/100;
				dst[bw*(bh*ko+ky)+kx]=sum;
			}
		}
	}
}
static int* c33_enc_pred(const short *buf, int bw0, int bw, int bh, int nbits, int *aux, int aw, int ah, int maxchannels, ArrayHandle preproc, ArrayHandle mainpart)
{
	Conv2dStage *c;
	int nc=maxchannels>>1;
	int *a1=aux, *a2=aux+aw*ah*nc;

	if(preproc->count!=2||mainpart->count!=8)
		LOG_ERROR("Outdated code");

	c=(Conv2dStage*)array_at(&preproc, 0);
	conv2d(a1, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a2);

	c=(Conv2dStage*)array_at(&preproc, 1);
	conv2d(a2, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a1);

	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ah-1;++ky)
		{
			for(int kx=0;kx<aw;++kx)
				a2[aw*(ah*kc+ky+2)+kx+1]=ky>=1&&ky<ah-1&&kx>=1&&kx<aw-1?buf[3*((aw-2)*ky+kx)+kc]<<(FIXED_POINT-nbits):0;
		}
	}
	c=(Conv2dStage*)array_at(&mainpart, 0);
	conv2d(a2, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 1, a1);
	
	c=(Conv2dStage*)array_at(&mainpart, 1);
	conv2d(a1, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a2);
	
	c=(Conv2dStage*)array_at(&mainpart, 2);
	conv2d(a2, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a1);
	
	c=(Conv2dStage*)array_at(&mainpart, 3);
	conv2d(a1, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a2);

	c=(Conv2dStage*)array_at(&mainpart, 4);
	conv2d(a2, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a1);
	
	c=(Conv2dStage*)array_at(&mainpart, 5);
	conv2d(a1, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a2);
	
	c=(Conv2dStage*)array_at(&mainpart, 6);
	conv2d(a2, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a1);
	
	c=(Conv2dStage*)array_at(&mainpart, 7);
	conv2d(a1, (int*)c->weight->data, (int*)c->bias->data, c->co, c->ci, c->kw, c->kh, aw, ah, 0, a2);

	return a2;
}
static void blit_meanconf(const int *result, int rw, int rh, int x1, int x2, int y1, int y2, int *buf_mean, int *buf_conf, int bw)
{
	int dx=x2-x1, dy=y2-y1;
	int conf_limit=(FIXED_ONE<<3)-1;
	for(int ky=0;ky<dy;++ky)
	{
		for(int kx=0;kx<dx;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int idx=3*(bw*(y1+ky)+x1+kx)+kc;
				if(kx+1<rw-1||ky+1<rh-1)
				{
					int val=result[rw*(rh*kc+ky+1)+kx+1];
					//val=CLAMP(-FIXED_ONE, val, FIXED_ONE);
					if(val<-FIXED_ONE)
						val=FIXED_ONE;
					if(val>FIXED_ONE)
						val=FIXED_ONE;
					buf_mean[idx]=val;

					val=result[rw*(rh*(kc+3)+ky+1)+kx+1];
					val=abs(val);
					if(val>conf_limit)
						val=conf_limit;
					buf_conf[idx]=val;
				}
				else
				{
					buf_mean[idx]=0;
					buf_conf[idx]=0;
				}
			}
		}
	}
}
static int powersof2[FIXED_POINT+4]=
{
	0x00100001, 0x00100001, 0x00100003, 0x00100006,
	0x0010000B, 0x00100016, 0x0010002C, 0x00100059,
	0x001000B1, 0x00100163, 0x001002C6, 0x0010058D,
	0x00100B1B, 0x0010163E, 0x00102C9A, 0x001059B1,
	0x0010B558, 0x001172B8, 0x001306FE, 0x0016A09E,
	0x00200000, 0x00400000, 0x01000000, 0x10000000,
};
//void twopower_init()
//{
//	for(int k=0;k<FIXED_POINT+4;++k)
//		powersof2[k]=(int)(pow(2., pow(2., (double)(FIXED_POINT-k))+FIXED_POINT)*FIXED_ONE);
//}
int twopower(int x)
{
	int product=FIXED_ONE;
	for(int k=0;k<FIXED_POINT+4;++k)
	{
		int bit=x>>k&1;
		if(bit)
			product=(int)((long long)product*powersof2[k]>>20);
	}
	return product;
}
static int error_func_p20(int x)
{
	//approximation 3 on Wikipedia
	const unsigned c[]=
	{
		0x120DD,//0x120DCCEB,//0.0705230784
		0x0AD30,//0x0AD2FE74,//0.0422820123
		0x025F9,//0x025F8DA3,//0.0092705272
		0x0009F,//0x0009F660,//0.0001520143
		0x00122,//0x00122007,//0.0002765672
		0x0002D,//0x0002D27E,//0.0000430638
	};
	int neg=x<0;
	unsigned long long x0=(unsigned long long)abs(x), res;//16.16 bits

	//if(x0>0x32A1F)//erf(3.16453508) = 0x0.FFFF8000003 ~= 1 in 16.16 bit
	if(x0>0x38F818)//erf(0x38F818>>20) = 0x0.FFFFF800004
		res=FIXED_ONE;
	else
	{
		res=(x0*c[5]>>FIXED_POINT)+c[4];
		res=(res*x0>>FIXED_POINT)+c[3];
		res=(res*x0>>FIXED_POINT)+c[2];
		res=(res*x0>>FIXED_POINT)+c[1];
		res=(res*x0>>FIXED_POINT)+c[0];
		res=(res*x0>>FIXED_POINT)+FIXED_ONE;

		res=(1LL<<(FIXED_POINT<<1)|res>>1)/res;//(ONE<<POINT)/res
		res=res*res>>FIXED_POINT;
		res=res*res>>FIXED_POINT;
		res=res*res>>FIXED_POINT;
		res=res*res>>FIXED_POINT;

		res=FIXED_ONE-res;
	}
	res^=-neg;
	res+=neg;
	return (int)res;
}
static int c33_calc_phi(int sym, int nbits, int mean, int conf)
{
	int x;

	x=(sym<<(20-nbits))-mean;
	x=(int)((long long)x*conf>>20);
	x=error_func_p20(x);

	x+=sym;
	//x=x<<nbits|sym;

	return x;
}
static void c33_enc3(int bw, int x, int y, int c, short nbits, const int *buf_mean, const int *buf_conf, const short *buf, unsigned long long *state, DList *list)
{
	static int ncalls=0;
	++ncalls;
	int idx=3*(bw*y+x)+c, half=(1<<(nbits-1));
	int mean=buf_mean[idx],
		lgconf=buf_conf[idx],
		sym=buf[idx];
	
	//if(sym<0||sym>=(1<<nbits))
	if(sym<-half||sym>=half)
		LOG_ERROR("Range error sym 0x%04X", sym);

	int conf=twopower(lgconf);

	//adaptive ANS
	int erf_start=c33_calc_phi(0         -half, nbits, mean, conf),//always -1.0
		erf_end  =c33_calc_phi((1<<nbits)-half, nbits, mean, conf),//always 1.0 (proof required)
		erf_curr =c33_calc_phi(sym            , nbits, mean, conf),
		erf_next =c33_calc_phi(sym+1          , nbits, mean, conf);

	if(erf_start>erf_curr||erf_curr>erf_next||erf_next>erf_end)
		LOG_ERROR("Out of range: start/curr/next/end: 0x%08X, 0x%08X, 0x%08X, 0x%08X", erf_start, erf_curr, erf_next, erf_end);

	long long den=(long long)erf_end-erf_start, CDF, freq;
	CDF =((long long)(erf_curr-erf_start)<<32|den>>1)/den;
	freq=((long long)(erf_next-erf_curr )<<32|den>>1)/den;

	if(CDF+freq>0x100000000)
		LOG_ERROR("CDF is not normalized CDF 0x%08llX freq 0x%08llX", CDF, freq);
	if(!freq)
		LOG_ERROR("Zero freq");

	if(*state>=((unsigned long long)freq<<32))//renorm
	{
		dlist_push_back(list, state, 4);
		*state>>=32;
	}

	*state=*state/freq<<32|(CDF+*state%freq);//update
}
static void c33_enc2(short *buf, int *buf_mean, int *buf_conf, int bw, int x1, int x2, int y1, int y2, int nbits, unsigned long long *state, DList *list)
{
	int dx=x2-x1, dy=y2-y1;
	for(int ky=y2-1;ky>=y1;--ky)
	{
		for(int kx=x2-1;kx>=x1;--kx)
		{
			for(int kc=2;kc>=0;--kc)
			{
				c33_enc3(bw, kx, ky, kc, nbits+(kc!=1), buf_mean, buf_conf, buf, state, list);
			}
		}
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
		buf[kd  ]=src[ks  ]-src[ks|1];//XGZ {9, 8, 9} bit
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
	int alen=aw*ah*3, *result;
	
	for(int ks=0;ks<nsizes-1;++ks)
	{
		int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
		int aw2=px+2, ah2=py+2;
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
					aux[aw2*(ah2*kc+ky+1)+kx+1]=curr-left;
				}
			}
		}
		result=c33_enc_pred(buf+3*px, bw, px, py, 10, aux, aw2, ah2, 24*2, ctx.preproc_x, ctx.mainpart_x);
		blit_meanconf(result, aw2, ah2, px, dx, 0, py, buf_mean, buf_conf, bw);

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
					aux[aw2*(ah2*kc+ky+1)+kx+1]=curr-top;
				}
			}
		}
		result=c33_enc_pred(buf+3*bw*py, bw, px, py, 10, aux, aw2, ah2, 24*2, ctx.preproc_y, ctx.mainpart_y);
		blit_meanconf(result, aw2, ah2, 0, px, py, dy, buf_mean, buf_conf, bw);
			
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
					aux[aw2*(ah2*kc+ky+1)+kx+1]=curr-left-top+topleft;
				}
			}
		}
		result=c33_enc_pred(buf+3*(bw*py+px), bw, px, py, 11, aux, aw2, ah2, 24*2, ctx.preproc_xy, ctx.mainpart_xy);
		blit_meanconf(result, aw2, ah2, px, dx, py, dy, buf_mean, buf_conf, bw);
	}
	blit_meanconf(0, 0, 0, 0, psizes[nsizes-1].w, 0, psizes[nsizes-1].h, buf_mean, buf_conf, bw);//bypass top-left subband
	free(temprow);
	free(aux);

	//encode with adaptive ANS
	DList list;
	dlist_init(&list, 1, 1024, 0);
	unsigned long long state=0x100000000;
	for(int ks=0;ks<nsizes-1;++ks)
	{
		int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
		c33_enc2(buf, buf_mean, buf_conf, bw, px, dx,  0, py,  9, &state, &list);
		c33_enc2(buf, buf_mean, buf_conf, bw,  0, px, py, dy,  9, &state, &list);
		c33_enc2(buf, buf_mean, buf_conf, bw, px, dx, py, dy, 10, &state, &list);
	}
	c33_enc2(buf, buf_mean, buf_conf, bw, 0, psizes[nsizes-1].w, 0, psizes[nsizes-1].h, 9, &state, &list);
	dlist_push_back(&list, &state, 8);

#if 1
	save_32bit("buf_mean.PNG",   buf_mean, bw, bh, 3, 0);
	save_32bit("buf_lgconf.PNG", buf_conf, bw, bh, 3, 0);
#endif

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);
	
	array_free(&sizes);
	free(buf_mean);
	free(buf_conf);
	free(buf);
	return 1;
}