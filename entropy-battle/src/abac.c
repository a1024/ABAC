#include"util.h"
#include<stdio.h>
#include<string.h>
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
static const char file[]=__FILE__;


	#define ABAC2_CONF_MSB_RELATION
	#define AC_MEASURE_PREDICTION


#define LOG_WINDOW_SIZE  16	//[2, 16]	do not change
#define LOG_CONFBOOST    14
#define WINDOW_SIZE      (1<<LOG_WINDOW_SIZE)
#define PROB_MAX         (WINDOW_SIZE-2)
#define PROB_HALF        (1<<(LOG_WINDOW_SIZE-1))
#define PROB_INIT        (PROB_HALF-1)
#define BOOST_POWER      4.
#define MIN_CONF         0.55
#define FAIL(RETTYPE, PROXYTYPE, MSG, ...)   return (RETTYPE)(PROXYTYPE)(LOG_ERROR(MSG, ##__VA_ARGS__), 0)
#define PROB_CLAMP(PROB) PROB=PROB<1?1:(PROB_MAX<PROB?PROB_MAX:PROB)
static const int magic_ac04='A'|'C'<<8|'0'<<16|'4'<<24;

#ifdef DEBUG_ABAC2
int examined_plane=9, examined_start=0, examined_end=100;
#endif
int abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud)
{
	DList list;
	const unsigned char *buffer;
	unsigned *header, *sizes, *conf;
	long long t1, t2;

	if(!src||!symcount||!bitdepth||!bytestride)
		FAIL(int, int, "abac4_encode(src=%p, symcount=%d, bitdepth=%d, stride=%d)", src, symcount, bitdepth, bytestride);
	t1=__rdtsc();

	dlist_init(&list, 1, 1024, 0);//bitdepth < 127 bit			for max bitdepth of 32: objpernode >= 260
	header=dlist_push_back(&list, 0, (1+bitdepth*2)*sizeof(int));//magic + sizes[bitdepth] + conf[bitdepth]		all must fit in list node
	header[0]=magic_ac04;
	sizes=header+1;//don't worry about alignment here
	conf=sizes+bitdepth;

	buffer=(const unsigned char*)src;

#ifdef AC_MEASURE_PREDICTION
	unsigned long long hitnum=0, hitden=0;//prediction efficiency
#endif

	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		size_t out_planestart=list.nobj;
		int bit_offset=(bitoffset+kp)>>3, bit_shift=(bitoffset+kp)&7;
		int bit_offset2=(kp+1)>>3, bit_shift2=(kp+1)&7;
		unsigned prob=PROB_HALF, prob_correct=PROB_HALF;//cheap weighted average predictor

		unsigned long long hitcount=1;
		
#ifdef PRINT_ABAC
		printf("\n");//
#endif
		for(int kb=0, kb2=0;kb<symcount;++kb, kb2+=bytestride)//analyze bitplane
		{
			int bit=buffer[kb2+bit_offset]>>bit_shift&1;
			int p0=prob-PROB_HALF;
			p0=((long long)p0*prob_correct>>16);
			p0=((long long)p0*p0>>16);
			p0+=PROB_HALF;
			//unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
			PROB_CLAMP(p0);
			int correct=bit^(p0>=PROB_HALF);
#ifdef PRINT_ABAC
			printf("%d", bit);//
#endif
			//if(kp==1)
			//	printf("%d", bit);//actual bits
			//	printf("%d", p0<PROB_HALF);//predicted bits
			//	printf("%d", !correct);//prediction error
			hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
		}
#ifdef PRINT_ABAC
		printf("\n");//
#endif
		conf[kp]=(int)hitcount;

		if(hitcount<symcount*MIN_CONF)//incompressible, bypass
		{
#ifdef PRINT_ABAC
			printf("BYPASS\n");//
#endif
			for(int kb=0, kb2=bit_offset;kb+7<symcount;kb+=8)
			{
				unsigned char frag=0;
				frag|= buffer[kb2]>>bit_shift&1, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<1, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<2, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<3, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<4, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<5, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<6, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<7, kb2+=bytestride;
				dlist_push_back1(&list, &frag);
			}
		}
		else
		{
			int hitratio_sure=(int)(0x10000*pow((double)hitcount/symcount, 1/BOOST_POWER)), hitratio_notsure=(int)(0x10000*pow((double)hitcount/symcount, BOOST_POWER));
			int hitratio_delta=hitratio_sure-hitratio_notsure;
		//	int hitratio_sure=int(0x10000*cbrt((double)hitcount/symcount)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)symcount*symcount*symcount));
		//	int hitratio_sure=int(0x10000*sqrt((double)hitcount/symcount)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)symcount*symcount));
			hitcount=(hitcount<<16)/symcount;

			//hitcount=unsigned(((unsigned long long)hitcount<<16)/symcount);
			//hitcount=abac2_normalize16(hitcount, logimsize);
			//hitcount*=invimsize;

			prob_correct=prob=PROB_HALF;

#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit0=0;
#endif
			
			unsigned start=0;
			unsigned long long range=0xFFFFFFFF;
			for(int kb=0, kb2=0;kb<symcount;kb2+=bytestride)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
			{
				int bit=buffer[kb2+bit_offset]>>bit_shift&1;
#ifdef ABAC2_CONF_MSB_RELATION
				int prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
			//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
				
				if(range<3)
				{
					dlist_push_back1(&list, (char*)&start+3);//big endian
					dlist_push_back1(&list, (char*)&start+2);
					dlist_push_back1(&list, (char*)&start+1);
					dlist_push_back1(&list, &start);

					start=0, range=0xFFFFFFFF;//because 1=0.9999...
				}
				
				//if(kp==7&&kb==5)
				//	int LOL_1=0;
				int p0=prob-PROB_HALF;
				p0=(int)((long long)p0*prob_correct>>16);
				p0=(int)((long long)p0*p0>>16);
				int sure=-(prevbit==prevbit0)&-(kp!=bitdepth-1);
				p0=(int)((long long)p0*(hitratio_notsure+(hitratio_delta&sure))>>16);
				//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
				//p0=(long long)p0*hitcount>>16;
				p0+=PROB_HALF;
				//if(prevbit!=prevbit0)
				//	p0=PROB_HALF;
				//	p0=0xFFFF-p0;

				//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

				//unsigned p0=(long long)(prob-PROB_HALF)*sqrthitcount>>16;
				//if(prevbit==prevbit0)
				//	p0=(long long)p0*hitcount>>16;
				//p0+=PROB_HALF;

				//int confboost=prevbit==prevbit0;
				//confboost-=!confboost;
				//confboost<<=LOG_CONFBOOST;
				//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(hitcount+confboost)>>16);

			//	unsigned p0=PROB_HALF+(int)((prob-PROB_HALF)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/symcount):(double)test_conf[kp]*test_conf[kp]/((double)symcount*symcount)));
			//	unsigned p0=prevbit==prevbit0?prob:PROB_HALF;
			//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*test_conf[kp]/symcount;
			//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
				PROB_CLAMP(p0);
				unsigned r2=(unsigned)(range*p0>>16);
				r2+=(r2==0)-(r2==range);
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("%6d %6d %d %08X+%08X %04X %08X %08X\n", kp, kb, bit, start, (int)range, p0, r2, start+r2);
#endif

				int correct=bit^(p0>=PROB_HALF);
			//	hitcount+=correct;
				prob=!bit<<15|prob>>1;
				prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
				prevbit0=prevbit;
#endif
#ifdef AC_MEASURE_PREDICTION
				hitnum+=correct, ++hitden;
#endif
				unsigned start0=start;
#ifdef PRINT_ABAC
				unsigned long long r0=range;
#endif
				if(bit)
				{
					++r2;
					start+=r2, range-=r2;
				}
				//	start=middle+1;
				else
					range=r2-1;
				//	end=middle-1;
#ifdef PRINT_ABAC
				if(kp==0)
					printf("%d %d bit %d p0 %04X %08X+%08X -> %08X+%08X\n", kp, kb, bit, p0, start0, r0, start, range);
#endif
				if(start<start0)//
				{
					FAIL(int, int, "AC OVERFLOW: start = %08X -> %08X, r2 = %08X", start0, start, r2);
					//printf("OVERFLOW\nstart = %08X -> %08X, r2 = %08X", start0, start, r2);
					//int k=0;
					//scanf_s("%d", &k);
				}
				++kb;
				
				while((start^(start+(unsigned)range))<0x1000000)//most significant byte has stabilized			zpaq 1.10
				{
#ifdef DEBUG_ABAC2
					if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
						printf("range 0x%08X byte-out 0x%02X\n", (int)range, start>>24);
#endif
					dlist_push_back1(&list, (char*)&start+3);
					start<<=8;
					range=range<<8|0xFF;
				}
			}
			//start+=range>>1;
			dlist_push_back1(&list, (char*)&start+3);//big endian
			dlist_push_back1(&list, (char*)&start+2);
			dlist_push_back1(&list, (char*)&start+1);
			dlist_push_back1(&list, &start);
		}
		if(loud)
			printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, conf[kp]-1, symcount, 100.*(conf[kp]-1)/symcount);
		sizes[kp]=(int)(list.nobj-out_planestart);
	}
	size_t offset0=*output?output[0]->count*output[0]->esize:0;
	dlist_appendtoarray(&list, output);
	t2=__rdtsc();

	if(loud)
	{
		int original_bitsize=symcount*bitdepth, compressed_bitsize=(int)(list.nobj)<<3;
		printf("AC encode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(original_bitsize>>3));
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize, (double)compressed_bitsize/symcount);
#ifdef AC_MEASURE_PREDICTION
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", (symcount+7)>>3);
		for(int k=0;k<bitdepth;++k)
			printf("%2d\t%5d\t%lf\n", bitdepth-1-k, sizes[k], (double)symcount/(sizes[k]<<3));
		
		printf("Preview:\n");
		int kprint=list.nobj<200?(int)list.nobj:200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", output[0]->data[offset0+k]&0xFF);
		printf("\n");
	}
	size_t csize=list.nobj;
	dlist_clear(&list);
	return (int)csize;
}
const void* abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud)//set the dst buffer to zero
{
	const unsigned char *data;
	unsigned char *buffer;
	int sizes[64], conf[64];
	long long t1, t2;

	if(!in_start||!in_end||!dst||!imsize||!bitdepth||!bytestride)
		FAIL(void*, size_t, "abac4_decode(in_start=%p, in_end=%p, dst=%p, imsize=%d, bitdepth=%d, stride=%d)", in_start, in_end, dst, imsize, bitdepth, bytestride);
	t1=__rdtsc();
	data=(const unsigned char*)in_start;
	buffer=(unsigned char*)dst;
	//memset(buffer, 0, imsize*bytestride);
	
#define PTR2OFFSET(PTR)		(int)((ptrdiff_t)(PTR)-(ptrdiff_t)in_start)
	int headercount=1+(bitdepth<<1);
	if((int)(headercount*sizeof(int))>=(unsigned char*)in_end-data)
		FAIL(void*, size_t, "File is %d bytes < %d", (int)((unsigned char*)in_end-data), headercount*sizeof(int));
	int magic;
	memcpy(&magic, data, sizeof(int));
	if(magic!=magic_ac04)
		FAIL(void*, size_t, "Invalid magic number 0x%08X, expected 0x%08X", magic, magic_ac04);
	data+=sizeof(int);
	memcpy(sizes, data, bitdepth*sizeof(int));
	data+=bitdepth*sizeof(int);
	memcpy(conf, data, bitdepth*sizeof(int));
	data+=bitdepth*sizeof(int);

#if 0
	printf("Compressed plane sizes:\n");
	for(int k=0;k<bitdepth;++k)
		printf("kp %d\t%d\n", k, sizes[k]);
#endif

	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop
	{
		int kp2=kp+bitoffset;
		int bit_offset=kp2>>3, bit_shift=kp2&7;
		int bit_offset2=(kp2+1)>>3, bit_shift2=(kp2+1)&7;
		int ncodes=sizes[kp];
		const unsigned char *plane_end=data+ncodes;
		//memcpy(&ncodes, sizes+kp*sizeof(int), sizeof(int));
		
		unsigned prob=PROB_HALF, prob_correct=PROB_HALF;
#if 1
		unsigned long long hitcount=conf[kp];
		//memcpy(&hitcount, conf+kp*sizeof(int), sizeof(int));
		if(hitcount<imsize*MIN_CONF)
		{
			if(plane_end>(const unsigned char*)in_end)
				FAIL(void*, size_t, "Decode error (bypass): plane end %d > input end %d", PTR2OFFSET(plane_end), PTR2OFFSET(in_end));
			for(int kb=0, kb2=0;kb<imsize;++kb, kb2+=bytestride)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				const unsigned char *b=data+byte_idx;
				//if((const unsigned char*)in_end<b)
				if(b>=plane_end)
					FAIL(void*, size_t, "Decode error (bypass): ptr %d > plane %d ~ %d", PTR2OFFSET(b), PTR2OFFSET(plane_end-ncodes), PTR2OFFSET(plane_end));
				int bit=*b>>bit_idx&1;
				buffer[kb2+bit_offset]|=bit<<bit_shift;
			//	buffer[kb]|=bit<<kp;
#ifdef ENABLE_GUIDE
				int original_bit=((unsigned char*)guide)[kb2+bit_offset]>>bit_shift&1;
				if(bit!=original_bit)
					FAIL(void*, size_t, "Decode error (bypass): data %d expected %d got %d", kp, original_bit, bit);
#endif
			}
			data=plane_end;
			continue;
		}
#ifdef ABAC2_CONF_MSB_RELATION
		int prevbit0=0;
#endif
		int hitratio_sure=(int)(0x10000*pow((double)hitcount/imsize, 1/BOOST_POWER)), hitratio_notsure=(int)(0x10000*pow((double)hitcount/imsize, BOOST_POWER));
		int hitratio_delta=hitratio_sure-hitratio_notsure;
		//int hitratio_sure=int(0x10000*cbrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)imsize*imsize*imsize));
		//int hitratio_sure=int(0x10000*sqrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)imsize*imsize));
		hitcount=(hitcount<<16)/imsize;
		//hitcount=unsigned(((unsigned long long)hitcount<<16)/imsize);
		//hitcount=abac2_normalize16(hitcount, logimsize);
		//hitcount*=invimsize;
#endif

		unsigned code, start=0;
		unsigned long long range=0xFFFFFFFF;

		data+=sizeof(int);
		if(data>plane_end)
			FAIL(void*, size_t, "Decode error: ptr %d > plane %d ~ %d, kp %d before loop", PTR2OFFSET(data), PTR2OFFSET(plane_end-ncodes), PTR2OFFSET(plane_end), kp);
		code=data[-4]<<24|data[-3]<<16|data[-2]<<8|data[-1];

		//printf("kp %d offset %d code %08X\n", kp, PTR2OFFSET(data), code);//

		for(int kb=0, kb2=0;kb<imsize;kb2+=bytestride)//bit-pixel loop
		{
			if(range<3)
			{
				data+=sizeof(int);
				if(data>plane_end)
					FAIL(void*, size_t, "Decode error: ptr %d > plane %d ~ %d, kp %d kb %d / %d", PTR2OFFSET(data), PTR2OFFSET(plane_end-ncodes), PTR2OFFSET(plane_end), kp, kb, imsize);
				code=data[-4]<<24|data[-3]<<16|data[-2]<<8|data[-1];
				start=0, range=0xFFFFFFFF;//because 1=0.9999...
			}
#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit=0;
			if(kp+1<bitdepth)
				prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
		//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
			//if(kp==7&&kb==5)
			//	int LOL_1=0;
			int p0=prob-PROB_HALF;
			p0=(long long)p0*prob_correct>>16;
			p0=(long long)p0*p0>>16;
			int sure=-(prevbit==prevbit0)&-(kp!=bitdepth-1);
			p0=(long long)p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
			//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
			//p0=(long long)p0*hitcount>>16;
			p0+=PROB_HALF;
			//if(prevbit!=prevbit0)
			//	p0=PROB_HALF;
			//	p0=0xFFFF-p0;

			//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

			//unsigned p0=(long long)(prob-PROB_HALF)*sqrthitcount>>16;
			//if(prevbit==prevbit0)
			//	p0=(long long)p0*hitcount>>16;
			//p0+=PROB_HALF;

			//int confboost=prevbit==prevbit0;
			//confboost-=!confboost;
			//confboost<<=LOG_CONFBOOST;
			//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(hitcount+confboost)>>16);

		//	unsigned p0=PROB_HALF+(int)((prob-PROB_HALF)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/imsize):(double)test_conf[kp]*test_conf[kp]/((double)imsize*imsize)));
		//	unsigned p0=prevbit==prevbit0?prob:PROB_HALF;
		//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*test_conf[kp]/imsize;
		//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
			PROB_CLAMP(p0);
			unsigned r2=(unsigned)(range*p0>>16);
			r2+=(r2==0)-(r2==range);
			unsigned middle=start+r2;
			int bit=code>middle;
#ifdef DEBUG_ABAC2
			if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
				printf("%6d %6d %d %08X+%08X %04X %08X %08X %08X\n", kp, kb, bit, start, (int)range, p0, r2, middle, code);
#endif
			
			int correct=bit^(p0>=PROB_HALF);
		//	hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
			prevbit0=prevbit;
#endif
			
#ifdef PRINT_ABAC
			auto start0=start;
			auto r0=range;
#endif
			if(bit)
			{
				++r2;
				start+=r2, range-=r2;
			}
			//	start=middle+1;
			else
				range=r2-1;
			//	end=middle-1;
#ifdef PRINT_ABAC
			if(kp==0)
				printf("%d %d bit %d p0 %04X %08X+%08X -> %08X+%08X c %08X\n", kp, kb, bit, p0, start0, r0, start, range, code);
#endif
			
			buffer[kb2+bit_offset]|=bit<<bit_shift;
		//	buffer[kb]|=bit<<kp;
#ifdef ENABLE_GUIDE
			int original_bit=((unsigned char*)guide)[kb2+bit_offset]>>bit_shift&1;
			if(bit!=original_bit)
				FAIL(void*, size_t, "Decode error (AC): data %d expected %d got %d", kp, original_bit, bit);
#endif
			++kb;
			
			while((start^(start+(unsigned)range))<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range 0x%08X byte-out 0x%02X\n", (int)range, code>>24);
#endif
				++data;
				if(data>plane_end)
					FAIL(void*, size_t, "Decode error: ptr %d > plane %d ~ %d, kp %d kb %d / %d", PTR2OFFSET(data), PTR2OFFSET(plane_end-ncodes), PTR2OFFSET(plane_end), kp, kb, imsize);
				code=code<<8|data[-1];
				start<<=8;
				range=range<<8|0xFF;

				//printf("[%d/%d]%02X %d-", PTR2OFFSET(data), PTR2OFFSET(plane_end), data[-1]&0xFF, kb);//
				//if(PTR2OFFSET(data)==160)//
				//	printf("kp %d kb %d code %08X\n", kp, kb, code);//
			}
		}
		if(data!=plane_end)
			FAIL(void*, size_t, "Decode error: end of plane: ptr %d, plane %d ~ %d, kp %d", PTR2OFFSET(data), PTR2OFFSET(plane_end-ncodes), PTR2OFFSET(plane_end), kp);
		//data+=ncodes;
		//cusize+=ncodes;
	}
	t2=__rdtsc();
#undef	PTR2OFFSET

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(imsize*bitdepth>>3));
	return data;
}