#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/stat.h>
uint8_t* file_load(const char *srcfn, ptrdiff_t *ret_size)
{
	struct stat info={0};
	int error=stat(srcfn, &info);
	if(error)
	{
		printf("Cannot stat \"%s\"", srcfn);
		return 0;
	}
	ptrdiff_t size=info.st_size;
	FILE *fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"", srcfn);
		return 0;
	}
	uint8_t *buf=(uint8_t*)malloc(size+16);
	if(!buf)
	{
		printf("Alloc error\n");
		return 0;
	}
	ptrdiff_t nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
	{
		printf("Read error\n");
		return 0;
	}
	fclose(fsrc);
	*ret_size=size;
	return buf;
}
int main(int argc, char **argv)
{
	if(argc!=3)
	{
		printf("Usage:  \"%s\"  file1.wav  file2.wav\n", argv[0]);
		return 1;
	}
	const char *fn1=argv[1], *fn2=argv[2];
	ptrdiff_t len1=0, len2=0;
	uint16_t *src1=(uint16_t*)file_load(fn1, &len1);
	uint16_t *src2=(uint16_t*)file_load(fn2, &len2);
	ptrdiff_t len=(len1<len2?len1:len2)/2;
	int64_t besterror=0;
	int bestshift=0;
	for(int ko=-256, it=0;ko<1024;++ko, ++it)
	{
		uint64_t error=0;
		for(int k=0;k<len;++k)
		{
			if((unsigned)(k+ko)<len)
			{
				int diff=(int16_t)src1[k]-(int16_t)src2[k+ko];
				error+=(uint64_t)diff*diff;
				//if((int64_t)error<0)//
				//{
				//	printf("diff %d\n", diff);
				//	printf("error %lld\n", error);
				//	printf("k %d\n", k);
				//	return 0;
				//}
			}
		}
		//printf("%lld\n", error);//
		if(!it||besterror>error)
		{
			//printf("^\n");//
			besterror=error;
			bestshift=ko;
		}
	}
	double rmse=sqrt((double)besterror/len);
	double psnr=-20*log10(rmse/0x7FFF);
	printf("RMSE %12.6lf  PSNR %12.6lf  phase shift %d\n", rmse, psnr, bestshift);
	free(src1);
	free(src2);
	return 0;
}