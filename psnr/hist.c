#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<sys/stat.h>
int hist1[0x10000]={0};
int hist0[256]={0};
void print_usage(const char *argv0)
{
	printf(
		"Usage:  %s  file  [2]\n"
		"  Optionally pass [2] for uint16 o0 entropy\n"
		, argv0
	);
}
int main(int argc, char **argv)
{
	if(argc!=2&&argc!=3)
	{
		print_usage(argv[0]);
		return 1;
	}
	const char *fn=argv[1];
	ptrdiff_t len=0;
	unsigned char *buf=0;

	int twobytes=0;
	if(argc==3)
	{
		if(!(argv[2][0]=='2'&&!argv[2][1]))
		{
			print_usage(argv[0]);
			return 1;
		}
		twobytes=1;
	}

	//load file
	{
		struct stat info={0};
		if(stat(fn, &info))
		{
			printf("Cannot open \"%s\"\n", fn);
			return 1;
		}
		len=info.st_size;
		if(!len)
		{
			printf("Empty file \"%s\"\n", fn);
			return 1;
		}
		FILE *fsrc=fopen(fn, "rb");
		if(!fsrc)
		{
			printf("Cannot open \"%s\"\n", fn);
			return 1;
		}
		buf=(unsigned char*)malloc(len+16);
		if(!buf)
		{
			fclose(fsrc);
			printf("Alloc error\n");
			return 1;
		}
		fread(buf, 1, len, fsrc);
		fclose(fsrc);
	}

	if(twobytes)
	{
		//accumulate histogram
		unsigned short *ptr=(unsigned short*)buf;
		len>>=1;
		for(int k=0;k<len;++k)
			++hist1[*ptr++];

		//calculate order-0 entropy
		double e0=0, norm=1./len;
		for(int ks=0;ks<0x10000;++ks)
		{
			int freq=hist1[ks];
			if(freq)
				e0-=freq*log2(freq*norm);
		}
		e0/=8;
		len<<=1;
		printf("Order-0 entropy: %12.2lf / %9td bytes  %12.6lf%%\n", e0, len, e0/len*100);
	}
	else
	{
		//accumulate histogram
		int prev=0;
		for(int k=0;k<len;++k)
		{
			unsigned char curr=buf[k];
			++hist0[curr];
			++hist1[prev<<8|curr];
			prev=buf[k];
		}

		//get peak
		int vmax=0;
		for(int ks=0;ks<256;++ks)
		{
			if(vmax<hist0[ks])
				vmax=hist0[ks];
		}

		//print histogram
		for(int ks=0;ks<256;++ks)
		{
			int freq=hist0[ks];
			printf("%3d %02X %8d", ks, ks, freq);
			if(freq)
			{
				int nstars=freq*96/vmax;
				printf(" ");
				for(int k=0;k<nstars;++k)
					printf("*");
			}
			printf("\n");
		}

		//calculate order-0 entropy
		double e0=0, norm=1./len;
		for(int ks=0;ks<256;++ks)
		{
			int freq=hist0[ks];
			if(freq)
				e0-=freq*log2(freq*norm);
		}
		e0/=8;
		printf("Order-0 entropy: %12.2lf / %9td bytes  %12.6lf%%\n", e0, len, e0/len*100);

		//calculate order-1 entropy
		double e1=0;
		for(int kW=0;kW<256;++kW)
		{
			int *curr_hist=hist1+((size_t)kW<<8);
			int sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=curr_hist[ks];
			if(sum)
			{
				norm=1./sum;
				for(int ks=0;ks<256;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
						e1-=freq*log2(freq*norm);
				}
			}
		}
		e1/=8;
		printf("Order-1 entropy: %12.2lf / %9td bytes  %12.6lf%%\n", e1, len, e1/len*100);
	}
	free(buf);
	return 0;
}