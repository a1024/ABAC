#include "e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

#if 0
void test1()
{
	//double w0=2.95, w1=3.05;
	double w0=4.5, w1=1;
	double lr=0.05;
	for(int it=0;it<100;++it)
	{
		double grad_w0=0, grad_w1=0;
		printf("iter %3d: w0=%g\tw1=%g\tmissed:", it, w0, w1);
		for(int n=-10;n<0;++n)
		{
			double c=w1*n+w0;
			if(c>0)//check negative integers
			{
				printf(" %d", n);
				++grad_w0;//delta is +ve for -ve numbers
				grad_w1+=n;
			}
		}
		for(int n=0;n<10;++n)
		{
			double c=w1*n+w0;
			if(c<0)//check positive integers
			{
				printf(" %d", n);
				--grad_w0;//delta is -ve for +ve numbers
				grad_w1-=n;
			}
		}
		if(!grad_w0)
			printf(" None");
		printf("\tgrad_w0=%g\tgrad_w1=%g", grad_w0, grad_w1);
		printf("\n");
		w0-=lr*grad_w0;
		w1-=lr*grad_w1;
	}
	printf("Done.\n");
	pause();
	exit(0);
}
#endif

int hist[256*4];

short jxl_Kodak_params[]=
{
	0x06E3,  0x1287,  0x1192,  0x0C0B,	 0x0007,  0x000E, -0x0180, -0x006B, -0x00DB, -0x0052,  0x0126,//01
	0x09B9,  0x0E20,  0x185F,  0x0BF3,	-0x0031, -0x004F,  0x003F, -0x003F,  0x009B, -0x0030, -0x0064,
	0x040B,  0x0F22,  0x0D82,  0x0BF8,	-0x0006,  0x0013, -0x011E, -0x0019, -0x00C6, -0x0027,  0x011A,

	0x04DF,  0x1309,  0x0EF3,  0x0C0B,	-0x000C, -0x0008, -0x00F5, -0x007B, -0x00DB, -0x003A,  0x0124,//02
	0x0835,  0x0F21,  0x185F,  0x0BF1,	-0x0031, -0x002F,  0x0073, -0x003F,  0x00C7, -0x0004, -0x00E3,
	0x009B,  0x0F22,  0x0A82,  0x0BA8,	 0x0008,  0x0006, -0x00BE,  0x001C, -0x00C6, -0x000B,  0x011A,
  
	0x04C0,  0x1369,  0x0FFB,  0x0B8D,	 0x001F,  0x0003, -0x01F5,  0x0007, -0x0103, -0x0039,  0x0124,//03
	0x0837,  0x0E21,  0x195F,  0x0BED,	-0x0040, -0x0067,  0x0053,  0x0181,  0x0133, -0x0004, -0x00D9,
	0x01CB,  0x0EE2,  0x0B96,  0x0981,	 0x0003,  0x0002, -0x015E,  0x006C, -0x00E8, -0x000A,  0x0119,
  
	0x04CF,  0x1367,  0x0CB0,  0x0CD1,	-0x0004,  0x0002, -0x00E5, -0x00C9, -0x0089, -0x000F,  0x012C,//04
	0x0876,  0x0FA5,  0x16E1,  0x0BEE,	-0x0078, -0x003F,  0x0053,  0x0028,  0x0037,  0x0054, -0x00C9,
	0x0177,  0x0ED2,  0x0A94,  0x0A43,	 0x0002,  0x0002, -0x015C, -0x0034, -0x00E0,  0x001E,  0x0138,
 
	0x04CF,  0x1361,  0x0CEA,  0x0AD2,	-0x0002, -0x0002, -0x00E5, -0x0069, -0x0141, -0x000B,  0x012C,//05
	0x0877,  0x0E65,  0x16E1,  0x0FD3,	-0x0090, -0x005F,  0x0027,  0x002A, -0x0013,  0x015A, -0x0085,
	0x01F7,  0x0E31,  0x0B94,  0x09AF,	-0x0016,  0x0002, -0x0057, -0x0013, -0x00D7,  0x0020,  0x0106,

	0x04CF,  0x134E,  0x0DAA,  0x0C0B,	 0x001E, -0x0006, -0x016F,  0x0060, -0x0161, -0x000B,  0x016E,//06
	0x08B7,  0x0E65,  0x1720,  0x0D2F,	-0x00AA, -0x0037,  0x0047,  0x007A,  0x011E, -0x004D,  0x007C,
	0x01F7,  0x0D31,  0x0C90,  0x09AF,	-0x0008,  0x0002, -0x00C0,  0x0062, -0x00F7,  0x000C,  0x0146,

	0x04CF,  0x134E,  0x0D2A,  0x0C4A,	 0x0002,  0x0002, -0x01AD, -0x0020, -0x0149, -0x001B,  0x016E,//07
	0x08F4,  0x0F65,  0x15A0,  0x0D2F,	-0x008A, -0x0013, -0x0019,  0x0034, -0x001A,  0x0131, -0x0019,
	0x0261,  0x0D30,  0x09AC,  0x0B36,	-0x0004,  0x0002, -0x00B7, -0x0004, -0x00B6,  0x000C,  0x0126,
 
	0x03CF,  0x1370,  0x0922,  0x0CDA,	 0x0003,  0x0004, -0x00DB, -0x002E, -0x00C9, -0x0043,  0x013F,//08
	0x07F4,  0x10E7,  0x159F,  0x0D2F,	-0x000A, -0x005B,  0x0081, -0x0104, -0x001A, -0x0062,  0x00F3,
	0x0260,  0x0E30,  0x07AC,  0x0B36,	-0x0004, -0x000E, -0x0083, -0x000B, -0x0092, -0x0021,  0x0126,

	0x03F2,  0x1370,  0x0913,  0x0C60,	 0x0002,  0x0024, -0x01DC,  0x003A, -0x0109, -0x0043,  0x014B,//09
	0x0877,  0x0F67,  0x157D,  0x0AD0,	-0x0022, -0x003F,  0x0087, -0x0080, -0x0012, -0x0062, -0x011F,
	0x02C4,  0x0E20,  0x07A4,  0x09A5,	-0x0004,  0x002A, -0x0181, -0x000F, -0x0112, -0x0021,  0x0155,

	0x0462,  0x1372,  0x07CD,  0x0B1E,	 0x0002,  0x0014, -0x020C,  0x0032, -0x012B, -0x0044,  0x0183,//10
	0x0A35,  0x0EA7,  0x13BC,  0x0BCF,	-0x0017, -0x003D,  0x0081, -0x002E, -0x0086,  0x00C6, -0x006A,
	0x03C4,  0x0F39,  0x0646,  0x0995,	-0x000B,  0x002A, -0x017E, -0x000B, -0x00CC, -0x0011,  0x0175,

	0x045C,  0x1336,  0x084E,  0x0C1D,	 0x0002,  0x0003, -0x01EC, -0x004A, -0x012B, -0x0044,  0x0233,//11
	0x0966,  0x0E9D,  0x147C,  0x084F,	-0x0058, -0x0032,  0x0081,  0x0007,  0x00DA,  0x0004, -0x0042,
	0x038C,  0x0F3F,  0x08A7,  0x0993,	-0x000B,  0x0002, -0x016A, -0x000B, -0x00D4, -0x000F,  0x0175,
  
	0x02EA,  0x1376,  0x099D,  0x0C0C,	 0x0003,  0x0011, -0x02E0,  0x006E, -0x012B, -0x0044,  0x020B,//12
	0x095B,  0x0E9D,  0x157B,  0x084B,	-0x0058, -0x002A, -0x006F,  0x0027,  0x015A, -0x0008, -0x0002,
	0x026E,  0x0EFF,  0x08B3,  0x09A0,	 0x0003,  0x001A, -0x01EA,  0x0071, -0x0154, -0x0007,  0x0174,

	0x03AA,  0x1377,  0x0A8F,  0x0A08,	 0x0003,  0x000D, -0x01DF,  0x0026, -0x01B9, -0x0028,  0x01FF,//13
	0x095B,  0x0E51,  0x157C,  0x0954,	-0x0080, -0x002A,  0x00A7,  0x0039,  0x015A,  0x001D, -0x013A,
	0x02EE,  0x0E7F,  0x08B4,  0x0693,	-0x001C, -0x0006, -0x012A, -0x0051, -0x00EC, -0x001F,  0x01DF,

	0x03CA,  0x11DB,  0x0A9F,  0x09F8,	 0x0003,  0x0001, -0x01DF,  0x002A, -0x017C, -0x002C,  0x01FE,//14
	0x0717,  0x0E51,  0x1585,  0x06E4,	-0x009F, -0x004A,  0x0054,  0x0039,  0x011A,  0x006B,  0x0006,
	0x02D4,  0x0E7F,  0x06E5,  0x09D3,	-0x0004,  0x0002, -0x0120, -0x0021, -0x0028,  0x0023,  0x0261,

	0x03CA,  0x11DB,  0x079D,  0x09F9,	 0x0002,  0x000B, -0x001B, -0x00D3, -0x00AD,  0x0030,  0x0116,//15
	0x0717,  0x0E53,  0x1543,  0x0A70,	 0x0004, -0x0002,  0x0034, -0x00F3, -0x000A,  0x00EB, -0x003B,
	0x02E2,  0x0E76,  0x0449,  0x08D3,	 0x0002,  0x0007, -0x0012, -0x0023, -0x007C,  0x0063,  0x0105,

	0x04CA,  0x115B,  0x09DD,  0x09F8,	 0x002E,  0x0002, -0x01A5,  0x0068, -0x014D, -0x0016,  0x0146,//16
	0x0717,  0x0E53,  0x1543,  0x0A6C,	-0x002C, -0x005E,  0x01B3,  0x0227,  0x0076,  0x007E,  0x014B,
	0x02E7,  0x0E16,  0x0857,  0x08CC,	-0x0006, -0x0002, -0x0126,  0x00D5, -0x00BD,  0x002F,  0x019D,
 
	0x070A,  0x115B,  0x08ED,  0x09F8,	 0x001A,  0x0010, -0x01E6, -0x0016, -0x012D, -0x0046,  0x00D6,//17
	0x07BD,  0x0E55,  0x12C3,  0x0A6C,	-0x0052, -0x0085,  0x01B2,  0x01E5,  0x006A,  0x0102, -0x01AF,
	0x036B,  0x0E56,  0x0853,  0x06AA,	-0x0008, -0x0024, -0x0106,  0x0053, -0x013D, -0x000A,  0x010A,
 
	0x06C6,  0x115B,  0x08ED,  0x08F4,	-0x000C, -0x0010, -0x0161, -0x0116, -0x00E9, -0x007E,  0x0154,//18
	0x08B1,  0x11D2,  0x12C0,  0x09EE,	-0x006A, -0x008C,  0x011A,  0x00DD,  0x00DA,  0x007E, -0x019F,
	0x0419,  0x0F56,  0x0803,  0x06AB,	-0x0026, -0x0020, -0x00EA, -0x0065, -0x00A0, -0x004A,  0x0148,
 
	0x0644,  0x1154,  0x08F0,  0x08F0,	 0x0004, -0x0008, -0x0165, -0x0016, -0x00EB, -0x007E,  0x0114,//19
	0x08AB,  0x1514,  0x12C0,  0x0788,	-0x003E, -0x0041,  0x00D2,  0x00B0,  0x00DA,  0x007C, -0x00CF,
	0x02FB,  0x0ED4,  0x0810,  0x06E7,	-0x0010, -0x0010, -0x012E,  0x0023, -0x00F0, -0x003A,  0x0144,

	0x05C4,  0x1153,  0x09F1,  0x077C,	 0x0008,  0x0006, -0x0165, -0x0018, -0x00D3, -0x003E,  0x011C,//20
	0x08AF,  0x1314,  0x141A,  0x08C7,	-0x0037, -0x0001,  0x00E6,  0x01CE,  0x00A2,  0x0136, -0x00B7,
	0x0501,  0x0DD6,  0x091D,  0x05E9,	 0x001E,  0x0002, -0x009E,  0x002C, -0x00CC,  0x0026,  0x0114,

	0x0584,  0x0E73,  0x09F1,  0x0884,	 0x0012,  0x0012, -0x01C8, -0x0008, -0x0114, -0x003B,  0x0121,//21
	0x08CF,  0x1310,  0x140A,  0x0437,	-0x0044, -0x0019,  0x00A6,  0x01CF,  0x0081,  0x0122,  0x0050,
	0x0488,  0x0B57,  0x093F,  0x06E1,	-0x001F, -0x0006, -0x010E, -0x0074, -0x00A4, -0x003A,  0x0128,

	0x0581,  0x0E6B,  0x08CF,  0x07E4,	 0x0003, -0x0012, -0x0146, -0x00E4, -0x00F8, -0x007B,  0x011F,//22
	0x08CF,  0x1498,  0x109A,  0x045D,	-0x003E, -0x0045,  0x00A6,  0x007F,  0x00A1,  0x0112, -0x00AA,
	0x0388,  0x0B97,  0x07BC,  0x06E0,	-0x000F, -0x0012, -0x00CE, -0x0081, -0x00A4, -0x003A,  0x0120,

	0x0589,  0x0EED,  0x0687,  0x0765,	 0x0003,  0x0001, -0x012E, -0x00B4, -0x00F8, -0x0035,  0x0119,//23
	0x08D0,  0x1461,  0x0E8A,  0x077F,	-0x003E, -0x0035,  0x0061,  0x018F,  0x004D,  0x017A, -0x006A,
	0x03F8,  0x0CCE,  0x06B4,  0x0700,	-0x0005, -0x0002, -0x010B, -0x0081, -0x00B5, -0x0038,  0x011D,

	0x045A,  0x0FF0,  0x0626,  0x0765,	 0x0003,  0x0001, -0x0277, -0x00AC, -0x00FC, -0x008B,  0x0199,//24
	0x090F,  0x147F,  0x0E8A,  0x063E,	-0x0046, -0x0059,  0x0059,  0x000F,  0x00DF,  0x0142,  0x0006,
	0x02F8,  0x0CD1,  0x06A4,  0x0700,	-0x0005, -0x0012, -0x010B, -0x007F, -0x00B5, -0x0030,  0x015E,
};
short jxl_CLIC30_params[]=
{
	  0x0B37,  0x130A,  0x121B,  0x0B7A,	-0x0023, -0x0032, -0x00AE, -0x0093, -0x0099,  0x0038,  0x00BF,//01
	  0x0D37,  0x1315,  0x181F,  0x0BF4,	-0x005C, -0x0073,  0x00A1, -0x0131,  0x00BD,  0x0026, -0x00BA,
	  0x06C2,  0x0EF2,  0x0C43,  0x0C14,	-0x0026, -0x001E, -0x0085, -0x0017, -0x0082,  0x0027,  0x00D4,

	  0x0BB6,  0x134A,  0x141B,  0x0AFA,	-0x0002, -0x00EA,  0x00D9, -0x0213, -0x0038,  0x0048,  0x0125,//02
	  0x0C37,  0x1518,  0x179F,  0x0BF4,	-0x0098, -0x0060,  0x0010,  0x0011,  0x0079,  0x00F0,  0x0001,
	  0x06C0,  0x0F02,  0x0D43,  0x0A12,	-0x0079, -0x0002, -0x0105, -0x0197,  0x00FF,  0x00A5,  0x013A,

	  0x0BB6,  0x134A,  0x141B,  0x0EFA,	-0x0007, -0x01DE,  0x00D8,  0x008D,  0x0006,  0x0041,  0x011E,//03
	  0x0BB7,  0x1518,  0x179F,  0x0D54,	-0x0007, -0x0330,  0x0021,  0x0004,  0x0008,  0x000A,  0x0107,
	  0x06C0,  0x0F04,  0x08BB,  0x0A12,	-0x0034,  0x0002, -0x0085, -0x0197,  0x01B9,  0x00A5,  0x0305,
 
	  0x0BB6,  0x13CA,  0x131B,  0x0DBA,	-0x001B, -0x003F,  0x00B0,  0x010D,  0x0010,  0x0021, -0x006B,//04
	  0x0BB7,  0x1618,  0x179F,  0x0D54,	-0x000F, -0x0033,  0x0070,  0x00B4, -0x0042,  0x000A, -0x007B,
	  0x06C4,  0x0F04,  0x08BB,  0x03CE,	 0x000F,  0x0006, -0x0085, -0x0197,  0x0176,  0x00A5,  0x0305,
  
	  0x0BB6,  0x13CA,  0x131B,  0x0DBA,	-0x001B, -0x003F,  0x00B0,  0x010D,  0x0010,  0x0021, -0x006B,//05
	  0x09B7,  0x1618,  0x179F,  0x10D4,	 0x0018,  0x0024,  0x0090,  0x016D, -0x001A,  0x0082, -0x007A,
	  0x06C4,  0x0F04,  0x08BB,  0x03CE,	 0x000F,  0x0006, -0x0085, -0x0197,  0x0176,  0x00A5,  0x0305,
  
	  0x0A97,  0x1208,  0x1363,  0x0FBB,	-0x0023, -0x001F,  0x00E8,  0x015D,  0x003F,  0x0080, -0x00EB,//06
	  0x08A8,  0x1618,  0x179E,  0x0F54,	-0x0020, -0x0012,  0x00A1,  0x016D, -0x0014,  0x00AD, -0x0086,
	  0x06C4,  0x0EF1,  0x0A3B, -0x00EC,	 0x0018, -0x0001, -0x0099, -0x0196,  0x0178,  0x007D,  0x0407,
  
	  0x08D7,  0x1388,  0x13C3,  0x0EBB,	-0x0029, -0x0013,  0x0088,  0x015D,  0x000E,  0x00AD, -0x0037,//07
	  0x09A8,  0x1618,  0x179E,  0x0C53,	-0x0080, -0x00CD,  0x009F,  0x0158,  0x005E,  0x007D,  0x0027,
	  0x06D2,  0x0E67,  0x0A17, -0x01AD,	 0x0019, -0x0019, -0x0119, -0x0196,  0x01F8, -0x0153,  0x0287,
  
	  0x0CD7,  0x1388,  0x1343,  0x0EBB,	-0x0069, -0x00BF,  0x0033,  0x0089, -0x002F,  0x00AF,  0x0049,//08
	  0x0A08,  0x1518,  0x179E,  0x0ED3,	-0x003E,  0x0013,  0x00A3,  0x0209,  0x005E,  0x0160, -0x00D9,
	  0x06DA,  0x0EA7,  0x0A17, -0x02B1,	-0x0040, -0x00BA, -0x0119, -0x0195,  0x01C7, -0x0573,  0x02A7,
 
	  0x0CD9,  0x164C,  0x12C3,  0x0EBB,	-0x007B, -0x003F,  0x002C,  0x00EB,  0x0039,  0x0096,  0x0045,//09
	  0x0980,  0x13F8,  0x17E2,  0x0ED4,	 0x0162, -0x003C, -0x0004,  0x01C1,  0x0017,  0x0160, -0x0339,
	  0x01CA,  0x0EA7,  0x0A1D, -0x01EF,	 0x0004,  0x0001, -0x011A, -0x0195,  0x0147, -0x0573,  0x02A7,
  
	  0x0C19,  0x1A0C,  0x1243,  0x0EBB,	-0x005E, -0x009F, -0x0008,  0x00AF,  0x0031,  0x003D,  0x00B3,//10
	  0x0980,  0x13F8,  0x17E2,  0x0ED4,	-0x007C, -0x015B,  0x01E8,  0x019C,  0x0098,  0x0117, -0x013D,
	  0x04C6,  0x0EA7,  0x0A1D,  0x0097,	-0x003D, -0x0080, -0x011A, -0x0195,  0x0147, -0x0573,  0x02A7,
 
	  0x0898,  0x1A0C,  0x1246,  0x0EBB,	-0x0060, -0x0025, -0x00AF,  0x003E, -0x002F,  0x003D,  0x012F,//11
	  0x0981,  0x13F8,  0x17E6,  0x09A9,	-0x005D, -0x008B,  0x0167,  0x019C,  0x0142,  0x0117, -0x0119,
	  0x0472,  0x0BB8,  0x0B5E,  0x003B,	-0x0058, -0x0051, -0x0122, -0x0195,  0x0147, -0x056D,  0x00E7,
	  
	  0x0708,  0x1A0C,  0x1246,  0x1003,	-0x002F, -0x0033, -0x003F, -0x008C,  0x0046,  0x0021,  0x01AF,//12
	  0x0A04,  0x1458,  0x14A4,  0x0869,	-0x0035, -0x004D,  0x0148,  0x01B8,  0x0032,  0x0109, -0x00E3,
	  0x0232,  0x0BBC,  0x0B56,  0x0182,	-0x0034, -0x0081, -0x0121, -0x01D5,  0x0145, -0x0515,  0x045F,
	  
	  0x0888,  0x1A0C,  0x14C6,  0x0FC2,	-0x0028,  0x0013, -0x007B,  0x0098, -0x004E,  0x002D,  0x00C7,//13
	  0x0B97,  0x1462,  0x149E,  0x0869,	-0x0015, -0x0023,  0x0168,  0x0228,  0x0033,  0x0107, -0x0127,
	  0x04D2,  0x0BC0,  0x0B3C, -0x0176,	-0x001E,  0x000E, -0x0019, -0x0259,  0x0155, -0x0514,  0x045F,
	  
	  0x087C,  0x1A19,  0x14C5,  0x101F,	-0x0010, -0x0009, -0x007B,  0x0018, -0x00AE,  0x002F,  0x00C6,//14
	  0x0BD7,  0x1523,  0x1116,  0x0A3C,	-0x0027, -0x000F,  0x0167,  0x0228,  0x0058,  0x00FD, -0x0107,
	  0x04D3,  0x0C8D,  0x0A3A, -0x015E,	-0x0036, -0x0013,  0x02E7, -0x0239,  0x0364, -0x0526,  0x041D,
	  
	  0x08E6,  0x1A1A,  0x15C5,  0x0BF2,	-0x0010, -0x0009, -0x0059,  0x002A, -0x004D,  0x002E,  0x0134,//15
	  0x0BD5,  0x1522,  0x1116,  0x06AC,	-0x0022, -0x0057,  0x0166,  0x01AC,  0x00F1,  0x00BD,  0x0139,
	  0x04A7,  0x0C8B,  0x0DAC, -0x0069,	-0x0038, -0x0034,  0x0327, -0x0229,  0x04E8, -0x0520,  0x0447,
	  
	  0x09E8,  0x1A1A,  0x15C5,  0x0BF2,	-0x00B3, -0x0026, -0x0055,  0x00AA,  0x00D7,  0x00AA,  0x0006,//16
	  0x09D5,  0x1522,  0x1217,  0x067C,	-0x0025, -0x00F2,  0x03F6,  0x01AC,  0x00D5,  0x0107,  0x0178,
	  0x0470,  0x0DD5,  0x0B92, -0x00AA,	-0x004B, -0x00B4,  0x0227, -0x0217,  0x04F0, -0x05A0,  0x04BF,
	  
	  0x0D68,  0x199A,  0x15C5,  0x0BF2,	-0x0052, -0x0010,  0x00A0,  0x0126,  0x00D6,  0x00D1, -0x0078,//17
	  0x09D5,  0x1522,  0x1117,  0x067C,	-0x0015,  0x0004,  0x0296,  0x0124,  0x00D5,  0x046B,  0x0112,
	  0x0392,  0x0E01,  0x0A72, -0x0026,	-0x002E, -0x001D,  0x0227, -0x0207,  0x0470, -0x0580,  0x044B,
	  
	  0x0DA8,  0x1ADA,  0x15DC,  0x0CE2,	-0x0052, -0x0002,  0x0000,  0x0126,  0x01D5,  0x017D, -0x0058,//18
	  0x09D5,  0x1523,  0x1119,  0x07BC,	-0x003D, -0x0044,  0x02D6,  0x01A4,  0x0218,  0x023D, -0x00BC,
	  0x02D2,  0x0DFF,  0x0BC2, -0x0225,	-0x0002, -0x000C, -0x0041, -0x0207,  0x0470, -0x0580,  0x03EB,
	  
	  0x0F08,  0x1ADA,  0x139C,  0x0C38,	-0x0060,  0x000D,  0x0102,  0x01C6,  0x0131,  0x0169, -0x009A,//19
	  0x09D3,  0x1525,  0x1118,  0x07CC,	-0x0022,  0x0028,  0x0256,  0x03A9,  0x00D1,  0x023D, -0x0214,
	  0x045B,  0x0CE7,  0x0B62, -0x01AD,	-0x003B,  0x0006, -0x00C0, -0x0288,  0x0356, -0x0580,  0x051F,
	  
	  0x0F08,  0x1BBA,  0x139C,  0x09B8,	-0x0027, -0x0072,  0x0099,  0x0062,  0x0119,  0x0064,  0x0006,//20
	  0x09CD,  0x163D,  0x0FE9,  0x0577,	-0x001D, -0x0036,  0x024E,  0x03A5,  0x0155,  0x016D, -0x01D3,
	  0x040D,  0x0C95,  0x0849,  0x0001,	-0x0016, -0x005D, -0x00C0, -0x0288,  0x0355, -0x0580,  0x051F,
	  
	  0x0E88,  0x1DBB,  0x129B,  0x09BA,	-0x0063, -0x0048, -0x0089, -0x000E,  0x0099,  0x0102,  0x0042,//21
	  0x0838,  0x163D,  0x0FE9,  0x03F7,	-0x0033, -0x004E,  0x024E,  0x00C8,  0x0259,  0x016B, -0x01D3,
	  0x04B3,  0x0D95,  0x086A, -0x00E0,	-0x0041, -0x00DB, -0x00B8, -0x0289,  0x00D5, -0x0620,  0x051F,
	  
	  0x0E88,  0x1D9B,  0x172B,  0x09BB,	-0x0047, -0x0022, -0x0087,  0x00B2,  0x0099,  0x003E,  0x00E5,//22
	  0x0838,  0x15F5,  0x134E,  0x03F7,	-0x0023, -0x008F,  0x00EA,  0x027D,  0x0257,  0x016B, -0x00F3,
	  0x04F2,  0x0D57,  0x0CE7, -0x0050,	-0x002F, -0x0055, -0x00B8, -0x0289,  0x00D5, -0x0640,  0x051F,
	  
	  0x0E8D,  0x1E19,  0x14EB,  0x0ADB,	-0x0005, -0x000C, -0x0061,  0x0098, -0x002B,  0x003E,  0x00F8,//23
	  0x0838,  0x15F5,  0x134C,  0x06F7,	 0x000D,  0x0049,  0x00EA,  0x027E,  0x0015,  0x00FB, -0x018A,
	  0x0672,  0x0D97,  0x0B17, -0x013E,	-0x002A, -0x000C, -0x0038, -0x0289,  0x0215, -0x05C0,  0x0563,
	  
	  0x0E8D,  0x1CD9,  0x14EB,  0x0E1B,	-0x0039,  0x000E,  0x003F,  0x0198,  0x005F,  0x00A1,  0x0059,//24
	  0x08FC,  0x15F6,  0x133D,  0x0AB2,	-0x0033, -0x0029,  0x00EF,  0x023E, -0x0033,  0x014B, -0x00FA,
	  0x064C,  0x0CB4,  0x0B3B, -0x0236,	-0x001B,  0x001C,  0x01CE, -0x048B,  0x01D5, -0x0604,  0x0517,
	  
	  0x0D6B,  0x1C99,  0x162B,  0x0C1C,	-0x002E, -0x0022,  0x0082,  0x0184,  0x0037,  0x0076, -0x008B,//25
	  0x09CD,  0x15F6,  0x129D,  0x0813,	-0x0027, -0x0013,  0x00D5,  0x0226,  0x003A,  0x00CB, -0x00F5,
	  0x054C,  0x0C94,  0x092B, -0x01F4,	-0x0002, -0x0003,  0x020D, -0x0493,  0x0366, -0x0664,  0x0377,
	  
	  0x0E6B,  0x1E19,  0x14AB,  0x0C1C,	-0x0069, -0x0031,  0x0066,  0x0135,  0x0037,  0x0177, -0x002B,//26
	  0x09CD,  0x15F6,  0x109D,  0x0893,	-0x001A, -0x0101,  0x0091,  0x0142, -0x00B3,  0x0181,  0x0015,
	  0x058C,  0x0C94,  0x078B, -0x0256,	-0x0046, -0x0035,  0x018D, -0x0412,  0x00E4, -0x0684,  0x0373,
	  
	  0x0E7B,  0x1DF0,  0x14AF,  0x0C1E,	-0x0003, -0x0001, -0x00CD, -0x004D, -0x003B,  0x00FB,  0x00E1,//27
	  0x08C5,  0x15F6,  0x109D,  0x0713,	-0x0001,  0x0000, -0x0040,  0x0116,  0x000F,  0x019C,  0x0097,
	  0x050D,  0x0CC0,  0x0A6D, -0x01C6,	-0x000B, -0x003D,  0x010D, -0x02D2,  0x0218, -0x0684,  0x0373,
	  
	  0x0E7B,  0x1DF0,  0x14AF,  0x0B9E,	-0x005F,  0x0007, -0x000D,  0x0137,  0x00AE,  0x011B,  0x000E,//28
	  0x08E5,  0x1576,  0x0F9D,  0x0913,	-0x004D,  0x0026, -0x0080,  0x0119,  0x0006,  0x015C,  0x0027,
	  0x03C3,  0x0CBD,  0x08ED, -0x01E6,	-0x003F, -0x0024, -0x00A3, -0x02D2, -0x0038, -0x0660,  0x0373,

	  0x0D84,  0x1DE7,  0x149F,  0x0B9A,	-0x004F,  0x0003,  0x006C,  0x0154,  0x00BC,  0x011B, -0x0072,//29
	  0x0B36,  0x1575,  0x0F1C,  0x08EF,	-0x0025,  0x0006,  0x008F,  0x0159,  0x004A,  0x0156, -0x00B9,
	  0x03E3,  0x0C7D,  0x08FD, -0x01BE,	-0x002F, -0x0017, -0x04FE, -0x029F, -0x013C, -0x065E,  0x02B0,
	
	  0x0C7C,  0x1DD7,  0x149D,  0x0D1A,	-0x0025, -0x01DA,  0x010C,  0x014C,  0x00F8,  0x010B,  0x01DF,//30
	  0x09B6,  0x1575,  0x0F1D,  0x08EF,	-0x0041, -0x00A0,  0x01AF,  0x010A,  0x010E,  0x0116,  0x0173,
	  0x01CF,  0x09EE,  0x090D,  0x0062,	-0x001D, -0x0101, -0x04FE, -0x029F, -0x013C, -0x065E,  0x02AF,
};

const char *g_extensions[]=
{
	"png",
	"jpg",
	"jpeg",
};
void batch_test(const char *path)
{
	int known_dataset=0;
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Start %s\n", g_buf);
	double t_start=time_ms();
	ArrayHandle filenames=get_filenames(path, g_extensions, COUNTOF(g_extensions), 1);
	if(!filenames)
	{
		printf("No images in \"%s\"\n", path);
		return;
	}
	long long
		count_PNG=0, count_JPEG=0,
		sum_cPNGsize=0, sum_cJPEGsize=0,
		sum_uPNGsize=0, sum_uJPEGsize=0,
		sum_testsize=0;
		//sum_test3size[2]={0};
#if 0	
	long long *hist=e10_start();//
	if(!hist)
	{
		exit(0);
		return;
	}
#endif
	{
		ArrayHandle path2=filter_path(path);
		if(!acme_stricmp(path2->data, "D:/ML/dataset-CLIC30/"))
			known_dataset=1;
		else if(!acme_stricmp(path2->data, "D:/ML/dataset-Kodak/"))
			known_dataset=2;
		array_free(&path2);
	}

	for(ptrdiff_t k=0;k<(ptrdiff_t)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);

		if(!fn)
		{
			LOG_ERROR("filename read error");
			continue;
		}

		ptrdiff_t formatsize=get_filesize(fn[0]->data);
		if(!formatsize||formatsize==-1)//skip non-images
			continue;

		int iw=0, ih=0, nch0=3, stride=4;
		long long cycles=__rdtsc();
		unsigned char *buf=image_load(fn[0]->data, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			printf("Cannot open \"%s\"\n", fn[0]->data);
			continue;
		}

		ptrdiff_t res=(ptrdiff_t)iw*ih, len=res*stride, usize=res*nch0;
		double ratio=(double)usize/formatsize;
#ifdef BATCHTEST_NO_B2
		printf("%3lld/%3lld  %.2lf%%\r", k+1, filenames->count, (k+1)*100./filenames->count);
#else
		printf("%3lld/%3lld  \"%s\"\tCR %lf (%lf BPP) Dec %lf CPB", k+1, filenames->count, fn[0]->data, ratio, 8/ratio, (double)cycles/usize);
#endif
		if(!acme_stricmp(fn[0]->data+fn[0]->count-3, "PNG"))
		{
			sum_cPNGsize+=formatsize;
			sum_uPNGsize+=usize;
			++count_PNG;
		}
		else//assumed
		{
			sum_cJPEGsize+=formatsize;
			sum_uJPEGsize+=usize;
			++count_JPEG;
		}
#ifndef BATCHTEST_NO_B2
		unsigned char *b2=(unsigned char*)malloc(len);
		if(!b2)
		{
			LOG_ERROR("Allocation error");
			return;
		}
		memset(b2, 0, len);
#endif
		
		//T34 ABAC + adaptive Bayesian inference
#if 1
		{
			ArrayHandle cdata=0;
			printf("\n");
#if 0
			if(known_dataset==1&&k<30)//pre-trained
				g_param_ptr=jxl_CLIC30_params+33*k;
			else if(known_dataset==2&&k<24)
				g_param_ptr=jxl_Kodak_params+33*k;
			else
			{
				memcpy(b2, buf, len);
				addbuf(b2, iw, ih, 3, 4, 128);
				colortransform_ycocb_fwd((char*)b2, iw, ih);
				pred_opt_opt_v6((char*)b2, iw, ih, 1);
				memset(b2, 0, len);
				pred_opt_printparam();
			}
#endif

			//printf("\nT35 (ABAC + context tree)\n");
			//printf("\nT39 Multiple estimators for all maps  WH %d*%d\n", iw, ih);
			//t35_encode(buf, iw, ih, &cdata, 1);
			t39_encode(buf, iw, ih, &cdata, 1);
			//t40_encode(buf, iw, ih, &cdata, 1);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!\n");

			//t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);
			t39_decode(cdata->data, cdata->count, iw, ih, b2, 1);
			//t40_decode(cdata->data, cdata->count, iw, ih, b2, 1);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "Test", 0);

			//printf("\nT34 (ABAC + adaptive Bayesian inference)\n");
			//t34_encode(buf, iw, ih, &cdata, 1);

			printf("\n");
		}
#endif

		//T29
#if 0
		{
			ArrayHandle cdata=0;
			double elapsed;

			elapsed=time_ms();
			cycles=__rdtsc();
			t29_encode(buf, iw, ih, &cdata, 0);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
			
			array_free(&cdata);
		}
#endif

		//test26: T16 with range coder
#if 0
		{
			int use_ans=0;
			ArrayHandle cdata=0;
			double elapsed;
			T26Params params[]=
			{
				{ 8, 26, 26, 26, 0xD3, 52},
				{23, 37, 37, 37, 0xD3, 52},
				{ 8, 26, 26, 26, 0xD3, 52},
			};
			printf("\nT26 (%s)\n", use_ans?"ANS":"AC");
			//printf("\nT26\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t26_encode(buf, iw, ih, params, use_ans, &cdata, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);

			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t26_decode(cdata->data, cdata->count, iw, ih, params, use_ans, b2, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Dec %11lf CPB  ", (double)cycles/usize);
			timedelta2str(0, 0, elapsed);
			printf("\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T26", 0);
			memset(b2, 0, len);
			printf("\n");
		}
#endif

		//T25: T16 optimizer
#if 0
		{
			int use_ans=0;
			ArrayHandle cdata=0;
			//int blockw[]={96, 96, 96}, blockh[]={96, 96, 96};
			int blockw[]={128, 128, 128}, blockh[]={128, 128, 128};
			double elapsed;
			printf("\nT25 (%s)\n", use_ans?"ANS":"AC");
			elapsed=time_ms();
			cycles=__rdtsc();
			t25_encode(buf, iw, ih, blockw, blockh, use_ans, &cdata, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			timedelta2str(0, 0, elapsed);
			
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");
		
			elapsed=time_ms();
			cycles=__rdtsc();
			t25_decode(cdata->data, cdata->count, iw, ih, blockw, blockh, use_ans, b2, 1);
			cycles=__rdtsc()-cycles;
			elapsed=time_ms()-elapsed;
			printf("Dec %11lf CPB  ", (double)cycles/usize);
			timedelta2str(0, 0, elapsed);
			printf("\n");

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, nch0, 4, "T25", 0);
			memset(b2, 0, len);
			printf("\n");
		}
#endif
		
		//test16 - THE BEST
#if 0
		{
			printf("\nT16\n");
#if 1
			memcpy(b2, buf, len);
			addbuf(b2, iw, ih, 3, 4, 128);
			colortransform_ycocb_fwd((char*)b2, iw, ih);
			pred_opt_opt_v3((char*)b2, iw, ih, 1);
			memset(b2, 0, len);
			pred_opt_printparam();
#endif

			ArrayHandle cdata=0;
			int alpha=0xD3E7,
				blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
				blockh[]={ 1,  1,  1},
				margin[]={26, 37, 26};

			cycles=__rdtsc();
			test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, 0);
			cycles=__rdtsc()-cycles;
			printf("Enc %lf CPB  CR %lf  csize %lld", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
			sum_testsize+=cdata->count;
			if((ptrdiff_t)cdata->count<formatsize)
				printf(" !!!");
			printf("\n");

			cycles=__rdtsc();
			test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
			cycles=__rdtsc()-cycles;
			printf("Dec %lf CPB\n", (double)cycles/usize);

			array_free(&cdata);
			compare_bufs_uint8(b2, buf, iw, ih, 3, 4, "T16", 0);
			memset(b2, 0, len);

			printf("\n");
		}
#endif

		//test16 estimate
#if 0
		apply_transforms_fwd(buf, iw, ih);
		double csize=test16_estimate_csize(buf, iw, ih, 32, 0.6f, 0);
		sum_testsize+=(long long)ceil(csize);
		printf("\tCR2 %f", usize/csize);
		if(csize<formatsize)
			printf(" !!!");
		printf("\n");
#endif

		//printf("\n");
		free(buf);
#ifndef BATCHTEST_NO_B2
		free(b2);
#endif
	}
	printf("Batch elapsed ");
	timedelta2str(0, 0, time_ms()-t_start);
	printf("\n");
#if 0
	e10_print(hist);
	free(hist);
#else
	ptrdiff_t totalusize=sum_uPNGsize+sum_uJPEGsize;
	if(totalusize)
	{
		printf("\nOn average:\n");
		printf("BMP     csize %9lld\n", totalusize);
		if(sum_cPNGsize)
			printf("PNG     csize %9lld  CR %lf  (%lld images)\n", sum_cPNGsize, (double)sum_uPNGsize/sum_cPNGsize, count_PNG);
		if(sum_cJPEGsize)
			printf("JPEG    csize %9lld  CR %lf  (%lld images)\n", sum_cJPEGsize, (double)sum_uJPEGsize/sum_cJPEGsize, count_JPEG);
		printf("test    csize %9lld  CR %lf\n", sum_testsize, (double)totalusize/sum_testsize);
		//printf("test3s  CR %lf\n", (double)totalusize/sum_test3size[0]);
		//printf("test3sd CR %lf\n", (double)totalusize/sum_test3size[1]);
	}
	else
		printf("\nNo valid images found\n");
#endif
	acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S");
	printf("Finish %s\n", g_buf);

	array_free(&filenames);

	printf("\nDone.\n");
	pause();
}
int main(int argc, char **argv)
{
	//const int LOL_1=(-1)/2;//0

	//int width=10,
	//	n=pyramid_getsize(width);
	//for(int k=0;k<n;++k)
	//	printf("%c", pyramid_getchar(width, k));
	//DCTtest();
	//DCTtest2();
	//test4();
	//test5();
	//test6();
	//test7();
	//test8();
	//test9();
	//test_swar();
	//system("cd");
	//init_vk();//
	//test1();

	printf("Entropy2\n");
#if 1
	long long cycles;
	int iw=0, ih=0, nch0=3,
		nch=4;
	size_t resolution=0, len=0;
	unsigned char *buf, *b2;
	const char *fn=0;
#ifdef _DEBUG
	//fn="C:/Projects/datasets/CLIC11-crop4-2.PNG";
	//fn="C:/Projects/datasets/CLIC11-small4.PNG";
	//fn="C:/Projects/datasets/dataset-CLIC30/11.png";
	//fn="C:/Projects/datasets/dataset-kodak";
	fn="C:/Projects/datasets/dataset-kodak/kodim13.png";

	//fn="D:/ML/dataset-CLIC30";
	//fn="D:/ML/dataset-kodak";
	//fn="D:/ML/dataset-CLIC30/16.png";//hardest noiseless CLIC30 image
	//fn="D:/ML/dataset-CLIC30/17.png";
	//fn="D:/ML/dataset-kodak/kodim13.png";
	//fn="D:/ML/dataset-kodak-small/13.PNG";
#endif
	if(fn||argc==2)
	{
		if(!fn)
			fn=argv[1];
		ptrdiff_t formatsize=get_filesize(fn);
		if(formatsize==-1)
		{
			LOG_ERROR("Cannot open \"%s\"", fn);
			return 0;
		}
		if(!formatsize)//path
		{
			batch_test(fn);
			return 0;
		}
		printf("Opening \"%s\"\n", fn);
		cycles=__rdtsc();
		buf=image_load(fn, &iw, &ih);
		cycles=__rdtsc()-cycles;
		if(!buf)
		{
			LOG_ERROR("Couldn't open \"%s\"", fn);
			return 0;
		}
		resolution=(size_t)iw*ih;
		len=resolution*nch;

		printf("Format Dec %lf CPB, ratio = %d * %d * %d / %lld = %lf\n", (double)cycles/(resolution*nch0), iw, ih, nch0, formatsize, (double)resolution*nch0/formatsize);
	}
	else if(argc==3)
	{
		const char *fn1=argv[1], *fn2=argv[2];
		int w2, h2;
		buf=image_load(fn1, &iw, &ih);
		b2 =image_load(fn2, &w2, &h2);
		if(!buf)
		{
			printf("Couldn't open %s\n", fn1);
			return 1;
		}
		if(!b2)
		{
			printf("Couldn't open %s\n", fn2);
			return 1;
		}
		if(iw!=w2||ih!=h2)
		{
			printf("Expected two images of SAME RESOLUTION. %dx%d != %dx%d\n", iw, ih, w2, h2);
			return 1;
		}
		ptrdiff_t formatsize=get_filesize(fn2);
		int res=iw*ih;
		long long sum[3]={0};
		for(int k=0;k<res;++k)
		{
			int dr=buf[k<<2  ]-b2[k<<2  ],
				dg=buf[k<<2|1]-b2[k<<2|1],
				db=buf[k<<2|2]-b2[k<<2|2];
			sum[0]+=dr*dr;
			sum[1]+=dg*dg;
			sum[2]+=db*db;
		}
		double rmse[]=
		{
			sqrt((double)sum[0]/res),
			sqrt((double)sum[1]/res),
			sqrt((double)sum[2]/res),
			sqrt((double)(sum[0]+sum[1]+sum[2])/(res*3)),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		double CR=res*3./formatsize;
		printf("T RMSE %lf PSNR %lf  CR %d/%d = %lf  BPP %lf\n", rmse[3], psnr[3], res*3, (int)formatsize, CR, 8/CR);
		printf("R RMSE %lf PSNR %lf\n", rmse[0], psnr[0]);
		printf("G RMSE %lf PSNR %lf\n", rmse[1], psnr[1]);
		printf("B RMSE %lf PSNR %lf\n", rmse[2], psnr[2]);
		return 0;
	}
	else
	{
		printf("Usage: e2.exe  file_or_path\n");
		pause();
		return 0;
#if 0
		iw=1920, ih=1080, nch0=3,//1080*1920*3	640*480		50		4*4*1
			nch=4;
		resolution=(size_t)iw*ih, len=resolution*nch;
		buf=(unsigned char*)malloc(len);
		if(!buf)
			return 0;
		//srand((unsigned)__rdtsc());
	
#ifdef UNIFORM
		printf("Generating test data (uniform)...\n");
		fill_uniform(buf, len);
#else
		int unibits=256;
		printf("Generating test data (%d bit binomial)...\n", unibits);
		fill_halfbinomial(buf, len, unibits);
#endif
#endif
	}

	if(nch0==3&&!buf[3])//set alpha
	{
		for(int k=3;k<len;k+=nch)
			buf[k]=0xFF;
	}

	b2=(unsigned char*)malloc(len);
	if(!b2)
		return 0;
	size_t usize=len*nch0>>2;

	printf("\n");
	
	ArrayHandle cdata=0;
	//const void *ptr, *end;
	
	int loud=0;

	//test16
#if 0
	{
		//debug_ptr=buf;//
		int besta=0, bestb=0, bestm=0, bestc=0;
		int it=0;
		for(int m=30;m<=30;++m)
		{
			for(int b=15;b<=15;++b)
			{
				for(int a=0xD3E7;a<=0xD3E7;++a, ++it)
				{
					int alpha=a,//(60<<16)/100;	0xA51F	0xD3DA	0xD3D6	0xD3EB
						bsize=b,//32 24 26 20
						margin=m;
					printf("T16 a 0x%04X b %3d m %3d ", alpha, bsize, margin);
					cycles=__rdtsc();
					test16_encode(buf, iw, ih, alpha, bsize, margin, &cdata, 1);
					cycles=__rdtsc()-cycles;
					printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

					if(!it||bestc>(int)cdata->count)
						besta=alpha, bestb=bsize, bestm=m, bestc=(int)cdata->count;

					cycles=__rdtsc();
					test16_decode(cdata->data, cdata->count, iw, ih, alpha, bsize, margin, b2);
					cycles=__rdtsc()-cycles;
					printf("Dec %11lf CPB ", (double)cycles/usize);

					array_free(&cdata);
					compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
					memset(b2, 0, len);
				}
				//printf("\n");
			}
		}
		if(it>1)
			printf("T16 best ABMS 0x%04X %2d %2d %d CR %lf\n", besta, bestb, bestm, bestc, (double)usize/bestc);
		else
			printf("\n");
	}
#endif


	//test16 codec with jxl predictor optimizer
#if 0
	{
		int alpha=0xD3E7,
			blockw[]={ 8, 23,  8},//best block for channels 0 & 2: 1x1
			blockh[]={ 1,  1,  1},
			margin[]={26, 37, 26};

#if 0
		int res=iw*ih;
		double step=0.001, CR0=0, CR, csize[3]={0};
		estimate_csize_from_transforms(buf, b2, iw, ih, csize);
		CR=res*3/(csize[0]+csize[1]+csize[2]);
		printf("%4d TRGB %lf [%lf %lf %lf]\n", 0, CR, res/csize[0], res/csize[1], res/csize[2]);

		for(int k=0;k<256;++k)
		{
			int idx=k%33;
			if(!(k+1)%33)
				step*=0.9;
			do
			{
				CR0=CR;
				jxlpred_params[idx]+=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);

			do
			{
				CR0=CR;
				jxlpred_params[idx]-=step;
				estimate_csize_from_transforms(buf, b2, iw, ih, csize);
				CR=res*3/(csize[0]+csize[1]+csize[2]);
				printf("%4d TRGB %lf [%lf %lf %lf]\n", k+1, CR, res/csize[0], res/csize[1], res/csize[2]);
			}
			while(CR>CR0);
		}
#endif

		//for(int k=0;k<3;++k)
		//	printf("%g\t%g\t%g\t%g\n", jxlpred_params[k<<2], jxlpred_params[k<<2|1], jxlpred_params[k<<2|2], jxlpred_params[k<<2|3]);
		//for(int k=12;k<33;k+=7)
		//	printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", jxlpred_params[k], jxlpred_params[k+1], jxlpred_params[k+2], jxlpred_params[k+3], jxlpred_params[k+4], jxlpred_params[k+5], jxlpred_params[k+6]);
		//for(int k=0;k<33;++k)
		//	printf("%3d  %lf\n", k, jxlpred_params[k]);
		//printf("\n");
		
		printf("T16\n");
		int bestcsizes[3]={0}, bestw[3]={0}, besth[3]={0}, bestm[3]={0};
		int it=0;
		//for(int bw=19;bw<25;++bw)//
		{
			//for(int m=31;m<48;++m)
			{
				int csizes[3];
				//blockw[1]=bw, blockh[1]=1, margin[1]=m;
				//blockw[1]=1+k, blockh[1]=1;
				//blockw[1]=4+k%10, blockh[1]=1+k/10;
			
				//blockw[0]=blockw[2]=blockw[1];
				//blockh[0]=blockh[2]=blockh[1];

				cycles=__rdtsc();
				test16_encode(buf, iw, ih, alpha, blockw, blockh, margin, &cdata, 1, csizes);
				cycles=__rdtsc()-cycles;
				printf("Enc %11lf CPB  CR %9lf  csize %lld ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);

				for(int kc=0;kc<3;++kc)
				{
					if(!it||bestcsizes[kc]>csizes[kc])
						bestcsizes[kc]=csizes[kc], bestw[kc]=blockw[kc], besth[kc]=blockh[kc], bestm[kc]=margin[kc];
				}

				cycles=__rdtsc();
				test16_decode(cdata->data, cdata->count, iw, ih, alpha, blockw, blockh, margin, b2);
				cycles=__rdtsc()-cycles;
				printf("Dec %11lf CPB ", (double)cycles/usize);

				array_free(&cdata);
				compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T16", 0);
				memset(b2, 0, len);
				printf("\n");
				++it;
			}
		}
		int res=iw*ih;
		printf("R %7d %lf %2dx%d  M %d\n", bestcsizes[0], (double)res/bestcsizes[0], bestw[0], besth[0], bestm[0]);
		printf("G %7d %lf %2dx%d  M %d\n", bestcsizes[1], (double)res/bestcsizes[1], bestw[1], besth[1], bestm[1]);
		printf("B %7d %lf %2dx%d  M %d\n", bestcsizes[2], (double)res/bestcsizes[2], bestw[2], besth[2], bestm[2]);
		printf("\n");
	}
#endif

	//test25: T16 optimizer
#if 0
	{
		int use_ans=0;
		printf("T25 (%s)\n", use_ans?"ANS":"AC");
		double elapsed;
		int blockw[]={128, 128, 128}, blockh[]={128, 128, 128};
		//int lbsizes[]=
		//{
		//	32, 32,
		//	32, 32,
		//	32, 32,
		//};
		//int sbsizes[]=//unused
		//{
		//	16, 16,
		//	16, 16,
		//	16, 16,
		//};
		elapsed=time_ms();
		cycles=__rdtsc();
		t25_encode(buf, iw, ih, blockw, blockh, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t25_decode(cdata->data, cdata->count, iw, ih, blockw, blockh, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T25", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//test26: T16 with range coder
#if 0
	{
		int use_ans=0;
		double elapsed;
		T26Params params[]=
		{
			{ 8, 26, 26, 26, 0xD3, 52},
			{23, 37, 37, 37, 0xD3, 52},
			{ 8, 26, 26, 26, 0xD3, 52},
		};
		printf("T26 (%s)\n", use_ans?"ANS":"AC");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t26_encode(buf, iw, ih, params, use_ans, &cdata, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Enc %11lf CPB  CR %9lf  csize %lld  ", (double)cycles/usize, (double)usize/cdata->count, cdata->count);
		timedelta2str(0, 0, elapsed);
		printf("\n");
		
		elapsed=time_ms();
		cycles=__rdtsc();
		t26_decode(cdata->data, cdata->count, iw, ih, params, use_ans, b2, 1);
		cycles=__rdtsc()-cycles;
		elapsed=time_ms()-elapsed;
		printf("Dec %11lf CPB  ", (double)cycles/usize);
		timedelta2str(0, 0, elapsed);
		printf("\n");

		array_free(&cdata);
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T26", 0);
		memset(b2, 0, len);
		printf("\n");
	}
#endif

	//T27: ABAC
#if 1
	//printf("T29 (ABAC + Bayesian inference)\n");
	//t29_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T30 (ABAC + Bayesian inference + predictor)\n");	//X
	//t30_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T31 (ABAC + adaptive Bayesian inference)\n");
	//t31_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T33 (ABAC + adaptive Bayesian inference with circular buffer)\n");	//X
	//t33_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);

	//printf("T34 (ABAC + adaptive Bayesian inference)\n");
	//t34_encode(buf, iw, ih, &cdata, 1);
	//array_free(&cdata);
	
	//old record
#if 0
	printf("T35 Entropy coding with context tree\n");
	//printf("T35 Combines spatial transform with entropy coding\n");
	t35_encode(buf, iw, ih, &cdata, 1);
	t35_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T35", 0);
	memset(b2, 0, len);
	printf("\n");
#endif
	
	//printf("T36 stretch & squish\n");
	//t36_encode(buf, iw, ih, &cdata, 1);
	//t36_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T36", 0);
	//memset(b2, 0, len);
	//printf("\n");
	
	//printf("T37 Fixed array as binary tree predictor\n");
	//t37_encode(buf, iw, ih, &cdata, 1);
	//t37_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T37", 0);
	//memset(b2, 0, len);
	//printf("\n");
	
	//printf("T38 Single simple bit predictor\n");
	//t38_encode(buf, iw, ih, &cdata, 1);
	//t38_decode(cdata->data, cdata->count, iw, ih, b2, 1);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T38", 0);
	//memset(b2, 0, len);
	//printf("\n");
	
	//record
#if 0
	printf("T39 Multiple estimators for all maps\n");
	t39_encode(buf, iw, ih, &cdata, 2);
	t39_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T39", 0);
	memset(b2, 0, len);
	printf("\n");
#endif
	
	//t40_encode(buf, iw, ih, &cdata, 2);	//X
	//t40_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T40", 0);
	//memset(b2, 0, len);
	//printf("\n");
	
	//t41_encode(buf, iw, ih, &cdata, 2);	//X
	//t41_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	//array_free(&cdata);
	//compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T41", 0);
	//memset(b2, 0, len);
	//printf("\n");
	
	t42_encode(buf, iw, ih, &cdata, 2);
	t42_decode(cdata->data, cdata->count, iw, ih, b2, 2);
	array_free(&cdata);
	compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "T42", 0);
	memset(b2, 0, len);
	printf("\n");
#endif

	//predict image
	apply_transforms_fwd(buf, iw, ih);
	//lodepng_encode_file("kodim21-YCoCgT-unplane.PNG", buf, iw, ih, LCT_RGBA, 8);//
	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//

	//colortransform_ycocb_fwd(buf, iw, ih);
	//save_channel(buf+1, iw, ih, 4, 0, "kodim13-YCoCb-jxl-luma.PNG");
#if 0
	printf("Predict image...\n");
	{
		ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 0);
		char *temp=(char*)malloc(MAXVAR(iw, ih));
		
		addbuf(buf, iw, ih, nch0, nch, 128);//unsigned char -> signed char
		
		//colortransform_ycocg_fwd((char*)buf, iw, ih);
		//colortransform_xgz_fwd((char*)buf, iw, ih);
		//colortransform_xyz_fwd((char*)buf, iw, ih);

		//char *b3=(char*)malloc(iw), *b4=(char*)malloc(iw);
		//if(!b3||!b4)
		//	return 0;
		//memcpy(b3, buf, iw);
		//dwt1d_squeeze_fwd(b3, iw, 1, b4);
		//dwt1d_squeeze_inv(b3, iw, 1, b4);
		//compare_bufs_uint8((unsigned char*)b3, buf, iw, 1, 1, 1, "squeeze row", 0);
		//free(b3);
		//free(b4);

#if 1
		memcpy(b2, buf, len);

		colortransform_ycocb_fwd((char*)b2, iw, ih);
		float jxlparams[33]=
		{
			 0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
			 0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
			 0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
		};
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 1);
		pred_jxl_apply((char*)b2, iw, ih, jxlparams, 0);

		//colortransform_xyz_fwd(b2, iw, ih);
		//colortransform_xyz_inv(b2, iw, ih);
		//for(int kc=0;kc<3;++kc)
		//{
		//	dwt2d_squeeze_fwd((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//	dwt2d_squeeze_inv((char*)b2+kc, (DWTSize*)sizes->data, 0, 2, 4, temp);
		//}
		compare_bufs_uint8(b2, buf, iw, ih, nch0, nch, "transform", 0);
		printf("\n");
#endif
		
		//for(int kc=0;kc<3;++kc)
		//	//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	dwt2d_squeeze_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, 2, 4, (char*)temp);
		//	//dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
		//	//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

		addbuf(buf, iw, ih, nch0, nch, 128);

		//save_DWT_int8("kodim21-squeeze-stage", buf, (DWTSize*)sizes->data, 2, 4);//
		//lodepng_encode_file("kodim21-cubic.PNG", buf, iw, ih, LCT_RGBA, 8);//

		array_free(&sizes);
		free(temp);
	}
	//squeeze_8bit_lossy(buf, iw, ih, nch0, nch);
//	image_pred(buf, iw, ih, nch0, nch);

	//lodepng_encode_file("kodim21-XGZ-diff2d.PNG", buf, iw, ih, LCT_RGBA, 8);//
#endif

	for(int kc=0;kc<nch0;++kc)
		calc_histogram(buf+kc, len, nch, hist+((size_t)kc<<8));
	//print_histogram(hist, 1);

	double entropy[6]={0};
	for(int kc=0;kc<nch0;++kc)
	{
		int freq;
		double p;
		for(int k=0;k<256;++k)
		{
			freq=hist[kc<<8|k];
			if(freq)
			{
				p=(double)freq/(len>>2);
				p*=0x10000-255;
				++p;
				p/=0x10000;
				entropy[kc]+=-p*log2(p);
			}

			//if(!kc)//
			//	printf("%3d %6d %lf\n", k, freq, p);//
		}
		printf("ch %d E = %lf / 8, optimal ratio = %lf\n", kc, entropy[kc], 8/entropy[kc]);
		entropy[4]+=entropy[kc];
	}
	entropy[4]/=nch0;
	printf("Av. E = %lf / 8, optimal ratio = %lf <- true limit\n", entropy[4], 8/entropy[4]);


	free(buf);
	free(b2);
#endif
	printf("Done.\n");
	pause();
	return 0;
}