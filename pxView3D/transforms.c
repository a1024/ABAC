#include"pxview3d.h"
#include<stdlib.h>
#include<math.h>
static const char file[]=__FILE__;

void addhalf(unsigned char *buf, int iw, int ih, int nch, int bytestride)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=128;
	}
}

void colortransform_xgz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xgz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xyz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		//r-=g;
		//g+=r>>1;
		//b-=g;
		//g+=b>>1;

		r-=b;
		g-=b;
		b+=(r+g)>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;

		//buf[k  ]=r-g+128;//XGZ
		//buf[k|1]=g;
		//buf[k|2]=b-g+128;

		//buf[k  ]=r;//RYZ
		//buf[k|1]=g-r+128;
		//buf[k|2]=b-r+128;

		//buf[k  ]=r-b+128;//XYBdash
		//buf[k|1]=g-b+128;
		//buf[k|2]=b;
	}
}
void colortransform_xyz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		//g-=b>>1;
		//b+=g;
		//g-=r>>1;
		//r+=g;

		b-=(r+g)>>1;
		g+=b;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=b;		//co = r-b			diff(r, b)
		b+=r>>1;	//(r+b)/2
		g-=b;		//cg = g-(r+b)/2	diff(g, av(r, b))
		b+=g>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		b-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocgt_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;		//diff(r, g)
		g+=r>>1;
		b-=g;		//diff(b, av(r, g))
		g+=b>>1;	//av(b, av(r, g))

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocgt_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}

void colortransform_exp_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char
				r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2], r2=0, g2=0, b2=0, r3=0, g3=0, b3=0;
			if(kx+1<iw)//right
				r2=buf[(idx+1)<<2], g2=buf[(idx+1)<<2|1], b2=buf[(idx+1)<<2|2];
			if(ky+1<ih)//below
				r3=buf[(idx+iw)<<2], g3=buf[(idx+iw)<<2|1], b3=buf[(idx+iw)<<2|2];

			if(abs(r2-g2)<abs(r2-b2)&&abs(r3-g3)<abs(r3-b3))
			{
				r-=g;
				g+=r>>1;
			}
			else if(abs(r2-g2)>abs(r2-b2)&&abs(r3-g3)>abs(r3-b3))
			{
				r-=b;
				b+=r>>1;
			}

			//if(abs(r2-g2-(r3-g3))<abs(r2-r3))
			//	r-=g;
			//else
			//	r-=(r2+r3)>>1;

			//r-=(g+g+r2+r3)>>2;

			//if(abs(r2-r3)<abs(g-b))
			//	r-=(r2+r3)>>1;
			//else
			//	r-=g;

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}
void colortransform_exp_inv(char *buf, int iw, int ih)
{
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx)
		{
			int idx=iw*ky+kx;
			char
				r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2], r2=0, g2=0, b2=0, r3=0, g3=0, b3=0;
			if(kx+1<iw)//right
				r2=buf[(idx+1)<<2], g2=buf[(idx+1)<<2|1], b2=buf[(idx+1)<<2|2];
			if(ky+1<ih)//below
				r3=buf[(idx+iw)<<2], g3=buf[(idx+iw)<<2|1], b3=buf[(idx+iw)<<2|2];
			
			if(abs(r2-g2)<abs(r2-b2)&&abs(r3-g3)<abs(r3-b3))
			{
				r+=g;
				g-=r>>1;
			}
			else if(abs(r2-g2)>abs(r2-b2)&&abs(r3-g3)>abs(r3-b3))
			{
				r+=b;
				b-=r>>1;
			}

			//if(abs(r2-g2-(r3-g3))<abs(r2-r3))
			//	r+=g;
			//else
			//	r+=(r2+r3)>>1;
			
			//r+=(g+g+r2+r3)>>2;

			//if(abs(r2-r3)<abs(g-b))
			//	r-=(r2+r3)>>1;
			//else
			//	r-=g;

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}

float lrct[]=
{
	-0.235968455672264100f, -0.005472748540341854f,  0.546232640743255600f,
	 0.042025841772556305f, -0.281145840883255000f, -0.041713438928127290f,
	-0.270333796739578250f,  0.382541775703430200f,  0.122593350708484650f,
	-0.654174447059631300f, -0.111721098423004150f,  0.049731694161891940f,
	-0.180651888251304630f, -0.665180742740631100f, -0.449275195598602300f,
	-0.634068489074707000f, -0.438444226980209350f, -0.079102315008640290f,
	-0.492376208305358900f, -0.011669249273836613f,  0.612039744853973400f,
	-0.140227079391479500f,  0.255208671092987060f, -0.405863404273986800f,
	-0.133800372481346130f, -0.049969632178545000f, -0.006972887087613344f,
	-0.498686790466308600f, -0.417643666267395000f,  0.129942402243614200f,
	-0.584400236606597900f, -0.273861020803451540f,  0.104215130209922790f,
	 0.297003507614135740f, -0.688771069049835200f,  0.330031186342239400f,
};
char lift(char v1, char v2, float *coeff)
{
	return (char)(127*(v1*coeff[0]/127+v2*coeff[1]/127+coeff[2]));
}
void colortransform_learned_fwd(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		r-=lift(g, b, lrct+3* 0);
		g-=lift(r, b, lrct+3* 1);
		b-=lift(r, g, lrct+3* 2);
		r+=lift(g, b, lrct+3* 3);
		g+=lift(r, b, lrct+3* 4);
		b+=lift(r, g, lrct+3* 5);
		r-=lift(g, b, lrct+3* 6);
		g-=lift(r, b, lrct+3* 7);
		b-=lift(r, g, lrct+3* 8);
		r+=lift(g, b, lrct+3* 9);
		g+=lift(r, b, lrct+3*10);
		b+=lift(r, g, lrct+3*11);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
void colortransform_learned_inv(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		b-=lift(r, g, lrct+3*11);
		g-=lift(r, b, lrct+3*10);
		r-=lift(g, b, lrct+3* 9);
		b+=lift(r, g, lrct+3* 8);
		g+=lift(r, b, lrct+3* 7);
		r+=lift(g, b, lrct+3* 6);
		b-=lift(r, g, lrct+3* 5);
		g-=lift(r, b, lrct+3* 4);
		r-=lift(g, b, lrct+3* 3);
		b+=lift(r, g, lrct+3* 2);
		g+=lift(r, b, lrct+3* 1);
		r+=lift(g, b, lrct+3* 0);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
#if 0
#if 1
float lrt_biases[]=
{
	-0.007593977730721235f,
	-0.20463228225708008f,
	0.048144787549972534f,
	-0.052841395139694214f,
	-0.229848250746727f,
	-0.10854046046733856f,
	-0.27017149329185486f,
	-0.28025728464126587f,
	0.14076033234596252f,
	0.11415114998817444f,
	-0.24117304384708405f,
	-0.22909477353096008f,
};
float lrt_c01[]=
{
	 0.21996481716632843f,
	 0.3447851538658142f ,
	 0.4309721887111664f ,
	-0.3914514482021332f ,
	 0.1072743684053421f ,
	 0.09292115271091461f,
	 0.09527700394392014f,
	-0.7720031142234802f ,
	-0.03891817852854729f,
	 0.18422403931617737f,
	 0.6045548319816589f ,
};
float lrt_c02[]=
{
	0.16338559985160828f,
	-0.21537575125694275f,
	-0.01003316044807434f,
	0.0806010365486145f,
	0.1516939401626587f,
	0.22873657941818237f,
	-0.13795271515846252f,
	0.13615399599075317f,
	-0.07778285443782806f,
	0.291059672832489f,
	0.22079020738601685f,
};
float lrt_c03[]=
{
	0.010352730751037598f,
	0.1270458996295929f,
	-0.18215587735176086f,
	-0.164699524641037f,
	0.23722952604293823f,
	0.27438366413116455f,
	-0.05265103280544281f,
	0.18550443649291992f,
	-0.21034197509288788f,
	-0.13812671601772308f,
	-0.015580296516418457f,
};
float lrt_c04[]=
{
	-0.06830897927284241f,
	-0.0477299690246582f,
	-0.2479163110256195f,
	-0.04657846689224243f,
	-0.0802721232175827f,
	-0.004736065864562988f,
	-0.07743765413761139f,
	-0.26771998405456543f,
	0.013397663831710815f,
	-0.23979492485523224f,
	0.12400856614112854f,
};
float lrt_c05[]=
{
	-0.04351922869682312f,
	-0.11087171733379364f,
	-0.11081892251968384f,
	0.301426351070404050f,
	0.090425968170166020f,
	0.171547293663024900f,
	0.166505038738250730f,
	0.086036622524261470f,
	-0.25268214941024780f,
	0.168854206800460820f,
	0.221817135810852050f,
};
float lrt_c06[]=
{
	0.216941952705383300f,
	0.145165383815765380f,
	-0.30066022276878357f,
	0.062808305025100710f,
	-0.11871317028999329f,
	-0.03821510076522827f,
	0.057682424783706665f,
	-0.15345700085163116f,
	-0.10404434800148010f,
	0.159864872694015500f,
	0.030692487955093384f,
};
float lrt_c07[]=
{
	-0.11159592866897583f,
	0.267521142959594700f,
	-0.04237860441207886f,
	-0.12383817136287689f,
	0.221899330615997310f,
	0.031740218400955200f,
	-0.11875738203525543f,
	0.109959483146667480f,
	-0.17683275043964386f,
	-0.12659740447998047f,
	-0.21873277425765990f,
};
float lrt_c08[]=
{
	-0.18569713830947876f,
	0.224896728992462160f,
	-0.28754058480262756f,
	0.111558407545089720f,
	-0.18325410783290863f,
	-0.27818745374679565f,
	0.241845905780792240f,
	0.179279416799545300f,
	-0.05782620608806610f,
	-0.10470718145370483f,
	0.128715872764587400f,
};
float lrt_c09[]=
{
	0.0736332833766937300f,
	-0.278686881065368650f,
	-0.062273472547531130f,
	0.2923494577407837000f,
	0.1254649162292480500f,
	0.0480349659919738800f,
	0.2714382410049438500f,
	-0.016680717468261720f,
	0.0569746196269989000f,
	-0.021338820457458496f,
	0.0394396781921386700f,
};
float lrt_c10[]=
{
	0.068611204624176030f,
	-0.12615610659122467f,
	-0.07715198397636414f,
	0.252162516117095950f,
	-0.02322089672088623f,
	0.062368571758270264f,
	0.138979107141494750f,
	0.039735138416290280f,
	0.171300053596496580f,
	-0.13656683266162872f,
	-0.07565732300281525f,
};
float lrt_c11[]=
{
	0.1174532473087310800f,
	0.0668951570987701400f,
	-0.071375116705894470f,
	-0.240586921572685240f,
	0.1582727432250976600f,
	0.1294519901275634800f,
	0.2438962459564209000f,
	-0.014519184827804565f,
	0.1907848119735717800f,
	-0.215079009532928470f,
	0.0485148429870605500f,
};
float lrt_c12[]=
{
	0.086784452199935910f,
	0.276734590530395500f,
	-0.24387010931968690f,
	0.102736681699752810f,
	-0.28351089358329773f,
	0.072889894247055050f,
	-0.18153974413871765f,
	0.022794693708419800f,
	-0.27390283346176150f,
	-0.27018001675605774f,
	0.102353841066360470f,
};
#endif
void learnedtransform_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			//if(!kx&&!ky)
			//	kx=0;

			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x01-=(char)(127*(lrt_c01[0]*x02/127+lrt_c01[1]*x03/127+lrt_c01[2]*x04/127+lrt_c01[3]*x05/127+lrt_c01[4]*x06/127+lrt_c01[5]*x07/127+lrt_c01[6]*x08/127+lrt_c01[7]*x09/127+lrt_c01[8]*x10/127+lrt_c01[9]*x11/127+lrt_c01[10]*x12/127+lrt_biases[ 0]));
			x02-=(char)(127*(lrt_c02[0]*x01/127+lrt_c02[1]*x03/127+lrt_c02[2]*x04/127+lrt_c02[3]*x05/127+lrt_c02[4]*x06/127+lrt_c02[5]*x07/127+lrt_c02[6]*x08/127+lrt_c02[7]*x09/127+lrt_c02[8]*x10/127+lrt_c02[9]*x11/127+lrt_c02[10]*x12/127+lrt_biases[ 1]));
			x03-=(char)(127*(lrt_c03[0]*x01/127+lrt_c03[1]*x02/127+lrt_c03[2]*x04/127+lrt_c03[3]*x05/127+lrt_c03[4]*x06/127+lrt_c03[5]*x07/127+lrt_c03[6]*x08/127+lrt_c03[7]*x09/127+lrt_c03[8]*x10/127+lrt_c03[9]*x11/127+lrt_c03[10]*x12/127+lrt_biases[ 2]));
			x04-=(char)(127*(lrt_c04[0]*x01/127+lrt_c04[1]*x02/127+lrt_c04[2]*x03/127+lrt_c04[3]*x05/127+lrt_c04[4]*x06/127+lrt_c04[5]*x07/127+lrt_c04[6]*x08/127+lrt_c04[7]*x09/127+lrt_c04[8]*x10/127+lrt_c04[9]*x11/127+lrt_c04[10]*x12/127+lrt_biases[ 3]));
			x05-=(char)(127*(lrt_c05[0]*x01/127+lrt_c05[1]*x02/127+lrt_c05[2]*x03/127+lrt_c05[3]*x04/127+lrt_c05[4]*x06/127+lrt_c05[5]*x07/127+lrt_c05[6]*x08/127+lrt_c05[7]*x09/127+lrt_c05[8]*x10/127+lrt_c05[9]*x11/127+lrt_c05[10]*x12/127+lrt_biases[ 4]));
			x06-=(char)(127*(lrt_c06[0]*x01/127+lrt_c06[1]*x02/127+lrt_c06[2]*x03/127+lrt_c06[3]*x04/127+lrt_c06[4]*x05/127+lrt_c06[5]*x07/127+lrt_c06[6]*x08/127+lrt_c06[7]*x09/127+lrt_c06[8]*x10/127+lrt_c06[9]*x11/127+lrt_c06[10]*x12/127+lrt_biases[ 5]));
			x07-=(char)(127*(lrt_c07[0]*x01/127+lrt_c07[1]*x02/127+lrt_c07[2]*x03/127+lrt_c07[3]*x04/127+lrt_c07[4]*x05/127+lrt_c07[5]*x06/127+lrt_c07[6]*x08/127+lrt_c07[7]*x09/127+lrt_c07[8]*x10/127+lrt_c07[9]*x11/127+lrt_c07[10]*x12/127+lrt_biases[ 6]));
			x08-=(char)(127*(lrt_c08[0]*x01/127+lrt_c08[1]*x02/127+lrt_c08[2]*x03/127+lrt_c08[3]*x04/127+lrt_c08[4]*x05/127+lrt_c08[5]*x06/127+lrt_c08[6]*x07/127+lrt_c08[7]*x09/127+lrt_c08[8]*x10/127+lrt_c08[9]*x11/127+lrt_c08[10]*x12/127+lrt_biases[ 7]));
			x09-=(char)(127*(lrt_c09[0]*x01/127+lrt_c09[1]*x02/127+lrt_c09[2]*x03/127+lrt_c09[3]*x04/127+lrt_c09[4]*x05/127+lrt_c09[5]*x06/127+lrt_c09[6]*x07/127+lrt_c09[7]*x08/127+lrt_c09[8]*x10/127+lrt_c09[9]*x11/127+lrt_c09[10]*x12/127+lrt_biases[ 8]));
			x10-=(char)(127*(lrt_c10[0]*x01/127+lrt_c10[1]*x02/127+lrt_c10[2]*x03/127+lrt_c10[3]*x04/127+lrt_c10[4]*x05/127+lrt_c10[5]*x06/127+lrt_c10[6]*x07/127+lrt_c10[7]*x08/127+lrt_c10[8]*x09/127+lrt_c10[9]*x11/127+lrt_c10[10]*x12/127+lrt_biases[ 9]));
			x11-=(char)(127*(lrt_c11[0]*x01/127+lrt_c11[1]*x02/127+lrt_c11[2]*x03/127+lrt_c11[3]*x04/127+lrt_c11[4]*x05/127+lrt_c11[5]*x06/127+lrt_c11[6]*x07/127+lrt_c11[7]*x08/127+lrt_c11[8]*x09/127+lrt_c11[9]*x10/127+lrt_c11[10]*x12/127+lrt_biases[10]));
			x12-=(char)(127*(lrt_c12[0]*x01/127+lrt_c12[1]*x02/127+lrt_c12[2]*x03/127+lrt_c12[3]*x04/127+lrt_c12[4]*x05/127+lrt_c12[5]*x06/127+lrt_c12[6]*x07/127+lrt_c12[7]*x08/127+lrt_c12[8]*x09/127+lrt_c12[9]*x10/127+lrt_c12[10]*x11/127+lrt_biases[11]));
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04;
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08;
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_fwd(buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);
}
void learnedtransform_inv(char *buf, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_inv(buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);

	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			//if(!kx&&!ky)
			//	kx=0;

			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x12+=(char)(127*(lrt_c12[0]*x01/127+lrt_c12[1]*x02/127+lrt_c12[2]*x03/127+lrt_c12[3]*x04/127+lrt_c12[4]*x05/127+lrt_c12[5]*x06/127+lrt_c12[6]*x07/127+lrt_c12[7]*x08/127+lrt_c12[8]*x09/127+lrt_c12[9]*x10/127+lrt_c12[10]*x11/127+lrt_biases[11]));
			x11+=(char)(127*(lrt_c11[0]*x01/127+lrt_c11[1]*x02/127+lrt_c11[2]*x03/127+lrt_c11[3]*x04/127+lrt_c11[4]*x05/127+lrt_c11[5]*x06/127+lrt_c11[6]*x07/127+lrt_c11[7]*x08/127+lrt_c11[8]*x09/127+lrt_c11[9]*x10/127+lrt_c11[10]*x12/127+lrt_biases[10]));
			x10+=(char)(127*(lrt_c10[0]*x01/127+lrt_c10[1]*x02/127+lrt_c10[2]*x03/127+lrt_c10[3]*x04/127+lrt_c10[4]*x05/127+lrt_c10[5]*x06/127+lrt_c10[6]*x07/127+lrt_c10[7]*x08/127+lrt_c10[8]*x09/127+lrt_c10[9]*x11/127+lrt_c10[10]*x12/127+lrt_biases[ 9]));
			x09+=(char)(127*(lrt_c09[0]*x01/127+lrt_c09[1]*x02/127+lrt_c09[2]*x03/127+lrt_c09[3]*x04/127+lrt_c09[4]*x05/127+lrt_c09[5]*x06/127+lrt_c09[6]*x07/127+lrt_c09[7]*x08/127+lrt_c09[8]*x10/127+lrt_c09[9]*x11/127+lrt_c09[10]*x12/127+lrt_biases[ 8]));
			x08+=(char)(127*(lrt_c08[0]*x01/127+lrt_c08[1]*x02/127+lrt_c08[2]*x03/127+lrt_c08[3]*x04/127+lrt_c08[4]*x05/127+lrt_c08[5]*x06/127+lrt_c08[6]*x07/127+lrt_c08[7]*x09/127+lrt_c08[8]*x10/127+lrt_c08[9]*x11/127+lrt_c08[10]*x12/127+lrt_biases[ 7]));
			x07+=(char)(127*(lrt_c07[0]*x01/127+lrt_c07[1]*x02/127+lrt_c07[2]*x03/127+lrt_c07[3]*x04/127+lrt_c07[4]*x05/127+lrt_c07[5]*x06/127+lrt_c07[6]*x08/127+lrt_c07[7]*x09/127+lrt_c07[8]*x10/127+lrt_c07[9]*x11/127+lrt_c07[10]*x12/127+lrt_biases[ 6]));
			x06+=(char)(127*(lrt_c06[0]*x01/127+lrt_c06[1]*x02/127+lrt_c06[2]*x03/127+lrt_c06[3]*x04/127+lrt_c06[4]*x05/127+lrt_c06[5]*x07/127+lrt_c06[6]*x08/127+lrt_c06[7]*x09/127+lrt_c06[8]*x10/127+lrt_c06[9]*x11/127+lrt_c06[10]*x12/127+lrt_biases[ 5]));
			x05+=(char)(127*(lrt_c05[0]*x01/127+lrt_c05[1]*x02/127+lrt_c05[2]*x03/127+lrt_c05[3]*x04/127+lrt_c05[4]*x06/127+lrt_c05[5]*x07/127+lrt_c05[6]*x08/127+lrt_c05[7]*x09/127+lrt_c05[8]*x10/127+lrt_c05[9]*x11/127+lrt_c05[10]*x12/127+lrt_biases[ 4]));
			x04+=(char)(127*(lrt_c04[0]*x01/127+lrt_c04[1]*x02/127+lrt_c04[2]*x03/127+lrt_c04[3]*x05/127+lrt_c04[4]*x06/127+lrt_c04[5]*x07/127+lrt_c04[6]*x08/127+lrt_c04[7]*x09/127+lrt_c04[8]*x10/127+lrt_c04[9]*x11/127+lrt_c04[10]*x12/127+lrt_biases[ 3]));
			x03+=(char)(127*(lrt_c03[0]*x01/127+lrt_c03[1]*x02/127+lrt_c03[2]*x04/127+lrt_c03[3]*x05/127+lrt_c03[4]*x06/127+lrt_c03[5]*x07/127+lrt_c03[6]*x08/127+lrt_c03[7]*x09/127+lrt_c03[8]*x10/127+lrt_c03[9]*x11/127+lrt_c03[10]*x12/127+lrt_biases[ 2]));
			x02+=(char)(127*(lrt_c02[0]*x01/127+lrt_c02[1]*x03/127+lrt_c02[2]*x04/127+lrt_c02[3]*x05/127+lrt_c02[4]*x06/127+lrt_c02[5]*x07/127+lrt_c02[6]*x08/127+lrt_c02[7]*x09/127+lrt_c02[8]*x10/127+lrt_c02[9]*x11/127+lrt_c02[10]*x12/127+lrt_biases[ 1]));
			x01+=(char)(127*(lrt_c01[0]*x02/127+lrt_c01[1]*x03/127+lrt_c01[2]*x04/127+lrt_c01[3]*x05/127+lrt_c01[4]*x06/127+lrt_c01[5]*x07/127+lrt_c01[6]*x08/127+lrt_c01[7]*x09/127+lrt_c01[8]*x10/127+lrt_c01[9]*x11/127+lrt_c01[10]*x12/127+lrt_biases[ 0]));
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04;
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08;
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
}
#endif

const int customparam_ct_w=2, customparam_ct_h=6, customparam_st_reach=2;
/*const double customparam_ct0[]=
{
	 1,   0,
	-0.5, 0,
	 0,   1,
	 0,   0,
	 0,   1,
	 0,   0,
};
const double customparam_st0[]=
{
	0,  0, 0, 0, 0,
	0, -1, 1, 0, 0,
	0,  1,
};*/
int customparam_sel=12;
double customparam_ct[12]={0};
double customparam_st[12]={0};
void customtransforms_resetparams()
{
	memset(customparam_ct, 0, sizeof(customparam_ct));
	memset(customparam_st, 0, sizeof(customparam_st));
	//memcpy(customparam_ct, customparam_ct0, sizeof(customparam_ct));
	//memcpy(customparam_st, customparam_st0, sizeof(customparam_st));
}
void colortransform_custom_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			char r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2];
			
			r-=(char)(customparam_ct[ 0]*g+customparam_ct[ 1]*b);
			g-=(char)(customparam_ct[ 2]*r+customparam_ct[ 3]*b);
			b-=(char)(customparam_ct[ 4]*r+customparam_ct[ 5]*g);
			r+=(char)(customparam_ct[ 6]*g+customparam_ct[ 7]*b);
			g+=(char)(customparam_ct[ 8]*r+customparam_ct[ 9]*b);
			b+=(char)(customparam_ct[10]*r+customparam_ct[11]*g);

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}
void colortransform_custom_inv(char *buf, int iw, int ih)
{
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx)
		{
			int idx=iw*ky+kx;
			char r=buf[idx<<2], g=buf[idx<<2|1], b=buf[idx<<2|2];
			
			b-=(char)(customparam_ct[10]*r+customparam_ct[11]*g);
			g-=(char)(customparam_ct[ 8]*r+customparam_ct[ 9]*b);
			r-=(char)(customparam_ct[ 6]*g+customparam_ct[ 7]*b);
			b+=(char)(customparam_ct[ 4]*r+customparam_ct[ 5]*g);
			g+=(char)(customparam_ct[ 2]*r+customparam_ct[ 3]*b);
			r+=(char)(customparam_ct[ 0]*g+customparam_ct[ 1]*b);

			buf[idx<<2  ]=r;
			buf[idx<<2|1]=g;
			buf[idx<<2|2]=b;
		}
	}
}


//spatial transforms
void image_differentiate(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				//if(kx||ky)
				//	sub-=128;
				buf[idx]-=sub;
			}
		}
	}
}
void image_integrate(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				//if(kx||ky)
				//	sub-=128;
				buf[idx]+=sub;
			}
		}
	}
}

void image_unplane(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					
					leftleft=kx-2>=0?buf[idx-(bytestride<<1)]:0,
					toptop=ky-2>=0?buf[idx-(rowlen<<1)]:0,
					topright=kx<iw-1&&ky?buf[idx-rowlen+bytestride]:0,

					pred;

				char xdelta=top-topleft, ydelta=left-topleft;
				if((xdelta>0)==(ydelta>0))
					pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
					//pred=topleft+((xdelta+ydelta)>>1);//average slope halved
				else
					pred=topleft+xdelta+ydelta;//average slope
#if 0
				char
					xdelta=((top-topleft)*3+(left-leftleft))>>2,
					ydelta=((left-topleft)*3+(top-toptop))>>2,
					xdelta2=top-topright;
				pred=(topleft+xdelta+ydelta+topright+xdelta2)>>1;
#endif

#if 0
				char vmax, vmin;
				if(top<left)
					vmin=top, vmax=left;
				else
					vmin=left, vmax=top;
#if 0
				//if(topleft>vmax&&vmin>topright)
				//	pred=top+left-topleft;
				//else
				if(topleft>vmax)
					pred=vmin;
				//else if(topleft<vmin&&vmax<topright)
				//{
				//	pred=top+left-topleft;
				//
				//	//pred=topleft+top-topleft+left-topleft;
				//
				//	//char xdelta=top-topleft, ydelta=left-topleft;
				//	//pred=topleft+xdelta+ydelta;
				//}
				else if(topleft<vmin)
					pred=vmax;
				else
					pred=top+left-topleft;
#endif

#if 1
				if(topleft>vmax)//choose steepest slope
					pred=vmin;
				else if(topleft<vmin)
					pred=vmax;
				else
					pred=left+top-topleft;//planar prediction (unplane)
#endif
#endif

				buf[idx]-=pred;
			}
		}
	}
}
void image_replane(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					
					leftleft=kx-2>=0?buf[idx-(bytestride<<1)]:0,
					toptop=ky-2>=0?buf[idx-(rowlen<<1)]:0,
					topright=kx<iw-1&&ky?buf[idx-rowlen+bytestride]:0,

					pred;

				char xdelta=top-topleft, ydelta=left-topleft;
				if((xdelta>0)==(ydelta>0))
					pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
					//pred=topleft+((xdelta+ydelta)>>1);//average slope halved
				else
					pred=topleft+xdelta+ydelta;//average slope
				
#if 0
				char
					xdelta=((top-topleft)*3+(left-leftleft))>>2,
					ydelta=((left-topleft)*3+(top-toptop))>>2,
					xdelta2=top-topright;
				pred=(topleft+xdelta+ydelta+topright+xdelta2)>>1;
#endif

#if 0
				char vmax, vmin;
				if(top<left)
					vmin=top, vmax=left;
				else
					vmin=left, vmax=top;
#if 0
				//if(topleft>vmax&&vmin>topright)
				//	pred=top+left-topleft;
				//else
				if(topleft>vmax)
					pred=vmin;
				//else if(topleft<vmin&&vmax<topright)
				//	pred=top+left-topleft;
				else if(topleft<vmin)
					pred=vmax;
				else
					pred=top+left-topleft;
#endif

#if 1
				if(topleft>vmax)
					pred=vmin;
				else if(topleft<vmin)
					pred=vmax;
				else
					pred=left+top-topleft;
#endif
#endif

				buf[idx]+=pred;
			}
		}
	}
}

void image_customst_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char comp[]=
				{
					kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
					kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)-bytestride]:0,
					         ky-2>=0?buf[idx-(rowlen<<1)]:0,
					kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+bytestride]:0,
					kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,

					kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
					kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride]:0,
					         ky-1>=0?buf[idx-rowlen]:0,
					kx+1<iw&&ky-1>=0?buf[idx-rowlen+bytestride]:0,
					kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,

					kx-2>=0?buf[idx-(bytestride<<1)]:0,
					kx-1>=0?buf[idx-bytestride]:0,
				};

				char pred=(char)(
					customparam_st[ 0]*comp[ 0]+
					customparam_st[ 1]*comp[ 1]+
					customparam_st[ 2]*comp[ 2]+
					customparam_st[ 3]*comp[ 3]+
					customparam_st[ 4]*comp[ 4]+
					customparam_st[ 5]*comp[ 5]+
					customparam_st[ 6]*comp[ 6]+
					customparam_st[ 7]*comp[ 7]+
					customparam_st[ 8]*comp[ 8]+
					customparam_st[ 9]*comp[ 9]+
					customparam_st[10]*comp[10]+
					customparam_st[11]*comp[11]
				);
				buf[idx]-=pred;
			}
		}
	}
}
void image_customst_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char comp[]=
				{
					kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
					kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)-bytestride]:0,
					         ky-2>=0?buf[idx-(rowlen<<1)]:0,
					kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+bytestride]:0,
					kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,

					kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
					kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride]:0,
					         ky-1>=0?buf[idx-rowlen]:0,
					kx+1<iw&&ky-1>=0?buf[idx-rowlen+bytestride]:0,
					kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,

					kx-2>=0?buf[idx-(bytestride<<1)]:0,
					kx-1>=0?buf[idx-bytestride]:0,
				};

				char pred=(char)(
					customparam_st[ 0]*comp[ 0]+
					customparam_st[ 1]*comp[ 1]+
					customparam_st[ 2]*comp[ 2]+
					customparam_st[ 3]*comp[ 3]+
					customparam_st[ 4]*comp[ 4]+
					customparam_st[ 5]*comp[ 5]+
					customparam_st[ 6]*comp[ 6]+
					customparam_st[ 7]*comp[ 7]+
					customparam_st[ 8]*comp[ 8]+
					customparam_st[ 9]*comp[ 9]+
					customparam_st[10]*comp[10]+
					customparam_st[11]*comp[11]
				);
				buf[idx]+=pred;
			}
		}
	}
}

ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	ArrayHandle sizes;
	int lw=floor_log2(iw), lh=floor_log2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	int nstages=0;
	DWTSize *p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=iw;
	p->h=ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=w2;
		p->h=h2;
	}
	return sizes;
}

//lifting-based 8bit lazy DWT
void dwt1d_lazy_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_lazy_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_lazy_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit Haar DWT
void dwt1d_haar_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd;++k)//const predictor (difference)
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd;++k)//update (O+(E-O)/2 = average)
		odd[k]+=even[k]>>1;


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_haar_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
		odd[k]-=even[k]>>1;
	
	for(int k=0;k<nodd;++k)//unpredict
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_haar_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit squeeze DWT
char smoothtendency(char B, char a, char n)
{
	char diff=0;
	if(B>=a&&a>=n)
	{
		diff=(4*B-3*n-a+6)/12;
		//      2C = a<<1 + diff - diff&1 <= 2B  so diff - diff&1 <= 2B - 2a
		//      2D = a<<1 - diff - diff&1 >= 2n  so diff + diff&1 <= 2a - 2n
		if(diff-(diff&1)>2*(B-a))
			diff=2*(B-a)+1;
		if(diff+(diff&1)>2*(a-n))
			diff=2*(a-n);
	}
	else if(B<=a&&a<=n)
	{
		diff=(4*B-3*n-a-6)/12;
		//      2C = a<<1 + diff + diff&1 >= 2B  so diff + diff&1 >= 2B - 2a
		//      2D = a<<1 - diff + diff&1 <= 2n  so diff - diff&1 >= 2a - 2n
		if(diff+(diff&1)<2*(B-a))
			diff=2*(B-a)-1;
		if(diff-(diff&1)<2*(a-n))
			diff=2*(a-n);
	}
	return diff;
}
void dwt1d_squeeze_fwd(char *buffer, int count, int stride, char *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		char
			o1=kd>0?buffer[ks-stride]:0,
			e =buffer[ks],
			o =(kd<<1)+1<count?buffer[ks+stride]:0,
			e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:0,//n-1-(idx-(n-1)) = ((n-1)<<1)-idx
			o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:0;

		e-=o;		//diff
		o+=e>>1;	//av
		e2-=o2;
		o2+=e2>>1;
		e-=smoothtendency(o1, o, o2);
		
		if(kd<nodd)
			odd[kd]=o;
		even[kd]=e;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_squeeze_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(int ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
	{
		char
			o1=ks>0?buffer[kd-stride]:0,
			e =even[ks],
			o=ks<nodd?odd[ks]:0,
			o2=ks+1<nodd?odd[ks+1]:0;

		e+=smoothtendency(o1, o, o2);
		o-=e>>1;
		e+=o;

		buffer[kd]=e;
		buffer[kd+stride]=o;
	}
}
void dwt2d_squeeze_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_squeeze_fwd(buffer+stride*kx, h2, rowlen, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dB.PNG", it);
	}
}
void dwt2d_squeeze_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_squeeze_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 5/3 DWT
void dwt1d_cdf53_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//linear predictor
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;

	//orginal CDF 5/3 DWT
#if 0
	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//linear predictor
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;
#endif


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf53_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf53_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf53_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	-0x1960C,	//-1.58613434342059f,	//alpha
	-0x00D90,	//-0.0529801185729f,	//beta
	 0x0E206,	// 0.8829110755309f,	//gamma
	 0x07189,	// 0.4435068520439f,	//delta
				// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_update(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unupdate(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
void dwt1d_cdf97_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[3]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

static unsigned qhist[256]={0};
void channel_entropy(unsigned char *buf, int resolution, int nch, int bytestride, float *cr, int *usage)
{
	double entropy[4]={0};
	memset(usage, 0, 4*sizeof(int));
	for(int kc=0;kc<nch;++kc)
	{
		memset(qhist, 0, 256*sizeof(unsigned));
		for(int k=0, end=resolution*bytestride;k<end;k+=bytestride)
		{
			unsigned char val=buf[k+kc];
			++qhist[val];
		}
		for(int ks=0;ks<256;++ks)
		{
			unsigned freq=qhist[ks];
			if(freq)
			{
				double p=(double)freq/resolution;
				entropy[kc]-=p*log2(p);
				++usage[kc];
			}
		}
		cr[kc]=(float)(8/entropy[kc]);
	}
	
	//calculate csize with joint histogram
	unsigned *h2=(unsigned*)malloc(0x1000000*sizeof(unsigned));
	if(!h2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(h2, 0, 0x1000000*sizeof(unsigned));
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		++h2[color];
	}
	for(int k=0;k<0x1000000;++k)
		usage[3]+=h2[k]!=0;
	double csize=0;
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		double p=(double)h2[color]/resolution;
		p*=0xFFFFFFFF;
		++p;
		p/=0x100000000;
		double bitsize=-log2(p);
		csize+=bitsize;
	}
	free(h2);
	csize/=8;
	cr[3]=(float)(resolution*3/csize);
}
void jointhistogram(unsigned char *buf, int resolution, int nbits, ArrayHandle *hist)
{
	int nlevels=1<<nbits, hsize=nlevels*nlevels*nlevels;
	if(*hist)
		array_free(hist);
	ARRAY_ALLOC(unsigned, *hist, 0, hsize, 0, 0);
	//unsigned *hist=(unsigned*)malloc(hsize*sizeof(unsigned));
	unsigned *h=(unsigned*)hist[0]->data;
	for(int k=0;k<resolution;++k)
	{
		unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
		int color=b<<(nbits<<1)|g<<nbits|r;

		++h[color];
	}

	//X  don't calculate csize from downsampled histogram
#if 0
	double csize=0;
	for(int k=0;k<resolution;++k)//calculate csize using Zipf's law
	{
		unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
		int color=b<<(nbits<<1)|g<<nbits|r;

		unsigned freq=h[color];
		double p=(double)freq/resolution;

		p*=0xFFFFFFFF;//guard against quantized zero
		++p;
		p/=0x100000000;

		csize-=log2(p);
	}
	csize/=nbits;
	float CR=(float)(resolution*3/csize);
#endif

#if 0
	float entropy=0;//also BPP
	for(int k=0;k<hsize;++k)//calculate entropy using Shannon's law
	{
		int freq=h[k];
		if(freq)
		{
			float p=(float)freq/resolution;
			entropy-=p*log2f(p);
		}
	}
	float CR=3*nbits/entropy;
#endif

	unsigned histmin=0;
	unsigned histmax=0;
	for(int k=0;k<hsize;++k)//get min & max
	{
		if(histmin>h[k])
			histmin=h[k];
		if(histmax<h[k])
			histmax=h[k];
	}

	unsigned char *h2=(unsigned char*)h;
	for(int k=0;k<hsize;++k)//normalize
		h2[k]=(h[k]<<8)/histmax;

	//for(int k=0;k<hsize;++k)
	//	h[k]=(unsigned)(((long long)h[k]<<32)/histmax);

	//for(int kz=0;kz<nlevels-1;++kz)
	//{
	//	for(int ky=0;ky<nlevels-1;++ky)
	//	{
	//		for(int kx=0;kx<nlevels-1;++kx)
	//		{
	//		}
	//	}
	//}

	//return CR;
}