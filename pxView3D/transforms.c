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

#if 0
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
#endif
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
double customparam_ct[12]={0}, customparam_st[12]={0};
int customparam_clamp[2]={-128, 127};
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

//learned predictor
#if 1
void mulmatvec_pd(double *dst, const double *matrix, const double *vector, int mw, int mh)
{
	for(int ky=0;ky<mh;++ky)
	{
		double sum=0;
		for(int kx=0;kx<mw;++kx)
			sum+=matrix[mw*ky+kx]*vector[kx];
		dst[ky]=sum;
	}
}
void mulvTmat_pd(double *dst, const double *vT, const double *matrix, int mw, int mh)
{
	for(int kx=0;kx<mw;++kx)
	{
		double sum=0;
		for(int ky=0;ky<mh;++ky)
			sum+=vT[ky]*matrix[mw*ky+kx];
		dst[kx]=sum;
	}
}
void mulvecs_pd(double *dst, const double *a, const double *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]*b[k];
}
void scalevec_pd(double *vec, int count, double factor)
{
	for(int k=0;k<count;++k)
		vec[k]*=factor;
}
void filternan_pd(double *vec, int count)
{
	for(int k=0;k<count;++k)
	{
		if(!isfinite(vec[k])||fabs(vec[k])>1e6)
			vec[k]=0;
	}
}
void addvecs_pd(double *dst, const double *a, const double *b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b[k];
}
void addvec1_pd(double *dst, const double *a, const double b, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=a[k]+b;
}
double vecsum_pd(const double *data, int count)
{
	double sum=0;
	for(int k=0;k<count;++k)
		sum+=data[k];
	return sum;
}
void negbuffer_pd(double *data, int count)
{
	for(int k=0;k<count;++k)
		data[k]=-data[k];
}
void subvecs_pd(double *dst, const double *pos, const double *neg, int count)
{
	for(int k=0;k<count;++k)
		dst[k]=pos[k]-neg[k];
}
void leakyReLU_pd(double *dst, const double *src, int count)
{
	for(int k=0;k<count;++k)
	{
		if(src[k]<0)
			dst[k]=src[k]*0.01;
		else
			dst[k]=src[k];
	}
}
void leakyReLUdash_pd(double *data, int count)
{
	for(int k=0;k<count;++k)
	{
		if(data[k]<0)
			data[k]=0.01;
		else
			data[k]=1;
	}
}
double calc_rmse_pd(const double *error, int count)
{
	double sum=0;
	for(int k=0;k<count;++k)
		sum+=error[k]*error[k];
	sum/=count;
	sum=sqrt(sum);
	return sum;
}
void fill_matrow_unchecked(const unsigned char *buf, int iw, int ih, int kc, int kx, int ky, double *mrow, double *target)
{
	int idx=iw*ky+kx;
	mrow[0]=buf[(idx-iw*2-2)<<2|kc];
	mrow[1]=buf[(idx-iw*2-1)<<2|kc];
	mrow[2]=buf[(idx-iw*2)<<2|kc];
	mrow[3]=buf[(idx-iw*2+1)<<2|kc];
	mrow[4]=buf[(idx-iw*2+2)<<2|kc];

	mrow[5]=buf[(idx-iw-2)<<2|kc];
	mrow[6]=buf[(idx-iw-1)<<2|kc];
	mrow[7]=buf[(idx-iw)<<2|kc];
	mrow[8]=buf[(idx-iw+1)<<2|kc];
	mrow[9]=buf[(idx-iw+2)<<2|kc];

	mrow[10]=buf[(idx-2)<<2|kc];
	mrow[11]=buf[(idx-1)<<2|kc];

	*target=buf[idx<<2|kc];
}
double leakyReLU1(double x)
{
	if(x<0)
		x*=0.01;
	return x;
}
#define OPT_N 12
#define OPT_B 12
static double g_mat[OPT_N*OPT_B], g_y[OPT_N], g_v[OPT_B], g_v2[OPT_B], g_grad[OPT_N+1];
double opt_causal_reach2(unsigned char *buf, int iw, int ih, int kc, double *x, double *bias, double lr, int test)
{
	addhalf(buf, iw, ih, 3, 4);
	int nupdates=0;
	double rmse=0;
	for(int ky=2, row=0;ky<ih-2;ky+=3)
	{
		for(int kx=2;kx<iw-4;kx+=5)
		{
			fill_matrow_unchecked(buf, iw, ih, kc, kx, ky, g_mat+row*12, g_y+row);
			row=(row+1)%12;
			if(!row)
			{
#if 0
				//forward		pred = leakyReLU(M x + b)		x & b are params, M & y are pixels
				mulmatvec_pd(g_v,  g_mat, x,  OPT_N, OPT_B);	//     M x
				addvec1_pd(g_v,  g_v, *bias,  OPT_N);			//v  = M x + b
				leakyReLU_pd(g_v2,  g_v,  OPT_N);				//     leakyReLU(M x + b)
				subvecs_pd(g_v2,  g_y, g_v2,  OPT_N);			//v2 = y - leakyReLU(M x + b)
				
				rmse+=calc_rmse_pd(g_v2, OPT_N);				//L  = sum i: (y[i] - leakyReLU(M x[i] + b))^2

				if(!test)
				{
					//backward
					leakyReLUdash_pd(g_v, OPT_N);				//v  = leakyReLU'(M x + b)
					negbuffer_pd(g_v, OPT_N);					//v  = -leakyReLU'(M x + b)

					mulvecs_pd(g_v2, g_v2, g_v, OPT_N);			//v2 = 2 * (y - leakyReLU(M x + b)) * -leakyReLU'(M x + b)
					mulvTmat_pd(g_grad, g_v2, g_mat, OPT_N, OPT_N);	//gx = v2T M
					g_grad[OPT_N]=vecsum_pd(g_v2, OPT_N);			//gb = sum(v2)
				
					scalevec_pd(g_grad, OPT_N+1, lr);

					subvecs_pd(x, x, g_grad, OPT_N);
					*bias-=g_grad[OPT_N];
				}
				++nupdates;
#endif
#if 1
				mulmatvec_pd(g_v, g_mat, x, OPT_N, OPT_N);
				subvecs_pd(g_v, g_y, g_v, OPT_N);
				rmse+=calc_rmse_pd(g_v, OPT_N);
				mulvTmat_pd(g_grad, g_v, g_mat, OPT_N, OPT_N);
				scalevec_pd(g_grad, OPT_N, lr);
				addvecs_pd(x, x, g_grad, OPT_N);
				filternan_pd(x, OPT_N);
				++nupdates;
#endif
			}
		}
	}
	addhalf(buf, iw, ih, 3, 4);
	rmse/=nupdates;
	return rmse;
}

#define O2_N 72
double customparam_hybrid[(7*3+3)*3*3]={0};//shape [3, 72]		maps 24 causal neighbors (72 values) -> 1 pixel (3 values)
double opt_causal_hybrid_r3(unsigned char *buf, int iw, int ih, double lr)
{
	double *x=customparam_hybrid;
	ArrayHandle tensors;
	ARRAY_ALLOC(double, tensors, 0, 0, 0, 0);
	int off_src=(int)ARRAY_APPEND_OFFSET(tensors, 0, O2_N, 1, 0);
	int off_grad=(int)ARRAY_APPEND_OFFSET(tensors, 0, O2_N*3, 1, 0);

	double *t=(double*)tensors->data;
	
	addhalf(buf, iw, ih, 3, 4);
	int nupdates=0;
	double rmse=0;
	for(int ky=3, row=0;ky<ih-3;ky+=4)
	{
		for(int kx=3;kx<iw-6;kx+=7)
		{
			int idx=iw*ky+kx;
			for(int ky2=0, kd=0;ky2<4;++ky2)//fetch causal neighbors
			{
				for(int kx2=0;kx2<7&&kd<O2_N;++kx2)
				{
					for(int kc=0;kc<3;++kc, ++kd)
						t[off_src+kd]=(char)buf[(idx+iw*(ky2-3)+kx2-3)<<2|kc];
				}
			}
			double v[3]={(char)buf[idx<<2], (char)buf[idx<<2|1], (char)buf[idx<<2|2]};

			for(int kc=0;kc<3;++kc)//v[3] = y[3] - x[3, 72] src.flatten[72]
			{
				for(int k=0;k<O2_N;++k)
					v[kc]-=x[kc*O2_N+k]*t[off_src+k];
			}

			rmse+=calc_rmse_pd(v, 3);

			for(int kc=0;kc<3;++kc)//grad = v[3] * -src.flatten[72]T		outer product
			{
				for(int k=0;k<O2_N;++k)
					t[off_grad+kc*O2_N+k]=v[kc]*-t[off_src+k];
			}
			
			scalevec_pd(t+off_grad, O2_N*3, lr);
			subvecs_pd(x, x, t+off_grad, O2_N*3);
			filternan_pd(x, O2_N*3);

			++nupdates;
		}
	}
	array_free(&tensors);
	rmse/=nupdates;
	addhalf(buf, iw, ih, 3, 4);
	return rmse;
}
void predict_hybrid(const char *buf, int iw, int kx, int ky, int idx, int rowlen, int *ret_pred)
{
	char cn[72];
	for(int ky2=0, kd=0;ky2<4;++ky2)//fetch causal neighbors
	{
		for(int kx2=0;kx2<7&&kd<O2_N;++kx2)
		{
			for(int kc=0;kc<3;++kc, ++kd)
				cn[kd]=ky2-3>=0&&BETWEEN_EXC(3, kx2, iw+3)>=0?buf[(idx+((iw*(ky2-3)+kx2-3)<<2))|kc]:0;
		}
	}
	for(int kc=0;kc<3;++kc)//v[3] = y[3] - x[3, 72] src.flatten[72]
	{
		double sum=0;
		for(int k=0;k<O2_N;++k)
			sum+=customparam_hybrid[kc*O2_N+k]*cn[k];
		sum=CLAMP(customparam_clamp[0], sum, customparam_clamp[1]);
		ret_pred[kc]=(int)round(sum);
	}
}
void pred_hybrid_fwd(char *buf, int iw, int ih)
{
	int rowlen=iw<<2;
	int idx=(iw*ih-1)<<2;
	for(int ky=ih-1;ky>=0;--ky)
	{
		for(int kx=iw-1;kx>=0;--kx, idx-=4)
		{
			int pred[3];
			predict_hybrid(buf, iw, kx, ky, idx, rowlen, pred);
			
			buf[idx  ]-=pred[0];
			buf[idx|1]-=pred[1];
			buf[idx|2]-=pred[2];
		}
	}
}
void pred_hybrid_inv(char *buf, int iw, int ih)
{
	int rowlen=iw<<2;
	int idx=0;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred[3];
			predict_hybrid(buf, iw, kx, ky, idx, rowlen, pred);

			buf[idx  ]+=pred[0];
			buf[idx|1]+=pred[1];
			buf[idx|2]+=pred[2];
		}
	}
}
#endif

#define PREDICT_LINEAR(A, B)		((A)<<1)-(B)
#define PREDICT_QUADRATIC(A, B, C)	((A)-3*(B)+3*(C))
#if 0
static char predict_linear(char prev2, char prev)
{
	int pred=(prev<<1)-prev2;
	//if(pred<-128)
	//	pred=-128;
	//if(pred>127)
	//	pred=127;
	return pred;

	//return prev2+((prev-prev2)<<1);
}
static int predict_quadratic(char prev3, char prev2, char prev)
{
	int pred=prev3-3*prev2+3*prev;
	//if(pred<-128)
	//	pred=-128;
	//if(pred>127)
	//	pred=127;
	return pred;
}
#endif
static int predict_simple(char topleft, char top, char left)
{
	int xdelta=top-topleft, ydelta=left-topleft, pred;
	if((xdelta>0)==(ydelta>0))
		pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
	else
		pred=topleft+xdelta+ydelta;//average slope
	return pred;
}
//int testhist[3]={0};//
int  predict_diff2d(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char cn[]=
	{
		kx-1>=0&&ky-1>=0?buf[idx-rowlen-bytestride]:0,
		         ky-1>=0?buf[idx-rowlen           ]:0,
		
		kx-1>=0?buf[idx-bytestride]:0,
	};
	int pred=cn[1]+cn[2]-cn[0];
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;

#if 0
	char cn[]=
	{
		kx-3>=0&&ky-3>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
		kx-2>=0&&ky-3>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-3>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-3>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-3>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-3>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
		kx+3<iw&&ky-3>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
		kx+3<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
		kx-3>=0&&ky-1>=0?buf[idx-rowlen- bytestride*3  ]:0,
		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,
		kx+3<iw&&ky-1>=0?buf[idx-rowlen+ bytestride*3  ]:0,
		
		kx-3>=0?buf[idx- bytestride*3  ]:0,
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};

	int predx0=cn[23], predy0=cn[17];
	int predx1=PREDICT_LINEAR(cn[22], cn[23]), predy1=PREDICT_LINEAR(cn[10], cn[17]);
	int predx2=PREDICT_QUADRATIC(cn[21], cn[22], cn[23]), predy2=PREDICT_QUADRATIC(cn[3], cn[10], cn[17]);
	int pred=(predx0+predy0)>>1;
	int best=abs(predx0-predy0),
		test=abs(predx1-predy1);
	if(best>test)
		best=test, pred=(predx1+predy1)>>1;
	test=abs(predx2-predy2);
	if(best>test)
		pred=(predx2+predy2)>>1;

	if(pred==(predx0+predy0)>>1)
		++testhist[0];
	else if(pred==(predx1+predy1)>>1)
		++testhist[1];
	else if(pred==(predx2+predy2)>>1)
		++testhist[2];

	if(pred<-128)
		pred=-128;
	if(pred>127)
		pred=127;

#endif
#if 0
	char cn[]=
	{
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,

		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,

		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	int error=
		predict_simple(cn[0], cn[1], cn[5])-cn[6]+
		predict_simple(cn[1], cn[2], cn[6])-cn[7]+
		predict_simple(cn[4], cn[3], cn[9])-cn[8]+
		predict_simple(cn[5], cn[6], cn[10])-cn[11];
	int pred=predict_simple(cn[7], cn[8], cn[11]);
	pred=pred-(pred*error>>9);//10
#endif

	//if(pred<-64)
	//	pred=-64;
	//if(pred>64)
	//	pred=64;
	//return pred;
}
void pred_diff2d_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char pred=predict_diff2d(buf, iw, kx, ky, idx, bytestride, rowlen);
#if 0
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=left+top-topleft;
				//if(kx||ky)
				//	pred-=128;
#endif
				buf[idx]-=pred;
			}
		}
	}
}
void pred_diff2d_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_diff2d(buf, iw, kx, ky, idx, bytestride, rowlen);
#if 0
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					pred=left+top-topleft;
				//if(kx||ky)
				//	pred-=128;
#endif
				buf[idx]+=pred;
			}
		}
	}
}

void pred_hpf_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
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
					pred=(left+top)>>1;

				pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				buf[idx]-=pred;
			}
		}
	}
}
void pred_hpf_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
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
					pred=(left+top)>>1;

				pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
				buf[idx]+=pred;
			}
		}
	}
}

int  predict_median(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char cn[]=
	{
	//	kx-3>=0&&ky-3>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
	//	kx-2>=0&&ky-3>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
	//	kx-1>=0&&ky-3>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
	//	         ky-3>=0?buf[idx-(rowlen<<1)                ]:0,
	//	kx+1<iw&&ky-3>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
	//	kx+2<iw&&ky-3>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-3>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
	//	kx-3>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride*3  ]:0,
		kx-2>=0&&ky-2>=0?buf[idx-(rowlen<<1)-(bytestride<<1)]:0,
		kx-1>=0&&ky-2>=0?buf[idx-(rowlen<<1)- bytestride    ]:0,
		         ky-2>=0?buf[idx-(rowlen<<1)                ]:0,
		kx+1<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride    ]:0,
		kx+2<iw&&ky-2>=0?buf[idx-(rowlen<<1)+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-2>=0?buf[idx-(rowlen<<1)+ bytestride*3  ]:0,
		
	//	kx-3>=0&&ky-1>=0?buf[idx-rowlen- bytestride*3  ]:0,
		kx-2>=0&&ky-1>=0?buf[idx-rowlen-(bytestride<<1)]:0,
		kx-1>=0&&ky-1>=0?buf[idx-rowlen- bytestride    ]:0,
		         ky-1>=0?buf[idx-rowlen                ]:0,
		kx+1<iw&&ky-1>=0?buf[idx-rowlen+ bytestride    ]:0,
		kx+2<iw&&ky-1>=0?buf[idx-rowlen+(bytestride<<1)]:0,
	//	kx+3<iw&&ky-1>=0?buf[idx-rowlen+ bytestride*3  ]:0,
		
	//	kx-3>=0?buf[idx- bytestride*3  ]:0,
		kx-2>=0?buf[idx-(bytestride<<1)]:0,
		kx-1>=0?buf[idx- bytestride    ]:0,
	};
	for(int k=1;k<COUNTOF(cn);++k)//insertion sort
	{
		char val=cn[k];
		int L=0, R=k-1;
		int idx=-1;
		while(L<=R)
		{
			int mid=(L+R)>>1;
			if(cn[mid]<val)
				L=mid+1;
			else if(cn[mid]>val)
				R=mid-1;
			else
			{
				idx=mid;
				break;
			}
		}
		if(idx==-1)
			idx=L+(L<k&&cn[L]<val);
		if(idx<k)
		{
			char temp=cn[k];
			memmove(cn+idx+1, cn+idx, (size_t)k-idx);
			cn[idx]=temp;
		}
	}
	int pred=(cn[(COUNTOF(cn)>>1)-1]+cn[COUNTOF(cn)>>1])>>1;

	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_median_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char pred=predict_median(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_median_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	//memset(testhist, 0, sizeof(testhist));//
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_median(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}

int  predict_grad(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
{
	char
		left=kx?buf[idx-bytestride]:0,
		top=ky?buf[idx-rowlen]:0,
		topleft=kx&&ky?buf[idx-rowlen-bytestride]:0;

	int pred;

	char vmax, vmin;
	if(top<left)
		vmin=top, vmax=left;
	else
		vmin=left, vmax=top;

	if(topleft>vmax)//choose steepest slope if both downward or upward
		pred=vmin;
	else if(topleft<vmin)
		pred=vmax;
	else
		pred=left+top-topleft;//planar prediction (unplane)

	//char xdelta=top-topleft, ydelta=left-topleft;
	//if((xdelta>0)==(ydelta>0))
	//	pred=topleft+(abs(xdelta)>abs(ydelta)?xdelta:ydelta);//take steepest slope once and stop, equivalent to original unplane
	//else
	//	pred=topleft+xdelta+ydelta;//average slope
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_grad_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_grad(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_grad_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_grad(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}

int  predict_custom(const char *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen)
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

		//leakyReLU1(//

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

		//+customparam_ct[11])//
	);
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return pred;
}
void pred_custom_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char pred=predict_custom(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_custom_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char pred=predict_custom(buf, iw, kx, ky, idx, bytestride, rowlen);

				buf[idx]+=pred;
			}
		}
	}
}


//DWTs
void dwt2d_grad_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	//for(int it=sizes_start;it<1;++it)
	{
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				char v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				//if((unsigned char)v[0]==0xFF||(unsigned char)v[1]==0xFF||(unsigned char)v[2]==0xFF||(unsigned char)v[3]==0xFF)
				//	v[0]=0xFF;

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]-=vmax;
				else if(v[0]>vmax)
					v[3]-=vmin;
				else
					v[3]-=v[1]+v[2]-v[0];

				v[2]-=v[0];
				v[1]-=v[0];
				v[0]+=(v[1]+v[2])>>2;

				buffer[ idx            *stride]=v[3];//grad
				buffer[(idx         +1)*stride]=v[2];//diffy
				buffer[(idx+sizes->w  )*stride]=v[1];//diffx
				buffer[(idx+sizes->w+1)*stride]=v[0];//av

				//buffer[ idx            <<2|kc]=v[0];//av
				//buffer[(idx         +1)<<2|kc]=v[1];//diffx
				//buffer[(idx+sizes->w  )<<2|kc]=v[2];//diffy
				//buffer[(idx+sizes->w+1)<<2|kc]=v[3];//grad
			}
		}
		dwt2d_lazy_fwd(buffer, sizes, it, it+2, stride, temp);
	}
}
void dwt2d_grad_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		dwt2d_lazy_inv(buffer, sizes, it, it+2, stride, temp);
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				char v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				v[0]-=(v[1]+v[2])>>2;
					
				v[1]+=v[0];
				v[2]+=v[0];

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]+=vmax;
				else if(v[0]>vmax)
					v[3]+=vmin;
				else
					v[3]+=v[1]+v[2]-v[0];

				buffer[ idx            *stride]=v[3];
				buffer[(idx         +1)*stride]=v[2];
				buffer[(idx+sizes->w  )*stride]=v[1];
				buffer[(idx+sizes->w+1)*stride]=v[0];
			}
		}
	}
}

void dwt1d_custom_fwd(char *buffer, int count, int stride, char *b2)
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


	even[0]-=(char)floor(odd[0]*(customparam_st[0]+customparam_st[5]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*customparam_st[0]+odd[k]*customparam_st[5]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(customparam_st[0]+customparam_st[5]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*customparam_st[1]+even[k+1]*customparam_st[6]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(customparam_st[1]+customparam_st[6]));


	even[0]-=(char)(odd[0]*(customparam_st[2]+customparam_st[7]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*customparam_st[2]+odd[k]*customparam_st[7]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(customparam_st[2]+customparam_st[7]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*customparam_st[3]+even[k+1]*customparam_st[8]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(customparam_st[3]+customparam_st[8]));


	even[0]-=(char)floor(odd[0]*(customparam_st[4]+customparam_st[9]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*customparam_st[4]+odd[k]*customparam_st[9]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(customparam_st[4]+customparam_st[9]));


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_custom_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	even[0]+=(char)floor(odd[0]*(customparam_st[4]+customparam_st[9]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*customparam_st[4]+odd[k]*customparam_st[9]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(customparam_st[4]+customparam_st[9]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*customparam_st[3]+even[k+1]*customparam_st[8]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(customparam_st[3]+customparam_st[8]));
	
	even[0]+=(char)floor(odd[0]*(customparam_st[2]+customparam_st[7]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*customparam_st[2]+odd[k]*customparam_st[7]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(customparam_st[2]+customparam_st[7]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*customparam_st[1]+even[k+1]*customparam_st[6]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(customparam_st[1]+customparam_st[6]));
	
	even[0]+=(char)floor(odd[0]*(customparam_st[0]+customparam_st[5]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*customparam_st[0]+odd[k]*customparam_st[5]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(customparam_st[0]+customparam_st[5]));


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_custom_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_custom_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_custom_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_custom_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_custom_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_custom_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

void dwt1d_exp_fwd(char *buffer, int count, int stride, char *b2)
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


	for(int k=0;k<nodd+extraeven;++k)//predict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			customparam_st[0]*(prev+next)+
			customparam_st[1]*(prev2+next2)+
			customparam_st[2]*(prev3+next3)+
			customparam_st[3]*(prev4+next4)+
			customparam_st[4]*(prev5+next5)
		));
		even[k]-=pred;
		//even[k]-=(9*(prev+next)+prev2+next2)>>4;
	}
	for(int k=0;k<nodd;++k)//update
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			customparam_st[5]*(prev+next)+
			customparam_st[6]*(prev2+next2)+
			customparam_st[7]*(prev3+next3)+
			customparam_st[8]*(prev4+next4)+
			customparam_st[9]*(prev5+next5)
		));
		odd[k]+=update;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_exp_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			customparam_st[5]*(prev+next)+
			customparam_st[6]*(prev2+next2)+
			customparam_st[7]*(prev3+next3)+
			customparam_st[8]*(prev4+next4)+
			customparam_st[9]*(prev5+next5)
		));
		odd[k]-=update;
	}
	for(int k=0;k<nodd+extraeven;++k)//unpredict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			customparam_st[0]*(prev+next)+
			customparam_st[1]*(prev2+next2)+
			customparam_st[2]*(prev3+next3)+
			customparam_st[3]*(prev4+next4)+
			customparam_st[4]*(prev5+next5)
		));
		even[k]+=pred;
		//even[k]+=(9*(prev+next)+prev2+next2)>>4;
	}


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_exp_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_exp_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_exp_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
}
void dwt2d_exp_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_exp_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_exp_inv(buffer+rowlen*ky, w2, stride, temp);
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
void dwt2d_lazy_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_lazy_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_haar_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_haar_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
		if(ks<nodd)
			buffer[kd+stride]=o;
	}
}
void dwt2d_squeeze_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_cdf53_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_cdf53_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	-0x1960C,	//-1.58613434342059f,	//alpha
	-0x00D90,	//-0.0529801185729f,	//beta
	 0x0E206,	// 0.8829110755309f,	//gamma
	 0x07189,	// 0.4435068520439f,	//delta
	 0x1264C,	// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
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
static void dwt1d_u8_unupdate(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_scale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=buf[k]*coeff>>16;
}
static void dwt1d_u8_unscale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=(buf[k]<<16)/coeff;
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

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[3]);
	//dwt1d_u8_scale(b2, count, cdf97_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	//dwt1d_u8_unscale(b2, count, cdf97_coeff[4]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
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
void dwt2d_cdf97_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

void dwt2d_dec_fwd(char *buffer, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	int w2=iw, h2=ih;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	free(temp);
}
void dwt2d_dec_inv(char *buffer, int iw, int ih)
{
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 3);
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[2].w;++kx)//vertical invDWT
			dwt1d_haar_inv(buffer+4*kx+kc, psizes[2].h, 4*iw, temp);
		for(int ky=0;ky<psizes[2].h;++ky)//horizontal invDWT
			dwt1d_haar_inv(buffer+4*iw*ky+kc, psizes[2].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[1].w;++kx)//vertical invDWT
			dwt1d_cdf53_inv(buffer+4*kx+kc, psizes[1].h, 4*iw, temp);
		for(int ky=0;ky<psizes[1].h;++ky)//horizontal invDWT
			dwt1d_cdf53_inv(buffer+4*iw*ky+kc, psizes[1].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[0].w;++kx)//vertical invDWT
			dwt1d_cdf97_inv(buffer+4*kx+kc, psizes[0].h, 4*iw, temp);
		for(int ky=0;ky<psizes[0].h;++ky)//horizontal invDWT
			dwt1d_cdf97_inv(buffer+4*iw*ky+kc, psizes[0].w, 4, temp);
	}
	free(temp);
	array_free(&sizes);
}


//DCTs

static void dct4_fwd_i8(char *x)
{
	x[3]=x[0]-x[3];
	x[0]-=x[3]>>1;
	x[2]=x[1]-x[2];
	x[1]-=x[2]>>1;
	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;
	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;
}
static void dct4_inv_i8(char *x)
{
	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];
	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];
	x[1]+=x[2]>>1;
	x[2]=x[1]-x[2];
	x[0]+=x[3]>>1;
	x[3]=x[0]-x[3];
}
void image_dct4_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[ kx>>2           ]=x[0];
				temp[(kx>>2)+(iw>>2)  ]=x[1];
				temp[(kx>>2)+(iw>>2)*2]=x[2];
				temp[(kx>>2)+(iw>>2)*3]=x[3];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[(ky>>2)          ]=x[0];
				temp[(ky>>2)+(ih>>2)  ]=x[1];
				temp[(ky>>2)+(ih>>2)*2]=x[2];
				temp[(ky>>2)+(ih>>2)*3]=x[3];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct4_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>2)          ],
					temp[(ky>>2)+(ih>>2)  ],
					temp[(ky>>2)+(ih>>2)*2],
					temp[(ky>>2)+(ih>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-3;kx+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>2)          ],
					temp[(kx>>2)+(iw>>2)  ],
					temp[(kx>>2)+(iw>>2)*2],
					temp[(kx>>2)+(iw>>2)*3],
				};

				dct4_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
}

static void dct8_fwd_i8(char *x)
{
	//binDCT-C7
	x[7]=x[0]-x[7];
	x[6]=x[1]-x[6];
	x[5]=x[2]-x[5];
	x[4]=x[3]-x[4];
	x[0]-=x[7]>>1;
	x[1]-=x[6]>>1;
	x[2]-=x[5]>>1;
	x[3]-=x[4]>>1;

	x[3]=x[0]-x[3];
	x[2]=x[1]-x[2];
	x[0]-=x[3]>>1;
	x[1]-=x[2]>>1;

	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;

	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;

	x[5]-=x[6]*13>>5;
	x[6]+=x[5]*11>>4;
	x[5]=(x[6]*15>>5)-x[5];

	x[5]=x[4]-x[5];
	x[6]=x[7]-x[6];
	x[4]-=x[5]>>1;
	x[7]-=x[6]>>1;

	x[4]=(x[7]*3>>4)-x[4];
	x[7]-=x[4]*3>>4;

	x[5]+=x[6]*11>>4;
	x[6]-=x[5]*15>>5;
}
static void dct8_inv_i8(char *x)
{
	//invBinDCT-C7
	x[6]+=x[5]*15>>5;
	x[5]-=x[6]*11>>4;

	x[7]+=x[4]*3>>4;
	x[4]=(x[7]*3>>4)-x[4];

	x[7]+=x[6]>>1;
	x[4]+=x[5]>>1;
	x[6]=x[7]-x[6];
	x[5]=x[4]-x[5];

	x[5]=(x[6]*15>>5)-x[5];
	x[6]-=x[5]*11>>4;
	x[5]+=x[6]*13>>5;

	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];

	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];

	x[1]+=x[2]>>1;
	x[0]+=x[3]>>1;
	x[2]=x[1]-x[2];
	x[3]=x[0]-x[3];

	x[3]+=x[4]>>1;
	x[2]+=x[5]>>1;
	x[1]+=x[6]>>1;
	x[0]+=x[7]>>1;
	x[4]=x[3]-x[4];
	x[5]=x[2]-x[5];
	x[6]=x[1]-x[6];
	x[7]=x[0]-x[7];
}
void image_dct8_fwd(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
					image[(idx+2)<<2|kc],
					image[(idx+3)<<2|kc],
					image[(idx+4)<<2|kc],
					image[(idx+5)<<2|kc],
					image[(idx+6)<<2|kc],
					image[(idx+7)<<2|kc],
				};

				//char y[8];
				//memcpy(y, x, 8);
				//dct8_fwd_i8(y);
				//dct8_inv_i8(y);
				//if(memcmp(x, y, 8))
				//	x[0]=y[0];

				dct8_fwd_i8(x);

				temp[ kx>>3           ]=x[0];
				temp[(kx>>3)+(iw>>3)  ]=x[1];
				temp[(kx>>3)+(iw>>3)*2]=x[2];
				temp[(kx>>3)+(iw>>3)*3]=x[3];
				temp[(kx>>3)+(iw>>3)*4]=x[4];
				temp[(kx>>3)+(iw>>3)*5]=x[5];
				temp[(kx>>3)+(iw>>3)*6]=x[6];
				temp[(kx>>3)+(iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
					image[(idx+iw*2)<<2|kc],
					image[(idx+iw*3)<<2|kc],
					image[(idx+iw*4)<<2|kc],
					image[(idx+iw*5)<<2|kc],
					image[(idx+iw*6)<<2|kc],
					image[(idx+iw*7)<<2|kc],
				};

				dct8_fwd_i8(x);

				temp[(ky>>3)          ]=x[0];
				temp[(ky>>3)+(ih>>3)  ]=x[1];
				temp[(ky>>3)+(ih>>3)*2]=x[2];
				temp[(ky>>3)+(ih>>3)*3]=x[3];
				temp[(ky>>3)+(ih>>3)*4]=x[4];
				temp[(ky>>3)+(ih>>3)*5]=x[5];
				temp[(ky>>3)+(ih>>3)*6]=x[6];
				temp[(ky>>3)+(ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
}
void image_dct8_inv(char *image, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(temp, 0, MAXVAR(iw, ih));
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih;++ky)
				temp[ky]=image[(iw*ky+kx)<<2|kc];
			for(int ky=0;ky<ih-7;ky+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(ky>>3)          ],
					temp[(ky>>3)+(ih>>3)  ],
					temp[(ky>>3)+(ih>>3)*2],
					temp[(ky>>3)+(ih>>3)*3],
					temp[(ky>>3)+(ih>>3)*4],
					temp[(ky>>3)+(ih>>3)*5],
					temp[(ky>>3)+(ih>>3)*6],
					temp[(ky>>3)+(ih>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+iw  )<<2|kc]=x[1];
				image[(idx+iw*2)<<2|kc]=x[2];
				image[(idx+iw*3)<<2|kc]=x[3];
				image[(idx+iw*4)<<2|kc]=x[4];
				image[(idx+iw*5)<<2|kc]=x[5];
				image[(idx+iw*6)<<2|kc]=x[6];
				image[(idx+iw*7)<<2|kc]=x[7];
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
				temp[kx]=image[(iw*ky+kx)<<2|kc];
			for(int kx=0;kx<iw-7;kx+=8)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					temp[(kx>>3)          ],
					temp[(kx>>3)+(iw>>3)  ],
					temp[(kx>>3)+(iw>>3)*2],
					temp[(kx>>3)+(iw>>3)*3],
					temp[(kx>>3)+(iw>>3)*4],
					temp[(kx>>3)+(iw>>3)*5],
					temp[(kx>>3)+(iw>>3)*6],
					temp[(kx>>3)+(iw>>3)*7],
				};

				dct8_inv_i8(x);
				
				image[ idx   <<2|kc]=x[0];
				image[(idx+1)<<2|kc]=x[1];
				image[(idx+2)<<2|kc]=x[2];
				image[(idx+3)<<2|kc]=x[3];
				image[(idx+4)<<2|kc]=x[4];
				image[(idx+5)<<2|kc]=x[5];
				image[(idx+6)<<2|kc]=x[6];
				image[(idx+7)<<2|kc]=x[7];
			}
		}
#endif
	}
	free(temp);
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
				p*=0xFF00;
				++p;
				p/=0x10000;
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