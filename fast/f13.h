#ifdef F13_ENC2
	int
		eNW	=rows[1][-1],
		eN	=rows[1][+0],
		eNE	=rows[1][+1],
		eW	=rows[0][-1];

	//const int ctx=0;

	int ctx1=QUANTIZE(eN);
	int ctx2=QUANTIZE(eW);
	int ctx3=QUANTIZE(eNW);
	int ctx4=QUANTIZE(eNE);
	//ctx2+=nclevels;
	//ctx3+=nclevels;
	ctx2+=nclevels;
	ctx3+=nclevels<<1;
	ctx4+=nclevels*3;

	//int aN=abs(eN), aW=abs(eW);
	//int alpha=(aN<<8)/(aN+aW+1);

	//int ctx=QUANTIZE(abs(eN-eW));

	//int ctx=QUANTIZE(eN+eW);
	//ctx=CLAMP(0, ctx, nclevels);

	//int ctx=QUANTIZE(eN);
	//ctx=3*ctx+QUANTIZE(eW);

	//ctx=3*ctx+QUANTIZE(eNW);
	//ctx=3*ctx+QUANTIZE(eNE);
//#ifdef DEBUG_CHECKS
//	if((unsigned)ctx>=3*3*3*3)
//		LOG_ERROR("ctx %d", ctx);
//#endif
	ctx1*=nqlevels+1;
	ctx2*=nqlevels+1;
	ctx3*=nqlevels+1;
	ctx4*=nqlevels+1;
	unsigned *curr_CDF1=CDF +ctx1;
	int     *curr_hist1=hist+ctx1;
	unsigned *curr_CDF2=CDF +ctx2;
	int     *curr_hist2=hist+ctx2;
	unsigned *curr_CDF3=CDF +ctx3;
	int     *curr_hist3=hist+ctx3;
	unsigned *curr_CDF4=CDF +ctx4;
	int     *curr_hist4=hist+ctx4;
	PROF(CTX);

	int delta=src->data[idx]-pred;
	rows[0][0]=(short)delta;
	delta+=nlevels>>1;
	delta&=nlevels-1;
	delta-=nlevels>>1;
	delta=delta<<1^-(delta<0);

	int token, bypass, nbits;
	token=quantize_pixel(delta, &bypass, &nbits);
	PROF(TOKEN);

	ac_enc_av4(ec, token, curr_CDF1, curr_CDF2, curr_CDF3, curr_CDF4, alpha, beta);
	//ac_enc_av2(ec, token, curr_CDF1, curr_CDF2, alpha);
	if(nbits)
		ac_enc_bypass(ec, bypass, 1<<nbits);
	PROF(CODE);
#ifdef DEBUG_CHECKS
	if(token>=nqlevels)
		LOG_ERROR("enc2: token %d/%d", token, nqlevels);
#endif

	//int updateA=(THREEWAY(curr_hist1[token]+curr_hist3[token], curr_hist2[token]+curr_hist4[token])<<14)+0x4000;
	//int updateB=(THREEWAY(curr_hist1[token]+curr_hist2[token], curr_hist3[token]+curr_hist4[token])<<14)+0x4000;
	//alpha+=(updateA-alpha)>>10;
	//beta +=(updateB-alpha)>>10;

	//int update=(THREEWAY(curr_hist1[token], curr_hist2[token])<<14)+0x4000;
	//alpha+=(update-alpha)>>10;

	++curr_hist1[token];
	++curr_hist1[nqlevels];
	if(curr_hist1[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist1, curr_CDF1, nqlevels);
		restale_hist(curr_hist1, nqlevels);
	}
	++curr_hist2[token];
	++curr_hist2[nqlevels];
	if(curr_hist2[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist2, curr_CDF2, nqlevels);
		restale_hist(curr_hist2, nqlevels);
	}
	++curr_hist3[token];
	++curr_hist3[nqlevels];
	if(curr_hist3[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist3, curr_CDF3, nqlevels);
		restale_hist(curr_hist3, nqlevels);
	}
	++curr_hist4[token];
	++curr_hist4[nqlevels];
	if(curr_hist4[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist4, curr_CDF4, nqlevels);
		restale_hist(curr_hist4, nqlevels);
	}
	PROF(STATS);
#endif
#ifdef F13_DEC2
#ifdef DEBUG_CHECKS
	static ctr=0;//
	++ctr;
#endif
	int
		eNW	=rows[1][-1],
		eN	=rows[1][+0],
		eNE	=rows[1][+1],
		eW	=rows[0][-1];
	
	//const int ctx=0;
	
	int ctx1=QUANTIZE(eN);
	int ctx2=QUANTIZE(eW);
	int ctx3=QUANTIZE(eNW);
	int ctx4=QUANTIZE(eNE);
	//ctx2+=nclevels;
	//ctx3+=nclevels;
	ctx2+=nclevels;
	ctx3+=nclevels<<1;
	ctx4+=nclevels*3;

	//int aN=abs(eN), aW=abs(eW);
	//int alpha=(aN<<8)/(aN+aW+1);

	//int ctx=QUANTIZE(abs(eN-eW));

	//int ctx=QUANTIZE(eN+eW);
	//ctx=CLAMP(0, ctx, nclevels);

	//int ctx=QUANTIZE(eN);
	//ctx=3*ctx+QUANTIZE(eW);

	//ctx=3*ctx+QUANTIZE(eNW);
	//ctx=3*ctx+QUANTIZE(eNE);
//#ifdef DEBUG_CHECKS
//	//if(ctr==15451)
//	//	printf("ctx %d  nb %d %d %d %d\n", ctx, eN, eW, eNW, eNE);
//#endif
	ctx1*=nqlevels+1;
	ctx2*=nqlevels+1;
	ctx3*=nqlevels+1;
	ctx4*=nqlevels+1;
	unsigned *curr_CDF1=CDF +ctx1;
	int     *curr_hist1=hist+ctx1;
	unsigned *curr_CDF2=CDF +ctx2;
	int     *curr_hist2=hist+ctx2;
	unsigned *curr_CDF3=CDF +ctx3;
	int     *curr_hist3=hist+ctx3;
	unsigned *curr_CDF4=CDF +ctx4;
	int     *curr_hist4=hist+ctx4;
	PROF(CTX);

	int token=ac_dec_packedsign_av4(ec, curr_CDF1, curr_CDF2, curr_CDF3, curr_CDF4, alpha, beta, nqlevels);
	//int token=ac_dec_packedsign_av2(ec, curr_CDF1, curr_CDF2, alpha, nqlevels);
	int bypass, nbits;
	PROF(CODE);
	int delta=token;
	if(delta>=(1<<CONFIG_EXP))
	{
		delta-=1<<CONFIG_EXP;
		int lsb=delta&((1<<CONFIG_LSB)-1);
		delta>>=CONFIG_LSB;
		int msb=delta&((1<<CONFIG_MSB)-1);
		delta>>=CONFIG_MSB;
		nbits=delta+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
		bypass=ac_dec_bypass(ec, 1<<nbits);
		delta=1;
		delta<<=CONFIG_MSB;
		delta|=msb;
		delta<<=nbits;
		delta|=bypass;
		delta<<=CONFIG_LSB;
		delta|=lsb;
	}
	delta=delta>>1^-(delta&1);
	delta+=pred;
	delta+=nlevels>>1;
	delta&=nlevels-1;
	delta-=nlevels>>1;
	rows[0][0]=(short)(delta-pred);
	PROF(TOKEN);
#ifdef DEBUG_CHECKS
	if(token>=nqlevels)
		LOG_ERROR("dec2: token %d/%d  ctr %d", token, nqlevels, ctr);
#endif
	
	//int updateA=(THREEWAY(curr_hist1[token]+curr_hist3[token], curr_hist2[token]+curr_hist4[token])<<14)+0x4000;
	//int updateB=(THREEWAY(curr_hist1[token]+curr_hist2[token], curr_hist3[token]+curr_hist4[token])<<14)+0x4000;
	//alpha+=(updateA-alpha)>>10;
	//beta +=(updateB-alpha)>>10;

	++curr_hist1[token];
	++curr_hist1[nqlevels];
	if(curr_hist1[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist1, curr_CDF1, nqlevels);
		restale_hist(curr_hist1, nqlevels);
	}
	++curr_hist2[token];
	++curr_hist2[nqlevels];
	if(curr_hist2[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist2, curr_CDF2, nqlevels);
		restale_hist(curr_hist2, nqlevels);
	}
	++curr_hist3[token];
	++curr_hist3[nqlevels];
	if(curr_hist3[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist3, curr_CDF3, nqlevels);
		restale_hist(curr_hist3, nqlevels);
	}
	++curr_hist4[token];
	++curr_hist4[nqlevels];
	if(curr_hist4[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist4, curr_CDF4, nqlevels);
		restale_hist(curr_hist4, nqlevels);
	}
	PROF(STATS);
#endif