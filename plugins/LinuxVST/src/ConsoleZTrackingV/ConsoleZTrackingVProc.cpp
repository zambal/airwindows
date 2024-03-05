/* ========================================
 *  ConsoleZTrackingV - ConsoleZTrackingV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTrackingV_H
#include "ConsoleZTrackingV.h"
#endif

#define BL 65536
#define BM 65535

void ConsoleZTrackingV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames)
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	bool enableHPF = B != 0.0;
	bool enableLPF = D != 1.0;
	bool enableComp = F != 0.0;
	bool enableChorus = J != 0.0;
	bool enableReverb = L != 0.0;
	bool enableuLaw = Q != 0.0;
	bool enableCreature = V != 0.5;
	bool enableToneSlant = X != 0.5;

	double sampleRate = getSampleRate() * 2.0;
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= sampleRate;

	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;
	//this is going to be 2 for 88.1 or 96k, 3 for silly people, 4 for 176 or 192k

	double KK;
	double norm;

	// BiquadOneHalf HPF

	double cutoff_hpf = ((B*B*B*0.9999)+0.0001)*0.499;
	if (cutoff_hpf < 0.0001) cutoff_hpf = 0.0001;

  double res_hpf = (C*C*C*29.99)+0.01;
	if (res_hpf < 0.0001) res_hpf = 0.0001;

	KK = tan(M_PI * cutoff_hpf);
	norm = 1.0 / (1.0 + KK / res_hpf + KK * KK);
	double b0_hpf = norm;
	double b1_hpf = -2.0 * b0_hpf;
	double b2_hpf = b0_hpf;
	double b3_hpf = 2.0 * (KK * KK - 1.0) * norm;
	double b4_hpf = (1.0 - KK / res_hpf + KK * KK) * norm;

	// BiquadOneHalf LPF

	double cutoff_lpf = ((D*D*D*0.9999)+0.0001)*0.499;
	if (cutoff_lpf < 0.0001) cutoff_lpf = 0.0001;

  double res_lpf = (E*E*E*29.99)+0.01;
	if (res_lpf < 0.0001) res_lpf = 0.0001;

	KK = tan(M_PI * cutoff_lpf);
	norm = 1.0 / (1.0 + KK / res_lpf + KK * KK);
	double b0_lpf = KK * KK * norm;
	double b1_lpf = 2.0 * b0_lpf;
	double b2_lpf = b0_lpf;
	double b3_lpf = 2.0 * (KK * KK - 1.0) * norm;
	double b4_lpf = (1.0 - KK / res_lpf + KK * KK) * norm;

	// Buttercomp2
	
	double bc_inputgain = pow(10.0,(F*14.0)/20.0);
	double bc_compfactor = 0.012 * (F / 135.0);
	double bc_output = G * 2.0;
	double bc_wet = H;
	//removed extra dry variable
	double bc_outputgain = bc_inputgain;
	bc_outputgain -= 1.0;
	bc_outputgain /= 1.5;
	bc_outputgain += 1.0;

	// StereoChorus

	if (chorus_cycle > cycleEnd-1) chorus_cycle = cycleEnd-1; //sanity check
	
	double speed = pow(0.32+(I/6),10);
	double depth = (J/60) / speed;
	double tupi = 3.141592653589793238 * 2.0;

	// ClearCoat

	double cc_subRate = 0.001 / overallscale;
	double cc_wet = L*2.0;
	double cc_dry = 2.0 - cc_wet;
	if (cc_wet > 1.0) cc_wet = 1.0;
	if (cc_wet < 0.0) cc_wet = 0.0;
	if (cc_dry > 1.0) cc_dry = 1.0;
	if (cc_dry < 0.0) cc_dry = 0.0;
	//this reverb makes 50% full dry AND full wet, not crossfaded.
	//that's so it can be on submixes without cutting back dry channel when adjusted:
	//unless you go super heavy, you are only adjusting the added verb loudness.

	// Channel9

	double ch9_localiirAmount = ch9_iirAmount / overallscale;
	double localthreshold = ch9_threshold; //we've learned not to try and adjust threshold for sample rate
	double density = M*2.0; //0-2
	if (density > 1.0) density = 1.0; //max out at full wet for Spiral aspect
	double nonLin = 5.0-density; //number is smaller for more intense, larger for more subtle
	double ch9_c = ch9_cutoff / sampleRate;
  double ch9_resA = 1.618033988749894848204586;
	double ch9_resB = 0.618033988749894848204586;
	
	KK = tan(M_PI * ch9_c);
	norm = 1.0 / (1.0 + KK / ch9_resA + KK * KK);
	double ch9_b0A = KK * KK * norm;
	double ch9_b1A = 2.0 * ch9_b0A;
	double ch9_b2A = ch9_b0A;
	double ch9_b3A = 2.0 * (KK * KK - 1.0) * norm;
	double ch9_b4A = (1.0 - KK / ch9_resA + KK * KK) * norm;
		
	KK = tan(M_PI * ch9_c);
	norm = 1.0 / (1.0 + KK / ch9_resB + KK * KK);
	double ch9_b0B = KK * KK * norm;
	double ch9_b1B = 2.0 * ch9_b0B;
	double ch9_b2B = ch9_b0B;
	double ch9_b3B = 2.0 * (KK * KK - 1.0) * norm;
	double ch9_b4B = (1.0 - KK / ch9_resB + KK * KK) * norm;

	// Focus

	double boost = pow(10.0,(M*12.0)/20.0);
	double fs_cutoff = (pow(M_E, freqX * O) * centerFreq)/sampleRate; //fixed frequency, 3.515775k
	double fs_res = pow(pow(N,3)*2,2)+0.0001; //resonance
	int mode = (int) ( P * 4.999 );
	double output = R;
	double wet = S;
	
	KK = tan(M_PI * fs_cutoff);
	norm = 1.0 / (1.0 + KK / fs_res + KK * KK);
	double f1 = KK / fs_res * norm;
	double f2 = -f1;
	double f3 = 2.0 * (KK * KK - 1.0) * norm;
	double f4 = (1.0 - KK / fs_res + KK * KK) * norm;

	// Creature

	double source = 1.0-pow(1.0-T,5);
	int stages = (pow(U,2)*32.0*sqrt(overallscale))+1;
	double cr_wet = (V*2.0)-1.0; //inv-dry-wet for highpass
	double cr_dry = 2.0-(V*2.0);
	if (cr_dry > 1.0) cr_dry = 1.0; //full dry for use with inv, to 0.0 at full wet

	// ToneSlant

	double ts_overallscale = (W*99.0)+1.0;
	double ts_applySlant = (X*2.0)-1.0;


	//count to f(gain) which will be 0. f(0) is x1
	for (int count = 0; count <= ts_overallscale; count+=4) {
		Vec4d os(ts_overallscale);
		((Vec4d(1.0) - (Vec4d(count, count + 1, count + 2, count + 3) / os)) / os).store_a(ts_f + count);
		//recalc the filter and don't change the buffer it'll apply to
	}

	// +/- 12dB gain trim
	double gain_trim = pow(M_E, freqX * Y) * 0.25;

	Vec4ui fpd; fpd.load(fpd_b);
	Vec4d t, drySample;

  while (--sampleFrames >= 0)
  {
		Vec4d inputSample(*in1, *in2, 0.0, 0.0);
		inputSample *= (input_gain * 2.0);
		inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		
		// Channel9
		{
			Vec2d biquad1A; biquad1A.load_a(ch9_biquad);
			Vec2d biquad2A; biquad2A.load_a(ch9_biquad + 2);	
			Vec2d biquad3A; biquad3A.load_a(ch9_biquad + 4);
			Vec2d biquad4A; biquad4A.load_a(ch9_biquad + 6);	

			Vec2d tempSample = ch9_b0A*inputSample.get_low()+ch9_b1A*biquad1A+ch9_b2A*biquad2A-ch9_b3A*biquad3A-ch9_b4A*biquad4A;
  		biquad2A = biquad1A; biquad1A = inputSample.get_low(); 
			Vec2d inputSampleLow = select(abs(tempSample)<1.18e-37, 0.0, tempSample);
			biquad4A = biquad3A; biquad3A = inputSampleLow; //DF1 

			tempSample = ch9_b0A*inputSample.get_high()+ch9_b1A*biquad1A+ch9_b2A*biquad2A-ch9_b3A*biquad3A-ch9_b4A*biquad4A;
  		biquad2A = biquad1A; biquad1A = inputSample.get_high(); 
			inputSample = concatenate2(inputSampleLow, select(abs(tempSample)<1.18e-37, 0.0, tempSample));
			biquad4A = biquad3A; biquad3A = inputSample.get_high(); //DF1 

			biquad1A.store_a(ch9_biquad);
			biquad2A.store_a(ch9_biquad + 2);	
			biquad3A.store_a(ch9_biquad + 4);
			biquad4A.store_a(ch9_biquad + 6);	
		}

	 inputSample += ((pow_const(inputSample,5u)/128.0) + (pow_const(inputSample,9u)/262144.0)) - ((pow_const(inputSample,3u)/8.0) + (pow_const(inputSample,7u)/4096.0));

		 // {
		 // 	t = abs(2.0-((inputSample+nonLin)/nonLin)); // dielectricScale
		 // 	Vec4d iir; iir.load_a(ch9_iirSamples);
		
		 // 	iir = select(abs(iir)<1.18e-37, 0.0, iir);
		 // 	iir = (iir * (1.0 - (ch9_localiirAmount * t))) + (inputSample * ch9_localiirAmount * t);
		 // 	inputSample -= iir;

		 // 	iir.store_a(ch9_iirSamples);
		 // }

		 if(enableHPF) {
		 	Vec4d biquadA; biquadA.load_a(biquad_hpf);
		 	Vec4d biquadB; biquadB.load_a(biquad_hpf + 4);

		 	t = inputSample * b0_hpf + biquadA;
		 	biquadA = (inputSample * b1_hpf) - (t * b3_hpf) + biquadB;
		 	biquadB = (inputSample * b2_hpf) - (t * b4_hpf);
		 	inputSample = min(max(t, -1.0), 1.0);

		 	biquadA.store_a(biquad_hpf);
		 	biquadB.store_a(biquad_hpf + 4);
		 }

		 if(enableLPF) {
		 	Vec4d biquadA; biquadA.load_a(biquad_lpf);
		 	Vec4d biquadB; biquadB.load_a(biquad_lpf + 4);

		 	t = inputSample * b0_lpf + biquadA;
		 	biquadA = (inputSample * b1_lpf) - (t * b3_lpf) + biquadB;
		 	biquadB = (inputSample * b2_lpf) - (t * b4_lpf);
		 	inputSample = min(max(t, -1.0), 1.0);

		 	biquadA.store_a(biquad_lpf);
		 	biquadB.store_a(biquad_lpf + 4);
		 }

		 inputSample += (pow_const(inputSample,3u)/4.0)+(pow_const(inputSample,5u)/8.0)+(pow_const(inputSample,7u)/16.0)+(pow_const(inputSample,9u)/32.0);
		// amplitude aspect

		 // Focus

		 drySample = inputSample;
		
		 inputSample = sin(inputSample);
		 //encode Console5: good cleanness
		 {
		 	Vec2d p1, p2; p1.load_a(figure); p2.load_a(figure + 2);
		 	Vec2d tempSample = (inputSample.get_low() * f1) + p1;
		 	p1 = -(tempSample * f3) + p2;
		 	p2 = (inputSample.get_low() * f2) - (tempSample * f4);
		 	Vec2d inputSampleLow = tempSample;
		
		 	tempSample = (inputSample.get_high() * f1) + p1;
		 	p1 = -(tempSample * f3) + p2;
		 	p2 = (inputSample.get_high() * f2) - (tempSample * f4);
		 	p1.store_a(figure); p2.store_a(figure + 2);
		 	inputSample = concatenate2(inputSampleLow, tempSample);
		 }		
		 inputSample = min(max(inputSample, -1.0), 1.0);
		 //without this, you can get a NaN condition where it spits out DC offset at full blast!

		 inputSample = Vec4d(asin(inputSample.extract(0)), asin(inputSample.extract(1)), asin(inputSample.extract(2)), asin(inputSample.extract(3)));
		 //decode Console5
		
		 Vec4d groundSample = drySample - inputSample; //set up UnBox

		 // uLawEncode

		 if(enableuLaw) {
		 	inputSample = min(max(inputSample, -1.0), 1.0);

		 	t = log(1.0+(255.0*abs(inputSample)));
		 	inputSample = if_mul(inputSample < 0, t, -1.0);
		 	inputSample /= log(Vec4d(256.0));
		 }

		 // StereoChorus

		 if(enableChorus) {
		 	chorus_cycle++;
		  if (chorus_cycle == cycleEnd) { //hit the end point and we do a chorus sample
		  	//assign working variables
		  	Vec2d airPrev; airPrev.load_a(airPrev_b);
		 	inputSample.get_high().store_a(airPrev_b);
		  	Vec2d airFactor = airPrev - inputSample.get_low();
		  	Vec2d airEven; airEven.load_a(airEven_b);
		  	Vec2d airOdd; airOdd.load_a(airOdd_b);
		 	airEven += airFactor; airOdd -= airFactor; airFactor = airEven;
		  	airOdd = ((airOdd - ((airOdd - airEven)*0.00390625)) * 0.9999);
		  	airEven = ((airEven - ((airEven - airOdd)*0.00390625)) * 0.9999);
		  	Vec2d tt = airFactor;
		  	airFactor = inputSample.get_low() - inputSample.get_high();
		 	airEven -= airFactor; airOdd += airFactor; airFactor = airOdd;
		  	airOdd = ((airOdd - ((airOdd - airEven)*0.00390625)) * 0.9999);
		  	airEven = ((airEven - ((airEven - airOdd)*0.00390625)) * 0.9999);
		 	airEven.store(airEven_b); airOdd.store(airOdd_b);
		 	inputSample += concatenate2(tt, airFactor);
		 	//air, compensates for loss of highs in flanger's interpolation
			
		   truncate_to_int32(inputSample*8388352.0).store(p_b + gcount);
		   Vec4i count(gcount, gcount + 1, gcount + 2, gcount + 3);
		   Vec2d prev_sweep; prev_sweep.load_a(sweep_b);
		   Vec4d sweep(prev_sweep, if_sub(prev_sweep + speed > tupi, prev_sweep + speed, tupi));
		   Vec4d offset = depth + (depth * sin(sweep));
		   count -= (truncate_to_int32(floor(offset)) * 2);
		   int count_b[4] __attribute__((aligned(32))); count.store_a(count_b);
		   Vec4i temp1(p_b[count_b[0] & BM], p_b[count_b[1] & BM], p_b[count_b[2] & BM], p_b[count_b[3] & BM]);
		   Vec4i temp2(p_b[(count_b[0] - 2) & BM], p_b[(count_b[1] - 2) & BM], p_b[(count_b[2] - 2) & BM], p_b[(count_b[3] - 2) & BM]);
		   Vec4i temp3(p_b[(count_b[0] - 4) & BM], p_b[(count_b[1] - 4) & BM], p_b[(count_b[2] - 4) & BM], p_b[(count_b[3] - 4) & BM]);
		 
		   temp1 = truncate_to_int32(to_double(temp1) * (1.0 - (offset - floor(offset))));
		   temp1 += temp2;
		   temp1 += truncate_to_int32(to_double(temp3) * (offset - floor(offset)));
		   temp1 -= ((temp1 - temp2)-(temp2-temp3)) / 50;
		 
		   sweep += speed;
		   sweep = if_sub(sweep > tupi, sweep, tupi);
		   sweep.get_high().store_a(sweep_b);
		 
		   gcount = (gcount + 4) & BM;
		   //still scrolling through the samples, remember

		   inputSample = to_double(temp1)/16776704.0;
		 
		  	if (cycleEnd == 4) {
		  		Vec2d l1; l1.load_a(lastRef + 16);
		  		Vec2d l2 = inputSample.get_low();
		  		Vec2d l3 = inputSample.get_high();

		  		l1.store_a(lastRef);
		  		((l1 + l1 + l1 + l2) * 0.25).store_a(lastRef + 2);
		  		((l1 + l2) / 2.0).store_a(lastRef + 4);
		  		((l1 + l2 + l2 + l2) * 0.25).store_a(lastRef + 6);
		  		l2.store_a(lastRef + 8);
		  		((l2 + l2 + l2 + l3) * 0.25).store_a(lastRef + 10);
		  		((l2 + l3) / 2.0).store_a(lastRef + 12);
		  		((l2 + l3 + l3 + l3) * 0.25).store_a(lastRef + 14);
		  		l3.store_a(lastRef + 16);
		  	}
		  	if (cycleEnd == 2) {
		  		Vec2d l1; l1.load_a(lastRef + 8);
		  		Vec2d l2 = inputSample.get_low();
		  		Vec2d l3 = inputSample.get_high();

		  		l1.store_a(lastRef);
		  		((l1 + l2) / 2.0).store_a(lastRef + 2);
		  		l2.store_a(lastRef + 4);
		  		((l2 + l3) / 2.0).store_a(lastRef + 6);
		  		l3.store_a(lastRef + 8);
		  	}
		  	if (cycleEnd == 1) {
		  		inputSample.store_a(lastRef);
		  	}
		  	chorus_cycle = 0; //reset
		  	inputSample.load_a(lastRef);
		  } else {
		  	inputSample.load_a(lastRef + (chorus_cycle * 4));
		  	//we are going through our references now
		  }
		 }
		 // ClearCoat
		
		 if(enableReverb) {
		 	Vec4d cc_drySample = inputSample;

		 	cc_cycle+=2;
		 	if (cc_cycle >= cycleEnd) { //hit the end point and we do a reverb sample
		 		Vec2d o1, o2, o3, o4, vv;
						
		 		vv = inputSample.get_low();
			
		 		o1.load_a(feedbackA); 
		 		o2.load_a(feedbackB); 
		 		o3.load_a(feedbackC); 
		 		o4.load_a(feedbackD); 

		 		o1 = (vv + (o1 * 0.04166666666)); 
		 		o2 = (vv + (o2 * 0.04166666666)); 
		 		o3 = (vv + (o3 * 0.04166666666)); 
		 		o4 = (vv + (o4 * 0.04166666666)); 

		 		o1.store_a(aA + (cc_gcount & 2047));
		 		o2.store_a(aB + (cc_gcount & 4095));
		 		o3.store_a(aC + (cc_gcount & 4095));
		 		o4.store_a(aD + (cc_gcount & 8191));

		 		o1 = Vec2d(aA[(cc_gcount - shortA) & 2047], aA[(cc_gcount + 1 - shortD) & 2047]);
		 		o2 = Vec2d(aB[(cc_gcount - shortB) & 4095], aB[(cc_gcount + 1 - shortH) & 4095]);
		 		o3 = Vec2d(aC[(cc_gcount - shortC) & 4095], aC[(cc_gcount + 1 - shortL) & 4095]);
		 		o4 = Vec2d(aD[(cc_gcount - shortD) & 8191], aD[(cc_gcount + 1 - shortP) & 8191]);

		 		(o1 - (o2 + o3 + o4)).store_a(aE + (cc_gcount & 4095));
		 		(o2 - (o1 + o3 + o4)).store_a(aF + (cc_gcount & 4095));
		 		(o3 - (o1 + o2 + o4)).store_a(aG + (cc_gcount & 4095));
		 		(o4 - (o1 + o2 + o3)).store_a(aH + (cc_gcount & 4095));

		 		o1 = Vec2d(aE[(cc_gcount - shortE) & 4095], aE[(cc_gcount + 1 - shortC) & 4095]);
		 		o2 = Vec2d(aF[(cc_gcount - shortF) & 4095], aF[(cc_gcount + 1 - shortG) & 4095]);
		 		o3 = Vec2d(aG[(cc_gcount - shortG) & 4095], aG[(cc_gcount + 1 - shortK) & 4095]);
		 		o4 = Vec2d(aH[(cc_gcount - shortH) & 4095], aH[(cc_gcount + 1 - shortO) & 4095]);

		 		(o1 - (o2 + o3 + o4)).store_a(aI + (cc_gcount & 4095));
		 		(o2 - (o1 + o3 + o4)).store_a(aJ + (cc_gcount & 8191));
		 		(o3 - (o1 + o2 + o4)).store_a(aK + (cc_gcount & 8191));
		 		(o4 - (o1 + o2 + o3)).store_a(aL + (cc_gcount & 8191));
	
		 		o1 = Vec2d(aI[(cc_gcount - shortI) & 4095], aI[(cc_gcount + 1 - shortB) & 4095]);
		 		o2 = Vec2d(aJ[(cc_gcount - shortJ) & 8191], aJ[(cc_gcount + 1 - shortF) & 8191]);
		 		o3 = Vec2d(aK[(cc_gcount - shortK) & 8191], aK[(cc_gcount + 1 - shortJ) & 8191]);
		 		o4 = Vec2d(aL[(cc_gcount - shortL) & 8191], aL[(cc_gcount + 1 - shortN) & 8191]);
	
		 		(o1 - (o2 + o3 + o4)).store_a(aM + (cc_gcount & 2047));
		 		(o2 - (o1 + o3 + o4)).store_a(aN + (cc_gcount & 8191));
		 		(o3 - (o1 + o2 + o4)).store_a(aO + (cc_gcount & 2047));
		 		(o4 - (o1 + o2 + o3)).store_a(aP + (cc_gcount & 8191));
	
		 		o1 = Vec2d(aM[(cc_gcount - shortM) & 2047], aM[(cc_gcount + 1 - shortA) & 2047]);
		 		o2 = Vec2d(aN[(cc_gcount - shortN) & 8191], aN[(cc_gcount + 1 - shortE) & 8191]);
		 		o3 = Vec2d(aO[(cc_gcount - shortO) & 2047], aO[(cc_gcount + 1 - shortI) & 2047]);
		 		o4 = Vec2d(aP[(cc_gcount - shortP) & 8191], aP[(cc_gcount + 1 - shortM) & 8191]);
	
		 		vv.load_a(prevMulchA);
		 		o1.store_a(prevMulchA);
		 		o1 = (o1 + o1 + o1 + vv)*0.25;

		 		(o1 - (o2 + o3 + o4)).store_a(feedbackA);
		 		(o2 - (o1 + o3 + o4)).store_a(feedbackB);
		 		(o3 - (o1 + o2 + o4)).store_a(feedbackC);
		 		(o4 - (o1 + o2 + o3)).store_a(feedbackD);

		 		vv = min(max((o1 + o2 + o3 + o4)*0.125, -1.0), 1.0);
		 		//and take the final combined sum of outputs, corrected for Householder gain

		 		cc_gcount = (cc_gcount + 2) & 8191;

		 		if (cycleEnd == 4) {
		 			Vec2d tt; tt.load_a(cc_lastRef + 4);
		 			((tt + vv) * 0.5).store_a(cc_lastRef);
		 			((tt + vv + vv + vv) * 0.25).store_a(cc_lastRef + 2);
		 			vv.store_a(cc_lastRef + 4);
		 			inputSample = concatenate2(tt, (tt + tt + tt + vv) * 0.25);
		 		}
		 		if (cycleEnd <= 2) {
		 			Vec2d tt; tt.load_a(cc_lastRef);
		 			vv.store_a(cc_lastRef);
		 			inputSample = concatenate2(tt, (tt + vv) * 0.5);
		 		}
		 		cc_cycle = 0; //reset
			
		 	} else {
		 		inputSample.load_a(cc_lastRef);
		 		//we are going through our references now
		 	}
	
		 	//begin SubTight section
		 	Vec4d subSample = inputSample * cc_subRate;
		 	Vec4d subA; subA.load_a(cc_sub_b);
		 	Vec4d subB; subB.load_a(cc_sub_b + 4);
		 	Vec4d subC; subC.load_a(cc_sub_b + 8);
		 	Vec4d subD; subD.load_a(cc_sub_b + 12);
		 	Vec4d scale = 0.5+abs(subSample*0.5);
		 	subSample = (subA+(sin(subA-subSample)*scale));
		 	subA = subSample*scale;
		 	scale = 0.5+abs(subSample*0.5);
		 	subSample = (subB+(sin(subB-subSample)*scale));
		 	subB = subSample*scale;
		 	scale = 0.5+abs(subSample*0.5);
		 	subSample = (subC+(sin(subC-subSample)*scale));
		 	subC = subSample*scale;
		 	scale = 0.5+abs(subSample*0.5);
		 	subSample = (subD+(sin(subD-subSample)*scale));
		 	subD = subSample*scale;
		 	subA.store_a(cc_sub_b);
		 	subB.store_a(cc_sub_b + 4);
		 	subC.store_a(cc_sub_b + 8);
		 	subD.store_a(cc_sub_b + 12);

		 	subSample = min(max(subSample, -0.25), 0.25);
		 	inputSample -= (subSample*16.0);
		 	//end SubTight section		

		 	if (cycleEnd > 1) {
		 		Vec2d tt; tt.load_a(tail); inputSample.get_high().store_a(tail);
		 		inputSample = (inputSample + concatenate2(tt, inputSample.get_low())) * 0.5;
		 	} //let's average only at elevated sample rates
	
		 	if (cc_wet < 1.0) {inputSample *= cc_wet;}
		 	if (cc_dry < 1.0) {cc_drySample *= cc_dry;}
		 	inputSample += cc_drySample;
		 	//this is our submix verb dry/wet: 0.5 is BOTH at FULL VOLUME
		 	//purpose is that, if you're adding verb, you're not altering other balances

		 }

		 inputSample *= boost; //now, focussed area gets cranked before distort
		
		 switch (mode)
		 {
		 	case 0: //Density
		 		inputSample = min(max(inputSample, -1.570796326794897), 1.570796326794897);
		 		//clip to 1.570796326794897 to reach maximum output
		 		inputSample = sin(inputSample);
		 		break;
		 	case 1: //Drive				
		 		inputSample -= (inputSample * (abs(inputSample) * 0.6) * (abs(inputSample) * 0.6));
		 		inputSample *= 1.6;
		 		break;
		 	case 2: //Spiral
		 		inputSample = min(max(inputSample, -1.2533141373155), 1.2533141373155);
		 		//clip to 1.2533141373155 to reach maximum output
		 		inputSample = sin(inputSample * abs(inputSample)) / select(abs(inputSample) == 0.0, Vec4d(1.0), abs(inputSample));
		 		break;
		 	case 3: //Mojo
		 		Vec4d mojo; mojo = pow(abs(inputSample),Vec4d(0.25));
		 		inputSample = select(mojo > 0.0, (sin(inputSample * mojo * M_PI * 0.5) / mojo) * 0.987654321, inputSample);
		 		//mojo is the one that flattens WAAAAY out very softly before wavefolding				
		 		break;
		 	case 4: //Dyno
		 		Vec4d dyno; dyno = pow_const(abs(inputSample),4);
		 		inputSample = select(dyno > 0.0, (sin(inputSample * dyno) / dyno) * 1.1654321, inputSample);
		 		//dyno is the one that tries to raise peak energy				
		 		break;
		 }				
		
		 if (output != 1.0) {
		 	inputSample *= output;
		 }
		
		 // uLawDecode

		 if(enableuLaw) {
		 	inputSample = min(max(inputSample, -1.0), 1.0);
		
		 	t = pow(Vec4d(256.0),abs(inputSample))-1.0;
		 	inputSample = if_mul(inputSample < 0, t, -1.0);
		 	inputSample /= 255.0;
		 }

		 inputSample += groundSample; //effectively UnBox
		
		 if (wet !=1.0) {
		 	inputSample = (inputSample * wet) + (drySample * (1.0-wet));
		 }


		 // Creature

		 if(enableCreature) {
		 	Vec4d cr_drySample = inputSample;
		
		 	Vec2d inputSample1 = inputSample.get_low(), inputSample2 = inputSample.get_high();
		 	Vec2d slew;
		 	for (int x = 0; x < stages; x++) {
		 		slew.load_a(cr_slew + (x * 2));
		 		inputSample1 = (slew+(sin(slew-inputSample1)*0.5))*source;
		 		slew = inputSample1 * 0.5;
		 		inputSample2 = (slew+(sin(slew-inputSample2)*0.5))*source;
		 		(0.5 * inputSample2).store_a(cr_slew + (x * 2));
		 	}

		 	inputSample = concatenate2(inputSample1, inputSample2);

		 	if (stages % 2 > 0) {
		 		inputSample = -inputSample;
		 	}
		
		 	inputSample = (cr_drySample * cr_dry) + (inputSample * cr_wet);
		 }

		 // ToneSlant
		 if(enableToneSlant) {
		 	inputSample.store_a(ts_b + ts_gcount);
		 	if(ts_gcount == 0) inputSample.store_a(ts_b + 256); 

		 	Vec4d acc = inputSample * ts_f[0];
		 	for (int count = 1; count < ts_overallscale; count++) {
		 		t.load(ts_b + ((ts_gcount - (count * 2)) & 255));
		 		acc += t * ts_f[count];
		 	}
		 	ts_gcount = (ts_gcount + 4) & 255;

		 	//we're gonna apply the total effect of all these calculations as a single subtract
		 	inputSample += ((inputSample - (acc*2.0)) * ts_applySlant);
		 	//our one math operation on the input data coming in
		 }

		 // ButterComp2
		
		 if(enableComp) {
		 	Vec4d controlpos; controlpos.load(controlpos_b);
		 	Vec4d controlneg; controlneg.load(controlneg_b);
		 	Vec4d targetpos; targetpos.load(targetpos_b);
		 	Vec4d targetneg; targetneg.load(targetneg_b);
		 	Vec4d lastOutput; lastOutput.load(lastOutput_b);

		 	Vec4d divisor = Vec4d(bc_compfactor) / (abs(lastOutput)+1.0);
		 	//this is slowing compressor recovery while output waveforms were high
		 	divisor /= overallscale;
		 	Vec4d remainder = divisor;
		 	divisor = 1.0 - divisor;
		 	//recalculate divisor every sample		

		 	Vec4d inputpos = max(inputSample + 1.0, 0);
		 	Vec4d outputpos = min(inputpos / 2.0, 1.0);
		 	inputpos *= inputpos;
		 	targetpos *= divisor;
		 	targetpos += (inputpos * remainder);
		 	Vec4d calcpos = pow_const((1.0/targetpos),2);

		 	Vec4d inputneg = max(-1.0 * inputSample + 1.0, 0);
		 	Vec4d outputneg = min(inputneg / 2.0, 1.0);
		 	inputneg *= inputneg;
		 	targetneg *= divisor;
		 	targetneg += (inputneg * remainder);
		 	Vec4d calcneg = pow_const((1.0/targetneg),2);
		 	//now we have mirrored targets for comp
		 	//outputpos and outputneg go from 0 to 1

		 	controlpos = select(inputSample > 0.0, (controlpos * divisor) + (calcpos * remainder), controlpos);
		 	controlneg = select(inputSample < 0.0, (controlneg * divisor) + (calcneg * remainder), controlneg);
		 	//this causes each of the four to update only when active and in the correct 'flip'

		 	t = (controlpos * outputpos) + (controlneg * outputneg);
		 	//this combines the sides according to flip, blending relative to the input value
		
		 	drySample = inputSample;
		 	inputSample *= bc_inputgain;
		 	inputSample /= bc_outputgain;
		 	inputSample *= t;
		
	
		 	if (bc_output != 1.0) {
		 		inputSample *= bc_output;
		 	}

		 	if (bc_wet !=1.0) {
		 		inputSample = (inputSample * bc_wet) + (drySample * (1.0-bc_wet));
		 	}
		
		 	lastOutput = inputSample;
		 	//we will make this factor respond to use of dry/wet

		 	controlpos.store(controlpos_b);
		 	controlneg.store(controlneg_b);
		 	targetpos.store(targetpos_b);
		 	targetneg.store(targetneg_b);
		 	lastOutput.store(lastOutput_b);
		 }

		 // Gain trim

		 inputSample = inputSample * gain_trim;

		 // ADClip

		 {
		 	double softness = 0.618033988749894848204586;
		 	Vec2d last; last.load_a(adcLastSample);
		 	Vec2d input1 = inputSample.get_low();
		 	Vec2d input2 = inputSample.get_high();

		 	last = select(last >= 0.5, select(input1 < 0.5, (0.5*softness) + (input1 * (1.0-softness)), 0.5), last);
		 	last = select(last <= -0.5, select(input1 > -0.5, (-0.5*softness) + (input1 * (1.0-softness)), -0.5), last);
		 	input1 = select(input1 > 0.5, select(last < 0.5, (0.5*softness) + (last * (1.0-softness)), 0.5), input1);
		 	input1 = select(input1 < -0.5, select(last > -0.5, (-0.5*softness) + (last * (1.0-softness)), -0.5), input1);
		  last = input1;
		 	last = select(last >= 0.5, select(input2 < 0.5, (0.5*softness) + (input2 * (1.0-softness)), 0.5), last);
		 	last = select(last <= -0.5, select(input2 > -0.5, (-0.5*softness) + (input2 * (1.0-softness)), -0.5), last);
		 	input2 = select(input2 > 0.5, select(last < 0.5, (0.5*softness) + (last * (1.0-softness)), 0.5), input2);
		 	input2 = select(input2 < -0.5, select(last > -0.5, (-0.5*softness) + (last * (1.0-softness)), -0.5), input2);

		 	last.store_a(adcLastSample);			

		 	inputSample = min(max(concatenate2(input1, input2), -0.5), 0.5);
		 	//final iron bar
		 }

		 // Channel9

		 {
		 	Vec2d biquad1B; biquad1B.load_a(ch9_biquad + 8);	
		 	Vec2d biquad2B; biquad2B.load_a(ch9_biquad + 10);	
		 	Vec2d biquad3B; biquad3B.load_a(ch9_biquad + 12);	
		 	Vec2d biquad4B; biquad4B.load_a(ch9_biquad + 14);	
		 	Vec2d lastSampleA; lastSampleA.load_a(ch9_lastSample);
		 	Vec2d lastSampleB; lastSampleB.load_a(ch9_lastSample + 2);
		 	Vec2d lastSampleC; lastSampleC.load_a(ch9_lastSample + 4);

		 	Vec2d clamp = (lastSampleB - lastSampleC) * 0.381966011250105;
		 	clamp -= (lastSampleA - lastSampleB) * 0.6180339887498948482045;
		 	clamp += inputSample.get_low() - lastSampleA; //regular slew clamping added
		
		 	lastSampleC = lastSampleB;
		 	lastSampleB = lastSampleA;
		 	lastSampleA = inputSample.get_low(); //now our output relates off lastSampleB
		
		 	Vec2d inputSampleLow = select(clamp > localthreshold, lastSampleB + localthreshold, inputSample.get_low());
		 	inputSampleLow = select(-clamp > localthreshold, lastSampleB - localthreshold, inputSampleLow);  
		
		 	lastSampleA = (lastSampleA*0.381966011250105)+(inputSampleLow*0.6180339887498948482045); //split the difference between raw and smoothed for buffer

		 	clamp = (lastSampleB - lastSampleC) * 0.381966011250105;
		 	clamp -= (lastSampleA - lastSampleB) * 0.6180339887498948482045;
		 	clamp += inputSample.get_high() - lastSampleA; //regular slew clamping added
		
		 	lastSampleC = lastSampleB;
		 	lastSampleB = lastSampleA;
		 	lastSampleA = inputSample.get_high(); //now our output relates off lastSampleB
		
		 	Vec2d inputSampleHigh = select(clamp > localthreshold, lastSampleB + localthreshold, inputSample.get_high());
		 	inputSampleHigh = select(-clamp > localthreshold, lastSampleB - localthreshold, inputSampleHigh);  
		
		 	lastSampleA = (lastSampleA*0.381966011250105)+(inputSampleHigh*0.6180339887498948482045); //split the difference between raw and smoothed for buffer
		 	inputSample = concatenate2(inputSampleLow, inputSampleHigh);
		
		 	Vec2d tempSample = ch9_b0B*inputSample.get_low()+ch9_b1B*biquad1B+ch9_b2B*biquad2B-ch9_b3B*biquad3B-ch9_b4B*biquad4B;
   		biquad2B = biquad1B; biquad1B = inputSample.get_low(); 
		 	inputSampleLow = select(abs(tempSample)<1.18e-37, 0.0, tempSample);
		 	biquad4B = biquad3B; biquad3B = inputSampleLow; //DF1 

		 	tempSample = ch9_b0B*inputSample.get_high()+ch9_b1B*biquad1B+ch9_b2B*biquad2B-ch9_b3B*biquad3B-ch9_b4B*biquad4B;
   		biquad2B = biquad1B; biquad1B = inputSample.get_high(); 
		 	inputSample = concatenate2(inputSampleLow, select(abs(tempSample)<1.18e-37, 0.0, tempSample));
		 	biquad4B = biquad3B; biquad3B = inputSample.get_high(); //DF1 

		 	biquad1B.store_a(ch9_biquad + 8);	
		 	biquad2B.store_a(ch9_biquad + 10);	
		 	biquad3B.store_a(ch9_biquad + 12);	
		 	biquad4B.store_a(ch9_biquad + 14);	
		 	lastSampleA.store_a(ch9_lastSample);
		 	lastSampleB.store_a(ch9_lastSample + 2);
		 	lastSampleC.store_a(ch9_lastSample + 4);
		 }		

		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (to_double(fpd) - 2147483647.0) * to_double(exp) * 5.5e-36l;
		//end 32 bit stereo floating point dither

		double result[2];
		inputSample.get_low().store_a(result);

		*out1 = result[0];
		*out2 = result[1];

		*in1++;
		*in2++;
		*out1++;
		*out2++;
	}
	fpd.store(fpd_b);
}

void ConsoleZTrackingV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
{
  double* in1  =  inputs[0];
  double* in2  =  inputs[1];
  double* out1 = outputs[0];
  double* out2 = outputs[1];

}
