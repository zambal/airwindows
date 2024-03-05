/* ========================================
 *  ConsoleZChannelV - ConsoleZChannelV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannelV_H
#include "ConsoleZChannelV.h"
#endif

void ConsoleZChannelV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames)
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	bool enableBiquad1 = D != 0.5;
	bool enableBiquad2 = H != 0.5;

	VstInt32 inFramesToProcess = sampleFrames; //vst doesn't give us this as a separate variable so we'll make it
	double sampleRate = getSampleRate() * 2.0;
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= sampleRate;

	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;

	double KK;
	double norm;

	// Hypersonic

	double cutoff = 25000.0 / sampleRate;
	if (cutoff > 0.45) cutoff = 0.45; //don't crash if run at 44.1k

	fixG[fix_freq] = fixE[fix_freq] = fixD[fix_freq] = fixA[fix_freq] = cutoff;
	
  fixA[fix_reso] = 4.46570214;
	fixD[fix_reso] = 0.70710678;
	fixE[fix_reso] = 0.59051105;
	fixG[fix_reso] = 0.50316379;
	
	KK = tan(M_PI * fixA[fix_freq]); //lowpass
	norm = 1.0 / (1.0 + KK / fixA[fix_reso] + KK * KK);
	fixA[fix_a0] = KK * KK * norm;
	fixA[fix_a1] = 2.0 * fixA[fix_a0];
	fixA[fix_a2] = fixA[fix_a0];
	fixA[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixA[fix_b2] = (1.0 - KK / fixA[fix_reso] + KK * KK) * norm;
	
	KK = tan(M_PI * fixD[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixD[fix_reso] + KK * KK);
	fixD[fix_a0] = KK * KK * norm;
	fixD[fix_a1] = 2.0 * fixD[fix_a0];
	fixD[fix_a2] = fixD[fix_a0];
	fixD[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixD[fix_b2] = (1.0 - KK / fixD[fix_reso] + KK * KK) * norm;
	
	KK = tan(M_PI * fixE[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixE[fix_reso] + KK * KK);
	fixE[fix_a0] = KK * KK * norm;
	fixE[fix_a1] = 2.0 * fixE[fix_a0];
	fixE[fix_a2] = fixE[fix_a0];
	fixE[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixE[fix_b2] = (1.0 - KK / fixE[fix_reso] + KK * KK) * norm;
	
	KK = tan(M_PI * fixG[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixG[fix_reso] + KK * KK);
	fixG[fix_a0] = KK * KK * norm;
	fixG[fix_a1] = 2.0 * fixG[fix_a0];
	fixG[fix_a2] = fixG[fix_a0];
	fixG[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixG[fix_b2] = (1.0 - KK / fixG[fix_reso] + KK * KK) * norm;

	// BiquadOneHalf 1

	int type_1 = ceil((A*3.999)+0.00001);

	double cutoff_1 = ((B*B*B*0.9999)+0.0001)*0.499;
	if (cutoff_1 < 0.0001) cutoff_1 = 0.0001;

  double res_1 = (C*C*C*29.99)+0.01;
	if (res_1 < 0.0001) res_1 = 0.0001;

	double wet_1 = (D*2.0)-1.0;

	double b0_1=0.0;
	double b1_1=0.0;
	double b2_1=0.0;
	double b3_1=0.0;
	double b4_1=0.0;

	if (type_1 == 1) { //lowpass
		KK = tan(M_PI * cutoff_1);
		norm = 1.0 / (1.0 + KK / res_1 + KK * KK);
		b0_1 = KK * KK * norm;
		b1_1 = 2.0 * b0_1;
		b2_1 = b0_1;
		b3_1 = 2.0 * (KK * KK - 1.0) * norm;
		b4_1 = (1.0 - KK / res_1 + KK * KK) * norm;
	}

	if (type_1 == 2) { //highpass
		KK = tan(M_PI * cutoff_1);
		norm = 1.0 / (1.0 + KK / res_1 + KK * KK);
		b0_1 = norm;
		b1_1 = -2.0 * b0_1;
		b2_1 = b0_1;
		b3_1 = 2.0 * (KK * KK - 1.0) * norm;
		b4_1 = (1.0 - KK / res_1 + KK * KK) * norm;
	}

	if (type_1 == 3) { //bandpass
		KK = tan(M_PI * cutoff_1);
		norm = 1.0 / (1.0 + KK / res_1 + KK * KK);
		b0_1 = KK / res_1 * norm;
		b1_1 = 0.0; //bandpass can simplify the biquad kernel: leave out this multiply
		b2_1 = -b0_1;
		b3_1 = 2.0 * (KK * KK - 1.0) * norm;
		b4_1 = (1.0 - KK / res_1 + KK * KK) * norm;
	}

	if (type_1 == 4) { //notch
		KK = tan(M_PI * cutoff_1);
		norm = 1.0 / (1.0 + KK / res_1 + KK * KK);
		b0_1 = (1.0 + KK * KK) * norm;
		b1_1 = 2.0 * (KK * KK - 1) * norm;
		b2_1 = b0_1;
		b3_1 = b1_1;
		b4_1 = (1.0 - KK / res_1 + KK * KK) * norm;
	}

	// BiquadOneHalf 2

	int type_2 = ceil((E*3.999)+0.00001);

	double cutoff_2 = ((F*F*F*0.9999)+0.0001)*0.499;
	if (cutoff_2 < 0.0001) cutoff_2 = 0.0001;

  double res_2 = (G*G*G*29.99)+0.01;
	if (res_2 < 0.0001) res_2 = 0.0001;

	double wet_2 = (H*2.0)-1.0;

	double b0_2=0.0;
	double b1_2=0.0;
	double b2_2=0.0;
	double b3_2=0.0;
	double b4_2=0.0;

	if (type_2 == 1) { //lowpass
		KK = tan(M_PI * cutoff_2);
		norm = 1.0 / (1.0 + KK / res_2 + KK * KK);
		b0_2 = KK * KK * norm;
		b1_2 = 2.0 * b0_2;
		b2_2 = b0_2;
		b3_2 = 2.0 * (KK * KK - 1.0) * norm;
		b4_2 = (1.0 - KK / res_2 + KK * KK) * norm;
	}

	if (type_2 == 2) { //highpass
		KK = tan(M_PI * cutoff_2);
		norm = 1.0 / (1.0 + KK / res_2 + KK * KK);
		b0_2 = norm;
		b1_2 = -2.0 * b0_2;
		b2_2 = b0_2;
		b3_2 = 2.0 * (KK * KK - 1.0) * norm;
		b4_2 = (1.0 - KK / res_2 + KK * KK) * norm;
	}

	if (type_2 == 3) { //bandpass
		KK = tan(M_PI * cutoff_2);
		norm = 1.0 / (1.0 + KK / res_2 + KK * KK);
		b0_2 = KK / res_2 * norm;
		b1_2 = 0.0; //bandpass can simplify the biquad kernel: leave out this multiply
		b2_2 = -b0_2;
		b3_2 = 2.0 * (KK * KK - 1.0) * norm;
		b4_2 = (1.0 - KK / res_2 + KK * KK) * norm;
	}

	if (type_2 == 4) { //notch
		KK = tan(M_PI * cutoff_2);
		norm = 1.0 / (1.0 + KK / res_2 + KK * KK);
		b0_2 = (1.0 + KK * KK) * norm;
		b1_2 = 2.0 * (KK * KK - 1) * norm;
		b2_2 = b0_2;
		b3_2 = b1_2;
		b4_2 = (1.0 - KK / res_2 + KK * KK) * norm;
	}

	// Desk 4

	double desk4_gain = (pow(I,2)*10)+0.0001;
	double gaintrim = (pow(I,2)*2)+1.0;
	double slewgain = (pow(J,3)*40)+0.0001;
	double prevslew = 0.105;
	double intensity = (pow(K,6)*15)+0.0001;
	double depthA = (pow(L,4)*940)+0.00001;
	int offsetA = (int)(depthA * overallscale);
	if (offsetA < 1) offsetA = 1;
	if (offsetA > 4880) offsetA = 4880;
	double balanceB = 0.0001;
	slewgain *= overallscale;
	prevslew *= overallscale;
	balanceB /= overallscale;
	double outputgain = M;
	double wet = N;
	double balanceA = 1.0 - balanceB;
	//removed extra dry variable

	// PurestDrive

	double pd_intensity = O;

	// ResEQ2

	//begin ResEQ2V Mid Boost
	double freqMPeak = pow(P+0.15,3);
	double amountMPeak = (Q*2.0)-1.0;
	double boostOrCut = amountMPeak >= 0 ? 1.0 : -1.0;
	amountMPeak = pow(amountMPeak,2);
	int maxMPeak = (amountMPeak*63.0)+1;
	if ((freqMPeak != prevfreqMPeak)||(amountMPeak != prevamountMPeak)) {
		for (int x = 0; x < maxMPeak; x+=4) {
			Vec4d xx(x, x + 1, x + 2, x + 3);
			select(
				(xx * freqMPeak) < M_PI_4,
				sin(xx * freqMPeak * 4) * freqMPeak*sin(((maxMPeak-x)/(double)maxMPeak)*M_PI_2),
				cos(xx*freqMPeak)*freqMPeak*sin(((maxMPeak-x)/(double)maxMPeak)*M_PI_2)
			).store_a(f + x);
		}
		prevfreqMPeak = freqMPeak; prevamountMPeak = amountMPeak;
	}//end ResEQ2V Mid Boost

	// ConsoleLA

	int limit = 4*cycleEnd;
	double divisor = 2.0/limit;

	double treble = (R*6.0)-3.0;
	midA = midB;
	midB = (S*6.0)-3.0;
	bassA = bassB;
	bassB = (T*6.0)-3.0;
	//these should stack to go up to -3.0 to 3.0
	if (treble < 0.0) treble /= 3.0;
	if (midB < 0.0) midB /= 3.0;
	if (bassB < 0.0) bassB /= 3.0;
	//and then become -1.0 to 3.0;
	//there will be successive sin/cos stages w. dry/wet in these
	double freqMid = 0.024;
	switch (cycleEnd)
	{
		case 1: //base sample rate, no change
			break;
		case 2: //96k tier
			freqMid = 0.012;
			break;
		case 3: //192k tier
			freqMid = 0.006;
			break;
	}

	gainA = gainB;
	gainB = V*2.0; //smoothed master fader from Z2 filters
	//BitShiftGain pre gain trim goes here

	double subTrim = 0.0011 / overallscale;

	// CStrip 2 comp

	double fpOld = 0.618033988749894848204586; //golden ratio!
	double fpNew = 1.0 - fpOld;

	//begin ButterComp
	double inputgain = (pow(W,4)*35)+1.0;
	double compoutgain = inputgain;
	compoutgain -= 1.0;
	compoutgain /= 1.2;
	compoutgain += 1.0;
	double comp_divisor = (0.008 * pow(X,2))+0.0004;
	//originally 0.012
	comp_divisor /= overallscale;
	double remainder = comp_divisor;
	comp_divisor = 1.0 - comp_divisor;
	bool engageComp = false;
	if (inputgain > 1.0) engageComp = true;
	//end ButterComp

	// Wider

	double densityside = (Y*2.0)-1.0;
	double densitymid = (Z*2.0)-1.0;
	//removed extra dry variable
	double offset = (densityside-densitymid)/2;
	if (offset > 0) offset = sin(offset);
	if (offset < 0) offset = -sin(-offset);
	offset = -(pow(offset,4) * 20 * overallscale);
	int near = (int)floor(fabs(offset));
	double farLevel = fabs(offset) - near;
	int far = near + 1;
	double nearLevel = 1.0 - farLevel;
	

	Vec4ui fpd; fpd.load(fpd_b);
	Vec4d g(gainL, gainR, gainL, gainR);
	Vec4d t, drySample;

  while (--sampleFrames >= 0)
  {
		Vec4d inputSample(*in1, *in2, 0.0, 0.0);
		inputSample *= 2.0;
		inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 

		// Hypersonic
		{
	    Vec2d sG1; sG1.load(fixG + fix_sL1);	
	    Vec2d sG2; sG2.load(fixG + fix_sL2);	
			Vec2d s0 = inputSample.get_low();
      Vec2d temp = (s0 * fixG[fix_a0]) + sG1;
      sG1 = (s0 * fixG[fix_a1]) - (temp * fixG[fix_b1]) + sG2;
      sG2 = (s0 * fixG[fix_a2]) - (temp * fixG[fix_b2]);
      s0 = temp; //fixed biquad filtering ultrasonics
			Vec2d s1 = inputSample.get_high();
      temp = (s1 * fixG[fix_a0]) + sG1;
      sG1 = (s1 * fixG[fix_a1]) - (temp * fixG[fix_b1]) + sG2;
      sG2 = (s1 * fixG[fix_a2]) - (temp * fixG[fix_b2]);
			inputSample = concatenate2(s0, temp);//fixed biquad filtering ultrasonics
		
	    sG1.store(fixG + fix_sL1);	
	    sG2.store(fixG + fix_sL2);	
		}


	  inputSample += ((pow_const(inputSample,5u)/128.0) + (pow_const(inputSample,9u)/262144.0)) - ((pow_const(inputSample,3u)/8.0) + (pow_const(inputSample,7u)/4096.0));

	  if(enableBiquad1) {
		 	Vec4d biquadA; biquadA.load_a(biquad_1);
		 	Vec4d biquadB; biquadB.load_a(biquad_1 + 4);

			drySample = inputSample;
		 	t = inputSample * b0_1 + biquadA;
		 	biquadA = (inputSample * b1_1) - (t * b3_1) + biquadB;
		 	biquadB = (inputSample * b2_1) - (t * b4_1);
		 	inputSample = min(max(t, -1.0), 1.0);

		 	biquadA.store_a(biquad_1);
		 	biquadB.store_a(biquad_1 + 4);

			if (wet_1 < 1.0) {
				inputSample = (inputSample*wet_1) + (drySample*(1.0-fabs(wet_1)));
				//inv/dry/wet lets us turn LP into HP and band into notch
			}
	  }

	  if(enableBiquad2) {
		 	Vec4d biquadA; biquadA.load_a(biquad_2);
		 	Vec4d biquadB; biquadB.load_a(biquad_2 + 4);

			drySample = inputSample;
		 	t = inputSample * b0_2 + biquadA;
		 	biquadA = (inputSample * b1_2) - (t * b3_2) + biquadB;
		 	biquadB = (inputSample * b2_2) - (t * b4_2);
		 	inputSample = min(max(t, -1.0), 1.0);

		 	biquadA.store_a(biquad_2);
		 	biquadB.store_a(biquad_2 + 4);

			if (wet_2 < 1.0) {
				inputSample = (inputSample*wet_2) + (drySample*(1.0-fabs(wet_2)));
				//inv/dry/wet lets us turn LP into HP and band into notch
			}
	  }

		// PurestDrive
		{
			Vec4d prev; prev.load_a(pd_previousSample);
			drySample = inputSample;
	
			inputSample = sin(inputSample);
			inputSample.store_a(pd_previousSample);
			//basic distortion factor
			t = (abs(prev + inputSample) / 2.0) * pd_intensity;
			//saturate less if previous sample was undistorted and low level, or if it was
			//inverse polarity. Lets through highs and brightness more.
			inputSample = (drySample * (1.0 - t)) + (inputSample * t);		
			//dry-wet control for intensity also has FM modulation to clean up highs
		}	
	
		 inputSample += (pow_const(inputSample,3u)/4.0)+(pow_const(inputSample,5u)/8.0)+(pow_const(inputSample,7u)/16.0)+(pow_const(inputSample,9u)/32.0);
		// amplitude aspect

		// Hypersonic
		{
	    Vec2d sE1; sE1.load(fixE + fix_sL1);	
	    Vec2d sE2; sE2.load(fixE + fix_sL2);	
			Vec2d s0 = inputSample.get_low();
      Vec2d temp = (s0 * fixE[fix_a0]) + sE1;
      sE1 = (s0 * fixE[fix_a1]) - (temp * fixE[fix_b1]) + sE2;
      sE2 = (s0 * fixE[fix_a2]) - (temp * fixE[fix_b2]);
      s0 = temp; //fixed biquad filtering ultrasonics
			Vec2d s1 = inputSample.get_high();
      temp = (s1 * fixE[fix_a0]) + sE1;
      sE1 = (s1 * fixE[fix_a1]) - (temp * fixE[fix_b1]) + sE2;
      sE2 = (s1 * fixE[fix_a2]) - (temp * fixE[fix_b2]);
			inputSample = concatenate2(s0, temp);//fixed biquad filtering ultrasonics
		
	    sE1.store(fixE + fix_sL1);	
	    sE2.store(fixE + fix_sL2);	
		}

		// Desk 4

		if(wet > 0.0) {
			Vec2d control; control.load(control_b);
			Vec2d lastSample; lastSample.load(lastSample_b);
			Vec2d lastOutSample; lastOutSample.load(lastOutSample_b);
			Vec2d lastSlew; lastSlew.load(lastSlew_b);
			Vec2d inputSample1, inputSample2;
			Vec4d a, b, c, d;

			if (gcount > 19600) {gcount = 9800;}

			drySample = inputSample;
			a = abs(inputSample)*intensity;
			a.store(d_b + gcount); a.store(d_b + gcount - 9800);
			b = a / offsetA;
			c.load(d_b + gcount - (offsetA * 2));
			b -= c / offsetA;
			b -= 0.000001;
			c = concatenate2(control + b.get_low(), control + b.get_low() + b.get_high()); // controls
			c = max(c, 0.0);
			d = max(select(c > 1.0, 1.0 - (c - 1.0), 1.0), 0.5); // clamp
			c = min(c, 1.0);
			control = c.get_high(); // new control
			a = ((1.0 - c) * 2.0) - 1.0; // thickness
			b = abs(a); // out
			c = min(abs(inputSample), 1.57079633); // bridgerectifier
			c = select(c > 0, sin(c), 1-cos(c));
			inputSample = select(inputSample > 0, (inputSample*(1-b))+(c*b), (inputSample*(1-b))-(c*b));
			inputSample *= d;
			a = concatenate2(lastSample, inputSample.get_low()); // lastSamples
			lastSample = inputSample.get_high(); // new lastSample
			b = inputSample - a; // slew
			c = abs(b*slewgain); // bridgerectifier
			c = select(c > 1.57079633, 1.0, sin(c));
			b = select(b > 0.0, c / slewgain, (c / slewgain) * -1.0);
			inputSample1 = (lastOutSample * balanceA) + (a.get_low() * balanceB) + b.get_low();
			lastOutSample = inputSample1;
			inputSample2 = (lastOutSample * balanceA) + (a.get_high() * balanceB) + b.get_high();
			lastOutSample = inputSample2;
			inputSample = concatenate2(inputSample1, inputSample2);
			d = min(abs(drySample + a), 1.0); // combSample
			c = concatenate2(lastSlew, b.get_low()); // lastSlews
			lastSlew = b.get_high(); // new lastSlew
			inputSample -= (c * d * prevslew);
			inputSample *= desk4_gain;
			c = abs(inputSample); // bridgerectifier
			c = select(c > 1.57079633, 1.0, sin(c));
			inputSample = select(inputSample > 0.0, c, c * -1.0);
			inputSample /= desk4_gain;
			inputSample *= gaintrim;

			if (outputgain != 1.0) {
				inputSample *= outputgain;
			}

			if (wet !=1.0) {
				inputSample = (inputSample * wet) + (drySample * (1.0-wet));
			}

			control.store(control_b);
			lastSample.store(lastSample_b);
			lastOutSample.store(lastOutSample_b);
			lastSlew.store(lastSlew_b);
		}

		// Hypersonic
		{
	    Vec2d sA1; sA1.load(fixA + fix_sL1);	
	    Vec2d sA2; sA2.load(fixA + fix_sL2);	
			Vec2d s0 = inputSample.get_low();
      Vec2d temp = (s0 * fixA[fix_a0]) + sA1;
      sA1 = (s0 * fixA[fix_a1]) - (temp * fixA[fix_b1]) + sA2;
      sA2 = (s0 * fixA[fix_a2]) - (temp * fixA[fix_b2]);
      s0 = temp; //fixed biquad filtering ultrasonics
			Vec2d s1 = inputSample.get_high();
      temp = (s1 * fixA[fix_a0]) + sA1;
      sA1 = (s1 * fixA[fix_a1]) - (temp * fixA[fix_b1]) + sA2;
      sA2 = (s1 * fixA[fix_a2]) - (temp * fixA[fix_b2]);
			inputSample = concatenate2(s0, temp);//fixed biquad filtering ultrasonics
		
	    sA1.store(fixA + fix_sL1);	
	    sA2.store(fixA + fix_sL2);	
		}

		// CStrip2 comp

		//begin compressor
		if (engageComp)
		{
			inputSample *= inputgain;
			
			Vec2d t; t.load_a(cs_avg); inputSample.get_high().store_a(cs_avg);
			Vec4d avg = concatenate2(t, inputSample.get_low());
			Vec4d input = max((inputSample * fpOld) + (avg * fpNew) + 1.0, 0.0);
			Vec4d outputpos = min(input * 0.5, 1.0);
			Vec4d target; target.load_a(cs_targetpos);
			input *= input;
			target *= comp_divisor;
			target += (input * remainder);
			Vec4d calcpos = pow_const((1.0/target),2);
			target.store_a(cs_targetpos);

			input = max((-inputSample * fpOld) + (-avg * fpNew) + 1.0, 0.0);
			Vec4d outputneg = min(input * 0.5, 1.0);
			target.load_a(cs_targetneg);
			input *= input;
			target *= comp_divisor;
			target += (input * remainder);
			Vec4d calcneg = pow_const((1.0/target),2);
			target.store_a(cs_targetneg);
			
			Vec4d controlpos; controlpos.load_a(cs_controlpos);
			Vec4d controlneg; controlneg.load_a(cs_controlneg);
			
			controlpos = select(inputSample > 0.0, (controlpos * comp_divisor) + (calcpos * remainder), controlpos);
			controlneg = select(inputSample < 0.0, (controlneg * comp_divisor) + (calcneg * remainder), controlneg);

			controlpos.store_a(cs_controlpos);
			controlneg.store_a(cs_controlneg);

			Vec4d totalmultiplier = (controlpos * outputpos) + (controlneg * outputneg);
			//this combines the sides according to flip, blending relative to the input value
			
			inputSample *= totalmultiplier;
			inputSample /= compoutgain;
		}
		//end compressor

		// ConsoleLA
		{
			double temp = (double)sampleFrames/inFramesToProcess;
			double gain = (gainA*temp)+(gainB*(1.0-temp));
			double mid = (midA*temp)+(midB*(1.0-temp));
			double bass = (bassA*temp)+(bassB*(1.0-temp));
		
			Vec4d b(0.0);
			//begin Hull2 Treble
			inputSample.store_a(hull + hullp);
			if(hullp == 0) inputSample.store_a(hull + 256);
			int x = hullp;
			while (x > hullp-limit) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(limit * 2)) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			b.store_a(hull + ((hullp - 40) & 255));
			if(((hullp - 40) & 255) == 0) b.store_a(hull + 256);
			x = hullp-40;
			b = 0.0;
			while (x > hullp-(40+limit)) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(40+(limit * 2))) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			b.store_a(hull + ((hullp - 80) & 255));
			if(((hullp - 80) & 255) == 0) b.store_a(hull + 256);
			x = hullp-80;
			b = 0.0;
			while (x > hullp-(80+limit)) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(80+(limit * 2))) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			t = inputSample - b; inputSample = b;
			hullp = (hullp + 4) & 255;

			//end Hull2 treble
		
			//begin Pear filter stages
			{
				//at this point 'bass' is actually still mid and bass
				Vec2d slew; Vec2d b0 = b.get_low(); Vec2d b1 = b.get_high();
				Vec2d p0, p1;
				p0.load_a(pearB);p1.load_a(pearB + 2);slew = ((b0 - p0) + p1)*freqMid*0.5;
				b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
				p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
				b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
				b1.store_a(pearB); slew.store_a(pearB + 2);

				p0.load_a(pearB + 4);p1.load_a(pearB + 6);slew = ((b0 - p0) + p1)*freqMid*0.5;
				b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
				p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
				b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
				b1.store_a(pearB + 4); slew.store_a(pearB + 6);

				p0.load_a(pearB + 8);p1.load_a(pearB + 10); slew = ((b0 - p0) + p1)*freqMid*0.5;
				b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
				p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
				b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
				b1.store_a(pearB + 8); slew.store_a(pearB + 10);

				p0.load_a(pearB + 12);p1.load_a(pearB + 14); slew = ((b0 - p0) + p1)*freqMid*0.5;
				b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
				p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
				b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
				b1.store_a(pearB + 12); slew.store_a(pearB + 14);

				p0.load_a(pearB + 16);p1.load_a(pearB + 18); slew = ((b0 - p0) + p1)*freqMid*0.5;
				b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
				p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
				b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
				b1.store_a(pearB + 16); slew.store_a(pearB + 18);

				b = concatenate2(b0, b1);
			}
			Vec4d m = inputSample - b;
			//we now have three bands out of hull and pear filters
		
			double w = 0.0; //filter into bands, apply the sin/cos to each band
			if (treble > 0.0) {
				w = treble > 1.0 ? 1.0 : treble;
				t = (t*(1.0-w)) + (sin(t*M_PI_2)*treble);
			}
			if (treble < 0.0) {
				t = min(max(t, -1.0), 1.0);
				w = -treble > 1.0 ? 1.0 : -treble;
				t = (t*(1.0-w))+if_mul(t < 0, ((1.0-cos(abs(t)*w))*(1.0-w)), -1.0);
			} //cosine stages for EQ or expansion
		
			m = min(max(m, -1.0), 1.0);
			if (mid > 0.0) {
				w = mid > 1.0 ? 1.0 : mid;
				m = (m*(1.0-w)) + (sin(m*M_PI_2)*mid);
			}
			if (mid < 0.0) {
				w = -mid > 1.0 ? 1.0 : -mid;
				m = (m*(1.0-w))+if_mul(m < 0, ((1.0-cos(abs(m)*w))*(1.0-w)), -1.0);
			} //cosine stages for EQ or expansion
		
			b = min(max(b, -1.0), 1.0);
			if (bass > 0.0) {
				w = bass > 1.0 ? 1.0 : bass;
				b = (b*(1.0-w)) + (sin(b*M_PI_2)*bass);
			}
			if (bass < 0.0) {
				w = -bass > 1.0 ? 1.0 : -bass;
				b = (b*(1.0-w))+if_mul(b < 0, ((1.0-cos(abs(b)*w))*(1.0-w)), -1.0);
			} //cosine stages for EQ or expansion
		
	
			inputSample = (b + m + t)*g*gain;
			//applies BitShiftPan pan section, and smoothed fader gain
		}

		//begin ResEQ2V Mid Boost
		{
			inputSample.store_a(mpk + mpc);
			if(mpc == 0) inputSample.store_a(mpk + 512);
			Vec4d midMPeak(0.0), m,k;
		
			for (int x = 0; x < maxMPeak; x++) {
				int y = x*cycleEnd;
				m = 0;
				for(int yy = y; yy < (y + cycleEnd); yy++) {
					int offset = (mpc - (yy*2)) & 511;
					k.load(mpk + offset);
					m += k;
				} 
				switch (cycleEnd)
				{
					case 1: m *= f[x];	break;
					case 2: m *= (f[x] * 0.5);	break;
					case 3: m *= (f[x] * 0.333);	break;
					case 4: m *= (f[x] * 0.25); break;
				}
				midMPeak += m;
			}
			mpc = (mpc + 4) & 511;
		
			inputSample = (midMPeak*(amountMPeak*boostOrCut))+((1.5-amountMPeak>1.0)?inputSample:inputSample*(1.5-amountMPeak));
		}
		//end ResEQ2V Mid Boost

		// Wider

		if(densityside != 0.0 || densitymid != 0.0) {
			//assign working variables
			Vec4d density(densitymid,densitymid,densityside,densityside);
		  Vec4d out = abs(density);
			Vec4d perm = permute4<0, 2, 0, 2>(inputSample);
			Vec4d ms = perm + (permute4<1, 3, 1, 3>(inputSample) * Vec4d(1, 1, -1, -1));
			Vec4d bridgerectifier = min(abs(ms)*1.57079633,1.57079633);
			bridgerectifier = select(density > 0, sin(bridgerectifier), 1-cos(bridgerectifier)) * out;
			ms = select(density != 0.0, (ms * (1.0-out)) + if_mul(ms < 0, bridgerectifier, -1.0), ms);
			//assign mid and side. Now, High Impact code

			if(offset != 0) {
				Vec2d t = offset > 0 ? ms.get_low() : ms.get_high();
				t.store_a(p + count); if(count == 0) t.store_a(p + 2048);
				Vec2d n; n.load(p + ((count - near) & 2047)); Vec2d f; f.load(p + ((count - far) & 2047)); 
				ms = offset > 0 ? concatenate2((n * nearLevel) + (f * farLevel), ms.get_high()) : concatenate2(ms.get_low(), (n * nearLevel) + (f * farLevel));
			}
			count = (count + 2) & 2047;
		
			inputSample = (permute4<0, 0, 1, 1>(ms) + (permute4<2, 2, 3, 3>(ms) * Vec4d(1, -1, 1, -1))) * 0.5;
		}

		//begin SubTight section
		{	
			Vec2d subA; subA.load_a(subA_b);
			Vec2d subB; subB.load_a(subB_b);
			Vec2d subC; subC.load_a(subC_b);
			Vec4d subSample = inputSample * subTrim;
			Vec2d s0 = subSample.get_low(); Vec2d s1 = subSample.get_high();
			Vec2d scale = 0.5+abs(s0*0.5);
			s0 = (subA+(sin(subA-s0)*scale));
			subA = s0*scale;
			scale = 0.5+abs(s0*0.5);
			s0 = (subB+(sin(subB-s0)*scale));
			subB = s0*scale;
			scale = 0.5+abs(s0*0.5);
			s0 = (subC+(sin(subC-s0)*scale));
			subC = s0*scale;

			scale = 0.5+abs(s1*0.5);
			s1 = (subA+(sin(subA-s1)*scale));
			subA = s1*scale;
			scale = 0.5+abs(s1*0.5);
			s1 = (subB+(sin(subB-s1)*scale));
			subB = s1*scale;
			scale = 0.5+abs(s1*0.5);
			s1 = (subC+(sin(subC-s1)*scale));
			subC = s1*scale;

			subSample = min(max(concatenate2(s0, s1), -0.25), 0.25) * 16.0;
			inputSample += subSample;

			subA.store_a(subA_b);
			subB.store_a(subB_b);
			subC.store_a(subC_b);
		}
		//end SubTight section		

		// Hypersonic
		{
	    Vec2d sD1; sD1.load(fixD + fix_sL1);	
	    Vec2d sD2; sD2.load(fixD + fix_sL2);	
			Vec2d s0 = inputSample.get_low();
      Vec2d temp = (s0 * fixD[fix_a0]) + sD1;
      sD1 = (s0 * fixD[fix_a1]) - (temp * fixD[fix_b1]) + sD2;
      sD2 = (s0 * fixD[fix_a2]) - (temp * fixD[fix_b2]);
      s0 = temp; //fixed biquad filtering ultrasonics
			Vec2d s1 = inputSample.get_high();
      temp = (s1 * fixD[fix_a0]) + sD1;
      sD1 = (s1 * fixD[fix_a1]) - (temp * fixD[fix_b1]) + sD2;
      sD2 = (s1 * fixD[fix_a2]) - (temp * fixD[fix_b2]);
			inputSample = concatenate2(s0, temp);//fixed biquad filtering ultrasonics
		
	    sD1.store(fixD + fix_sL1);	
	    sD2.store(fixD + fix_sL2);	
		}

		//begin Console7 Channel processing
		inputSample = min(max(inputSample, -1.097), 1.097);
		Vec4d absSample = abs(inputSample);
		inputSample = ((sin(inputSample*absSample)/select(absSample == 0.0, 1, absSample))*0.8)+(sin(inputSample)*0.2);
		//this is a version of Spiral blended 80/20 with regular Density.
		//It's blending between two different harmonics in the overtones of the algorithm


		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (to_double(fpd) - 2147483647.0) * to_double(exp) * 5.5e-36l;
		//end 32 bit stereo floating point dither

		double result[2];
		inputSample.get_low().store_a(result);

		*out1 = result[0];
		*out2 = result[1];

		in1++;
		in2++;
		out1++;
		out2++;
	}
	fpd.store(fpd_b);
}

void ConsoleZChannelV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
{
  double* in1  =  inputs[0];
  double* in2  =  inputs[1];
  double* out1 = outputs[0];
  double* out2 = outputs[1];

			// if(flip) {
			// 	//begin 64 bit stereo floating point dither
			// 	//int expon; frexp((double)inputSampleL, &expon);
			// 	fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
			// 	//inputSampleL += ((double(fpdL)-uint32_t(0x7fffffff)) * 1.1e-44l * pow(2,expon+62));
			// 	//frexp((double)inputSampleR, &expon);
			// 	fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;
			// 	//inputSampleR += ((double(fpdR)-uint32_t(0x7fffffff)) * 1.1e-44l * pow(2,expon+62));
			// 	//end 64 bit stereo floating point dither

			// 	*out1 = inputSampleL;
			// 	*out2 = inputSampleR;

			// 	in1++;
			// 	in2++;
			// 	out1++;
			// 	out2++;
			// 	--sampleFrames;
	  //   }
	
			// flip = !flip;
   //  }
}
