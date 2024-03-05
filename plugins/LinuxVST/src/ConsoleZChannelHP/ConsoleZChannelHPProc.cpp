/* ========================================
 *  ConsoleZChannelHPHP - ConsoleZChannelHPHP.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannelHPHP_H
#include "ConsoleZChannelHPHP.h"
#endif

void ConsoleZChannelHPHP::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames)
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
	double cutoff;
	double res;

	// Hypersonic

	double cutoff = 25000.0 / sampleRate;
	if (cutoff > 0.45) cutoff = 0.45; //don't crash if run at 44.1k

	
  double resA = 4.46570214;
	double resD = 0.70710678;
	double resE = 0.59051105;
	double resG = 0.50316379;
	
	KK = tan(M_PI * cutoff); //lowpass
	norm = 1.0 / (1.0 + KK / resA + KK * KK);
	fixA[fix_a0] = KK * KK * norm;
	fixA[fix_a1] = 2.0 * fixA[fix_a0];
	fixA[fix_a2] = fixA[fix_a0];
	fixA[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixA[fix_b2] = (1.0 - KK / resA + KK * KK) * norm;
	
	KK = tan(M_PI * cutoff);
	norm = 1.0 / (1.0 + KK / resD + KK * KK);
	fixD[fix_a0] = KK * KK * norm;
	fixD[fix_a1] = 2.0 * fixD[fix_a0];
	fixD[fix_a2] = fixD[fix_a0];
	fixD[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixD[fix_b2] = (1.0 - KK / resD + KK * KK) * norm;
	
	KK = tan(M_PI * cutoff);
	norm = 1.0 / (1.0 + KK / resE + KK * KK);
	fixE[fix_a0] = KK * KK * norm;
	fixE[fix_a1] = 2.0 * fixE[fix_a0];
	fixE[fix_a2] = fixE[fix_a0];
	fixE[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixE[fix_b2] = (1.0 - KK / resE + KK * KK) * norm;
	
	KK = tan(M_PI * cutoff);
	norm = 1.0 / (1.0 + KK / resG + KK * KK);
	fixG[fix_a0] = KK * KK * norm;
	fixG[fix_a1] = 2.0 * fixG[fix_a0];
	fixG[fix_a2] = fixG[fix_a0];
	fixG[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixG[fix_b2] = (1.0 - KK / resG + KK * KK) * norm;

	// BiquadOneHalf 1

	int type_1 = ceil((A*3.999)+0.00001);

	cutoff = ((B*B*B*0.9999)+0.0001)*0.499;
	if (cutoff < 0.0001) cutoff = 0.0001;

    res = (C*C*C*29.99)+0.01;
	if (res < 0.0001) res = 0.0001;

	double wet_1 = (D*2.0)-1.0;

	if (type_1 == 1) { //lowpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_1[0] = KK * KK * norm;
		biquad_1[1] = 2.0 * biquad_1[2];
		biquad_1[2] = biquad_1[2];
		biquad_1[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_1[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_1 == 2) { //highpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_1[0] = norm;
		biquad_1[1] = -2.0 * biquad_1[2];
		biquad_1[2] = biquad_1[2];
		biquad_1[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_1[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_1 == 3) { //bandpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_1[0] = KK / res * norm;
		biquad_1[1] = 0.0; //bandpass can simplify the biquad_1_ kernel: leave out this multiply
		biquad_1[2] = -biquad_1[2];
		biquad_1[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_1[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_1 == 4) { //notch
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_1[0] = (1.0 + KK * KK) * norm;
		biquad_1[1] = 2.0 * (KK * KK - 1) * norm;
		biquad_1[2] = biquad_1[2];
		biquad_1[3] = biquad_1[3];
		biquad_1[4] = (1.0 - KK / res + KK * KK) * norm;
	}

	// BiquadOneHalf 2

	int type_2 = ceil((E*3.999)+0.00001);

	cutoff = ((F*F*F*0.9999)+0.0001)*0.499;
	if (cutoff < 0.0001) cutoff = 0.0001;

  res = (G*G*G*29.99)+0.01;
	if (res < 0.0001) res = 0.0001;

	double wet_2 = (H*2.0)-1.0;

	if (type_2 == 1) { //lowpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_2[0] = KK * KK * norm;
		biquad_2[1] = 2.0 * biquad_2[2];
		biquad_2[2] = biquad_2[2];
		biquad_2[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_2[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_2 == 2) { //highpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_2[0] = norm;
		biquad_2[1] = -2.0 * biquad_2[2];
		biquad_2[2] = biquad_2[2];
		biquad_2[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_2[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_2 == 3) { //bandpass
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_2[0] = KK / res * norm;
		biquad_2[1] = 0.0; //bandpass can simplify the biquad_2_ kernel: leave out this multiply
		biquad_2[2] = -biquad_2[2];
		biquad_2[3] = 2.0 * (KK * KK - 1.0) * norm;
		biquad_2[4] = (1.0 - KK / res + KK * KK) * norm;
	}
	
	if (type_2 == 4) { //notch
		KK = tan(M_PI * cutoff);
		norm = 1.0 / (1.0 + KK / res + KK * KK);
		biquad_2[0] = (1.0 + KK * KK) * norm;
		biquad_2[1] = 2.0 * (KK * KK - 1) * norm;
		biquad_2[2] = biquad_2[2];
		biquad_2[3] = biquad_2[3];
		biquad_2[4] = (1.0 - KK / res + KK * KK) * norm;
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
	//removed extra dry variable

	double clampL;
	double clampR;
	double thicknessL;
	double thicknessR;
	double out;
	double balanceA = 1.0 - balanceB;
	double bridgerectifier;
	double slewL;
	double slewR;
	double combSampleL;
	double combSampleR;

	// PurestDrive

	double pd_intensity = O;
	double pd_apply;

	// ResEQ2

	//begin ResEQ2 Mid Boost
	double freqMPeak = pow(P+0.15,3);
	double amountMPeak = pow(Q,2);
	int maxMPeak = (amountMPeak*63.0)+1;
	if ((freqMPeak != prevfreqMPeak)||(amountMPeak != prevamountMPeak)) {
		for (int x = 0; x < maxMPeak; x++) {
			if (((double)x*freqMPeak) < M_PI_4) f[x] = sin(((double)x*freqMPeak)*4.0)*freqMPeak*sin(((double)(maxMPeak-x)/(double)maxMPeak)*M_PI_2);
			else f[x] = cos((double)x*freqMPeak)*freqMPeak*sin(((double)(maxMPeak-x)/(double)maxMPeak)*M_PI_2);
		}
		prevfreqMPeak = freqMPeak; prevamountMPeak = amountMPeak;
	}//end ResEQ2 Mid Boost

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


	int bitshiftL = 0;
	int bitshiftR = 0;
	double panControl = (U*2.0)-1.0; //-1.0 to 1.0
	double panAttenuation = (1.0-fabs(panControl));
	int panBits = 20; //start centered
	if (panAttenuation > 0.0) panBits = floor(1.0 / panAttenuation);
	if (panControl > 0.25) bitshiftL += panBits;
	if (panControl < -0.25) bitshiftR += panBits;
	if (bitshiftL < 0) bitshiftL = 0; if (bitshiftL > 17) bitshiftL = 17;
	if (bitshiftR < 0) bitshiftR = 0; if (bitshiftR > 17) bitshiftR = 17;
	double gainL = 1.0;
	double gainR = 1.0;
	switch (bitshiftL)
	{
		case 17: gainL = 0.0; break;
		case 16: gainL = 0.0000152587890625; break;
		case 15: gainL = 0.000030517578125; break;
		case 14: gainL = 0.00006103515625; break;
		case 13: gainL = 0.0001220703125; break;
		case 12: gainL = 0.000244140625; break;
		case 11: gainL = 0.00048828125; break;
		case 10: gainL = 0.0009765625; break;
		case 9: gainL = 0.001953125; break;
		case 8: gainL = 0.00390625; break;
		case 7: gainL = 0.0078125; break;
		case 6: gainL = 0.015625; break;
		case 5: gainL = 0.03125; break;
		case 4: gainL = 0.0625; break;
		case 3: gainL = 0.125; break;
		case 2: gainL = 0.25; break;
		case 1: gainL = 0.5; break;
		case 0: break;
	}
	switch (bitshiftR)
	{
		case 17: gainR = 0.0; break;
		case 16: gainR = 0.0000152587890625; break;
		case 15: gainR = 0.000030517578125; break;
		case 14: gainR = 0.00006103515625; break;
		case 13: gainR = 0.0001220703125; break;
		case 12: gainR = 0.000244140625; break;
		case 11: gainR = 0.00048828125; break;
		case 10: gainR = 0.0009765625; break;
		case 9: gainR = 0.001953125; break;
		case 8: gainR = 0.00390625; break;
		case 7: gainR = 0.0078125; break;
		case 6: gainR = 0.015625; break;
		case 5: gainR = 0.03125; break;
		case 4: gainR = 0.0625; break;
		case 3: gainR = 0.125; break;
		case 2: gainR = 0.25; break;
		case 1: gainR = 0.5; break;
		case 0: break;
	}

	gainA = gainB;
	gainB = V*2.0; //smoothed master fader from Z2 filters
	//BitShiftGain pre gain trim goes here

	double subTrim = 0.0011 / overallscale;

	// CStrip 2 comp

	double fpOld = 0.618033988749894848204586; //golden ratio!
	double fpNew = 1.0 - fpOld;

	//begin ButterComp
	double inputpos;
	double inputneg;
	double calcpos;
	double calcneg;
	double outputpos;
	double outputneg;
	double totalmultiplier;
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

	double mid;
	double side;
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
	

	double inputSample[4];
	double drySample[4];
	double temp[4];

  while (sampleFrames > 0)
  {

		inputSample[0] = fpd[0] * 1.18e-17;
		inputSample[1] = fpd[1] * 1.18e-17;
		inputSample[3] = *in1 * 2.0;
		inputSample[4] = *in2 * 2.0;

		if (fabs(inputSample[3])<1.18e-23) inputSample[3] = fpd[0] * 1.18e-17;
		if (fabs(inputSample[4])<1.18e-23) inputSample[4] = fpd[1] * 1.18e-17;

		// Hypersonic

		for(int i = 0; i < 4; ++i) {
      temp[i] = (inputSample[i] * fixG[fix_a0]) + fixG[fix_sL1 + (i % 2)];
      fixG[fix_sL1 + (i % 2)] = (inputSample[i] * fixG[fix_a1]) - (temp[i] * fixG[fix_b1]) + fixG[fix_sL2 + (i % 2)];
      fixG[fix_sL2 + (i % 2)] = (inputSample[i] * fixG[fix_a2]) - (temp[i] * fixG[fix_b2]);
      inputSample[i] = temp[i]; //fixed biquad filtering ultrasonics

			inputSample[i] += ((pow(inputSample[i],5)/128.0) + (pow(inputSample[i],9)/262144.0)) - ((pow(inputSample[i],3)/8.0) + (pow(inputSample[i],7)/4096.0));

			if(enableBiquad1) {
				drySample[i] = inputSample[i];

				temp = (inputSample[i] * biquad_1[0]) + biquad_1[5 + i];
				biquad_1[5 + i] = (inputSample[i] * biquad_1[1]) - (temp * biquad_1[3]) + biquad_1[9 + 1];
				biquad_1[9 + i] = (inputSample[i] * biquad_1[2]) - (temp * biquad_1[4]);
				inputSample[i] = temp;

				if (wet_1 < 1.0) {
					inputSample[i] = (inputSample[i]*wet_1) + (drySample[i]*(1.0-fabs(wet_1)));
					//inv/dry/wet lets us turn LP into HP and band into notch
				}
			}

			if(enableBiquad2) {
				drySample[i] = inputSample[i];

				temp = (inputSample[i] * biquad_2[0]) + biquad_2[5 + i];
				biquad_2[5 + i] = (inputSample[i] * biquad_2[1]) - (temp * biquad_2[3]) + biquad_2[9 + i];
				biquad_2[9 + i] = (inputSample[i] * biquad_2[2]) - (temp * biquad_2[4]);
				inputSample[i] = temp;

				if (wet_2 < 1.0) {
					inputSample[i] = (inputSample[i]*wet_2) + (drySample[i]*(2.0-fabs(wet_2)));
					//inv/dry/wet lets us turn LP into HP and band into notch
				}
			}

			// PurestDrive

			drySample[i] = inputSample[i];
	
			inputSample[i] = sin(inputSample[i]);
			//basic distortion factor
			pd_apply = (fabs(pd_previousSample[i % 2] + inputSample[i]) / 2.0) * pd_intensity;
			//saturate less if previous sample was undistorted and low level, or if it was
			//inverse polarity. [i]ets through highs and brightness more.
			inputSample[i] = (drySample[i] * (1.0 - pd_apply)) + (inputSample[i] * pd_apply);		
			//dry-wet control for intensity also has FM modulation to clean up highs
			pd_previousSample[i % 2] = sin(drySample[i]);
			//apply the sine while storing previous sample
	
	
			inputSample[i] += (pow(inputSample[i],3)/4.0)+(pow(inputSample[i],5)/8.0)+(pow(inputSample[i],7)/16.0)+(pow(inputSample[i],9)/32.0);

			// Hypersonic

			temp = (inputSampleL * fixE[fix_a0]) + fixE[fix_sL1];
			fixE[fix_sL1] = (inputSampleL * fixE[fix_a1]) - (temp * fixE[fix_b1]) + fixE[fix_sL2];
			fixE[fix_sL2] = (inputSampleL * fixE[fix_a2]) - (temp * fixE[fix_b2]);
			inputSampleL = temp; //fixed biquad filtering ultrasonics
			temp = (inputSampleR * fixE[fix_a0]) + fixE[fix_sR1];
			fixE[fix_sR1] = (inputSampleR * fixE[fix_a1]) - (temp * fixE[fix_b1]) + fixE[fix_sR2];
			fixE[fix_sR2] = (inputSampleR * fixE[fix_a2]) - (temp * fixE[fix_b2]);
			inputSampleR = temp; //fixed biquad filtering ultrasonics
		}

		// Desk 4

		if(wet > 0.0) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			if (gcount < 0 || gcount > 4900) {gcount = 4900;}

			dL[gcount+4900] = dL[gcount] = fabs(inputSampleL)*intensity;
			controlL += (dL[gcount] / offsetA);
			controlL -= (dL[gcount+offsetA] / offsetA);
			controlL -= 0.000001;
			clampL = 1;
			if (controlL < 0) {controlL = 0;}
			if (controlL > 1) {clampL -= (controlL - 1); controlL = 1;}
			if (clampL < 0.5) {clampL = 0.5;}

			dR[gcount+4900] = dR[gcount] = fabs(inputSampleR)*intensity;
			controlR += (dR[gcount] / offsetA);
			controlR -= (dR[gcount+offsetA] / offsetA);
			controlR -= 0.000001;
			clampR = 1;
			if (controlR < 0) {controlR = 0;}
			if (controlR > 1) {clampR -= (controlR - 1); controlR = 1;}
			if (clampR < 0.5) {clampR = 0.5;}


			gcount--;
			//control = 0 to 1
			thicknessL = ((1.0 - controlL) * 2.0) - 1.0;
			thicknessR = ((1.0 - controlR) * 2.0) - 1.0;

			out = fabs(thicknessL);
			bridgerectifier = fabs(inputSampleL);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.57079633;
			//max value for sine function
			if (thicknessL > 0) bridgerectifier = sin(bridgerectifier);
			else bridgerectifier = 1-cos(bridgerectifier);
			//produce either boosted or starved version
			if (inputSampleL > 0) inputSampleL = (inputSampleL*(1-out))+(bridgerectifier*out);
			else inputSampleL = (inputSampleL*(1-out))-(bridgerectifier*out);
			//blend according to density control

			out = fabs(thicknessR);
			bridgerectifier = fabs(inputSampleR);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.57079633;
			//max value for sine function
			if (thicknessR > 0) bridgerectifier = sin(bridgerectifier);
			else bridgerectifier = 1-cos(bridgerectifier);
			//produce either boosted or starved version
			if (inputSampleR > 0) inputSampleR = (inputSampleR*(1-out))+(bridgerectifier*out);
			else inputSampleR = (inputSampleR*(1-out))-(bridgerectifier*out);
			//blend according to density control

			inputSampleL *= clampL;
			inputSampleR *= clampR;

			slewL = inputSampleL - lastSampleL;
			lastSampleL = inputSampleL;
			//Set up direct reference for slew

			slewR = inputSampleR - lastSampleR;
			lastSampleR = inputSampleR;
			//Set up direct reference for slew

			bridgerectifier = fabs(slewL*slewgain);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (slewL > 0) slewL = bridgerectifier/slewgain;
			else slewL = -(bridgerectifier/slewgain);

			bridgerectifier = fabs(slewR*slewgain);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (slewR > 0) slewR = bridgerectifier/slewgain;
			else slewR = -(bridgerectifier/slewgain);

			inputSampleL = (lastOutSampleL*balanceA) + (lastSampleL*balanceB) + slewL;
			//go from last slewed, but include some raw values
			lastOutSampleL = inputSampleL;
			//Set up slewed reference

			inputSampleR = (lastOutSampleR*balanceA) + (lastSampleR*balanceB) + slewR;
			//go from last slewed, but include some raw values
			lastOutSampleR = inputSampleR;
			//Set up slewed reference

			combSampleL = fabs(drySampleL*lastSampleL);
			if (combSampleL > 1.0) combSampleL = 1.0;
			//bailout for very high input gains

			combSampleR = fabs(drySampleR*lastSampleR);
			if (combSampleR > 1.0) combSampleR = 1.0;
			//bailout for very high input gains

			inputSampleL -= (lastSlewL * combSampleL * prevslew);
			lastSlewL = slewL;
			//slew interaction with previous slew

			inputSampleR -= (lastSlewR * combSampleR * prevslew);
			lastSlewR = slewR;
			//slew interaction with previous slew

			inputSampleL *= desk4_gain;
			bridgerectifier = fabs(inputSampleL);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleL > 0) inputSampleL = bridgerectifier;
			else inputSampleL = -bridgerectifier;
			//drive section
			inputSampleL /= desk4_gain;
			inputSampleL *= gaintrim;
			//end of Desk section

			inputSampleR *= desk4_gain;
			bridgerectifier = fabs(inputSampleR);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleR > 0) inputSampleR = bridgerectifier;
			else inputSampleR = -bridgerectifier;
			//drive section
			inputSampleR /= desk4_gain;
			inputSampleR *= gaintrim;
			//end of Desk section

			if (outputgain != 1.0) {
				inputSampleL *= outputgain;
				inputSampleR *= outputgain;
			}

			if (wet !=1.0) {
				inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0-wet));
				inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0-wet));
			}
		}

		// Hypersonic

		temp = (inputSampleL * fixA[fix_a0]) + fixA[fix_sL1];
		fixA[fix_sL1] = (inputSampleL * fixA[fix_a1]) - (temp * fixA[fix_b1]) + fixA[fix_sL2];
		fixA[fix_sL2] = (inputSampleL * fixA[fix_a2]) - (temp * fixA[fix_b2]);
		inputSampleL = temp; //fixed biquad filtering ultrasonics
		temp = (inputSampleR * fixA[fix_a0]) + fixA[fix_sR1];
		fixA[fix_sR1] = (inputSampleR * fixA[fix_a1]) - (temp * fixA[fix_b1]) + fixA[fix_sR2];
		fixA[fix_sR2] = (inputSampleR * fixA[fix_a2]) - (temp * fixA[fix_b2]);
		inputSampleR = temp; //fixed biquad filtering ultrasonics

		// CStrip2 comp

		//begin compressor
		if (engageComp)
		{
			//begin L
			inputSampleL *= inputgain;

			inputpos = (inputSampleL * fpOld) + (avgLA * fpNew) + 1.0;
			avgLA = inputSampleL;

			if (inputpos < 0.0) inputpos = 0.0;
			outputpos = inputpos / 2.0;
			if (outputpos > 1.0) outputpos = 1.0;
			inputpos *= inputpos;
			targetposL *= comp_divisor;
			targetposL += (inputpos * remainder);
			calcpos = pow((1.0/targetposL),2);

			inputneg = (-inputSampleL * fpOld) + (nvgLA * fpNew) + 1.0;
			nvgLA = -inputSampleL;

			if (inputneg < 0.0) inputneg = 0.0;
			outputneg = inputneg / 2.0;
			if (outputneg > 1.0) outputneg = 1.0;
			inputneg *= inputneg;
			targetnegL *= comp_divisor;
			targetnegL += (inputneg * remainder);
			calcneg = pow((1.0/targetnegL),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1

			if (inputSampleL > 0)
			{ //working on pos
				if (true == flip)
				{
					controlAposL *= comp_divisor;
					controlAposL += (calcpos*remainder);

				}
				else
				{
					controlBposL *= comp_divisor;
					controlBposL += (calcpos*remainder);
				}
			}
			else
			{ //working on neg
				if (true == flip)
				{
					controlAnegL *= comp_divisor;
					controlAnegL += (calcneg*remainder);
				}
				else
				{
					controlBnegL *= comp_divisor;
					controlBnegL += (calcneg*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'

			if (true == flip)
			{totalmultiplier = (controlAposL * outputpos) + (controlAnegL * outputneg);}
			else
			{totalmultiplier = (controlBposL * outputpos) + (controlBnegL * outputneg);}
			//this combines the sides according to flip, blending relative to the input value

			inputSampleL *= totalmultiplier;
			inputSampleL /= compoutgain;
			//end L

			//begin R
			inputSampleR *= inputgain;

			inputpos = (inputSampleR * fpOld) + (avgRA * fpNew) + 1.0;
			avgRA = inputSampleR;

			if (inputpos < 0.0) inputpos = 0.0;
			outputpos = inputpos / 2.0;
			if (outputpos > 1.0) outputpos = 1.0;
			inputpos *= inputpos;
			targetposR *= comp_divisor;
			targetposR += (inputpos * remainder);
			calcpos = pow((1.0/targetposR),2);

			inputneg = (-inputSampleR * fpOld) + (nvgRA * fpNew) + 1.0;
			nvgRA = -inputSampleR;

			if (inputneg < 0.0) inputneg = 0.0;
			outputneg = inputneg / 2.0;
			if (outputneg > 1.0) outputneg = 1.0;
			inputneg *= inputneg;
			targetnegR *= comp_divisor;
			targetnegR += (inputneg * remainder);
			calcneg = pow((1.0/targetnegR),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1

			if (inputSampleR > 0)
			{ //working on pos
				if (true == flip)
				{
					controlAposR *= comp_divisor;
					controlAposR += (calcpos*remainder);

				}
				else
				{
					controlBposR *= comp_divisor;
					controlBposR += (calcpos*remainder);
				}
			}
			else
			{ //working on neg
				if (true == flip)
				{
					controlAnegR *= comp_divisor;
					controlAnegR += (calcneg*remainder);
				}
				else
				{
					controlBnegR *= comp_divisor;
					controlBnegR += (calcneg*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'

			if (true == flip)
			{totalmultiplier = (controlAposR * outputpos) + (controlAnegR * outputneg);}
			else
			{totalmultiplier = (controlBposR * outputpos) + (controlBnegR * outputneg);}
			//this combines the sides according to flip, blending relative to the input value

			inputSampleR *= totalmultiplier;
			inputSampleR /= compoutgain;
			//end R
		}
		//end compressor

		// ConsoleLA

		temp = (double)sampleFrames/inFramesToProcess;
		double gain = (gainA*temp)+(gainB*(1.0-temp));
		double mid = (midA*temp)+(midB*(1.0-temp));
		double bass = (bassA*temp)+(bassB*(1.0-temp));

		//begin Hull2 Treble
		hullp--; if (hullp < 0) hullp += 60;
		hullL[hullp] = hullL[hullp+60] = inputSampleL;
		hullR[hullp] = hullR[hullp+60] = inputSampleR;

		int x = hullp;
		double bassL = 0.0;
		double bassR = 0.0;
		while (x < hullp+(limit/2)) {
			bassL += hullL[x] * divisor;
			bassR += hullR[x] * divisor;
			x++;
		}
		bassL += bassL * 0.125;
		bassR += bassR * 0.125;
		while (x < hullp+limit) {
			bassL -= hullL[x] * 0.125 * divisor;
			bassR -= hullR[x] * 0.125 * divisor;
			x++;
		}
		hullL[hullp+20] = hullL[hullp+80] = bassL;
		hullR[hullp+20] = hullR[hullp+80] = bassR;
		x = hullp+20;
		bassL = bassR = 0.0;
		while (x < hullp+20+(limit/2)) {
			bassL += hullL[x] * divisor;
			bassR += hullR[x] * divisor;
			x++;
		}
		bassL += bassL * 0.125;
		bassR += bassR * 0.125;
		while (x < hullp+20+limit) {
			bassL -= hullL[x] * 0.125 * divisor;
			bassR -= hullR[x] * 0.125 * divisor;
			x++;
		}
		hullL[hullp+40] = hullL[hullp+100] = bassL;
		hullR[hullp+40] = hullR[hullp+100] = bassR;
		x = hullp+40;
		bassL = bassR = 0.0;
		while (x < hullp+40+(limit/2)) {
			bassL += hullL[x] * divisor;
			bassR += hullR[x] * divisor;
			x++;
		}
		bassL += bassL * 0.125;
		bassR += bassR * 0.125;
		while (x < hullp+40+limit) {
			bassL -= hullL[x] * 0.125 * divisor;
			bassR -= hullR[x] * 0.125 * divisor;
			x++;
		}
		double trebleL = inputSampleL - bassL; inputSampleL = bassL;
		double trebleR = inputSampleR - bassR; inputSampleR = bassR;
		//end Hull2 treble

		//begin Pear filter stages
		//at this point 'bass' is actually still mid and bass
		double slew = ((bassL - pearB[0]) + pearB[1])*freqMid*0.5;
		pearB[0] = bassL = (freqMid * bassL) + ((1.0-freqMid) * (pearB[0] + pearB[1]));
		pearB[1] = slew; slew = ((bassR - pearB[2]) + pearB[3])*freqMid*0.5;
		pearB[2] = bassR = (freqMid * bassR) + ((1.0-freqMid) * (pearB[2] + pearB[3]));
		pearB[3] = slew; slew = ((bassL - pearB[4]) + pearB[5])*freqMid*0.5;
		pearB[4] = bassL = (freqMid * bassL) + ((1.0-freqMid) * (pearB[4] + pearB[5]));
		pearB[5] = slew; slew = ((bassR - pearB[6]) + pearB[7])*freqMid*0.5;
		pearB[6] = bassR = (freqMid * bassR) + ((1.0-freqMid) * (pearB[6] + pearB[7]));
		pearB[7] = slew; slew = ((bassL - pearB[8]) + pearB[9])*freqMid*0.5;
		pearB[8] = bassL = (freqMid * bassL) + ((1.0-freqMid) * (pearB[8] + pearB[9]));
		pearB[9] = slew; slew = ((bassR - pearB[10]) + pearB[11])*freqMid*0.5;
		pearB[10] = bassR = (freqMid * bassR) + ((1.0-freqMid) * (pearB[10] + pearB[11]));
		pearB[11] = slew; slew = ((bassL - pearB[12]) + pearB[13])*freqMid*0.5;
		pearB[12] = bassL = (freqMid * bassL) + ((1.0-freqMid) * (pearB[12] + pearB[13]));
		pearB[13] = slew; slew = ((bassR - pearB[14]) + pearB[15])*freqMid*0.5;
		pearB[14] = bassR = (freqMid * bassR) + ((1.0-freqMid) * (pearB[14] + pearB[15]));
		pearB[15] = slew; slew = ((bassL - pearB[16]) + pearB[17])*freqMid*0.5;
		pearB[16] = bassL = (freqMid * bassL) + ((1.0-freqMid) * (pearB[16] + pearB[17]));
		pearB[17] = slew; slew = ((bassR - pearB[18]) + pearB[19])*freqMid*0.5;
		pearB[18] = bassR = (freqMid * bassR) + ((1.0-freqMid) * (pearB[18] + pearB[19]));
		pearB[19] = slew;
		double midL = inputSampleL - bassL;
		double midR = inputSampleR - bassR;
		//we now have three bands out of hull and pear filters

		double w = 0.0; //filter into bands, apply the sin/cos to each band
		if (treble > 0.0) {
			w = treble; if (w > 1.0) w = 1.0;
			trebleL = (trebleL*(1.0-w)) + (sin(trebleL*M_PI_2)*treble);
			trebleR = (trebleR*(1.0-w)) + (sin(trebleR*M_PI_2)*treble);
		}
		if (treble < 0.0) {
			if (trebleL > 1.0) trebleL = 1.0; if (trebleL < -1.0) trebleL = -1.0;
			if (trebleR > 1.0) trebleR = 1.0; if (trebleR < -1.0) trebleR = -1.0;
			w = -treble; if (w > 1.0) w = 1.0;
			if (trebleL > 0) trebleL = (trebleL*(1.0-w))+((1.0-cos(trebleL*w))*(1.0-w));
			else trebleL = (trebleL*(1.0-w))+((-1.0+cos(-trebleL*w))*(1.0-w));
			if (trebleR > 0) trebleR = (trebleR*(1.0-w))+((1.0-cos(trebleR*w))*(1.0-w));
			else trebleR = (trebleR*(1.0-w))+((-1.0+cos(-trebleR*w))*(1.0-w));
		} //cosine stages for EQ or expansion

		if (midL > 1.0) midL = 1.0; if (midL < -1.0) midL = -1.0;
		if (midR > 1.0) midR = 1.0; if (midR < -1.0) midR = -1.0;
		if (mid > 0.0) {
			w = mid; if (w > 1.0) w = 1.0;
			midL = (midL*(1.0-w)) + (sin(midL*M_PI_2)*mid);
			midR = (midR*(1.0-w)) + (sin(midR*M_PI_2)*mid);
		}
		if (mid < 0.0) {
			w = -mid; if (w > 1.0) w = 1.0;
			if (midL > 0) midL = (midL*(1.0-w))+((1.0-cos(midL*w))*(1.0-w));
			else midL = (midL*(1.0-w))+((-1.0+cos(-midL*w))*(1.0-w));
			if (midR > 0) midR = (midR*(1.0-w))+((1.0-cos(midR*w))*(1.0-w));
			else midR = (midR*(1.0-w))+((-1.0+cos(-midR*w))*(1.0-w));
		} //cosine stages for EQ or expansion

		if (bassL > 1.0) bassL = 1.0; if (bassL < -1.0) bassL = -1.0;
		if (bassR > 1.0) bassR = 1.0; if (bassR < -1.0) bassR = -1.0;
		if (bass > 0.0) {
			w = bass; if (w > 1.0) w = 1.0;
			bassL = (bassL*(1.0-w)) + (sin(bassL*M_PI_2)*bass);
			bassR = (bassR*(1.0-w)) + (sin(bassR*M_PI_2)*bass);
		}
		if (bass < 0.0) {
			w = -bass; if (w > 1.0) w = 1.0;
			if (bassL > 0) bassL = (bassL*(1.0-w))+((1.0-cos(bassL*w))*(1.0-w));
			else bassL = (bassL*(1.0-w))+((-1.0+cos(-bassL*w))*(1.0-w));
			if (bassR > 0) bassR = (bassR*(1.0-w))+((1.0-cos(bassR*w))*(1.0-w));
			else bassR = (bassR*(1.0-w))+((-1.0+cos(-bassR*w))*(1.0-w));
		} //cosine stages for EQ or expansion

		inputSampleL = (bassL + midL + trebleL)*gainL*gain;
		inputSampleR = (bassR + midR + trebleR)*gainR*gain;
		//applies BitShiftPan pan section, and smoothed fader gain

		//begin ResEQ2 Mid Boost
		mpc++; if (mpc < 1 || mpc > 2001) mpc = 1;
		mpkL[mpc] = inputSampleL;
		mpkR[mpc] = inputSampleR;
		double midMPeakL = 0.0;
		double midMPeakR = 0.0;
		for (int x = 0; x < maxMPeak; x++) {
			int y = x*cycleEnd;
			switch (cycleEnd)
			{
				case 1:
					midMPeakL += (mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x]);
					midMPeakR += (mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x]); break;
				case 2:
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.5);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.5); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.5);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.5); break;
				case 3:
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.333); break;
				case 4:
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25); y--;
					midMPeakL += ((mpkL[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25);
					midMPeakR += ((mpkR[(mpc-y)+((mpc-y < 1)?2001:0)] * f[x])*0.25); //break
			}
		}
		inputSampleL = (midMPeakL*amountMPeak)+((1.5-amountMPeak>1.0)?inputSampleL:inputSampleL*(1.5-amountMPeak));
		inputSampleR = (midMPeakR*amountMPeak)+((1.5-amountMPeak>1.0)?inputSampleR:inputSampleR*(1.5-amountMPeak));
		//end ResEQ2 Mid Boost

		// Wider

		if(densityside != 0.0 || densitymid != 0.0) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;
			//assign working variables		
			mid = inputSampleL + inputSampleR;
			side = inputSampleL - inputSampleR;
			//assign mid and side. Now, High Impact code
	
			if (densityside != 0.0)
			{
				temp = fabs(densityside);
				bridgerectifier = fabs(side)*1.57079633;
				if (bridgerectifier > 1.57079633) bridgerectifier = 1.57079633;
				//max value for sine function
				if (densityside > 0) bridgerectifier = sin(bridgerectifier);
				else bridgerectifier = 1-cos(bridgerectifier);
				//produce either boosted or starved version
				if (side > 0) side = (side*(1-temp))+(bridgerectifier*temp);
				else side = (side*(1-temp))-(bridgerectifier*temp);
				//blend according to density control
			}
	
			if (densitymid != 0.0)
			{
				temp = fabs(densitymid);
				bridgerectifier = fabs(mid)*1.57079633;
				if (bridgerectifier > 1.57079633) bridgerectifier = 1.57079633;
				//max value for sine function
				if (densitymid > 0) bridgerectifier = sin(bridgerectifier);
				else bridgerectifier = 1-cos(bridgerectifier);
				//produce either boosted or starved version
				if (mid > 0) mid = (mid*(1-temp))+(bridgerectifier*temp);
				else mid = (mid*(1-temp))-(bridgerectifier*temp);
				//blend according to density control
			}
	
			if (count < 1 || count > 2048) {count = 2048;}
			if (offset > 0)
			{
				p[count+2048] = p[count] = mid;
				mid = p[count+near]*nearLevel;
				mid += p[count+far]*farLevel;
			}
	
			if (offset < 0)
			{
				p[count+2048] = p[count] = side;
				side = p[count+near]*nearLevel;
				side += p[count+far]*farLevel;
			}
			count -= 1;
	
			inputSampleL = (mid+side) * 0.5;
			inputSampleR = (mid-side) * 0.5;
		}

		//begin SubTight section
		double subSampleL = inputSampleL * subTrim;
		double subSampleR = inputSampleR * subTrim;

		double scale = 0.5+fabs(subSampleL*0.5);
		subSampleL = (subAL+(sin(subAL-subSampleL)*scale));
		subAL = subSampleL*scale;
		scale = 0.5+fabs(subSampleR*0.5);
		subSampleR = (subAR+(sin(subAR-subSampleR)*scale));
		subAR = subSampleR*scale;
		scale = 0.5+fabs(subSampleL*0.5);
		subSampleL = (subBL+(sin(subBL-subSampleL)*scale));
		subBL = subSampleL*scale;
		scale = 0.5+fabs(subSampleR*0.5);
		subSampleR = (subBR+(sin(subBR-subSampleR)*scale));
		subBR = subSampleR*scale;
		scale = 0.5+fabs(subSampleL*0.5);
		subSampleL = (subCL+(sin(subCL-subSampleL)*scale));
		subCL = subSampleL*scale;
		scale = 0.5+fabs(subSampleR*0.5);
		subSampleR = (subCR+(sin(subCR-subSampleR)*scale));
		subCR = subSampleR*scale;
		if (subSampleL > 0.25) subSampleL = 0.25;
		if (subSampleL < -0.25) subSampleL = -0.25;
		if (subSampleR > 0.25) subSampleR = 0.25;
		if (subSampleR < -0.25) subSampleR = -0.25;
		inputSampleL += (subSampleL*16.0);
		inputSampleR += (subSampleR*16.0);
		//end SubTight section

		// Hypersonic

		temp = (inputSampleL * fixD[fix_a0]) + fixD[fix_sL1];
		fixD[fix_sL1] = (inputSampleL * fixD[fix_a1]) - (temp * fixD[fix_b1]) + fixD[fix_sL2];
		fixD[fix_sL2] = (inputSampleL * fixD[fix_a2]) - (temp * fixD[fix_b2]);
		inputSampleL = temp; //fixed biquad filtering ultrasonics
		temp = (inputSampleR * fixD[fix_a0]) + fixD[fix_sR1];
		fixD[fix_sR1] = (inputSampleR * fixD[fix_a1]) - (temp * fixD[fix_b1]) + fixD[fix_sR2];
		fixD[fix_sR2] = (inputSampleR * fixD[fix_a2]) - (temp * fixD[fix_b2]);
		inputSampleR = temp; //fixed biquad filtering ultrasonics

		//begin Console7 Channel processing
		if (inputSampleL > 1.097) inputSampleL = 1.097;
		if (inputSampleL < -1.097) inputSampleL = -1.097;
		if (inputSampleR > 1.097) inputSampleR = 1.097;
		if (inputSampleR < -1.097) inputSampleR = -1.097;
		inputSampleL = ((sin(inputSampleL*fabs(inputSampleL))/((fabs(inputSampleL) == 0.0) ?1:fabs(inputSampleL)))*0.8)+(sin(inputSampleL)*0.2);
		inputSampleR = ((sin(inputSampleR*fabs(inputSampleR))/((fabs(inputSampleR) == 0.0) ?1:fabs(inputSampleR)))*0.8)+(sin(inputSampleR)*0.2);
		//this is a version of Spiral blended 80/20 with regular Density.
		//It's blending between two different harmonics in the overtones of the algorithm

		if(flip) {
			//begin 32 bit stereo floating point dither

			int expon; frexpf((float)inputSampleL, &expon);
			fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
			inputSampleL += ((double(fpdL)-uint32_t(0x7fffffff)) * 5.5e-36l * pow(2,expon+62));
			frexpf((float)inputSampleR, &expon);
			fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;
			inputSampleR += ((double(fpdR)-uint32_t(0x7fffffff)) * 5.5e-36l * pow(2,expon+62));
			//end 32 bit stereo floating point dither

			*out1 = inputSampleL;
			*out2 = inputSampleR;

			*in1++;
			*in2++;
			*out1++;
			*out2++;
			--sampleFrames;
    }

		flip = !flip;
  }
}

void ConsoleZChannelHPHP::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
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
