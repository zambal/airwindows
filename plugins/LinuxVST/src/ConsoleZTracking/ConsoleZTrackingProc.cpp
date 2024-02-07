/* ========================================
 *  ConsoleZTracking - ConsoleZTracking.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTracking_H
#include "ConsoleZTracking.h"
#endif

void ConsoleZTracking::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames)
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	bool enableHPF = N > 0.0;
	bool enableLPF = O > 0.0;
	bool enableuLaw = P > 0.0;

	double drySampleL;
	double drySampleR;

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	double KK;
	double norm;

	// BitshiftGain

	double input_gain = 1.0;
	switch ((int)(A * 32)-16)
	{
		case -16: input_gain = 0.0000152587890625; break;
		case -15: input_gain = 0.000030517578125; break;
		case -14: input_gain = 0.00006103515625; break;
		case -13: input_gain = 0.0001220703125; break;
		case -12: input_gain = 0.000244140625; break;
		case -11: input_gain = 0.00048828125; break;
		case -10: input_gain = 0.0009765625; break;
		case -9: input_gain = 0.001953125; break;
		case -8: input_gain = 0.00390625; break;
		case -7: input_gain = 0.0078125; break;
		case -6: input_gain = 0.015625; break;
		case -5: input_gain = 0.03125; break;
		case -4: input_gain = 0.0625; break;
		case -3: input_gain = 0.125; break;
		case -2: input_gain = 0.25; break;
		case -1: input_gain = 0.5; break;
		case 0: input_gain = 1.0; break;
		case 1: input_gain = 2.0; break;
		case 2: input_gain = 4.0; break;
		case 3: input_gain = 8.0; break;
		case 4: input_gain = 16.0; break;
		case 5: input_gain = 32.0; break;
		case 6: input_gain = 64.0; break;
		case 7: input_gain = 128.0; break;
		case 8: input_gain = 256.0; break;
		case 9: input_gain = 512.0; break;
		case 10: input_gain = 1024.0; break;
		case 11: input_gain = 2048.0; break;
		case 12: input_gain = 4096.0; break;
		case 13: input_gain = 8192.0; break;
		case 14: input_gain = 16384.0; break;
		case 15: input_gain = 32768.0; break;
		case 16: input_gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// Hypersonic

	double cutoff = 25000.0 / getSampleRate();
	if (cutoff > 0.45) cutoff = 0.45; //don't crash if run at 44.1k

	fixE[fix_freq] = fixD[fix_freq] = fixC[fix_freq] = fixB[fix_freq] = fixA[fix_freq] = cutoff;

  fixA[fix_reso] = 4.46570214;
	fixB[fix_reso] = 1.51387132;
	fixC[fix_reso] = 0.93979296;
	fixD[fix_reso] = 0.70710678;
	fixE[fix_reso] = 0.59051105;

	KK = tan(M_PI * fixA[fix_freq]); //lowpass
	norm = 1.0 / (1.0 + KK / fixA[fix_reso] + KK * KK);
	fixA[fix_a0] = KK * KK * norm;
	fixA[fix_a1] = 2.0 * fixA[fix_a0];
	fixA[fix_a2] = fixA[fix_a0];
	fixA[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixA[fix_b2] = (1.0 - KK / fixA[fix_reso] + KK * KK) * norm;

	KK = tan(M_PI * fixB[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixB[fix_reso] + KK * KK);
	fixB[fix_a0] = KK * KK * norm;
	fixB[fix_a1] = 2.0 * fixB[fix_a0];
	fixB[fix_a2] = fixB[fix_a0];
	fixB[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixB[fix_b2] = (1.0 - KK / fixB[fix_reso] + KK * KK) * norm;

	KK = tan(M_PI * fixC[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixC[fix_reso] + KK * KK);
	fixC[fix_a0] = KK * KK * norm;
	fixC[fix_a1] = 2.0 * fixC[fix_a0];
	fixC[fix_a2] = fixC[fix_a0];
	fixC[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixC[fix_b2] = (1.0 - KK / fixC[fix_reso] + KK * KK) * norm;

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

	// Interstage

	double firstStage = 0.381966011250105 / overallscale;
	double iirAmount = 0.00295 / overallscale;
	double threshold = 0.381966011250105;

	// BiquadOneHalf HPF

	biquad_hpf_AL[0] = ((B*B*B*0.9999)+0.0001)*0.499;
	if (biquad_hpf_AL[0] < 0.0001) biquad_hpf_AL[0] = 0.0001;
	
    biquad_hpf_AL[1] = (C*C*C*29.99)+0.01;
	if (biquad_hpf_AL[1] < 0.0001) biquad_hpf_AL[1] = 0.0001;

	KK = tan(M_PI * biquad_hpf_AL[0]);
	norm = 1.0 / (1.0 + KK / biquad_hpf_AL[1] + KK * KK);
	biquad_hpf_AL[2] = norm;
	biquad_hpf_AL[3] = -2.0 * biquad_hpf_AL[2];
	biquad_hpf_AL[4] = biquad_hpf_AL[2];
	biquad_hpf_AL[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquad_hpf_AL[6] = (1.0 - KK / biquad_hpf_AL[1] + KK * KK) * norm;

	for (int x = 0; x < 7; x++) {biquad_hpf_AR[x] = biquad_hpf_BL[x] = biquad_hpf_BR[x] = biquad_hpf_AL[x];}

	// BiquadOneHalf LPF

	biquad_lpf_AL[0] = ((D*D*D*0.9999)+0.0001)*0.499;
	if (biquad_lpf_AL[0] < 0.0001) biquad_lpf_AL[0] = 0.0001;
	
    biquad_lpf_AL[1] = (E*E*E*29.99)+0.01;
	if (biquad_lpf_AL[1] < 0.0001) biquad_lpf_AL[1] = 0.0001;

	KK = tan(M_PI * biquad_lpf_AL[0]);
	norm = 1.0 / (1.0 + KK / biquad_lpf_AL[1] + KK * KK);
	biquad_lpf_AL[2] = KK * KK * norm;
	biquad_lpf_AL[3] = 2.0 * biquad_lpf_AL[2];
	biquad_lpf_AL[4] = biquad_lpf_AL[2];
	biquad_lpf_AL[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquad_lpf_AL[6] = (1.0 - KK / biquad_lpf_AL[1] + KK * KK) * norm;

	for (int x = 0; x < 7; x++) {biquad_lpf_AR[x] = biquad_lpf_BL[x] = biquad_lpf_BR[x] = biquad_lpf_AL[x];}

	// Tape

	double inputgain = pow(10.0,((F-0.5)*24.0)/20.0);
	double bumpgain = G*0.1;
	double HeadBumpFreq = 0.12/overallscale;
	double softness = 0.618033988749894848204586;
	double RollAmount = (1.0 - softness) / overallscale;
	//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
	//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
	biquadAL[0] = biquadBL[0] = biquadAR[0] = biquadBR[0] = 0.0072/overallscale;
	biquadAL[1] = biquadBL[1] = biquadAR[1] = biquadBR[1] = 0.0009;
	KK = tan(M_PI * biquadBR[0]);
	norm = 1.0 / (1.0 + KK / biquadBR[1] + KK * KK);
	biquadAL[2] = biquadBL[2] = biquadAR[2] = biquadBR[2] = KK / biquadBR[1] * norm;
	biquadAL[4] = biquadBL[4] = biquadAR[4] = biquadBR[4] = -biquadBR[2];
	biquadAL[5] = biquadBL[5] = biquadAR[5] = biquadBR[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquadAL[6] = biquadBL[6] = biquadAR[6] = biquadBR[6] = (1.0 - KK / biquadBR[1] + KK * KK) * norm;

	biquadCL[0] = biquadDL[0] = biquadCR[0] = biquadDR[0] = 0.032/overallscale;
	biquadCL[1] = biquadDL[1] = biquadCR[1] = biquadDR[1] = 0.0007;
	KK = tan(M_PI * biquadDR[0]);
	norm = 1.0 / (1.0 + KK / biquadDR[1] + KK * KK);
	biquadCL[2] = biquadDL[2] = biquadCR[2] = biquadDR[2] = KK / biquadDR[1] * norm;
	biquadCL[4] = biquadDL[4] = biquadCR[4] = biquadDR[4] = -biquadDR[2];
	biquadCL[5] = biquadDL[5] = biquadCR[5] = biquadDR[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquadCL[6] = biquadDL[6] = biquadCR[6] = biquadDR[6] = (1.0 - KK / biquadDR[1] + KK * KK) * norm;

	// BitshiftGain

	double trim_gain = 1.0;
	switch ((int)(H * 32)-16)
	{
		case -16: trim_gain = 0.0000152587890625; break;
		case -15: trim_gain = 0.000030517578125; break;
		case -14: trim_gain = 0.00006103515625; break;
		case -13: trim_gain = 0.0001220703125; break;
		case -12: trim_gain = 0.000244140625; break;
		case -11: trim_gain = 0.00048828125; break;
		case -10: trim_gain = 0.0009765625; break;
		case -9: trim_gain = 0.001953125; break;
		case -8: trim_gain = 0.00390625; break;
		case -7: trim_gain = 0.0078125; break;
		case -6: trim_gain = 0.015625; break;
		case -5: trim_gain = 0.03125; break;
		case -4: trim_gain = 0.0625; break;
		case -3: trim_gain = 0.125; break;
		case -2: trim_gain = 0.25; break;
		case -1: trim_gain = 0.5; break;
		case 0: trim_gain = 1.0; break;
		case 1: trim_gain = 2.0; break;
		case 2: trim_gain = 4.0; break;
		case 3: trim_gain = 8.0; break;
		case 4: trim_gain = 16.0; break;
		case 5: trim_gain = 32.0; break;
		case 6: trim_gain = 64.0; break;
		case 7: trim_gain = 128.0; break;
		case 8: trim_gain = 256.0; break;
		case 9: trim_gain = 512.0; break;
		case 10: trim_gain = 1024.0; break;
		case 11: trim_gain = 2048.0; break;
		case 12: trim_gain = 4096.0; break;
		case 13: trim_gain = 8192.0; break;
		case 14: trim_gain = 16384.0; break;
		case 15: trim_gain = 32768.0; break;
		case 16: trim_gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// Creature

	double source = 1.0-pow(1.0-I,5);
	int stages = (pow(J,2)*32.0*sqrt(overallscale))+1;
	double wet = (K*2.0)-1.0; //inv-dry-wet for highpass
	double dry = 2.0-(K*2.0);
	if (dry > 1.0) dry = 1.0; //full dry for use with inv, to 0.0 at full wet


	// ToneSlant

	double ts_inputSampleL;
	double ts_inputSampleR;
	double ts_correctionSampleL;
	double ts_correctionSampleR;
	double ts_accumulatorSampleL;
	double ts_accumulatorSampleR;
	double ts_drySampleL;
	double ts_drySampleR;
	double ts_overallscale = (L*99.0)+1.0;
	double ts_applySlant = (M*2.0)-1.0;
	
	
	ts_f[0] = 1.0 / ts_overallscale;
	//count to f(gain) which will be 0. f(0) is x1
	for (int count = 1; count < 102; count++) {
		if (count <= ts_overallscale) {
			ts_f[count] = (1.0 - (count / ts_overallscale)) / ts_overallscale;
			//recalc the filter and don't change the buffer it'll apply to
		} else {
			ts_bL[count] = 0.0; //blank the unused buffer so when we return to it, no pops
			ts_bR[count] = 0.0; //blank the unused buffer so when we return to it, no pops
		}
	}

  while (--sampleFrames >= 0)
  {

		double inputSampleL = *in1 * input_gain;
		double inputSampleR = *in2 * input_gain;
		double outSample;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// Hypersonic

		outSample = (inputSampleL * fixA[fix_a0]) + fixA[fix_sL1];
		fixA[fix_sL1] = (inputSampleL * fixA[fix_a1]) - (outSample * fixA[fix_b1]) + fixA[fix_sL2];
		fixA[fix_sL2] = (inputSampleL * fixA[fix_a2]) - (outSample * fixA[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixA[fix_a0]) + fixA[fix_sR1];
		fixA[fix_sR1] = (inputSampleR * fixA[fix_a1]) - (outSample * fixA[fix_b1]) + fixA[fix_sR2];
		fixA[fix_sR2] = (inputSampleR * fixA[fix_a2]) - (outSample * fixA[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// Tape

		double HighsSampleL = 0.0;
		double HighsSampleR = 0.0;
		double NonHighsSampleL = 0.0;
		double NonHighsSampleR = 0.0;
		double tempSample;

		if (flip)
		{
			// Interstage
			
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			iirSampleAL = (iirSampleAL * (1 - firstStage)) + (inputSampleL * firstStage); inputSampleL = iirSampleAL;
			iirSampleCL = (iirSampleCL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleCL;
			iirSampleEL = (iirSampleEL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleEL;
			inputSampleL = drySampleL - inputSampleL;
			//make highpass
			if (inputSampleL - iirSampleAL > threshold) inputSampleL = iirSampleAL + threshold;
			if (inputSampleL - iirSampleAL < -threshold) inputSampleL = iirSampleAL - threshold;
			//slew limit against lowpassed reference point

			iirSampleAR = (iirSampleAR * (1 - firstStage)) + (inputSampleR * firstStage); inputSampleR = iirSampleAR;
			iirSampleCR = (iirSampleCR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleCR;
			iirSampleER = (iirSampleER * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleER;
			inputSampleR = drySampleR - inputSampleR;
			//make highpass
			if (inputSampleR - iirSampleAR > threshold) inputSampleR = iirSampleAR + threshold;
			if (inputSampleR - iirSampleAR < -threshold) inputSampleR = iirSampleAR - threshold;
			//slew limit against lowpassed reference point

			if(enableHPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_hpf_AL[2]) + biquad_hpf_AL[7];
				biquad_hpf_AL[7] = (inputSampleL * biquad_hpf_AL[3]) - (tempSample * biquad_hpf_AL[5]) + biquad_hpf_AL[8];
				biquad_hpf_AL[8] = (inputSampleL * biquad_hpf_AL[4]) - (tempSample * biquad_hpf_AL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_hpf_AR[2]) + biquad_hpf_AR[7];
				biquad_hpf_AR[7] = (inputSampleR * biquad_hpf_AR[3]) - (tempSample * biquad_hpf_AR[5]) + biquad_hpf_AR[8];
				biquad_hpf_AR[8] = (inputSampleR * biquad_hpf_AR[4]) - (tempSample * biquad_hpf_AR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}

			if(enableLPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_lpf_AL[2]) + biquad_lpf_AL[7];
				biquad_lpf_AL[7] = (inputSampleL * biquad_lpf_AL[3]) - (tempSample * biquad_lpf_AL[5]) + biquad_lpf_AL[8];
				biquad_lpf_AL[8] = (inputSampleL * biquad_lpf_AL[4]) - (tempSample * biquad_lpf_AL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_lpf_AR[2]) + biquad_lpf_AR[7];
				biquad_lpf_AR[7] = (inputSampleR * biquad_lpf_AR[3]) - (tempSample * biquad_lpf_AR[5]) + biquad_lpf_AR[8];
				biquad_lpf_AR[8] = (inputSampleR * biquad_lpf_AR[4]) - (tempSample * biquad_lpf_AR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}
			
			// Tape
			
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			iirMidRollerAL = (iirMidRollerAL * (1.0 - RollAmount)) + (inputSampleL * RollAmount);
			iirMidRollerAR = (iirMidRollerAR * (1.0 - RollAmount)) + (inputSampleR * RollAmount);
			HighsSampleL = inputSampleL - iirMidRollerAL;
			HighsSampleR = inputSampleR - iirMidRollerAR;
			NonHighsSampleL = iirMidRollerAL;
			NonHighsSampleR = iirMidRollerAR;

			iirHeadBumpAL += (inputSampleL * 0.05);
			iirHeadBumpAR += (inputSampleR * 0.05);
			iirHeadBumpAL -= (iirHeadBumpAL * iirHeadBumpAL * iirHeadBumpAL * HeadBumpFreq);
			iirHeadBumpAR -= (iirHeadBumpAR * iirHeadBumpAR * iirHeadBumpAR * HeadBumpFreq);
			iirHeadBumpAL = sin(iirHeadBumpAL);
			iirHeadBumpAR = sin(iirHeadBumpAR);

			tempSample = (iirHeadBumpAL * biquadAL[2]) + biquadAL[7];
			biquadAL[7] = (iirHeadBumpAL * biquadAL[3]) - (tempSample * biquadAL[5]) + biquadAL[8];
			biquadAL[8] = (iirHeadBumpAL * biquadAL[4]) - (tempSample * biquadAL[6]);
			iirHeadBumpAL = tempSample; //interleaved biquad
			if (iirHeadBumpAL > 1.0) iirHeadBumpAL = 1.0;
			if (iirHeadBumpAL < -1.0) iirHeadBumpAL = -1.0;
			iirHeadBumpAL = asin(iirHeadBumpAL);

			tempSample = (iirHeadBumpAR * biquadAR[2]) + biquadAR[7];
			biquadAR[7] = (iirHeadBumpAR * biquadAR[3]) - (tempSample * biquadAR[5]) + biquadAR[8];
			biquadAR[8] = (iirHeadBumpAR * biquadAR[4]) - (tempSample * biquadAR[6]);
			iirHeadBumpAR = tempSample; //interleaved biquad
			if (iirHeadBumpAR > 1.0) iirHeadBumpAR = 1.0;
			if (iirHeadBumpAR < -1.0) iirHeadBumpAR = -1.0;
			iirHeadBumpAR = asin(iirHeadBumpAR);

			inputSampleL = sin(inputSampleL);
			tempSample = (inputSampleL * biquadCL[2]) + biquadCL[7];
			biquadCL[7] = (inputSampleL * biquadCL[3]) - (tempSample * biquadCL[5]) + biquadCL[8];
			biquadCL[8] = (inputSampleL * biquadCL[4]) - (tempSample * biquadCL[6]);
			inputSampleL = tempSample; //interleaved biquad
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
			inputSampleL = asin(inputSampleL);

			inputSampleR = sin(inputSampleR);
			tempSample = (inputSampleR * biquadCR[2]) + biquadCR[7];
			biquadCR[7] = (inputSampleR * biquadCR[3]) - (tempSample * biquadCR[5]) + biquadCR[8];
			biquadCR[8] = (inputSampleR * biquadCR[4]) - (tempSample * biquadCR[6]);
			inputSampleR = tempSample; //interleaved biquad
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
			inputSampleR = asin(inputSampleR);
		} else {
			// Interstage
			
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			iirSampleBL = (iirSampleBL * (1 - firstStage)) + (inputSampleL * firstStage); inputSampleL = iirSampleBL;
			iirSampleDL = (iirSampleDL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleDL;
			iirSampleFL = (iirSampleFL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleFL;
			inputSampleL = drySampleL - inputSampleL;
			//make highpass
			if (inputSampleL - iirSampleBL > threshold) inputSampleL = iirSampleBL + threshold;
			if (inputSampleL - iirSampleBL < -threshold) inputSampleL = iirSampleBL - threshold;
			//slew limit against lowpassed reference point

			iirSampleBR = (iirSampleBR * (1 - firstStage)) + (inputSampleR * firstStage); inputSampleR = iirSampleBR;
			iirSampleDR = (iirSampleDR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleDR;
			iirSampleFR = (iirSampleFR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleFR;
			inputSampleR = drySampleR - inputSampleR;
			//make highpass
			if (inputSampleR - iirSampleBR > threshold) inputSampleR = iirSampleBR + threshold;
			if (inputSampleR - iirSampleBR < -threshold) inputSampleR = iirSampleBR - threshold;
			//slew limit against lowpassed reference point

			if(enableHPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_hpf_BL[2]) + biquad_hpf_BL[7];
				biquad_hpf_BL[7] = (inputSampleL * biquad_hpf_BL[3]) - (tempSample * biquad_hpf_BL[5]) + biquad_hpf_BL[8];
				biquad_hpf_BL[8] = (inputSampleL * biquad_hpf_BL[4]) - (tempSample * biquad_hpf_BL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_hpf_BR[2]) + biquad_hpf_BR[7];
				biquad_hpf_BR[7] = (inputSampleR * biquad_hpf_BR[3]) - (tempSample * biquad_hpf_BR[5]) + biquad_hpf_BR[8];
				biquad_hpf_BR[8] = (inputSampleR * biquad_hpf_BR[4]) - (tempSample * biquad_hpf_BR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}

			if(enableLPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_lpf_BL[2]) + biquad_lpf_BL[7];
				biquad_lpf_BL[7] = (inputSampleL * biquad_lpf_BL[3]) - (tempSample * biquad_lpf_BL[5]) + biquad_lpf_BL[8];
				biquad_lpf_BL[8] = (inputSampleL * biquad_lpf_BL[4]) - (tempSample * biquad_lpf_BL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_lpf_BR[2]) + biquad_lpf_BR[7];
				biquad_lpf_BR[7] = (inputSampleR * biquad_lpf_BR[3]) - (tempSample * biquad_lpf_BR[5]) + biquad_lpf_BR[8];
				biquad_lpf_BR[8] = (inputSampleR * biquad_lpf_BR[4]) - (tempSample * biquad_lpf_BR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}
			
			// Tape
			
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			iirMidRollerBL = (iirMidRollerBL * (1.0 - RollAmount)) + (inputSampleL * RollAmount);
			iirMidRollerBR = (iirMidRollerBR * (1.0 - RollAmount)) + (inputSampleR * RollAmount);
			HighsSampleL = inputSampleL - iirMidRollerBL;
			HighsSampleR = inputSampleR - iirMidRollerBR;
			NonHighsSampleL = iirMidRollerBL;
			NonHighsSampleR = iirMidRollerBR;

			iirHeadBumpBL += (inputSampleL * 0.05);
			iirHeadBumpBR += (inputSampleR * 0.05);
			iirHeadBumpBL -= (iirHeadBumpBL * iirHeadBumpBL * iirHeadBumpBL * HeadBumpFreq);
			iirHeadBumpBR -= (iirHeadBumpBR * iirHeadBumpBR * iirHeadBumpBR * HeadBumpFreq);
			iirHeadBumpBL = sin(iirHeadBumpBL);
			iirHeadBumpBR = sin(iirHeadBumpBR);

			tempSample = (iirHeadBumpBL * biquadBL[2]) + biquadBL[7];
			biquadBL[7] = (iirHeadBumpBL * biquadBL[3]) - (tempSample * biquadBL[5]) + biquadBL[8];
			biquadBL[8] = (iirHeadBumpBL * biquadBL[4]) - (tempSample * biquadBL[6]);
			iirHeadBumpBL = tempSample; //interleaved biquad
			if (iirHeadBumpBL > 1.0) iirHeadBumpBL = 1.0;
			if (iirHeadBumpBL < -1.0) iirHeadBumpBL = -1.0;
			iirHeadBumpBL = asin(iirHeadBumpBL);

			tempSample = (iirHeadBumpBR * biquadBR[2]) + biquadBR[7];
			biquadBR[7] = (iirHeadBumpBR * biquadBR[3]) - (tempSample * biquadBR[5]) + biquadBR[8];
			biquadBR[8] = (iirHeadBumpBR * biquadBR[4]) - (tempSample * biquadBR[6]);
			iirHeadBumpBR = tempSample; //interleaved biquad
			if (iirHeadBumpBR > 1.0) iirHeadBumpBR = 1.0;
			if (iirHeadBumpBR < -1.0) iirHeadBumpBR = -1.0;
			iirHeadBumpBR = asin(iirHeadBumpBR);

			inputSampleL = sin(inputSampleL);
			tempSample = (inputSampleL * biquadDL[2]) + biquadDL[7];
			biquadDL[7] = (inputSampleL * biquadDL[3]) - (tempSample * biquadDL[5]) + biquadDL[8];
			biquadDL[8] = (inputSampleL * biquadDL[4]) - (tempSample * biquadDL[6]);
			inputSampleL = tempSample; //interleaved biquad
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
			inputSampleL = asin(inputSampleL);

			inputSampleR = sin(inputSampleR);
			tempSample = (inputSampleR * biquadDR[2]) + biquadDR[7];
			biquadDR[7] = (inputSampleR * biquadDR[3]) - (tempSample * biquadDR[5]) + biquadDR[8];
			biquadDR[8] = (inputSampleR * biquadDR[4]) - (tempSample * biquadDR[6]);
			inputSampleR = tempSample; //interleaved biquad
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
			inputSampleR = asin(inputSampleR);
		}
		flip = !flip;

		// Hypersonic

		outSample = (inputSampleL * fixB[fix_a0]) + fixB[fix_sL1];
		fixB[fix_sL1] = (inputSampleL * fixB[fix_a1]) - (outSample * fixB[fix_b1]) + fixB[fix_sL2];
		fixB[fix_sL2] = (inputSampleL * fixB[fix_a2]) - (outSample * fixB[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixB[fix_a0]) + fixB[fix_sR1];
		fixB[fix_sR1] = (inputSampleR * fixB[fix_a1]) - (outSample * fixB[fix_b1]) + fixB[fix_sR2];
		fixB[fix_sR2] = (inputSampleR * fixB[fix_a2]) - (outSample * fixB[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		double groundSampleL = drySampleL - inputSampleL; //set up UnBox
		double groundSampleR = drySampleR - inputSampleR; //set up UnBox

		if (inputgain != 1.0) {
			inputSampleL *= inputgain;
			inputSampleR *= inputgain;
		} //gain boost inside UnBox: do not boost fringe audio

		double applySoften = fabs(HighsSampleL)*1.57079633;
		if (applySoften > 1.57079633) applySoften = 1.57079633;
		applySoften = 1-cos(applySoften);
		if (HighsSampleL > 0) inputSampleL -= applySoften;
		if (HighsSampleL < 0) inputSampleL += applySoften;
		//apply Soften depending on polarity
		applySoften = fabs(HighsSampleR)*1.57079633;
		if (applySoften > 1.57079633) applySoften = 1.57079633;
		applySoften = 1-cos(applySoften);
		if (HighsSampleR > 0) inputSampleR -= applySoften;
		if (HighsSampleR < 0) inputSampleR += applySoften;
		//apply Soften depending on polarity

		if (inputSampleL > 1.2533141373155) inputSampleL = 1.2533141373155;
		if (inputSampleL < -1.2533141373155) inputSampleL = -1.2533141373155;
		//clip to 1.2533141373155 to reach maximum output
		inputSampleL = sin(inputSampleL * fabs(inputSampleL)) / ((fabs(inputSampleL) == 0.0) ?1:fabs(inputSampleL));
		//Spiral, for cleanest most optimal tape effect
		if (inputSampleR > 1.2533141373155) inputSampleR = 1.2533141373155;
		if (inputSampleR < -1.2533141373155) inputSampleR = -1.2533141373155;
		//clip to 1.2533141373155 to reach maximum output
		inputSampleR = sin(inputSampleR * fabs(inputSampleR)) / ((fabs(inputSampleR) == 0.0) ?1:fabs(inputSampleR));
		//Spiral, for cleanest most optimal tape effect

		double suppress = (1.0-fabs(inputSampleL)) * 0.00013;
		if (iirHeadBumpAL > suppress) iirHeadBumpAL -= suppress;
		if (iirHeadBumpAL < -suppress) iirHeadBumpAL += suppress;
		if (iirHeadBumpBL > suppress) iirHeadBumpBL -= suppress;
		if (iirHeadBumpBL < -suppress) iirHeadBumpBL += suppress;
		//restrain resonant quality of head bump algorithm
		suppress = (1.0-fabs(inputSampleR)) * 0.00013;
		if (iirHeadBumpAR > suppress) iirHeadBumpAR -= suppress;
		if (iirHeadBumpAR < -suppress) iirHeadBumpAR += suppress;
		if (iirHeadBumpBR > suppress) iirHeadBumpBR -= suppress;
		if (iirHeadBumpBR < -suppress) iirHeadBumpBR += suppress;
		//restrain resonant quality of head bump algorithm

		inputSampleL += groundSampleL; //apply UnBox processing
		inputSampleR += groundSampleR; //apply UnBox processing

		inputSampleL += ((iirHeadBumpAL + iirHeadBumpBL) * bumpgain);//and head bump
		inputSampleR += ((iirHeadBumpAR + iirHeadBumpBR) * bumpgain);//and head bump

		// Hypersonic

		outSample = (inputSampleL * fixC[fix_a0]) + fixC[fix_sL1];
		fixC[fix_sL1] = (inputSampleL * fixC[fix_a1]) - (outSample * fixC[fix_b1]) + fixC[fix_sL2];
		fixC[fix_sL2] = (inputSampleL * fixC[fix_a2]) - (outSample * fixC[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixC[fix_a0]) + fixC[fix_sR1];
		fixC[fix_sR1] = (inputSampleR * fixC[fix_a1]) - (outSample * fixC[fix_b1]) + fixC[fix_sR2];
		fixC[fix_sR2] = (inputSampleR * fixC[fix_a2]) - (outSample * fixC[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// uLawEncode

		if(enableuLaw) {
			static int noisesourceL = 0;
			static int noisesourceR = 850010;
			int residue;
			double applyresidue;
		
			noisesourceL = noisesourceL % 1700021; noisesourceL++;
			residue = noisesourceL * noisesourceL;
			residue = residue % 170003; residue *= residue;
			residue = residue % 17011; residue *= residue;
			residue = residue % 1709; residue *= residue;
			residue = residue % 173; residue *= residue;
			residue = residue % 17;
			applyresidue = residue;
			applyresidue *= 0.00000001;
			applyresidue *= 0.00000001;
			inputSampleL += applyresidue;
			if (inputSampleL<1.2e-38 && -inputSampleL<1.2e-38) {
				inputSampleL -= applyresidue;
			}
		
			noisesourceR = noisesourceR % 1700021; noisesourceR++;
			residue = noisesourceR * noisesourceR;
			residue = residue % 170003; residue *= residue;
			residue = residue % 17011; residue *= residue;
			residue = residue % 1709; residue *= residue;
			residue = residue % 173; residue *= residue;
			residue = residue % 17;
			applyresidue = residue;
			applyresidue *= 0.00000001;
			applyresidue *= 0.00000001;
			inputSampleR += applyresidue;
			if (inputSampleR<1.2e-38 && -inputSampleR<1.2e-38) {
				inputSampleR -= applyresidue;
			}
			//for live air, we always apply the dither noise. Then, if our result is 
			//effectively digital black, we'll subtract it auLawEncode. We want a 'air' hiss

			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
		
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
		
			if (inputSampleL > 0) inputSampleL = log(1.0+(255*fabs(inputSampleL))) / log(256);
			if (inputSampleL < 0) inputSampleL = -log(1.0+(255*fabs(inputSampleL))) / log(256);
		
			if (inputSampleR > 0) inputSampleR = log(1.0+(255*fabs(inputSampleR))) / log(256);
			if (inputSampleR < 0) inputSampleR = -log(1.0+(255*fabs(inputSampleR))) / log(256);
		}
	
		// Creature

		if(wet != 0) {
			double cr_drySampleL = inputSampleL;
			double cr_drySampleR = inputSampleR;
		
			for (int x = 0; x < stages; x++) {
				inputSampleL = (slewL[x]+(sin(slewL[x]-inputSampleL)*0.5))*source;
				slewL[x] = inputSampleL*0.5;
				inputSampleR = (slewR[x]+(sin(slewR[x]-inputSampleR)*0.5))*source;
				slewR[x] = inputSampleR*0.5;
			}
			if (stages % 2 > 0) {
				inputSampleL = -inputSampleL;
				inputSampleR = -inputSampleR;
			}
		
			inputSampleL *= wet;
			inputSampleR *= wet;
			cr_drySampleL *= dry;
			cr_drySampleR *= dry;
			inputSampleL += cr_drySampleL;
			inputSampleR += cr_drySampleR;
		}

		// uLawDecode

		if(enableuLaw) {
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
		
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
		
			if (inputSampleL > 0) inputSampleL = (pow(256,fabs(inputSampleL))-1.0) / 255;
			if (inputSampleL < 0) inputSampleL = -(pow(256,fabs(inputSampleL))-1.0) / 255;
		
			if (inputSampleR > 0) inputSampleR = (pow(256,fabs(inputSampleR))-1.0) / 255;
			if (inputSampleR < 0) inputSampleR = -(pow(256,fabs(inputSampleR))-1.0) / 255;
		}

		// BitshiftGain

		inputSampleL *= trim_gain;
		inputSampleR *= trim_gain;

		// Hypersonic

		outSample = (inputSampleL * fixD[fix_a0]) + fixD[fix_sL1];
		fixD[fix_sL1] = (inputSampleL * fixD[fix_a1]) - (outSample * fixD[fix_b1]) + fixD[fix_sL2];
		fixD[fix_sL2] = (inputSampleL * fixD[fix_a2]) - (outSample * fixD[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixD[fix_a0]) + fixD[fix_sR1];
		fixD[fix_sR1] = (inputSampleR * fixD[fix_a1]) - (outSample * fixD[fix_b1]) + fixD[fix_sR2];
		fixD[fix_sR2] = (inputSampleR * fixD[fix_a2]) - (outSample * fixD[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// ToneSlant
		
		for (int count = ts_overallscale; count >= 0; count--) {
			ts_bL[count+1] = ts_bL[count];
			ts_bR[count+1] = ts_bR[count];
		}
		
		ts_bL[0] = ts_accumulatorSampleL = ts_drySampleL = inputSampleL;
		ts_bR[0] = ts_accumulatorSampleR = ts_drySampleR = inputSampleR;
		
		ts_accumulatorSampleL *= ts_f[0];
		ts_accumulatorSampleR *= ts_f[0];

		for (int count = 1; count < ts_overallscale; count++) {
			ts_accumulatorSampleL += (ts_bL[count] * ts_f[count]);
			ts_accumulatorSampleR += (ts_bR[count] * ts_f[count]);
		}
		
		ts_correctionSampleL = inputSampleL - (ts_accumulatorSampleL*2.0);
		ts_correctionSampleR = inputSampleR - (ts_accumulatorSampleR*2.0);
		//we're gonna apply the total effect of all these calculations as a single subtract
		
		inputSampleL += (ts_correctionSampleL * ts_applySlant);
		inputSampleR += (ts_correctionSampleR * ts_applySlant);
		//our one math operation on the input data coming in

		// ADClip
		
		if (lastSampleL >= 0.5)
		{
			if (inputSampleL < 0.5) lastSampleL = ((0.5*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = 0.5;
		}

		if (lastSampleL <= -0.5)
		{
			if (inputSampleL > -0.5) lastSampleL = ((-0.5*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = -0.5;
		}

		if (inputSampleL > 0.5)
		{
			if (lastSampleL < 0.5) inputSampleL = ((0.5*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = 0.5;
		}

		if (inputSampleL < -0.5)
		{
			if (lastSampleL > -0.5) inputSampleL = ((-0.5*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = -0.5;
		}
		lastSampleL = inputSampleL; //end ADClip L


		if (lastSampleR >= 0.5)
		{
			if (inputSampleR < 0.5) lastSampleR = ((0.5*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = 0.5;
		}

		if (lastSampleR <= -0.5)
		{
			if (inputSampleR > -0.5) lastSampleR = ((-0.5*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = -0.5;
		}

		if (inputSampleR > 0.5)
		{
			if (lastSampleR < 0.5) inputSampleR = ((0.5*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = 0.5;
		}

		if (inputSampleR < -0.5)
		{
			if (lastSampleR > -0.5) inputSampleR = ((-0.5*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = -0.5;
		}
		lastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.5) inputSampleL = 0.5;
		if (inputSampleL < -0.5) inputSampleL = -0.5;
		//final iron bar
		if (inputSampleR > 0.5) inputSampleR = 0.5;
		if (inputSampleR < -0.5) inputSampleR = -0.5;
		//final iron bar

		// Hypersonic

		outSample = (inputSampleL * fixE[fix_a0]) + fixE[fix_sL1];
		fixE[fix_sL1] = (inputSampleL * fixE[fix_a1]) - (outSample * fixE[fix_b1]) + fixE[fix_sL2];
		fixE[fix_sL2] = (inputSampleL * fixE[fix_a2]) - (outSample * fixE[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixE[fix_a0]) + fixE[fix_sR1];
		fixE[fix_sR1] = (inputSampleR * fixE[fix_a1]) - (outSample * fixE[fix_b1]) + fixE[fix_sR2];
		fixE[fix_sR2] = (inputSampleR * fixE[fix_a2]) - (outSample * fixE[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics
		
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
    }
}

void ConsoleZTracking::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
{
  double* in1  =  inputs[0];
  double* in2  =  inputs[1];
  double* out1 = outputs[0];
  double* out2 = outputs[1];

	bool enableHPF = N > 0.0;
	bool enableLPF = O > 0.0;
	bool enableuLaw = P > 0.0;

	double drySampleL;
	double drySampleR;

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	double KK;
	double norm;

	// BitshiftGain

	double input_gain = 1.0;
	switch ((int)(A * 32)-16)
	{
		case -16: input_gain = 0.0000152587890625; break;
		case -15: input_gain = 0.000030517578125; break;
		case -14: input_gain = 0.00006103515625; break;
		case -13: input_gain = 0.0001220703125; break;
		case -12: input_gain = 0.000244140625; break;
		case -11: input_gain = 0.00048828125; break;
		case -10: input_gain = 0.0009765625; break;
		case -9: input_gain = 0.001953125; break;
		case -8: input_gain = 0.00390625; break;
		case -7: input_gain = 0.0078125; break;
		case -6: input_gain = 0.015625; break;
		case -5: input_gain = 0.03125; break;
		case -4: input_gain = 0.0625; break;
		case -3: input_gain = 0.125; break;
		case -2: input_gain = 0.25; break;
		case -1: input_gain = 0.5; break;
		case 0: input_gain = 1.0; break;
		case 1: input_gain = 2.0; break;
		case 2: input_gain = 4.0; break;
		case 3: input_gain = 8.0; break;
		case 4: input_gain = 16.0; break;
		case 5: input_gain = 32.0; break;
		case 6: input_gain = 64.0; break;
		case 7: input_gain = 128.0; break;
		case 8: input_gain = 256.0; break;
		case 9: input_gain = 512.0; break;
		case 10: input_gain = 1024.0; break;
		case 11: input_gain = 2048.0; break;
		case 12: input_gain = 4096.0; break;
		case 13: input_gain = 8192.0; break;
		case 14: input_gain = 16384.0; break;
		case 15: input_gain = 32768.0; break;
		case 16: input_gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// Hypersonic

	double cutoff = 25000.0 / getSampleRate();
	if (cutoff > 0.45) cutoff = 0.45; //don't crash if run at 44.1k

	fixE[fix_freq] = fixD[fix_freq] = fixC[fix_freq] = fixB[fix_freq] = fixA[fix_freq] = cutoff;

  fixA[fix_reso] = 4.46570214;
	fixB[fix_reso] = 1.51387132;
	fixC[fix_reso] = 0.93979296;
	fixD[fix_reso] = 0.70710678;
	fixE[fix_reso] = 0.59051105;

	KK = tan(M_PI * fixA[fix_freq]); //lowpass
	norm = 1.0 / (1.0 + KK / fixA[fix_reso] + KK * KK);
	fixA[fix_a0] = KK * KK * norm;
	fixA[fix_a1] = 2.0 * fixA[fix_a0];
	fixA[fix_a2] = fixA[fix_a0];
	fixA[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixA[fix_b2] = (1.0 - KK / fixA[fix_reso] + KK * KK) * norm;

	KK = tan(M_PI * fixB[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixB[fix_reso] + KK * KK);
	fixB[fix_a0] = KK * KK * norm;
	fixB[fix_a1] = 2.0 * fixB[fix_a0];
	fixB[fix_a2] = fixB[fix_a0];
	fixB[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixB[fix_b2] = (1.0 - KK / fixB[fix_reso] + KK * KK) * norm;

	KK = tan(M_PI * fixC[fix_freq]);
	norm = 1.0 / (1.0 + KK / fixC[fix_reso] + KK * KK);
	fixC[fix_a0] = KK * KK * norm;
	fixC[fix_a1] = 2.0 * fixC[fix_a0];
	fixC[fix_a2] = fixC[fix_a0];
	fixC[fix_b1] = 2.0 * (KK * KK - 1.0) * norm;
	fixC[fix_b2] = (1.0 - KK / fixC[fix_reso] + KK * KK) * norm;

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

	// Interstage

	double firstStage = 0.381966011250105 / overallscale;
	double iirAmount = 0.00295 / overallscale;
	double threshold = 0.381966011250105;

	// BiquadOneHalf HPF

	biquad_hpf_AL[0] = ((B*B*B*0.9999)+0.0001)*0.499;
	if (biquad_hpf_AL[0] < 0.0001) biquad_hpf_AL[0] = 0.0001;
	
    biquad_hpf_AL[1] = (C*C*C*29.99)+0.01;
	if (biquad_hpf_AL[1] < 0.0001) biquad_hpf_AL[1] = 0.0001;

	KK = tan(M_PI * biquad_hpf_AL[0]);
	norm = 1.0 / (1.0 + KK / biquad_hpf_AL[1] + KK * KK);
	biquad_hpf_AL[2] = norm;
	biquad_hpf_AL[3] = -2.0 * biquad_hpf_AL[2];
	biquad_hpf_AL[4] = biquad_hpf_AL[2];
	biquad_hpf_AL[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquad_hpf_AL[6] = (1.0 - KK / biquad_hpf_AL[1] + KK * KK) * norm;

	for (int x = 0; x < 7; x++) {biquad_hpf_AR[x] = biquad_hpf_BL[x] = biquad_hpf_BR[x] = biquad_hpf_AL[x];}

	// BiquadOneHalf LPF

	biquad_lpf_AL[0] = ((D*D*D*0.9999)+0.0001)*0.499;
	if (biquad_lpf_AL[0] < 0.0001) biquad_lpf_AL[0] = 0.0001;
	
    biquad_lpf_AL[1] = (E*E*E*29.99)+0.01;
	if (biquad_lpf_AL[1] < 0.0001) biquad_lpf_AL[1] = 0.0001;

	KK = tan(M_PI * biquad_lpf_AL[0]);
	norm = 1.0 / (1.0 + KK / biquad_lpf_AL[1] + KK * KK);
	biquad_lpf_AL[2] = KK * KK * norm;
	biquad_lpf_AL[3] = 2.0 * biquad_lpf_AL[2];
	biquad_lpf_AL[4] = biquad_lpf_AL[2];
	biquad_lpf_AL[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquad_lpf_AL[6] = (1.0 - KK / biquad_lpf_AL[1] + KK * KK) * norm;

	for (int x = 0; x < 7; x++) {biquad_lpf_AR[x] = biquad_lpf_BL[x] = biquad_lpf_BR[x] = biquad_lpf_AL[x];}

	// Tape

	double inputgain = pow(10.0,((F-0.5)*24.0)/20.0);
	double bumpgain = G*0.1;
	double HeadBumpFreq = 0.12/overallscale;
	double softness = 0.618033988749894848204586;
	double RollAmount = (1.0 - softness) / overallscale;
	//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
	//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
	biquadAL[0] = biquadBL[0] = biquadAR[0] = biquadBR[0] = 0.0072/overallscale;
	biquadAL[1] = biquadBL[1] = biquadAR[1] = biquadBR[1] = 0.0009;
	KK = tan(M_PI * biquadBR[0]);
	norm = 1.0 / (1.0 + KK / biquadBR[1] + KK * KK);
	biquadAL[2] = biquadBL[2] = biquadAR[2] = biquadBR[2] = KK / biquadBR[1] * norm;
	biquadAL[4] = biquadBL[4] = biquadAR[4] = biquadBR[4] = -biquadBR[2];
	biquadAL[5] = biquadBL[5] = biquadAR[5] = biquadBR[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquadAL[6] = biquadBL[6] = biquadAR[6] = biquadBR[6] = (1.0 - KK / biquadBR[1] + KK * KK) * norm;

	biquadCL[0] = biquadDL[0] = biquadCR[0] = biquadDR[0] = 0.032/overallscale;
	biquadCL[1] = biquadDL[1] = biquadCR[1] = biquadDR[1] = 0.0007;
	KK = tan(M_PI * biquadDR[0]);
	norm = 1.0 / (1.0 + KK / biquadDR[1] + KK * KK);
	biquadCL[2] = biquadDL[2] = biquadCR[2] = biquadDR[2] = KK / biquadDR[1] * norm;
	biquadCL[4] = biquadDL[4] = biquadCR[4] = biquadDR[4] = -biquadDR[2];
	biquadCL[5] = biquadDL[5] = biquadCR[5] = biquadDR[5] = 2.0 * (KK * KK - 1.0) * norm;
	biquadCL[6] = biquadDL[6] = biquadCR[6] = biquadDR[6] = (1.0 - KK / biquadDR[1] + KK * KK) * norm;

	// BitshiftGain

	double trim_gain = 1.0;
	switch ((int)(H * 32)-16)
	{
		case -16: trim_gain = 0.0000152587890625; break;
		case -15: trim_gain = 0.000030517578125; break;
		case -14: trim_gain = 0.00006103515625; break;
		case -13: trim_gain = 0.0001220703125; break;
		case -12: trim_gain = 0.000244140625; break;
		case -11: trim_gain = 0.00048828125; break;
		case -10: trim_gain = 0.0009765625; break;
		case -9: trim_gain = 0.001953125; break;
		case -8: trim_gain = 0.00390625; break;
		case -7: trim_gain = 0.0078125; break;
		case -6: trim_gain = 0.015625; break;
		case -5: trim_gain = 0.03125; break;
		case -4: trim_gain = 0.0625; break;
		case -3: trim_gain = 0.125; break;
		case -2: trim_gain = 0.25; break;
		case -1: trim_gain = 0.5; break;
		case 0: trim_gain = 1.0; break;
		case 1: trim_gain = 2.0; break;
		case 2: trim_gain = 4.0; break;
		case 3: trim_gain = 8.0; break;
		case 4: trim_gain = 16.0; break;
		case 5: trim_gain = 32.0; break;
		case 6: trim_gain = 64.0; break;
		case 7: trim_gain = 128.0; break;
		case 8: trim_gain = 256.0; break;
		case 9: trim_gain = 512.0; break;
		case 10: trim_gain = 1024.0; break;
		case 11: trim_gain = 2048.0; break;
		case 12: trim_gain = 4096.0; break;
		case 13: trim_gain = 8192.0; break;
		case 14: trim_gain = 16384.0; break;
		case 15: trim_gain = 32768.0; break;
		case 16: trim_gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// Creature

	double source = 1.0-pow(1.0-I,5);
	int stages = (pow(J,2)*32.0*sqrt(overallscale))+1;
	double wet = (K*2.0)-1.0; //inv-dry-wet for highpass
	double dry = 2.0-(K*2.0);
	if (dry > 1.0) dry = 1.0; //full dry for use with inv, to 0.0 at full wet


	// ToneSlant

	double ts_inputSampleL;
	double ts_inputSampleR;
	double ts_correctionSampleL;
	double ts_correctionSampleR;
	double ts_accumulatorSampleL;
	double ts_accumulatorSampleR;
	double ts_drySampleL;
	double ts_drySampleR;
	double ts_overallscale = (L*99.0)+1.0;
	double ts_applySlant = (M*2.0)-1.0;
	
	
	ts_f[0] = 1.0 / ts_overallscale;
	//count to f(gain) which will be 0. f(0) is x1
	for (int count = 1; count < 102; count++) {
		if (count <= ts_overallscale) {
			ts_f[count] = (1.0 - (count / ts_overallscale)) / ts_overallscale;
			//recalc the filter and don't change the buffer it'll apply to
		} else {
			ts_bL[count] = 0.0; //blank the unused buffer so when we return to it, no pops
			ts_bR[count] = 0.0; //blank the unused buffer so when we return to it, no pops
		}
	}

  while (--sampleFrames >= 0)
  {

		double inputSampleL = *in1 * input_gain;
		double inputSampleR = *in2 * input_gain;
		double outSample;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// Hypersonic

		outSample = (inputSampleL * fixA[fix_a0]) + fixA[fix_sL1];
		fixA[fix_sL1] = (inputSampleL * fixA[fix_a1]) - (outSample * fixA[fix_b1]) + fixA[fix_sL2];
		fixA[fix_sL2] = (inputSampleL * fixA[fix_a2]) - (outSample * fixA[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixA[fix_a0]) + fixA[fix_sR1];
		fixA[fix_sR1] = (inputSampleR * fixA[fix_a1]) - (outSample * fixA[fix_b1]) + fixA[fix_sR2];
		fixA[fix_sR2] = (inputSampleR * fixA[fix_a2]) - (outSample * fixA[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// Tape

		drySampleL = inputSampleL;
		drySampleR = inputSampleR;

		double HighsSampleL = 0.0;
		double HighsSampleR = 0.0;
		double NonHighsSampleL = 0.0;
		double NonHighsSampleR = 0.0;
		double tempSample;

		if (flip)
		{
			// Interstage
			
			iirSampleAL = (iirSampleAL * (1 - firstStage)) + (inputSampleL * firstStage); inputSampleL = iirSampleAL;
			iirSampleCL = (iirSampleCL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleCL;
			iirSampleEL = (iirSampleEL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleEL;
			inputSampleL = drySampleL - inputSampleL;
			//make highpass
			if (inputSampleL - iirSampleAL > threshold) inputSampleL = iirSampleAL + threshold;
			if (inputSampleL - iirSampleAL < -threshold) inputSampleL = iirSampleAL - threshold;
			//slew limit against lowpassed reference point

			iirSampleAR = (iirSampleAR * (1 - firstStage)) + (inputSampleR * firstStage); inputSampleR = iirSampleAR;
			iirSampleCR = (iirSampleCR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleCR;
			iirSampleER = (iirSampleER * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleER;
			inputSampleR = drySampleR - inputSampleR;
			//make highpass
			if (inputSampleR - iirSampleAR > threshold) inputSampleR = iirSampleAR + threshold;
			if (inputSampleR - iirSampleAR < -threshold) inputSampleR = iirSampleAR - threshold;
			//slew limit against lowpassed reference point

			if(enableHPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_hpf_AL[2]) + biquad_hpf_AL[7];
				biquad_hpf_AL[7] = (inputSampleL * biquad_hpf_AL[3]) - (tempSample * biquad_hpf_AL[5]) + biquad_hpf_AL[8];
				biquad_hpf_AL[8] = (inputSampleL * biquad_hpf_AL[4]) - (tempSample * biquad_hpf_AL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_hpf_AR[2]) + biquad_hpf_AR[7];
				biquad_hpf_AR[7] = (inputSampleR * biquad_hpf_AR[3]) - (tempSample * biquad_hpf_AR[5]) + biquad_hpf_AR[8];
				biquad_hpf_AR[8] = (inputSampleR * biquad_hpf_AR[4]) - (tempSample * biquad_hpf_AR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}

			if(enableLPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_lpf_AL[2]) + biquad_lpf_AL[7];
				biquad_lpf_AL[7] = (inputSampleL * biquad_lpf_AL[3]) - (tempSample * biquad_lpf_AL[5]) + biquad_lpf_AL[8];
				biquad_lpf_AL[8] = (inputSampleL * biquad_lpf_AL[4]) - (tempSample * biquad_lpf_AL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_lpf_AR[2]) + biquad_lpf_AR[7];
				biquad_lpf_AR[7] = (inputSampleR * biquad_lpf_AR[3]) - (tempSample * biquad_lpf_AR[5]) + biquad_lpf_AR[8];
				biquad_lpf_AR[8] = (inputSampleR * biquad_lpf_AR[4]) - (tempSample * biquad_lpf_AR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}
			
			// Tape
			
			iirMidRollerAL = (iirMidRollerAL * (1.0 - RollAmount)) + (inputSampleL * RollAmount);
			iirMidRollerAR = (iirMidRollerAR * (1.0 - RollAmount)) + (inputSampleR * RollAmount);
			HighsSampleL = inputSampleL - iirMidRollerAL;
			HighsSampleR = inputSampleR - iirMidRollerAR;
			NonHighsSampleL = iirMidRollerAL;
			NonHighsSampleR = iirMidRollerAR;

			iirHeadBumpAL += (inputSampleL * 0.05);
			iirHeadBumpAR += (inputSampleR * 0.05);
			iirHeadBumpAL -= (iirHeadBumpAL * iirHeadBumpAL * iirHeadBumpAL * HeadBumpFreq);
			iirHeadBumpAR -= (iirHeadBumpAR * iirHeadBumpAR * iirHeadBumpAR * HeadBumpFreq);
			iirHeadBumpAL = sin(iirHeadBumpAL);
			iirHeadBumpAR = sin(iirHeadBumpAR);

			tempSample = (iirHeadBumpAL * biquadAL[2]) + biquadAL[7];
			biquadAL[7] = (iirHeadBumpAL * biquadAL[3]) - (tempSample * biquadAL[5]) + biquadAL[8];
			biquadAL[8] = (iirHeadBumpAL * biquadAL[4]) - (tempSample * biquadAL[6]);
			iirHeadBumpAL = tempSample; //interleaved biquad
			if (iirHeadBumpAL > 1.0) iirHeadBumpAL = 1.0;
			if (iirHeadBumpAL < -1.0) iirHeadBumpAL = -1.0;
			iirHeadBumpAL = asin(iirHeadBumpAL);

			tempSample = (iirHeadBumpAR * biquadAR[2]) + biquadAR[7];
			biquadAR[7] = (iirHeadBumpAR * biquadAR[3]) - (tempSample * biquadAR[5]) + biquadAR[8];
			biquadAR[8] = (iirHeadBumpAR * biquadAR[4]) - (tempSample * biquadAR[6]);
			iirHeadBumpAR = tempSample; //interleaved biquad
			if (iirHeadBumpAR > 1.0) iirHeadBumpAR = 1.0;
			if (iirHeadBumpAR < -1.0) iirHeadBumpAR = -1.0;
			iirHeadBumpAR = asin(iirHeadBumpAR);

			inputSampleL = sin(inputSampleL);
			tempSample = (inputSampleL * biquadCL[2]) + biquadCL[7];
			biquadCL[7] = (inputSampleL * biquadCL[3]) - (tempSample * biquadCL[5]) + biquadCL[8];
			biquadCL[8] = (inputSampleL * biquadCL[4]) - (tempSample * biquadCL[6]);
			inputSampleL = tempSample; //interleaved biquad
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
			inputSampleL = asin(inputSampleL);

			inputSampleR = sin(inputSampleR);
			tempSample = (inputSampleR * biquadCR[2]) + biquadCR[7];
			biquadCR[7] = (inputSampleR * biquadCR[3]) - (tempSample * biquadCR[5]) + biquadCR[8];
			biquadCR[8] = (inputSampleR * biquadCR[4]) - (tempSample * biquadCR[6]);
			inputSampleR = tempSample; //interleaved biquad
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
			inputSampleR = asin(inputSampleR);
		} else {
			// Interstage
			
			iirSampleBL = (iirSampleBL * (1 - firstStage)) + (inputSampleL * firstStage); inputSampleL = iirSampleBL;
			iirSampleDL = (iirSampleDL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleDL;
			iirSampleFL = (iirSampleFL * (1 - iirAmount)) + (inputSampleL * iirAmount); inputSampleL = iirSampleFL;
			inputSampleL = drySampleL - inputSampleL;
			//make highpass
			if (inputSampleL - iirSampleBL > threshold) inputSampleL = iirSampleBL + threshold;
			if (inputSampleL - iirSampleBL < -threshold) inputSampleL = iirSampleBL - threshold;
			//slew limit against lowpassed reference point

			iirSampleBR = (iirSampleBR * (1 - firstStage)) + (inputSampleR * firstStage); inputSampleR = iirSampleBR;
			iirSampleDR = (iirSampleDR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleDR;
			iirSampleFR = (iirSampleFR * (1 - iirAmount)) + (inputSampleR * iirAmount); inputSampleR = iirSampleFR;
			inputSampleR = drySampleR - inputSampleR;
			//make highpass
			if (inputSampleR - iirSampleBR > threshold) inputSampleR = iirSampleBR + threshold;
			if (inputSampleR - iirSampleBR < -threshold) inputSampleR = iirSampleBR - threshold;
			//slew limit against lowpassed reference point

			if(enableHPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_hpf_BL[2]) + biquad_hpf_BL[7];
				biquad_hpf_BL[7] = (inputSampleL * biquad_hpf_BL[3]) - (tempSample * biquad_hpf_BL[5]) + biquad_hpf_BL[8];
				biquad_hpf_BL[8] = (inputSampleL * biquad_hpf_BL[4]) - (tempSample * biquad_hpf_BL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_hpf_BR[2]) + biquad_hpf_BR[7];
				biquad_hpf_BR[7] = (inputSampleR * biquad_hpf_BR[3]) - (tempSample * biquad_hpf_BR[5]) + biquad_hpf_BR[8];
				biquad_hpf_BR[8] = (inputSampleR * biquad_hpf_BR[4]) - (tempSample * biquad_hpf_BR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}

			if(enableLPF) {
				//inputSampleL = sin(inputSampleL);
				inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
				tempSample = (inputSampleL * biquad_lpf_BL[2]) + biquad_lpf_BL[7];
				biquad_lpf_BL[7] = (inputSampleL * biquad_lpf_BL[3]) - (tempSample * biquad_lpf_BL[5]) + biquad_lpf_BL[8];
				biquad_lpf_BL[8] = (inputSampleL * biquad_lpf_BL[4]) - (tempSample * biquad_lpf_BL[6]);
				inputSampleL = tempSample;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				//inputSampleL = asin(inputSampleL);
				inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);

				//inputSampleR = sin(inputSampleR);
				inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));
				tempSample = (inputSampleR * biquad_lpf_BR[2]) + biquad_lpf_BR[7];
				biquad_lpf_BR[7] = (inputSampleR * biquad_lpf_BR[3]) - (tempSample * biquad_lpf_BR[5]) + biquad_lpf_BR[8];
				biquad_lpf_BR[8] = (inputSampleR * biquad_lpf_BR[4]) - (tempSample * biquad_lpf_BR[6]);
				inputSampleR = tempSample;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				//inputSampleR = asin(inputSampleR);
				inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);
			}
			
			// Tape
			
			iirMidRollerBL = (iirMidRollerBL * (1.0 - RollAmount)) + (inputSampleL * RollAmount);
			iirMidRollerBR = (iirMidRollerBR * (1.0 - RollAmount)) + (inputSampleR * RollAmount);
			HighsSampleL = inputSampleL - iirMidRollerBL;
			HighsSampleR = inputSampleR - iirMidRollerBR;
			NonHighsSampleL = iirMidRollerBL;
			NonHighsSampleR = iirMidRollerBR;

			iirHeadBumpBL += (inputSampleL * 0.05);
			iirHeadBumpBR += (inputSampleR * 0.05);
			iirHeadBumpBL -= (iirHeadBumpBL * iirHeadBumpBL * iirHeadBumpBL * HeadBumpFreq);
			iirHeadBumpBR -= (iirHeadBumpBR * iirHeadBumpBR * iirHeadBumpBR * HeadBumpFreq);
			iirHeadBumpBL = sin(iirHeadBumpBL);
			iirHeadBumpBR = sin(iirHeadBumpBR);

			tempSample = (iirHeadBumpBL * biquadBL[2]) + biquadBL[7];
			biquadBL[7] = (iirHeadBumpBL * biquadBL[3]) - (tempSample * biquadBL[5]) + biquadBL[8];
			biquadBL[8] = (iirHeadBumpBL * biquadBL[4]) - (tempSample * biquadBL[6]);
			iirHeadBumpBL = tempSample; //interleaved biquad
			if (iirHeadBumpBL > 1.0) iirHeadBumpBL = 1.0;
			if (iirHeadBumpBL < -1.0) iirHeadBumpBL = -1.0;
			iirHeadBumpBL = asin(iirHeadBumpBL);

			tempSample = (iirHeadBumpBR * biquadBR[2]) + biquadBR[7];
			biquadBR[7] = (iirHeadBumpBR * biquadBR[3]) - (tempSample * biquadBR[5]) + biquadBR[8];
			biquadBR[8] = (iirHeadBumpBR * biquadBR[4]) - (tempSample * biquadBR[6]);
			iirHeadBumpBR = tempSample; //interleaved biquad
			if (iirHeadBumpBR > 1.0) iirHeadBumpBR = 1.0;
			if (iirHeadBumpBR < -1.0) iirHeadBumpBR = -1.0;
			iirHeadBumpBR = asin(iirHeadBumpBR);

			inputSampleL = sin(inputSampleL);
			tempSample = (inputSampleL * biquadDL[2]) + biquadDL[7];
			biquadDL[7] = (inputSampleL * biquadDL[3]) - (tempSample * biquadDL[5]) + biquadDL[8];
			biquadDL[8] = (inputSampleL * biquadDL[4]) - (tempSample * biquadDL[6]);
			inputSampleL = tempSample; //interleaved biquad
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
			inputSampleL = asin(inputSampleL);

			inputSampleR = sin(inputSampleR);
			tempSample = (inputSampleR * biquadDR[2]) + biquadDR[7];
			biquadDR[7] = (inputSampleR * biquadDR[3]) - (tempSample * biquadDR[5]) + biquadDR[8];
			biquadDR[8] = (inputSampleR * biquadDR[4]) - (tempSample * biquadDR[6]);
			inputSampleR = tempSample; //interleaved biquad
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
			inputSampleR = asin(inputSampleR);
		}
		flip = !flip;

		// Hypersonic

		outSample = (inputSampleL * fixB[fix_a0]) + fixB[fix_sL1];
		fixB[fix_sL1] = (inputSampleL * fixB[fix_a1]) - (outSample * fixB[fix_b1]) + fixB[fix_sL2];
		fixB[fix_sL2] = (inputSampleL * fixB[fix_a2]) - (outSample * fixB[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixB[fix_a0]) + fixB[fix_sR1];
		fixB[fix_sR1] = (inputSampleR * fixB[fix_a1]) - (outSample * fixB[fix_b1]) + fixB[fix_sR2];
		fixB[fix_sR2] = (inputSampleR * fixB[fix_a2]) - (outSample * fixB[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		double groundSampleL = drySampleL - inputSampleL; //set up UnBox
		double groundSampleR = drySampleR - inputSampleR; //set up UnBox

		if (inputgain != 1.0) {
			inputSampleL *= inputgain;
			inputSampleR *= inputgain;
		} //gain boost inside UnBox: do not boost fringe audio

		double applySoften = fabs(HighsSampleL)*1.57079633;
		if (applySoften > 1.57079633) applySoften = 1.57079633;
		applySoften = 1-cos(applySoften);
		if (HighsSampleL > 0) inputSampleL -= applySoften;
		if (HighsSampleL < 0) inputSampleL += applySoften;
		//apply Soften depending on polarity
		applySoften = fabs(HighsSampleR)*1.57079633;
		if (applySoften > 1.57079633) applySoften = 1.57079633;
		applySoften = 1-cos(applySoften);
		if (HighsSampleR > 0) inputSampleR -= applySoften;
		if (HighsSampleR < 0) inputSampleR += applySoften;
		//apply Soften depending on polarity

		if (inputSampleL > 1.2533141373155) inputSampleL = 1.2533141373155;
		if (inputSampleL < -1.2533141373155) inputSampleL = -1.2533141373155;
		//clip to 1.2533141373155 to reach maximum output
		inputSampleL = sin(inputSampleL * fabs(inputSampleL)) / ((fabs(inputSampleL) == 0.0) ?1:fabs(inputSampleL));
		//Spiral, for cleanest most optimal tape effect
		if (inputSampleR > 1.2533141373155) inputSampleR = 1.2533141373155;
		if (inputSampleR < -1.2533141373155) inputSampleR = -1.2533141373155;
		//clip to 1.2533141373155 to reach maximum output
		inputSampleR = sin(inputSampleR * fabs(inputSampleR)) / ((fabs(inputSampleR) == 0.0) ?1:fabs(inputSampleR));
		//Spiral, for cleanest most optimal tape effect

		double suppress = (1.0-fabs(inputSampleL)) * 0.00013;
		if (iirHeadBumpAL > suppress) iirHeadBumpAL -= suppress;
		if (iirHeadBumpAL < -suppress) iirHeadBumpAL += suppress;
		if (iirHeadBumpBL > suppress) iirHeadBumpBL -= suppress;
		if (iirHeadBumpBL < -suppress) iirHeadBumpBL += suppress;
		//restrain resonant quality of head bump algorithm
		suppress = (1.0-fabs(inputSampleR)) * 0.00013;
		if (iirHeadBumpAR > suppress) iirHeadBumpAR -= suppress;
		if (iirHeadBumpAR < -suppress) iirHeadBumpAR += suppress;
		if (iirHeadBumpBR > suppress) iirHeadBumpBR -= suppress;
		if (iirHeadBumpBR < -suppress) iirHeadBumpBR += suppress;
		//restrain resonant quality of head bump algorithm

		inputSampleL += groundSampleL; //apply UnBox processing
		inputSampleR += groundSampleR; //apply UnBox processing

		inputSampleL += ((iirHeadBumpAL + iirHeadBumpBL) * bumpgain);//and head bump
		inputSampleR += ((iirHeadBumpAR + iirHeadBumpBR) * bumpgain);//and head bump

		// Hypersonic

		outSample = (inputSampleL * fixC[fix_a0]) + fixC[fix_sL1];
		fixC[fix_sL1] = (inputSampleL * fixC[fix_a1]) - (outSample * fixC[fix_b1]) + fixC[fix_sL2];
		fixC[fix_sL2] = (inputSampleL * fixC[fix_a2]) - (outSample * fixC[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixC[fix_a0]) + fixC[fix_sR1];
		fixC[fix_sR1] = (inputSampleR * fixC[fix_a1]) - (outSample * fixC[fix_b1]) + fixC[fix_sR2];
		fixC[fix_sR2] = (inputSampleR * fixC[fix_a2]) - (outSample * fixC[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// uLawEncode

		if(enableuLaw) {
			static int noisesourceL = 0;
			static int noisesourceR = 850010;
			int residue;
			double applyresidue;
		
			noisesourceL = noisesourceL % 1700021; noisesourceL++;
			residue = noisesourceL * noisesourceL;
			residue = residue % 170003; residue *= residue;
			residue = residue % 17011; residue *= residue;
			residue = residue % 1709; residue *= residue;
			residue = residue % 173; residue *= residue;
			residue = residue % 17;
			applyresidue = residue;
			applyresidue *= 0.00000001;
			applyresidue *= 0.00000001;
			inputSampleL += applyresidue;
			if (inputSampleL<1.2e-38 && -inputSampleL<1.2e-38) {
				inputSampleL -= applyresidue;
			}
		
			noisesourceR = noisesourceR % 1700021; noisesourceR++;
			residue = noisesourceR * noisesourceR;
			residue = residue % 170003; residue *= residue;
			residue = residue % 17011; residue *= residue;
			residue = residue % 1709; residue *= residue;
			residue = residue % 173; residue *= residue;
			residue = residue % 17;
			applyresidue = residue;
			applyresidue *= 0.00000001;
			applyresidue *= 0.00000001;
			inputSampleR += applyresidue;
			if (inputSampleR<1.2e-38 && -inputSampleR<1.2e-38) {
				inputSampleR -= applyresidue;
			}
			//for live air, we always apply the dither noise. Then, if our result is 
			//effectively digital black, we'll subtract it auLawEncode. We want a 'air' hiss

			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
		
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
		
			if (inputSampleL > 0) inputSampleL = log(1.0+(255*fabs(inputSampleL))) / log(256);
			if (inputSampleL < 0) inputSampleL = -log(1.0+(255*fabs(inputSampleL))) / log(256);
		
			if (inputSampleR > 0) inputSampleR = log(1.0+(255*fabs(inputSampleR))) / log(256);
			if (inputSampleR < 0) inputSampleR = -log(1.0+(255*fabs(inputSampleR))) / log(256);
		}
	
		// Creature

		if(wet != 0) {
			double cr_drySampleL = inputSampleL;
			double cr_drySampleR = inputSampleR;
		
			for (int x = 0; x < stages; x++) {
				inputSampleL = (slewL[x]+(sin(slewL[x]-inputSampleL)*0.5))*source;
				slewL[x] = inputSampleL*0.5;
				inputSampleR = (slewR[x]+(sin(slewR[x]-inputSampleR)*0.5))*source;
				slewR[x] = inputSampleR*0.5;
			}
			if (stages % 2 > 0) {
				inputSampleL = -inputSampleL;
				inputSampleR = -inputSampleR;
			}
		
			inputSampleL *= wet;
			inputSampleR *= wet;
			cr_drySampleL *= dry;
			cr_drySampleR *= dry;
			inputSampleL += cr_drySampleL;
			inputSampleR += cr_drySampleR;
		}

		// uLawDecode

		if(enableuLaw) {
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;
		
			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;
		
			if (inputSampleL > 0) inputSampleL = (pow(256,fabs(inputSampleL))-1.0) / 255;
			if (inputSampleL < 0) inputSampleL = -(pow(256,fabs(inputSampleL))-1.0) / 255;
		
			if (inputSampleR > 0) inputSampleR = (pow(256,fabs(inputSampleR))-1.0) / 255;
			if (inputSampleR < 0) inputSampleR = -(pow(256,fabs(inputSampleR))-1.0) / 255;
		}

		// BitshiftGain

		inputSampleL *= trim_gain;
		inputSampleR *= trim_gain;

		// Hypersonic

		outSample = (inputSampleL * fixD[fix_a0]) + fixD[fix_sL1];
		fixD[fix_sL1] = (inputSampleL * fixD[fix_a1]) - (outSample * fixD[fix_b1]) + fixD[fix_sL2];
		fixD[fix_sL2] = (inputSampleL * fixD[fix_a2]) - (outSample * fixD[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixD[fix_a0]) + fixD[fix_sR1];
		fixD[fix_sR1] = (inputSampleR * fixD[fix_a1]) - (outSample * fixD[fix_b1]) + fixD[fix_sR2];
		fixD[fix_sR2] = (inputSampleR * fixD[fix_a2]) - (outSample * fixD[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics

		// ToneSlant
		
		for (int count = ts_overallscale; count >= 0; count--) {
			ts_bL[count+1] = ts_bL[count];
			ts_bR[count+1] = ts_bR[count];
		}
		
		ts_bL[0] = ts_accumulatorSampleL = ts_drySampleL = inputSampleL;
		ts_bR[0] = ts_accumulatorSampleR = ts_drySampleR = inputSampleR;
		
		ts_accumulatorSampleL *= ts_f[0];
		ts_accumulatorSampleR *= ts_f[0];

		for (int count = 1; count < ts_overallscale; count++) {
			ts_accumulatorSampleL += (ts_bL[count] * ts_f[count]);
			ts_accumulatorSampleR += (ts_bR[count] * ts_f[count]);
		}
		
		ts_correctionSampleL = inputSampleL - (ts_accumulatorSampleL*2.0);
		ts_correctionSampleR = inputSampleR - (ts_accumulatorSampleR*2.0);
		//we're gonna apply the total effect of all these calculations as a single subtract
		
		inputSampleL += (ts_correctionSampleL * ts_applySlant);
		inputSampleR += (ts_correctionSampleR * ts_applySlant);
		//our one math operation on the input data coming in

		// ADClip
		
		if (lastSampleL >= 0.5)
		{
			if (inputSampleL < 0.5) lastSampleL = ((0.5*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = 0.5;
		}

		if (lastSampleL <= -0.5)
		{
			if (inputSampleL > -0.5) lastSampleL = ((-0.5*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = -0.5;
		}

		if (inputSampleL > 0.5)
		{
			if (lastSampleL < 0.5) inputSampleL = ((0.5*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = 0.5;
		}

		if (inputSampleL < -0.5)
		{
			if (lastSampleL > -0.5) inputSampleL = ((-0.5*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = -0.5;
		}
		lastSampleL = inputSampleL; //end ADClip L


		if (lastSampleR >= 0.5)
		{
			if (inputSampleR < 0.5) lastSampleR = ((0.5*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = 0.5;
		}

		if (lastSampleR <= -0.5)
		{
			if (inputSampleR > -0.5) lastSampleR = ((-0.5*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = -0.5;
		}

		if (inputSampleR > 0.5)
		{
			if (lastSampleR < 0.5) inputSampleR = ((0.5*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = 0.5;
		}

		if (inputSampleR < -0.5)
		{
			if (lastSampleR > -0.5) inputSampleR = ((-0.5*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = -0.5;
		}
		lastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.5) inputSampleL = 0.5;
		if (inputSampleL < -0.5) inputSampleL = -0.5;
		//final iron bar
		if (inputSampleR > 0.5) inputSampleR = 0.5;
		if (inputSampleR < -0.5) inputSampleR = -0.5;
		//final iron bar

		// Hypersonic

		outSample = (inputSampleL * fixE[fix_a0]) + fixE[fix_sL1];
		fixE[fix_sL1] = (inputSampleL * fixE[fix_a1]) - (outSample * fixE[fix_b1]) + fixE[fix_sL2];
		fixE[fix_sL2] = (inputSampleL * fixE[fix_a2]) - (outSample * fixE[fix_b2]);
		inputSampleL = outSample; //fixed biquad filtering ultrasonics
		outSample = (inputSampleR * fixE[fix_a0]) + fixE[fix_sR1];
		fixE[fix_sR1] = (inputSampleR * fixE[fix_a1]) - (outSample * fixE[fix_b1]) + fixE[fix_sR2];
		fixE[fix_sR2] = (inputSampleR * fixE[fix_a2]) - (outSample * fixE[fix_b2]);
		inputSampleR = outSample; //fixed biquad filtering ultrasonics
			//begin 64 bit stereo floating point dither
		//int expon; frexp((double)inputSampleL, &expon);
		fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
		//inputSampleL += ((double(fpdL)-uint32_t(0x7fffffff)) * 1.1e-44l * pow(2,expon+62));
		//frexp((double)inputSampleR, &expon);
		fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;
		//inputSampleR += ((double(fpdR)-uint32_t(0x7fffffff)) * 1.1e-44l * pow(2,expon+62));
		//end 64 bit stereo floating point dither

		*out1 = inputSampleL;
		*out2 = inputSampleR;

		*in1++;
		*in2++;
		*out1++;
		*out2++;
    }
}
