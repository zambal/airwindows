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

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	// BitshiftGain

	int bitshiftGain = (A * 32)-16;
	double gain = 1.0;
	switch (bitshiftGain)
	{
		case -16: gain = 0.0000152587890625; break;
		case -15: gain = 0.000030517578125; break;
		case -14: gain = 0.00006103515625; break;
		case -13: gain = 0.0001220703125; break;
		case -12: gain = 0.000244140625; break;
		case -11: gain = 0.00048828125; break;
		case -10: gain = 0.0009765625; break;
		case -9: gain = 0.001953125; break;
		case -8: gain = 0.00390625; break;
		case -7: gain = 0.0078125; break;
		case -6: gain = 0.015625; break;
		case -5: gain = 0.03125; break;
		case -4: gain = 0.0625; break;
		case -3: gain = 0.125; break;
		case -2: gain = 0.25; break;
		case -1: gain = 0.5; break;
		case 0: gain = 1.0; break;
		case 1: gain = 2.0; break;
		case 2: gain = 4.0; break;
		case 3: gain = 8.0; break;
		case 4: gain = 16.0; break;
		case 5: gain = 32.0; break;
		case 6: gain = 64.0; break;
		case 7: gain = 128.0; break;
		case 8: gain = 256.0; break;
		case 9: gain = 512.0; break;
		case 10: gain = 1024.0; break;
		case 11: gain = 2048.0; break;
		case 12: gain = 4096.0; break;
		case 13: gain = 8192.0; break;
		case 14: gain = 16384.0; break;
		case 15: gain = 32768.0; break;
		case 16: gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// UltrasonicLite

	biquadA[0] = 24000.0 / getSampleRate();
	if (getSampleRate() < 88000.0) biquadA[0] = 21000.0 / getSampleRate();
    biquadA[1] = 0.70710678;

	double KK = tan(M_PI * biquadA[0]); //lowpass
	double u_norm = 1.0 / (1.0 + KK / biquadA[1] + KK * KK);
	biquadA[2] = KK * KK * u_norm;
	biquadA[3] = 2.0 * biquadA[2];
	biquadA[4] = biquadA[2];
	biquadA[5] = 2.0 * (KK * KK - 1.0) * u_norm;
	biquadA[6] = (1.0 - KK / biquadA[1] + KK * KK) * u_norm;

	// Interstage

	double firstStage = 0.381966011250105 / overallscale;
	double iirAmount = 0.00295 / overallscale;
	double threshold = 0.381966011250105;

	// Baxandall2

	double trebleGain = pow(10.0,((B*48.0)-24.0)/20.0);
	double trebleFreq = (4410.0*trebleGain)/getSampleRate();
	if (trebleFreq > 0.45) trebleFreq = 0.45;
	trebleAL[0] = trebleBL[0] = trebleAR[0] = trebleBR[0] = trebleFreq;
	double bassGain = pow(10.0,((C*48.0)-24.0)/20.0);
	double bassFreq = pow(10.0,-((C*48.0)-24.0)/20.0);
	bassFreq = (8820.0*bassFreq)/getSampleRate();
	if (bassFreq > 0.45) bassFreq = 0.45;
	bassAL[0] = bassBL[0] = bassAR[0] = bassBR[0] = bassFreq;
    trebleAL[1] = trebleBL[1] = trebleAR[1] = trebleBR[1] = 0.4;
    bassAL[1] = bassBL[1] = bassAR[1] = bassBR[1] = 0.2;

	double b_K = tan(M_PI * trebleAL[0]);
	double b_norm = 1.0 / (1.0 + b_K / trebleAL[1] + b_K * b_K);
	trebleBL[2] = trebleAL[2] = trebleBR[2] = trebleAR[2] = b_K * b_K * b_norm;
	trebleBL[3] = trebleAL[3] = trebleBR[3] = trebleAR[3] = 2.0 * trebleAL[2];
	trebleBL[4] = trebleAL[4] = trebleBR[4] = trebleAR[4] = trebleAL[2];
	trebleBL[5] = trebleAL[5] = trebleBR[5] = trebleAR[5] = 2.0 * (b_K * b_K - 1.0) * b_norm;
	trebleBL[6] = trebleAL[6] = trebleBR[6] = trebleAR[6] = (1.0 - b_K / trebleAL[1] + b_K * b_K) * b_norm;

	b_K = tan(M_PI * bassAL[0]);
	b_norm = 1.0 / (1.0 + b_K / bassAL[1] + b_K * b_K);
	bassBL[2] = bassAL[2] = bassBR[2] = bassAR[2] = b_K * b_K * b_norm;
	bassBL[3] = bassAL[3] = bassBR[3] = bassAR[3] = 2.0 * bassAL[2];
	bassBL[4] = bassAL[4] = bassBR[4] = bassAR[4] = bassAL[2];
	bassBL[5] = bassAL[5] = bassBR[5] = bassAR[5] = 2.0 * (b_K * b_K - 1.0) * b_norm;
	bassBL[6] = bassAL[6] = bassBR[6] = bassAR[6] = (1.0 - b_K / bassAL[1] + b_K * b_K) * b_norm;

	// Tape

	double inputgain = pow(10.0,((D-0.5)*24.0)/20.0);
	double bumpgain = E*0.1;
	double HeadBumpFreq = 0.12/overallscale;
	double softness = 0.618033988749894848204586;
	double RollAmount = (1.0 - softness) / overallscale;
	//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
	//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
	biquadAL[0] = biquadBL[0] = biquadAR[0] = biquadBR[0] = 0.0072/overallscale;
	biquadAL[1] = biquadBL[1] = biquadAR[1] = biquadBR[1] = 0.0009;
	double K = tan(M_PI * biquadBR[0]);
	double norm = 1.0 / (1.0 + K / biquadBR[1] + K * K);
	biquadAL[2] = biquadBL[2] = biquadAR[2] = biquadBR[2] = K / biquadBR[1] * norm;
	biquadAL[4] = biquadBL[4] = biquadAR[4] = biquadBR[4] = -biquadBR[2];
	biquadAL[5] = biquadBL[5] = biquadAR[5] = biquadBR[5] = 2.0 * (K * K - 1.0) * norm;
	biquadAL[6] = biquadBL[6] = biquadAR[6] = biquadBR[6] = (1.0 - K / biquadBR[1] + K * K) * norm;

	biquadCL[0] = biquadDL[0] = biquadCR[0] = biquadDR[0] = 0.032/overallscale;
	biquadCL[1] = biquadDL[1] = biquadCR[1] = biquadDR[1] = 0.0007;
	K = tan(M_PI * biquadDR[0]);
	norm = 1.0 / (1.0 + K / biquadDR[1] + K * K);
	biquadCL[2] = biquadDL[2] = biquadCR[2] = biquadDR[2] = K / biquadDR[1] * norm;
	biquadCL[4] = biquadDL[4] = biquadCR[4] = biquadDR[4] = -biquadDR[2];
	biquadCL[5] = biquadDL[5] = biquadCR[5] = biquadDR[5] = 2.0 * (K * K - 1.0) * norm;
	biquadCL[6] = biquadDL[6] = biquadCR[6] = biquadDR[6] = (1.0 - K / biquadDR[1] + K * K) * norm;

    while (--sampleFrames >= 0)
    {

		// BitshiftGain

		double inputSampleL = *in1 * gain;
		double inputSampleR = *in2 * gain;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// UltrasonicLite

		double outSampleL = biquadA[2]*inputSampleL+biquadA[3]*biquadA[7]+biquadA[4]*biquadA[8]-biquadA[5]*biquadA[9]-biquadA[6]*biquadA[10];
		biquadA[8] = biquadA[7]; biquadA[7] = inputSampleL; inputSampleL = outSampleL; biquadA[10] = biquadA[9]; biquadA[9] = inputSampleL; //DF1 left

		double outSampleR = biquadA[2]*inputSampleR+biquadA[3]*biquadA[11]+biquadA[4]*biquadA[12]-biquadA[5]*biquadA[13]-biquadA[6]*biquadA[14];
		biquadA[12] = biquadA[11]; biquadA[11] = inputSampleR; inputSampleR = outSampleR; biquadA[14] = biquadA[13]; biquadA[13] = inputSampleR; //DF1 right

		// Interstage

		double drySampleL = inputSampleL;
		double drySampleR = inputSampleR;

		if (i_flip) {
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
		} else {
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
		}
		i_flip = !i_flip;
		i_lastSampleL = inputSampleL;
		i_lastSampleR = inputSampleR;

		// Baxandall2

		double trebleSampleL;
		double bassSampleL;
		double trebleSampleR;
		double bassSampleR;

		if (flip)
		{
			trebleSampleL = (inputSampleL * trebleAL[2]) + trebleAL[7];
			trebleAL[7] = (inputSampleL * trebleAL[3]) - (trebleSampleL * trebleAL[5]) + trebleAL[8];
			trebleAL[8] = (inputSampleL * trebleAL[4]) - (trebleSampleL * trebleAL[6]);
			trebleSampleL = inputSampleL - trebleSampleL;

			bassSampleL = (inputSampleL * bassAL[2]) + bassAL[7];
			bassAL[7] = (inputSampleL * bassAL[3]) - (bassSampleL * bassAL[5]) + bassAL[8];
			bassAL[8] = (inputSampleL * bassAL[4]) - (bassSampleL * bassAL[6]);

			trebleSampleR = (inputSampleR * trebleAR[2]) + trebleAR[7];
			trebleAR[7] = (inputSampleR * trebleAR[3]) - (trebleSampleR * trebleAR[5]) + trebleAR[8];
			trebleAR[8] = (inputSampleR * trebleAR[4]) - (trebleSampleR * trebleAR[6]);
			trebleSampleR = inputSampleR - trebleSampleR;

			bassSampleR = (inputSampleR * bassAR[2]) + bassAR[7];
			bassAR[7] = (inputSampleR * bassAR[3]) - (bassSampleR * bassAR[5]) + bassAR[8];
			bassAR[8] = (inputSampleR * bassAR[4]) - (bassSampleR * bassAR[6]);
		}
		else
		{
			trebleSampleL = (inputSampleL * trebleBL[2]) + trebleBL[7];
			trebleBL[7] = (inputSampleL * trebleBL[3]) - (trebleSampleL * trebleBL[5]) + trebleBL[8];
			trebleBL[8] = (inputSampleL * trebleBL[4]) - (trebleSampleL * trebleBL[6]);
			trebleSampleL = inputSampleL - trebleSampleL;

			bassSampleL = (inputSampleL * bassBL[2]) + bassBL[7];
			bassBL[7] = (inputSampleL * bassBL[3]) - (bassSampleL * bassBL[5]) + bassBL[8];
			bassBL[8] = (inputSampleL * bassBL[4]) - (bassSampleL * bassBL[6]);

			trebleSampleR = (inputSampleR * trebleBR[2]) + trebleBR[7];
			trebleBR[7] = (inputSampleR * trebleBR[3]) - (trebleSampleR * trebleBR[5]) + trebleBR[8];
			trebleBR[8] = (inputSampleR * trebleBR[4]) - (trebleSampleR * trebleBR[6]);
			trebleSampleR = inputSampleR - trebleSampleR;

			bassSampleR = (inputSampleR * bassBR[2]) + bassBR[7];
			bassBR[7] = (inputSampleR * bassBR[3]) - (bassSampleR * bassBR[5]) + bassBR[8];
			bassBR[8] = (inputSampleR * bassBR[4]) - (bassSampleR * bassBR[6]);
		}

		trebleSampleL *= trebleGain;
		bassSampleL *= bassGain;
		inputSampleL = bassSampleL + trebleSampleL; //interleaved biquad
		trebleSampleR *= trebleGain;
		bassSampleR *= bassGain;
		inputSampleR = bassSampleR + trebleSampleR; //interleaved biquad


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

		if (lastSampleL >= 0.333333333333333333)
		{
			if (inputSampleL < 0.333333333333333333) lastSampleL = ((0.333333333333333333*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = 0.333333333333333333;
		}

		if (lastSampleL <= -0.333333333333333333)
		{
			if (inputSampleL > -0.333333333333333333) lastSampleL = ((-0.333333333333333333*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = -0.333333333333333333;
		}

		if (inputSampleL > 0.333333333333333333)
		{
			if (lastSampleL < 0.333333333333333333) inputSampleL = ((0.333333333333333333*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = 0.333333333333333333;
		}

		if (inputSampleL < -0.333333333333333333)
		{
			if (lastSampleL > -0.333333333333333333) inputSampleL = ((-0.333333333333333333*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = -0.333333333333333333;
		}
		lastSampleL = inputSampleL; //end ADClip L


		if (lastSampleR >= 0.333333333333333333)
		{
			if (inputSampleR < 0.333333333333333333) lastSampleR = ((0.333333333333333333*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = 0.333333333333333333;
		}

		if (lastSampleR <= -0.333333333333333333)
		{
			if (inputSampleR > -0.333333333333333333) lastSampleR = ((-0.333333333333333333*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = -0.333333333333333333;
		}

		if (inputSampleR > 0.333333333333333333)
		{
			if (lastSampleR < 0.333333333333333333) inputSampleR = ((0.333333333333333333*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = 0.333333333333333333;
		}

		if (inputSampleR < -0.333333333333333333)
		{
			if (lastSampleR > -0.333333333333333333) inputSampleR = ((-0.333333333333333333*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = -0.333333333333333333;
		}
		lastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.333333333333333333) inputSampleL = 0.333333333333333333;
		if (inputSampleL < -0.333333333333333333) inputSampleL = -0.333333333333333333;
		//final iron bar
		if (inputSampleR > 0.333333333333333333) inputSampleR = 0.333333333333333333;
		if (inputSampleR < -0.333333333333333333) inputSampleR = -0.333333333333333333;
		//final iron bar

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

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	// BitshiftGain

	int bitshiftGain = (A * 32)-16;
	double gain = 1.0;
	switch (bitshiftGain)
	{
		case -16: gain = 0.0000152587890625; break;
		case -15: gain = 0.000030517578125; break;
		case -14: gain = 0.00006103515625; break;
		case -13: gain = 0.0001220703125; break;
		case -12: gain = 0.000244140625; break;
		case -11: gain = 0.00048828125; break;
		case -10: gain = 0.0009765625; break;
		case -9: gain = 0.001953125; break;
		case -8: gain = 0.00390625; break;
		case -7: gain = 0.0078125; break;
		case -6: gain = 0.015625; break;
		case -5: gain = 0.03125; break;
		case -4: gain = 0.0625; break;
		case -3: gain = 0.125; break;
		case -2: gain = 0.25; break;
		case -1: gain = 0.5; break;
		case 0: gain = 1.0; break;
		case 1: gain = 2.0; break;
		case 2: gain = 4.0; break;
		case 3: gain = 8.0; break;
		case 4: gain = 16.0; break;
		case 5: gain = 32.0; break;
		case 6: gain = 64.0; break;
		case 7: gain = 128.0; break;
		case 8: gain = 256.0; break;
		case 9: gain = 512.0; break;
		case 10: gain = 1024.0; break;
		case 11: gain = 2048.0; break;
		case 12: gain = 4096.0; break;
		case 13: gain = 8192.0; break;
		case 14: gain = 16384.0; break;
		case 15: gain = 32768.0; break;
		case 16: gain = 65536.0; break;
	}
	//we are directly punching in the gain values rather than calculating them

	// UltrasonicLite

	biquadA[0] = 24000.0 / getSampleRate();
	if (getSampleRate() < 88000.0) biquadA[0] = 21000.0 / getSampleRate();
    biquadA[1] = 0.70710678;

	double KK = tan(M_PI * biquadA[0]); //lowpass
	double u_norm = 1.0 / (1.0 + KK / biquadA[1] + KK * KK);
	biquadA[2] = KK * KK * u_norm;
	biquadA[3] = 2.0 * biquadA[2];
	biquadA[4] = biquadA[2];
	biquadA[5] = 2.0 * (KK * KK - 1.0) * u_norm;
	biquadA[6] = (1.0 - KK / biquadA[1] + KK * KK) * u_norm;

	// Interstage

	double firstStage = 0.381966011250105 / overallscale;
	double iirAmount = 0.00295 / overallscale;
	double threshold = 0.381966011250105;

	// Baxandall2

	double trebleGain = pow(10.0,((B*48.0)-24.0)/20.0);
	double trebleFreq = (4410.0*trebleGain)/getSampleRate();
	if (trebleFreq > 0.45) trebleFreq = 0.45;
	trebleAL[0] = trebleBL[0] = trebleAR[0] = trebleBR[0] = trebleFreq;
	double bassGain = pow(10.0,((C*48.0)-24.0)/20.0);
	double bassFreq = pow(10.0,-((C*48.0)-24.0)/20.0);
	bassFreq = (8820.0*bassFreq)/getSampleRate();
	if (bassFreq > 0.45) bassFreq = 0.45;
	bassAL[0] = bassBL[0] = bassAR[0] = bassBR[0] = bassFreq;
    trebleAL[1] = trebleBL[1] = trebleAR[1] = trebleBR[1] = 0.4;
    bassAL[1] = bassBL[1] = bassAR[1] = bassBR[1] = 0.2;

	double b_K = tan(M_PI * trebleAL[0]);
	double b_norm = 1.0 / (1.0 + b_K / trebleAL[1] + b_K * b_K);
	trebleBL[2] = trebleAL[2] = trebleBR[2] = trebleAR[2] = b_K * b_K * b_norm;
	trebleBL[3] = trebleAL[3] = trebleBR[3] = trebleAR[3] = 2.0 * trebleAL[2];
	trebleBL[4] = trebleAL[4] = trebleBR[4] = trebleAR[4] = trebleAL[2];
	trebleBL[5] = trebleAL[5] = trebleBR[5] = trebleAR[5] = 2.0 * (b_K * b_K - 1.0) * b_norm;
	trebleBL[6] = trebleAL[6] = trebleBR[6] = trebleAR[6] = (1.0 - b_K / trebleAL[1] + b_K * b_K) * b_norm;

	b_K = tan(M_PI * bassAL[0]);
	b_norm = 1.0 / (1.0 + b_K / bassAL[1] + b_K * b_K);
	bassBL[2] = bassAL[2] = bassBR[2] = bassAR[2] = b_K * b_K * b_norm;
	bassBL[3] = bassAL[3] = bassBR[3] = bassAR[3] = 2.0 * bassAL[2];
	bassBL[4] = bassAL[4] = bassBR[4] = bassAR[4] = bassAL[2];
	bassBL[5] = bassAL[5] = bassBR[5] = bassAR[5] = 2.0 * (b_K * b_K - 1.0) * b_norm;
	bassBL[6] = bassAL[6] = bassBR[6] = bassAR[6] = (1.0 - b_K / bassAL[1] + b_K * b_K) * b_norm;

	// Tape

	double inputgain = pow(10.0,((D-0.5)*24.0)/20.0);
	double bumpgain = E*0.1;
	double HeadBumpFreq = 0.12/overallscale;
	double softness = 0.618033988749894848204586;
	double RollAmount = (1.0 - softness) / overallscale;
	//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
	//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
	biquadAL[0] = biquadBL[0] = biquadAR[0] = biquadBR[0] = 0.0072/overallscale;
	biquadAL[1] = biquadBL[1] = biquadAR[1] = biquadBR[1] = 0.0009;
	double K = tan(M_PI * biquadBR[0]);
	double norm = 1.0 / (1.0 + K / biquadBR[1] + K * K);
	biquadAL[2] = biquadBL[2] = biquadAR[2] = biquadBR[2] = K / biquadBR[1] * norm;
	biquadAL[4] = biquadBL[4] = biquadAR[4] = biquadBR[4] = -biquadBR[2];
	biquadAL[5] = biquadBL[5] = biquadAR[5] = biquadBR[5] = 2.0 * (K * K - 1.0) * norm;
	biquadAL[6] = biquadBL[6] = biquadAR[6] = biquadBR[6] = (1.0 - K / biquadBR[1] + K * K) * norm;

	biquadCL[0] = biquadDL[0] = biquadCR[0] = biquadDR[0] = 0.032/overallscale;
	biquadCL[1] = biquadDL[1] = biquadCR[1] = biquadDR[1] = 0.0007;
	K = tan(M_PI * biquadDR[0]);
	norm = 1.0 / (1.0 + K / biquadDR[1] + K * K);
	biquadCL[2] = biquadDL[2] = biquadCR[2] = biquadDR[2] = K / biquadDR[1] * norm;
	biquadCL[4] = biquadDL[4] = biquadCR[4] = biquadDR[4] = -biquadDR[2];
	biquadCL[5] = biquadDL[5] = biquadCR[5] = biquadDR[5] = 2.0 * (K * K - 1.0) * norm;
	biquadCL[6] = biquadDL[6] = biquadCR[6] = biquadDR[6] = (1.0 - K / biquadDR[1] + K * K) * norm;

    while (--sampleFrames >= 0)
    {

		// BitshiftGain

		double inputSampleL = *in1 * gain;
		double inputSampleR = *in2 * gain;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// UltrasonicLite

		double outSampleL = biquadA[2]*inputSampleL+biquadA[3]*biquadA[7]+biquadA[4]*biquadA[8]-biquadA[5]*biquadA[9]-biquadA[6]*biquadA[10];
		biquadA[8] = biquadA[7]; biquadA[7] = inputSampleL; inputSampleL = outSampleL; biquadA[10] = biquadA[9]; biquadA[9] = inputSampleL; //DF1 left

		double outSampleR = biquadA[2]*inputSampleR+biquadA[3]*biquadA[11]+biquadA[4]*biquadA[12]-biquadA[5]*biquadA[13]-biquadA[6]*biquadA[14];
		biquadA[12] = biquadA[11]; biquadA[11] = inputSampleR; inputSampleR = outSampleR; biquadA[14] = biquadA[13]; biquadA[13] = inputSampleR; //DF1 right

		// Interstage

		double drySampleL = inputSampleL;
		double drySampleR = inputSampleR;

		if (i_flip) {
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
		} else {
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
		}
		i_flip = !i_flip;
		i_lastSampleL = inputSampleL;
		i_lastSampleR = inputSampleR;

		// Baxandall2

		double trebleSampleL;
		double bassSampleL;
		double trebleSampleR;
		double bassSampleR;

		if (flip)
		{
			trebleSampleL = (inputSampleL * trebleAL[2]) + trebleAL[7];
			trebleAL[7] = (inputSampleL * trebleAL[3]) - (trebleSampleL * trebleAL[5]) + trebleAL[8];
			trebleAL[8] = (inputSampleL * trebleAL[4]) - (trebleSampleL * trebleAL[6]);
			trebleSampleL = inputSampleL - trebleSampleL;

			bassSampleL = (inputSampleL * bassAL[2]) + bassAL[7];
			bassAL[7] = (inputSampleL * bassAL[3]) - (bassSampleL * bassAL[5]) + bassAL[8];
			bassAL[8] = (inputSampleL * bassAL[4]) - (bassSampleL * bassAL[6]);

			trebleSampleR = (inputSampleR * trebleAR[2]) + trebleAR[7];
			trebleAR[7] = (inputSampleR * trebleAR[3]) - (trebleSampleR * trebleAR[5]) + trebleAR[8];
			trebleAR[8] = (inputSampleR * trebleAR[4]) - (trebleSampleR * trebleAR[6]);
			trebleSampleR = inputSampleR - trebleSampleR;

			bassSampleR = (inputSampleR * bassAR[2]) + bassAR[7];
			bassAR[7] = (inputSampleR * bassAR[3]) - (bassSampleR * bassAR[5]) + bassAR[8];
			bassAR[8] = (inputSampleR * bassAR[4]) - (bassSampleR * bassAR[6]);
		}
		else
		{
			trebleSampleL = (inputSampleL * trebleBL[2]) + trebleBL[7];
			trebleBL[7] = (inputSampleL * trebleBL[3]) - (trebleSampleL * trebleBL[5]) + trebleBL[8];
			trebleBL[8] = (inputSampleL * trebleBL[4]) - (trebleSampleL * trebleBL[6]);
			trebleSampleL = inputSampleL - trebleSampleL;

			bassSampleL = (inputSampleL * bassBL[2]) + bassBL[7];
			bassBL[7] = (inputSampleL * bassBL[3]) - (bassSampleL * bassBL[5]) + bassBL[8];
			bassBL[8] = (inputSampleL * bassBL[4]) - (bassSampleL * bassBL[6]);

			trebleSampleR = (inputSampleR * trebleBR[2]) + trebleBR[7];
			trebleBR[7] = (inputSampleR * trebleBR[3]) - (trebleSampleR * trebleBR[5]) + trebleBR[8];
			trebleBR[8] = (inputSampleR * trebleBR[4]) - (trebleSampleR * trebleBR[6]);
			trebleSampleR = inputSampleR - trebleSampleR;

			bassSampleR = (inputSampleR * bassBR[2]) + bassBR[7];
			bassBR[7] = (inputSampleR * bassBR[3]) - (bassSampleR * bassBR[5]) + bassBR[8];
			bassBR[8] = (inputSampleR * bassBR[4]) - (bassSampleR * bassBR[6]);
		}

		trebleSampleL *= trebleGain;
		bassSampleL *= bassGain;
		inputSampleL = bassSampleL + trebleSampleL; //interleaved biquad
		trebleSampleR *= trebleGain;
		bassSampleR *= bassGain;
		inputSampleR = bassSampleR + trebleSampleR; //interleaved biquad


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

		if (lastSampleL >= 0.333333333333333333)
		{
			if (inputSampleL < 0.333333333333333333) lastSampleL = ((0.333333333333333333*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = 0.333333333333333333;
		}

		if (lastSampleL <= -0.333333333333333333)
		{
			if (inputSampleL > -0.333333333333333333) lastSampleL = ((-0.333333333333333333*softness) + (inputSampleL * (1.0-softness)));
			else lastSampleL = -0.333333333333333333;
		}

		if (inputSampleL > 0.333333333333333333)
		{
			if (lastSampleL < 0.333333333333333333) inputSampleL = ((0.333333333333333333*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = 0.333333333333333333;
		}

		if (inputSampleL < -0.333333333333333333)
		{
			if (lastSampleL > -0.333333333333333333) inputSampleL = ((-0.333333333333333333*softness) + (lastSampleL * (1.0-softness)));
			else inputSampleL = -0.333333333333333333;
		}
		lastSampleL = inputSampleL; //end ADClip L


		if (lastSampleR >= 0.333333333333333333)
		{
			if (inputSampleR < 0.333333333333333333) lastSampleR = ((0.333333333333333333*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = 0.333333333333333333;
		}

		if (lastSampleR <= -0.333333333333333333)
		{
			if (inputSampleR > -0.333333333333333333) lastSampleR = ((-0.333333333333333333*softness) + (inputSampleR * (1.0-softness)));
			else lastSampleR = -0.333333333333333333;
		}

		if (inputSampleR > 0.333333333333333333)
		{
			if (lastSampleR < 0.333333333333333333) inputSampleR = ((0.333333333333333333*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = 0.333333333333333333;
		}

		if (inputSampleR < -0.333333333333333333)
		{
			if (lastSampleR > -0.333333333333333333) inputSampleR = ((-0.333333333333333333*softness) + (lastSampleR * (1.0-softness)));
			else inputSampleR = -0.333333333333333333;
		}
		lastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.333333333333333333) inputSampleL = 0.333333333333333333;
		if (inputSampleL < -0.333333333333333333) inputSampleL = -0.333333333333333333;
		//final iron bar
		if (inputSampleR > 0.333333333333333333) inputSampleR = 0.333333333333333333;
		if (inputSampleR < -0.333333333333333333) inputSampleR = -0.333333333333333333;
		//final iron bar

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
