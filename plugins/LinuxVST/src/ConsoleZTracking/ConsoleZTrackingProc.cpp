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

	bool enableHPF = B != 0.0;
	bool enableLPF = D != 1.0;
	bool enableComp = F != 0.0;
	bool enableChorus = J != 0.0;
	bool enableReverb = L != 0.0;
	bool enableuLaw = Q != 0.0;
	bool enableCreature = V != 0.5;
	bool enableToneSlant = X != 0.5;

	double inputSampleL;
	double inputSampleR;
	double drySampleL;
	double drySampleR;
	double temp;	

	double sampleRate = getSampleRate() * 2.0;
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= sampleRate;

	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;
	//this is going to be 2 for 88.1 or 96k, 3 for silly people, 4 for 176 or 192k

	double outSample;
	double KK;
	double norm;

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

	if (cc_cycle > cycleEnd-1) cc_cycle = cycleEnd-1; //sanity check
	
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
	ch9_biquadB[0] = ch9_biquadA[0] = ch9_cutoff / sampleRate;
    ch9_biquadA[1] = 1.618033988749894848204586;
	ch9_biquadB[1] = 0.618033988749894848204586;
	
	KK = tan(M_PI * ch9_biquadA[0]); //lowpass
	norm = 1.0 / (1.0 + KK / ch9_biquadA[1] + KK * KK);
	ch9_biquadA[2] = KK * KK * norm;
	ch9_biquadA[3] = 2.0 * ch9_biquadA[2];
	ch9_biquadA[4] = ch9_biquadA[2];
	ch9_biquadA[5] = 2.0 * (KK * KK - 1.0) * norm;
	ch9_biquadA[6] = (1.0 - KK / ch9_biquadA[1] + KK * KK) * norm;
	
	KK = tan(M_PI * ch9_biquadA[0]);
	norm = 1.0 / (1.0 + KK / ch9_biquadB[1] + KK * KK);
	ch9_biquadB[2] = KK * KK * norm;
	ch9_biquadB[3] = 2.0 * ch9_biquadB[2];
	ch9_biquadB[4] = ch9_biquadB[2];
	ch9_biquadB[5] = 2.0 * (KK * KK - 1.0) * norm;
	ch9_biquadB[6] = (1.0 - KK / ch9_biquadB[1] + KK * KK) * norm;

	// Focus

	double boost = pow(10.0,(M*12.0)/20.0);
	figureL[0] = figureR[0] = (pow(M_E, freqX * O) * centerFreq)/sampleRate; //fixed frequency, 3.515775k
	figureL[1] = figureR[1] = pow(pow(N,3)*2,2)+0.0001; //resonance
	int mode = (int) ( P * 4.999 );
	double output = R;
	double wet = S;
	
	KK = tan(M_PI * figureR[0]);
	norm = 1.0 / (1.0 + KK / figureR[1] + KK * KK);
	figureL[2] = figureR[2] = KK / figureR[1] * norm;
	figureL[4] = figureR[4] = -figureR[2];
	figureL[5] = figureR[5] = 2.0 * (KK * KK - 1.0) * norm;
	figureL[6] = figureR[6] = (1.0 - KK / figureR[1] + KK * KK) * norm;

	// Creature

	double source = 1.0-pow(1.0-T,5);
	int stages = (pow(U,2)*32.0*sqrt(overallscale))+1;
	double cr_wet = (V*2.0)-1.0; //inv-dry-wet for highpass
	double cr_dry = 2.0-(V*2.0);
	if (cr_dry > 1.0) cr_dry = 1.0; //full dry for use with inv, to 0.0 at full wet

	// ToneSlant

	double ts_correctionSampleL;
	double ts_correctionSampleR;
	double ts_accumulatorSampleL;
	double ts_accumulatorSampleR;
	double ts_overallscale = (W*99.0)+1.0;
	double ts_applySlant = (X*2.0)-1.0;


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

	// +/- 12dB gain trim
	double gain_trim = pow(M_E, freqX * Y) * 0.25;

	while (sampleFrames > 0)
  {

		if(flip) {
			inputSampleL = fpdL * 1.18e-17;
			inputSampleR = fpdR * 1.18e-17;
		}
		else {
			inputSampleL = *in1 * input_gain * 2.0;
			inputSampleR = *in2 * input_gain * 2.0;

			if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
			if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;
		}

		// Channel9
		temp = ch9_biquadA[2]*inputSampleL+ch9_biquadA[3]*ch9_biquadA[7]+ch9_biquadA[4]*ch9_biquadA[8]-ch9_biquadA[5]*ch9_biquadA[9]-ch9_biquadA[6]*ch9_biquadA[10];
		ch9_biquadA[8] = ch9_biquadA[7]; ch9_biquadA[7] = inputSampleL; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleL = temp;
		ch9_biquadA[10] = ch9_biquadA[9]; ch9_biquadA[9] = inputSampleL; //DF1 left
		temp = ch9_biquadA[2]*inputSampleR+ch9_biquadA[3]*ch9_biquadA[11]+ch9_biquadA[4]*ch9_biquadA[12]-ch9_biquadA[5]*ch9_biquadA[13]-ch9_biquadA[6]*ch9_biquadA[14];
		ch9_biquadA[12] = ch9_biquadA[11]; ch9_biquadA[11] = inputSampleR; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleR = temp;
		ch9_biquadA[14] = ch9_biquadA[13]; ch9_biquadA[13] = inputSampleR; //DF1 right

		//inputSampleL = sin(inputSampleL);
		inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
		//inputSampleR = sin(inputSampleR);
		inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));

		double dielectricScaleL = fabs(2.0-((inputSampleL+nonLin)/nonLin));
		double dielectricScaleR = fabs(2.0-((inputSampleR+nonLin)/nonLin));

		if (flip)
		{

			if (fabs(ch9_iirSampleLA)<1.18e-37) ch9_iirSampleLA = 0.0; 
			ch9_iirSampleLA = (ch9_iirSampleLA * (1.0 - (ch9_localiirAmount * dielectricScaleL))) + (inputSampleL * ch9_localiirAmount * dielectricScaleL);
			inputSampleL = inputSampleL - ch9_iirSampleLA;
			if (fabs(ch9_iirSampleRA)<1.18e-37) ch9_iirSampleRA = 0.0; 
			ch9_iirSampleRA = (ch9_iirSampleRA * (1.0 - (ch9_localiirAmount * dielectricScaleR))) + (inputSampleR * ch9_localiirAmount * dielectricScaleR);
			inputSampleR = inputSampleR - ch9_iirSampleRA;

			if(enableHPF) {
				temp = (inputSampleL * biquad_hpf_AL[2]) + biquad_hpf_AL[7];
				biquad_hpf_AL[7] = (inputSampleL * biquad_hpf_AL[3]) - (temp * biquad_hpf_AL[5]) + biquad_hpf_AL[8];
				biquad_hpf_AL[8] = (inputSampleL * biquad_hpf_AL[4]) - (temp * biquad_hpf_AL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_hpf_AR[2]) + biquad_hpf_AR[7];
				biquad_hpf_AR[7] = (inputSampleR * biquad_hpf_AR[3]) - (temp * biquad_hpf_AR[5]) + biquad_hpf_AR[8];
				biquad_hpf_AR[8] = (inputSampleR * biquad_hpf_AR[4]) - (temp * biquad_hpf_AR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}

			if(enableLPF) {
				temp = (inputSampleL * biquad_lpf_AL[2]) + biquad_lpf_AL[7];
				biquad_lpf_AL[7] = (inputSampleL * biquad_lpf_AL[3]) - (temp * biquad_lpf_AL[5]) + biquad_lpf_AL[8];
				biquad_lpf_AL[8] = (inputSampleL * biquad_lpf_AL[4]) - (temp * biquad_lpf_AL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_lpf_AR[2]) + biquad_lpf_AR[7];
				biquad_lpf_AR[7] = (inputSampleR * biquad_lpf_AR[3]) - (temp * biquad_lpf_AR[5]) + biquad_lpf_AR[8];
				biquad_lpf_AR[8] = (inputSampleR * biquad_lpf_AR[4]) - (temp * biquad_lpf_AR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}
		} else {

			if (fabs(ch9_iirSampleLB)<1.18e-37) ch9_iirSampleLB = 0.0; 
			ch9_iirSampleLB = (ch9_iirSampleLB * (1.0 - (ch9_localiirAmount * dielectricScaleL))) + (inputSampleL * ch9_localiirAmount * dielectricScaleL);
			inputSampleL = inputSampleL - ch9_iirSampleLB;
			if (fabs(ch9_iirSampleRB)<1.18e-37) ch9_iirSampleRB = 0.0; 
			ch9_iirSampleRB = (ch9_iirSampleRB * (1.0 - (ch9_localiirAmount * dielectricScaleR))) + (inputSampleR * ch9_localiirAmount * dielectricScaleR);
			inputSampleR = inputSampleR - ch9_iirSampleRB;

			if(enableHPF) {
				temp = (inputSampleL * biquad_hpf_BL[2]) + biquad_hpf_BL[7];
				biquad_hpf_BL[7] = (inputSampleL * biquad_hpf_BL[3]) - (temp * biquad_hpf_BL[5]) + biquad_hpf_BL[8];
				biquad_hpf_BL[8] = (inputSampleL * biquad_hpf_BL[4]) - (temp * biquad_hpf_BL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_hpf_BR[2]) + biquad_hpf_BR[7];
				biquad_hpf_BR[7] = (inputSampleR * biquad_hpf_BR[3]) - (temp * biquad_hpf_BR[5]) + biquad_hpf_BR[8];
				biquad_hpf_BR[8] = (inputSampleR * biquad_hpf_BR[4]) - (temp * biquad_hpf_BR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}

			if(enableLPF) {
				temp = (inputSampleL * biquad_lpf_BL[2]) + biquad_lpf_BL[7];
				biquad_lpf_BL[7] = (inputSampleL * biquad_lpf_BL[3]) - (temp * biquad_lpf_BL[5]) + biquad_lpf_BL[8];
				biquad_lpf_BL[8] = (inputSampleL * biquad_lpf_BL[4]) - (temp * biquad_lpf_BL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_lpf_BR[2]) + biquad_lpf_BR[7];
				biquad_lpf_BR[7] = (inputSampleR * biquad_lpf_BR[3]) - (temp * biquad_lpf_BR[5]) + biquad_lpf_BR[8];
				biquad_lpf_BR[8] = (inputSampleR * biquad_lpf_BR[4]) - (temp * biquad_lpf_BR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}
		}

		//inputSampleL = asin(inputSampleL);
		inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);
		//inputSampleR = asin(inputSampleR);
		inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);

		// Focus

		drySampleL = inputSampleL;
		drySampleR = inputSampleR;
		
		inputSampleL = sin(inputSampleL);
		inputSampleR = sin(inputSampleR);
		//encode Console5: good cleanness
		
		temp = (inputSampleL * figureL[2]) + figureL[7];
		figureL[7] = -(temp * figureL[5]) + figureL[8];
		figureL[8] = (inputSampleL * figureL[4]) - (temp * figureL[6]);
		inputSampleL = temp;
		
		temp = (inputSampleR * figureR[2]) + figureR[7];
		figureR[7] = -(temp * figureR[5]) + figureR[8];
		figureR[8] = (inputSampleR * figureR[4]) - (temp * figureR[6]);
		inputSampleR = temp;
				
		if (inputSampleL > 1.0) inputSampleL = 1.0;
		if (inputSampleL < -1.0) inputSampleL = -1.0;
		if (inputSampleR > 1.0) inputSampleR = 1.0;
		if (inputSampleR < -1.0) inputSampleR = -1.0;
		//without this, you can get a NaN condition where it spits out DC offset at full blast!
		inputSampleL = asin(inputSampleL);
		inputSampleR = asin(inputSampleR);
		//decode Console5
		
		double groundSampleL = drySampleL - inputSampleL; //set up UnBox
		double groundSampleR = drySampleR - inputSampleR; //set up UnBox

		// uLawEncode

		if(enableuLaw) {
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;

			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;

			if (inputSampleL > 0) inputSampleL = log(1.0+(255*fabs(inputSampleL))) / log(256);
			if (inputSampleL < 0) inputSampleL = -log(1.0+(255*fabs(inputSampleL))) / log(256);

			if (inputSampleR > 0) inputSampleR = log(1.0+(255*fabs(inputSampleR))) / log(256);
			if (inputSampleR < 0) inputSampleR = -log(1.0+(255*fabs(inputSampleR))) / log(256);
		}

		// StereoChorus

		if(enableChorus) {
			chorus_cycle++;
			if (chorus_cycle == cycleEnd) { //hit the end point and we do a chorus sample
				//assign working variables
				airFactorL = airPrevL - inputSampleL;
				if (chorus_flip) {airEvenL += airFactorL; airOddL -= airFactorL; airFactorL = airEvenL;}
				else {airOddL += airFactorL; airEvenL -= airFactorL; airFactorL = airOddL;}
				airOddL = (airOddL - ((airOddL - airEvenL)/256.0)) / 1.0001;
				airEvenL = (airEvenL - ((airEvenL - airOddL)/256.0)) / 1.0001;
				airPrevL = inputSampleL;
				inputSampleL += airFactorL;
				//left
				airFactorR = airPrevR - inputSampleR;
				if (chorus_flip) {airEvenR += airFactorR; airOddR -= airFactorR; airFactorR = airEvenR;}
				else {airOddR += airFactorR; airEvenR -= airFactorR; airFactorR = airOddR;}
				airOddR = (airOddR - ((airOddR - airEvenR)/256.0)) / 1.0001;
				airEvenR = (airEvenR - ((airEvenR - airOddR)/256.0)) / 1.0001;
				airPrevR = inputSampleR;
				inputSampleR += airFactorR;
				//right
				chorus_flip = !chorus_flip;
				//air, compensates for loss of highs in flanger's interpolation
			
				int tempL = 0;
				int tempR = 0;
				if (gcount < 1 || gcount > 32760) {gcount = 32760;}
				int count = gcount;
				pL[count+32760] = pL[count] = (int)(inputSampleL*8388352.0);
				//double buffer -8388352 to 8388352 is equal to 24 bit linear space
				double offset = depth + (depth * sin(sweepL));
				count += (int)floor(offset);
				tempL += (int)(pL[count] * (1-(offset-floor(offset)))); //less as value moves away from .0
				tempL += pL[count+1]; //we can assume always using this in one way or another?
				tempL += (int)(pL[count+2] * (offset-floor(offset))); //greater as value moves away from .0
				tempL -= (int)(((pL[count]-pL[count+1])-(pL[count+1]-pL[count+2]))/50); //interpolation hacks 'r us
				//left
			
				count = gcount;
				pR[count+32760] = pR[count] = (int)(inputSampleR*8388352.0);
				//double buffer -8388352 to 8388352 is equal to 24 bit linear space
				offset = depth + (depth * sin(sweepR));
				count += (int)floor(offset);
				tempR += (int)(pR[count] * (1-(offset-floor(offset)))); //less as value moves away from .0
				tempR += pR[count+1]; //we can assume always using this in one way or another?
				tempR += (int)(pR[count+2] * (offset-floor(offset))); //greater as value moves away from .0
				tempR -= (int)(((pR[count]-pR[count+1])-(pR[count+1]-pR[count+2]))/50); //interpolation hacks 'r us
				//right
			
				sweepL += speed;
				sweepR += speed;
				if (sweepL > tupi){sweepL -= tupi;}
				if (sweepR > tupi){sweepR -= tupi;}
				gcount--;
				//still scrolling through the samples, remember
			
				inputSampleL = ((double)(tempL/16776704.0));
				inputSampleR = ((double)(tempR/16776704.0));
				if (cycleEnd == 4) {
					lastRefL[0] = lastRefL[4]; //start from previous last
					lastRefL[2] = (lastRefL[0] + inputSampleL)/2; //half
					lastRefL[1] = (lastRefL[0] + lastRefL[2])/2; //one quarter
					lastRefL[3] = (lastRefL[2] + inputSampleL)/2; //three quarters
					lastRefL[4] = inputSampleL; //full
					lastRefR[0] = lastRefR[4]; //start from previous last
					lastRefR[2] = (lastRefR[0] + inputSampleR)/2; //half
					lastRefR[1] = (lastRefR[0] + lastRefR[2])/2; //one quarter
					lastRefR[3] = (lastRefR[2] + inputSampleR)/2; //three quarters
					lastRefR[4] = inputSampleR; //full
				}
				if (cycleEnd == 3) {
					lastRefL[0] = lastRefL[3]; //start from previous last
					lastRefL[2] = (lastRefL[0]+lastRefL[0]+inputSampleL)/3; //third
					lastRefL[1] = (lastRefL[0]+inputSampleL+inputSampleL)/3; //two thirds
					lastRefL[3] = inputSampleL; //full
					lastRefR[0] = lastRefR[3]; //start from previous last
					lastRefR[2] = (lastRefR[0]+lastRefR[0]+inputSampleR)/3; //third
					lastRefR[1] = (lastRefR[0]+inputSampleR+inputSampleR)/3; //two thirds
					lastRefR[3] = inputSampleR; //full
				}
				if (cycleEnd == 2) {
					lastRefL[0] = lastRefL[2]; //start from previous last
					lastRefL[1] = (lastRefL[0] + inputSampleL)/2; //half
					lastRefL[2] = inputSampleL; //full
					lastRefR[0] = lastRefR[2]; //start from previous last
					lastRefR[1] = (lastRefR[0] + inputSampleR)/2; //half
					lastRefR[2] = inputSampleR; //full
				}
				if (cycleEnd == 1) {
					lastRefL[0] = inputSampleL;
					lastRefR[0] = inputSampleR;
				}
				chorus_cycle = 0; //reset
				inputSampleL = lastRefL[chorus_cycle];
				inputSampleR = lastRefR[chorus_cycle];
			} else {
				inputSampleL = lastRefL[chorus_cycle];
				inputSampleR = lastRefR[chorus_cycle];
				//we are going through our references now
			}
		}

		// ClearCoat
		
		if(enableReverb) {
			double cc_drySampleL = inputSampleL;
			double cc_drySampleR = inputSampleR;
		
			cc_cycle++;
			if (cc_cycle == cycleEnd) { //hit the end point and we do a reverb sample
				aAL[countAL] = inputSampleL + (feedbackAL * 0.04166666666);
				aBL[countBL] = inputSampleL + (feedbackBL * 0.04166666666);
				aCL[countCL] = inputSampleL + (feedbackCL * 0.04166666666);
				aDL[countDL] = inputSampleL + (feedbackDL * 0.04166666666);
			
				aDR[countDR] = inputSampleR + (feedbackDR * 0.04166666666);
				aHR[countHR] = inputSampleR + (feedbackHR * 0.04166666666);
				aLR[countLR] = inputSampleR + (feedbackLR * 0.04166666666);
				aPR[countPR] = inputSampleR + (feedbackPR * 0.04166666666);
				//exactly halfway between infinite sustain at 0.0625
				//and 6dB down, almost no regen at 0.03125
				//means roughly half the results work as BitShiftGain
			
				countAL++; if (countAL < 0 || countAL > shortA) countAL = 0;
				countBL++; if (countBL < 0 || countBL > shortB) countBL = 0;
				countCL++; if (countCL < 0 || countCL > shortC) countCL = 0;
				countDL++; if (countDL < 0 || countDL > shortD) countDL = 0;
			
				countDR++; if (countDR < 0 || countDR > shortD) countDR = 0;
				countHR++; if (countHR < 0 || countHR > shortH) countHR = 0;
				countLR++; if (countLR < 0 || countLR > shortL) countLR = 0;
				countPR++; if (countPR < 0 || countPR > shortP) countPR = 0;
			
				double outAL = aAL[countAL-((countAL > shortA)?shortA+1:0)];
				double outBL = aBL[countBL-((countBL > shortB)?shortB+1:0)];
				double outCL = aCL[countCL-((countCL > shortC)?shortC+1:0)];
				double outDL = aDL[countDL-((countDL > shortD)?shortD+1:0)];
			
				double outDR = aDR[countDR-((countDR > shortD)?shortD+1:0)];
				double outHR = aHR[countHR-((countHR > shortH)?shortH+1:0)];
				double outLR = aLR[countLR-((countLR > shortL)?shortL+1:0)];
				double outPR = aPR[countPR-((countPR > shortP)?shortP+1:0)];
			
				aEL[countEL] = outAL - (outBL + outCL + outDL);
				aFL[countFL] = outBL - (outAL + outCL + outDL);
				aGL[countGL] = outCL - (outAL + outBL + outDL);
				aHL[countHL] = outDL - (outAL + outBL + outCL);
			
				aCR[countCR] = outDR - (outHR + outLR + outPR);
				aGR[countGR] = outHR - (outDR + outLR + outPR);
				aKR[countKR] = outLR - (outDR + outHR + outPR);
				aOR[countOR] = outPR - (outDR + outHR + outLR);
			
				countEL++; if (countEL < 0 || countEL > shortE) countEL = 0;
				countFL++; if (countFL < 0 || countFL > shortF) countFL = 0;
				countGL++; if (countGL < 0 || countGL > shortG) countGL = 0;
				countHL++; if (countHL < 0 || countHL > shortH) countHL = 0;
			
				countCR++; if (countCR < 0 || countCR > shortC) countCR = 0;
				countGR++; if (countGR < 0 || countGR > shortG) countGR = 0;
				countKR++; if (countKR < 0 || countKR > shortK) countKR = 0;
				countOR++; if (countOR < 0 || countOR > shortO) countOR = 0;
			
				double outEL = aEL[countEL-((countEL > shortE)?shortE+1:0)];
				double outFL = aFL[countFL-((countFL > shortF)?shortF+1:0)];
				double outGL = aGL[countGL-((countGL > shortG)?shortG+1:0)];
				double outHL = aHL[countHL-((countHL > shortH)?shortH+1:0)];
			
				double outCR = aCR[countCR-((countCR > shortC)?shortC+1:0)];
				double outGR = aGR[countGR-((countGR > shortG)?shortG+1:0)];
				double outKR = aKR[countKR-((countKR > shortK)?shortK+1:0)];
				double outOR = aOR[countOR-((countOR > shortO)?shortO+1:0)];
			
				aIL[countIL] = outEL - (outFL + outGL + outHL);
				aJL[countJL] = outFL - (outEL + outGL + outHL);
				aKL[countKL] = outGL - (outEL + outFL + outHL);
				aLL[countLL] = outHL - (outEL + outFL + outGL);
			
				aBR[countBR] = outCR - (outGR + outKR + outOR);
				aFR[countFR] = outGR - (outCR + outKR + outOR);
				aJR[countJR] = outKR - (outCR + outGR + outOR);
				aNR[countNR] = outOR - (outCR + outGR + outKR);
			
				countIL++; if (countIL < 0 || countIL > shortI) countIL = 0;
				countJL++; if (countJL < 0 || countJL > shortJ) countJL = 0;
				countKL++; if (countKL < 0 || countKL > shortK) countKL = 0;
				countLL++; if (countLL < 0 || countLL > shortL) countLL = 0;
			
				countBR++; if (countBR < 0 || countBR > shortB) countBR = 0;
				countFR++; if (countFR < 0 || countFR > shortF) countFR = 0;
				countJR++; if (countJR < 0 || countJR > shortJ) countJR = 0;
				countNR++; if (countNR < 0 || countNR > shortN) countNR = 0;
			
				double outIL = aIL[countIL-((countIL > shortI)?shortI+1:0)];
				double outJL = aJL[countJL-((countJL > shortJ)?shortJ+1:0)];
				double outKL = aKL[countKL-((countKL > shortK)?shortK+1:0)];
				double outLL = aLL[countLL-((countLL > shortL)?shortL+1:0)];
			
				double outBR = aBR[countBR-((countBR > shortB)?shortB+1:0)];
				double outFR = aFR[countFR-((countFR > shortF)?shortF+1:0)];
				double outJR = aJR[countJR-((countJR > shortJ)?shortJ+1:0)];
				double outNR = aNR[countNR-((countNR > shortN)?shortN+1:0)];
			
				aML[countML] = outIL - (outJL + outKL + outLL);
				aNL[countNL] = outJL - (outIL + outKL + outLL);
				aOL[countOL] = outKL - (outIL + outJL + outLL);
				aPL[countPL] = outLL - (outIL + outJL + outKL);
			
				aAR[countAR] = outBR - (outFR + outJR + outNR);
				aER[countER] = outFR - (outBR + outJR + outNR);
				aIR[countIR] = outJR - (outBR + outFR + outNR);
				aMR[countMR] = outNR - (outBR + outFR + outJR);
			
				countML++; if (countML < 0 || countML > shortM) countML = 0;
				countNL++; if (countNL < 0 || countNL > shortN) countNL = 0;
				countOL++; if (countOL < 0 || countOL > shortO) countOL = 0;
				countPL++; if (countPL < 0 || countPL > shortP) countPL = 0;
			
				countAR++; if (countAR < 0 || countAR > shortA) countAR = 0;
				countER++; if (countER < 0 || countER > shortE) countER = 0;
				countIR++; if (countIR < 0 || countIR > shortI) countIR = 0;
				countMR++; if (countMR < 0 || countMR > shortM) countMR = 0;
			
				double outML = aML[countML-((countML > shortM)?shortM+1:0)];
				double outNL = aNL[countNL-((countNL > shortN)?shortN+1:0)];
				double outOL = aOL[countOL-((countOL > shortO)?shortO+1:0)];
				double outPL = aPL[countPL-((countPL > shortP)?shortP+1:0)];
			
				double outAR = aAR[countAR-((countAR > shortA)?shortA+1:0)];
				double outER = aER[countER-((countER > shortE)?shortE+1:0)];
				double outIR = aIR[countIR-((countIR > shortI)?shortI+1:0)];
				double outMR = aMR[countMR-((countMR > shortM)?shortM+1:0)];			
				double outSample = (outML + outML + outML + prevMulchAL)*0.25;
				prevMulchAL = outML; outML = outSample;
				outSample = (outAR + outAR + outAR + prevMulchAR)*0.25;
				prevMulchAR = outAR; outAR = outSample;
			
				feedbackAL = outML - (outNL + outOL + outPL);
				feedbackDR = outAR - (outER + outIR + outMR);
				feedbackBL = outNL - (outML + outOL + outPL);
				feedbackHR = outER - (outAR + outIR + outMR);
				feedbackCL = outOL - (outML + outNL + outPL);
				feedbackLR = outIR - (outAR + outER + outMR);
				feedbackDL = outPL - (outML + outNL + outOL);
				feedbackPR = outMR - (outAR + outER + outIR);
				//which we need to feed back into the input again, a bit
			
				inputSampleL = (outML + outNL + outOL + outPL)/8.0;
				inputSampleR = (outAR + outER + outIR + outMR)/8.0;
				//and take the final combined sum of outputs, corrected for Householder gain
			
				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			
				if (cycleEnd == 4) {
					cc_LastRefL[0] = cc_LastRefL[4]; //start from previous last
					cc_LastRefL[2] = (cc_LastRefL[0] + inputSampleL)/2; //half
					cc_LastRefL[1] = (cc_LastRefL[0] + cc_LastRefL[2])/2; //one quarter
					cc_LastRefL[3] = (cc_LastRefL[2] + inputSampleL)/2; //three quarters
					cc_LastRefL[4] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[4]; //start from previous last
					cc_LastRefR[2] = (cc_LastRefR[0] + inputSampleR)/2; //half
					cc_LastRefR[1] = (cc_LastRefR[0] + cc_LastRefR[2])/2; //one quarter
					cc_LastRefR[3] = (cc_LastRefR[2] + inputSampleR)/2; //three quarters
					cc_LastRefR[4] = inputSampleR; //full
				}
				if (cycleEnd == 3) {
					cc_LastRefL[0] = cc_LastRefL[3]; //start from previous last
					cc_LastRefL[2] = (cc_LastRefL[0]+cc_LastRefL[0]+inputSampleL)/3; //third
					cc_LastRefL[1] = (cc_LastRefL[0]+inputSampleL+inputSampleL)/3; //two thirds
					cc_LastRefL[3] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[3]; //start from previous last
					cc_LastRefR[2] = (cc_LastRefR[0]+cc_LastRefR[0]+inputSampleR)/3; //third
					cc_LastRefR[1] = (cc_LastRefR[0]+inputSampleR+inputSampleR)/3; //two thirds
					cc_LastRefR[3] = inputSampleR; //full
				}
				if (cycleEnd == 2) {
					cc_LastRefL[0] = cc_LastRefL[2]; //start from previous last
					cc_LastRefL[1] = (cc_LastRefL[0] + inputSampleL)/2; //half
					cc_LastRefL[2] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[2]; //start from previous last
					cc_LastRefR[1] = (cc_LastRefR[0] + inputSampleR)/2; //half
					cc_LastRefR[2] = inputSampleR; //full
				}
				if (cycleEnd == 1) {
					cc_LastRefL[0] = inputSampleL;
					cc_LastRefR[0] = inputSampleR;
				}
				cc_cycle = 0; //reset
				inputSampleL = cc_LastRefL[cc_cycle];
				inputSampleR = cc_LastRefR[cc_cycle];
			} else {
				inputSampleL = cc_LastRefL[cc_cycle];
				inputSampleR = cc_LastRefR[cc_cycle];
				//we are going through our references now
			}
		
			//begin SubTight section
			double subSampleL = inputSampleL * cc_subRate;
			double subSampleR = inputSampleR * cc_subRate;
		
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
			scale = 0.5+fabs(subSampleL*0.5);
			subSampleL = (subDL+(sin(subDL-subSampleL)*scale));
			subDL = subSampleL*scale;
			scale = 0.5+fabs(subSampleR*0.5);
			subSampleR = (subDR+(sin(subDR-subSampleR)*scale));
			subDR = subSampleR*scale;
			if (subSampleL > 0.25) subSampleL = 0.25;
			if (subSampleL < -0.25) subSampleL = -0.25;
			if (subSampleR > 0.25) subSampleR = 0.25;
			if (subSampleR < -0.25) subSampleR = -0.25;
			inputSampleL -= (subSampleL*16.0);
			inputSampleR -= (subSampleR*16.0);
			//end SubTight section		
		
			if (cycleEnd > 1) {
				temp = (inputSampleL + tailL)*0.5;
				tailL = inputSampleL; inputSampleL = temp;
				temp = (inputSampleR + tailR)*0.5;
				tailR = inputSampleR; inputSampleR = temp;
			} //let's average only at elevated sample rates
		
			if (cc_wet < 1.0) {inputSampleL *= cc_wet; inputSampleR *= cc_wet;}
			if (cc_dry < 1.0) {cc_drySampleL *= cc_dry; cc_drySampleR *= cc_dry;}
			inputSampleL += cc_drySampleL; inputSampleR += cc_drySampleR;
			//this is our submix verb dry/wet: 0.5 is BOTH at FULL VOLUME
			//purpose is that, if you're adding verb, you're not altering other balances
		}

		inputSampleL *= boost; //now, focussed area gets cranked before distort
		inputSampleR *= boost; //now, focussed area gets cranked before distort
		
		switch (mode)
		{
			case 0: //Density
				if (inputSampleL > 1.570796326794897) inputSampleL = 1.570796326794897;
				if (inputSampleL < -1.570796326794897) inputSampleL = -1.570796326794897;
				if (inputSampleR > 1.570796326794897) inputSampleR = 1.570796326794897;
				if (inputSampleR < -1.570796326794897) inputSampleR = -1.570796326794897;
				//clip to 1.570796326794897 to reach maximum output
				inputSampleL = sin(inputSampleL);
				inputSampleR = sin(inputSampleR);
				break;
			case 1: //Drive				
				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				inputSampleL -= (inputSampleL * (fabs(inputSampleL) * 0.6) * (fabs(inputSampleL) * 0.6));
				inputSampleR -= (inputSampleR * (fabs(inputSampleR) * 0.6) * (fabs(inputSampleR) * 0.6));
				inputSampleL *= 1.6;
				inputSampleR *= 1.6;
				break;
			case 2: //Spiral
				if (inputSampleL > 1.2533141373155) inputSampleL = 1.2533141373155;
				if (inputSampleL < -1.2533141373155) inputSampleL = -1.2533141373155;
				if (inputSampleR > 1.2533141373155) inputSampleR = 1.2533141373155;
				if (inputSampleR < -1.2533141373155) inputSampleR = -1.2533141373155;
				//clip to 1.2533141373155 to reach maximum output
				inputSampleL = sin(inputSampleL * fabs(inputSampleL)) / ((fabs(inputSampleL) == 0.0) ?1:fabs(inputSampleL));
				inputSampleR = sin(inputSampleR * fabs(inputSampleR)) / ((fabs(inputSampleR) == 0.0) ?1:fabs(inputSampleR));
				break;
			case 3: //Mojo
				temp = pow(fabs(inputSampleL),0.25);
				if (temp > 0.0) inputSampleL = (sin(inputSampleL * temp * M_PI * 0.5) / temp) * 0.987654321;
				temp = pow(fabs(inputSampleR),0.25);
				if (temp > 0.0) inputSampleR = (sin(inputSampleR * temp * M_PI * 0.5) / temp) * 0.987654321;
				//mojo is the one that flattens WAAAAY out very softly before wavefolding				
				break;
			case 4: //Dyno
				temp = pow(fabs(inputSampleL),4);
				if (temp > 0.0) inputSampleL = (sin(inputSampleL * temp) / temp) * 1.1654321;
				temp = pow(fabs(inputSampleR),4);
				if (temp > 0.0) inputSampleR = (sin(inputSampleR * temp) / temp) * 1.1654321;
				//dyno is the one that tries to raise peak energy				
				break;
		}				
		
		if (output != 1.0) {
			inputSampleL *= output;
			inputSampleR *= output;
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

		inputSampleL += groundSampleL; //effectively UnBox
		inputSampleR += groundSampleR; //effectively UnBox
		
		if (wet !=1.0) {
			inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0-wet));
			inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0-wet));
		}


		// Creature

		if(enableCreature) {
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

			inputSampleL *= cr_wet;
			inputSampleR *= cr_wet;
			drySampleL *= cr_dry;
			drySampleR *= cr_dry;
			inputSampleL += drySampleL;
			inputSampleR += drySampleR;
		}

		// ToneSlant
		if(enableToneSlant) {
			for (int count = ts_overallscale; count >= 0; count--) {
				ts_bL[count+1] = ts_bL[count];
				ts_bR[count+1] = ts_bR[count];
			}

			ts_bL[0] = ts_accumulatorSampleL = inputSampleL;
			ts_bR[0] = ts_accumulatorSampleR = inputSampleR;

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
		}

		// ButterComp2
		
		if(enableComp) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;
		
			inputSampleL *= bc_inputgain;
			inputSampleR *= bc_inputgain;
		
			double divisor = bc_compfactor / (1.0+fabs(lastOutputL));
			//this is slowing compressor recovery while output waveforms were high
			divisor /= overallscale;
			double remainder = divisor;
			divisor = 1.0 - divisor;
			//recalculate divisor every sample		
		
			double inputposL = inputSampleL + 1.0;
			if (inputposL < 0.0) inputposL = 0.0;
			double outputposL = inputposL / 2.0;
			if (outputposL > 1.0) outputposL = 1.0;		
			inputposL *= inputposL;
			targetposL *= divisor;
			targetposL += (inputposL * remainder);
			double calcposL = pow((1.0/targetposL),2);
		
			double inputnegL = (-inputSampleL) + 1.0;
			if (inputnegL < 0.0) inputnegL = 0.0;
			double outputnegL = inputnegL / 2.0;
			if (outputnegL > 1.0) outputnegL = 1.0;		
			inputnegL *= inputnegL;
			targetnegL *= divisor;
			targetnegL += (inputnegL * remainder);
			double calcnegL = pow((1.0/targetnegL),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1
		
			if (inputSampleL > 0)
			{ //working on pos
				if (flip)
				{
					controlAposL *= divisor;
					controlAposL += (calcposL*remainder);
				
				}
				else
				{
					controlBposL *= divisor;
					controlBposL += (calcposL*remainder);
				}
			}
			else
			{ //working on neg
				if (flip)
				{
					controlAnegL *= divisor;
					controlAnegL += (calcnegL*remainder);
				}
				else
				{
					controlBnegL *= divisor;
					controlBnegL += (calcnegL*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'
		
			divisor = bc_compfactor / (1.0+fabs(lastOutputR));
			//this is slowing compressor recovery while output waveforms were high
			divisor /= overallscale;
			remainder = divisor;
			divisor = 1.0 - divisor;
			//recalculate divisor every sample		
		
			double inputposR = inputSampleR + 1.0;
			if (inputposR < 0.0) inputposR = 0.0;
			double outputposR = inputposR / 2.0;
			if (outputposR > 1.0) outputposR = 1.0;		
			inputposR *= inputposR;
			targetposR *= divisor;
			targetposR += (inputposR * remainder);
			double calcposR = pow((1.0/targetposR),2);
		
			double inputnegR = (-inputSampleR) + 1.0;
			if (inputnegR < 0.0) inputnegR = 0.0;
			double outputnegR = inputnegR / 2.0;
			if (outputnegR > 1.0) outputnegR = 1.0;		
			inputnegR *= inputnegR;
			targetnegR *= divisor;
			targetnegR += (inputnegR * remainder);
			double calcnegR = pow((1.0/targetnegR),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1
		
			if (inputSampleR > 0)
			{ //working on pos
				if (flip)
				{
					controlAposR *= divisor;
					controlAposR += (calcposR*remainder);
				
				}
				else
				{
					controlBposR *= divisor;
					controlBposR += (calcposR*remainder);
				}
			}
			else
			{ //working on neg
				if (flip)
				{
					controlAnegR *= divisor;
					controlAnegR += (calcnegR*remainder);
				}
				else
				{
					controlBnegR *= divisor;
					controlBnegR += (calcnegR*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'
		
			double totalmultiplierL;
			double totalmultiplierR;
			if (flip)
			{
				totalmultiplierL = (controlAposL * outputposL) + (controlAnegL * outputnegL);
				totalmultiplierR = (controlAposR * outputposR) + (controlAnegR * outputnegR);
			}
			else
			{
				totalmultiplierL = (controlBposL * outputposL) + (controlBnegL * outputnegL);
				totalmultiplierR = (controlBposR * outputposR) + (controlBnegR * outputnegR);
			}
			//this combines the sides according to flip, blending relative to the input value
		
			inputSampleL *= totalmultiplierL;
			inputSampleL /= bc_outputgain;
		
			inputSampleR *= totalmultiplierR;
			inputSampleR /= bc_outputgain;
		
			if (bc_output != 1.0) {
				inputSampleL *= bc_output;
				inputSampleR *= bc_output;
			}

			if (bc_wet !=1.0) {
				inputSampleL = (inputSampleL * bc_wet) + (drySampleL * (1.0-bc_wet));
				inputSampleR = (inputSampleR * bc_wet) + (drySampleR * (1.0-bc_wet));
			}
		
			lastOutputL = inputSampleL;
			lastOutputR = inputSampleR;
			//we will make this factor respond to use of dry/wet
		}

		// Gain trim

		inputSampleL = inputSampleL * gain_trim;
		inputSampleR = inputSampleR * gain_trim;

		// ADClip

		double softness = 0.618033988749894848204586;

		if (adcLastSampleL >= 0.5)
		{
			if (inputSampleL < 0.5) adcLastSampleL = ((0.5*softness) + (inputSampleL * (1.0-softness)));
			else adcLastSampleL = 0.5;
		}

		if (adcLastSampleL <= -0.5)
		{
			if (inputSampleL > -0.5) adcLastSampleL = ((-0.5*softness) + (inputSampleL * (1.0-softness)));
			else adcLastSampleL = -0.5;
		}

		if (inputSampleL > 0.5)
		{
			if (adcLastSampleL < 0.5) inputSampleL = ((0.5*softness) + (adcLastSampleL * (1.0-softness)));
			else inputSampleL = 0.5;
		}

		if (inputSampleL < -0.5)
		{
			if (adcLastSampleL > -0.5) inputSampleL = ((-0.5*softness) + (adcLastSampleL * (1.0-softness)));
			else inputSampleL = -0.5;
		}
		adcLastSampleL = inputSampleL; //end ADClip L


		if (adcLastSampleR >= 0.5)
		{
			if (inputSampleR < 0.5) adcLastSampleR = ((0.5*softness) + (inputSampleR * (1.0-softness)));
			else adcLastSampleR = 0.5;
		}

		if (adcLastSampleR <= -0.5)
		{
			if (inputSampleR > -0.5) adcLastSampleR = ((-0.5*softness) + (inputSampleR * (1.0-softness)));
			else adcLastSampleR = -0.5;
		}

		if (inputSampleR > 0.5)
		{
			if (adcLastSampleR < 0.5) inputSampleR = ((0.5*softness) + (adcLastSampleR * (1.0-softness)));
			else inputSampleR = 0.5;
		}

		if (inputSampleR < -0.5)
		{
			if (adcLastSampleR > -0.5) inputSampleR = ((-0.5*softness) + (adcLastSampleR * (1.0-softness)));
			else inputSampleR = -0.5;
		}
		adcLastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.5) inputSampleL = 0.5;
		if (inputSampleL < -0.5) inputSampleL = -0.5;
		//final iron bar
		if (inputSampleR > 0.5) inputSampleR = 0.5;
		if (inputSampleR < -0.5) inputSampleR = -0.5;
		//final iron bar

		// Channel9

		//begin L
		temp = (ch9_lastSampleBL - ch9_lastSampleCL) * 0.381966011250105;
		temp -= (ch9_lastSampleAL - ch9_lastSampleBL) * 0.6180339887498948482045;
		temp += inputSampleL - ch9_lastSampleAL; //regular slew clamping added
		
		ch9_lastSampleCL = ch9_lastSampleBL;
		ch9_lastSampleBL = ch9_lastSampleAL;
		ch9_lastSampleAL = inputSampleL; //now our output relates off ch9_lastSampleB
		
		if (temp > localthreshold)
			inputSampleL = ch9_lastSampleBL + localthreshold;
		if (-temp > localthreshold)
			inputSampleL = ch9_lastSampleBL - localthreshold;
		
		ch9_lastSampleAL = (ch9_lastSampleAL*0.381966011250105)+(inputSampleL*0.6180339887498948482045); //split the difference between raw and smoothed for buffer
		//end L
		
		//begin R
		temp = (ch9_lastSampleBR - ch9_lastSampleCR) * 0.381966011250105;
		temp -= (ch9_lastSampleAR - ch9_lastSampleBR) * 0.6180339887498948482045;
		temp += inputSampleR - ch9_lastSampleAR; //regular slew clamping added
		
		ch9_lastSampleCR = ch9_lastSampleBR;
		ch9_lastSampleBR = ch9_lastSampleAR;
		ch9_lastSampleAR = inputSampleR; //now our output relates off ch9_lastSampleB
		
		if (temp > localthreshold)
			inputSampleR = ch9_lastSampleBR + localthreshold;
		if (-temp > localthreshold)
			inputSampleR = ch9_lastSampleBR - localthreshold;
		
		ch9_lastSampleAR = (ch9_lastSampleAR*0.381966011250105)+(inputSampleR*0.6180339887498948482045); //split the difference between raw and smoothed for buffer
		//end R
		
		temp = ch9_biquadB[2]*inputSampleL+ch9_biquadB[3]*ch9_biquadB[7]+ch9_biquadB[4]*ch9_biquadB[8]-ch9_biquadB[5]*ch9_biquadB[9]-ch9_biquadB[6]*ch9_biquadB[10];
		ch9_biquadB[8] = ch9_biquadB[7]; ch9_biquadB[7] = inputSampleL; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleL = temp;
		ch9_biquadB[10] = ch9_biquadB[9]; ch9_biquadB[9] = inputSampleL; //DF1 left
		temp = ch9_biquadB[2]*inputSampleR+ch9_biquadB[3]*ch9_biquadB[11]+ch9_biquadB[4]*ch9_biquadB[12]-ch9_biquadB[5]*ch9_biquadB[13]-ch9_biquadB[6]*ch9_biquadB[14];
		ch9_biquadB[12] = ch9_biquadB[11]; ch9_biquadB[11] = inputSampleR; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleR = temp;
		ch9_biquadB[14] = ch9_biquadB[13]; ch9_biquadB[13] = inputSampleR; //DF1 right

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

void ConsoleZTracking::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
{
  double* in1  =  inputs[0];
  double* in2  =  inputs[1];
  double* out1 = outputs[0];
  double* out2 = outputs[1];

	bool enableHPF = B != 0.0;
	bool enableLPF = D != 1.0;
	bool enableComp = F != 0.0;
	bool enableChorus = J != 0.0;
	bool enableReverb = L != 0.0;
	bool enableuLaw = Q != 0.0;
	bool enableCreature = V != 0.5;
	bool enableToneSlant = X != 0.5;

	double inputSampleL;
	double inputSampleR;
	double drySampleL;
	double drySampleR;
	double temp;	

	double sampleRate = getSampleRate() * 2.0;
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= sampleRate;

	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;
	//this is going to be 2 for 88.1 or 96k, 3 for silly people, 4 for 176 or 192k

	double outSample;
	double KK;
	double norm;

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

	if (cc_cycle > cycleEnd-1) cc_cycle = cycleEnd-1; //sanity check
	
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
	ch9_biquadB[0] = ch9_biquadA[0] = ch9_cutoff / sampleRate;
    ch9_biquadA[1] = 1.618033988749894848204586;
	ch9_biquadB[1] = 0.618033988749894848204586;
	
	KK = tan(M_PI * ch9_biquadA[0]); //lowpass
	norm = 1.0 / (1.0 + KK / ch9_biquadA[1] + KK * KK);
	ch9_biquadA[2] = KK * KK * norm;
	ch9_biquadA[3] = 2.0 * ch9_biquadA[2];
	ch9_biquadA[4] = ch9_biquadA[2];
	ch9_biquadA[5] = 2.0 * (KK * KK - 1.0) * norm;
	ch9_biquadA[6] = (1.0 - KK / ch9_biquadA[1] + KK * KK) * norm;
	
	KK = tan(M_PI * ch9_biquadA[0]);
	norm = 1.0 / (1.0 + KK / ch9_biquadB[1] + KK * KK);
	ch9_biquadB[2] = KK * KK * norm;
	ch9_biquadB[3] = 2.0 * ch9_biquadB[2];
	ch9_biquadB[4] = ch9_biquadB[2];
	ch9_biquadB[5] = 2.0 * (KK * KK - 1.0) * norm;
	ch9_biquadB[6] = (1.0 - KK / ch9_biquadB[1] + KK * KK) * norm;

	// Focus

	double boost = pow(10.0,(M*12.0)/20.0);
	figureL[0] = figureR[0] = (pow(M_E, freqX * O) * centerFreq)/sampleRate; //fixed frequency, 3.515775k
	figureL[1] = figureR[1] = pow(pow(N,3)*2,2)+0.0001; //resonance
	int mode = (int) ( P * 4.999 );
	double output = R;
	double wet = S;
	
	KK = tan(M_PI * figureR[0]);
	norm = 1.0 / (1.0 + KK / figureR[1] + KK * KK);
	figureL[2] = figureR[2] = KK / figureR[1] * norm;
	figureL[4] = figureR[4] = -figureR[2];
	figureL[5] = figureR[5] = 2.0 * (KK * KK - 1.0) * norm;
	figureL[6] = figureR[6] = (1.0 - KK / figureR[1] + KK * KK) * norm;

	// Creature

	double source = 1.0-pow(1.0-T,5);
	int stages = (pow(U,2)*32.0*sqrt(overallscale))+1;
	double cr_wet = (V*2.0)-1.0; //inv-dry-wet for highpass
	double cr_dry = 2.0-(V*2.0);
	if (cr_dry > 1.0) cr_dry = 1.0; //full dry for use with inv, to 0.0 at full wet

	// ToneSlant

	double ts_correctionSampleL;
	double ts_correctionSampleR;
	double ts_accumulatorSampleL;
	double ts_accumulatorSampleR;
	double ts_overallscale = (W*99.0)+1.0;
	double ts_applySlant = (X*2.0)-1.0;


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

	// +/- 12dB gain trim
	double gain_trim = pow(M_E, freqX * Y) * 0.25;

	while (sampleFrames > 0)
  {

		if(flip) {
			inputSampleL = fpdL * 1.18e-17;
			inputSampleR = fpdR * 1.18e-17;
		}
		else {
			inputSampleL = *in1 * input_gain * 2.0;
			inputSampleR = *in2 * input_gain * 2.0;

			if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
			if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;
		}

		// Channel9
		temp = ch9_biquadA[2]*inputSampleL+ch9_biquadA[3]*ch9_biquadA[7]+ch9_biquadA[4]*ch9_biquadA[8]-ch9_biquadA[5]*ch9_biquadA[9]-ch9_biquadA[6]*ch9_biquadA[10];
		ch9_biquadA[8] = ch9_biquadA[7]; ch9_biquadA[7] = inputSampleL; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleL = temp;
		ch9_biquadA[10] = ch9_biquadA[9]; ch9_biquadA[9] = inputSampleL; //DF1 left
		temp = ch9_biquadA[2]*inputSampleR+ch9_biquadA[3]*ch9_biquadA[11]+ch9_biquadA[4]*ch9_biquadA[12]-ch9_biquadA[5]*ch9_biquadA[13]-ch9_biquadA[6]*ch9_biquadA[14];
		ch9_biquadA[12] = ch9_biquadA[11]; ch9_biquadA[11] = inputSampleR; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleR = temp;
		ch9_biquadA[14] = ch9_biquadA[13]; ch9_biquadA[13] = inputSampleR; //DF1 right

		//inputSampleL = sin(inputSampleL);
		inputSampleL += ((pow(inputSampleL,5)/128.0) + (pow(inputSampleL,9)/262144.0)) - ((pow(inputSampleL,3)/8.0) + (pow(inputSampleL,7)/4096.0));
		//inputSampleR = sin(inputSampleR);
		inputSampleR += ((pow(inputSampleR,5)/128.0) + (pow(inputSampleR,9)/262144.0)) - ((pow(inputSampleR,3)/8.0) + (pow(inputSampleR,7)/4096.0));

		double dielectricScaleL = fabs(2.0-((inputSampleL+nonLin)/nonLin));
		double dielectricScaleR = fabs(2.0-((inputSampleR+nonLin)/nonLin));

		if (flip)
		{

			if (fabs(ch9_iirSampleLA)<1.18e-37) ch9_iirSampleLA = 0.0; 
			ch9_iirSampleLA = (ch9_iirSampleLA * (1.0 - (ch9_localiirAmount * dielectricScaleL))) + (inputSampleL * ch9_localiirAmount * dielectricScaleL);
			inputSampleL = inputSampleL - ch9_iirSampleLA;
			if (fabs(ch9_iirSampleRA)<1.18e-37) ch9_iirSampleRA = 0.0; 
			ch9_iirSampleRA = (ch9_iirSampleRA * (1.0 - (ch9_localiirAmount * dielectricScaleR))) + (inputSampleR * ch9_localiirAmount * dielectricScaleR);
			inputSampleR = inputSampleR - ch9_iirSampleRA;

			if(enableHPF) {
				temp = (inputSampleL * biquad_hpf_AL[2]) + biquad_hpf_AL[7];
				biquad_hpf_AL[7] = (inputSampleL * biquad_hpf_AL[3]) - (temp * biquad_hpf_AL[5]) + biquad_hpf_AL[8];
				biquad_hpf_AL[8] = (inputSampleL * biquad_hpf_AL[4]) - (temp * biquad_hpf_AL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_hpf_AR[2]) + biquad_hpf_AR[7];
				biquad_hpf_AR[7] = (inputSampleR * biquad_hpf_AR[3]) - (temp * biquad_hpf_AR[5]) + biquad_hpf_AR[8];
				biquad_hpf_AR[8] = (inputSampleR * biquad_hpf_AR[4]) - (temp * biquad_hpf_AR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}

			if(enableLPF) {
				temp = (inputSampleL * biquad_lpf_AL[2]) + biquad_lpf_AL[7];
				biquad_lpf_AL[7] = (inputSampleL * biquad_lpf_AL[3]) - (temp * biquad_lpf_AL[5]) + biquad_lpf_AL[8];
				biquad_lpf_AL[8] = (inputSampleL * biquad_lpf_AL[4]) - (temp * biquad_lpf_AL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_lpf_AR[2]) + biquad_lpf_AR[7];
				biquad_lpf_AR[7] = (inputSampleR * biquad_lpf_AR[3]) - (temp * biquad_lpf_AR[5]) + biquad_lpf_AR[8];
				biquad_lpf_AR[8] = (inputSampleR * biquad_lpf_AR[4]) - (temp * biquad_lpf_AR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}
		} else {

			if (fabs(ch9_iirSampleLB)<1.18e-37) ch9_iirSampleLB = 0.0; 
			ch9_iirSampleLB = (ch9_iirSampleLB * (1.0 - (ch9_localiirAmount * dielectricScaleL))) + (inputSampleL * ch9_localiirAmount * dielectricScaleL);
			inputSampleL = inputSampleL - ch9_iirSampleLB;
			if (fabs(ch9_iirSampleRB)<1.18e-37) ch9_iirSampleRB = 0.0; 
			ch9_iirSampleRB = (ch9_iirSampleRB * (1.0 - (ch9_localiirAmount * dielectricScaleR))) + (inputSampleR * ch9_localiirAmount * dielectricScaleR);
			inputSampleR = inputSampleR - ch9_iirSampleRB;

			if(enableHPF) {
				temp = (inputSampleL * biquad_hpf_BL[2]) + biquad_hpf_BL[7];
				biquad_hpf_BL[7] = (inputSampleL * biquad_hpf_BL[3]) - (temp * biquad_hpf_BL[5]) + biquad_hpf_BL[8];
				biquad_hpf_BL[8] = (inputSampleL * biquad_hpf_BL[4]) - (temp * biquad_hpf_BL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_hpf_BR[2]) + biquad_hpf_BR[7];
				biquad_hpf_BR[7] = (inputSampleR * biquad_hpf_BR[3]) - (temp * biquad_hpf_BR[5]) + biquad_hpf_BR[8];
				biquad_hpf_BR[8] = (inputSampleR * biquad_hpf_BR[4]) - (temp * biquad_hpf_BR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}

			if(enableLPF) {
				temp = (inputSampleL * biquad_lpf_BL[2]) + biquad_lpf_BL[7];
				biquad_lpf_BL[7] = (inputSampleL * biquad_lpf_BL[3]) - (temp * biquad_lpf_BL[5]) + biquad_lpf_BL[8];
				biquad_lpf_BL[8] = (inputSampleL * biquad_lpf_BL[4]) - (temp * biquad_lpf_BL[6]);
				inputSampleL = temp;

				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;

				temp = (inputSampleR * biquad_lpf_BR[2]) + biquad_lpf_BR[7];
				biquad_lpf_BR[7] = (inputSampleR * biquad_lpf_BR[3]) - (temp * biquad_lpf_BR[5]) + biquad_lpf_BR[8];
				biquad_lpf_BR[8] = (inputSampleR * biquad_lpf_BR[4]) - (temp * biquad_lpf_BR[6]);
				inputSampleR = temp;

				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			}
		}

		//inputSampleL = asin(inputSampleL);
		inputSampleL += (pow(inputSampleL,3)/4.0)+(pow(inputSampleL,5)/8.0)+(pow(inputSampleL,7)/16.0)+(pow(inputSampleL,9)/32.0);
		//inputSampleR = asin(inputSampleR);
		inputSampleR += (pow(inputSampleR,3)/4.0)+(pow(inputSampleR,5)/8.0)+(pow(inputSampleR,7)/16.0)+(pow(inputSampleR,9)/32.0);

		// ButterComp2
		
		if(enableComp) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;
		
			inputSampleL *= bc_inputgain;
			inputSampleR *= bc_inputgain;
		
			double divisor = bc_compfactor / (1.0+fabs(lastOutputL));
			//this is slowing compressor recovery while output waveforms were high
			divisor /= overallscale;
			double remainder = divisor;
			divisor = 1.0 - divisor;
			//recalculate divisor every sample		
		
			double inputposL = inputSampleL + 1.0;
			if (inputposL < 0.0) inputposL = 0.0;
			double outputposL = inputposL / 2.0;
			if (outputposL > 1.0) outputposL = 1.0;		
			inputposL *= inputposL;
			targetposL *= divisor;
			targetposL += (inputposL * remainder);
			double calcposL = pow((1.0/targetposL),2);
		
			double inputnegL = (-inputSampleL) + 1.0;
			if (inputnegL < 0.0) inputnegL = 0.0;
			double outputnegL = inputnegL / 2.0;
			if (outputnegL > 1.0) outputnegL = 1.0;		
			inputnegL *= inputnegL;
			targetnegL *= divisor;
			targetnegL += (inputnegL * remainder);
			double calcnegL = pow((1.0/targetnegL),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1
		
			if (inputSampleL > 0)
			{ //working on pos
				if (flip)
				{
					controlAposL *= divisor;
					controlAposL += (calcposL*remainder);
				
				}
				else
				{
					controlBposL *= divisor;
					controlBposL += (calcposL*remainder);
				}
			}
			else
			{ //working on neg
				if (flip)
				{
					controlAnegL *= divisor;
					controlAnegL += (calcnegL*remainder);
				}
				else
				{
					controlBnegL *= divisor;
					controlBnegL += (calcnegL*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'
		
			divisor = bc_compfactor / (1.0+fabs(lastOutputR));
			//this is slowing compressor recovery while output waveforms were high
			divisor /= overallscale;
			remainder = divisor;
			divisor = 1.0 - divisor;
			//recalculate divisor every sample		
		
			double inputposR = inputSampleR + 1.0;
			if (inputposR < 0.0) inputposR = 0.0;
			double outputposR = inputposR / 2.0;
			if (outputposR > 1.0) outputposR = 1.0;		
			inputposR *= inputposR;
			targetposR *= divisor;
			targetposR += (inputposR * remainder);
			double calcposR = pow((1.0/targetposR),2);
		
			double inputnegR = (-inputSampleR) + 1.0;
			if (inputnegR < 0.0) inputnegR = 0.0;
			double outputnegR = inputnegR / 2.0;
			if (outputnegR > 1.0) outputnegR = 1.0;		
			inputnegR *= inputnegR;
			targetnegR *= divisor;
			targetnegR += (inputnegR * remainder);
			double calcnegR = pow((1.0/targetnegR),2);
			//now we have mirrored targets for comp
			//outputpos and outputneg go from 0 to 1
		
			if (inputSampleR > 0)
			{ //working on pos
				if (flip)
				{
					controlAposR *= divisor;
					controlAposR += (calcposR*remainder);
				
				}
				else
				{
					controlBposR *= divisor;
					controlBposR += (calcposR*remainder);
				}
			}
			else
			{ //working on neg
				if (flip)
				{
					controlAnegR *= divisor;
					controlAnegR += (calcnegR*remainder);
				}
				else
				{
					controlBnegR *= divisor;
					controlBnegR += (calcnegR*remainder);
				}
			}
			//this causes each of the four to update only when active and in the correct 'flip'
		
			double totalmultiplierL;
			double totalmultiplierR;
			if (flip)
			{
				totalmultiplierL = (controlAposL * outputposL) + (controlAnegL * outputnegL);
				totalmultiplierR = (controlAposR * outputposR) + (controlAnegR * outputnegR);
			}
			else
			{
				totalmultiplierL = (controlBposL * outputposL) + (controlBnegL * outputnegL);
				totalmultiplierR = (controlBposR * outputposR) + (controlBnegR * outputnegR);
			}
			//this combines the sides according to flip, blending relative to the input value
		
			inputSampleL *= totalmultiplierL;
			inputSampleL /= bc_outputgain;
		
			inputSampleR *= totalmultiplierR;
			inputSampleR /= bc_outputgain;
		
			if (bc_output != 1.0) {
				inputSampleL *= bc_output;
				inputSampleR *= bc_output;
			}

			if (bc_wet !=1.0) {
				inputSampleL = (inputSampleL * bc_wet) + (drySampleL * (1.0-bc_wet));
				inputSampleR = (inputSampleR * bc_wet) + (drySampleR * (1.0-bc_wet));
			}
		
			lastOutputL = inputSampleL;
			lastOutputR = inputSampleR;
			//we will make this factor respond to use of dry/wet
		}

		// Focus

		drySampleL = inputSampleL;
		drySampleR = inputSampleR;
		
		inputSampleL = sin(inputSampleL);
		inputSampleR = sin(inputSampleR);
		//encode Console5: good cleanness
		
		temp = (inputSampleL * figureL[2]) + figureL[7];
		figureL[7] = -(temp * figureL[5]) + figureL[8];
		figureL[8] = (inputSampleL * figureL[4]) - (temp * figureL[6]);
		inputSampleL = temp;
		
		temp = (inputSampleR * figureR[2]) + figureR[7];
		figureR[7] = -(temp * figureR[5]) + figureR[8];
		figureR[8] = (inputSampleR * figureR[4]) - (temp * figureR[6]);
		inputSampleR = temp;
				
		if (inputSampleL > 1.0) inputSampleL = 1.0;
		if (inputSampleL < -1.0) inputSampleL = -1.0;
		if (inputSampleR > 1.0) inputSampleR = 1.0;
		if (inputSampleR < -1.0) inputSampleR = -1.0;
		//without this, you can get a NaN condition where it spits out DC offset at full blast!
		inputSampleL = asin(inputSampleL);
		inputSampleR = asin(inputSampleR);
		//decode Console5
		
		double groundSampleL = drySampleL - inputSampleL; //set up UnBox
		double groundSampleR = drySampleR - inputSampleR; //set up UnBox

		// uLawEncode

		if(enableuLaw) {
			if (inputSampleL > 1.0) inputSampleL = 1.0;
			if (inputSampleL < -1.0) inputSampleL = -1.0;

			if (inputSampleR > 1.0) inputSampleR = 1.0;
			if (inputSampleR < -1.0) inputSampleR = -1.0;

			if (inputSampleL > 0) inputSampleL = log(1.0+(255*fabs(inputSampleL))) / log(256);
			if (inputSampleL < 0) inputSampleL = -log(1.0+(255*fabs(inputSampleL))) / log(256);

			if (inputSampleR > 0) inputSampleR = log(1.0+(255*fabs(inputSampleR))) / log(256);
			if (inputSampleR < 0) inputSampleR = -log(1.0+(255*fabs(inputSampleR))) / log(256);
		}

		// StereoChorus

		if(enableChorus) {
			chorus_cycle++;
			if (chorus_cycle == cycleEnd) { //hit the end point and we do a chorus sample
				//assign working variables
				airFactorL = airPrevL - inputSampleL;
				if (chorus_flip) {airEvenL += airFactorL; airOddL -= airFactorL; airFactorL = airEvenL;}
				else {airOddL += airFactorL; airEvenL -= airFactorL; airFactorL = airOddL;}
				airOddL = (airOddL - ((airOddL - airEvenL)/256.0)) / 1.0001;
				airEvenL = (airEvenL - ((airEvenL - airOddL)/256.0)) / 1.0001;
				airPrevL = inputSampleL;
				inputSampleL += airFactorL;
				//left
				airFactorR = airPrevR - inputSampleR;
				if (chorus_flip) {airEvenR += airFactorR; airOddR -= airFactorR; airFactorR = airEvenR;}
				else {airOddR += airFactorR; airEvenR -= airFactorR; airFactorR = airOddR;}
				airOddR = (airOddR - ((airOddR - airEvenR)/256.0)) / 1.0001;
				airEvenR = (airEvenR - ((airEvenR - airOddR)/256.0)) / 1.0001;
				airPrevR = inputSampleR;
				inputSampleR += airFactorR;
				//right
				chorus_flip = !chorus_flip;
				//air, compensates for loss of highs in flanger's interpolation
			
				int tempL = 0;
				int tempR = 0;
				if (gcount < 1 || gcount > 32760) {gcount = 32760;}
				int count = gcount;
				pL[count+32760] = pL[count] = (int)(inputSampleL*8388352.0);
				//double buffer -8388352 to 8388352 is equal to 24 bit linear space
				double offset = depth + (depth * sin(sweepL));
				count += (int)floor(offset);
				tempL += (int)(pL[count] * (1-(offset-floor(offset)))); //less as value moves away from .0
				tempL += pL[count+1]; //we can assume always using this in one way or another?
				tempL += (int)(pL[count+2] * (offset-floor(offset))); //greater as value moves away from .0
				tempL -= (int)(((pL[count]-pL[count+1])-(pL[count+1]-pL[count+2]))/50); //interpolation hacks 'r us
				//left
			
				count = gcount;
				pR[count+32760] = pR[count] = (int)(inputSampleR*8388352.0);
				//double buffer -8388352 to 8388352 is equal to 24 bit linear space
				offset = depth + (depth * sin(sweepR));
				count += (int)floor(offset);
				tempR += (int)(pR[count] * (1-(offset-floor(offset)))); //less as value moves away from .0
				tempR += pR[count+1]; //we can assume always using this in one way or another?
				tempR += (int)(pR[count+2] * (offset-floor(offset))); //greater as value moves away from .0
				tempR -= (int)(((pR[count]-pR[count+1])-(pR[count+1]-pR[count+2]))/50); //interpolation hacks 'r us
				//right
			
				sweepL += speed;
				sweepR += speed;
				if (sweepL > tupi){sweepL -= tupi;}
				if (sweepR > tupi){sweepR -= tupi;}
				gcount--;
				//still scrolling through the samples, remember
			
				inputSampleL = ((double)(tempL/16776704.0));
				inputSampleR = ((double)(tempR/16776704.0));
				if (cycleEnd == 4) {
					lastRefL[0] = lastRefL[4]; //start from previous last
					lastRefL[2] = (lastRefL[0] + inputSampleL)/2; //half
					lastRefL[1] = (lastRefL[0] + lastRefL[2])/2; //one quarter
					lastRefL[3] = (lastRefL[2] + inputSampleL)/2; //three quarters
					lastRefL[4] = inputSampleL; //full
					lastRefR[0] = lastRefR[4]; //start from previous last
					lastRefR[2] = (lastRefR[0] + inputSampleR)/2; //half
					lastRefR[1] = (lastRefR[0] + lastRefR[2])/2; //one quarter
					lastRefR[3] = (lastRefR[2] + inputSampleR)/2; //three quarters
					lastRefR[4] = inputSampleR; //full
				}
				if (cycleEnd == 3) {
					lastRefL[0] = lastRefL[3]; //start from previous last
					lastRefL[2] = (lastRefL[0]+lastRefL[0]+inputSampleL)/3; //third
					lastRefL[1] = (lastRefL[0]+inputSampleL+inputSampleL)/3; //two thirds
					lastRefL[3] = inputSampleL; //full
					lastRefR[0] = lastRefR[3]; //start from previous last
					lastRefR[2] = (lastRefR[0]+lastRefR[0]+inputSampleR)/3; //third
					lastRefR[1] = (lastRefR[0]+inputSampleR+inputSampleR)/3; //two thirds
					lastRefR[3] = inputSampleR; //full
				}
				if (cycleEnd == 2) {
					lastRefL[0] = lastRefL[2]; //start from previous last
					lastRefL[1] = (lastRefL[0] + inputSampleL)/2; //half
					lastRefL[2] = inputSampleL; //full
					lastRefR[0] = lastRefR[2]; //start from previous last
					lastRefR[1] = (lastRefR[0] + inputSampleR)/2; //half
					lastRefR[2] = inputSampleR; //full
				}
				if (cycleEnd == 1) {
					lastRefL[0] = inputSampleL;
					lastRefR[0] = inputSampleR;
				}
				chorus_cycle = 0; //reset
				inputSampleL = lastRefL[chorus_cycle];
				inputSampleR = lastRefR[chorus_cycle];
			} else {
				inputSampleL = lastRefL[chorus_cycle];
				inputSampleR = lastRefR[chorus_cycle];
				//we are going through our references now
			}
		}

		// ClearCoat
		
		if(enableReverb) {
			double cc_drySampleL = inputSampleL;
			double cc_drySampleR = inputSampleR;
		
			cc_cycle++;
			if (cc_cycle == cycleEnd) { //hit the end point and we do a reverb sample
				aAL[countAL] = inputSampleL + (feedbackAL * 0.04166666666);
				aBL[countBL] = inputSampleL + (feedbackBL * 0.04166666666);
				aCL[countCL] = inputSampleL + (feedbackCL * 0.04166666666);
				aDL[countDL] = inputSampleL + (feedbackDL * 0.04166666666);
			
				aDR[countDR] = inputSampleR + (feedbackDR * 0.04166666666);
				aHR[countHR] = inputSampleR + (feedbackHR * 0.04166666666);
				aLR[countLR] = inputSampleR + (feedbackLR * 0.04166666666);
				aPR[countPR] = inputSampleR + (feedbackPR * 0.04166666666);
				//exactly halfway between infinite sustain at 0.0625
				//and 6dB down, almost no regen at 0.03125
				//means roughly half the results work as BitShiftGain
			
				countAL++; if (countAL < 0 || countAL > shortA) countAL = 0;
				countBL++; if (countBL < 0 || countBL > shortB) countBL = 0;
				countCL++; if (countCL < 0 || countCL > shortC) countCL = 0;
				countDL++; if (countDL < 0 || countDL > shortD) countDL = 0;
			
				countDR++; if (countDR < 0 || countDR > shortD) countDR = 0;
				countHR++; if (countHR < 0 || countHR > shortH) countHR = 0;
				countLR++; if (countLR < 0 || countLR > shortL) countLR = 0;
				countPR++; if (countPR < 0 || countPR > shortP) countPR = 0;
			
				double outAL = aAL[countAL-((countAL > shortA)?shortA+1:0)];
				double outBL = aBL[countBL-((countBL > shortB)?shortB+1:0)];
				double outCL = aCL[countCL-((countCL > shortC)?shortC+1:0)];
				double outDL = aDL[countDL-((countDL > shortD)?shortD+1:0)];
			
				double outDR = aDR[countDR-((countDR > shortD)?shortD+1:0)];
				double outHR = aHR[countHR-((countHR > shortH)?shortH+1:0)];
				double outLR = aLR[countLR-((countLR > shortL)?shortL+1:0)];
				double outPR = aPR[countPR-((countPR > shortP)?shortP+1:0)];
			
				aEL[countEL] = outAL - (outBL + outCL + outDL);
				aFL[countFL] = outBL - (outAL + outCL + outDL);
				aGL[countGL] = outCL - (outAL + outBL + outDL);
				aHL[countHL] = outDL - (outAL + outBL + outCL);
			
				aCR[countCR] = outDR - (outHR + outLR + outPR);
				aGR[countGR] = outHR - (outDR + outLR + outPR);
				aKR[countKR] = outLR - (outDR + outHR + outPR);
				aOR[countOR] = outPR - (outDR + outHR + outLR);
			
				countEL++; if (countEL < 0 || countEL > shortE) countEL = 0;
				countFL++; if (countFL < 0 || countFL > shortF) countFL = 0;
				countGL++; if (countGL < 0 || countGL > shortG) countGL = 0;
				countHL++; if (countHL < 0 || countHL > shortH) countHL = 0;
			
				countCR++; if (countCR < 0 || countCR > shortC) countCR = 0;
				countGR++; if (countGR < 0 || countGR > shortG) countGR = 0;
				countKR++; if (countKR < 0 || countKR > shortK) countKR = 0;
				countOR++; if (countOR < 0 || countOR > shortO) countOR = 0;
			
				double outEL = aEL[countEL-((countEL > shortE)?shortE+1:0)];
				double outFL = aFL[countFL-((countFL > shortF)?shortF+1:0)];
				double outGL = aGL[countGL-((countGL > shortG)?shortG+1:0)];
				double outHL = aHL[countHL-((countHL > shortH)?shortH+1:0)];
			
				double outCR = aCR[countCR-((countCR > shortC)?shortC+1:0)];
				double outGR = aGR[countGR-((countGR > shortG)?shortG+1:0)];
				double outKR = aKR[countKR-((countKR > shortK)?shortK+1:0)];
				double outOR = aOR[countOR-((countOR > shortO)?shortO+1:0)];
			
				aIL[countIL] = outEL - (outFL + outGL + outHL);
				aJL[countJL] = outFL - (outEL + outGL + outHL);
				aKL[countKL] = outGL - (outEL + outFL + outHL);
				aLL[countLL] = outHL - (outEL + outFL + outGL);
			
				aBR[countBR] = outCR - (outGR + outKR + outOR);
				aFR[countFR] = outGR - (outCR + outKR + outOR);
				aJR[countJR] = outKR - (outCR + outGR + outOR);
				aNR[countNR] = outOR - (outCR + outGR + outKR);
			
				countIL++; if (countIL < 0 || countIL > shortI) countIL = 0;
				countJL++; if (countJL < 0 || countJL > shortJ) countJL = 0;
				countKL++; if (countKL < 0 || countKL > shortK) countKL = 0;
				countLL++; if (countLL < 0 || countLL > shortL) countLL = 0;
			
				countBR++; if (countBR < 0 || countBR > shortB) countBR = 0;
				countFR++; if (countFR < 0 || countFR > shortF) countFR = 0;
				countJR++; if (countJR < 0 || countJR > shortJ) countJR = 0;
				countNR++; if (countNR < 0 || countNR > shortN) countNR = 0;
			
				double outIL = aIL[countIL-((countIL > shortI)?shortI+1:0)];
				double outJL = aJL[countJL-((countJL > shortJ)?shortJ+1:0)];
				double outKL = aKL[countKL-((countKL > shortK)?shortK+1:0)];
				double outLL = aLL[countLL-((countLL > shortL)?shortL+1:0)];
			
				double outBR = aBR[countBR-((countBR > shortB)?shortB+1:0)];
				double outFR = aFR[countFR-((countFR > shortF)?shortF+1:0)];
				double outJR = aJR[countJR-((countJR > shortJ)?shortJ+1:0)];
				double outNR = aNR[countNR-((countNR > shortN)?shortN+1:0)];
			
				aML[countML] = outIL - (outJL + outKL + outLL);
				aNL[countNL] = outJL - (outIL + outKL + outLL);
				aOL[countOL] = outKL - (outIL + outJL + outLL);
				aPL[countPL] = outLL - (outIL + outJL + outKL);
			
				aAR[countAR] = outBR - (outFR + outJR + outNR);
				aER[countER] = outFR - (outBR + outJR + outNR);
				aIR[countIR] = outJR - (outBR + outFR + outNR);
				aMR[countMR] = outNR - (outBR + outFR + outJR);
			
				countML++; if (countML < 0 || countML > shortM) countML = 0;
				countNL++; if (countNL < 0 || countNL > shortN) countNL = 0;
				countOL++; if (countOL < 0 || countOL > shortO) countOL = 0;
				countPL++; if (countPL < 0 || countPL > shortP) countPL = 0;
			
				countAR++; if (countAR < 0 || countAR > shortA) countAR = 0;
				countER++; if (countER < 0 || countER > shortE) countER = 0;
				countIR++; if (countIR < 0 || countIR > shortI) countIR = 0;
				countMR++; if (countMR < 0 || countMR > shortM) countMR = 0;
			
				double outML = aML[countML-((countML > shortM)?shortM+1:0)];
				double outNL = aNL[countNL-((countNL > shortN)?shortN+1:0)];
				double outOL = aOL[countOL-((countOL > shortO)?shortO+1:0)];
				double outPL = aPL[countPL-((countPL > shortP)?shortP+1:0)];
			
				double outAR = aAR[countAR-((countAR > shortA)?shortA+1:0)];
				double outER = aER[countER-((countER > shortE)?shortE+1:0)];
				double outIR = aIR[countIR-((countIR > shortI)?shortI+1:0)];
				double outMR = aMR[countMR-((countMR > shortM)?shortM+1:0)];			
				double outSample = (outML + outML + outML + prevMulchAL)*0.25;
				prevMulchAL = outML; outML = outSample;
				outSample = (outAR + outAR + outAR + prevMulchAR)*0.25;
				prevMulchAR = outAR; outAR = outSample;
			
				feedbackAL = outML - (outNL + outOL + outPL);
				feedbackDR = outAR - (outER + outIR + outMR);
				feedbackBL = outNL - (outML + outOL + outPL);
				feedbackHR = outER - (outAR + outIR + outMR);
				feedbackCL = outOL - (outML + outNL + outPL);
				feedbackLR = outIR - (outAR + outER + outMR);
				feedbackDL = outPL - (outML + outNL + outOL);
				feedbackPR = outMR - (outAR + outER + outIR);
				//which we need to feed back into the input again, a bit
			
				inputSampleL = (outML + outNL + outOL + outPL)/8.0;
				inputSampleR = (outAR + outER + outIR + outMR)/8.0;
				//and take the final combined sum of outputs, corrected for Householder gain
			
				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
			
				if (cycleEnd == 4) {
					cc_LastRefL[0] = cc_LastRefL[4]; //start from previous last
					cc_LastRefL[2] = (cc_LastRefL[0] + inputSampleL)/2; //half
					cc_LastRefL[1] = (cc_LastRefL[0] + cc_LastRefL[2])/2; //one quarter
					cc_LastRefL[3] = (cc_LastRefL[2] + inputSampleL)/2; //three quarters
					cc_LastRefL[4] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[4]; //start from previous last
					cc_LastRefR[2] = (cc_LastRefR[0] + inputSampleR)/2; //half
					cc_LastRefR[1] = (cc_LastRefR[0] + cc_LastRefR[2])/2; //one quarter
					cc_LastRefR[3] = (cc_LastRefR[2] + inputSampleR)/2; //three quarters
					cc_LastRefR[4] = inputSampleR; //full
				}
				if (cycleEnd == 3) {
					cc_LastRefL[0] = cc_LastRefL[3]; //start from previous last
					cc_LastRefL[2] = (cc_LastRefL[0]+cc_LastRefL[0]+inputSampleL)/3; //third
					cc_LastRefL[1] = (cc_LastRefL[0]+inputSampleL+inputSampleL)/3; //two thirds
					cc_LastRefL[3] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[3]; //start from previous last
					cc_LastRefR[2] = (cc_LastRefR[0]+cc_LastRefR[0]+inputSampleR)/3; //third
					cc_LastRefR[1] = (cc_LastRefR[0]+inputSampleR+inputSampleR)/3; //two thirds
					cc_LastRefR[3] = inputSampleR; //full
				}
				if (cycleEnd == 2) {
					cc_LastRefL[0] = cc_LastRefL[2]; //start from previous last
					cc_LastRefL[1] = (cc_LastRefL[0] + inputSampleL)/2; //half
					cc_LastRefL[2] = inputSampleL; //full
					cc_LastRefR[0] = cc_LastRefR[2]; //start from previous last
					cc_LastRefR[1] = (cc_LastRefR[0] + inputSampleR)/2; //half
					cc_LastRefR[2] = inputSampleR; //full
				}
				if (cycleEnd == 1) {
					cc_LastRefL[0] = inputSampleL;
					cc_LastRefR[0] = inputSampleR;
				}
				cc_cycle = 0; //reset
				inputSampleL = cc_LastRefL[cc_cycle];
				inputSampleR = cc_LastRefR[cc_cycle];
			} else {
				inputSampleL = cc_LastRefL[cc_cycle];
				inputSampleR = cc_LastRefR[cc_cycle];
				//we are going through our references now
			}
		
			//begin SubTight section
			double subSampleL = inputSampleL * cc_subRate;
			double subSampleR = inputSampleR * cc_subRate;
		
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
			scale = 0.5+fabs(subSampleL*0.5);
			subSampleL = (subDL+(sin(subDL-subSampleL)*scale));
			subDL = subSampleL*scale;
			scale = 0.5+fabs(subSampleR*0.5);
			subSampleR = (subDR+(sin(subDR-subSampleR)*scale));
			subDR = subSampleR*scale;
			if (subSampleL > 0.25) subSampleL = 0.25;
			if (subSampleL < -0.25) subSampleL = -0.25;
			if (subSampleR > 0.25) subSampleR = 0.25;
			if (subSampleR < -0.25) subSampleR = -0.25;
			inputSampleL -= (subSampleL*16.0);
			inputSampleR -= (subSampleR*16.0);
			//end SubTight section		
		
			if (cycleEnd > 1) {
				temp = (inputSampleL + tailL)*0.5;
				tailL = inputSampleL; inputSampleL = temp;
				temp = (inputSampleR + tailR)*0.5;
				tailR = inputSampleR; inputSampleR = temp;
			} //let's average only at elevated sample rates
		
			if (cc_wet < 1.0) {inputSampleL *= cc_wet; inputSampleR *= cc_wet;}
			if (cc_dry < 1.0) {cc_drySampleL *= cc_dry; cc_drySampleR *= cc_dry;}
			inputSampleL += cc_drySampleL; inputSampleR += cc_drySampleR;
			//this is our submix verb dry/wet: 0.5 is BOTH at FULL VOLUME
			//purpose is that, if you're adding verb, you're not altering other balances
		}

		inputSampleL *= boost; //now, focussed area gets cranked before distort
		inputSampleR *= boost; //now, focussed area gets cranked before distort
		
		switch (mode)
		{
			case 0: //Density
				if (inputSampleL > 1.570796326794897) inputSampleL = 1.570796326794897;
				if (inputSampleL < -1.570796326794897) inputSampleL = -1.570796326794897;
				if (inputSampleR > 1.570796326794897) inputSampleR = 1.570796326794897;
				if (inputSampleR < -1.570796326794897) inputSampleR = -1.570796326794897;
				//clip to 1.570796326794897 to reach maximum output
				inputSampleL = sin(inputSampleL);
				inputSampleR = sin(inputSampleR);
				break;
			case 1: //Drive				
				if (inputSampleL > 1.0) inputSampleL = 1.0;
				if (inputSampleL < -1.0) inputSampleL = -1.0;
				if (inputSampleR > 1.0) inputSampleR = 1.0;
				if (inputSampleR < -1.0) inputSampleR = -1.0;
				inputSampleL -= (inputSampleL * (fabs(inputSampleL) * 0.6) * (fabs(inputSampleL) * 0.6));
				inputSampleR -= (inputSampleR * (fabs(inputSampleR) * 0.6) * (fabs(inputSampleR) * 0.6));
				inputSampleL *= 1.6;
				inputSampleR *= 1.6;
				break;
			case 2: //Spiral
				if (inputSampleL > 1.2533141373155) inputSampleL = 1.2533141373155;
				if (inputSampleL < -1.2533141373155) inputSampleL = -1.2533141373155;
				if (inputSampleR > 1.2533141373155) inputSampleR = 1.2533141373155;
				if (inputSampleR < -1.2533141373155) inputSampleR = -1.2533141373155;
				//clip to 1.2533141373155 to reach maximum output
				inputSampleL = sin(inputSampleL * fabs(inputSampleL)) / ((fabs(inputSampleL) == 0.0) ?1:fabs(inputSampleL));
				inputSampleR = sin(inputSampleR * fabs(inputSampleR)) / ((fabs(inputSampleR) == 0.0) ?1:fabs(inputSampleR));
				break;
			case 3: //Mojo
				temp = pow(fabs(inputSampleL),0.25);
				if (temp > 0.0) inputSampleL = (sin(inputSampleL * temp * M_PI * 0.5) / temp) * 0.987654321;
				temp = pow(fabs(inputSampleR),0.25);
				if (temp > 0.0) inputSampleR = (sin(inputSampleR * temp * M_PI * 0.5) / temp) * 0.987654321;
				//mojo is the one that flattens WAAAAY out very softly before wavefolding				
				break;
			case 4: //Dyno
				temp = pow(fabs(inputSampleL),4);
				if (temp > 0.0) inputSampleL = (sin(inputSampleL * temp) / temp) * 1.1654321;
				temp = pow(fabs(inputSampleR),4);
				if (temp > 0.0) inputSampleR = (sin(inputSampleR * temp) / temp) * 1.1654321;
				//dyno is the one that tries to raise peak energy				
				break;
		}				
		
		if (output != 1.0) {
			inputSampleL *= output;
			inputSampleR *= output;
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

		inputSampleL += groundSampleL; //effectively UnBox
		inputSampleR += groundSampleR; //effectively UnBox
		
		if (wet !=1.0) {
			inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0-wet));
			inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0-wet));
		}


		// Creature

		if(enableCreature) {
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

			inputSampleL *= cr_wet;
			inputSampleR *= cr_wet;
			drySampleL *= cr_dry;
			drySampleR *= cr_dry;
			inputSampleL += drySampleL;
			inputSampleR += drySampleR;
		}

		// ToneSlant
		if(enableToneSlant) {
			for (int count = ts_overallscale; count >= 0; count--) {
				ts_bL[count+1] = ts_bL[count];
				ts_bR[count+1] = ts_bR[count];
			}

			ts_bL[0] = ts_accumulatorSampleL = inputSampleL;
			ts_bR[0] = ts_accumulatorSampleR = inputSampleR;

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
		}

		// Gain trim

		inputSampleL = inputSampleL * gain_trim;
		inputSampleR = inputSampleR * gain_trim;

		// ADClip

		double softness = 0.618033988749894848204586;

		if (adcLastSampleL >= 0.5)
		{
			if (inputSampleL < 0.5) adcLastSampleL = ((0.5*softness) + (inputSampleL * (1.0-softness)));
			else adcLastSampleL = 0.5;
		}

		if (adcLastSampleL <= -0.5)
		{
			if (inputSampleL > -0.5) adcLastSampleL = ((-0.5*softness) + (inputSampleL * (1.0-softness)));
			else adcLastSampleL = -0.5;
		}

		if (inputSampleL > 0.5)
		{
			if (adcLastSampleL < 0.5) inputSampleL = ((0.5*softness) + (adcLastSampleL * (1.0-softness)));
			else inputSampleL = 0.5;
		}

		if (inputSampleL < -0.5)
		{
			if (adcLastSampleL > -0.5) inputSampleL = ((-0.5*softness) + (adcLastSampleL * (1.0-softness)));
			else inputSampleL = -0.5;
		}
		adcLastSampleL = inputSampleL; //end ADClip L


		if (adcLastSampleR >= 0.5)
		{
			if (inputSampleR < 0.5) adcLastSampleR = ((0.5*softness) + (inputSampleR * (1.0-softness)));
			else adcLastSampleR = 0.5;
		}

		if (adcLastSampleR <= -0.5)
		{
			if (inputSampleR > -0.5) adcLastSampleR = ((-0.5*softness) + (inputSampleR * (1.0-softness)));
			else adcLastSampleR = -0.5;
		}

		if (inputSampleR > 0.5)
		{
			if (adcLastSampleR < 0.5) inputSampleR = ((0.5*softness) + (adcLastSampleR * (1.0-softness)));
			else inputSampleR = 0.5;
		}

		if (inputSampleR < -0.5)
		{
			if (adcLastSampleR > -0.5) inputSampleR = ((-0.5*softness) + (adcLastSampleR * (1.0-softness)));
			else inputSampleR = -0.5;
		}
		adcLastSampleR = inputSampleR; //end ADClip R

		if (inputSampleL > 0.5) inputSampleL = 0.5;
		if (inputSampleL < -0.5) inputSampleL = -0.5;
		//final iron bar
		if (inputSampleR > 0.5) inputSampleR = 0.5;
		if (inputSampleR < -0.5) inputSampleR = -0.5;
		//final iron bar

		// Channel9

		//begin L
		temp = (ch9_lastSampleBL - ch9_lastSampleCL) * 0.381966011250105;
		temp -= (ch9_lastSampleAL - ch9_lastSampleBL) * 0.6180339887498948482045;
		temp += inputSampleL - ch9_lastSampleAL; //regular slew clamping added
		
		ch9_lastSampleCL = ch9_lastSampleBL;
		ch9_lastSampleBL = ch9_lastSampleAL;
		ch9_lastSampleAL = inputSampleL; //now our output relates off ch9_lastSampleB
		
		if (temp > localthreshold)
			inputSampleL = ch9_lastSampleBL + localthreshold;
		if (-temp > localthreshold)
			inputSampleL = ch9_lastSampleBL - localthreshold;
		
		ch9_lastSampleAL = (ch9_lastSampleAL*0.381966011250105)+(inputSampleL*0.6180339887498948482045); //split the difference between raw and smoothed for buffer
		//end L
		
		//begin R
		temp = (ch9_lastSampleBR - ch9_lastSampleCR) * 0.381966011250105;
		temp -= (ch9_lastSampleAR - ch9_lastSampleBR) * 0.6180339887498948482045;
		temp += inputSampleR - ch9_lastSampleAR; //regular slew clamping added
		
		ch9_lastSampleCR = ch9_lastSampleBR;
		ch9_lastSampleBR = ch9_lastSampleAR;
		ch9_lastSampleAR = inputSampleR; //now our output relates off ch9_lastSampleB
		
		if (temp > localthreshold)
			inputSampleR = ch9_lastSampleBR + localthreshold;
		if (-temp > localthreshold)
			inputSampleR = ch9_lastSampleBR - localthreshold;
		
		ch9_lastSampleAR = (ch9_lastSampleAR*0.381966011250105)+(inputSampleR*0.6180339887498948482045); //split the difference between raw and smoothed for buffer
		//end R
		
		temp = ch9_biquadB[2]*inputSampleL+ch9_biquadB[3]*ch9_biquadB[7]+ch9_biquadB[4]*ch9_biquadB[8]-ch9_biquadB[5]*ch9_biquadB[9]-ch9_biquadB[6]*ch9_biquadB[10];
		ch9_biquadB[8] = ch9_biquadB[7]; ch9_biquadB[7] = inputSampleL; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleL = temp;
		ch9_biquadB[10] = ch9_biquadB[9]; ch9_biquadB[9] = inputSampleL; //DF1 left
		temp = ch9_biquadB[2]*inputSampleR+ch9_biquadB[3]*ch9_biquadB[11]+ch9_biquadB[4]*ch9_biquadB[12]-ch9_biquadB[5]*ch9_biquadB[13]-ch9_biquadB[6]*ch9_biquadB[14];
		ch9_biquadB[12] = ch9_biquadB[11]; ch9_biquadB[11] = inputSampleR; if (fabs(temp)<1.18e-37) temp = 0.0; inputSampleR = temp;
		ch9_biquadB[14] = ch9_biquadB[13]; ch9_biquadB[13] = inputSampleR; //DF1 right

		if(flip) {
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
			--sampleFrames;
    }

		flip = !flip;
  }
}
