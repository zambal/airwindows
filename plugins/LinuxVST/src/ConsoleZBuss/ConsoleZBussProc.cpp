/* ========================================
 *  ConsoleZBuss - ConsoleZBuss.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZBuss_H
#include "ConsoleZBuss.h"
#endif

void ConsoleZBuss::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames)
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	VstInt32 inFramesToProcess = sampleFrames; //vst doesn't give us this as a separate variable so we'll make it
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	// ConsoleLABuss

	gainA = gainB;
	gainB = sqrt(A); //smoothed master fader from Z2 filters
	//this will be applied three times: this is to make the various tone alterations
	//hit differently at different master fader drive levels.
	//in particular, backing off the master fader tightens the super lows
	//but opens up the modified Sinew, because more of the attentuation happens before
	//you even get to slew clipping :) and if the fader is not active, it bypasses completely.

	double threshSinew = 0.718/overallscale;
	double subTrim = 0.0011 / overallscale;

	// Desk4

	double gain = (pow(B,2)*10)+0.0001;
	double gaintrim = (pow(B,2)*2)+1.0;
	double slewgain = (pow(C,3)*40)+0.0001;
	double prevslew = 0.105;
	double intensity = (pow(D,6)*15)+0.0001;
	double depthA = (pow(E,4)*940)+0.00001;
	int offsetA = (int)(depthA * overallscale);
	if (offsetA < 1) offsetA = 1;
	if (offsetA > 4880) offsetA = 4880;
	double balanceB = 0.0001;
	slewgain *= overallscale;
	prevslew *= overallscale;
	balanceB /= overallscale;
	double desk4_outputgain = F;
	double wet = G;
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
	double drySampleL;
	double drySampleR;
	double inputSampleL;
	double inputSampleR;

	// Buttercomp2

	double inputgain = pow(10.0,(H*14.0)/20.0);
	double compfactor = 0.012 * (H / 135.0);
	double output = I * 2.0;
	double butter_wet = J;
	//removed extra dry variable
	double outputgain = inputgain;
	outputgain -= 1.0;
	outputgain /= 1.5;
	outputgain += 1.0;

	// VariMu

	double threshold = 1.001 - (1.0-pow(1.0-K,3));
	double muMakeupGain = sqrt(1.0 / threshold);
	muMakeupGain = (muMakeupGain + sqrt(muMakeupGain))/2.0;
	muMakeupGain = sqrt(muMakeupGain);
	double outGain = sqrt(muMakeupGain);
	//gain settings around threshold
	double release = pow((1.15-L),5)*32768.0;
	release /= overallscale;
	double fastest = sqrt(release);
	//speed settings around release
	double coefficient;
	double mu_output = outGain * M;
	double mu_wet = N;
	double squaredSampleL;
	double squaredSampleR;

	// BussColors4

	const int maxConvolutionBufferSize = (int)(34.0 * overallscale); //we won't use more of the buffer than we have to
	for (int count = 0; count < 34; count++) c[count] = (int)(count * overallscale); //assign conv taps
	double bc_drySampleL;
	double bc_drySampleR;
	double applyconvL;
	double applyconvR;
	int bc_offsetA = 3;
	double dynamicconvL = 3.0;
	double dynamicconvR = 3.0;
	double bc_gain = 0.436;
	double bc_outgain = 1.0;
	double bc_bridgerectifier;

	double bc_inputSampleL;
	double bc_inputSampleR;

	int console = (int)( O * 7.999 )+1; //the AU used a 1-8 index, will just keep it
	switch (console)
	{
		case 1: bc_offsetA = 4; bc_gain = g[1]; bc_outgain = outg[1]; break; //Dark (Focusrite)
		case 2: bc_offsetA = 3; bc_gain = g[2]; bc_outgain = outg[2]; break; //Rock (SSL)
		case 3: bc_offsetA = 5; bc_gain = g[3]; bc_outgain = outg[3]; break; //Lush (Neve)
		case 4: bc_offsetA = 8; bc_gain = g[4]; bc_outgain = outg[4]; break; //Vibe (Elation)
		case 5: bc_offsetA = 5; bc_gain = g[5]; bc_outgain = outg[5]; break; //Holo (Precision 8)
		case 6: bc_offsetA = 7; bc_gain = g[6]; bc_outgain = outg[6]; break; //Punch (API)
		case 7: bc_offsetA = 7; bc_gain = g[7]; bc_outgain = outg[7]; break; //Steel (Calibre)
		case 8: bc_offsetA = 6; bc_gain = g[8]; bc_outgain = outg[8]; break; //Tube (Manley)
	}
	bc_offsetA = (int)(bc_offsetA * overallscale); //we extend the sag buffer too, at high sample rates
	if (bc_offsetA < 3) bc_offsetA = 3; //are we getting divides by zero?
	bc_gain *= pow(10.0,((P * 36.0)-18.0)/14.0); //add adjustment factor
	bc_outgain *= pow(10.0,(((Q * 36.0)-18.0)+3.3)/14.0); //add adjustment factor
	double bc_wet = R;
	//removed extra dry variable


    while (--sampleFrames >= 0)
    {
		double inputSampleL = *in1;
		double inputSampleR = *in2;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// Desk4

		if(wet > 0.0) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			if (desk4_gcount < 0 || desk4_gcount > 4900) {desk4_gcount = 4900;}

			desk4_dL[desk4_gcount+4900] = desk4_dL[desk4_gcount] = fabs(inputSampleL)*intensity;
			desk4_controlL += (desk4_dL[desk4_gcount] / offsetA);
			desk4_controlL -= (desk4_dL[desk4_gcount+offsetA] / offsetA);
			desk4_controlL -= 0.000001;
			clampL = 1;
			if (desk4_controlL < 0) {desk4_controlL = 0;}
			if (desk4_controlL > 1) {clampL -= (desk4_controlL - 1); desk4_controlL = 1;}
			if (clampL < 0.5) {clampL = 0.5;}

			desk4_dR[desk4_gcount+4900] = desk4_dR[desk4_gcount] = fabs(inputSampleR)*intensity;
			desk4_controlR += (desk4_dR[desk4_gcount] / offsetA);
			desk4_controlR -= (desk4_dR[desk4_gcount+offsetA] / offsetA);
			desk4_controlR -= 0.000001;
			clampR = 1;
			if (desk4_controlR < 0) {desk4_controlR = 0;}
			if (desk4_controlR > 1) {clampR -= (desk4_controlR - 1); desk4_controlR = 1;}
			if (clampR < 0.5) {clampR = 0.5;}


			desk4_gcount--;
			//control = 0 to 1
			thicknessL = ((1.0 - desk4_controlL) * 2.0) - 1.0;
			thicknessR = ((1.0 - desk4_controlR) * 2.0) - 1.0;

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

			inputSampleL *= gain;
			bridgerectifier = fabs(inputSampleL);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleL > 0) inputSampleL = bridgerectifier;
			else inputSampleL = -bridgerectifier;
			//drive section
			inputSampleL /= gain;
			inputSampleL *= gaintrim;
			//end of Desk section

			inputSampleR *= gain;
			bridgerectifier = fabs(inputSampleR);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleR > 0) inputSampleR = bridgerectifier;
			else inputSampleR = -bridgerectifier;
			//drive section
			inputSampleR /= gain;
			inputSampleR *= gaintrim;
			//end of Desk section

			if (desk4_outputgain != 1.0) {
				inputSampleL *= desk4_outputgain;
				inputSampleR *= desk4_outputgain;
			}

			if (wet !=1.0) {
				inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0-wet));
				inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0-wet));
			}
		}

		// ConsoleLABuss

		double temp = (double)sampleFrames/inFramesToProcess;
		double la_gain = (gainA*temp)+(gainB*(1.0-temp));
		//setting up smoothed master fader

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
		inputSampleL -= (subSampleL*16.0);
		inputSampleR -= (subSampleR*16.0);
		//end SubTight section

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//subtight is always fully engaged: tighten response when restraining full console

		//begin ConsoleZeroBuss which is the one we choose for ConsoleZ
		if (inputSampleL > 2.8) inputSampleL = 2.8;
		if (inputSampleL < -2.8) inputSampleL = -2.8;
		if (inputSampleL > 0.0) inputSampleL = (inputSampleL*2.0)/(3.0-inputSampleL);
		else inputSampleL = -(inputSampleL*-2.0)/(3.0+inputSampleL);

		if (inputSampleR > 2.8) inputSampleR = 2.8;
		if (inputSampleR < -2.8) inputSampleR = -2.8;
		if (inputSampleR > 0.0) inputSampleR = (inputSampleR*2.0)/(3.0-inputSampleR);
		else inputSampleR = -(inputSampleR*-2.0)/(3.0+inputSampleR);
		//ConsoleZero Buss

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//after C0Buss but before EverySlew: allow highs to come out a bit more
		//when pulling back master fader. Less drive equals more open

		temp = inputSampleL;
		double clamp = inputSampleL - lastSinewL;
		if (lastSinewL > 1.0) lastSinewL = 1.0;
		if (lastSinewL < -1.0) lastSinewL = -1.0;
		double sinew = threshSinew * cos(lastSinewL);
		if (clamp > sinew) temp = lastSinewL + sinew;
		if (-clamp > sinew) temp = lastSinewL - sinew;
		inputSampleL = lastSinewL = temp;
		temp = inputSampleR;
		clamp = inputSampleR - lastSinewR;
		if (lastSinewR > 1.0) lastSinewR = 1.0;
		if (lastSinewR < -1.0) lastSinewR = -1.0;
		sinew = threshSinew * cos(lastSinewR);
		if (clamp > sinew) temp = lastSinewR + sinew;
		if (-clamp > sinew) temp = lastSinewR - sinew;
		inputSampleR = lastSinewR = temp;

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//after EverySlew fades the total output sound: least change in tone here.

		// Buttercomp2

		if(butter_wet > 0.0) {

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
			//effectively digital black, we'll subtract it aButterComp2. We want a 'air' hiss
			double butter_drySampleL = inputSampleL;
			double butter_drySampleR = inputSampleR;

			inputSampleL *= inputgain;
			inputSampleR *= inputgain;

			double divisor = compfactor / (1.0+fabs(lastOutputL));
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
				if (butter_flip)
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
				if (butter_flip)
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
			//this causes each of the four to update only when active and in the correct 'butter_flip'

			divisor = compfactor / (1.0+fabs(lastOutputR));
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
				if (butter_flip)
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
				if (butter_flip)
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
			//this causes each of the four to update only when active and in the correct 'butter_flip'

			double totalmultiplierL;
			double totalmultiplierR;
			if (butter_flip)
			{
				totalmultiplierL = (controlAposL * outputposL) + (controlAnegL * outputnegL);
				totalmultiplierR = (controlAposR * outputposR) + (controlAnegR * outputnegR);
			}
			else
			{
				totalmultiplierL = (controlBposL * outputposL) + (controlBnegL * outputnegL);
				totalmultiplierR = (controlBposR * outputposR) + (controlBnegR * outputnegR);
			}
			//this combines the sides according to butter_flip, blending relative to the input value

			inputSampleL *= totalmultiplierL;
			inputSampleL /= outputgain;

			inputSampleR *= totalmultiplierR;
			inputSampleR /= outputgain;

			if (output != 1.0) {
				inputSampleL *= output;
				inputSampleR *= output;
			}

			if (butter_wet !=1.0) {
				inputSampleL = (inputSampleL * butter_wet) + (butter_drySampleL * (1.0-butter_wet));
				inputSampleR = (inputSampleR * butter_wet) + (butter_drySampleR * (1.0-butter_wet));
			}

			lastOutputL = inputSampleL;
			lastOutputR = inputSampleR;
			//we will make this factor respond to use of dry/wet

			butter_flip = !butter_flip;
		}

		// VariMu
		if(mu_wet > 0.0) {
			static int mu_noisesourceL = 0;
			static int mu_noisesourceR = 850010;
			int mu_residue;
			double mu_applyresidue;

			mu_noisesourceL = mu_noisesourceL % 1700021; mu_noisesourceL++;
			mu_residue = mu_noisesourceL * mu_noisesourceL;
			mu_residue = mu_residue % 170003; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17011; mu_residue *= mu_residue;
			mu_residue = mu_residue % 1709; mu_residue *= mu_residue;
			mu_residue = mu_residue % 173; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17;
			mu_applyresidue = mu_residue;
			mu_applyresidue *= 0.00000001;
			mu_applyresidue *= 0.00000001;
			inputSampleL += mu_applyresidue;
			if (inputSampleL<1.2e-38 && -inputSampleL<1.2e-38) {
				inputSampleL -= mu_applyresidue;
			}

			mu_noisesourceR = mu_noisesourceR % 1700021; mu_noisesourceR++;
			mu_residue = mu_noisesourceR * mu_noisesourceR;
			mu_residue = mu_residue % 170003; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17011; mu_residue *= mu_residue;
			mu_residue = mu_residue % 1709; mu_residue *= mu_residue;
			mu_residue = mu_residue % 173; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17;
			mu_applyresidue = mu_residue;
			mu_applyresidue *= 0.00000001;
			mu_applyresidue *= 0.00000001;
			inputSampleR += mu_applyresidue;
			if (inputSampleR<1.2e-38 && -inputSampleR<1.2e-38) {
				inputSampleR -= mu_applyresidue;
			}
			//for live air, we always apply the dither noise. Then, if our result is
			//effectively digital black, we'll subtract it aVariMu. We want a 'air' hiss
			double mu_drySampleL = inputSampleL;
			double mu_drySampleR = inputSampleR;

			if (fabs(inputSampleL) > fabs(previousL)) squaredSampleL = previousL * previousL;
			else squaredSampleL = inputSampleL * inputSampleL;
			previousL = inputSampleL;
			inputSampleL *= muMakeupGain;

			if (fabs(inputSampleR) > fabs(previousR)) squaredSampleR = previousR * previousR;
			else squaredSampleR = inputSampleR * inputSampleR;
			previousR = inputSampleR;
			inputSampleR *= muMakeupGain;

			//adjust coefficients for L
			if (flip)
			{
				if (fabs(squaredSampleL) > threshold)
				{
					muVaryL = threshold / fabs(squaredSampleL);
					muAttackL = sqrt(fabs(muSpeedAL));
					muCoefficientAL = muCoefficientAL * (muAttackL-1.0);
					if (muVaryL < threshold)
					{
						muCoefficientAL = muCoefficientAL + threshold;
					}
					else
					{
						muCoefficientAL = muCoefficientAL + muVaryL;
					}
					muCoefficientAL = muCoefficientAL / muAttackL;
				}
				else
				{
					muCoefficientAL = muCoefficientAL * ((muSpeedAL * muSpeedAL)-1.0);
					muCoefficientAL = muCoefficientAL + 1.0;
					muCoefficientAL = muCoefficientAL / (muSpeedAL * muSpeedAL);
				}
				muNewSpeedL = muSpeedAL * (muSpeedAL-1);
				muNewSpeedL = muNewSpeedL + fabs(squaredSampleL*release)+fastest;
				muSpeedAL = muNewSpeedL / muSpeedAL;
			}
			else
			{
				if (fabs(squaredSampleL) > threshold)
				{
					muVaryL = threshold / fabs(squaredSampleL);
					muAttackL = sqrt(fabs(muSpeedBL));
					muCoefficientBL = muCoefficientBL * (muAttackL-1);
					if (muVaryL < threshold)
					{
						muCoefficientBL = muCoefficientBL + threshold;
					}
					else
					{
						muCoefficientBL = muCoefficientBL + muVaryL;
					}
					muCoefficientBL = muCoefficientBL / muAttackL;
				}
				else
				{
					muCoefficientBL = muCoefficientBL * ((muSpeedBL * muSpeedBL)-1.0);
					muCoefficientBL = muCoefficientBL + 1.0;
					muCoefficientBL = muCoefficientBL / (muSpeedBL * muSpeedBL);
				}
				muNewSpeedL = muSpeedBL * (muSpeedBL-1);
				muNewSpeedL = muNewSpeedL + fabs(squaredSampleL*release)+fastest;
				muSpeedBL = muNewSpeedL / muSpeedBL;
			}
			//got coefficients, adjusted speeds for L

			//adjust coefficients for R
			if (flip)
			{
				if (fabs(squaredSampleR) > threshold)
				{
					muVaryR = threshold / fabs(squaredSampleR);
					muAttackR = sqrt(fabs(muSpeedAR));
					muCoefficientAR = muCoefficientAR * (muAttackR-1.0);
					if (muVaryR < threshold)
					{
						muCoefficientAR = muCoefficientAR + threshold;
					}
					else
					{
						muCoefficientAR = muCoefficientAR + muVaryR;
					}
					muCoefficientAR = muCoefficientAR / muAttackR;
				}
				else
				{
					muCoefficientAR = muCoefficientAR * ((muSpeedAR * muSpeedAR)-1.0);
					muCoefficientAR = muCoefficientAR + 1.0;
					muCoefficientAR = muCoefficientAR / (muSpeedAR * muSpeedAR);
				}
				muNewSpeedR = muSpeedAR * (muSpeedAR-1);
				muNewSpeedR = muNewSpeedR + fabs(squaredSampleR*release)+fastest;
				muSpeedAR = muNewSpeedR / muSpeedAR;
			}
			else
			{
				if (fabs(squaredSampleR) > threshold)
				{
					muVaryR = threshold / fabs(squaredSampleR);
					muAttackR = sqrt(fabs(muSpeedBR));
					muCoefficientBR = muCoefficientBR * (muAttackR-1);
					if (muVaryR < threshold)
					{
						muCoefficientBR = muCoefficientBR + threshold;
					}
					else
					{
						muCoefficientBR = muCoefficientBR + muVaryR;
					}
					muCoefficientBR = muCoefficientBR / muAttackR;
				}
				else
				{
					muCoefficientBR = muCoefficientBR * ((muSpeedBR * muSpeedBR)-1.0);
					muCoefficientBR = muCoefficientBR + 1.0;
					muCoefficientBR = muCoefficientBR / (muSpeedBR * muSpeedBR);
				}
				muNewSpeedR = muSpeedBR * (muSpeedBR-1);
				muNewSpeedR = muNewSpeedR + fabs(squaredSampleR*release)+fastest;
				muSpeedBR = muNewSpeedR / muSpeedBR;
			}
			//got coefficients, adjusted speeds for R

			if (flip)
			{
				coefficient = (muCoefficientAL + pow(muCoefficientAL,2))/2.0;
				inputSampleL *= coefficient;
				coefficient = (muCoefficientAR + pow(muCoefficientAR,2))/2.0;
				inputSampleR *= coefficient;
			}
			else
			{
				coefficient = (muCoefficientBL + pow(muCoefficientBL,2))/2.0;
				inputSampleL *= coefficient;
				coefficient = (muCoefficientBR + pow(muCoefficientBR,2))/2.0;
				inputSampleR *= coefficient;
			}
			//applied compression with vari-vari-µ-µ-µ-µ-µ-µ-is-the-kitten-song o/~
			//applied gain correction to control output level- tends to constrain sound rather than inflate it
			flip = !flip;

			if (mu_output < 1.0) {
				inputSampleL *= mu_output;
				inputSampleR *= mu_output;
			}
			if (mu_wet < 1.0) {
				inputSampleL = (mu_drySampleL * (1.0-mu_wet)) + (inputSampleL * mu_wet);
				inputSampleR = (mu_drySampleR * (1.0-mu_wet)) + (inputSampleR * mu_wet);
			}
			//nice little output stage template: if we have another scale of floating point
			//number, we really don't want to meaninglessly multiply that by 1.0.
		}
		// BussColors4

		if(bc_wet > 0.0) {

			bc_drySampleL = inputSampleL;
			bc_drySampleR = inputSampleR;

			if (bc_gain != 1.0) {
				inputSampleL *= bc_gain;
				inputSampleR *= bc_gain;
			}


			bc_bridgerectifier = fabs(inputSampleL);
			slowdynL *= 0.999;
			slowdynL += bc_bridgerectifier;
			if (slowdynL > 1.5) {slowdynL = 1.5;}
			//before the iron bar- fry console for crazy behavior
			dynamicconvL = 2.5 + slowdynL;

			if (bc_bridgerectifier > 1.57079633) bc_bridgerectifier = 1.0;
			else bc_bridgerectifier = sin(bc_bridgerectifier);
			if (inputSampleL > 0) inputSampleL = bc_bridgerectifier;
			else inputSampleL = -bc_bridgerectifier;
			//end pre saturation stage L

			bc_bridgerectifier = fabs(inputSampleR);
			slowdynR *= 0.999;
			slowdynR += bc_bridgerectifier;
			if (slowdynR > 1.5) {slowdynR = 1.5;}
			//before the iron bar- fry console for crazy behavior
			dynamicconvR = 2.5 + slowdynR;

			if (bc_bridgerectifier > 1.57079633) bc_bridgerectifier = 1.0;
			else bc_bridgerectifier = sin(bc_bridgerectifier);
			if (inputSampleR > 0) inputSampleR = bc_bridgerectifier;
			else inputSampleR = -bc_bridgerectifier;
			//end pre saturation stage R

			if (gcount < 0 || gcount > 44) gcount = 44;

			dL[gcount+44] = dL[gcount] = fabs(inputSampleL);
			controlL += (dL[gcount] / bc_offsetA);
			controlL -= (dL[gcount+bc_offsetA] / bc_offsetA);
			controlL -= 0.000001;
			if (controlL < 0) controlL = 0;
			if (controlL > 100) controlL = 100;
			applyconvL = (controlL / bc_offsetA) * dynamicconvL;
			//now we have a 'sag' style average to apply to the conv, L

			dR[gcount+44] = dR[gcount] = fabs(inputSampleR);
			controlR += (dR[gcount] / bc_offsetA);
			controlR -= (dR[gcount+bc_offsetA] / bc_offsetA);
			controlR -= 0.000001;
			if (controlR < 0) controlR = 0;
			if (controlR > 100) controlR = 100;
			applyconvR = (controlR / bc_offsetA) * dynamicconvR;
			//now we have a 'sag' style average to apply to the conv, R

			gcount--;

			//now the convolution
			for (int count = maxConvolutionBufferSize; count > 0; --count) {bL[count] = bL[count-1];} //was 173
			//we're only doing assigns, and it saves us an add inside the convolution calculation
			//therefore, we'll just assign everything one step along and have our buffer that way.
			bL[0] = inputSampleL;

			for (int count = maxConvolutionBufferSize; count > 0; --count) {bR[count] = bR[count-1];} //was 173
			//we're only doing assigns, and it saves us an add inside the convolution calculation
			//therefore, we'll just assign everything one step along and have our buffer that way.
			bR[0] = inputSampleR;
			//The reason these are separate is, since it's just a big assign-fest, it's possible compilers can
			//come up with a clever way to do it. Interleaving the samples might make it enough more complicated
			//that the compiler wouldn't know to do 'em in batches or whatever. Just a thought, profiling would
			//be the correct way to find out this (or indeed, whether doing another add insode the convolutions would
			//be the best bet,

			//The convolutions!

			switch (console)
			{
				case 1:
					//begin Cider (Focusrite) MCI
					inputSampleL += (bL[c[1]] * (0.61283288942201319  + (0.00024011410669522*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.24036380659761222  - (0.00020789518206241*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.09104669761717916  + (0.00012829642741548*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.02378290768554025  - (0.00017673646470440*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.02832818490275965  - (0.00013536187747384*applyconvL)));
					inputSampleL += (bL[c[6]] * (0.03268797679215937  + (0.00015035126653359*applyconvL)));
					inputSampleL -= (bL[c[7]] * (0.04024464202655586  - (0.00015034923056735*applyconvL)));
					inputSampleL += (bL[c[8]] * (0.01864890074318696  + (0.00014513281680642*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.01632731954100322  - (0.00015509089075614*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.00318907090555589  - (0.00014784812076550*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.00208573465221869  - (0.00015350520779465*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.00907033901519614  - (0.00015442964157250*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00199458794148013  - (0.00015595640046297*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.00705979153201755  - (0.00015730069418051*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00429488975412722  - (0.00015743697943505*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.00497724878704936  - (0.00016014760011861*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00506059305562353  - (0.00016194824072466*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00483432223285621  - (0.00016329050124225*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00495100420886005  - (0.00016297509798749*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00489319520555115  - (0.00016472839684661*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00489177657970308  - (0.00016791875866630*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00487900894707044  - (0.00016755993898534*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00486234009335561  - (0.00016968157345446*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00485737490288736  - (0.00017180713324431*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00484106070563455  - (0.00017251073661092*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00483219429408410  - (0.00017321683790891*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00482013597437550  - (0.00017392186866488*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00480949628051497  - (0.00017569098775602*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00479992055604049  - (0.00017746046369449*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00478750757986987  - (0.00017745630047554*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00477828651185740  - (0.00017958043287604*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00476906544384494  - (0.00018170456527653*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00475700712413634  - (0.00018099144598088*applyconvL)));
					//end Cider (Focusrite) MCI
					break;
				case 2:
					//begin Rock (SSL) conv
					inputSampleL += (bL[c[1]] * (0.67887916185274055  + (0.00068787552301086*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.25671050678827934  + (0.00017691749454490*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.15135839896615280  + (0.00007481480365043*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.11813512969090802  + (0.00005191138121359*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.08329104347166429  + (0.00001871054659794*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.07663817456103936  + (0.00002751359071705*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.05477586152148759  + (0.00000744843212679*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.05547314737187786  + (0.00001025289931145*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.03822948356540711  - (0.00000249791561457*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.04199383340841713  - (0.00000067328840674*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.02695796542339694  - (0.00000796704606548*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.03228715059431878  - (0.00000579711816722*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.01846929689819187  - (0.00000984017804950*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.02528050435045951  - (0.00000701189792484*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.01207844846859765  - (0.00001522630289356*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01894464378378515  - (0.00001205456372080*applyconvL)));
					inputSampleL += (bL[c[17]] * (0.00667804407593324  - (0.00001343604283817*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.01408418045473130  - (0.00001246443581504*applyconvL)));
					inputSampleL += (bL[c[19]] * (0.00228696509481569  - (0.00001506764046927*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.01006277891348454  - (0.00000970723079112*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00132368373546377  + (0.00001188847238761*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00676615715578373  - (0.00001209129844861*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00426288438418556  + (0.00001286836943559*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00408897698639688  - (0.00001102542567911*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00662040619382751  + (0.00001206328529063*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00196101294183599  - (0.00000950703614981*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00845620581010342  + (0.00001279970295678*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00032595215043616  - (0.00000920518241371*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00982957737435458  + (0.00001177745362317*applyconvL)));
					inputSampleL += (bL[c[30]] * (0.00086920573760513  + (0.00000913758382404*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.01079020871452061  + (0.00000900750153697*applyconvL)));
					inputSampleL += (bL[c[32]] * (0.00167613606334460  + (0.00000732769151038*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.01138050011044332  + (0.00000946908207442*applyconvL)));
					//end Rock (SSL) conv
					break;
				case 3:
					//begin Lush (Neve) conv
					inputSampleL += (bL[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvL)));
					inputSampleL += (bL[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvL)));
					inputSampleL += (bL[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvL)));
					inputSampleL -= (bL[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvL)));
					inputSampleL += (bL[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvL)));
					inputSampleL += (bL[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvL)));
					inputSampleL += (bL[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvL)));
					inputSampleL += (bL[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvL)));
					inputSampleL += (bL[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvL)));
					inputSampleL += (bL[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvL)));
					inputSampleL += (bL[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvL)));
					inputSampleL += (bL[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvL)));
					inputSampleL += (bL[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvL)));
					inputSampleL += (bL[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvL)));
					inputSampleL += (bL[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvL)));
					inputSampleL += (bL[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvL)));
					inputSampleL += (bL[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvL)));
					inputSampleL += (bL[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvL)));
					inputSampleL += (bL[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvL)));
					inputSampleL += (bL[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvL)));
					inputSampleL += (bL[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvL)));
					inputSampleL += (bL[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvL)));
					inputSampleL += (bL[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvL)));
					//end Lush (Neve) conv
					break;
				case 4:
					//begin Elation (LA2A) vibe
					inputSampleL -= (bL[c[1]] * (0.25867935358656502  - (0.00045755657070112*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.11509367290253694  - (0.00017494270657228*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.06709853575891785  - (0.00058913102597723*applyconvL)));
					inputSampleL += (bL[c[4]] * (0.01871006356851681  - (0.00003387358004645*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00794797957360465  - (0.00044224784691203*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.01956921817394220  - (0.00006718936750076*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01682120257195205  + (0.00032857446292230*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03401069039824205  - (0.00013634182872897*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.02369950268232634  + (0.00023112685751657*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.03477071178117132  - (0.00018029792231600*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.02024369717958201  + (0.00017337813374202*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.02819087729102172  - (0.00021438538665420*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.01147946743141303  + (0.00014424066034649*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.01894777011468867  - (0.00021549146262408*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.00301370330346873  + (0.00013527460148394*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01067147835815486  - (0.00020960689910868*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00402715397506384  - (0.00014421582712470*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00502221703392005  - (0.00019805767015024*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00808788533308497  - (0.00016095444141931*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00232696588842683  - (0.00018384470981829*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00943950821324531  - (0.00017098987347593*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00193709517200834  - (0.00018151995939591*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00899713952612659  - (0.00017385835059948*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00280584331659089  - (0.00017742164162470*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00780381001954970  - (0.00018002500755708*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00400370310490333  - (0.00017471691087957*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00661527728186928  - (0.00018137323370347*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00496545526864518  - (0.00017681872601767*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00580728820997532  - (0.00018186220389790*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00549309984725666  - (0.00017722985399075*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00542194777529239  - (0.00018486900185338*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00565992080998939  - (0.00018005824393118*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00532121562846656  - (0.00018643189636216*applyconvL)));
					//end Elation (LA2A)
					break;
				case 5:
					//begin Precious (Precision 8) Holo
					inputSampleL += (bL[c[1]] * (0.59188440274551890  - (0.00008361469668405*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.24439750948076133  + (0.00002651678396848*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.14109876103205621  - (0.00000840487181372*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.10053507128157971  + (0.00001768100964598*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.05859287880626238  - (0.00000361398065989*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.04337406889823660  + (0.00000735941182117*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01589900680531097  + (0.00000207347387987*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.01087234854973281  + (0.00000732123412029*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.00845782429679176  - (0.00000133058605071*applyconvL)));
					inputSampleL += (bL[c[10]] * (0.00662278586618295  - (0.00000424594730611*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.02000592193760155  + (0.00000632896879068*applyconvL)));
					inputSampleL += (bL[c[12]] * (0.01321157777167565  - (0.00001421171592570*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.02249955362988238  + (0.00000163937127317*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.01196492077581504  - (0.00000535385220676*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.01905917427000097  + (0.00000121672882030*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00761909482108073  - (0.00000326242895115*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01362744780256239  + (0.00000359274216003*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00200183122683721  - (0.00000089207452791*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00833042637239315  + (0.00000946767677294*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00258481175207224  - (0.00000087429351464*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00459744479712244  - (0.00000049519758701*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00534277030993820  + (0.00000397547847155*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00272332919605675  + (0.00000040077229097*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00637243782359372  - (0.00000139419072176*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00233001590327504  + (0.00000420129915747*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00623296727793041  + (0.00000019010664856*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00276177096376805  + (0.00000580301901385*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00559184754866264  + (0.00000080597287792*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00343180144395919  - (0.00000243701142085*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00493325428861701  + (0.00000300985740900*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00396140827680823  - (0.00000051459681789*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00448497879902493  + (0.00000744412841743*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00425146888772076  - (0.00000082346016542*applyconvL)));
					//end Precious (Precision 8) Holo
					break;
				case 6:
					//begin Punch (API) conv
					inputSampleL += (bL[c[1]] * (0.09299870608542582  - (0.00009582362368873*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.11947847710741009  - (0.00004500891602770*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.09071606264761795  + (0.00005639498984741*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.08561982770836980  - (0.00004964855606916*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.06440549220820363  + (0.00002428052139507*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.05987991812840746  + (0.00000101867082290*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.03980233135839382  + (0.00003312430049041*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03648402630896925  - (0.00002116186381142*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.01826860869525248  + (0.00003115110025396*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.01723968622495364  - (0.00002450634121718*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.00187588812316724  + (0.00002838206198968*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.00381796423957237  - (0.00003155815499462*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00852092214496733  - (0.00001702651162392*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.00315560292270588  + (0.00002547861676047*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.01258630914496868  - (0.00004555319243213*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00536435648963575  + (0.00001812393657101*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01272975658159178  - (0.00004103775306121*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00403818975172755  + (0.00003764615492871*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.01042617366897483  - (0.00003605210426041*applyconvL)));
					inputSampleL += (bL[c[20]] * (0.00126599583390057  + (0.00004305458668852*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00747876207688339  - (0.00003731207018977*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00149873689175324  - (0.00005086601800791*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00503221309488033  - (0.00003636086782783*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00342998224655821  - (0.00004103091180506*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00355585977903117  - (0.00003698982145400*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00437201792934817  - (0.00002720235666939*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00299217874451556  - (0.00004446954727956*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00457924652487249  - (0.00003859065778860*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00298182934892027  - (0.00002064710931733*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00438838441540584  - (0.00005223008424866*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00323984218794705  - (0.00003397987535887*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00407693981307314  - (0.00003935772436894*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00350435348467321  - (0.00005525463935338*applyconvL)));
					//end Punch (API) conv
					break;
				case 7:
					//begin Calibre (?) steel
					inputSampleL -= (bL[c[1]] * (0.23505923670562212  - (0.00028312859289245*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.08188436704577637  - (0.00008817721351341*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.05075798481700617  - (0.00018817166632483*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.00455811821873093  + (0.00001922902995296*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00027610521433660  - (0.00013252525469291*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.03529246280346626  - (0.00002772989223299*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01784111585586136  + (0.00010230276997291*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.04394950700298298  - (0.00005910607126944*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.01990770780547606  + (0.00007640328340556*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.04073629569741782  - (0.00007712327117090*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.01349648572795252  + (0.00005959130575917*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.03191590248003717  - (0.00008418000575151*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.00348795527924766  + (0.00005489156318238*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.02198496281481767  - (0.00008471601187581*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00504771152505089  - (0.00005525060587917*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01391075698598491  - (0.00007929630732607*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01142762504081717  - (0.00005967036737742*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00893541815021255  - (0.00007535697758141*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.01459704973464936  - (0.00005969199602841*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00694755135226282  - (0.00006930127097865*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.01516695630808575  - (0.00006365800069826*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00705917318113651  - (0.00006497209096539*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.01420501209177591  - (0.00006555654576113*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00815905656808701  - (0.00006105622534761*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.01274326525552961  - (0.00006542652857017*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00937146927845488  - (0.00006051267868722*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.01146573981165209  - (0.00006381511607749*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.01021294359409007  - (0.00005930397856398*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.01065217095323532  - (0.00006371505438319*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.01058751196699751  - (0.00006042857480233*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.01026557827762401  - (0.00006007776163871*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.01060929183604604  - (0.00006114703012726*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.01014533525058528  - (0.00005963567932887*applyconvL)));
					//end Calibre (?)
					break;
				case 8:
					//begin Tube (Manley) conv
					inputSampleL -= (bL[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvL)));
					//end Tube (Manley) conv
					break;
			}

			switch (console)
			{
				case 1:
					//begin Cider (Focusrite) MCI
					inputSampleR += (bR[c[1]] * (0.61283288942201319  + (0.00024011410669522*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.24036380659761222  - (0.00020789518206241*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.09104669761717916  + (0.00012829642741548*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.02378290768554025  - (0.00017673646470440*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.02832818490275965  - (0.00013536187747384*applyconvR)));
					inputSampleR += (bR[c[6]] * (0.03268797679215937  + (0.00015035126653359*applyconvR)));
					inputSampleR -= (bR[c[7]] * (0.04024464202655586  - (0.00015034923056735*applyconvR)));
					inputSampleR += (bR[c[8]] * (0.01864890074318696  + (0.00014513281680642*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.01632731954100322  - (0.00015509089075614*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.00318907090555589  - (0.00014784812076550*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.00208573465221869  - (0.00015350520779465*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.00907033901519614  - (0.00015442964157250*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00199458794148013  - (0.00015595640046297*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.00705979153201755  - (0.00015730069418051*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00429488975412722  - (0.00015743697943505*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.00497724878704936  - (0.00016014760011861*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00506059305562353  - (0.00016194824072466*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00483432223285621  - (0.00016329050124225*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00495100420886005  - (0.00016297509798749*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00489319520555115  - (0.00016472839684661*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00489177657970308  - (0.00016791875866630*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00487900894707044  - (0.00016755993898534*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00486234009335561  - (0.00016968157345446*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00485737490288736  - (0.00017180713324431*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00484106070563455  - (0.00017251073661092*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00483219429408410  - (0.00017321683790891*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00482013597437550  - (0.00017392186866488*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00480949628051497  - (0.00017569098775602*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00479992055604049  - (0.00017746046369449*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00478750757986987  - (0.00017745630047554*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00477828651185740  - (0.00017958043287604*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00476906544384494  - (0.00018170456527653*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00475700712413634  - (0.00018099144598088*applyconvR)));
					//end Cider (Focusrite) MCI
					break;
				case 2:
					//begin Rock (SSL) conv
					inputSampleR += (bR[c[1]] * (0.67887916185274055  + (0.00068787552301086*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.25671050678827934  + (0.00017691749454490*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.15135839896615280  + (0.00007481480365043*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.11813512969090802  + (0.00005191138121359*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.08329104347166429  + (0.00001871054659794*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.07663817456103936  + (0.00002751359071705*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.05477586152148759  + (0.00000744843212679*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.05547314737187786  + (0.00001025289931145*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.03822948356540711  - (0.00000249791561457*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.04199383340841713  - (0.00000067328840674*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.02695796542339694  - (0.00000796704606548*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.03228715059431878  - (0.00000579711816722*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.01846929689819187  - (0.00000984017804950*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.02528050435045951  - (0.00000701189792484*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.01207844846859765  - (0.00001522630289356*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01894464378378515  - (0.00001205456372080*applyconvR)));
					inputSampleR += (bR[c[17]] * (0.00667804407593324  - (0.00001343604283817*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.01408418045473130  - (0.00001246443581504*applyconvR)));
					inputSampleR += (bR[c[19]] * (0.00228696509481569  - (0.00001506764046927*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.01006277891348454  - (0.00000970723079112*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00132368373546377  + (0.00001188847238761*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00676615715578373  - (0.00001209129844861*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00426288438418556  + (0.00001286836943559*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00408897698639688  - (0.00001102542567911*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00662040619382751  + (0.00001206328529063*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00196101294183599  - (0.00000950703614981*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00845620581010342  + (0.00001279970295678*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00032595215043616  - (0.00000920518241371*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00982957737435458  + (0.00001177745362317*applyconvR)));
					inputSampleR += (bR[c[30]] * (0.00086920573760513  + (0.00000913758382404*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.01079020871452061  + (0.00000900750153697*applyconvR)));
					inputSampleR += (bR[c[32]] * (0.00167613606334460  + (0.00000732769151038*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.01138050011044332  + (0.00000946908207442*applyconvR)));
					//end Rock (SSL) conv
					break;
				case 3:
					//begin Lush (Neve) conv
					inputSampleR += (bR[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvR)));
					inputSampleR += (bR[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvR)));
					inputSampleR += (bR[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvR)));
					inputSampleR -= (bR[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvR)));
					inputSampleR += (bR[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvR)));
					inputSampleR += (bR[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvR)));
					inputSampleR += (bR[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvR)));
					inputSampleR += (bR[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvR)));
					inputSampleR += (bR[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvR)));
					inputSampleR += (bR[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvR)));
					inputSampleR += (bR[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvR)));
					inputSampleR += (bR[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvR)));
					inputSampleR += (bR[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvR)));
					inputSampleR += (bR[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvR)));
					inputSampleR += (bR[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvR)));
					inputSampleR += (bR[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvR)));
					inputSampleR += (bR[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvR)));
					inputSampleR += (bR[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvR)));
					inputSampleR += (bR[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvR)));
					inputSampleR += (bR[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvR)));
					inputSampleR += (bR[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvR)));
					inputSampleR += (bR[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvR)));
					inputSampleR += (bR[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvR)));
					//end Lush (Neve) conv
					break;
				case 4:
					//begin Elation (LA2A) vibe
					inputSampleR -= (bR[c[1]] * (0.25867935358656502  - (0.00045755657070112*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.11509367290253694  - (0.00017494270657228*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.06709853575891785  - (0.00058913102597723*applyconvR)));
					inputSampleR += (bR[c[4]] * (0.01871006356851681  - (0.00003387358004645*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00794797957360465  - (0.00044224784691203*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.01956921817394220  - (0.00006718936750076*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01682120257195205  + (0.00032857446292230*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03401069039824205  - (0.00013634182872897*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.02369950268232634  + (0.00023112685751657*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.03477071178117132  - (0.00018029792231600*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.02024369717958201  + (0.00017337813374202*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.02819087729102172  - (0.00021438538665420*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.01147946743141303  + (0.00014424066034649*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.01894777011468867  - (0.00021549146262408*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.00301370330346873  + (0.00013527460148394*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01067147835815486  - (0.00020960689910868*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00402715397506384  - (0.00014421582712470*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00502221703392005  - (0.00019805767015024*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00808788533308497  - (0.00016095444141931*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00232696588842683  - (0.00018384470981829*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00943950821324531  - (0.00017098987347593*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00193709517200834  - (0.00018151995939591*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00899713952612659  - (0.00017385835059948*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00280584331659089  - (0.00017742164162470*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00780381001954970  - (0.00018002500755708*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00400370310490333  - (0.00017471691087957*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00661527728186928  - (0.00018137323370347*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00496545526864518  - (0.00017681872601767*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00580728820997532  - (0.00018186220389790*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00549309984725666  - (0.00017722985399075*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00542194777529239  - (0.00018486900185338*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00565992080998939  - (0.00018005824393118*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00532121562846656  - (0.00018643189636216*applyconvR)));
					//end Elation (LA2A)
					break;
				case 5:
					//begin Precious (Precision 8) Holo
					inputSampleR += (bR[c[1]] * (0.59188440274551890  - (0.00008361469668405*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.24439750948076133  + (0.00002651678396848*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.14109876103205621  - (0.00000840487181372*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.10053507128157971  + (0.00001768100964598*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.05859287880626238  - (0.00000361398065989*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.04337406889823660  + (0.00000735941182117*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01589900680531097  + (0.00000207347387987*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.01087234854973281  + (0.00000732123412029*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.00845782429679176  - (0.00000133058605071*applyconvR)));
					inputSampleR += (bR[c[10]] * (0.00662278586618295  - (0.00000424594730611*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.02000592193760155  + (0.00000632896879068*applyconvR)));
					inputSampleR += (bR[c[12]] * (0.01321157777167565  - (0.00001421171592570*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.02249955362988238  + (0.00000163937127317*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.01196492077581504  - (0.00000535385220676*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.01905917427000097  + (0.00000121672882030*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00761909482108073  - (0.00000326242895115*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01362744780256239  + (0.00000359274216003*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00200183122683721  - (0.00000089207452791*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00833042637239315  + (0.00000946767677294*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00258481175207224  - (0.00000087429351464*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00459744479712244  - (0.00000049519758701*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00534277030993820  + (0.00000397547847155*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00272332919605675  + (0.00000040077229097*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00637243782359372  - (0.00000139419072176*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00233001590327504  + (0.00000420129915747*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00623296727793041  + (0.00000019010664856*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00276177096376805  + (0.00000580301901385*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00559184754866264  + (0.00000080597287792*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00343180144395919  - (0.00000243701142085*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00493325428861701  + (0.00000300985740900*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00396140827680823  - (0.00000051459681789*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00448497879902493  + (0.00000744412841743*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00425146888772076  - (0.00000082346016542*applyconvR)));
					//end Precious (Precision 8) Holo
					break;
				case 6:
					//begin Punch (API) conv
					inputSampleR += (bR[c[1]] * (0.09299870608542582  - (0.00009582362368873*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.11947847710741009  - (0.00004500891602770*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.09071606264761795  + (0.00005639498984741*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.08561982770836980  - (0.00004964855606916*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.06440549220820363  + (0.00002428052139507*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.05987991812840746  + (0.00000101867082290*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.03980233135839382  + (0.00003312430049041*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03648402630896925  - (0.00002116186381142*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.01826860869525248  + (0.00003115110025396*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.01723968622495364  - (0.00002450634121718*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.00187588812316724  + (0.00002838206198968*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.00381796423957237  - (0.00003155815499462*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00852092214496733  - (0.00001702651162392*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.00315560292270588  + (0.00002547861676047*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.01258630914496868  - (0.00004555319243213*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00536435648963575  + (0.00001812393657101*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01272975658159178  - (0.00004103775306121*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00403818975172755  + (0.00003764615492871*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.01042617366897483  - (0.00003605210426041*applyconvR)));
					inputSampleR += (bR[c[20]] * (0.00126599583390057  + (0.00004305458668852*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00747876207688339  - (0.00003731207018977*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00149873689175324  - (0.00005086601800791*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00503221309488033  - (0.00003636086782783*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00342998224655821  - (0.00004103091180506*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00355585977903117  - (0.00003698982145400*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00437201792934817  - (0.00002720235666939*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00299217874451556  - (0.00004446954727956*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00457924652487249  - (0.00003859065778860*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00298182934892027  - (0.00002064710931733*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00438838441540584  - (0.00005223008424866*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00323984218794705  - (0.00003397987535887*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00407693981307314  - (0.00003935772436894*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00350435348467321  - (0.00005525463935338*applyconvR)));
					//end Punch (API) conv
					break;
				case 7:
					//begin Calibre (?) steel
					inputSampleR -= (bR[c[1]] * (0.23505923670562212  - (0.00028312859289245*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.08188436704577637  - (0.00008817721351341*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.05075798481700617  - (0.00018817166632483*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.00455811821873093  + (0.00001922902995296*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00027610521433660  - (0.00013252525469291*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.03529246280346626  - (0.00002772989223299*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01784111585586136  + (0.00010230276997291*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.04394950700298298  - (0.00005910607126944*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.01990770780547606  + (0.00007640328340556*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.04073629569741782  - (0.00007712327117090*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.01349648572795252  + (0.00005959130575917*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.03191590248003717  - (0.00008418000575151*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.00348795527924766  + (0.00005489156318238*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.02198496281481767  - (0.00008471601187581*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00504771152505089  - (0.00005525060587917*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01391075698598491  - (0.00007929630732607*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01142762504081717  - (0.00005967036737742*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00893541815021255  - (0.00007535697758141*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.01459704973464936  - (0.00005969199602841*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00694755135226282  - (0.00006930127097865*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.01516695630808575  - (0.00006365800069826*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00705917318113651  - (0.00006497209096539*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.01420501209177591  - (0.00006555654576113*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00815905656808701  - (0.00006105622534761*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.01274326525552961  - (0.00006542652857017*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00937146927845488  - (0.00006051267868722*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.01146573981165209  - (0.00006381511607749*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.01021294359409007  - (0.00005930397856398*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.01065217095323532  - (0.00006371505438319*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.01058751196699751  - (0.00006042857480233*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.01026557827762401  - (0.00006007776163871*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.01060929183604604  - (0.00006114703012726*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.01014533525058528  - (0.00005963567932887*applyconvR)));
					//end Calibre (?)
					break;
				case 8:
					//begin Tube (Manley) conv
					inputSampleR -= (bR[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvR)));
					//end Tube (Manley) conv
					break;
			}

			bc_bridgerectifier = fabs(inputSampleL);
			bc_bridgerectifier = 1.0-cos(bc_bridgerectifier);
			if (inputSampleL > 0) inputSampleL -= bc_bridgerectifier;
			else inputSampleL += bc_bridgerectifier;

			bc_bridgerectifier = fabs(inputSampleR);
			bc_bridgerectifier = 1.0-cos(bc_bridgerectifier);
			if (inputSampleR > 0) inputSampleR -= bc_bridgerectifier;
			else inputSampleR += bc_bridgerectifier;


			if (bc_outgain != 1.0) {
				inputSampleL *= bc_outgain;
				inputSampleR *= bc_outgain;
			}

			if (bc_wet !=1.0) {
				inputSampleL = (inputSampleL * bc_wet) + (bc_drySampleL * (1.0-bc_wet));
				inputSampleR = (inputSampleR * bc_wet) + (bc_drySampleR * (1.0-bc_wet));
			}
		}

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

		in1++;
		in2++;
		out1++;
		out2++;
    }
}

void ConsoleZBuss::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames)
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

	VstInt32 inFramesToProcess = sampleFrames; //vst doesn't give us this as a separate variable so we'll make it
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	// ConsoleLABuss

	gainA = gainB;
	gainB = sqrt(A); //smoothed master fader from Z2 filters
	//this will be applied three times: this is to make the various tone alterations
	//hit differently at different master fader drive levels.
	//in particular, backing off the master fader tightens the super lows
	//but opens up the modified Sinew, because more of the attentuation happens before
	//you even get to slew clipping :) and if the fader is not active, it bypasses completely.

	double threshSinew = 0.718/overallscale;
	double subTrim = 0.0011 / overallscale;

	// Desk4

	double gain = (pow(B,2)*10)+0.0001;
	double gaintrim = (pow(B,2)*2)+1.0;
	double slewgain = (pow(C,3)*40)+0.0001;
	double prevslew = 0.105;
	double intensity = (pow(D,6)*15)+0.0001;
	double depthA = (pow(E,4)*940)+0.00001;
	int offsetA = (int)(depthA * overallscale);
	if (offsetA < 1) offsetA = 1;
	if (offsetA > 4880) offsetA = 4880;
	double balanceB = 0.0001;
	slewgain *= overallscale;
	prevslew *= overallscale;
	balanceB /= overallscale;
	double desk4_outputgain = F;
	double wet = G;
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
	double drySampleL;
	double drySampleR;
	double inputSampleL;
	double inputSampleR;

	// Buttercomp2

	double inputgain = pow(10.0,(H*14.0)/20.0);
	double compfactor = 0.012 * (H / 135.0);
	double output = I * 2.0;
	double butter_wet = J;
	//removed extra dry variable
	double outputgain = inputgain;
	outputgain -= 1.0;
	outputgain /= 1.5;
	outputgain += 1.0;

	// VariMu

	double threshold = 1.001 - (1.0-pow(1.0-K,3));
	double muMakeupGain = sqrt(1.0 / threshold);
	muMakeupGain = (muMakeupGain + sqrt(muMakeupGain))/2.0;
	muMakeupGain = sqrt(muMakeupGain);
	double outGain = sqrt(muMakeupGain);
	//gain settings around threshold
	double release = pow((1.15-L),5)*32768.0;
	release /= overallscale;
	double fastest = sqrt(release);
	//speed settings around release
	double coefficient;
	double mu_output = outGain * M;
	double mu_wet = N;
	double squaredSampleL;
	double squaredSampleR;

	// BussColors4

	const int maxConvolutionBufferSize = (int)(34.0 * overallscale); //we won't use more of the buffer than we have to
	for (int count = 0; count < 34; count++) c[count] = (int)(count * overallscale); //assign conv taps
	double bc_drySampleL;
	double bc_drySampleR;
	double applyconvL;
	double applyconvR;
	int bc_offsetA = 3;
	double dynamicconvL = 3.0;
	double dynamicconvR = 3.0;
	double bc_gain = 0.436;
	double bc_outgain = 1.0;
	double bc_bridgerectifier;

	double bc_inputSampleL;
	double bc_inputSampleR;

	int console = (int)( O * 7.999 )+1; //the AU used a 1-8 index, will just keep it
	switch (console)
	{
		case 1: bc_offsetA = 4; bc_gain = g[1]; bc_outgain = outg[1]; break; //Dark (Focusrite)
		case 2: bc_offsetA = 3; bc_gain = g[2]; bc_outgain = outg[2]; break; //Rock (SSL)
		case 3: bc_offsetA = 5; bc_gain = g[3]; bc_outgain = outg[3]; break; //Lush (Neve)
		case 4: bc_offsetA = 8; bc_gain = g[4]; bc_outgain = outg[4]; break; //Vibe (Elation)
		case 5: bc_offsetA = 5; bc_gain = g[5]; bc_outgain = outg[5]; break; //Holo (Precision 8)
		case 6: bc_offsetA = 7; bc_gain = g[6]; bc_outgain = outg[6]; break; //Punch (API)
		case 7: bc_offsetA = 7; bc_gain = g[7]; bc_outgain = outg[7]; break; //Steel (Calibre)
		case 8: bc_offsetA = 6; bc_gain = g[8]; bc_outgain = outg[8]; break; //Tube (Manley)
	}
	bc_offsetA = (int)(bc_offsetA * overallscale); //we extend the sag buffer too, at high sample rates
	if (bc_offsetA < 3) bc_offsetA = 3; //are we getting divides by zero?
	bc_gain *= pow(10.0,((P * 36.0)-18.0)/14.0); //add adjustment factor
	bc_outgain *= pow(10.0,(((Q * 36.0)-18.0)+3.3)/14.0); //add adjustment factor
	double bc_wet = R;
	//removed extra dry variable


    while (--sampleFrames >= 0)
    {
		double inputSampleL = *in1;
		double inputSampleR = *in2;

		if (fabs(inputSampleL)<1.18e-23) inputSampleL = fpdL * 1.18e-17;
		if (fabs(inputSampleR)<1.18e-23) inputSampleR = fpdR * 1.18e-17;

		// Desk4

		if(wet > 0.0) {
			drySampleL = inputSampleL;
			drySampleR = inputSampleR;

			if (desk4_gcount < 0 || desk4_gcount > 4900) {desk4_gcount = 4900;}

			desk4_dL[desk4_gcount+4900] = desk4_dL[desk4_gcount] = fabs(inputSampleL)*intensity;
			desk4_controlL += (desk4_dL[desk4_gcount] / offsetA);
			desk4_controlL -= (desk4_dL[desk4_gcount+offsetA] / offsetA);
			desk4_controlL -= 0.000001;
			clampL = 1;
			if (desk4_controlL < 0) {desk4_controlL = 0;}
			if (desk4_controlL > 1) {clampL -= (desk4_controlL - 1); desk4_controlL = 1;}
			if (clampL < 0.5) {clampL = 0.5;}

			desk4_dR[desk4_gcount+4900] = desk4_dR[desk4_gcount] = fabs(inputSampleR)*intensity;
			desk4_controlR += (desk4_dR[desk4_gcount] / offsetA);
			desk4_controlR -= (desk4_dR[desk4_gcount+offsetA] / offsetA);
			desk4_controlR -= 0.000001;
			clampR = 1;
			if (desk4_controlR < 0) {desk4_controlR = 0;}
			if (desk4_controlR > 1) {clampR -= (desk4_controlR - 1); desk4_controlR = 1;}
			if (clampR < 0.5) {clampR = 0.5;}


			desk4_gcount--;
			//control = 0 to 1
			thicknessL = ((1.0 - desk4_controlL) * 2.0) - 1.0;
			thicknessR = ((1.0 - desk4_controlR) * 2.0) - 1.0;

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

			inputSampleL *= gain;
			bridgerectifier = fabs(inputSampleL);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleL > 0) inputSampleL = bridgerectifier;
			else inputSampleL = -bridgerectifier;
			//drive section
			inputSampleL /= gain;
			inputSampleL *= gaintrim;
			//end of Desk section

			inputSampleR *= gain;
			bridgerectifier = fabs(inputSampleR);
			if (bridgerectifier > 1.57079633) bridgerectifier = 1.0;
			else bridgerectifier = sin(bridgerectifier);
			if (inputSampleR > 0) inputSampleR = bridgerectifier;
			else inputSampleR = -bridgerectifier;
			//drive section
			inputSampleR /= gain;
			inputSampleR *= gaintrim;
			//end of Desk section

			if (desk4_outputgain != 1.0) {
				inputSampleL *= desk4_outputgain;
				inputSampleR *= desk4_outputgain;
			}

			if (wet !=1.0) {
				inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0-wet));
				inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0-wet));
			}
		}

		// ConsoleLABuss

		double temp = (double)sampleFrames/inFramesToProcess;
		double la_gain = (gainA*temp)+(gainB*(1.0-temp));
		//setting up smoothed master fader

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
		inputSampleL -= (subSampleL*16.0);
		inputSampleR -= (subSampleR*16.0);
		//end SubTight section

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//subtight is always fully engaged: tighten response when restraining full console

		//begin ConsoleZeroBuss which is the one we choose for ConsoleZ
		if (inputSampleL > 2.8) inputSampleL = 2.8;
		if (inputSampleL < -2.8) inputSampleL = -2.8;
		if (inputSampleL > 0.0) inputSampleL = (inputSampleL*2.0)/(3.0-inputSampleL);
		else inputSampleL = -(inputSampleL*-2.0)/(3.0+inputSampleL);

		if (inputSampleR > 2.8) inputSampleR = 2.8;
		if (inputSampleR < -2.8) inputSampleR = -2.8;
		if (inputSampleR > 0.0) inputSampleR = (inputSampleR*2.0)/(3.0-inputSampleR);
		else inputSampleR = -(inputSampleR*-2.0)/(3.0+inputSampleR);
		//ConsoleZero Buss

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//after C0Buss but before EverySlew: allow highs to come out a bit more
		//when pulling back master fader. Less drive equals more open

		temp = inputSampleL;
		double clamp = inputSampleL - lastSinewL;
		if (lastSinewL > 1.0) lastSinewL = 1.0;
		if (lastSinewL < -1.0) lastSinewL = -1.0;
		double sinew = threshSinew * cos(lastSinewL);
		if (clamp > sinew) temp = lastSinewL + sinew;
		if (-clamp > sinew) temp = lastSinewL - sinew;
		inputSampleL = lastSinewL = temp;
		temp = inputSampleR;
		clamp = inputSampleR - lastSinewR;
		if (lastSinewR > 1.0) lastSinewR = 1.0;
		if (lastSinewR < -1.0) lastSinewR = -1.0;
		sinew = threshSinew * cos(lastSinewR);
		if (clamp > sinew) temp = lastSinewR + sinew;
		if (-clamp > sinew) temp = lastSinewR - sinew;
		inputSampleR = lastSinewR = temp;

		if (la_gain < 1.0) {
			inputSampleL *= la_gain;
			inputSampleR *= la_gain;
		} //if using the master fader, we are going to attenuate three places.
		//after EverySlew fades the total output sound: least change in tone here.

		// Buttercomp2

		if(butter_wet > 0.0) {

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
			//effectively digital black, we'll subtract it aButterComp2. We want a 'air' hiss
			double butter_drySampleL = inputSampleL;
			double butter_drySampleR = inputSampleR;

			inputSampleL *= inputgain;
			inputSampleR *= inputgain;

			double divisor = compfactor / (1.0+fabs(lastOutputL));
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
				if (butter_flip)
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
				if (butter_flip)
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
			//this causes each of the four to update only when active and in the correct 'butter_flip'

			divisor = compfactor / (1.0+fabs(lastOutputR));
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
				if (butter_flip)
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
				if (butter_flip)
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
			//this causes each of the four to update only when active and in the correct 'butter_flip'

			double totalmultiplierL;
			double totalmultiplierR;
			if (butter_flip)
			{
				totalmultiplierL = (controlAposL * outputposL) + (controlAnegL * outputnegL);
				totalmultiplierR = (controlAposR * outputposR) + (controlAnegR * outputnegR);
			}
			else
			{
				totalmultiplierL = (controlBposL * outputposL) + (controlBnegL * outputnegL);
				totalmultiplierR = (controlBposR * outputposR) + (controlBnegR * outputnegR);
			}
			//this combines the sides according to butter_flip, blending relative to the input value

			inputSampleL *= totalmultiplierL;
			inputSampleL /= outputgain;

			inputSampleR *= totalmultiplierR;
			inputSampleR /= outputgain;

			if (output != 1.0) {
				inputSampleL *= output;
				inputSampleR *= output;
			}

			if (butter_wet !=1.0) {
				inputSampleL = (inputSampleL * butter_wet) + (butter_drySampleL * (1.0-butter_wet));
				inputSampleR = (inputSampleR * butter_wet) + (butter_drySampleR * (1.0-butter_wet));
			}

			lastOutputL = inputSampleL;
			lastOutputR = inputSampleR;
			//we will make this factor respond to use of dry/wet

			butter_flip = !butter_flip;
		}

		// VariMu
		if(mu_wet > 0.0) {
			static int mu_noisesourceL = 0;
			static int mu_noisesourceR = 850010;
			int mu_residue;
			double mu_applyresidue;

			mu_noisesourceL = mu_noisesourceL % 1700021; mu_noisesourceL++;
			mu_residue = mu_noisesourceL * mu_noisesourceL;
			mu_residue = mu_residue % 170003; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17011; mu_residue *= mu_residue;
			mu_residue = mu_residue % 1709; mu_residue *= mu_residue;
			mu_residue = mu_residue % 173; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17;
			mu_applyresidue = mu_residue;
			mu_applyresidue *= 0.00000001;
			mu_applyresidue *= 0.00000001;
			inputSampleL += mu_applyresidue;
			if (inputSampleL<1.2e-38 && -inputSampleL<1.2e-38) {
				inputSampleL -= mu_applyresidue;
			}

			mu_noisesourceR = mu_noisesourceR % 1700021; mu_noisesourceR++;
			mu_residue = mu_noisesourceR * mu_noisesourceR;
			mu_residue = mu_residue % 170003; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17011; mu_residue *= mu_residue;
			mu_residue = mu_residue % 1709; mu_residue *= mu_residue;
			mu_residue = mu_residue % 173; mu_residue *= mu_residue;
			mu_residue = mu_residue % 17;
			mu_applyresidue = mu_residue;
			mu_applyresidue *= 0.00000001;
			mu_applyresidue *= 0.00000001;
			inputSampleR += mu_applyresidue;
			if (inputSampleR<1.2e-38 && -inputSampleR<1.2e-38) {
				inputSampleR -= mu_applyresidue;
			}
			//for live air, we always apply the dither noise. Then, if our result is
			//effectively digital black, we'll subtract it aVariMu. We want a 'air' hiss
			double mu_drySampleL = inputSampleL;
			double mu_drySampleR = inputSampleR;

			if (fabs(inputSampleL) > fabs(previousL)) squaredSampleL = previousL * previousL;
			else squaredSampleL = inputSampleL * inputSampleL;
			previousL = inputSampleL;
			inputSampleL *= muMakeupGain;

			if (fabs(inputSampleR) > fabs(previousR)) squaredSampleR = previousR * previousR;
			else squaredSampleR = inputSampleR * inputSampleR;
			previousR = inputSampleR;
			inputSampleR *= muMakeupGain;

			//adjust coefficients for L
			if (flip)
			{
				if (fabs(squaredSampleL) > threshold)
				{
					muVaryL = threshold / fabs(squaredSampleL);
					muAttackL = sqrt(fabs(muSpeedAL));
					muCoefficientAL = muCoefficientAL * (muAttackL-1.0);
					if (muVaryL < threshold)
					{
						muCoefficientAL = muCoefficientAL + threshold;
					}
					else
					{
						muCoefficientAL = muCoefficientAL + muVaryL;
					}
					muCoefficientAL = muCoefficientAL / muAttackL;
				}
				else
				{
					muCoefficientAL = muCoefficientAL * ((muSpeedAL * muSpeedAL)-1.0);
					muCoefficientAL = muCoefficientAL + 1.0;
					muCoefficientAL = muCoefficientAL / (muSpeedAL * muSpeedAL);
				}
				muNewSpeedL = muSpeedAL * (muSpeedAL-1);
				muNewSpeedL = muNewSpeedL + fabs(squaredSampleL*release)+fastest;
				muSpeedAL = muNewSpeedL / muSpeedAL;
			}
			else
			{
				if (fabs(squaredSampleL) > threshold)
				{
					muVaryL = threshold / fabs(squaredSampleL);
					muAttackL = sqrt(fabs(muSpeedBL));
					muCoefficientBL = muCoefficientBL * (muAttackL-1);
					if (muVaryL < threshold)
					{
						muCoefficientBL = muCoefficientBL + threshold;
					}
					else
					{
						muCoefficientBL = muCoefficientBL + muVaryL;
					}
					muCoefficientBL = muCoefficientBL / muAttackL;
				}
				else
				{
					muCoefficientBL = muCoefficientBL * ((muSpeedBL * muSpeedBL)-1.0);
					muCoefficientBL = muCoefficientBL + 1.0;
					muCoefficientBL = muCoefficientBL / (muSpeedBL * muSpeedBL);
				}
				muNewSpeedL = muSpeedBL * (muSpeedBL-1);
				muNewSpeedL = muNewSpeedL + fabs(squaredSampleL*release)+fastest;
				muSpeedBL = muNewSpeedL / muSpeedBL;
			}
			//got coefficients, adjusted speeds for L

			//adjust coefficients for R
			if (flip)
			{
				if (fabs(squaredSampleR) > threshold)
				{
					muVaryR = threshold / fabs(squaredSampleR);
					muAttackR = sqrt(fabs(muSpeedAR));
					muCoefficientAR = muCoefficientAR * (muAttackR-1.0);
					if (muVaryR < threshold)
					{
						muCoefficientAR = muCoefficientAR + threshold;
					}
					else
					{
						muCoefficientAR = muCoefficientAR + muVaryR;
					}
					muCoefficientAR = muCoefficientAR / muAttackR;
				}
				else
				{
					muCoefficientAR = muCoefficientAR * ((muSpeedAR * muSpeedAR)-1.0);
					muCoefficientAR = muCoefficientAR + 1.0;
					muCoefficientAR = muCoefficientAR / (muSpeedAR * muSpeedAR);
				}
				muNewSpeedR = muSpeedAR * (muSpeedAR-1);
				muNewSpeedR = muNewSpeedR + fabs(squaredSampleR*release)+fastest;
				muSpeedAR = muNewSpeedR / muSpeedAR;
			}
			else
			{
				if (fabs(squaredSampleR) > threshold)
				{
					muVaryR = threshold / fabs(squaredSampleR);
					muAttackR = sqrt(fabs(muSpeedBR));
					muCoefficientBR = muCoefficientBR * (muAttackR-1);
					if (muVaryR < threshold)
					{
						muCoefficientBR = muCoefficientBR + threshold;
					}
					else
					{
						muCoefficientBR = muCoefficientBR + muVaryR;
					}
					muCoefficientBR = muCoefficientBR / muAttackR;
				}
				else
				{
					muCoefficientBR = muCoefficientBR * ((muSpeedBR * muSpeedBR)-1.0);
					muCoefficientBR = muCoefficientBR + 1.0;
					muCoefficientBR = muCoefficientBR / (muSpeedBR * muSpeedBR);
				}
				muNewSpeedR = muSpeedBR * (muSpeedBR-1);
				muNewSpeedR = muNewSpeedR + fabs(squaredSampleR*release)+fastest;
				muSpeedBR = muNewSpeedR / muSpeedBR;
			}
			//got coefficients, adjusted speeds for R

			if (flip)
			{
				coefficient = (muCoefficientAL + pow(muCoefficientAL,2))/2.0;
				inputSampleL *= coefficient;
				coefficient = (muCoefficientAR + pow(muCoefficientAR,2))/2.0;
				inputSampleR *= coefficient;
			}
			else
			{
				coefficient = (muCoefficientBL + pow(muCoefficientBL,2))/2.0;
				inputSampleL *= coefficient;
				coefficient = (muCoefficientBR + pow(muCoefficientBR,2))/2.0;
				inputSampleR *= coefficient;
			}
			//applied compression with vari-vari-µ-µ-µ-µ-µ-µ-is-the-kitten-song o/~
			//applied gain correction to control output level- tends to constrain sound rather than inflate it
			flip = !flip;

			if (mu_output < 1.0) {
				inputSampleL *= mu_output;
				inputSampleR *= mu_output;
			}
			if (mu_wet < 1.0) {
				inputSampleL = (mu_drySampleL * (1.0-mu_wet)) + (inputSampleL * mu_wet);
				inputSampleR = (mu_drySampleR * (1.0-mu_wet)) + (inputSampleR * mu_wet);
			}
			//nice little output stage template: if we have another scale of floating point
			//number, we really don't want to meaninglessly multiply that by 1.0.
		}
		// BussColors4

		if(bc_wet > 0.0) {

			bc_drySampleL = inputSampleL;
			bc_drySampleR = inputSampleR;

			if (bc_gain != 1.0) {
				inputSampleL *= bc_gain;
				inputSampleR *= bc_gain;
			}


			bc_bridgerectifier = fabs(inputSampleL);
			slowdynL *= 0.999;
			slowdynL += bc_bridgerectifier;
			if (slowdynL > 1.5) {slowdynL = 1.5;}
			//before the iron bar- fry console for crazy behavior
			dynamicconvL = 2.5 + slowdynL;

			if (bc_bridgerectifier > 1.57079633) bc_bridgerectifier = 1.0;
			else bc_bridgerectifier = sin(bc_bridgerectifier);
			if (inputSampleL > 0) inputSampleL = bc_bridgerectifier;
			else inputSampleL = -bc_bridgerectifier;
			//end pre saturation stage L

			bc_bridgerectifier = fabs(inputSampleR);
			slowdynR *= 0.999;
			slowdynR += bc_bridgerectifier;
			if (slowdynR > 1.5) {slowdynR = 1.5;}
			//before the iron bar- fry console for crazy behavior
			dynamicconvR = 2.5 + slowdynR;

			if (bc_bridgerectifier > 1.57079633) bc_bridgerectifier = 1.0;
			else bc_bridgerectifier = sin(bc_bridgerectifier);
			if (inputSampleR > 0) inputSampleR = bc_bridgerectifier;
			else inputSampleR = -bc_bridgerectifier;
			//end pre saturation stage R

			if (gcount < 0 || gcount > 44) gcount = 44;

			dL[gcount+44] = dL[gcount] = fabs(inputSampleL);
			controlL += (dL[gcount] / bc_offsetA);
			controlL -= (dL[gcount+bc_offsetA] / bc_offsetA);
			controlL -= 0.000001;
			if (controlL < 0) controlL = 0;
			if (controlL > 100) controlL = 100;
			applyconvL = (controlL / bc_offsetA) * dynamicconvL;
			//now we have a 'sag' style average to apply to the conv, L

			dR[gcount+44] = dR[gcount] = fabs(inputSampleR);
			controlR += (dR[gcount] / bc_offsetA);
			controlR -= (dR[gcount+bc_offsetA] / bc_offsetA);
			controlR -= 0.000001;
			if (controlR < 0) controlR = 0;
			if (controlR > 100) controlR = 100;
			applyconvR = (controlR / bc_offsetA) * dynamicconvR;
			//now we have a 'sag' style average to apply to the conv, R

			gcount--;

			//now the convolution
			for (int count = maxConvolutionBufferSize; count > 0; --count) {bL[count] = bL[count-1];} //was 173
			//we're only doing assigns, and it saves us an add inside the convolution calculation
			//therefore, we'll just assign everything one step along and have our buffer that way.
			bL[0] = inputSampleL;

			for (int count = maxConvolutionBufferSize; count > 0; --count) {bR[count] = bR[count-1];} //was 173
			//we're only doing assigns, and it saves us an add inside the convolution calculation
			//therefore, we'll just assign everything one step along and have our buffer that way.
			bR[0] = inputSampleR;
			//The reason these are separate is, since it's just a big assign-fest, it's possible compilers can
			//come up with a clever way to do it. Interleaving the samples might make it enough more complicated
			//that the compiler wouldn't know to do 'em in batches or whatever. Just a thought, profiling would
			//be the correct way to find out this (or indeed, whether doing another add insode the convolutions would
			//be the best bet,

			//The convolutions!

			switch (console)
			{
				case 1:
					//begin Cider (Focusrite) MCI
					inputSampleL += (bL[c[1]] * (0.61283288942201319  + (0.00024011410669522*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.24036380659761222  - (0.00020789518206241*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.09104669761717916  + (0.00012829642741548*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.02378290768554025  - (0.00017673646470440*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.02832818490275965  - (0.00013536187747384*applyconvL)));
					inputSampleL += (bL[c[6]] * (0.03268797679215937  + (0.00015035126653359*applyconvL)));
					inputSampleL -= (bL[c[7]] * (0.04024464202655586  - (0.00015034923056735*applyconvL)));
					inputSampleL += (bL[c[8]] * (0.01864890074318696  + (0.00014513281680642*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.01632731954100322  - (0.00015509089075614*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.00318907090555589  - (0.00014784812076550*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.00208573465221869  - (0.00015350520779465*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.00907033901519614  - (0.00015442964157250*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00199458794148013  - (0.00015595640046297*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.00705979153201755  - (0.00015730069418051*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00429488975412722  - (0.00015743697943505*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.00497724878704936  - (0.00016014760011861*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00506059305562353  - (0.00016194824072466*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00483432223285621  - (0.00016329050124225*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00495100420886005  - (0.00016297509798749*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00489319520555115  - (0.00016472839684661*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00489177657970308  - (0.00016791875866630*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00487900894707044  - (0.00016755993898534*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00486234009335561  - (0.00016968157345446*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00485737490288736  - (0.00017180713324431*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00484106070563455  - (0.00017251073661092*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00483219429408410  - (0.00017321683790891*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00482013597437550  - (0.00017392186866488*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00480949628051497  - (0.00017569098775602*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00479992055604049  - (0.00017746046369449*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00478750757986987  - (0.00017745630047554*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00477828651185740  - (0.00017958043287604*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00476906544384494  - (0.00018170456527653*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00475700712413634  - (0.00018099144598088*applyconvL)));
					//end Cider (Focusrite) MCI
					break;
				case 2:
					//begin Rock (SSL) conv
					inputSampleL += (bL[c[1]] * (0.67887916185274055  + (0.00068787552301086*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.25671050678827934  + (0.00017691749454490*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.15135839896615280  + (0.00007481480365043*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.11813512969090802  + (0.00005191138121359*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.08329104347166429  + (0.00001871054659794*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.07663817456103936  + (0.00002751359071705*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.05477586152148759  + (0.00000744843212679*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.05547314737187786  + (0.00001025289931145*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.03822948356540711  - (0.00000249791561457*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.04199383340841713  - (0.00000067328840674*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.02695796542339694  - (0.00000796704606548*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.03228715059431878  - (0.00000579711816722*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.01846929689819187  - (0.00000984017804950*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.02528050435045951  - (0.00000701189792484*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.01207844846859765  - (0.00001522630289356*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01894464378378515  - (0.00001205456372080*applyconvL)));
					inputSampleL += (bL[c[17]] * (0.00667804407593324  - (0.00001343604283817*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.01408418045473130  - (0.00001246443581504*applyconvL)));
					inputSampleL += (bL[c[19]] * (0.00228696509481569  - (0.00001506764046927*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.01006277891348454  - (0.00000970723079112*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00132368373546377  + (0.00001188847238761*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00676615715578373  - (0.00001209129844861*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00426288438418556  + (0.00001286836943559*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00408897698639688  - (0.00001102542567911*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00662040619382751  + (0.00001206328529063*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00196101294183599  - (0.00000950703614981*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00845620581010342  + (0.00001279970295678*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00032595215043616  - (0.00000920518241371*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00982957737435458  + (0.00001177745362317*applyconvL)));
					inputSampleL += (bL[c[30]] * (0.00086920573760513  + (0.00000913758382404*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.01079020871452061  + (0.00000900750153697*applyconvL)));
					inputSampleL += (bL[c[32]] * (0.00167613606334460  + (0.00000732769151038*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.01138050011044332  + (0.00000946908207442*applyconvL)));
					//end Rock (SSL) conv
					break;
				case 3:
					//begin Lush (Neve) conv
					inputSampleL += (bL[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvL)));
					inputSampleL += (bL[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvL)));
					inputSampleL += (bL[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvL)));
					inputSampleL -= (bL[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvL)));
					inputSampleL += (bL[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvL)));
					inputSampleL += (bL[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvL)));
					inputSampleL += (bL[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvL)));
					inputSampleL += (bL[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvL)));
					inputSampleL += (bL[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvL)));
					inputSampleL += (bL[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvL)));
					inputSampleL += (bL[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvL)));
					inputSampleL += (bL[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvL)));
					inputSampleL += (bL[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvL)));
					inputSampleL += (bL[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvL)));
					inputSampleL += (bL[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvL)));
					inputSampleL += (bL[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvL)));
					inputSampleL += (bL[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvL)));
					inputSampleL += (bL[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvL)));
					inputSampleL += (bL[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvL)));
					inputSampleL += (bL[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvL)));
					inputSampleL += (bL[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvL)));
					inputSampleL += (bL[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvL)));
					inputSampleL += (bL[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvL)));
					//end Lush (Neve) conv
					break;
				case 4:
					//begin Elation (LA2A) vibe
					inputSampleL -= (bL[c[1]] * (0.25867935358656502  - (0.00045755657070112*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.11509367290253694  - (0.00017494270657228*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.06709853575891785  - (0.00058913102597723*applyconvL)));
					inputSampleL += (bL[c[4]] * (0.01871006356851681  - (0.00003387358004645*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00794797957360465  - (0.00044224784691203*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.01956921817394220  - (0.00006718936750076*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01682120257195205  + (0.00032857446292230*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03401069039824205  - (0.00013634182872897*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.02369950268232634  + (0.00023112685751657*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.03477071178117132  - (0.00018029792231600*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.02024369717958201  + (0.00017337813374202*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.02819087729102172  - (0.00021438538665420*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.01147946743141303  + (0.00014424066034649*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.01894777011468867  - (0.00021549146262408*applyconvL)));
					inputSampleL += (bL[c[15]] * (0.00301370330346873  + (0.00013527460148394*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01067147835815486  - (0.00020960689910868*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00402715397506384  - (0.00014421582712470*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00502221703392005  - (0.00019805767015024*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00808788533308497  - (0.00016095444141931*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00232696588842683  - (0.00018384470981829*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00943950821324531  - (0.00017098987347593*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00193709517200834  - (0.00018151995939591*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00899713952612659  - (0.00017385835059948*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00280584331659089  - (0.00017742164162470*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00780381001954970  - (0.00018002500755708*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00400370310490333  - (0.00017471691087957*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00661527728186928  - (0.00018137323370347*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00496545526864518  - (0.00017681872601767*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00580728820997532  - (0.00018186220389790*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00549309984725666  - (0.00017722985399075*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00542194777529239  - (0.00018486900185338*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00565992080998939  - (0.00018005824393118*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00532121562846656  - (0.00018643189636216*applyconvL)));
					//end Elation (LA2A)
					break;
				case 5:
					//begin Precious (Precision 8) Holo
					inputSampleL += (bL[c[1]] * (0.59188440274551890  - (0.00008361469668405*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.24439750948076133  + (0.00002651678396848*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.14109876103205621  - (0.00000840487181372*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.10053507128157971  + (0.00001768100964598*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.05859287880626238  - (0.00000361398065989*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.04337406889823660  + (0.00000735941182117*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01589900680531097  + (0.00000207347387987*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.01087234854973281  + (0.00000732123412029*applyconvL)));
					inputSampleL -= (bL[c[9]] * (0.00845782429679176  - (0.00000133058605071*applyconvL)));
					inputSampleL += (bL[c[10]] * (0.00662278586618295  - (0.00000424594730611*applyconvL)));
					inputSampleL -= (bL[c[11]] * (0.02000592193760155  + (0.00000632896879068*applyconvL)));
					inputSampleL += (bL[c[12]] * (0.01321157777167565  - (0.00001421171592570*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.02249955362988238  + (0.00000163937127317*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.01196492077581504  - (0.00000535385220676*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.01905917427000097  + (0.00000121672882030*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00761909482108073  - (0.00000326242895115*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01362744780256239  + (0.00000359274216003*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00200183122683721  - (0.00000089207452791*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00833042637239315  + (0.00000946767677294*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00258481175207224  - (0.00000087429351464*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00459744479712244  - (0.00000049519758701*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00534277030993820  + (0.00000397547847155*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00272332919605675  + (0.00000040077229097*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00637243782359372  - (0.00000139419072176*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00233001590327504  + (0.00000420129915747*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00623296727793041  + (0.00000019010664856*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00276177096376805  + (0.00000580301901385*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00559184754866264  + (0.00000080597287792*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00343180144395919  - (0.00000243701142085*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00493325428861701  + (0.00000300985740900*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00396140827680823  - (0.00000051459681789*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00448497879902493  + (0.00000744412841743*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00425146888772076  - (0.00000082346016542*applyconvL)));
					//end Precious (Precision 8) Holo
					break;
				case 6:
					//begin Punch (API) conv
					inputSampleL += (bL[c[1]] * (0.09299870608542582  - (0.00009582362368873*applyconvL)));
					inputSampleL -= (bL[c[2]] * (0.11947847710741009  - (0.00004500891602770*applyconvL)));
					inputSampleL += (bL[c[3]] * (0.09071606264761795  + (0.00005639498984741*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.08561982770836980  - (0.00004964855606916*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.06440549220820363  + (0.00002428052139507*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.05987991812840746  + (0.00000101867082290*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.03980233135839382  + (0.00003312430049041*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03648402630896925  - (0.00002116186381142*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.01826860869525248  + (0.00003115110025396*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.01723968622495364  - (0.00002450634121718*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.00187588812316724  + (0.00002838206198968*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.00381796423957237  - (0.00003155815499462*applyconvL)));
					inputSampleL -= (bL[c[13]] * (0.00852092214496733  - (0.00001702651162392*applyconvL)));
					inputSampleL += (bL[c[14]] * (0.00315560292270588  + (0.00002547861676047*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.01258630914496868  - (0.00004555319243213*applyconvL)));
					inputSampleL += (bL[c[16]] * (0.00536435648963575  + (0.00001812393657101*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01272975658159178  - (0.00004103775306121*applyconvL)));
					inputSampleL += (bL[c[18]] * (0.00403818975172755  + (0.00003764615492871*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.01042617366897483  - (0.00003605210426041*applyconvL)));
					inputSampleL += (bL[c[20]] * (0.00126599583390057  + (0.00004305458668852*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.00747876207688339  - (0.00003731207018977*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00149873689175324  - (0.00005086601800791*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00503221309488033  - (0.00003636086782783*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00342998224655821  - (0.00004103091180506*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00355585977903117  - (0.00003698982145400*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00437201792934817  - (0.00002720235666939*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00299217874451556  - (0.00004446954727956*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00457924652487249  - (0.00003859065778860*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00298182934892027  - (0.00002064710931733*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00438838441540584  - (0.00005223008424866*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00323984218794705  - (0.00003397987535887*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00407693981307314  - (0.00003935772436894*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00350435348467321  - (0.00005525463935338*applyconvL)));
					//end Punch (API) conv
					break;
				case 7:
					//begin Calibre (?) steel
					inputSampleL -= (bL[c[1]] * (0.23505923670562212  - (0.00028312859289245*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.08188436704577637  - (0.00008817721351341*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.05075798481700617  - (0.00018817166632483*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.00455811821873093  + (0.00001922902995296*applyconvL)));
					inputSampleL -= (bL[c[5]] * (0.00027610521433660  - (0.00013252525469291*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.03529246280346626  - (0.00002772989223299*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.01784111585586136  + (0.00010230276997291*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.04394950700298298  - (0.00005910607126944*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.01990770780547606  + (0.00007640328340556*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.04073629569741782  - (0.00007712327117090*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.01349648572795252  + (0.00005959130575917*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.03191590248003717  - (0.00008418000575151*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.00348795527924766  + (0.00005489156318238*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.02198496281481767  - (0.00008471601187581*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00504771152505089  - (0.00005525060587917*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.01391075698598491  - (0.00007929630732607*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.01142762504081717  - (0.00005967036737742*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00893541815021255  - (0.00007535697758141*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.01459704973464936  - (0.00005969199602841*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00694755135226282  - (0.00006930127097865*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.01516695630808575  - (0.00006365800069826*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00705917318113651  - (0.00006497209096539*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.01420501209177591  - (0.00006555654576113*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00815905656808701  - (0.00006105622534761*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.01274326525552961  - (0.00006542652857017*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00937146927845488  - (0.00006051267868722*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.01146573981165209  - (0.00006381511607749*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.01021294359409007  - (0.00005930397856398*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.01065217095323532  - (0.00006371505438319*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.01058751196699751  - (0.00006042857480233*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.01026557827762401  - (0.00006007776163871*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.01060929183604604  - (0.00006114703012726*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.01014533525058528  - (0.00005963567932887*applyconvL)));
					//end Calibre (?)
					break;
				case 8:
					//begin Tube (Manley) conv
					inputSampleL -= (bL[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvL)));
					inputSampleL += (bL[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvL)));
					inputSampleL -= (bL[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvL)));
					inputSampleL -= (bL[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvL)));
					inputSampleL += (bL[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvL)));
					inputSampleL -= (bL[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvL)));
					inputSampleL += (bL[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvL)));
					inputSampleL -= (bL[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvL)));
					inputSampleL += (bL[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvL)));
					inputSampleL -= (bL[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvL)));
					inputSampleL += (bL[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvL)));
					inputSampleL -= (bL[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvL)));
					inputSampleL += (bL[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvL)));
					inputSampleL -= (bL[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvL)));
					inputSampleL -= (bL[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvL)));
					inputSampleL -= (bL[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvL)));
					inputSampleL -= (bL[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvL)));
					inputSampleL -= (bL[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvL)));
					inputSampleL -= (bL[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvL)));
					inputSampleL -= (bL[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvL)));
					inputSampleL -= (bL[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvL)));
					inputSampleL -= (bL[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvL)));
					inputSampleL -= (bL[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvL)));
					inputSampleL -= (bL[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvL)));
					inputSampleL -= (bL[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvL)));
					inputSampleL -= (bL[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvL)));
					inputSampleL -= (bL[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvL)));
					inputSampleL -= (bL[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvL)));
					inputSampleL -= (bL[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvL)));
					inputSampleL -= (bL[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvL)));
					inputSampleL -= (bL[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvL)));
					inputSampleL -= (bL[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvL)));
					inputSampleL -= (bL[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvL)));
					//end Tube (Manley) conv
					break;
			}

			switch (console)
			{
				case 1:
					//begin Cider (Focusrite) MCI
					inputSampleR += (bR[c[1]] * (0.61283288942201319  + (0.00024011410669522*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.24036380659761222  - (0.00020789518206241*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.09104669761717916  + (0.00012829642741548*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.02378290768554025  - (0.00017673646470440*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.02832818490275965  - (0.00013536187747384*applyconvR)));
					inputSampleR += (bR[c[6]] * (0.03268797679215937  + (0.00015035126653359*applyconvR)));
					inputSampleR -= (bR[c[7]] * (0.04024464202655586  - (0.00015034923056735*applyconvR)));
					inputSampleR += (bR[c[8]] * (0.01864890074318696  + (0.00014513281680642*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.01632731954100322  - (0.00015509089075614*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.00318907090555589  - (0.00014784812076550*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.00208573465221869  - (0.00015350520779465*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.00907033901519614  - (0.00015442964157250*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00199458794148013  - (0.00015595640046297*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.00705979153201755  - (0.00015730069418051*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00429488975412722  - (0.00015743697943505*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.00497724878704936  - (0.00016014760011861*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00506059305562353  - (0.00016194824072466*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00483432223285621  - (0.00016329050124225*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00495100420886005  - (0.00016297509798749*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00489319520555115  - (0.00016472839684661*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00489177657970308  - (0.00016791875866630*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00487900894707044  - (0.00016755993898534*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00486234009335561  - (0.00016968157345446*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00485737490288736  - (0.00017180713324431*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00484106070563455  - (0.00017251073661092*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00483219429408410  - (0.00017321683790891*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00482013597437550  - (0.00017392186866488*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00480949628051497  - (0.00017569098775602*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00479992055604049  - (0.00017746046369449*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00478750757986987  - (0.00017745630047554*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00477828651185740  - (0.00017958043287604*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00476906544384494  - (0.00018170456527653*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00475700712413634  - (0.00018099144598088*applyconvR)));
					//end Cider (Focusrite) MCI
					break;
				case 2:
					//begin Rock (SSL) conv
					inputSampleR += (bR[c[1]] * (0.67887916185274055  + (0.00068787552301086*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.25671050678827934  + (0.00017691749454490*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.15135839896615280  + (0.00007481480365043*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.11813512969090802  + (0.00005191138121359*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.08329104347166429  + (0.00001871054659794*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.07663817456103936  + (0.00002751359071705*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.05477586152148759  + (0.00000744843212679*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.05547314737187786  + (0.00001025289931145*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.03822948356540711  - (0.00000249791561457*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.04199383340841713  - (0.00000067328840674*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.02695796542339694  - (0.00000796704606548*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.03228715059431878  - (0.00000579711816722*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.01846929689819187  - (0.00000984017804950*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.02528050435045951  - (0.00000701189792484*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.01207844846859765  - (0.00001522630289356*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01894464378378515  - (0.00001205456372080*applyconvR)));
					inputSampleR += (bR[c[17]] * (0.00667804407593324  - (0.00001343604283817*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.01408418045473130  - (0.00001246443581504*applyconvR)));
					inputSampleR += (bR[c[19]] * (0.00228696509481569  - (0.00001506764046927*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.01006277891348454  - (0.00000970723079112*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00132368373546377  + (0.00001188847238761*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00676615715578373  - (0.00001209129844861*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00426288438418556  + (0.00001286836943559*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00408897698639688  - (0.00001102542567911*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00662040619382751  + (0.00001206328529063*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00196101294183599  - (0.00000950703614981*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00845620581010342  + (0.00001279970295678*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00032595215043616  - (0.00000920518241371*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00982957737435458  + (0.00001177745362317*applyconvR)));
					inputSampleR += (bR[c[30]] * (0.00086920573760513  + (0.00000913758382404*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.01079020871452061  + (0.00000900750153697*applyconvR)));
					inputSampleR += (bR[c[32]] * (0.00167613606334460  + (0.00000732769151038*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.01138050011044332  + (0.00000946908207442*applyconvR)));
					//end Rock (SSL) conv
					break;
				case 3:
					//begin Lush (Neve) conv
					inputSampleR += (bR[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvR)));
					inputSampleR += (bR[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvR)));
					inputSampleR += (bR[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvR)));
					inputSampleR -= (bR[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvR)));
					inputSampleR += (bR[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvR)));
					inputSampleR += (bR[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvR)));
					inputSampleR += (bR[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvR)));
					inputSampleR += (bR[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvR)));
					inputSampleR += (bR[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvR)));
					inputSampleR += (bR[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvR)));
					inputSampleR += (bR[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvR)));
					inputSampleR += (bR[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvR)));
					inputSampleR += (bR[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvR)));
					inputSampleR += (bR[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvR)));
					inputSampleR += (bR[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvR)));
					inputSampleR += (bR[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvR)));
					inputSampleR += (bR[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvR)));
					inputSampleR += (bR[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvR)));
					inputSampleR += (bR[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvR)));
					inputSampleR += (bR[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvR)));
					inputSampleR += (bR[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvR)));
					inputSampleR += (bR[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvR)));
					inputSampleR += (bR[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvR)));
					//end Lush (Neve) conv
					break;
				case 4:
					//begin Elation (LA2A) vibe
					inputSampleR -= (bR[c[1]] * (0.25867935358656502  - (0.00045755657070112*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.11509367290253694  - (0.00017494270657228*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.06709853575891785  - (0.00058913102597723*applyconvR)));
					inputSampleR += (bR[c[4]] * (0.01871006356851681  - (0.00003387358004645*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00794797957360465  - (0.00044224784691203*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.01956921817394220  - (0.00006718936750076*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01682120257195205  + (0.00032857446292230*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03401069039824205  - (0.00013634182872897*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.02369950268232634  + (0.00023112685751657*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.03477071178117132  - (0.00018029792231600*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.02024369717958201  + (0.00017337813374202*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.02819087729102172  - (0.00021438538665420*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.01147946743141303  + (0.00014424066034649*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.01894777011468867  - (0.00021549146262408*applyconvR)));
					inputSampleR += (bR[c[15]] * (0.00301370330346873  + (0.00013527460148394*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01067147835815486  - (0.00020960689910868*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00402715397506384  - (0.00014421582712470*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00502221703392005  - (0.00019805767015024*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00808788533308497  - (0.00016095444141931*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00232696588842683  - (0.00018384470981829*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00943950821324531  - (0.00017098987347593*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00193709517200834  - (0.00018151995939591*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00899713952612659  - (0.00017385835059948*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00280584331659089  - (0.00017742164162470*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00780381001954970  - (0.00018002500755708*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00400370310490333  - (0.00017471691087957*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00661527728186928  - (0.00018137323370347*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00496545526864518  - (0.00017681872601767*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00580728820997532  - (0.00018186220389790*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00549309984725666  - (0.00017722985399075*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00542194777529239  - (0.00018486900185338*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00565992080998939  - (0.00018005824393118*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00532121562846656  - (0.00018643189636216*applyconvR)));
					//end Elation (LA2A)
					break;
				case 5:
					//begin Precious (Precision 8) Holo
					inputSampleR += (bR[c[1]] * (0.59188440274551890  - (0.00008361469668405*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.24439750948076133  + (0.00002651678396848*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.14109876103205621  - (0.00000840487181372*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.10053507128157971  + (0.00001768100964598*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.05859287880626238  - (0.00000361398065989*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.04337406889823660  + (0.00000735941182117*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01589900680531097  + (0.00000207347387987*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.01087234854973281  + (0.00000732123412029*applyconvR)));
					inputSampleR -= (bR[c[9]] * (0.00845782429679176  - (0.00000133058605071*applyconvR)));
					inputSampleR += (bR[c[10]] * (0.00662278586618295  - (0.00000424594730611*applyconvR)));
					inputSampleR -= (bR[c[11]] * (0.02000592193760155  + (0.00000632896879068*applyconvR)));
					inputSampleR += (bR[c[12]] * (0.01321157777167565  - (0.00001421171592570*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.02249955362988238  + (0.00000163937127317*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.01196492077581504  - (0.00000535385220676*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.01905917427000097  + (0.00000121672882030*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00761909482108073  - (0.00000326242895115*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01362744780256239  + (0.00000359274216003*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00200183122683721  - (0.00000089207452791*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00833042637239315  + (0.00000946767677294*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00258481175207224  - (0.00000087429351464*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00459744479712244  - (0.00000049519758701*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00534277030993820  + (0.00000397547847155*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00272332919605675  + (0.00000040077229097*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00637243782359372  - (0.00000139419072176*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00233001590327504  + (0.00000420129915747*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00623296727793041  + (0.00000019010664856*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00276177096376805  + (0.00000580301901385*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00559184754866264  + (0.00000080597287792*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00343180144395919  - (0.00000243701142085*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00493325428861701  + (0.00000300985740900*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00396140827680823  - (0.00000051459681789*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00448497879902493  + (0.00000744412841743*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00425146888772076  - (0.00000082346016542*applyconvR)));
					//end Precious (Precision 8) Holo
					break;
				case 6:
					//begin Punch (API) conv
					inputSampleR += (bR[c[1]] * (0.09299870608542582  - (0.00009582362368873*applyconvR)));
					inputSampleR -= (bR[c[2]] * (0.11947847710741009  - (0.00004500891602770*applyconvR)));
					inputSampleR += (bR[c[3]] * (0.09071606264761795  + (0.00005639498984741*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.08561982770836980  - (0.00004964855606916*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.06440549220820363  + (0.00002428052139507*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.05987991812840746  + (0.00000101867082290*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.03980233135839382  + (0.00003312430049041*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03648402630896925  - (0.00002116186381142*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.01826860869525248  + (0.00003115110025396*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.01723968622495364  - (0.00002450634121718*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.00187588812316724  + (0.00002838206198968*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.00381796423957237  - (0.00003155815499462*applyconvR)));
					inputSampleR -= (bR[c[13]] * (0.00852092214496733  - (0.00001702651162392*applyconvR)));
					inputSampleR += (bR[c[14]] * (0.00315560292270588  + (0.00002547861676047*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.01258630914496868  - (0.00004555319243213*applyconvR)));
					inputSampleR += (bR[c[16]] * (0.00536435648963575  + (0.00001812393657101*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01272975658159178  - (0.00004103775306121*applyconvR)));
					inputSampleR += (bR[c[18]] * (0.00403818975172755  + (0.00003764615492871*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.01042617366897483  - (0.00003605210426041*applyconvR)));
					inputSampleR += (bR[c[20]] * (0.00126599583390057  + (0.00004305458668852*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.00747876207688339  - (0.00003731207018977*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00149873689175324  - (0.00005086601800791*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00503221309488033  - (0.00003636086782783*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00342998224655821  - (0.00004103091180506*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00355585977903117  - (0.00003698982145400*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00437201792934817  - (0.00002720235666939*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00299217874451556  - (0.00004446954727956*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00457924652487249  - (0.00003859065778860*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00298182934892027  - (0.00002064710931733*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00438838441540584  - (0.00005223008424866*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00323984218794705  - (0.00003397987535887*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00407693981307314  - (0.00003935772436894*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00350435348467321  - (0.00005525463935338*applyconvR)));
					//end Punch (API) conv
					break;
				case 7:
					//begin Calibre (?) steel
					inputSampleR -= (bR[c[1]] * (0.23505923670562212  - (0.00028312859289245*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.08188436704577637  - (0.00008817721351341*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.05075798481700617  - (0.00018817166632483*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.00455811821873093  + (0.00001922902995296*applyconvR)));
					inputSampleR -= (bR[c[5]] * (0.00027610521433660  - (0.00013252525469291*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.03529246280346626  - (0.00002772989223299*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.01784111585586136  + (0.00010230276997291*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.04394950700298298  - (0.00005910607126944*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.01990770780547606  + (0.00007640328340556*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.04073629569741782  - (0.00007712327117090*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.01349648572795252  + (0.00005959130575917*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.03191590248003717  - (0.00008418000575151*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.00348795527924766  + (0.00005489156318238*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.02198496281481767  - (0.00008471601187581*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00504771152505089  - (0.00005525060587917*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.01391075698598491  - (0.00007929630732607*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.01142762504081717  - (0.00005967036737742*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00893541815021255  - (0.00007535697758141*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.01459704973464936  - (0.00005969199602841*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00694755135226282  - (0.00006930127097865*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.01516695630808575  - (0.00006365800069826*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00705917318113651  - (0.00006497209096539*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.01420501209177591  - (0.00006555654576113*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00815905656808701  - (0.00006105622534761*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.01274326525552961  - (0.00006542652857017*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00937146927845488  - (0.00006051267868722*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.01146573981165209  - (0.00006381511607749*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.01021294359409007  - (0.00005930397856398*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.01065217095323532  - (0.00006371505438319*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.01058751196699751  - (0.00006042857480233*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.01026557827762401  - (0.00006007776163871*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.01060929183604604  - (0.00006114703012726*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.01014533525058528  - (0.00005963567932887*applyconvR)));
					//end Calibre (?)
					break;
				case 8:
					//begin Tube (Manley) conv
					inputSampleR -= (bR[c[1]] * (0.20641602693167951  - (0.00078952185394898*applyconvR)));
					inputSampleR += (bR[c[2]] * (0.07601816702459827  + (0.00022786334179951*applyconvR)));
					inputSampleR -= (bR[c[3]] * (0.03929765560019285  - (0.00054517993246352*applyconvR)));
					inputSampleR -= (bR[c[4]] * (0.00298333157711103  - (0.00033083756545638*applyconvR)));
					inputSampleR += (bR[c[5]] * (0.00724006282304610  + (0.00045483683460812*applyconvR)));
					inputSampleR -= (bR[c[6]] * (0.03073108963506036  - (0.00038190060537423*applyconvR)));
					inputSampleR += (bR[c[7]] * (0.02332434692533051  + (0.00040347288688932*applyconvR)));
					inputSampleR -= (bR[c[8]] * (0.03792606869061214  - (0.00039673687335892*applyconvR)));
					inputSampleR += (bR[c[9]] * (0.02437059376675688  + (0.00037221210539535*applyconvR)));
					inputSampleR -= (bR[c[10]] * (0.03416764311979521  - (0.00040314850796953*applyconvR)));
					inputSampleR += (bR[c[11]] * (0.01761669868102127  + (0.00035989484330131*applyconvR)));
					inputSampleR -= (bR[c[12]] * (0.02538237753523052  - (0.00040149119125394*applyconvR)));
					inputSampleR += (bR[c[13]] * (0.00770737340728377  + (0.00035462118723555*applyconvR)));
					inputSampleR -= (bR[c[14]] * (0.01580706228482803  - (0.00037563141307594*applyconvR)));
					inputSampleR -= (bR[c[15]] * (0.00055119240005586  - (0.00035409299268971*applyconvR)));
					inputSampleR -= (bR[c[16]] * (0.00818552143438768  - (0.00036507661042180*applyconvR)));
					inputSampleR -= (bR[c[17]] * (0.00661842703548304  - (0.00034550528559056*applyconvR)));
					inputSampleR -= (bR[c[18]] * (0.00362447476272098  - (0.00035553012761240*applyconvR)));
					inputSampleR -= (bR[c[19]] * (0.00957098027225745  - (0.00034091691045338*applyconvR)));
					inputSampleR -= (bR[c[20]] * (0.00193621774016660  - (0.00034554529131668*applyconvR)));
					inputSampleR -= (bR[c[21]] * (0.01005433027357935  - (0.00033878223153845*applyconvR)));
					inputSampleR -= (bR[c[22]] * (0.00221712428802004  - (0.00033481410137711*applyconvR)));
					inputSampleR -= (bR[c[23]] * (0.00911255639207995  - (0.00033263425232666*applyconvR)));
					inputSampleR -= (bR[c[24]] * (0.00339667169034909  - (0.00032634428038430*applyconvR)));
					inputSampleR -= (bR[c[25]] * (0.00774096948249924  - (0.00032599868802996*applyconvR)));
					inputSampleR -= (bR[c[26]] * (0.00463907626773794  - (0.00032131993173361*applyconvR)));
					inputSampleR -= (bR[c[27]] * (0.00658222997260378  - (0.00032014977430211*applyconvR)));
					inputSampleR -= (bR[c[28]] * (0.00550347079924993  - (0.00031557153256653*applyconvR)));
					inputSampleR -= (bR[c[29]] * (0.00588754981375325  - (0.00032041307242303*applyconvR)));
					inputSampleR -= (bR[c[30]] * (0.00590293898419892  - (0.00030457857428714*applyconvR)));
					inputSampleR -= (bR[c[31]] * (0.00558952010441800  - (0.00030448053548086*applyconvR)));
					inputSampleR -= (bR[c[32]] * (0.00598183557634295  - (0.00030715064323181*applyconvR)));
					inputSampleR -= (bR[c[33]] * (0.00555223929714115  - (0.00030319367948553*applyconvR)));
					//end Tube (Manley) conv
					break;
			}

			bc_bridgerectifier = fabs(inputSampleL);
			bc_bridgerectifier = 1.0-cos(bc_bridgerectifier);
			if (inputSampleL > 0) inputSampleL -= bc_bridgerectifier;
			else inputSampleL += bc_bridgerectifier;

			bc_bridgerectifier = fabs(inputSampleR);
			bc_bridgerectifier = 1.0-cos(bc_bridgerectifier);
			if (inputSampleR > 0) inputSampleR -= bc_bridgerectifier;
			else inputSampleR += bc_bridgerectifier;


			if (bc_outgain != 1.0) {
				inputSampleL *= bc_outgain;
				inputSampleR *= bc_outgain;
			}

			if (bc_wet !=1.0) {
				inputSampleL = (inputSampleL * bc_wet) + (bc_drySampleL * (1.0-bc_wet));
				inputSampleR = (inputSampleR * bc_wet) + (bc_drySampleR * (1.0-bc_wet));
			}
		}

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

		in1++;
		in2++;
		out1++;
		out2++;
    }
}
