/* ========================================
 *  ButterComp2V - ButterComp2V.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ButterComp2V_H
#include "ButterComp2V.h"
#endif

void ButterComp2V::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	double inputgain = pow(10.0,(A*14.0)/20.0);
	double compfactor = 0.012 * (A / 135.0);
	double output = B * 2.0;
	double wet = C;
	//removed extra dry variable
	double outputgain = inputgain;
	outputgain -= 1.0;
	outputgain /= 1.5;
	outputgain += 1.0;

	Vec4ui fpd; fpd.load(fpd_b);
	Vec4d controlpos; controlpos.load(controlpos_b);
	Vec4d controlneg; controlneg.load(controlneg_b);
	Vec4d targetpos; targetpos.load(targetpos_b);
	Vec4d targetneg; targetneg.load(targetneg_b);
	Vec4d lastOutput; lastOutput.load(lastOutput_b);
	
    
    while (sampleFrames >= 0)
    {
		Vec4d	totalmultiplier;
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		{
			Vec4d divisor = Vec4d(compfactor) / (abs(lastOutput)+1.0);
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

			totalmultiplier = (controlpos * outputpos) + (controlneg * outputneg);
			//this combines the sides according to flip, blending relative to the input value
		}
		
		Vec4d drySample = inputSample;
		inputSample *= inputgain;
		inputSample /= outputgain;
		inputSample *= totalmultiplier;
		
	
		if (output != 1.0) {
			inputSample *= output;
		}

		if (wet !=1.0) {
			inputSample = (inputSample * wet) + (drySample * (1.0-wet));
		}
		
		lastOutput = inputSample;
		//we will make this factor respond to use of dry/wet
		
		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (to_double(fpd) - 2147483647.0) * to_double(exp) * 5.5e-36l;
		//end 32 bit stereo floating point dither

		double result[4];
		inputSample.store(result);
		*out1 = result[0];
		*out2 = result[1];
		*(out1 + 1) = result[2];
		*(out2 + 1) = result[3];

		in1+= 2;
		in2+= 2;
		out1+= 2;
		out2+= 2;

		sampleFrames -= 2;
	}

	fpd.store(fpd_b);
	controlpos.store(controlpos_b);
	controlneg.store(controlneg_b);
	targetpos.store(targetpos_b);
	targetneg.store(targetneg_b);
	lastOutput.store(lastOutput_b);
}

void ButterComp2V::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
