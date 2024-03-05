/* ========================================
 *  Channel9V - Channel9V.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __Channel9V_H
#include "Channel9V.h"
#endif

void Channel9V::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();	
	double localiirAmount = iirAmount / overallscale;
	double localthreshold = threshold; //we've learned not to try and adjust threshold for sample rate
	double density = B*2.0; //0-2
	double phattity = density - 1.0;
	if (density > 1.0) density = 1.0; //max out at full wet for Spiral aspect
	if (phattity < 0.0) phattity = 0.0; //
	double nonLin = 5.0-density; //number is smaller for more intense, larger for more subtle
	double c = cutoff / getSampleRate();
  double resA = 1.618033988749894848204586;
	double resB = 0.618033988749894848204586;
	
	double K = tan(M_PI * c);
	double norm = 1.0 / (1.0 + K / resA + K * K);
	double b0A = K * K * norm;
	double b1A = 2.0 * b0A;
	double b2A = b0A;
	double b3A = 2.0 * (K * K - 1.0) * norm;
	double b4A = (1.0 - K / resA + K * K) * norm;
		
	K = tan(M_PI * c);
	norm = 1.0 / (1.0 + K / resB + K * K);
	double b0B = K * K * norm;
	double b1B = 2.0 * b0B;
	double b2B = b0B;
	double b3B = 2.0 * (K * K - 1.0) * norm;
	double b4B = (1.0 - K / resB + K * K) * norm;
		
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec2d biquad1A; biquad1A.load_a(biquad);
	Vec2d biquad2A; biquad2A.load_a(biquad + 2);	
	Vec2d biquad3A; biquad3A.load_a(biquad + 4);
	Vec2d biquad4A; biquad4A.load_a(biquad + 6);	
	Vec2d biquad1B; biquad1B.load_a(biquad + 8);	
	Vec2d biquad2B; biquad2B.load_a(biquad + 10);	
	Vec2d biquad3B; biquad3B.load_a(biquad + 12);	
	Vec2d biquad4B; biquad4B.load_a(biquad + 14);	
	Vec2d lastSampleA; lastSampleA.load_a(lastSample);
	Vec2d lastSampleB; lastSampleB.load_a(lastSample + 2);
	Vec2d lastSampleC; lastSampleC.load_a(lastSample + 4);
	Vec4d iir; iir.load_a(iirSamples);

  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		
		if (cutoff < 0.49999) {
			Vec2d tempSample = b0A*inputSample.get_low()+b1A*biquad1A+b2A*biquad2A-b3A*biquad3A-b4A*biquad4A;
  		biquad2A = biquad1A; biquad1A = inputSample.get_low(); 
			Vec2d inputSampleLow = select(abs(tempSample)<1.18e-37, 0.0, tempSample);
			biquad4A = biquad3A; biquad3A = inputSampleLow; //DF1 

			tempSample = b0A*inputSample.get_high()+b1A*biquad1A+b2A*biquad2A-b3A*biquad3A-b4A*biquad4A;
  		biquad2A = biquad1A; biquad1A = inputSample.get_high(); 
			inputSample = concatenate2(inputSampleLow, select(abs(tempSample)<1.18e-37, 0.0, tempSample));
			biquad4A = biquad3A; biquad3A = inputSample.get_high(); //DF1 
		}		
		
		{
			Vec4d dielectricScale = abs(2.0-((inputSample+nonLin)/nonLin));
		

			iir = select(abs(iir)<1.18e-37, 0.0, iir);
			iir = (iir * (1.0 - (localiirAmount * dielectricScale))) + (inputSample * localiirAmount * dielectricScale);
			inputSample -= iir;
			//highpass section
		}
		{
			Vec4d drySample = inputSample;

			inputSample = min(max(inputSample, -1.0), 1.0);
		
			Vec4d phatSample = sin(inputSample * 1.57079633);
			inputSample *= 1.2533141373155;
			//clip to 1.2533141373155 to reach maximum output, or 1.57079633 for pure sine 'phat' version
		
			Vec4d distSample = sin(inputSample * abs(inputSample)) / select(abs(inputSample) == 0.0, Vec4d(1.0), abs(inputSample));
		
			inputSample = distSample; //purest form is full Spiral
			if (density < 1.0) inputSample = (drySample*(1-density))+(distSample*density); //fade Spiral aspect
			if (phattity > 0.0) inputSample = (inputSample*(1-phattity))+(phatSample*phattity); //apply original Density on top
		}
		{
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
		}
		if (C < 1.0) {
			inputSample *= C;
		}
		
		if (cutoff < 0.49999) {
			Vec2d tempSample = b0B*inputSample.get_low()+b1B*biquad1B+b2B*biquad2B-b3B*biquad3B-b4B*biquad4B;
  		biquad2B = biquad1B; biquad1B = inputSample.get_low(); 
			Vec2d inputSampleLow = select(abs(tempSample)<1.18e-37, 0.0, tempSample);
			biquad4B = biquad3B; biquad3B = inputSampleLow; //DF1 

			tempSample = b0B*inputSample.get_high()+b1B*biquad1B+b2B*biquad2B-b3B*biquad3B-b4B*biquad4B;
  		biquad2B = biquad1B; biquad1B = inputSample.get_high(); 
			inputSample = concatenate2(inputSampleLow, select(abs(tempSample)<1.18e-37, 0.0, tempSample));
			biquad4B = biquad3B; biquad3B = inputSample.get_high(); //DF1 
		}		
		
		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (to_double(fpd) - 2147483647.0) * to_double(exp) * 5.5e-36l;
		//end 32 bit stereo floating point dither

		double result[4];
		inputSample.store_a(result);
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

	fpd.store_a(fpd_b);
	biquad1A.store_a(biquad);
	biquad2A.store_a(biquad + 2);	
	biquad3A.store_a(biquad + 4);
	biquad4A.store_a(biquad + 6);	
	biquad1B.store_a(biquad + 8);	
	biquad2B.store_a(biquad + 10);	
	biquad3B.store_a(biquad + 12);	
	biquad4B.store_a(biquad + 14);	
	lastSampleA.store_a(lastSample);
	lastSampleB.store_a(lastSample + 2);
	lastSampleC.store_a(lastSample + 4);
	iir.store_a(iirSamples);
}

void Channel9V::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
