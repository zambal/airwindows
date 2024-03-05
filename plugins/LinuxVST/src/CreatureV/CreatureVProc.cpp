/* ========================================
 *  CreatureV - CreatureV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __CreatureV_H
#include "CreatureV.h"
#endif

void CreatureV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();
	
	double source = 1.0-pow(1.0-A,5);
	int stages = (pow(B,2)*32.0*sqrt(overallscale))+1;
	double wet = (C*2.0)-1.0; //inv-dry-wet for highpass
	double dry = 2.0-(C*2.0);
	if (dry > 1.0) dry = 1.0; //full dry for use with inv, to 0.0 at full wet
	
	Vec4ui fpd; fpd.load_a(fpd_b);
	while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 

		Vec4d drySample = inputSample;
		
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
		
		inputSample = (drySample * dry) + (inputSample * wet);
		
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
}

void CreatureV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
