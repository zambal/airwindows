/* ========================================
 *  SubTightV - SubTightV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __SubTightV_H
#include "SubTightV.h"
#endif

void SubTightV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();
	int subStages = pow(B,2)*16.0;
	if (subStages < 1) subStages = 1;
	double subTrim = pow((A*0.3)+(pow(B,2)*0.2),subStages)/overallscale;
	//to use this as an analog modeler for restricting digital lows, find set values that still show bass
	//Note that this is best used sparingly, on the 'not enough subtraction' side of the node.

	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d fpd_d = to_double(fpd);

  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, fpd_d * 1.18e-17, inputSample); 
		
		//you want subStages and subTrim to be hardcoded values when embedding this into something else
		//then it only needs the sub[] array, and to have it initialized to 0.0
		
		//begin SubTightV section
		Vec4d subSample = inputSample * subTrim;
		Vec2d subSample1 = subSample.get_low(), subSample2 = subSample.get_high();

		for (int x = 0; x < subStages; x++) {
			Vec2d scale = 0.5+abs(subSample1*0.5);
			Vec2d s; s.load_a(sub + (x * 2));
			subSample1 = (s+(sin(s-subSample1)*scale));
			s = subSample1*scale;

			scale = 0.5+abs(subSample2*0.5);
			subSample2 = (s+(sin(s-subSample2)*scale));
			(subSample2*scale).store_a(sub + (x * 2));
		}

		subSample = concatenate2(subSample1, subSample2);

		if (subStages % 2 > 0) {
			subSample = -subSample;
		}

		inputSample -= (min(max(subSample, -0.25), 0.25)*16.0);
		//end SubTightV section
		
		//cut the level WAY down, then the modified Creature code blows up subs.
		//the adjustment of scale destabilizes the routine to blow up more DC.
		//this is boosted by 24dB and subtracted from the dry signal
		
		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		fpd_d = to_double(fpd);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (fpd_d - 2147483647.0) * to_double(exp) * 5.5e-36l;
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

void SubTightV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];
	
}
