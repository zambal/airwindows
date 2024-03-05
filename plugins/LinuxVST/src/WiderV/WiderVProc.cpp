/* ========================================
 *  WiderV - WiderV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __WiderV_H
#include "WiderV.h"
#endif

void WiderV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();

	double densityside = (A*2.0)-1.0;
	double densitymid = (B*2.0)-1.0;
	double wet = C;
	//removed extra dry variable
	wet *= 0.5; //we make mid-side by adding/subtracting both channels into each channel
	//and that's why we gotta divide it by 2: otherwise everything's doubled. So, premultiply it to save an extra 'math'
	double offset = (densityside-densitymid)/2;
	if (offset > 0) offset = sin(offset);
	if (offset < 0) offset = -sin(-offset);
	offset = -(pow(offset,4) * 20 * overallscale);
	int near = (int)floor(fabs(offset));
	double farLevel = fabs(offset) - near;
	int far = near + 1;
	double nearLevel = 1.0 - farLevel;
	double bridgerectifier;
	//interpolating the sample
    
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d density(densitymid,densitymid,densityside,densityside);
  Vec4d out = abs(density);
	while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 

		Vec4d drySample = inputSample;

		//assign working variables
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
		
		inputSample = (drySample * (1.0-wet)) + ((permute4<0, 0, 1, 1>(ms) + (permute4<2, 2, 3, 3>(ms) * Vec4d(1, -1, 1, -1))) * wet * 0.5);
		
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

void WiderV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
