/* ========================================
 *  ToneSlantV - ToneSlantV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ToneSlantV_H
#include "ToneSlantV.h"
#endif

void ToneSlantV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];


	double overallscale = (A*99.0)+1.0;
	double applySlant = (B*2.0)-1.0;
	
	
	//count to f(gain) which will be 0. f(0) is x1
	for (int count = 0; count <= overallscale; count+=4) {
		Vec4d os(overallscale);
		((Vec4d(1.0) - (Vec4d(count, count + 1, count + 2, count + 3) / os)) / os).store_a(ts_f + count);
		//recalc the filter and don't change the buffer it'll apply to
	}
		
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d fpd_d = to_double(fpd);

  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, fpd_d * 1.18e-17, inputSample); 

		inputSample.store_a(ts_b + gcount);
		if(gcount == 0) inputSample.store_a(ts_b + 256); 

		Vec4d acc = inputSample * ts_f[0];
		Vec4d b;
		for (int count = 1; count < overallscale; count++) {
			b.load(ts_b + ((gcount - (count * 2)) & 255));
			acc += b * ts_f[count];
		}
		//we're gonna apply the total effect of all these calculations as a single subtract
		inputSample += ((inputSample - (acc*2.0)) * applySlant);
		//our one math operation on the input data coming in

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
		gcount = (gcount + 4) & 255;
	}

	fpd.store_a(fpd_b);
}

void ToneSlantV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];


}
