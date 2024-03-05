/* ========================================
 *  ResEQ2V - ResEQ2V.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ResEQ2V_H
#include "ResEQ2V.h"
#endif

void ResEQ2V::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();
	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;
	//this is going to be 2 for 88.1 or 96k, 3 for silly people, 4 for 176 or 192k
	
	//begin ResEQ2V Mid Boost
	double freqMPeak = pow(A+0.15,3);
	double amountMPeak = pow(B,2);
	int maxMPeak = (amountMPeak*63.0)+1;
	// if ((freqMPeak != prevfreqMPeak)||(amountMPeak != prevamountMPeak)) {
		for (int x = 0; x < maxMPeak; x++) {
			if (((double)x*freqMPeak) < M_PI_4) f[x] = sin(((double)x*freqMPeak)*4.0)*freqMPeak*sin(((double)(maxMPeak-x)/(double)maxMPeak)*M_PI_2);
			else f[x] = cos((double)x*freqMPeak)*freqMPeak*sin(((double)(maxMPeak-x)/(double)maxMPeak)*M_PI_2);
		}
		// prevfreqMPeak = freqMPeak; prevamountMPeak = amountMPeak;
	// }//end ResEQ2V Mid Boost
    
	Vec4ui fpd; fpd.load_a(fpd_b);

  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		
		//begin ResEQ2V Mid Boost
		inputSample.store_a(mpk + mpc);
		if(mpc == 0) inputSample.store_a(mpk + 512);
		Vec4d midMPeak(0.0), m,k;
		
		for (int x = 0; x < maxMPeak; x++) {
			int y = x*cycleEnd;
			m = 0;
			for(int yy = y; yy < (y + cycleEnd); yy++) {
				int offset = (mpc - (yy*2)) & 511;
				k.load(mpk + offset);
				m += k;
			} 
			switch (cycleEnd)
			{
				case 1: m *= f[x];	break;
				case 2: m *= (f[x] * 0.5);	break;
				case 3: m *= (f[x] * 0.333);	break;
				case 4: m *= (f[x] * 0.25); break;
			}
			midMPeak += m;
		}
		mpc = (mpc + 4) & 511;
		
		inputSample = (midMPeak*amountMPeak)+((1.5-amountMPeak>1.0)?inputSample:inputSample*(1.5-amountMPeak));
		//end ResEQ2V Mid Boost

	
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

void ResEQ2V::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
