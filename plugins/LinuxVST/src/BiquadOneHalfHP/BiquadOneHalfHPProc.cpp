/* ========================================
 *  BiquadOneHalfHP - BiquadOneHalfHP.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __BiquadOneHalfHP_H
#include "BiquadOneHalfHP.h"
#endif


void BiquadOneHalfHP::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];


		double b0;
		double b1;
		double b2;
		double b3;
		double b4;

		int type = ceil((A*3.999)+0.00001);

		double cutoff = ((B*B*B*0.9999)+0.0001)*0.499;
		if (cutoff < 0.0001) cutoff = 0.0001;

	  double res = (C*C*C*29.99)+0.01;
		if (res < 0.0001) res = 0.0001;

		//biquad contains these values:
		//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
		//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
		//[2] is a0 but you need distinct ones for additional biquad instances so it's here
		//[3] is a1 but you need distinct ones for additional biquad instances so it's here
		//[4] is a2 but you need distinct ones for additional biquad instances so it's here
		//[5] is b1 but you need distinct ones for additional biquad instances so it's here
		//[6] is b2 but you need distinct ones for additional biquad instances so it's here
		//[7] is stored delayed sample (freq and res are stored so you can move them sample by sample)
		//[8] is stored delayed sample (you have to include the coefficient making code if you do that)

		//to build a dedicated filter, rename 'biquad' to whatever the new filter is, then
		//put this code either within the sample buffer (for smoothly modulating freq or res)
		//or in this 'read the controls' area (for letting you change freq and res with controls)
		//or in 'reset' if the freq and res are absolutely fixed (use GetSampleRate to define freq)

		if (type == 1) { //lowpass
			double K = tan(M_PI * cutoff);
			double norm = 1.0 / (1.0 + K / res + K * K);
			b0 = K * K * norm;
			b1 = 2.0 * b0;
			b2 = b0;
			b3 = 2.0 * (K * K - 1.0) * norm;
			b4 = (1.0 - K / res + K * K) * norm;
		}

		if (type == 2) { //highpass
			double K = tan(M_PI * cutoff);
			double norm = 1.0 / (1.0 + K / res + K * K);
			b0 = norm;
			b1 = -2.0 * b0;
			b2 = b0;
			b3 = 2.0 * (K * K - 1.0) * norm;
			b4 = (1.0 - K / res + K * K) * norm;
		}

		if (type == 3) { //bandpass
			double K = tan(M_PI * cutoff);
			double norm = 1.0 / (1.0 + K / res + K * K);
			b0 = K / res * norm;
			b1 = 0.0; //bandpass can simplify the biquad kernel: leave out this multiply
			b2 = -b0;
			b3 = 2.0 * (K * K - 1.0) * norm;
			b4 = (1.0 - K / res + K * K) * norm;
		}

		if (type == 4) { //notch
			double K = tan(M_PI * cutoff);
			double norm = 1.0 / (1.0 + K / res + K * K);
			b0 = (1.0 + K * K) * norm;
			b1 = 2.0 * (K * K - 1) * norm;
			b2 = b0;
			b3 = b1;
			b4 = (1.0 - K / res + K * K) * norm;
		}

		Vec4d wet = Vec4d((D*2.0)-1.0);

		Vec4d biquadA; biquadA.load(biquadPrev);
		Vec4d biquadB; biquadB.load(biquadPrev + 4);
		Vec4ui fpd; fpd.load(fpd_b);

    while (sampleFrames > 0)
    {
			Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
			inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
			Vec4d drySample = inputSample;
		
			inputSample += ((pow_const(inputSample,5)/128.0) + (pow_const(inputSample,9)/262144.0)) - ((pow_const(inputSample,3)/8.0) + (pow_const(inputSample,7)/4096.0));
			//encode Console5: good cleanness

			Vec4d tempSample = inputSample * b0 + biquadA;
			biquadA = (inputSample * b1) - (tempSample * b3) + biquadB;
			biquadB = (inputSample * b2) - (tempSample * b4);
			inputSample = min(max(tempSample, -1.0), 1.0);
			//without this, you can get a NaN condition where it spits out DC offset at full blast!

			inputSample += (pow_const(inputSample,3)/4.0)+(pow_const(inputSample,5)/8.0)+(pow_const(inputSample,7)/16.0)+(pow_const(inputSample,9)/32.0);
			//amplitude aspect

			inputSample = (inputSample*wet) + (drySample*(Vec4d(1.0)-abs(wet)));
			//inv/dry/wet lets us turn [i]P into HP and band into notch
	
		
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

		biquadA.store(biquadPrev);
		biquadB.store(biquadPrev + 4);
		fpd.store(fpd_b);
}

void BiquadOneHalfHP::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];
	

}
