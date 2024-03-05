/* ========================================
 *  CSCompV - CSCompV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __CSCompV_H
#include "CSCompV.h"
#endif

void CSCompV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];
	
	double overallscale = 1.0;
	overallscale /= 44100.0;
	double compscale = overallscale;
	overallscale = getSampleRate();
	compscale = compscale * overallscale;
	//compscale is the one that's 1 or something like 2.2 for 96K rates
	double fpOld = 0.618033988749894848204586; //golden ratio!
	double fpNew = 1.0 - fpOld;
	
	//begin ButterComp
	double inputgain = (pow(A,4)*35)+1.0;
	double compoutgain = inputgain;
	compoutgain -= 1.0;
	compoutgain /= 1.2;
	compoutgain += 1.0;
	double divisor = (0.008 * pow(B,2))+0.0004;
	//originally 0.012
	divisor /= compscale;
	double remainder = divisor;
	divisor = 1.0 - divisor;
	bool engageComp = false;
	if (inputgain > 1.0) engageComp = true;
	//end ButterComp

	double outputgain = C*3.0; //0-2
	
	Vec4ui fpd; fpd.load_a(fpd_b);
	while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		
		//begin compressor
		if (engageComp)
		{
			inputSample *= inputgain;
			
			Vec2d t; t.load_a(cs_avg); inputSample.get_high().store_a(cs_avg);
			Vec4d avg = concatenate2(t, inputSample.get_low());
			Vec4d input = max((inputSample * fpOld) + (avg * fpNew) + 1.0, 0.0);
			Vec4d outputpos = min(input * 0.5, 1.0);
			Vec4d target; target.load_a(cs_targetpos);
			input *= input;
			target *= divisor;
			target += (input * remainder);
			Vec4d calcpos = pow_const((1.0/target),2);
			target.store_a(cs_targetpos);

			input = max((-inputSample * fpOld) + (-avg * fpNew) + 1.0, 0.0);
			Vec4d outputneg = min(input * 0.5, 1.0);
			target.load_a(cs_targetneg);
			input *= input;
			target *= divisor;
			target += (input * remainder);
			Vec4d calcneg = pow_const((1.0/target),2);
			target.store_a(cs_targetneg);
			
			Vec4d controlpos; controlpos.load_a(cs_controlpos);
			Vec4d controlneg; controlneg.load_a(cs_controlneg);
			
			controlpos = select(inputSample > 0.0, (controlpos * divisor) + (calcpos * remainder), controlpos);
			controlneg = select(inputSample < 0.0, (controlneg * divisor) + (calcneg * remainder), controlneg);

			controlpos.store_a(cs_controlpos);
			controlneg.store_a(cs_controlneg);

			Vec4d totalmultiplier = (controlpos * outputpos) + (controlneg * outputneg);
			//this combines the sides according to flip, blending relative to the input value
			
			inputSample *= totalmultiplier;
			inputSample /= compoutgain;

			//built in output trim and dry/wet if desired
			if (outputgain != 1.0) {
				inputSample *= outputgain;
			}
		}
		//end compressor

		
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

void CSCompV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];
	
}
