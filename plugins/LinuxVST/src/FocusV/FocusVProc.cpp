/* ========================================
 *  FocusV - FocusV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __FocusV_H
#include "FocusV.h"
#endif

void FocusV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];
	
	//[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
	//[1] is resonance, 0.7071 is Butterworth. Also can't be zero
	double boost = pow(10.0,(A*12.0)/20.0);
	double cutoff = 3515.775/getSampleRate(); //fixed frequency, 3.515775k
	double res = pow(pow(B,3)*2,2)+0.0001; //resonance
	int mode = (int) ( C * 4.999 );
	double output = D;
	double wet = E;
	
	double K = tan(M_PI * cutoff);
	double norm = 1.0 / (1.0 + K / res + K * K);
	double f1 = K / res * norm;
	double f2 = -f1;
	double f3 = 2.0 * (K * K - 1.0) * norm;
	double f4 = (1.0 - K / res + K * K) * norm;
		
	Vec4ui fpd; fpd.load_a(fpd_b);

  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		Vec4d drySample = inputSample;
		
		inputSample = sin(inputSample);
		//encode Console5: good cleanness
		{
			Vec2d p1, p2; p1.load_a(figure); p2.load_a(figure + 2);
			Vec2d tempSample = (inputSample.get_low() * f1) + p1;
			p1 = -(tempSample * f3) + p2;
			p2 = (inputSample.get_low() * f2) - (tempSample * f4);
			Vec2d inputSampleLow = tempSample;
		
			tempSample = (inputSample.get_high() * f1) + p1;
			p1 = -(tempSample * f3) + p2;
			p2 = (inputSample.get_high() * f2) - (tempSample * f4);
			p1.store_a(figure); p2.store_a(figure + 2);
			inputSample = concatenate2(inputSampleLow, tempSample);
		}		
		inputSample = min(max(inputSample, -1.0), 1.0);
		//without this, you can get a NaN condition where it spits out DC offset at full blast!

		inputSample = Vec4d(asin(inputSample.extract(0)), asin(inputSample.extract(1)), asin(inputSample.extract(2)), asin(inputSample.extract(3)));
		//decode Console5
		
		Vec4d groundSample = drySample - inputSample; //set up UnBox
		inputSample *= boost; //now, focussed area gets cranked before distort
		
		switch (mode)
		{
			case 0: //Density
				inputSample = min(max(inputSample, -1.570796326794897), 1.570796326794897);
				//clip to 1.570796326794897 to reach maximum output
				inputSample = sin(inputSample);
				break;
			case 1: //Drive				
				inputSample -= (inputSample * (abs(inputSample) * 0.6) * (abs(inputSample) * 0.6));
				inputSample *= 1.6;
				break;
			case 2: //Spiral
				inputSample = min(max(inputSample, -1.2533141373155), 1.2533141373155);
				//clip to 1.2533141373155 to reach maximum output
				inputSample = sin(inputSample * abs(inputSample)) / select(abs(inputSample) == 0.0, Vec4d(1.0), abs(inputSample));
				break;
			case 3: //Mojo
				Vec4d mojo; mojo = pow(abs(inputSample),Vec4d(0.25));
				inputSample = select(mojo > 0.0, (sin(inputSample * mojo * M_PI * 0.5) / mojo) * 0.987654321, inputSample);
				//mojo is the one that flattens WAAAAY out very softly before wavefolding				
				break;
			case 4: //Dyno
				Vec4d dyno; dyno = pow_const(abs(inputSample),4);
				inputSample = select(dyno > 0.0, (sin(inputSample * dyno) / dyno) * 1.1654321, inputSample);
				//dyno is the one that tries to raise peak energy				
				break;
		}				
		
		if (output != 1.0) {
			inputSample *= output;
		}
		
		inputSample += groundSample; //effectively UnBox
		
		if (wet !=1.0) {
			inputSample = (inputSample * wet) + (drySample * (1.0-wet));
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
}

void FocusV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];
	
}
