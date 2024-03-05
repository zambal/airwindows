/* ========================================
 *  Desk4V - Desk4V.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __Desk4V_H
#include "Desk4V.h"
#endif

void Desk4V::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();
	
	double gain = (pow(A,2)*10)+0.0001;
	double gaintrim = (pow(A,2)*2)+1.0;
	double slewgain = (pow(B,3)*40)+0.0001;	
	double prevslew = 0.105;
	double intensity = (pow(C,6)*15)+0.0001;
	double depthA = (pow(D,4)*940)+0.00001;
	int offsetA = (int)(depthA * overallscale);
	if (offsetA < 1) offsetA = 1;
	if (offsetA > 4880) offsetA = 4880;
	double balanceB = 0.0001;	
	slewgain *= overallscale;
	prevslew *= overallscale;
	balanceB /= overallscale;
	double outputgain = E;
	double wet = F;
	double balanceA = 1.0 - balanceB;
	//removed extra dry variable
	
	Vec4ui fpd; fpd.load(fpd_b);
	Vec2d control; control.load(control_b);
	Vec2d lastSample; lastSample.load(lastSample_b);
	Vec2d lastOutSample; lastOutSample.load(lastOutSample_b);
	Vec2d lastSlew; lastSlew.load(lastSlew_b);

	while (sampleFrames >= 0)
	{
		if (gcount > 19600) {gcount = 9800;}

		Vec2d inputSample1, inputSample2;
		Vec4d a, b, c, d, drySample;
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		drySample = inputSample;
		a = abs(inputSample)*intensity;
		a.store(d_b + gcount); a.store(d_b + gcount - 9800);
		b = a / offsetA;
		c.load(d_b + gcount - (offsetA * 2));
		b -= c / offsetA;
		b -= 0.000001;
		c = concatenate2(control + b.get_low(), control + b.get_low() + b.get_high()); // controls
		c = max(c, 0.0);
		d = max(select(c > 1.0, 1.0 - (c - 1.0), 1.0), 0.5); // clamp
		c = min(c, 1.0);
		control = c.get_high(); // new control
		a = ((1.0 - c) * 2.0) - 1.0; // thickness
		b = abs(a); // out
		c = min(abs(inputSample), 1.57079633); // bridgerectifier
		c = select(c > 0, sin(c), 1-cos(c));
		inputSample = select(inputSample > 0, (inputSample*(1-b))+(c*b), (inputSample*(1-b))-(c*b));
		inputSample *= d;
		a = concatenate2(lastSample, inputSample.get_low()); // lastSamples
		lastSample = inputSample.get_high(); // new lastSample
		b = inputSample - a; // slew
		c = abs(b*slewgain); // bridgerectifier
		c = select(c > 1.57079633, 1.0, sin(c));
		b = select(b > 0.0, c / slewgain, (c / slewgain) * -1.0);
		inputSample1 = (lastOutSample * balanceA) + (a.get_low() * balanceB) + b.get_low();
		lastOutSample = inputSample1;
		inputSample2 = (lastOutSample * balanceA) + (a.get_high() * balanceB) + b.get_high();
		lastOutSample = inputSample2;
		inputSample = concatenate2(inputSample1, inputSample2);
		d = min(abs(drySample + a), 1.0); // combSample
		c = concatenate2(lastSlew, b.get_low()); // lastSlews
		lastSlew = b.get_high(); // new lastSlew
		inputSample -= (c * d * prevslew);
		inputSample *= gain;
		c = abs(inputSample); // bridgerectifier
		c = select(c > 1.57079633, 1.0, sin(c));
		inputSample = select(inputSample > 0.0, c, c * -1.0);
		inputSample /= gain;
		inputSample *= gaintrim;

		if (outputgain != 1.0) {
			inputSample *= outputgain;
		}

		if (wet !=1.0) {
			inputSample = (inputSample * wet) + (drySample * (1.0-wet));
		}

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
		gcount +=4;
	}

	fpd.store(fpd_b);
	control.store(control_b);
	lastSample.store(lastSample_b);
	lastOutSample.store(lastOutSample_b);
	lastSlew.store(lastSlew_b);
}

void Desk4V::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
