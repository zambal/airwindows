/* ========================================
 *  ConsoleLAChannelV - ConsoleLAChannelV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleLAChannelV_H
#include "ConsoleLAChannelV.h"
#endif

void ConsoleLAChannelV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
  float* in1  =  inputs[0];
  float* in2  =  inputs[1];
  float* out1 = outputs[0];
  float* out2 = outputs[1];

	VstInt32 inFramesToProcess = sampleFrames; //vst doesn't give us this as a separate variable so we'll make it
	double overallscale = 1.0;
	overallscale /= 44100.0;
	overallscale *= getSampleRate();
	int cycleEnd = floor(overallscale);
	if (cycleEnd < 1) cycleEnd = 1;
	if (cycleEnd > 4) cycleEnd = 4;
	
	int limit = 4*cycleEnd;
	double divisor = 2.0/limit;
	
	treble = (A*6.0)-3.0;
	midA = midB;
	midB = (B*6.0)-3.0;
	bassA = bassB;
	bassB = (C*6.0)-3.0;
	//these should stack to go up to -3.0 to 3.0
	if (treble < 0.0) treble /= 3.0;
	if (midB < 0.0) midB /= 3.0;
	if (bassB < 0.0) bassB /= 3.0;
	//and then become -1.0 to 3.0;
	//there will be successive sin/cos stages w. dry/wet in these

	double freqMid = 0.024;
	switch (cycleEnd)
	{
		case 1: //base sample rate, no change
			break;
		case 2: //96k tier
			freqMid = 0.012;
			break;
		case 3: //192k tier
			freqMid = 0.006;
			break;
	}
	
	gainA = gainB;
	gainB = E*2.0; //smoothed master fader from Z2 filters
	//BitShiftGain pre gain trim goes here
	double subTrim = 0.0011 / overallscale;
	
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d g(gainL, gainR, gainL, gainR);
  while (sampleFrames > 0)
  {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		
		double temp = (double)sampleFrames/inFramesToProcess;
		double gain = (gainA*temp)+(gainB*(1.0-temp));
		double mid = (midA*temp)+(midB*(1.0-temp));
		double bass = (bassA*temp)+(bassB*(1.0-temp));
		
		Vec4d b(0.0), t;
		//begin Hull2 Treble
		{
			inputSample.store_a(hull + hullp);
			if(hullp == 0) inputSample.store_a(hull + 256);
			int x = hullp;
			while (x > hullp-limit) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(limit * 2)) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			b.store_a(hull + ((hullp - 40) & 255));
			if(((hullp - 40) & 255) == 0) b.store_a(hull + 256);
			x = hullp-40;
			b = 0.0;
			while (x > hullp-(40+limit)) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(40+(limit * 2))) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			b.store_a(hull + ((hullp - 80) & 255));
			if(((hullp - 80) & 255) == 0) b.store_a(hull + 256);
			x = hullp-80;
			b = 0.0;
			while (x > hullp-(80+limit)) {
				t.load_a(hull + (x & 255));
				b += t * divisor;
				x-=2;
			}
			b += b * 0.125;
			while (x > hullp-(80+(limit * 2))) {
				t.load_a(hull + (x & 255));
				b -= t * 0.125 * divisor;
				x-=2;
			}
			// t = inputSample - b; inputSample = t;
			inputSample = inputSample - b;
		}
		hullp = (hullp + 4) & 255;

		//end Hull2 treble
		
		// //begin Pear filter stages
		// {
		// 	//at this point 'bass' is actually still mid and bass
		// 	Vec2d slew; Vec2d b0 = b.get_low(); Vec2d b1 = b.get_high();
		// 	Vec2d p0, p1;
		// 	p0.load_a(pearB);p1.load_a(pearB + 2);slew = ((b0 - p0) + p1)*freqMid*0.5;
		// 	b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
		// 	p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
		// 	b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
		// 	b1.store_a(pearB); slew.store_a(pearB + 2);

		// 	p0.load_a(pearB + 4);p1.load_a(pearB + 6);slew = ((b0 - p0) + p1)*freqMid*0.5;
		// 	b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
		// 	p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
		// 	b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
		// 	b1.store_a(pearB + 4); slew.store_a(pearB + 6);

		// 	p0.load_a(pearB + 8);p1.load_a(pearB + 10); slew = ((b0 - p0) + p1)*freqMid*0.5;
		// 	b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
		// 	p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
		// 	b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
		// 	b1.store_a(pearB + 8); slew.store_a(pearB + 10);

		// 	p0.load_a(pearB + 12);p1.load_a(pearB + 14); slew = ((b0 - p0) + p1)*freqMid*0.5;
		// 	b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
		// 	p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
		// 	b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
		// 	b1.store_a(pearB + 12); slew.store_a(pearB + 14);

		// 	p0.load_a(pearB + 16);p1.load_a(pearB + 18); slew = ((b0 - p0) + p1)*freqMid*0.5;
		// 	b0 = (freqMid * b0) + ((1.0-freqMid) * (p0 + p1));
		// 	p0 = b0; p1 = slew; slew = ((b1 - p0) + p1)*freqMid*0.5;
		// 	b1 = (freqMid * b1) + ((1.0-freqMid) * (p0 + p1));
		// 	b1.store_a(pearB + 16); slew.store_a(pearB + 18);

		// 	b = concatenate2(b0, b1);
		// }
		// Vec4d m = inputSample - b;
		// //we now have three bands out of hull and pear filters
		
		// double w = 0.0; //filter into bands, apply the sin/cos to each band
		// if (treble > 0.0) {
		// 	w = treble > 1.0 ? 1.0 : treble;
		// 	t = (t*(1.0-w)) + (sin(t*M_PI_2)*treble);
		// }
		// if (treble < 0.0) {
		// 	t = min(max(t, -1.0), 1.0);
		// 	w = -treble > 1.0 ? 1.0 : -treble;
		// 	t = (t*(1.0-w))+if_mul(t < 0, ((1.0-cos(abs(t)*w))*(1.0-w)), -1.0);
		// } //cosine stages for EQ or expansion
		
		// m = min(max(m, -1.0), 1.0);
		// if (mid > 0.0) {
		// 	w = mid > 1.0 ? 1.0 : mid;
		// 	m = (m*(1.0-w)) + (sin(m*M_PI_2)*mid);
		// }
		// if (mid < 0.0) {
		// 	w = -mid > 1.0 ? 1.0 : -mid;
		// 	m = (m*(1.0-w))+if_mul(m < 0, ((1.0-cos(abs(m)*w))*(1.0-w)), -1.0);
		// } //cosine stages for EQ or expansion
		
		// b = min(max(b, -1.0), 1.0);
		// if (bass > 0.0) {
		// 	w = bass > 1.0 ? 1.0 : bass;
		// 	b = (b*(1.0-w)) + (sin(b*M_PI_2)*bass);
		// }
		// if (bass < 0.0) {
		// 	w = -bass > 1.0 ? 1.0 : -bass;
		// 	b = (b*(1.0-w))+if_mul(b < 0, ((1.0-cos(abs(b)*w))*(1.0-w)), -1.0);
		// } //cosine stages for EQ or expansion
		
	
		// inputSample = (b + m + t)*g*gain;
		// //applies BitShiftPan pan section, and smoothed fader gain
		
		// //begin SubTight section
		// {	
		// 	Vec2d subA; subA.load_a(subA_b);
		// 	Vec2d subB; subB.load_a(subB_b);
		// 	Vec2d subC; subC.load_a(subC_b);
		// 	Vec4d subSample = inputSample * subTrim;
		// 	Vec2d s0 = subSample.get_low(); Vec2d s1 = subSample.get_high();
		// 	Vec2d scale = 0.5+abs(s0*0.5);
		// 	s0 = (subA+(sin(subA-s0)*scale));
		// 	subA = s0*scale;
		// 	scale = 0.5+abs(s0*0.5);
		// 	s0 = (subB+(sin(subB-s0)*scale));
		// 	subB = s0*scale;
		// 	scale = 0.5+abs(s0*0.5);
		// 	s0 = (subC+(sin(subC-s0)*scale));
		// 	subC = s0*scale;

		// 	scale = 0.5+abs(s1*0.5);
		// 	s1 = (subA+(sin(subA-s1)*scale));
		// 	subA = s1*scale;
		// 	scale = 0.5+abs(s1*0.5);
		// 	s1 = (subB+(sin(subB-s1)*scale));
		// 	subB = s1*scale;
		// 	scale = 0.5+abs(s1*0.5);
		// 	s1 = (subC+(sin(subC-s1)*scale));
		// 	subC = s1*scale;

		// 	subSample = min(max(concatenate2(s0, s1), -0.25), 0.25) * 16.0;
		// 	inputSample += subSample;

		// 	subA.store_a(subA_b);
		// 	subB.store_a(subB_b);
		// 	subC.store_a(subC_b);
		// }
		// //end SubTight section		

		//begin Console7 Channel processing
		inputSample = min(max(inputSample, -1.097), 1.097);
		Vec4d absSample = abs(inputSample);
		inputSample = ((sin(inputSample*absSample)/select(absSample == 0.0, 1, absSample))*0.8)+(sin(inputSample)*0.2);
		//this is a version of Spiral blended 80/20 with regular Density.
		//It's blending between two different harmonics in the overtones of the algorithm
		
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

void ConsoleLAChannelV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
