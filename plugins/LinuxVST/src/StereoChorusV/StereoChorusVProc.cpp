/* ========================================
 *  StereoChorusV - StereoChorusV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __StereoChorusV_H
#include "StereoChorusV.h"
#endif

#define BL 65536
#define BM 65535
void StereoChorusV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
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
	if (cycle > cycleEnd-1) cycle = cycleEnd-1; //sanity check
	
	double speed = pow(0.32+(A/6),10);
	double depth = (B/60) / speed;
	double tupi = 3.141592653589793238 * 2.0;
	
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d fpd_d = to_double(fpd);
    while (sampleFrames > 0)
    {
	  Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
	  inputSample = select(abs(inputSample) < 1.18e-23, fpd_d * 1.18e-17, inputSample); 
		
		 cycle++;
		 if (cycle == cycleEnd) { //hit the end point and we do a chorus sample
		 	//assign working variables
		 	Vec2d airPrev; airPrev.load_a(airPrev_b);
			inputSample.get_high().store_a(airPrev_b);
		 	Vec2d airFactor = airPrev - inputSample.get_low();
		 	Vec2d airEven; airEven.load_a(airEven_b);
		 	Vec2d airOdd; airOdd.load_a(airOdd_b);
			airEven += airFactor; airOdd -= airFactor; airFactor = airEven;
		 	airOdd = ((airOdd - ((airOdd - airEven)*0.00390625)) * 0.9999);
		 	airEven = ((airEven - ((airEven - airOdd)*0.00390625)) * 0.9999);
		 	Vec2d t = airFactor;
		 	airFactor = inputSample.get_low() - inputSample.get_high();
			airEven -= airFactor; airOdd += airFactor; airFactor = airOdd;
		 	airOdd = ((airOdd - ((airOdd - airEven)*0.00390625)) * 0.9999);
		 	airEven = ((airEven - ((airEven - airOdd)*0.00390625)) * 0.9999);
			airEven.store(airEven_b); airOdd.store(airOdd_b);
			inputSample += concatenate2(t, airFactor);
			//air, compensates for loss of highs in flanger's interpolation
			
		  truncate_to_int32(inputSample*8388352.0).store(p_b + gcount);
		  Vec4i count(gcount, gcount + 1, gcount + 2, gcount + 3);
		  Vec2d prev_sweep; prev_sweep.load_a(sweep_b);
		  Vec4d sweep(prev_sweep, if_sub(prev_sweep + speed > tupi, prev_sweep + speed, tupi));
		  Vec4d offset = depth + (depth * sin(sweep));
		  count -= (truncate_to_int32(floor(offset)) * 2);
		  int count_b[4] __attribute__((aligned(32))); count.store_a(count_b);
		  Vec4i temp1(p_b[count_b[0] & BM], p_b[count_b[1] & BM], p_b[count_b[2] & BM], p_b[count_b[3] & BM]);
		  Vec4i temp2(p_b[(count_b[0] - 2) & BM], p_b[(count_b[1] - 2) & BM], p_b[(count_b[2] - 2) & BM], p_b[(count_b[3] - 2) & BM]);
		  Vec4i temp3(p_b[(count_b[0] - 4) & BM], p_b[(count_b[1] - 4) & BM], p_b[(count_b[2] - 4) & BM], p_b[(count_b[3] - 4) & BM]);
		 
		  temp1 = truncate_to_int32(to_double(temp1) * (1.0 - (offset - floor(offset))));
		  temp1 += temp2;
		  temp1 += truncate_to_int32(to_double(temp3) * (offset - floor(offset)));
		  temp1 -= ((temp1 - temp2)-(temp2-temp3)) / 50;
		 
		  sweep += speed;
		  sweep = if_sub(sweep > tupi, sweep, tupi);
		  sweep.get_high().store_a(sweep_b);
		 
		  gcount = (gcount + 4) & BM;
		  //still scrolling through the samples, remember

		  inputSample = to_double(temp1)/16776704.0;
		 
		 	if (cycleEnd == 4) {
		 		Vec2d l1; l1.load_a(lastRef + 16);
		 		Vec2d l2 = inputSample.get_low();
		 		Vec2d l3 = inputSample.get_high();

		 		l1.store_a(lastRef);
		 		((l1 + l1 + l1 + l2) * 0.25).store_a(lastRef + 2);
		 		((l1 + l2) / 2.0).store_a(lastRef + 4);
		 		((l1 + l2 + l2 + l2) * 0.25).store_a(lastRef + 6);
		 		l2.store_a(lastRef + 8);
		 		((l2 + l2 + l2 + l3) * 0.25).store_a(lastRef + 10);
		 		((l2 + l3) / 2.0).store_a(lastRef + 12);
		 		((l2 + l3 + l3 + l3) * 0.25).store_a(lastRef + 14);
		 		l3.store_a(lastRef + 16);
		 	}
		 	if (cycleEnd == 3) {
		 		Vec2d l1; l1.load_a(lastRef + 12);
		 		Vec2d l2 = inputSample.get_low();
		 		Vec2d l3 = inputSample.get_high();

		 		l1.store_a(lastRef);
		 		((l1 + l1 + l2) / 3.0).store_a(lastRef + 2);
		 		((l1 + l2 + l2) / 3.0).store_a(lastRef + 4);
		 		l2.store_a(lastRef + 6);
		 		((l2 + l2 + l3) / 3.0).store_a(lastRef + 8);
		 		((l2 + l3 + l3) / 3.0).store_a(lastRef + 10);
		 		l3.store_a(lastRef + 12);
		 	}
		 	if (cycleEnd == 2) {
		 		Vec2d l1; l1.load_a(lastRef + 8);
		 		Vec2d l2 = inputSample.get_low();
		 		Vec2d l3 = inputSample.get_high();

		 		l1.store_a(lastRef);
		 		((l1 + l2) / 2.0).store_a(lastRef + 2);
		 		l2.store_a(lastRef + 4);
		 		((l2 + l3) / 2.0).store_a(lastRef + 6);
		 		l3.store_a(lastRef + 8);
		 	}
		 	if (cycleEnd == 1) {
		 		inputSample.store_a(lastRef);
		 	}
		 	cycle = 0; //reset
		 	inputSample.load_a(lastRef);
		 } else {
		 	inputSample.load_a(lastRef + (cycle * 4));
		 	//we are going through our references now
		 }
		
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

void StereoChorusV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
