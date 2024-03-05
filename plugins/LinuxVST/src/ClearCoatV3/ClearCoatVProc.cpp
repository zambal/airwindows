/* ========================================
 *  ClearCoatV - ClearCoatV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ClearCoatV_H
#include "ClearCoatV.h"
#endif

void ClearCoatV::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
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
	// if (cycle > ((cycleEnd-1) * 4)) cycle = (cycleEnd * 4)-4; //sanity check
	 
	double subRate = 0.002 / overallscale;
	double wet = B*2.0;
	double dry = 2.0 - wet;
	if (wet > 1.0) wet = 1.0;
	if (wet < 0.0) wet = 0.0;
	if (dry > 1.0) dry = 1.0;
	if (dry < 0.0) dry = 0.0;
	//this reverb makes 50% full dry AND full wet, not crossfaded.
	//that's so it can be on submixes without cutting back dry channel when adjusted:
	//unless you go super heavy, you are only adjusting the added verb loudness.
	
	double result[4] __attribute__((aligned(32)));
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d v, fpd_d = to_double(fpd);

  while (sampleFrames > 0)
  {
		Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
		inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
		Vec4d drySample = inputSample;

		cycle+=2;
		if (cycle >= cycleEnd) { //hit the end point and we do a reverb sample
			Vec2d o1, o2, o3, o4, vv;
						
			vv = inputSample.get_low();
			
			o1.load_a(feedbackA); 
			o2.load_a(feedbackB); 
			o3.load_a(feedbackC); 
			o4.load_a(feedbackD); 

			o1 = (vv + (o1 * 0.04166666666)); 
			o2 = (vv + (o2 * 0.04166666666)); 
			o3 = (vv + (o3 * 0.04166666666)); 
			o4 = (vv + (o4 * 0.04166666666)); 

			o1.store_a(aA + (gcount & 2047));
			o2.store_a(aB + (gcount & 4095));
			o3.store_a(aC + (gcount & 4095));
			o4.store_a(aD + (gcount & 8191));

			o1 = Vec2d(aA[(gcount - shortA) & 2047], aA[(gcount + 1 - shortD) & 2047]);
			o2 = Vec2d(aB[(gcount - shortB) & 4095], aB[(gcount + 1 - shortH) & 4095]);
			o3 = Vec2d(aC[(gcount - shortC) & 4095], aC[(gcount + 1 - shortL) & 4095]);
			o4 = Vec2d(aD[(gcount - shortD) & 8191], aD[(gcount + 1 - shortP) & 8191]);

			(o1 - (o2 + o3 + o4)).store_a(aE + (gcount & 4095));
			(o2 - (o1 + o3 + o4)).store_a(aF + (gcount & 4095));
			(o3 - (o1 + o2 + o4)).store_a(aG + (gcount & 4095));
			(o4 - (o1 + o2 + o3)).store_a(aH + (gcount & 4095));

			o1 = Vec2d(aE[(gcount - shortE) & 4095], aE[(gcount + 1 - shortC) & 4095]);
			o2 = Vec2d(aF[(gcount - shortF) & 4095], aF[(gcount + 1 - shortG) & 4095]);
			o3 = Vec2d(aG[(gcount - shortG) & 4095], aG[(gcount + 1 - shortK) & 4095]);
			o4 = Vec2d(aH[(gcount - shortH) & 4095], aH[(gcount + 1 - shortO) & 4095]);

			(o1 - (o2 + o3 + o4)).store_a(aI + (gcount & 4095));
			(o2 - (o1 + o3 + o4)).store_a(aJ + (gcount & 8191));
			(o3 - (o1 + o2 + o4)).store_a(aK + (gcount & 8191));
			(o4 - (o1 + o2 + o3)).store_a(aL + (gcount & 8191));
	
			o1 = Vec2d(aI[(gcount - shortI) & 4095], aI[(gcount + 1 - shortB) & 4095]);
			o2 = Vec2d(aJ[(gcount - shortJ) & 8191], aJ[(gcount + 1 - shortF) & 8191]);
			o3 = Vec2d(aK[(gcount - shortK) & 8191], aK[(gcount + 1 - shortJ) & 8191]);
			o4 = Vec2d(aL[(gcount - shortL) & 8191], aL[(gcount + 1 - shortN) & 8191]);
	
			(o1 - (o2 + o3 + o4)).store_a(aM + (gcount & 2047));
			(o2 - (o1 + o3 + o4)).store_a(aN + (gcount & 8191));
			(o3 - (o1 + o2 + o4)).store_a(aO + (gcount & 2047));
			(o4 - (o1 + o2 + o3)).store_a(aP + (gcount & 8191));
	
			o1 = Vec2d(aM[(gcount - shortM) & 2047], aM[(gcount + 1 - shortA) & 2047]);
			o2 = Vec2d(aN[(gcount - shortN) & 8191], aN[(gcount + 1 - shortE) & 8191]);
			o3 = Vec2d(aO[(gcount - shortO) & 2047], aO[(gcount + 1 - shortI) & 2047]);
			o4 = Vec2d(aP[(gcount - shortP) & 8191], aP[(gcount + 1 - shortM) & 8191]);
	
			vv.load_a(prevMulchA);
			o1.store_a(prevMulchA);
			o1 = (o1 + o1 + o1 + vv)*0.25;

			(o1 - (o2 + o3 + o4)).store_a(feedbackA);
			(o2 - (o1 + o3 + o4)).store_a(feedbackB);
			(o3 - (o1 + o2 + o4)).store_a(feedbackC);
			(o4 - (o1 + o2 + o3)).store_a(feedbackD);

			vv = min(max((o1 + o2 + o3 + o4)*0.125, -1.0), 1.0);
			//and take the final combined sum of outputs, corrected for Householder gain

			gcount = (gcount + 2) & 8191;

			if (cycleEnd == 4) {
				Vec2d tt; tt.load_a(lastRef + 4);
				((tt + vv) * 0.5).store_a(lastRef);
				((tt + vv + vv + vv) * 0.25).store_a(lastRef + 2);
				vv.store_a(lastRef + 4);
				inputSample = concatenate2(tt, (tt + tt + tt + vv) * 0.25);
			}
			if (cycleEnd <= 2) {
				Vec2d tt; tt.load_a(lastRef);
				vv.store_a(lastRef);
				inputSample = concatenate2(tt, (tt + vv) * 0.5);
			}
			cycle = 0; //reset
			
		} else {
			inputSample.load_a(lastRef);
			//we are going through our references now
		}
	
		//begin SubTight section
		Vec4d subSample = inputSample * subRate;
		Vec4d subA; subA.load_a(cc_sub_b);
		Vec4d subB; subB.load_a(cc_sub_b + 4);
		Vec4d subC; subC.load_a(cc_sub_b + 8);
		Vec4d subD; subD.load_a(cc_sub_b + 12);
		Vec4d scale = 0.5+abs(subSample*0.5);
		subSample = (subA+(sin(subA-subSample)*scale));
		subA = subSample*scale;
		scale = 0.5+abs(subSample*0.5);
		subSample = (subB+(sin(subB-subSample)*scale));
		subB = subSample*scale;
		scale = 0.5+abs(subSample*0.5);
		subSample = (subC+(sin(subC-subSample)*scale));
		subC = subSample*scale;
		scale = 0.5+abs(subSample*0.5);
		subSample = (subD+(sin(subD-subSample)*scale));
		subD = subSample*scale;
		subA.store_a(cc_sub_b);
		subB.store_a(cc_sub_b + 4);
		subC.store_a(cc_sub_b + 8);
		subD.store_a(cc_sub_b + 12);

		subSample = min(max(subSample, -0.25), 0.25);
		inputSample -= (subSample*16.0);
		//end SubTight section		

		if (cycleEnd > 1) {
			Vec2d tt; tt.load_a(tail); inputSample.get_high().store_a(tail);
			inputSample = (inputSample + concatenate2(tt, inputSample.get_low())) * 0.5;
		} //let's average only at elevated sample rates
	
		if (wet < 1.0) {inputSample *= wet;}
		if (dry < 1.0) {drySample *= dry;}
		inputSample += drySample;
		//this is our submix verb dry/wet: 0.5 is BOTH at FULL VOLUME
		//purpose is that, if you're adding verb, you're not altering other balances
	
		//begin 32 bit stereo floating point dither
		fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
		fpd_d = to_double(fpd);
		Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample)) + 62);
		inputSample += (fpd_d - 2147483647.0) * to_double(exp) * 5.5e-36l;
		//end 32 bit stereo floating point dither

		inputSample.store_a(result);
		*out1 = result[0];
		*(out1 + 1) = result[2];
		*out2 = result[1];
		*(out2 + 1) = result[3];

		in1+= 2;
		in2+= 2;
		out1+= 2;
		out2+= 2;

		sampleFrames -= 2;
  }
	fpd.store_a(fpd_b);
}

void ClearCoatV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
