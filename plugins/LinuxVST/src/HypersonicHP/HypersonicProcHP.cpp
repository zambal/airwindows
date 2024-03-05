/* ========================================
 *  HypersonicHP - HypersonicHP.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __HypersonicHP_H
#include "HypersonicHP.h"
#endif

void HypersonicHP::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

		Vec4ui fpd; fpd.load(fpd_b);
    Vec2d sA1; sA1.load(fixA + fix_sL1);	
    Vec2d sA2; sA2.load(fixA + fix_sL2);	
    Vec2d sB1; sB1.load(fixB + fix_sL1);	
    Vec2d sB2; sB2.load(fixB + fix_sL2);	
    Vec2d sC1; sC1.load(fixC + fix_sL1);	
    Vec2d sC2; sC2.load(fixC + fix_sL2);	
    Vec2d sD1; sD1.load(fixD + fix_sL1);	
    Vec2d sD2; sD2.load(fixD + fix_sL2);	
    Vec2d sE1; sE1.load(fixE + fix_sL1);	
    Vec2d sE2; sE2.load(fixE + fix_sL2);	
    Vec2d sF1; sF1.load(fixF + fix_sL1);	
    Vec2d sF2; sF2.load(fixF + fix_sL2);	
    Vec2d sG1; sG1.load(fixG + fix_sL1);	
    Vec2d sG2; sG2.load(fixG + fix_sL2);	

    while (sampleFrames >= 0)
    {
        Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
        inputSample = select(abs(inputSample) < 1.18e-23, to_double(fpd) * 1.18e-17, inputSample); 
        Vec2d innerSample[2]; innerSample[0] = inputSample.get_low(); innerSample[1] = inputSample.get_high();
        for(int i = 0; i < 2; ++i) {
          Vec2d temp = (innerSample[i] * fixA[fix_a0]) + sA1;
          sA1 = (innerSample[i] * fixA[fix_a1]) - (temp * fixA[fix_b1]) + sA2;
          sA2 = (innerSample[i] * fixA[fix_a2]) - (temp * fixA[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics
	
          temp = (innerSample[i] * fixB[fix_a0]) + sB1;
          sB1 = (innerSample[i] * fixB[fix_a1]) - (temp * fixB[fix_b1]) + sB2;
          sB2 = (innerSample[i] * fixB[fix_a2]) - (temp * fixB[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics

          temp = (innerSample[i] * fixC[fix_a0]) + sC1;
          sC1 = (innerSample[i] * fixC[fix_a1]) - (temp * fixC[fix_b1]) + sC2;
          sC2 = (innerSample[i] * fixC[fix_a2]) - (temp * fixC[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics

          temp = (innerSample[i] * fixD[fix_a0]) + sD1;
          sD1 = (innerSample[i] * fixD[fix_a1]) - (temp * fixD[fix_b1]) + sD2;
          sD2 = (innerSample[i] * fixD[fix_a2]) - (temp * fixD[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics

          temp = (innerSample[i] * fixE[fix_a0]) + sE1;
          sE1 = (innerSample[i] * fixE[fix_a1]) - (temp * fixE[fix_b1]) + sE2;
          sE2 = (innerSample[i] * fixE[fix_a2]) - (temp * fixE[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics

          temp = (innerSample[i] * fixF[fix_a0]) + sF1;
          sF1 = (innerSample[i] * fixF[fix_a1]) - (temp * fixF[fix_b1]) + sF2;
          sF2 = (innerSample[i] * fixF[fix_a2]) - (temp * fixF[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics

          temp = (innerSample[i] * fixG[fix_a0]) + sG1;
          sG1 = (innerSample[i] * fixG[fix_a1]) - (temp * fixG[fix_b1]) + sG2;
          sG2 = (innerSample[i] * fixG[fix_a2]) - (temp * fixG[fix_b2]);
          innerSample[i] = temp; //fixed biquad filtering ultrasonics
    		}

        inputSample = concatenate2(innerSample[0], innerSample[1]);
		
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

    fpd.store(fpd_b);
    sA1.store(fixA + fix_sL1);	
    sA2.store(fixA + fix_sL2);	
    sB1.store(fixB + fix_sL1);	
    sB2.store(fixB + fix_sL2);	
    sC1.store(fixC + fix_sL1);	
    sC2.store(fixC + fix_sL2);	
    sD1.store(fixD + fix_sL1);	
    sD2.store(fixD + fix_sL2);	
    sE1.store(fixE + fix_sL1);	
    sE2.store(fixE + fix_sL2);	
    sF1.store(fixF + fix_sL1);	
    sF2.store(fixF + fix_sL2);	
    sG1.store(fixG + fix_sL1);	
    sG2.store(fixG + fix_sL2);	
}

void HypersonicHP::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];
	
}
