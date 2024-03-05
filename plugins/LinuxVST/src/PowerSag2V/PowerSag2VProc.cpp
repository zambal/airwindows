/* ========================================
 *  PowerSag2V - PowerSag2V.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __PowerSag2V_H
#include "PowerSag2V.h"
#endif

void PowerSag2V::processReplacing(float **inputs, float **outputs, VstInt32 sampleFrames) 
{
    float* in1  =  inputs[0];
    float* in2  =  inputs[1];
    float* out1 = outputs[0];
    float* out2 = outputs[1];

    double depth = pow(A,4);
    int offset = (int)(depth * 16383) + 1;	
    double sqrt_offset = sqrt(offset);
    double wet = (B*2.0)-1.0;
    
    while (sampleFrames >= 0)
    {
        Vec4ui fpd_a = fpd;
        fpd_a = fpd_a ^ (fpd_a << 13); fpd_a = fpd_a ^ (fpd_a << 17); fpd_a = fpd_a ^ (fpd_a << 5);
        Vec4d fpd_mult = to_double(blend4<0, 1, 4, 5>(fpd, fpd_a)) * 1.18e-17;
        Vec4d inputSample(*in1, *in2, *(in1 + 1), *(in2 + 1));
        inputSample = select(abs(inputSample) < 1.18e-23, fpd_mult, inputSample); 
        Vec4d drySample = inputSample;

    		if (gcount > 32766) {gcount = 0;}		

    		Vec4d v = abs(inputSample);
        v.store(d + gcount);
    		control += v;
        if(gcount-(offset*2) < 3)
            v.load(d + gcount-(offset*2) + 32770);
        else if(gcount-(offset*2) < 0) {
            Vec2d v1; v1.load(d + gcount-(offset*2) + 32770);
            Vec2d v2; v2.load(d + gcount-(offset*2) + 2);
            v = concatenate2(v1, v2);
        }
        else
            v.load(d + gcount-(offset*2));

        gcount+=4;
        control -= v;
        control = max(min(control, offset), 0.0);
    		v = inputSample * (control / sqrt_offset);
        v = inputSample / select(v == 0.0, 1.0, v);
        v = max(min(v, 1.0), 0.0);
    		inputSample *= v;
        v = drySample - inputSample;
        if(wet < 0.0) drySample *= (wet+1.0);
    		inputSample = drySample; // - (v * wet);
		
		
        //begin 32 bit stereo floating point dither
        Vec4d exp2(2.0);
        Vec4i expon = exponent(to_float(inputSample)) + 62;
        Vec4d mult = pow(exp2, to_double(expon)) * 5.5e-36l;
        fpd = fpd_a ^ (fpd_a << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
        inputSample += (to_double(blend4<0, 1, 4, 5>(fpd_a, fpd)) - Vec4d(2147483647.0)) * mult;

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
}

void PowerSag2V::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
                    
