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
	 
	int shorts1[8] __attribute__((aligned(32))) = {shortA + 3, shortB + 3, shortC + 3, shortD + 3, shortD + 3, shortH + 3, shortL + 3, shortP + 3};
	int shorts2[8] __attribute__((aligned(32))) = {shortE + 3, shortF + 3, shortG + 3, shortH + 3, shortC + 3, shortG + 3, shortK + 3, shortO + 3};
	int shorts3[8] __attribute__((aligned(32))) = {shortI + 3, shortJ + 3, shortK + 3, shortL + 3, shortB + 3, shortF + 3, shortJ + 3, shortN + 3};
	int shorts4[8] __attribute__((aligned(32))) = {shortM + 3, shortN + 3, shortO + 3, shortP + 3, shortA + 3, shortE + 3, shortI + 3, shortM + 3};
	int load_c[8] __attribute__((aligned(32))), copy_c[8] __attribute__((aligned(32))), copy_s[8] __attribute__((aligned(32)));

	double subRate = 0.001 / overallscale;
	double wet = B*2.0;
	double dry = 2.0 - wet;
	if (wet > 1.0) wet = 1.0;
	if (wet < 0.0) wet = 0.0;
	if (dry > 1.0) dry = 1.0;
	if (dry < 0.0) dry = 0.0;
	//this reverb makes 50% full dry AND full wet, not crossfaded.
	//that's so it can be on submixes without cutting back dry channel when adjusted:
	//unless you go super heavy, you are only adjusting the added verb loudness.
	
	Vec4ui fpd; fpd.load_a(fpd_b);
	Vec4d inputSample[2], v;
	double drySample_b[8] __attribute__((aligned(32)));

  while (sampleFrames > 0)
  {
		Vec4d l(*in1, *(in1 + 1), *(in1 + 2), *(in1 + 3));
		l = select(abs(l) < 1.18e-23, to_double(fpd) * 1.18e-17, l); 
		Vec4d r(*in2, *(in2 + 1), *(in2 + 2), *(in2 + 3));
		r = select(abs(r) < 1.18e-23, to_double(fpd) * 1.18e-17, r); 
		blend4<0,4,1,5>(l, r).store_a(drySample_b);
		blend4<2,6,3,7>(l, r).store_a(drySample_b + 4);
		cycle+=4;
		if (cycle == (cycleEnd * 4)) { //hit the end point and we do a reverb sample
			Vec4d o1, o2, o3, o4;
			Vec8i c, s, ks, t;
			
			v.load_a(feedbackAL); (l + (v * 0.04166666666)).store(aAL + counters1[0]);
			v.load_a(feedbackBL); (l + (v * 0.04166666666)).store(aBL + counters1[1]);
			v.load_a(feedbackCL); (l + (v * 0.04166666666)).store(aCL + counters1[2]);
			v.load_a(feedbackDL); (l + (v * 0.04166666666)).store(aDL + counters1[3]);

			v.load_a(feedbackDR); (r + (v * 0.04166666666)).store(aDR + counters1[4]);
			v.load_a(feedbackHR); (r + (v * 0.04166666666)).store(aHR + counters1[5]);
			v.load_a(feedbackLR); (r + (v * 0.04166666666)).store(aLR + counters1[6]);
			v.load_a(feedbackPR); (r + (v * 0.04166666666)).store(aPR + counters1[7]);
	
			c.load_a(counters1);
			s.load_a(shorts1);
			select((c + 4) > s, c - s + 3, c + 4).store_a(load_c);
			t = select(c > (s - 3), max(s + 1, c), -1); t.store_a(copy_c);
			select((c > (s - 3)), (c + 4 - t) * sizeof(double), -1).store_a(copy_s);
			select((c + 4) > (s + 3), c - s + 3, c + 4).store_a(counters1);

			if(copy_c[0] != -1) memcpy(aEL + copy_c[0] - (shorts2[0] + 1), aEL + copy_c[0], copy_s[0]);
			if(copy_c[1] != -1) memcpy(aFL + copy_c[1] - (shorts2[1] + 1), aFL + copy_c[1], copy_s[1]);
			if(copy_c[2] != -1) memcpy(aGL + copy_c[2] - (shorts2[2] + 1), aGL + copy_c[2], copy_s[2]);
			if(copy_c[3] != -1) memcpy(aHL + copy_c[3] - (shorts2[3] + 1), aHL + copy_c[3], copy_s[3]);

			if(copy_c[4] != -1) memcpy(aCR + copy_c[4] - (shorts2[4] + 1), aCR + copy_c[4], copy_s[4]);
			if(copy_c[5] != -1) memcpy(aGR + copy_c[5] - (shorts2[5] + 1), aGR + copy_c[5], copy_s[5]);
			if(copy_c[6] != -1) memcpy(aKR + copy_c[6] - (shorts2[6] + 1), aKR + copy_c[6], copy_s[6]);
			if(copy_c[7] != -1) memcpy(aOR + copy_c[7] - (shorts2[7] + 1), aOR + copy_c[7], copy_s[7]);
		

			o1.load(aAL + load_c[0]);
			o2.load(aBL + load_c[1]);
			o3.load(aCL + load_c[2]);
			o4.load(aDL + load_c[3]);
	
			(o1 - (o2 + o3 + o4)).store(aEL + counters2[0]);
			(o2 - (o1 + o3 + o4)).store(aFL + counters2[1]);
			(o3 - (o1 + o2 + o4)).store(aGL + counters2[2]);
			(o4 - (o1 + o2 + o3)).store(aHL + counters2[3]);

			o1.load(aDR + load_c[4]);
			o2.load(aHR + load_c[5]);
			o3.load(aLR + load_c[6]);
			o4.load(aPR + load_c[7]);

			(o1 - (o2 + o3 + o4)).store(aCR + counters2[4]);
			(o2 - (o1 + o3 + o4)).store(aGR + counters2[5]);
			(o3 - (o1 + o2 + o4)).store(aKR + counters2[6]);
			(o4 - (o1 + o2 + o3)).store(aOR + counters2[7]);

			c.load_a(counters2);
			s.load_a(shorts2);
			select((c + 4) > s, c - s + 3, c + 4).store_a(load_c);
			t = select(c > (s - 3), max(s + 1, c), -1); t.store_a(copy_c);
			select((c > (s - 3)), (c + 4 - t) * sizeof(double), -1).store_a(copy_s);
			select((c + 4) > (s + 3), c - s + 3, c + 4).store_a(counters2);
		
			if(copy_c[0] != -1) memcpy(aEL + copy_c[0] - (shorts2[0] + 1), aEL + copy_c[0], copy_s[0]);
			if(copy_c[1] != -1) memcpy(aFL + copy_c[1] - (shorts2[1] + 1), aFL + copy_c[1], copy_s[1]);
			if(copy_c[2] != -1) memcpy(aGL + copy_c[2] - (shorts2[2] + 1), aGL + copy_c[2], copy_s[2]);
			if(copy_c[3] != -1) memcpy(aHL + copy_c[3] - (shorts2[3] + 1), aHL + copy_c[3], copy_s[3]);

			if(copy_c[4] != -1) memcpy(aCR + copy_c[4] - (shorts2[4] + 1), aCR + copy_c[4], copy_s[4]);
			if(copy_c[5] != -1) memcpy(aGR + copy_c[5] - (shorts2[5] + 1), aGR + copy_c[5], copy_s[5]);
			if(copy_c[6] != -1) memcpy(aKR + copy_c[6] - (shorts2[6] + 1), aKR + copy_c[6], copy_s[6]);
			if(copy_c[7] != -1) memcpy(aOR + copy_c[7] - (shorts2[7] + 1), aOR + copy_c[7], copy_s[7]);
		
			o1.load(aEL + load_c[0]);
			o2.load(aFL + load_c[1]);
			o3.load(aGL + load_c[2]);
			o4.load(aHL + load_c[3]);
	
			(o1 - (o2 + o3 + o4)).store(aIL + counters3[0]);
			(o2 - (o1 + o3 + o4)).store(aJL + counters3[1]);
			(o3 - (o1 + o2 + o4)).store(aKL + counters3[2]);
			(o4 - (o1 + o2 + o3)).store(aLL + counters3[3]);
		
			o1.load(aCR + load_c[4]);
			o2.load(aGR + load_c[5]);
			o3.load(aKR + load_c[6]);
			o4.load(aOR + load_c[7]);

			(o1 - (o2 + o3 + o4)).store(aBR + counters3[4]);
			(o2 - (o1 + o3 + o4)).store(aFR + counters3[5]);
			(o3 - (o1 + o2 + o4)).store(aJR + counters3[6]);
			(o4 - (o1 + o2 + o3)).store(aNR + counters3[7]);
		
			c.load_a(counters3);
			s.load_a(shorts3);
			select((c + 4) > s, c - s + 3, c + 4).store_a(load_c);
			t = select(c > (s - 3), max(s + 1, c), -1); t.store_a(copy_c);
			select((c > (s - 3)), (c + 4 - t) * sizeof(double), -1).store_a(copy_s);
			select((c + 4) > (s + 3), c - s + 3, c + 4).store_a(counters3);

			if(copy_c[0] != -1) memcpy(aIL + copy_c[0] - (shorts3[0] + 1), aIL + copy_c[0], copy_s[0]);
			if(copy_c[1] != -1) memcpy(aJL + copy_c[1] - (shorts3[1] + 1), aJL + copy_c[1], copy_s[1]);
			if(copy_c[2] != -1) memcpy(aKL + copy_c[2] - (shorts3[2] + 1), aKL + copy_c[2], copy_s[2]);
			if(copy_c[3] != -1) memcpy(aLL + copy_c[3] - (shorts3[3] + 1), aLL + copy_c[3], copy_s[3]);

			if(copy_c[4] != -1) memcpy(aBR + copy_c[4] - (shorts3[4] + 1), aBR + copy_c[4], copy_s[4]);
			if(copy_c[5] != -1) memcpy(aFR + copy_c[5] - (shorts3[5] + 1), aFR + copy_c[5], copy_s[5]);
			if(copy_c[6] != -1) memcpy(aJR + copy_c[6] - (shorts3[6] + 1), aJR + copy_c[6], copy_s[6]);
			if(copy_c[7] != -1) memcpy(aNR + copy_c[7] - (shorts3[7] + 1), aNR + copy_c[7], copy_s[7]);
		

			o1.load(aIL + load_c[0]);
			o2.load(aJL + load_c[1]);
			o3.load(aKL + load_c[2]);
			o4.load(aLL + load_c[3]);
		
			(o1 - (o2 + o3 + o4)).store(aML + counters4[0]);
			(o2 - (o1 + o3 + o4)).store(aNL + counters4[1]);
			(o3 - (o1 + o2 + o4)).store(aOL + counters4[2]);
			(o4 - (o1 + o2 + o3)).store(aPL + counters4[3]);
		
			o1.load(aBR + load_c[4]);
			o2.load(aFR + load_c[5]);
			o3.load(aJR + load_c[6]);
			o4.load(aNR + load_c[7]);
		
			(o1 - (o2 + o3 + o4)).store(aAR + counters4[4]);
			(o2 - (o1 + o3 + o4)).store(aER + counters4[5]);
			(o3 - (o1 + o2 + o4)).store(aIR + counters4[6]);
			(o4 - (o1 + o2 + o3)).store(aMR + counters4[7]);
		
			c.load_a(counters4);
			s.load_a(shorts4);
			select((c + 4) > s, c - s + 3, c + 4).store_a(load_c);
			t = select(c > (s - 3), max(s + 1, c), -1); t.store_a(copy_c);
			select((c > (s - 3)), (c + 4 - t) * sizeof(double), -1).store_a(copy_s);
			select((c + 4) > (s + 3), c - s + 3, c + 4).store_a(counters4);

			if(copy_c[0] != -1) memcpy(aML + copy_c[0] - (shorts4[0] + 1), aML + copy_c[0], copy_s[0]);
			if(copy_c[1] != -1) memcpy(aNL + copy_c[1] - (shorts4[1] + 1), aNL + copy_c[1], copy_s[1]);
			if(copy_c[2] != -1) memcpy(aOL + copy_c[2] - (shorts4[2] + 1), aOL + copy_c[2], copy_s[2]);
			if(copy_c[3] != -1) memcpy(aPL + copy_c[3] - (shorts4[3] + 1), aPL + copy_c[3], copy_s[3]);

			if(copy_c[4] != -1) memcpy(aAR + copy_c[4] - (shorts4[4] + 1), aAR + copy_c[4], copy_s[4]);
			if(copy_c[5] != -1) memcpy(aER + copy_c[5] - (shorts4[5] + 1), aER + copy_c[5], copy_s[5]);
			if(copy_c[6] != -1) memcpy(aIR + copy_c[6] - (shorts4[6] + 1), aIR + copy_c[6], copy_s[6]);
			if(copy_c[7] != -1) memcpy(aMR + copy_c[7] - (shorts4[7] + 1), aMR + copy_c[7], copy_s[7]);

			o1.load(aML + load_c[0]);
			o2.load(aNL + load_c[1]);
			o3.load(aOL + load_c[2]);
			o4.load(aPL + load_c[3]);
		
			v = permute4<3, 0, 1, 2>(o1);
			double x = v.extract(0);
			v.insert(0, prevMulchAL);
			prevMulchAL = x;
			o1 = (o1 + o1 + o1 + v)*0.25;

			(o1 - (o2 + o3 + o4)).store(feedbackAL);
			(o2 - (o1 + o3 + o4)).store(feedbackBL);
			(o3 - (o1 + o2 + o4)).store(feedbackCL);
			(o4 - (o1 + o2 + o3)).store(feedbackDL);

			l = min(max((o1 + o2 + o3 + o4)/8.0, -1.0), 1.0);

			o1.load(aAR + load_c[4]);
			o2.load(aER + load_c[5]);
			o3.load(aIR + load_c[6]);
			o4.load(aMR + load_c[7]);			

			v = permute4<3, 0, 1, 2>(o1);
			x = v.extract(0);
			v.insert(0, prevMulchAR);
			prevMulchAR = x;
			o1 = (o1 + o1 + o1 + v)*0.25;
		
			(o1 - (o2 + o3 + o4)).store_a(feedbackDR);
			(o2 - (o1 + o3 + o4)).store_a(feedbackHR);
			(o3 - (o1 + o2 + o4)).store_a(feedbackLR);
			(o4 - (o1 + o2 + o3)).store_a(feedbackPR);
			//which we need to feed back into the input again, a bit
		
			r = min(max((o1 + o2 + o3 + o4)/8.0, -1.0), 1.0);
			// l = o1;
			// r = o1;
			//and take the final combined sum of outputs, corrected for Householder gain

			if (cycleEnd == 4) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[32]);
				lastRef[32] = x;
				scatter<0, 8, 16, 24>(v, lastRef);
				scatter<2, 10, 18, 26>((v + v + v + l) / 4.0, lastRef);
				scatter<4, 12, 20, 28>((v + v + l + l) / 4.0, lastRef);
				scatter<6, 14, 22, 30>((v + l + l + l) / 4.0, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[33]);
				lastRef[33] = x;
				scatter<1, 9, 17, 25>(v, lastRef);
				scatter<3, 11, 19, 27>((v + v + v + l) / 4.0, lastRef);
				scatter<5, 13, 21, 29>((v + v + l + l) / 4.0, lastRef);
				scatter<7, 15, 23, 31>((v + l + l + l) / 4.0, lastRef);
			}
			if (cycleEnd == 3) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[24]);
				lastRef[24] = x;
				scatter<0, 6, 12, 18>(v, lastRef);
				scatter<2, 8, 14, 20>((v + v + l) / 3.0, lastRef);
				scatter<4, 10, 16, 22>((v + l + l) / 3.0, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[25]);
				lastRef[25] = x;
				scatter<1, 7, 13, 19>(v, lastRef);
				scatter<3, 9, 15, 21>((v + v + r) / 3.0, lastRef);
				scatter<5, 11, 17, 23>((v + r + r) / 3.0, lastRef);
			}
			if (cycleEnd == 2) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[16]);
				lastRef[16] = x;
				scatter<0, 4, 8, 12>(v, lastRef);
				scatter<2, 6, 10, 14>((v + l) / 2.0, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[17]);
				lastRef[17] = x;
				scatter<1, 5, 9, 13>(v, lastRef);
				scatter<3, 7, 11, 15>((v + r) / 2.0, lastRef);
			}
			if (cycleEnd == 1) {
				scatter<0, 2, 4, 6>(l, lastRef);
				scatter<1, 3, 5, 7>(r, lastRef);
			}
			cycle = 0; //reset
			inputSample[0].load(lastRef);
			inputSample[1].load(lastRef + 4);
		} else {
			inputSample[0].load(lastRef + cycle);
			inputSample[1].load(lastRef + cycle + 4);
			//we are going through our references now
		}
	
		Vec4d drySample[2]; drySample[0].load_a(drySample_b); drySample[1].load_a(drySample_b + 4);
		for(int i = 0; i < 2; i ++) {
			{
				//begin SubTight section
				Vec2d subA; subA.load_a(subA_b);
				Vec2d subB; subB.load_a(subB_b);
				Vec2d subC; subC.load_a(subC_b);
				Vec2d subD; subD.load_a(subD_b);
				Vec4d subSample = inputSample[i] * subRate;
				Vec2d subSample2[2]; subSample2[0] = subSample.get_low(); subSample2[1] = subSample.get_high();
				for(int i = 0; i < 2; i++) {
					Vec2d scale = 0.5+abs(subSample2[i]*0.5);
					subSample2[i] = (subA+(sin(subA-subSample2[i])*scale));
					subA = subSample2[i]*scale;
					scale = 0.5+abs(subSample2[i]*0.5);
					subSample2[i] = (subB+(sin(subB-subSample2[i])*scale));
					subB = subSample2[i]*scale;
					scale = 0.5+abs(subSample2[i]*0.5);
					subSample2[i] = (subC+(sin(subC-subSample2[i])*scale));
					subC = subSample2[i]*scale;
					scale = 0.5+abs(subSample2[i]*0.5);
					subSample2[i] = (subD+(sin(subD-subSample2[i])*scale));
					subD = subSample2[i]*scale;
				}
				subSample = concatenate2(subSample2[0], subSample2[1]);
		
				subA.store_a(subA_b);
				subB.store_a(subB_b);
				subC.store_a(subC_b);
				subD.store_a(subD_b);
				subSample = min(max(subSample, -0.25), 0.25);
				inputSample[i] -= (subSample*16.0);
				//end SubTight section		
			}

			if (cycleEnd > 1) {
				Vec2d a = inputSample[i].get_low();
				Vec2d b = inputSample[i].get_high();
				double x[2] __attribute__((aligned(32))); b.store(x); b.load(tail); tail[0]=x[0];tail[1]=x[1];
				v = concatenate2(b, a);
				inputSample[i] = (inputSample[i] + v)*0.5;
			} //let's average only at elevated sample rates
		
			if (wet < 1.0) {inputSample[i] *= wet;}
			if (dry < 1.0) {drySample[i] *= dry;}
			inputSample[i] += drySample[i];
			//this is our submix verb dry/wet: 0.5 is BOTH at FULL VOLUME
			//purpose is that, if you're adding verb, you're not altering other balances
		
			//begin 32 bit stereo floating point dither
			fpd = fpd ^ (fpd << 13); fpd = fpd ^ (fpd << 17); fpd = fpd ^ (fpd << 5);
			Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample[i])) + 62);
			inputSample[i] += (to_double(fpd) - 2147483647.0) * to_double(exp) * 5.5e-36l;
			//end 32 bit stereo floating point dither

			double result[4] __attribute__((aligned(32)));
			inputSample[i].store_a(result);
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
	fpd.store_a(fpd_b);
}

void ClearCoatV::processDoubleReplacing(double **inputs, double **outputs, VstInt32 sampleFrames) 
{
    double* in1  =  inputs[0];
    double* in2  =  inputs[1];
    double* out1 = outputs[0];
    double* out2 = outputs[1];

}
