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
	 
	int shorts1[8] __attribute__((aligned(32))) = {shortA, shortB, shortC, shortD, shortD, shortH, shortL, shortP};
	int shorts2[8] __attribute__((aligned(32))) = {shortE, shortF, shortG, shortH, shortC, shortG, shortK, shortO};
	int shorts3[8] __attribute__((aligned(32))) = {shortI, shortJ, shortK, shortL, shortB, shortF, shortJ, shortN};
	int shorts4[8] __attribute__((aligned(32))) = {shortM, shortN, shortO, shortP, shortA, shortE, shortI, shortM};
	int kshorts1[8] __attribute__((aligned(32))) = {kshortA, kshortB, kshortC, kshortD, kshortD, kshortH, kshortL, kshortP};
	int kshorts2[8] __attribute__((aligned(32))) = {kshortE, kshortF, kshortG, kshortH, kshortC, kshortG, kshortK, kshortO};
	int kshorts3[8] __attribute__((aligned(32))) = {kshortI, kshortJ, kshortK, kshortL, kshortB, kshortF, kshortJ, kshortN};
	int kshorts4[8] __attribute__((aligned(32))) = {kshortM, kshortN, kshortO, kshortP, kshortA, kshortE, kshortI, kshortM};
	int load_c[8] __attribute__((aligned(32)));

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
	Vec4d inputSample[2], v, fpd_d = to_double(fpd);
	Vec4d drySample[2];

  while (sampleFrames > 0)
  {
		Vec4d l(*in1, *(in1 + 1), *(in1 + 2), *(in1 + 3));
		l = select(abs(l) < 1.18e-23, fpd_d * 1.18e-17, l); 
		Vec4d r(*in2, *(in2 + 1), *(in2 + 2), *(in2 + 3));
		r = select(abs(r) < 1.18e-23, fpd_d * 1.18e-17, r); 
		drySample[0] = blend4<0,4,1,5>(l, r);
		drySample[1] = blend4<2,6,3,7>(l, r);
		cycle+=4;
		if (cycle == (cycleEnd * 4)) { //hit the end point and we do a reverb sample
			Vec4d o1, o2, o3, o4;
			Vec8i c, s, ks, t;
			
			v.load_a(feedbackAL); v = (l + (v * 0.04166666666)); v.store_a(aAL + counters1[0]);
			if(counters1[0] == 0) v.store_a(aEL + kshorts1[0]);
			v.load_a(feedbackBL); v = (l + (v * 0.04166666666)); v.store_a(aBL + counters1[1]);
			if(counters1[1] == 0) v.store_a(aFL + kshorts1[1]);
			v.load_a(feedbackCL); v = (l + (v * 0.04166666666)); v.store_a(aCL + counters1[2]);
			if(counters1[2] == 0) v.store_a(aGL + kshorts1[2]);
			v.load_a(feedbackDL); v = (l + (v * 0.04166666666)); v.store_a(aDL + counters1[3]);
			if(counters1[3] == 0) v.store_a(aHL + kshorts1[3]);

			v.load_a(feedbackDR); v = (r + (v * 0.04166666666)); v.store_a(aDR + counters1[4]);
			if(counters1[4] == 0) v.store_a(aCR + kshorts1[4]);
			v.load_a(feedbackHR); v = (r + (v * 0.04166666666)); v.store_a(aHR + counters1[5]);
			if(counters1[5] == 0) v.store_a(aGR + kshorts1[5]);
			v.load_a(feedbackLR); v = (r + (v * 0.04166666666)); v.store_a(aLR + counters1[6]);
			if(counters1[6] == 0) v.store_a(aKR + kshorts1[6]);
			v.load_a(feedbackPR); v = (r + (v * 0.04166666666)); v.store_a(aPR + counters1[7]);
			if(counters1[7] == 0) v.store_a(aOR + kshorts1[7]);
			
			c.load_a(counters1);
			s.load_a(shorts1);
			ks.load_a(kshorts1);
			t = c - s; if_add(t < 0, t, ks).store_a(load_c);
			t = c + 4; select(t > ks, 0, t).store_a(counters1);

			o1.load(aAL + load_c[0]);
			o2.load(aBL + load_c[1]);
			o3.load(aCL + load_c[2]);
			o4.load(aDL + load_c[3]);
	
			v = (o1 - (o2 + o3 + o4)); v.store_a(aEL + counters2[0]);
			if(counters2[0] == 0) v.store_a(aEL + kshorts2[0]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aFL + counters2[1]);
			if(counters2[1] == 0) v.store_a(aFL + kshorts2[1]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aGL + counters2[2]);
			if(counters2[2] == 0) v.store_a(aGL + kshorts2[2]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aHL + counters2[3]);
			if(counters2[3] == 0) v.store_a(aHL + kshorts2[3]);

			o1.load(aDR + load_c[4]);
			o2.load(aHR + load_c[5]);
			o3.load(aLR + load_c[6]);
			o4.load(aPR + load_c[7]);

			v = (o1 - (o2 + o3 + o4)); v.store_a(aCR + counters2[4]);
			if(counters2[4] == 0) v.store_a(aCR + kshorts2[4]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aGR + counters2[5]);
			if(counters2[5] == 0) v.store_a(aGR + kshorts2[5]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aKR + counters2[6]);
			if(counters2[6] == 0) v.store_a(aKR + kshorts2[6]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aOR + counters2[7]);
			if(counters2[7] == 0) v.store_a(aOR + kshorts2[7]);

			c.load_a(counters2);
			s.load_a(shorts2);
			ks.load_a(kshorts2);
			t = c - s; if_add(t < 0, t, ks).store_a(load_c);
			t = c + 4; select(t > ks, 0, t).store_a(counters2);
		
			o1.load(aEL + load_c[0]);
			o2.load(aFL + load_c[1]);
			o3.load(aGL + load_c[2]);
			o4.load(aHL + load_c[3]);
	
			v = (o1 - (o2 + o3 + o4)); v.store_a(aIL + counters3[0]);
			if(counters3[0] == 0) v.store_a(aIL + kshorts3[0]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aJL + counters3[1]);
			if(counters3[1] == 0) v.store_a(aJL + kshorts3[1]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aKL + counters3[2]);
			if(counters3[2] == 0) v.store_a(aKL + kshorts3[2]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aLL + counters3[3]);
			if(counters3[3] == 0) v.store_a(aLL + kshorts3[3]);
		
			o1.load(aCR + load_c[4]);
			o2.load(aGR + load_c[5]);
			o3.load(aKR + load_c[6]);
			o4.load(aOR + load_c[7]);

			v = (o1 - (o2 + o3 + o4)); v.store_a(aBR + counters3[4]);
			if(counters3[4] == 0) v.store_a(aBR + kshorts3[4]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aFR + counters3[5]);
			if(counters3[5] == 0) v.store_a(aFR + kshorts3[5]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aJR + counters3[6]);
			if(counters3[6] == 0) v.store_a(aJR + kshorts3[6]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aNR + counters3[7]);
			if(counters3[7] == 0) v.store_a(aNR + kshorts3[7]);
		
			c.load_a(counters3);
			s.load_a(shorts3);
			ks.load_a(kshorts3);
			t = c - s; if_add(t < 0, t, ks).store_a(load_c);
			t = c + 4; select(t > ks, 0, t).store_a(counters3);

			o1.load(aIL + load_c[0]);
			o2.load(aJL + load_c[1]);
			o3.load(aKL + load_c[2]);
			o4.load(aLL + load_c[3]);
		
			v = (o1 - (o2 + o3 + o4)); v.store_a(aML + counters4[0]);
			if(counters4[0] == 0) v.store_a(aML + kshorts4[0]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aNL + counters4[1]);
			if(counters4[1] == 0) v.store_a(aNL + kshorts4[1]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aOL + counters4[2]);
			if(counters4[2] == 0) v.store_a(aOL + kshorts4[2]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aPL + counters4[3]);
			if(counters4[3] == 0) v.store_a(aPL + kshorts4[3]);
		
			o1.load(aBR + load_c[4]);
			o2.load(aFR + load_c[5]);
			o3.load(aJR + load_c[6]);
			o4.load(aNR + load_c[7]);
		
			v = (o1 - (o2 + o3 + o4)); v.store_a(aAR + counters4[4]);
			if(counters4[4] == 0) v.store_a(aAR + kshorts4[4]);
			v = (o2 - (o1 + o3 + o4)); v.store_a(aER + counters4[5]);
			if(counters4[5] == 0) v.store_a(aER + kshorts4[5]);
			v = (o3 - (o1 + o2 + o4)); v.store_a(aIR + counters4[6]);
			if(counters4[6] == 0) v.store_a(aIR + kshorts4[6]);
			v = (o4 - (o1 + o2 + o3)); v.store_a(aMR + counters4[7]);
			if(counters4[7] == 0) v.store_a(aMR + kshorts4[7]);
		
			c.load_a(counters4);
			s.load_a(shorts4);
			ks.load_a(kshorts4);
			t = c - s; if_add(t < 0, t, ks).store_a(load_c);
			t = c + 4; select(t > ks, 0, t).store_a(counters4);

			o1.load(aML + load_c[0]);
			o2.load(aNL + load_c[1]);
			o3.load(aOL + load_c[2]);
			o4.load(aPL + load_c[3]);
		
			v = permute4<3, 0, 1, 2>(o1);
			double x = v.extract(0);
			v.insert(0, prevMulchAL);
			prevMulchAL = x;
			o1 = (o1 + o1 + o1 + v)*0.25;

			(o1 - (o2 + o3 + o4)).store_a(feedbackAL);
			(o2 - (o1 + o3 + o4)).store_a(feedbackBL);
			(o3 - (o1 + o2 + o4)).store_a(feedbackCL);
			(o4 - (o1 + o2 + o3)).store_a(feedbackDL);

			l = min(max((o1 + o2 + o3 + o4)*0.125, -1.0), 1.0);

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
		
			r = min(max((o1 + o2 + o3 + o4)*0.125, -1.0), 1.0);
			//and take the final combined sum of outputs, corrected for Householder gain

			if (cycleEnd == 4) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[32]);
				lastRef[32] = x;
				scatter<0, 8, 16, 24>(v, lastRef);
				scatter<2, 10, 18, 26>((v + v + v + l) * 0.25, lastRef);
				scatter<4, 12, 20, 28>((v + v + l + l) * 0.25, lastRef);
				scatter<6, 14, 22, 30>((v + l + l + l) * 0.25, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[33]);
				lastRef[33] = x;
				scatter<1, 9, 17, 25>(v, lastRef);
				scatter<3, 11, 19, 27>((v + v + v + r) * 0.25, lastRef);
				scatter<5, 13, 21, 29>((v + v + r + r) * 0.25, lastRef);
				scatter<7, 15, 23, 31>((v + r + r + r) * 0.25, lastRef);
			}
			if (cycleEnd == 3) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[24]);
				lastRef[24] = x;
				scatter<0, 6, 12, 18>(v, lastRef);
				scatter<2, 8, 14, 20>((v + v + l) * 0.333, lastRef);
				scatter<4, 10, 16, 22>((v + l + l) * 0.333, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[25]);
				lastRef[25] = x;
				scatter<1, 7, 13, 19>(v, lastRef);
				scatter<3, 9, 15, 21>((v + v + r) * 0.333, lastRef);
				scatter<5, 11, 17, 23>((v + r + r) * 0.333, lastRef);
			}
			if (cycleEnd == 2) {
				v = permute4<3, 0, 1, 2>(l);
				double x = v.extract(0);
				v.insert(0, lastRef[16]);
				lastRef[16] = x;
				scatter<0, 4, 8, 12>(v, lastRef);
				scatter<2, 6, 10, 14>((v + l) * 0.5, lastRef);

				v = permute4<3, 0, 1, 2>(r);
				x = v.extract(0);
				v.insert(0, lastRef[17]);
				lastRef[17] = x;
				scatter<1, 5, 9, 13>(v, lastRef);
				scatter<3, 7, 11, 15>((v + r) * 0.5, lastRef);
			}
			if (cycleEnd == 1) {
				scatter<0, 2, 4, 6>(l, lastRef);
				scatter<1, 3, 5, 7>(r, lastRef);
			}
			cycle = 0; //reset
			inputSample[0].load_a(lastRef);
			inputSample[1].load_a(lastRef + 4);
		} else {
			inputSample[0].load_a(lastRef + cycle);
			inputSample[1].load_a(lastRef + cycle + 4);
			//we are going through our references now
		}
	
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
				double x[2] __attribute__((aligned(32))); b.store_a(x); b.load_a(tail); tail[0]=x[0];tail[1]=x[1];
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
			fpd_d = to_double(fpd);
			Vec4q exp = Vec4q(1) << extend(exponent(to_float(inputSample[i])) + 62);
			inputSample[i] += (fpd_d - 2147483647.0) * to_double(exp) * 5.5e-36l;
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
