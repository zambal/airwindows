/* ========================================
 *  ConsoleZTracking - ConsoleZTracking.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) 2011 __MyCompanyName__, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTracking_H
#define __ConsoleZTracking_H

#ifndef __audioeffect__
#include "audioeffectx.h"
#endif

#include <set>
#include <string>
#include <math.h>

enum {
	kParamA = 0,
	kParamB = 1,
	kParamC = 2,
	kParamD = 3,
	kParamE = 4,
	kParamF = 5,
	kParamG = 6,
	kParamH = 7,
	kParamI = 8,
	kParamJ = 9,
	kParamK = 10,
	kParamL = 11,
	kParamM = 12,
	kParamN = 13,
	kParamO = 14,
	kParamP = 15,
	kParamQ = 16,
	kParamR = 17,
	kParamS = 18,
	kParamT = 19,
	kParamU = 20,
	kParamV = 21,
	kParamW = 22,
	kParamX = 23,
	kParamY = 24,
	kParamZ = 25,
  kNumParameters = 26
}; //

// ClearCoat
const int kshortA = 350;
const int kshortB = 1710;
const int kshortC = 1610;
const int kshortD = 835;
const int kshortE = 700;
const int kshortF = 1260;
const int kshortG = 1110;
const int kshortH = 1768;
const int kshortI = 280;
const int kshortJ = 2645;
const int kshortK = 1410;
const int kshortL = 1175;
const int kshortM = 12;
const int kshortN = 3110;
const int kshortO = 120;
const int kshortP = 2370;

// Focus
const double centerFreq = 3515.775 / 4.0;

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'cztr';    //Change this to what the AU identity is!

class ConsoleZTracking :
    public AudioEffectX
{
public:
    ConsoleZTracking(audioMasterCallback audioMaster);
    ~ConsoleZTracking();
    virtual bool getEffectName(char* name);                       // The plug-in name
    virtual VstPlugCategory getPlugCategory();                    // The general category for the plug-in
    virtual bool getProductString(char* text);                    // This is a unique plug-in string provided by Steinberg
    virtual bool getVendorString(char* text);                     // Vendor info
    virtual VstInt32 getVendorVersion();                          // Version number
    virtual void processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames);
    virtual void processDoubleReplacing (double** inputs, double** outputs, VstInt32 sampleFrames);
    virtual void getProgramName(char *name);                      // read the name from the host
    virtual void setProgramName(char *name);                      // changes the name of the preset displayed in the host
  	virtual VstInt32 getChunk (void** data, bool isPreset);
  	virtual VstInt32 setChunk (void* data, VstInt32 byteSize, bool isPreset);
    virtual float getParameter(VstInt32 index);                   // get the parameter value at the specified index
    virtual void setParameter(VstInt32 index, float value);       // set the parameter at index to value
    virtual void getParameterLabel(VstInt32 index, char *text);  // label for the parameter (eg dB)
    virtual void getParameterName(VstInt32 index, char *text);    // name of the parameter
    virtual void getParameterDisplay(VstInt32 index, char *text); // text description of the current value
    virtual VstInt32 canDo(char *text);
private:
    char _programName[kVstMaxProgNameLen + 1];
    std::set< std::string > _canDo;


  	// BitshiftGain

  	double input_gain;

    // BiquadOneHalf HPF

  	double biquad_hpf_AL[9];
  	double biquad_hpf_AR[9];
  	double biquad_hpf_BL[9];
  	double biquad_hpf_BR[9];

    // BiquadOneHalf LPF

  	double biquad_lpf_AL[9];
  	double biquad_lpf_AR[9];
  	double biquad_lpf_BL[9];
  	double biquad_lpf_BR[9];

    // ButterComp2

  	double controlAposL;
  	double controlAnegL;
  	double controlBposL;
  	double controlBnegL;
  	double targetposL;
  	double targetnegL;
  	double lastOutputL;
  	double controlAposR;
  	double controlAnegR;
  	double controlBposR;
  	double controlBnegR;
  	double targetposR;
  	double targetnegR;
  	double lastOutputR;

    // StereoChorus

  	int pL[65536];
  	int pR[65536];
  	double sweepL;
  	double sweepR;
  	int	gcount;
  	double airPrevL;
  	double airEvenL;
  	double airOddL;
  	double airFactorL;
  	double airPrevR;
  	double airEvenR;
  	double airOddR;
  	double airFactorR;
	
  	double lastRefL[7];
  	double lastRefR[7];
  	int chorus_cycle;
    bool chorus_flip;

    // ClearCoat

  	double aAL[kshortA+5];
  	double aBL[kshortB+5];
  	double aCL[kshortC+5];
  	double aDL[kshortD+5];
  	double aEL[kshortE+5];
  	double aFL[kshortF+5];
  	double aGL[kshortG+5];
  	double aHL[kshortH+5];
  	double aIL[kshortI+5];
  	double aJL[kshortJ+5];
  	double aKL[kshortK+5];
  	double aLL[kshortL+5];
  	double aML[kshortM+5];
  	double aNL[kshortN+5];
  	double aOL[kshortO+5];
  	double aPL[kshortP+5];
	
  	double aAR[kshortA+5];
  	double aBR[kshortB+5];
  	double aCR[kshortC+5];
  	double aDR[kshortD+5];
  	double aER[kshortE+5];
  	double aFR[kshortF+5];
  	double aGR[kshortG+5];
  	double aHR[kshortH+5];
  	double aIR[kshortI+5];
  	double aJR[kshortJ+5];
  	double aKR[kshortK+5];
  	double aLR[kshortL+5];
  	double aMR[kshortM+5];
  	double aNR[kshortN+5];
  	double aOR[kshortO+5];
  	double aPR[kshortP+5];
	
  	double feedbackAL;
  	double feedbackBL;
  	double feedbackCL;
  	double feedbackDL;
	
  	double feedbackDR;
  	double feedbackHR;
  	double feedbackLR;
  	double feedbackPR;
	
  	double previousAL;
  	double previousBL;
  	double previousCL;
  	double previousDL;
  	double previousEL;
	
  	double cc_LastRefL[7];
	
  	double previousAR;
  	double previousBR;
  	double previousCR;
  	double previousDR;
  	double previousER;
	
  	double cc_LastRefR[7];
	
  	int countAL;
  	int countBL;
  	int countCL;
  	int countDL;
  	int countEL;
  	int countFL;
  	int countGL;
  	int countHL;
  	int countIL;
  	int countJL;
  	int countKL;
  	int countLL;		
  	int countML;		
  	int countNL;		
  	int countOL;		
  	int countPL;		
	
  	int countAR;
  	int countBR;
  	int countCR;
  	int countDR;
  	int countER;
  	int countFR;
  	int countGR;
  	int countHR;
  	int countIR;
  	int countJR;
  	int countKR;
  	int countLR;		
  	int countMR;		
  	int countNR;		
  	int countOR;		
  	int countPR;	
	
  	int cc_cycle;
	
  	double prevMulchAL;
  	double prevMulchAR;
	
  	double tailL;
  	double tailR;
	
  	int shortA;
  	int shortB;
  	int shortC;
  	int shortD;
  	int shortE;
  	int shortF;
  	int shortG;
  	int shortH;
  	int shortI;
  	int shortJ;
  	int shortK;
  	int shortL;
  	int shortM;
  	int shortN;
  	int shortO;
  	int shortP;
	
  	int clearcoat;
  	int prevclearcoat;
	
  	double subAL;
  	double subAR;
  	double subBL;
  	double subBR;
  	double subCL;
  	double subCR;
  	double subDL;
  	double subDR;	

    // Channel9
  
  	double ch9_iirSampleLA;
  	double ch9_iirSampleRA;
  	double ch9_iirSampleLB;
  	double ch9_iirSampleRB;
  	double ch9_lastSampleAL;
  	double ch9_lastSampleBL;
  	double ch9_lastSampleCL;
  	double ch9_lastSampleAR;
  	double ch9_lastSampleBR;
  	double ch9_lastSampleCR;
  	double ch9_biquadA[15];
  	double ch9_biquadB[15];
  	double ch9_iirAmount;
  	double ch9_threshold;
  	double ch9_cutoff;

    // Focus

  	double figureL[9];
  	double figureR[9];

    double freqX;

    // Creature

  	double slewL[102]; //probably worth just using a number here
  	double slewR[102]; //probably worth just using a number here

    // ToneSlant

  	double ts_bL[102];
  	double ts_bR[102];
  	double ts_f[102];

    bool flip;
    double adcLastSampleL;
    double adcLastSampleR;
  	uint32_t fpdL;
  	uint32_t fpdR;
  	//default stuff

    // Parameters

    float A;
    float B;
    float C;
    float D;
    float E;
    float F;
    float G;
    float H;
    float I;
    float J;
    float K;
    float L;
    float M;
    float N;
    float O;
    float P;
    float Q;
    float R;
    float S;
    float T;
    float U;
    float V;
    float W;
    float X;
    float Y;
    float Z;
};

#endif
