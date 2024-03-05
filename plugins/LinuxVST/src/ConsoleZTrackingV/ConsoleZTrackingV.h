/* ========================================
 *  ConsoleZTrackingV - ConsoleZTrackingV.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) 2011 __MyCompanyName__, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTrackingV_H
#define __ConsoleZTrackingV_H

#ifndef __audioeffect__
#include "audioeffectx.h"
#endif

#define MAX_VECTOR_SIZE 256

#include <set>
#include <string>
#include <math.h>
#include "vectorclass.h"

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


// Focus
const double centerFreq = 3515.775 / 4.0;

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'cztv';    //Change this to what the AU identity is!

class ConsoleZTrackingV :
    public AudioEffectX
{
public:
    ConsoleZTrackingV(audioMasterCallback audioMaster);
    ~ConsoleZTrackingV();
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

    double biquad_hpf[8] __attribute__((aligned(32)));

    // BiquadOneHalf LPF

  	double biquad_lpf[8] __attribute__((aligned(32)));

    // ButterComp2

  	double controlpos_b[4] __attribute__((aligned(32)));
  	double controlneg_b[4] __attribute__((aligned(32)));
  	double targetpos_b[4] __attribute__((aligned(32)));
  	double targetneg_b[4] __attribute__((aligned(32)));
  	double lastOutput_b[4] __attribute__((aligned(32)));

    // StereoChorus

  	int p_b[65536] __attribute__((aligned(32)));
  	double sweep_b[2] __attribute__((aligned(32)));
  	int	gcount;
  	double airPrev_b[2] __attribute__((aligned(32)));
  	double airEven_b[2] __attribute__((aligned(32)));
  	double airOdd_b[2] __attribute__((aligned(32)));
	
  	double lastRef[20] __attribute__((aligned(32)));
  	int chorus_cycle;
    bool chorus_flip;

    // ClearCoat

  	double aA[2048] __attribute__((aligned(32)));
  	double aB[4096] __attribute__((aligned(32)));
  	double aC[4096] __attribute__((aligned(32)));
  	double aD[8192] __attribute__((aligned(32)));
  	double aE[4096] __attribute__((aligned(32)));
  	double aF[4096] __attribute__((aligned(32)));
  	double aG[4096] __attribute__((aligned(32)));
  	double aH[4096] __attribute__((aligned(32)));
  	double aI[4096] __attribute__((aligned(32)));
  	double aJ[8192] __attribute__((aligned(32)));
  	double aK[8192] __attribute__((aligned(32)));
  	double aL[8192] __attribute__((aligned(32)));
  	double aM[2048] __attribute__((aligned(32)));
  	double aN[8192] __attribute__((aligned(32)));
  	double aO[2048] __attribute__((aligned(32)));
  	double aP[8192] __attribute__((aligned(32)));
	
	
  	double feedbackA[2] __attribute__((aligned(32)));
  	double feedbackB[2] __attribute__((aligned(32)));
  	double feedbackC[2] __attribute__((aligned(32)));
  	double feedbackD[2] __attribute__((aligned(32)));
	
  	int cc_cycle;
    int cc_gcount;
	
  	double cc_lastRef[8] __attribute__((aligned(32)));
  	double prevMulchA[2] __attribute__((aligned(32)));
  	double tail[2] __attribute__((aligned(32)));
	
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
	
  	int prevclearcoat;
	
  	double cc_sub_b[16] __attribute__((aligned(32)));

    // Channel9
  
  	double ch9_iirSamples[4] __attribute__((aligned(32)));
  	double ch9_lastSample[6] __attribute__((aligned(32)));
  	double ch9_biquad[16] __attribute__((aligned(32)));
  	double ch9_iirAmount;
  	double ch9_threshold;
  	double ch9_cutoff;

    // Focus

  	double figure[4] __attribute__((aligned(32)));

    double freqX;

    // Creature

    double cr_slew[256] __attribute__((aligned(32))); //probably worth just using a number here

    // ToneSlant

  	double ts_b[260] __attribute__((aligned(32)));
  	double ts_f[128] __attribute__((aligned(32)));		
    int ts_gcount;

    bool flip;
    double adcLastSample[2] __attribute__((aligned(32)));
  	uint32_t fpd_b[4] __attribute__((aligned(32)));
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
