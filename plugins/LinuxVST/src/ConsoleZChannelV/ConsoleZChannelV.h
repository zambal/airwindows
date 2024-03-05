/* ========================================
 *  ConsoleZChannelV - ConsoleZChannelV.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannelV_H
#define __ConsoleZChannelV_H

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

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'czv';    //Change this to what the AU identity is!

class ConsoleZChannelV :
    public AudioEffectX
{
public:
    ConsoleZChannelV(audioMasterCallback audioMaster);
    ~ConsoleZChannelV();
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

	uint32_t fpd_b[4] __attribute__((aligned(32)));
	//default stuff

  // Hypersonic

	enum {
		fix_sL1,
		fix_sR1,
		fix_sL2,
		fix_sr2,
		fix_freq,
		fix_reso,
		fix_a0,
		fix_a1,
		fix_a2,
		fix_b1,
		fix_b2,
		fix_total
	}; //fixed frequency biquad filter for ultrasonics, stereo

	double fixA[fix_total];
	double fixD[fix_total];
	double fixE[fix_total];
	double fixG[fix_total];

  // BiquadOneHalf 1

  double biquad_1[8] __attribute__((aligned(32)));

  // BiquadOneHalf 2

	double biquad_2[8] __attribute__((aligned(32)));

  // Desk4

	double d_b[20000] __attribute__((aligned(32)));
	double control_b[2] __attribute__((aligned(32)));
	double lastSample_b[2] __attribute__((aligned(32)));
	double lastOutSample_b[2] __attribute__((aligned(32)));
	double lastSlew_b[2] __attribute__((aligned(32)));

	int gcount;

  // PurestDrive
	double pd_previousSample[4] __attribute__((aligned(32)));

  // ResEQ2

	double mpk[1024] __attribute__((aligned(32)));
	double f[64] __attribute__((aligned(32)));
	double prevfreqMPeak;
	double prevamountMPeak;
	int mpc;

  // ConsoleLA

	double subA_b[2] __attribute__((aligned(32)));
	double subB_b[2] __attribute__((aligned(32)));
	double subC_b[2] __attribute__((aligned(32)));
	double hull[260] __attribute__((aligned(32)));	
	int hullp;
	double pearB[20] __attribute__((aligned(32)));
  double treble;
	double midA;
	double midB;
	double bassA;
	double bassB;
	double gainA;
	double gainB; //smoothed master fader for channel, from Z2 series filter code	
	double gainL;
	double gainR;

  // CStrip2 comp

	//begin ButterComp
	double cs_controlpos[4] __attribute__((aligned(32)));
	double cs_controlneg[4] __attribute__((aligned(32)));
	double cs_targetpos[4] __attribute__((aligned(32)));
	double cs_targetneg[4] __attribute__((aligned(32)));
	double cs_avg[2] __attribute__((aligned(32)));
	double cs_nvg[2] __attribute__((aligned(32)));
	//end ButterComp

  // Wider

	double p[2050] __attribute__((aligned(32)));
	int count;
  

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
