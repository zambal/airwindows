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
  kNumParameters = 16
}; //

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

  // Hypersonic

	enum {
		fix_freq,
		fix_reso,
		fix_a0,
		fix_a1,
		fix_a2,
		fix_b1,
		fix_b2,
		fix_sL1,
		fix_sL2,
		fix_sR1,
		fix_sR2,
		fix_total
	}; //fixed frequency biquad filter for ultrasonics, stereo

	double fixA[fix_total];
	double fixB[fix_total];
	double fixC[fix_total];
	double fixD[fix_total];
	double fixE[fix_total];

  // Interstage

	double iirSampleAL;
	double iirSampleBL;
	double iirSampleCL;
	double iirSampleDL;
	double iirSampleEL;
	double iirSampleFL;
	double iirSampleAR;
	double iirSampleBR;
	double iirSampleCR;
	double iirSampleDR;
	double iirSampleER;
	double iirSampleFR;

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

  // Tape

	double iirMidRollerAL;
	double iirMidRollerBL;
	double iirHeadBumpAL;
	double iirHeadBumpBL;

	double iirMidRollerAR;
	double iirMidRollerBR;
	double iirHeadBumpAR;
	double iirHeadBumpBR;

	double biquadAL[9];
	double biquadBL[9];
	double biquadCL[9];
	double biquadDL[9];

	double biquadAR[9];
	double biquadBR[9];
	double biquadCR[9];
	double biquadDR[9];
	bool flip;

	double lastSampleL;
	double lastSampleR;

  // Creature

	double slewL[102]; //probably worth just using a number here
	double slewR[102]; //probably worth just using a number here
  
  // ToneSlant

	double ts_bL[102];
	double ts_bR[102];
	double ts_f[102];		

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
};

#endif
