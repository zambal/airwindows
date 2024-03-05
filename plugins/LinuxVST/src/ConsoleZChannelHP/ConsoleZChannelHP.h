/* ========================================
 *  ConsoleZChannelHP - ConsoleZChannelHP.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannelHP_H
#define __ConsoleZChannelHP_H

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

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'czc';    //Change this to what the AU identity is!

class ConsoleZChannelHP :
    public AudioEffectX
{
public:
    ConsoleZChannelHP(audioMasterCallback audioMaster);
    ~ConsoleZChannelHP();
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

	uint32_t fpd[2];
	//default stuff

  // Hypersonic

	enum {
		fix_a0,
		fix_a1,
		fix_a2,
		fix_b1,
		fix_b2,
		fix_sL1,
		fix_sR1,
		fix_sL2,
		fix_sR2,
		fix_total
	}; //fixed frequency biquad filter for ultrasonics, stereo

	double fixA[fix_total];
	double fixD[fix_total];
	double fixE[fix_total];
	double fixG[fix_total];

  // BiquadOneHalf 1

	double biquad_1[13];

  // BiquadOneHalf 2

	double biquad_2[13];

  // Desk4

	double d[10000][2];
	double control[2];
	double lastSample[2];
	double lastOutSample[2];
	double lastSlew[2];

	int gcount;

  // PurestDrive
	double pd_previousSample[2];

  // Res2

	double mpk[2005][2];
	double f[66];
	double prevfreqMPeak;
	double prevamountMPeak;
	int mpc;

  // ConsoleLA

	double subA[2];
	double subB[2];
	double subC[2];
	double hull[225][2];
	int hullp;
	double pearB[22];
	double midA;
	double midB;
	double bassA;
	double bassB;
	double gainA;
	double gainB; //smoothed master fader for channel, from Z2 series filter code

  // CStrip2 comp

	//begin ButterComp
	double controlApos[2];
	double controlAneg[2];
	double controlBpos[2];
	double controlBneg[2];
	double targetpos[2];
	double targetneg[2];
	double avgA[2];
	double avgB[2];
	double nvgA[2];
	double nvgB[2];

	bool flip;
	//end ButterComp

  // Wider

	double p[4099];
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
