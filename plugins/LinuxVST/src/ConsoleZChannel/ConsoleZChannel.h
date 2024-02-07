/* ========================================
 *  ConsoleZChannel - ConsoleZChannel.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannel_H
#define __ConsoleZChannel_H

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
  kNumParameters = 17
}; //

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'czc';    //Change this to what the AU identity is!

class ConsoleZChannel :
    public AudioEffectX
{
public:
    ConsoleZChannel(audioMasterCallback audioMaster);
    ~ConsoleZChannel();
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

	uint32_t fpdL;
	uint32_t fpdR;
	//default stuff

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

	double fixF[fix_total];
	double fixG[fix_total];

  // Desk4

	double dL[10000];
	double controlL;
	double lastSampleL;
	double lastOutSampleL;
	double lastSlewL;

	double dR[10000];
	double controlR;
	double lastSampleR;
	double lastOutSampleR;
	double lastSlewR;

	int gcount;

  // Res2

	double mpkL[2005];
	double mpkR[2005];
	double f[66];
	double prevfreqMPeak;
	double prevamountMPeak;
	int mpc;

  // ConsoleLA

	double subAL;
	double subAR;
	double subBL;
	double subBR;
	double subCL;
	double subCR;
	double hullL[225];
	double hullR[225];
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
	double controlAposL;
	double controlAnegL;
	double controlBposL;
	double controlBnegL;
	double targetposL;
	double targetnegL;
	double avgLA;
	double avgLB;
	double nvgLA;
	double nvgLB;

	double controlAposR;
	double controlAnegR;
	double controlBposR;
	double controlBnegR;
	double targetposR;
	double targetnegR;
	double avgRA;
	double avgRB;
	double nvgRA;
	double nvgRB;
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
  float Q; //parameters. Always 0-1, and we scale/alter them elsewhere.

};

#endif
