/* ========================================
 *  ConsoleZBuss - ConsoleZBuss.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZBuss_H
#define __ConsoleZBuss_H

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
  kNumParameters = 18
}; //

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'czb';    //Change this to what the AU identity is!

class ConsoleZBuss :
    public AudioEffectX
{
public:
    ConsoleZBuss(audioMasterCallback audioMaster);
    ~ConsoleZBuss();
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

  // ConsoleLABuss
	double lastSinewL;
	double lastSinewR;

	double subAL;
	double subAR;
	double subBL;
	double subBR;
	double subCL;
	double subCR;

	double gainA;
	double gainB; //smoothed master fader for channel, from Z2 series filter code

  // Desk4
	double desk4_dL[10000];
	double desk4_controlL;
	double lastSampleL;
	double lastOutSampleL;
	double lastSlewL;

	double desk4_dR[10000];
	double desk4_controlR;
	double lastSampleR;
	double lastOutSampleR;
	double lastSlewR;

	int desk4_gcount;

  // Buttercomp2

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
	bool butter_flip;

  // VariMu

	double muVaryL;
	double muAttackL;
	double muNewSpeedL;
	double muSpeedAL;
	double muSpeedBL;
	double muCoefficientAL;
	double muCoefficientBL;
	double previousL;

	double muVaryR;
	double muAttackR;
	double muNewSpeedR;
	double muSpeedAR;
	double muSpeedBR;
	double muCoefficientAR;
	double muCoefficientBR;
	double previousR;
	bool flip;

  // BussColors4

	double bL[175]; //full buffer for high sample rates. Scales to 192K
	double bR[175]; //full buffer for high sample rates. Scales to 192K
	double dL[100]; //buffer for calculating sag as it relates to the dynamic impulse synthesis. To 192K.
	double dR[100]; //buffer for calculating sag as it relates to the dynamic impulse synthesis. To 192K.
	int c[35]; //just the number of taps we use, doesn't have to scale
	double g[9]; //console model
	double outg[9]; //console model
	double controlL;
	double controlR;
	double slowdynL;
	double slowdynR;
	int gcount;

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
};

#endif
