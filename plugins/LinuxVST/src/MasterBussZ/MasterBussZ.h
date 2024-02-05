/* ========================================
 *  MasterBussZ - MasterBussZ.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) 2011 __MyCompanyName__, Airwindows uses the MIT license
 * ======================================== */

#ifndef __MasterBussZ_H
#define __MasterBussZ_H

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
  kNumParameters = 5
}; //

enum {
    kDKAD,
    kDKCD,
    kPEAK,
    kSLEW,
    kSUBS,
    kMONO,
    kSIDE,
    kVINYL,
    kAURAT,
    kMONORAT,
    kMONOLAT,
    kPHONE,
    kCANSA,
    kCANSB,
    kCANSC,
    kCANSD,
    kTRICK,
};


const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'mbz';    //Change this to what the AU identity is!

class MasterBussZ :
    public AudioEffectX
{
public:
    MasterBussZ(audioMasterCallback audioMaster);
    ~MasterBussZ();
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

  // ADClip7

	double ad_lastSampleL;
	double ad_lastSampleR;
	float ad_bL[22200];
	float ad_bR[22200];
	int gcount;
	double lowsL;
	double lowsR;
	double iirLowsAL;
	double iirLowsAR;
	double iirLowsBL;
	double iirLowsBR;
	double refclipL;
	double refclipR;

  // Monitoring2

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
	};
	double biquad[fix_total];
	//Bandpasses

	float darkSampleL[100];
	float darkSampleR[100];

	double aL[1503], bL[1503], cL[1503], dL[1503];
	double aR[1503], bR[1503], cR[1503], dR[1503];
	int ax, bx, cx, dx;
	//PeaksOnly
	double lastSampleL, lastSampleR;
	//SlewOnly
	double iirSampleAL, iirSampleBL, iirSampleCL, iirSampleDL, iirSampleEL, iirSampleFL, iirSampleGL;
	double iirSampleHL, iirSampleIL, iirSampleJL, iirSampleKL, iirSampleLL, iirSampleML, iirSampleNL, iirSampleOL, iirSamplePL;
	double iirSampleQL, iirSampleRL, iirSampleSL;
	double iirSampleTL, iirSampleUL, iirSampleVL;
	double iirSampleWL, iirSampleXL, iirSampleYL, iirSampleZL;

	double iirSampleAR, iirSampleBR, iirSampleCR, iirSampleDR, iirSampleER, iirSampleFR, iirSampleGR;
	double iirSampleHR, iirSampleIR, iirSampleJR, iirSampleKR, iirSampleLR, iirSampleMR, iirSampleNR, iirSampleOR, iirSamplePR;
	double iirSampleQR, iirSampleRR, iirSampleSR;
	double iirSampleTR, iirSampleUR, iirSampleVR;
	double iirSampleWR, iirSampleXR, iirSampleYR, iirSampleZR; // o/`
	//SubsOnly

	uint32_t fpdL;
	uint32_t fpdR;
	//default stuff

  // Parameters

  float A;
  float B;
  float C;
  float D;
  float E;
};

#endif
