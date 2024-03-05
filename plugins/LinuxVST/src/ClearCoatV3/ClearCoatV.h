/* ========================================
 *  ClearCoatV - ClearCoatV.h
 *  Created 8/12/11 by SPIAdmin 
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ClearCoatV_H
#define __ClearCoatV_H

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
  kNumParameters = 2
}; //

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'clc3';    //Change this to what the AU identity is!

class ClearCoatV : 
    public AudioEffectX 
{
public:
    ClearCoatV(audioMasterCallback audioMaster);
    ~ClearCoatV();
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
	
	double lastRef[8] __attribute__((aligned(32)));

	int cycle;
  int gcount;	
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
		
    float A;
    float B;
};

#endif
