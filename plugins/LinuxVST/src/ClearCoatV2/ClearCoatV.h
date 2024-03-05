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

const int kshortA = 352;
const int kshortB = 1712;
const int kshortC = 1612;
const int kshortD = 836;
const int kshortE = 700;
const int kshortF = 1260;
const int kshortG = 1112;
const int kshortH = 1768;
const int kshortI = 280;
const int kshortJ = 2648;
const int kshortK = 1412;
const int kshortL = 1176;
const int kshortM = 12;
const int kshortN = 3112;
const int kshortO = 120;
const int kshortP = 2372;

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'clcv';    //Change this to what the AU identity is!

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

	double aAL[kshortA+4] __attribute__((aligned(32)));
	double aBL[kshortB+4] __attribute__((aligned(32)));
	double aCL[kshortC+4] __attribute__((aligned(32)));
	double aDL[kshortD+4] __attribute__((aligned(32)));
	double aEL[kshortE+4] __attribute__((aligned(32)));
	double aFL[kshortF+4] __attribute__((aligned(32)));
	double aGL[kshortG+4] __attribute__((aligned(32)));
	double aHL[kshortH+4] __attribute__((aligned(32)));
	double aIL[kshortI+4] __attribute__((aligned(32)));
	double aJL[kshortJ+4] __attribute__((aligned(32)));
	double aKL[kshortK+4] __attribute__((aligned(32)));
	double aLL[kshortL+4] __attribute__((aligned(32)));
	double aML[kshortM+4] __attribute__((aligned(32)));
	double aNL[kshortN+4] __attribute__((aligned(32)));
	double aOL[kshortO+4] __attribute__((aligned(32)));
	double aPL[kshortP+4] __attribute__((aligned(32)));
	
	double aAR[kshortA+4] __attribute__((aligned(32)));
	double aBR[kshortB+4] __attribute__((aligned(32)));
	double aCR[kshortC+4] __attribute__((aligned(32)));
	double aDR[kshortD+4] __attribute__((aligned(32)));
	double aER[kshortE+4] __attribute__((aligned(32)));
	double aFR[kshortF+4] __attribute__((aligned(32)));
	double aGR[kshortG+4] __attribute__((aligned(32)));
	double aHR[kshortH+4] __attribute__((aligned(32)));
	double aIR[kshortI+4] __attribute__((aligned(32)));
	double aJR[kshortJ+4] __attribute__((aligned(32)));
	double aKR[kshortK+4] __attribute__((aligned(32)));
	double aLR[kshortL+4] __attribute__((aligned(32)));
	double aMR[kshortM+4] __attribute__((aligned(32)));
	double aNR[kshortN+4] __attribute__((aligned(32)));
	double aOR[kshortO+4] __attribute__((aligned(32)));
	double aPR[kshortP+4] __attribute__((aligned(32)));
	
	double feedbackAL[4] __attribute__((aligned(32)));
	double feedbackBL[4] __attribute__((aligned(32)));
	double feedbackCL[4] __attribute__((aligned(32)));
	double feedbackDL[4] __attribute__((aligned(32)));
	double feedbackDR[4] __attribute__((aligned(32)));
	double feedbackHR[4] __attribute__((aligned(32)));
	double feedbackLR[4] __attribute__((aligned(32)));
	double feedbackPR[4] __attribute__((aligned(32)));
	
	double lastRef[34] __attribute__((aligned(32)));

  int counters1[8] __attribute__((aligned(32)));
  int counters2[8] __attribute__((aligned(32)));
  int counters3[8] __attribute__((aligned(32)));
  int counters4[8] __attribute__((aligned(32)));
	
	int cycle;
	
	double prevMulchAL;
	double prevMulchAR;
	
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
	
	double subA_b[2] __attribute__((aligned(32)));
	double subB_b[2] __attribute__((aligned(32)));
	double subC_b[2] __attribute__((aligned(32)));
	double subD_b[2] __attribute__((aligned(32)));
		
    float A;
    float B;
};

#endif
