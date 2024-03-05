/* ========================================
 *  CSCompV - CSCompV.h
 *  Created 8/12/11 by SPIAdmin 
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __CSCompV_H
#define __CSCompV_H

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
	kNumParameters = 3
}; //

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'cscv';    //Change this to what the AU identity is!

class CSCompV : 
    public AudioEffectX 
{
public:
    CSCompV(audioMasterCallback audioMaster);
    ~CSCompV();
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
	
	//begin ButterComp
	double cs_controlpos[4] __attribute__((aligned(32)));
	double cs_controlneg[4] __attribute__((aligned(32)));
	double cs_targetpos[4] __attribute__((aligned(32)));
	double cs_targetneg[4] __attribute__((aligned(32)));
	double cs_avg[2] __attribute__((aligned(32)));
	double cs_nvg[2] __attribute__((aligned(32)));
	//end ButterComp
	//flip is already covered in EQ
	
    float A;
    float B;
    float C;
	
};

#endif
