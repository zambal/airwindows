/* ========================================
 *  Console8BussOut - Console8BussOut.h
 *  Created 8/12/11 by SPIAdmin 
 *  Copyright (c) 2011 __MyCompanyName__, All rights reserved
 * ======================================== */

#ifndef __Console8BussOut_H
#define __Console8BussOut_H

#ifndef __audioeffect__
#include "audioeffectx.h"
#endif

#include <set>
#include <string>
#include <math.h>

enum {
	kParamA = 0,
  kNumParameters = 1
}; //

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'c8bo';    //Change this to what the AU identity is!

class Console8BussOut : 
    public AudioEffectX 
{
public:
    Console8BussOut(audioMasterCallback audioMaster);
    ~Console8BussOut();
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
    
	double inTrimA;
	double inTrimB;
	bool hsr;
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
	double fix[fix_total];
	
	double lastSampleL;
	double intermediateL[18];
	bool wasPosClipL;
	bool wasNegClipL; //ClipOnly2
	double NSOddL;
	double NSEvenL;
	double prevShapeL;
	double darkSampleL[100];
	
	double lastSampleR;
	double intermediateR[18];
	bool wasPosClipR;
	bool wasNegClipR;
	double NSOddR;
	double NSEvenR;
	double prevShapeR; //Ten Nines
	double darkSampleR[100]; //Dark
	
	int spacing;
	bool flip;
	int depth;
	uint32_t fpdL;
	uint32_t fpdR;
	//default stuff
    float A;
};

#endif
