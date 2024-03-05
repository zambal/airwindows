/* ========================================
 *  Discontinuity - Discontinuity.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __Discontinuity_H
#include "Discontinuity.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new Discontinuity(audioMaster);}

Discontinuity::Discontinuity(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.5;

	for(int count = 0; count < predelay+2; count++) {
		dBaL[count] = 0.0;
		dBaR[count] = 0.0;
		dBbL[count] = 0.0;
		dBbR[count] = 0.0;
		dBcL[count] = 0.0;
		dBcR[count] = 0.0;
		dBdL[count] = 0.0;
		dBdR[count] = 0.0;
	}
	dBaX = 1;
	dBbX = 1;
	dBcX = 1;
	dBdX = 1;
	
	fpdL = 1.0; while (fpdL < 16386) fpdL = rand()*UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand()*UINT32_MAX;
	//this is reset: values being initialized only once. Startup values, whatever they are.
	
    _canDo.insert("plugAsChannelInsert"); // plug-in can be used as a channel insert effect.
    _canDo.insert("plugAsSend"); // plug-in can be used as a send effect.
    _canDo.insert("x2in2out"); 
    setNumInputs(kNumInputs);
    setNumOutputs(kNumOutputs);
    setUniqueID(kUniqueId);
    canProcessReplacing();     // supports output replacing
    canDoubleReplacing();      // supports double precision processing
	programsAreChunks(true);
    vst_strncpy (_programName, "Default", kVstMaxProgNameLen); // default program name
}

Discontinuity::~Discontinuity() {}
VstInt32 Discontinuity::getVendorVersion () {return 1000;}
void Discontinuity::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void Discontinuity::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 Discontinuity::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	chunkData[0] = A;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you 
	 started with. */
	
	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 Discontinuity::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void Discontinuity::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float Discontinuity::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void Discontinuity::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Top dB", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void Discontinuity::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (80.0+(A*40.0), text, kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void Discontinuity::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "dB", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 Discontinuity::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool Discontinuity::getEffectName(char* name) {
    vst_strncpy(name, "Discontinuity", kVstMaxProductStrLen); return true;
}

VstPlugCategory Discontinuity::getPlugCategory() {return kPlugCategEffect;}

bool Discontinuity::getProductString(char* text) {
  	vst_strncpy (text, "airwindows Discontinuity", kVstMaxProductStrLen); return true;
}

bool Discontinuity::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
