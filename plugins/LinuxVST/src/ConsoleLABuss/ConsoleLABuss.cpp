/* ========================================
 *  ConsoleLABuss - ConsoleLABuss.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleLABuss_H
#include "ConsoleLABuss.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ConsoleLABuss(audioMaster);}

ConsoleLABuss::ConsoleLABuss(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 1.0;

	lastSinewL = lastSinewR = 0.0;
	subAL = subAR = subBL = subBR = subCL = subCR = 0.0;
	gainA = gainB = 1.0;

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

ConsoleLABuss::~ConsoleLABuss() {}
VstInt32 ConsoleLABuss::getVendorVersion () {return 1000;}
void ConsoleLABuss::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ConsoleLABuss::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ConsoleLABuss::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	chunkData[0] = A;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you 
	 started with. */
	
	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 ConsoleLABuss::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void ConsoleLABuss::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float ConsoleLABuss::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ConsoleLABuss::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Master", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ConsoleLABuss::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (A, text, kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ConsoleLABuss::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ConsoleLABuss::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ConsoleLABuss::getEffectName(char* name) {
    vst_strncpy(name, "ConsoleLABuss", kVstMaxProductStrLen); return true;
}

VstPlugCategory ConsoleLABuss::getPlugCategory() {return kPlugCategEffect;}

bool ConsoleLABuss::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ConsoleLABuss", kVstMaxProductStrLen); return true;
}

bool ConsoleLABuss::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
