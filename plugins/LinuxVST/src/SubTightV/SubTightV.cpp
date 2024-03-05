/* ========================================
 *  SubTightV - SubTightV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __SubTightV_H
#include "SubTightV.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new SubTightV(audioMaster);}

SubTightV::SubTightV(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.3;
	B = 0.5;
	for (int x = 0; x < 32; x++) {
		sub[x] = 0.0;
	}
	
	fpd_b[0] = 1.0; while (fpd_b[0] < 16386) fpd_b[0] = rand()*UINT32_MAX;
	fpd_b[1] = 1.0; while (fpd_b[1] < 16386) fpd_b[1]= rand()*UINT32_MAX;
  fpd_b[2] = fpd_b[0]; fpd_b[2] ^= fpd_b[2] << 13; fpd_b[2] ^= fpd_b[2] >> 17; fpd_b[2] ^= fpd_b[2] << 5; 
  fpd_b[3] = fpd_b[1]; fpd_b[3] ^= fpd_b[3] << 13; fpd_b[3] ^= fpd_b[3] >> 17; fpd_b[3] ^= fpd_b[3] << 5;
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

SubTightV::~SubTightV() {}
VstInt32 SubTightV::getVendorVersion () {return 1000;}
void SubTightV::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void SubTightV::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 SubTightV::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	chunkData[0] = A;
	chunkData[1] = B;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you 
	 started with. */
	
	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 SubTightV::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	B = pinParameter(chunkData[1]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void SubTightV::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        case kParamB: B = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float SubTightV::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        case kParamB: return B; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void SubTightV::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Trim", kVstMaxParamStrLen); break;
		case kParamB: vst_strncpy (text, "Steep", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void SubTightV::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (A, text, kVstMaxParamStrLen); break;
        case kParamB: float2string (B, text, kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void SubTightV::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 SubTightV::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool SubTightV::getEffectName(char* name) {
    vst_strncpy(name, "SubTightV", kVstMaxProductStrLen); return true;
}

VstPlugCategory SubTightV::getPlugCategory() {return kPlugCategEffect;}

bool SubTightV::getProductString(char* text) {
  	vst_strncpy (text, "airwindows SubTightV", kVstMaxProductStrLen); return true;
}

bool SubTightV::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
