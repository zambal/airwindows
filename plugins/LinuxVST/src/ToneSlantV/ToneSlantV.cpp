/* ========================================
 *  ToneSlantV - ToneSlantV.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ToneSlantV_H
#include "ToneSlantV.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ToneSlantV(audioMaster);}

ToneSlantV::ToneSlantV(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.0;
	B = 0.0;
	for(int count = 0; count < 128; count++) {ts_f[count] = 0.0;}
	for(int count = 0; count < 260; count++) {ts_b[count] = 0.0;}
  gcount = 0;
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

ToneSlantV::~ToneSlantV() {}
VstInt32 ToneSlantV::getVendorVersion () {return 1000;}
void ToneSlantV::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ToneSlantV::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ToneSlantV::getChunk (void** data, bool isPreset)
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

VstInt32 ToneSlantV::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	B = pinParameter(chunkData[1]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void ToneSlantV::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        case kParamB: B = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float ToneSlantV::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        case kParamB: return B; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ToneSlantV::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Voicing", kVstMaxParamStrLen); break;
		case kParamB: vst_strncpy (text, "Highs", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ToneSlantV::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string ((A*99.0)+1.0, text, kVstMaxParamStrLen); break;
        case kParamB: float2string ((B*2.0)-1.0, text, kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ToneSlantV::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "taps", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, " ", kVstMaxParamStrLen); break; //the percent
        default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ToneSlantV::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ToneSlantV::getEffectName(char* name) {
    vst_strncpy(name, "ToneSlantV", kVstMaxProductStrLen); return true;
}

VstPlugCategory ToneSlantV::getPlugCategory() {return kPlugCategEffect;}

bool ToneSlantV::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ToneSlantV", kVstMaxProductStrLen); return true;
}

bool ToneSlantV::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
