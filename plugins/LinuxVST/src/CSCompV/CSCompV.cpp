/* ========================================
 *  CSCompV - CSCompV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __CSCompV_H
#include "CSCompV.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new CSCompV(audioMaster);}

CSCompV::CSCompV(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.0; //Compres 0-1
	B = 0.0; //CompSpd 0-1
	C = 0.33; //Output
		
	//begin ButterComp
	cs_controlpos[0] = 1.0;
	cs_controlpos[1] = 1.0;
	cs_controlpos[2] = 1.0;
	cs_controlpos[3] = 1.0;
	cs_controlneg[0] = 1.0;
	cs_controlneg[1] = 1.0;
	cs_controlneg[2] = 1.0;
	cs_controlneg[3] = 1.0;
	cs_targetpos[0] = 1.0;
	cs_targetpos[1] = 1.0;
	cs_targetpos[2] = 1.0;
	cs_targetpos[3] = 1.0;
	cs_targetneg[0] = 1.0;
	cs_targetneg[1] = 1.0;
	cs_targetneg[2] = 1.0;
	cs_targetneg[3] = 1.0;
	cs_avg[0] = cs_avg[1] = 0.0;
	//end ButterComp	
	
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

CSCompV::~CSCompV() {}
VstInt32 CSCompV::getVendorVersion () {return 1000;}
void CSCompV::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void CSCompV::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 CSCompV::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	chunkData[0] = A;
	chunkData[1] = B;
	chunkData[2] = C;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you 
	 started with. */
	
	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 CSCompV::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	B = pinParameter(chunkData[1]);
	C = pinParameter(chunkData[2]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void CSCompV::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        case kParamB: B = value; break;
        case kParamC: C = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float CSCompV::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        case kParamB: return B; break;
        case kParamC: return C; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void CSCompV::getParameterName(VstInt32 index, char *text) {
    switch (index) {
		case kParamA: vst_strncpy (text, "Compres", kVstMaxParamStrLen); break;
		case kParamB: vst_strncpy (text, "CompSpd", kVstMaxParamStrLen); break;
		case kParamC: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void CSCompV::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (A, text, kVstMaxParamStrLen); break; //Gate 0-1
        case kParamB: float2string (B, text, kVstMaxParamStrLen); break; //Compres 0-1
        case kParamC: float2string (C, text, kVstMaxParamStrLen); break; //CompSpd 0-1
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void CSCompV::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamC: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 CSCompV::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool CSCompV::getEffectName(char* name) {
    vst_strncpy(name, "CSCompV", kVstMaxProductStrLen); return true;
}

VstPlugCategory CSCompV::getPlugCategory() {return kPlugCategEffect;}

bool CSCompV::getProductString(char* text) {
  	vst_strncpy (text, "airwindows CSCompV", kVstMaxProductStrLen); return true;
}

bool CSCompV::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
