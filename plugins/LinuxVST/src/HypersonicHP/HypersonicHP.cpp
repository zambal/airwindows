/* ========================================
 *  HypersonicHP - HypersonicHP.h
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __HypersonicHP_H
#include "HypersonicHP.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new HypersonicHP(audioMaster);}

HypersonicHP::HypersonicHP(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	for (int x = 0; x < fix_total; x++) {
		fixA[x] = 0.0;
		fixB[x] = 0.0;
		fixC[x] = 0.0;
		fixD[x] = 0.0;
		fixE[x] = 0.0;
		fixF[x] = 0.0;
		fixG[x] = 0.0;
	}

  fixA[fix_reso] = 4.46570214;
	fixB[fix_reso] = 1.51387132;
	fixC[fix_reso] = 0.93979296;
	fixD[fix_reso] = 0.70710678;
	fixE[fix_reso] = 0.59051105;
	fixF[fix_reso] = 0.52972649;
	fixG[fix_reso] = 0.50316379;

	fpd_b[0] = 1.0; while (fpd_b[0] < 16386) fpd_b[0] = rand()*UINT32_MAX;
	fpd_b[1] = 1.0; while (fpd_b[1] < 16386) fpd_b[1]= rand()*UINT32_MAX;
  fpd_b[2] = fpd_b[0]; fpd_b[2] ^= fpd_b[2] << 13; fpd_b[2] ^= fpd_b[2] >> 17; fpd_b[2] ^= fpd_b[2] << 5; 
  fpd_b[3] = fpd_b[1]; fpd_b[3] ^= fpd_b[3] << 13; fpd_b[3] ^= fpd_b[3] >> 17; fpd_b[3] ^= fpd_b[3] << 5;
	//this is reset: values being initialized only once. Startup values, whatever they are.
	
    setSampleRate(getSampleRate());
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

HypersonicHP::~HypersonicHP() {}
VstInt32 HypersonicHP::getVendorVersion () {return 1000;}
void HypersonicHP::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void HypersonicHP::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 HypersonicHP::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you 
	 started with. */
	
	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 HypersonicHP::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void HypersonicHP::setParameter(VstInt32 index, float value) {
    switch (index) {
		default: throw; // unknown parameter, shouldn't happen!
    }
}

float HypersonicHP::getParameter(VstInt32 index) {
    switch (index) {
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void HypersonicHP::setSampleRate(float sampleRate) {
  this->sampleRate = sampleRate;

  double cutoff = 25000.0 / sampleRate;
  if (cutoff > 0.49) cutoff = 0.49; //don't crash if run at 44.1k

  fixG[fix_freq] = fixF[fix_freq] = fixE[fix_freq] = fixD[fix_freq] = fixC[fix_freq] = fixB[fix_freq] = fixA[fix_freq] = cutoff;


  double K = tan(M_PI * fixA[fix_freq]); //lowpass
  double norm = 1.0 / (1.0 + K / fixA[fix_reso] + K * K);
  fixA[fix_a0] = K * K * norm;
  fixA[fix_a1] = 2.0 * fixA[fix_a0];
  fixA[fix_a2] = fixA[fix_a0];
  fixA[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixA[fix_b2] = (1.0 - K / fixA[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixB[fix_freq]);
  norm = 1.0 / (1.0 + K / fixB[fix_reso] + K * K);
  fixB[fix_a0] = K * K * norm;
  fixB[fix_a1] = 2.0 * fixB[fix_a0];
  fixB[fix_a2] = fixB[fix_a0];
  fixB[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixB[fix_b2] = (1.0 - K / fixB[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixC[fix_freq]);
  norm = 1.0 / (1.0 + K / fixC[fix_reso] + K * K);
  fixC[fix_a0] = K * K * norm;
  fixC[fix_a1] = 2.0 * fixC[fix_a0];
  fixC[fix_a2] = fixC[fix_a0];
  fixC[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixC[fix_b2] = (1.0 - K / fixC[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixD[fix_freq]);
  norm = 1.0 / (1.0 + K / fixD[fix_reso] + K * K);
  fixD[fix_a0] = K * K * norm;
  fixD[fix_a1] = 2.0 * fixD[fix_a0];
  fixD[fix_a2] = fixD[fix_a0];
  fixD[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixD[fix_b2] = (1.0 - K / fixD[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixE[fix_freq]);
  norm = 1.0 / (1.0 + K / fixE[fix_reso] + K * K);
  fixE[fix_a0] = K * K * norm;
  fixE[fix_a1] = 2.0 * fixE[fix_a0];
  fixE[fix_a2] = fixE[fix_a0];
  fixE[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixE[fix_b2] = (1.0 - K / fixE[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixF[fix_freq]);
  norm = 1.0 / (1.0 + K / fixF[fix_reso] + K * K);
  fixF[fix_a0] = K * K * norm;
  fixF[fix_a1] = 2.0 * fixF[fix_a0];
  fixF[fix_a2] = fixF[fix_a0];
  fixF[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixF[fix_b2] = (1.0 - K / fixF[fix_reso] + K * K) * norm;

  K = tan(M_PI * fixG[fix_freq]);
  norm = 1.0 / (1.0 + K / fixG[fix_reso] + K * K);
  fixG[fix_a0] = K * K * norm;
  fixG[fix_a1] = 2.0 * fixG[fix_a0];
  fixG[fix_a2] = fixG[fix_a0];
  fixG[fix_b1] = 2.0 * (K * K - 1.0) * norm;
  fixG[fix_b2] = (1.0 - K / fixG[fix_reso] + K * K) * norm;
}

void HypersonicHP::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void HypersonicHP::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void HypersonicHP::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 HypersonicHP::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool HypersonicHP::getEffectName(char* name) {
    vst_strncpy(name, "HypersonicHP", kVstMaxProductStrLen); return true;
}

VstPlugCategory HypersonicHP::getPlugCategory() {return kPlugCategEffect;}

bool HypersonicHP::getProductString(char* text) {
  	vst_strncpy (text, "airwindows HypersonicHP", kVstMaxProductStrLen); return true;
}

bool HypersonicHP::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
