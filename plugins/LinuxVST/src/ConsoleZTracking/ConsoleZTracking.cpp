/* ========================================
 *  ConsoleZTracking - ConsoleZTracking.cpp
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTracking_H
#include "ConsoleZTracking.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ConsoleZTracking(audioMaster);}

ConsoleZTracking::ConsoleZTracking(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.5;
	B = 0.0;
	C = 0.0;
	D = 0.0;
	E = 0.0;
	F = 0.5;
	G = 0.5;
	H = 0.5;
	I = 0.0;
	J = 0.0;
	K = 0.5;
	L = 0.0;
	M = 0.5;
	N = 0.0;
	O = 0.0;
	P = 0.0;

  // Hypersonic

	for (int x = 0; x < fix_total; x++) {
		fixA[x] = 0.0;
		fixB[x] = 0.0;
		fixC[x] = 0.0;
		fixD[x] = 0.0;
		fixE[x] = 0.0;
	}

  // Interstage

	iirSampleAL = 0.0;
	iirSampleBL = 0.0;
	iirSampleCL = 0.0;
	iirSampleDL = 0.0;
	iirSampleEL = 0.0;
	iirSampleFL = 0.0;
	iirSampleAR = 0.0;
	iirSampleBR = 0.0;
	iirSampleCR = 0.0;
	iirSampleDR = 0.0;
	iirSampleER = 0.0;
	iirSampleFR = 0.0;

  // BiquadOneHalf HPF

	for (int x = 0; x < 9; x++) {biquad_hpf_AL[x] = 0.0; biquad_hpf_AR[x] = 0.0; biquad_hpf_BL[x] = 0.0; biquad_hpf_BR[x] = 0.0;}
  
  // BiquadOneHalf LPF

	for (int x = 0; x < 9; x++) {biquad_lpf_AL[x] = 0.0; biquad_lpf_AR[x] = 0.0; biquad_lpf_BL[x] = 0.0; biquad_lpf_BR[x] = 0.0;}
  
  // Tape

	iirMidRollerAL = 0.0;
	iirMidRollerBL = 0.0;
	iirHeadBumpAL = 0.0;
	iirHeadBumpBL = 0.0;
	iirMidRollerAR = 0.0;
	iirMidRollerBR = 0.0;
	iirHeadBumpAR = 0.0;
	iirHeadBumpBR = 0.0;
	for (int x = 0; x < 9; x++) {
		biquadAL[x] = 0.0;biquadBL[x] = 0.0;biquadCL[x] = 0.0;biquadDL[x] = 0.0;
		biquadAR[x] = 0.0;biquadBR[x] = 0.0;biquadCR[x] = 0.0;biquadDR[x] = 0.0;
	}
	flip = false;
	lastSampleL = 0.0;
	lastSampleR = 0.0;

  // Creature

	for (int x = 0; x < 101; x++) {
		slewL[x] = 0.0;
		slewR[x] = 0.0;
	}

  // ToneSlant

	for(int count = 0; count < 102; count++) {ts_bL[count] = 0.0; ts_bR[count] = 0.0; ts_f[count] = 0.0;}


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

ConsoleZTracking::~ConsoleZTracking() {}
VstInt32 ConsoleZTracking::getVendorVersion () {return 1000;}
void ConsoleZTracking::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ConsoleZTracking::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ConsoleZTracking::getChunk (void** data, bool isPreset)
{
	float *chunkData = (float *)calloc(kNumParameters, sizeof(float));
	chunkData[0] = A;
	chunkData[1] = B;
	chunkData[2] = C;
	chunkData[3] = D;
	chunkData[4] = E;
	chunkData[5] = F;
	chunkData[6] = G;
	chunkData[7] = H;
	chunkData[8] = I;
	chunkData[9] = J;
  chunkData[10] = K;
	chunkData[11] = L;
	chunkData[12] = M;
	chunkData[13] = N;
  chunkData[14] = O;
  chunkData[15] = P;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you
	 started with. */

	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 ConsoleZTracking::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	B = pinParameter(chunkData[1]);
	C = pinParameter(chunkData[2]);
	D = pinParameter(chunkData[3]);
	E = pinParameter(chunkData[4]);
	F = pinParameter(chunkData[5]);
	G = pinParameter(chunkData[6]);
	H = pinParameter(chunkData[7]);
	I = pinParameter(chunkData[8]);
	J = pinParameter(chunkData[9]);
	K = pinParameter(chunkData[10]);
	L = pinParameter(chunkData[11]);
	M = pinParameter(chunkData[12]);
	N = pinParameter(chunkData[13]);
	O = pinParameter(chunkData[14]);
	P = pinParameter(chunkData[15]);
	/* We're ignoring byteSize as we found it to be a filthy liar */

	/* calculate any other fields you need here - you could copy in
	 code from setParameter() here. */
	return 0;
}

void ConsoleZTracking::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        case kParamB: B = value; break;
        case kParamC: C = value; break;
        case kParamD: D = value; break;
        case kParamE: E = value; break;
        case kParamF: F = value; break;
        case kParamG: G = value; break;
        case kParamH: H = value; break;
        case kParamI: I = value; break;
        case kParamJ: J = value; break;
        case kParamK: K = value; break;
        case kParamL: L = value; break;
        case kParamM: M = value; break;
        case kParamN: N = value; break;
        case kParamO: O = value; break;
        case kParamP: P = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float ConsoleZTracking::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        case kParamB: return B; break;
        case kParamC: return C; break;
        case kParamD: return D; break;
        case kParamE: return E; break;
        case kParamF: return F; break;
        case kParamG: return G; break;
        case kParamH: return H; break;
        case kParamI: return I; break;
        case kParamJ: return J; break;
        case kParamK: return K; break;
        case kParamL: return L; break;
        case kParamM: return M; break;
        case kParamN: return N; break;
        case kParamO: return O; break;
        case kParamP: return P; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ConsoleZTracking::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Input", kVstMaxParamStrLen); break;
    		case kParamB: vst_strncpy (text, "HPF Freq", kVstMaxParamStrLen); break;
    		case kParamC: vst_strncpy (text, "HPF Q", kVstMaxParamStrLen); break;
    		case kParamD: vst_strncpy (text, "LPF Freq", kVstMaxParamStrLen); break;
    		case kParamE: vst_strncpy (text, "LPF Q", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "Slam", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "Bump", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "Trim", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "Drive", kVstMaxParamStrLen); break;
    		case kParamJ: vst_strncpy (text, "Depth", kVstMaxParamStrLen); break;
    		case kParamK: vst_strncpy (text, "Inv/Wet", kVstMaxParamStrLen); break;
        case kParamL: vst_strncpy (text, "Voicing", kVstMaxParamStrLen); break;
    		case kParamM: vst_strncpy (text, "Highs", kVstMaxParamStrLen); break;
        case kParamN: vst_strncpy (text, "HPF", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "LPF", kVstMaxParamStrLen); break;
        case kParamP: vst_strncpy (text, "uLaw", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ConsoleZTracking::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: int2string ((VstInt32)((A * 32)-16), text, kVstMaxParamStrLen); break;
        case kParamB: float2string ((B*B*B*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamC: float2string ((C*C*C*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamD: float2string ((D*D*D*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamE: float2string ((E*E*E*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamF: float2string ((F-0.5)*24.0, text, kVstMaxParamStrLen); break;
        case kParamG: float2string (G, text, kVstMaxParamStrLen); break;
        case kParamH: int2string ((VstInt32)((H * 32)-16), text, kVstMaxParamStrLen); break;
        case kParamI: float2string (I, text, kVstMaxParamStrLen); break;
        case kParamJ: float2string (J, text, kVstMaxParamStrLen); break;
        case kParamK: float2string (K, text, kVstMaxParamStrLen); break;
        case kParamL: float2string ((L*99.0)+1.0, text, kVstMaxParamStrLen); break;
        case kParamM: float2string ((M*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamN: if(N == 0.0) {vst_strncpy (text, "Off", kVstMaxParamStrLen);} else {vst_strncpy (text, "On", kVstMaxParamStrLen);} break;
        case kParamO: if(O == 0.0) {vst_strncpy (text, "Off", kVstMaxParamStrLen);} else {vst_strncpy (text, "On", kVstMaxParamStrLen);} break;
        case kParamP: if(P == 0.0) {vst_strncpy (text, "Off", kVstMaxParamStrLen);} else {vst_strncpy (text, "On", kVstMaxParamStrLen);} break;
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ConsoleZTracking::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "bits", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamC: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamD: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamE: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "dB", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "bits", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamJ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamK: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamL: vst_strncpy (text, "taps", kVstMaxParamStrLen); break;
        case kParamM: vst_strncpy (text, " ", kVstMaxParamStrLen); break; //the percent
        case kParamN: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamP: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ConsoleZTracking::canDo(char *text)
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ConsoleZTracking::getEffectName(char* name) {
    vst_strncpy(name, "ConsoleZTracking", kVstMaxProductStrLen); return true;
}

VstPlugCategory ConsoleZTracking::getPlugCategory() {return kPlugCategEffect;}

bool ConsoleZTracking::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ConsoleZTracking", kVstMaxProductStrLen); return true;
}

bool ConsoleZTracking::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
