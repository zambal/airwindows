/* ========================================
 *  ConsoleZChannel - ConsoleZChannel.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZChannel_H
#include "ConsoleZChannel.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ConsoleZChannel(audioMaster);}

ConsoleZChannel::ConsoleZChannel(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
  // Hypersonic

	for (int x = 0; x < fix_total; x++) {
		fixA[x] = 0.0;
		fixD[x] = 0.0;
		fixE[x] = 0.0;
		fixG[x] = 0.0;
	}

  // Biquad 1

	A = 0.0;
	B = 0.0;
  C = 0.0;
  D = 0.5;

	for (int x = 0; x < 9; x++) {biquad_1_AL[x] = 0.0; biquad_1_AR[x] = 0.0; biquad_1_BL[x] = 0.0; biquad_1_BR[x] = 0.0;}

  // Biquad 2
	E = 0.0;
	F = 0.0;
	G = 0.0;
	H = 0.5;

	for (int x = 0; x < 9; x++) {biquad_2_AL[x] = 0.0; biquad_2_AR[x] = 0.0; biquad_2_BL[x] = 0.0; biquad_2_BR[x] = 0.0;}

  // Desk 4

	I = 0.0;
	J = 0.0;
	K = 0.0;
	L = 0.0;
	M = 0.0;
	N = 0.0;

	for(int count = 0; count < 9999; count++) {dL[count] = 0; dR[count] = 0;}
	controlL = 0;
	lastSampleL = 0.0;
	lastOutSampleL = 0.0;
	lastSlewL = 0.0;
	controlR = 0;
	lastSampleR = 0.0;
	lastOutSampleR = 0.0;
	lastSlewR = 0.0;
	gcount = 0;

  // PurestDrive

  O = 0.0;

	pd_previousSampleL = 0.0;
	pd_previousSampleR = 0.0;

  // Res2

  P = 0.5;
	Q = 0.0;

	for(int count = 0; count < 2004; count++) {mpkL[count] = 0.0; mpkR[count] = 0.0;}
	for(int count = 0; count < 65; count++) {f[count] = 0.0;}
	prevfreqMPeak = -1;
	prevamountMPeak = -1;
	mpc = 1;

  // ConsoleLA

  R = 0.5;
	S = 0.5;
	T = 0.5;
	U = 0.5;
	V = 0.5;

	for(int count = 0; count < 222; count++) {hullL[count] = 0.0; hullR[count] = 0.0;}
	hullp = 1;
	for (int x = 0; x < 21; x++) pearB[x] = 0.0;
	subAL = subAR = subBL = subBR = subCL = subCR = 0.0;
	midA = midB = 0.0;
	bassA = bassB = 0.0;
	gainA = gainB = 1.0;

	// CStrip2 comp

	W = 0.0; //Compres 0-1
	X = 0.0; //CompSpd 0-1

	//begin ButterComp
	controlAposL = 1.0;
	controlAnegL = 1.0;
	controlBposL = 1.0;
	controlBnegL = 1.0;
	targetposL = 1.0;
	targetnegL = 1.0;
	avgLA = avgLB = 0.0;
	nvgLA = nvgLB = 0.0;

	controlAposR = 1.0;
	controlAnegR = 1.0;
	controlBposR = 1.0;
	controlBnegR = 1.0;
	targetposR = 1.0;
	targetnegR = 1.0;
	avgRA = avgRB = 0.0;
	nvgRA = nvgRB = 0.0;
  flip = false;
	//end ButterComp

  // Wider

	Y = 0.5;
	Z = 0.5;

	for(int fcount = 0; fcount < 4098; fcount++) {p[fcount] = 0.0;}
	count = 0;

  // FP-Dither

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

ConsoleZChannel::~ConsoleZChannel() {}
VstInt32 ConsoleZChannel::getVendorVersion () {return 1000;}
void ConsoleZChannel::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ConsoleZChannel::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ConsoleZChannel::getChunk (void** data, bool isPreset)
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
	chunkData[16] = Q;
	chunkData[17] = R;
	chunkData[18] = S;
	chunkData[19] = T;
	chunkData[20] = U;
	chunkData[21] = V;
	chunkData[22] = W;
	chunkData[23] = X;
	chunkData[24] = Y;
	chunkData[25] = Z;
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you
	 started with. */

	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 ConsoleZChannel::setChunk (void* data, VstInt32 byteSize, bool isPreset)
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
	Q = pinParameter(chunkData[16]);
	R = pinParameter(chunkData[17]);
	S = pinParameter(chunkData[18]);
	T = pinParameter(chunkData[19]);
	U = pinParameter(chunkData[20]);
	V = pinParameter(chunkData[21]);
	W = pinParameter(chunkData[22]);
	X = pinParameter(chunkData[23]);
	Y = pinParameter(chunkData[24]);
	Z = pinParameter(chunkData[25]);
	/* We're ignoring byteSize as we found it to be a filthy liar */

	/* calculate any other fields you need here - you could copy in
	 code from setParameter() here. */
	return 0;
}

void ConsoleZChannel::setParameter(VstInt32 index, float value) {
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
        case kParamQ: Q = value; break;
        case kParamR: R = value; break;
        case kParamS: S = value; break;
        case kParamT: T = value; break;
        case kParamU: U = value; break;
        case kParamV: V = value; break;
        case kParamW: W = value; break;
        case kParamX: X = value; break;
        case kParamY: Y = value; break;
        case kParamZ: Z = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float ConsoleZChannel::getParameter(VstInt32 index) {
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
        case kParamQ: return Q; break;
        case kParamR: return R; break;
        case kParamS: return S; break;
        case kParamT: return T; break;
        case kParamU: return U; break;
        case kParamV: return V; break;
        case kParamW: return W; break;
        case kParamX: return X; break;
        case kParamY: return Y; break;
        case kParamZ: return Z; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ConsoleZChannel::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "F1 Type", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "F1 Freq", kVstMaxParamStrLen); break;
        case kParamC: vst_strncpy (text, "F1 Res", kVstMaxParamStrLen); break;
        case kParamD: vst_strncpy (text, "F1 D/W", kVstMaxParamStrLen); break;
        case kParamE: vst_strncpy (text, "F2 Type", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "F2 Freq", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "F2 Res", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "F2 D/W", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "Overdrive", kVstMaxParamStrLen); break;
    		case kParamJ: vst_strncpy (text, "Hi Choke", kVstMaxParamStrLen); break;
    		case kParamK: vst_strncpy (text, "PowerSag", kVstMaxParamStrLen); break;
    		case kParamL: vst_strncpy (text, "Frequency", kVstMaxParamStrLen); break;
    		case kParamM: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
    		case kParamN: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "PDrive", kVstMaxParamStrLen); break;
        case kParamP: vst_strncpy (text, "MSweep", kVstMaxParamStrLen); break;
    		case kParamQ: vst_strncpy (text, "MBoost", kVstMaxParamStrLen); break;
        case kParamR: vst_strncpy (text, "Treble", kVstMaxParamStrLen); break;
    		case kParamS: vst_strncpy (text, "Mid", kVstMaxParamStrLen); break;
    		case kParamT: vst_strncpy (text, "Bass", kVstMaxParamStrLen); break;
    		case kParamU: vst_strncpy (text, "Pan", kVstMaxParamStrLen); break;
    		case kParamV: vst_strncpy (text, "Fader", kVstMaxParamStrLen); break;
    		case kParamW: vst_strncpy (text, "Compres", kVstMaxParamStrLen); break;
    		case kParamX: vst_strncpy (text, "CompSpd", kVstMaxParamStrLen); break;
        case kParamY: vst_strncpy (text, "Width", kVstMaxParamStrLen); break;
    		case kParamZ: vst_strncpy (text, "Center", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ConsoleZChannel::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string ((float)ceil((A*3.999)+0.00001), text, kVstMaxParamStrLen); break;
        case kParamB: float2string ((B*B*B*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamC: float2string ((C*C*C*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamD: float2string ((D*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamE: float2string ((float)ceil((E*3.999)+0.00001), text, kVstMaxParamStrLen); break;
        case kParamF: float2string ((F*F*F*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamG: float2string ((G*G*G*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamH: float2string ((H*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamI: float2string (I, text, kVstMaxParamStrLen); break;
        case kParamJ: float2string (J, text, kVstMaxParamStrLen); break;
        case kParamK: float2string (K, text, kVstMaxParamStrLen); break;
        case kParamL: float2string (L, text, kVstMaxParamStrLen); break;
        case kParamM: float2string (M, text, kVstMaxParamStrLen); break;
        case kParamN: float2string (N, text, kVstMaxParamStrLen); break;
        case kParamO: float2string (O, text, kVstMaxParamStrLen); break;
        case kParamP: float2string (P, text, kVstMaxParamStrLen); break;
        case kParamQ: float2string (Q, text, kVstMaxParamStrLen); break;
        case kParamR: float2string (R, text, kVstMaxParamStrLen); break;
        case kParamS: float2string (S, text, kVstMaxParamStrLen); break;
        case kParamT: float2string (T, text, kVstMaxParamStrLen); break;
        case kParamU: float2string (U, text, kVstMaxParamStrLen); break;
        case kParamV: float2string (V, text, kVstMaxParamStrLen); break;
        case kParamW: float2string (W, text, kVstMaxParamStrLen); break;
        case kParamX: float2string (X, text, kVstMaxParamStrLen); break;
        case kParamY: float2string ((Y*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamZ: float2string ((Z*2.0)-1.0, text, kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ConsoleZChannel::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamC: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamD: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamE: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamJ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamK: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamL: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamM: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamN: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamP: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamQ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamR: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamS: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamT: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamU: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamV: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamW: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamX: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamY: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamZ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ConsoleZChannel::canDo(char *text)
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ConsoleZChannel::getEffectName(char* name) {
    vst_strncpy(name, "ConsoleZChannel", kVstMaxProductStrLen); return true;
}

VstPlugCategory ConsoleZChannel::getPlugCategory() {return kPlugCategEffect;}

bool ConsoleZChannel::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ConsoleZChannel", kVstMaxProductStrLen); return true;
}

bool ConsoleZChannel::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
