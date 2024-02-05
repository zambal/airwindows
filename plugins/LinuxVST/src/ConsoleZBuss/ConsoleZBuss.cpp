/* ========================================
 *  ConsoleZBuss - ConsoleZBuss.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZBuss_H
#include "ConsoleZBuss.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ConsoleZBuss(audioMaster);}

ConsoleZBuss::ConsoleZBuss(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
//LABuss
	A = 1.0;
// Desk4
	B = 0.27;
	C = 0.18;
	D = 0.26;
	E = 0.54;
	F = 0.84;
	G = 0.0;
// Buttercomp2
	H = 0.0;
	I = 0.5;
	J = 0.0;
// VariMu
	K = 0.0;
	L = 0.5;
	M = 1.0;
	N = 0.0;
// BussColors4
	O = 0.0;
	P = 0.5;
	Q = 0.5;
	R = 1.0;


// ConsoleLABuss

	lastSinewL = lastSinewR = 0.0;
	subAL = subAR = subBL = subBR = subCL = subCR = 0.0;
	gainA = gainB = 1.0;

// Desk4

	for(int count = 0; count < 9999; count++) {desk4_dL[count] = 0; desk4_dR[count] = 0;}
	desk4_controlL = 0;
	lastSampleL = 0.0;
	lastOutSampleL = 0.0;
	lastSlewL = 0.0;
	desk4_controlR = 0;
	lastSampleR = 0.0;
	lastOutSampleR = 0.0;
	lastSlewR = 0.0;
	desk4_gcount = 0;

// Buttercomp2

	controlAposL = 1.0;
	controlAnegL = 1.0;
	controlBposL = 1.0;
	controlBnegL = 1.0;
	targetposL = 1.0;
	targetnegL = 1.0;
	lastOutputL = 0.0;

	controlAposR = 1.0;
	controlAnegR = 1.0;
	controlBposR = 1.0;
	controlBnegR = 1.0;
	targetposR = 1.0;
	targetnegR = 1.0;
	lastOutputR = 0.0;

	butter_flip = false;

// VariMu

	muSpeedAL = 10000;
	muSpeedBL = 10000;
	muCoefficientAL = 1;
	muCoefficientBL = 1;
	muVaryL = 1;
	previousL = 0.0;

	muSpeedAR = 10000;
	muSpeedBR = 10000;
	muCoefficientAR = 1;
	muCoefficientBR = 1;
	muVaryR = 1;
	previousR = 0.0;
	flip = false;

// BussColors4

	for (int count = 0; count < 174; count++) {bL[count] = 0; bR[count] = 0;}
	for (int count = 0; count < 99; count++) {dL[count] = 0; dR[count] = 0;}
	for (int count = 0; count < 34; count++) c[count] = count; //initial setup for 44.1K
	g[1] = pow(10.0, -5.2 / 14.0); //dark
	g[2] = pow(10.0, -6.2 / 14.0); //rock
	g[3] = pow(10.0, -2.9 / 14.0); //lush
	g[4] = pow(10.0, -1.1 / 14.0); //vibe
	g[5] = pow(10.0, -5.1 / 14.0); //holo
	g[6] = pow(10.0, -3.6 / 14.0); //punch
	g[7] = pow(10.0, -2.3 / 14.0); //steel
	g[8] = pow(10.0, -2.9 / 14.0); //tube
	//preset gains for models
	outg[1] = pow(10.0, -0.3 / 14.0); //dark
	outg[2] = pow(10.0, 0.5 / 14.0); //rock
	outg[3] = pow(10.0, -0.7 / 14.0); //lush
	outg[4] = pow(10.0, -0.6 / 14.0); //vibe
	outg[5] = pow(10.0, -0.2 / 14.0); //holo
	outg[6] = pow(10.0, 0.3 / 14.0); //punch
	outg[7] = pow(10.0, 0.1 / 14.0); //steel
	outg[8] = pow(10.0, 0.9 / 14.0); //tube
	//preset gains for models
	controlL = 0;
	controlR = 0;
	slowdynL = 0;
	slowdynR = 0;
	gcount = 0;


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

ConsoleZBuss::~ConsoleZBuss() {}
VstInt32 ConsoleZBuss::getVendorVersion () {return 1000;}
void ConsoleZBuss::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ConsoleZBuss::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ConsoleZBuss::getChunk (void** data, bool isPreset)
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
	/* Note: The way this is set up, it will break if you manage to save settings on an Intel
	 machine and load them on a PPC Mac. However, it's fine if you stick to the machine you
	 started with. */

	*data = chunkData;
	return kNumParameters * sizeof(float);
}

VstInt32 ConsoleZBuss::setChunk (void* data, VstInt32 byteSize, bool isPreset)
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
	/* We're ignoring byteSize as we found it to be a filthy liar */

	/* calculate any other fields you need here - you could copy in
	 code from setParameter() here. */
	return 0;
}

void ConsoleZBuss::setParameter(VstInt32 index, float value) {
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
        default: throw; // unknown parameter, shouldn't happen!
    }
}

float ConsoleZBuss::getParameter(VstInt32 index) {
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
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ConsoleZBuss::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Master", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "Overdrive", kVstMaxParamStrLen); break;
        case kParamC: vst_strncpy (text, "Hi Choke", kVstMaxParamStrLen); break;
        case kParamD: vst_strncpy (text, "PowerSag", kVstMaxParamStrLen); break;
        case kParamE: vst_strncpy (text, "Frequency", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "Compress", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        case kParamJ: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamK: vst_strncpy (text, "Intensity", kVstMaxParamStrLen); break;
        case kParamL: vst_strncpy (text, "Speed", kVstMaxParamStrLen); break;
        case kParamM: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        case kParamN: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "Color", kVstMaxParamStrLen); break;
        case kParamP: vst_strncpy (text, "Input", kVstMaxParamStrLen); break;
        case kParamQ: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        case kParamR: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ConsoleZBuss::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (A, text, kVstMaxParamStrLen); break;
        case kParamB: float2string (B, text, kVstMaxParamStrLen); break;
        case kParamC: float2string (C, text, kVstMaxParamStrLen); break;
        case kParamD: float2string (D, text, kVstMaxParamStrLen); break;
        case kParamE: float2string (E, text, kVstMaxParamStrLen); break;
        case kParamF: float2string (F, text, kVstMaxParamStrLen); break;
        case kParamG: float2string (G, text, kVstMaxParamStrLen); break;
        case kParamH: float2string (H, text, kVstMaxParamStrLen); break;
        case kParamI: float2string (2.0 * I, text, kVstMaxParamStrLen); break;
        case kParamJ: float2string (J, text, kVstMaxParamStrLen); break;
        case kParamK: float2string (K, text, kVstMaxParamStrLen); break;
        case kParamL: float2string (L, text, kVstMaxParamStrLen); break;
        case kParamM: float2string (M, text, kVstMaxParamStrLen); break;
        case kParamN: float2string (N, text, kVstMaxParamStrLen); break;
        case kParamO: switch((VstInt32)( O * 7.999 )) //0 to almost edge of # of params
		{case 0: vst_strncpy (text, "Dark", kVstMaxParamStrLen); break;
			case 1: vst_strncpy (text, "Rock", kVstMaxParamStrLen); break;
			case 2: vst_strncpy (text, "Lush", kVstMaxParamStrLen); break;
			case 3: vst_strncpy (text, "Vibe", kVstMaxParamStrLen); break;
			case 4: vst_strncpy (text, "Holo", kVstMaxParamStrLen); break;
			case 5: vst_strncpy (text, "Punch", kVstMaxParamStrLen); break;
			case 6: vst_strncpy (text, "Steel", kVstMaxParamStrLen); break;
			case 7: vst_strncpy (text, "Tube", kVstMaxParamStrLen); break;
			default: break; // unknown parameter, shouldn't happen!
		} break; //completed A 'popup' parameter, exit
        case kParamP: float2string ((P * 36.0)-18.0, text, kVstMaxParamStrLen); break;
        case kParamQ: float2string ((Q * 36.0)-18.0, text, kVstMaxParamStrLen); break;
        case kParamR: float2string (R, text, kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ConsoleZBuss::getParameterLabel(VstInt32 index, char *text) {
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
        case kParamP: vst_strncpy (text, "db", kVstMaxParamStrLen); break;
        case kParamQ: vst_strncpy (text, "db", kVstMaxParamStrLen); break;
        case kParamR: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ConsoleZBuss::canDo(char *text)
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ConsoleZBuss::getEffectName(char* name) {
    vst_strncpy(name, "ConsoleZBuss", kVstMaxProductStrLen); return true;
}

VstPlugCategory ConsoleZBuss::getPlugCategory() {return kPlugCategEffect;}

bool ConsoleZBuss::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ConsoleZBuss", kVstMaxProductStrLen); return true;
}

bool ConsoleZBuss::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
