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
	A = 0.5;
  // Desk4
	B = 0.0;
	C = 0.0;
	D = 0.0;
	E = 0.0;
	F = 0.0;
	G = 0.0;
  // Buttercomp2
	H = 0.0;
	I = 0.5;
	J = 0.0;
  // VariMu
	K = 0.0;
	L = 0.5;
	M = 0.0;
	N = 0.0;
  // Biquad1
	O = 0.0;
	P = 0.0;
	Q = 0.0;
	R = 0.5;
  // Biquad2
	S = 0.0;
	T = 0.0;
	U = 0.0;
	V = 0.5;
  // BussColors4
	W = 0.0;
	X = 0.5;
	Y = 0.5;
	Z = 0.0;
	// ADClip7
	Z0 = 0.0;
	Z1 = 0.5;
	Z2 = 0.5;
	Z3 = 0.0;
		
	// Monitoring2

	Z4 = 0.0;
	
  // Hypersonic

	for (int x = 0; x < fix_total; x++) {
		fixB[x] = 0.0;
		fixC[x] = 0.0;
		fixF[x] = 0.0;
	}

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

  // Biquad 1

	for (int x = 0; x < 9; x++) {biquad_1_AL[x] = 0.0; biquad_1_AR[x] = 0.0; biquad_1_BL[x] = 0.0; biquad_1_BR[x] = 0.0;}

  // Biquad 2

	for (int x = 0; x < 9; x++) {biquad_2_AL[x] = 0.0; biquad_2_AR[x] = 0.0; biquad_2_BL[x] = 0.0; biquad_2_BR[x] = 0.0;}

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

  // ADClip7

	ad_lastSampleL = 0.0;
	ad_lastSampleR = 0.0;
	for(int count = 0; count < 22199; count++) {ad_bL[count] = 0; ad_bR[count] = 0;}
	gcount = 0;
	lowsL = 0;
	lowsR = 0;
	refclipL = 0.89;
	refclipR = 0.89;
	iirLowsAL = 0.0;
	iirLowsAR = 0.0;
	iirLowsBL = 0.0;
	iirLowsBR = 0.0;

  // Monitoring2

	for(int count = 0; count < 99; count++) {
		darkSampleL[count] = 0;
		darkSampleR[count] = 0;
	}

	for(int count = 0; count < 1502; count++) {
		m_aL[count] = 0.0; m_bL[count] = 0.0; m_cL[count] = 0.0; m_dL[count] = 0.0;
		m_aR[count] = 0.0; m_bR[count] = 0.0; m_cR[count] = 0.0; m_dR[count] = 0.0;
	}
	ax = 1; bx = 1; cx = 1; dx = 1;
	//PeaksOnly
	m_lastSampleL = 0.0; m_lastSampleR = 0.0;
	//SlewOnly
	iirSampleAL = 0.0; iirSampleBL = 0.0; iirSampleCL = 0.0; iirSampleDL = 0.0; iirSampleEL = 0.0; iirSampleFL = 0.0; iirSampleGL = 0.0;
	iirSampleHL = 0.0; iirSampleIL = 0.0; iirSampleJL = 0.0; iirSampleKL = 0.0; iirSampleLL = 0.0; iirSampleML = 0.0; iirSampleNL = 0.0; iirSampleOL = 0.0; iirSamplePL = 0.0;
	iirSampleQL = 0.0; iirSampleRL = 0.0; iirSampleSL = 0.0;
	iirSampleTL = 0.0; iirSampleUL = 0.0; iirSampleVL = 0.0;
	iirSampleWL = 0.0; iirSampleXL = 0.0; iirSampleYL = 0.0; iirSampleZL = 0.0;

	iirSampleAR = 0.0; iirSampleBR = 0.0; iirSampleCR = 0.0; iirSampleDR = 0.0; iirSampleER = 0.0; iirSampleFR = 0.0; iirSampleGR = 0.0;
	iirSampleHR = 0.0; iirSampleIR = 0.0; iirSampleJR = 0.0; iirSampleKR = 0.0; iirSampleLR = 0.0; iirSampleMR = 0.0; iirSampleNR = 0.0; iirSampleOR = 0.0; iirSamplePR = 0.0;
	iirSampleQR = 0.0; iirSampleRR = 0.0; iirSampleSR = 0.0;
	iirSampleTR = 0.0; iirSampleUR = 0.0; iirSampleVR = 0.0;
	iirSampleWR = 0.0; iirSampleXR = 0.0; iirSampleYR = 0.0; iirSampleZR = 0.0; // o/`
	//SubsOnly
	for (int x = 0; x < fix_total; x++) {biquad[x] = 0.0;}
	//Bandpasses

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
	chunkData[0] =  A;
	chunkData[1] =  B;
	chunkData[2] =  C;
	chunkData[3] =  D;
	chunkData[4] =  E;
	chunkData[5] =  F;
	chunkData[6] =  G;
	chunkData[7] =  H;
	chunkData[8] =  I;
	chunkData[9] =  J;
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
	chunkData[26] = Z0;
	chunkData[27] = Z1;
	chunkData[28] = Z2;
	chunkData[29] = Z3;
	chunkData[30] = Z4;
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
	S = pinParameter(chunkData[18]);
	T = pinParameter(chunkData[19]);
	U = pinParameter(chunkData[20]);
	V = pinParameter(chunkData[21]);
	W = pinParameter(chunkData[22]);
	X = pinParameter(chunkData[23]);
	Y = pinParameter(chunkData[24]);
	Z = pinParameter(chunkData[25]);
	Z0 = pinParameter(chunkData[26]);
	Z1 = pinParameter(chunkData[27]);
	Z2 = pinParameter(chunkData[28]);
	Z3 = pinParameter(chunkData[29]);
	Z4 = pinParameter(chunkData[30]);
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
        case kParamS: S = value; break;
        case kParamT: T = value; break;
        case kParamU: U = value; break;
        case kParamV: V = value; break;
        case kParamW: W = value; break;
        case kParamX: X = value; break;
        case kParamY: Y = value; break;
        case kParamZ: Z = value; break;
        case kParamZ0: Z0 = value; break;
        case kParamZ1: Z1 = value; break;
        case kParamZ2: Z2 = value; break;
        case kParamZ3: Z3 = value; break;
        case kParamZ4: Z4 = value; break;
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
        case kParamS: return S; break;
        case kParamT: return T; break;
        case kParamU: return U; break;
        case kParamV: return V; break;
        case kParamW: return W; break;
        case kParamX: return X; break;
        case kParamY: return Y; break;
        case kParamZ: return Z; break;
        case kParamZ0: return Z0; break;
        case kParamZ1: return Z1; break;
        case kParamZ2: return Z2; break;
        case kParamZ3: return Z3; break;
        case kParamZ4: return Z4; break;
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
        case kParamO: vst_strncpy (text, "F1 Type", kVstMaxParamStrLen); break;
				case kParamP: vst_strncpy (text, "F1 Freq", kVstMaxParamStrLen); break;
				case kParamQ: vst_strncpy (text, "F1 Q", kVstMaxParamStrLen); break;
				case kParamR: vst_strncpy (text, "F1 Inv/Wet", kVstMaxParamStrLen); break;
        case kParamS: vst_strncpy (text, "F2 Type", kVstMaxParamStrLen); break;
				case kParamT: vst_strncpy (text, "F2 Freq", kVstMaxParamStrLen); break;
				case kParamU: vst_strncpy (text, "F2 Q", kVstMaxParamStrLen); break;
				case kParamV: vst_strncpy (text, "F2 Inv/Wet", kVstMaxParamStrLen); break;
        case kParamW: vst_strncpy (text, "Color", kVstMaxParamStrLen); break;
        case kParamX: vst_strncpy (text, "Input", kVstMaxParamStrLen); break;
        case kParamY: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
        case kParamZ: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamZ0: vst_strncpy (text, "Boost", kVstMaxParamStrLen); break;
    		case kParamZ1: vst_strncpy (text, "Soften", kVstMaxParamStrLen); break;
    		case kParamZ2: vst_strncpy (text, "Enhance", kVstMaxParamStrLen); break;
    		case kParamZ3: vst_strncpy (text, "Mode", kVstMaxParamStrLen); break;
        case kParamZ4: vst_strncpy (text, "Monitor", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ConsoleZBuss::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: float2string (2.0 * A, text, kVstMaxParamStrLen); break;
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
        case kParamO: float2string ((float)ceil((O*3.999)+0.00001), text, kVstMaxParamStrLen); break;
        case kParamP: float2string ((P*P*P*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamQ: float2string ((Q*Q*Q*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamR: float2string ((R*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamS: float2string ((float)ceil((S*3.999)+0.00001), text, kVstMaxParamStrLen); break;
        case kParamT: float2string ((T*T*T*0.9999)+0.0001, text, kVstMaxParamStrLen); break;
        case kParamU: float2string ((U*U*U*29.99)+0.01, text, kVstMaxParamStrLen); break;
        case kParamV: float2string ((V*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamW: switch((VstInt32)( W * 7.999 )) //0 to almost edge of # of params
				{ case 0: vst_strncpy (text, "Dark", kVstMaxParamStrLen); break;
					case 1: vst_strncpy (text, "Rock", kVstMaxParamStrLen); break;
					case 2: vst_strncpy (text, "Lush", kVstMaxParamStrLen); break;
					case 3: vst_strncpy (text, "Vibe", kVstMaxParamStrLen); break;
					case 4: vst_strncpy (text, "Holo", kVstMaxParamStrLen); break;
					case 5: vst_strncpy (text, "Punch", kVstMaxParamStrLen); break;
					case 6: vst_strncpy (text, "Steel", kVstMaxParamStrLen); break;
					case 7: vst_strncpy (text, "Tube", kVstMaxParamStrLen); break;
					default: break; // unknown parameter, shouldn't happen!
				} break; //completed A 'popup' parameter, exit
        case kParamX: float2string ((X * 36.0)-18.0, text, kVstMaxParamStrLen); break;
        case kParamY: float2string ((Y * 36.0)-18.0, text, kVstMaxParamStrLen); break;
        case kParamZ: float2string (Z, text, kVstMaxParamStrLen); break;
        case kParamZ0: float2string (Z0*18.0, text, kVstMaxParamStrLen); break;
        case kParamZ1: float2string (Z1, text, kVstMaxParamStrLen); break;
        case kParamZ2: float2string (Z2, text, kVstMaxParamStrLen); break;
        case kParamZ3: switch((VstInt32)( Z3 * 2.999 )) //0 to almost edge of # of params
				{ case 0: vst_strncpy (text, "Normal", kVstMaxParamStrLen); break;
				  case 1: vst_strncpy (text, "Atten", kVstMaxParamStrLen); break;
				  case 2: vst_strncpy (text, "Clips", kVstMaxParamStrLen); break;
				  default: break; // unknown parameter, shouldn't happen!
				} break;
        case kParamZ4: switch((VstInt32)( Z4 * 16.999 )) //0 to almost edge of # of params
				{	case kDKAD: vst_strncpy (text, "Out24", kVstMaxParamStrLen); break;
					case kDKCD: vst_strncpy (text, "Out16", kVstMaxParamStrLen); break;
					case kPEAK: vst_strncpy (text, "Peaks", kVstMaxParamStrLen); break;
					case kSLEW: vst_strncpy (text, "Slew", kVstMaxParamStrLen); break;
					case kSUBS: vst_strncpy (text, "Subs", kVstMaxParamStrLen); break;
					case kMONO: vst_strncpy (text, "Mono", kVstMaxParamStrLen); break;
					case kSIDE: vst_strncpy (text, "Side", kVstMaxParamStrLen); break;
					case kVINYL: vst_strncpy (text, "Vinyl", kVstMaxParamStrLen); break;
					case kAURAT: vst_strncpy (text, "Aurat", kVstMaxParamStrLen); break;
					case kMONORAT: vst_strncpy (text, "MonoRat", kVstMaxParamStrLen); break;
					case kMONOLAT: vst_strncpy (text, "MonoLat", kVstMaxParamStrLen); break;
					case kPHONE: vst_strncpy (text, "Phone", kVstMaxParamStrLen); break;
					case kCANSA: vst_strncpy (text, "Cans A", kVstMaxParamStrLen); break;
					case kCANSB: vst_strncpy (text, "Cans B", kVstMaxParamStrLen); break;
					case kCANSC: vst_strncpy (text, "Cans C", kVstMaxParamStrLen); break;
					case kCANSD: vst_strncpy (text, "Cans D", kVstMaxParamStrLen); break;
					case kTRICK: vst_strncpy (text, "V Trick", kVstMaxParamStrLen); break;
					default: break; // unknown parameter, shouldn't happen!
				} break;
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
        case kParamP: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamQ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamR: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamS: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamT: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamU: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamV: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamW: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamX: vst_strncpy (text, "db", kVstMaxParamStrLen); break;
        case kParamY: vst_strncpy (text, "db", kVstMaxParamStrLen); break;
        case kParamZ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamZ0: vst_strncpy (text, "dB", kVstMaxParamStrLen); break;
        case kParamZ1: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamZ2: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamZ3: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamZ4: vst_strncpy (text, "", kVstMaxParamStrLen); break;
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
