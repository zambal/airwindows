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
	D = 1.0;
	E = 0.0;
	F = 0.0;
	G = 0.5;
	H = 1.0;
	I = 0.0;
	J = 0.0;
	K = 0.5;
	L = 0.0;
	M = 0.0;
	N = 0.5;
	O = 0.5;
	P = 0.5;
	Q = 0.0;
	R = 1.0;
	S = 0.0;
	T = 0.26;
	U = 0.26;
	V = 0.5;
	W = 0.0;
	X = 0.5;
	Y = 0.5;
	Z = 0.0;


	// BitshiftGain

	input_gain = 1.0;

  // BiquadOneHalf HPF

	for (int x = 0; x < 9; x++) {biquad_hpf_AL[x] = 0.0; biquad_hpf_AR[x] = 0.0; biquad_hpf_BL[x] = 0.0; biquad_hpf_BR[x] = 0.0;}

  // BiquadOneHalf LPF

	for (int x = 0; x < 9; x++) {biquad_lpf_AL[x] = 0.0; biquad_lpf_AR[x] = 0.0; biquad_lpf_BL[x] = 0.0; biquad_lpf_BR[x] = 0.0;}

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

  // StereoChorus

	for(int count = 0; count < 65535; count++) {pL[count] = 0;pR[count] = 0;}
	sweepL = 3.141592653589793238 / 2.7;
	sweepR = 3.141592653589793238;
	gcount = 0;
	airPrevL = 0.0;
	airEvenL = 0.0;
	airOddL = 0.0;
	airFactorL = 0.0;
	airPrevR = 0.0;
	airEvenR = 0.0;
	airOddR = 0.0;
	airFactorR = 0.0;
	
	for(int count = 0; count < 6; count++) {lastRefL[count] = 0.0;lastRefR[count] = 0.0;}
	chorus_cycle = 0;
  chorus_flip = false;

  // ClearCoat

	for(int count = 0; count < kshortA+2; count++) {aAL[count] = 0.0; aAR[count] = 0.0;}
	for(int count = 0; count < kshortB+2; count++) {aBL[count] = 0.0; aBR[count] = 0.0;}
	for(int count = 0; count < kshortC+2; count++) {aCL[count] = 0.0; aCR[count] = 0.0;}
	for(int count = 0; count < kshortD+2; count++) {aDL[count] = 0.0; aDR[count] = 0.0;}
	for(int count = 0; count < kshortE+2; count++) {aEL[count] = 0.0; aER[count] = 0.0;}
	for(int count = 0; count < kshortF+2; count++) {aFL[count] = 0.0; aFR[count] = 0.0;}
	for(int count = 0; count < kshortG+2; count++) {aGL[count] = 0.0; aGR[count] = 0.0;}
	for(int count = 0; count < kshortH+2; count++) {aHL[count] = 0.0; aHR[count] = 0.0;}
	for(int count = 0; count < kshortI+2; count++) {aIL[count] = 0.0; aIR[count] = 0.0;}
	for(int count = 0; count < kshortJ+2; count++) {aJL[count] = 0.0; aJR[count] = 0.0;}
	for(int count = 0; count < kshortK+2; count++) {aKL[count] = 0.0; aKR[count] = 0.0;}
	for(int count = 0; count < kshortL+2; count++) {aLL[count] = 0.0; aLR[count] = 0.0;}
	for(int count = 0; count < kshortM+2; count++) {aML[count] = 0.0; aMR[count] = 0.0;}
	for(int count = 0; count < kshortN+2; count++) {aNL[count] = 0.0; aNR[count] = 0.0;}
	for(int count = 0; count < kshortO+2; count++) {aOL[count] = 0.0; aOR[count] = 0.0;}
	for(int count = 0; count < kshortP+2; count++) {aPL[count] = 0.0; aPR[count] = 0.0;}
	
	feedbackAL = 0.0;
	feedbackBL = 0.0;
	feedbackCL = 0.0;
	feedbackDL = 0.0;
	
	previousAL = 0.0;
	previousBL = 0.0;
	previousCL = 0.0;
	previousDL = 0.0;
	previousEL = 0.0;
	
	feedbackDR = 0.0;
	feedbackHR = 0.0;
	feedbackLR = 0.0;
	feedbackPR = 0.0;
	
	previousAR = 0.0;
	previousBR = 0.0;
	previousCR = 0.0;
	previousDR = 0.0;
	previousER = 0.0;
	
	prevMulchAL = 0.0;
	prevMulchAR = 0.0;
	
	tailL = 0.0;
	tailR = 0.0;
	
	for(int count = 0; count < 6; count++) {cc_LastRefL[count] = 0.0; cc_LastRefR[count] = 0.0;}
	
	countAL = 1;
	countBL = 1;
	countCL = 1;
	countDL = 1;	
	countEL = 1;
	countFL = 1;
	countGL = 1;
	countHL = 1;
	countIL = 1;
	countJL = 1;
	countKL = 1;
	countLL = 1;
	countML = 1;
	countNL = 1;
	countOL = 1;
	countPL = 1;
	
	countAR = 1;
	countBR = 1;
	countCR = 1;
	countDR = 1;	
	countER = 1;
	countFR = 1;
	countGR = 1;
	countHR = 1;
	countIR = 1;
	countJR = 1;
	countKR = 1;
	countLR = 1;
	countMR = 1;
	countNR = 1;
	countOR = 1;
	countPR = 1;
	
	cc_cycle = 0;
	
	shortA = 336;
	shortB = 1660;
	shortC = 386;
	shortD = 623;
	shortE = 693;
	shortF = 1079;
	shortG = 891;
	shortH = 1574;
	shortI = 24;
	shortJ = 2641;
	shortK = 1239;
	shortL = 775;
	shortM = 11;
	shortN = 3104;
	shortO = 55;
	shortP = 2366;
	prevclearcoat = -1;
	
	subAL = subAR = subBL = subBR = subCL = subCR = subDL = subDR = 0.0;

  // Channel9

	for (int x = 0; x < 15; x++) {ch9_biquadA[x] = 0.0; ch9_biquadB[x] = 0.0;}
	ch9_iirSampleLA = 0.0;
	ch9_iirSampleRA = 0.0;
	ch9_iirSampleLB = 0.0;
	ch9_iirSampleRB = 0.0;
	ch9_lastSampleAL = ch9_lastSampleBL = ch9_lastSampleCL = 0.0;
	ch9_lastSampleAR = ch9_lastSampleBR = ch9_lastSampleCR = 0.0;
	ch9_iirAmount = 0.005832;
	ch9_threshold = 0.33362176; //instantiating with Neve values
	ch9_cutoff = 28811.0;
  
  // Focus

	for (int x = 0; x < 9; x++) {figureL[x] = 0.0;figureR[x] = 0.0;}
  freqX = 2.0 * log((16.0 - 4.0) / (4.0 - 1.0));

  // Creature

	for (int x = 0; x < 101; x++) {
		slewL[x] = 0.0;
		slewR[x] = 0.0;
	}

  // ToneSlant

	for(int count = 0; count < 102; count++) {ts_bL[count] = 0.0; ts_bR[count] = 0.0; ts_f[count] = 0.0;}

  adcLastSampleL = 0.0;
  adcLastSampleR = 0.0;

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

  	// BitshiftGain

  	switch ((int)(A * 32)-16)
  	{
  		case -16: input_gain = 0.0000152587890625; break;
  		case -15: input_gain = 0.000030517578125; break;
  		case -14: input_gain = 0.00006103515625; break;
  		case -13: input_gain = 0.0001220703125; break;
  		case -12: input_gain = 0.000244140625; break;
  		case -11: input_gain = 0.00048828125; break;
  		case -10: input_gain = 0.0009765625; break;
  		case -9: input_gain = 0.001953125; break;
  		case -8: input_gain = 0.00390625; break;
  		case -7: input_gain = 0.0078125; break;
  		case -6: input_gain = 0.015625; break;
  		case -5: input_gain = 0.03125; break;
  		case -4: input_gain = 0.0625; break;
  		case -3: input_gain = 0.125; break;
  		case -2: input_gain = 0.25; break;
  		case -1: input_gain = 0.5; break;
  		case 0: input_gain = 1.0; break;
  		case 1: input_gain = 2.0; break;
  		case 2: input_gain = 4.0; break;
  		case 3: input_gain = 8.0; break;
  		case 4: input_gain = 16.0; break;
  		case 5: input_gain = 32.0; break;
  		case 6: input_gain = 64.0; break;
  		case 7: input_gain = 128.0; break;
  		case 8: input_gain = 256.0; break;
  		case 9: input_gain = 512.0; break;
  		case 10: input_gain = 1024.0; break;
  		case 11: input_gain = 2048.0; break;
  		case 12: input_gain = 4096.0; break;
  		case 13: input_gain = 8192.0; break;
  		case 14: input_gain = 16384.0; break;
  		case 15: input_gain = 32768.0; break;
  		case 16: input_gain = 65536.0; break;
  	}
  	//we are directly punching in the gain values rather than calculating them

    // ClearCoat

  	clearcoat = (int)(K*16.999);
	
  	if (clearcoat != prevclearcoat) {
  		for(int count = 0; count < kshortA+2; count++) {aAL[count] = 0.0; aAR[count] = 0.0;}
  		for(int count = 0; count < kshortB+2; count++) {aBL[count] = 0.0; aBR[count] = 0.0;}
  		for(int count = 0; count < kshortC+2; count++) {aCL[count] = 0.0; aCR[count] = 0.0;}
  		for(int count = 0; count < kshortD+2; count++) {aDL[count] = 0.0; aDR[count] = 0.0;}
  		for(int count = 0; count < kshortE+2; count++) {aEL[count] = 0.0; aER[count] = 0.0;}
  		for(int count = 0; count < kshortF+2; count++) {aFL[count] = 0.0; aFR[count] = 0.0;}
  		for(int count = 0; count < kshortG+2; count++) {aGL[count] = 0.0; aGR[count] = 0.0;}
  		for(int count = 0; count < kshortH+2; count++) {aHL[count] = 0.0; aHR[count] = 0.0;}
  		for(int count = 0; count < kshortI+2; count++) {aIL[count] = 0.0; aIR[count] = 0.0;}
  		for(int count = 0; count < kshortJ+2; count++) {aJL[count] = 0.0; aJR[count] = 0.0;}
  		for(int count = 0; count < kshortK+2; count++) {aKL[count] = 0.0; aKR[count] = 0.0;}
  		for(int count = 0; count < kshortL+2; count++) {aLL[count] = 0.0; aLR[count] = 0.0;}
  		for(int count = 0; count < kshortM+2; count++) {aML[count] = 0.0; aMR[count] = 0.0;}
  		for(int count = 0; count < kshortN+2; count++) {aNL[count] = 0.0; aNR[count] = 0.0;}
  		for(int count = 0; count < kshortO+2; count++) {aOL[count] = 0.0; aOR[count] = 0.0;}
  		for(int count = 0; count < kshortP+2; count++) {aPL[count] = 0.0; aPR[count] = 0.0;}		
  		countAL = 1;
  		countBL = 1;
  		countCL = 1;
  		countDL = 1;	
  		countEL = 1;
  		countFL = 1;
  		countGL = 1;
  		countHL = 1;
  		countIL = 1;
  		countJL = 1;
  		countKL = 1;
  		countLL = 1;
  		countML = 1;
  		countNL = 1;
  		countOL = 1;
  		countPL = 1;
		
  		countAR = 1;
  		countBR = 1;
  		countCR = 1;
  		countDR = 1;	
  		countER = 1;
  		countFR = 1;
  		countGR = 1;
  		countHR = 1;
  		countIR = 1;
  		countJR = 1;
  		countKR = 1;
  		countLR = 1;
  		countMR = 1;
  		countNR = 1;
  		countOR = 1;
  		countPR = 1;
  		switch (clearcoat)
  		{
  			case 0:
  				shortA = 65; shortB = 124; shortC = 83; shortD = 180; shortE = 200; shortF = 291; shortG = 108; shortH = 189; shortI = 73; shortJ = 410; shortK = 479; shortL = 310; shortM = 11; shortN = 928; shortO = 23; shortP = 654; break; //5 to 51 ms, 96 seat room. Scarcity, 1 in 125324
  				//Short96
  			case 1:
  				shortA = 114; shortB = 205; shortC = 498; shortD = 195; shortE = 205; shortF = 318; shortG = 143; shortH = 254; shortI = 64; shortJ = 721; shortK = 512; shortL = 324; shortM = 11; shortN = 782; shortO = 26; shortP = 394; break; //7 to 52 ms, 107 seat club. Scarcity, 1 in 65537
  				//Short107
  			case 2:
  				shortA = 118; shortB = 272; shortC = 292; shortD = 145; shortE = 200; shortF = 241; shortG = 204; shortH = 504; shortI = 50; shortJ = 678; shortK = 424; shortL = 412; shortM = 11; shortN = 1124; shortO = 47; shortP = 766; break; //8 to 58 ms, 135 seat club. Scarcity, 1 in 196272
  				//Short135
  			case 3:
  				shortA = 19; shortB = 474; shortC = 301; shortD = 275; shortE = 260; shortF = 321; shortG = 371; shortH = 571; shortI = 50; shortJ = 410; shortK = 697; shortL = 414; shortM = 11; shortN = 986; shortO = 47; shortP = 522; break; //7 to 61 ms, 143 seat club. Scarcity, 1 in 165738
  				//Short143
  			case 4:
  				shortA = 112; shortB = 387; shortC = 452; shortD = 289; shortE = 173; shortF = 476; shortG = 321; shortH = 593; shortI = 73; shortJ = 343; shortK = 829; shortL = 91; shortM = 11; shortN = 1055; shortO = 43; shortP = 862; break; //8 to 66 ms, 166 seat club. Scarcity, 1 in 158437
  				//Short166
  			case 5:
  				shortA = 60; shortB = 368; shortC = 295; shortD = 272; shortE = 210; shortF = 284; shortG = 326; shortH = 830; shortI = 125; shortJ = 236; shortK = 737; shortL = 486; shortM = 11; shortN = 1178; shortO = 75; shortP = 902; break; //9 to 70 ms, 189 seat club. Scarcity, 1 in 94790
  				//Short189
  			case 6:
  				shortA = 73; shortB = 311; shortC = 472; shortD = 251; shortE = 134; shortF = 509; shortG = 393; shortH = 591; shortI = 124; shortJ = 1070; shortK = 340; shortL = 525; shortM = 11; shortN = 1367; shortO = 75; shortP = 816; break; //7 to 79 ms, 225 seat club. Scarcity, 1 in 257803
  				//Short225
  			case 7:
  				shortA = 159; shortB = 518; shortC = 514; shortD = 165; shortE = 275; shortF = 494; shortG = 296; shortH = 667; shortI = 75; shortJ = 1101; shortK = 116; shortL = 414; shortM = 11; shortN = 1261; shortO = 79; shortP = 998; break; //11 to 80 ms, 252 seat club. Scarcity, 1 in 175192
  				//Short252
  			case 8:
  				shortA = 41; shortB = 741; shortC = 274; shortD = 59; shortE = 306; shortF = 332; shortG = 291; shortH = 767; shortI = 42; shortJ = 881; shortK = 959; shortL = 422; shortM = 11; shortN = 1237; shortO = 45; shortP = 958; break; //8 to 83 ms, 255 seat club. Scarcity, 1 in 185708
  				//Short255
  			case 9:
  				shortA = 251; shortB = 437; shortC = 783; shortD = 189; shortE = 130; shortF = 272; shortG = 244; shortH = 761; shortI = 128; shortJ = 1190; shortK = 320; shortL = 491; shortM = 11; shortN = 1409; shortO = 58; shortP = 455; break; //10 to 93 ms, 323 seat club. Scarcity, 1 in 334044
  				//Short323
  			case 10:
  				shortA = 316; shortB = 510; shortC = 1087; shortD = 349; shortE = 359; shortF = 74; shortG = 79; shortH = 1269; shortI = 34; shortJ = 693; shortK = 749; shortL = 511; shortM = 11; shortN = 1751; shortO = 93; shortP = 403; break; //9 to 110 ms, 427 seat theater. Scarcity, 1 in 200715
  				//Short427
  			case 11:
  				shortA = 254; shortB = 651; shortC = 845; shortD = 316; shortE = 373; shortF = 267; shortG = 182; shortH = 857; shortI = 215; shortJ = 1535; shortK = 1127; shortL = 315; shortM = 11; shortN = 1649; shortO = 97; shortP = 829; break; //15 to 110 ms, 470 seat theater. Scarcity, 1 in 362673
  				//Short470
  			case 12:
  				shortA = 113; shortB = 101; shortC = 673; shortD = 357; shortE = 340; shortF = 229; shortG = 278; shortH = 1008; shortI = 265; shortJ = 1890; shortK = 155; shortL = 267; shortM = 11; shortN = 2233; shortO = 116; shortP = 600; break; //11 to 131 ms, 606 seat theater. Scarcity, 1 in 238058
  				//Short606
  			case 13:
  				shortA = 218; shortB = 1058; shortC = 862; shortD = 505; shortE = 297; shortF = 580; shortG = 532; shortH = 1387; shortI = 120; shortJ = 576; shortK = 1409; shortL = 473; shortM = 11; shortN = 1991; shortO = 76; shortP = 685; break; //14 to 132 ms, 643 seat theater. Scarcity, 1 in 193432
  				//Short643
  			case 14:
  				shortA = 78; shortB = 760; shortC = 982; shortD = 528; shortE = 445; shortF = 1128; shortG = 130; shortH = 708; shortI = 22; shortJ = 2144; shortK = 354; shortL = 1169; shortM = 11; shortN = 2782; shortO = 58; shortP = 1515; break; //5 to 159 ms, 809 seat hall. Scarcity, 1 in 212274
  				//Short809
  			case 15:
  				shortA = 330; shortB = 107; shortC = 1110; shortD = 371; shortE = 620; shortF = 143; shortG = 1014; shortH = 1763; shortI = 184; shortJ = 2068; shortK = 1406; shortL = 595; shortM = 11; shortN = 2639; shortO = 33; shortP = 1594; break; //10 to 171 ms, 984 seat hall. Scarcity, 1 in 226499
  				//Short984
  			case 16:
  			default:
  				shortA = 336; shortB = 1660; shortC = 386; shortD = 623; shortE = 693; shortF = 1079; shortG = 891; shortH = 1574; shortI = 24; shortJ = 2641; shortK = 1239; shortL = 775; shortM = 11; shortN = 3104; shortO = 55; shortP = 2366; break; //24 to 203 ms, 1541 seat hall. Scarcity, 1 in 275025
  				//Short1541
  		}
  		prevclearcoat = clearcoat;
  	}

    // Channel 9

  	switch((VstInt32)( Z * 4.999 ))
  	{  
  		case 0: ch9_iirAmount = 0.005832; ch9_threshold = 0.33362176; ch9_cutoff = 28811.0; break; //Neve
  		case 1: ch9_iirAmount = 0.004096; ch9_threshold = 0.59969536; ch9_cutoff = 27216.0; break; //API
  		case 2: ch9_iirAmount = 0.004913; ch9_threshold = 0.84934656; ch9_cutoff = 23011.0; break; //SSL
  		case 3: ch9_iirAmount = 0.009216; ch9_threshold = 0.149; ch9_cutoff = 18544.0; break; //Teac
  		case 4: ch9_iirAmount = 0.011449; ch9_threshold = 0.092; ch9_cutoff = 19748.0; break; //Mackie
  		default: break; //should not happen
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

void ConsoleZTracking::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Input", kVstMaxParamStrLen); break;
    		case kParamB: vst_strncpy (text, "HPF Freq", kVstMaxParamStrLen); break;
    		case kParamC: vst_strncpy (text, "HPF Q", kVstMaxParamStrLen); break;
    		case kParamD: vst_strncpy (text, "LPF Freq", kVstMaxParamStrLen); break;
    		case kParamE: vst_strncpy (text, "LPF Q", kVstMaxParamStrLen); break;
        case kParamF: vst_strncpy (text, "Compress", kVstMaxParamStrLen); break;
    		case kParamG: vst_strncpy (text, "CompOut", kVstMaxParamStrLen); break;
    		case kParamH: vst_strncpy (text, "CompD/W", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "ChSpeed", kVstMaxParamStrLen); break;
        case kParamJ: vst_strncpy (text, "ChDepth", kVstMaxParamStrLen); break;
        case kParamK: vst_strncpy (text, "RevSelect", kVstMaxParamStrLen); break;
    		case kParamL: vst_strncpy (text, "RevD/W", kVstMaxParamStrLen); break;
        case kParamM: vst_strncpy (text, "Boost", kVstMaxParamStrLen); break;
    		case kParamN: vst_strncpy (text, "Focus", kVstMaxParamStrLen); break;
    		case kParamO: vst_strncpy (text, "CenterF", kVstMaxParamStrLen); break;
    		case kParamP: vst_strncpy (text, "Mode", kVstMaxParamStrLen); break;
        case kParamQ: vst_strncpy (text, "uLaw", kVstMaxParamStrLen); break;
    		case kParamR: vst_strncpy (text, "Output", kVstMaxParamStrLen); break;
    		case kParamS: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        case kParamT: vst_strncpy (text, "Drive", kVstMaxParamStrLen); break;
    		case kParamU: vst_strncpy (text, "Depth", kVstMaxParamStrLen); break;
    		case kParamV: vst_strncpy (text, "Inv/Wet", kVstMaxParamStrLen); break;
        case kParamW: vst_strncpy (text, "Voicing", kVstMaxParamStrLen); break;
    		case kParamX: vst_strncpy (text, "Highs", kVstMaxParamStrLen); break;
    		case kParamY: vst_strncpy (text, "Trim", kVstMaxParamStrLen); break;
    		case kParamZ: vst_strncpy (text, "Console", kVstMaxParamStrLen); break;
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
        case kParamF: float2string (F, text, kVstMaxParamStrLen); break;
        case kParamG: float2string (G * 2.0, text, kVstMaxParamStrLen); break;
        case kParamH: float2string (H, text, kVstMaxParamStrLen); break;
        case kParamI: float2string (I, text, kVstMaxParamStrLen); break;
        case kParamJ: int2string (J, text, kVstMaxParamStrLen); break;
        case kParamK: int2string ((VstInt32)( K * 16.999 ), text, kVstMaxParamStrLen); break;
        case kParamL: float2string (L, text, kVstMaxParamStrLen); break;
        case kParamM: float2string (M*12.0, text, kVstMaxParamStrLen); break;
        case kParamN: float2string (N, text, kVstMaxParamStrLen); break;
        case kParamO: float2string (pow(M_E, freqX * O) * centerFreq, text, kVstMaxParamStrLen); break;
        case kParamP: switch((VstInt32)( P * 4.999 )) //0 to almost edge of # of params
    		{
    			case 0: vst_strncpy (text, "Density", kVstMaxParamStrLen); break;
    			case 1: vst_strncpy (text, "Drive", kVstMaxParamStrLen); break;
    			case 2: vst_strncpy (text, "Spiral", kVstMaxParamStrLen); break;
    			case 3: vst_strncpy (text, "Mojo", kVstMaxParamStrLen); break;
    			case 4: vst_strncpy (text, "Dyno", kVstMaxParamStrLen); break;
    			default: break; // unknown parameter, shouldn't happen!
    		} break;
        case kParamQ: if(Q == 0.0) {vst_strncpy (text, "Off", kVstMaxParamStrLen);} else {vst_strncpy (text, "On", kVstMaxParamStrLen);} break;
        case kParamS: float2string (S, text, kVstMaxParamStrLen); break;
        case kParamT: float2string (T, text, kVstMaxParamStrLen); break;
        case kParamU: float2string (U, text, kVstMaxParamStrLen); break;
        case kParamV: float2string (V, text, kVstMaxParamStrLen); break;
        case kParamW: float2string ((W*99.0)+1.0, text, kVstMaxParamStrLen); break;
        case kParamX: float2string ((X*2.0)-1.0, text, kVstMaxParamStrLen); break;
        case kParamY: float2string (pow(M_E, freqX * Y) * 0.25, text, kVstMaxParamStrLen); break;
        case kParamZ: switch((VstInt32)( Z * 4.999 )) //0 to almost edge of # of params
    		{	case 0: vst_strncpy (text, "Neve", kVstMaxParamStrLen); break;
    			case 1: vst_strncpy (text, "API", kVstMaxParamStrLen); break;
    			case 2: vst_strncpy (text, "SSL", kVstMaxParamStrLen); break;
    			case 3: vst_strncpy (text, "Teac", kVstMaxParamStrLen); break;
    			case 4: vst_strncpy (text, "Mackie", kVstMaxParamStrLen); break;
    			default: break; // unknown parameter, shouldn't happen!
    		} break; //completed consoletype 'popup' parameter, exit
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
        case kParamF: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamG: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamH: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamI: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamJ: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamK: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamL: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamM: vst_strncpy (text, "dB", kVstMaxParamStrLen); break;
        case kParamN: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamO: vst_strncpy (text, "Hz", kVstMaxParamStrLen); break;
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
