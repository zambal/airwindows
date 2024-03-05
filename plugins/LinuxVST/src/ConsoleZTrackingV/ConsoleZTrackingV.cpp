/* ========================================
 *  ConsoleZTrackingV - ConsoleZTrackingV.cpp
 *  Copyright (c) 2016 airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZTrackingV_H
#include "ConsoleZTrackingV.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ConsoleZTrackingV(audioMaster);}

ConsoleZTrackingV::ConsoleZTrackingV(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.5;
	B = 0.0;
	C = 0.0;
	D = 0.0;
	E = 0.0;
	F = 0.0;
	G = 0.5;
	H = 0.0;
	I = 0.0;
	J = 0.0;
	K = 0.0;
	L = 0.0;
	M = 0.0;
	N = 0.5;
	O = 0.5;
	P = 0.5;
	Q = 0.0;
	R = 0.0;
	S = 0.0;
	T = 0.0;
	U = 0.0;
	V = 0.5;
	W = 0.0;
	X = 0.5;
	Y = 0.5;
	Z = 0.0;


	// BitshiftGain

	input_gain = 1.0;

  // BiquadOneHalf HPF

	for (int x = 0; x < 8; x++) {biquad_hpf[x] = 0.0;}

  // BiquadOneHalf LPF

	for (int x = 0; x < 8; x++) {biquad_lpf[x] = 0.0;}

  // Buttercomp2

  for(int i=0; i<4; i++) {
  	controlpos_b[i] = 1.0;
  	controlneg_b[i] = 1.0;
  	targetpos_b[i] = 1.0;
  	targetneg_b[i] = 1.0;
  	lastOutput_b[i] = 0.0;
  }

  // StereoChorus

	for(int count = 0; count < 65536; count++) {p_b[count] = 0;}
	sweep_b[0] = 3.141592653589793238 / 2.7;
	sweep_b[1] = 3.141592653589793238;
	gcount = 0;
	airEven_b[0] = airEven_b[1] = 0.0;
	airOdd_b[0] = airOdd_b[1] = 0.0;
	airPrev_b[0] = airPrev_b[1] = 0.0;
	
	for(int count = 0; count < 20; count++) {lastRef[count] = 0.0;}
	chorus_cycle = 0;
  chorus_flip = false;

  // ClearCoat

	for(int count = 0; count < 8192; count++) {
  	aA[count & 2047] = 0.0;
  	aB[count & 4095] = 0.0;
  	aC[count & 4095] = 0.0;
  	aD[count & 8191] = 0.0;
  	aE[count & 4095] = 0.0;
  	aF[count & 4095] = 0.0;
  	aG[count & 4095] = 0.0;
  	aH[count & 4095] = 0.0;
  	aI[count & 4095] = 0.0;
  	aJ[count & 8191] = 0.0;
  	aK[count & 8191] = 0.0;
  	aL[count & 8191] = 0.0;
  	aM[count & 2047] = 0.0;
  	aN[count & 8191] = 0.0;
  	aO[count & 2047] = 0.0;
  	aP[count & 8191] = 0.0;
	}
	
  for(int count = 0; count < 2; count++) {
  	feedbackA[count] = 0.0;
  	feedbackB[count] = 0.0;
  	feedbackC[count] = 0.0;
  	feedbackD[count] = 0.0;
  }
	
	for(int count = 0; count < 8; count++) {cc_lastRef[count] = 0.0;}

	prevMulchA[0] = 0.0;
	prevMulchA[1] = 0.0;

  cc_cycle = 0;
  cc_gcount = 0;
	
	shortA = 2 * 336;
	shortB = 2 * 1660;
	shortC = 2 * 386;
	shortD = 2 * 623;
	shortE = 2 * 693;
	shortF = 2 * 1079;
	shortG = 2 * 891;
	shortH = 2 * 1574;
	shortI = 2 * 24;
	shortJ = 2 * 2641;
	shortK = 2 * 1239;
	shortL = 2 * 775;
	shortM = 2 * 11;
	shortN = 2 * 3104;
	shortO = 2 * 55;
	shortP = 2 * 2366;
	prevclearcoat = -1;
	
	tail[0] = tail[1] = 0.0;
	for(int count = 0; count < 16; count++) { cc_sub_b[count] = 0.0; }

  // Channel9

	for (int x = 0; x < 16; x++) {ch9_biquad[x] = 0.0;}
  ch9_iirSamples[0] = 0.0; ch9_iirSamples[1] = 0.0; ch9_iirSamples[2] = 0.0; ch9_iirSamples[3] = 0.0;
	ch9_lastSample[0] = ch9_lastSample[1] = ch9_lastSample[2] = ch9_lastSample[3] = ch9_lastSample[4] = ch9_lastSample[5] = 0.0;
	ch9_iirAmount = 0.005832;
	ch9_threshold = 0.33362176; //instantiating with Neve values
	ch9_cutoff = 28811.0;
  
  // Focus

	for (int x = 0; x < 4; x++) {figure[x] = 0.0;}
  freqX = 2.0 * log((16.0 - 4.0) / (4.0 - 1.0));

  // Creature

	for (int x = 0; x < 256; x++) {
		cr_slew[x] = 0.0;
	}

  // ToneSlant

	for(int count = 0; count < 128; count++) {ts_f[count] = 0.0;}
	for(int count = 0; count < 260; count++) {ts_b[count] = 0.0;}
  ts_gcount = 0;

  adcLastSample[0] = 0.0;
  adcLastSample[1] = 0.0;

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

ConsoleZTrackingV::~ConsoleZTrackingV() {}
VstInt32 ConsoleZTrackingV::getVendorVersion () {return 1000;}
void ConsoleZTrackingV::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ConsoleZTrackingV::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ConsoleZTrackingV::getChunk (void** data, bool isPreset)
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

VstInt32 ConsoleZTrackingV::setChunk (void* data, VstInt32 byteSize, bool isPreset)
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

void ConsoleZTrackingV::setParameter(VstInt32 index, float value) {
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

  	int clearcoat = (int)(K*16.999);
	
  	if (clearcoat != prevclearcoat) {
      for(int count = 0; count < 8192; count++) {
      	aA[count & 2047] = 0.0;
      	aB[count & 4095] = 0.0;
      	aC[count & 4095] = 0.0;
      	aD[count & 8191] = 0.0;
      	aE[count & 4095] = 0.0;
      	aF[count & 4095] = 0.0;
      	aG[count & 4095] = 0.0;
      	aH[count & 4095] = 0.0;
      	aI[count & 4095] = 0.0;
      	aJ[count & 8191] = 0.0;
      	aK[count & 8191] = 0.0;
      	aL[count & 8191] = 0.0;
      	aM[count & 2047] = 0.0;
      	aN[count & 8191] = 0.0;
      	aO[count & 2047] = 0.0;
      	aP[count & 8191] = 0.0;
      }

  		switch (clearcoat)
  		{
  			case 0:
  				shortA = 2 * 65; shortB = 2 * 124; shortC = 2 * 83; shortD = 2 * 180; shortE = 2 * 200; shortF = 2 * 291; shortG = 2 * 108; shortH = 2 * 189; shortI = 2 * 73; shortJ = 2 * 410; shortK = 2 * 479; shortL = 2 * 310; shortM = 2 * 11; shortN = 2 * 928; shortO = 2 * 23; shortP = 2 * 654; break; //5 to 51 ms, 96 seat room. Scarcity, 1 in 125324
  				//Short96
  			case 1:
  				shortA = 2 * 114; shortB = 2 * 205; shortC = 2 * 498; shortD = 2 * 195; shortE = 2 * 205; shortF = 2 * 318; shortG = 2 * 143; shortH = 2 * 254; shortI = 2 * 64; shortJ = 2 * 721; shortK = 2 * 512; shortL = 2 * 324; shortM = 2 * 11; shortN = 2 * 782; shortO = 2 * 26; shortP = 2 * 394; break; //7 to 52 ms, 107 seat club. Scarcity, 1 in 65537
  				//Short107
  			case 2:
  				shortA = 2 * 118; shortB = 2 * 272; shortC = 2 * 292; shortD = 2 * 145; shortE = 2 * 200; shortF = 2 * 241; shortG = 2 * 204; shortH = 2 * 504; shortI = 2 * 50; shortJ = 2 * 678; shortK = 2 * 424; shortL = 2 * 412; shortM = 2 * 11; shortN = 2 * 1124; shortO = 2 * 47; shortP = 2 * 766; break; //8 to 58 ms, 135 seat club. Scarcity, 1 in 196272
  				//Short135
  			case 3:
  				shortA = 2 * 19; shortB = 2 * 474; shortC = 2 * 301; shortD = 2 * 275; shortE = 2 * 260; shortF = 2 * 321; shortG = 2 * 371; shortH = 2 * 571; shortI = 2 * 50; shortJ = 2 * 410; shortK = 2 * 697; shortL = 2 * 414; shortM = 2 * 11; shortN = 2 * 986; shortO = 2 * 47; shortP = 2 * 522; break; //7 to 61 ms, 143 seat club. Scarcity, 1 in 165738
  				//Short143
  			case 4:
  				shortA = 2 * 112; shortB = 2 * 387; shortC = 2 * 452; shortD = 2 * 289; shortE = 2 * 173; shortF = 2 * 476; shortG = 2 * 321; shortH = 2 * 593; shortI = 2 * 73; shortJ = 2 * 343; shortK = 2 * 829; shortL = 2 * 91; shortM = 2 * 11; shortN = 2 * 1055; shortO = 2 * 43; shortP = 2 * 862; break; //8 to 66 ms, 166 seat club. Scarcity, 1 in 158437
  				//Short166
  			case 5:
  				shortA = 2 * 60; shortB = 2 * 368; shortC = 2 * 295; shortD = 2 * 272; shortE = 2 * 210; shortF = 2 * 284; shortG = 2 * 326; shortH = 2 * 830; shortI = 2 * 125; shortJ = 2 * 236; shortK = 2 * 737; shortL = 2 * 486; shortM = 2 * 11; shortN = 2 * 1178; shortO = 2 * 75; shortP = 2 * 902; break; //9 to 70 ms, 189 seat club. Scarcity, 1 in 94790
  				//Short189
  			case 6:
  				shortA = 2 * 73; shortB = 2 * 311; shortC = 2 * 472; shortD = 2 * 251; shortE = 2 * 134; shortF = 2 * 509; shortG = 2 * 393; shortH = 2 * 591; shortI = 2 * 124; shortJ = 2 * 1070; shortK = 2 * 340; shortL = 2 * 525; shortM = 2 * 11; shortN = 2 * 1367; shortO = 2 * 75; shortP = 2 * 816; break; //7 to 79 ms, 225 seat club. Scarcity, 1 in 257803
  				//Short225
  			case 7:
  				shortA = 2 * 159; shortB = 2 * 518; shortC = 2 * 514; shortD = 2 * 165; shortE = 2 * 275; shortF = 2 * 494; shortG = 2 * 296; shortH = 2 * 667; shortI = 2 * 75; shortJ = 2 * 1101; shortK = 2 * 116; shortL = 2 * 414; shortM = 2 * 11; shortN = 2 * 1261; shortO = 2 * 79; shortP = 2 * 998; break; //11 to 80 ms, 252 seat club. Scarcity, 1 in 175192
  				//Short252
  			case 8:
  				shortA = 2 * 41; shortB = 2 * 741; shortC = 2 * 274; shortD = 2 * 59; shortE = 2 * 306; shortF = 2 * 332; shortG = 2 * 291; shortH = 2 * 767; shortI = 2 * 42; shortJ = 2 * 881; shortK = 2 * 959; shortL = 2 * 422; shortM = 2 * 11; shortN = 2 * 1237; shortO = 2 * 45; shortP = 2 * 958; break; //8 to 83 ms, 255 seat club. Scarcity, 1 in 185708
  				//Short255
  			case 9:
  				shortA = 2 * 251; shortB = 2 * 437; shortC = 2 * 783; shortD = 2 * 189; shortE = 2 * 130; shortF = 2 * 272; shortG = 2 * 244; shortH = 2 * 761; shortI = 2 * 128; shortJ = 2 * 1190; shortK = 2 * 320; shortL = 2 * 491; shortM = 2 * 11; shortN = 2 * 1409; shortO = 2 * 58; shortP = 2 * 455; break; //10 to 93 ms, 323 seat club. Scarcity, 1 in 334044
  				//Short323
  			case 10:
  				shortA = 2 * 316; shortB = 2 * 510; shortC = 2 * 1087; shortD = 2 * 349; shortE = 2 * 359; shortF = 2 * 74; shortG = 2 * 79; shortH = 2 * 1269; shortI = 2 * 34; shortJ = 2 * 693; shortK = 2 * 749; shortL = 2 * 511; shortM = 2 * 11; shortN = 2 * 1751; shortO = 2 * 93; shortP = 2 * 403; break; //9 to 110 ms, 427 seat theater. Scarcity, 1 in 200715
  				//Short427
  			case 11:
  				shortA = 2 * 254; shortB = 2 * 651; shortC = 2 * 845; shortD = 2 * 316; shortE = 2 * 373; shortF = 2 * 267; shortG = 2 * 182; shortH = 2 * 857; shortI = 2 * 215; shortJ = 2 * 1535; shortK = 2 * 1127; shortL = 2 * 315; shortM = 2 * 11; shortN = 2 * 1649; shortO = 2 * 97; shortP = 2 * 829; break; //15 to 110 ms, 470 seat theater. Scarcity, 1 in 362673
  				//Short470
  			case 12:
  				shortA = 2 * 113; shortB = 2 * 101; shortC = 2 * 673; shortD = 2 * 357; shortE = 2 * 340; shortF = 2 * 229; shortG = 2 * 278; shortH = 2 * 1008; shortI = 2 * 265; shortJ = 2 * 1890; shortK = 2 * 155; shortL = 2 * 267; shortM = 2 * 11; shortN = 2 * 2233; shortO = 2 * 116; shortP = 2 * 600; break; //11 to 131 ms, 606 seat theater. Scarcity, 1 in 238058
  				//Short606
  			case 13:
  				shortA = 2 * 218; shortB = 2 * 1058; shortC = 2 * 862; shortD = 2 * 505; shortE = 2 * 297; shortF = 2 * 580; shortG = 2 * 532; shortH = 2 * 1387; shortI = 2 * 120; shortJ = 2 * 576; shortK = 2 * 1409; shortL = 2 * 473; shortM = 2 * 11; shortN = 2 * 1991; shortO = 2 * 76; shortP = 2 * 685; break; //14 to 132 ms, 643 seat theater. Scarcity, 1 in 193432
  				//Short643
  			case 14:
  				shortA = 2 * 78; shortB = 2 * 760; shortC = 2 * 982; shortD = 2 * 528; shortE = 2 * 445; shortF = 2 * 1128; shortG = 2 * 130; shortH = 2 * 708; shortI = 2 * 22; shortJ = 2 * 2144; shortK = 2 * 354; shortL = 2 * 1169; shortM = 2 * 11; shortN = 2 * 2782; shortO = 2 * 58; shortP = 2 * 1515; break; //5 to 159 ms, 809 seat hall. Scarcity, 1 in 212274
  				//Short809
  			case 15:
  				shortA = 2 * 330; shortB = 2 * 107; shortC = 2 * 1110; shortD = 2 * 371; shortE = 2 * 620; shortF = 2 * 143; shortG = 2 * 1014; shortH = 2 * 1763; shortI = 2 * 184; shortJ = 2 * 2068; shortK = 2 * 1406; shortL = 2 * 595; shortM = 2 * 11; shortN = 2 * 2639; shortO = 2 * 33; shortP = 2 * 1594; break; //10 to 171 ms, 984 seat hall. Scarcity, 1 in 226499
  				//Short984
  			case 16:
  			default:
  				shortA = 2 * 336; shortB = 2 * 1660; shortC = 2 * 386; shortD = 2 * 623; shortE = 2 * 693; shortF = 2 * 1079; shortG = 2 * 891; shortH = 2 * 1574; shortI = 2 * 24; shortJ = 2 * 2641; shortK = 2 * 1239; shortL = 2 * 775; shortM = 2 * 11; shortN = 2 * 3104; shortO = 2 * 55; shortP = 2 * 2366; break; //24 to 203 ms, 1541 seat hall. Scarcity, 1 in 275025
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

float ConsoleZTrackingV::getParameter(VstInt32 index) {
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

void ConsoleZTrackingV::getParameterName(VstInt32 index, char *text) {
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

void ConsoleZTrackingV::getParameterDisplay(VstInt32 index, char *text) {
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

void ConsoleZTrackingV::getParameterLabel(VstInt32 index, char *text) {
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

VstInt32 ConsoleZTrackingV::canDo(char *text)
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ConsoleZTrackingV::getEffectName(char* name) {
    vst_strncpy(name, "ConsoleZTrackingV", kVstMaxProductStrLen); return true;
}

VstPlugCategory ConsoleZTrackingV::getPlugCategory() {return kPlugCategEffect;}

bool ConsoleZTrackingV::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ConsoleZTrackingV", kVstMaxProductStrLen); return true;
}

bool ConsoleZTrackingV::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
