/* ========================================
 *  ClearCoatV - ClearCoatV.h
 *  Copyright (c) airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ClearCoatV_H
#include "ClearCoatV.h"
#endif

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {return new ClearCoatV(audioMaster);}

ClearCoatV::ClearCoatV(audioMasterCallback audioMaster) :
    AudioEffectX(audioMaster, kNumPrograms, kNumParameters)
{
	A = 0.5;
	B = 1.0;
	
	for(int count = 0; count < 8192; count++) {
		aA[count] = 0.0;
		aB[count] = 0.0;
		aC[count] = 0.0;
		aD[count] = 0.0;
		aE[count] = 0.0;
		aF[count] = 0.0;
		aG[count] = 0.0;
		aH[count] = 0.0;
		aI[count] = 0.0;
		aJ[count] = 0.0;
		aK[count] = 0.0;
		aL[count] = 0.0;
		aM[count] = 0.0;
		aN[count] = 0.0;
		aO[count] = 0.0;
		aP[count] = 0.0;
	}
	
  for(int count = 0; count < 2; count++) {
  	feedbackA[count] = 0.0;
  	feedbackB[count] = 0.0;
  	feedbackC[count] = 0.0;
  	feedbackD[count] = 0.0;
  }
	
	for(int count = 0; count < 34; count++) {lastRef[count] = 0.0;}

	prevMulchA[0] = 0.0;
	prevMulchA[1] = 0.0;
	cycle = 0;
	gcount = 0;
	
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

ClearCoatV::~ClearCoatV() {}
VstInt32 ClearCoatV::getVendorVersion () {return 1000;}
void ClearCoatV::setProgramName(char *name) {vst_strncpy (_programName, name, kVstMaxProgNameLen);}
void ClearCoatV::getProgramName(char *name) {vst_strncpy (name, _programName, kVstMaxProgNameLen);}
//airwindows likes to ignore this stuff. Make your own programs, and make a different plugin rather than
//trying to do versioning and preventing people from using older versions. Maybe they like the old one!

static float pinParameter(float data)
{
	if (data < 0.0f) return 0.0f;
	if (data > 1.0f) return 1.0f;
	return data;
}

VstInt32 ClearCoatV::getChunk (void** data, bool isPreset)
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

VstInt32 ClearCoatV::setChunk (void* data, VstInt32 byteSize, bool isPreset)
{	
	float *chunkData = (float *)data;
	A = pinParameter(chunkData[0]);
	B = pinParameter(chunkData[1]);
	/* We're ignoring byteSize as we found it to be a filthy liar */
	
	/* calculate any other fields you need here - you could copy in 
	 code from setParameter() here. */
	return 0;
}

void ClearCoatV::setParameter(VstInt32 index, float value) {
    switch (index) {
        case kParamA: A = value; break;
        case kParamB: B = value; break;
        default: throw; // unknown parameter, shouldn't happen!
    }

	int clearcoat = (int)(A*16.999);
	
	if (clearcoat != prevclearcoat) {
		for(int count = 0; count < 8192; count++) {
			aA[count] = 0.0;
			aB[count] = 0.0;
			aC[count] = 0.0;
			aD[count] = 0.0;
			aE[count] = 0.0;
			aF[count] = 0.0;
			aG[count] = 0.0;
			aH[count] = 0.0;
			aI[count] = 0.0;
			aJ[count] = 0.0;
			aK[count] = 0.0;
			aL[count] = 0.0;
			aM[count] = 0.0;
			aN[count] = 0.0;
			aO[count] = 0.0;
			aP[count] = 0.0;
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
}

float ClearCoatV::getParameter(VstInt32 index) {
    switch (index) {
        case kParamA: return A; break;
        case kParamB: return B; break;
        default: break; // unknown parameter, shouldn't happen!
    } return 0.0; //we only need to update the relevant name, this is simple to manage
}

void ClearCoatV::getParameterName(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "Select", kVstMaxParamStrLen); break;
		case kParamB: vst_strncpy (text, "Dry/Wet", kVstMaxParamStrLen); break;
        default: break; // unknown parameter, shouldn't happen!
    } //this is our labels for displaying in the VST host
}

void ClearCoatV::getParameterDisplay(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: int2string ((VstInt32)( A * 16.999 ), text, kVstMaxParamStrLen); break;
        case kParamB: float2string (B, text, kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
	} //this displays the values and handles 'popups' where it's discrete choices
}

void ClearCoatV::getParameterLabel(VstInt32 index, char *text) {
    switch (index) {
        case kParamA: vst_strncpy (text, "", kVstMaxParamStrLen); break;
        case kParamB: vst_strncpy (text, "", kVstMaxParamStrLen); break;
		default: break; // unknown parameter, shouldn't happen!
    }
}

VstInt32 ClearCoatV::canDo(char *text) 
{ return (_canDo.find(text) == _canDo.end()) ? -1: 1; } // 1 = yes, -1 = no, 0 = don't know

bool ClearCoatV::getEffectName(char* name) {
    vst_strncpy(name, "ClearCoatV3", kVstMaxProductStrLen); return true;
}

VstPlugCategory ClearCoatV::getPlugCategory() {return kPlugCategEffect;}

bool ClearCoatV::getProductString(char* text) {
  	vst_strncpy (text, "airwindows ClearCoatV3", kVstMaxProductStrLen); return true;
}

bool ClearCoatV::getVendorString(char* text) {
  	vst_strncpy (text, "airwindows", kVstMaxVendorStrLen); return true;
}
