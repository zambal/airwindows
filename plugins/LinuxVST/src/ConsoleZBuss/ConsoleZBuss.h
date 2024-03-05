/* ========================================
 *  ConsoleZBuss - ConsoleZBuss.h
 *  Created 8/12/11 by SPIAdmin
 *  Copyright (c) Airwindows, Airwindows uses the MIT license
 * ======================================== */

#ifndef __ConsoleZBuss_H
#define __ConsoleZBuss_H

#ifndef __audioeffect__
#include "audioeffectx.h"
#endif

#include <set>
#include <string>
#include <math.h>

enum {
	kParamA = 0,
	kParamB = 1,
	kParamC = 2,
	kParamD = 3,
	kParamE = 4,
	kParamF = 5,
	kParamG = 6,
	kParamH = 7,
	kParamI = 8,
	kParamJ = 9,
	kParamK = 10,
	kParamL = 11,
	kParamM = 12,
	kParamN = 13,
	kParamO = 14,
	kParamP = 15,
	kParamQ = 16,
	kParamR = 17,
	kParamS = 18,
	kParamT = 19,
	kParamU = 20,
	kParamV = 21,
	kParamW = 22,
	kParamX = 23,
	kParamY = 24,
	kParamZ = 25,
	kParamZ0 = 26,
	kParamZ1 = 27,
	kParamZ2 = 28,
	kParamZ3 = 29,
	kParamZ4 = 30,
  kNumParameters = 31
}; //

enum {
    kDKAD,
    kDKCD,
    kPEAK,
    kSLEW,
    kSUBS,
    kMONO,
    kSIDE,
    kVINYL,
    kAURAT,
    kMONORAT,
    kMONOLAT,
    kPHONE,
    kCANSA,
    kCANSB,
    kCANSC,
    kCANSD,
    kTRICK,
};

const int kNumPrograms = 0;
const int kNumInputs = 2;
const int kNumOutputs = 2;
const unsigned long kUniqueId = 'czb';    //Change this to what the AU identity is!

class ConsoleZBuss :
    public AudioEffectX
{
public:
    ConsoleZBuss(audioMasterCallback audioMaster);
    ~ConsoleZBuss();
    virtual bool getEffectName(char* name);                       // The plug-in name
    virtual VstPlugCategory getPlugCategory();                    // The general category for the plug-in
    virtual bool getProductString(char* text);                    // This is a unique plug-in string provided by Steinberg
    virtual bool getVendorString(char* text);                     // Vendor info
    virtual VstInt32 getVendorVersion();                          // Version number
    virtual void processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames);
    virtual void processDoubleReplacing (double** inputs, double** outputs, VstInt32 sampleFrames);
    virtual void getProgramName(char *name);                      // read the name from the host
    virtual void setProgramName(char *name);                      // changes the name of the preset displayed in the host
	virtual VstInt32 getChunk (void** data, bool isPreset);
	virtual VstInt32 setChunk (void* data, VstInt32 byteSize, bool isPreset);
    virtual float getParameter(VstInt32 index);                   // get the parameter value at the specified index
    virtual void setParameter(VstInt32 index, float value);       // set the parameter at index to value
    virtual void getParameterLabel(VstInt32 index, char *text);  // label for the parameter (eg dB)
    virtual void getParameterName(VstInt32 index, char *text);    // name of the parameter
    virtual void getParameterDisplay(VstInt32 index, char *text); // text description of the current value
    virtual VstInt32 canDo(char *text);
private:
    char _programName[kVstMaxProgNameLen + 1];
    std::set< std::string > _canDo;

	uint32_t fpdL;
	uint32_t fpdR;
	//default stuff

  // Hypersonic

	enum {
		fix_freq,
		fix_reso,
		fix_a0,
		fix_a1,
		fix_a2,
		fix_b1,
		fix_b2,
		fix_sL1,
		fix_sL2,
		fix_sR1,
		fix_sR2,
		fix_total
	}; //fixed frequency biquad filter for ultrasonics, stereo

	double fixB[fix_total];
	double fixC[fix_total];
	double fixF[fix_total];

  // ConsoleLABuss
	double lastSinewL;
	double lastSinewR;

	double subAL;
	double subAR;
	double subBL;
	double subBR;
	double subCL;
	double subCR;

	double gainA;
	double gainB; //smoothed master fader for channel, from Z2 series filter code

  // Desk4
	double desk4_dL[10000];
	double desk4_controlL;
	double lastSampleL;
	double lastOutSampleL;
	double lastSlewL;

	double desk4_dR[10000];
	double desk4_controlR;
	double lastSampleR;
	double lastOutSampleR;
	double lastSlewR;

	int desk4_gcount;

  // Buttercomp2

	double controlAposL;
	double controlAnegL;
	double controlBposL;
	double controlBnegL;
	double targetposL;
	double targetnegL;
	double lastOutputL;
	double controlAposR;
	double controlAnegR;
	double controlBposR;
	double controlBnegR;
	double targetposR;
	double targetnegR;
	double lastOutputR;
	bool butter_flip;

  // VariMu

	double muVaryL;
	double muAttackL;
	double muNewSpeedL;
	double muSpeedAL;
	double muSpeedBL;
	double muCoefficientAL;
	double muCoefficientBL;
	double previousL;

	double muVaryR;
	double muAttackR;
	double muNewSpeedR;
	double muSpeedAR;
	double muSpeedBR;
	double muCoefficientAR;
	double muCoefficientBR;
	double previousR;
	bool flip;

  // BiquadOneHalf 1

	double biquad_1_AL[9];
	double biquad_1_AR[9];
	double biquad_1_BL[9];
	double biquad_1_BR[9];

  // BiquadOneHalf 2

	double biquad_2_AL[9];
	double biquad_2_AR[9];
	double biquad_2_BL[9];
	double biquad_2_BR[9];

  // BussColors4

	double bL[175]; //full buffer for high sample rates. Scales to 192K
	double bR[175]; //full buffer for high sample rates. Scales to 192K
	double dL[100]; //buffer for calculating sag as it relates to the dynamic impulse synthesis. To 192K.
	double dR[100]; //buffer for calculating sag as it relates to the dynamic impulse synthesis. To 192K.
	int c[35]; //just the number of taps we use, doesn't have to scale
	double g[9]; //console model
	double outg[9]; //console model
	double controlL;
	double controlR;
	double slowdynL;
	double slowdynR;
	int gcount;

  // ADClip7

	double ad_lastSampleL;
	double ad_lastSampleR;
	float ad_bL[22200];
	float ad_bR[22200];
	int ad_gcount;
	double lowsL;
	double lowsR;
	double iirLowsAL;
	double iirLowsAR;
	double iirLowsBL;
	double iirLowsBR;
	double refclipL;
	double refclipR;

  // Monitoring2

	double biquad[fix_total];
	//Bandpasses

	float darkSampleL[100];
	float darkSampleR[100];

	double m_aL[1503], m_bL[1503], m_cL[1503], m_dL[1503];
	double m_aR[1503], m_bR[1503], m_cR[1503], m_dR[1503];
	int ax, bx, cx, dx;
	//PeaksOnly
	double m_lastSampleL, m_lastSampleR;
	//SlewOnly
	double iirSampleAL, iirSampleBL, iirSampleCL, iirSampleDL, iirSampleEL, iirSampleFL, iirSampleGL;
	double iirSampleHL, iirSampleIL, iirSampleJL, iirSampleKL, iirSampleLL, iirSampleML, iirSampleNL, iirSampleOL, iirSamplePL;
	double iirSampleQL, iirSampleRL, iirSampleSL;
	double iirSampleTL, iirSampleUL, iirSampleVL;
	double iirSampleWL, iirSampleXL, iirSampleYL, iirSampleZL;

	double iirSampleAR, iirSampleBR, iirSampleCR, iirSampleDR, iirSampleER, iirSampleFR, iirSampleGR;
	double iirSampleHR, iirSampleIR, iirSampleJR, iirSampleKR, iirSampleLR, iirSampleMR, iirSampleNR, iirSampleOR, iirSamplePR;
	double iirSampleQR, iirSampleRR, iirSampleSR;
	double iirSampleTR, iirSampleUR, iirSampleVR;
	double iirSampleWR, iirSampleXR, iirSampleYR, iirSampleZR; // o/`
	//SubsOnly

  // Parameters

  float A;
  float B;
  float C;
  float D;
  float E;
  float F;
  float G;
  float H;
  float I;
  float J;
  float K;
  float L;
  float M;
  float N;
  float O;
  float P;
  float Q;
  float R;
  float S;
  float T;
  float U;
  float V;
  float W;
  float X;
  float Y;
  float Z;
  float Z0;
  float Z1;
  float Z2;
  float Z3;
  float Z4;
};

#endif
