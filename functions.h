#ifndef __functions_h__
#define __functions_h__

// #include "histoFunc.h"

#include <iostream>

#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSpectrum.h"

const unsigned int nS = 1024;
const float GS_s = 5.12;



//void CountDarkCounts(float* vals);

float DummyFunc();

//std::vector<std::pair<unsigned int,unsigned int> > CheckDarkCounts(float* vals);
//void AnalyzeDarkCounts(float* vals, unsigned int& pedSMin, unsigned int& pedSMax, unsigned int& ampSMin, unsigned int& ampNS);

std::pair<std::vector<float>*, std::vector<float>*> GetInvSubWaveform(std::vector<float> *times, std::vector<float> *amps, float& delay);

float GetPedestal (std::vector<float> *times, std::vector<float> *amps, const unsigned int& s1, const unsigned int& s2);
std::pair<float,float> GetAmplitude(std::vector<float> *times, std::vector<float> *amps, const unsigned int& s1, const unsigned int& s2, float& polarity, const unsigned int& nAve);

float GetIntegral (std::vector<float> *times, std::vector<float> *amps, const unsigned int& s1, const unsigned int& s2, float& polarity, const unsigned int& minS, const unsigned int& maxS);

std::pair<float, TF1*> GetTime (std::vector<float> *times, std::vector<float> *amps, const unsigned int& s1, const unsigned int& s2, float& polarity, const float& th, 
								    const unsigned int& minS, const unsigned int& maxS,
								    const unsigned int& nSBef, const unsigned int& nSAft,
                                                                    const bool& isFraction, const bool& deleteFit=true, 
								    const unsigned int& nAveAmp=1);

TGraph* PointsAroundTh(std::vector<float>* times, std::vector<float>* amps, float& polarity, const unsigned int& sBef, const unsigned int& nSBef, const unsigned int& nSAft);

double riseSingleDecay(double* x, double* par);
double riseSingleDecaySub(double* x, double* par);
double riseSingleDecayNoIRF(double* x, double* par);
double riseSingleDecayNoIRFSub(double* x, double* par);


#endif
