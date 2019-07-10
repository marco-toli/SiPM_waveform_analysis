// g++ -Wall -o pulses_analyzer.exe pulses_analyzer.cpp functions.cpp `root-config --cflags --glibs` -lSpectrum

#include "VectorSmallestInterval.hh"
#include "VectorSmallestInterval.cc"
// #include "fitUtils.h"
// #include "setTDRStyle.h"
// #include "ConfigParser.h"

#include "functions.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>    // std::sort
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"

#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"

#include "TSpline.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TSpectrum.h"
#include <TROOT.h>


using namespace std;


int main(int argc, char** argv)
{

  TApplication* theApp = new TApplication("App", &argc, argv);  TApplication a(“a”, 0, 0)

  
  std::string file_name  = "./root_files/";  
  file_name += std::string(argv[1]);
  std::cout << "input file name: " << file_name << std::endl; 
  
  std::string output_file_name  = "./output/";  
  output_file_name += std::string(argv[2]);
  std::cout << "output file name: " << output_file_name << std::endl; 

  float bias = 0.;
//   if (argc>3) bias = std::atof(argv[3]);
  
  TLegend * leg;
  
  //defining input tree
  int index;
  int NCH=2;
  
  std::vector<float> *t_time[NCH];
  std::vector<float> *t_amp[NCH];
  
  for(int iCh = 0; iCh<NCH;iCh++)
  {
      t_time[iCh]   = new std::vector<float>();
      t_amp[iCh]    = new std::vector<float>();
  }  
  
  TFile * inputFile = new TFile(file_name.c_str() , "READ");
  TTree * waveTree = (TTree*) inputFile->Get("waveTree");

  waveTree -> SetBranchAddress("index",   &index);
  waveTree -> SetBranchAddress("t_time_l",&t_time[0]);
  waveTree -> SetBranchAddress("t_amp_l", &t_amp[0]);
  waveTree -> SetBranchAddress("t_time",  &t_time[1]);
  waveTree -> SetBranchAddress("t_amp",   &t_amp[1]);
  
//   waveTree -> SetBranchAddress("t_time_l",     &t_time_l);
  
  std::cout << "read file, initializing..."<<std::endl;
  
//   std::map<int,float> ped;
//   std::map<int,float> amp;
//   std::map<int,float> time_linFit;
//   std::map<int,float> duration_linFit;
    
  unsigned int pedSMin[NCH];  
  unsigned int pedSMax[NCH];
  pedSMin[0] = 5;
  pedSMax[0] = 40;
  pedSMin[1] = 270;
  pedSMax[1] = 305;
//   pedSMin[1] = 5;
//   pedSMax[1] = 40;
  
  unsigned int pedSMin_sub[NCH];
  unsigned int pedSMax_sub[NCH];
  pedSMin_sub[0] = 5;
  pedSMax_sub[0] = 40;
  pedSMin_sub[1] = 245;
  pedSMax_sub[1] = 285;
  
  
  unsigned int nSBef[NCH];
  unsigned int nSAft[NCH];
  nSBef[0] = 2;
  nSAft[0] = 2;  
  nSBef[1] = 2;
  nSAft[1] = 3;
  
  
  float delay = 500;    // in ps for subtracted waveform
    
  unsigned int nSBef_sub[NCH];
  unsigned int nSAft_sub[NCH];
  nSBef_sub[0] = 2;
  nSAft_sub[0] = 2;  
  nSBef_sub[1] = 2;
  nSAft_sub[1] = 3;
  
  
  //samples to integrate signal
    
  unsigned int minS[NCH];
  unsigned int maxS[NCH];
  //-1 for full gate integration after pedestal
//   minS[0] = -1;
//   maxS[0] = -1;
//   minS[1] = -1;
//   maxS[1] = -1
  
  //for CTR
  minS[0] = 40;
  maxS[0] = 120;  
  minS[1] = 300;
  maxS[1] = 380;
  
  /*
  //for XTALK
  minS[0] = 100;
  maxS[0] = 220;  
  minS[1] = 100;
  maxS[1] = 200;
  */
  
  std::cout << "ok til here" << std::endl;
  
  waveTree->GetEntry(10);  
  
  float lowestT = t_time[0]->at(0);
  
  std::cout << "lowest Time =  " << lowestT << std::endl;
  
  for (int iCh = 0; iCh < NCH; iCh++)
  {
    std::cout << "minS[" << iCh << "] = " << minS[iCh] << " :: maxS = " << maxS[iCh] << std::endl;
  }
  
  
  //defining histos and graphs
  int NPULSES_EX = 4;
  TGraphErrors* gPulse[NPULSES_EX][NCH];
  TGraphErrors* gPulseSub[NPULSES_EX][NCH];
  
  
//   TF1 * fitPulseLED[NPULSES_EX][NCH];
  TF1 * fitPulseCFD[NPULSES_EX][NCH];
  TF1 * fitPulseCFD_sub[NPULSES_EX][NCH];
  
  float polarity[NCH];
  polarity[0] = 1.;
  polarity[1] = -1.;
  
  
  int NTH = 18;
  float stepTh = 0.05;
  float cfdTh[NTH];
//   float cfdTh[NCH];
  for (int iTh = 0; iTh<NTH; iTh++)
  {
    cfdTh[iTh] = 0.1+iTh*stepTh;
  }
  
  TH1F * hPed[NCH];
  TH1F * hPedSample[NCH];
  TH1F * hPedSampleSub[NCH];
  TH1F * hAmp[NCH];
  TH1F * hInt[NCH];
  
  TH1F * hTimeCFD[NCH][NTH];
  TH1F * hTimeLED[NCH][NTH];
  TH1F * hTimeCFD_sub[NCH][NTH];
  
  TH1F * hSlopeCFD[NCH][NTH];
//   TH1F * hSlopeLED[NCH][NTH];
  TH1F * hSlopeCFD_sub[NCH][NTH];
  
  TH1F * hTimeTemplateFit = new TH1F ("hTimeTemplateFit", "hTimeTemplateFit", 20000, 0, 50);
  TH1F * hSlopeTemplateFit = new TH1F ("hSlopeTemplateFit", "hSlopeTemplateFit", 20000, -10, 20);
  
  TH1F * hTimeTemplateFit_sub = new TH1F ("hTimeTemplateFit_sub", "hTimeTemplateFit_sub", 20000, 0, 50);
  TH1F * hSlopeTemplateFit_sub = new TH1F ("hSlopeTemplateFit_sub", "hSlopeTemplateFit_sub", 20000, -10, 20);
  
  TH1F * hChisquareLin     = new TH1F ("hChisquareLin", "hChisquareLin", 5000, 0, 1);
  TH1F * hChisquareLinSub  = new TH1F ("hChisquareLinSub", "hChisquareLinSub", 5000, 0, 1);
  
  TH1F * hChisquareTemp    = new TH1F ("hChisquareTemp", "hChisquareTemp", 5000, 0, 1);
  TH1F * hChisquareTempSub = new TH1F ("hChisquareTempSub", "hChisquareTempSub", 5000, 0, 1);
  
  
  TProfile * pAmpTimeRawLin  = new TProfile ("pAmpTimeRawLin", "pAmpTimeRawLin", 10000, -0.5, 5);
  TProfile * pAmpTimeSubLin  = new TProfile ("pAmpTimeSubLin", "pAmpTimeSubLin", 10000, -0.5, 5);
  TProfile * pAmpTimeRawTemp = new TProfile ("pAmpTimeRawTemp", "pAmpTimeRawTemp", 10000, -0.5, 5);
  TProfile * pAmpTimeSubTemp = new TProfile ("pAmpTimeSubTemp", "pAmpTimeSubTemp", 10000, -0.5, 5);
  
  
  for (int iCh = 0; iCh < NCH; iCh++)
  {
      hPed[iCh] = new TH1F (Form("hPed_%d", iCh), Form("hPed_%d", iCh), 2000, -0.5, 0.5);
      hPedSample[iCh] = new TH1F (Form("hPedSample_%d", iCh), Form("hPedSample_%d", iCh), 1000, -0.5, 0.5);
      hPedSampleSub[iCh] = new TH1F (Form("hPedSampleSub_%d", iCh), Form("hPedSampleSub_%d", iCh), 1000, -0.5, 0.5);
      hAmp[iCh] = new TH1F (Form("hAmp_%d", iCh), Form("hAmp_%d", iCh), 10000, -0.5, 5);
      hInt[iCh] = new TH1F (Form("hInt_%d", iCh), Form("hInt_%d", iCh), 10000, -5, 50);
      
//       hTimeCFD[iCh] = new TH1F (Form("hTimeCFD_%d", iCh), Form("hTimeCFD_%d", iCh), 20000, 0, 50);
//       hTimeLED[iCh] = new TH1F (Form("hTimeLED_%d", iCh), Form("hTimeLED_%d", iCh), 20000, 0, 50);
//       hTimeCFD_sub[iCh] = new TH1F (Form("hTimeCFD_sub_%d", iCh), Form("hTimeCFD_sub_%d", iCh), 20000, 0, 50);
      for (int iTh = 0; iTh<NTH; iTh++)
      {
        hTimeCFD[iCh][iTh]     = new TH1F (Form("hTimeCFD_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hTimeCFD_%d_th_%.2f", iCh, cfdTh[iTh]), 20000, 0, 50);
        hTimeLED[iCh][iTh]     = new TH1F (Form("hTimeLED_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hTimeLED_%d_th_%.2f", iCh, cfdTh[iTh]), 20000, 0, 50);
        hTimeCFD_sub[iCh][iTh] = new TH1F (Form("hTimeCFD_sub_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hTimeCFD_sub_%d_th_%.2f", iCh, cfdTh[iTh]), 20000, 0, 50);
        
        hSlopeCFD[iCh][iTh]     = new TH1F (Form("hSlopeCFD_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hSlopeCFD_%d_th_%.2f", iCh, cfdTh[iTh]), 4000, -10, 20);
//         hSlopeLED[iCh][iTh]     = new TH1F (Form("hSlopeLED_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hSlopeLED_%d_th_%.2f", iCh, cfdTh[iTh]), 4000, -10, 20);
        hSlopeCFD_sub[iCh][iTh] = new TH1F (Form("hSlopeCFD_sub_%d_th_%.2f", iCh, cfdTh[iTh]), Form("hSlopeCFD_sub_%d_th_%.2f", iCh, cfdTh[iTh]), 4000, -10, 20);
      
      }            
      
      for (int iPulse = 0; iPulse< NPULSES_EX; iPulse++)
      {        
        gPulse[iPulse][iCh] = new TGraphErrors();
        gPulseSub[iPulse][iCh] = new TGraphErrors();
      }
  }
  
//   TH1F * hCTR = new TH1F ("hCTR", "hCTR", 20000, 0, 50);
//   TH1F * hCTR_LED = new TH1F ("hCTR_LED", "hCTR_LED", 10000, 0, 50);
      
  

//   int NSAMPLES = 80000;
  
//   float time_led[NCH][NTH];
  float time_cfd[NCH][NTH];
  float time_cfd_sub[NCH][NTH];
  

  
//   ledTh[0] = 1.5;
//   ledTh[1] = 0.05;
  
  
  
  int nEntries = waveTree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  
//   nEntries = 10000;
  
  
  
  
  
  //**************************************************************************
  //                cycle 1 --> find single photoelectron value
  //**************************************************************************
  
  for (int iEntry = 0; iEntry < nEntries; iEntry++)
  {
      waveTree->GetEntry(iEntry);            
      
//       std::cout << "vector size = " << t_time[0]->size() << std::endl;
           
      if (t_time[0]->size()<10) continue;
      
//       std::cout << "iEntry: " << iEntry << std::endl;

      for (int iCh = 0; iCh <NCH; iCh++)
      {
          
        //filling pedestal for single samples in baseline gate
        for(unsigned int sIt = pedSMin[iCh]; sIt < pedSMax[iCh]; ++sIt)
        {
          hPedSample[iCh]->Fill(t_amp[iCh]->at(sIt));                    
        }
            
        float ped = GetPedestal(t_time[iCh], t_amp[iCh], pedSMin[iCh], pedSMax[iCh]);
        hPed[iCh]->Fill(ped);
        
        std::pair<float, float> pulse_amp = GetAmplitude(t_time[iCh], t_amp[iCh], pedSMin[iCh], pedSMax[iCh], polarity[iCh], 4);        
        hAmp[iCh]->Fill(pulse_amp.first);
        
        float integral = GetIntegral(t_time[iCh], t_amp[iCh], pedSMin[iCh], pedSMax[iCh], polarity[iCh], minS[iCh], maxS[iCh]);
        hInt[iCh]->Fill(integral);
  
      }
  }
  
  
  //drawing first plots
  float ped_noise[NCH];
  
  TCanvas * cPedSample[NCH];
  for (int iCh = 0; iCh<NCH; iCh++)
  {
      cPedSample[iCh] = new TCanvas (Form("cPedSample_%d", iCh), Form("cPedSample_%d", iCh), 500, 500);
      cPedSample[iCh]->cd();
      hPedSample[iCh]->Draw();
      
      hPedSample[iCh]->GetXaxis()->SetTitle("Time [s]");
      hPedSample[iCh]->GetYaxis()->SetTitle("Counts");
      
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -10, 10);
      hPedSample[iCh]->Fit(fitGaus, "QR");
      ped_noise[iCh] = fitGaus->GetParameter(2)*1e3;
  
      std::cout << "ped_noise [" << iCh << "] = " << ped_noise[iCh] << " mV " << std::endl;
  }
  /*
  TCanvas * cPed[NCH];
  for (int iCh = 0; iCh<NCH; iCh++)
  {
      cPed[iCh] = new TCanvas (Form("cPed_%d", iCh), Form("cPed_%d", iCh), 500, 500);
      cPed[iCh]->cd();
      hPed[iCh]->Draw();      
      hPed[iCh]->GetXaxis()->SetTitle("Time [s]");
      hPed[iCh]->GetYaxis()->SetTitle("Counts");
  }
  */
  TCanvas * cAmp[NCH];
  for (int iCh = 0; iCh<NCH; iCh++)
  {
      cAmp[iCh] = new TCanvas (Form("cAmp_%d", iCh), Form("cAmp_%d", iCh), 500, 500);
      cAmp[iCh]->cd();
      hAmp[iCh]->Draw();
      hAmp[iCh]->GetXaxis()->SetTitle("Max amplitude [mV]");
      hAmp[iCh]->GetYaxis()->SetTitle("Counts");
  }
  
    
  TCanvas * cInt[NCH];
  for (int iCh = 1; iCh<NCH; iCh++)
  {
      cInt[iCh] = new TCanvas (Form("cInt_%d", iCh), Form("cInt_%d", iCh), 500, 500);
      cInt[iCh]->cd();
      hInt[iCh]->Draw();
      hInt[iCh]->GetXaxis()->SetTitle("Signal integral [s*V]");
      hInt[iCh]->GetYaxis()->SetTitle("Counts");
  }
  
  
  
  
  
  
  
  //**************************************************************************
  //            cycle 2 --> apply event selection around NPHE
  //**************************************************************************
  

  
  //defining SiPM template waveform fit
  TF1 * fitRiseDecay = new TF1 ("fitRiseDecay", riseSingleDecay, 23.0e-9, 27.5e-9, 6);
  
    //template function parameters
  float myPar[6];
  if (bias == 67)
  {
      myPar[0] = 2.83829e-10;
      myPar[1] = 2.77329e-10;
      myPar[2] = 6.75821e-10;
      myPar[5] = 2.22547e-10;      
  }
  else if (bias == 68)
  {
      myPar[0] = 1.54741e-09;
      myPar[1] = 3.30974e-10;
      myPar[2] = 4.48079e-10;
      myPar[5] = 2.6e-10;      
  }
  else if (bias == 69)
  {
      myPar[0] = 3.24955e-09;
      myPar[1] = 2.85118e-10;
      myPar[2] = 5.14006e-10;
      myPar[5] = 1.93242e-10;
  }
  else if (bias == 70)
  {
      myPar[0] = 1.20466e-08;
      myPar[1] = 5.23e-9;
      myPar[2] = 2.34564e-10;
      myPar[5] = 9.3213e-11;      
  }
  else if (bias == 71)
  {
      myPar[0] = 1.08392e-08;
      myPar[1] = 4.40786e-09;
      myPar[2] = 4.11933e-10;      
      myPar[5] = 1.82012e-10;      
  }
  else if (bias == 72)
  {
      myPar[0] = 1.08392e-08;
      myPar[1] = 4.40786e-10;
      myPar[2] = 4.11933e-10;      
      myPar[5] = 1.82012e-10;   
  }  
  else
  {
      myPar[0] = 3.24955e-09;
      myPar[1] = 2.85118e-10;
      myPar[2] = 5.14006e-10;
      myPar[5] = 1.93242e-10;
  }
  
  fitRiseDecay->SetParameter(0, myPar[0]);       //normalization
  fitRiseDecay->SetParLimits(0, 0, 30.23501e-06);  
  fitRiseDecay->FixParameter(3, 0);
  fitRiseDecay->SetParameter(4, 2.51113e-08);
  
  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72) fitRiseDecay->SetParameter(1, myPar[1]);   //rise time
//   else                                                                                  
      fitRiseDecay->SetParLimits(1, 0.11686e-10, 100.31686e-09);   //rise time

  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72) fitRiseDecay->SetParameter(2, myPar[2]);  //decay time
//   else                                                                                  
      fitRiseDecay->SetParLimits(2, 0.5e-10, 51e-10);   //decay time    
  
  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72) fitRiseDecay->SetParameter(5, myPar[5]);        
//   else                                                                                  
      fitRiseDecay->SetParLimits(5, 1e-11, 20e-11);  
  

  
  
  TF1 * fitRiseDecaySub = new TF1 ("fitRiseDecaySub", riseSingleDecaySub, 24.0e-9, 27.0e-9, 6);
//   TF1 * fitRiseDecaySub = new TF1 ("fitRiseDecaySub", riseSingleDecay, 25.0e-9, 26.0e-9, 6);
  
  fitRiseDecaySub->SetParameter(0, myPar[0]);       //normalization
  fitRiseDecaySub->SetParLimits(0, 1e-08, 30e-06);  
  fitRiseDecaySub->FixParameter(3, delay*1e-12);
  fitRiseDecaySub->SetParameter(4, 2.545e-08);
  fitRiseDecaySub->SetParLimits(4, 2.5e-08, 2.6e-08);
  
    
  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72) fitRiseDecaySub->FixParameter(1, myPar[1]);   //rise time
//   else                                                                                  
      
//       fitRiseDecaySub->SetParLimits(1, 0.11686e-10, 10.31686e-09);   //rise time
  
  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72) fitRiseDecaySub->FixParameter(2, myPar[2]);  //decay time
//   else                                                                                  
//       fitRiseDecaySub->SetParLimits(2, 0.5e-11, 5e-9);   //decay time
      
  if (bias == 67 || bias == 68 || bias == 69 || bias == 70 || bias == 71 || bias == 72)  fitRiseDecaySub->SetParameter(5, myPar[5]);        
//    else                                                                                 
       fitRiseDecaySub->SetParLimits(5, 5e-11, 12e-11);  
  
//   TF1 * fitRiseDecaySub = new TF1 ("fitRiseDecaySub", "gaus", 23.0e-9, 27.50e-9);
//   fitRiseDecaySub->SetParameters(1, 25e-9, 0.5e-9);
  
  TGraphErrors * gAverageWaveform    = new TGraphErrors();
  TGraphErrors * gAverageWaveformSub = new TGraphErrors();
  std::vector<float> *amp_avg        = new std::vector<float>();  
  std::vector<float> *amp_avg_sub    = new std::vector<float>();  
  
  
  std::vector<double> *sigma_eff_raw_lin[NCH][NTH];  
  std::vector<double> *sigma_eff_sub_lin[NCH][NTH];
  std::vector<double> *sigma_eff_raw_tmp[NCH];
  std::vector<double> *sigma_eff_sub_tmp[NCH];
  
  for (int iCh = 0; iCh < NCH; iCh++)
  {
      
      sigma_eff_raw_tmp[iCh] = new std::vector<double> ();
      sigma_eff_sub_tmp[iCh] = new std::vector<double> ();
      
      for (int iTh = 0; iTh<NTH; iTh++)
      {          
          sigma_eff_raw_lin[iCh][iTh] = new std::vector<double> ();
          sigma_eff_sub_lin[iCh][iTh] = new std::vector<double> ();
                    
      }
  }
  
  
  int selPulse = 0;
  int nEfficiency = 0;
//   int NPHE_max = 1.;
  
  //Gundi subtraction
  std::vector<float> *subTimes;// = new std::vector<float>();
  std::vector<float> *subAmps;//  = new std::vector<float>();
  
      
  float ChiSquareCut = 1;
  float ChiSquareCutTemp = 1;
      
  int maxEvents = 9999;
  if (nEntries>maxEvents) nEntries = maxEvents;
  
  std::pair<float, float> sipm_amp;
  std::pair<float, float> sipm_amp_sub;
  
  for (int iEntry = 0; iEntry < nEntries; iEntry++)
  {
      waveTree->GetEntry(iEntry);            
      if (iEntry%10 == 0) std::cout << "processing entry: " << iEntry << "\r" << std::flush;
//       std::cout << "vector size = " << t_time[0]->size() << std::endl;
      if (t_time[1]->size()<10) continue;
           
//       std::cout << "iEntry: " << iEntry << std::endl;
      float ped_sipm        = GetPedestal(t_time[1], t_amp[1], pedSMin[1], pedSMax[1]);
//       float integral_sipm   = GetIntegral(t_time[1], t_amp[1], pedSMin[1], pedSMax[1], polarity[1], minS[1], maxS[1]);
      sipm_amp = GetAmplitude(t_time[1], t_amp[1], pedSMin[1], pedSMax[1], polarity[1], 6);        
      
      float time_zero = GetTime(t_time[0], t_amp[0], pedSMin[0], pedSMax[0], polarity[0], 0.55*hAmp[0]->GetMean(), minS[0], maxS[0], nSBef[0], nSAft[0], false, false, 4).first*1e9;
      float meanSiPM_amp = hAmp[1]->GetMean();
      
      if (sipm_amp.first<meanSiPM_amp*0.9 || sipm_amp.first>meanSiPM_amp*1.1)  continue;
//       if ( integral_sipm<(NPHE_min*single - single_sigma*sqrt(NPHE_min)) || integral_sipm>(NPHE_min*single + single_sigma*sqrt(NPHE_min))  ) continue;                                                         
    
      

      //scanning channels
      for (int iCh = 0; iCh <NCH; iCh++)
      {
                 
        std::pair<std::vector<float>*, std::vector<float>*> subWaveform = GetInvSubWaveform(t_time[iCh], t_amp[iCh], delay);
        subTimes = subWaveform.first;
        subAmps  = subWaveform.second;
        
        float ped_sipm_sub        = GetPedestal(subTimes, subAmps, pedSMin[1], pedSMax[1]);
        sipm_amp_sub = GetAmplitude(subTimes, subAmps, pedSMin[1], pedSMax[1], polarity[1], 4);        
        //filling pedestal for single samples in baseline gate
        for(unsigned int sIt = pedSMin[iCh]; sIt < pedSMax[iCh]; ++sIt)
        {
          hPedSampleSub[iCh]->Fill(subAmps->at(sIt));                    
        }
      
        
        std::pair<float, TF1*> time_info_cfd[NTH];
        std::pair<float, TF1*> time_info_cfd_sub[NTH];
        
        for (int iTh = 0; iTh <NTH; iTh++)
        {
            //CFD on standard pulse
            time_info_cfd[iTh] = GetTime(t_time[iCh], t_amp[iCh], pedSMin[iCh], pedSMax[iCh], polarity[iCh], cfdTh[iTh], minS[iCh], maxS[iCh], nSBef[iCh], nSAft[iCh], true, false, 6);
//         std::pair<float, TF1*> time_info_led = GetTime(t_time[iCh], t_amp[iCh], pedSMin[iCh], pedSMax[iCh], polarity[iCh], ledTh[iCh], nSBef[iCh], nSAft[iCh], false, false);        
            time_cfd[iCh][iTh] = time_info_cfd[iTh].first * 1e9;
//         time_led[iCh] = time_info_led.first * 1e9;
            
            //CFD on substracted waveform
            time_info_cfd_sub[iTh] = GetTime(subTimes, subAmps, pedSMin_sub[iCh], pedSMax_sub[iCh], polarity[iCh], cfdTh[iTh], minS[iCh], maxS[iCh], nSBef_sub[iCh], nSAft_sub[iCh], true, false, 4);        
            time_cfd_sub[iCh][iTh] = time_info_cfd_sub[iTh].first * 1e9;
            
            
//         std::cout << "time_cfd[" << iCh << "] = " << time_cfd[iCh] << " :: time_led[" << iCh << "] = " << time_led[iCh] << " :: integral_sipm = " << integral_sipm << std::endl;
//         std::cout << "time_cfd[" << iCh << "] = " << time_cfd[iCh] << " :: integral_sipm = " << integral_sipm << std::endl;

            if (time_cfd[iCh][iTh]>0 && time_cfd_sub[1][iTh]>0)
            {
            
                        
                //LED
//             TF1 * fitPulse_temp_led = time_info_led.second;
//             hTimeLED[iCh]->Fill(time_led[iCh]);
//             hSlopeLED[iCh] ->Fill(fitPulse_temp_led->GetParameter(1)*1e-9);

                //CFD
                TF1 * fitPulse_temp_cfd = time_info_cfd[iTh].second;       
            
            
//                 hChisquareLin[iTh]->Fill(fitPulse_temp_cfd[iTh]->GetChisquare());// / (fitPulse_temp_cfd->GetNDF()-1));
                if (fitPulse_temp_cfd->GetChisquare() <ChiSquareCut) // /(fitPulse_temp_cfd->GetNDF() ) -1)<ChiSquareCut) 
                {
                    if (iCh!=0) 
                    {
                        hTimeCFD[iCh][iTh] -> Fill(time_cfd[iCh][iTh]-time_zero);
                        sigma_eff_raw_lin[iCh][iTh]->push_back(time_cfd[iCh][iTh]-time_zero);
                        
                    }
                    else
                    {
                        hTimeCFD[iCh][iTh] -> Fill(time_cfd[iCh][iTh]);
                        sigma_eff_raw_lin[iCh][iTh]->push_back(time_cfd[iCh][iTh]);
                    }
                    
                    hSlopeCFD[iCh][iTh] ->Fill(fitPulse_temp_cfd->GetParameter(1)*1e-9);            
                    if (iCh == 1 && iTh == NTH/2) pAmpTimeRawLin->Fill(sipm_amp.first, time_cfd[iCh][iTh]);                
                }
                        
                //CFD sub
                TF1 * fitPulse_temp_cfd_sub = time_info_cfd_sub[iTh].second;        
            
            
//                 hChisquareLinSub[iTh]->Fill(fitPulse_temp_cfd_sub->GetChisquare());// / (fitPulse_temp_cfd_sub->GetNDF()-1));
                if (fitPulse_temp_cfd_sub->GetChisquare()<ChiSquareCut) // / (fitPulse_temp_cfd_sub->GetNDF())-1)<ChiSquareCut) 
                {
                    if (iCh!=0) 
                    {
                        hTimeCFD_sub[iCh][iTh] -> Fill(time_cfd_sub[iCh][iTh]-time_zero);
                        sigma_eff_sub_lin[iCh][iTh]->push_back(time_cfd_sub[iCh][iTh]-time_zero);
                    }
                    else
                    {
                        hTimeCFD_sub[iCh][iTh] -> Fill(time_cfd_sub[iCh][iTh]);
                        sigma_eff_sub_lin[iCh][iTh]->push_back(time_cfd_sub[iCh][iTh]);
                    }
                    
                    
                    hSlopeCFD_sub[iCh][iTh] -> Fill(fitPulse_temp_cfd_sub->GetParameter(1)*1e-9);            
                    if (iCh == 1 && iTh == NTH/2) pAmpTimeSubLin->Fill(sipm_amp_sub.first, time_cfd_sub[iCh][iTh]);
                }      
            }
        }
            
            
        //fill test waveforms
        if (time_cfd[iCh][NTH/2]>0 && time_cfd_sub[1][NTH/2]>0)
        {
            if (selPulse< NPULSES_EX) 
            {
                
                for (unsigned int iSample = 0; iSample < t_time[iCh]->size(); iSample++)
                {                    
                    gPulse[selPulse][iCh]->SetPoint(iSample,  t_time[iCh]->at(iSample), t_amp[iCh]->at(iSample)*polarity[iCh]);                              
//                     gPulse[selPulse][iCh]->SetPointError(iSample,  0, ped_noise[iCh]/1e3);                              
//                     std::cout << "selectedPulse[" << iCh << "] = " << selPulse << " :: time[" << iSample << "] = " << t_time[iCh]->at(iSample) << " :: amp = " << t_amp[iCh]->at(iSample) << std::endl;
                }
                fitPulseCFD[selPulse][iCh] = time_info_cfd[NTH/2].second;                
                
                
                for (int i = 0; i<1; i++) gPulse[selPulse][iCh]->Fit(fitRiseDecay, "QR");
                                
                
                for (unsigned int iSample = 0; iSample < subTimes->size(); iSample++)
                {
                    gPulseSub[selPulse][iCh]->SetPoint(iSample, subTimes->at(iSample), subAmps->at(iSample)*polarity[iCh]);  
//                     gPulseSub[selPulse][iCh]->SetPointError(iSample,  0, ped_noise[iCh]/1e3);                              
//                    std::cout << "selectedPulse[" << iCh << "] = " << selPulse << " :: time[" << iSample << "] = " << t_time[iCh]->at(iSample) << " :: amp = " << t_amp[iCh]->at(iSample) << std::endl;
                }
                
                fitPulseCFD_sub[selPulse][iCh] = time_info_cfd_sub[NTH/2].second;                
                                
                fitRiseDecaySub->SetLineColor(kCyan+1);
                for (int i = 0; i<3; i++) gPulseSub[selPulse][iCh]->Fit(fitRiseDecaySub, "QR");
                
                
                if (iCh==1) selPulse++;                            
                //             fitPulseLED[iEntry][iCh] = time_info_led.second;
            }
        }
        
        
        //fill average waveforms
        if (iCh == 1)
        {
            //integrating average SiPM waveform
//             std::cout << "filling average waveforms" << std::endl;
            for(unsigned int sIt = 0; sIt < t_amp[1]->size(); ++sIt)
            {
                
                if (nEfficiency == 0) amp_avg    ->push_back(t_amp[1]->at(sIt) - ped_sipm);                    
                else                  amp_avg    ->at(sIt)+=(t_amp[1]->at(sIt) - ped_sipm);
            }
//             std::cout << "filling average waveform subtracted" << std::endl;
            for(unsigned int sIt = 0; sIt < subTimes->size(); ++sIt)
            {
                if (nEfficiency == 0) amp_avg_sub    ->push_back(subAmps->at(sIt) - ped_sipm_sub);                    
                else                  amp_avg_sub    ->at(sIt)+=(subAmps->at(sIt) - ped_sipm_sub);
            }
            
            nEfficiency++;
        }
        
        
        
        //*********************************************
        //try template fit on SiPM waveform
        //*********************************************
        
        TGraphErrors * gPulseSingle = new TGraphErrors();
        for (unsigned int iSample = 0; iSample < t_time[1]->size(); iSample++)
        {             
              gPulseSingle->SetPoint(iSample,  t_time[1]->at(iSample), t_amp[1]->at(iSample)*polarity[1]);                              
//               gPulseSingle->SetPointError(iSample, 0, ped_noise[1]/1.0e3);                              
//                     std::cout << "selectedPulse[" << iCh << "] = " << selPulse << " :: time[" << iSample << "] = " << t_time[iCh]->at(iSample) << " :: amp = " << t_amp[iCh]->at(iSample) << std::endl;
        }
        for (int i = 0; i<1; i++) gPulseSingle->Fit(fitRiseDecay, "QR");
        hChisquareTemp->Fill(fitRiseDecay->GetChisquare());// / (fitRiseDecay->GetNDF()-1));
//         std::cout << "time from template fit = " << fitRiseDecay->GetParameter(4) << std::endl;
        
        
        double maxT  = fitRiseDecay->GetMaximumX();
        double zeroT = fitRiseDecay->GetParameter(4);
//         double cfdT  = (maxT+zeroT)/2.;
        float cfdT  = fitRiseDecay->GetX(maxT*0.5, zeroT, maxT, 1e-12);
        double slope = (fitRiseDecay->Eval(cfdT+10.0e-12) - fitRiseDecay->Eval(cfdT-10.0e-12)) / 4.0e-12*1e-9;
//         std::cout << "maxT = " << maxT << " :: zeroT = " << zeroT << " :: cfdT = " << cfdT << " :: A_up = " << fitRiseDecay->Eval(cfdT+2.0e-12) << " :: A_low = " << fitRiseDecay->Eval(cfdT-2.0e-12) << " :: Delta_A = " <<fitRiseDecay->Eval(cfdT+2.0e-12)-fitRiseDecay->Eval(cfdT-2.0e-12) << " ::  slope = " << slope << std::endl;
        if (fitRiseDecay->GetChisquare()<ChiSquareCutTemp) // / (fitRiseDecay->GetNDF()-1)<ChiSquareCut) 
        {
            hTimeTemplateFit->Fill(fitRiseDecay->GetParameter(4)/1e-9);        
            sigma_eff_raw_tmp[iCh]->push_back(fitRiseDecay->GetParameter(4)/1e-9);
            hSlopeTemplateFit->Fill(slope);
            pAmpTimeRawTemp->Fill(fitRiseDecay->GetMaximum(), fitRiseDecay->GetParameter(4)/1e-9);
        }
        gPulseSingle->Delete();
        
        
        
        
        //*********************************************
        //try template fit on SiPM waveform subtracted
        //*********************************************
        
        TGraphErrors * gPulseSingleSub = new TGraphErrors();
        for (unsigned int iSample = 0; iSample < subTimes->size(); iSample++)
        {             
              gPulseSingleSub->SetPoint(iSample,  subTimes->at(iSample), subAmps->at(iSample)*polarity[1]);            
//               gPulseSingleSub->SetPointError(iSample, 0, ped_noise[1]/1.0e3);                              
//                     std::cout << "selectedPulse[" << iCh << "] = " << selPulse << " :: time[" << iSample << "] = " << t_time[iCh]->at(iSample) << " :: amp = " << t_amp[iCh]->at(iSample) << std::endl;
        }

        for (int i = 0; i<1; i++) gPulseSingleSub->Fit(fitRiseDecaySub, "QR");
        
        hChisquareTempSub->Fill(fitRiseDecaySub->GetChisquare()); // / (fitRiseDecaySub->GetNDF()-1));
//         std::cout << "time from template fit = " << fitRiseDecay->GetParameter(4) << std::endl;
        
//         hTimeTemplateFit_sub->Fill(fitRiseDecaySub->GetParameter(1)/1e-9);        
        
        maxT  = fitRiseDecaySub->GetMaximumX();
        zeroT = fitRiseDecaySub->GetParameter(4);
//         zeroT = fitRiseDecaySub->GetParameter(1);
//         double cfdT  = (maxT+zeroT)/2.;
        cfdT  = fitRiseDecaySub->GetX(maxT*0.5, zeroT, maxT, 1e-12);
        slope = (fitRiseDecaySub->Eval(cfdT+10.0e-12) - fitRiseDecaySub->Eval(cfdT-10.0e-12)) / 4.0e-12*1e-9;
//         std::cout << "maxT = " << maxT << " :: zeroT = " << zeroT << " :: cfdT = " << cfdT << " :: A_up = " << fitRiseDecay->Eval(cfdT+2.0e-12) << " :: A_low = " << fitRiseDecay->Eval(cfdT-2.0e-12) << " :: Delta_A = " <<fitRiseDecay->Eval(cfdT+2.0e-12)-fitRiseDecay->Eval(cfdT-2.0e-12) << " ::  slope = " << slope << std::endl;
        
        if (fitRiseDecaySub->GetChisquare()<ChiSquareCutTemp && slope > 0.1) // / (fitRiseDecaySub->GetNDF()-1)<ChiSquareCut) 
        {
            hSlopeTemplateFit_sub->Fill(slope);         
            hTimeTemplateFit_sub->Fill(fitRiseDecaySub->GetParameter(4)/1e-9);        
            sigma_eff_sub_tmp[iCh]->push_back(fitRiseDecaySub->GetParameter(4)/1e-9);
            pAmpTimeSubTemp->Fill(fitRiseDecaySub->GetMaximum(), fitRiseDecaySub->GetParameter(4)/1e-9);
                   
        }
        gPulseSingleSub->Delete(); 
//         std::cout << "ch[" << iCh << "] --> time_cfd = " << time_cfd[iCh] << " :: time_led = " << time_led[iCh] << std::endl;
      }      
      
      
      
//       hCTR->Fill(time_cfd[1]-time_cfd[0]);
//       hCTR->Fill(time_cfd[1][0]);
//       hCTR_LED->Fill(time_led[1]-time_led[0]);
      
  }
  std::cout << std::endl;
  
  
  //average SiPM waveform
  waveTree->GetEntry(10);            
  std::cout << "selected Waveforms: " << nEfficiency << std::endl;
  
  for(unsigned int sIt = 0; sIt < amp_avg->size(); ++sIt)
  {            
      gAverageWaveform   ->SetPoint(sIt, t_time[1]->at(sIt),  -amp_avg      ->at(sIt)/nEfficiency);      
  }
  
  for(unsigned int sIt = 0; sIt < amp_avg_sub->size(); ++sIt)
  {            
      gAverageWaveformSub->SetPoint(sIt, subTimes->at(sIt),  -amp_avg_sub  ->at(sIt)/nEfficiency);
  }
  
  TCanvas * cAverageWaveform = new TCanvas ("cAverageWaveform", "cAverageWaveform", 600, 500);
  cAverageWaveform->cd();
  gAverageWaveform->GetXaxis()->SetTitle("Time [s]");
  gAverageWaveform->GetYaxis()->SetTitle("Amplitude [V]");
  gAverageWaveform->GetXaxis()->SetTitleSize(0.05);
  gAverageWaveform->GetYaxis()->SetTitleSize(0.05);
  gAverageWaveform->GetXaxis()->SetTitleOffset(0.95);
  gAverageWaveform->GetYaxis()->SetTitleOffset(0.95);
  gAverageWaveform->GetXaxis()->SetRangeUser(20e-9, 40e-9);
  gAverageWaveform->GetYaxis()->SetRangeUser(-0.05, 2.5);
  gAverageWaveform->SetMarkerStyle(7);
  gAverageWaveform->SetLineWidth(2);
  gAverageWaveform->Draw("ALP");
  
//   TF1 * fitRiseDecayAve = new TF1 ("fitRiseDecayAve", riseSingleDecayNoIRF, -5.0e-9, 12.0e-9, 5);
  TF1 * fitRiseDecayAve = new TF1 ("fitRiseDecayAve", riseSingleDecay, 15.0e-9, 45.0e-9, 6);
  fitRiseDecayAve->SetParameter(0, fitRiseDecay->GetParameter(0));     //normalization
  fitRiseDecayAve->SetParameter(1, fitRiseDecay->GetParameter(1));	//rise time
  fitRiseDecayAve->SetParameter(2, fitRiseDecay->GetParameter(2));	//decay time  
  fitRiseDecayAve->FixParameter(3, 0);	//electronic delay
  fitRiseDecayAve->SetParameter(4, fitRiseDecay->GetParameter(4));	//delay
  
  fitRiseDecayAve->SetParameter(5, 2.6e-10);	//impulse response function
//    fitRiseDecay->SetParameter(5, 0);	//impulse response function
    
  for (int i = 0; i<3; i++) gAverageWaveform->Fit(fitRiseDecayAve, "R");

  
  
  
  TF1 * fitRiseDecayAveSub = new TF1 ("fitRiseDecayAveSub", riseSingleDecay, 24.0e-9, 26.5e-9, 6);
  fitRiseDecayAveSub->SetParameter(0, fitRiseDecaySub->GetParameter(0));     //normalization
  fitRiseDecayAveSub->SetParameter(1, fitRiseDecaySub->GetParameter(1));	//rise time
  fitRiseDecayAveSub->SetParameter(2, fitRiseDecaySub->GetParameter(2));	//decay time  
  fitRiseDecayAveSub->FixParameter(3, fitRiseDecaySub->GetParameter(3));	//electronic delay
  fitRiseDecayAveSub->SetParameter(4, fitRiseDecaySub->GetParameter(4));	//delay
  fitRiseDecayAveSub->SetParameter(5, fitRiseDecaySub->GetParameter(5));	//impulse response function
  
//
  gAverageWaveformSub->SetLineColor(kBlue+1);
  gAverageWaveformSub->SetMarkerColor(kBlue+1);
  gAverageWaveformSub->SetMarkerStyle(7);
  gAverageWaveformSub->SetLineWidth(2);
  gAverageWaveformSub->Draw("same LP");
  
//   fitRiseDecayAveSub->SetParameter(2, fitRiseDecaySub->GetParameter(2)/3);	//decay time    
  fitRiseDecayAveSub->SetLineColor(kCyan);
  for (int i = 0; i<3; i++) gAverageWaveformSub->Fit(fitRiseDecayAveSub, "R");
  
  leg = new TLegend(0.56,0.66,0.88,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  
  leg->AddEntry(gAverageWaveform, "Raw waveform", "lpe");   
  leg->AddEntry(gAverageWaveformSub, "Processed waveform", "lpe");   
  leg->Draw("same");
  
  
  
  //drawing stuff
  
  float ped_noise_sub[NCH];
  
  TCanvas * cPedSampleSub[NCH];
  for (int iCh = 1; iCh<NCH; iCh++)
  {
      cPedSampleSub[iCh] = new TCanvas (Form("cPedSampleSub_%d", iCh), Form("cPedSampleSub_%d", iCh), 500, 500);
      cPedSampleSub[iCh]->cd();
      hPedSampleSub[iCh]->Draw();
      
      hPedSampleSub[iCh]->GetXaxis()->SetTitle("Time [s]");
      hPedSampleSub[iCh]->GetYaxis()->SetTitle("Counts");
      
      TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -10, 10);
      hPedSampleSub[iCh]->Fit(fitGaus, "QR");
      ped_noise_sub[iCh] = fitGaus->GetParameter(2)*1e3;
  
      std::cout << "ped_noise_sub [" << iCh << "] = " << ped_noise_sub[iCh] << " mV " << std::endl;
  }
  
  
  
  TCanvas * cChisquares = new TCanvas ("cChisquares", "cChisquares", 600, 600);
  cChisquares->cd();
  hChisquareTemp->SetLineColor(kRed+1);
  hChisquareTemp->GetXaxis()->SetTitle("ChiSquare / NDF");
  hChisquareTemp->GetYaxis()->SetTitle("Counts");
  hChisquareTemp->GetXaxis()->SetRangeUser(0, 0.03);
  hChisquareTemp->Draw();
  
  hChisquareTempSub->SetLineColor(kCyan+1);
  hChisquareTempSub->Draw("same");
  
  hChisquareLin->SetLineColor(kOrange+1);
  hChisquareLin->Draw("same");
  
  hChisquareLinSub->SetLineColor(kBlue+1);
  hChisquareLinSub->Draw("same");
  
  TLine * linCutChiSquare = new TLine(ChiSquareCut, 0, ChiSquareCut, 1000);
  linCutChiSquare->SetLineColor(kBlack);
  linCutChiSquare->SetLineWidth(2);
  linCutChiSquare->SetLineStyle(7);
  linCutChiSquare->Draw("same");
  
  gPad->SetLogy();
  
  leg = new TLegend(0.45,0.66,0.88,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  
  leg->AddEntry(hChisquareLin, "Raw waveform: linear fit", "lpe");   
  leg->AddEntry(hChisquareLinSub, "Processed waveform: linear fit", "lpe");   
  leg->AddEntry(hChisquareTemp, "Raw waveform: template fit", "lpe");   
  leg->AddEntry(hChisquareTempSub, "Processed waveform: template fit", "lpe");   
  leg->Draw("same");
  
  
  
  TGraphErrors* gSlope_vs_Th[NCH];
  TGraphErrors* gSlopeSub_vs_Th[NCH];
  TGraphErrors* gTimeCFD_vs_Th[NCH];
  TGraphErrors* gTimeCFDSub_vs_Th[NCH];
  
  TGraphErrors* gTimeEffCFD_vs_Th[NCH];
  TGraphErrors* gTimeEffCFDSub_vs_Th[NCH];
  
  for (int iCh = 0; iCh < NCH; iCh++)
  {
      gSlope_vs_Th[iCh]      = new TGraphErrors();
      gSlopeSub_vs_Th[iCh]   = new TGraphErrors();
      gTimeCFD_vs_Th[iCh]    = new TGraphErrors();
      gTimeCFDSub_vs_Th[iCh] = new TGraphErrors();
      gTimeEffCFD_vs_Th[iCh]    = new TGraphErrors();
      gTimeEffCFDSub_vs_Th[iCh] = new TGraphErrors();
  }
  
  //calculating dV/dt
  std::cout << "drawing slope... " << std::endl;
  float slope_cfd[NCH][NTH];
//   float slope_led[NCH][NTH];
  float slope_cfd_sub[NCH][NTH];
  
  TCanvas * cSlope[NCH];
  for (int iCh = 0; iCh<NCH; iCh++)
  {
      cSlope[iCh] = new TCanvas (Form("cSlope_%d", iCh), Form("cSlope_%d", iCh), 500, 500);
      cSlope[iCh]->cd();
      hSlopeCFD[iCh][0]->Draw();
//       hSlopeLED[iCh]->SetLineColor(kGreen+2);
//       hSlopeLED[iCh]->Draw("same");
      for (int iTh = 0; iTh<NTH; iTh++)
      {        
        hSlopeCFD[iCh][iTh]->Draw("same");
        hSlopeCFD[iCh][iTh]->SetLineColor(iTh+1);        
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -5, 10);
        hSlopeCFD[iCh][iTh]->Fit(fitGaus, "QR");
        slope_cfd[iCh][iTh] = fitGaus->GetParameter(1);
        std::cout << "slope CFD [" << iCh << "] for th [" << cfdTh[iTh] << "] = " << slope_cfd[iCh][iTh] << " mV/ns" << std::endl;
        gSlope_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], slope_cfd[iCh][iTh]);
      }
  }
  
  TCanvas * cSlopeSub[NCH];
  for (int iCh = 1; iCh<NCH; iCh++)
  {
      cSlopeSub[iCh] = new TCanvas (Form("cSlopeSub_%d", iCh), Form("cSlopeSub_%d", iCh), 500, 500);
      cSlopeSub[iCh]->cd();
      hSlopeCFD_sub[iCh][0]->Draw();
      
      for (int iTh = 0; iTh<NTH; iTh++)
      {
        hSlopeCFD_sub[iCh][iTh]->Draw("same");
        hSlopeCFD_sub[iCh][iTh]->SetLineColor(iTh+1);
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -5, 10);
        hSlopeCFD_sub[iCh][iTh]->Fit(fitGaus, "QR");
        slope_cfd_sub[iCh][iTh] = fitGaus->GetParameter(1);
        std::cout << "slope CFD _sub [" << iCh << "] for th [" << cfdTh[iTh] << "] = " << slope_cfd_sub[iCh][iTh] << " mV/ns" << std::endl;
        gSlopeSub_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], slope_cfd_sub[iCh][iTh]);
      }
  }
  
  TCanvas * cSlopes_vsTh = new TCanvas ("cSlopes_vsTh", "cSlopes_vsTh", 500, 500);
  cSlopes_vsTh->cd();
  gSlope_vs_Th[1]->Draw("ALPE");
  gSlopeSub_vs_Th[1]->SetLineColor(kBlue+1);
  gSlopeSub_vs_Th[1]->Draw("same LPE");
  
  
  float sigma_cfd[NCH][NTH];
  float sigma_cfd_sub[NCH][NTH];
  
  TCanvas * cTime[NCH];
  for (int iCh = 0; iCh<NCH; iCh++)
  {
      cTime[iCh] = new TCanvas (Form("cTime_%d", iCh), Form("cTime_%d", iCh), 500, 500);
      cTime[iCh]->cd();
      
      hTimeCFD[iCh][0]->Draw();                
      hTimeCFD[iCh][0]->GetXaxis()->SetRangeUser(hTimeCFD[iCh][0]->GetMean()*0.5, hTimeCFD[iCh][0]->GetMean()*1.5);
//       hTimeLED[iCh]->SetLineColor(kGreen+2);
//       hTimeLED[iCh]->Draw("same");
      for (int iTh = 0; iTh< NTH; iTh++)
      {
        for (int iR = 0; iR < 5; iR++)
        {
            if (hTimeCFD[iCh][iTh]->GetBinContent(hTimeCFD[iCh][iTh]->GetMaximumBin() ) <30) hTimeCFD[iCh][iTh]->Rebin(2);
            else break;
        }
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", hTimeCFD[iCh][iTh]->GetMean()-hTimeCFD[iCh][iTh]->GetRMS(), hTimeCFD[iCh][iTh]->GetMean()+hTimeCFD[iCh][iTh]->GetRMS()*2);
//         fitGaus->SetParameters(10,  hTimeCFD[iCh][iTh]->GetMean(), hTimeCFD[iCh][iTh]->GetRMS()/2);
        hTimeCFD[iCh][iTh]->Fit(fitGaus, "QR");
        sigma_cfd[iCh][iTh] = fitGaus->GetParameter(2);
      
//       hTimeLED[iCh]->Fit(fitGaus, "QR");
//       float sigma_led = fitGaus->GetParameter(2);
  
//       std::cout << "CTR CFD (sigma) [" << iCh << "] = " << sigma_cfd*1e3 << " ps :: CTR LED (sigma) = " << sigma_led*1e3 << " ps :: CTR intrinsic noise (sigma) = " << ped_noise[iCh]/slope_cfd[iCh] << " ps "  << std::endl;
        if (slope_cfd[iCh][iTh]!=0) std::cout << "CTR CFD (sigma) [" << iCh << "]  for th[" << cfdTh[iTh] << "] = " << sigma_cfd[iCh][iTh]*1e3 << " ps :: CTR intrinsic noise (sigma) = " << ped_noise[iCh]/slope_cfd[iCh][iTh] << " ps "  << " :: CTR without noise = " << sqrt(pow(sigma_cfd[iCh][iTh]*1e3,2) - pow(ped_noise[iCh]/slope_cfd[iCh][iTh],2)) << std::endl;
        gTimeCFD_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], sigma_cfd[iCh][iTh]*1e3);
        gTimeCFD_vs_Th[iCh]->SetPointError(iTh, 0, fitGaus->GetParError(2)*1e3);
//         std::cout << "here ..." <<std::endl;
        
        double sigma_eff = sigma_cfd[iCh][iTh];
        if (iCh != 0) sigma_eff = FindSmallestInterval(hTimeCFD[iCh][iTh]->GetMean(), hTimeCFD[iCh][iTh]->GetMeanError(), 0, 2000, sigma_eff_raw_lin[iCh][iTh], 0.68, 0)/2.;
        
//         std::cout << "here ..." <<std::endl;
        gTimeEffCFD_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], sigma_eff*1e3);
      }      
  }
  
  
  TCanvas * cTimeSub[NCH];
  for (int iCh = 1; iCh<NCH; iCh++)
  {
      cTimeSub[iCh] = new TCanvas (Form("cTimeSub_%d", iCh), Form("cTimeSub_%d", iCh), 500, 500);
      cTimeSub[iCh]->cd();
      
      hTimeCFD_sub[iCh][0]->Draw();
      hTimeCFD_sub[iCh][0]->GetXaxis()->SetRangeUser(hTimeCFD_sub[iCh][0]->GetMean()*0.5, hTimeCFD_sub[iCh][0]->GetMean()*1.5);
//       hTimeLED[iCh]->SetLineColor(kGreen+2);
//       hTimeLED[iCh]->Draw("same");
      for (int iTh = 0; iTh< NTH; iTh++)
      {
        for (int iR = 0; iR < 5; iR++)
        {
            if (hTimeCFD_sub[iCh][iTh]->GetBinContent(hTimeCFD_sub[iCh][iTh]->GetMaximumBin() ) <30) hTimeCFD_sub[iCh][iTh]->Rebin(2);
            else break;
        }
        TF1 * fitGaus = new TF1 ("fitGaus", "gaus", hTimeCFD_sub[iCh][iTh]->GetMean()-hTimeCFD_sub[iCh][iTh]->GetRMS(), hTimeCFD_sub[iCh][iTh]->GetMean()+hTimeCFD_sub[iCh][iTh]->GetRMS()*2);
//         fitGaus->SetParameters(10,  hTimeCFD_sub[iCh][iTh]->GetMean(), hTimeCFD_sub[iCh][iTh]->GetRMS()/2);
        hTimeCFD_sub[iCh][iTh]->Fit(fitGaus, "QR");
        sigma_cfd_sub[iCh][iTh] = fitGaus->GetParameter(2);
      
//       hTimeLED[iCh]->Fit(fitGaus, "QR");
//       float sigma_led = fitGaus->GetParameter(2);
  
//       std::cout << "CTR CFD (sigma) [" << iCh << "] = " << sigma_cfd*1e3 << " ps :: CTR LED (sigma) = " << sigma_led*1e3 << " ps :: CTR intrinsic noise (sigma) = " << ped_noise[iCh]/slope_cfd[iCh] << " ps "  << std::endl;
        if (slope_cfd_sub[iCh][iTh]!=0) std::cout << "CTR CFD sub (sigma) [" << iCh << "] for th[" << cfdTh[iTh] << "] = " << sigma_cfd_sub[iCh][iTh]*1e3 << " ps :: CTR intrinsic noise (sigma) = " << ped_noise_sub[iCh]/slope_cfd_sub[iCh][iTh] << " ps "  << " :: CTR sub without noise = " << sqrt(pow(sigma_cfd_sub[iCh][iTh]*1e3,2) - pow(ped_noise_sub[iCh]/slope_cfd_sub[iCh][iTh],2)) << std::endl;
        gTimeCFDSub_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], sigma_cfd_sub[iCh][iTh]*1e3);
        gTimeCFDSub_vs_Th[iCh]->SetPointError(iTh, 0, fitGaus->GetParError(2)*1e3);

                
        double sigma_eff = FindSmallestInterval(hTimeCFD_sub[iCh][iTh]->GetMean(), hTimeCFD_sub[iCh][iTh]->GetMeanError(), 0, 2000, sigma_eff_sub_lin[iCh][iTh], 0.68, 0)/2.;
        gTimeEffCFDSub_vs_Th[iCh]->SetPoint(iTh, cfdTh[iTh], sigma_eff*1e3);
      }            
      
  }
  
  
  TCanvas * cTimeCFD_vsTh = new TCanvas ("cTimeCFD_vsTh", "cTimeCFD_vsTh", 500, 500);
  cTimeCFD_vsTh->cd();
  gTimeCFD_vs_Th[1]->Draw("ALPE");
  gTimeCFD_vs_Th[1]->GetXaxis()->SetTitle("CFD threshold");
  gTimeCFD_vs_Th[1]->GetYaxis()->SetTitle("#sigma_{t}");
  
  
  gTimeCFDSub_vs_Th[1]->SetLineColor(kBlue+1);
  gTimeCFDSub_vs_Th[1]->Draw("same LPE");
  
  
  gTimeEffCFD_vs_Th[1]->SetLineColor(kBlack);
  gTimeEffCFD_vs_Th[1]->SetLineStyle(7);
  gTimeEffCFD_vs_Th[1]->Draw("same LPE");
  
  gTimeEffCFDSub_vs_Th[1]->SetLineColor(kBlue+1);
  gTimeEffCFDSub_vs_Th[1]->SetLineStyle(7);
  gTimeEffCFDSub_vs_Th[1]->Draw("same LPE");
  
  gTimeCFD_vs_Th[0]->SetLineColor(kRed+1);
  gTimeCFD_vs_Th[0]->Draw("same LPE");
  
  
    

  
  TCanvas * cTimeTemplateFit = new TCanvas (Form("cTimeTemplateFit"), Form("cTimeTemplateFit"), 500, 500);
  cTimeTemplateFit->cd();
  hTimeTemplateFit->Draw();
  hTimeTemplateFit->GetXaxis()->SetRangeUser(hTimeTemplateFit->GetMean()*0.5, hTimeTemplateFit->GetMean()*1.5);
      
  TF1 * fitGausTemp = new TF1 ("fitGausTemp", "gaus", hTimeTemplateFit->GetMean()-hTimeTemplateFit->GetRMS()*2, hTimeTemplateFit->GetMean()+hTimeTemplateFit->GetRMS()*2);
  hTimeTemplateFit->Fit(fitGausTemp, "QR");
  float sigma_template = fitGausTemp->GetParameter(2);
  
  
  TCanvas * cSlopeTemplateFit = new TCanvas (Form("cSlopeTemplateFit"), Form("cSlopeTemplateFit"), 500, 500);
  cSlopeTemplateFit->cd();
  hSlopeTemplateFit->Draw();
      
  fitGausTemp = new TF1 ("fitGausTemp", "gaus", hSlopeTemplateFit->GetMean()-hSlopeTemplateFit->GetRMS()*2, hSlopeTemplateFit->GetMean()+hSlopeTemplateFit->GetRMS()*2);
  hSlopeTemplateFit->Fit(fitGausTemp, "QR");
  float slope_template = fitGausTemp->GetParameter(1);
  

  std::cout << "CTR template fit (sigma) = " << sigma_template*1e3 <<" ps :: CTR intrinsic noise (sigma) temp = " << ped_noise[1]/slope_template << " ps "  << " :: CTR temp without noise = " << sqrt(pow(sigma_template*1e3,2) - pow(ped_noise[1]/slope_template,2)) <<  std::endl;
      

  TCanvas * cTimeTemplateFitSub = new TCanvas (Form("cTimeTemplateFitSub"), Form("cTimeTemplateFitSub"), 500, 500);
  cTimeTemplateFitSub->cd();
  hTimeTemplateFit_sub->Draw();
  hTimeTemplateFit_sub->GetXaxis()->SetRangeUser(hTimeTemplateFit_sub->GetMean()*0.5, hTimeTemplateFit_sub->GetMean()*1.5);
      
  fitGausTemp = new TF1 ("fitGausTemp", "gaus", hTimeTemplateFit_sub->GetMean()-hTimeTemplateFit_sub->GetRMS()*2, hTimeTemplateFit_sub->GetMean()+hTimeTemplateFit_sub->GetRMS()*2);
  hTimeTemplateFit_sub->Fit(fitGausTemp, "QR");
  float sigma_template_sub = fitGausTemp->GetParameter(2);
  
  
  TCanvas * cSlopeTemplateFitSub = new TCanvas (Form("cSlopeTemplateFitSub"), Form("cSlopeTemplateFitSub"), 500, 500);
  cSlopeTemplateFitSub->cd();
  hSlopeTemplateFit_sub->Draw();
      
  fitGausTemp = new TF1 ("fitGausTemp", "gaus", hSlopeTemplateFit_sub->GetMean()-hSlopeTemplateFit_sub->GetRMS()*2, hSlopeTemplateFit_sub->GetMean()+hSlopeTemplateFit_sub->GetRMS()*2);
  hSlopeTemplateFit_sub->Fit(fitGausTemp, "QR");
  float slope_template_sub = fitGausTemp->GetParameter(1);
  

  std::cout << "CTR sub template fit (sigma) = " << sigma_template_sub*1e3 <<" ps :: CTR intrinsic noise (sigma) temp sub = " << ped_noise[1]/slope_template_sub << " ps "  << " :: CTR sub temp without noise = " << sqrt(pow(sigma_template_sub*1e3,2) - pow(ped_noise[1]/slope_template_sub,2)) <<  std::endl;  
  
  
  
  TCanvas * cAmpTimeRawLin = new TCanvas ("cAmpTimeRawLin", "cAmpTimeRawLin", 600, 500);
  cAmpTimeRawLin->cd();
  pAmpTimeRawLin->Draw();
  pAmpTimeRawLin->GetXaxis()->SetTitle("amplitude raw lin");
  pAmpTimeRawLin->GetYaxis()->SetTitle("time [ns]");
  
  
  TCanvas * cAmpTimeSubLin = new TCanvas ("cAmpTimeSubLin", "cAmpTimeSubLin", 600, 500);
  cAmpTimeSubLin->cd();
  pAmpTimeSubLin->Draw();
  pAmpTimeSubLin->GetXaxis()->SetTitle("amplitude raw lin");
  pAmpTimeSubLin->GetYaxis()->SetTitle("time [ns]");
  
  TCanvas * cAmpTimeRawTemp = new TCanvas ("cAmpTimeRawTemp", "cAmpTimeRawTemp", 600, 500);
  cAmpTimeRawTemp->cd();
  pAmpTimeRawTemp->Draw();
  pAmpTimeRawTemp->GetXaxis()->SetTitle("amplitude raw lin");
  pAmpTimeRawTemp->GetYaxis()->SetTitle("time [ns]");
  
  TCanvas * cAmpTimeSubTemp = new TCanvas ("cAmpTimeSubTemp", "cAmpTimeSubTemp", 600, 500);
  cAmpTimeSubTemp->cd();
  pAmpTimeSubTemp->Draw();
  pAmpTimeSubTemp->GetXaxis()->SetTitle("amplitude raw lin");
  pAmpTimeSubTemp->GetYaxis()->SetTitle("time [ns]");
  
    
/*
      
  TCanvas * cCTR = new TCanvas ("cCTR", "cCTR", 500, 500);
  cCTR->cd();
  hCTR->Draw();
  hCTR->GetXaxis()->SetRangeUser(15, 30);
  
  TF1 * fitGaus = new TF1 ("fitGaus", "gaus", hCTR->GetMean()-hCTR->GetRMS()*2, hCTR->GetMean()+hCTR->GetRMS()*2 );
  hCTR->Fit(fitGaus, "QR");
  
  std::cout << "CTR (sigma) = " << fitGaus->GetParameter(2)*1e3 << " ps :: CTR (CTR) = " << fitGaus->GetParameter(2)*1e3*2.355 << std::endl;
  
  leg = new TLegend(0.56,0.66,0.88,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  
//   leg->AddEntry(gCTR_vs_BS_L, Form("left only"), "lpe");     
//   leg->AddEntry(gCTR_vs_BS_R, Form("right only"), "lpe");     
//   leg->AddEntry(gCTR_vs_BS,   Form("average"), "lpe");     
  leg->Draw("same");
  
//   hCTR_LED->SetLineColor(kGreen+1);
//   hCTR_LED->Draw("same");
  
  gPad->SetLogy();
  */
  
  
  TCanvas * cPulse[NPULSES_EX];
  
//   TLine * line_led_th[NCH];
//   TLine * line_cfd_th[NCH];
//   TLine * line_led_time[NCH];
//   TLine * line_cfd_time[NCH];
  
  TLine * int_gate_low[NCH];
  TLine * int_gate_up[NCH];
  
  for (int iPulse = 0; iPulse < NPULSES_EX; iPulse++)
  {
      
      cPulse[iPulse] = new TCanvas (Form("cPulse_%d", iPulse), Form("cPulse_%d", iPulse), 600, 500);
      cPulse[iPulse]->cd();

//       std::cout << "drawing single pulse SiPM... " << iPulse << std::endl;
      
      gPulse[iPulse][1]->SetTitle("Irradiated 3x3 mm^{2} SiPM - DCR = 20 GHz (Single event waveform)");
      gPulse[iPulse][1]->SetLineColor(kGreen+1);
      gPulse[iPulse][1]->SetLineWidth(2);
      gPulse[iPulse][1]->SetMarkerColor(kGreen+1);      
      gPulse[iPulse][1]->SetMarkerStyle(6);
      gPulse[iPulse][1]->Draw("ALP");
      gPulse[iPulse][1]->GetYaxis()->SetRangeUser(-0.05, 2.5);
      gPulse[iPulse][1]->GetXaxis()->SetRangeUser(0e-9, 40e-9);
      gPulse[iPulse][1]->GetXaxis()->SetTitle("Time [s]");
      gPulse[iPulse][1]->GetYaxis()->SetTitle("Amplitude [V]");
      gPulse[iPulse][1]->GetXaxis()->SetTitleSize(0.05);
      gPulse[iPulse][1]->GetYaxis()->SetTitleSize(0.05);
      gPulse[iPulse][1]->GetXaxis()->SetTitleOffset(0.92);
      gPulse[iPulse][1]->GetYaxis()->SetTitleOffset(0.95);
      
      
      
//       std::cout << "drawing single pulse LASER... " << iPulse << std::endl;
      
//       gPulse[iPulse][0]->SetMarkerStyle(20);
//       gPulse[iPulse][0]->Draw("same LPE");
      
      
//       std::cout << "drawing Sub SiPM pulse... " << iPulse << std::endl;
      gPulseSub[iPulse][1]->SetLineWidth(2);
      gPulseSub[iPulse][1]->Draw("same LP");
      gPulseSub[iPulse][1]->SetLineColor(kBlue+1);
      gPulseSub[iPulse][1]->SetMarkerColor(kBlue+1);
      gPulseSub[iPulse][1]->SetMarkerStyle(7);
      
      
//       fitPulseLED[iPulse][0]->SetLineColor(kRed);
//       fitPulseLED[iPulse][0]->Draw("same");
//       fitPulseLED[iPulse][1]->SetLineColor(kYellow+2);
//       fitPulseLED[iPulse][1]->Draw("same");
      
      /*
      std::cout << "drawing fit 1... " << iPulse << std::endl;
      if (fitPulseCFD[iPulse][0]!= NULL)
      {
          fitPulseCFD[iPulse][0]->SetLineColor(kMagenta);
          fitPulseCFD[iPulse][0]->Draw("same");
      }
      */
//       std::cout << "drawing fit 2... " << iPulse << std::endl;
      
      if (fitPulseCFD[iPulse][1]!= NULL)
      {
          fitPulseCFD[iPulse][1]->SetLineColor(kOrange+1);
          fitPulseCFD[iPulse][1]->Draw("same");
      }
//       std::cout << "drawing fit 3... " << iPulse << std::endl;
      if (fitPulseCFD_sub[iPulse][1]!= NULL)
      {
          fitPulseCFD_sub[iPulse][1]->SetLineColor(kCyan+1);
          fitPulseCFD_sub[iPulse][1]->Draw("same");
      }
      
//       float minTrange = -9.0e-9, maxTrange = 45.0e-9; 
      
//       std::cout << "done drawing waveforms..." << std::endl;

      double sampling = 100e-12;
      
      for (int iCh = 0; iCh < NCH; iCh++)
      {
//           line_cfd_th[iCh] = new TLine (minTrange, cfdTh[iCh], maxTrange, cfdTh[iCh]);
//           line_led_th[iCh] = new TLine (minTrange, ledTh[iCh], maxTrange, ledTh[iCh]);
//           line_led_th[iCh]->Draw("same");
          
//           int_gate_low[iCh] = new TLine ((float) sampling*(minS[iCh]-50), 0, (float) sampling*(minS[iCh]-50), 2);
//           int_gate_up[iCh] = new TLine ((float) sampling*(maxS[iCh]-50), 0, (float) sampling*(maxS[iCh]-50), 2);
          
          double minT = sampling*minS[iCh]+lowestT;
          double maxT = sampling*maxS[iCh]+lowestT;
          
          int_gate_low[iCh] = new TLine (minT, 0, minT, 2);
          int_gate_up[iCh] = new TLine (maxT, 0, maxT, 2);
          
//           std::cout << "minLine[" << iCh << "] : " << minT << " :: maxLine = " << maxT << std::endl;
          
          int_gate_low[iCh]->SetLineColor(iCh+1);
          int_gate_up[iCh]->SetLineColor(iCh+1);
          
          int_gate_low[iCh]->Draw("same");
          int_gate_up[iCh]->Draw("same");
      }
      
      
      leg = new TLegend(0.56,0.66,0.88,0.88,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);  
  
      leg->AddEntry(gPulse[iPulse][1], "Raw waveform", "lpe");   
      leg->AddEntry(gPulseSub[iPulse][1], "Processed waveform", "lpe");   
      leg->AddEntry(fitPulseCFD[iPulse][1], "Raw: leading edge fit", "lpe");   
      leg->AddEntry(fitPulseCFD_sub[iPulse][1], "Processed: leading edge fit", "lpe");   
      leg->AddEntry(fitRiseDecay, "Raw: template fit", "lpe");
      leg->AddEntry(fitRiseDecaySub, "Processed: template fit", "lpe");
      
      leg->Draw("same");
      
//       gPad->SetLogy();
  }

  std::cout << "ready to write output file " << std::endl;
//   
  
  TFile * outputFile = new TFile(output_file_name.c_str(), "RECREATE"); //DCR_S12572-015C-2percent-5VOV-LED0.root
//   TFile * outputFile = new TFile(Form("./output/out_laser2_irr_%d.root", bias), "RECREATE");
  outputFile->cd();
  
  for (int iCh = 0; iCh <NCH; iCh++)
  {
      
//     hTimeCFD[iCh]->Write();
//     hSlopeCFD[iCh]->Write();
    hInt[iCh]           ->Write();
    hAmp[iCh]           ->Write();
    hPedSample[iCh]     ->Write();
    hPedSampleSub[iCh]  ->Write();
     
    gSlope_vs_Th[iCh]      ->SetName(Form("gSlope_vs_Th_%d", iCh)      );
    gSlopeSub_vs_Th[iCh]   ->SetName(Form("gSlopeSub_vs_Th_%d", iCh)   );
    gTimeCFD_vs_Th[iCh]    ->SetName(Form("gTimeCFD_vs_Th_%d_", iCh)   );
    gTimeCFDSub_vs_Th[iCh] ->SetName(Form("gTimeCFDSub_vs_Th_%d", iCh) );
    gTimeEffCFD_vs_Th[iCh]    ->SetName(Form("gTimeEffCFD_vs_Th_%d_", iCh)   );
    gTimeEffCFDSub_vs_Th[iCh] ->SetName(Form("gTimeEffCFDSub_vs_Th_%d", iCh) );
    
    gSlope_vs_Th[iCh]      ->Write();
    gSlopeSub_vs_Th[iCh]   ->Write();
    gTimeCFD_vs_Th[iCh]    ->Write();
    gTimeCFDSub_vs_Th[iCh] ->Write();
    gTimeEffCFD_vs_Th[iCh]    ->Write();
    gTimeEffCFDSub_vs_Th[iCh] ->Write();
    
//     hTimeCFD_sub[iCh]->Write();
//     hSlopeCFD_sub[iCh]->Write();
    
    for (int iPulse = 0; iPulse < NPULSES_EX; iPulse++)
    {
        gPulse[iPulse][iCh]->SetName(Form("gPulse_%d_%d", iPulse, iCh) );
        gPulseSub[iPulse][iCh]->SetName(Form("gPulseSub_%d_%d", iPulse, iCh) );
        
        gPulse[iPulse][iCh]->Write();        
        gPulseSub[iPulse][iCh]->Write();
    }
    
  }
  
  hChisquareLin->Write();
  hChisquareLinSub->Write();
  hChisquareTemp->Write();
  hChisquareTempSub->Write();
  
  hTimeTemplateFit->Write();
  hSlopeTemplateFit->Write();
  
  hTimeTemplateFit_sub->Write();
  hSlopeTemplateFit_sub->Write();
  
//   gPeaks->SetName("gPeaks");
//   gPeaks->Write();
  
  outputFile->Write();
  outputFile->Close();
  
  std::cout << "done writing output! -->  " << output_file_name.c_str() << std::endl;
 
  return 0;
//   std::exit(0);
//    theApp->Run();
  
  
}









