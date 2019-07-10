#include "functions.h"
#include "TSpectrum.h"



float DummyFunc()
{
     return 1;
}


std::pair<std::vector<float>*,std::vector<float>*>  GetInvSubWaveform(std::vector<float>* times, std::vector<float>* amps, float& delay)
{
    
   std::vector<float> *new_times = new std::vector<float> ();
   std::vector<float> *new_amps  = new std::vector<float> ();
   
   float samplingRate = 100;    //ps
   unsigned int delayIt = (int) delay/samplingRate;
//    std::cout << "delay = " << delay << " :: delayIt = " << delayIt << std::endl;
    
  for(unsigned int sIt = 0; sIt < times->size(); ++sIt)
  {
      if (sIt > delayIt)
      {        
        new_times->push_back(times->at(sIt));
        new_amps->push_back( amps->at(sIt) - amps->at(sIt-delayIt) );
//         std::cout << "new_amp[" << sIt << "] = " <<  amps->at(sIt) - amps->at(sIt-delayIt)  << " :: delay = " << delay << " :: delayIt = " << delayIt << std::endl;
      }      
  }
  
  std::pair<std::vector<float>*,std::vector<float>*> ret(new_times, new_amps);
  return ret;
}

float GetPedestal(std::vector<float>* times, std::vector<float>* amps, const unsigned int& s1, const unsigned int& s2)
{
  float ped = 0.;
  
  for(unsigned int sIt = s1; sIt < s2; ++sIt)
  {
    ped += amps->at(sIt);
  }
  
  return 1. * ped / (s2-s1);
}



std::pair<float,float>  GetAmplitude(std::vector<float>* times, std::vector<float>* amps, const unsigned int& s1, const unsigned int& s2, float& polarity, const unsigned int& nAve)
{
  float ped = GetPedestal(times, amps, s1, s2)*polarity;
  
  float max_amp = -999999999.;    
  float time_max = -9999999;    
  int mySit = 0; //number of samples used for averaging the max
  
  //find absolute max amp point
  for(unsigned int sIt = s2+6; sIt < times->size()-6; ++sIt)
  {
      if (amps->at(sIt)*polarity > max_amp )
      {
        max_amp = amps->at(sIt)*polarity;
        time_max = times->at(sIt);
        mySit = sIt;
      }
  }
  
  //average across neighboring points  
  
  int nPointsAverage = nAve;
  max_amp = 0;
  for (int iPoint = 0; iPoint < nPointsAverage; iPoint++)
  {
      max_amp += amps->at(mySit - (nPointsAverage-1)/2 + iPoint)*polarity;
  }
  max_amp /= nPointsAverage;
  
  
  
  //gaussian fit on neighboring points
  /*
  TGraph * g = new TGraph();
  int nPointsAverage = nAve;
  max_amp = 0;
  for (int iPoint = 0; iPoint < nPointsAverage; iPoint++)
  {
      g->SetPoint(iPoint, iPoint, amps->at(mySit - (nPointsAverage-1)/2 + iPoint)*polarity);      
  }
  TF1 * fitGaus = new TF1 ("fitGaus", "gaus", -6, 6); 
  fitGaus->SetParameters(amps->at(mySit)*polarity, mySit, nPointsAverage);
  g->Fit(fitGaus, "QR");
  max_amp = fitGaus->GetMaximum();
  delete g;
  delete fitGaus;
  */
  
  std::pair<float,float> ret(max_amp-ped,time_max);
//   g->Delete();
//   fitGaus->Delete();

  return ret;
}


float GetIntegral(std::vector<float>* times, std::vector<float>* amps, const unsigned int& s1, const unsigned int& s2, float& polarity, const unsigned int& minS, const unsigned int& maxS)
{
  float ped = GetPedestal(times, amps, s1, s2);    
  float integral = 0;
  unsigned int myMinS = minS;
  unsigned int myMaxS = maxS;
  
  if (minS == (unsigned int) -1) myMinS = s2;
  if (maxS == (unsigned int) -1) myMaxS = times->size();
  
  for(unsigned int sIt = myMinS; sIt < myMaxS; ++sIt)
  {
      integral += (amps->at(sIt) - ped)*polarity;
  }
    
  return integral;
}


double riseSingleDecay(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double t = x[0];
  
  double rho 		= par[0];
  double tau_r 		= par[1];
  double tau_d 		= par[2];
  double delta_m 	= par[3];
  double theta	 	= par[4];
  double sigma_rf 	= par[5];
  

  double f = (1/(2*(tau_d - tau_r)) * exp((2*tau_d*(delta_m + theta - t) + pow(sigma_rf,2)) / (2*pow(tau_d,2))) *(1-erf( (tau_d*(delta_m+theta-t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_d) )))
	   - (1/(2*(tau_d - tau_r)) * exp((2*tau_r*(delta_m + theta - t) + pow(sigma_rf,2)) / (2*pow(tau_r,2))) *(1-erf( (tau_r*(delta_m+theta-t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_r) )));
  
  return rho*f;
}



double riseSingleDecaySub(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double t = x[0];
  
  double rho 		= par[0];
  double tau_r 		= par[1];
  double tau_d 		= par[2];
  double delta_m 	= par[3];
  double theta	 	= par[4];
  double sigma_rf 	= par[5];  
  

  double f = (1/(2*(tau_d - tau_r)) * exp((2*tau_d*(theta - t) + pow(sigma_rf,2)) / (2*pow(tau_d,2))) *(1-erf( (tau_d*(theta-t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_d) )))
	   - (1/(2*(tau_d - tau_r)) * exp((2*tau_r*(theta - t) + pow(sigma_rf,2)) / (2*pow(tau_r,2))) *(1-erf( (tau_r*(theta-t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_r) )))
           - ( (1/(2*(tau_d - tau_r)) * exp((2*tau_d*(delta_m + theta - t) + pow(sigma_rf,2)) / (2*pow(tau_d,2))) *(1-erf( (tau_d*(delta_m+theta -t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_d) )))
	   - (1/(2*(tau_d - tau_r)) * exp((2*tau_r*(delta_m + theta - t) + pow(sigma_rf,2)) / (2*pow(tau_r,2))) *(1-erf( (tau_r*(delta_m+theta -t) + pow(sigma_rf,2))  / (sqrt(2)*sigma_rf*tau_r) ))) )        ;
  
  return rho*f;
}



double riseSingleDecayNoIRF(double* x, double* par)
{
  
  double t = x[0];
  
  double rho 		= par[0];
  double tau_r 		= par[1];
  double tau_d 		= par[2];
  double delta_m 	= par[3];
  double theta	 	= par[4];


  if (t<theta+delta_m)
  {
    return 0;
  }
  
  else
  {
    return rho * (exp((delta_m + theta -t) / tau_d) - exp((delta_m + theta - t) / tau_r)) / (tau_d - tau_r);
  }
  
  
}

double riseSingleDecayNoIRFSub(double* x, double* par)
{
  
  double t = x[0];
  
  double rho 		= par[0];
  double tau_r 		= par[1];
  double tau_d 		= par[2];
  double delta_m 	= par[3];
  double theta	 	= par[4];


  if (t<theta+delta_m)
  {
    return 0;
  }
  
  else
  {
    return rho * (exp((delta_m + theta -t) / tau_d) - exp((delta_m + theta - t) / tau_r)) / (tau_d - tau_r);
  }
  
  
}

/*
std::pair<float,float> GetAmplitudeSquare(float* vals, const unsigned int& s1, const unsigned int& s2, const float& SqTH, const unsigned int& nS, const bool& isNegative)
{
  float ped = GetPedestal(vals,s1,s2);
  float threshold = SqTH*1000;
  int sMin = 0;
  
//   threshold = 100;
  
  for(unsigned int sIt = 0; sIt < 1024; ++sIt)
  {
     if( isNegative)
     {
//        std::cout << "vals = " << vals[sIt] << " :: ped = "  << ped << " :: threshold = " << threshold << std::endl;
       if( fabs(vals[sIt]-ped) > threshold )
       {
         sMin = sIt;         
         break;
       }
     }
  }
  
  sMin+=10;
  
  float amp = GetPedestal(vals,sMin,sMin+nS);
  //std::cout << "sMin = " << sMin << " :: amp = " << amp << std::endl;
  
  std::pair<float,float> ret(ped,fabs(amp-ped));
  return ret;
}
*/


std::pair<float, TF1*> GetTime(std::vector<float>* times, std::vector<float>* amps, const unsigned int& s1, const unsigned int& s2,  
                                                                float& polarity, const float& th, 
                                                                const unsigned int& minS, const unsigned int& maxS,
                                                                const unsigned int& nSBef, const unsigned int& nSAft,
                                                                const bool& isFraction, const bool& deleteFit, const unsigned int& nAveAmp)
{
  
    
  float ped = GetPedestal(times, amps, s1, s2)*polarity;    
  float max_amp = GetAmplitude(times, amps, s1, s2, polarity, nAveAmp).first;
    
    
  float threshold = th;
  if( isFraction )  threshold = max_amp * th;
  
  
  unsigned int sBef = 0;
  
  unsigned int myMinS = minS;
  unsigned int myMaxS = maxS;
  
  if (minS == (unsigned int) -1) myMinS = s2;
  if (maxS == (unsigned int) -1) myMaxS = times->size();
  
  for(unsigned int sIt = myMinS; sIt < myMaxS; ++sIt)
  {
     if( (amps->at(sIt)*polarity - ped) > threshold )
     {
         sBef = std::max(0, int(sIt-1) );
         break;
     }           
   }
   
//    std::cout << "sBef1 = " << sBef1 << " :: sBef2 = " << sBef2 << std::endl;
  
  std::pair<float,TF1*> ret;
  
  
  if( sBef > 0 )
  {
    TGraph* g = PointsAroundTh(times, amps, polarity, sBef, nSBef, nSAft);
    g -> SetMarkerColor(kBlue+1);
    g -> SetMarkerStyle(21);
    g -> GetXaxis()->SetLimits(0,1024);
    g -> GetYaxis()->SetLimits(0,1);
    g -> Fit("pol1","QR");
    //g -> Draw("ALPE");
    
    TF1* fitFunc = (TF1*)( g->GetFunction("pol1") );
    fitFunc -> SetLineColor(kMagenta);
    fitFunc -> SetLineWidth(1);
    
    float time = -1.;
    time = ( (ped+threshold*polarity) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1);
    
    ret.first  = time;
    ret.second = fitFunc;
    
    if( deleteFit )
    {
      delete fitFunc;
      delete g;
    }
  }
  
  
  return ret;
}


TGraph* PointsAroundTh(std::vector<float>* times, std::vector<float>* amps, float& polarity, const unsigned int& sBef, const unsigned int& nSBef, const unsigned int& nSAft)
{
  // n-point interpolation: nSBef points before threshold, BSAft points after threshold
  TGraph* g = new TGraph();

  int point = 0;
  for(unsigned int iS = 0; iS < nSBef; ++iS)
  {
    unsigned int sample = sBef - (nSBef-iS-1);        
    g -> SetPoint(point, times->at(sample), amps->at(sample)*polarity);
    ++point;
  }
  for(unsigned int iS = 0; iS < nSAft; ++iS)
  {
    unsigned int sample = sBef + 1 +iS;        
    g -> SetPoint(point, times->at(sample), amps->at(sample)*polarity);
    ++point;
  }    
  
  return g;
  delete g;
  
}

