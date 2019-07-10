#include "VectorSmallestInterval.hh"


double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  double delta = 999999.;
  
//   if (vals->size()>0)
  {
    std::sort(vals->begin(),vals->end());

    unsigned int nPoints = vals->size();
    unsigned int maxPoints = (unsigned int)(fraction * nPoints);

//   unsigned int minPoint = 0;
//   unsigned int maxPoint = 0;
    
    for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
    {
        double tmpMin = vals->at(point);
        double tmpMax = vals->at(point+maxPoints-1);
        if( tmpMax-tmpMin < delta )
        {
            delta = tmpMax - tmpMin;
            min = tmpMin;
            max = tmpMax;
//       minPoint = point;
//       maxPoint = point + maxPoints - 1;
        }
    }
  }

  
  return delta;
  
}
