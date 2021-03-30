#ifndef CARBON_ELASTIC_FORM_FACTOR_HH
#define CARBON_ELASTIC_FORM_FACTOR_HH

// Carbon elastic form factors
// - Based on Sum of Gaussians fit 
// - Atomic data: 
//   - Elem C  
//   - Z    6  
//   - A    12 
//   - RMS  2.469      
//   - RP   1.2 
// - Reference: http://discovery.phys.virginia.edu/research/groups/ncd/download.html     

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

#define C_SIZE 11 

class CarbonElasticFF: public ElasticFormFactor { 

   private:
      double fR_rms; 
      double fR[C_SIZE],fQch[C_SIZE],fQmag[C_SIZE];

      double SumOfGaussians(int type,double q); 

   public: 
      CarbonElasticFF(int target=RC::Target::kCarbon);
      ~CarbonElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
