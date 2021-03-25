#ifndef AMROUN_ELASTIC_FORM_FACTOR_HH
#define AMROUN_ELASTIC_FORM_FACTOR_HH

// 3He or 3H elastic form factors
// - Reference: Nucl. Phys. A, 579 596 (1994)    

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

#define AM_SIZE 12 

class AmrounElasticFF: public ElasticFormFactor { 

   private:
      double fR[AM_SIZE],fQch[AM_SIZE],fQmag[AM_SIZE];

      double SumOfGaussians(int type,double q); 

   public: 
      AmrounElasticFF(int target=RC::Target::kHe3);
      ~AmrounElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
