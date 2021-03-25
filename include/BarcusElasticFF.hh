#ifndef BARCUS_ELASTIC_FORM_FACTOR_HH
#define BARCUS_ELASTIC_FORM_FACTOR_HH

// 3He elastic form factors
// - Reference: S. Barcus Ph.D. Thesis (College of William and Mary, 2019)     

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

#define SB_SIZE 12 

class BarcusElasticFF: public ElasticFormFactor { 

   private:
      double fR[SB_SIZE],fQch[SB_SIZE],fQmag[SB_SIZE];

      double SumOfGaussians(int type,double q); 

   public: 
      BarcusElasticFF(int target=RC::Target::kHe3);
      ~BarcusElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
