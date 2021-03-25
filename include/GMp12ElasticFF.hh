#ifndef GMP12_ELASTIC_FORM_FACTOR_HH
#define GMP12_ELASTIC_FORM_FACTOR_HH

// elastic form factors for the proton 
// - fit from GMp(12) preprint  

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

class GMp12ElasticFF: public ElasticFormFactor { 

   private:

   public: 
      GMp12ElasticFF();
      ~GMp12ElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
