#ifndef STEIN_ELASTIC_FORM_FACTOR_HH
#define STEIN_ELASTIC_FORM_FACTOR_HH

// elastic form factors for Al, Cu, Au, C,... 
// - References: Phys. Rev. D 12, 1884 (1975)
//              Rev. Mod. Phys. 28, 214 (1956)      

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

class SteinElasticFF: public ElasticFormFactor { 

   private:

   public: 
      SteinElasticFF(int target=RC::Target::kProton);
      ~SteinElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
