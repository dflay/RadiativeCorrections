#ifndef DIPOLE_ELASTIC_FORM_FACTOR_HH
#define DIPOLE_ELASTIC_FORM_FACTOR_HH

// dipole model of the elastic form factor for the proton
// - see Phys. Rev. D 12, 1884 for details (A7 and A8)   

#include <cstdlib>
#include <iostream>

#include "RCConstants.hh"
#include "ElasticFormFactor.hh"

class DipoleElasticFF: public ElasticFormFactor { 

   private:

   public: 
      DipoleElasticFF();
      ~DipoleElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
