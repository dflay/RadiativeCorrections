#ifndef ELASTIC_FORM_FACTOR_HH
#define ELASTIC_FORM_FACTOR_HH

// Elastic form factor class
// Notes: - Abstract base class 
//        - User provides a model deriving from this to determine GE and GM as a function of Q2   

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "RCConstants.hh"

class ElasticFormFactor {

   protected:
      int fType;     // target type 
      double fA,fZ;  // target A and Z 
      double fMT;    // target mass

   public:
      ElasticFormFactor()  {};
      ~ElasticFormFactor() {};

      void SetTarget(int type); 

      virtual double GetGE(double Q2)=0;
      virtual double GetGM(double Q2)=0;

};

#endif 
