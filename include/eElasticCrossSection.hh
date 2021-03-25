#ifndef eELASTIC_CROSS_SECTION_HH
#define eELASTIC_CROSS_SECTION_HH

// eElasticCrossSection class
// Notes: - Kinematic variables: beam energy (Es) and scattering angle (th) is in degrees.   
//        - Target variables: Z, A 
//        - The resulting cross section is for elastic electron scattering  

#include <cstdlib> 
#include <iostream>
#include <iomanip> 
#include <cmath>

#include "RCConstants.hh"
#include "Kinematics.hh"
#include "ElasticFormFactor.hh"
#include "eInclusiveCrossSection.hh"

class eElasticCrossSection: public eInclusiveCrossSection  {

   private: 
      ElasticFormFactor *fFormFactor; 	

   public: 
      eElasticCrossSection(double Z=1,double A=1);
      ~eElasticCrossSection();

      void SetTargetParameters(double Z,double A){ fZ = Z; fA = A; }
      
      void SetFormFactor(ElasticFormFactor *ff) { fFormFactor = ff; }
       
      double GetBornXS();

};

#endif
