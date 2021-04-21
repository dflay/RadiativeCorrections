#ifndef RADIATIVE_CORRECTIONS_TARGET_HH
#define RADIATIVE_CORRECTIONS_TARGET_HH

// target material for radiative corrections 

#include <cstdlib> 
#include <iostream>
#include <cmath> 

#include "RCConstants.hh"

class Target { 
   private: 
      std::string fName;
      int fType;           // target type (see RCConstants.hh for choices) 
      double fA;           // atomic mass [g/mol]
      double fZ;           // atomic number (number of protons) 
      double fN;           // number of neutrons 
      double fMT;          // target mass [GeV] 
      double fRho;         // density [g/cm^3]  
      double fTb,fTa;      // thickness before and after scattering [#X0]   

   public: 
      Target(int type=-1,double A=0,double Z=0,double rho=0);
      ~Target();

      void Print(); 
      void SetA(double a)     { fA   = a; fMT = fA*RC::Constants::proton_mass; } 
      void SetZ(double z)     { fZ   = z;   } 
      void SetN(double n)     { fN   = n;   } 
      void SetRho(double rho) { fRho = rho; }
      void SetMass(double m)  { fMT  = m;   } 
      void SetRadiationLengths(double tb,double ta) { fTb = tb; fTa = ta; }  

      double GetA()     const { return fA;   }  
      double GetZ()     const { return fZ;   }  
      double GetN()     const { return fN;   }  
      double GetRho()   const { return fRho; }
      double GetMass()  const { return fMT;  }  
      double GetTb()    const { return fTb;  }  
      double GetTa()    const { return fTa;  } 
}; 

#endif 
