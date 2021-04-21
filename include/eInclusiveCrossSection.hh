#ifndef eINCLUSIVE_CROSS_SECTION_HH
#define eINCLUSIVE_CROSS_SECTION_HH

// eInclusiveCrossSection class
// Notes: - Abstract base class 
//        - Variables: beam energy (Es) and scattered energy (Ep) are in GeV, 
//          scattering angle (th) is in degrees.  A is in g/mol. 
//        - Needed values: Z, A, Es, Ep, th 
//        - The resulting cross section is for inclusive inelastic electron scattering  

#include <cstdlib> 
#include <iostream>
#include <iomanip> 
#include <cmath>

#include "RCConstants.hh"
#include "Target.hh"

class eInclusiveCrossSection {

   protected: 
      double fZ,fA;
      double fEs,fEp,fTh;
      
      Target *fTgt;  

      void Init();

   public: 
      eInclusiveCrossSection();
      ~eInclusiveCrossSection();
 
      void SetTarget(Target *t); 
      void SetKinematicVariables(double Es,double Ep,double th); 

      void SetZ(double v)  { fZ  = v;    }
      void SetA(double v)  { fA  = v;    }
      void SetEs(double v) { fEs = v;    } 
      void SetEp(double v) { fEp = v;    }
      void SetTh(double v) { fTh = v;    }
     
      double GetMottXS(double,double);
      double GetZ()  const { return fZ;  }
      double GetA()  const { return fA;  }
      double GetEs() const { return fEs; }
      double GetEp() const { return fEp; }
      double GetTh() const { return fTh; }

      virtual double GetBornXS()=0;

};

#endif
