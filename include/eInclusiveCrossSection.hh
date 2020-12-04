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

#define ELECTRON_MASS 5.11e-4 // in GeV
#define PROTON_MASS 0.938     // in GeV 
#define PION_MASS 0.140       // in GeV 
#define PI 3.14159265359
#define ALPHA 1./137. 
#define DEG_TO_RAD PI/180.
#define HBAR_C 624.4197  // in GeV*nb^(1/2) 

using namespace std; 

class eInclusiveCrossSection{

   protected: 
      double fZ,fA;
      double fEs,fEp,fTh; 	

      void Init();

   public: 
      eInclusiveCrossSection();
      ~eInclusiveCrossSection();

      void SetZ(double v){fZ = v;}
      void SetA(double v){fA = v;}
      void SetEs(double v){fEs = v;} 
      void SetEp(double v){fEp = v;}
      void SetTh(double v){fTh = v;} 

      double GetMottXS(double,double);
      double GetZ(){return fZ;}
      double GetA(){return fA;}
      double GetEs(){return fEs;}
      double GetEp(){return fEp;}
      double GetTh(){return fTh;}

      virtual double GetBornXS()=0;

};

#endif
