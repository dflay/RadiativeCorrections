#ifndef INCLUSIVE_E_RADIATIVE_CORRECTIONS_HH
#define INCLUSIVE_E_RADIATIVE_CORRECTIONS_HH

// a class for applying radiative effects to cross section models, 
// or unfolding radiative effects from experimental cross section data 

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "constants.hh"
#include "Kinematics.hh"
#include "eInclusiveCrossSection.hh" 
#include "ElasticFormFactor.hh"

namespace RC {
   enum thrType_t { kElastic, kPion }; 
} 

class RadiativeCorrections {

   private:
      int fVerbosity; 

      RC::thrType_t fThreshold;  // integration threshold: elastic or pion 

      double fZ,fA;
      double fDeltaE;
      double fMT;
      double fEs,fEp,fThDeg;
      double fR,fCFACT;
      double fTa,fTb,fT,fEta,fXi,fb;

      void CalculateB();
      void CalculateXi();
      void CalculateEta();
      void CalculateR();
      void CalculateCFACT();

      double CalculateEsIntegral();
      double CalculateEpIntegral();

      // For the Radiate method
      double GetTr(double);           // Tr(Q2): changes PER INTEGRATION BIN 
      double GetFTilde(double);       // FTilde(Q2): changes PER INTEGRATION BIN 
      double GetPhi(double);          // phi(v), v = arbitrary value  
      double GetEsMin(double);        // EsMin(Ep) 
      double GetEpMax(double);        // EpMax(Es)                        
      double GetSpence(double);       // Spence(x), x = arbitrary value 

      double EsIntegrand(const double);
      double EpIntegrand(const double);
      double Integrate(double (RadiativeCorrections::*)(const double),double,double,double,int);
      double AdaptiveSimpsonAux(double (RadiativeCorrections::*)(const double),double,double,double,double,double,double,double,int);

      // elastic tail 
      double GetF_soft(); 
      double GetWs(double,double,double); 
      double GetWp(double,double,double);
      double sigma_el(double); 
      double ElasticTail_sigmaEx_Integrand(const double);    

      eInclusiveCrossSection *fInclXS;
      ElasticFormFactor *fFormFactor;

   public:
      RadiativeCorrections();
      ~RadiativeCorrections();

      void Init();
      void Print();
     
      // use these for testing only 
      void SetTargetVariables(double Z,double A)       { fZ = Z; fA = A; fMT = A*proton_mass; }  
      void SetKinematicVariables(double Es,double Ep,double th); 
      void CalculateVariables();  // compute various variables when Es, Ep, th change

      void SetIntegrationThreshold(RC::thrType_t t)    { fThreshold = t; } 

      void SetVerbosity(int v)  { fVerbosity = v; } 

      void SetTb(double tb) { fTb = tb; }
      void SetTa(double ta) { fTa = ta; }
      void SetRadiationLengths(double tb,double ta)    { fTb = tb; fTa = ta; fT = tb + ta; } 

      void SetCrossSection(eInclusiveCrossSection *XS) { fInclXS     = XS; }
      void SetFormFactor(ElasticFormFactor *ff)        { fFormFactor = ff; }

      double Radiate();

      double ElasticTail_peakApprox();  
      double ElasticTail_exact();  

      double ElasticTail_sigmaP();    
      double ElasticTail_sigmaB();    
      double ElasticTail_sigmaEx();    

}; 

#endif 
