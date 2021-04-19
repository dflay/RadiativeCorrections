#ifndef RADIATIVE_CORRECTIONS_HH
#define RADIATIVE_CORRECTIONS_HH

// a class for applying radiative effects to cross section models, 
// or unfolding radiative effects from experimental cross section data 

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector> 

#include "RCConstants.hh"
#include "Kinematics.hh"
#include "eInclusiveCrossSection.hh" 
#include "ElasticFormFactor.hh"

namespace RC {
   enum thrType_t  { kElastic, kPion }; 
   enum unitType_t { kMicrobarnPerGeVPerSr, kNanobarnPerGeVPerSr, kPicobarnPerGeVPerSr }; 
} 

class RadiativeCorrections {

   private:
      int fVerbosity;
      int fNumIter; 
      bool fElasticTail,fElasticApprox; 

      RC::thrType_t fThreshold;  // integration threshold: elastic or pion
      RC::unitType_t fUnit;  

      double fZ,fA;
      double fDeltaE;
      double fMT;
      double fEs,fEp,fThDeg;
      double fR,fCFACT;
      double fTa,fTb,fT,fEta,fXi,fb;
      double fIterThresh; 

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
      double GetWs(double,double,double); 
      double GetWp(double,double,double);
      double sigma_el(double Es); 
      double sigma_el_tilde(double Es); 
      // double ElasticTail_sigmaEx_Integrand(const double);   

      // elastic peak kinematic factors 
      double GetX(double Es,double th); 
      double GetRho(double Es,double th); 
      double GetEta(double Es,double th); 
      double GetEta_MTS(double Es,double th); 

      eInclusiveCrossSection *fInclXS;
      ElasticFormFactor *fFormFactor;

   public:
      RadiativeCorrections();
      ~RadiativeCorrections();

      void Init();
      void Print();
  
      void ElasticTailEnable(bool v=true)              { fElasticTail   = v; } 
      void UseElasticTailApprox(bool v=true)           { fElasticApprox = v; } 

      void SetUnits(RC::unitType_t u)                  { fUnit = u; };
      void SetIntegrationThreshold(RC::thrType_t t)    { fThreshold = t; } 
      void SetVerbosity(int v)                         { fVerbosity = v; } 

      void SetTb(double tb) { fTb = tb; }
      void SetTa(double ta) { fTa = ta; }
      void SetRadiationLengths(double tb,double ta)    { fTb = tb; fTa = ta; fT = tb + ta; } 

      void SetCrossSection(eInclusiveCrossSection *XS) { fInclXS     = XS; }
      void SetFormFactor(ElasticFormFactor *ff)        { fFormFactor = ff; }

      // settings for unfolding born xs 
      void SetNumberOfIterations(int i)                { fNumIter    = i;   }  // number of iterations 
      void SetUnfoldingTolerance(double thr)           { fIterThresh = thr; }  // convergence threshold  

      int Unfold(double Es,double th,std::vector<double> Ep,std::vector<double> xsr,std::vector<double> &xsb); 

      double Radiate();
      double ElasticTail_peakApprox();  
      double ElasticTail_exact();  
     
      // use these for testing only 
      void SetTargetVariables(double Z,double A)       { fZ = Z; fA = A; fMT = A*RC::Constants::proton_mass; }  
      void SetKinematicVariables(double Es,double Ep,double th); 
      void CalculateVariables();  // compute various variables when Es, Ep, th change
      double GetF_soft(); 
    
      // for testing.  will be private once things are finalized 
      double ElasticTail_sigmaP();    
      double ElasticTail_sigmaB();    
      double ElasticTail_sigmaEx();   
      double ElasticTail_sigmaEx_Integrand(const double);   

      // Mo & Tsai 
      double ElasticPeak_Delta_MTS();  
      double ElasticPeak_Z0_MTS();  
      double ElasticPeak_Z1_MTS();  
      double ElasticPeak_Z2_MTS();  

      // Meister and Yennie  
      double ElasticPeak_Delta_MY();  
      double ElasticPeak_Z0_MY();  
      double ElasticPeak_Z1_MY();  
      double ElasticPeak_Z2_MY();  
      
      // Maximon and Tjon  
      double ElasticPeak_Delta_MTJ_E1THDE(double Es,double th,double deltaE);  
      double ElasticPeak_Z0_MTJ(double Es,double th,double deltaE);  
      double ElasticPeak_Z1_MTJ(double Es,double th,double deltaE);  
      double ElasticPeak_Z2_MTJ(double Es,double th,double deltaE);  
      double ElasticPeak_DeltaEl_MTJ();  

}; 

#endif 
