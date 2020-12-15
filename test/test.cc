// test radiative corrections 

#include <cstdlib>
#include <iostream>
#include <vector>

#include "ElasticFormFactor.hh"
#include "DipoleElasticFF.hh"
#include "RadiativeCorrections.hh"

#include "CSVManager.hh"

int main(){

   double Es=0;
   std::cout << "Enter beam energy (1, 5, 20 GeV): ";
   std::cin >> Es;

   // from Phys. Rev. D 12, 1884 (1975) 
   double tb    = 5.459E-3;
   double ta    = 11.129E-3;
   double th    = 5.;    // in deg 
   double EpMin = 0.100; // in GeV
   double EpMax = 4.901; // in GeV 
   double Z     = 1;
   double A     = 1;

   int ES_INT = (int)Es;

   if( ES_INT==1 ){
      EpMin = 0.100;
      EpMax = 0.99;
   }else if( ES_INT==5 ) {
      EpMin = 0.100;
      EpMax = 4.901;
   }else if( ES_INT==20 ){
      EpMin = 1;
      EpMax = 18.6;
   }

   double xMin = EpMin*0.8;
   double xMax = EpMax*1.1;

   DipoleElasticFF *myDip = new DipoleElasticFF();

   RadiativeCorrections *myRC = new RadiativeCorrections();
   myRC->SetFormFactor(myDip);
   myRC->SetTargetVariables(Z,A);
   // myRC->SetRadiationLengths(tb,ta);
   // myRC->SetVerbosity(1);  

   char msg[200]; 

   const int N = 100;
   double min=1E-7,max=0;
   std::vector<double> EP,ET,ETA;
   double Ep=0,el_tail=0,el_tail_approx=0;
   double step = (EpMax - EpMin)/( (double)N );
   for(int i=0;i<N;i++){
      Ep = EpMin + ( (double)i )*step;
      myRC->SetKinematicVariables(Es,Ep,th);
      el_tail        = MUB_PER_GEV_SR*myRC->ElasticTail_exact();      // need to convert to mub/GeV/sr! 
      el_tail_approx = MUB_PER_GEV_SR*myRC->ElasticTail_peakApprox(); // need to convert to mub/GeV/sr!
      if(el_tail<min)        min = el_tail;
      if(el_tail_approx<min) min = el_tail_approx;
      if(el_tail>max)        max = el_tail;
      if(el_tail_approx<max) max = el_tail_approx;
      EP.push_back(Ep);
      ET.push_back(el_tail);
      ETA.push_back(el_tail_approx);
      sprintf(msg,"Es = %.3lf GeV, Ep = %.3lf GeV, th = %.1lf deg, el_tail = %.3E, el_tail_approx = %.3E",
              Es,Ep,th,el_tail,el_tail_approx);
      std::cout << msg << std::endl; 
   }

   const int NCOL = 3;
   const int NROW = EP.size();

   char outpath[200]; 
   sprintf(outpath,"./output/rc-test_Es-%d.csv",ES_INT); 

   util_df::CSVManager *myCSV = new util_df::CSVManager();
   myCSV->InitTable(NROW,NCOL);
   myCSV->SetHeader("Ep(GeV),tail(mub/GeV/sr),tail_approx(mub/GeV/sr)"); 
   myCSV->SetColumn<double>(0,EP ); 
   myCSV->SetColumn<double>(1,ET ); 
   myCSV->SetColumn<double>(2,ETA);
   myCSV->WriteFile(outpath);  

   delete myDip;
   delete myRC;
   delete myCSV;

   return 0;
}
