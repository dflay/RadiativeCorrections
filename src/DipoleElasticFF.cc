#include "DipoleElasticFF.hh"
//______________________________________________________________________________
DipoleElasticFF::DipoleElasticFF(){

}
//______________________________________________________________________________
DipoleElasticFF::~DipoleElasticFF(){

}
//______________________________________________________________________________
double DipoleElasticFF::GetGE(double Q2){
   // Phys. Rev. D 12, 1884 for details (A7 and A8)   
   // assumes Q2 is in GeV^2 
   const int N = 6;
   double H[N] = {1.0007,1.01807,1.05584,0.836380,0.6864584,0.672830}; 
   double P=0;
   double prod_sum=1;
   double q = sqrt(Q2); 
   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	 if(i!=j) prod_sum *= (q-j)/(i-j);  
      }
      P += H[i]*prod_sum;
      prod_sum = 1;
   }
   double GE = P/pow(1. + Q2/0.71,2.); 
   return GE; 
}
//______________________________________________________________________________
double DipoleElasticFF::GetGM(double Q2){
   double GE = GetGE(Q2); 
   double GM = (1. + RC::Constants::kappa_p)*GE;
   return GM; 
}
