#include "GMp12ElasticFF.hh"
//______________________________________________________________________________
GMp12ElasticFF::GMp12ElasticFF(){

}
//______________________________________________________________________________
GMp12ElasticFF::~GMp12ElasticFF(){

}
//______________________________________________________________________________
double GMp12ElasticFF::GetGE(double Q2){
   double GE=0;

   double tau = Q2/(4.*RC::Constants::proton_mass*RC::Constants::proton_mass); 

   const int NA  = 1; 
   double a[NA]  = {-0.21}; 
   // double da[NA] = {0.064}; 
   const int NB  = 3;  
   double b[NB]  = {11.36,18.6,-1.59};    
   // double db[NB] = {0.45 ,2.88, 0.75};   

   double num = 1.; 
   for(int i=0;i<NA;i++) num += a[i]*pow(tau,i+1); 
   double den = 1; 
   for(int i=0;i<NB;i++) den += b[i]*pow(tau,i+1);

   if(den!=0) GE = num/den;
   return GE; 
}
//______________________________________________________________________________
double GMp12ElasticFF::GetGM(double Q2){
   double GM=0;

   double tau = Q2/(4.*RC::Constants::proton_mass*RC::Constants::proton_mass); 

   const int NA  = 1; 
   double a[NA]  = {0.048}; 
   // double da[NA] = {0.042}; 
   const int NB  = 3;  
   double b[NB]  = {10.35,20.6,0.31};    
   // double db[NB] = {0.24,0.31,1.34};   

   double num = 1.; 
   for(int i=0;i<NA;i++) num += a[i]*pow(tau,i+1); 
   double den = 1; 
   for(int i=0;i<NB;i++) den += b[i]*pow(tau,i+1);

   if(den!=0) GM = RC::Constants::mu_p*num/den;
   return GM; 
}
