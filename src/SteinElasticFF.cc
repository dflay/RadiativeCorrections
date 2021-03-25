#include "SteinElasticFF.hh"
//______________________________________________________________________________
SteinElasticFF::SteinElasticFF(int type){
   SetTarget(type);
}
//______________________________________________________________________________
SteinElasticFF::~SteinElasticFF(){

}
//______________________________________________________________________________
double SteinElasticFF::GetGE(double Q2){
   // Q2 assumed in GeV^2 
   double tau   = Q2/(4.*fMT*fMT);        // unitless  
   // q => sqrt(Q2) in fm^-1 
   double q_eV  = sqrt(Q2)*1E+9;
   double q_fm  = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1
   double q2_fm = pow(q_fm,2.); 
   double b     = 2.4;                       // [fm] 
   double c     = 1.07*pow(fA,1./3.);        // [fm]
   // compute F(q) 
   double num   = exp(-(1./6.)*q2_fm*b*b);
   double den   = 1. + (1./6.)*q2_fm*c*c; 
   double F     = num/den;
   // compute W_2  
   double W_2   = fZ*fZ*F*F;
   // compute G_E  
   // W_2 = (GE^2 + tau*GM^2)/(1+tau) => GE = sqrt( (1+tau)W_2 ) since GM = 0 
   double GE = sqrt( (1.+tau)*W_2 ); 
   return GE; 
}
//______________________________________________________________________________
double SteinElasticFF::GetGM(double Q2){
   // W1 = tau*GM^2 = 0 => GM = 0
   return 0; 
}
