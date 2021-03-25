#include "eElasticCrossSection.hh"
//______________________________________________________________________________
eElasticCrossSection::eElasticCrossSection(double Z,double A){
   fZ  = Z;
   fA  = A;
}
//______________________________________________________________________________
eElasticCrossSection::~eElasticCrossSection(){

}
//______________________________________________________________________________
double eElasticCrossSection::GetBornXS(){
   // Es = incident electron energy
   // th = scattered electron angle
   double MT    = fA*RC::Constants::proton_mass; 
   double thr   = fTh*RC::Constants::deg_to_rad;
   double TAN   = tan(thr/2.); 
   double TAN2  = TAN*TAN;
   // kinematics and form factors 
   fEp = Kinematics::GetEp_Elastic(fEs,fTh,MT);  // scattered electron energy 
   double Q2    = Kinematics::GetQ2(fEs,fEp,fTh);
   double tau   = Q2/(4.*MT*MT);
   double GE    = fFormFactor->GetGE(Q2);
   double GM    = fFormFactor->GetGM(Q2);
   double W1    = tau*GM*GM;
   double W2    = (GE*GE + tau*GM*GM)/(1+tau);
   // Mott cross section
   double xs_mott = GetMottXS(fEs,fTh); // in nb  
   // term 2 
   double T1=0;
   double T1_num = fEp;
   double T1_den = fEs;
   if(T1_den!=0) T1 = T1_num/T1_den;
   // term 3 
   double T2     = W2 + 2.*TAN2*W1;
   double xs_el  = xs_mott*T1*T2;
   return xs_el;
}
