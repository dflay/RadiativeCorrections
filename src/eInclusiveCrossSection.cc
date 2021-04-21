#include "eInclusiveCrossSection.hh"
//________________________________________________________________________
eInclusiveCrossSection::eInclusiveCrossSection(){
   Init();
}
//________________________________________________________________________
eInclusiveCrossSection::~eInclusiveCrossSection(){

}
//________________________________________________________________________
void eInclusiveCrossSection::Init(){
   fEs = 0;
   fEp = 0;
   fTh = 0;
}
//________________________________________________________________________
void eInclusiveCrossSection::SetTarget(Target *t){
   fTgt = t;
   fZ   = fTgt->GetZ();  
   fA   = fTgt->GetA(); 
}
//________________________________________________________________________
void eInclusiveCrossSection::SetKinematicVariables(double Es,double Ep,double th){
   fEs  = Es; 
   fEp  = Ep; 
   fTh  = th; 
}
//________________________________________________________________________
double eInclusiveCrossSection::GetMottXS(double Es,double th){
   double thr  = th*RC::Constants::deg_to_rad;
   double COS  = cos(thr/2.0);
   double SIN  = sin(thr/2.0);
   double SIN2 = SIN*SIN;
   double T    = RC::Constants::HBAR_C*RC::Constants::alpha*COS/(2.0*SIN2);
   double Mott = T*T/(Es*Es);
   return Mott;  // in nb
}
