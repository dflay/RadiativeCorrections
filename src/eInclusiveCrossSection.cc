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
double eInclusiveCrossSection::GetMottXS(double Es,double th){
   double thr  = th*deg_to_rad;
   double COS  = cos(thr/2.0);
   double SIN  = sin(thr/2.0);
   double SIN2 = SIN*SIN;
   double T    = HBAR_C*alpha*COS/(2.0*SIN2);
   double Mott = T*T/(Es*Es);
   return Mott;  // in nb
}
