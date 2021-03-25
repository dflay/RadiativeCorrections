#include "AmrounElasticFF.hh"
//______________________________________________________________________________
AmrounElasticFF::AmrounElasticFF(int type){

   // in fm 
   fR[0]  = 0.1; 
   fR[1]  = 0.5;
   fR[2]  = 0.9;
   fR[3]  = 1.3;
   fR[4]  = 1.6;
   fR[5]  = 2.0;
   fR[6]  = 2.4;
   fR[7]  = 2.9;
   fR[8]  = 3.4;
   fR[9]  = 4.0;
   fR[10] = 4.6;
   fR[11] = 5.2; 

   if(type==RC::Target::kHe3){
      // electric 
      fQch[0] = 0.027614;
      fQch[1] = 0.170847;
      fQch[2] = 0.219805;
      fQch[3] = 0.170486;
      fQch[4] = 0.134453;
      fQch[5] = 0.100953;
      fQch[6] = 0.074310;
      fQch[7] = 0.053970;
      fQch[8] = 0.023689;
      fQch[9] = 0.017502;
      fQch[10] = 0.002034;
      fQch[11] = 0.004338; 
      // magnetic 
      fQmag[0] = 0.059785; 
      fQmag[1] = 0.138368;
      fQmag[2] = 0.281326;
      fQmag[3] = 0.000037;
      fQmag[4] = 0.289808;
      fQmag[5] = 0.019056;
      fQmag[6] = 0.114825;
      fQmag[7] = 0.042296;
      fQmag[8] = 0.028345;
      fQmag[9] = 0.018312;
      fQmag[10] = 0.007843;
      fQmag[11] = 0.000000;
   }else if(type==RC::Target::kH3){
      // electric 
      fQch[0] = 0.054706;  
      fQch[1] = 0.172505;  
      fQch[2] = 0.313852;  
      fQch[3] = 0.072056;  
      fQch[4] = 0.225333;  
      fQch[5] = 0.020849;  
      fQch[6] = 0.097374;  
      fQch[7] = 0.022273;  
      fQch[8] = 0.011933;  
      fQch[9] = 0.009121;  
      fQch[10] = 0.000000;  
      fQch[11] = 0.000000;  
      // magnetic 
      fQmag[0]  = 0.075234; 
      fQmag[1]  = 0.164700;
      fQmag[2]  = 0.273033;
      fQmag[3]  = 0.037591;
      fQmag[4]  = 0.252089;
      fQmag[5]  = 0.027036;
      fQmag[6]  = 0.098445;
      fQmag[7]  = 0.040160;
      fQmag[8]  = 0.016696;
      fQmag[9]  = 0.015077;
      fQmag[10] = 0.000000;
      fQmag[11] = 0.000000; 
   }


}
//______________________________________________________________________________
AmrounElasticFF::~AmrounElasticFF(){

}
//______________________________________________________________________________
double AmrounElasticFF::GetGE(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GE   = SumOfGaussians(0,q_fm);
   return GE; 
}
//______________________________________________________________________________
double AmrounElasticFF::GetGM(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GM   = SumOfGaussians(1,q_fm);
   return GM; 
}
//______________________________________________________________________________
double AmrounElasticFF::SumOfGaussians(int type,double q){
   // sum of gaussians (SOG) fit
   // q => sqrt(Q2) in fm^-1 
   double gamma = sqrt(2./3.)*0.8; // in fermi; based on the RMS radius of the proton 
   double T1    = exp(-0.5*q*q*gamma*gamma); 
   double a=0,b=0,sum=0,qr=0,rg=0,FF=0;
   for(int i=0;i<AM_SIZE;i++){
      qr = q*fR[i];
      rg = fR[i]/gamma;
      if(type==0) FF = fQch[i];  
      if(type==1) FF = fQmag[i];  
      a = FF/( 1. + 2.*rg*rg ); 
      b = cos(qr) + 2.*rg*rg*(sin(qr)/qr);  
      sum += a*b; 
   }
   // ensure it's positive
   double res_sq = pow(T1*sum,2.); 
   double res    = sqrt(res_sq);  
   return T1*sum; 
}
