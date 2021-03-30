#include "CarbonElasticFF.hh"
//______________________________________________________________________________
CarbonElasticFF::CarbonElasticFF(int type){
   fR_rms = 2.469; // rms radius in fm  
   // in fm 
   fR[0]  = 0.0; 
   fR[1]  = 0.4;
   fR[2]  = 1.0;
   fR[3]  = 1.3;
   fR[4]  = 1.7;
   fR[5]  = 2.3;
   fR[6]  = 2.7;
   fR[7]  = 3.5;
   fR[8]  = 4.3;
   fR[9]  = 5.4;
   fR[10] = 6.7;

   for(int i=0;i<C_SIZE;i++){
      fQch[i]  = 0;
      fQmag[i] = 0;
   }

   if(type==RC::Target::kCarbon){
      // electric 
      fQch[0]  = 0.016690;
      fQch[1]  = 0.050325;
      fQch[2]  = 0.128621;
      fQch[3]  = 0.180515;
      fQch[4]  = 0.219097;
      fQch[5]  = 0.278416;
      fQch[6]  = 0.058779;
      fQch[7]  = 0.057817;
      fQch[8]  = 0.007739;
      fQch[9]  = 0.002001;
      fQch[10] = 0.000007;
   }else{
      std::cout << "[CarbonElasticFF]: Invalid target type! " << std::endl;
   }
}
//______________________________________________________________________________
CarbonElasticFF::~CarbonElasticFF(){

}
//______________________________________________________________________________
double CarbonElasticFF::GetGE(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GE   = SumOfGaussians(0,q_fm);
   return GE; 
}
//______________________________________________________________________________
double CarbonElasticFF::GetGM(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GM   = SumOfGaussians(1,q_fm);
   return GM; 
}
//______________________________________________________________________________
double CarbonElasticFF::SumOfGaussians(int type,double q){
   // sum of gaussians (SOG) fit
   // q => sqrt(Q2) in fm^-1 
   // double gamma = sqrt(2./3.)*0.8; // in fermi; based on the RMS radius of proton 
   char msg[200];  
   double gamma = sqrt(2./3.)*fR_rms; // in fermi; based on the RMS radius of C atom 
   double T1    = exp(-0.5*q*q*gamma*gamma); 
   double a=0,b=0,sum=0,qr=0,rg=0,FF=0;
   for(int i=0;i<C_SIZE;i++){
      qr = q*fR[i];
      rg = fR[i]/gamma;
      if(type==0) FF = fQch[i];  
      if(type==1) FF = fQmag[i];  
      a = FF/( 1. + 2.*rg*rg ); 
      if(qr==0){
	 b = cos(qr);
      }else{
	 b = cos(qr) + 2.*rg*rg*(sin(qr)/qr); 
      } 
      sum += a*b; 
      // sprintf(msg,"q = %.3E, r = %.3lf, a = %.3E, b = %.3E, sum = %.3E",q,fR[i],a,b,sum); 
      // std::cout << msg << std::endl;
   }
   // ensure its positive
   double res_sq = pow(T1*sum,2.); 
   double res    = sqrt(res_sq);  
   return T1*sum; 
}
