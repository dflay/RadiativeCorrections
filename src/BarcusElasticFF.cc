#include "BarcusElasticFF.hh"
//______________________________________________________________________________
BarcusElasticFF::BarcusElasticFF(int type){
   // in fm 
   fR[0]  = 0.3; 
   fR[1]  = 0.7;
   fR[2]  = 0.9;
   fR[3]  = 1.1;
   fR[4]  = 1.5;
   fR[5]  = 1.6;
   fR[6]  = 2.2;
   fR[7]  = 2.7;
   fR[8]  = 3.3;
   fR[9]  = 4.2;
   fR[10] = 4.3;
   fR[11] = 4.8; 

   if(type==RC::Target::kHe3){
      // electric 
      fQch[0]  = 0.0996392;
      fQch[1]  = 0.214304;
      fQch[2]  = 0.0199385;
      fQch[3]  = 0.195676;
      fQch[4]  = 0.0785533;
      fQch[5]  = 0.167223;
      fQch[6]  = 0.126926;
      fQch[7]  = 0.0549379;
      fQch[8]  = 0.0401401;
      fQch[9]  = 0.0100803;
      fQch[10] = 0.0007217;
      fQch[11] = 4.98962E-12;
      // magnetic 
      fQmag[0]  = 0.159649; 
      fQmag[1]  = 0.0316168;
      fQmag[2]  = 0.277843;
      fQmag[3]  = 0.0364955;
      fQmag[4]  = 0.0329718;
      fQmag[5]  = 0.233469;
      fQmag[6]  = 0.117059;
      fQmag[7]  = 0.0581085;
      fQmag[8]  = 0.0485212;
      fQmag[9]  = 1.77602E-12;
      fQmag[10] = 0.0240927;
      fQmag[11] = 8.94934E-12;
   }else{
      for(int i=0;i<SB_SIZE;i++){
	 fQch[i] = 0;
	 fQmag[i] = 0;
      }
   }

}
//______________________________________________________________________________
BarcusElasticFF::~BarcusElasticFF(){

}
//______________________________________________________________________________
double BarcusElasticFF::GetGE(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GE   = SumOfGaussians(0,q_fm);
   return GE; 
}
//______________________________________________________________________________
double BarcusElasticFF::GetGM(double Q2){
   // Q2 assumed in GeV^2 
   // q => sqrt(Q2) in fm^-1 
   double q_eV = sqrt(Q2)*1E+9;
   double q_fm = (q_eV/RC::Constants::hc_meter)*1E-15;  // [eV][1/(eV*m)][10^-15 m]/[fm] = fm^-1 
   double GM   = SumOfGaussians(1,q_fm);
   return GM; 
}
//______________________________________________________________________________
double BarcusElasticFF::SumOfGaussians(int type,double q){
   // sum of gaussians (SOG) fit
   // q => sqrt(Q2) in fm^-1 
   double gamma = sqrt(2./3.)*0.8; // in fermi; based on the RMS radius of the proton 
   double T1    = exp(-0.5*q*q*gamma*gamma); 
   double a=0,b=0,sum=0,qr=0,rg=0,FF=0;
   for(int i=0;i<SB_SIZE;i++){
      qr = q*fR[i];
      rg = fR[i]/gamma;
      if(type==0) FF = fQch[i];  
      if(type==1) FF = fQmag[i];  
      a = FF/( 1. + 2.*rg*rg ); 
      b = cos(qr) + 2.*rg*rg*(sin(qr)/qr);  
      sum += a*b; 
   }
   // ensure its positive
   double res_sq = pow(T1*sum,2.); 
   double res    = sqrt(res_sq);  
   return T1*sum; 
}
