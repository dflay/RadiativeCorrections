#include "RadiativeCorrections.hh"
//________________________________________________________________________
RadiativeCorrections::RadiativeCorrections(){
   Init();
}
//________________________________________________________________________
RadiativeCorrections::~RadiativeCorrections(){

}
//________________________________________________________________________
void RadiativeCorrections::Init(){
   fDeltaE   = 0.01;           // in GeV
   fMT       = 0;
   fZ        = 0;
   fA        = 0;
   fb        = 0;
   fXi       = 0;
   fEta      = 0;
   fTa       = 0;
   fTb       = 0;
   fT        = 0;
   fThDeg    = 0;
   fEs       = 0;
   fEp       = 0;
   fR        = 0;
   fCFACT    = 0;
   fMT       = 0;
   fThreshold = RC::kPion; 
}
//_____________________________________________________________________________________________
void RadiativeCorrections::SetKinematicVariables(double Es,double Ep,double thDeg){
   // set important variables 
   fEs    = Es; 
   fEp    = Ep; 
   fThDeg = thDeg; 
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateVariables(){
   // update variables that depend on Es, Ep, th, Ta, Tb  
   CalculateEta();
   CalculateB();
   CalculateXi();
   CalculateR();
   CalculateCFACT();
}
//_____________________________________________________________________________________________
double RadiativeCorrections::Radiate(){
   // set important variables 
   fZ     = fInclXS->GetZ();
   fA     = fInclXS->GetA();
   fEs    = fInclXS->GetEs();
   fEp    = fInclXS->GetEp();
   fThDeg = fInclXS->GetTh();
   fMT    = fA*proton_mass;      // set the target mass 
   fT     = fTa + fTb;

   if( (fTa==0)||(fTb==0) ){
      std::cout << "[RadiativeCorrections::Radiate]: WARNING! Radiation lengths are zero! " << std::endl;
   }

   CalculateVariables(); 

   double BornXS = fInclXS->GetBornXS();
   double AnsEs  = CalculateEsIntegral();
   double AnsEp  = CalculateEpIntegral();
   double RadXS  = fCFACT*BornXS + AnsEs + AnsEp;

   return RadXS;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetPhi(double v){
   double phi = 1.0 - v + (3.0/4.0)*v*v;
   return phi;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetTr(double Q2){
   // General terms
   double M2 = electron_mass*electron_mass;
   // Individual terms
   double T1 = (1.0/fb)*(alpha/PI);
   double T2 = log(Q2/M2) - 1.0;
   // Put it all together 
   double Tr = T1*T2;
   return Tr;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetFTilde(double Q2){
   // General terms
   double M2     = electron_mass*electron_mass;
   double PI2    = PI*PI;
   double thr    = fThDeg*deg_to_rad;
   double COS    = cos(thr/2.0);
   double COS2   = COS*COS;
   double SPENCE = GetSpence(COS2);
   // Individual terms 
   double T1     = 1.0 + 0.5772*fb*fT;
   double T2     = (2.0*alpha/PI)*( (-14.0/9.0) + (13.0/12.0)*log(Q2/M2) );
   double T3     = (-1.0)*(alpha/(2.0*PI))*log( pow(fEs/fEp,2.0) );
   double T4     = (alpha/PI)*( (PI2/6.0) - SPENCE );
   // Put it all together
   double FTilde = T1+T2+T3+T4;
   return FTilde;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetEsMin(double Ep){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Ep;
      denom = 1.0 - (2.0*Ep/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = pion_mass*pion_mass + 2.*proton_mass*pion_mass + 2.*proton_mass*Ep;
      denom = 2.*proton_mass - (4.0*Ep)*SIN2;
   }

   double EsMin = num/denom;
   return EsMin;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetEpMax(double Es){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Es;
      denom = 1.0 + (2.0*Es/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = 2.*proton_mass*Es - 2.*proton_mass*pion_mass - pion_mass*pion_mass;
      denom = 2.*proton_mass + (4.0*Es)*SIN2;
   }

   double EpMax = num/denom;
   return EpMax;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetSpence(double x){
   // converted from radcor.f: 
   double num=0,denom=0,Index=0;
   double PI2 = PI*PI;
   double ans = (PI2/6.0) - log(x)*log(1.-x);

   for(int i=0;i<50;i++){
      Index   = (double)i + 1.0;
      num     = pow(x,i+1.0);
      denom   = pow(Index,2.);
      if(denom>0){
	 ans -= num/denom;
      }else{
	 ans -= 0;
      }
   }

   return ans;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateEta(){
   double Z23   = pow(fZ,-2.0/3.0);
   double Z13   = pow(fZ,-1.0/3.0);
   double num   = log(1440.0*Z23);
   double denom = log(183.0*Z13);
   fEta         = num/denom;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateB(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = 1.0;
   double T2  = (1.0/9.0)*( (fZ+1.0)/(fZ+fEta) );
   double T3  = 1.0/log(183.0*Z13);
   fb         = (4.0/3.0)*(T1 + T2*T3);
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateXi(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = PI*electron_mass/(2.0*alpha);
   double T2  = fT/( (fZ+fEta)*log(183.0*Z13) );
   fXi        = T1*T2;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateR(){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;
   double num   = fMT + 2.0*fEs*SIN2;
   double denom = fMT - 2.0*fEp*SIN2;
   fR           = num/denom;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateCFACT(){
   // General terms 
   double Q2     = Kinematics::GetQ2(fEs,fEp,fThDeg);
   double Tr     = GetTr(Q2);
   double FTilde = GetFTilde(Q2);
   // First term
   double Term1  = fR*fDeltaE/fEs;
   double Exp1   = fb*(fTb+Tr);
   double T1     = pow(Term1,Exp1);
   // Second term 
   double Term2  = fDeltaE/fEp;
   double Exp2   = fb*(fTa+Tr);
   double T2     = pow(Term2,Exp2);
   // Third term
   double num    = fXi/fDeltaE;
   double denom  = 1.0 - fb*(fTa+fTb+2.0*Tr);
   double T3     = 1.0 - num/denom;
   // Put it all together 
   fCFACT       = FTilde*T1*T2*T3;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::EsIntegrand(const double EsPrime){
   // general terms
   double Q2       = Kinematics::GetQ2(EsPrime,fEp,fThDeg);
   double FTilde   = GetFTilde(Q2);
   double Tr       = GetTr(Q2);
   double dEs      = fEs-EsPrime;
   double v        = dEs/fEs;
   double phi      = GetPhi(v);
   fInclXS->SetEs(EsPrime);
   fInclXS->SetEp(fEp);
   fInclXS->SetTh(fThDeg);
   double Sig      = fInclXS->GetBornXS();
   if(Sig!=Sig){
      std::cout << "[RadiativeCorrections::EsIntegrand]: Invalid cross section! " << std::endl;
      exit(1);
   }
   double SigTilde = FTilde*Sig;
   // first term 
   double Term1    = dEs/(fEp*fR);
   double Exp1     = fb*(fTa+Tr);
   double T1       = pow(Term1,Exp1);
   // second term
   double Term2    = dEs/fEs;
   double Exp2     = fb*(fTb+Tr);
   double T2       = pow(Term2,Exp2);
   // third term 
   double T3       = fb*( ((fTb+Tr)/dEs)*phi + fXi/(2.0*pow(dEs,2.0)) );
   // put it all together
   double FES      = T1*T2*T3*SigTilde;

   return FES;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::EpIntegrand(const double EpPrime){
   // general terms 
   double Q2       = Kinematics::GetQ2(fEs,EpPrime,fThDeg);
   double Tr       = GetTr(Q2);
   double FTilde   = GetFTilde(Q2);
   double dEp      = EpPrime-fEp;
   double v        = dEp/EpPrime;
   double phi      = GetPhi(v);
   fInclXS->SetEs(fEs);
   fInclXS->SetEp(EpPrime);
   fInclXS->SetTh(fThDeg);
   double Sig      = fInclXS->GetBornXS();
   if(Sig!=Sig){
      std::cout << "[RadiativeCorrections::EpIntegrand]: Invalid cross section! " << std::endl;
      exit(1);
   }
   double SigTilde = FTilde*Sig;
   // first term 
   double Term1    = dEp/(EpPrime);
   double Exp1     = fb*(fTa+Tr);
   double T1       = pow(Term1,Exp1);
   // second term
   double Term2    = (dEp*fR)/fEs;
   double Exp2     = fb*(fTb+Tr);
   double T2       = pow(Term2,Exp2);
   // third term 
   double T3       = fb*( ((fTa+Tr)/dEp)*phi + fXi/(2.0*pow(dEp,2.0)) );
   // put it all together
   double FEP      = T1*T2*T3*SigTilde;

   return FEP;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::CalculateEsIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = GetEsMin(fEp);
   double max     = fEs - fR*fDeltaE;
   double AnsEs   = Integrate(&RadiativeCorrections::EsIntegrand,min,max,epsilon,depth);
   return AnsEs;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::CalculateEpIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = fEp + fDeltaE;
   double max     = GetEpMax(fEs);
   double AnsEp   = Integrate(&RadiativeCorrections::EpIntegrand,min,max,epsilon,depth);
   return AnsEp;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::ElasticTail_exact(){
   // Elastic radiative tail using the exact formalism 
   // Phys. Rev. D 12, 1884 (A62)
   CalculateVariables(); 
   double Rt       = 1.;  // FIXME: account for radiation from the target.  Eq A59--A61.  
   double sigma_ex = ElasticTail_sigmaEx();
   double sigma_b  = ElasticTail_sigmaB();  
   double fsoft    = GetF_soft();
   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::ElasticTail_exact]: " << std::endl;
      std::cout << "Es = " << fEs << ", Ep = " << fEp << ", th = " << fThDeg << std::endl;
      std::cout << "sigma_ex = " << sigma_ex << " " 
	        << "sigma_b = "  << sigma_b  << " " 
	        << "F_soft = "   << fsoft << std::endl;
   }
   double el_tail  = fsoft*(sigma_ex*Rt + sigma_b); 
   return el_tail; 
}
//___________________________________________________________________________________
double RadiativeCorrections::ElasticTail_peakApprox(){
   // Elastic radiative tail, peaking approximation 
   // Phys. Rev. D 12, 1884 (A63) 
   CalculateVariables(); 
   double sigma_p = ElasticTail_sigmaP(); 
   double sigma_b = ElasticTail_sigmaB();
   double fsoft   = GetF_soft();  
   double el_tail = fsoft*(sigma_p + sigma_b);
   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::ElasticTail_peakApprox]: " << std::endl;
      std::cout << "Es = " << fEs << ", Ep = " << fEp << ", th = " << fThDeg << std::endl;
      std::cout << "sigma_p = " << sigma_p  << " "
                << "sigma_b = " << sigma_b  << " "
                << "F_soft = "  << fsoft << std::endl;
   }
 
   return el_tail; 
}
//_____________________________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaEx(){
   // Elastic radiative tail using the exact formalism 
   // Phys. Rev. D 12, (A24)
   int depth = 20; 
   double epsilon = 1e-10;
   double min = -1; 
   double max =  1; 
   double Ans = Integrate(&RadiativeCorrections::ElasticTail_sigmaEx_Integrand,min,max,epsilon,depth);
   // scale factor 
   double sf=0;
   double sf_num = pow(alpha,3)*fEp; 
   double sf_den = 2.*PI*fEs;
   if(sf_den!=0) sf = sf_num/sf_den; 
   return sf*Ans;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaEx_Integrand(const double cos_thk){
   // Elastic radiative tail using the exact formalism 
   // Phys. Rev. D 12, (A24--41)
   // 4-vector definitions
   // s = 4-momentum of incident electron (Es,vec(s))  
   // p = 4-momentum of scattered electron (Ep,vec(p))
   // t = 4-momentum of target (MT,0) 
   // k = 4-momentum of real photon emitted (w,vec(k))  
   // u = s + t - p 
   // scattering angle theta
   double thr     = fThDeg*deg_to_rad; 
   double COS     = cos(thr);
   // vector calcs  
   double s_mag   = sqrt(fEs*fEs - electron_mass*electron_mass); // FIXME: is this right?  
   double p_mag   = sqrt(fEp*fEp - electron_mass*electron_mass); // FIXME: is this right?  
   double s_dot_p = fEs*fEp - p_mag*s_mag*COS;
   double u_0     = fEs + fMT - fEp;
   double u_sq    = 2.*electron_mass*electron_mass + fMT*fMT - 2.*s_dot_p + 2.*fMT*(fEs-fEp); 
   double u_mag   = sqrt(u_0*u_0-u_sq);  
   // theta_k (angle between u and k)  
   double sin_thk = sqrt(1. - cos_thk*cos_thk);  
   // theta_p 
   double cos_thp = (s_mag*COS - p_mag)/u_mag;  
   double sin_thp = sqrt(1. - cos_thp*cos_thp);  
   // theta_s 
   double cos_ths = (s_mag - p_mag*COS)/u_mag;  
   // other variables 
   double w       = 0.5*(u_sq - fMT*fMT)/(u_0  - u_mag*cos_thk);  
   double q_sq    = 2.*electron_mass*electron_mass - 2.*s_dot_p - 2.*w*(fEs-fEp) + 2.*w*u_mag*cos_thk;
   double a       = w*(fEp - p_mag*cos_thp*cos_thk);
   double a_pr    = w*(fEs - s_mag*cos_ths*cos_thk);
   double b_pr    = (-1.)*w*p_mag*sin_thp*sin_thk;
   double v       = 1./(a_pr - a);
   double x       = sqrt( a*a - b_pr*b_pr );  
   double y       = sqrt( a_pr*a_pr - b_pr*b_pr ); 
   // form factors (NOTE: -q2 = Q2!)
   double FTilde   = GetFTilde(-q_sq);   // FIXME: If we have the proton, is this set to 1? 
   double tau      = -q_sq/(4.*fMT*fMT); 
   double GE       = fFormFactor->GetGE(-q_sq);   
   double GM       = fFormFactor->GetGM(-q_sq);  
   double W1       = tau*GM*GM;
   double W2       = (GE*GE + tau*GM*GM)/(1+tau);
   double W1_tilde = FTilde*W1; 
   double W2_tilde = FTilde*W2; 
   // integrand terms
   // T0 multiplies everything 
   double T0=0;
   double T0_num = 2.*fMT*w; 
   double T0_den = q_sq*q_sq*(u_0 - u_mag*cos_thk); 
   if(T0_den!=0) T0 = T0_num/T0_den;
   // T1 is scaled by W2_tilde 
   double T1a=0;
   double T1a_num = (-1.)*(a*electron_mass*electron_mass)*(2.*fEs*(fEp + w) + q_sq/2. );  
   double T1a_den = pow(x,3);
   if(T1a_den!=0) T1a = T1a_num/T1a_den; 
   double T1b=0;
   double T1b_num = (-1.)*(a_pr*electron_mass*electron_mass)*(2.*fEp*(fEs - w) + q_sq/2. );  
   double T1b_den = pow(y,3);
   if(T1b_den!=0) T1b = T1b_num/T1b_den;
   double T1c = -2.;
   double T1d_sf=0;                // note this is rewritten to be nicer in code
   double T1d_sf_num = 2.*v*(y-x);   
   double T1d_sf_den = x*y; 
   if(T1d_sf_den!=0) T1d_sf = T1d_sf_num/T1d_sf_den;
   double T1d = T1d_sf*(electron_mass*electron_mass*(s_dot_p - w*w) + s_dot_p*(2.*fEs*fEp - s_dot_p + w*(fEs-fEp)) ); 
   double T1e=0; 
   double T1e_num = 2.*(fEs*fEp + fEs*w + fEp*fEp) + q_sq/2. - s_dot_p - electron_mass*electron_mass;
   double T1e_den = x; 
   if(T1e_den!=0) T1e = T1e_num/T1e_den; 
   double T1f=0; 
   double T1f_num = 2.*(fEs*fEp - fEp*w + fEs*fEs) + q_sq/2. - s_dot_p - electron_mass*electron_mass;
   double T1f_den = y; 
   if(T1f_den!=0) T1f = (-1.)*T1f_num/T1f_den; // NOTE the minus sign!  
   double T1 = W2_tilde*(T1a + T1b + T1c + T1d + T1e + T1f);
   // T2 is scaled by W1_tilde
   double T2a_sf=0; 
   double T2a_sf_num = a*pow(y,3) + a_pr*pow(x,3);
   double T2a_sf_den = pow(x,3)*pow(y,3);
   if(T2a_sf_den!=0) T2a_sf = T2a_sf_num/T2a_sf_den;
   double T2a = T2a_sf*electron_mass*electron_mass*(2.*electron_mass*electron_mass + q_sq);
   double T2b = 4.; 
   double T2c=0;
   double T2c_num = 4.*v*(y-x)*s_dot_p*(s_dot_p - 2.*electron_mass*electron_mass);  
   double T2c_den = x*y;  
   if(T2c_den!=0) T2c = T2c_num/T2c_den;  
   double T2d=0;
   double T2d_num = (y-x)*( 2.*s_dot_p + 2.*electron_mass*electron_mass - q_sq); 
   double T2d_den = x*y; 
   if(T2d_den!=0) T2d = T2d_num/T2d_den;  
   double T2 = W1_tilde*(T2a + T2b + T2c + T2d); 
   double val = T0*(T1 + T2);
   if(fVerbosity>1){
      std::cout << "[RadiativeCorrections::ElasticTail_sigmaEx_Integrand]: " << std::endl;
      std::cout << "Es = " << fEs << " Ep = " << fEp << " th = " << fThDeg << " Q2 = " << -q_sq << std::endl; 
      std::cout << "W1_tilde = " << W1_tilde << " W2_tilde = " << W2_tilde << " " 
	 << "T0 = " << T0 << " T1 = " << T1 << " T2 = " << T2 << " el_tail = " << val << std::endl; 
      std::cout << "----------------------" << std::endl; 
   } 
   return val;
}

//___________________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaB(){
   // real bremsstrahlung and ionization loss  
   // looks like the angle-peaking approximation (c.f., sigmaP below) 
   // Phys. Rev. D 12, 1884 (A49) 
   double thr    = fThDeg*deg_to_rad; 
   double SIN    = sin(thr/2.);
   double SIN2   = SIN*SIN;
   double ws     = GetWs(fEs,fEp,fThDeg);  
   double wp     = GetWp(fEs,fEp,fThDeg);
   double vp     = wp/(fEp+wp);  
   double vs     = ws/fEs;  
   double Q2     = Kinematics::GetQ2(fEs,fEp,fThDeg); 
   // first term  
   double T1=0;
   double T1_num = fMT + 2.*(fEs-ws)*SIN2;
   double T1_den = fMT - 2.*fEp*SIN2; 
   if(T1_den!=0) T1 = T1_num/T1_den; 
   // second term  
   double FTilde = GetFTilde(Q2);
   double T2a=0; 
   double T2a_num = fb*fTb*GetPhi(vs);
   double T2a_den = ws; 
   if(T2a_den!=0) T2a = T2a_num/T2a_den;   
   double T2b=0; 
   double T2b_num = fXi;
   double T2b_den = 2.*ws*ws; 
   if(T2b_den!=0) T2b = T2b_num/T2b_den;   
   double T2_sf = FTilde*sigma_el(fEs-ws);
   double T2    = T2_sf*(T2a + T2b);  
   // third term
   double T3a=0;
   double T3a_num = fb*fTa*GetPhi(vp); 
   double T3a_den = wp; 
   if(T3a_den!=0) T3a = T3a_num/T3a_den;
   double T3b=0; 
   double T3b_num = fXi;
   double T3b_den = 2.*wp*wp; 
   if(T3b_den!=0) T3b = T3b_num/T3b_den;   
   double T3_sf = FTilde*sigma_el(fEs);  
   double T3    = T3_sf*(T3a + T3b);  
   // put it together 
   double val = T1*T2 + T3; 
   return val;  
}
//___________________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaP(){
   // Elastic radiative tail using the angle-peaking approximation
   // Phys. Rev. D. 12, 1884 (A56)  
   double thr    = fThDeg*deg_to_rad; 
   double SIN    = sin(thr/2.);
   double SIN2   = SIN*SIN;
   double ws     = GetWs(fEs,fEp,fThDeg);  
   double wp     = GetWp(fEs,fEp,fThDeg);
   double vp     = wp/(fEp+wp);  
   double vs     = ws/fEs;  
   double Q2     = Kinematics::GetQ2(fEs,fEp,fThDeg); 
   double Tr     = GetTr(Q2); 
   // first term  
   double T1=0;
   double T1_num = fMT + 2.*(fEs-ws)*SIN2;
   double T1_den = fMT - 2.*fEp*SIN2; 
   if(T1_den!=0) T1 = T1_num/T1_den; 
   // second term  
   double FTilde = GetFTilde(Q2);
   double T2_sf  = FTilde*sigma_el(fEs-ws);
   double T2_num = fb*Tr*GetPhi(vs); 
   double T2_den = ws;
   double T2=0; 
   if(T2_den!=0) T2 = T2_sf*T2_num/T2_den;
   // third term 
   double T3_sf  = FTilde*sigma_el(fEs);
   double T3_num = fb*Tr*GetPhi(vp); 
   double T3_den = wp;
   double T3=0; 
   if(T3_den!=0) T3 = T3_sf*T3_num/T3_den;
   // put it together 
   double val = T1*T2 + T3; 
   return val;  
}
//___________________________________________________________________________________
double RadiativeCorrections::sigma_el(double Es){
   // elastic cross section
   // Phys.Rev.D 12,1884 (A13)
   double thr   = fThDeg*deg_to_rad; 
   double COS   = cos(thr/2.); 
   double COS2  = COS*COS; 
   double SIN   = sin(thr/2.); 
   double SIN2  = SIN*SIN;
   double TAN2=0;
   if(COS2!=0) TAN2 = SIN2/COS2; 
   double Ep    = Es/(1 + (2.*Es/fMT)*SIN2);
   double Q2    = Kinematics::GetQ2(Es,Ep,fThDeg); 
   double tau   = Q2/(4.*fMT*fMT); 
   double GE    = fFormFactor->GetGE(Q2);   
   double GM    = fFormFactor->GetGM(Q2);  
   double W1    = tau*GM*GM;
   double W2    = (GE*GE + tau*GM*GM)/(1+tau);
   // term 1 
   double T1=0;
   double T1_num = alpha*alpha*COS2; 
   double T1_den = 4.*Es*Es*SIN2*SIN2; 
   if(T1_den!=0) T1 = T1_num/T1_den;  
   // term 2 
   double T2=0;
   double T2_num = Ep;  
   double T2_den = Es;
   if(T2_den!=0) T2 = T2_num/T2_den; 
   // term 3 
   double T3     = W2 + 2.*TAN2*W1;  
   double xs_el  = T1*T2*T3; 
   return xs_el; 
}
//___________________________________________________________________________________
double RadiativeCorrections::GetF_soft(){
   // Multiple-photon correction
   // Phys.Rev.D 12,1884 (A58)
   double ws    = GetWs(fEs,fEp,fThDeg);  
   double wp    = GetWp(fEs,fEp,fThDeg);  
   double Q2    = Kinematics::GetQ2(fEs,fEp,fThDeg); 
   double Tr    = GetTr(Q2); 
   double arg1  = fb*(fTb+Tr); 
   double arg2  = fb*(fTa+Tr);
   double Fsoft = pow(ws/fEs,arg1)*pow(wp/(fEp+wp),arg2);
   return Fsoft;  
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetWs(double Es,double Ep,double th){
   double thr    = th*deg_to_rad; 
   double SIN    = sin(thr/2.);
   double SIN2   = SIN*SIN; 
   double T1     = Es; 
   double T2_num = Ep;
   double T2_den = 1. - (2.*Ep/fMT)*SIN2;
   double T2=0; 
   if(T2_den!=0) T2 = T2_num/T2_den;  
   double ws = T1 - T2; 
   return ws;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetWp(double Es,double Ep,double th){
   double thr    = th*deg_to_rad; 
   double SIN    = sin(thr/2.);
   double SIN2   = SIN*SIN; 
   double T2     = Ep; 
   double T1_num = Es; 
   double T1_den = 1. + (2.*Es/fMT)*SIN2;
   double T1=0;
   if(T1_den!=0) T1 = T1_num/T1_den;  
   double wp = T1 - T2; 
   return wp;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::Integrate(double (RadiativeCorrections::*f)(const double),double A,double B,double epsilon,int Depth){
   // Adaptive Simpson's Rule
   double C   = (A + B)/2.0;
   double H   = B - A;
   double fa  = (this->*f)(A);
   double fb  = (this->*f)(B);
   double fc  = (this->*f)(C);
   double S   = (H/6.0)*(fa + 4.0*fc + fb);
   double ans = AdaptiveSimpsonAux(f,A,B,epsilon,S,fa,fb,fc,Depth);
   return ans;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::AdaptiveSimpsonAux(double (RadiativeCorrections::*f)(const double),
      double A,double B,double epsilon,
      double S,double fa,double fb,double fc,int bottom){
   // Recursive auxiliary function for AdaptiveSimpson() function
   double C      = (A + B)/2.0;
   double H      = B - A;
   double D      = (A + C)/2.0;
   double E      = (C + B)/2.0;
   double fd     = (this->*f)(D);
   double fe     = (this->*f)(E);
   double Sleft  = (H/12.0)*(fa + 4.0*fd + fc);
   double Sright = (H/12.0)*(fc + 4.0*fe + fb);
   double S2     = Sleft + Sright;
   if( (bottom <= 0) || (fabs(S2 - S) <= 15.0*epsilon) ){
      return S2 + (S2 - S)/15;
   }
   double arg = AdaptiveSimpsonAux(f,A,C,epsilon/2.0,Sleft, fa,fc,fd,bottom-1) +
      AdaptiveSimpsonAux(f,C,B,epsilon/2.0,Sright,fc,fb,fe,bottom-1);
   return arg;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::Print(){
   std::cout << "------------------------------------"              << std::endl;
   std::cout << "Radiative correction quantities: "                 << std::endl;
   std::cout << "DeltaE = " << std::fixed      << std::setprecision(4) << fDeltaE << " [GeV]"  << std::endl;
   std::cout << "Constants for given thicknesses: "                 << std::endl;
   std::cout << "Tb     = " << std::scientific << std::setprecision(4) << fTb     << " [#X0]"  << std::endl;
   std::cout << "Ta     = " << std::scientific << std::setprecision(4) << fTa     << " [#X0]"  << std::endl;
   std::cout << "eta    = " << std::scientific << std::setprecision(4) << fEta    << " [-]"    << std::endl;
   std::cout << "b      = " << std::fixed      << std::setprecision(4) << fb      << " [-]"    << std::endl;
   std::cout << "xi     = " << std::scientific << std::setprecision(4) << fXi     << " [GeV]"  << std::endl;
   std::cout << "Values that change for each (Es,Ep): "   << std::endl;
   std::cout << "Es         = " << std::fixed      << std::setprecision(4) << fEs     << " [GeV]"   << std::endl;
   std::cout << "Ep         = " << std::fixed      << std::setprecision(4) << fEp     << " [GeV]"   << std::endl;
   std::cout << "R          = " << std::fixed      << std::setprecision(4) << fR      << " [-]"     << std::endl;
   std::cout << "CFACT      = " << std::scientific << std::setprecision(4) << fCFACT  << " [-]"     << std::endl;
}
