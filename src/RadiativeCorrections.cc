#include "RadiativeCorrections.hh"
//______________________________________________________________________________
RadiativeCorrections::RadiativeCorrections(){
   Init();
}
//______________________________________________________________________________
RadiativeCorrections::~RadiativeCorrections(){

}
//______________________________________________________________________________
void RadiativeCorrections::Init(){
   fDeltaE        = 0.01;           // in GeV
   fMT            = 0;
   fZ             = 0;
   fA             = 0;
   fb             = 0;
   fXi            = 0;
   fEta           = 0;
   fTa            = 0;
   fTb            = 0;
   fT             = 0;
   fThDeg         = 0;
   fEs            = 0;
   fEp            = 0;
   fR             = 0;
   fCFACT         = 0;
   fMT            = 0;
   fNumIter       = 4;                         // number of iterations for unfolding born xs 
   fIterThresh    = 1E-3;                      // tolerance threshold for unfolding born xs  
   fThreshold     = RC::kPion;
   fUnit          = RC::kMicrobarnPerGeVPerSr; // mub/GeV/sr 
   fElasticTail   = false;
   fElasticApprox = false;  
}
//______________________________________________________________________________
void RadiativeCorrections::SetKinematicVariables(double Es,double Ep,double thDeg){
   // set important variables 
   // do this directly for the elastic tail 
   fEs    = Es; 
   fEp    = Ep; 
   fThDeg = thDeg;
   // if we have a XS model, set that too 
   if(fInclXS!=NULL){
      fInclXS->SetEs(Es);  
      fInclXS->SetEp(Ep);  
      fInclXS->SetTh(thDeg);  
   }
}
//______________________________________________________________________________
void RadiativeCorrections::CalculateVariables(){
   // update variables that depend on Es, Ep, th, Z, A, MT, Ta, Tb  
   // target details  
   if(fInclXS!=NULL){
      fZ  = fInclXS->GetZ();
      fA  = fInclXS->GetA();
      fMT = fA*RC::Constants::proton_mass; 
   }else{
      std::cout << "[RadiativeCorrections::CalculateVariables]: NO TARGET SET!" << std::endl; 
      std::cout << "A = " << fA << ", Z = " << fZ << ", MT = " << fMT << std::endl;
   } 
   CalculateEta();
   CalculateB();
   CalculateXi();
   CalculateR();
   CalculateCFACT();
}
//______________________________________________________________________________
double RadiativeCorrections::Radiate(){
   // set important variables
   // target details  
   // fZ     = fInclXS->GetZ();
   // fA     = fInclXS->GetA();
   // fMT    = fA*RC::Constants::proton_mass;  
   // kinematics  
   fEs    = fInclXS->GetEs();
   fEp    = fInclXS->GetEp();
   fThDeg = fInclXS->GetTh();
   fT     = fTa + fTb;

   // if( (fTa==0)||(fTb==0) ){
   //    std::cout << "[RadiativeCorrections::Radiate]: WARNING! Radiation lengths are zero! " << std::endl;
   // }

   CalculateVariables(); 

   double BornXS = fInclXS->GetBornXS();
   double AnsEs  = CalculateEsIntegral();
   double AnsEp  = CalculateEpIntegral();
   double RadXS  = fCFACT*BornXS + AnsEs + AnsEp;

   double el_tail=0;
   if(fElasticTail){
      if(fElasticApprox){
	 el_tail = ElasticTail_peakApprox(); 
      }else{
	 el_tail = ElasticTail_exact(); 
      }
      RadXS += el_tail;   
   }

   return RadXS;
}
//______________________________________________________________________________
int RadiativeCorrections::Unfold(double Es,double th,std::vector<double> Ep,std::vector<double> xsr,std::vector<double> &xsb){
   // unfold the BORN cross section from an input spectrum of RADIATED cross section data 
   // input: experimental/radiated cross section specrum
   // - Es  = beam energy [GeV]
   // - th  = scattered e- angle [deg]  
   // - Ep  = vector of scattered e- energies [GeV]
   // - xsr = vector of experimental/radiated cross sections [user-defined units; e.g., mub/GeV/sr, nb/GeV/sr] 
   // output: BORN cross section  
   // - xsb = vector of born cross sections [user-defined units; e.g., mub/GeV/sr, nb/GeV/sr] 

   // FIXME: Need to call the correct cross section in the integrands! (i.e., the unfolded ones since this is iterative...) 
   // - Requires using a RADIATED model => build the radiated spectra first and store to an interpolation object.  
   // - Update the unfolded spectra on each iteration and propagate to interpolation object  

   const int NPTS = Ep.size(); 

   bool stopCalc=false;
   int k=0;  // track the number of iterations 

   double arg=0,num=0,den=0,sum=0,conv=0; 
   double esInt=0,epInt=0,mott=0,xs_red=0,xs_cor=0;

   std::vector<double> v;                  // to determine convergence 
   std::vector<double> xsz;                // corrected cross section as a function of Ep 
   std::vector< std::vector<double> > xsi; // corrected cross sections on each iteration 

   // fill with the input data 
   xsi.push_back(xsr); 

   char msg[200]; 

   for(int i=1;i<=fNumIter;i++){
      // loop over number of iterations
      if(fVerbosity>1) std::cout << "[RadiativeCorrections::Unfold]: Iteration = " << i << std::endl;
      for(int j=0;j<NPTS;j++){
	 // loop over Ep bins 
         SetKinematicVariables(Es,Ep[j],th);
	 CalculateVariables();
	 esInt = CalculateEsIntegral(); // note: EsMin,EsMax are computed inside this function for the (Es,Ep,th) bin 
	 epInt = CalculateEpIntegral(); // note: EpMin,EpMax are computed inside this function for the (Es,Ep,th) bin
         mott  = fInclXS->GetMottXS(Es,th);  
         // construct reduced cross section 
         // xs_red = xsr[j]/mott;
         xs_red = xsi[i-1][j]/mott;
         xs_cor = (xs_red - (esInt+epInt)/mott)/fCFACT;  // corrected cross section 
         // save the corrected cross section 
         xsz.push_back(xs_cor);
         sprintf(msg,"   Ep = %.3lf, xs_red = %.3lf, esInt = %.3lf, epInt = %.3lf, CFACT = %.3lf, xs_cor = %.3lf",
                 Ep[j],xs_red,esInt,epInt,fCFACT,xs_cor); 
         if(fVerbosity>1) std::cout << msg << std::endl; 
      }
      // done with all points; save result of iteration 
      xsi.push_back(xsz); 
      // check for convergence
      if(i>=2){ 
	 for(int j=0;j<NPTS;j++){
	    num  = xsi[i][j]-xsi[i-1][j];
	    den  = xsi[i][j]+xsi[i-1][j];
            arg  = fabs(num)/den;
	    sum += arg;

	 }
         conv = sum/( (double)NPTS );
	 if(conv<=fIterThresh){
	    stopCalc = true;
	 } 
         sum = 0; 
      }
      // set up for next iteration
      xsz.clear(); 
      k++; // iterate the counter of number of iterations performed  
      // check to see if we're done  
      if(stopCalc){
	 sprintf(msg,"[RadiativeCorrections::Unfold]: Reached threshold of %.3lf with delta = %.3lf",fIterThresh,conv); 
	 std::cout << msg << std::endl;
	 break;
      } 
   }

   std::cout << "[RadiativeCorrections::Unfold]: Calculation complete! Number of iterations: " << k << std::endl;
   // save last iteration to the output vector xsb 
   for(int j=0;j<NPTS;j++){
      xsb.push_back(xsi[k-1][j]); 
   }
 
   return 0;  
}
//______________________________________________________________________________
double RadiativeCorrections::GetPhi(double v){
   double phi = 1.0 - v + (3.0/4.0)*v*v;
   return phi;
}
//______________________________________________________________________________
double RadiativeCorrections::GetTr(double Q2){
   // General terms
   double me = RC::Constants::electron_mass;
   double M2 = me*me;
   // Individual terms
   double T1 = (1.0/fb)*(RC::Constants::alpha/RC::Constants::PI);
   double T2 = log(Q2/M2) - 1.0;
   // Put it all together 
   double Tr = T1*T2;
   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::GetTr]: Q2 = " 
	        << Q2 << ", b = " << fb << ", T1 = " << T1 << ", T2 = " << T2 << std::endl;
   }
   return Tr;
}
//______________________________________________________________________________
double RadiativeCorrections::GetFTilde(double q2){
   // Phys. Rev. D 12, 1884 (1975), eq A44 
   // NOTE: q2 != Q2 here! Q2 depends on fEs, fEp, fThDeg.  q2 is computed as needed.  
   // General terms
   double alpha  = RC::Constants::alpha; 
   double M2     = RC::Constants::electron_mass*RC::Constants::electron_mass;
   double PI     = RC::Constants::PI; 
   double PI2    = PI*PI;
   double thr    = fThDeg*RC::Constants::deg_to_rad;
   double COS    = cos(thr/2.0);
   double COS2   = COS*COS;
   double SPENCE = GetSpence(COS2);
   // Individual terms 
   double T1     = 1.0 + 0.5772*fb*fT;
   double T2     = (2.0*alpha/PI)*( (-14.0/9.0) + (13.0/12.0)*log(q2/M2) );
   double T3     = (-1.0)*(alpha/(2.0*PI))*log(fEs/fEp)*log(fEs/fEp);
   double T4     = (alpha/PI)*( (PI2/6.0) - SPENCE );
   // Put it all together
   double FTilde = T1+T2+T3+T4;
   return FTilde;
}
//______________________________________________________________________________
double RadiativeCorrections::GetEsMin(double Ep){
   double thr   = fThDeg*RC::Constants::deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;
   double Mp    = RC::Constants::proton_mass;
   double Mpi   = RC::Constants::pion_mass;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Ep;
      denom = 1.0 - (2.0*Ep/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = Mpi*Mpi + 2.*Mp*Mpi + 2.*Mp*Ep;
      denom = 2.*Mp - (4.0*Ep)*SIN2;
   }

   double EsMin = num/denom;
   return EsMin;
}
//______________________________________________________________________________
double RadiativeCorrections::GetEpMax(double Es){
   double thr   = fThDeg*RC::Constants::deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;
   double Mp    = RC::Constants::proton_mass;
   double Mpi   = RC::Constants::pion_mass;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Es;
      denom = 1.0 + (2.0*Es/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = 2.*Mp*Es - 2.*Mp*Mpi - Mpi*Mpi;
      denom = 2.*Mp + (4.0*Es)*SIN2;
   }

   double EpMax = num/denom;
   return EpMax;
}
//______________________________________________________________________________
double RadiativeCorrections::GetSpence(double x){
   // converted from radcor.f: 
   double num=0,denom=0,Index=0;
   double PI  = RC::Constants::PI; 
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
//______________________________________________________________________________
void RadiativeCorrections::CalculateEta(){
   double Z23   = pow(fZ,-2.0/3.0);
   double Z13   = pow(fZ,-1.0/3.0);
   double num   = log(1440.0*Z23);
   double denom = log(183.0*Z13);
   fEta         = num/denom;
}
//______________________________________________________________________________
void RadiativeCorrections::CalculateB(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = 1.0;
   double T2  = (1.0/9.0)*( (fZ+1.0)/(fZ+fEta) );
   double T3  = 1.0/log(183.0*Z13);
   fb         = (4.0/3.0)*(T1 + T2*T3);
}
//______________________________________________________________________________
void RadiativeCorrections::CalculateXi(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = RC::Constants::PI*RC::Constants::electron_mass/(2.0*RC::Constants::alpha);
   double T2  = fT/( (fZ+fEta)*log(183.0*Z13) );
   fXi        = T1*T2;
}
//______________________________________________________________________________
void RadiativeCorrections::CalculateR(){
   double thr   = fThDeg*RC::Constants::deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;
   double num   = fMT + 2.0*fEs*SIN2;
   double denom = fMT - 2.0*fEp*SIN2;
   fR           = num/denom;
}
//______________________________________________________________________________
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
//______________________________________________________________________________
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
      return 0;
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
//______________________________________________________________________________
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
      return 0;
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
//______________________________________________________________________________
double RadiativeCorrections::CalculateEsIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = GetEsMin(fEp);
   double max     = fEs - fR*fDeltaE;
   double AnsEs   = Integrate(&RadiativeCorrections::EsIntegrand,min,max,epsilon,depth);
   return AnsEs;
}
//______________________________________________________________________________
double RadiativeCorrections::CalculateEpIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = fEp + fDeltaE;
   double max     = GetEpMax(fEs);
   double AnsEp   = Integrate(&RadiativeCorrections::EpIntegrand,min,max,epsilon,depth);
   return AnsEp;
}
//______________________________________________________________________________
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

   double units = 1; 
   if(fUnit==RC::kMicrobarnPerGeVPerSr) units = RC::Constants::MUB_PER_GEV_SR;
   if(fUnit==RC::kNanobarnPerGeVPerSr ) units = RC::Constants::NB_PER_GEV_SR;
   if(fUnit==RC::kPicobarnPerGeVPerSr ) units = RC::Constants::PB_PER_GEV_SR;

   double el_tail  = units*fsoft*(sigma_ex*Rt + sigma_b); 
   return el_tail; 
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticTail_peakApprox(){
   // Elastic radiative tail, peaking approximation 
   // Phys. Rev. D 12, 1884 (A63) 
   CalculateVariables(); 
   double sigma_p = ElasticTail_sigmaP(); 
   double sigma_b = ElasticTail_sigmaB();
   double fsoft   = GetF_soft();  

   double units = 1; 
   if(fUnit==RC::kMicrobarnPerGeVPerSr) units = RC::Constants::MUB_PER_GEV_SR;
   if(fUnit==RC::kNanobarnPerGeVPerSr ) units = RC::Constants::NB_PER_GEV_SR;
   if(fUnit==RC::kPicobarnPerGeVPerSr ) units = RC::Constants::PB_PER_GEV_SR;

   double el_tail = units*fsoft*(sigma_p + sigma_b);

   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::ElasticTail_peakApprox]: " << std::endl;
      std::cout << "Es = " << fEs << ", Ep = " << fEp << ", th = " << fThDeg << std::endl;
      std::cout << "sigma_p = " << sigma_p  << " "
                << "sigma_b = " << sigma_b  << " "
                << "F_soft = "  << fsoft << std::endl;
   }
 
   return el_tail; 
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaEx(){
   // Elastic radiative tail using the exact formalism 
   // Phys. Rev. D 12, 1884 (A24)
   int depth = 20; 
   double epsilon = 1e-10;
   double min = -1; 
   double max =  1; 
   double Ans = Integrate(&RadiativeCorrections::ElasticTail_sigmaEx_Integrand,min,max,epsilon,depth);
   // scale factor 
   double sf=0;
   double sf_num = pow(RC::Constants::alpha,3)*fEp; 
   double sf_den = 2.*RC::Constants::PI*fEs;
   if(sf_den!=0) sf = sf_num/sf_den; 
   return sf*Ans;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaEx_Integrand(const double cos_thk){
   // Elastic radiative tail using the exact formalism 
   // Phys. Rev. D 12, 1884 (A24--41)
   // 4-vector definitions
   // s = 4-momentum of incident electron (Es,vec(s))  
   // p = 4-momentum of scattered electron (Ep,vec(p))
   // t = 4-momentum of target (MT,0) 
   // k = 4-momentum of real photon emitted (w,vec(k))  
   // u = s + t - p 
   double me      = RC::Constants::electron_mass; 
   // scattering angle theta
   double thr     = fThDeg*RC::Constants::deg_to_rad; 
   double COS     = cos(thr);
   // vector calcs  
   double s_mag   = sqrt(fEs*fEs - me*me); // FIXME: is this right?  
   double p_mag   = sqrt(fEp*fEp - me*me); // FIXME: is this right?  
   double s_dot_p = fEs*fEp - p_mag*s_mag*COS;
   double u_0     = fEs + fMT - fEp;
   double u_sq    = 2.*me*me + fMT*fMT - 2.*s_dot_p + 2.*fMT*(fEs-fEp); 
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
   double q_sq    = 2.*me*me - 2.*s_dot_p - 2.*w*(fEs-fEp) + 2.*w*u_mag*cos_thk;
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
   double T1a_num = (-1.)*(a*me*me)*(2.*fEs*(fEp + w) + q_sq/2. );  
   double T1a_den = pow(x,3);
   if(T1a_den!=0) T1a = T1a_num/T1a_den; 
   double T1b=0;
   double T1b_num = (-1.)*(a_pr*me*me)*(2.*fEp*(fEs - w) + q_sq/2. );  
   double T1b_den = pow(y,3);
   if(T1b_den!=0) T1b = T1b_num/T1b_den;
   double T1c = -2.;
   double T1d_sf=0;                // note this is rewritten to be nicer in code
   double T1d_sf_num = 2.*v*(y-x);   
   double T1d_sf_den = x*y; 
   if(T1d_sf_den!=0) T1d_sf = T1d_sf_num/T1d_sf_den;
   double T1d = T1d_sf*(me*me*(s_dot_p - w*w) + s_dot_p*(2.*fEs*fEp - s_dot_p + w*(fEs-fEp)) ); 
   double T1e=0; 
   double T1e_num = 2.*(fEs*fEp + fEs*w + fEp*fEp) + q_sq/2. - s_dot_p - me*me;
   double T1e_den = x; 
   if(T1e_den!=0) T1e = T1e_num/T1e_den; 
   double T1f=0; 
   double T1f_num = 2.*(fEs*fEp - fEp*w + fEs*fEs) + q_sq/2. - s_dot_p - me*me;
   double T1f_den = y; 
   if(T1f_den!=0) T1f = (-1.)*T1f_num/T1f_den; // NOTE the minus sign!  
   double T1 = W2_tilde*(T1a + T1b + T1c + T1d + T1e + T1f);
   // T2 is scaled by W1_tilde
   double T2a_sf=0; 
   double T2a_sf_num = a*pow(y,3) + a_pr*pow(x,3);
   double T2a_sf_den = pow(x,3)*pow(y,3);
   if(T2a_sf_den!=0) T2a_sf = T2a_sf_num/T2a_sf_den;
   double T2a = T2a_sf*me*me*(2.*me*me + q_sq);
   double T2b = 4.; 
   double T2c=0;
   double T2c_num = 4.*v*(y-x)*s_dot_p*(s_dot_p - 2.*me*me);  
   double T2c_den = x*y;  
   if(T2c_den!=0) T2c = T2c_num/T2c_den;  
   double T2d=0;
   double T2d_num = (y-x)*( 2.*s_dot_p + 2.*me*me - q_sq); 
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
//______________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaB(){
   // real bremsstrahlung and ionization loss  
   // looks like the angle-peaking approximation (c.f., sigmaP below) 
   // Phys. Rev. D 12, 1884 (A49) 
   double thr    = fThDeg*RC::Constants::deg_to_rad; 
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
   double T2a=0; 
   double T2a_num = fb*fTb*GetPhi(vs);
   double T2a_den = ws; 
   if(T2a_den!=0) T2a = T2a_num/T2a_den;   
   double T2b=0; 
   double T2b_num = fXi;
   double T2b_den = 2.*ws*ws; 
   if(T2b_den!=0) T2b = T2b_num/T2b_den;  
   double T2_sf = sigma_el_tilde(fEs-ws);
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
   double T3_sf = sigma_el_tilde(fEs);  
   double T3    = T3_sf*(T3a + T3b);  
   // put it together 
   double val = T1*T2 + T3; 
   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::ElasticTail_sigmaB]: Es = " 
                << fEs << ", Ep = " << fEp 
                << ", T1 = " << T1 << ", T2 = " << T2 << " , T3 = " << T3 << std::endl;
   } 
   return val;  
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticTail_sigmaP(){
   // Elastic radiative tail using the angle-peaking approximation
   // Phys. Rev. D. 12, 1884 (A56)  
   double thr    = fThDeg*RC::Constants::deg_to_rad; 
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
   double T2_sf  = sigma_el_tilde(fEs-ws);
   double T2_num = fb*Tr*GetPhi(vs); 
   double T2_den = ws;
   double T2=0; 
   if(T2_den!=0) T2 = T2_sf*T2_num/T2_den;
   // third term 
   double T3_sf  = sigma_el_tilde(fEs);
   double T3_num = fb*Tr*GetPhi(vp); 
   double T3_den = wp;
   double T3=0; 
   if(T3_den!=0) T3 = T3_sf*T3_num/T3_den;
   // put it together 
   double val = T1*T2 + T3;
   if(fVerbosity>0){
      std::cout << "[RadiativeCorrections::ElasticTail_sigmaP]: Es = " 
                << fEs << ", Ep = " << fEp << ", Q2 = " << Q2 << ", ws = " << ws << ", wp = " << wp << ", Tr = " << Tr
                << ", T1 = " << T1 << ", T2 = " << T2 << " , T3 = " << T3 << std::endl;
   } 
   return val;  
}
//______________________________________________________________________________
double RadiativeCorrections::sigma_el_tilde(double Es){
   // Phys. Rev. D. 12, 1884 (A55) 
   double Ep        = Kinematics::GetEp_Elastic(Es,fThDeg,fMT); 
   double Q2        = Kinematics::GetQ2(Es,Ep,fThDeg);
   double FTilde    = GetFTilde(Q2); 
   double sigmaEl   = sigma_el(Es);
   double sig_tilde = FTilde*sigmaEl;
   if(fVerbosity>0){
      std::cout << "RadiativeCorrections::sigma_el_tilde]: Es = " 
                << Es << ", Ep = " << Ep << ", F_tilde = " << FTilde << ", sig_el = " << sigmaEl << std::endl; 
   }
   return sig_tilde;  
}
//______________________________________________________________________________
double RadiativeCorrections::sigma_el(double Es){
   // elastic cross section
   // Phys.Rev.D 12,1884 (A13)
   double thr   = fThDeg*RC::Constants::deg_to_rad; 
   double COS   = cos(thr/2.); 
   double COS2  = COS*COS; 
   double SIN   = sin(thr/2.); 
   double SIN2  = SIN*SIN;
   double TAN2=0;
   if(COS2!=0) TAN2 = SIN2/COS2; 
   double Ep    = Kinematics::GetEp_Elastic(Es,fThDeg,fMT); 
   double Q2    = Kinematics::GetQ2(Es,Ep,fThDeg); 
   double tau   = Q2/(4.*fMT*fMT); 
   double GE    = fFormFactor->GetGE(Q2);   
   double GM    = fFormFactor->GetGM(Q2);  
   double W1    = tau*GM*GM;
   double W2    = (GE*GE + tau*GM*GM)/(1+tau);
   // term 1 
   double T1=0;
   double alpha  = RC::Constants::alpha; 
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
//______________________________________________________________________________
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
//______________________________________________________________________________
double RadiativeCorrections::GetWs(double Es,double Ep,double th){
   double thr    = th*RC::Constants::deg_to_rad; 
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
//______________________________________________________________________________
double RadiativeCorrections::GetWp(double Es,double Ep,double th){
   double thr    = th*RC::Constants::deg_to_rad; 
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
//______________________________________________________________________________
double RadiativeCorrections::GetRho(double Es,double th){
   // for elastic radiative corrections
   double Ep     = Kinematics::GetEp_Elastic(Es,th,fMT); 
   double Q2     = Kinematics::GetQ2(Es,Ep,th);
   double omega  = Q2/(2.*fMT);  
   double rho_sq = Q2 + 4.*fMT*fMT + 4.*fMT*omega; 
   return sqrt(rho_sq);  
}
//______________________________________________________________________________
double RadiativeCorrections::GetX(double Es,double th){
   // for elastic radiative corrections
   // Note: NOT x Bjorken! 
   double Ep  = Kinematics::GetEp_Elastic(Es,th,fMT); 
   double Q2  = Kinematics::GetQ2(Es,Ep,th);
   double q   = sqrt(Q2); 
   double rho = GetRho(Es,th); 
   double num = pow(rho+q,2); 
   double den = 4.*fMT*fMT; 
   double x   = num/den; 
   return x;  
}
//______________________________________________________________________________
double RadiativeCorrections::GetEta(double Es,double th){
   // lab system recoil factor
   double thr = th*RC::Constants::deg_to_rad; 
   double COS = cos(thr); 
   double eta = 1. + (Es/fMT)*(1.-COS); 
   return eta; 
}
//______________________________________________________________________________
double RadiativeCorrections::GetEta_MTS(double Es,double Ep){
   // from Mo & Tsai, Rev Mod Phys 41 205 (1969), pg 208  
   return Es/Ep;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Delta_MTS(){
   // radiative correction to elastic peak
   // from Mo & Tsai, Rev Mod Phys 41 205 (1969), eq II.6 
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z0_MTS(){
   // radiative correction to elastic peak
   // from Mo & Tsai, Rev Mod Phys 41 205 (1969), eq II.6 
   // Z^(0) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z1_MTS(){
   // radiative correction to elastic peak
   // from Mo & Tsai, Rev Mod Phys 41 205 (1969), eq II.6
   // Z^(1) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z2_MTS(){
   // radiative correction to elastic peak
   // from Mo & Tsai, Rev Mod Phys 41 205 (1969), eq II.6
   // Z^(2) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Delta_MY(){
   // radiative correction to elastic peak
   // from Meister and Yennie 
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z0_MY(){
   // radiative correction to elastic peak
   // from Meister and Yennie 
   // Z^(0) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z1_MY(){
   // radiative correction to elastic peak
   // from Meister and Yennie 
   // Z^(1) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z2_MY(){
   // radiative correction to elastic peak
   // from Meister and Yennie 
   // Z^(2) term  
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Delta_MTJ_E1THDE(double Es,double th,double deltaE){
   // radiative correction to elastic peak
   // from Maximon and Tjon, Phys. Rev. C 62, 054320 (2000), eq 5.2
   // inputs: 
   // - Es     = incident electron energy 
   // - th     = scattered electron angle
   // - deltaE = electron detector acceptance in lab frame 
   double T_z0  = ElasticPeak_Z0_MTJ(Es,th,deltaE);  
   double T_z1  = ElasticPeak_Z1_MTJ(Es,th,deltaE);  
   double T_z2  = ElasticPeak_Z2_MTJ(Es,th,deltaE);  
   double T_del = ElasticPeak_DeltaEl_MTJ(); 
   double res   = T_z0 + fZ*T_z1 + (fZ*fZ)*T_z2 + T_del; 
   return res;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z0_MTJ(double Es,double th,double deltaE){
   // radiative correction to elastic peak
   // from Maximon and Tjon, Phys. Rev. C 62, 054320 (2000), eq 5.2
   // Z^(0) term  
   // inputs: 
   // - Es     = incident electron energy 
   // - th     = scattered electron angle
   // - deltaE = electron detector acceptance in lab frame 
   // general terms
   double E1     = Es; 
   double E3     = Kinematics::GetEp_Elastic(E1,th,fMT);
   double DE     = deltaE;    
   double thr    = th*RC::Constants::deg_to_rad;
   double COS    = cos(thr/2.); 
   double COS2   = COS*COS;
   double SPENCE = GetSpence(COS2);  
   double eta    = GetEta(E1,th);
   double Q2     = Kinematics::GetQ2(E1,E3,th);  
   double m2     = RC::Constants::electron_mass*RC::Constants::electron_mass; 
   double M2     = fMT*fMT;
   // construct Z0 term
   double PI    = RC::Constants::PI;   
   double alpha = RC::Constants::alpha; 
   double sf    = alpha/PI; 
   double T1    = (13./6.)*log(Q2/m2) - (28./9.) - (log(Q2/m2) - 1)*log( 4.*E1*E3/pow(2*eta*DE,2.) ) 
                - 0.5*log(eta)*log(eta) + SPENCE - PI*PI/6.;
   double res   = sf*T1;
   return res;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z1_MTJ(double Es,double th,double deltaE){
   // radiative correction to elastic peak
   // from Maximon and Tjon, Phys. Rev. C 62, 054320 (2000), eq 5.2 
   // Z^(1) term 
   // inputs: 
   // - Es     = incident electron energy 
   // - th     = scattered electron angle
   // - deltaE = electron detector acceptance in lab frame 
   // general terms
   double E1      = Es; 
   double E3      = Kinematics::GetEp_Elastic(E1,th,fMT);
   double DE      = deltaE;    
   double eta     = GetEta(E1,th);
   double x       = GetX(E1,th);
   double arg1    = 1. - eta/x; 
   double SPENCE1 = GetSpence(arg1);  
   double arg2    = 1. - 1./(eta*x); 
   double SPENCE2 = GetSpence(arg2);
   double Q2      = Kinematics::GetQ2(E1,E3,th);  
   double m2      = RC::Constants::electron_mass*RC::Constants::electron_mass; 
   double M2      = fMT*fMT;
   // construct Z1 term  
   double PI      = RC::Constants::PI;    
   double alpha   = RC::Constants::alpha;  
   double sf      = 2.*alpha/PI;
   double T1      = (-1.)*log(eta)*log( Q2*x/pow(2.*eta*DE,2.) ) + SPENCE1 - SPENCE2; 
   double res     = sf*T1;  
   return res;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_Z2_MTJ(double Es,double th,double deltaE){
   // radiative correction to elastic peak
   // from Maximon and Tjon, Phys. Rev. C 62, 054320 (2000), eq 5.2 
   // Z^(2) term  
   // inputs: 
   // - Es     = incident electron energy 
   // - th     = scattered electron angle
   // - deltaE = electron detector acceptance in lab frame 
   // general terms
   double E1      = Es; 
   double E3      = Kinematics::GetEp_Elastic(E1,th,fMT);
   double DE      = deltaE;   
   double Q2      = Kinematics::GetQ2(E1,E3,th);  
   double M2      = fMT*fMT;
   double omega   = Q2/(2.*fMT); 
   double E4      = fMT + omega; 
   double p4      = sqrt(E4*E4 - M2); 
   // double beta4   = sqrt(1. - M2/(E4*E4));  
   double eta     = GetEta(E1,th);
   double x       = GetX(E1,th);
   double rho     = GetRho(E1,th); 
   double arg1    = 1. - 1./(x*x); 
   double SPENCE1 = GetSpence(arg1);  
   double arg2    = -1./x; 
   double SPENCE2 = GetSpence(arg2);
   // construct Z2 term 
   double PI     = RC::Constants::PI;     
   double alpha  = RC::Constants::alpha;     
   double sf     = alpha/PI;
   double T1     = (E4/p4)*( -0.5*log(x)*log(x) - log(x)*log(rho*rho/M2) + log(x) ); 
   double T2     = (-1.)*( (E4/p4)*log(x) - 1.)*log( M2/pow(2.*eta*DE,2.) );
   double T3     = 1.;
   double T4     = (E4/p4)*( (-1.)*SPENCE1 + 2.*SPENCE2+ PI*PI/6. ); 
   double res    = sf*(T1+T2+T3+T4);  
   return res;
}
//______________________________________________________________________________
double RadiativeCorrections::ElasticPeak_DeltaEl_MTJ(){
   // radiative correction to elastic peak
   // from Maximon and Tjon, Phys. Rev. C 62, 054320 (2000), eq 3.37 
   // delta_el term, accounts for nucleon size effects 
   return 0;
}
//______________________________________________________________________________
double RadiativeCorrections::Integrate(double (RadiativeCorrections::*f)(const double),
      double A,double B,double epsilon,int Depth){
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
//______________________________________________________________________________
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
//______________________________________________________________________________
void RadiativeCorrections::Print(){
   std::cout << "------------------------------------"              << std::endl;
   std::cout << "Settings: " << std::endl;
   std::cout << "Elastic tail enabled:    " << fElasticTail    << std::endl; 
   std::cout << "Elastic tail exact form: " << !fElasticApprox << std::endl;    
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
