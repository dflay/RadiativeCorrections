#ifndef RC_CONSTANTS_HH
#define RC_CONSTANTS_HH

namespace RC {
   namespace Target {
      enum targType_t { kProton, kNeutron, kDeuteron, kH3, kHe3, kCarbon, kAluminum, kCopper, kIron, kGold };
   } //::Target

   namespace Constants { 
      const double HBAR_C = 624.4197;  // in GeV*nb^(1/2)

      const double hc_micrometer = 1.23984193;            // h*c in eV*um 
      const double hc_meter      = hc_micrometer*1E-6;    // h*c in eV*m 

      const double mu_p=2.793; 
      const double mu_n=-1.913;  

      const double sin2thetaW=0.232;

      const double gv_e=-0.5+2.0*sin2thetaW;
      const double gA_e=-0.5;

      // u
      const double gv_u=0.5-4.0/3.0*sin2thetaW;
      const double gA_u=0.5;

      //d
      const double gv_d=-0.5+2.0/3.0*sin2thetaW;
      const double gA_d=-0.5;

      //s
      const double gv_s=-0.5+2.0/3.0*sin2thetaW;
      const double gA_s=-0.5;

      //charge
      const double Q_e=-1;
      const double Q_u=2.0/3.0;
      const double Q_d=-1.0/3.0;
      const double Q_s=-1.0/3.0;


      const double GF=1.166389E-5;       //GeV^-2
      const double alpha=1.0/137.0359895;
      const double PI=3.14159265;

      // mass 
      const double electron_mass=0.511E-3;  // GeV
      const double proton_mass=0.9383;      // GeV
      const double pion_mass=0.140;         // GeV

      const double deg_to_rad = PI/180.0;

      const double kappa_p = 1.7927; 

      // cross section (in cm^2)  
      const double megabarn	= 1E-18;   
      const double kilobarn	= 1E-21; 
      const double barn	        = 1E-24; 
      const double millibarn	= 1E-27; 
      const double microbarn	= 1E-30; 
      const double nanobarn	= 1E-33; 
      const double picobarn	= 1E-36; 
      const double femtobarn	= 1E-39; 
      const double attobarn	= 1E-42; 
      const double zeptobarn	= 1E-45; 
      const double yoctobarn	= 1E-48;

      // conversion factors  
      const double MUB_PER_GEV_SR = 389.44E+0;  
      const double NB_PER_GEV_SR  = 389.44E+3;  
      const double PB_PER_GEV_SR  = 389.44E+6; 

   } //::Constants 
} //::RC 

#endif
