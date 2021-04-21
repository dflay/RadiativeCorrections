#include "Target.hh"
//______________________________________________________________________________
Target::Target(int type,double A,double Z,double rho){

   fType = type;

   if(fType==RC::Target::kProton){
      fA    = 1;
      fZ    = 1;
      fName = "Proton"; 
   }else if(fType==RC::Target::kNeutron){
      fA    = 1;
      fZ    = 0;
      fName = "Neutron";
   }else if(fType==RC::Target::kDeuteron){
      fA    = 2.01410177811;  
      fZ    = 1;  
      fName = "Deuteron";
   }else if(fType==RC::Target::kHe3){
      fA    = 3.0160293;  
      fZ    = 2; 
      fName = "Helium-3";  
   }else if(fType==RC::Target::kCarbon){
      fA    = 12.011;
      fZ    = 6;
      fName = "Carbon"; 
   }else if(fType==RC::Target::kAluminum){
      fA    = 26.9815385;
      fZ    = 13;
      fName = "Aluminum"; 
   } else if(fType==RC::Target::kCopper){
      fA    = 63.546;
      fZ    = 29;
      fName = "Copper"; 
   }else if(fType==RC::Target::kIron){
      fA    = 55.845;
      fZ    = 26;
      fName = "Iron"; 
   }else if(fType==RC::Target::kGold){
      fA    = 196.966569;
      fZ    = 79;
      fName = "Gold"; 
   }

   fMT = RC::Constants::proton_mass*fA;
   fN  = (int)(fA - fZ);

}
//______________________________________________________________________________
Target::~Target(){

}
//______________________________________________________________________________
void Target::Print(){
   char msg[200]; 
   std::cout << "----------------------------------------" << std::endl; 
   std::cout << "Target Parameters" << std::endl;
   std::cout << "Name: " << fName << std::endl;
   sprintf(msg,"A = %.2lf g/mol\nZ = %.0lf\nN = %0lf\nrho = %.3lf g/cm^3\nM = %.3lf GeV",fA,fZ,fN,fRho,fMT);
   std::cout << msg << std::endl; 
   sprintf(msg,"Tb = %.3lf #X0\nTa = %.3lf #X0",fTb,fTa);
   std::cout << msg << std::endl;
   std::cout << "----------------------------------------" << std::endl; 
}

