#include "Target.hh"
//______________________________________________________________________________
Target::Target(int type,double A,double Z,double rho){

   fType = type;

   if(fType==RC::Target::kProton){
      fA = 1;
      fZ = 1;
      fName = "Proton"; 
   }else if(fType==RC::Target::kNeutron){
      fA = 1;
      fZ = 0;
      fName = "Neutron"; 
   }else if(fType==RC::Target::kCarbon){
      fA  = 12.011;
      fZ  = 6;
      fName = "Carbon"; 
   }else if(fType==RC::Target::kAluminum){
      fA = 26.9815385;
      fZ = 13;
      fName = "Aluminum"; 
   } else if(fType==RC::Target::kCopper){
      fA = 63.546;
      fZ = 29;
      fName = "Copper"; 
   }else if(fType==RC::Target::kIron){
      fA = 55.845;
      fZ = 26;
      fName = "Iron"; 
   }else if(fType==RC::Target::kGold){
      fA = 196.966569;
      fZ = 79;
      fName = "Gold"; 
   }
   fMT = RC::Constants::proton_mass*fA;

}
//______________________________________________________________________________
Target::~Target(){

}
//______________________________________________________________________________
void Target::Print(){
   
   std::cout << "----------------------------------------" << std::endl; 
   std::cout << "Target Parameters" << std::endl;
   std::cout << "Name: " << fName << std::endl;

   char msg[200]; 
   sprintf(msg,"A = %.2lf g/mol, Z = %.2lf, rho = %.3lf g/cm^3, M = %.3lf GeV",fA,fZ,fRho,fMT);
   std::cout << msg << std::endl; 
   sprintf(msg,"Tb = %.3lf #X0, %.3lf #X0",fTb,fTa);
   std::cout << msg << std::endl;
   std::cout << "----------------------------------------" << std::endl; 
    
}

