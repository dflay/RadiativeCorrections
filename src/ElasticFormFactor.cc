#include "ElasticFormFactor.hh"
//______________________________________________________________________________
void ElasticFormFactor::SetTarget(int type){
   fType = type;
   if(fType==RC::Target::kProton){
      fA = 1; 
      fZ = 1; 
   }else if(fType==RC::Target::kNeutron){
      fA = 1; 
      fZ = 0; 
   }else if(fType==RC::Target::kCarbon){
      fA  = 12.011;
      fZ  = 6;
   }else if(fType==RC::Target::kAluminum){
      fA = 26.9815385;
      fZ = 13;
   } else if(fType==RC::Target::kCopper){
      fA = 63.546;
      fZ = 29; 
   }else if(fType==RC::Target::kIron){
      fA = 55.845;
      fZ = 26; 
   }else if(fType==RC::Target::kGold){
      fA = 196.966569;
      fZ = 79;  
   }
   fMT = RC::Constants::proton_mass*fA;
}
