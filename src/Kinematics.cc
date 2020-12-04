#include "../include/Kinematics.hh"
//______________________________________________________________________________
namespace Kinematics {
   //________________________________________________________________________
   double GetQ2(double Es,double Ep,double th){
      double thr  = th*deg_to_rad;
      double SIN  = sin(thr/2.);
      double SIN2 = SIN*SIN;
      double Q2   = 4.*Es*Ep*SIN2;
      return Q2;
   }
   //________________________________________________________________________
   double GetW(double Es,double Ep,double th,double M){
      double Nu = Es-Ep;
      double Q2 = GetQ2(Es,Ep,th);
      double W2 = M*M + 2.*M*Nu - Q2;
      return sqrt(W2);
   }
   //________________________________________________________________________
   double GetXbj(double Es,double Ep,double th,double M){
      double Q2    = GetQ2(Es,Ep,th); 
      double Nu    = Es - Ep;
      double num   = Q2;
      double denom = 2.0*M*Nu;
      double x = -1; 
      if(denom!=0){
	 x = num/denom;
      }else{
	 std::cout << "[Kinematics::GetXbj]: ERROR! denominator = 0!" << std::endl;
      }
      return x;
   }
   //________________________________________________________________________
   double GetEpsilon(double Es,double Ep,double th){
      double Nu    = Es - Ep; 
      double Q2    = GetQ2(Es,Ep,th); 
      double thr   = th*deg_to_rad;
      double TAN   = tan(thr/2.0);
      double TAN2  = TAN*TAN;
      double num   = 1.0;
      double denom = 1.0 + 2.0*( 1.0 + Nu*Nu/Q2 )*TAN2;
      double eps   = -1;
      if(denom!=0){
	 eps = num/denom;
      }else{
	 std::cout << "[Kinematics::GetEpsilon]: ERROR! denominator = 0!" << std::endl;
      }
      return eps;
   }
}
