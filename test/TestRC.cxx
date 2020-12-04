// simple environment to test code out 

#include <cstdlib>
#include <iostream>
#include <vector>

#include "UtilDFGraph.hh"
#include "CSVManager.hh"

#include "./src/DipoleElasticFF.C"
#include "./src/GMp12ElasticFF.C"
#include "./src/Kinematics.C"
#include "./src/RadiativeCorrections.C"

#include "TStopwatch.h"

TGraph *GetGraph(util_df::CSVManager *data,std::string xAxis,std::string yAxis); 
TGraph *GetError(const char *yAxis,util_df::CSVManager *MS,TGraph *rc); 

int TestRC(){

   double Es=0;
   std::cout << "Enter beam energy (1, 5, 20 GeV): ";
   std::cin >> Es; 

   // from Phys. Rev. D 12, 1884 (1975) 
   double tb    = 5.459E-3; 
   double ta    = 11.129E-3; 
   double th    = 5.;    // in deg 
   double EpMin = 0.100; // in GeV
   double EpMax = 4.901; // in GeV 
   double Z  = 1;  
   double A  = 1;  

   int ES_INT = (int)Es; 

   if( ES_INT==1 ){
      EpMin = 0.100; 
      EpMax = 0.98; 
   }else if( ES_INT==5 ) {
      EpMin = 0.100; 
      EpMax = 4.901; 
   }else if( ES_INT==20 ){
      EpMin = 1; 
      EpMax = 18.6; 
   }

   double xMin = EpMin*0.8;  
   double xMax = EpMax*1.1;  

   char inpath[20];    
   sprintf(inpath,"./input/mo-and-tsai_Es-%d.csv",(int)Es);

   // note: Mo & Tsai is in mub/GeV/sr
   util_df::CSVManager *myCSV = new util_df::CSVManager(); 
   myCSV->ReadFile(inpath,true);
   // myCSV->Print();

   TGraph *gMSe = GetGraph(myCSV,"Ep","exact"    ); 
   TGraph *gMSa = GetGraph(myCSV,"Ep","ms-approx"); 

   util_df::Graph::SetParameters(gMSe,20,kBlack); 
   util_df::Graph::SetParameters(gMSa,21,kRed  ); 

   DipoleElasticFF *myDip = new DipoleElasticFF(); 

   RadiativeCorrections *myRC = new RadiativeCorrections(); 
   myRC->SetFormFactor(myDip);
   // myRC->SetRadiationLengths(tb,ta);
   myRC->SetTargetVariables(Z,A);
   // myRC->SetVerbosity(1);  

   TStopwatch *sw = new TStopwatch(); 
   double tStart=0,tStop=0,dt=0; 

   double min=1E-7,max=0;

   const int N = 100; 
   std::vector<double> EP,ET,ETA; 
   double Ep=0,el_tail=0,el_tail_approx=0;
   double step = (EpMax - EpMin)/( (double)N ); 
   for(int i=0;i<N;i++){
      Ep = EpMin + ( (double)i )*step;
      myRC->SetKinematicVariables(Es,Ep,th); 
      el_tail        = MUB_PER_GEV_SR*myRC->ElasticTail_exact();      // need to convert to mub/GeV/sr! 
      el_tail_approx = MUB_PER_GEV_SR*myRC->ElasticTail_peakApprox(); // need to convert to mub/GeV/sr!
      if(el_tail<min)        min = el_tail;  
      if(el_tail_approx<min) min = el_tail_approx;  
      if(el_tail>max)        max = el_tail;  
      if(el_tail_approx<max) max = el_tail_approx;  
      EP.push_back(Ep);  
      ET.push_back(el_tail);  
      ETA.push_back(el_tail_approx);
      // std::cout << Form("Es = %.3lf GeV, Ep = %.3lf GeV, th = %.1lf deg, el_tail = %.3E, el_tail_approx = %.3E",Es,Ep,th,el_tail,el_tail_approx) << std::endl; 
   }

   max *= 1.2; 

   TGraph *gExact  = util_df::Graph::GetTGraph(EP,ET ); 
   TGraph *gApprox = util_df::Graph::GetTGraph(EP,ETA);
   util_df::Graph::SetParameters(gExact ,20,kBlack,0.5,3);  
   util_df::Graph::SetParameters(gApprox,20,kRed  ,0.5,3);
   gExact->SetLineStyle(2);  
   gApprox->SetLineStyle(2);  

   TGraph *gErr    = GetError("exact"    ,myCSV,gExact ); 
   TGraph *gErrApp = GetError("ms-approx",myCSV,gApprox); 
   util_df::Graph::SetParameters(gErr   ,20,kBlack,1,3);  
   util_df::Graph::SetParameters(gErrApp,20,kRed  ,1,3);

   TMultiGraph *mg = new TMultiGraph(); 
   mg->Add(gExact ,"l"); 
   mg->Add(gApprox,"l");
   mg->Add(gMSa   ,"p");
   mg->Add(gMSe   ,"p");

   TMultiGraph *mgErr = new TMultiGraph(); 
   mgErr->Add(gErr   ,"lp"); 
   mgErr->Add(gErrApp,"lp");

   TString Title      = Form("Elastic Tail Test (E_{s} = %.1lf GeV, #theta = %.1lf#circ)",Es,th);
   TString xAxisTitle = Form("E_{p} (GeV)");
   TString yAxisTitle = Form("#sigma_{el} (#mub/GeV/sr)");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gExact ,"Exact" ,"l"); 
   L->AddEntry(gApprox,"Approx","l");
   L->AddEntry(gMSe   ,"Mo & Tsai (exact)" ,"p");  
   L->AddEntry(gMSa   ,"Mo & Tsai (approx)","p");  
 
   TCanvas *c1 = new TCanvas("c1","Elastic Tail Test",1000,800); 
   c1->Divide(1,2); 

   c1->cd(1);
   mg->Draw("a");
   util_df::Graph::SetLabels(mg,Title,xAxisTitle,yAxisTitle); 
   util_df::Graph::SetLabelSizes(mg,0.05,0.06); 
   mg->GetXaxis()->SetLimits(xMin,xMax); 
   mg->GetYaxis()->SetRangeUser(min,max); 
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgErr->Draw("a");
   util_df::Graph::SetLabels(mgErr,"Error Relative to Mo & Tsai",xAxisTitle,"(DF-MS)/MS (%)"); 
   util_df::Graph::SetLabelSizes(mgErr,0.05,0.06); 
   mgErr->GetXaxis()->SetLimits(xMin,xMax); 
   // mgErr->GetYaxis()->SetRangeUser(0,100); 
   mgErr->Draw("a");
   c1->Update();
 
   delete myDip;
   delete myRC; 
   delete myCSV;  

   return 0;
}
//______________________________________________________________________________
TGraph *GetError(const char *yAxis,util_df::CSVManager *MS,TGraph *rc){
   // define error as: (my_calc-ms)/ms 
   // my_calc = my code 
   // ms      = mo & tsai
   std::vector<double> x,y,ERR; 
   MS->GetColumn_byName<double>("Ep" ,x);  
   MS->GetColumn_byName<double>(yAxis,y); 

   double f=0,err=0;
   const int N = x.size();
   for(int i=0;i<N;i++){
      f = rc->Eval(x[i]);
      err = 100.*(f-y[i])/y[i]; 
      ERR.push_back(err);
      // std::cout << Form("Ep = %.3lf GeV, err = %.3lf",x[i],err) << std::endl; 
   }
   TGraph *g = util_df::Graph::GetTGraph(x,ERR);
   return g; 
}
//______________________________________________________________________________
TGraph *GetGraph(util_df::CSVManager *data,std::string xAxis,std::string yAxis){
   std::vector<double> x,y; 
   data->GetColumn_byName<double>(xAxis,x);
   data->GetColumn_byName<double>(yAxis,y);
   TGraph *g = util_df::Graph::GetTGraph(x,y);
   return g;
} 
