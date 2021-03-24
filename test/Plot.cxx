// simple environment to test code out
// NOTE: This relies on the util_df library (compiles into ROOT)   

#include <cstdlib>
#include <iostream>
#include <vector>

#include "TString.h"
#include "TLine.h"

#include "UtilDFGraph.hh"
#include "CSVManager.hh"

TGraph *GetGraph(util_df::CSVManager *data,std::string xAxis,std::string yAxis); 
TGraphErrors *GetTGraphErrors(const char *xAxis,bool isBorn,util_df::CSVManager *data); 

TGraph *GetError(const char *yAxis,util_df::CSVManager *MS,TGraph *rc); 
TGraph *GetDiff(const char *yAxis,util_df::CSVManager *data,TGraph *model); 

int Plot(){

   double Es=0;
   std::cout << "Enter beam energy (1, 5, 20 GeV): ";
   std::cin >> Es; 

   // from Phys. Rev. D 12, 1884 (1975) 
   // double tb = 5.459E-3; 
   // double ta = 11.129E-3; 

   double th = 5.;    // in deg 
   double EpMin=0,EpMax=0;
   int ES_INT = (int)Es; 

   if( ES_INT==1 ){
      EpMin = 0.100; 
      EpMax = 0.99; 
   }else if( ES_INT==5 ) {
      EpMin = 0.100; 
      EpMax = 4.901; 
   }else if( ES_INT==20 ){
      EpMin = 1; 
      EpMax = 18.6; 
   }

   double xMin = EpMin*0.8;  
   double xMax = EpMax*1.1;  

   // note: Mo & Tsai is in mub/GeV/sr
   char inpath[200];    
   sprintf(inpath,"./input/mo-and-tsai_Es-%d.csv",ES_INT);
   util_df::CSVManager *myCSV = new util_df::CSVManager(); 
   myCSV->ReadFile(inpath,true);
   // myCSV->Print();

   // note: output is in mub/GeV/sr
   char inpath_rc[200];
   sprintf(inpath_rc,"./output/rc-test_Es-%d.csv",ES_INT); 
   util_df::CSVManager *myRC = new util_df::CSVManager(); 
   myRC->ReadFile(inpath_rc,true);
   // myRC->Print();

   TGraph *gMSe = GetGraph(myCSV,"Ep","exact"    ); 
   TGraph *gMSa = GetGraph(myCSV,"Ep","ms-approx"); 

   TGraph *gRCe = GetGraph(myRC,"Ep(GeV)","tail(mub/GeV/sr)"       ); 
   TGraph *gRCa = GetGraph(myRC,"Ep(GeV)","tail_approx(mub/GeV/sr)");

   util_df::Graph::SetParameters(gMSe,20,kBlack); 
   util_df::Graph::SetParameters(gMSa,21,kRed  ); 

   util_df::Graph::SetParameters(gRCe,20,kBlack,0.5,3); 
   util_df::Graph::SetParameters(gRCa,21,kRed  ,0.5,3); 
   gRCe->SetLineStyle(2);  
   gRCa->SetLineStyle(2);  

   TGraph *gErr    = GetError("exact"    ,myCSV,gRCe); 
   TGraph *gErrApp = GetError("ms-approx",myCSV,gRCa); 
   util_df::Graph::SetParameters(gErr   ,20,kBlack,1,3);  
   util_df::Graph::SetParameters(gErrApp,20,kRed  ,1,3);

   TMultiGraph *mg = new TMultiGraph(); 
   mg->Add(gRCe,"l"); 
   mg->Add(gRCa,"l");
   mg->Add(gMSa,"p");
   mg->Add(gMSe,"p");

   TMultiGraph *mgErr = new TMultiGraph(); 
   mgErr->Add(gErr   ,"lp"); 
   mgErr->Add(gErrApp,"lp");

   TString Title      = Form("Elastic Tail Test (E_{s} = %.1lf GeV, #theta = %.1lf#circ)",Es,th);
   TString xAxisTitle = Form("E_{p} (GeV)");
   TString yAxisTitle = Form("#sigma_{el} (#mub/GeV/sr)");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gRCe,"Exact"             ,"l"); 
   L->AddEntry(gRCa,"Approx"            ,"l");
   L->AddEntry(gMSe,"Mo & Tsai (exact)" ,"p");  
   L->AddEntry(gMSa,"Mo & Tsai (approx)","p");  

   TLine *zeroLine = new TLine(xMin,0,xMax,0); 
 
   TCanvas *c1 = new TCanvas("c1","Elastic Tail Test",1000,800); 
   c1->Divide(1,2); 

   c1->cd(1);
   gPad->SetLogy(); 
   mg->Draw("a");
   util_df::Graph::SetLabels(mg,Title,xAxisTitle,yAxisTitle); 
   util_df::Graph::SetLabelSizes(mg,0.05,0.06); 
   mg->GetXaxis()->SetLimits(xMin,xMax); 
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update();

   c1->cd(2);
   mgErr->Draw("a");
   util_df::Graph::SetLabels(mgErr,"Error Relative to Mo & Tsai",xAxisTitle,"(DF-MS)/MS (%)"); 
   util_df::Graph::SetLabelSizes(mgErr,0.05,0.06); 
   mgErr->GetXaxis()->SetLimits(xMin,xMax); 
   mgErr->GetYaxis()->SetRangeUser(-40,40); 
   mgErr->Draw("a");
   zeroLine->Draw("same"); 
   c1->Update();
 
   delete myCSV;  
   delete myRC; 

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(const char *xAxis,bool isBorn,util_df::CSVManager *data){
   std::vector<double> x,y,ey,stat,syst;
   data->GetColumn_byName<double>(xAxis,x);
   if(isBorn){ 
      data->GetColumn_byName<double>("xs_born" ,y); 
      data->GetColumn_byName<double>("xsb_stat",stat); 
      data->GetColumn_byName<double>("xsb_syst",syst);
   }else{
      data->GetColumn_byName<double>("xs_rad"  ,y); 
      data->GetColumn_byName<double>("xsr_stat",stat); 
      data->GetColumn_byName<double>("xsr_syst",syst);
   } 

   double arg=0;
   const int N = x.size();
   for(int i=0;i<N;i++){
      arg = TMath::Sqrt( stat[i]*stat[i] + syst[i]*syst[i] );
      ey.push_back(arg);  
   }

   TGraphErrors *g = util_df::Graph::GetTGraphErrors(x,y,ey); 
   return g;
}
//______________________________________________________________________________
TGraph *GetDiff(const char *yAxis,util_df::CSVManager *data,TGraph *model){
   // define error as (model-data)/data 
   std::vector<double> x,y,ERR; 
   data->GetColumn_byName<double>("Ep" ,x);  
   data->GetColumn_byName<double>(yAxis,y);
   double f=0,err=0;
   const int N = x.size();
   for(int i=0;i<N;i++){
      f   = model->Eval(x[i]);
      err = 100.*(f-y[i])/y[i]; 
      ERR.push_back(err);
      // std::cout << Form("Ep = %.3lf GeV, err = %.3lf",x[i],err) << std::endl; 
   }
   TGraph *g = util_df::Graph::GetTGraph(x,ERR);
   return g; 
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

   const int N = y.size(); 
   for(int i=0;i<N;i++) y[i] *= 1E+3; // mub -> nb 

   TGraph *g = util_df::Graph::GetTGraph(x,y);
   return g;
} 
