// simple environment to test code out 

#include <cstdlib>
#include <iostream>
#include <vector>

#include "UtilDFGraph.hh"

#include "./src/DipoleElasticFF.C"
#include "./src/GMp12ElasticFF.C"
// #include "./src/RadiativeCorrections.C"

int TestFF(){

   DipoleElasticFF *myDip = new DipoleElasticFF(); 
   GMp12ElasticFF *myFit  = new GMp12ElasticFF(); 

   std::vector<double> Q2,GE,GM,GE12,GM12,R,R12; 

   int NPTS = 200; 
   double min=0,max=30; 
   double step = (max-min)/( (double)NPTS );
   double ge=0,gm=0,iQ2=0; 
   for(int i=0;i<NPTS;i++){
      iQ2 = min + ( (double)i )*step; 
      // dipole 
      ge = myDip->GetGE(iQ2); 
      gm = myDip->GetGM(iQ2);
      Q2.push_back(iQ2); 
      GE.push_back(ge);  
      GM.push_back(gm); 
      R.push_back(mu_p*ge/gm);  
      // GMp12 fit 
      ge = myFit->GetGE(iQ2); 
      gm = myFit->GetGM(iQ2);
      GE12.push_back(ge);  
      GM12.push_back(gm);  
      R12.push_back(mu_p*ge/gm);  
   } 

   TGraph *Ge = util_df::Graph::GetTGraph(Q2,GE); 
   TGraph *Gm = util_df::Graph::GetTGraph(Q2,GM);
   TGraph *Gr = util_df::Graph::GetTGraph(Q2,R);

   TGraph *Ge12 = util_df::Graph::GetTGraph(Q2,GE12); 
   TGraph *Gm12 = util_df::Graph::GetTGraph(Q2,GM12);
   TGraph *Gr12 = util_df::Graph::GetTGraph(Q2,R12);

   util_df::Graph::SetParameters(Ge,20,kRed ,0.5,3);  
   util_df::Graph::SetParameters(Gm,20,kBlue,0.5,3);  
   util_df::Graph::SetParameters(Gr,20,kBlack,0.5,3);  

   util_df::Graph::SetParameters(Ge12,20,kRed ,0.5,3);  
   util_df::Graph::SetParameters(Gm12,20,kBlue,0.5,3);  
   util_df::Graph::SetParameters(Gr12,20,kBlack,0.5,3);  

   Ge12->SetLineStyle(2); 
   Gm12->SetLineStyle(2); 
   Gr12->SetLineStyle(2); 

   TMultiGraph *mge = new TMultiGraph();
   mge->Add(Ge  ,"l"); 
   mge->Add(Ge12,"l"); 

   TMultiGraph *mgm = new TMultiGraph();
   mgm->Add(Gm  ,"l"); 
   mgm->Add(Gm12,"l"); 

   TMultiGraph *mgr = new TMultiGraph();
   mgr->Add(Gr  ,"l"); 
   mgr->Add(Gr12,"l"); 

   double xSize = 0.05; 
   double ySize = 0.07; 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(Ge  ,"Dipole"   ,"l"); 
   L->AddEntry(Ge12,"GMp12 Fit","l"); 

   TCanvas *c1 = new TCanvas("c1","GE and GM",800,800);
   c1->Divide(1,3); 

   c1->cd(1);
   mge->Draw("a");
   util_df::Graph::SetLabels(mge,"G_{E}","Q^{2} (GeV^{2})","G_{E}");  
   util_df::Graph::SetLabelSizes(mge,xSize,ySize); 
   mge->Draw("a");
   L->Draw("same");
   c1->Update();

   c1->cd(2);
   mgm->Draw("a");
   util_df::Graph::SetLabels(mgm,"G_{M}","Q^{2} (GeV^{2})","G_{M}");  
   util_df::Graph::SetLabelSizes(mgm,xSize,ySize); 
   mgm->Draw("a");
   c1->Update();

   c1->cd(3);
   mgr->Draw("a");
   util_df::Graph::SetLabels(mgr,"#mu_{p}G_{E}/G_{M}","Q^{2} (GeV^{2})","#mu_{p}G_{E}/G_{M}");  
   util_df::Graph::SetLabelSizes(mgr,xSize,ySize); 
   mgr->Draw("a");
   c1->Update();

   delete myDip; 
   delete myFit; 

   return 0;
} 
