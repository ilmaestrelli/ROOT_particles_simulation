#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

#include <iostream>

int main(){

    gStyle->SetOptFit(111);
    gStyle->SetOptStat(210);

    TFile *file = TFile::Open("analysis.root");
    TFile *flMC_file = new TFile("local.root");
        
    auto c1bis = new TCanvas();
    auto c2bis = new TCanvas();
    c1bis->Divide(2, 2);
    c2bis->Divide(1, 3);
       
    double pi = M_PI;
    
    TH1F *h[12];
    h[0] = (TH1F *)file->Get("hTypes");
    h[1] = (TH1F *)file->Get("hAzimuthal");
    h[2] = (TH1F *)file->Get("hPolar");
    h[3] = (TH1F *)file->Get("hP");
    h[4] = (TH1F *)file->Get("hTransverseP");
    h[5] = (TH1F *)file->Get("hEnergy");
    h[6] = (TH1F *)file->Get("hInvMass");
    h[7] = (TH1F *)file->Get("hInvMassDis");
    h[8] = (TH1F *)file->Get("hInvMassCon");
    h[9] = (TH1F *)file->Get("hInvMassPiKDis");
    h[10] = (TH1F *)file->Get("hInvMassPiKCon");
    h[11] = (TH1F *)file->Get("hInvMassDec");
    
    flMC_file->GetList()->Write();
    
  
    //settings axis
    h[0]->GetXaxis()->SetTitle("Popolazioni");
    h[1]->GetXaxis()->SetTitle("Azimuthal angle (rad)");
    h[2]->GetXaxis()->SetTitle("Polar angle (rad)");
    h[3]->GetXaxis()->SetTitle("Impulse module (GeV)");
    h[4]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    h[5]->GetXaxis()->SetTitle("Invariant Mass(GeV/c^2)");
    h[6]->GetXaxis()->SetTitle("Invariant Mass(GeV/c^2)");
    
    for (int i=0; i<7; ++i) 
     {
        h[i]->GetYaxis()->SetTitle("Occurrences");
     }
    
    //c1
    c1bis->cd(1);
    h[0]->Draw();
    c1bis->cd(4);
    h[1]->Draw();
    c1bis->cd(3);
    h[2]->Draw();
    c1bis->cd(2);
    h[3]->Draw();
    
    //c2
    c2bis->cd(1);
    h[11]->Draw();
    c2bis->cd(2);
    h13->Draw();
    c2bis->cd(3);
    h14->Draw();

    c1bis->Print("simulation1.pdf");
    c2bis->Print("simulation2.pdf");
    
}
