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
        
    auto c1bis = new TCanvas();
    auto c2bis = new TCanvas();
    c1bis->Divide(2, 2);
    c2bis->Divide(1, 3);
    
    TH1F *h[14];
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
    h[12] = (TH1F *)file->Get("hSubCharge");
    h[13] = (TH1F *)file->Get("hSubPiK");
    
    //cosmetic
    h[0]->SetFillColor(kBlue);
    
    h[0]->GetXaxis()->SetTitle("Particle types");
    h[1]->GetXaxis()->SetTitle("Azimuthal angle (rad)");
    h[2]->GetXaxis()->SetTitle("Polar angle (rad)");
    h[3]->GetXaxis()->SetTitle("Impulse (GeV)");
    h[11]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    h[12]->GetXaxis()->SetTitle("Invariant Mass(GeV/c^2)");
    h[13]->GetXaxis()->SetTitle("Invariant Mass(GeV/c^2)");
    
    for (int i=0; i<14; ++i) 
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
    h[12]->Draw();
    c2bis->cd(3);
    h[13]->Draw();
     
    //print final canvas
    c1bis->Print("canvasRelazione1.pdf");
    c1bis->Print("canvasRelazione1.C");
    c1bis->Print("canvasRelazione1.root");
    c2bis->Print("canvasRelazione2.pdf");
    c2bis->Print("canvasRelazione2.C"); 
    c2bis->Print("canvasRelazione2.root"); 
}
