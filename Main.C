#include <string>
#include <cmath>
#include <iostream>

#include "Particle.h"

#include "TRandom.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TROOT.h"

int main() 
{

gROOT->LoadMacro("ParticleType.cxx");
gROOT->LoadMacro("ResonanceType.cxx");
gROOT->LoadMacro("Particle.cxx");

Particle particles[120];
int resonanceDecay = 0;  //counter for extra particles due to K* decay
double pi = M_PI;
 
 Particle::AddParticleType("Pion+", 0.13957, 1, 0);
 Particle::AddParticleType("Pion-", 0.13957, -1, 0);
 Particle::AddParticleType("Kaon+", 0.49367, 1, 0);
 Particle::AddParticleType("Kaon-", 0.49367, -1, 0);
 Particle::AddParticleType("Proton+", 0.93827, 1, 0);
 Particle::AddParticleType("Proton-", 0.93827, -1, 0);
 Particle::AddParticleType("Kaon*", 0.899166, 1, 0.05);
  
double phi, theta=0;
int nGen=1000;

//root file
TFile *file = new TFile("analysis.root", "RECREATE");

//histograms
auto h1 = new TH1F("h1", "Particle type distribution", 7, 0, 7);  
auto h2 = new TH1F("h2", "Azimuthal angle distribution (Phi)", 500, 0, 2*pi);
auto h3 = new TH1F("h3", "Polar angle distribution (Theta)", 250, 0, pi);
auto h4 = new TH1F("h4", "Impulse", 800, 0, 5.5);     
auto h5 = new TH1F("h5", "Transverse impulse", 600, 0, 4);  
auto h6 = new TH1F("h6", "Energy", 600, 0, 5);  
auto h7 = new TH1F("h7", "Invariant masses", 400, 0, 5); //all particles
auto h8 = new TH1F("h8", "Invariant masses (different charge)", 400, 0, 5);
auto h9 = new TH1F("h9", "Invariant masses (same charge)", 400, 0, 5); 
auto h10 = new TH1F("h10", "Invariant masse P+/K-, P-/K+", 400, 0, 5); //different sign combination: Pion+/Kaon-, Pion-/Kaon+ 
auto h11 = new TH1F("h11", "Invariant masse P+/K+, P-/K-", 400, 0, 5); //same sign combination: Pion+/Kaon+, Pion-/Kaon-
auto h12 = new TH1F("h12", "Invariant masse K* decay", 400, 0.5, 1.2); //invariant masses among particles from resonance decay

auto h10bis = new TH1F("h10bis", "Invarian masse P+/K-, P-/K+", 5000, 0.5, 1.5); 
auto h11bis = new TH1F("h11bis", "Invariant masse P+/K+, P-/K-", 5000, 0.5, 1.5);
auto h12bis = new TH1F("h12bis", "Invariant masse K* Decay", 5000, 0.5, 1.5);  


gRandom->SetSeed();

for (int i=0; i<nGen; ++i)
{
  for (int j=0; j<100; ++j) 
  {
   
   //phi
   phi = gRandom->Uniform(0, 2*pi);
   h2->Fill(phi);
   
   //theta
   theta = gRandom->Uniform(0, TMath::Pi());
   h3->Fill(theta);
   
   //impulse
   double P = gRandom->Exp(1);
   h4->Fill(P);
   
   double Px = TMath::Sin(theta) * TMath::Cos(phi);
   double Py = TMath::Sin(theta) * TMath::Sin(phi);
   double Pz = TMath::Cos(theta);
   
   particles[j].SetP(Px, Py, Pz);
   
   //transverse impulse
   double transverseP = sqrt(Px*Px + Py*Py);
   h5->Fill(transverseP);
   
   
   double x = 0;
   x = gRandom->Rndm();
   
   if (x<0.4) {particles[j].SetParticleID("Pion+"); h1->Fill(0);}

   else if (x<0.8) {particles[j].SetParticleID("Pion-"); h1->Fill(1);}
           
   else if (x<0.85) {particles[j].SetParticleID("Kaon+"); h1->Fill(2);}
                    
   else if (x<0.9) {particles[j].SetParticleID("Kaon-"); h1->Fill(3);}
                  
   else if (x<0.945) {particles[j].SetParticleID("Proton+"); h1->Fill(4);}
                    
   else if (x<0.99) {particles[j].SetParticleID("Proton-"); h1->Fill(5);} 
   
   //K* -> P+ + K-
   else if (x<0.995) 
   {
    particles[j].SetParticleID("Kaon*");
    h1->Fill(6);
   
    particles[100+resonanceDecay].SetParticleID("Pion+");
    particles[100+resonanceDecay+1].SetParticleID("Kaon-");
    particles[j].Decay2body(particles[100+resonanceDecay], particles[100+resonanceDecay+1]); 
   
    resonanceDecay ++;
    resonanceDecay ++; 
   }
   
   //K* -> P- + K+
   else 
   {
    particles[j].SetParticleID("Kaon*");
    h1->Fill(6);
   
    particles[100+resonanceDecay].SetParticleID("Pion-");
    particles[100+resonanceDecay+1].SetParticleID("Kaon+");
    particles[j].Decay2body(particles[100+resonanceDecay], particles[100+resonanceDecay+1]); 
   
    resonanceDecay ++;
    resonanceDecay ++; 
   }
   
   
   //energy
   double E = particles[j].GetEnergy();
   h6->Fill(E);
   }
   
    for (int z = 0; z<100+resonanceDecay-1; ++z)
        {
            for (int n = z + 1; n<100+resonanceDecay; n++)
            {
                h7->Fill(particles[z].InvMass(particles[n])); 


                //different charge sign 
                if (((particles[z].GetIndex() == 0 ||
                   particles[z].GetIndex() == 2 ||
                   particles[z].GetIndex() == 4) && (particles[n].GetIndex() == 1 ||
                   particles[n].GetIndex() == 3 ||
                   particles[n].GetIndex() == 5)) ||
                
                ((particles[z].GetIndex() == 1 ||
                  particles[z].GetIndex() == 3 ||
                  particles[z].GetIndex() == 5) && (particles[n].GetIndex() == 0 ||
                  particles[n].GetIndex() == 2 ||
                  particles[n].GetIndex() == 4)))
                {
                    h8->Fill(particles[z].InvMass(particles[n]));
                }


                //same charge sign
                if (((particles[z].GetIndex() == 0 ||
                      particles[z].GetIndex() == 2 ||
                      particles[z].GetIndex() == 4) && (particles[n].GetIndex() == 0 ||
                      particles[n].GetIndex() == 2 ||
                      particles[n].GetIndex() == 4)) ||
                      
                    ((particles[z].GetIndex() == 1 ||
                      particles[z].GetIndex() == 3 ||
                      particles[z].GetIndex() == 5) && (particles[n].GetIndex() == 1 ||
                      particles[n].GetIndex() == 3 || 
                      particles[n].GetIndex() == 5)))
                {
                    h9->Fill(particles[z].InvMass(particles[n]));
                }

                //combination P+/K-, P-/K+ 
                if (((particles[z].GetIndex() == 0) && (particles[n].GetIndex() == 3)) ||
                    ((particles[z].GetIndex() == 1) && (particles[n].GetIndex() == 2)))
                {
                    h10->Fill(particles[z].InvMass(particles[n]));
                    h10bis->Fill(particles[z].InvMass(particles[n]));
                }

                //combination P+/K+, P-/K- 
                if (((particles[z].GetIndex() == 0) && (particles[n].GetIndex() == 2)) ||
                    ((particles[z].GetIndex() == 1) && (particles[n].GetIndex() == 3)))
                  {
                    h11->Fill(particles[z].InvMass(particles[n]));
                    h11bis->Fill(particles[z].InvMass(particles[n]));
                  }
                 }
               }
                
                if (resonanceDecay != 0)
                {
                   for (int j = 0; j < resonanceDecay; j += 2)
                   { 
                     h12->Fill(particles[100 + j].InvMass(particles[100 + 1 + j]));
                     h12bis->Fill(particles[100 + j].InvMass(particles[100 + 1 + j]));
                   }
                }   
              }  //end generation of 10000 events
   

    
    //check bin content and error of particle type distibution (h1)
    for(int i=1; i<=7; i++) 
    {
     std::cout<<"BIN CONTENT("<<i<<"):\t"<<h1->GetBinContent(i)<<'\n';
     std::cout<<"BIN ERROR("<<i<<"):\t"<<h1->GetBinError(i)<<'\n';
    }
   
    
    //substractions
    auto histoSubPiK= new TH1F("histoSubPiK","Subctraction Invariant mass Pion/Kaon",500,0,3);
    histoSubPiK->Add(h10, h11, 1, -1);
    auto histoSub= new TH1F("histoSub","Subctraction invariant mass of all particles",500,0,3);
    histoSub->Add(h8, h9, 1, -1);
    
    //cosmetic
    h1->GetXaxis()->SetTitle("PARTICLE TYPES");
    h1->GetYaxis()->SetTitle("ENTRIES");
    h2->GetXaxis()->SetTitle("AZIMUTAL ANGLE (rad)");
    h2->GetYaxis()->SetTitle("COUNTS");
    h3->GetXaxis()->SetTitle("POLAR ANGLE (rad)");
    h3->GetYaxis()->SetTitle("COUNTS");
    h4->GetXaxis()->SetTitle("MOMENTUM (GeV/c)");
    h4->GetYaxis()->SetTitle("COUNTS");
    h12->GetXaxis()->SetTitle("INVARIANT MASS (GeV/c^2)");
    h12->GetYaxis()->SetTitle("COUNTS (1MeV/c^2)");
    histoSub->GetXaxis()->SetTitle("INVARIANT MASS (GeV/c^2)");
    histoSub->GetYaxis()->SetTitle("COUNTS (1MeV/c^2)");
    histoSubPiK->GetXaxis()->SetTitle("INVARIANT MASS (GeV/c^2)");
    histoSubPiK->GetYaxis()->SetTitle("COUNTS (1MeV/c^2)");
    
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetYaxis()->SetTitleOffset(1.5);
    h2->GetXaxis()->SetTitleOffset(1.1);
    h2->GetYaxis()->SetTitleOffset(1.5);
    h3->GetXaxis()->SetTitleOffset(1.1);
    h3->GetYaxis()->SetTitleOffset(1.5);
    h4->GetXaxis()->SetTitleOffset(1.1);
    h4->GetYaxis()->SetTitleOffset(1.5);
    h12->GetXaxis()->SetTitleOffset(1.1);
    h12->GetYaxis()->SetTitleOffset(1.2);
    histoSub->GetXaxis()->SetTitleOffset(1.1);
    histoSub->GetYaxis()->SetTitleOffset(1.2);
    histoSubPiK->GetXaxis()->SetTitleOffset(1.1);
    histoSubPiK->GetYaxis()->SetTitleOffset(1.2);
    
    
    //Fit  
    gStyle->SetOptFit(1);
    h1->SetFillColor(kBlue);
    h4->SetLineWidth(1.8);
    
    auto fP = new TF1("fP","expo"); //exponential
    fP->SetLineColor(kRed); fP->SetLineWidth(0.3);
    h4->Fit(fP,"","",0,10);
    gStyle->SetOptFit(111);
    
    auto fUniform = new TF1("fUniform","pol0"); //uniform
    fUniform->SetLineColor(kRed); fUniform->SetLineWidth(0.5);
    h2->Fit(fUniform,"","",0,2*M_PI);
    h3->Fit(fUniform,"","",0,M_PI);
    h7->SetLineWidth(1.8);
    
    auto fGaus = new TF1("fGaus", "gaus");  //gaussian
    fGaus->SetLineColor(kRed); fGaus->SetLineWidth(0.3);
    h7->Fit("fGaus","","",0,3);
    histoSubPiK->Fit("fGaus","","",0,3);
    histoSub->Fit("fGaus","","",0,3);
    
    
   //draw and canvas
   auto c1 = new TCanvas("c1","c1",10,20,500,600);
   c1->Divide(3,2);
   c1->cd(1);
   h1->Draw();
   c1->cd(2);
   h4->Draw();
   c1->cd(3);
   h3->Draw();
   c1->cd(4);
   h2->Draw();
   c1-> Print("c1.pdf");
   
   auto c2 = new TCanvas("c2","c2",10,20,500,600);
   c2->Divide(3,2);
   c2->cd(1);
   h7->Draw();
   c2->cd(2);
   histoSub->Draw();
   c2->cd(3);
   histoSubPiK->Draw();
   c2-> Print("c2.pdf"); 
   
   auto c3 = new TCanvas("c3","c3",10,20,500,600);
   c3->Divide(2, 1);
   c3->cd(1);
   c3->cd(2);
   c3->Print("c3.pdf");
   
   auto c5 = new TCanvas("c5", "Impulse Fit",10,20,500,600);
   c5->Divide(2, 3);
   c5->cd(1);
   h10->Draw();
   c5->cd(2);
   h11->Draw();
   c5->cd(3);
   h12bis->Draw();
   c5->cd(4);
 
 
    
    //saving on file
    file->Write();
    file->Close();
}
