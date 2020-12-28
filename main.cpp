
#include "particle.h"
#include <iostream>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TRandom.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TStyle.h"

int main()
{
    gROOT->LoadMacro("particletype.cpp");
    gROOT->LoadMacro("resonancetype.cpp");
    gROOT->LoadMacro("particle.cpp");
    
    //root file
    TFile *file = new TFile("analysis.root", "RECREATE");
    const int nGen = 100000; //number of generations   
    double phi, theta, impulse;
    double pi = M_PI; 
    
   //adding particle to the array
    Particle::AddParticleType("Pion+", 0.13957, 1);          
    Particle::AddParticleType("Pion-", 0.13957, -1);          
    Particle::AddParticleType("Kaon+", 0.49367, 1);          
    Particle::AddParticleType("Kaon-", 0.49367, -1);          
    Particle::AddParticleType("Proton+", 0.93827, 1);  
    Particle::AddParticleType("Proton-", 0.93827, -1);   
    Particle::AddParticleType("K*", 0.89166, 0, 0.050); 
    
    auto hTypes = new TH1F("hTypes", "PARTICLE TYPES DISTRIBUTION", 8, 0, 8);                           
    auto hAzimuthal = new TH1F("hAzimuthal", "AZIMUTHAL ANGLE DISTRIBUTION (Phi)", 1000, 0, 2*pi); 
    auto hPolar = new TH1F("hPolar", "POLAR ANGLE DISTRIBUTION (Theta)", 1000, 0, pi);         
    auto hP = new TH1F("hP", "IMPULSE DISTRIBUTION", 1000, 0, 5.5);                        
    auto hTransverseP = new TH1F("hTransverseP", "TRANSVERSE IMPULSE DISTRIBUTION", 1000, 0, 4);                
    auto hEnergy = new TH1F("hEnergy", "ENERGY DISTRIBUTION", 1000, 0, 6);                            
    auto hInvMass = new TH1F("hInvMass", "INVARIANT MASS DISTRIBUTION", 1000, 0, 3);             
    auto hInvMassDis = new TH1F("hInvMassDis", "INVARIANT MASS DISTRIBUTION OF PARTICLES WITH OPPOSITE CHARGE", 500, 0, 3); 
    auto hInvMassCon = new TH1F("hInvMassCon", "INVARIANT MASS DISTRIBUTION OF PARTICLES WITH CONCORDE CHARGE", 500, 0, 3);      
    auto hInvMassPiKDis = new TH1F("hInvMassPiKDis", "INVARIANT MASS DISTRIBUTION P+/K-, P-/K+", 500, 0, 3);         
    auto hInvMassPiKCon = new TH1F("hInvMassPiKCon", "INVARIANT MASS DISTRIBUTION P+/K+, P-/K-", 500, 0, 3);              
    auto hInvMassDec = new TH1F("hInvMassDec ", "INVARIANT MASS DISTRIBUTION (K*)", 500, 0, 3);   
    hInvMassDec->Sumw2();
    
    gStyle->SetOptFit(111);        

    gRandom->SetSeed();

    for (int i = 0; i<nGen; i++)
    { 
          
    Particle particles[120];
    int nResonance = 0;  
    
        for (int j=0; j<100; j++)
        { 
            phi = gRandom->Uniform(2*pi);
            hAzimuthal->Fill(phi);
            
            theta = gRandom->Uniform(pi);
            hPolar->Fill(theta);
            
            impulse = gRandom->Exp(1);
            hP->Fill(impulse);

            double Px, Py, Pz;
            Px = impulse*sin(theta)*cos(phi);
            Py = impulse*sin(theta)*sin(phi);
            Pz = impulse*cos(theta);
            particles[j].SetImpulse(Px,Py,Pz);

            double TransverseImpulse = sqrt(pow(Px,2) + pow(Py,2));
            hTransverseP->Fill(TransverseImpulse);

            //particle abundancies
            double p = gRandom->Uniform(1);
            if (p < 0.40) { particles[j].SetParticleID("Pion+"); hTypes->Fill(0);}
                    
            else if (p < 0.80)
            {particles[j].SetParticleID("Pion-"); hTypes->Fill(1);}
            
            else if (p < 0.85)
            { particles[j].SetParticleID("Kaon+"); hTypes->Fill(2);}
            
            else if (p < 0.90)
            {particles[j].SetParticleID("Kaon-"); hTypes->Fill(3);}
            
            else if (p < 0.945)
            {particles[j].SetParticleID("Proton+"); hTypes->Fill(4);}
            
            else if (p < 0.99)
            {particles[j].SetParticleID("Proton-"); hTypes->Fill(5);}
            
            else if (p < 0.995)
            { 
                particles[j].SetParticleID("K*"); hTypes->Fill(6);
                particles[100 + nResonance].SetParticleID("Pion+");
                particles[100 + nResonance + 1].SetParticleID("Kaon-");
                particles[j].Decay2body(particles[100 + nResonance], particles[100 + nResonance + 1]);
                nResonance++;
                nResonance++;
            }
            else
            { 
                particles[j].SetParticleID("K*"); hTypes->Fill(6);
                particles[100 + nResonance].SetParticleID("Pion-");
                particles[100 + nResonance + 1].SetParticleID("Kaon+");
                particles[j].Decay2body(particles[100 + nResonance], particles[100 + nResonance + 1]);
                nResonance++;
                nResonance++;
            }

            double energy = particles[j].GetEnergy();
            hEnergy->Fill(energy);
        }


        for (int k=0; k<100+nResonance-1; k++)
        {
            for (int n=k+1; n<100 + nResonance; n++)
            {
                hInvMass->Fill(particles[k].GetInvMass(particles[n])); 

                //opposite charge combination
                if (((particles[k].GetIndex() == 0 || 
                      particles[k].GetIndex() == 2 ||
                      particles[k].GetIndex() == 4) && (particles[n].GetIndex() == 1 || 
                      particles[n].GetIndex() == 3 || 
                      particles[n].GetIndex() == 5)) ||
                      
                    ((particles[k].GetIndex() == 1 ||
                      particles[k].GetIndex() == 3 ||
                      particles[k].GetIndex() == 5) && (particles[n].GetIndex() == 0 || 
                      particles[n].GetIndex() == 2 || 
                      particles[n].GetIndex() == 4)))
                {
                    hInvMassDis->Fill(particles[k].GetInvMass(particles[n]));
                }

                //concorde charge combination
                if (((particles[k].GetIndex() == 0 || 
                      particles[k].GetIndex() == 2 ||
                      particles[k].GetIndex() == 4) && (particles[n].GetIndex() == 0 || 
                      particles[n].GetIndex() == 2 ||
                      particles[n].GetIndex() == 4)) ||
                      
                    ((particles[k].GetIndex() == 1 || 
                      particles[k].GetIndex() == 3 ||
                      particles[k].GetIndex() == 5) && (particles[n].GetIndex() == 1 || 
                      particles[n].GetIndex() == 3 || 
                      particles[n].GetIndex() == 5)))
                {
                    hInvMassCon->Fill(particles[k].GetInvMass(particles[n]));
                }

                //P+/K-, P-/K+
                if (((particles[k].GetIndex() == 0) && (particles[n].GetIndex() == 3)) ||
                    ((particles[k].GetIndex() == 1) && (particles[n].GetIndex() == 2)))
                {
                    hInvMassPiKDis->Fill(particles[k].GetInvMass(particles[n]));
                }

                //P+/K+, P-/K-
                if (((particles[k].GetIndex() == 0) && (particles[n].GetIndex() == 2)) ||
                    ((particles[k].GetIndex() == 1) && (particles[n].GetIndex() == 3)))
                {
                    hInvMassPiKCon->Fill(particles[k].GetInvMass(particles[n]));
                }
            }
        }

        if (nResonance!=0)
        {
            for (int j=0; j<nResonance; j+=2)
            {
                hInvMassDec->Fill(particles[100+j].GetInvMass(particles[100+1+j]));
            }
        }
    } //end of 10000 events generation
    
    //canvas
    //c1
    auto c1 = new TCanvas();
    c1->Divide(3, 2);
    c1->cd(1);
    hTypes->Draw();
    c1->cd(2);
    hAzimuthal->Draw();
    c1->cd(3);
    hPolar->Draw();
    c1->cd(4);
    hP->Draw();
    c1->cd(5);
    hTransverseP->Draw();
    c1->cd(6);
    hEnergy->Draw();
    
    //c2
    auto c2 = new TCanvas();
    c2->Divide(3, 2);
    c1->cd(1);
    hTypes->Draw();
    c1->cd(2);
    hAzimuthal->Draw();
    c1->cd(3);
    hPolar->Draw();
    c1->cd(4);
    hP->Draw();
    c1->cd(5);
    hTransverseP->Draw();
    c1->cd(6);
    hEnergy->Draw();
    c2->cd(1);
    hInvMass->Draw();
    c2->cd(2);
    hInvMassDis->Draw();
    c2->cd(3);
    hInvMassCon->Draw();
    c2->cd(4);
    hInvMassPiKDis->Draw();
    c2->cd(5);
    hInvMassPiKCon->Draw();
    
    //resonance fit
    auto fitRes = new TF1("K* FIT", "gaus", 0.60, 1.30); //gaussian 
    hInvMassDec->Fit(fitRes);

    c2->cd(6);
    hInvMassDec->Draw();
    
    c1->Print("c1.pdf");
    c1->Print("c1.C");
    c2->Print("c2.pdf");
    c2->Print("c2.C");

    //Cheking bins contents and errors
    const double total= 100*nGen;
    std::cout << "Pion+   are " << hTypes->GetBinContent(1)/total<< " ± " << hTypes->GetBinError(1)/total<< ",  0.40 expected"<<'\n';
    std::cout << "Pion-   are " << hTypes->GetBinContent(2)/total<< " ± " << hTypes->GetBinError(2)/total<< ",  0.40 expected"<<'\n';
    std::cout << "Kaon+   are " << hTypes->GetBinContent(3)/total<< " ± " << hTypes->GetBinError(3)/total << ", 0.05 expected"<<'\n';
    std::cout << "Kaon-   are " << hTypes->GetBinContent(4)/total << " ± " << hTypes->GetBinError(4)/total << ",0.05 expected"<<'\n';
    std::cout << "Proton+ are " << hTypes->GetBinContent(5)/total<< " ± " << hTypes->GetBinError(5)/total<< ",  0.045 expected"<<'\n';
    std::cout << "Proton- are " << hTypes->GetBinContent(6)/total<< " ± " << hTypes->GetBinError(6)/total<< ",  0.045 expected"<<'\n';
    std::cout << "K*      are " << hTypes->GetBinContent(7)/total<< " ± " << hTypes->GetBinError(7)/total<< ",  0.01 expected"<<'\n';
    
    //c3
    //phi and theta fit
    TCanvas *c3 = new TCanvas("c3", "ANGLES FIT");
    c3->Divide(2, 1);
    auto fitPhi = new TF1("AZIMUTHAL ANGLE FIT", "pol0", 0, 2 * pi); //pol0
    auto fitTheta = new TF1("POLAR ANGLE FIT", "pol0", 0, pi);
    c3->cd(1);
    hAzimuthal->Fit(fitPhi);
    c3->cd(2);
    hPolar->Fit(fitTheta);
    c3->Print("c3.pdf");
    c3->Print("c3.C");

    //c4
    //impulse fit
    TCanvas *c4 = new TCanvas("c4", "IMPULSE FIT");
    auto fitImpulse = new TF1("IMPULSE FIT", "expo", 0, 5.5); //exponential
    c4->cd(1);
    hP->Fit(fitImpulse);
    c4->Print("c4.pdf");
    c4->Print("c4.C");
    
    //c5
    TCanvas *c5 = new TCanvas("c5", "IMPULSE FIT");
    c5->Divide(2, 3);
    c5->cd(1);
    hInvMassPiKDis->Draw();
    c5->cd(2);
    hInvMassPiKCon->Draw();
    c5->cd(3);
    hInvMassDec->Draw();
    c5->cd(4);

    auto fitDecay = new TF1("K* FIT", "gaus", 0.60, 1.30); //gaussian
    hInvMassDec->Fit(fitDecay);

    //subctractions
    auto hSubPiK = new TH1F("hSubPiK", "INVARIANT MASS DISTRIBUTION (disc./conc. Pions-Kaons)", 500, 0, 3);
    hSubPiK->Sumw2();
    hSubPiK->Add(hInvMassPiKDis, hInvMassPiKCon, 1, -1);
    auto fitSubPiK = new TF1("Fit subtraction Pion Kaon", "gaus", 0.6, 3); //gaussian
    hSubPiK->Fit(fitSubPiK);
    hSubPiK->SetEntries(hSubPiK->Integral());
     
    c5->cd(5);
  
    auto hSubCharge = new TH1F("hSubCharge", "INVARIANT MASS DISTRIBUTION (disc./conc. charge)", 500, 0, 3);
    hSubCharge->Sumw2();
    hSubCharge->Add(hInvMassDis, hInvMassCon, 1, -1);
    auto fitSubCharge= new TF1("Fit subtraction invariant masses (charge)", "gaus", 0.6, 3); //gaussian
    hSubCharge->Fit(fitSubCharge);
    hSubCharge->SetEntries(hSubCharge->Integral());

    c5->Print("c5.pdf");
    c5->Print("c5.C");
    
    //file writing
    file->Write();
    file->Close();
    
    std::cout << "Simulation ended";
} 
