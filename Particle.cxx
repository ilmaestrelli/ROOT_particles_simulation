#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include <cstdlib>
#include "Particle.h"

int Particle::fNParticleType_=0;

std::array<ParticleType*, fMaxNumParticleType> Particle::fParticleType_;

Particle::Particle(char* name, double Px, double Py, double Pz)
                                 : fPx_(Px), fPy_(Py), fPz_(Pz) 
                            {fIndex_ = FindParticle(name);};
                            
Particle::Particle() : fIndex_(-1), fPx_(0), fPy_(0), fPz_(0) {};
                            
int Particle::GetIndex() const {return fIndex_;}

/* void Particle::SetIndex(int fIndex)
{
  if (FindParticle(fIndex) != 10)
  {
    fIndex_ = fIndex;
  }
}; */

void Particle::SetParticleID(char* Name)
{
  if (FindParticle(Name) != 10)
  {
    fIndex_ = FindParticle(Name);
  }
}; 

/* int Particle::FindParticle(int fIndex)
{
    for (int i = 0; i < fNParticleType_; i++)
    {
        if (fParticleType_[i]->GetIndex()==fIndex)
        {
            return i;
        }
    }
    std::cout << fIndex << " not found" << '\n';
    return fNParticleType_ + 1;
};  */
 
int Particle::FindParticle(char* Name)
{ 
    for (int i = 0; i < fNParticleType_; i++)
    {
        if (fParticleType_[i]->GetName()==Name)
        {
            return i;
        }
        
        else 
        {
         std::cout << Name << " not found" << '\n';
        }
    }
    
    return fNParticleType_ + 1;
};

void Particle::AddParticleType(char* name, double mass, int charge, double width) 
{
  if (fNParticleType_ < fMaxNumParticleType && FindParticle(name) == 10) {
  
    if (width==0)  {
      fParticleType_[fNParticleType_] = new ParticleType(name, mass, charge);
      ++fNParticleType_;
    }
    else
    {
      fParticleType_[fNParticleType_] = new ResonanceType(name, mass, charge, width);
      ++fNParticleType_;
    }
  }
};


void Particle::PrintParticleID() const 
{
  std::cout<< "Name: " << GetName() << '\n';
  std::cout<< "Index: " << GetIndex() << '\n';
  std::cout<<"Impulse: P=("<<GetPx()<<", "<<GetPy()<<", "<<GetPz()<<")"<<'\n';
};

void Particle::PrintParticle()
{  
    for (int i = 0; i < fNParticleType_; ++i) {
      std::cout <<"                     "<<'\n'
                <<"Name: "<< fParticleType_[i]->GetName()<<'\n'
                <<"Mass: "<< fParticleType_[i]->GetMass()<<"GeV/c^2"<<'\n'
                <<"Charge: "<< fParticleType_[i]->GetCharge()<<'\n'
                <<"--------------------"<<'\n';
     }
};


double Particle::GetPx() const {return fPx_;}      
double Particle::GetPy() const {return fPy_;} 
double Particle::GetPz() const {return fPz_;} 


const char* Particle::GetName() const { return fParticleType_[fIndex_]->GetName(); }

double Particle::GetMass() const { return fParticleType_[fIndex_]->GetMass(); }

int Particle::GetCharge() const { return fParticleType_[fIndex_]->GetCharge(); }


double Particle::GetEnergy() const {
    double P = sqrt(pow(fPx_, 2) + pow(fPy_, 2) + pow(fPz_, 2));
 
    double E = sqrt(pow(GetMass(), 2) + pow(P, 2));
 
    return E;
};


double Particle::InvMass(Particle &p)  
{
    double P1 = sqrt(pow(fPx_, 2) + pow(fPy_, 2) + pow(fPz_, 2));
 
    double P2 = sqrt(pow(p.GetPx(), 2) + pow(p.GetPy(), 2) + pow(p.GetPz(), 2));
 
    double M_inv = sqrt(pow((GetEnergy() + p.GetEnergy()), 2) + pow((P1 + P2), 2));
    return M_inv;
};

void Particle::SetP(double Px, double Py, double Pz) 
{
    fPx_ = Px;
    fPy_ = Py;
    fPz_ = Pz;
};


Particle::~Particle() 
{
  for (int i = 0; i != Particle::fNParticleType_; i++) {
      delete Particle::fParticleType_[i];
    }
};          
                           
                           
//////////DECAY PROCESS//////////////////////////

int Particle::Decay2body(Particle &dau1,Particle &dau2) {
  if(GetMass() == 0.0){
    std::cout<<"Decayment cannot be preformed if mass is zero"<<'\n';
    return 1;
  }
  
  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if(fIndex_ > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleType_[fIndex_]->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    std::cout<<"Decayment cannot be preformed because mass is too low in this channel"<<'\n';
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx_*fPx_ + fPy_*fPy_ + fPz_*fPz_ + massMot*massMot);

  double bx = fPx_/energy;
  double by = fPy_/energy;
  double bz = fPz_/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz)
{

  double energy = GetEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx_ + by*fPy_ + bz*fPz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx_ += gamma2*bp*bx + gamma*bx*energy;
  fPy_ += gamma2*bp*by + gamma*by*energy;
  fPz_ += gamma2*bp*bz + gamma*bz*energy;
}           

                  
