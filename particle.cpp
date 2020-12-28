#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

int Particle::fNParticleType_ = 0;
ParticleType *Particle::fParticleType_[fMaxNumParticleType_];

int Particle::GetIndex() const { return fIndex_; };

Particle::Particle(std::string name, double fPx, double fPy, double fPz) : fPx_{fPx}, fPy_{fPy}, fPz_{fPz}
{
    fIndex_ = FindParticle(name);
};

int Particle::FindParticle(std::string name)
{
    for (int i = 0; i < fNParticleType_; i++)
    {
        if (fParticleType_[i]->GetName() == name)
        {
            return i;
        }
    }
   // std::cout << name << " not found" << '\n';
    return fNParticleType_ + 1;
};

/*
  int Particle::FindParticle(fIndex_)
   {
    for (int i = 0; i < fNParticleType_; i++)
    {
        if (fParticleType_[i]->GetIndex() == index)
        {
            return i;
        }
    }
    std::cout << name << " not found" << '\n';
    return fNParticleType_ + 1;
    };
 */  

void Particle::AddParticleType(std::string name, double mass, int charge, double width)
{

    int index = FindParticle(name);

    if (width == 0)
    {
        if (index < fNParticleType_)
        {
            std::cout << "Particle " << name << " already in the array\n";
            return;
        }

        if (fNParticleType_ < fMaxNumParticleType_)
        {
            fParticleType_[fNParticleType_] = new ParticleType(name, mass, charge);
            fNParticleType_++;
        }
        else
        {
            std::cout << "Error: The array is full!"<<'\n';
        }
    }
    else
    {

        if (index < fNParticleType_)
        {
            std::cout << "Error: the particle is already in the array"<<'\n';
            return;
        }

        if (fNParticleType_ < fMaxNumParticleType_)
        {
            fParticleType_[fNParticleType_] = new ResonanceType(name, mass, charge, width);
            fNParticleType_++;
        }
        else
        {
            std::cout << "Error: The array is full!"<<'\n';
        }
    }
};


double Particle::GetPx() const { return fPx_; };
double Particle::GetPy() const { return fPy_; };
double Particle::GetPz() const { return fPz_; };
double Particle::GetImpulse() const { return sqrt(fPx_*fPx_ + fPy_*fPy_ + fPz_*fPz_); };

double Particle::GetMass() const
{
    int i = this->fIndex_;
    return fParticleType_[i]->GetMass();
};


double Particle::GetEnergy() const
{ 
    double energy = sqrt(pow(this->GetMass(), 2) + pow(this->GetImpulse(), 2));
    return energy;
};

double Particle::GetInvMass(Particle &p)
{
    double totalEnergy = pow(GetEnergy() + p.GetEnergy(), 2);
    double totalImpulse = pow(GetPx() + p.GetPx(), 2) + pow(GetPy() + p.GetPy(), 2) + pow(GetPz() + p.GetPz(), 2);
    double invMass = sqrt(totalEnergy-totalImpulse);
    return invMass;
};


void Particle::PrintParticle() const
{
    std::cout << "The particle is a " << fParticleType_[GetIndex()]->GetName() << '\n';
    std::cout << "Index is " << GetIndex() << "\n";
    std::cout << "Its impulse is:\n x: " << GetPx() << "\n y: " << GetPy()
              << "\n z: " << GetPz() << "\n\n";
};

void Particle::SetIndex(int index)
{
    if (index < fNParticleType_)
    {
        fIndex_ = index;
    }
    else
    {
        std::cout << "Error: the particle doesn't exist"<<'\n';
    }
};

void Particle::SetParticleID(std::string name)
{
    int index = FindParticle(name);
    if (index < fNParticleType_)
    {
        fIndex_ = index;
    }
    else
    {
        std::cout << "Error: the particle doesn't exist"<<'\n';
    }
};

void Particle::SetImpulse(double Px, double Py, double Pz)
{
    fPx_ = Px;
    fPy_ = Py;
    fPz_ = Pz;
};

void Particle::PrintArray()
{
    std::cout << "There are " << fNParticleType_ << " particles"<<'\n';
    for (int i = 0; i < fNParticleType_; i++)
    {
        fParticleType_[i]->Print();

    }
};


void Particle::Boost(double bx, double by, double bz)
{

    double energy = GetEnergy();

    //Boost this Lorentz vector
    double b2 = bx * bx + by * by + bz * bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx * fPx_ + by * fPy_ + bz * fPz_;
    double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

    fPx_ += gamma2 * bp * bx + gamma * bx * energy;
    fPy_ += gamma2 * bp * by + gamma * by * energy;
    fPz_ += gamma2 * bp * bz + gamma * bz * energy;
}
int Particle::Decay2body(Particle &dau1, Particle &dau2) const
{
    if (GetMass() == 0.0)
    {
        printf("Decayment cannot be preformed if mass is zero\n");
        return 1;
    }

    double massMot = GetMass();
    double massDau1 = dau1.GetMass();
    double massDau2 = dau2.GetMass();

    if (fIndex_ > -1)
    { // add width effect

        // gaussian random numbers

        float x1, x2, w, y1, y2;

        double invnum = 1. / RAND_MAX;
        do
        {
            x1 = 2.0 * rand() * invnum - 1.0;
            x2 = 2.0 * rand() * invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;
        y2 = x2 * w;

        massMot += fParticleType_[fIndex_]->GetWidth() * y1;
    }

    if (massMot < massDau1 + massDau2)
    {
        printf("Decayment cannot be preformed because mass is too low in this channel\n");
        return 2;
    }

    double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

    double norm = 2 * M_PI / RAND_MAX;

    double phi = rand() * norm;
    double theta = rand() * norm * 0.5 - M_PI / 2.;
    dau1.SetImpulse(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
    dau2.SetImpulse(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

    double energy = sqrt(fPx_ * fPx_ + fPy_ * fPy_ + fPz_ * fPz_ + massMot * massMot);

    double bx = fPx_ / energy;
    double by = fPy_ / energy;
    double bz = fPz_ / energy;

    dau1.Boost(bx, by, bz);
    dau2.Boost(bx, by, bz);

    return 0;
}



