#include "particletype.h"

#include <iostream>

ParticleType::ParticleType(std::string name, double mass, int charge) : fName_{name}, fMass_{mass}, fCharge_{charge} {};

std::string ParticleType::GetName() const { return fName_; };
double ParticleType::GetMass() const { return fMass_; };
int ParticleType::GetCharge() const { return fCharge_; };
double ParticleType::GetWidth() const { return 0; };

void ParticleType::Print() const
{
    std::cout<< "The particle is a " << fName_ <<'\n';
    std::cout<< "Particle mass is " << fMass_ << "GeV/c^2"<<'\n';
    std::cout<< "Particle charge is" << fCharge_ <<'\n';
}


