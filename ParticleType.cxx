#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include "ParticleType.h"

ParticleType::ParticleType (const char* Name,const double Mass, const int Charge)  
                 : fName_(Name), fMass_(Mass), fCharge_(Charge) {};
                 
const char* ParticleType::GetName() const { return fName_; }
double ParticleType::GetMass() const { return fMass_; }
int ParticleType::GetCharge() const { return fCharge_; }
double ParticleType::GetWidth() const { return 0;}

void ParticleType::Print() const {

    std::cout << "Name: " << fName_ << '\n'
              << "Mass: " << fMass_ << '\n'
              << "Charge: " << fCharge_ << '\n';
}
